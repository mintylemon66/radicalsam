#!/usr/bin/env python3
"""Are the standard cross-validation protocols calibrated on taxonomically structured data?

A protocol is CALIBRATED if, on data containing no real signal, it returns AUC ~ 0.5.
Signal-free data is built by permuting labels WITHIN each genus: this destroys any true
sequence-to-label relationship while preserving the taxonomic structure of the labels exactly.

Protocols compared, in increasing order of how careful the field considers them:
  A random split,        pooled AUC   - the naive default
  B identity-grouped,    pooled AUC   - CD-HIT clusters as CV groups (the usual leakage control)
  C genus-grouped,       pooled AUC   - hold out whole genera (the stricter fix people reach for)
  D genus-grouped, within-genus AUC   - hold out genera AND score inside each genus, then average

If C is inflated while D is calibrated, then grouping the SPLIT is not sufficient: a pooled
readout re-imports the very confound the split removed, because between-genus differences in
label prevalence are still rankable.

Features are amino-acid composition + log length, deliberately trivial, so that nothing here
depends on a language model.

Run: mlenv/bin/python scripts/cv_protocol_calibration.py
"""
import csv, os, re, subprocess, sys
from collections import Counter, defaultdict
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AA = 'ACDEFGHIKLMNPQRSTVWY'
NPERM = 30
MIN_POS = 60


def is_anchor(r):
    try:
        return (abs(int(r["n_start"]) - int(r["q_start"])) <= 1
                and abs(int(r["n_end"]) - int(r["q_end"])) <= 1)
    except (TypeError, ValueError, KeyError):
        return False


def load_tasks():
    sys.path.insert(0, f'{ROOT}/scripts')
    import nacs_ranking as N
    anch, genus = N.load_anchor_map(), N.load_genus()
    seqs, k, buf = {}, None, []
    for line in open(f'{ROOT}/data/mqn_targets.fasta'):
        if line.startswith('>'):
            if k: seqs[k] = ''.join(buf)
            m = re.match(r'^>(?:sp|tr)\|([^|]+)\|', line); k = m.group(1) if m else None; buf = []
        else: buf.append(line.strip())
    if k: seqs[k] = ''.join(buf)

    hits, anchors = defaultdict(set), set()
    for line in open(f'{ROOT}/audit/pfam_wide.domtbl'):
        if line.startswith('#'): continue
        p = line.split()
        if len(p) < 23: continue
        acc = p[4].split('.')[0] if p[4] != '-' else p[3]
        d = p[-1].split('|')
        if len(d) < 4: continue
        q, grp = d[0], d[1]
        if grp != 'mqn' or q not in seqs or q not in genus: continue
        try:
            gi = int(d[2].split('=')[1]); gap = int(d[3].split('=')[1])
        except (IndexError, ValueError): continue
        anchors.add(q)
        # exclude the anchor gene itself, by coordinates AND by offset zero
        if gap > 8000 or abs(gi) > 10 or gi == 0 or anch.get((q, p[0]), False): continue
        hits[acc].add(q)

    anchors = sorted(anchors)
    idx = {a: i for i, a in enumerate(anchors)}
    X = np.array([[seqs[a].count(c) / max(len(seqs[a]), 1) for c in AA] + [np.log(len(seqs[a]))]
                  for a in anchors])
    gen = np.array([genus[a] for a in anchors])
    tasks = []
    for acc, s in hits.items():
        y = np.zeros(len(anchors), int)
        for q in s:
            y[idx[q]] = 1
        if MIN_POS <= y.sum() <= len(y) - MIN_POS:
            # need genera that actually contain both classes, or nothing is measurable
            ok = sum(1 for g in set(gen)
                     if (gen == g).sum() >= 20 and 0 < y[gen == g].sum() < (gen == g).sum())
            if ok >= 3:
                tasks.append((acc, y, ok))
    tasks.sort(key=lambda t: -t[2])
    return X, gen, anchors, seqs, tasks


def clusters(anchors, seqs, ident=0.5):
    fa, out = f'{ROOT}/data/plmbench/_anch.faa', f'{ROOT}/data/plmbench/_anch50'
    os.makedirs(f'{ROOT}/data/plmbench', exist_ok=True)
    with open(fa, 'w') as fh:
        for a in anchors: fh.write(f'>{a}\n{seqs[a]}\n')
    if not os.path.exists(out + '.clstr'):
        subprocess.run(['cd-hit', '-i', fa, '-o', out, '-c', str(ident), '-n', '2', '-d', '0',
                        '-M', '2000'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    cl, cur = {}, None
    for line in open(out + '.clstr'):
        if line.startswith('>Cluster'): cur = line.split()[1]
        else:
            m = re.search(r'>(\S+?)\.\.\.', line)
            if m: cl[m.group(1)] = cur
    return np.array([cl.get(a, a) for a in anchors])


def fit_predict(X, y, tr, te):
    if len(set(y[tr])) < 2: return None
    sc = StandardScaler().fit(X[tr])
    m = LogisticRegression(max_iter=2000, C=0.1, class_weight='balanced')
    m.fit(sc.transform(X[tr]), y[tr])
    return m.predict_proba(sc.transform(X[te]))[:, 1]


def protocol(X, y, gen, clu, mode, seed=0):
    n = len(y); oof = np.full(n, np.nan)
    if mode == 'A':   splits = StratifiedKFold(5, shuffle=True, random_state=seed).split(X, y)
    elif mode == 'B': splits = StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(X, y, clu)
    else:             splits = StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(X, y, gen)
    for tr, te in splits:
        p = fit_predict(X, y, tr, te)
        if p is not None: oof[te] = p
    ok = ~np.isnan(oof)
    if len(set(y[ok])) < 2: return np.nan
    if mode != 'D':
        return roc_auc_score(y[ok], oof[ok])
    aucs, w = [], []
    for g in set(gen):
        m = ok & (gen == g)
        if m.sum() >= 20 and 0 < y[m].sum() < m.sum():
            aucs.append(roc_auc_score(y[m], oof[m])); w.append(m.sum())
    return float(np.average(aucs, weights=w)) if aucs else np.nan


def main():
    X, gen, anchors, seqs, tasks = load_tasks()
    clu = clusters(anchors, seqs)
    print(f"anchors={len(anchors)}  genera={len(set(gen))}  cd-hit@50 clusters={len(set(clu))}")
    print(f"admissible tasks={len(tasks)}; using top {min(8,len(tasks))}\n")
    names = {'A':'random split,      pooled AUC','B':'identity-grouped,  pooled AUC',
             'C':'genus-grouped,     pooled AUC','D':'genus-grouped, within-genus AUC'}
    agg = {m: {'null': [], 'real': []} for m in 'ABCD'}
    rng = np.random.default_rng(0)
    for acc, y, ok in tasks[:8]:
        for m in 'ABCD':
            agg[m]['real'].append(protocol(X, y, gen, clu, m))
        for _ in range(NPERM):
            ys = y.copy()
            for g in set(gen):
                i = gen == g
                ys[i] = rng.permutation(y[i])          # signal-free, genus structure preserved
            for m in 'ABCD':
                agg[m]['null'].append(protocol(X, ys, gen, clu, m))
        print(f"  task {acc}: pos={y.sum()}/{len(y)}, {ok} mixed genera", flush=True)

    print(f"\n{'protocol':<34}{'null median':>12}{'null p95':>10}{'FPR@.55':>9}{'real median':>13}")
    print('-' * 79)
    rows = []
    for m in 'ABCD':
        nl = np.array([v for v in agg[m]['null'] if not np.isnan(v)])
        rl = np.array([v for v in agg[m]['real'] if not np.isnan(v)])
        fpr = float((nl > 0.55).mean())
        print(f"{names[m]:<34}{np.median(nl):>12.3f}{np.percentile(nl,95):>10.3f}"
              f"{fpr:>9.1%}{np.median(rl):>13.3f}")
        rows.append({'protocol': names[m], 'null_median': round(float(np.median(nl)),4),
                     'null_p95': round(float(np.percentile(nl,95)),4),
                     'fpr_at_0.55': round(fpr,4), 'real_median': round(float(np.median(rl)),4)})
    with open(f'{ROOT}/audit/cv_protocol_calibration.tsv','w',newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\nA calibrated protocol has null median ~0.50 and low FPR while keeping real median high.")
    print("wrote audit/cv_protocol_calibration.tsv")


if __name__ == '__main__':
    main()
