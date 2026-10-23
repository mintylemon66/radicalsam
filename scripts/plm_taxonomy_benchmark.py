#!/usr/bin/env python3
"""How much of a protein language model embedding is taxonomy, and can it be removed?

Protein language models are increasingly used as general-purpose protein representations in
comparative genomics. If their embeddings encode the SPECIES a protein came from as strongly as
they encode what the protein DOES, then any cross-genome analysis built on them inherits a
phylogenetic confound, and ordinary permutation nulls will not detect it.

This benchmark measures both quantities on the same proteins, for several representations:

    TAXONOMY signal : balanced accuracy predicting genus (and order) from the embedding
    FUNCTION signal : balanced accuracy predicting Pfam family from the embedding

A useful correction drives TAXONOMY down while leaving FUNCTION intact. Both are reported
against a deliberately trivial baseline, 20 amino-acid frequencies plus log length, because a
learned representation that cannot beat residue counting has not earned its runtime.

All cross-validation is grouped by CD-HIT sequence-identity cluster, so near-identical
homologs never span train and test.

Run: mlenv/bin/python scripts/plm_taxonomy_benchmark.py [--models esm2_t30_150M_UR50D,...]
"""
import argparse, csv, json, os, subprocess, sys, time
from collections import Counter
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REC = f'{ROOT}/data/plmbench/records.jsonl'
AA = 'ACDEFGHIKLMNPQRSTVWY'
MIN_PER_CLASS = 30


def load(min_genus=MIN_PER_CLASS):
    rows = [json.loads(l) for l in open(REC) if l.strip()]
    gc = Counter(r['genus'] for r in rows)
    rows = [r for r in rows if gc[r['genus']] >= min_genus]
    fc = Counter(r['pfam'] for r in rows)
    rows = [r for r in rows if fc[r['pfam']] >= min_genus]
    return rows


def cluster(rows, ident=0.7):
    """CD-HIT clusters as CV groups so homologs cannot span folds."""
    fa, out = f'{ROOT}/data/plmbench/_all.faa', f'{ROOT}/data/plmbench/_cl{int(ident*100)}'
    with open(fa, 'w') as fh:
        for r in rows:
            fh.write(f">{r['acc']}\n{r['seq']}\n")
    if not os.path.exists(out + '.clstr'):
        n = 5 if ident >= 0.7 else 4
        subprocess.run(['cd-hit', '-i', fa, '-o', out, '-c', str(ident), '-n', str(n),
                        '-d', '0', '-M', '2000'], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    import re
    cl, cur = {}, None
    for line in open(out + '.clstr'):
        if line.startswith('>Cluster'):
            cur = line.split()[1]
        else:
            m = re.search(r'>(\S+?)\.\.\.', line)
            if m:
                cl[m.group(1)] = cur
    return cl


def order_of(lineage):
    for t in lineage:
        if t.endswith('ales'):
            return t
    return lineage[-1] if lineage else 'unknown'


def composition(rows):
    X = []
    for r in rows:
        s = r['seq']; n = max(len(s), 1)
        X.append([s.count(c) / n for c in AA] + [np.log(n)])
    return np.array(X)


def embed(rows, model_name):
    npy = f"{ROOT}/data/plmbench/emb_{model_name}.npy"
    if os.path.exists(npy):
        return np.load(npy), None
    import torch, esm
    dev = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model, alphabet = getattr(esm.pretrained, model_name)()
    model = model.eval().to(dev)
    bc = alphabet.get_batch_converter(); L = model.num_layers
    t0 = time.time(); out = []
    B = 4 if '650M' in model_name else 8
    for i in range(0, len(rows), B):
        ch = [(r['acc'], r['seq'][:1022]) for r in rows[i:i+B]]
        _, _, tok = bc(ch)
        with torch.no_grad():
            rep = model(tok.to(dev), repr_layers=[L])['representations'][L]
        for j, (_, s) in enumerate(ch):
            out.append(rep[j, 1:min(len(s), 1022)+1].mean(0).float().cpu().numpy())
        if i % 800 == 0:
            print(f"    {i+len(ch)}/{len(rows)}", flush=True)
    dt = time.time() - t0
    E = np.vstack(out); np.save(npy, E)
    print(f"    {model_name}: {len(rows)} seqs in {dt:.0f}s on {dev}", flush=True)
    return E, dt


def bacc(X, y, groups, seed=0, folds=5):
    o = []
    for tr, te in StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed).split(X, y, groups):
        if len(set(y[tr])) < 2:
            continue
        sc = StandardScaler().fit(X[tr])
        m = LogisticRegression(max_iter=3000, C=0.1, class_weight='balanced')
        m.fit(sc.transform(X[tr]), y[tr])
        o.append(balanced_accuracy_score(y[te], m.predict(sc.transform(X[te]))))
    return float(np.mean(o)) if o else float('nan')


def center_within(M, labels, minn=3):
    M = M.copy(); gm = M.mean(0)
    lab = np.asarray(labels)
    for u in set(lab):
        i = lab == u
        M[i] = M[i] - M[i].mean(0) if i.sum() >= minn else M[i] - gm
    return M


def rm_pcs(M, k):
    M = M - M.mean(0)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    return M - (U[:, :k] * S[:k]) @ Vt[:k]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--models', default='esm2_t12_35M_UR50D,esm2_t30_150M_UR50D')
    a = ap.parse_args()

    rows = load()
    genus = np.array([r['genus'] for r in rows])
    fam = np.array([r['pfam'] for r in rows])
    order = np.array([order_of(r['lineage']) for r in rows])
    cl = cluster(rows)
    grp = np.array([cl.get(r['acc'], r['acc']) for r in rows])
    print(f"n={len(rows)}  families={len(set(fam))}  genera={len(set(genus))}  "
          f"orders={len(set(order))}  clusters={len(set(grp))}")
    print(f"chance: genus {1/len(set(genus)):.3f}  order {1/len(set(order)):.3f}  family {1/len(set(fam)):.3f}\n")

    res = []
    def row(tag, X):
        g, o, f = bacc(X, genus, grp), bacc(X, order, grp), bacc(X, fam, grp)
        print(f"  {tag:<34}{g:>9.3f}{o:>9.3f}{f:>9.3f}")
        res.append({'representation': tag, 'genus': round(g,4), 'order': round(o,4), 'family': round(f,4)})

    print(f"  {'representation':<34}{'GENUS':>9}{'ORDER':>9}{'FAMILY':>9}")
    print('  ' + '-' * 61)
    row('amino-acid composition (trivial)', composition(rows))
    for mn in a.models.split(','):
        short = mn.replace('esm2_', '').replace('_UR50D', '')
        E, _ = embed(rows, mn)
        row(f'{short} raw', E)
        row(f'{short} + within-genus centring', center_within(E, genus))
        row(f'{short} + within-order centring', center_within(E, order))
        row(f'{short} + top-10 PCs removed', rm_pcs(E, 10))

    with open(f'{ROOT}/audit/plm_taxonomy_benchmark.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=['representation','genus','order','family'])
        w.writeheader(); w.writerows(res)
    print("\nwrote audit/plm_taxonomy_benchmark.tsv")


if __name__ == '__main__':
    main()
