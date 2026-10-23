#!/usr/bin/env python3
"""Phylogeny-controlled, MSA-free co-evolution screen.

For each neighbour family F of MqnE, ask whether F's sequence variation tracks MqnE's
variation ACROSS GENOMES, after removing the component explained by species phylogeny.

    D_MqnE, D_F  : pairwise cosine distance between per-dimension z-scored ESM-2 embeddings
    D_tax        : taxonomic distance from shared NCBI lineage depth
    statistic    : partial Mantel r(D_MqnE, D_F | D_tax), permutation p

The z-scoring is not cosmetic. Mean-pooled ESM-2 embeddings of a conserved enzyme family are
dominated by a large shared component; without standardisation the usable dynamic range collapses
(mean cosine distance 0.0041 for MqnE) and the statistic loses power. See coevo_report.md.

DUF4229 is the candidate. The other families are the control panel: they occupy the same genomes
and the same neighbourhood window, so any statistic that cannot separate DUF4229 from them is
reporting genome co-membership rather than functional coupling.

Run: mlenv/bin/python scripts/coevo_panel.py
"""
import csv, json, os, sys, time
import numpy as np
from scipy.spatial.distance import pdist, squareform

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NPERM = 9999


def read_fasta(p):
    d, k, buf = {}, None, []
    for line in open(p):
        if line.startswith('>'):
            if k: d[k] = ''.join(buf)
            k = line[1:].strip().split()[0]; buf = []
        else: buf.append(line.strip())
    if k: d[k] = ''.join(buf)
    return d


def zscore(M):
    return (M - M.mean(0)) / (M.std(0) + 1e-8)


def partial_mantel(A, B, dt, nperm=NPERM, seed=1):
    d1, d2 = pdist(zscore(A), 'cosine'), pdist(zscore(B), 'cosine')
    r12 = np.corrcoef(d1, d2)[0, 1]
    r13 = np.corrcoef(d1, dt)[0, 1]
    r23 = np.corrcoef(d2, dt)[0, 1]
    den = np.sqrt((1 - r13**2) * (1 - r23**2))
    pr = (r12 - r13 * r23) / den
    X = squareform(d1); n = A.shape[0]; rg = np.random.default_rng(seed); c = 0
    for _ in range(nperm):
        p = rg.permutation(n)
        dp = squareform(X[np.ix_(p, p)], checks=False)
        a12 = np.corrcoef(dp, d2)[0, 1]; a13 = np.corrcoef(dp, dt)[0, 1]
        if (a12 - a13 * r23) / np.sqrt((1 - a13**2) * (1 - r23**2)) >= pr:
            c += 1
    return r12, pr, (c + 1) / (nperm + 1), r13, r23


def tax_matrix(lineages):
    n = len(lineages); D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = lineages[i], lineages[j]; c = 0
            for x, y in zip(a, b):
                if x == y: c += 1
                else: break
            D[i, j] = D[j, i] = 1 - c / max(len(a), len(b))
    return squareform(D, checks=False)


def embed_cached(fa, npy):
    if not os.path.exists(npy):
        sys.path.insert(0, os.path.join(ROOT, 'scripts'))
        import plm_embed
        plm_embed.main(fa, npy)
    names = open(npy.replace('.npy', '.names.txt')).read().split('\n')
    return np.load(npy), names


def main():
    T = json.load(open(f'{ROOT}/data/coevo/taxonomy.json'))
    acc_tax, lin = T['acc_tax'], T['lineage']
    idx = list(csv.DictReader(open(f'{ROOT}/data/coevo/pairs_index.tsv'), delimiter='\t'))
    mqn_all, mqn_names = embed_cached(f'{ROOT}/data/coevo/pairs_mqn.faa', f'{ROOT}/data/coevo/emb_mqn.npy')
    pos = {a: i for i, a in enumerate(mqn_names)}

    fams = [('PF14012', 'DUF4229 (candidate)', f'{ROOT}/data/coevo/pairs_duf.faa', None)]
    cdir = f'{ROOT}/data/coevo/controls'
    labels = {'PF01037': 'AsnC_trans_reg', 'PF02441': 'Flavoprotein', 'PF00583': 'Acetyltransf_1',
              'PF01040': 'UbiA', 'PF01977': 'UbiD'}
    for f in sorted(labels):
        fa = f'{cdir}/{f}.faa'
        if os.path.exists(fa) and os.path.exists(f'{cdir}/{f}.map.json'):
            fams.append((f, labels[f] + ' (control)', fa, f'{cdir}/{f}.map.json'))

    print(f"{'family':<10}{'label':<24}{'n':>5}{'raw r':>9}{'partial r':>11}{'p':>10}{'r_tax':>8}")
    print('-' * 77)
    out = []
    for acc, label, fa, mapf in fams:
        seqs = read_fasta(fa)
        if mapf:
            m = json.load(open(mapf))
            keep = [(a, pid) for a, pid in m.items() if a in pos and lin.get(acc_tax.get(a, ''))
                    and (pid in seqs or pid.split('.')[0] in seqs)]
        else:
            keep = [(r['mqn_acc'], r['duf_protein_id']) for r in idx
                    if r['mqn_acc'] in pos and lin.get(acc_tax.get(r['mqn_acc'], ''))]
        if len(keep) < 40:
            print(f"{acc:<10}{label:<24}{len(keep):>5}   skipped (n too small)"); continue
        keep.sort()
        A = mqn_all[[pos[a] for a, _ in keep]]
        tmp_fa = f'{cdir}/_{acc}_sub.faa'
        with open(tmp_fa, 'w') as fh:
            for a, pid in keep:
                s = seqs.get(pid) or seqs.get(pid.split('.')[0])
                fh.write(f">{a}\n{s}\n")
        B, _ = embed_cached(tmp_fa, f'{cdir}/_{acc}_sub.npy')
        dt = tax_matrix([lin[acc_tax[a]] for a, _ in keep])
        r12, pr, p, r13, r23 = partial_mantel(A, B, dt)
        star = '  <==' if acc == 'PF14012' else ''
        print(f"{acc:<10}{label:<24}{len(keep):>5}{r12:>9.4f}{pr:>11.4f}{p:>10.4f}{r23:>8.3f}{star}")
        out.append({'pfam': acc, 'label': label, 'n': len(keep), 'raw_r': round(r12, 4),
                    'partial_r': round(pr, 4), 'p': p, 'r_family_vs_taxonomy': round(r23, 3)})
    with open(f'{ROOT}/audit/coevo_panel.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    print(f"\nwrote audit/coevo_panel.tsv")


if __name__ == '__main__':
    main()
