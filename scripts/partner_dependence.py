#!/usr/bin/env python3
"""Does a radical SAM enzyme's own sequence predict whether it has an accessory partner?

Label: for each MqnE enzyme, 1 if a DUF4229 gene lies in its canonical neighbourhood window,
0 otherwise. Features: mean-pooled ESM-2 embedding of the MqnE protein ONLY. The partner's own
sequence is never shown to the model, so the task is not trivially circular.

The threat is taxonomy: DUF4229 is concentrated in certain lineages, so a model could score well
by recognising the genus rather than any functional feature. Every evaluation here is therefore
leave-one-genus-out, and is reported beside two controls:

  genus-prior baseline : predict each held-out genus using only the training-set positive rate,
                         i.e. the best a taxonomy-only model could do. Beating this is the
                         entire claim.
  label permutation    : shuffle labels within genus. Must collapse to chance or the CV leaks.

Run: mlenv/bin/python scripts/partner_dependence.py
"""
import csv, os, sys, numpy as np
from collections import Counter
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load():
    rows = list(csv.DictReader(open(f'{ROOT}/data/coevo/mqne_all.labels.tsv'), delimiter='\t'))
    names = open(f'{ROOT}/data/coevo/emb_mqne_all.names.txt').read().split('\n')
    X = np.load(f'{ROOT}/data/coevo/emb_mqne_all.npy')
    pos = {a: i for i, a in enumerate(names)}
    keep = [r for r in rows if r['acc'] in pos]
    return (X[[pos[r['acc']] for r in keep]],
            np.array([int(r['label']) for r in keep]),
            np.array([r['genus'] for r in keep]))


def grouped_auc(X, y, g, min_per_genus=8, seed=0):
    """Leave-one-genus-out. Only genera with both classes present can be scored."""
    cnt = Counter(g)
    big = [k for k, v in cnt.items() if v >= min_per_genus]
    mask = np.isin(g, big)
    X, y, g = X[mask], y[mask], g[mask]
    logo = LeaveOneGroupOut()
    oof = np.full(len(y), np.nan)
    prior = np.full(len(y), np.nan)
    for tr, te in logo.split(X, y, g):
        if y[tr].sum() < 5 or y[tr].sum() == len(tr):
            continue
        sc = StandardScaler().fit(X[tr])
        clf = LogisticRegression(max_iter=3000, C=0.05, class_weight='balanced')
        clf.fit(sc.transform(X[tr]), y[tr])
        oof[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
        prior[te] = y[tr].mean()          # taxonomy-blind training prior
    ok = ~np.isnan(oof)
    return oof, prior, ok, y, g


def main():
    X, y, g = load()
    print(f"n = {len(y)}   positives = {y.sum()} ({100*y.mean():.1f}%)   genera = {len(set(g))}")

    oof, prior, ok, y2, g2 = grouped_auc(X, y, g)
    n_ev = ok.sum()
    genera_ev = len(set(g2[ok]))
    auc = roc_auc_score(y2[ok], oof[ok])
    print(f"\nevaluated on {n_ev} enzymes across {genera_ev} held-out genera")
    print(f"  ESM-2 sequence model, leave-one-genus-out AUC = {auc:.4f}")

    # control 1: best a taxonomy-only model could do on the same folds
    try:
        auc_prior = roc_auc_score(y2[ok], prior[ok])
    except ValueError:
        auc_prior = float('nan')
    print(f"  genus-prior baseline (taxonomy only)      = {auc_prior:.4f}   <- must be beaten")

    # control 2: labels shuffled within genus
    rng = np.random.default_rng(0); ys = y.copy()
    for gen in set(g):
        m = g == gen
        ys[m] = rng.permutation(y[m])
    oof_s, _, ok_s, ys2, _ = grouped_auc(X, ys, g)
    print(f"  label permuted within genus               = {roc_auc_score(ys2[ok_s], oof_s[ok_s]):.4f}   <- must be ~0.50")

    # per-genus breakdown, so a single dominant genus cannot carry the result
    print("\n  per-genus AUC (genera with both classes and n >= 15):")
    for gen in sorted(set(g2[ok])):
        m = ok & (g2 == gen)
        if m.sum() < 15 or len(set(y2[m])) < 2:
            continue
        print(f"    {gen:<22}{m.sum():>5}  pos={y2[m].sum():>4}  AUC={roc_auc_score(y2[m], oof[m]):.3f}")


if __name__ == '__main__':
    main()
