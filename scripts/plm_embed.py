#!/usr/bin/env python3
"""Mean-pooled ESM-2 embeddings for a FASTA, saved as .npy. Times itself."""
import sys, time, numpy as np, torch, esm

def main(fa, out):
    names, seqs = [], []
    k = None; buf = []
    for line in open(fa):
        if line.startswith('>'):
            if k: names.append(k); seqs.append(''.join(buf))
            k = line[1:].strip().split()[0]; buf = []
        else: buf.append(line.strip())
    if k: names.append(k); seqs.append(''.join(buf))

    dev = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
    model = model.eval().to(dev)
    bc = alphabet.get_batch_converter()
    L = model.num_layers

    t0 = time.time(); embs = []
    B = 8
    for i in range(0, len(seqs), B):
        chunk = [(n, s[:1022]) for n, s in zip(names[i:i+B], seqs[i:i+B])]
        _, _, toks = bc(chunk)
        toks = toks.to(dev)
        with torch.no_grad():
            rep = model(toks, repr_layers=[L])['representations'][L]
        for j, (n, s) in enumerate(chunk):
            n_res = min(len(s), 1022)
            embs.append(rep[j, 1:n_res+1].mean(0).float().cpu().numpy())
        if i % 80 == 0:
            print(f"  {i+len(chunk)}/{len(seqs)}", flush=True)
    dt = time.time() - t0
    E = np.vstack(embs)
    np.save(out, E)
    with open(out.replace('.npy', '.names.txt'), 'w') as fh:
        fh.write('\n'.join(names))
    print(f"device={dev}  n={len(seqs)}  dim={E.shape[1]}  embed_seconds={dt:.1f}")
    return dt

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
