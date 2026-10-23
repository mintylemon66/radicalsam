#!/usr/bin/env python3
"""Fetch a taxonomy- and family-labelled bacterial protein benchmark from UniProt.

Each record is (sequence, Pfam family, genus, NCBI lineage): exactly the label pair needed to
measure how much of a protein language model embedding encodes taxonomy versus function.
Uses curl, because bare urllib stalls against this endpoint on this machine.
"""
import json, os, subprocess, sys, time

FAMILIES = ["PF04055","PF00005","PF00072","PF01037","PF02441","PF00583","PF01040","PF01977",
            "PF02621","PF14012","PF07293","PF08812","PF07315","PF02681","PF00375","PF13396",
            "PF00106","PF00561","PF00196","PF01047","PF00392","PF02518","PF00534","PF00990",
            "PF00701","PF00445","PF00296","PF00135","PF00126","PF00589"]
CAP = int(os.environ.get("CAP", "250"))
OUT = "data/plmbench/records.jsonl"


def fetch(pfam):
    q = (f"xref:pfam-{pfam}%20AND%20taxonomy_id:2%20AND%20length:%5B60%20TO%20700%5D")
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
           f"&fields=accession,sequence,organism_name&format=json&size={CAP}")
    tmp = f"/tmp/_pf_{pfam}.json"
    for attempt in range(3):
        r = subprocess.run(["curl", "-s", "--max-time", "120", "-A", "Mozilla/5.0 research",
                            url, "-o", tmp], capture_output=True)
        if r.returncode == 0 and os.path.exists(tmp) and os.path.getsize(tmp) > 200:
            break
        time.sleep(4)
    else:
        return []
    try:
        data = json.load(open(tmp))
    except Exception:
        return []
    out = []
    for it in data.get("results", []):
        seq = (it.get("sequence") or {}).get("value")
        org = it.get("organism") or {}
        name, lin = org.get("scientificName", ""), list(org.get("lineage") or [])
        if not seq or not name or not lin:
            continue
        out.append({"acc": it["primaryAccession"], "pfam": pfam, "seq": seq,
                    "organism": name, "genus": name.split()[0], "lineage": lin})
    os.remove(tmp)
    return out


def main():
    os.makedirs("data/plmbench", exist_ok=True)
    total = 0
    with open(OUT, "w") as fh:
        for pf in FAMILIES:
            rows = fetch(pf)
            for r in rows:
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            total += len(rows)
            g = len({r["genus"] for r in rows})
            print(f"{pf}: {len(rows):>4} proteins, {g:>3} genera   (total {total})", flush=True)
            time.sleep(0.3)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
