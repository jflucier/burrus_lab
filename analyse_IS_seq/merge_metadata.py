#!/usr/bin/env python3
"""Usage: merge_meta.py samples.tsv contig_assignment.tsv > contig_assignment_meta.tsv
Joint les métadonnées de samples.tsv (colonne 'sample') sur la colonne 'barcode'."""
import csv, sys

samples_f, table_f = sys.argv[1], sys.argv[2]

TARGET = {"pxo38km": "pOX38Km", "pox38km": "pOX38Km",
          "ecoli_chr": "chromosome", "chromosome": "chromosome"}
META = [("is_type", "type_IS"), ("exp_target", "cible_IS"), ("box", "box_localisation"),
        ("bio_rep", "N_biol"), ("clone", "clone"), ("passage", "passage")]
OLD = {"Type d'IS", "Cible de l'IS"}   # colonnes ajoutées à la main : remplacées

with open(samples_f, encoding="utf-8-sig") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    missing = [c for _, c in META if c not in rd.fieldnames] + \
              (["sample"] if "sample" not in rd.fieldnames else [])
    if missing:
        sys.exit(f"Colonnes absentes de {samples_f}: {missing} (trouvées: {rd.fieldnames})")
    meta = {r["sample"].strip(): r for r in rd}

with open(table_f, encoding="utf-8-sig") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    cols = [c for c in rd.fieldnames if c not in OLD and c != "barcode"]
    print("\t".join(["barcode"] + [n for n, _ in META] + cols + ["on_target"]))
    unknown = set()
    for r in rd:
        bc = r["barcode"].strip()
        m = meta.get(bc)
        if m is None:
            unknown.add(bc)
        g = lambda c: ((m or {}).get(c, "") or "").strip() or "NA"
        vals = [g(c) for _, c in META]
        vals[1] = TARGET.get(vals[1].lower(), vals[1])        # exp_target normalisée
        asg = r.get("assignation", "NA")
        on = "NA" if vals[1] == "NA" else ("yes" if vals[1] == asg else "no")
        print("\t".join([bc] + vals + [r[c] for c in cols] + [on]))
    if unknown:
        print("Barcodes absents de samples.tsv: " + ", ".join(sorted(unknown)), file=sys.stderr)