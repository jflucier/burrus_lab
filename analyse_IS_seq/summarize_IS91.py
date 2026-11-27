#!/usr/bin/env python3
"""summarize_IS91.py -- table par clone + fréquences par type d'IS / passage (+ présence de Mu).
Usage: summarize_IS91.py analysis/IS91_insertions.tsv analysis/samples.tsv analysis/ [analysis/IS91_mu.tsv]
Écrit: IS91_per_clone.tsv, IS91_freq_by_passage.tsv, IS91_freq_by_rep.tsv"""
import csv, sys, os
from collections import defaultdict, Counter
ev_f, samp_f, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
mu_f = sys.argv[4] if len(sys.argv) > 4 else None
TMAP = {"pxo38km": "pOX38Km", "ecoli_chr": "chromosome"}
COMPLETE = {"simple", "tandem", "OET", "OET_candidate", "rearranged_junction", "repeat_site", "ori_side_donor", "donor_both_sides", "mu_junction", "other_junction"}

samples = {}
for r in csv.DictReader(open(samp_f, encoding="utf-8-sig"), delimiter="\t"):
    samples[r["sample"].strip()] = r
events = defaultdict(list)
for r in csv.DictReader(open(ev_f), delimiter="\t"):
    if r["event_class"] in COMPLETE and r["copy_index"] == "1":      # une ligne par événement (réseau)
        events[r["barcode"]].append(r)

mu = defaultdict(list)
if mu_f and os.path.exists(mu_f):
    for r in csv.DictReader(open(mu_f), delimiter="\t"): mu[r["barcode"]].append(r)

clones = []
for bc, s in sorted(samples.items()):
    tgt = TMAP.get(s["cible_IS"].lower(), s["cible_IS"])
    ev = events.get(bc, [])
    on = [e for e in ev if e["on_target"] == "yes"] or ev        # priorité aux événements sur la cible attendue
    on.sort(key=lambda e: -int(e["n_copies"]))
    e = on[0] if on else None
    rec = dict(barcode=bc, is_type=s["type_IS"], exp_target=tgt, passage=s["passage"],
               bio_rep=s["N_biol"], clone=s["clone"], n_events=len(ev),
               n_events_on_target=sum(x["on_target"] == "yes" for x in ev))
    if e:
        rec.update(event_class=e["event_class"], oet=e["oet"], n_copies=e["n_copies"], strand=e["strand"],
                   contig=e["contig"], arr_start=e["arr_start"], arr_end=e["arr_end"],
                   len_vs_expected=e["len_vs_expected"], target_4nt=e["target_4nt"],
                   array_open=e["array_open"], left_map=e["left_map"], right_map=e["right_map"], ref_gap=e["ref_gap"], ref_gap_seq=e["ref_gap_seq"],
                   on_target=e["on_target"])
    else:
        rec.update(event_class="none_detected", oet="NA", n_copies=0)
    m = mu.get(bc, [])
    states = sorted({x["state"] for x in m})
    rec.update(mu_present="yes" if m else "no", mu_state="+".join(states) if m else "NA", mu_n=len(m),
               mu_site=";".join(f"{x['site_ref']}:{x['site_start']}-{x['site_end']}" for x in m if x["site_ref"] != "NA") or "NA",
               mu_tsd=";".join(x["tsd_seq"] for x in m if x["tsd_seq"] != "NA") or "NA",
               mu_inverted_bp=max([int(x["inverted_bp"]) for x in m] or [0]),
               mu_snps=";".join(sorted({x.get("snp_profile", "NA") for x in m})) or "NA")
    rec["mu_dist_to_array"] = e["mu_dist"] if e else "NA"
    clones.append(rec)

cols = ["barcode", "is_type", "exp_target", "passage", "bio_rep", "clone", "n_events", "n_events_on_target",
        "event_class", "oet", "n_copies", "strand", "contig", "arr_start", "arr_end", "len_vs_expected", "array_open",
        "target_4nt", "left_map", "right_map", "ref_gap", "ref_gap_seq", "on_target",
        "mu_present", "mu_state", "mu_n", "mu_site", "mu_tsd", "mu_inverted_bp", "mu_snps", "mu_dist_to_array"]
with open(os.path.join(outdir, "IS91_per_clone.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in clones: fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")

def freq(keys, fn):
    g = defaultdict(list)
    for r in clones: g[tuple(r[k] for k in keys)].append(r)
    with open(os.path.join(outdir, fn), "w") as fh:
        fh.write("\t".join(keys + ["n_clones", "n_none", "n_simple", "n_OET", "n_OET_candidate", "n_tandem",
                                   "n_other", "pct_OET", "pct_tandem", "mean_copies_if_insertion", "n_mu", "pct_mu"]) + "\n")
        for k, rs in sorted(g.items()):
            n = len(rs); cl = Counter(r["event_class"] for r in rs)
            ins = [r for r in rs if r["event_class"] != "none_detected"]
            tand = sum(int(r["n_copies"]) > 1 for r in ins)
            oet = sum(r["oet"] == "yes" for r in ins)
            mc = sum(int(r["n_copies"]) for r in ins) / len(ins) if ins else 0
            other = len(ins) - cl["simple"] - cl["OET"] - cl["OET_candidate"] - cl["tandem"]
            fh.write("\t".join(list(k) + [str(x) for x in (n, cl["none_detected"], cl["simple"], cl["OET"],
                     cl["OET_candidate"], cl["tandem"], other)] +
                     [f"{100*oet/n:.1f}", f"{100*tand/n:.1f}", f"{mc:.2f}",
                     str(sum(r["mu_present"] == "yes" for r in rs)), f"{100*sum(r['mu_present'] == 'yes' for r in rs)/n:.1f}"]) + "\n")
freq(["is_type", "exp_target", "passage"], "IS91_freq_by_passage.tsv")
freq(["is_type", "exp_target", "passage", "bio_rep"], "IS91_freq_by_rep.tsv")
