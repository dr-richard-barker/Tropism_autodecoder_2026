#!/usr/bin/env python3
"""
Build per-study DESeq2 inputs from downloaded OSDR raw counts + harmonized metadata.

For each <counts>/OSD-<id>_raw.csv, writes into <de-dir>:
  OSD-<id>_counts_for_de.csv  (integer counts, gene x sample)
  OSD-<id>_conditions.csv     (sample, condition in {Ground Control, Space Flight})
which Code/run_de.R then consumes. Sample -> condition is resolved from
Data/harmonized_metadata.tsv by exact sample_id, then GSM, then a GC/FLT name token.

Usage:
  python Code/build_de_inputs.py --counts bulk/counts --meta Data/harmonized_metadata.tsv --de-dir bulk/de_results
"""
import argparse, glob, os, re
import pandas as pd

VALID = ("Space Flight", "Ground Control")


def gsm(s):
    m = re.search(r"GSM\d+", str(s))
    return m.group(0) if m else None


def token_cond(name):
    n = name.upper()
    if re.search(r"(^|_)(FLT|SF|SPACEFLIGHT|FLIGHT)(_|$)", n):
        return "Space Flight"
    if re.search(r"(^|_)(GC|GRC|GROUND)(_|$)", n):
        return "Ground Control"
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", default="bulk/counts", help="dir with OSD-<id>_raw.csv")
    ap.add_argument("--meta", default="Data/harmonized_metadata.tsv")
    ap.add_argument("--de-dir", default="bulk/de_results", help="output dir for DE inputs")
    args = ap.parse_args()
    os.makedirs(args.de_dir, exist_ok=True)

    m = pd.read_csv(args.meta, sep="\t")
    m["GSM"] = m["sample_id"].apply(gsm)
    sid2cond, gsm2cond = {}, {}
    for _, r in m.iterrows():
        if r["spaceflight_condition"] in VALID:
            sid2cond.setdefault(str(r["sample_id"]), r["spaceflight_condition"])
    for _, r in m.dropna(subset=["GSM"]).iterrows():
        if r["spaceflight_condition"] in VALID:
            gsm2cond.setdefault(r["GSM"], r["spaceflight_condition"])

    summary = []
    for cf in sorted(glob.glob(os.path.join(args.counts, "OSD-*_raw.csv"))):
        sid = re.search(r"(OSD-\d+)", os.path.basename(cf)).group(1)
        df = pd.read_csv(cf)
        gene_col, samples = df.columns[0], list(df.columns[1:])
        rows = [(s, sid2cond.get(s) or (gsm2cond.get(gsm(s)) if gsm(s) else None) or token_cond(s))
                for s in samples]
        cond_df = pd.DataFrame(rows, columns=["sample", "condition"]).dropna(subset=["condition"])
        nflt = int((cond_df.condition == "Space Flight").sum())
        ngc = int((cond_df.condition == "Ground Control").sum())
        if nflt < 2 or ngc < 2:
            summary.append((sid, len(samples), nflt, ngc, "SKIP <2/group")); continue
        keep = cond_df["sample"].tolist()
        out = df[[gene_col] + keep].copy()
        out.columns = ["gene"] + keep
        for c in keep:
            out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0).round().astype(int)
        out.to_csv(os.path.join(args.de_dir, f"{sid}_counts_for_de.csv"), index=False)
        cond_df.to_csv(os.path.join(args.de_dir, f"{sid}_conditions.csv"), index=False)
        summary.append((sid, len(samples), nflt, ngc, "ok"))

    print(f"{'study':10} {'cols':>5} {'FLT':>4} {'GC':>4}  status")
    for sid, n, f, g, st in summary:
        print(f"{sid:10} {n:5} {f:4} {g:4}  {st}")
    print(f"\nprepared {sum(1 for x in summary if x[4] == 'ok')} studies for DE")


if __name__ == "__main__":
    main()
