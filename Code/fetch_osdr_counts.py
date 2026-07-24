#!/usr/bin/env python3
"""
Re-acquire bulk RNA-seq raw-count matrices from NASA OSDR.

Part of the reproduce-from-public-data path (see Code/REACQUIRE_BULK.md): the original
count-acquisition step was not committed with the pipeline, so this rebuilds it directly
from the live OSDR file API. For each OSD study it finds the study's *unnormalized* count
matrix (RSEM preferred, then STAR) and downloads it to <out>/OSD-<id>_raw.csv.

Usage:
  python Code/fetch_osdr_counts.py --out bulk/counts
  python Code/fetch_osdr_counts.py --out bulk/counts --studies 37 120 321
"""
import argparse, json, os, re, urllib.request

HOST = "https://osdr.nasa.gov"
# RNA-seq OSD studies with both Space Flight + Ground Control in harmonized_metadata.tsv.
DEFAULT_STUDIES = [37, 38, 120, 193, 217, 218, 219, 223, 281, 314, 321, 427, 437, 522, 678]


def _get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": "tropism-reacquire"})
    return urllib.request.urlopen(req, timeout=timeout)


def find_counts_url(osd):
    """Return (file_name, remote_url) of the preferred unnormalized count matrix, or (None, None)."""
    data = json.load(_get(f"{HOST}/osdr/data/osd/files/{osd}"))
    pairs = []

    def walk(o):
        if isinstance(o, dict):
            if "file_name" in o and "remote_url" in o:
                pairs.append((o["file_name"], o["remote_url"]))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(data)
    prefs = [
        lambda n: "RSEM_Unnormalized_Counts_GLbulkRNAseq.csv" in n and "rRNArm" not in n,
        lambda n: "STAR_Unnormalized_Counts_GLbulkRNAseq.csv" in n and "rRNArm" not in n,
        lambda n: re.search(r"Unnormalized_Counts.*\.csv$", n) and "rRNArm" not in n,
        lambda n: re.search(r"Unnormalized_Counts.*\.csv$", n),
    ]
    for pref in prefs:
        hits = [(n, u) for n, u in pairs if pref(n)]
        if hits:
            return hits[0]
    return None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="bulk/counts", help="output directory for OSD-<id>_raw.csv files")
    ap.add_argument("--studies", type=int, nargs="+", default=DEFAULT_STUDIES,
                    help="OSD study ids (default: the 15 RNA-seq SF+GC studies)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ok, fail = [], []
    for osd in args.studies:
        dest = os.path.join(args.out, f"OSD-{osd}_raw.csv")
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            print(f"OSD-{osd}: already have"); ok.append(osd); continue
        try:
            name, url = find_counts_url(osd)
            if not url:
                print(f"OSD-{osd}: NO unnormalized count file found"); fail.append(osd); continue
            with _get(HOST + url) as r, open(dest, "wb") as f:
                f.write(r.read())
            print(f"OSD-{osd}: {name}  ({os.path.getsize(dest)/1e6:.1f} MB)"); ok.append(osd)
        except Exception as e:
            print(f"OSD-{osd}: ERROR {e}"); fail.append(osd)
    print(f"\nDONE. ok={len(ok)} {ok}  fail={fail}")


if __name__ == "__main__":
    main()
