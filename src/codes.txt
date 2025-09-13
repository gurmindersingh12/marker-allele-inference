import argparse
import pandas as pd
from collections import Counter

BASES = set("ACGT")
IUPAC = {
    "A":{"A"}, "C":{"C"}, "G":{"G"}, "T":{"T"},
    "R":{"A","G"}, "Y":{"C","T"}, "S":{"G","C"},
    "W":{"A","T"}, "K":{"G","T"}, "M":{"A","C"},
    "B":{"C","G","T"}, "D":{"A","G","T"}, "H":{"A","C","T"},
    "V":{"A","C","G"}, "N":set()
}
# You said you'll remove HH; we also ignore it just in case.
MISSING = {"", ".", "-", "NA", "N", "NN", "HH"}

def parse_cell(cell):
    """Return (alleles_set, is_homozygote?)."""
    if pd.isna(cell): return set(), None
    s = str(cell).strip().upper()
    if s in MISSING: return set(), None

    # Two-letter encodings like AA, AT, TT
    if len(s) == 2 and s[0] in BASES and s[1] in BASES:
        if s[0] == s[1]:
            return {s[0]}, True
        return {s[0], s[1]}, False

    # Single-letter IUPAC
    if len(s) == 1 and s in IUPAC:
        return IUPAC[s] & BASES, False

    return set(), None

def infer_alleles(row_vals):
    """Prefer homozygotes. If not enough, consider all calls."""
    hom = Counter()
    anyc = Counter()

    for v in row_vals:
        alleles, is_hom = parse_cell(v)
        for b in alleles:
            anyc[b] += 1
            if is_hom:
                hom[b] += 1

    uniq_any = [b for b,c in anyc.items() if c > 0]
    n_unique = len(uniq_any)
    tri_plus = n_unique > 2

    # zero or mono
    if n_unique == 0:
        return "", 0, False
    if n_unique == 1:
        b = sorted(uniq_any)[0]
        return f"{b}/{b}", 1, False

    # counts to choose top two (prefer homozygotes if we saw >=2 bases as homozygotes)
    if sum(c > 0 for c in hom.values()) >= 2:
        counts = hom
    else:
        counts = anyc

    top2 = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:2]
    alleles = "/".join(sorted([top2[0][0], top2[1][0]]))
    return alleles, n_unique, tri_plus

def main():
    ap = argparse.ArgumentParser(description="Infer per-marker allele pair from a genotype matrix and output a tidy CSV with QC flags.")
    ap.add_argument("--in", dest="infile", required=True, help="Input CSV file")
    ap.add_argument("--out", dest="outfile", default="marker_with_alleles.csv", help="Output CSV file")
    ap.add_argument("--marker-col", default=None, help="Existing marker ID column name (optional)")
    ap.add_argument("--drop-multi", action="store_true", help="Drop markers with >2 alleles")
    ap.add_argument("--no-clean-headers", dest="clean_headers", action="store_false",
                    help="Don't trim headers at semicolons")
    ap.set_defaults(clean_headers=True)
    args = ap.parse_args()

    # Load as strings to avoid accidental type coercion
    df = pd.read_csv(args.infile, dtype=str, keep_default_na=False)

    # Optionally clean headers like "CLOUSEAU;IonCode_0701" -> "CLOUSEAU"
    if args.clean_headers:
        df.columns = [c.split(";")[0] if isinstance(c, str) else c for c in df.columns]

    # Figure out metadata vs genotype columns
    meta_cols = []
    if args.marker_col and args.marker_col in df.columns:
        meta_cols = [args.marker_col]
    else:
        # If no marker ID column present, create one
        if not any(str(c).lower() in {"marker", "id", "rs#"} for c in df.columns):
            df.insert(0, "Marker", [f"marker_{i+1:05d}" for i in range(len(df))])
            meta_cols = ["Marker"]
        else:
            # Use whichever is present
            for c in df.columns:
                if str(c).lower() in {"marker", "id", "rs#"}:
                    meta_cols = [c]; break

    geno_cols = [c for c in df.columns if c not in meta_cols]

    # Infer
    alle_info = df[geno_cols].apply(lambda r: infer_alleles(r.values), axis=1)
    df.insert(1, "alleles", alle_info.map(lambda t: t[0]))       # e.g., "A/T"
    df.insert(2, "n_alleles", alle_info.map(lambda t: t[1]))     # number of distinct bases seen
    df.insert(3, "multi_allelic", alle_info.map(lambda t: "YES" if t[2] else ""))

    if args.drop_multi:
        df = df[df["multi_allelic"] != "YES"].copy()

    df.to_csv(args.outfile, index=False)
    # Small summary
    total = len(df)
    multi = (df["multi_allelic"] == "YES").sum() if "multi_allelic" in df else 0
    print(f"Wrote {args.outfile} | rows: {total} | multi-allelic flagged: {multi}")

if __name__ == "__main__":
    main()
