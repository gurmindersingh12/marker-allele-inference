# Code Explanation (Design & Logic)

This document explains how the script infers per-marker alleles and why certain design choices were made.

## Overview

For each marker (row) across all samples (columns), the script determines which **two alleles** best represent that locus and reports them as `A/T` (or `C/C` if monomorphic). It also computes simple QC metadata: the number of distinct bases (`n_alleles`) and whether the site appears tri-/multi-allelic.

## Parsing a genotype cell: `parse_cell`

- Accepts typical two-letter encodings (`AA`, `AT`, `TT`).
  - `AA`, `CC`, `GG`, `TT` are treated as **homozygotes**.
  - `AT`, `AG`, `CG`, etc., are **heterozygotes**.
- Accepts IUPAC single-letter codes:
  - R(A/G), Y(C/T), S(G/C), W(A/T), K(G/T), M(A/C), B(C/G/T), D(A/G/T), H(A/C/T), V(A/C/G), N(missing).
  - Single letters are treated as **not explicitly homozygous** by default (conservative across vendor formats).
- Recognizes missing tokens: `""`, `.`, `-`, `NA`, `N`, `NN`, `HH`.

Returns a tuple `(allele_set, is_homozygote)`:
- `allele_set` is a set of bases among `{A,C,G,T}` found/encoded in the cell (possibly empty).
- `is_homozygote` is `True` for explicit `AA/CC/..`, `False` for heterozygous or ambiguous/IUPAC, and `None` for missing/unknown.

## Deciding alleles per marker: `infer_alleles`

We use two counters:
- `hom`: counts alleles seen in **homozygote** calls only.
- `anyc`: counts alleles seen in **any** call (homo+hetero+IUPAC).

Algorithm:
1. Scan all sample calls for the marker, update both counters.
2. Let `uniq_any` be the list of alleles with non-zero counts in `anyc`; `n_unique = len(uniq_any)`.
3. If `n_unique == 0`, return `("", 0, False)` (no calls).
4. If `n_unique == 1`, return `("X/X", 1, False)` for the single allele `X` (monomorphic).
5. Otherwise (`n_unique >= 2`):
   - If we observed **two or more distinct alleles as homozygotes**, choose the top two by **homozygote counts** (more robust).
   - Else, fall back to **overall counts** in `anyc`.
   - Sort the two bases lexicographically and format as `A/T`.
6. Flag `tri_plus = (n_unique > 2)` for QA/QC.

This approach avoids letting ambiguous heterozygote/IUPAC noise incorrectly introduce a third allele when there is solid homozygote evidence for the biallelic pair.

## CLI & I/O

- Reads input CSV as **strings** to avoid type coercion.
- Optionally trims sample headers at semicolons (e.g., `SampleA;IonCode_0701` → `SampleA`) to remove barcode tails.
- If there is no explicit marker ID column (`Marker`, `ID`, or `rs#`), it auto-inserts a `Marker` column.
- Writes output CSV with `alleles`, `n_alleles`, `multi_allelic` prepended after the ID column.

## Edge cases & notes

- **Single-letter base calls (A/C/G/T)**: treated as not explicitly homozygous by default. If your pipeline guarantees that single letters denote homozygotes, change `parse_cell` accordingly.
- **Tri-allelic sites**: These are flagged. Use `--drop-multi` to exclude them for biallelic-only pipelines.
- **Missing-heavy markers**: With very sparse data, the script may return empty alleles and `n_alleles=0`.
- **Tie-breaking**: When counts tie, lexicographic order determines the top two bases (stable and deterministic).

## Possible extensions

- Add per-marker call rate and heterozygosity %.
- Add a minor-allele frequency (MAF) threshold to ignore rare third alleles.
- Parse slash-delimited genotypes like `A/T` or `C|T`.
- Output per-allele counts to a sidecar table for auditability.

