# marker-allele-infer

Infer the **allele pair** per SNP/marker from a genotype matrix and export a tidy CSV with QC flags.

- Input: CSV with rows = markers, columns = samples (genotypes).
- Output: Adds three columns up front:
  - `alleles` → e.g., `A/T` or `C/C`
  - `n_alleles` → count of distinct bases seen across samples for that marker
  - `multi_allelic` → `YES` if >2 distinct bases were seen

> Designed for quick QA/QC and for downstream tools (e.g., PLINK/GWAS/Shiny) that expect biallelic markers.

---

## Features

- Prefers **homozygote evidence** to determine the allele pair; falls back to all calls when needed.
- Understands **IUPAC single-letter** codes (R, Y, S, W, K, M, B, D, H, V, N).
- Handles common **missing** tokens: `""`, `.`, `-`, `NA`, `N`, `NN`, `HH`.
- Flags **tri-/multi-allelic** markers and optionally **drops** them (`--drop-multi`).
- Cleans **sample headers** by trimming text after `;` (disable with `--no-clean-headers`).
- Auto-creates a **Marker** ID column if the input has none.

## Install

```bash
git clone https://github.com/<YOUR-USERNAME>/marker-allele-infer.git
cd marker-allele-infer
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Quickstart

```bash
python src/infer_marker_alleles.py --in examples/example_genotypes.csv --out marker_with_alleles.csv
```

With tri-allelic removal:

```bash
python src/infer_marker_alleles.py --in examples/example_genotypes.csv --out marker_with_alleles_biallelic.csv --drop-multi
```

## Input expectations

- Genotype cells can be:
  - Two-letter genotypes: `AA`, `AT`, `TT`, etc. (explicitly mark homozygotes like `AA`, `CC`).
  - IUPAC single-letter codes: `R (A/G)`, `Y (C/T)`, `S (G/C)`, `W (A/T)`, `K (G/T)`, `M (A/C)`, `B (C/G/T)`, `D (A/G/T)`, `H (A/C/T)`, `V (A/C/G)`, `N (missing)`.
  - Missing: `""`, `.`, `-`, `NA`, `N`, `NN`, `HH`.
- **Note:** Single-letter bases `A/C/G/T` are treated as **not explicitly homozygous** by default (conservative). If your data always uses single letters for homozygotes, see customization below.

## Output columns

- `Marker` (added if missing): unique ID if none provided.
- `alleles`: the inferred allele pair (`A/T` or `C/C` for monomorphic).
- `n_alleles`: number of distinct bases observed across samples (0, 1, 2, or >2).
- `multi_allelic`: `YES` if `n_alleles > 2`, else blank.
- Followed by your original genotype columns.

## CLI options

```
--in            Path to input CSV (required)
--out           Path to output CSV (default: marker_with_alleles.csv)
--marker-col    Name of existing marker ID column (optional)
--drop-multi    Drop markers with >2 alleles
--no-clean-headers   Do NOT trim sample headers at semicolons
```

## Customization

Common tweaks (edit `src/infer_marker_alleles.py`):

- **Count single-letter A/C/G/T as homozygotes:** in `parse_cell`, when `len(s)==1 and s in IUPAC`, you may special-case `s in {'A','C','G','T'}` to return `( {s}, True )`.
- **Ignore rare third alleles:** apply a frequency threshold before deciding `n_alleles` and top two alleles (e.g., ignore alleles with count < 2).
- **Accept slash genotypes (`A/T`):** add a branch in `parse_cell` to split on `/` and handle length-3 strings like `A/T`.

## Example data

See `examples/example_genotypes.csv` for a toy matrix and try the quickstart command above.

## Why this tool

We often receive genotype matrices without an explicit "allele" descriptor per marker. This script inspects calls across samples to infer the allele pair (favoring homozygotes) and flags multi-allelic sites—making the table cleaner for QA/QC and biallelic-only pipelines.

## License

MIT — see [LICENSE](LICENSE).

## Citation / Acknowledgment

If you use this in a paper or report, you can cite the repo URL and mention: *“marker-allele-infer: a small utility to infer per-marker allele pairs and flag multi-allelic sites from genotype matrices.”*
