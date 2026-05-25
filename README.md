# GC-MS Metabolic Profiling Tool

A self-contained, **fully data-driven** Python tool that turns a folder
of Shimadzu *QUALITATIVE ANALYSIS REPORT* PDFs into a complete
metabolic-profiling deliverable: cleaned identifications, abundance
matrix, alpha-diversity statistics, PCA ordination, compound-class
heat-maps, shared / unique metabolite tables, and two publication-ready
Word reports — for **any** sample matrix, **any** set of conditions,
**any** number of days / fractions and **any** metabolites.

The same single command works for plant-stress experiments, food-
spoilage studies, microbiome work, environmental samples, drug-
metabolism profiles, fermentation broths, etc. — anything that
produces a Shimadzu QGD report in the standard "QUALITATIVE ANALYSIS
REPORT" format.

---

## Table of contents

1. [What you get](#1-what-you-get)
2. [Requirements & installation](#2-requirements--installation)
3. [Quick start (5 minutes)](#3-quick-start-5-minutes)
4. [Project layout](#4-project-layout)
5. [Input files](#5-input-files)
6. [The pipeline, step by step](#6-the-pipeline-step-by-step)
7. [Output files](#7-output-files)
8. [Customisation](#8-customisation)
9. [CLI reference](#9-cli-reference)
10. [Robustness & failure modes](#10-robustness--failure-modes)
11. [Performance & caching](#11-performance--caching)
12. [Troubleshooting](#12-troubleshooting)
13. [Extending the tool](#13-extending-the-tool)
14. [How it works internally](#14-how-it-works-internally)
15. [License & citation](#15-license--citation)

---

## 1. What you get

Every run produces, automatically:

- **Re-identification audit** — biologically implausible NIST library
  picks (drugs, ergot alkaloids, scintillator dyes, halogenated
  synthetics, triglycerides whose retention index is incompatible with
  the oven program) replaced with the best plausible alternative from
  the hit list. Trimethylsilyl / TBDMS suffixes stripped when no
  derivatisation step was performed.
- **PubChem enrichment** — every distinct compound name is queried
  against the free PubChem PUG REST API to obtain canonical IUPAC
  names, CIDs, InChIKeys, molecular formulae and exact masses.
  Stereoisomers / salt forms / tautomers of the same compound are
  merged using the InChIKey-14 connectivity hash.
- **Quality cleanup** — drops invalid / NaN / negative rows, solvent-
  front peaks, low-similarity matches; merges within-sample
  co-elutions of the same compound; normalises spelling variants;
  optionally drops dataset-wide singletons.
- **Metabolite × group abundance matrix** — one row per metabolite,
  one column per `(condition, day)` group; intracellular / extracellular
  (or any other) fractions are merged at the group level.
- **Alpha-diversity statistics** — richness (S), Shannon (H'),
  Simpson (1−D), Pielou evenness (J), total signal per group.
- **Principal-component ordination** — log-transformed, standardised,
  2-component PCA with variance explained shown on the axes.
- **Compound-class diversity** — stacked-bar chart and abundance heat-
  map across ~30 biochemical classes (n-alkanes, fatty alcohols, free
  fatty acids, fatty acid amides, mono/di-glycerides, phenolics,
  diketopiperazines, pyrazines, indoles, sterols, terpenes, etc.).
- **Top-40 metabolite heat-map** — log-scaled abundance of the most
  informative metabolites.
- **Shared / unique tables** — core metabolome (compounds in *all*
  conditions) and the strict per-condition unique sets, with
  abundance values.
- **Two Word reports** — a full narrative report with embedded
  figures + tables, and a focused side-by-side unique-metabolite
  comparison table.
- **Full audit CSVs** — every transformation is reversible and
  documented (`peaks_raw.csv` → `peaks_reidentified.csv` →
  `peaks_cleaned.csv` → `metabolite_abundance_matrix.csv`).

All Word-report text is generated from the data — no hardcoded
biological claims about your samples. If you supply an optional
`data/project.json`, your own narrative paragraphs (sample source,
extraction protocol, condition descriptions) are embedded under the
relevant headings.

---

## 2. Requirements & installation

### System

- **Python 3.10 or newer** (tested on 3.10–3.13).
- Works on **Windows, macOS and Linux**.
- ~200 MB of free disk for the dependencies and output artefacts.
- An internet connection is recommended for PubChem enrichment
  (~30–60 seconds the first time, instant on every re-run thanks to
  the on-disk cache). The tool also runs fully offline with `--offline`.

### Install

```bash
git clone <this-repo>           # or copy the GCMS/ folder
cd GCMS
python -m venv .venv            # optional but recommended
. .venv/Scripts/activate        # Windows
# . .venv/bin/activate          # macOS / Linux
pip install -r requirements.txt
```

### Verify

```bash
python gcms_tool.py --help
```

Should print the CLI help.

---

## 3. Quick start (5 minutes)

```text
1. Drop your Shimadzu QUALITATIVE ANALYSIS REPORT PDFs into:
        data/PDFs/

2. Run:
        python gcms_tool.py

   The first run inspects the filenames and writes a starter
   metadata template at:
        data/samples_template.csv

   It will exit cleanly and ask you to review that file.

3. Open data/samples_template.csv in Excel / a text editor.
   Fill in any blank `condition` / `day` / `fraction` cells and verify
   the auto-inferred ones. Save it as:
        data/samples.csv

4. (Optional) Personalise the Word-report text by editing:
        data/project.json
   (a starter `data/project_template.json` is auto-generated on first
   run; rename it to project.json once filled in).

5. (Optional) Adjust quality thresholds by editing:
        data/quality.json
   (a starter `data/quality_template.json` is auto-generated on first
   run; rename it to quality.json once edited).

6. Run again:
        python gcms_tool.py

   All outputs land under output/ — including two Word reports,
   six 300-DPI PNG figures and ~20 CSV tables.
```

---

## 4. Project layout

```text
GCMS/
├── data/                              # YOUR inputs (you create / edit)
│   ├── PDFs/                          # drop Shimadzu QGD PDFs here
│   ├── samples.csv                    # required (auto-template available)
│   ├── conditions.csv                 # optional — filename -> condition map
│   ├── project.json                   # optional — Word-report personalisation
│   ├── quality.json                   # optional — QC thresholds
│   └── pubchem_cache.json             # auto-created on first online run
│
├── gcms/                              # the package (no editing needed)
│   ├── __init__.py
│   ├── config.py                      # knowledge base (extend for new matrices)
│   ├── parser.py                      # PDF parsing (parallel + cache)
│   ├── reidentify.py                  # candidate-walking re-identification
│   ├── enrich.py                      # PubChem PUG-REST enrichment
│   ├── quality.py                     # error filters + dedup + co-elution merge
│   ├── metadata.py                    # filename heuristics + samples.csv loader
│   ├── project_info.py                # optional user narrative loader
│   ├── profile.py                     # abundance matrix + diversity stats
│   ├── figures.py                     # all six figures
│   ├── report.py                      # two Word-document writers
│   └── pipeline.py                    # the orchestrator
│
├── output/                            # GENERATED automatically (recreated each run)
│   ├── peaks_raw.csv                  # parser cache
│   ├── hits_raw.csv                   # parser cache
│   ├── peaks_reidentified.csv         # re-identification audit trail
│   ├── reidentification_corrections.csv
│   ├── pubchem_enrichment.csv         # NIST name -> PubChem CID/IUPAC/InChIKey/...
│   ├── pubchem_name_remap.csv         # NIST -> PubChem canonical name renames
│   ├── synonym_collapses.csv          # NIST -> normalised name collapses
│   ├── quality_report.csv             # per-step row counts of the QC pipeline
│   ├── peaks_cleaned.csv              # final cleaned peak table
│   ├── metabolite_abundance_matrix.csv
│   ├── metabolite_presence.csv
│   ├── diversity_indices.csv
│   ├── class_abundance.csv
│   ├── class_richness.csv
│   ├── set_core_metabolites.csv
│   ├── set_unique_to_<condition>.csv  # one per condition
│   ├── set_all_summary.csv
│   ├── lists/
│   │   ├── top40_most_abundant.csv
│   │   ├── top5_per_group.csv
│   │   └── unique_<condition>.csv     # one per condition
│   ├── figures/                       # six 300-DPI PNGs
│   │   ├── 01_PCA.png
│   │   ├── 02_class_diversity_bar.png
│   │   ├── 03_alpha_diversity.png
│   │   ├── 04_heatmap_top40_metabolites.png
│   │   ├── 05_heatmap_classes.png
│   │   └── 06_unique_vs_core.png
│   ├── GCMS_Metabolic_Profiling_Report.docx
│   └── Unique_metabolites_table.docx
│
├── gcms_tool.py                       # single command-line entrypoint
├── requirements.txt
├── README.md                          # this file
└── .gitignore
```

---

## 5. Input files

### 5.1  `data/PDFs/*.pdf` (required)

Place one Shimadzu *QUALITATIVE ANALYSIS REPORT* PDF per GC-MS
injection. Filenames become sample identifiers — keep them descriptive
(for example `MyOrg-Day1-Polar.pdf` becomes sample `MyOrg-Day1-Polar`).

### 5.2  `data/samples.csv` (required after first run)

A plain CSV mapping each sample to its experimental design. Required
columns:

| Column      | Type    | Description |
|-------------|---------|-------------|
| `sample`    | string  | Must equal the PDF filename without the `.pdf` extension |
| `condition` | string  | Free-form label (e.g. `Control`, `Drought`, `Cd`) |
| `day`       | integer | Incubation / time-point in days |
| `fraction`  | string  | Free-form label (e.g. `Intra`, `Extra`, `Polar`, `Apolar`); samples with the same `(condition, day)` and different `fraction` are merged by area-% summation |

Example:

```csv
sample,condition,day,fraction
MyOrg-1-Polar,Control,1,Polar
MyOrg-1-Apolar,Control,1,Apolar
MyOrg-15-Polar,DroughtStress,15,Polar
MyOrg-15-Apolar,DroughtStress,15,Apolar
```

The first run auto-generates `data/samples_template.csv` from
filename heuristics; review it, fix anything wrong, save as
`samples.csv`.

### 5.3  `data/conditions.csv` (optional)

Helps the auto-template detect your `condition` and `day` from
filenames. Useful if you have a regular naming pattern. Three columns:

```csv
token,condition,default_day
AS15,As,15
CDR10,Cd,10
G,Control,1
```

If a filename starts with `<token>` (case-insensitive), that condition
is assigned. The day is taken from the filename if a number is present;
otherwise `default_day` is used. Without this file, the auto-template
just leaves `condition` blank for unknown filenames.

### 5.4  `data/project.json` (optional)

Personalises the Word report with your own narrative paragraphs. A
starter `project_template.json` is auto-created on first run.
Schema:

```json
{
  "title":         "<report title>",
  "subtitle":      "<one-line subtitle>",
  "sample_source": "<one paragraph describing where the samples come from>",
  "extraction":    "<one paragraph describing the extraction protocol>",
  "fraction_note": "<comment on how fractions are merged>",
  "conditions": {
    "Control":       "<one-line description>",
    "DroughtStress": "<one-line description>"
  },
  "notes":         "<free-form notes appended at the end>"
}
```

If `project.json` is absent, the report uses neutral generic text.

### 5.5  `data/quality.json` (optional)

Tunes the quality-control thresholds. A starter
`quality_template.json` is auto-created on first run.

| Key | Type | Default | Effect |
|-----|------|---------|--------|
| `min_rt`             | float | `2.0`  | Drop peaks with retention time below this (in min) — removes solvent-front peaks |
| `si_min`             | int   | `50`   | Drop peaks whose final similarity index is below this (0–100). Below ~50 the library match is too weak to trust |
| `trace_min`          | float | `0.0`  | Drop peaks with `area_pct` below this. Default 0 = off; set to e.g. `0.05` to remove integration-noise peaks |
| `drop_singletons`    | bool  | `false`| If true, drop compounds detected in only one sample with very low abundance |
| `singleton_max_area` | float | `0.05` | Threshold used by `drop_singletons` |
| `merge_coelutions`   | bool  | `true` | Merge multiple peaks of the same compound within one sample |

---

## 6. The pipeline, step by step

```
┌─────────────────────────────────────────────────────────────────────────────┐
│   data/PDFs/*.pdf                                                           │
└─────────────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
   [1/5] PDF parsing  ───────────────────────────────────────────► peaks_raw.csv
         pdfplumber, parallel processing, on-disk cache              hits_raw.csv
                              │
                              ▼
   [2/5] Re-identification  ─────────────────────────────────► peaks_reidentified.csv
         candidate-walking engine                              reidentification_corrections.csv
                              │
                              ▼
   [2a/5] PubChem enrichment  ──────────────────────────────► pubchem_enrichment.csv
          name → CID → IUPACName / InChIKey / formula / mass   pubchem_name_remap.csv
                              │
                              ▼
   [2b/5] Quality cleanup  ─────────────────────────────────► peaks_cleaned.csv
          invalid rows / solvent-front / low SI / trace        quality_report.csv
          synonym normalisation                                synonym_collapses.csv
          InChIKey-based co-elution merging
                              │
                              ▼
   [3/5] Profiling  ─────────────────────────────────────► metabolite_abundance_matrix.csv
         abundance / presence / class / set tables           diversity_indices.csv
                                                             class_abundance.csv
                                                             class_richness.csv
                                                             metabolite_presence.csv
                                                             set_core_metabolites.csv
                                                             set_unique_to_<cond>.csv
                                                             lists/top40_most_abundant.csv
                                                             lists/top5_per_group.csv
                                                             lists/unique_<cond>.csv
                              │
                              ▼
   [4/5] Figures  ─────────────────────────────────────────► figures/01_PCA.png
                                                             figures/02_class_diversity_bar.png
                                                             figures/03_alpha_diversity.png
                                                             figures/04_heatmap_top40_metabolites.png
                                                             figures/05_heatmap_classes.png
                                                             figures/06_unique_vs_core.png
                              │
                              ▼
   [5/5] Word reports ─────────────────────────────────────► GCMS_Metabolic_Profiling_Report.docx
                                                             Unique_metabolites_table.docx
```

### 6.1  Step 1 — PDF parsing

`gcms/parser.py` opens every `*.pdf` in `data/PDFs/`, extracts the
peak table (Peak#, R.Time, area, area%, height, height%, A/H, mark,
system Hit-#1 name) AND the library appendix (every Hit#:Entry:
Library: SI:formula:CAS:MolWeight:RetIndex:CompName for up to five
alternative hits per peak). PDFs are processed in parallel
(up to 8 worker processes). Results are cached as
`output/peaks_raw.csv` and `output/hits_raw.csv`; subsequent runs
skip parsing if no PDF is newer than the cache.

### 6.2  Step 2 — Re-identification (candidate walking)

`gcms/reidentify.py` walks the alternative library hits per peak and
selects the first plausible candidate using lists in `gcms/config.py`:

- **Implausible substrings** — names that are essentially never
  produced biologically: pharmaceuticals (Floxuridine, Debrisoquine,
  Ergotamine, Bromocriptine), scintillator dyes, halogenated
  synthetics, triglycerides whose NIST retention index is impossible
  in this oven program, etc. These are skipped.
- **Contaminant substrings** — lab and instrument artefacts: cyclic
  siloxanes (column bleed), phthalates and other plasticisers,
  perfluorinated derivatising-agent residues. These are tagged as
  contaminants and excluded from biological diversity statistics.
- **TMS / TBDMS suffix stripping** — when the protocol does not
  include silylation, library entries annotated as
  `, TMS derivative` are converted to their underlying compound.

The original system pick, the chosen final pick, the SI scores, the
chosen Hit number and a textual reason are preserved for every peak.

### 6.3  Step 2a — PubChem enrichment

`gcms/enrich.py` sends every distinct re-identified name to PubChem's
free PUG REST API:

- `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<NAME>/cids/JSON`
  → CID
- `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/<CID>/property/IUPACName,InChIKey,MolecularFormula,ExactMass/JSON`
  → properties

Returned values:

- `cid` — PubChem identifier (used for hyperlinks if you want them)
- `canonical_name` — the IUPAC name (replaces the NIST notation)
- `inchikey` — 27-character compound hash
- `inchikey14` — first 14 characters; identical for stereoisomers,
  salt forms and tautomers, so this is the *canonical compound key*
- `formula`, `exact_mass` — cross-checks against NIST values

Lookups are cached in `data/pubchem_cache.json`. Network failures
(rate limit, no internet) degrade gracefully — any compound that
cannot be matched simply keeps its original NIST name and a
lower-cased name is used as the deduplication key. Use `--offline`
to skip enrichment entirely.

### 6.4  Step 2b — Quality cleanup and deduplication

`gcms/quality.py` applies, in order:

1. Drop rows with missing or non-positive `area_pct`, `rt`, `area`
   or `height`.
2. Drop solvent-front peaks below `min_rt` (default 2.0 min).
3. Drop low-similarity peaks below `si_min` (default 50).
4. Optionally drop trace peaks below `trace_min`.
5. Normalise compound-name spelling and stereochemistry on the few
   compounds PubChem could not match (`Oleic Acid, (Z)-` →
   `Oleic acid`, `n-Hexadecanoic acid` → `Hexadecanoic acid`).
6. Merge within-sample co-elutions on `compound_key` — when the same
   compound appears in multiple peaks of the same sample (broadened
   peaks, integration splits, stereoisomers detected separately), the
   rows are summed: total `area`, sum `area_pct`, max `final_si`,
   median `rt`.
7. Optionally drop dataset-wide singletons.

Per-step counts go to `output/quality_report.csv`.

### 6.5  Step 3 — Profiling

`gcms/profile.py` builds:

- `metabolite_abundance_matrix.csv` — metabolites (rows) × groups
  (columns), values are sum-of-area-% across every fraction in that
  group. Contaminants and unidentified compounds are excluded.
- `metabolite_presence.csv` — same shape, presence/absence + a
  `#groups` column counting how many groups each metabolite is in.
- `diversity_indices.csv` — per group: richness `S`, Shannon `H'`,
  Gini–Simpson `1−D`, Pielou evenness `J`, total signal.
- `class_abundance.csv` and `class_richness.csv` — by compound class
  (n-alkane, fatty alcohol, free fatty acid, fatty acid amide,
  mono/di-glyceride, diketopiperazine, etc.).
- `set_core_metabolites.csv` — compounds present in *all* conditions.
- `set_unique_to_<condition>.csv` — compounds detected in only one
  condition (across all its days and fractions).
- `lists/top40_most_abundant.csv` — top 40 by maximum across-group
  abundance, with `argmax_group` and `#groups_present`.
- `lists/top5_per_group.csv` — top 5 metabolites within each group.

### 6.6  Step 4 — Figures

`gcms/figures.py` renders six 300-DPI PNGs:

| File | Plot | Description |
|------|------|-------------|
| `01_PCA.png` | PCA | Log-transformed standardised metabolite × group matrix; PC1 vs PC2 with variance explained |
| `02_class_diversity_bar.png` | Stacked bar | Distinct metabolites per compound class, per group |
| `03_alpha_diversity.png` | Bar triplet | Richness, Shannon and Simpson per group |
| `04_heatmap_top40_metabolites.png` | Heat-map | Top-40 metabolites by max across-group abundance, log10 colour scale |
| `05_heatmap_classes.png` | Heat-map | Compound-class abundance per group, log10 colour scale, raw annotations |
| `06_unique_vs_core.png` | Bar | Number of metabolites unique to each condition + the shared core |

### 6.7  Step 5 — Word reports

`gcms/report.py` writes two `.docx` files:

- `GCMS_Metabolic_Profiling_Report.docx` — full narrative report with
  ~10 chapters, all six figures embedded, ~10 tables, four annexes
  (per-condition uniques, top-40, top-5-per-group, files index).
  Every sentence is built from the actual numbers in the data.
- `Unique_metabolites_table.docx` — focused, side-by-side comparison
  of the per-condition unique metabolites with abundance values.

---

## 7. Output files

### CSV tables under `output/`

| File | Rows | Columns | Description |
|------|------|---------|-------------|
| `peaks_raw.csv` | 1 per peak | sample, peak#, RT, area, area%, ..., system_name | Direct parse of the PDF peak table |
| `hits_raw.csv` | up to 5 per peak | sample, line, hit#, entry, library, SI, formula, CAS, MolWeight, RetIndex, CompName | All NIST library hits per peak |
| `peaks_reidentified.csv` | 1 per peak | + system_si, final_name, final_si, hit_chosen, reidentification, is_*  | Audit trail of the re-identification |
| `reidentification_corrections.csv` | 1 per (system_name, final_name) pair | + occurrences, mean_area_pct | Aggregated correction table |
| `pubchem_enrichment.csv` | 1 per distinct compound name | nist_name, cid, canonical_name, inchikey, inchikey14, formula, exact_mass, source | PubChem lookup results |
| `pubchem_name_remap.csv` | 1 per remap | nist_name, canonical_name | Pre→post-PubChem name renames |
| `synonym_collapses.csv` | 1 per collapse | original_name, normalised_name | Final-name normalisation |
| `quality_report.csv` | 1 per QC step | metric, count | Counts at each cleanup step |
| `peaks_cleaned.csv` | 1 per cleaned row | (same as `peaks_reidentified.csv` + canonical_name, inchikey14, …) | Final cleaned peak table |
| `metabolite_abundance_matrix.csv` | 1 per metabolite | (one column per group) | Sum-of-area% matrix |
| `metabolite_presence.csv` | 1 per metabolite | (one column per group) + #groups + groups_present + mean/max area_pct | Presence / absence + summary |
| `diversity_indices.csv` | 1 per group | richness, shannon, simpson, evenness, total_signal | Alpha-diversity indices |
| `class_abundance.csv` | 1 per class | (one column per group) | Class abundance |
| `class_richness.csv` | 1 per class | (one column per group) | Distinct compounds per class |
| `set_core_metabolites.csv` | 1 per metabolite | (one column per condition) + #conditions | Core compounds |
| `set_unique_to_<cond>.csv` | 1 per metabolite | (one column per condition) | Per-condition unique compounds |
| `lists/top40_most_abundant.csv` | 40 | max_area_pct, mean_area_pct, #groups_present, argmax_group | Top-40 by maximum |
| `lists/top5_per_group.csv` | 5 × #groups | group, rank, metabolite, area_pct | Top-5 per group |
| `lists/unique_<cond>.csv` | 1 per metabolite | abundance_area_pct | Per-condition unique table |

### Figures under `output/figures/`

Six 300-DPI PNGs (see §6.6).

### Word reports under `output/`

Two `.docx` files (see §6.7).

---

## 8. Customisation

The tool is designed so most users never need to edit Python. Three
optional input files cover the most common customisations:

| File | What it does | When to use |
|------|--------------|-------------|
| `data/conditions.csv`     | Filename → (condition, default day) mapping | Your filenames follow a regular pattern (e.g. `AS15-Intra`) |
| `data/quality.json`       | QC thresholds (min RT, min SI, trace cut, …) | Different solvent cut, different SI tolerance, want to remove trace peaks |
| `data/project.json`       | Personalises Word-report narrative text | You want the report to mention your organism / extraction protocol |

If you need to extend the *knowledge base* (new contaminant patterns,
new biochemical class), edit `gcms/config.py` directly — see §13.

---

## 9. CLI reference

```text
python gcms_tool.py [--root PATH] [--unique-conditions A,B,C] [--offline]

  --root PATH                Project root containing data/ and output/.
                             Default: the folder containing gcms_tool.py.

  --unique-conditions A,B,C  Restrict the side-by-side unique-metabolite
                             table to specific conditions (e.g. heavy-
                             metal stresses only, excluding a substrate
                             group). Default: all conditions in
                             samples.csv.

  --offline                  Skip PubChem enrichment entirely. The tool
                             still produces every output but without
                             InChIKey-based deduplication, and uses
                             NIST names verbatim.

  -h, --help                 Show this help.
```

Examples:

```bash
# Standard run
python gcms_tool.py

# Run on a different project folder
python gcms_tool.py --root D:\Studies\Plant_Drought_2025

# Compare only stress conditions (drop the substrate-control group)
python gcms_tool.py --unique-conditions Control,Drought,Heat

# Fully offline (no internet access)
python gcms_tool.py --offline
```

---

## 10. Robustness & failure modes

The pipeline is designed to do something sensible in every realistic
edge case:

| Scenario | Behaviour |
|----------|-----------|
| Missing `samples.csv` | Auto-writes `samples_template.csv` and asks you to review |
| Filenames the heuristics can't parse | Template row gets blank `condition`; tool refuses to proceed until filled in |
| Extra PDF not in `samples.csv` | Warns and skips |
| `samples.csv` row with no matching PDF | Warns and ignores |
| Corrupted / unreadable PDF | Logged as `FAILED <name>: <error>`; pipeline continues with the rest |
| Different number of conditions / days / fractions | Group order, palette, PCA, heat-maps, diversity tables and Word text all adapt automatically |
| Single-condition / single-group data | Figures requiring ≥2 groups are skipped silently |
| Word file already open in Word | Tool writes to `…__new.docx` instead of crashing |
| No internet / PubChem rate limit | Auto-retries 503/429 up to 3 times with exponential backoff; remaining failures keep NIST names |
| No internet for the entire run | Use `--offline`; tool still produces all artefacts |
| Very large dataset (>1000 PDFs) | Parallel parser keeps memory bounded; cache means re-runs are fast |
| Re-run after editing `samples.csv` only | PDF parser cache reused; total run time becomes seconds |
| Unicode in compound names | UTF-8 stdout enforced; Word docs handle Unicode natively |

---

## 11. Performance & caching

Three levels of cache make re-runs fast:

1. **PDF parser cache** — `output/peaks_raw.csv` + `output/hits_raw.csv`.
   Reused unless a PDF in `data/PDFs/` is newer than these files.
   Cold parse of 18 PDFs takes ~5 minutes on one core; parallelised
   it takes ~1 minute. Re-run with cache: <2 seconds.

2. **PubChem enrichment cache** — `data/pubchem_cache.json`.
   Reused for every compound the tool has seen before. Cold cache:
   ~30–60 seconds for 200 compounds (5 RPS rate limit). Warm cache:
   instant.

3. **Quality / profiling / figures** — recomputed every run, but the
   underlying tables are small (~1000 peaks → ~200 metabolites)
   so the entire downstream pipeline finishes in <10 seconds.

---

## 12. Troubleshooting

| Problem | Cause / fix |
|---------|-------------|
| `No PDF files found in data/PDFs` | Drop your PDFs into `data/PDFs/` (not anywhere else). |
| `samples.csv is missing required columns` | Required columns are exactly `sample,condition,day,fraction`. Edit and re-save. |
| `samples.csv has N row(s) with non-integer 'day'` | The `day` column must be an integer (1, 7, 15 …). Fix any text values. |
| `samples.csv has N row(s) with blank 'condition'` | Fill in every `condition` cell — the auto-template often leaves them blank when it can't infer a token. |
| Word report shows blank "Conditions" subsection | Either no `data/project.json`, or none of its `conditions:{}` entries match the labels in `samples.csv`. Match the labels exactly. |
| PubChem step is slow on first run | This is a one-time cost. The cache makes every re-run instant. Use `--offline` to skip it entirely. |
| Some compounds show `source: no_match` | Either the NIST name is too obscure for PubChem, or the API was busy. Re-run later — failed lookups are not cached, so the tool will try them again. |
| `PermissionError` writing the `.docx` | The Word file is currently open in Microsoft Word. Tool falls back to `…__new.docx`. Close Word and re-run. |
| PCA figure is empty / one-dimensional | Fewer than 2 groups in the data. PCA needs ≥2 groups. |

---

## 13. Extending the tool

All knowledge bases are in `gcms/config.py` and are pure Python lists
/ regex tuples — extend them in place to handle new sample matrices.

| Knowledge base | What it does | When to extend |
|----------------|--------------|----------------|
| `IMPLAUSIBLE_SUBSTRINGS` | Names rejected from the candidate hit walk | A new biology-implausible NIST entry pollutes your data |
| `CONTAMINANT_SUBSTRINGS` | Names tagged as contaminants and excluded from diversity | New column-bleed pattern, new plastic additive, etc. |
| `CLASS_RULES` | Regex rules that bin compounds into biochemical classes | Your matrix has compound classes not yet covered |
| `SYNONYM_MAP` | Manual NIST → canonical-name remaps | A specific NIST notation you want collapsed |

The pipeline, parser and report are orthogonal to these lists, so
you can ship a customised `config.py` to colleagues without modifying
any other file.

To add a new figure, write a function in `gcms/figures.py` and call it
from `gcms/pipeline.py`. To add a new section to the Word report,
write a helper in `gcms/report.py`.

---

## 14. How it works internally

### Module dependencies

```
gcms_tool.py  ──► gcms.pipeline.run_pipeline
                        │
                        ├──► gcms.metadata           (samples.csv, conditions.csv)
                        ├──► gcms.parser             (PDF -> peaks/hits)
                        ├──► gcms.reidentify         (candidate-walking + classifier)
                        │     └──► gcms.config       (knowledge base)
                        ├──► gcms.enrich             (PubChem PUG REST)
                        ├──► gcms.quality            (filters + dedup + co-elution merge)
                        ├──► gcms.profile            (matrix + diversity + sets)
                        ├──► gcms.figures            (matplotlib + seaborn)
                        ├──► gcms.project_info       (optional user narrative)
                        └──► gcms.report             (python-docx writers)
```

### Data flow

```
peaks_raw  +  hits_raw            (Step 1, parser cache)
        │
        ▼
peaks_reidentified                 (Step 2, with audit trail)
        │
        ▼
peaks_reidentified  +  pubchem_enrichment    (Step 2a, joined by name)
        │
        ▼
peaks_cleaned                       (Step 2b, deduped + filtered)
        │
        ▼
metabolite × group abundance        (Step 3)
        │
        ▼
diversity / class / set tables      (Step 3, all derived)
        │
        ▼
figures (PCA, bars, heat-maps)      (Step 4)
        │
        ▼
Word reports                        (Step 5)
```

### Key design choices

- **Plain CSVs at every step** — every transformation is reversible
  and inspectable in Excel. No proprietary formats, no binary blobs.
- **Knowledge separated from logic** — `gcms/config.py` is the only
  file you need to edit to handle a new matrix.
- **Generic narrative** — `gcms/report.py` builds every sentence from
  the data; biology-specific language only appears if you supply
  `data/project.json`.
- **Optional online enrichment** — PubChem is a strong booster, but
  the tool works fully offline.
- **Three-level caching** — parser, PubChem and quality results all
  cache on disk so iteration is fast.

---

## 15. License & citation

This tool was developed for academic / research use. PubChem data is
public-domain under the U.S. Government Open Data policy. NIST library
data belongs to NIST and your local GCMSsolution license.

If you use the tool in a publication, please cite:

- **PubChem PUG REST** —
  Kim S. *et al.* PubChem 2023 update. *Nucleic Acids Research* 51(D1),
  D1373-D1380, 2023.

(A formal citation entry for this tool will be added once a permanent
identifier is issued.)

---

*Generated artefacts (CSVs, figures, Word reports) are reproduced on
every run from the contents of `data/`. Delete `output/` any time to
force a fresh build.*
