"""
GC-MS Metabolic Profiling Tool — single-command entrypoint.

Usage
-----
1. Place your Shimadzu QUALITATIVE ANALYSIS REPORT PDFs into:
       data/PDFs/

2. Run:
       python gcms_tool.py

3. On first run the tool inspects the filenames in data/PDFs and writes
   a metadata template at data/samples_template.csv. Open it, verify the
   auto-inferred (condition, day, fraction) values, then save it as
   data/samples.csv.

4. Run again:
       python gcms_tool.py

   The pipeline parses every PDF, re-identifies misclassified peaks,
   builds the abundance matrix, computes diversity statistics, runs PCA,
   renders all figures and produces two Word reports under output/:

      output/GCMS_Metabolic_Profiling_Report.docx   <- full narrative report
      output/Unique_metabolites_table.docx          <- focused unique table

   plus all CSV tables under output/  and high-resolution PNG figures
   under output/figures/.

Optional CLI flags
------------------
  --root PATH                    Project root (default: this script's folder)
  --unique-conditions A,B,C      Comma-separated list of condition labels to
                                 use for the side-by-side unique table.
                                 Default = all conditions in samples.csv.
"""
from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

# Force UTF-8 stdout/stderr on Windows so Unicode column labels never crash
# the print pipeline regardless of the terminal's default encoding.
try:
    if hasattr(sys.stdout, "buffer"):
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                       errors="replace", line_buffering=True)
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8",
                                       errors="replace", line_buffering=True)
except Exception:
    pass

# allow running this script directly without `pip install`
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gcms.pipeline import run_pipeline  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="gcms_tool",
        description="GC-MS metabolic profiling, re-identification and "
                    "diversity analysis from Shimadzu PDFs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument("--root", type=Path,
                   default=Path(__file__).resolve().parent,
                   help="Project root (must contain data/PDFs).")
    p.add_argument("--unique-conditions", type=str, default=None,
                   help="Comma-separated condition labels to include in the "
                        "side-by-side unique-metabolite table.")
    p.add_argument("--offline", action="store_true",
                   help="Skip PubChem enrichment (no internet calls). The "
                        "tool still produces all outputs but without "
                        "InChIKey-based deduplication.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cu = (None if args.unique_conditions is None
          else [c.strip() for c in args.unique_conditions.split(",")
                if c.strip()])
    try:
        run_pipeline(args.root, conditions_for_unique=cu,
                     online=not args.offline)
    except SystemExit as e:
        # gracefully reported above (e.g. samples_template generated)
        print(str(e))
        return 0
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
