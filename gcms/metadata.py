"""
Sample metadata management.

Handles:
  - Auto-generation of `data/samples_template.csv` from filename heuristics
    when `data/samples.csv` is missing.
  - Loading and validating an existing `data/samples.csv`.
  - Cross-checking metadata against the PDFs actually present (warns
    about missing or extra files).

Filename heuristics
-------------------
A condition prefix is recognised in the filename stem when the stem
starts with a token (case-insensitive). The token-to-condition mapping
is read from `data/conditions.csv` if present. If that file is absent,
the tool simply leaves the `condition` column blank in the auto-
generated template and the user fills it in manually.

`data/conditions.csv` has three columns:

    token,condition,default_day
    AS15,As,15
    CDR10,Cd,10
    G,Control,1

Tokens are checked in the order they appear in the file (longest /
most-specific first is recommended). Day inference still tries to
extract a number from the filename first; `default_day` is only used
when no number can be inferred.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

# Each entry is (filename token, default condition label, default day).
# The token is matched at the START of the filename stem (case-insensitive).
# `default_day` is used only if no number can be inferred from the stem.
# Order matters — list more specific tokens first.
#
# This is the LAST-RESORT fallback used only when neither
# `data/conditions.csv` nor an existing `data/samples.csv` is present.
# It is intentionally empty so the tool stays generic. Users supply
# their own token map via `data/conditions.csv` or just edit the
# auto-generated `samples_template.csv` directly.
DEFAULT_CONDITION_TOKENS: list[tuple[str, str, int | None]] = []


def load_condition_tokens(data_dir: Path) -> list[tuple[str, str, int | None]]:
    """Load filename heuristics from data/conditions.csv if present."""
    cf = data_dir / "conditions.csv"
    if not cf.exists():
        return list(DEFAULT_CONDITION_TOKENS)
    df = pd.read_csv(cf)
    required = {"token", "condition"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"{cf} must have at least columns {sorted(required)} "
            f"(plus an optional 'default_day' column).")
    tokens: list[tuple[str, str, int | None]] = []
    for _, r in df.iterrows():
        tok = str(r["token"]).strip()
        cond = str(r["condition"]).strip()
        day = None
        if "default_day" in df.columns:
            v = r["default_day"]
            if pd.notna(v) and str(v).strip() != "":
                try:
                    day = int(v)
                except Exception:
                    day = None
        if tok and cond:
            tokens.append((tok, cond, day))
    return tokens


def write_conditions_template(data_dir: Path) -> Path:
    """Write a template conditions.csv if absent, with a few empty rows."""
    cf = data_dir / "conditions.csv"
    if cf.exists():
        return cf
    # purely illustrative — user fills in their own
    pd.DataFrame([
        {"token": "", "condition": "", "default_day": ""},
    ]).to_csv(cf, index=False)
    return cf

INTRA_RE  = re.compile(r"intra", re.I)
EXTRA_RE  = re.compile(r"extra", re.I)
DAY_RE    = re.compile(r"(\d+)\s*[Dd]")              # "3D", "10d"
TRAIL_NUM = re.compile(r"(\d+)\s*$")                 # trailing number


def _infer_one(stem: str,
               tokens: list[tuple[str, str, int | None]]) -> dict:
    """Return a metadata dict inferred from one filename stem.
    Missing fields are returned as empty strings so the user can fill
    them in."""
    s = stem.strip()
    fraction = ("Intra" if INTRA_RE.search(s)
                else "Extra" if EXTRA_RE.search(s) else "")

    # strip the trailing intra/extra suffix to get the base condition+day
    base = re.sub(r"[-_\s]?(intra|extra)\s*$", "", s, flags=re.I).strip(" -_")

    cond = ""
    day: int | None = None

    for tok, lbl, default_day in tokens:
        if base.upper().startswith(tok.upper()):
            cond = lbl
            tail = base[len(tok):].strip(" -_")
            m = DAY_RE.search(tail) or DAY_RE.search(base)
            if m:
                day = int(m.group(1))
            else:
                m2 = TRAIL_NUM.search(tail)
                if m2:
                    day = int(m2.group(1))
                elif default_day is not None:
                    day = default_day
            break

    if cond == "" and day is None:
        m = re.search(r"(\d+)", base)
        if m:
            day = int(m.group(1))

    return dict(
        sample=stem,
        condition=cond,
        day=("" if day is None else day),
        fraction=fraction,
    )


def infer_metadata(pdf_files: list[Path],
                    tokens: list[tuple[str, str, int | None]] | None = None
                    ) -> pd.DataFrame:
    if tokens is None:
        tokens = list(DEFAULT_CONDITION_TOKENS)
    rows = [_infer_one(f.stem, tokens) for f in pdf_files]
    df = pd.DataFrame(rows, columns=["sample", "condition", "day", "fraction"])
    return df.sort_values("sample").reset_index(drop=True)


def write_template(pdf_files: list[Path], out_csv: Path,
                    tokens: list[tuple[str, str, int | None]] | None = None
                    ) -> pd.DataFrame:
    """Write the auto-inferred metadata template to disk."""
    df = infer_metadata(pdf_files, tokens=tokens)
    df.to_csv(out_csv, index=False)
    return df


def load_metadata(samples_csv: Path) -> pd.DataFrame:
    """Load and validate a user-supplied samples.csv.

    Required columns: sample, condition, day, fraction.
    `day` must be an integer. `fraction` defaults to 'Unknown' if blank.
    """
    df = pd.read_csv(samples_csv)
    required = {"sample", "condition", "day", "fraction"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"samples.csv is missing required columns: {sorted(missing)}. "
            f"Required columns are exactly: {sorted(required)}.")

    # clean whitespace
    for col in ("sample", "condition", "fraction"):
        df[col] = df[col].astype(str).str.strip()
    df["fraction"] = df["fraction"].replace({"": "Unknown", "nan": "Unknown"})

    # coerce day to int with helpful error
    bad = df[~df["day"].apply(lambda x: pd.notna(x)
                                          and str(x).strip().lstrip("-").isdigit())]
    if not bad.empty:
        raise ValueError(
            f"samples.csv has {len(bad)} row(s) with non-integer 'day' "
            f"values. Offending rows:\n{bad}")
    df["day"] = df["day"].astype(int)

    # warn / fail on missing condition
    blank = df[df["condition"] == ""]
    if not blank.empty:
        raise ValueError(
            f"samples.csv has {len(blank)} row(s) with blank 'condition'. "
            f"Please fill in the condition label for every sample.\n"
            f"Offending rows:\n{blank[['sample','condition','day']]}")

    # warn on duplicates
    dups = df[df.duplicated("sample", keep=False)]
    if not dups.empty:
        raise ValueError(
            f"samples.csv has duplicate 'sample' entries:\n"
            f"{dups.sort_values('sample')}")

    df["group"] = df["condition"] + "_d" + df["day"].astype(str)
    return df


def cross_check(meta: pd.DataFrame, pdf_files: list[Path]) -> tuple[list[str], list[str]]:
    """Return (missing, extra) — samples in metadata without a PDF, and
    PDFs without a metadata row."""
    pdf_stems = {f.stem for f in pdf_files}
    meta_stems = set(meta["sample"])
    missing = sorted(meta_stems - pdf_stems)
    extra = sorted(pdf_stems - meta_stems)
    return missing, extra
