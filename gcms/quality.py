"""
Data quality, deduplication and error filtering.

Applied AFTER re-identification, BEFORE the abundance matrix is built.

Steps (in order):

1. drop_invalid_rows
   - drop rows with area_pct ≤ 0 or NaN
   - drop rows with retention time, height or area NaN

2. drop_solvent_front
   - drop rows with rt < min_rt (default 2.0 min) — the solvent cut was
     at 1.5 min in the supplied method, so anything before this is
     residual solvent / matrix front. Configurable via QualityConfig.

3. drop_low_si
   - drop rows with final_si < si_min (default 50). The Shimadzu /
     NIST similarity index is a 0–100 score; below ~50 the assignment
     is little better than random. Configurable.

4. drop_trace_peaks
   - drop rows with area_pct < trace_min (default 0.0; off by default
     because trace peaks can be biologically real markers).

5. normalize_compound_names
   - collapse common spelling / stereochemistry variants
       'Oleic Acid, (Z)-'  -> 'Oleic acid'
       'n-Hexadecanoic acid' -> 'Hexadecanoic acid'
       'Palmitic Acid'     -> 'Palmitic acid'
   - case-fold the first letter of the molecule's name only when an
     all-lowercase variant already exists.

6. merge_within_sample_coelutions
   - if the same compound is reported at multiple peaks in the same
     sample (co-elution), KEEP one row per (sample, compound) and
     SUM the area / area_pct values, taking the median rt and the
     maximum SI. This fixes the "Heneicosane appears 9 times in
     T20-3D-Extra" artefact while preserving the total signal.

7. drop_dataset_singletons (optional, default off)
   - drop compounds that appear in only one sample with very low
     abundance (likely noise).

A QualityReport object captures how many rows were filtered at each
step so the report can document them.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class QualityConfig:
    min_rt: float = 2.0          # min retention time (min)
    si_min: int = 50             # min similarity index (0-100)
    trace_min: float = 0.0       # min area_pct threshold
    drop_singletons: bool = False
    singleton_max_area: float = 0.05
    merge_coelutions: bool = True


@dataclass
class QualityReport:
    n_in: int = 0
    n_out: int = 0
    invalid_rows: int = 0
    solvent_front: int = 0
    low_si: int = 0
    trace: int = 0
    coelutions_merged: int = 0
    singletons: int = 0
    synonym_collapses: dict[str, str] = field(default_factory=dict)

    def as_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame([
            ("rows_in",            self.n_in),
            ("rows_out",           self.n_out),
            ("invalid_rows",       self.invalid_rows),
            ("solvent_front_rows", self.solvent_front),
            ("low_si_rows",        self.low_si),
            ("trace_rows",         self.trace),
            ("coelutions_merged",  self.coelutions_merged),
            ("singletons_dropped", self.singletons),
            ("synonym_collapses",  len(self.synonym_collapses)),
        ], columns=["metric", "count"]).set_index("metric")


# ---------------------------------------------------------------------------
# Compound-name normalisation
# ---------------------------------------------------------------------------
_STEREO_RE = re.compile(
    r",?\s*\((?:[ZE]|R|S|\+|-|cis|trans|alpha|beta|"
    r"\.alpha\.|\.beta\.|\.gamma\.|\.delta\.)\)-?", re.I)
_PREFIX_N_RE = re.compile(r"^\s*n-(?=[A-Z])")
_TRIM_TRAILING = re.compile(r"[,\s]+$")
_COLLAPSE_WS = re.compile(r"\s{2,}")


def _normalize_name(name: str) -> str:
    """Best-effort canonicalisation of a compound name.

    Conservative — preserves chemistry while flattening the noisiest
    inconsistencies (case of common acid names, stereochemistry tokens,
    prefix 'n-', trailing punctuation, repeated whitespace)."""
    if not isinstance(name, str):
        return name
    n = name.strip()
    # strip stereochemistry tokens like ', (Z)-', '(R)-', '(.alpha.)-'
    n = _STEREO_RE.sub("", n)
    n = _TRIM_TRAILING.sub("", n)
    n = _COLLAPSE_WS.sub(" ", n).strip()
    # strip leading 'n-' for normal-alkyl prefixes
    n = _PREFIX_N_RE.sub("", n)
    # case-fold common acid-name variations
    repls = [
        ("Palmitic Acid", "Palmitic acid"),
        ("Oleic Acid",    "Oleic acid"),
        ("Stearic Acid",  "Stearic acid"),
        ("Linoleic Acid", "Linoleic acid"),
        ("Myristic Acid", "Myristic acid"),
        ("Lauric Acid",   "Lauric acid"),
        ("Behenic Acid",  "Behenic acid"),
    ]
    for a, b in repls:
        if n == a:
            n = b
    return n


def normalize_compound_names(reid: pd.DataFrame
                              ) -> tuple[pd.DataFrame, dict[str, str]]:
    out = reid.copy()
    raw = out["final_name"].astype(str)
    out["final_name"] = raw.map(_normalize_name)
    changed = {a: b for a, b in zip(raw, out["final_name"]) if a != b}
    return out, changed


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------
def drop_invalid_rows(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    n0 = len(df)
    bad = (df["area_pct"].isna() | (df["area_pct"] <= 0)
           | df["rt"].isna() | df["area"].isna())
    out = df[~bad].copy()
    return out, n0 - len(out)


def drop_solvent_front(df: pd.DataFrame, min_rt: float
                        ) -> tuple[pd.DataFrame, int]:
    n0 = len(df)
    out = df[df["rt"] >= min_rt].copy()
    return out, n0 - len(out)


def drop_low_si(df: pd.DataFrame, si_min: int) -> tuple[pd.DataFrame, int]:
    n0 = len(df)
    out = df[df["final_si"].fillna(0) >= si_min].copy()
    return out, n0 - len(out)


def drop_trace_peaks(df: pd.DataFrame, trace_min: float
                      ) -> tuple[pd.DataFrame, int]:
    if trace_min <= 0:
        return df, 0
    n0 = len(df)
    out = df[df["area_pct"] >= trace_min].copy()
    return out, n0 - len(out)


def merge_within_sample_coelutions(df: pd.DataFrame
                                    ) -> tuple[pd.DataFrame, int]:
    """Sum area of multiple peaks of the same compound in the same sample.

    Uses ``compound_key`` (PubChem InChIKey-14 when available, else the
    normalised final_name) so that stereo-isomers and salt forms of the
    same compound are merged correctly.
    """
    if df.empty:
        return df, 0

    key_col = "compound_key" if "compound_key" in df.columns else "final_name"

    # Identify duplicates so we know how many merges occurred
    grp = df.groupby(["sample", key_col])
    n_extra = (grp.size() - 1).clip(lower=0).sum()
    if n_extra == 0:
        return df, 0

    # When merging on compound_key, the final_name we keep is the most
    # frequent name among the merged peaks (so the user sees the
    # PubChem-canonical form when present).
    aggs = {
        "final_name": lambda s: s.value_counts().index[0],
        "peak":       "first",
        "rt":         "median",
        "it":         "min",
        "ft":         "max",
        "area":       "sum",
        "area_pct":   "sum",
        "height":     "max",
        "height_pct": "sum",
        "ah":         "mean",
        "mark":       "first",
        "system_name": "first",
        "system_si":  "max",
        "final_formula": "first",
        "final_mw":   "first",
        "final_si":   "max",
        "final_retindex": "first",
        "hit_chosen": "first",
        "reidentification": "first",
        "is_system_implausible": "first",
        "is_system_contaminant": "first",
        "final_class": "first",
        "is_contaminant_final": "first",
        "is_unidentified": "first",
        "canonical_name": "first",
        "cid":            "first",
        "inchikey":       "first",
        "inchikey14":     "first",
        "formula":        "first",
        "exact_mass":     "first",
    }
    extra_cols = [c for c in ("condition", "day", "fraction", "group")
                   if c in df.columns]
    for c in extra_cols:
        aggs[c] = "first"

    keep_aggs = {k: v for k, v in aggs.items() if k in df.columns}
    out = (df.groupby(["sample", key_col], as_index=False)
              .agg(keep_aggs))
    return out, int(n_extra)


def drop_dataset_singletons(df: pd.DataFrame, max_area: float
                             ) -> tuple[pd.DataFrame, int]:
    """Drop compounds that appear in exactly ONE sample with low abundance."""
    counts = df.groupby("final_name")["sample"].nunique()
    abund = df.groupby("final_name")["area_pct"].max()
    singletons = counts[(counts == 1) & (abund < max_area)].index
    n0 = len(df)
    out = df[~df["final_name"].isin(singletons)].copy()
    return out, n0 - len(out)


# ---------------------------------------------------------------------------
# Master pipeline
# ---------------------------------------------------------------------------
def clean(reid: pd.DataFrame, cfg: QualityConfig | None = None
           ) -> tuple[pd.DataFrame, QualityReport]:
    cfg = cfg or QualityConfig()
    rep = QualityReport(n_in=len(reid))

    df, n = drop_invalid_rows(reid)
    rep.invalid_rows = n

    df, n = drop_solvent_front(df, cfg.min_rt)
    rep.solvent_front = n

    df, n = drop_low_si(df, cfg.si_min)
    rep.low_si = n

    df, n = drop_trace_peaks(df, cfg.trace_min)
    rep.trace = n

    df, syn = normalize_compound_names(df)
    rep.synonym_collapses = syn

    if cfg.merge_coelutions:
        df, n = merge_within_sample_coelutions(df)
        rep.coelutions_merged = n

    if cfg.drop_singletons:
        df, n = drop_dataset_singletons(df, cfg.singleton_max_area)
        rep.singletons = n

    rep.n_out = len(df)
    return df.reset_index(drop=True), rep


def load_quality_config(data_dir: Path) -> QualityConfig:
    """Load optional data/quality.json. Use defaults if absent."""
    import json
    p = data_dir / "quality.json"
    if not p.exists():
        return QualityConfig()
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
    except Exception as e:  # noqa: BLE001
        print(f"WARNING: failed to parse {p}: {e}. Using defaults.")
        return QualityConfig()
    return QualityConfig(**{k: v for k, v in d.items()
                             if k in QualityConfig.__dataclass_fields__})


def write_quality_template(data_dir: Path) -> Path:
    """Write a quality_template.json with default values, if absent."""
    import json
    p = data_dir / "quality_template.json"
    if p.exists():
        return p
    cfg = QualityConfig()
    p.write_text(json.dumps({
        "min_rt":              cfg.min_rt,
        "si_min":              cfg.si_min,
        "trace_min":           cfg.trace_min,
        "drop_singletons":     cfg.drop_singletons,
        "singleton_max_area":  cfg.singleton_max_area,
        "merge_coelutions":    cfg.merge_coelutions,
    }, indent=2), encoding="utf-8")
    return p
