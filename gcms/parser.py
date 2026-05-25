"""
PDF parser for Shimadzu QUALITATIVE ANALYSIS reports.

Extracts:
  * Peak table  -> peak#, R.Time, area, area%, height, height%, mark, name
  * Library hits -> Hit#, Entry, Library, SI, Formula, CAS, MW, RetIndex,
                     CompName (one row per Hit per peak)
"""
from __future__ import annotations

import os
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import pdfplumber

PEAK_RE = re.compile(
    r"^\s*(?P<peak>\d+)\s+"
    r"(?P<rt>\d+\.\d+)\s+(?P<it>\d+\.\d+)\s+(?P<ft>\d+\.\d+)\s+"
    r"(?P<area>\d+)\s+(?P<area_pct>\d+\.\d+)\s+"
    r"(?P<height>\d+)\s+(?P<height_pct>\d+\.\d+)\s+"
    r"(?P<ah>\d+\.\d+)\s+"
    r"(?P<rest>.+)$"
)
LINE_HEADER_RE = re.compile(r"Line#:(\d+)\s+R\.Time:(\d+\.\d+)\(Scan#:(\d+)\)")
HIT_HEADER_RE  = re.compile(r"Hit#:(\d+)\s+Entry:(\d+)\s+Library:(\S+)")
SI_RE = re.compile(
    r"SI:(\d+)\s+Formula:(\S+)\s+CAS:(\S+)\s+MolWeight:(\d+)\s+RetIndex:(\d+)"
)
COMP_RE = re.compile(r"CompName:(.+?)(?:\$\$|$)")

_NON_PEAK_PREFIXES = (
    "Peak#", "Library", "<<", "Line#", "Hit#", "Total", "===",
    "Chromatogram", "[", "BasePeak", "BG Mode", "RawMode",
    "Group ", "MassPeaks", "100", "80", "60", "40", "20",
)


def _parse_peak_table(text: str) -> pd.DataFrame:
    rows = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = PEAK_RE.match(lines[i])
        if not m:
            i += 1
            continue
        rest = m.group("rest").strip()
        mark = ""
        if rest.startswith(("V ", "MI ", "M ")):
            tok = rest.split(maxsplit=1)
            mark, rest = tok[0], tok[1] if len(tok) > 1 else ""
        # join wrapped name continuations
        j = i + 1
        while (
            j < len(lines)
            and lines[j].strip()
            and not PEAK_RE.match(lines[j])
            and not lines[j].lstrip().startswith(_NON_PEAK_PREFIXES)
            and not re.match(r"^\s*\d+\s+\d", lines[j])
        ):
            rest = (rest + " " + lines[j].strip()).strip()
            j += 1
        rows.append(dict(
            peak=int(m.group("peak")),
            rt=float(m.group("rt")),
            it=float(m.group("it")),
            ft=float(m.group("ft")),
            area=int(m.group("area")),
            area_pct=float(m.group("area_pct")),
            height=int(m.group("height")),
            height_pct=float(m.group("height_pct")),
            ah=float(m.group("ah")),
            mark=mark,
            system_name=rest.strip(),
        ))
        i = j
    return pd.DataFrame(rows)


def _parse_library_hits(text: str) -> pd.DataFrame:
    rows: list[dict] = []
    blocks = re.split(r"<<\s*Target\s*>>", text)
    current_line, current_rt = None, None
    for blk in blocks:
        mh = LINE_HEADER_RE.search(blk)
        if mh:
            current_line = int(mh.group(1))
            current_rt = float(mh.group(2))
        for hm in HIT_HEADER_RE.finditer(blk):
            hit_no = int(hm.group(1))
            tail = blk[hm.end(): hm.end() + 1500]
            sim = SI_RE.search(tail)
            cmp_ = COMP_RE.search(tail)
            if not sim or not cmp_:
                continue
            comp = cmp_.group(1).strip()
            primary = comp.split("$$")[0].strip().rstrip(",").strip()
            rows.append(dict(
                line=current_line,
                rt=current_rt,
                hit=hit_no,
                entry=int(hm.group(2)),
                library=hm.group(3),
                si=int(sim.group(1)),
                formula=sim.group(2),
                cas=sim.group(3),
                molweight=int(sim.group(4)),
                retindex=int(sim.group(5)),
                compname=primary,
                compname_full=comp[:300],
            ))
    return pd.DataFrame(rows)


def parse_pdf(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (peaks_df, hits_df) for one Shimadzu PDF."""
    with pdfplumber.open(str(path)) as pdf:
        text = "\n".join((p.extract_text() or "") for p in pdf.pages)
    text = text.replace("\u00a0", " ")
    return _parse_peak_table(text), _parse_library_hits(text)


def _parse_one(path_str: str) -> tuple[str, pd.DataFrame, pd.DataFrame, str | None]:
    """Worker for parallel parsing. Returns (sample, peaks, hits, error)."""
    path = Path(path_str)
    sample = path.stem
    try:
        peaks, hits = parse_pdf(path)
        peaks["sample"] = sample
        hits["sample"] = sample
        return sample, peaks, hits, None
    except Exception as e:  # noqa: BLE001
        return sample, pd.DataFrame(), pd.DataFrame(), str(e)


def parse_pdf_directory(pdf_dir: Path,
                         *, n_workers: int | None = None
                         ) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Parse every *.pdf in pdf_dir in parallel.

    Returns (peaks_all, hits_all, file_log). The number of worker
    processes defaults to min(8, cpu_count), but can be overridden.
    """
    files = sorted(pdf_dir.glob("*.pdf"))
    if not files:
        return pd.DataFrame(), pd.DataFrame(), []

    if n_workers is None:
        try:
            cpu = os.cpu_count() or 1
        except Exception:
            cpu = 1
        n_workers = max(1, min(8, cpu - 1, len(files)))

    peaks_all, hits_all, log = [], [], []
    if n_workers <= 1:
        # Serial fallback
        for f in files:
            sample, p, h, err = _parse_one(str(f))
            if err:
                log.append(f"FAILED {sample}: {err}")
                continue
            peaks_all.append(p)
            hits_all.append(h)
            log.append(f"{sample:<25s}  peaks={len(p):4d}  hits={len(h):5d}")
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            for sample, p, h, err in ex.map(_parse_one,
                                             [str(f) for f in files]):
                if err:
                    log.append(f"FAILED {sample}: {err}")
                    continue
                peaks_all.append(p)
                hits_all.append(h)
                log.append(f"{sample:<25s}  peaks={len(p):4d}  "
                           f"hits={len(h):5d}")

    if not peaks_all:
        return pd.DataFrame(), pd.DataFrame(), log
    return (pd.concat(peaks_all, ignore_index=True),
            pd.concat(hits_all, ignore_index=True),
            log)
