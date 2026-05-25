"""
PubChem-based compound enrichment and InChIKey deduplication.

For every distinct compound name in the re-identified peak list, the
PubChem PUG REST API is queried (free, no API key required) to obtain:

    - canonical IUPAC name        (preferred report name)
    - PubChem CID
    - InChIKey                    (used as the canonical compound ID)
    - molecular formula           (cross-checks the NIST formula)
    - exact mass                  (cross-checks the NIST molweight)

The first 14 characters of the InChIKey (the "connectivity hash") are
identical for stereoisomers, salt forms and tautomers of the same
compound. We therefore deduplicate on this prefix, so that

    'Oleic acid'                  -> CID 445639 -> InChIKey RZJQGNCSTQAWON-...
    'Oleic Acid, (Z)-'            -> CID 445639 -> InChIKey RZJQGNCSTQAWON-...
    '9-Octadecenoic acid (Z)-'    -> CID 445639 -> InChIKey RZJQGNCSTQAWON-...

all collapse to one compound.

The implementation is robust to network failures: every successful
lookup is cached on disk (data/pubchem_cache.json) so re-runs are
essentially instantaneous, and any compound for which no PubChem hit
is returned (or the API call fails) is simply kept under its NIST name.
This means the tool ALWAYS produces a complete output, with or without
an internet connection.
"""
from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional
from urllib.parse import quote

import pandas as pd

try:
    import requests
    _HAS_REQUESTS = True
except Exception:  # noqa: BLE001
    _HAS_REQUESTS = False


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PUBCHEM_TIMEOUT = 8.0
PUBCHEM_THROTTLE = 0.20  # max 5 requests / sec per PUG REST policy
HEADERS = {"User-Agent": "GCMS-Profiling-Tool/1.0 (research)"}


@dataclass
class EnrichmentRecord:
    name_query: str
    cid: Optional[int] = None
    canonical_name: Optional[str] = None
    inchikey: Optional[str] = None
    formula: Optional[str] = None
    exact_mass: Optional[float] = None
    error: Optional[str] = None


def _clean_query(name: str) -> str:
    """Trim NIST-style annotations that confuse PubChem name search."""
    n = str(name).strip()
    n = re.sub(r",\s*\d*TMS derivative\s*$", "", n, flags=re.I)
    n = re.sub(r",\s*\d*TBDMS derivative\s*$", "", n, flags=re.I)
    # remove trailing ', (Z)-' / ', (R)-' / '.alpha.' style stereo
    n = re.sub(r",?\s*\((?:[ZE]|R|S|cis|trans|\+|-)\)-?\s*$", "", n, flags=re.I)
    n = re.sub(r"\.\s*(?:alpha|beta|gamma|delta)\s*\.", "", n, flags=re.I)
    n = re.sub(r"\s+", " ", n).strip()
    return n


def _http_get(session: Any, url: str, *, max_retries: int = 3) -> Any:
    """GET with automatic retry on 429 / 503 (PubChem rate-limit / busy)."""
    last_err = None
    for attempt in range(max_retries):
        try:
            r = session.get(url, timeout=PUBCHEM_TIMEOUT, headers=HEADERS)
            if r.status_code in (429, 503):
                # back off exponentially: 0.5, 1.0, 2.0 s
                time.sleep(0.5 * (2 ** attempt))
                last_err = f"HTTP {r.status_code}"
                continue
            return r
        except Exception as e:  # noqa: BLE001
            last_err = str(e)
            time.sleep(0.5 * (2 ** attempt))
    raise RuntimeError(last_err or "PubChem request failed")


def _pubchem_lookup(name: str, session: Any) -> EnrichmentRecord:
    """Query PubChem for one compound. Returns best-effort record."""
    rec = EnrichmentRecord(name_query=name)
    if not _HAS_REQUESTS:
        rec.error = "requests not installed"
        return rec
    q = _clean_query(name)
    if not q:
        rec.error = "empty query"
        return rec

    # Step 1: name -> CID (with retry)
    try:
        url = f"{PUBCHEM_BASE}/name/{quote(q, safe='')}/cids/JSON"
        r = _http_get(session, url)
        if r.status_code == 404:
            rec.error = "no CID"
            return rec
        r.raise_for_status()
        cids = (r.json().get("IdentifierList", {}) or {}).get("CID", [])
        if not cids:
            rec.error = "no CID"
            return rec
        rec.cid = int(cids[0])
    except Exception as e:  # noqa: BLE001
        rec.error = f"name->CID: {e}"
        return rec

    # Step 2: CID -> properties (with retry)
    try:
        url = (f"{PUBCHEM_BASE}/cid/{rec.cid}/property/"
               f"IUPACName,InChIKey,MolecularFormula,ExactMass/JSON")
        r = _http_get(session, url)
        r.raise_for_status()
        props = (r.json().get("PropertyTable", {}) or {}).get("Properties", [])
        if props:
            p = props[0]
            rec.canonical_name = p.get("IUPACName")
            rec.inchikey = p.get("InChIKey")
            rec.formula = p.get("MolecularFormula")
            try:
                rec.exact_mass = (float(p["ExactMass"])
                                   if "ExactMass" in p else None)
            except Exception:
                rec.exact_mass = None
    except Exception as e:  # noqa: BLE001
        rec.error = f"CID->props: {e}"
    return rec


def _load_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _save_cache(path: Path, data: dict[str, dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def enrich_compounds(names: list[str], cache_path: Path,
                      *, online: bool = True, verbose: bool = True
                      ) -> pd.DataFrame:
    """Enrich a list of NIST names with PubChem properties.

    Returns a DataFrame indexed by the original NIST name with columns:
        cid, canonical_name, inchikey, inchikey14, formula, exact_mass,
        source ('cache', 'pubchem', 'offline', 'no_match', 'error: ...')
    """
    cache = _load_cache(cache_path)
    out_rows: list[dict[str, Any]] = []

    session = None
    if online and _HAS_REQUESTS:
        session = requests.Session()

    n_to_query = sum(1 for n in set(names) if n not in cache)
    if verbose and online and n_to_query:
        print(f"   PubChem enrichment: querying {n_to_query} new compound(s) "
              "(cached results reused)")

    distinct = list(dict.fromkeys(names))  # preserve order, drop dups
    for i, name in enumerate(distinct, 1):
        if name in cache:
            rec = cache[name]
            rec.setdefault("source", "cache")
        elif not online or session is None:
            rec = dict(name_query=name, cid=None, canonical_name=None,
                       inchikey=None, formula=None, exact_mass=None,
                       error=None, source="offline")
        else:
            r = _pubchem_lookup(name, session)
            rec = r.__dict__.copy()
            rec["source"] = ("pubchem" if r.inchikey
                              else f"no_match: {r.error or '-'}")
            cache[name] = {k: v for k, v in rec.items() if k != "source"}
            time.sleep(PUBCHEM_THROTTLE)
            if verbose and i % 25 == 0:
                print(f"     {i} / {len(distinct)} ({rec.get('source', '?')})")
        out_rows.append({
            "nist_name":      name,
            "cid":            rec.get("cid"),
            "canonical_name": rec.get("canonical_name"),
            "inchikey":       rec.get("inchikey"),
            "inchikey14":     (rec.get("inchikey")[:14]
                                if rec.get("inchikey") else None),
            "formula":        rec.get("formula"),
            "exact_mass":     rec.get("exact_mass"),
            "source":         rec.get("source", ""),
            "error":          rec.get("error"),
        })

    if online and session is not None:
        _save_cache(cache_path, cache)
        session.close()

    return pd.DataFrame(out_rows).set_index("nist_name")


def apply_enrichment(reid: pd.DataFrame, enrich_df: pd.DataFrame
                      ) -> tuple[pd.DataFrame, dict[str, str]]:
    """Apply PubChem-based deduplication to a re-identified peak table.

    Two new columns are added:
        canonical_name -- IUPAC name from PubChem if available, else NIST name
        inchikey14     -- the 14-char connectivity hash, used as the unique
                          compound key

    The 'final_name' column is REPLACED by the canonical_name so that
    stereoisomer / salt aliases are merged in every down-stream stat.
    """
    out = reid.copy()
    enrich_df = enrich_df.reset_index()

    # Build one row per (nist_name -> canonical_name, inchikey14)
    map_df = (enrich_df[["nist_name", "canonical_name",
                         "cid", "inchikey", "inchikey14",
                         "formula", "exact_mass"]]
              .drop_duplicates("nist_name"))

    # When PubChem returns a canonical name, prefer it
    out = out.merge(map_df, left_on="final_name", right_on="nist_name",
                    how="left").drop(columns=["nist_name"])

    # final_name choice: keep PubChem canonical when present, else NIST
    new_name = out["canonical_name"].where(
        out["canonical_name"].notna() & (out["canonical_name"] != ""),
        out["final_name"]
    )
    name_changes = {
        a: b for a, b in zip(out["final_name"], new_name)
        if isinstance(a, str) and isinstance(b, str) and a != b
    }
    out["final_name"] = new_name

    # Group key for downstream dedup. When inchikey is missing, fall
    # back to a normalised lower-case name so peaks of unmatched
    # compounds still merge correctly within a sample.
    out["compound_key"] = out["inchikey14"].where(
        out["inchikey14"].notna() & (out["inchikey14"] != ""),
        out["final_name"].astype(str).str.lower().str.strip()
    )

    return out, name_changes
