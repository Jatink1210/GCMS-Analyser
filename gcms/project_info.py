"""
Optional, user-supplied project metadata.

If `data/project.json` exists, its values are used to personalise the
Word report (title, subtitle, sample-source description, extraction
protocol, condition descriptions, footer). If it is absent, fully
generic placeholder text is used so the report still works for any
matrix and any set of conditions without any user input.

Recognised keys (all optional):

{
    "title":        "<report title>",
    "subtitle":     "<one-line subtitle>",
    "sample_source": "<one-paragraph description of where samples come from>",
    "extraction":    "<extraction / sample-prep description>",
    "fraction_note": "<comment on how Intra/Extra are merged, etc.>",
    "conditions": {
        "<condition_label>": "<one-line description of this condition>"
    },
    "notes":         "<free-form notes appended at the end>"
}

When a key is missing the report falls back to a neutral default that
references only the data itself (number of samples, conditions, days,
metabolites, etc.).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

DEFAULT_PROJECT_INFO: dict[str, Any] = {
    "title": "GC-MS Metabolic Profiling Report",
    "subtitle": ("Re-identification of NIST library mismatches, metabolite "
                 "abundance matrix, alpha-diversity, PCA, compound-class "
                 "statistics, and shared / unique metabolite tables"),
    "sample_source": "",
    "extraction": "",
    "fraction_note": ("If multiple fractions are present (for example "
                       "intracellular and extracellular extracts of the same "
                       "biological replicate), they are merged at the "
                       "(condition, day) level by summation of the area-% "
                       "values."),
    "conditions": {},
    "notes": "",
}


def _example_project_info(conditions: list[str]) -> dict[str, Any]:
    cd = {c: f"<one-line description of condition '{c}'>" for c in conditions}
    return {
        "title": "GC-MS Metabolic Profiling Report",
        "subtitle": ("Re-identification, abundance, diversity, PCA and "
                      "shared / unique metabolite analysis"),
        "sample_source": ("<one paragraph describing the biological "
                           "source of the samples — organism, growth "
                           "conditions, harvesting, etc.>"),
        "extraction": ("<one paragraph describing the extraction and "
                        "sample-preparation protocol — solvents, "
                        "concentration, derivatisation if any, "
                        "redissolution solvent, injection volume, etc.>"),
        "fraction_note": ("If your dataset contains multiple fractions "
                           "per biological sample (e.g. intracellular and "
                           "extracellular), they are merged at the "
                           "(condition, day) level by summation."),
        "conditions": cd,
        "notes": "",
    }


def load_project_info(data_dir: Path,
                       conditions: list[str] | None = None) -> dict[str, Any]:
    pj = data_dir / "project.json"
    info = dict(DEFAULT_PROJECT_INFO)
    if pj.exists():
        try:
            user = json.loads(pj.read_text(encoding="utf-8"))
            if isinstance(user, dict):
                # shallow merge — user's keys override defaults
                info.update({k: v for k, v in user.items()
                             if v not in (None, "", {})})
        except Exception as e:  # noqa: BLE001
            print(f"WARNING: failed to parse {pj}: {e}. Using default text.")
    return info


def write_project_template(data_dir: Path, conditions: list[str]) -> Path:
    """Write a project_template.json the user can fill in."""
    pj = data_dir / "project_template.json"
    if pj.exists():
        return pj
    pj.write_text(json.dumps(_example_project_info(conditions), indent=2),
                  encoding="utf-8")
    return pj
