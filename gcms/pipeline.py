"""
End-to-end orchestrator for the GC-MS analysis tool.

Single entry-point — `run_pipeline(root)` — does everything.
Designed to be robust to:
    * missing or malformed PDFs (parser logs the failure and continues),
    * any number / set of conditions, days and fractions,
    * filenames that the inference heuristics can't parse (they appear
      with blank condition in samples_template.csv for the user to fix),
    * samples.csv that references missing PDFs or omits some PDFs,
    * single-condition or single-group experiments (figures are skipped
      gracefully when ordination is not meaningful),
    * peaks where no plausible alternative exists in the hit list.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import enrich as enrich_mod
from . import figures as fig_mod
from . import metadata as meta_mod
from . import parser as parser_mod
from . import profile as prof
from . import project_info as proj_mod
from . import quality as qc_mod
from . import reidentify as reid_mod
from . import report as report_mod


# ---------------------------------------------------------------------------
def _ensure_dirs(root: Path) -> tuple[Path, Path, Path, Path, Path]:
    data_dir = root / "data"
    pdf_dir = data_dir / "PDFs"
    out_dir = root / "output"
    fig_dir = out_dir / "figures"
    list_dir = out_dir / "lists"
    for d in (data_dir, pdf_dir, out_dir, fig_dir, list_dir):
        d.mkdir(parents=True, exist_ok=True)
    return data_dir, pdf_dir, out_dir, fig_dir, list_dir


def _ordered_groups(meta: pd.DataFrame) -> list[str]:
    """Build a stable group order: by condition first-seen, then day asc."""
    seen = list(dict.fromkeys(meta["condition"].tolist()))
    groups: list[str] = []
    for c in seen:
        days = sorted(meta.loc[meta["condition"] == c, "day"].unique().tolist())
        for d in days:
            groups.append(f"{c}_d{d}")
    return groups


def _print_banner(text: str) -> None:
    print(f"\n{'-' * 78}\n {text}\n{'-' * 78}")


# ---------------------------------------------------------------------------
def run_pipeline(
    root: Path,
    *,
    conditions_for_unique: list[str] | None = None,
    verbose: bool = True,
    online: bool = True,
) -> dict[str, Path]:
    """Run the full analysis. Returns a dict of key output paths."""
    data_dir, pdf_dir, out_dir, fig_dir, list_dir = _ensure_dirs(root)

    # =====================================================================
    # 1. Discover PDFs
    # =====================================================================
    pdf_files = sorted(pdf_dir.glob("*.pdf"))
    if not pdf_files:
        raise FileNotFoundError(
            f"No PDF files found in {pdf_dir}.\n"
            f"  -> Place Shimadzu QUALITATIVE ANALYSIS REPORT PDFs in this "
            f"folder and re-run.")
    if verbose:
        _print_banner(f"Found {len(pdf_files)} PDF(s) in {pdf_dir}")

    # =====================================================================
    # 2. Metadata: auto-template if missing, validate if present
    # =====================================================================
    samples_csv = data_dir / "samples.csv"
    template_csv = data_dir / "samples_template.csv"

    # Load filename-token heuristics from data/conditions.csv if present.
    cond_tokens = meta_mod.load_condition_tokens(data_dir)

    if not samples_csv.exists():
        meta_mod.write_template(pdf_files, template_csv, tokens=cond_tokens)
        token_hint = ""
        if not cond_tokens:
            token_hint = ("\nTip: if your filenames follow a regular pattern "
                          "(e.g. 'AS15-Intra' = condition AS, day 15, "
                          "fraction Intra), you can speed up the template "
                          "by listing your filename tokens in "
                          f"{data_dir / 'conditions.csv'}: a CSV with "
                          "columns 'token,condition,default_day'.\n")
        msg = (f"\nNo samples.csv found at:\n  {samples_csv}\n\n"
               f"A template has been written to:\n  {template_csv}\n\n"
               "Open it, verify the auto-inferred condition / day / "
               "fraction for every sample (blank cells need filling), "
               "save it as samples.csv in the same folder, and re-run "
               "the tool." + token_hint)
        raise SystemExit(msg)

    meta = meta_mod.load_metadata(samples_csv)
    missing, extra = meta_mod.cross_check(meta, pdf_files)
    if missing:
        print(f"WARNING: samples.csv references {len(missing)} file(s) "
              f"that are not in {pdf_dir}: {missing[:5]}"
              f"{'...' if len(missing) > 5 else ''}")
        # ignore missing rows — restrict to PDFs that exist
        meta = meta[meta["sample"].isin({f.stem for f in pdf_files})]
    if extra:
        print(f"WARNING: {len(extra)} PDF(s) not listed in samples.csv "
              f"will be ignored: {extra[:5]}"
              f"{'...' if len(extra) > 5 else ''}")

    if meta.empty:
        raise RuntimeError(
            "After cross-checking samples.csv against the PDFs in "
            f"{pdf_dir}, no valid samples remain. Please verify the "
            "filenames in samples.csv match the PDF names exactly "
            "(without the .pdf extension).")

    group_order = _ordered_groups(meta)
    conditions_present = list(dict.fromkeys(meta["condition"].tolist()))

    if conditions_for_unique is None:
        conditions_for_unique = conditions_present
    else:
        # validate that user-specified conditions exist
        unknown = [c for c in conditions_for_unique
                   if c not in conditions_present]
        if unknown:
            print(f"WARNING: --unique-conditions includes labels not "
                  f"present in samples.csv: {unknown}. Falling back to "
                  f"all conditions.")
            conditions_for_unique = conditions_present

    if verbose:
        print(f"Conditions in metadata: {conditions_present}")
        print(f"Groups (condition x day): {group_order}")

    # Optional user-supplied narrative text
    project_info = proj_mod.load_project_info(data_dir,
                                                conditions=conditions_present)
    # Always offer a starter project_template.json the user can fill in
    proj_mod.write_project_template(data_dir, conditions_present)

    # =====================================================================
    # 3. Parse PDFs (with cache to skip re-parsing on re-runs)
    # =====================================================================
    if verbose:
        _print_banner(f"[1/5] Parsing {len(pdf_files)} PDF(s)")
    pdf_files_used = [f for f in pdf_files if f.stem in set(meta["sample"])]
    if not pdf_files_used:
        raise RuntimeError(
            "No PDFs match the sample names listed in samples.csv. "
            "Sample names must equal the PDF filename without the .pdf "
            "extension.")

    cache_peaks = out_dir / "peaks_raw.csv"
    cache_hits = out_dir / "hits_raw.csv"
    cache_valid = False
    if cache_peaks.exists() and cache_hits.exists():
        cache_mtime = min(cache_peaks.stat().st_mtime,
                          cache_hits.stat().st_mtime)
        newest_pdf = max(f.stat().st_mtime for f in pdf_files_used)
        cache_valid = cache_mtime > newest_pdf

    if cache_valid:
        if verbose:
            print(f"   Reusing parser cache (no PDF newer than cache):")
            print(f"     {cache_peaks}")
            print(f"     {cache_hits}")
        peaks = pd.read_csv(cache_peaks)
        hits = pd.read_csv(cache_hits)
        parser_log = [f"{s:<25s}  peaks={(peaks['sample']==s).sum():4d}  "
                      f"hits={(hits['sample']==s).sum():5d}"
                      for s in sorted(peaks["sample"].unique())]
        if verbose:
            for line in parser_log:
                print(" ", line)
    else:
        peaks, hits, parser_log = parser_mod.parse_pdf_directory(pdf_dir)
        if verbose:
            for line in parser_log:
                print(" ", line)
        if not peaks.empty:
            peaks.to_csv(cache_peaks, index=False)
            hits.to_csv(cache_hits, index=False)

    if peaks.empty:
        raise RuntimeError(
            "PDF parsing produced 0 peaks. Verify the PDFs are Shimadzu "
            "QUALITATIVE ANALYSIS REPORT files (the format produced by "
            "GCMSsolution).")

    peaks = peaks[peaks["sample"].isin(meta["sample"])]
    hits = hits[hits["sample"].isin(meta["sample"])]

    samples_with_no_peaks = (set(meta["sample"]) - set(peaks["sample"]))
    if samples_with_no_peaks:
        print(f"WARNING: {len(samples_with_no_peaks)} sample(s) yielded no "
              f"peaks during parsing: {sorted(samples_with_no_peaks)}")

    # =====================================================================
    # 4. Re-identification
    # =====================================================================
    if verbose:
        _print_banner(f"[2/5] Re-identifying {len(peaks)} peaks")
    reid = reid_mod.reidentify(peaks, hits)
    reid_with_meta = reid.merge(
        meta[["sample", "condition", "day", "fraction", "group"]],
        on="sample", how="left")
    reid_with_meta.to_csv(out_dir / "peaks_reidentified.csv", index=False)

    n_kept = (reid_with_meta["reidentification"] == "kept_system_pick").sum()
    n_total = len(reid_with_meta)
    n_changed = n_total - n_kept
    if verbose:
        print(f"   Kept system pick : {n_kept}")
        print(f"   Re-assigned      : {n_changed} "
              f"({100 * n_changed / max(1, n_total):.1f} %)")

    corr = reid_mod.correction_summary(reid_with_meta)
    if not corr.empty:
        corr.to_csv(out_dir / "reidentification_corrections.csv", index=False)

    # =====================================================================
    # 4a. PubChem enrichment + InChIKey-based deduplication
    # =====================================================================
    if verbose:
        _print_banner("[2a/5] PubChem enrichment (canonical names "
                      "& InChIKey deduplication)")
    cache_path = data_dir / "pubchem_cache.json"
    distinct_names = sorted(set(reid_with_meta["final_name"].dropna()
                                  .astype(str).tolist()))
    enrich_df = enrich_mod.enrich_compounds(
        distinct_names, cache_path, online=online, verbose=verbose)
    enrich_df.to_csv(out_dir / "pubchem_enrichment.csv")

    # join PubChem props onto reid + replace final_name with PubChem
    # canonical and add compound_key (InChIKey-14 or normalised name)
    reid_enriched, name_changes = enrich_mod.apply_enrichment(
        reid_with_meta, enrich_df)

    if name_changes:
        pd.DataFrame(name_changes.items(),
                     columns=["nist_name", "canonical_name"]
                     ).drop_duplicates().to_csv(
                         out_dir / "pubchem_name_remap.csv", index=False)

    if verbose:
        n_with_inchi = enrich_df["inchikey"].notna().sum()
        n_total = len(enrich_df)
        print(f"   PubChem matches      : {n_with_inchi} / {n_total} "
              f"({100*n_with_inchi/max(1,n_total):.1f} %)")
        print(f"   Name remaps applied  : {len(name_changes)}")

    # =====================================================================
    # 4b. Data quality cleanup (deduplication, error filtering)
    # =====================================================================
    qc_cfg = qc_mod.load_quality_config(data_dir)
    qc_mod.write_quality_template(data_dir)
    if verbose:
        _print_banner("[2b/5] Data quality cleanup and deduplication")
    reid_clean, qrep = qc_mod.clean(reid_enriched, qc_cfg)
    qrep_df = qrep.as_dataframe()
    qrep_df.to_csv(out_dir / "quality_report.csv")
    # save the synonym-collapse log too for full audit
    if qrep.synonym_collapses:
        pd.DataFrame(
            list(qrep.synonym_collapses.items()),
            columns=["original_name", "normalised_name"]
        ).drop_duplicates().to_csv(
            out_dir / "synonym_collapses.csv", index=False)
    if verbose:
        for label, n in (("Invalid (NaN/<=0)",  qrep.invalid_rows),
                         ("Solvent front",      qrep.solvent_front),
                         ("Low-SI noise",       qrep.low_si),
                         ("Trace peaks",        qrep.trace),
                         ("Co-elutions merged", qrep.coelutions_merged),
                         ("Singletons dropped", qrep.singletons)):
            print(f"   {label:<22s}: {n}")
        print(f"   Synonym collapses     : {len(qrep.synonym_collapses)}")
        print(f"   Rows in -> out        : {qrep.n_in} -> {qrep.n_out}")

    # =====================================================================
    # 5. Profiling
    # =====================================================================
    if verbose:
        _print_banner("[3/5] Building abundance matrix and diversity stats")
    # Persist the cleaned table for full traceability
    reid_clean.to_csv(out_dir / "peaks_cleaned.csv", index=False)

    mat = prof.build_abundance_matrix(reid_clean, meta, group_order)
    mat.to_csv(out_dir / "metabolite_abundance_matrix.csv")

    if mat.empty:
        raise RuntimeError(
            "Every peak was filtered out (Unidentified or Contaminant). "
            "Cannot build the abundance matrix. Check that your re-"
            "identification rules are not too aggressive for this "
            "dataset.")

    div = prof.diversity_indices(mat)
    div.to_csv(out_dir / "diversity_indices.csv")

    class_ab = prof.class_abundance(reid_clean, meta, group_order)
    class_ab.to_csv(out_dir / "class_abundance.csv")

    class_rich = prof.class_richness(reid_clean, meta, group_order)
    class_rich.to_csv(out_dir / "class_richness.csv")

    presence = prof.shared_unique(mat)
    presence.to_csv(out_dir / "metabolite_presence.csv")

    buckets = prof.core_vs_unique(mat, conditions_for_unique, group_order)
    for k, v in buckets.items():
        v.to_csv(out_dir / f"set_{k}.csv")

    top40 = prof.top_abundant_overall(mat, top_n=40)
    top40.to_csv(list_dir / "top40_most_abundant.csv")

    top5 = prof.top_n_per_group(mat, n=5)
    top5.to_csv(list_dir / "top5_per_group.csv", index=False)

    cond_unique = prof.cond_unique_3way(mat, conditions_for_unique, group_order)
    for c, df in cond_unique.items():
        df.to_csv(list_dir / f"unique_{c}.csv", index_label="metabolite")

    if verbose:
        print(f"   Distinct metabolites : {mat.shape[0]}")
        print(f"   Groups               : {mat.shape[1]} -> {group_order}")
        print(f"   Conditions for unique table: {conditions_for_unique}")
        print(f"   Core (shared by all conditions): "
              f"{len(buckets.get('core_metabolites', []))}")
        for c in conditions_for_unique:
            print(f"   Unique to {c:<12s}: "
                  f"{len(cond_unique.get(c, pd.DataFrame()))}")

    # =====================================================================
    # 6. Figures
    # =====================================================================
    if verbose:
        _print_banner(f"[4/5] Rendering figures into {fig_dir}")
    fig_mod.fig_pca(mat, fig_dir / "01_PCA.png")
    fig_mod.fig_class_diversity_bar(class_rich,
                                     fig_dir / "02_class_diversity_bar.png")
    fig_mod.fig_alpha_diversity(div, fig_dir / "03_alpha_diversity.png")
    fig_mod.fig_heatmap_top(mat, fig_dir / "04_heatmap_top40_metabolites.png",
                             top_n=40)
    fig_mod.fig_heatmap_classes(class_ab, fig_dir / "05_heatmap_classes.png")
    fig_mod.fig_unique_vs_core(buckets, conditions_for_unique,
                                fig_dir / "06_unique_vs_core.png")

    # =====================================================================
    # 7. Word reports
    # =====================================================================
    if verbose:
        _print_banner(f"[5/5] Writing Word reports")
    main_docx = out_dir / "GCMS_Metabolic_Profiling_Report.docx"
    try:
        report_mod.write_main_report(
            out_path=main_docx,
            figures_dir=fig_dir,
            sample_meta=meta,
            parser_log=parser_log,
            reid=reid_clean,
            corrections=corr,
            mat=mat,
            div=div,
            class_ab=class_ab,
            class_rich=class_rich,
            buckets=buckets,
            top40=top40,
            top5=top5,
            cond_unique=cond_unique,
            conditions=conditions_for_unique,
            group_order=group_order,
            project_info=project_info,
            quality_report=qrep,
        )
    except PermissionError as e:
        alt = out_dir / "GCMS_Metabolic_Profiling_Report__new.docx"
        print(f"WARNING: cannot overwrite the open Word file ({e}); "
              f"writing to {alt} instead.")
        report_mod.write_main_report(
            out_path=alt,
            figures_dir=fig_dir,
            sample_meta=meta,
            parser_log=parser_log,
            reid=reid_clean,
            corrections=corr,
            mat=mat,
            div=div,
            class_ab=class_ab,
            class_rich=class_rich,
            buckets=buckets,
            top40=top40,
            top5=top5,
            cond_unique=cond_unique,
            conditions=conditions_for_unique,
            group_order=group_order,
            project_info=project_info,
            quality_report=qrep,
        )
        main_docx = alt

    unique_docx = out_dir / "Unique_metabolites_table.docx"
    try:
        report_mod.write_unique_table(
            out_path=unique_docx,
            cond_unique=cond_unique,
            conditions=conditions_for_unique,
            project_info=project_info,
        )
    except PermissionError as e:
        alt = out_dir / "Unique_metabolites_table__new.docx"
        print(f"WARNING: cannot overwrite the open Word file ({e}); "
              f"writing to {alt} instead.")
        report_mod.write_unique_table(
            out_path=alt,
            cond_unique=cond_unique,
            conditions=conditions_for_unique,
            project_info=project_info,
        )
        unique_docx = alt

    if verbose:
        _print_banner("Done")
        print(f"   Main report  : {main_docx}")
        print(f"   Unique table : {unique_docx}")
        print(f"   Figures      : {fig_dir}")
        print(f"   CSV tables   : {out_dir}")

    return {
        "main_report": main_docx,
        "unique_table": unique_docx,
        "figures_dir": fig_dir,
        "output_dir": out_dir,
    }
