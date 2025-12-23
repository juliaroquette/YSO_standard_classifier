#! /usr/bin/env python
#
# @juliaroquette 22 December 2025:
# This is an example script to query VizieR SED service for a list of sources,
#
import os
import numpy as np
import pandas as pd
from tqdm import trange

from sed.vizierSED import Source
from sed.config_blacklist import load_black_list, apply_blacklist
from sed.classifyYSO import YSOClassifier   


def main(
    input_file,                 # CSV file with list of sources
    output_dir="./output/",      # directory for saving outputs
    ID_col="Internal_ID",        # column name for source ID in input_file
    ra_col="RA",                # column name for RA in input_file
    dec_col="DE",               # column name for DEC in input_file
    search_radius=3.0,          # search radius in arcseconds
    raw_cache_dir=None,         # None -> no raw file saved
    processed_dir="processed",  # sub-directory for exporting processed SEDs
    force_download=False,       # if True, always query VizieR even if raw file exists
    save_processed=True,        # if True, save processed SEDs to disk

    # NEW: alpha_ir outputs
    save_alpha_table=True,
    alpha_outfile="alpha_ir.csv",
    fit_prefix="_fit",
    make_plots=False,
    plot_dir="plot",

    # Blacklist config (Option A)
    blacklist_yaml=None,        # path to YAML file; if None, defaults to <repo>/input/black_list.yaml
    blacklist_keys=None,        # e.g. ["outdated_surveys"] or ["outdated_surveys","weird_surveys"]; None -> apply all YAML keys
):
    # location of input file
    coordspath = os.path.abspath(input_file)
    # location for exporting seds
    out_root = os.path.abspath(output_dir)
    # location of raw seds if any
    out_raw_cache = (
        os.path.join(out_root, raw_cache_dir)
        if raw_cache_dir is not None
        else None
    )
    # processed output directory
    out_processed = os.path.join(out_root, processed_dir)

    # NEW: plot + alpha output paths
    out_plot = os.path.join(out_root, plot_dir)
    out_alpha = os.path.join(out_root, alpha_outfile)

    if out_raw_cache is not None:
        os.makedirs(out_raw_cache, exist_ok=True)

    if save_processed:
        os.makedirs(out_processed, exist_ok=True)

    if make_plots:
        os.makedirs(out_plot, exist_ok=True)

    # read coordinates
    coords = pd.read_csv(coordspath)

    required_cols = [ID_col, ra_col, dec_col]
    for col in required_cols:
        if col not in coords.columns:
            raise ValueError("Missing required column in input_file: %s" % col)

    # NEW: initialize output columns
    if save_alpha_table:
        coords["alpha_ir" + fit_prefix] = np.nan
        coords["intercept" + fit_prefix] = np.nan
        coords["wl_min" + fit_prefix] = np.nan
        coords["wl_max" + fit_prefix] = np.nan
        coords["n_points" + fit_prefix] = 0
        coords["class" + fit_prefix] = ""

    # load survey black list config
    if blacklist_yaml is None:
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        blacklist_yaml = os.path.join(repo_root, "input", "black_list.yaml")

    blacklist_config = load_black_list(blacklist_yaml)

    for i in trange(len(coords)):
        ra = coords.loc[i, ra_col]
        dec = coords.loc[i, dec_col]
        int_id = coords.loc[i, ID_col]

        # ---- Query VizieR SED ----
        star = Source(
            ra,
            dec,
            search_radius,
            ID=int_id,
            raw_cache_dir=out_raw_cache,
            force_download=force_download,
        )
        if star.sed is None:
            continue

        apply_blacklist(star, blacklist_config, keys=blacklist_keys)

        # ---- Optional post-processing steps ----
        star.keep_wavelength_range(.5, 50.)
        star.resolveMultipleDetections()
        star.mergeRepeatedPoints(tol_arcsec=0.1, flux_tol_frac=1e-3)
        star.aggregateFilter()

        # ---- Alpha_IR classification (NEW) ----
        if save_alpha_table:
            clf = YSOClassifier(star, lower=2.0, upper=24.0)
            res = clf.run()

            coords.loc[i, "alpha_ir" + fit_prefix] = res["alpha_ir"]
            coords.loc[i, "intercept" + fit_prefix] = res["intercept"]
            coords.loc[i, "wl_min" + fit_prefix] = res["wl_min"]
            coords.loc[i, "wl_max" + fit_prefix] = res["wl_max"]
            coords.loc[i, "n_points" + fit_prefix] = res["n_points"]
            coords.loc[i, "class" + fit_prefix] = res["class"]

            if make_plots:
                plot_name = str(int_id) + ".png"
                clf.plot(savepath=os.path.join(out_plot, plot_name))

        # -- Write processed SED to disk --
        if save_processed:
            star.save(out_processed, suffix="_processed")

    # NEW: save augmented coordinates table
    if save_alpha_table:
        coords.to_csv(out_alpha, index=False)


if __name__ == "__main__":
    main(
        "./nemesis_example_alphair.csv",
        output_dir="./output/",
        ID_col="Internal_ID",
        ra_col="RA",
        dec_col="DE",
        search_radius=3.0,
        raw_cache_dir="raw",
        processed_dir="processed",
        force_download=False,
        save_processed=True,

        # NEW:
        save_alpha_table=True,
        alpha_outfile="alpha_ir.csv",
        fit_prefix="_fit",
        make_plots=True,
        plot_dir="plot",

        # blacklist_yaml=None means: <repo_root>/input/black_list.yaml
        blacklist_yaml=None,
        blacklist_keys=["outdated_surveys"],
    )
