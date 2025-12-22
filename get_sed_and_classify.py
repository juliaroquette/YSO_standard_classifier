
#! /usr/bin/env python


from sedan.vizierSED import Source, load_black_list
import os
import pandas as pd
from tqdm import trange

def main(
    input_file , # file with list of sources to query
    output_dir="./SED/", # directory for saving outputs
    ID_col="Internal_ID",          # column name for source ID in input_file
    ra_col="RA",         # column name for RA in input_file
    dec_col="DE",        # column name for DEC in input_file
    search_radius=3.0,   # search radius in arcseconds
    raw_cache_dir=None,        # None → no raw file saved
    processed_dir="processed", # sub-directory for exporting processed SEDs
    force_download=False, # if True, always query VizieR even if raw file exists
    save_processed=True,   # if True, save processed SEDs to disk
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

    if out_raw_cache is not None:
        os.makedirs(out_raw_cache, exist_ok=True)

    if save_processed:
        os.makedirs(out_processed, exist_ok=True)

    # read coordinates
    coords = pd.read_csv(coordspath)

    # load survey black list config
    catalog_dir = os.path.dirname(os.path.abspath(input_file))
    blacklist_config = load_black_list(catalog_dir)

    for i in trange(len(coords)):
        ra = coords.loc[i, ra_col]
        dec = coords.loc[i, dec_col]
        int_id = coords.loc[i, ID_col]
        # 
        # ---- Query VizieR SED ----
        # Just instanciate it will:
        # 1. query VizieR-SED and retrieve the data
        # 2. Add field for distance to central coordinate
        # 3. do basic processing (drop invalid freq, zero error, negative flux)
        # 4. convert units to lambda [microns], flux [Jy], eflux [Jy] 
        star = Source(
            ra,
            dec,
            search_radius,
            ID=int_id,
            raw_cache_dir=out_raw_cache,
            force_download=force_download,
            black_list=blacklist_config,
        )
        # ---- Optional post-processing steps ----
        # __ To be commented/uncommented as needed __
        if star.sed is None:
            continue

        #
        # -- Target wavelength-range selection --
        #
        # - Keeps only measurements within a desired wavelength range
        #
        star.keep_wavelength_range(2., 24.)
        #
        # - Alternatively, drop measurements within a given range
        # for example, in Orion we dropped most optical data to 
        # avoid including light-curves
        # 
        # star.drop_wavelength_range(.4, 1.0)
        #
        # -- Confine to smaller radius --
        #
        # star.confine(1.0)  # arcseconds
        # 
        # -- Drop Filter --
        # To drop all data in a given filter
        #
        #  star.dropFilter('W4')
        #
        # -- Drop Catalog --
        #
        # to drop all data from a given catalog
        # 
        # star.dropCatalog('II/349/apass9')
        #
        # -- Dealing with duplicates/repeats --
        #
        # - Deal with multiple detections
        # If multiple sources are detected from the same survey in the same filter,
        # only the closest one to the central coordinate is kept,
        # and it is flagged with 'mult=True' in the output.
        #
        star.resolveMultipleDetections()
        #
        # - Merge repeated data points
        # Data points are considered repeated if they have the same fitler name,
        # coordinate (down to tol_arcsec), and flux (down to flux_tol_frac)
        #
        star.mergeRepeatedPoints(tol_arcsec=0.1, flux_tol_frac=1e-3)
        # 
        # - Aggregate all measurements by filter
        # If multiple measurements exist in the same filter, they are combined
        # into a single point using weighted average.
        # 
        star.aggregateFilter()
        #
        # -- Remove known problematic surveys --
        #
        # Surveys listed in the YAML config file are removed here
        # the name of the YAML file at the same location as the input catalog
        #
        star.remove_outdated()
        star.remove_weird()
        star.remove_custom()
        #
        # -- Write processed SED to disk --
        if save_processed:
            star.save(out_processed, suffix="_processed")

if __name__ == "__main__":
    main('./nemesis_example_alphair.csv',
         output_dir="SED", # directory for saving outputs
         ID_col="Internal_ID",          # column name for source ID in input_file
         ra_col="RA",         # column name for RA in input_file
         dec_col="DE",        # column name for DEC in input_file
         search_radius=3.0,   # search radius in arcseconds
         raw_cache_dir='raw',        # None → no raw file saved
         processed_dir="processed", # sub-directory for exporting processed SEDs
         force_download=False, # if True, always query VizieR even if raw file exists
         save_processed=True,   # if True, save processed SEDs to disk)
    )
    exit(0)