#! /usr/bin/env python
"""
@juliaroquette 4 December 2025
This is a generalised script to query VizieR SED service for a list of sources
and produce SEDs.

This is based on the original code by David Hernandez (https://github.com/Starvexx/SED_query)
which was tailored for the Horizon 2020 project NEMESIS (https://nemesis.konkoly.hu/).

Basic class `source` queries VizieR SED for a given coordinate and search radius, retrieves the data,
processes it to remove invalid points, and saves it as a csv file.

For example, to query a single source:

ra = 83.86937 
dec = -4.79072
radius = 2.0 # arcseconds
ID = "V1547Ori"
raw_path = "./SED/raw_test/"
star = source(ra, dec, radius, raw_path, ID=ID)

This will save a file `./SED/raw_test/V1547Ori.csv` with the SED data queried from VizieR SED, 
with basic processing to remove invalid points, and converted to wavelength in microns, and flux density in Jansky.

"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from tqdm import tqdm
import astropy.constants as const
import yaml
from astropy.constants import c as c_light

class Source:
    def __init__(self, ra, dec, radius,
                  ID=None, 
                  raw_cache_dir=None,
                  force_download=False,
                  black_list=None):
        """
        Class that instanciate a single SED source around (ra, dec) by 
        querying a given coordinate using VizieR SED service 
        (https://vizier.cds.unistra.fr/vizier/sed/assets/sedapi.gml)

        Behavior:
        - If raw_cache_dir is None:
            * Always query VizieR
            * Never load or save raw SED files automatically
        - If raw_cache_dir is provided:
            * If force_download is False and a raw cache file exists:
                - Load the cached raw SED
            * Otherwise:
                - Query VizieR
                - Automatically save the raw SED to raw_cache_dir        

        Parameters
        ----------
        ra, dec : float
            ICRS coordinates in degrees.
        radius : float
            Search query radius in arcseconds.
        raw_path : str
            Directory for outputting raw CSV files.
        ID : int or str, optional
            Internal source identifier. If None, use 'ra_dec' as filename.
        raw_cache_dir : str or None
            Directory where raw SED cache files are stored.  
            If None, now raw file is downloaded.
        force_download : bool, optional
            If True, always query VizieR even if a cached raw SED exists.
        """

        self.ra = ra
        self.dec = dec
        self.r = radius
        self.r_confine = None 
        self.internal_id = ID
        self.raw_cache_dir = raw_cache_dir
        self.force_download = force_download
        self.black_list = black_list

        # define central coordinate for query
        self.centerCoord = SkyCoord(
            ra=ra * u.degree,
            dec=dec * u.degree,
            frame="icrs"
        )

        # initialize SED attributes
        self.sed = None

        # define target part for the vizier url
        self.__target = f"{ra},{dec}"


        # gets the name of the file
        raw_file = self._rawCacheFile()

        if (
            raw_file is not None
            and os.path.exists(raw_file)
            and not self.force_download
        ):
            if self._loadSedFromFile(raw_file):
                return
                # if the file is loaded, we are skipping all below 
        self.__getData()

        if self.sed is None:
            return

        self.__toDataframe()
        self.__computeCoordDistance()  
        # cary out basic processing
        self.__dropInvalidFrequency()
        self.__dropZeroError()
        self.__dropNegativeFlux()
        self.__convertData()

        # save raw SED to cache
        if self.raw_cache_dir is not None and self.sed is not None:
            self.save(self.raw_cache_dir, suffix="_raw")



    def _baseFilename(self):
        """
        Define base filename for the source SED file.
        """
        if self.internal_id is None:
            return f"{self.ra:.5f}_{self.dec:.5f}_r{self.r:.2f}arcsec"
        return str(self.internal_id)

    def _rawCacheFile(self):
        """
        Define raw cache filename for the source SED file.
        """
        if self.raw_cache_dir is None:
            return None
        return os.path.join(
            self.raw_cache_dir,
            f"{self._baseFilename()}_raw.csv"
        )

    def _loadSedFromFile(self, path):
        """
        Load SED data from a CSV file.
        """
        try:
            self.sed = pd.read_csv(path)
            return True
        except Exception as e:
            tqdm.write(f"Error! Could not read {path}: {e}")
            self.sed = None
            return False

    def __getData(self):
        """
        Query VizieR for the SED
        """
        url = (
            "https://vizier.cds.unistra.fr/"
            f"viz-bin/sed?-c={self.__target}&-c.rs={self.r}"
        )        
        try:
            # instructions of SED rely on astropy Tables
            sed = Table.read(url)
        except Exception as e:
            tqdm.write(f"Error! VizieR SED retrieval failed for {self.__target}: {e}")
            self.sed = None
            return
                
        self.sed = sed


    def save(self, out_dir, *, suffix):
        """
        Save the current state of the SED to disk.

        Parameters
        ----------
        out_dir : str
            Output directory.
        suffix : str
            Suffix to append to the filename (e.g. '_raw', '_processed').
        """
        if self.sed is None:
            return

        os.makedirs(out_dir, exist_ok=True)

        fname = f"{self._baseFilename()}{suffix}.csv"
        out_file = os.path.join(out_dir, fname)

        self.sed.to_csv(out_file, index=False)


    def __toDataframe(self):
        """
        Convert internal SED Astropy Table 
        to a pandas DataFrame."""
        if self.sed is None:
            return
        if isinstance(self.sed, Table):
            self.sed = self.sed.to_pandas()\

    def __dropZeroError(self):
        """
        Drop rows with null negative or NaN flux error. 
        """
        if self.sed is None:
            return

        if 'sed_eflux' not in self.sed.columns:
            tqdm.write(f"[ERROR] sed_eflux not found for source {self.internal_id}")
            return

        e = self.sed['sed_eflux']
        self.sed = self.sed[(e > 0) & np.isfinite(e)]
        # check if no datapoint was left
        if len(self.sed) == 0:
            self.sed = None

    def __dropInvalidFrequency(self):
        """
        Sanity check to make sure alls ed_freq returned by VizieR SED are positive.
        """
        if self.sed is None:
            return

        if 'sed_freq' not in self.sed.columns:
            tqdm.write(f"Error! sed_freq missing for source {self.internal_id}.")
            self.sed = None
            return

        f = self.sed['sed_freq']
        self.sed = self.sed[np.isfinite(f) & (f > 0)]

        if len(self.sed) == 0:
            self.sed = None
    
    def __computeCoordDistance(self):
        """
        Estimate how far from central coordinate each measurement is.
        Add separation in arcseconds as 'dist_to_center' column in the SED file.
        """
        if self.sed is None:
            return
        # sanity check if VizieR SED is not bugging coordinates
        if '_RAJ2000' not in self.sed.columns or '_DEJ2000' not in self.sed.columns:
            return
        try:
            # gets separation
            dist = self.centerCoord.separation(SkyCoord(
                ra=self.sed['_RAJ2000'].values * u.deg,
                dec=self.sed['_DEJ2000'].values * u.deg,
                frame='icrs'
            )).arcsec
            self.sed['dist_to_center'] = dist

        except Exception as e:
            tqdm.write(f"Warn! Could not compute measurement coordinates "
                    f"for source {self.internal_id}: {e}")

    def __dropNegativeFlux(self):
        """
        Makes sure no bugged flux values came from VizieR

        """
        if self.sed is None:
            return

        if 'sed_flux' not in self.sed.columns:
            tqdm.write(f"Error! sed_flux missing for source {self.internal_id}")
            self.sed = None
            return

        f = self.sed['sed_flux']
        self.sed = self.sed[np.isfinite(f) & (f >= 0)]

        if len(self.sed) == 0:
            self.sed = None

    def __convertData(self):
        """
        VizieR SED provides 
        sed_freq [GHz] -> 1e9 u.Hz
        sed_flux [Jy] -> 1 u.Jy
        sed_eflux [Jy] -> 1 u.Jy,

        Here I use astropy units to convert these to 
        lambdas [micrometers]
        while keeping
        flux [Jy]
        eflux [Jy]

        At this tage, data is sorted by wavelength.
        """
        if self.sed is None:
            return

        try:
            freq = self.sed["sed_freq"].to_numpy() * u.GHz
            lamb = (const.c / freq).to(u.micrometer).value
            flux  = (self.sed['sed_flux'].to_numpy()  * u.Jansky).to(u.Jansky).value
            eflux = (self.sed['sed_eflux'].to_numpy() * u.Jansky).to(u.Jansky).value

            # Save converted columns
            self.sed['sed_lambda'] = lamb
            self.sed['sed_flux']   = flux
            self.sed['sed_eflux']  = eflux

            # Sort by wavelength increasing
            self.sed = self.sed.sort_values(by="sed_lambda", ascending=True)

        except Exception as e:
            tqdm.write(f"Error! Unit conversion failed for source {self.internal_id}: {e}")
            self.sed = None

    def keep_wavelength_range(self, lower, upper):
        """
        Keep only SED data within a given wavelength range.

        All measurements with wavelengths outside the interval
        ``[lower, upper]`` are discarded.

        Parameters
        ----------
        lower : float
            Lower bound of the wavelength range (microns).
        upper : float
            Upper bound of the wavelength range (microns).

        Notes
        -----
        - Wavelengths are assumed to be stored in the column ``sed_lambda``
        in units of microns.
        - This operation is irreversible: data outside the specified range
        are permanently removed from the SED.
        """
        if self.sed is None:
            return

        if lower > upper:
            raise ValueError("lower must be <= upper")

        mask = (
            (self.sed["sed_lambda"] >= lower) &
            (self.sed["sed_lambda"] <= upper)
        )

        self.sed = self.sed[mask].reset_index(drop=True)

    def drop_wavelength_range(self, lower, upper):
        """
        Remove SED data within a given wavelength range.

        All measurements with wavelengths inside the interval
        ``[lower, upper]`` are discarded; data outside the range
        are preserved.

        Parameters
        ----------
        lower : float
            Lower bound of the wavelength interval to remove (microns).
        upper : float
            Upper bound of the wavelength interval to remove (microns).

        Notes
        -----
        - Wavelengths are assumed to be stored in the column ``sed_lambda``
        in units of microns.
        - This is useful for masking known problematic wavelength regions
        (e.g. strong atmospheric bands or instrumental artifacts).
        """
        if self.sed is None:
            return

        if lower > upper:
            raise ValueError("lower must be <= upper")

        mask = (
            (self.sed["sed_lambda"] < lower) |
            (self.sed["sed_lambda"] > upper)
        )

        self.sed = self.sed[mask].reset_index(drop=True)

    def confine(self, r):
        """
        Reduce the searching radius to a radius r (arcseconds)
        around the central coordinate.
        """
        self.r_confine = r

        if self.sed is None:
            return
        if "dist_to_center" not in self.sed.columns:
            self.__computeCoordDistance()
            if self.sed is None or "dist_to_center" not in self.sed.columns:
                return


    def dropCatalog(self, tabname):
        """
        Drop rows belonging to a specific VizieR table (_tabname).
        """
        if self.sed is None:
            return
        if '_tabname' not in self.sed.columns:
            return

        self.sed = self.sed[self.sed['_tabname'] != tabname].reset_index(drop=True)

    def dropFilter(self, sed_filter):
        """
        Drop rows with a specific sed_filter.
        """
        if self.sed is None:
            return
        if 'sed_filter' not in self.sed.columns:
            return
        # Accept single string or a collection
        if isinstance(tabname, str):
            mask = self.sed["_tabname"] != tabname
        else:
            # tabname is iterable: drop any in the list/set
            tabset = {str(t) for t in tabname if t is not None}
            if not tabset:
                return
            mask = ~self.sed["_tabname"].astype(str).isin(tabset)

        self.sed = self.sed.loc[mask].reset_index(drop=True)

    def mergeRepeatedPoints(self, tol_arcsec=0.01, flux_tol_frac=1e-5):
        """
        Identify and merge repeated measurements:
        To be considered the same, they must have the same filter name (sed_filter),
        nearly identical RA/DEC (down to a tolerance tol_arcsec), 
        nearly identical flux (down to a tolerance flux_tol_frac).

        Groups of repeated measurements are merged into a single row with:
        - flux: inverse-variance weighted mean
        - flux uncertainty: propagated from inverse-variance weighting
        - metadata row: the measurement with the smallest reported flux uncertainty

        All contributing survey catalog names (``_tabname``) are stored in the
        output column ``aggregateTabName``.


        Parameters
        ----------
        tol_arcsec : float.
        flux_tol_frac : float
            Relative flux tolerance: |f1 - f2| < flux_tol_frac * max(f1,f2)
        """
        if self.sed is None or len(self.sed) == 0:
            return

        df = self.sed.copy()

        # --- bin size in degrees ---
        bin_size = tol_arcsec / 3600.0

        # --- build bins ---
        df["_ra_bin"]  = np.floor(df["_RAJ2000"] / bin_size).astype(int)
        df["_dec_bin"] = np.floor(df["_DEJ2000"] / bin_size).astype(int)

        merged_rows = []

        # --- group first by filter, then by sky bin ---
        for (_, ra_bin, dec_bin), g in df.groupby(
            ["sed_filter", "_ra_bin", "_dec_bin"]
        ):
            if len(g) == 1:
                row = g.iloc[0].copy()
                row["aggregateTabName"] = (
                    g["_tabname"].iloc[0] if "_tabname" in g.columns else ""
                )
                merged_rows.append(row)
                continue

            # small group → safe to do O(k²)
            coords = SkyCoord(
                ra=g["_RAJ2000"].values * u.deg,
                dec=g["_DEJ2000"].values * u.deg,
            )

            used = np.zeros(len(g), dtype=bool)

            for i in range(len(g)):
                if used[i]:
                    continue

                idx = [i]

                for j in range(i + 1, len(g)):
                    if used[j]:
                        continue

                    # separation check
                    if coords[i].separation(coords[j]).arcsec > tol_arcsec:
                        continue

                    # flux similarity
                    f1 = g.iloc[i]["sed_flux"]
                    f2 = g.iloc[j]["sed_flux"]
                    if abs(f1 - f2) > flux_tol_frac * max(abs(f1), abs(f2)):
                        continue

                    idx.append(j)
                    used[j] = True

                group = g.iloc[idx]

                # --- weighted merge ---
                f = group["sed_flux"].to_numpy()
                e = group["sed_eflux"].to_numpy()
                w = 1.0 / e**2

                fbar = np.sum(f * w) / np.sum(w)
                ebar = np.sqrt(1.0 / np.sum(w))

                rep = group.iloc[np.argmin(e)].copy()
                rep["sed_flux"]  = fbar
                rep["sed_eflux"] = ebar

                if "_tabname" in group.columns:
                    rep["aggregateTabName"] = ",".join(group["_tabname"].unique())
                else:
                    rep["aggregateTabName"] = ""

                merged_rows.append(rep)

        new_df = pd.DataFrame(merged_rows).reset_index(drop=True)

        # cleanup
        new_df = new_df.drop(columns=["_ra_bin", "_dec_bin"], errors="ignore")

        self.sed = new_df

    def aggregateFilter(self):
        """
        Aggregate all measurements that have the same filter.

        This was modified from the original version, where smae filters
        where aggregated by the simple mean. 

        Assumes that '_mergeRepeatedPoints()' was already applied and that
        no true duplicates remain.

        For each filter:
        - compute weighted mean flux and propagated uncertainty
        - pick representative row = smallest error
        - store all contributing _tabname values in 'aggregateTabName'
        - return one row per sed_filter
        """

        if self.sed is None or len(self.sed) == 0:
            return
        if 'sed_filter' not in self.sed.columns:
            return

        df = self.sed
        aggregated = []

        for filt, g in df.groupby('sed_filter'):

            # if only one datapoint exists → keep it unchanged
            if len(g) == 1:
                row = g.iloc[0].copy()

                # add aggregateTabName for consistency
                if '_tabname' in g.columns:
                    row['aggregateTabName'] = g['_tabname'].iloc[0]
                else:
                    row['aggregateTabName'] = ""

                aggregated.append(row)
                continue

            # Weighted mean flux
            f = g['sed_flux'].to_numpy()
            e = g['sed_eflux'].to_numpy()
            w = 1.0 / e**2

            fbar = np.sum(f * w) / np.sum(w)
            ebar = np.sqrt(1.0 / np.sum(w))

            # Representative row = smallest error
            rep = g.iloc[np.argmin(e)].copy()
            rep['sed_flux']  = fbar
            rep['sed_eflux'] = ebar

            # Collect all contributing catalogs
            if '_tabname' in g.columns:
                rep['aggregateTabName'] = ",".join(g['_tabname'].unique())
            else:
                rep['aggregateTabName'] = ""

            aggregated.append(rep)

        # Build final dataframe
        self.sed = pd.DataFrame(aggregated).reset_index(drop=True)

    def resolveMultipleDetections(self):
        """
        Detect and resolve multiple detections from the same survey (_tabname)
        in the same filter (sed_filter).

        If multiple rows have the same (_tabname, sed_filter), the one closest
        to the central coordinate is kept, and it receives a flag 'mult=True'.
        All others are removed.

        This is useful to clean sources with possible visual pair and flag them. 
        """

        if self.sed is None:
            return

        required = ['_tabname', 'sed_filter', '_RAJ2000', '_DEJ2000']
        if not all(c in self.sed.columns for c in required):
            return

        df = self.sed.copy()

        # Build a SkyCoord array for separation computation
        coords = SkyCoord(
            ra=df['_RAJ2000'].values * u.deg,
            dec=df['_DEJ2000'].values * u.deg,
            frame='icrs'
        )

        # Initialize 'mult' flag
        if 'mult' not in df.columns:
            df['mult'] = False

        to_remove = []

        # Group by (tabname, filter)
        for (tab, filt), g in df.groupby(['_tabname', 'sed_filter']):

            if len(g) <= 1:
                continue  # nothing to resolve

            # Compute separations from central coordinate
            group_coords = SkyCoord(
                ra=g['_RAJ2000'].values * u.deg,
                dec=g['_DEJ2000'].values * u.deg,
                frame='icrs'
            )
            sep = self.centerCoord.separation(group_coords).arcsec

            # Select the closest measurement
            keep_idx = g.index[np.argmin(sep)]

            # Flag the kept one
            df.loc[keep_idx, 'mult'] = True

            # Mark others for removal
            for idx in g.index:
                if idx != keep_idx:
                    to_remove.append(idx)

        # Remove non-closest detections
        df = df.drop(to_remove).reset_index(drop=True)

        self.sed = df

    def get_plot_arrays(self, *, drop_nonfinite=True):
        """
        At any given state of Source, this will return the current version of the SED relevant fields:
        sed_lambda, sed_flux, sed_eflux <- returns the three as numpy arrays.

        """
        if self.sed is None or len(self.sed) == 0:
            return np.array([]), np.array([]), np.array([])

        df = self.sed

        # Required columns
        required = ("sed_lambda", "sed_flux", "sed_eflux")
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise KeyError(f"Missing required SED columns: {missing}")

        lam = df["sed_lambda"].to_numpy(dtype=float)
        flux = df["sed_flux"].to_numpy(dtype=float)
        eflux = df["sed_eflux"].to_numpy(dtype=float)

        if drop_nonfinite:
            m = (
                np.isfinite(lam) & (lam > 0) &
                np.isfinite(flux) & (flux >= 0) &
                np.isfinite(eflux) & (eflux > 0)
            )
            lam, flux, eflux = lam[m], flux[m], eflux[m]

        # Sort by wavelength
        if lam.size > 1:
            idx = np.argsort(lam)
            lam, flux, eflux = lam[idx], flux[idx], eflux[idx]

        return lam, flux, eflux

    def plot_sed(self, ax=None, *, use_confine=True, y="flux", show=True, **kwargs):
        """
        Makes a quick plot of a SED

        Parameters
        ----------
        y : {"flux", "nuFnu"}
            flux -> plot sed_flux (Jy)
            nuFnu -> plot nu * Fnu
        """

        lam, flux, eflux = self.get_plot_arrays(use_confine=use_confine)

        if ax is None:
            fig, ax = plt.subplots(figsize=(6.5, 4.0), dpi=150)

        if lam.size == 0:
            ax.text(0.5, 0.5, "No SED data", ha="center", va="center", transform=ax.transAxes)
            if show:
                plt.show()
            return ax

        if y == "nuFnu":
            # nu = c / lambda
            nu = (c_light.value) / (lam * 1e-6)  # Hz (lam in um)
            yy = nu * flux
            yerr = nu * eflux if np.all(np.isfinite(eflux)) else None
            ax.set_ylabel(r"$\nu F_\nu$ (Hz·Jy)")
        else:
            yy = flux
            yerr = eflux if np.all(np.isfinite(eflux)) else None
            ax.set_ylabel(r"$F_\nu$ (Jy)")

        ax.errorbar(lam, yy, yerr=yerr, fmt="o", ms=4, lw=1, capsize=2, **kwargs)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\lambda\ (\mu m)$")

        title = str(self.internal_id) if self.internal_id is not None else f"{self.ra:.5f}, {self.dec:.5f}"
        ax.set_title(title)

        if show:
            plt.show()
        return ax