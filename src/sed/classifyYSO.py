"""
@juliaroquette 22 December 2025:

This is an adaptation of David Hernandez code: https://github.com/Starvexx/standard_classification
YSOClassifier has now been tailored to be fed directly a Source instance with the data received from VizieR-SED.

    - "sed_lambda"  (micrometers)
    - "sed_flux"    (Jy)   = F_nu
    - "sed_eflux"   (Jy)   = sigma(F_nu)

Alternatively, if one wants to use external data, the source_from_arrays function can be used to mimic the Source class.
src = source_from_arrays(lam, flux, eflux, ID="V1547Ori")

YSOClassifier can estimate YSO infrared spectral index (alpha_IR). It assumes:
-  data is provided as F_nu in Jy at given wavelengths in micrometers.
- converts F_nu (Jy) -> nu*F_nu (erg/s/cm^2) using astropy units
- fits log10(nuFnu) vs log10(lambda) with a 2-step procedure 
- returns alpha_ir, intercept, wl_min/max, n_points, class label
- allows custom classification schemes using a dictionary with rules
- can plot the SED + fit

"""

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.constants import c

from scipy import optimize
from scipy.optimize import Bounds


#
# Using classification rules adopted in Roquette et al. 2025, A&A, 702, A63 Table 1. 
#
DEFAULT_RULES = {
    "0/I": lambda a: np.isfinite(a) and (a > 0.3),
    "flat": lambda a: np.isfinite(a) and (-0.3 < a < 0.3),
    "II": lambda a: np.isfinite(a) and (-1.6 < a < -0.3),
    "III thin disk": lambda a: np.isfinite(a) and (-2.5 < a < -1.6),
    "III no disk / MS": lambda a: np.isfinite(a) and (a < -2.5),
}


def apply_rules(alpha, rules):
    for label, pred in rules.items():
        if pred(alpha):
            return label
    return "not classified"


# -----------------------------------------------------------------------------
# Unit-safe conversion: 
# -----------------------------------------------------------------------------
def jy_to_nuFnu(lambda_um, fnu_jy):
    """ 
    Converts flux density into flux
    Jy (F_nu) -> nu*F_nu (erg/s/cm^2)
    """
    lam = np.asarray(lambda_um, dtype=float)
    fnu = np.asarray(fnu_jy, dtype=float)

    out = np.full(lam.shape, np.nan) * (u.erg / (u.s * u.cm**2))

    good = np.isfinite(lam) & np.isfinite(fnu) & (lam > 0)
    if np.any(good):
        lam_q = lam[good] * u.um
        fnu_q = fnu[good] * u.Jy

        nu = (c / lam_q).to(u.Hz)
        nuFnu = (fnu_q.to(u.erg / (u.s * u.cm**2 * u.Hz)) * nu).to(u.erg / (u.s * u.cm**2))
        out[good] = nuFnu

    return out


def jyerr_to_nuFnu_err(lambda_um, efnu_jy):
    """
    propagates flux density errors to nu*F_nu errors
    """
    lam = np.asarray(lambda_um, dtype=float)
    efnu = np.asarray(efnu_jy, dtype=float)

    out = np.full(lam.shape, np.nan) * (u.erg / (u.s * u.cm**2))

    good = np.isfinite(lam) & np.isfinite(efnu) & (lam > 0)
    if np.any(good):
        lam_q = lam[good] * u.um
        efnu_q = efnu[good] * u.Jy

        nu = (c / lam_q).to(u.Hz)
        enuFnu = (efnu_q.to(u.erg / (u.s * u.cm**2 * u.Hz)) * nu).to(u.erg / (u.s * u.cm**2))
        out[good] = enuFnu

    return out


def _line(x, k, d):
    """
    Simple Linear model to be used with curve_fit
    """
    return k * x + d


def _powerlaw(lam_um, k, d):
    """
    Power-law model to be used with curve_fit
    """
    lam_um = np.asarray(lam_um, dtype=float)
    return 10.0 ** (k * np.log10(lam_um) + d)


def fit_alpha(lambda_um, nuFnu, enuFnu, lower=2.0, upper=24.0, min_points=2, wl_span_threshold_um=2.0):
    """
    Fit the infrared spectral index alpha_IR using a 2-step fitting procedure. 
    Fits first a linear model in log-log space, then refines with a power-law fit. 
    (Makes sure things won't go crazy with bounds from the first fit.)
    """
    lam = np.asarray(lambda_um, dtype=float)
    f = np.asarray(nuFnu, dtype=float)
    ef = np.asarray(enuFnu, dtype=float)

    mask = (
        np.isfinite(lam) & np.isfinite(f) & np.isfinite(ef) &
        (lam > 0) & (f > 0) &
        (lam > lower) & (lam < upper)
    )

    n = int(np.count_nonzero(mask))
    if n < int(min_points):
        return np.nan, np.nan, np.nan, np.nan, n

    wl_min = float(np.min(lam[mask]))
    wl_max = float(np.max(lam[mask]))

    if (wl_max - wl_min) < float(wl_span_threshold_um):
        return np.nan, np.nan, np.nan, np.nan, n

    x = np.log10(lam[mask])
    y = np.log10(f[mask])

    try:
        popt, pcov = optimize.curve_fit(_line, x, y)

        if (pcov is None) or (pcov.shape != (2, 2)):
            se_k, se_d = 5.0, 5.0
        else:
            se_k = float(np.sqrt(pcov[0, 0])) if np.isfinite(pcov[0, 0]) else 5.0
            se_d = float(np.sqrt(pcov[1, 1])) if np.isfinite(pcov[1, 1]) else 5.0

        bounds = Bounds(
            (popt[0] - se_k, popt[1] - se_d),
            (popt[0] + se_k, popt[1] + se_d),
        )

        popt2, _ = optimize.curve_fit(_powerlaw, lam[mask], f[mask], bounds=bounds)

        alpha = float(popt2[0])
        intercept = float(popt2[1])

        return alpha, intercept, wl_min, wl_max, n

    except Exception:
        return np.nan, np.nan, np.nan, np.nan, n


class YSOClassifier(object):
    """
    Class to classify YSOs based on their infrared spectral index (alpha_IR).
    """
    def __init__(self, source, # should be a vizierSED Source instance 
                 lower=2.0, # minimum wavelength in micrometers for fit
                 upper=24.0, # maximum wavelength in micrometers for fit
                 min_points=2, # minimum number of points required for fit
                 wl_span_threshold_um=2.0, # minimum wavelength span required for fit
                 rules=None, # classification rules to apply -  If not provided, use DEFAULT_RULES
                 ):
        self.source = source
        self.lower = float(lower)
        self.upper = float(upper)
        self.min_points = int(min_points)
        self.wl_span_threshold_um = float(wl_span_threshold_um)
        self.rules = DEFAULT_RULES if rules is None else rules

        if source.sed is None:
            raise ValueError("YSOClassifier: source.sed is None (no SED data).")

        # Source.sed is a pandas DataFrame with columns:
        #   sed_lambda [um], sed_flux [Jy], sed_eflux [Jy]
        self.lambda_um = source.sed["sed_lambda"].to_numpy(dtype=float)
        self.fnu_jy = source.sed["sed_flux"].to_numpy(dtype=float)
        self.efnu_jy = source.sed["sed_eflux"].to_numpy(dtype=float)

        # unit-safe conversion to nuFnu (Quantity)
        self.nuFnu_q = jy_to_nuFnu(self.lambda_um, self.fnu_jy)
        self.enuFnu_q = jyerr_to_nuFnu_err(self.lambda_um, self.efnu_jy)

        # numeric arrays for scipy
        self.nuFnu = self.nuFnu_q.value
        self.enuFnu = self.enuFnu_q.value

        # outputs
        self.alpha_ir = np.nan
        self.intercept = np.nan
        self.wl_min = np.nan
        self.wl_max = np.nan
        self.n_points = 0
        self.yso_class = "not classified"

    def run(self):
        a, d, wlmin, wlmax, n = fit_alpha(
            self.lambda_um,
            self.nuFnu,
            self.enuFnu,
            lower=self.lower,
            upper=self.upper,
            min_points=self.min_points,
            wl_span_threshold_um=self.wl_span_threshold_um,
        )

        self.alpha_ir = a
        self.intercept = d
        self.wl_min = wlmin
        self.wl_max = wlmax
        self.n_points = int(n)
        self.yso_class = apply_rules(self.alpha_ir, self.rules)

        return {
            "alpha_ir": self.alpha_ir,
            "intercept": self.intercept,
            "wl_min": self.wl_min,
            "wl_max": self.wl_max,
            "n_points": self.n_points,
            "class": self.yso_class,
        }

    def plot(self, ax=None, savepath=None, dpi=150, show_fit=True):
        if ax is None:
            fig, ax = plt.subplots(figsize=(4, 4), dpi=dpi)
        else:
            fig = ax.figure

        mask = (
            np.isfinite(self.lambda_um) & np.isfinite(self.nuFnu) & np.isfinite(self.enuFnu) &
            (self.lambda_um > 0) & (self.nuFnu > 0) &
            (self.lambda_um > self.lower) & (self.lambda_um < self.upper)
        )

        ax.errorbar(
            self.lambda_um, self.nuFnu, yerr=self.enuFnu,
            fmt=".", ms=3.5, elinewidth=0.75, capsize=2
        )

        if np.any(mask):
            ax.scatter(self.lambda_um[mask], self.nuFnu[mask], s=10)

        if show_fit and np.isfinite(self.alpha_ir) and np.isfinite(self.intercept):
            wl_grid = np.logspace(np.log10(0.5), np.log10(1000), 200)
            ax.plot(wl_grid, _powerlaw(wl_grid, self.alpha_ir, self.intercept), ls="--", lw=1)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\lambda\ (\mu m)$")
        ax.set_ylabel(r"$\nu F_{\nu}\ (erg\ s^{-1}\ cm^{-2})$")

        ax.set_title("alpha={0:.2f} | {1} | N={2}".format(self.alpha_ir, self.yso_class, self.n_points))

        ax.axvspan(0.5, self.lower, alpha=0.1)
        ax.axvspan(self.lower, self.upper, alpha=0.15)
        ax.axvspan(self.upper, 1e3, alpha=0.1)

        if savepath is not None:
            fig.savefig(savepath, dpi=dpi, bbox_inches="tight")

        return ax

class SimpleSource(object):
    """
    This is a class that allows to mimic a Source class with external input data. 
    It will recreate the minimal expected structure for a Source from vizierSED
    """
    def __init__(self, sed, ID=None, ra=None, dec=None):
        self.sed = sed
        self.ID = ID
        self.ra = ra
        self.dec = dec


def source_from_arrays(sed_lambda, sed_flux, sed_eflux, ID=None, ra=None, dec=None):
    """
    Build a minimal Source-like object from three arrays.

    Parameters
    ----------
    sed_lambda : array-like, micrometers
    sed_flux   : array-like, Jy
    sed_eflux  : array-like, Jy

    Returns
    -------
    SimpleSource instance mimicking a Source from vizierSED.py
    """
    lam = np.asarray(sed_lambda, dtype=float)
    flx = np.asarray(sed_flux, dtype=float)
    eflx = np.asarray(sed_eflux, dtype=float)

    if not (lam.shape == flx.shape == eflx.shape):
        raise ValueError("source_from_arrays: sed_lambda, sed_flux, sed_eflux must have the same shape.")

    sed = pd.DataFrame({
        "sed_lambda": lam,
        "sed_flux": flx,
        "sed_eflux": eflx,
    })

    return SimpleSource(sed, ID=ID, ra=ra, dec=dec)
