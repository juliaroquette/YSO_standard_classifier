# Get SEDs and classify them

**@juliaroquette 22 December 2025:** This little package put together a series of tools developed during the [Horizon 2020 project NEMESIS project](https://nemesis.konkoly.hu/) and allows to query the [VizieR database](https://vizier.cds.unistra.fr) for a given sky coordinate, download and process the relevant photometric data, and use this data to infer standard for Young Stellar Objects. 

The tools here were built on other tools put together by [David Hernandez](https://github.com/Starvexx) and Gabor Marton:
-  [SED_query](https://github.com/Starvexx/SED_query)
- [stanard_classification](https://github.com/Starvexx/standard_classification)

Other relevant references:
- VizieR SED service: https://vizier.cds.unistra.fr/vizier/sed/
- VizieR SED API documentation: https://vizier.cds.unistra.fr/vizier/sed/assets/sedapi.gml
- VizieR Meta-filter description: https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=METAfltr 
  Photometric filter names
- VizieR META-phot description: https://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=METAphot
  Photometric systems details (zero-points, etc)
- Roquette et al. 2025, A&A, 702, A63 https://www.aanda.org/articles/aa/full_html/2025/10/aa53588-24/aa53588-24.html

# Installation

To install, run from the package-root folder

```bash
python -m pip install -e .
```

# Quickstart (query one source)

Basic class `Source` queries VizieR SED for a given coordinate and search radius, retrieves the data,
processes it to remove invalid points, and saves it as a csv file.

```python
from sed.vizierSED import Source
```

For example, let's say we want to query the T Tauri star [V1547 Ori](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=V1547Ori&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id) using a 2" searching radius. 

```python
ra = 83.86937
dec = -4.79072
radius = 2.0  # arcseconds
ID = "V1547Ori"
```

If we want to simply export the raw output SED (no processing yet), to the following directory:
```python
raw_path = "./output/raw"  # where to save the raw downloaded SED table (optional)
```

Then:

```python
star = Source(
    ra,
    dec,
    radius,
    ID=ID,
    raw_cache_dir=raw_path,  # set to None if you don't want to save raw files
    force_download=True  
)
```


will save a file `./output/raw/V1547Ori.csv` with the SED data queried from VizieR SED, with basic processing to remove invalid points, and converted to wavelength in microns, and flux density in Jansky. Here `force_download=True` ensures that a raw SED will be downloaded from VizieR-SED even if a file with the same name already exists in `./output/raw/`. Note that the `ID` field here is passed as a labelling tool only. If only the coordinate is passed, name convention will follow `<RA>_<Dec>_r<radius>arcsec`, where the example above would be `83.86937_-4.79072_r2.00arcsec_processed.csv`

Next, a list of extra processing steps is available:


References:

"""