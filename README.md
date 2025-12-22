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

Next, a list of extra processing steps is available:


References:

"""