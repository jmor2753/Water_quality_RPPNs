# Curitiba water chemistry | R replication
Use R 4.5 or newer. Open the extracted R_release folder as the working directory in RStudio. Install the dependencies once:
```r
install.packages(c("vegan", "xml2", "zip", "ggplot2", "ggrepel", "patchwork", "sf", "geobr"))
source("replicate.R", encoding="UTF-8")
```
Command-line equivalent after installation: Rscript replicate.R

The script performs preprocessing, PCA, site PERMANOVA, dispersion testing, the seven sensitivity scenarios, sample- and site-level Spearman correlations, landscape regressions and RDA. It generates Tables 1, 2 and S2, all three figures (PDF and PNG), sessionInfo.txt and the tracked Word manuscript under results/. Statistics and figures are implemented entirely in R. xml2 and zip are used only to assemble the Word output. The manuscript input is a review template containing the original figure artwork as layout references and values from the earlier numerical cross-check. Use the regenerated Word output for the R version. This script has not been executed in the preparation environment because that environment had no R installation.

data/chemistry.csv contains the 36 original sample records and units, including censoring flags and the negative blank-corrected Al measurement. data/sites.csv contains the 12 representative coordinates and the land-cover percentages for 600 m buffers. All analyses start from those two files. The two compressed integer matrices fix the 9,999 permutation draws and are read directly by R. The chemistry test uses sites in alphabetical order and samples in their original within-site order. Landscape tests use site-file order. The land-cover predictors are supplied as extracted from the original spatial workflow; the MapBiomas raster is not required for the models or site map.

Values below LOQ are assigned 0.001 µg/L in the primary analysis. Sensitivity substitutes 0, 25%, 50% and 100% of each element's observed minimum quantified concentration, then separately omits Cr, Zn and both. The empirical substitution levels are not estimates of laboratory limits.

## Figure design
The design_reference directory contains the three unmodified PNGs extracted from the submitted manuscript. figures.R follows that manuscript layout. Figure 1 uses pale neighbourhood polygons and a dark municipal outline, a bottom legend, labelled points, Marco Zero, and 600 m circles. Figure 2 pairs the PCA and correlation heatmap. Figure 3 has two labelled regression panels and a shared site-type legend. The plotting script uses the original colours and orients PC2 with positive Zn loading, matching the submitted biplot. Signed slopes and plotted site scores follow the same convention. Figures are regenerated from the revised data. Exact label positions may change. On first run, Figure 1 downloads the two geobr boundary layers and caches them in data/curitiba_boundaries.rds. Internet access is required for that first fetch. The original PNGs are design references and retain their original analytical inputs.
