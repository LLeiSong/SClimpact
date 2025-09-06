# SClimpact

Explainable Artificial Intelligence Reveals Spatially Divergent Effects of Global Change on Mammals.

## Overview

This study uses SHapley Additive exPlanations (SHAP) to assess how climate and land cover drivers influence future mammal distributions across different global change scenarios. We quantify both directional changes (favoring/disfavoring turnovers) and gradual changes (favorability shifts) to evaluate species-specific spatial risks and opportunities under global change.

The repository includes:

- Scripts for preprocessing species and environmental data
- Code for fitting species distribution models (SDMs)
- SHAP-based analysis and visualization
- Reproducible figures and tables

## Repository structure

```text
├── docs/                    # Workflows and technical reports.
│   ├── figures/             # Scripts to generate manuscript figures
├── R/                       # All analysis scripts
│   ├── figures/             # Scripts to generate manuscript figures
├── Schedulers/              # Supportive HPC schedule scripts for R scripts.
├── data/
│   └── generation_length/   # Generation length data of species
│   └── IUCN/                # Species IUCN range polygons
│   └── IUCN_AOH_100m/       # The refined AOH of species
│   └── occurrences/         # Species occurrence data for modeling
│   └── terr-ecoregions-TNC/ # TNC ecoregion polygons used in the analysis
│   └── variables            # Environmental variables used for modeling
├── results/
│   ├── sdm/                 # All modeling results
│   └── species_analysis/    # Species-level analysis
│   └── climate_change/      # Global-level analysis aggregating all species
│   └── figures/             # All visual figures
│   └── tables/              # All summarized tables to support results
└── README.md
```

>Note: Large input and output files (e.g., occurrence maps, SHAP matrices) are not tracked in this repo but are shared via [OSF](https://doi.org/10.17605/OSF.IO/7MDJR).

## Demo

Please check the [demo](https://htmlpreview.github.io/?https://github.com/LLeiSong/SClimpact/blob/main/demo/demo.html) for the method used in this study.

## Start

Please check documents: [docs/workflow.Rmd](docs/workflow.Rmd) or its [online page](https://htmlpreview.github.io/?https://github.com/LLeiSong/SClimpact/blob/main/docs/workflow.html) to start the analysis.

## Session info

All analyses were conducted in R. The full R session information, including package versions and system details, is saved in [session_info.txt](session_info.txt) for reproducibility.

*NOTE:* due to a bug in CRAN version of package spatialEco, users need to install the [development version](https://github.com/jeffreyevans/spatialEco) of this package.

## Acknowledgement

This work is supported by project “BoCP-Implementation: BioFI- Biodiversity Forecasting Initiative to Understand Population, Community and Ecosystem Function Under Global Change” funded by a grant to A.E.F. from the National Science Foundation (award number: 2416164).
