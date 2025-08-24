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

>Note: Large input and output files (e.g., occurrence maps, SHAP matrices) are not tracked in this repo but can be shared via Zenodo.

## Start

Please check documents: `docs/workflow.Rmd` to start the analysis.

## Acknowledgement

This work is supported by project “BoCP-Implementation: BioFI- Biodiversity Forecasting Initiative to Understand Population, Community and Ecosystem Function Under Global Change” funded by a grant to A.E.F. from the National Science Foundation (award number: 2416164).