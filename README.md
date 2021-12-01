# gmc_spawning_migration

This repository contains all the code required to analyse acoustic telemetry data from Giant Mud Crab (*Scylla serrata*) tagged in the Clarence and Kalang River. Data can be accessed via the [IMOS animal tracking database](https://animaltracking.aodn.org.au/).

This work is currently being reviewed for publication in [Estuaries and Coasts](https://www.springer.com/journal/12237).

The following briefly describes the workflow:
- `rsp.R` is the script used to recreate the most likely trajectories of tagged crabs as they moved throughout the estuary (using the R package 'RSP'), this process is fully described in [Niella et al., 2020](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13484), while the [RSP github](https://github.com/YuriNiella/RSP) provides some useful vignettes.
- `migration_metrics.R` is the main script, in this some data wrangling and the main analysis and plotting is done.
- `oceanic_crabs.R` creates a summary of the compiled oceanic detections obtained from the [IMOS animal tracking database](https://animaltracking.aodn.org.au/) and receivers maintained under the [NSW Shark Management Strategy](https://www.sharksmart.nsw.gov.au/).
- site_maps.R` creates the site map for the study (Figure 1).
- `wq_plots.R` produces Supplementary Figure 1, which is a timeseries of all the environmental data collected during the project.
