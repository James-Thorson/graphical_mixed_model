# Code and data to reproduce results from Thorson (2025)

Thorson, J. T. (2025). Generalized graphical mixed models connect ecological theory with widely used statistical models. EcoEvoRxiv. https://doi.org/10.32942/X2963M

# Contents

### Primary scripts
* `abundance_at_age.R`:  R script to reproduce proportional abundance-at-age prediction using path coefficients along cohort, age, and/or year
* `diffusion_enhanced.R`:  R script to reproduce diffusion-enhanced spatio-temporal model on a gridded discretization of 2D space, using path coefficients along time, space, and space-time
* `mammal_traits.R`:  R script to reproduce to reproduce phylogenetic structural equation model (PSEM) for mammal body size, metabolic rate and range size while estimating different stabilizing selection for each trait
* `Fig_1.R`:  R script to produce visual illustration of path matrices

### Helper scripts
* `trait_functions.R`:  R functions sourced by `mammal_traits.R`, using RTMB to implement the PSEM
* `make_dsem_ram.R`:  R functions sourced by `diffusion_enhanced.R`, parsing arrow-and-lag notation to construct path matrices

### Results
* `2025-05-13`:  Directory containing all output from primary scripts (see above), including all figures from submission

### Data
Directory `data` containing:
* `GOA_Rex_8_2021_run7.dat`:  The Stock Synthesis data file for Gulf of Alaska rex sole ([here](https://github.com/noaa-afsc/goa_rex/blob/main/runs/2025_cie_review/2021_accepted_model_inputs/)), where lines 206-222 contain proportional abundance at age for ages (columns) and years (rows) of available sampling, and see `abundance_at_age.R` for further extraction and processing;
* `PanTHERIA_1-0_WR05_Aug2008.txt`: Mammal species traits downloaded from PanTHERIA ([here](https://esapubs.org/archive/ecol/E090/184/metadata.htm))
* `VertTree_mammals.tre`:  A dated phylogeny for mammals downloaded from VertLife ([here](https://vertlife.org/phylosubsets/))
