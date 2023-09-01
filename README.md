

# A repository for:

Fidino, M et al.  Gentrification drives patterns of alpha and beta diversity in cities

## Links to different parts of the readme file

1. [What's in this repository?](#whats-in-this-repository)
2. [The working directory](#the-working-directory)
3. [city codes](#city-codes)
4. [The data folder (`./data`)](#the-data-folder-data)
5. [The figures folder (`./figures`)](#the-figures-folder-figures)
6. [The JAGS folder (`./JAGS`)](#the-jags-folder-jags)
7. [The mcmc output folder (`./mcmc_output`)](#the-mcmc-output-folder-mcmc_output)
8. [The R folder (`./R`)](#the-r-folder-r)


## What's in this repository?

This repository stores all of the data and code used to query U.S. Census data,
fit a multi-city, multi-species occupancy model, estimate alpha and beta diversity at sites with uncertainty, and then finally fit separate alpha and beta
diversity meta-analytic models. The folder organization seperates the data (`./data`),figures from the manuscript (`./figures`), JAGS model (`./JAGS`), the mcmc outputs from the model we fitted to the data (`./mcmc_outputs`), the  R code (`./R`), and the markdown file used to create the supplemental material (`./supplemental`).

This document here serves as a road map that describes all of the files present in this repository.


[Back to table of contents ⤒](#a-repository-for)


### The working directory

Aside from the aforementioned folders, the working directory here stores the `.gitignore` file for this repository, this README file (`README.md`) and the `.Rproj` file (for if you are using RStudio, `uwin-gentrification.Rproj`).

[Back to table of contents ⤒](#a-repository-for)

### City codes

When there is data specific to a given city we use 4 letter abbreviations.
These abbreviations are:

| Code  |       City         |
|-------|--------------------|
|`ahga` | Athens, GA         |
|`boma` | Bay Area, CA       |
|`chil` | Boston, MA         |
|`deco` | Chicago, IL        |
|`hote` | Denver, CO         |
|`inin` | Houston, TX        |
|`ioio` | Indianapolis, IN   |
|`jams` | Iowa City, IA      |
|`lrar` | Jackson, MS        |
|`mawi` | Little Rock, AR    |
|`mela` | Madison, WI        |
|`naca` | Metro LA, CA       |
|`baca` | Washington D.C.    |
|`phaz` | Phoenix, AZ        |
|`poor` | Portland, OR       |
|`rony` | Rochester, NY      |
|`safl` | Sanford, FL        |
|`scut` | Salt Lake City, UT |
|`sewa` | Seattle, WA        |
|`slmo` | Saint Louis, MO    |
|`tawa` | Tacoma, WA         |
|`uril` | Urbana, IL         |
|`wide` | Wilmington, DE     |

[Back to table of contents ⤒](#a-repository-for)

### The data folder (`./data`)

This folder has 4 files and 2 sub-folders.

|**`./data/detection_data.csv`** | The camera trap data used in our analysis.
This csv file has 576,000 rows and 9 columns. 

[Back to table of contents ⤒](#a-repository-for)




### The figures folder (`./figures`)

This folder houses some of the raw figures I generated in `R` (which I cleaned up using Inkscape), as well as other figures that were publication ready. All of the extra "cleaning" I needed to do was related to the maps I had made (there was too much spacing between images).


|File|Description|
|---|---|
|**`./figures/figure_1.pdf`** | Figure 1 in the manuscript, publication ready, which is an example of the modeling framework.|
|**`./figures/figure_1.png`** | Figure 1 as a png so I could use it in this README.md.|
|**`./figures/figure_1.svg`** | Figure 1 as a scaleable vector graphic.|
|**`./figures/figure_2.pdf`** | Figure 2 in the manuscript, publication ready, which shows each species distribution, their conflict potential, and where they are most likely to come in to conflict with humans.|
|**`./figures/figure_2.svg`** | Figure 2 as a scaleable vector graphic.|
|**`./figures/figure_3.tiff`** | Figure 3 in the manuscript, publication ready, which shows the slope terms from the model.|
|**`./figures/figure_4.tiff`** | Figure 4 in the manuscript, publication ready, which shows the conflict potential response for different covariates.|
|**`./figures/rough_1.svg`** | Output from R that was used to generate figure 1.|
|**`./figures/spatial_variables.jpeg`** | Plot of the spatial variables used in the model, used in this README.md.|
|**`./figures/supl_coyote.png`** | Supplemental figure for coyote, which shows their spatiotemporal correlation in occupancy across the 12 seasons of sampling.|
|**`./figures/supl_coyote.svg`** | Same, but as a scaleable vector graphic (output from R).|
|**`./figures/supl_opossum.png`** | Supplemental figure for Virginia opossum, which shows their spatiotemporal correlation in occupancy across the 12 seasons of sampling.|
|**`./figures/supl_opossum.svg`** | Same, but as a scaleable vector graphic (output from R).|
|**`./figures/supl_raccoon.png`** | Supplemental figure for raccoon, which shows their spatiotemporal correlation in occupancy across the 12 seasons of sampling.|
|**`./figures/supl_raccoon.svg`** | Same, but as a scaleable vector graphic (output from R).|
|**`./figures/supl_urb_map.png`** | Supplemental figure for the two urban intensity metrics and the variables that were used to create them.|
|**`./figures/supl_urb_map.svg`** | Same, but as a scaleable vector graphic (output from R).|


[Back to table of contents ⤒](#a-repository-for)

### The JAGS folder (`./JAGS`)

This folder houses one script, which is the `JAGS` model that is fit to the data, which is titled `dynamic_integrated_occupancy_gam.R`. The code is commented out within the model pretty well (I had to write it in a way that was more difficult to read because it ran far faster that way).

[Back to table of contents ⤒](#a-repository-for)

### The mcmc output folder (`./mcmc_output`)

This folder houses the posterior distributions from the models we fit to the data of coyote (`./mcmc_output/coyote_model.RDS`), opossum (`./mcmc_output/opossum_model.RDS`), and  raccoon (`./mcmc_output/raccoon_model.RDS`), which are stored as RDS files. I used `run.jags` to fit the models, so they are `runjags` objects. If you want to grab the posterior distribution from the models it can be done as:

```R
my_model <- readRDS("./mcmc_output/coyote_model.RDS)

my_mcmc <- do.call("rbind", my_model$mcmc)
```

This folder also has a summary of all the parameters as a pdf (which was the supplemental material for the manuscript). There is an R markdown file (`./mcmc_output/Appendix_S1.Rmd`) and the associated PDF that was knitted from it (`./mcmc_output/Appendix_S1.pdf`).

#### The diagnostic plots sub-folder (`./mcmc_output/diagnostic_plots)

This folder has three sub-folders (one for coyote, one for opossum, and one for raccoon). Inside of each of these are the traceplots for each model parameter. They are generated when `fit_models.R` is run.  See `./JAGS/dynamic_integrated_occupancy_gam.R` for where each parameter fits into the model.


#### The validation sub-folder (`./mcmc_output/validation`)

This sub-folder contains the posterior distributions from the models fit to only a portion of the human-wildife conflict data so that we could do some model validation for the coyote (`./mcmc_output/valdiation/coyote_validation_model.RDS`), opossum (`./mcmc_output/valdiation/opossum_validation_model.RDS`), and raccoon (`./mcmc_output/valdiation/raccoon_validation_model.RDS`). There are two other RDS files, which are created via `./R/validate_model.R`, and include `ROC_scores.RDS` and `model_auc.RDS`. The former is a list object of length three (one element for each species) that contains the ROC scores for each potential threshold and sampling period. The latter is also a list object of length three (one element for each species) that summarises the AUC scores. Each element is also a list and contains the global AUC across the whole study (e.g., `model_auc$coyote$auc_global`), for each sampling period (e.g.,`model_auc$coyote$auc_year`), as well as the true positive rate (tpr), false positive rate (fpr), and accuracy (acc) for each threshold and sampling period. We used 41 threshold values, evenly spaced between 0 and 1.

[Back to table of contents ⤒](#a-repository-for)


### The parameter summary folder (`./parameter_summary`)

This folder houses three files that contain a summary of of parameters estimated by the model for coyote (`./parameter_summary/coyote_parameter_estimates.csv`), opossum (`./parameter_summary/opossum_parameter_estimates.csv`), and raccoon (`./parameter_summary/raccoon_parameter_estimates.csv`). Each file contains the median estimate of each parameter, 95% credible interval, effective sample size, and Gelman-Rubin diagnostic. The file `./parameter_summary/metadata.csv` contains information to link the parameter names (as specified in the `JAGS` script) to the model explanation in the manuscript.

[Back to table of contents ⤒](#a-repository-for)

### The R folder (`./R`)

This folder has 12 R scripts, plus an additional sub-folder titled `./R/functions` that contains 2 seperate scripts of different functions for running the analysis or plotting the results.

| File                                 | Description                                                                                                                                                                                                                                                                                                                                                            |
|--------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **`./R/Appendix_S2.Rmd`**          | A markdown file to create Appendix S2 for this manuscript. It houses the model validation results. |
| **`./R/clean_conflicts.R`**          | This script goes through the raw conflict files in `./data/conflict_raw` and uses regular expressions to query nuisance wildlife reports where the species ID may be uncertain. It then loops through the descriptions of each offending report and prints out the description tied to the report (and allows the user to determine if the report should be censored). |
| **`./R/figure1.R`**                  | The script I wrote to create the first figure in the manuscript.                                                                                                                                                                                                                                                                                                       |
|**`./R/fit_models_validation.R`**                  | The script used to fit the models to a subset of the data for validation purposes. Similar to `fit_models_validation.R`.                                                                                                                                                                                                                                                                                                     |
| **`./R/format_data_for_analysis.R`** | The main workhorse of this project. All of the data is pulled in with this script and it gets formatted and made ready to fit into a JAGS model. The main outputs from this are `data_list`, which is the list object supplied to `JAGS`, as well as the `my_inits()` function, which supplies initial values to JAGS. It's commented out as well.                     |
| **`./R/format_data_for_validation.R`** | Similar to `./R/format_data_for_analysis.R` expect it subsets the data before fitting the model.                     |
|**`./R/geocode_conflicts.R`**        | Once `./R/clean_conflicts.R` has been ran, this is the script that will take the block-level address and geocode it to latitude and longitude.                                                                                                                                                                                                                         |
| **`./R/plots_models.R`**             | The script I wrote to make the rest of the figures for this manuscript. `./R/format_data_for_analysis.R` still gets ran as it provides a lot of the spatial data we need for making maps.                                                                                                                                                                              |
| **`./R/sourcer.R`**                  | A script that just sources all of the functions in `./R/functions`. That way I can just call `source(./R/sourcer.R)` instead of sourcing each file in that folder.                                                                                                                                                                                                     |
| **`./summarise_high_res_lulc.R`**    | The spatial environmental data we used had sub-meter resolution, which we needed to aggregate up to 500 meters for this analysis. This script does that aggregration and prepares the environmental data layers for statistical analysis.                                                                                                                              |
 |**`./URB_PCA_plot.R`**    | Used to create a rough draft of Appendix S1: Fig S1|
 |**`./validate_model.R`**    | Calculates the AUC from the posterior distributions created via `./R/fit_models_validation.R`|

#### The functions sub-folder (`./R/functions`)

This sub-folder houses a bunch of custom functions I wrote to streamline this project. I've commented out all of the functions and their arguments in the scripts themselves, so I leave those who are interested in what functions are there to explore the scripts themselves. The two scripts here are`./R/functions/plot_utility.R`, which houses a lot of the functions I used to make the figures for the project, and `./R/functions/utility_script.R`, which was a catchall for non-plot related functions.


[Back to table of contents ⤒](#a-repository-for)




