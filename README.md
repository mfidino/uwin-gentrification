

# A repository for:

Fidino, M et al.  Gentrification drives patterns of alpha and beta diversity in cities

## Links to different parts of the readme file

1. [What's in this repository?](#whats-in-this-repository)
2. [The working directory](#the-working-directory)
3. [city codes](#city-codes)
4. [The data folder (`./data`)](#the-data-folder-data)
5. [The JAGS folder(`./JAGS`)](#the-jags-folder-jags)
6. [The mcmc output folder (`./mcmc_output`)](#the-mcmc-output-folder-mcmc_output)
7. [The mcmc plots folder(`./mcmc_plots`)](#the-mcmc-plots-folder-mcmc-plots)
5. [The figures folder (`./figures`)](#the-figures-folder-figures)
6. 

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

This folder has 3 files and 2 sub-folders.

|**`./data/detection_data.csv`** | The camera trap data used in our analysis.
This csv file has 576,000 rows and 9 columns. 

| Column  | Date type | Explanation                                                                                                                                                                                                   |
|---------|-----------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Species | `factor`  | The species associated to this data point                                                                                                                                                                     |
| Season  | `factor`  | A seasonal code for the season the data comes from. It combines the first two letters of the season and the last two digits of the year. Seasonal codes are SP = Spring, SU = Summer, FA = Fall, WI = Winter. |
| Site    | `factor`  | The site associated to this data point                                                                                                                                                                        |
| City    | `factor`  | The city abbreviation associated to this data point                                                                                                                                                           |
| Long    | `numeric` | The longitude of the site associated to this data point (WGS 84)                                                                                                                                              |
| Lat     | `numeric` | The latitude of the site associated to this data point (WGS 84)                                                                                                                                               |
| Crs     | `integer` | The coordinate reference system code for the site coordinates                                                                                                                                                 |
| Y       | `integer` | The number of days the species was detected during a given season                                                                                                                                             |
| J       | `integer` | The number of days the camera trap was operational during a given season. If `J == 0` then no sampling occurred.                                                                                              |


|**`./data/gentri_all_coordinates.csv`** | The coordinates for all sites in our analysis. This dataset has 1000 rows and 4 columns.

| Column | Date type | Explanation                                                      |
|--------|-----------|------------------------------------------------------------------|
| Site   | `factor`  | The site associated to this data point                           |
| City   | `factor`  | The city abbreviation associated to this data point              |
| Long   | `numeric` | The longitude of the site associated to this data point (WGS 84) |
| Lat    | `numeric` | The latitude of the site associated to this data point (WGS 84)  |


|**`./data/species_in_analysis.csv`** | The common names of the species analyzed in our occupancy model. This dataset has 21 rows and 1 column.

| Column  | Date type | Explanation      |
|---------|-----------|------------------|
| species | `factor`  | The species code |


|**`./data/census_data/`** | This sub-folder contains all of the ancillary files generated to classify gentrification across Census tracts in each city. The files themselves are not, overall, necessary to run the entire analysis. This sub-folder contains one folder for each of the three metrics used to classify gentrification (educational attainment, median income, and race), as well as many `RDS` files that are summaries of the Census data (e.g., `./data/census_data/gent_tracts.rds` is a spatial data.frame of which sites are and are not gentrified). Explanations of each object here are not necessary, but to see how they were all created you can look at `./supplemental/supplemental.Rmd`.

|**`./data/cleaned_data/covariates`** | This sub-folder contains the covariates used in the occupancy analysis, as well as one file used to make supplemental maps.

|**`./data/cleaned_data/covariates/census_tract_imperv_cover.csv`** | Used to make supplemental maps. This csv has 7816 rows and 5 columns.

| Column  | Date type    | Explanation                                                                         |
|---------|--------------|-------------------------------------------------------------------------------------|
| city    | `factor`     | The four letter code for a city                                                     |
| gent    | `boolean`    | Whether the Census tract was gentrified (1) or not (0)                              |
| vuln    | `boolean`    | Whether the Census tract was vulnerable to gentfrication (1) or not (0)             |
| mean_11 | `proportion` | The proportion of impervious cover in a Census tract in 2011. Ranges from 0 to 100. |
| mean_19 | `proportion` | The proportion of impervious cover in a Census tract in 2019. Ranges from 0 to 100. |


|**`./data/cleaned_data/covariates/dist_city_to_species_edge.csv`** | The distance of a city from the edge of a speecies known range, where negative values mean a city is within a species range and positive values mean a city is outside a species range. This dataset contains 21 rows and 24 columns.

| Column     | Date type | Explanation                                                                                                                                                              |
|------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| species    | `factor`  | The species used included in the occupancy model                                                                                                                         |
| City codes | `numeric` | The distance of a respective city, which is represented by the column name (columns 2 through 24 contain city name codes) from the edge of a species known distribution. Distances are in kilometers.|

|**`./data/cleaned_data/covariates/site_covariates.csv`** | The site specific covariates used in our analyses. This dataset contains 999 rows and 5 columns. Note that this dataset has 1 fewer row than `./data/gent_all_coordinates.csv`. This is because one site in Chicago had no one living nearby and so was dropped from the analysis.


| Column      | Date type    | Explanation                                                                         |
|-------------|--------------|-------------------------------------------------------------------------------------|
| City        | `factor`     | The four letter code for a city                                                     |
| Site        | `numeric`    | The site associated to this data point                                              |
| gentrifying | `boolean`    | Whether the site was gentrified (1) or not (0)                                      |
| vulnerable  | `boolean`    | Whether the site was vulnerable to gentrification (1) or not (0)                    |
| mean_19     | `proportion` | The proportion of impervious cover in a Census tract in 2019. Ranges from 0 to 100. |

[Back to table of contents ⤒](#a-repository-for)


### The JAGS folder (`./JAGS`)

This folder contains the `JAGS` scripts used to fit each model.

|**`./JAGS/beta_model_collapsed_norm.R`** | The beta diversity generlized dissimilarity model. Uses half-normal priors for positive slope terms.

|**`./JAGS/impute_alpha.R`** | The alpha diversity model.

|**`./JAGS/multi_scale_occupancy.R`** | The multi-species multi-city occupancy model.


[Back to table of contents ⤒](#a-repository-for)


### The mcmc_output folder (`./mcmc_output`)

This contains two folders, one for the alpha diversity analysis results (`./mcmc_output/alpha_output`) and another for the beta diversity analysis results (`./mcmc_output/beta_output`). However, the `RDS` files that contain the posteriors are too large and so are not stored on GitHub. 

|**`./mcmc_output/alpha_output`** | This sub-folder is blank, but stores the `RDS` file that contains the posterior from the alpha diversity analysis (which is added to the `.gitignore` file because it is large).


|**`./mcmc_output/beta_output`** | This sub-folder contains some of the objects needed to fit the generalized dissimilarity model. 

|**`./mcmc_output/beta_output/dmat_site_ids_collapsed.csv`** | This dataset links the sites together that are used to generate the pairwise dissimilarity metric. It contains  27,757 rows and 5 columns.

| Column   | Date type | Explanation                                                                                   |
|----------|-----------|-----------------------------------------------------------------------------------------------|
| siteA    | `factor`  | The name of one of the two sites                                                              |
| siteA_id | `numeric` | A numeric identifier for site A                                                               |
| siteB    | `factor`  | Whether the site was gentrified (1) or not (0)                                                |
| siteB_id | `numeric` | A numeric identifier for site B                                                               |
| City_id  | `numeric` | A numeric identifier for each city, the levels are just the city names sorted alphabetically  | 

|**`./mcmc_output/beta_output/knots.csv`** | This dataset has the knots to used to generate the I-splines for the generalized dissimilarity model. It contains 46 rows and 6 columns.

| Column    | Date type | Explanation                                                                                       |
|-----------|-----------|---------------------------------------------------------------------------------------------------|
| covariate | `factor`  | The covariate the knot is associated to. Either geographic distance or mean_19 (impervious cover) |
| nspline   | `numeric` | The number of splines used                                                                        |
| min       | `numeric` | The minimum knot                                                                                  |
| median    | `numeric` | The median knot                                                                                   |
| max       | `numeric` | The maximum knot                                                                                  |
| City      | `factor`  | The city code                                                                                     |


|**`./mcmc_output/beta_output/site_splines.csv`** | The I-splines used to fit the generalized dissimilarity model. This dataset contains 27,757 rows and 8 columns.

| Column | Date type | Explanation                           |
|--------|-----------|---------------------------------------|
| siteA  | `factor`  | The first site                        |
| siteB  | `factor`  | The second site                       |
| X      | `numeric` | The first geographic distance spline  |
| X.1    | `numeric` | The second geographic distance spline |
| X.2    | `numeric` | The third geographic distance spline  |
| X.3    | `numeric` | The first impervious cover spline     |
| X.4    | `numeric` | The second impervious cover spline    |
| X.5    | `numeric` | The third impervious cover spline     |



[Back to table of contents ⤒](#a-repository-for)


### The mcmc plots folder(`./mcmc_plots`)

This folder is intentionally left blank, and stored temporary files for the traceplots of model parameters. 

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




