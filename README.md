[![DOI](https://zenodo.org/badge/419362229.svg)](https://zenodo.org/doi/10.5281/zenodo.10413280)

# A repository for:

Fidino, M et al. 2024.  Gentrification drives patterns of alpha and beta diversity in cities. PNAS.

[Link to study here](https://www.pnas.org/doi/10.1073/pnas.2318596121)


## Links to different parts of the readme file

1.  [What's in this repository?](#whats-in-this-repository)
2.  [The working directory](#the-working-directory)
3.  [city codes](#city-codes)
4.  [The data folder (`./data`)](#the-data-folder-data)
5.  [The JAGS folder(`./JAGS`)](#the-jags-folder-jags)
6.  [The mcmc output folder (`./mcmc_output`)](#the-mcmc-output-folder-mcmc_output)
7.  [The mcmc plots folder(`./mcmc_plots`)](#the-mcmc-plots-folder-mcmc_plots)
8.  [The plots folder(`./plots`)](#the-plots-folder-plots)
9.  [The R folder (`./R`)](#the-r-folder-r)
10. [The results folder (`./results`)](#the-results-folder-results)
11. [The supplemental folder (`./supplemental`)](#the-supplemental-folder-supplemental)


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


### The mcmc output folder (`./mcmc_output`)

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


### The mcmc plots folder (`./mcmc_plots`)

This folder is intentionally left blank, and stored temporary files for the traceplots of model parameters. 

[Back to table of contents ⤒](#a-repository-for)



### The plots folder (`./plots`)

This folder houses some of the raw figures I generated in `R` (which I cleaned up using Inkscape), as well as other figures that were publication ready. All of the extra "cleaning" I needed to do was related to the maps I had made (there was too much spacing between images).


|File|Description|
|---|---|
|**`./plots/alpha_beta_expected.tiff`** | Figure 5 in the manuscript.|
|**`./plots/alpha_beta_results.tiff`** | Figure 3 in the manuscript.|
|**`./plots/alpha_examples.tiff`** | Variation in alpha diversity as a function of gentrification and impervious cover. Used in presentations as an example.|
|**`./plots/beta_examples.tiff`** | Variation in beta diversity as a function of gentrification and impervious cover. Used in presentations as an example.|
|**`./plots/figure_1.tiff`** | Figure 1 in the manuscript.|
|**`./plots/gentrification_map.tiff`** | Figure 2 in the manuscript.|
|**`./plots/gentrification_map_supp_finished.svg`** | Figure 2 in the manuscript, but with names added to it (for supplemental material). Stored as an svg file so I could edit it in Inkscape.|
|**`./plots/occupancy_result.tiff`** | Figure 4 in the manuscript.|

[Back to table of contents ⤒](#a-repository-for)


### The R folder (`./R`)

This folder has 25 R scripts, plus an additional sub-folder titled `./R/functions` that contains 7 seperate scripts of different functions for running the analysis or plotting the results.

Because there are so many scripts in here, I have ordered them here in such a way that if someone was interested in completely running the analysis, they could do so.


#### Step 1. Compile site coordinates

| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
|**./R/compile_coordinates.R**| Read in detection data and get a unique set of sites for each city. Save output. | `dplyr`|

#### Step 2. Query Census data

| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
| **./R/pull_census.R** | Query the Census data in each city across a few different years. | `dplyr`, `tidycensus`, `sf` |
|**./R/summarise_census.R**| Summarise the Census data after it has been queried. | `sf`, `dplyr`, `uwinspatialtools`|
|**./R/census_functions.R**| Helpful functions to be used within `summarise_census.R`. Each function is also explained with `summarise_census.R` | `sf`, `dplyr`, `uwinspatialtools`|

NOTE: `uwinspatialtools` is an R package I deleveloped, which essentially has some wrapper functions for `sf` and `raster`. It can be found at www.github/com/mfidino/uwinspatialtools.

#### Step 3. Pull impervious cover data

| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
| **./R/pull_imperv.R** | Query impervious cover in each city across a few different years. | `dplyr`, `raster`, `sf`, `FedData` |

NOTE: When I wrote these scripts `raster` was still usable. I'm not rewriting them now that it will soon be retired. 

#### Step 4. Pull species range data

| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
| **./R/pull_ranges.R** | Pulls IUCN range data and calculates the distance to the edge of a species range froma given city. The great lakes get filled in as well, or else all Chicago distances would just be to Lake Michigan. | `dplyr`, `sf`, `mapview`, `smoothr`|
|**.R/functions/common_to_binomial.R** | Converts common name for species to scientific name so that it can be queried.| Only base R used | 

#### Step 5. Quantify gentrification

Since I wanted to have the outputs from this in supplemental material, all data and code to classify locations as gentrified can be found in `./supplemental/supplemental.Rmd` in the `Additional Gentrification Metrics` section.

#### Step 6. Fit the multi-city, multi-species occupancy model


| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
|**./R/prep_data_occupancy.R**| Get's ran via `./R/fit_occupancy_model.R`, but it pulls in all the necessary data and gets it ready for analysis. | `dplyr`, `runjags`|
|**./R/fit_occupancy_model.R**| Fits the occupancy model and saves the output| `dplyr`, `runjags`|

#### Step 7. Estimate alpha and beta diversity with uncertainty

| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
| **./R/simulate_latent.R** | Given an input (`alpha` or `beta`), use the output from the occupancy model to estimate species presence on the landscape and then calculate alpha or beta diversity | `dplyr`, `runjags`|
| **./R/sample_z.R** | Used in the other script, does all the sampling for the species incidence matrix given the posterior. | `dplyr`, `runjags`|


#### Step 8. Fit alpha diversity model

| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
|**./R/fit_alpha.R** | Fits alpha diversity model and saves output. | `dplyr`, `runjags`|

#### Step 9. Fit beta diversity model

| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
|**./R/fit_beta.R** | Fits beta diversity model and saves output. | `dplyr`, `runjags`, `gtools`|

#### Step 10. Summarise alpha and beta diversity models and plot results
| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
| **./R/alpha_beta_functions.R**| Utility functions for summarizing posteriors and plotting | Only base R used|
| **./R/plot_alpha.R** | Make some figures related to alpha diversity analysis. None are used in the manuscript, but I retained the code here as I've used them in presentations. |`runjags`, `bbplot` | **./R/gdm_functions.R** | Utility functions for the beta diversity analysis | `vegan`, `bbplot`, `scales` |
| **./R/plot_beta.R** |Similar to `./R/plot_alpha.R`, these are some code for some extra figures used in presentations and the like. | `bbplot`, `sf`, `runjags`, `dplyr` |
| **./R/alpha_beta_plots.R** | The actual script to generate most of the analysis related figures in the manuscript. | `sf`, `runjags`, `dplyr`, `cli`, `MCMCvis`, `pals`, `bbplot`, `scales` |
| **./R/make_supp_maps.R** | Make city maps for supplemental | `pander`, `dplyr`, `bbplot`, `viridis`, `raster`, `sf`, `runjags`, `prettymapr`|


NOTE: `bbplot` is a plotting package developed by Juniper Simonis that I have collaborated on. It can be found at https://github.com/dapperstats/bbplot. 


#### Step 11. Summarise occupancy model and plot results


| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
| **./R/occ_model_unpack.R**| Summarise the occupancy model. | `dplyr`, `data.table`|
| **./R/plot_occupancy.R** | Make the occupancy model figure in the manuscript | `runjags`, `dplyr`, `bbplot` |


#### Step 12. Some additional descriptive data calculations

| File                  | Description                                                      | Packages required           |
|-----------------------|------------------------------------------------------------------|-----------------------------|
| **./R/income_across_gentrification.R**| Calculate differences in income across gentrified and not gentrified sites. | `sf`, `dplyr`, `uwinspatialtools`|
| **./R/greenspace_by_census_tract.R**| Calculate differences in greenspace across gentrified and not gentrified Census tracts. | `sf`, `dplyr`, `uwinspatialtools`|
| **./R/imperv_by_census_tract.R**| Calculate differences in income across gentrified and not gentrified Census tracts. | `sf`, `dplyr`, `uwinspatialtools`|


#### The functions sub-folder (`./R/functions`)

This sub-folder houses some of custom functions I wrote to streamline this project. 


[Back to table of contents ⤒](#a-repository-for)


### The results folder (`./results`)

This folder stores the MCMC output from many of the models. While `./mcmc_output` stores the results from the alpha and beta diversity models, the alpha and beta estimates (with associated uncertainty) are stored here. Some other RDS files can also be stored here (e.g., occupancy model output), but have been added to the `.gitignore` because they are way too big.


|**`./results/alpha_for_stage_two_collapsed.csv`** | The alpha diversity estimates as well as the associated uncertainty. This csv has 999 rows and 6 columns.

| Column      | Date type    | Explanation                                                                                |
|-------------|--------------|--------------------------------------------------------------------------------------------|
| Site        | `factor`     | The site associated to this data point                                                     |
| City        | `factor`     | The city abbreviation associated to this data point                                        |
| gentrifying | `boolean`    | Whether the site was gentrified (TRUE) or not (FALSE)                                      |
| mean_19     | `proportion` | The proportion of impervious cover within 1000 km of a site in 2019. Ranges from 0 to 100. |
| mu          | `numeric`    | The estimated species richness at a site across all survey periods                         |
| sd          | `numeric`    | The estimated variation in species richness estimates at a site across all survey periods  |

|**`./results/beta_summary_for_analysis_collapsed_vegan.csv`** | The beta diversity estimates as well as the associated uncertainty. This csv has 27,757 rows and 6 columns.


| Column   | Date type    | Explanation                                                                                                                                                                          |
|----------|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| loc      | `factor`     | The row location, used to line up with the covariate data                                                                                                                            |
| mu_beta  | `proportion` | The estimated pairwise dissimilarity metric between a pair of sites. Ranges between 0 and 1, where 0 is identical communities and 1 is completely different communities.             |
| var_beta | `numeric`    | The estimate variation in the pairwise dissimilarity metric.                                                                                                                         |
| mu_rich  | `numeric`    | Estimated species richness across the pair of sites.                                                                                                                                 |
| var_rich | `numeric`    | Estimated variation in species richness across a pair of sites.                                                                                                                      |
| na_count | `numeric`    | The number of MCMC samples for a specific pair that dissimilarity could not be calculated. Note that we were always able to make at least 2500 estimates for each pair of locations. |

|**`./results/city_gamma_samples.csv`** | MCMC samples of gamma diversity estimates for each city (i.e., the total number of the 21 species estimated to be present at the city scale). This csv has 10,000 rows (MCMC samples) and 23 columns (one for each city).


| Column         | Date type | Explanation                                                                                                                           |
|----------------|-----------|---------------------------------------------------------------------------------------------------------------------------------------|
| V1 through V23 | `integer` | MCMC samples of gamma diversity across each city. Cities are ordered alphabetically based on their 4 letter abbreviation (see README) |


[Back to table of contents ⤒](#a-repository-for)


### The Supplemental folder (`./supplemental`)

This folder contains most of the necessary data to construct the entire supplemental material for this manuscript within a R markdown file.

|**`./supplemental/supp_maps`** | sub-folder that contains 23 maps, one for each city, that are used in the supplemental material.


|**`./supplemental/occupancy_summary.RDS`** | The median estimate of every model parameter in the multi-city, multi-species occupancy model as well as uncertainty (90% and 95% CI).

|**`./supplemental/science_logo.png`** | Used to place at the top of the supplemental material

|**`./supplemental/supplemental.Rmd`** | The R markdown file of the supplemental material.

|**`./supplemental/Fidino_2023_supplemental.pdf`** | A rendered PDF file of the supplemental material.


[Back to table of contents ⤒](#a-repository-for)
