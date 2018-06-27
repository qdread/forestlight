# Forestlight code pipeline

Last updated by QDR on 25 June 2018

## File locations:

- Code is on the repo `github.com/qdread/forestlight`, referred to as *GitHub* in file paths below.
- Data is on google drive at `ForestLight/data`.
- Any intermediate figures are on google drive at `ForestLight/figs`.
- Function fitting takes a lot of computing time/power so is on MSU's cluster. The huge files with the raw MCMC samples are kept there and only summary files are put on google drive.

## Pipeline

The pipeline is organized by script. I tried to summarize what each script does, and say what input files are required and what output files each one makes.

### Step 1: Quality control

`GitHub/forestlight/code/complete_workflows/bciqc.r`

*Inputs required*
 
- Condit's raw BCI data (all in `google_drive/ForestLight/data/BCI_raw/bcidata`)
- Wood specific gravity and taper parameters for BCI trees (both in `google_drive/ForestLight/data/BCI_raw/bci_taper`), provided by KC Cushman
- Location of quadrats that are "young" or secondary forest, provided by Nadja (`google_drive/ForestLight/data/BCI_light/habitats2.txt`)

*What the script does*

- Remove tree ferns and strangler figs
- Correct dbh for tapering and buttressed species
- Recalculate biomass allometry using corrected dbh and using the wood specific gravity values for each species
- Assign trees to either being in the 42.84-ha study area, the "young" area, or within 20 m of the edge
 
*Outputs generated*

- R data object `google_drive/ForestLight/data/BCI_raw/bcidata/bciqcrun.R`

### Step 2: Lots of data processing

`GitHub/forestlight/code/complete_workflows/workflow_june_newFGs.r`

*Inputs required*

- R data object `bciqcrun.R` from previous step
- Nadja's functional group classifications (`google_drive/ForestLight/data/Ruger/fgroups_dynamics_new.txt`)

*What the script does*

- Exclude the young-forest trees and edge trees from the dataset
- Calculate annual biomass increment from the 5-year biomass increments
- Remove outliers from the biomass increment dataset that appear to be errors
- Convert all the values to the right units, and get per-hectare values by dividing by the study area where needed
- Split the data into separate datasets for each functional group
- Calculate all allometries (light received, height, crown area)
- Do all logarithmic binning, both for all trees together each census and for each functional group separately each census (bin abundance, total growth, and individual/total light)
- Save raw and binned data
- Make some figures of the binned data for temporary visualization

*Outputs generated*

- Raw data R object
- Binned data R object
- Individual CSVs for each set of binned data

### Step 3: Function fitting

The basic procedure here is to use Stan to fit the models. That requires a few things. First, you specify the model in a `.stan` file, which Stan translates into C++ code, then compiles into an executable file. You also need the input data in text format, which can be made in R using a so-called "Rdump." Then, since the models take a lot of time and power to run, we put everything on MSU's cluster. One shell script contains the instructions for the different calls to the Stan programs, and a text file of qsubs runs that script multiple times to do calls for the different functional form/guild/year combinations. Finally, a number of other scripts are used to get information out of the large CSVs that Stan generates with the sampler's output. 

The four `.stan` scripts with the model code for the four different models are all in `GitHub/forestlight/stan`:

- `model_ppow_withlik.stan`
- `model_pexp_withlik.stan`
- `model_wpow_withlik.stan`
- `model_wexp_withlik.stan`

The `.stan` script with the model code for the von Bertalanffy fit is `GitHub/forestlight/stan/vonb.stan`.

All the shell scripts and R scripts for working with the Stan scripts are in `GitHub/forestlight/stan/final_workflow/`.

#### Step 3a: Create data objects usable by Stan

The Rdumps are done with the scripts `datadump_densprod.r` and `datadump_light.r` for the main models and the light von Bertalanffy models, respectively.

*Inputs required*

The raw data R objects from step 2

*Outputs generated*

R objects called Rdumps that Stan can use to load the data


#### Step 3b: Fit models with Stan

The script `fitproduction.sh` calls the Stan programs to fit the density-production models, if you input which model, which guild, and which year, along with the number of sampling and warm-up iterations. A list of `qsub` calls that will run this script remotely are in the document `cmdstanqsubs.sh`. Cmdstan is the name of Stan's command line interface. The script `fitlight.sh` does the same for the light models.

*Inputs required*

- Rdumps from step 3a
- The compiled Stan models

*Outputs generated*

Each model generates a huge CSV with all the MCMC samples for the parameters.

#### Step 3c: Extract summary information from raw Stan output

*density-production models*: The script `extraction_functions_productionfits.r` has the source code for 5 functions. The first four functions extract the following from a model fit:

- Parameter values and credible intervals
- Fitted values and credible intervals (and prediction intervals)
- Fitted log slopes and credible intervals
- Bayesian R-squared and credible intervals

The final function calls all these functions as well as running LOOIC on the model fit.

The functions in the above script are run on each model in parallel in the script `extract_ci_productionfits.r`. Since that script is run in parallel, you next need to combine all the output into CSV files using `combinestanoutput.r`.

*light models*: The script `extract_ci_lightfits.r` gets all needed information out of the light model fits. It calculates median and quantile values for the model parameters and for fitted values for the curves, as well as for the maximum log slope. Lastly a separate part of the script manually calculates the Bayesian R-squared, using the correct method endorsed by Gelman. This function isn't needed to be run in parallel because it's a lot quicker and runs in a few seconds for every model separately.

*Inputs required (for both density-production models and light models)*

The CSV files with the MCMC samples

*Outputs generated*

First, R objects for each model, then the combining script loads those R objects and writes them into CSVs with summary information from all the models

### Step 4: Plot and analyze results of function fitting

#### Create plotting data for density-production models

`GitHub/forestlight/code/data_export/createplotdata_final.r`

*Inputs required*

- The raw data R object created in step 2 
- The fitted model values created in step 3c (for density-production models)

*Outputs generated*

- Several CSVs with observed and predicted values in separate files

#### Create plotting data for light models

`GitHub/forestlight/code/data_export/createlightplotdata_final.r`

*Inputs required*

- The raw data R object created in step 2 (for light & growth per unit area)
- The fitted model values created in step 3c (for light models)

*Outputs generated*

- Several CSVs with observed and predicted values in separate files

#### Draw density-production plots

`GitHub/forestlight/code/plotting/plottingcodemodelfits.r`

Uses the observed and predicted CSVs to draw many different plots based on the density and production models.

#### Draw light plots

`GitHub/forestlight/code/plotting/plottingcodelightmodelfits.r`

Uses the observed and predicted CSVs to draw many different plots based on the light models.