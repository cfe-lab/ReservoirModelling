# Instructions for performing our analysis

First, it should be noted that we've typically run our analyses in an
*interactive* R console using the `source()` command rather than as 
self-contained scripts called from a command line.  In most cases we
imagine that the scripts *would* work in this fashion with minimal 
tweaking, but we preferred to do it this way so that we could perform
the slow dynamical system simulations only once and then "build upon 
the calculations" by simply `source()`-ing the followup script, which
would typically perform the plotting or further analyses.

## Input data

Our scripts currently look up the input files from hard-coded locations.
If your code is in 

```
[project root]/ReservoirModelling
```

then place the CSV files `BrooksData.csv` and `BrooksDataVLs.csv` in

```
[project root]/data
```

## `example_sampled_distributions.r`

This produces "composition plots" with only the empirical distributions
for two specific cases:

 - N133M
 - Z634F

## "Typical behaviour" analyses and plots

These analyses are independent of any participant data and purely rely
on our dynamical system with "typical" conditions.

As the dynamical system is relatively slow to run, we perform that calculation
in its own script: `typical_pvl_demonstration_ode.r`.  Running this script will
produce an `.RData` file that contains all of the variables defined in the
script; this is so that you can later load this file to do other analyses
rather than having to rerun the dynamical system simulation.

To perform these analyses, first either run this script in an interactive
R console using `source()` to set up all the variables needed for subsequent
analyses, or if you've already done this then you can simply load the `.RData`
file it produced.  The various follow-up analyses that build on this script are:

### `typical_pvl_demonstration_plots.r`

This produces the plots of "typical pVL curves" where ART is initiated after either
3 years or 7 years of infection, along with corresponding plots of the reservoir
composition in each of these settings.

### `particpant_vl_vs_typical.r`

This produces plots of each individual's pVL curve with the "typical" curve overlaid.

### `dynamical_system_all_compartments.r`

This produces plots of all of the compartments in the dynamical system.

## Dynamical system-based reservoir analyses

These analyses all use the participants' measured VL data (which have some "grafted"
data for the acute phase) to inform the dynamical system which models their infection
and creation of their latent reservoir.

First, either run `brooks_data_ode.r`, which will perform all the simulations and generate
`brooks_data_ode.RData`, which will contain the results; or, if you've already done that,
simply load the `.RData` file if needed.

The analysis scripts below can then be run.

### `ode_based_compute_estimates.r`

This computes the likelihoods and the MLEs if you want to directly look at them without
generating any further plots.  Most of this analysis is repeated in the scripts that actually
generate plots, but this one also performs a sanity check on the data to make sure
that the infection dates and ART initiation dates are cogent between the viral
load data and the integration data.

### `ode_based_likelihood_plots.r`

This produces plots of the likelihoods used to estimate the decay rate, with the
MLEs and their error bounds marked where possible.

### `ode_based_composition_plots.r`

This produces the plots of reservoir composition, with the empirical distribution
plotted as a bar graph and with modelled distributions shown with stepped lines.
Only one of them has a legend, as they are all meant to be presented on a single page.

### `ode_based_half_life_plots.r`

This generates plots of the participants' half-lives, using estimates drawn from the
above analyses.

### Analyses performed using integration data with duplicate proviruses retained

Each of the above scripts also has an analogue that performs the same analysis
on data where duplicate proviruses are retained:

 - `ode_based_with_duplicates_compute_estimates.r`
 - `ode_based_with_duplicates_likelihood_plots.r`
 - `ode_based_with_duplicates_composition_plots.r`
 - `ode_based_with_duplicates_half_life_plots.r`

## Poisson GLM-based decay rate estimates

These analyses do not use the dynamical system at all, and instead infer a decay rate
(and therefore a half-life) using a Poisson GLM fitted directly to the estimated 
integration dates (under various different binning schemes).

### `model_free_half_life_plots.r`

In this script, integration dates are binned together by year prior to ART initiation.
This means that the "earliest" year is typically abbreviated, as the person was not 
infected for the entirety of that year.

### `model_free_no_binning_half_life_plots.r`

This is like the previous script, but with "no binning" on integration dates: that is,
the data is binned by day, which is the limit of resolution of our integration date
data.

### Analyses performed using integration data with duplicate proviruses retained

Each of the above scripts also has an analogue that performs the same analysis
on data where duplicate proviruses are retained:

 - `model_free_with_duplicates_half_life_plots.r`
 - `model_free_with_duplicates_no_binning_half_life_plots.r`
