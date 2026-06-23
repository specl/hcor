# hCor: Hierarchical Correlation & Attenuation Correction

hCor is an R utility designed to estimate the latent correlation between two experimental tasks while accounting for measurement error and trial-by-trial variance. It provides both Bayesian hierarchical modeling (via JAGS) and classical psychometric corrections (Spearman's correction for attenuation).

## Reference

Rouder, JN (2026) Hierarchical-Model Correlations Should Replace Sample Correlations in Experimental Tasks.

## Help

jrouder@uci.edu 

## Why use hCor?

When calculating correlations between experimental tasks, trial-level noise attenuates (shrinks) the observed correlation coefficient. hCor solves this by modeling the underlying true scores hierarchically, giving you an unattenuated estimate of the true relationship.

It supports two designs:
- Model 1: Simple baseline/mean comparison between two tasks.
- Model 2: Within-subject experimental contrast comparisons (e.g., Condition 2 minus Condition 1) across two tasks.

## Prerequisites

Before running hCor, you must have JAGS (Just Another Gibbs Sampler) installed on your machine.
Download JAGS here: https://mcmc-jags.sourceforge.io/

Then, install the required R packages in your R console:
```R
install.packages(c('R2jags', 'coda'))
```

## Data Structure Requirements

Your input data frame must be in long format and contain the following column names exactly:

- sub: Integer. Subject ID. Must be contiguous integers starting from 1 (e.g., 1, 2, 3...).
- task: Integer. Task ID. Must be exactly 1 or 2.
- y: Numeric. The dependent variable / trial-level score.
- cond: Integer. (Optional) Condition ID for Model 2. Must be exactly 1 or 2.

## User-Facing Functions Reference

### hcorRun(dat, prior, chains, iter, burnin)
Runs the Bayesian hierarchical correlation model via JAGS. It automatically detects whether to run Model 1 or Model 2 based on whether a cond column exists in your data.

- Arguments:
  - dat: Your long-format data frame.
  - prior: A list containing hyper-parameters for the priors.
  - chains: Number of MCMC chains to run (default = 4).
  - iter: Total number of MCMC iterations per chain (default = 1000).
  - burnin: Number of initial iterations to discard as burn-in (default = 100).
- Returns: A jags object containing the full MCMC posterior distributions.

### hcorDiagnostic(chains)
Calculates convergence diagnostics for the key parameters using the Gelman-Rubin statistic.

- Arguments:
  - chians: The jags object returned by hcorRun().
- Returns: A gelman.diag object. Look for values close to 1.00 to confirm convergence.

### hcorConv(dat)
Computes conventional, non-hierarchical statistical metrics for comparison.

- Arguments:
  - dat: Your long-format data frame.
- Returns: A list containing:
  - cortest: The standard Pearson correlation (cor.test) on observed subject means/contrasts.
  - reliability: The estimated reliability of each task.
  - spearman: The classical Spearman correction for attenuation.

### hcorPlot(chains, dat)
Generates a density histogram of the posterior distribution for the correlation coefficient. It overlays the standard Pearson correlation (dot and line at the top), the classical Spearman correction (asterisk), and the Bayesian credible interval (blue line) for easy visual comparison.
- Arguments:
  - chains: The jags object returned by hcorRun().
  - dat: Your long-format data frame.
- Returns: Generates a base R plot.

### hcorRho(chains)
A quick helper utility to extract the raw vector of MCMC posterior samples for the correlation coefficient.
- Arguments:
  - chains: The jags object returned by hcorRun().
- Returns: A numeric vector of posterior samples.

### hcorBF(chains,interval)
Calculates a Savage-Dickey density ratio approximation of the Bayes Factor for a specified hypothesis interval versus a uniform prior.
- Arguments:
  - chains: The jags object returned by hcorRun().
  - interval: A numeric vector of length 2 defining the hypothesis window (e.g., c(0, 1) to test a positive correlation).
- Returns: A list containing the boundary interval, the calculated Bayes Factor (bf), sample counts, and a logical flag indicating if the count hit an extreme boundary.

## Quick Start Example

### 1. Load the functions
Source the file containing your hCor functions:
```R
source('hCor.R')
```

### 2. Define your Priors
JAGS requires a list of hyper-parameters for the priors. Here is a standard baseline template:
```R
myPrior <- list(
  tau.scale=c(175,175),  # trial noise / measurement error
  mu.m=c(60,60),         # slope means
  mu.sd=c(30,30),        # slope sd
  sig.lower=c(10,10),    # lower bound for true score SD
  sig.upper=c(100,100),  # upper bound for true score SD
  alpha.m=c(800,800),    # intercept means (Model 2 only)
  alpha.sd=c(1000,1000)) # intercept standard deviations ( Model 2 only)
```

### 3. Execution Pipeline
```R
fit <- hcorRun(dat = myData, prior = myPriors, chains = 4, iter = 2000, burnin = 500)
hcorDiagnostic(fit)
conventional_results <- hcorConv(myData)
bf_positive <- hcorBF(fit, interval = c(0, 1))
hcorPlot(fit, myData)
```
