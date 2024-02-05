Bayesian semiparametric modeling of spatially referenced multistate
current status data
================

Implementation of modeling framework proposed in the paper, Das, S.,
Chae, M., Pati, D., and Bandyopadhyay, D. (2024+) “Bayesian
semiparametric modeling of spatially referenced multistate current
status data”.

## Overview

Assessment of multistate disease progression is commonplace in
biomedical research, such as, in periodontal disease (PD). However, the
presence of current status endpoints, where subjects’ teeth
transitioning through various PD states are observed at a single random
inspection time, complicates the inferential framework. In addition,
these endpoints are clustered (teeth within subjects), and spatially
associated i.e., proximal teeth share similar PD status compared to
distal ones.

We propose a Bayesian semiparametric accelerated failure time model with
(spatial) random effects, and flexible errors that follow a Dirichlet
process mixture of Gaussians. For elegant interpretation, we formulate a
monotone single index model, whose (unknown) link function is estimated
via a novel integrated basis expansion with constrained Gaussian process
on the basis coefficients. Efficient posterior computing is achieved
using elliptical slice sampling and fast circulant embedding techniques.

## Main Functions

- **data_generate.R** - Functions to generate clustered multistate
  currentstatus data. The functions `data_model` and
  `data_misspecification` generate data from the correctly specified and
  misspecified models respectively.

- **SPMSM.cpp** - Functions to implement our proposed Bayesian modeling
  framework. The function `GP_MSM` performs posterior sampling of
  parameters of our proposed model.

- **SOP_TP.R** - Functions to calculate the state occupation probability
  (SOP) and transition probability (TP) of the disease states over time.
  The functions `get_sop` and `get_tp` generate the SOP and TP,
  respectively, using posterior samples of parameters obtained from
  `GP_MSM`.

## Additional Functions

- **singleIndex.h** - Functions to fit a Bayesian monotone single index
  model using constrained Gaussian process and Bernstein polynomial
  basis functions.

- **circulantembedding.R**, **inv_chol.cpp** - Functions to estimate a
  monotone function using a constrained Gaussian process basis
  expansion, by employing the algorithm proposed by [Ray et
  al. (2020)](https://link.springer.com/article/10.1007/s11222-020-09922-0).
  Code is taken from the repository,
  [raypallavi/BNP-Computations](https://github.com/raypallavi/BNP-Computations).

- **DPM.h** - Functions to fit a Dirichlet process mixture of Gaussians,
  using the blocked Gibbs sampler proposed by [Ishwaran & James
  (2001)](https://www.tandfonline.com/doi/abs/10.1198/016214501750332758).

- **rTruncDist.h** - Functions to generate samples from truncated
  distributions.

- **LogLik.h** - Functions to evaluate the log-likelihood of the fitted
  models.

- **RcppBasic.h** - Functions to basic vector and matrix operations in
  Rcpp.
