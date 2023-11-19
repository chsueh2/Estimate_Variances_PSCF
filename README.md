# Estimate Variances of Model Parameters Using Perturbed SSE Curve Fitting (PSCF) Method

In this work, we implement Perturbed SSE Curve Fitting (PSCF) algorithm and compare it with the other common statistics methods including nonlinear regression models, bootstrap confidence intervals and the delta normality method.

[Project report](https://rpubs.com/clh2021/1114624)

Key features:

- Variance Estimates
- Nonlinear Regression
- Bootstrap confidence Intervals
- Delta Normality Method

R packages used:

- `conflicted`: tools to deal with R package conflicts
- `here`: enables easy file referencing and builds file paths in a OS-independent way
- `scales`: formats and labels scales nicely for better visualization
- `skimr`: provides summary statistics about variables in data frames, tibbles, data tables and vectors
- `glue`: embeds and evaluates R expressions into strings to be printed as messages
- `gt`: pretty-looking tables
- `GGally`: extends 'ggplot2' by combining geometric objects with transformed data
- `ggpubr`: provides easy-to-use functions for creating and customizing ggplot2-based publication ready plots
- `broom`: summarizes key information about statistical objects in tidy tibbles
- `nlshelper`: summarizes, tests, and plots non-linear regression models fits 
- `invgamma`: standard distribution functions for the inverse gamma distribution
- `tidyverse`: includes collections of useful packages like `dplyr` (data manipulation), `tidyr` (tidying data),  `ggplots` (creating graphs), etc.
- `zealot`: provides a `%<-%` operator to perform multiple, unpacking, and destructing assignment 

## Project Report

[Project report](https://rpubs.com/clh2021/1114624) ([Github Markdown](./z_main.md))([PDF](./z_main.pdf))([Quarto Markdown](./z_main.qmd)).

[Presentation slides](./Project_Presentation.pdf)

The analysis results with all theoretical backgrounds and math derivations are included.

Chien-Lan Hsueh (chienlan.hsueh at gmail.com)

## Overview and Project Goal

Perturbed SSE Curve Fitting (PSCF) Method proposed by Prof. Robert Hayes (Department of Nuclear Engineering, NCSU) is used to estimate variances of model parameters. 

- Implement Perturbed SSE Curve Fitting (PSCF) algorithm in R
- Is this method good in estimating variances of a model parameter?
- How does it compared to other statistical methods?


## Part 1 - Introduction 

Background and research questions:

1. how good the algorithm is in term of determining the variances of physical model parameters?
1. if there is a solid statistics ground to support and backup the validity of this method?
1. if yes in (2), can it be improved and further generalized to any physical models?
1. if no in (2), how well does it estimate as an approximation approach?

## Part 2 - Perturbed SSE Curve Fitting (PSCF) Method

Implementation of the PSCF algorithm.

## Part 3 - Study case: Attenuation decay data

A physical model "attenuation decay" is used to generate simulated measurement data with a predetermined randomness in the model parameter. We use PSCF method to estimate the variances of the model parameter and check how good PSCF method is.

## Part 4 - Comparison to Other Statistical Methods

In addition to the proposed algorithm Perturbed SSE Curve Fitting (PSCF) Method, we are going to use the following statistical methods to estimate the variance of an unknown parameter of interest.

1. Best Fit Model using Nonlinear Regression Model
1. Bootstrapping CI
   1. Parametric bootstrap
   1. Non-parametric bootstrap
1. Delta Method Normality

## Part 5 - Conclusion

After performing a extensive study on a simulated data set to test how well the algorithm performs in the estimating the parameter’s standard deviation, we summarize and compare the results.
