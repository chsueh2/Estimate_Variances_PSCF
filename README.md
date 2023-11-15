# Estimate Variances of Model Parameters Using Perturbed SSE Curve Fitting (PSCF) Method

Perturbed SSE Curve Fitting (PSCF) Method proposed by Prof. Robert Hayes (Department of Nuclear Engineering, NCSU) is used to estimate variances of model parameters. In this work, we implement  PSCF algorithm  and compare it with the other common statistics methods including nonlinear regression models, bootstrap confidence intervals and the delta normality method.

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

Consider three different hospitals and each hospital has patients that end up with infections. We are interested whether or not the multinomials are homogenous across the hospitals.

Conduct chi-square tests for homogeneity using R. Also manually calculate the LRT statistic, the Pearson Chi-square statistic, the critical value, and find approximate p-values for the hypotheses using both test statistics.

- Determine how well the asymptotic rejection region performs at controlling $\alpha$
- Determine the power of the asymptotic test when comparing certain alternative situations

## Part 1 - Data Example

Create matrix for the given hospital data and conduct chi-square tests for homogeneity. Implement the  test in R and verify the result with R build-in test function `chisq.test()`.

## Part 2 - Derivation

Derive the likelihood ratio test by considering the generic case of comparing $J$ independent multinomial distributions, each with $I$ categories. 

## Part 3 - Simulation

The Pearson statistic can be derived as a Taylor series approximation to the LRT. Investigate the $\alpha$ control of the Pearson chi-square test and its power. 

Consider the following three different setups:

1. Two multinomial case only, where each multinomial has three categories
1. All combinations of four sample sizes for each multinomial (16 total cases)
1. Three different probabilities that may generate a particular multinomial

## Part 4 - Summarization 

Analyze and summarize the results using faceted line plots for easy comparison.
