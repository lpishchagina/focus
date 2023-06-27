<a id="top"></a>
#  focus Vignette
### Liudmila Pishchagina
### June 27, 2023

## Quick Start

ATTENTION: for this package it is necessary to install the library of C++  "Qhull" (see [this link](https://github.com/qhull/qhull))  

` focus ` focus is an R package written in Rcpp/C++ and developed to detect change using the Multidimensional Fast Online Changepoint Detection via Functional Pruning CUSUM statistics  method (MdFOCuS) in `p`-variate time series of length `n` from Exponential Natural family (Gaussian and Poisson Models).
A pruning strategy of our method is based  by building of covex hull on the set of potential candidates.

We present a basic use of the main functions of the `focus` package. 

We install the package from Github:

```r
#devtools::install_github("lpishchagina/focus")
library(focus)
```

## The function generate_ts 

The `generate_ts` is the generation of data (Gaussian or Poisson distribution) of dimension `p` with a given values of means and changes.

`type`  is the distribution (`gauss` or `poisson`).

`n`  is the time series length.

`p`  is the time series dimension.

`changes` is the changepoint vector that gives the last index of each segment.

The last element of `changes` is always less than to the length of time series.

By default, `changes = NULL` (for the data without changes). 

`means` is the matrix of successive means for the `p`-variate time series.

By default, `means = matrix(0, ncol = 1, nrow = p)`  for the Gaussian model and `means = matrix(1, ncol = 1, nrow = p)` for the Poisson model (for the data without changes). 

The length of each matrix row is equal to the length of `changes` plus one.

`noise` is a variance of the time series. By default, `noise = 1`(for Gaussian model).

```r

#parameters

set.seed(21)
N <- 100
P <- 2
Change <- N %/% 2
theta0 <-  rep(0, P)

```

```r

#Data generation

##the time series with one change

ts_gauss <- generate_ts(type = "gauss", p = P, n = N, changes = Change, means = matrix(c(theta0, theta0 + 5), nrow = P))
ts_poisson <- generate_ts(type = "poisson", p = P, n = N, changes = Change,  means = matrix(c(theta0 + 1, theta0 + 5), nrow = P))


##the time series with one change

ts_gauss0 <- generate_ts(type = "gauss", p = P, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = P))
ts_poisson0 <- generate_ts(type = "poisson", p = P, n = N, changes = NULL, means = matrix(1, ncol = 1, nrow = P))

```



## The function getChangePoints

The ` getChangePoints ` function returns the result of the segmentation of MdFOCuS -method.

` method ` is the type of the MdFOCuS method (`FOCuS0` or `method = FOCuS`). If `method = FOCuS0`(by default) the mean of the first segment is constant (zero).

` cost ` is the type of the cost function (` gauss `(by default) or `poisson `). 

` common_difference_step ` is the difference parameter for the steps of for the convex hull generation.

` common_ratio_step ` is the ratio parameter for the steps of for the convex hull generation.

` first_step_qhull ` is the first time for the convex hull generation (by default, for Gaussian model `t = first_step_qhull+p `, for Poisson model ` t = first_step_qhull + p + 10`).

` cand_nb ` is  the candidate number at the last iteration.

` opt_changes ` is the logical parameter (if ` opt_changes = TRUE ` get the optimal change at each iteration).

` opt_costs ` is the logical parameter (if ` opt_costs = TRUE ` get the optimal cost at each iteration).

` cands ` is the logical parameter (if ` cands = TRUE ` get the candidates at the last iteration).



```r


CPD <- list()

CPD[[1]] <- getChangePoints(ts_gauss, cands = TRUE, 
                                common_difference_step=1, 
                                common_ratio_step=2, 
                                first_step_qhull=2)

CPD[[2]] <- getChangePoints(ts_gauss, method ='FOCuS', cands = TRUE, 
                                common_difference_step=1, 
                                common_ratio_step=2, 
                                first_step_qhull=2)
                                
CPD[[3]] <- getChangePoints(ts_gauss0, cands = TRUE, 
                                common_difference_step=1, 
                                common_ratio_step=2, 
                                first_step_qhull=2)

CPD[[4]] <- getChangePoints(ts_gauss0, method ='FOCuS', cands = TRUE, 
                                common_difference_step=1, 
                                common_ratio_step=2, 
                                first_step_qhull=2)
                                
#Poisson model   

CPD[[5]] <- getChangePoints(ts_poisson, cost ="poisson", cands = TRUE, 
                                common_difference_step=1, 
                                common_ratio_step=2, 
                                first_step_qhull=2)

CPD[[6]] <- getChangePoints(ts_poisson, cost ="poisson", method ='FOCuS', cands = TRUE, 
                                common_difference_step=1, 
                                common_ratio_step=2, 
                                first_step_qhull=2)
                                
CPD[[7]] <- getChangePoints(ts_poisson0, cost ="poisson", cands = TRUE, 
                                common_difference_step=1, 
                                common_ratio_step=2, 
                                first_step_qhull=2)

CPD[[8]] <- getChangePoints(ts_poisson0, cost ="poisson", method ='FOCuS', cands = TRUE, 
                                common_difference_step=1, 
                                common_ratio_step=2, 
                                first_step_qhull=2)                               
                                
```

```r
#Change in 50

#Gaussian Model : 
#FOCuS0

CPD[[1]]$change
#[1] 50

CPD[[1]]$nb_candidates
#[1] 31

CPD[[1]]$candidates
# [1]  1  2  3  4  5 11 22 29 31 32 33 34 36 44 46 49 50 51 53 54 82 85 86 89 91 92 94 95 96 97 99

#FOCuS
CPD[[2]]$change
#[1] 50

CPD[[2]]$nb_candidates
#[1] 31

CPD[[2]]$candidates
# [1]  1  2  3  4  5 11 22 29 31 32 33 34 36 44 46 49 50 51 53 54 82 85 86 89 91 92 94 95 96 97 99

```

```r
#Change in 50

#Poisson Model : 
#FOCuS0

CPD[[5]]$change
#[1] 50

CPD[[5]]$nb_candidates
#[1] 27

CPD[[5]]$candidates
# [1] 1  2  3  4  5 12 14 16 23 24 33 45 48 49 50 51 52 54 55 72 89 90 94 95 96 98 99

#FOCuS
CPD[[6]]$change
#[1] 50

CPD[[6]]$nb_candidates
#[1] 27

CPD[[6]]$candidates
# [1] 1  2  3  4  5 12 14 16 23 24 33 45 48 49 50 51 52 54 55 72 89 90 94 95 96 98 99

```

`change` is the changepoint index.
`nb_candidates`  is  the candidate number at the last iteration.

`candidates` is the vector of the candidates at the last iteration.

[Back to Top](#top)
