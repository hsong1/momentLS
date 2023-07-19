
<!-- README.md is generated from README.Rmd. Please edit that file -->

# momentLS

<!-- badges: start -->
<!-- badges: end -->

The goal of momentLS is to implement Moment LS estimators in “Efficient
shape-constrained inference for the autocovariance sequence from a
reversible Markov chain” <https://arxiv.org/abs/2207.12705>

## Installation

You can install the development version of momentLS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hsong1/momentLS")
```

## Example

We provide a basic usage example of momentLS package using a simulated
AR(1) chain:

``` r
library(momentLS)

# generate an AR(1) chain of length M with rho =.5
set.seed(12345); M=10000
chainParams = list(type="AR", rho = 0.5, M = M)
ch = generateChain(chainParams)

# compute the empirical autocovariances
r=autocov(ch$x)

# fit MomentLS
delta_hat = tune_delta(r,method="ft",seq_type="average") # tune delta
m =SR1(r,delta = delta_hat[2]) # fit
```

### Asymptotic variance estimator comparison

``` r
# compute the MomentLS asymptotic variance estimator
MomentLS_avar = asympVariance(weights = m$weights, support = m$support) 

# comparison with  BM, OLBM, init-seq estimators
result = 
  c("Truth"=ch$varTruth,
    "MomentLS"=MomentLS_avar,
    "BM"=M*mcmcse::mcse(ch$x,method="bm")$se^2, 
    "OLBM"=M*mcmcse::mcse(ch$x,method="obm")$se^2,
    "Init-pos"=mcmc::initseq(ch$x)$var.pos,
    "Init-dec"=mcmc::initseq(ch$x)$var.dec,
    "Init-conv"=mcmc::initseq(ch$x)$var.con
  )

knitr::kable(t(result))
```

| Truth | MomentLS |      BM |   OLBM | Init-pos | Init-dec | Init-conv |
|------:|---------:|--------:|-------:|---------:|---------:|----------:|
|     4 | 3.970975 | 4.37627 | 4.3261 | 3.926424 | 3.888211 |  3.884543 |

### Comparison of empirical and MomentLS autocovariances

``` r
# first 100 estimated autocovariances
estimated_autocov = computeMoments(support = m$support, weights = m$weights,M = 100)
plot(0:99, r[1:100], col="lightblue", xlab="lags", ylab = "autocovariances", 
     ylim= c(min(r[1:100], estimated_autocov),max(r[1:100], estimated_autocov)))
points(0:99, estimated_autocov,pch=3)
```

<img src="man/figures/README-autocov_seq-1.png" width="50%" height="50%" style="display: block; margin: auto;" />
