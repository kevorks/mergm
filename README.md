# Iterative Estimation of Mixed Exponential Random Graph Models with Nodal Random Effects

## Overview

This R package provides a function for fitting Exponential Random Graph Models (ERGMs) with nodal random effects. The estimation process combines the stepping algorithm for fixed parameters with pseudo-likelihood for nodal random effects. It is particularly useful for modeling network data with unobserved heterogeneity in the nodes.

## Installation

To install this package, you can use the `devtools` library in R:

```         
devtools::install_github("kevorks/mergm")
```

## Usage

The `mergm` function can be used to estimate mixed ERGMs with nodal random effects. Here are the main arguments:

-   **formula**: An R formula object representing the network model. The formula should be of the form `<network> ~ <model terms>`, where `<network>` is a network object, and `<model terms>` are ERGM terms. The term `sociality(1:nodes)` is automatically added to the formula as an offset.

-   **iter**: The number of iterations for the estimation process. Typically, 10 iterations are sufficient.

## Example

```         
library(mergm) 
data("zach", package = "ergm.count")

set.seed(2410)

mod <- mergm(zach ~ edges + kstar(2), 10)

gof(mod$model)

plot(gof(mod.model))

summary(mod$model)
```

## Prerequisites

```         
library(ergm)
library(network)
library(mgcv)
```

## Author

**Sevag Kevork** - *Author/Data Scientist* - [me\@sevagkevork.net](https://github.com/kevorks)

## References

For more details on the algorithm and methodology, you can refer to the following publication:

Kevork, S. and Kauermann, G. (2021). Iterative Estimation of Mixed Exponential Random Graph Models with Nodal Random Effects. Network Science 9, 478--498. <https://doi.org/10.1017/nws.2021.22>.

## License

This project is licensed under the [MIT License](https://chat.openai.com/c/LICENSE).

## Contribution

Contributions to the project are welcome. Feel free to submit issues or pull requests on the GitHub repository.
