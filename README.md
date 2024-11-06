# rOpenScPCA

This repository contains the `rOpenScPCA` package, which offers  utility functions to support single-cell RNAseq analysis in the [OpenScPCA project](https://openscpca.readthedocs.io/en/latest/).

## Installation

`rOpenScPCA` can either be installed with the `remotes` package, or with `renv` if you need to track it in an `renv.lock` file:

```r
# Install the package with remotes
remotes::install_github("AlexsLemonade/rOpenScPCA")

# Install the package with renv
renv::install("AlexsLemonade/rOpenScPCA")
# You can then add to a renv.lock file with renv::snapshot()
```

<!--
## Usage

[FORTHCOMING]

We have compiled example notebooks using this package:

- The `OpenScPCA-analysis` module `hello-clusters` demonstrates how to perform and evaluate clustering with `rOpenScPCA`
-->
