---
output:
  github_document:
    html_preview: false
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
options(tibble.print_min = 5, tibble.print_max = 5)
```

# dplyr <img src="man/figures/logo.png" align="right" />

## Overview


This package provides Gaussien Process with histogram intersection kernel and variance approximations
The package is design to provide the highest speed. It only encomprises a very limites number of functions and is supposed to be light weight.

## Installation

```{r, eval = FALSE}
# The easiest way to get gpHist is:
install.packages("gpHist")


# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("dennisthemenace2/gpHist")
```

If you encounter a clear bug, please file a minimal reproducible example.

## Usage

```{r, message = FALSE}
library(gpHist)


```

