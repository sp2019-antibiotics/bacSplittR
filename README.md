# bacSplittR

R package estimating ECOFF based on the Iterative Statistical Method

## Overview

This package is designed to find the Epidemiological Cutoff Value, ECOFF, of pairings of bacteria and antibiotics. Its use is to split bacteria into a resistant and a sensitive group for each antibiotic.
The algorithm implemented in the R package `bacSplittR` is an implementation for zone diameter data based on the Iterative Statistical Method introduced by Turnidge, Kahlmeter and Kronvall in 2006.


## Installation

Install the package using `devtools`.

```{r, eval = FALSE}

# install.packages("devtools")

devtools::install_github("sp2019-antibiotics/bacSplittR")

```
