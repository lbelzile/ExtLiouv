# ExtLiouv
## Miscellaneous functions for `Extremal attractors of Liouville copulas`


Some utilities for pairwise composite likelihood estimation of unstructured Husler-Reiss and extremal Student models (both are experimental until further notice).

- `fmvcpot` performs pairwise composite likelihood estimation and returns a list with estimated values. It handles missing values in the observation matrix.
- `compositemat`: Godambe information matrix and the variability matrix can be obtained by means of a nonparametric boostrap using the function .

Also included are exponent measures for the scaled Dirichlet extreme value distribution, 
- `fscore` for gradient score estimation for the class of scaled Dirichlet models, 
- `specdens` for the spectral density
- `isAbove`, a function that returns a logical matrix indicating whether observations exceed marginally a given threshold vector


### Installation 
The package is currently under development and is NOT available from CRAN. 
The current syntax is subject to changes and some functions may be ported to the package `mev` in due time.

To install from Github, run

```R
devtools::install_github("lbelzile/ExtLiouv")
```

after installing `devtools`.
