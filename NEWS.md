# PAGFL 1.1.3 (development version)

* Added a progress bar displaying the percentage of convergence to `pagfl` and `tv_pagfl` if `verbose` is selected
* Slight changes in the formatting of the `print` method for `summary` of `pagfl` and `grouped_plm` objects
* Fix of missing regressor name in the output of `pagfl` and `grouped_plm` when having a single regressor and passing an empty formula
* Deprecation of `DGP` argument in `sim_tv_DGP`
* Bugfix when passing a prescribed coefficient matrix and selecting `dynamic = TRUE` in `sim_DGP`
* Bugfix to preclude explosive processes when selecting `dynamic = TRUE` in `sim_tv_DGP`
* Added seeds to all function examples
* Support for RcppArmadillo 14.4.0

# PAGFL 1.1.2 (stable version)

* Fixed backwards compatibility issue with the generic functions `fitted` and `resid`
* Bugfix when passing index variables and an empty formula `y ~ .` for `tv_pagfl` and `grouped_tv_plm`
* Improved documentation
* Changed x-axis label in the plot produced by calling `summary()` of a `tvpagfl` object
* Bugfix in the plot produced by calling `fitted()` and passing a character time-index variable
* Added the current algorithm iteration and tuning parameter as a progress counter to the console for `pagfl` and `tv_pagfl` if option `verbose` is selected
* Changed the IC selecting the best fitting tuning parameter for the `tv_pagfl` and `grouped_tv_plm` procedures to include the logarithmic mean squared error as the fitness term

# PAGFL 1.1.1

* Introduction of `grouped_plm` and `grouped_tv_plm` to estimate grouped (time-varying) panel data models given an exogenous group structure
* Remove warm-starts when iterating across different tuning parameters
* Small efficiency upgrades
* More extensive unit testing
* Improved documentation
* Small updates to generic methods

# PAGFL 1.1.0

* Introduction of the time-varying *PAGFL* `tv_pagfl`
* Introduction of s3 methods (`summary()`, `coef()`, `fitted()`, `resid()`) for the output of `pagfl` and `tv_pagfl`
* Change of the user interface to formula objects
* Implementation of unit testing
* Renamed functions to be consistently snake case
* Support for unordered and/ or unbalanced panel data sets via `index`
* Possibility to estimate a mix of time-constant and time-varying coefficients in the same panel data model (`const_coef`)
* Added row and column names to the estimation output
* Maximum within-group heterogeneity `group_tol` set to a machine inaccuracy value
* Improved documentation, checks, and error/warning messages
* Enabled 64bit C++ compiler to support large datasets
* Improvements to efficiency and numerical stability
* Bugfixes

# PAGFL 1.0.1

* Patch to fix freak bug in the M1Mac build

# PAGFL 1.0.0

* Initial CRAN submission.
