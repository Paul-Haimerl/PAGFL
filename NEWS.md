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
