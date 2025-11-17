# PAGFL v. 1.1.4

## R CMD check results

0 errors \| 0 warnings \| 2 notes

The two notes are:
(1) unable to verify current time
(2) tv_pagfl: no visible binding for global variable 'y'
    tv_pagfl: no visible binding for global variable 'p'
    tv_pagfl: no visible binding for global variable 'N'
    Undefined global functions or variables:
      N p y
      
The latter note stems from the legacy function `tv_pagfl`, which was renamed to `fuse_time` in the current version. 
To avoid hard breaks with older versions, `tv_pagfl` is now a wrapper of `fuse_time`. Subsequently, all arguments of `tv_pagfl` are identical to, and passed on to `fuse_time`.
Therefore, argument defaults of `tv_pagfl` which include 'y', 'p', and 'N' are only evaluated within `fuse_time`, triggering the note.

## Github CI tests passed:

-   windows-latest (devel)
-   windows-latest (oldrel)
-   windows-latest (release)
-   macOS-latest (release)
-   ubuntu-latest (devel)
-   ubuntu-latest (oldrel-1)
-   ubuntu-latest (release)

## Unit-test coverage

96%
