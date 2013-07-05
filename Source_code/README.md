### compile source code

require [Rtools](http://cran.r-project.org/bin/windows/Rtools/)
#### 1. Windows
use compile.bat

#### 2. Linux
use commend: 
```
R CMD SHLIB -c --preclean FMM-HMM.c nrutil.c
```
