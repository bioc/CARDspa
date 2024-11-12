# CARDspa

## Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics 
We developed a statistical method for spatially informed cell type deconvolution 
for spatial transcriptomics. Briefly,CARD is a reference-based 
deconvolution method that estimates cell type composition in spatial 
transcriptomics based on cell type specific expression information obtained 
from a reference scRNA-seq data. A key feature of CARD is its ability to 
accommodate spatial correlation in the cell type composition across tissue 
locations, enabling accurate and spatially informed cell type deconvolution 
as well as refined spatial map construction. CARD relies on an efficient 
optimization algorithm for constrained maximum likelihood estimation and 
is scalable to spatial transcriptomics with tens of thousands of spatial 
locations and tens of thousands of genes. CARD is implemented as an 
open-source R package, freely available at www.xzlab.org/software.html. 

Installation
------------
You can install the released version of CARDspa using BioConductor.

``` r
# install BiocManager if necessary
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CARDspa")

```

Or you can install using devtools in R.
``` r
# install devtools if necessary
install.packages('devtools')

# install the CARDspa package
devtools::install_github('https://github.com/YMa-lab/CARDspa')

# load package
library(CARDspa)

```

# Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure 
to raise issues with a detailed and reproducible exmple and also please 
provide the output of your sessionInfo() in R! 

How to cite `CARDspa`
-------------------
Ma, Y., Zhou, X. Spatially informed cell-type deconvolution for spatial 
transcriptomics. Nat Biotechnol 40, 1349–1359 (2022). 
https://doi.org/10.1038/s41587-022-01273-7

How to use `CARDspa`
-------------------
Please see vignettes.
``` r
browseVignettes("CARDspa")
```
