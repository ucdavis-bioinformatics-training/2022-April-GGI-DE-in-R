## Please complete the following before the course

* Install [R](https://urldefense.com/v3/__https:/cran.r-project.org/__;!!LQC6Cpwp!7awQqIxGimWN1BaG-W6M3P8S76YbSpCHSDrodTM_CWzHvtq0NAqYZ_4r_trAsYtkXOE$) and [RStudio](https://urldefense.com/v3/__https:/www.rstudio.com/products/rstudio/download/__;!!LQC6Cpwp!7awQqIxGimWN1BaG-W6M3P8S76YbSpCHSDrodTM_CWzHvtq0NAqYZ_4r_trAEp5qrJc$).  If you already have R installed, make sure you have an R version 4.0.0 or later.  Due to substantive changes in R between versions 3 and 4, some code in the course materials will not work on earlier versions. 

* Install necessary R packages by cutting and pasting the lines below into your RStudio console _one at a time_ (if you get a prompt about installing from source, say "no"):

```
install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("topGO")
BiocManager::install("KEGGREST")
BiocManager::install("gplots")
BiocManager::install("pathview")
BiocManager::install("Rgraphviz")
BiocManager::install("RColorBrewer")
```

Test the installations by cutting and pasting these lines into your RStudio console (don't worry if you get messages about an object being masked):

```
library(edgeR)
library(org.Mm.eg.db)
library(topGO)
library(KEGGREST)
library(gplots)
library(pathview)
library(Rgraphviz)
library(RColorBrewer)
```


* If you are new to R, please work through the introductory R materials [here](https://ucdavis-bioinformatics-training.github.io/2022_February_Introduction_to_R_for_Bioinformatics/index.html)  
