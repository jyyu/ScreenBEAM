###ScreenBEAM: 

[ScreenBEAM](https://github.com/jyyu/ScreenBEAM) is an [R](http://www.r-project.org) package to do a gene-level meta-anlaysis of high-throughput functional genomics RNAi or CRISPR screening data. Both microarray and NGS data are supported. f you find this package useful, please cite the following paper:

Yu J, Silva, J, Califano A. ScreenBEAM: a Novel Meta-Analysis Algorithm for Functional Genomics Screens via Bayesian Hierarchical Modeling. Bioinformatics (In revision), 2015.


#### Installation

Install R/qtlcharts from its
[GitHub repository](https://github.com/jyyu/ScreenBEAM). You first need to [devtools](https://github.com/hadley/devtools) package.

```r
install.packages(c("devtools"))
```

Then install R/qtlcharts using the `install_github` function in the
[devtools](https://github.com/hadley/devtools) package. (With
`build_vignettes=TRUE`, the vignettes will be built and installed.)

```r
library(devtools)
install_github("kbroman/qtlcharts", build_vignettes=TRUE)
```

#### Example use

Try the following example, which creates an interactive chart with LOD
curves linked to estimated QTL effects.

```r
library(ScreenBEAM)

```

#### Licenses

Licensed under the [MIT license](LICENSE). ([More information here](http://en.wikipedia.org/wiki/MIT_License).)


