###ScreenBEAM

[ScreenBEAM](https://github.com/jyyu/ScreenBEAM) is an [R](http://www.r-project.org) package to do gene-level meta-anlaysis of high-throughput functional genomics RNAi or CRISPR screening data. Both microarray and NGS data are supported. f you find this package useful, please cite the following paper:

**Yu J, Silva, J, Califano A. ScreenBEAM: a Novel Meta-Analysis Algorithm for Functional Genomics Screens via Bayesian Hierarchical Modeling. _Bioinformatics_ (In revision), 2015.**


#### Installation

Install ScreenBEAM from its
[GitHub repository](https://github.com/jyyu/ScreenBEAM). You first need to [devtools](https://github.com/hadley/devtools) package.

```r
install.packages(c("devtools"))
```

Install required Biobase package from Bioconductor.
```r
source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
```

Then install R/ScreenBEAM using the `install_github` function in the
[devtools](https://github.com/hadley/devtools) package.

```r
library(devtools)
install_github("jyyu/ScreenBEAM")
```

#### Input file format

Input file format for ScreenBEAM is tab-separated text file with "sh/sgRNA_Id", "gene id or symbol" as the first two columns, followed by the samples. Examples are below:

- [NGS count data] (https://github.com/jyyu/ScreenBEAM/blob/master/inst/extdata/NGS.example.tsv)
- [Microarray log2-transformed data] (https://github.com/jyyu/ScreenBEAM/blob/master/inst/extdata/microarray.example.tsv)


#### Example use

Try the following example to apply ScreenBEAM on NGS or microarry data. A data.frame with all statistical results is returned.

```r
library(ScreenBEAM)

###NGS data
r<-ScreenBEAM(

  ###input format
  input.file=system.file("extdata", "NGS.example.tsv", package = "ScreenBEAM")#tab-separted file
  ,
  control.samples=c('T0_A','T0_B','T0_C')#column names of control samples
  ,
  case.samples=c('T16_A','T16_B','T16_C')#column names of case/treated samples
  ,
  control.groupname='T0'#name your control group
  ,
  case.groupname='T16'#name your case group
  ,

  ###data pre-processing
  data.type='NGS'#data type
  ,
  do.normalization=TRUE
  ,
  filterLowCount=TRUE
  ,
  filterBy = 'control'
  ,
  count.cutoff=4
  ,

  ###Bayesian computing
  nitt=1500,#number of MCMC iterations, use small number here for testing, please use larger number in real data, 15000 is default
  burnin=500#number of burnin in MCMC sampling, 5000 is default

  )


###microarray data
r<-ScreenBEAM(

  ###input format
  input.file=system.file("extdata", "microarray.example.tsv", package = "ScreenBEAM")#tab-separted file
  ,
  control.samples=c('T0_A','T0_B','T0_C')#column names of control samples
  ,
  case.samples=c('T16_A','T16_B','T16_C')#column names of case/treated samples
  ,
  control.groupname='T0'#name your control group
  ,
  case.groupname='T16'#name your case group
  ,

  ###data pre-processing
  data.type='microarray'#data type
  ,
  do.normalization=FALSE#assuming the microarry data is normalized
  ,

  ###Bayesian computing
  nitt=1500,#number of MCMC iterations, use small number here for testing, please use larger number in real data, 15000 is default
  burnin=500#number of burnin in MCMC sampling, 5000 is default

)

head(r)

###save your results
write.csv(r,file=file.path('results.ScreenBEAM.csv'),row.names=FALSE,na='')

```

#### Google discussion group

https://groups.google.com/d/forum/screenbeam


#### Licenses

Licensed under the [MIT license](LICENCE). ([More information here](http://en.wikipedia.org/wiki/MIT_License).)

