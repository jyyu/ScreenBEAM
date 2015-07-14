#ScreenBEAM Running examples

library(ScreenBEAM)

#NGS data
r<-ScreenBEAM(

  ###input format
  input.file=system.file("extdata", "NGS.example.tsv", package = "ScreenBEAM")
  ,
  control.samples=c('T0_A','T0_B','T0_C')
  ,
  case.samples=c('T16_A','T16_B','T16_C')
  ,
  control.groupname='T0'
  ,
  case.groupname='T16'
  ,

  ###data pre-processing
  data.type='NGS'
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
  input.file=system.file("extdata", "microarray.example.tsv", package = "ScreenBEAM")
  ,
  control.samples=c('T0_A','T0_B','T0_C')
  ,
  case.samples=c('T16_A','T16_B','T16_C')
  ,
  control.groupname='T0'
  ,
  case.groupname='T16'
  ,

  ###data pre-processing
  data.type='microarray'
  ,
  do.normalization=FALSE
  ,

  ###Bayesian computing
  nitt=1500,#number of MCMC iterations, use small number here for testing, please use larger number in real data, 15000 is default
  burnin=500#number of burnin in MCMC sampling, 5000 is default

)

head(r)


###save your results
write.csv(r,file=file.path('results.ScreenBEAM.csv'),row.names=FALSE,na='')

