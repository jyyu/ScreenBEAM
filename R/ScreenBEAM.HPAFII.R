# FileName: T14vsT0.shRNALevel.R
# 
# Author: Jiyang Yu
# Time: 3:17:32 PM
# Date: Jul 19, 2011
#
# Description: 
###############################################################################


source('comm.inc.R')

source(file.path(inc.path,'ScreenBEAM.inc.R'))

source(file.path(inc.path,'ScreenBEAM_plot.inc.R'))

comm.path<-'Achilles/V2.4.3/selectedLines/HPAFII'

	
####read data
eset<-generateEset(
		input.file=file.path(data.path,comm.path,'HPAFII.countPerM.txt')
		,
		control.samples=c('pDNA.p0.10025','pDNA.p0.10025_10056','pDNA.p0.110625','pDNA.p0.110709','pDNA.p0.110709_10056','pDNA.p0.54_10051')
		,
		case.samples=c('HPAFII_PANCREAS.p5.A','HPAFII_PANCREAS.p5.B','HPAFII_PANCREAS.p5.C','HPAFII_PANCREAS.p5.D')
		,
		control.groupname='T0',
		case.groupname='HPAFII.T5',
		gene.columnId=2
)

table(pData(eset)$condition)
table(pData(eset)$group)


#QC
QC.expresso(eset,outputdir = file.path(data.path,comm.path,'QC1'),intgroup='group',do.logtransform=TRUE)


#correlation of replicates
####no filtering
plotWithinGroupCorrelation(
		eset,
		plot.dir=file.path(data.path,comm.path,'QC2.replicates'),
		group.column='group'
		,
		filterLowCount=FALSE
)

Sys.time()
t.start<-Sys.time()

DR.df<-DRAgeneLevel(
		eset=eset
		,
		do.log2=TRUE
		,
		do.restand=TRUE
		,
		filterLowCount=FALSE
		,
		filterBy = 'control'
		,
		count.cutoff=32
)

t.end<-Sys.time()

#time
cat(t.end-t.start,'\n')
# 4.487811 


dim(eset)
# Features  Samples 
#    46874       10 


nrow(subset(fData(eset),!duplicated(gene)))
# [1] 11359


write.csv(DR.df,file=file.path(data.path,comm.path,'DR.ScreenBEAM.csv'),row.names=FALSE,na='')

