# Author: Jiyang Yu
###############################################################################
#####
#diff representation analyis at gene level
#eset
#condition: treatment group names (must be in pData(eset)$condition)
#ctrl: ctrl grou names (paired with condition, must be in pData(eset)$condition)
#do.restand: whether do restand in the analysis
#family: model to use, gaussian, poisson
#estimation.method: Bayesian or MLE
DRAgeneLevel<-function(eset,data.type=c('microarray','NGS'),do.normalization=FALSE,filterLowCount=TRUE,filterBy='control',count.cutoff=4,nitt=15000,burnin=5000,...){

  if(missing(data.type)){
    data.type<-ifelse(all(exprs(eset)>0) & is.integer(exprs(eset)),'NGS','microarray')
  }

  data.type<-match.arg(data.type,c('microarray','NGS'))
  cat(paste('data.type:',data.type,'data\n'))

  #do.log2: whether to do log2 transformation for the data
  do.log2<-ifelse(data.type=='NGS',TRUE,FALSE)

  ###normalization
  if(do.normalization){
    if(data.type=='NGS'){
      cat('normalization: scale normlaization to NGS count data!\n')
      pseudoCount<-ifelse(all(exprs(eset)>0),0,1)
      total<-ifelse(all(apply(exprs(eset),2,sum)<1e6),1e6,NULL)
  exprs(eset)<-as.matrix(normalize.scale(exprs(eset),total=total,pseudoCount = pseudoCount))
    }

    if(data.type=='microarray'){
#quantile normalization
      cat('normalization: quantile normlaization to microarray log2(intensity) data!\n')
      exprs(eset)<-normalize.quantile(exprs(eset), ties=TRUE)
    }
  }


	if(sum(names(fData(eset))=='gene')==0){
		if(sum(names(fData(eset))=='geneSymbol')==1){
			names(fData(eset))[names(fData(eset))=='geneSymbol']<-'gene'
		}else{
			stop('gene column doesn\'t exist!\n')
		}
	}

  DR<-data.frame(gene=unique(fData(eset)$gene))

  condition<-'treatment'
	ctrl<-'control'

	if(length(condition)!=length(ctrl))
		stop('condion and ctrl have diff length!\n')

	if(!all(c(condition,ctrl)%in%pData(eset)$condition))
		stop('conditon, ctrl are not all in pData(eset)$condition!\n')

	#add n.shRNA column
	n.shRNAs<-table(fData(eset)$gene)
	n.shRNAs<-data.frame(gene=names(n.shRNAs),n.sh_sgRNAs.raw=as.integer(n.shRNAs))
	DR<-merge(n.shRNAs,DR,by='gene')


	eset.sel<-eset

  if(data.type=='NGS'){
	if(filterLowCount){
		sel<-grep(filterBy,pData(eset)$condition)
		if(length(sel)>=1){
			eset.sel<-eset[apply(exprs(eset[,sel]),1,median)>=count.cutoff,]
		}
	}
}

#i=1

	for(i in 1:length(condition)){

		cat(condition[i],'vs.', ctrl[i], 'in processs...\n')

		d.eset<-eset.sel[,pData(eset.sel)$condition %in% c(ctrl[i],condition[i])]

		condition.tag<-pData(d.eset)$group[pData(d.eset)$condition%in%condition[i]][1]
		ctrl.tag<-pData(d.eset)$group[pData(d.eset)$condition%in%ctrl[i]][1]

		if(do.log2){
			#take log2 for count data
			exprs(d.eset)<-log2(exprs(d.eset))
		}

		comp<-factor(gsub(condition[i],1,gsub(ctrl[i],0,pData(d.eset)$condition)))

		d<-data.frame(gene=as.character(fData(d.eset)$gene),exprs(d.eset),stringsAsFactors=FALSE)

		dr<-plyr::ddply(d,'gene','combRowEvid.2grps',comp=comp,family=gaussian,method='Bayesian',nitt=nitt,burnin=burnin,thin=1,logTransformed=do.log2,restand=FALSE,pooling='partial',...)

		if(do.log2){
			dr$AveSignal<-round(2^dr$AveSignal,0)
		}else{
			dr$AveSignal<-round(dr$AveSignal,0)
		}

		#FDR.partial<-fdrtool(dr$pval.partial,'pvalue',plot=FALSE, verbose=FALSE)$qval
		FDR.BH.partial<-p.adjust(dr$pval.partial,'BH')
		partial.id<-grep('pval.partial',names(dr))
		dr<-data.frame(dr[,1:partial.id],FDR.BH.partial=FDR.BH.partial,dr[,(partial.id+1):ncol(dr)])

		names(dr)<-gsub('.partial','',names(dr))

		m<-min(dr$pval[dr$pval>0])[1]
		dr$pval[dr$pval==0]<-ifelse(m<1e-7,m,1e-7)
		dr$z<-sign(dr$t)*qnorm(dr$pval/2,lower.tail = FALSE)

		#transform FC to log2FC
		dr$FC<-sign(dr$FC)*log2(abs(dr$FC))
		names(dr)<-gsub('^FC','log2FC',names(dr))

		#ctrl.short<-ifelse(is.null(ctrl.tag),gsub('.*(.*)\\.(.*)','\\2',ctrl[i]),ctrl.tag[i])
		names(dr)[-1]<-paste(names(dr)[-1],'.',sum(comp==1),'VS',sum(comp==0),'.',condition.tag,'_VS_',ctrl.tag,sep='')

		#DR<-merge(DR,dr,by.x='gene',by.y='gene',all.x=TRUE)
		DR<-merge(DR,dr,by.x='gene',by.y='gene')
	}


	names(DR)<-gsub('^FDR.BH','FDR',gsub('^sd','B.sd',gsub('^coef','B',gsub('n.levels','n.sh_sgRNAs.passFilter',names(DR)))))

	col.sel<-c('n.sh_sgRNAs.passFilter','B','z','pval','FDR','B.sd')

	condition<-gsub('z.','',grep('^z',names(DR),value=T));condition

	col.pre<-c(
			###annotation columns
			#names(DR)[1:(grep('log2FC',names(DR))[1]-1)]
		'gene'
		,
		'n.sh_sgRNAs.raw'
		,
			#col.sel
			#DR columns
			paste(rep(col.sel,each=length(condition)),rep(condition,length(col.sel)),sep='.')
	)

	DR<-DR[,
			c(col.pre
							#,
							#setdiff(names(DR),col.pre)
					)
			]

			DR<-DR[order(DR[,grep('^B\\.',names(DR))]),]

			DR

}


#d: dataframe as below: (first column is the annotation)
#			geneSymbol  DMSO_1 DMSO_2 DMSO_3    DEX_1    DEX_2    DEX_3
#P0112072    PRKAR2B 0.50099 1.2108 1.0524 -0.34881 -0.13441 -0.87112
#P0112073    PRKAR2B 1.84579 2.0356 2.6025  1.62954  1.88281  1.29604
#comp as below:
# [1] 0 0 0 1 1 1
# Levels: 0 1
#d<-d[1:10,]
#' @export
combRowEvid.2grps<-function(d,comp,
		family=gaussian
		,method=c('MLE','Bayesian'),pooling=c('partial'),

		n.iter=1000
		,
		prior.V.scale=0.02
		,
		prior.R.nu=1
		,
prior.G.nu=2
,
nitt = 13000
,
burnin =3000
,
thin=10
,
restand=FALSE
,
logTransformed=TRUE
,
log.base=2
,
average.method=c('geometric')
,
pseudoCount=0
){

	if(!all(comp %in% c(1,0))){
		stop('comp only takes 1 or 0 !!! \n')
	}

	if(missing(pooling))
		pooling<-c('full','no','partial')
	else
		pooling<-tolower(pooling)

	pooling<-match.arg(pooling,several.ok=T)


	if (is.character(family))
		family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family))
		family <- family()
	if (!family$family %in% c('gaussian','binomial', 'poisson')) {
		print(family)
		stop("Only Gaussian Poisson, and Binomial are supported!!! \n")
	}

	if(missing(method))
		method<-'Bayesian'
	if(grepl('Bayes',method,ignore.case = T)) method<-'Bayesian'
	if(grepl('MLE|MaxLikelihood',method,ignore.case = T)) method<-'MLE'
	method<-match.arg(method)

	#d<-d[1,]
	#cat(unique(d[,1]),'\t',dim(d),'\n')
	d<-data.frame(t(d[,-1]))
	#cat(dim(d),'\n')


	dat<-data.frame(
			response=
					c(unlist(d[comp==levels(comp)[1],]),unlist(d[comp==levels(comp)[2],]))
			,
			treatment=
					factor(c(rep(levels(comp)[1],sum(comp==levels(comp)[1])*ncol(d)),rep(levels(comp)[2],sum(comp==levels(comp)[2])*ncol(d))))
			,
			probe=
					factor(c(rep(colnames(d),each=sum(comp==levels(comp)[1])),rep(colnames(d),each=sum(comp==levels(comp)[2]))))
	)


#calculate FC
	FC.val<-FC(dat$response,dat$treatment,logTransformed=logTransformed,log.base=log.base,average.method='geometric',pseudoCount=pseudoCount)

	AveSignal<-mean(dat$response)
	n.levels<-nlevels(dat$probe)


	rs<-c(FC=FC.val,AveSignal=AveSignal,n.levels=as.integer(n.levels))

	#cat(rs,'\n')

	if('arithmetic' %in% average.method){
		FC.ari<-FC(dat$response,dat$treatment,logTransformed=logTransformed,log.base=log.base,average.method='arithmetic',pseudoCount=0)
		rs<-c(FC.ari_raw=FC.ari,rs)
	}


#re-standarize the input data
	if(restand & sd(dat$response)>0)
		dat$response<-0.5*(dat$response-mean(dat$response))/sd(dat$response)

	#MLE approach to estimate parameters
	if(method=='MLE'){
		if(family$family=='gaussian'){
			#comlete pooling
			if('full'%in%pooling){
				M.full<-glm(response ~ treatment, data=dat)
				sum.tmp<-summary(M.full)
				rs.full<-c(
						coef.full=sum.tmp$coef[2,1],
						se.full=sum.tmp$coef[2,2],
						t.full=sum.tmp$coef[2,3],
						pval.full=sum.tmp$coef[2,4],
						#z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
						df.full=sum.tmp$df.residual,
						AIC.full=sum.tmp$aic,
						BIC.full=AIC(M.full,k=log(nrow(dat))),
						Dev.full=sum.tmp$deviance
				)
				rs<-c(rs,rs.full)
			}

			#no pooling
			if('no'%in%pooling){
				if(n.levels>1)
					M.no<-glm(response ~ treatment + probe, data=dat)
				else
					M.no<-glm(response ~ treatment, data=dat)

				sum.tmp<-summary(M.no)
				rs.no<-c(
						coef.no=sum.tmp$coef[2,1],
						se.no=sum.tmp$coef[2,2],
						t.no=sum.tmp$coef[2,3],
						pval.no=sum.tmp$coef[2,4],
						#z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
						df.no=sum.tmp$df.residual,
						AIC.no=sum.tmp$aic,
						BIC.no=AIC(M.no,k=log(nrow(dat))),
						Dev.no=sum.tmp$deviance
				)
				rs<-c(rs,rs.no)
			}


			#partial pooling with multilevel model of varing slopes and intercepts
			if('partial'%in%pooling){
				if(n.levels>1){
					#consider sd==0 case
					if(sd(dat$response)==0)
					{
						dat$response<-rnorm(nrow(dat),mean(dat$response),sd(dat$response)+0.001)
					}
					M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat)
					sum.tmp<-summary(M.partial)
					t.partial<-sum.tmp$coef[2,3]
					df.partial<-as.numeric(	sum.tmp$devcomp$dims['n']-	sum.tmp$devcomp$dims['p']-	sum.tmp$devcomp$dims['q'])
					pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
					z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
					#########################################
					#Another way to calculate pvalue
					#########################################
					#M.null<-glmer(response ~ (treatment + 1 | probe), data=dat)
					#anova(M.partial, M.null)
					rs.partial<-c(
							coef.partial=sum.tmp$coef[2,1],
							se.partial=sum.tmp$coef[2,2],
							t.partial=t.partial,
							pval.partial=pval.partial,
							z.partial=z.partial,
							df.partial=df.partial,
							AIC.partial=as.numeric(sum.tmp$AICtab[1]),
							BIC.partial=as.numeric(sum.tmp$AICtab[2]),
							REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
							logLik=-as.numeric(sum.tmp$AICtab[3]),
							Dev.partial=-as.numeric(sum.tmp$AICtab[4])
					)
				}else{
					M.partial<-glm(response ~ treatment, data=dat)
					sum.tmp<-summary(M.partial)
					t.partial<-sum.tmp$coef[2,3]
					df.partial<-sum.tmp$df.residual
					pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
					z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
					rs.partial<-c(
							coef.partial=sum.tmp$coef[2,1],
							se.partial=sum.tmp$coef[2,2],
							t.partial=t.partial,
							pval.partial=pval.partial,
							z.partial=z.partial,
							df.partial=df.partial,
							AIC.partial=sum.tmp$aic,
							BIC.partial=AIC(M.partial,k=log(nrow(dat))),
							REMLDev.partial=sum.tmp$deviance,
							logLik=logLik(M.partial),
							Dev.partial=sum.tmp$deviance
					)
				}
				rs<-c(rs,rs.partial)
			}
		}else if(family$family=='binomial'){
			if('full'%in%pooling){
				M.full<-glm(treatment ~ response, data=dat, family=family)
				sum.tmp<-summary(M.full)
				rs.full<-c(
						coef.full=sum.tmp$coef[2,1],
						se.full=sum.tmp$coef[2,2],
						z.full=sum.tmp$coef[2,3],
						pval.full=sum.tmp$coef[2,4],
						df.full=sum.tmp$df.residual,
						AIC.full=sum.tmp$aic,
						BIC.full=AIC(M.full,k=log(nrow(dat))),
						Dev.full=sum.tmp$deviance
				)
				rs<-c(rs,rs.full)
			}
			if('no'%in%pooling){
				if(n.levels>1)
					M.no<-glm(treatment ~ response + probe, data=dat, family=family)
				else
					M.no<-glm(treatment ~ response, data=dat, family=family)

				sum.tmp<-summary(M.no)
				rs.no<-c(
						coef.no=sum.tmp$coef[2,1],
						se.no=sum.tmp$coef[2,2],
						pval.no=sum.tmp$coef[2,4],
						z.no=sum.tmp$coef[2,3],
						df.no=sum.tmp$df.residual,
						AIC.no=sum.tmp$aic,
						BIC.no=AIC(M.no,k=log(nrow(dat))),
						Dev.no=sum.tmp$deviance
				)
				rs<-c(rs,rs.no)
			}

			if('partial'%in%pooling){
				if(n.levels>1){
					M.partial<-glmer(treatment ~ response + (response + 1 | probe), data=dat, family=family)
					sum.tmp<-summary(M.partial)
					t.partial<-sum.tmp$coef[2,3]
					df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
					pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
					z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
					rs.partial<-c(
							coef.partial=sum.tmp$coef[2,1],
							se.partial=sum.tmp$coef[2,2],
							t.partial=t.partial,
							pval.partial=pval.partial,
							z.partial=z.partial,
							df.partial=df.partial,
							AIC.partial=as.numeric(sum.tmp$AICtab[1]),
							BIC.partial=as.numeric(sum.tmp$AICtab[2]),
							Dev.partial=as.numeric(sum.tmp$AICtab[4])
					)
				}else{
					M.partial<-glm(treatment ~ response, data=dat, family=family)
					sum.tmp<-summary(M.partial)
					t.partial<-sum.tmp$coef[2,3]
					df.partial<-sum.tmp$df.residual
					pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
					z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
					rs.partial<-c(
							coef.partial=sum.tmp$coef[2,1],
							se.partial=sum.tmp$coef[2,2],
							t.partial=t.partial,
							pval.partial=pval.partial,
							z.partial=z.partial,
							df.partial=df.partial,
							AIC.partial=sum.tmp$aic,
							BIC.partial=AIC(M.partial,k=log(nrow(dat))),
							Dev.partial=sum.tmp$deviance
					)
				}
				rs<-c(rs,rs.partial)
			}


		}else if(family$family=='poisson'){
			#comlete pooling
			if('full'%in%pooling){
				M.full<-glm(response ~ treatment, data=dat,family='poisson')
				sum.tmp<-summary(M.full)
				rs.full<-c(
						coef.full=sum.tmp$coef[2,1],
						se.full=sum.tmp$coef[2,2],
						t.full=sum.tmp$coef[2,3],
						pval.full=sum.tmp$coef[2,4],
						#z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
						df.full=sum.tmp$df.residual,
						AIC.full=sum.tmp$aic,
						BIC.full=AIC(M.full,k=log(nrow(dat))),
						Dev.full=sum.tmp$deviance
				)
				rs<-c(rs,rs.full)
			}

			#no pooling
			if('no'%in%pooling){
				if(n.levels>1)
					M.no<-glm(response ~ treatment + probe, data=dat,family='poisson')
				else
					M.no<-glm(response ~ treatment, data=dat,family='poisson')

				sum.tmp<-summary(M.no)
				rs.no<-c(
						coef.no=sum.tmp$coef[2,1],
						se.no=sum.tmp$coef[2,2],
						t.no=sum.tmp$coef[2,3],
						pval.no=sum.tmp$coef[2,4],
						#z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
						df.no=sum.tmp$df.residual,
						AIC.no=sum.tmp$aic,
						BIC.no=AIC(M.no,k=log(nrow(dat))),
						Dev.no=sum.tmp$deviance
				)
				rs<-c(rs,rs.no)
			}

			#partial pooling with multilevel model of varing slopes and intercepts
			if('partial'%in%pooling){
				if(n.levels>1){
					M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat,family='poisson')
					sum.tmp<-summary(M.partial)
					t.partial<-sum.tmp$coef[2,3]
					df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
					pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
					z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
					rs.partial<-c(
							coef.partial=sum.tmp$coef[2,1],
							se.partial=sum.tmp$coef[2,2],
							t.partial=t.partial,
							pval.partial=pval.partial,
							z.partial=z.partial,
							df.partial=df.partial,
							AIC.partial=as.numeric(sum.tmp$AICtab[1]),
							BIC.partial=as.numeric(sum.tmp$AICtab[2]),
							#REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
							logLik=-as.numeric(sum.tmp$AICtab[3]),
							Dev.partial=-as.numeric(sum.tmp$AICtab[4])
					)
				}else{
					M.partial<-glm(response ~ treatment, data=dat,family='poisson')
					sum.tmp<-summary(M.partial)
					t.partial<-sum.tmp$coef[2,3]
					df.partial<-sum.tmp$df.residual
					pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
					z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
					rs.partial<-c(
							coef.partial=sum.tmp$coef[2,1],
							se.partial=sum.tmp$coef[2,2],
							t.partial=t.partial,
							pval.partial=pval.partial,
							z.partial=z.partial,
							df.partial=df.partial,
							AIC.partial=sum.tmp$aic,
							BIC.partial=AIC(M.partial,k=log(nrow(dat))),
							#REMLDev.partial=sum.tmp$deviance,
							logLik=logLik(M.partial),
							Dev.partial=sum.tmp$deviance
					)
				}
				rs<-c(rs,rs.partial)
			}

		}
		else{
			stop('Only liner model with gaussian or poisson distrn and binomial family model are supported !!! \n')
		}
	}

	#Bayesian approach
	else if(method=='Bayesian'){
	  if('full'%in%pooling | 'no'%in%pooling){
	  require(arm)
	  }else{
	  require(MCMCglmm)
}
		if(family$family=='gaussian'){
			if('full'%in%pooling){

        #comlete pooling
				M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter)
				sum.tmp<-summary(M.full)

				rs.full<-c(
						coef.full=sum.tmp$coef[2,1],
						se.full=sum.tmp$coef[2,2],
						t.full=sum.tmp$coef[2,3],
						pval.full=sum.tmp$coef[2,4],
						#z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
						df.full=sum.tmp$df.residual-sum.tmp$df[1],
						AIC.full=sum.tmp$aic,
						BIC.full=AIC(M.full,k=log(nrow(dat))),
						Dev.full=sum.tmp$deviance
				)
				rs<-c(rs,rs.full)
			}

			if('no'%in%pooling){
				#no pooling
				if(n.levels>1)
					M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter)
				else
					M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter)
				sum.tmp<-summary(M.no)

				rs.no<-c(
						coef.no=sum.tmp$coef[2,1],
						se.no=sum.tmp$coef[2,2],
						t.no=sum.tmp$coef[2,3],
						pval.no=sum.tmp$coef[2,4],
						#z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
						df.no=sum.tmp$df.residual-sum.tmp$df[1],
						AIC.no=sum.tmp$aic,
						BIC.no=AIC(M.no,k=log(nrow(dat))),
						Dev.no=sum.tmp$deviance
				)
				rs<-c(rs,rs.no)

			}

			if('partial'%in%pooling){

        #partial pooling with multilevel model of varing slopes and intercepts
				prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))

				if(n.levels>1){
					M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
					df.partial<-nrow(dat)-(n.levels+1)*2
					#if(df.partial<2)
					#	df.partial<-2
				}else{
					prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu))

					M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)

					df.partial<-nrow(dat)-2
					#if(df.partial<2)
					#	df.partial<-2
					#cat(n.levels,df.partial,'\n')
				}
				sum.tmp<-summary(M.partial)
				t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
				pval.partial<-2*pt(-abs(t.partial),df=df.partial)

				if(is.na(pval.partial))
					pval.partial<-sum.tmp$sol[2,5]
				z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
				rs.partial<-c(
						coef.partial=sum.tmp$sol[2,1]
				,
						'l-95% CI.partial'=sum.tmp$sol[2,2],
						'u-95% CI.partial'=sum.tmp$sol[2,3],
						t.partial=t.partial,
						pval.partial=pval.partial,
						z.partial=z.partial,
						df.partial=df.partial,
						#pMCMC
						pvalMCMC.partial=sum.tmp$sol[2,5],
						#z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
						#effective sample size
						n.effet.partial=sum.tmp$sol[2,4],
						n.MCMC.partial=sum.tmp$cstats[3],
						#Standard Deviation
						sd.partial=sd(M.partial$Sol[,2]),
						#only DIC saved, set AIC, BIC as DIC
						DIC.partial=sum.tmp$DIC,
						Dev.partial=mean(M.partial$Deviance)
				)
				rs<-c(rs,rs.partial)
			}

		}else if(family$family=='binomial'){

			if(family$link=='logit')
				prior.scale<-2.5
			else if(family$link=='probit')
				prior.scale<-2.5*1.6

			if('full'%in%pooling){
				M.full<-bayesglm(treatment ~ response, data=dat, family=family, n.iter = n.iter, prior.scale = prior.scale)
				sum.tmp<-summary(M.full)
				rs.full<-c(
						coef.full=sum.tmp$coef[2,1],
						se.full=sum.tmp$coef[2,2],
						z.full=sum.tmp$coef[2,3],
						pval.full=sum.tmp$coef[2,4],
						#a bug in bayeglm for summary, fix: total obs - rank of the model
						df.full=sum.tmp$df.residual-sum.tmp$df[1],
						AIC.full=sum.tmp$aic,
						BIC.full=AIC(M.full,k=log(nrow(dat))),
						Dev.full=sum.tmp$deviance
				)
				rs<-c(rs,rs.full)
			}

			if('no'%in%pooling){
				if(n.levels>1)
					M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
				else
					M.no<-bayesglm(treatment ~ response, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)

				sum.tmp<-summary(M.no)
				rs.no<-c(
						coef.no=sum.tmp$coef[2,1],
						se.no=sum.tmp$coef[2,2],
						pval.no=sum.tmp$coef[2,4],
						z.no=sum.tmp$coef[2,3],
						df.no=sum.tmp$df.residual-sum.tmp$df[1],
						AIC.no=sum.tmp$aic,
						BIC.no=AIC(M.no,k=log(nrow(dat))),
						Dev.no=sum.tmp$deviance
				)
				rs<-c(rs,rs.no)
			}
			if('partial'%in%pooling){
				#categorial=logit, ordinal=probit
				if(family$link=='logit'){
					glmm.family<-'categorial'
					if(n.levels>1){
						prior<-list(R = list(V = prior.V.scale, nu=n.levels), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
					}else{
						prior<-list(R = list(V = prior.V.scale, nu=n.levels))
					}
				}
				else if(family$link=='probit'){
					glmm.family<-'ordinal'
					if(n.levels>1){
						prior<-list(R = list(V = prior.V.scale, nu=n.levels+1), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
					}else{
						prior<-list(R = list(V = prior.V.scale, nu=n.levels+1))
					}
				}
				else
					stop('For multilevel model with Binomial family and Bayeisan method, only logit and probit model are supported!!! \n')

#			coef.scale<-ifelse(family$link=='logit', 2.5^2, (2.5*1.6)^2)
#			intercept.scale<-ifelse(family$link=='logit', 10^2, (10*1.6)^2)
#			scale<-var(dat$treatment)/var(dat$response)
				################################# prior to be modified #################################
				#alpha.mu=rep(0,2),alpha.V=diag(2)*25^2)))
				#########################################################################################
				sd.threshold<-prior.scale*4
				#coef.sign<- -sign(logFC)
#				while(sd.threshold>prior.scale*2){
				if(n.levels>1){
					M.partial<-MCMCglmm(treatment~response, random=~idh(1+response):probe,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
					df.partial<-nrow(dat)-(n.levels+1)*2
				}else{
					M.partial<-MCMCglmm(treatment~response,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
					df.partial<-nrow(dat)-2
				}
				sum.tmp<-summary(M.partial)
				t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
				pval.partial<-2*pt(-abs(t.partial),df=df.partial)
				z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
				rs.partial<-c(
						coef.partial=sum.tmp$sol[2,1],
						'l-95% CI.partial'=sum.tmp$sol[2,2],
						'u-95% CI.partial'=sum.tmp$sol[2,3],
						t.partial=t.partial,
						pval.partial=pval.partial,
						z.partial=z.partial,
						df.partial=df.partial,
						#pMCMC
						pvalMCMC.partial=sum.tmp$sol[2,5],
#						z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
						#effective sample size
						n.effet.partial=sum.tmp$sol[2,4],
						n.MCMC.partial=sum.tmp$cstats[3],
						#Standard Deviation
						sd.partial=sd(M.partial$Sol[,2]),
						#only DIC saved, set AIC, BIC as DIC
						DIC.partial=sum.tmp$DIC,
						Dev.partial=mean(M.partial$Deviance)
				)
				#coef.sign<-sign(rs.partial['coef.partial'])
				sd.threshold<-as.double(rs.partial['sd.partial'])
#				}
				rs<-c(rs,rs.partial)
			}
		}else if(family$family=='poisson'){
			if('full'%in%pooling){
				#comlete pooling
				M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter,family='poisson')
				sum.tmp<-summary(M.full)
				rs.full<-c(
						coef.full=sum.tmp$coef[2,1],
						se.full=sum.tmp$coef[2,2],
						t.full=sum.tmp$coef[2,3],
						pval.full=sum.tmp$coef[2,4],
						#z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
						df.full=sum.tmp$df.residual-sum.tmp$df[1],
						AIC.full=sum.tmp$aic,
						BIC.full=AIC(M.full,k=log(nrow(dat))),
						Dev.full=sum.tmp$deviance
				)
				rs<-c(rs,rs.full)
			}

			if('no'%in%pooling){
				#no pooling
				if(n.levels>1)
					M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter,family='poisson')
				else
					M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter,family='poisson')
				sum.tmp<-summary(M.no)
				rs.no<-c(
						coef.no=sum.tmp$coef[2,1],
						se.no=sum.tmp$coef[2,2],
						t.no=sum.tmp$coef[2,3],
						pval.no=sum.tmp$coef[2,4],
						#z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
						df.no=sum.tmp$df.residual-sum.tmp$df[1],
						AIC.no=sum.tmp$aic,
						BIC.no=AIC(M.no,k=log(nrow(dat))),
						Dev.no=sum.tmp$deviance
				)
				rs<-c(rs,rs.no)
			}

			if('partial'%in%pooling){
				#partial pooling with multilevel model of varing slopes and intercepts
				prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))

				if(n.levels>1){
					M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
					df.partial<-nrow(dat)-(n.levels+1)*2
				}
				else{
					prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu))
					M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
					df.partial<-nrow(dat)-2
				}
				sum.tmp<-summary(M.partial)
				t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
				pval.partial<-2*pt(-abs(t.partial),df=df.partial)
				z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
				rs.partial<-c(
						coef.partial=sum.tmp$sol[2,1],
						'l-95% CI.partial'=sum.tmp$sol[2,2],
						'u-95% CI.partial'=sum.tmp$sol[2,3],
						t.partial=t.partial,
						pval.partial=pval.partial,
						z.partial=z.partial,
						df.partial=df.partial,
						#pMCMC
						pvalMCMC.partial=sum.tmp$sol[2,5],
#					z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
						n.effet.partial=sum.tmp$sol[2,4],
						n.MCMC.partial=sum.tmp$cstats[3],
						#Standard Deviation
						sd.partial=sd(M.partial$Sol[,2]),
						#only DIC saved, set AIC, BIC as DIC
						DIC.partial=sum.tmp$DIC,
						Dev.partial=mean(M.partial$Deviance)
				)
				rs<-c(rs,rs.partial)
			}
		}else{
			stop('Only liner model and binomial family model are supported !!! \n')
		}
	}

	#taking care of extreme cases where se or sd is zero
	if(is.na(rs[grep('^t|^z',names(rs))]) & rs[grep('^se|^sd',names(rs))]==0){
		rs[grep('^t|^z',names(rs))]<-0
		rs[grep('^pval',names(rs))]<-1
	}

	rs

}




#combine Pvalues by Fisher or Stouffer's method
#input:
#dat, a dataframe containing pvalue.cols
#sign.cols: column labels indicating sign of stat, if NULL, use absolute values to calculate overall stat for Stouffer's method
#' @export
combinePvalues<-function(dat,pvalue.cols,sign.cols=NULL,FC.cols=NULL,logTransformed=FALSE,log.base=2,method=c('Stouffer','Fisher'),twosided=TRUE,signed=TRUE,byRow=FALSE){
	if(!is.data.frame(dat) & !is.matrix(dat))
		stop('dat must be in data.frame or matrix format!!! \n')

	dat.pvals<-as.matrix(dat[,pvalue.cols])

	if(missing(method))
		method<-'Stouffer'

	if(!is.null(sign.cols)){
		signs<-sign(as.matrix(dat[,sign.cols]))
		signs[signs==0]<-1
		if(ncol(signs) != ncol(dat.pvals)){
			stop('sign.cols must have the same length with pvalue.cols !!!')
		}
	}else{
		signs<-matrix(1,ncol=ncol(dat.pvals),nrow=nrow(dat.pvals))
	}
	if(byRow){
		rs0<-apply(dat.pvals*signs,1,combinePvalVector,method=method,twosided=twosided,signed=signed)
		rs<-t(rs0)
		colnames(rs)<-c('comStat','comPval')
		nPvals<-unlist(apply(dat.pvals,1,f<-function(v){
							v<-as.vector(v);sum(!(is.na(v)|is.null(v)))
						}))
		if(!is.null(FC.cols)){
			FCs<-as.data.frame(dat[,FC.cols])
			if(logTransformed){
				logFC<-apply(FCs,1,mean,na.rm = TRUE)
			}else{
				FCs[FCs==0]<-1
				logFC<-apply(sign(FCs)*log(abs(FCs),base=log.base),1,mean,na.rm = TRUE)
				logFC<-sign(logFC)*(log.base^abs(logFC))
			}
			rs<-data.frame(nPvals=nPvals,FC=logFC,rs)
		}else{
			rs<-data.frame(nPvals=nPvals,rs)
		}
	}else{
		rs0<-apply(dat.pvals*signs,2,combinePvalVector,method=method,twosided=twosided,signed=signed)
		rs<-as.vector(rs0)
		if(ncol(rs0)==1)
			names(rs)<-c('comStat','comPval')
		else
			names(rs)<-paste(rep(c('comStat','comPval'),ncol(rs0)),rep(colnames(rs0),each=2),sep='.')
		if(!is.null(FC.cols)){
#		FC<-apply(as.data.frame(dat[,FC.cols]),2,mean,na.rm = TRUE)
			FCs<-as.data.frame(dat[,FC.cols])
			if(logTransformed){
				logFC<-apply(FCs,2,mean,na.rm = TRUE)
			}else{
				FCs[FCs==0]<-1
				logFC<-apply(sign(FCs)*log(abs(FCs),base=log.base),2,mean,na.rm = TRUE)
				logFC<-sign(logFC)*(log.base^abs(logFC))
			}
			if(length(logFC)==1)
				names(logFC)<-'FC'
			rs<-c(nPvals=nrow(dat.pvals),logFC,rs)
		}else{
			rs<-c(nPvals=nrow(dat.pvals),rs)
		}
	}
	rs
}

#Debugging
#pvals<-c(0.21,-.05)
#pvals<-c(runif(3,0,0.5),runif(3,-0.15,0)) ;pvals
#combinePvalVector(pvals,method='Fisher',signed=T,twosided = T)
#combinePvalVector(pvals,signed=T, twosided = T)
#combinePvalVector(c(1,1),method='Fisher')
#pvals<-c(1,1)
#combine one pvalues vector, the signed version
# For Fisher's method: if twosided, pvalues are transformed into single
combinePvalVector<-function(pvals,method=c('Stouffer','Fisher'),signed=TRUE,twosided=TRUE){

	#remove NA pvalues
	pvals<-pvals[!is.na(pvals) & !is.null(pvals)]

	if(sum(is.na(pvals))>=1){
		stat<-NA
		pval<-NA
	}else{
		if(twosided & (sum(pvals>1 | pvals< -1)>=1))
			stop('pvalues must between 0 and 1!\n')
		if(!twosided & (sum(pvals>0.5 | pvals< -0.5)>=1))
			stop('One-sided pvalues must between 0 and 0.5!\n')

		if(missing(method))
			method<-'Stouffer'

		if(!signed){
			pvals<-abs(pvals)
		}

		signs<-sign(pvals)
		signs[signs==0]<-1

		if(grepl('Fisher',method,ignore.case = TRUE)){
			if(twosided & signed){
				neg.pvals<-pos.pvals<-abs(pvals)/2
				pos.pvals[signs<0]<-1-pos.pvals[signs<0]
				neg.pvals[signs>0]<-1-neg.pvals[signs>0]
			}else{
				neg.pvals<-pos.pvals<-abs(pvals)
			}

			pvals<-c(1,-1)*c(pchisq(-2*sum(log(as.numeric(pos.pvals))),df=2*length(pvals),lower.tail = FALSE)/2,pchisq(-2*sum(log(as.numeric(neg.pvals))),df=2*length(pvals),lower.tail = FALSE)/2)

			pval<-min(abs(pvals))[1]
			#if two pvals are equal, pick up the first one
			stat <- sign(pvals[abs(pvals)==pval])[1]*qnorm(pval,lower.tail=F)[1]
			pval<-2*pval
		}
		else if(grepl('Stou',method,ignore.case = TRUE)){
			if(twosided){
				zs<-signs*qnorm(abs(pvals)/2,lower.tail=FALSE)
				stat<-sum(zs)/sqrt(length(zs))
				pval<-2*pnorm(abs(stat),lower.tail=FALSE)
			}
			else{
				zs<-signs*qnorm(abs(pvals),lower.tail=FALSE)
				stat<-sum(zs)/sqrt(length(zs))
				pval<-pnorm(abs(stat),lower.tail=FALSE)
			}
		}
		else{
			stop('Only \"Fisher\" or \"Stouffer\" method is supported!!!\n')
		}
	}
	return(c(stat=stat,pvalue=pval))
}



#fold change function
#positive: class1/class0
#negative: class0/class1
FC <- function(x,cl,logTransformed=TRUE,log.base=2,average.method=c('geometric','arithmetic'),pseudoCount=0){
	x.class0 <- x[(cl == 0)]+pseudoCount
	x.class1 <- x[(cl == 1)]+pseudoCount
	if(missing(average.method))
		average.method<-'geometric'
	if(logTransformed){
		if(is.na(log.base)|log.base<0)
			stop('You must specify log.bsae !\n')
		logFC<-mean(x.class1)-mean(x.class0)
		FC.val<-sign(logFC)*log.base^abs(logFC)
	}else{
		logFC<-ifelse(average.method=='arithmetic',log(mean(x.class1))-log(mean(x.class0)),mean(log(x.class1)-mean(log(x.class0))))
		FC.val<-sign(logFC)*exp(abs(logFC))
	}
	FC.val[FC.val==0 | is.na(FC.val)]<-1
	FC.val
}



##scale normalization
#' @export
normalize.scale<-function(d,total=NULL,pseudoCount=1){

	if(!is.data.frame(d)) d<-data.frame(d)

	if(!all(d>0)) d<-d+pseudoCount
	s<-apply(d,2,sum)
	m<-ifelse(is.null(total),as.integer(mean(s)),as.integer(total))
	options(digits=2+nchar(m))
	fac<-m/s
	for(i in 1:length(s)){
		d[,i]<-round(d[,i]*fac[i],0)
	}
	if(!all(d>0)) d<-d+pseudoCount
	d

}


normalize.quantile<-function (M, ties = TRUE)
{
  n <- dim(M)
  if (is.null(n))
    return(M)
  if (n[2] == 1)
    return(M)
  O <- S <- array(, n)
  nobs <- rep(n[1], n[2])
  i <- (0:(n[1] - 1))/(n[1] - 1)
  for (j in 1:n[2]) {
    Si <- sort(M[, j], method = "quick", index.return = TRUE)
    nobsj <- length(Si$x)
    if (nobsj < n[1]) {
      nobs[j] <- nobsj
      isna <- is.na(M[, j])
      S[, j] <- approx((0:(nobsj - 1))/(nobsj - 1), Si$x,
                       i, ties = "ordered")$y
      O[!isna, j] <- ((1:n[1])[!isna])[Si$ix]
    }
    else {
      S[, j] <- Si$x
      O[, j] <- Si$ix
    }
  }
  m <- rowMeans(S)
  for (j in 1:n[2]) {
    if (ties)
      r <- rank(M[, j])
    if (nobs[j] < n[1]) {
      isna <- is.na(M[, j])
      if (ties)
        M[!isna, j] <- approx(i, m, (r[!isna] - 1)/(nobs[j] -
                                                      1), ties = "ordered")$y
      else M[O[!isna, j], j] <- approx(i, m, (0:(nobs[j] -
                                                   1))/(nobs[j] - 1), ties = "ordered")$y
    }
    else {
      if (ties)
        M[, j] <- approx(i, m, (r - 1)/(n[1] - 1), ties = "ordered")$y
      else M[O[, j], j] <- m
    }
  }
  M
}


###input.file: a tabular-separated file with the following columns: shRNA/sgRNA, gene, control_1, control_2, case_1, case_2, etc
#control.samples: column names of control samples
#case.samples: column names of case samples
#gene.columnId: column number for gene symbol, only 1 or 2, default is 2
#' @export
generateEset<-function(input.file,control.samples,case.samples,control.groupname='control',case.groupname='treatment',gene.columnId=2){

	require(Biobase)

	d<-read.table(input.file,sep='\t',header=T)

	gene.column<-as.integer(gene.columnId)
	if(!gene.column%in%1:2) stop ('gene.columnId must be 1 or 2!\n')
	id.column<-setdiff(1:2,gene.column)

	#id.column<-(1:2)[apply(apply(d[,1:2],2,duplicated),2,sum)==0]
	#gene.column<-setdiff(1:2,id.column)

	fd<-data.frame(d[,c(id.column,gene.column)],row.names=d[,id.column])
	names(fd)<-c('id','gene')

	if(is.character(c(control.samples,case.samples))){
		s.ctrl<-setdiff(control.samples,names(d))
		if(length(s.ctrl)>0) stop(paste('NO control samples: ',s.ctrl,'\n',collapse=', '))
		s.case<-setdiff(case.samples,names(d))
		if(length(s.case)>0) stop(paste('NO case samples: ',s.case,'\n',collapse=', '))
	}

	d1<-data.frame(d[,c(control.samples,case.samples)],row.names=d[,id.column])

	pd<-data.frame(sampleName=names(d1),group=c(rep(control.groupname,length(control.samples)),rep(case.groupname,length(case.samples))),condition=c(rep('control',length(control.samples)),rep('treatment',length(case.samples))),row.names=names(d1))

	eset<-new("ExpressionSet",phenoData = new("AnnotatedDataFrame",pd),featureData=new("AnnotatedDataFrame",fd),exprs=as.matrix(d1))

	eset

}



#diff representation analyis at shRNA level
#eset
#condition: treatment group names (must be in pData(eset)$condition)
#ctrl: ctrl grou names (paired with condition, must be in pData(eset)$condition)
#do.log2: whether to do log2 transformation for the data
#do.restand: whether do restand in the analysis
#family: model to use, gaussian, poisson
#estimation.method: Bayesian or MLE
DRAshRNALevel<-function(eset,condition,ctrl,ctrl.tag=NULL,do.log2=TRUE,do.restand=TRUE,filterLowCount=TRUE,filterBy='T0',count.cutoff=32,family=gaussian,estimation.method='Bayesian'){

  if(sum(names(fData(eset))=='gene')==0){
    if(sum(names(fData(eset))=='geneSymbol')==1){
      names(fData(eset))[names(fData(eset))=='geneSymbol']<-'gene'
    }else{
      stop('gene doesn\'t exist!\n')
    }
  }

  names(fData(eset))[1]<-'shId'
  DR<-fData(eset)

  if(length(condition)!=length(ctrl))
    stop('condion and ctrl have diff length!\n')

  if(!all(c(condition,ctrl)%in%pData(eset)$condition))
    stop('conditon, ctrl are not all in pData(eset)$condition!\n')

  eset.sel<-eset

  if(filterLowCount){
    sel<-grep(filterBy,pData(eset)$condition)
    if(length(sel)>1){
      eset.sel<-eset[apply(exprs(eset[,sel]),1,median)>=count.cutoff,]
    }
  }

  #i=1

  for(i in 1:length(condition)){

    cat(i,':',condition[i],'vs.', ctrl[i], 'in processs...\n')

    d.eset<-eset.sel[,pData(eset.sel)$condition %in% c(ctrl[i],condition[i])]

    if(do.log2){
      #take log2 for count data
      exprs(d.eset)<-log2(exprs(d.eset))
    }

    comp<-factor(gsub(condition[i],1,gsub(ctrl[i],0,pData(d.eset)$condition)))
    table(comp)
    d<-data.frame(shId=as.character(fData(d.eset)$shId),exprs(d.eset),stringsAsFactors=FALSE)

    #dr<-ddply(d,'shId','combRowEvid.2grps',comp=comp,family=poisson,method='Bayesian',n.iter=5000,nitt=25000,burnin=5000,thin=1,pooling=c('full'),logTransformed=FALSE,restand=FALSE,pseudoCount=0)

    dr<-plyr::ddply(d,'shId','combRowEvid.2grps',comp=comp,family=family,method=estimation.method,n.iter=5000,nitt=25000,burnin=5000,thin=1,pooling=c('full'),logTransformed=do.log2,restand=do.restand,pseudoCount=0)
    dr<-dr[,!names(dr)%in%c('n.levels')]

    if(do.log2){
      dr$AveSignal<-round(2^dr$AveSignal,0)
    }else{
      dr$AveSignal<-round(dr$AveSignal,0)
    }

    #FDR.full<-fdrtool(dr$pval.full,'pvalue',plot=FALSE, verbose=FALSE)$qval
    FDR.BH.full<-p.adjust(dr$pval.full,'BH')
    full.id<-grep('pval.full',names(dr))
    dr<-data.frame(dr[,1:full.id],FDR.BH.full=FDR.BH.full,dr[,(full.id+1):ncol(dr)])

    names(dr)<-gsub('.full','',names(dr))

    m<-min(dr$pval[dr$pval>0])[1]
    dr$pval[dr$pval==0]<-ifelse(m<1e-7,m,1e-7)
    dr$z<-sign(dr$t)*qnorm(dr$pval/2,lower.tail = FALSE)


    #median-based FC
    if(do.log2){
      FC.m<-data.frame(shId=as.character(fData(d.eset)$shId),medianFC=apply(data.frame(exprs(d.eset)[,comp==1]),1,median)-apply(data.frame(exprs(d.eset)[,comp==0]),1,median))
    }else{
      FC.m<-data.frame(shId=as.character(fData(d.eset)$shId),medianFC=log2(apply(data.frame(exprs(d.eset)[,comp==1]),1,median)/apply(data.frame(exprs(d.eset)[,comp==0]),1,median)))
    }

    FC.m$medianFC<-sign(FC.m$medianFC)*2^abs(FC.m$medianFC)
    FC.m$medianFC[FC.m$medianFC==0]<-1
    dr<-merge(dr,FC.m,by='shId')


    ctrl.short<-ifelse(is.null(ctrl.tag),gsub('.*(.*)\\.(.*)','\\2',ctrl[i]),ctrl.tag[i])
    names(dr)[-1]<-paste(names(dr)[-1],'.',sum(comp==1),'VS',sum(comp==0),'.',condition[i],'_VS_',ctrl.short,sep='')

    DR<-merge(DR,dr,by='shId',all.x=TRUE)

  }

  #add n.shRNA column
  n.shRNAs<-table(DR$gene)
  n.shRNAs<-data.frame(gene=names(n.shRNAs),n.shRNAs=as.integer(n.shRNAs))
  DR<-merge(DR,n.shRNAs,by='gene',all.x=TRUE)

  DR<-DR[,c(2,1,3,ncol(DR),4:(ncol(DR)-1))]

  condition<-gsub('z.','',grep('^z',names(DR),value=T));condition

  col.pre<-c(
    ###annotation columns
    names(DR)[1:(grep('^FC',names(DR))[1]-1)]
    ,
    #DR columns
    paste(rep(c('FC','pval','FDR.BH','z','AveSignal','medianFC'),each=length(condition)),rep(condition,6),sep='.')
  )


  fData(eset)<-data.frame(
    DR[
      as.integer(sapply(featureNames(eset),match,DR$shId))
      ,
      c(col.pre,
        setdiff(names(DR),col.pre)
      )
      ],row.names=featureNames(eset))

  eset

}
