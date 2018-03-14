library('stringr')
library(BayesFactor)

root="../../share//hedgeData//"

readDat=function(fileString,topLen)
{
	cmd=paste("ls -1 ",root,fileString, " >temp",sep='')
	system(cmd)
	filenames=as.vector(read.table('temp'))

	dat=NULL
	nFiles=dim(filenames)[1]
	subLab=sesLab=1:nFiles
	for (n in 1:nFiles){
		string=filenames[n,1]
		sesLab[n]=as.integer(str_sub(string,start = -5, end= -5 ))
		subLab[n]=ifelse(str_length(string)==topLen,
			as.integer(str_sub(string,start=9,end=10)),
			as.integer(str_sub(string,start=9,end=9)))
		partial0=read.csv(head=F,
			paste(root,fileString,'//',as.character(string),sep=''))
		colnames(partial0)=
			c('blk','trl','stim','cond','acc','rt')
		sub=rep(subLab[n],length(partial0[,1]))
		session=rep(sesLab[n],length(partial0[,1]))
		partial1=cbind(sub,session,partial0)
		dat=rbind(dat,partial1)
	}
	return(dat)
}

cleanData=function(dat,rtBounds=c(.25,1.5),badSub=NULL,first=25,return=10,badCond=NULL,acc=T)
{
	bad1= dat$rt<rtBounds[1] | dat$rt>rtBounds[2]
	bad2 = dat$sub %in% badSub
	bad3 = dat$blk==1 & dat$trl < (first+1)
	bad4 = dat$blk>1 & dat$trl< (return+1)
	bad5 = dat$cond %in% badCond
	bad6=NULL
	if (acc) bad6 = !dat$acc
	bad= bad1 | bad2 | bad3 | bad4 |bad5 | bad6	
	return(dat[!bad,])
}

designEst=function(dat){
	dat$sub=as.integer(as.factor(dat$sub))	
	I=max(dat$sub)
	N=dim(dat)[1]
	condVal=c(-1/2,1/2)
	sesVal=c(-1,1)
	Xsubsession=matrix(0,nrow=N,ncol=2*I)
	Xsubsession[cbind(1:N,dat$sub+I*(dat$session-1))]=1
	Xcond=matrix(0,nrow=N,ncol=2)
	Xcond[cbind(1:N,dat$session)]=condVal[dat$cond]
	Xsubcond=matrix(0,nrow=N,ncol=I)
	Xsubcond[cbind(1:N,dat$sub)]=condVal[dat$cond]
	Xint=matrix(0,nrow=N,ncol=I)
	Xint[cbind(1:N,dat$sub)]=condVal[dat$cond]*sesVal[dat$session]
	return(cbind(Xsubsession,Xcond,Xsubcond,Xint))}

est=function(dat,rScale){
	dat$sub=as.integer(as.factor(dat$sub))
	I=max(dat$sub)
	X=designEst(dat)
	out=nWayAOV(dat$rt,X,rep(0:3,c(2*I,2,I,I)),rscale=rScale,posterior=T)
	return(out)}


est.effect=function(out){	
	top=dim(out)[2]
	I=(top-8)/4
	meanIndex=1+2*I+1:2
	mainIndex=3+2*I+(1:I)
	intIndex=3+3*I+(1:I)
	effect1=out[,meanIndex[1]]+out[,mainIndex]-out[,intIndex]
	effect2=out[,meanIndex[2]]+out[,mainIndex]+out[,intIndex]
	effect=cbind(colMeans(effect1),colMeans(effect2))
	return(effect)}	

est.cor=function(out){
	top=dim(out)[2]
	I=(top-8)/4
	meanIndex=1+2*I+1:2
	mainIndex=3+2*I+(1:I)
	intIndex=3+3*I+(1:I)
	effect1=out[,meanIndex[1]]+out[,mainIndex]-out[,intIndex]
	effect2=out[,meanIndex[2]]+out[,mainIndex]+out[,intIndex]
	effect=array(dim=c(dim(effect1),2))
	effect[,,1]=effect1
	effect[,,2]=effect2
	samp=apply(effect,1,cor)[2,]
	rho=(out[,top-1]-out[,top])/(out[,top-1]+out[,top])
	return(cbind(samp,rho))
}


design1Session = function(dat){ 
	I=max(dat$sub)
	condVal=c(-.5,.5)
	N=dim(dat)[1]
	Xcond=condVal[dat$cond]
	Xsubcond=matrix(0,nrow=N,ncol=I)
	Xsubcond[cbind(1:N,dat$sub)]=condVal[dat$cond]
	Xsub=matrix(0,nrow=N,ncol=I)
	Xsub[cbind(1:N,dat$sub)]=1
	X=cbind(Xsub,Xcond,Xsubcond)
	return(X)}
	
est1Session=function(dat,rScale){
	dat$sub=as.integer(as.factor(dat$sub))
	I=max(dat$sub)
	X=design1Session(dat)
	return(nWayAOV(dat$rt,X,rep(0:2,c(I,1,I)),rscale=rScale,posterior=T))}


designBF=function(dat){
	dat$sub=as.integer(as.factor(dat$sub))	
	I=max(dat$sub)
	N=dim(dat)[1]
	condVal=c(-1/2,1/2)
	sesVal=c(-1,1)
	Xsubsession=matrix(0,nrow=N,ncol=2*I)
	Xsubsession[cbind(1:N,dat$sub+I*(dat$session-1))]=1
	Xcond=matrix(0,nrow=N,ncol=2)
	Xcond[cbind(1:N,dat$session)]=condVal[dat$cond]
	Xsubcond=matrix(0,nrow=N,ncol=I)
	Xsubcond[cbind(1:N,dat$sub)]=condVal[dat$cond]
	Xint=matrix(0,nrow=N,ncol=I)
	Xint[cbind(1:N,dat$sub)]=condVal[dat$cond]*sesVal[dat$session]
	Xrand=matrix(0,nrow=N,ncol=2*I)
	Xrand[cbind(1:N,dat$sub+I*(dat$session-1))]=condVal[dat$cond]
	X0=cbind(Xsubsession)
	X1=cbind(Xsubsession,Xcond)
	X2=cbind(Xsubsession,Xcond,Xsubcond)
	X3=cbind(Xsubsession,Xcond,Xrand)
	X4=cbind(Xsubsession,Xcond,Xsubcond,Xint)
	return(list(X0=X0,X1=X1,X2=X2,X3=X3,X4=X4))
	}

bf=function(dat,rScale){
	b=1:5	
	names(b)=paste('Mod',0:4,sep='')
	dat$sub=as.integer(as.factor(dat$sub))
	dat$cond=1+dat$cond/2
	I=max(dat$sub)
	X=designBF(dat)
	b[1]=nWayAOV(dat$rt,X$X0,rep(0,2*I),rscale=rScale[1])$bf
	b[2]=nWayAOV(dat$rt,X$X1,rep(0:1,c(2*I,2)),rscale=rScale[1:2])$bf
	b[3]=nWayAOV(dat$rt,X$X2,rep(0:2,c(2*I,2,I)),rscale=rScale[c(1,2,3)])$bf
	b[4]=nWayAOV(dat$rt,X$X3,rep(0:2,c(2*I,2,2*I)),rscale=rScale[c(1,2,3)])$bf
	b[5]=nWayAOV(dat$rt,X$X4,rep(0:3,c(2*I,2,I,I)),rscale=rScale[c(1,2,4,5)])$bf
	return(b)}



