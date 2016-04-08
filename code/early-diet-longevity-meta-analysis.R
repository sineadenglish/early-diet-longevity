## meta-analysis of effect of early diet manipulation on lifespan across taxa
library(MCMCglmm); library(lme4)
source("http://www.math.mcmaster.ca/bolker/R/misc/coefplot_new.R")

# read in data

setwd("~//Dropbox//work//meta-analysis//DevAge//analysis//")

datlong <- read.csv("DataMeanLongevity.csv")
datHR	<- read.csv("DataLogHR.csv")


# clarify certain data points
datlong$ExptLifeStage[datlong$ExptLifeStage=="PreNatal-PostNatal"] <- "Prenatal"  # Check analysis robust to this being classed as post-natal (they do - see below)
datlong$ExptLifeStage <- droplevels(datlong$ExptLifeStage)

datHR$AdultDiet[datHR$AdultDiet=="O"] <- "C"
datHR$ExptLifeStage[datHR$ExptLifeStage=="PreNatal-PostNatal"] <- "Prenatal"  # Check analysis robust to this being classed as post-natal (they do - see below)
datHR$ExptLifeStage <- droplevels(datHR$ExptLifeStage)


datlong$Sex = factor(datlong$Sex,levels(datlong$Sex)[c(2,3,1)])
datHR$Sex = factor(datHR$Sex,levels(datHR$Sex)[c(2,3,1)])


##########################
## ADD EFFECT SIZES
##########################

################################
# funtion to calcuate d

d.func<-function(m1,m2,sd1,sd2,n1,n2){
	spool<-sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
	d<-(m1-m2)/spool
	return(d)
}

Vd.func<-function(n1,n2,d){
	Vd<-(n1+n2)/(n1*n2)+d^2/(2*(n1+n2-2))
	return(Vd)
}
###################################



# d1: (C-E); do offspring live shorter if receive poor nutrition in early life

D1<-d.func(datlong$MeanE,datlong$MeanC,datlong$SD_E,datlong$SD_C,datlong$NStartExpt,datlong$NStartControl)

VD1<-Vd.func(datlong$NStartExpt,datlong$NStartControl,D1)


datlong$d = D1
datlong$vd = VD1


## merge datasets and look at correlation between the effect sizes

datmerge <- merge(datlong[,c("EffectID","Phylum","Sex","ExptLifeStage","d","vd")],datHR[,c("EffectID","lnHR","varlnHR")],by="EffectID")
# not very tight correlation but seems sensible



## add studies to publication bias plot
usedat <- datlong #[!(dat$StudyNo %in% c(28,100)),]

wm1<-weighted.mean(usedat$d,1/usedat$vd,na.rm=T)

par(mfrow=c(1,1),bty="l")
# all data
plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (d)"),	
	ylab="precision (1/s.e.)",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=c(0,10),xlim=c(-2,2),main="All data")
	axis(1,at=seq(-2,2,1))
	axis(2,at=seq(0,10,5))

par(mfrow=c(2,3))

# examine data: adult diet, quality/quantity, sex, pre-natal or post-natal, invertebrate/vertebrate, catch-up

# 1 SEX

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (d)"),	
	ylab="precision (1/s.e.)",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=c(0,10),xlim=c(-2,2),main="Sex")
	axis(1,at=seq(-2,2,1))
	axis(2,at=seq(0,10,5))

plotdat <- usedat[usedat$Sex=="F",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$Sex=="M",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("M","F"),col=c("blue","red"),pch=16)

## no obvious effect

# 2 ADULT DIET

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (d)"),	
	ylab="precision (1/s.e.)",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=c(0,10),xlim=c(-2,2),main="Post-manipulation diet")
	axis(1,at=seq(-2,2,1))
	axis(2,at=seq(0,10,5))

plotdat <- usedat[usedat$AdultDiet=="R",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$AdultDiet=="C",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("Control","Deprived"),col=c("blue","red"),pch=16)

## no obvious effect

# 3 QUANTITY OR QUALITY

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (d)"),	
	ylab="precision (1/s.e.)",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=c(0,10),xlim=c(-2,2),main="Diet manipulation type")
	axis(1,at=seq(-2,2,1))
	axis(2,at=seq(0,10,5))

plotdat <- usedat[usedat$ManipType=="Quantity",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$ManipType=="Quality",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("Quality","Quantity"),col=c("blue","red"),pch=16)



## also mixed

# 4 PRENATAL OR POSTNATAL

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (d)"),	
	ylab="precision (1/s.e.)",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=c(0,10),xlim=c(-2,2),main="Pre or post natal")
	axis(1,at=seq(-2,2,1))
	axis(2,at=seq(0,10,5))

plotdat <- usedat[usedat$ExptLifeStage=="Prenatal",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$ExptLifeStage=="Postnatal",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("Postnatal","Prenatal"),col=c("blue","red"),pch=16)




# 5 VERTEBRATE OR INVERTEBRATE

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (d)"),	
	ylab="precision (1/s.e.)",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=c(0,10),xlim=c(-2,2),main="Vertebrate/Invertebrate")
	axis(1,at=seq(-2,2,1))
	axis(2,at=seq(0,10,5))

plotdat <- usedat[usedat$Phylum=="Vertebrate",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$Phylum=="Invertebrate",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("Invertebrate","Vertebrate"),col=c("blue","red"),pch=16)




# 6 CATCH UP GROWTH MEASURED

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (d)"),	
	ylab="precision (1/s.e.)",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=c(0,10),xlim=c(-2,2),main="Catch up")
	axis(1,at=seq(-2,2,1))
	axis(2,at=seq(0,10,5))

plotdat <- usedat[usedat$CatchUp==1,]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$CatchUp==0,]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)

legend("topright",legend=c("NoCatchup","Catchup"),col=c("blue","red"),pch=16)


# SAME AGAIN FOR LOG HR

## add studies to publication bias plot
usedat <- datHR #[!(dat$StudyNo %in% c(28,100)),]

usedat$d <- usedat$lnHR
usedat$vd <- usedat$varlnHR

wm1<-weighted.mean(usedat$d,1/usedat$vd,na.rm=T)
use.ylim <- c(0,12)
use.xlim <- c(-1.5,1.5)
par(mfrow=c(1,1),bty="l")
# all data
plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (lnHR)"),	
	ylab="precision (1/var(lnHR))",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=use.ylim,xlim=use.xlim,main="All data")
	axis(1,at=seq(use.xlim[1],use.xlim[2],0.5))
	axis(2,at=seq(use.ylim[1],use.ylim[2],2))

par(mfrow=c(2,3))

# examine data: adult diet, quality/quantity, sex, pre-natal or post-natal, invertebrate/vertebrate, catch-up

# 1 SEX

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (lnHR)"),	
	ylab="precision (1/var(lnHR))",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=use.ylim,xlim=use.xlim,main="Sex")
	axis(1,at=seq(use.xlim[1],use.xlim[2],0.5))
	axis(2,at=seq(use.ylim[1],use.ylim[2],2))

plotdat <- usedat[usedat$Sex=="F",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$Sex=="M",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("M","F"),col=c("blue","red"),pch=16)

## no obvious effect

# 2 ADULT DIET

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (lnHR)"),	
	ylab="precision (1/var(lnHR))",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=use.ylim,xlim=use.xlim,main="Post-manipulation diet")
	axis(1,at=seq(use.xlim[1],use.xlim[2],0.5))
	axis(2,at=seq(use.ylim[1],use.ylim[2],2))

plotdat <- usedat[usedat$AdultDiet=="R",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$AdultDiet=="C",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("Control","Deprived"),col=c("blue","red"),pch=16)

## no obvious effect

# 3 QUANTITY OR QUALITY

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (lnHR)"),	
	ylab="precision (1/var(lnHR))",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=use.ylim,xlim=use.xlim,main="Diet manipulation type")
	axis(1,at=seq(use.xlim[1],use.xlim[2],0.5))
	axis(2,at=seq(use.ylim[1],use.ylim[2],2))

plotdat <- usedat[usedat$ManipType=="Quantity",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$ManipType=="Quality",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("Quality","Quantity"),col=c("blue","red"),pch=16)



## also mixed

# 4 PRENATAL OR POSTNATAL

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (lnHR)"),	
	ylab="precision (1/var(lnHR))",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=use.ylim,xlim=use.xlim,main="Pre or post natal")
	axis(1,at=seq(use.xlim[1],use.xlim[2],0.5))
	axis(2,at=seq(use.ylim[1],use.ylim[2],2))

plotdat <- usedat[usedat$ExptLifeStage=="Prenatal",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$ExptLifeStage=="Postnatal",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("Postnatal","Prenatal"),col=c("blue","red"),pch=16)




# 5 VERTEBRATE OR INVERTEBRATE

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (lnHR)"),	
	ylab="precision (1/var(lnHR))",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=use.ylim,xlim=use.xlim,main="Vertebrate/Invertebrate")
	axis(1,at=seq(use.xlim[1],use.xlim[2],0.5))
	axis(2,at=seq(use.ylim[1],use.ylim[2],2))

plotdat <- usedat[usedat$Phylum=="Vertebrate",]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$Phylum=="Invertebrate",]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)
legend("topright",legend=c("Invertebrate","Vertebrate"),col=c("blue","red"),pch=16)



# 6 CATCH UP GROWTH MEASURED

plot(usedat$d,1/sqrt(usedat$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression("effect size (lnHR)"),	
	ylab="precision (1/var(lnHR))",cex.lab=1,cex.axis=1,bty="l",axes=FALSE,ylim=use.ylim,xlim=use.xlim,main="Catch up")
	axis(1,at=seq(use.xlim[1],use.xlim[2],0.5))
	axis(2,at=seq(use.ylim[1],use.ylim[2],2))

plotdat <- usedat[usedat$CatchUp==1,]
points(plotdat$d,1/sqrt(plotdat$vd),col="red",pch=16,cex=1.5)

plotdat <- usedat[usedat$CatchUp==0,]
points(plotdat$d,1/sqrt(plotdat$vd),col="blue",pch=16,cex=1.5)

legend("topright",legend=c("NoCatchup","Catchup"),col=c("blue","red"),pch=16)



### START ANALYSIS: 
setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datlong/")

dat = datlong

## DO THE ANALYSIS
prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# intercept-only model
mod0.1 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.2 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.3 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

# full model including all terms of interest						
mod1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.2 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.3 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

save(mod0.1,file="mod0.1.Rdata")
save(mod0.2,file="mod0.2.Rdata")
save(mod0.3,file="mod0.3.Rdata")

save(mod1.1,file="mod1.1.Rdata")
save(mod1.2,file="mod1.2.Rdata")
save(mod1.3,file="mod1.3.Rdata")

# if run, then load: 
load("mod0.1.Rdata"); load("mod0.2.Rdata"); load("mod0.3.Rdata")
load("mod1.1.Rdata"); load("mod1.2.Rdata"); load("mod1.3.Rdata")



### Model diagnostics: convergence, effects of parameters etc...

alltest <- c()

name=datlong

for(mod in c(0,1))
	{
		m1 <- get(paste("mod",mod,".1",sep=""))
		m2 <- get(paste("mod",mod,".2",sep=""))
		m3 <- get(paste("mod",mod,".3",sep=""))
		
		x1 <- gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
		x2 <- gelman.diag(list(m1$VCV[,c(1,3)],m2$VCV[,c(1,3)],m3$VCV[,c(1,3)]))
		x3 <- gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))
		
		allx <- c(unlist(x1),unlist(x2),unlist(x3))
		alltest <- rbind(alltest,allx)
		
	}

colnames(alltest)<-c("sol.est","sol.CI","vcv1.est","vcv2.est","vcv1.CI","vcv2.CI","multi.est","dev.est","dev.CI")
max(alltest) # all potential scale reduction factor less than 1.01


#### Output posterior modes and HPD intervals on each parameter in intercept and full models

x<-"--------------------------"

# fixed effects
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output1 <- posterior.mode(thismod$Sol)
	output2 <- HPDinterval(thismod$Sol[])
	output3 <- posterior.mode(thismod$VCV)
	print(paste(name,mod,".1",sep=""))
	print(output1)
	print(output2)
	print(output3)
	print(x)
}
		
# just HPD intervals for text
output <- c()
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output2 <- HPDinterval(thismod$Sol[])
	output2<-data.frame(output2)
	output2$modname<-paste(name,mod,".1",sep="")
	output <- rbind(output, output2)
}

write.csv(output, "HPDintervals.csv", quote=FALSE)


#########################
## FOR I^2
#########################
modnames <- c("mod")

modsummary <- c()
for(i in 1:length(modnames)) for(type in 0:1)
{
	thismod <- get(as.character(paste(modnames[i],type,".1",sep="")))
	thisdat	<- dat
	
	thisDIC <- thismod$DIC
	
	W <- 1/thisdat$vd
	
	s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) 
	
	I2<-100*(thismod$VCV[,"StudyNo"])/
	(thismod$VCV[,"StudyNo"]+thismod$VCV[,"units"]+s2)

	UnitEff<-100*(thismod$VCV[,"units"])/
	(thismod$VCV[,"StudyNo"]+thismod$VCV[,"units"]+s2)	
		
	temp <- cbind(as.character(paste(modnames[i],type,".1",sep="")),thisDIC,posterior.mode(I2),posterior.mode(UnitEff))
	
	modsummary <- rbind(modsummary, temp)

}

modsummary <- data.frame(modsummary) # ignore error message
colnames(modsummary) <- c("ModelName","DIC","I2","Resid")
write.csv(modsummary, "ModSummary.csv",quote=F,row.names=F)




#########################
## PUBLICATION BIAS
#########################
wm<-weighted.mean(dat$d,1/dat$vd,na.rm=T)


par(mfrow=c(1,1),bty="l",oma=c(2,5,1,1))


plot(dat$d, 1/sqrt(dat$vd), abline(v=c(0,wm),lty=c(2,1),col="black"),col="grey30",xlab="effect size (std. mean)",
	ylab="precision (1/s.e.)",cex.lab=1.5,cex.axis=1.2,bty="l",axes=FALSE,xlim=c(-2,2),ylim=c(0,20))
#	axis(1,at=seq(-40,40,10),cex.axis=1.3)
#	axis(2,at=seq(0,175,25),cex.axis=1.3)
	



# try: text(xpd=TRUE, min(x[,1]), 0.46, adj=0, "a", font=1, family="serif", cex=1)

#########################
## FOR EGGER PLOT
#########################

# re-run three full models using pr=TRUE

prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# standardized means
eg.mod<-MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=T,prior=prior1)
save(eg.mod, file="eg.mod.Rdata")

thismod <- eg.mod

thismod$Random$formula<-update(thismod$Random$formula, ~.+leg(mev, -1, FALSE):units)

#Pred<-predict(thismod, marginal=~leg(mev, -1, FALSE):units) ## not work... 

W<-1/dat$vd; Resp<-dat$d; thisdat<-dat; ZVar<-dat$vd;

Prec<-sqrt(W) # precision
ES<-Resp-predict(thismod, marginal=~leg(mev, -1, FALSE):units) # residuals
zES<-ES*Prec

prior1<-list(R=list(V=1E-10,nu=-1))

e0.1<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.1,file="E1.Rdata")

e0.2<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.2,file="E2.Rdata")

e0.3<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.3,file="E3.Rdata")

gelman.diag(list(e0.1$Sol[,1],e0.2$Sol[,1],e0.3$Sol[,1]))
gelman.diag(list(e0.1$VCV[,1],e0.2$VCV[,1],e0.3$VCV[,1])) # changed from DAH code (c(1,3))
gelman.diag(list(e0.1$Deviance,e0.2$Deviance,e0.3$Deviance))

# fixed effects
posterior.mode(e0.1$Sol[,1:2])
summary(e0.1$Sol[,1:2])
HPDinterval(e0.1$Sol[,1:2])

# random effects
posterior.mode(e0.1$VCV)
summary(e0.1$VCV)

# HPD interval for b0 (intercept) overlaps 0 for both effect-size models => no publication bias? 


### do all again for logHR data

setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datHR/")

# change to new folder

dat = datHR
dat$d = dat$lnHR
dat$vd = dat$varlnHR

## DO THE ANALYSIS
prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# intercept-only model
mod0.1 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.2 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.3 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

# full model including all terms of interest						
mod1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.2 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.3 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

save(mod0.1,file="mod0.1.Rdata")
save(mod0.2,file="mod0.2.Rdata")
save(mod0.3,file="mod0.3.Rdata")

save(mod1.1,file="mod1.1.Rdata")
save(mod1.2,file="mod1.2.Rdata")
save(mod1.3,file="mod1.3.Rdata")


# if run, then load: 
load("mod0.1.Rdata"); load("mod0.2.Rdata"); load("mod0.3.Rdata")
load("mod1.1.Rdata"); load("mod1.2.Rdata"); load("mod1.3.Rdata")



### Model diagnostics: convergence, effects of parameters etc...

alltest <- c()

#for(name in c("carry","match")) for(mod in c(2,3)) to check additional models
for(mod in c(0,1))
	{
		m1 <- get(paste("mod",mod,".1",sep=""))
		m2 <- get(paste("mod",mod,".2",sep=""))
		m3 <- get(paste("mod",mod,".3",sep=""))
		
		x1 <- gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
		x2 <- gelman.diag(list(m1$VCV[,c(1,3)],m2$VCV[,c(1,3)],m3$VCV[,c(1,3)]))
		x3 <- gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))
		
		allx <- c(unlist(x1),unlist(x2),unlist(x3))
		alltest <- rbind(alltest,allx)
		
	}

colnames(alltest)<-c("sol.est","sol.CI","vcv1.est","vcv2.est","vcv1.CI","vcv2.CI","multi.est","dev.est","dev.CI")
max(alltest) # all potential scale reduction factor less than 1.01


#### Output posterior modes and HPD intervals on each parameter in intercept and full models

x<-"--------------------------"

# fixed effects
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output1 <- posterior.mode(thismod$Sol)
	output2 <- HPDinterval(thismod$Sol[])
	output3 <- posterior.mode(thismod$VCV)
	print(paste(name,mod,".1",sep=""))
	print(output1)
	print(output2)
	print(output3)
	print(x)
}
		
# just HPD intervals for text
output <- c()
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output2 <- HPDinterval(thismod$Sol[])
	output2<-data.frame(output2)
	output2$modname<-paste(name,mod,".1",sep="")
	output <- rbind(output, output2)
}

write.csv(output, "HPDintervals.csv", quote=FALSE)


#########################
## FOR I^2
#########################
modnames <- c("mod")

modsummary <- c()
for(i in 1:length(modnames)) for(type in 0:1)
{
	thismod <- get(as.character(paste(modnames[i],type,".1",sep="")))
	thisdat	<- dat
	
	thisDIC <- thismod$DIC
	
	W <- 1/thisdat$vd
	
	s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) 
	
	I2<-100*(thismod$VCV[,"StudyNo"])/
	(thismod$VCV[,"StudyNo"]+thismod$VCV[,"units"]+s2)

	UnitEff<-100*(thismod$VCV[,"units"])/
	(thismod$VCV[,"StudyNo"]+thismod$VCV[,"units"]+s2)	
		
	temp <- cbind(as.character(paste(modnames[i],type,".1",sep="")),thisDIC,posterior.mode(I2),posterior.mode(UnitEff))
	
	modsummary <- rbind(modsummary, temp)

}

modsummary <- data.frame(modsummary) # ignore error message
colnames(modsummary) <- c("ModelName","DIC","I2","Resid")
write.csv(modsummary, "ModSummary.csv",quote=F,row.names=F)




#########################
## PUBLICATION BIAS
#########################
wm<-weighted.mean(dat$d,1/dat$vd,na.rm=T)


par(mfrow=c(1,1),bty="l",oma=c(2,5,1,1))


plot(dat$d, 1/sqrt(dat$vd), abline(v=c(0,wm),lty=c(2,1),col="black"),col="grey30",xlab="effect size (std. mean)",
	ylab="precision (1/s.e.)",cex.lab=1.5,cex.axis=1.2,bty="l",axes=FALSE,xlim=c(-2,2),ylim=c(0,20))
#	axis(1,at=seq(-40,40,10),cex.axis=1.3)
#	axis(2,at=seq(0,175,25),cex.axis=1.3)
	



# try: text(xpd=TRUE, min(x[,1]), 0.46, adj=0, "a", font=1, family="serif", cex=1)

#########################
## FOR EGGER PLOT
#########################

# re-run three full models using pr=TRUE

prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# standardized means
eg.mod<-MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=T,prior=prior1)
save(eg.mod, file="eg.mod.Rdata")

thismod <- eg.mod

thismod$Random$formula<-update(thismod$Random$formula, ~.+leg(mev, -1, FALSE):units)

#Pred<-predict(thismod, marginal=~leg(mev, -1, FALSE):units) ## not work... 

W<-1/dat$vd; Resp<-dat$d; thisdat<-dat; ZVar<-dat$vd;

Prec<-sqrt(W) # precision
ES<-Resp-predict(thismod, marginal=~leg(mev, -1, FALSE):units) # residuals
zES<-ES*Prec

prior1<-list(R=list(V=1E-10,nu=-1))

e0.1<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.1,file="E1.Rdata")

e0.2<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.2,file="E2.Rdata")

e0.3<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.3,file="E3.Rdata")

gelman.diag(list(e0.1$Sol[,1],e0.2$Sol[,1],e0.3$Sol[,1]))
gelman.diag(list(e0.1$VCV[,1],e0.2$VCV[,1],e0.3$VCV[,1])) # changed from DAH code (c(1,3))
gelman.diag(list(e0.1$Deviance,e0.2$Deviance,e0.3$Deviance))

# fixed effects
posterior.mode(e0.1$Sol[,1:2])
summary(e0.1$Sol[,1:2])
HPDinterval(e0.1$Sol[,1:2])

# random effects
posterior.mode(e0.1$VCV)
summary(e0.1$VCV)

# HPD interval for b0 (intercept) overlaps 0 => no publication bias? 





## CHECK ROBUSTNESS TO KEY STUDIES

# Pre- and post-natal driven by Ozanne study? 

outlier1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=datlong[datlong$StudyNo!=14,],
			family="gaussian",
			mev=datlong$vd[datlong$StudyNo!=14], 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)







##################
# Compare results when have Species as random term rather than Study
##################


setwd("~/Documents/work/meta-analysis/DevAge/analysis/results_datlong_species/")

dat = datlong

## DO THE ANALYSIS
prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# intercept-only model
mod0.1 <- MCMCglmm(d ~ 1,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.2 <- MCMCglmm(d ~ 1,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.3 <- MCMCglmm(d ~ 1,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

# full model including all terms of interest						
mod1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.2 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.3 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

save(mod0.1,file="mod0.1.Rdata")
save(mod0.2,file="mod0.2.Rdata")
save(mod0.3,file="mod0.3.Rdata")

save(mod1.1,file="mod1.1.Rdata")
save(mod1.2,file="mod1.2.Rdata")
save(mod1.3,file="mod1.3.Rdata")

# if run, then load: 
load("mod0.1.Rdata"); load("mod0.2.Rdata"); load("mod0.3.Rdata")
load("mod1.1.Rdata"); load("mod1.2.Rdata"); load("mod1.3.Rdata")



### Model diagnostics: convergence, effects of parameters etc...

alltest <- c()

name=datlong

for(mod in c(0,1))
	{
		m1 <- get(paste("mod",mod,".1",sep=""))
		m2 <- get(paste("mod",mod,".2",sep=""))
		m3 <- get(paste("mod",mod,".3",sep=""))
		
		x1 <- gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
		x2 <- gelman.diag(list(m1$VCV[,c(1,3)],m2$VCV[,c(1,3)],m3$VCV[,c(1,3)]))
		x3 <- gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))
		
		allx <- c(unlist(x1),unlist(x2),unlist(x3))
		alltest <- rbind(alltest,allx)
		
	}

colnames(alltest)<-c("sol.est","sol.CI","vcv1.est","vcv2.est","vcv1.CI","vcv2.CI","multi.est","dev.est","dev.CI")
max(alltest) # all potential scale reduction factor less than 1.01


#### Output posterior modes and HPD intervals on each parameter in intercept and full models

x<-"--------------------------"

# fixed effects
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output1 <- posterior.mode(thismod$Sol)
	output2 <- HPDinterval(thismod$Sol[])
	output3 <- posterior.mode(thismod$VCV)
	print(paste(name,mod,".1",sep=""))
	print(output1)
	print(output2)
	print(output3)
	print(x)
}
		
# just HPD intervals for text
output <- c()
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output2 <- HPDinterval(thismod$Sol[])
	output2<-data.frame(output2)
	output2$modname<-paste(name,mod,".1",sep="")
	output <- rbind(output, output2)
}

write.csv(output, "HPDintervals.csv", quote=FALSE)


#########################
## FOR I^2
#########################
modnames <- c("mod")

modsummary <- c()
for(i in 1:length(modnames)) for(type in 0:1)
{
	thismod <- get(as.character(paste(modnames[i],type,".1",sep="")))
	thisdat	<- dat
	
	thisDIC <- thismod$DIC
	
	W <- 1/thisdat$vd
	
	s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) 
	
	I2<-100*(thismod$VCV[,"Species"])/
	(thismod$VCV[,"Species"]+thismod$VCV[,"units"]+s2)

	UnitEff<-100*(thismod$VCV[,"units"])/
	(thismod$VCV[,"Species"]+thismod$VCV[,"units"]+s2)	
		
	temp <- cbind(as.character(paste(modnames[i],type,".1",sep="")),thisDIC,posterior.mode(I2),posterior.mode(UnitEff))
	
	modsummary <- rbind(modsummary, temp)

}

modsummary <- data.frame(modsummary) # ignore error message
colnames(modsummary) <- c("ModelName","DIC","I2","Resid")
write.csv(modsummary, "ModSummary.csv",quote=F,row.names=F)




#########################
## PUBLICATION BIAS
#########################
wm<-weighted.mean(dat$d,1/dat$vd,na.rm=T)


par(mfrow=c(1,1),bty="l",oma=c(2,5,1,1))


plot(dat$d, 1/sqrt(dat$vd), abline(v=c(0,wm),lty=c(2,1),col="black"),col="grey30",xlab="effect size (std. mean)",
	ylab="precision (1/s.e.)",cex.lab=1.5,cex.axis=1.2,bty="l",axes=FALSE,xlim=c(-2,2),ylim=c(0,20))
#	axis(1,at=seq(-40,40,10),cex.axis=1.3)
#	axis(2,at=seq(0,175,25),cex.axis=1.3)
	



# try: text(xpd=TRUE, min(x[,1]), 0.46, adj=0, "a", font=1, family="serif", cex=1)

#########################
## FOR EGGER PLOT
#########################

# re-run three full models using pr=TRUE

prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# standardized means
eg.mod<-MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=T,prior=prior1)
save(eg.mod, file="eg.mod.Rdata")

thismod <- eg.mod

thismod$Random$formula<-update(thismod$Random$formula, ~.+leg(mev, -1, FALSE):units)

#Pred<-predict(thismod, marginal=~leg(mev, -1, FALSE):units) ## not work... 

W<-1/dat$vd; Resp<-dat$d; thisdat<-dat; ZVar<-dat$vd;

Prec<-sqrt(W) # precision
ES<-Resp-predict(thismod, marginal=~leg(mev, -1, FALSE):units) # residuals
zES<-ES*Prec

prior1<-list(R=list(V=1E-10,nu=-1))

e0.1<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.1,file="E1.Rdata")

e0.2<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.2,file="E2.Rdata")

e0.3<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.3,file="E3.Rdata")

gelman.diag(list(e0.1$Sol[,1],e0.2$Sol[,1],e0.3$Sol[,1]))
gelman.diag(list(e0.1$VCV[,1],e0.2$VCV[,1],e0.3$VCV[,1])) # changed from DAH code (c(1,3))
gelman.diag(list(e0.1$Deviance,e0.2$Deviance,e0.3$Deviance))

# fixed effects
posterior.mode(e0.1$Sol[,1:2])
summary(e0.1$Sol[,1:2])
HPDinterval(e0.1$Sol[,1:2])

# random effects
posterior.mode(e0.1$VCV)
summary(e0.1$VCV)

# HPD interval for b0 (intercept) overlaps 0 for both effect-size models => no publication bias? 


### do all again for logHR data

setwd("~/Documents/work/meta-analysis/DevAge/analysis/results_datHR_species/")

# change to new folder

dat = datHR
dat$d = dat$lnHR
dat$vd = dat$varlnHR

## DO THE ANALYSIS
prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# intercept-only model
mod0.1 <- MCMCglmm(d ~ 1,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.2 <- MCMCglmm(d ~ 1,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.3 <- MCMCglmm(d ~ 1,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

# full model including all terms of interest						
mod1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.2 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.3 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

save(mod0.1,file="mod0.1.Rdata")
save(mod0.2,file="mod0.2.Rdata")
save(mod0.3,file="mod0.3.Rdata")

save(mod1.1,file="mod1.1.Rdata")
save(mod1.2,file="mod1.2.Rdata")
save(mod1.3,file="mod1.3.Rdata")


# if run, then load: 
load("mod0.1.Rdata"); load("mod0.2.Rdata"); load("mod0.3.Rdata")
load("mod1.1.Rdata"); load("mod1.2.Rdata"); load("mod1.3.Rdata")



### Model diagnostics: convergence, effects of parameters etc...

alltest <- c()

#for(name in c("carry","match")) for(mod in c(2,3)) to check additional models
for(mod in c(0,1))
	{
		m1 <- get(paste("mod",mod,".1",sep=""))
		m2 <- get(paste("mod",mod,".2",sep=""))
		m3 <- get(paste("mod",mod,".3",sep=""))
		
		x1 <- gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
		x2 <- gelman.diag(list(m1$VCV[,c(1,3)],m2$VCV[,c(1,3)],m3$VCV[,c(1,3)]))
		x3 <- gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))
		
		allx <- c(unlist(x1),unlist(x2),unlist(x3))
		alltest <- rbind(alltest,allx)
		
	}

colnames(alltest)<-c("sol.est","sol.CI","vcv1.est","vcv2.est","vcv1.CI","vcv2.CI","multi.est","dev.est","dev.CI")
max(alltest) # all potential scale reduction factor less than 1.01


#### Output posterior modes and HPD intervals on each parameter in intercept and full models

x<-"--------------------------"

# fixed effects
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output1 <- posterior.mode(thismod$Sol)
	output2 <- HPDinterval(thismod$Sol[])
	output3 <- posterior.mode(thismod$VCV)
	print(paste(name,mod,".1",sep=""))
	print(output1)
	print(output2)
	print(output3)
	print(x)
}
		
# just HPD intervals for text
output <- c()
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output2 <- HPDinterval(thismod$Sol[])
	output2<-data.frame(output2)
	output2$modname<-paste(name,mod,".1",sep="")
	output <- rbind(output, output2)
}

write.csv(output, "HPDintervals.csv", quote=FALSE)


#########################
## FOR I^2
#########################
modnames <- c("mod")

modsummary <- c()
for(i in 1:length(modnames)) for(type in 0:1)
{
	thismod <- get(as.character(paste(modnames[i],type,".1",sep="")))
	thisdat	<- dat
	
	thisDIC <- thismod$DIC
	
	W <- 1/thisdat$vd
	
	s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) 
	
	I2<-100*(thismod$VCV[,"Species"])/
	(thismod$VCV[,"Species"]+thismod$VCV[,"units"]+s2)

	UnitEff<-100*(thismod$VCV[,"units"])/
	(thismod$VCV[,"Species"]+thismod$VCV[,"units"]+s2)	
		
	temp <- cbind(as.character(paste(modnames[i],type,".1",sep="")),thisDIC,posterior.mode(I2),posterior.mode(UnitEff))
	
	modsummary <- rbind(modsummary, temp)

}

modsummary <- data.frame(modsummary) # ignore error message
colnames(modsummary) <- c("ModelName","DIC","I2","Resid")
write.csv(modsummary, "ModSummary.csv",quote=F,row.names=F)




#########################
## PUBLICATION BIAS
#########################
wm<-weighted.mean(dat$d,1/dat$vd,na.rm=T)


par(mfrow=c(1,1),bty="l",oma=c(2,5,1,1))


plot(dat$d, 1/sqrt(dat$vd), abline(v=c(0,wm),lty=c(2,1),col="black"),col="grey30",xlab="effect size (std. mean)",
	ylab="precision (1/s.e.)",cex.lab=1.5,cex.axis=1.2,bty="l",axes=FALSE,xlim=c(-2,2),ylim=c(0,20))
#	axis(1,at=seq(-40,40,10),cex.axis=1.3)
#	axis(2,at=seq(0,175,25),cex.axis=1.3)
	



# try: text(xpd=TRUE, min(x[,1]), 0.46, adj=0, "a", font=1, family="serif", cex=1)

#########################
## FOR EGGER PLOT
#########################

# re-run three full models using pr=TRUE

prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# standardized means
eg.mod<-MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~Species,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=T,prior=prior1)
save(eg.mod, file="eg.mod.Rdata")

thismod <- eg.mod

thismod$Random$formula<-update(thismod$Random$formula, ~.+leg(mev, -1, FALSE):units)

#Pred<-predict(thismod, marginal=~leg(mev, -1, FALSE):units) ## not work... 

W<-1/dat$vd; Resp<-dat$d; thisdat<-dat; ZVar<-dat$vd;

Prec<-sqrt(W) # precision
ES<-Resp-predict(thismod, marginal=~leg(mev, -1, FALSE):units) # residuals
zES<-ES*Prec

prior1<-list(R=list(V=1E-10,nu=-1))

e0.1<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.1,file="E1.Rdata")

e0.2<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.2,file="E2.Rdata")

e0.3<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.3,file="E3.Rdata")

gelman.diag(list(e0.1$Sol[,1],e0.2$Sol[,1],e0.3$Sol[,1]))
gelman.diag(list(e0.1$VCV[,1],e0.2$VCV[,1],e0.3$VCV[,1])) # changed from DAH code (c(1,3))
gelman.diag(list(e0.1$Deviance,e0.2$Deviance,e0.3$Deviance))

# fixed effects
posterior.mode(e0.1$Sol[,1:2])
summary(e0.1$Sol[,1:2])
HPDinterval(e0.1$Sol[,1:2])

# random effects
posterior.mode(e0.1$VCV)
summary(e0.1$VCV)

# HPD interval for b0 (intercept) overlaps 0 => no publication bias? 



#############################
####### check if behaves differently if study number considered a factor rather than integer

##
### START ANALYSIS: 
setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datlong_studyfactor/")

dat = datlong
dat$StudyNo = as.factor(dat$StudyNo)
class(dat$StudyNo)

## DO THE ANALYSIS
prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# intercept-only model
mod0.1 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.2 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
mod0.3 <- MCMCglmm(d ~ 1,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

# full model including all terms of interest						
mod1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.2 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
mod1.3 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

save(mod0.1,file="mod0.1.Rdata")
save(mod0.2,file="mod0.2.Rdata")
save(mod0.3,file="mod0.3.Rdata")

save(mod1.1,file="mod1.1.Rdata")
save(mod1.2,file="mod1.2.Rdata")
save(mod1.3,file="mod1.3.Rdata")

# if run, then load: 
load("mod0.1.Rdata"); load("mod0.2.Rdata"); load("mod0.3.Rdata")
load("mod1.1.Rdata"); load("mod1.2.Rdata"); load("mod1.3.Rdata")



### Model diagnostics: convergence, effects of parameters etc...

alltest <- c()

name=datlong

for(mod in c(0,1))
	{
		m1 <- get(paste("mod",mod,".1",sep=""))
		m2 <- get(paste("mod",mod,".2",sep=""))
		m3 <- get(paste("mod",mod,".3",sep=""))
		
		x1 <- gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
		x2 <- gelman.diag(list(m1$VCV[,c(1,3)],m2$VCV[,c(1,3)],m3$VCV[,c(1,3)]))
		x3 <- gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))
		
		allx <- c(unlist(x1),unlist(x2),unlist(x3))
		alltest <- rbind(alltest,allx)
		
	}

colnames(alltest)<-c("sol.est","sol.CI","vcv1.est","vcv2.est","vcv1.CI","vcv2.CI","multi.est","dev.est","dev.CI")
max(alltest) # all potential scale reduction factor less than 1.01


#### Output posterior modes and HPD intervals on each parameter in intercept and full models

x<-"--------------------------"

# fixed effects
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output1 <- posterior.mode(thismod$Sol)
	output2 <- HPDinterval(thismod$Sol[])
	output3 <- posterior.mode(thismod$VCV)
	print(paste(name,mod,".1",sep=""))
	print(output1)
	print(output2)
	print(output3)
	print(x)
}
		
# just HPD intervals for text
output <- c()
for(name in c("mod")) for (mod in 0:1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output2 <- HPDinterval(thismod$Sol[])
	output2<-data.frame(output2)
	output2$modname<-paste(name,mod,".1",sep="")
	output <- rbind(output, output2)
}

write.csv(output, "HPDintervals.csv", quote=FALSE)


#########################
## FOR I^2
#########################
modnames <- c("mod")

modsummary <- c()
for(i in 1:length(modnames)) for(type in 0:1)
{
	thismod <- get(as.character(paste(modnames[i],type,".1",sep="")))
	thisdat	<- dat
	
	thisDIC <- thismod$DIC
	
	W <- 1/thisdat$vd
	
	s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) 
	
	I2<-100*(thismod$VCV[,"StudyNo"])/
	(thismod$VCV[,"StudyNo"]+thismod$VCV[,"units"]+s2)

	UnitEff<-100*(thismod$VCV[,"units"])/
	(thismod$VCV[,"StudyNo"]+thismod$VCV[,"units"]+s2)	
		
	temp <- cbind(as.character(paste(modnames[i],type,".1",sep="")),thisDIC,posterior.mode(I2),posterior.mode(UnitEff))
	
	modsummary <- rbind(modsummary, temp)

}

modsummary <- data.frame(modsummary) # ignore error message
colnames(modsummary) <- c("ModelName","DIC","I2","Resid")
write.csv(modsummary, "ModSummary.csv",quote=F,row.names=F)




#########################
## PUBLICATION BIAS
#########################
wm<-weighted.mean(dat$d,1/dat$vd,na.rm=T)


par(mfrow=c(1,1),bty="l",oma=c(2,5,1,1))


plot(dat$d, 1/sqrt(dat$vd), abline(v=c(0,wm),lty=c(2,1),col="black"),col="grey30",xlab="effect size (std. mean)",
	ylab="precision (1/s.e.)",cex.lab=1.5,cex.axis=1.2,bty="l",axes=FALSE,xlim=c(-2,2),ylim=c(0,20))
#	axis(1,at=seq(-40,40,10),cex.axis=1.3)
#	axis(2,at=seq(0,175,25),cex.axis=1.3)
	



# try: text(xpd=TRUE, min(x[,1]), 0.46, adj=0, "a", font=1, family="serif", cex=1)

#########################
## FOR EGGER PLOT
#########################

# re-run three full models using pr=TRUE

prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# standardized means
eg.mod<-MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=T,prior=prior1)
save(eg.mod, file="eg.mod.Rdata")

thismod <- eg.mod

thismod$Random$formula<-update(thismod$Random$formula, ~.+leg(mev, -1, FALSE):units)

#Pred<-predict(thismod, marginal=~leg(mev, -1, FALSE):units) ## not work... 

W<-1/dat$vd; Resp<-dat$d; thisdat<-dat; ZVar<-dat$vd;

Prec<-sqrt(W) # precision
ES<-Resp-predict(thismod, marginal=~leg(mev, -1, FALSE):units) # residuals
zES<-ES*Prec

prior1<-list(R=list(V=1E-10,nu=-1))

e0.1<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.1,file="E1.Rdata")

e0.2<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.2,file="E2.Rdata")

e0.3<-MCMCglmm(zES~I(1/sqrt(ZVar)),family="gaussian",data=thisdat,verbose=F,nitt=50000,thin=25,burnin=25000)
save(e0.3,file="E3.Rdata")

gelman.diag(list(e0.1$Sol[,1],e0.2$Sol[,1],e0.3$Sol[,1]))
gelman.diag(list(e0.1$VCV[,1],e0.2$VCV[,1],e0.3$VCV[,1])) # changed from DAH code (c(1,3))
gelman.diag(list(e0.1$Deviance,e0.2$Deviance,e0.3$Deviance))

# fixed effects
posterior.mode(e0.1$Sol[,1:2])
summary(e0.1$Sol[,1:2])
HPDinterval(e0.1$Sol[,1:2])

# random effects
posterior.mode(e0.1$VCV)
summary(e0.1$VCV)

# HPD interval for b0 (intercept) overlaps 0 for both effect-size models => no publication bias? 

#### no difference, so will stick with original models for now. 



### SENSITIVITY ANALYSES: 

# 1. ONLY THOSE STUDIES WHICH APPEAR IN BOTH DATASETS

############# 
# work out discrepancy between different analysis

setwd("~/Documents/work/meta-analysis/DevAge/analysis/results_merged_datasets/")

# try run models on restricted dataset for both d and lnHR

datmerge <- merge(datlong[,c("EffectID","StudyNo","Phylum","Sex","ExptLifeStage","CatchUp","ManipType","AdultDiet","d","vd")],datHR[,c("EffectID","lnHR","varlnHR")],by="EffectID")

dat=datmerge

# full model including all terms of interest						
dmod1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
dmod1.2 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
dmod1.3 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

# full model including all terms of interest						
HRmod1.1 <- MCMCglmm(lnHR ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$varlnHR, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
HRmod1.2 <- MCMCglmm(lnHR ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$varlnHR, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

HRmod1.3 <- MCMCglmm(lnHR ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$varlnHR, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)




### Model diagnostics: convergence, effects of parameters etc...

alltest <- c()

for(name in c("d","HR"))
{
	for(mod in c(1))
	{
		m1 <- get(paste(name,"mod",mod,".1",sep=""))
		m2 <- get(paste(name,"mod",mod,".2",sep=""))
		m3 <- get(paste(name,"mod",mod,".3",sep=""))
		
		x1 <- gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
		x2 <- gelman.diag(list(m1$VCV[,c(1,3)],m2$VCV[,c(1,3)],m3$VCV[,c(1,3)]))
		x3 <- gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))
		
		allx <- c(unlist(x1),unlist(x2),unlist(x3))
		alltest <- rbind(alltest,allx)
		
	}

}

colnames(alltest)<-c("sol.est","sol.CI","vcv1.est","vcv2.est","vcv1.CI","vcv2.CI","multi.est","dev.est","dev.CI")
max(alltest) # all potential scale reduction factor less than 1.01

		
# just HPD intervals for text
output <- c()
for(name in c("dmod","HRmod")) for (mod in 1)
{
	thismod <- get(paste(name,mod,".1",sep=""))
	output2 <- HPDinterval(thismod$Sol[])
	output2<-data.frame(output2)
	output2$modname<-paste(name,mod,".1",sep="")
	output <- rbind(output, output2)
}

write.csv(output, "HPDintervals.csv", quote=FALSE)


#########################
## FOR I^2
#########################
modnames <- c("dmod","HRmod")

modsummary <- c()
for(i in 1:length(modnames)) for(type in 1)
{
	thismod <- get(as.character(paste(modnames[i],type,".1",sep="")))
	thisdat	<- dat
	
	thisDIC <- thismod$DIC
	
	W <- 1/thisdat$vd
	
	s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) 
	
	I2<-100*(thismod$VCV[,"StudyNo"])/
	(thismod$VCV[,"StudyNo"]+thismod$VCV[,"units"]+s2)

	UnitEff<-100*(thismod$VCV[,"units"])/
	(thismod$VCV[,"StudyNo"]+thismod$VCV[,"units"]+s2)	
		
	temp <- cbind(as.character(paste(modnames[i],type,".1",sep="")),thisDIC,posterior.mode(I2),posterior.mode(UnitEff))
	
	modsummary <- rbind(modsummary, temp)

}

modsummary <- data.frame(modsummary) # ignore error message
colnames(modsummary) <- c("ModelName","DIC","I2","Resid")
write.csv(modsummary, "ModSummary.csv",quote=F,row.names=F)




# 2. WORK OUT WHICH ARE INFLUENTIAL STUDIES IN LONGEVITY ANALYSIS AND REMOVE

setwd("~/Documents/work/meta-analysis/DevAge/analysis/results_datlong/")

sensOutput = data.frame(x=rep(NA,10))

dat = datlong
for(i in 2:length(unique(dat$StudyNo)))
{
	excludeStudy = unique(dat$StudyNo)[i]
	usedat = dat[dat$StudyNo!=excludeStudy,]
	
	sens_mod <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=usedat,
			family="gaussian",
			mev=usedat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)
			
	output2 <- HPDinterval(sens_mod$Sol[])
	output2<-data.frame(output2)
	output2 <- rbind(rep(excludeStudy,2),output2)
	
	sensOutput <- cbind(sensOutput, output2)
}		
		
write.csv(sensOutput,"sensitivity_datlong.csv",quote=F)



########## do the plot for the paper

# Panel plot (2,2) of funnel points and confidence intervals

# import models

setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datlong/")

load("mod0.1.Rdata"); load("mod1.1.Rdata")
dmod0.1 <- mod0.1
dmod1.1 <- mod1.1

setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datHR/")

load("mod0.1.Rdata"); load("mod1.1.Rdata")
HRmod0.1 <- mod0.1
HRmod1.1 <- mod1.1

vnames <- c("Adult diet (R)","Expt. Stage (prenatal)", "Manip. type (quantity)", "Sex (both)", "Sex (M)", "Catch-up (Y)", "Catch-up (N)", "Vert/invert (Vert)", "Intercept")
vnames_null <- rep("",9)

wm1<-weighted.mean(datlong$d,1/datlong$vd,na.rm=T)
wm2<-weighted.mean(datHR$lnHR,1/datHR$varlnHR,na.rm=T)


mtext("precision (1/s.e.)", side=2, outer=T)	
mtext("effect size",side=1,outer=T,line=-0.2,at=0.17)
mtext("(standardized mean)",side=1,outer=T,line=1.3,at=0.17)
	#,xlab="effect size (std. mean)"
	#ylab=""
	
plot(Data$D1,1/sqrt(Data$VD1), abline(v=c(0,wm2),lty=c(2,1),col="black"),col="grey30",xlab="",	
	ylab="",cex.lab=1.5,cex.axis=1.2,bty="l",axes=FALSE,ylim=c(0,30),xlim=c(-6,6))
	mtext("(b)",line=0, adj=-0.25)
	axis(1,at=seq(-6,6,2),cex.axis=1.3)
	axis(2,at=seq(0,30,5),cex.axis=1.3)
mtext("effect size",side=1,outer=T,line=-0.2,at=0.505)
mtext(expression(paste("(Hedge's ",italic(d)," - matching)")),side=1,outer=T,line=1.3,at=0.505)


quartz()
par(mfrow=c(2,2))

# funnel plot - HR
plot(datHR$lnHR, 1/sqrt(datHR$varlnHR), abline(v=c(0,wm2),lty=c(2,1),col="black"),col="grey30",xlab="effect size (ln(HR))",
	ylab="precision (1/SE)",cex.axis=1.2,bty="l",axes=FALSE,xlim=c(-1.5,1.5),ylim=c(0,15))
	axis(1,at=seq(-1.5,1.5,0.5),cex.axis=1)
	axis(2,at=seq(0,15,5),cex.axis=1)
mtext("(a)", side=3, line=1, at=-2)

# coef plot - HR
coefplot(HRmod1.1, intercept=TRUE, col=c("white",rep("black",10)),
	xlim=c(-1.5,1), xlab=expression(paste("effect size (ln(HR))",sep="")),cex.var=1, top.axis=FALSE,main="",varnames=rev(vnames)) # varnames=" ", 
coefplot(HRmod0.1, intercept=TRUE, add=TRUE, col="grey20", pch.pts=15, cex.pts=1.6)
mtext("(b)", side=3, line=1, at=-2)

# funnel plot - datlong
plot(datlong$d, 1/sqrt(datlong$vd), abline(v=c(0,wm1),lty=c(2,1),col="black"),col="grey30",xlab=expression(paste("effect size (Hedge's ", italic(d),")",sep="")),
	ylab="precision (1/SE)",cex.axis=1.2,bty="l",axes=FALSE,xlim=c(-1.5,1.5),ylim=c(0,15))
	axis(1,at=seq(-1.5,1.5,0.5),cex.axis=1)
	axis(2,at=seq(0,15,5),cex.axis=1)
mtext("(c)", side=3, line=1, at=-2)

# coef plot - datlong
coefplot(dmod1.1, intercept=TRUE, col=c("white",rep("black",10)),
	xlim=c(-1.5,1), xlab=expression(paste("effect size (Hedge's ", italic(d),")",sep="")),cex.var=1, top.axis=FALSE, main="",varnames=rev(vnames)) # varnames=" ", 
coefplot(dmod0.1, intercept=TRUE, add=TRUE, col="grey20", pch.pts=15, cex.pts=1.6)
mtext("(d)", side=3, line=1, at=-2)


## for appendix
par(mfrow=c(1,1))

plot(datmerge$lnHR,datmerge$d,xlab="ln(HR)",ylab=expression(paste("Hedge's ",italic(d))))



# for paper, create table of all studies

datall0 <- merge(datlong[,c("EffectID","StudyNo","Author","Year","Journal","Species","Phylum","d")],datHR[,c("EffectID","StudyNo","Author","Year","Journal","Species","Phylum","lnHR")],by=c("EffectID"),all.x=T,all.y=T)

datall <- merge(datlong[,c("EffectID","StudyNo","Author","Year","Journal","Species","Phylum","d")],datHR[,c("EffectID","StudyNo","Author","Year","Journal","Species","Phylum","lnHR")],by=c("EffectID","StudyNo","Author","Year","Journal","Species","Phylum"),all.x=T,all.y=T)


summdat<-ddply(datall,c("StudyNo","Author","Year","Journal","Species","Phylum"),summarize,d=length(na.omit(d)),lnHR=length(na.omit(lnHR)),n=length(StudyNo))

summdat <- summdat[with(summdat, order(Phylum,Year)),]




## try to do counter-shaded funnel plot
library(metafor)

par(mfrow=c(1,2))
res <- rma(yi=d, vi=vd, data=datlong, slab=StudyNo) # slab=paste(Author,Year,sep="-")) # slab=paste(Author,Year,sep="-")
funnel(res, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0,yaxis="seinv")
text(res$yi+runif(1,min=-0.1,max=-0.05),1/sqrt(res$vi),res$slab,cex=.8)

res2 <- rma(yi=lnHR, vi=varlnHR, data=datHR, slab=StudyNo)
funnel(res2, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0,yaxis="seinv")
text(res2$yi+runif(1,min=0.05,max=0.1),1/sqrt(res2$vi),res2$slab,cex=.8)


library(metafor)
 
### to save as png file
png(filename="contour_enhanced_funnel_plot.png",
    res=92, width=680, height=600, type="cairo")
 
### decrease margins so the full space is used
par(mar=c(5,4,1,2))
 
### load BCG vaccine data
data(dat.bcg)
 
### fit random-effects model
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
           slab=paste(author, year, sep=", "), method="REML")
 
### create contour enhanced funnel plot (with funnel centered at 0)
funnel(res, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)

dev.off()


head(dat.bcg)

?rma

## try with package 'meta' instead? 




###### calculate R2_mcmcglmm_marginal and R2_mcmcglmm_conditional for all models in table 

# import models

setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datlong/")

load("mod0.1.Rdata"); load("mod1.1.Rdata")
dmod0.1 <- mod0.1
dmod1.1 <- mod1.1

setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datHR/")

load("mod0.1.Rdata"); load("mod1.1.Rdata")
HRmod0.1 <- mod0.1
HRmod1.1 <- mod1.1

setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datlong_species/") # with species as random term

load("mod0.1.Rdata"); load("mod1.1.Rdata")
dmod0.1_sp <- mod0.1
dmod1.1_sp <- mod1.1

setwd("~/Dropbox/work/meta-analysis/DevAge/analysis/results_datHR_species/") # with species as random term

load("mod0.1.Rdata"); load("mod1.1.Rdata")
HRmod0.1_sp <- mod0.1
HRmod1.1_sp <- mod1.1


# get R2 with credible intervals (see Nakagawa & Schielzeth code)

mmF <- HRmod1.1_sp # cycle through for each model: dmod0.1, dmod1.1, HRmod0.1, HRmod1.1, dmod0.1_sp, dmod1.1_sp, HRmod0.1_sp, HRmod1.1_sp

vmVarF<-numeric(1000)
for(i in 1:1000){
Var<-var(as.vector(mmF$Sol[i,] %*% t(mmF$X)))
vmVarF[i]<-Var}

R2m<-vmVarF/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])
mean(R2m)
posterior.mode(R2m)
HPDinterval(R2m)

# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV[,1]+mmF$VCV[,2])/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])
mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c)


## check analyses similar if re-categorise the two ambiguous pre- versus post-natal studies: 

# read in data

setwd("~//Dropbox//work//meta-analysis//DevAge//analysis//")

datlong2 <- read.csv("DataMeanLongevity.csv")
datHR2	<- read.csv("DataLogHR.csv")

# clarify certain data points
datlong2$ExptLifeStage[datlong$ExptLifeStage=="PreNatal-PostNatal"] <- "Postnatal"  
datlong2$ExptLifeStage <- droplevels(datlong$ExptLifeStage)

datHR2$AdultDiet[datHR$AdultDiet=="O"] <- "C"
datHR2$ExptLifeStage[datHR$ExptLifeStage=="PreNatal-PostNatal"] <- "Postnatal"  
datHR2$ExptLifeStage <- droplevels(datHR$ExptLifeStage)

datlong2$Sex = factor(datlong2$Sex,levels(datlong2$Sex)[c(2,3,1)])
datHR2$Sex = factor(datHR2$Sex,levels(datHR2$Sex)[c(2,3,1)])

##########################
## ADD EFFECT SIZES
##########################


# d1: (C-E); do offspring live shorter if receive poor nutrition in early life

D1<-d.func(datlong2$MeanE,datlong2$MeanC,datlong2$SD_E, datlong2$SD_C,datlong2$NStartExpt,datlong2$NStartControl)

VD1<-Vd.func(datlong2$NStartExpt,datlong2$NStartControl,D1)


datlong2$d = D1
datlong2$vd = VD1


dat = datlong2
dat$StudyNo = as.factor(dat$StudyNo)
class(dat$StudyNo)

## DO THE ANALYSIS
prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# intercept-only model
mod1.1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

summary(mod1.1.1) # similar to dmod1.1


dat = datHR2
dat$d = dat$lnHR
dat$vd = dat$varlnHR

dat$StudyNo = as.factor(dat$StudyNo)
class(dat$StudyNo)

## DO THE ANALYSIS
prior1<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

# intercept-only model
mod2.1.1 <- MCMCglmm(d ~ Phylum + CatchUp + Sex + ManipType + ExptLifeStage + AdultDiet,
			random=~StudyNo,
			data=dat,
			family="gaussian",
			mev=dat$vd, 										
			nitt=400000,thin=100,burnin=300000,pr=F,prior=prior1,verbose=F)

summary(mod2.1.1) # similar to HRmod1.1







