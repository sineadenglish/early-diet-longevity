## meta-analysis of effect of early diet manipulation on lifespan across taxa
library(MCMCglmm); library(lme4)
source("http://www.math.mcmaster.ca/bolker/R/misc/coefplot_new.R")

# read in data

# setwd("...") -> set working directory to folder where files are located

datlong <- read.csv("DataMeanLongevity.csv")
datHR	<- read.csv("DataLogHR.csv")


# clarify certain data points and order factor levels
datlong$ExptLifeStage[datlong$ExptLifeStage=="PreNatal-PostNatal"] <- "Prenatal"  # Check analysis robust to this study being classed as post-natal (they do - see below)
datlong$ExptLifeStage <- droplevels(datlong$ExptLifeStage)

datHR$AdultDiet[datHR$AdultDiet=="O"] <- "C"
datHR$ExptLifeStage[datHR$ExptLifeStage=="PreNatal-PostNatal"] <- "Prenatal"  # Check analysis robust to this study being classed as post-natal (they do - see below)
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


### START ANALYSIS: 
# setwd("~/.../results_datlong/") => set working directory to folder in which will hold results for longevity analysis

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

# HPD interval for b0 (intercept) overlaps 0 for both effect-size models 


### do all again for logHR data

# setwd("~/.../results_datHR/") # change to new folder

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

# HPD interval for b0 (intercept) overlaps 0 => no publication bias




### CREATE PLOTS FOR PAPER:

library(metafor)

tiff(filename = "Fig1.tiff",
     width = 480, height = 480, units = "px", pointsize = 12)

par(mfrow=c(2,2))

# funnel plot - HR
res <- rma(yi=lnHR, vi=varlnHR, data=datHR) #, slab=StudyNo)
funnel(res, level=c(90, 95, 99), ylim=c(1,13),xlim=c(-1.5,1.5),shade=c("white", "gray", "darkgray"), refline=0,
	yaxis="seinv",xlab="effect size (ln(HR))",ylab="precision (1/SE)",digits=c(1,0))
#text(res$yi+runif(1,min=-0.1,max=-0.05),1/sqrt(res$vi),res$slab,cex=.8)
mtext("(a)", side=3, line=1, at=-2)

# coef plot - HR
coefplot(HRmod1.1, intercept=TRUE, col=c("white",rep("black",10)),
	xlim=c(-1.5,1), xlab=expression(paste("effect size (ln(HR))",sep="")),cex.var=1, top.axis=FALSE,main="",varnames=rev(vnames)) # varnames=" ", 
coefplot(HRmod0.1, intercept=TRUE, add=TRUE, col="grey20", pch.pts=15, cex.pts=1.6)
mtext("(b)", side=3, line=1, at=-2)

# funnel plot - datlong
res2 <- rma(yi=d, vi=vd, data=datlong) 
funnel(res2, level=c(90, 95, 99), ylim=c(1,13),xlim=c(-1.5,1.5),shade=c("white", "gray", "darkgray"), refline=0,
	yaxis="seinv",xlab=expression(paste("effect size (Hedge's ", italic(d),")",sep="")),ylab="precision (1/SE)",digits=c(1,0))
#text(res2$yi+runif(1,min=0.05,max=0.1),1/sqrt(res2$vi),res2$slab,cex=.8)
mtext("(c)", side=3, line=1, at=-2)

# coef plot - datlong
coefplot(dmod1.1, intercept=TRUE, col=c("white",rep("black",10)),
	xlim=c(-1.5,1), xlab=expression(paste("effect size (Hedge's ", italic(d),")",sep="")),cex.var=1, top.axis=FALSE, main="",varnames=rev(vnames)) # varnames=" ", 
coefplot(dmod0.1, intercept=TRUE, add=TRUE, col="grey20", pch.pts=15, cex.pts=1.6)
mtext("(d)", side=3, line=1, at=-2)

dev.off()



### SENSITIVITY ANALYSES: WORK OUT WHICH ARE INFLUENTIAL STUDIES IN LONGEVITY ANALYSIS

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


# get R2 with credible intervals (see Nakagawa & Schielzeth)

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
