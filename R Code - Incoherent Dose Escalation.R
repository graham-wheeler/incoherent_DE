###############################################
# Wheeler, GM (2018): "Incoherent dose-escalation in phase I trials using the escalation with overdose control approach",
# Statistical Papers, 59, 801-811 (2018). DOI: 10.1007/s00362-016-0790-7
# First published 24/06/2016
#
# R Code
#
# Graham Wheeler
###############################################

rm(list=ls())

####################
# Install packages #
####################

install.packages("R2OpenBUGS")
library(R2OpenBUGS)

#############
# FUnctions #
#############

ptox.fn<-function(rho0,gamma,dose.vec,theta){
  plogis((1/(gamma - dose.vec[1]))*(gamma*qlogis(rho0) - dose.vec[1]*qlogis(theta)+(qlogis(theta)-qlogis(rho0))*dose.vec))
  }

coherence.trial.fn<-function(start.seed, numsim, numpat, alpha, strong.prior=FALSE, sdose, ptox0, theta, start=1, doseskip=TRUE, dose.inc=NULL, permissible.alphas=seq(minalpha,maxalpha,by=0.01), n.chains=1, n.iter=10000, n.burnin=5000, n.thin=4){
  rounded.dosemat <- outmat <- array(NA, dim=c(numsim, numpat-1, length(permissible.alphas)))
  incomat<-matrix(NA, nrow=numsim, ncol = numpat-1)
  l <- length(sdose)
  Xmin <- min(sdose)
  Xmax <- max(sdose)
  x.mat <- y.mat <- tox.mat <- matrix(rep(NA, numsim*numpat), nrow = numsim)
  mtd.dists <- matrix(rep(NA, numsim*5), ncol = 5)
  seed.fail <- NA
  seed.count <- start.seed
  for(i in 1:numsim){
    # First patient treated at minimum dose (dose 1), with extremely low prob of tox (ideally)
    alpha.vec<-NULL
    y.vec<-1
    while(y.vec==1){
      x.vec<-Xmin
      N<-1
      next.vec<-start
      set.seed(seed.count)
      y.vec <- sum( runif(1) < ptox0[next.vec] )
      if(y.vec==1){
        seed.count<-seed.count+1
        seed.fail<-c(seed.fail, seed.count-1)
       }
      }
    while(N<=numpat){
      data1<-list(Xmin=Xmin, Xmax=Xmax, theta=theta, K=N, Y=c(y.vec,NA), X=c(x.vec,NA))
      if(strong.prior==FALSE){
        brugs<-bugs.fit(model=model1, data=data1, params=c("gamma"), n.chains=n.chains, n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin, program="openbugs", seed=NULL)
        }else{
        brugs<-bugs.fit(model=model1strong, data=data1, params=c("gamma"), n.chains=n.chains, n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin, program="openbugs", seed=NULL)
        }
      t<-unlist(lapply(1:n.chains,function(z) brugs[[z]][,1]))
      ewoc.dose<-quantile(t, alpha)
      # If all patients treated, return posterior distribution of gamma and do not compute future doses
      if(N==numpat){
        break
        }
      # Find dose "closest" to the "alpha"th quantile found in ewoc.dose
      next.dose<-which.min((ewoc.dose-sdose)^2)
      # Update record of recommended doses and records of DLTs/non-DLTs
      ifelse(doseskip==TRUE, x.vec<-c(x.vec,sdose[next.dose]),
             ifelse(sdose[next.dose] > x.vec[N], x.vec<-c(x.vec, sdose[which(sdose==x.vec[N])+dose.inc]),
                    x.vec<-c(x.vec, sdose[next.dose])))
      next.dose<-which(sdose==x.vec[N+1])
      y.j <- sum( runif(1) < ptox0[next.dose] )        
      y.vec<-c(y.vec, y.j)
      next.vec<-c(next.vec, next.dose)
      data1<-list(Xmin=Xmin, Xmax=Xmax, theta=theta, K=N+1, Y=c(y.vec[1:N],1,NA), X=c(x.vec,NA))
      if(strong.prior==FALSE){
        brugs<-bugs.fit(model=model1, data=data1, params=c("gamma"), n.chains=n.chains, n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin, program="openbugs", seed=NULL)
        }else{
        brugs<-bugs.fit(model=model1strong, data=data1, params=c("gamma"), n.chains=n.chains, n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin, program="openbugs", seed=NULL)
        }
      t<-unlist(lapply(1:n.chains,function(z) brugs[[z]][,1]))
      # array of Ntrials by Npatients by different feasibility bounds to choose next dose - used to determine incomat
      outmat[i,N,]<-round(quantile(t,permissible.alphas),0)
      rounded.dosemat[i,N,]<-sdose[sapply(1:length(outmat[i,N,]), function(z) which.min((outmat[i,N,z]-sdose)^2))]
      # matrix of Ntrials by Npatients - entries are smallest increment at which incoherence occurs - index of permissible.alphas
      incomat[i,N]<-ifelse(max(outmat[i,N,])<x.vec[N+1],NA,min(which(outmat[i,N,]>=x.vec[N+1])))
      # Add next patient and repeat
      N<-N+1
      }
    x.mat[i,]<-x.vec
    y.mat[i,]<-y.vec
    tox.mat[i,]<-ptox0[next.vec]
    mtd.dists[i,]<-c(mean(t),median(t),sd(t),quantile(t,0.025),quantile(t,0.975))
    cat(paste("Simulation",i,"complete\n"))
    seed.count<-seed.count+1
    }
  seed.fail<-seed.fail[-1]
  return(list(x.mat=x.mat, y.mat=y.mat, tox.mat=tox.mat, mtd.dists=mtd.dists, seed.fail=seed.fail, 
              outmat=outmat, incomat=incomat, rounded.dosemat=rounded.dosemat))
  }

scenalphas.fn<-function(scenario, alphas, ntrials, npats){
  obj<-get(paste("scen",scenario,"",sep=""))
  newinco<-matrix(NA,nrow=ntrials,ncol=npats-1)
  for(i in 1:ntrials){
    for(N in 1:(npats-1)){
      newinco[i,N]<-ifelse(max(obj$rounded.dosemat[i,N,])<obj$x.mat[i,N+1], NA, max(which(obj$rounded.dosemat[i,N,]<=obj$x.mat[i,N+1]))+1)
    }
  }
  matrix(alphas[t(newinco)],nrow=ntrials,byrow=T)
}

####################################
# Set up Simulations - Section 3.3 #
####################################

theta <- 1/3
dose.FU <- seq(140, 425, by = 1)
ptox <- ptox.fn(0.08, 300, dose.FU, theta)
num.scens <- 6
nsims <- 100
maxpat <- 40
start.seed <- 804
minalpha <- 0.25
maxalpha <- 0.50
alphas <- seq(minalpha, maxalpha, by = 0.01)
n.chains <- 2
n.iter <- 10000
n.burnin <- 5000
n.thin <- 2

###################
# Run Simulations #
###################

model1<-function(){
  for (i in 1:K){
    Y[i]~dbern(p[i])
    logit(p[i])<- (1/(gamma - Xmin))*(gamma*logit(rho0) - Xmin*logit(theta)+(logit(theta)-logit(rho0))*X[i])
    }
  gamma ~ dunif(Xmin, Xmax)
  rho0 ~ dunif(0,theta)
  }

# Original
scen1<-coherence.trial.fn(start.seed = start.seed, numsim = nsims, numpat = maxpat, alpha = minalpha, strong.prior=FALSE,
                          sdose=dose.FU, ptox0=ptox, theta=theta, start=1, doseskip=TRUE, permissible.alphas=alphas, 
                          n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin)

# 20 dose levels (every 15mg/m2)
scen2<-coherence.trial.fn(start.seed = start.seed, numsim = nsims, numpat = maxpat, alpha = minalpha, strong.prior=FALSE, 
                          sdose=dose.FU[seq(1,length(dose.FU),length=20)], ptox0=ptox[seq(1,length(dose.FU),length=20)],
                          theta=theta, start=1, doseskip=TRUE, permissible.alphas=alphas, n.chains=n.chains, 
                          n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin)

# 16 dose levels (every 19mg/m2)
scen3<-coherence.trial.fn(start.seed = start.seed, numsim = nsims, numpat = maxpat, alpha = minalpha, strong.prior=FALSE, 
                          sdose=dose.FU[seq(1,length(dose.FU),length=16)], ptox0=ptox[seq(1,length(dose.FU),length=16)],
                          theta=theta, start=1, doseskip=TRUE, permissible.alphas=alphas, 
                          n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin)

# 6 dose levels (every 57mg/m2)
scen4<-coherence.trial.fn(start.seed = start.seed, numsim = nsims, numpat = maxpat, alpha = minalpha, strong.prior=FALSE, 
                          sdose=dose.FU[seq(1,length(dose.FU),length=6)], ptox0=ptox[seq(1,length(dose.FU),length=6)],
                          theta=theta, start=1, doseskip=TRUE, permissible.alphas=alphas,
                          n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin)

# Strong skew prior

model1strong<-function(){
  for (i in 1:K){
    Y[i]~dbern(p[i])
    logit(p[i])<- (1/(gamma - Xmin))*(gamma*logit(rho0) - Xmin*logit(theta)+(logit(theta)-logit(rho0))*X[i])
    }
  gamma<-Xmin + (Xmax-Xmin)*mubeta
  mubeta ~ dbeta(3,7)
  rho0<-theta*rhomu
  rhomu~dbeta(7,3)
  }

scen5<-coherence.trial.fn(start.seed = start.seed, numsim = nsims, numpat = maxpat, alpha = minalpha, strong.prior=TRUE,
                          sdose=dose.FU, ptox0=ptox, theta=theta, start=1, doseskip=TRUE, permissible.alphas=alphas,
                          n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin)

# Strong centered prior

model1strong<-function(){
  for (i in 1:K){
    Y[i]~dbern(p[i])
    logit(p[i])<- (1/(gamma - Xmin))*(gamma*logit(rho0) - Xmin*logit(theta)+(logit(theta)-logit(rho0))*X[i])
    }
  gamma<-Xmin + (Xmax-Xmin)*mubeta
  mubeta ~ dbeta(5,5)
  rho0<-theta*rhomu
  rhomu~dbeta(5,5)
  }

scen6<-coherence.trial.fn(start.seed = start.seed, numsim = nsims, numpat = maxpat, alpha = minalpha, strong.prior=TRUE,
                          sdose=dose.FU, ptox0=ptox, theta=theta, start=1, doseskip=TRUE, permissible.alphas=alphas, 
                          n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin)

#################
# Plot Figure 1 #
#################

for(i in 1:num.scens){
  assign(paste0("scen",i,"alphas"), scenalphas.fn(i,nsims,npats=maxpat))
  assign(paste0("mean.inco",i), apply(get(paste0("scen",i,"alphas")),2,mean, na.rm=T))
  assign(paste0("med.inco",i), apply(get(paste0("scen",i,"alphas")),2,quantile, na.rm=T, probs=0.5))
  assign(paste0("quant",i,"low"), apply(get(paste0("scen",i,"alphas")),2,quantile, na.rm=T, probs=0.025))
  assign(paste0("quant",i,"hi"), apply(get(paste0("scen",i,"alphas")),2,quantile, na.rm=T, probs=0.975))
  assign(paste0("min",i), apply(get(paste0("scen",i,"alphas")),2,min, na.rm=T))
  assign(paste0("max",i), apply(get(paste0("scen",i,"alphas")),2,max, na.rm=T))
  }

cuts<-seq(10,maxpat,by=10)

par(plt=c(0.11,0.46,0.75,0.95))
plot(2:maxpat, mean.inco1,type="l",col="black",xlab="", ylab=expression(paste(alpha[n+1]^"min")), ylim=c(minalpha,maxalpha), main="Scenario 1", las=1, bty="L")
lines(2:maxpat, quant1low, lty=2, col="black")
lines(2:maxpat, quant1hi, lty=2, col="black")
min_feas_excess<-!is.na(scen1alphas)
min_feas_excess_n<-apply(min_feas_excess[,cuts - 1],2,sum)
for(i in 1:length(cuts)){
  mtext(paste0("(N = ",min_feas_excess_n[i],")"), 1, cex=0.75, line=2, at=cuts[i])  
 }
par(new=TRUE,plt=c(0.6,0.95,0.75,0.95))
plot(2:maxpat, mean.inco2, type="l", col="black", xlab="", ylab="", ylim=c(minalpha,maxalpha), main="Scenario 2", las=1, bty="L")
lines(2:maxpat, quant2low, lty=2, col="black")
lines(2:maxpat, quant2hi, lty=2, col="black")
min_feas_excess<-!is.na(scen2alphas)
min_feas_excess_n<-apply(min_feas_excess[,cuts - 1],2,sum)
for(i in 1:length(cuts)){
  mtext(paste0("(N = ",min_feas_excess_n[i],")"), 1, cex=0.75, line=2, at=cuts[i])  
 }
par(new=TRUE,plt=c(0.11,0.46,0.45,0.65))
plot(2:maxpat, mean.inco3, type="l", col="black",xlab="", ylab=expression(paste(alpha[n+1]^"min")), ylim=c(minalpha,maxalpha), main="Scenario 3", las=1, bty="L")
lines(2:maxpat, quant3low, lty=2, col="black")
lines(2:maxpat, quant3hi, lty=2, col="black")
min_feas_excess<-!is.na(scen3alphas)
min_feas_excess_n<-apply(min_feas_excess[,cuts - 1],2,sum)
for(i in 1:length(cuts)){
  mtext(paste0("(N = ",min_feas_excess_n[i],")"), 1, cex=0.75, line=2, at=cuts[i])  
 }
par(new=TRUE,plt=c(0.6,0.95,0.45,0.65))
plot(2:maxpat, mean.inco4, type="l", col="black",xlab="", ylab="", ylim=c(minalpha,maxalpha), main="Scenario 4", las=1, bty="L")
lines(2:maxpat, quant4low, lty=2, col="black")
lines(2:maxpat, quant4hi, lty=2, col="black")
min_feas_excess<-!is.na(scen4alphas)
min_feas_excess_n<-apply(min_feas_excess[,cuts - 1],2,sum)
for(i in 1:length(cuts)){
  mtext(paste0("(N = ",min_feas_excess_n[i],")"), 1, cex=0.75, line=2, at=cuts[i])  
 }
par(new=TRUE,plt=c(0.11,0.46,0.15,0.35))
plot(2:maxpat, mean.inco5, type="l", col="black", xlab="", ylab=expression(paste(alpha[n+1]^"min")), ylim=c(minalpha,maxalpha), main="Scenario 5", las=1, bty="L")
lines(2:maxpat, quant5low, lty=2, col="black")
lines(2:maxpat, quant5hi, lty=2, col="black")
min_feas_excess<-!is.na(scen5alphas)
min_feas_excess_n<-apply(min_feas_excess[,cuts - 1],2,sum)
for(i in 1:length(cuts)){
  mtext(paste0("(N = ",min_feas_excess_n[i],")"), 1, cex=0.75, line=2, at=cuts[i])  
 }
mtext("Patient", 1, cex=1, line=3.25, at=23)
par(new=TRUE,plt=c(0.6,0.95,0.15,0.35))
plot(2:maxpat, mean.inco6, type="l", col="black", xlab="", ylab="", ylim=c(minalpha,maxalpha), main="Scenario 6", las=1, bty="L")
lines(2:maxpat, quant6low, lty=2, col="black")
lines(2:maxpat, quant6hi, lty=2, col="black")
min_feas_excess<-!is.na(scen6alphas)
min_feas_excess_n<-apply(min_feas_excess[,cuts - 1],2,sum)
for(i in 1:length(cuts)){
  mtext(paste0("(N = ",min_feas_excess_n[i],")"), 1, cex=0.75, line=2, at=cuts[i])  
 }
mtext("Patient", 1, cex=1, line=3.25, at=23)
par(new=TRUE,plt=c(0.11,0.95,0.02,0.08))
plot(0,0,bty="n", xaxt="n",yaxt="n", type="n", xlim=c(-1,1),ylim=c(-1,1), xlab="", ylab="")
legend(-0.7,0.3,legend=c("Mean", "95% Credible Interval"), lty=c(1,2), horiz=TRUE)


###########################################################
# Code to produce results similar to those in for Table 1 #
###########################################################

############
# Function #
############

ewoc.seq.fn<-function(Xmin, Xmax, theta, N, y.vec, x.vec, minalpha, maxalpha, n.chains, n.iter, n.burnin, n.thin){
  ntrials<-length(y.vec)/maxpat
  outmat<-array(NA, dim = c(ntrials, N-1, length(seq(minalpha, maxalpha, by=0.01))))
  incomat<-matrix(NA, nrow=ntrials, ncol=N-1)
  for(i in 1:ntrials){
    for(j in 1:(N-1)){
      if(ntrials == 1){
        yvec_trial <- y.vec
        xvec_trial <- x.vec
      }else{
        yvec_trial <- y.vec[i,]
        xvec_trial <- x.vec[i,]
      }
      data1<-list(Xmin=Xmin, Xmax=Xmax, theta=theta, K=j+1, Y=c(yvec_trial[1:j],1,NA), X=c(xvec_trial[1:(j+1)],NA))
      brugs<-bugs.fit(model = model1, data=data1,  n.chains = n.chains, params=c("gamma"), n.burnin = n.burnin, n.iter = n.iter, 
                      n.thin = n.thin, program = "openbugs", seed = 1:n.chains)
      # array of Ntrials by Npatients by different feasibility bounds to choose next dose after "possibles" - used to determine incomat
      outmat[i,j,]<-round(quantile(unlist(lapply(1:n.chains,function(z) brugs[[z]][,1])),seq(minalpha,maxalpha,by=0.01)),0)
      # matrix of Ntrials by Npatients - entries are smallest increment at which incoherence occurs - index of seq(minalpha,maxalpha,by=0.01)
      incomat[i,j]<-ifelse(max(outmat[i,j,])<xvec_trial[j+1],NA,min(which(outmat[i,j,]>=xvec_trial[j+1])))
      }
    }
  return(list(incomat=incomat, outmat=outmat))
  }

############
# Run Code #
############

pick_trial<-7
# Choose "pick_trial" from, for example, previous Scenario 1 simulations...
foo<-ewoc.seq.fn(Xmin = min(dose.FU), Xmax = max(dose.FU), theta, N=maxpat, scen1$y.mat[pick_trial, ], 
                 scen1$x.mat[pick_trial, ], minalpha=minalpha, maxalpha=maxalpha, 
                 n.chains = n.chains, n.thin = n.thin, n.burnin = n.burnin, n.iter = n.iter)

final_out<-rbind(scen1$x.mat[pick_trial,], scen1$y.mat[pick_trial,], c(NA, alphas[foo$incomat[1, ]]))
row.names(final_out)<-c("Dose", "DLT Outcome", "alpha_min")
final_out

#######
# END #
#######