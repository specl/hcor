require(R2jags)
require(coda)

# ==============================================================================
# 1. USER-FACING API FUNCTIONS
# ==============================================================================

hcorRun=function(dat,prior,chains=4,iter=1000,burnin=100){
  hcor.check(dat) # Automatically checks structure before crashing JAGS
  contrast=hcor.contrastCheck(dat)
  if(contrast){
    out=hcor.m2Run(dat,prior,chains,iter,burnin)
  } else {
    out=hcor.m1Run(dat,prior,chains,iter,burnin)
  }
  return(out)
}

hcorDiagnostic=function(chains){
  mcmc_fit <- as.mcmc(chains)
  mcmc_sel <- mcmc_fit[, c("rho", "sigma[1]", "sigma[2]", "pTau[1]", "pTau[2]"), drop = FALSE]
  return(gelman.diag(mcmc_sel, multivariate = TRUE))
}

hcorConv=function(dat){
  hcor.check(dat) # Automatically validates the dataset
  contrast=hcor.contrastCheck(dat)
  if (contrast){
    m0=tapply(dat$y,list(dat$sub,dat$task,dat$cond),mean)
    m=m0[,,2]-m0[,,1]
    reliability=hcor.reliability2(dat)
    spearman=hcor.spearman2(dat)
  } else {
    m=tapply(dat$y,list(dat$sub,dat$task),mean)
    reliability=hcor.reliability1(dat)
    spearman=hcor.spearman1(dat)
  }
  conventional=list(
    "cortest"=cor.test(m[,1],m[,2]),
    "reliability"=reliability,
    "spearman"=spearman)
  return(conventional)
}

hcorPlot=function(chains,dat){
  conv=hcorConv(dat)
  pars=chains$BUGSoutput$sims.list
  h=hist(pars$rho,breaks=seq(-1,1,.05),plot=F)
  top=max(h$density)
  hist(pars$rho,breaks=seq(-1,1,.05),prob=T,xlim=c(-1,1),ylim=c(0,top*1.1),   
       xlab="Correlation Coefficient",ylab="Density",
       border=NA,main="",col="lightblue")
  lines(c(-1,-1,1,1),c(0,.5,.5,.0))
  lines(conv$cortest$conf.int,rep(top*1.05,2),col='darkred')
  points(conv$cortest$estimate,top*1.05,pch=19,cex=1.3,col='darkred')
  points(conv$spearman,top*1.05,pch=8,col='darkred',cex=1.3)
  q=quantile(pars$rho,c(.025,.975))
  lines(q,rep(top*.5,2),col='darkblue')
  points(mean(pars$rho),top*.5,pch=19,cex=1.3,col='darkblue')
}

hcorRho=function(chains) chains$BUGSoutput$sims.list$rho

hcorBF=function(chains,interval){
  rho=hcorRho(chains)
  prior=.5*diff(interval)
  count=sum(interval[1]<rho & rho <=interval[2])
  post=(count+1)/(length(rho)+2)
  extremeCount = (count==length(rho) | count == 0)
  out=list(
    "interval"=interval,
    "bf"=post/prior,
    "extreme"=extremeCount,
    "count"=count,
    "N"=length(rho))
  return(out)
}

# ==============================================================================
# 2. DATA VALIDATION AND UTILITIES
# ==============================================================================

hcor.check=function(dat){
  if (!inherits(dat, "data.frame")) stop("Input data 'dat' must be a data frame.")
  required <- c("sub", "task", "y")
  for (cn in required){
    if(!(cn %in% names(dat))) stop(paste("Data frame is missing required column:", cn))
  }
  if (!all(dat$task %in% c(1, 2))) stop("The 'task' column must contain only 1s and 2s.")
}

hcor.contrastCheck=function(dat){
  has_cond <- 'cond' %in% names(dat)
  if (has_cond && !all(dat$cond %in% c(1, 2))) {
    stop("The 'cond' column must contain only 1s and 2s.")
  }
  return(has_cond)
}

# ==============================================================================
# 3. CLASSICAL RELIABILITY & ATTENUATION ESTIMATION
# ==============================================================================

hcor.reliability2=function(dat){
  I=max(dat$sub)
  m=tapply(dat$y,list(dat$sub,dat$task,dat$cond),mean)
  eff=m[,,2]-m[,,1]
  pred=m[cbind(dat$sub,dat$task,dat$cond)]
  sqrerr=(dat$y-pred)^2
  ss=tapply(sqrerr,dat$task,sum)
  N=table(dat$task)
  tauEst=ss/(N-2*I)
  L=tapply(dat$y,list(dat$sub,dat$task,dat$cond),length)
  p0=apply(1/L,2,sum)/I
  Vd=apply(eff,2,var)
  return(1-(tauEst*p0/Vd))
}

hcor.reliability1=function(dat){
  I=max(dat$sub)
  m=tapply(dat$y,list(dat$sub,dat$task),mean)
  pred=m[cbind(dat$sub,dat$task)]
  sqrerr=(dat$y-pred)^2
  ss=tapply(sqrerr,dat$task,sum)
  N=table(dat$task)
  tauEst=ss/(N-I)
  L=tapply(dat$y,list(dat$sub,dat$task),length)
  p0=apply(1/L,2,sum)/I
  Vd=apply(m,2,var)
  return(1-(tauEst*p0/Vd))
}

hcor.spearman2=function(dat){
  r=hcor.reliability2(dat)
  m=tapply(dat$y,list(dat$sub,dat$task,dat$cond),mean)
  eff=m[,,2]-m[,,1]
  obs=cor(eff)[1,2]
  return(obs/sqrt(prod(r)))
}

hcor.spearman1=function(dat){
  r=hcor.reliability1(dat)
  m=tapply(dat$y,list(dat$sub,dat$task),mean)
  obs=cor(m)[1,2]
  return(obs/sqrt(prod(r)))
}

# ==============================================================================
# 4. JAGS MODEL DEFINITIONS & INTERNAL RUNNERS
# ==============================================================================

# Model 1: One condition per task
hcor.m1Def = "
model{
  for (i in 1:I){
    for (j in 1:2){
    ybar[i,j] ~ dnorm(theta[i,j],L[i,j]*pTau[j])
    V[i,j] ~ dgamma((L[i,j]-1)/2,.5*(L[i,j]-1)*pTau[j])}}
  for (i in 1:I){  
    theta[i,1:2] ~ dmnorm(mu,Omega)}
  for (j in 1:2){
    mu[j] ~dnorm(mu.m[j],1/mu.sd[j]^2)
    sigma[j] ~ dunif(sig.lower[j],sig.upper[j])
    pTau[j] ~ dgamma(.5,.5*tau.scale[j]^2)}
  a=1/(1-rho^2)
  Omega[1,1] <- a/sigma[1]^2
  Omega[2,2] <- a/sigma[2]^2
  Omega[1,2] <- -(rho*a) /(sigma[1] * sigma[2])
  Omega[2,1] <- Omega[1,2]
  rho ~ dunif(-1, 1)  
}"

# Model 2: Two conditions per task (Experimental Contrast)
hcor.m2Def = "
model{
  for (i in 1:I){
    for (j in 1:2){
      for (k in 1:2){
        ybar[i,j,k] ~ dnorm(alpha[i,j]+x[k]*theta[i,j],L[i,j,k]*pTau[j])
        V[i,j,k] ~ dgamma((L[i,j,k]-1)/2,.5*(L[i,j,k]-1)*pTau[j])}}}
  for (i in 1:I){
    for (j in 1:2){ 
      alpha[i,j] ~ dnorm(alpha.m[j],1/alpha.sd[j]^2)}
    theta[i,1:2] ~ dmnorm(mu,Omega)}
  for (j in 1:2){
    mu[j] ~dnorm(mu.m[j],1/mu.sd[j]^2)
    sigma[j] ~ dunif(sig.lower[j],sig.upper[j])
    pTau[j] ~ dgamma(.5,.5*tau.scale[j]^2)}
  a=1/(1-rho^2)
  Omega[1,1] <- a/sigma[1]^2
  Omega[2,2] <- a/sigma[2]^2
  Omega[1,2] <- -(rho*a) /(sigma[1] * sigma[2])
  Omega[2,1] <- Omega[1,2]
  rho ~ dunif(-1, 1)  
}"

hcor.m1Run=function(dat,prior,chains,iter,burnin){
  L=table(dat$sub,dat$task)
  ybar=tapply(dat$y,list(dat$sub,dat$task),mean)
  V=tapply(dat$y,list(dat$sub,dat$task),var)
  setup=list(
    L=L,
    I=max(dat$sub),
    ybar=ybar,
    V=V)
  pars=c("theta","pTau","mu","sigma","rho")
  out=jags(data=c(setup,prior),parameters=pars, 
           model.file = textConnection(hcor.m1Def), 
           n.chains=chains,
           n.iter=iter,
           n.burnin=burnin)
  return(out)}

hcor.m2Run=function(dat,prior,chains,iter,burnin){
  L=table(dat$sub,dat$task,dat$cond)
  ybar=tapply(dat$y,list(dat$sub,dat$task,dat$cond),mean)
  V=tapply(dat$y,list(dat$sub,dat$task,dat$cond),var)
  setup=list(
    L=L,
    I=max(dat$sub),
    x=c(-.5,.5),
    ybar=ybar,
    V=V)
  pars=c("alpha","theta","pTau","mu","sigma","rho")
  out=jags(data=c(setup,prior),parameters=pars, 
           model.file = textConnection(hcor.m2Def), 
           n.chains=chains,
           n.iter=iter,
           n.burnin=burnin)
  return(out)}