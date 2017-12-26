library(truncnorm)
library(Iso)
library(splines)
library(locpol)
# x : exposure
# y : response, 0 or 1
# x is sorted in an ascending order and y is sorted accordingly
# isotonic regression estimate
pavaest=function(x,y){
  pxy=isoreg(x,y)$yf
  return(pxy)
}

# n : sample size
# ph: fitted prob.
# bw: bandwidth
# op is the quantity we are using to choose the optimal bandwidth
rnfun=function(y,n,ph,bw){
  sig2=mean(diff(y)^2)/2  
  op=(sum(y-ph)+sig2*1.5/bw)/n  
  return(op)
  
}
# xeval : where monotone local linear fitted probability is evaluated 
# the function returns the L and LI estimates
lifun=function(x,y,xeval,bw){
  
  ll=locLinSmootherC(x, y, xeval, bw, EpaK)$beta0  # local linear estimate
  
  ll[ll>1]=1  # truncation
  ll[ll<0]=0 
  
  count=c()
  for(i in 1:length(xeval)){
    count[i]=sum(x==xeval[i])
  }
  
  llp=rep(ll, count)
  
  lip=llp
  lirn=rn=rnfun(y,n,llp,bw)
  
  # if the ll estimate is not monotone, we apply isotonic regression to get the LI estimate
  
  if(any(diff(llp)<0)){
    lip=pavaest(x,llp)
    lip[lip>1]=1
    lip[lip<0]=0
    
    lirn=rnfun(y,n,lip,bw)
  }
  
  return(list(llp=llp,lip=lip,rn=rn,lirn=lirn))
}

# bws: grid points of bandwidth
# lws: length of bws

lep=function(x,y,xeval,n,bws,lws){
  
  lirn=rn=rep(1e6,lws) 
  
  for(i in 1:lws){
    
    obj=lifun(x,y,xeval,bws[i])
    rn[i]=obj$rn
    lirn[i]=obj$lirn
  }
  # choose the optimal bandwidth at which 'rn' and 'lirn' are minimized 
  
  orn=order(rn)[1]
  ob=bws[orn]
  lle=lifun(x,y,xeval,ob)$llp
  
  olp=order(lirn)[1]
  obli=bws[olp]
  lie=lifun(x,y,xeval,obli)$lip
  
  bh=c(ob,obli)
  
  return(list(lle=lle,lie=lie,bh=bh))
}

lifun2=function(x,y,xeval,bw){
  
  ll=locLinSmootherC(x, y, xeval, bw, EpaK)$beta0  # local linear estimate
  
  ll[ll>1]=1  # truncation
  ll[ll<0]=0      
  
  li=ll
  lirn=rn=rnfun(y,n,ll,bw)
  
  # if the ll estimate is not monotone, we apply isotonic regression to get the LI estimate
  
  if(any(diff(ll)<0)){
    li=pavaest(x,ll)
    li[li>1]=1
    li[li<0]=0
    
    lirn=rnfun(y,n,li,bw)
  }
  
  return(list(ll=ll,li=li,rn=rn,lirn=lirn))
}

# bws: grid points of bandwidth
# lws: length of bws

lep2=function(x,y,n,bws,lws){
  
  ll=li=matrix(0,n,lws)      
  lirn=rn=rep(1e6,lws) 
  
  for(i in 1:lws){
    
    obj=lifun2(x,y,x,bws[i])
    ll[,i]=obj$ll
    li[,i]=obj$li
    rn[i]=obj$rn
    lirn[i]=obj$lirn
  }
  # choose the optimal bandwidth at which 'rn' and 'lirn' are minimized 
  
  orn=order(rn)[1]
  ob=bws[orn]
  lle=ll[,orn]
  
  olp=order(lirn)[1]
  obli=bws[olp]
  lie=li[,olp]
  
  bh=c(ob,obli)
  
  return(list(lle=lle,lie=lie,bh=bh))
}

maf.fun=function(y,probx){
  probxx=(probx-probx[1])/probx
  maf=sum(probxx[y==1])/sum(y==1)
  return(maf)
}
## Qin and Leung estimator
# Evaluate semiparametric likelihood in Qin and Leung

ql.loglikelihood.for.optim=function(param,y,d){
  alpha=param[1];
  beta=param[2];
  lambdastar=param[3];
  p=param[4];
  m0=sum((d==0)[y==0]);
  n0=sum((d==0)[y==1]);
  N1=length(y)-m0-n0;
  m1=sum(y==0)-m0;
  n1=sum(y==1)-n0;
  # Take log transformation of positive values
  xpos=(d[d>0 & y==1]);
  tpos=(d[d>0]);
  feverpos=y[d>0];
  nu=lambdastar*(1/N1)*sum((exp(alpha+beta*xpos))/((1-lambdastar)+lambdastar*exp(alpha+beta*xpos)));
  l2=-sum(log(1+nu*(exp(alpha+beta*tpos)-1)))+sum(log((1-lambdastar)+lambdastar*exp(alpha+beta*xpos)));
  l1=m0*log(p)+m1*log(1-p)+n0*log(p*(1-lambdastar)/(1-lambdastar*p))+n1*log(1-p*(1-lambdastar)/(1-lambdastar*p));
  l1+l2;
}

# Function for taking in values and computing the ql estimator using
# optim
ql.optim=function(y,d){
  m0=sum((d==0)[y==0]);
  n0=sum((d==0)[y==1]);
  N1=length(y)-m0-n0;
  m1=sum(y==0)-m0;
  n1=sum(y==1)-n0;
  # Take log transformation of positive values
  xpos=(d[d>0 & y==1]);
  tpos=(d[d>0]);
  feverpos=y[d>0];
  # Starting values
  # Find starting value for p and lambda based on binomial likelihood
  parasite.prevalence=(d>0);
  lambdahat.binomial=(mean(parasite.prevalence[y==1])-mean(parasite.prevalence[y==0]))/(1-mean(parasite.prevalence[y==0]));
  phat.binomial=1-mean(parasite.prevalence[y==0]);
  p.curr=phat.binomial;
  lambdastar.curr=lambdahat.binomial/(1-phat.binomial*(1-lambdahat.binomial));
  # Assume normal distribution on logged positive values to estimate starting values for alpha, beta
  mu1=mean(tpos[feverpos==0]);
  sigmasq=var(tpos[feverpos==0]);
  mu2=(mean(tpos[feverpos==1])-(1-lambdastar.curr)*mu1)/lambdastar.curr;
  alpha.curr=(mu1^2-mu2^2)/(2*sigmasq);
  beta.curr=(mu2-mu1)/sigmasq;
  param=c(alpha.curr,beta.curr,lambdastar.curr,p.curr)
  
  temp.optim=optim(param,ql.loglikelihood.for.optim,method="L-BFGS-B",lower=c(-Inf,-Inf,.01,.01),upper=c(Inf,Inf,.99,.99),control=list(fnscale=-1),y=y,d=d);
  alphahat=temp.optim$par[1];
  betahat=temp.optim$par[2];
  lambdastarhat=temp.optim$par[3]
  phat=temp.optim$par[4];
  lambdahat=(lambdastarhat-lambdastarhat*phat)/(1-lambdastarhat*phat);
  list(lambdahat=lambdahat,alphahat=alphahat,betahat=betahat,lambdastarhat=lambdastarhat,phat=phat);
}

expit=function(x){
  val=exp(x)/(1+exp(x))
  val
}

####################################
simnum=10
maf.logi=maf.power=maf.nonpara=est.sp=rep(NA, simnum)

t1=Sys.time();
for(k in 1:simnum){
  n=500
  true.p=0.2
  #true.lamb=0.5
  beta=0.8
  
  p.mi=0.25
  p.nmi=0.2
  
  y.mi=c(rep(0, n*(1-p.mi)), rep(1, n*p.mi))
  y.nmi=c(rep(0, n*(1-p.nmi)*(1-p.mi)), rep(1, n*p.nmi*(1-p.mi)), rep(0, n*(1-p.nmi)*p.mi), rep(1, n*p.nmi*p.mi))
  y.obs=y.mi+y.nmi>0
  
  #rand.num=rbinom(n, size=1, prob=true.p)
  d.no.nmi=d.cur=d.obs=rep(NA, n)
  u=runif(n,0,1)
  for(i in 1:n){
    if(y.mi[i]==0){
      d.no.nmi[i]=ifelse(u[i]<true.p, 0, rtruncnorm(1, a=0, b=Inf, mean=5, sd=5))
      #d.no.nmi[i]=ifelse(u[i]<true.p, 0, rpois(1, lambda=5))
    }else{
      d.no.nmi[i]=rtruncnorm(1, a=0, b=Inf, mean=10, sd=5)
      #d.no.nmi[i]=rpois(1, lambda=10)
      #d.no.nmi[i]=rgamma(1, shape=5, rate=0.5)
      #d.no.nmi[i]=ifelse(rand.num2==0, rgamma(1, shape=6, rate=0.5), rgamma(1, shape=5, rate=0.5))
    }
    
    if(y.mi[i]==0 & y.nmi[i]==1){
      d.cur[i]=beta*d.no.nmi[i]
    }else{
      d.cur[i]=d.no.nmi[i]
    }
    d.obs[i]=rpois(1, lambda=d.cur[i])
    #d.obs[i]=d.cur[i]
  }
  

  ########################
  #### Semiparametric ####
  #### Qin and Leung  ####
  ########################
  
  phat=sum(y.obs[d.obs==0])/sum(d.obs==0)
  qlopt=try(ql.optim(y=y.obs, d=d.obs))
  if(isTRUE(class(qlopt)=="try-error")){ next} else{maf.semipara=qlopt$lambdahat}
  est.sp[k]=(maf.semipara - maf.semipara*mean(y.obs==1))/(1 - maf.semipara*mean(y.obs==1))

    
  ###
  ## Power
  
  powermodelest=function(y,d,lowerlimit.tau=.001,upperlimit.tau=5){
    
    devfunc=function(tau,y,d){
      d.tau=d^tau;
      tempmodel=glm(y~d.tau,family=binomial);
      deviance(tempmodel);
    }
    
    tauhat=optimize(devfunc,c(lowerlimit.tau,upperlimit.tau),y=y,d=d)$minimum;
    d.tauhat=d^tauhat;
    tempmodel=glm(y~d.tauhat,family=binomial);
    alphahat=coef(tempmodel)[1];
    betahat=coef(tempmodel)[2];
    list(alphahat=alphahat,betahat=betahat,tauhat=tauhat);
  }
  
  mypower=try(powermodelest(y=y.obs, d=d.obs))
  if(isTRUE(class(mypower)=="try-error")){ next} else{alphahat=as.numeric(mypower[1])
                                                      betahat=as.numeric(mypower[2])
                                                      tauhat=as.numeric(mypower[3])}
  
  
  sum.maf.func=function(alpha, beta, tau, d, y){
    prob.logi=(exp(alpha+beta*(d^tau)))/(1+exp(alpha+beta*(d^tau)))
    prob.logi[is.na(prob.logi)]=1
    p.nmi.e=exp(alpha)/(1+exp(alpha))
    maf=(prob.logi-p.nmi.e)/prob.logi
    sum.maf=sum(maf[y==1])
    final.maf=sum.maf/sum(y==1)
    return(final.maf)
  }
  
  maf.power[k]=sum.maf.func(alphahat, betahat, tauhat, d=d.obs, y=y.obs)
  
  
  ## Nonparametric
  log.dobs=log(1+d.obs)
  #log.dobs=d.obs
  order.log.dobs=log.dobs[order(d.obs)]
  order.y=y.obs[order(d.obs)]
  order.dobs=d.obs[order(d.obs)]
  bws=seq(1.5,3,by=0.1)
  lws=length(bws)
  
  log.leo=try(lep2(order.log.dobs, order.y, n, bws, lws))
  if(isTRUE(class(log.leo)=="try-error")){ next} else{maf.nonpara[k]=maf.fun(order.y, log.leo$lie)}
  

  #print(k)
  if(k%%10==0)cat("..", k)
}
t2=Sys.time();
t2-t1

maff.mat=cbind(maf.power, est.sp, maf.nonpara)
apply(maff.mat, 2, mean)

#save(maff.mat, file="res.RData")


