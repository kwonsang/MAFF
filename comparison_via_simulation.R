library(splines)


g.alpha1=function(alpha, mat){
  alpha0=exp(alpha[1])/(1+exp(alpha[1]))
  new.vec=mat%*%alpha[-1]
  phi.alpha=log(sum(exp(new.vec)))
  func.g=exp(new.vec-phi.alpha)
  func.g=as.vector(func.g)
  return(c(alpha0, (1-alpha0)*func.g))
}
g.alpha2=function(alpha,eta,d,mat){
  new.d=(d-mean(d))/sd(d)
  #new.d=d
  new.vec=mat%*%alpha+eta*new.d
  phi.alpha=log(sum(exp(new.vec)))
  func.g=exp(new.vec-phi.alpha)
  func.g=as.vector(func.g)
  return(c(0,func.g))
}

pois.err=function(x, d){
  n=length(x)
  m=length(d)
  mat=matrix(NA, nrow=n, ncol=m)
  for(i in 1:n){
    #mat[i,]=dpois(x[i], lambda=(exp(d)-1))
    mat[i,]=dpois(x[i], lambda=d)
  }
  return(mat)
}


freq.count=function(data, x, d){
  dobs=data[,1]; yobs=data[,2]
  n=length(x)
  freq0=freq1=rep(NA, n)
  for(i in 1:(n-1)){
    freq0[i]=sum(dobs[yobs==0]>=x[i] & dobs[yobs==0]<x[i+1])
    freq1[i]=sum(dobs[yobs==1]>=x[i] & dobs[yobs==1]<x[i+1])
  }
  freq0[n]=sum(dobs[yobs==0]>=x[n])
  freq1[n]=sum(dobs[yobs==1]>=x[n])
  freq.mat=cbind(freq0, freq1)
  return(freq.mat)
}

penal.poisson.fke=function(par,x,d,freq.mat,mat, fke, c0){
  pois.err.mat=pois.err(x=x, d=d)
  pois.err.mat.fke=pois.err(x=x, d=fke*d)
  
  m=length(par)
  lamb=exp(par[1])/(1+exp(par[1]))
  alpha=par[2:(m-1)]
  eta=par[m]
  
  g1=g.alpha1(alpha=alpha, mat=mat)
  g2=g.alpha2(alpha=alpha[-1], eta=eta, d=d[-1], mat=mat)
  trans.g1=pois.err.mat%*%g1
  trans.g2=pois.err.mat%*%g2
  trans.g1star=pois.err.mat.fke%*%g1
  
  freq0=freq.mat[,1]; freq1=freq.mat[,2]
  likeli.y0=sum(freq0*log(trans.g1))
  likeli.y1=sum(freq1*log(trans.g1star*(1-lamb)+trans.g2*lamb))
  likeli=likeli.y0+likeli.y1
  return(-likeli+c0*sqrt(sum(alpha[-1]^2)))
}




#####################################################
simnum=1
est.maff=penal.est.maff=rep(NA, simnum)

t1=Sys.time();
for(k in 1:simnum){
  library(truncnorm)
  n=500
  q=0.2
  true.p=q
  #true.lamb=0.5
  beta=1
  
  p.mi=0.25
  p.nmi=0.2
  
  y.mi=c(rep(0, n*(1-p.mi)), rep(1, n*p.mi))
  y.nmi=c(rep(0, n*(1-p.nmi)*(1-p.mi)), rep(1, n*p.nmi*(1-p.mi)), rep(0, n*(1-p.nmi)*p.mi), rep(1, n*p.nmi*p.mi))
  y.obs=y.mi+y.nmi>0
  
  d.no.nmi=d.cur=d.obs=rep(NA, n)
  u=runif(n,0,1)
  for(i in 1:n){
    if(y.mi[i]==0){
      d.no.nmi[i]=ifelse(u[i]<true.p, 0, rtruncnorm(1, a=0, b=Inf, mean=5, sd=5))
      #d.no.nmi[i]=ifelse(u[i]<true.p, 0, rexp(1, rate=1/100))
    }else{
      d.no.nmi[i]=rtruncnorm(1, a=0, b=Inf, mean=10, sd=5)
      #d.no.nmi[i]=rexp(1, rate=1/200)
    }
    
    if(y.mi[i]==0 & y.nmi[i]==1){
      d.cur[i]=beta*d.no.nmi[i]
    }else{
      d.cur[i]=d.no.nmi[i]
    }
    #d.obs[i]=rnbinom(1, mu=d.cur[i], size=10)
    d.obs[i]=rpois(1, lambda=d.cur[i])
  }
  ####
  x.grid=sort(unique(d.obs))
  theta.space=seq(0,20,by=0.5)
  
  d.grid=theta.space
  
  freq.dobs=freq.count(data=cbind(d.obs, y.obs), x=x.grid, d=d.grid)
  
  df=5
  q.mat=ns(d.grid[-1], df=df-1)
  for(j in 1:length(q.mat[1,])){
    q.mat[,j]=(q.mat[,j]-mean(q.mat[,j]))/sd(q.mat[,j])
  }
  
  ###
  fke=beta
  penal.est.par=nlminb(start=rep(0,df+2), penal.poisson.fke, lower=rep(-100,df+2), upper=rep(100,df+2), x=x.grid, d=d.grid, freq.mat=freq.dobs, mat=q.mat, fke=fke, c0=50)
  est.par=nlminb(start=rep(0,df+2), penal.poisson.fke, lower=rep(-100,df+2), upper=rep(100,df+2), x=x.grid, d=d.grid, freq.mat=freq.dobs, mat=q.mat, fke=fke, c0=0)
  
  penal.lambdastar.est=exp(penal.est.par$par[1])/(1+exp(penal.est.par$par[1]))
  penal.lambda.est=penal.lambdastar.est*(1-mean(y.obs==1))/(1-mean(y.obs==1)*penal.lambdastar.est)
  lambdastar.est=exp(est.par$par[1])/(1+exp(est.par$par[1]))
  lambda.est=lambdastar.est*(1-mean(y.obs==1))/(1-mean(y.obs==1)*lambdastar.est)
  est.maff[k]=lambda.est
  penal.est.maff[k]=penal.lambda.est
  print(k)
  #if(k%%10==0)cat("..", k)
}

mean(est.maff); sd(est.maff)
mean(penal.est.maff); sd(penal.est.maff)

