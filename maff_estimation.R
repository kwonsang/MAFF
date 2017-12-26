library(splines)


g.alpha1=function(alpha, mat){
  alpha0=alpha[1]
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
  freq0[n]=sum(dobs[yobs==0]>x[n])
  freq1[n]=sum(dobs[yobs==1]>x[n])
  freq.mat=cbind(freq0, freq1)
  return(freq.mat)
}

penal.poisson.fke=function(par,x,d,freq.mat,mat, fke, c0){
  pois.err.mat=pois.err(x=x, d=d)
  pois.err.mat.fke=pois.err(x=x, d=d*fke)
  
  m=length(par)
  lamb=par[1]
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

  return(-likeli+c0*sum(abs(alpha[-1])))
}

##########################
# Data generation process:
# Sample size = 200
# P(Y^{mi}=1)=0.25 # Y^{mi} = \mathcal{M}
# P(Y^{nmi}=1)=0.2 # Y^{nmi} = \mathcal{N}
# Y^{obs}=Y^{mi}+Y^{nmi}
# P(D^{no.nmi}=0 | Y^{mi}=0) = 0.8 (q)
# Size of fever killing effect (beta)= 0.2 which means 80% fever killing effect
# True MAFF = 0.5
#############################
dataset=read.csv("simulated_dataset.csv")[,-1]
n=length(dataset[,1])

attach(dataset)
####
x.grid=sort(unique(d.obs))
theta.space=seq(0, 30, by=1)
theta.space.sub=theta.space[-1]
d.grid=theta.space

## df = number of parameters in g(x)
df=5

freq.dobs=freq.count(data=cbind(d.obs, y.obs), x=x.grid, d=d.grid)

q.mat=ns(d.grid[-1], df=df-1)
for(j in 1:length(q.mat[1,])){
  q.mat[,j]=(q.mat[,j]-mean(q.mat[,j]))/sd(q.mat[,j])
}

## When the size of fever killing effect (beta) is known. True beta is 0.2.  
t1=Sys.time()
penal.est.par=optim(rep(0.5, df+2),penal.poisson.fke, method="L-BFGS-B", lower=c(0.01, 0.01, rep(-Inf, df)), upper=c(0.99, 0.99, rep(Inf, df)), control=list(trace=TRUE), x=x.grid, d=d.grid, freq.mat=freq.dobs, mat=q.mat, fke=0.2, c0=5)

t2=Sys.time()
t2-t1
penal.est.par$par; # this output has a length of df+2. The first component is an estimate of lambdastar and the second component is an estimate of proportion of zero (i.e., q). 

## Plot the densities g1 and g2. 
plot(d.grid, g.alpha1(alpha=penal.est.par$par[2:(df+1)], mat=q.mat), type="l", ylim=c(0,1), xlab="Parasite Density", ylab="")
lines(d.grid, g.alpha2(alpha=penal.est.par$par[3:(df+1)], eta=penal.est.par$par[(df+2)], d=d.grid[-1], mat=q.mat), col="blue")

## Adjustment step for estimating MAFF(lambda) from lambstar and p.
lambstar=penal.est.par$par[1]
prob.yobs=mean(y.obs==1)
est.lamb=(lambstar*(1-prob.yobs))/(1-lambstar*prob.yobs)
est.lamb # MAFF estimate (the true MAFF is 0.5)

############
## Sensitivity analysis
## when beta is not known. 
## We investigate a range of beta. Here, we assume that beta is in (0.05, 1)
fke.vec=seq(0.05,1,by=0.05)
est.lamb=rep(NA, length(fke.vec))

for(k in 1:length(fke.vec)){
  penal.est.par=optim(rep(0.5, df+2),penal.poisson.fke, method="L-BFGS-B", lower=c(0.01, 0.01, rep(-Inf, df)), upper=c(0.99, 0.99, rep(Inf, df)), control=list(trace=TRUE), x=x.grid, d=d.grid, freq.mat=freq.dobs, mat=q.mat, fke=fke.vec[k], c0=5)
  
  lambstar=penal.est.par$par[1]
  prob.yobs=mean(y.obs==1)
  est.lamb[k]=(lambstar*(1-prob.yobs))/(1-lambstar*prob.yobs)
  print(k)
}

plot(100*(1-fke.vec), est.lamb, type="l", xlab="fever killing effect (%)", ylab="MAFF", ylim=c(0.4, 0.6))
abline(h=0.5, lty=2)
