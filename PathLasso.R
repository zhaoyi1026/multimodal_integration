#########################################
# soft-thresholding function
soft.thred<-function(mu,omega)
{
  if(length(mu)==1)
  {
    return(max(abs(mu)-omega,0)*sign(mu))
  }else
  {
    return(pmax(abs(mu)-omega,rep(0,length(mu)))*sign(mu))
  }
}

# solution lemma for pathway lasso
pathlasso.sol<-function(lambda,omega,phi1,phi2,mu1,mu2)
{
  if(lambda==0)
  {
    a<-soft.thred(mu1,omega)/phi1
    b<-soft.thred(mu2,omega)/phi2
  }else
  {
    C1<-((phi2*mu1-lambda*mu2)>omega*(phi2-lambda))&((phi1*mu2-lambda*mu1)>omega*(phi1-lambda))
    C2<-((phi2*mu1+lambda*mu2)>omega*(phi2-lambda))&((phi1*mu2+lambda*mu1)<(-omega*(phi1-lambda)))
    C3<-((phi2*mu1+lambda*mu2)<(-omega*(phi2-lambda)))&((phi1*mu2+lambda*mu1)>omega*(phi1-lambda))
    C4<-((phi2*mu1-lambda*mu2)<(-omega*(phi2-lambda)))&((phi1*mu2-lambda*mu1)<(-omega*(phi1-lambda)))
    C5<-((abs(mu1)>omega)&((phi1*abs(mu2)-lambda*abs(mu1)<=omega*(phi1-lambda))))
    C6<-((abs(mu2)>omega)&((phi2*abs(mu1)-lambda*abs(mu2)<=omega*(phi2-lambda))))
    
    x1<-phi2*(mu1-omega)-lambda*(mu2-omega)
    x2<-phi2*(mu1-omega)+lambda*(mu2+omega)
    x3<-phi2*(mu1+omega)+lambda*(mu2-omega)
    x4<-phi2*(mu1+omega)-lambda*(mu2+omega)
    
    y1<-phi1*(mu2-omega)-lambda*(mu1-omega)
    y2<-phi1*(mu2+omega)+lambda*(mu1-omega)
    y3<-phi1*(mu2-omega)+lambda*(mu1+omega)
    y4<-phi1*(mu2+omega)-lambda*(mu1+omega)
    
    de<-phi1*phi2-lambda^2
    
    if(C1)
    {
      a<-x1/de
      b<-y1/de
    }else
      if(C2)
      {
        a<-x2/de
        b<-y2/de
      }else
        if(C3)
        {
          a<-x3/de
          b<-y3/de
        }else
          if(C4)
          {
            a<-x4/de
            b<-y4/de
          }else
            if(C5)
            {
              a<-(abs(mu1)-omega)*sign(mu1)/phi1
              b<-0
            }else
              if(C6)
              {
                a<-0
                b<-(abs(mu2)-omega)*sign(mu2)/phi2
              }else
              {
                a<-0
                b<-0
              }
  }
  
  return(list(a=a,b=b))
}
#########################################

#########################################
# two-block pathway lasso
library("glmnet")

# log-likelihood function
log.Lik<-function(X,M1,M2,Y,beta,theta,zeta,pi,Lambda,delta)
{
  l1<-sum(diag(t(M1-X%*%beta)%*%(M1-X%*%beta)))*(-1/2)
  l2<-sum(diag(t(M2-X%*%zeta-M1%*%Lambda)%*%(M2-X%*%zeta-M1%*%Lambda)))*(-1/2)
  l3<-(t(Y-X*delta-M1%*%theta-M2%*%pi)%*%(Y-X*delta-M1%*%theta-M2%*%pi))[1,1]*(-1/2)
  
  re<-data.frame(M1=l1,M2=l2,Y=l3,sum=l1+l2+l3)
  rownames(re)<-"logLik"
  
  return(re)
}

# objective function
obj.func<-function(X,M1,M2,Y,beta,theta,zeta,pi,Lambda,delta,kappa1,kappa2,kappa3,kappa4,nu1,nu2,mu1,mu2)
{
  l1<-sum(diag(t(M1-X%*%beta)%*%(M1-X%*%beta)))
  l2<-sum(diag(t(M2-X%*%zeta-M1%*%Lambda)%*%(M2-X%*%zeta-M1%*%Lambda)))
  l3<-(t(Y-X*delta-M1%*%theta-M2%*%pi)%*%(Y-X*delta-M1%*%theta-M2%*%pi))[1,1]
  
  l<-l1+l2+l3
  
  p1<-kappa1*(sum(abs(beta)*abs(theta))+sum(nu1*(beta^2))+sum(nu1*(theta^2)))
  p2<-kappa2*(sum(abs(zeta)*abs(pi))+sum(nu1*(zeta^2))+sum(nu1*(pi^2)))
  p3<-kappa3*(sum(abs(Lambda)))+kappa4*abs(delta)
  p4<-mu1*(sum(abs(beta))+sum(abs(theta)))+mu2*(sum(abs(zeta))+sum(abs(pi)))
  
  pen<-p1+p2+p3+p4
  
  return(l/2+pen)
}

# data standardized
pathlasso.2b.std<-function(X,M1,M2,Y,kappa1,kappa2,kappa3,kappa4,nu1=1,nu2=1,mu1=0,mu2=0,rho=1,max.itr=10000,tol=1e-6,rho.increase=FALSE,trace=FALSE,
                           beta0=NULL,theta0=NULL,zeta0=NULL,pi0=NULL,Lambda0=NULL)
{
  n<-length(X)
  p1<-ncol(M1)
  p2<-ncol(M2)
  
  if(trace)
  {
    beta.trace<-NULL
    theta.trace<-NULL
    zeta.trace<-NULL
    pi.trace<-NULL
    
    td.beta.trace<-NULL
    td.theta.trace<-NULL
    td.zeta.trace<-NULL
    td.pi.trace<-NULL
    
    Lambda.trace<-NULL
    delta.trace<-NULL
    
    tau1.trace<-NULL
    tau2.trace<-NULL
    tau3.trace<-NULL
    tau4.trace<-NULL
    
    objfunc<-NULL
  }
  
  if(rho.increase)
  {
    rho0<-rho
  }else
  {
    rho0<-0
  }
  
  # set initial values
  if(is.null(beta0))
  {
    td.beta0=beta0<-matrix(rep(0,p1),nrow=1)
  }else
  {
    td.beta0<-beta0
  }
  if(is.null(theta0))
  {
    td.theta0=theta0<-rep(0,p1)  
  }else
  {
    td.theta0<-theta0
  }
  if(is.null(zeta0))
  {
    td.zeta0=zeta0<-matrix(rep(0,p2),nrow=1)  
  }else
  {
    td.zeta0<-zeta0
  }
  if(is.null(pi0))
  {
    td.pi0=pi0<-rep(0,p2) 
  }else
  {
    td.pi0<-pi0
  }
  
  delta0<-0
  if(is.null(Lambda0))
  {
    Lambda0<-matrix(0,p1,p2) 
  }
  
  tau1=tau2<-rep(0,p1)
  tau3=tau4<-rep(0,p2)
  
  s<-0
  diff<-100
  
  time<-system.time(
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update beta
    beta.new<-(t(X)%*%M1-t(tau1)+rho*td.beta0)/((t(X)%*%X)[1,1]+rho)
    # update theta
    theta.new<-solve(t(M1)%*%M1+rho*diag(rep(1,p1)))%*%(t(M1)%*%(Y-X*delta0-M2%*%pi0)-tau2+rho*td.theta0)
    # update zeta
    zeta.new<-(t(X)%*%(M2-M1%*%Lambda0)-t(tau3)+rho*td.zeta0)/((t(X)%*%X)[1,1]+rho)
    # update pi
    pi.new<-solve(t(M2)%*%M2+rho*diag(rep(1,p2)))%*%(t(M2)%*%(Y-X*delta0-M1%*%theta.new)-tau4+rho*td.pi0)
    
    # update delta
    delta.new<-soft.thred((t(X)%*%(Y-M1%*%theta.new-M2%*%pi.new))[1,1],kappa4)/((t(X)%*%X)[1,1])
    
    # update Lambda
    Lambda.new<-matrix(NA,p1,p2)
    for(j in 1:p2)
    {
      fit.tmp<-glmnet(x=M1,y=M2[,j]-X*zeta.new[j],family="gaussian",lambda=kappa3,standardize=FALSE)
      Lambda.new[,j]<-coef(fit.tmp)[-1]
    }
    
    # update tilde.beta and tilde.theta
    td.beta.new<-matrix(rep(NA,p1),nrow=1)
    td.theta.new<-rep(NA,p1)
    for(k in 1:p1)
    {
      sol.tmp<-pathlasso.sol(lambda=kappa1,omega=mu1,phi1=2*kappa1*nu1+rho,phi2=2*kappa1*nu1+rho,mu1=tau1[k]+rho*beta.new[k],mu2=tau2[k]+rho*theta.new[k])
      td.beta.new[k]<-sol.tmp$a
      td.theta.new[k]<-sol.tmp$b
    }
    
    # update tilde.zeta and tilde.pi
    td.zeta.new<-matrix(rep(NA,p2),nrow=1)
    td.pi.new<-rep(NA,p2)
    for(j in 1:p2)
    {
      sol.tmp<-pathlasso.sol(lambda=kappa2,omega=mu2,phi1=2*kappa2*nu2+rho,phi2=2*kappa2*nu2+rho,mu1=tau3[j]+rho*zeta.new[j],mu2=tau4[j]+rho*pi.new[j])
      td.zeta.new[j]<-sol.tmp$a
      td.pi.new[j]<-sol.tmp$b
    }
    
    # update tau1, tau2, tau3, tau4
    tau1<-tau1+rho*c(beta.new-td.beta.new)
    tau2<-tau2+rho*c(theta.new-td.theta.new)
    tau3<-tau3+rho*c(zeta.new-td.zeta.new)
    tau4<-tau4+rho*c(pi.new-td.pi.new)
    
    if(trace)
    {
      beta.trace<-cbind(beta.trace,c(beta.new))
      theta.trace<-cbind(theta.trace,c(theta.new))
      zeta.trace<-cbind(zeta.trace,c(zeta.new))
      pi.trace<-cbind(pi.trace,c(pi.new))
      
      td.beta.trace<-cbind(td.beta.trace,c(td.beta.new))
      td.theta.trace<-cbind(td.theta.trace,c(td.theta.new))
      td.zeta.trace<-cbind(td.zeta.trace,c(td.zeta.new))
      td.pi.trace<-cbind(td.pi.trace,c(td.pi.new))
      
      Lambda.trace<-cbind(Lambda.trace,c(Lambda.new))
      delta.trace<-c(delta.trace,delta.new)
      
      tau1.trace<-cbind(tau1.trace,tau1)
      tau2.trace<-cbind(tau2.trace,tau2)
      tau3.trace<-cbind(tau3.trace,tau3)
      tau4.trace<-cbind(tau4.trace,tau4)
      
      objfunc<-c(objfunc,obj.func(X,M1,M2,Y,td.beta.new,td.theta.new,td.zeta.new,td.pi.new,Lambda.new,delta.new,kappa1,kappa2,kappa3,kappa4,nu1,nu2,mu1,mu2))
    }
    
    diff.beta<-max(abs(beta.new-beta0))
    diff.theta<-max(abs(theta.new-theta0))
    diff.zeta<-max(abs(zeta.new-zeta0))
    diff.pi<-max(abs(pi.new-pi0))
    diff.Lambda<-max(abs(Lambda.new-Lambda0))
    diff.delta<-abs(delta.new-delta0)
    
    diff<-max(c(diff.beta,diff.theta,diff.zeta,diff.pi,diff.Lambda,diff.delta))
    
    # print(data.frame(diff=diff,beta=diff.beta,theta=diff.theta,zeta=diff.zeta,pi=diff.pi,Lambda=diff.Lambda,delta=diff.delta,obj=objfunc[s]))
    
    beta0<-beta.new
    theta0<-theta.new
    zeta0<-zeta.new
    pi0<-pi.new
    Lambda0<-Lambda.new
    delta0<-delta.new
    td.beta0<-td.beta.new
    td.theta0<-td.theta.new
    td.zeta0<-td.zeta.new
    td.pi0<-td.pi.new
    
    rho<-rho+rho0
  })
  
  if(s>max.itr)
  {
    warning("Method does not converge!")
  }
  
  constraint1=constraint2<-matrix(NA,1,4)
  colnames(constraint1)=colnames(constraint2)<-c("beta=td.beta","theta=td.theta","zeta=td.zeta","pi=td.pi")
  constraint1[1,1]<-(max(abs(beta.new-td.beta.new))<tol)
  constraint1[1,2]<-(max(abs(theta.new-td.theta.new))<tol)
  constraint1[1,3]<-(max(abs(zeta.new-td.zeta.new))<tol)
  constraint1[1,4]<-(max(abs(pi.new-td.pi.new))<tol)
  constraint2[1,1]<-max(abs(beta.new-td.beta.new))
  constraint2[1,2]<-max(abs(theta.new-td.theta.new))
  constraint2[1,3]<-max(abs(zeta.new-td.zeta.new))
  constraint2[1,4]<-max(abs(pi.new-td.pi.new))
  constraint<-cbind(data.frame(t(constraint1)),data.frame(t(constraint2)))
  colnames(constraint)<-c("Satisfied","value")
  
  beta.est<-td.beta.new
  colnames(beta.est)<-paste0("M1.",1:p1)
  rownames(beta.est)<-"X"
  theta.est<-matrix(td.theta.new,ncol=1)
  rownames(theta.est)<-paste0("M1.",1:p1)
  colnames(theta.est)<-"Y"
  zeta.est<-td.zeta.new
  colnames(zeta.est)<-paste0("M2.",1:p2)
  rownames(zeta.est)<-"X"
  pi.est<-matrix(td.pi.new,ncol=1)
  rownames(pi.est)<-paste0("M2.",1:p2)
  colnames(pi.est)<-"Y"
  Lambda.est<-Lambda.new
  rownames(Lambda.est)<-paste0("M1.",1:p1)
  colnames(Lambda.est)<-paste0("M2.",1:p2)
  
  net.mat<-matrix(NA,1+p1+p2,p1+p2+1)
  rownames(net.mat)<-c("X",paste0("M1.",1:p1),paste0("M2.",1:p2))
  colnames(net.mat)<-c(paste0("M1.",1:p1),paste0("M2.",1:p2),"Y")
  net.mat[1,]<-c(beta.est,zeta.est,delta.new)
  net.mat[2:(p1+1),(p1+1):(p1+p2+1)]<-cbind(Lambda.est,theta.est)
  net.mat[(1+p1+1):(1+p1+p2),(p1+p2+1)]<-pi.est
  
  IE.M1.est<-c(beta.est)*theta.est
  names(IE.M1.est)<-paste0("M1.",1:p1)
  IE.M2.est<-c(zeta.est)*pi.est
  names(IE.M2.est)<-paste0("M2.",1:p2)
  IE.M1M2.est<-matrix(NA,p1,p2)
  rownames(IE.M1M2.est)<-paste0("M1.",1:p1)
  colnames(IE.M1M2.est)<-paste0("M2.",1:p2)
  for(k in 1:p1)
  {
    for(j in 1:p2)
    {
      IE.M1M2.est[k,j]<-beta.est[k]*Lambda.est[k,j]*pi.est[j]
    }
  }
  
  if(trace)
  {
    re<-list(beta=beta.est,theta=theta.est,zeta=zeta.est,pi=pi.est,Lambda=Lambda.est,delta=delta.new,para.mat=net.mat,IE.M1=IE.M1.est,IE.M2=IE.M2.est,IE.M1M2=IE.M1M2.est,
             logLik=log.Lik(X,M1,M2,Y,td.beta.new,td.theta.new,td.zeta.new,td.pi.new,Lambda.new,delta.new),converge=(s<=max.itr),constraint=constraint,time=time,
             beta.trace=beta.trace,theta.trace=theta.trace,zeta.trace=zeta.trace,pi.trace=pi.trace,
             td.beta.trace=td.beta.trace,td.theta.trace=td.theta.trace,td.zeta.trace=td.zeta.trace,td.pi.trace=td.pi.trace,
             Lambda.trace=Lambda.trace,delta.trace=delta.trace,objfunc=objfunc)
  }else
  {
    re<-list(beta=beta.est,theta=theta.est,zeta=zeta.est,pi=pi.est,Lambda=Lambda.est,delta=delta.new,para.mat=net.mat,
             logLik=log.Lik(X,M1,M2,Y,td.beta.new,td.theta.new,td.zeta.new,td.pi.new,Lambda.new,delta.new),converge=(s<=max.itr),constraint=constraint,time=time)
  }
  
  return(re)
}

pathlasso.2b<-function(X,M1,M2,Y,kappa1,kappa2,kappa3,kappa4,nu1=1,nu2=1,mu1=0,mu2=0,rho=1,standardize=TRUE,max.itr=10000,tol=1e-6,rho.increase=FALSE,trace=FALSE,
                       beta0=NULL,theta0=NULL,zeta0=NULL,pi0=NULL,Lambda0=NULL)
{
  n<-length(X)
  p1<-ncol(M1)
  p2<-ncol(M2)
  
  if(standardize)
  {
    # standardize data
    X.sd<-sd(X)
    M1.sd<-apply(M1,2,sd)
    M2.sd<-apply(M2,2,sd)
    Y.sd<-sd(Y)
    
    X.std<-scale(X,center=TRUE,scale=TRUE)
    M1.std<-scale(M1,center=TRUE,scale=TRUE)
    M2.std<-scale(M2,center=TRUE,scale=TRUE)
    Y.std<-scale(Y,center=TRUE,scale=TRUE)
    
    re.std<-pathlasso.2b.std(X.std,M1.std,M2.std,Y.std,kappa1=kappa1,kappa2=kappa2,kappa3=kappa3,kappa4=kappa4,nu1=nu1,nu2=nu2,mu1=mu1,mu2=mu2,rho=rho,max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=trace,
                             beta0=beta0,theta0=theta0,zeta0=zeta0,pi0=pi0,Lambda0=Lambda0)
    
    beta.est<-re.std$beta*(M1.sd/X.sd)
    theta.est<-re.std$theta*(Y.sd/M1.sd)
    zeta.est<-re.std$zeta*(M2.sd/X.sd)
    pi.est<-re.std$pi*(Y.sd/M2.sd)
    Lambda.est<-re.std$Lambda
    for(j in 1:p2)
    {
      Lambda.est[,j]<-re.std$Lambda[,j]*(M2.sd[j]/M1.sd)
    }
    delta.est<-re.std$delta*(Y.sd/X.sd)
    
    net.mat<-matrix(NA,1+p1+p2,p1+p2+1)
    rownames(net.mat)<-c("X",paste0("M1.",1:p1),paste0("M2.",1:p2))
    colnames(net.mat)<-c(paste0("M1.",1:p1),paste0("M2.",1:p2),"Y")
    net.mat[1,]<-c(beta.est,zeta.est,delta.est)
    net.mat[2:(p1+1),(p1+1):(p1+p2+1)]<-cbind(Lambda.est,theta.est)
    net.mat[(1+p1+1):(1+p1+p2),(p1+p2+1)]<-pi.est
    
    IE.M1.est<-c(c(beta.est)*theta.est)
    names(IE.M1.est)<-paste0("M1.",1:p1)
    IE.M2.est<-c(c(zeta.est)*pi.est)
    names(IE.M2.est)<-paste0("M2.",1:p2)
    IE.M1M2.est<-matrix(NA,p1,p2)
    rownames(IE.M1M2.est)<-paste0("M1.",1:p1)
    colnames(IE.M1M2.est)<-paste0("M2.",1:p2)
    for(k in 1:p1)
    {
      for(j in 1:p2)
      {
        IE.M1M2.est[k,j]<-beta.est[k]*Lambda.est[k,j]*pi.est[j]
      }
    }
    
    re<-list(beta=beta.est,theta=theta.est,zeta=zeta.est,pi=pi.est,Lambda=Lambda.est,delta=delta.est,para.mat=net.mat,IE.M1=IE.M1.est,IE.M2=IE.M2.est,IE.M1M2=IE.M1M2.est,
             logLik=log.Lik(X,M1,M2,Y,beta.est,theta.est,zeta.est,pi.est,Lambda.est,delta.est),converge=re.std$converge,constraint=re.std$constraint,time=re.std$time,
             out.scaled=re.std)
  }else
  {
    re<-pathlasso.2b.std(X,M1,M2,Y,kappa1=kappa1,kappa2=kappa2,kappa3=kappa3,kappa4=kappa4,nu1=nu1,nu2=nu2,mu1=mu1,mu2=mu2,rho=rho,max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=trace,
                         beta0=beta0,theta0=theta0,zeta0=zeta0,pi0=pi0,Lambda0=Lambda0)
  }
  
  return(re)
}
#########################################