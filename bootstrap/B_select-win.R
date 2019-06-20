#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###               B_select: Function for obtaining bootstrap CI for Fgrid function
###                               V1
###                           Juan Santos
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

B_select_win<-function(X,Start,modelo,lp,n.min,fix_par,B,ncores,exportfuns){
  
cat("This process is using", ncores, "cores from your computer \n", "be pacient and wait for the results" )


  require(doParallel)

  require(foreach)

  require(reshape2)

  require(plyr)
  
  cl <- makeCluster(ncores)
  
  registerDoParallel(cl)
  
  
  r <- foreach(i=1:B,.export=exportfuns) %dopar% {
    
    x=if(i==1) f_selB(X,resampling=F) else {f_selB(X,resampling=T)}
    
    F_select(x=x,lp=lp,modelo=modelo,n.min=n.min,Start=Start,fix_par=fix_par)
    
  }
  
  stopCluster(cl)
#[1] "fun"           "Betas"         "nbetas"        "l"             "L"             "p.phi"         "selective_ind" "phi_check"     "phi"           "PHI"           "modelhood"     "aic"           "aicc"         

#Extract pseudo-results
b.phi<-t(ldply(r,function(x) x$PHI))

b.Betas<-ldply(r,function(x) x$Betas)

#b.Ind<-ldply(r,function(x) x$selective_ind)

#original data
fun<-r[[1]]$link

l<-r[[1]]$l

L<-r[[1]]$L

p.phi<-r[[1]]$p.phi

#Estimations
aic<-r[[1]]$aic

aicc<-r[[1]]$aicc

Betas.ci<-round(t(apply(b.Betas,2,function(x) f_bci(x,q=95))),3)

#Ind.ci<-round(t(apply(b.Ind[,-1],2,function(x) f_bci(x,q=95))),3)

phi.ci<-t(apply(b.phi,1,function(x) f_bci(x,q=95)))


return(list(fun=fun,B=B,l=l,L=L,p.phi=p.phi,aic=aic,aicc=aicc,Betas.ci=Betas.ci,phi.ci=phi.ci,betas=b.Betas)) }



















