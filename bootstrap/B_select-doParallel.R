#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###               B-ecc: Function for obtaining bootstrap CI for F_emg function
###                               V1
###                           Juan Santos
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

B_ecc.dp<-function(X,lp,B,ncores,poli=3,exportfuns,nstart=10){
  
  cat("This process is using", ncores, "cores from your computer \n", "be pacient and wait for the results" )

  require(doParallel)
  require(foreach)
  
  cl <- makeCluster(ncores)
  
  registerDoParallel(cl)



r <- foreach(i=1:B,.export=exportfuns) %dopar% {
  
  if(i==1) x<-f_bwhv(X,resampling=F) else  x<-f_bwhv(X,resampling=T)
  
  F_ecc(x,poli=poli,lp=lp,n.min=1,nstart=nstart,Check=FALSE)
  
}

stopCluster(cl)

#Extract pseudo-results
b.sl<-t(ldply(r,function(x) x$sl))

#original data
fun<-r[[1]]$fun
l<-r[[1]]$l
L<-r[[1]]$L
p.spanel<-r[[1]]$p.spanel


#Estimations
aic<-r[[1]]$aic
aicc<-r[[1]]$aicc
Betas.ci<-t(apply(b.Betas,2,function(x) F_bci(x,q=95)))
sl.ci<-t(apply(b.sl,1,function(x) F_bci(x,q=95)))

return(list(fun=fun,B=B,l=l,L=L,p.spanel=p.spanel,aic=aic,aicc=aicc,Betas.ci=Betas.ci,sl.ci=sl.ci))}












