
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###
###            F_select-check: Fit statistic and Information values 
###                for paired gear selective models class SELECT 
###                          Main function (V3)
###                        Juan Santos - 10.2014
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Arguments:

 # dat: data modeled, consisting of:

#'data.frame':	141 obs. of  4 variables:
# $ l   : num  18.5 19.5 19.8 20.2 20.5 ...
# $ ncd1: num  0 1 1 1 3 2 5 5 1 6 ...
# $ ncd2: num  1 0 0 1 0 0 0 0 1 3 ...

 # mod: Model from class selR::select models.


# Dependencies:

#No dependencies


F_check<-function(dat,mod){
  
       
                            nbetas<-mod[["nbetas"]]
                            
                            modelhood<-mod[["modelhood"]]
                            
                            phi<-mod[["phi_check"]]
  
                            #data
                            
                            l<-dat$l
                            
                            nl<-length(l)
                            
                            
                            ncd1<-dat$ncd1
                            
                            ncd2<-dat$ncd2
                            
                            N<-ncd1+ncd2
                            
                            ncd1.hat<-N*phi
                            

                              
                              if( any((N==0)&(ncd1>0 | ncd1.hat>0)) ) stop("Wrong data for function Check")
                            
                              if( any((ncd1.hat == 0 & N >0) | (ncd1.hat == N & ncd1 < N)) )  stop("Impossibility in function Check")
                            
                              p=ifelse(N>0,ncd1/N,0.5)
                              
                              phat=ifelse(N>0,ncd1.hat/N,0.5)
                              
                              #Deviance Residuals
                              resDeviance0<-ifelse(ncd1>0,p*log(p/phat),0)+ifelse(ncd1<N,(1-p)*log((1-p)/(1-phat)),0)
                              
                              resDeviance1<-sqrt(2*N*resDeviance0)
                              
                              sign=ifelse(ncd1>=ncd1.hat,1,-1)
                              
                              resDeviance<-sign*resDeviance1
                              
                              #Deviance
                              Deviance=round(sum(resDeviance^2),4)
                              
                              pvalue=1-pchisq(Deviance,df=(nl-nbetas))
                              
                              #Chi square value
                              #Pearson Residual
                              resPearson=sqrt(N)*ifelse(is.na((p-phat)/(sqrt(phat*(1-phat)))), 0, (p-phat)/(sqrt(phat*(1-phat))))
                              
                              Chi=round(sum(resPearson^2),4)
                              
                              # AIC's
                              aic=2*(nbetas)+2*modelhood
                              
                              aicc=aic+2*nbetas*(nbetas+1)/(nl-nbetas-1)
                              
                              
                              
                              if(pvalue<0.0001) pvalue=">0.001" else pvalue=pvalue
  
  return(list(resPearson=resPearson,resDeviance=resDeviance,Deviance=Deviance,AIC=aic,AICc=aicc,PearsonChi=Chi,pvalue=pvalue,df=c((nl-nbetas),length(N)))) }

