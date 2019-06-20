#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###      F_select-funs: SELECT structural model for paired gear selective analysis 
###             V1: selectivity mdoels, MLE and auxiliary functions
###                           Juan Santos
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#  functions available 

# - Logit, probit  Gompertz and Richards

##Dependencies:

#RNo dependencies


# Selectivity functions ---------------------------------------------------



f_1<-function(par,l){
  
                    eta<-exp(log(9)*((l-par[2])/par[3]))
  
                    rl<-eta/(1+eta)
      
                    phi<-par[1]*rl/( (1-par[1]) + par[1]*rl)
  
                    return(phi)
                    
                    }

#probit

f_2<-function(par,l){
  
                     rl<-pnorm(1.349*(l-par[2])/par[3])
  
                     phi<-par[1]*rl/( (1-par[1]) + par[1]*rl)
  
                     return(phi)  
                    
                    }


#gompertz

f_3<-function(par,l){
  
                     rl<-exp(-exp(-(0.365+1.573*(l-par[2])/par[3])))
  
                     phi<-par[1]*rl/( (1-par[1]) + par[1]*rl)
  
                     return(phi)  
        
                    }


#richard


f_4<-function(par,l){
  
                    eta<-exp(log(0.5^par[4]/(1-0.5^par[4]))+(log(0.75^par[4]/(1-0.75^par[4])) - log(0.25^par[4]/(1-0.25^par[4])))*(l-par[2])/par[3])
  
                    rl<-eta/(1+eta)
  
                    rl<-rl^(1/par[4])
  
                    phi<-par[1]*rl/( (1-par[1]) + par[1]*rl)
  
                    return(phi)  
                  
                    }






# select model with ML target  --------------------------------------------



#logit

F_1<-function(par,x,...){
  
  l<-x$l
  
  ncd1<-x$ncd1
  
  ncd2<-x$ncd2
  
  qcd1<-x$qcd1
  
  qcd2<-x$qcd2
  
  qt<-qcd1/qcd2


  phi<-f_1(par,l)

  
  return(-sum(ncd1*log(qt*phi/(qt*phi+(1-phi)))+ ncd2*log((1-phi)/(qt*phi+((1-phi)))))) 

                        }

#probit

F_2<-function(par,x,...){
  
  l<-x$l
  
  ncd1<-x$ncd1
  
  ncd2<-x$ncd2
  
  qcd1<-x$qcd1
  
  qcd2<-x$qcd2
  
  qt<-unique(qcd1/qcd2)
  
  phi<-f_2(par,l)
  
  return(-sum(ncd1*log(qt*phi/(qt*phi+(1-phi)))+ ncd2*log((1-phi)/(qt*phi+((1-phi)))))) 
  
                        }

#gompertz

F_3<-function(par,x,...){
  
  l<-x$l
  
  ncd1<-x$ncd1
  
  ncd2<-x$ncd2
  
  qcd1<-x$qcd1
  
  qcd2<-x$qcd2
  
  qt<-unique(qcd1/qcd2)
  
  phi<-f_1(par,l)
  
  return(-sum(ncd1*log(qt*phi/(qt*phi+(1-phi)))+ ncd2*log((1-phi)/(qt*phi+((1-phi)))))) 
        
                      }

#richards

F_4<-function(par,x,...){
  
  l<-x$l
  
  ncd1<-x$ncd1
  
  ncd2<-x$ncd2
  
  qcd1<-x$qcd1
  
  qcd2<-x$qcd2
  
  qt<-unique(qcd1/qcd2)
  
  phi<-f_4(par,l)
  
 return(-sum(ncd1*log(qt*phi/(qt*phi+(1-phi)))+ ncd2*log((1-phi)/(qt*phi+((1-phi)))))) 

                      }




# F for predictions  ------------------------------------------------------



f_phi.hat<-function(modelo,par.hat,l,qr){
  
      phi.hat<-if(modelo==1){
        
                              phi<- f_1(par=par.hat,l=l)
         
                              phi<-qr*phi/(qr*phi+(1-phi))
      
                              } else{
                                  
                                 if(modelo==2){
                
                                               phi<- f_2(par=par.hat,l=l)
                
                                               phi<-qr*phi/(qr*phi+(1-phi))
                
                                                                           } else{
                 
                                                                                   if(modelo==3){
                                                                                     
                                                                                                 phi<- f_3(par=par.hat,l=l)
                                                                                     
                                                                                                 phi<-qr*phi/(qr*phi+(1-phi))
                                                                                              
                                                                                                  }  else{
                                                                                                    
                                                                                                          phi<- f_4(par=par.hat,l=l)
                                                                                                    
                                                                                                          phi<-qr*phi/(qr*phi+(1-phi))
                                                                                                           }
                                                                                        }
                                                   }
                   return(phi.hat) }




# F to identify and add name of size selection function  ------------------




funnome<-function(modelo){
  
                          nome<-if(modelo==1){"logit"} else{
      
                                                            if(modelo==2){"probit"} else{
                                                              
                                                                                          if(modelo==3){"gompertz"}  else  {"richards"}
                                                                                              
                                                                                         }
   
                                                             }
        return(nome)
}

  
  
  

# F to calculate selectivity indicators -----------------------------------




f_si<-function(ref_length,l,test,pop){
  
                                 
                                       v_ref<-ifelse(l<ref_length, "below","above") 
    
                                       tab.ind<-aggregate(cbind(test,pop), by=list(v_ref),FUN=sum)
    
                                       tab.ind<-rbind(tab.ind,data.frame(Group.1="total",test=sum(tab.ind[,2]),pop=sum(tab.ind[,3])))
    
                                       sel_ind<-c(ref_length,round(100*tab.ind[,2]/(tab.ind[,3]),2),round(100*tab.ind[2,2]/tab.ind[3,2],2))
    
                                       names(sel_ind)<-c("ref","nPa","nPb","nPt","nDr")
    
                                       return(sel_ind) 
    
                                      }
  

