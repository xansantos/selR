#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###
###      F_select: SELECT structural model for paired gear selective analysis 
###                          Main function (V3)
###                        Juan Santos - 10.2014
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###Arguments:

#x=data frame with the following str:

            #'data.frame':	144 obs. of  3 variables:
            #  $ l   : num  18.2 19.2 19.8 20 20.2 ...
            #$ ncd1: num  1 1 1 1 1 2 2 1 2 5 ...
            #$ ncd2: num  0 0 0 0 0 0 0 0 1 0 ...
            #$ qcd1: num  0 0 0 0 0 0 0 0 1 0 ...
            #$ qcd2: num  0 0 0 0 0 0 0 0 1 0 ...

##Dependencies:

#R_select::S_select-funs


# Main Function -----------------------------------------------------------
#x
#lp=c(10,100)
#modelo=1
#ref_length=35
#n.min=1
#Start=Start
#fix_sp=c()

F_select<-function(x,lp,modelo,ref_length=NULL,n.min=1,Start,fix_par=NULL){
  
  require(plyr)
  
  # original data:
  
  # arrange the data for the model ----------------------
  
  # 'data.frame':	342 obs. of  11 variables:
  #  $ l   : num  18 18.5 18.8 19.5 19.8 ...
  #  $ ncd1: num  2 1 1 3 1 1 5 3 2 5 ...
  #  $ ncd2: num  0 0 0 0 0 0 0 0 0 0 ...
  #  $ qcd1: num  0 0 0 0 0 0 0 0 0 0 ...
  #  $ qcd2: num  0 0 0 0 0 0 0 0 0 0 ...

 modelos<-c("F_1","F_2","F_3","F_4") 
  
  
  # subset the data ------------------------------------------------------------------------
  
  x<-x[order(x$l),]
  
  x<-x[x$ncd1+x$ncd2>=n.min,]
  
  # To vectorize the data ---------------------------------------------------
  
  l=x$l
  
  ncd1=as.integer(x$ncd1)
  
  ncd2=as.integer(x$ncd2)
  
  qcd1=x$qcd1
  
  qcd2=x$qcd2
  
  qr<-unique(qcd1 / qcd2)
  
  pop<-ncd2 / qcd2
  
  test<-ncd1 / qcd1
  
  ntotal=test + pop
  
  p.phi<-test / (ntotal) 
  
  #Length Grid
  L<-seq(lp[1],lp[2],.25)
  
  res<-rep(NA,8)
  
  names(res)<-c("sp","l50","sr","d","logLik","nbetas","aic","aicc")


  
  

# Starting values ---------------------------------------------------------

if(!is.null(fix_par)) Start[,fix_par[[1]]]<-fix_par[[2]]  
  
  par<-as.numeric(Start[1,])
  
  lpar<-as.numeric(Start[2,])
  
  upar<-as.numeric(Start[3,])

  

# Models ------------------------------------------------------------------


  
  if(modelo==0){
    
      mod.l<- ldply(modelos,function(i){
      
      mod<-tryCatch(suppressWarnings(nlminb(par,eval(parse(text=i)), lower=lpar,upper=upar,x=x, control=list(abs.tol=0.0000000000000001))),error = function(e) {cat("")})   
      
      if(is.null(mod)){resultado=res} else {
        
        if(is.null(mod$par)|all(mod$par==par)|all(is.na(mod$par))|mod$objective==0) {resultado=res} else{
          
          parametros=mod$par
          
          ll=mod$objective
          
          nbetas=length(parametros)-length(which(parametros==par))
          
          aic.j=2*nbetas+2*ll
          
          aicc.j=aic.j+2*nbetas*(nbetas+1)/(nrow(x)-nbetas-1)
          
          aic.f<-round(aic.j,3)
          
          aicc.f<-round(aicc.j,3)
          
          resultado<-c(parametros,ll,nbetas,aic.f,aicc.f)
          
          names(resultado)<-names(res)
                                                                                                    }
                      }
      
      return(resultado)                 }  )
    
      mod_sel_table<-cbind(model=1:4,rl_fun=c("logit","probit","gompertz","richards"), mod.l)
      
      mod_sel_table<- mod_sel_table[order( mod_sel_table$aicc,decreasing=F),]
      
      out.res<-mod_sel_table
  
  }else{
      
    
    mod<-tryCatch(suppressWarnings(nlminb(par,eval(parse(text=modelos[modelo])), lower=lpar,upper=upar,x=x, control=list(abs.tol=0.0000000000000001))),error = function(e) {cat("")})
    
    if(!is.null(mod)){
      
      if(is.null(mod$par)|all(mod$par==par)|all(is.na(mod$par))|mod$objective==0) resultado=res else{
        
        parametros=mod$par
        
        nbetas=length(parametros)-length(which(parametros==par))
        
        ll=mod$objective
        
        nbetas=length(parametros)-length(which(parametros==par))
        
        aic.j=2*nbetas+2*ll
        
        aicc.j=aic.j+2*nbetas*(nbetas+1)/(nrow(x)-nbetas-1)
        
        aic.f<-round(aic.j,3)
        
        aicc.f<-round(aicc.j,3)
        
        resultado<-c(parametros,ll,nbetas,aic.f,aicc.f)
        
        names(resultado)<-names(res)
        
      }}
    
    
    par.hat<-as.numeric(parametros)
    
    phi_check<-f_phi.hat(modelo=modelo,par.hat=par.hat,l=l,qr=qr)
    
    phi<-f_phi.hat(modelo=modelo,par.hat=par.hat,l=l,qr=1)
    
    PHI<-f_phi.hat(modelo=modelo,par.hat=par.hat,l=L,qr=1)
    
    names(parametros)<-c("split","l50","sr","delta")
    
    
    if(is.null(ref_length)) selective_ind=NULL else selective_ind=f_si(ref_length,l,test,pop)
    
    
    out.res<- list(fun=funnome(modelo),
                   
                    Betas=round(parametros,4),
                   
                    nbetas=nbetas,
                   
                    l=l,
                   
                    L=L,
                   
                    p.phi=p.phi,
                   
                    selective_ind=selective_ind,
                   
                    phi_check=as.numeric(phi_check),
                   
                    phi=as.numeric(phi),
                    
                    PHI=as.numeric(PHI),
                   
                    modelhood=as.numeric(ll),
                   
                    aic=as.numeric(aic.f),
                   
                    aicc=as.numeric(aicc.f))  
    
    
    }
    

          


  return(out.res)}
