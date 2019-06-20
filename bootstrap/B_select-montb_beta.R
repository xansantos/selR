#########################################                                      #
#   TO INVESTIGATE SOURCES OF VARIATION
#       B_selec for resampling          #
#         on selectivity data           #
#########################################





# Functions for two-nested blocks montecarlo resampling method -----------------------------------------------------------

#Within haul variation function(nested in f_bwhv)
      
      # Arguments:
      # U= individual haul data from Ltabs$Lw
#U<-X.1[X.1$haul==1,]

f_blockb<-function(U){
  
  require(plyr)
  
  haul=unique(U$haul)
  qcd2<-if(length(unique(U$qcd2))==1) unique(U$qcd2) else unique(U$qcd2)[1]
  qcd1<-if(length(unique(U$qcd1))==1) unique(U$qcd1) else unique(U$qcd1)[1]
  
  ncd2<-with(U,rep(l,ncd2))
  ncd1<-with(U,rep(l,ncd1))
  
  ncd2.b<-sample(ncd2,length(ncd2),replace=T)
  ncd1.b<-sample(ncd1,length(ncd1),replace=T)
  
  t.ncd2<-table(ncd2.b)
  t.ncd1<-table(ncd1.b)
  
  df.ncd2<-data.frame(l=type.convert(as.character(names(t.ncd2))),ncd2=as.vector(t.ncd2))
  df.ncd1<-data.frame(l=type.convert(as.character(names(t.ncd1))),ncd1=as.vector(t.ncd1))
  
  #df.cc$cc<-as.integer(df.cc$cc/qcd2)
  #df.cd$cd<-as.integer(df.cd$cd/qcd1)
  
  df.wh<-join(df.ncd2,df.ncd1,by="l",type="full")
  df.wh<-df.wh[order(df.wh$l),]
  df.wh$haul<-haul
  df.wh<-df.wh[,c(1,3,2)]
  df.wh[is.na(df.wh)]<-0
  df.wh$ncd2=df.wh$ncd2/qcd2
  df.wh$ncd1=df.wh$ncd1/qcd1
  
  return(df.wh)}

#Full resampling function:
      
      # Arguments:
      # X=data set from Ltabs$Lw

f_blockB<-function(X,resampling=T){
  
  
    
    names(X)<-c("haul","l","ncd1","ncd2","qcd1","qcd2") 
    

  
  #haul vector 
  
  haul<-unique(X$haul)  
  if(resampling==FALSE) h.b=haul else {h.b<-sample(haul,length(haul),replace=T)}
  
  #to define haul vector and bootstrap iterations
  
  b1<-1:length(h.b)
  
  #first resampling stage: between haul variation
  
  X.1<-ldply(b1,function(i){x<-X[X$haul==h.b[i],];x$haul<-i;return(x)})
  
  
  
  #second resampling stage: within haul variation
  if(resampling==FALSE) {X.2<-X.1;X.2$ncd2<-X.2$ncd2/X.2$qcd2;X.2$ncd1<-X.2$ncd1/X.2$qcd1}else{X.2<-ddply(X.1,.(haul),f_blockb)}
  
  #pool data
  X.bw<-aggregate(cbind(ncd1,ncd2)~l,data=X.2,sum)
  X.bw$qcd2<-1
  X.bw$qcd1<-1
  
  X.bw<-data.frame(h=1,X.bw)
  return(X.bw)}


# Functions for two-nested levels montecarlo resampling method -----------------------------------------------------------

# resampling function for individual samples
# Dependency: NONE

f_selb<-function(X){
  
  nr<-nrow(X)
  
  tot=X$ncd1+X$ncd2
  
  p=X$ncd1/tot
  
  ncd1_b<-sapply(1:nr, function(i){rbinom(1,tot[i],p[i])})
  
  ncd2_b<-tot-ncd1_b
  
  X$ncd1<-ncd1_b
  
  X$ncd2<-ncd2_b

  return(X)}


## Main function
# Dependency: f_resfish

f_selB<-function(X,resampling=F){
  
  #Arguments
  
  #X: data frame with the following structure:
  
  #'data.frame':	235 obs. of  6 variables:
  #  $ h   : int  2 2 2 2 2 2 2 2 2 5 ...
  #$ l   : int  21 24 25 26 28 30 31 34 35 16 ...
  #$ ncd1: int  0 0 0 0 0 1 0 0 1 0 ...
  #$ ncd2: int  1 1 1 1 1 2 1 1 0 1 ...
  #$ qcd1: int  1 1 1 1 1 1 1 1 1 1 ...
  #$ qcd2: int  1 1 1 1 1 1 1 1 1 1 ...
  
  
  #resampling: Character (NA,full","haul","sample")
  
  names(X)<-c("haul","l","ncd1","ncd2","qcd1","qcd2") 
    
  
  if(is.na(resampling)==T) {
    
                        X_b<-X
    
                        X_b$ncd1<-X_b$ncd1/X_b$qcd1
    
                        X_b$ncd2<-X_b$ncd2/X_b$qcd2
    
                          } else { 
    
                                   #haul vector 
                                   haul<-unique(X$haul)
                          
                                   h.b<-sample(haul,length(haul),replace=T)
  
                                   b1<-1:length(h.b)
    
                                   # First level: resampling hauls  
                                   X_temp<-ldply(b1,function(i){
                          
                                                                    x<-X[X$haul==h.b[i],]
                                
                                                                    x$haul<-i
                                
                                                                    return(x)
                                                                    
                                                                    }
                                                    )
  
                                  # second levels: resampling measured fish
                                  X_b<-f_selb(X_temp)
    
                                  X_b$ncd1<-X_b$ncd1/X_b$qcd1
    
                                  X_b$ncd2<-X_b$ncd2/X_b$qcd2
    
                                  }
  
  #pooling data over pseudo-hauls
  X_out<-aggregate(cbind(ncd1,ncd2)~l,data=X_b,sum)
  
  # Set raising factors to 1.0 (samples already sampled)
  X_out$qcd1<-1
  
  X_out$qcd2<-1

  X_out<-data.frame(h=1,X_out)
  
  return(X_out)}



# Functions for BOOTSTRAP inference -----------------------------------------------------------


##f_eci: Efron-bootstrap confidence intervals function 

    # Arguments:
    # x= vector b={1,...,B} from the montecarlo resampling scheme
    # q= quantile*100
    

    f_bci<-function(x,q){
  xm<-x[1]
  xci<-quantile(x,probs=c((1-q/100)/2,1-((1-q/100)/2)))
  res<-c(xci[1],Mean.Val=xm,xci[2])
  return(res)}

##f_tbci: t-bootstrap confidence intervals function 



# Others -----------------------------------------------------------

#f_getclogitCI: estimate CIs for clogit curves (dependent of f_bci). 

    # Arguments:
    # x= vector with the pseudo-curve estimations 1,...,b
    # q= quantile


f_getclogitCI<-function(X,q){
  
  rl<-function(L,l50,sr){
    eta<-exp(log(9)*(L-l50)/sr)  
    return(eta/(1+eta))
  }
  
  f_curve<-function(L,x){
    C=unlist(x[1])
    l50<-unlist(x[2])
    sr<-unlist(x[3])
    
    eg<-C*(1-rl(L,l50,sr))
    rg<-1-eg
    curve<-.5*rg/(.5*rg+.5)
    
    #curve<-C*(1-rl(L,l50,sr))
    
    return(curve)}
  
  
  L<-seq(0,100,by=.5)
  
  
  
  df.matcurve<-apply(X[[1]],1,function(x){f_curve(L,x)})
  
  df.curve<-data.frame(t(apply(df.matcurve,1,function(x){f_bci(x,q=q)})))
  df.curve$L<-L
  
  return(list(df.curve,df.matcurve))
  
  
}


#f_getclogitCI: estimate CIs for clogit curves (dependent of B_crates). 

# Arguments:
# x= object [[2]] from  B_crates object


f_bracketsCI<-function(x){
  
  x<-t(x)
  vx<-apply(x,1,function(y){paste(round(y[2],2),"(",round(y[1],2),"-",round(y[3],2),")",sep="")})
  return(vx) 
  
}



# Bootstrap estimation of catch rates in test codend relative to the reference codend

#Arguments:
# X= Ltabs list
# sp= Species names
# setup= gear setup
# mrs= minimun reference size: if mrs = 0, then rates above and below mrs not estimated
# B = number of bootstrap iterations


B_crates<-function(X,sp,setup,mrs,B){
  require(plyr)
  
  if(!exists("f_whv")==TRUE | !exists("f_bwhv")==TRUE) stop("Montecarlo functions not loaded in the current Working space:\nPlease run the lines in the source section.")
  
  
  X<-f_semdat(X,sp=sp,setup=setup)  
  
  
  vb<-1:B
  
  BOOT<-ldply(vb,function(b){
    if(b==1) X.b=f_bwhv(X,resampling=F) else X.b<-f_bwhv(X,resampling=T)
    #X.b<-X.b[X.b$test+X.b$ref>1,]
    
    vtot<-apply(X.b[,c(2,3)],2,sum)
    crt<-as.numeric(100*vtot[1]/vtot[2])
    
    if(mrs==0){ 
      RES<-c(crt,crt,crt)
      names(RES)<-c("total","above","below")
      
    }else{
      
      X.b$mrs<-ifelse(X.b$l<mrs,"below","above")      
      #rate above below
      RES<-c(crt,ddply(X.b,.(mrs),function(x){v<-apply(x[,c(2,3)],2,sum);cr<-as.numeric(100*v[1]/v[2]);return(cr)})[,2])
      names(RES)<-c("total","above","below")
      
    }
    
    
    
    
    cat("\n iteration=", b)
    return(RES)})
  
  CI<-apply(BOOT,2,function(x){f_bci(x,q=95)})
  
  RES<-list(BOOT,CI)
  
  return(RES)}






