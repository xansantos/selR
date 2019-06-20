#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###      F_select: SELECT structural model for paired gear selective analysis 
###                           User routine
###                           Juan Santos - 10.2014
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#input data:

#X<-read.csv("/home/santos/Documents/7_SOFTWARE/selR/R_select/X_bacoma_on.csv",sep=",")
#f_bwhv(X,side=NULL,resampling=T)

load("/home/santos/Documents/7_SOFTWARE/selR/R_select/.RData")


  
  
str(x) # x is an example data!!!

#'data.frame':	106 obs. of  6 variables:
#  $ h   : num  1 1 1 1 1 1 1 1 1 1 ...
#$ l   : num  23 25 25.5 26 26.5 27 27.5 28 28.5 29 ...
#$ ncd1: int  0 0 0 0 0 0 0 0 0 0 ...
#$ ncd2: int  1 3 2 2 2 3 2 1 3 4 ...
#$ qcd1: num  0.0816 0.0816 0.0816 0.0816 0.0816 ...
#$ qcd2: num  0.049 0.049 0.049 0.049 0.049 ...

Start<-read.csv("/home/santos/Documents/7_SOFTWARE/selR/R_select/start_values.csv")[,-1]


# Fixed split:
model_selection_fixedSplit<-F_select(x=x,lp=c(10,100),modelo=0,n.min=1,Start=Start,fix_par=list("sp",.5))

# estimated split:
model_selection_EstimatedSplit<<-F_select(x=x,lp=c(10,100),modelo=0,n.min=1,Start=Start,fix_par=NULL)


# fixed or estimated split??
model_selection_fixedSplit[1,10]
model_selection_EstimatedSplit[1,10]


# Estimate best model (probit and fixed split)
best_model<-F_select(x=x,lp=c(10,100),modelo=2,ref_length=35,n.min=1,Start=Start,fix_par=list("sp",.5))

# Diagnosis statistics /residuals from the best model:
F_check(dat=x,mod=best_model)

# Figure with the average curve:
P_select(mod=best_model,nome="",Dir=NULL)



# Here starts the estimation of selectivity for the BACOMA_ON -----------------------------


Start_on<-read.csv("/home/santos/Documents/2_CRUISES/2019/SOLEA/SO759/9_software/start_on.csv",row.names=1)

COD_on<-read.csv("/home/santos/Documents/7_SOFTWARE/selR/R_select/X_bacoma_on.csv")

cod_on<-f_bwhv(X=X_on,resampling=F) # pool the data over hauls

# model selection procedure within models with fixed and estimated split:
ms_cod_on_fixedSplit<-F_select(x=cod_on,lp=c(10,100),modelo=0,n.min=1,Start=Start_on,fix_par=list("sp",0.5))
ms_cod_on_estimSplit<-F_select(x=cod_on,lp=c(10,100),modelo=0,n.min=1,Start=Start_on,fix_par=NULL)

# model selection between models with fixed and estimated split
ms_cod_on_fixedSplit[1,10]
ms_cod_on_estimSplit[1,10]

# Estimation of the best model for cod (bacoma_on)
mod_cod_on<-F_select(x=cod_on,lp=c(10,100),modelo=1,ref_length=35,n.min=1,Start=Start_on,fix_par=list("sp",0.5))


F_check(dat=cod_on,mod=mod_cod_on)

P_select(mod=mod_cod_on,nome="",Dir=NULL)



MOD_cod_on<-B_select(X=COD_on,modelo=1,ref_length=35,lp=c(10,100),n.min=1,fix_par=list("sp",.5),B=10,ncores=2)








