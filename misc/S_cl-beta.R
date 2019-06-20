#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###      F_select: SELECT structural model for paired gear selective analysis 
###                           User routine
###                           Juan Santos - 10.2014
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#input data:

X<-read.csv("/home/santos/Documents/7_SOFTWARE/selR/R_select/X_bacoma_on.csv",sep=",")

load("/home/santos/Documents/7_SOFTWARE/selR/R_select/.RData")

f_bwhv(X,side=NULL,resampling=T)
  
  
str(x)

#'data.frame':	106 obs. of  6 variables:
#  $ h   : num  1 1 1 1 1 1 1 1 1 1 ...
#$ l   : num  23 25 25.5 26 26.5 27 27.5 28 28.5 29 ...
#$ ncd1: int  0 0 0 0 0 0 0 0 0 0 ...
#$ ncd2: int  1 3 2 2 2 3 2 1 3 4 ...
#$ qcd1: num  0.0816 0.0816 0.0816 0.0816 0.0816 ...
#$ qcd2: num  0.049 0.049 0.049 0.049 0.049 ...

Start<-read.csv("/home/santos/Documents/7_SOFTWARE/selR/R_select/start_values.csv")[,-1]


model_selection<-F_select(x=x,lp=c(10,100),modelo=0,n.min=1,Start=Start,fix_par=list("sp",.5))

best_model<-F_select(x=x,lp=c(10,100),modelo=2,ref_length=35,n.min=1,Start=Start)

F_check(dat=x,mod=best_model)

P_select(mod=best_model,nome="",Dir=NULL)


#Bootstrap

Start_on<-read.csv("/home/santos/Documents/2_CRUISES/2019/SOLEA/SO759/9_software/start_on.csv",row.names=1)

X_on<-read.csv("/home/santos/Documents/7_SOFTWARE/selR/R_select/X_bacoma_on.csv")

x_on<-f_bwhv(X=X_on,resampling=F)

model_selection2<-F_select(x=x_on,lp=c(10,100),modelo=0,n.min=1,Start=Start_on,fix_par=list("sp",0.5))

best_model<-F_select(x=x_on,lp=c(10,100),modelo=1,ref_length=35,n.min=1,Start=Start_on)

best_model2<-F_select(x=x_on,lp=c(10,100),modelo=1,ref_length=35,n.min=1,Start=Start_on,fix_par=list("sp",0.5))

F_check(dat=x_on,mod=best_model2)

P_select(mod=best_model2,nome="",Dir=NULL)
lines(f_cl(par=c( 0.99000, 41.88703, 12.18668, 25.00000 , 5.00000),best_model2$L)~best_model2$L,
      col=1,lwd=2)


B_select(X=X_on,modelo=1,ref_length=35,lp=c(10,65),n.min=1,fix_par=list("sp",.5),B=10,ncores=2)




# bacoma off --------------------------------------------------------------

Start_off<-read.csv("/home/santos/Documents/2_CRUISES/2019/SOLEA/SO759/9_software/start_off.csv",row.names=1)

X_off<-read.csv("/home/santos/Documents/2_CRUISES/2019/SOLEA/SO759/3_dataR/COD_bacoma_off.csv",sep=" ")


x_off<-f_dualB(X_off, resampling = F)


ms_fixedSP_off<-F_select(x=x_off,lp=c(10,100),modelo=0,n.min=1,Start=Start_off,fix_par=list("sp",.5))

ms_off<-F_select(x=x_off,lp=c(10,100),modelo=0,n.min=1,Start=Start_off,fix_par=NULL)

model_selection_split<-F_select(x=x_off,lp=c(10,100),modelo=0,n.min=1,Start=Start_off,fix_par=NULL)

mod_fixedSP_off<-F_select(x=x_off,lp=c(10,100),modelo=1,n.min=1,Start=Start_off,fix_par=list("sp",.5))
mod_off<-F_select(x=x_off,lp=c(10,100),modelo=1,n.min=1,Start=Start_off,fix_par=NULL)

F_check(dat=x,mod=best_model)

P_select(mod=mod_fixedSP_off,nome="",Dir=NULL)



# Benchmark boot ----------------------------------------------------------
B=2000

cont1a<-numeric(B)
cont2a<-numeric(B)

cont1b<-numeric(B)
cont2b<-numeric(B)

for(i in 1:B){
  
  cont1a[i]<-  as.numeric(apply(f_dualB(X_on, resampling = T),2,sum)[3])
  cont2a[i]<-as.numeric(apply(f_dualB(X_on, resampling = T),2,sum)[4])
  
  cont1b[i]<-  as.numeric(apply(f_bwhv(X_on, resampling = T),2,sum)[3])
  cont2b[i]<-as.numeric(apply(f_bwhv(X_on, resampling = T),2,sum)[4])
  
  
  
}


rbind(c(mu_ncd1a=mean(cont1a),mu_ncd1b=mean(cont1b),mu_ncd2a=mean(cont2a),mu_ncd2b=mean(cont2b)),
c(er_ncd1a=sd(cont1a),er_ncd1b=sd(cont1b),er_ncd2a=sd(cont2a),er_ncd2b=sd(cont2b)))




