
#Bfile=B_FG;lp=c(25,100);Cut=TRUE;medida="mm";txt="Analysis conducted on fresh measurments";nome.file="FiberGlassGrid";nome.plot="Fiber Glass Grid (Danish Drop Shape)";Dir=Dir.dr

Bp_grid<-function(Bfile,lp,Cut=TRUE,medida="cm",txt="",nome.file="BootPlot",nome.plot="BootPlot",Dir=NULL){
  
  em<-Bfile[[1]];B<-Bfile[[2]];l<-Bfile[[3]];L<-Bfile[[4]];p.rl<-Bfile[[5]];rl.ci<-Bfile[[9]]

#define layout matrix
  
  if(ncol(Bfile$betas)==2){
    
    nf <- matrix(c(1,2,3,3,3,3), 3, 2, byrow=TRUE)
    
  } else if(ncol(Bfile$betas)==3){
    
    nf <- matrix(c(1,2,3,4,4,4,4,4,4), 3, 3, byrow=TRUE)
  } else if(ncol(Bfile$betas)==4) {
    
    nf <- matrix(c(1,2,3,4,5,5,5,5), 4, 2, byrow=TRUE)
  }
  
  
  betas_df<-Bfile$betas
  betas_ci<-Bfile$Betas.ci
  
  
  

 
  if(Cut==T){
    rl.ci<-rl.ci[L>=min(l) & L<=max(l),]
    
    p.rl<-p.rl[L>=min(l) & L<=max(l)]
    L<-L[L>=min(l) & L<=max(l)]
  }
  
  
if(!is.null(Dir)){
          setwd(Dir)  
          pdf(paste(nome.file,".pdf",sep=""),width=12, height=9)
                  } 
            
         
  layout(nf,respect=F)
  
  for(i in 1:ncol(betas_df)){
    
    par(cex=1,cex.axis=.8,cex.main=1.2,cex.lab=.8,mar= c(5, 5, 4, .5))
    
    par_ci<-paste(betas_ci[i,2], " (",betas_ci[i,1],"-",betas_ci[i,3],")",sep="")
    
    hist(betas_df[,i],main=paste(toupper(names(betas_df)[i]),par_ci,sep="="),xlab="",col=i+1,las=1)
    
    abline(v=Bfile$Betas.ci[i,2],lwd=3,col=1,lty=2)
    
    # legend("topleft",paste("   Mean Value",par_ci,sep="\n"),bty="n",cex=1.3)
    
    
  }

  
  
          par(cex=1.5,cex.axis=1.5,cex.main=1.5,cex.lab=1.5,mar= c(5, 5, 3, .5))
          plot(rl.ci[,2]~L ,type = 'n',bty="n",ylim=c(0,1),xlim=lp,ylab="Retention probability r(l)",xlab=paste("Length"," (",medida,")",sep=""), cex=1.5,
               main=nome.plot,las=1)
          # add fill
          #points(Pt~L,pch=21,col=colors()[185],bg=colors()[234])
          polygon(c(rev(L), L), c(rev(rl.ci[,3]), rl.ci[,1]), col = rgb(.9, .1, .1,alpha=.3), border = NA)
          lines(rl.ci[,2]~ L,lty=1,lwd=3,col=1)
          lines(rl.ci[,1]~ L,lty = 2,lwd=2,col=rgb(.9, .1, .1,alpha=.3))
          lines(rl.ci[,3]~ L,lty = 2,lwd=2,col=rgb(.9, .1, .1,alpha=.3))
          points(p.rl~L,pch=21,col=rgb(.9, .1, .1,alpha=.3),cex=1.5,bg=rgb(.5, .5, .5,alpha=.8))
          mtext(paste("model ",em,sep=":"), side=3, adj=0, line=0,cex=1.5,col="darkgray")
          mtext(paste("Bootstrap CIs: ","(",B,"R",")",sep=""), side=3, adj=1, line=0,cex=1.5,col="darkgray")
          mtext(txt,side=1,line=-1,adj=.9,col=2,cex=.8)
          #abline(h=.5,col="darkgreen",lwd=1.5)
          
          if(!is.null(Dir)){dev.off()}  
          
          }
            
            
# Compare CIs -------------------------------------------------------------


F_ciComp<-function(Bmod1,Bmod2,gears,xL=c(0,100),nome.plot=NULL, Cor=F,posLenda=NULL){
  
  
  # prelude -----------------------------------------------------------------  
  gear1<-gears[1]
  gear2<-gears[2]
  
  if(is.null(nome.plot)) {nome.plot<-paste(gear1,gear2,sep=" Vs. ")}
  
  par(mfrow=c(1,1),las=1,cex.lab=1.5,cex.axis=1.5)
  
  with(Bmod1,plot(p.rl~L,type = 'n',bty="n",main=nome.plot,ylim=c(0,1),xlim=xL,ylab="Retention probability",
                  xlab="Length (mm)",col=1))
  #abline(h=0.5,lty=2,lwd=1.2,col=colors()[185])
  if(Cor==T){
    
    col1=rgb(red=0.15, green=0, blue=0.72, alpha=.5, maxColorValue = 1)
    col2=rgb(red=0.54, green=0.72, blue=0, alpha=.5, maxColorValue = 1) 
    
  }else{
   
    col1=rgb(red=0.80, green=0.80, blue=0.80, alpha=.7, maxColorValue = 1)
    col2=rgb(red=0.20, green=0.20, blue=0.20, alpha=.7, maxColorValue = 1)
    
  }

  
  rng1<-with(Bmod1,range(l))
  rng2<-with(Bmod2,range(l))
  dat.rates1<-as.data.frame(with(Bmod1,rl.ci[L>=rng1[1] & L<=rng1[2],]))
  names(dat.rates1)<-c("lci","m","uci")
  dat.rates1$l<-with(Bmod1,L[L>=rng1[1] & L<=rng1[2]])
  dat.rates2<-as.data.frame(with(Bmod2,rl.ci[L>=rng2[1] & L<=rng2[2],]))
  names(dat.rates2)<-c("lci","m","uci")
  dat.rates2$l<-with(Bmod2,L[L>=rng2[1] & L<=rng2[2]])
  
  
  
  with(dat.rates1,polygon(c(rev(l), l), c(rev(uci), lci), col = col1, border = 1,lty=1,lwd=2))
  with(dat.rates2,polygon(c(rev(l), l), c(rev(uci), lci), col = col2, border = 1,lty=2,lwd=2))
  
  if(!is.null(posLenda)){
    legend(posLenda,c(gear1,gear2),fill=c(col1,col2),bty="n",lwd=2,lty=c(1,2),cex=2) 
  }
  
  
}




F_npolygon<-function(x,xlim,medida="cm",pos1,pos2,nome,Dir=NULL){
  
  x$tot<-x[,2]+x[,3]

  ymx<-max(x$tot,na.rm = T)
  yst=if(ymx<10) 10 else {
    if(ymx>9 & ymx<25) 25 else{
      if(ymx>24 & ymx<50) 50  else{
        if(ymx>49 & ymx<75) 75 else{
          if(ymx>74 & ymx<100) 100 else{
            if(ymx>99 & ymx<150) 150  else{
              if(ymx>149 & ymx<300) 300  else{
                if(ymx>299 & ymx<500) 500  else{
                  if(ymx>499 & ymx<750) 750  else{
                    if(ymx>749 & ymx<1000) 1000 else{ 
                    
                    if(ymx>999 & ymx<1500) 1500  else{4000
                    }}}}}}}}}}}
  
  f_yax<-function(yst){ 
    stp=if(yst<10) 2 else {
      if(yst>=10 & yst<26) 5 else{
        if(yst>=26 & yst<75) 10  else{
          if(yst>=75 & yst<100) 25 else{
            if(yst>=100 & yst<150) 25 else{
              if(yst>=150 & yst<300) 50  else{
                if(yst>=300& yst<500) 100  else{ 250
                }}}}}}}
    
    return(stp)}
  
  if(!is.null(Dir)){
    
    setwd(Dir)  
    pdf(paste(nome,".pdf",sep=""),width=10, height=10)
    par(cex=3,cex.axis=3,cex.main=3,cex.lab=2.5,mar= c(5, 5, 4, 2),mfrow=c(1,3))
    
    
    #plot comparing total catch with codends 1 and 2 
    plot(tot~l,data=x,type="n",lwd=3,yaxt="n",lty=1,bty="n",xlab=paste("Length"," (",medida,")",sep=""),ylab="Catch numbers",
         ylim=c(0,yst),xlim=xlim)
    axis(2, at = seq(0, yst, by = f_yax(yst)),line=0)
    polygon(c(rev(x$l), x$l), c(rev(x$tot), rep(0,length(x$tot))),
            col=rgb(200/255, 200/255, 200/255,0.8), border = NA)
    lines(cd1~l,data=x,type="l",lwd=3)
    legend(pos1[1],pos1[2],"Total catch",fill=rgb(200/255, 200/255, 200/255,0.8),bty="n",cex=1.5)
    legend(pos2[1],pos2[2],"Codend 1",lty=1,bty="n",cex=3,lwd=1.5)  
   
    dev.off()
    
  } else {
    
    par(cex=1.5,cex.axis=1.5,cex.main=1.5,cex.lab=1.25,mar= c(5, 5, 4, 2))
    #plot comparing total catch with codends 1 and 2 
    plot(tot~l,data=x,type="n",lwd=3,yaxt="n",lty=1,bty="n",xlab=paste("Length"," (",medida,")",sep=""),ylab="Catch numbers",
         ylim=c(0,yst),xlim=xlim)
    axis(2, at = seq(0, yst, by = f_yax(yst)),line=0)
    polygon(c(rev(x$l), x$l), c(rev(x$tot), rep(0,length(x$tot))),
            col=rgb(200/255, 200/255, 200/255,0.8), border = NA)
    lines(cd1~l,data=x,type="l",lwd=3)
    legend(pos1[1],pos1[2],"    Total catch",fill=rgb(200/255, 200/255, 200/255,0.8),bty="n",cex=1)
    legend(pos2[1],pos2[2],"Codend 1",lty=1,bty="n",cex=1,lwd=1.5)  
    
    
  }
  
}




            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
          
