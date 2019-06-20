
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###
###            F_select-plot: Visualization of the average curve 
###                predicted by paired-gear selective models class SELECT 
###                          Main function (V3)
###                        Juan Santos - 10.2014
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Arguments:

  # mod: Model from class selR::select models.

  # nome: Main name of the plot

  # Dir if not NULL, figure will be sunk to the given directory


# Dependencies:

#No dependencies


#names(mod1)

#[1] "fun"       "Betas"     "l"         "p.phi"     "L"         "phi"       "modelhood" "aic"       "aicc" 


P_select<-function(mod,nome,Dir){
  
  fun_sel<-mod[["fun"]]

  par_sel<-paste("Split= ",round(mod[["Betas"]]["split"],2)," L50= ", round(mod[["Betas"]]["l50"],2)," SR= ", round(mod[["Betas"]]["sr"],2),sep="")
    

if(!is.null(Dir)){
  
  setwd(Dir)  
  
  pdf(paste(nome,".pdf",sep=""),width=12, height=9)
  
  par(cex=2,cex.axis=1.5,cex.main=1.5,cex.lab=1.5,mar= c(5, 4, 4, 1),mfrow=c(1,1))
  
  with(mod, plot(phi~l,type="n",bty="n",col=2,ylim=c(0,1),xlim=range(L),ylab="Catch sharing ",xlab="length (cm)"))
  
  abline(h=.5,lty=3,lwd=2,col="darkgreen")
  
  with(mod, points(p.phi~l,pch=21,col=2,cex=3,bg="darkgrey"))
  
  with(mod, lines(phi~l,type="l",lwd=3,col=2))
  
  mtext(par_sel, 3, line=-0.2,cex=2,col="darkgrey")
  
  mtext("Equal catch sharing", 3, line=-17.5,cex=1,col="darkgreen",adj=0)
  
  mtext(paste("model:", fun_sel,sep=""),3,line=0,cex=1.5,adj=0.05,padj=2,col="darkgrey",outer=T)
  
  dev.off()
  
} else {
  
  
  par(cex=2,cex.axis=1.5,cex.main=1.5,cex.lab=1.5,mar= c(5, 5, 4, 1),mfrow=c(1,1))
  
  with(mod, plot(phi~l,type="n",bty="n",col=2,ylim=c(0,1),xlim=range(l),ylab= expression(paste("Catch comparison ", phi, "(l)")),
                 xlab="length (cm)"))
  
  abline(h=.5,lty=3,lwd=2,col="darkgreen")
  
  with(mod, points(p.phi~l,pch=21,col=2,cex=3,bg="darkgrey"))
  
  with(mod, lines(phi~l,type="l",lwd=3,col=2))
  
  mtext(par_sel, 3, line=-0.2,cex=2,col="darkgrey")
  
  mtext("Equal catch sharing", 3, line=-17.5,cex=1,col="darkgreen",adj=0)
  
  mtext(paste("model:", fun_sel,sep=""),3,line=0,cex=1.5,adj=0.05,padj=2,col="darkgrey",outer=T)

  
}}










