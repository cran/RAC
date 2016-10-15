#' Postprocess the Clam population bioenergetic model results
#'
#' @param userpath the path where the working folder is located
#' @param output the output list of the
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param N the number of individuals
#'
#' @return  a list containing the clam weights, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Clam_pop_post<-function(userpath,output,times,Dates,N) {

  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end

  # Extracts results from output list
  Wd_stat=output[[1]]
  Ww_stat=output[[2]]
  L_stat=output[[3]]
  A_stat=output[[4]]
  C_stat=output[[5]]
  fgT=output[[6]]
  frT=output[[7]]

cat('Data post-processing\n')
cat('\n')

# Adjusts results acoording with integration extremes
# now day 1 coincides with ti

WdSave=Wd_stat[,ti:tf]
WwSave=Ww_stat[,ti:tf]
LSave=L_stat[,ti:tf]

ASave=A_stat[,ti:tf]
CSave=C_stat[,ti:tf]

fgT=fgT[ti:tf]
frT=frT[ti:tf]

N=N[ti:tf]

output=list(WdSave,WwSave,LSave,ASave,CSave,fgT,frT,N)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1) # create a dates vector to plot results
currentpath=getwd()

# Plot weight
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//Dryweight.jpeg")
jpeg(filepath,800,600)
ub=WdSave[1,]+WdSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(WdSave[1,]-WdSave[2,])){
lb[i]=max(WdSave[1,i]-WdSave[2,i],0)
}
maxub=max(WdSave[1,]+WdSave[2,])
plot(days,WdSave[1,],ylab="Weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,WdSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot length
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//Length.jpeg")
jpeg(filepath,800,600)
ub=LSave[1,]+LSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(LSave[1,]-LSave[2,])){
  lb[i]=max(LSave[1,i]-LSave[2,i],0)
}
maxub=max(LSave[1,]+LSave[2,])
plot(days,LSave[1,],ylab="Length (mm)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,LSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//Tfun.jpeg")
jpeg(filepath,800,600)
ub=max(max(fgT),max(frT))
plot(days,fgT,ylab="Temperature limitation functions",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,frT,col="blue")
legend("topright",c("Anabolism limitation","Catabolism limitation"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//metabolism.jpeg")
jpeg(filepath,800,600)
Aub=ASave[1,]+ASave[2,]
Cub=CSave[1,]+CSave[2,]
Alb=as.matrix(matrix(0,nrow=length(Aub),ncol=1))
Clb=as.matrix(matrix(0,nrow=length(Cub),ncol=1))
for (i in 1:length(ASave[1,]-ASave[2,])){
  Alb[i]=max(ASave[1,i]-ASave[2,i],0)
  Clb[i]=max(CSave[1,i]-CSave[2,i],0)
}
maxub=max(Aub,Cub)
plot(days,ASave[1,],ylab="Metabolic rates (J/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(Alb,rev(Aub)),col="grey75",border=FALSE)
lines(days,ASave[1,],lwd=2,col="red")
polygon(c(days,rev(days)),c(Clb,rev(Cub)),col="grey75",border=FALSE)
lines(days,CSave[1,],lwd=2,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
legend("topleft",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
dev.off()

# plot population dynamics
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//Population.jpeg")
jpeg(filepath,800,600)
plot(days, N, ylab="Individuals", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//Dryweight.csv")
write.csv(t(WdSave),filepath)

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//Wetweight.csv")
write.csv(t(WwSave),filepath)

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//Length.csv")
write.csv(t(LSave),filepath)

return(output)

}
