#' Postprocess the Mussel indivual bioenergetic model results
#'
#' @param userpath the path where the working folder is located
#' @param output output list containing the output of the RK solver
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param CS the commercial size of Mussel
#'
#' @return a list containing the weights of the mussel, the excreted CNP, the mussel CNP, temperature limitation functions, metabolic rates, oxygen consumption
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Mussel_ind_post<-function(userpath,output,times,Dates,CS) {

cat('Data post-processing\n')
cat('\n')

ti=times[1]           # Integration beginning
tf=times[2]           # Integration end

# Extracts outputs from the output list
W=output[[1]]
pfec=output[[2]]
fec=output[[3]]
comp=output[[4]]
tfun=output[[5]]
metab=output[[6]]
cons=output[[7]]

# Adjusts results acoording with integration extremes
# now day 1 coincides with ti
weightSave=t(W[,ti:tf])
pfecSave=pfec[ti:tf,]
fecSave=fec[ti:tf,]
compSave=comp[ti:tf,]
tfunSave=tfun[ti:tf,]
metabSave=metab[ti:tf,]
consSave=cons[ti:tf]

foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(weightSave)
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  daysToSize="Not reaching the commercial size"
}else{  daysToSize <- min(NonNAindex)
}
daysToSize<-as.list(daysToSize)

output=list(weightSave,fecSave,compSave,tfunSave,metabSave,consSave,daysToSize)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1) # create a dates vector to plot results
currentpath=getwd()

# Plot weight
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//dryweight.jpeg")
jpeg(filepath,800,600)
plot(days,weightSave[,1],ylab="Dry weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red")
lines(days,weightSave[,2],col="blue")
lines(days,weightSave[,3],col="black")
legend("topleft",c("Somatic tissue","Gonadic tissue","Total"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot length
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//length.jpeg")
jpeg(filepath,800,600)
plot(days,weightSave[,5],ylab="Length (cm)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot pseudofecies
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//pseudofaeces.jpeg")
jpeg(filepath,800,600)
ub=max(max(pfecSave[,1]),max(pfecSave[,2]),max(pfecSave[,3]))
plot(days,pfecSave[,1],ylab="pseudofaeces production (g/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,pfecSave[,2],col="blue")
lines(days,pfecSave[,3],col="black")
legend("topleft",c("Excreted C","Excreted N","excreted P"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot fecies
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//faeces.jpeg")
jpeg(filepath,800,600)
ub=max(max(fecSave[,1]),max(fecSave[,2]),max(fecSave[,3]))
plot(days,fecSave[,1],ylab="faeces production (g/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,fecSave[,2],col="blue")
lines(days,fecSave[,3],col="black")
legend("topleft",c("Excreted C","Excreted N","excreted P"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot CNP mytilus
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//composition.jpeg")
jpeg(filepath,800,600)
ub=max(max(compSave[,1]),max(compSave[,2]),max(compSave[,3]))
plot(days,compSave[,1],ylab="CNP mytilus (g)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,compSave[,2],col="blue")
lines(days,compSave[,3],col="black")
legend("topleft",c("C","N","P"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//T_limitation.jpeg")
jpeg(filepath,800,600)
ub=max(max(tfunSave[,1]),max(tfunSave[,2]))
plot(days,tfunSave[,1],ylab="Temperature limitation functions",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,tfunSave[,2],col="blue")
legend("topright",c("Anabolism limitation","Catabolism limitation"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//metabolism.jpeg")
jpeg(filepath,800,600)
ub=max(max(metabSave[,1]),max(metabSave[,2]))
plot(days,metabSave[,1],ylab="Metabolic rate (J/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,metabSave[,2],col="blue")
legend("topright",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot O2 consumption
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//O2consumption.jpeg")
jpeg(filepath,800,600)
plot(days,consSave,ylab="O2 consumption (g/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//weight.csv")
write.csv(weightSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//faeces.csv")
write.csv(fecSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//pseudofaeces.csv")
write.csv(pfecSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//CNPcontent.csv")
write.csv(compSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//O2consumption.csv")
write.csv(consSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//Days_to_comercial_size.csv")
write.csv(daysToSize,filepath)

return(output)

}
