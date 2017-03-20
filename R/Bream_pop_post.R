#' Postprocess the Bream population bioenergetic model results
#'
#' @param userpath the path where the working folder is located
#' @param output output list containing the output of the RK solver
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param N the number of individuals
#' @param CS the commercial size of Seabream
#'
#' @return a list containing the fish weight, proteines, lipids and carbohydrates wasted or produced with excretions, potential and actual ingestion rates, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Bream_pop_post<-function(userpath,output,times,Dates,N,CS) {

cat('Data post-processing\n')
cat('\n')

ti=times[1]           # Integration beginning
tf=times[2]           # Integration end

# Extracts outputs from the output list
W_stat=output[[1]]
Pexc_stat=output[[2]]
Lexc_stat=output[[3]]
Cexc_stat=output[[4]]
ingestion_stat=output[[5]]
Pwst_stat=output[[6]]
Lwst_stat=output[[7]]
Cwst_stat=output[[8]]
A_stat=output[[9]]
C_stat=output[[10]]
fgT=output[[11]]
frT=output[[12]]


# Adjusts results acoording with integration extremes
# now day 1 coincides with ti
weightSave=W_stat[ti:tf,]

PexcSave=Pexc_stat[ti:tf,]
LexcSave=Lexc_stat[ti:tf,]
CexcSave=Cexc_stat[ti:tf,]

ingestionSave=ingestion_stat[ti:tf,]

PwstSave=Pwst_stat[ti:tf,]
LwstSave=Lwst_stat[ti:tf,]
CwstSave=Cwst_stat[ti:tf,]

ASave=A_stat[ti:tf,]
CSave=C_stat[ti:tf,]

fgT=fgT[ti:tf]
frT=frT[ti:tf]

N=N[ti:tf]

# Days to commercial size

# Lower bound
foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(weightSave[,1]-weightSave[,2])
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  Lb_daysToSize="Not reaching the commercial size"
}else{  Lb_daysToSize <- min(NonNAindex)
}

# Mean
foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(weightSave[,1])
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  Mean_daysToSize="Not reaching the commercial size"
}else{  Mean_daysToSize <- min(NonNAindex)
}

# Upper bound
foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(weightSave[,1]+weightSave[,2])
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  Ub_daysToSize="Not reaching the commercial size"
}else{  Ub_daysToSize <- min(NonNAindex)
}

# List containing days to size
daysToSize<-as.list(cbind(Ub_daysToSize,Mean_daysToSize,Lb_daysToSize))

output=list(weightSave,PexcSave,LexcSave,CexcSave,ingestionSave,PwstSave,LwstSave,CwstSave,ASave,CSave,fgT,frT,N, daysToSize)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1) # create a dates vector to plot results
currentpath=getwd()

# Plot weight
filepath=paste0(userpath,"/Bream_population/Outputs/Out_plots//weight.jpeg")
jpeg(filepath,800,600)
ub=weightSave[,1]+weightSave[,2]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(weightSave[,1]-weightSave[,2])){
lb[i]=max(weightSave[i,1]-weightSave[i,2],0)
}
maxub=max(weightSave[,1]+weightSave[,2])
plot(days,weightSave[,1],ylab="Weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,weightSave[,1],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot excretion
filepath=paste0(userpath,"/Bream_population/Outputs/Out_plots//excretion.jpeg")
jpeg(filepath,800,600)
Lub=LexcSave[,1]+LexcSave[,2]
Pub=PexcSave[,1]+PexcSave[,2]
Cub=CexcSave[,1]+CexcSave[,2]
Llb=as.matrix(matrix(0,nrow=length(Lub),ncol=1))
Plb=as.matrix(matrix(0,nrow=length(Pub),ncol=1))
Clb=as.matrix(matrix(0,nrow=length(Cub),ncol=1))
for (i in 1:length(LexcSave[,1]-LexcSave[,2])){
  Llb[i]=max(LexcSave[i,1]-LexcSave[i,2],0)
  Plb[i]=max(PexcSave[i,1]-PexcSave[i,2],0)
  Clb[i]=max(CexcSave[i,1]-CexcSave[i,2],0)
}
maxub=max(Lub,Pub,Cub)
plot(days,LexcSave[,1],ylab="Excreted quantities (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(Llb,rev(Lub)),col="grey75",border=FALSE)
lines(days,LexcSave[,1],lwd=2,col="red")
polygon(c(days,rev(days)),c(Plb,rev(Pub)),col="grey75",border=FALSE)
lines(days,PexcSave[,1],lwd=2,col="green")
polygon(c(days,rev(days)),c(Clb,rev(Cub)),col="grey75",border=FALSE)
lines(days,CexcSave[,1],lwd=2,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
legend("topleft",c("Excreted Proteins","Excreted Lipids","excreted Carbohydrates"),fill=c("red","green","blue"))
dev.off()


# plot wasted food
filepath=paste0(userpath,"/Bream_population/Outputs/Out_plots//waste.jpeg")
jpeg(filepath,800,600)
Lub=LwstSave[,1]+LwstSave[,2]
Pub=PwstSave[,1]+PwstSave[,2]
Cub=CwstSave[,1]+CwstSave[,2]
Llb=as.matrix(matrix(0,nrow=length(Lub),ncol=1))
Plb=as.matrix(matrix(0,nrow=length(Pub),ncol=1))
Clb=as.matrix(matrix(0,nrow=length(Cub),ncol=1))
for (i in 1:length(LwstSave[,1]-LwstSave[,2])){
  Llb[i]=max(LwstSave[i,1]-LwstSave[i,2],0)
  Plb[i]=max(PwstSave[i,1]-PwstSave[i,2],0)
  Clb[i]=max(CwstSave[i,1]-CwstSave[i,2],0)
}
maxub=max(Lub,Pub,Cub)
plot(days,LwstSave[,1],ylab="Excreted quantities (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(Llb,rev(Lub)),col="grey75",border=FALSE)
lines(days,LwstSave[,1],lwd=2,col="red")
polygon(c(days,rev(days)),c(Plb,rev(Pub)),col="grey75",border=FALSE)
lines(days,PwstSave[,1],lwd=2,col="green")
polygon(c(days,rev(days)),c(Clb,rev(Cub)),col="grey75",border=FALSE)
lines(days,CwstSave[,1],lwd=2,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
legend("topleft",c("Proteins to waste","Lipids to waste","Carbohydrates to waste"),fill=c("red","green","blue"))
dev.off()

# plot ingested food
filepath=paste0(userpath,"/Bream_population/Outputs/Out_plots//ingestion.jpeg")
jpeg(filepath,800,600)
ub=ingestionSave[,1]+ingestionSave[,2]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(ingestionSave[,1]-ingestionSave[,2])){
  lb[i]=max(ingestionSave[i,1]-ingestionSave[i,2],0)
}
maxub=max(ingestionSave[,1]+ingestionSave[,2])
plot(days,ingestionSave[,1],ylab="Ingested food (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,ingestionSave[,1],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
filepath=paste0(userpath,"/Bream_population/Outputs/Out_plots//Tfun.jpeg")
jpeg(filepath,800,600)
ub=max(max(fgT),max(frT))
plot(days,fgT,ylab="Temperature limitation functions",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,frT,col="blue")
legend("topright",c("Anabolism limitation","Catabolism limitation"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/Bream_population/Outputs/Out_plots//metabolism.jpeg")
jpeg(filepath,800,600)
Aub=ASave[,1]+ASave[,2]
Cub=CSave[,1]+CSave[,2]
Alb=as.matrix(matrix(0,nrow=length(Aub),ncol=1))
Clb=as.matrix(matrix(0,nrow=length(Cub),ncol=1))
for (i in 1:length(LwstSave[,1]-LwstSave[,2])){
  Alb[i]=max(ASave[i,1]-ASave[i,2],0)
  Clb[i]=max(CSave[i,1]-CSave[i,2],0)
}
maxub=max(Aub,Cub)
plot(days,ASave[,1],ylab="Metabolic rates (J/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(Alb,rev(Aub)),col="grey75",border=FALSE)
lines(days,ASave[,1],lwd=2,col="red")
polygon(c(days,rev(days)),c(Clb,rev(Cub)),col="grey75",border=FALSE)
lines(days,CSave[,1],lwd=2,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
legend("topleft",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
dev.off()

# plot population dynamics
filepath=paste0(userpath,"/Bream_population/Outputs/Out_plots//Population.jpeg")
jpeg(filepath,800,600)
plot(days, N, ylab="Individuals", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//weight.csv")
write.csv(weightSave,filepath)

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//Pexc.csv")
write.csv(PexcSave,filepath)

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//Lexc.csv")
write.csv(LexcSave,filepath)

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//Cexc.csv")
write.csv(CexcSave,filepath)

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//Pwst.csv")
write.csv(PwstSave,filepath)

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//Lwst.csv")
write.csv(LwstSave,filepath)

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//Cwst.csv")
write.csv(CwstSave,filepath)

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//ing.csv")
write.csv(ingestionSave,filepath)

filepath=paste0(userpath,"/Bream_population/Outputs/Out_csv//Days_to_comercial_size.csv")
write.csv(daysToSize,filepath)

return(output)

}
