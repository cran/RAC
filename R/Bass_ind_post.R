#' Postprocess the Bass indivual bioenergetic model results
#'
#' @param userpath the path where the working folder is located
#' @param output the output list of the
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#'
#' @return a list containing the fish weight, proteines, lipids and carbohydrates wasted or produced with excretions, potential and actual ingestion rates, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Bass_ind_post<-function(userpath,output,times,Dates) {

cat('Data post-processing\n')
cat('\n')

ti=times[1]           # Integration beginning
tf=times[2]           # Integration end

# Extracts outputs from the output list
weight=unlist(output[1])
exc=output[[2]]
wst=output[[3]]
ing=unlist(output[4])
ingvero=unlist(output[5])
Tfun=output[[6]]
metab=output[[7]]

# Adjusts results acoording with integration extremes
# now day 1 coincides with ti
weightSave=weight[(ti+1):tf]
excSave=exc[(ti+1):tf,]
wstSave=wst[(ti+1):tf,]
ingSave=ing[(ti+1):tf]
ingveroSave=ingvero[(ti+1):tf]
TfunSave=Tfun[(ti+1):tf,]
metabSave=metab[(ti+1):tf,]
output=list(weightSave,excSave,ingSave,ingveroSave,wstSave,metabSave,TfunSave)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti) # create a dates vector to plot results

# Plot weight
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//weight.jpeg")
jpeg(filepath,800,600)
plot(days,weightSave,ylab="Weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot excretion
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//excretion.jpeg")
jpeg(filepath,800,600)
ub=max(max(excSave[,1]),max(excSave[,2]),max(excSave[,3]))
plot(days,excSave[,1],ylab="Excreted quantities (g/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,excSave[,2],col="blue")
lines(days,excSave[,3],col="black")
legend("topleft",c("Excreted Proteins","Excreted Lipids","excreted Carbohydrates"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot wasted food
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//waste.jpeg")
jpeg(filepath,800,600)
ub=max(max(wstSave[,1]),max(wstSave[,2]),max(wstSave[,3]))
plot(days,wstSave[,1],ylab="Quantities to waste (g/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,wstSave[,2],col="blue")
lines(days,wstSave[,3],col="black")
legend("topleft",c("Proteins to waste","Lipids to waste","Carbohydrates to waste"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot ingested food
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//ingestion.jpeg")
jpeg(filepath,800,600)
plot(days,ingveroSave,ylab="Ingested food (g)",xlab=" ",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//T_limitation.jpeg")
jpeg(filepath,800,600)
ub=max(max(TfunSave[,1]),max(TfunSave[,2]))
plot(days,TfunSave[,1],ylab="Temperature limitation functions",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,TfunSave[,2],col="blue")
legend("topright",c("Anabolism limitation","Catabolism limitation"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//metabolism.jpeg")
jpeg(filepath,800,600)
ub=max(max(metabSave[,1]),max(metabSave[,2]))
plot(days,metabSave[,1],ylab="Metabolic rate (J/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,metabSave[,2],col="blue")
legend("topright",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//weight.csv")
write.csv(weightSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//excretion.csv")
write.csv(excSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//waste.csv")
write.csv(wstSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//potential_ingestion.csv")
write.csv(ingSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//actual_ingestion.csv")
write.csv(ingveroSave,filepath)

return(output)
}

