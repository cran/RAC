% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clam_ind_dataloader.R
\name{Clam_ind_dataloader}
\alias{Clam_ind_dataloader}
\title{Function that loads forcings data for Clam individual bioenergetic model and performs the interpolation}
\usage{
Clam_ind_dataloader(userpath)
}
\arguments{
\item{userpath}{the path where folder containing model inputs and outputs is located}
}
\value{
a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l], particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
}
\description{
Function that loads forcings data for Clam individual bioenergetic model and performs the interpolation
}
