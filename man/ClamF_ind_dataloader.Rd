% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClamF_ind_dataloader.R
\name{ClamF_ind_dataloader}
\alias{ClamF_ind_dataloader}
\title{Function that loads forcings data for Clam individual bioenergetic model (alternative version) and performs the interpolation}
\usage{
ClamF_ind_dataloader(userpath)
}
\arguments{
\item{userpath}{the path where folder containing model inputs and outputs is located}
}
\value{
a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3]
}
\description{
Function that loads forcings data for Clam individual bioenergetic model (alternative version) and performs the interpolation
}
