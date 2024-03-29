% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClamF_ind_main.R
\name{ClamF_ind_main}
\alias{ClamF_ind_main}
\title{Clam bioenergetic individual model (alternative version)}
\usage{
ClamF_ind_main(userpath, forcings)
}
\arguments{
\item{userpath}{the path where the working folder is located}

\item{forcings}{a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3]}
}
\value{
A list containing model outputs: weights, temperature limitation functions and metabolic rates
}
\description{
Clam bioenergetic individual model (alternative version)
}
