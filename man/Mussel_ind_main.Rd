% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mussel_ind_main.R
\name{Mussel_ind_main}
\alias{Mussel_ind_main}
\title{Mussel bioenergetic individual model}
\usage{
Mussel_ind_main(userpath, forcings)
}
\arguments{
\item{userpath}{the path where the working folder is located}

\item{forcings}{a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l] and its characterization in terms of C/P and N/P molar ratios, particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]}
}
\value{
A list containing model outputs: weight, length mussel CNP, pseudofaeces CNP production, faeces CNP production, temperature limitation functions, metabolic rates and oxygen consumption
}
\description{
Solves the bioenergetic balance for Mussel
}
