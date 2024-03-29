% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bass_ind_RKsolver.R
\name{Bass_ind_RKsolver}
\alias{Bass_ind_RKsolver}
\title{Solves the Seabass bioenergetic balance with a 4th order Runge Kutta method}
\usage{
Bass_ind_RKsolver(Param, Temperature, G, Food, IC, times)
}
\arguments{
\item{Param}{vector containing all metabolic parameters}

\item{Temperature}{water temperature forcing time series}

\item{G}{food entering the cage time series}

\item{Food}{food characterization (Proteins, Lipids, Carbohydrates)}

\item{IC}{initial conditions}

\item{times}{vector containing integration extremes and integration timestep}
}
\value{
a list containing the fish weight, proteines, lipids and carbohydrates wasted or produced with excretions, potential and actual ingestion rates, temperature limitation functions and metabolic rates
}
\description{
Solves the Seabass bioenergetic balance with a 4th order Runge Kutta method
}
