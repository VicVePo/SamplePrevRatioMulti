% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{logistica_estudio_transversal}
\alias{logistica_estudio_transversal}
\title{Logística del Estudio Transversal}
\usage{
logistica_estudio_transversal(
  n_final,
  tasa_rechazo,
  tasa_elegibilidad,
  sujetos_por_dia,
  dias_laborables_mes
)
}
\arguments{
\item{n_final}{Tamaño muestral final calculado}

\item{tasa_rechazo}{Tasa de rechazo esperada (proporción entre 0 y 1)}

\item{tasa_elegibilidad}{Tasa de elegibilidad esperada (proporción entre 0 y 1)}

\item{sujetos_por_dia}{Número de sujetos/registros que pueden ser procesados por día}

\item{dias_laborables_mes}{Número de días laborables por mes}
}
\value{
Una lista con los cálculos logísticos del estudio
}
\description{
Logística del Estudio Transversal
}
