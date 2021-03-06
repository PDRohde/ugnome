% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diseaseIDs.R
\name{getID}
\alias{getID}
\title{Identify disease cases within the UK Biobank.}
\usage{
getID(ICD = NULL, SR = NULL)
}
\arguments{
\item{ICD}{a vector of ICD-10 codes}

\item{SR}{a vector of self-reported disease codes.}
}
\description{
In the UK Biobank (UKBB) diseases are reported based on either in-hospital records (ICD-10 classification system; field id 41270),
or based on self-reported diseases (field id 20001 or 20002).

The reported diseases are listed in large matrices where many of the cells are filled with NA's.
This utility function enables easy extraction of disease cases based on a vector of ICD-10 codes and/or a vector of self-reported diseases.

This function requires two matrices 'self' and 'icd10', which are matrices with self-reported and ICD-10 diagnoses,
where each row correspond to one individual and the columns correspond to the array number (see UKBB showcase for details).
}
\author{
Palle Duun Rohde
}
