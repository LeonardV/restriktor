\name{Burns}
\alias{Burns}
\docType{data}
\title{ 
Relation between the response variable PTSS and gender, age, TBSA, guilt and anger.
}
\description{
Simulated dataset based on the original model parameters. The original 
data are based on two cohort studies in children from 0 to 4 and 8 
to 18 years old with burns and their mother. 
}

\usage{data(Burns)}
\format{
  A data frame of 278 observations of 4 variables.
  \describe{
    \item{\code{PTSS}}{post-traumatic stress symptoms}
    \item{\code{gender}}{gender}
    \item{\code{age}}{age in years}
    \item{\code{TBSA}}{estimated percentage total body surface area affected by second and third degree burns}
    \item{\code{guilt}}{parental guilt feelings in relation to the burn event}
    \item{\code{anger}}{parental anger feelings in relation to the burn event}
  }
}
\references{
Bakker A, Van der Heijden PG, Van Son MJ, Van Loey NE. Course of traumatic stress reactions in couples after a burn event to their young child. Health Psychology 2013; 10(32):1076-1083, doi:10.1037/a0033983.

Egberts MR, van de Schoot R, Boekelaar A, Hendrickx H, Geenen R, NEE V. Child and adolescent internalizing and externalizing problems 12 months postburn: the potential role of preburn functioning, parental posttraumatic stress, and informant bias. Child \& Adolescent Psychiatry 2016; 25:791-803.
}
\examples{
head(Burns)
}
