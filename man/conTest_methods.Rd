\name{conTest-methods}
\alias{conTest-methods}
\alias{iht-methods}
\alias{print.conTest}

\title{Methods for conTest}
\description{The print() function prints the output of objects of class \code{conTest}. 
}

\usage{

\method{print}{conTest}(x, digits = max(3, getOption("digits") - 2), \dots)

}

\arguments{
  \item{x}{an object of class \code{conTest}.} 
  
  \item{digits}{the number of significant digits to use when printing.}
  
  \item{\dots}{no additional arguments for now.}
}
  
  
  \examples{
  # unrestricted linear model for ages (in months) at which an 
  # infant starts to walk alone.
  
  # prepare data
  DATA <- subset(ZelazoKolb1972, Group != "Control")
  
  # fit unrestricted linear model
  fit_lm <- lm(Age ~ -1 + Group, data = DATA)
  
  # restricted linear model with restrictions that the walking 
  # exercises would not have a negative effect of increasing the 
  # mean age at which a child starts to walk. 
  fit_con <- restriktor(fit_lm, constraints = "GroupActive  < GroupPassive; 
                                               GroupPassive < GroupNo")
  
  fit_con  
}
  