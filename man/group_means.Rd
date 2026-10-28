\name{group_means}
\alias{group_means}
\alias{print.goric_group_means}
\title{Group means in terms of the parameters of an lm object (for goric hypotheses)}
\description{
Expresses the group means (cell means) of a fitted \code{lm} (or \code{aov},
\code{glm}) object as linear combinations of its regression coefficients and
gives them a short notation (e.g., \code{mu_12}). These names can be used
directly in the \code{hypotheses} of \code{\link{goric}}: the definitions are
added automatically. The function can also be used to obtain the syntax
(\code{mu_12 := ...}) which can be copied into the hypothesis text.

This is especially useful when the model contains an intercept and/or more
than one grouping variable (factor), since the group means are then not
parameters of the model themselves.
}
\usage{
group_means(object, labels = c("number", "abbrev", "level"),
            prefix = "mu", sep = "_", digits = 10)

\method{print}{goric_group_means}(x, \dots)
}
\arguments{
  \item{object}{a fitted object of class \code{lm}, \code{mlm}, \code{aov}
    or \code{glm}.}
  \item{labels}{notation of the subscripts (one per grouping variable):
    \code{"number"} (default) uses the level number of each grouping
    variable (e.g., \code{mu_12} for level 1 of the first and level 2 of the
    second grouping variable; if a grouping variable has 10 or more levels,
    the subscripts are separated by \code{sep}, e.g., \code{mu_1_12});
    \code{"abbrev"} uses the shortest unique abbreviation of the level names
    (e.g., \code{mu_O_0.5}), as obtained by
    \code{abbreviate(levels(<factor>), minlength = 1)} (which may remove
    lower-case vowels; e.g., 'Melon' may become 'Ml'); levels that are numbers
    (e.g., 0.5) are used as is; use \code{group_means(model, labels = "abbrev")}
    to see the resulting names; \code{"level"} uses the full level names
    (e.g., \code{mu_OJ_0.5}). In \code{goric}, all three notations can be
    used.}
  \item{prefix}{character; the name of the group means. Default \code{"mu"}.}
  \item{sep}{character; separator used in the names. Default \code{"_"}.}
  \item{digits}{number of digits used to round the weights of the
    coefficients.}
  \item{x}{an object of class \code{goric_group_means}.}
  \item{\dots}{not used.}
}
\details{
Grouping variables (factors) are: factors (incl. ordered factors), character
and logical variables, and numeric variables taking on exactly two distinct
values (e.g., 0/1 dummy variables). All combinations of the levels of these
variables are formed. Other (continuous) predictors are not used to form the
groups: they are set to 0, so the group means are then the adjusted group
means at a value of 0 for the continuous predictors (center these predictors
to obtain the adjusted means at their average value).

The group means are computed from the model matrix, so they are correct for
any coding (contrasts), with or without intercept and with or without
interactions. For a main-effects model (no interaction), the group means are
the means implied by that model (not the observed cell means).

The coefficient names in the definitions follow the restriktor notation,
that is, \code{(Intercept)} becomes \code{.Intercept.} and \code{:} becomes
\code{.}.

For an \code{mlm} (multivariate lm, e.g., \code{lm(cbind(y1, y2) ~ g)}), the
group means are given per response variable and named like the parameters,
that is, \code{<response>.<name>} (e.g., \code{y1.mu_12}, cf. the parameter
\code{y1.g2}).

For a \code{glm} with a non-identity link, the "group means" are on the
scale of the linear predictor (e.g., log-odds in a logistic regression).

In \code{goric}, for an \code{lm} object, each name of a group mean that is
used in a (character) hypothesis gets its definition added automatically,
unless the user defines that name him/herself (via \code{:=}) or it is a
parameter name of the model.
}
\value{
An object of class \code{goric_group_means}, a list with:
  \item{table}{data frame with the name of each group mean, (for an mlm) the
    response variable, the levels of the grouping variables, the definition in
    terms of the coefficients, and the number of observations (\code{n}).}
  \item{syntax}{named character vector with the definitions in the form
    \code{"mu_12 := .Intercept. + dose1"}.}
  \item{L}{matrix with the weights of the coefficients for each group mean
    (group means = \code{L} times the coefficients, for an mlm stacked per
    response variable).}
  \item{group_vars}{the names of the grouping variables.}
  \item{continuous_vars}{the names of the other (continuous) predictors.}
  \item{intercept}{logical; does the model contain an intercept.}
}
\author{Rebecca M. Kuiper}
\seealso{\code{\link{goric}}}
\examples{
ToothGrowth2 <- ToothGrowth
ToothGrowth2$dose <- factor(ToothGrowth2$dose)

# model with intercept and two (interacting) factors
fit <- lm(len ~ supp * dose, data = ToothGrowth2)

# notation and definitions of the 2 x 3 group means
group_means(fit)

# use the group means directly in the hypotheses
H1 <- "mu_11 < mu_12 < mu_13; mu_21 < mu_22 < mu_23"
goric(fit, hypotheses = list(H1 = H1))

# the same, but with the level names as notation
group_means(fit, labels = "level")
goric(fit, hypotheses = list(H1 = "mu_OJ_0.5 < mu_OJ_1 < mu_OJ_2"))

# multivariate (mlm): group means per response variable
ToothGrowth2$len2 <- ToothGrowth2$len + rnorm(nrow(ToothGrowth2))
fit_mlm <- lm(cbind(len, len2) ~ supp * dose, data = ToothGrowth2)
group_means(fit_mlm)
goric(fit_mlm, hypotheses = list(H1 = "len.mu_11 < len.mu_12; len2.mu_11 < len2.mu_12"))
}
