restriktor is a free, open source R package for linear equality and inequality 
constrained statistical estimation, inference and evaluation for linear models. More 
information can be found online at http://restriktor.org.


##### Install R #####
Restriktor is implemented as an R package. This means that before installing
restriktor, you should have installed a recent version (>= 3.0.0) of R. You can
download the latest version of R from the
[R-project website](http://www.r-project.org).

##### Install Graphical User Interface (GUI) #####
R is a command line driven program. This means that it does not have a graphical
user interface (GUI). Luckily, there are many good GUI's to make life easier, for
example [Rstudio](http://rstudio.org), [R Commander](http://www.rcommander.com/)
and [RKWard](https://rkward.kde.org/).

##### Install restriktor #####
Once you have installed R, the next step is to install restriktor. This can be
done by typing in R:

`install.packages("restriktor", dependencies = TRUE)`

To check if the installation was successful, you can load the restriktor package
and try for example:

`library(restriktor)`

`# construct constraint syntax based on the factor level names` 

`constraints <- 'GroupActive  < GroupPassive
                GroupPassive < GroupControl
                GroupControl < GroupNo'
`

Fit the unrestricted linear model, where "Age" is the response
variable and `"Group"` a factor with four treatment groups.

`fit.ANOVA <- lm(Age ~ -1 + Group, data = ZelazoKolb1972)`

`# fit the restricted model`
`restr.ANOVA <- restriktor(fit.ANOVA, constraints = constraints)
`

`# summary of the restricted parameter estimates`
`summary(restr.ANOVA)
`

If you can see the output, everything is set up and ready.

#### Getting Started #####
Now that you have installed restriktor, 
