library(restriktor)
fit_ANOVA <- lm(Anger ~ -1 + Group, data = AngerManagement)

myConstraints_cin <- ' GroupNo < GroupPhysical
                       GroupPhysical < GroupBehavioral
                       GroupBehavioral == GroupBoth '

restr_ANOVA <- restriktor(fit_ANOVA, constraints = myConstraints_cin)

out <- iht(restr_ANOVA, boot = "parametric", R=9)
out
out <- iht(restr_ANOVA, boot = "model.based", R=9)
out


myConstraints_ceq <- ' GroupNo == GroupPhysical
                       GroupPhysical == GroupBehavioral
                       GroupBehavioral == GroupBoth '

restr_ANOVA <- restriktor(fit_ANOVA, constraints = myConstraints_ceq)

out <- iht(restr_ANOVA, boot = "parametric", R=9)
out
out <- iht(restr_ANOVA, boot = "model.based", R=9)
out
