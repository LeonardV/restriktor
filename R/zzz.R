.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("This is ",paste(pkgname, version))
    #packageStartupMessage(pkgname, " is BETA software! Please report any bugs.")
    packageStartupMessage("Please report any bugs to info@restriktor.org")
}

utils::globalVariables(c("Group", "Value", "green", "reset", 
                         "Group_hypo_comparison", "Group_pop_values", 
                         "sample_value", "x", "y"))