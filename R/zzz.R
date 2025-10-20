.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("This is ",paste(pkgname, version))
    packageStartupMessage("Please report any bugs to info@restriktor.org")
}

utils::globalVariables(c("green", "reset", "Value", 
                         "Group_hypo_comparison", "Group_pop_values", 
                         "sample_value", "percentile_value",
                         "percentile_label", "Group"
                         ))