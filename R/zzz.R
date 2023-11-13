.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("This is ",paste(pkgname, version))
    #packageStartupMessage(pkgname, " is BETA software! Please report any bugs.")
    packageStartupMessage("Please report any bugs.")
}




# restriktorStartupMessage <- function() {
#   msg <- c(paste0("
#   
# ██████╗ ███████╗███████╗████████╗██████╗ ██╗██╗  ██╗████████╗ ██████╗ ██████╗ 
# ██╔══██╗██╔════╝██╔════╝╚══██╔══╝██╔══██╗██║██║ ██╔╝╚══██╔══╝██╔═══██╗██╔══██╗
# ██████╔╝█████╗  ███████╗   ██║   ██████╔╝██║█████╔╝    ██║   ██║   ██║██████╔╝
# ██╔══██╗██╔══╝  ╚════██║   ██║   ██╔══██╗██║██╔═██╗    ██║   ██║   ██║██╔══██╗
# ██║  ██║███████╗███████║   ██║   ██║  ██║██║██║  ██╗   ██║   ╚██████╔╝██║  ██║
# ╚═╝  ╚═╝╚══════╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝ 
# version ", packageVersion("restriktor"), " (restriktor.org)"))
#   return(msg)
# }
  

# .onAttach <- function(lib, pkg) {
#   # startup message
#   msg <- restriktorStartupMessage()
#   if(!interactive())
#     msg[1] <- paste("Package 'restriktor' version", packageVersion("restriktor"))
#   packageStartupMessage(msg)      
#   invisible()
# }