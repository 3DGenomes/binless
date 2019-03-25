.onAttach <- function(libname, pkgname) {
   packageStartupMessage(paste0("Binless version ",as.character(packageVersion("binless")),"\n"))
}

