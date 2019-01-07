# TODO: remove?
# .onLoad <- function(libname, pkgname) 
# {
#   library.dynam("ppgmmga", pkgname, libname)
# }
# 
# .onUnload <- function (lib)
# {
#   library.dynam.unload("ppgmmga", lib)
# }

ppgmmgaStartupMessage <- function()
{

  # Startup message obtained as 
# > figlet -f slant ppgmmga
  msg <- c(paste0(
"    ____  ____  ____ _____ ___  ____ ___  ____ _____ _
   / __ \\/ __ \\/ __ `/ __ `__ \\/ __ `__ \\/ __ `/ __ `/
  / /_/ / /_/ / /_/ / / / / / / / / / / / /_/ / /_/ / 
 / .___/ .___/\\__, /_/ /_/ /_/_/ /_/ /_/\\__, /\\__,_/  
/_/   /_/    /____/                    /____/  version ", packageVersion("ppgmmga")))

# Startup message obtained as 
# > figlet -f smslant ppgmmga
  msg <- c(paste0(
"   ___  ___  ___ ___ _  __ _  ___ ____ _
  / _ \\/ _ \\/ _ `/  ' \\/  ' \\/ _ `/ _ `/
 / .__/ .__/\\_, /_/_/_/_/_/_/\\_, /\\_,_/ 
/_/  /_/   /___/            /___/       version ", 
packageVersion("ppgmmga")))
    
# "\nType 'citation(\"ppgmmga\")' for citing this R package in publications."
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # unlock .ppgmmga variable allowing its modification
  unlockBinding(".ppgmmga", asNamespace("ppgmmga"))
  # startup message
  msg <- ppgmmgaStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'ppgmmga' version", packageVersion("ppgmmga"))
  packageStartupMessage(msg)
  invisible()
}


