.First.lib <- function(lib, pkg)
  library.dynam("sjedrp", pkg, lib)

.Last.lib <- function (libpath) {
  ## Run when the package is being unloaded.  This allows us to test
  ## packages within same session when dynlib is updated.
  library.dynam.unload("sjedrp", libpath)
}
