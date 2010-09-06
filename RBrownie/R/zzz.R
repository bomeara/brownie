
# use this with a namespace
.onLoad <- function(lib, pkg) {
    require(phylobase)
    require(methods) # not sure if this is needed
    phylobase.options(singleton="ok")
}
