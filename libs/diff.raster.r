diff.raster <- function(r) {
    rout = r[[-1]]
    for (i in 1:nlayers(rout)) 
        rout[[i]] = rout[[i]] - r[[i]]
    return(rout)
}
