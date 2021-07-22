library(ncdf4)
process.jules.file <- function(file, level, varName) {
    print(file)
    

    nc = try(nc_open(file))
    if (class(nc) == "try-error") return(NULL)
    vars = names(nc$var)
	
    if (all(vars != varName)) {
	noFileWarning(c(), varName)
	return(NULL)
    }
	
    getVar <- function(var) {
	var = nc$var[[which(vars == var)]]
	dat = ncvar_get( nc, var)
	return(dat)
    }
	
    dat = getVar(varName)
    lat = getVar("latitude")
    lon = getVar("longitude")
    tim = getVar("time_bounds")
	
    l = length(lat)
	
    multiLayer <- function(mn, leveli = level) {
	mdat = dat[, leveli, mn]
	if (!is.null(dim(mdat)))
	    mdat = apply(mdat,1 , sum)
	return(mdat)
    }
	
    singleLayer <- function(mn) dat[, mn]
	
    monthizeData <- function(mn, FUN, ...) {
	mdat = FUN(mn, ...)
	r = rasterFromXYZ(cbind(lon, lat, mdat))
	return(r)
    }
    
    if (length(dim(dat)) == 2 && ncol(dat) == 12)
        r = layer.apply(1:12, monthizeData, singleLayer)
    else if (length(dim(dat)) == 3) {
	if (varName != modFacVar) {
	    openWeightLayer <- function(fracLevel) {
	        frac = process.jules.file(file, fracLevel, modFacVar)
	        r = layer.apply(1:12, monthizeData, multiLayer, fracLevel)
	        frac * r 
	    }
	    ri = mapply(openWeightLayer, 1:13)            
	    r = ri[[1]]
            for (i in ri[-1]) r = r + i		
            
	} else r = layer.apply(1:12, monthizeData, multiLayer)
    } else r = sum(rasterFromXYZ(cbind(lon, lat, dat))[[level]])
        	
    return(r)
}
