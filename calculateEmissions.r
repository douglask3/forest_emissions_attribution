source("cfg.r")

dir = 'data/modInputs/'

obsFile = 'obs_fireEmissions-20032020.nc'

years = 2003:2019
nms =  paste(rep(years, each = 12), 1:12, c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
             sep = '-')

obs = brick(paste0(dir, obsFile))
#browser()
#names(obs) = nms


writeRaster.gitInfo.time <- function(r, file, ...) {
    nms = names(r)
    nms = substr(nms, 2, 99)
    nms = gsub('.', '-', nms, fixed = TRUE)
    r = setZ(r, as.Date(nms), name = "Date")
    
    writeRaster.gitInfo(r, file = file, overwrite = TRUE, ...)
    #nc = nc_open(file, write = TRUE)
   
    #var = ncvar_def( 'Time', '', list(nc$dim[[3]]), NaN)
    #nc = ncvar_add(nc, var)	
    #ncvar_put(nc, nc$var[[length(nc$var)]], names(r))
    #nc_close(nc)
}
runningAv.raster <- function(dat, N = 12, aFUN = function(i) sum(i), outFile = NULL)  {
    if (!is.null(outFile))  memSafeFile.initialise()
    FUN <- function(mn) {
        print(mn)
        out = aFUN(dat[[(mn-N+1):mn]])
        if (!is.null(outFile))  out = writeRaster(out, file = memSafeFile(), overwrite = TRUE)
        print(filename(out))
        return(out)
    }
    
    outs = layer.apply(N:nlayers(dat), FUN) # 
    names(outs) = names(dat)[(N-round(N/2)):(nlayers(dat)-round(N/2))] #
    if (!is.null(outFile)) {
        outs = writeRaster.gitInfo.time(outs, outFile)
    
        memSafeFile.remove()
    }
    return(outs)
}

obs12 = runningAv.raster(obs, outFile = "outputs/totalFireEmissions-withHumans.nc")


