datConvert <- function(dat, scale, annualAver = FALSE, trendC = FALSE, trendT = FALSE,
                       Seasonal = FALSE, modality = FALSE, cover = FALSE, 
                       modEG = NULL, ...) {
    RS <- function(dat) {
        if (!is.null(modEG)) dat = raster::resample(dat, modEG)
        dat[is.na(modEG)] = NaN
        return(dat)    }
    if (annualAver) {
        dat = mean(dat) * scale
        dat = RS(dat)
    } else if (trendC) {        
        dat = RS(dat)
        dat = makeTrendCoe(dat * scale, ...)
    } else if (trendT) {     
        dat = RS(dat)
        dat = makeTrendCoe(dat * scale, TRUE, ...)
    } else if (Seasonal) {
        dat = PhaseConcMod(dat, modEG)
    } else if (cover) {
        dat = sum(dat)
        dat = RS(dat)
    }
    return(dat)
}
