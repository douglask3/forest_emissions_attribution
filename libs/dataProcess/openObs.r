openObs <- function(obs, Layers, scale, modEG, ..., layer = NULL) { 
    if (!exists("extent") || is.null(extent)) extent = c(-180, 180, -90, 90)
    
    if (is.list(Layers)) tLayers = paste(sapply(Layers, paste0, collapse = '_'), collapse = '-')
        else tLayers = range(Layers)
    tempFile = paste(c('temp/', filename.noPath(obs, TRUE),
                     tLayers, ..., scale, extent, '.nc'),
                     collapse = '-')
    cat("\nopening:", obs)
    cat("\n\tinto:", tempFile)
    
    if (file.exists(tempFile)) obs = brick(tempFile)
    else {
        
        openFile <- function(file) {
            obs = brick(file)
        
            openLayer <- function(Layer) {
                obs = obs[[Layer]]
                obs = convert_pacific_centric_2_regular(obs)     
                     
                obs = datConvert(obs, scale, modEG = modEG[[1]], ...)
                #obs[is.na(modEG)] = NaN
            
                return(obs)
            }
            if (is.list(Layers)) obs = layer.apply(Layers, openLayer) else
                obs = openLayer(Layers)
            return(obs)
        }
        if (length(obs) == 1) obs = openFile(obs) else obs = layer.apply(obs, openFile)
        
        writeRaster(obs, file = tempFile, overwrite = TRUE)
    }
    if (!is.null(layer)) obs = obs[[layer]]
    if (nlayers(obs) == 1) obs = obs[[1]]
    return(obs)
}
