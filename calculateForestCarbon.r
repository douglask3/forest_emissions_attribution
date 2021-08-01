source("cfg.r")

dir = 'data/modInputs/'

obsFile = 'obs_ABG_carbon.nc'

obsYrs = c(2010, 2017, 2018)
modYrs = 2000:2019

modFiles_opt = c('mod_withHuman_vegCarbon_cv.nc', 'mod_withHuman_DPM_cs.nc')
modFile_correct = list(c(modFiles_opt, 'mod_withHuman_deadCarbon_cs_gb.nc'),
                       c('mod_noHuman_vegCarbon_cv.nc', 'mod_noHuman_DPM_cs.nc',
                         'mod_noHuman_deadCarbon_cs_gb.nc'))
tbase = 500


yrC = sapply(obsYrs, function(i) which (i == modYrs))
yr = 1:length(modYrs)

brickDIR <- function(file) brick(paste0(dir, '/', file))

obs = brickDIR(obsFile)
mods = lapply(modFiles_opt, brickDIR)
mask = which(!is.na(sum(obs + mods[[1]][[1]]))[])

correct = correct0 =  matrix(0, ncol = tbase, nrow = nlayers(mods[[1]]))
CorrectForGridCell <- function(i, id, modVeg, modDPM, returnRes = TRUE) {
    idf = tbase*ceiling(id/tbase)
    idi = id - idf + tbase
    
    tfile = paste0("temp/correct_ABGcarbon_withDMP", idf, ".csv")
    
    if (file.exists(tfile)) {
        if (idf == id || id == length(mask)) {  
            print(tfile)          
            if (returnRes) return(read.csv(tfile)[,-1]) else return(invisible())
        } else {
            return(invisible())
        }
    }
    
    M = modVeg[i][1,] + modDPM[i][1,]; O = obs[i][1,]    
    O[O<0] = 0; M[M<0] = 0
    if (all(O==O[1]) && all(M[yrC] == O[1])) C = rep(0, length(yr))
    else { 
        M = log(M + 0.00000001)
        O = log(O + 0.00000001)
        if (all(O==O[1]) && all(M[yrC] == M[yrC[1]])) {
            C = rep(M[yrC[1]] - O[1], length(yr))
        } else {
            D = M[yrC] - O

            offset = mean(D)
            D = D - offset
            sdev = sd(D)
            D = D /sdev

            if (length(yrC) < 4) {
                yrC = c(yrC + 0.0001, yrC - 0.0001)
                D = c(D, D)
            }
            fit = smooth.spline(yrC, D)

            C = predict(fit, yr)[[2]]*sdev+offset
            if (any(is.na(C))) browser()
        }
        correct[,idi] = C
    }
    if (idf == id || id == length(mask)) {
        write.csv(correct, file = tfile)
        print(tfile)
        print(length(mask))
        if (returnRes) return(correct) else return(invisible())
    }
    correct <<- correct
    return(invisible())
}

index = 1:length(mask)#572000#
out = mapply(CorrectForGridCell, mask[index], index,
             MoreArgs = list(mods[[1]], mods[[2]],  returnRes = FALSE))
out = mapply(CorrectForGridCell, mask[index], index, MoreArgs = list(mods[[1]]))
out = out[!sapply(out, is.null)]
corrV= do.call(cbind, out)

corrected = mods[[1]][[1]]
corrected[!is.na(corrected)] = 0
#corrected = log(mods[[1]][[yrC[1]]]) - log(obs[[1]])
#corrected[corrected < -20] = -20

rasterizeCorrection <- function(i) {
    tfile = paste0("temp/calForestCorrecetedVar-", i, ".nc")
    if (file.exists(tfile)) return(raster(tfile))
    corrected[mask[index]] = as.numeric(corrV[i,])
    corrected[corrected>200] = 200 
    corrected = writeRaster(corrected, file = tfile)
    return(corrected)
}

corrected = layer.apply(1:nrow(corrV), rasterizeCorrection)

corrTotCarbon <- function(i, id = 1:length(mods), mods, nme) {
    tfile = paste0(c("temp/JULES_tot_carbon_corrected-", nme, '-', i, '-', id, ".nc"),
                   collapse = '-')
    if (file.exists(tfile)) return(raster(tfile))
    print(i)
    
    totCarbon = sum(layer.apply(mods[id], function(mod) mod[[i]]))
    
    totCarbon[totCarbon<0] = 0
    totCarbon = log(totCarbon + 0.00000001)
    
    totCarbon = totCarbon - corrected[[i]]
    totCarbon = exp(totCarbon)
    totCarbon = writeRaster(totCarbon, tfile, overwrite = TRUE)
    
    return(totCarbon)
}

forMods <- function(mods, nm) {
    outs = lapply(1:3,
                 function(i) layer.apply(1:nlayers(corrected), corrTotCarbon, i, mods, nm))
    out3 = layer.apply(1:nlayers(corrected), corrTotCarbon, 1:3, mods, nm)
    nms = paste(modYrs, 6, 15,  sep = '-')
    out3 = setZ(out3, as.Date(nms), 'Date')    
    writeRaster.gitInfo(out3, paste0("outputs/ForestCarbon-", nm, '.nc'), overwrite = TRUE)
    return(out3)
}

modss = lapply(modFile_correct, lapply, brickDIR)
outs = mapply(forMods, modss, c("withHumans", "noHumans"))

