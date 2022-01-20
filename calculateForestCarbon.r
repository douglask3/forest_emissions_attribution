source("cfg.r")

dir = 'data/modInputs/'

obsFile = 'obs_ABG_carbon.nc'

obsYrs = c(2010, 2017, 2018)
modYrs = 2000:2019

modFiles_opt = c('mod_withHuman_vegCarbon_cv.nc', 'mod_withHuman_DPM_cs.nc')

files = paste0(c('ForestCarbon_c_veg_1-2-3-4-5', 'NoneForestCarbon_c_veg_7-8-9-10-11-12-13', 
                 'deadCarbon_cs_gb_1', 'DPM_cs_1'), '.nc')

modFile_correct = list(paste0('mod_withHuman_', files),
                       paste0('mod_noHuman_', files))

modFiles_opt = modFile_correct[[1]][[1]]
tbase = 500

yrC = sapply(obsYrs, function(i) which (i == modYrs))
yr = 1:length(modYrs)

brickDIR <- function(file) brick(paste0(dir, '/', file))

obs = brickDIR(obsFile)
mods = brickDIR(modFiles_opt)
mask = which(!is.na(sum(obs + mods[[1]]))[])

correct = correct0 =  matrix(NaN, ncol = tbase, nrow = nlayers(mods))
CorrectForGridCell <- function(i, id, modVeg, returnRes = TRUE) {
    idf = tbase*ceiling(id/tbase)
    idi = id - idf + tbase
    
    tfile = paste0("temp/correct_ABGcarbon_ForestOnly", idf, ".csv")
    
    if (file.exists(tfile)) {
        if (idf == id || id == length(mask)) {  
            print(tfile)          
            if (returnRes) return(read.csv(tfile)[,-1]) else return(invisible())
        } else {
            return(invisible())
        }
    }
    
    M = modVeg[i][1,]; O = obs[i][1,]    
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

index = 1:length(mask)#272000#
out = mapply(CorrectForGridCell, mask[index], index,
             MoreArgs = list(mods,  returnRes = FALSE))
out = mapply(CorrectForGridCell, mask[index], index, MoreArgs = list(mods))

out = out[!sapply(out, is.null)]
out[[length(out)]] =  out[[length(out)]][,!is.na(out[[length(out)]][1,])]
corrV= do.call(cbind, out)

corrected = mods[[1]]
corrected[!is.na(corrected)] = 0
#corrected = log(mods[[1]][[yrC[1]]]) - log(obs[[1]])
#corrected[corrected < -20] = -20

rasterizeCorrection <- function(i) {
    
    tfile = paste0("temp/calForestCorrecetedVar-liveVeg", i, ".nc")
    print(tfile)
    if (file.exists(tfile)) return(raster(tfile))
    
    corrected[mask[index]] = as.numeric(corrV[i,])
    
    corrected = writeRaster(corrected, file = tfile)
    return(corrected)
}

corrected = layer.apply(1:nrow(corrV), rasterizeCorrection)

corrTotCarbon <- function(i, id = 1:length(mods), mods, nme) {
    tfile = paste0(c("temp/JULES_tot_carbon_corrected-liveVeg", nme, '-', i, '-', id, ".nc"),
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

forMods <- function(mods, nm, fnames) {
    outs = lapply(1:3,
                 function(i) layer.apply(1:nlayers(corrected), corrTotCarbon, i, mods, nm))
    out3 = layer.apply(1:nlayers(corrected), corrTotCarbon, 1:3, mods, nm)
    nms = paste(modYrs, 6, 15,  sep = '-')
    
    
    writeRaster.gitInfo.Dates <- function(r, ...) {
       
        r = setZ(r, as.Date(nms), 'Date')  
        writeRaster.gitInfo(r, ..., overwrite = TRUE)
    }
    fnames = paste0('outputs/cci_corrected-', sapply(fnames, filename.noPath))
    mapply(writeRaster.gitInfo.Dates, outs, fnames)

    writeRaster.gitInfo.Dates(out3, paste0("outputs/ForestCarbon-", nm, '.nc'))
    return(out3)
}

modss = lapply(modFile_correct, lapply, brickDIR)
outs = mapply(forMods, modss, c("withHumans", "noHumans"), modFile_correct)

