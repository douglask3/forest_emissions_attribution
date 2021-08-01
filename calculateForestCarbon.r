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

Mobs = mean(obs[[-1]])
Mobs = Mobs[!is.na(Mobs)]
Mobs = Mobs[Mobs>10 & Mobs < 400]

png("figs/obsCarbon.png", height = 5.5, width = 5, res = 300, units = 'in')
layout(matrix(1:5, nrow = 5), heights = c(2.5, 1, 1, 1, 0.3))
par(mar = c(5.5, 0, 0, 0), oma = rep(0.5, 4))

mc = max(hist(Mobs, 100, yaxt = 'n', xlab = '', ylab = '', main = '')$counts)
Cobs = knots(ecdf(Mobs))
y = seq(0, mc, length.out = length(Cobs))

lines(Cobs, y, lwd = 2)
minP = Cobs[which.max(diff(Cobs[0:300000]))]
lines(c(minP, minP), c(0, 9E9), lty = 2)
mtext.units(side = 1, line = 3, 'Above ground carbon (kg/~m2~)')

cols = c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45',
         '#006d2c','#00441b')

limits =  c(0, 10, 20, 50, 100, 50, 200, 250, 300)
plotFUN <- function(r, title) {
    plotStandardMap(r, cols = cols, limits = limits)
    contour(r, levels = minP, drawlabels=FALSE, add = TRUE, col = '#c51b7d')
    mtext(side = 3, title)
}
par(mar = rep(0, 4))
mapply(plotFUN, layers2list(obs), obsYrs)
StandardLegend(cols, limits, obs[[1]], extend_max = TRUE, add = FALSE)
dev.off()

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

#corrected = mods[[1]][[1]]
#corrected[!is.na(corrected)] = 0
corrected = log(mods[[1]][[yrC[1]]]) - log(obs[[1]])
corrected[corrected < -20] = -20
     

rasterizeCorrection <- function(i) {
    corrected[mask[index]] = as.numeric(corrV[i,])
    corrected[corrected>200] = 200 
    corrected
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

#modss = lapply(modFile_correct, lapply, brickDIR)

forMods <- function(mods, nm) {
    outs = lapply(1:3,
                 function(i) layer.apply(1:nlayers(corrected), corrTotCarbon, i, mods, nm))
    out3 = layer.apply(1:nlayers(corrected), corrTotCarbon, 1:3, mods, nm)
    writeRaster.gitInfo(out3, paste0("outputs/ForestCarbon-", nm, '.nc', overwrite = TRUE))
    return(out3)
}
outs = mapply(forMods, modss, c("withHumans", "noHumans"))

