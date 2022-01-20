source("cfg.r")
load('data/CountryInfo.Rd')
load('outputs/minP.Rd')
load("outputs/forestDef.Rd")

yearsLinePlot = 2010:2020

TCthreshold = 0.05#'Olson'

addSoilC = TRUE

brickQ <- function(file) {
    print(file)
    dat = brick(file)
    nms = names(dat)
    #browser()
    if (max(dat[[1]][], na.rm = TRUE)<1000) return(dat)
    qX = quantile(dat[[nlayers(dat)-1]], 0.95)
    
    dat = layer.apply(dat, function(i) {i[i>qX] = qX; i})
    names(dat) = nms
    dat
}
if (F){
FC = brickQ('outputs/ForestCarbon-withHumans.nc')
EM = brickQ('outputs/totalFireEmissions-withHumans.nc')
TC = brickQ('outputs/cci_corrected-mod_withHuman_ForestCarbon_c_veg_1-2-3-4-5.nc')
NC = brickQ("outputs/cci_corrected-mod_withHuman_NoneForestCarbon_c_veg_7-8-9-10-11-12-13.nc")
DC = brickQ("outputs/cci_corrected-mod_withHuman_deadCarbon_cs_gb_1.nc")
RP = brickQ("outputs/cci_corrected-mod_noHuman_DPM_cs_1.nc")
}
if (is.numeric(TCthreshold)) { 
    TC = raster('../jules_benchmarking/data/TreeCover.nc')/0.8
    TC = convert_pacific_centric_2_regular(TC)
    TC = raster::resample(TC, FC[[1]])
    mask =  TC < TCthreshold
} else if (TCthreshold == 'minP') {
    mask = brick("data/modInputs/obs_ABG_carbon.nc")[[1]] < minP
} else if (TCthreshold == 'Olson') {
    mask = !forestBiome
} else {
    browser()
}
ctry = raster::resample(ctry, mask)

colsFC = c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476',
           '#41ab5d','#238b45','#006d2c','#00441b')
dcolsFC = c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5',
            '#c7eae5','#80cdc1','#35978f','#01665e','#003c30')
limitsFC = c(0, 10, 100, 200, 300, 400, 500) * 10
limitsFC = c(0, 50, 100, 150, 200, 250, 300, 350) * 100/2
dlimitsFC = c(-100, -50, -20, -10, -5, -2, -1, 1, 2, 5, 10, 20, 50, 100) * 10

colsEM = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
           '#fc4e2a','#e31a1c','#bd0026','#800026')

limitsEM = c(0, 0.1, 0.5, 1, 5, 10, 50, 100)
limitsEM = c(0, 5, 10, 20, 50, 100, 200, 500, 1000)

FCscale = 100/2
EMscale = 1000

carbCols = c("#BBFF99", "#99DD77", "#77AA44")
carbCols = c('#d9f0a3','#addd8e','#78c679','#41ab5d')

if (F) {

plot3ColMaps <- function(fname, rs, rscales, limitss, colss,   titles) {
    png(paste0("figs/", fname, ".png"), res = 300, units = 'in', height = 5*6.3/5.3, width = 7.2*3/2)
    layout(cbind(1:7, 8:14, 15:21), heights = c(rep(1, 6), 0.3))
    par(mar = rep(0, 4), oma = c(0, 0, 2, 0))

    plotMaps <- function(r, cols, limits, scale = 1, title, years =  2015:2019) {
        mask = raster::crop(mask, r)
        r[mask] = NaN
        r = r * scale
        plotStandardMap(mean(r), 'Mean', limits = limits, cols = cols)
        mtext(title, side = 3)
        mapply(plotStandardMap, layers2list(r), years,
               MoreArgs = list(limits = limits, cols = cols))
        StandardLegend(cols, limits, r[[1]], extend_max = TRUE, units = 'gC ~m-2~',
                      extend_min = c(FALSE, TRUE)[1 + (limits[1] < 0)])
    }
    #browser()
    mapply(plotMaps, rs, colss, limitss, rscales, titles, list(2015:2019, '', ''))
    #plotMaps(FC4plot[[-1]], colsFC, limitsFC, FCscale)
    #plotMaps(dFC4plot, dcolsFC, dlimitsFC, FCscale, years = '')
    #plotMaps(EM4plot, colsEM, limitsEM, EMscale, years = '')
    
    #mtext(, outer = TRUE, adj = 0.25)
    #mtext(, outer = TRUE, adj = 0.5)
    #mtext( , outer = TRUE, adj = 0.75)
    
    dev.off()
}
FC4plot = TC[[(nlayers(TC)-5):nlayers(TC)]]
dFC4plot = diff.raster(FC4plot)
EM4plot = EM[[rev(seq(nlayers(EM), 1, by = -12)[1:5])]]
plot3ColMaps('totFosC_and_EM', list(FC4plot[[-1]], dFC4plot, EM4plot), 
                               list(FCscale, FCscale, EMscale),
                               list(limitsFC, dlimitsFC, limitsEM),
                               list(colsFC, dcolsFC, colsEM),
                               c('Forest carbon', 'Change in forest carbon', 'Fire emission'))

#FCs = 
}

forCountry <- function(id, name) {
    #if (id == 5) browser()
    mask = (!mask) & abs(ctry - id)<0.5
    maskS = mask
    maskS[is.na(maskS)] = 0
    maskS = which(maskS[]==1)
    totCarb <- function(i, mask, marea) {
        print("yay")
        sum(i[maskS] * marea, na.rm = TRUE)
    }
    totsCarbs <-function(r, scale) {
        
        yrs = substr(names(r), 2, 5)
        tst = which(sapply(yrs, function(i) any(i == yearsLinePlot)))
               
        r = r[[tst]]
        
        date = as.numeric(substr(names(r), 2, 5)) + as.numeric(substr(names(r), 7, 8))/12 + 
               as.numeric(substr(names(r), 10, 11))/365
        
        mask = raster::crop(mask, extent(r[[1]]))
        marea = raster::area(mask)[mask]
                
        out = scale * unlist(layer.apply(r, totCarb, mask, marea))/1000000000
        return(cbind(date, out))
    }

    tfile = paste0("temp/countryEmission-3", TCthreshold, '-', id, ".Rd")
    #cfile = paste0("outputs/countryEmission", TCthreshold, '-', name, ".csv")
    if (file.exists(tfile)) load(tfile)
    else {
        em = totsCarbs(EM, EMscale)
        #em = sapply(1:10, function(yr) sum(em[seq((yr-1)*12 +1, 12*yr),2])) 
        fc = totsCarbs(FC, FCscale)
        nc = totsCarbs(NC, FCscale)
        nc[2,2] = nc[1,2]
        dc = totsCarbs(DC, FCscale)
        rp = totsCarbs(RP, FCscale)
        
        tcc = fc[,2] - nc[,2] - dc[,2]
        #out = cbind(fc, tcc, nc[,2], rp[,2], dc[,2] - rp[,2], em)
        #`=colnames(out) = c('date', 'Total carbon', 'Tree live carbon', 'GrassShrub live carbon', 'DPM', 
        #                  'Dead/soil carbon', 'Fire emissions')
        #write.csv(out, file = cfile)
        save(fc, nc, dc, rp, em, file = tfile)
    }
    
    #dev.new()
    
    addBar <- function(xy, dx = 0.5, col = "#99DD77", ...) 
        polygon(xy[1] + dx * c(-1, 1, 1, -1, -1), xy[2] * c(0, 0, 1, 1, 0), col = col, ...)
    
    nc0 = nc
    tc = fc; tc[,2] = tc[,2] - nc[,2] - dc[,2]   

    nc[,2] = nc[,2] + dc[,2]
    sc = dc; sc[,2] = sc[,2] - rp[,2]

    if (!addSoilC) {
        fc[,2] = fc[,2] - sc[,2]
        nc[,2] = nc[,2] - sc[,2]
        dc[,2] = dc[,2] - sc[,2]
        sc[,2] = 0
    } 

    plot(range(fc[1:9,1]), c(0, 1.1*max(fc[1:9,2])), type = 'n', xlab = '', ylab = '', yaxs = 'i')

    fc1 = fc
    apply(fc, 1, addBar, 0.5, carbCols[1])
    fc1[-c(1, 8, 9),2] = 0
    
    apply(fc1, 1, addBar, 0.5, NULL, 10)
    apply(nc, 1, addBar, 0.5, carbCols[2])
    apply(dc, 1, addBar, 0.5, carbCols[3])
    apply(sc, 1, addBar, 0.5, carbCols[4])
    
    fc[,2] = fc[,2] + nc[,2] + dc[,2] + sc[,2]
    dfc = cbind(fc[-1,1] - diff(fc[,1])/2, -diff(fc[,2]))
    if (pcEmission) {
        em = em[sapply(fc[,1], function(i) which.min(abs(i - em[,1]))),]
        em[,2] = 100*em[,2]/fc[,2]
    }
    
    #ems = sapply(list(fc, tc, nc0, rp, sc) function(i) em[,2]/i[,2])
    rEM = range(c(dfc[,2], em[,2]))
    if (pcEmission) rEM = range(em[,2])
    scaleFun <- function(y)  par("usr")[4] * (y - rEM[1])/(1.1*(rEM[2] - rEM[1]))
    
    addLine <- function(xy, col) {
        x = xy[,1]
        y = xy[,2]
        y = scaleFun(y)
        lines(x, y, lwd = 2, col = col)
    }
    addLine(em, 'red')
    if (!pcEmission) addLine(dfc, 'blue')
     
    for (i in 1:5) {
        labels = round(seq(rEM[1], rEM[2], length.out = 6), i)
        labels = unique(labels)
        if (length(labels) >= 5) break
    }
    at = scaleFun(labels)
    axis(side = 4, labels = labels, at = at)

    mtext(name, side = 3, line = -1.3, font = 2)
}

FUN <- function(nme) {
png(nme, height = 7.2*3/3, width = 7.2*3/2, res = 300,  units = 'in')
    par(mfrow = c(3, 3), mar = rep(2, 4), oma = c(0, 2,0, 2))
    mapply(forCountry, 1:length(cntryNames), names(cntryNames))
    mtext(outer = TRUE, side = 2, 'Forest carbon (Pg)', line = 0.3)
    if(pcEmission) lab = 'Fire emissions (% forest carbon)'
    else lab = 'Fire emissions/carbon loss (Pg/yr)'
mtext(outer = TRUE, side = 4, lab, line = 0.3)

    plot.new()
    labs = c('Tree carbon', 'Grass/Shrub carbon', 'DPM')
    lty =  c(0, 0, 0, 0, 1, 1); lwd = c(0, 0, 0, 0, 2, 2); cex = c(3, 3, 3, 3, 0, 0)
    pch =  c(15,15, 15, 15, NaN, NaN)
    if (addSoilC) {
        labs = c(labs, 'Soil carbon') 
    } else {
        pch = pch[-1];  lty = lty[-1]; lwd = lwd[-1]; cex = cex[-1];
        carbCols = head(carbCols, -1)
    }

    legend('left', c(labs, 'Forest carbon loss', 'fire emissions'), 
           pch = pch, col = c(carbCols, 'blue', 'red'), lty = lty, lwd = lwd, pt.cex  = cex)
dev.off()
}
pcEmission = FALSE
#FUN("figs/CountryTS.png")
addSoilC = FALSE
FUN("figs/CountryTS-noSoilC.png")##
pcEmission = TRUE

addSoilC = TRUE
FUN("figs/CountryTS-pcEm.png")
addSoilC = FALSE
FUN("figs/CountryTS-noSoilC-pcEm.png")


