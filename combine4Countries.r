source("cfg.r")
load('data/CountryInfo.Rd')
load('outputs/minP.Rd')

yearsLinePlot = 2010:2020

TCthreshold = 0.5

FC = brick('outputs/ForestCarbon-withHumans.nc')
EM = brick('outputs/totalFireEmissions-withHumans.nc')

mask = brick("data/modInputs/obs_ABG_carbon.nc")[[1]] < minP
TC = raster('../jules_benchmarking/data/TreeCover.nc')/0.8
TC = convert_pacific_centric_2_regular(TC)
TC = raster::resample(TC, FC[[1]])
mask =  TC < TCthreshold
if (F) {
colsFC = c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476',
           '#41ab5d','#238b45','#006d2c','#00441b')

limitsFC = c(0, 10, 100, 200, 300, 400, 500) * 10

colsEM = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
           '#fc4e2a','#e31a1c','#bd0026','#800026')

limitsEM = c(0, 0.1, 0.5, 1, 5, 10, 50, 100)

FCscale = 10
EMscale = 1000

FC4plot = FC[[(nlayers(FC)-4):nlayers(FC)]]
EM4plot = EM[[rev(seq(nlayers(EM), 1, by = -12)[1:5])]]

png("figs/carbonYrMaps.png", res = 300, units = 'in', height = 5, width = 7.2)
    layout(cbind(1:6, 7:12), heights = c(rep(1, 5), 0.3))
    par(mar = rep(0, 4), oma = c(0, 0, 2, 0))

    plotMaps <- function(r, cols, limits, scale = 1, years =  2016:2020) {
        mask = raster::crop(mask, r)
        r[mask] = NaN
        r = r * scale
        mapply(plotStandardMap, layers2list(r), years,
               MoreArgs = list(limits = limits, cols = cols))
        StandardLegend(cols, limits, r[[1]], extend_max = TRUE, units = 'gC ~m-2~')
    }

    plotMaps(FC4plot, colsFC, limitsFC, FCscale)
    plotMaps(EM4plot, colsEM, limitsEM, EMscale, years = '')
    mtext('Biomass carbon', outer = TRUE, adj = 0.33)
    mtext('Fire emission' , outer = TRUE, adj = 0.67)
    
dev.off()
}
forCountry <- function(id, name) {
    #if (id == 5) browser()
    mask = (!mask) & ctry == id

    totCarb <- function(i, mask, marea) {
        print("yay")
        sum(i[mask] * marea, na.rm = TRUE)
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

    tfile = paste0("temp/countryEmission-", TCthreshold, '-', id, ".Rd")
    if (file.exists(tfile)) load(tfile)
    else {
        fc = totsCarbs(FC, FCscale)
        em = totsCarbs(EM, EMscale)
        save(fc, em, file = tfile)
    }
    plot(range(fc[,1]), c(0, 1.1*max(fc[,2])), type = 'n', xlab = '', ylab = '', yaxs = 'i')

    addBar <- function(xy, dx = 0.5) 
        polygon(xy[1] + dx * c(-1, 1, 1, -1, -1), xy[2] * c(0, 0, 1, 1, 0), col = "#99DD77")
    
    apply(fc, 1, addBar)
    dfc = cbind(fc[-1,1] - diff(fc[,1])/2, -diff(fc[,2]))

    rEM = range(c(dfc[,2], em[,2]))
    #rEM = range(em[,2])
    scaleFun <- function(y)  par("usr")[4] * (y - rEM[1])/(1.1*(rEM[2] - rEM[1]))
    
    addLine <- function(xy, col) {
        x = xy[,1]
        y = xy[,2]
        y = scaleFun(y)
        lines(x, y, lwd = 2, col = col)
    }
    addLine(em, 'red')
    addLine(dfc, 'blue')

    for (i in 1:5) {
        labels = round(seq(rEM[1], rEM[2], length.out = 6), i)
        labels = unique(labels)
        if (length(labels) >= 5) break
    }
    at = scaleFun(labels)
    axis(side = 4, labels = labels, at = at)

    mtext(name, side = 3, line = -1.3, font = 2)
}

png("figs/CountryTS.png", height = 7.2, width = 7.2, res = 300,  units = 'in')
    par(mfrow = c(3, 2), mar = rep(2, 4), oma = c(0, 2,0, 2))
    mapply(forCountry, 1:length(cntryNames), names(cntryNames))
    mtext(outer = TRUE, side = 2, 'Forest carbon (Pg)', line = 0.3)
    mtext(outer = TRUE, side = 4, 'Fire emissions/carbon loss (Pg/yr)', line = 0.3)

    plot.new()
    legend('left', c('Forest carbon', 'Forest carbon\nloss', 'fire emissions'), 
           pch = c(15, NaN, NaN), col = c('#99DD77', 'blue', 'red'), 
           lty = c(0, 1, 1), lwd = c(0, 2, 2), pt.cex = c(3, 0, 0))
dev.off()
