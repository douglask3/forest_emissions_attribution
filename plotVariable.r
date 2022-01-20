source("cfg.r")
options(error=recover)
graphics.off()
titles = c("Burnt area", "Surface temp", "Runoff", "Soil moisture", "Agriculature", "Vegetation carbon", "Soil carbon")

filess = list(
burntArea = list("JULES-ES" = "data/modInputs/mod_withHuman_burnt_area_gb_burnt_area_gb_1.nc",
                 "ConFire" = "../savanna_fire_feedback_test/data/LimFIRE/outputs/ensemble_156/LimFIRE_fire-rw.nc",
                 "GFED4s.1 Obs." = "../savanna_fire_feedback_test/data/LimFIRE/outputs/fire2000-2014.nc"),
tas = list("JULES-ES" = "data/modInputs/mod_withHuman_t1p5m_gb_t1p5m_gb_1.nc"),
runoff = list("JULES-ES" = "data/modInputs/mod_withHuman_runoff_runoff_1.nc"),
smc = list("JULES-ES" = "data/modInputs/mod_withHuman_smc_tot_smc_tot_1.nc"),
Agri = list(Pasture = "../savanna_fire_feedback_test/data/LimFIRE/outputs/pasture2000-2014.nc",
          Cropland ='../savanna_fire_feedback_test/data/LimFIRE/outputs/cropland2000-2014.nc'),
vegCarbon  = list("CCI carbon" = 'data/modInputs/obs_ABG_carbon.nc',
                  "JULES-ES"   = "data/modInputs/mod_withHuman_vegCarbon_cv.nc"),
soilCarbon = list("JULES-ES"   = "data/modInputs/mod_withHuman_cs_gb.nc"))

extent = c(-125, 155, -35, 30)

colss = list(burntArea = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026'),
             tas = rev(c('#d73027','#f46d43','#fdae61','#fee090',
                         '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4')),
            runoff = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6',
                       '#4292c6','#2171b5','#08519c','#08306b'),
            smc = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                     '#3690c0','#0570b0','#045a8d','#023858'),
             Agri       = c('#f7fcfd','#e5f5f9','#ccece6','#99d8c9',
                            '#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
             vegCarbon  = c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb',
                           '#41b6c4','#1d91c0','#225ea8','#253494','#081d58'),
             soilCarbon = c('#ffffe5','#fff7bc','#fee391','#fec44f',
                            '#fe9929','#ec7014','#cc4c02','#993404','#662506'))

limitss = list(burntArea = c(0.1, 1, 2, 5, 10, 20, 50),
               tas = c(16, 18, 20, 22, 24, 26, 28, 30),
               runoff = c(0, 100, 200, 400, 600, 800, 1000, 1500), 
               smc = c(0, 100, 200, 400, 600, 800, 1000, 1200),
               Agri = c(0, 1, 2, 5, 10, 20, 50, 60, 80),
               vegCarbon  = c(0, 40, 80, 120, 160, 200, 240, 280, 320),
               soilCarbon = c(0, 40, 80, 120, 160, 200, 240, 280, 320))

unitss = c('%', 'C', 'kg/m2/yr', 'kg/m2', '%', 'MgC/Ha', 'MgC/Ha')

openPlot <- function(vname, files, cols, limits, title, units) {
    print(vname)
    rs = lapply(files, function(file) mean(brick(file)[[1:12]]))
    rs = lapply(rs, crop, extent)
    if (vname ==  "vegCarbon") rs[[2]] = rs[[2]]*12
    if (vname == "soilCarbon") rs[[1]] = rs[[1]]*12
    if (vname == "Agri" || vname == "burntArea") rs = lapply(rs, function(i) i*100)
    if (vname == "tas") rs = lapply(rs, function(i) i-273.15)
    if (vname == "runoff") rs = lapply(rs, function(i) i*60*60*24*360)
    if (vname == "burntArea") {
        rs[[1]] = rs[[1]]*60*60*24*30
        rs = lapply(rs, function(i) i*12)
    }
    
    plotR <- function(name, r) {
        print(name)
        fname = paste("figs/map", vname, name, ".png", sep = '-')
        png(fname, height = 2.5, width = 5, res = 300, units = 'in')
            layout(rbind(c(1, 1, 4), 
                         c(3, 2, 4),
                         c(3, 2, 2), 
                         c(3, 5, 5)), width = c(2.3,2.4, 1.8), height = c(1.6, 0.4, 1, 0.5))
            par(mar = rep(0.2,4))
            pFUN <- function(x_range = extent[1:2], y_range = extent[3:4], ...) {
                print("plotting")
                plotStandardMap(r, '', limits = limits, cols = cols, 
                                x_range = x_range, y_range = y_range, ...)
                polyx = c(x_range[c(1, 2, 2, 1, 1)], -180, -180, 180, 180, -180)
                polyy = c(y_range[c(1, 1, 2, 2, 1)], -90, 90, 90, -90, -90)
                polygon(polyx, polyy, col = "white", border = NA)
            }
            pFUN()
            pFUN(x_range = c(92, 150), y_range = c(-11, 10))
            pFUN(x_range = c(-95, -30), y_range = c(-30, 12.5))
            pFUN(x_range = c(31, 42.5), y_range = c(-5, 5), cntrThinkness = 2)
            
            
            StandardLegend(cols = cols, limits =  limits, dat = r, units = units)
            mtext(side = 3, paste(title, "-", name), outer = TRUE, line = -1.2, adj = 0.1)
        dev.off()
    }
    mapply(plotR, names(files), rs)
}

mapply(openPlot, names(filess), filess, colss, limitss, titles, unitss)
