source("cfg.r")

obs = brick('data/modInputs/obs_ABG_carbon.nc')

obsYrs = c(2010, 2017, 2018)

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

revN = gitVersionNumber()
save(minP, revN, file = 'outputs/minP.Rd')



