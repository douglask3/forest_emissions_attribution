source("cfg.r") 

library(rgdal)

cntryNames = c(Bolivia = 'Bolivia', Brazil = 'Brazil',
               "DR Congo" = 'Democratic Republic of the Congo',
               Indonesia = 'Indonesia', Paraguay = 'Paraguay')


FC_file = 'outputs/ForestCarbon-withHumans.nc'
ctry_file = 'data/WB_countries_Admin0_10m/'

FC = brick(FC_file)

ctryS = readOGR(ctry_file)
ctry0 = ctry = rasterize(ctryS, FC[[1]])

id = sapply(cntryNames, function(i) which(grepl(i, ctryS[[5]])))

ctry[!is.na(ctry)] = 0
for (i in 1:length(id)) ctry[ctry0 == id[i]] = i

writeRaster.gitInfo(ctry0, "outputs/AllCountries.nc", overwrite = TRUE)
writeRaster.gitInfo(ctry , "outputs/HighlightCountries.nc", overwrite = TRUE)

save(ctry, ctry0, ctryS, id, cntryNames, file = 'data/CountryInfo.Rd')
