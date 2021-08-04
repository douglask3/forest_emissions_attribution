source("cfg.r")

extent = c(-180, 180, -35, 35)
years = 2003:2019
mod_years = 2000:2019
out_dir = "data/modInputs/"

#############
## gfas    ##
#############
if (F) {
file_emissions_gfas = 'data/raw/cams_gfas_co2fire_2003-2020-timestep_monthly-.nc'

gfas = brick(file_emissions_gfas)
gfas = raster::crop(gfas, extent)

mnths = c(paste0('0', 1:9), 10:12)    
nms =  paste0(mnths, '-', rep(years, each = 12))
gfas = gfas[[1:length(nms)]]
names(gfas) = nms

fname = paste0(out_dir, 'obs_fireEmissions.nc')
gfas = writeRaster(gfas, fname, overwrite = TRUE)
}
###########
## GFED  ##
###########
if (F) {
gfed_dir = "data/gfed4.1s_nc/"
files = list.files(gfed_dir)

openResample <- function(file) {
    print(file)
    browser()
    dat = brick(paste0(gfed_dir, '/', file))
    dat = raster::crop(dat, extent)
    dat = raster::resample(dat, gfas[[1]])
    dat[dat<0] = 0
    names(dat) = years[1:nlayers(dat)]

    fname = paste0(out_dir, "obs_", file)
    writeRaster(dat, fname, overwrite = TRUE)
    return(dat)
}

dats = sapply(files, openResample)
}
#############
## biomass ##
#############
dir = "data/raw/biomass_cci/"
files = list.files(dir, full.names = TRUE)
dat = stack(files)
dat = raster::resample(dat, gfas[[1]])
yr = sapply(files, function(i) strsplit(strsplit(i, '100m-')[[1]][2], '-')[[1]][1])
names(dat) = paste0('X', yr)
dat[dat<0] = 0

fname = paste0(out_dir, 'obs_ABG_carbon.nc')
dat = writeRaster(dat, fname, overwrite = TRUE)

###########
## jules ##
###########
jules_dir = "/hpc/data/d05/cburton/jules_output/"
jules_run = c(withHuman = "u-cb972", noHuman = "u-cb972_NatIgn")

sec2year = 60*60*24*365
jules_vars_scale = c(sec2year, sec2year, sec2year, 1, 1)

forVarMod <- function(var, vname, scale, mod, name) {    
    dat = openMod(mod, jules_dir, var, mod_years, scale, stream = jules_pattern)
    dat = raster::resample(dat, gfas[[1]])
    if (nlayers(dat) == length(mod_years)) {
        names(dat) = mod_years
    } else {
        mnths = c(paste0('0', 1:9), 10:12)    
        names(dat) =  paste0(mnths, '-', rep(mod_years, each = 12))
    }
    
    fname = paste0(out_dir, "mod_", name, '_', vname, '_', var, ".nc")
    print(fname)
    dat = writeRaster(dat, fname, overwrite = TRUE)
}

RunJules <- function() 
    mapply(function(i, j) mapply(forVarMod, jules_vars, names, jules_vars_scale, i, j), 
          jules_run, names(jules_run))


jules_pattern = 'S2.Annual'
jules_vars = c("cs", "cs_gb", "cv")
names = c("DPM", "deadCarbon", "vegCarbon")
RunJules()

jules_pattern = 'S3.ilamb'
jules_vars = c("veg_c_fire_emission_gb", "burnt_carbon_dpm", "burnt_carbon_rpm")
names = jules_vars



