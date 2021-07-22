source("cfg.r")

extent = c(-180, 180, -35, 35)
years = 2003:2020
out_dir = "data/modInputs/"

#############
## gfas    ##
#############

file_emissions_gfas = "data/raw/cams_gfas_co2fire_2003-2020-timestep_annual-.nc"
dir_biomass_cci = "data/raw/biomass_cci/"

gfas = brick(file_emissions_gfas)
gfas = raster::crop(gfas, extent)
names(gfas) = years

fname = paste0(out_dir, 'obs_fireEmissions.nc')
gfas = writeRaster(gfas, fname, overwrite = TRUE)

###########
## GFED  ##
###########
gfed_dir = "data/raw/GFED4.1s/"

files = list.files(gfed_dir)

browser()

#############
## biomass ##
#############
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
jules_pattern = 'S3.Annual'
jules_vars = c("veg_c_fire_emission_gb", "burnt_carbon_dpm", "burnt_carbon_rpm", "cs_gb", "cv")

sec2year = 60*60*24*365
jules_vars_scale = c(sec2year, sec2year, sec2year, 1, 1)


forVarMod <- function(var, scale, mod) {
    
    dat = openMod(mod, jules_dir, var, years, scale, stream = jules_pattern)
    dat = raster::resample(dat, gfas[[1]])
    fname = paste0(out_dir, "mod_", var, ".nc")
    dat = writeRaster(dat, fname, overwrite = TRUE)
}


lapply(jules_run, function(i) mapply(forVarMod, jules_vars, jules_vars_scale, i))

