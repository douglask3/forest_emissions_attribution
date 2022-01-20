source("cfg.r")

extent = c(-180, 180, -35, 35)
years = 2003:2020
mod_years = 2000:2019
out_dir = "data/modInputs/"

TC = raster('../jules_benchmarking/data/TreeCover.nc')
TC = convert_pacific_centric_2_regular(TC)

grab_cache = TRUE

gfas_fname = paste0(c(out_dir, 'obs_fireEmissions-', range(years), '.nc'), collapse = '')
file_emissions_gfas = 'data/raw/cams_gfas_co2fire_2003-2020-timestep_monthly-.nc'
#############
## gfas    ##
#############
if (file.exists(gfas_fname) && grab_cache) {
    gfas = brick(gfas_fname) 
} else {
    print("processing GFAS")
    gfas = brick(file_emissions_gfas)

    memSafeFile.initialise("temp/")
    gfas = gfas0 = layer.apply(gfas, memSafeFile.crop, extent, overwrite = TRUE)

    nms =  as.Date(paste(rep(years, each = 12), 1:12, 15, sep = '-'))
    gfas = gfas[[1:length(nms)]]
    gfas = setZ(gfas, nms, 'Date')

    gfas = writeRaster.gitInfo(gfas, gfas_fname, overwrite = TRUE)
}

###########
## GFED  ##
###########
if (T) {
gfed_dir = "data/gfed4.1s_nc/"
files = list.files(gfed_dir)

openResample <- function(file) {
    fname = paste0(out_dir, "obs_", file)
    
    if (file.exists(fname) && grab_cache) return(brick(fname))
    print(file)
    
    dat = brick(paste0(gfed_dir, '/', file))
    dat = raster::crop(dat, extent)
    dat = raster::resample(dat, gfas[[1]])
    dat[dat<0] = 0
    names(dat) = years[1:nlayers(dat)]    
    writeRaster(dat, fname, overwrite = TRUE)
    return(dat)
}

dats = sapply(files[1:2], openResample)
}

#############
## biomass ##
#############
dir = "data/raw/biomass_cci/"
files = list.files(dir, full.names = TRUE)
dat = stack(files)
dat = raster::resample(dat, gfas[[1]])
yr = sapply(files, function(i) strsplit(strsplit(i, '100m-')[[1]][2], '-')[[1]][1])

nms = paste(yr, 6, 15, sep  = '-')
dat[dat<0] = 0

setZ(dat, nms, 'Date')

fname = paste0(out_dir, 'obs_ABG_carbon.nc')
dat = writeRaster(dat, fname, overwrite = TRUE)

###########
## jules ##
###########
jules_dir = "/hpc/data/d05/cburton/jules_output/"
jules_run = c(withHuman = "u-cb972", noHuman = "u-cb972_NatIgn")
modFacVar = 'frac'
sec2year = 60*60*24*365

forVarMod <- function(var, levels, vname, scale, correct = NULL, mod, name) {   
    forLevel <- function(level, v) {
        print(level)
        openMod(mod, jules_dir, v, mod_years, scale, levels = level, stream = jules_pattern)
    }
    if (length(levels) > 1) {    
        
        dat  = lapply(levels, forLevel, var)
        frac = lapply(levels, forLevel, modFacVar)
        dat = mapply('*', dat, frac)
        
        dat = sumList.raster(dat)
        if (!is.null(correct)) {
            fracS = sumList.raster(frac)
            if (nlayers(correct) == 1) fracS = mean(fracS)
            dat = dat * correct/fracS
        }
    } else {
        dat = forLevel(levels, var)
    }

    dat = raster::resample(dat, gfas[[1]])
    if (nlayers(dat) == length(mod_years)) {
        names(dat) = mod_years
        setZ(dat, as.Date(paste(mod_years, 6, 15, sep = '-')), 'Date')   
    } else {
        mnths = c(paste0('0', 1:9), 10:12)    
        names(dat) =  paste0(mnths, '-', rep(mod_years, each = 12))
    }
    fname = paste0(levels, collapse = '-')
    fname = paste0(out_dir, "mod_", name, '_', vname, '_', var, '_', fname, ".nc")
    print(fname)
    if (is.na(max(dat[[1]][], na.rm = TRUE))) browser()
    dat = writeRaster.gitInfo(dat, fname, overwrite = TRUE)
}

RunJules <- function() 
    mapply(function(i, j) mapply(forVarMod, jules_vars, levels, names, jules_vars_scale, 
                                 corrects, i, j), 
           jules_run, names(jules_run))



jules_pattern = 'S2.Annual'

jules_vars = c("c_veg", "c_veg", "cs_gb", "cs")
names = paste0(c("ForestCarbon", "NoneForestCarbon", "deadCarbon", "DPM"), "-noLU")
jules_vars_scale = c(1, 1, 1, 1)
levels = list(1:5, 7:13, 1, 1)
corrects = list(NULL, NULL, NULL, NULL)
#RunJules()

jules_pattern = 'S3.Annual'
names = paste0(c("ForestCarbon", "NoneForestCarbon", "deadCarbon", "DPM"), "-LU")
#RunJules()

jules_pattern = 'S2.Annual'

jules_vars = c("c_veg", "c_veg", "cs_gb", "cs")
names = c("ForestCarbon", "NoneForestCarbon", "deadCarbon", "DPM")
jules_vars_scale = c(1, 1, 1, 1)
levels = list(1:5, 7:13, 1, 1)
corrects = list(TC, 1 - TC, NULL, NULL)
#RunJules()



jules_vars = c("cs", "cs_gb", "cv")
names = c("DPM", "deadCarbon", "vegCarbon")
jules_vars_scale = c(1, 1, 1)
levels = c(1, 1, 1)
corrects = list(NULL, NULL, NULL)
#RunJules()

jules_pattern = 'S3.ilamb'
jules_vars = c("burnt_area_gb")
jules_vars_scale = c(sec2year)
levels = c(1)
names = jules_vars
RunJules()

jules_pattern = 'S3.ilamb'
jules_vars = c("veg_c_fire_emission_gb", "burnt_carbon_dpm", "burnt_carbon_rpm")
jules_vars_scale = c(sec2year, sec2year, sec2year)
levels = c(1, 1, 1)
names = jules_vars
#RunJules()

jules_pattern = 'S3.ilamb'
jules_vars = c("t1p5m_gb", "smc_tot", "runoff")
jules_vars_scale = c(1, 1, sec2year)
levels = c(1, 1, 1)
names = jules_vars
#RunJules()

jules_pattern = 'S3.Monthly'
jules_vars = c("fch4_wetl")
jules_vars_scale = c(1)
levels = c(1)
names = jules_vars
#RunJules()


jules_pattern = 'S3.ilamb'
jules_vars = c("frac")
jules_vars_scale = c(1)
levels = list(1:5, 7:13)
names = c("TreeCover", "nonTreeCover")
#RunJules()

jules_pattern = 'S2.ilamb'
jules_vars = c("frac")
jules_vars_scale = c(1)
levels = list(1:5, 7:13)
names = c("TreeCover_noLU", "nonTreeCover_noLU")
#RunJules()



jules_pattern = 'S3.ilamb'
jules_vars = c("npp_gb", "resp_s_to_atmos_gb")
jules_vars_scale = c(sec2year, sec2year)
levels = c(1, 1)
names = c("npp_gb-LU", "resp_s_to_atmos_gb-LU")
#RunJules()

jules_pattern = 'S2.ilamb'
names = c("npp_gb-noLU", "resp_s_to_atmos_gb-noLU")
#RunJules()


