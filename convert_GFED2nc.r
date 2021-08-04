################################################################################
## cfg                                                                        ##
################################################################################
## Libraries etc
source('cfg.r')
library(rhdf5)
    
dir = "data/raw/GFED4.1s/"

files = list.files(dir, full.names = TRUE)
files = files[grepl('GFED4.1s_', files)]

extent = c(-180, 180, -90, 90)


openHDFandConvert2Nc <- function(month, variable, file) {
        print(month)
	if (month < 10) month = paste('0', month, sep = '')
	print(file)
        print(variable)
	layer = paste('emissions', month, 'partitioning', variable, sep = '/')
       
	dat = h5read(file, layer)
	dat = raster(t(dat))
         
	#dat = aggregate(dat)
	extent(dat) = extent(extent)
	H5close()
	
	#writeRaster(dat, file = memSafeFile(), overwrite = TRUE)
	
	fDate = file.info(file)$ctime
	fDate = as.character(fDate)
	names(fDate) = file	
	fileDate <<- c(fileDate, fDate)
	
	return(dat)
}

convertVar <- function(var) {
    fileDate <<- c()
    if (substr(var, 1, 1) == "C") files = files[!grepl('beta', files)]

    forFile <- function(file) {
        dat = layer.apply(1:12, openHDFandConvert2Nc, var, file)
        #if (abs(max.raster(dat) - 1) < 0.00001) dat = mean(dat) else dat = sum(dat)
        return(dat)
    }
    dat = layer.apply(files, forFile)	
    
    yrs =  sapply(files, function(file) substr(strsplit(file, '1s_')[[1]][2], 1, 4))
    mnths = c(paste0('0', 1:9), 10:12)
    names(dat) = paste0(mnths, '-', rep(yrs, each = 12))
   
    names(fileDate) = paste(names(fileDate), 'obtained on')
				 
    comment = list('Data from GFEDv4.1s' = 	'Raw data file list on data/gfed/file_list.txt',
		   'Data obtained on' = fileDate)
    #dat = dat[[7:(nlayers(dat) - 6)]]
    fname = paste0("data/gfed4.1s_nc/", var, ".nc")        
    
    writeRaster.gitInfo.time(dat, fname,
                        comment = comment, overwrite = TRUE)
    # browser()
}

varnames = c("C_AGRI", "C_BORF", "C_DEFO", "C_PEAT", "C_SAVA", "C_TEMF", "DM_AGRI", "DM_BORF",
             "DM_DEFO", "DM_PEAT", "DM_SAVA", "DM_TEMF")
sapply(varnames, convertVar)
					
#memSafeFile.remove()

