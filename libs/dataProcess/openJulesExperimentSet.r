openJulesExperimentSet <- function(name, JULES_control, JULES_experiments, 
                                   ..., grab_cache = TRUE) { 
    temp_file = paste0('temp/plot_jules_diff_from_control_new', name, '.Rd')
    if (file.exists(temp_file) & grab_cache) {
	load(temp_file)
    } else {        
	Jules_TC_control = openJulesTree(JULES_control, JULES_dir, ...)
	Jules_PFTs       = openJulesTree(JULES_control, JULES_dir, splitPFTs = TRUE, ...) 

    Jules_BA_control = openJulesTree(JULES_control, JULES_dir, 1,
                                 "burnt_area_gb", yrs = yrs) *   60 * 60 * 24 * 365 
	Jules_BA_exp     = layer.apply(JULES_experiments, openJulesTree, JULES_dir, 1,
                                 "burnt_area_gb", yrs = yrs) *   60 * 60 * 24 * 365 
                         
	Jules_BA_control = raster::resample(Jules_BA_control, Jules_TC_control)
	Jules_BA_exp     = raster::resample(Jules_BA_exp    , Jules_TC_control)
	Jules_TC_exp     = layer.apply(JULES_experiments, openJulesTree, JULES_dir, ...)
		
	Jules_dout = (Jules_TC_control - Jules_TC_exp)/Jules_TC_control
	Jules_dout =  squeeze(Jules_dout, 500)
	
	save(Jules_TC_control, Jules_BA_control, Jules_BA_exp, Jules_TC_exp, Jules_dout, Jules_PFTs,
             Jules_PFTs, file = temp_file)
    }
    return(list(Jules_TC_control, Jules_TC_exp, Jules_BA_control, Jules_BA_exp, Jules_dout, Jules_PFTs))
}
