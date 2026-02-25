source(file.path(Sys.getenv("MODULESHOME"), "init/R"))
module("load", "netcdf-c")
#install.packages('ncdf4')
library(ncdf4)

scratch_dir <- ''

extension <- 'CESM/d651076/ppe/cam_ppe/rerun_PPE_250/PD/PD_timeseries/'
extension2 <- 'CESM_FSNT/d651076/ppe/cam_ppe/rerun_PPE_250/PD/PD_timeseries/'

sub_folders <- list.files(paste0(scratch_dir, extension))

further_extension <- c('/atm/hist/')
i = 1

file <- list.files(paste0(scratch_dir, extension, sub_folders[i], further_extension))
FLNT <- nc_open(filename =paste0(scratch_dir, extension, 
                                 sub_folders[i], further_extension,
                                 file), 
                write = F)
lat <- FLNT$dim$lat$vals
long <- FLNT$dim$lon$vals
time <- FLNT$dim$time$vals
init_vals <- (sin(pi / 180 * lat) - sin(pi / 180 *dplyr::lag(lat))) 
init_vals[1] <- 0
init_vals[length(init_vals) + 1] <- 0

area_lat <- (init_vals + dplyr::lag(init_vals))/2
area_lat <- area_lat[-1]

area_lat <- area_lat/sum(area_lat)

area_weights <- matrix(nrow=  length(long), ncol = length(lat), 
                       area_lat, byrow = T)
area_weights <- area_weights/sum(area_weights)
area_weights_all <- array(area_weights, dim = c(dim(area_weights), length(time)))
RESTOM_vals <- list()

for (i in 1:length(sub_folders)) {
  file <- list.files(paste0(scratch_dir, extension, sub_folders[i], further_extension))
  FLNT <- nc_open(filename =paste0(scratch_dir, extension, 
                                   sub_folders[i], further_extension,
                                   file), 
                  write = F)
  # FLNT_vals[[i]] <- apply(ncvar_get(FLNT, 'FLNT'), c(1,2), mean)
  
  file_FSNT <- list.files(paste0(scratch_dir, extension2, sub_folders[i], further_extension))
  
  FSNT <- nc_open(filename =paste0(scratch_dir, extension2, 
                                   sub_folders[i], further_extension,
                                   file_FSNT), 
                  write = F)
  # FSNT_vals[[i]] <- apply(ncvar_get(FSNT, 'FSNT'), c(1,2), mean)
  RESTOM_vals[[i]] <- apply(area_weights_all*(ncvar_get(FSNT, 'FSNT') - 
                                            ncvar_get(FLNT, 'FLNT')), c(3), sum)
  nc_close(FSNT)
  nc_close(FLNT)
  print(i)
}


hist(sapply(RESTOM_vals, mean))
summary(sapply(RESTOM_vals, mean))
#W/m2
RESTOM_vector <- sapply(RESTOM_vals, mean)

save(RESTOM_vals,RESTOM_vector, sub_folders, # lat, long, time, area_weights_all,
     file = 'data/CESM_data.RData')

pars_folder <- 'CESM_pars/'

list.files('CESM_pars/')
par_file <- nc_open(paste0(pars_folder, 'parameter_262_w_control_fix.nc'), write = F)

variables <- list()
for (i in 1:length(par_file$var)) {
  variables[[i]] <- ncvar_get(par_file, names(par_file$var)[i])
}
names(variables) <-  names(par_file$var)
variables_df <- as.data.frame(variables)
nc_close(par_file)


save(variables_df, 
     file = 'data/CESM_vars.RData')




