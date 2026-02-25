

setwd('Autotuning-NGD/')
files <- list.files('data')
library(ncdf4)
restom_file <- nc_open('data/lat_lon_10yr_180x360_ANN.nc')


e3sm_parameters <- ncvar_get(restom_file, 'lhs')
e3sm_parameters <- t(e3sm_parameters)
colnames(e3sm_parameters) <- ncvar_get(restom_file, 'x')


FLNT <- ncvar_get(restom_file, 'FLNT')
FSNT <- ncvar_get(restom_file, 'FSNT')
area <- ncvar_get(restom_file, 'area')
area_norm <-  area/sum(area) * 250
E3SM_RESTOM <- apply(area_norm * (FSNT - FLNT),3, sum)

PRECC <- ncvar_get(restom_file, 'PRECC')
PRECL <- ncvar_get(restom_file, 'PRECL')

PRECC <-  apply(area_norm *PRECC,3, sum)
PRECL <-  apply(area_norm *PRECL,3, sum)

save(E3SM_RESTOM, PRECC, PRECL, e3sm_parameters, file = 'data/E3SM_RESTOM.RData')



