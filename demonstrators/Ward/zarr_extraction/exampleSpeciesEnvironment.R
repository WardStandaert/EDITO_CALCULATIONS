# load functions, but also cached stac catalog named stacCatalog
source("editoTools.R")
options("outputdebug"=c('L','M'))

#the cached stacCatalog is called 'EDITOSTAC'
load(file = "./data-raw/editostacv2.par")



# the file to process

acf <- S3FileSystem$create(
  anonymous = T,
  scheme = "https",
  endpoint_override = "s3.waw3-1.cloudferro.com"
)

eurobis=arrow::open_dataset(acf$path("emodnet/biology/eurobis_occurence_data/eurobisgeoparquet/eurobis_no_partition_sorted.parquet" ))
Abra_alba = eurobis |> filter(aphiaidaccepted==141433 ) |> collect() 
pts= Abra_alba %>% select(Latitude=latitude, Longitude=longitude)








#the requested timestep resolution of the dataset in milliseconds
#in this case we work with daily data (1day = 24*3600*1000 = 8640000 milliseconds)
#
timeSteps=c()

#lookups without a time value will fail for now.. needs to change because not all layers have a time component
#for those that do, the user should specify a time, or we need statistical or climatology data products, like bio-oracle
# if(!"Time" %in% colnames(pts))
#   pts$Time = as.POSIXct( "2020-01-01" ) 


#the requested parameters, names in the stac catalogue table field 'par' .. see unique(stacCatalogue$par) for a list 
parameters= c('thetao','elevation','Substrate')


#check if they all exist
for ( parameter in parameters) {
  if(! parameter %in% unique(EDITOSTAC$par) ) 
  { dbl("Unknown parameter ", parameter) 
    
  }
}


#add verbose= anything to get additional info on the positions ( par_x, par_y, par_z ) and time (par_t) found in the zarr files
enhanced_DF=enhanceDF(inputPoints = pts , requestedParameters = parameters, requestedTimeSteps = timeSteps, stacCatalogue = EDITOSTAC, verbose="on")







