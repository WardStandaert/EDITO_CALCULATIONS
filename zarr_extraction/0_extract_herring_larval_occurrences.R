# load functions, but also cached stac catalog named stacCatalog
source("editoTools.R")
options("outputdebug"=c('L','M'))

#the cached stacCatalog is called 'EDITOSTAC'
load(file = "./data-raw/editostacv2.par")

library(arrow)

# the file to process

acf <- S3FileSystem$create(
  anonymous = T,
  scheme = "https",
  endpoint_override = "s3.waw3-1.cloudferro.com"
)

eurobis = arrow::open_dataset(acf$path("emodnet/biology/eurobis_occurence_data/eurobisgeoparquet/eurobis_no_partition_sorted.parquet" ))
Clupea_harengus = eurobis |> 
  filter(aphiaidaccepted==126417, datasetid==4423,
         longitude > -12, longitude < 10,
         latitude > 48, latitude < 62,
         observationdate >= as.POSIXct("2000-01-01"),
         observationdate <= as.POSIXct("2020-12-31")) |> 
  collect() 

glimpse(Clupea_harengus)

pts = Clupea_harengus %>% select(Latitude=latitude, 
                                 Longitude=longitude, 
                                 Time=observationdate)

glimpse(pts)
