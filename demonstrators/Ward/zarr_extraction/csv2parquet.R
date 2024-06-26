###script used for conversion of large csv to parquet
# # /usr/bin/time --verbose R CMD BATCH --no-save "csv2parquet.R" "csv2parquet.R_$(date +"%F").Rout" &

library(arrow)
library(geoarrow)
# library(sfarrow)
library(dplyr)

library(RPostgreSQL)


directory="~/parquet"
if (!dir.exists(directory)) {dir.create(directory)}
setwd(directory)

csv_file="obis_mrgid_export.csv"



####obis_mrgid collection from aphia mongodb####
##mongodb collection (obis_mrgid) to csv (25.8GB) to parquet (8.5GB)
#1) export a mongodb collection to csv, using studio3T (codewise looping and appending csv creates errors atm)
#2) read csv and write parquet with this script

##set schema

#all columns as string
# a=arrow::open_dataset("obis_mrgid_export.csv", format = "csv")
# s <- schema(
#   purrr::map(names(a), ~Field$create(name = .x, type = string()))
# )
# a=arrow::open_dataset("obis_mrgid_export.csv", format = "csv", schema = s, skip=1)
# arrow::write_dataset(a, "obis_mrgid", format = "parquet")

#set fitting type
a=arrow::open_dataset(csv_file, format = "csv")
s=a$schema
s[["eventDate"]]=Field$create("eventDate", string()) #eventDate as string, because problem with timezones: Error: Invalid: In CSV column #3: Row #20447321: CSV conversion error to timestamp[s, tz=UTC]: expected a zone offset in '2009-05-03T01:27:14'. If these timestamps are in local time, parse them as timestamps without timezone, then call assume_timezone.
s[["aphiaID"]]=Field$create("aphiaID", string()) #not all aphiaID have integers
a=arrow::open_dataset("obis_mrgid_export.csv", format = "csv", schema = s, skip=1)
# h=dplyr::slice_head(a, n=100000)

arrow::write_dataset(a, "obis_mrgid_types", format = "parquet") #uses 10GB RAM, 3-4 threads, 5min 
arrow::write_dataset(a, "obis_mrgid_partitions", format = "parquet", partitioning = "mrgid_mr")

# a=arrow::read_parquet("obis_mrgid/part-0.parquet", as_data_frame = F) #fully read parquet into memory
a=arrow::open_dataset("obis_mrgid") #open with binding
ap=arrow::open_dataset("obis_mrgid_partitions")
 

# #tests
# system.time({dplyr::distinct(ap, mrgid_mr) %>% collect()})
# # user  system elapsed
# # 64.220  13.723  37.530
# 
# system.time({m<<-dplyr::distinct(a, mrgid_mr) %>% collect()})
# # user  system elapsed
# # 12.504   1.269   2.057
# 
# 
# system.time({
#   x <<- dplyr::filter(ap, mrgid_mr %in% c("26565", "25650", "25653", "25663", "26564", "25647", "25652",
#                                           "26560", "25549", "25411")) %>% collect()
# })
# # user  system elapsed
# # 73.544   8.767   6.937 -> partitioning helps here
# 
# system.time({
#   x <<- dplyr::filter(a, mrgid_mr %in% c("26565", "25650", "25653", "25663", "26564", "25647", "25652",
#                                          "26560", "25549", "25411")) %>% collect()
# })
# # user  system elapsed
# # 123.331  10.022  11.535
# 
# # tibble [6,815,942 × 16]
# 
# x=a %>%
#   slice_head(n = 100000) %>%
#   dplyr::count(mrgid_mr, aphiaID) %>%
#   collect()
# 
# x=a %>%
#   dplyr::count(mrgid_mr, mrgid_iho, mrgid_meow, mrgid_eez) %>%
#   collect()
# x




####eurobis from postgres####
csv_file="eurobis.csv"

#query get postgres eurobis
pgdrv <- dbDriver(drvName = "PostgreSQL")
con <-DBI::dbConnect(pgdrv,
                     dbname="dataportal",
                     host="pgetn.vliz.be", port=5432,
                     user = '',
                     password = '')

#5min
data.table::fwrite(DBI::dbGetQuery(con,
  "SELECT
  obs.occurrenceid_unique AS occurrenceid,
  obs.datasetid,
  obs.observationdate,
  obs.longitude,
  obs.latitude,
  obs.scientificname,
  obs.aphiaid,
  aph.scientificname_accepted,
  obs.aphiaidaccepted,
  obs.taxonrank
  FROM eurobis.observations obs
  LEFT JOIN eurobis.aphia_matrix aph ON obs.aphiaid = aph.aphiaid"
), file=csv_file, sep=",", na="", row.names=FALSE, col.names=TRUE )





##partition on species: aphiaidaccepted, faster to filter on species, although max amount of species not tested, bbox query is really slow here, so making multiple datasets with different partition 
output_parquet_file="eurobis_species.parquet"
unlink(output_parquet_file,recursive=T)

a=arrow::open_dataset(csv_file, format = "csv", read_options = list(use_threads=T, block_size = 2048576L))
a=dplyr::arrange(a,aphiaidaccepted)
# arrow::write_dataset(a, output_parquet_file, format = "parquet", partitioning = c("latitude","longitude"), max_partitions = 500000) #stops at 16393 files

#workaround for Error: Invalid: Fragment would be written into 1082 partitions. This exceeds the maximum of 1024 -> increasing max_partitions option stops as 16k files
a=dplyr::collect(a)
x=split(a,a$aphiaidaccepted)
for (n in 1:length(x)) {
  if(n %% 1000==0) print(n)
  df=as.data.frame(x[n])
  colnames(df)=gsub("X.+\\.","",colnames(df))
  arrow::write_dataset(df, output_parquet_file, format = "parquet", partitioning = "aphiaidaccepted")
}




##partition on date: year, month, day
output_parquet_file="eurobis_date.parquet"
unlink(output_parquet_file,recursive=T)

a=arrow::open_dataset(csv_file, format = "csv", read_options = list(use_threads=T, block_size = 2048576L))
a=dplyr::arrange(a,observationdate) |> 
  dplyr::mutate(
    year=format(observationdate, format="%Y"),
    month=format(observationdate, format="%m"),
    day=format(observationdate, format="%d")) |>
  dplyr::collect()
# arrow::write_dataset(a, output_parquet_file, format = "parquet", partitioning = c("year","month","day"), max_partitions = 500000) #stops at 52001 files

x=split(a,a$observationdate)
for (n in 1:length(x)) {
  if(n %% 1000==0) print(n)
  df=as.data.frame(x[n])
  colnames(df)=gsub("X.+\\.","",colnames(df))
  arrow::write_dataset(df, output_parquet_file, format = "parquet",  partitioning = c("year","month","day"))
}

# ep=arrow::open_dataset(output_parquet_file)
# e=arrow::open_dataset("eurobis_no_partition.parquet/")
# 
# x=dplyr::filter(ep,year==2015) |>
#   dplyr::collect()
# 
# x2=dplyr::filter(e,
#   observationdate >= as.POSIXct("2015-01-01",tz="UTC"),
#   observationdate < as.POSIXct("2016-01-01",tz="UTC")
#   ) |>
#   dplyr::collect()



##geospatial index: this is not available at the moment and is being developed, tried to create one myself with rtree, but my implementation speeds up only small bbox queries (transform rtree to 1 id column, partition on this and filter on max 60 ids, otherwise it's slower) and you get approximate data (too much), working this out takes too much time, so for now sorting data on geometry column with no partition seems the quickest 
# https://github.com/opengeospatial/geoparquet/issues/13
# https://cloudnativegeo.org/blog/2023/10/the-admin-partitioned-geoparquet-distribution/


#sort on geometry column -> speeds up bbox query
output_parquet_file="eurobis_spatial.parquet"
unlink(output_parquet_file,recursive=T)

a=arrow::open_dataset(csv_file, format = "csv", read_options = list(use_threads=T, block_size = 2048576L))
a_sf=dplyr::collect(a) |> 
  sf::st_as_sf(coords = c("longitude","latitude"), crs=4326, remove=FALSE) |>
  dplyr::arrange(geometry)
a=sf::st_drop_geometry(a_sf)
arrow::write_dataset(a, output_parquet_file, format = "parquet")

#geoparquet with geometry column, but slower for bbox query, and can't find way to filter by polygon?
output_parquet_file="eurobis_spatial_geoparquet.parquet"
unlink(output_parquet_file,recursive=T)

geoarrow::write_geoparquet(a_sf, output_parquet_file) 

# aspe=arrow::open_dataset("eurobis_species.parquet/")
# ad=arrow::open_dataset("eurobis_date.parquet/")
# aspa=arrow::open_dataset("eurobis_spatial.parquet/")
# aspag=arrow::open_dataset("eurobis_spatial_geoparquet.parquet")


