library(rstac)
library(purrr)
source("demonstrators/Ward/zarr_extraction/editoTools.R")

#1. Query collections ----
# Define STAC endpoint URL
stac_endpoint_url <- 'https://catalog.dive.edito.eu/'

# Perform a request to get all collections in the catalog
collections <- stac(stac_endpoint_url) %>%
  collections() %>%
  get_request()

# Print collections
print(collections)

f <- function(col, name) {
  collection_list <- col$collections
  matches <- keep(collection_list, \(col) {grepl(name, col$id, ignore.case = TRUE) |  grepl(name, col$title, ignore.case = TRUE)})
  return(matches)
}


f_timechunked <- function(col, name) {
  collection_list <- col$collections
  matches <- keep(collection_list, \(col) {grepl(name, col$id, ignore.case = TRUE) |  grepl(name, col$title, ignore.case = TRUE)})
  return(matches)
}

# alternative
# eunis <- collects$collections %>%
#   keep(~ str_detect(tolower(.x$title), "eunis")) %>%
#   .[[1]]

# Filtered collections in id & title
t <- f(collections, "broad")[[1]]


col_seabed_habitat <- f(collections, "broad")[[1]]
col_seabed_habitat

#2. Create link and extract raster ----
link <- rstac::assets_url(col_seabed_habitat, "Zarr")
#look at what data sets are in there
names(gdalinfo(link)$arrays)

#still uses EditoTools
dsn <- toGDAL(link)
sdsn=sprintf('%s:/%s',dsn,"marine_strategy_framework_directive_benthic_broad_habitat_type")
sdsn
r <- rast(sdsn)
crs(r) <- "epsg:4326"
plot(r)
plot(flip(r))




# staging endpoint ----
stac_endpoint_url <- 'https://catalog.staging.edito.eu'

collects <- stac_endpoint_url %>%
  stac() %>%
  collections %>%
  get_request()

eunis <- collects$collections %>%
  keep(~ str_detect(tolower(.x$title), "eunis")) %>%
  .[[1]]

glimpse(eunis)

q <- stac_search(q = stac(stac_endpoint_url),
            collections = "emodnet-eunis_seabed_habitat_class_2019") %>%

  get_request()

href <- keep(eunis$links, ~ str_detect(.x$rel, 'child|item')) %>%
  .[[1]] %>%
  .[["href"]]
href
retryJson(href)


seab <- collects$collections %>%
  keep(~ str_detect(tolower(.x$title), "habitat")) %>%
  .[[6]]

seab



# endpoint_copernicus <- "https://s3.waw3-1.cloudferro.com/mdl-metadata/metadata/catalog.stac.json"
# cop_links <- retryJson(endpoint_copernicus)$links

# new endpoint ----
stac_endpoint_url <- 'https://catalog.dive.edito.eu/'

# q <- stac_search(q = stac(stac_endpoint_url), limit = 5000) |> get_request()
# 
# map_chr(q$features, \(x) x$collection)

#does not work - query for correct product / collection
# collects <- stac_endpoint_url %>%
#   stac() %>%
#   get_request()
# 
# eunis <- collects$collections %>%
#   keep(~ str_detect(tolower(.x$title), "seabed_energy")) %>%
#   .[[1]]

#request layer
q <- stac_search(q = stac(stac_endpoint_url),
                 collections = "emodnet-eunis_seabed_habitat_class_2019") %>%
  get_request()

#create link and extract raster
link <- rstac::assets_url(q, "Zarr")
#look at what data sets are in there
names(gdalinfo(link)$arrays)

#still uses EditoTools
dsn <- toGDAL(link)
sdsn=sprintf('%s:/%s',dsn,"marine_strategy_framework_directive_benthic_broad_habitat_type")
sdsn
r <- rast(sdsn)
crs(r) <- "epsg:4326"
plot(r)
plot(flip(r))

#still uses EditoTools
dsn <- toGDAL(link)
sdsn=sprintf('%s:/%s',dsn,"seabed_energy")
sdsn
r <- rast(sdsn)
crs(r) <- "epsg:4326"
plot(flip(r))

#COPERNICUS layer
# does not work --> how to find collection id?
q <- stac_search(q = stac(stac_endpoint_url),
                 collections = "climate_forecast-sea_water_speed") %>%
  get_request()

#create link and extract raster
link <- rstac::assets_url(q, "Zarr")
#look at what data sets are in there
names(gdalinfo(link)$arrays)

#still uses EditoTools
dsn <- toGDAL(link)
sdsn=sprintf('%s:/%s',dsn,"...")
sdsn
r <- rast(sdsn)
crs(r) <- "epsg:4326"
plot(r)
plot(flip(r))








library(stars)
r <- read_stars(sdsn,driver = "zarr")
plot(r)
# Extract values as a matrix
values_matrix <- as.matrix(r[[1]])

# Flip the matrix along the y-axis
flipped_matrix <- values_matrix[, ncol(values_matrix):1]

# Recreate the stars object with the flipped matrix
raster_flipped <- st_as_stars(flipped_matrix)

# Copy the dimensions from the original raster to the flipped raster
st_dimensions(raster_flipped) <- st_dimensions(r)

# Plot to check the result
plot(r, main = "Original Raster")
plot(raster_flipped, main = "Flipped Raster")

# Extract values as a matrix
values_matrix <- as.matrix(raster)

# Flip the matrix along the y-axis
flipped_matrix <- values_matrix[nrow(values_matrix):1, ]

# Ensure ncol and nrow return numeric values
ncol_raster <- ncol(raster)
nrow_raster <- nrow(raster)

# Create an empty raster with the same extent and resolution
raster_flipped <- rast(ncols=ncol_raster, nrows=nrow_raster, ext=ext(raster), crs=crs(raster))

# Assign the flipped values back to the new raster
values(raster_flipped) <- flipped_matrix





#try different categorical layer

stac_endpoint_url <- 'https://catalog.dive.edito.eu/'
collects <- stac_endpoint_url %>%
  stac() %>%
  get_request()

eunis <- collects$collections %>%
  keep(~ str_detect(tolower(.x$title), "analysis"))

q <- stac_search(q = stac(stac_endpoint_url),
                 collections = "sea_water_potential_temperature") %>%
  get_request()

link <- q$features[[1]]$assets$Zarr$href
dsn <- toGDAL(link)
sdsn=sprintf('%s:/%s',dsn,"eunis_seabed_habitat_class_2019")
sdsn
rast(sdsn)

s3 <- paws.storage::s3(credentials = list(anonymous = "logical"))
read_zarr_array(link, s3_client = s3)

url <- rstac::assets_url(q, "Zarr")
link = toGDAL(url, par = "eunis_seabed_habitat_class_2019")
rast(link)

#try copernicus

stac_endpoint_url <- 'https://catalog.dive.edito.eu/'
collects <- stac_endpoint_url %>%
  stac() %>%
  get_request()

eunis <- collects$collections %>%
  keep(~ str_detect(tolower(.x$title), "analysis"))

q <- stac_search(q = stac(stac_endpoint_url),
                 collections = "climate_forecast-temperature_of_analysis_of_sea_water") %>%
  get_request()


link <- q$features[[1]]$assets$Zarr$href
dsn <- toGDAL(link)
sdsn=sprintf('%s:/%s',dsn,"eunis_seabed_habitat_class_2019")
sdsn
rast(sdsn)

s3 <- paws.storage::s3(credentials = list(anonymous = "logical"))
read_zarr_array(link, s3_client = s3)


#### tests ----
b <- st_bbox(c(xmin = 16.1, xmax = 16.6, ymax = 48.6, ymin = 47.9))

q <- stac(stac_endpoint_url) %>%
  stac_search(collections = "emodnet-eunis_seabed_habitat_class_2019") %>%
  get_request()

lapply(q$features, \(x) names(x$properties))
ext_filter()


req <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") %>%
  stac_search(limit = 5)

t <- req %>% ext_filter(
  collection %in% c("landsat-c2-l2", "sentinel-2-l2a") &&
    datetime > "2019-01-01" &&
    datetime < "2019-06-01") %>%
  post_request()



req2 <- stac(stac_endpoint_url) %>%
  stac_search(bbox = c(-10, 10, 62, 68) ,
              datetime = "2020-01-01T00:00:00Z/2020-12-31T00:00:00Z",
              limit = 2000) %>% 
  get_request()

filtered_list <- req2$features %>%
  keep(~any(str_detect(.x$assets,"timeChunked.zarr")))

map(filtered_list, \(x) x$properties$`moi:short_name`)

filtered_list <- filtered_list %>%
  keep(~str_detect(.x$properties$`moi:short_name`,"vo"))


map_chr(filtered_list, "collection")

l <-filtered_list[[1]]$assets$`arco-geo-series`$href
dsn = toGDAL(l)
sdsn = sprintf('%s:/%s',dsn,"vo")
sdsn=sprintf('%s:%s',sdsn,1-1)
sdsn=sprintf('%s:%s',sdsn,50-1)

r = rast(sdsn)
plot(r)


req <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") %>%
  stac_search(limit = 5)
req2 <- stac(stac_endpoint_url, force_version = "1.0.0") %>%
  stac_search(limit = 5)

# Equal operator '=' with collection property
req %>% ext_filter(collection == "sentinel-2-l2a") %>% post_request()

req2 %>% ext_filter(collection == "dfff7b6b-c6dd-5178-bc30-986c057325d0") %>% post_request()
# req2 <- stac(stac_endpoint_url) %>%
#   stac_search() %>%
#   ext_filter(t_intersects(datetime, interval("1985-07-16T05:32:00Z",
#                                              "1985-07-24T16:50:35Z"))
#              ) %>%
#   post_request(encode = "form")



rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
  rstac::ext_filter(
    collection == "landsat-c2-l2" &&
      t_intersects(datetime, )
  )

t2 <- req2 %>% 
  ext_filter(collection == "emodnet-eunis_seabed_habitat_class_2019") %>% post_request()


stac(stac_endpoint_url) %>%
  ext_filter(collection == "emodnet-eunis_seabed_habitat_class_2019") %>% post_request()


req2 %>% 
  ext_filter(t_intersects(datetime, )) %>% post_request()

req2 %>% 
  stac_search(datetime %in% "2020-01-01/2020-01-31") %>% post_request()



### rstac tutorial ----
library("sf")
library("terra")

stac_endpoint_url <- 'https://catalog.dive.edito.eu/'
stac_source <- stac(stac_endpoint_url)
coll_q <- stac_source |>
  collections() |>
  get_request() 

stac_search(q = stac_source, 
            collections = "eunis_seabed_habitat_class_2019",
            limit = 999)
b <- st_bbox(c(xmin = 16.1, xmax = 16.6, ymax = 48.6, ymin = 47.9))
t <- stac_search(q = stac_source, 
            collections = "emodnet-eunis_seabed_habitat_class_2019",
            limit = 999,
            bbox = b) %>% get_request()

stac_search(q = stac_source, 
            collections = "emodnet-eunis_seabed_habitat_class_2019",
            limit = 999,
            bbox = b) %>% get_request()

signed_stac_query <- items_sign(t, sign_planetary_computer())
signed_stac_query
#does not work
assets_download(signed_stac_query, output_dir = "demonstrators/Ward/zarr_extraction/SAVE/")

h_url <- assets_url(t, "Zarr")
make_vsicurl_url <- function(base_url) {
  paste0(
    "/vsicurl", 
    "?pc_url_signing=yes",
    "&pc_collection=usgs-lcmap-conus-v13",
    "&url=",
    base_url
  )
}
lcpri_url <- make_vsicurl_url(h_url)



stac_source |>
  stac_search(datetime = "2000-01-01/2020-12-31") |> 
  ext_filter(collection == "dfff7b6b-c6dd-5178-bc30-986c057325d0") |>
  post_request()


stac_source |>
  stac_search() |>
  get_request()

stac_source |>
  stac_search(collections = "emodnet-eunis_seabed_habitat_class_2019") |>
  get_request() |>
  items_fields(field = "properties")
a <-  stac_source |>
  stac_search(collections = "emodnet-eunis_seabed_habitat_class_2019") |>
  get_request()
a <-  stac_source |>
  stac_search(collections = "emodnet-eunis_seabed_habitat_class_2019") |>
  get_request() |>
  items_filter(str_detect(properties$productIdentifier, "habitat"))

stac_source |>
  collections("emodnet-eunis_seabed_habitat_class_2019") |>
  queryables() |>
  get_request()


map(a$features, \(x) x$properties$productIdentifier)

# --> only queries the first 500 collections




stac_source |>
  stac_search(limit = ) |>
  get_request() 

a

s_obj <- stac("https://brazildatacube.dpi.inpe.br/stac/")
b <- s_obj |>
  stac_search(
    collections = "CB4_64_16D_STK-1",
    datetime = "2019-01-01/2019-12-31",
    limit = 100) |> 
  post_request()

b <- s_obj |> 
  stac_search(collections = c("CB4_64_16D_STK", "S2-16D-2")) |>
  get_request()

s_obj |>
  stac_search(
    collections = "CB4_64_16D_STK-1",
    datetime = "2019-01-01/2019-12-31",
    limit = 100) |> 
  post_request() |>
  items_fields(field = "properties")
s_obj |>
  stac_search(
    collections = c("CB4_64_16D_STK", "S2-16D-2")) |> 
  get_request() |>
  items_filter(properties$`eo:cloud_cover` < 10)


# |>
 #  items_filter(properties$`eo:cloud_cover` < 10)

