library(rstac)
library(purrr)
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

eunis$links %>%
  keep(~ .x)
c('child','item')

href <- keep(eunis$links, ~ str_detect(.x$rel, 'child|item')) %>%
  .[[1]] %>%
  .[["href"]] %>% items_fetch()
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
collects <- stac_endpoint_url %>%
  stac() %>%
  get_request()

eunis <- collects$links %>%
  keep(~ str_detect(tolower(.x$title), "thetao")) %>%
  .[[1]]

q <- stac_search(q = stac(stac_endpoint_url),
                 collections = "emodnet-eunis_seabed_habitat_class_2019") %>%
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