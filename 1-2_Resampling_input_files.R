library(ncdf4)
library(raster)
library(tidyverse)
library(terra)

# resample files to resolution of biological data (10 NM x 10 NM) ----
## netcdf files ----

setwd("/home/onyxia/work/BAR/DATA")

nc_chl <- ncdf4::nc_open("1.DOWNLOAD/environmental_variables/chlorophyll_concentration.nc")
v <- nc_chl$var[[1]]
size <- v$varsize
dims <- v$ndims
nt <- size[dims]                  # length of time dimension
lat <- nc_chl$dim$latitude$vals   # latitude position
lon <- nc_chl$dim$longitude$vals  # longitude position
t_index <- tibble(index = seq(1:size[4]),
                  year = rep(c(2000:2020), each = 12),
                  month = rep(seq(1:12),21))

# read chl variable
r<-list()

for (i in 1:nt) {
  start <- rep(1,dims)         # begin with start=(1,1,...,1)
  start[dims] <- i             # change to start=(1,1,...,i) to read timestep i
  count <- size                # begin with count=(nx,ny,...,nt), reads entire var
  count[dims] <- 1             # change to count=(nx,ny,...,1) to read 1 tstep
  
  dt <- ncvar_get(nc_chl, varid = 'chl', start = start, count = count)
  # convert to raster
  # transpose the raster to have correct orientation
  r[i]<-t(flip(raster(dt), direction = "x"))
}

# create layer stack with time dimension
st_chl <-stack(r)
extent(st_chl) <- extent(c(range(lon), range(lat)))
crs(st_chl) <- CRS("+init=epsg:4326")

# remove edges of raster that have NA values
r_ex <- st_chl[[1]]

r1NaM <- is.na(as.matrix(r_ex))
colNotNA <- which(colSums(r1NaM) != nrow(r_ex))
rowNotNA <- which(rowSums(r1NaM) != ncol(r_ex))
r_extent <- extent(r_ex, 
                   rowNotNA[1], rowNotNA[length(rowNotNA)],
                   colNotNA[1], colNotNA[length(colNotNA)])
r_ex <- crop(r_ex, r_extent)

layout(matrix(1:2, nrow=1))
plot(st_chl[[1]])
plot(r_ex)

# change spatial resolution to 10 x 10 NM (because this is the resolution of the biological data)
r_proj <- projectRaster(r_ex, crs = crs("+init=epsg:3034"))
# 1 nautical mile is 1852 m, so 10 NM is 18520 m
r_empty <- raster(nrows = 1000, ncols = 1000, xmn = extent(r_proj)[1], 
                       xmx = extent(r_proj)[1] + 18520 * 1000, ymn = extent(r_proj)[3], ymx = extent(r_proj)[3] + 18520 * 1000, crs = "+init=epsg:3034")
r_empty <- raster(xmn = extent(r_proj)[1], xmx = extent(r_proj)[2], 
                  ymn = extent(r_proj)[3], ymx = extent(r_proj)[4],
                  res = 18520,
                  crs = "+init=epsg:3034")

r_res <- resample(r_proj, r_empty, method = "bilinear")
r_res <- projectRaster(r_res, crs = crs("+init=epsg:4326"))

r1NaM <- is.na(as.matrix(r_res))
colNotNA <- which(colSums(r1NaM) != nrow(r_res))
rowNotNA <- which(rowSums(r1NaM) != ncol(r_res))
r_extent <- extent(r_res, 
                   rowNotNA[1], rowNotNA[length(rowNotNA)],
                   colNotNA[1], colNotNA[length(colNotNA)])
r_res <- crop(r_res, r_extent)


desired_nc_vars <- c("thetao", "so", "max_v", "o2", "chl", "phyc", "zooc", "zeu")
key_nc_vars <- c("SST", "SSS", "max_SSV", "DO", "Chl", "Phyto", "ZooPl", "EuphD")

#resample netcdf files
fls <- list.files("1.DOWNLOAD/environmental_variables/", pattern = ".nc$", full.names = "TRUE")
for (f in fls) {
  nc <- ncdf4::nc_open(f)
  nc_v <- which(names(nc$var) %in% desired_nc_vars)
  
  for (v_n in nc_v) {
    v <- nc$var[[v_n]]
    v_name <- v$name
    v_name_out <- key_nc_vars[which(desired_nc_vars == v_name)]
    size <- v$varsize
    dims <- v$ndims
    nt <- size[dims]       # length of time dimension
    lat <- nc$dim$lat$vals # latitude position
    lon <- nc$dim$lon$vals # longitude position
    
    # read the variable
    r <-list()
    
    for (i in 1:nt) {
      start <- rep(1,dims)         # begin with start=(1,1,...,1)
      start[dims] <- i             # change to start=(1,1,...,i) to read timestep i
      count <- size                # begin with count=(nx,ny,...,nt), reads entire var
      count[dims] <- 1             # change to count=(nx,ny,...,1) to read 1 tstep
      
      dt <- ncvar_get(nc, varid = v_name, start = start, count = count)
      # convert to raster
      # transpose the raster to have correct orientation
      r[i]<-t(flip(raster(dt), direction = "x"))
    }
    
    # create layer stack with time dimension
    st<-stack(r)
    extent(st) <- extent(c(range(lon), range(lat)))
    crs(st) <- CRS("+init=epsg:4326")

    #reproject to coarse resolution
    st_res <- resample(st, r_res)
    names(st_res) <- paste(t_index$year, t_index$month, sep = "_")
    
    raster::writeRaster(st_res, paste("2.PREPROCESSED/environmental_variables/",v_name_out, ".grd", sep = ""), format = "raster", overwrite = TRUE)

    print(paste("processed", v_name))
  }
}

## tif files (numerical - bilinear) ----
fls <- list.files("1.DOWNLOAD/environmental_variables/", pattern = ".tif$", full.names = "TRUE")
for (f in fls[1]) {
  r_f <- raster(f)
  r_name <- f %>% str_remove(".tif$") %>% str_remove("1.DOWNLOAD/environmental_variables/")
  
  r_res <- resample(r_f, r_res)
  
  raster::writeRaster(r_res, paste("2.PREPROCESSED/environmental_variables/",r_name, ".grd", sep = ""), format = "raster", overwrite = TRUE)

  print(paste("processed", r_name))
}

## tif files (categorical - nearest neighbour) ----
fls <- list.files("1.DOWNLOAD/environmental_variables/", pattern = ".tif$", full.names = "TRUE")
for (f in fls[c(2,3)]) {
  r_f <- raster(f)
  r_name <- f %>% str_remove(".tif$") %>% str_remove("1.DOWNLOAD/environmental_variables/")
  
  r_res <- resample(r_f, r_res, method = "ngb")
  
  raster::writeRaster(r_res, paste("2.PREPROCESSED/environmental_variables/", r_name, ".grd", sep = ""), format = "raster", overwrite = TRUE)

  print(paste("processed", r_name))
}

unique(values(r_res))

#resample grd - windfarms
fls <- list.files("1.DOWNLOAD/environmental_variables/", pattern = ".grd$", full.names = "TRUE")
for (f in fls) {
  st_f <- stack(f)
  st_name <- f %>% str_remove(".grd$") %>% str_remove("1.DOWNLOAD/environmental_variables/")
  
  st_res <- resample(st_f, r_res, method = "ngb")
  
  raster::writeRaster(st_res, paste("2.PREPROCESSED/environmental_variables/", st_name, ".grd", sep = ""), format = "raster", overwrite = TRUE)

  print(paste("processed", st_name))
}
plot(st_res$X2020_12)
