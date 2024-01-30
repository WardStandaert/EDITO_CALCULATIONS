library(ncdf4)
library(lubridate)
library(tictoc)
library(dplyr)

# ONLY WORKS FOR R VERSION OF 3.6.1 OR LOWER !!
#also need installation of older version of Rtools
USERNAME <- "wstandaert1"    # input your username
PASSWORD <- "FPrYpgCygTr4r2" # input your password

server = "@my.cmems-du.eu";                # Copernicus Marine server
datasetID = "cmems_mod_glo_phy_my_0.083_P1D-m";  # datasetID

# creates the OPeNDAP url
url <- paste ("https://",USERNAME, ":", PASSWORD,"@",server,"/thredds/dodsC/",datasetID, sep = "")

# Open the connection
ds <- nc_open(url)   
print(ds)


## Get dataset's dimensions 
# Longitude
lon <- ncvar_get(ds, "longitude")
nlon <- dim(lon)
# Latitude
lat <- ncvar_get(ds, "latitude")
nlat <- dim(lat)
# Depth
depth <- ncvar_get(ds, "depth")
ndepth <- dim(depth)

# Check dimensions
print(c(nlon,nlat,ndepth))

# Time
time<-ncvar_get(ds,"time")
nt <- dim(time)
t_units <- ncatt_get(ds, "time", "units")
t_units
# convert time -- split the time units string into fields
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) + dhours(time)      
date

## Define the parameters and ranges for subset
#Bounding box
x <- c(-20, 30)          # longitude
y <- c(40, 65)           # latitude
z <- c(9, 10)

head(cbind(lon, lat, depth), n = 10)

# Function to get the indices from the ranges
btw <- function(data, num){
  c(min(which(between(data,num[1],num[2]))),
    max(which(between(data,num[1],num[2]))))
}

# Starting indices
lon_indices <- btw(data = lon, num = x)
lat_indices <- btw(data = lat, num = y)
depth_indices <- btw(data = depth, num = z)

# Count number of indices to extract along each dimension
lon_range <- lon_indices[-1] - lon_indices[1]+1
lat_range <- lat_indices[-1] - lat_indices[1]+1
depth_range <- depth_indices[-1] - depth_indices[1]+1

# Start and Count vectors
offset <- c(lon_indices[1], lat_indices[1], depth_indices[1])    #lon,lat,depth
count <- c(lon_range, lat_range, depth_range)                    #dimension var


##### download per year, aggregate monthly and combine north and east #####
for (yr in 2000:2020) {
  tic(yr)
  
  t_start <- ymd(paste(yr, "-01-01", sep = ""))
  t_end <- ymd(paste(yr, "-12-31", sep = ""))
  
  t_start_ind <- which(year(date) == yr & month(date) == 1 &  day(date) == 1)
  t_end_ind <- which(year(date) == yr & month(date) == 12 &  day(date) == 31)
  
  t_count_ind <- t_end_ind - t_start_ind + 1
  # Get subsetted variables   
  var_uo <- ncvar_get(ds,"uo", start = c(offset,t_start_ind), count = c(count,t_count_ind))
  var_vo <- ncvar_get(ds,"vo", start = c(offset,t_start_ind), count = c(count,t_count_ind))
  
  time_yr <- ncvar_get(ds,"time", start = t_start_ind, count = t_count_ind)
  date_yr <- ymd(t_dstr) + dhours(time_yr)
  
  for (m in 1:12) {
    #subset per month
    d1 <- which(date(date_yr) == paste(yr, m, "01", sep = "-"))
    days_in_month <- sum(month(date_yr) == m)
    d2 <- d1 - 1 + days_in_month
    
    var_uo_m <- abs(var_uo[,,d1:d2])
    var_vo_m <- abs(var_vo[,,d1:d2])
    out_uo <- apply(var_uo_m, c(1,2), max, na.rm = T)  ## we calculate the maximum current
    out_vo <- apply(var_vo_m, c(1,2), max, na.rm = T)  ## we calculate the maximum current
    out <- pmax(out_uo, out_vo, na.rm = T)
    print(m)
    assign(paste("values", yr, m, sep = "-"), out)
  }
  toc()
}

# Function to select a range from a list
select_range <- function(lst, indices) {
  sub_lst <- lst[[1]]  # Get the list
  return(sub_lst[indices[1]:indices[2]])
}

##create output NetCDF file with monthly aggregates
lat_list <- list(ds[["dim"]][["latitude"]][["vals"]])
lon_list <- list(ds[["dim"]][["longitude"]][["vals"]])
time_list <- list(ds[["dim"]][["time"]][["vals"]])

# Extract the dimension range using the indices
range_lat <- select_range(lat_list, lat_indices)          #lat_indices
range_lon <- select_range(lon_list, lon_indices)          #lon_indices
range_time <- as.numeric(difftime(seq(dmy("01/01/2000"),  #define hours since 1950-01-01
                                      dmy("31/12/2020"), 
                                      by = "month"), 
                                  ymd("1950-01-01"), 
                                  unit = "hours"))

currents_array <- array(c(`values-2000-1`, `values-2000-2`, `values-2000-3`, `values-2000-4`, `values-2000-5`, `values-2000-6`, 
                        `values-2000-7`, `values-2000-8`, `values-2000-9`, `values-2000-10`, `values-2000-11`, `values-2000-12`, 
                        `values-2001-1`, `values-2001-2`, `values-2001-3`, `values-2001-4`, `values-2001-5`, `values-2001-6`, 
                        `values-2001-7`, `values-2001-8`, `values-2001-9`, `values-2001-10`, `values-2001-11`, `values-2001-12`,
                        `values-2002-1`, `values-2002-2`, `values-2002-3`, `values-2002-4`, `values-2002-5`, `values-2002-6`, 
                        `values-2002-7`, `values-2002-8`, `values-2002-9`, `values-2002-10`, `values-2002-11`, `values-2002-12`,
                        `values-2003-1`, `values-2003-2`, `values-2003-3`, `values-2003-4`, `values-2003-5`, `values-2003-6`, 
                        `values-2003-7`, `values-2003-8`, `values-2003-9`, `values-2003-10`, `values-2003-11`, `values-2003-12`,
                        `values-2004-1`, `values-2004-2`, `values-2004-3`, `values-2004-4`, `values-2004-5`, `values-2004-6`, 
                        `values-2004-7`, `values-2004-8`, `values-2004-9`, `values-2004-10`, `values-2004-11`, `values-2004-12`,
                        `values-2005-1`, `values-2005-2`, `values-2005-3`, `values-2005-4`, `values-2005-5`, `values-2005-6`, 
                        `values-2005-7`, `values-2005-8`, `values-2005-9`, `values-2005-10`, `values-2005-11`, `values-2005-12`,
                        `values-2006-1`, `values-2006-2`, `values-2006-3`, `values-2006-4`, `values-2006-5`, `values-2006-6`, 
                        `values-2006-7`, `values-2006-8`, `values-2006-9`, `values-2006-10`, `values-2006-11`, `values-2006-12`,
                        `values-2007-1`, `values-2007-2`, `values-2007-3`, `values-2007-4`, `values-2007-5`, `values-2007-6`, 
                        `values-2007-7`, `values-2007-8`, `values-2007-9`, `values-2007-10`, `values-2007-11`, `values-2007-12`,
                        `values-2008-1`, `values-2008-2`, `values-2008-3`, `values-2008-4`, `values-2008-5`, `values-2008-6`, 
                        `values-2008-7`, `values-2008-8`, `values-2008-9`, `values-2008-10`, `values-2008-11`, `values-2008-12`,
                        `values-2009-1`, `values-2009-2`, `values-2009-3`, `values-2009-4`, `values-2009-5`, `values-2009-6`, 
                        `values-2009-7`, `values-2009-8`, `values-2009-9`, `values-2009-10`, `values-2009-11`, `values-2009-12`,
                        `values-2010-1`, `values-2010-2`, `values-2010-3`, `values-2010-4`, `values-2010-5`, `values-2010-6`, 
                        `values-2010-7`, `values-2010-8`, `values-2010-9`, `values-2010-10`, `values-2010-11`, `values-2010-12`,
                        `values-2011-1`, `values-2011-2`, `values-2011-3`, `values-2011-4`, `values-2011-5`, `values-2011-6`, 
                        `values-2011-7`, `values-2011-8`, `values-2011-9`, `values-2011-10`, `values-2011-11`, `values-2011-12`,
                        `values-2012-1`, `values-2012-2`, `values-2012-3`, `values-2012-4`, `values-2012-5`, `values-2012-6`, 
                        `values-2012-7`, `values-2012-8`, `values-2012-9`, `values-2012-10`, `values-2012-11`, `values-2012-12`,
                        `values-2013-1`, `values-2013-2`, `values-2013-3`, `values-2013-4`, `values-2013-5`, `values-2013-6`, 
                        `values-2013-7`, `values-2013-8`, `values-2013-9`, `values-2013-10`, `values-2013-11`, `values-2013-12`,
                        `values-2014-1`, `values-2014-2`, `values-2014-3`, `values-2014-4`, `values-2014-5`, `values-2014-6`, 
                        `values-2014-7`, `values-2014-8`, `values-2014-9`, `values-2014-10`, `values-2014-11`, `values-2014-12`,
                        `values-2015-1`, `values-2015-2`, `values-2015-3`, `values-2015-4`, `values-2015-5`, `values-2015-6`, 
                        `values-2015-7`, `values-2015-8`, `values-2015-9`, `values-2015-10`, `values-2015-11`, `values-2015-12`,
                        `values-2016-1`, `values-2016-2`, `values-2016-3`, `values-2016-4`, `values-2016-5`, `values-2016-6`, 
                        `values-2016-7`, `values-2016-8`, `values-2016-9`, `values-2016-10`, `values-2016-11`, `values-2016-12`,
                        `values-2017-1`, `values-2017-2`, `values-2017-3`, `values-2017-4`, `values-2017-5`, `values-2017-6`, 
                        `values-2017-7`, `values-2017-8`, `values-2017-9`, `values-2017-10`, `values-2017-11`, `values-2017-12`,
                        `values-2018-1`, `values-2018-2`, `values-2018-3`, `values-2018-4`, `values-2018-5`, `values-2018-6`, 
                        `values-2018-7`, `values-2018-8`, `values-2018-9`, `values-2018-10`, `values-2018-11`, `values-2018-12`,
                        `values-2019-1`, `values-2019-2`, `values-2019-3`, `values-2019-4`, `values-2019-5`, `values-2019-6`, 
                        `values-2019-7`, `values-2019-8`, `values-2019-9`, `values-2019-10`, `values-2019-11`, `values-2019-12`,
                        `values-2020-1`, `values-2020-2`, `values-2020-3`, `values-2020-4`, `values-2020-5`, `values-2020-6`, 
                        `values-2020-7`, `values-2020-8`, `values-2020-9`, `values-2020-10`, `values-2020-11`, `values-2020-12`),
                      dim = c(length(range_lon),length(range_lat),length(range_time)))

currents_array[which(currents_array == - Inf)] <- -32767

# Define the dimensions
dim_lon <- ncdim_def("lon", "degrees_east", range_lon)
dim_lat <- ncdim_def("lat", "degrees_north", range_lat)
dim_time <- ncdim_def("time", "hours since 1950-01-01 00:00:00", range_time)

# Define the dimensions of the variables Zeu and Zooc
count <- c(lon_range, lat_range, length(range_time)) #dimension Zeu & Zooc

var_currents <- ncvar_def("max_v", "m/s", list(dim_lon, dim_lat, dim_time), -32767, longname = "maximum sea water velocity", prec = "float")

setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Input-output_files/0.Preprocessed files")
# Create a new NetCDF file object
ncnew_currents <- nc_create("max_currents.nc", var_currents)
ncvar_put(ncnew_currents, var_currents, currents_array, start = NA, count = NA)
nc_close(ncnew_currents)
