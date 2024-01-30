library(ncdf4)
library(lubridate)
library(tictoc)
library(dplyr)

# ONLY WORKS FOR R VERSION OF 3.6.1 OR LOWER !!
#also need installation of older version of Rtools
USERNAME <- "wstandaert1"    # input your username
PASSWORD <- "FPrYpgCygTr4r2" # input your password

server = "@my.cmems-du.eu";                # Copernicus Marine server
datasetID = "cmems_mod_glo_phy_my_0.083deg_P1D_m";  # datasetID

# creates the OPeNDAP url
url <- paste ("https://",USERNAME, ":", PASSWORD,"@",server,"/thredds/dodsC/",datasetID, sep = "")

# Open the connection
ds <- nc_open(url)   
print(ds)


## Get dataset's dimensions 
# Longitude
lon <- ncvar_get(ds, "lon")
nlon <- dim(lon)
# Latitude
lat <- ncvar_get(ds, "lat")
nlat <- dim(lat)
# Depth
depth <- ncvar_get(ds, "depth")
ndepth <- dim(depth)

# Check dimensions
print(c(nlon,nlat))

# Time
time<-ncvar_get(ds,"time")
nt <- dim(time)
t_units <- ncatt_get(ds, "time", "units")
t_units
# convert time -- split the time units string into fields
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) + ddays(time)      
date


## Define the parameters and ranges for subset
#Bounding box
x <- c(-5.77, 31.15)               # longitude
y <- c(52.2997, 66.7372)                # latitude
z <- c(9, 10)

# Function to get the indices from the ranges
btw <- function(data, num){
  c(min(which(num<=data)), max(which(num>=data)))
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
offset <- c(lon_indices[1], lat_indices[1], depth_indices[1])    #lon,lat,depth,time
count <- c(lon_range, lat_range, depth_range) #dimension CHL


##### download 2000-2020 directly, without aggregation #####
t_start_ind <- which(year(date) == 2000 & month(date) == 1)
t_end_ind <- which(year(date) == 2020 & month(date) == 12)
t_count_ind <- t_end_ind - t_start_ind + 1

var_o2 <- ncvar_get(ds,"o2", start = c(offset,t_start_ind), count = c(count,t_count_ind))
var_chl <- ncvar_get(ds,"chl", start = c(offset,t_start_ind), count = c(count,t_count_ind))

time_yr <- ncvar_get(ds,"time", start = t_start_ind, count = t_count_ind)
date_yr <- ymd(t_dstr) + ddays(time_yr)

##create output NetCDF file with monthly aggregates
lat_list <- list(ds[["dim"]][["lat"]][["vals"]])
lon_list <- list(ds[["dim"]][["lon"]][["vals"]])
time_list <- list(ds[["dim"]][["time"]][["vals"]])

# Function to select a range from a list
select_range <- function(lst, indices) {
  sub_lst <- lst[[1]]  # Get the list
  return(sub_lst[indices[1]:indices[2]])
}

# Extract the dimension range using the indices
range_lat <- select_range(lat_list, lat_indices)    #lat_indices
range_lon <- select_range(lon_list, lon_indices)   #lon_indices
range_time <- as.numeric(difftime(seq(dmy("01/01/2000"),  #define hours since 1950-01-01
                                      dmy("31/12/2020"), 
                                      by = "month"), 
                                  ymd("1950-01-01"), 
                                  unit = "hours"))

# Define the dimensions of the variables
dimo2  <- ncdim_def("o2", "mmol m-3", count)
dimchl <- ncdim_def("chl", "mg m-3", count)
dim_lon <- ncdim_def("longitude", "degrees_east", range_lon)
dim_lat <- ncdim_def("latitude", "degrees_north", range_lat)
dim_time <- ncdim_def("time", "days since 1900-01-01 00:00:00", range_time)


v_o2 <- ncvar_def("o2", "mmol m-3", list(dim_lon, dim_lat, dim_time), -999, longname = "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water", prec = "float")
v_chl <- ncvar_def("chl", "mg m-3", list(dim_lon, dim_lat, dim_time), -999, longname = "mass_concentration_of_chlorophyll_a_in_sea_water", prec = "float")

# Create a new NetCDF file object
ncnew_o2<- nc_create("dissolved_oxygen_Baltic.nc", v_o2)
ncvar_put(ncnew_o2, v_o2, var_o2, start = NA, count = NA)
nc_close(ncnew_o2)

ncnew_chl <- nc_create("chlorophyll_Baltic.nc", v_chl)
ncvar_put(ncnew_chl, v_chl, var_chl, start = NA, count = NA)
nc_close(ncnew_chl)




##### download per year and aggregate monthly #####
for (yr in 2000:2020) {
  tic(yr)
  
  t_start <- ymd(paste(yr, "-01-01", sep = ""))
  t_end <- ymd(paste(yr, "-12-31", sep = ""))
  
  t_start_ind <- which(year(date) == yr & month(date) == 1 &  day(date) == 1)
  t_end_ind <- which(year(date) == yr & month(date) == 12 &  day(date) == 31)
  
  t_count_ind <- t_end_ind - t_start_ind + 1
  # Get subsetted variable   
  var <- ncvar_get(ds,"chl", start = c(offset,t_start_ind), count = c(count,t_count_ind))
  time_yr <- ncvar_get(ds,"time", start = t_start_ind, count = t_count_ind)
  date_yr <- ymd(t_dstr) + dseconds(time_yr)

  for (m in 1:12) {
    #subset per month
    d1 <- yday(paste(yr, m, "01", sep = "-"))
    d2 <- d1 - 1 + days_in_month(paste(yr, m, "01", sep = "-"))
    
    var_m <- var[,,d1:d2]
    out <- apply(var_m, c(1,2), mean, na.rm = TRUE)
    print(m)
    assign(paste("zooc_values", yr, m, sep = "-"), out)
  }
  toc()
}

#download per year and aggregate monthly
for (yr in 2000:2020) {
  tic(yr)
  
  t_start <- ymd(paste(yr, "-01-01", sep = ""))
  t_end <- ymd(paste(yr, "-12-31", sep = ""))
  
  t_start_ind <- which(year(date) == yr & month(date) == 1 &  day(date) == 1)
  t_end_ind <- which(year(date) == yr & month(date) == 12 &  day(date) == 31)
  
  t_count_ind <- t_end_ind - t_start_ind + 1
  # Get subsetted variable   
  var <- ncvar_get(ds,"zeu", start = c(offset,t_start_ind), count = c(count,t_count_ind))
  time_yr <- ncvar_get(ds,"time", start = t_start_ind, count = t_count_ind)
  date_yr <- ymd(t_dstr) + dseconds(time_yr)
  
  for (m in 1:12) {
    #subset per month
    d1 <- yday(paste(yr, m, "01", sep = "-"))
    d2 <- d1 - 1 + days_in_month(paste(yr, m, "01", sep = "-"))
    
    var_m <- var[,,d1:d2]
    out <- apply(var_m, c(1,2), mean, na.rm = TRUE)
    print(m)
    assign(paste("zeu_values", yr, m, sep = "-"), out)
  }
  toc()
}

# Function to select a range from a list
select_range <- function(lst, indices) {
  sub_lst <- lst[[1]]  # Get the list
  return(sub_lst[indices[1]:indices[2]])
}

# Extract the dimension range using the indices
range_lat <- select_range(lat_list, lat_indices)    #lat_indices
range_lon <- select_range(lon_list, lon_indices)   #lon_indices
range_time <- as.numeric(difftime(seq(dmy("01/01/2000"),  #define hours since 1950-01-01
                                      dmy("31/12/2020"), 
                                      by = "month"), 
                                  ymd("1950-01-01"), 
                                  unit = "hours"))

zooc_array <- array(c(`zooc_values-2000-1`, `zooc_values-2000-2`, `zooc_values-2000-3`, `zooc_values-2000-4`, `zooc_values-2000-5`, `zooc_values-2000-6`, 
                             `zooc_values-2000-7`, `zooc_values-2000-8`, `zooc_values-2000-9`, `zooc_values-2000-10`, `zooc_values-2000-11`, `zooc_values-2000-12`, 
                             `zooc_values-2001-1`, `zooc_values-2001-2`, `zooc_values-2001-3`, `zooc_values-2001-4`, `zooc_values-2001-5`, `zooc_values-2001-6`, 
                             `zooc_values-2001-7`, `zooc_values-2001-8`, `zooc_values-2001-9`, `zooc_values-2001-10`, `zooc_values-2001-11`, `zooc_values-2001-12`,
                             `zooc_values-2002-1`, `zooc_values-2002-2`, `zooc_values-2002-3`, `zooc_values-2002-4`, `zooc_values-2002-5`, `zooc_values-2002-6`, 
                             `zooc_values-2002-7`, `zooc_values-2002-8`, `zooc_values-2002-9`, `zooc_values-2002-10`, `zooc_values-2002-11`, `zooc_values-2002-12`,
                             `zooc_values-2003-1`, `zooc_values-2003-2`, `zooc_values-2003-3`, `zooc_values-2003-4`, `zooc_values-2003-5`, `zooc_values-2003-6`, 
                             `zooc_values-2003-7`, `zooc_values-2003-8`, `zooc_values-2003-9`, `zooc_values-2003-10`, `zooc_values-2003-11`, `zooc_values-2003-12`,
                             `zooc_values-2004-1`, `zooc_values-2004-2`, `zooc_values-2004-3`, `zooc_values-2004-4`, `zooc_values-2004-5`, `zooc_values-2004-6`, 
                             `zooc_values-2004-7`, `zooc_values-2004-8`, `zooc_values-2004-9`, `zooc_values-2004-10`, `zooc_values-2004-11`, `zooc_values-2004-12`,
                             `zooc_values-2005-1`, `zooc_values-2005-2`, `zooc_values-2005-3`, `zooc_values-2005-4`, `zooc_values-2005-5`, `zooc_values-2005-6`, 
                             `zooc_values-2005-7`, `zooc_values-2005-8`, `zooc_values-2005-9`, `zooc_values-2005-10`, `zooc_values-2005-11`, `zooc_values-2005-12`,
                             `zooc_values-2006-1`, `zooc_values-2006-2`, `zooc_values-2006-3`, `zooc_values-2006-4`, `zooc_values-2006-5`, `zooc_values-2006-6`, 
                             `zooc_values-2006-7`, `zooc_values-2006-8`, `zooc_values-2006-9`, `zooc_values-2006-10`, `zooc_values-2006-11`, `zooc_values-2006-12`,
                             `zooc_values-2007-1`, `zooc_values-2007-2`, `zooc_values-2007-3`, `zooc_values-2007-4`, `zooc_values-2007-5`, `zooc_values-2007-6`, 
                             `zooc_values-2007-7`, `zooc_values-2007-8`, `zooc_values-2007-9`, `zooc_values-2007-10`, `zooc_values-2007-11`, `zooc_values-2007-12`,
                             `zooc_values-2008-1`, `zooc_values-2008-2`, `zooc_values-2008-3`, `zooc_values-2008-4`, `zooc_values-2008-5`, `zooc_values-2008-6`, 
                             `zooc_values-2008-7`, `zooc_values-2008-8`, `zooc_values-2008-9`, `zooc_values-2008-10`, `zooc_values-2008-11`, `zooc_values-2008-12`,
                             `zooc_values-2009-1`, `zooc_values-2009-2`, `zooc_values-2009-3`, `zooc_values-2009-4`, `zooc_values-2009-5`, `zooc_values-2009-6`, 
                             `zooc_values-2009-7`, `zooc_values-2009-8`, `zooc_values-2009-9`, `zooc_values-2009-10`, `zooc_values-2009-11`, `zooc_values-2009-12`,
                             `zooc_values-2010-1`, `zooc_values-2010-2`, `zooc_values-2010-3`, `zooc_values-2010-4`, `zooc_values-2010-5`, `zooc_values-2010-6`, 
                             `zooc_values-2010-7`, `zooc_values-2010-8`, `zooc_values-2010-9`, `zooc_values-2010-10`, `zooc_values-2010-11`, `zooc_values-2010-12`,
                             `zooc_values-2011-1`, `zooc_values-2011-2`, `zooc_values-2011-3`, `zooc_values-2011-4`, `zooc_values-2011-5`, `zooc_values-2011-6`, 
                             `zooc_values-2011-7`, `zooc_values-2011-8`, `zooc_values-2011-9`, `zooc_values-2011-10`, `zooc_values-2011-11`, `zooc_values-2011-12`,
                             `zooc_values-2012-1`, `zooc_values-2012-2`, `zooc_values-2012-3`, `zooc_values-2012-4`, `zooc_values-2012-5`, `zooc_values-2012-6`, 
                             `zooc_values-2012-7`, `zooc_values-2012-8`, `zooc_values-2012-9`, `zooc_values-2012-10`, `zooc_values-2012-11`, `zooc_values-2012-12`,
                             `zooc_values-2013-1`, `zooc_values-2013-2`, `zooc_values-2013-3`, `zooc_values-2013-4`, `zooc_values-2013-5`, `zooc_values-2013-6`, 
                             `zooc_values-2013-7`, `zooc_values-2013-8`, `zooc_values-2013-9`, `zooc_values-2013-10`, `zooc_values-2013-11`, `zooc_values-2013-12`,
                             `zooc_values-2014-1`, `zooc_values-2014-2`, `zooc_values-2014-3`, `zooc_values-2014-4`, `zooc_values-2014-5`, `zooc_values-2014-6`, 
                             `zooc_values-2014-7`, `zooc_values-2014-8`, `zooc_values-2014-9`, `zooc_values-2014-10`, `zooc_values-2014-11`, `zooc_values-2014-12`,
                             `zooc_values-2015-1`, `zooc_values-2015-2`, `zooc_values-2015-3`, `zooc_values-2015-4`, `zooc_values-2015-5`, `zooc_values-2015-6`, 
                             `zooc_values-2015-7`, `zooc_values-2015-8`, `zooc_values-2015-9`, `zooc_values-2015-10`, `zooc_values-2015-11`, `zooc_values-2015-12`,
                             `zooc_values-2016-1`, `zooc_values-2016-2`, `zooc_values-2016-3`, `zooc_values-2016-4`, `zooc_values-2016-5`, `zooc_values-2016-6`, 
                             `zooc_values-2016-7`, `zooc_values-2016-8`, `zooc_values-2016-9`, `zooc_values-2016-10`, `zooc_values-2016-11`, `zooc_values-2016-12`,
                             `zooc_values-2017-1`, `zooc_values-2017-2`, `zooc_values-2017-3`, `zooc_values-2017-4`, `zooc_values-2017-5`, `zooc_values-2017-6`, 
                             `zooc_values-2017-7`, `zooc_values-2017-8`, `zooc_values-2017-9`, `zooc_values-2017-10`, `zooc_values-2017-11`, `zooc_values-2017-12`,
                             `zooc_values-2018-1`, `zooc_values-2018-2`, `zooc_values-2018-3`, `zooc_values-2018-4`, `zooc_values-2018-5`, `zooc_values-2018-6`, 
                             `zooc_values-2018-7`, `zooc_values-2018-8`, `zooc_values-2018-9`, `zooc_values-2018-10`, `zooc_values-2018-11`, `zooc_values-2018-12`,
                             `zooc_values-2019-1`, `zooc_values-2019-2`, `zooc_values-2019-3`, `zooc_values-2019-4`, `zooc_values-2019-5`, `zooc_values-2019-6`, 
                             `zooc_values-2019-7`, `zooc_values-2019-8`, `zooc_values-2019-9`, `zooc_values-2019-10`, `zooc_values-2019-11`, `zooc_values-2019-12`,
                             `zooc_values-2020-1`, `zooc_values-2020-2`, `zooc_values-2020-3`, `zooc_values-2020-4`, `zooc_values-2020-5`, `zooc_values-2020-6`, 
                             `zooc_values-2020-7`, `zooc_values-2020-8`, `zooc_values-2020-9`, `zooc_values-2020-10`, `zooc_values-2020-11`, `zooc_values-2020-12`),
                    dim = c(length(range_lon),length(range_lat),length(range_time)))

zeu_array <- array(c(`zeu_values-2000-1`, `zeu_values-2000-2`, `zeu_values-2000-3`, `zeu_values-2000-4`, `zeu_values-2000-5`, `zeu_values-2000-6`, 
                      `zeu_values-2000-7`, `zeu_values-2000-8`, `zeu_values-2000-9`, `zeu_values-2000-10`, `zeu_values-2000-11`, `zeu_values-2000-12`, 
                      `zeu_values-2001-1`, `zeu_values-2001-2`, `zeu_values-2001-3`, `zeu_values-2001-4`, `zeu_values-2001-5`, `zeu_values-2001-6`, 
                      `zeu_values-2001-7`, `zeu_values-2001-8`, `zeu_values-2001-9`, `zeu_values-2001-10`, `zeu_values-2001-11`, `zeu_values-2001-12`,
                      `zeu_values-2002-1`, `zeu_values-2002-2`, `zeu_values-2002-3`, `zeu_values-2002-4`, `zeu_values-2002-5`, `zeu_values-2002-6`, 
                      `zeu_values-2002-7`, `zeu_values-2002-8`, `zeu_values-2002-9`, `zeu_values-2002-10`, `zeu_values-2002-11`, `zeu_values-2002-12`,
                      `zeu_values-2003-1`, `zeu_values-2003-2`, `zeu_values-2003-3`, `zeu_values-2003-4`, `zeu_values-2003-5`, `zeu_values-2003-6`, 
                      `zeu_values-2003-7`, `zeu_values-2003-8`, `zeu_values-2003-9`, `zeu_values-2003-10`, `zeu_values-2003-11`, `zeu_values-2003-12`,
                      `zeu_values-2004-1`, `zeu_values-2004-2`, `zeu_values-2004-3`, `zeu_values-2004-4`, `zeu_values-2004-5`, `zeu_values-2004-6`, 
                      `zeu_values-2004-7`, `zeu_values-2004-8`, `zeu_values-2004-9`, `zeu_values-2004-10`, `zeu_values-2004-11`, `zeu_values-2004-12`,
                      `zeu_values-2005-1`, `zeu_values-2005-2`, `zeu_values-2005-3`, `zeu_values-2005-4`, `zeu_values-2005-5`, `zeu_values-2005-6`, 
                      `zeu_values-2005-7`, `zeu_values-2005-8`, `zeu_values-2005-9`, `zeu_values-2005-10`, `zeu_values-2005-11`, `zeu_values-2005-12`,
                      `zeu_values-2006-1`, `zeu_values-2006-2`, `zeu_values-2006-3`, `zeu_values-2006-4`, `zeu_values-2006-5`, `zeu_values-2006-6`, 
                      `zeu_values-2006-7`, `zeu_values-2006-8`, `zeu_values-2006-9`, `zeu_values-2006-10`, `zeu_values-2006-11`, `zeu_values-2006-12`,
                      `zeu_values-2007-1`, `zeu_values-2007-2`, `zeu_values-2007-3`, `zeu_values-2007-4`, `zeu_values-2007-5`, `zeu_values-2007-6`, 
                      `zeu_values-2007-7`, `zeu_values-2007-8`, `zeu_values-2007-9`, `zeu_values-2007-10`, `zeu_values-2007-11`, `zeu_values-2007-12`,
                      `zeu_values-2008-1`, `zeu_values-2008-2`, `zeu_values-2008-3`, `zeu_values-2008-4`, `zeu_values-2008-5`, `zeu_values-2008-6`, 
                      `zeu_values-2008-7`, `zeu_values-2008-8`, `zeu_values-2008-9`, `zeu_values-2008-10`, `zeu_values-2008-11`, `zeu_values-2008-12`,
                      `zeu_values-2009-1`, `zeu_values-2009-2`, `zeu_values-2009-3`, `zeu_values-2009-4`, `zeu_values-2009-5`, `zeu_values-2009-6`, 
                      `zeu_values-2009-7`, `zeu_values-2009-8`, `zeu_values-2009-9`, `zeu_values-2009-10`, `zeu_values-2009-11`, `zeu_values-2009-12`,
                      `zeu_values-2010-1`, `zeu_values-2010-2`, `zeu_values-2010-3`, `zeu_values-2010-4`, `zeu_values-2010-5`, `zeu_values-2010-6`, 
                      `zeu_values-2010-7`, `zeu_values-2010-8`, `zeu_values-2010-9`, `zeu_values-2010-10`, `zeu_values-2010-11`, `zeu_values-2010-12`,
                      `zeu_values-2011-1`, `zeu_values-2011-2`, `zeu_values-2011-3`, `zeu_values-2011-4`, `zeu_values-2011-5`, `zeu_values-2011-6`, 
                      `zeu_values-2011-7`, `zeu_values-2011-8`, `zeu_values-2011-9`, `zeu_values-2011-10`, `zeu_values-2011-11`, `zeu_values-2011-12`,
                      `zeu_values-2012-1`, `zeu_values-2012-2`, `zeu_values-2012-3`, `zeu_values-2012-4`, `zeu_values-2012-5`, `zeu_values-2012-6`, 
                      `zeu_values-2012-7`, `zeu_values-2012-8`, `zeu_values-2012-9`, `zeu_values-2012-10`, `zeu_values-2012-11`, `zeu_values-2012-12`,
                      `zeu_values-2013-1`, `zeu_values-2013-2`, `zeu_values-2013-3`, `zeu_values-2013-4`, `zeu_values-2013-5`, `zeu_values-2013-6`, 
                      `zeu_values-2013-7`, `zeu_values-2013-8`, `zeu_values-2013-9`, `zeu_values-2013-10`, `zeu_values-2013-11`, `zeu_values-2013-12`,
                      `zeu_values-2014-1`, `zeu_values-2014-2`, `zeu_values-2014-3`, `zeu_values-2014-4`, `zeu_values-2014-5`, `zeu_values-2014-6`, 
                      `zeu_values-2014-7`, `zeu_values-2014-8`, `zeu_values-2014-9`, `zeu_values-2014-10`, `zeu_values-2014-11`, `zeu_values-2014-12`,
                      `zeu_values-2015-1`, `zeu_values-2015-2`, `zeu_values-2015-3`, `zeu_values-2015-4`, `zeu_values-2015-5`, `zeu_values-2015-6`, 
                      `zeu_values-2015-7`, `zeu_values-2015-8`, `zeu_values-2015-9`, `zeu_values-2015-10`, `zeu_values-2015-11`, `zeu_values-2015-12`,
                      `zeu_values-2016-1`, `zeu_values-2016-2`, `zeu_values-2016-3`, `zeu_values-2016-4`, `zeu_values-2016-5`, `zeu_values-2016-6`, 
                      `zeu_values-2016-7`, `zeu_values-2016-8`, `zeu_values-2016-9`, `zeu_values-2016-10`, `zeu_values-2016-11`, `zeu_values-2016-12`,
                      `zeu_values-2017-1`, `zeu_values-2017-2`, `zeu_values-2017-3`, `zeu_values-2017-4`, `zeu_values-2017-5`, `zeu_values-2017-6`, 
                      `zeu_values-2017-7`, `zeu_values-2017-8`, `zeu_values-2017-9`, `zeu_values-2017-10`, `zeu_values-2017-11`, `zeu_values-2017-12`,
                      `zeu_values-2018-1`, `zeu_values-2018-2`, `zeu_values-2018-3`, `zeu_values-2018-4`, `zeu_values-2018-5`, `zeu_values-2018-6`, 
                      `zeu_values-2018-7`, `zeu_values-2018-8`, `zeu_values-2018-9`, `zeu_values-2018-10`, `zeu_values-2018-11`, `zeu_values-2018-12`,
                      `zeu_values-2019-1`, `zeu_values-2019-2`, `zeu_values-2019-3`, `zeu_values-2019-4`, `zeu_values-2019-5`, `zeu_values-2019-6`, 
                      `zeu_values-2019-7`, `zeu_values-2019-8`, `zeu_values-2019-9`, `zeu_values-2019-10`, `zeu_values-2019-11`, `zeu_values-2019-12`,
                      `zeu_values-2020-1`, `zeu_values-2020-2`, `zeu_values-2020-3`, `zeu_values-2020-4`, `zeu_values-2020-5`, `zeu_values-2020-6`, 
                      `zeu_values-2020-7`, `zeu_values-2020-8`, `zeu_values-2020-9`, `zeu_values-2020-10`, `zeu_values-2020-11`, `zeu_values-2020-12`),
                    dim = c(length(range_lon),length(range_lat),length(range_time)))

##create output NetCDF file with monthly aggregates
lat_list <- list(ds[["dim"]][["latitude"]][["vals"]])
lon_list <- list(ds[["dim"]][["longitude"]][["vals"]])
time_list <- list(ds[["dim"]][["time"]][["vals"]])

# Define the dimensions
dim_lon <- ncdim_def("longitude", "degrees_east", range_lon)
dim_lat <- ncdim_def("latitude", "degrees_north", range_lat)
dim_time <- ncdim_def("time", "hours since 1950-01-01 00:00:00", range_time)

# Define the dimensions of the variables Zeu and Zooc
count <- c(lon_range, lat_range, length(range_time)) #dimension Zeu & Zooc

dimzoo <- ncdim_def("zooc", "g m-2", count)
dimzeu <- ncdim_def("zeu", "m", count)

var_zoo <- ncvar_def("zooc", "g m-2", list(dim_lon, dim_lat, dim_time), -32767, longname = "mass_content_of_zooplankton_expressed_as_carbon_in_sea_water", prec = "float")
var_zeu <- ncvar_def("zeu", "m", list(dim_lon, dim_lat, dim_time), -32767, longname = "euphotic_zone_depth", prec = "float")

# Create a new NetCDF file object
ncnew_zoo <- nc_create("zooplankton_concentration.nc", var_zoo)
ncvar_put(ncnew_zoo, var_zoo, zooc_array, start = NA, count = NA)
nc_close(ncnew_zoo)


ncnew_zeu <- nc_create("euphotic_zone_depth.nc", var_zeu)

# Write the chlorophyll data
ncvar_put(ncnew_zeu, var_zeu, zeu_array, start = NA, count = NA)

# Close the netCDF to save it
nc_close(ncnew_zeu)
