# load functions, but also cached stac catalog named stacCatalog
source("editoTools.R")
options("outputdebug"=c('L','M'))

#the cached stacCatalog is called 'EDITOSTAC'
load(file = "./data-raw/editostacv2.par")



# the file to process
datafile="./data-raw/extract_test.csv"

pts=read.delim(datafile, sep=",") 



#the requested timestep resolution of the dataset in milliseconds
#in this case we work with daily data (1day = 24*3600*1000 = 8640000 milliseconds)
timeSteps=c(86400000)
timeSteps=c()


#the requested parameters, names in the stac catalogue table field 'par' .. see unique(stacCatalogue$par) for a list 
parameters= c('Substrate')


#parameters= c('Substrate','elevation','thetao','so','phyc','zooc','Energy')


#check if they all exist
for ( parameter in parameters) {
  if(! parameter %in% unique(EDITOSTAC$par) ) 
    { dbl("Unknown parameter ", parameter) 
      
  }
}


#add verbose= anything to get additional info on the positions ( par_x, par_y, par_z ) and time (par_t) found in the zarr files
enhanced_DF=enhanceDF(inputPoints = pts , requestedParameters = parameters, requestedTimeSteps = timeSteps, stacCatalogue = EDITOSTAC, verbose="on")


#check if position and time 
#View(enhanced_DF %>% select( Latitude,Longitude, any_of(c('Date','Time')) , ends_with(c('_y','_x','_t'))))
View(enhanced_DF %>% select(  -ends_with(c('_y','_x','_t'))))



