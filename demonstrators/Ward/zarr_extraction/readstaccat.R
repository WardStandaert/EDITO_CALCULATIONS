#reads the stac catalogue
library(foreach)
library(doParallel)

source("editoTools.R")
options("outputdebug"=c('M'))

cf <- S3FileSystem$create(
  access_key = '915692b659c44f77ad21330649ebfe8e',
  secret_key = 'c3b758f0fb09469d8f6f4911440a2738',
  scheme = "https",
  endpoint_override = "s3.waw3-1.cloudferro.com"
)

#stac files behind the url are hierarchical, recursive loop through the links called 'child' or 'item'
getChildrenLinks<-function(url="",allowed=c('child','item'),step=1 )
{
  
  if(is.na(url)) return()
  if(url=="") return()
  
  
  
  result=""
  
  try({ result= retryJson(url) })
  
  if( class(result)=="character") {dbl('ends here');return()}
  
  links=result$links
  
  if(is.null(links) ||  class(links)!='data.frame') {dbl('end here');return()}
  
  children= links  %>% dplyr::filter(rel %in% allowed)
  dbl('step ',step, 'links', nrow(children), url )
  
  if(nrow(children)==0) return(url)
  
  refs = stringr::str_replace(url,tail(unlist(stringr::str_split(url,stringr::fixed('/'))),n=1),stringr::str_remove(children$href,stringr::fixed("./")))
  
  
  all=c()
  i=1
  all <- foreach::foreach(i=1:length(refs) , .combine='c') %dopar% 
    {    getChildrenLinks(url=refs[i], step=step+1) }
  
  
  
  
  
  
  return(all)
  
}

#collect all info needed to decide on best zarr file to use from the stac file
#this function is very specific for the way the edito data lake is organized
#at extract time some additional info is collected from the zarr meta data
#the cmems and edito catalogues are slightly different, the url decides where to look for info in the file 
parseStacJson<-function(url)
{
  
  cat=retryJson(url)
  if(length(names(cat$assets)) == 0) return()
  
  
  thisstacdf=dplyr::tibble("asset"="","par"="","standardname"="","href"="","reftype"="","ori"="","latmin"=0,"latmax"=0,'lonmin'=0,'lonmax'=0,'latstep'=0,'lonstep'=0,'timemin'=0,'timemax'=0,'timestep'=0, 'timeunits'= "" 
                           , 'levels'=0, 'timetype'="", 'start_datetime'='','end_datetime'="",'title'="", "categories"=0 ,"catalogue"="", "chunktype"="")
  
  
  for( a in names(cat$assets)) 
  {  
    step3=cat$assets[[a]]
    vars=names(step3$viewVariables)
    
    
    dbl(cat$assets[[a]]$type)
    
    
    if(stringr::str_detect(url,"emodnet")) # &  str_detect(cat$assets[[a]]$type,"zarr")) # of ?  %in% c("application/zarr" , "application/vnd+zarr"))
    { 
      
      
      
      # try({ 
      thisstacdf=dplyr::add_row(thisstacdf,
                                "asset"=a,
                                "par"=coalesce(cat$properties$standardName,cat$properties$dataset_variable_name, cat$properties$native_dataset_parameter, cat$collection, 'missing' ), # needs translation after standardizing
                                "standardname"=coalesce(cat$properties$standardName,cat$properties$dataset_variable_name, cat$properties$native_dataset_parameter, cat$collection, 'missing' ),
                                "href"=cat$assets[[a]]$href,
                                "reftype"=cat$assets[[a]]$type,
                                "title"=ifelse(is.null(cat$properties$title), cat$id , cat$properties$title ) ,
                                "ori"=url,
                                "catalogue"="EMODNET",
                                "chunktype"=ifelse(str_detect(cat$assets[[a]]$type,"zarr"),"chunked","native"),
                                "latmin"=cat$properties$`cube:dimensions`$latitude$extent[1],
                                "latmax"=cat$properties$`cube:dimensions`$latitude$extent[2],
                                
                                
                                "lonmin"=cat$properties$`cube:dimensions`$longitude$extent[1],
                                "lonmax"=cat$properties$`cube:dimensions`$longitude$extent[2],
                                
                                
                                "start_datetime"=cat$properties$start_datetime,
                                "end_datetime"=cat$properties$end_datetime ,
                                
                                "categories"= ifelse(length(cat$properties$categorical_encodings)>0,length(cat$properties$categorical_encodings), 
                                                     ifelse(str_detect(cat$assets[[a]]$href,"scale") | str_detect(cat$assets[[a]]$href,"eunis") , 1, 0) ) 
      )
      #or try this
      if(is.na(thisstacdf$lonmax[2] )) {
        thisstacdf$lonmax[2]=cat$bbox[3]
        thisstacdf$latmax[2]=cat$bbox[4]
        thisstacdf$lonmin[2]=cat$bbox[1]
        thisstacdf$latmin[2]=cat$bbox[2]
      }
      # })  
    }  
    else if(stringr::str_detect(url,"cmems") &   !is.null(vars))
    {
      for(g in 1:length(vars)){
        try({ thisstacdf=dplyr::add_row(thisstacdf,
                                        "asset"=a,
                                        "par"=vars[g],
                                        "standardname"=cat$properties$`cube:variables`[[vars[g]]]$standardName,
                                        "href"=cat$assets[[a]]$href,
                                        "title"=ifelse(is.null(cat$properties$title), cat$id , cat$properties$title ) ,
                                        "ori"=url,
                                        "catalogue"="CMEMS",
                                        "chunktype"=a,
                                        "latmin"=cat$assets[[a]]$viewDims$latitude$coords$min,
                                        "latmax"=cat$assets[[a]]$viewDims$latitude$coords$max,
                                        "latstep"=cat$assets[[a]]$viewDims$latitude$coords$step,
                                        
                                        "lonmin"=cat$assets[[a]]$viewDims$longitude$coords$min,
                                        "lonmax"=cat$assets[[a]]$viewDims$longitude$coords$max,
                                        "lonstep"=cat$assets[[a]]$viewDims$longitude$coords$step,
                                        
                                        "timemin"=cat$assets[[a]]$viewDims$time$coords$min,
                                        "timemax"=cat$assets[[a]]$viewDims$time$coords$max,
                                        "timestep"=cat$assets[[a]]$viewDims$time$coords$step,
                                        "timetype"=cat$assets[[a]]$viewDims$time$coords$type,
                                        
                                        "timeunits"=cat$assets[[a]]$viewDims$time$units,
                                        
                                        "start_datetime"=cat$properties$start_datetime,
                                        "end_datetime"=cat$properties$end_datetime,
                                        
                                        "levels"=cat$assets[[a]]$viewDims$elevation$len
                                        
        )})
        
        dbl(a,vars[g])
      }
    }
    
    thisstacdf=dplyr::filter(thisstacdf, par!='')
    
  }
  
  
  return(thisstacdf)
}

loadStacCatalogue<-function(catalogue=c('https://s3.waw3-1.cloudferro.com/emodnet/stac/catalog.json',"https://s3.waw3-1.cloudferro.com/mdl-metadata/metadata/catalog.stac.json"))
  # or "https://s3.waw3-1.cloudferro.com/mdl-metadata/metadata/catalog.stac.json"
  # or is it https://s3.waw3-1.cloudferro.com/emodnet/stac/catalog.json
  # "https://edito-infra.dev.emodnet.eu/emodnetstacdev/catalog.json"
{
  doParallel::registerDoParallel(16)
  dbl("workers:",getDoParWorkers())
  
  dbl("getting links")
  links=c()
  
  for(c in catalogue)
    links=c(links, getChildrenLinks(c))
  
  #or test like this :  links=c('https://s3.waw3-1.cloudferro.com/emodnet/stac/elevation/climate_forecast-geoid_height_above_reference_ellipsoid/geoid_height_above_reference_ellipsoid_1816-01-01T00:00:00_2022-07-11T00:00:00_-36.0_25.0_43.0_85.0/geoid_height_above_reference_ellipsoid_1816-01-01T00:00:00_2022-07-11T00:00:00_-36.0_25.0_43.0_85.0.json')
  
  dbl(length(links), " assets")
  options("outputdebug"=c('M'))
  
  
  
  if(length(links)==0) stop("no links in url")
  
  Estacdf=tibble()
  Estacdf <- foreach(l=1:length(links) ,.options.silent=F, .combine='CombineR') %dopar%     
    { 
      
      dbl( l, links[l] )
      tble=tibble()
      try({
        tble=parseStacJson(links[l])
      }) 
      dbl( nrow(tble) )
      
      return(tble)
    }
  dbl("got " , nrow(Estacdf) , " datasets ")
  
  doParallel::stopImplicitCluster()
  
  return(Estacdf) # save later to EDITOSTAC.par
}

#kan beter.., parallel enz

ns = loadStacCatalogue()

arrow::write_dataset(ns, cf$path("emodnet/edito_stac_cache.parquet"), format="parquet")





arraysInZarr<-function(fn){
  r=jsonlite::fromJSON(gdal_utils("mdiminfo", toGDAL(fn), quiet = T))
  if(!is.null(r$arrays)) return(length(r$arrays))
  
  return(0)
}  

