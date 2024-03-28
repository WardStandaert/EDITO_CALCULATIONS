#reads the stac catalogue
library(foreach)
library(doParallel)

source("helpers.R")
options("outputdebug"=c('M'))

getChildrenLinks<-function(url="https://edito-infra.dev.emodnet.eu/emodnetstacdev/catalog.json",allowed=c('child','item'),step=1 )
{
  
 
  
  dbl('step',step)
  links= retryJson(url)$links
  if(is.null(links) ||  class(links)!='data.frame') {dbl('end here');return(paste('dont',url))}
  
  children= links  %>% dplyr::filter(rel %in% allowed)
  refs = stringr::str_replace(url,tail(unlist(stringr::str_split(url,stringr::fixed('/'))),n=1),stringr::str_remove(children$href,stringr::fixed("./")))
  #dbl('to do ',length(refs))
  
  if(length(refs)==0) return(url)
  
  
all=c()
# 
#    for(r in refs) {
#      newlinks=getChildrenLinks(url=r, step=step+1)
#   #   dbl('new links', length(newlinks))
#      if(length(newlinks)>0) all=c(all,newlinks)
#     
#    }
  
  all= foreach(r=1:length(refs) ,.options.silent=F, .combine='append') %dopar%    
    {
      newlinks=getChildrenLinks(url=refs[r], step=step+1)
         
      
    }
   return(all)
  
}

parseStacJson<-function(url)
{
  
  cat=retryJson(url)
  if(length(names(cat$assets)) == 0) return()
  
  thisstacdf=dplyr::tibble("asset"="","par"="","standardname"="","href"="","reftype"="","ori"="","latmin"=0,"latmax"=0,'lonmin'=0,'lonmax'=0,'latstep'=0,'lonstep'=0,'timemin'=0,'timemax'=0,'timestep'=0, 'timeunits'= "" 
                       , 'levels'=0, 'timetype'="", 'start_datetime'='','end_datetime'="",'title'="", "categories"=0 ,"catalogue"="", "chunktype"="")
  #the case with cmems zarr files
  for( a in names(cat$assets)) 
  { 
    step3=cat$assets[[a]]
    vars=names(step3$viewVariables)
    
    if(stringr::str_detect(url,"emodnet") &  cat$assets[[a]]$type== "application/zarr")
    {
      try({ 
        thisstacdf=dplyr::add_row(thisstacdf,
                                  "asset"=a,
                                  "par"=cat$properties$dataset_variable_name, # needs translation after standardizing
                                  "standardname"=cat$properties$dataset_variable_name,
                                  "href"=cat$assets[[a]]$href,
                                  "reftype"=cat$assets[[a]]$type,
                                  "title"=ifelse(is.null(cat$properties$title), cat$id , cat$properties$title ) ,
                                  "ori"=url,
                                  "catalogue"="EMODNET",
                                  "chunktype"="chunked",
                                  "latmin"=cat$properties$`cube:dimensions`$latitude$extent[1],
                                  "latmax"=cat$properties$`cube:dimensions`$latitude$extent[2],
                                  
                                  
                                  "lonmin"=cat$properties$`cube:dimensions`$longitude$extent[1],
                                  "lonmax"=cat$properties$`cube:dimensions`$longitude$extent[2],
                                  
                                  
                                  "start_datetime"=cat$properties$start_datetime,
                                  "end_datetime"=cat$properties$end_datetime ,
                                  
                                  "categories"= length(cat$properties$categorical_encodings)
        )
        #or try this
        if(is.na(thisstacdf$lonmax[2] )) {
          thisstacdf$lonmax[2]=cat$bbox[3]
          thisstacdf$latmax[2]=cat$bbox[4]
          thisstacdf$lonmin[2]=cat$bbox[1]
          thisstacdf$latmin[2]=cat$bbox[2]
        }
      })  
    }  
    else if(!is.null(vars))
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
   thisstacdf=dplyr::filter(thisstacdf, asset != "")
    
  }
 return(thisstacdf)
}

loadEMODNETStacCatalogue<-function()
{
  registerDoParallel(16)
  dbl("workers:",getDoParWorkers())
  
  dbl("getting links")
  options("outputdebug"=c('silent'))
  Estacdf=tibble()
  elinks=getChildrenLinks("https://edito-infra.dev.emodnet.eu/emodnetstacdev/catalog.json")
  links=getChildrenLinks("https://s3.waw3-1.cloudferro.com/mdl-metadata/metadata/catalog.stac.json")
  links=c(elinks,links)
  options("outputdebug"=c('M'))
  
  dbl(length(links), "links")
  
  if(length(links)==0) stop("no links in url")

Estacdf= foreach(l=1:length(links) ,.options.silent=F, .combine='CombineR') %dopar%     
  {
    tble=parseStacJson(links[l])
    dbl(nrow(tble) )
    return(tble)
  }
  dbl("got " , nrow(Estacdf) , " datasets ")
  return(Estacdf) # save later to EDITOSTAC.par
}

#kan beter.., parallel enz
testAssets<-function(stacdf)
{
for(i in 1:nrow(stacdf)){
  r=gdal_utils("mdiminfo", toGDAL(stacdf$href[i]), quiet = T)
  if(!stacdf$par[i] %in% names(jsonlite::fromJSON(r)$arrays)  ) print(paste0(stacdf$href[i],  " not found")) 
  cat('.')
  }  
}

