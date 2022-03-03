#' Read and clean historical data
#'
#' @return
#' @export
#'

read_clean_hist_data <- function(){
  
  dat <- read.csv(here::here("data", "presences", "FinalCetaceanBIOTHistPresence.csv"),stringsAsFactors = TRUE, header = T)
  dat <- dat[c("Lon", "Lat")]
  dat <- dat[!duplicated(dat),]
  dat <- dat[ which(dat$Lon < 84.75), ]
  
  return(dat)
}



#' Read and clean modern data
#'
#' @return
#' @export
#'

read_clean_mod_data <- function(){
  
  dat <- read.csv(here::here("data", "presences", "FinalCetaceanBIOTModPresence.csv"),stringsAsFactors = TRUE, header = T)
  dat <- dat[c("Lon", "Lat")]
  dat <- dat[!duplicated(dat),]
  dat <- dat[ which(dat$Lon < 84.75), ]
  
  return(dat)
}



#' Read rasters of bathymetrically derived predictors
#'
#' @param var 
#'
#' @return
#' @export
#'

read_bathy_data <- function(var){
  
  var <- raster::raster(here::here("data", "environment", paste0(var, ".tif")))
  
  return(var)
}


