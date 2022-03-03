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
  dat <- dat[which(dat$Lon < 84.75), ]
  #remove erroneous observation
  dat <- dat[which(dat$Lat > -33.34), ]

  return(dat)
}





#' Read transects data
#'
#' @param var
#'
#' @return
#' @export
#'

read_transect_data <- function(var){

  if (var == "OA"){
      spf <- "OA_line.shp"
  }
  if (var == "REMMOA"){
    spf <- "Transects_realises.shp"
  }
  if (var == "NOAA"){
    spf <- "NOAAline.shp"
  }

  trans = rgdal::readOGR(dsn = here::here("data", "effortlines", spf))

  return(trans)
}



#' Read rasters of bathymetrically derived predictors
#'
#' @param var
#'
#' @return
#' @export
#'

read_bathy_data <- function(var){

dat <- raster::raster(here::here("data", "environment", paste0(var, ".tif")))

return(dat)
}


#' Remove duplicate occurrences (ie those that share a grid cell in the predictor variable raster)
#'
#' @param modelstack
#' @param occ
#'
#' @return
#' @export
#'

remove_duplicates_occ <- function(modelstack, occ){

  occscells <- raster::extract(modelstack[[1]], occ, cellnumbers = TRUE)
  occscelldups <- duplicated(occscells[,1])
  occ <- occ[!occscelldups,]

  return(occ)
}



#' Extract bathymetric data at point (presence or pseudo-absence point)
#'
#' @param modelstack
#' @param type
#'
#' @return
#' @export
#'

extract_bathy_at_point <- function(modelstack, type){

dat <- na.omit(cbind(type, raster::extract(modelstack, type)))

return(dat)
}





#' Get a given number of random points
#'
#' @param modelstack
#' @param ext
#' @param occ
#' @param nb
#'
#' @return
#' @export
#'

get_random_points <- function(modelstack, ext, occ, nb){

  occ <- dismo::randomPoints(mask = modelstack, ext = ext, p = occ, nb)
  colnames(occ) <- c("Lon", "Lat")
  return(occ)

}



#' Define equally spaced points along survey transects
#'
#' @param trans
#' @param length
#'
#' @return
#' @export
#'

define_points_along_survey_transects <- function(trans, length){

  trans_sf <- sf::st_as_sf(trans)
  sum <- units::drop_units(sum(sf::st_length(trans_sf))) / 1000
  nbpts <- sum / length
  trans_pts <- sf::st_sample(trans_sf, size = round(nbpts) , type="regular")
  trans_pts = trans_pts[!sf::st_is_empty(trans_pts), drop=FALSE] #drop empty geometries
  trans_pts <- as(trans_pts, "Spatial")

  return(trans_pts)

}



#' Create background points from survey points (on the basis of bias file)
#'
#' @param modelstack
#' @param pts
#'
#' @return
#' @export
#'

create_background_points_from_survey <- function(modelstack, pts, nb){

  # create a raster file on the basis of survey points density
  biasfile <- modelstack[["depth"]]
  biasfile[] <- 0
  tab <- table(raster::cellFromXY(biasfile, pts))
  biasfile[as.numeric(names(tab))] <- tab

  #select background points on the basis of the bias file, keeping the ratio of occurrence points to background points available in historical data
  ratio_occ_to_bg <- (nrow(SPMod) * nb)/nrow(SPHis)
  bg <- raster::xyFromCell(biasfile, sample(raster::ncell(biasfile),
                                            ratio_occ_to_bg,
                                            prob = raster::values(biasfile),
                                            size = ratio_occ_to_bg))
  colnames(bg) <- c("Lon", "Lat")

  return(bg)

}
