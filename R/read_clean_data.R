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



#' Project survey transects to projection
#'
#' @param trans
#' @param proj
#'
#' @return
#' @export
#'

project_transect <- function(trans, proj){

  trans_proj <- sp::spTransform(trans, proj)

  return(trans_proj)
}





#' Read logbook historical data
#'
#' @return
#' @export
#'

read_logbook_data <- function(){

  dat <- readxl::read_xlsx(here::here("data", "whaling_logbooks", "AmericanOffshoreWhalingLogbookData", "aowl_download_20180104.xlsx"), col_names = T, skip = 1) %>%
    dplyr::mutate(Encounter = as.factor(Encounter)) %>%
    dplyr::mutate(Species = as.factor(Species)) %>%
    dplyr::mutate(Source = as.factor(Source)) -> dat2

  return(dat2)

}




#' Select logbook historical data
#'
#' @param dat
#' @param ext
#'
#' @return
#' @export
#'

select_logbook_data <- function(dat, ext){

  dat %>%
    #select non encounters
    dplyr::filter(Encounter == "NoEnc") %>%
    #add origin column
    dplyr::mutate(origin = "no_encounter") -> dat_noenc

  dat %>%
    #select non sperm whale presences
    dplyr::filter(Species != "Sperm" | is.na(Species)) %>%
    #add origin column
    dplyr::mutate(origin = "no_sperm_whale_sighting") -> dat_nonsw

  #extract lat/lon from study region extent
  lonmin <- ext[1]
  lonmax <- ext[2]
  latmin <- ext[3]
  latmax <- ext[4]

  #bind
  rbind(dat_noenc, dat_nonsw) %>%
    #restrict to study region
    dplyr::filter(Lat > latmin & Lat < latmax) %>%
    dplyr::filter(Lon > lonmin & Lon < lonmax) %>%
    #selection
    dplyr::select(Lon, Lat, origin) -> datnew

return(datnew)

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




#' Get a given number of random points from logbooks
#'
#' @param dat
#' @param n
#'
#' @return
#' @export
#'

get_random_points_from_logbooks <- function(dat, n){

  r <- dplyr::slice_sample(dat, n= n, replace = FALSE)

  r <- as.data.frame(r)

  return(r)

}




#' First, we need to define a function, which accepts two-column matrix of coordinates and lenghts from and to, which clip the segment.
#' Distances are always calculated from 0, therefore resulting segment will start at distance start from the beginning of the line and end at distance
#' end from the beginning.
#' from http://rstudio-pubs-static.s3.amazonaws.com/10685_1f7266d60db7432486517a111c76ac8b.html
#'
#' @param coords
#' @param from
#' @param to

CreateSegment <- function(coords, from, to) {

  distance <- 0
  coordsOut <- c()
  biggerThanFrom <- F
  for (i in 1:(nrow(coords) - 1)) {
    d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i +
                                                                             1, 2])^2)
    distance <- distance + d
    if (!biggerThanFrom && (distance > from)) {
      w <- 1 - (distance - from)/d
      x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
      y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
      coordsOut <- rbind(coordsOut, c(x, y))
      biggerThanFrom <- T
    }
    if (biggerThanFrom) {
      if (distance > to) {
        w <- 1 - (distance - to)/d
        x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
        y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
        coordsOut <- rbind(coordsOut, c(x, y))
        break
      }
      coordsOut <- rbind(coordsOut, c(coords[i + 1, 1], coords[i + 1,
                                                               2]))
    }
  }
  return(coordsOut)

}


#' Now, lets say we want to cut the given line to several segments with given length. Another option is to give number of desired parts (instead of their length)
#' as an argument.
#' from http://rstudio-pubs-static.s3.amazonaws.com/10685_1f7266d60db7432486517a111c76ac8b.html
#' @param coords
#' @param length
#' @param n.parts
#'
CreateSegments <- function(coords, length = 0, n.parts = 0) {

  stopifnot((length > 0 || n.parts > 0))
  # calculate total length line
  total_length <- 0
  for (i in 1:(nrow(coords) - 1)) {
    d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i +
                                                                             1, 2])^2)
    total_length <- total_length + d
  }

  # calculate stationing of segments
  if (length > 0) {
    stationing <- c(seq(from = 0, to = total_length, by = length), total_length)
  } else {
    stationing <- c(seq(from = 0, to = total_length, length.out = n.parts),
                    total_length)
  }

  # calculate segments and store the in list
  newlines <- list()
  for (i in 1:(length(stationing) - 1)) {
    newlines[[i]] <- CreateSegment(coords, stationing[i], stationing[i +
                                                                       1])
  }
  return(newlines)

}






#' Since the actual length of line is rarely a multiple of given length, last segment is shorter. We can, however, very simply merge it with penultimate segment,
#' if we wish.
#' from http://rstudio-pubs-static.s3.amazonaws.com/10685_1f7266d60db7432486517a111c76ac8b.html
#'
#' @param lst

MergeLast <- function(lst) {

  l <- length(lst)
  lst[[l - 1]] <- rbind(lst[[l - 1]], lst[[l]])
  lst <- lst[1:(l - 1)]
  return(lst)

}









#' Now, lets create SpatialLines and make function to perform segmentation on SpatialLines. Attributes of lines, even in case they are given as
#' SpatialLinesDataFrame, are not kept. If desired length of segments is bigger than actual length of original line, the original line is returned.
#' from http://rstudio-pubs-static.s3.amazonaws.com/10685_1f7266d60db7432486517a111c76ac8b.html
#'
#' @param sl
#' @param length
#' @param n.parts
#' @param merge.last

SegmentSpatialLines <- function(sl, length = 0, n.parts = 0, merge.last = FALSE) {

  stopifnot((length > 0 || n.parts > 0))
  id <- 0
  newlines <- list()
  sl <- as(sl, "SpatialLines")
  for (lines in sl@lines) {
    for (line in lines@Lines) {
      crds <- line@coords
      # create segments
      segments <- CreateSegments(coords = crds, length, n.parts)
      if (merge.last && length(segments) > 1) {
        # in case there is only one segment, merging would result into error
        segments <- MergeLast(segments)
      }
      # transform segments to lineslist for SpatialLines object
      for (segment in segments) {
        newlines <- c(newlines, Lines(list(Line(unlist(segment))), ID = as.character(id)))
        id <- id + 1
      }
    }
  }

  return(SpatialLines(newlines))

}





#' get observed absence points from transect segments
#'
#' @param seg
#' @param occ
#'
#' @return
#' @export
#'

get_observed_absence_points <- function(seg, occ){

  #transform segments to sf
  seg_sf <- sf::st_as_sf(seg)

  #transform occurrence to sf
  occ_sf <- sf::st_as_sf(x = occ,
                           coords = c("Lon", "Lat"),
                           crs = "+proj=longlat +datum=WGS84 +no_defs")

  #add id to segments
  seg_sf %>%
    dplyr::mutate(id = seq(1, nrow(seg_sf))) -> seg_sf

  #create a small buffer to allow intersection
  occ_sf_buf <- sf::st_buffer(occ_sf, 1000)

  #intersect segments and occurrences
  int <- sf::st_intersection(seg_sf, occ_sf_buf)

  #eliminate segments with intersection
  seg_sf %>%
    dplyr::filter(!id %in% int$id) %>%
    dplyr::select(-id) -> seg_sf_no_int

  #get centroids, ie observed absence points
  seg_sf_no_int_cen <- sf::st_centroid(seg_sf_no_int)

  return(seg_sf_no_int_cen)

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



#' Create background points from survey segments centroids (on the basis of bias file)
#'
#' @param modelstack
#' @param seg
#'
#' @return
#' @export
#'

create_background_points_from_survey <- function(modelstack, segs, nb){

  #transform segs to sf
  segs <- sf::st_as_sf(segs)

  #get centroids
  segs <- sf::st_centroid(segs)

  # create a raster file on the basis of survey points density
  biasfile <- modelstack[["depth"]]
  biasfile[] <- 0
  tab <- table(raster::cellFromXY(biasfile, sf::st_coordinates(segs)))
  biasfile[as.numeric(names(tab))] <- tab

  #select background points on the basis of the bias file
  bg <- raster::xyFromCell(biasfile, sample(raster::ncell(biasfile),
                                            prob = raster::values(biasfile),
                                            replace = FALSE,
                                            size = nb))
  colnames(bg) <- c("Lon", "Lat")

  bg <- as.data.frame(bg)

  return(bg)

}
