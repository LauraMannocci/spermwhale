
# load all functions in R
devtools::load_all()

####################  PREPARE ENVIRONMENT DATA

# read rasters of bathymetrically derived predictors - each has an ecological rational
depth.m <- read_bathy_data("depth")
slope.per <- read_bathy_data("slope")
di_1000.km <- read_bathy_data("di_1000m")
di_seaM.km <- read_bathy_data("di_seaM")
di_spRid.km <- read_bathy_data("di_spRid")

# stack all predictor rasters into a single stack
modelStack <- raster::stack(c(depth.m, slope.per, di_1000.km, di_seaM.km, di_spRid.km), quick=TRUE)

# reduce resolution by a factor of 10 to decrease processing time
modelStack10 <- raster::aggregate(modelStack, fact = 10)

# check correlations
cor <- raster::layerStats(modelStack10, 'pearson', na.rm = TRUE)

# stack again after removing collinear predictors (none here)
modelStack <- raster::stack(c(depth.m, slope.per, di_1000.km, di_seaM.km, di_spRid.km), quick=TRUE)

# reduce resolution by a factor of 20 to decrease processing time
modelStack <- raster::aggregate(modelStack, fact = 20)



# get multidimensional extrapolation extent - the next two lines took a few minutes to run on 50 threads
load(here::here("make_for_server.RData"))

df_extraMod <- get_extra_extent(c("depth",  "slope",    "di_1000m",  "di_seaM", "di_spRid"), modelStack, BgMod.z, "modern")
df_extraHis <- get_extra_extent(c("depth",  "slope",    "di_1000m",  "di_seaM", "di_spRid"), modelStack, BgHis.z, "historical")

save(df_extraMod, df_extraHis, file = here::here("make_for_pc.RData"))
