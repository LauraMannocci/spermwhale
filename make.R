#initialize see for random sampling
set.seed(1)

# load all functions in R
devtools::load_all()



#save / load rdata
#load(here::here("make.RData"))
#save(df_predMod, df_extraMod, df_predHis, df_extraHis, BgMod.z, BgHis.z, predHis, predMod, SPMod.x, SPHis.x, opt.seqMod, opt.seqHis, mod.seqHis, mod.seqMod, file = here::here("make.RData"))






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
write.csv(cor, here::here("outputs", "raster_correlations.csv"))

# stack again after removing collinear predictors (none here)
modelStack <- raster::stack(c(depth.m, slope.per, di_1000.km, di_seaM.km, di_spRid.km), quick=TRUE)

# reduce resolution by a factor of 20 to decrease processing time
modelStack <- raster::aggregate(modelStack, fact = 20)






####################  READ HISTORICAL AND MODERN PRESENCE DATA

# read and clean historical data
SPHis <- read_clean_hist_data()

# read and clean modern data
SPMod <- read_clean_mod_data()

# Removing occurrences that have the same coordinates is good practice to avoid pseudoreplication.
SPHis <- SPHis[!duplicated(SPHis),]
SPMod <- SPMod[!duplicated(SPMod),]







####################  READ TRANSECTS MODERN DATA

oa <- read_transect_data("OA")
remmoa <- read_transect_data("REMMOA")
noaa <- read_transect_data("NOAA")

# project transects according to custom projection
custom_proj <- sp::CRS("+proj=aeqd +lat_0=0 +lon_0=60 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
oa_proj <- project_transect(oa, custom_proj)
remmoa_proj <- project_transect(remmoa, custom_proj)
noaa_proj <- project_transect(noaa, custom_proj)







####################  READ LOGBOOK HISTORICAL DATA

#read logbook data
log <- read_logbook_data()

#clean logbook data
log <- select_logbook_data(log, raster::extent(depth.m))










####################  MAKE QUICK MAPS

# map sperm whale occurrences
jpeg(here::here("outputs", "occurrences.jpeg"), width = 1100, height = 550)
par(mfrow=c(1,2), mar=c(2,2,2,6))
raster::plot(depth.m, main="Historic", legend = FALSE)
maps::map('world', fill=T, col= "grey", add=TRUE)
points(SPHis,col="black", pch=20)
raster::plot(depth.m, main="Modern")
maps::map('world',fill=T , col= "grey", add=TRUE)
points(SPMod,col="black", pch=20)
dev.off()


# map rasters
p1 <- rasterVis::levelplot(modelStack[[1]], col.regions = rev(terrain.colors(255)), cuts=254, margin=FALSE,
                           colorkey = list(title = "[m]",
                                           height=1, width=1.4,
                                           labels = list(cex = 1)),
                           xlab = "", ylab = "", main = names(modelStack[[1]]))

p2 <- rasterVis::levelplot(modelStack[[2]], col.regions = rev(terrain.colors(255)), cuts=254, margin=FALSE,
                           colorkey = list(title = "[%]",
                                           height=1, width=1.4,
                                           labels = list(cex = 1)),
                           xlab = "", ylab = "", main = names(modelStack[[2]]))

p3 <- rasterVis::levelplot(modelStack[[3]], col.regions = rev(terrain.colors(255)), cuts=254, margin=FALSE,
                           colorkey = list(title = "[km]",
                                           height=1, width=1.4,
                                           labels = list(cex = 1)),
                           xlab = "", ylab = "", main = names(modelStack[[3]]))

p4 <- rasterVis::levelplot(modelStack[[4]], col.regions = rev(terrain.colors(255)), cuts=254, margin=FALSE,
                           colorkey = list(title = "[km]",
                                           height=1, width=1.4,
                                           labels = list(cex = 1)),
                           xlab = "", ylab = "", main = names(modelStack[[4]]))

p5 <- rasterVis::levelplot(modelStack[[5]], col.regions = rev(terrain.colors(255)), cuts=254, margin=FALSE,
                           colorkey = list(title = "[km]",
                                           height=1, width=1.4,
                                           labels = list(cex = 1)),
                           xlab = "", ylab = "", main = names(modelStack[[5]]))

# multiplot
jpeg(here::here("outputs", "rasters.jpeg"), width = 1340, height = 960)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol=3)
dev.off()







####################  PREPARATION OF OCCURRENCES

# Let's now remove occurrences that are cell duplicates -- these are
# occurrences that share a grid cell in the predictor variable rasters.

SPMod <- remove_duplicates_occ(modelStack, SPMod)
SPHis <- remove_duplicates_occ(modelStack, SPHis)


# Extract the predictor variable values at our occurrence localities.
SPHis.z <- extract_bathy_at_point(modelStack, SPHis)
SPMod.z <- extract_bathy_at_point(modelStack, SPMod)


# Violin plots of predictors for occurrences
SPHis.z$type <- factor(c('Historical'))
SPMod.z$type <- factor(c('Modern'))
Mod.Hist <- rbind(SPHis.z, SPMod.z)

depthP <- predictor_violin_occurrences(Mod.Hist, "depth", "Seabed depth (m)")
slopeP <- predictor_violin_occurrences(Mod.Hist, "slope", "Seabed slope (%)")
di_spRidP <- predictor_violin_occurrences(Mod.Hist, "di_spRid", "Distance to spreading ridge (km)")
di_seaMP <- predictor_violin_occurrences(Mod.Hist, "di_seaM", "Distance to seamount (km)")
di_1000mP <- predictor_violin_occurrences(Mod.Hist, "di_1000m", "Distance to 1,000 m contour (km)")



# multiplot
jpeg(here::here("outputs", "predictors_violin_occurrences.jpeg"), width = 1550, height = 1200)
cowplot::ggdraw() +
  cowplot::draw_plot(depthP, 0.005, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(slopeP, 0.34, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(di_1000mP, 0.67, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(di_seaMP, 0.005, 0.1, 0.32, 0.4) +
  cowplot::draw_plot(di_spRidP, 0.34, 0.1, 0.32, 0.4) +
  cowplot::draw_text("(a)", x = 0.03, y = 0.92, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.35, y = 0.92, size = 25, fontface = "italic") +
  cowplot::draw_text("(c)", x = 0.69, y = 0.92, size = 25, fontface = "italic") +
  cowplot::draw_text("(d)", x = 0.03, y = 0.47, size = 25, fontface = "italic") +
  cowplot::draw_text("(e)", x = 0.35, y = 0.47, size = 25, fontface = "italic")
dev.off()








####################  SELECTION OF BACKGROUND POINTS (n=10000)


#-------- Historical data: get random background points from whaling logbooks
BgHis <- get_random_points_from_logbooks(log, 10000)


#-------- Modern data: generate background points from observed absences on transects

# segment transects into 10km pieces
library(sp) #to ensure SegmentSpatialLines is working
oa_proj_seg <- SegmentSpatialLines(sl = oa_proj, length = 10000, merge.last = TRUE)
remmoa_proj_seg <- SegmentSpatialLines(sl = remmoa_proj, length = 10000, merge.last = TRUE)
noaa_proj_seg <- SegmentSpatialLines(sl = noaa_proj, length = 10000, merge.last = TRUE)
proj4string(oa_proj_seg) <- custom_proj
proj4string(remmoa_proj_seg) <- custom_proj
proj4string(noaa_proj_seg) <- custom_proj


# reproject to geographic projection
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs")
oa_seg <- project_transect(oa_proj_seg, proj)
remmoa_seg <- project_transect(remmoa_proj_seg, proj)
noaa_seg <- project_transect(noaa_proj_seg, proj)


# merge all segments
sp.lines <- list()
sp.lines[[1]] <- oa_seg
sp.lines[[2]] <- remmoa_seg
sp.lines[[3]] <- noaa_seg
seg <- do.call(rbind, sp.lines)


# create a raster file on the basis of survey segment centroids density (on the basis of bias file)
#(use 1/12 resolution for bias file)
BgMod <- create_background_points_from_survey(modelStack10, seg, 10000)



# write background points
write.csv(BgMod, here::here("outputs", "background_points_modern.csv"), row.names = F)
write.csv(BgHis, here::here("outputs", "background_points_historial.csv"), row.names = F)



# once the randomPoints have been generated, we can batch-extract all variables from the raster stack
BgMod.z <- extract_bathy_at_point(modelStack, BgMod)
BgHis.z <- extract_bathy_at_point(modelStack, BgHis[,c("Lon", "Lat")])



# map background points
jpeg(here::here("outputs", "background_points.jpeg"), width = 1500, height = 750)
par(mfrow=c(1,2), mar=c(2,2,2,6))
raster::plot(depth.m, main = "Historical (n = 10000)")
maps::map('world', fill=T , col= "grey", add=TRUE)
points(BgHis, col="blue", pch=20, cex=0.1)
raster::plot(depth.m, main = "Modern (n = 10000)")
maps::map('world', fill=T , col= "grey", add=TRUE)
points(BgMod, col="blue", pch=20, cex=0.1)
dev.off()



## remove unnecessary objects

rm(abs, cor, depth.m, di_1000.km, di_guy.km, di_rid.km, di_seaM.km, di_spRid.km, di_tren.km, di_trough.km,
   noaa, noaa_proj, noaa_proj_seg, noaa_seg, oa, oa_proj, oa_proj_seg, oa_seg, remmoa, remmoa_proj, remmoa_proj_seg,
   remmoa_seg, seg, slope.per, sp.lines)









#################### MODEL EVALUATION


# do the checkerboard partitioning
cbMod <- checkerboard_partitioning(SPMod, SPMod.z, BgMod, BgMod.z, modelStack)
cbHis <- checkerboard_partitioning(SPHis, SPHis.z, BgHis, BgHis.z, modelStack)


# model evaluation with user partition
SPMod.x <- evaluate_model(SPMod,
                          modelStack,
                          BgMod,
                          list(fc = c("L", "Q"), rm = 1), #only linear fc / impose low penalty on model complexity
                          list(occs.grp = cbMod$occs.grp, bg.grp = cbMod$bg.grp),
                          "modern")

SPHis.x <- evaluate_model(SPHis,
                          modelStack,
                          BgHis[,c("Lon", "Lat")],
                          list(fc = c("L", "Q"), rm = 1), #only linear fc /  impose low penalty on model complexity
                          list(occs.grp = cbHis$occs.grp, bg.grp = cbHis$bg.grp),
                          "historical")







####################  MODEL SELECTION AND OPTIMISATION

# select best model
opt.seqMod <- select_best_model(SPMod.x, "Modern")
opt.seqHis <- select_best_model(SPHis.x, "Historical")


# write best model metrics
best_model_metrics <- rbind(opt.seqHis, opt.seqMod)
best_model_metrics$type <- factor(c('Historical', 'Modern'))
write.csv(best_model_metrics, here::here("outputs", "metrics_best_models.csv"))


# We can select a single model from the ENMevaluation object using the tune.args of our optimal model
mod.seqHis <- ENMeval::eval.models(SPHis.x)[[opt.seqHis$tune.args]]
mod.seqMod <- ENMeval::eval.models(SPMod.x)[[opt.seqMod$tune.args]]









####################  MODEL COEFFICIENTS AND PARTIAL PLOTS

# Here are the non-zero coefficients in our model.
coefHis <- tibble::enframe(mod.seqHis$betas)
coefMod <- tibble::enframe(mod.seqMod$betas)

coefHis$type <- factor(c('Historical'))
coefMod$type <- factor(c('Modern'))

Coef <- rbind(coefHis, coefMod)
write.csv(Coef, here::here("outputs", "coefficients.csv"))



# plot coefficients
coefsp <- make_coefficients_plot(Coef, c(-0.0015, 0.0015))



# barplot version
coefBar <- make_coefficients_barplot(Coef, c(-0.0015, 0.0015))





# partial plots

# And these are the marginal response curves for the predictor variables with non-zero
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).


# on separate plots
jpeg(here::here("outputs", "partial_plot_historical.jpeg"), width = 1000, height = 700)
plot_partial_curves(mod.seqHis, type = "cloglog", lwd = 5, col = "#F8766D")
dev.off()

jpeg(here::here("outputs", "partial_plot_modern.jpeg"), width = 1000, height = 700)
plot_partial_curves(mod.seqMod, type = "cloglog", lwd = 5, col = "#00BA38")
dev.off()



# together on the same plot

dat_his = plot_partial_curves(mod.seqHis, type = "cloglog")
dat_mod = plot_partial_curves(mod.seqMod, type = "cloglog")


p1 <- plot_partial_curves_together(dat_his, dat_mod, "depth")
p2 <- plot_partial_curves_together(dat_his, dat_mod, "slope")
p3 <- plot_partial_curves_together(dat_his, dat_mod, "di_1000m")
p4 <- plot_partial_curves_together(dat_his, dat_mod, "di_seaM")
p5 <- plot_partial_curves_together(dat_his, dat_mod, "di_spRid")

# multiplot
jpeg(here::here("outputs", "partial_plots.jpeg"), width = 1000, height = 700)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol=3)
dev.off()




# multiplot partial plots and barplot
jpeg(here::here("outputs", "partial_plots_barplot.jpeg"), width = 1400, height = 1400)
cowplot::ggdraw() +
  cowplot::draw_plot(coefBar, 0.08, 0.52, 0.9, 0.47) +
  cowplot::draw_plot(p1, 0.05, 0.25, 0.25, 0.25) +
  cowplot::draw_plot(p2, 0.38, 0.25, 0.25, 0.25) +
  cowplot::draw_plot(p3, 0.71, 0.25, 0.25, 0.25) +
  cowplot::draw_plot(p4, 0.05, 0, 0.25, 0.25) +
  cowplot::draw_plot(p5, 0.38, 0, 0.25, 0.25) +
  cowplot::draw_text("(a)", x = 0.08, y = 0.97, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.04, y = 0.48, size = 25, fontface = "italic") +
  cowplot::draw_text("(c)", x = 0.37, y = 0.48, size = 25, fontface = "italic") +
  cowplot::draw_text("(d)", x = 0.7, y = 0.48, size = 25, fontface = "italic") +
  cowplot::draw_text("(e)", x = 0.04, y = 0.23, size = 25, fontface = "italic") +
  cowplot::draw_text("(f)", x = 0.37, y = 0.23, size = 25, fontface = "italic")
dev.off()








####################  MODEL PREDICTIONS


# make predictions from selected model
predMod <- ENMeval::eval.predictions(SPMod.x)[[opt.seqMod$tune.args]]
predHis <- ENMeval::eval.predictions(SPHis.x)[[opt.seqHis$tune.args]]


# get countries and mpa for plotting
wio <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") #world map with 50m resolution
mpa_sf <- sf::read_sf(here::here("data", "mpa_sf2", "mpa_sf.shp"), stringsAsFactors=TRUE)


# convert prediction rasters to dataframes for plotting
df_predHis <- convert_predictions_to_df(predHis, "historical")
df_predMod <- convert_predictions_to_df(predMod, "modern")


# bind prediction dataframes
df_pred <- rbind(df_predHis, df_predMod)


# plot predictions
gMod <- plot_predictions(wio, df_predMod, "modern")
gHis <- plot_predictions(wio, df_predHis, "historical")


# plot predictions residuals
res <- plot_predictions_residuals(wio, predHis, predMod)









####################  VISUALISE EXTRAPOLATION EXTENT

# get multidimensional extrapolation extent - the next two lines of code take a few minutes to run on 50 threads on a server
# but they will crash a laptop computer
# so here we load df_extraMod and df_extraHis objects generated on a server
load(("make_for_pc.RData"))
# df_extraMod <- get_extra_extent(c("depth",  "slope",    "di_1000m",  "di_seaM", "di_spRid"), modelStack, BgMod.z, "modern")
# df_extraHis <- get_extra_extent(c("depth",  "slope",    "di_1000m",  "di_seaM", "di_spRid"), modelStack, BgHis.z, "historical")


# Plot predictions with extrapolation extent
gMod <- plot_predictions_with_extra(wio, df_predMod, "Modern", df_extraMod)
gHis <- plot_predictions_with_extra(wio, df_predHis, "Historical", df_extraHis)



# Plot predictions with extrapolation extent and mpas
gMod <- plot_predictions_with_extra_mpas(wio, df_predMod, "Modern", df_extraMod, mpa_sf)
gHis <- plot_predictions_with_extra_mpas(wio, df_predHis, "Historical", df_extraHis, mpa_sf)




# plot predictions residuals with extrapolation extent
res <- plot_predictions_residuals_with_extra(wio, predHis, predMod, df_extraMod, df_extraHis)




####################  PREDICTIONS COMPARISON AND RESIDUAL MPA EFFECT



# compare predictions (removing extrapolation zones)
compare_predictions_extra(predHis, df_extraHis, predMod, df_extraMod)



# test for "residual" mpa effect with kruskall test (removing extrapolation zones)
df_mpaHis <- test_residual_mpa_effect_extra(predHis,  df_extraHis, mpa_sf, "Historical")
df_mpaMod <- test_residual_mpa_effect_extra(predMod, df_extraMod, mpa_sf, "Modern")



# density plot of prediction distributions for all predictions vs predictions in mpa (removing extrapolation zones)
df_mpa <- rbind(df_mpaHis, df_mpaMod)
df_mpa$type <- factor(df_mpa$type , levels=c("region", "mpa"))

df_mpa <- df_mpa %>%
  dplyr::mutate(model = forcats::fct_relevel(model, levels = "Modern", "Historical"))

dens <- ggplot2::ggplot(data = df_mpa, ggplot2::aes(x = value, y = model, fill = model, alpha = type, linetype = type, scale = 0.9)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = c("#F8766D","#00BA38")) +
  ggplot2::scale_alpha_manual(values = c(0.6, .1)) +
  ggplot2::xlab('Habitat suitability') +
  ggplot2::ylab("") +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "none",
                 axis.text = ggplot2::element_text(size = 18),
                 axis.title = ggplot2::element_text(size = 18)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2)

ggplot2::ggsave(here::here("outputs", "density_plot_all_models.png"), dens, width = 9, height = 7)




# Calculate median predictions in and out of mpas (removing extrapolation zones)
calculate_median_predictions_in_out_mpa_extra(predHis, df_extraHis, mpa_sf)
calculate_median_predictions_in_out_mpa_extra(predMod, df_extraMod, mpa_sf)

# Calculate median predictions in mpas vs in all region (removing extrapolation zones)
calculate_median_predictions_in_mpa_all_region_extra(predHis, df_extraHis, mpa_sf)
calculate_median_predictions_in_mpa_all_region_extra(predMod, df_extraMod, mpa_sf)


# multiplot
jpeg(here::here("outputs", "summary_plot1.jpeg"), width = 1320, height = 960)
cowplot::ggdraw() +
  cowplot::draw_plot(gHis, 0.02, 0.52, 0.46, 0.46) +
  cowplot::draw_plot(dens, 0.46, 0.52, 0.46, 0.46) +
  cowplot::draw_plot(gMod, 0.02, 0.02, 0.46, 0.46) +
  cowplot::draw_plot(res, 0.46, 0.02, 0.46, 0.46) +
  cowplot::draw_text("(a)", x = 0.06, y = 0.97, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.48, y = 0.97, size = 25, fontface = "italic") +
  cowplot::draw_text("(c)", x = 0.06, y = 0.47, size = 25, fontface = "italic") +
  cowplot::draw_text("(d)", x = 0.48, y = 0.47, size = 25, fontface = "italic")
dev.off()







####################  PREDICTIONS WITHIN SPECIFIC EEZS


# read eez shapefile
eez <- sf::st_read(here::here("data", "eez", "World_EEZ_v11_20191118_gpkg", "eez_v11.gpkg"))



# Clip and plot predictions in eez and add violin plot to compare them
plot_predictions_in_eez(eez, "Seychelles", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Madagascar", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Sri Lanka", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Mauritius", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Chagos", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Maldives", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Reunion", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Comoros", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Oman", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Yemen", predHis, predMod, SPMod, SPHis, wio)
plot_predictions_in_eez(eez, "Somali", predHis, predMod, SPMod, SPHis, wio)


# Clip and plot predictions in eez and add violin plot to compare them (removing extrapolation zones)
plot_predictions_in_eez_extra(eez, "Seychelles", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Madagascar", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Sri Lanka", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Mauritius", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Chagos", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Maldives", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Reunion", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Comoros", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Oman", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Yemen", predHis, predMod, df_extraHis, df_extraMod, wio)
plot_predictions_in_eez_extra(eez, "Somali", predHis, predMod, df_extraHis, df_extraMod, wio)









####################  PREDICTIONS WITHIN EEZS VS HIGH SEAS


# Clip and plot predictions in high seas and add violin plot to compare them
plot_predictions_in_high_seas(eez, predHis, predMod, wio)


# Clip and plot predictions in high seas and add violin plot to compare them  (removing extrapolation zones)
plot_predictions_in_high_seas_extra(eez, predHis, predMod, df_extraHis, df_extraMod, wio)


# Clip and plot predictions in all eezs and add violin plot to compare them  (removing extrapolation zones)
plot_predictions_in_eezs_extra(eez, predHis, predMod, df_extraHis, df_extraMod, wio)


# Clip and plot predictions in high seas versus eezs with violin plot (removing extrapolation zones)
plot_predictions_high_seas_vs_eezs_extra_violin(eez, predHis, predMod, df_extraHis, df_extraMod, wio)

