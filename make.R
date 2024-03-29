#initialize see for random sampling
set.seed(1)

# load all functions in R
devtools::load_all()








####################  PREPARE ENVIRONMENT DATA -------------------------------------

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








####################  READ TRANSECTS MODERN DATA -------------------------------------

oa <- read_transect_data("OA")
remmoa <- read_transect_data("REMMOA")
noaa <- read_transect_data("NOAA")

# project transects according to custom projection
custom_proj <- sp::CRS("+proj=aeqd +lat_0=0 +lon_0=60 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
oa_proj <- project_transect(oa, custom_proj)
remmoa_proj <- project_transect(remmoa, custom_proj)
noaa_proj <- project_transect(noaa, custom_proj)







####################  READ LOGBOOK HISTORICAL DATA -------------------------------------

#read logbook data
logbook <- read_logbook_data()

#clean logbook effort data
log <- select_logbook_effort_data(logbook, raster::extent(depth.m))







####################  READ HISTORICAL AND MODERN OCCURRENCE DATA -------------------------------------

# read and clean historical occurrence data
SPHis <- select_logbook_occurrence_data(logbook, raster::extent(depth.m))
write.csv(SPHis, here::here("outputs", "occurrence_points_historical.csv"), row.names = F)

# read and clean modern occurrence data
SPMod <- read_clean_mod_data()

dim(SPMod) #219
dim(SPHis) #602



# Removing occurrences that have the same coordinates is good practice to avoid pseudoreplication.
SPHis <- SPHis[!duplicated(SPHis),]
SPMod <- SPMod[!duplicated(SPMod),]

dim(SPMod)#219
dim(SPHis)#543






####################  MAKE QUICK MAPS -------------------------------------

# map sperm whale occurrences
jpeg(here::here("outputs", "occurrences.jpeg"), width = 1100, height = 550)
par(mfrow=c(1,2), mar=c(2,2,2,6))
raster::plot(depth.m, main="Historical (n=543)", legend = FALSE)
maps::map('world', fill=T, col= "grey", add=TRUE)
points(SPHis,col="black", pch=20)
raster::plot(depth.m, main="Modern (n=219)")
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







####################  PREPARATION OF OCCURRENCES -------------------------------------

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
jpeg(here::here("outputs", "figure3.jpeg"), width = 1550, height = 1200)
cowplot::ggdraw() +
  cowplot::draw_plot(depthP, 0.005, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(slopeP, 0.34, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(di_1000mP, 0.67, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(di_seaMP, 0.005, 0.1, 0.32, 0.4) +
  cowplot::draw_plot(di_spRidP, 0.34, 0.1, 0.32, 0.4) +
  cowplot::draw_text("(a)", x = 0.03, y = 0.94, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.35, y = 0.94, size = 25, fontface = "italic") +
  cowplot::draw_text("(c)", x = 0.69, y = 0.94, size = 25, fontface = "italic") +
  cowplot::draw_text("(d)", x = 0.03, y = 0.49, size = 25, fontface = "italic") +
  cowplot::draw_text("(e)", x = 0.35, y = 0.49, size = 25, fontface = "italic")
dev.off()


pdf(here::here("outputs", "figure3.pdf"), width = 15, height = 12.5)
cowplot::ggdraw() +
  cowplot::draw_plot(depthP, 0.005, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(slopeP, 0.34, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(di_1000mP, 0.67, 0.55, 0.32, 0.4) +
  cowplot::draw_plot(di_seaMP, 0.005, 0.1, 0.32, 0.4) +
  cowplot::draw_plot(di_spRidP, 0.34, 0.1, 0.32, 0.4) +
  cowplot::draw_text("(a)", x = 0.03, y = 0.94, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.35, y = 0.94, size = 25, fontface = "italic") +
  cowplot::draw_text("(c)", x = 0.69, y = 0.94, size = 25, fontface = "italic") +
  cowplot::draw_text("(d)", x = 0.03, y = 0.49, size = 25, fontface = "italic") +
  cowplot::draw_text("(e)", x = 0.35, y = 0.49, size = 25, fontface = "italic")
dev.off()






####################  SELECTION OF BACKGROUND POINTS  -------------------------------------


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









#################### MODEL EVALUATION  -------------------------------------


# do the checkerboard partitioning
cbMod <- checkerboard_partitioning(SPMod, SPMod.z, BgMod, BgMod.z, modelStack)
cbHis <- checkerboard_partitioning(SPHis, SPHis.z, BgHis[,c("Lon", "Lat")], BgHis.z, modelStack)


# model evaluation with user partition
SPMod.x <- evaluate_model(SPMod,
                          modelStack,
                          BgMod,
                          list(fc = c("L", "Q"), rm = 1), #linear and quadratic - impose low penalty on model complexity
                          list(occs.grp = cbMod$occs.grp, bg.grp = cbMod$bg.grp),
                          "modern")

SPHis.x <- evaluate_model(SPHis,
                          modelStack,
                          BgHis[,c("Lon", "Lat")],
                          list(fc = c("L", "Q"), rm = 1), #linear and quadratic - impose low penalty on model complexity
                          list(occs.grp = cbHis$occs.grp, bg.grp = cbHis$bg.grp),
                          "historical")







####################  MODEL SELECTION AND OPTIMISATION  -------------------------------------

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









####################  MODEL COEFFICIENTS AND PARTIAL PLOTS  -------------------------------------

# Here are the non-zero coefficients in our model.
coefHis <- tibble::enframe(mod.seqHis$betas)
coefMod <- tibble::enframe(mod.seqMod$betas)

coefHis$type <- factor(c('Historical'))
coefMod$type <- factor(c('Modern'))

Coef <- rbind(coefHis, coefMod)
write.csv(Coef, here::here("outputs", "coefficients.csv"))



# plot coefficients
coefsp <- make_coefficients_plot(Coef, c(-0.004, 0.002))



# barplot version
coefBar <- make_coefficients_barplot(Coef, c(-0.004, 0.002))





# partial plots

# And these are the marginal response curves for the predictor variables with non-zero
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).


# on separate plots
jpeg(here::here("outputs", "partial_plot_historical.jpeg"), width = 1000, height = 700)
plot_partial_curves(mod.seqHis, type = "cloglog", lwd = 5, col = "#e9a3c9")
dev.off()

jpeg(here::here("outputs", "partial_plot_modern.jpeg"), width = 1000, height = 700)
plot_partial_curves(mod.seqMod, type = "cloglog", lwd = 5, col = "#a1d76a")
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
jpeg(here::here("outputs", "figure4.jpeg"), width = 1400, height = 1400)
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

pdf(here::here("outputs", "figure4.pdf"), width = 18, height = 12)
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









####################  VISUALISE EXTRAPOLATION EXTENT  -------------------------------------

# get multidimensional extrapolation extent - IMPORTANT the next two lines of code take a few minutes to run on 50 threads on a server
# but they will crash a laptop computer
df_extraMod <- get_extra_extent(c("depth",  "slope",    "di_1000m",  "di_seaM", "di_spRid"), modelStack, BgMod.z, "modern")
df_extraHis <- get_extra_extent(c("depth",  "slope",    "di_1000m",  "di_seaM", "di_spRid"), modelStack, BgHis.z, "historical")


# Plot predictions with extrapolation extent
gMod <- plot_predictions_with_extra(wio, df_predMod, "Modern", df_extraMod)
gHis <- plot_predictions_with_extra(wio, df_predHis, "Historical", df_extraHis)



# Plot predictions with extrapolation extent and mpas
gMod <- plot_predictions_with_extra_mpas(wio, df_predMod, "Modern", df_extraMod, mpa_sf)
gHis <- plot_predictions_with_extra_mpas(wio, df_predHis, "Historical", df_extraHis, mpa_sf)




# plot predictions residuals with extrapolation extent
res <- plot_predictions_residuals_with_extra(wio, predHis, predMod, df_extraMod, df_extraHis)




####################  PREDICTIONS COMPARISON AND RESIDUAL MPA EFFECT -------------------------------------



# compare predictions (removing extrapolation zones)
compare_predictions_extra(predHis, df_extraHis, predMod, df_extraMod)



# test for "residual" mpa effect (mpa vs whole region) with kruskal test (removing extrapolation zones)
df_mpaHis <- test_residual_mpa_effect_mpa_vs_region_extra(predHis,  df_extraHis, mpa_sf, "Historical")
df_mpaMod <- test_residual_mpa_effect_mpa_vs_region_extra(predMod, df_extraMod, mpa_sf, "Modern")



# density plot of prediction distributions for predictions in mpa vs whole region (removing extrapolation zones)
df_mpa <- rbind(df_mpaHis, df_mpaMod)

df_mpa <- df_mpa %>%
  dplyr::mutate(model = forcats::fct_relevel(model, levels = "Modern", "Historical"))

dens_mpa_region <- ggplot2::ggplot(data = df_mpa, ggplot2::aes(x = value, y = model, fill = model, alpha = type, linetype = type, scale = 0.9)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = c("#e9a3c9","#a1d76a")) +
  ggplot2::scale_alpha_manual(values = c(1, .2)) +
  ggplot2::xlab('Habitat suitability') +
  ggplot2::ylab("") +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "none",
                 axis.text = ggplot2::element_text(size = 18),
                 axis.title = ggplot2::element_text(size = 18)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2)

ggplot2::ggsave(here::here("outputs", "density_plot_all_models_mpa_vs_region.png"), dens_mpa_region, width = 9, height = 7)




# test for "residual" mpa effect (mpa vs out mpa) with kruskal test (removing extrapolation zones)
df_mpaHis <- test_residual_mpa_effect_mpa_vs_outmpa_extra(predHis,  df_extraHis, mpa_sf, "Historical")
df_mpaMod <- test_residual_mpa_effect_mpa_vs_outmpa_extra(predMod, df_extraMod, mpa_sf, "Modern")



# density plot of prediction distributions for predictions in mpa vs out mpa (removing extrapolation zones)
df_mpa <- rbind(df_mpaHis, df_mpaMod)
df_mpa$type <- factor(df_mpa$type , levels=c("out", "mpa"))

df_mpa <- df_mpa %>%
  dplyr::mutate(model = forcats::fct_relevel(model, levels = "Modern", "Historical"))

dens_mpa_out <- ggplot2::ggplot(data = df_mpa, ggplot2::aes(x = value, y = model, fill = model, alpha = type, linetype = type, scale = 0.9)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = c("#e9a3c9","#a1d76a")) +
  ggplot2::scale_alpha_manual(values = c(1, .2)) +
  ggplot2::xlab('Habitat suitability') +
  ggplot2::ylab("") +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "none",
                 axis.text = ggplot2::element_text(size = 18),
                 axis.title = ggplot2::element_text(size = 18)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2)

ggplot2::ggsave(here::here("outputs", "density_plot_all_models_mpa_vs_outmpa.png"), dens_mpa_out, width = 9, height = 7)




# Calculate median predictions in and out of mpas (removing extrapolation zones)
calculate_median_predictions_mpa_vs_outmpa_extra(predHis, df_extraHis, mpa_sf)
calculate_median_predictions_mpa_vs_outmpa_extra(predMod, df_extraMod, mpa_sf)



# Calculate median predictions in mpas vs in whole region (removing extrapolation zones)
calculate_median_predictions_mpa_vs_region_extra(predHis, df_extraHis, mpa_sf)
calculate_median_predictions_mpa_vs_region_extra(predMod, df_extraMod, mpa_sf)






#################### REFUGEE SPECIES EFFECT -------------------------------------


# select threshold indicating highly suitable habitat (using max SSS following Liu et al 2013)
sssMod <- select_sss_threshold(SPMod, BgMod.z, modelStack)
sssMod #0.4346055
sssHis <- select_sss_threshold(SPHis, BgHis.z, modelStack)
sssHis #0.5369416


# Barplot of nb of cells above thresholds in and out of MPAs for historical and modern predictions
barplot_inout_mpas <- barplot_predictions_in_out_mpas_above_threshold(predMod, predHis, sssMod, sssHis, mpa_sf)




# Plot predictions above sss threshold with extrapolation extent and mpas
gMod_above <- plot_predictions_with_extra_mpas_above_threshold(wio, df_predMod, "Modern", df_extraMod, mpa_sf, sssMod)
gHis_above <- plot_predictions_with_extra_mpas_above_threshold(wio, df_predHis, "Historical", df_extraHis, mpa_sf, sssHis)




# Map high suitability predictions that have become low suitability and low suitability predictions that have become high suitability
gLostGained <- plot_predictions_lost_gained(predMod, predHis, sssMod, sssHis, df_extraMod, df_extraHis, mpa_sf)


# Barplot of cell numbers for lost and gained habitat in and out of MPAs
barplotLostGained <- barplot_cells_lost_gained_in_out_mpa(predMod, predHis, sssMod, sssHis, mpa_sf)


# Mask distance to coast with lost/gained predictions
di_coast <- read_bathy_data("di_coast")
di_coast_masked_lost <- mask_dist_to_coast_predictions_lost(di_coast, predMod, predHis, sssMod, sssHis, df_extraMod, df_extraHis)
di_coast_masked_gained <- mask_dist_to_coast_predictions_gained(di_coast, predMod, predHis, sssMod, sssHis, df_extraMod, df_extraHis)


#density plot
# convert to dataframe for plotting
dflost <- as(di_coast_masked_lost, "SpatialPixelsDataFrame")
dflost <- as.data.frame(dflost)
dflost$model <- "Lost"

dfgained <- as(di_coast_masked_gained, "SpatialPixelsDataFrame")
dfgained <- as.data.frame(dfgained)
dfgained$model <- "Gained"

df <- rbind(dflost, dfgained)
df$model <- factor(df$model , levels=c("Lost", "Gained"))

densLostGained <- ggplot2::ggplot(data = df, ggplot2::aes(x = di_coast, y = model, fill = model, scale = 0.5)) +
      ggridges::geom_density_ridges(panel_scaling = FALSE) +
      ggplot2::scale_fill_manual(values = c("#fdae61", "#2b83ba")) +
      ggplot2::scale_x_continuous(breaks = seq(0, 2800, 400)) +
      ggplot2::xlab('Distance to coast (km)') +
      ggplot2::ylab("") +
      ggplot2::theme_light() +
      ggplot2::theme(legend.position = "none",
                     axis.text = ggplot2::element_text(size = 18),
                     axis.title = ggplot2::element_text(size = 18)) +
      ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2)

ggplot2::ggsave(here::here("outputs", "density_plot_dist_to_coast_predictions_lost_gained.png"), densLostGained, width = 9, height = 7)



# multiplot figures

jpeg(here::here("outputs", "figure5.jpeg"), width = 1320, height = 960)
cowplot::ggdraw() +
  cowplot::draw_plot(gHis, 0.02, 0.52, 0.46, 0.46) +
  cowplot::draw_plot(gMod, 0.46, 0.52, 0.46, 0.46) +
  cowplot::draw_plot(gHis_above, 0.02, 0.02, 0.46, 0.46) +
  cowplot::draw_plot(gMod_above, 0.46, 0.02, 0.46, 0.46) +
  cowplot::draw_text("(a)", x = 0.06, y = 0.97, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.48, y = 0.97, size = 25, fontface = "italic") +
  cowplot::draw_text("(c)", x = 0.06, y = 0.47, size = 25, fontface = "italic") +
  cowplot::draw_text("(d)", x = 0.48, y = 0.47, size = 25, fontface = "italic")
dev.off()

pdf(here::here("outputs", "figure5.pdf"), width = 16, height = 12)
cowplot::ggdraw() +
  cowplot::draw_plot(gHis, 0.02, 0.52, 0.46, 0.46) +
  cowplot::draw_plot(gMod, 0.46, 0.52, 0.46, 0.46) +
  cowplot::draw_plot(gHis_above, 0.02, 0.02, 0.46, 0.46) +
  cowplot::draw_plot(gMod_above, 0.46, 0.02, 0.46, 0.46) +
  cowplot::draw_text("(a)", x = 0.06, y = 0.97, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.48, y = 0.97, size = 25, fontface = "italic") +
  cowplot::draw_text("(c)", x = 0.06, y = 0.47, size = 25, fontface = "italic") +
  cowplot::draw_text("(d)", x = 0.48, y = 0.47, size = 25, fontface = "italic")
dev.off()

jpeg(here::here("outputs", "figure6.jpeg"), width = 1320, height = 600)
cowplot::ggdraw() +
  cowplot::draw_plot(dens_mpa_out, 0.02, 0.12, 0.46, 0.76) +
  cowplot::draw_plot(barplot_inout_mpas, 0.50, 0.12, 0.46, 0.76) +
  cowplot::draw_text("(a)", x = 0.06, y = 0.87, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.52, y = 0.87, size = 25, fontface = "italic")
dev.off()

pdf(here::here("outputs", "figure6.pdf"), width = 14, height = 9)
cowplot::ggdraw() +
  cowplot::draw_plot(dens_mpa_out, 0.02, 0.12, 0.46, 0.76) +
  cowplot::draw_plot(barplot_inout_mpas, 0.50, 0.12, 0.46, 0.76) +
  cowplot::draw_text("(a)", x = 0.06, y = 0.87, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.52, y = 0.87, size = 25, fontface = "italic")
dev.off()

jpeg(here::here("outputs", "figure7.jpeg"), width = 1320, height = 600)
cowplot::ggdraw() +
  cowplot::draw_plot(gLostGained, 0.02, 0.12, 0.46, 0.76) +
  cowplot::draw_plot(densLostGained, 0.45, 0.10, 0.46, 0.76) +
  cowplot::draw_text("(a)", x = 0.04, y = 0.87, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.48, y = 0.87, size = 25, fontface = "italic")
dev.off()

pdf(here::here("outputs", "figure7.pdf"), width = 16, height = 9)
cowplot::ggdraw() +
  cowplot::draw_plot(gLostGained, 0.04, 0.12, 0.46, 0.76) +
  cowplot::draw_plot(densLostGained, 0.51, 0.10, 0.46, 0.76) +
  cowplot::draw_text("(a)", x = 0.04, y = 0.87, size = 25, fontface = "italic") +
  cowplot::draw_text("(b)", x = 0.53, y = 0.87, size = 25, fontface = "italic")
dev.off()



#################### SENSIVITY ANALYSIS (ie sub-sampling historical occurrences to sample size of modern occurrences)  -------------------------------------

his_coef_sens_all <- list()
his_metrics_sens_all <- list()
his_pred_all <- raster::stack()

for (i in 1:10){

  print(paste("sensitivity step", i))

  ######  Random sample of historical occurrences

  indices <- sample(1:nrow(SPHis), size = nrow(SPMod), replace = F)
  SPHis_sens <- SPHis[indices,]



  ######  PREPARATION OF OCCURRENCES

  # Let's now remove occurrences that are cell duplicates -- these are
  # occurrences that share a grid cell in the predictor variable rasters.
  SPHis_sens <- remove_duplicates_occ(modelStack, SPHis_sens)

  # Extract the predictor variable values at our occurrence localities.
  SPHis_sens.z <- extract_bathy_at_point(modelStack, SPHis_sens)



  ######  MODEL EVALUATION

  # do the checkerboard partitioning
  cbHis_sens <- checkerboard_partitioning(SPHis_sens, SPHis_sens.z, BgHis[,c("Lon", "Lat")], BgHis.z, modelStack)

  # model evaluation with user partition
  SPHis_sens.x <- evaluate_model(SPHis_sens,
                            modelStack,
                            BgHis[,c("Lon", "Lat")],
                            list(fc = c("L", "Q"), rm = 1), #linear and quadratic - impose low penalty on model complexity
                            list(occs.grp = cbHis_sens$occs.grp, bg.grp = cbHis_sens$bg.grp),
                            "historical",
                            sensitivity = TRUE)



  ######  MODEL SELECTION AND OPTIMISATION

  # select best model
  his_metrics_sens <- select_best_model(SPHis_sens.x, "Historical", sensitivity = TRUE)
  his_metrics_sens$type <- factor('Historical')
  his_metrics_sens$step <- factor(i)

  # fill in metrics list
  his_metrics_sens_all[[i]] <- his_metrics_sens

  # We can select a single model from the ENMevaluation object using the tune.args of our optimal model
  his_mod_sens <- ENMeval::eval.models(SPHis_sens.x)[[his_metrics_sens$tune.args]]



  ######  MODEL COEFFICIENTS AND PARTIAL PLOTS

  # Here are the non-zero coefficients in our model.
  his_coef_sens <- tibble::enframe(his_mod_sens$betas)
  his_coef_sens$type <- factor('Historical')
  his_coef_sens$step <- factor(i)

  # fill in coef list
  his_coef_sens_all[[i]] <- his_coef_sens



  ######  MODEL PREDICTIONS

  # make predictions from selected model
  predHis_sens <- ENMeval::eval.predictions(SPHis_sens.x)[[his_metrics_sens$tune.args]]

  # stack prediction rasters
  his_pred_all <- raster::stack(his_pred_all, predHis_sens)

}




#### convert lists to dataframes and write csv

his_coef_sens_all <- do.call(rbind.data.frame, his_coef_sens_all)
write.csv(his_coef_sens_all, here::here("outputs", "historical_coefficients_sensitivity.csv"))

his_metrics_sens_all <- do.call(rbind.data.frame, his_metrics_sens_all)
write.csv(his_metrics_sens_all, here::here("outputs", "historical_metrics_best_models_sensitivity.csv"))

mean(his_metrics_sens_all$auc.val.avg) #0.6508068
sd(his_metrics_sens_all$auc.val.avg) #0.009309328




### make coefficient barplot

# preprocessing

his_coef_sens_all %>%
  dplyr::mutate(name = as.factor(name)) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(mean_value = mean(value),
                   sd = sd(value)) %>%
  dplyr::mutate(type = 'Historical') %>%
  dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename('value' = mean_value) -> his_coef_sens_all2

coefMod %>%
  dplyr::mutate(name = as.factor(name)) %>%
  dplyr::mutate(sd = 0) %>%
  dplyr::select(-type) %>%
  dplyr::mutate(type = 'Modern') %>%
  dplyr::mutate(type = as.factor(type)) -> coefMod2

coefs <- rbind(his_coef_sens_all2, coefMod2)

make_coefficients_barplot_sensitivity(coefs, c(-0.004, 0.002))






### plot mean and sd historical predictions with mpa and extrap extent

# calculate mean and sd predictions

his_pred_mean <- raster::calc(his_pred_all, fun = mean)
his_pred_sd <- raster::calc(his_pred_all, fun = sd)

# plot

plot_mean_sd_historical_predictions_with_extra_mpas(wio, his_pred_mean, his_pred_sd, df_extraHis, mpa_sf)
