

# load all functions in R
devtools::load_all()



#reload rdata
#load(here::here("make.RData"))


####################  READ PRESENCE DATA

# read and clean historical data
SPHis <- read_clean_hist_data()

# read and clean modern data
SPMod <- read_clean_mod_data()

# Removing occurrences that have the same coordinates is good practice to
# avoid pseudoreplication.
SPHis <- SPHis[!duplicated(SPHis),]
SPMod <- SPMod[!duplicated(SPMod),]







####################  READ TRANSECTS DATA (MODERN DATA)

oa <- read_transect_data("OA")
remmoa <- read_transect_data("REMMOA")
noaa <- read_transect_data("NOAA")

#project transects according to custom projection
custom_proj <- sp::CRS("+proj=aeqd +lat_0=0 +lon_0=60 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
oa_proj <- project_transect(oa, custom_proj)
remmoa_proj <- project_transect(remmoa, custom_proj)
noaa_proj <- project_transect(noaa, custom_proj)




####################  PREPARE ENVIRONMENT DATA

# read rasters of bathymetrically derived predictors - each has an ecological rational
depth.m <- read_bathy_data("depth")
slope.per <- read_bathy_data("slope")
di_1000.km <- read_bathy_data("di_1000m")
di_seaM.km <- read_bathy_data("di_seaM")
# di_guy.km <- read_bathy_data("di_guy")
# di_trough.km <- read_bathy_data("di_Trough")
# di_tren.km <- read_bathy_data("di_Tren")
# di_rid.km <- read_bathy_data("di_rid")
di_spRid.km <- read_bathy_data("di_spRid")

# stack all predictor rasters into a single stack
modelStack <- raster::stack(c(depth.m, slope.per, di_1000.km, di_seaM.km, di_spRid.km), quick=TRUE)

# reduce resolution by a factor of 10 to decrease processing time
modelStack <- raster::aggregate(modelStack, fact = 10)

# check correlations
cor <- raster::layerStats(modelStack, 'pearson', na.rm = TRUE)
write.csv(cor, here::here("outputs", "raster_correlations.csv"))

# stack again after removing
# di_rid.km due to cor >.7 with di_seaM -> di_rid removed
modelStack <- raster::stack(c(depth.m, slope.per, di_1000.km, di_seaM.km, di_spRid.km), quick=TRUE)

# reduce again resolution by a factor of 10 to decrease processing time
modelStack <- raster::aggregate(modelStack, fact = 20)










####################  QUICK MAPS

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

jpeg(here::here("outputs", "rasters.jpeg"), width = 1340, height = 960)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol=3)
dev.off()






####################  SELECTION OF BACKGROUND POINTS

## Let's now remove occurrences that are cell duplicates -- these are
# occurrences that share a grid cell in the predictor variable rasters.

SPMod <- remove_duplicates_occ(modelStack, SPMod)
SPHis <- remove_duplicates_occ(modelStack, SPHis)


# Extract the predictor variable values at our occurrence localities.

SPHis.z <- extract_bathy_at_point(modelStack, SPHis)
SPMod.z <- extract_bathy_at_point(modelStack, SPMod)


# define model extent with raster
model.extent <- raster::extent(26.25, 84.75, -40.25, 25.25) #Extent defined by the West Indian Ocean


#Historical data: get random background points in whole region
BgHis <- get_random_points(modelStack, model.extent, SPHis, 5000)



#Modern data: generate background points from observed absences on transects

#segment transects into 10km pieces

library(sp) #to ensure SegmentSpatialLines is working
oa_proj_seg <- SegmentSpatialLines(sl = oa_proj, length = 10000, merge.last = TRUE)
remmoa_proj_seg <- SegmentSpatialLines(sl = remmoa_proj, length = 10000, merge.last = TRUE)
noaa_proj_seg <- SegmentSpatialLines(sl = noaa_proj, length = 10000, merge.last = TRUE)
proj4string(oa_proj_seg) <- custom_proj
proj4string(remmoa_proj_seg) <- custom_proj
proj4string(noaa_proj_seg) <- custom_proj



#reproject to geographic projection

proj <- CRS("+proj=longlat +datum=WGS84 +no_defs")
oa_seg <- project_transect(oa_proj_seg, proj)
remmoa_seg <- project_transect(remmoa_proj_seg, proj)
noaa_seg <- project_transect(noaa_proj_seg, proj)

# plot(noaa_seg, col = rep(c(1, 2), length.out = length(noaa_seg)), axes = T)
# plot(remmoa_seg, col = rep(c(1, 2), length.out = length(remmoa_seg)), axes = T, add = T)
# plot(oa_seg, col = rep(c(1, 2), length.out = length(oa_seg)), axes = T, add = T)
# points(SPMod, col=2)



#merge all segments

sp.lines <- list()
sp.lines[[1]] <- oa_seg
sp.lines[[2]] <- remmoa_seg
sp.lines[[3]] <- noaa_seg
seg <- do.call(rbind, sp.lines)


#get observed absence points from transect segments

abs <- get_observed_absence_points(seg, SPMod)

plot(abs, cex = 0.1)
points(SPMod, col=2, cex = 0.1)


# create a raster file on the basis of survey points density
BgMod <- create_background_points_from_survey(modelStack, abs, 5000)




#once the randomPoints have been generated, we can batch-extract all variables from the raster stack
BgMod.z <- extract_bathy_at_point(modelStack, BgMod)
BgHis.z <- extract_bathy_at_point(modelStack, BgHis)



# map background points
jpeg(here::here("outputs", "background_points.jpeg"), width = 1500, height = 750)
par(mfrow=c(1,2), mar=c(2,2,2,6))
raster::plot(depth.m, main = "Historical")
maps::map('world', fill=T , col= "grey", add=TRUE)
points(BgHis, col="blue", pch=20, cex=0.1)
raster::plot(depth.m, main = "Modern")
maps::map('world', fill=T , col= "grey", add=TRUE)
points(BgMod, col="blue", pch=20, cex=0.1)
dev.off()

# write background points
write.csv(BgMod, here::here("outputs", "background_points_modern.csv"), row.names = F)
write.csv(BgHis, here::here("outputs", "background_points_historial.csv"), row.names = F)





## remove unnecessary objects

rm(abs, cor, depth.m, di_1000.km, di_guy.km, di_rid.km, di_seaM.km, di_spRid.km, di_tren.km, di_trough.km,
   noaa, noaa_proj, noaa_proj_seg, noaa_seg, oa, oa_proj, oa_proj_seg, oa_seg, remmoa, remmoa_proj, remmoa_proj_seg,
   remmoa_seg, seg, slope.per, sp.lines)





#################### MODEL EVALUATION


### do the checkerboard partitioning

cbMod <- checkerboard_partitioning(SPMod, SPMod.z, BgMod, BgMod.z, modelStack)
cbHis <- checkerboard_partitioning(SPHis, SPHis.z, BgHis, BgHis.z, modelStack)



### model evaluation with user partition

SPMod.x <- evaluate_model(SPMod,
                          modelStack,
                          BgMod,
                          list(fc = c("L"), rm = 1), #only linear fc / impose low penalty on model complexity
                          list(occs.grp = cbMod$occs.grp, bg.grp = cbMod$bg.grp),
                          "modern")

SPHis.x <- evaluate_model(SPHis,
                          modelStack,
                          BgHis,
                          list(fc = c("L"), rm = 1), #only linear fc /  impose low penalty on model complexity
                          list(occs.grp = cbHis$occs.grp, bg.grp = cbHis$bg.grp),
                          "historical")





####################  MODEL SELECTION AND OPTIMISATION

#select best model

opt.seqMod <- select_best_model(SPMod.x)
opt.seqHis <- select_best_model(SPHis.x)

#write model metrics

Model_metrics <- rbind(opt.seqHis, opt.seqMod)
Model_metrics$type <- factor(c('Historical', 'Modern'))
write.csv(Model_metrics, here::here("outputs", "model_metrics.csv"))

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.

mod.seqHis <- ENMeval::eval.models(SPHis.x)[[opt.seqHis$tune.args]]
mod.seqMod <- ENMeval::eval.models(SPMod.x)[[opt.seqMod$tune.args]]



# Here are the non-zero coefficients in our model.

coefHis <- tibble::enframe(mod.seqHis$betas)
coefMod <- tibble::enframe(mod.seqMod$betas)

coefHis$type <- factor(c('Historical'))
coefMod$type <- factor(c('Modern'))

Coef <- rbind(coefHis, coefMod)
write.csv(Coef, here::here("outputs", "coefficients.csv"))

# plot coefficients

coefsp <- ggplot2::ggplot(Coef, ggplot2::aes(x=name, y=value, color = type)) +
  ggplot2::geom_point() +
  ggforce::facet_zoom(ylim = c(-0.0043, 0.002)) +
  ggplot2::theme_light() +
  ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="right", text = ggplot2::element_text(size=16),
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::xlab("Predictors") +
  ggplot2::ylab("Coefficient") +
  ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "black", size=.5) +
  ggplot2::scale_colour_manual(values = c("#00BA38", "#F8766D"), guide = ggplot2::guide_legend(ncol = 1, direction = "horizontal",
                                                                                                          label.position="left", label.hjust = 0.5, label.vjust = 0,
                                                                                                          label.theme = ggplot2::element_text(angle = 0)))

ggplot2::ggsave(here::here("outputs", "coefficients_plot.png"), coefsp, width = 10, height = 7)

# barplot version

coefBar <- ggplot2::ggplot(Coef, ggplot2::aes(x=name, y=value, fill = type)) +
  ggplot2::geom_bar(position="dodge", stat="identity") +
  ggforce::facet_zoom(ylim = c(-0.0043, 0.002)) +
  ggplot2::theme_light() +
  ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="right", text = ggplot2::element_text(size=16),
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::xlab("Predictors") +
  ggplot2::ylab("Coefficient") +
  ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "black", size=.5) +
  ggplot2::scale_fill_manual(values = c("#00BA38", "#F8766D"), guide = ggplot2::guide_legend(ncol = 1, direction = "horizontal",
                                                                                                        label.position="left", label.hjust = 0.5, label.vjust = 0,
                                                                                                        label.theme = ggplot2::element_text(angle = 0)))

ggplot2::ggsave(here::here("outputs", "coefficients_barplot.png"), coefBar, width = 10, height = 7)












# And these are the marginal response curves for the predictor variables with non-zero
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

# curves on separate plots
jpeg(here::here("outputs", "partial_plot_historical.jpeg"), width = 1000, height = 700)
plot_partial_curves(mod.seqHis, type = "cloglog", lwd = 5, col = "#F8766D")
dev.off()

jpeg(here::here("outputs", "partial_plot_modern.jpeg"), width = 1000, height = 700)
plot_partial_curves(mod.seqMod, type = "cloglog", lwd = 5, col = "#00BA38")
dev.off()



# curves on the same plot

dat_his = plot_partial_curves(mod.seqHis, type = "cloglog")
dat_mod = plot_partial_curves(mod.seqMod, type = "cloglog")

p1 <- ggplot2::ggplot(dat_his$depth, ggplot2::aes(y = pred, x = depth)) +
        ggplot2::geom_line(color = "#00BA38", size = 2) +
        ggplot2::geom_line(data = dat_mod$depth, color = "#F8766D",  size = 2) +
        ggplot2::ylim(0,1) +
        ggplot2::ylab("predictions") +
       ggplot2::scale_x_continuous(limits=c(min(dat_his$depth, dat_mod$depth), 0)) +
        ggplot2::theme_bw()

p2 <- ggplot2::ggplot(dat_his$slope, ggplot2::aes(y = pred, x = slope)) +
        ggplot2::geom_line(color = "#00BA38", size = 2) +
        ggplot2::geom_line(data = dat_mod$slope, color = "#F8766D",  size = 2) +
        ggplot2::ylim(0,1) +
        ggplot2::ylab("predictions") +
        ggplot2::scale_x_continuous(limits=c(0, max(dat_his$slope, dat_mod$slope))) +
        ggplot2::theme_bw()

p3 <- ggplot2::ggplot(dat_his$di_1000m, ggplot2::aes(y = pred, x = di_1000m)) +
        ggplot2::geom_line(color = "#00BA38", size = 2) +
        ggplot2::geom_line(data = dat_mod$di_1000m, color = "#F8766D",  size = 2) +
        ggplot2::ylim(0,1) +
        ggplot2::ylab("predictions") +
        ggplot2::scale_x_continuous(limits=c(0, max(dat_his$di_1000m, dat_mod$di_1000m))) +
        ggplot2::theme_bw()

p4 <- ggplot2::ggplot(dat_his$di_seaM, ggplot2::aes(y = pred, x = di_seaM)) +
        ggplot2::geom_line(color = "#00BA38", size = 2) +
        ggplot2::geom_line(data = dat_mod$di_seaM, color = "#F8766D",  size = 2) +
        ggplot2::ylim(0,1) +
        ggplot2::ylab("predictions") +
        ggplot2::scale_x_continuous(limits=c(0, max(dat_his$di_seaM, dat_mod$di_seaM))) +
        ggplot2::theme_bw()

p5 <- ggplot2::ggplot(dat_his$di_spRid, ggplot2::aes(y = pred, x = di_spRid)) +
        ggplot2::geom_line(color = "#00BA38", size = 2) +
        ggplot2::geom_line(data = dat_mod$di_spRid, color = "#F8766D",  size = 2) +
        ggplot2::ylim(0,1) +
        ggplot2::ylab("predictions") +
        ggplot2::scale_x_continuous(limits=c(0, max(dat_his$di_spRid, dat_mod$di_spRid))) +
        ggplot2::theme_bw()


jpeg(here::here("outputs", "partial_plots.jpeg"), width = 1000, height = 700)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, coefBar, ncol=3)
dev.off()







# make predictions from selected model
predMod <- ENMeval::eval.predictions(SPMod.x)[[opt.seqMod$tune.args]]
predHis <- ENMeval::eval.predictions(SPHis.x)[[opt.seqHis$tune.args]]












#comparison of  predictions

#Pearsons's correlation coefficients between predictions

predStack <- raster::stack(c(predHis, predMod))
names(predStack) <- c("historical", "modern")
jnk <- raster::layerStats(predStack, 'pearson', na.rm=T)
corr_matrix <- jnk$'pearson correlation coefficient'
corr_matrix






#plot predictions and histogram

### get contries for plotting
wio <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") #world map with 50m resolution


### convert prediction rasters to dataframes for plotting

df_predHis <- convert_predictions_to_df(predHis, "historical")
df_predMod <- convert_predictions_to_df(predMod, "modern")


# bind prediction dataframes
df_pred <- rbind(df_predHis, df_predMod)




### MPA values for residual MPA testing

## read and plot mpas
mpa_sf <- sf::read_sf(here::here("data", "mpa_sf2", "mpa_sf.shp"), stringsAsFactors=TRUE)


### test for "residual" mpa effect with kruskall test

df_mpaHis <- test_residual_mpa_effect(predHis, df_predHis, mpa_sf, "historical")
df_mpaMod <- test_residual_mpa_effect(predMod, df_predMod, mpa_sf, "modern")

df_mpa <- rbind(df_mpaHis, df_mpaMod)
df_mpa$type <- factor(df_mpa$type , levels=c("region", "mpa"))
summary(df_mpa)

### wilcox test (non-parametric alternative to paired t-test used to compare paired data)
pairwise.wilcox.test(df_mpa$value, df_mpa$type, p.adjust.method = 'BH')



### density plot of prediction distributions

df_mpa <- df_mpa %>%
  dplyr::mutate(model = forcats::fct_relevel(model, levels = "modern", "historical"))

d <- ggplot2::ggplot(data = df_mpa, ggplot2::aes(x = value, y = model, fill = model, alpha = type, linetype = type, scale = 0.9)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = c("#F8766D","#00BA38", "#619CFF")) +
  ggplot2::scale_alpha_manual(values = c(0.6, .1)) +
  ggplot2::xlab('Habitat suitability') +
  ggplot2::ylab("") +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "none", axis.text = ggplot2::element_text(size=12)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2)

ggplot2::ggsave(here::here("outputs", "density_plot_all_models.png"), d, width = 9, height = 7)





### plot predictions

gMod <- plot_predictions(wio, df_predMod, "modern")
gHis <- plot_predictions(wio, df_predHis, "historical")


### plot predictions difference

d <- plot_predictions_difference(wio, df_predHis, df_predMod)




### get multidimensional extrapolation extent - the next two lines took a few minutes to run on 50 threads

df_extraMod <- get_extra_extent(c("depth",  "slope",    "di_1000m",  "di_seaM", "di_spRid"), modelStack, BgMod.z, "modern")
df_extraHis <- get_extra_extent(c("depth",  "slope",    "di_1000m",  "di_seaM", "di_spRid"), modelStack, BgHis.z, "historical")



# Plot predictions with extrapolation extent

gMod2 <- plot_predictions_with_extra(wio, df_predMod, "modern", df_extraMod)
gHis2 <- plot_predictions_with_extra(wio, df_predHis, "historical", df_extraHis)



# Plot predictions with extrapolation extent and mpas

gMod2 <- plot_predictions_with_extra_mpas(wio, df_predMod, "modern", df_extraMod, mpa_sf)
gHis2 <- plot_predictions_with_extra_mpas(wio, df_predHis, "historical", df_extraHis, mpa_sf)



### save objects
save(df_predMod, df_extraMod, df_predHis, df_extraHis, BgMod.z, BgHis.z, predHis, predMod, SPMod.x, SPHis.x, opt.seqMod, opt.seqHis, file = here::here("make.RData"))


#------------------------------------------------------------------------------------------------------


#read eez shapefile
eez <- sf::st_read(here::here("data", "eez", "World_EEZ_v11_20191118_gpkg", "eez_v11.gpkg"))



#Clip and plot predictions in eez and add violin plot to compare them
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


plot_predictions_in_high_seas(eez, predHis, predMod, wio)
plot_predictions_in_high_seas_extra(eez, predHis, predMod, df_extraHis, df_extraMod, wio)

plot_predictions_in_eezs_extra(eez, predHis, predMod, df_extraHis, df_extraMod, wio)



##extract predictions along transect (need to mask extrapolation zones)
#https://rdrr.io/cran/inlmisc/man/ExtractAlongTransect.html

x <- c(20, 90)
y <- c(-10, -10)
transect <- sp::SpatialLines(list(sp::Lines(sp::Line(cbind(x,y)), ID="a")))
sp::proj4string(l) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

t <- inlmisc::ExtractAlongTransect(transect, predMod)



segs <- inlmisc::ExtractAlongTransect(transect, predMod)
for (i in seq_along(segs)) points(segs[[i]])

dev.new()
xlab <- "Distance along transect"
ylab <- "Raster value"
xlim <- range(vapply(segs, function(i) {
  range(i@data[, "dist"])
}, c(0, 0)))
ylim <- range(vapply(segs, function(i) {
  range(i@data[, "fc.L_rm.1"], na.rm = TRUE)
}, c(0, 0)))
inlmisc::PlotGraph(NA, xlab = xlab, ylab = ylab,
          xlim = xlim, ylim = ylim, type = "n")
cols <- inlmisc::GetColors(length(segs), scheme = "bright")
for (i in seq_along(segs))
  lines(segs[[i]]@data[, c("dist", "fc.L_rm.1")],
        col = cols[i], lwd = 2)
coords <- sp::coordinates(transect)
n <- length(transect)
d <- cumsum(c(0, as.matrix(dist((coords)))[cbind(1:(n - 1), 2:n)]))
abline(v = d, lty = 2)
mtext(sprintf("(%d, %d)", coords[1, 1], coords[1, 2]),
      line = -1, adj = 0, cex = 0.7)
mtext(sprintf("(%d, %d)", coords[n, 1], coords[n, 2]),
      line = -1, adj = 1, cex = 0.7)

