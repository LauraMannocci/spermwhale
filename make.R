

# load all functions in R
devtools::load_all()


####################  READ PRESENCE DATA

# read and clean historical data
SPHis <- read_clean_hist_data()

# read and clean modern data
SPMod <- read_clean_mod_data()

# bind historical and modern data
SPCom <- rbind(SPHis, SPMod)

# Removing occurrences that have the same coordinates is good practice to
# avoid pseudoreplication.
SPHis <- SPHis[!duplicated(SPHis),]
SPMod <- SPMod[!duplicated(SPMod),]
SPCom <- SPCom[!duplicated(SPCom),]





####################  READ TRANSECTS DATA (MODERN DATA)

oa <- read_transect_data("OA")
remmoa <- read_transect_data("REMMOA")
noaa <- read_transect_data("NOAA")

# maps::map('world', fill=T, col= "grey", xlim = c(30, 90), ylim = c(-25, 30))
# sp::plot(oa, add = T)
# sp::plot(remmoa, add = T, col = 2)
# sp::plot(noaa, add = T, col = 3)






####################  PREPARE ENVIRONMENT DATA

# read rasters of bathymetrically derived predictors
depth.m <- read_bathy_data("depth")
slope.per <- read_bathy_data("slope")
di_200.km <- read_bathy_data("di_200m")
di_1000.km <- read_bathy_data("di_1000m")
di_2500.km <- read_bathy_data("di_2500m")
di_5000.km <- read_bathy_data("di_5000m")
di_coast.km <- read_bathy_data("di_coast")
di_seaM.km <- read_bathy_data("di_seaM")
di_guy.km <- read_bathy_data("di_guy")
di_trough.km <- read_bathy_data("di_Trough")
di_tren.km <- read_bathy_data("di_Tren")
di_rid.km <- read_bathy_data("di_rid")
di_spRid.km <- read_bathy_data("di_spRid")

# stack all predictor rasters into a single stack
modelStack <- raster::stack(c(depth.m, slope.per, di_200.km, di_1000.km, di_2500.km, di_5000.km, di_coast.km, di_seaM.km, di_guy.km, di_trough.km, di_tren.km, di_rid.km, di_spRid.km), quick=TRUE)

# reduce resolution by a factor of 10 to decrease processing time
modelStack <- raster::aggregate(modelStack, fact = 10)

# check correlations
cor <- raster::layerStats(modelStack, 'pearson', na.rm = TRUE)
write.csv(cor, here::here("outputs", "raster_correlations.csv"))

# stack again after removing
# di_200.km due to cor >.7 with di_coast
# di_rid.km due to cor >.7 with di_seaM
# di_1000.km due to cor >.6 with di_coast
# di_5000.km due to cor >.6 with depth
# modelStack <- raster::stack(c(depth.m, slope.per, di_1000.km, di_2500.km, di_5000.km, di_coast.km, di_seaM.km, di_guy.km, di_trough.km, di_tren.km, di_spRid.km), quick=TRUE)
modelStack <- raster::stack(c(depth.m, slope.per, di_2500.km, di_coast.km, di_seaM.km, di_guy.km, di_trough.km, di_tren.km, di_spRid.km), quick=TRUE)

# reduce again resolution by a factor of 10 to decrease processing time
modelStack <- raster::aggregate(modelStack, fact = 10)










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
jpeg(here::here("outputs", "rasters.jpeg"), width = 960, height = 960)
par(mar=c(0,0,1,0))
raster::plot(modelStack)
dev.off()









####################  SELECTION OF BACKGROUND POINTS

## Let's now remove occurrences that are cell duplicates -- these are
# occurrences that share a grid cell in the predictor variable rasters.

SPCom <- remove_duplicates_occ(modelStack, SPCom)
SPMod <- remove_duplicates_occ(modelStack, SPMod)
SPHis <- remove_duplicates_occ(modelStack, SPHis)


# Extract the predictor variable values at our occurrence localities.

SPCom.z <- extract_bathy_at_point(modelStack, SPCom)
SPHis.z <- extract_bathy_at_point(modelStack, SPHis)
SPMod.z <- extract_bathy_at_point(modelStack, SPMod)




# define model extent with raster
model.extent <- raster::extent(26.25, 84.75, -40.25, 25.25) #Extent defined by the West Indian Ocean

#then, we take 10000 random background points
#these are our background points, our absence points
set.seed(0)

#Combined data
BgCom <- get_random_points(modelStack, model.extent, SPCom, 10000)

#Historical data (same as combined)
BgHis <- get_random_points(modelStack, model.extent, SPHis, 10000)



#Modern data: generate background points from transects

#define points spaced by 50km along transect lines for each survey
oa_pts <- define_points_along_survey_transects(oa, 50)
remmoa_pts <- define_points_along_survey_transects(remmoa, 50)
noaa_pts <- define_points_along_survey_transects(noaa, 50)

#merge all survey pts
pts = rbind(oa_pts, remmoa_pts, noaa_pts)
pts  = data.frame(Lon = sp::coordinates(pts)[,1], Lat = sp::coordinates(pts)[,2])

# create a raster file on the basis of survey points density
BgMod <- create_background_points_from_survey(modelStack, pts, 10000)


#once the randomPoints have been generated, we can batch-extract all variables from the raster stack
BgCom.z <- extract_bathy_at_point(modelStack, BgCom)
BgHis.z <- extract_bathy_at_point(modelStack, BgHis)
BgMod.z <- extract_bathy_at_point(modelStack, BgMod)


# map background points
jpeg(here::here("outputs", "background_points.jpeg"), width = 1500, height = 750)
par(mfrow=c(2,2), mar=c(2,2,2,6))
raster::plot(depth.m, main = "Historical")
maps::map('world', fill=T , col= "grey", add=TRUE)
points(BgHis, col="blue", pch=20, cex=0.1)
raster::plot(depth.m, main = "Modern")
maps::map('world', fill=T , col= "grey", add=TRUE)
points(BgMod, col="blue", pch=20, cex=0.1)
raster::plot(depth.m, main = "Combined")
maps::map('world', fill=T , col= "grey", add=TRUE)
points(BgCom, col="blue", pch=20, cex=0.1)
dev.off()

# write background points
write.csv(BgMod, here::here("outputs", "background_points_modern.csv"), row.names = F)
write.csv(BgHis, here::here("outputs", "background_points_historial.csv"), row.names = F)
write.csv(BgCom, here::here("outputs", "background_points_combined.csv"), row.names = F)








#################### MODEL EVALUATION


### do the checkerboard partitioning

cbCom <- checkerboard_partitioning(SPCom, SPCom.z, BgCom, BgCom.z, modelStack)
cbMod <- checkerboard_partitioning(SPMod, SPMod.z, BgMod, BgMod.z, modelStack)
cbHis <- checkerboard_partitioning(SPHis, SPHis.z, BgHis, BgHis.z, modelStack)



# create rasters limited to range of variables in background points
modelStack_BgMod <- limit_rasters_to_background_points(modelStack, BgMod.z, "modern")
modelStack_BgHis <- limit_rasters_to_background_points(modelStack, BgHis.z, "historical")





### model evaluation with user partition
SPCom.x <- evaluate_model(SPCom,
                          modelStack,
                          BgCom,
                          list(fc = c("L"), rm = 1:5), #only linear fc
                          list(occs.grp = cbCom$occs.grp, bg.grp = cbCom$bg.grp),
                          "combined")


#modelStack_BgMod <- raster::dropLayer(modelStack_BgMod, 5)

SPMod.x <- evaluate_model(SPMod,
                          modelStack,
                          # modelStack_BgMod, # use limited rasters for modern
                          BgMod,
                          list(fc = c("L"), rm = 1:5), #only linear fc
                          list(occs.grp = cbMod$occs.grp, bg.grp = cbMod$bg.grp),
                          "modern")

SPHis.x <- evaluate_model(SPHis,
                          modelStack,
                          BgHis,
                          list(fc = c("L"), rm = 1:5), #only linear fc
                          list(occs.grp = cbHis$occs.grp, bg.grp = cbHis$bg.grp),
                          "historical")





####################  MODEL SELECTION AND OPTIMISATION

#select best model

opt.seqCom <- select_best_model(SPCom.x)
opt.seqMod <- select_best_model(SPMod.x)
opt.seqHis <- select_best_model(SPHis.x)

#write model metrics

Model_metrics <- rbind(opt.seqCom, opt.seqHis, opt.seqMod)
Model_metrics$type <- factor(c('Combined', 'Historical', 'Modern'))
write.csv(Model_metrics, here::here("outputs", "model_metrics.csv"))

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.

mod.seqCom <- ENMeval::eval.models(SPCom.x)[[opt.seqCom$tune.args]]
mod.seqHis <- ENMeval::eval.models(SPHis.x)[[opt.seqHis$tune.args]]
mod.seqMod <- ENMeval::eval.models(SPMod.x)[[opt.seqMod$tune.args]]



# Here are the non-zero coefficients in our model.

coefCom <- tibble::enframe(mod.seqCom$betas)
coefHis <- tibble::enframe(mod.seqHis$betas)
coefMod <- tibble::enframe(mod.seqMod$betas)

coefCom$type <- factor(c('Combined'))
coefHis$type <- factor(c('Historical'))
coefMod$type <- factor(c('Modern'))

Coef <- rbind(coefCom, coefHis, coefMod)
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
  ggplot2::scale_colour_manual(values = c("#619CFF", "#00BA38", "#F8766D"), guide = ggplot2::guide_legend(ncol = 1, direction = "horizontal",
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
  ggplot2::scale_fill_manual(values = c("#619CFF", "#00BA38", "#F8766D"), guide = ggplot2::guide_legend(ncol = 1, direction = "horizontal",
                             label.position="left", label.hjust = 0.5, label.vjust = 0,
                             label.theme = ggplot2::element_text(angle = 0)))

ggplot2::ggsave(here::here("outputs", "coefficients_barplot.png"), coefBar, width = 10, height = 7)






# predictor contributions (requires to rebuild a maxent model)
contribMod <- get_predictors_contributions(modelStack, SPMod, BgMod, "modern")
contribHis <- get_predictors_contributions(modelStack, SPHis, BgHis, "historical")
contribCom <- get_predictors_contributions(modelStack, SPCom, BgCom, "combined")



# predictor contributions - all models
dat <- rbind(contribMod, contribHis, contribCom)

g <- ggplot2::ggplot(dat, ggplot2::aes(x = reorder(var, contrib), y = contrib, fill = type)) +
  ggplot2::geom_bar(stat="identity", position="dodge") +
  ggplot2::coord_flip() +
  ggplot2::ggtitle("all models") +
  ggplot2::ylab("Contributions (%)") +
  ggplot2::xlab("") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

ggplot2::ggsave(here::here("outputs", "predictor_contributions_all_models.png"), g, width = 9, height = 7)




# predictor importance (requires to rebuild a maxent model)
impMod <- get_predictors_importance(modelStack, SPMod, BgMod, "modern")
impHis <- get_predictors_importance(modelStack, SPHis, BgHis, "historical")
impCom <- get_predictors_importance(modelStack, SPCom, BgCom, "combined")

# predictor importance - all models
dat <- rbind(impMod, impHis, impCom)

g <- ggplot2::ggplot(dat, ggplot2::aes(x = reorder(var, contrib), y = contrib, fill = type)) +
  ggplot2::geom_bar(stat="identity", position="dodge") +
  ggplot2::coord_flip() +
  ggplot2::ggtitle("all models") +
  ggplot2::ylab("Importance (%)") +
  ggplot2::xlab("") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

ggplot2::ggsave(here::here("outputs", "predictor_importance_all_models.png"), g, width = 9, height = 7)









# And these are the marginal response curves for the predictor variables with non-zero
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

#par(mfrow=c(2,11), mar=c(0.1,0.1,0.1,0.1))

jpeg(here::here("outputs", "partial_plot_combined.jpeg"), width = 1000, height = 700)
plot(mod.seqCom, type = "cloglog")
dev.off()

jpeg(here::here("outputs", "partial_plot_historical.jpeg"), width = 1000, height = 700)
plot(mod.seqHis, type = "cloglog")
dev.off()

jpeg(here::here("outputs", "partial_plot_modern.jpeg"), width = 1000, height = 700)
plot(mod.seqMod, type = "cloglog")
dev.off()




# make predictions from selected model
predMod <- ENMeval::eval.predictions(SPMod.x)[[opt.seqMod$tune.args]] #limited raster
predHis <- ENMeval::eval.predictions(SPHis.x)[[opt.seqHis$tune.args]]
predCom <- ENMeval::eval.predictions(SPCom.x)[[opt.seqCom$tune.args]]












#comparison of  predictions

#Pearsons's correlation coefficients between predictions

predStack <- raster::stack(c(predCom, predHis, predMod))
names(predStack) <- c("combined", "historical", "modern")
jnk <- raster::layerStats(predStack, 'pearson', na.rm=T)
corr_matrix <- jnk$'pearson correlation coefficient'
corr_matrix






#plot predictions and histogram

### get contries for plotting
wio <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") #world map with 50m resolution


### convert prediction rasters to dataframes for plotting

df_predCom <- convert_predictions_to_df(predCom, "combined")
df_predHis <- convert_predictions_to_df(predHis, "historical")
df_predMod <- convert_predictions_to_df(predMod, "modern")


# bind prediction dataframes
df_pred <- rbind(df_predCom, df_predHis, df_predMod)




### MPA values for residual MPA testing

## read and plot mpas
mpa_sf <- sf::read_sf(here::here("data", "mpa_sf2", "mpa_sf.shp"), stringsAsFactors=TRUE)

ggMPA <- ggMPA +
  ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
  ggplot2::theme(legend.position = "none", axis.title = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = 'white'))



### test for "residual" mpa effect

df_mpaCom <- test_residual_mpa_effect(predCom, df_predCom, mpa_sf, "combined")
df_mpaHis <- test_residual_mpa_effect(predHis, df_predHis, mpa_sf, "historical")
df_mpaMod <- test_residual_mpa_effect(predMod, df_predMod, mpa_sf, "modern")

df_mpa <- rbind(df_mpaCom, df_mpaHis, df_mpaMod)
df_mpa$type <- factor(df_mpa$type , levels=c("region", "mpa"))
summary(df_mpa)

### wilcox test
pairwise.wilcox.test(df_mpa$value, df_mpa$type, p.adjust.method = 'BH')



### density plot

df_mpa <- df_mpa %>%
  dplyr::mutate(model = forcats::fct_relevel(model, levels = "modern", "historical", "combined"))

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
gCom <- plot_predictions(wio, df_predCom, "combined")





#plot residuals with respect to combined model

ggHis <- plot_residuals_wrt_combined(predHis, predCom, wio, "historical")
ggMod <- plot_residuals_wrt_combined(predMod, predCom, wio, "modern")





### summary plot

sum <- cowplot::ggdraw() +
  cowplot::draw_plot(d, 0, .66, .5, .33)+
  cowplot::draw_plot(gCom, .46, .66, .55, .33) +
  cowplot::draw_plot(gHis, 0, .33, .5, .33) +
  cowplot::draw_plot(gMod, .45, .33, .5, .33) +
  cowplot::draw_plot(ggHis, 0, 0, .5, .33) +
  cowplot::draw_plot(ggMod, .46, 0, .55, .33) +
  cowplot::draw_plot_label(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), c(0, 0.5, 0, 0.5, 0, .5), c(1, 1, .69, .69, .36, .36), font = "italic", size = 15)
sum


jpeg(here::here("outputs", "summary_plot.jpeg"), width = 2300, height = 3000, res = 300)
sum
dev.off()











##comparison of predictor values between his and mod occurrences

SPHis.z$type <- factor(c('Historical'))
SPMod.z$type <- factor(c('Modern'))

Mod.Hist <- rbind(SPHis.z, SPMod.z)
summary(Mod.Hist)


## violin plot with significance based on kruskal test
depthP <- predictor_violin(Mod.Hist, "depth", "Seabed depth (m)")
slopeP <- predictor_violin(Mod.Hist, "slope", "Seabed slope (%)")
di_TroughP <- predictor_violin(Mod.Hist, "di_Trough", "Distance to trough (km)")
di_TrenP <- predictor_violin(Mod.Hist, "di_Tren", "Distance to trench (km)")
di_spRidP <- predictor_violin(Mod.Hist, "di_spRid", "Distance to spreading ridge (km)")
di_seaMP <- predictor_violin(Mod.Hist, "di_seaM", "Distance to seamount (km)")
di_guyP <- predictor_violin(Mod.Hist, "di_guy", "Distance to guyot (km)")
di_coastP <- predictor_violin(Mod.Hist, "di_coast", "Distance to coast (km)")
di_5000mP <- predictor_violin(Mod.Hist, "di_5000m", "Distance to 5,000 m contour (km)")
di_2500mP <- predictor_violin(Mod.Hist, "di_2500m", "Distance to 2,500 m contour (km)")
di_1000mP <- predictor_violin(Mod.Hist, "di_1000m", "Distance to 1,000 m contour (km)")


p <- ggpubr::ggarrange(depthP, slopeP, di_TroughP,
                       di_1000mP, di_2500mP, di_5000mP,
                       di_coastP, di_guyP, di_seaMP, di_spRidP, di_TrenP,
                       ncol = 3, nrow = 4, labels = list("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)"), font.label = list(face = "italic", size =14))

jpeg(here::here("outputs", "predictors_violin_all.jpeg"), width = 2300, height = 3180, res = 300)
p
dev.off()










#plot habitat maps >10th percentile and turn into 1s to map consistent habitat
#select 10th percentiles of predictions from the optimal models, look for overlap. and difference


predCqfoc <-predC
qC <- quantile(predCqfoc, 0.9)
predCqfoc[predCqfoc < qC] <- 0
predCqfoc[predCqfoc > qC] <- 1
predCqfoc[predCqfoc < qC] <- NA
plot(predCqfoc)

predHqfoc <-predH
qH <- quantile(predHqfoc, 0.9)
predHqfoc[predHqfoc < qH] <- 0
predHqfoc[predHqfoc > qH] <- 1
predHqfoc[predHqfoc < qH] <- NA
plot(predHqfoc)

predMqfoc <-predM
qM <- quantile(predMqfoc, 0.9)
predMqfoc[predMqfoc < qM] <- 0
predMqfoc[predMqfoc > qM] <- 1
predMqfoc[predMqfoc < qM] <- NA


df_predCqfoc<- as(predCqfoc, "SpatialPixelsDataFrame")
df_predCqfoc <- as.data.frame(df_predCqfoc)
colnames(df_predCqfoc) <- c("value", "x", "y")

df_predHqfoc<- as(predHqfoc, "SpatialPixelsDataFrame")
df_predHqfoc <- as.data.frame(df_predHqfoc)
colnames(df_predHqfoc) <- c("value", "x", "y")

df_predMqfoc<- as(predMqfoc, "SpatialPixelsDataFrame")
df_predMqfoc <- as.data.frame(df_predMqfoc)
colnames(df_predMqfoc) <- c("value", "x", "y")




## plot range maps

ggCqfoc <-ggplot() + geom_sf(data = wio, fill ="black", colour = "black")+
  geom_tile(data=df_predCqfoc, alpha =0.6, aes(x=x, y=y), fill= "#619CFF")+ coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE)+ ylab("")+ xlab("")+ theme_light()+
  theme(legend.position = "none", axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
ggCqfoc

ggHqfoc <-ggplot() + geom_sf(data = wio, fill ="black", colour = "black")+
  geom_tile(data=df_predHqfoc, alpha =0.6, aes(x=x, y=y), fill= "#00BA38")+ theme_light()+  coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE)+ylab("")+xlab("")+
  theme(legend.position = "none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
ggHqfoc

ggMqfoc <-ggplot() + geom_sf(data = wio, fill ="black", colour = "black")+
  geom_tile(data=df_predMqfoc, alpha =0.6,aes(x=x, y=y), fill = "#F8766D")+theme_light()+coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE)+ylab("")+xlab("")+
  theme(legend.position = "none", axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
ggMqfoc





#test whether spermwhales are a refugee species
#need to see whether high quality habitat in the present and historical are marginal or suboptimal quality in the  combined model

##Start with high quality habitat in the modern

predMref <-predM
qM <- quantile(predMref, 0.9)
predMref[predMref < qM] <- NA
plot(predMref)
#and historical

predHref <-predH
qH <- quantile(predHref, 0.9)
predHref[predHref < qH] <- NA
plot(predHref)

##generate high quality habitat in the combined (top 10th percentile) for reference

predCref <-predC
qC <- quantile(predCref, 0.9)
predCref[predCref < qC] <- NA
df_predCref<- as(predCref, "SpatialPixelsDataFrame")
df_predCref <- as.data.frame(df_predCref)
colnames(df_predCref) <- c("value", "x", "y")
df_predCref$model <- factor(c('Combined'))

summary(df_predCref)



#For the locations of high quality habitat in the combined, extract values in Historical model

CHref <- mask(predH, predCref)
df_CHref<- as(CHref, "SpatialPixelsDataFrame")
df_CHref <- as.data.frame(df_CHref)
colnames(df_CHref) <- c("value", "x", "y")

#For the locations of high quality habitat in the combined, extract values in Modern model

CMref <- mask(predM, predCref)
df_CMref<- as(CMref, "SpatialPixelsDataFrame")
df_CMref <- as.data.frame(df_CMref)
colnames(df_CMref) <- c("value", "x", "y")

CHref <- mask(predC, predHref)
df_CHref<- as(CHref, "SpatialPixelsDataFrame")
df_CHref <- as.data.frame(df_CHref)
colnames(df_CHref) <- c("value", "x", "y")

#combined
df_CMref$model <- factor(c('Modern'))
df_CHref$model <- factor(c('Historical'))

summary(df_predCref)
summary(df_CHref)

df_CHMref<-rbind(df_predCref, df_CMref, df_CHref)

summary(df_CHMref)

kruskal.test(value ~ model, data = df_CHMref)

pairwise.wilcox.test(df_CHMref$value, df_CHMref$model,
                     p.adjust.method = "BH")


#df_predCreg <-subset(df_predC, select = -c(value,x,y))
#df_predCreg$model <- factor(c('Combined ref'))

df_CHMref <- df_CHMref %>%
  mutate(model = fct_relevel(model, levels = "Combined","Historical", "Modern"))

CHMrefDen <- ggplot(data = df_CHMref, aes(x = value, y=model, fill = model, scale =0.9)) +
  geom_violin()+geom_boxplot(width =.2) +
  scale_fill_manual(values = c( "#619CFF", "#00BA38","#F8766D"))+
  xlab('Habitat suitability')+ylab("") +
  theme_light()+ theme(legend.position = "none", axis.text=element_text(size=12))+coord_flip()
  #stat_density_ridges(quantile_lines = TRUE, quantiles = 2)
CHMrefDen



HMbp <- ggdraw()+ draw_plot(ggCqfoc, 0.03, 0, .3, 1)+
  draw_plot(ggHqfoc, .36, 0, .3, 1) +
  draw_plot(ggMqfoc, .69, 0, .3, 1) +
  draw_plot_label(c("(a)", "(b)", "(c)"), c(0, 0.33, 0.66), c(1, 1, 1), font = "italic", size = 15)

jpeg("HMbp.jpeg", width = 2500, height =950, res = 300)
HMbp
dev.off()

### refugee species plotting


summary(df_predCref)
summary(df_CMref)
summary(df_CHref)

gg_predCref <-ggplot() + geom_sf(data = wio, fill = 'black', colour = 'black')+
  geom_tile(data=df_predCref, aes(x=x, y=y, fill=value), alpha=0.8)+
  scale_fill_viridis_c(limits = c(0.2, 1), option = "viridis")+theme_light()+
  coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE)+ylab("") + xlab('')+
  theme(legend.position = "none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
#gg_predCref

gg_CHref <-ggplot() + geom_sf(data = wio, fill = 'black', colour = 'black')+
  geom_tile(data=df_CHref, aes(x=x, y=y, fill=value), alpha=0.8)+
  scale_fill_viridis_c(limits = c(0.2, 1), option = "viridis")+theme_light()+
  coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE)+ylab("") + xlab('')+
  theme(legend.position ="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
#gg_CHref

gg_CMref <-ggplot() + geom_sf(data = wio, fill = 'black', colour = 'black')+
  geom_tile(data=df_CMref, aes(x=x, y=y, fill=value), alpha=0.8)+
  scale_fill_viridis_c(limits = c(0.2, 1), option = "viridis")+theme_light()+
  coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE)+ylab("") + xlab('')+
  labs(fill='Habitat \nsuitability') + theme(axis.text.x=element_blank(),
                                             axis.text.y=element_blank(),
                                             axis.ticks.x=element_blank(),
                                             axis.ticks.y=element_blank(),
                                             legend.position="right",
                                           legend.justification="left",
                                           legend.margin=margin(0,0,0,0),
                                           legend.box.margin=margin(-6,-10,-6,-10))
gg_predCref

Refugee <- ggdraw()+ draw_plot(gg_predCref, 0.02, .5, .3, .45)+
  draw_plot(gg_CHref, .28, .5, .4, .45) +
  draw_plot(gg_CMref, .55, .5, .5, .45) +
  draw_plot(CHMrefDen, 0, 0, 0.95, .5) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"), c(0, 0.3, 0.6, 0), c(1, 1, 1, .5), font = "italic", size = 15)

jpeg("Refugee.jpeg", width = 2300, height =1550, res = 300)
Refugee
dev.off()

RangeRef <- ggdraw() +
  draw_plot(ggCqfoc, 0.02, 0.66, .3, .3)+
  draw_plot(ggHqfoc, 0.325, 0.66, .3, .3)+
  draw_plot(ggMqfoc, 0.6, 0.66, .3, .3)+
  draw_plot(gg_predCref, 0.02, .33, .3, .3)+
  draw_plot(gg_CHref, .28, .33, .4, .3) +
  draw_plot(gg_CMref, .56, .33, .5, .3) +
  draw_plot(CHMrefDen, 0, 0, 0.9, .3) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"), c(0, 0.32, 0.615, 0, 0.32, 0.615, 0), c(1, 1, 1, .66, .66, .66, .33), font = "italic", size = 15)

jpeg("RangeRef.jpeg", width = 2300, height =2250, res = 300)
RangeRef
dev.off()

