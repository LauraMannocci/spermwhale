
#' Partition into checkerboard (spatial folding with k=4)
#'
#' @param occ
#' @param occ.z
#' @param bg
#' @param bg.z
#' @param modelstack
#'
#' @return
#' @export
#'

checkerboard_partitioning <- function (occ, occ.z, bg, bg.z, modelstack){

  #checkerboard partitioning using spatial folding (k=4) for combined data
  #checkerboard methods subdivide geographic space equally but do not ensure a balanced number of occurrence localities in each bin
  #aggregation.factor specifies the number of grids cells to aggregate when making the underlying checkerboard pattern.

  cb <- ENMeval::get.checkerboard2(occ, modelstack, bg, aggregation.factor=c(5,5))

  #plot checkerboard partitioning

  print(ENMeval::evalplot.grps(pts = occ, pts.grp = cb$occs.grp, envs = modelstack)) #occurrences
  print(ENMeval::evalplot.grps(pts = bg, pts.grp = cb$bg.grp, envs = modelstack)) #background

  # Look at environment similarity between occurrence partitions and background partitions
  # Higher positive values indicate increasing similarity, while higher negative values indicate dissimilarity.

  #print(ENMeval::evalplot.envSim.hist(sim.type = "mess", occs.z = occ.z,
  #                                    bg.z = bg.z, occs.grp = cb$occs.grp, bg.grp = cb$bg.grp))

  # Look at environment similarity between occurrence partitions and all raster cells

  #print(ENMeval::evalplot.envSim.map(sim.type = "mess", envs = modelstack, occs.z = occ.z,
  #                                  bg.z = bg.z, occs.grp = cb$occs.grp, bg.grp = cb$bg.grp, bb.buf = 7))

  return(cb)

}








#' Evaluate niche model with user partition
#'
#' @param occ
#' @param modelstack
#' @param bg
#' @param args
#' @param user.grp
#' @param type
#'
#' @return
#' @export
#'

evaluate_model <- function(occ, modelstack, bg, args, user.grp, type) {

  # # Define a custom function that implements a performance metric not included in ENMeval.
  # # The function should have a single argument "vars", which is a list that includes the
  # # data most performance metrics should require -- the total list of these data can be found
  # # here: ?ENMevaluate. Make sure you return a data frame that specifies the names you want to
  # # see in the results tables.

  proc <- function(vars) {

    proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
    out <- data.frame(proc_auc_ratio = proc$pROC_summary[1],
                      proc_pval = proc$pROC_summary[2], row.names = NULL)

    return(out)

  }

  # evaluation with user partition

  mod <- ENMeval::ENMevaluate(occs = occ,
                       envs = modelstack,
                       bg = bg,
                       tune.args = args,
                       algorithm = 'maxnet',
                       partitions = "user",
                       user.grp = user.grp,
                       parallel = TRUE,
                       numCores = 6,
                       user.eval = proc)

  #Visualizing evaluation results

  jpeg(here::here("outputs", paste0("model_evaluation_", type, ".jpeg")), width = 1000, height = 700)
  print(ENMeval::evalplot.stats(e = mod, stats = "auc.val", color = "fc", x.var = "rm"))
  dev.off()

  return(mod)


}





#' select best model
#'
#' @param mod
#'
#' @return
#' @export
#'

select_best_model <- function(mod, type){

  # get overall results
  res <- ENMeval::eval.results(mod)

  # save metrics
  write.csv(res, here::here("outputs", paste0("metrics_ranking_", type,"_model.csv")))

  # select models maximizing auc
  opt <- res %>%
    dplyr::filter(auc.val.avg == max(auc.val.avg))

  return(opt)

}







#' make coefficients plot
#'
#' @param df_coef
#' @param ylim
#'
#' @return
#' @export
#'

make_coefficients_plot <- function(df_coef, ylim){

  p <- ggplot2::ggplot(df_coef, ggplot2::aes(x=name, y=value, color = type)) +
    ggplot2::geom_point() +
    ggforce::facet_zoom(ylim = ylim) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="right", text = ggplot2::element_text(size=16),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::xlab("Predictors") +
    ggplot2::ylab("Coefficient") +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "black", size=.5) +
    ggplot2::scale_colour_manual(values = c("#00BA38", "#F8766D"), guide = ggplot2::guide_legend(ncol = 1, direction = "horizontal",
                                                                                                 label.position="left", label.hjust = 0.5, label.vjust = 0,
                                                                                                 label.theme = ggplot2::element_text(angle = 0)))

  ggplot2::ggsave(here::here("outputs", "coefficients_plot.png"), p, width = 10, height = 7)

  return(p)

}









#' make coefficients barplot
#'
#' @param df_coef
#' @param ylim
#'
#' @return
#' @export
#'

make_coefficients_barplot <- function(df_coef, ylim){

  p <- ggplot2::ggplot(df_coef, ggplot2::aes(x=name, y=value, fill = type, width = 0.8)) +
    ggplot2::geom_bar(position="dodge", stat="identity") +
    ggforce::facet_zoom(ylim = ylim, zoom.size = 1) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size=20),
                   legend.position = "right",
                   legend.text = ggplot2::element_text(size=24),
                   axis.text.y = ggplot2::element_text(size=16),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 18)) +
    ggplot2::ylab("Coefficients") +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "black", size=.5) +
    ggplot2::scale_fill_manual(values = c("#00BA38", "#F8766D"), guide = ggplot2::guide_legend(ncol = 1, direction = "horizontal",
                                                                                               label.position="right"))

  ggplot2::ggsave(here::here("outputs", "coefficients_barplot.png"), p, width = 10, height = 7)

  return(p)

}






#' plot partial curves (modified function maxnet::plot.maxnet)
#'
#' @param x
#' @param vars
#' @param common.scale
#' @param type
#' @param ylab
#' @param plot
#' @param ...
#'
#' @return
#' @export
#'

plot_partial_curves <- function (x, vars = names(x$samplemeans), common.scale = T,
                                 type = c("link", "exponential", "cloglog", "logistic"),
                                 ylab = NULL, plot = TRUE, ...){

  type <- match.arg(type)
  nc <- ceiling(sqrt(length(vars)))
  nr <- ceiling(length(vars)/nc)

  opar = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  graphics::par(mfrow = c(nr, nc), mar = c(5, 5, 4, 2) + 0.1)

  ylim = NULL
  if (common.scale && (type == "link" || type == "exponential")) {
    vals <- do.call(c, lapply(vars, function(v) maxnet::response.plot(x,
                                                                      v, type, plot = F)$pred))
    ylim = c(min(vals), max(vals))
  }
  if (type == "cloglog" || type == "logistic")
    ylim = c(0, 1)


  if (plot) {
    for (v in vars) {
      maxnet::response.plot(x, v, type, ylim = ylim, ylab = ylab, ...)
    }
  }


  res = list()
  for (v in vars) {
    res[[v]] = maxnet::response.plot(x, v, type, ylim = ylim, plot = FALSE)
  }

  invisible(res)

}




#' plot partial curves together
#'
#'
#' @param dathis
#' @param datmod
#' @param var_name
#'
#' @return
#' @export
#'

plot_partial_curves_together <- function (dathis, datmod, var_name){

  # x limits
  if (var_name == "depth"){
    min <- min(dat_his[[var_name]], dat_mod[[var_name]])
    max <- 0
  }else{
    min <- 0
    max <- max(dat_his[[var_name]], dat_mod[[var_name]])
  }

  p <- ggplot2::ggplot(dat_his[[var_name]], ggplot2::aes(y = pred, x = get(var_name))) +
    ggplot2::geom_line(color = "#00BA38", size = 2) +
    ggplot2::geom_line(data = dat_mod[[var_name]], color = "#F8766D",  size = 2) +
    ggplot2::ylim(0,1) +
    ggplot2::ylab("Predictions") +
    ggplot2::xlab(var_name) +
    ggplot2::scale_x_continuous(limits = c(min, max)) +
    ggplot2::theme(text = ggplot2::element_text(size=20),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey90"),
                   axis.line = ggplot2::element_line(color="black"))

  ggplot2::ggsave(here::here("outputs", paste0("partial_plots_together_", var_name, ".png")), p, width = 9, height = 7)

  return(p)

}






#' Get predictor contributions
#'
#' @param modelstack
#' @param occ
#' @param bg
#' @param type
#'
#' @return
#' @export
#'

get_predictors_contributions <- function(modelstack, occ, bg, type){

  # get predictor contributions (requires to rebuild a maxent model)
  mx <- dismo::maxent(x = modelstack,
                      p = as.data.frame(occ),
                      a = bg,
                      path = here::here("outputs", "maxent"))

  contrib <- mx@results[grep("contribution", row.names(mx@results))]
  var <- names(modelStack)

  dat <- data.frame(var, contrib, type = type)

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = reorder(var, contrib), y = contrib)) +
          ggplot2::geom_bar(stat="identity") +
          ggplot2::coord_flip() +
          ggplot2::ggtitle(paste(type, "model")) +
          ggplot2::ylab("Contributions (%)") +
          ggplot2::xlab("") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(here::here("outputs", paste0("predictor_contributions_", type, ".png")), g, width = 9, height = 7)

  return(dat)

}





#' Get predictor importance
#'
#' @param modelstack
#' @param occ
#' @param bg
#' @param type
#'
#' @return
#' @export
#'

get_predictors_importance <- function(modelstack, occ, bg, type){

  # get predictor contributions (requires to rebuild a maxent model)
  mx <- dismo::maxent(x = modelstack,
                      p = as.data.frame(occ),
                      a = bg,
                      path = here::here("outputs", "maxent"))

  contrib <- mx@results[grep("importance", row.names(mx@results))]
  var <- names(modelStack)

  dat <- data.frame(var, contrib, type = type)

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = reorder(var, contrib), y = contrib)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(paste(type, "model")) +
    ggplot2::ylab("Importance (%)") +
    ggplot2::xlab("") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(here::here("outputs", paste0("predictor_importance_", type, ".png")), g, width = 9, height = 7)

  return(dat)

}




#' Make violin plot of predictors for occurrences
#'
#' @param modhis
#' @param var
#' @param ylabel
#'
#' @return
#' @export
#'

predictor_violin_occurrences <- function(modhis, var, ylabel) {

  #kruskal test
  pval = kruskal.test(get(var) ~ type, data = modhis)$p.value

  #significance level
  if(pval > 0.05) {signif <- "ns"}
  if(pval <= 0.05 & pval > 0.01) {signif <- "*"}
  if(pval <= 0.01 & pval > 0.001) {signif <- "**"}
  if(pval <= 0.001) {signif <- "***"}

  p <- ggplot2::ggplot(modhis, ggplot2::aes(x = type, y = get(var), color = type)) +
    ggplot2::geom_violin() +
    ggplot2::geom_jitter(shape = 16, size = 1, position = ggplot2::position_jitter(0.2), alpha = 0.4) +
    ggplot2::geom_boxplot(width = 0.1, fill = "white") +
    ggplot2::stat_summary(fun = mean, color = "black", geom="text", size = 5, ggplot2::aes(label=round(..y.., digits=0))) +
    ggplot2::ggtitle(signif)+
    ggplot2::xlab("") +
    ggplot2::ylab(ylabel) +
    ggplot2::theme_light() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   legend.position = "none",
                   axis.text = ggplot2::element_text(size = 20),
                   axis.title = ggplot2::element_text(size = 20)) +
    ggplot2::scale_color_manual(values = c("#00BA38", "#F8766D"))

  ggplot2::ggsave(here::here("outputs", paste0("predictors_violin_occurrences_", var, ".png")), p, width = 9, height = 7)

  return(p)

}




#' Limit rasters to predictor range of background points
#'
#' @param modelstack
#' @param bg
#' @param type
#'
#' @return
#' @export
#'

limit_rasters_to_background_points <- function(modelstack, bg, type){

  modelStack_limited <- modelstack

  depth_rg <- range(as.data.frame(bg)["depth"], na.rm = TRUE)
  slope_rg <- range(as.data.frame(bg)["slope"], na.rm = TRUE)
  # di_1000m_rg <- range(as.data.frame(bg)["di_1000m"], na.rm = TRUE)
  di_2500m_rg <- range(as.data.frame(bg)["di_2500m"], na.rm = TRUE)
  # di_5000m_rg <- range(as.data.frame(bg)["di_5000m"], na.rm = TRUE)
  di_coast_rg <- range(as.data.frame(bg)["di_coast"], na.rm = TRUE)
  di_seaM_rg <- range(as.data.frame(bg)["di_seaM"], na.rm = TRUE)
  di_guy_rg <- range(as.data.frame(bg)["di_guy"], na.rm = TRUE)
  di_Trough_rg <- range(as.data.frame(bg)["di_Trough"], na.rm = TRUE)
  di_Tren_rg <- range(as.data.frame(bg)["di_Tren"], na.rm = TRUE)
  di_spRid_rg <- range(as.data.frame(bg)["di_spRid"], na.rm = TRUE)

  modelStack_limited[["depth"]] <- raster::clamp(modelStack_limited[["depth"]], lower = depth_rg[1], upper = depth_rg[2], useValues = FALSE)
  modelStack_limited[["slope"]] <- raster::clamp(modelStack_limited[["slope"]], lower = slope_rg[1], upper = slope_rg[2], useValues = FALSE)
  # modelStack_limited[["di_1000m"]] <- raster::clamp(modelStack_limited[["di_1000m"]], lower = di_1000m_rg[1], upper = di_1000m_rg[2], useValues = FALSE)
  modelStack_limited[["di_2500m"]] <- raster::clamp(modelStack_limited[["di_2500m"]], lower = di_2500m_rg[1], upper = di_2500m_rg[2], useValues = FALSE)
  # modelStack_limited[["di_5000m"]] <- raster::clamp(modelStack_limited[["di_5000m"]], lower = di_5000m_rg[1], upper = di_5000m_rg[2], useValues = FALSE)
  modelStack_limited[["di_coast"]] <- raster::clamp(modelStack_limited[["di_coast"]], lower = di_coast_rg[1], upper = di_coast_rg[2], useValues = FALSE)
  modelStack_limited[["di_seaM"]] <- raster::clamp(modelStack_limited[["di_seaM"]], lower = di_seaM_rg[1], upper = di_seaM_rg[2], useValues = FALSE)
  modelStack_limited[["di_guy"]] <- raster::clamp(modelStack_limited[["di_guy"]], lower = di_guy_rg[1], upper = di_guy_rg[2], useValues = FALSE)
  modelStack_limited[["di_Trough"]] <- raster::clamp(modelStack_limited[["di_Trough"]], lower = di_Trough_rg[1], upper = di_Trough_rg[2], useValues = FALSE)
  modelStack_limited[["di_Tren"]] <- raster::clamp(modelStack_limited[["di_Tren"]], lower = di_Tren_rg[1], upper = di_Tren_rg[2], useValues = FALSE)
  modelStack_limited[["di_spRid"]] <- raster::clamp(modelStack_limited[["di_spRid"]], lower = di_spRid_rg[1], upper = di_spRid_rg[2], useValues = FALSE)

  # map rasters limited to background of modern model
  jpeg(here::here("outputs", paste0("rasters_limited_background_", type, ".jpeg")), width = 960, height = 960)
  par(mar=c(0,0,1,0))
  raster::plot(modelStack_limited)
  dev.off()

  return(modelStack_limited)

}






#' Helper function used by get_extra_extent to facilitate computations of convex hulls and Gower's distances
#'
#' @param calibration_data
#' @param test_data
#' @param var_name
#' @param choice
#'
#' @return
#' @export
#'

make_cfact <- function (calibration_data, test_data, var_name, choice) {


  ## helper function to rescale predictor
  rescale2 <- function (ynew, y) {

    return ((ynew - mean (y, na.rm = TRUE)) / (sd (y, na.rm = TRUE)))

  }

  ## standardize new data to predict from
  # this simplifies computation A LOT!
  make_X <- function (calibration_data, test_data, var_name){
    X <- sapply(var_name,
                function (k) { rescale2 (ynew = test_data[, k],
                                         y = calibration_data[, k]
                )}
    )
    X <- as.data.frame (X)
    names (X) <- var_name
    return (X)
  }
  # compute counterfactuals
  cfact <- WhatIf::whatif (formula = NULL,
                           data = make_X (calibration_data = calibration_data, test_data = calibration_data, var_name),
                           cfact = make_X (calibration_data = calibration_data, test_data = test_data, var_name),
                           mc.cores = parallel::detectCores()/1.3,
                           choice = choice)

  # function to get the % of interpolations
  prop_cfact <- function (cfact) {
    return (round (100 * length (which (cfact$in.hull == TRUE)) / length (cfact$in.hull), 2))
  }

  return (list (cfact = cfact, interpolation = prop_cfact (cfact = cfact)))

}








#' get multidimensional extrapolation extent using background points as reference set and prediction region as test set
#'
#' @param variables
#' @param bg.z
#' @param type
#'
#' @return
#' @export
#'

get_extra_extent <- function(variables, modelstack, bg.z, type) {

  #define calibration datasets
  cal <- bg.z[, variables]

  #define test datasets from the raster stack
  test <- do.call(cbind, lapply(1:length(variables), function(i) raster::rasterToPoints(modelstack[[variables[i]]])))
  test <- as.data.frame(test)
  test <- test[, c("x", "y", variables)]

  #calculate extent based on gower distance
  cfact1 <- make_cfact(calibration_data = na.omit(cal), test_data = na.omit(test), var_name = variables, choice = "both")

  #assign results of cfact to test data (interpolation 1, extrapolation -1)
  test_na <- na.omit(test)
  test$cfact1 <- NA
  test$cfact1[as.numeric (row.names(test_na))] <- ifelse (cfact1$cfact$in.hull, 1, -1)

  #output map of interpolation vs extrapolation
  g <- ggplot2::ggplot() +
    ggplot2::geom_tile (data = test, ggplot2::aes(x = x, y = y, fill = factor(cfact1))) +
    ggplot2::xlab("Easting") +
    ggplot2::ylab("Northing") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "grey50"), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=10),
          axis.title = ggplot2::element_text(size=6),
          axis.text = ggplot2::element_text(size=6)) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_manual (name="", values= c("midnightblue","yellow1"), labels=c("Extra","Inter"), breaks=c("-1","1"), limits=c("-1","1")) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8),
          legend.title = ggplot2::element_text(size=8))

  ggplot2::ggsave(here::here("outputs", paste0("map_inter_vs_extra_", type, ".png")), g, width = 9, height = 7)

  return(test)

}








#' Convert predictions to dataframe
#'
#' @param pred
#' @param type
#'
#' @return
#' @export
#'

convert_predictions_to_df <- function(pred, type){

  df_pred <- as(pred, "SpatialPixelsDataFrame")
  df_pred <- as.data.frame(df_pred)
  colnames(df_pred) <- c("value", "x", "y")
  df_pred$type <- factor(c(type))

  return(df_pred)

}





#' Plot predictions
#'
#' @param wio
#' @param df_pred
#' @param type
#'
#' @return
#' @export
#'

plot_predictions <- function(wio, df_pred, type){

  g <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = wio) +
    ggplot2::geom_tile(data = df_pred, ggplot2::aes(x = x, y = y, fill = value), alpha=0.8) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab(paste(type, "model")) +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))

  ggplot2::ggsave(here::here("outputs", paste0("map_predictions_", type, ".png")), g, width = 9, height = 7)

  return(g)

}




#' Plot predictions residuals
#'
#' @param wio
#' @param df_pred
#' @param type
#'
#' @return
#' @export
#'

plot_predictions_residuals <- function(wio, pred_his, pred_mod){

  #get residuals of linear model
  model <- lm(raster::getValues(pred_his) ~ raster::getValues(pred_mod), na.action=na.exclude)
  res <- pred_his
  res[] <- residuals.lm(model)

  #get r squared
  print("r-squared :")
  print(summary(model)$r.squared)

  #convert to dataframe for plotting
  df_res <- as(res, "SpatialPixelsDataFrame")
  df_res <- as.data.frame(df_res)
  colnames(df_res) <- c("value", "x", "y")

  g <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = wio) +
    ggplot2::geom_tile(data = df_res, ggplot2::aes(x = x, y = y, fill = value), alpha=0.8) +
    ggplot2::scale_fill_viridis_c(limits = c(-1.2, 0.6), option = "magma")+
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("Residuals") +
    ggplot2::labs(fill = 'Residuals')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 17),
                   legend.title = ggplot2::element_text(size = 17),
                   axis.title = ggplot2::element_text(size = 18),
                   axis.text = ggplot2::element_text(size = 12),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))

  ggplot2::ggsave(here::here("outputs", "map_predictions_residuals.png"), g, width = 9, height = 7)

  return(g)

}





#' Plot predictions with extrapolation extent
#'
#' @param wio
#' @param df_pred
#' @param type
#'
#' @return
#' @export
#'

plot_predictions_with_extra <- function(wio, df_pred, type, df_test){

  #select extrapolation data points
  df_test <- subset(df_test, cfact1 == -1)

  g <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_test, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    #add countries
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab(paste(type, "model")) +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA))

  ggplot2::ggsave(here::here("outputs", paste0("map_predictions_with_extra_", type, ".png")), g, width = 9, height = 7)

  return(g)

}




#' Plot predictions with extrapolation extent and mpas
#'
#' @param wio
#' @param df_pred
#' @param type
#'
#' @return
#' @export
#'

plot_predictions_with_extra_mpas <- function(wio, df_pred, type, df_test, mpa){

  #select extrapolation data points
  df_test <- subset(df_test, cfact1 == -1)

  g <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_test, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    #add mpas
    ggplot2::geom_sf(data = mpa, color = "white", fill = NA, size = 0.6) +
    #add countries
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab(paste(type, "model")) +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 17),
                   legend.title = ggplot2::element_text(size = 17),
                   axis.title = ggplot2::element_text(size = 18),
                   axis.text = ggplot2::element_text(size = 10),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA))

  ggplot2::ggsave(here::here("outputs", paste0("map_predictions_with_extra_", type, "_mpas.png")), g, width = 9, height = 7)

  return(g)

}




#' Plot predictions residuals with extrapolation extent
#'
#' @param pred_his
#' @param pred_mod
#' @param df_his
#' @param df_mod
#' @param wio
#'
#' @return
#' @export
#'

plot_predictions_residuals_with_extra <- function(wio, pred_his, pred_mod, df_his, df_mod){

  #get residuals of linear model
  model <- lm(raster::getValues(pred_his) ~ raster::getValues(pred_mod), na.action=na.exclude)
  res <- pred_his
  res[] <- residuals.lm(model)

  #get r squared
  print("r-squared :")
  print(summary(model)$r.squared)

  #select extrapolation data points
  df_his <- subset(df_his, cfact1 == -1)
  df_mod <- subset(df_mod, cfact1 == -1)

  #convert to dataframe for plotting
  df_res <- as(res, "SpatialPixelsDataFrame")
  df_res <- as.data.frame(df_res)
  colnames(df_res) <- c("value", "x", "y")

  g <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_res, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_his, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::geom_tile(data = df_mod, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::scale_fill_viridis_c(limits = c(-1.2, 0.6), option = "magma")+
    #add countries
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("Residuals") +
    ggplot2::labs(fill = 'Residuals')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 17),
                   legend.title = ggplot2::element_text(size = 17),
                   axis.title = ggplot2::element_text(size = 18),
                   axis.text = ggplot2::element_text(size = 10),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA))

  ggplot2::ggsave(here::here("outputs", "map_predictions_residuals_with_extra.png"), g, width = 9, height = 7)

  return(g)

}



#' compare predictions (removing extrapolation zones)
#'
#' @param pred_his
#' @param df_extra_his
#' @param pred_mod
#' @param df_extra_mod
#'
#' @return
#' @export
#'

compare_predictions_extra <- function(pred_his, df_extra_his, pred_mod, df_extra_mod){

  ### historical

  # convert extrapolation dataframe to raster
  extra_his <- raster::rasterFromXYZ(df_extra_his[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")

  # mask predictions raster with extrapolation (-1 value) raster and convert to dataframe
  extra_his[extra_his == -1] <- NA
  pred_his <- raster::mask(pred_his, extra_his)
  df_pred_his <- as(pred_his, "SpatialPixelsDataFrame")
  df_pred_his <- as.data.frame(df_pred_his)
  colnames(df_pred_his) <- c("value", "x", "y")


  ### modern

  # convert extrapolation dataframe to raster
  extra_mod <- raster::rasterFromXYZ(df_extra_mod[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")

  # mask predictions raster with extrapolation (-1 value) raster and convert to dataframe
  extra_mod[extra_mod == -1] <- NA
  pred_mod <- raster::mask(pred_mod, extra_mod)
  df_pred_mod <- as(pred_mod, "SpatialPixelsDataFrame")
  df_pred_mod <- as.data.frame(df_pred_mod)
  colnames(df_pred_mod) <- c("value", "x", "y")


  ### bind dataframes
  df_pred_his$type <- factor(c('his'))
  df_pred_mod$type <- factor(c('mod'))
  df <- rbind(df_pred_his, df_pred_mod)


  ### calculating medians
  print("median historical predictions:")
  print(median(df_pred_his$value))
  print("median modern predictions:")
  print(median(df_pred_mod$value))


  ### kruskall tests of predictions in mpa vs all predictions
  # non-parametric method for testing whether samples originate from the same distribution
  print(kruskal.test(value ~ type, data = df))



  ### pearsons's correlation coefficients
  predStack <- raster::stack(c(pred_his, pred_mod))
  names(predStack) <- c("historical", "modern")
  jnk <- raster::layerStats(predStack, 'pearson', na.rm=T)
  corr_matrix <- jnk$'pearson correlation coefficient'
  print("pearsons's correlation coefficients")
  print(corr_matrix)




}




#' test for "residual" mpa effect with kruskall test
#'
#' @param pred
#' @param df_pred
#' @param mpa_sf
#' @param name
#'
#' @return
#' @export
#'

test_residual_mpa_effect <- function(pred, df_pred, mpa_sf, name){

  # mask predictions with mpa and convert to dataframe
  mpamod <- raster::mask(pred, mpa_sf)
  df_mpamod <- as(mpamod, "SpatialPixelsDataFrame")
  df_mpamod <- as.data.frame(df_mpamod)
  colnames(df_mpamod) <- c("value", "x", "y")

  ### kruskall tests of predictions in mpa vs all predictions
  # non-parametric method for testing whether samples originate from the same distribution
  df_mpamod$type <- factor(c('mpa'))
  df_pred$type <- factor(c('region'))
  df_mpamod_pred <- rbind(df_pred, df_mpamod)
  df_mpamod_pred$model <- factor(name)

  print(kruskal.test(value ~ type, data = df_mpamod_pred))

  return(df_mpamod_pred)

}




#' test for "residual" mpa effect with kruskall test (removing extrapolation zones)
#'
#' @param pred
#' @param mpa_sf
#' @param name
#' @param df_extra
#'
#' @return
#' @export
#'

test_residual_mpa_effect_extra <- function(pred, df_extra, mpa_sf, name){

  # convert extrapolation dataframe to raster
  extra <- raster::rasterFromXYZ(df_extra[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")

  # mask predictions raster with extrapolation (-1 value) raster and convert to dataframe
  extra[extra == -1] <- NA
  pred <- raster::mask(pred, extra)
  df_pred <- as(pred, "SpatialPixelsDataFrame")
  df_pred <- as.data.frame(df_pred)
  colnames(df_pred) <- c("value", "x", "y")

  # mask predictions with mpa and convert to dataframe
  pred_mpa <- raster::mask(pred, mpa_sf)
  df_pred_mpa <- as(pred_mpa, "SpatialPixelsDataFrame")
  df_pred_mpa <- as.data.frame(df_pred_mpa)
  colnames(df_pred_mpa) <- c("value", "x", "y")

  #bind dataframes
  df_pred_mpa$type <- factor(c('mpa'))
  df_pred$type <- factor(c('region'))
  df <- rbind(df_pred, df_pred_mpa)
  df$model <- factor(name)

  ### kruskall tests of predictions in mpa vs all predictions
  # non-parametric method for testing whether samples originate from the same distribution
  print(kruskal.test(value ~ type, data = df))

  return(df)
}





#' Calculate median predictions in and out of mpas (removing extrapolation zones)
#'
#' @param pred
#' @param df_extra
#' @param mpa_sf
#'
#' @return
#' @export
#'

calculate_median_predictions_in_out_mpa_extra <- function(pred, df_extra, mpa_sf){

  # convert extrapolation dataframe to raster
  extra <- raster::rasterFromXYZ(df_extra[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")

  # mask predictions raster with extrapolation (-1 value) raster and convert to dataframe
  extra[extra == -1] <- NA
  pred <- raster::mask(pred, extra)
  df_pred <- as(pred, "SpatialPixelsDataFrame")
  df_pred <- as.data.frame(df_pred)
  colnames(df_pred) <- c("value", "x", "y")

  # mask predictions with mpa and convert to dataframe (predictions in mpa)
  pred_inmpa <- raster::mask(pred, mpa_sf)
  df_pred_inmpa <- as(pred_inmpa, "SpatialPixelsDataFrame")
  df_pred_inmpa <- as.data.frame(df_pred_inmpa)
  colnames(df_pred_inmpa) <- c("value", "x", "y")

  # inverse mask predictions with mpa and convert to dataframe (predictions outside mpa)
  pred_outmpa <- raster::mask(pred, mpa_sf, inverse = TRUE)
  df_pred_outmpa <- as(pred_outmpa, "SpatialPixelsDataFrame")
  df_pred_outmpa <- as.data.frame(df_pred_outmpa)
  colnames(df_pred_outmpa) <- c("value", "x", "y")

  #calculating medians in and out mpa
  print("median in mpa :")
  print(median(df_pred_inmpa$value))
  print("median oustide mpa :")
  print(median(df_pred_outmpa$value))

}



#' Calculate median predictions in mpas vs in all region (removing extrapolation zones)
#'
#' @param pred
#' @param df_extra
#' @param mpa_sf
#'
#' @return
#' @export
#'

calculate_median_predictions_in_mpa_all_region_extra <- function(pred, df_extra, mpa_sf){

  # convert extrapolation dataframe to raster
  extra <- raster::rasterFromXYZ(df_extra[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")

  # mask predictions raster with extrapolation (-1 value) raster and convert to dataframe
  extra[extra == -1] <- NA
  pred <- raster::mask(pred, extra)
  df_pred <- as(pred, "SpatialPixelsDataFrame")
  df_pred <- as.data.frame(df_pred)
  colnames(df_pred) <- c("value", "x", "y")

  # mask predictions with mpa and convert to dataframe (predictions in mpa)
  pred_inmpa <- raster::mask(pred, mpa_sf)
  df_pred_inmpa <- as(pred_inmpa, "SpatialPixelsDataFrame")
  df_pred_inmpa <- as.data.frame(df_pred_inmpa)
  colnames(df_pred_inmpa) <- c("value", "x", "y")

  #calculating medians in and out mpa
  print("median in mpa :")
  print(median(df_pred_inmpa$value))
  print("median all region :")
  print(median(df_pred$value))

}




#' Clip and plot predictions in eez and add violin plot to compare them
#'
#' @param eez_sh ...
#' @param eez_name ...
#' @param pred_his ...
#' @param pred_mod ...
#' @param wio ...
#' @param obs_mod
#' @param obs_his
#'
#' @return
#' @export
#'

plot_predictions_in_eez <- function(eez_sh, eez_name, pred_his, pred_mod, obs_mod, obs_his, wio){

  if (eez_name == "Seychelles"){names_from_sh <- "Seychellois Exclusive Economic Zone"}
  if (eez_name == "Madagascar"){names_from_sh <- "Madagascan Exclusive Economic Zone"  }
  if (eez_name == "Mauritius"){names_from_sh <- "Mauritian Exclusive Economic Zone"}
  if (eez_name == "Chagos"){names_from_sh <- "Chagos Archipelago Exclusive Economic Zone"}
  if (eez_name == "Maldives"){names_from_sh <- "Maldives Exclusive Economic Zone"}
  if (eez_name == "Sri Lanka"){names_from_sh <- "Sri Lankan Exclusive Economic Zone"}
  if (eez_name == "Reunion"){names_from_sh <- "Réunion Exclusive Economic Zone"}
  if (eez_name == "Comoros"){names_from_sh <- "Comoran Exclusive Economic Zone"}
  if (eez_name == "Oman"){names_from_sh <- "Omani Exclusive Economic Zone"  }
  if (eez_name == "Yemen"){names_from_sh <- "Yemeni Exclusive Economic Zone" }
  if (eez_name == "Somali"){names_from_sh <-"Somali Exclusive Economic Zone" }



  #select eez polygon
  eez_sh %>%
    dplyr::filter(GEONAME == names_from_sh) -> eez

  #mask rasters with eez polygon
  pred_his <- raster::mask(pred_his, eez)
  pred_mod <- raster::mask(pred_mod, eez)

  ###Historical

  #convert to dataframe
  type <- "Historical"
  df_pred_his <- as(pred_his, "SpatialPixelsDataFrame")
  df_pred_his <- as.data.frame(df_pred_his)
  colnames(df_pred_his) <- c("value", "x", "y")
  df_pred_his$type <- factor(c(type))

  #make map
  gHis <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_his, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_point(data = obs_his, ggplot2::aes(x = Lon, y = Lat, alpha = 0.6)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(raster::extent(eez)[1], raster::extent(eez)[2]), ylim = c(raster::extent(eez)[3], raster::extent(eez)[4]), expand = FALSE) +
    ggspatial::fixed_plot_aspect(ratio = 1) +
    ggplot2::ylab("") +
    ggplot2::xlab("Historical") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))


  ###Modern

  #convert to dataframe
  type <- "Modern"
  df_pred_mod <- as(pred_mod, "SpatialPixelsDataFrame")
  df_pred_mod <- as.data.frame(df_pred_mod)
  colnames(df_pred_mod) <- c("value", "x", "y")
  df_pred_mod$type <- factor(c(type))


  #make map
  gMod <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_mod, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_point(data = obs_mod, ggplot2::aes(x = Lon, y = Lat, alpha = 0.6)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(raster::extent(eez)[1], raster::extent(eez)[2]), ylim = c(raster::extent(eez)[3], raster::extent(eez)[4]), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("Modern") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))



  #violin plot
  pred <- rbind(df_pred_his, df_pred_mod)
  violin <- ggplot2::ggplot(pred, ggplot2::aes(x = type, y = value, color = type)) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1, fill = "white") +
    ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("#00BA38", "#F8766D"))


  if (eez_name == "Yemen"){
    g <- cowplot::ggdraw() +
      cowplot::draw_plot(gHis, x = 0, y = 0.1, height = .7, width = .31) +
      cowplot::draw_plot(gMod, x = .35, y = 0.1, height = .7, width = .31) +
      cowplot::draw_plot(violin, x = .69, y = 0.1, height = .7, width = .3) +
      cowplot::draw_label(eez_name, x = 0.5,  y = 0.92,  size = 26)
  }else{
    g <- cowplot::ggdraw() +
      cowplot::draw_plot(gHis, x = -0.18, y = 0.1, height = .7, width = .7) +
      cowplot::draw_plot(gMod, x = .16, y = 0.1, height = .7, width = .7) +
      cowplot::draw_plot(violin, x = .69, y = 0.1, height = .7, width = .3) +
      cowplot::draw_label(eez_name, x = 0.5,  y = 0.92,  size = 26)
  }

  ggplot2::ggsave(here::here("outputs", paste0("plot_predictions_in_", eez_name, ".jpeg")), g, width = 15.5, height = 5, bg = "white")


}





#' Clip and plot predictions in eez and add violin plot to compare them (removing extrapolation zones)
#' @param eez_sh ...
#' @param eez_name ...
#' @param pred_his ...
#' @param pred_mod ...
#' @param wio ...
#' @param df_extra_his
#' @param df_extra_mod
#'
#' @return
#' @export
#'

plot_predictions_in_eez_extra <- function(eez_sh, eez_name, pred_his, pred_mod, df_extra_his, df_extra_mod, wio){

  if (eez_name == "Seychelles"){names_from_sh <- "Seychellois Exclusive Economic Zone"}
  if (eez_name == "Madagascar"){names_from_sh <- "Madagascan Exclusive Economic Zone"  }
  if (eez_name == "Mauritius"){names_from_sh <- "Mauritian Exclusive Economic Zone"}
  if (eez_name == "Chagos"){names_from_sh <- "Chagos Archipelago Exclusive Economic Zone"}
  if (eez_name == "Maldives"){names_from_sh <- "Maldives Exclusive Economic Zone"}
  if (eez_name == "Sri Lanka"){names_from_sh <- "Sri Lankan Exclusive Economic Zone"}
  if (eez_name == "Reunion"){names_from_sh <- "Réunion Exclusive Economic Zone"}
  if (eez_name == "Comoros"){names_from_sh <- "Comoran Exclusive Economic Zone"}
  if (eez_name == "Oman"){names_from_sh <- "Omani Exclusive Economic Zone"  }
  if (eez_name == "Yemen"){names_from_sh <- "Yemeni Exclusive Economic Zone" }
  if (eez_name == "Somali"){names_from_sh <-"Somali Exclusive Economic Zone" }


  #select eez polygon
  eez_sh %>%
    dplyr::filter(GEONAME == names_from_sh) -> eez


  #mask prediction rasters with eez polygon
  pred_his <- raster::mask(pred_his, eez)
  pred_mod <- raster::mask(pred_mod, eez)


  #mask extrapolation rasters with eez polygon
  extra_his <- raster::rasterFromXYZ(df_extra_his[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")
  extra_mod <- raster::rasterFromXYZ(df_extra_mod[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")
  extra_his <- raster::mask(extra_his, eez)
  extra_mod <- raster::mask(extra_mod, eez)


  #convert extrapolation rasters to dataframes
  type <- "Historical"
  df_extra_his <- as(extra_his, "SpatialPixelsDataFrame")
  df_extra_his <- as.data.frame(df_extra_his)
  colnames(df_extra_his) <- c("value", "x", "y")
  df_extra_his$type <- factor(c(type))

  type <- "Modern"
  df_extra_mod <- as(extra_mod, "SpatialPixelsDataFrame")
  df_extra_mod <- as.data.frame(df_extra_mod)
  colnames(df_extra_mod) <- c("value", "x", "y")
  df_extra_mod$type <- factor(c(type))


  #select extrapolation data points (-1)
  df_extra_his <- subset(df_extra_his, value == -1)
  df_extra_mod <- subset(df_extra_mod, value == -1)


  #convert predictions to dataframes
  type <- "Historical"
  df_pred_his <- as(pred_his, "SpatialPixelsDataFrame")
  df_pred_his <- as.data.frame(df_pred_his)
  colnames(df_pred_his) <- c("value", "x", "y")
  df_pred_his$type <- factor(c(type))

  type <- "Modern"
  df_pred_mod <- as(pred_mod, "SpatialPixelsDataFrame")
  df_pred_mod <- as.data.frame(df_pred_mod)
  colnames(df_pred_mod) <- c("value", "x", "y")
  df_pred_mod$type <- factor(c(type))



  ######Historical map

  gHis <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_his, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_his, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(raster::extent(eez)[1], raster::extent(eez)[2]), ylim = c(raster::extent(eez)[3], raster::extent(eez)[4]), expand = FALSE) +
    ggspatial::fixed_plot_aspect(ratio = 1) +
    ggplot2::ylab("") +
    ggplot2::xlab("Historical") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))


  ######Modern map

  gMod <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_mod, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_mod, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(raster::extent(eez)[1], raster::extent(eez)[2]), ylim = c(raster::extent(eez)[3], raster::extent(eez)[4]), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("Modern") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))



  ######violin plot

  #union of extrapolation rasters (first need to extend to extent of predictions)
  newextent = raster::extent(pred_mod)
  extra_his_extend = raster::extend(extra_his, newextent)
  extra_mod_extend = raster::extend(extra_mod, newextent)
  s <- raster::stack(extra_his_extend, extra_mod_extend)
  extra <- raster::calc(s, base::sum, na.rm=TRUE)
  extra[extra == -2] <- NA #-2 means there is extrapolation for the historical and modern rasters
  extra[extra == 0] <- NA #0 means there is extrapolation for one of the two rasters
  #this leaves 2 pixels where there is no extrapolation


  #mask predictions to extrapolation raster (set NA to prediction pixels)
  pred_his <- raster::mask(pred_his, extra)
  pred_mod <- raster::mask(pred_mod, extra)


  #convert predictions to dataframes
  type <- "Historical"
  df_pred_his <- as(pred_his, "SpatialPixelsDataFrame")
  df_pred_his <- as.data.frame(df_pred_his)
  colnames(df_pred_his) <- c("value", "x", "y")
  df_pred_his$type <- factor(c(type))

  type <- "Modern"
  df_pred_mod <- as(pred_mod, "SpatialPixelsDataFrame")
  df_pred_mod <- as.data.frame(df_pred_mod)
  colnames(df_pred_mod) <- c("value", "x", "y")
  df_pred_mod$type <- factor(c(type))


  #bind
  pred <- rbind(df_pred_his, df_pred_mod)


  violin <- ggplot2::ggplot(pred, ggplot2::aes(x = type, y = value, color = type)) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1, fill = "white") +
    ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("#00BA38", "#F8766D"))


  if (eez_name == "Yemen"){
    g <- cowplot::ggdraw() +
      cowplot::draw_plot(gHis, x = 0, y = 0.1, height = .7, width = .31) +
      cowplot::draw_plot(gMod, x = .35, y = 0.1, height = .7, width = .31) +
      cowplot::draw_plot(violin, x = .69, y = 0.1, height = .7, width = .3) +
      cowplot::draw_label(eez_name, x = 0.5,  y = 0.92,  size = 26)
  }else{
    g <- cowplot::ggdraw() +
      cowplot::draw_plot(gHis, x = -0.18, y = 0.1, height = .7, width = .7) +
      cowplot::draw_plot(gMod, x = .16, y = 0.1, height = .7, width = .7) +
      cowplot::draw_plot(violin, x = .69, y = 0.1, height = .7, width = .3) +
      cowplot::draw_label(eez_name, x = 0.5,  y = 0.92,  size = 26)
  }

  ggplot2::ggsave(here::here("outputs", paste0("plot_predictions_in_", eez_name, "_extra.jpeg")), g, width = 15.5, height = 5, bg = "white")


}






#' Clip and plot predictions in high seas and add violin plot to compare them
#'
#' @param eez_sh
#' @param pred_his
#' @param pred_mod
#' @param wio
#'
#' @return
#' @export
#'

plot_predictions_in_high_seas <- function(eez_sh, pred_his, pred_mod, wio){

  ###Historical

  #inverse mask predictions with eez to extract high seas
  pred_his2 <- raster::mask(pred_his, eez, inverse=T)

  #remove remaining pixels along coastlines
  wio2 <- sf::as_Spatial(wio)
  wio3 <- raster::buffer(wio2, width = 0.5)
  wio4 <- raster::crop(wio3, raster::extent(pred_his))
  pred_his3 <- raster::mask(pred_his2, wio4, inverse=T) #these are the predictions extracted in high seas

  #convert to dataframe
  type <- "Historical"
  df_pred_his <- as(pred_his3, "SpatialPixelsDataFrame")
  df_pred_his <- as.data.frame(df_pred_his)
  colnames(df_pred_his) <- c("value", "x", "y")
  df_pred_his$type <- factor(c(type))

  #make map
  gHis <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_his, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(raster::extent(pred_his)[1], raster::extent(pred_his)[2]), ylim = c(raster::extent(pred_his)[3], raster::extent(pred_his)[4]), expand = FALSE) +
    ggspatial::fixed_plot_aspect(ratio = 1) +
    ggplot2::ylab("") +
    ggplot2::xlab("Historical") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))




  ###Modern

  #inverse mask predictions with eez to extract high seas
  pred_mod2 <- raster::mask(pred_mod, eez, inverse=T)

  #remove remaining pixels along coastlines
  pred_mod3 <- raster::mask(pred_mod2, wio4, inverse=T) #these are the predictions extracted in high seas

  #convert to dataframe
  type <- "Modern"
  df_pred_mod <- as(pred_mod3, "SpatialPixelsDataFrame")
  df_pred_mod <- as.data.frame(df_pred_mod)
  colnames(df_pred_mod) <- c("value", "x", "y")
  df_pred_mod$type <- factor(c(type))


  #make map
  gMod <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_mod, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(raster::extent(pred_mod)[1], raster::extent(pred_mod)[2]), ylim = c(raster::extent(pred_mod)[3], raster::extent(pred_mod)[4]), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("Modern") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))



  #violin plot
  pred <- rbind(df_pred_his, df_pred_mod)

  violin <- ggplot2::ggplot(pred, ggplot2::aes(x = type, y = value, color = type)) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1, fill = "white") +
    ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("#00BA38", "#F8766D"))


  g <- cowplot::ggdraw() +
    cowplot::draw_plot(gHis, x = -0.18, y = 0.1, height = .7, width = .7) +
    cowplot::draw_plot(gMod, x = .16, y = 0.1, height = .7, width = .7) +
    cowplot::draw_plot(violin, x = .69, y = 0.1, height = .7, width = .3) +
    cowplot::draw_label("High seas", x = 0.5,  y = 0.92,  size = 26)

  ggplot2::ggsave(here::here("outputs", "plot_predictions_in_high_seas.jpeg"), g, width = 15.5, height = 5, bg = "white")

}





#' Clip and plot predictions in high seas and add violin plot to compare them (removing extrapolation zones)
#'
#' @param eez_sh
#' @param pred_his
#' @param pred_mod
#' @param wio
#'
#' @return
#' @export
#'

plot_predictions_in_high_seas_extra <- function(eez_sh, pred_his, pred_mod, df_extra_his, df_extra_mod, wio){


  #convert extrapolation df to rasters
  extra_his <- raster::rasterFromXYZ(df_extra_his[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")
  extra_mod <- raster::rasterFromXYZ(df_extra_mod[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")


  #inverse mask predictions with eez to extract high seas
  pred_his2 <- raster::mask(pred_his, eez, inverse=T)
  pred_mod2 <- raster::mask(pred_mod, eez, inverse=T)


  #inverse mask extrapolations with eez to extract high seas
  extra_his2 <- raster::mask(extra_his, eez, inverse=T)
  extra_mod2 <- raster::mask(extra_mod, eez, inverse=T)


  #remove remaining pixels along coastlines for predictions
  wio2 <- sf::as_Spatial(wio)
  wio3 <- raster::buffer(wio2, width = 0.5)
  wio4 <- raster::crop(wio3, raster::extent(pred_his))
  pred_his3 <- raster::mask(pred_his2, wio4, inverse=T) #these are the predictions extracted in high seas
  pred_mod3 <- raster::mask(pred_mod2, wio4, inverse=T) #these are the predictions extracted in high seas


  #remove remaining pixels along coastlines for extrapolations
  extra_his3 <- raster::mask(extra_his2, wio4, inverse=T) #these are the extrapolations extracted in high seas
  extra_mod3 <- raster::mask(extra_mod2, wio4, inverse=T) #these are the extrapolations extracted in high seas



  #convert prediction rasters to dataframes
  type <- "Historical"
  df_pred_his <- as(pred_his3, "SpatialPixelsDataFrame")
  df_pred_his <- as.data.frame(df_pred_his)
  colnames(df_pred_his) <- c("value", "x", "y")
  df_pred_his$type <- factor(c(type))

  type <- "Modern"
  df_pred_mod <- as(pred_mod3, "SpatialPixelsDataFrame")
  df_pred_mod <- as.data.frame(df_pred_mod)
  colnames(df_pred_mod) <- c("value", "x", "y")
  df_pred_mod$type <- factor(c(type))



  #convert extrapolation rasters to dataframes
  type <- "Historical"
  df_extra_his <- as(extra_his3, "SpatialPixelsDataFrame")
  df_extra_his <- as.data.frame(df_extra_his)
  colnames(df_extra_his) <- c("value", "x", "y")
  df_extra_his$type <- factor(c(type))

  type <- "Modern"
  df_extra_mod <- as(extra_mod3, "SpatialPixelsDataFrame")
  df_extra_mod <- as.data.frame(df_extra_mod)
  colnames(df_extra_mod) <- c("value", "x", "y")
  df_extra_mod$type <- factor(c(type))



  #select extrapolation data points (-1)
  df_extra_his <- subset(df_extra_his, value == -1)
  df_extra_mod <- subset(df_extra_mod, value == -1)




  #make map modern
  gHis <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_his, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_his, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    ggplot2::coord_sf(xlim = c(raster::extent(pred_his)[1], raster::extent(pred_his)[2]), ylim = c(raster::extent(pred_his)[3], raster::extent(pred_his)[4]), expand = FALSE) +
    ggspatial::fixed_plot_aspect(ratio = 1) +
    ggplot2::ylab("") +
    ggplot2::xlab("Historical") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))



  #make map historical
  gMod <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_mod, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_mod, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ggplot2::coord_sf(xlim = c(raster::extent(pred_mod)[1], raster::extent(pred_mod)[2]), ylim = c(raster::extent(pred_mod)[3], raster::extent(pred_mod)[4]), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("Modern") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))




  #violin plot
  pred <- rbind(df_pred_his, df_pred_mod)

  violin <- ggplot2::ggplot(pred, ggplot2::aes(x = type, y = value, color = type)) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1, fill = "white") +
    ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("#00BA38", "#F8766D"))


  g <- cowplot::ggdraw() +
    cowplot::draw_plot(gHis, x = -0.18, y = 0.1, height = .7, width = .7) +
    cowplot::draw_plot(gMod, x = .16, y = 0.1, height = .7, width = .7) +
    cowplot::draw_plot(violin, x = .69, y = 0.1, height = .7, width = .3) +
    cowplot::draw_label("High seas", x = 0.5,  y = 0.92,  size = 26)

  ggplot2::ggsave(here::here("outputs", "plot_predictions_in_high_seas_extra.jpeg"), g, width = 15.5, height = 5, bg = "white")

}





#' Clip and plot predictions in all eezs and add violin plot to compare them (removing extrapolation zones)
#'
#' @param eez_sh
#' @param pred_his
#' @param pred_mod
#' @param wio
#'
#' @return
#' @export
#'

plot_predictions_in_eezs_extra <- function(eez_sh, pred_his, pred_mod, df_extra_his, df_extra_mod, wio){


  #convert extrapolation df to rasters
  extra_his <- raster::rasterFromXYZ(df_extra_his[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")
  extra_mod <- raster::rasterFromXYZ(df_extra_mod[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")


  #mask predictions with eez to extract eez
  pred_his2 <- raster::mask(pred_his, eez)
  pred_mod2 <- raster::mask(pred_mod, eez)


  #mask extrapolations with eez to extract eezs
  extra_his2 <- raster::mask(extra_his, eez)
  extra_mod2 <- raster::mask(extra_mod, eez)


  #convert prediction rasters to dataframes
  type <- "Historical"
  df_pred_his <- as(pred_his2, "SpatialPixelsDataFrame")
  df_pred_his <- as.data.frame(df_pred_his)
  colnames(df_pred_his) <- c("value", "x", "y")
  df_pred_his$type <- factor(c(type))

  type <- "Modern"
  df_pred_mod <- as(pred_mod2, "SpatialPixelsDataFrame")
  df_pred_mod <- as.data.frame(df_pred_mod)
  colnames(df_pred_mod) <- c("value", "x", "y")
  df_pred_mod$type <- factor(c(type))



  #convert extrapolation rasters to dataframes
  type <- "Historical"
  df_extra_his <- as(extra_his2, "SpatialPixelsDataFrame")
  df_extra_his <- as.data.frame(df_extra_his)
  colnames(df_extra_his) <- c("value", "x", "y")
  df_extra_his$type <- factor(c(type))

  type <- "Modern"
  df_extra_mod <- as(extra_mod2, "SpatialPixelsDataFrame")
  df_extra_mod <- as.data.frame(df_extra_mod)
  colnames(df_extra_mod) <- c("value", "x", "y")
  df_extra_mod$type <- factor(c(type))



  #select extrapolation data points (-1)
  df_extra_his <- subset(df_extra_his, value == -1)
  df_extra_mod <- subset(df_extra_mod, value == -1)




  #make map modern
  gHis <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_his, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_his, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    ggplot2::coord_sf(xlim = c(raster::extent(pred_his)[1], raster::extent(pred_his)[2]), ylim = c(raster::extent(pred_his)[3], raster::extent(pred_his)[4]), expand = FALSE) +
    ggspatial::fixed_plot_aspect(ratio = 1) +
    ggplot2::ylab("") +
    ggplot2::xlab("Historical") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))



  #make map historical
  gMod <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_mod, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_mod, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ggplot2::coord_sf(xlim = c(raster::extent(pred_mod)[1], raster::extent(pred_mod)[2]), ylim = c(raster::extent(pred_mod)[3], raster::extent(pred_mod)[4]), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("Modern") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))




  #violin plot
  pred <- rbind(df_pred_his, df_pred_mod)

  violin <- ggplot2::ggplot(pred, ggplot2::aes(x = type, y = value, color = type)) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1, fill = "white") +
    ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("#00BA38", "#F8766D"))


  g <- cowplot::ggdraw() +
    cowplot::draw_plot(gHis, x = -0.18, y = 0.1, height = .7, width = .7) +
    cowplot::draw_plot(gMod, x = .16, y = 0.1, height = .7, width = .7) +
    cowplot::draw_plot(violin, x = .69, y = 0.1, height = .7, width = .3) +
    cowplot::draw_label("All eezs", x = 0.5,  y = 0.92,  size = 26)

  ggplot2::ggsave(here::here("outputs", "plot_predictions_in_eezs_extra.jpeg"), g, width = 15.5, height = 5, bg = "white")

}




#' Clip and plot predictions in high seas versus eezs with violin plot (removing extrapolation zones)
#'
#' @param eez_sh
#' @param pred_his
#' @param pred_mod
#' @param wio
#'
#' @return
#' @export
#'

plot_predictions_high_seas_vs_eezs_extra_violin <- function(eez_sh, pred_his, pred_mod, df_extra_his, df_extra_mod, wio){

  ########################## high seas ##########################

  #convert extrapolation df to rasters
  extra_his <- raster::rasterFromXYZ(df_extra_his[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")
  extra_mod <- raster::rasterFromXYZ(df_extra_mod[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")


  #inverse mask predictions with eez to extract high seas
  pred_his_highsea <- raster::mask(pred_his, eez, inverse=T)
  pred_mod_highsea <- raster::mask(pred_mod, eez, inverse=T)


  #inverse mask extrapolations with eez to extract high seas
  extra_his_highsea <- raster::mask(extra_his, eez, inverse=T)
  extra_mod_highsea <- raster::mask(extra_mod, eez, inverse=T)


  #remove remaining pixels along coastlines for predictions
  wio2 <- sf::as_Spatial(wio)
  wio3 <- raster::buffer(wio2, width = 0.5)
  wio4 <- raster::crop(wio3, raster::extent(pred_his))
  pred_his_highsea <- raster::mask(pred_his_highsea, wio4, inverse=T) #these are the predictions extracted in high seas
  pred_mod_highsea <- raster::mask(pred_mod_highsea, wio4, inverse=T) #these are the predictions extracted in high seas


  #remove remaining pixels along coastlines for extrapolations
  extra_his_highsea <- raster::mask(extra_his_highsea, wio4, inverse=T) #these are the extrapolations extracted in high seas
  extra_mod_highsea <- raster::mask(extra_mod_highsea, wio4, inverse=T) #these are the extrapolations extracted in high seas



  #convert prediction rasters to dataframes
  type <- "Historical hs"
  df_pred_his_highsea <- as(pred_his_highsea, "SpatialPixelsDataFrame")
  df_pred_his_highsea <- as.data.frame(df_pred_his_highsea)
  colnames(df_pred_his_highsea) <- c("value", "x", "y")
  df_pred_his_highsea$type <- factor(c(type))

  type <- "Modern hs"
  df_pred_mod_highsea <- as(pred_mod_highsea, "SpatialPixelsDataFrame")
  df_pred_mod_highsea <- as.data.frame(df_pred_mod_highsea)
  colnames(df_pred_mod_highsea) <- c("value", "x", "y")
  df_pred_mod_highsea$type <- factor(c(type))



  #convert extrapolation rasters to dataframes
  type <- "Historical hs"
  df_extra_his_highsea <- as(extra_his_highsea, "SpatialPixelsDataFrame")
  df_extra_his_highsea <- as.data.frame(df_extra_his_highsea)
  colnames(df_extra_his_highsea) <- c("value", "x", "y")
  df_extra_his_highsea$type <- factor(c(type))

  type <- "Modern hs"
  df_extra_mod_highsea <- as(extra_mod_highsea, "SpatialPixelsDataFrame")
  df_extra_mod_highsea <- as.data.frame(df_extra_mod_highsea)
  colnames(df_extra_mod_highsea) <- c("value", "x", "y")
  df_extra_mod_highsea$type <- factor(c(type))



  #select extrapolation data points (-1)
  df_extra_his_highsea <- subset(df_extra_his_highsea, value == -1)
  df_extra_mod_highsea <- subset(df_extra_mod_highsea, value == -1)




  #make map modern
  gHis_highsea <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_his_highsea, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_his_highsea, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(raster::extent(pred_his_highsea)[1], raster::extent(pred_his_highsea)[2]), ylim = c(raster::extent(pred_his_highsea)[3], raster::extent(pred_his_highsea)[4]), expand = FALSE) +
    ggspatial::fixed_plot_aspect(ratio = 1) +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))



  #make map historical
  gMod_highsea <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_mod_highsea, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_mod_highsea, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ggplot2::coord_sf(xlim = c(raster::extent(pred_mod_highsea)[1], raster::extent(pred_mod_highsea)[2]), ylim = c(raster::extent(pred_mod_highsea)[3], raster::extent(pred_mod_highsea)[4]), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))




  ########################## eez ##########################

  #mask predictions with eez to extract eez
  pred_his_eez <- raster::mask(pred_his, eez)
  pred_mod_eez <- raster::mask(pred_mod, eez)


  #mask extrapolations with eez to extract eezs
  extra_his_eez <- raster::mask(extra_his, eez)
  extra_mod_eez <- raster::mask(extra_mod, eez)


  #convert prediction rasters to dataframes
  type <- "Historical eez"
  df_pred_his_eez <- as(pred_his_eez, "SpatialPixelsDataFrame")
  df_pred_his_eez <- as.data.frame(df_pred_his_eez)
  colnames(df_pred_his_eez) <- c("value", "x", "y")
  df_pred_his_eez$type <- factor(c(type))

  type <- "Modern eez"
  df_pred_mod_eez <- as(pred_mod_eez, "SpatialPixelsDataFrame")
  df_pred_mod_eez <- as.data.frame(df_pred_mod_eez)
  colnames(df_pred_mod_eez) <- c("value", "x", "y")
  df_pred_mod_eez$type <- factor(c(type))



  #convert extrapolation rasters to dataframes
  type <- "Historical eez"
  df_extra_his_eez <- as(extra_his_eez, "SpatialPixelsDataFrame")
  df_extra_his_eez <- as.data.frame(df_extra_his_eez)
  colnames(df_extra_his_eez) <- c("value", "x", "y")
  df_extra_his_eez$type <- factor(c(type))

  type <- "Modern eez"
  df_extra_mod_eez <- as(extra_mod_eez, "SpatialPixelsDataFrame")
  df_extra_mod_eez <- as.data.frame(df_extra_mod_eez)
  colnames(df_extra_mod_eez) <- c("value", "x", "y")
  df_extra_mod_eez$type <- factor(c(type))



  #select extrapolation data points (-1)
  df_extra_his_eez <- subset(df_extra_his_eez, value == -1)
  df_extra_mod_eez <- subset(df_extra_mod_eez, value == -1)




  #make map modern
  gHis_eez <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_his_eez, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_his_eez, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(raster::extent(pred_his_eez)[1], raster::extent(pred_his_eez)[2]), ylim = c(raster::extent(pred_his_eez)[3], raster::extent(pred_his_eez)[4]), expand = FALSE) +
    ggspatial::fixed_plot_aspect(ratio = 1) +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))



  #make map historical
  gMod_eez <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_mod_eez, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra_mod_eez, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ggplot2::coord_sf(xlim = c(raster::extent(pred_mod_eez)[1], raster::extent(pred_mod_eez)[2]), ylim = c(raster::extent(pred_mod_eez)[3], raster::extent(pred_mod_eez)[4]), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("") +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                   axis.title = ggplot2::element_text(size=16, face = "bold"),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))


  ############################## combined violin plot ######################################

  pred_highsea <- rbind(df_pred_his_highsea, df_pred_mod_highsea)
  pred_eez <- rbind(df_pred_his_eez, df_pred_mod_eez)

  cat("highseas to eez ratio in historical model", median(df_pred_his_highsea$value) / median(df_pred_his_eez$value), sep="\n")
  cat("highseas to eez ratio in modern model", median(df_pred_mod_highsea$value) / median(df_pred_mod_eez$value), sep="\n")

  violin <- ggplot2::ggplot(pred_highsea, ggplot2::aes(x = type, y = value, color = type)) +
    #violin plot for high seas
    ggplot2::geom_violin(width = 0.8) +
    ggplot2::geom_boxplot(width = 0.1, fill = "white") +
    #violin plot for eez
    ggplot2::geom_violin(data = pred_eez, ggplot2::aes(x = type, y = value, color = type), width = 0.8) +
    ggplot2::geom_boxplot(data = pred_eez, width = 0.1, fill = "white") +
    ggplot2::ylab("Habitat suitability") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.text.x = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(size = 18),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::annotate(geom = "text", x = 1.5, y = -0.02, size = 10, label = "Historical", color = "#00BA38") +
    ggplot2::annotate(geom = "text", x = 3.5, y = -0.02, size = 10, label = "Modern", color = "#F8766D") +
    ggplot2::annotate(geom = "text", x = 1, y = 1.025, size = 8, label = "All eezs") +
    ggplot2::annotate(geom = "text", x = 3, y = 1.025, size = 8, label = "All eezs") +
    ggplot2::annotate(geom = "text", x = 2, y = 1.025, size = 8, label = "High seas") +
    ggplot2::annotate(geom = "text", x = 4, y = 1.025, size = 8, label = "High seas") +
    ggplot2::scale_color_manual(values = c(rep("#00BA38", 2), rep("#F8766D", 2)))



  # multiplot
  jpeg(here::here("outputs", "plot_predictions_high_seas_vs_eezs_extra_violin.jpeg"), width = 1500, height = 900)
  cowplot::ggdraw() +
    cowplot::draw_plot(gHis_eez, x = -0.08, y = 0.47, height = .449, width = .449) +
    cowplot::draw_plot(gHis_highsea, x = 0.2, y = 0.47, height = .449, width = .449) +
    cowplot::draw_plot(gMod_eez, x = - 0.08, y = 0.02, height = .449, width = .449) +
    cowplot::draw_plot(gMod_highsea, x = 0.2, y = 0.02, height = .449, width = .449) +
    cowplot::draw_plot(violin, x = 0.565, y = 0.145, height = .7, width = .43) +
    cowplot::draw_label("High seas", x = 0.425,  y = 0.96,  size = 26) +
    cowplot::draw_label("All eezs", x = 0.15,  y = 0.96,  size = 26) +
    cowplot::draw_label("Historical", x = 0.01,  y = 0.71,  size = 26, angle = 90) +
    cowplot::draw_label("Modern", x = 0.01,  y = 0.27,  size = 26, angle = 90) +
    cowplot::draw_text("(a)", x = 0.01, y = 0.91, size = 25, fontface = "italic") +
    cowplot::draw_text("(b)", x = 0.29, y = 0.91, size = 25, fontface = "italic") +
    cowplot::draw_text("(c)", x = 0.01, y = 0.46, size = 25, fontface = "italic") +
    cowplot::draw_text("(d)", x = 0.29, y = 0.46, size = 25, fontface = "italic") +
    cowplot::draw_text("(e)", x = 0.58, y = 0.825, size = 25, fontface = "italic")
  dev.off()


}



#' Select threshold indicating highly suitable habitat (using max SSS following Liu et al 2013)
#'
#' @param occ
#' @param bg
#' @param modelstack
#'
#' @return
#' @export
#'

select_sss_threshold <- function(occ, bg.z, modelstack){

  # Create SWD object
  data <- SDMtune::prepareSWD(species = "SW", p = occ, a = bg.z[c("Lon", "Lat")],
                              env = modelstack)

  # Train model (needed to get threshold)
  model <- SDMtune::train(method = "Maxnet", data = data, fc = "l", reg = 1, iter = 500)

  # Get Maximum training sensitivity plus specificity threshold (sss)
  t <- SDMtune::thresholds(model, type = "logistic")

  sss <- t[3,2]

  return(sss)

}




#' Plot predictions above sss threshold with extrapolation extent and mpas
#'
#' @param wio
#' @param df_pred
#' @param type
#' @param df_test
#' @param mpa
#' @param sss
#'
#' @return
#' @export
#'

plot_predictions_with_extra_mpas_above_threshold <- function(wio, df_pred, type, df_test, mpa, sss){

  #select extrapolation data points
  df_test <- subset(df_test, cfact1 == -1)

  #select predictions above sss threshold
  df_pred_sss <- subset(df_pred, df_pred$value >= sss)

  g <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_pred_sss, ggplot2::aes(x = x, y = y, fill = value)) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_test, ggplot2::aes(x = x, y = y), fill = "grey30") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    #add mpas
    ggplot2::geom_sf(data = mpa, color = "black", fill = NA, size = 0.6) +
    #add countries
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab(paste(type, "model")) +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 17),
                   legend.title = ggplot2::element_text(size = 17),
                   axis.title = ggplot2::element_text(size = 18),
                   axis.text = ggplot2::element_text(size = 10),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA))

  ggplot2::ggsave(here::here("outputs", paste0("map_predictions_with_extra_", type, "_mpas_above_threshold.png")), g, width = 9, height = 7)

  return(g)

}






#' Mask distance to coast with predictions above threshold
#'
#' @param pred
#' @param df_extra
#' @param sss
#' @param dcoast
#'
#' @return
#' @export
#'

mask_dist_to_coast_predictions_above_threshold <- function(dcoast, pred, df_extra, sss){

  # set NA to predictions below sss
  pred <- pred
  pred[pred < sss] <- NA

  #convert extrapolation df to raster
  r_extra <- raster::rasterFromXYZ(df_extra[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")

  #mask predictions that are in extrapolation zone
  pred <- raster::mask(pred, r_extra, maskvalue = -1)

  #get area of predictions (nb of cells) above threshold
  print("area of predictions (nb of cells) above threshold")
  print(raster::ncell(pred[!is.na(pred)]))

  #mask distance to coast with these predictions
  dcoast <- raster::aggregate(dcoast, fact = 20)
  dcoast_masked <- raster::mask(dcoast, pred)

  raster::plot(dcoast_masked)

return(dcoast_masked)
}






#' Map high suitability predictions that have become low suitability and low suitability predictions that have become high suitability
#'
#' @param pred_mod
#' @param pred_his
#' @param sss_mod
#' @param sss_his
#' @param df_extra_mod
#' @param df_extra_his
#'
#' @return
#' @export
#'

plot_predictions_lost_gained <- function(pred_mod, pred_his, sss_mod, sss_his, df_extra_mod, df_extra_his, mpa){

  ##### get lost areas (high quality predictions in historical that have become low quality in modern)

  # set NA to historical predictions below sss to get high quality predictions
  pred_his1 <- pred_his
  pred_his1[pred_his1 < sss_his] <- NA

  # set NA to modern predictions above sss to get low quality predictions
  pred_mod1 <- pred_mod
  pred_mod1[pred_mod1 > sss_mod] <- NA

  # mask to get high quality historical area that have become low quality
  lost <- raster::mask(pred_his1, pred_mod1)
  lost[!is.na(lost)] <- 1

  # convert to dataframe for plotting
  df_lost <- as(lost, "SpatialPixelsDataFrame")
  df_lost <- as.data.frame(df_lost)
  colnames(df_lost) <- c("value", "x", "y")


  ##### get gained areas (low quality predictions in historical that have become high quality in modern)

  # set NA to historical predictions above sss to get low quality predictions
  pred_his1 <- pred_his
  pred_his1[pred_his1 >= sss_his] <- NA

  # set NA to modern predictions below sss to get high quality predictions
  pred_mod1 <- pred_mod
  pred_mod1[pred_mod1 <= sss_mod] <- NA

  # mask to get low quality historical area that have become high quality
  gained <- raster::mask(pred_his1, pred_mod1)
  gained[!is.na(gained)] <- 1

  # convert to dataframe for plotting
  df_gained <- as(gained, "SpatialPixelsDataFrame")
  df_gained <- as.data.frame(df_gained)
  colnames(df_gained) <- c("value", "x", "y")

  ##### prepare modern and historical extrapolation dataframe
  df_extra_mod <- subset(df_extra_mod, cfact1 == -1)
  df_extra_his <- subset(df_extra_his, cfact1 == -1)
  df_extra <- rbind(df_extra_mod, df_extra_his)

  g <- ggplot2::ggplot() +
    #lost
    ggplot2::geom_tile(data = df_lost, ggplot2::aes(x = x, y = y, fill = "lost")) +
    #gained
    ggplot2::geom_tile(data = df_gained, ggplot2::aes(x = x, y = y, fill = "gained")) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_extra, ggplot2::aes(x = x, y = y), fill = "grey30") +
    #add mpas
    ggplot2::geom_sf(data = mpa, color = "black", fill = NA, size = 0.6) +
    #add countries
    ggplot2::geom_sf(data = wio, color = "white", fill = "grey80", size = 0.2) +
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::scale_fill_manual(values = c("lost" = "orange", "gained" ="blue"),
                               labels = c("Lost high \nsuitability", "Gained high \nsuitability"),
                               name = "") +
    ggplot2::theme(legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 14),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 10),
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA))

  ggplot2::ggsave(here::here("outputs", "map_predictions_lost_gained.png"), g, width = 9, height = 7)

  return(g)

}











#' Mask distance to coast with lost predictions
#'
#' @param dcoast
#' @param pred_mod
#' @param pred_his
#' @param sss_mod
#' @param sss_his
#' @param df_extra_mod
#' @param df_extra_his
#'
#' @return
#' @export
#'

mask_dist_to_coast_predictions_lost <- function(dcoast, pred_mod, pred_his, sss_mod, sss_his, df_extra_mod, df_extra_his){


  ##### prepare modern and historical extrapolation dataframe
  df_extra_mod <- subset(df_extra_mod, cfact1 == -1)
  df_extra_his <- subset(df_extra_his, cfact1 == -1)
  df_extra <- rbind(df_extra_mod, df_extra_his)

  #convert extrapolation df to raster
  r_extra <- raster::rasterFromXYZ(df_extra[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")


  ##### get lost areas (high quality predictions in historical that have become low quality in modern)

  # set NA to historical predictions below sss to get high quality predictions
  pred_his1 <- pred_his
  pred_his1[pred_his1 < sss_his] <- NA

  # set NA to modern predictions above sss to get low quality predictions
  pred_mod1 <- pred_mod
  pred_mod1[pred_mod1 > sss_mod] <- NA

  # mask to get high quality historical area that have become low quality
  lost <- raster::mask(pred_his1, pred_mod1)
  lost[!is.na(lost)] <- 1

  # mask lost predictions that are in extrapolation zone
  lost <- raster::mask(lost, r_extra, maskvalue = -1)

  # get area of predictions (nb of cells) above threshold
  print("area of predictions for lost locations (nb of cells)")
  print(raster::ncell(lost[!is.na(lost)]))

  # mask distance to coast with these predictions
  dcoast <- raster::aggregate(dcoast, fact = 20)
  dcoast_masked <- raster::mask(dcoast, lost)

  raster::plot(dcoast_masked)

  return(dcoast_masked)

}







#' Mask distance to coast with gained predictions
#'
#' @param dcoast
#' @param pred_mod
#' @param pred_his
#' @param sss_mod
#' @param sss_his
#' @param df_extra_mod
#' @param df_extra_his
#'
#' @return
#' @export
#'

mask_dist_to_coast_predictions_gained <- function(dcoast, pred_mod, pred_his, sss_mod, sss_his, df_extra_mod, df_extra_his){


  ##### prepare modern and historical extrapolation dataframe
  df_extra_mod <- subset(df_extra_mod, cfact1 == -1)
  df_extra_his <- subset(df_extra_his, cfact1 == -1)
  df_extra <- rbind(df_extra_mod, df_extra_his)

  #convert extrapolation df to raster
  r_extra <- raster::rasterFromXYZ(df_extra[,c("x", "y", "cfact1")], crs = "+proj=longlat +datum=WGS84 +no_defs")


  ##### get gained areas (high quality predictions in historical that have become low quality in modern)

  # set NA to historical predictions above sss to get low quality predictions
  pred_his1 <- pred_his
  pred_his1[pred_his1 >= sss_his] <- NA

  # set NA to modern predictions below sss to get high quality predictions
  pred_mod1 <- pred_mod
  pred_mod1[pred_mod1 <= sss_mod] <- NA

  # mask to get low quality historical area that have become high quality
  gained <- raster::mask(pred_his1, pred_mod1)
  gained[!is.na(gained)] <- 1

  # mask gained predictions that are in extrapolation zone
  gained <- raster::mask(gained, r_extra, maskvalue = -1)

  # get area of predictions (nb of cells) above threshold
  print("area of predictions for gained locations (nb of cells)")
  print(raster::ncell(gained[!is.na(gained)]))

  # mask distance to coast with these predictions
  dcoast <- raster::aggregate(dcoast, fact = 20)
  dcoast_masked <- raster::mask(dcoast, gained)

  raster::plot(dcoast_masked)

  return(dcoast_masked)

}
