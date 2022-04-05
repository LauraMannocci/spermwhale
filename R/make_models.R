
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

select_best_model <- function(mod){

  # get overall results
  res <- ENMeval::eval.results(mod)

  # select models minimizing OP and maximizing auc
  opt <- res %>%
    dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
    dplyr::filter(auc.val.avg == max(auc.val.avg))

  return(opt)

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




#' Make violin plot of predictors for occurrences and background points
#'
#' @param modhis
#' @param var
#' @param ylabel
#'
#' @return
#' @export
#'

predictor_violin <- function(modhis, var, ylabel) {

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
    ggplot2::stat_summary(fun = mean, color = "black", geom="text", ggplot2::aes(label=round(..y.., digits=0))) +
    ggplot2::ggtitle(signif)+
    ggplot2::ylab(ylabel) +
    ggplot2::theme_light() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("#00BA38", "#F8766D"))

  ggplot2::ggsave(here::here("outputs", paste0("predictors_violin_", var, ".png")), p, width = 9, height = 7)

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
                           mc.cores = parallel::detectCores()/1.5,
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







#' Plot residuals with respect to combined model
#'
#' @param pred
#' @param predcom
#' @param wio
#' @param type
#'
#' @return
#' @export
#'

plot_residuals_wrt_combined <- function(pred, predcom, wio, type){

  model <- lm(raster::getValues(pred) ~ raster::getValues(predcom), na.action=na.exclude)
  res <- pred
  res[] <- residuals.lm(model)
  df_res <- as(res, "SpatialPixelsDataFrame")
  df_res <- as.data.frame(df_res)
  colnames(df_res) <- c("value", "x", "y")

  gg <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = wio) +
    ggplot2::geom_tile(data = df_res, ggplot2::aes(x = x, y = y, fill = value), alpha=0.8) +
    ggplot2::scale_fill_viridis_c(option = "inferno", limit = c(-0.34,0.71)) +
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("") +
    ggplot2::xlab(paste(type, 'residuals')) +
    ggplot2::theme(legend.position ="none")

  ggplot2::ggsave(here::here("outputs", paste0("residuals_", type, ".png")), gg, width = 9, height = 7)

  return(gg)

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




#' Plot predictions difference
#'
#' @param wio
#' @param df_pred
#' @param type
#'
#' @return
#' @export
#'

plot_predictions_difference <- function(wio, df_predHis, df_predMod){

  df_predMod$diff = df_predMod$value - df_predHis$value

  g <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = wio) +
    ggplot2::geom_tile(data = df_predMod, ggplot2::aes(x = x, y = y, fill = diff), alpha=0.8) +
    ggplot2::scale_fill_viridis_c(limits = c(-0.5, 1), option = "magma")+
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab("") +
    ggplot2::labs(fill = 'Difference \nmodern - \nhistorical')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))

  ggplot2::ggsave(here::here("outputs", "map_predictions_difference.png"), g, width = 9, height = 7)

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
    ggplot2::geom_sf(data = wio) +
    ggplot2::geom_tile(data = df_pred, ggplot2::aes(x = x, y = y, fill = value), alpha=0.8) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_test, ggplot2::aes(x = x, y = y), fill = "black", alpha=0.5) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab(paste(type, "model")) +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))

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
    ggplot2::geom_sf(data = wio) +
    ggplot2::geom_tile(data = df_pred, ggplot2::aes(x = x, y = y, fill = value), alpha=0.8) +
    #ask extrapolation mask
    ggplot2::geom_tile(data = df_test, ggplot2::aes(x = x, y = y), fill = "black", alpha=0.5) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
    #add mpas
    ggplot2::geom_sf(data = mpa, fill = "white", alpha = 0.1) +
    ggplot2::coord_sf(xlim = c(26, 85), ylim = c(-40, 25), expand = FALSE) +
    ggplot2::ylab("")+
    ggplot2::xlab(paste(type, "model")) +
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))

  ggplot2::ggsave(here::here("outputs", paste0("map_predictions_with_extra_", type, "_mpas.png")), g, width = 9, height = 7)

  return(g)

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






#' Clip and plot predictions in eez and add violin plot to compare them
#'
#' @param eez_sh ...
#' @param eez_name ...
#' @param pred_Historical ...
#' @param pred_Modern ...
#' @param wio ...
#'
#' @return
#' @export
#'

plot_predictions_in_eez <- function(eez_sh, eez_name, pred_Historical, pred_Modern, obs_mod, obs_his, wio){

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
  predHis <- raster::mask(pred_Historical, eez)
  predMod <- raster::mask(pred_Modern, eez)

  ###Historical

  #convert to dataframe
  type <- "Historical"
  df_predHis <- as(predHis, "SpatialPixelsDataFrame")
  df_predHis <- as.data.frame(df_predHis)
  colnames(df_predHis) <- c("value", "x", "y")
  df_predHis$type <- factor(c(type))

  #make map
  gHis <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = wio) +
    ggplot2::geom_tile(data = df_predHis, ggplot2::aes(x = x, y = y, fill = value), alpha=0.8) +
    ggplot2::geom_point(data = obs_his, ggplot2::aes(x = Lon, y = Lat, alpha = 0.6)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis")+
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
  df_predMod <- as(predMod, "SpatialPixelsDataFrame")
  df_predMod <- as.data.frame(df_predMod)
  colnames(df_predMod) <- c("value", "x", "y")
  df_predMod$type <- factor(c(type))


  #make map
  gMod <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = wio) +
    ggplot2::geom_tile(data = df_predMod, ggplot2::aes(x = x, y = y, fill = value), alpha=0.8) +
    ggplot2::geom_point(data = obs_mod, ggplot2::aes(x = Lon, y = Lat, alpha = 0.6)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
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
  pred <- rbind(df_predHis, df_predMod)
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

