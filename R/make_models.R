
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
  di_1000m_rg <- range(as.data.frame(bg)["di_1000m"], na.rm = TRUE)
  di_2500m_rg <- range(as.data.frame(bg)["di_2500m"], na.rm = TRUE)
  di_5000m_rg <- range(as.data.frame(bg)["di_5000m"], na.rm = TRUE)
  di_coast_rg <- range(as.data.frame(bg)["di_coast"], na.rm = TRUE)
  di_seaM_rg <- range(as.data.frame(bg)["di_seaM"], na.rm = TRUE)
  di_guy_rg <- range(as.data.frame(bg)["di_guy"], na.rm = TRUE)
  di_Trough_rg <- range(as.data.frame(bg)["di_Trough"], na.rm = TRUE)
  di_Tren_rg <- range(as.data.frame(bg)["di_Tren"], na.rm = TRUE)
  di_spRid_rg <- range(as.data.frame(bg)["di_spRid"], na.rm = TRUE)

  modelStack_limited[["depth"]] <- raster::clamp(modelStack_limited[["depth"]], lower = depth_rg[1], upper = depth_rg[2], useValues = FALSE)
  modelStack_limited[["slope"]] <- raster::clamp(modelStack_limited[["slope"]], lower = slope_rg[1], upper = slope_rg[2], useValues = FALSE)
  modelStack_limited[["di_1000m"]] <- raster::clamp(modelStack_limited[["di_1000m"]], lower = di_1000m_rg[1], upper = di_1000m_rg[2], useValues = FALSE)
  modelStack_limited[["di_2500m"]] <- raster::clamp(modelStack_limited[["di_2500m"]], lower = di_2500m_rg[1], upper = di_2500m_rg[2], useValues = FALSE)
  modelStack_limited[["di_5000m"]] <- raster::clamp(modelStack_limited[["di_5000m"]], lower = di_5000m_rg[1], upper = di_5000m_rg[2], useValues = FALSE)
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
    ggplot2::xlab("Combined model")+
    ggplot2::labs(fill = 'Habitat \nsuitability')+
    ggplot2::theme(legend.position = "right",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(-6,-10,-6,-10))

  ggplot2::ggsave(here::here("outputs", paste0("predictions_", type, ".png")), g, width = 9, height = 7)

  return(g)

}




#' test for "residual" mpa effect
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
  df_mpamod$type <- factor(c('mpa'))
  df_pred$type <- factor(c('region'))
  df_mpamod_pred <- rbind(df_pred, df_mpamod)
  df_mpamod_pred$model <- factor(name)

  print(kruskal.test(value ~ type, data = df_mpamod_pred))

  return(df_mpamod_pred)

}

