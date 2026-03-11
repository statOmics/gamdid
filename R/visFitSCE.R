#' Visualise the fit of one feature
#'
#' This functions visualizes the histogram and the fitted GAMs used for the
#' distributional analysis. If there is a significant difference in distribution,
#' the shaded regions pinpoint regions that drive the significant difference.
#'
#'
#' @param sce Single cell experiment object
#' @param feature The name of feature (or the number of the row of the feature)
#' that should be visualized
#' @param plot hist or diff, specify the plot you want
#' @return Interpretable visualisation
#' @export
visFitSCE <- function(sce, feature, plot = "hist", shade = TRUE){
  covar <- metadata(sce)$covariate
  bins <- metadata(sce)$bins
  comps <-  metadata(sce)$comparisons
  if ("partialResiduals" %in% names(assays(sce))){
    assay_dat <- assays(sce)$partialResiduals[feature,]
  } else {
    assay_dat <- assay(sce)[feature,]
  }

  dat <- createCounts(x = assay_dat, colDat = colData(sce),
                           covariate = covar, bins = bins, minObs = 0)

  fit_gam <- mgcv::gam(
    counts ~ covariate + s(midsC, by = covariate) + offset(carrier),
    data = dat,
    family = poisson(),
    method = "REML")

  plot_list <- vector(mode='list', length=length(comps))
  names(plot_list) <- comps
  for (comp in comps){
    comp_levels <- gsub(covar,"",unlist(strsplit(comp, "_VS_")))
    hist_dat <- dat %>% filter(covariate %in% comp_levels)
  # Create difference smoother and confidence bands
  ## create an even grid
  grid <- seq(min(hist_dat$midsC), max(hist_dat$midsC), length.out = 100)
  ## create model matrices for the difference smoother
  newdatA <- data.frame(midsC = grid, carrier = hist_dat$carrier[1], covariate = comp_levels[1])
  newdatB <- data.frame(midsC = grid, carrier = hist_dat$carrier[bins+1], covariate = comp_levels[2])
  xA <- predict(fit_gam, newdata = newdatA, type = "lpmatrix")
  xB <- predict(fit_gam, newdata = newdatB, type = "lpmatrix")
  xdiff <- xB - xA
  ## create the difference smoother
  diff <- as.vector(xdiff %*% fit_gam$coefficients)
  ## construct the confidence band
  sediff <- sqrt(rowSums(xdiff %*% fit_gam$Vp * xdiff))
  alpha <- metadata(sce)$adj_alpha/2 # two-sided test


  # lwr <- diff + qnorm(alpha/qr(xdiff)$rank)*sediff
  # upr <- diff - qnorm(alpha/qr(xdiff)$rank)*sediff
  lwr <- diff + qnorm(1-(1-alpha)^(1/max(summary(fit_gam)$edf)))*sediff
  upr <- diff - qnorm(1-(1-alpha)^(1/max(summary(fit_gam)$edf)))*sediff


  mids_gray <- grid[which(!((lwr <= 0) & (upr >= 0)))]
  grid_bw <- grid[2]- grid[1]
  if (length(mids_gray)>0){
    shade_df <- data.frame(xmin = mids_gray - grid_bw/2, xmax = mids_gray + grid_bw/2)
  }
  # Create plot
  newdat <- data.frame(midsC = rep(seq(min(hist_dat$midsC),max(hist_dat$midsC),length.out = 100), each = 2),
                       covariate = comp_levels,
                       #carrier = rep(c(hist_dat$'offset(carrier)'[1], hist_dat$'offset(carrier)'[bins+1]), each = 100))
                       carrier = 0)
  newdat$fit <- predict(fit_gam, newdata = newdat, type = "response")
  colnames(newdat)[which(colnames(newdat) == "covariate")] <- covar
  colnames(hist_dat)[which(colnames(hist_dat) == "covariate")] <- covar
  if (plot == "hist"){
    ggplot() +
      geom_col(aes(x = midsC, y = counts/exp(carrier), fill = .data[[covar]]), alpha = 0.3, position = "identity", colour = "black",data = hist_dat) +
      scale_fill_manual(values = c("salmon", "skyblue")) +
      geom_line(data = newdat, mapping = aes(x = midsC, y = fit, colour = .data[[covar]]), linetype = 1, linewidth = 1.2) +
      labs(x = "midpoints", y = "density", title = comp) +
      theme_bw() +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) -> plot_list[[comp]]

    if ( (length(mids_gray) > 0) & shade ){
      plot_list[[comp]] +
        geom_rect(data = shade_df, aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = Inf, inherit.aes = F,
                  fill = "gray", alpha = 0.3) + theme_bw() +
        theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) -> plot_list[[comp]]
    }

  } else if (plot == "diff"){
    # Combine data into one dataframe
    hist_dat_combined <- hist_dat |>
      mutate(panel = "density",
             x = midsC,
             y = counts/exp(carrier),
             lwr = NA,
             upr = NA) |>
      select(x,y,.data[[covar]],lwr,upr, panel)

    df <- data.frame(grid = grid, diff = diff, lwr = lwr, upr = upr)
    diff_dat_combined <- df |>
      mutate(panel = "log(difference)",
             x = grid,
             y = diff,
             "{covar}" := NA) |>
      select(x,y,.data[[covar]],lwr,upr, panel)

    combined_data <- rbind(diff_dat_combined, hist_dat_combined)

    ggplot() +
      # histogram layer
      geom_col(data = cbind(hist_dat, panel = "density"),
               aes(x = midsC, y = counts/exp(carrier), fill = .data[[covar]]),
               alpha = 0.3, position = "identity", colour = "black") +
      geom_line(data = cbind(newdat, panel = "density"),
                mapping = aes(x = midsC, y = fit, colour = .data[[covar]]),
                linetype = 1, linewidth = 1.2) +
      # difference layer
      geom_line(data = cbind(df, panel = "log(difference)"),
                aes(x = grid, y = diff)) +
      geom_line(data = cbind(df, panel = "log(difference)"),
                aes(x = grid, y = lwr), linetype = "dashed") +
      geom_line(data = cbind(df, panel = "log(difference)"),
                aes(x = grid, y = upr), linetype = "dashed") +
      geom_hline(data = data.frame(panel = "log(difference)"), aes(yintercept = 0),
                 linetype = "dotted", color = "red") +
      facet_wrap(~panel, ncol = 1, scales = "free_y", strip.position = "left") +
      theme_bw() +
      labs(x = "x") +
      ggtitle(comp) +
      theme(strip.placement = "outside",
            strip.background = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "top",
            plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) -> histdiff_plot

    if ( (length(mids_gray) > 0) & shade ){
      histdiff_plot +
        geom_rect(data = cbind(shade_df, panel = "density"), aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = Inf, inherit.aes = F,
                  fill = "gray", alpha = 0.3) +
        geom_rect(data = cbind(shade_df, panel = "log(difference)"), aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf, inherit.aes = F,
                  fill = "gray", alpha = 0.3) -> histdiff_plot
      }
    plot_list[[comp]] <- histdiff_plot
    }
  }
  suppressWarnings(print(plot_list))
  invisible(plot_list)
}
