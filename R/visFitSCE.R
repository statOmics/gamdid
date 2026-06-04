#' @importFrom MultiAssayExperiment metadata
#' @importFrom dplyr filter
#' @importFrom stats predict qnorm
#' @importFrom ggplot2 ggplot geom_col aes geom_line labs theme_bw theme element_text element_line element_blank geom_rect geom_hline facet_wrap ggtitle scale_color_manual scale_fill_manual %+replace%
#' @importFrom grid unit
#' @importFrom utils globalVariables
#' @importFrom mgcv gam

#' @include createCounts.R

utils::globalVariables(c("covariate", "midsC", "carrier", ".data", "fit", "xmin", "xmax"))

visFitSCE <- function(sce, feature, plot = "hist", shade = TRUE, alpha = NULL){
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

  fit_gam <- gam(
    counts ~ covariate + s(midsC, by = covariate) + offset(carrier),
    data = dat,
    family = poisson(),
    method = "REML")

  plot_list <- vector(mode='list', length=length(comps))
  names(plot_list) <- comps
  for (comp in comps){
    comp_levels <- gsub(covar,"",unlist(strsplit(comp, "_VS_")))
    hist_dat <- dat |> filter(covariate %in% comp_levels)
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

  if (is.null(alpha)){
    alpha <- metadata(sce)$adj_alpha/2 # two-sided test
  } else {
    alpha <- alpha/2 # two-sided test
  }


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

  base_size   <- 10
  base_family <- "sans"

  PALETTE_GAMDID <- c(
    "#D55E00", "#56B4E9"
  )

  scale_color_gamdid <- function(...) {
    scale_color_manual(values = PALETTE_GAMDID, ...)
  }

  scale_fill_gamdid <- function(...) {
    scale_fill_manual(values = PALETTE_GAMDID, ...)
  }

  if (plot == "hist"){
    ggplot() +
      geom_col(aes(x = midsC, y = counts/exp(carrier), fill = .data[[covar]]), alpha = 0.4, position = "identity", colour = "black",data = hist_dat) +
      scale_fill_gamdid() +
      geom_line(data = newdat, mapping = aes(x = midsC, y = fit, colour = .data[[covar]]), linetype = 1, linewidth = 1.2) +
      scale_color_gamdid() +
      labs(x = "midpoints", y = "density", title = comp) +
      theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            # Axes
            axis.line = element_line(linewidth = 0.2, color = "black"),
            axis.ticks = element_line(linewidth = 0.2, color = "black"),
            axis.text        = element_text(size = base_size, color = "black"),
            axis.title       = element_text(size = base_size, color = "black"),
            # Legend
            legend.title     = element_text(size = base_size, face = "bold"),
            legend.text      = element_text(size = base_size),
            legend.key.size  = unit(0.2, "cm"),
            legend.background = element_blank(),
            legend.position  = "top",
            legend.direction = "horizontal") -> plot_list[[comp]]

    if ( (length(mids_gray) > 0) & shade ){
      plot_list[[comp]] +
        geom_rect(data = shade_df, aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = Inf, inherit.aes = F,
                  fill = "gray", alpha = 0.3) -> plot_list[[comp]]
    }

  } else if (plot == "diff"){
    df <- data.frame(grid = grid, diff = diff, lwr = lwr, upr = upr)

    ggplot() +
      # histogram layer
      geom_col(data = cbind(hist_dat, panel = "density"),
               aes(x = midsC, y = counts/exp(carrier), fill = .data[[covar]]),
               alpha = 0.4, position = "identity", colour = "black") +
      scale_fill_gamdid() +
      geom_line(data = cbind(newdat, panel = "density"),
                mapping = aes(x = midsC, y = fit, colour = .data[[covar]]),
                linetype = 1, linewidth = 1.2) +
      scale_color_gamdid() +
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
      labs(x = "x") +
      ggtitle(comp) +
      theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            # Axes
            axis.line = element_line(linewidth = 0.2, color = "black"),
            axis.ticks = element_line(linewidth = 0.2, color = "black"),
            axis.text        = element_text(size = base_size, color = "black"),
            axis.title       = element_text(size = base_size, color = "black"),
            axis.title.y = element_blank(),
            # Legend
            legend.title     = element_text(size = base_size, face = "bold"),
            legend.text      = element_text(size = base_size),
            legend.key.size  = unit(0.2, "cm"),
            legend.background = element_blank(),
            legend.position  = "top",
            legend.direction = "horizontal",
            # Facets
            strip.background = element_blank(),
            strip.text       = element_text(size = base_size, color = "black"),
            strip.placement = "outside") -> histdiff_plot

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
