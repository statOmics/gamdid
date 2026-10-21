#' @importFrom S4Vectors metadata
#' @importFrom dplyr filter
#' @importFrom stats predict qnorm
#' @importFrom ggplot2 ggplot geom_col aes geom_line labs theme_bw theme
#' @importFrom ggplot2 element_text element_line element_blank geom_rect
#' @importFrom ggplot2 geom_hline facet_wrap ggtitle scale_color_manual
#' @importFrom ggplot2 scale_fill_manual %+replace%
#' @importFrom grid unit
#' @importFrom utils globalVariables
#' @importFrom mgcv gam

#' @include createCounts.R

utils::globalVariables(c(
    "covariate", "midsC", "carrier", ".data",
    "fit", "xmin", "xmax", "grid", "lwr", "upr"
))

PALETTE_GAMDID <- c("#D55E00", "#56B4E9")

visFitSCE <- function(sce, feature, plot = "hist", shade = TRUE,
    alpha = NULL) {
    if (!is.character(plot) || length(plot) != 1 ||
        !plot %in% c("hist", "diff")) {
        stop("'plot' must be \"hist\" or \"diff\"")
    }
    covar <- metadata(sce)$covariate
    bins <- metadata(sce)$bins
    comps <- metadata(sce)$comparisons

    feature_fit <- fitFeature(sce, feature, covar, bins)
    dat <- feature_fit$dat
    fit_gam <- feature_fit$fit_gam

    plot_list <- vector(mode = "list", length = length(comps))
    names(plot_list) <- comps

    ## Two-sided significance level, the same for every comparison
    alpha_half <- visFitAlpha(sce, alpha) / 2

    for (comp in comps) {
        ## Strip the covariate name from the start of both level names
        comp_levels <- substring(
            strsplit(comp, "_VS_", fixed = TRUE)[[1]],
            nchar(covar) + 1
        )
        hist_dat <- dat |> filter(covariate %in% comp_levels)

        band <- diffSmoother(fit_gam, hist_dat, comp_levels, bins, alpha_half)
        newdat <- densityCurves(fit_gam, hist_dat, comp_levels)
        colnames(newdat)[colnames(newdat) == "covariate"] <- covar
        colnames(hist_dat)[colnames(hist_dat) == "covariate"] <- covar

        if (plot == "hist") {
            plot_list[[comp]] <- histPlot(
                hist_dat, newdat, band, covar, comp, shade
            )
        } else {
            plot_list[[comp]] <- diffPlot(
                hist_dat, newdat, band, covar, comp, shade
            )
        }
    }

    suppressWarnings(print(plot_list))
    invisible(plot_list)
}

## Histogram counts and Poisson GAM of one feature, on the partial residuals
## when a nuisance was specified in compDistr
fitFeature <- function(sce, feature, covar, bins) {
    if ("partialResiduals" %in% names(assays(sce))) {
        assay_dat <- assays(sce)$partialResiduals[feature, ]
    } else {
        assay_dat <- assay(sce)[feature, ]
    }

    dat <- createCounts(
        x = assay_dat, colDat = colData(sce),
        covariate = covar, bins = bins, minObs = 0
    )

    fit_gam <- gam(
        counts ~ covariate + s(midsC, by = covariate) + offset(carrier),
        data = dat, family = poisson(), method = "REML"
    )
    list(dat = dat, fit_gam = fit_gam)
}

## Significance level of the confidence bands: alpha when given, else the
## BH threshold of compDistr, or the most stringent BH level when no feature
## is significant
visFitAlpha <- function(sce, alpha) {
    if (is.null(alpha)) {
        alpha <- metadata(sce)$adj_alpha
        if (is.na(alpha)) {
            ## No significant feature: use the most stringent BH level
            alpha_fdr <- metadata(sce)$alpha
            n_tests <- metadata(sce)$n_tests
            if (is.null(alpha_fdr) || is.null(n_tests) || n_tests == 0) {
                stop(
                    "No feature is significant and the BH threshold ",
                    "cannot be derived: please specify 'alpha'."
                )
            }
            alpha <- alpha_fdr / n_tests
            message(
                "No feature is significant at FDR ", alpha_fdr,
                "; confidence bands use the most stringent BH level ",
                "alpha / ", n_tests, " = ", signif(alpha, 2), "."
            )
        }
    }
    alpha
}

## Difference smoother of two groups with its confidence band on a fine
## grid, and the regions where the band excludes zero (NULL if none)
diffSmoother <- function(fit_gam, hist_dat, comp_levels, bins, alpha_half) {
    grid <- seq(min(hist_dat$midsC), max(hist_dat$midsC), length.out = 100)
    newdatA <- data.frame(
        midsC = grid,
        carrier = hist_dat$carrier[1],
        covariate = comp_levels[1]
    )
    newdatB <- data.frame(
        midsC = grid,
        carrier = hist_dat$carrier[bins + 1],
        covariate = comp_levels[2]
    )
    xA <- predict(fit_gam, newdata = newdatA, type = "lpmatrix")
    xB <- predict(fit_gam, newdata = newdatB, type = "lpmatrix")
    xdiff <- xB - xA
    diff <- as.vector(xdiff %*% fit_gam$coefficients)
    sediff <- sqrt(rowSums(xdiff %*% fit_gam$Vp * xdiff))

    ## Sidak-corrected pointwise significance threshold
    crit <- qnorm(1 - (1 - alpha_half)^(1 / max(summary(fit_gam)$edf)))
    lwr <- diff + crit * sediff
    upr <- diff - crit * sediff

    mids_gray <- grid[!((lwr <= 0) & (upr >= 0))]
    grid_bw <- grid[2] - grid[1]
    shade_df <- NULL
    if (length(mids_gray) > 0) {
        shade_df <- data.frame(
            xmin = mids_gray - grid_bw / 2,
            xmax = mids_gray + grid_bw / 2
        )
    }
    list(grid = grid, diff = diff, lwr = lwr, upr = upr, shade_df = shade_df)
}

## Fitted GAM densities of two groups (carrier = 0 gives the density scale)
densityCurves <- function(fit_gam, hist_dat, comp_levels) {
    newdat <- data.frame(
        midsC = rep(seq(min(hist_dat$midsC), max(hist_dat$midsC),
            length.out = 100
        ), each = 2),
        covariate = comp_levels,
        carrier = 0
    )
    newdat$fit <- predict(fit_gam, newdata = newdat, type = "response")
    newdat
}

gamdidTheme <- function(base_size = 10, base_family = "sans") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            axis.line = element_line(linewidth = 0.2, color = "black"),
            axis.ticks = element_line(linewidth = 0.2, color = "black"),
            axis.text = element_text(size = base_size, color = "black"),
            axis.title = element_text(size = base_size, color = "black"),
            legend.title = element_text(size = base_size, face = "bold"),
            legend.text = element_text(size = base_size),
            legend.key.size = unit(0.2, "cm"),
            legend.background = element_blank(),
            legend.position = "top",
            legend.direction = "horizontal"
        )
}

## Layers of the histograms and the fitted GAM densities
densityLayers <- function(hist_dat, newdat, covar) {
    list(
        geom_col(
            data = hist_dat,
            aes(
                x = midsC, y = counts / exp(carrier),
                fill = .data[[covar]]
            ),
            alpha = 0.4, position = "identity", colour = "black"
        ),
        scale_fill_manual(values = PALETTE_GAMDID),
        geom_line(
            data = newdat,
            mapping = aes(x = midsC, y = fit, colour = .data[[covar]]),
            linetype = 1, linewidth = 1.2
        ),
        scale_color_manual(values = PALETTE_GAMDID)
    )
}

## Grey rectangles over the regions where the groups differ
shadeLayer <- function(shade_df, ymin) {
    geom_rect(
        data = shade_df,
        aes(xmin = xmin, xmax = xmax),
        ymin = ymin, ymax = Inf, inherit.aes = FALSE,
        fill = "gray", alpha = 0.3
    )
}

## Histograms and fitted GAM densities, with shaded differential regions
histPlot <- function(hist_dat, newdat, band, covar, comp, shade) {
    p <- ggplot() +
        densityLayers(hist_dat, newdat, covar) +
        labs(x = "midpoints", y = "density", title = comp) +
        gamdidTheme()

    if (!is.null(band$shade_df) && shade) {
        p <- p + shadeLayer(band$shade_df, ymin = 0)
    }
    p
}

## Layers of the difference smoother and its confidence band
diffLayers <- function(band) {
    df_plot <- data.frame(
        grid = band$grid, diff = band$diff, lwr = band$lwr, upr = band$upr,
        panel = "log(difference)"
    )
    list(
        geom_line(data = df_plot, aes(x = grid, y = diff)),
        geom_line(data = df_plot, aes(x = grid, y = lwr), linetype = "dashed"),
        geom_line(data = df_plot, aes(x = grid, y = upr), linetype = "dashed"),
        geom_hline(
            data = data.frame(panel = "log(difference)"),
            aes(yintercept = 0),
            linetype = "dotted", color = "red"
        )
    )
}

## Two panels: histograms with GAM densities, and the difference smoother
## with its confidence band
diffPlot <- function(hist_dat, newdat, band, covar, comp, shade) {
    p <- ggplot() +
        densityLayers(
            cbind(hist_dat, panel = "density"),
            cbind(newdat, panel = "density"), covar
        ) +
        diffLayers(band) +
        facet_wrap(~panel,
            ncol = 1, scales = "free_y",
            strip.position = "left"
        ) +
        labs(x = "x") +
        ggtitle(comp) +
        gamdidTheme() %+replace%
        theme(
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 10, color = "black"),
            strip.placement = "outside"
        )

    if (!is.null(band$shade_df) && shade) {
        p <- p +
            shadeLayer(cbind(band$shade_df, panel = "density"), ymin = 0) +
            shadeLayer(
                cbind(band$shade_df, panel = "log(difference)"),
                ymin = -Inf
            )
    }
    p
}
