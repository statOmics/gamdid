#' @importFrom MSqRob getContrast getFitMethod getCoef

topFeaturesOmnibus <-
  function (models, contrasts, adjust.method = "BH", sort = TRUE,
            alpha = 1)
  {

    qrL <- qr(contrasts)
    Lh <- qr.Q(qrL)[,1:qrL$rank]
    logFCs <- apply(contrasts, 2, function(models, L) sapply(models, getContrast, L = L), models = models)



    full_coef <- rownames(contrasts)
    full_varcov <- matrix(0, nrow = length(full_coef), ncol = length(full_coef),
                           dimnames = list(full_coef, full_coef)
    )
    Fstat <- sapply(models,
                    function(model, Lh, r){
                      if (getFitMethod(model)=="fitError")
                        return(NA) else {
                          cfs <- getCoef(model)
                          present_cfs <- intersect(names(cfs), full_coef)
                          varcov <- full_varcov
                          ## fill out missing coefficients to allow matrix multiplication
                          varcov[present_cfs,present_cfs] <- as.matrix(getVcovUnscaled(model))
                          cfs[setdiff(full_coef, names(cfs))] <- 0
                          cfs[full_coef]

                          est <- cfs %*% Lh
                          varEst <- t(Lh)%*%(varcov * getVarPosterior(model))%*%Lh
                          if (inherits(try(solve(varEst), silent = TRUE), "try-error")) return(NA)
                          Fstat <- (est %*% solve(varEst) %*% t(est) )/ r
                          return(Fstat)
                        }
                    },
                    Lh = Lh,
                    r = qrL$rank)
    df <- sapply(models, getDfPosterior)
    pval <- pf(Fstat, qrL$rank, df, lower.tail=FALSE) #pchisq(Fstat*qrL$rank, df = qrL$rank, lower.tail = F)
    adjPval <- p.adjust(pval, method = adjust.method)
    out <- data.frame(logFCs, df1 = qrL$rank, df2 = df, Fstat, pval, adjPval)
    if (sort) {
      if (alpha < 1) {
        ids <- adjPval < alpha
        out <- na.exclude(out[ids, ])
      }
      return(out[order(out$pval), ])
    }
    else {
      return(out)
    }
  }

# setGeneric("omnibusTest",function(object,...) standardGeneric("omnibusTest"))
#
# setMethod("omnibusTest","QFeatures",
#           function(object,
#                    i,
#                    contrasts,
#                    adjust.method="BH",
#                    modelColumn="msqrobModels",
#                    resultsColumnName="omnibusTest",
#                    overwrite=FALSE){
#             if (is.null(object[[i]])) stop(paste0("QFeatures object does not contain an assay with the name ",i))
#             if(!(modelColumn %in% colnames(rowData(object[[i]])))) stop(paste0("There is no column named \'", modelColumn,"\' with stored models of an msqrob fit in the rowData of assay ",i,"of the QFeatures object."))
#             if(is.null(colnames(contrasts)) & resultsColumnName=="") resultsColumnName<-"omnibusTest"
#             if((resultsColumnName %in% colnames(rowData(object[[i]])))&!overwrite) stop(paste0("There is/are already column(s) named \'", resultsColumnName, "\' in the rowData of assay ",i," of the QFeatures object, set the argument overwrite=TRUE to replace the column(s) with the new results or use another name for the argument resultsColumnName"))
#             rowData(object[[i]])[[resultsColumnName]] <- topFeaturesOmnibus(rowData(object[[i]])[,modelColumn], contrasts = contrasts, adjust.method = adjust.method, sort = FALSE, alpha = 1)
#             return(object)
#           })
