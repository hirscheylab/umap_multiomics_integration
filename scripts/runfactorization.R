##### runfactoization runs all the considered multi-omics factorization
### the required inputs are:
### "folder" corresponding to the path to the folder where the input files are contained, in the idea that all omics matrices are organized inside a unique folder
### "file.names" corresponding to a vector containing the names of all the omics files
### "num.factors" containing the number of factors in which we desire to decompose the matrices
### "sep=" "" corresponding to the separator used in the omics files required to properly read them
### "single.cell" indicating if the data are single cell data. In this case the filtering of the data will be more intense in respect to other data types
### "filtering"

### the input files need to be log2 transformed before running the analysis

library(tidyverse)
library(RGCCA)
library(IntNMF)
library(omicade4)
library(MOFA2)
library(r.jive)
library(tensorBSS)
source("tICA.R")
library(ubmi)
# remotes::install_version("Matrix", version = "1.6-1")

runfactorization <- function(folder, 
                             file.names, 
                             num.factors,
                             n_components_conc = 2,
                             min_pts = 5,
                             compute_features = TRUE,
                             sep = " ", 
                             filtering = "none",
                             methods = c("RGCCA", "MCIA", "MOFA", "intNMF", "JIVE", "tICA", "UBMI")
                             ) {
  factorizations <- list()
  t <- 1
  method <- numeric(0)
  
  num.factors <- as.numeric(num.factors)
  
  ##creating list of omics
  omics <- list()
  for (i in 1:length(file.names)){
    omics[[i]] <- as.matrix(read.table(paste(folder, file.names[i], sep = "/"), sep = sep, row.names = 1, header = T))
  }
  
  ##restricting to common samples and filtering
  samples <- colnames(omics[[1]])
  for (j in 1:length(omics)) {
    samples <- intersect(samples,colnames(omics[[j]]))
  }
  for (j in 1:length(omics)) {
    omics[[j]] <- omics[[j]][,samples]
    if (filtering != "none") {
      x <- apply(omics[[j]], 1, sd)
      x <- as.matrix(sort(x, decreasing = T))
      w <- which(x > 0)
      if (filtering == "stringent") {
        selected <- rownames(x)[1:min(w[length(w)],5000)]
      } else {
        selected <- rownames(x)[1:min(w[length(w)],6000)]
      }
      m <- match(rownames(omics[[j]]), selected)
      w <- which(!is.na(m))
      omics[[j]] <- omics[[j]][w,]
    } else {
      omics[[j]] <- omics[[j]][, which(apply(omics[[j]], 2, sd) > 0)]
    }
  }
  
  omics_pos <- list()
  for(j in 1:length(omics)) {
    if (min(omics[[j]]) < 0){
      omics_pos[[j]] <- omics[[j]] + abs(min(omics[[j]]))
    } else {
      omics_pos[[j]] <- omics[[j]]
    }
    omics_pos[[j]] <- omics_pos[[j]]/max(omics_pos[[j]])
  }
  
  if ("RGCCA" %in% methods) { ### RGCCA 
    factorizations_RGCCA <- rgcca(lapply(omics, function(x) t(x)), ncomp = rep(num.factors, length(omics)), 
                                  scheme = "centroid", scale = TRUE, init = "svd", bias = TRUE, tol = 1e-08, verbose = F)
    factors_rgcca <- as.matrix(factorizations_RGCCA$Y[[1]])
    metagenes_rgcca <- list()
    for(j in 1:length(omics)){
      metagenes_rgcca[[j]] <- as.matrix(factorizations_RGCCA$a[[j]])
      rownames(metagenes_rgcca[[j]]) <- rownames(omics[[j]])
      colnames(metagenes_rgcca[[j]]) <- 1:num.factors
    }
    
    factorizations[[t]] <- list(factors_rgcca, metagenes_rgcca)
    t <- t + 1
    method <- c(method, "RGCCA")
    print("Done RGCCA!") 
  }
  
  if ("MCIA" %in% methods) { ### MCIA
    factorizations_mcia <- mcia(omics_pos, cia.nf = num.factors)
    factors_mcia <- as.matrix(factorizations_mcia$mcoa$SynVar)
    metagenes_mcia <- list()
    for (j in 1:length(omics)) {
      metagenes_mcia[[j]] <- as.matrix(factorizations_mcia$mcoa$axis[1:dim(omics[[j]])[1],])
      rownames(metagenes_mcia[[j]]) <- rownames(omics[[j]])
      colnames(metagenes_mcia[[j]]) <- 1:num.factors
    }
    
    factorizations[[t]] <- list(factors_mcia, metagenes_mcia)
    t <- t + 1
    method <- c(method, "MCIA")
    print("Done MCIA!") 
  }
  
  if ("MOFA" %in% methods) { ### MOFA
    MOFAobject <- MOFA2::create_mofa(omics)
    DataOptions <- MOFA2::get_default_data_options(MOFAobject)
    ModelOptions <- MOFA2::get_default_model_options(MOFAobject)
    ModelOptions$num_factors <- num.factors
    TrainOptions <- MOFA2::get_default_training_options(MOFAobject)
    
    MOFAobject <- MOFA2::prepare_mofa(
      MOFAobject,
      data_options = DataOptions,
      model_options = ModelOptions,
      training_options = TrainOptions
    )
    
    MOFAobject <- MOFA2::run_mofa(MOFAobject, use_basilisk = TRUE)
    metagenes_mofa <- MOFA2::get_weights(MOFAobject)
    factors_mofa <- MOFA2::get_factors(MOFAobject)
    
    factorizations[[t]] <- list(factors_mofa$group1, metagenes_mofa)
    t <- t + 1
    method <- c(method, "MOFA")
    print("Done MOFA!") 
  }
  
  if ("intNMF" %in% methods) { ### intNMF
    factorizations_intnmf <- nmf.mnnals(dat = lapply(omics_pos, function(x) t(x)), k = num.factors)
    factors_intNMF <- as.matrix(factorizations_intnmf$W)
    colnames(factors_intNMF) <- 1:num.factors
    metagenes_intNMF <- list()
    for (j in 1:length(omics)) {
      metagenes_intNMF[[j]] <- t(factorizations_intnmf$H[[j]])
      rownames(metagenes_intNMF[[j]]) <- rownames(omics[[j]])
      colnames(metagenes_intNMF[[j]]) <- 1:num.factors
    }
    
    factorizations[[t]] <- list(factors_intNMF, metagenes_intNMF)
    t <- t + 1
    method <- c(method, "intNMF")
    print("Done intNMF!") 
  }
  
  if ("JIVE" %in% methods) { ### JIVE
    factorizations_jive <- jive(omics, rankJ = num.factors, rankA = rep(num.factors, length(omics)), 
                                method = "given", conv = "default", maxiter = 100, showProgress = FALSE)
    rankJV <- factorizations_jive$rankJ
    rankIV.v <- factorizations_jive$rankA
    J <- numeric(0)
    ng <- 0
    metagenes_jive <- list()
    for (j in 1:length(omics)) {
      J <- rbind(J,factorizations_jive$joint[[j]])
      ng <- c(ng, dim(factorizations_jive$joint[[j]])[1])
    }
    svd.o <- svd(J)
    jV <- svd.o$v %*% diag(svd.o$d)
    for (j in 1:length(omics)) {
      metagenes_jive[[j]] <- svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV] ###error in dimension
      rownames(metagenes_jive[[j]]) <- rownames(omics[[j]])
      colnames(metagenes_jive[[j]]) <- 1:num.factors
    }
    factors_jive <- jV[,1:rankJV]
    rownames(factors_jive) <- colnames(omics[[1]])
    colnames(factors_jive) <- 1:num.factors
    
    factorizations[[t]] <- list(factors_jive, metagenes_jive)
    t <- t + 1
    method <- c(method, "JIVE")
    print("Done JIVE!")
  }
  
  if ("tICA" %in% methods) { ### tICA
    omics_tensor <- list()
    for (j in 1:length(omics)) {
      omics_tensor[[j]] <- cor(omics[[j]], method = "spearman")
    }
    
    S <- vector(length = dim(omics[[1]])[2]*dim(omics[[1]])[2]*length(omics))
    dim(S) <- c(length(omics), dim(omics[[1]])[2], dim(omics[[1]])[2])
    for (j in 1:length(omics)) {
      S[j,,] <- t(omics_tensor[[j]])
    }
    tICA <- DoTICA(S, num.factors, method = "FOBI")
    factors_tica <- as.matrix(tICA$signals)
    rownames(factors_tica) <- colnames(omics[[1]])
    metagenes_tica <- list()
    for (j in 1:length(omics)) {
      metagenes_tica[[j]] <- omics[[j]] %*% ginv(t(tICA$signals))
      rownames(metagenes_tica[[j]]) <- rownames(omics[[j]])
      colnames(metagenes_tica[[j]]) <- 1:num.factors
    }
    
    factorizations[[t]] <- list(factors_tica, metagenes_tica)
    t <- t + 1
    method <- c(method, "tICA")
    print("Done tICA!")
  }
  
  if ("UBMI" %in% methods) { ### UBMI

    ubmi_object <- ubmi(omics, 
                        umap_params = list(n_components = num.factors),
                        umap_params_conc = list(n_components = n_components_conc),
                        compute_features = compute_features, 
                        samples_in_rows = FALSE, 
                        # xgboost_params = list(lambda = 0, eta = 0.8, gamma = 10),
                        min_pts = min_pts)

    factorizations[[t]] <- list(ubmi_object@factors[, colnames(ubmi_object@factors) != "clust"], 
                                ubmi_object@metagenes)
    t <- t + 1
    method <- c(method, "UBMI")
    print("Done UBMI!")
  }
  
  ## OUTPUT
  out <- list(factorizations = factorizations, method = method)
  
  if ("intNMF" %in% methods) {
    out <- append(out, list(intNMF.clusters = as.matrix(factorizations_intnmf$clusters)))
  }
  if ("UBMI" %in% methods) {
    out <- append(out, list(UBMI.clusters = as.matrix(ubmi_object@factors$clust)))
  }
  
  return(out)
}

