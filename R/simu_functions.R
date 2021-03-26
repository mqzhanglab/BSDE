#' this function returns the parameters for changing mean and variance of the zero-inflated negative binomial(ZINB) distribution
#' @param mu the original parameter mu of ZINB distribtuion
#' @param theta the original parameter theta of ZINB distribtuion, the distribution will be more close to zero-inflated poission with larger parameter size.
#' @param drop the original parameter dropout_rate of of ZINB distribution
#' @param r_m the targeted fold change of mean
#' @param v_m the targeted fold change of variance

#' @return two number for parameter mu and theta of modified ZINB distribution

#' @note For numerical meaningful calculation, is is required that r_m/r_v <1+mu/theta

#' @example
#' calc_zinb_param(2,3,r_m=1.3)
#' calc_zinb_param(2,3,r_v=1.5)

calc_zinb_param = function(mu,
                           theta,
                           drop = 0,
                           r_m = 1,
                           r_v = 1) {
  #mu=mean
  #theta=overdispersion
  mu2 = r_m * mu
  theta2 = theta * mu * r_m / (mu * r_v * r_m + (r_v * r_m - 1) * theta +
                                 (r_v - 1) * r_m * mu * drop * theta)
  if (theta2 < 0) {
    stop("negative theta2\n")
  }
  return(c(mu2, theta2))
}

#' This function returns the overdispersion parameters(size3) for the single ZINB distribution ZINB(mu_3, size_3, pi).
#' This distribution has the same dropout_rate as the mixture of two ZINB with equal 50%-50% proportions
#' The two ZINB have the same dropout_rate, same overdispersion (size) and the same distance of expectation from the mu_3,
#' i.e. ZINB(mu_3(1+t),size,pi_3) and ZINB(mu_3(1-t),size,pi_3)

#' @param size parameter theta of the zinb model, the distribution will be more close to zero-inflated poission with larger parameter size.
#' @param parameter t as the numeric 0<=t<=1
#' @return size3 The parameter theta in the notation
#' @example calc_multimod_param(25,t=0.3)
calc_multimod_param = function(size, t = 0.5) {
  size3 = size / (size * t ^ 2 + t ^ 2 + 1)
  size3
}

#' This function works for the basic parameter generation of ZINB models from the given reference data.
#' @discription The reference data is as the format of the output of DCA.
#'  The reference data includes the gene x cell matrices of the mean, dipersion and dropout rate
#' The reference data also includes a meta parameter, which is the dataframe for cell features.
#' The feature named as "individual" is necessary for individual information.
#' The the feature "RIN" is also necessary if parameter RIN_adj=TRUE.
#' @param t_mean a genexcell matrix, represent the mu parameter of ZINB models from reference data. DCA can generate it.
#' @param t_disp a genexcell matrix, represent the overdispersion parameter of ZINB models from reference data. DCA can generate it.
#' @param t_drop a genexcell matrix, represent the dropout_rate parameter of ZINB models from reference data. DCA can generate it.
#' @param nTotal numbers of genes needs to be simulated
#' @param nTotal numbers of subjects needs to be simulated
#' @return a list includes two objects sample_ctrl and RNA.simu
#' @return sample_ctrl is a gene x individual x param array, which gives the parameter "mean", "dispersion", and "dropout" for simulate a ZINB distribution.
#' It also gives the residual standard deviations of means.
#'
#'
#' @references Eraslan, G., Simon, L. M., Mircea, M., Mueller, N. S., & Theis, F. J. (2019). Single-cell RNA-seq denoising using a deep count autoencoder. Nature communications, 10(1), 1-14.
#'
simu_base_param=function(t_mean,t_disp,t_drop,t_meta,nTotal=30,nall=40, RIN_adj=TRUE){
  #cell mean sum adj
  cell_mean_sum = colSums(t_mean)
  summary(cell_mean_sum)
  t_mean_adj = t(t(t_mean) * (10000 / cell_mean_sum))
  ############### collect sample information ###################

  col_info   = strsplit(colnames(t_mean), split = "_")
  sample_ids = sapply(col_info, function(x)
    x[2])
  table(sample_ids)

  tapply_mean <-
    function(x) {
      tapply(x, sample_ids, function(x)
        return(mean(x, na.rm = TRUE)))
    }
  tapply_median <-
    function(x) {
      tapply(x, sample_ids, function(x)
        return(median(x, na.rm = TRUE)))
    }
  tapply_sd <-
    function(x) {
      tapply(x, sample_ids, function(x)
        return(sd(x, na.rm = TRUE)))
    }
  tapply_1 <- function(x) {
    tapply(x, sample_ids, function(x)
      x[1])
  }

  mid_mean = t(apply(t_mean, 1, tapply_median))
  mid_disp = t(apply(t_disp, 1, tapply_median))
  mid_drop = t(apply(t_drop, 1, tapply_median))
  mid_mean_adj = t(apply(t_mean_adj, 1, tapply_median))
  sample_log_mean_sd = t(apply(t_mean, 1, function(x)
    tapply_sd(log(as.numeric(
      x
    )))))
  sample_meta = t_meta[match(colnames(mid_mean), t_meta$individual), ]
  mid_drop[mid_drop == 0] = min(mid_drop[mid_drop > 0]) / 2 #/sum(mid_drop==0)#adjust for the situation that mid-drop==0

  sample_log_mean = log(mid_mean)
  sample_log_disp = log(mid_disp)
  sample_logit_drop = log((mid_drop) / (1 - mid_drop))
  sample_log_mean_adj = log(mid_mean_adj)

  # log residual sd
  xmt = model.matrix( ~ scale(log(cell_mean_sum)))
  dim(xmt)
  xmt[1:2, ]
  log_t_mean_adj   = log(t_mean_adj)
  log_t_mean_resid = matrix(NA, nrow = nrow(t_mean), ncol = ncol(t_mean))

  coef = matrix(NA, nrow = nrow(log_t_mean_adj), ncol = 2)
  for (i in 1:nrow(log_t_mean_adj)) {
    yi = log_t_mean_adj[i, ]
    li = lm.fit(x = xmt, y = yi)
    coef[i, ] = li$coefficients
    log_t_mean_resid[i, ] = li$residuals
  }
  summary(coef)
  sample_logresid_mean_sd = t(apply(log_t_mean_resid, 1, function(x)
    tapply_sd(x)))
  sample_loglogresid_mean_sd = log(sample_logresid_mean_sd)

  # log transformed cell readdepth
  log_cell_mean_sum = log(cell_mean_sum)
  sample_log_cell_mean_sum_sd = tapply_sd(log_cell_mean_sum)
  sample_log_cell_mean_sum_mean = tapply_mean(log_cell_mean_sum)

  sample_meta = cbind(sample_meta,
                      sample_log_cell_mean_sum_mean,
                      sample_log_cell_mean_sum_sd)

  dim(sample_log_mean)
  dim(sample_log_disp)
  dim(sample_logit_drop)
  dim(sample_log_mean_sd)
  dim(sample_loglogresid_mean_sd)
  dim(sample_log_mean_adj)

  sample_log_mean[1:2, 1:5]
  sample_log_mean_sd[1:2, 1:5]
  sample_loglogresid_mean_sd[1:2, 1:5]

  sample_log_cell_readdepth_mean = sample_meta$sample_log_cell_mean_sum_mean
  sample_log_cell_readdepth_sd = sample_meta$sample_log_cell_mean_sum_sd


  #characterize the parameters
  sample_ctrl = array(dim = c(nTotal, nall, 4),
                      dimnames = list(
                        paste0("gene", 1:nTotal),
                        paste0("ind", 1:nall),
                        c("mean", "dispersion", "dropout", "resid_mean_sd")
                      ))

  ind_strength = 0.5 #from 0 to 1, use this to adjust individual mean expression strength in the simulation.
  for (ig in 1:nTotal) {
    #sample_data=cbind(c(sample_log_mean[ig,]),c(sample_log_disp[ig,]),c(sample_logit_drop[ig,]),c(sample_log_mean_sd[ig,]))
    sample_data = cbind(
      c(sample_log_mean_adj[ig, ]),
      c(sample_log_disp[ig, ]),
      c(sample_logit_drop[ig, ]),
      c(sample_loglogresid_mean_sd[ig, ])
    ) #version2 with adjusted residuals

    sample_data_reg = apply(sample_data, 2, mean)

    sample_data2 = sample_data * ind_strength + matrix(rep(sample_data_reg, each =
                                                             nrow(sample_data)), nrow = nrow(sample_data)) * (1 - ind_strength)
    sample_data_mean = apply(sample_data2, 2, function(x)
      mean(x, na.rm = TRUE))

    cov_matrix = matrix(NA,
                        nrow = ncol(sample_data),
                        ncol = ncol(sample_data))
    cor_test = cov_matrix
    for (i in 1:ncol(sample_data)) {
      for (j in 1:ncol(sample_data)) {
        cov_matrix[i, j] = cov(sample_data[, i], sample_data[, j])
        cor_test[i, j] = cor.test(sample_data[, i], sample_data[, j])$p.value
      }
    }



    #adjust for RIN
    RIN.simu=NA
    if(RIN_adj==TRUE){
      log_mean_ig = sample_data2[, 1]
      lmi  = lm(log_mean_ig ~ sample_meta$RIN)
      beta = lmi$coefficients
      # add some extra variance for the mean parameter
      e1 = rnorm(nall, mean = 0, sd = sqrt(cov_matrix[1, 1]))
      RIN.simu = rnorm(nall)
      log_mean_ig_simu = beta[1] + beta[2] * RIN.simu + e1
      for (j in 1:nall) {
        sample_data_mean_j    = sample_data_mean
        sample_data_mean_j[1] = log_mean_ig_simu[j]

        sample_ctrl[ig, j,]  = exp(MASS::mvrnorm(1, mu = sample_data_mean_j, Sigma = cov_matrix))
      }
    }
    if(RIN_adj==FALSE){
      sample_ctrl[ig,,]=exp(MASS::mvrnorm(nall, mu = sample_data_mean_j, Sigma = cov_matrix,empirical = TRUE))
    }
  }

  sample_ctrl[, , "dropout"] = sample_ctrl[, , "dropout"] / (1 + sample_ctrl[, , "dropout"]) #logit trans
  sample_ctrl[, , "resid_mean_sd"] = exp(sample_ctrl[, , "resid_mean_sd"]) #double-log trans

  return(list(sample_ctrl=sample_ctrl,RIN.simu=RIN.simu))
}
