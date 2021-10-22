#' Given a ZINB(\eqn{\mu}, \eqn{\theta}, \eqn{z}), returns ZINB(\eqn{\mu'}, \eqn{\theta'}, \eqn{z})
#' such that the mean changes by a factor of \eqn{r_m} and the variance changes
#' by a factor of \eqn{v_m}.
#' @param mu parameter \eqn{\mu} of ZINB distribution.
#' @param theta parameter \eqn{\theta} of ZINB distribution.
#' @param drop dropout rate \eqn{z} of ZINB distribution.
#' @param r_m targeted fold change of mean
#' @param r_v targeted fold change of variance
#'
#' @return \eqn{\mu'} and \eqn{\theta'} of the modified ZINB distribution
#'
#' @note It is required that \eqn{r_m / r_v < 1 + \mu / \theta}.
#' @export
#' @examples
#' calc_zinb_param(2,3,r_m=1.3)
#' calc_zinb_param(2,3,r_v=1.5)
#'
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


#' For a given mixture \eqn{0.5 ZINB(\mu(1-t), \theta, z) + 0.5 ZINB(\mu(1+t), \theta, z)},
#' return \eqn{ZINB(\mu, \theta', z)} such that both mean and variance match.
#'
#' @param size parameter \eqn{\theta} of the ZINB model.
#' @param t Parameter \eqn{t} (between 0 and 1).
#' @return Parameter \eqn{\theta'} of the output ZINB model.
#' @export
#'
calc_multimod_param = function(size, t = 0.5) {
  size3 = size / (size * t ^ 2 + t ^ 2 + 1)
  size3
}

#' Given parameters of a reference dataset, calculate the ZINB parameters for
#' simulating a new dataset of a different size.
#' @description The basic parameters of ZINB(\eqn{\mu}, \eqn{\theta}, \eqn{z}) from the reference dataset can be
#' estimated with DCA (Deep Count Autoencoder).
#' @param t_mean mean parameter matrix \eqn{\mu} of the reference dataset (gene x cell).
#' @param t_disp dispersion parameter matrix \eqn{\theta} of the reference dataset (gene x cell).
#' @param t_drop dropout parameter matrix \eqn{z} of the reference dataset (gene x cell).
#' @param t_meta Data frame (nrow = cell) of meta information.
#'   It must contain a column \code{individual} recording from which subject each cell is collected.
#'   If \code{RIN_adj=TRUE}, then it must also contain column \code{RIN} for RNA integrity number.
#' @param nTotal numbers of genes to simulate
#' @param RIN_adj If \code{TRUE}, then \code{t_meta} must contain column \code{RIN} for RNA integrity number.
#'   (default: FALSE)
#' @param nall numbers of individuals to simulate
#' @return a list containing:
#'  \itemize{
#'  \item \code{sample_ctrl}: gene x individual x param array,
#'    which gives the parameter \eqn{\mu}, \eqn{\theta} and \eqn{z} of ZINB for simulation.
#'  \item \code{RNA.simu}: residual standard deviations of means.
#'  }
#' @export
#' @references Eraslan, G., Simon, L. M., Mircea, M., Mueller, N. S., & Theis, F. J. (2019). Single-cell RNA-seq denoising using a deep count Autoencoder. Nature Communications, 10(1), 1-14.
#'
simu_base_param=function(t_mean, t_disp, t_drop, t_meta, nTotal=30, nall=40, RIN_adj=FALSE){
  #cell mean sum adj
  cell_mean_sum = colSums(t_mean)
  summary(cell_mean_sum)
  t_mean_adj = t(t(t_mean) * (10000 / cell_mean_sum))
  sample_ids=t_meta$individual

  tapply_mean <-
    function(x) {
      tapply(x, sample_ids, function(x)
        return(mean(x, na.rm = TRUE)))
    }
  tapply_median <-
    function(x) {
      tapply(x, sample_ids, function(x)
        return(stats::median(x, na.rm = TRUE)))
    }
  tapply_sd <-
    function(x) {
      tapply(x, sample_ids, function(x)
        return(stats::sd(x, na.rm = TRUE)))
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
  xmt = stats::model.matrix( ~ scale(log(cell_mean_sum)))
  dim(xmt)
  xmt[1:2, ]
  log_t_mean_adj   = log(t_mean_adj)
  log_t_mean_resid = matrix(NA, nrow = nrow(t_mean), ncol = ncol(t_mean))

  coef = matrix(NA, nrow = nrow(log_t_mean_adj), ncol = 2)
  for (i in 1:nrow(log_t_mean_adj)) {
    yi = log_t_mean_adj[i, ]
    li = stats::lm.fit(x = xmt, y = yi)
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
        cov_matrix[i, j] = stats::cov(sample_data[, i], sample_data[, j])
        cor_test[i, j] = stats::cor.test(sample_data[, i], sample_data[, j])$p.value
      }
    }



    #adjust for RIN
    RIN.simu=NA
    if(RIN_adj==TRUE){
      log_mean_ig = sample_data2[, 1]
      lmi  = stats::lm(log_mean_ig ~ sample_meta$RIN)
      beta = lmi$coefficients
      # add some extra variance for the mean parameter
      e1 = stats::rnorm(nall, mean = 0, sd = sqrt(cov_matrix[1, 1]))
      RIN.simu = stats::rnorm(nall)
      log_mean_ig_simu = beta[1] + beta[2] * RIN.simu + e1
      for (j in 1:nall) {
        sample_data_mean_j    = sample_data_mean
        sample_data_mean_j[1] = log_mean_ig_simu[j]

        sample_ctrl[ig, j,]  = exp(MASS::mvrnorm(1, mu = sample_data_mean_j, Sigma = cov_matrix))
      }
    }
    if(RIN_adj==FALSE){
      sample_ctrl[ig,,]=exp(MASS::mvrnorm(nall, mu = sample_data_mean, Sigma = cov_matrix,empirical = TRUE))
    }
  }

  sample_ctrl[, , "dropout"] = sample_ctrl[, , "dropout"] / (1 + sample_ctrl[, , "dropout"]) #logit trans
  sample_ctrl[, , "resid_mean_sd"] = exp(sample_ctrl[, , "resid_mean_sd"]) #double-log trans

  return(list(sample_ctrl=sample_ctrl,RIN.simu=RIN.simu))
}
