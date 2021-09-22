# This is the example of the simulation experiment with given parameters in the paper.

##################################################################
#
#                    Part I simulation
#
##################################################################

library("abind")
library("emdbook")
library("MASS")
library("DESeq2")
library("moments")
library("MAST")
library("lme4")
library("BSDE")

#senarios
nMean = 3 # of genes shows mean differences
nVar = 3  # of genes shows variance differences
nMult = 3 # of genes shows multimodality differences
nDP = 3   # of genes shows different proportion for mixture models
nGeneBlank = 18
nTotal = nMean + nVar + nMult + nDP + nGeneBlank
ncase = 20   # of case subjects
nctrl = 20   # of ctrl subjects
ncell = 100  # of cells per individual
nall = ncase + nctrl

#simulation of size factors
r_mean = 1.2
r_var = 1.2
r_mult = 0.6
r_dp = 0.2


i_mean = 1:nMean
i_var = (nMean + 1):(nMean + nVar)
i_mult = (nMean + nVar + 1):(nMean + nVar + nMult)
i_dp = (nMean + nVar + nMult + 1):(nMean + nVar + nMult + nDP)
i_blank = (nMean + nVar + nMult + nDP + 1):nTotal

#randomlized settings, switching cases and control samples to make sure the library sizes the same
case_modify_flag = rbinom(nTotal, 1, 0.5)
i_case_modify = which(case_modify_flag == 1)

#first case, then ctrl
i_case = 1:(ncase * ncell)
i_ctrl = (ncase * ncell + 1):(nall * ncell)

# ############ Simulation Data Preparation ###############################

file_tag = paste0(r_mean, "_", r_var, "_", r_mult, "_", r_dp)
dir = "../data/"
setwd(dir)

RIN_adj = TRUE

#readin data
t_mean = readRDS(paste0(dir, "autism/dca_PFCL23_100x10perc/mean.rds"))
t_disp = readRDS(paste0(dir, "autism/dca_PFCL23_100x10perc/dispersion.rds"))
t_drop = readRDS(paste0(dir, "autism/dca_PFCL23_100x10perc/dropout.rds"))
t_meta = readRDS(paste0(dir, "autism/meta_PFCL23_100x10perc.rds"))

dim(t_mean)
t_mean[1:2, 1:5]

dim(t_disp)
t_disp[1:2, 1:5]

dim(t_drop)
t_drop[1:2, 1:5]

table(t_meta$individual)
t_meta$RIN = t_meta$RNA.Integrity.Number

###########simulate basic parameters of ZINB model################
res = simu_base_param(t_mean, t_disp, t_drop, t_meta, nTotal, nall, RIN_adj =
                        RIN_adj)
sample_ctrl = res$sample_ctrl
RIN.simu = res$RIN.simu

random_idx_gene = sample.int(nTotal)
random_idx_sam  = sample.int(nall)

sample_ctrl = sample_ctrl[random_idx_gene, random_idx_sam,]

##########  calculate parameters for the case data #################

# sample gene index for genes differential expressed by mean or variance.
special_index = sample.int(nTotal, (nMean + nVar + nMult + nDP))
mean_index    = as.numeric(special_index[i_mean])
var_index     = as.numeric(special_index[i_var])
mult_index     = as.numeric(special_index[i_mult])
dp_index     = as.numeric(special_index[i_dp])
# label and save the DE index information.
de.mean = rep(0, nTotal)
de.var  = rep(0, nTotal)
de.mult  = rep(0, nTotal)
de.dp  = rep(0, nTotal)
de.mean[mean_index] = 1
de.var[var_index]   = 1
de.mult[mult_index]   = 1
de.dp[dp_index]   = 1
# To make sure all parameters are non-negative, we do some transformation
r_mean2 = r_mean
r_var2  = r_var

if (r_mean < 1) {
  r_mean2 = 1 / r_mean
}

if (r_var < 1) {
  r_var2 = 1 / r_var
}

# modify parameters for cases mean and var
# now r_m (r_mean2) is smaller than 1. We need to modify gene expression in
# x proportion of genes by fold r_m and (1-x) proportion of genes by
# fold of 1/r_m, so that x(1 - r_m) = (1-x)(1/r_m - 1)
# x (1 - r_m + 1/r_m -1) = (1/r_m -1) => x = 1/(1 + r_m)


sample_param_case = abind::abind(sample_ctrl[, 1:ncase,], sample_ctrl[, 1:ncase,], along =
                                   4)#gene x ind x (mean,overdispersion, droupout) x #of mixture models(2)

sample_param_ctrl = abind::abind(sample_ctrl[, (ncase + 1):nall,], sample_ctrl[, (ncase +
                                                                                    1):nall,], along = 4)


for (i in mean_index) {
  #distr1==distr2, case!=ctrl
  for (j in 1:ncase) {
    x = sample_param_case[i, j, , 1]
    sample_param_case[i, j, 1:2, 1] = calc_zinb_param(
      mu = x[1],
      theta = x[2],
      drop = x[3],
      r_m = r_mean2,
      r_v = 1
    )
    sample_param_case[i, j, 1:2, 2] = sample_param_case[i, j, 1:2, 1]
  }
}
for (i in var_index) {
  #distr1==distr2, case!=ctrl
  for (j in 1:ncase) {
    x = sample_param_case[i, j, , 1]
    sample_param_case[i, j, 1:2, 1] = calc_zinb_param(
      mu = x[1],
      theta = x[2],
      drop = x[3],
      r_m = 1,
      r_v = r_var2
    )
    sample_param_case[i, j, 1:2, 2] = sample_param_case[i, j, 1:2, 1]
  }
}

for (i in mult_index) {
  for (j in 1:ncase) {
    #modify mean
    sample_param_case[i, j, 1, 1] = sample_param_case[i, j, 1, 1] + sample_param_case[i, j, 1, 1] *
      (r_mult)
    sample_param_case[i, j, 1, 2] = sample_param_case[i, j, 1, 2] - sample_param_case[i, j, 1, 2] *
      (r_mult)
    #modify disp
    sample_param_case[i, j, 2, 1] = sample_param_case[i, j, 2, 1] * 2
    sample_param_case[i, j, 2, 2] = sample_param_case[i, j, 2, 1]
  }
  for (j in 1:nctrl) {
    sample_param_ctrl[i, j, 2, 1] = calc_multimod_param(size = (sample_param_ctrl[i, j, 2, 1] *
                                                                  2), t = r_mult)
    sample_param_ctrl[i, j, 2, 2] = sample_param_ctrl[i, j, 2, 1]
  }
}

for (i in dp_index) {
  for (j in 1:ncase) {
    x = sample_param_case[i, j, , 1]
    sample_param_case[i, j, 1:2, 1] = calc_zinb_param(
      mu = x[1],
      theta = x[2],
      drop = x[3],
      r_m = 4,
      r_v = 1
    ) #here we set r_m=4
  }
  for (j in 1:nctrl) {
    x = sample_param_ctrl[i, j, , 2]
    sample_param_ctrl[i, j, 1:2, 2] = calc_zinb_param(
      mu = x[1],
      theta = x[2],
      drop = x[3],
      r_m = 4,
      r_v = 1
    ) #here we set r_m=4
  }
}

temp = NA
temp = sample_param_case[i_case_modify, , ,]
sample_param_case[i_case_modify, , ,] = sample_param_ctrl[i_case_modify, , ,]
sample_param_ctrl[i_case_modify, , ,] = temp



# simulate scRNAseq based on zinb parameters of cases and controls:
sim_matrix = matrix(nrow = nTotal, ncol = nall * ncell)
sim_param = array(dim = c(nTotal, nall * ncell, 3))
date()
for (i in 1:nall) {
  idx_i = ((i - 1) * ncell + 1):(i * ncell)
  for (k in 1:ncell) {
    cell_distr_flag = rbinom(1, 1, 0.5) + 1 #indicate the minor prop for simulation
    sub_sample_param_case = sample_param_case[, , , cell_distr_flag] #equally separated for mult as design
    sub_sample_param_ctrl = sample_param_ctrl[, , , cell_distr_flag] #equally separated for mult as design
    cell_distr_flag_dp_add = (cell_distr_flag - 1) * rbinom(1, 1, prob =
                                                              (r_dp * 2)) + 1
    sub_sample_param_case[which(de.dp == 1), ,] = sample_param_case[which(de.dp ==
                                                                            1), , , cell_distr_flag_dp_add] #equally separated for dp as design
    sub_sample_param_ctrl[which(de.dp == 1), ,] = sample_param_ctrl[which(de.dp ==
                                                                            1), , , cell_distr_flag_dp_add] #equally separated for dp as design

    if (i > ncase) {
      mean_i = sub_sample_param_ctrl[, (i - ncase), 1]
      disp_i = sub_sample_param_ctrl[, (i - ncase), 2]
      drop_i = sub_sample_param_ctrl[, (i - ncase), 3]
      sample_mean_sd_i = sub_sample_param_ctrl[, (i - ncase), 4]
    } else{
      mean_i = sub_sample_param_case[, i, 1]
      disp_i = sub_sample_param_case[, i, 2]
      drop_i = sub_sample_param_case[, i, 3]
      sample_mean_sd_i = sub_sample_param_case[, i, 4]
    }

    sample_mean_k = exp(rnorm(nTotal, log(mean_i), sample_mean_sd_i))
    for (ig in 1:nTotal) {
      cur_count = emdbook::rzinbinom(1, sample_mean_k[ig], disp_i[ig], drop_i[ig])
      if (is.na(cur_count)) {
        print(c(i, ik, ig, sample_mean_k[ig], disp_i[ig], drop_i[ig]))
        while (is.na(cur_count)) {
          cur_count = emdbook::rzinbinom(1, sample_mean_k[ig], disp_i[ig], drop_i[ig])
          if (!is.na(cur_count)) {
            print(c(i, ik, ig, cur_count))
          }
        }
      }
      sim_matrix[ig, idx_i[k]] = cur_count
      sim_param[ig, idx_i[k],] = c(sample_mean_k[ig], disp_i[ig], drop_i[ig])
    }
  }
  #print(i)
}
date()

dim(sim_matrix)
sim_matrix[1:4, 1:4]

table(c(sim_matrix) == 0)
table(c(sim_matrix) == 0) / (nrow(sim_matrix) * ncol(sim_matrix))



# scatter plot 2:check simulation
par(mfrow = c(2, 5))

for (ig in which(de.mean + de.var + de.mult + de.dp == 0)[1:2]) {
  for (i in c(1:5, (nall - 4):nall)) {
    idx_i = ((i - 1) * ncell + 1):(i * ncell)
    if (i <= ncase) {
      cur_col = "red"
      cur_label = paste0("non-DE gene ", ig, ", case ")
    }
    if (i >= ncase) {
      cur_col = "blue"
      cur_label = paste0("non-DE gene ", ig, ", ctrl ")
    }
    hist(sim_matrix[ig, idx_i],
         col = cur_col,
         main = cur_label,
         xlab = "count")
  }
}

for (ig in which(de.mean == 1)[1:2]) {
  for (i in c(1:5, (nall - 4):nall)) {
    idx_i = ((i - 1) * ncell + 1):(i * ncell)
    if (i <= ncase) {
      cur_col = "red"
      cur_label = paste0("Mean-DE gene ", ig, ", case ")
    }
    if (i >= ncase) {
      cur_col = "blue"
      cur_label = paste0("Mean-DE gene ", ig, ", ctrl ")
    }
    hist(sim_matrix[ig, idx_i],
         col = cur_col,
         main = cur_label,
         xlab = "count")
  }
}

for (ig in which(de.var == 1)[1:2]) {
  for (i in c(1:5, (nall - 4):nall)) {
    idx_i = ((i - 1) * ncell + 1):(i * ncell)
    if (i <= ncase) {
      cur_col = "red"
      cur_label = paste0("Var-DE gene ", ig, ", case ")
    }
    if (i >= ncase) {
      cur_col = "blue"
      cur_label = paste0("Var-DE gene ", ig, ", ctrl ")
    }
    hist(sim_matrix[ig, idx_i],
         col = cur_col,
         main = cur_label,
         xlab = "count")
  }
}

for (ig in which(de.mult == 1)[1:2]) {
  for (i in c(1:5, (nall - 4):nall)) {
    idx_i = ((i - 1) * ncell + 1):(i * ncell)
    if (i <= ncase) {
      cur_col = "red"
      cur_label = paste0("Mult-DE gene ", ig, ", case ")
    }
    if (i >= ncase) {
      cur_col = "blue"
      cur_label = paste0("Mult-DE gene ", ig, ", ctrl ")
    }
    hist(sim_matrix[ig, idx_i],
         col = cur_col,
         main = cur_label,
         xlab = "count")
  }
}

for (ig in which(de.dp == 1)[1:2]) {
  for (i in c(1:5, (nall - 4):nall)) {
    idx_i = ((i - 1) * ncell + 1):(i * ncell)
    if (i <= ncase) {
      cur_col = "red"
      cur_label = paste0("DP-DE gene ", ig, ", case ")
    }
    if (i >= ncase) {
      cur_col = "blue"
      cur_label = paste0("DP-DE gene ", ig, ", ctrl ")
    }
    hist(sim_matrix[ig, idx_i],
         col = cur_col,
         main = cur_label,
         xlab = "count")
  }
}


####################### Meta information collection #################

#the phenotype and individual information of simulated samples.
phenotype = c(rep(0, nctrl * ncell), rep(1, ncase * ncell))
individual = paste0("ind", c(rep(1:nall, each = ncell)))

#Count info for matrix
cell_id = paste0("cell", 1:(nall * ncell))
gene_id = paste0("gene", 1:nTotal)

rownames(sim_matrix) = gene_id
colnames(sim_matrix) = cell_id

dimnames(sim_param) = list(gene_id, cell_id, c("mean", "overdisp", "dropout"))

#Cell info for meta
cellsum = apply(sim_matrix, 2, sum)
genesum = apply(sim_matrix, 1, sum)
CDR  = apply(sim_matrix > 0, 2, sum) / nrow(sim_matrix)
meta = data.frame(cell_id,
                  individual,
                  phenotype,
                  cellsum,
                  CDR,
                  stringsAsFactors = FALSE)

meta$RIN = rep(RIN.simu, each = ncell) #if RNA_adj==FALSE in simu_base_param, meta$RIN=NA

dim(meta)
meta[1:2,]

######### Some basic stat & Preparation for Bulk RNAseq analysis ############
library(DESeq2)
# individual level info
cur_individual = unique(meta$individual)
cell_num = matrix(ncol = 1, nrow = length(cur_individual))
rownames(cell_num) = cur_individual
colnames(cell_num) = "cell_num"

read_depth = matrix(ncol = 1, nrow = length(cur_individual))
rownames(read_depth) = cur_individual
colnames(read_depth) = "read_depth"

phenotype_ind = matrix(ncol = 1, nrow = length(cur_individual))
rownames(phenotype_ind) = cur_individual
colnames(phenotype_ind) = "phenotype"

zero_rate_ind = matrix(nrow = nrow(sim_matrix),
                       ncol = length(cur_individual))
rownames(zero_rate_ind) = rownames(sim_matrix)
colnames(zero_rate_ind) = cur_individual

sim_matrix_bulk = matrix(nrow = nrow(sim_matrix),
                         ncol = length(cur_individual))
rownames(sim_matrix_bulk) = rownames(sim_matrix)
colnames(sim_matrix_bulk) = cur_individual

for (i_ind in 1:length(cur_individual)) {
  cur_ind = cur_individual[i_ind]
  #fit org
  cur_ind_m = sim_matrix[, meta$individual == cur_ind]
  cell_num[i_ind]   = ncol(cur_ind_m)
  read_depth[i_ind] = sum(cur_ind_m, na.rm = TRUE) / cell_num[i_ind] * 1000
  phenotype_ind[i_ind] = meta$phenotype[meta$individual == cur_ind][1]

  zero_rate_ind[, i_ind] = rowSums(cur_ind_m == 0, na.rm = TRUE) / cell_num[i_ind]
  sim_matrix_bulk[, i_ind] = rowSums(cur_ind_m, na.rm = TRUE)
}

tapply(read_depth, phenotype_ind, summary)

########################Other cell level info################################
library("MAST")
library("lme4")
cell_readdepth = apply(sim_matrix, 2, sum)
cell_count_quantile = apply(sim_matrix, 2, quantile)

#almost no need to adjust library size.
# saveRDS(de.mean,
#         paste0(dir, "simulated/sim_de.mean_", file_tag, ".rds"))
# saveRDS(de.var, paste0(dir, "simulated/sim_de.var_", file_tag, ".rds"))
# saveRDS(de.mult,
#         paste0(dir, "simulated/sim_de.mult_", file_tag, ".rds"))
# saveRDS(de.dp, paste0(dir, "simulated/sim_de.dp_", file_tag, ".rds"))
#
# saveRDS(read_depth,
#         paste0(dir, "simulated/sim_ind_readdepth_", file_tag, ".rds"))
# saveRDS(phenotype_ind,
#         paste0(dir, "simulated/sim_ind_phenotye_", file_tag, ".rds"))
# saveRDS(zero_rate_ind,
#         paste0(dir, "simulated/sim_ind_zero_rate_", file_tag, ".rds"))
#
# saveRDS(sim_param,
#         paste0(dir, "simulated/sim_param_", file_tag, ".rds"))
# saveRDS(sim_matrix,
#         paste0(dir, "simulated/sim_matrix_", file_tag, ".rds"))
# saveRDS(meta, paste0(dir, "simulated/sim_meta_", file_tag, ".rds"))
# saveRDS(sim_matrix_bulk,
#         paste0(dir, "simulated/sim_matrix_bulk_", file_tag, ".rds"))
#
# saveRDS(cell_readdepth,
#         paste0(dir, "simulated/sim_cell_readdepth_", file_tag, ".rds"))
# saveRDS(
#   cell_count_quantile,
#   paste0(dir, "simulated/sim_cell_count_quantile_", file_tag, ".rds")
# )
#
# write.csv(sim_matrix,
#           paste0(dir, "simulated/sim_matrix_", file_tag, ".csv"))


##################################################################
#
#                    Part II simulation
#
##################################################################

###DESeq2######

# We calculate bulk information by summing up raw counts of
# all cells(of certain cluster) of an individual within a genes

cur_info = meta[, c("individual", "phenotype")]
cur_info = unique(cur_info)
rownames(cur_info) = cur_info$individual
cur_info$phenotype = as.factor(cur_info$phenotype)

# object construction
dds = DESeqDataSetFromMatrix(countData = sim_matrix_bulk,
                             colData = cur_info,
                             design = ~ phenotype)

# observed pvalue calculation
dds = DESeq(dds)
deseq_pval = results(dds)$pvalue


#####MAST######

#we take the method of current usage of MAST,
#this code care with the results of simulation data, and calculate the MAST


# the input of MAST analysis can be many format including matrix and
# SingleCellAssay.

# input:
# (1)the log-transformed expression count matrix,
#     with each column represents the cells and each rows represents the genes.
# (2)the meta data, including cell information and individual information.

# we get the p-values based on the Hurdle model ("H" model)

sim_matrix_log = log2(1 + sim_matrix) #log transformed data

dim(sim_matrix_log)
sim_matrix_log[1:2, 1:5]

cell_id = colnames(sim_matrix_log)   #get the cell id from the data
gene_id = rownames(sim_matrix_log)   #get the gene id from the data

fData = data.frame(primerid = gene_id)
cData = data.frame(wellKey = cell_id)

diagnosis = as.character(meta$phenotype) #

diagnosis2 = matrix("Control", ncol = 1, nrow = length(diagnosis))
diagnosis2[which(diagnosis == "1")] = "Case"
diagnosis = as.factor(diagnosis2)

sca = MAST::FromMatrix(sim_matrix_log, cData, fData)
colData(sca)$cngeneson = as.numeric(meta$CDR)
colData(sca)$diagnosis = as.factor(diagnosis)
colData(sca)$ind = as.factor(meta$individual)

colData(sca)

print(date())
print(gc())
b1 = tryCatch(
  MAST::zlm(
    formula = ~ diagnosis + (1 |
                               ind),
    sca = sca,
    method = 'glmer',
    ebayes = FALSE,
    parallel = TRUE
  ),
  error = function(e) {
    NA
  }
)
print(date())
print(gc())

lrt1 = tryCatch(
  MAST::lrTest(b1, "diagnosis"),
  error = function(e) {
    NA
  }
)
MAST_pval1 = tryCatch(
  apply(lrt1, 1, function(x) {
    x[3, 3]
  }),
  error = function(e) {
    NA
  }
)
print(date())
print(gc())

####### BSDE ##########

perm_num = 500
lambda_fold = 1.5
unif_round_unit = 0.05
library("reticulate")
library("doRNG")
library("doParallel")
library("ggplot2")
library("emdbook")
library("Rcpp")
library("Barycenter")
library("BSDE")

op_pval = matrix(ncol = 1, nrow = nrow(sim_matrix_log))


print(date())
print(gc())

for (i_g in 1:nrow(sim_matrix_log)) {
  cur_sim = sim_matrix_log[i_g, ]
  cur_ind = meta$individual
  cur_pheno = phenotype

  op_pval[i_g] = tryCatch({
    cal_w2_pval(
      count_per_gene = cur_sim,
      meta_individual = cur_ind,
      meta_phenotype = cur_pheno,
      perm_num = perm_num,
      unif_round_unit = 0.5
    )[[1]]
  }, error = function(e) {
    NA
  })
}
print(date())
print(gc())
rownames(op_pval) = gene_id

##################################################################
#
#              Part III simulation result analysis
#
##################################################################
