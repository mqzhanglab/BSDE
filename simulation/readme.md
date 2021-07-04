
# Introduction

This markdown is for individual case-control study based on Barycenter and Wasserstein distance estimation on the gene differential expression (DE) analysis of certain cell type based on simulated scRNAseq data from zero-inflated negative binomial (ZINB) distribution and compared the result with the [MAST](https://www.bioconductor.org/packages/release/bioc/html/MAST.html) packages(for scRNAseq) and the [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) packages(for bulk RNAseq)


# Initial Setup

R packages required in this pipeline:
```{r load_libraries,eval=runcode}
library("MASS")
library("emdbook")
library("moments")
library("MAST")
library("lme4")
library("DESeq2")
library("abind")

source("~/Dropbox/Barycenter/BSDEpackage/BSDE/R/simu_functions.R")
```

# Simulation

Here we present a demo version of our simulation here, which are minimized and based on the DCA output of the 10% autism dataset with top 100 expressed genes.

In this simulation example, we will generate 30 genes from a particular type of cell of 20 case subjects and 20 control subjects, with each subjects have 100 cells. We simulate the basic parameters, i.e., the mean, dispersion and dropout parameter from the distribution of the real reference scRNAseq database from paper [Single-cell genomics identifies cell type–specific molecular changes in autism](https://science.sciencemag.org/content/364/6441/685.abstract).


Within this 30 genes, we simulate 4 different type of DEs by adding on particular ratios to the basic parameters to get expected fold changes.

(1). mean DE: 10% genes are DE in mean but not in var, the genes are differential expressed in mean between cases and controls. We implement this setting by simply change the parameter /mu and /dispersion to suit the expected fold changes in mean with the function calc_zinb_param.

(2).var DE: 10% genes are DE in var but not mean, the genes are differential expressed in mean between cases and controls. In practice, we simulate 150 genes changes in cases, and another 150 genes changes in control. We implement this setting by simply change the parameter /mu and /dispersion to suit the expected fold changes in mean with the function calc_zinb_param.

(3).multimodality DE: 10% genes are DE due to Multimodality. Suppose Z have equaliy probability to falls into two NB/ZINB distribution Z1 from ZINB(mu1,size1,drop1), Z2 from ZINB(mu2,size2,drop2). Suppose another variable Z3 from ZINB(mu3,size3,drop3), who has the same mean and variance as Z1. Then it is not too difficult to calculate the mu3,size3,drop3 with known m1, size1, drop1 and m2, size2, drop2. We implement this setting by supporting the m1=m3+t, m2=m3-t, while the t equals to given proportion of m3(eg: t=0.2*m3). To avoid generating negative parameters. We will set the basic parameters as the parameters as the distribution of Z2, and then generate Z1 and Z3 according to t.

(4).disp DE: 10% genes are DE in dispersion.As the one of the parameter to describe ZINB model, dispersion majorly describes the shape of ZINB distribution. According to bulk gene analysis [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), the dispersion is related to mean when the expression level is not very high. The size factor is varied from r_p = 0.6, 0.7, 0.8, 0.9.

## Basic Settings

Senario settings in our analysis are as follows:


```{r basic_settings,eval=runcode}
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
```

## Simulation References Data Preparation

The reference data for simulation is from the paper [Single-cell genomics identifies cell type–specific molecular changes in autism](https://science.sciencemag.org/content/364/6441/685.abstract). Some details of their analysis is provided [here](https://science.sciencemag.org/content/suppl/2019/05/15/364.6441.685.DC1). Specifically, we do some basic data cleaning and processing. Then we take the first 3K genes who expressed in the most cells, and we randomly take 10 percentage of the single cells among all cells in the dataset. After this step, the data was denoised with [DCA](https://github.com/theislab/dca) package, a software to denoise the scRNAseq counts with autoencoder, an artificial neural network used to learn efficient data codings in an unsupervised manner.The input rawM3k10.csv includes a CSV/TSV-formatted raw count matrix with genes in rows and cells in columns.

```{bash dca_comments,eval=FALSE}
dca input.csv res_dca_input/
```

The output of DCA is an estimation of the expression of each gene and each cell in an ZINB distribution, which includes the estimation of the mean, dropout probabilities and dispersion for each cell and gene as an the ZINB distribution. In our simulation, we will borrow the parameters from the ouput of DCA.

Note: DCA itself includes the size factor of each cells in its pipelines.Thus the input of DCA should be the raw count data without any normalization or log transformation.

```{r load_expression_data,eval=runcode}
file_tag = paste0(r_mean, "_", r_var, "_", r_mult, "_", r_dp)
dir = "~/Dropbox/Barycenter/BSDEpackage/BSDE/data/"
setwd(dir)

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
t_meta$RIN=t_meta$RNA.Integrity.Number

```

## Simulate Parameters of ZINB Model
We use the function simu_base_param to generate the basic parameters of ZINB models for simulation.The output of simu_base_param, sample_ctrl, is a gene x individual x param array, which gives the parameter "mean", "dispersion", and "dropout" for simulate a ZINB distribution. It also gives the residual standard deviations of means.

```{r zinb_param,eval=runcode}
RIN_adj=TRUE
res=simu_base_param(t_mean,t_disp,t_drop,t_meta,nTotal,nall, RIN_adj=RIN_adj)
sample_ctrl=res$sample_ctrl
RIN.simu=res$RIN.simu

random_idx_gene = sample.int(nTotal)
random_idx_sam  = sample.int(nall)

sample_ctrl = sample_ctrl[random_idx_gene, random_idx_sam, ]

```



Once we got the basic parameters for simulation, we add size factors by introducing the differences between cases and controls.
```{r adjust_param,eval=runcode}

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

#gene x ind x (mean,overdispersion, droupout) x #of mixture models(2)
sample_param_case = abind::abind(sample_ctrl[, 1:ncase, ], 
                                 sample_ctrl[, 1:ncase, ], along = 4)

sample_param_ctrl = abind::abind(sample_ctrl[, (ncase + 1):nall, ], 
                                 sample_ctrl[, (ncase + 1):nall, ], along = 4)




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
temp = sample_param_case[i_case_modify, , , ]
sample_param_case[i_case_modify, , , ] = sample_param_ctrl[i_case_modify, , , ]
sample_param_ctrl[i_case_modify, , , ] = temp

```

## Simulate Counts from Parameters of ZINB Model
Once got the parameters, we simulate the counts of each cell and each genes from zero inflated negative binomial models.
```{r param2count,eval=runcode}
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
    sub_sample_param_case[which(de.dp == 1), , ] = sample_param_case[which(de.dp ==
                                                                             1), , , cell_distr_flag_dp_add] #equally separated for dp as design
    sub_sample_param_ctrl[which(de.dp == 1), , ] = sample_param_ctrl[which(de.dp ==
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
      sim_param[ig, idx_i[k], ] = c(sample_mean_k[ig], disp_i[ig], drop_i[ig])
    }
  }
  #print(i)
}
date()

dim(sim_matrix)
sim_matrix[1:4, 1:4]

table(c(sim_matrix) == 0)
table(c(sim_matrix) == 0) / (nrow(sim_matrix) * ncol(sim_matrix))

```

We also simulate the meta information and collect some basic statistics of the counts we simulated.

```{r meta_simulation,eval=runcode}
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
meta[1:2, ]

# Some basic stat & Preparation for Bulk RNAseq analysis 

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

#Other cell level info
cell_readdepth = apply(sim_matrix, 2, sum)
cell_count_quantile = apply(sim_matrix, 2, quantile)

```

# Simulated Differential Expressed Analysis
We calculate the p-values of the simulated genes through the BSDE method as well as the [MAST](https://www.bioconductor.org/packages/release/bioc/html/MAST.html) (with a widely used protocol mixed-effect model "glmer" & Likelihood-ratio Test) and the [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

## DESeq2

For DESeq2 analysis, we construct the bulk RNAseq situation though suming up the counts per gene, per cluster, per individual, per 1000 cells.

```{r deseq2,eval=runcode}

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
```

## MAST

We take the method of current usage of MAST, we get the p-values based on the Hurdle model ("H" model)
```{r MAST,eval=runcode}

sim_matrix_log = log2(1 + sim_matrix) #log transformed data

dim(sim_matrix_log)
sim_matrix_log[1:2, 1:5]

cell_id = colnames(sim_matrix_log)   #get the cell id from the data
gene_id = rownames(sim_matrix_log)   #get the gene id from the data

fData = data.frame(primerid = gene_id)
cData = data.frame(wellKey = cell_id)

diagnosis = as.character(meta$phenotype) #

diagnosis2=matrix("Control",ncol=1,nrow=length(diagnosis))
diagnosis2[which(diagnosis == "1")] = "Case"
diagnosis= as.factor(diagnosis2)

sca = MAST::FromMatrix(sim_matrix_log, cData, fData)
colData(sca)$cngeneson = as.numeric(meta$CDR)
colData(sca)$diagnosis = as.factor(diagnosis)
colData(sca)$ind = as.factor(meta$individual)

colData(sca)

print(date())
print(gc())
b1 = tryCatch(MAST::zlm(formula = ~ diagnosis + ( 1 | ind ), sca = sca, method = 'glmer', ebayes = FALSE, parallel = TRUE), error = function(e) {NA} )
print(date())
print(gc())

lrt1 = tryCatch(MAST::lrTest(b1, "diagnosis"), error = function(e) {NA} )
MAST_pval1 = tryCatch(apply(lrt1, 1, function(x){x[3,3]}), error = function(e) {NA} )
print(date())
print(gc())


```

## BSDE
We also calculate the p-values using the BSDE as follows.

```{r BSDE,eval=runcode}
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

use_python("/nas/longleaf/apps/python/3.5.1/bin/python3")
source("~/Dropbox/Barycenter/BSDEpackage/BSDE/R/op_functions.R")
py_run_file("~/Dropbox/Barycenter/BSDEpackage/op_functions.py")
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

```

# Further Analysis and Method Evaluation
Once we have the p-values, we can do further analysis on the method evaluation.For example, we can calculate the proportion of p-values < 0.05. We can also draw similar plot for power, where power is calculated as the proportion of p-values < 0.05 in permuted data. 

