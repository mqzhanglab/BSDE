decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#' Calculate the BSDE p-value for a SingleCellExperiment object.
#'
#' @param sce A SingleCellExperiment object.
#'   See \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param perm_num The number of permutations for computing Monte Carlo p-value.
#' @param unif_round_unit Bin width for empirical histogram (default: 0.2).
#' @param weight Weight vector (non-negative, length=number of subjects) for subjects.
#'   Can be used for inverse probability weighting to control for confounding.
#'   Default: 1 for equal weights.
#' @param shrink If \code{TRUE}, the range with zero count will be removed to
#'   cope with sparsity and speed up the computation. Use \code{TRUE} if
#'   extremely high expressions are present. Default: \code{FALSE}.
#' @return a list of p-values
#' @export
BSDE <- function(sce, perm_num=200, unif_round_unit=0.2, weight=1, shrink=FALSE) {
  sce <- scater::logNormCounts(sce)
  count_log <- SingleCellExperiment::logcounts(sce)
  cur_ind <- SingleCellExperiment::colData(sce)$individual
  cur_condition <- SingleCellExperiment::colData(sce)$condition
  gene_id <- SingleCellExperiment::rowData(sce)$gene_id
  op_pval <- plyr::laply(1:nrow(count_log),
                        function(i_g) {
                          cur_count <- count_log[i_g, ]
                          if(sum(cur_count)==0){
                            return(1)
                          } else {
                          return(tryCatch({
                            BSDE::cal_w2_pval(
                              count_per_gene = cur_count,
                              meta_individual = cur_ind,
                              meta_phenotype = cur_condition,
                              perm_num = perm_num,
                              unif_round_unit = unif_round_unit,
                              weight = weight,
                              shrink = shrink
                            )[[1]]
                          }, error = function(e) {NA}))
                        }}, .progress = "text")
  rownames(op_pval) = gene_id
  op_pval
}
#' This function calculates the BSDE Barycenter p-value of a given gene.
#'
#' @param count_per_gene Gene count. Vector of length n (n is number of cells).
#' @param meta_individual Label of individuals (length n).
#' @param meta_phenotype 0/1 vector of length n (0=control, 1=case).
#' @param perm_num Number of random permutations for computing p-value.
#' @param unif_round_unit Bin width for empirical histogram (default: 0.2).
#' @param weight Weight vector (non-negative, length=number of subjects) for subjects.
#'   Can be used for inverse probability weighting to control for confounding.
#'   Default: 1 for equal weights.
#' @param shrink If \code{TRUE}, the range with zero count will be removed to
#'   cope with sparsity and speed up the computation. Use \code{TRUE} if
#'   extremely high expressions are present. Default: \code{FALSE}.
#' @return a list contains:
#' \itemize{
#'   \item pval: p-value
#'   \item case_bc_ob: density of the aggregated Barycenter distribution of the cases
#'   \item ctrl_bc_ob: density of the aggregated Barycenter distribution of the controls
#' }
#' @export
cal_w2_pval <- function(count_per_gene, meta_individual, meta_phenotype,
                        perm_num=200, unif_round_unit=0.2,
                        weight=1, shrink=FALSE){
  cur_individual=unique(meta_individual)
  phenotype=meta_phenotype[match(cur_individual,meta_individual)]

  n=length(cur_individual)

  #count_per_gene=round(count_per_gene,unif_round_unit)
  cur_range=seq(min(count_per_gene),(max(count_per_gene)+unif_round_unit),by=unif_round_unit) #cur range +1
  unif_round_unit_digit=min(decimalplaces(unif_round_unit),5)
  cur_range=floor(cur_range*10^unif_round_unit_digit)/10^unif_round_unit_digit
  cur_range=c(cur_range, (max(cur_range)+cur_range[2]-cur_range[1])) #cur range expand +2
  len=length(cur_range)-1
  cur_burden=matrix((rep(cur_range[1:len],times=len)-rep(cur_range[1:len],each=len))^2,len,len)
  #cur_burden=cur_burden/(10^unif_round_unit)
  cur_distr=matrix(0,nrow=len,ncol=n) #row for ranges col for individuals
  #background_image=marix(0,nrow=len,ncol=(len-1))
  colnames(cur_distr)=cur_individual
  rownames(cur_distr)=cur_range[1:len]

  for(i_n in 1:n){
    cur_ind=cur_individual[i_n]
    cur_ind_count=count_per_gene[meta_individual==cur_ind]
    cur_freq=graphics::hist(cur_ind_count,breaks=cur_range,plot=FALSE)$counts
    cur_distr[,i_n]=cur_freq/sum(cur_freq)
  }

  if(shrink){ #however, typically, we don't shrink, since it may just by chance that there are no numbers between
    #remove zero counts, to fasten the calculation.
    shrink_index=which(apply(cur_distr,1,sum)>0)
    cur_burden=cur_burden[shrink_index,shrink_index]
    cur_distr=cur_distr[shrink_index,]
    cur_range=cur_range[shrink_index]
    len=length(cur_range)-1
  }

  `%dorng%` <- doRNG::`%dorng%`
  `%dopar%` <- foreach::`%dopar%`
  w2_res=foreach::foreach(ip=0:perm_num) %dorng% {
  #for(ip in 0:10) {
    #print(ip)
    if(ip>0){
      cur_phenotype=phenotype[sample.int(n,n)]
    }
    if(ip==0){
      cur_phenotype=phenotype
    }
    case_distr=cur_distr[,cur_phenotype==1]
    ctrl_distr=cur_distr[,cur_phenotype==0]

    w_case=1
    w_ctrl=1
    if(length(weight)!=1){
      stopifnot("weight length must match # of subjects"=length(weight)==length(cur_phenotype))
      w_case=weight[cur_phenotype==1]
      w_ctrl=weight[cur_phenotype==0]
      w_case =  w_case / sum(w_case)
      w_ctrl =  w_ctrl / sum(w_ctrl)
    }

    py <- reticulate::import_main()
    py$case_wass=NULL
    py$ctrl_wass=NULL
    py$case_distr=reticulate::r_to_py(case_distr)
    py$ctrl_distr=reticulate::r_to_py(ctrl_distr)
    py$w_case=reticulate::r_to_py(w_case)
    py$w_ctrl=reticulate::r_to_py(w_ctrl)

    reticulate::py_run_string("case_wass = cal_bary_wass(case_distr,w=w_case)")
    reticulate::py_run_string("ctrl_wass = cal_bary_wass(ctrl_distr,w=w_ctrl)")
    case_bc=py$case_wass
    ctrl_bc=py$ctrl_wass
    names(case_bc)=cur_range[1:len]
    names(ctrl_bc)=cur_range[1:len]
    #calculate w2

    w2=Barycenter::Greenkhorn(as.matrix(case_bc),(as.matrix(ctrl_bc)),costm=abs(cur_burden), lambda=0.01)$Distance
    #w2
    list(w2,case_bc,ctrl_bc)
  }
  #w2_ob=w2_res[[1]]
  #w2_perm=sapply(w2_res,unlist)[-1]
  w2_ob=w2_res[[1]][[1]]
  w2_perm=sapply(w2_res[-1],function(x)unlist(x[[1]]))
  case_bc_ob=w2_res[[1]][[2]]
  ctrl_bc_ob=w2_res[[1]][[3]]

  pval=mean(w2_perm>=w2_ob,na.rm=TRUE)
  return(list(pval,case_bc_ob,ctrl_bc_ob))
}

#' Calculate p-value for a given gene based on bulk RNAseq.
#' The expressions are represented as one numerical value per gene, per individual.
#'
#' @param cur_ind_count vector of gene counts (length = number of individuals).
#' @param cur_individual individual annotations (length = number of individuals).
#' @param phenotype 0/1 individual labels (length = number of individuals, 0=control, 1=case).
#' @param perm_num number of permutations for computing p-value (integer).
#' @param unif_round_unit Bin width for empirical histogram (default: 0.2).
#' @return a list contains:
#' \itemize{
#' \item pval: p-value
#' \item case_distr: density of the aggregated Barycenter distribution of the cases
#' \item ctrl_distr: density of the aggregated Barycenter distribution of the controls
#' }
#' @export
#' @examples
#' \dontrun{
#' count_per_gene=c(rpois(6,10),c(rpois(6,4)))
#' meta_individual=paste0("ind",rep(1:12,each=1))
#' meta_phenotype=c(rep(1,6),rep(0,6))
#' cal_w2_bulk_pval(count_per_gene,meta_individual,meta_phenotype)
#' }
#'
cal_w2_bulk_pval=function(cur_ind_count,cur_individual,phenotype,
                          perm_num=500,unif_round_unit=0.2){
  n=length(cur_individual)
  cur_range=seq(min(cur_ind_count),(max(cur_ind_count)+unif_round_unit),by=unif_round_unit) #cur range +1
  unif_round_unit_digit=min(decimalplaces(unif_round_unit),5)
  cur_range=floor(cur_range*10^unif_round_unit_digit)/10^unif_round_unit_digit
  cur_range=c(cur_range, (max(cur_range)+cur_range[2]-cur_range[1])) #cur range expand +2
  len=length(cur_range)-1
  cur_burden=matrix((rep(cur_range[1:len],times=len)-rep(cur_range[1:len],each=len))^2,len,len)

  `%dorng%` <- doRNG::`%dorng%`
  `%dopar%` <- foreach::`%dopar%`
  w2_res=foreach::foreach(ip=0:perm_num) %dorng% {
    if(ip>0){
      cur_phenotype=phenotype[sample.int(n,n)]
    }
    if(ip==0){
      cur_phenotype=phenotype
    }
    case_distr=graphics::hist(cur_ind_count[cur_phenotype==1],breaks=cur_range,plot=FALSE)$counts
    ctrl_distr=graphics::hist(cur_ind_count[cur_phenotype==0],breaks=cur_range,plot=FALSE)$counts
    names(case_distr)=cur_range[1:len]
    names(ctrl_distr)=cur_range[1:len]
    w2=Barycenter::Greenkhorn(as.matrix(case_distr),(as.matrix(ctrl_distr)),costm=abs(cur_burden), lambda=0.1)$Distance
    list(w2,case_distr,ctrl_distr)
  }
  w2_ob=w2_res[[1]][[1]]
  w2_perm=sapply(w2_res[-1],function(x)unlist(x[[1]]))
  case_distr=w2_res[[1]][[2]]
  ctrl_distr=w2_res[[1]][[3]]

  pval=mean(w2_perm>=w2_ob,na.rm=TRUE)
  return(list(pval,case_distr,ctrl_distr))
}


