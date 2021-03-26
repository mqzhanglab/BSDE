decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#' This function calculates the w2_pval based on BSDE for a given gene
#' @param count_per_gene is a len x 1 matrix,length=cell num.  It represents gene count.
#' @param meta_individual is a len x 1 vector, length=cell num. It represents individual labels.
#' @param meta_individual is a len x 1 vector, length=cell num. It is indicator labels (contains 0,1) for cases (1) and controls(0)
#' @param perm_num is a ingeter for times of permutation.
#' @param merge_method provides the support package for barycenter calculation,
#' the default setting is based on the python package 'ot'(merge_method="python"),
#' if python doesn't work, the R package merge_method="R" is also available
#' @param unif_round_unit is the unit length for frequency calculation.
#' @param l_fold is a positive number, represents the fold-change to the parameter lamda for merge_method R, the lambda is the Non-negative regularization parameter (for small lambda the Barycenter is close to the EMD).
#' @param weight is a len x 1 vector, length=subjects num. It represents the importance for subjects.
#' @param shrink is a bool value. If True, the range with zero density will be removed, for more robust performance to extreme values and reduce the computational burden. Choose TRUE if there are some extereme high expressions.
#' @return a list contains:
#' pval a p-value, from the based on Monte Carlo Permutation Procedure.
#' case_bc_ob a vectors contains the density of case subjects distribution of Barycenter
#' ctrl_bc_ob a vectors contains the density of ctrl subjects distribution of Barycenter
#'
#' @export
#' @examples
#' library("Barycenter")
#' library("Rcpp")
#' library("reticulate")
#' library("doRNG")
#' library("doParallel")
#' py_run_file("../inst/op_functions.py")
#'
#' count_per_gene=c(rpois(60,6),c(rpois(60,4)))
#' meta_individual=paste0("ind",rep(1:12,each=10))
#' meta_phenotype=c(rep(1,60),rep(0,60))
#' cal_w2_pval(count_per_gene,meta_individual,meta_phenotype )
#' cal_w2_pval(count_per_gene,meta_individual,meta_phenotype,unif_round_unit = 1 )

cal_w2_pval=function(count_per_gene,meta_individual,meta_phenotype,perm_num=200,unif_round_unit=0.2,l_fold=1,merge_method="python",weight=1,shrink=FALSE){ #num_method="unif", "empirical"
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
    cur_freq=hist(cur_ind_count,breaks=cur_range,plot=FALSE)$counts
    cur_distr[,i_n]=cur_freq/sum(cur_freq)
  }

  if(shrink){ #however, typically, we don't shrink, since it may just by chance that there are no numbers between
    shrink_index=which(apply(cur_distr,1,sum)>0)  #remove zero counts, to fastern the calculation.
    cur_burden=cur_burden[shrink_index,shrink_index]
    cur_distr=cur_distr[shrink_index,]
  }

  w2_res=foreach(ip=0:perm_num) %dorng% {
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
      w_case=weight[cur_phenotype==1]
      w_ctrl=weight[cur_phenotype==0]
    }

    #1. R method
    if(merge_method=="R"){
      #case_list=lapply(seq_len(ncol(case_distr)), function(i) cbind(as.matrix(case_distr[,i]),background_image))
      #ctrl_list=lapply(seq_len(ncol(ctrl_distr)), function(i) cbind(as.matrix(ctrl_distr[,i]),background_image))
      case_list=lapply(seq_len(ncol(case_distr)), function(i) as.matrix(case_distr[,i]))
      ctrl_list=lapply(seq_len(ncol(ctrl_distr)), function(i) as.matrix(ctrl_distr[,i]))
      cur_costm=cur_burden
      cur_costm=cur_costm/max(cur_costm)
      case_bc=WaBarycenter2(case_list,lambda_fold=l_fold,costm = cur_costm)
      ctrl_bc=WaBarycenter2(ctrl_list,lambda_fold=l_fold,costm = cur_costm)
    }

    #2. python POT default method
    if(merge_method=="python"){
      #method 1 python
      py$case_wass=NULL
      py$ctrl_wass=NULL
      py$case_distr=r_to_py(case_distr)
      py$ctrl_distr=r_to_py(ctrl_distr)
      py$w_case=r_to_py(w_case)
      py$w_ctrl=r_to_py(w_ctrl)

      py_run_string("case_wass = cal_bary_wass(case_distr,w=w_case)")
      py_run_string("ctrl_wass = cal_bary_wass(ctrl_distr,w=w_ctrl)")
      case_bc=py$case_wass
      ctrl_bc=py$ctrl_wass
    }


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




#' This function calculates the w2_pval based on BSDE for a given gene, based on bulk RNAseq.
#' The expression are represented as one numerical value per gene, per individual.
#'
#' @param count_per_gene: len x 1 matrix,length=ind num,  gene count
#' @param meta_individual: len x 1 vector, length=ind num, individual labels
#' @param meta_individual: len x 1 vector, length=ind num, indicator labels (contains 0,1)
#' @param perm_num:ingeter, permutation time
#' @return a list contains:
#' pval a p-value, from the based on Monte Carlo Permutation Procedure.
#' case_distr a vectors contains the density of case subjects distribution of Barycenter
#' ctrl_distr a vectors contains the density of ctrl subjects distribution of Barycenter
#'
#' @export
#' @examples
# count_per_gene=c(rpois(6,60),c(rpois(6,4)))
# cur_individual=paste0("ind",rep(1:12,each=1))
# cur_phenotype=c(rep(1,6),rep(0,6))
# cal_w2_bulk_pval(count_per_gene,meta_individual,meta_phenotype )

cal_w2_bulk_pval=function(cur_ind_count,cur_individual,phenotype,perm_num=500,unif_round_unit=0.5){ #num_method="unif", "empirical"
  n=length(cur_individual)
  #cur_ind_count=round(cur_ind_count,unif_round_unit)
  cur_range=seq(min(cur_ind_count),(max(cur_ind_count)+unif_round_unit),by=unif_round_unit) #cur range +1
  unif_round_unit_digit=min(decimalplaces(unif_round_unit),5)
  cur_range=floor(cur_range*10^unif_round_unit_digit)/10^unif_round_unit_digit
  cur_range=c(cur_range, (max(cur_range)+cur_range[2]-cur_range[1])) #cur range expand +2
  len=length(cur_range)-1
  cur_burden=matrix((rep(cur_range[1:len],times=len)-rep(cur_range[1:len],each=len))^2,len,len)


  w2_res=foreach(ip=0:perm_num) %dorng% {
    if(ip>0){
      cur_phenotype=phenotype[sample.int(n,n)]
    }
    if(ip==0){
      cur_phenotype=phenotype
    }

    case_distr=hist(cur_ind_count[cur_phenotype==1],breaks=cur_range,plot=FALSE)$counts
    ctrl_distr=hist(cur_ind_count[cur_phenotype==0],breaks=cur_range,plot=FALSE)$counts


    names(case_distr)=cur_range[1:len]
    names(ctrl_distr)=cur_range[1:len]

    #calculate w2
    w2=Greenkhorn(as.matrix(case_distr),(as.matrix(ctrl_distr)),costm=abs(cur_burden), lambda=0.01)$Distance
    #w2
    list(w2,case_distr,ctrl_distr)
  }
  #w2_ob=w2_res[[1]]
  #w2_perm=sapply(w2_res,unlist)[-1]
  w2_ob=w2_res[[1]][[1]]
  w2_perm=sapply(w2_res,unlist)[-1]
  case_distr_ob=w2_res[[1]][[2]]
  ctrl_distr_ob=w2_res[[1]][[3]]

  pval=mean(w2_perm>=w2_ob,na.rm=TRUE)

  return(list(pval,case_distr,ctrl_distr))
}



#' This function is modified for 1-dimensional data from WaBarycenter function of R package 'Barycenter', and used for the option "R" of merge_method at function cal_w2_pval
#' A list of matrices satisfying the prerequisites described above.
#' @param images The input vector, which could be the log-transformed and normalized expression counts.
#' @param maxIter Maximum number of iterations.
#' @param lambda Non-negative regularization parameter (for large lambda the regularized Barycenter is close to its true counterpart). If FALSE the algorithm uses a lambda depending on costm.
#' @param costm A matrix of pairwise distances between the locations. If FALSE the algorithm uses the usual euclidean distance matrix on a [0,1]x[0,1] equidistant pixel grid.
#' @param lambda_fold A modification for the default lambda adjustment.
#' @return Value The Barycenter of the input, represented by a  vector.

#' @references
#' Cuturi, M.: Fast Computation of Wasserstein Barycenters, Proceedings of the International Conference on Machine Learning, Beijing, China, 2014


WaBarycenter2=function (images, maxIter = 10, lambda = FALSE, costm = FALSE,lambda_fold=1) {
  #time <- proc.time()
  if (is.list(images) == FALSE) {
    stop("The images have to be passed as a list each entry representing a matrix!")
  }
  if (length(unique(lapply(images, dim))) == 1) {
    dimension <- dim(images[[1]])
  }
  if (length(unique(lapply(images, dim))) != 1) {
    stop("Dimensions of the images are not equal!")
  }
  if (is.matrix(costm) == FALSE) {
    n <- dimension[1] * dimension[2]
    coord1 <- seq(0, 1, length.out = dimension[2])
    coord2 <- rev(seq(0, 1, length.out = dimension[1]))
    coordinates <- expand.grid(coord1, coord2)
    costm <- as.matrix(dist(coordinates, diag = TRUE, upper = TRUE))
  }
  if (is.matrix(costm) != FALSE) {
    n <- dimension[1] * dimension[2]
    if (identical(dim(costm), rep(n, 2)) == FALSE) {
      print(costm)
      stop("Dimension of the cost matrix is not compatible with the given images!")
    }
  }
  if (lambda == FALSE) {
    lambda <- 60/median(costm)*lambda_fold
  }
  a_tild <- rep(1/n, n)
  a_hat <- rep(1/n, n)
  t_0 <- 2
  t <- t_0
  for (i in 1:maxIter) {
    beta <- (t + 1)/2
    a <- (1 - 1/beta) * a_hat + (1/beta) * a_tild
    ALPHA <- 0
    for (j in 1:length(images)) {
      b=Barycenter:::Subgradient(a, t(images[[j]]), costm, lambda)
      b[is.na(b)]=0
      ALPHA <-  b+ ALPHA
    }
    ALPHA <- (1/length(images)) * ALPHA
    a_tild <- a_tild * exp(-(t_0) * beta * ALPHA)
    a_tild <- a_tild/sum(a_tild)
    a_hat <- (1 - 1/beta) * a_hat + (1/beta) * a_tild
    t <- t + 1
  }
  # if (length(unique(dimension)) == 1) {
  a.temp <- matrix(a, dimension[1], dimension[2], byrow = TRUE)
  # a.temp <- a.temp[, nrow(a.temp):1]
  a.temp <- a.temp[, ncol(a.temp):1]
  # }
  #print(proc.time() - time)
  return(a.temp)
}





