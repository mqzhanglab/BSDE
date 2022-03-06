#' Meta information
"meta"

#' \code{SingleCellExperiment::SingleCellExperiment} object
"sce"

#' count data (row: gene, column: cell)
"sim_matrix"

#' indicator vector (1 for differential expression in proportion).
"de.dp"

#' indicator vector (1 for differential expression in mean)
"de.mean"


#' indicator vector (1 for differential expression in multimodality)
"de.mult"

#' indicator vector (1 for differential expression in variance)
"de.var"

#' matrix for bulk RNA data (gene x individual)
"sim_matrix_bulk"
