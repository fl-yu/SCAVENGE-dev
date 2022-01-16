#' @title get_sigcell_simple
#'
#' @description a function to determine statistically trait-enriched cell by permutation test
#'
#' @param knn_sparse_mat a sparse matrix used for network propagation, which indicates the adjacent matrix (m x m, where m is the cell number) of cell-to-cell network (M-kNN graph)
#' @param topidx a logical vector indicating seed cells (TRUE) and non-seed cells (FALSE) with length of m, where m is the cell number. The length and position are corresponding to knn_sparse_mat
#' @param permutation_times an integer describe times of permutation test for each cell
#' @param true_cell_significance a numeric value between 0-1 indicating the significant threshold used to determine statistically trait-enriched cell
#' @param out_rda path and name for resulting r.data
#' @param mycores how many cores used for permutation test
#'
#' @return a dataframe with two columns, the first one is seed index (the same with input), the second one is significance from permutation test. In addition, an R data is generated, which contains this dataframe and another dataframe depicting the network propagation score for each cell at each permuation
#' @export
#'
#' @import Matrix
#' @import parallel
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#' get_sigcell_simple(knn_sparse_mat=mutualknn30,
#' topidx=seed_p0.05,
#' permutation_times=1000,
#' true_cell_significance=0.05,
#' out_rda="true_cell_df.rda",
#'  mycores=4)}
#'
get_sigcell_simple <- function(knn_sparse_mat=mutualknn30,
                               topidx=seed_p0.05,
                               permutation_times=1000,
                               true_cell_significance=0.05,
                               out_rda="true_cell_df.rda",
                               mycores=4){
  if(permutation_times<100){
    warning("Permutation times less than 100")
  }
  cell_mat <- data.frame(cell=1:nrow(knn_sparse_mat), degree=colSums(knn_sparse_mat))
  cell_table <- data.frame(table(cell_mat$degree))

  seed_mat_top  <- data.frame(seed=which(topidx), degree=colSums(knn_sparse_mat[, topidx]))
  # summary(seed_mat_top $degree)
  seed_table_top <- data.frame(table(seed_mat_top$degree))
  xx_top <- tapply(cell_mat[, 1], cell_mat[, 2], list)
  xx2_top <- xx_top[names(xx_top) %in% seed_table_top$Var1]
  # permutation_score_top <- data.frame(matrix(nrow=nrow(knn_sparse_mat), ncol=1000))
  permutation_score_top <- mclapply(1:permutation_times, mc.cores = mycores, function(i){
    sampled_cellid <- xx2_top %>%
      mapply(sample, ., seed_table_top$Freq) %>%
      unlist %>%
      sort
    xx <- randomWalk_sparse(intM=knn_sparse_mat, queryGenes=as.numeric(sampled_cellid), gamma=0.05)
    if (i %% 100 == 0) {print(i)}
    return(xx)
  }
  )
  names(permutation_score_top) <- paste0("permutation_", 1:permutation_times)
  permutation_score_top <- as.data.frame(permutation_score_top)

  permutation_df_top <- data.frame(matrix(nrow=nrow(knn_sparse_mat), ncol=permutation_times))

  permutation_df_top <- apply(permutation_score_top, 2, function(x) { temp <- x<=topseed_trs; return(temp) } )
  message("cells passed 0.001 threshold: ", (sum(rowSums(permutation_df_top)>=0.001*permutation_times)*100)/nrow(permutation_df_top))
  message("cells passed 0.01 threshold: ", (sum(rowSums(permutation_df_top)>=0.01*permutation_times)*100)/nrow(permutation_df_top))
  message("cells passed 0.05 threshold: ", (sum(rowSums(permutation_df_top)>=0.05*permutation_times)*100)/nrow(permutation_df_top))
  true_cell_top_idx <- rowSums(permutation_df_top)>=(true_cell_significance*permutation_times)
  message("Your threshold: ", true_cell_significance)
  message("Fold of true cell over seed: ", sum(true_cell_top_idx)/sum(topidx)) # fold of true cell over seed
  message("How many propertion of seed were true cells: ", sum(true_cell_top_idx & topidx)*100/sum(topidx)) # how many propertion of seed were true cells
  message("How many propertion true cell over all cells: ", (sum(true_cell_top_idx)*100)/nrow(permutation_df_top)) # how many propertion true cell over all cells
  message("How many propertion of seed over all cells: ", (sum(topidx)*100)/nrow(permutation_df_top)) # how many propertion of seed over all cells

  true_cell_top_filter_idx <-
    ture_cell_df <- data.frame(topidx,
                               true_cell_top_idx)
  save(ture_cell_df, permutation_score_top, file=out_rda)
  return(ture_cell_df)
}
