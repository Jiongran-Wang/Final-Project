#' @title Calculate the representative cluster label from posterior samples
#'
#' @param C A matrix (N * S) containing posterior samples of labels
#' @param burnin The number of samples will be discarded
#' @param K The predefined number of clusters
#'
#' @return A vector containing the representative label
#' @export
#'
#' @examples representative_cluster(C = matrix(sample(10, 5000, replace = TRUE), nrow = 10), burnin = 250, K = 10)
representative_cluster <- function(C, burnin, K){
  N <- nrow(C)
  S <- ncol(C)
  p <- matrix(NA, nrow = N, ncol = K)
  for (i in 1:K){
    p[,i] <- apply(C[, (burnin+1):S],
                   MARGIN = 1,
                   FUN = function(x){return(sum(x == i))}) / ncol(C[, (burnin+1):S])
  }
  predicted_label <- apply(p, MARGIN = 1, FUN = function(x){return(which.max(x))})
  return(predicted_label)
}
