# helper functions for DINA

#' sample.split
#'
#' This function splits weighted samples into folds with (approximately) equal weight sum.
#'
#' @param weight weight vector
#' @param k number of folds. Default is 2.
#'
#' @return a matrix of samples. Each column stands for a fold. Each sample is repeated for \code{weight} times.
#' @export
sample.split = function(weight, k){
  weight.sum = sum(weight)
  n = length(weight)
  weight = c(weight, 0)
  weight.order = order(weight, decreasing = T)
  weight = weight[weight.order]

  membership = matrix(c(weight.order[weight.sort > 1], rep(n+1, 2 * k - sum(weight.sort > 1))), nrows = 2 * k)
  membership[seq(k+1, 2 *k), ] = membership[seq(2 *k, k+1), ]
  weight.sum = apply(membership, 1, function(x){return(sum(weight[x]))})
  weight.sum.rest = c(rep(ceiling(n/k), n %% k), rep(floor(n/k), k - n %% k)) - weight.sum
  weight.order[cumsum]

  result = list()
  for(i in 1:k){
    result[[i]] = rep(weight.order, weight.order + 1)
  }
  result[[]]
  membership = matrix(c(weight.order[weight.sort > 1], rep(n+1, 2 * k - sum(weight.sort > 1))), nrows = 2 * k)
  result = list()
  # weight >= 2
  weight.order = order(weight, decreasing = T)
  weight.order = sort(weight[weight.order], decreasing = T);
  # weight == 1

  # weight == 0

}

#' asymp.variance
#'
#' This function estimates the variance matrix of DINA by the asymptotic variance matrix with reasonably good nuisance function estimators and regularity conditions
#'
#' @param a structure?
#'
#' @return the variance matrix
#' @export
asymp.variance = function(){
  # TODO
  result$Sigma1 = t(cbind(1, X)) %*% diag((data.HTE$resW)^2) %*% cbind(1, X)/length(Y) # add inference ...
  result$Sigma2 = t(cbind(1, X)) %*% diag((data.HTE$resW * glm.model$residuals)^2) %*% cbind(1, X)/length(Y)

  result = list(); result$estimator = record[[1]]
  result$Sigma1 = apply(array(unlist(Sigma1), dim = c(d+1,d+1,k)), c(1,2), mean)
  result$Sigma2 = apply(array(unlist(Sigma2), dim = c(d+1,d+1,k)), c(1,2), mean)
  result$Sigma = solve(result$Sigma1) %*% result$Sigma2 %*% solve(result$Sigma1)
}




