#' DINA
#'
#' This function adapts the predict.lm for DINA.
#' @method predict DINA
#' @importFrom stats predict
#'
#' @param object
#' @param newdata the data frame in which to look for variables with which to predict. If omitted, the fitted values are used. Default is NULL.
#' @param interval: "none", "confidence", "prediction?". Type of interval calculation. If omitted, point estimates are provided. Default is NULL.
#' @param level confidence level. Default is 0.95.
#' @param sd.method "bootstrap", "approximation".
#'
#' @export
predict.DINA = function(object, ..., interval = "none", sd.method = "bootstrap"){

}
