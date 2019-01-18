#' meanTruncPareto
#'
#' Calculate the mean of a truncated Pareto distribution?
#' TODO: description and title needed
#'
#' @param k TODO: description needed
#' @param lower TODO: description needed
#' @param upper TODO: description needed
#' @param alpha TODO: description needed
#'
#' @return TODO: description needed
#'
#' @export
meanTruncPareto <- function(k, lower, upper, alpha) {
  k * lower^k * (upper^(1 - k) - alpha^(1 - k)) / ((1 - k) * (1 - (alpha/upper)^k))
}
