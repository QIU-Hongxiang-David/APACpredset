#' @importFrom assertthat assert_that is.string has_name noNA is.count is.flag is.number
#' @importFrom DescTools BinomCI
#' @import SuperLearner


#' @title APACpredset package
#' 
#' @description
#' In this package, we implement methods to calculate Asymptotically Probably Approximately Correct (APAC) prediction sets under unknown covariate shift. Inputs are
#' \itemize{
#' \item{Observations from source population}{This includes both covariates \eqn{X}{X} and outcome/dependent variable/label \eqn{Y}{Y}}
#' \item{Observations from target population}{This includes covariates \eqn{X}{X} while \eqn{Y}{Y} may be missing}
#' \item{A pretrained scoring function}{This function \eqn{s}{s} takes in values of \eqn{(X,Y)}{(X,Y)} and outputs a score of how likely it is to observe \eqn{Y}{Y} given observing \eqn{X}{X}. This function should be pretrained from an independent data set that is typically drawn from the source population, and may be trained with the aid of data from the target population with missing \eqn{Y}{Y}}
#' \item{Candidate thresholds}{A finite set of candidate thresholds. (More info below)}
#' }
#' We consider prediction sets of the form
#' \deqn{C_\tau(x)=\{y: s(x,y) \ge \tau\}}{C_tau(x)={y: s(x,y) >= tau}}
#' where \eqn{\tau}{tau} is a threshold.
#' 
#' A prediction set \eqn{\hat{C}_n}{hat(C)_n} estimated from data is said to be APAC if
#' \deqn{P(P(Y \notin \hat{C}_n(X)) \le \epsilon | \hat{C}_n) \ge 1-\delta-o(1)}{P(P(Y notin hat(C)_n(X)) <= epsilon | hat(C)_n)>=1-delta-o(1)}
#' as sample size \eqn{n}{n} tends to infinity. Here, \eqn{(X,Y)}{(X,Y)} is an new draw from the target population that is independent of \eqn{\hat{C}_n}{hat(C)_n}, \eqn{\epsilon}{epsilon} is a desired upper bound on the true coverage error, and \eqn{1-\delta}{1-delta} is a desired confidence level of the statement on the coverage error of \eqn{\hat{C}_n}{hat(C)_n}. The methods in this package aim at selecting a threshold from candidate thresholds that leads to an APAC prediction set that is as small as possible (namely the threshold is as large as possible).
#' 
#' When either \eqn{Y}{Y} is observed in data from the target population or the likelihood ratio between the target population and source population is known, other methods that produce PAC prediction sets may be preferrable to those in this package because PAC is a finite-sample guarantee while APAC is large-sample property.
#' 
#' @docType package
#' @name APACpredset-package
NULL
