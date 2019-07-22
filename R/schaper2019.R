
#' Ratings for 96 objects
#'
#' An item set that was used in experiments by Schaper, Kuhlmann and
#' Bayen (2019; in press).
#'
#' @format A data frame with 96 rows and 7 variables
#' \describe{
#'   \item{item}{The name of the object (in German)}
#'   \item{room}{The room in which the item is typically found; can be
#'       'kitchen' or 'bathroom'}
#'   \item{rating_consistent}{How expected would it
#'       be to find the \code{item} in the typical \code{room}}
#'   \item{rating_inconsistent}{How expected would it
#'       be to find the \code{item} in the atypical \code{room}}
#'   \item{syllables}{The number of syllables in the object name}
#'   \item{frequency}{A value indicating the relative frequency of the
#'       object name in German language}
#'   \item{list}{Represents the set affiliation of the \code{item} as
#'       realized in experiments by Schaper et al.}
#' }
#'
#' @source
#' Courteously provided by Marie Lusia Schaper and Ute Bayen.
#'
#' @examples
#'
#' head(schaper2019)
#' features <- schaper2019[, 3:6]
#'
#' # Optimize the variance criterion
#' # (tends to maximize similarity in feature means)
#' anticlusters <- anticlustering(
#'   features,
#'   K = 3,
#'   objective = "variance",
#'   categories = schaper2019$room,
#'   method = "exchange"
#' )
#'
#' # Means are quite similar across sets:
#' by(features, anticlusters, function(x) round(colMeans(x), 2))
#' # Check differences in standard deviations:
#' by(features, anticlusters, function(x) round(apply(x, 2, sd), 2))
#' # Room is balanced between the three sets:
#' table(Room = schaper2019$room, Set = anticlusters)
#'
#' # Optimize the distance criterion
#' # (tends to maximize similarity in feature means)
#' ac_dist <- anticlustering(
#'   features,
#'   K = 3,
#'   objective = "distance",
#'   categories = schaper2019$room,
#'   method = "exchange"
#' )
#' # With the distance criterion, means tend to be less similar,
#' # but standard deviations tend to be more similar:
#' by(features, ac_dist, function(x) round(colMeans(x), 2))
#' by(features, ac_dist, function(x) round(apply(x, 2, sd), 2))
#'
#'
#' ## Min-max anticlustering
#' # To make one variable dissimilar between sets, we use the
#' # function argument `iv` (= independent variable)
#'
#' anticlusters <- anticlustering(
#'   schaper2019[, 3:5],
#'   K = 2,
#'   preclustering = TRUE,
#'   iv = schaper2019[, 6]
#' )
#' # The similarity of the features in the min-max application is usually
#' # improved by activating preclustering, but this decreases the dissimilarity
#' # between the sets on the independent variable.
#'
#' # Check out feature means:
#' by(schaper2019[, 3:6], anticlusters, function(x) round(colMeans(x), 2))
#'
#'
#' @references
#'
#' Schaper, M. L., Kuhlmann, B. G., & Bayen, U. J. (2019). Metamemory
#' expectancy illusion and schema-consistent guessing in source monitoring.
#' Journal of Experimental Psychology: Learning, Memory, and Cognition,
#' 45, 470-496. http://dx.doi.org/10.1037/xlm0000602
#'
#' Schaper, M. L., Kuhlmann, B. G., & Bayen, U. J. (in press).
#' Metacognitive expectancy effects in source monitoring: Beliefs,
#' in-the-moment experiences, or both? Journal of Memory and Language.
#'
"schaper2019"
