
#' Ratings for 96 words
#'
#' A stimulus set that was used in experiments by Schaper, Kuhlmann and
#' Bayen (2019a; 2019b). The item pool consists of 96 German words. 
#' Each word represents an object that is either 
#' typically found in a bathroom or in a kitchen. 
#'
#' @format A data frame with 96 rows and 7 variables
#' \describe{
#'   \item{item}{The name of an object (in German)}
#'   \item{room}{The room in which the item is typically found; can be
#'       'kitchen' or 'bathroom'}
#'   \item{rating_consistent}{How expected would it
#'       be to find the \code{item} in the typical \code{room}}
#'   \item{rating_inconsistent}{How expected would it
#'       be to find the \code{item} in the atypical \code{room}}
#'   \item{syllables}{The number of syllables in the object name}
#'   \item{frequency}{A value indicating the relative frequency of the
#'       object name in German language (lower values indicate higher
#'       frequency)}
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
#' # Maximize the diversity criterion
#' ac_dist <- anticlustering(
#'   features,
#'   K = 3,
#'   objective = "diversity",
#'   categories = schaper2019$room,
#'   method = "exchange"
#' )
#' # With the distance criterion, means tend to be less similar,
#' # but standard deviations tend to be more similar:
#' by(features, ac_dist, function(x) round(colMeans(x), 2))
#' by(features, ac_dist, function(x) round(apply(x, 2, sd), 2))
#'
#'
#' @references
#'
#' Schaper, M. L., Kuhlmann, B. G., & Bayen, U. J. (2019a). Metacognitive expectancy
#' effects in source monitoring: Beliefs, in-the-moment experiences, or both? Journal
#' of Memory and Language, 107, 95–110. https://doi.org/10.1016/j.jml.2019.03.009
#'
#' Schaper, M. L., Kuhlmann, B. G., & Bayen, U. J. (2019b). Metamemory expectancy illusion
#' and schema-consistent guessing in source monitoring. Journal of Experimental
#' Psychology: Learning, Memory, and Cognition, 45, 470.
#' https://doi.org/10.1037/xlm0000602
#'
#'
"schaper2019"
