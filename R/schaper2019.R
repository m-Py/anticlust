
#' Consistency rating for 96 German words
#'
#' A dataset containing average consistency and inconsistency ratings,
#' number of syllables and frequency for 96 objects.
#'
#' @format A data frame with 96 rows and 7 variables
#' \describe{
#'   \item{Item}{The object (in German language)}
#'   \item{Room}{1 = the item is typically found in a kitchen;
#'     2 = the item is typically found in a bathroom}
#'   \item{ConsistencyRating}{How expected would it
#'      be to find the \code{Item} in the typical \code{Room}}
#'   \item{InconsistencyRating}{How expected would it
#'      be to find the \code{Item} in the atypical \code{Room}}
#'   \item{Syllables}{The number of syllables of the word}
#'   \item{Frequency}{A value indicating the relative frequency of the word}
#'   \item{List}{Represents the set affiliation of the \code{Item} as
#'       realized in experiments by Schaper, Kuhlmann and Bayen
#'       (2019; in press)}
#' }
#'
#' @source
#' Courteously provided by Marie Lusia Schaper and Ute Bayen.
#'
#' @examples
#'
#' features <- schaper2019[, 3:6]
#' head(features)
#' # This command scales all features to have a similar range:
#' features <- apply(features, 2, function(x) x / max(x))
#' head(features)
#'
#' # Optimize the variance criterion
#' # (tends to maximize similarity in feature means)
#' anticlusters <- anticlustering(
#'   features,
#'   K = 3,
#'   objective = "variance",
#'   parallelize = TRUE,
#'   categories = schaper2019$Room,
#'   seed = 451
#' )
#'
#' # Means are quite similar across sets:
#' by(schaper2019[, 3:6], anticlusters, function(x) round(colMeans(x), 2))
#' # Room is balanced between the three sets:
#' table(Room = schaper2019$Room, Set = anticlusters)
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
