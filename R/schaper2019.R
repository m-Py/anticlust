
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
#'      be to find the item in the typical room
#'      (see\code{Room})}
#'   \item{InconsistencyRating}{How expected would it
#'      be to find the item in the atypical room
#'      (see\code{Room})}
#'   \item{Syllables}{The number of syllables}
#'   \item{Frequency}{A value indicating the relative frequency of a word}
#'   \item{List}{The group the item was affiliated with.}
#' }
#'
#' @source
#' Courteously provided by Marie Lusia Schaper and Ute Bayen.
#'
#' @examples
#'
#' features <- schaper2019[, 3:6]
#' # This scales all features to have a similar range in their values:
#' features <- apply(features, 2, function(x) x / max(x))
#'
#' # Optimize the variance criterion
#' start <- Sys.time()
#' anticlusters <- anticlustering(
#'   features,
#'   K = 3,
#'   objective = "variance",
#'   parallelize = TRUE,
#'   categories = schaper2019$Room,
#'   seed = 451
#' )
#' Sys.time() - start
#' # Room is balanced between the three sets:
#' table(schaper2019$Room, anticlusters)
#' # Means are quite similar across sets:
#' by(schaper2019[, 3:6], anticlusters, function(x) round(colMeans(x), 2))
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
