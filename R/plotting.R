
#' Visualize a clustering of two features
#'
#' @param features A data.frame or matrix representing the features that
#'     are used. Must have two columns.
#' @param clustering A vector representing the clustering
#' @param connect_clusters Boolean. Connect the groups through lines?
#'     Useful to illustrate a graph structure.
#' @param show_axes Boolean, display values on the x and y-axis? Defaults
#'     to `FALSE`.
#' @param xlab The label for the x-axis
#' @param ylab The label for the y-axis
#' @param cols The coloring of the groups. Does not need to be passed.
#'     If it is passed, it needs to be a character vector of the same
#'     length as there are clusters, and each element of the vector is a
#'     color.
#' @param pch A numeric vector representing the symbol used to plot
#'     data points (see \code{\link{par}}). Should either be of length 1,
#'     then the same symbol is used for all plots, or should have the
#'     same length as there are different clusters. In the latter case,
#'     each cluster has a different plotting symbol.
#' @param main The title of the plot
#' @param cex The size of the plotting symbols, see \code{\link{par}}
#' @param cex.axis The size of the values on the axes
#' @param cex.axis The size of the labels of the axes
#'
#' @export
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics lines par plot
#'
#' @details In most cases, the argument `clustering` is a vector
#'   returned by one of the functions \code{\link{anticlustering}} or
#'   \code{\link{clustering}}. However, the plotting function can also
#'   be used to plot the results of other cluster functions such as
#'   \code{\link{kmeans}}.
#'
#' @examples
#'
#' n_elements <- 90
#' features <- matrix(runif(n_elements * 2), ncol = 2)
#' n_groups <- 3
#' clusters <- clustering(features, n_groups)
#' anticlusters <- anticlustering(features, n_groups, method = "sampling")
#' par(mfrow = c(1, 2))
#' plot_clusters(features, clusters, pch = c(15, 16, 17))
#' plot_clusters(features, anticlusters, pch = c(15, 16, 17))
#'
plot_clusters <- function(features, clustering, connect_clusters = FALSE,
                          show_axes = FALSE, xlab = NULL, ylab = NULL,
                          col = NULL, pch = 19, main = "", cex = 1.2,
                          cex.axis = 1.2, cex.lab = 1.2) {
  if (ncol(features) != 2)
    stop("Can only plot two features")
  if (nrow(features) != length(clustering))
    stop("There must be a cluster for each element - nrow(features) != length(clustering)")

  ## Just in case anybody passes something weird (i.e., cluster numbers
  ## do not start at 1):
  clustering <- as.numeric(factor(clustering))

  ## More input handling:
  n_cliques <- length(unique(clustering))
  if (length(pch) == n_cliques)
    pch <- pch[clustering]
  else if (length(pch) == 1)
    pch <- pch # :-)
  else
    stop("Length of argument pch must be 1 the number of different anticlusters")


  ## Specify axis labels if they cannot be assumed from data matrix or argument
  if (is.null(xlab))
    xlab <- ifelse(is.null(colnames(features)[1]), "Feature 1", colnames(features)[1])
  if (is.null(ylab))
    ylab <- ifelse(is.null(colnames(features)[2]), "Feature 2", colnames(features)[2])

  x <- features[, 1]
  y <- features[, 2]
  if (is.null(col) == TRUE)
    col <- rainbow(n_cliques)
  ## Set colors for cliques
  if (length(col) != n_cliques)
    stop("If you specify colors, you have to specify as many colors as there are clusters")
  col <- col[clustering]
  def_mar <- c(5, 4, 4, 2) + 0.1
  par(mar = def_mar + c(0, 1, 0, 0))
  axt <- "n"
  if (show_axes == TRUE)
    axt <- "s"
  plot(x, y, las = 1, cex.axis = cex.axis, cex.lab = cex.lab,
       col = col, xlab = xlab, ylab = ylab, cex = cex,
       xaxt = axt, yaxt = axt, pch = pch, main = main)
  if (connect_clusters)
    draw_all_cliques(x, y, clustering, cols = col)
  par(mar = def_mar)
}

# Note: length(cols) must be == length(x)
draw_all_cliques <- function(x, y, assignment, lwd = 1.5, lty = 1, cols) {
  n_cliques <- length(unique(assignment))
  for (i in 1:n_cliques) {
    draw_clique(x[assignment == i], y[assignment == i],
                col = cols[assignment == i], lwd = lwd,
                lty = lty)
  }
}

# draw edges connecting all elements of a clique
draw_clique <- function(x, y, col, lwd = 1.5, lty = 1) {
  connections <- connections_within_clique(x, y)
  cords <- data.frame(x = x, y = y)
  cords <- cords[c(t(connections)), ]
  startpoints <- seq(1, nrow(cords) - 1, 2)
  for (i in startpoints) {
    start <- i
    end   <- i + 1
    lines(c(cords$x[start], cords$x[end]),
          c(cords$y[start], cords$y[end]),
            col = col, lwd = lwd, lty = lty)
  }
}

connections_within_clique <- function(x, y) {
  nitems <- length(x)
  connections <- expand.grid(1:nitems, 1:nitems)
  connections <- connections[connections[, 1] < connections[, 2], ]
  return(connections)
}
