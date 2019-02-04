
#' Visualize a clustering of two features
#'
#' @param features A data.frame or matrix representing the features that
#'     are used. Must have two columns.
#' @param clustering A vector representing the clustering
#' @param connect_clusters Boolean. Connect the groups through lines?
#'     Useful to illustrate a graph structure.
#' @param xlab The label for the x-axis
#' @param ylab The label for the y-axis
#' @param cols The coloring of the groups. Does not need to be passed.
#'     If it is passed, it needs to be a character vector of the same
#'     length as there are groups, and each element of the vector is a
#'     color.
#' @param pch The symbol used for data points, see ?par
#' @param main The title of the plot
#' @param cex The size of the plotting symbols, see ?par
#' @param cex.axis The size of the values on the axes
#' @param cex.axis The size of the labels of the axes
#'
#' @export
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics lines par plot
#'
plot_clusters <- function(features, clustering, connect_clusters = FALSE,
                          xlab = NULL, ylab = NULL, cols = NULL,
                          pch = 19, main = "", cex = 1.2,
                          cex.axis = 1.2, cex.lab = 1.2) {

  if (ncol(features) != 2)
    stop("Can only plot two features")

  ## Specify axis labels if they cannot be assumed from data matrix or argument
  if (is.null(xlab))
    xlab <- ifelse(is.null(colnames(features)[1]), "Feature 1", colnames(features)[1])
  if (is.null(ylab))
    ylab <- ifelse(is.null(colnames(features)[2]), "Feature 2", colnames(features)[2])

  x <- features[, 1]
  y <- features[, 2]
  n_cliques <- length(unique(clustering))
  if (is.null(cols) == TRUE)
    cols <- rainbow(n_cliques)
  ## Set colors for cliques
  if (length(cols) != n_cliques)
    stop("If you specify colors, you have to specify as many colors as there are clusters")
  cols <- cols[clustering]
  def_mar <- c(5, 4, 4, 2) + 0.1
  par(mar = def_mar + c(0, 1, 0, 0))
  plot(x, y, las = 1, cex.axis = cex.axis, cex.lab = cex.lab,
       col = cols, xlab = xlab, ylab = ylab, cex = cex,
       xaxt = "n", yaxt = "n", pch = pch, main = main)
  if (connect_clusters)
    draw_all_cliques(x, y, clustering, cols = cols)
  par(mar = def_mar)
}

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
