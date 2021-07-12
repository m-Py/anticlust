
#' Visualize a cluster analysis
#'
#' @param features A data.frame or matrix representing the features that
#'     are plotted. Must have two columns.
#' @param clusters A vector representing the clustering
#' @param within_connection Boolean. Connect the elements within each
#'     clusters through lines? Useful to illustrate a graph structure.
#' @param between_connection Boolean. Connect the elements between each
#'     clusters through lines? Useful to illustrate a graph structure.
#'     (This argument only works for two clusters).
#' @param illustrate_variance Boolean. Illustrate the variance criterion
#'     in the plot?
#' @param show_axes Boolean, display values on the x and y-axis? Defaults
#'     to `FALSE`.
#' @param xlab The label for the x-axis
#' @param ylab The label for the y-axis
#' @param xlim The limits for the x-axis
#' @param ylim The limits for the y-axis
#' @param main The title of the plot
#' @param cex The size of the plotting symbols, see \code{\link{par}}
#' @param cex.axis The size of the values on the axes
#' @param cex.lab The size of the labels of the axes
#' @param lwd The width of the lines connecting elements.
#' @param lty The line type of the lines connecting elements
#'    (see \code{\link{par}}).
#' @param frame.plot a logical indicating whether a box should be drawn
#'    around the plot.
#' @param cex_centroid The size of the cluster center symbol (has an
#'    effect only if \code{illustrate_variance} is \code{TRUE})
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @importFrom grDevices palette rainbow
#' @importFrom graphics lines par plot points
#'
#' @details In most cases, the argument \code{clusters} is a vector
#'   returned by one of the functions \code{\link{anticlustering}},
#'   \code{\link{balanced_clustering}} or \code{\link{matching}}. 
#'   However, the plotting function can also be used to plot the results 
#'   of other cluster functions such as \code{\link{kmeans}}. This function
#'   is usually just used to get a fast impression of the results of an 
#'   (anti)clustering assignment, but limited in its functionality. 
#'   It is useful for depicting the intra-cluster connections using 
#'   argument \code{within_connection}.
#'
#' @examples
#'
#' N <- 15
#' features <- matrix(runif(N * 2), ncol = 2)
#' K <- 3
#' clusters <- balanced_clustering(features, K = K)
#' anticlusters <- anticlustering(features, K = K)
#' user_par <- par("mfrow")
#' par(mfrow = c(1, 2))
#' plot_clusters(features, clusters, main = "Cluster editing", within_connection = TRUE)
#' plot_clusters(features, anticlusters, main = "Anticluster editing", within_connection = TRUE)
#' par(mfrow = user_par)
#' 
plot_clusters <- function(features, clusters, within_connection = FALSE,
                          between_connection = FALSE, illustrate_variance = FALSE,
                          show_axes = FALSE, xlab = NULL, ylab = NULL,
                          xlim = NULL, ylim = NULL,
                          main = "", cex = 1.2,
                          cex.axis = 1.2, cex.lab = 1.2, lwd = 1.5,
                          lty = 2, frame.plot = FALSE, cex_centroid = 2) {
  features <- as.matrix(features)
  if (!argument_exists(xlim)) {
    xlim <- c(min(features[, 1]), max(features[, 1]))
  }
  if (!argument_exists(ylim)) {
    ylim <- c(min(features[, 2]), max(features[, 2]))
  }
  if (ncol(features) != 2)
    stop("Can only plot two features")
  if (nrow(features) != length(clusters))
    stop("There must be a cluster for each element - nrow(features) != length(clustering)")

  ## Just in case anybody passes something weird (i.e., cluster numbers
  ## do not start at 1):
  clusters <- to_numeric(clusters)
  # Replace NAs in clustering vector
  clusters[is.na(clusters)] <- seq(
    max(clusters, na.rm = TRUE) + 1, length.out = sum(is.na(clusters))
  )

  ## Specify axis labels if they cannot be assumed from data matrix or argument
  if (is.null(xlab))
    xlab <- ifelse(is.null(colnames(features)[1]), "Feature 1", colnames(features)[1])
  if (is.null(ylab))
    ylab <- ifelse(is.null(colnames(features)[2]), "Feature 2", colnames(features)[2])

  x <- features[, 1]
  y <- features[, 2]

  ## Set colors for clusters
  K <- length(unique(clusters))
  if (K < 8) {
    col <- palette()[2:(K+1)]
    pch <- c(15:18, 3, 4, 8)[1:K]
    pch <- pch[clusters]
  } else {
    col <- rainbow(K)
    pch <- 19
  }

  col <- col[clusters]
  
  # Draw plot
  axt <- "n"
  if (show_axes == TRUE)
    axt <- "s"
  plot(x, y, las = 1, cex.axis = cex.axis, cex.lab = cex.lab,
       col = col, xlab = xlab, ylab = ylab, cex = cex, bg = col,
       xaxt = axt, yaxt = axt, pch = pch, main = main,
       frame.plot = frame.plot, xlim = xlim, ylim = ylim)
  ## Draw graph structure
  if (within_connection == TRUE) {
    draw_all_cliques(x, y, clusters, cols = col, lwd = lwd, lty = lty)
  }
  if (between_connection == TRUE) {
    if (K == 2) {
      draw_between_cliques(x, y, clusters, lwd = lwd, lty = lty)
      ## redraw points so that the lines do notlay the points
      points(x, y, cex = cex, pch = pch, col = col)
    } else {
      warning("Connections between elements of different clusters were ",
              "not drawn because there were more than two clusters.")
    }
  }
  if (illustrate_variance) {
    illustrate_variance(features, clusters, col, lwd, lty, cex_centroid)
  }
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
  if (length(x) == 1) {
    return(invisible(NULL))
  }
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

## Function to illustrate distances between elements of different cliques
## (does not work in general, but only for K = 2)
draw_between_cliques <- function(x, y, assignment, lwd = 1.5,
                                 lty = 2, col = "darkgrey") {
  selfx <- x[assignment == 1] ## first cluster
  selfy <- y[assignment == 1]
  otherx <- x[assignment != 1] ## second cluster
  othery <- y[assignment != 1]
  for (j in seq_along(selfx)) {
    for (k in seq_along(otherx)) {
      lines(c(selfx[j], otherx[k]), c(selfy[j], othery[k]),
            col = col, lwd = lwd, lty = lty)
    }
  }
}

illustrate_variance <- function(features, clusters, cols, lwd, lty, cex_centroid) {
  K <- length(unique(clusters))
  centers <- cluster_centers(features, clusters)
  # only K colors; in the order implied by "clusters" so that the color
  # scheme is consistent with the points
  cols <- unique(cols[order(clusters)])
  for (i in 1:K) {
    feats <- features[clusters == i, ]
    draw_between_cliques(
      c(feats[, 1], centers[i, 1]),
      c(feats[, 2], centers[i, 2]),
      c(rep(1, nrow(feats)), 2),
      lwd = lwd,
      col = cols[i],
      lty = lty
    )
  }
  points(centers, cex = cex_centroid, pch = 24,
         col = "black", bg = cols, lwd = cex_centroid + 1)
}
