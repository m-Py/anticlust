% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objective-diversity.R
\name{diversity_objective}
\alias{diversity_objective}
\title{(Anti)cluster editing "diversity" objective}
\usage{
diversity_objective(x, clusters)
}
\arguments{
\item{x}{The data input. Can be one of two structures: (1) A data matrix
where rows correspond to elements and columns correspond to
features (a single numeric feature can be passed as a vector). (2)
An N x N matrix dissimilarity matrix; can be an object of class
\code{dist} (e.g., returned by \code{\link{dist}} or
\code{\link{as.dist}}) or a \code{matrix} where the entries of
the upper and lower triangular matrix represent the pairwise
dissimilarities.}

\item{clusters}{A vector representing (anti)clusters (e.g.,
returned by \code{\link{anticlustering}}).}
}
\value{
The cluster editing objective
}
\description{
Compute the diversity for a given clustering.
}
\details{
The objective function used in (anti)cluster editing is the
diversity, i.e., the sum of the pairwise distances between elements
within the same groups. When the input \code{x} is a feature
matrix, the Euclidean distance is computed as the basic distance
unit of this objective.
}
\examples{

data(iris)
distances <- dist(iris[1:60, -5])
## Clustering
clusters <- balanced_clustering(distances, K = 3)
# This is low:
diversity_objective(distances, clusters)
## Anticlustering
anticlusters <- anticlustering(distances, K = 3)
# This is higher:
diversity_objective(distances, anticlusters)

}
\references{
Brusco, M. J., Cradit, J. D., & Steinley, D. (2020). Combining
diversity and dispersion criteria for anticlustering: A bicriterion
approach. British Journal of Mathematical and Statistical
Psychology, 73, 275-396. https://doi.org/10.1111/bmsp.12186

Papenberg, M., & Klau, G. W. (2021). Using anticlustering to partition 
data sets into equivalent parts. Psychological Methods, 26(2), 
161–174. https://doi.org/10.1037/met0000301.
}
\author{
Martin Papenberg \email{martin.papenberg@hhu.de}
}
