
# test obj function for multiple split variables

oasis <- read.csv("https://raw.githubusercontent.com/aenneb/OASIS-beauty/master/means_per_image.csv")

design <- c(2, 3)
iv <- c("Valence_mean", "Arousal_mean")
covariate <- "beauty_mean"
obj_funs <- make_obj_function(oasis, covariate, iv, design)

levels <- expand.grid(lapply(design, function(x) 1:x))

cl <- balanced_clustering(oasis[, "Valence_mean"], K = prod(design))

obj_funs(cl, oasis)

obj_value_distance(cl, oasis[, "beauty_mean", drop = FALSE]) - (
  obj_value_distance(levels[, 1][cl], oasis[, "Valence_mean", drop = FALSE]) +
  obj_value_distance(levels[, 2][cl], oasis[, "Arousal_mean", drop = FALSE])
)

