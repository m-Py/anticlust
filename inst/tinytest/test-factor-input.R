
library("anticlust")
library("tinytest")

# Test that using factor variables as input now works

gr <- anticlustering(iris, K = 3, standardize = TRUE)
mean_sd_tab(iris[-5], gr)
table(gr, iris$Species)

gr <- anticlustering(iris, K = 3, standardize = TRUE, method = "3phase")
mean_sd_tab(iris[-5], gr)
table(gr, iris$Species)

gr <- anticlustering(iris, K = 5, standardize = TRUE, method = "3phase", objective = "kplus")
mean_sd_tab(iris[-5], gr)
table(gr, iris$Species)

gr <- anticlustering(iris$Species, K = 5, method = "3phase",  objective = "variance")
table(gr, iris$Species)
expect_true(all(table(gr, iris$Species) == 10))
