library("anticlust")


# Test feature input:
data <- schaper2019[, 2:6]
data[1:3,1] <- NA
tt <- anticlustering(data, K = 3, standardize = TRUE, objective = "kplus")
mean_sd_tab(data[-1], tt)
table(tt, data$room, useNA = "ifany")

# Test distance input
dists <- dist(anticlust:::get_anticlustering_features(data, standardize = TRUE, objective = "kplus"))
sum(is.na(dists))
tt <- anticlustering(dists, K = 3)
mean_sd_tab(data[-1], tt)
table(tt, data$room, useNA = "ifany")

# Test distance input (but NAs in distances, which is not allowed!)
dists <- dist(schaper2019[, 3:6])
dists[1] <- NA
expect_error(
  anticlustering(dists, K = 3), 
  pattern = "No NA allowed in distance matrix input."
)

# Test distance input (but NAs in distances in matrix form, which is not allowed!)
dists <- dist(schaper2019[, 3:6])
dists[1] <- NA
expect_error(
  anticlustering(as.matrix(dists), K = 3), 
  pattern = "No NA allowed in distance matrix input."
)

