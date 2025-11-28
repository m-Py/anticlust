
library(anticlust)
library(tinytest)


table(brunel2025$target_word_emotionality)

blocks <- brunel2025$target_word_emotionality

features <- brunel2025[, c("valence_target_word", "arousal_target_word")]

groups <- anticlustering(features, K = 2, objective = "kplus", blocks = blocks, standardize = TRUE, method = "local-maximum")
tab <- table(blocks, groups)

# is condition blocked?
expect_true(all(abs(tab[ ,1] - tab[, 2]) <= 1))

# non blocked assignment
groups2 <- anticlustering(features, K = 2, objective = "kplus", categories = blocks, standardize = TRUE, method = "local-maximum")

# verify balance, overall and within blocks, and compare to non-blocked assignment:

# overall
mean_sd_tab(features, groups)
mean_sd_tab(features, groups2)

# within blocks
mean_sd_tab(features[blocks == "Neutral", ], groups[blocks == "Neutral"])
mean_sd_tab(features[blocks == "Neutral", ], groups2[blocks == "Neutral"])

mean_sd_tab(features[blocks == "Positive", ], groups[blocks == "Positive"])
mean_sd_tab(features[blocks == "Positive", ], groups2[blocks == "Positive"])

mean_sd_tab(features[blocks == "Negative", ], groups[blocks == "Negative"])
mean_sd_tab(features[blocks == "Negative", ], groups2[blocks == "Negative"])


# useful results with distance matrix?
groups <- anticlustering(dist(features), K = 2, blocks = blocks, method = "local-maximum")
mean_sd_tab(features[blocks == "Negative", ], groups[blocks == "Negative"])
groups <- anticlustering(dist(features)^2, K = 2, blocks = blocks, method = "local-maximum")
mean_sd_tab(features[blocks == "Negative", ], groups[blocks == "Negative"])
# yep, seems to work with distance input

## Test wrapper function anticlustering() with different input specifications
anticlustering(features, K = 2, objective = "diversity", blocks = blocks, standardize = FALSE)

anticlustering(features, K = 2, objective = "diversity", blocks = blocks)

anticlustering(features, K = 2, objective = "variance", blocks = blocks)

anticlustering(features, K = 2, objective = "kplus", blocks = blocks, categories = brunel2025$sentence_emotionality)

## use other (random) data input

N <- 1000
blocks <- sample(1:10, size = N, replace = TRUE)
foo <- anticlustering(1:1000, K = 10, blocks = blocks, objective = "variance")

tab <- table(blocks, foo)
expect_true(all(abs(tab[ ,1] - tab[, 2]) <= 1))



## verify that some combinations of arguments throw errors
N <- 10
K <- 2
data <- rnorm(N)
must_link <- c(1, 1, 1, 2, 2, 2, 3, 4, 5, 6)
blocks <- rep(1:K, each = N/K)
expect_error(
  anticlustering(data, K = 2, blocks = blocks, must_link = must_link),
  pattern = "must-link"
)

expect_error(
  anticlustering(data, K = 2, objective = "kplus", blocks = blocks, method ="3phase"), 
  pattern = "exchange"
)
expect_error(
  anticlustering(data, K = 2, objective = "kplus", blocks = blocks, preclustering = TRUE),
  pattern = "preclustering"
)

expect_error(
  anticlustering(data, K = 2, objective = "kplus", blocks = blocks, repetitions = 10), 
  pattern = "repetitions"
)

# use cannot-link constraints
N <- 100
data <- rnorm(N)
blocks <- rep(1:K, each = N/K)
cannot_link <- cbind(1, 2)
tt <- anticlustering(data, K = 2, blocks = blocks, cannot_link = cannot_link, objective = "variance")
expect_true(tt[1] != tt[2])
tapply(data, tt, mean) |> round(2)


cannot_link <- cbind(c(1, 2, 3), c(2, 3, 1))
tt <- anticlustering(data, K = 3, blocks = blocks, cannot_link = cannot_link, objective = "variance")
expect_true(tt[1] != tt[2])
expect_true(tt[2] != tt[3])
expect_true(tt[1] != tt[3])
tapply(data, tt, mean) |> round(2)
