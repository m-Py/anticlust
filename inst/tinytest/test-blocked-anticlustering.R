
library(anticlust)
library(tinytest)


table(brunel2025$target_word_emotionality)

blocks <- brunel2025$target_word_emotionality

features <- brunel2025[, c("valence_target_word", "arousal_target_word")]

groups <- anticlustering(features, K = 2, objective = "kplus", blocks = blocks, standardize = TRUE, method = "exchange")
tab <- table(blocks, groups)

# is condition blocked?
expect_true(all(abs(tab[ ,1] - tab[, 2]) <= 1))

# non blocked assignment
groups2 <- anticlustering(features, K = 2, objective = "kplus", categories = blocks, standardize = TRUE, method = "exchange")

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


## Test wrapper function anticlustering() with different input specifications
anticlustering(features, K = 2, objective = "diversity", blocks = blocks, standardize = FALSE)

anticlustering(features, K = 2, objective = "diversity", blocks = blocks, method = "local-maximum")

anticlustering(features, K = 2, objective = "diversity", blocks = blocks, method = "local-maximum", repetitions = 10)

groups <- anticlustering(features, K = 2, objective = "kplus", blocks = blocks, method = "local-maximum", repetitions = 10, preclustering = TRUE)

foo <- anticlustering(features, K = 2, objective = "kplus", blocks = blocks, method ="3phase")

tab <- table(blocks, foo)
expect_true(all(abs(tab[ ,1] - tab[, 2]) <= 1))

anticlustering(features, K = 2, objective = "kplus", blocks = blocks, categories = brunel2025$sentence_emotionality)

## use other (random) data input

N <- 1000
blocks <- sample(1:10, size = N, replace = TRUE)
foo <- anticlustering(1:1000, K = 10, blocks = blocks, objective = "variance")

tab <- table(blocks, foo)
expect_true(all(abs(tab[ ,1] - tab[, 2]) <= 1))

