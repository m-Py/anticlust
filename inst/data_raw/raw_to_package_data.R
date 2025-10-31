
# See COPYRIGHT file for source of the data, and LICENSE file.

library(anticlust)
library(readxl)
library(usethis)

tt <- read_excel("./inst/data_raw/Indices_PhraseEM.xlsx")

tt$target_word_emotionality <- factor(tt$`Word Emotionality`)
tt$sentence_emotionality <- factor(tt$`Sentence Emotionality`)
tt$percentage_target_word <- tt$`% of occurrence target word`
tt$valence_target_word <- tt$`Valence Target Word`
tt$arousal_target_word <- tt$`Arousal Target Word`
tt$target_word <- tt$`Target Word`
tt$sentence <- tt$Sentence

brunel2025 <- subset(as.data.frame(tt), select = c(target_word, sentence, target_word_emotionality:arousal_target_word))

usethis::use_data(brunel2025, overwrite = TRUE)
