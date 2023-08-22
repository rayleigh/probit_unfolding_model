library(tidyverse)

#setwd()

#Download data from https://mqscores.lsa.umich.edu/replication.php
load("data_files/mqData2021.Rda")

mqVotes <- t(mqData[, 1:(ncol(mqData) - 3)])
mqTime <- mqData$time
pos_inds <- c(39)
neg_inds <- c(12, 29, 9)
neg_year_inds <- c(list(1:29), list(1:24), list(1))
save(mqVotes, mqTime, pos_inds, neg_inds, neg_year_inds,
     file = "data_files/mq_supreme_court_vote_info_2021.Rdata")


