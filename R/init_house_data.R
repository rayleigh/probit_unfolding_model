library(tidyverse)
library(wnominate)

#setwd()

#Download code and data from 10.1214/21-AOAS1454SUPPB
source("circ_house_code/source/read_kh2.R")
source("circ_house_code/source/ymat_spit.R")

for (i in 100:115) {
  house_votes_info <- ymat_spit(i, T)
  house_votes_m <- house_votes_info[[1]]
  rownames(house_votes_m) <- house_votes_info[[2]]
  remove_ind <- 
    which(apply(house_votes_m, 2, function(col) {length(unique(col[!is.na(col)])) == 1}))
  house_votes_m <- house_votes_m[, -remove_ind]
  save(house_votes_m, file = paste(
    "data_files/processed_house_votes_", i, ".Rdata", sep = ""))
}

# Download Voteview votes and member info from
# https://voteview.com/data for 116th and 117th House
party_affiliation <- function(party_code) {
  
  affiliation = party_code
  
  affiliation[party_code == 100] = "D"
  affiliation[party_code == 200] = "R"
  affiliation[party_code == 328] = "I"
  
  return(affiliation)
}

translate_cast_code <- function(cast_code) {
  
  vote = rep(NA, length(cast_code))
  vote[cast_code <= 3] = 1
  vote[cast_code <= 6 & cast_code > 3] = 0
  return(vote)
}

get_non_na_val <- function(vote) {
  if (any(!is.na(vote))) {
    return(vote[!is.na(vote)])
  }
  return(NA)
}

house_116_vote_m <- read_csv("data_files/H116_votes.csv")
house_116_members <- read_csv("data_files/H116_members.csv")
house_116_members <- house_116_members %>% separate(bioname, c("last_name"), sep = ",")
house_116_members <- house_116_members %>% mutate(party_abbrev = party_affiliation(party_code))
house_116_members$member_name <- apply(house_116_members, 1, function(row) {
  paste(row["last_name"], " (",
        row["party_abbrev"], " ",
        row["state_abbrev"], "-", as.numeric(row["district_code"]), ")", sep = "")
})
house_116_member_info <- house_116_members %>% select(icpsr, member_name) 

preprocessed_house_votes_m <- 
  house_116_vote_m %>% select(rollnumber:cast_code) %>% 
  mutate(vote = translate_cast_code(cast_code)) %>%
  left_join(house_116_member_info, by = "icpsr") %>%
  select(rollnumber, vote, member_name) %>%
  pivot_wider(names_from = rollnumber, values_from = vote) 
amash_inds <- grep("AMASH", preprocessed_house_votes_m$member_name)
preprocessed_house_votes_m[amash_inds[1], -1] <-
  t(apply(preprocessed_house_votes_m[amash_inds, -1], 2, get_non_na_val))
preprocessed_house_votes_m <-
  preprocessed_house_votes_m[-amash_inds[2],]

mitchell_inds <- grep("MITCHELL", preprocessed_house_votes_m$member_name)
preprocessed_house_votes_m[mitchell_inds[1], -1] <-
  t(apply(preprocessed_house_votes_m[mitchell_inds, -1], 2, get_non_na_val))
preprocessed_house_votes_m <-
  preprocessed_house_votes_m[-mitchell_inds[2],]

absent_members <- which(rowMeans(is.na(preprocessed_house_votes_m[,-1])) >= 0.4)
preprocessed_house_votes_m <- preprocessed_house_votes_m[-absent_members,]
house_votes_m <- data.matrix(preprocessed_house_votes_m[,-1])
rownames(house_votes_m) <- preprocessed_house_votes_m$member_name
unanimous_votes <- which(apply(house_votes_m, 2, function(col) {length(table(col)) == 1}))
house_votes_m <- house_votes_m[, -unanimous_votes]
save(house_votes_m, file = "data_files/processed_house_votes_116.Rdata")

house_117_vote_m <- read_csv("data_files/H117_votes.csv")
house_117_members <- read_csv("data_files/H117_members.csv")
house_117_members <- house_117_members %>% separate(bioname, c("last_name"), sep = ",")
house_117_members <- house_117_members %>% mutate(party_abbrev = party_affiliation(party_code))
house_117_members$member_name <- apply(house_117_members, 1, function(row) {
  # print(as.numeric(row["party_code"]))
  paste(row["last_name"], " (",
        row["party_abbrev"], " ",
        row["state_abbrev"], "-", as.numeric(row["district_code"]), ")", sep = "")
})
house_117_member_info <- house_117_members %>% select(icpsr, member_name) 

preprocessed_house_votes_m <- 
  house_117_vote_m %>% select(rollnumber:cast_code) %>% 
  mutate(vote = translate_cast_code(cast_code)) %>%
  left_join(house_117_member_info, by = "icpsr") %>%
  select(rollnumber, vote, member_name) %>%
  pivot_wider(names_from = rollnumber, values_from = vote) 
absent_members <- which(rowMeans(is.na(preprocessed_house_votes_m[,-1])) >= 0.4)
preprocessed_house_votes_m <- preprocessed_house_votes_m[-absent_members,]
house_votes_m <- data.matrix(preprocessed_house_votes_m[,-1])
rownames(house_votes_m) <- preprocessed_house_votes_m$member_name
unanimous_votes <- which(apply(house_votes_m, 2, function(col) {length(table(col)) == 1}))
house_votes_m <- house_votes_m[, -unanimous_votes]
save(house_votes_m, file = "data_files/processed_house_votes_117.Rdata")


