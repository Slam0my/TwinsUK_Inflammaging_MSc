#Creating the pedigree file

setwd("/users/k25046756/Twins_Project/Rsolar/solar_analysis/data/")


#Loading the dataset containing the information
load("/users/k25046756/Twins_Project/raw_olink_ip/data/dataset.RData")

#Loading library
library(dplyr)


dataset_clean <- dataset %>%
  mutate(across(c(id, fid), ~as.numeric(as.character(.)))) %>%
  filter(`QC Warning` == "Pass") 

#Creating a new data frame with chosen columns
pedigree_solar <-dataset_clean %>%
  group_by(fid) %>%
  filter(n() == 2) %>% #No singletons
  mutate(mztwin = if_else(first(Z) == "MZ", first(fid), 0)) %>%
  ungroup() %>%
  mutate(
    id = id,
    fa = 100000 + fid, #father
    mo = 200000 + fid, #mother
    sex = 2,     #All females
    famid = fid, #Family ID used to link twins 
  ) %>%
  #Reordering the columns 
  select(id, fa, mo, sex, famid, mztwin)


write.csv(pedigree_solar, "twinsuk_pedigree.csv", row.names = FALSE, quote = FALSE)

#table(table(pedigree_solar$famid))
#table(pedigree_solar$mztwin)


