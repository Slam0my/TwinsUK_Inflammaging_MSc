#Phenotype Creation

setwd("/users/k25046756/Twins_Project/Rsolar/solar_analysis/data/")


#Loading the dataset containing the information
load("/users/k25046756/Twins_Project/raw_olink_ip/data/dataset.RData")

#load("/scratch/users/k25046756/Twins_Project/raw_olink_ip/results/Olink_associations.RData")

#To make sure the ids in pedigree.csv matches the ids here (removing the singletons)
ped <- read.csv("/users/k25046756/Twins_Project/Rsolar/solar_analysis/data/twinsuk_pedigree.csv")
ped_ids <- ped$id


#Loading library
library(dplyr)

#Filtering for QC and no twinless twin
dataset_clean <- dataset %>%
  mutate(across(c(id, fid), ~as.numeric(as.character(.)))) %>%
  filter(`QC Warning` == "Pass") %>%
  group_by(fid) %>%
  filter(n() == 2) %>% 
  ungroup()

#Adding season as a covariate but before, needs numerical conversion 

dataset_final <- dataset_clean %>%
  mutate(
    season_num = case_when(
      season == "Win" ~ 1,
      season == "Spr" ~ 2,
      season == "Sum" ~ 3,
      season == "Aut" ~ 4,
      TRUE ~ NA_real_
    ),
    plate_num = as.numeric(as.factor(`Plate ID`))
  )

#Use columns numbers for all IP and protein
phenotype_solar <- dataset_final %>%
  filter(id %in% ped_ids) %>%  
  select(id, Age, BMI, season_num, plate_num , Freezing_time, 24:42393)


#Inverse normal distribution 
#Only on the traits columns

trait_cols <- names(phenotype_solar)[7:ncol(phenotype_solar)]

phenotype_solar[trait_cols] <- lapply(phenotype_solar[trait_cols], function(x) {
  #Skipping NA values
  if(all(is.na(x))) return(x)
  
  #Blom style 
  n <- sum(!is.na(x))
  r <- rank(x, na.last = "keep", ties.method = "random")
  qnorm((r - 0.5) / n)
})

#need to change formatting for solar
#Note: need to be careful as name mismatches in the dataset and tr data frames with this phenotype csv
colnames(phenotype_solar) <- gsub("-| ", "_", colnames(phenotype_solar))

write.csv(phenotype_solar, "twinsuk_phenotypes.csv", row.names = FALSE, na = "", quote = FALSE)
