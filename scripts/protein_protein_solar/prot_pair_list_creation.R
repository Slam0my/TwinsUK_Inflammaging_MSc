#Generating a list containing all of the protein-protein pairs to be used for automation within the SOLAR heritability analysis.

setwd("/users/k25046756/Twins_Project/scripts/protein_protein_solar/")

#Loading the tr data frame that contains all  the Protein names
load("/users/k25046756/Twins_Project/data/raw_olink_ip/results/Olink_associations.RData")

#Library
library(dplyr)


proteins <- tr %>%
  pull(Protein) %>%
  unique() %>%
  gsub("-| ", "_", .)

#Make sure no overlapping/redundancy of a x b , b x a, a x a
pair_matrix <- combn(proteins, 2)
pair_list <- as.data.frame(t(pair_matrix))

colnames(pair_list) <- c("Protein1", "Protein2")

#Saving it as a csv
write.csv(pair_list, "/users/k25046756/Twins_Project/data/protein_protein_solar/prot_pair_list.csv", row.names = FALSE, quote = FALSE)

#should have 2080 pairs as (65 x 64)/2