#Generating a list containing all of the IP pairs to be used for automation within the SOLAR heritability analysis.
#For this, we only want to use those that were genetically significant 0.1 h2r hits (to save computational and analysis time)

setwd("/users/k25046756/Twins_Project/Rsolar/solar_analysis/scripts/")

#Loading the saved data frame that contains all the IP  names 
load("/users/k25046756/Twins_Project/Rsolar/solar_analysis/scripts/solar_analysis_filtering_results.RData")

#Library
library(dplyr)

#only IPs
unique_ips <- unique(genetic_0.1_hits$trait1)
unique_ips <- gsub("-| ", "_", unique_ips)

pair_matrix <- combn(unique_ips, 2)
pair_list <- as.data.frame(t(pair_matrix))

colnames(pair_list) <- c("IP1", "IP2")
# Check the number now
nrow(pair_list)



#Saving it as a csv
write.csv(pair_list, "/users/k25046756/Twins_Project/Rsolar/ip_ip_solar/gen_0.1/data/gen_pair_list.csv", row.names = FALSE, quote = FALSE)

nrow(pair_list)
