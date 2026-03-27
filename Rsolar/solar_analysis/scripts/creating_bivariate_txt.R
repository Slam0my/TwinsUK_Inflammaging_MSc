#Generating text files that contains all the solar and traits to be submitted all at once
#By using the pair list, loop the solar command texts to run through all 

#Import the pair list to be used 
#18629 pairs
pairs <- read.csv("/users/k25046756/Twins_Project/Rsolar/solar_analysis/data/pair_list.csv")  

#As this is a heavy job, split the pairs into chunks to run parallel jobs
#Trying 200 pairs per chunk
chunks<- 200
chunk_size <- ceiling(nrow(pairs) / chunks)

chunk_id <- 1

for (start in seq(1, nrow(pairs), by = chunk_size)) {
  
  end <- min(start + chunk_size - 1, nrow(pairs))
  
  
  #Selecting the pairs for the chunks
 # start <- (chunk - 1) * chunk_size + 1
  #end   <- min(chunk * chunk_size, nrow(pairs))
  
  chunk_pairs <- pairs[start:end, ]
  
  batch_file <- paste0("/users/k25046756/Twins_Project/Rsolar/solar_analysis/scripts/solar_chunks/bivariate_chunk_"
                       , chunk_id, ".txt")

#This is what will be written for SOLAR at the beginning
cat("load pedigree /users/k25046756/Twins_Project/Rsolar/solar_analysis/data/twinsuk_pedigree.csv\n", file=batch_file)
cat("load phenotypes /users/k25046756/Twins_Project/Rsolar/solar_analysis/data/twinsuk_phenotypes.csv\n", file=batch_file, append=TRUE)

#Looping the model commands with all pairs as inputs

for(i in 1:nrow(chunk_pairs)){
  
  trait1 <- chunk_pairs$IP[i]
  trait2 <- chunk_pairs$Protein[i]
  
  #Need model new to  let SOLAR know that it is separate for each pair
  if(i > 1){
    cat("model new\n", file=batch_file, append=TRUE) 
  }

  #Testing the two new pair traits
  cat(sprintf("trait %s %s\n", trait1, trait2), file=batch_file, append=TRUE)
  
  #Setting covariates 
  cat("covariate Age BMI season_num plate_num Freezing_time\n", file=batch_file, append=TRUE)
  
  #For bivariate analysis, need -testrhog for genetic and -testrhoe for environmental 
  cat("polygenic -testrhog -testrhoe\n\n", file=batch_file, append=TRUE)
}
chunk_id <- chunk_id + 1
}
