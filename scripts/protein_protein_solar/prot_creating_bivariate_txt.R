#Generating text files that contains all the solar and traits to be submitted all at once
#By using the pair list, loop the solar command texts to run through all 

#Import the pair list to be used 
pairs <- read.csv("/users/k25046756/Twins_Project/data/protein_protein_solar/prot_pair_list.csv")  


#Make a directory for our chunks
dir.create("/users/k25046756/Twins_Project/scripts/protein_protein_solar/solar_chunks", recursive = TRUE, showWarnings = FALSE)

#As this is a heavy job, split the pairs into chunks to run parallel jobs
#Trying 200 chunks
chunks<- 200
chunk_size <- ceiling(nrow(pairs) / chunks)

chunk_id <- 1

for (start in seq(1, nrow(pairs), by = chunk_size)) {
  
  end <- min(start + chunk_size - 1, nrow(pairs))
  
  
  #Selecting the pairs for the chunks
  # start <- (chunk - 1) * chunk_size + 1
  #end   <- min(chunk * chunk_size, nrow(pairs))
  
  chunk_pairs <- pairs[start:end, ]
  
  batch_file <- paste0("/users/k25046756/Twins_Project/scripts/protein_protein_solar/solar_chunks/bivariate_chunk_"
                       , chunk_id, ".txt")
  
  #This is what will be written for SOLAR at the beginning
  cat("load pedigree /users/k25046756/Twins_Project/data/solar_analysis/twinsuk_pedigree.csv\n", file=batch_file)
  cat("load phenotypes /users/k25046756/Twins_Project/data/solar_analysis/twinsuk_phenotypes.csv\n", file=batch_file, append=TRUE)
  
  #Looping the model commands with all pairs as inputs
  
  for(i in 1:nrow(chunk_pairs)){
    
    trait1 <- chunk_pairs$Protein1[i]
    trait2 <- chunk_pairs$Protein2[i]
    
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
