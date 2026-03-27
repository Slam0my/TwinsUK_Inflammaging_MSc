#Collating all solar outputs into one file for easier readability and accessibility 
  
#Load libraries 
library(stringr)
library(dplyr)

#Set working directory
setwd("/users/k25046756/Twins_Project/Rsolar/solar_analysis/data")

#In /results there are separate directories for each job_chunk_
job_chunks <- list.files("/users/k25046756/Twins_Project/Rsolar/solar_analysis/results",
                      pattern ="job_chunk_", full.names=TRUE)

#Save as data frame (empty to store results)
summary <- data.frame()

#Loop for obtaining all of the information in each job_chunk_
for(job in job_chunks){
  
  #Using the directory names to identify the pairs in each analysis
  pairs <- list.dirs(job, full.names=TRUE, recursive = FALSE)
  
  #Loop for each pair per directory
  for(pair_dir in pairs){
    
    #In each pair directory, finds polygeneic.out file where the SOLAR results are
    file <- file.path(pair_dir, "polygenic.out")
    
    #Reads the file to extract information
    #Check the file exists
    if(file.exists(file)){
      
      #Reads the file output
      lines <- readLines(file)
      
      pair_name <- basename(pair_dir)
      parts <- unlist(str_split(pair_name, "\\."))
      trait1 <- parts[1]         
      trait2 <- parts[length(parts)] 
      
      #Heritability extraction
      #Finding the lines
      h2_lines <- lines[grepl("H2r\\(.+\\) is [0-9.Ee+-]+", lines)]
      
      if(length(h2_lines) >= 2) {
        
        #IP trait1 heritability
        #(?<=) means lookbehind
        #\\( means in the bracket 
        #. means any character 
        #+ means one or more 
        #? stop at first match
        #.+? takes the matching
        #(?=\\) means lookahead stop before the close bracket
      

        
        h2r1 <- as.numeric(str_extract(h2_lines[1], "(?<=is )[0-9.Ee+]+")) #These takes the numerical values after h2r is _
        
        h2r1_se_line <- lines[grepl("H2r.*Std\\. Error:", lines)]
        h2r1_se <- if(length(h2r1_se_line) > 0) as.numeric(gsub(".*Std\\. Error:\\s*([0-9.Ee+-]+).*", "\\1", h2r1_se_line[1])) else NA

        #Protein trait 2 heritability
        h2r2 <- as.numeric(str_extract(h2_lines[2], "(?<=is )[0-9.Ee+-]+"))
        
        h2r2_se_line <- if(length(h2r1_se_line) > 1) h2r1_se_line[2] else NA
        h2r2_se <- if(!is.na(h2r2_se_line)) as.numeric(gsub(".*Std\\. Error:\\s*([0-9.Ee+-]+).*", "\\1", h2r2_se_line)) else NA
        
      } else next
      
      #RhoE environmental correlations
      rhoE_line <- lines[grep("RhoE is", lines)]
      
      rhoE <- if(length(rhoE_line) > 0)
        as.numeric(str_extract(rhoE_line[1], "(?<=RhoE is )-?[0-9.Ee+]+")) else NA
      
      rhoE_p <- if(length(rhoE_line) > 0)
        as.numeric(str_extract(rhoE_line[1], "(?<=p = )[0-9.Ee+]+")) else NA
      
      
      rhoE_se_line <- lines[grepl("RhoE.*Std\\. Error:", lines)]
      rhoE_se <- if(length(rhoE_se_line) > 0) as.numeric(gsub(".*Std\\. Error:\\s*([0-9.Ee+-]+).*", "\\1", rhoE_se_line[1])) else NA
      
      #RhoG genetic correlations
      rhoG_line <- lines[grep("RhoG is", lines)]
      
      rhoG <- if(length(rhoG_line) > 0)
        as.numeric(str_extract(rhoG_line[1], "(?<=RhoG is )-?[0-9.Ee+]+")) else NA
      

      rhoG_se_line <- lines[grepl("RhoG.*Std\\. Error:", lines)]
      rhoG_se <- if(length(rhoG_se_line) > 0) as.numeric(gsub(".*Std\\. Error:\\s*([0-9.Ee+-]+).*", "\\1", rhoG_se_line[1])) else NA
      
      
      rhoG_p0_line <- lines[grep("RhoG different from zero", lines)] #if p<0.05 then genetic relationship
      
      
      rhoG_p0 <- if(length(rhoG_p0_line) > 0) 
        as.numeric(str_extract(rhoG_p0_line[1], "([0-9.]+([eE][-+]?[0-9]+)?)$")) else NA
      
      rhoG_p1_line <- lines[grep("RhoG different from (1\\.0|-1\\.0)", lines)]
      rhoG_p1 <- if(length(rhoG_p1_line) > 0) {
          as.numeric(str_extract(rhoG_p1_line[1], "([0-9.]+([eE][-+]?[0-9]+)?)$"))
        } else NA
      
      
      #RhoP phenotypic correlation 
      rhoP_line <- lines[grep("Derived Estimate of RhoP is", lines)]
      
      rhoP <- if(length(rhoP_line) > 0)
        as.numeric(str_extract(rhoP_line[1], "-?[0-9.Ee+]+$")) else NA

      #Binding into one table 
      summary <- rbind(summary,
                       data.frame(trait1, trait2, h2r1, h2r1_se, h2r2, h2r2_se, rhoE, rhoE_p, rhoE_se, rhoG,rhoG_p0, rhoG_p1, rhoG_se, rhoP, stringsAsFactors = FALSE))
    }
  }
}

#Saving the results to a csv file
write.table(summary, "summary_heritability.csv", sep="\t", row.names=FALSE, quote=FALSE)

#For viewing table for all 
summary_results<- read.delim("summary_heritability.csv")


#Showing the NAs
na_rows <- which(rowSums(is.na(summary_results)) > 0)
summary_results[na_rows, ]
#There is 1 NA


#Finding if any pairs didn't collate
missing_pairs <- c()

for(job in job_chunks){
  
  pairs <- list.dirs(job, full.names=TRUE, recursive=FALSE)
  
  for(pair_dir in pairs){
    
    file <- file.path(pair_dir, "polygenic.out")
    
    if(file.exists(file)){
      
      lines <- readLines(file)
      
      h2_lines <- lines[grepl("H2r\\(.+\\) is", lines)]
      
      if(length(h2_lines) < 2){
        missing_pairs <- c(missing_pairs, pair_dir)
      }
      
    } else {
      missing_pairs <- c(missing_pairs, pair_dir)
    }
  }
}

length(missing_pairs)
missing_pairs

for(pair_dir in missing_pairs){
  file <- file.path(pair_dir, "polygenic.out")
  if(file.exists(file)){
    lines <- readLines(file)
    h2_check <- grep("H2r", lines, ignore.case = TRUE, value = TRUE)
    cat("Found", length(h2_check), "H2r lines in:", basename(pair_dir), "\n")
  }
}


