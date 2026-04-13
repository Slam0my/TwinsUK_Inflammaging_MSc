#libraries
.libPaths(c("~/Rlibs/4.3",
            "/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/r-4.3.0-xhenh7qr6yqtm5wvocw23qtch4vcq3ap/rlib/R/library"))

library(plyr) #data manipulation (older dplyr)
library(dplyr) #data manipulation (modern)
library(tidyr) #part of tidyverse
library(purrr) #functional programming with map()
#library(tidyverse)
library(lme4) #linear mixed effects model for modelling hierarchical data
library(magrittr) #pipe %>%
library(lubridate) # date handling

#command line arguments
#capture only custom arguments that pass to a script
#access inputs directly via index without system flags
#just gets the protein name without the messages
args <- commandArgs(trailingOnly = T)
print(args) #prints what protein is running

#runs the external r scripts
source('~/Twins_Project/scripts/lmer_testing/copy_helper_functions/copy.data.normalisation.R')
source('~/Twins_Project/scripts/lmer_testing/copy_helper_functions/copy.return.lmer.coef.R')
#load the main data frame
load('~/Twins_Project/data/raw_olink_ip/data/dataset.RData')

#protein passed as command line argument
prot=args[1]
#empty data frame to store results for all IP against the protein
sum.table = data.frame()

#loop all columns starting with "IMMUNO_" (each IP)
#For each immune cell subset/IP, run mixed models predicted protein using cell, adjusting for covariates
#loop through all immune cell types
for ( ip in colnames(dataset)[grep("IMMUNO_", colnames(dataset))]) {

  #temporary dataset for each immune cell + protein
  tmp <- dataset %>%
    #select certain columns
        dplyr::select(iid, all_of(c(prot, ip)), Age, BMI, fid) %>%
    #remove extreme outliers from protein (>4SD)
    mutate_at(vars(all_of(c(prot))), ~remove.outliers(.,4)) %>%
    #misses NA values
    na.omit() %>%
    #renames columns proteins IL and immune cell IP
    dplyr::rename(c("IL" = all_of(prot), "IP" = all_of(ip)))

  N = nrow(tmp) #number of valid samples after filtering

  #rank based inverse normal transformation?
  #formula models
  #(2|fid) = random intercept per twin pair - correlation within twins
  frm1 <- formula(paste("inverse.normal(IL) ~ Age + BMI + IP + (1|fid)"))
  frm2 <- formula(paste("inverse.normal(IL) ~ Age + IP + (1|fid)"))
  frm3 <- formula(paste("inverse.normal(IL) ~ IP + (1|fid)"))

  #runs the linear mixed models using lmer
  #fixed effects(beta,se,p) and random effects
  res1 <- lmer(frm1, data=tmp)
  res2 <- lmer(frm2, data=tmp)
  res3 <- lmer(frm3, data=tmp)

  #extract coefficients for IP
  #extracts beta, se, p for IP
  #beta effect size
  #unlist converts to numeric vector
  r1 = unlist(return.coeff.lmer(res1, 'IP'))
  r2 = unlist(return.coeff.lmer(res2, 'IP'))
  r3 = unlist(return.coeff.lmer(res3, 'IP'))

  #store results in table summary
  #results table of all proteins and immune cells associations
  sum.table = rbind(sum.table, data.frame(Protein=prot,
                                          IP=ip,
                                          n=N,
                                          beta=r1[1],
                                          se=r1[2],
                                          p=r1[3],
                                          beta1=r2[1],
                                          se1=r2[2],
                                          p1=r2[3],
                                          beta0=r3[1],
                                          se0=r3[2],
                                          p0=r3[3])
  )
}

#Each protein gets own results file
write.csv(sum.table, paste0("~/Twins_Project/data/lmer_testing/associations_test/", prot, "_ips_associations.csv"), row.names = F)
