require(lme4)

# Return regression coefficients from a lmer object

# for the twins 
#defines function. X is fitted model, var is list of variable names to extract (age, bmi, ip)
return.coeff.lmer <- function(x, var) {

  # empty to store results for each variable
 l = list()
  #if failed, fills in with NA and moves onto next... doesn't crash
 if (class(x)[1] != "lmerMod") { 

  for(i in 1:length(var)) {
   l[[i]] <- c(NA, NA, NA)
  }

 } else { 
  #extracts variance covariance matrix
  vc <- vcov(x, useScale = FALSE)
  #fixed effects(beta coefficients)
  b <- fixef(x)
  #standard error - diagnonal and squares them
  se <- sqrt(diag(vc))
  #z score is how many sd the beta is from zero
  z <- b / sqrt(diag(vc))
  #p value calculations 
  #manual calculation using normal distribution
  P <- 2 * (1 - pnorm(abs(z)))
  names(se) <- names(P)
  #searches model results for specific variable name to make sure it's the right row
  #packs beta, se, p value into a single vector for that variable 
  for(i in 1:length(var)) {
   l[[i]] <- c(b[grep(var[[i]], names(b))], se[grep(var[[i]], names(se))], P[grep(var[[i]], names(P))])
  }

 }

 names(l) <- unlist(var)
 return(l)

}


# Return regression coefficients from a lm object

#for standard linear models (no random effects/twins)
return.coeff.lm <- function(x, var) {
 l = list() 
 #check if actually standard linear model 
 if (class(x) != "lm") { 
 
  for(i in 1:length(var)) {
   l[[i]] <- c(NA, NA, NA)
  }

 } else { 
  for(i in 1:length(var)) {
   l[[i]] <- summary(x)$coefficients[grep(var[[i]], rownames(summary(x)$coefficients)), c(1,2,4)]
  }
 }
 
 names(l) <- unlist(var)
 return(l)
}


