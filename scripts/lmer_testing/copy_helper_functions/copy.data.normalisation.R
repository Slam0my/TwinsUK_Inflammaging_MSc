inverse.normal <- function(data)
# Data inverse-normalisation
# @author Niccolò Rossi
# 
# Args:
#	data		 : vector (can be non-numeric)

# Output:
# 	inverse-normalised numeric vector

{
	v <- as.numeric(as.character(data))
	v <- round(qnorm((rank(v, na.last="keep")-0.5)/sum(!is.na(v))), digits=2)
	v
}


remove.outliers <- function(data, t)
# remove outliers from
# a numeric vector
#
# @author Niccolò Rossi
# 
# Args:
#	data		 : vector (can be non-numeric)
#	t		 : sd threshold 
# Output:
# 	numeric vector, outliers are set to NA

{
	v <- as.numeric(as.character(data))
	
	m <- mean(v, na.rm = T)
	s <- sd(v, na.rm = T) * t
	v <- sapply(v, function(x) {ifelse(x > m+s | x < m-s, NA, x)})
	
	v
}

