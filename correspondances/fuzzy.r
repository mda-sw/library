fuzzy <- function(x)
 {
  
   # Fuzzy or piecewise linear coding of a vector of values. 

   if (!is.vector(x)) stop("Only vectors handled currently as input.\n")
   plin.x <- cbind(x, x)

   lo    <- quantile(x, .33)
   hi    <- quantile(x, .67)

   plin.x[x >= hi, 1] <- 1   
   plin.x[x >= hi, 2] <- 0   

   plin.x[x <= lo, 1] <- 0
   plin.x[x <= lo, 2] <- 1

   plin.x[x > lo & x < hi, 1] <- (x[x > lo & x < hi] - lo) / (hi - lo)
   plin.x[x > lo & x < hi, 2] <- (hi - x[x > lo & x < hi]) / (hi - lo)

   plin.x

 }

