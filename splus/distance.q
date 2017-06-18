distance <- function(x, metric="euclidean")
{

# Function to use, and enhance, current S+ function 'dist' 
# Also to attach appropriate attributes, for subsequent calls to hier. 
# clust. functions, etc.
# Author: F. Murtagh, May 1992

    
   dst <- dist(x, metric)

   
#  'dist' provides attr. 'Size'.   Now, we'll add a few more.
   attr(dst, "class")       <- "dissdat"
#  Note: potential pitfall in following: input data is assumed located in
#  particular position in call arguments.  Fix up later.
   attr(dst, "origdata")   <- parse(text=match.call())[2]
   attr(dst, "metric")     <- metric

   dst

}


