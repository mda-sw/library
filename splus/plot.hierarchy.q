plot.hierarchy <- function(x, ...)
{
# Function to implement 'plot' method for objects of class 'hierarchy'
# Author: F. Murtagh, May 1992


         if (!inherits(x, "hierarchy"))
                  stop("Not legitimate hierarchy")


#        'plclust' doesn't like its main argument to have some 'attr's, so we
#        temporarily remove them

         cl                       <- attr(x, "class")
         cal                      <- attr(x, "call")
         wh                       <- attr(x, "criterion")
         ori                      <- attr(x, "origdata")
         ord                      <- attr(x, "origdist")
         attr(x, "class")         <- NULL
         attr(x, "call")          <- NULL
         attr(x, "criterion")     <- NULL
         attr(x, "origdata")      <- NULL
         attr(x, "origdist")      <- NULL
       
         plclust(x, ...)

         attr(x, "class")         <- cl
         attr(x, "call")          <- cal
         attr(x, "criterion")     <- wh
         attr(x, "origdata")      <- ori
         attr(x, "origdist")      <- ord

         invisible(x)
}


