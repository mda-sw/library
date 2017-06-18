print.hierarchy <- function(x, ...)
{
# Function to implement 'print' method for objects of class 'hierarchy'

             if (!inherits(x, "hierarchy"))
                     stop("Not legitimate hierarchy")

             cat("\nSequence of merges:\n")
             print(x$merge)

             cat("\nList of node heights:\n")
             print(x$height)

             cat("\nList of 'horizontal' ordering of observations:\n")
             print(x$order)

             invisible(x)
}

