print.clustering <- function(obj, ...)
{
# Function to implement 'print' method for objects of class 'clustering'

          if (!inherits(obj, "clustering"))
                  stop("Not legitimate clustering")

          n                       <- length(obj)
          ng                      <- length(unique(obj))
          partn                   <- obj
          attr(partn, "class")    <- NULL
          attr(partn, "call")     <- NULL
          attr(partn, "origdata") <- NULL
          cal                     <- attr(obj, "call")

          cat("Number of observations: ", n, "\n")
          cat("Number of clusters:     ", ng, "\n")
          cat("Partition:              ", "\n", partn, "\n")
          cat("Original data:          ", attr(obj, "origdata"), "\n")
          cat("Call which set this up: ", cal, "\n")

          invisible(obj)
}

