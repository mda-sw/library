summary.hierarchy <- function(object, ...)
{
# Method `summary' for objects of class `hierarchy'

         if (!inherits(object, "hierarchy"))
                stop("Not legitimate hierarchy")

         obj              <- list(call = object$call)
         obj$type         <- "\nHierarchical clustering:\n"
         obj$call         <- attr(object, "call")
         obj$nobs         <- length(object$order)
         obj$criterion    <- attr(object, "criterion")
         obj$bas          <- attr(object, "origdata")
         obj$basd         <- attr(object, "origdist")
         class(obj)  <- "summary.hierarchy"

         cat("\n",obj$type,"\n")
         cat(obj$call,"\n")
         cat("Number of observations: ", obj$nobs, "\n")
         cat("Method used:            ", obj$criterion, "\n")
         cat("Original data:          ", obj$bas, "\n")
         cat("Dissimilarity used:     ", obj$basd, "\n")

         invisible(obj)
}

