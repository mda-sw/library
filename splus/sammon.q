sammon <- function(a, p=2, maxit=100, tol=0.05, alpha=0.3, diagnostics=F)
{

# Sammon (non-metric multidimensional scaling) mapping.  
# Author: F. Murtagh, May 1992

   n <- nrow(a)
   m <- ncol(a)
   storage.mode(a) <- "single"
   b <- matrix(runif(n*p), n, p)
   storage.mode(b) <- "single"
   dstar <- as.vector(dist(a))
   dd    <- as.vector(dist(b))
   ndis  <- length(dstar)
   iter  <- 0
   err   <- 1.e+5
   diag  <- 0
   if (diagnostics) {
#     Note that most diagnostics are from within the Fortran program!
      cat("Mapping error\n")
      diag <- 1
   }

   redmap <- .Fortran("sammon",
             n     = as.integer(n),
             m     = as.integer(m),
             p     = as.integer(p),
             a     = as.matrix(a),
             b     = as.matrix(b),
             ndis  = as.integer(ndis),
             dstar = as.single(dstar),
             dd    = as.single(dd),
             alpha = as.single(alpha),
             maxit = as.integer(maxit),
             diag  = as.integer(diag),
             iter  = as.integer(iter),
             tol   = as.single(tol),
             err   = as.single(err))

   cat("Number of iterations: ",redmap$iter,"\n")
   cat("Mapping error:        ",redmap$err,"\n")

   rproj <- redmap$b
   res   <- list(rproj = rproj)
   attr(res,"class") <- "reddim"

   res

}


