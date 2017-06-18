partition <- function(a, g, iter.max=15, option="mindst", diagnostics=FALSE,
                      showslice=FALSE)
{
#
#          Carry out partitioning; produce object of class 'clustering'
#          Author: F. Murtagh, May 1992

#          `g' is either: (i) number of clusters, in which case `option'
#          will be used (default option: "mindst"; alternative: "exch"; for
#          small datasets, former is better behaved - latter gives empty
#          clusters, and Inf of NaN center and compactness values);
#          or (ii) `g' is set of initial guesses for the cluster centers.

#          `a' is an object of class `hierarchy' or of class `coordata' or
#          of class `dissdata'.


# HERE IS THE BIG DIVIDE... between current S+ routine 'kmeans', and 
# alternative routine 'exch' and 'mindst'.  But first, in accordance with
# class of `a'.

if (inherits(a, "hierarchy")) {
           n <- length(a$order)
           if (g != 0) ht <- n - g + 1
           if (g < 1 || g > n) stop("Invalid number of clusters requested\n")
           htreal <- a$height[ht-1]
           if (diagnostics == TRUE) {
              cat("Partition with number of clusters = ",g,"\n")
              cat("Corresponds to dendrogram height  = ",
                  round(as.single(htreal),2),"\n")
            }
  
# Draw assoc. dotted line thru dendrogram, assuming window open -- and
# assuming that dendrogram is already plotted.

          if (length(dev.list()) > 0 && showslice) abline(htreal,0.0,lty=2)


# Following inelegant removing, and replacing, of attributes from 'xdata' is
# due to use of 'cutree'.  Latter screams if there are attributes present.
               cl                           <- attr(a, "class")
               cal                          <- attr(a, "call")
               wh                           <- attr(a, "criterion")
               ori                          <- attr(a, "origdata")
               ord                          <- attr(a, "origdist")
               attr(a, "class")         <- NULL
               attr(a, "call")          <- NULL
               attr(a, "criterion")     <- NULL
               attr(a, "origdata")      <- NULL
               attr(a, "origdist")      <- NULL

               partition   <- cutree(a, g, ht)

               attr(a, "class")         <- cl
               attr(a, "call")          <- cal
               attr(a, "criterion")     <- wh
               attr(a, "origdata")      <- ori
               attr(a, "origdist")      <- ord

#              Ok, we now have what we came for...
#              attr(partition, "height")   <- NULL
               attr(partition, "class")    <- "clustering"
               attr(partition, "call")     <- deparse(sys.call())
               attr(partition, "origdata") <- ori
               attr(partition, "origdist") <- ord

}

if (inherits(a, "reddim")) a <- a$rproj

if (!inherits(a, "hierarchy")) {
# Assume that `a' is of class `coordata' or of class `dissdata'.



           if (is.matrix(g))
{
# HARTIGAN & WONG ALGORITHM
           partition <- kmeans(a, g, iter.max)

           attr(partition, "class")      <-  "clustering"
           attr(partition, "call")       <-  deparse(sys.call())
#Note:  we assume argument 'a' (i/p data) is at partic. posn. in arg. list.
           attr(partition, "origdata")   <-  parse(text=match.call())[2]
           attr(partition, "origdist")   <-  "Euclidean distance"
}
           else
{
# SPAETH ALGORITHMS (MURTAGH IMPLEMENTATION)
           ng <- g
           n <- nrow(a)
           m <- ncol(a)
           gpcen <- matrix(0.0, ng, m)
           storage.mode(gpcen) <- "single"
#          Check for valid number of groups
           if (ng <= 1) stop("ng (number of groups) must be > 1.")
           if (ng >= n) 
                    stop("ng (number of groups) must be < number of objects.")
#          Assign objects arbitrarily to groups.  Generate n uniform values,
#          and scale these up by just less than ng (= number of groups).  
#          'Floor'ing this gives an assignment in the range 0..ng-1 to each
#          object, so add 1 to this.
           memgp <- floor(runif(n)*(ng-0.00001))+1
#          Lowest acceptable number of objects assigned to a group:
           ng0 <- 1
#          Initial value of error return indicator:
           ierr <- 0
#          Now, go for it ...
           if (option == "exch") 
              partit <- .Fortran("exch",
                              a = as.single(a),      # input data
                              n = as.integer(n),
                              m = as.integer(m),
                              memgp = as.integer(memgp),  # group memberships
                              ngp0 = integer(1),     # smallest accept. cardin.
                              numgp = integer(ng),   # num. objects per group
                              gpcen = as.matrix(gpcen),# gp. cent. coords.
                              ng = as.integer(ng),   # number of gps. requested
                              comp = single(ng),     # group compactness vals.
                              ctot = single(1),      # sum of comp values
                              iter = integer(1),     # iterns. to converge
                              iter.max = as.integer(iter.max),
                              ierr = integer(1))     # error indicator
        else   
           partit <- .Fortran("mindst",
                              a = as.single(a),      # input data
                              n = as.integer(n),
                              m = as.integer(m),
                              memgp = as.integer(memgp),  # group memberships
                              ngp0 = integer(1),     # smallest accept. cardin.
                              numgp = integer(ng),   # num. objects per group
                              gpcen = as.matrix(gpcen),# gp. cen. coords.
                              ng = as.integer(ng),   # number of gps. requested
                              comp = single(ng),     # group compactness vals.
                              ctot = single(1),      # sum of comp values
                              iter = integer(1),     # iterns. to converge
                              iter.max = as.integer(iter.max),
                              ierr = integer(1))     # error indicator
           if (partit$ierr==1) stop("Invalid group number.")
#          Above error due to an object's group assignment being <1 or >ng,
#          in routine gmeans, called from mindst or from exch, which are the 
#          two principal routines stored in file part.f.  Could something
#          have gone wrong with the initial, arbitrary assignment of group
#          memberships at the beginning of this routine?  Hard to imagine...
           if (partit$ierr==2) 
            stop("A gp. has too few mems. - reduce # of gps. and try again.")


           if (partit$numgp[1] == 0) {
            cat("Warning: 1 or more classes are empty.\n")
            cat("Note that $centers values consequently contain NaN values.\n")
           } 

           partition <- partit$memgp

#          Resequence the cluster nos.
#          print(partition)
           tmp0 <- unique(partition)
           tmp  <- partition
           for (ich in 1:length(tmp0)) {
               tmp[partition==tmp0[ich]] <- ich
           }
           partition <- tmp
#          print(partition)
           if (diagnostics) {
           cat("Within sum of squares:", partit$comp,"\n")
           cat("Cluster cardinalities:", partit$numgp,"\n")
           cat("Cluster centers (displayed horizontally for each cluster):", 
              "\n")
              prmatrix(partit$gpcen)    
           }
   }

   attr(partition, "class")      <-  "clustering"
   attr(partition, "call")       <-  deparse(sys.call())
#Note:  we assume argument 'a' (i/p data) is at partic. posn. in arg. list.
   attr(partition, "origdata")   <-  parse(text=match.call())[2]
   attr(partition, "origdist")   <-  "Euclidean distance"

}


   partition

}

