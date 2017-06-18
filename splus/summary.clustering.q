summary.clustering <- function(group, ...)
{

      if (!inherits(group,"clustering")) stop("Not legitimate clustering\n")
      ng <- length(unique(group))

      n <- length(group)
      cat("Partition with number of clusters = ",ng,"\n")

      data     <- get(attr(group, "origdata"))
      if (is.matrix(data))    m   <- ncol(data)
      if (!is.matrix(data))   m   <- 1


      if (m > 1) {
#     Produce and output statistics on clusters
      for (i in 1:ng) {
               gdata <- data[group==i,]
               if (is.matrix(gdata)) {
                  xmean <- apply(gdata, 2, mean)
                  xvar  <- apply(gdata, 2, var)
                  xmax  <- apply(gdata, 2, max)
                  xmin  <- apply(gdata, 2, min)
               }
               if (!is.matrix(gdata)) {
                  xmean <- gdata
                  xvar  <- rep(0.0, length(gdata))
                  xmax  <- gdata
                  xmin  <- gdata
               }
               cat(" \n")
               cat("Cluster: ",i,"\n")
               alldat     <- rbind(round(xmean,4),round(xvar,4),
                                   round(xmax,4), round(xmin,4) )
              dimnames(alldat)[[1]]<- c("Means","Variances","Maxima","Minima")
#            dimnames(alldat)[[1]]<- list("Means","Variances","Maxima","Minima")
               prmatrix(alldat)
          }
      }

    if (m == 1) {
       cat("Either 1-dim. input data, or diss. input.\n")
       cat("Nothing available to summarize.  Sorry.\n")
    }

    invisible(group)
}

           
