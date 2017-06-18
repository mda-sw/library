plot.clustering <- function(group, rep="pc", ...)
{


# Function to implement 'plot' method for objects of class 'clustering'
# Author: F. Murtagh, May 1992



      if (!inherits(group,"clustering")) stop("Not legitimate clustering\n")
      ng <- length(unique(group))



#     Note that following won't effect boxplots for unidim. data
      if (rep != "pc" && rep != "pa") 
      stop("Unacceptable representation requested.\n")



      n <- length(group)
      cat("Partition with number of clusters = ",ng,"\n")

         data     <- eval(attr(group, "origdata"))


      if (!inherits(data,"reddim")) {
         datanam  <- attr(group, "origdata")
         out      <- 
         paste("We will use raw or diss. data, or eval. expr.:",datanam,"\n")
         cat(out)
      }
      if (inherits(data,"reddim")) {
         data     <- data$rproj
         datanam  <- attr(group, "origdata")
         out      <-
         paste("We will use rproj component of reddim object:",datanam,"\n")
         cat(out)
      }


# Remark: in following, 'data' may in fact be half-matrix of dissimilarities
      if (is.matrix(data))    m   <- ncol(data)
      if (!is.matrix(data))   m   <- 1
#     Now check out situation for dissimilarity input
      dst    <- FALSE
      if (!is.matrix(data)) {
         lngth  <- n*(n-1)/2
         if (lngth == length(data))   dst <- TRUE
      }



#   Look after new window.
#   cat(" ","\n")
    cat("------------------------------------------------------------------",
        "\n")
    cat("Opening new window to show partition. Position using mouse button.\n")
    motif()




#   'text' screams when 'group' has attributes.  So temporarily remove them.
    cl    <- attr(group, "class")
    cal   <- attr(group, "call")
    ori   <- attr(group, "origdata")
    attr(group, "class")    <- NULL
    attr(group, "call")     <- NULL
    attr(group, "origdata") <- NULL




#   Case 1: unidimensional data
    if (!dst && m == 1) boxplot(split(data,group), xlab="Classes", ylab=
              "Values of unidimensional variable", main=
              "Boxplot of variable's values by class")  




#   Case 2: multidimensional data, first two principal components
    if (rep == "pc" && m > 2) {
            cat("First two principal components used to display partition.\n")
            pc <- pca(data)
            titl <- paste("Data set ",ori,
                      ".  Partition into ", ng," clusters.")
            plot(pc$rproj[,1],pc$rproj[,2],xlab="Principal component 1",
                      ylab="Principal component 2",main=titl,type="n")
            text(pc$rproj[,1],pc$rproj[,2],group)
            plaxes(pc$rproj[,1],pc$rproj[,2])
            }


    if (rep == "pc" && m == 2) {
            cat("Two coordinates used to display partition.\n")
            titl <- paste("Data set ",ori,
                      ".  Partition into ", ng," clusters.")
            plot(data[,1], data[,2],xlab="Coordinate 1",
                      ylab="Coordinate 2",main=titl,type="n")
            text(data[,1],data[,2],group)
            plaxes(data[,1],data[,2])
            }

    if (rep == "pc" && dst) {
            cat("First two principal coordinates used to display partition.\n")

#           Temporarily remove all attributes, to use old-S fn. 'cmdscale'
#           But note: 'Size' attr. must be present!!
            attda     <- attr(data, "class")
            attor     <- attr(data, "origdata")
            attme     <- attr(data, "metric")
            attr(data, "class")    <- NULL
            attr(data, "origdata") <- NULL
            attr(data, "metric")   <- NULL

            pc <- cmdscale(data)

            attr(data, "class")    <- attda
            attr(data, "origdata") <- attor
            attr(data, "metric")   <- attme
            
            titl <- paste("Data set ", attor,
                     ".  Partition into ", ng," clusters.")
            plot(pc[,1], pc[,2], xlab="Principal coordinate 1",
                     ylab="Principal coordinate 2",main=titl,type="n")
            text(pc[,1], pc[,2], group)
            plaxes(pc$rproj[,1],pc$rproj[,2])
            }




#   Case 3: multidimensional data, all-pairs plot
    if (!dst && rep == "pa" && m > 2) {
    pairs(data, glab=group, panel=function(x, y,glab) { points(x,y,type="n");
                                 text(x,y,glab) } )
    }




    cat("Delete plot window? [y/Y]\n")
    resp <- scan("",what=character(),1)
    if (length(resp)==0 || resp == "y" || resp == "Y") 
        dev.off()
    else
    cat("Note - plot window is current ('dev.off()' to delete).\n")




    if (ng <= 4 && m <= 4 && !dst) {
    cat("No. classes and no. vbes. both <= 4: boxplots follow.\n")
    cat("Opening new window to show partition. Position using mouse button.\n")
    motif()
    par(mfrow=c(m,1))
    for (vbe in 1:m) {
        gdata  <- data[,vbe]
        vbelab <- paste("Variable ",vbe) 
        boxplot(split(gdata,group), ylab=vbelab)
    }
    title(sub="Clusters")
    cat("Delete plot window? [y/Y]\n")
    resp <- scan("",what=character(),1)
    if (length(resp)==0 || resp == "y" || resp == "Y") 
        dev.off()
    else
        cat("Note - plot window is current ('dev.off()' to delete).\n")
    }


 

    attr(group, "class")    <- cl
    attr(group, "call")     <- cal
    attr(group, "origdata") <- ori



    invisible(group)


}


           


