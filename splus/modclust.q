modclust <- function(x, signif = rep(0, dim(x)[2]), method = "S*", 
                 noise = F, scale = rep(1,dim(x)[2]), 
                 shape = c(1, rep(0.2, dim(x)[2] - 1)), 
                 workspace = (dim(x)[1] * (dim(x)[1] - 1))/2 + 10 * dim(x)[1])
{

    hier <- mclust(x, signif, method, noise, scale, shape, workspace)

    hierar <- hier$tree
    hierar$awe <- hier$awe  


    attr(hierar,"class")         <- "hierarchy"
    attr(hierar,"call")          <- deparse(sys.call())
    attr(hierar,"criterion")     <- method
#Note: we assume arg. 'a' (i/p data/dist.) is at particular posn. in par. list.
    attr(hierar,"origdata")      <- parse(text=match.call())[2]
    attr(hierar, "origdist")     <- "weighted Euclidean distance"

 
    hierar

}

