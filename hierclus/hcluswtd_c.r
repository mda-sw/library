HierClust <- function(arr, wts)
{

# FM, Jan. 2004

# Example of use.
# dyn.load("C:/TEMP/hierclust.dll") 
# arr <- matrix(c(1,4,2,3,0,2,1,2,4,2,3,2),nrow=4,ncol=3)
# wts <- as.single(c(0.7, 0.3, 0.5, 0.8))

# NOTE: By careful with input data.  read.table produces an object of 
# type list.  Convert this to a matrix (e.g. y <- as.matrix(x) ).

order <- as.integer(rep(0,nrow(arr)))
ia <- as.integer(rep(0,nrow(arr)-1))
ib <- as.integer(rep(0,nrow(arr)-1))
height <- as.single(rep(0,nrow(arr)-1))

output <-  .C("HierClust", as.integer(order),
                          as.integer(ia),
                          as.integer(ib),
                          as.single(height),
                          as.single(arr),
                          as.single(wts),
                          as.integer(nrow(arr)),
                          as.integer(ncol(arr)))


retlist <- list(merge=cbind(output[[2]], output[[3]]),
                height=output[[4]],order=output[[1]])

retlist 

}                        
