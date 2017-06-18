facor <- function(a, nhier=4, nf=2) {

# R PROGRAM facor PURPOSE
# Carry out a Correspondence Analysis of input matrix, a.
# Consider the best factor space of dimensionality, nf.
# Carry out a hierarchical clustering on row projections, and on 
#   column projections in the factor space of dimensionaly nf. 
#   In both cases, use the Correspondence Analysis row and column 
#   masses as weights.  
# Choose a partition for both hierarchical clusterings, with nhier clusters.
# Determine the coordinates of the clusters (on either rows, or columns) 
#   in the factor space of dimensionality nf.
# Author: F. Murtagh, 2010 Jan. 17.  f m u r t a g h _at_ a c m _dot_ o r g 

# INPUT DATASET
# File casa2.prn characterizes the 77 successive scenes of the film 
# Casablanca by 13 characterics related to the characters and location,
# defined in the following labels:
# collabels <- c("No","Int","Ext","Day","Night","Action","Dialog","Rick",
#  "Ilsa","Renault","Strasser","Laszlo","Minor")
# "No" is sequence number of scene. "Int", "Ext": interior, exterior.
# "Minor" is minor character.  

# EXAMPLE OF USE
# cas <- matrix(scan("~/casablanca/casa2.prn"),nrow=77,ncol=13,byrow=T)
# casc <- ca(cas[,-1])   # Col 1 is just seq nos.
# casfacor <- facor(cas[,-1], nhier=8)  # Same CA is repeated -note same evals.
# Set up plot:
# plot(casc$rproj[,1], casc$rproj[,2], type="n")
# Plot the row cluster centers:
# nhier <- 8
# text(casfacor$centersr[,1], casfacor$centersr[,2],1:nhier)
# Plot rows/observations:
# text(casc$rproj[,1], casc$rproj[,2], 1:77, col="red")
# Plot the col.	       cluster	       centers:
# text(casfacor$centersc[,1], casfacor$centersc[,2],1:10,col="blue")
# Plot cols./attributes:
# text(casc$cproj[,1], casc$cproj[,2], 1:13, col="green")

# PARAMETERS
#  nf is number of factor to use for the clustering.  
#  If there are linear dependencies, then one cannot use all 
#  factors which are, resulting from the dependency, undefined.
#  0 eigenvalues express this.  Recommend: use nf = 2 or other low number.

# ASSUMPTIONS
# Not checked: a, input data matrix, must contain non-negative values.
# Marginal sums of a must not contain a 0 value.  Else zero-divides occur.

# REQUIRED 
#  Require:   ca.r         Correspondence Analysis.
#  Require:   hcluswtd.r   Hierarchical Clustering, using weights.
#  Require:   facor-new.r  facor, Factor coordinates of clusters.

cat("nhier =", nhier, "\n")
   k <- sum(a)
   fIJ <- a/k
   fI <- apply(fIJ, 1, sum)
   fJ <- apply(fIJ, 2, sum)

   cac <- ca(fIJ)

   nclr <- nrow(fIJ)
   nclc <- ncol(fIJ)
   rproj <- cac$rproj[,1:nf]
   cproj <- cac$cproj[,1:nf]

   hclr <- hierclust(rproj, fI)
   hclc <- hierclust(cproj, fJ)

labsr     <- 1:nclr
labsallr  <- dimnames(fIJ)[[1]]
# membersr  <- cutree(hclr, nclr) # I.e. fine partition; reproduces data exactly.
membersr  <- cutree(hclr, nhier)  # NEW

labsc     <- 1:nclc
labsallc  <- dimnames(fIJ)[[2]]
# membersc  <- cutree(hclc, nclc) # I.e. fine partition; reproduces data exactly.
membersc  <- cutree(hclc, nhier)    # NEW

centersr  <- NULL
#for (k in 1:nclr) {
for (k in 1:nhier) {      # NEW
  if (length(fI[membersr==k]) > 1) 
    centersr <- rbind(centersr, apply(rproj[membersr==k,],2,sum))
  if (length(fI[membersr==k]) == 1)    # Card of cluster = 1 => bypass summing.
    centersr <- rbind(centersr, rproj[membersr==k,])
  labsr[k] <-  list(labsallr[membersr==k])
}
#  NEW  Averaging (in equiwtd. Eucl. factor space)  the cluster coords.
for (k in 1:nhier) {
   centersr[k,]  <- centersr[k,]/length(fI[membersr==k])   
}

centersc  <- NULL
#for (k in 1:nclc) {
for (k in 1:nhier) {    #  NEW
  if (length(fJ[membersc==k]) > 1) 
    centersc <- rbind(centersc, apply(cproj[membersc==k,],2,sum))
  if (length(fJ[membersc==k]) == 1)    # Card of cluster = 1 => bypass summing.
    centersc <- rbind(centersc, cproj[membersc==k,])
  labsc[k] <-  list(labsallc[membersc==k])
}

#  NEW  Averaging (in equiwtd. Eucl. factor space)  the cluster coords.
for (k in 1:nhier) {
   centersc[k,]  <- centersc[k,]/length(fJ[membersr==k])   
}

list(centersr=centersr, cluslabrow=labsr, centersc=centersc, cluslabcol=labsc) 

}

