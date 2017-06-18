facor <- function(kIJ, nclr=nrow(kIJ), nclc=ncol(kIJ), nf=2) {

# Example of reading input data.
# parm <- read.table("c:/parmenides.dat")
# kIJ <- parm[, -3:-4]

# Definition of total mass, frequencies, profile vectors on I and on J:
k       <- sum(kIJ)
fIJ     <- kIJ/k
fJI     <- t(fIJ)
fI      <- apply(fIJ, 1, sum)
fJ      <- apply(fIJ, 2, sum)
fJsupI  <- sweep(fIJ, 1, fI, FUN="/")
fIsupJ  <- sweep(fIJ, 2, fJ, FUN="/")

# Carry out Correspondence Analysis. Diagonalization: definition of factors.
# Following three lines yield matrix to be diagonalized (SVD):
# s_jj' = sum_i (fij * fij') / (fi * sqrt(fj) * sqrt(fj'))
s <- as.matrix(t(fJsupI)) %*% as.matrix(fIJ)
s1 <- sweep(s, 1, sqrt(fJ), FUN="/")
s2 <- sweep(s1, 2, sqrt(fJ), FUN="/")
sres <- eigen(s2)
# Eigenvectors divided rowwise by sqrt(fJ):
evectors <- sweep(sres$vectors, 1, sqrt(fJ), FUN="/")
# Projections on factors of rows and columns
# Note: first column of rproj is trivially 1-valued.
# Note: Observations x factors.  Read projections with factors 1, 2, ... from
# cols. 2, 3, ...
rproj <- as.matrix(fJsupI) %*% evectors
temp  <- as.matrix(s2) %*% sres$vectors
# Following divides rowwise by sqrt(fJ) and columnwise by sqrt(eigenvalues):
# Note: first column of cproj is trivially 1-valued.
# Note: Variables x factors.  Read projectsion with factors 1, 2, ... from
# cols; 2, 3, ...
cproj <- sweep ( sweep(temp,1,sqrt(fJ),FUN="/"), 2,sqrt(sres$values),FUN="/")

# Create hierarchical clusterings on observations (rows) and variables 
# (columns), based on factor projections.  Use weights of obs. or vbes.  
# We could limit the number of factors used here, e.g. rproj[,2:nlimit]
hclr <- hierclust( rproj[,-1], fI)
hclc <- hierclust( cproj[,-1], fJ)

labsr     <- 1:nclr
labsallr  <- dimnames(fIJ)[[1]]
membersr  <- cutree(hclr, nclr)

labsc     <- 1:nclc
labsallc  <- dimnames(fIJ)[[2]]
membersc  <- cutree(hclc, nclc)

centersr  <- NULL
for (k in 1:nclr) {
  if (length(fI[membersr==k]) > 1) 
    centersr <- rbind(centersr, apply(rproj[membersr==k,-1],2,sum))
  if (length(fI[membersr==k]) == 1)    # Card of cluster = 1 => bypass summing.
    centersr <- rbind(centersr, rproj[membersr==k,-1])
  labsr[k] <-  list(labsallr[membersr==k])
}

centersc  <- NULL
for (k in 1:nclc) {
  if (length(fJ[membersc==k]) > 1) 
    centersc <- rbind(centersc, apply(cproj[membersc==k,-1],2,sum))
  if (length(fJ[membersc==k]) == 1)    # Card of cluster = 1 => bypass summing.
    centersc <- rbind(centersc, cproj[membersc==k,-1])
  labsc[k] <-  list(labsallc[membersc==k])
}

list(centersr=centersr, cluslabrow=labsr, centersc=centersc, cluslabcol=labsc) 

}
