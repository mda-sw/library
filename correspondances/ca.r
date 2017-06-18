ca <- function(xtab) {

# Correspondence analysis of principal table.  List returned with 
# projections, correlations, and contributions of rows (observations),
# and columns (attributes).  Eigenvalues are output to display device.
# FM, 2003/12.

tot <- sum(xtab)
fIJ <- xtab/tot
fI <- apply(fIJ, 1, sum)
fJ <- apply(fIJ, 2, sum)
if (length(fI[fI <= 0]) >= 1) cat("Note and check: fI terms le 0. <-- 1.\n")
if (length(fJ[fJ <= 0]) >= 1) cat("Note and check: fJ terms le 0. <-- 1.\n")
# Note on following: any positive value is reqd. for mass; distance will be
# 0.  Hence inertia contribution will be 0.  FM, 2008/7.
fI[fI <= 0] <- 1
fJ[fJ <= 0] <- 1
fJsupI <- sweep(fIJ, 1, fI, FUN="/")
fIsupJ <- sweep(fIJ, 2, fJ, FUN="/")
s <- as.matrix(t(fJsupI)) %*% as.matrix(fIJ)
s1 <- sweep(s, 1, sqrt(fJ), FUN="/")
s2 <- sweep(s1, 2, sqrt(fJ), FUN="/")
# In following s2 is symmetric.  However due to precision S-Plus didn't 
# find it to be symmetric.  And function eigen in S-Plus uses a different
# normalization for the non-symmetric case (in the case of some data)!  
sres <- eigen(s2,symmetric=T)
sres$values[sres$values < 1.0e-8] <- 0.0
cat("Eigenvalues follow (trivial first eigenvalue removed).\n")
cat(sres$values[-1], "\n")
cat("Eigenvalue rate, in thousandths.\n")
tot <- sum(sres$values[-1])
cat(1000*sres$values[-1]/tot,"\n")
# Eigenvectors divided rowwise by sqrt(fJ):
evectors <- sweep(sres$vectors, 1, sqrt(fJ), FUN="/")

# PROJECTIONS ON FACTORS OF ROWS AND COLUMNS
rproj <- as.matrix(fJsupI) %*% evectors
temp  <- as.matrix(s2) %*% sres$vectors
# Following divides rowwise by sqrt(fJ) and columnwise by sqrt(eigenvalues):
# Note: first column of cproj is trivially 1-valued.
# NOTE: VBESxFACTORS. READ PROJS WITH FACTORS 1,2,... FROM COLS 2,3,...
cproj <- sweep(sweep(temp,1,sqrt(fJ),FUN="/"),2,sqrt(sres$values),FUN="/")

# CONTRIBUTIONS TO FACTORS BY ROWS AND COLUMNS
# Contributions: mass times projection distance squared.
temp <- sweep( rproj^2, 1, fI, FUN="*")
# Normalize such that sum of contributions for a factor equals 1.
sumCtrF <- apply(temp, 2, sum)
# Note: Obs. x factors. Read cntrs. with factors 1,2,... from cols. 2,3,...
rcntr <- sweep(temp, 2, sumCtrF, FUN="/")
temp <- sweep( cproj^2, 1, fJ, FUN="*")
sumCtrF <- apply(temp, 2, sum)
# Note: Vbs. x factors. Read cntrs. with factors 1,2,... from cols. 2,3,...
ccntr <- sweep(temp, 2, sumCtrF, FUN="/")

# CORRELATIONS WITH FACTORS BY ROWS AND COLUMNS
# dstsq(i) = sum_j 1/fj (fj^i - fj)^2
temp <- sweep(fJsupI, 2, fJ, "-")
dstsq <- apply( sweep( temp^2, 2, fJ, "/"), 1, sum)
# NOTE: Obs. x factors. Read corrs. with factors 1,2,... from cols. 2,3,...
rcorr <- sweep(rproj^2, 1, dstsq, FUN="/")
temp <- sweep(fIsupJ, 1, fI, "-")
dstsq <- apply( sweep( temp^2, 1, fI, "/"), 2, sum)
# NOTE: Vbs. x factors. Read corrs. with factors 1,2,... from cols. 2,3,...
ccorr <- sweep(cproj^2, 1, dstsq, "/")

# Value of this function on return: list containing projections, correlations,
# and contributions for rows (observations), and for columns (variables).
# In all cases, allow for first trivial first eigenvector.
list(rproj=rproj[,-1], rcorr=rcorr[,-1], rcntr=rcntr[,-1],
     cproj=cproj[,-1], ccorr=ccorr[,-1], ccntr=ccntr[,-1])

}

