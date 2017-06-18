caSuppCol <- function(xtab, csupp) {

# Projections of supplementary columns.  One or more 
# supplementary columns can be used.  
# Example of use on gobelets data (25 x 6):
# suppl <- caSuppCol(gobelets[,-3], gobelets[,3])
# This uses all except column 3 as principal; and 
# column 3 as supplementary.
# FM, 2003/12.

tot <- sum(xtab)
fIJ <- xtab/tot
fI <- apply(fIJ, 1, sum)
fJ <- apply(fIJ, 2, sum)
if (length(fI[fI <= 0]) >= 1) cat("Note and check: fI terms le 0. <-- 1.\n")
if (length(fJ[fJ <= 0]) >= 1) cat("Note and check: fJ terms le 0. <-- 1.\n")
# Note on following: any positive value is reqd. for mass; distance will be
# 0.  Hence inertia contribution will be 0.   FM, 2008/7.
fI[fI <= 0] <- 1
fJ[fJ <= 0] <- 1
fJsupI <- sweep(fIJ, 1, fI, FUN="/")
fIsupJ <- sweep(fIJ, 2, fJ, FUN="/")
s <- as.matrix(t(fJsupI)) %*% as.matrix(fIJ)
s1 <- sweep(s, 1, sqrt(fJ), FUN="/")
s2 <- sweep(s1, 2, sqrt(fJ), FUN="/")
sres <- eigen(s2)
sres$values[sres$values < 1.0e-8] <- 0.0
# cat("Eigenvalues follow (trivial first eigenvalue removed).\n")
# cat(sres$values[-1], "\n")
# cat("Eigenvalue rate, in thousandths.\n")
# tot <- sum(sres$values[-1])
# cat(1000*sres$values[-1]/tot,"\n")
# Eigenvectors divided rowwise by sqrt(fJ):
evectors <- sweep(sres$vectors, 1, sqrt(fJ), FUN="/")
rproj <- as.matrix(fJsupI) %*% evectors

# Note: we must coerce csupp to matrix type, which 
# propagates to csuppIJ
csuppIJ <- as.matrix(csupp)/tot
if (ncol(csuppIJ) > 1) csuppJ <- apply(csuppIJ, 2, sum)
if (ncol(csuppIJ) == 1) csuppJ <- sum(csuppIJ)
csuppproj <- t(csuppIJ) %*% rproj
temp <- csuppproj
# Divide rows by mass; and then cols. by sqrt of evals.
csuppproj <- sweep ( sweep(temp,1,csuppJ,FUN="/"),2,
     sqrt(sres$values),FUN="/")

# Value of this function on return: table of projections,
# rows = set of supplementary columns; columns = set of factors.
# (More than 1 supplementary column => labels will be retained.)
# Adjust for trivial factor.
csuppproj[,-1]

}