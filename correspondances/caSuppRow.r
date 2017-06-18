caSuppRow <- function(xtab, rsupp) {

# Projections of supplementary rows.  One or more supplementary rows
# can be used.  Example of use:
# x <- read.table("c:/mandible77s.dat")                 # Read data
# xca <- ca(x[1:77,])                                   # Corr. analysis
# xcar <- caSuppRow(x[1:77,], x[78:86,])                # Suppl. rows
# plot(c(xca$rproj[,1],xca$cproj[,1]),                  # Prepare plot
#       c(xca$rproj[,2],xca$cproj[,2]),
#       type="n", xlab="Factor 1 (50.0% of inertia)",
#       ylab="Factor 2 (21.0% of inertia)")
# text(xca$rproj[,1],xca$rproj[,2],dimnames(x)[[1]])    # Plot prin. rows
# text(xca$cproj[,1],xca$cproj[,2],dimnames(x)[[2]],font=4) # Plot cols.
# text(xcar[,1],xcar[,2],dimnames(xcar)[[1]],font=3)    # Plot supp. rows
# title("77 mandibles, 9 supplementary (groups), crossed by 9 attributes")
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
# rproj <- as.matrix(fJsupI) %*% evectors
temp  <- as.matrix(s2) %*% sres$vectors
# Following divides rowwise by sqrt(fJ) and 
# columnwise by sqrt(eigenvalues):
# Note: first column of cproj is trivially 1-valued.
cproj <- sweep ( sweep(temp,1,sqrt(fJ),FUN="/"), 2,
                 sqrt(sres$values),FUN="/")

# Note: we must coerce rsupp to matrix type, which 
# propagates to rsuppIJ
rsuppIJ <- as.matrix(rsupp)/tot
if (nrow(rsuppIJ) > 1) rsuppI <- apply(rsuppIJ, 1, sum)
if (nrow(rsuppIJ) == 1) rsuppI <- sum(rsuppIJ)

rsuppproj <- rsuppIJ %*% cproj 
temp <- rsuppproj
# Divide cols. by mass; and then rows. by sqrt of evals.
rsuppproj <- sweep ( sweep(temp,1,rsuppI,FUN="/"),2,
                     sqrt(sres$values),FUN="/")

# Value of this function on return: table of projections,
# rows = set of supplementary rows; columns = set of factors.
# (More than 1 supplementary rows => labels will be retained.)
# Adjust for trivial factor.
rsuppproj[,-1]

}

