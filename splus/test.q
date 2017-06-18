# Test or example ...
# Do not source: instead, cut and paste



# Form 150x4 matrix of iris data (available as a frame in S-Plus):
iris.dat <- rbind(iris[,,1], iris[,,2], iris[,,3])

iris.pca.corr <- pca(iris.dat)
iris.pca.cov  <- pca(iris.dat, 2)
iris.hc.mvar  <- hierclust(iris.dat)
iris.hc.slnk  <- hierclust(iris.dat, 2)

summary(iris.hc.mvar)
summary(iris.hc.slnk)

motif()
par(mfrow=c(2,2))
lbls <- c(rep("*",50), rep("+",50), rep("o",50))

# Window 1
plot(iris.pca.corr$rproj[,1], iris.pca.corr$rproj[,2], xlab="PC1", ylab="PC2",
    type="n", main="PCA of correlations")
text(iris.pca.corr$rproj[,1], iris.pca.corr$rproj[,2], lbls)
plaxes(iris.pca.corr$rproj[,1], iris.pca.corr$rproj[,2])

# Window 2
plot(iris.pca.cov$rproj[,1], iris.pca.cov$rproj[,2], xlab="PC1", ylab="PC2",
    type="n", main="PCA of covariancess")
text(iris.pca.cov$rproj[,1], iris.pca.cov$rproj[,2], lbls)
plaxes(iris.pca.cov$rproj[,1], iris.pca.cov$rproj[,2])

# Window 3
plot(iris.hc.mvar, labels=lbls)
title("HC/min. var. criterion")

# Window 4
plot(iris.hc.slnk, labels=lbls)
title("HC/slink criterion")

# New window for plot, new window for boxplots
plot(partition(iris.hc.mvar, 3))

# New window for plot, new window for boxplots
plot(partition(iris.hc.slnk, 3))

# New window for plot, new window for boxplots
set.seed(17)
plot(partition(iris.dat, 3))

# Window 1 - Sammon non-metric MDS mapping
set.seed(19)
mds <- sammon(iris.dat)
# Slloooowwww...
# No. iterations = 22, mapping error = 0.0401686
plot(mds$rproj[,1], mds$rproj[,2], xlab="Coord. 1", ylab="Coord. 2", type="n")
text(mds$rproj[,1], mds$rproj[,2], lbls)
plaxes(mds$rproj[,1], mds$rproj[,2])
title("Sammon MDS")

# CA
iris.fuzz <- cbind(flou(iris.dat[,1]), flou(iris.dat[,2]), flou(iris.dat[,3]),
  flou(iris.dat[,4]))
iris.ca  <- ca(iris.fuzz)

# Window 2
plot(iris.ca$rproj[,1], iris.ca$rproj[,2], xlab="Factor 1", ylab="Factor 2",
  type="n")
text(iris.ca$rproj[,1], iris.ca$rproj[,2], lbls)
title("CA of rows")
plaxes(iris.ca$rproj[,1], iris.ca$rproj[,2])

# Window 3
lbls2 <- c("1H", "1L", "2H", "2L", "3H", "3L", "4H", "4L")
plot(iris.ca$cproj[,1], iris.ca$cproj[,2], xlab="Factor 1", ylab="Factor 2",
  type="n", main="CA of variables")
text(iris.ca$cproj[,1], iris.ca$cproj[,2], lbls2)
plaxes(iris.ca$cproj[,1], iris.ca$cproj[,2])


# Full window
par(mfrow=c(1,1))
# Sllllooooooooowwwwwwwwww... for 150 observations
hierclust(iris.dat, movie=T)






