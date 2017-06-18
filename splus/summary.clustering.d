.BG
.FN summary.clustering
.TL
Summary of a Clustering Object
.CS
summary(x, ...)
.RA
.AG x
object of class `clustering'.  Such an object is produced by function
`partition'.
.SH VALUE
the partition, as input, is returned invisibly.
.SH DETAILS
Means, variances, maxima and minima are produced for each cluster's members.
If the original input was dissimilarity data, or if the original coordinate
data was one-dimensional, then nothing can be summarized.
.SA
summary, partition, plot, print
.EX
# Produce a partition with 4 clusters.  Input data set is `a'.
clresult1 <- partition(a, ng=4)
# Next, a hierarchical clustering.
clhier    <- hierclust(a)
# Derive the 4-cluster partition from this.
clresult2 <- partition(clhier, gp=4)
# See how the results compare.  One can use `plot' or `print' or ...
summary(clresult1)
summary(clresult2)
.KW multivariate
.KW cluster
.WR

