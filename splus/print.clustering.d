.BG
.FN print.clustering
.TL
Print a Clustering Object
.CS
print(x, ...)
.RA
.AG x
object of class `clustering'.  Such an object is produced by function
`partition'.
.SH VALUE
the partition, as input, is returned invisibly.
.SH DETAILS
Pretty-print some of the attributes of a clustering object.  Included are:
number of observations on which the partition is based, number of clusters,
original data set name, function call which produced the clustering object,
etc.
.SA
print, summary, partition, plot
.EX
# Produce a partition with 4 clusters.  Input data set is `a'.
clresult1 <- partition(a, ng=4)
# Next, a hierarchical clustering.
clhier    <- hierclust(a)
# Derive the 4-cluster partition from this.
clresult2 <- partition(clhier, ng=4)
# See how the results compare.  One can use `plot' or `summary' or ...
print(clresult1)
print(clresult2)
.KW multivariate
.KW cluster
.WR

