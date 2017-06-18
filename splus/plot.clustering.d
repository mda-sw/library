.BG
.FN plot.clustering
.TL
Plot a Clustering Object
.CS
plot(x, rep="pc", ...)
.RA
.AG x
object of class `clustering'.  Such an object is produced by function
`partition'.
.AG rep
preferred representation of the partition: acceptable values are `pc'
(default) 
for principal plane of a principal components analysis, or - in the case
of dissimilarity input data - principal coordinates analysis; and `pa' for
all pairwise plots.  Note however that this function tries to judge whether
or not these representations are appropriate (cf. below).
.SH VALUE
the partition, as input, is returned invisibly.
.SH SIDE EFFECTS
a plot window is opened and a plot (see below) is produced.  
.SH DETAILS
User choice as regards representation may be overridden on the grounds of
(1) whether the input data set is one-dimensional or multidimensional; and
(2) whether the input data was coordinate data or dissimilarity data.  The
latter is indicated using attributes associated with the clustering
object passed to this function.  Note that we do not test for the existence
of such data sets: they are assumed to be present in a working directory on
the search path.
.PP
The options implemented are as follows.
.PP
(1) If the original data is coordinate data, and is one-dimensional, then
box-plots are used in the plotting.
.PP
(2) If the original data is coordinate data, and is multidimensional, 
and `rep = pc' is specified, then the cluster sequence numbers are plotted
in the principal plane of a principal components analysis.
.PP
(3) If the original data is dissimilarity data, and is multidimensional, and
`rep = pc' is specified, then the dissimilarities are processed by 
principal coordinates analysis (or classical multidimensional scaling; 
S function `cmdscale') and the cluster sequence numbers of the clusters'
members are plotted in the principal plane.  Note that principal coordinates
analysis is based on a Euclidean distance perspective: in a future version
of this function it is intended to also have non-metric multidimensional
scaling available as an option.
.PP
(4) If the original data is coordinate data, and is multidimensional, and if
`rep = pa' is specified, then a set of all pairwise plots is produced,
with labeling based on cluster sequence numbers.
.PP
(5) If the original data is coordinate data, and if the number of clusters in
the input partition is less than or equal to 4, and if the dimensionality of
the input data is less than or equal to 4, then box plots of the clusters'
members are produced.  This 5th option may be offered, under the appropriate
circumstances, in addition to other options.
.PP
(6) If the original data has resulted from a dimensionality-reduction 
technique (i.e. of class `"reddim"', produced by `pca', `ca' or `sammon'),
then the `rproj' component is used instead of coordinate data.
.PP
In all cases the user is prompted as regards closing the plot window.  If the
user so wishes, the plot window remains active.
.SA
plot, partition, print, summary
.EX
# Produce a partition with 4 clusters.  Input data set is `a'.
clresult1 <- partition(a, ng=4)
# Next, a hierarchical clustering.
clhier    <- hierclust(a)
# Derive the 4-cluster partition from this.
clresult2 <- partition(clhier, 4)
# See how the results compare.  (Different information is available from
# `print' and `summary'.)
plot(clresult1)
plot(clresult2)
# Next a Sammon mapping.
mds <- sammon(iris.var, tol=0.05, maxit=200)
mds.part <- partition(mds, 3)
plot(mds.part)
# In the foregoing, there are some evaluation/recursive call problems if one 
# does: plot(partition(sammon(iris.var, tol=0.05, maxit=200), 3))  
# To be corrected in a future version...
.KW multivariate
.KW clustering
.WR
