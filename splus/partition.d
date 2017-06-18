.BG
.FN partition
.TL
Partitioning by Iterative Optimization; or through Slicing a Dendrogram
.DN 
Returns cluster memberships.  If input is a hierarchy, then the dendrogram
is cut to yield the requested number of clusters in the partition.  Otherwise,
in the case of a data matrix, the user specifies either a requested number of 
clusters; or - implicitly - as many clusters as there are rows in an initial 
guess at the cluster centers.  Output is an object of class `clustering'.  
Incorporates S function `kmeans'.
.CS
partition(a, ng, iter.max=15, option="mindst", diagnostics=F)
.sp
or
.sp
partition(a, centers, iter.max=15, diagnostics=F)
.sp
or
.sp
partition(a, ng, showslice=F, diagnostics=F)
.RA
.AG a
in the first two of the three examples of usage above, matrix of multivariate 
data.
Each row corresponds to an observation, and each column corresponds to 
a variable.  If `a' is of class `"reddim"' (as yielded by principal components
analysis, `pca', or correspondence analysis, `ca', or Sammon mapping, 
`sammon'), then the `rproj' component is used.
In the third of the three examples of usage above, a hierarchy
object as produced by `hierclust' or `modclust'.  
.AG ng
Number of groups or clusters.
.AG centers
matrix of initial guesses for the cluster centers.  Each row represents a
cluster center, and thus `centers' must have the same number of columns as
`a'.  The number of rows in `centers' is the number of clusters that will be
formed.
.OA
.AG iter.max
maximum number of iterations.
.AG option
Either `mindst' (default) or `exch'.  Options for the Spaeth (1985) 
algorithm. In the former case, the variances of
the groups are optimized by assigning the objects to groups such that they
are minimally distant from group centers.  In the latter case, the variances
are optimized by exchanging the objects between groups such that they are
minimally distant from group centers. 
.AG diagnostics
FALSE (default) implies that cluster cardinalities, cluster center coordinates,
and cluster compactness values will not be output to the command window.  
.AG showslice
whether or not the partition will be indicated with a dotted line drawn through
the dendrogram.  It this parameter is TRUE, then it is assumed that the
dendrogram has been plotted in the active plot window.
.RT
in the case of the Hartigan algorithm (where initial cluster centers 
are provided) a list with the following components:
.RC cluster
vector of integers, ranging from `1' to `ng' (or number of rows in 
`centers'), with length the same as the number of rows of `a'.  
The `i'th value indicates the group in which the
`i'th observation belongs.
.RC centers
matrix containing the coordinates of the final group centers.
.RC withinss
vector of length `ng' or number of rows in `centers'.  The `i'th value 
gives the within group sum of squares for the `i'th group.
.RC size
vector of length `ng', or the number of rows in `centers'.  
The `i'th value gives the number of observations in group `i'.
.PP
In the case of the Spaeth algorithm (where the number of clusters is given),
or in the case of the slicing of a dendrogram, the first of these, only, is
returned.  
.PP
In all cases, the returned object has a number of attributes:
.Co class ,
the class of the object, `"clustering"'.
.Co call ,
the calling statement which produced the object.
.Co origdata ,
the input data for this function.
.Co origdist ,
the distance measure used, currently `"Euclidean distance"'.
.SH METHOD
Consider first the clustering of coordinate data.
The object is to find a partition of the observations with `ng' (or the 
number of rows in `centers') groups that
minimizes the sum of `withinss' values.  To actually guarantee the 
minimum would be computationally
infeasible in many settings; this function finds a local minimum, that is,
a solution such that there is no single switch of an observation from one
group to another group that will decrease the objective criterion.  In the
case of the Spaeth (1985) algorithms, the local
minimum arrived at is dependent on the initial arbitrary assignment of 
observations to groups.  This arbitrary assignment is carried out using a
random number generator.  Subsequent executions of the `partition' function
will make use of different initial arbitrary assignments.  The function
`partition' should therefore be run a number of times, and the best result
kept.
.PP
The initial data is not standardized: it may be advisable to divide all column
values by the range of the column's values, or to standardize in some other 
way.
.PP
Note that this routine is best used by specifying a small number of groups;
specifying a large number of groups may cause many empty classes.  When this
happens, group centroids and compactness values may be undefined.  This may
also happen when there are very large values in the input data set.  
Cluster memberships, though, can still be used.  
.PP
Sample timings for the Spaeth (1985) algorithms (accessed by specifying the 
desired number of groups): 33000 rows, 4 columns, 3 groups took about 40 
secs. on a Sun SPARCstation 1; 50 groups took about 140 secs.  
.PP
When deciding on the number of clusters, Hartigan (1975, pp. 90-91) suggests
the following rough rule of thumb.  If `k' is the result of `partition' with
`k' groups, and `kplus1' is the result with `k+1' groups, then it is 
justifiable to add the extra group when 
.sp
(sum(k$withinss)/sum(kplus1$withinss)-1)*(nrow(a)-k+1)
.sp
is greater than 10.
.PP
When a partition is being obtained from a dendrogram, the procedure is
straightforward: a slice is made, orthogonal to the cluster criterion value
axis, and the clusters read off.
.SH REFERENCES
Spaeth, H. (1985).  
.ul
Cluster Dissection and Analysis: Theory, Fortran Programs, Examples.
Ellis Horwood, Chichester.
.sp
Hartigan, J.A. (1975).
.ul
Clustering Algorithms.
Wiley, New York.
.sp
Hartigan, J.A. and Wong, M.A. (1979).
`A k-means clustering algorithm',
.ul
Applied Statistics,
vol. 28, pp. 100-108.
.SA
`kmeans' (incorporated into this function); `hierclust'
(hierarchical clustering routines, incorporating `hclust'); `modclust'
(hierarchical clustering, incorporating `mclust').  To 
display results, `summary', `plot', etc.
.EX
# Firstly, specifying the desired number of groups.
ng <- 3
pp <- partition(a,ng)
# Plot the results in the plane comprising the first two columns of `a'
x <- a[,1]       # for convenience
y <- a[,2]
plot(x, y, type="n")     # set up plot scales, etc.
points(x[pp==1], y[pp==1], pch="*")
points(x[pp==2], y[pp==2], pch=".")
points(x[pp==3], y[pp==3], pch="o")
# Plot the results in the plane of the first two principal components
# using routine prcomp
pra <- prcomp(a)$x[,1:2]
plot(pra, type="n")      # set up plot scales, etc.
points(pra[pp==1,], pch="*")
points(pra[pp==2,], pch="+")
points(pra[pp==3,], pch="o")
.sp
# Secondly, specifying guesses at cluster centers.
irismean <- t(apply(iris, c(2, 3), 'mean'))
x <- rbind(iris[,,1], iris[,,2], iris[,,3])
km <- partition(x, irismean)
wrong <- km$cluster!=rep(1:3, c(50, 50, 50))
spin(x, highlight=wrong)
plot(x[,2], x[,3], type="n")
text(x[!wrong, 2], x[!wrong, 3], km$cluster)
# identify cluster membership that is correct
points(x[wrong, 2], x[wrong, 3], pch=15)
# boxes for points in error
title(main="K-Means Clustering of the Iris Data")
.KW clustering
.KW partition
.KW multivariate
.KW optimization
.WR

