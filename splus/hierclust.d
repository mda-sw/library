.BG
.FN hierclust
.TL
Hierarchical Clustering
.DN
Agglomerates the closest pair of observations, and replaces them with
a cluster.  Agglomerates the next closest pair of observations or clusters.
Continues through `n-1' agglomerations, assuming there are `n' observations.
The output is an object of class `hierarchy'.
.CS
hierclust(a, method=1, bign = F, diagnostics=F, sim=" ",
             movie=F, option="thru",alpha=1.0,repres="chull", 
             show="change")
.sp
or
.sp
hierclust(dist, method="compact", sim=)
.RA
.sp
First variant:
.sp
.AG a
data matrix to be analyzed, the rows representing observations and the
columns variables; or, alternatively, lower-diagonal dissimilarity matrix in
vector or list format.  In latter case, the number of values is `n(n-1)/2' 
where `n' is the number of observations.
.OA
.sp
First variant:
.sp
.AG method
First possibility: integer taking the value between 1 and 7.  These methods 
correspond to:
Ward's minimum variance or error sum of squares method (1); single linkage
or nearest neighbor method (2); complete linkage or diameter (3); average
linkage, group average, or UPGMA method (4); McQuitty's or WPGMA method
(5); median, Gower's or WPGMC method (6); and centroid or UPGMC method (7).
If the `movie' option is used, only `method's 1 to 4 are supported.
Second possibility: string `compact' (default; complete linkage), 
`average' (average linkage), or `connected' (single linkage).
.AG sim
structure giving similarities rather than dissimilarities.  This can either
be a symmetric matrix or a vector with a `Size' attribute.  Only one of 
dissimilarities or similarities may be specified (or, of course, as an
alternative, a data matrix.
.AG bign
is `n' big?  If storage is problemsome, a different implementation of the 
Ward criterion may be tried.  This determines dissimilarities on the fly, and
hence requires `O(n)' storage.
.AG diagnostics
show, or suppress, some information directed to the command window. 
Suppression is useful if functional programming capabilities are availed of,
leading to lazy evaluation; this often results in multiple evaluations, 
with consequent repetition of diagnostic information.
.AG movie
whether interactive display of construction of hierarchy is wanted.  This
necessitates an open plot window.  It does not produce an output object.  
All remaining arguments apply only in the case of the `movie=T'.
.AG option
(Only for `movie'.) "thru" (default) or "prompt": Go right thru all n-1 
agglomerations, or ask
the user whether to continue at each agglomeration?  Here is what we
recommend: On the first occasion, go right through.  Check out the cluster
criterion values reported on in the command window.  Then using the 'prompt'
option, go as far as the desired partition, and print it out.
.AG alpha
(Only for `movie'.) 1.0 (default) or 0.0 < alpha <= 1.0.
Only used by the minimum variance agglomerative criterion. Alpha downweights
the variance component on the principal axis of the cluster, and thereby
redefines the cluster criterion such that linear clusters are preferred.  A
value of alpha = 1.0 gives exactly the traditional Ward's minimum variance
method.  See Murtagh and Raftery (1984) for more on this altered criterion.
.AG repres
(Only for `movie'.) repres = "chull" (default), "lines", "gauss":
Representation of clusters: For input data which has dimensionality greater
than 2, only "chull" is allowed (since the representations were found to be
too problematic, given that the clustering takes place in m-dimensions, and
the the representation is necessarily a 2-dimensional projection of this).
Representation "chull" = convex hull, or straight line in case
of collinear points.  (We test for collinearity using principal
components analysis.)
Representation "lines" = lines of best fit (least squares fit is used,
and the extremities are defined by the intersection of the min. and max. x
values on this LS fit line). "gauss" = 1-sigma circles.  Side-effect of
latter: warning messages indicate whenever circles overlap plot boundaries.
.AG show
(Only for `movie'.) "change" or "all":
In this method's present implementation, the clusters at each stage are
`highlighted' by means of a different symbol.  Option show = "change" has
this highlighting changed at each agglomeration, so that only the most recent
agglomerands are shown.  Option show = "all" does not turn off the
highlighting in this way.
.RA
.sp
Second variant (which re-wraps S function `hclust', adding appropriate 
attributes, including class `hierarchy': cf. help file on `hclust'):
.sp

Exactly one of dist or sim must be specified.
.sp
.AG dist
a distance structure or distance matrix.  Normally this
will be the result of the function `distance', but it can be any
data of the form returned by `distance', or  a  full,  symmetric
matrix.  Missing values are not allowed.
.OA
.AG method
character string giving  the  clustering  method.   The
three  methods  currently implemented are "average", "con-
nected" (single linkage) and "compact" (complete linkage).
(The first three characters of the method are sufficient.)
.AG sim=
structure giving similarities  rather  than  distances.
This  can  either be a symmetric matrix or a vector with a
"Size" attribute.  Missing values are not allowed.
.RT
When `movie' is not specified in the first variant, this is a 
list describing the hierarchical clustering of the observations.
The `movie' option is considerably slower than alternative options, due
to the processing involved.  Hence its use is recommended only on small
data sets (less than 100 observations).
.RC merge
an `n-1' by 2 matrix. Row `i' of `merge' describes the merging of clusters
at step `i' of the clustering.  If an element `j' in the row is negative,
then observation `-j' was merged at this stage.  If `j' is positive 
then the merge
was with the cluster formed at the (earlier) stage `j' of the algorithm. Thus
negative entries in `merge' indicate agglomerations of singletons, and 
positive entries indicate agglomerations of non-singletons.
.RC height 
a set of `n-1' non-decreasing real values.  The
clustering "height": that is, the value of the criterion associated with
the clustering `method' for the particular agglomeration.
.RC order
a vector giving the permutation of the original observations suitable for 
plotting,
in the sense that a cluster plot using this ordering and matrix `merge' will
not have crossings of the branches.
.PP
In hierarchical cluster displays, a decision is needed at each merge to
specify which subtree should go on the left and which on the right.  Since, 
for 
`n' observations there are `n-1' merges, there are `2^(n-1)' possible
orderings for the leaves in a cluster tree, or dendrogram.  The default
algorithm in `hc' is to order the subtree so that the tighter cluster is on
the left (the last, i.e. most recent, merge of the left subtree is at a lower 
value than the last
merge of the right subtree).  Observations are the tightest clusters possible,
and merges involving two observations place them in order by their 
observation sequence number.
.PP
Attributes associated with the returned hierarchy include: the class,
`"hierarchy"'; the original data name; and the metric used (the latter
being implemented, so far, mainly for Euclidean distance); the call which
created the object; and the agglomerative clustering criterion used.
.SH NOTE
The squared Euclidean distance is used in the `movie' case, and in the case
when a data matrix is input. Note that the `movie' option returns 
invisibly.   
.SH METHOD
In the case of all methods when `bign = F' (default), 
the Lance-Williams dissimilarity update formula is used to allow the 7 methods
to be implemented.   
.PP
In the direct data matrix case, 
a list of nearest neighbors is maintained, which 
allows the closest pair of clusters or observations to be determined by
scanning this list to find the pair of agglomerands on each iteration.
Hence `O(n)' operations are required for each of `n-1' iterations, i.e.
`O(n^2)' computation time in total.  The 
initial setting up of the `n' nearest neighbors also requires `O(n^2)' time.
The processing to be carried out following each agglomeration is `O(n)'.
In total, all methods implemented in `hc' are of time complexity `O(n^2)'.
Unless `bign = T' (see below), storage 
complexity is `O(n^2)'.  Typical times on a Sun SPARCstation
1 for 400x4 and 800x4 arrays: 12 and 48 secs.
.PP
Ward's minimum variance method aims at finding compact, spherical clusters.
The complete link method finds similar clusters. The single link method,
which is closely related to the minimal spanning tree, adopts a 
`friends of friends' clustering strategy.  
The other methods can be regarded as aiming
for clusters with characteristics somewhere between the single and complete
link methods. 
.PP
In the case of the Ward method, and `bign = T', the nearest-neighbor 
chain algorithm is used, which gives an `O(n^2)'
computation time algorithm. The input matrix is stored, but not the
dissimilarities produced.  Computation time is marginally inferior to the
implementation of Ward's method when `bign = F' (default).  Sample 
times for a 400x4 and an 800x4
array on a Sun SPARCstation 1: 12.5 and 50 secs.
.SH REFERENCES
A small sample: 
.sp
B. Everitt, 
.ul
Cluster Analysis
Heinemann Educ. Books, London 1974; 
.sp
P.H.A. Sneath and R.R. Sokal, 
.ul
Numerical Taxonomy
Freeman, San Francisco, 1973; 
.sp
M.R. Anderberg, 
.ul
Cluster Analysis for Applications
Academic Press, New York, 1973;
.sp
A.D. Gordon, 
.ul
Classification
Chapman and Hall, London, 1981; and
.sp
F. Murtagh, 
.ul
Multidimensional Clustering Algorithms
COMPSTAT Lectures 4,
Physica-Verlag, Wuerzburg, 1985 (for algorithmic details of algorithms used).
.SA
`plot', `summary', `print' to look at hierarchical clustering results.
'distances' to produce (dis)similarity input.  `members' may be used to 
obtain cluster assignments of all observations.  Function `hclust' is
included in `hierclust', with additing of appropriate attributes.  In a
similar way, function `distance' supersedes `dist'.
.EX
# Default Ward criterion
h <- hierclust(a)
plot(h)
# Real-time display of clustering
hierclust(a, movie=T)
# Single link criterion
h <- hierclust(distances(a), method=connected)
# Get cluster assignments corresponding to 5-cluster partition:
k <- members(h)
k$class[,4]
# Large `n'
h <- hierclust(a, bign=T)
# An example of the incorporated `hclust' function, which allows processing
# of dissimilarity data:
hierclust(distance(a))
.KW cluster
.KW multivariate
.WR
