.BG
.FN distance
.TL
Distance Matrix Calculation
.DN
Returns a distance structure that represents all of the pairwise distances
between objects in the data.
The choices for the metric are currently restricted to 
`"euclidean"', `"maximum"', `"manhattan"', and `"binary"'.
.CS
distance(x, metric="euclidean")
.RA
.AG x
matrix (typically a data matrix).  The distances computed
will be among the rows of `x'.  Missing values (`NA's) are
allowed.
.OA
.AG metric
character string specifying the distance metric to be used.
The currently available options are `"euclidean"',
`"maximum"', `"manhattan"', and `"binary"'.
Euclidean distances are root sum-of-squares of differences,
`"maximum"' is the maximum difference, `"manhattan"' is the sum
of absolute differences, and `"binary"' is the proportion
of non-zeros that two vectors do not have in common (the number 
of occurrences of a zero and a one, or a one and a zero
divided by the number of times at least one vector has a one).
.RT
the distances among
the rows of `x'.  Since there are many distances
and since the result of `dist' is typically an argument to
`hclust' or `cmdscale', a vector is returned, rather than a symmetric matrix.
For `i' less than `j', the distance between row `i' and row `j' is 
element `nrow(x)*(i-1) - i*(i-1)/2 + j-i' of the result.
.PP
The returned object has a number of attributes:
.Co Size ,
giving
the number of objects, that is, `nrow(x)'.
The length of the vector that is returned is `nrow(x)*(nrow(x)-1)/2',
that is, it is of order `nrow(x)' squared.
.Co class ,
indicating the type of data returned.  Currently only `"disddat"' is
supported, indicating dissimilarity data.  Note that this latter is taken
here as a class, encompassing similarities and distances as well.
.Co origdata ,
giving the name of the original input data.
.Co metric ,
giving the option chosen for the `metric' parameter.
.DT
Missing values in a row of `x' are not included in any distances
involving that row.
If the metric is `"euclidean"' and `ng' is the number of columns in which
no missing values occur for the given rows, 
then the distance returned is `sqrt(ncol(x)/ng)' times the Euclidean
distance between the two vectors of length `ng' 
shortened to exclude `NA's.
The rule is similar for the `"manhattan"' metric, except that the coefficient
is `ncol(x)/ng'.
The `"binary"' metric excludes columns in which either row has an `NA'.
If all values for a particular distance are
excluded, the distance is `NA'.
.SH NOTE
This function incorporates the S function `dist'.  The only difference in the
present implementation is additional set of attributes associated with the
function's output.
.PP
If the columns of a matrix are in different units, it is usually advisable
to scale the matrix before using `distance'.  A column that is much more 
variable than the others will dominate the distance measure.
.SH BACKGROUND
Distance measures are used in cluster analysis and in multidimensional scaling.
The choice of metric may have a large impact.
.SH REFERENCES
Everitt, B. (1980).
.ul
Cluster Analysis (second edition). 
Halsted, New York.
.sp
Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979).
.ul
Multivariate Analysis.
Academic Press, London.
.SA
`cmdscale', `hierclust', `hclust', `scale'.
.EX
dist(x,"max") # distances among rows by maximum
dist(t(x)) # distances among cols in Euclidean metric

# below is a function that converts a distance structure to a matrix
dist2full <- function(dis)
{
        n <- attr(dis, "Size")
        full <- matrix(0, n, n)
        for(rowi in seq(n - 1))
                for(colj in seq(from = rowi + 1, to = n)) {
                        full[rowi, colj] <- full[colj, rowi] <- dis[n * (rowi -
                                1) - (rowi * (rowi - 1))/2 + colj - rowi]
                }
        full
}
.KW multivariate
.KW cluster
.WR

