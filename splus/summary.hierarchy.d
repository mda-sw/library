.BG
.FN summary.hierarchy
.TL
Summary Method for Objects of Class Hierarchy
.CS
print(summary(x), ...)
.RA
.AG x
object of class `hierarchy'.  Such an object is produced by the function
`hierclust' (which incorporates S function 'hclust') and `modclust'
(which incorporates `mclust'). 
.SH VALUE
the input parameter is returned invisibly.
.SH DETAILS
Pretty-print some of the summary information concerning a hierarchy object.
Included are:
number of observations on which the partition is based, agglomerative
method used,
original data set name, and dissimilarity used.
.SA
print, summary, hierclust, modclust, plot
.EX
summary(hierclust(a))
.KW multivariate
.KW cluster
.WR

