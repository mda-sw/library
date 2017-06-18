.BG
.FN print.hierarchy
.TL
Print a Hierarcy Object
.CS
print(x, ...)
.RA
.AG x
object of class `hierarchy'.  Such an object is produced by function
`hierclust' (which incorporates function `hclust') and function
`modclust' (which incorporates function `mclust').
.SH VALUE
the hierarchy, as input, is returned invisibly.
.SH DETAILS
Pretty-print some of the attributes of a hierarchy object.  Included are:
sequence of merges, list of node heights, and "horizontal" ordering of 
observations.
.SA
print, summary, hierclust, modclust, plot
.EX
# Produce a hierarchical clustering and print some of its attributes.
clhier    <- hierclust(a)
print(clhier)
.KW multivariate
.KW cluster
.WR

