.BG
.FN plot.hierarchy
.TL
Plot a Hierarchy Object
.CS
plot(x, ...)
.RA
.AG x
object of class `hierarchy'.  Such an object is produced by the function
`hierclust' (which incorporates S function `hclust') and `modclust'
(which incorporates `mclust').
.SH VALUE
the hierarchy, as input, is returned invisibly.
.SH SIDE EFFECTS
the dendrogram is plotted on an already-opened, active plot window.
.SH DETAILS
The S function `plclust' is used, together with attributes of the hierarchy
being properly handled.
.SA
plot, plclust, hierclust, modclust
.EX
# Produce a hierarchical clustering and plot it.  
# Note: the plot window must be open.
clhier    <- hierclust(a)
motif()
plot(clhier)
.KW multivariate
.KW cluster
.WR
