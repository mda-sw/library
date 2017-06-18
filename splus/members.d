.BG
.FN members
.TL
Cluster Memberships 
.DN
From a succession of agglomerations produced by a hierarchical routine
in function `hierclust' or `hclust', determine cluster assignments of all 
objects, at all levels of the hierarchy.
.CS
members(a)
.RA
.AG a
output produced by functions `hierclust' (or 'hclust') or `modclust'
(or `mclust').
.RT
matrix of dimensions `n' by `n-2' giving cluster assignments to the 2, 3, ...
n-1 cluster partitions of the hierarchy.  This corresponds to levels `n-2'
to 2.  The observations correspond to the rows, and are in
sequence, `1, 2, ...'.  Column `j' specifies which of the `2, 3, ... j-1'
clusters each observation is associated with at the level `j-1' of the 
agglomerative process.  The two clusters which merge in going from level 
`j-1' to `j-2' are replaced in the `j-2'nd column with the lower of the
two cluster sequence numbers.  
.SH NOTE
.PP
The time requirement of 'members' is `O(n^3)'.  Sample time for Sun 
SPARCstation 1, with 600 observations: 303 secs.
.SA
Functions `hierclust', `hclust', `modclust', or `mclust' produce the 
array containing the 
sequence of agglomerations which can be used as input for `members'.
.EX
h <- hierclust(a, method=2)
k <- members(h)
.KW cluster
.KW multivariate
.WR
