.BG
.FN modclust
.TL
Model-Based Hierarchical Clustering
.DN
Carries out a model-based hierarchical clustering.  Refer to help for
`mclust'.  This function packages the latter with appropriate attributes 
for plot, summary and other methods.
.CS
modclust(a)
.RA
.AG a
matrix.  Cf. help for `mclust' for other parameters.
.RT
an object of class `hierarchy' (including `merge', `height' and `order'
components). In addition: `awe' component.
.PP
The returned object has a number of attributes:
.Co class ,
`hierarchy'
.Co call ,
calling statement
.Co criterion ,
agglomerative method used
.Co origdata
original input data set name
.Co origdist
currently set to "weighted Euclidean distance"
.SA
`mclust', `mreloc', `mclass'
.KW multivariate
.KW cluster
.WR

