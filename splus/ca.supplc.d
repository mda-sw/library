.BG
.FN supplc
.TL
Supplementary Columns in Correspondence Analysis
.DN
Using the results of a correspondence analysis, project new columns into the
factor space.
.CS
supplc(a, ca.res)
.RA
.AG a
data matrix to be projected.  Must have same number of rows as matrix which
was initially input to the correspondence analysis.
.AG ca.res
the output of a correspondence analysis.  The following components of this
object are used: `evals', `rproj' and `cproj'.
.RT
a matrix, 
projections of the columns of `a' on the correspondence analysis factors.
.SH REFERENCES
See function `ca'.
.SA
Correspondence analysis: `ca'. 
Supplementary rows and columns: `supplr', `supplc'.  Initial data coding:
`flou', `logique'.  Other functions producing objects of class "reddim":
`pca', `sammon'.  Other related functions: `prcomp', `cancor', `cmdscale'.
Plotting tool: `plaxes'.
.EX
cares <- ca(logarray)
newproj <- supplc(newcols, cares)
# plot of first and second factors, and of supplementary columns:
plot(cares$cproj[,1], cares$cproj[,2],type="n")
text(cares$cproj[,1], cares$cproj[,2])
points(newproj[,1], newproj[,2])
# Place additional axes through x=0 and y=0:
plaxes(cares$cproj[,1], cares$cproj[,2])
.KW multivariate
.KW algebra
.WR

