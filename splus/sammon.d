.BG
.FN sammon
.TL
Sammon's Non-Linear Mapping
.DN
Finds a new, reduced-dimensionality, coordinate system for multivariate data 
such that the an error criterion between distances in the given space, 
and distances in the result space, is minimized.
.CS
sammon(a, p=2, maxit=100, tol=0.05, alpha=0.3, diagnostics=F)
.RA
.AG a
input matrix, coordinate data. No missing values.
.OA
.AG p
dimensionality of output space.
.AG tol
tolerance on error criterion.
.AG alpha
step size for gradient descent optimization.
.AG diagnostics
whether or not error value is output at each step of the iterative 
optimization.
.RT
projections in the new space.
.RC rproj
matrix of projections of the row observations, as yielded by the nonlinear
mapping algorithm.
.SH NOTE
It may be useful to run this routine a number of times and to keep the
result yielding the smallest error.  This mapping error, and the number of
iterations required for convergence, are output to the command window.
The object returned by this function is of class `"reddim"'.
.SH REFERENCES
W. Siedlecki, K. Siedlecka and J. Sklansky, 
An overview of mapping techniques for exploratory pattern analysis, 
.ul
Pattern Recognition,
21, 411-429, 1988.
.sp
J.W. Sammon, 
A nonlinear mapping for data structure analysis,
.ul
IEEE Trans. Computers, C-18, 401-409, 1969.
.EX
# Assume the 150x4 iris data matrix is in `iris.var'.  
mds <- sammon(iris.var, tol=0.05, maxit=200)
# Now plot observations 1-50, 51-100 and 101-150 distinctively;
# add a set of axes through x=0 and y=0.
plot(mds$rproj[,1], mds$rproj[,2], type="n", xlab="Axis 1", ylab="Axis 2",
main="2-d Sammon mapping of iris data")
points(mds$rproj[1:50,1], mds$rproj[1:50,2], pch="*")
points(mds$rproj[51:100,1], mds$rproj[51:100,2], pch="+")
points(mds$rproj[101:150,1], mds$rproj[101:150,2], pch="o")
plaxes(mds$rproj[,1], mds$rproj[,2])
.KW multivariate
.WR


