.BG
.FN ca
.TL
Correspondence Analysis
.DN
Finds a new coordinate system for multivariate data such that the first
coordinate has maximal inertia, the second coordinate has maximal inertia
subject to being orthogonal to the first, etc.  Compared to Principal 
Components Analysis, each row and column point has an associated mass 
(related to the row or column totals); and the chi-squared distance takes the 
place of the Euclidean distance.  The issue of how to code the input data
is important: this takes the place of input data transformation in PCA.
.CS
ca(a)
.RA
.AG a
data matrix to be decomposed, the rows representing observations and the
columns variables.
.AG nf
number of factors or axes to be sought; default 7.  
.RT
list of class `"reddim"' describing the correspondence analysis:
.RC rproj
projections of row points on the factors.
.RC cproj
projections of column points on the factors.
.RC evals
eigenvalues associated with the new factors. These provide figures of merit for
the "inertia explained" by the factors.  They are usually quoted in terms
of percentage of the total, or in terms of cumulative percentage of the total.
.RC evecs
definition of the factors in terms of the original variables. 
The first column is the linear combination of columns of 
`a' defining the first factor, etc.  
.RC rcntr
contributions of observations to the factors.  The contributions are mass times
projection (on the factor) squared.  Since contributions take account of the 
mass, they more accurately indicate influential observations for the 
interpretation of the factor, compared to the projections alone.
.RC ccntr
contributions of variables to the factors. See above remark concerning 
row contributions.
.SH NOTE
Very small negative eigenvalues, if they arise, are an artifact of the SVD
algorithm used, and are to be treated as zero.
.SH METHOD
A singular value decomposition is carried out.  
.SH BACKGROUND
Correspondence analysis defines the axis which provides the best fit
to both the row points and the column points.  A second axis is determined
which best fits the data subject to being orthogonal to the first.  Third and
subsequent axes are similarly found.  Best fit is in the least squares sense,
relative to the chi-squared distance.  This can be viewed as a weighted
Euclidean distance between `profiles'.
.PP
The question of `coding' of input data is an important one.  For instance,
in a matrix of scores, one might wish to adjoin extra columns to the input
matrix such that both the initial score, and the maximum score minus it,
are included in the observation's set of values.  Note that this has the
effect that all row masses are equal.  Hence the variables alone are 
differentially weighted.  This is known as `doubling' the observations.
In the case of binary data, such coding is known as 
`complete disjunctive form'.  
.PP
Other forms of input data for which correspondence analysis can be used include
frequencies, or contingency-type data.  In this case, the totaled chi-squared
distances of all (row or column) points from the origin is the familiar
chi-squared statistic. Hence the graphical output of correspondence analysis
allows assessment of departure from a null hypothesis of no dependence of 
rows and columns.
.PP
Supplementary rows or columns are projected into the factor space, after
carrying out a correspondence analysis.  That is to say, such row or 
column profiles are assumed to have zero mass, and their projections are
to be found under such an assumption.  Functions `supplr' and `supplc' may
be used for this purpose.  Supplementary rows or columns are of a different
nature compared to the basis data analyzed (e.g. sex in the context of a 
questionnaire); or they are rows or columns which, one suspects, would 
untowardly influence the definition of the factors.
.SH REFERENCES
Extensive works of J.-P. Benzecri including
.ul
Correspondence Analysis Handbook
Marcel Dekker, Basel, 1992.
.sp
M.J. Greenacre,
.ul
Theory and Applications of Correspondence Analysis
Academic Press, New York, 1984.
.sp
L. Lebart, A. Morineau and K.M. Warwick,
.ul
Multivariate Descriptive Statistical Analysis
Wiley, New York, 1984.
.sp
S. Nishisato, 
.ul
Analysis of Categorical Data: Dual Scaling and Its Applications
University of Toronto Press, Toronto, 1980.
.sp
(An extensive annotated bibliography is to be found in Greenacre.)
.SA
Supplementary rows and columns: `supplr', `supplc'.  Initial data coding:
`flou', `logique'.  Other functions producing objects of class "reddim":
`pca', `sammon'.  Other related functions: `prcomp', `cancor', `cmdscale'.
Plotting tool: `plaxes'.
.EX
# correspondence analysis of the breakfast cereal data, 
# in complete disjunctive form:
bfpos <- t(cereal.attitude)
bfneg <- max(bfpos) - bfpos
bfposneg <- cbind(bfpos, bfneg)
corr <- ca(bfposneg)
# plot of first and second factors
plot(corr$rproj[,1], corr$rproj[,2],type="n")
text(corr$rproj[,1], corr$rproj[,2], labels=dimnames(bfposneg)[[1]])
# Place additional axes through x=0 and y=0:
plaxes(corr$rproj[,1], corr$rproj[,2])
# check of row contributions
corr$rcntr
#
# Fuzzy coding of input variables, `a', `b', `c':
a.fuzz <- flou(a)
b.fuzz <- flou(b)
c.fuzz <- flou(c)
newdata <- cbind(a.fuzz, b.fuzz, c.fuzz)
ca.newdata <- ca(newdata)
.KW multivariate
.KW algebra
.WR

