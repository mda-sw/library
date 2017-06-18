.BG
.FN flou
.TL
Fuzzy Coding (3-Way)
.DN
Simple fuzzy, or piecewise linear, coding of a vector: each value in the 
vector is replaced by a 1 (if it is above or equal to the 67th quantile),
by a 0 (if it is below or equal to the 33rd quantile), and by a linearly
interpolated value between 0 and 1 (if it lies between the 33rd and 67th
quantiles).  
.CS
flou(a)
.RA
.AG a
real-valued vector, with no missing values.
.RT
matrix of `length(a)' rows, and two columns.  The first column contains the
fuzzily coded values of `a', and the second column contains their 
complements.  Hence each row of this returned matrix necessarily sums to 1.
.SH BACKGROUND
This form of coding is suitable for a subsequent correspondence analysis.
When all variable have been fuzzily (or logically) coded, the row masses 
(proportional to the row sums) are identical.  
.SH REFERENCES
J.-P. Benzecri
.ul
Correspondence Analysis Handbook
Marcel Dekker, Basel, 1992.
.sp
F.J. Gallego,
Codage flou en analyse des correspondances,
.ul
Les Cahiers de l'Analyse des Donnees
vol. VII, 413-430, 1982
.SA
`logique', `ca', `supplr', `supplc'.  
.EX
# Fuzzy coding of input variables, `a', `b', `c':
a.fuzz <- flou(a)
b.fuzz <- flou(b)
c.fuzz <- flou(c)
newdata <- cbind(a.fuzz, b.fuzz, c.fuzz)
ca.newdata <- ca(newdata)
.KW multivariate
.KW algebra
.WR

