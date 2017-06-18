# 'source' and 'dyn.load' everything in Cluster Analysis Suite

cat("Installing Cluster Analysis Suite...\n")
cat("F. Murtagh, January 1993.  Email: fmurtagh@eso.org\n")
cat("F. Murtagh, January 1999.  Email: fmurtagh@acm.org\n")

cat(" \n")
cat("(1) We will compile eight Fortran programs in this directory (f77 -c).\n")
cat("(2) We will source S-Plus script files (extension .q).\n")
cat("(3) We will dyn.load all six object files.\n")
cat("(4) We will delete *all* object files (extension .o).\n") 
cat("(5) We will attempt to create directory /.Data/.Help (ineffective if it already exists).\n")
cat("(6) We will put all help files into the latter directory.\n")
cat("Okay?  If you wish, say no, and make changes in install.q.\n")
cat("Continue [y/Y]  ")
doit <- F
resp <- scan("", what=character(),1)
if (length(resp)==0 || resp=="y" || resp =="Y") doit <- T

if (doit == F) stop("Okay, we're stopping.\n")

cat("Compiling all Fortran programs...\n")
!f77 -c bea.f
!f77 -c hc.f
!f77 -c hcmovie.f
!f77 -c pca.f
!f77 -c members.f
!f77 -c partition.f
!f77 -c sammon.f
!f77 -c ca.f

cat("Sourcing all S-Plus scripts...\n")
source("bea.q")
source("distance.q")
source("pca.q")
source("plaxes.q")
source("members.q")
source("partition.q")
source("modclust.q")
source("plot.clustering.q")
source("print.clustering.q")
source("summary.clustering.q")
source("plot.hierarchy.q")
source("print.hierarchy.q")
source("summary.hierarchy.q")
source("sammon.q")
source("ca.q")
source("ca.supplr.q")
source("ca.supplc.q")
source("flou.q")
source("logique.q")
source("hc.q")

cat("Dyn.load'ing all object files...\n")
dyn.load("bea.o")
dyn.load("hc.o")
dyn.load("hcmovie.o")
dyn.load("pca.o")
dyn.load("members.o")
dyn.load("partition.o")
dyn.load("sammon.o")
dyn.load("ca.o")

NEW:
dyn.open("bea.o")
dyn.open("hc.o")
dyn.open("hcmovie.o")
dyn.open("pca.o")
dyn.open("members.o")
dyn.open("partition.o")
dyn.open("sammon.o")
dyn.open("ca.o")

# cat("Done.  Deleting all object files...\n")
# !rm *.o

cat("Now copy help files to ~/.Data/.Help\n")
cat("Attempt to create this directory in case it does not already exist.\n")
!mkdir ~/.Data/.Help
!cp      bea.d        ~/.Data/.Help/bea
!cp      distance.d   ~/.Data/.Help/distance
!cp      hierclust.d  ~/.Data/.Help/hierclust
!cp      members.d    ~/.Data/.Help/members
!cp      modclust.d   ~/.Data/.Help/modclust
!cp      partition.d  ~/.Data/.Help/partition
!cp      pca.d        ~/.Data/.Help/pca
!cp      plot.clustering.d   ~/.Data/.Help/plot.clustering
!cp      plot.hierarchy.d    ~/.Data/.Help/plot.hierarchy
!cp      print.clustering.d  ~/.Data/.Help/print.clustering
!cp      print.hierarchy.d   ~/.Data/.Help/print.hierarchy
!cp      summary.clustering.d ~/.Data/.Help/summary.clustering
!cp      summary.hierarchy.d  ~/.Data/.Help/summary.hierarchy
!cp      sammon.d       ~/.Data/.Help/sammon
!cp      ca.d           ~/.Data/.Help/ca
!cp      ca.supplr.d    ~/.Data/.Help/supplr
!cp      ca.supplc.d    ~/.Data/.Help/supplc
!cp      flou.d         ~/.Data/.Help/flou
!cp      logique.d      ~/.Data/.Help/logique

cat("Done.  \n")

