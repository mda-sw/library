import java.text.*;
import java.io.*;

/**  
 * Carry out pairwise agglomerations.  For n items, therefore there
 * are n-1 agglomerations.  Represent the cluster labels is an nxn 
 * cluster label matrix.  Column no. n will be the singleton labels, 
 * 1 to n.  Column no. n-1 will have n-1 unique values (or label sequence
 * numbers).  Column no. n-2 will have n-2 unique values.  Column no. 1 
 * will have the value 1 only, implying that all n items are in one 
 * cluster.
 * <p>
 * ClustMat is our agglomeration "engine".  It looks after labeling 
 * only, and is independent of any agglomerative clustering criterion.
 * <p>
 * Other utility methods:
 * <p>
 * Dissim ... calculate dissimilarity matrix <br>
 * getNNs ... get nearest neighbors and associated
 *            nearest neighbor dissimilarities <br>
 * getSpaces ... helping in output formating <br>
 * printMatrix ... print matrix of doubles, or integers <br>
 * printVect ... print vector of doubles, or integers
 * <p> 
 * main does the following:
 * <ol>
 * <li> Reads data. (Format: integer row, column dimensions, followed
 *      by matrix values, read row-wise.)  No preprocessing on input. 
 * <li> Calculate pairwise dissimilarities, and determines nearest neighbors 
 *       and corresponding dissimilarities.  (Squared Euclidean distance 
 *       used.) 
 * <li> Determines the closest nearest neighbors.
 * <li> Carries out an agglomeration in ClustMat.
 * <li> Updates the  pairwise dissimilarity matrix, and then, on the basis
 *  of this, the nearest neighbors, and the nearest neighbor 
 *  dissimilarities.
 * <li> Repeats while no. of clusters is greater than 2. 
 * </ol>
 * Constant MAXVAL is used in, resp., dissimilarities and 
 * nearest neighbor dissimilarities, to indicate when items are processed and 
 * no longer exist as singletons.  <br>
 * Note also how flag = 1 denotes an active observation, and flag = 0
 * denotes an inactive one (since it has been agglomerated).  It is not 
 * necessary to use the flag since exceptionally high dissimilarities 
 * will signify inactive observations.  However the use of flag is 
 * helpful computationally.  
 * <p>
 * Step 5 here determines the agglomerative clustering criterion.  We are
 * currently using the minimum variance method.  It is indicated
 * in the code where to change to use other agglomerative criteria. <p>
 * Output cluster labels using original sequence numbers.  The ordering
 * of observations is not such that a dendrogram can be directly 
 * constructed.  However the cluster labels do allow selections and
 * further inter and intra cluster processing of the input data.  <p>
 * Example of use: <p>
 * <tt> javac HCL.java </tt> <br>
 * <tt> java HCL <a href="../iris.dat">iris.dat</a> > 
 *           <a href="../hcloutput.txt">hcloutput.txt</a> </tt> <p>
 * First version: 1999 Nov. <br>
 * Version: 2002 Oct. 26 <br>
 * Author: F. Murtagh, f.murtagh@qub.ac.uk
 * @version 2002 Oct. 26
 * @author F. Murtagh, f.murtagh@qub.ac.uk
 */

public class HCL
{
    public static final double MAXVAL = 1.0e12;

    /** 
     * Method Dissim, calculates dissimilarity n x n array
     * @param nrow integer row dimension
     * @param ncol integer column dimension
     * @param A floating row/column matrix
     * @return Adiss floating n x n dissimilarity array 
     */
    public static double[][] Dissim(int nrow, int ncol, 
                  double[] mass, double[][] A)
    {
    double[][] Adiss = new double[nrow][nrow];

       // Initialize to 0.  Caters for the diagonal which stays 0.
       for (int i1 = 0; i1 < nrow; i1++)  {
         for (int i2 = 0; i2 < nrow; i2++) Adiss[i1][i2] = 0.0;
       }
       for (int i1 = 0; i1 < nrow; i1++)  {
	   // All dissimilarity we are dealing with are symmetric
	   // so just calculate for half the array and fill in later
	   for (int i2 = 0; i2 < i1; i2++) {
	     for (int j = 0; j < ncol; j++) {
	       // We are happy with the squared dissimilarity
	       // 0.5 term since this equals 
               // (clcard[i1]*clcard[i2])/(clcard[i1]+clcard[i2]) 
               // for minimum variance criterion
	       Adiss[i1][i2] += 0.5*Math.pow(A[i1][j] - A[i2][j], 2.0); }
	       // Adiss[i1][i2] += ((mass[i1]*mass[i2])/(mass[i1]+mass[i2]))*
	       //   Math.pow( A[i1][j] - A[i2][j], 2.0 ); 
  	       // Fill the other half of the array 
               Adiss[i2][i1] = Adiss[i1][i2];
	   }
        }
        return Adiss;
    } // Dissim 
           

    /**
     * Method getNNs, determine NNs and NN dissimilarities
     * @param nrow row dimension or number of observations (input)
     * @param flag =1 for active observation, = 0 for inactive one (input) 
     * @param diss dissimilarity matrix (input)
     * @param nn nearest neighbor sequence number (calculated)
     * @param nndiss nearest neigbor dissimilarity (calculated)
     */
   public static void getNNs(int nrow, int[] flag, double[][] diss,
           int[] nn, double[]nndiss)
    {
       int minobs;  
       double mindist; 
       for (int i1 = 0; i1 < nrow; i1++)  {
           if (flag[i1] == 1) {
              minobs = -1;
              mindist = MAXVAL;
              for (int i2 = 0; i2 < nrow; i2++)   {
	         if ((diss[i1][i2] < mindist) && (i1 != i2)) {
                    mindist = diss[i1][i2];
                    minobs = i2;
                  }
              }
              nn[i1] = minobs + 1;
              nndiss[i1] = mindist;
           }
       }
       // Return type void => no return stmt
    } // getNNs


    /**
     * Method ClustMat, updates cluster structure matrix following 
     * an agglomeration
     * @param nrow row dimension or number of observations (input)
     * @param clusters list of agglomerations, stored as array of 
     *        pairs of cluster sequence numbers (input, and updated)
     * @param clust1 first agglomerand (input)
     * @param clust2 second agglomerand (input)
     * @param ncl number of clusters remaining (input)
     */
    public static void ClustMat(int nrow, int[][] clusters, 
           int clust1, int clust2, int ncl)
    {
        // If either clust* is not initialized, then we must init. clusters
	if ((clust1 == 0) || (clust2 == 0)) {
           for (int j = 0; j < nrow; j++) {
               for (int i = 0; i < nrow; i++) {
                   clusters[i][j] = 0;
	       }
	   }
	   // Adjust for 0-sequencing in extreme right-hand col values.
           for (int i = 0; i < nrow; i++) clusters[i][ncl-1] = i + 1;
           return;
	}

        // For some agglomeration, we are told that label clust1 and 
        // label clust2, among all labels in col. ncl-1 (ncl ranges over
        // 0 thru nrow-1) are to be agglomerated
	// We have arranged that order is such that clust1 < clust2
	// Note: clust1 and clust2 are 1-sequenced and not 0-sequenced 

        int ncl1; 
        ncl1 = ncl - 1; 
        for (int i = 0; i < nrow; i++) {
            clusters[i][ncl1] = clusters[i][ncl];
            if (clusters[i][ncl1] == clust2) clusters[i][ncl1] = clust1;
        }
        // Return type void => no return stmt
    } // ClustMat
 

    // Little method for helping in output formating  
    public static String getSpaces(int n) {
    StringBuffer sb = new StringBuffer(n);
    for (int i = 0; i < n; i++) sb.append(' ');
    return sb.toString();
    }


    /**
     * Method for printing a double float matrix  <br>
     * Based on ER Harold, "Java I/O", O'Reilly, around p. 473.
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values, double 
     * @param d display precision, number of decimal places
     * @param w display precision, total width of floating value
     */
    public static void printMatrix(int n1, int n2, double[][] m, int d, int w)
    {
        // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      myFormat.setMaximumIntegerDigits(d);
      myFormat.setMaximumFractionDigits(d);
      myFormat.setMinimumFractionDigits(d);
      for (int i=0; i<n1; i++)
          {
              // Print each row, elements separated by spaces
              for (int j=0; j<n2; j++) {
                  String valString = myFormat.format(
                       m[i][j], new StringBuffer(), fp).toString();
                  valString = getSpaces(w - fp.getEndIndex()) + valString;
                  System.out.print(valString);
                  }
              // Start a new line at the end of a row
              System.out.println();
          }
      // Leave a gap after the entire matrix
      System.out.println();
    }


    /**
     * Method for printing an integer matrix  <br>
     * Based on ER Harold, "Java I/O", O'Reilly, around p. 473.
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values 
     * @param d display precision, number of decimal places
     * @param w display precision, total width of floating value
     */
    public static void printMatrix(int n1, int n2, int[][] m, int d, int w)
    {
        // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      myFormat.setMaximumIntegerDigits(d);
      for (int i=0; i<n1; i++)   {
              // Print each row, elements separated by spaces
              for (int j=0; j<n2; j++)  {
                  String valString = myFormat.format(
                       m[i][j], new StringBuffer(), fp).toString();
		  // 4 character locations per integer number 
                  valString = getSpaces(w - fp.getEndIndex()) + valString;
                  System.out.print(valString);
                  }
              // Start a new line at the end of a row
              System.out.println();
          }
      // Leave a gap after the entire matrix
      System.out.println();
    }  // printMatrix


    /**
     * Method printVect for printing a double float vector <br>
     * Based on ER Harold, "Java I/O", O'Reilly, around p. 473.
     * @param m input vector of length m.length
     * @param d display precision, number of decimal places
     * @param w display precision, total width of floating value
     */
    public static void printVect(double[] m, int d, int w)
    {
        // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      myFormat.setMaximumIntegerDigits(d);
      myFormat.setMaximumFractionDigits(d);
      myFormat.setMinimumFractionDigits(d);
      int len = m.length;
      for (int i=0; i<len; i++)  {
             String valString = myFormat.format(
                   m[i], new StringBuffer(), fp).toString();
             valString = getSpaces(w - fp.getEndIndex()) + valString;
                   System.out.print(valString);
          }
          // Start a new line at the row end
          System.out.println();
      // Leave a gap after the entire vector
      System.out.println();
    }  // printVect


    /**
     * Method printVect for printing an integer vector <br>
     * Based on ER Harold, "Java I/O", O'Reilly, around p. 473.
     * @param m input vector of length m.length
     * @param d display precision, number of decimal places
     * @param w display precision, total width of floating value
     */
    public static void printVect(int[] m, int d, int w)
    {
        // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      myFormat.setMaximumIntegerDigits(d);
      int len = m.length;
      for (int i=0; i<len; i++)  {
             String valString = myFormat.format(
                   m[i], new StringBuffer(), fp).toString();
             valString = getSpaces(w - fp.getEndIndex()) + valString;
                   System.out.print(valString);
          }
          // Start a new line at the row end
          System.out.println();
      // Leave a gap after the entire vector
      System.out.println();
    } // printVect


    // The main method contains the body of the program
    public static void main(String[] argv)
    {  
   PrintStream out = System.out;

   try{
     if (argv.length == 0) {
        System.out.println(" Syntax: java HCL infile.dat ");
        System.out.println(" Input file format: ");
        System.out.println(" Line 1: integer no. rows, no. cols.");
        System.out.println(" Successive lines: matrix values, floating");
        System.out.println(" Read in row-wise");
        System.exit (1);
     }
     String filname = argv[0]; 
     System.out.println (" Input file name: " + filname);

     // Open the matrix file
     FileInputStream is  = new FileInputStream(filname);
     BufferedReader  bis = new BufferedReader(new InputStreamReader(is));
     StreamTokenizer st  = new StreamTokenizer(bis);

     // Row and column sizes, read in first
     st.nextToken();  int nrow = (int)st.nval;
     st.nextToken();  int ncol = (int)st.nval;

     System.out.println(" No. of rows, nrow = " + nrow);
     System.out.println(" No. of cols, ncol = " + ncol);
     
     // Input array, values to be read in successively, float
     double[][] indat = new double[nrow][ncol]; 
     double inval;

     // New read in input array values, successively
     System.out.println
        (" Input data sample follows as a check, first 4 values.");
     for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
             st.nextToken();
             inval = (double)st.nval;
             indat[i][j] = inval;
             if (i < 2 && j < 2) {
                System.out.println(" value = " + inval);
             }
        }
     }  
     System.out.println();

     //-------------------------------------------------------------------

     // Some definitions and initializations
        int[][] clusters = new int[nrow][nrow];  // cluster label matrix
        int[] nn = new int[nrow];                // cluster nearest neigh's
        int[] flag = new int[nrow];              // flag of active obs
        double[] nndiss = new double[nrow];      // nearest neigh diss's
        double[] clcard = new double[nrow];      // cluster cardinality
	double[] mass = new double[nrow];        // observation masses
	double[] cpoids = new double[ncol];      // column weights
        int minobs;
        double mindist;  
        int ncl;   // Initial number of clusters 
        ncl = nrow; 

        for (int i = 0; i < nrow; i++) {
            flag[i] = 1;
            clcard[i] = 1.0;
            mass[i] = 1.0; 
        }
        for (int j = 0; j < ncol; j++) {
            cpoids[j] = 0.0; 
        }

        // We may wish to normalize out input data, 
	// and create mass and cpoids, resp. row and col. masses
	// In the following: row/col. masses are row/col sums div. by total
	// double tot = 0.0; 
        // for (int i = 0; i < nrow; i++) {
        //   for (int j = 0; j < ncol; j++) {
        //      mass[i] += indat[i][j];
	//	cpoids[j] += indat[i][j]; 
        //      tot += indat[i][j];
	//    }
	// }
	// And then we det. fij <-  xij/fi*sqrt(fj)
	// where fi and fj are, resp., row and col. masses
        // for (int j = 0; j < ncol; j++) {
	//  cpoids[j] = Math.sqrt(cpoids[j]/tot);
	// }
	// for (int i = 0; i < nrow; i++) {
	//  mass[i] /= tot; 
	//  for (int j = 0; j < ncol; j++) {
	//	indat[i][j] /= (tot*mass[i]*cpoids[j]);
	//    }
	// }
	//System.out.println("Normalized array to be analyzed:");
        //printMatrix(nrow, ncol, indat, 4, 10);
	//System.out.println("Row masses:");
        //printVect(mass, 4, 10);

	// Alternative normalizations include:
	// Centre row vectors to zero mean, and reduce to unit std. dev.
	// We'll use no normalization.  We'll take masses as equal to 
	// cluster cardinalities.  

        // Get dissimilarities
        double[][] diss = new double[nrow][nrow];
        diss = Dissim(nrow, ncol, mass, indat); 
        // Print it out
        System.out.println("Dissimilarity matrix for analysis:");
        printMatrix(nrow, nrow, diss, 4, 10);

	// Get nearest neighbors and nearest neighbor dissimilarities
        getNNs(nrow, flag, diss, nn, nndiss);
	//  System.out.println("Nearest neighbors:");
	//  printVect(nn, 4, 10); 
        //  System.out.println("Nearest neighbors dissimilarities:");
        //  printVect(nndiss, 4, 10);

        // Get closest neighbors, using nndiss, followed by nn
        int clust1 = 0;
        int clust2 = 0;
        int cl1    = 0;
        int cl2    = 0;
	// Call to initialize
        ClustMat(nrow, clusters, clust1, clust2, ncl);
        //  System.out.println("No. of clusters:  " + ncl);  
	//  System.out.println(" Cluster label matrix:");
        //  printMatrix(nrow, nrow, clusters, 4, 10);


	// Loop to carry out series of nrow-1 agglomerations
        do 
	    {

        minobs = -1;
        mindist = MAXVAL;
	    for (int i = 0; i < nrow; i++) {
                if (flag[i] == 1) {
                   if (nndiss[i] < mindist) {
                      mindist = nndiss[i];
                      minobs = i;
		   } 
                }
	    } 
	    // minobs is one cluster label, the other is nn[minobs]
	    // Adjust for 0-sequencing
  	    if (minobs < nn[minobs]) {
               clust1 = minobs + 1;
               clust2 = nn[minobs];
            }
  	    if (minobs > nn[minobs]) {
               clust2 = minobs + 1;
               clust1 = nn[minobs];
            }
	     //    Now we have clust1 < clust2, and we'll agglomerate
	    System.out.println(" clus#1: "  + clust1 +
                               ";  clus#2: "  + clust2 +
                               ";  new card: " + (clcard[clust1-1]+
                                                  clcard[clust2-1]) + 
                               "; # clus left: " + ncl +
                               "; mindiss: " + mindist);

	    // Now we will carry out an agglomeration
            ncl = ncl - 1; 
            ClustMat(nrow, clusters, clust1, clust2, ncl);
	    // System.out.println("#clusters left: " + ncl +
            //                 "; cluster label matrix: ");
            // printMatrix(nrow, nrow, clusters, 4, 10);

	    // Update the following: diss, nndiss
            // Strategy:
            // nn[clust2] ceases to exist; similarly nndiss[clust2]
            // ... for all occurrences of clust2
            // nn[clust1] must be updated, as must nndiss[clust1]
            // Only other change is for any nn[i] such that nn[i] =
            // ... clust1 or clust2; this must be updated.

            // First, update diss
           
            cl1 = clust1 - 1;
            cl2 = clust2 - 1;
            for (int i = 0; i < nrow; i++) {
		// The following test is important:
		if ( (i != cl1) && (i != cl2) && (flag[i] == 1) ) {
 		   // Minimum variance criterion
                   diss[cl1][i] = 
                   (mass[cl1]+mass[i])/(mass[cl1]+mass[cl2]+mass[i])
                   *diss[cl1][i]   +
                   (mass[cl2]+mass[i])/(mass[cl1]+mass[cl2]+mass[i])
                   *diss[cl2][i]   -
                   (mass[i])/(mass[cl1]+mass[cl2]+mass[i])
                   *diss[cl1][cl2] ;
		   // For other alternative agglomeration criteria,
		   // e.g. single or complete linkage, see the update
		   // formulas in F. Murtagh, "Multidimensional 
		   // Clustering Algorithms", Physica-Verlag, 1985, p. 68.
		// ( (clcard[cl1]*clcard[cl2])/(clcard[cl1]+clcard[cl2]) )*
                //   (diss[cl1][i] - diss[cl2][i]);
                diss[i][cl1] = diss[cl1][i];
		}
	    }
            clcard[cl1] = clcard[cl1] + clcard[cl2]; 
	    mass[cl1] = mass[cl1] + mass[cl2];
	    // Cluster label clust2 is knocked out
            for (int i = 0; i < nrow; i++) {
                diss[cl2][i] = MAXVAL;
                diss[i][cl2] = diss[cl2][i];
                flag[cl2] = 0; 
                nndiss[cl2] = MAXVAL; 
		mass[cl2] = 0; 
	    }


	// Get nearest neighbors and nearest neighbor dissimilarities
        getNNs(nrow, flag, diss, nn, nndiss);
	// System.out.println("Nearest neighbors of items 1, 2, ... n:");
	// printVect(nn, 4, 10); 
        // System.out.println 
        //   ("Corresponding nearest neighbors dissimilarities:");
        // printVect(nndiss, 4, 10);

	    }  // End of agglomerations loop 
            while (ncl > 1);

	// For convenience, transpose the cluster labels
        int[][] tclusters = new int[nrow][nrow]; 
        for (int i1 = 0; i1 < nrow; i1++) {
            for (int i2 = 0; i2 < nrow; i2++) {
                tclusters[i2][i1] = clusters[i1][i2];
            }
        }

	// Row 1: just one cluster
        // Row 2: two clusters
        // Row 3: three clusters
        // ...
        // Row nrow: nrow clusters
	printMatrix(nrow, nrow, tclusters, 4, 4);

  } // End of try

  catch (IOException e) {
     out.println("error: " + e);
     System.exit(1);
  }
  
    } // End of main  

} // End of class HCL
 





