import java.io.*;
import java.util.*;
import java.text.*;
import Jama.*;

/**
 * Partition - k-means partitioning using exchange method <p>
 * Uses: PCAcorr - Principal Components Analysis on correlations <p>
 * Jama Matrix class package, "JAMA: A Java Matrix Package" is used.<br>
 * See: http://math.nist.gov/javanumerics/jama <p>
 * Example of use: <p>
 * <tt> javac partition.java </tt> <br>
 * <tt> java partition <a href="../iris0.dat">iris0.dat</a> > 
 * <a href="partoutput.txt">partoutput.txt</a></tt> <p>
 * Format of input data set: 
 * <ul>
 * <li> integer row and column dimensions, integer number of clusers,  
 * <li> followed by floating values 
 * <li> which are read row-wise.
 * </ul>
 * Outputs produced:
 * <ul>
 * <li> Echo of input file name, input dimensions, k, sample of data
 * <li> Variable means and standard deviations
 * </ul>
 * Changes in 2003 Nov. 
 * <ul>
 * <li> Corrections in inSort routine (edn. 1 of Numerical Recipes errors!)
 * <li> Clusters prevented from having cardinality of 0 or 1 (zero divide
 * caused).  In this case: abort.
 * </ul>
 * 
 * Version: 2003 Nov. 14 <br>
 * Author: F. Murtagh, f.murtagh@qub.ac.uk
 * @version 2003 Nov. 14 
 * @author F. Murtagh, f.murtagh@qub.ac.uk
 */
public class partition {

public static final double MAXVAL = 1.0e12;
    public static final double R = 0.999;        // See Spaeth, p. 103

public static void main (String argv[])
{
  PrintStream out = System.out;

  try{

     //-------------------------------------------------------------------
     if (argv.length == 0) {
        System.out.println(" Syntax: java partition infile.dat ");
        System.out.println(" Input file format: ");
        System.out.println(" Line 1: integer no. rows, no. cols., no clus.");
        System.out.println(" Successive lines: matrix values, floating");
        System.out.println(" Read in row-wise");
        System.exit (1);
     }
     System.out.println (" K-Means using Spaeth's exchange algorithm,");
     System.out.println (" initialized using first PC projections.");
     String filname = argv[0]; 
     System.out.println (" Input file name: " + filname);

     // Open the matrix file
     FileInputStream is  = new FileInputStream(filname);
     BufferedReader  bis = new BufferedReader(new InputStreamReader(is));
     StreamTokenizer st  = new StreamTokenizer(bis);

     // Row and column sizes, read in first
     st.nextToken();  int n = (int)st.nval;
     st.nextToken();  int m = (int)st.nval;
     st.nextToken();  int k = (int)st.nval;

     System.out.println(" No. of rows, n     = " + n);
     System.out.println(" No. of columns, m  = " + m);
     System.out.println(" No. of clusters, k = " + k);
     
     // Input array, values to be read in successively, float
     double[][] indat = new double[n][m]; 
     double inval;

     // New read in input array values, successively
     for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
             st.nextToken();
             inval = (double)st.nval;
             indat[i][j] = inval;
        }
     }  

     // Data preprocessing - standardization and determining correlations 
     double[][] indatstd = Standardize(n, m, indat);
     // use Jama matrix class  
     Matrix X = new Matrix(indatstd);
     // Sums of squares and cross-products matrix
     Matrix Xprime = X.transpose();
     Matrix SSCP   = Xprime.times(X);
     // Eigen decomposition
     EigenvalueDecomposition evaldec = SSCP.eig();
     Matrix evecs = evaldec.getV();

     // evecs contains the cols. ordered right to left
     // Evecs is the more natural order with cols. ordered left to right
     // So to repeat: leftmost col. of Evecs is assoc with largest Evals
     // Evecs ordered from left to right
     // reverse order of Matrix evecs into Matrix Evecs
     double[][] tempold = evecs.getArray();
     double[][] tempnew = new double[m][m]; 
     for (int j1 = 0; j1 < m; j1++) {
         for (int j2 = 0; j2 < m; j2++) {
             tempnew[j1][j2] = tempold[j1][m - j2 - 1];
	 }
     }
     Matrix Evecs = new Matrix(tempnew);
     // Evecs.print(10,4);

     // Projections - for rows
     // Row projections in new space, X U  Dims: (n x m) x (m x m)
     Matrix rowproj = X.times(Evecs);
     // rowproj.print(10,4);

     //-------------------------------------------------------------------
     //  Now k-means partitioning using exchange method

     double[] tempvec = new double[n];
     double[] rproj = new double[n]; 
     tempold = rowproj.getArray();
     for (int i = 0; i < n; i++) {
         tempvec[i] = tempold[i][0];
	 rproj[i] = tempold[i][0];
     }
     inSort(tempvec);  // Sort in place, modifying values.
     // printVect(tempvec, 4, 10);

     // First choose breakpoints
     int[] breakindexes = new int[k+1];
     double[] breakpoints = new double[k+1];
     breakpoints[0] = tempvec[0] - 0.1; //Offset so that strict ">" used later 
     breakpoints[k] = tempvec[n-1];
     for (int i = 1; i < k; i++) {
	 breakindexes[i] = i*n/k;
	 breakpoints[i] = tempvec[breakindexes[i]];
     }
     // printVect(breakpoints, 4, 10);

     double[][] gpmeans = new double[k][m];
     double[] cardinality = new double[k]; 
     int[] assignment = new int[n];
     for (int i = 0; i < k-1; i++) {
         cardinality[i] = 0.0;
	 for (int j = 0; j < m-1; j++) {
	     gpmeans[i][j] = 0.0;
	 }
     }
 
     for (int icl = 0; icl < k; icl++) {
	 // lo, hi (resp.) are breakpoints[i], breakpoints[i+1] 
	 // cluster = icl, with values 0, 1, 2, ... k-1
	 for (int i = 0; i < n; i++) {
	     // Not terribly efficient - fix later
 	     if ( (rproj[i] > breakpoints[icl]) && 
                (rproj[i] <= breakpoints[icl+1]) ) {
		 assignment[i] = icl;
		 cardinality[icl] += 1.0;
		 for (int j = 0; j < m; j++) {
		     gpmeans[icl][j] += indat[i][j];
		 }
	     }
	 }
	 for (int j=0; j<m; j++) gpmeans[icl][j] /= cardinality[icl];
     }
     for (int i = 0; i < k; i++) {
         System.out.println
         (" Cluster number, cardinality = " + i + " " + cardinality[i]);
         if (cardinality[i] <= 1) {
            System.out.println
	    (" A cluster has cardinality <= 1; stopping.");
	    // Card = 1 in expressions below, with card - 1 will cause 
	    // zero divide.  So use this as a condition to halt.     
	    System.exit(1);
	 }   
     }

     System.out.println("Initial cluster means:");
     printMatrix(k, m, gpmeans, 4, 10);
     // printVect(assignment, 4, 10);

     //-------------------------------------------------------------------
     // Having an initial partition, with group mean vectors, 
     // cardinalities, and asssignments, we now proceed...


     double compactness = 0.0;
     for (int i = 0; i < n; i++) {
	 for (int j = 0; j < m; j++) {
	     compactness += Math.pow( indat[i][j] - 
				      gpmeans[assignment[i]][j], 2.0); 
	 }
     }
     System.out.println("Compactness = " + compactness);

     int epochmax = 15;
     int epoch = 0; 
     double oldcompactness = MAXVAL;
     double prevcontrib = 0.0; 
     double newcontrib = 0.0;  
     double bestnewcontrib; 
     int bestclus;
     int oldclus = 0;
     double eps = 0.001; 
     double dissim = 0.0; 
     int nochange = 1; 

     while ( (nochange > 0) && (compactness < (oldcompactness-eps) ) && 
         (epoch <= epochmax) ) {
         nochange = 0; 
	 epoch += 1;      
         System.out.println("EPOCH NUMBER IS = " + epoch);

	 for (int i = 0; i < n; i++) {

	     bestnewcontrib = MAXVAL; 
	     bestclus = 0; 

	     for (int c = 0; c < k; c++) {

                 dissim = 0.0;
                 for (int j = 0; j < m; j++) 
                     dissim += Math.pow(indat[i][j] - gpmeans[c][j], 2.0);

		 if (assignment[i] == c) {  // Same cluster
		     prevcontrib = 
                           (cardinality[c]/(cardinality[c]-1.0))*dissim;
		 }
		 else  { // Potentially new cluster
		     newcontrib = 
                           (cardinality[c]/(cardinality[c]+1.0))*dissim; 
		     if (newcontrib < bestnewcontrib) {
		         bestnewcontrib = newcontrib; 
		         bestclus = c;
		     } 
		 }
	     } // End of c loop 
	     

	     // Is bestnewcontrib better than prevcontrib?	      
	     if (bestnewcontrib < R*prevcontrib) {
		 // A change is in store
                 nochange += 1;

		 oldcompactness = compactness;
		 compactness = compactness + bestnewcontrib - prevcontrib; 
                 //System.out.println("Old compactness = " + oldcompactness);
                 System.out.println("New compactness = " + compactness);


		 oldclus = assignment[i];
		 // bestclus is the new cluster
		 for (int j = 0; j < m; j++) {
		   gpmeans[oldclus][j] = ( 1.0/(cardinality[oldclus]-1.0) )
		       * ( cardinality[oldclus] * gpmeans[oldclus][j]
			   - indat[i][j] );
		   gpmeans[bestclus][j] = ( 1.0/(cardinality[bestclus]+1.0) )
		       * ( cardinality[bestclus] * gpmeans[bestclus][j]
			   + indat[i][j] );
		 }

                 cardinality[oldclus] -= 1.0;
                 if (cardinality[oldclus] <= 1.0) {
		     System.out.println
                       ("Cluster has become too small: stopping."); 
		     System.out.println
                       ("Try a smaller number of clusters instead."); 
		     System.exit(1);
		 }
                 cardinality[bestclus] += 1.0;
                 assignment[i] = bestclus;

		 // Check on progress
		 // printVect(cardinality, 4, 10);

	     }
	 

	 } // End of i loop 
     } // End of while loop

     System.out.println();
     System.out.println("Cluster centers:");
     printMatrix(k, m, gpmeans, 4, 10);
     System.out.println("Final cluster cardinalities:");
     printVect(cardinality, 4, 10);
     System.out.println("Assignments:");
     printVect(assignment, 1, 2);
  

     //-------------------------------------------------------------------
     // That's it.

  } // End of try loop

  catch (IOException e) {
     out.println("error: " + e);
     System.exit(1);
  }

} // end of main

    //-------------------------------------------------------------------  
    // Little method for helping in output formating  
    public static String getSpaces(int n) {

    StringBuffer sb = new StringBuffer(n);
    for (int i = 0; i < n; i++) sb.append(' ');
    return sb.toString();
    } // getSpaces
   
    //-------------------------------------------------------------------
    /**
     * Method for printing a matrix  <br>
     * Based on ER Harold, "Java I/O", O'Reilly, around p. 473.
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values 
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
              for (int j=0; j<n2; j++)
		  // Following unfortunately doesn't format at all
		  //                  System.out.print(m[i][j] + "  ");
                  {
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
    } // printMatrix

    //-------------------------------------------------------------------
    /**
     * Method printVect for printing a vector <br>
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
      for (int i=0; i<len; i++)
	  {
	      // Following would be nice, but doesn't format adequately
	      //                  System.out.print(m[i] + "  ");
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

    //-------------------------------------------------------------------
    /**
     * Method printVect for printing a vector <br>
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
      for (int i=0; i<len; i++)
	  {
	      // Following would be nice, but doesn't format adequately
	      //                  System.out.print(m[i] + "  ");
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

    //-------------------------------------------------------------------
    /**
     * Method for standardizing the input data <p>
     * Note the formalas used (since these vary between implementations):<br>
     * reduction: (vect - meanvect)/sqrt(nrow)*colstdev <br>
     * colstdev:  sum_cols ((vect - meanvect)^2/nrow) <br>
     * if colstdev is close to 0, then set it to 1.
     * @param nrow number of rows in input matrix
     * @param ncol number of columns in input matrix
     * @param A input matrix values
     */
    public static double[][] Standardize(int nrow, int ncol, double[][] A)
    {
    double[] colmeans = new double[ncol];
    double[] colstdevs = new double[ncol];
    // Adat will contain the standardized data and will be returned
    double[][] Adat = new double[nrow][ncol];
    double[] tempcol = new double[nrow];
    double tot; 

    // Determine means and standard deviations of variables/columns
    for (int j=0; j<ncol; j++)
        {
         tot = 0.0;
         for (int i=0; i<nrow; i++)
             {
                tempcol[i] = A[i][j];
                tot += tempcol[i];
             }

	 // For this col, det mean
         colmeans[j] = tot/(double)nrow;
         for (int i=0; i<nrow; i++) {
                colstdevs[j] += Math.pow(tempcol[i]-colmeans[j], 2.0);
             }
             colstdevs[j] = Math.sqrt(colstdevs[j]/((double)nrow));
             if (colstdevs[j] < 0.0001) { colstdevs[j] = 1.0; }
	}

        System.out.println(" Variable means:");
        printVect(colmeans, 4, 8);
        System.out.println(" Variable standard deviations:");
        printVect(colstdevs, 4, 8);

	// Now ceter to zero mean, and reduce to unit standard deviation
    for (int j=0; j<ncol; j++)
        {
         for (int i=0; i<nrow; i++)
             {
               Adat[i][j] = (A[i][j] - colmeans[j])/
                 (Math.sqrt((double)nrow)*colstdevs[j]);
             }
        }
    return Adat;
    } // Standardize
    //-------------------------------------------------------------------


    //-------------------------------------------------------------------
    /**
     * Method to sort a double vector, by inefficient straight insertion
     * See section 8.1, p. 243, of Numerical Recipes in C.
     * Two corrections, FM, 2003/11/14
     * @param invect input data to be sorted, sorted in place
     */
   private static void inSort (double[] invect) {
        double a;
        int i; 

        //for (int j = 2; j < invect.length; j++) {   CORRECTED! 
        for (int j = 1; j < invect.length; j++) {
            a = invect[j]; 
            i = j - 1;
            //while (i > 0 && invect[i] > a) {        CORRECTED!
            while (i >= 0 && invect[i] > a) {
                invect[i+1] = invect[i]; 
                i--; 
	    }
            invect[i+1] = a; 
	}
	// Return type void
   }  // inSort



}  // PCAcorr


