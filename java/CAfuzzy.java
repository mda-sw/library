import java.io.*;
import java.util.*;
import java.text.*;
import Jama.*;

/**
 * CAfuzzy - Correspondence Analysis on Fuzzy Input <p>
 * Jama Matrix class package, "JAMA: A Java Matrix Package" is used.<br>
 * See: http://math.nist.gov/javanumerics/jama <p>
 * Example of use: <p>
 * <tt> javac CAfuzzy.java </tt> <br>
 * <tt> java CAfuzzy <a href="../iris.dat">iris.dat</a> > 
 * <a href="pcaoutput.txt">caoutput.txt</a></tt> <p>
 * Format of input data set: 
 * <ul>
 * <li> integer row and column dimensions, 
 * <li> followed by floating values 
 * <li> which are read row-wise.
 * </ul>
 * Outputs produced: 
 * <ul>
 * <li> Echo of input file name, input dimensions, sample of data
 * <li> Matrix to be diagonalized
 * <li> Eigenvectors
 * <li> Eigenvalues and as cumulative percentages
 * <li> Row projections in new factor space
 * <li> Column projections in new factor space
 * <li> Row contributions to factors
 * <li> Column contributions to factors
 * </ul>
 * Version: 2002 Aug. 16 <br>
 * Author: F. Murtagh, f.murtagh@qub.ac.uk
 * @version 2002 Aug 16
 * @author F. Murtagh, f.murtagh@qub.ac.uk
 */
public class CAfuzzy{

public static final double EPS = 1.0e-8;

public static void main (String argv[])
{
  PrintStream out = System.out;

  try{

     //-------------------------------------------------------------------
     if (argv.length == 0) {
        System.out.println(" Syntax: java CAfuzzy infile.dat ");
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
     st.nextToken();  int n = (int)st.nval;
     st.nextToken();  int m = (int)st.nval;

     System.out.println(" No. of rows, n = " + n);
     System.out.println(" No. of cols, m = " + m);
     
     // Input array, values to be read in successively, float
     double[][] indat = new double[n][m]; 
     double inval;

     // New read in input array values, successively
     System.out.println
        (" Input data sample follows as a check, first 4 values.");
     for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
             st.nextToken();
             inval = (double)st.nval;
             indat[i][j] = inval;
             if (i < 2 && j < 2) {
                System.out.println(" value = " + inval);
             }
        }
     }  
     System.out.println();

     // Some output formating parameters
     int collo = 2;   // start col for output printing
     int colhi = 5;   // last col for output printing
     int nplac = 10;  // no. of places in output float values
     int ndec = 4;    // no. of dec. places in output

     //-------------------------------------------------------------------
     int mnew;
     // Data preprocessing - fuzzify 
     double[][] indatstd = Fuzzify(n, m, indat);
     mnew = 2*m;         // fuzzification has expanded col dim by 2
     // Data to be analyzed is now indatstd, with dims n x 2m, or n x mnew

     double rowsums[] = new double[n];
     double colsums[] = new double[mnew];
     double total = 0.0; 

     // Row sums and overall total 
     for (int i = 0; i < n; i++) {
         rowsums[i] = 0.0;
         for (int j = 0; j < mnew; j++) {
             rowsums[i] += indatstd[i][j];
             total += indatstd[i][j];
         }
     }         

     // Col sums 
     for (int j = 0; j < mnew; j++) {
         colsums[j] = 0.0;
         for (int i = 0; i <n; i++) {
             colsums[j] += indatstd[i][j];
	 }
     }

     // Finalize nomalization to provide masses by dividing by total
     for (int i = 0; i < n; i++) rowsums[i] /= total;
     for (int j = 0; j < mnew; j++) colsums[j] /= total;
     for (int i = 0; i < n; i++) {
	 for (int j = 0; j < mnew; j++) {
             indatstd[i][j] /= total;
	 }
     }

     // Form matrix to be analyzed (diagonalized)
     double CP[][] = new double[mnew][mnew]; //cross-products, e.g. Burt table
     for (int j1 = 0; j1 < mnew; j1++) {
	 for (int j2 = 0; j2 < mnew; j2++) {
             CP[j1][j2] = 0.0;
             for (int i = 0; i < n; i++) {
                 CP[j1][j2] += indatstd[i][j1] *
                               indatstd[i][j2] /
                               (rowsums[i] * 
                               Math.sqrt(colsums[j1]*colsums[j2])); 
	     }
	 }
     }

     // use Jama matrix class  
     Matrix cp = new Matrix(CP);

     // Print out cp matrix 
     System.out.println 
            (" Matrix to be diagonalized");
     cp.print(6,2);

     //-------------------------------------------------------------------
     // Eigen decomposition
     EigenvalueDecomposition evaldec = cp.eig();
     Matrix evecs = evaldec.getV();
     double[] evals = evaldec.getRealEigenvalues();

     // print out eigenvectors
     System.out.println 
        (" Eigenvectors (leftmost col <--> largest eval, first always = 1):");
     // evecs contains the cols. ordered right to left
     // Evecs is the more natural order with cols. ordered left to right
     // So to repeat: leftmost col. of Evecs is assoc with largest Evals
     // Evals and Evecs ordered from left to right

     double tot = 0.0; 
     for (int j = 0; j < evals.length; j++)  {
         tot += evals[j]; 
     }

     // reverse order of evals into Evals
     double[] Evals = new double[mnew];
     for (int j = 0; j < mnew; j++) {
         Evals[j] = evals[mnew - j - 1];
     }
     // reverse order of Matrix evecs into Matrix Evecs
     double[][] tempold = evecs.getArray();
     double[][] Evex = new double[mnew][mnew]; 
     for (int j1 = 0; j1 < mnew; j1++) {
         for (int j2 = 0; j2 < mnew; j2++) {
             Evex[j1][j2] = tempold[j1][mnew - j2 - 1] / 
                                    Math.sqrt(colsums[j1]);
	 }
     }
     Matrix Evecs = new Matrix(Evex);
     Evecs.print(10,4);
     // So as a JAMA Matrix, we have Evecs, and 
     // as a multidim. array of dims mnew x mnew, we have Evex

     System.out.println();
     System.out.println 
        (" E-vals. (first always = 1) and as cumul. perc. (first excl.):"); 
     double runningtotal = 0.0;
     double[] percentevals = new double[mnew];
     // low index in following = 1 to exclude first trivial eval.
     percentevals[0] = 0.0; 
     for (int j = 1; j < Evals.length; j++) {
         percentevals[j] = runningtotal + 100.0*Evals[j]/(tot-1.0);
         runningtotal    = percentevals[j];
     }
     printVect(Evals, 4, 10);
     printVect(percentevals, 4, 10);

     //-------------------------------------------------------------------
     // Projections - row, and col
     // Row projections in new space, X U  Dims: (n x m) x (m x m)
     System.out.println();
     System.out.println
       (" Row projections in new factor space. Start, end columns: " +
          collo + " " + colhi);
     double rowproj[][] = new double[n][mnew]; 
     for (int i = 0; i < n; i++) {
         for (int j1 = 0; j1 < mnew; j1++) {
             rowproj[i][j1] = 0.0; 
             for (int j2 = 0; j2 < mnew; j2++) {
                rowproj[i][j1] += indatstd[i][j2]*Evex[j2][j1];
	     }
             if (rowsums[i] >= EPS) rowproj[i][j1] /= rowsums[i];
             if (rowsums[i] < EPS) rowproj[i][j1] = 0.0;
	 }
     }
     // Print from cols collo to colhi; 
     // precision: nplac chars in total, ndec dec places
     printMatrix(n, mnew, rowproj, collo, colhi, ndec, nplac);

     System.out.println ();
     System.out.println 
          (" Column projections in new factor space. Start, end columns: " +
          collo + " " + colhi);
     double colproj[][] = new double[mnew][mnew];
     for (int j1 = 0; j1 < mnew; j1++) {
         for (int j2 = 0; j2 < mnew; j2++) {
             colproj[j1][j2] = 0.0;
             for (int j3 = 0; j3 < mnew; j3++) {
                 colproj[j1][j2] += CP[j1][j3] * Evex[j3][j2] *
                                    Math.sqrt(colsums[j3]);
	     }
             if (colsums[j1] >= EPS && Evals[j2] >= EPS) 
                colproj[j1][j2] /= Math.sqrt(Evals[j2]*colsums[j1]);
             if (colsums[j1] < EPS && Evals[j2] < EPS) 
                colproj[j1][j2] = 0.0; 
	 }
     }
     // Print from cols collo to colhi; 
     // precision: nplac chars in total, ndec dec places
     printMatrix(mnew, mnew, colproj, collo, colhi, ndec, nplac);

     System.out.println ();
     System.out.println 
         (" Row contributions to factors. Start, end columns: " +
          collo + " " + colhi);
     double rowcntr[][] = new double[n][mnew]; 
     double rowconColsum; 
     for (int j = 0; j < mnew; j++) {
         rowconColsum = 0.0;
         for (int i = 0; i < n; i++) {
             rowcntr[i][j] = rowsums[i]*Math.pow(rowproj[i][j], 2.0);
             rowconColsum += rowcntr[i][j];
	 }
         // Normalize so that sum of contributions for a factor equals 1
         for (int i = 0; i < n; i++) {
             if (rowconColsum > EPS) rowcntr[i][j] /= rowconColsum;
             if (rowconColsum <= EPS) rowcntr[i][j] = 0.0; 
	 }
     }
     // Print from cols collo to colhi; 
     // precision: nplac chars in total, ndec dec places
     printMatrix(n, mnew, rowcntr, collo, colhi, ndec, nplac); 

     System.out.println ();
     System.out.println 
          (" Column contributions to factors. Start, end columns: " +
          collo + " " + colhi);
     double colcntr[][] = new double[mnew][mnew]; 
     double colconColsum; 
     for (int j1 = 0; j1 < mnew; j1++) {
         colconColsum = 0.0;
         for (int j2 = 0; j2 < mnew; j2++) {
             colcntr[j2][j1] = colsums[j2]*Math.pow(colproj[j2][j1], 2.0);  
             colconColsum +=  colcntr[j2][j1];
         }
	 // Normalize so that sum of contributions for a factor sum to 1
         for (int j2 = 0; j2 < mnew; j2++) {
             if (colconColsum > EPS) colcntr[j2][j1] /= colconColsum;
             if (colconColsum <= EPS) colcntr[j2][j1] = 0.0;
	 }
     }
     // Print from cols collo to colhi; 
     // precision: nplac chars in total, ndec dec places
     printMatrix(mnew, mnew, colcntr, collo, colhi, ndec, nplac); 



     //-------------------------------------------------------------------
     // That's it.

  }

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
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values
     * @param collo start column for printing 
     * @param colhi end column for printing 
     * @param d display precision, number of decimal places
     * @param w display precision, total width of floating value
     */
    public static void printMatrix(int n1, int n2, double[][] m, 
           int collo, int colhi, int d, int w)
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
              for (int j = (collo-1); j < colhi; j++)
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
     * Method for fuzzyifying the input data <p>
     * @param nrow number of rows in input matrix
     * @param ncol number of columns in input matrix
     * @param A input matrix values
     * @param ncolnew new number of columns in transformed data
     */
    public static double[][] Fuzzify(int nrow, int ncol, double[][] A)
    {
    // Adat will contain the fuzzified data and will be returned
    double[][] Adat = new double[nrow][2*ncol];

     double[] coltemp = new double[nrow]; 
     double perc25, perc50, perc75;       // percentiles
     int medloc, p25loc, p75loc;          // percentile offsets
     medloc = nrow/2; 
     p25loc = nrow/4;
     p75loc = (3*nrow)/4;

     for (int j = 0; j < ncol; j++) {
	 // for each col
         for (int i = 0; i < nrow; i++) coltemp[i] = A[i][j];
         inSort(coltemp);   // sort in place, modifying values
         perc50 = coltemp[medloc];
         perc25 = coltemp[p25loc];
         perc75 = coltemp[p75loc];
         // System.out.println(perc50 + " " + 
         //                    perc25 + " " +
         //                    perc75);

         for (int i = 0; i < nrow; i++) {
             if (A[i][j] > perc75) {
                 Adat[i][2*j] = 0;
                 Adat[i][2*j+1] = 1; 
	     }
             if (A[i][j] < perc25) {
                 Adat[i][2*j] = 1;
                 Adat[i][2*j+1] = 0;
	     }
             if (A[i][j] <= perc75 && A[i][j] >= perc25) {
                 Adat[i][2*j] = (A[i][j] - perc25)/(perc75 - perc25);
                 Adat[i][2*j+1] = 1.0 - Adat[i][2*j];
	     }		    
	 }
     }

    return Adat;
    } // Fuzzify
    //-------------------------------------------------------------------


    /**
     * Method to sort a double vector, by inefficient straight insertion
     * See section 8.1, p. 243, of Numerical Recipes in C.
     * @param invect input data to be sorted, sorted in place
     */
   private static void inSort (double[] invect) {
        double a;
        int i; 

        for (int j = 2; j < invect.length; j++) {
            a = invect[j]; 
            i = j - 1;
            while (i > 0 && invect[i] > a) {
                invect[i+1] = invect[i]; 
                i--; 
	    }
            invect[i+1] = a; 
	}
	// Return type void
   }  // inSort



}  // CAfuzzy


