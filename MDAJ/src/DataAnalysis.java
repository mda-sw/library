
import java.io.*;
import java.util.*;
import java.text.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.table.*;
// import com.macfaq.io.*;

/**
 * DataAnalysis
 <P>
 DataAnalysisJava is a first version of a data analysis package for 
 Java.  It provides functionality for correspondence analysis, 
 plots of factor projections, hierarchical classification, and 
 mutual interpretations of factors and clusters.  It is comprised
 of the following classes: 
 <OL>
 <LI> DataAnalyis, which parses input data sets, manages the 
 eigen-processing and related correspondence analysis,
 followed by the plots, hierarchical clustering, and interpretations.
 <LI> CorrAnal, carrying out most of the correspondence analysis.
 <LI> HierClass, carrying out the hierarchical clustering.
 <LI> Interpret, carrying out the factor vs. cluster interpretations.
 <LI> JAMA classes: CholeskyDecomposition, EigenvalueDecomposition,
 LUDecomposition, Matrix, Maths, QRDecomposition, SingularValueDecomposition.
 (JAMA: A Java Matrix Package, math.nist.gov/javanumerics/jama)
 <LI> Imatrix, Header, ManageOutput, PrintUtilities, ResultCA, 
 ResultHC, ResultSVD, ResultInterpret
 </OL>
 Usage:
 <P>
 In command prompt window on a Windows NT/2000/XP system, or a window
 on a Unix or Linux system, we will assume that you will either
 <OL>
 <LI> compile all source files by doing <i>javac -classpath . *.java</i>
 (assuming the Java compiler, <i>javac</i>, is installed on your system),
 and/or followed by
 <LI> run: <i>java -classpath . DataAnalysis</i>
 </OL>
 If your <i>classpath</i> is already configured it can be omitted in the 
 above.  Note: this code has not yet been tested on Macintosh OS X, but 
 should run, modulo the uncommenting of the import near the start of the
 source files.
 <P>
 @author Fionn Murtagh, fmurtagh@acm.org
 @version 0.0, 2003-March-15
 */

public class DataAnalysis extends JFrame{

    public static final double EPS = 1.0e-8;
    public static final double SMALL = -1.0e10;
    public static final double MAXVAL = 1.0e12;

    public static void main (String argv[])
    {
        PrintStream out = System.out;

	try {

	    /* Step 1 -------------------------------------------------
	     * Manage window which will provide input file 
	     * -------------------------------------------------------- */

           // Create a file chooser that opens up as an Open dialog
           // Open the input file
	    String filname = "test.dat";           // Initialize arbitrarily
	    String pathname = " ";                 // Initialize

           JFrame parent = new JFrame();
           JFileChooser chooser = new JFileChooser();
           chooser.setSize(250,300);
           int option = chooser.showOpenDialog(parent);
           if (option == JFileChooser.APPROVE_OPTION) {
	       File f = chooser.getSelectedFile();
	       filname = f.getName();
	       pathname = f.getPath(); 
	       System.out.println("You chose file " + filname); 
	       System.out.println("With full path " + pathname); 
           }
           else
	   { 
	       System.out.println("You cancelled - no input file - stopping.");
	       System.exit(1);
	   }
           String filename = filname.toLowerCase();
           if (filename.endsWith(".dat")) System.out.println("File okay!");

           FileInputStream is  = new FileInputStream(pathname);
           BufferedReader  bis = new BufferedReader(new InputStreamReader(is));
           StreamTokenizer st  = new StreamTokenizer(bis);

     /* Step 2 -------------------------------------------------
      * Create and manage a number of windows, which will contain 
      * various phases of the output.
      * -------------------------------------------------------- */

     String outTitle = " ";  // To be used in window header.
     int fontSize = 0;       // Usually: 14 pt; for plots: 9 pt.
     boolean close = false;  // Whether or not closing windows exits all.

     // Output text window for correspondence analysis results
     JTextArea outext1 = new JTextArea();
     // outTitle = "R\u00e9sultats de l'Analyse des Correspondances"; 
     outTitle = "Results of Correspondance Analysis"; 
     fontSize = 14;
     createWindow(outext1, outTitle, fontSize, close); 

     // Output text window for correspondence analysis plot 1 vs. 2 result
     JTextArea outext1f12 = new JTextArea();
     // outTitle = "Plan des facteurs 1 et 2"; 
     outTitle = "Plane of factors 1 and 2"; 
     fontSize = 9;
     createWindow(outext1f12, outTitle, fontSize, close); 

     // Output text window for correspondence analysis plot 1 vs. 3 result
     JTextArea outext1f13 = new JTextArea();
     // outTitle = "Plan des facteurs 1 et 3"; 
     outTitle = "Plane of factors 1 and 3"; 
     fontSize = 9;
     createWindow(outext1f13, outTitle, fontSize, close); 

     // Output text window for correspondence analysis plot 2 vs. 3 result
     JTextArea outext1f23 = new JTextArea();
     // outTitle = "Plan des facteurs 2 et 3"; 
     outTitle = "Plane of factors 2 and 3"; 
     fontSize = 9;
     createWindow(outext1f23, outTitle, fontSize, close); 

     // Output text window for hierarchical clustering results
     JTextArea outext2 = new JTextArea();
     // outTitle = "R\u00e9sultats de la Classification Hi\u00e9rarchique"; 
     outTitle = "Results of Hierarchical Classification"; 
     fontSize = 14;
     createWindow(outext2, outTitle, fontSize, close); 

     // Output text window for factors vs. clusters interpretation results
     JTextArea outext3 = new JTextArea();
     // outTitle = "R\u00e9sultats de l'Analyse FACOR";
     outTitle = "Results of FACOR Interpretational Analysis";
     fontSize = 14;
     createWindow(outext3, outTitle, fontSize, close); 

     // Output text window for variables vs. clusters interpretation results
     JTextArea outext4 = new JTextArea();
     // outTitle = "R\u00e9sultats de l'Analyse VACOR";
     outTitle = "Results of VACOR Interpretational Analysis";
     fontSize = 14;
     createWindow(outext4, outTitle, fontSize, close); 

     // Output text window for logging phases of processing
     JTextArea outext5 = new JTextArea();
     // outTitle = "Fen\u00eatre de suivie: l'analyse du fichier " + filename;
     outTitle = "Control Window: Analysis of File " + filename;
     fontSize = 14;
     close = true; 
     createWindow(outext5, outTitle, fontSize, close); 


     /* Step 3 --------------------------------------------------
      * Some output into a control window, for information.
      * --------------------------------------------------------- */

     Date date = new Date(); 
     outext5.append(date.toString() + "\n");
     // outext5.append("Fichier d'entr\u00e9e : " + filename + "\n");
     outext5.append("Input file : " + filename + "\n");
     // outext5.append
     //  ("Les 5 fen\u00eatres de sortie (initialement surimpos\u00e9es): \n");
     // outext5.append
     //  ("1. R\u00e9sultat de l'Analyse des Correspondances. \n");
     // outext5.append
     //  ("2. R\u00e9sultat de la Classification Hi\u00e9rarchique. \n");
     // outext5.append("3. R\u00e9sultat de l'Analyse FACOR \n");
     // outext5.append("4. R\u00e9sultat de l'Analyse VACOR \n");
     // outext5.append("5. Fen\u00eatre de suive et d'information. \n");
     // outext5.append("Fermer celle-ci pour terminer le programme. \n");
     // outext5.append
     //   ("Col et copier \u00e0 partir de chaque fen\u00eatre. \n \n");
     // outext5.append("Lecture des donn\u00e9es ... \n"); 
     outext5.append
        ("The 5 output windows (initially superimposed): \n");
     outext5.append("1. Results of Correspondence Analysis. \n");
     outext5.append("2. Results of Hierarchical Classification. \n");
     outext5.append("3. Results of FACOR Analysis. \n");
     outext5.append("4. Results of VACOR Analysis. \n");
     outext5.append("5. Control Window. (Close to Terminate Program.) \n");
     outext5.append("Cut and paste from any window. \n \n");
     outext5.append("Inputing data ... \n"); 

     /* Step 4a ----------------------------------------------------
      * Read - parse and store - input data file.
      * First - get data file header information, which consists of:
      * title, nrows, ncols, whether or not suppl. elements, 
      * no. factors, no. clusters for analysis purposes, followed 
      * by all column labels.
      * --------------------------------------------------------- */

     String title = null;   // title 
     int norig = 0;   // (Original) number of rows
     int morig = 0;   // (Original) number of columns
     int nfactors = 0;      // Number of factors to be output by default
     int nclusters = 4;     // Number of clusters to be analyzed 
     String[] collabsOrig = null;  // Column labels
     boolean supR = false;  // Supplementary rows present?
     boolean supC = false;  // Supplementary columns present?
     boolean supRC = false; // Supplementary rows and columns present?
                            // Warning: both not allowed simultaneously!
     boolean supNO = false; // Default, to allow for any misspelled string
     int n = 0, m = 0;      // If we have supplementary rows or columns, 
                            // then n or m will equal norig+1 or morig+1
     String suppIndicator = null;  // 'SUPR' or 'SUPC' or (default) 'SUPNO'

     // Object dataparams of class Header will store the header information.
     Header dataparams = new Header(norig, morig, nfactors, nclusters,
     title, suppIndicator, collabsOrig); // Initializing...

     getHeader(dataparams, bis, st, outext5);

     title         = dataparams.getTitle(); 
     norig         = dataparams.getNrows();
     morig         = dataparams.getNcols();
     nfactors      = dataparams.getNfactors();
     nclusters     = dataparams.getNclusters();
     collabsOrig   = dataparams.getColumnIdentifiers();
     suppIndicator = dataparams.getSuppIndicator();

     // Note: we can have 
     // - no indicator row, or               SUPNO or SUPC
     // - no indicator column, or            SUPNO or SUPR
     // - an indicator row, or               SUPR
     // - an indicator column.               SUPC
     // We do not allow currently: both an ind. row and an ind. column. 

     // Temporarily, in this program, to read in data, we expand the 
     // numbers of rows and/or of columns.  SUPR implies an extra column
     // of indicator 1s and 0s.  SUPC implies an extra row of indicator 
     // 1s and 0s.  SUPNO implies no change in the number of rows or cols.
 
     n = norig; 
     m = morig; 

     // To avoid any problems with upper/lower case, use only lower case.
     String suppInd = suppIndicator.toLowerCase();   
     if (suppInd.equals("supr")) { 
        supR = true;
	m++; 
     }
     if (suppInd.equals("supc")) {
        supC = true;
	n++; 
     }
     if (suppInd.equals("supno")) supNO = true;     // No change in n, m.

     // A quick test on validity of input.
     if (supR == false && supC == false && supRC == false && supNO == false) {
	 System.out.println
	     ("We were expecting indicator of supplementary elements");
	 System.out.println("following numbers of rows, columns.");
	 System.out.println
             ("Use one only of: SUPR, SUPC, (default) SUPNO.");
	 System.exit(1);
     }

     /* Step 4b ----------------------------------------------------
      * Read - parse and store - input data file.
      * Secondly, rest of essential input data storage needed - 
      * row labels, data array.
      * --------------------------------------------------------- */

     // Input array, values to be read in successively, float
     double[][] indat = new double[n][m]; 
     double inval;
     // Row labels
     String[] rowlabsOrig = new String[n];

     // Read in input data array row-wise, each row preceded by a label 
     // (character identifier), knowing number of rows and of columns.

     getArray(n, m, st, indat, rowlabsOrig); 
     //    public static void getArray(int n, int m, StreamTokenizer st,
     //    double[][] indat, String[] rowlabs) { 

     outext5.append
        ("Input file parsed ... \n");
     // ("Fichier d'entr\u00e9e lu, syntatiquement analys\u00e9 ... \n");

     outext5.append("\n");
     // outext5.append("Donn\u00e9es d'entr\u00e9ee: \n");
     outext5.append("Input data: \n");
     // Some definitions for handling output formating
     NumberFormat myFormat = NumberFormat.getNumberInstance();
     FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
     // Following suppresses e.g. comma in 1,000 = 1000 for English locale.
     myFormat.setGroupingUsed(false);
     // In the next few lines we say how we want numbers to appear.
     // max integer digits = max no. of digits before dec. point
     // max fraction digits = max no. of digits following dec. point
     // min fraction digits = min no. of digits following dec. point
     myFormat.setMaximumIntegerDigits(4);
     // We will set min and max nos. of digits following dec. point to 0
     myFormat.setMaximumFractionDigits(2);
     myFormat.setMinimumFractionDigits(2);
     for (int i = 0; i < n; i++) {
	 for (int j = 0; j < m; j++) {
             String valString = myFormat.format(
                       indat[i][j], new StringBuffer(), fp).toString();
	     // With a max field width of 6, we pack spaces before val.
             valString = getSpaces(8 - fp.getEndIndex()) + valString;
             outext5.append(valString);
	 }
	 outext5.append("\n");
     }
     outext5.append("\n");

     // outext5.append("D\u00e9roulement de l'analyse ... \n");
     outext5.append("Proceeding with analysis ... \n");


     /* Step 5 ------------------------------------------------------
      * Having out input data, organize our principal data array,
      * and our supplementary rows (if present), or our supplementary
      * columns (if present).
      * ------------------------------------------------------------- */ 

     // indat contains our data, and is of dims n x m 
     // n = norig + 1 if supC || supRC
     // m = morig + 1 if supR || supRC

     int suppRowCount = 0; 
     int mainRowCount = 0; 
     int suppColCount = 0; 
     int mainColCount = 0; 

     // The case of no supplementary rows nor supplementary columns
     if (supNO) {
	// In this case n = norig, m = morig -- check this
	if ( (n != norig) && (m != morig) ) {
	    System.out.println
		("There is a problem: we should have n=norig and m=morig.");
	    System.out.println("Aborting.");
	    System.exit(1);
	}
	// Following just for the record
	mainRowCount = n;
	mainColCount = m; 
	suppRowCount = 0; 
	suppColCount = 0; 
        outext5.append
               ("Supplementary rows: none present.\n"); 
	//     ("Lignes suppl\u00e9mentaires: non pr\u00e9sentes.\n"); 
        outext5.append
               ("Supplementary columns: none present.\n"); 
        //     ("Colonnes suppl\u00e9mentaires: non pr\u00e9sentes.\n"); 
     }  // End of: if (supNO) 

     // Check the case of supplemenatary rows AND supplementary columns
     if (supRC) {
	 System.out.println("Sorry: the case of simultaneously both");
	 System.out.println("supplementary rows and supplementary columns");
	 System.out.println("is not currently supported.");
	 System.exit(1);
     }

     int indx = 0; 

     // Check the case of supplementary rows
     if (supR) {
	 // Some checks 
	 if (n != norig) {
	     System.out.println
             ("Problem: we have supplementary rows but n not equal to norig.");
	     System.out.println("Aborting.");
	     System.exit(1); 
	 }
	 if (m != morig + 1) {
            System.out.println
	 ("Problem: we have supplementary rows but m not equal to morig+1."); 
	    System.out.println("Aborting.");
	    System.exit(1);
	 }
	 mainRowCount = 0; 
	 suppRowCount = 0; 
         for (int i = 0; i < n; i++) { // Here, n,m and not norig,morig used.
	     // Check last column: is it binary?  Implies: 1 = main, 0 = suppl.
	     if (indat[i][m-1] == 0) suppRowCount++;
             if (indat[i][m-1] == 1) mainRowCount++;
         }
         if ( suppRowCount + mainRowCount != n) {    // And here n = norig
             System.out.println
	     ("Problem: with suppl. rows we are finding an inconsistency.");
	     System.out.println
	     ("We find: suppRowCount + mainRowCount not equal to n.");
	     System.out.println("Aborting.");
	     System.exit(1);
         }
	 mainColCount = morig; 
	 suppColCount = 0; 
         outext5.append
            ("Supplementary rows: present.\n"); 
	 // ("Lignes suppl\u00e9mentaires: pr\u00e9sentes.\n"); 
         outext5.append 
           ("Number of supplementary rows: " + suppRowCount + "\n");
	 // ("Nombre de lignes suppl\u00e9mentaires: " + suppRowCount + "\n");
     }  // end of: if (supR) 

     // Check the case of supplementary columns
     if (supC) {
	 // Some checks 
	 if (m != morig) {
	     System.out.println
             ("Problem: we have supplementary cols but m not equal to morig.");
	     System.out.println("Aborting.");
	     System.exit(1); 
	 }
	 if (n != norig + 1) {
            System.out.println
	 ("Problem: we have supplemenatary cols but n not equal to norig+1."); 
	    System.out.println("Aborting.");
	    System.exit(1);
	 }
	 mainColCount = 0; 
	 suppColCount = 0; 
         for (int j = 0; j < m; j++) { // Here, n,m and not norig,morig used.
	     // Check last row: is it binary?  Implies: 1 = main, 0 = suppl.
	     if (indat[n-1][j] == 0) suppColCount++;
             if (indat[n-1][j] == 1) mainColCount++;
         }
         if ( suppColCount + mainColCount != m) {
             System.out.println
	     ("Problem: with suppl. cols. we are finding an inconsistency.");
	     System.out.println
	     ("We find: suppColCount + mainColCount not equal to m.");
	     System.out.println("Aborting.");
	     System.exit(1);
         }
	 mainRowCount = norig; 
	 suppRowCount = 0; 
         outext5.append
          ("Supplementary columns: present.\n"); 
         // ("Colonnes suppl\u00e9mentaires: pr\u00e9sentes.\n"); 
         outext5.append 
          ("Number of supplementary columns: present: " + suppColCount + "\n");
        // ("Nombre de colonnes suppl\u00e9mentaires: " + suppColCount + "\n");
     }  // end of: if (supC) 


     double[][] indatstd   = new double[mainRowCount][mainColCount];
     String[] rowlabs      = new String[mainRowCount];
     String[] collabs      = new String[mainColCount];
     double[][] suppRowDat = new double[suppRowCount][mainColCount]; 
     double[][] suppColDat = new double[mainRowCount][suppColCount]; 
     String[] suppRowlabs  = new String[suppRowCount];
     String[] suppCollabs  = new String[suppColCount]; 
      
     // By default, take 1st main*Count values from *labsOrig
     for (int i = 0; i < mainRowCount; i++) {
	 rowlabs[i] = rowlabsOrig[i]; 
     }
     for (int j = 0; j < mainColCount; j++) {
	 collabs[j] = collabsOrig[j]; 
     }

     if (supNO) {
        indatstd = Copy(mainRowCount, mainColCount, indat);
	for (int i = 0; i < n; i++) {
	    rowlabs[i] = rowlabsOrig[i];
	}
	for (int j = 0; j < m; j++) {
	    collabs[j] = collabsOrig[j]; 
	}
     }

     if (supR) { 
	 indx = 0; 
	 for (int i = 0; i < n; i++) {
	     if (indat[i][m-1] == 1) {  // Note: check for 1 => principal row.
 	         rowlabs[indx] = rowlabsOrig[i]; 
		 for (int j = 0; j < m-1; j++) { // Omit last col = supp ind
		     indatstd[indx][j] = indat[i][j];
		 }
		 indx++;
	     }
	 }
	 indx = 0; 
	 for (int i = 0; i < n; i++) {
	     if (indat[i][m-1] == 0) {   // Note: check for 0 => suppl. row.
		 suppRowlabs[indx] = rowlabsOrig[i];  
		 for (int j = 0; j < m-1; j++) {  // Omit last col = supp ind
		     suppRowDat[indx][j] = indat[i][j];
		 }
		 indx++;
	     }
	 }
     }  // End of: if (supR) 

     if (supC) {
	 indx = 0; 
	 for (int j = 0; j < m; j++) {
	     if (indat[n-1][j] == 1) {  // Note: check for 1 => principal row.
	         collabs[indx] = collabsOrig[j]; 
		 for (int i = 0; i < n-1; i++) {  // Omit last row = supp ind
		     indatstd[i][indx] = indat[i][j];
		 }
		 indx++;
	     }
	 }
	 indx = 0; 
	 for (int j = 0; j < m; j++) {   // m = mainColCount = morig
	     if (indat[n-1][j] == 0) {   // Note: check for 0 => suppl. row. 
		 suppCollabs[indx] = collabsOrig[j]; 
		 for (int i = 0; i < n-1; i++) {  // Omit last row = supp ind
		     suppColDat[i][indx] = indat[i][j];
		 }
		 indx++;
	     }
	 }
     } // End of: if (supC) 


     // For convenience put n, m <-- mainRowCount, mainColCount
     n = mainRowCount; 
     m = mainColCount; 


     /* Step 6 ------------------------------------------------------
      * Finally, principal array.
      * ------------------------------------------------------------- */
       
     double[] rowsums = new double[n];
     double[] colsums = new double[m];
     double total = 0.0; 

     // Row sums and overall total 
     for (int i = 0; i < n; i++) {
         rowsums[i] = 0.0;
         for (int j = 0; j < m; j++) {
             rowsums[i] += indatstd[i][j];
             total += indatstd[i][j];
         }  
     }         

     // Col sums 
     for (int j = 0; j < m; j++) {
         colsums[j] = 0.0;
         for (int i = 0; i <n; i++) colsums[j] += indatstd[i][j];
     }

     // Finalize normalization to provide masses by dividing by total
     for (int i = 0; i < n; i++) rowsums[i] /= total;
     for (int j = 0; j < m; j++) colsums[j] /= total;
     for (int i = 0; i < n; i++) {
	 for (int j = 0; j < m; j++) indatstd[i][j] /= total;
     }


     /* Step 7 ------------------------------------------------------
      * Set up "principal", an object of class Imatrix, and 
      * potentially, either one of: "supplementaryRows" or 
      * "supplementaryColumns", each an object of class Imatrix.
      * ------------------------------------------------------------- */ 

     // Now set up "principal" object of class Imatrix:
     // this will be the input data for principal rows and columns.
     // Parameters for Imatrix: 
     // data, nrows, ncols, row-masses, col-masses, grand total,
     // row-labels, col-labels. 
     Imatrix principal = new Imatrix(indatstd, n, m, title, filename, 
			 rowsums, colsums, total, rowlabs, collabs); 

     // Case of supplementary rows
     // Note: only allow either suppl. rows or suppl. cols. in one analysis

     double[] rowsumsSuppR = new double[suppRowCount];
     Imatrix supplementaryRows = principal;   // Arbitrary initialization

     if (supR) {
	// Determine row sums
	for (int i = 0; i < suppRowCount; i++) {
	    rowsumsSuppR[i] = 0.0;
	    for (int j = 0; j < m; j++) 
		rowsumsSuppR[i] += suppRowDat[i][j]; 
	}
	for (int i = 0; i < suppRowCount; i++) rowsumsSuppR[i] /= total;
	for (int i = 0; i < suppRowCount; i++) {
	    for (int j = 0; j < m; j++) suppRowDat[i][j] /= total;
	}

	// Set up supplementaryRows object of class Imatrix: 
	// this will be used as input data for supplementary rows. 
	// Parameters for Imatrix: 
	// data, nrows, ncols, row-masses, col-masses, grand total,
	// row-labels, col-labels.
	Imatrix xsupplementaryRows = new Imatrix(suppRowDat, suppRowCount,
						m, title, filename, 
                                                rowsumsSuppR, 
						colsums, total, suppRowlabs, 
						collabs);  
	supplementaryRows = xsupplementaryRows; 
     }  // End of: if (supR)


     // Case of supplementary columns
     // Note: only allow either suppl. rows or suppl. cols. in one analysis

     double[] colsumsSuppC = new double[suppColCount];

     if (supC) {
	// Determine column sums
	for (int j = 0; j < suppColCount; j++) {
	    colsumsSuppC[j] = 0.0;
	    for (int i = 0; i < n; i++) 
		colsumsSuppC[j] += suppColDat[i][j]; 
	}
	for (int j = 0; j < suppColCount; j++) colsumsSuppC[j] /= total;
	for (int j = 0; j < suppColCount; j++) {
	    for (int i = 0; i < n; i++) suppColDat[i][j] /= total;
	}

	// Set up supplementaryColumns object of class Imatrix: 
	// this will be used as input data for supplementary columns.  
	// Parameters for Imatrix: 
	// data, nrows, ncols, row-masses, col-masses, grand total,
	// row-labels, col-labels. 
	Imatrix supplementaryColumns = new Imatrix(suppColDat, n, suppColCount,
  				   title, filename, rowsums, 
                                   colsumsSuppC, total, rowlabs, suppCollabs);
     } // End of: if (supC)

     // Some information output
     outext5.append
         ("Preprocessing of the data ... \n"); 
     //  ("Pr\u00e9traitement sur les donn\u00e9es ... \n"); 

       
     /* Step 8 --------------------------------------------------------
      * Into correspondence analysis now.
      * --------------------------------------------------------------- */

     // ---------------------------------------------------------------
     // Make use of the various methods in the Correspondence Analysis
     // class, viz. class CorrAnal
     // These methods are: diagonalization, projections, contributions,
     // correlations, and formatprint; and for supplementary elements:
     // projectionsSupp, contributionsSupp, correlationsSupp, and 
     // formatprintSupp.
     // ---------------------------------------------------------------

     CorrAnal MyCA; 
     MyCA = new CorrAnal();

     double[] Evals = new double[m];               // Eigenvalues
     double[][] Evex = new double[m][m];           // Eigenvectors
     double[][] CP = new double[m][m];             // Cross-products 
     double[] rate = new double[m];                // Rates of inertia
     double trce;                                  // trace

     trce = MyCA.diagonalization(principal, CP, outext1, Evals, Evex, rate);
     outext5.append
        ("Diagonalization complete ...\n"); 
     // ("Diagonalisation accomplie ...\n"); 

     // We collect all eigen-analysis results in an object of class ResultSVD
     ResultSVD resEigenAnalysis = new ResultSVD(Evex, Evals, CP, rate, trce);

     double[][] rowproj = new double[n][m];        // row projections
     double[][] colproj = new double[m][m];        // column projections

     MyCA.projections(principal, resEigenAnalysis, rowproj, colproj);
     outext5.append
        ("Projections calculated ... \n"); 
     // ("Projections calcul\u00e9es ... \n"); 

     double[][] rowcntr = new double[n][m];        // row contributions
     double[][] colcntr = new double[m][m];        // column contributins
     MyCA.contributions(principal, rowproj, colproj, rowcntr, colcntr);
     outext5.append
        ("Contributions calculated ... \n"); 
     // ("Contributions calcul\u00e9es ... \n"); 

     double[][] rowcorr = new double[n][m];        // row correlations 
     double[][] colcorr = new double[m][m];        // column correlations 
     MyCA.correlations(principal, rowproj, colproj, rowcorr, colcorr); 
     outext5.append
        ("Correlations calculated ... \n"); 
     // ("Corr\u00e9lations calcul\u00e9es ... \n"); 

     if ( (m-1) < nfactors) nfactors = m -1; 
     double[] cinertia = new double[m];    // relative inertia
     MyCA.formatprint (principal, nfactors, rowproj, colproj, cinertia, 
                       rowcntr, colcntr, rowcorr, colcorr, trce, outext1); 
     outext5.append
         ("Results output to screen ... \n"); 
     //  ("R\u00e9sultats imprim\u00e9s \u00e0 l'\u00e9cran ... \n"); 
	
     //---------------------------------------------------SUPPL ROWS
     // Projections on suppl. rows
     double[][] rowprojS = new double[suppRowCount][m];     
     double[][] rowcorrS = new double[suppRowCount][m];     
     double[][] rowcntrS = new double[suppRowCount][m];     
	
     if (supR) {
	 // Case of supplementary rows
    
         MyCA.projectionsSupp(m, supplementaryRows, Evals, colproj, rowprojS);
         outext5.append
            ("Projections of supplementary rows calculated ... \n"); 
         // ("Projections des lignes suppl. calcul\u00e9es ... \n"); 
    
	 // Correlations on suppl. rows
         MyCA.correlationsSupp(m, supplementaryRows, rowprojS, rowcorrS); 
         outext5.append
            ("Correlations of supplementary rows calculated ... \n"); 
         // ("Corr\u00e9lations des lignes suppl. calcul\u00e9es ... \n"); 

	 // Contributions on suppl. rows
         MyCA.contributionsSupp(m, supplementaryRows, rowprojS, 
                                     rowcntrS, Evals);
         outext5.append
            ("Contributions of supplementary rows calculated ... \n"); 
         // ("Contributions des lignes suppl. calcul\u00e9es ... \n"); 
     
	 outext1.append
            ("Supplementary rows: \n");
         // ("Lignes suppl\u00e9mentaires: \n");
         MyCA.formatprintSupp (m, nfactors, supplementaryRows, 
                       rowprojS, 
                       rowcntrS, rowcorrS, trce,
                       outext1); 

     }
	                            
     // -------------------------------------------------------SUPPL COLS
     // Projections on suppl. columns
     double[][] colprojS = new double[suppColCount][m]; 
     double[][] colcorrS = new double[suppColCount][m]; 
     double[][] colcntrS = new double[suppColCount][m]; 

     if (supC) {
	 // Supplementary columns
	 // We have: suppColDat double[mainRowCount][suppColCount]
	 // and      suppCollabs String[suppColcount] 

	 // First determine transpose of suppl. col. data
	 double[][] suppColDatTrans = new double[suppColCount][n];
	 for (int i = 0; i < n; i++) {
	     for (int j = 0; j < suppColCount; j++) {
		 suppColDatTrans[j][i] = suppColDat[i][j];
	     }
	 }
	 Imatrix supplementaryColumnsTransposed = new Imatrix(suppColDatTrans,
                 suppColCount, n, title, filename,colsumsSuppC, rowsums, 
                 total, suppCollabs, rowlabs);

	 MyCA.projectionsSupp(m, supplementaryColumnsTransposed, 
               Evals, rowproj, colprojS);
         outext5.append
            ("Projections of supplementary columns calculated ... \n"); 
         // ("Projections des colonnes suppl. calcul\u00e9es ... \n"); 

	 // Correlations on suppl. columns
	 MyCA.correlationsSupp(m, supplementaryColumnsTransposed, 
	         colprojS, colcorrS); 
         outext5.append
             ("Correlations of supplementary columns calculated ... \n"); 
         //  ("Corr\u00e9lations des colonnes suppl. calcul\u00e9es ... \n"); 

	 MyCA.contributionsSupp(m, supplementaryColumnsTransposed, colprojS, 
	                           colcntrS, Evals);
         outext5.append
             ("Contributions of supplementary columns calculated ... \n"); 
         //  ("Contributions des lignes suppl. calcul\u00e9es ... \n"); 

	 outext1.append
             ("Supplementary columns: \n");
         //  ("Colonnes suppl\u00e9mentaires: \n");

         MyCA.formatprintSupp (m, nfactors, supplementaryColumnsTransposed, 
	               colprojS, 
		       colcntrS, colcorrS, trce,
                       outext1); 

         outext5.append
       ("Projections on supplementary elements calculated ... \n"); 
    // ("Projections sur les \u00e9l\u00e9ments supp. calcul\u00e9es ... \n"); 

     }
       
     /* Step 8 ---------------------------------------------------------
      * Collect correspondence analysis results into one object.
      * ---------------------------------------------------------------- */

     String[] factorlabs = new String[m];
     for (int j = 0; j < m; j++) factorlabs[j] = String.valueOf(j);
     double[] dummy1 = new double[n];
     for (int i = 0; i < n; i++) dummy1[i] = 1.0;
     double[] dummy2 = new double[m];
     for (int j = 0; j < m; j++) dummy2[j] = 1.0;
     Imatrix irowproj = new Imatrix(rowproj, n, m, title, filename,
			rowsums, dummy2, total, rowlabs, factorlabs); 
     Imatrix icolproj = new Imatrix(colproj, m, m, title, filename,
			colsums, dummy1, total, collabs, factorlabs); 
     Imatrix irowcorr = new Imatrix(rowcorr, n, m, title, filename,
			rowsums, dummy2, total, rowlabs, factorlabs); 
     Imatrix icolcorr = new Imatrix(colcorr, m, m, title, filename,
			colsums, dummy1, total, collabs, factorlabs); 
     Imatrix irowcntr = new Imatrix(rowcntr, n, m, title, filename,
			rowsums, dummy2, total, rowlabs, factorlabs); 
     Imatrix icolcntr = new Imatrix(colcntr, m, m, title, filename,
			colsums, dummy1, total, collabs, factorlabs); 
     Imatrix irowprojS = new Imatrix(rowprojS, suppRowCount, m, title, 
        filename, rowsumsSuppR, dummy2, total, suppRowlabs, factorlabs); 
     Imatrix icolprojS = new Imatrix(colprojS, suppColCount, m, title, 
        filename, colsumsSuppC, dummy1, total, suppCollabs, factorlabs); 
     Imatrix irowcorrS = new Imatrix(rowcorrS, suppRowCount, m, title, 
        filename, rowsumsSuppR, dummy2, total, suppRowlabs, factorlabs); 
     Imatrix icolcorrS = new Imatrix(colcorrS, suppColCount, m, title, 
        filename, colsumsSuppC, dummy1, total, suppCollabs, factorlabs); 
     Imatrix irowcntrS = new Imatrix(rowcntrS, suppRowCount, m, title, 
        filename, rowsumsSuppR, dummy2, total, suppRowlabs, factorlabs); 
     Imatrix icolcntrS = new Imatrix(colcntrS, suppColCount, m, title, 
        filename, colsumsSuppC, dummy1, total, suppCollabs, factorlabs); 

     ResultCA myCAresult = new ResultCA( irowproj, icolproj, irowcorr,
					 icolcorr, irowcntr, icolcntr,
					 irowprojS, icolprojS, irowcorrS,
					 icolcorrS, irowcntrS, icolcntrS);


     /* Step 8 ---------------------------------------------------------
      * Produce and output plots.
      * ---------------------------------------------------------------- */

     // ----------------------------------------------------------------
     // Plot factor planes... Factors 1 and 2, 1 and 3, 2 and 3.
     // ----------------------------------------------------------------

     // To be plotted - observation projections and variable projections
     int nRows = n + suppRowCount; 
     int nCols = m + suppColCount; 

     double[] obsx = new double[nRows];
     double[] obsy = new double[nRows]; 
     double[] vbex = new double[nCols];
     double[] vbey = new double[nCols]; 
     double tauxx, tauxy; 
     int axisx, axisy; 
     String[] allrowlabs = new String[nRows];
     String[] allcollabs = new String[nCols]; 

     // ----------------------------------------------------------------
     // Factors 1 and 2
     // ----------------------------------------------------------------
     axisx = 1; 
     axisy = 2; 

     preparePlotData (n, suppRowCount, m, suppColCount, 
		      rowproj, colproj, rowprojS, colprojS, 
                      obsx, obsy, vbex, vbey,  
		      rowlabs, collabs, suppRowlabs, suppCollabs, 
                      allrowlabs, allcollabs, axisx, axisy); 
 
     tauxx = rate[axisx];  // rate of inertia explained by axis
     tauxy = rate[axisy];  // rate of inertia explained by axis 
     plotGraph (nRows, nCols, obsx, obsy, vbex, vbey, allrowlabs, allcollabs, 
               axisx, axisy, tauxx, tauxy, outext1f12);
     outext5.append
        ("Plane of projections on the first and second factors ... \n");
     // ("Plan des premiers 1er, 2e facteurs trac\u00e9 ... \n");

     // ----------------------------------------------------------------
     // Factors 1 and 3
     // ----------------------------------------------------------------
     axisx = 1; 
     axisy = 3; 

     preparePlotData (n, suppRowCount, m, suppColCount, 
		      rowproj, colproj, rowprojS, colprojS, 
                      obsx, obsy, vbex, vbey,  
		      rowlabs, collabs, suppRowlabs, suppCollabs, 
                      allrowlabs, allcollabs, axisx, axisy); 
 
     tauxx = rate[axisx];  // rate of inertia explained by axis
     tauxy = rate[axisy];  // rate of inertia explained by axis 
     plotGraph (n, m, obsx, obsy, vbex, vbey, rowlabs, collabs, 
               axisx, axisy, tauxx, tauxy, outext1f13);
     outext5.append
        ("Plane of projections on the first and third factors ... \n");
     // ("Plan des premiers 1er, 3e facteurs trac\u00e9 ... \n");

     // ----------------------------------------------------------------
     // Factors 2 and 3
     // ----------------------------------------------------------------
     axisx = 2; 
     axisy = 3; 

     preparePlotData (n, suppRowCount, m, suppColCount, 
		      rowproj, colproj, rowprojS, colprojS, 
                      obsx, obsy, vbex, vbey,  
		      rowlabs, collabs, suppRowlabs, suppCollabs, 
                      allrowlabs, allcollabs, axisx, axisy); 
 
     tauxx = rate[axisx];  // rate of inertia explained by axis
     tauxy = rate[axisy];  // rate of inertia explained by axis 
     plotGraph (n, m, obsx, obsy, vbex, vbey, rowlabs, collabs, 
               axisx, axisy, tauxx, tauxy, outext1f23);
     outext5.append
        ("Plane of projections on the second and third factors ... \n");
     // ("Plan des premiers 2e, 3e facteurs trac\u00e9 ... \n");



     /* Step 9 -----------------------------------------------------------
      * Hierarchical clustering.   The work is done in methods 
      * reciprocalNearestNeighbors and printClustering.
      * ------------------------------------------------------------------ */

     //-------------------------------------------------------------------
     // Next stage: CAH (hierarchical clustering, agglomerative 
     // clustering, classification ascendante hierarchique) follows now.
     //-------------------------------------------------------------------

     // Working on factor projections implies normalization: 
     // remember that we carried out the correspondence analysis on 
     // profiles in the chi-squared metric - not on the raw data. 
     // The factor space is implicitly Euclidean.  Which helps in display.

     // nfactors must be <= m.  If not the case, we will force it to be.
     if (nfactors > m) nfactors = m; 
     // if nfactors is 0 or negative, we will force it to be min(m, 7).
     if (nfactors < 1) nfactors = Math.min(m, 7); 
 

     double[][] indat2clstr = new double[n][nfactors];
     for (int i = 0; i < n; i++) {
	 // Note: we avoid the j=0 trivial factor consisting of all 1s.
	 // In case we have one column too few, assuming we are using 
	 // all factors, since we omitted the trivial 1-valued factor - 
	 // we will fill the last column of indat2clstr with 0s.  
	 indat2clstr[i][nfactors-1] = 0.0;
	 if ( (nfactors+1) <= m) {
	    for (int j = 1; j < (nfactors+1); j++) {
	        // Coordinates on each factor:
	        indat2clstr[i][j-1] = rowproj[i][j];
	    }
         }
	 else {
	    for (int j = 1; j < nfactors; j++) {
	        // Coordinates on each factor:
	        indat2clstr[i][j-1] = rowproj[i][j];
	    }
	 }
     }

     outext2.append 
           ("Inertias and level indices are calculated in\n");
           // ("Les inerties et indices de niveau sont calcul\u00e9s dans\n");
     outext2.append
           ("the space definied by the " + nfactors + 
                " axes used for the CAH. \n"); 
           //  ("l'espace endendr\u00e9 par les " + nfactors + 
           //    " axes utilis\u00e9s pour la CAH. \n"); 

     HierClass HC;  
     HC = new HierClass();

     // Recall: there are n-1 levels in a hierarchy, or n if we 
     // include the set of all singletons.  
     int[][] clusters = new int[n][n];        // cluster label matrix
     int[][] clusters2 = new int[n][n]; 
     // Acl, Bcl convention: cluster has seqno. = nrow+1, ... 2*nrow-1
     int[] Acl = new int[n-1];                // aine, elder
     int[] Bcl = new int[n-1];                // benjamin, younger
     double[] levels = new double[n-1];       // cluster levels 
     double[] card = new double[n-1];         // cluster cardinalities
     int[] order = new int[n];                // dendrogram-consistent
                                              // order of obs. 
     double totint;                           // total inertia
     double[] mass = new double[n];           // masses  
     // Mass values will be changed following agglomerations.  
     for (int i = 0; i < n; i++) mass[i] = rowsums[i]; 

     HC.reciprocalNearestNeighbors
        (n, nfactors, indat2clstr, mass, outext2, clusters, 
	 clusters2, Acl, Bcl, levels, card);
     outext5.append
        ("Classification built ...\n"); 
     // ("Classification b\u00e2tie ...\n"); 

     totint = HC.printClustering(n, levels, outext2, clusters, clusters2, 
         card, Acl, Bcl, rowlabs, order); 
     outext5.append
	 ("Results output to screen ... \n"); 
     //  ("R\u00e9sultats imprim\u00e9s \u00e0 l'\u00e9cran ... \n"); 



     /* Step 10 -------------------------------------------------------
      * Interpretion of factors in terms of clusters.
      * Method FacorVacor is used.  
      * --------------------------------------------------------------- */

     //----------------------------------------------------------------
     // Next - on to FACOR and VACOR interpretation.
     //----------------------------------------------------------------

     // Note: Selected no. of clusters for analysis is: nclusters

     Interpret interprt; 
     interprt = new Interpret(); 
     interprt.FacorVacor(n, m, nclusters, nfactors, rowsums, colsums, 
		       indatstd, cinertia, order, clusters, Acl, Bcl, 
                       Evex, Evals, totint, collabs, outext3, outext4); 
     outext5.append
	 ("Results of FACOR, and of VACOR on J, output ... \n");  
      // ("R\u00e9sultats de FACOR, et de VACOR sur J, sortis ... \n");  

     outext5.append
         ("Analysis complete.\n"); 
     //  ("Analyse termin\u00e9e.\n"); 
 
        }


	/* Step 11 -----------------------------------------------------
	 * Catch some exceptions and do some (limited) explaining of them.
	 * ------------------------------------------------------------- */

	catch (java.io.FileNotFoundException e) {
            // java.io.FileNotFoundException potentially thrown by:
            // FileInputStream is = new FileInputStream(pathname);
	    System.err.println("Error: " + e); 
	    System.err.println("Check for existence of file.");
	}
	catch (java.io.IOException e) {
            // java.io.IOException potentially thrown by:
            // title = bis.readLine();  Or by:
            // st.nextToken();     (Various occurrences.)
	    System.err.println("Error: " + e); 
	    System.err.println
              ("Check for numeric/char format problem in input data."); 
	}
	catch (java.lang.ArrayIndexOutOfBoundsException e) {
	    System.err.println("Error: " + e); 
	    System.err.println
              ("Check validity of the data in the input file.");
	}
	catch (Exception e) {
	    System.err.println("Error: " + e); 
	}
	finally {
	    System.err.println("All processing completed."); 
	    System.err.println("(c) f.murtagh@qub.ac.uk");
	}

} // end of main

    //-------------------------------------------------------------------  

    /** 
     * Little method for helping in output formating <br>
     * @param n     total length of field to be padded with blanks
     */

    public static String getSpaces(int n) {

    StringBuffer sb = new StringBuffer(n);
    for (int i = 0; i < n; i++) sb.append(' ');
    return sb.toString();
    } // getSpaces
   
    //-------------------------------------------------------------------

    /**
     * Method for printing a matrix, values scaled up by scaling factor <br>
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values
     * @param collo start column for printing 
     * @param colhi end column for printing 
     */

    public static void printMatrix(int n1, int n2, double[][] m, 
           int collo, int colhi, double scaling, JTextArea outext)
    {
      // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      // Following suppresses e.g. comma in 1,000 = 1000 for English locale.
      myFormat.setGroupingUsed(false);

      if (n2 < colhi) colhi = n2;   // Case of small no. of variables. 

      int w; 
      double maxval; 
      maxval = SMALL;
      for (int i =0; i < n1; i++) {
          for (int j =0; j < n2; j++) {
              // Scale up values by 'scaling' factor
              if (m[i][j] > maxval) maxval = m[i][j];
          }
      }
      
      String temp = myFormat.format(maxval*scaling, 
                      new StringBuffer(), fp).toString();
      // Output display field width of each number: extra 2 to 
      // account for possible minus sign, and 1 space.
      w = fp.getEndIndex() - fp.getBeginIndex() + 2;
      // System.out.println(" Output field width = " + w);

      // In the next few lines we say how we want numbers to appear.
      // max integer digits = max no. of digits before dec. point
      // max fraction digits = max no. of digits following dec. point
      // min fraction digits = min no. of digits following dec. point
      myFormat.setMaximumIntegerDigits(w);
      // We will set min and max nos. of digits following dec. point to 0
      myFormat.setMaximumFractionDigits(0);
      myFormat.setMinimumFractionDigits(0);
      for (int i=0; i<n1; i++)
	  {
	      // Print each row, elements separated by spaces
              for (int j = (collo-1); j < colhi; j++)
		  // Following unfortunately doesn't format at all
		  //                  System.out.print(m[i][j] + "  ");
                  {
                  String valString = myFormat.format(
                       scaling*m[i][j], new StringBuffer(), fp).toString();
		  // With a max field width of w, we pack spaces before val.
                  valString = getSpaces(w - fp.getEndIndex()) + valString;
                  outext.append(valString);
                  }
          }
      outext.append("\n");
    } // printMatrix


    //-------------------------------------------------------------------

    /**
     * Method for printing all, values scaled up by scaling factor <br>
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values
     */

    public static void printAll(int n1, int n2, double[][] m, 
			   String[] labs, double scaling, JTextArea outext)
    {
      // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      // Following suppresses e.g. comma in 1,000 = 1000 for English locale.
      myFormat.setGroupingUsed(false);

      int w; 
      double maxval; 
      maxval = SMALL;

      for (int i =0; i < n1; i++) {
          for (int j =0; j < n2; j++) {
              // Scale up values by 'scaling' factor
              if (m[i][j] > maxval) maxval = m[i][j];
          }
      }

      String temp = myFormat.format(maxval*scaling, 
                      new StringBuffer(), fp).toString();
      // Output display field width of each number: extra 2 to 
      // account for possible minus sign, and 1 space.
      w = fp.getEndIndex() - fp.getBeginIndex() + 2;
      // System.out.println(" Output field width = " + w);

      // In the next few lines we say how we want numbers to appear.
      // max integer digits = max no. of digits before dec. point
      // max fraction digits = max no. of digits following dec. point
      // min fraction digits = min no. of digits following dec. point
      myFormat.setMaximumIntegerDigits(w);
      // We will set min and max nos. of digits following dec. point to 0
      myFormat.setMaximumFractionDigits(0);
      myFormat.setMinimumFractionDigits(0);
      for (int i=0; i<n1; i++)
	  {
	      // First handle cols for labels, QLT, PDS, INR:
	      String myString = labs[i];
	      myString = getSpaces(4 - myString.length()) + myString; 
              outext.append("|" + myString + "|");
	      for (int j = 0; j < 3; j++) {
                  String valString = myFormat.format(
                       scaling*m[i][j], new StringBuffer(), fp).toString();
		  // With a max field width of w, we pack spaces before val.
                  valString = getSpaces(4 - fp.getEndIndex()) + valString;
                  outext.append(valString);
	      }

	      // Print each row, elements separated by spaces
              for (int j = 3; j < n2; j++)
                  {
		      // Scaling is in thousandths; but note if +ve only.
		      w = 4;      // All cntr and corr are positive.
		      if (j/3 == (double)j/3.0) {  // Here, will handle proj.
			  outext.append("|");
			  w = 5;                   // Allow for -ve vals.
		      }
                      String valString = myFormat.format(
                          scaling*m[i][j], new StringBuffer(), fp).toString();
		      // With max field width of w, we pack spaces before val.
                      valString = getSpaces(w - fp.getEndIndex()) + valString;
                      outext.append(valString);
                  }
	      outext.append("|");
              // Start a new line at the end of a row
	      outext.append("\n");
          }
    } // printAll


    //-------------------------------------------------------------------

    /**
     * Method for printing all, values scaled up by scaling factor <br>
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values
     */

    public static void printAll(int n1, int n2, double[][] m, 
		  int[][] cluslabs, double scaling, JTextArea outext)
    {
      // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      // Following suppresses e.g. comma in 1,000 = 1000 for English locale.
      myFormat.setGroupingUsed(false);

      int w; 
      double maxval; 
      maxval = SMALL;

      for (int i =0; i < n1; i++) {
          for (int j =0; j < n2; j++) {
              // Scale up values by 'scaling' factor
              if (m[i][j] > maxval) maxval = m[i][j];
          }
      }

      String temp = myFormat.format(maxval*scaling, 
                      new StringBuffer(), fp).toString();
      // Output display field width of each number: extra 2 to 
      // account for possible minus sign, and 1 space.
      w = fp.getEndIndex() - fp.getBeginIndex() + 2;
      // System.out.println(" Output field width = " + w);

      // In the next few lines we say how we want numbers to appear.
      // max integer digits = max no. of digits before dec. point
      // max fraction digits = max no. of digits following dec. point
      // min fraction digits = min no. of digits following dec. point
      myFormat.setMaximumIntegerDigits(w);
      // We will set min and max nos. of digits following dec. point to 0
      myFormat.setMaximumFractionDigits(0);
      myFormat.setMinimumFractionDigits(0);
      String myString;
      for (int i=0; i<n1; i++)
	  {
	      // First cluslabs
	      for (int j = 0; j < 3; j++) {
                  myString = myFormat.format(
                       cluslabs[i][j], new StringBuffer(), fp).toString();
		  // With a max field width of 5, we pack spaces before val.
                  myString = getSpaces(5 - fp.getEndIndex()) + myString;
                  outext.append(myString);
	      }
	      outext.append("|");
	      for (int j = 0; j < 3; j++) {
                  myString = myFormat.format(
                       scaling*m[i][j], new StringBuffer(), fp).toString();
		  // With a max field width of w, we pack spaces before val.
                  myString = getSpaces(4 - fp.getEndIndex()) + myString;
                  outext.append(myString);
	      }

	      // Print each row, elements separated by spaces
              for (int j = 5; j < n2; j++)
                  {
		      // Scaling is in thousandths; but note if +ve only.
		      w = 4;      // All cntr and corr are positive.
		      if ((j-2)/3 == (double)(j-2)/3.0) {  //Here, handle proj.
			  outext.append("|");      // New factor divider. 
			  w = 5;                   // Allow for -ve vals.
		      }
                      myString = myFormat.format(
                          scaling*m[i][j], new StringBuffer(), fp).toString();
		      // With max field width of w, we pack spaces before val.
                      myString = getSpaces(w - fp.getEndIndex()) + myString;
                      outext.append(myString);
                  }
	      outext.append("|");
              outext.append("\n");
          }
    } // printAll


    //-------------------------------------------------------------------

    /**
     * Method printVect for printing a vector <br>
     * @param m input vector of length m.length
     * @param collo start column number, usually 1 (= 2nd, not 0 = 1st)
     * @param colhi end column number, usually 7 by convention
     * @param abswidth display precision, fixed constant
     * @param scaling scaling factor (e.g. 100 or 1000)
     * @param label description label of row of values
     */

    public static void printVect(double[] m, int collo, int colhi, 
           int abswidth, double scaling, String label, JTextArea outext)
    {
      // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      // Following suppresses e.g. comma in 1,000 = 1000 for English locale.
      myFormat.setGroupingUsed(false);

      int w; 
      double maxval; 
      maxval = SMALL;

      // Test colhi against m.length
      // Our indexing runs from 1 since we ignore 0/trivial eigenvalue
      if ( (m.length-1) < colhi) colhi = m.length-1; 

      for (int i = 0; i < m.length; i++) {
          // Scale up values by 'scaling' factor
          if (m[i] > maxval) maxval = m[i];
      }

      String temp = myFormat.format(maxval*scaling, 
                      new StringBuffer(), fp).toString();
      // Output display field width of each number: extra 2 to 
      // account for possible minus sign, and 1 space.
      w = fp.getEndIndex() - fp.getBeginIndex() + 2;
      // System.out.println(" Output field width = " + w);
      // abswidth is our preferred field (or numeric token) width,
      // but if our fields require more space, we will grant this: 
      if ( (w-1) > abswidth) abswidth = w; 

      // In the next few lines we say how we want numbers to appear.
      // max integer digits = max no. of digits before dec. point
      // max fraction digits = max no. of digits following dec. point
      // min fraction digits = min no. of digits following dec. point
      myFormat.setMaximumIntegerDigits(w);
      myFormat.setMaximumFractionDigits(0);
      myFormat.setMinimumFractionDigits(0);


      outext.append(label);
      for (int i = collo; i <= colhi; i++)
	  {
             String valString = myFormat.format(
                   m[i]*scaling, new StringBuffer(), fp).toString();
	     // Replace abswidth with w in the following for adaptive spacing.
             valString = getSpaces(abswidth - fp.getEndIndex()) + valString;
		   outext.append(valString);
          }
      // Flush the line at the row end
      outext.append("\n");
    } // printVect


    //-------------------------------------------------------------------

    /**
     * Method printVect for printing a vector <br>
     * @param m input vector of length m.length
     * @param collo start column number, usually 1 (= 2nd, not 0 = 1st)
     * @param colhi end column number, usually 7 by convention
     * @param abswidth display precision, fixed constant
     * @param scaling scaling factor (e.g. 100 or 1000)
     * @param label description label of row of values
     */

    public static void printVect(int[] m, int collo, int colhi, 
           int abswidth, int scaling, String label, JTextArea outext)
    {
      // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      // Following suppresses e.g. comma in 1,000 = 1000 for English locale.
      myFormat.setGroupingUsed(false);

      int w; 
      int maxval; 
      maxval = 0;

      // Test colhi against m.length
      // Our indexing runs from 1 since we ignore 0/trivial eigenvalue
      if ( (m.length-1) < colhi) colhi = m.length-1; 

      for (int i = collo; i <= colhi; i++) {
          // Scale up values by 'scaling' factor
          if (m[i] > maxval) maxval = m[i];
      }

      String temp = myFormat.format(maxval*scaling, 
                      new StringBuffer(), fp).toString();
      // Output display field width of each number: extra 2 to 
      // account for possible minus sign, and 1 space.
      w = fp.getEndIndex() - fp.getBeginIndex() + 2;
      // System.out.println(" Output field width = " + w);
      // abswidth is our preferred field (or numeric token) width,
      // but if our fields require more space, we will grant this: 
      if ( (w-1) > abswidth) abswidth = w; 

      // In the next few lines we say how we want numbers to appear.
      // max integer digits = max no. of digits before dec. point
      // max fraction digits = max no. of digits following dec. point
      // min fraction digits = min no. of digits following dec. point
      myFormat.setMaximumIntegerDigits(w);
      // myFormat.setMaximumFractionDigits(0);
      // myFormat.setMinimumFractionDigits(0);


      outext.append(label);
      for (int i = collo; i <= colhi; i++)
	  {
             String valString = myFormat.format(
                   m[i]*scaling, new StringBuffer(), fp).toString();
	     // Replace abswidth with w in the following for adaptive spacing.
             valString = getSpaces(abswidth - fp.getEndIndex()) + valString;
		   outext.append(valString);
          }
      // Flush the line at the row end
      outext.append("\n");
    } // printVect


    //-------------------------------------------------------------------

    /**
     * Method printVect for printing a vector <br>
     * @param m input vector of length m.length
     * @param abswidth display precision, fixed constant
     * @param labels vector of row labels
     */

    public static void printVect(int[] m, int abswidth, String[] labels,
           String label, JTextArea outext)
    {
      // Some definitions for handling output formating
      NumberFormat myFormat = NumberFormat.getNumberInstance();
      FieldPosition fp = new FieldPosition(NumberFormat.INTEGER_FIELD);
      // Following suppresses e.g. comma in 1,000 = 1000 for English locale.
      myFormat.setGroupingUsed(false);

      int w; 
      // Test colhi against m.length
      int n;
      n = labels.length;
      int collo = 0;
      int colhi; 
      colhi = m.length;

      int maxval= -10000;
      for (int i = collo; i < colhi; i++) {
          if (m[i] > maxval) maxval = m[i];
      }
      String temp = myFormat.format(maxval, 
                      new StringBuffer(), fp).toString();
      // Output display field width of each number.
      w = fp.getEndIndex() - fp.getBeginIndex();
      // System.out.println(" Output field width = " + w);
      // abswidth is our preferred field (or numeric token) width,
      // but if our fields require more space, we will grant this: 
      if ( (w-1) > abswidth) abswidth = w; 

      // In the next few lines we say how we want numbers to appear.
      // max integer digits = max no. of digits before dec. point
      // max fraction digits = max no. of digits following dec. point
      // min fraction digits = min no. of digits following dec. point
      myFormat.setMaximumIntegerDigits(w);
      // myFormat.setMaximumFractionDigits(0);
      // myFormat.setMinimumFractionDigits(0);

      outext.append(label);
      String valString; 
      for (int i = collo; i < colhi; i++)
	  {
	     if (m[i] <= n) {
             valString = getSpaces(abswidth - labels[m[i]-1].length() ) + 
                         labels[m[i]-1];
	     }
	     else {
             valString = myFormat.format(
                   m[i], new StringBuffer(), fp).toString();
  	           valString = getSpaces(abswidth-fp.getEndIndex())+valString;
	     }             
	     // Replace abswidth with w in the following for adaptive spacing.
		   outext.append(valString);
          }
      // Flush the line at the row end
      outext.append("\n");
    } // printVect


    //-------------------------------------------------------------------

    /**
     * Method to sort a double vector, by inefficient straight insertion.<br>
     * See overloaded version taking integer array input for discussion <br>
     * of origin of this code.
     * @param invect input data to be sorted, sorted in place
     */

   public static void inSort (double[] invect) {
        double a;
        int i; 

        for (int j = 1; j < invect.length; j++) {
            a = invect[j]; 
            i = j - 1;
            while (i >= 0 && invect[i] > a) {
                invect[i+1] = invect[i]; 
                i--; 
	    }
            invect[i+1] = a; 
	}
	// Return type void
   }  // inSort


    //-------------------------------------------------------------------

    /**
     * Method to sort a double vector, by inefficient straight insertion.
     * <P>
     * See section 8.1, p 243, of Numerical Recipes in C. <br>
     * Note following difference: Numerical Recipes uses an array <br>
     * with elements 1..n whereas we use 1..length-1 <br>
     * Consequently the main loop here is <br>
     *         for (int j = 1; j < invect.length; j++) {
     * Was:    for (int j = 2; j < invect.length; j++) {
     * and the while loop here is <br>
     *         while (i >= 0 && invect[i] > a) {
     * Was:    while (i > 0 && invect[i] > a) {
     * This sort routine is O(n^2) and therefore very inefficient.
     * @param invect input data to be sorted, sorted in place
     */

   public static void inSort (int[] invect) {
        int a;
        int i; 

        for (int j = 1; j < invect.length; j++) {
            a = invect[j]; 
            i = j - 1;
            while (i >= 0 && invect[i] > a) {
                invect[i+1] = invect[i]; 
                i--; 
	    }
            invect[i+1] = a; 
	}
	// Return type void
   }  // inSort


    //-------------------------------------------------------------------

    /**
     * Method for straight copying of the input data <br>
     * @param nrow number of rows in input matrix
     * @param ncol number of columns in input matrix
     * @param A input matrix values
     * @param ncolnew new number of columns in transformed data
     */

    public static double[][] Copy(int nrow, int ncol, double[][] A)
    {
    // Adat will contain the copy of the data and will be returned
    double[][] Adat = new double[nrow][ncol];

     for (int j = 0; j < ncol; j++) {
         for (int i = 0; i < nrow; i++) {
	     Adat[i][j] = A[i][j]; 
	 }
     }

    return Adat;
    } // Copy


    //-------------------------------------------------------------------

    /**
     * Method for printing a double float matrix.  <P>
     * Based on ER Harold, "Java I/O", O'Reilly, around p 473. <br>
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values, double 
     * @param d display precision, number of decimal places
     * @param w display precision, total width of floating value
     */

    public static void printMatrix(int n1, int n2, double[][] m, int d, int w,
               JTextArea outext)
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
                  outext.append(valString);
                  }
              // Start a new line at the end of a row
	      outext.append("\n");
          }
    } // End of: printMatrix


    //--------------------------------------------------------------------

    /**
     * Method for printing an integer matrix  <P>
     * Based on ER Harold, "Java I/O", O'Reilly, around p 473. <br>
     * @param n1 row dimension of matrix
     * @param n2 column dimension of matrix
     * @param m input matrix values 
     * @param d display precision, number of decimal places
     * @param w display precision, total width of floating value
     */

    public static void printMatrix(int n1, int n2, int[][] m, int d, int w,
                 JTextArea outext)
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
                  outext.append(valString);
                  }
              // Start a new line at the end of a row
	      outext.append("\n");
          }
    }  // printMatrix


    //--------------------------------------------------------------------

    /**
     * Method printVect for printing a double float vector. <P>
     * Based on ER Harold, "Java I/O", O'Reilly, around p 473.  <br>
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
    }  // printVect


    //-------------------------------------------------------------------

    /**
     * Method printVect for printing an integer vector. <P>
     * Based on ER Harold, "Java I/O", O'Reilly, around p 473. <br>
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
    } // printVect

    //--------------------------------------------------------------------
    /**
     * Method getCoordInf for getting factor projections/loadings,
     * correlations and contributions. <br>
     * @param nrow
     * @param ncol
     * @param proj
     * @param corr
     * @param ctrn
     * @param data
     * @param evecs
     * @param evals
     * @param evals 
     * @param center
     * @param wts
     * @param cards
     */

    public static void getCoordInf(int nrow, int ncol, 
	   double[][] proj, double[][] corr, double[][] ctrn, 
	   double[][] data, double[][] evecs, 
	   double[] evals, double[] center, double[] wts, int[] cards)

    {
	// Falpha(c), projections of cluster c on factor alpha
	// or coordinate of c on the axis alpha.  In absolution value
	// Falpha(c) measures the distance which separates the centre 
	// of the cloud and the projection of the centre c onto the
	// axis alpha.  [Handbook, p. 609]
        for (int i = 0; i < nrow; i++) {
            for (int j1 = 0; j1 < ncol; j1++) {
                proj[i][j1] = 0.0; 
                for (int j2 = 0; j2 < ncol; j2++) 
                    proj[i][j1] += data[i][j2]*evecs[j2][j1];
	     }
        }
	// System.out.println("Cluster projections on factors. F #:");
	// printMatrix(nrow, ncol, proj, 4,4);
	// System.out.println();

	// CORalpha(c), sometimes denoted CO2alpha(c), squared cosine
        // of the angle made by the radius (O,c) with the axis alpha:
	// CORalpha(c) = (Falpha(c))^2 / ||O-c||^2
	// This is also the share of the factor alpha in the explanation
	// of the distance of c from the centre of gravity, and this 
	// quantity is necessarily less than or equal to 1, and decreases
	// from 1 to 0 as the angle made by (O,c) with the axis alpha
	// increases from 0 to 90 deg.  
        // [Correspondence Analysis Handbook, p. 609]
	double dstsq; 
	for (int i = 0; i < nrow; i++) {
	    dstsq = 0.0;
	    for (int j2 = 0; j2 < ncol; j2++)  //j2 ranges over variables
		// ||c-O||^2 is defined as: \sum_j 1/f_j 
		// ( q_kj/q_k - f_j)^2   where f_J = center
		dstsq += Math.pow((data[i][j2]-center[j2]), 2.0)/center[j2]; 
	    for (int j = 0; j < ncol; j++)   {  // j ranges over factors
		corr[i][j] = Math.pow(proj[i][j], 2.0)/dstsq;
	    }
	}
	// System.out.println("Cluster correlations, CO2:");
	// printMatrix(nrow, ncol, corr, 4, 4);

	// CTRalpha(c), the inertia of c wrt O = g projected onto the 
	// axis alpha, divided by the total inertia of the cloud projected
	// onto the axis alpha: 
	// CTRalpha(c) = (PDS(c).(Falpha(c))^2)/lambda_alpha
	for (int i = 0; i < nrow; i++) {
	    for (int j = 0; j < ncol; j++) {
		ctrn[i][j] = wts[i]*Math.pow(proj[i][j],2.0)/evals[j];
	    }
	}
	// System.out.println();
	// System.out.println("Cluster contributions, CTR:");
	// printMatrix(nrow, ncol, ctrn, 4, 4);
    } // getCoordInf


    //--------------------------------------------------------------------

    /**
     *  Method plotGraph 
     * @param n
     * @param m
     * @param obsx
     * @param obsy
     * @param vbex
     * @param vbey
     * @param rowlabs
     * @param collabs
     * @param axisx
     * @param axisy
     * @param tauxx
     * @param tauxy
     * @param outext
     */

    public static void plotGraph(int n, int m, double[] obsx, double[] obsy,
		  double[] vbex, double[] vbey, 
                  String[] rowlabs, String[] collabs, 
                  int axisx, int axisy, 
                  double tauxx, double tauxy, JTextArea outext)
    {

     // Set up display formating for floating values.
     NumberFormat nf = NumberFormat.getInstance();
     nf.setGroupingUsed(false);
     nf.setMaximumFractionDigits(4);
     nf.setMinimumFractionDigits(4);


     int maxLenRowlabs = 0;
     for (int i = 0; i < n; i++) {
	 if (rowlabs[i].length() > maxLenRowlabs) 
            maxLenRowlabs = rowlabs[i].length(); 
     }

     int maxLenCollabs = 0; 
     for (int j = 0; j < m; j++) {
	 if (collabs[j].length() > maxLenCollabs) 
            maxLenCollabs = collabs[j].length(); 
     }

     //     System.out.println("Max row and column label sizes are " +
     //			maxLenRowlabs + ", " + maxLenCollabs); 
 

     // Approach: define a 51x105 array initially of blanks, 
     // effective plot area is 50x100/104 (latter: 100 with label to 104).
     // Here we assume aspect ratios of 1 for x (factor 1) and 0.5 for y
     // (factor 2).

     // | chars at extremities, 103 other values.
     String plotlin0 = 
	 "|-------------------------------------------------------------------------------------------------------|";
     String plotline = 
	 "|                                                                                                       |";
     String[] plot = new String[105]; 
     for (int k = 0; k < 52; k++) {
     	 plot[k] = plotline; 
     } 

     double xymax, xymin, xmax, xmin, ymax, ymin;
     int xnull, ynull;
     xymin = MAXVAL; 
     xymax = SMALL;
     xmin = MAXVAL; 
     ymin = MAXVAL;
     xmax = SMALL; 
     ymax = SMALL; 
     for (int i = 0;  i < n; i++) {
	 if (obsx[i] < xymin) xymin = obsx[i];
	 if (obsx[i] > xymax) xymax = obsx[i];
	 if (obsy[i] < xymin) xymin = obsy[i];
	 if (obsy[i] > xymax) xymax = obsy[i];
	 if (obsx[i] < xmin) xmin = obsx[i];
	 if (obsx[i] > xmax) xmax = obsx[i];
	 if (obsy[i] < ymin) ymin = obsy[i];
	 if (obsy[i] > ymax) ymax = obsy[i]; 
     }
     for (int j = 0; j < m; j++) {
         if (vbex[j] < xymin) xymin = vbex[j];
         if (vbex[j] > xymax) xymax = vbex[j];
         if (vbey[j] < xymin) xymin = vbey[j];
         if (vbey[j] > xymax) xymax = vbey[j];
	 if (vbex[j] < xmin) xmin = vbex[j];
	 if (vbex[j] > xmax) xmax = vbex[j];
	 if (vbey[j] < ymin) ymin = vbey[j];
	 if (vbey[j] > ymax) ymax = vbey[j]; 
     }

     int[] intabs = new int[n];
     int[] intord = new int[n]; 
     for (int i = 0; i < n; i++) {
	 // min value = 1, max value = 100
	 // 4-char label at 100 => use locns 100, 101, 102, 103.
	 intabs[i] = (int)(99.0*(obsx[i]-xymin)/(xymax-xymin)+1.0);
	 // On ordinate, min value = 1, max value = 50
	 intord[i] = (int)(49.0*(obsy[i]-xymin)/(xymax-xymin)+1.0);
     }
     xnull = (int)(99.0*(0.0 - xymin)/(xymax-xymin)+1.0);
     ynull = (int)(49.0*(0.0 - xymin)/(xymax-xymin)+1.0);
     plot[ynull] = plotlin0;
     String plotrow, bit1, bit2, myString; 
     plotrow = plotline;   // Initialize
     for (int k = 0; k < 52; k++) {
	 plotrow = plot[k]; 
         plot[k] = plotrow.substring(0,xnull) + "|" + 
	     plotrow.substring(xnull+1, 105); 
     }


     // Now, assume rowlabs[i] is 4-char label.  New: < maxLenRowlabs chars.
     // Ordinate defines where to place vertically: plotrow = plot[intord[i]] 
     // Abscissa key location is intabs[i]
     // Take plotrow.substring(0, intabs[i]-1)  
     // + rowlabs[i]
     // + plotrow.substring(intabs[i]+4, 104);

     for (int i = 0; i < n; i++) {
	 myString = getSpaces(maxLenRowlabs - rowlabs[i].length()) + 
                              rowlabs[i]; 
	 plotrow = plot[intord[i]];
	 bit1 = plotrow.substring(0, intabs[i]);
	 bit2 = plotrow.substring(intabs[i] + maxLenRowlabs, 105);
	 plot[intord[i]] = bit1 + myString + bit2;
     }

     // Repeat for variables
     int[] colintabs = new int[m];
     int[] colintord = new int[m]; 
     // Note: scale using the same min and max projections as for rows.
     for (int j = 0; j < m; j++) {
	 // min value = 1, max value = 100
	 colintabs[j] = (int)(99.0*(vbex[j]-xymin)/(xymax-xymin)+1.0);
	 // On ordinate, min value = 1, max value = 50
	 colintord[j] = (int)(49.0*(vbey[j]-xymin)/(xymax-xymin)+1.0);
     }
     for (int j = 0; j < m; j++) {
	 myString = getSpaces(maxLenCollabs - collabs[j].length()) + 
                              collabs[j]; 
	 plotrow = plot[colintord[j]];
	 bit1 = plotrow.substring(0, colintabs[j]);
	 bit2 = plotrow.substring(colintabs[j]+maxLenCollabs, 105);
	 plot[colintord[j]] = bit1 + myString + bit2;
     }

     outext.append
                   ("Horizontal axis: " + axisx + "; min = " + 
		// ("Axe horizontal: " + axisx + "; min = " + 
		   nf.format(xmin) + "; max = " + nf.format(xmax) +
		   "; rate = " + nf.format(tauxx) + "\n");
                // "; taux = " + nf.format(tauxx) + "\n");
     outext.append
                   ("Vertical axis: " + axisy + "; min = " + 
		// ("Axe vertical: " + axisy + "; min = " + 
		   nf.format(ymin) + "; max = " + nf.format(ymax) +
		   "; rate = " + nf.format(tauxy) + "\n");
                //  "; taux = " + nf.format(tauxy) + "\n");
     outext.append
 //("A noter: des \u00e9l\u00e9ments peuvent etre non repr\u00e9sent\u00e9s ou partiellement cach\u00e9s.\n");
  ("Note: elements can be partially or fully overlapping.\n");
     outext.append(plotlin0 + "\n");
     for (int k = 0; k < 52; k++) {
	 outext.append(plot[k] + "\n");
     }
     outext.append(plotlin0 + "\n"); 


    } // plotGraph



    //--------------------------------------------------------------------

    /**
     *  Method preparePlotData
     * @param n
     * @param suppRowCount
     * @param m
     * @param suppColCount
     * @param rowproj
     * @param colproj
     * @param rowprojS
     * @param colprojS
     * @param obsx
     * @param obsy
     * @param rowlabs
     * @param collabs
     * @param suppRowlabs
     * @param suppCollabs
     * @param allrowlabs
     * @param allcollabs
     * @param axisx
     * @param axisy
     */

    public static void preparePlotData (int n, int suppRowCount, 
		     int m, int suppColCount,
		     double[][] rowproj, double[][] colproj, 
		     double[][] rowprojS, double[][] colprojS, 
		     double[] obsx, double[] obsy,
		     double[] vbex, double[] vbey,  
		     String[] rowlabs, String[] collabs, 
		     String[] suppRowlabs, String[] suppCollabs, 
		     String[] allrowlabs, String[] allcollabs, 
					int axisx, int axisy)
    {

         for (int i = 0;  i < n; i++) {     // observation projections
	     // Note: coordinate 0 is the trivial unit-valued factor,
	     // so we deal with factors numbered 1, 2, etc. 
	     obsx[i] = rowproj[i][axisx];
	     obsy[i] = rowproj[i][axisy];
	     allrowlabs[i] = rowlabs[i]; 
         }
         for (int i = 0; i < suppRowCount; i++) {  // suppl. obs. projns. 
	     obsx[n+i] = rowprojS[i][axisx]; 
	     obsy[n+i] = rowprojS[i][axisy];
	     allrowlabs[n+i] = suppRowlabs[i]; 
         }
         for (int j = 0;  j < m; j++) {   // variable projections
	     // Note: coordinate 0 is the trivial unit-valued factor.
	     vbex[j] = colproj[j][axisx];
	     vbey[j] = colproj[j][axisy];
	     allcollabs[j] = collabs[j]; 
         }
         for (int j = 0;  j < suppColCount; j++) { // suppl. vbe. projections
	     vbex[m+j] = colproj[j][axisx];
	     vbey[m+j] = colproj[j][axisy];
	     allcollabs[m+j] = suppCollabs[j]; 
         }

    } // End of method preparePlotData


    /** 
     * Method createWindow
     * @param outext     name of text area
     * @param outTitle   title for window header
     * @param fontSize   font size to use in the window (one size allowed)
     * @param close      allow close of window to terminate program?
     */

    public static void createWindow (JTextArea outext, String outTitle, 
				     int fontSize, boolean close) {

       // Output text window and layout as a JFrame
       outext.setFont(new Font("Monospaced", Font.PLAIN, fontSize));
       JScrollPane scroll = new JScrollPane(outext);
       JFrame frame = new JFrame(outTitle);
       if (close) frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
       frame.getContentPane().add(scroll);
       frame.setSize(400,200);
       frame.setJMenuBar(new ManageOutput(outext));        
       frame.setVisible(true);

    } // End of: createWindow


    /** 
     * Method getHeader
     * @param dataparams    header elements, object of class Header
     * @param bis           buffered reader
     * @param st            stream tokenizer
     * @param outext5       JTextArea text window
     */

    public static void getHeader ( Header dataparams,  
        BufferedReader bis, StreamTokenizer st, JTextArea outext5) {

	// Header information:
        // title, nrows, ncols, whether or not suppl. elements, 
        // no. factors, no. clusters for analysis purposes, followed 
        // by all column labels.


	try {            // Need to catch exceptions when reading data.

	    String title = null;          // Title
	    int norig = 0;                // No. rows of data  
	    int morig = 0;                // No. cols. of data
	    String suppIndicator = null;  // Whether or not suppl. elts.
	    int nfactors = 0;             // No. factors (for output)
	    int nclusters = 0;            // No. clusters (for interpretation)
	    String[] collabsOrig = null;  // Col. labels

            // Title is read until end-of-line control character is found.
            title = bis.readLine();

            // Row and column sizes, read in first
            st.nextToken();  norig = (int)st.nval;   // (original) no. rows
            st.nextToken();  morig = (int)st.nval;   // (original) no. cols

            // In addition to (original) nos. of rows, cols., we can have an 
            // indicator row or col. for, resp., suppl. cols. and suppl. rows.
            // Now look for indication of supplementary rows/columns.
            // Indication is the string: 
            // SUPR (suppl. rows), SUPC (suppl. cols.), SUPRC (both suppl. rows
            // and cols.), and by default SUPNO (no suppl. rows nor cols.).
            // Note: even with SUPR, SUPC or SUPRC, we will check for zero such
            // suppl. elements being present.
            switch (st.nextToken()) {
                  case StreamTokenizer.TT_NUMBER : 
	               System.out.println
	               ("We were expecting indicator of suppl. elements.");
	               System.out.println
                       ("This should follow nos. of rows, columns.");
	               System.out.println
                       ("Use one of: SUPR, SUPC, (default) SUPNO.");
	               System.exit(1);
	          case StreamTokenizer.TT_WORD : 
	               // Okay
            }
            suppIndicator = st.sval; 

            int n, m; 
            n = norig; 
            m = morig;
	    // Force upper case to allow for any use of lower case. 
            if ("SUPR".equals(suppIndicator.toUpperCase())) m++;
            if ("SUPC".equals(suppIndicator.toUpperCase())) n++; 

            // Column labels, which we first define and then populate.
            collabsOrig = new String[m]; 

            for (int j = 0; j < m; j++) {
	        // We now look for 
	        // EITHER nfactors AND nclusters AND set of column labels
	        // OR set of column labels only. 
	        switch (st.nextToken()) { 
	     	     case StreamTokenizer.TT_WORD :  
                           // We are reading the col. labels
		           collabsOrig[j] = st.sval;
                           if (nfactors == 0) nfactors = 7;  // set as default
	     		   break; 
	             case StreamTokenizer.TT_NUMBER : 
                           // We have found a value for nfactors
	                   if (j > 0) {
	     			  // Trouble! We have found > 1 numeric value
                            System.out.println("Expecting column labels.");
	     		       System.exit(1);
	     		   }
	                  nfactors = (int)st.nval; 
			  j = j-1;  // We still need to get all col labels
			  // Read a second numeric value, before proceeding
			  // to col. labels:
			  st.nextToken();
			  nclusters = (int)st.nval; 
	     		   break;
		}
	    }

            // Output some information.
            outext5.append
               ("Number of factors for outputing: " + nfactors + "\n");
            // ("Nombre de facteurs pour les sorties: " + nfactors + "\n");
            outext5.append
               ("Note: calculations are carried out in the \n");
            // ("Remarque: les calculs sont effectu\u00e9s dans \n");
            outext5.append
               ("factor space of full rank. \n"); 
            // ("l'espace factoriel de rang entier. \n"); 

            dataparams.title = title;
            dataparams.nrows = norig;
            dataparams.ncols = morig;
            dataparams.nfactors = nfactors;
            dataparams.nclusters = nclusters;
            dataparams.suppIndicator = suppIndicator;
            dataparams.columnIdentifiers = collabsOrig;

	}   // End of: try

        catch (java.io.IOException e) {
	    // Potentially thrown by: st.nextToken(); 
	    System.err.println("Error in getHeader: " + e); 
	    System.err.println
               ("Check for numeric/char problem in input data."); 
	}

   } // End of method: getHeader


    /** 
     * Method getArray
     * @param n         number of rows
     * @param m         number of columns
     * @param st        stream tokenizer
     * @returns         input data, as an object of class Imatrix
     */

    public static void getArray(int n, int m, StreamTokenizer st,
           double[][] indat, String[] rowlabs) { 

       double inval;

       // New read in input array values, successively
       // We include indicator row and/or col (i.e. indicating 
       // supplementary element) if they are present.
       // Warning: later we will allow only the case either of 
       // supplementary rows, or of supplementary columns. 

       try {
           for (int i = 0; i < n; i++) {
	       // New row, so first pick up row label (string) token.  
	       switch (st.nextToken()) {
 	          case StreamTokenizer.TT_NUMBER :
		  System.out.println("Expecting row label.");
                  System.exit(1);
                  break;
	       case StreamTokenizer.TT_WORD : 
                  rowlabs[i] = st.sval; // Token is a word, as wanted.
		  break; 
	       }
               for (int j = 0; j < m; j++) {
	          // Pick up # variables (cols.) numeric tokens.
	          // Examples of acceptable annotations following a value:
	          // Big value.  Uncertain  Estimate
	          // NOT acceptable is an exclamation mark: Big!  
                  switch (st.nextToken()) {
	             case StreamTokenizer.TT_NUMBER : inval = (double)st.nval;
                                              indat[i][j] = inval;
					      break;
	            case StreamTokenizer.TT_WORD :   System.out.println
					      ("Annotation: " + st.sval + 
					      " (location i,j: " + i + " " + 
                                              j + ")");
                                	      // Don't miss out next time:
		                              j = j - 1;
					      break;  
		  }
	       }
	   }
        } // End of: try

        catch (java.io.IOException e) {
	    // Potentially thrown by: st.nextToken(); 
	    System.err.println("Error in getArray: " + e); 
	    System.err.println
               ("Check for numeric/char problem in input data."); 
	}

   } // End of method: getArray



}  // End of DataAnalysis class


