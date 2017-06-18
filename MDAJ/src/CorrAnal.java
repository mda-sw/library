import java.io.*;
import java.util.*;
import java.text.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
// import com.macfaq.io.*;

/**
 * CorrAnal
 <P>
 DataAnalysisJava is a first version of a data analysis package for 
 Java.  The CorrAnal class carries out eigen-reduction, and determines
 factor rates of inertia, qualities, projections, correlations and 
 contributions, for main elements (rows, columns) and supplementary 
 elements.
 <P>
 @author Fionn Murtagh, fmurtagh@acm.org
 @version 0.0, 2003-March-15
 */

public class CorrAnal {

public static final double EPS = 1.0e-8;
public static final double SMALL = -1.0e10;
public static final double MAXVAL = 1.0e12;


    /* Diagonalization or eigenvector/eigenvalue decomposition
    @param principal    input matrix to be analyzed.
    @param CP           sums of squares and cross-products matrix to be 
                        diagonalized.
    @param Evals        eigenvalues.
    @param Evex         eigenvectors.
    @param rate         rate of inertia associated with eigenvectors.
    @param outext1      output window for display of results. 
    */

    double diagonalization (Imatrix principal, double[][] CP,  
      JTextArea outext1, double[] Evals, double[][] Evex, double[] rate)   {

	int nrow         = principal.getRowDimension();
	int ncol         = principal.getColumnDimension();
	String title     = principal.getTitle();
        String filename  = principal.getFilename();
	double[] rowmass = principal.getRowMass();
	double[] colmass = principal.getColumnMass();
	double[][] data  = principal.getArray(); 

        // Form matrix of cross-products to be analyzed, i.e. diagonalized.
        //double[][] CP = new double[m][m]; //cross-products, e.g. Burt table.
        for (int j1 = 0; j1 < ncol; j1++) {
	    for (int j2 = 0; j2 < ncol; j2++) {
                CP[j1][j2] = 0.0;
                for (int i = 0; i < nrow; i++) 
                    CP[j1][j2] += ( data[i][j1] * data[i][j2] ) /
                      ( rowmass[i] * Math.sqrt(colmass[j1]*colmass[j2]) ); 
	    }
        }

        // We will use Jama matrix class because it has the methods needed.  
        Matrix cp = new Matrix(CP);

        // Eigen decomposition
        EigenvalueDecomposition evaldec = cp.eig();
        Matrix evecs = evaldec.getV();
        double[] evals = evaldec.getRealEigenvalues();

        // Trace is adjusted by a value 1.0 because always in CA, 
        // the first eigenvalue is trivially 1-valued.
        double trce = cp.trace() - 1.0;

        // evecs contains the cols. ordered right to left.
        // Evecs is the more natural order with cols. ordered left to right.
        // So to repeat: leftmost col. of Evecs is assoc'd. with largest Evals.
        // Evals and Evecs ordered from left to right.
        double tot = 0.0; 
        for (int j = 0; j < evals.length; j++)  tot += evals[j]; 

        // Reverse order of evals into Evals.
        for (int j = 0; j < ncol; j++) Evals[j] = evals[ncol - j - 1];

        // Reverse order of Matrix evecs into Matrix Evecs.
        double[][] tempold = evecs.getArray();
        for (int j1 = 0; j1 < ncol; j1++) {
            for (int j2 = 0; j2 < ncol; j2++) 
                Evex[j1][j2] = tempold[j1][ncol - j2 - 1]/
                               Math.sqrt(colmass[j1]);
        }
        Matrix Evecs = new Matrix(Evex);
        // So as a Jama Matrix, we have Evecs, and 
        // as a multidim. array of dims ncol x ncol, we have Evex.

        double runningtotal = 0.0;
        double[] percentevals = new double[ncol];
        // Low index in following = 1 to exclude first trivial eval.
        percentevals[0] = 0.0; 
        rate[0] = 0.0;
        for (int j = 1; j < Evals.length; j++) {
            percentevals[j] = runningtotal + 100.0*Evals[j]/(tot-1.0);
	    rate[j] = Evals[j]/trce;
            runningtotal    = percentevals[j];
        }

	// Output printing follows.

        // First input data filename, and title.
        outext1.append (filename + "\n");
        outext1.append (title + "\n");

        // Some printing.  Next: trace, rank, lambda, taux, cumul.

        // Set up display formating for trace.  Floating 'trce' is 
        // properly formatted as 'num'.
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        nf.setMaximumFractionDigits(4);
        nf.setMinimumFractionDigits(4);
        String num = nf.format(trce);
        outext1.append("trace  :  " + num + "\n");
        //System.out.flush();

        outext1.append
           ("Lambda, rates and cumulative values are scaled by 10^-4 \n");
        // ("Les lambdas, taux et cumuls sont compt\u00e9s en 10^-4 \n");
        // Print from 2nd eval., since CA has trivial unitary 1st eval.
        int collo = 1;            // We will ignore trivial eigenvalue 1.
        int colhi = 7;            // With nrow > ncol, 
	                          // we have evals 0,1,2,...,m-1 
	// Note: above will be overruled in method DataAnalysis.printVect 
	// if we are dealing with a small no. (less than 7) of variables. 
        int abswidth = 6;         // Absolute field (or numeric token) spacing.
        // Ranks will be manually handled.
        outext1.append
            ("rank   :     1     2     3     4     5     6     7\n");
        //  ("rang   :     1     2     3     4     5     6     7\n");
        // Penultimate parameter is scaling.
        String lab;
        lab = "lambda :";
        DataAnalysis.printVect
           (Evals, collo, colhi, abswidth, 10000.0, lab, outext1);
        lab = "rates  :";
        // lab = "taux   :";
        DataAnalysis.printVect
           (rate, collo, colhi, abswidth, 10000.0, lab, outext1);
        lab = "cumul  :";
        DataAnalysis.printVect
                 (percentevals, collo, colhi, abswidth, 100.0, lab, outext1);


	// That completes the printing.

	return trce; 

    }   // end of: diagonalization method in class CorrAnal



    /** Projections on factors
    @param principal        matrix containing input data
    @param CP               sums of squares and cross-products matrix
    @param Evex             eigenvectors
    @param Evals            eigenvalues
    @param rowproj          row projections
    @param colproj          column projections
    */

    void projections (Imatrix principal, ResultSVD resEigenAnalysis, 
		      double[][] rowproj, double[][] colproj) {
	// double[][] CP, double[][] Evex, double[] Evals, 

	int n            = principal.getRowDimension();
	int m            = principal.getColumnDimension();
	double[] rowmass = principal.getRowMass();
	double[] colmass = principal.getColumnMass();
	double[][] data  = principal.getArray(); 
	double[][] Evex  = resEigenAnalysis.getEvectors();
	double[] Evals   = resEigenAnalysis.getEvalues();
	double[][] CP    = resEigenAnalysis.getCP();

        // Projections on factors - row, and column
        // Row projections in new space, X U  Dims: (n x m) x (m x m)
        for (int i = 0; i < n; i++) {
            for (int j1 = 0; j1 < m; j1++) {
                rowproj[i][j1] = 0.0; 
                for (int j2 = 0; j2 < m; j2++) 
                    rowproj[i][j1] += data[i][j2]*Evex[j2][j1];
                if (rowmass[i] >= EPS) rowproj[i][j1] /= rowmass[i];
                if (rowmass[i] < EPS)  rowproj[i][j1] = 0.0;
	    }
        }
        for (int j1 = 0; j1 < m; j1++) {
            for (int j2 = 0; j2 < m; j2++) {
                colproj[j1][j2] = 0.0;
                for (int j3 = 0; j3 < m; j3++) {
                    colproj[j1][j2] += CP[j1][j3] * Evex[j3][j2] *
                                       Math.sqrt(colmass[j3]);
	        }
                if (colmass[j1] >= EPS && Evals[j2] >= EPS) 
                   colproj[j1][j2] /= Math.sqrt(Evals[j2]*colmass[j1]);
                if (colmass[j1] < EPS && Evals[j2] < EPS) 
                   colproj[j1][j2] = 0.0; 
	    }
        }
    
    }  // End of projections method in CorrAnal class



    /** method projectionsSupp: determine projections on factors of 
    *   supplemenatary elements.
    @param m         dimensionality of original data space
    @param datSupp   data: on the basis of this, determine projections
    @param Evals     eigenvalues
    @param proj      full space row projections
    @param projSupp  projections to be determined
    */

    void projectionsSupp (int m, Imatrix IdatSupp, double[] Evals, 
                  double[][] proj, double[][] projSupp) {

        int nSupp        = IdatSupp.getRowDimension();
        int n            = IdatSupp.getColumnDimension();
        double[] sumSupp = IdatSupp.getRowMass();
	double[][] datSupp = IdatSupp.getArray();

        for (int i = 0; i < nSupp; i++) {
            for (int j1 = 0; j1 < m; j1++) {
                projSupp[i][j1] = 0.0; 
                for (int j2 = 0; j2 < n; j2++) {
                    projSupp[i][j1] += datSupp[i][j2]*proj[j2][j1];
	        }
	    }
        }
        for (int i = 0; i < nSupp; i++) {
            for (int j1 = 0; j1 < m; j1++) {
		if (sumSupp[i] >= EPS) projSupp[i][j1] /= sumSupp[i]; 
                if (sumSupp[i] < EPS) projSupp[i][j1] = 0.0;
		if (Evals[j1] >= EPS) 
                              projSupp[i][j1] /= Math.sqrt(Evals[j1]); 
		if (Evals[j1] < EPS) projSupp[i][j1] = 0.0; 
	    }	     
	}

    }  // End of projectionsSupp



    /** Row and column contributions to the factors.
    @param principal    input data 
    @param rowproj      row projectsions
    @param colproj      column projections
    @param rowcntr      row contributions
    @param colcntr      column contributions
    */

    void contributions (Imatrix principal,  
                        double[][] rowproj, double[][] colproj, 
                        double[][] rowcntr, double[][] colcntr) {

	int n            = principal.getRowDimension();
	int m            = principal.getColumnDimension();
	double[] rowmass = principal.getRowMass();
	double[] colmass = principal.getColumnMass();

        // Contributions to factors - row, and column
        double rowconColsum; 
        for (int j = 0; j < m; j++) {
            rowconColsum = 0.0;
            for (int i = 0; i < n; i++) {
                rowcntr[i][j] = rowmass[i]*Math.pow(rowproj[i][j], 2.0);
                rowconColsum += rowcntr[i][j];
	    }
            // Normalize so that sum of contributions for a factor equals 1
            for (int i = 0; i < n; i++) {
                if (rowconColsum > EPS) rowcntr[i][j] /= rowconColsum;
                if (rowconColsum <= EPS) rowcntr[i][j] = 0.0; 
	    }
        }
        double colconColsum; 
        for (int j1 = 0; j1 < m; j1++) {
            colconColsum = 0.0;
            for (int j2 = 0; j2 < m; j2++) {
                colcntr[j2][j1] = colmass[j2]*Math.pow(colproj[j2][j1], 2.0);  
                colconColsum +=  colcntr[j2][j1];
            }
	    // Normalize so that sum of contributions for a factor sum to 1
            for (int j2 = 0; j2 < m; j2++) {
                if (colconColsum > EPS) colcntr[j2][j1] /= colconColsum;
                if (colconColsum <= EPS) colcntr[j2][j1] = 0.0;
	    }
        }
    }  // End of contributions method in class CorrAnal



    /** Determine contributions of supplementary elements to factors
    @param dat          Imatrix containing data
    @param projSupp     projections 
    @param cntrSupp     contributions
    @param Evals        eigenvalues
    */

    void contributionsSupp (int m, Imatrix dat, double[][] projSupp, 
                        double[][] cntrSupp, 
                        double[] Evals) {
 
	int n            = dat.getRowDimension();
	double[] rowmass = dat.getRowMass();

	// See formula, p. 304 of Handbook
	 for (int j = 0; j < m; j++) {
	      for (int i = 0; i < n; i++) {
	          cntrSupp[i][j] = rowmass[i]*Math.pow(projSupp[i][j], 2.0);
	      }
	     for (int i = 0; i < n; i++) {
	         if (Evals[j] > EPS) cntrSupp[i][j] /= Evals[j];
	         if (Evals[j] <= EPS) cntrSupp[i][j] = 0.0; 
	     }
	 }

    }  // End of contributionsSupp method in class CorrAnal



    /** Correlations of rows and columns with factors.
    @param principal    input data
    @param rowproj      row projections
    @param colproj      column projections
    @param rowcorr      row correlations
    @param colcorr      column correlations
    */

    void correlations (Imatrix principal,
		       double[][] rowproj, double[][] colproj,
		       double[][] rowcorr, double[][] colcorr) {

	int n            = principal.getRowDimension();
	int m            = principal.getColumnDimension();
	double[] rowmass = principal.getRowMass();
	double[] colmass = principal.getColumnMass();
	double[][] data  = principal.getArray(); 

        // Correlations with factors - rows and column.
        // Return to this later: we may want to restrict to the 
	// calculation of 7 factors.
        // We're assuming here that n > mnew >= 7

        // First rows. 
        double distsq; 
        for (int i = 0; i < n; i++) {
	    distsq = 0.0;
            for (int j = 0; j < m; j++) {
                distsq += 
                     Math.pow((data[i][j]/rowmass[i] - colmass[j]),2.0) /
		     colmass[j]; 
	    }
	    for (int j = 0; j < m; j++) {
                rowcorr[i][j] = Math.pow(rowproj[i][j], 2.0)/distsq; 
	    }
        }

        // Now columns.
        // Return to this later: we may want to restrict to 7 factors.
        // We're assuming here that n > mnew >= 7
        for (int j1 = 0; j1 < m; j1++) {
	    distsq = 0.0;
            for (int i = 0; i < n; i++) {
                distsq += 
                 Math.pow((data[i][j1]/colmass[j1] - rowmass[i]),2.0) /
		 rowmass[i]; 
	    }
	    for (int j2 = 0; j2 < m; j2++) {
                colcorr[j1][j2] = Math.pow(colproj[j1][j2], 2.0)/distsq; 
	    }
        }
    }   // End of correlations method in class CorrAnal



    /** Correlations of supplementary elements on factors
    @param m        dimensionality of input data
    @param datSupp  Imatrix with data
    @param projSupp projections 
    @param corrSupp correlations to be determined
    */

    void correlationsSupp (int m, 
                       Imatrix datSupp,
		       double[][] projSupp, 
		       double[][] corrSupp) {

	int nSupp        = datSupp.getRowDimension();
	int n            = datSupp.getColumnDimension();
	double[] rowmass = datSupp.getRowMass();
	double[] colmass = datSupp.getColumnMass();
	double[][] data  = datSupp.getArray(); 



        double distsq; 
        for (int i = 0; i < nSupp; i++) {
	    distsq = 0.0;
            for (int j = 0; j < n; j++) {
                distsq += 
                     Math.pow((data[i][j]/rowmass[i] - 
                     colmass[j]),2.0) / colmass[j]; 
	    }
	    for (int j = 0; j < m; j++) {
                corrSupp[i][j] = Math.pow(projSupp[i][j], 2.0)/distsq; 
	    }
        }

    }   // End of correlationsSupp method in class CorrAnal




    /** Format and print output for: projections, correlations, contributions.
    @param principla    input data
    @param nfactors     number of factors wanted
    @param rowproj      row projections
    @param colproj      column projections
    @param cinertia     column inertias
    @param rowcntr      row contributions
    @param colcntr      column contributions
    @param rowcorr      row correlations
    @param colcorr      column correlations
    @param trce         trace
    @param outext1      output window
    */

    void formatprint (Imatrix principal, int nfactors, 
		      double[][] rowproj, double[][] colproj,
		      double[] cinertia,
		      double[][] rowcntr, double[][] colcntr,
		      double[][] rowcorr, double[][] colcorr,
		      double trce, JTextArea outext1) {

	int n            = principal.getRowDimension();
	int m            = principal.getColumnDimension();
	double[] rowmass = principal.getRowMass();
	double[] colmass = principal.getColumnMass();
	double[][] data  = principal.getArray(); 
	String[] rowlabs = principal.getRowIdentifier();
	String[] collabs = principal.getColumnIdentifier();
	String title = principal.getTitle();
	String filname = principal.getFilename();

	// Note: projections, correlations and contributions are
	// determined for the full-rank factor spaces.  However,
	// information on "nfactors" factors is output.

	// First determine inertias of observations (rows) and 
	// variables (columns), and qualities of representation 
	// in the new factor space.
        double[] rinertia = new double[n];       // relative inertia
        double[] rquality = new double[n];       // quality 
        double[] cquality = new double[m]; 
        for (int i = 0; i < n; i++) {
            rinertia[i] = 0.0; 
	    rquality[i] = 0.0; 
            for (int j = 0; j < m; j++) {
                rinertia[i] += 
                  Math.pow( (data[i][j]/rowmass[i]-colmass[j]),2.0 )/
                  colmass[j]; 
	    }
	    rinertia[i] = rowmass[i]*rinertia[i]/trce;
	    for (int k = 1; k <= nfactors; k++) {
	        rquality[i] += rowcorr[i][k]; 
	    }
        }

        for (int j = 0; j < m; j++) {
	    cinertia[j] = 0.0;
	    cquality[j] = 0.0; 
	    for (int i = 0; i < n; i++) {
	        cinertia[j] += 
                   Math.pow( (data[i][j]/colmass[j]-rowmass[i]),2.0 )/
                   rowmass[i];
	    }
	    cinertia[j] = colmass[j]*cinertia[j]/trce;
	    for (int k = 1; k <= nfactors; k++) {
	        cquality[j] += colcorr[j][k];
	    }
        }

        // Collect all row results.  
        int mbig = 3*nfactors + 3;   
        double[][] rowres = new double[n][mbig];
        for (int i = 0; i < n; i++) {
            rowres[i][0] = rquality[i];                 // QLT
	    //         rowres[i][1] = rowcntr[i][0];    // PDS
	    rowres[i][1] = rowmass[i];                  // PDS
            rowres[i][2] = rinertia[i];                 // INR 
            for (int k = 1; k < ((mbig-3)/3)+1; k++) {
                rowres[i][3*k] = rowproj[i][k];         // F#1, F#2, ...
                rowres[i][3*k+1] = rowcorr[i][k];       // CO2
                rowres[i][3*k+2] = rowcntr[i][k];       // CTR
	    }
        }
	// Define following for number formatting in labels
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        String myTitle = "|IDNI| QLT WTS INR|"; 
        // String myTitle = "|SIGI| QLT PDS INR|"; 
        for (int k = 1; k < ((mbig-3)/3)+1; k++) {
            if (nfactors < 10) {
	       // Check that sigle digit k is output as one digit in all cases.
               nf.setMaximumIntegerDigits(1);
               myTitle = myTitle + "  F#" + nf.format(k) +" CO2 CTR|"; 
	    }
            if (nfactors > 9) {
	       // Assumed is that k is never more than a 2-digit number.
               nf.setMaximumIntegerDigits(2);
               myTitle = myTitle + "  F" + nf.format(k) +" CO2 CTR|"; 
	    }
         }
	// Following, and in formatprintSupp, was: 5*mbig-3
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");
        outext1.append(myTitle + "\n");
        DataAnalysis.printAll(n, mbig, rowres, rowlabs, 1000.0, outext1); 
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");

        // Collect all column results.
        double[][] colres = new double[m][mbig];
        for (int j = 0; j < m; j++) {
            colres[j][0] = cquality[j];              // QLT
            colres[j][1] = colcntr[j][0];            // PDS
            colres[j][2] = cinertia[j];              // INR 
            for (int k = 1; k < ((mbig-3)/3)+1; k++) {
                colres[j][3*k] = colproj[j][k];      // F#1, F#2, ...
                colres[j][3*k+1] = colcorr[j][k];    // CO2
                colres[j][3*k+2] = colcntr[j][k];    // CTR
	    }
        }
        myTitle = "|IDNJ| QLT WTS INR|"; 
        // myTitle = "|SIGJ| QLT PDS INR|"; 
        for (int k = 1; k < ((mbig-3)/3)+1; k++) {
            if (nfactors < 10) {
	       // Check that sigle digit k is output as one digit in all cases.
               nf.setMaximumIntegerDigits(1);
               myTitle = myTitle + "  F#" + nf.format(k) +" CO2 CTR|"; 
	    }
            if (nfactors > 9) {
	       // Assumed is that k is never more than a 2-digit number.
               nf.setMaximumIntegerDigits(2);
               myTitle = myTitle + "  F" + nf.format(k) +" CO2 CTR|"; 
	    }
        }
        outext1.append(myTitle + "\n");
        DataAnalysis.printAll(m, mbig, colres, collabs, 1000.0, outext1); 
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");
    }  // End of method formatprint in class CorrAnal
 



    /** Format and print factor information for supplementary elements.
     * @param m
     * @param nfactors
     * @param datSupp
     * @param projSupp
     * @param cntrSupp
     * @param corrSupp
     * @param trce
     * @param outext1
    */

    void formatprintSupp (int m, int nfactors, 
                      Imatrix datSupp, 
		      double[][] projSupp, 
		      double[][] cntrSupp, 
		      double[][] corrSupp, 
		      double trce, 
		      JTextArea outext1) {

	int nSupp = datSupp.getRowDimension();
	int n     = datSupp.getColumnDimension();
	double[] sumSupp = datSupp.getRowMass();
	double[] sumC = datSupp.getColumnMass();
	String title = datSupp.getTitle();
	String filname = datSupp.getFilename();
	String[] labSupp = datSupp.getRowIdentifier();
	double[][] data = datSupp.getArray();

	// Note: projections, correlations and contributions are
	// determined for the full-rank factor spaces.  However,
	// information on "nfactors" factors is output.

	// First determine inertias 
	// and qualities of representation 
	// in the new factor space.
        double[] inertia = new double[nSupp];       // relative inertia
        double[] quality = new double[nSupp];       // quality 

        for (int i = 0; i < nSupp; i++) {
            inertia[i] = 0.0; 
	    quality[i] = 0.0; 
            for (int j = 0; j < n; j++) {
                inertia[i] += 
                  Math.pow( (data[i][j]/sumSupp[i]-sumC[j]),2.0 )/
                  sumC[j]; 
	    }
	    inertia[i] = sumSupp[i]*inertia[i]/trce;
	    for (int k = 1; k <= nfactors; k++) {
	        quality[i] += corrSupp[i][k]; 
	    }
        }

        // =================================================================
        // Collect all row results.
        int mbig = 3*nfactors + 3;   
        double[][] res = new double[nSupp][mbig];
        for (int i = 0; i < nSupp; i++) {
            res[i][0] = quality[i];                   // QLT
	    res[i][1] = sumSupp[i];                   // PDS or WTS
            res[i][2] = inertia[i];                   // INR 
            for (int k = 1; k < ((mbig-3)/3)+1; k++) {
                res[i][3*k] = projSupp[i][k];         // F#1, F#2, ...
                res[i][3*k+1] = corrSupp[i][k];       // CO2
                res[i][3*k+2] = cntrSupp[i][k];       // CTR
	    }
        }
	// Define following for number formatting in labels
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        String myTitle = "|IDNI| QLT WTS INR|"; 
        // String myTitle = "|SIGI| QLT PDS INR|"; 
        for (int k = 1; k < ((mbig-3)/3)+1; k++) {
            if (nfactors < 10) {
	       // Check that sigle digit k is output as one digit in all cases.
               nf.setMaximumIntegerDigits(1);
               myTitle = myTitle + "  F#" + nf.format(k) +" CO2 CTR|"; 
	    }
            if (nfactors > 9) {
	       // Assumed is that k is never more than a 2-digit number.
               nf.setMaximumIntegerDigits(2);
               myTitle = myTitle + "  F" + nf.format(k) +" CO2 CTR|"; 
	    }
         }
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");
        outext1.append(myTitle + "\n");
        DataAnalysis.printAll(nSupp, mbig, res, labSupp, 1000.0, outext1); 
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");

    }  // End of method formatprintSupp in class CorrAnal
 


    //======================================================================

    // Following are deprecated overloaded versions of the foregoing methods.

    // ---------------------------------------------------------------------
    double diagonalization (int n, int mnew,  
                          double[][] indatstd, double[][] CP, 
                          JTextArea outext1,
			  double[] rowsums, double[] colsums,  
			  double[] Evals, double[][] Evex, double[] rate,
                          String title, String filname, 
                          boolean commandline)
    {

        // Form matrix to be analyzed (diagonalized)
        //double[][] CP = new double[mnew][mnew]; //x-products, e.g. Burt table
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

        // Use Jama matrix class.  
        Matrix cp = new Matrix(CP);
        // For information: print out cp matrix.
        // System.out.println(" Matrix to be diagonalized");
        // cp.print(6,2);

        //-------------------------------------------------------------------
        // Eigen decomposition
        EigenvalueDecomposition evaldec = cp.eig();
        Matrix evecs = evaldec.getV();
        double[] evals = evaldec.getRealEigenvalues();
        // Trace is adjusted by a value 1.0 because always in CA, 
        // the first eigenvalue is trivially 1-valued.
        double trce = cp.trace() - 1.0;

        // evecs contains the cols. ordered right to left
        // Evecs is the more natural order with cols. ordered left to right
        // So to repeat: leftmost col. of Evecs is assoc with largest Evals
        // Evals and Evecs ordered from left to right

        double tot = 0.0; 
        for (int j = 0; j < evals.length; j++)  {
            tot += evals[j]; 
        }
        // reverse order of evals into Evals
        // double[] Evals = new double[mnew];  // passed as parameter 
        for (int j = 0; j < mnew; j++) {
            Evals[j] = evals[mnew - j - 1];
        }
        // reverse order of Matrix evecs into Matrix Evecs
        double[][] tempold = evecs.getArray();
        // double[][] Evex = new double[mnew][mnew];  // Passed as parameter
        for (int j1 = 0; j1 < mnew; j1++) {
            for (int j2 = 0; j2 < mnew; j2++) {
                Evex[j1][j2] = tempold[j1][mnew - j2 - 1] / 
                               Math.sqrt(colsums[j1]);
	    }
        }
        Matrix Evecs = new Matrix(Evex);
        //     Evecs.print(10,4);
        // So as a JAMA Matrix, we have Evecs, and 
        // as a multidim. array of dims mnew x mnew, we have Evex

        double runningtotal = 0.0;
        double[] percentevals = new double[mnew];
        // double[] rate = new double[mnew];  // Passed as parameter
        // low index in following = 1 to exclude first trivial eval.
        percentevals[0] = 0.0; 
        rate[0] = 0.0;
        for (int j = 1; j < Evals.length; j++) {
            percentevals[j] = runningtotal + 100.0*Evals[j]/(tot-1.0);
	    rate[j] = Evals[j]/trce;
            runningtotal    = percentevals[j];
        }

        // ==================================================================
	// Output printing follows.
	// Two possibilities - to the command prompt window, if 
	// the boolean variable commandline is true; and in any case
	// to a JTextArea panel.

        // Some printing.  First input data filename, and title.
        if (commandline) System.out.println (filname);
        if (commandline) System.out.println(title);
        outext1.append (filname + "\n");
        outext1.append (title + "\n");
        // Some printing.  Next: trace, rank, lambda, taux, cumul.

        // Set up display formating for trace.  Floating 'trce' is 
        // properly format as 'num'.
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        nf.setMaximumFractionDigits(4);
        nf.setMinimumFractionDigits(4);
        String num = nf.format(trce);
        if (commandline) System.out.print("trace  :  " + num);
        outext1.append("trace  :  " + num + "\n");
        System.out.flush();
        if (commandline) System.out.println();

        if (commandline) 
        System.out.println("Les lambdas, taux et cumuls sont comptes en e-4");
        outext1.append
           ("Les lambdas, taux et cumuls sont compt\u00e9s en 10^-4 \n");
        // Print from 2nd eval., since CA has trivial unitary 1st eval.
        int collo = 1;  // We will ignore trivial eigenvalue 1.
        int colhi = 7;  // With n > mnew, we have evals 0,1,2,...,mnew-1 
	// Note: above will be overruled in method DataAnalysis.printVect 
	// if we are dealing with a small no. (less than 7) of variables. 
        int abswidth = 6;       // Absolute field (or numeric token) spacing.
        // Ranks will be manually handled.
        if (commandline) System.out.println
                   ("rang   :     1     2     3     4     5     6     7");
        outext1.append("rang   :     1     2     3     4     5     6     7\n");
        // Penultimate parameter is scaling.
        String lab;
        lab = "lambda :";
        DataAnalysis.printVect
           (Evals, collo, colhi, abswidth, 10000.0, lab, outext1);
        lab = "taux   :";
        DataAnalysis.printVect
           (rate, collo, colhi, abswidth, 10000.0, lab, outext1);
        lab = "cumul  :";
        DataAnalysis.printVect
                 (percentevals, collo, colhi, abswidth, 100.0, lab, outext1);

	// That completes the printing.
        // ==================================================================

	return trce; 
    }   // end of: diagonalization method in class CorrAnal
    // ---------------------------------------------------------------------


    // ---------------------------------------------------------------------
    void projections (int n, int mnew, double[][] indatstd,
		      double[][] CP,  
		      double[][] Evex, double[] Evals, 
		      double[][] rowproj, double[][] colproj) {

        double[] rowsums = new double[n]; 
        double[] colsums = new double[mnew]; 

	for (int i = 0; i < n; i++) {
	    rowsums[i] = 0.0; 
	    for (int j = 0; j < mnew; j++) {
		rowsums[i] += indatstd[i][j]; 
	    }
	}
	for (int j = 0; j < mnew; j++) {
	    colsums[j] = 0.0; 
	    for (int i = 0; i < n; i++) {
		colsums[j] += indatstd[i][j]; 
	    }
	}

        //-------------------------------------------------------------------
        // Projections on factors - row, and column
        // Row projections in new space, X U  Dims: (n x m) x (m x m)
        System.out.println();
        // double[][] rowproj = new double[n][mnew];  // Passed as parameter
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
        // double[][] colproj = new double[mnew][mnew];  // Passed as parameter
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
    
    }  // End of projections method in CorrAnal class
   // ---------------------------------------------------------------------


    // ---------------------------------------------------------------------
    // Alternative method for row projections, yielding 
    // identical output to method projectionsSupp
    void projectionsSuppx (int nSupp, int m, double[][] datSupp,
		      double[][] Evex, double[] sumSupp, 
                      double[][] projSupp) {

        for (int i = 0; i < nSupp; i++) {
            for (int j1 = 0; j1 < m; j1++) {
                projSupp[i][j1] = 0.0; 
                for (int j2 = 0; j2 < m; j2++) {
                    projSupp[i][j1] += datSupp[i][j2]*Evex[j2][j1];
	        }
                if (sumSupp[i] >= EPS) projSupp[i][j1] /= 
                    sumSupp[i];
                if (sumSupp[i] < EPS) projSupp[i][j1] = 0.0;
	    }
        }

    }  // End of projectionsSuppx


    // ---------------------------------------------------------------------
    void projectionsSupp (int nSupp, int n, int m, double[][] datSupp,
		      double[] Evals, double[][] proj, double[] sumSupp, 
                      double[][] projSupp) {

        for (int i = 0; i < nSupp; i++) {
            for (int j1 = 0; j1 < m; j1++) {
                projSupp[i][j1] = 0.0; 
                for (int j2 = 0; j2 < n; j2++) {
                    projSupp[i][j1] += datSupp[i][j2]*proj[j2][j1];
	        }
	    }
        }
        for (int i = 0; i < nSupp; i++) {
            for (int j1 = 0; j1 < m; j1++) {
		if (sumSupp[i] >= EPS) projSupp[i][j1] /= sumSupp[i]; 
                if (sumSupp[i] < EPS) projSupp[i][j1] = 0.0;
		if (Evals[j1] >= EPS) 
                              projSupp[i][j1] /= Math.sqrt(Evals[j1]); 
		if (Evals[j1] < EPS) projSupp[i][j1] = 0.0; 
	    }	     
	}

    }  // End of projectionsSupp
    // ---------------------------------------------------------------------


    // ---------------------------------------------------------------------
    void contributions (int n, int mnew, 
                        double[][] rowproj, double[][] colproj, 
                        double[] rowsums, double[] colsums, 
                        double[][] rowcntr, double[][] colcntr) {

        // Contributions to factors - row, and column
        // double[][] rowcntr = new double[n][mnew];  // Passed as parameter
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
        // double[][] colcntr = new double[mnew][mnew];  // Passed as parameter
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
    }  // End of contributions method in class CorrAnal
    // ---------------------------------------------------------------------


    // ---------------------------------------------------------------------
    void contributionsSupp (int nSupp, int m, double[][] projSupp, 
                        double[] sumSupp, double[][] cntrSupp, 
                        double[] Evals) {
 
	// See formula, p. 304 of Handbook
	 for (int j = 0; j < m; j++) {
	      for (int i = 0; i < nSupp; i++) {
	          cntrSupp[i][j] = sumSupp[i]*Math.pow(projSupp[i][j], 2.0);
	      }
	     for (int i = 0; i < nSupp; i++) {
	         if (Evals[j] > EPS) cntrSupp[i][j] /= Evals[j];
	         if (Evals[j] <= EPS) cntrSupp[i][j] = 0.0; 
	     }
	 }

    }  // End of contributionsSupp method in class CorrAnal
    // ---------------------------------------------------------------------


    // ---------------------------------------------------------------------
    void correlations (int n, int mnew, double[][] indatstd,
		       double[] rowsums, double[] colsums,
		       double[][] rowproj, double[][] colproj,
		       double[][] rowcorr, double[][] colcorr) {

        // Correlations with factors - rows and column.
        // Return to this later: we may want to restrict to the 
	// calculation of 7 factors.
        // We're assuming here that n > mnew >= 7
        // First rows. 
        // double[][] rowcorr = new double[n][mnew];   // Passed as parameter
        double distsq; 
        for (int i = 0; i < n; i++) {
	    distsq = 0.0;
            for (int j = 0; j < mnew; j++) {
                distsq += 
                     Math.pow((indatstd[i][j]/rowsums[i] - colsums[j]),2.0) /
		     colsums[j]; 
	    }
	    for (int j = 0; j < mnew; j++) {
                rowcorr[i][j] = Math.pow(rowproj[i][j], 2.0)/distsq; 
	    }
        }
        // Now columns.
        // Return to this later: we may want to restrict to 7 factors.
        // We're assuming here that n > mnew >= 7
        // double[][] colcorr = new double[mnew][mnew]; // Passed as parameter
        for (int j1 = 0; j1 < mnew; j1++) {
	    distsq = 0.0;
            for (int i = 0; i < n; i++) {
                distsq += 
                 Math.pow((indatstd[i][j1]/colsums[j1] - rowsums[i]),2.0) /
		 rowsums[i]; 
	    }
	    for (int j2 = 0; j2 < mnew; j2++) {
                colcorr[j1][j2] = Math.pow(colproj[j1][j2], 2.0)/distsq; 
	    }
        }
    }   // End of correlations method in class CorrAnal
    // ---------------------------------------------------------------------


    // ---------------------------------------------------------------------
    void correlationsSupp (int nSupp, int n, int m, 
                       double[][] datSupp,
		       double[] sumSupp, 
		       double[][] projSupp, 
		       double[][] corrSupp, double[] sumC) {

        double distsq; 
        for (int i = 0; i < nSupp; i++) {
	    distsq = 0.0;
            for (int j = 0; j < n; j++) {
                distsq += 
                     Math.pow((datSupp[i][j]/sumSupp[i] - 
                     sumC[j]),2.0) / sumC[j]; 
	    }
	    for (int j = 0; j < m; j++) {
                corrSupp[i][j] = Math.pow(projSupp[i][j], 2.0)/distsq; 
	    }
        }

    }   // End of correlationsSupp method in class CorrAnal
    // --------------------------------------------------------------------- 

    // ---------------------------------------------------------------------
    void formatprint (int n, int mnew, int nfactors, double[][] indatstd, 
		      double[][] rowproj, double[][] colproj,
		      double[] rowsums, double[] colsums, double[] cinertia,
		      double[][] rowcntr, double[][] colcntr,
		      double[][] rowcorr, double[][] colcorr,
		      String[] rowlabs, String[] collabs, double trce, 
		      String title, String filname,
		      JTextArea outext1, boolean commandline) {

	// Note: projections, correlations and contributions are
	// determined for the full-rank factor spaces.  However,
	// information on "nfactors" factors is output.

	// First determine inertias of observations (rows) and 
	// variables (columns), and qualities of representation 
	// in the new factor space.
        double[] rinertia = new double[n];       // relative inertia
        // double[] cinertia = new double[mnew];    // relative inertia
        double[] rquality = new double[n];       // quality 
        double[] cquality = new double[mnew]; 
        for (int i = 0; i < n; i++) {
            rinertia[i] = 0.0; 
	    rquality[i] = 0.0; 
            for (int j = 0; j < mnew; j++) {
                rinertia[i] += 
                  Math.pow( (indatstd[i][j]/rowsums[i]-colsums[j]),2.0 )/
                  colsums[j]; 
	    }
	    rinertia[i] = rowsums[i]*rinertia[i]/trce;
	    for (int k = 1; k <= nfactors; k++) {
	        rquality[i] += rowcorr[i][k]; 
	    }
        }

        for (int j = 0; j < mnew; j++) {
	    cinertia[j] = 0.0;
	    cquality[j] = 0.0; 
	    for (int i = 0; i < n; i++) {
	        cinertia[j] += 
                   Math.pow( (indatstd[i][j]/colsums[j]-rowsums[i]),2.0 )/
                   rowsums[i];
	    }
	    cinertia[j] = colsums[j]*cinertia[j]/trce;
	    for (int k = 1; k <= nfactors; k++) {
	        cquality[j] += colcorr[j][k];
	    }
        }

        // =================================================================
        // Collect all row results.
        int mbig = 3*nfactors + 3;   
        double[][] rowres = new double[n][mbig];
        for (int i = 0; i < n; i++) {
            rowres[i][0] = rquality[i];                 // QLT
	    //         rowres[i][1] = rowcntr[i][0];    // PDS
	    rowres[i][1] = rowsums[i];                  // PDS
            rowres[i][2] = rinertia[i];                 // INR 
            for (int k = 1; k < ((mbig-3)/3)+1; k++) {
                rowres[i][3*k] = rowproj[i][k];         // F#1, F#2, ...
                rowres[i][3*k+1] = rowcorr[i][k];       // CO2
                rowres[i][3*k+2] = rowcntr[i][k];       // CTR
	    }
        }
	// Define following for number formatting in labels
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        String myTitle = "|SIGI| QLT PDS INR|"; 
        for (int k = 1; k < ((mbig-3)/3)+1; k++) {
            if (nfactors < 10) {
	       // Check that sigle digit k is output as one digit in all cases.
               nf.setMaximumIntegerDigits(1);
               myTitle = myTitle + "  F#" + nf.format(k) +" CO2 CTR|"; 
	    }
            if (nfactors > 9) {
	       // Assumed is that k is never more than a 2-digit number.
               nf.setMaximumIntegerDigits(2);
               myTitle = myTitle + "  F" + nf.format(k) +" CO2 CTR|"; 
	    }
         }
	// Following, and in formatprintSupp, was: 5*mbig-3
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");
        if (commandline) System.out.println(myTitle);
        outext1.append(myTitle + "\n");
        DataAnalysis.printAll(n, mbig, rowres, rowlabs, 1000.0, outext1); 
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");

        // Collect all column results.
        double[][] colres = new double[mnew][mbig];
        for (int j = 0; j < mnew; j++) {
            colres[j][0] = cquality[j];              // QLT
            colres[j][1] = colcntr[j][0];            // PDS
            colres[j][2] = cinertia[j];              // INR 
            for (int k = 1; k < ((mbig-3)/3)+1; k++) {
                colres[j][3*k] = colproj[j][k];      // F#1, F#2, ...
                colres[j][3*k+1] = colcorr[j][k];    // CO2
                colres[j][3*k+2] = colcntr[j][k];    // CTR
	    }
        }
        myTitle = "|SIGJ| QLT PDS INR|"; 
        for (int k = 1; k < ((mbig-3)/3)+1; k++) {
            if (nfactors < 10) {
	       // Check that sigle digit k is output as one digit in all cases.
               nf.setMaximumIntegerDigits(1);
               myTitle = myTitle + "  F#" + nf.format(k) +" CO2 CTR|"; 
	    }
            if (nfactors > 9) {
	       // Assumed is that k is never more than a 2-digit number.
               nf.setMaximumIntegerDigits(2);
               myTitle = myTitle + "  F" + nf.format(k) +" CO2 CTR|"; 
	    }
        }
        if (commandline) System.out.println(myTitle);
        outext1.append(myTitle + "\n");
        DataAnalysis.printAll(mnew, mbig, colres, collabs, 1000.0, outext1); 
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");
    }  // End of method formatprint in class CorrAnal
    // ---------------------------------------------------------------------

    // ---------------------------------------------------------------------
    void formatprintSupp (int nSupp, int n, int m, int nfactors, 
                      double[][] datSupp, 
		      double[][] projSupp, 
		      double[] sumSupp, 
                      double[] sumC, 
		      double[][] cntrSupp, 
		      double[][] corrSupp, 
		      String[] labSupp, double trce, 
		      String title, String filname,
		      JTextArea outext1, boolean commandline) {

	// Note: projections, correlations and contributions are
	// determined for the full-rank factor spaces.  However,
	// information on "nfactors" factors is output.

	// First determine inertias 
	// and qualities of representation 
	// in the new factor space.
        double[] inertia = new double[nSupp];       // relative inertia
        double[] quality = new double[nSupp];       // quality 

        for (int i = 0; i < nSupp; i++) {
            inertia[i] = 0.0; 
	    quality[i] = 0.0; 
            for (int j = 0; j < n; j++) {
                inertia[i] += 
                  Math.pow( (datSupp[i][j]/sumSupp[i]-sumC[j]),2.0 )/
                  sumC[j]; 
	    }
	    inertia[i] = sumSupp[i]*inertia[i]/trce;
	    for (int k = 1; k <= nfactors; k++) {
	        quality[i] += corrSupp[i][k]; 
	    }
        }

	// double tot = 0.0; 
	// for (int i = 0; i < nSupp; i++) {
	//    tot += sumSupp[i]; 
	// }
	// for (int i = 0; i < nSupp; i++) {
	//    sumSupp[i] /= tot; 
	// }


        // =================================================================
        // Collect all row results.
        int mbig = 3*nfactors + 3;   
        double[][] res = new double[nSupp][mbig];
        for (int i = 0; i < nSupp; i++) {
            res[i][0] = quality[i];                   // QLT
	    res[i][1] = sumSupp[i];                   // PDS
            res[i][2] = inertia[i];                   // INR 
            for (int k = 1; k < ((mbig-3)/3)+1; k++) {
                res[i][3*k] = projSupp[i][k];         // F#1, F#2, ...
                res[i][3*k+1] = corrSupp[i][k];       // CO2
                res[i][3*k+2] = cntrSupp[i][k];       // CTR
	    }
        }
	// Define following for number formatting in labels
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        String myTitle = "|SIGI| QLT PDS INR|"; 
        for (int k = 1; k < ((mbig-3)/3)+1; k++) {
            if (nfactors < 10) {
	       // Check that sigle digit k is output as one digit in all cases.
               nf.setMaximumIntegerDigits(1);
               myTitle = myTitle + "  F#" + nf.format(k) +" CO2 CTR|"; 
	    }
            if (nfactors > 9) {
	       // Assumed is that k is never more than a 2-digit number.
               nf.setMaximumIntegerDigits(2);
               myTitle = myTitle + "  F" + nf.format(k) +" CO2 CTR|"; 
	    }
         }
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");
        if (commandline) System.out.println(myTitle);
        outext1.append(myTitle + "\n");
        DataAnalysis.printAll(nSupp, mbig, res, labSupp, 1000.0, outext1); 
        for (int j = 0; j < (5*mbig-1); j++) outext1.append("_");
        outext1.append("\n");

    }  // End of method formatprintSupp in class CorrAnal
    // ---------------------------------------------------------------------


}  // End of class CorrAnal


