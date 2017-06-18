import java.io.*;
import java.util.*;
import java.text.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
// import com.macfaq.io.*;

/**
 * HierClass
 <P>
 DataAnalysisJava is a first version of a data analysis package for 
 Java.  The HierClass class carries out hierarchical clustering using
 a reciprocal nearest neighbors algorithm.  The minimal variance 
 agglomerative criterion is used. 
 @author Fionn Murtagh, fmurtagh@acm.org
 @version 0.0, 2003-March-15
 */

public class HierClass {

public static final double EPS = 1.0e-8;
public static final double SMALL = -1.0e10;
public static final double MAXVAL = 1.0e12;

    void reciprocalNearestNeighbors (int nrow, int ncol, 
          double[][] indat2clstr, 
          double[] mass, JTextArea outext2,
	  int[][] clusters, int[][] clusters2, int[] Acl, int[] Bcl, 
          double[] levels, double[] card) 
        {

		  // Algorithm:
		  // - get all pairwise dissimilarities
		  // - get nearest neighbors and NN dissimilarities
		  // - loop to carry our series of agglomerations
		  // - do this by finding smallest dissimilarity
		  // - then update nearest neighbor list and NN dissims.
		  // - we also update the cluster label matrix,
		  //   following each agglomeration
 
        // Some definitions and initializations
	// int[][] clusters = new int[nrow][nrow];  // cluster label matrix
	// int[][] clusters2 = new int[nrow][nrow]; 
        int[] nn = new int[nrow];                // cluster nearest neigh's
        int[] flag = new int[nrow];              // flag of active obs
        double[] nndiss = new double[nrow];      // nearest neigh diss's
        double[] clcard = new double[nrow];      // cluster cardinality
	// A, B convention: cluster has seqno. = min. subcluster seqno.
	int[] A = new int[nrow-1];               // aine, elder
	int[] B = new int[nrow-1];               // benjamin, younger
	// Following are passed as parameters:
	// double[] levels = new double[nrow-1]; // cluster levels 
	// double[] card = new double[nrow-1];   // cluster cardinalities
	// Acl, Bcl convention: cluster has seqno. = nrow+1, ... 2*nrow-1
	// int[] Acl = new int[nrow-1];             // aine, elder
	// int[] Bcl = new int[nrow-1];             // benjamin, younger
        int ncl;   
        ncl = nrow;                              // Initial no. of clusters  

        for (int i = 0; i < nrow-1; i++) {       // nrow-1 elements
            flag[i] = 1;
            clcard[i] = 1.0; 
	    levels[i] = 0.0;
	    card[i] = 0.0;
	    A[i] = 0;
	    B[i] = 0;
	    Acl[i] = 0;
	    Bcl[i] = 0;
        }
	flag[nrow-1] = 1;                        // nrow'th value
	clcard[nrow-1] = 1.0;

        // Get dissimilarities to start the analysis on.
        double[][] diss = new double[nrow][nrow];
        diss = Dissim(nrow, ncol, mass, indat2clstr);  // Dissim method
	// Get nearest neighbors and nearest neighbor dissimilarities.
        getNNs(nrow, flag, diss, nn, nndiss);          // getNNs method
	//  System.out.println("Nearest neighbors:");
	//  DataAnalysis.printVect(nn, 4, 10); 
        //  System.out.println("Nearest neighbors dissimilarities:");
        //  DataAnalysis.printVect(nndiss, 4, 10);

        // Get closest neighbors, using nndiss, followed by nn.
        int clust1 = 0;
        int clust2 = 0;
        int cl1    = 0;
        int cl2    = 0;
	// Call to initialize
        ClustMat(nrow, clusters, clusters2, clust1, clust2, ncl, Acl, Bcl);
        //  System.out.println("No. of clusters:  " + ncl);  
	//  System.out.println(" Cluster label matrix:");
        //  DataAnalysis.printMatrix(nrow, nrow, clusters, 4, 10);

        int minobs;
        double mindist;  
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
	    }  // end of i-loop 
	    // minobs is one cluster label, the other is nn[minobs].
	    // Adjust for 0-sequencing: indices have to be 0-sequenced,
	    // but for our convenience we want cluster numbers to be 
	    // 1-sequenced.
  	    if (minobs < nn[minobs]) {
               clust1 = minobs + 1;
               clust2 = nn[minobs];
            }
  	    if (minobs > nn[minobs]) {
               clust2 = minobs + 1;
               clust1 = nn[minobs];
            }
	    // Now we have clust1 < clust2, and we'll agglomerate.
	    // System.out.println(" clus#1: "  + clust1 +
	    //                    ";  clus#2: "  + clust2 +
            //                    ";  new card: " + (clcard[clust1-1]+
            //                                       clcard[clust2-1]) + 
            //                    "; # clus left: " + ncl +
            //                    "; mindiss: " + mindist);
	    // We will carry out the agglomeration.
            ncl = ncl - 1; 
	    levels[ncl-1] = mindist; 
	    card[ncl-1] = (double)clcard[clust1-1] + (double)clcard[clust2-1];
	    A[ncl-1] = clust1;
	    B[ncl-1] = clust2;

	    // Update cluster label matrix following the agglomeration.
            ClustMat(nrow, clusters, clusters2, clust1, clust2, ncl, Acl,Bcl);
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
		if ( (i != cl1) && (i != cl2) && (flag[i] == 1) ) {
		// Minimum variance criterion 
		// (This is the only criterion used here, since it is 
		// coherent with the Euclidean factor space metric used.)
		    diss[cl1][i] = 
		    ((mass[cl1]+mass[i])/(mass[cl1]+mass[cl2]+mass[i]))
		                *diss[cl1][i]   +
		    ((mass[cl2]+mass[i])/(mass[cl1]+mass[cl2]+mass[i]))
		                *diss[cl2][i]   -
		    ((mass[i])/(mass[cl1]+mass[cl2]+mass[i]))
		                *diss[cl1][cl2] ;
                    diss[i][cl1] = diss[cl1][i];
		}
	    }
	    // Updates of cardinality and mass:
            clcard[cl1] = clcard[cl1] + clcard[cl2];
	    mass[cl1] = mass[cl1] + mass[cl2];  
	    // Cluster label clust2 is knocked out:
            flag[cl2] = 0; 
            nndiss[cl2] = MAXVAL; 
	    mass[cl2] = 0.0;  clcard[cl2] = 0.0;
            for (int i = 0; i < nrow; i++) {
                diss[cl2][i] = MAXVAL;
                diss[i][cl2] = diss[cl2][i];
	    }

 	    // Update nearest neighbors and nearest neighbor dissimilarities.
            getNNs(nrow, flag, diss, nn, nndiss);
	    // System.out.println("Nearest neighbors of items 1, 2, ... n:");
	    // DataAnalysis.printVect(nn, 4, 10); 
            // System.out.println 
            //   ("Corresponding nearest neighbors dissimilarities:");
            // DataAnalysis.printVect(nndiss, 4, 10);

	    } // End of agglomerations loop 
            while (ncl > 1);
	// This completes the loop to carry out the series of 
	// nrow-1 agglomerations.

    }  // End of method reciprocalNearestNeighbors in class HierClass


    double printClustering (int nrow, double[] levels, JTextArea outext2,
    int[][] clusters, int[][] clusters2, double[] card, int[] Acl, int[] Bcl, 
                          String[] rowlabs, int[] order) 
        {

        // Prepare the level indices for output 
	double totint = 0.0;
	for (int i = 0; i < nrow-1; i++) {
	    totint += levels[i];
	}
        // Set up display formating for totint.  
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        nf.setMaximumFractionDigits(4);
        nf.setMinimumFractionDigits(4);
	String num = nf.format(totint);
	// num (above) is formatted; totint (below) is identical but 
	// unformatted.
	// System.out.println("The sum of level indices is " + totint);
	outext2.append
           ("The sum of level indices is " + num + "\n");
        // ("La somme des indices de niveau est " + num + "\n");
	outext2.append
           ("Rates T are scaled by 10^-4 \n");
        // ("Les taux T sont compt\u00e9s en 10^-4 \n");
	for (int i = 0; i < nrow-1; i++) {
	    levels[i] /= totint;
	}
	double[] seqnums = new double[nrow-1];
	for (int i = 0; i < nrow-1; i++) {
	    seqnums[i] = 2.0*(double)nrow - 1.0 - (double)i;
	}
	int collo = 0;
	int colhi = nrow - 2;
	int abswidth = 5;
	String lab = "c   ";
	DataAnalysis.printVect
               (seqnums, collo, colhi, abswidth, 1.0, lab, outext2);
        lab = "car ";
	DataAnalysis.printVect
               (card, collo, colhi, abswidth, 1.0, lab, outext2);
	lab = "T   ";
        DataAnalysis.printVect
               (levels, collo, colhi, abswidth, 10000.0, lab, outext2);
	// The following can be picked up from method 
	// reciprocalNearestNeigbors if of interest.
	// They provide an alt. sequencing of the non-singleton clusters.
        // lab = "A   ";
        // printVect( A, collo, colhi, abswidth, 1, lab, outext2);
        // lab = "B   ";
        // printVect( B, collo, colhi, abswidth, 1, lab, outext2);
	abswidth = 5; 
	lab = "A   ";
	DataAnalysis.printVect( Acl, abswidth, rowlabs, lab, outext2);	
	lab = "B   ";
        DataAnalysis.printVect( Bcl, abswidth, rowlabs, lab, outext2);

	// Prepare the observation labels for outputing
	String myTitle = "";
	String valString; 
        for (int i = 0; i < nrow; i++) {
	    valString = rowlabs[i];
	    abswidth = 4;
	    if (valString.length() < abswidth) 
		valString = DataAnalysis.getSpaces
                    (abswidth - valString.length() - 1) + valString;
	    myTitle = myTitle + " " + valString;
	}
	// Transpose the cluster labels array
	// First "clusters" array
        int[][] tclusters = new int[nrow][nrow]; 
        for (int i1 = 0; i1 < nrow; i1++) {
            for (int i2 = 0; i2 < nrow; i2++) 
                      tclusters[i2][i1] = clusters[i1][i2];
        }
	// We now have: 
	// Row 1: just one cluster
        // Row 2: two clusters
        // Row 3: three clusters
        // ...
        // Row nrow: nrow clusters

	// We consider three representations of the hierarchy, based on 
	// observation sequence numbers and cluster labels.
	// We retain output representation 3. 
	// Representations 1 and 2 are useful in determining rep. 3.

	// Output representation 1: cluster sequence number is the 
	// smallest consistuent observation sequence number. 
	// So cluster with all observations is: 1.
	// System.out.println();
	// System.out.println
        // ("Output representation 1 
	//         (cluster seqno. is smallest obs. seqno.):");
	// outext2.append (
        // "Output representation 1 
	//         (cluster seqno. is smallest obs. seqno.): \n");
	// DataAnalysis.printMatrix(nrow, nrow, tclusters, 4, 4, outext2);
	// System.out.println(myTitle);
	// outext2.append(myTitle + "\n");

	// Next seqnos. of clusters up to 2*nrow-1 using "clusters2" array
        for (int i1 = 0; i1 < nrow; i1++) {
            for (int i2 = 0; i2 < nrow; i2++)
                tclusters[i2][i1] = clusters2[i1][i2];
        }

	// Ouput representation 2: cluster labels use sequence numbers 
	// from n+1 through 2*n-1.  Order of singletons is the given 
	// input order. 
	// System.out.println();
	// System.out.println
        //   ("Output representation 2 (cluster seqnos. up to 2*n-1):");
	// outext2.append
        //   ("Output representation 2 (cluster seqnos. up to 2*n-1): \n");
	// DataAnalysis.printMatrix(nrow, nrow, tclusters, 4, 4, outext2);
	// System.out.println(myTitle);
	// outext2.append(myTitle + "\n");

	// Get dendrogram-consistent order of individuals
	// int[] order = new int[nrow];   // Passed as parameter
	int elements = 0;             // Number of elements so far in order
	for (int i1 = 0; i1 < nrow; i1++) order[i1] = 0;
	// Proceed downwards back through sequence of agglomerations,
	// expanding cluster sequence labels at each level.
	order[0] = Acl[0];  // First we use level 0, top non-trivial level
	order[1] = Bcl[0]; 
	elements = 2;       // Level 0 gives us two elements 
	for (int i1 = 1; i1 < nrow-1; i1++) {        // Now, from level 1
	    for (int i2 = 0; i2 < elements; i2++) {  // Find level in order
		if (order[i2] == (2*nrow-1-i1)) {    // Expand for 2 new elts
		    for (int i3 = elements; i3 > i2; i3--) {  // Make room
			order[i3] = order[i3-1];
		    }                                // End of making room
		    order[i2] = Acl[i1];             // Now insert 2 new elts
		    order[i2+1] = Bcl[i1];
		    elements += 1;                   // One extra locn. 
		}                                    // Continue expanding
	    }                                        // Keep searching order
	}                                     // Move on to next level down
	// Reorder assignment array
	for (int i1 = 0; i1 < nrow; i1++) {
	    for (int i2 = 0; i2 < nrow; i2++) {
		clusters[i1][i2] = tclusters[i1][order[i2]-1];
	    }
	}
	// Reorder observation labels 
	myTitle = "";
	//	String valString; 
        for (int i = 0; i < nrow; i++) {
	    valString = rowlabs[order[i]-1];
	    abswidth = 4;
	    if (valString.length() < abswidth) 
		valString = DataAnalysis.getSpaces
                    (abswidth - valString.length() - 1) + valString;
	    myTitle = myTitle + " " + valString;
	}

	// Output representation 3: clusters are labeled n+1 to 2*n-1,
	// and singletons are ordered such that there are no cross-overs
	// in the dendrogram.
	// System.out.println
	// ("Output representation 3 (dendrogram order of observations):");
	// outext2.append
        // ("Output representation 3 (dendrogram order of observations): \n");
	outext2.append 
           ("\nRepresentation of the hierarchy. \n");
        // ("\nRepr\u00e9sentation de la hi\u00e9rarchie. \n");
        outext2.append
           ("Classes are labeled n+1 to 2*n-1.\n");
        // ("Les classes sont libell\u00e9es n+1 jusqu'\u00e0 2*n-1.\n");
	// 4, 4 correspond to number of digits, and space allowed per value.
	// Note: labels assumed to be 4-character.  
	// Problems will arise if 2*n - 1 becomes 4-valued (i.e. 1000).  
	// I.e. n up to and including 500 is okay.
	DataAnalysis.printMatrix(nrow, nrow, clusters, 4, 4, outext2);
	outext2.append(myTitle + "\n");

	// Method returns total inertia: 
	return totint; 

    }  // End of method printClustering in class HierClass


    //--------------------------------------------------------------------
    /** 
     * Method Dissim, calculates dissimilarity n x n array
     * @param nrow integer row dimension
     * @param ncol integer column dimension
     * @param mass double row masses
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
	        Adiss[i1][i2] += (mass[i1]*mass[i2])/(mass[i1]+mass[i2])*
	                       Math.pow(A[i1][j] - A[i2][j], 2.0); }
	   // Adiss[i1][i2] += 0.5*Math.pow(A[i1][j] - A[i2][j], 2.0); }
	   // Fill the other half of the array 
           Adiss[i2][i1] = Adiss[i1][i2];
         }
       }
       return Adiss;
    } // Dissim 


    //-----------------------------------------------------------------

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



    //------------------------------------------------------------------

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
           int[][] clusters2, int clust1, int clust2, int ncl, 
           int[] Acl, int[] Bcl)
    {
        // If either clust* is not initialized, then we must init. clusters
	if ((clust1 == 0) || (clust2 == 0)) {
           for (int j = 0; j < nrow; j++) {
               for (int i = 0; i < nrow; i++) {
                   clusters[i][j] = 0;
		   clusters2[i][j] = 0;
	       }
	   }
	   // Adjust for 0-sequencing in extreme right-hand col values.
           for (int i = 0; i < nrow; i++) {
               clusters[i][ncl-1] = i + 1;
               clusters2[i][ncl-1] = i + 1;
	   }
           return;
	}

        // For some agglomeration, we are told that label clust1 and 
        // label clust2, among all labels in col. ncl-1 (ncl ranges over
        // 0 thru nrow-1) are to be agglomerated
	// We have arranged that order is such that clust1 < clust2
	// Note: clust1 and clust2 are 1-sequenced and not 0-sequenced 

        int ncl1; 
        ncl1 = ncl - 1; 
	int newcl = 2*nrow -1 -ncl1;
        for (int i = 0; i < nrow; i++) {
	    // Start by replicating previous level:
            clusters[i][ncl1] = clusters[i][ncl];
	    clusters2[i][ncl1] = clusters2[i][ncl];
	    // In clusters2, relabel both new clusters, and control 
	    // this using the clusters seq. nos.  BEFORE updates to clusters.
            if (clusters[i][ncl1] == clust2) {
		// Get old label from clusters2 into Bcl before update.
		  Bcl[ncl1] = clusters2[i][ncl1];
                  clusters2[i][ncl1] = newcl;
	    }
            if (clusters[i][ncl1] == clust1) { 
		// Get old label from clusters2 into Acl before update.
		Acl[ncl1] = clusters2[i][ncl1];
                clusters2[i][ncl1] = newcl;
	    }

	    // In clusters, new cluster clust1 is fine; but relabel clust2:
            if (clusters[i][ncl1] == clust2) clusters[i][ncl1] = clust1;
        }

	// By convention we will require Acl < Bcl 
	int temp = 0;
        for (int i = 0; i < nrow-1; i++) {
	    if (Acl[i] > Bcl[i]) {
		temp = Acl[i];
		Acl[i] = Bcl[i];
		Bcl[i] = temp;
	    }
	}

        // Return type void => no return stmt
    } // ClustMat
    //------------------------------------------------------------------       


}  // End of class HierClass

