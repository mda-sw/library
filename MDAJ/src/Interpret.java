import java.io.*;
import java.util.*;
import java.text.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
// import com.macfaq.io.*;

/**
 * Interpret
 <P>
 DataAnalysisJava is a first version of a data analysis package for 
 Java.  Class Interpret provides interpretation of factors in terms
 of clusters (FACOR) and variables (VACOR). 
 <P>
 @author Fionn Murtagh, fmurtagh@acm.org
 @version 0.0, 2003-March-15
 */


public class Interpret {


    void FacorVacor (int nrow, int mnew, int nclusters, 
		     int nfactors, double[] rowsums, double[] colsums, 
                     double[][] indatstd, double[] cinertia,
		     int[] order, int[][] clusters, 
		     int[] Acl, int[] Bcl, double[][] Evex, double[] Evals,
                     double totint, String[] collabs, 
		     JTextArea outext3, JTextArea outext4) 
     {


	// --------------------------------------------------------------
	if (nclusters < 2) {
	    outext3.append
		("Unacceptably small number of clusters; setting to 2.");
	    nclusters = 2;
	}
	if (nclusters > nrow-1) {
	    outext3.append
		("Unacceptably large number of clusters; setting to " +
		 (nrow-1) + ".");
	    nclusters = nrow-1;
	}

	// --------------------------------------------------------------
	// nclusters is the selected no. of clusters for analysis
	// Row 0 of table clusters corresponds to 1 cluster
	// Row 1 of table clusters corresponds to 2 clusters
	// Etc.
	// nclusters = 4 (e.g.) implies that we analyze row nclusters-1
	// We will use an array of dimensions 2*nclusters-1 which will
	// related to 
	// firstly - top nclusters-1 clusters/nodes in the hierarchy, and
	// secondly - the nclusters components of the specified partition.

	double[] origin = new double[mnew];
	// Origin is fJ:
	for (int j = 0; j < mnew; j++) origin[j] = colsums[j];


	// --------------------------------------------------------------
	// Analysis of top nodes 
	// --------------------------------------------------------------

        // Analysis for top nodes in hierarchical tree 
        // 2*n-1, 2*n-2, 2*n-3, ... 2*n-(nclusters-1)
        // Number of these top nodes: nclusters-1
        // We don't need to care about the order of the observations
        double[][] topnodes = new double[nclusters-1][mnew];
        double[] topnodeswts = new double[nclusters-1];
        int[] topnodescard = new int[nclusters-1]; 
        for (int k = 0; k < nclusters-1; k++) {
	    topnodeswts[k] = 0.0;
	    topnodescard[k] = 0;
	    for (int j = 0; j < mnew; j++) {
	        topnodes[k][j] = 0.0;
	    }
        }
        int klust;
        for (int k = 0; k < nclusters-1; k++) {
	    klust = 2*nrow - 1 - k;            // node seq. number
	    for (int i = 0; i < nrow; i++) {  
	        if (clusters[k][i] == klust) {
                   for (int j = 0; j < mnew; j++) {
                   topnodes[k][j] += indatstd[order[i]-1][j];
		   }
		topnodeswts[k] += rowsums[order[i]-1];
		topnodescard[k]++;
		}
	    }
	}
        for (int k = 0; k < nclusters-1; k++) {
	    for (int j = 0; j < mnew; j++) {
	        topnodes[k][j] /= topnodeswts[k];
	    }
        }

        double[][] topnodesproj = new double[nclusters-1][mnew];
        double[][] topnodescor = new double[nclusters-1][mnew];
        double[][] topnodesctr = new double[nclusters-1][mnew];
        DataAnalysis.getCoordInf
                (nclusters-1, mnew, topnodesproj, topnodescor, topnodesctr,
		 topnodes, Evex, Evals, origin, topnodeswts, topnodescard);
        // We need to zero topnodescor of the trivial partition:
        for (int j = 0; j < mnew; j++) topnodescor[0][j] = 0.0; 

	// --------------------------------------------------------------
	// Formatting of output for top nodes 
	// --------------------------------------------------------------

        // Collect all class results in array topres - "top of tree results"
        int mbigcl = 3*nfactors + 5;   
        int cl = 0; int clA = 0; int clB = 0;
        // Now build up big array for outputing:
        double[][] topres = new double[nclusters-1][mbigcl];
        int[][] toplabs = new int[nclusters-1][3]; 
        for (int k = 0; k < nclusters-1; k++) {
	    // QLT has still to be completed - return to this later.
            topres[k][0] = 1.0;                         // QLT
	    topres[k][1] = topnodeswts[k];              // PDS or WTS
	    double dst = 0.0;                           // Det. INR 
	    for (int j = 0; j < mnew; j++)
                dst += Math.pow((topnodes[k][j]-origin[j]), 2.0)/origin[j];
            topres[k][2] = topnodeswts[k]* dst / totint;  // INR 
            for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
                topres[k][3*j+2] = topnodesproj[k][j];      // F#1, F#2, ...
                topres[k][3*j+3] = topnodescor[k][j];       // CO2
                topres[k][3*j+4] = topnodesctr[k][j];       // CTR
	    }
	    for (int k2 = 0; k2 < nclusters-1; k2++) {
                toplabs[k2][0] =  2*nrow - 1 - k2;   
                toplabs[k2][1] = Acl[k2];
                toplabs[k2][2] = Bcl[k2];
	    }
        }


	// --------------------------------------------------------------
	// Analysis of clusters of partition
	// --------------------------------------------------------------

        // Analysis for clusters of selected partition

	// Determine the nclusters labels
	int[] temprow = new int[nrow];
	int[] clusnums = new int[nclusters];
	for (int i = 0; i < nrow; i++) temprow[i] = clusters[nclusters-1][i];

        DataAnalysis.inSort(temprow);

	clusnums[0] = temprow[0];   // First cluster number = first entry
	int foundsofar = 1;  
	for (int i = 1; i < nrow; i++) {
	    if (temprow[i] != temprow[i-1]) {
		clusnums[foundsofar] = temprow[i];
		foundsofar += 1;
	    }
	}


	int[] cluscard = new int[nclusters];
	for (int k = 0; k < nclusters; k++) {
	    cluscard[k] = 0;
	    for (int i = 0; i < nrow; i++) {
		if (clusters[nclusters-1][i] == clusnums[k]) {
		    // System.out.println("Obs " + order[i] + 
                    // " is in cluster " + clusnums[k]);
		    cluscard[k]++; 
		}
	    }
	}

        double[][] clusmeans = new double[nclusters][mnew];
	double[] cluspoids = new double[nclusters];
	for (int k = 0; k < nclusters; k++) {
	    cluspoids[k] = 0.0;
	    for (int j = 0; j < mnew; j++) {
                clusmeans[k][j] = 0.0;
	    }
	}
	for (int k = 0; k < nclusters; k++) {
	        for (int i = 0; i < nrow; i++) {
		    if (clusters[nclusters-1][i] == clusnums[k]) {
			// Obs index is given by: order[i]-1
	                for (int j = 0; j < mnew; j++)
                            clusmeans[k][j] += 
                            indatstd[order[i]-1][j];
	            cluspoids[k] += rowsums[order[i]-1];
		    }
		}
	}
	for (int k = 0; k < nclusters; k++) {
	    for (int j = 0; j < mnew; j++) clusmeans[k][j] /= cluspoids[k];
	}

        double[][] clusproj = new double[nclusters][mnew]; 
	double[][] cluscor  = new double[nclusters][mnew];
	double[][] clusctr  = new double[nclusters][mnew];
	DataAnalysis.getCoordInf(nclusters, mnew, clusproj, cluscor, clusctr, 
		    clusmeans, Evex, Evals, origin, cluspoids, cluscard); 

	//for (int i1 = 0; i1 < nclusters; i1++) {
	//    for (int i2 = 0; i2 < 3; i2++) {
	//	System.out.println("i1, i2, clusproj-i1-i2 = " + i1 + 
	//			   " " + i2 + " " + clusproj[i1][i2]); 
	//    }
	//}


	// --------------------------------------------------------------
	// Formatting of output for clusters of partition 
	// --------------------------------------------------------------

        // Collect all class results in array clusres
        // int mbigcl = 3*nfactors + 5;       // No change
        cl = 0; clA = 0; clB = 0;
        // Now build up big array for outputing:
        double[][] clusres = new double[nclusters][mbigcl];
        int[][] cluslabs = new int[nclusters][3]; 
        for (int k = 0; k < nclusters; k++) {
	    // QLT has still to be completed - return to this later.
            clusres[k][0] = 1.0;                         // QLT
	    clusres[k][1] = cluspoids[k];                // PDS or WTS
	    double dst = 0.0;                            // Det. INR 
	    for (int j = 0; j < mnew; j++)
                dst += Math.pow((clusmeans[k][j]-origin[j]), 2.0)/origin[j];
            clusres[k][2] = cluspoids[k]* dst / totint;  // INR 
            for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
                clusres[k][3*j+2] = clusproj[k][j];      // F#1, F#2, ...
                clusres[k][3*j+3] = cluscor[k][j];       // CO2
                clusres[k][3*j+4] = clusctr[k][j];       // CTR
	    }
	    if (clusnums[k] > nrow) {  // Case of cl != singleton.
	        cluslabs[k][0] = clusnums[k];
	        cluslabs[k][1] = Acl[2*nrow-clusnums[k]-1];
	        cluslabs[k][2] = Bcl[2*nrow-clusnums[k]-1];
	    }
	    else {   // Case of cl = singleton, so no elder and younger.
	        cluslabs[k][0] = clusnums[k];
	        cluslabs[k][1] = 0;
	        cluslabs[k][2] = 0;
	    }
        }
	// For number formatting in labels 
        NumberFormat nf = NumberFormat.getInstance();
        nf.setGroupingUsed(false);
        // String myTitle = "|CLAS AINE BNJM| QLT PDS INR|"; 
        String myTitle = "|CLAS ELDR YNGR| QLT WTS INR|"; 
        for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
            if (nfactors < 10) {
               nf.setMaximumIntegerDigits(1);
               myTitle = myTitle + "  F#" + nf.format(j) +" CO2 CTR|"; 
	    }
            if (nfactors > 9) {
               nf.setMaximumIntegerDigits(2);
               myTitle = myTitle + "  F" + nf.format(j) +" CO2 CTR|"; 
	    }
        } 
        for (int j = 0; j < (5*mbigcl-3); j++) outext3.append("_");
        outext3.append("\n");
        outext3.append(myTitle + "\n");
        for (int j = 0; j < (5*mbigcl-3); j++) outext3.append("_");
        outext3.append("\n");
        outext3.append
             ("Representation on the factorial axes of the " + 
	//   ("Repr\u00e9sentation sur les axes factoriels des " + 
             (nclusters-1) + 
             " selected nodes \n");
        //   " noeuds choisis \n");
        DataAnalysis.printAll
             (nclusters-1, mbigcl, topres, toplabs, 1000.0, outext3); 
        outext3.append
             ("Representation on the factorial axes of the " + 
	      // ("Repr\u00e9sentation sur les axes fact des " + 
             nclusters + 
             " classes of the selected partition \n");
	     // " classes de la partition choisie \n");
        DataAnalysis.printAll
             (nclusters, mbigcl, clusres, cluslabs, 1000.0, outext3); 
        for (int j = 0; j < (5*mbigcl-3); j++) outext3.append("_");
        outext3.append("\n");

	// --------------------------------------------------------------
	// Analysis of dipoles
	// --------------------------------------------------------------

        // Dipole analysis

        double[][] selnodes = new double[3*nclusters-3][mnew];
        double[] selnodeswts = new double[3*nclusters-3];
        int[] selnodescard = new int[3*nclusters-3]; 
        for (int k = 0; k < 3*nclusters-3; k++) {
	    selnodeswts[k] = 0.0;
	    selnodescard[k] = 0;
	    for (int j = 0; j < mnew; j++) {
	        selnodes[k][j] = 0.0;
	    }
        }
        int k3 = 0;
        for (int k = 0; k < nclusters-1; k++) {
	    for (int k2 = 0; k2 < 3; k2++) {
	        k3 = 3*k + k2;
	        klust = toplabs[k][k2];            // node seq. number

		// klust ranges over 
		// cluster no., aine no., benjamin no. 
		// for top clusters, listed in that order.  

		if (klust <= nrow) {
		    for (int j = 0; j < mnew; j++) {
                        selnodes[k3][j] = indatstd[klust-1][j]; 
		    }
                    selnodeswts[k3] = rowsums[klust-1];
                    selnodescard[k3]++;
		}  // End of: klust <= nrow (klust singleton).

		if (klust > nrow) { 
	            for (int i = 0; i < nrow; i++) {  
	                if (clusters[2*nrow-klust-1][i] == klust) {
                           for (int j = 0; j < mnew; j++) {
                               selnodes[k3][j] += indatstd[order[i]-1][j];
		           }
		           selnodeswts[k3] += rowsums[order[i]-1];
		           selnodescard[k3]++;
	                 }
                    }
		}  // End of: klust > nrow (klust non-singleton).

	    }
        }
        for (int k = 0; k < 3*nclusters-3; k++) {
	    for (int j = 0; j < mnew; j++) {
	        selnodes[k][j] /= selnodeswts[k];
	    }
        }

        // printMatrix(3*nclusters-3, mnew, selnodes, 4, 4);
        double[][] selproj = new double[3*nclusters-3][mnew];

        for (int k = 0; k < 3*nclusters-3; k++) {
            for (int j = 0; j < mnew; j++) {
                selproj[k][j] = 0.0; 
                for (int j2 = 0; j2 < mnew; j2++) 
                    selproj[k][j] += selnodes[k][j2]*Evex[j2][j];
	     }
        }

        // Our storage scheme: nodes n are 0, 3, 6, ... 
        // i.e. 3*k for k = 0, 1, 2, ... and 
        // corresponding A[n] and B[n] are resp. 3*k+1 and 3*k+2 : 
        // 0: 1 - 2
        // 3: 4 - 5
        // 6: 7 - 8

        double[][] seldiff = new double[nclusters-1][mnew];
        double[][] seldprj = new double[nclusters-1][mnew];
        double[][] seldctd = new double[nclusters-1][mnew];
        double[] selQLD = new double[nclusters-1]; 
        double[] selPDS = new double[nclusters-1];
        double[] selIND = new double[nclusters-1];
        int[][] selcl = new int[nclusters-1][3];
        double dst; 

        // Pick out node seq. no. together with the node's A[n] and B[n]:
        for (int k = 0; k < nclusters-1; k++) {
	    for (int k2 = 0; k2 < 3; k2++) {
	        k3 = 3*k + k2;
	        selcl[k][k2] = toplabs[k][k2];            // node seq. number
	    }
        }

        for (int k = 0; k < nclusters-1; k++) {
	    // Aagh! Weeks of work because following init was before k-loop!!
	    // So dst was being accumulated for each node.   Aaaagh!!!
            dst = 0.0; 
	    selQLD[k] = 0.0;
	    selIND[k] = 0.0;

	    for (int j = 0; j < mnew; j++) {
	        // Here we are determining distance in the space of variables:
	        dst += Math.pow((selnodes[3*k+1][j] - selnodes[3*k+2][j]),2.0)/
                     colsums[j];
	        //  Equally acceptable alternative is to use all factors:
	        //  dst += Math.pow((selproj[3*k+1][j]-selproj[3*k+2][j]),2.0);
	    }

	    for (int j = 0; j < mnew; j++) {

	        // seldiff is: Dalpha(n) = Falpha(A[n]) - Falpha(B[n])
	        // i.e. difference between projections of A and B on factors:
	        seldiff[k][j] = selproj[3*k+1][j] - selproj[3*k+2][j];

	        // seldprj is CODalpha(n) = 
	        // (Falpha(A[n]) - Falpha(B[n]))^2 / || A[n] - B[n] ||^2.  
	        // Numerator is squared projection difference term.
	        // Denominator has been calculated above before this loop. 
	        seldprj[k][j] = 
		    Math.pow((selproj[3*k+1][j] - selproj[3*k+2][j]), 2.0)/dst;

	        // selproj[3*k][j] is proj of n
	        // selproj[3*k+1][j] is proj of A[n],mass is selnodeswts[3*k+1]
	        // selproj[3*k+2][j] is proj of B[n],mass is selnodeswts[3*k+2]
	        // In following, determine internal inertia of dipole, i.e.
	        // inertia of the system of its two extremities, A[n] and B[n],
	        // projected onto factor alpha, with unchanged masses,
	        // where inertia is w.r.t. centre of gravity n projected also.
	        // Normalize by total inertia of cloud projected onto alpha, 
	        // i.e. lambda_alpha.
	        seldctd[k][j] = 
                    (selnodeswts[3*k+1]*
                     Math.pow(selproj[3*k+1][j]-selproj[3*k][j],2.0) +
                     selnodeswts[3*k+2]*
                     Math.pow(selproj[3*k+2][j]-selproj[3*k][j],2.0))/
  	             Evals[j];

	        selQLD[k] += seldprj[k][j]; 
	        selPDS[k] = selnodeswts[3*k];   // Weight of class n.
	        // selIND is internal inertia of dipole, divided by 
		// total inertia:
	        selIND[k] += seldctd[k][j]*Evals[j]/totint;
	    } 
	}


	// --------------------------------------------------------------
	// Formatting of output of analysis of dipoles
	// --------------------------------------------------------------

        // Collect all class results in array clusres
        // int mbigcl = 3*nfactors + 5;       // No change
        // Now build up big array for outputing:
        // double[][] clusres = new double[nclusters][mbigcl];  // No change
        // int[][] cluslabs = new int[nclusters][3]; // No change
        for (int k = 0; k < nclusters-1; k++) {
            clusres[k][0] = selQLD[k];                   // QLT  
	    clusres[k][1] = selPDS[k];                   // PDS or WTS
	    clusres[k][2] = selIND[k];                   // INR 
            for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
                clusres[k][3*j+2] = seldiff[k][j];       // D#1, D#2, ...
                clusres[k][3*j+3] = seldprj[k][j];       // COD
                clusres[k][3*j+4] = seldctd[k][j];       // CTD
	    }
            cluslabs[k][0] = selcl[k][0];
	    cluslabs[k][1] = selcl[k][1];
	    cluslabs[k][2] = selcl[k][2];
        }
        // myTitle = "|CDIP AINE BNJM| QLT PDS IND|"; 
        myTitle = "|CDIP ELDR YNGR| QLT WTS IND|"; 
        for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
            if (nfactors < 10) {
               nf.setMaximumIntegerDigits(1);
               myTitle = myTitle + "  D#" + nf.format(j) +" COD CTD|"; 
	    }
            if (nfactors > 9) {
               nf.setMaximumIntegerDigits(2);
               myTitle = myTitle + "  D" + nf.format(j) +" COD CTD|"; 
	    }
        }
        for (int j = 0; j < (5*mbigcl-3); j++) outext3.append("_");
        outext3.append("\n");
        outext3.append(myTitle + "\n");
        for (int j = 0; j < (5*mbigcl-3); j++) outext3.append("_");

        outext3.append("\n");
        outext3.append
            ("Representation on the factorial axes of the " + 
            (nclusters-1) + " selected dipoles \n");
        //  ("Repr\u00e9sentation sur les axes factoriels des " + 
        //  (nclusters-1) + " dip\u00f4les choisis \n");
        DataAnalysis.printAll
            (nclusters-1, mbigcl, clusres, toplabs, 1000.0, outext3); 
        for (int j = 0; j < (5*mbigcl-3); j++) outext3.append("_");
        outext3.append("\n");


     //-------------------------------------------------------------------
     // NEXT... VACOR
     //-------------------------------------------------------------------


	// --------------------------------------------------------------
	// Analysis of top nodes 
	// --------------------------------------------------------------

     // Correlations for top nodes
     for (int k = 0; k < nclusters-1; k++) {  // k ranges over nodes
	 dst = 0.0;
         for (int j = 0; j < mnew; j++)       // j ranges over variables
	     dst += Math.pow(topnodes[k][j]-colsums[j], 2.0)/colsums[j]; 
	 for (int j = 0; j < mnew; j++)   {  
	     topnodescor[k][j] = Math.pow(topnodes[k][j]-colsums[j], 2.0)/
                                                       (colsums[j]*dst);
	    }
     }
     // We need to zero topnodescor of the trivial partition:
     for (int j = 0; j < mnew; j++) topnodescor[0][j] = 0.0; 

     // Correlations for clusters of partition
     for (int k = 0; k < nclusters; k++) {  // k ranges over nodes
	 dst = 0.0;
         for (int j = 0; j < mnew; j++)     // j ranges over variables
	     dst += Math.pow(clusmeans[k][j]-colsums[j], 2.0)/colsums[j]; 
	 for (int j = 0; j < mnew; j++)   {  
	     cluscor[k][j] = Math.pow(clusmeans[k][j]-colsums[j], 2.0)/
                                                       (colsums[j]*dst);
	    }
     }

     // Contributions for top nodes
     for (int k = 0; k < nclusters-1; k++) {  // k ranges over nodes
	 for (int j = 0; j < mnew; j++) 
	     topnodesctr[k][j] = topnodeswts[k]*
                             Math.pow(topnodes[k][j]-colsums[j],2.0)/
                             (colsums[j]*cinertia[j]*totint);
     }

     // Contributions for clusters of partition
     for (int k = 0; k < nclusters; k++) {  // k ranges over clusters
	 for (int j = 0; j < mnew; j++) 
	     clusctr[k][j] = cluspoids[k]*
                             Math.pow(clusmeans[k][j]-colsums[j],2.0)/
                             (colsums[j]*cinertia[j]*totint);
     }


	// --------------------------------------------------------------
	// Formatting of output of analysis of top nodes 
	// --------------------------------------------------------------

     // Format for output printing - first results for top nodes 
     // Collect all class results in array topres - "top of tree results"
     // int mbigcl = 3*nfactors + 5;
     // Here we must have mbigcl = 3*mnew + 5;    
     // int cl = 0; int clA = 0; int clB = 0;
     // Now build up big array for outputing:
     // double[][] topres = new double[nclusters-1][mbigcl];
     // int[][] toplabs = new int[nclusters-1][3]; 

     mbigcl = 3*mnew + 5; 
     double[][] topres2 = new double[nclusters-1][mbigcl]; 

     for (int k = 0; k < nclusters-1; k++) {
	 // QLT has still to be completed - return to this later.
         topres2[k][0] = 1.0;                         // QLT
	 topres2[k][1] = topnodeswts[k];              // PDS or WTS
	 dst = 0.0;                            // Det. INR 
	 for (int j = 0; j < mnew; j++)
             dst += Math.pow((topnodes[k][j]-origin[j]), 2.0)/origin[j];
         topres2[k][2] = topnodeswts[k]* dst / totint;  // INR 
         for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
             topres2[k][3*j+2] = topnodes[k][j-1];      // V#1, V#2, ...
             topres2[k][3*j+3] = topnodescor[k][j-1];       // CO2
             topres2[k][3*j+4] = topnodesctr[k][j-1];       // CTR
	 }
	 for (int k2 = 0; k2 < nclusters-1; k2++) {
             toplabs[k2][0] =  2*nrow - 1 - k2;   
             toplabs[k2][1] = Acl[k2];
             toplabs[k2][2] = Bcl[k2];
	 }
     }

	// --------------------------------------------------------------
	// Analysis of clusters of partition
	// --------------------------------------------------------------

     // Format for output printing - second, clusters of partition 
     // Collect all class results in array clusres
     // int mbigcl = 3*nfactors + 5;       // To be changed to 3*mnew 
     // cl = 0; clA = 0; clB = 0;
     // Now build up big array for outputing:
     // double[][] clusres = new double[nclusters][mbigcl];
     // int[][] cluslabs = new int[nclusters][3]; 

     double[][] clusres2 = new double[nclusters][mbigcl]; 

     for (int k = 0; k < nclusters; k++) {
	 // QLT has still to be completed - return to this later.
         clusres2[k][0] = 1.0;                         // QLT  NO QLT
	 // WE WILL WANT TO SKIP THE ABOVE VALUE.
	 clusres2[k][1] = cluspoids[k];                // PDS or WTS
	 dst = 0.0;                            // Det. INR 
	 for (int j = 0; j < mnew; j++)
             dst += Math.pow((clusmeans[k][j]-origin[j]), 2.0)/origin[j];
         clusres2[k][2] = cluspoids[k]* dst / totint;  // INR 
         for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
             clusres2[k][3*j+2] = clusmeans[k][j-1];      // Vbe#1, Vbe#2, ...
             clusres2[k][3*j+3] = cluscor[k][j-1];       // CO2
             clusres2[k][3*j+4] = clusctr[k][j-1];       // CTR
	 }
	 if (clusnums[k] > nrow) {  // Case of cl != singleton.
	        cluslabs[k][0] = clusnums[k];
	        cluslabs[k][1] = Acl[2*nrow-clusnums[k]-1];
	        cluslabs[k][2] = Bcl[2*nrow-clusnums[k]-1];
	 }
	 else {   // Case of cl = singleton, so no elder and younger.
	        cluslabs[k][0] = clusnums[k];
	        cluslabs[k][1] = 0;
	        cluslabs[k][2] = 0;
	 }
     }


	// --------------------------------------------------------------
	// Formatting of output of clusters of partition
	// --------------------------------------------------------------

     myTitle = "|CLAS ELDR YNGR| QLT WTS INR|"; 
     // myTitle = "|CLAS AINE BNJM| QLT PDS INR|"; 
     String myString;
     for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
         myString = collabs[j-1];
         myString = DataAnalysis.getSpaces(4 - myString.length()) + myString;
         myTitle = myTitle + " " + myString + " CO2 CTR|"; 
     }  
     for (int j = 0; j < (5*mbigcl-3); j++) outext4.append("_");
     outext4.append("\n");
     outext4.append(myTitle + "\n");
     for (int j = 0; j < (5*mbigcl-3); j++) outext4.append("_");
     outext4.append("\n");
     outext4.append
            ("Representation of the " + 
            (nclusters-1) + " nodes on the set J \n");
     //     ("Repr\u00e9sentation des " + 
     //     (nclusters-1) + " noeuds sur l'ensemble J \n");
     DataAnalysis.printAll
            (nclusters-1, mbigcl, topres2, toplabs, 1000.0, outext4); 
     outext4.append
             ("Representation of the " + 
     //      ("Repr\u00e9sentation des " + 
             nclusters + " classes on the set J \n");
     //      nclusters + " classes sur l'ensemble J \n");
     DataAnalysis.printAll
            (nclusters, mbigcl, clusres2, cluslabs, 1000.0, outext4); 
     for (int j = 0; j < (5*mbigcl-3); j++) outext4.append("_");
     outext4.append("\n");

	// --------------------------------------------------------------
	// Analysis of dipoles 
	// --------------------------------------------------------------

     // Next analysis of dipoles

      for (int k = 0; k < nclusters-1; k++) {
	  dst = 0.0;
	  for (int j = 0; j < mnew; j++) {
	      // Here we are determining distance in the space of variables:
	      dst += Math.pow((selnodes[3*k+1][j] - selnodes[3*k+2][j]),2.0)/
                  colsums[j];
	  }
	  for (int j = 0; j < mnew; j++) {
	      // seldiff is: D(n) = pr_vbe1(A[n]) - pr_vbe1(B[n])
	      // i.e. difference between projections of A and B on variables:
	      seldiff[k][j] = selnodes[3*k+1][j] - selnodes[3*k+2][j];
	      seldprj[k][j] = 
		 Math.pow((selnodes[3*k+1][j] - selnodes[3*k+2][j]), 2.0)/
                             (colsums[j]*dst);
	      seldctd[k][j] = 
                 (selnodeswts[3*k+1]*
                 Math.pow(selnodes[3*k+1][j]-selnodes[3*k][j],2.0)/
                        (colsums[j]*cinertia[j]*totint))+
                 (selnodeswts[3*k+2]*
                 Math.pow(selnodes[3*k+2][j]-selnodes[3*k][j],2.0)/
                        (colsums[j]*cinertia[j]*totint));
	  }
      }

	// --------------------------------------------------------------
	// Formatting of output of analysis of dipoles
	// --------------------------------------------------------------

      // Format for output printing
      // Collect all class results in array clusres2
      // int mbigcl = 3*nfactors + 5;       // No change
      // Now build up big array for outputing:
      // double[][] clusres = new double[nclusters][mbigcl];  // No change
      // int[][] cluslabs = new int[nclusters][3]; // No change
      for (int k = 0; k < nclusters-1; k++) {
          clusres2[k][0] = selQLD[k];                   // QLT  
	  clusres2[k][1] = selPDS[k];                   // PDS  or WTS
	  clusres2[k][2] = selIND[k];                   // INR 
          for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
              clusres2[k][3*j+2] = seldiff[k][j-1];       // V#1, V#2, ...
              clusres2[k][3*j+3] = seldprj[k][j-1];       // COD
              clusres2[k][3*j+4] = seldctd[k][j-1];       // CTD
	  }
          cluslabs[k][0] = selcl[k][0];
	  cluslabs[k][1] = selcl[k][1];
	  cluslabs[k][2] = selcl[k][2];
      }
      myTitle = "|CDIP ELDR YNGR| QLT WTS IND|"; 
      // myTitle = "|CDIP AINE BNJM| QLT PDS IND|"; 
      for (int j = 1; j < ((mbigcl-3)/3)+1; j++) {
         myString = collabs[j-1];
         myString = DataAnalysis.getSpaces(4 - myString.length()) + myString;
         myTitle = myTitle + " " + myString + " COD CTD|"; 
      }
      for (int j = 0; j < (5*mbigcl-3); j++) outext4.append("_");
      outext4.append("\n");
      outext4.append(myTitle + "\n");
      for (int j = 0; j < (5*mbigcl-3); j++) outext4.append("_");
      outext4.append("\n");
      outext4.append
         ("Representation of the " + 
         (nclusters-1) + " dipoles on the set J \n");
      //  ("Repr\u00e9sentation des " + 
      //  (nclusters-1) + " dip\u00f4les sur l'ensemble J \n");
      DataAnalysis.printAll
         (nclusters-1, mbigcl, clusres2, toplabs, 1000.0, outext4);
      for (int j = 0; j < (5*mbigcl-3); j++) outext4.append("_");
      outext4.append("\n");
	
     } // End of method FacorVacor in class Interpret

}  // End of class Interpret
