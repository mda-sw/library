import java.io.*;
import java.util.*;
import java.text.*;

/**
 * Header
 <P>
 DataAnalysisJava is a first version of a data analysis package for 
 Java.  The Header class handles data and operations related to the
 header information in the input data file. 
 @author Fionn Murtagh, fmurtagh@acm.org
 @version 0.0, 2003-March-15
 */


public class Header {

    String title;                   // Title of analysis
    int nrows;                      // Number of rows
    int ncols;                      // Number of columns
    String suppIndicator;           // 'SUPC' or 'SUPR' or (default) 'SUPNO' 
    int nfactors;                   // Number of factors
    int nclusters;                  // Number of clusters
    String[] columnIdentifiers;     // Column labels

    public Header (int nrows, int ncols, int nfactors, int nclusters,
                    String title, String suppIndicator, 
                    String[] columnIdentifiers)
    {
	this.nrows = nrows; 
	this.ncols = ncols; 
	this.nfactors = nfactors;
	this.nclusters = nclusters;
	this.suppIndicator = suppIndicator; 
	this.title = title;
	this.columnIdentifiers = columnIdentifiers;
    }


    /* -------------------
      Public Methods 
     * ------------------- */

    public String getTitle () {
	 return title;
    }

    public int getNrows () {
	 return nrows;
    }

    public int getNcols () {
	return ncols; 
    }

    public int getNfactors () {
	return nfactors; 
    }

    public int getNclusters () {
	return nclusters; 
    }

    public String getSuppIndicator() {
	return suppIndicator; 
    }

    public String[] getColumnIdentifiers () {
	return columnIdentifiers; 
    }

}


