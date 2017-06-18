import java.io.*;
import java.util.*;
import java.text.*;

/** 
 * resultCA - result of Correspondence Analysis
 */ 

public class ResultCA {

    Imatrix projectionsRows;
    Imatrix projectionsColumns;
    Imatrix correlationsRows;
    Imatrix correlationsColumns;
    Imatrix contributionsRows;
    Imatrix contributionsColumns;
    Imatrix projectionsRowsSupp;
    Imatrix projectionsColumnsSupp;
    Imatrix correlationsRowsSupp;
    Imatrix correlationsColumnsSupp;
    Imatrix contributionsRowsSupp;
    Imatrix contributionsColumnsSupp;

    public ResultCA ( Imatrix projectionsRows,
                      Imatrix projectionsColumns,
		      Imatrix correlationsRows,
		      Imatrix correlationsColumns,
		      Imatrix contributionsRows,
		      Imatrix contributionsColumns,
		      Imatrix projectionsRowsSupp,
		      Imatrix projectionsColumnsSupp,
		      Imatrix correlationsRowsSupp,
		      Imatrix correlationsColumnsSupp,
		      Imatrix contributionsRowsSupp,
		      Imatrix contributionsColumnsSupp )
    {

	this.projectionsRows = projectionsRows;
	this.projectionsColumns = projectionsColumns;
	this.correlationsRows = correlationsRows;
	this.correlationsColumns = correlationsColumns;
	this.contributionsRows = contributionsRows;
	this.contributionsColumns = contributionsColumns;
	this.projectionsRowsSupp = projectionsRowsSupp;
	this.projectionsColumnsSupp = projectionsColumnsSupp;
	this.correlationsRowsSupp = correlationsRowsSupp;
	this.correlationsColumnsSupp = correlationsColumnsSupp;
	this.contributionsRowsSupp = contributionsRowsSupp;
	this.contributionsColumnsSupp = contributionsColumnsSupp;

    }

    /** 
     * Public methods
     */

    public Imatrix getProjectionsRows() {
	return projectionsRows; 
    }

    public Imatrix getProjectionsColumns() {
	return projectionsColumns; 
    }

    public Imatrix getCorrelationsRows() {
	return correlationsRows; 
    }

    public Imatrix getCorrelationsColumns() {
	return correlationsColumns; 
    }

    public Imatrix getContributionsRows() {
	return contributionsRows; 
    }

    public Imatrix getContributionsColumns() {
	return contributionsColumns; 
    }

    public Imatrix getProjectionsRowsSupp() {
	return projectionsRowsSupp; 
    }

    public Imatrix getProjectionsColumnsSupp() {
	return projectionsColumnsSupp; 
    }

    public Imatrix getCorrelationsRowsSupp() {
	return correlationsRowsSupp; 
    }

    public Imatrix getCorrelationsColumnsSupp() {
	return correlationsColumnsSupp; 
    }

    public Imatrix getContributionsRowsSupp() {
	return contributionsRowsSupp; 
    }

    public Imatrix getContributionsColumnsSupp() {
	return contributionsColumnsSupp; 
    }

} // End of class: ResultCA


