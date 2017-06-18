import java.io.*;
import java.util.*;
import java.text.*;

/** 
 * resultSVD - result of eigen-reduction
 */

public class ResultSVD {

    double[][] evectors;
    double[] evalues;
    double[][] CP;
    double[] rates; 
    double trace;

    public ResultSVD ( double[][] evectors, double[] evalues,
		       double[][] CP, double[] rates, double trace)
    {

	this.evectors = evectors;
	this.evalues = evalues;
	this.CP = CP;
	this.rates = rates; 
	this.trace = trace; 

    }

    /** 
     * Public methods
     */

	public double[][] getEvectors() {
	    return evectors; 
	}

	public double[] getEvalues() {
	    return evalues; 
	}

	public double[][] getCP() {
	    return CP; 
	}

	public double[] getRates() {
	    return rates; 
	}

	public double getTrace() {
	    return trace; 
	}

    }
