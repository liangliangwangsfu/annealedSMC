package pty.smc.models;

import pty.ObservationDimensions;

public class BrownianModel implements ObservationDimensions {
	final public int nsites;
	final public double varianceScale;
	public BrownianModel (int sites, double varianceScale) {
		this.nsites = sites;
		this.varianceScale = varianceScale;
	}
//	 public BrownianModel (int sites, double variance) {
//	    nsites = sites;
//	    variances = new double[nsites];
//	    for (int i =0; i < nsites;i++)
//	      variances[i] = variance;
//	  }
	
	public int nSites (){
		return nsites;
	}
	
	// Not meaningful for Brownian
	public int nCharacter (int site) {
		return 2;
	}
}
