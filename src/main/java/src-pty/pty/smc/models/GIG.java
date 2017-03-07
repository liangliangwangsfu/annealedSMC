package pty.smc.models;

public class GIG {
	// GIG (a,b,p)
	
	public static double GIGLognorm (double a, double b, double p) {
		double lognorm = Math.log(2) + VarianceMarginalUtils.logModifiedBesselSecondKind(p, Math.sqrt(a*b)) 
		- 0.5 * p * (Math.log(a) - Math.log(b));
		return lognorm;
	}
	
	public static double GIGapproxLognorm (double a, double b, double p) {
		double lognorm = 0.5 * ( Math.log (4*Math.PI) - Math.log(2*p+1)) 
							+ p * (Math.log(2) + Math.log(p) - Math.log(a) - 1);
		return lognorm;
		
	}

}