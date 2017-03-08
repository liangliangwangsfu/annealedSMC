package pty.smc.models;

import nuts.math.TrapezoidLogSpaceIntegrator;
import nuts.util.MathUtils;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.special.Gamma;


public class Tet {
	
	public Tet () {
	}
	
	
	
	public double logGammaDensity (double x, double a, double b ) {
		double y = (a-1)*Math.log(x) - b*x + a * Math.log(b) - Gamma.logGamma(a);
		return y;
	}
	
	public double logNormalDensity (double x, double a, double b) {
		double y = - Math.pow(x-a,2) /(2*b) - 0.5 * Math.log (2*Math.PI*b);
		return y;
	}
	
	public void compareDensityApproximations (double a, double b, double m) {
		double[] x = new double[1000];
		for (int i  = 0; i < x.length; i++)
			x[i] = i*0.05;
		
		for (int i = 0 ; i < x.length; i++) {
			double n = (m-2)/2;
			double ygamma = logGammaDensity ((n + x[i]*Math.sqrt(n))/b, n, 0.5*b) + 0.5 * (Math.log(n) - Math.log(b*b));
			double ynormal = logNormalDensity (x[i], 0, 1) ;
			
			
			System.out.println ("gamma = " + ygamma + "\tnormal = " + ynormal);
		}
		
	}
	

	public double numericalIntegrate (double a, double b, double m, double c) {	
		//LegendreGaussIntegrator lg = new LegendreGaussIntegrator (5, 100);
		TrapezoidLogSpaceIntegrator trap = new TrapezoidLogSpaceIntegrator(
												new GIGLogDensity(b,a,0.5*m-1));
		double i;
		try {
			//i = Math.log(t.integrate(new GIGDensity(b,a,0.5*m-1), 0, 1/c)) ;
			i = trap.integrate (0, 1/c) ;
			m = MathUtils.nint(m);
			System.out.println ("p=" + (0.5*m-1) + "m=" + m);
			double lognorm = Math.log(2) + VarianceMarginalUtils.logModifiedBesselSecondKind(0.5*m-1, Math.sqrt(a*b)) 
							- 0.25 * (m-2) * (Math.log(b) - Math.log(a));
			i = i  - lognorm;
		} catch (Exception e) {
			System.out.println ("error");
			e.printStackTrace ();
			i = -1;
		}
		return i;				
	}

	public double analyticalIntegrateNormal (double a, double b, double m, double c) {
		double result = 0;
		NormalDistributionImpl normal = new NormalDistributionImpl((m-2)/b, Math.sqrt(2*(m-2))/b );
		
		try {
			result += Math.log (normal.cumulativeProbability (0,1/c));		
		} catch (Exception e ) {
			e.printStackTrace ();
		}		
		return result;
	}
	
	
	public double analyticalIntegrateGamma (double a, double b, double m, double c) {
		double result;
		GammaDistributionImpl gamma = new GammaDistributionImpl(0.5*m-1, 1/(0.5*b));
		try {
			result = Math.log (gamma.cumulativeProbability (0,1/c));		
			
		} catch (Exception e ) {
			e.printStackTrace ();
			result = -1;
		}
		
		return result;
	}
	
	
	public void compareIntegralApproximations (double a, double b, double m) {
		double[] c = new double[10];
		for (int i  = 0; i < c.length; i++)
			c[i] = (i+1)*0.05;
		
		for (int i = 0 ; i < c.length; i++) {
			//double numerical = numericalIntegrate ( a, b, m, c[i]);
			double numerical = 0;
			double gamma = analyticalIntegrateGamma (a, b, m, c[i]);
			double normal = analyticalIntegrateNormal (a, b, m, c[i]);
			double normal0 = analyticalIntegrateNormal (a, b, m, 0);
			
			double normal2 = normalcdf2 ((1/c[i] - (m-2)/b)/Math.sqrt(2*(m-2))/b);
			double normal20 = normalcdf2 ((- (m-2)/b)/Math.sqrt(2*(m-2))/b);
			normal2 = Math.log(normal2 - normal20);
			
			double zorder = logGammaDensity(1/c[i], .5*(m-2), 0.5*b) + (1/c[i]);
			
			System.out.println ("c = " + c[i] + 
						"\tNumerical = " + numerical + 
						"\tGamma approximation = " + gamma +
						"\tNormal approximation = " + normal2 +
						"\tzorder = " + zorder);
		}
		
	}
	
	private class GIGDensity implements UnivariateRealFunction {
		double a,b,m;	
		public GIGDensity (double a, double b, double m) {
			this.a = a;
			this.b = b;
			this.m = m;
		}
		
		public double value (double x) {
			return Math.exp (-b/x - a*x) * Math.pow (x,m-1);
			
		}
	}
	
	private class GIGLogDensity implements UnivariateRealFunction {
		double a,b,m;	
		public GIGLogDensity (double a, double b, double m) {
			this.a = a;
			this.b = b;
			this.m = m;
		}
		
		public double value (double x) {
			return -b/x - a*x  + (m-1)*Math.log(x);
			
		}
	}

	public static double normalcdf1 (double x) {
		double p0 = 1.253314137;
		double a = 0.212023887;
		double b = 0.282455120;
		double x2 = x*x;
		
		double y = Math.pow(p0*x, 2)	+ Math.exp(-0.5*x2)*Math.sqrt(1+b*x2)/(1+a*x2);
		y = Math.exp(-0.5*x2) - Math.sqrt(y);
		y = p0 + y/x;
		
		return y;
		
	}
	
	public static  double normalcdf2 (double x) {
		double y = 0.07056*Math.pow(x, 3) +1.5976*x;
		y = 1/(1+Math.exp(-y));		
		return y;

	}

	
	public static  double normalcdf3 (double x) {
		double y = 1.702*x;
		y = 1/(1+Math.exp(y));		
		return y;

	}
	
	
	public static void main (String args[]) {
		System.out.println ("Running Test");
		try {
		
				Tet t = new Tet ();
				t.compareIntegralApproximations(1, 1, 999);
			//	t.compareDensityApproximations(1, 1, 999);

		} catch (Exception e) {
				e.printStackTrace();			
		}
	}
}
