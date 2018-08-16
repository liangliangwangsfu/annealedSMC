package evolmodel;

import java.util.ArrayList;


public interface EvolutionParameters {
	public double[] getParameters();
	public void setParameters(double[] para);	 


	public static class GTR implements EvolutionParameters
	{		  
//		private double[] subsRates;         // six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
//		private double[] statFreqs;         // stationary state frequencies. pi_A, pi_C, pi_G, pi_T
//		private double alpha=0.5;           // shape parameter in the Gamma distribution
//		private double pInv=0.2;            // the proportion of invariant sites 
        private double[] para = new double[12];
        
        public GTR(double[] para)
        {
        	 this.para = para; 
        }
        
		@Override
		public double[] getParameters() {
			return para;
		}

		@Override
		public void setParameters(double[] para) {
			// para include: six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT, stationary state frequencies. pi_A, pi_C, pi_G, pi_T, shape parameter in the Gamma distribution,  the proportion of invariant sites 
            this.para = para; 
		}
	}

	
	public static class K2P implements EvolutionParameters
	{		  
        private double kappa;
        
        public K2P(double kappa)
        {
        	 this.kappa = kappa; 
        }
        
		@Override
		public double[] getParameters() {
			return new double[] {kappa};
		}

		@Override
		public void setParameters(double[] para) { 
            this.kappa = para[0]; 
		}
	}

	
}


