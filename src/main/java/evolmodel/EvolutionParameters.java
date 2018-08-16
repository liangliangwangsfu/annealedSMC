package evolmodel;

import java.util.Random;

import fig.prob.Dirichlet;
import nuts.math.Sampling;

public interface EvolutionParameters {
	public double[] getParameters();
	public void setParameters(double[] para);	 
	public void sampleFromPrior(Random rand);

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

		@Override
		public void sampleFromPrior(Random rand) {
			double[] subsRates=Dirichlet.sample(rand, new double[]{10,10,10,10,10,10});  
			// six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
			double[] statFreqs=Dirichlet.sample(rand, new double[]{10,10,10,10});  
			// stationary state frequencies. pi_A, pi_C, pi_G, pi_T				
			double alpha=Sampling.nextDouble(rand,0.1,0.9);       // shape parameter in the Gamma distribution
			double pInv=0;        // the proportion of invariant sites
			para = new double[] {subsRates[0],subsRates[1],subsRates[2],subsRates[3],subsRates[4],subsRates[5],statFreqs[0],statFreqs[1],statFreqs[2],statFreqs[3],alpha,pInv};
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

		@Override
		public void sampleFromPrior(Random rand) {
			kappa = Sampling.nextDouble(rand,1.0/5,5);   
			
		}
	}

	
}


