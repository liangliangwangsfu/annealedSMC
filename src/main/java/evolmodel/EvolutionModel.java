package evolmodel;

import pty.smc.models.CTMC;

public enum EvolutionModel {

	GTRIGammaCTMC{

		@Override
		public CTMC instantiateCTMC(EvolutionParameters evolPara, int nSites) {
			// para include: six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT, stationary state frequencies. pi_A, pi_C, pi_G, pi_T, shape parameter in the Gamma distribution,  the proportion of invariant sites
			double[] para = evolPara.getParameters();
			double[] subsRates = new double[]{para[0],para[1],para[2],para[3],para[4],para[5]};     
			double[] statFreqs = new double[]{para[6],para[7],para[8],para[9]};         // stationary state frequencies. pi_A, pi_C, pi_G, pi_T
			double alpha = para[10];           // shape parameter in the Gamma distribution
			double pInv = para[11];            // the proportion of invariant sites 
			return  new CTMC.GTRIGammaCTMC(statFreqs,subsRates,4,nSites,alpha,4,pInv);
		}					
	},
	K2P{

		@Override
		public CTMC instantiateCTMC(EvolutionParameters evolPara, int nSites) {
			double[] para = evolPara.getParameters();
			double csmc_trans2tranv = para[0];						
			return CTMC.SimpleCTMC.dnaCTMC(nSites,csmc_trans2tranv);
		}							
	};

	public abstract CTMC instantiateCTMC(EvolutionParameters evolPara, int nSites);
}


