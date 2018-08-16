package evolmodel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import evolmodel.EvolutionParameterProposalDistribution.Util.Trans2tranvProposal;
import nuts.math.Sampling;
import pty.UnrootedTree;
import fig.basic.ListUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.prob.Dirichlet;

public interface EvolutionParameterProposalDistribution
{
	/**
	 * 
	 * @param current
	 * @param rand
	 * @return (proposedState, log proposal-ratio), where proposal-ratio =
	 *         Q(old|new)/Q(new|old)
	 */
	public Pair<EvolutionParameters, Double> propose(EvolutionParameters current, Random rand);

	public String description();

	public static class Options {
		@Option 
		public double a_alpha = 1.5; 		
		@Option
		public boolean useK2PProposal = true;
	}


	public static class Util {
		public static final EvolutionParameterProposalDistribution.Options _defaultProposalDistributionOptions = new EvolutionParameterProposalDistribution.Options();

		public static List<EvolutionParameterProposalDistribution> proposalList(Options options,
				UnrootedTree nct, Random rand) {
			List<EvolutionParameterProposalDistribution> result = new ArrayList<EvolutionParameterProposalDistribution>();
			if (options.useK2PProposal) result.add(new K2PProposal(options.a_alpha));		
			Collections.shuffle(result, rand);
			return result;
		}

		
		public static class Trans2tranvProposal{
			private final double a; // higher will have lower accept rate		
		    private Random rand;
			public Trans2tranvProposal(double a, Random rand) {
				if (a <= 1)
					throw new RuntimeException();
				this.a = a;
				this.rand=rand;			
			}

			public  Pair<Double,Double> propose(double currentTrans2tranv) {
				double lambda=2*Math.log(a);
				double rvUnif = Sampling.nextDouble(rand, 0, 1);
				double m  = Math.exp(lambda*(rvUnif-0.5));		
				final double newTrans2tranv = m * currentTrans2tranv;
				return Pair.makePair(newTrans2tranv, Math.log(m));						
			}
		}

		
		
		public static double[] proposeAlpha(Random rand, double alpha, double a_alpha, double low, double high){

			double [] result=new double[2]; 
			double scale=0,proposedAlpha=Double.MAX_VALUE; 		
			while(proposedAlpha<low || proposedAlpha>high){
				scale=Sampling.nextDouble(rand, 1.0/a_alpha, a_alpha);
				proposedAlpha=scale*alpha;		
			}		
			result[0]=scale; 
			result[1]=proposedAlpha;
			return result; 
		}

		public static double proposePInv(Random rand, double pInv, double a_pInv, double low, double high){		
			double proposedPInv=Double.MAX_VALUE; 		
			while(proposedPInv<low || proposedPInv>high){
				proposedPInv=Sampling.nextDouble(rand, Math.max(0, pInv-a_pInv), Math.min(1, pInv+a_pInv));
			}				
			return proposedPInv; 
		}

		public static double[] proposeFromDirichlet(Random rand, double a, double[] rates){		
			double[] alphas = new double[rates.length];
			for (int i = 0; i < rates.length; i++)alphas[i]=a*rates[i];
			double[] result = Dirichlet.sample(rand, alphas); 
			return result;
		}

		public static double logProposal(double scale, double[] rates){
			double[] alphas = new double[rates.length];
			for (int i = 0; i < rates.length; i++)alphas[i]=scale*rates[i];	      	        	    		
			return Dirichlet.logProb(alphas, ListUtils.sum(alphas), rates);
		}

	}


	public static class K2PProposal implements EvolutionParameterProposalDistribution {		
		private final double a; 
		public K2PProposal(double a)
		{
			this.a = a; 
		}
		@Override
		public String description() {
			return "K2PProposal";
		}

		@Override
		public Pair<EvolutionParameters, Double> propose(EvolutionParameters current, Random rand) {
            double currentTrans2tranv = current.getParameters()[0];  
      		Trans2tranvProposal kappaProposal=new Trans2tranvProposal(a,rand);
      		Pair<Double,Double> proposed=kappaProposal.propose(currentTrans2tranv);		
      		double proposedTrans2tranv=proposed.getFirst();		      		
      		return Pair.makePair(new EvolutionParameters.K2P(proposedTrans2tranv), proposed.getSecond());
		}
	}

	

}
