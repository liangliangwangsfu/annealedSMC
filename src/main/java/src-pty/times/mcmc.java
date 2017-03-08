package times;


import fig.basic.Pair;

import java.util.Random;


public class mcmc {
	State state;
	Random rand;
	int iteration ;
	
	public mcmc (State state, Random rand) {
		this.state = state;
		this.rand = rand;
		this.iteration = 0 ;
	}
	
	public State sample ()  {
		Pair<State, Double> proposed = state.propose (rand);
		State newState = proposed.getFirst();
		double q = proposed.getSecond ();
		q += newState.getLogDensity() - state.getLogDensity();
		q = Math.exp(q);
		q  = q < 1? q:1;
		if ( rand.nextDouble() < q ) 
			state = newState;
		iteration ++;
		return state;
	}
	
}