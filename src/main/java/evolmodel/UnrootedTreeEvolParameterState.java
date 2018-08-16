package evolmodel;

import pty.mcmc.UnrootedTreeState;

public class UnrootedTreeEvolParameterState {
	
	private UnrootedTreeState unrootedTreeState; 
	private EvolutionParameters evolParameter;
	
	public UnrootedTreeEvolParameterState(UnrootedTreeState unrootedTreeState, EvolutionParameters evolParameter)
	{
		this.unrootedTreeState = unrootedTreeState; 
		this.evolParameter = evolParameter;  
	}
	
	public UnrootedTreeState getUnrootedTreeState() {
		return unrootedTreeState;
	}
	public void setUnrootedTreeState(UnrootedTreeState unrootedTreeState) {
		this.unrootedTreeState = unrootedTreeState;
	}
	public EvolutionParameters getEvolParameter() {
		return evolParameter;
	}
	public void setEvolParameter(EvolutionParameters evolParameter) {
		this.evolParameter = evolParameter;
	} 
	
	

}
