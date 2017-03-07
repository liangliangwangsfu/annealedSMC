package pty.smc;

import pty.io.Dataset;
import pty.smc.models.CTMC;

public class PartialCoalescentState4BackForwardKernel extends
		PartialCoalescentState {
	private PartialCoalescentState4BackForwardKernel parent = null;
	private double deltaOld = 0;

	public double getDeltaOld() {
		return deltaOld;
	}

	public void setDeltaOld(double deltaOld) {
		this.deltaOld = deltaOld;
	}

	public PartialCoalescentState4BackForwardKernel(PartialCoalescentState pcs,
			PartialCoalescentState4BackForwardKernel parent, double deltaOld) {
		super(pcs);
		this.parent = parent;
		this.deltaOld = deltaOld;
	}
	
	public static PartialCoalescentState initFastState(Dataset data, CTMC ctmc,
			boolean isClock) {
		PartialCoalescentState.initFastState(false, data, ctmc, isClock);
		return initFastState(false, data, ctmc, isClock);
	}

	public PartialCoalescentState4BackForwardKernel parentState() {
		return parent;
	}

	public void setParent(PartialCoalescentState4BackForwardKernel parent) {
		this.parent = parent;
	}

	public boolean hasParent()
{
		return (this.parent != null);
			}
}
