package times;

import java.util.*;
import fig.basic.Pair;

public interface State {
	public double getLogDensity();
	public Pair<State, Double>	propose (Random rand);
}
