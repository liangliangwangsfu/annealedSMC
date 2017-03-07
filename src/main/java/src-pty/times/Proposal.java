package times;

import java.util.Random;
import fig.basic.Pair;


public interface Proposal {
	public Pair<State,Double> propose (State current, Random rand);

}
