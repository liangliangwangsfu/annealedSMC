package pty.smc.exp;
import java.io.*;
import java.util.*;

import fig.basic.Pair;

import nuts.io.IO;
import nuts.math.Plot2D;
import nuts.math.Plot2D.Plot2DOptions;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class CreateGraph implements Runnable
{

  /**
   * @param args
   */
  public static void main(String[] args)
  {
    IO.run(args, new CreateGraph());
  }

  @Override
  public void run()
  {
    Plot2DOptions options = new Plot2DOptions();
    options.xlogscale = true;
    String prefix = "graph-data-";
    for (int graphIdx = 1; graphIdx <=2 ; graphIdx++)
    {
      Plot2D cur = new Plot2D(options);
      
      for (String method : Arrays.asList("smc", "mcmc"))
      {
        List<Pair<Double,Double>> curList = list();
        for (String line : IO.i(new File(prefix + graphIdx + ".txt." + method)))
          if (!line.matches("^//s*$"))
          {
            String [] fields = line.split("\\s+");
            curList.add(Pair.makePair(Double.parseDouble(fields[0]), Double.parseDouble(fields[1])));
          }
        cur.addSeries(curList, true, method.toUpperCase()); 
      }
      cur.savePlot(new File("../graphics/graph" + (graphIdx == 1 ? "-med" : "") + ".pdf"));
    }
  }

}
