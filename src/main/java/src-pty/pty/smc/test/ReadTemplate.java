package pty.smc.test;
import java.io.*;
import java.util.*;

import fig.basic.IOUtils;
import nuts.io.IO;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class ReadTemplate
{

  /**
   * @param args
   */
  public static void main(String[] args)
  {
    Random rand = new Random();
    for (int i = 0; i < 4; i++)
    {
      PrintWriter out = IOUtils.openOutHard(args[0] + ".rand" + i );
      for (String line : IO.i(args[0]))
        out.println(line.replaceAll("rand", "" + rand.nextDouble()));
      out.close();
    }
  }

}
