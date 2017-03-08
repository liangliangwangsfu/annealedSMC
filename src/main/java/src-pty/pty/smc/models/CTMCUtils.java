package pty.smc.models;
import java.io.*;
import java.util.*;

import nuts.io.IO;
import nuts.lang.StringUtils;
import nuts.math.RateMtxUtils;

import fig.basic.IOUtils;
import fig.exec.Execution;
import goblin.HLParams;

public class CTMCUtils
{
  public static void saveInExec(CTMC ctmc, String fileNamePrefix)
  {
    String seriFile = Execution.getFile(fileNamePrefix + ".CTMC"),
            txtFile = Execution.getFile(fileNamePrefix + ".CTMC.txt");
    serialize(ctmc, new File(seriFile));
    IO.writeToDisk(txtFile, ctmc.toString());
  }
  public static void serialize(CTMC ctmc, File file)
  {
    try
    {
      ObjectOutputStream oos = IOUtils.openBinOut(
          file);
      oos.writeObject(ctmc);
      oos.close();
    } 
    catch (IOException e) { throw new RuntimeException(e); }
  }
  public static CTMC unSerialize(File file)
  {
    try
    {
      ObjectInputStream ois = IOUtils.openBinIn(file);
      return (CTMC) ois.readObject();
    } catch (Exception e) { throw new RuntimeException(e); }
  }
}
