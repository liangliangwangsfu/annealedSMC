package pty.smc.models;
import java.io.*;
import java.util.*;


import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class CachedEigenDecomp implements Serializable
{
  private static final long serialVersionUID = 1L;
  private final Matrix V, Vinv, D;
  private final double [] imag,real;
  public CachedEigenDecomp(EigenvalueDecomposition ed)
  {
    V = ed.getV();
    Vinv =  V.inverse();
    D = ed.getD();
    this.imag = ed.getImagEigenvalues();
    this.real = ed.getRealEigenvalues();
  }
  public Matrix getV()
  {
    return V;
  }
  public Matrix getVinv()
  {
    return Vinv;
  }
  public Matrix getD()
  {
    return D;
  }
  public double [] getImagEigenvalues()
  {
    return imag;
  }
  public double [] getRealEigenvalues()
  {
    return real;
  }
}
