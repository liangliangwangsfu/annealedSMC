package pty.smc.models;
import java.io.*;
import java.util.*;

import nuts.util.CollUtils;

public class ProductModel implements LikelihoodModelCalculator
{
  private final List<LikelihoodModelCalculator> subModels;
  public ProductModel(List<LikelihoodModelCalculator> subModels)
  {
    this.subModels = CollUtils.archive(subModels);
  }
  public LikelihoodModelCalculator combine(LikelihoodModelCalculator _node1,
      LikelihoodModelCalculator _node2, double delta1, double delta2, boolean doNotBuildCache)
  {
    ProductModel
      node1 = (ProductModel) _node1,
      node2 = (ProductModel) _node2;
    List<LikelihoodModelCalculator> subModels = CollUtils.list();
    final int size = node1.subModels.size();
    if (size != node2.subModels.size()) 
      throw new RuntimeException();
    for (int i = 0; i < size; i++)
      subModels.add(node1.subModels.get(i).combine(
          node1.subModels.get(i), node2.subModels.get(i),
          delta1, delta2, doNotBuildCache));
    return new ProductModel(subModels);
  }
  public double extendLogLikelihood(double delta)
  {
    double sum = 0.0;
    for (LikelihoodModelCalculator lmc : subModels)
      sum += lmc.extendLogLikelihood(delta);
    return sum;
  }
  public boolean isReversible()
  {
    for (LikelihoodModelCalculator lmc : subModels)
      if (!lmc.isReversible())
        return false;
    return true;
  }
  public double logLikelihood()
  {
    double sum = 0.0;
    for (LikelihoodModelCalculator lmc : subModels)
      sum += lmc.logLikelihood();
    return sum;
  }
  public double peekCoalescedLogLikelihood(LikelihoodModelCalculator _node1,
      LikelihoodModelCalculator _node2, double delta1, double delta2)
  {
    double sum = 0.0;
    ProductModel
      node1 = (ProductModel) _node1,
      node2 = (ProductModel) _node2;
    final int size = node1.subModels.size();
    if (size != node2.subModels.size()) 
      throw new RuntimeException();
    for (int i = 0; i < size; i++)
      sum += node1.subModels.get(i).peekCoalescedLogLikelihood(
          node1.subModels.get(i), node2.subModels.get(i),
          delta1, delta2);
    return sum;
  }
}
