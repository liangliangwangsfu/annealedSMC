package pty.learn;
import java.io.*;
import java.util.*;

import nuts.math.RateMtxUtils;
import Jama.Matrix;

public class Estimators
{
  public static double[][] getGeneralRateMatrixMLE(Matrix suffStats)
  {
    return getGeneralRateMatrixMLE(suffStats, Double.POSITIVE_INFINITY);
  }
  public static double[][] getGeneralRateMatrixMLE(Matrix suffStats, double maxValue)
  {
    final int nChars = suffStats.getColumnDimension();
    double [][] result = new double[nChars][nChars];
    for (int a = 0; a < nChars; a++)
      for (int b = 0; b < nChars; b++)
        if (a!=b)
          result[a][b] = Math.min(maxValue, suffStats.get(a,b) / suffStats.get(a,a));
    RateMtxUtils.fillRateMatrixDiagonalEntries(result);
    return result;
  }
  
//
//public static double[][] getMLERateMatrix(Matrix suffStats)
//{
//  final int nChars = suffStats.getColumnDimension();
////  double [] statDistn = new double[nChars];
////  for (int i = 0; i < nChars; i++)
////    statDistn[i] = suffStats.get(i,i);
////  NumUtils.normalize(statDistn);
//  double [][] result = new double[nChars][nChars];
//  for (int row = 0; row < nChars; row++)
//  {
//    double sum = 0;
//    for (int col = 0; col < nChars; col++)
//      if (col != row)
//        sum += result[row][col] = statDistn[col] * 
//          ((suffStats.get(row,col)+suffStats.get(col,row)) /
//           (suffStats.get(col,col)*statDistn[row] + suffStats.get(row,row)*statDistn[col]));
//    result[row][row] = -sum;
//  }
//  return result;
//}
}
