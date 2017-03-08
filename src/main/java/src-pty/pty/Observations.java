package pty;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import nuts.util.Arbre;
import pty.smc.models.CTMC;

public interface Observations extends ObservationDimensions
{
  Map<Taxon,double[][]> observations();
}
