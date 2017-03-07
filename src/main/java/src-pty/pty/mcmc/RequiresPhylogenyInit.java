package pty.mcmc;
import java.io.*;
import java.util.*;

import conifer.Phylogeny;

import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;

@Deprecated
public interface RequiresPhylogenyInit
{
  public void init(Phylogeny phylogeny);
}
