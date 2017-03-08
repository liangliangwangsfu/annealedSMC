package times;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import nuts.io.IO;
import nuts.util.Arbre;
import nuts.util.Tree;
import ma.newick.NewickParser;

public class GetClades {
	
	private HashMap<Set<String>, String > cladesMap;

	public GetClades (String file)  {
		try  {
			String treeString =  IO.f2s(file);
			NewickParser np = new NewickParser (treeString);
			Tree<String> tree= np.parse();
			Arbre<String> ar = Arbre.tree2Arbre(tree);
			Map <Arbre<String>, Set<String>> leavesMap = Arbre.leavesMap(ar);
			cladesMap = new HashMap < Set<String>, String> ();
		
		for (Arbre<String> a:leavesMap.keySet()){
			String label = a.getContents();
			cladesMap.put(leavesMap.get(a),label);
		}
		
		BufferedWriter bw = new BufferedWriter (new FileWriter ("clades.txt"));

		for (Set<String> a:cladesMap.keySet()) {
			
			String s  = a.toString();
			s = s.replaceAll(" ", "");
			
			bw.write (s + "\t" + cladesMap.get(a) + "\n");
		}
		bw.close ();

		} catch (Exception e) { 
			e.printStackTrace();
		}
		
	}
	
	public static void main (String args[]) {
		new GetClades (args[0]);
	}
}
