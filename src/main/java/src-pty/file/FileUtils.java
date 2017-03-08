package file;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;


public class FileUtils {

	
	public static HashMap<String,String> readMap (String file) throws FileNotFoundException , IOException {
		HashMap <String, String> tmp = new HashMap <String, String> ();
		
		
			BufferedReader br = new BufferedReader (new FileReader (file));
			String s;
			while  ( (s = br.readLine()) !=null) {
				StringTokenizer tok = new StringTokenizer (s);
				String name = tok.nextToken();
				String val = tok.nextToken();
				tmp.put(name, val);
				
			}
		
		return tmp;
	}
}
