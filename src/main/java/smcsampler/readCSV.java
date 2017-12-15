package smcsampler;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Random;
import java.util.Scanner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class readCSV {
	
	public static List<Double> DeterministicTem(String fileName){
		File file = new File(fileName);
		List<Double> temp = new ArrayList<Double>();

		try {
			Scanner inputStream = new Scanner(file);
			while(inputStream.hasNext()) {
				String data = inputStream.next();
				temp.add(Double.parseDouble(data));
			}
			inputStream.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return temp;

	}
	



}
