package de.mpi_cbg.revant.factorize;

import java.util.*;
import java.io.*;


public class GetQualities {
	
	/**
	 * 
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		int i;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE),10000);
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)!='I') {
				str=br.readLine();
				continue;
			}
			i=str.indexOf(' '); i=str.indexOf(' ',i+1);
			System.out.println(str.substring(i+1));
			str=br.readLine();
		}
		br.close();
	}

}