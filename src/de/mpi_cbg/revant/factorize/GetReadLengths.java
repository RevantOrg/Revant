package de.mpi_cbg.revant.factorize;

import java.util.*;
import java.io.*;


public class GetReadLengths {
	
	/**
	 * 
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		int i, j;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE),10000);
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)!='L') {
				str=br.readLine();
				continue;
			}
			i=str.indexOf(' '); i=str.indexOf(' ',i+1);
			j=str.indexOf(' ',i+1);
			System.out.println((Integer.parseInt(str.substring(j+1))-Integer.parseInt(str.substring(i+1,j)))+"");
			str=br.readLine();
		}
		br.close();
	}

}