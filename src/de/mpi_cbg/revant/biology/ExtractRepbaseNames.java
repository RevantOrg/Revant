package de.mpi_cbg.revant.biology;

import java.io.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;


public class ExtractRepbaseNames {
	
	/**
	 * 
	 */
	public static void main(String[] args) throws IOException {
		final int MAX_LENGTH_OF_WRONG_LABEL = 16;
		final String INPUT_FILE = args[0];
		int i, j;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE),IO.BUFFER_SIZE);
		str=br.readLine();
		while (str!=null) {
			if (str.length()==0 || str.charAt(0)!='>') {
				str=br.readLine();
				continue;
			}
			i=str.indexOf("\t");
			if (i<0) i=0;
			j=str.indexOf("\t",i+1);
			if (j<0) j=Math.min(str.length(),MAX_LENGTH_OF_WRONG_LABEL);
			System.out.println(str.substring(i+1,j));
			str=br.readLine();
		}
		br.close();
	}

}