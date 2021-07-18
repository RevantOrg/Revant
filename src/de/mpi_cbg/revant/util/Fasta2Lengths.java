package de.mpi_cbg.revant.util;

import java.io.*;

/**
 * Prints the length of every Fasta string in a multi-Fasta file.
 */
public class Fasta2Lengths {
	
	/**
	 * 
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		int length;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE),IO.BUFFER_SIZE);
		str=br.readLine(); length=0;
		while (str!=null) {
			if (str.length()==0 || str.charAt(0)=='>') {
				if (length>0) System.out.println(length);
				length=0;
				str=br.readLine();
				continue;
			}
			length+=str.length();
			str=br.readLine();
		}
		if (length>0) System.out.println(length);
		br.close();
	}

}