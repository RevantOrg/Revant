package de.mpi_cbg.revant.factorize;

import java.io.*;

/**
 * Transforms the quality file produced by DBdump into a file containing just the sequence 
 * of phred scores of every read.
 */
public class DBdump2Phred {
	
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		
		char c;
		int i, p;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE));
		i=0; str=br.readLine();
		while (str!=null) {
			c=str.charAt(0);
			if (c!='I') {
				str=br.readLine();
				continue;
			}
			p=str.indexOf(" ");
			p=str.indexOf(" ",p+1);
			System.out.println(str.substring(p+1));
			str=br.readLine();
		}
		br.close();
	}
	
}
