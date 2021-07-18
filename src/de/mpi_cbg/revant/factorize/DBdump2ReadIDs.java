package de.mpi_cbg.revant.factorize;

import java.io.*;

/**
 * Extracts just the read IDs and makes them zero-based.
 */
public class DBdump2ReadIDs {
	
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		
		char c;
		int p;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE));
		str=br.readLine();
		while (str!=null) {
			c=str.charAt(0);
			if (c!='R') {
				str=br.readLine();
				continue;
			}
			System.out.println((Integer.parseInt(str.substring(2))-1)+"");
			str=br.readLine();
		}
		br.close();
	}
	
}
