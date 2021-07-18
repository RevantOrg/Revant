package de.mpi_cbg.revant.factorize;

import java.io.*;

/**
 * 
 */
public class DBdump2ReadLengths {
	
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		
		char c;
		int i, p, q;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE));
		i=0; str=br.readLine();
		while (str!=null) {
			c=str.charAt(0);
			if (c!='L') {
				str=br.readLine();
				continue;
			}
			p=str.indexOf(" ");
			p=str.indexOf(" ",p+1);
			q=str.indexOf(" ",p+1);
			System.out.println(Integer.parseInt(str.substring(q+1))-Integer.parseInt(str.substring(p+1,q)));
			str=br.readLine();
		}
		br.close();
	}
	
}
