package de.mpi_cbg.revant.factorize;

import java.io.*;

/**
 * Transforms a track file produced by DBdump (e.g. a dust track or a tandem track) into a
 * file containing just the sequence of intervals of every read (or an empty line if a 
 * read contains no interval in the track).
 */
public class DBdump2Track {
	
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		
		char c;
		int i, p;
		int read, lastRead;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE));
		str=br.readLine(); lastRead=-1;
		while (str!=null) {
			c=str.charAt(0);
			if (c!='R') {
				str=br.readLine();
				continue;
			}
			p=str.indexOf(" ");
			read=Integer.parseInt(str.substring(p+1));
			for (i=lastRead+1; i<read-1; i++) System.out.println();
			lastRead=read;
			do {
				str=br.readLine();
				c=str.charAt(0);
			} while (c!='T');
			p=str.indexOf(" "); p=str.indexOf(" ",p+1);
			if (p>=0) System.out.println(str.substring(p+1));
			else System.out.println();
			str=br.readLine();
		}
		br.close();
	}
	
}
