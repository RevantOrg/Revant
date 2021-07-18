package de.mpi_cbg.revant.util;

import java.io.*;

public class Fasta2LengthsHeader {


	public static void main(String[] args) throws IOException {
		final String FASTA_FILE = args[0];
		String str;
		
		BufferedReader br = new BufferedReader(new FileReader(FASTA_FILE));
		str=br.readLine();
		while (str!=null) {
			System.out.println(str.substring(str.indexOf("_")+1));
			str=br.readLine();
			str=br.readLine();
		}
		br.close();
	}


}