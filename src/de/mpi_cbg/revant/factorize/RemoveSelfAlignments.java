package de.mpi_cbg.revant.factorize;

import java.io.*;
import de.mpi_cbg.revant.util.IO;


public class RemoveSelfAlignments {

	
	/**
	 * 
	 */
	public static void main(String[] args) throws IOException {
		String header1, header2, str;
		BufferedReader br;

		br = new BufferedReader(new FileReader(args[0]),IO.BUFFER_SIZE);
		header1=br.readLine(); System.out.println(header1);
		header2=br.readLine(); System.out.println(header2);
		str=br.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			if (Alignments.readA!=Alignments.readB) System.out.println(str);
			str=br.readLine();
		}
		br.close();
	}
	

}