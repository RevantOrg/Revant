package de.mpi_cbg.revant.factorize;

import java.io.*;
import de.mpi_cbg.revant.util.IO;


public class Test {
		
	/**
	 * java Test ../simulations/drosophila/LAshow-6.txt ../simulations/drosophila/connection-6.txt
	 * java Test ../simulations/drosophila/LAshow-all.txt ../simulations/drosophila/connection-all.txt
	 */
	public static void main(String[] args) throws IOException {
		final BufferedReader alignmentsFile = new BufferedReader(new FileReader(args[0]),IO.BUFFER_SIZE);
		final BufferedReader connectionFile = new BufferedReader(new FileReader(args[1]),IO.BUFFER_SIZE);
		int line;
		String str1, str2;
		
		str1=alignmentsFile.readLine();
		str1=alignmentsFile.readLine();  // Skipping the first two lines
		str1=alignmentsFile.readLine();
		str2=connectionFile.readLine();
		line=3;
		while (str1!=null) {
			if (str2.length()==0) {
				str1=alignmentsFile.readLine();
				str2=connectionFile.readLine();
				line++;
				continue;
			}
			Alignments.readAlignmentFile(str1);
			Factorize.readIntervalConnectionFile(str2);
			if ( Alignments.startA<=Factorize.start-500 || 
			     Alignments.endA>=Factorize.end+500 ||
				 Alignments.endA<=Factorize.start-100 ||
				 Alignments.startA>=Factorize.end+100
			   ) {
				System.err.println("ERROR in line "+line+" (starting from one) of the connection file. ALIGNMENT: startA="+Alignments.startA+" endA="+Alignments.endA+" start="+Factorize.start+" end="+Factorize.end+" implied="+Factorize.implied);
				System.err.println(str1);
				System.err.println("CONNECTION FILE:");
				System.err.println(str2);
			}
			str1=alignmentsFile.readLine();
			str2=connectionFile.readLine();
			line++;
		}
	}
	
}