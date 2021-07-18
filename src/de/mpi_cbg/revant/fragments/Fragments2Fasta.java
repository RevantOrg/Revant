package de.mpi_cbg.revant.fragments;

import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class Fragments2Fasta {
	
	/**
	 * Transforms a fragments file into a FASTA file to be fed to $fasta2DB$.
	 */
	public static void main(String[] args) throws IOException {
		final String FRAGMENTS_FILE = args[0];
		final String DATABASE_FILE = args[1];
		final String DBSHOW_COMMAND = args[2];
		final String OUTPUT_FILE = args[3];
		final String SEPARATOR = ",";
		int read;
		String str;
		BufferedReader br, brPrime;
		BufferedWriter bw;
		Process proc;
		String[] tokens;
		
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE),IO.BUFFER_SIZE);
		br = new BufferedReader(new FileReader(FRAGMENTS_FILE),IO.BUFFER_SIZE);
		str=br.readLine();
		while (str!=null) {
			tokens=str.split(SEPARATOR);
			read=Integer.parseInt(tokens[0])+1;
			proc=Runtime.getRuntime().exec(DBSHOW_COMMAND+" -w"+Math.POSITIVE_INFINITY+" "+DATABASE_FILE+" "+read);
			brPrime = new BufferedReader(new InputStreamReader(proc.getInputStream()),IO.BUFFER_SIZE);
			bw.write(brPrime.readLine()+"\n");  // FASTA header
			bw.write(brPrime.readLine().substring(Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2])+1)+"\n");
			brPrime.close();
			str=br.readLine();
		}
		br.close(); bw.close();
	}

}