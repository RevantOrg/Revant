package de.mpi_cbg.revant.factorize;

import java.io.*;


public class BuildCoverageHistograms {

	
	/**
	 * 
	 */
	public static void main(String[] args) throws IOException {
		String INPUT_FILE = args[0];
		String OUTPUT_FILE = args[1];
		int LONGEST_READ = Integer.parseInt(args[2]);
		int N_READS = Integer.parseInt(args[3]);
		String LENGTHS_FILE = args[4];
		int[] tmpArray = new int[LONGEST_READ];
		
		Reads.nReads=N_READS;
		Reads.loadReadLengths(LENGTHS_FILE);
		System.err.println("Read lengths loaded");
		Alignments.getCoverageHistograms(INPUT_FILE,OUTPUT_FILE,Reads.readLengths,tmpArray);
	}
	

}