package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Breaks reads at low-quality regions, and builds the data structures for translating
 * alignments from the original reads to the broken reads and vice versa.
 * The number of new reads is printed to STDOUT.
 */
public class BreakReads1 {
		
	public static void main(String[] args) throws IOException {
		final int N_READS = Integer.parseInt(args[0]);
		final String READ_LENGTHS_FILE = args[1];
		final String QUALITIES_FILE = args[2];
		final String QUALITY_THRESHOLDS_FILE = args[3];
		final int MIN_LOW_QUALITY_LENGTH = Integer.parseInt(args[4]);
		final String OLD2NEW_FILE = args[5];
		final String NEW2OLD_FILE = args[6];
		final String READ_LENGTHS_FILE_NEW = args[7];
		
		final int nReads_new = Reads.breakReads(null,N_READS,READ_LENGTHS_FILE,0,QUALITIES_FILE,QUALITY_THRESHOLDS_FILE,MIN_LOW_QUALITY_LENGTH,OLD2NEW_FILE,null);		
		Reads.breakReads_serialize(N_READS,nReads_new,null,NEW2OLD_FILE,READ_LENGTHS_FILE_NEW);
		System.out.println(nReads_new+"");
	}

}