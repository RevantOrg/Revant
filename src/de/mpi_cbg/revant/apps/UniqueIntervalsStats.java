package de.mpi_cbg.revant.apps;

import java.io.*;

/**
 * Prints stats on unique intervals.
 */
public class UniqueIntervalsStats {
	
	public static void main(String[] args) throws IOException {
		final String UNIQUE_INTERVALS_FILE = args[0];
		final String BOUNDARIES_FILE = args[1];
		final String OUTPUT_FILE = args[2];
		
		final int MAX_BLOCKS_PER_READ = 100;  // Arbitrary
		RepeatAlphabet.printUniqueIntervalsStats(UNIQUE_INTERVALS_FILE,BOUNDARIES_FILE,MAX_BLOCKS_PER_READ,OUTPUT_FILE);
	}

}