package de.mpi_cbg.revant.apps;

import java.io.*;

/**
 * Collects all tandem intervals from every translated read.
 */
public class CollectTandems {
	
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String TRANSLATED_FILE = args[1];
		final String BOUNDARIES_FILE = args[2];
		final String READ_LENGTHS_FILE = args[3];
		final String REPEAT_LENGTHS_FILE = args[4];
		final int N_REPEATS = Integer.parseInt(args[5]);
		final String TANDEMS_FILE = args[6];  // Output file
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		final long nIntervals = RepeatAlphabet.getTandemIntervals(TRANSLATED_FILE,BOUNDARIES_FILE,READ_LENGTHS_FILE,TANDEMS_FILE);
		System.err.println(nIntervals+" tandem intervals found");
	}

}