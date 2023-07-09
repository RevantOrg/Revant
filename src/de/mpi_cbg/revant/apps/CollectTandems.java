package de.mpi_cbg.revant.apps;

import java.io.*;

/**
 * Collects all tandem intervals from every translated read.
 */
public class CollectTandems {
	
	public static void main(String[] args) throws IOException {
		final boolean GET_PERIODIC_TANDEMS = Integer.parseInt(args[0])==1;
		final boolean GET_NONPERIODIC_TANDEMS = Integer.parseInt(args[1])==1;
		final int NONPERIODIC_MODE = Integer.parseInt(args[2]);
		final String ALPHABET_FILE = args[3];
		final String TRANSLATED_FILE = args[4];
		final String BOUNDARIES_FILE = args[5];
		final String READ_LENGTHS_FILE = args[6];
		final String REPEAT_LENGTHS_FILE = args[7];
		final int N_REPEATS = Integer.parseInt(args[8]);
		final String TANDEMS_FILE = args[9];  // Output file
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		final long nIntervals = RepeatAlphabet.getTandemIntervals(GET_PERIODIC_TANDEMS,GET_NONPERIODIC_TANDEMS,NONPERIODIC_MODE,TRANSLATED_FILE,BOUNDARIES_FILE,READ_LENGTHS_FILE,TANDEMS_FILE);
		System.err.println(nIntervals+" tandem intervals found");
	}

}