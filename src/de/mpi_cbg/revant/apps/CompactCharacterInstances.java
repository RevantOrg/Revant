package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;


public class CompactCharacterInstances {
	
	public static void main(String[] args) throws IOException {
		final String INSTANCES_FILE = args[0];
		final int N_REPEATS = Integer.parseInt(args[1]);
		final String REPEAT_LENGTHS_FILE = args[2];
		final String OUTPUT_FILE = args[3];
		
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		RepeatAlphabet.deserializeAlphabet(INSTANCES_FILE,2);
		if (RepeatAlphabet.lastAlphabet>=0) { 
			RepeatAlphabet.compactInstances();
			RepeatAlphabet.closeAlphabetByRC();
			// $compactInstances()$ might be needed again: see $closeAlphabetByRC()$.
			RepeatAlphabet.compactInstances();
		}
		RepeatAlphabet.serializeAlphabet(OUTPUT_FILE);
	}


}