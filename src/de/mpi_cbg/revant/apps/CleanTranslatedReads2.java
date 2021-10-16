package de.mpi_cbg.revant.apps;

import java.io.*;

/**
 * Removes from the alphabet all rare characters, and adds to the alphabet all the new
 * unique characters induced by the removal.
 */
public class CleanTranslatedReads2 {
	
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String ALPHABET_COUNTS_FILE = args[1];
		final int N_NEW_CHARACTERS = Integer.parseInt(args[2]);
		final String NEW_CHARACTERS_FILE = args[3];
		final int MIN_FREQUENCY = Integer.parseInt(args[4]);
		final String NEW_UNIQUE_FILE = args[5];  // Sorted in decreasing order
		final String OUTPUT_FILE = args[6];
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE);
		RepeatAlphabet.loadAlphabetCount(ALPHABET_COUNTS_FILE,RepeatAlphabet.lastAlphabet+1);
		RepeatAlphabet.cleanTranslatedRead_updateAlphabet(N_NEW_CHARACTERS,NEW_CHARACTERS_FILE,MIN_FREQUENCY);
		BufferedReader br = new BufferedReader(new FileReader(NEW_UNIQUE_FILE));
		RepeatAlphabet.maxOpenLength_unique=Integer.parseInt(br.readLine());
		br.close();
		RepeatAlphabet.serializeAlphabet(OUTPUT_FILE);		
	}

}