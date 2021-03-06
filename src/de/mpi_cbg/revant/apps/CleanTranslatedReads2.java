package de.mpi_cbg.revant.apps;

import java.io.*;

/**
 * Removes from the alphabet all rare characters, and adds to the alphabet all the new
 * unique characters induced by the removal.
 */
public class CleanTranslatedReads2 {
	/**
	 * @param args 5: keep periodic characters untouched (1/0).
	 */
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String ALPHABET_COUNTS_FILE = args[1];
		final int N_NEW_CHARACTERS = Integer.parseInt(args[2]);
		final String NEW_CHARACTERS_FILE = args[3];
		final int MIN_FREQUENCY = Integer.parseInt(args[4]);
		final boolean KEEP_PERIODIC = Integer.parseInt(args[5])==1;
		final String NEW_UNIQUE_FILE = args[6];  // Sorted in decreasing order
		final String NEW_ALPHABET_FILE = args[7];
		final String OLD2NEW_FILE = args[8];  // Map oldNonUnique -> newNonUnique
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.loadAlphabetCount(ALPHABET_COUNTS_FILE,RepeatAlphabet.lastAlphabet+1);
		int[] old2new = RepeatAlphabet.cleanTranslatedRead_updateAlphabet(N_NEW_CHARACTERS,NEW_CHARACTERS_FILE,MIN_FREQUENCY,KEEP_PERIODIC);
		
		BufferedReader br = new BufferedReader(new FileReader(NEW_UNIQUE_FILE));
		RepeatAlphabet.maxOpenLength_unique=Integer.parseInt(br.readLine());
		br.close();
		
		RepeatAlphabet.serializeAlphabet(NEW_ALPHABET_FILE);
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(OLD2NEW_FILE));
		for (int i=0; i<old2new.length; i++) bw.write(old2new[i]+"\n");
		bw.close();
	}

}