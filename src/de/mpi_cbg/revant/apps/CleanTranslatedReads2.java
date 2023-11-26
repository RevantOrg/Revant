package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;

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
        final int N_READS = Integer.parseInt(args[2]);
        final int AVG_READ_LENGTH = Integer.parseInt(args[3]);
        final int SPANNING_BPS = Integer.parseInt(args[4]);
        final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[5]);  // Read-repeat
        final long GENOME_LENGTH = Long.parseLong(args[6]);  // Of one haplotype
        final int N_HAPLOTYPES = Integer.parseInt(args[7]);
		final int N_NEW_CHARACTERS = Integer.parseInt(args[8]);
		final String NEW_CHARACTERS_FILE = args[9];
		final boolean KEEP_PERIODIC = Integer.parseInt(args[10])==1;
		final String NEW_UNIQUE_FILE = args[11];  // Sorted in decreasing order
		final String NEW_ALPHABET_FILE = args[12];
		final String OLD2NEW_FILE = args[13];  // Map oldNonUnique -> newNonUnique
		
        final double SIGNIFICANCE_LEVEL = 0.05;  // Conventional
        final int MIN_MISSING_LENGTH = IO.quantum;  // Arbitrary
        
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.loadAlphabetCount(ALPHABET_COUNTS_FILE,RepeatAlphabet.lastAlphabet+1);
		int[] old2new = RepeatAlphabet.cleanTranslatedRead_updateAlphabet(N_NEW_CHARACTERS,NEW_CHARACTERS_FILE,KEEP_PERIODIC,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,MIN_ALIGNMENT_LENGTH,MIN_MISSING_LENGTH,SIGNIFICANCE_LEVEL);
		
		BufferedReader br = new BufferedReader(new FileReader(NEW_UNIQUE_FILE));
		RepeatAlphabet.maxOpenLength_unique=Integer.parseInt(br.readLine());
		br.close();
		
		RepeatAlphabet.serializeAlphabet(NEW_ALPHABET_FILE);
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(OLD2NEW_FILE));
		for (int i=0; i<old2new.length; i++) bw.write(old2new[i]+"\n");
		bw.close();
	}

}