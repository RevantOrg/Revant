package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Makes sure that the result of $FixEndBlocks$ is consistent on the two sides of every
 * long low-quality interval if the reads were broken.
 */
public class BreakReads3 {
	
	public static void main(String[] args) throws IOException {
		final boolean MODE = Integer.parseInt(args[0])==1;  // 1=replacement; 0=insertion.
		final int LENGTH_TOLERANCE = Integer.parseInt(args[1]);
		final String NEW2OLD_FILE = args[2];
		final int N_READS_NEW = Integer.parseInt(args[3]);
		final String READ_LENGTHS_FILE_NEW = args[4];
		final String ALPHABET_FILE = args[5];
		final String TRANSLATED_FILE = args[6];
		final String TRANSLATED_FILE_DISAMBIGUATED = args[7];
		final String BOUNDARIES_FILE = args[8];
		final String OUTPUT_FILE = args[9];
		
		Reads.nReads=N_READS_NEW;
		Reads.loadReadLengths(READ_LENGTHS_FILE_NEW);
		Reads.breakReads_new2old_deserialize(N_READS_NEW,NEW2OLD_FILE);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.breakReads_checkDisambiguation(MODE,LENGTH_TOLERANCE,TRANSLATED_FILE,TRANSLATED_FILE_DISAMBIGUATED,BOUNDARIES_FILE,OUTPUT_FILE);
	}

}