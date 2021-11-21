package de.mpi_cbg.revant.apps;

import java.io.*;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.util.IO;

/**
 * Writes to the output instructions on which alignments to keep or to remove.
 * Loose filtering (MODE=0) discards only alignments such that the sequence of both 
 * intervals is likely to occur multiple times in the genome. Tight filtering (MODE=1) 
 * keeps only alignments such that the sequence of both intervals is likely to occur just 
 * once in the genome. MODE=2 additionally requires the intervals to cover matching 
 * characters of the repeat alphabet.
 */
public class FilterAlignments {
	
	public static void main(String[] args) throws IOException {
		final String ALIGNMENTS_FILE = args[0];
		final int N_READS = Integer.parseInt(args[1]);
		final String READ_LENGTHS_FILE = args[2];
		final String READ_IDS_FILE = args[3];
		final String TRANSLATED_READS_FILE = args[4];
		final String TRANSLATED_BOUNDARIES_FILE = args[5];
		final String FULLY_UNIQUE_FILE = args[6];
		final int N_FULLY_UNIQUE = Integer.parseInt(args[7]);
		final String FULLY_CONTAINED_FILE = args[8];
		final int N_FULLY_CONTAINED = Integer.parseInt(args[9]);
		final String UNIQUE_INTERVALS_FILE = args[10];
		final int MODE = Integer.parseInt(args[11]);
		final String ALPHABET_FILE = args[12];  // Discarded if MODE==0
		final String OUTPUT_FILE = args[13];
		
		final int MIN_INTERSECTION = IO.quantum<<1;  // Arbitrary
		
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		RepeatAlphabet.loadReadsFully(FULLY_UNIQUE_FILE,N_FULLY_UNIQUE,FULLY_CONTAINED_FILE,N_FULLY_CONTAINED);
		RepeatAlphabet.loadBlueIntervals(UNIQUE_INTERVALS_FILE);
		if (MODE==0) {
			RepeatAlphabet.loadAllBoundaries(TRANSLATED_READS_FILE,true,false,TRANSLATED_BOUNDARIES_FILE);
			RepeatAlphabet.filterAlignments_loose(ALIGNMENTS_FILE,OUTPUT_FILE,MIN_INTERSECTION);
		}
		else {
			RepeatAlphabet.loadAllBoundaries(TRANSLATED_READS_FILE,false,true,TRANSLATED_BOUNDARIES_FILE);
			RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
			RepeatAlphabet.filterAlignments_tight(ALIGNMENTS_FILE,OUTPUT_FILE,MODE==1?false:true,MIN_INTERSECTION);
		}
	}

}