package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Reads;

/**
 * Collects all the non-repetitive blocks that are adjacent to a tandem, and tries to  
 * assign a repeat character to them. The program prints to STDOUT: 1 if there is no
 * tandem spacer; 2 if there are tandem spacers but none of them aligns to a repeat;
 * 3 if most tandem spacers do not have a solution after propagation; 0 otherwise.
 *
 * Remark: this program processes all reads and all alignments sequentially. This is done
 * just for simplicity. Collecting tandem spacers could be done in parallel on different
 * chunks of reads, but propagating solutions needs all tandem spacers in memory.
 */
public class FixTandemSpacers1 {
	
	/**
	 * @param args 12 TRUE=use all non-repetitive blocks in a read.
	 */
	public static void main(String[] args) throws IOException {
		final int N_READS = Integer.parseInt(args[0]);
		final String READ_IDS_FILE = args[1];
		final String READ_LENGTHS_FILE = args[2];
		final String ALPHABET_FILE = args[3];
		final String TRANSLATED_READS_CHARACTERS_FILE = args[4];
		final String TRANSLATED_READS_BOUNDARIES_FILE = args[5];
		final String FULLY_UNIQUE_FILE = args[6];
		final int N_FULLY_UNIQUE = Integer.parseInt(args[7]);
		final String FULLY_CONTAINED_FILE = args[8];
		final int N_FULLY_CONTAINED = Integer.parseInt(args[9]);
		final String READ_READ_ALIGNMENTS_FILE = args[10];
		final String TANDEMS_FILE = args[11];
		final int NONREPETITIVE_BLOCKS_MODE = Integer.parseInt(args[12]);
		final String OUTPUT_FILE = args[13];
		
		final int DISTANCE_THRESHOLD = IO.quantum;
        final int DISTANCE_THRESHOLD_CONSISTENCY = IO.quantum<<2;  // Arbitrary. Using a large threshold for consistency is ok, since it is just checking whether a spacer is assigned very different substrings of the same repeat.
		
		int maxBlockLength;
		int[] tmpArray = new int[2];
		
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		maxBlockLength=RepeatAlphabet.loadAllBoundaries(TRANSLATED_READS_CHARACTERS_FILE,true,true,TRANSLATED_READS_BOUNDARIES_FILE);
		RepeatAlphabet.loadReadsFully(FULLY_UNIQUE_FILE,N_FULLY_UNIQUE,FULLY_CONTAINED_FILE,N_FULLY_CONTAINED);
		RepeatAlphabet.loadTandemIntervals(TANDEMS_FILE,N_READS);
		RepeatAlphabet.loadTandemSpacers(NONREPETITIVE_BLOCKS_MODE);
		RepeatAlphabet.loadTandemSpacers_blocks(READ_READ_ALIGNMENTS_FILE,DISTANCE_THRESHOLD,NONREPETITIVE_BLOCKS_MODE,tmpArray);
		if (RepeatAlphabet.lastSpacer==-1) { System.out.println("1"); return; }
		RepeatAlphabet.loadFullyContainedTranslation(TRANSLATED_READS_CHARACTERS_FILE,N_FULLY_CONTAINED);
		if (RepeatAlphabet.loadTandemSpacerNeighbors(READ_READ_ALIGNMENTS_FILE,NONREPETITIVE_BLOCKS_MODE,tmpArray)==0) { System.out.println("2"); return; }
		if (!RepeatAlphabet.propagateSolutions(DISTANCE_THRESHOLD_CONSISTENCY)) { System.out.println("3"); return; }
		RepeatAlphabet.serializeSpacers(OUTPUT_FILE);
		System.out.println("0");
	}

}