package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Collects short non-periodic spacers between periodic endpoints, and assigns arbitrary
 * breakpoints to them. The program prints "1" to STDOUT if spacers are already correct
 * and the spacer correction pipeline should stop, "0" otherwise.
 *
 * Remark: this program processes all reads and all alignments sequentially. This is done
 * just for simplicity. Collecting spacers could be done in parallel on different chunks 
 * of reads, but assigning breakpoints needs all spacers in memory.
 */
public class FixPeriodicEndpoints1 {
	
	public static void main(String[] args) throws IOException {
		final int MAX_SPACER_LENGTH = Integer.parseInt(args[0]);
		final int MIN_ALIGNMENT_LENGTH_READ_REPEAT = Integer.parseInt(args[1]);
		final int N_READS = Integer.parseInt(args[2]);
		final String READ_IDS_FILE = args[3];
		final String READ_LENGTHS_FILE = args[4];
		final String ALPHABET_FILE = args[5];
		final String TRANSLATED_READS_CHARACTERS_FILE = args[6];
		final String TRANSLATED_READS_BOUNDARIES_FILE = args[7];
		final String FULLY_UNIQUE_FILE = args[8];
		final int N_FULLY_UNIQUE = Integer.parseInt(args[9]);
		final String FULLY_CONTAINED_FILE = args[10];
		final int N_FULLY_CONTAINED = Integer.parseInt(args[11]);
		final String READ_READ_ALIGNMENTS_FILE = args[12];
		final int AVG_READ_LENGTH = Integer.parseInt(args[13]);
		final long GENOME_LENGTH = Long.parseLong(args[14]);
		final String OUTPUT_FILE = args[15];
		
        final double SIGNIFICANCE_LEVEL = 0.05;  // Conventional
        
		int maxBlockLength;
		int[] tmpArray1, tmpArray2, tmpArray3;
		
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		maxBlockLength=RepeatAlphabet.loadAllBoundaries(TRANSLATED_READS_CHARACTERS_FILE,true,true,TRANSLATED_READS_BOUNDARIES_FILE);
		tmpArray1 = new int[maxBlockLength];
		tmpArray2 = new int[maxBlockLength];
		tmpArray3 = new int[maxBlockLength];
		RepeatAlphabet.loadReadsFully(FULLY_UNIQUE_FILE,N_FULLY_UNIQUE,FULLY_CONTAINED_FILE,N_FULLY_CONTAINED);
		RepeatAlphabet.loadSpacers(MAX_SPACER_LENGTH);
		if (RepeatAlphabet.lastSpacer==-1) { System.out.println("1"); return; }
		RepeatAlphabet.loadSpacerNeighbors(READ_READ_ALIGNMENTS_FILE,MIN_ALIGNMENT_LENGTH_READ_REPEAT,tmpArray1,tmpArray2,tmpArray3);
		if (RepeatAlphabet.getSpacerGraphStatistics(N_READS,AVG_READ_LENGTH,GENOME_LENGTH,SIGNIFICANCE_LEVEL,false)) { System.out.println("1"); return; }        
		RepeatAlphabet.assignBreakpoints();
		RepeatAlphabet.serializeSpacers(OUTPUT_FILE);
		System.out.println("0");
	}

}