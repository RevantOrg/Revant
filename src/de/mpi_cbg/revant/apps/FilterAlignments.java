package de.mpi_cbg.revant.apps;

import java.io.*;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;

/**
 * Writes to the output instructions on which alignments to keep or to remove.
 * Loose filtering (MODE=0) discards only alignments such that the sequence of both 
 * intervals is likely to occur multiple times in the genome. Tight filtering (MODE=1) 
 * keeps only alignments such that the sequence of both intervals is likely to occur just 
 * once in the genome. MODE=2 additionally requires the intervals to cover matching 
 * characters of the repeat alphabet.
 *
 * Remark: an obvious limitation of character-based alignment filtering is that it cannot 
 * achieve sub-character resolution. E.g. even if a 3-mer is unique in the genome, the
 * filtering would discard an alignment that spans both adjacencies of the 3-mer, without 
 * spanning the whole 3-mer. This is reasonable, since the first/last substrings are 
 * partial and could match other characters (e.g. longer substrings of the same repeat).
 *
 * Remark: this program does not filter by error rate. This basic check should have been
 * already performed upstream.
 *
 * Remark: if there are no low-quality regions, all local alignments (i.e. not suffix-
 * prefix) are likely repeat-induced, so one could discard them upstream and feed just
 * suffix-prefix alignments to this program. Local alignments that are not filtered out by
 * the program might be due to: (1) the repeat database not containing some repeats in the
 * dataset (their occurrences might get tagged as non-repetitive); (2) rare repeat 
 * characters wrongly recoded as non-repetitive characters by our pipeline; (3) repeated 
 * k-mers with borderline frequency considered unique by our pipeline; (4) the heuristics 
 * of the aligner failing to continue an alignment; (4) wrong read correction/patching at 
 * previous stages of the assembly pipeline; (5) chimeric reads.
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
		final int MIN_ALIGNMENT_LENGTH_READ_READ = Integer.parseInt(args[14]);
		
		final int MIN_INTERSECTION = MIN_ALIGNMENT_LENGTH_READ_READ>>2;  // Arbitrary
		long[][] stats;
		
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.loadReadsFully(FULLY_UNIQUE_FILE,N_FULLY_UNIQUE,FULLY_CONTAINED_FILE,N_FULLY_CONTAINED);
		RepeatAlphabet.loadBlueIntervals(UNIQUE_INTERVALS_FILE);
		stats = new long[3][3];
		Math.set(stats,0);
		if (MODE==0) {
			RepeatAlphabet.loadAllBoundaries(TRANSLATED_READS_FILE,true,false,TRANSLATED_BOUNDARIES_FILE);
			RepeatAlphabet.filterAlignments_loose(ALIGNMENTS_FILE,OUTPUT_FILE,MIN_INTERSECTION,stats);
		}
		else {
			RepeatAlphabet.loadAllBoundaries(TRANSLATED_READS_FILE,false,true,TRANSLATED_BOUNDARIES_FILE);
			RepeatAlphabet.filterAlignments_tight(ALIGNMENTS_FILE,OUTPUT_FILE,MODE==1?false:true,MIN_INTERSECTION,stats);
		}
		System.err.println("All alignments:  (input, output)");
		System.err.println("Suffix/prefix overlaps: \t"+stats[0][0]+" ("+stats[2][0]+") -> "+stats[1][0]+" ("+(100*((double)(stats[1][0]-stats[0][0]))/stats[0][0])+"%)");
		System.err.println("Local substring: \t\t"+stats[0][1]+" ("+stats[2][1]+") -> "+stats[1][1]+" ("+(100*((double)(stats[1][1]-stats[0][1]))/stats[0][1])+"%)");
		System.err.println("Full containment/identity: \t"+stats[0][2]+" ("+stats[2][2]+") -> "+stats[1][2]+" ("+(100*((double)(stats[1][2]-stats[0][2]))/stats[0][2])+"%)");
	}

}