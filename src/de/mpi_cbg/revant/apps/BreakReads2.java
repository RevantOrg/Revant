package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Translates alignments from the original reads to the broken reads.
 * The alignments file might contain just a chunk of the entire set of alignments.
 */
public class BreakReads2 {
	
	public static void main(String[] args) throws IOException {
		final int N_READS_OLD = Integer.parseInt(args[0]);
		final String READ_LENGTHS_FILE = args[1];
		final String OLD2NEW_FILE = args[2];
		final boolean TRANSLATE_B = Integer.parseInt(args[3])==1;
		final boolean WRITE_HEADER = Integer.parseInt(args[4])==1;
		final String READ_READ_ALIGNMENTS = args[5];
		final String READ_READ_ALIGNMENTS_NEW = args[6];
			
		Reads.nReads=N_READS_OLD;
		Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.breakReads_old2new_deserialize(N_READS_OLD,OLD2NEW_FILE);
		RepeatAlphabet.breakReads_translateAlignments(READ_READ_ALIGNMENTS,TRANSLATE_B,WRITE_HEADER,READ_READ_ALIGNMENTS_NEW);
	}

}