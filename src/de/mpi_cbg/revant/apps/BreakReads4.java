package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Converts the bitvector built by $FilterAlignments$ (which is assumed to mark alignments
 * on the broken reads produced by $BreakReads$), to a bitvector that marks the original 
 * alignments on the original reads with low-quality regions.
 */
public class BreakReads4 {
	
	public static void main(String[] args) throws IOException {
		final String BITVECTOR_FILE_NEW = args[0];
		final String ALIGNMENTS_FILE_NEW = args[1];
		final String NEW2OLD_FILE = args[2];
		final int NREADS_NEW = Integer.parseInt(args[3]);
		final String OLD2NEW_FILE = args[4];
		final int NREADS_OLD = Integer.parseInt(args[5]);
		final String BITVECTOR_FILE_OLD = args[6];
		final String ALIGNMENTS_FILE_OLD = args[7];
		
		Reads.breakReads_old2new_deserialize(NREADS_OLD,OLD2NEW_FILE);
		Reads.breakReads_new2old_deserialize(NREADS_NEW,NEW2OLD_FILE);
		RepeatAlphabet.breakReads_translateBitvector(BITVECTOR_FILE_NEW,ALIGNMENTS_FILE_NEW,BITVECTOR_FILE_OLD,ALIGNMENTS_FILE_OLD);
	}

}