package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Converts the bitvectors built by $FilterAlignments$ (which are assumed to mark
 * alignments on the broken reads produced by $BreakReads$), to bitvectors that marks the
 * original alignments on the original reads that contained low-quality regions.
 */
public class BreakReads5 {
	
	public static void main(String[] args) throws IOException {
		final String BITVECTOR_FILE_NEW = args[0];
		final String TANDEM_BITVECTOR_FILE_NEW = args[1];
		final String ALIGNMENTS_FILE_NEW = args[2];
		final String NEW2OLD_FILE = args[3];
		final int NREADS_NEW = Integer.parseInt(args[4]);
		final String OLD2NEW_FILE = args[5];
		final int NREADS_OLD = Integer.parseInt(args[6]);
		final String ALIGNMENTS_FILE_OLD = args[7];
		final String READ_LENGTHS_FILE_OLD = args[8];
		final String BITVECTOR_FILE_OLD = args[9];
		final String TANDEM_BITVECTOR_FILE_OLD = args[10];
		
		Reads.nReads=NREADS_OLD;
		Reads.loadReadLengths(READ_LENGTHS_FILE_OLD);
		Reads.breakReads_old2new_deserialize(NREADS_OLD,OLD2NEW_FILE);
		Reads.breakReads_new2old_deserialize(NREADS_NEW,NEW2OLD_FILE);
		RepeatAlphabet.breakReads_translateBitvector(BITVECTOR_FILE_NEW,TANDEM_BITVECTOR_FILE_NEW,ALIGNMENTS_FILE_NEW,BITVECTOR_FILE_OLD,TANDEM_BITVECTOR_FILE_OLD,ALIGNMENTS_FILE_OLD);
	}

}