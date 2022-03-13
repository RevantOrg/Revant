package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Splits the bitvector built by $FilterAlignments$ (which is assumed to mark alignments
 * on the broken reads), and the alignments on the broken reads themselves, into chunks
 * that correspond to those of the alignments file of the unbroken reads.
 */
public class BreakReads4 {
	/**
	 * @param args 7 file containing the last readA (zero-based) in every chunk of
	 * alignments of the unbroken reads.
	 */
	public static void main(String[] args) throws IOException {
		final String ALIGNMENTS_FILE_NEW = args[0];
		final String BITVECTOR_FILE_NEW = args[1].equalsIgnoreCase("null")?null:args[1];
		final int NREADS_OLD = Integer.parseInt(args[2]);
		final String OLD2NEW_FILE = args[3];
		final String ALIGNMENTS_PREFIX_NEW = args[4];
		final String BIVECTOR_PREFIX_NEW = args[5];
		final int N_CHUNKS = Integer.parseInt(args[6]);
		final String LAST_READA_FILE = args[7];
		
		int i;
		String str;
		BufferedReader br;
		int[] lastReadA_old;
		
		lastReadA_old = new int[N_CHUNKS];
		br = new BufferedReader(new FileReader(LAST_READA_FILE));
		str=br.readLine(); i=-1;
		while (str!=null && str.length()!=0) {
			lastReadA_old[++i]=Integer.parseInt(str)+1;  // The file is zero-based
			str=br.readLine();
		}
		br.close();
		Reads.breakReads_old2new_deserialize(NREADS_OLD,OLD2NEW_FILE);
		RepeatAlphabet.breakReads_splitAlignments(lastReadA_old,N_CHUNKS,ALIGNMENTS_FILE_NEW,BITVECTOR_FILE_NEW,ALIGNMENTS_PREFIX_NEW,BIVECTOR_PREFIX_NEW);
	}

}