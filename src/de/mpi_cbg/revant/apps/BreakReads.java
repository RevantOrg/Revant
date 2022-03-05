package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Assume that the input reads are CLR and that they contain long low-quality regions. The 
 * program breaks each read at every one of its long low-quality regions, and translates
 * read-read and read-repeat alignments accordingly. This transformation allows handling 
 * low-quality regions correctly while translating the reads in the alphabet of repeats:
 * specifically, the whole translation pipeline can be run as a black box on the broken
 * reads.
 *
 * Remark: assume that we just run the repeat translation pipeline as-is on reads with 
 * long low-quality regions. If a low-quality region $L$ falls inside a repetitive region
 * $RLS$, the repetitive substrings $R$ and $S$ would be modeled as closed characters and,
 * since low-quality regions are randomly distributed, such characters would likely have
 * low frequency in the dataset and they would likely be filtered out of the alphabet. 
 * $R$ and $S$ would thus be converted into non-repetitive substrings, and some repeat-
 * induced alignments that contain them might be preserved because of this.
 * 
 * Remark: one might think of modeling $R,S$ above as open characters. Open characters 
 * might be compressed away during alphabet construction, so $R,S$ might get replaced with 
 * closed characters of much bigger length. The translated read might then become $UVW$,
 * where $U,V$ are closed characters (possibly the same closed character) with $|U|>|R|$,
 * $|W|>|S|$, and where $V$ is non-repetitive. Any k-mer built on such a translation would 
 * not reflect the real structure of the read. One might assume that such a k-mer would 
 * likely be globally infrequent and would thus be removed, but it might be kept if it
 * happens to be frequent.
 *
 * Remark: breaking the read at $L$ above would still allow to model $R,S$ as half-open,
 * while at the same time taking advantage of the existing determinization procedures if 
 * $R,S$ happen to match several repeats, and forbidding any k-mer to contain $L$. 
 *
 * Remark: breaking reads is acceptable in this context, since we don't need their
 * contiguity information (as we do in assembly): we just need to translate reads in the 
 * alphabet of repeats and to mark unique substrings in the recoded alphabet, and such 
 * unique substrings cannot contain a low-quality region anyway.
 * 
 * Remark: the program prints the number of new reads to STDOUT.
 */
public class BreakReads {
		
	public static void main(String[] args) throws IOException {
		final int N_READS = Integer.parseInt(args[0]);
		final String READ_LENGTHS_FILE = args[1];
		final String QUALITIES_FILE = args[2];
		final String QUALITY_THRESHOLDS_FILE = args[3];
		final int MIN_LOW_QUALITY_LENGTH = Integer.parseInt(args[4]);
		final String READ_READ_ALIGNMENTS = args[5];
		final String READ_REPEAT_ALIGNMENTS = args[6];
		final String OLD2NEW_FILE = args[7];
		final String NEW2OLD_FILE = args[8];
		final String READ_LENGTHS_FILE_NEW = args[9];
		final String READ_READ_ALIGNMENTS_NEW = args[10];
		final String READ_REPEAT_ALIGNMENTS_NEW = args[11];
		
		final int nReads_new = Reads.breakReads(null,N_READS,READ_LENGTHS_FILE,0,QUALITIES_FILE,QUALITY_THRESHOLDS_FILE,MIN_LOW_QUALITY_LENGTH,null);		
		Reads.breakReads_serialize(N_READS,nReads_new,OLD2NEW_FILE,NEW2OLD_FILE,READ_LENGTHS_FILE_NEW);
		RepeatAlphabet.breakReads_translateAlignments(READ_READ_ALIGNMENTS,true,READ_READ_ALIGNMENTS_NEW);
		RepeatAlphabet.breakReads_translateAlignments(READ_REPEAT_ALIGNMENTS,false,READ_REPEAT_ALIGNMENTS_NEW);
		System.out.println(nReads_new+"");
	}

}