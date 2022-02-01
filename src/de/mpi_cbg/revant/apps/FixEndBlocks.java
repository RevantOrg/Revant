package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Tries to disambiguate the first and the last block of every translated read using the 
 * sorrounding context of length $k$.
 */
public class FixEndBlocks {
	private static final int MAX_AMBIGUITY_HISTOGRAM = 20;  // Arbitrary
	
	/**
	 * Remark: the program writes $X,Y,Z$ to STDOUT, where X is the number of read ends
	 * that have been disambiguated, Y is the max number of read ends that could have been
	 * possibly disambiguated, and Z is the total number of reads.
	 *
	 * @param args 
	 * 1: input file of translated reads with ambiguous starts/ends;
	 * 5: output file of translated reads with some starts/ends disambiguated.
	 */
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String OLD_TRANSLATED_FILE = args[1];
		final String KMERS_FILE = args[2];
		final int K = Integer.parseInt(args[3]);
		final boolean TIGHT_MODE = Integer.parseInt(args[4])==1;
		final String NEW_TRANSLATED_FILE = args[5];
		final String STATS_FILE = args[6];
		
		int i;
		int nReads, nFixed, last;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		RepeatAlphabet.Kmer newKmer, context;
		HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer> kmers;
		int[] out, tmpArray1, tmpArray2, tmpArray3;
		int[][] ambiguityHistogram;
		
		// Loading k-mers and alphabet
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		kmers = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
		br = new BufferedReader(new FileReader(KMERS_FILE));
		str=br.readLine();
		while (str!=null) {
			newKmer = new RepeatAlphabet.Kmer(str,K);
			kmers.put(newKmer,newKmer);
			str=br.readLine();
		}
		br.close();
		
		// Fixing end blocks
		context = new RepeatAlphabet.Kmer();
		context.sequence = new int[K];
		out = new int[] {0,0};
		tmpArray1 = new int[3000/*Arbitrary*/];
		tmpArray2 = new int[K]; tmpArray3 = new int[(K)<<1];
		ambiguityHistogram = new int[2][MAX_AMBIGUITY_HISTOGRAM+1];
		Arrays.fill(ambiguityHistogram[0],0); Arrays.fill(ambiguityHistogram[1],0);
		br = new BufferedReader(new FileReader(OLD_TRANSLATED_FILE));
		bw = new BufferedWriter(new FileWriter(NEW_TRANSLATED_FILE));
		str=br.readLine(); nReads=0;
		while (str!=null) {
			nReads++;
			RepeatAlphabet.fixEndBlocks(str,K,kmers,TIGHT_MODE,context,tmpArray1,tmpArray2,tmpArray3,bw,out,ambiguityHistogram);
			str=br.readLine();
		}
		br.close(); bw.close();
		bw = new BufferedWriter(new FileWriter(STATS_FILE));
		bw.write(out[0]+","+out[1]+","+nReads+"\n");
		bw.close();
		System.err.println("Distribution of endpoints ambiguity:  (endblocks, interior blocks)");
		last=-1;
		for (i=0; i<=MAX_AMBIGUITY_HISTOGRAM; i++) {
			if (ambiguityHistogram[0][i]!=0 || ambiguityHistogram[1][i]!=0) last=i;
		}
		for (i=0; i<=last; i++) System.err.println(i+","+ambiguityHistogram[0][i]+","+ambiguityHistogram[1][i]);
	}

}