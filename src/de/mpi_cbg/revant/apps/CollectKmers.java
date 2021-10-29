package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Collects all distinct k-mers from a chunk of translated reads, and prints them sorted
 * in output.
 */
public class CollectKmers {
	
	public static void main(String[] args) throws IOException {
		final int K = Integer.parseInt(args[0]);
		final String TRANSLATED_FILE = args[1];
		final String ALPHABET_FILE = args[2];
		final int UNIQUE_MODE = Integer.parseInt(args[3]);  // See $RepeatAlphabet.isValidWindow()$
		final int OPEN_MODE = Integer.parseInt(args[4]);
		final int MULTI_MODE = Integer.parseInt(args[5]);
		final String OUTPUT_FILE = args[6];
		
		int i;
		int row, nKmers;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		RepeatAlphabet.Kmer tmpKmer = new RepeatAlphabet.Kmer();
		HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer> kmers;
		int[] tmpArray2 = new int[K];
		int[] tmpArray3 = new int[(K)<<1];
		RepeatAlphabet.Kmer[] keys;
		
		// Collecting k-mers
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		kmers = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
		br = new BufferedReader(new FileReader(TRANSLATED_FILE));
		str=br.readLine(); row=0;
		while (str!=null) {
			if (row%100000==0) System.err.println("Processed "+row+" reads, "+kmers.size()+" distinct "+K+"-mers.");
			RepeatAlphabet.getKmers(str,K,UNIQUE_MODE,OPEN_MODE,MULTI_MODE,kmers,tmpKmer,tmpArray2,tmpArray3);
			str=br.readLine(); row++;
		}
		br.close(); nKmers=kmers.size();
		System.err.println(nKmers+" total distinct "+K+"-mers (uniqueMode="+UNIQUE_MODE+", openMode="+OPEN_MODE+", multiMode="+MULTI_MODE+")");
		
		// Serializing sorted k-mers
		keys = new RepeatAlphabet.Kmer[nKmers];
		kmers.keySet().toArray(keys);
		if (nKmers>1) Arrays.sort(keys,0,nKmers);
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		for (i=0; i<nKmers; i++) bw.write(keys[i].toString()+","+keys[i].count+"\n");
		bw.close();
	}

}