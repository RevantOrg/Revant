package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.HashMap;
import java.util.Set;
import java.util.Iterator;

/**
 * Prints the unsorted list of all distinct (k-1)-mers from a given list of k-mers.
 * See $RepeatAlphabet.getKMinusOneMers()$ for details. Output format:
 * $sequence,count,previousCharacter,nextCharacter$.
 */
public class GetKMinusOneMers {
	
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String KMERS_FILE = args[1];
		final int K = Integer.parseInt(args[2]);
		final String K_MINUS_ONE_MERS_FILE = args[3];
		
		String str;
		RepeatAlphabet.Kmer tmpKmer;
		Iterator<RepeatAlphabet.Kmer> iterator;
		HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer> kmers, kMinusOneMers;
		BufferedReader br;
		BufferedWriter bw;
		int[] tmpArray;
		
		// Loading k-mers and alphabet
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		kmers = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
		br = new BufferedReader(new FileReader(KMERS_FILE));
		str=br.readLine();
		while (str!=null) {
			tmpKmer = new RepeatAlphabet.Kmer(str,K);
			kmers.put(tmpKmer,tmpKmer);
			str=br.readLine();
		}
		br.close();
		
		// Building (k-1)-mers
		tmpKmer = new RepeatAlphabet.Kmer(); tmpKmer.sequence = new int[K-1];
		tmpArray = new int[(K-1)<<1];
		kMinusOneMers = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>(kmers.size());
		RepeatAlphabet.getKMinusOneMers(kmers,kMinusOneMers,K,tmpKmer,tmpArray);
		kmers.clear(); kmers=null;
		
		// Printing the unsorted output
		bw = new BufferedWriter(new FileWriter(K_MINUS_ONE_MERS_FILE));
		iterator=kMinusOneMers.keySet().iterator();
		while (iterator.hasNext()) {
			tmpKmer=iterator.next();
			bw.write(tmpKmer.toString()+","+tmpKmer.count+","+tmpKmer.previousCharacter+","+tmpKmer.nextCharacter+"\n");
		}
		bw.close();
	}

}