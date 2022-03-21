package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Collects all distinct k-mers from a chunk of translated reads, and prints them sorted
 * in output. If a list of intervals is provided, the program avoids every window of a 
 * sequence of blocks that contains an interval.
 */
public class CollectKmers {
	
	public static void main(String[] args) throws IOException {
		final int K = Integer.parseInt(args[0]);
		final String TRANSLATED_FILE = args[1];
		final String ALPHABET_FILE = args[2];
		final int UNIQUE_MODE = Integer.parseInt(args[3]);  // See $RepeatAlphabet.isValidWindow()$
		final int MULTI_MODE = Integer.parseInt(args[4]);
		final String AVOIDED_INTERVALS_FILE = args[5];  // NULL to discard it
		final String KMERS_FILE = args[6];  // Output file
		
		boolean INTERVALS_FILE_EXISTS = !AVOIDED_INTERVALS_FILE.equalsIgnoreCase("null");
		
		int i;
		int row, nKmers, lastAvoidedInterval;
		String str1, str2;
		BufferedReader br1, br2;
		BufferedWriter bw;
		RepeatAlphabet.Kmer tmpKmer = new RepeatAlphabet.Kmer();
		HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer> kmers, tmpMap;
		int[] avoidedIntervals;
		int[] tmpArray2 = new int[K];
		int[] tmpArray3 = new int[(K)<<1];
		String[] tokens;
		RepeatAlphabet.Kmer[] keys;
		
		// Collecting k-mers
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.kmerPool_init(K);
		kmers = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
		tmpMap = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
		br1 = new BufferedReader(new FileReader(TRANSLATED_FILE));
		if (INTERVALS_FILE_EXISTS) {
			br2 = new BufferedReader(new FileReader(AVOIDED_INTERVALS_FILE));
			avoidedIntervals = new int[10];  // Arbitrary, multiple of 2.
		}
		else { br2=null; avoidedIntervals=null; }
		str1=br1.readLine(); str2=INTERVALS_FILE_EXISTS?br2.readLine():null; row=0;
		while (str1!=null) {
			if (row%100000==0) System.err.println("Processed "+row+" reads, "+kmers.size()+" distinct "+K+"-mers.");
			if (INTERVALS_FILE_EXISTS) {
				if (str2.length()==0) lastAvoidedInterval=-1;
				else {
					tokens=str2.split(","); lastAvoidedInterval=tokens.length-1;
					if (lastAvoidedInterval>=avoidedIntervals.length) avoidedIntervals = new int[lastAvoidedInterval+1];
					for (i=0; i<=lastAvoidedInterval; i++) avoidedIntervals[i]=Integer.parseInt(tokens[i]);
				}
			}
			else lastAvoidedInterval=-1;
			RepeatAlphabet.getKmers(str1,K,UNIQUE_MODE,MULTI_MODE,kmers,null,avoidedIntervals,lastAvoidedInterval,-1,tmpKmer,tmpArray2,tmpArray3,tmpMap);
			str1=br1.readLine(); str2=INTERVALS_FILE_EXISTS?br2.readLine():null; row++;
		}
		br1.close(); nKmers=kmers.size();
		if (INTERVALS_FILE_EXISTS) br2.close();
		System.err.println(nKmers+" total distinct "+K+"-mers (uniqueMode="+UNIQUE_MODE+", multiMode="+MULTI_MODE+")");
		
		// Serializing sorted k-mers
		keys = new RepeatAlphabet.Kmer[nKmers];
		kmers.keySet().toArray(keys);
		if (nKmers>1) Arrays.sort(keys,0,nKmers);
		bw = new BufferedWriter(new FileWriter(KMERS_FILE));
		for (i=0; i<nKmers; i++) bw.write(keys[i].toString()+","+keys[i].count+","+keys[i].sameReadCount+"\n");
		bw.close();
	}

}