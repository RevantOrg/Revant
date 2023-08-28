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
	/**
     * @param args 0: 0 (enumerating) or 1 (counting): see $RepeatAlphabet.getKmers()$.
     */
	public static void main(String[] args) throws IOException {
        final int MODE = Integer.parseInt(args[0]);
		final int K = Integer.parseInt(args[1]);
		final String TRANSLATED_FILE = args[2];
		final String BOUNDARIES_FILE = args[3];
		final String READ_LENGTHS_FILE = args[4];
		final String ALPHABET_FILE = args[5];
		final String AVOIDED_INTERVALS_FILE = args[6];  // NULL to discard it
        final String KMERS_FILE_INPUT = args[7];  // NULL to discard it
		final String KMERS_FILE_OUTPUT = args[8];  // Output file
		
		boolean INTERVALS_FILE_EXISTS = !AVOIDED_INTERVALS_FILE.equalsIgnoreCase("null");
		
		int i;
		int row, nKmers, lastAvoidedInterval, readLength;
		String str1, str2, str3, str4;
		BufferedReader br1, br2, br3, br4;
		BufferedWriter bw;
		RepeatAlphabet.Character tmpChar = new RepeatAlphabet.Character();
		RepeatAlphabet.Kmer tmpKmer = new RepeatAlphabet.Kmer();
		HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer> kmers, tmpMap;
		int[] avoidedIntervals;
		int[] tmpArray2 = new int[K];
		int[] tmpArray3 = new int[(K)<<1];
		String[] tokens;
		RepeatAlphabet.Kmer[] keys;
		
		// Collecting k-mers
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
        if (!KMERS_FILE_INPUT.equalsIgnoreCase("null")) kmers=deserializeKmers(KMERS_FILE_INPUT,K);
        else kmers = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
		RepeatAlphabet.kmerPool_init(K);
		tmpMap = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
		br1 = new BufferedReader(new FileReader(TRANSLATED_FILE));
		if (INTERVALS_FILE_EXISTS) {
			br2 = new BufferedReader(new FileReader(AVOIDED_INTERVALS_FILE));
			avoidedIntervals = new int[10];  // Arbitrary, multiple of 2.
		}
		else { br2=null; avoidedIntervals=null; }
		br3 = new BufferedReader(new FileReader(BOUNDARIES_FILE));
		br4 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		str1=br1.readLine(); str2=INTERVALS_FILE_EXISTS?br2.readLine():null; 
		str3=br3.readLine(); str4=br4.readLine(); row=0;
		while (str1!=null) {
			if (row%100000==0) System.err.println("Processed "+row+" reads, "+kmers.size()+" distinct "+K+"-mers.");
			if (INTERVALS_FILE_EXISTS) {
				if (str2.length()==0) lastAvoidedInterval=-1;
				else {
					tokens=str2.split(RepeatAlphabet.SEPARATOR_MINOR+""); lastAvoidedInterval=tokens.length-1;
					if (lastAvoidedInterval>=avoidedIntervals.length) avoidedIntervals = new int[lastAvoidedInterval+1];
					for (i=0; i<=lastAvoidedInterval; i++) avoidedIntervals[i]=Integer.parseInt(tokens[i]);
				}
			}
			else lastAvoidedInterval=-1;
			RepeatAlphabet.loadBoundaries(str3);
			readLength=Integer.parseInt(str4);
			RepeatAlphabet.getKmers(MODE,str1,K,kmers,null,avoidedIntervals,lastAvoidedInterval,readLength,-1/*argument not used*/,-1/*argument not used*/,-1/*argument not used*/,-1/*argument not used*/,-1/*argument not used*/,-1/*argument not used*/,RepeatAlphabet.boundaries,-1/*argument not used*/,-1/*argument not used*/,-1.0/*argument not used*/,tmpKmer,tmpArray2,tmpArray3,tmpMap,tmpChar);
			str1=br1.readLine(); str2=INTERVALS_FILE_EXISTS?br2.readLine():null; 
			str3=br3.readLine(); str4=br4.readLine(); row++;
		}
		br1.close(); br3.close(); br4.close(); nKmers=kmers.size();
		if (INTERVALS_FILE_EXISTS) br2.close();
		System.err.println(nKmers+" total distinct "+K+"-mers");
		
		// Serializing sorted k-mers
		keys = new RepeatAlphabet.Kmer[nKmers];
		kmers.keySet().toArray(keys);
		if (nKmers>1) Arrays.sort(keys,0,nKmers);
		bw = new BufferedWriter(new FileWriter(KMERS_FILE_OUTPUT));
		for (i=0; i<nKmers; i++) {
            bw.write(keys[i].toString());
            if (MODE==1) bw.write((RepeatAlphabet.SEPARATOR_MINOR+"")+keys[i].count+(RepeatAlphabet.SEPARATOR_MINOR+"")+keys[i].countPartial+(RepeatAlphabet.SEPARATOR_MINOR+"")+keys[i].sameReadCount);
            bw.newLine();
        }
		bw.close();
	}
    
    
    private static final HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer> deserializeKmers(String path, int k) throws IOException {
        String str;
        BufferedReader br;
        RepeatAlphabet.Kmer kmer;
        HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer> out;
        
        out = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();
        while (str!=null) {
            kmer = new RepeatAlphabet.Kmer(str,k);
			out.put(kmer,kmer);
            str=br.readLine();
        }
        br.close();
        return out;
    }
    
}