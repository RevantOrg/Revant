package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Let a shortest unique intervals file contain, for every translated read, the list of
 * all intervals such that: (1) every interval is the occurrence of a unique h-mer for all
 * h < k; (2) intervals are not contained in one another. The program adds to every such
 * list all the occurrences of unique k-mers that do not contain an interval already in
 * the list.
 *
 * Remark: when we mark a unique k-mer in a read, we don't propagate this information to
 * the reads that align to it, so the aligned region of another read might not be marked
 * as unique, and it might even have completely different characters and block boundaries.
 */
public class GetShortestUniqueIntervals {
	
	public static void main(String[] args) throws IOException {
		final int K = Integer.parseInt(args[0]);
		final String TRANSLATED_FILE = args[1];
		final String BOUNDARIES_FILE = args[2];
		final String READ_LENGTHS_FILE = args[3];
		final String ALPHABET_FILE = args[4];
		final String UNIQUE_KMERS_FILE = args[5];
		final int HAPLOTYPE_COVERAGE = Integer.parseInt(args[6]);
		final int IDENTITY_THRESHOLD = Integer.parseInt(args[7]);
		final int DISTANCE_THRESHOLD = Integer.parseInt(args[8]);
		final double CHARACTER_FRACTION = Double.parseDouble(args[9]);
		final String OLD_INTERVALS_FILE = args[10];  // NULL to discard it
		final String NEW_INTERVALS_FILE = args[11];  // Output
		
		boolean OLD_INTERVALS_FILE_EXISTS = !OLD_INTERVALS_FILE.equalsIgnoreCase("null");
		
		int i, p;
		int row, nBlocks, nPairs, lastUniqueInterval, readLength;
		String str1, str2, str3, str4;
		BufferedReader br1, br2, br3, br4;
		BufferedWriter bw;
		RepeatAlphabet.Character tmpChar = new RepeatAlphabet.Character();
		RepeatAlphabet.Kmer tmpKmer = new RepeatAlphabet.Kmer();
		RepeatAlphabet.Kmer newKmer;
		HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer> kmers;
		int[] uniqueIntervals;
		int[] tmpArray2 = new int[K];
		int[] tmpArray3 = new int[(K)<<1];
		String[] tokens;
		Pair[] pairs;
		
		// Loading unique k-mers
		kmers = new HashMap<RepeatAlphabet.Kmer,RepeatAlphabet.Kmer>();
		br1 = new BufferedReader(new FileReader(UNIQUE_KMERS_FILE));
		str1=br1.readLine();
		while (str1!=null) {
			newKmer = new RepeatAlphabet.Kmer(str1,K);
			kmers.put(newKmer,newKmer);
			str1=br1.readLine();
		}
		br1.close();
		
		// Building the new intervals file
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		br1 = new BufferedReader(new FileReader(TRANSLATED_FILE));
		uniqueIntervals = new int[30];  // Arbitrary, multiple of 3.
		pairs = new Pair[10];  // Arbitrary
		for (i=0; i<pairs.length; i++) pairs[i] = new Pair();
		if (OLD_INTERVALS_FILE_EXISTS) br2 = new BufferedReader(new FileReader(OLD_INTERVALS_FILE));
		else br2=null;
		br3 = new BufferedReader(new FileReader(BOUNDARIES_FILE));
		br4 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		bw = new BufferedWriter(new FileWriter(NEW_INTERVALS_FILE));
		str1=br1.readLine(); str2=OLD_INTERVALS_FILE_EXISTS?br2.readLine():null; 
		str3=br3.readLine(); str4=br4.readLine(); row=0;
		while (str1!=null) {
			if (row%100000==0) System.err.println("Processed "+row+" reads");
			nBlocks=0; p=str1.indexOf(RepeatAlphabet.SEPARATOR_MAJOR+"");
			while (p>=0) {
				nBlocks++;
				p=str1.indexOf(RepeatAlphabet.SEPARATOR_MAJOR+"",p+1);
			}
			nBlocks++;
			if (uniqueIntervals.length<nBlocks*3) uniqueIntervals = new int[nBlocks*3];
			if (OLD_INTERVALS_FILE_EXISTS) {
				if (str2.length()==0) lastUniqueInterval=-1;
				else {
					tokens=str2.split(RepeatAlphabet.SEPARATOR_MINOR+""); lastUniqueInterval=tokens.length-1;
					for (i=0; i<=lastUniqueInterval; i++) uniqueIntervals[i]=Integer.parseInt(tokens[i]);
				}
			}
			else lastUniqueInterval=-1;
			RepeatAlphabet.loadBoundaries(str3);
			readLength=Integer.parseInt(str4);
			lastUniqueInterval=RepeatAlphabet.getKmers(str1,K,null,kmers,uniqueIntervals,lastUniqueInterval,HAPLOTYPE_COVERAGE,readLength,RepeatAlphabet.boundaries,IDENTITY_THRESHOLD,DISTANCE_THRESHOLD,CHARACTER_FRACTION,tmpKmer,tmpArray2,tmpArray3,null,tmpChar);
			if (lastUniqueInterval>0) {
				nPairs=(lastUniqueInterval+1)/3;
				if (pairs.length<nPairs) {
					Pair[] newArray = new Pair[nPairs];
					System.arraycopy(pairs,0,newArray,0,pairs.length);
					for (i=pairs.length; i<nPairs; i++) newArray[i] = new Pair();
					pairs=newArray;
				}
				for (i=0; i<nPairs; i++) pairs[i].set(uniqueIntervals[3*i],uniqueIntervals[3*i+1],uniqueIntervals[3*i+2]);
				if (nPairs>1) Arrays.sort(pairs,0,nPairs);
				bw.write(pairs[0].position+(RepeatAlphabet.SEPARATOR_MINOR+"")+pairs[0].length+(RepeatAlphabet.SEPARATOR_MINOR+"")+pairs[0].nHaplotypes);
				for (i=1; i<nPairs; i++) bw.write((RepeatAlphabet.SEPARATOR_MINOR+"")+pairs[i].position+(RepeatAlphabet.SEPARATOR_MINOR+"")+pairs[i].length+(RepeatAlphabet.SEPARATOR_MINOR+"")+pairs[i].nHaplotypes);
			}
			bw.newLine();
			str1=br1.readLine(); str2=OLD_INTERVALS_FILE_EXISTS?br2.readLine():null; 
			str3=br3.readLine(); str4=br4.readLine(); row++;
		}
		if (OLD_INTERVALS_FILE_EXISTS) br2.close();
		br1.close(); br3.close(); br4.close(); bw.close();
	}


	private static class Pair implements Comparable {
		public int position, length, nHaplotypes;
		
		public Pair() { }
		
		public void set(int p, int l, int n) { this.position=p; this.length=l; this.nHaplotypes=n; }
		
		public int compareTo(Object other) {
			Pair otherPair = (Pair)other;
			if (position<otherPair.position) return -1;
			else if (position>otherPair.position) return 1;
			return 0;
		}
		
		public boolean equals(Object other) {
			Pair otherPair = (Pair)other;
			return position==otherPair.position && length==otherPair.length && nHaplotypes==otherPair.nHaplotypes;
		}
	}

}