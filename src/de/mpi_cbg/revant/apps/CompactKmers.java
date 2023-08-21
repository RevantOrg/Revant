package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;

/**
 * Given a sorted file that contains the k-mers extracted by several threads, the program
 * sums up all the counts of the same k-mer, it discards k-mers whose total count fails a 
 * significance test at every possible ploidy of the k-mer, or that occur multiple times
 * inside the same read, and it prints a histogram of total counts for all k-mers.
 */
public class CompactKmers {
    /**
     * @param args 10 TRUE=keep every k-mer that passes a one-sided significance test in
     * the model where it occurs on just one haplotype.
     */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		final int K = Integer.parseInt(args[1]);
        final long GENOME_LENGTH = Long.parseLong(args[2]);  // Of one haplotype
        final int N_HAPLOTYPES = Integer.parseInt(args[3]);
		final int N_READS = Integer.parseInt(args[4]);
        final int AVG_READ_LENGTH = Integer.parseInt(args[5]);
        final int SPANNING_BPS = Integer.parseInt(args[6]);
        final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[7]);  // Read-repeat
		final boolean DISCARD_SAME_READ_KMERS = Integer.parseInt(args[8])==1;
        final String ALPHABET_FILE = args[9];
        final boolean KEEP_ALL_FREQUENT = Integer.parseInt(args[10])==1;
        final int MAX_HISTOGRAM_COUNT = Integer.parseInt(args[11]);
        final String OUTPUT_FILE_KMERS = args[12];
		final String OUTPUT_FILE_HISTOGRAM = args[13].equalsIgnoreCase("null")?null:args[13];
		
        final double SIGNIFICANCE_LEVEL = 0.05;  // Conventional
        final String SEPARATOR = ",";
        
		boolean equal;
		int i;
		int count, countPartial, previousCount, previousCountPartial, sameReadCount, previousSameReadCount;
		String str;
		BufferedReader br;
		BufferedWriter bw;
        RepeatAlphabet.Kmer kmer = new RepeatAlphabet.Kmer();
		int[] previous, current, tmpArray;
		long[] histogram;
		String[] tokens;
        
        RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		br = new BufferedReader(new FileReader(INPUT_FILE));
		str=br.readLine();
		if (str==null) {
			br.close();
			System.err.println("ERROR: empty file "+INPUT_FILE);
			System.exit(1);
		}
		if (OUTPUT_FILE_HISTOGRAM!=null) {
			histogram = new long[MAX_HISTOGRAM_COUNT+1];
			Math.set(histogram,MAX_HISTOGRAM_COUNT,0);
		}
		else histogram=null;
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_KMERS));
		previous = new int[K];
		tokens=str.split(SEPARATOR);
		for (i=0; i<K; i++) previous[i]=Integer.parseInt(tokens[i]);
		previousCount=Integer.parseInt(tokens[K]);
        previousCountPartial=Integer.parseInt(tokens[K+1]);
		previousSameReadCount=Integer.parseInt(tokens[K+2]);
		current = new int[K];
		str=br.readLine();
		while (str!=null) {
			tokens=str.split(SEPARATOR);
			equal=true;
			for (i=0; i<K; i++) {
				current[i]=Integer.parseInt(tokens[i]);
				if (current[i]!=previous[i]) equal=false;
			}
			count=Integer.parseInt(tokens[K]);
            countPartial=Integer.parseInt(tokens[K+1]);
			sameReadCount=Integer.parseInt(tokens[K+2]);
			if (equal) {
				previousCount+=count;
                previousCountPartial+=countPartial;
				previousSameReadCount=Math.max(previousSameReadCount,sameReadCount);
			}
			else {
                kmer.set(previous,K,previousCount,previousCountPartial,previousSameReadCount);
                if ( (DISCARD_SAME_READ_KMERS?previousSameReadCount==1:true) && 
                     (KEEP_ALL_FREQUENT?kmer.isFrequent(K,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,SIGNIFICANCE_LEVEL):(kmer.isUnique(K,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,MIN_ALIGNMENT_LENGTH,SIGNIFICANCE_LEVEL)!=-1))
                   ) {
					for (i=0; i<K; i++) bw.write(previous[i]+SEPARATOR);
					bw.write(previousCount+SEPARATOR+previousCountPartial+SEPARATOR+previousSameReadCount+"\n");
				}
				if (OUTPUT_FILE_HISTOGRAM!=null) histogram[previousCount+previousCountPartial>MAX_HISTOGRAM_COUNT?MAX_HISTOGRAM_COUNT:previousCount+previousCountPartial]++;
				tmpArray=previous; previous=current; current=tmpArray;
				previousCount=count; previousCountPartial=countPartial; previousSameReadCount=sameReadCount;
			}
			str=br.readLine();
		}
		br.close();
        kmer.set(previous,K,previousCount,previousCountPartial,previousSameReadCount);
		if ( (DISCARD_SAME_READ_KMERS?previousSameReadCount==1:true) && 
             (KEEP_ALL_FREQUENT?kmer.isFrequent(K,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,SIGNIFICANCE_LEVEL):(kmer.isUnique(K,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,MIN_ALIGNMENT_LENGTH,SIGNIFICANCE_LEVEL)!=-1))
           ) {
			for (i=0; i<K; i++) bw.write(previous[i]+SEPARATOR);
			bw.write(previousCount+SEPARATOR+previousCountPartial+SEPARATOR+previousSameReadCount+"\n");
		}
		bw.close();
		if (OUTPUT_FILE_HISTOGRAM!=null) {
			histogram[previousCount+previousCountPartial>MAX_HISTOGRAM_COUNT?MAX_HISTOGRAM_COUNT:previousCount+previousCountPartial]++;
			bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_HISTOGRAM));
			for (i=0; i<=MAX_HISTOGRAM_COUNT; i++) bw.write(i+","+histogram[i]+"\n");
			bw.close();
		}
	}

}