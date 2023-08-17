package de.mpi_cbg.revant.apps;

import java.io.*;
import org.apache.commons.statistics.distribution.PoissonDistribution;

import de.mpi_cbg.revant.util.Math;

/**
 * Given a sorted file that contains the k-mers extracted by several threads, the program
 * sums up all the counts of the same k-mer, it discards k-mers whose total count fails a 
 * significance test at every possible ploidy of the k-mer, or that occur multiple times
 * inside the same read, and it prints a histogram of total counts for all k-mers.
 */
public class CompactKmers {
    /**
     * @param args 9 TRUE=keep every k-mer that passes a one-sided significance test in
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
		final boolean DISCARD_SAME_READ_KMERS = Integer.parseInt(args[7])==1;
        final String ALPHABET_FILE = args[8];
        final boolean KEEP_ALL_FREQUENT = Integer.parseInt(args[9])==1;
        final int MAX_HISTOGRAM_COUNT = Integer.parseInt(args[10]);
        final String OUTPUT_FILE_KMERS = args[11];
		final String OUTPUT_FILE_HISTOGRAM = args[12].equalsIgnoreCase("null")?null:args[12];
		
        final double SIGNIFICANCE_LEVEL = 0.05;  // Conventional
        
		boolean equal;
		int i;
		int count, previousCount, sameReadCount, previousSameReadCount;
		String str;
		BufferedReader br;
		BufferedWriter bw;
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
		tokens=str.split(",");
		for (i=0; i<K; i++) previous[i]=Integer.parseInt(tokens[i]);
		previousCount=Integer.parseInt(tokens[K]);
		previousSameReadCount=Integer.parseInt(tokens[K+1]);
		current = new int[K];
		str=br.readLine();
		while (str!=null) {
			tokens=str.split(",");
			equal=true;
			for (i=0; i<K; i++) {
				current[i]=Integer.parseInt(tokens[i]);
				if (current[i]!=previous[i]) equal=false;
			}
			count=Integer.parseInt(tokens[K]);
			sameReadCount=Integer.parseInt(tokens[K+1]);
			if (equal) {
				previousCount+=count;
				previousSameReadCount=Math.max(previousSameReadCount,sameReadCount);
			}
			else {
                if ( (DISCARD_SAME_READ_KMERS?previousSameReadCount==1:true) && 
                     (KEEP_ALL_FREQUENT?isFrequent(previousCount,previous,K,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,SIGNIFICANCE_LEVEL):(isUnique(previousCount,previous,K,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,SIGNIFICANCE_LEVEL)!=-1))
                   ) {
					for (i=0; i<K; i++) bw.write(previous[i]+",");
					bw.write(previousCount+"\n");
				}
				if (OUTPUT_FILE_HISTOGRAM!=null) histogram[previousCount>MAX_HISTOGRAM_COUNT?MAX_HISTOGRAM_COUNT:previousCount]++;
				tmpArray=previous; previous=current; current=tmpArray;
				previousCount=count; previousSameReadCount=sameReadCount;
			}
			str=br.readLine();
		}
		br.close();
		if ( (DISCARD_SAME_READ_KMERS?previousSameReadCount==1:true) && 
             (KEEP_ALL_FREQUENT?isFrequent(previousCount,previous,K,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,SIGNIFICANCE_LEVEL):(isUnique(previousCount,previous,K,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,SIGNIFICANCE_LEVEL)!=-1))
           ) {
			for (i=0; i<K; i++) bw.write(previous[i]+",");
			bw.write(previousCount+"\n");
		}
		bw.close();
		if (OUTPUT_FILE_HISTOGRAM!=null) {
			histogram[previousCount>MAX_HISTOGRAM_COUNT?MAX_HISTOGRAM_COUNT:previousCount]++;
			bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_HISTOGRAM));
			for (i=0; i<=MAX_HISTOGRAM_COUNT; i++) bw.write(i+","+histogram[i]+"\n");
			bw.close();
		}
	}
    
    
    /**
     * @param nReads in the entire dataset;
     * @param spanningBps basepairs before and after $kmer$ for it to be considered
     * observed in a read;
     * @param genomeLength of one haplotype;
     * @return -1 if the p-value of $count$ fails the two-sided significance test for
     * every ploidy of $kmer$ (assuming a Poisson distribution for every possible ploidy);
     * this happens if $kmer$ is noise or a repeat; otherwise, the smallest ploidy of
     * $kmer$ for which the significance test does not fail.
     */
    private static final int isUnique(int count, int[] kmer, int k, int nReads, int avgReadLength, int spanningBps, long genomeLength, int nHaplotypes, double significanceLevel) {
        int i, length;
        double p;
        PoissonDistribution distribution;
        
        length=0;
        for (i=0; i<k; i++) RepeatAlphabet.alphabet[kmer[i]].getLength();
        final double base = ((double)(avgReadLength-((spanningBps)<<1)-length))/genomeLength;
        final int quantum = nReads/nHaplotypes;
        for (i=0; i<nHaplotypes; i++) {
            distribution=PoissonDistribution.of(base*quantum*(i+1));
            p=distribution.cumulativeProbability(count);
            p=2.0*Math.min(p,1.0-p);
            if (p>significanceLevel) return i+1;
        }
        return -1;
    }
    
    
    /**
     * Just a one-sided test on the model with one haplotype.
     */
    private static final boolean isFrequent(int count, int[] kmer, int k, int nReads, int avgReadLength, int spanningBps, long genomeLength, int nHaplotypes, double significanceLevel) {
        int i, length;
        PoissonDistribution distribution;
        
        length=0;
        for (i=0; i<k; i++) RepeatAlphabet.alphabet[kmer[i]].getLength();
        final double base = ((double)(avgReadLength-((spanningBps)<<1)-length))/genomeLength;
        final int quantum = nReads/nHaplotypes;
        distribution=PoissonDistribution.of(base*quantum);
        return distribution.cumulativeProbability(count)>significanceLevel;
    }

}