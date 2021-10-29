package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;

/**
 * Given a sorted file that contains the k-mers extracted by several threads, the program
 * sums up all counts of the same k-mer, it discards k-mers whose total count is outside a 
 * given range, and it prints a histogram of total counts for all k-mers.
 */
public class CompactKmers {
	
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		final int K = Integer.parseInt(args[1]);
		final int MIN_COUNT = Integer.parseInt(args[2]);
		final int MAX_COUNT = Integer.parseInt(args[3]);
		final String OUTPUT_FILE_KMERS = args[4];
		final int MAX_HISTOGRAM_COUNT = Integer.parseInt(args[5]);
		final String OUTPUT_FILE_HISTOGRAM = args[6];
		
		boolean equal;
		int i;
		long count, previousCount;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		int[] previous, current, tmpArray;
		long[] histogram;
		String[] tokens;
		
		br = new BufferedReader(new FileReader(INPUT_FILE));
		str=br.readLine();
		if (str==null) {
			br.close();
			System.err.println("ERROR: empty file "+INPUT_FILE);
			System.exit(1);
		}
		histogram = new long[MAX_HISTOGRAM_COUNT+1];
		Math.set(histogram,MAX_HISTOGRAM_COUNT,0);
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_KMERS));
		previous = new int[K];
		tokens=str.split(",");
		for (i=0; i<K; i++) previous[i]=Integer.parseInt(tokens[i]);
		previousCount=Integer.parseInt(tokens[K]);
		current = new int[K];
		str=br.readLine();
		while (str!=null) {
			tokens=str.split(",");
			equal=true;
			for (i=0; i<K; i++) {
				current[i]=Integer.parseInt(tokens[i]);
				if (current[i]!=previous[i]) equal=false;
			}
			count=Long.parseLong(tokens[K]);
			if (equal) previousCount+=count;
			else {
				if (previousCount>=MIN_COUNT && previousCount<=MAX_COUNT) {
					for (i=0; i<K; i++) bw.write(previous[i]+",");
					bw.write(previousCount+"\n");
				}
				histogram[previousCount>MAX_HISTOGRAM_COUNT?MAX_HISTOGRAM_COUNT:(int)previousCount]++;
				tmpArray=previous; previous=current; current=tmpArray;
				previousCount=count;
			}
			str=br.readLine();
		}
		br.close();
		if (previousCount>=MIN_COUNT && previousCount<=MAX_COUNT) {
			for (i=0; i<K; i++) bw.write(previous[i]+",");
			bw.write(previousCount+"\n");
		}
		bw.close();
		histogram[previousCount>MAX_HISTOGRAM_COUNT?MAX_HISTOGRAM_COUNT:(int)previousCount]++;
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_HISTOGRAM));
		for (i=0; i<=MAX_HISTOGRAM_COUNT; i++) bw.write(i+","+histogram[i]+"\n");
		bw.close();
	}

}