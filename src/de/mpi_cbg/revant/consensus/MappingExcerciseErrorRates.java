package de.mpi_cbg.revant.consensus;

import java.io.*;
import java.util.*;
import de.mpi_cbg.revant.util.Math;

/**
 * Used in the experiment in which we mapped my lungfish modules to 1x of the reads.
 * Just computes the distribution of error rates from an input file that contains just 
 * those, filtered out form a LAshow file.
 */
public class MappingExcerciseErrorRates {
		
	/**
	 * java MappingExcerciseErrorRates lashowFile nBins
	 */
	public static void main(String[] args) throws Exception {
		final String INPUT_FILE = args[0];
		final int N_BINS = Integer.parseInt(args[1]);
		final double QUANTUM = 1.0/N_BINS;
		int i;
		double value;
		String str;
		BufferedReader br;
		int[] histogram = new int[N_BINS+1];
		
		// Building histogram
		Math.set(histogram,histogram.length-1,0);
		br = new BufferedReader(new FileReader(INPUT_FILE));
		str=br.readLine();
		while (str!=null) {
			value=Double.parseDouble(str)/100;
			histogram[(int)(value/QUANTUM)]++;
			str=br.readLine();
		}
		br.close(); br=null;
		
		// Printing
		for (i=0; i<histogram.length; i++) System.out.println((i*QUANTUM)+","+histogram[i]);
		//for (i=histogram.length-2; i>=0; i--) histogram[i]+=histogram[i+1];
		for (i=1; i<histogram.length; i++) histogram[i]+=histogram[i-1];
		for (i=0; i<histogram.length; i++) System.err.println((i*QUANTUM)+","+histogram[i]);
	}

}