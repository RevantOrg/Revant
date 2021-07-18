package de.mpi_cbg.revant.factorize;

import java.io.*;

/**
 * This is just a command-line version of $Reads.qualityHistogram()$: see that procedure
 * for details.
 */
public class QualityHistogram  {

	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		final int MAX_QUALITY = Integer.parseInt(args[1]);
		
		long[] histogram = new long[MAX_QUALITY+1];
		Reads.qualityHistogram(INPUT_FILE,histogram);
		for (int i=0; i<=MAX_QUALITY; i++) System.out.println(i+","+histogram[i]);
	}

}


