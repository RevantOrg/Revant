package de.mpi_cbg.revant.apps;

import java.util.Random;
import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Builds a random quality track with long low-quality intervals of fixed length.
 */
public class BuildRandomQualities {
	/**
	 * @param lowQualityMode every long low-quality region: 1=replaces a substring with
	 * a random substring of approx. the same length; 0=is an insertion;
	 * @param args 3: 0=uses the convention of DASqv and DBdump; 1=uses the convention of
     * DAmar.
	 */
	public static void main(String[] args) throws IOException {
		final boolean LOW_QUALITY_MODE = Integer.parseInt(args[0])==1;
		final String READ_LENGTHS_FILE = args[1];
		final int LOW_QUALITY_LENGTH = Integer.parseInt(args[2]);
		final double LOW_QUALITY_PROBABILITY = Double.parseDouble(args[3]);
		final byte QUALITY_MODE = Byte.parseByte(args[4]);
		final String QUALITY_THRESHOLDS_FILE = args[5];
		final String OUTPUT_FILE = args[6];
		
		final Random random = new Random();
		Reads.loadQualities_thresholds(QUALITY_THRESHOLDS_FILE);
		Reads.breakReads_buildRandomQualities(LOW_QUALITY_MODE,READ_LENGTHS_FILE,OUTPUT_FILE,LOW_QUALITY_LENGTH,LOW_QUALITY_PROBABILITY,QUALITY_MODE,random);
	}

}