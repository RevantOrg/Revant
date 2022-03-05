package de.mpi_cbg.revant.apps;

import java.util.Random;
import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Builds a random quality track with long low-quality intervals of fixed length.
 */
public class BuildRandomQualities {
	/**
	 * @param args 3: 0=uses the convention of DASqv and DBdump; 1=uses the convention of
     * DAmar.
	 */
	public static void main(String[] args) throws IOException {
		final String READ_LENGTHS_FILE = args[0];
		final int LOW_QUALITY_LENGTH = Integer.parseInt(args[1]);
		final double LOW_QUALITY_PROBABILITY = Double.parseDouble(args[2]);
		final byte MODE = Byte.parseByte(args[3]);
		final String QUALITY_THRESHOLDS_FILE = args[4];
		final String OUTPUT_FILE = args[5];
		
		final Random random = new Random();
		Reads.loadQualities_thresholds(QUALITY_THRESHOLDS_FILE);
		Reads.breakReads_buildRandomQualities(READ_LENGTHS_FILE,OUTPUT_FILE,LOW_QUALITY_LENGTH,LOW_QUALITY_PROBABILITY,MODE,random);
	}

}