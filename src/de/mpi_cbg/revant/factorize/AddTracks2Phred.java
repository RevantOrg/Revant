package de.mpi_cbg.revant.factorize;

import java.io.*;
import de.mpi_cbg.revant.util.IO;

/**
 * Just a command-line version of $Reads.addTracksToPhred()$: see that procedure for 
 * details.
 */
public class AddTracks2Phred {
	
	public static void main(String[] args) throws IOException {
		final byte MIN_BAD_QUALITY = Byte.parseByte(args[0]);
		final String PHRED_FILE = args[1];
		final String QUALITY_THRESHOLDS_FILE = args[2];
		final String DUST_FILE = args[3];
		final String TANDEM_FILE = args[4];
		final String OUTPUT_FILE = args[5];
		
		IO.initialize();
		Reads.loadQualities_thresholds(QUALITY_THRESHOLDS_FILE);
		Reads.addTracksToPhred(PHRED_FILE,MIN_BAD_QUALITY,QUALITY_THRESHOLDS_FILE,DUST_FILE,TANDEM_FILE,OUTPUT_FILE);
	}
	
}
