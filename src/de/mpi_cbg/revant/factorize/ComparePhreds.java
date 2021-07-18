package de.mpi_cbg.revant.factorize;

import java.io.*;
import de.mpi_cbg.revant.util.IO;

/**
 * Just a command-line version of $Reads.comparePhred2track()$ and  
 * $Reads.compareTrack2track()$ that compares a dust track and a tandem track. See those
 * procedures for details.
 */
public class ComparePhreds {
	
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE_1 = args[0];
		final String INPUT_FILE_2 = args[1];
		final String DUST_FILE = args[2];
		final String TANDEM_FILE = args[3];
		
		IO.initialize();
		Reads.QUALITY_SPACING = 100;
		final int minAlignmentLength = 1500;
		final int identityThreshold = 100;
		
		Reads.comparePhreds(INPUT_FILE_1,INPUT_FILE_2,51,true);
		Reads.comparePhreds(INPUT_FILE_1,INPUT_FILE_2,45,false);
		Reads.comparePhreds(INPUT_FILE_1,INPUT_FILE_2,40,false);
		Reads.comparePhreds(INPUT_FILE_1,INPUT_FILE_2,35,false);
		Reads.comparePhreds(INPUT_FILE_1,INPUT_FILE_2,30,false);
		Reads.comparePhreds(INPUT_FILE_1,INPUT_FILE_2,25,false);
		
		System.out.println();
		System.err.println("COMPARING QUALITY FILE 1 TO DUST TRACK:  (merge threshold="+minAlignmentLength+")");
		Reads.comparePhred2track(INPUT_FILE_1,DUST_FILE,51,minAlignmentLength);
		System.out.println();
		System.err.println("COMPARING QUALITY FILE 1 TO TANDEM TRACK:  (merge threshold="+identityThreshold+")");
		Reads.comparePhred2track(INPUT_FILE_1,TANDEM_FILE,51,identityThreshold);
		System.out.println();
		System.err.println("COMPARING DUST TRACK (1) TO TANDEM TRACK (2):");
		Reads.compareTrack2track(DUST_FILE,TANDEM_FILE,minAlignmentLength,identityThreshold);
		
	}
	
}
