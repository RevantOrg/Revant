package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Intervals;


public class GetCharacterInstances {
	
	public static void main(String[] args) throws IOException {
		final int N_READS = Integer.parseInt(args[0]);
		final String READ_IDS_FILE = args[1];
		final String READ_LENGTHS_FILE = args[2];
		final int N_REPEATS = Integer.parseInt(args[3]);
		final String REPEAT_LENGTHS_FILE = args[4];
		final String ISPERIODIC_FILE = args[5];
		final String ALIGNMENTS_FILE = args[6];
		final double MAX_ERROR = Double.parseDouble(args[7]);
		final String INSTANCES_FILE = args[8];
		
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,Reads.nReads);
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		RepeatAlphabet.loadIsPeriodic(ISPERIODIC_FILE,N_REPEATS);
		RepeatAlphabet.collectCharacterInstances(ALIGNMENTS_FILE,MAX_ERROR,IO.quantum,IO.quantum<<1,INSTANCES_FILE);
		System.err.println("Histogram of recoded sequence lengths:");
		RepeatAlphabet.printSequenceLengths();
	}


}