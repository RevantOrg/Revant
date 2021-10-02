package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Intervals;


public class BuildRepeatAlphabet {
	
	public static void main(String[] args) throws IOException {
		final int N_READS = Integer.parseInt(args[0]);
		final String READ_IDS_FILE = args[1];
		final String READ_LENGTHS_FILE = args[2];
		final int N_REPEATS = Integer.parseInt(args[3]);
		final String REPEAT_LENGTHS_FILE = args[4];
		final String ISPERIODIC_FILE = args[5];
		final String ALIGNMENTS_FILE = args[6];
		final double MAX_ERROR = Double.parseDouble(args[7]);
		final String ALPHABET_FILE = args[8];
		
		RepeatAlphabet.initialize(N_READS,READ_IDS_FILE,READ_LENGTHS_FILE,N_REPEATS,REPEAT_LENGTHS_FILE,ISPERIODIC_FILE);
		RepeatAlphabet.buildAlphabet(ALIGNMENTS_FILE,MAX_ERROR,IO.quantum,IO.quantum<<1);
		RepeatAlphabet.serializeAlphabet(ALPHABET_FILE);
	}


}