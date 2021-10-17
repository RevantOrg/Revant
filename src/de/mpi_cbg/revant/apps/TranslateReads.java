package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Reads;


public class TranslateReads {
	
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
		final int LAST_TRANSLATED_READ = Integer.parseInt(args[9]);
		final String TRANSLATED_FILE_CHARACTERS = args[10];
		final String TRANSLATED_FILE_BOUNDARIES = args[11];
		
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,Reads.nReads);
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		RepeatAlphabet.loadIsPeriodic(ISPERIODIC_FILE,N_REPEATS);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.translateReads(ALIGNMENTS_FILE,LAST_TRANSLATED_READ,MAX_ERROR,IO.quantum,TRANSLATED_FILE_CHARACTERS,TRANSLATED_FILE_BOUNDARIES);
	}

}