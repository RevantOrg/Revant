package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Intervals;


public class TranslateReads {
	
	public static void main(String[] args) throws IOException {
		final int N_READS = Integer.parseInt(args[0]);
		final String READ_IDS_FILE = args[1];
		final String READ_LENGTHS_FILE = args[2];
		final String ALIGNMENTS_FILE = args[3];
		final double MAX_ERROR = Double.parseDouble(args[4]);
		final String ALPHABET_FILE = args[5];
		final int ALPHABET_SIZE = Integer.parseInt(args[6]);
		final String TRANSLATED_FILE = args[7];
		
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,Reads.nReads);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,ALPHABET_SIZE);
		RepeatAlphabet.translateReads(ALIGNMENTS_FILE,MAX_ERROR,IO.quantum,IO.quantum<<1,TRANSLATED_FILE);
	}


}