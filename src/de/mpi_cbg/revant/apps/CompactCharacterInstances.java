package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;


public class CompactCharacterInstances {
	
	public static void main(String[] args) throws IOException {
		final String INSTANCES_FILE = args[0];
		final String OUTPUT_FILE = args[1];
		
		RepeatAlphabet.deserializeAlphabet(INSTANCES_FILE);
		RepeatAlphabet.compactInstances(IO.quantum,IO.quantum<<1);
		RepeatAlphabet.serializeAlphabet(OUTPUT_FILE);
	}


}