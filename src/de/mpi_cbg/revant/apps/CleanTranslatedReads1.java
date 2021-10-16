package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Reads;

/**
 * Collects new unique characters induced by removing rare characters.
 * Prints to STDOUT the new value of $maxOpenLength_unique$.
 */
public class CleanTranslatedReads1 {
	
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String ALPHABET_COUNTS_FILE = args[1];
		final int N_READS = Integer.parseInt(args[2]);
		final String READ_IDS_FILE = args[3];
		final String READ_LENGTHS_FILE = args[4];
		final String TRANSLATED_READS_CHARACTERS_FILE = args[5];
		final String TRANSLATED_READS_BOUNDARIES_FILE = args[6];
		final int MIN_FREQUENCY = Integer.parseInt(args[7]);
		final String OUTPUT_FILE = args[8];
		
		int i;
		String str1, str2;
		RepeatAlphabet.Character tmpChar = new RepeatAlphabet.Character();
		BufferedReader br1, br2;
		BufferedWriter bw;
		
		Reads.nReads=N_READS;
		Reads.loadReadIDs(READ_IDS_FILE,Reads.nReads);
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE);
		RepeatAlphabet.loadAlphabetCount(ALPHABET_COUNTS_FILE,RepeatAlphabet.lastAlphabet+1);		
		br1 = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
		br2 = new BufferedReader(new FileReader(TRANSLATED_READS_BOUNDARIES_FILE));
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		i=0; str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null) {
			RepeatAlphabet.cleanTranslatedRead_collectCharacterInstances(str1,str2,Reads.readLengths[i],MIN_FREQUENCY,IO.quantum,bw,tmpChar);
			i++; str1=br1.readLine(); str2=br2.readLine();
		}
		br1.close(); br2.close(); bw.close();
		System.out.println(RepeatAlphabet.maxOpenLength_unique+"");
	}

}