package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Reads;

/**
 * Updates the translation of every read to use the alphabet built by
 * $CleanTranslatedReads2$.
 */
public class CleanTranslatedReads3 {
	
	public static void main(String[] args) throws IOException {
		final int N_READS = Integer.parseInt(args[0]);
		final String READ_IDS_FILE = args[1];
		final String READ_LENGTHS_FILE = args[2];
		final String ALPHABET_FILE_OLD = args[3];
		final String ALPHABET_COUNTS_FILE_OLD = args[4];
		final String TRANSLATED_FILE_CHARACTERS_OLD = args[5];
		final String TRANSLATED_FILE_BOUNDARIES_OLD = args[6];
		final int MIN_FREQUENCY = Integer.parseInt(args[7]);
		final String ALPHABET_FILE_NEW = args[8];
		final String OLD2NEW_FILE = args[9];
		final String TRANSLATED_FILE_CHARACTERS_NEW = args[10];
		final String TRANSLATED_FILE_BOUNDARIES_NEW = args[11];		
		
		int i;
		int lastUnique_old, lastPeriodic_old, lastAlphabet_old;
		int lastUnique_new, lastPeriodic_new, lastAlphabet_new;
		String str1, str2;
		RepeatAlphabet.Character tmpChar = new RepeatAlphabet.Character();
		BufferedReader br1, br2;
		BufferedWriter bw1, bw2;
		int[] old2new;
		RepeatAlphabet.Character[] oldAlphabet, newAlphabet;
		
		// Loading old alphabet
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_OLD,2);
		RepeatAlphabet.loadAlphabetCount(ALPHABET_COUNTS_FILE_OLD,RepeatAlphabet.lastAlphabet+1);
		oldAlphabet=RepeatAlphabet.alphabet;
		lastUnique_old=RepeatAlphabet.lastUnique;
		lastPeriodic_old=RepeatAlphabet.lastPeriodic;
		lastAlphabet_old=RepeatAlphabet.lastAlphabet;
		br1 = new BufferedReader(new FileReader(OLD2NEW_FILE));
		old2new = new int[lastAlphabet_old-lastUnique_old];
		for (i=0; i<old2new.length; i++) old2new[i]=Integer.parseInt(br1.readLine());
		br1.close();
		
		// Updating the translation
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_NEW,2);
		newAlphabet=RepeatAlphabet.alphabet;
		lastUnique_new=RepeatAlphabet.lastUnique;
		lastPeriodic_new=RepeatAlphabet.lastPeriodic;
		lastAlphabet_new=RepeatAlphabet.lastAlphabet;
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,Reads.nReads);
		br1 = new BufferedReader(new FileReader(TRANSLATED_FILE_CHARACTERS_OLD));
		br2 = new BufferedReader(new FileReader(TRANSLATED_FILE_BOUNDARIES_OLD));
		bw1 = new BufferedWriter(new FileWriter(TRANSLATED_FILE_CHARACTERS_NEW));
		bw2 = new BufferedWriter(new FileWriter(TRANSLATED_FILE_BOUNDARIES_NEW));
		i=0; str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null) {		
			RepeatAlphabet.cleanTranslatedRead_updateTranslation(str1,str2,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,old2new,Reads.readLengths[i],MIN_FREQUENCY,IO.quantum,bw1,bw2,tmpChar);
			bw1.newLine(); bw2.newLine();
			i++; str1=br1.readLine(); str2=br2.readLine();
		}
		br1.close(); br2.close(); bw1.close(); bw2.close();
	}

}