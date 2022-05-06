package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;

/**
 * Uses spacers, and the new alphabet induced by them, to update the translation of every
 * read. This is designed to work on a chunk of reads.
 */
public class FixPeriodicEndpoints3 {

	public static void main(String[] args) throws IOException {
		final int MAX_SPACER_LENGTH = Integer.parseInt(args[0]);
		final String SPACERS_FILE = args[1];  // Of a chunk of reads
		final int N_SPACERS = Integer.parseInt(args[2]);
		final String ALPHABET_FILE_OLD = args[3];  // Of all reads
		final String ALPHABET_FILE_NEW = args[4];  // Of all reads
		final String READ_IDS_FILE = args[5];  // Of a chunk of reads
		final String READ_LENGTHS_FILE = args[6];  // Of a chunk of reads
		final String READ2CHARACTERS_FILE_OLD = args[7];  // Of a chunk of reads
		final String READ2BOUNDARIES_FILE_OLD = args[8];  // Of a chunk of reads
		final String READ2CHARACTERS_FILE_NEW = args[9];  // Of a chunk of reads
		final String READ2BOUNDARIES_FILE_NEW = args[10];  // Of a chunk of reads
		
		int i, j;
		int sum, lastUnique_old, lastPeriodic_old, lastAlphabet_old, lastUnique_new, lastPeriodic_new, lastAlphabet_new;
		String str1, str2, str3, str4;
		BufferedReader br1, br2, br3, br4;
		BufferedWriter bw1, bw2;
		boolean[] used;
		int[] histogram;
		RepeatAlphabet.Character[] oldAlphabet, newAlphabet;
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_OLD,2);
		oldAlphabet=RepeatAlphabet.alphabet;
		lastUnique_old=RepeatAlphabet.lastUnique;
		lastPeriodic_old=RepeatAlphabet.lastPeriodic;
		lastAlphabet_old=RepeatAlphabet.lastAlphabet;
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_NEW,2);
		newAlphabet=RepeatAlphabet.alphabet;
		lastUnique_new=RepeatAlphabet.lastUnique;
		lastPeriodic_new=RepeatAlphabet.lastPeriodic;
		lastAlphabet_new=RepeatAlphabet.lastAlphabet;
		RepeatAlphabet.deserializeSpacers(SPACERS_FILE,N_SPACERS);
		histogram = new int[11];  // Arbitrary
		Math.set(histogram,0,histogram.length-1);
		br1 = new BufferedReader(new FileReader(READ_IDS_FILE));
		br2 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		br3 = new BufferedReader(new FileReader(READ2CHARACTERS_FILE_OLD));
		br4 = new BufferedReader(new FileReader(READ2BOUNDARIES_FILE_OLD));
		bw1 = new BufferedWriter(new FileWriter(READ2CHARACTERS_FILE_NEW));
		bw2 = new BufferedWriter(new FileWriter(READ2BOUNDARIES_FILE_NEW));
		i=0; j=0;
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); str4=br4.readLine();
		while (str1!=null) {
			j=RepeatAlphabet.fixPeriodicEndpoints_updateTranslation(Integer.parseInt(str1),Integer.parseInt(str2),j,MAX_SPACER_LENGTH,str3,str4,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,bw1,bw2,histogram);
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); str4=br4.readLine();
		}
		br1.close(); br2.close(); br3.close(); br4.close(); bw1.close(); bw2.close();
		sum=0;
		for (i=0; i<histogram.length; i++) sum+=histogram[i];
		System.err.println("Fixed "+sum+" spacers. Histogram of fixed spacers lengths:");
		for (i=0; i<histogram.length; i++) System.err.println((i*IO.quantum)+": "+histogram[i]);
	}

}