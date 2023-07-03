package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;

/**
 * Uses tandem spacers, and the new alphabet induced by them, to update the translation of
 * every read. This is designed to work on a chunk of reads.
 */
public class FixTandemSpacers3 {

	public static void main(String[] args) throws IOException {
		final String SPACERS_FILE = args[0];  // Of a chunk of reads
		final int N_SPACERS = Integer.parseInt(args[1]);
		final String ALPHABET_FILE_OLD = args[2];  // Of all reads
		final String ALPHABET_FILE_NEW = args[3];  // Of all reads
		final String READ_IDS_FILE = args[4];  // Of a chunk of reads
		final String READ_LENGTHS_FILE = args[5];  // Of a chunk of reads
		final String READ2CHARACTERS_FILE_OLD = args[6];  // Of a chunk of reads
		final String READ2BOUNDARIES_FILE_OLD = args[7];  // Of a chunk of reads
		final String READ2CHARACTERS_FILE_NEW = args[8];  // Of a chunk of reads
		final String TANDEMS_FILE = args[9];  // Of a chunk of reads
		final String REPEAT_LENGTHS_FILE = args[10];
		final int N_REPEATS = Integer.parseInt(args[11]);
		
		final int DISTANCE_THRESHOLD = IO.quantum;
		
		int i, j;
		int sum, lastUnique_old, lastPeriodic_old, lastAlphabet_old, lastUnique_new, lastPeriodic_new, lastAlphabet_new, nBlocks;
		String str1, str2, str3, str4, str5;
		BufferedReader br1, br2, br3, br4, br5;
		BufferedWriter bw;
		RepeatAlphabet.Character tmpCharacter = new RepeatAlphabet.Character();
		int[] histogram;
		RepeatAlphabet.Character[] oldAlphabet, newAlphabet;
		
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		
		// Loading alphabets. The order of the two deserializations matters.
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_NEW,2);
		newAlphabet=RepeatAlphabet.alphabet;
		lastUnique_new=RepeatAlphabet.lastUnique;
		lastPeriodic_new=RepeatAlphabet.lastPeriodic;
		lastAlphabet_new=RepeatAlphabet.lastAlphabet;
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_OLD,2);
		oldAlphabet=RepeatAlphabet.alphabet;
		lastUnique_old=RepeatAlphabet.lastUnique;
		lastPeriodic_old=RepeatAlphabet.lastPeriodic;
		lastAlphabet_old=RepeatAlphabet.lastAlphabet;
		
		RepeatAlphabet.deserializeSpacers(SPACERS_FILE,N_SPACERS);
		histogram = new int[21];  // Arbitrary
		Math.set(histogram,0,histogram.length-1);
		br1 = new BufferedReader(new FileReader(READ_IDS_FILE));
		br2 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		br3 = new BufferedReader(new FileReader(READ2CHARACTERS_FILE_OLD));
		br4 = new BufferedReader(new FileReader(READ2BOUNDARIES_FILE_OLD));
		br5 = new BufferedReader(new FileReader(TANDEMS_FILE));
		bw = new BufferedWriter(new FileWriter(READ2CHARACTERS_FILE_NEW));
		i=0; j=0;
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); str4=br4.readLine(); str5=br5.readLine();
		while (str1!=null) {
			j=RepeatAlphabet.tandemSpacers_updateTranslation(Integer.parseInt(str1),Integer.parseInt(str2),j,str3,str4,str5,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,bw,DISTANCE_THRESHOLD,histogram,tmpCharacter);
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); str4=br4.readLine(); str5=br5.readLine();
		}
		br1.close(); br2.close(); br3.close(); br4.close(); br5.close(); bw.close();
		sum=0;
		for (i=0; i<histogram.length; i++) sum+=histogram[i];
		System.err.println("Fixed "+sum+" tandem spacers. Histogram of fixed tandem spacers lengths:");
		for (i=0; i<histogram.length; i++) System.err.println((i*IO.quantum)+": "+histogram[i]);
	}

}