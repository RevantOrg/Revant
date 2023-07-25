package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;

/**
 * Uses the new alphabet induced by concatenation, to update the translation of every
 * read. This is designed to work on a chunk of reads.
 */
public class ConcatenateBlocks2 {

	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE_OLD = args[0];  // Of all reads
		final String ALPHABET_FILE_NEW = args[1];  // Of all reads
		final String READ_LENGTHS_FILE = args[2];  // Of a chunk of reads
		final String REPEAT_LENGTHS_FILE = args[3];
		final int N_REPEATS = Integer.parseInt(args[4]);
		final String READ2CHARACTERS_FILE_OLD = args[5];  // Of a chunk of reads
		final String READ2BOUNDARIES_FILE_OLD = args[6];  // Of a chunk of reads
		final int DISTANCE_THRESHOLD = Integer.parseInt(args[7]);
		final String READ2CHARACTERS_FILE_NEW = args[8];  // Of a chunk of reads
		final String READ2BOUNDARIES_FILE_NEW = args[9];  // Of a chunk of reads
		final String FULLYCONTAINED_FILE_NEW = args[10];  // Of a chunk of reads
		
		final int CAPACITY = 100;  // Arbitrary
		
		int i;
		int lastUnique_old, lastPeriodic_old, lastAlphabet_old, lastUnique_new, lastPeriodic_new, lastAlphabet_new, nBlocks;
		int nBlocks_total, nBlocks_concatenated;
		String str1, str2, str3;
		BufferedReader br1, br2, br3;
		BufferedWriter bw1, bw2, bw3;
		RepeatAlphabet.Character tmpCharacter = new RepeatAlphabet.Character();
		boolean[] tmpBoolean1, tmpBoolean2;
		int[] tmpArray = new int[2];
		int[] stats = new int[2];
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
		
		tmpBoolean1 = new boolean[CAPACITY]; 
		tmpBoolean2 = new boolean[CAPACITY];
		br1 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		br2 = new BufferedReader(new FileReader(READ2CHARACTERS_FILE_OLD));
		br3 = new BufferedReader(new FileReader(READ2BOUNDARIES_FILE_OLD));
		bw1 = new BufferedWriter(new FileWriter(READ2CHARACTERS_FILE_NEW));
		bw2 = new BufferedWriter(new FileWriter(READ2BOUNDARIES_FILE_NEW));
		bw3 = new BufferedWriter(new FileWriter(FULLYCONTAINED_FILE_NEW));
		nBlocks_total=0; nBlocks_concatenated=0;
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); i=0;
		while (str1!=null) {
            if (str2.length()>0) {
    			nBlocks=1+((str3.length()+1)>>1);  // Loose upper bound
    			if (tmpBoolean1.length<nBlocks) {
    				tmpBoolean1 = new boolean[nBlocks];
    				tmpBoolean2 = new boolean[nBlocks];
    			}
    			RepeatAlphabet.concatenateBlocks_updateTranslation(i,Integer.parseInt(str1),str2,str3,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,bw1,bw2,bw3,DISTANCE_THRESHOLD,stats,tmpCharacter,tmpArray,tmpBoolean1,tmpBoolean2);
    			nBlocks_concatenated+=stats[0]; nBlocks_total+=stats[1];
            }
            else { bw1.newLine(); bw2.newLine(); }
			if (i%10000==0) System.err.println("Processed "+i+" reads");
			i++;
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
		}
		br1.close(); br2.close(); br3.close(); bw1.close(); bw2.close(); bw3.close();
		System.err.println(nBlocks_concatenated+" blocks involved in a concatenation out of "+nBlocks_total+" total blocks ("+((100.0*nBlocks_concatenated)/nBlocks_total)+"%).");
	}

}