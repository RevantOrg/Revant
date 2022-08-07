package de.mpi_cbg.revant.apps;

import de.mpi_cbg.revant.util.IO;
import java.io.*;

/**
 * Adds to some translated blocks the additional characters that result from wobbling its
 * characters.
 * 
 * This is designed to work on a chunk of reads.
 */
public class Wobble {

	public static void main(String[] args) throws IOException {
		final String TRANSLATED_READS_CHARACTERS_FILE = args[0];  // Of a chunk of reads
		final int WOBBLE_LENGTH = Integer.parseInt(args[1]);
		final String ALPHABET_FILE_OLD = args[2];  // Of all reads
		final String ALPHABET_FILE_NEW = args[3];  // Of all reads
		final String ALPHABET_FILE_OLD2NEW = args[4];  // Of all reads
		final String REPEAT_LENGTHS_FILE = args[5];
		final int N_REPEATS = Integer.parseInt(args[6]);
		final String OUTPUT_FILE = args[7];  // Of a chunk of reads
		
		int i;
		int nBlocks, lastUnique_new, lastPeriodic_new, lastAlphabet_new;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		int[] old2new, tmpArray1, tmpArray2, tmpArray3;
		RepeatAlphabet.Character[] alphabet_new;
		
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_NEW,2);
		alphabet_new=RepeatAlphabet.alphabet; lastUnique_new=RepeatAlphabet.lastUnique; lastPeriodic_new=RepeatAlphabet.lastPeriodic; lastAlphabet_new=RepeatAlphabet.lastAlphabet;
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_OLD,2);
		old2new = new int[lastAlphabet_new+1];
		br = new BufferedReader(new FileReader(ALPHABET_FILE_OLD2NEW));
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) old2new[i]=Integer.parseInt(br.readLine());
		br.close();
		tmpArray1 = new int[100];  // Arbitrary
		tmpArray2 = new int[RepeatAlphabet.lastAlphabet+1];
		tmpArray3 = new int[] {0,0};
		br = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		str=br.readLine();
		while (str!=null) {
			nBlocks=1+((str.length()+1)>>1);  // Loose upper bound
			if (tmpArray1.length<nBlocks) tmpArray1 = new int[nBlocks];
			RepeatAlphabet.wobble(str,WOBBLE_LENGTH,IO.quantum,old2new,alphabet_new,lastUnique_new,lastPeriodic_new,lastAlphabet_new,bw,tmpArray1,tmpArray2,tmpArray3);		
			str=br.readLine();
		}
		br.close(); bw.close();
		System.err.println("Applied wobbling to "+tmpArray3[0]+" blocks out of "+tmpArray3[1]+" total ("+((100.0*tmpArray3[0])/tmpArray3[1])+"%)");
	}

}