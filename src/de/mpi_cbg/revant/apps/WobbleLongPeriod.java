package de.mpi_cbg.revant.apps;

import de.mpi_cbg.revant.util.IO;
import java.io.*;

/**
 * Like $Wobble.java$.
 */
public class WobbleLongPeriod {

	public static void main(String[] args) throws IOException {
		final String TRANSLATED_READS_CHARACTERS_FILE = args[0];  // Of a chunk of reads
		final int WOBBLE_LENGTH = Integer.parseInt(args[1]);
		final String ALPHABET_FILE_OLD = args[2];  // Of all reads
		final String ALPHABET_FILE_NEW = args[3];  // Of all reads
		final String ALPHABET_FILE_OLD2NEW = args[4];  // Of all reads
		final String REPEAT_LENGTHS_FILE = args[5];
		final int N_REPEATS = Integer.parseInt(args[6]);
        final String TANDEMS_FILE = args[7];  // Of a chunk of reads. Non-periodic only.
		final String OUTPUT_FILE = args[8];  // Of a chunk of reads
		
		int i;
		int nBlocks, lastUnique_new, lastPeriodic_new, lastAlphabet_new;
		String str1, str2;
		BufferedReader br1, br2;
		BufferedWriter bw;
		int[] old2new, tmpArray1, tmpArray2, tmpArray3;
		RepeatAlphabet.Character[] alphabet_new;
		
		final int MULTIPLICITY = (1+WOBBLE_LENGTH/IO.quantum)*3;
		final int MAX_NEWCHARS_PER_CHAR = MULTIPLICITY*MULTIPLICITY;
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_NEW,2);
		alphabet_new=RepeatAlphabet.alphabet; lastUnique_new=RepeatAlphabet.lastUnique; lastPeriodic_new=RepeatAlphabet.lastPeriodic; lastAlphabet_new=RepeatAlphabet.lastAlphabet;
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_OLD,2);
		old2new = new int[lastAlphabet_new+1];
		br1 = new BufferedReader(new FileReader(ALPHABET_FILE_OLD2NEW));
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) old2new[i]=Integer.parseInt(br1.readLine());
		br1.close();
		tmpArray1 = new int[100];  // Arbitrary
		tmpArray2 = new int[(RepeatAlphabet.lastAlphabet+1)*MAX_NEWCHARS_PER_CHAR];
		tmpArray3 = new int[] {0,0};
		br1 = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
        br2 = new BufferedReader(new FileReader(TANDEMS_FILE));
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null) {
			nBlocks=1+((str1.length()+1)>>1);  // Loose upper bound
			if (tmpArray1.length<nBlocks) tmpArray1 = new int[nBlocks];
            RepeatAlphabet.wobble_longPeriod(str1,str2,WOBBLE_LENGTH,IO.quantum,old2new,alphabet_new,lastUnique_new,lastPeriodic_new,lastAlphabet_new,bw,tmpArray1,tmpArray2,tmpArray3);
			str1=br1.readLine(); str2=br2.readLine();
		}
		br1.close(); br2.close(); bw.close();
		System.err.println("Applied wobbling to "+tmpArray3[0]+" blocks out of "+tmpArray3[1]+" total ("+((100.0*tmpArray3[0])/tmpArray3[1])+"%)");
	}

}