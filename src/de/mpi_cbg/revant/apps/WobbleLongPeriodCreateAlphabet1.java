package de.mpi_cbg.revant.apps;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import java.io.*;

/**
 * Marks all characters of the alphabet that are adjacent to a long-period tandem or that
 * belong to a long-period tandem in some translation, and prints the corresponding
 * bitvector in output.
 *
 * This is designed to work on a chunk of reads.
 */
public class WobbleLongPeriodCreateAlphabet1 {

	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];  // Of all reads
		final String TRANSLATED_READS_CHARACTERS_FILE = args[1];  // Of a chunk of reads
        final String TANDEMS_FILE = args[2];  // Of a chunk of reads. Non-periodic only.
        final String READ_LENGTHS_FILE = args[3];  // Of a chunk of reads
		final String OUTPUT_FILE = args[4];
		
		int i;
		int nBlocks, nFlags;
		String str1, str2, str3;
        RepeatAlphabet.Character tmpCharacter;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		boolean[] flags;
		int[] tmpArray;
		
        tmpCharacter = new RepeatAlphabet.Character();
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		flags = new boolean[RepeatAlphabet.lastAlphabet+1];
		Math.set(flags,RepeatAlphabet.lastAlphabet,false);
		if (RepeatAlphabet.lastAlphabet>RepeatAlphabet.lastPeriodic) {
			tmpArray = new int[100];  // Arbitrary
			br1 = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
            br2 = new BufferedReader(new FileReader(TANDEMS_FILE));
            br3 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
            i=0;
			while (str1!=null) {
				nBlocks=1+((str1.length()+1)>>1);  // Loose upper bound
				if (tmpArray.length<nBlocks) tmpArray = new int[nBlocks];
                RepeatAlphabet.wobble_longPeriod_markAlphabet(str1,str2,Integer.parseInt(str3),flags,tmpCharacter,tmpArray);                
				i++;
				if (i%10000==0) System.err.println("Processed "+i+" reads");
				str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
			}
			br1.close(); br2.close(); br3.close();
		}
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) bw.write(flags[i]?"1\n":"0\n");
		bw.close();
	}

}