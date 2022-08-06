package de.mpi_cbg.revant.apps;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import java.io.*;

/**
 * Marks all characters of the alphabet that are adjacent to periodic characters in some 
 * translation, and prints the corresponding bitvector in output.
 *
 * This is designed to work on a chunk of reads.
 */
public class WobbleCreateAlphabet1 {

	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];  // Of all reads
		final String TRANSLATED_READS_CHARACTERS_FILE = args[1];  // Of a chunk of reads
		final String OUTPUT_FILE = args[2];
		
		int i;
		int nBlocks, nFlags;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		boolean[] flags;
		int[] tmpArray;
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		flags = new boolean[RepeatAlphabet.lastAlphabet+1];
		Math.set(flags,RepeatAlphabet.lastAlphabet,false);
		tmpArray = new int[100];  // Arbitrary
		br = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
		str=br.readLine(); i=0;
		while (str!=null) {
			nBlocks=1+((str.length()+1)>>1);  // Loose upper bound
			if (tmpArray.length<nBlocks) tmpArray = new int[nBlocks];
			RepeatAlphabet.wobble_markAlphabet(str,flags,tmpArray);
			i++;
			if (i%10000==0) System.err.println("Processed "+i+" reads");
			str=br.readLine();
		}
		br.close();
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) bw.write(flags[i]?"1\n":"0\n");
		bw.close();
	}

}