package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;

/**
 * Concatenates characters in adjacent non-periodic blocks, that are also adjacent on
 * their repeat. This is designed to work on a chunk of reads.
 */
public class ConcatenateBlocks1 {
	/**
	 * @param args
	 * 
	 * 6: prints to this file all instances of new characters created by concatenation,
	 * and all instances of old characters that are still used after concatenation. The
	 * output is not sorted and may contain duplicates.
	 */
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];  // Of all reads
		final String READ_LENGTHS_FILE = args[1];  // Of a chunk of reads
		final String TRANSLATED_READS_CHARACTERS_FILE = args[2];  // Of a chunk of reads
		final String TRANSLATED_READS_BOUNDARIES_FILE = args[3];  // Of a chunk of reads
		final String REPEAT_LENGTHS_FILE = args[4];
		final int N_REPEATS = Integer.parseInt(args[5]);
		final int DISTANCE_THRESHOLD = Integer.parseInt(args[6]);
		final String OUTPUT_FILE = args[7];  // Of a chunk of reads
		
		final int CAPACITY = 100;  // Arbitrary
		final int QUANTUM = IO.quantum;
		
		int i;
		int nBlocks;
		String str1, str2, str3;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		RepeatAlphabet.Character tmpCharacter = new RepeatAlphabet.Character();
		boolean[] used, tmpBoolean1, tmpBoolean2, tmpBoolean3;
		int[] tmpArray = new int[2];
		
		RepeatAlphabet.loadRepeatLengths(REPEAT_LENGTHS_FILE,N_REPEATS);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		used = new boolean[RepeatAlphabet.lastAlphabet+2];
		Math.set(used,used.length-1,false);
		br1 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		br2 = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
		br3 = new BufferedReader(new FileReader(TRANSLATED_READS_BOUNDARIES_FILE));
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		tmpBoolean1 = new boolean[CAPACITY]; 
		tmpBoolean2 = new boolean[CAPACITY];
		tmpBoolean3 = new boolean[CAPACITY];
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); i=0;
		while (str1!=null) {
            if (str2.length()>0) {
    			nBlocks=1+((str3.length()+1)>>1);  // Loose upper bound
    			if (tmpBoolean1.length<nBlocks) {
    				tmpBoolean1 = new boolean[nBlocks];
    				tmpBoolean2 = new boolean[nBlocks];
    				tmpBoolean3 = new boolean[nBlocks];
    			}
    			RepeatAlphabet.concatenateBlocks_collectCharacterInstances(str2,str3,Integer.parseInt(str1),DISTANCE_THRESHOLD,QUANTUM,used,bw,tmpCharacter,tmpBoolean1,tmpBoolean2,tmpBoolean3,tmpArray);
            }
			if (i%10000==0) System.err.println("Processed "+i+" reads");
			i++;
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
		}
		br1.close(); br2.close(); br3.close();
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) {
			if (used[i]) {
				bw.write(RepeatAlphabet.alphabet[i].toString());
				bw.newLine();
			}
		}
		bw.close();
	}

}