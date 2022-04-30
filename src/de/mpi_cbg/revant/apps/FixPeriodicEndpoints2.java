package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;

/**
 * Uses spacers to update the alphabet. This is designed to work on a chunk of reads.
 */
public class FixPeriodicEndpoints2 {
	/**
	 * @param args
	 * 
	 * 8: prints to this file all instances of new characters created by spacers, and all 
	 * instances of old characters that are still used after spacers have been applied.
	 * The output is not sorted and may contain duplicates.
	 */
	public static void main(String[] args) throws IOException {
		final int MAX_SPACER_LENGTH = Integer.parseInt(args[0]);
		final String SPACERS_FILE = args[1];  // Of a chunk of reads
		final int N_SPACERS = Integer.parseInt(args[2]);
		final String ALPHABET_FILE = args[3];  // Of all reads
		final String READ_IDS_FILE = args[4];  // Of a chunk of reads
		final String READ_LENGTHS_FILE = args[5];  // Of a chunk of reads
		final String TRANSLATED_READS_CHARACTERS_FILE = args[6];  // Of a chunk of reads
		final String TRANSLATED_READS_BOUNDARIES_FILE = args[7];  // Of a chunk of reads
		final String OUTPUT_FILE = args[8];  // Of a chunk of reads
		final String OUTPUT_FILE_UNIQUE = args[9];  // Of a chunk of reads
		
		final int CAPACITY = 100;  // Arbitrary
		
		int i, j;
		int nBlocks;
		String str1, str2, str3, str4;
		BufferedReader br1, br2, br3, br4;
		BufferedWriter bw;
		boolean[] used, isBlockPeriodic, isBlockNonperiodic;
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.deserializeSpacers(SPACERS_FILE,N_SPACERS);
		used = new boolean[RepeatAlphabet.lastAlphabet+2];
		Math.set(used,used.length-1,false);
		isBlockPeriodic = new boolean[CAPACITY];
		isBlockNonperiodic = new boolean[CAPACITY];
		br1 = new BufferedReader(new FileReader(READ_IDS_FILE));
		br2 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		br3 = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
		br4 = new BufferedReader(new FileReader(TRANSLATED_READS_BOUNDARIES_FILE));
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		i=0; j=0;
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); str4=br4.readLine();
		while (str1!=null) {
			nBlocks=1+((str3.length()+1)>>1);  // Loose upper bound
			if (isBlockPeriodic.length<nBlocks) {
				isBlockPeriodic = new boolean[nBlocks];
				isBlockNonperiodic = new boolean[nBlocks];
			}
			j=RepeatAlphabet.fixPeriodicEndpoints_collectCharacterInstances(Integer.parseInt(str1),j,str3,str4,Integer.parseInt(str2),MAX_SPACER_LENGTH,bw,used,isBlockPeriodic,isBlockNonperiodic);
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); str4=br4.readLine();
		}
		br1.close(); br2.close(); br3.close(); br4.close();
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) {
			if (used[i]) {
				bw.write(RepeatAlphabet.alphabet[i].toString());
				bw.newLine();
			}
		}
		bw.close();
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_UNIQUE));
		bw.write(RepeatAlphabet.maxOpenLength_unique+"");
		bw.close();
	}

}