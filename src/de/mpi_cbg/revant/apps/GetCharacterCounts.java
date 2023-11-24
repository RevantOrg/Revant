package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;


/**
 * Counts the full and partial frequency of every character in any orientation.
 */
public class GetCharacterCounts {
	
	public static void main(String[] args) throws IOException {
		final String TRANSLATED_FILE = args[0];
        final String BOUNDARIES_FILE = args[1];
        final String READ_LENGTHS_FILE = args[2];
		final String ALPHABET_FILE = args[3];
		final String COUNTS_FILE = args[4];
		final String HISTOGRAM_FILE = args[5];
		
		final int MAX_HISTOGRAM_FREQUENCY = 100000;
		
		int i, j;
		int row;	
		String str1, str2, str3;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		boolean[] marked;
		long[][] characterCount;
		long[][] characterHistogram;
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		
		// Collecting character counts
		characterCount = new long[RepeatAlphabet.lastAlphabet+2][2];
		Math.set(characterCount,0);
		br1 = new BufferedReader(new FileReader(TRANSLATED_FILE));
        br2 = new BufferedReader(new FileReader(BOUNDARIES_FILE));
        br3 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); row=0;
		while (str1!=null) {
			if (row%10000==0) System.err.println("Processed "+row+" reads");
			RepeatAlphabet.incrementCharacterCounts(str1,str2,Integer.parseInt(str3),characterCount);
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); row++;
		}
		br1.close(); br2.close(); br3.close();
		marked = new boolean[RepeatAlphabet.lastAlphabet+1];
		RepeatAlphabet.symmetrizeCharacterCounts(characterCount,marked);
		bw = new BufferedWriter(new FileWriter(COUNTS_FILE));
		for (i=0; i<characterCount.length; i++) bw.write(characterCount[i][0]+","+characterCount[i][1]+"\n");
		bw.close();
		
		// Building count histogram
		characterHistogram=RepeatAlphabet.getCharacterHistogram(characterCount,2,marked,MAX_HISTOGRAM_FREQUENCY);
		bw = new BufferedWriter(new FileWriter(HISTOGRAM_FILE));
		bw.write("# Histogram of symmetrized character frequencies (full+partial) \n");
		for (i=0; i<characterHistogram.length; i++) {
			bw.write(i+","+characterHistogram[i][0]);
			for (j=1; j<characterHistogram[i].length; j++) bw.write(","+characterHistogram[i][j]);
			bw.newLine();
		}
		bw.close();
	}

}