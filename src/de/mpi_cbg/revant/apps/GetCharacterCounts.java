package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.Math;


/**
 * Counts the frequency of every character in both orientations
 */
public class GetCharacterCounts {
	
	public static void main(String[] args) throws IOException {
		final String TRANSLATED_FILE = args[0];
		final String ALPHABET_FILE = args[1];
		final String COUNTS_FILE = args[2];
		final String HISTOGRAM_FILE = args[3];
		
		final int MAX_HISTOGRAM_FREQUENCY = 100000;
		
		int i, j;
		int row;	
		String str;
		BufferedReader br;
		BufferedWriter bw;
		boolean[] marked;
		long[] characterCount;
		long[][] characterHistogram;
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		
		// Collecting character counts
		characterCount = new long[RepeatAlphabet.lastAlphabet+2];
		Math.set(characterCount,characterCount.length-1,0);
		br = new BufferedReader(new FileReader(TRANSLATED_FILE));
		str=br.readLine(); row=0;
		while (str!=null) {
			if (row%10000==0) System.err.println("Processed "+row+" reads");
			RepeatAlphabet.incrementCharacterCounts(str,characterCount);
			str=br.readLine(); row++;
		}
		br.close();
		marked = new boolean[RepeatAlphabet.lastAlphabet+1];
		RepeatAlphabet.symmetrizeCharacterCounts(characterCount,marked);
		bw = new BufferedWriter(new FileWriter(COUNTS_FILE));
		for (i=0; i<characterCount.length; i++) bw.write(characterCount[i]+"\n");
		bw.close();
		
		// Building count histogram
		characterHistogram=RepeatAlphabet.getCharacterHistogram(characterCount,marked,MAX_HISTOGRAM_FREQUENCY);
		bw = new BufferedWriter(new FileWriter(HISTOGRAM_FILE));
		bw.write("# Histogram of symmetrized character frequencies \n");
		for (i=0; i<characterHistogram.length; i++) {
			bw.write(i+","+characterHistogram[i][0]);
			for (j=1; j<characterHistogram[i].length; j++) bw.write(","+characterHistogram[i][j]);
			bw.newLine();
		}
		bw.close();
	}

}