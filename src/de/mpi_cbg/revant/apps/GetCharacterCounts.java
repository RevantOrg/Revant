package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Hashtable;

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
		bw = new BufferedWriter(new FileWriter(COUNTS_FILE));
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) bw.write(characterCount[i]+"\n");
		bw.close();
		
		// Building count histogram
		characterHistogram=getCharacterHistogram(characterCount,MAX_HISTOGRAM_FREQUENCY);
		bw = new BufferedWriter(new FileWriter(HISTOGRAM_FILE));
		bw.write("# Histogram of character frequencies \n");
		for (i=0; i<characterHistogram.length; i++) {
			bw.write(i+","+characterHistogram[i][0]);
			for (j=1; j<characterHistogram[i].length; j++) bw.write(","+characterHistogram[i][j]);
			bw.newLine();
		}
		bw.close();
	}
	
	
	/**
	 * Builds a histogram of $characterCount$ values. Rows: counts. Columns:
	 * 0: unique, open; 
	 * 1: unique, closed;
	 * 2: periodic, open;
	 * 3: periodic, closed; 
	 * 4: nonperiodic, open;
	 * 5: nonperiodic, closed.
	 */
	private static final long[][] getCharacterHistogram(long[] characterCount, int maxFrequency) throws IOException {
		int i;
		int nElements;
		long count;
		RepeatAlphabet.Character character;
		String key;
		Long value;
		Hashtable<String,Long> hashTable;
		Long[] values;
		long[][] characterHistogram;

		// Unique, open and closed.
		characterHistogram = new long[maxFrequency+1][3<<1];
		Math.set(characterHistogram,0);
		for (i=0; i<=RepeatAlphabet.lastUnique; i++) {
			count=characterCount[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][RepeatAlphabet.alphabet[i].isOpen()?0:1]++;
		}
		
		// Periodic, open.
		hashTable = new Hashtable<String,Long>();
		for (i=RepeatAlphabet.lastUnique+1; i<=RepeatAlphabet.lastPeriodic; i++) {
			character=RepeatAlphabet.alphabet[i];
			if (!character.isOpen()) continue;
			key=character.repeat+","+character.length;
			value=hashTable.get(key);
			count=characterCount[i];
			if (value==null) hashTable.put(key,Long.valueOf(count));
			else hashTable.put(key,Long.valueOf(value.longValue()+count));
		}
		values = new Long[hashTable.size()];
		hashTable.values().toArray(values);
		for (i=0; i<values.length; i++) {
			count=values[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][2]++;
		}
		
		// Periodic, closed.
		hashTable.clear();
		for (i=RepeatAlphabet.lastUnique+1; i<=RepeatAlphabet.lastPeriodic; i++) {
			character=RepeatAlphabet.alphabet[i];
			if (character.isOpen()) continue;
			key=character.repeat+","+character.length;
			value=hashTable.get(key);
			count=characterCount[i];
			if (value==null) hashTable.put(key,Long.valueOf(count));
			else hashTable.put(key,Long.valueOf(value.longValue()+count));
		}
		nElements=hashTable.size();
		if (values.length<nElements) values = new Long[nElements];
		hashTable.values().toArray(values);
		for (i=0; i<nElements; i++) {
			count=values[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][3]++;
		}
		
		// Nonperiodic, open.
		hashTable.clear();
		for (i=RepeatAlphabet.lastPeriodic+1; i<=RepeatAlphabet.lastAlphabet; i++) {
			character=RepeatAlphabet.alphabet[i];
			if (!character.isOpen()) continue;
			key=character.repeat+","+character.start+","+character.end;
			value=hashTable.get(key);
			count=characterCount[i];
			if (value==null) hashTable.put(key,Long.valueOf(count));
			else hashTable.put(key,Long.valueOf(value.longValue()+count));
		}
		nElements=hashTable.size();
		if (values.length<nElements) values = new Long[nElements];
		hashTable.values().toArray(values);
		for (i=0; i<nElements; i++) {
			count=values[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][4]++;
		}
		
		// Nonperiodic, closed.
		hashTable.clear();
		for (i=RepeatAlphabet.lastPeriodic+1; i<=RepeatAlphabet.lastAlphabet; i++) {
			character=RepeatAlphabet.alphabet[i];
			if (character.isOpen()) continue;
			key=character.repeat+","+character.start+","+character.end;
			value=hashTable.get(key);
			count=characterCount[i];
			if (value==null) hashTable.put(key,Long.valueOf(count));
			else hashTable.put(key,Long.valueOf(value.longValue()+count));
		}
		nElements=hashTable.size();
		if (values.length<nElements) values = new Long[nElements];
		hashTable.values().toArray(values);
		for (i=0; i<nElements; i++) {
			count=values[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][5]++;
		}
		
		return characterHistogram;
	}

}