package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Hashtable;

import de.mpi_cbg.revant.util.Math;


public class CollectKmers {
	
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String TRANSLATED_FILE = args[1];
		final int K = Integer.parseInt(args[2]);
		final String OUTPUT_DIR = args[3];
		
		final int MAX_FREQUENCY = 100000;
		
		int i, j;
		int row, length, nKmers, nElements;
		long count, max;
		Long value;
		String str, key;
		StringBuilder sb;
		BufferedReader br;
		BufferedWriter bw;
		RepeatAlphabet.Character character;
		Hashtable<String,Long> hashTable;
		long[] characterCount, histogram;
		long[][] characterHistogram;
		Long[] values;
		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE);
		
		System.err.println("Collecting character statistics...");		
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
		// Unique, open and closed.
		characterHistogram = new long[MAX_FREQUENCY+1][3<<1];
		Math.set(characterHistogram,0);
		for (i=0; i<=RepeatAlphabet.lastUnique; i++) {
			count=characterCount[i];
			if (count>MAX_FREQUENCY) count=MAX_FREQUENCY;
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
			if (count>MAX_FREQUENCY) count=MAX_FREQUENCY;
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
			if (count>MAX_FREQUENCY) count=MAX_FREQUENCY;
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
			if (count>MAX_FREQUENCY) count=MAX_FREQUENCY;
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
			if (count>MAX_FREQUENCY) count=MAX_FREQUENCY;
			characterHistogram[(int)count][5]++;
		}
		// Outputting
		bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/histogram-characters.txt"));
		bw.write("# Histogram of character frequencies \n");
		for (i=0; i<characterHistogram.length; i++) {
			bw.write(i+","+characterHistogram[i][0]);
			for (j=1; j<characterHistogram[i].length; j++) bw.write(","+characterHistogram[i][j]);
			bw.newLine();
		}
		bw.close();
		bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/alphabet-counts.txt"));
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) bw.write(characterCount[i]+"\n");
		bw.close();		
		
		
		
		
		
		
		
		
		
		
System.exit(1);		
		

		System.err.println("Collecting "+K+"-mers...");
		sb = new StringBuilder();
		hashTable = new Hashtable<String,Long>();
		br = new BufferedReader(new FileReader(TRANSLATED_FILE));
		str=br.readLine(); row=0;
		while (str!=null) {
			if (row%1000==0) System.err.println("Processed "+row+" reads, "+hashTable.size()+" "+K+"-mers.");
			RepeatAlphabet.getKmers(str,K,hashTable,sb);
			str=br.readLine(); row++;
		}
		br.close(); nKmers=hashTable.size();
		System.err.println("DONE  "+nKmers+" distinct "+K+"-mers");
		
		System.err.println("Computing frequency histogram...");
		values = new Long[nKmers];
		hashTable.values().toArray(values);
		length=values.length; max=0;
		for (i=0; i<length; i++) {
			count=values[i].longValue();
			if (count>max) max=count;
		}
		histogram = new long[(int)max+1];
		Math.set(histogram,(int)max,0);
		for (i=0; i<length; i++) histogram[(int)(values[i].longValue())]++;
		System.err.println(K+"-mer histogram:");
		for (i=0; i<=max; i++) System.err.println(i+","+histogram[i]);
	}

}