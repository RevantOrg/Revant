package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Hashtable;

import de.mpi_cbg.revant.util.Math;


public class CollectKmers {
	
	public static void main(String[] args) throws IOException {
		final String TRANSLATED_FILE = args[0];
		final int K = Integer.parseInt(args[1]);
		
		int i, j;
		int row, max, length, nKmers;
		String str;
		StringBuilder sb;
		BufferedReader br;
		Hashtable<String,Integer> kmers;
		int[] histogram;
		Integer[] values;
		
		
		System.err.println("Collecting "+K+"-mers...");
		sb = new StringBuilder();
		kmers = new Hashtable<String,Integer>();
		br = new BufferedReader(new FileReader(TRANSLATED_FILE));
		str=br.readLine(); row=0;
		while (str!=null) {
			if (row%1000==0) System.err.println("Processed "+row+" reads, "+kmers.size()+" "+K+"-mers.");
			RepeatAlphabet.loadTranslatedRead(str,0,K,kmers,sb);
			str=br.readLine(); row++;
		}
		br.close(); nKmers=kmers.size();
		System.err.println("DONE  "+nKmers+" distinct "+K+"-mers");
		
		System.err.println("Computing frequency histogram...");
		values = new Integer[nKmers];
		kmers.values().toArray(values);
		length=values.length; max=0;
		for (i=0; i<length; i++) {
			j=values[i].intValue();
			max=Math.max(max,j);
		}
		histogram = new int[max+1];
		Math.set(histogram,max,0);
		for (i=0; i<length; i++) histogram[values[i].intValue()]++;
		System.err.println(K+"-mer histogram:");
		for (i=0; i<=max; i++) System.err.println(i+","+histogram[i]);
	}

}