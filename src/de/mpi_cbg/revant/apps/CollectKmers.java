package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Hashtable;

import de.mpi_cbg.revant.util.Math;


public class CollectKmers {
	
	public static void main(String[] args) throws IOException {
		final int K = Integer.parseInt(args[0]);
		final String TRANSLATED_FILE = args[1];
		final String ALPHABET_FILE = args[2];
		final int UNIQUE_MODE = Integer.parseInt(args[3]);
		final boolean OPEN_MODE = Integer.parseInt(args[4])==1;
		final boolean MULTI_MODE = Integer.parseInt(args[5])==1;
		final int MAX_HISTOGRAM_FREQUENCY = Integer.parseInt(args[6]);
		final String OUTPUT_DIR = args[7];
		
		final String HISTOGRAM_FILE = OUTPUT_DIR+"/histogram-k"+K+".txt";
		
		int i;
		int row, nKmers, intCount;
		long count;
		Long value;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		RepeatAlphabet.Kmer tmpKmer = new RepeatAlphabet.Kmer();
		Hashtable<RepeatAlphabet.Kmer,Long> kmers;
		int[] tmpArray2 = new int[K];
		int[] tmpArray3 = new int[2*K];
		int[] histogram;
		Long[] values;
		
		System.err.println("Collecting "+K+"-mers...");
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		kmers = new Hashtable<RepeatAlphabet.Kmer,Long>();
		br = new BufferedReader(new FileReader(TRANSLATED_FILE));
		str=br.readLine(); row=0;
		while (str!=null) {
			if (row%10000==0) System.err.println("Processed "+row+" reads, "+kmers.size()+" distinct "+K+"-mers.");
			RepeatAlphabet.getKmers(str,K,UNIQUE_MODE,OPEN_MODE,MULTI_MODE,kmers,tmpKmer,tmpArray2,tmpArray3);
			str=br.readLine(); row++;
		}
		br.close(); nKmers=kmers.size();
		System.err.println("DONE  "+nKmers+" total distinct "+K+"-mers (uniqueMode="+UNIQUE_MODE+", openMode="+OPEN_MODE+", multiMode="+MULTI_MODE+")");
		
		System.err.println("Computing frequency histogram... ");
		values = new Long[nKmers];
		kmers.values().toArray(values);
		histogram = new int[MAX_HISTOGRAM_FREQUENCY+1];
		Math.set(histogram,MAX_HISTOGRAM_FREQUENCY,0);
		for (i=0; i<nKmers; i++) {
			count=values[i].longValue();
			intCount=count<=MAX_HISTOGRAM_FREQUENCY?((int)count):MAX_HISTOGRAM_FREQUENCY;
			histogram[intCount]++;
		}
		bw = new BufferedWriter(new FileWriter(HISTOGRAM_FILE));
		for (i=0; i<=MAX_HISTOGRAM_FREQUENCY; i++) bw.write(i+","+histogram[i]+"\n");
		bw.close();
		System.err.println("DONE");
	}

}