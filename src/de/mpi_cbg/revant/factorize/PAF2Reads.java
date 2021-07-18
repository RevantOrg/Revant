package de.mpi_cbg.revant.factorize;

import java.util.*;
import java.io.*;

/**
 * Given a multi-FASTA file, and a mapping from miniasm read IDs to integers, the program 
 * creates one file per read in the mapping, with all the sequence on a single line.
 */
public class PAF2Reads {

	
	public static void main(String[] args) throws IOException {
		final String OLDNAMES_NEWNAMES_FILE = args[0];
		final int N_READS = Integer.parseInt(args[1]);
		final String FASTA_FILE = args[2];
		final String OUTPUT_DIR = args[3];
		
		int p;
		String str, key;
		HashMap<String,Integer> sequenceNames;
		BufferedReader br;
		BufferedWriter bw;
		
		System.err.print("Loading names map... ");
		sequenceNames = new HashMap<String,Integer>(N_READS);
		br = new BufferedReader(new FileReader(OLDNAMES_NEWNAMES_FILE));
		str=br.readLine(); 
		while (str!=null) {
			p=str.indexOf(",");
			sequenceNames.put(str.substring(0,p),Integer.valueOf(Integer.parseInt(str.substring(p+1))));
			str=br.readLine();
		}
		br.close();
		System.err.println(" DONE");
		
		System.err.print("Building reads... ");
		br = new BufferedReader(new FileReader(FASTA_FILE));
		bw=null; str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)=='>') {
				if (bw!=null) { bw.close(); bw=null; }
				key=str.substring(1);
				if (sequenceNames.containsKey(key)) bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/read"+sequenceNames.get(key).intValue()+".txt"));
			}
			else if (bw!=null) bw.write(str);
			str=br.readLine();
		}
		if (bw!=null) bw.close();
		br.close();
		System.err.println(" DONE");
	}
	

}