package de.mpi_cbg.revant.apps;

import java.io.*;


public class MergeAlphabetHeaders {
	
	public static void main(String[] args) throws IOException {
		final String HEADERS_DIR = args[0];
		final String HEADERS_PREFIX = args[1];
		final String OUTPUT_FILE = args[2];
		
		int i, j;
		int nNonrepetitive, nPeriodic, nCharacters, maxOpenLength_nonperiodic;
		final int headersPrefixLength = HEADERS_PREFIX.length();
		String str;
		BufferedReader br;
		BufferedWriter bw;
		String[] files, tokens;
		
		files=(new File(HEADERS_DIR)).list();
		nNonrepetitive=0; nPeriodic=0; nCharacters=0; maxOpenLength_nonperiodic=0;
		for (i=0; i<files.length; i++) {
			if (files[i].length()>headersPrefixLength && files[i].substring(0,headersPrefixLength).equals(HEADERS_PREFIX)) {
				br = new BufferedReader(new FileReader(HEADERS_DIR+"/"+files[i]));
				str=br.readLine();
				br.close();
				tokens=str.split(RepeatAlphabet.SEPARATOR_MINOR);
				j=Integer.parseInt(tokens[0]);
				nNonrepetitive+=j+1;
				j=Integer.parseInt(tokens[1]);
				nPeriodic+=j+1;
				j=Integer.parseInt(tokens[2]);
				nCharacters+=j+1;
				j=Integer.parseInt(tokens[3]);
				if (j>maxOpenLength_nonperiodic) maxOpenLength_nonperiodic=j;
			}
		}
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		bw.write((nNonrepetitive-1)+RepeatAlphabet.SEPARATOR_MINOR+(nNonrepetitive+nPeriodic-1)+RepeatAlphabet.SEPARATOR_MINOR+(nCharacters-1)+RepeatAlphabet.SEPARATOR_MINOR+maxOpenLength_nonperiodic);
		bw.newLine(); bw.close();
	}

}