package de.mpi_cbg.revant.apps;

import java.io.*;


public class MergeAlphabetHeaders {
	
	public static void main(String[] args) throws IOException {
		final String HEADERS_DIR = args[0];
		final String HEADERS_PREFIX = args[1];
		final String OUTPUT_FILE = args[2];
		
		int i, p, q, r;
		int nUnique, nPeriodic, nCharacters, maxOpenLength_nonperiodic;
		final int headersPrefixLength = HEADERS_PREFIX.length();
		String str;
		BufferedReader br;
		BufferedWriter bw;
		String[] files, tokens;
		
		files=(new File(HEADERS_DIR)).list();
		nUnique=0; nPeriodic=0; nCharacters=0; maxOpenLength_nonperiodic=0;
		for (i=0; i<files.length; i++) {
			if (files[i].length()>headersPrefixLength && files[i].substring(0,headersPrefixLength).equals(HEADERS_PREFIX)) {
				br = new BufferedReader(new FileReader(HEADERS_DIR+"/"+files[i]));
				str=br.readLine();
				br.close();
				tokens=str.split(RepeatAlphabet.SEPARATOR_MINOR);
				p=Integer.parseInt(tokens[0]);
				q=Integer.parseInt(tokens[1]);
				r=Integer.parseInt(tokens[2]);
				nUnique+=p+1; nPeriodic+=q-p; nCharacters+=r+1;
				maxOpenLength_nonperiodic=Math.max(maxOpenLength_nonperiodic,Integer.parseInt(tokens[3]));
			}
		}
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		bw.write((nUnique-1)+RepeatAlphabet.SEPARATOR_MINOR+(nUnique+nPeriodic-1)+RepeatAlphabet.SEPARATOR_MINOR+(nCharacters-1)+RepeatAlphabet.SEPARATOR_MINOR+maxOpenLength_nonperiodic);
		bw.newLine(); bw.close();
	}

}