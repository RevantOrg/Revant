package de.mpi_cbg.revant.apps;

import java.io.*;


public class MergeAlphabetHeaders {
	
	public static void main(String[] args) throws IOException {
		final String HEADERS_DIR = args[0];
		final String HEADERS_PREFIX = args[1];
		final String HEADERS_PREFIX_UNIQUE = args[2];
		final String OUTPUT_FILE = args[3];
		
		int i, p, q, r;
		int nUnique, nPeriodic, nCharacters, maxOpenLength_unique;
		final int headersPrefixLength = HEADERS_PREFIX.length();
		final int headersPrefixLengthUnique = HEADERS_PREFIX_UNIQUE.length();
		String str;
		BufferedReader br;
		BufferedWriter bw;
		String[] files, tokens;
		
		files=(new File(HEADERS_DIR)).list();
		nUnique=0; nPeriodic=0; nCharacters=0; maxOpenLength_unique=0;
		for (i=0; i<files.length; i++) {
			if (files[i].length()>headersPrefixLength && files[i].substring(0,headersPrefixLength).equals(HEADERS_PREFIX)) {
				br = new BufferedReader(new FileReader(HEADERS_DIR+"/"+files[i]));
				str=br.readLine();
				br.close();
				tokens=str.split(RepeatAlphabet.SEPARATOR_MINOR+"");
				p=Integer.parseInt(tokens[0]);
				q=Integer.parseInt(tokens[1]);
				r=Integer.parseInt(tokens[2]);
				nUnique+=p+1; nPeriodic+=q-p; nCharacters+=r+1;
			}
			else if (files[i].length()>headersPrefixLengthUnique && files[i].substring(0,headersPrefixLengthUnique).equals(HEADERS_PREFIX_UNIQUE)) {
				br = new BufferedReader(new FileReader(HEADERS_DIR+"/"+files[i]));
				str=br.readLine();
				br.close();
System.err.println("VITTU> read "+str+" from file "+files[i]);				
				maxOpenLength_unique=Math.max(maxOpenLength_unique,Integer.parseInt(str));
			}
		}
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		bw.write((nUnique-1)+(RepeatAlphabet.SEPARATOR_MINOR+"")+(nUnique+nPeriodic-1)+(RepeatAlphabet.SEPARATOR_MINOR+"")+(nCharacters-1)+(RepeatAlphabet.SEPARATOR_MINOR+"")+maxOpenLength_unique);
		bw.newLine(); bw.close();
	}

}