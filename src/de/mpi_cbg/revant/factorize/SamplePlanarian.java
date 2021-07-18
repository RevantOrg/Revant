package de.mpi_cbg.revant.factorize;

import java.util.*;
import java.io.*;


public class SamplePlanarian {
	
	/**
	 * Given a directory with a number of FASTA files, collects reads of total length
	 * approx. equal to the length of the planarian genome, using only reads with distinct
	 * ID up to the second slash.
	 */
	public static void main(String[] args) throws IOException {
		final long GENOME_LENGTH = 8000000000L;  // We want 10x coverage
		final String FASTA_SUFFIX = ".fasta";
		final String READS_DIR = args[0];
		final String OUTPUT_FILE = args[1];
		final int MIN_READ_LENGTH = Integer.parseInt(args[2]);
		int i, p;
		long currentSize;
		String str, id, header;
		BufferedReader br;
		BufferedWriter bw;
		HashSet<String> ids;
		StringBuilder sb = new StringBuilder();
		String[] files;
		
		ids = new HashSet<String>();
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE),10000);
		currentSize=0;
		files = new File(READS_DIR).list();
		for (i=0; i<files.length; i++) {
			if (files[i].length()<FASTA_SUFFIX.length() || !files[i].substring(files[i].length()-FASTA_SUFFIX.length()).equalsIgnoreCase(FASTA_SUFFIX)) continue;
			br = new BufferedReader(new FileReader(READS_DIR+"/"+files[i]),10000);
			str=br.readLine();
			do {
				// $str$ is always a header here
				header=str;
				p=str.indexOf('/'); p=str.indexOf('/',p+1); id=str.substring(1,p);
				sb.delete(0,sb.length());
				str=br.readLine();
				while (str!=null && str.charAt(0)!='>') {
					sb.append(str);
					str=br.readLine();
				}
				if (!ids.contains(id) && sb.length()>=MIN_READ_LENGTH) {
					ids.add(id);
					bw.write(header+"\n");
					bw.write(sb+"\n");
					currentSize+=sb.length();
					if (currentSize>=GENOME_LENGTH) break;
				}
			} while (str!=null);
			br.close();
			if (currentSize>=GENOME_LENGTH) break;
		}
		bw.close();
	}

}