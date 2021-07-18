package de.mpi_cbg.revant.biology;

import java.io.*;


public class PrintHistograms {

	/**
	 * java PrintHistograms copia ./edgeHistograms.txt 1> histogram-copia.txt 2> histogram-copia-clusters.txt
	 */
	public static void main(String[] args) throws IOException {
		final String REPBASE_LABEL = args[0].trim();
		final String HISTOGRAM_FILE = args[1];
		final String SEPARATOR = " -> ";
		int i;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(HISTOGRAM_FILE));
		str=br.readLine();
		while (str!=null) {
			i=str.indexOf(SEPARATOR);
			if (i<0) {
				str=br.readLine();
				continue;
			}
			if (str.substring(i+SEPARATOR.length()).trim().equalsIgnoreCase(REPBASE_LABEL)) {
				System.err.println(str.substring(0,i).trim());
				str=br.readLine();
				System.out.println(str);
			}
			str=br.readLine();
		}
	}

}