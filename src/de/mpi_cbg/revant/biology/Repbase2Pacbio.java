package de.mpi_cbg.revant.biology;

import java.io.*;
import de.mpi_cbg.revant.util.IO;


public class Repbase2Pacbio {
	
	/**
	 * Puts the headers of a RepBase file into PacBio format
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		final String OUTPUT_FILE = args[1];
		int i;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE),IO.BUFFER_SIZE);
		br = new BufferedReader(new FileReader(INPUT_FILE),IO.BUFFER_SIZE);
		str=br.readLine(); i=0;
		while (str!=null) {
			if (str.length()==0) bw.write("\n");
			else {
				if (str.charAt(0)=='>') {
					// Random string that happens to comply with the PacBio format
					bw.write(">L416/"+(i++)+"/0_16506 RQ=0.873 \n");
				}
				else bw.write(str+"\n");
			}
			str=br.readLine();
		}
		br.close(); bw.close();
	}

}