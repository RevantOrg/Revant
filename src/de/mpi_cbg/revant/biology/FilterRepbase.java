package de.mpi_cbg.revant.biology;

import java.io.*;
import de.mpi_cbg.revant.util.IO;


public class FilterRepbase {
	
	/**
	 * Only sequences of $args[0]$ whose header contains $args[1]$ (case ignored) are 
	 * output to file $args[2]$.
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		final String KEYWORD = args[1].toLowerCase();
		final String OUTPUT_FILE = args[2];
		boolean selected;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE),IO.BUFFER_SIZE);
		br = new BufferedReader(new FileReader(INPUT_FILE),IO.BUFFER_SIZE);
		str=br.readLine();
		selected=false;
		while (str!=null) {
			if (str.length()==0) {
				str=br.readLine();
				continue;
			}
			if (str.charAt(0)=='>') {
				if (str.toLowerCase().indexOf(KEYWORD)>=0) {
					selected=true;
					bw.write(str+"\n");
				}
				else selected=false;
			}
			else if (selected) bw.write(str+"\n");
			str=br.readLine();
		}
		br.close(); bw.close();
	}

}