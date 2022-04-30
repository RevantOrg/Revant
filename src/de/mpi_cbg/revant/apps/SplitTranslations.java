package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Splits a pair of translation files (characters and boundaries) in chunks specified in a
 * given last read file.
 */
public class SplitTranslations {
	
	public static void main(String[] args) throws IOException {
		final String READ_IDS_FILE = args[0];  // 0-based
		final String CHARACTERS_FILE = args[1];
		final String BOUNDARIES_FILE = args[2];
		final String LAST_READ_FILE = args[3];  // 0-based, last block included.
		final String CHARACTERS_FILE_PREFIX = args[4];
		final String BOUNDARIES_FILE_PREFIX = args[5];
		
		int i;
		int read, nBlocks;
		String str1, str2, str3;
		int[] lastRead;
		BufferedReader br1, br2, br3;
		BufferedWriter bw1, bw2;
		
		// Loading last read file
		nBlocks=0;
		br1 = new BufferedReader(new FileReader(LAST_READ_FILE));
		str1=br1.readLine();
		while (str1!=null) {
			nBlocks++;
			str1=br1.readLine();
		}
		br1.close();
		lastRead = new int[nBlocks];
		br1 = new BufferedReader(new FileReader(LAST_READ_FILE));
		for (i=0; i<nBlocks; i++) lastRead[i]=Integer.parseInt(br1.readLine());
		br1.close();
		
		// Splitting
		bw1 = new BufferedWriter(new FileWriter(CHARACTERS_FILE_PREFIX+"0.txt"));
		bw2 = new BufferedWriter(new FileWriter(BOUNDARIES_FILE_PREFIX+"0.txt"));
		br1 = new BufferedReader(new FileReader(READ_IDS_FILE));
		br2 = new BufferedReader(new FileReader(CHARACTERS_FILE));
		br3 = new BufferedReader(new FileReader(BOUNDARIES_FILE));
		i=0; str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
		while (str1!=null) {
			read=Integer.parseInt(str1);
			if (read>lastRead[i]) {
				bw1.close(); bw2.close(); i++;
				bw1 = new BufferedWriter(new FileWriter(CHARACTERS_FILE_PREFIX+i+".txt"));
				bw2 = new BufferedWriter(new FileWriter(BOUNDARIES_FILE_PREFIX+i+".txt"));
			}
			bw1.write(str2); bw1.newLine(); bw2.write(str3); bw2.newLine();
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
		}
		bw1.close(); bw2.close(); br1.close(); br2.close(); br3.close();
	}

}