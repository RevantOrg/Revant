package de.mpi_cbg.revant.apps;

import java.io.*;

/**
 * Splits a spacers file, a read IDs file and a read lengths file into blocks, based on a
 * file that specifies the last read ID of each block.
 */
public class SplitSpacers {
	/**
	 * @param args 2: use $null$ to avoid splitting spacers.
	 */
	public static void main(String[] args) throws IOException {
		final String LAST_READ_FILE = args[0];  // Of every block (last block included).
		final int N_BLOCKS = Integer.parseInt(args[1]);
		final String SPACERS_FILE = args[2];
		final String READ_IDS_FILE = args[3];
		final String READ_LENGTHS_FILE = args[4];
		final String OUTPUT_PREFIX = args[5];
		
		int i;
		String str1, str2;
		BufferedReader br1, br2;
		BufferedWriter bw1, bw2;
		int[] lastRead;

		lastRead = new int[N_BLOCKS];
		br1 = new BufferedReader(new FileReader(LAST_READ_FILE));
		for (i=0; i<N_BLOCKS; i++) lastRead[i]=Integer.parseInt(br1.readLine());
		br1.close();
		
		// Serializing spacers
		if (!SPACERS_FILE.equalsIgnoreCase("null")) {
			i=0;
			bw1 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"0.txt"));
			br1 = new BufferedReader(new FileReader(SPACERS_FILE));
			str1=br1.readLine();
			while (str1!=null) {
				if (Integer.parseInt(str1.substring(0,str1.indexOf(RepeatAlphabet.SEPARATOR_MINOR)))>lastRead[i]) {
					bw1.close();
					i++;
					bw1 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+i+".txt"));
				}
				bw1.write(str1); bw1.newLine();
				str1=br1.readLine();
			}
			bw1.close(); br1.close();
		}
		
		// Serializing read IDs and lengths
		i=0;
		bw1 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"ids-0.txt"));
		bw2 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"lengths-0.txt"));
		br1 = new BufferedReader(new FileReader(READ_IDS_FILE));
		br2 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null) {
			if (Integer.parseInt(str1)>lastRead[i]) {
				bw1.close(); bw2.close();
				i++;
				bw1 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"ids-"+i+".txt"));
				bw2 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"lengths-"+i+".txt"));
			}
			bw1.write(str1); bw1.newLine();
			bw2.write(str2); bw2.newLine();
			str1=br1.readLine(); str2=br2.readLine();
		}
		bw1.close(); bw2.close(); br1.close(); br2.close();
	}

}