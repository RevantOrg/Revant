package de.mpi_cbg.revant.factorize;

import java.io.*;
import de.mpi_cbg.revant.util.Math;


public class Damar2Phred  {

	public static void main(String[] args) throws IOException {
		final String DAMAR_FILE = args[0];
		final String READ_LENGTHS_FILE = args[1];
		final String READ_IDS_FILE = args[2];
		final int N_READS = Integer.parseInt(args[3]);
		final String QUALITY_THRESHOLD_FILE = args[4];

		int i, j, p;
		int length, read, previousRead;
		String str;
		BufferedReader br;
		
		Reads.loadQualities_thresholds(QUALITY_THRESHOLD_FILE);
		Reads.nReads=N_READS;
		Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		br = new BufferedReader(new FileReader(DAMAR_FILE));
		str=br.readLine(); previousRead=Reads.readIDs[0]-1;
		while (str!=null) {
			p=str.indexOf(' ');
			read=Integer.parseInt(str.substring(0,p))-1;  // Read IDs are one-based
			for (i=previousRead+1; i<read; i++) {
				length=Math.ceil(Reads.getReadLength(i),Reads.QUALITY_SPACING);
				for (j=0; j<length; j++) System.out.print('Y');
				System.out.println();
			}
			p=str.indexOf(' ',p+1);
			System.out.println(str.substring(p+1)+"\"");  // Adding a good value at the end
			previousRead=read;
			str=br.readLine();
		}
		br.close();
		read=Reads.readAtIndex(Reads.nReads-1);
		for (i=previousRead+1; i<=read; i++) {
			length=Math.ceil(Reads.getReadLength(i),Reads.QUALITY_SPACING);
			for (j=0; j<length; j++) System.out.print('Y');
			System.out.println();
		}
	}

}


