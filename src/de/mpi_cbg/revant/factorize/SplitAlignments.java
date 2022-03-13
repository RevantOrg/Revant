package de.mpi_cbg.revant.factorize;

import java.io.*;
import de.mpi_cbg.revant.util.IO;


/**
 * Splits an alignment file in approx. equal-sized parts, aligned with readA boundaries.
 */
public class SplitAlignments {

	public static void main(String[] args) throws IOException {
		final int N_ALIGNMENTS = Integer.parseInt(args[0]);
		final int N_PIECES = Integer.parseInt(args[1]);
		final String INPUT_FILE = args[2];
		final String OUTPUT_PREFIX = args[3];
		final String LAST_READA_FILE = args[4];  // Last readA (zero-based) in each piece
		
		int p, q;
		int fileID, nAlignmentsInFile, quantum, readA, previousReadA;
		String header1, header2, str;
		BufferedReader br;
		BufferedWriter outputFile = null;
		BufferedWriter lastReadsFile;
		
		quantum=N_ALIGNMENTS/N_PIECES;
		br = new BufferedReader(new FileReader(INPUT_FILE),IO.BUFFER_SIZE);
		header1=br.readLine();
		header2=br.readLine();  // Skipping the first two lines
		lastReadsFile = new BufferedWriter(new FileWriter(LAST_READA_FILE));
		previousReadA=-1; fileID=-1; nAlignmentsInFile=0;
		str=br.readLine();
		while (str!=null) {
			p=0;
			while (p<str.length() && str.charAt(p)==' ') p++;
			q=p+1;
			while (q<str.length() && str.charAt(q)!=' ') q++;
			readA=Alignments.intSubstring(str,p,q)-1;
			if (previousReadA==-1 || (readA!=previousReadA && nAlignmentsInFile>=quantum)) {
				if (outputFile!=null) {
					outputFile.close();
					System.out.println(nAlignmentsInFile+" alignments");
				}
				if (previousReadA!=-1) lastReadsFile.write(previousReadA+"\n");
				fileID++;
				outputFile = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+fileID+".txt"));
				outputFile.write(header1); outputFile.newLine();
				outputFile.write(header2); outputFile.newLine();
				nAlignmentsInFile=0;
			}
			outputFile.write(str); outputFile.newLine(); nAlignmentsInFile++;
			previousReadA=readA;
			str=br.readLine();
		}
		if (previousReadA!=-1) lastReadsFile.write(previousReadA+"\n");
		System.out.println(nAlignmentsInFile+" alignments");
		outputFile.close(); lastReadsFile.close();
		br.close();
	}

}