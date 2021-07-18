package de.mpi_cbg.revant.factorize;

import java.io.*;


public class SplitTracks {

	/**
	 *  
	 */
	public static void main(String[] args) throws IOException {
		int i;
		int readID, fileID, nPieces;
		String str, trackFilePrefix, lastReadsFile, trackFile;
		BufferedReader br;
		BufferedWriter outputFile = null;
		int[] lastReads;   // Last readID (zero-based) in each fragment.

		trackFilePrefix=args[0];
		nPieces=Integer.parseInt(args[1]);
		lastReadsFile=args[2];
		
		lastReads = new int[nPieces-1];
		br = new BufferedReader(new FileReader(lastReadsFile));
		for (i=0; i<nPieces-1; i++) lastReads[i]=Integer.parseInt(br.readLine());
		br.close();
		
		br = new BufferedReader(new FileReader(trackFilePrefix+".txt"));
		fileID=0;
		outputFile = new BufferedWriter(new FileWriter(trackFilePrefix+"-"+fileID+".txt"));
		br.readLine(); br.readLine(); br.readLine(); br.readLine();  // Skipping header
		readID=0; str=br.readLine(); 
		while (str!=null) {
			i=str.indexOf(" "); 
			// Removing line prefix
			if (i<0) outputFile.write("\n");
			else {
				i=str.indexOf(" ",i+1);
				if (i<0) outputFile.write("\n");
				else outputFile.write(str.substring(i+1)+"\n");
			}
			if (fileID<nPieces-1 && readID==lastReads[fileID]) {
				outputFile.close();
				fileID++;
				outputFile = new BufferedWriter(new FileWriter(trackFilePrefix+"-"+fileID+".txt"));
			}
			str=br.readLine(); readID++;
		}
		outputFile.close();
		br.close();
	}

}