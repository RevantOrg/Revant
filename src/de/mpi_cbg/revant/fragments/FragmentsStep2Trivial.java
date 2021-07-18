package de.mpi_cbg.revant.fragments;

import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Factorize;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Alignments;


/**
 * Simplified version of $FragmentsStep2$ that keeps the existing reference instead of 
 * creating a new one. The program copies the input fragments to the output, filtering 
 * them using the connection file built by $FragmentsStep1$ and the alignments of 
 * fragments to the reference.
 *
 * Remark: the program does not copy the existing reference sequence to the output; this
 * must be done by the caller.
 */
public class FragmentsStep2Trivial {
	
	/**
	 * @param args 
	 * 0: number of fragments (including those in the old reference sequence);
	 * 1: directory containing alignments of fragments (to fragments and reference);
	 * 4: output dir containing the new fragments, i.e. only the fragments that align 
	 * to the reference, with the boundaries specified in the connection file.
	 */
	public static void main(String[] args) throws IOException {
		Reads.nReads=Integer.parseInt(args[0]);
		final String ALIGNMENTS_DIR = args[1];
		final String BASIN_DESCRIPTOR_ID = args[2];
		final String FRAGMENTS_STRINGS_FILE = args[3];
		final String NEW_FRAGMENTS_DIR = args[4];
		
		final String READ_LENGTHS_FILE = ALIGNMENTS_DIR+"/reads-lengths.txt";
		final String FF_ALIGNMENTS_PATH = ALIGNMENTS_DIR+"/LAshow.txt";
		final String FR_ALIGNMENTS_PATH = ALIGNMENTS_DIR+"/"+IO.FRAGMENTS_REFERENCE_FILE+".txt";
		final String CONNECTION_FILE = ALIGNMENTS_DIR+"/connection.txt";
		final String NEW_FRAGMENTS_FILE = NEW_FRAGMENTS_DIR+"/fragments-"+BASIN_DESCRIPTOR_ID+".txt";
		final String NEW_FRAGMENTS_LENGTHS = NEW_FRAGMENTS_DIR+"/fragments-"+BASIN_DESCRIPTOR_ID+"-lengths.txt";
		
		int i;
		int length, readA, header;
		String str1, str2;
		BufferedReader br1, br2;
		BufferedWriter bw1, bw2;
		boolean[] aligns, inReference;
		int[][] newCoordinates;
		
		Reads.loadReadIDs(0,Reads.nReads-1);
		Reads.loadReadLengths(READ_LENGTHS_FILE);
		
		// Loading $aligns$ and $inReference$.
		aligns = new boolean[Reads.nReads];
		Math.set(aligns,Reads.nReads-1,false);
		inReference = new boolean[Reads.nReads];
		Math.set(inReference,Reads.nReads-1,false);
		br1 = new BufferedReader(new FileReader(FR_ALIGNMENTS_PATH));
		str1=br1.readLine(); str1=br1.readLine();  // Skipping header
		str1=br1.readLine();
		while (str1!=null) {
			Alignments.readAlignmentFile(str1);
			readA=Alignments.readA-1;
			aligns[readA]=true;
			if (Alignments.diffs<=FragmentsStep1.TOO_FEW_DIFFS) inReference[readA]=true;
			str1=br1.readLine();
		}
		br1.close();
		
		// Loading new fragment coordinates
		newCoordinates = new int[Reads.nReads][2];
		Math.set(newCoordinates,-1);
		br1 = new BufferedReader(new FileReader(FF_ALIGNMENTS_PATH));
		br2 = new BufferedReader(new FileReader(CONNECTION_FILE));
		str1=br1.readLine(); str1=br1.readLine();  // Skipping header
		str1=br1.readLine(); 
		str2=br2.readLine();
		while (str1!=null) {
			Alignments.readAlignmentFile(str1);
			if (str2.length()!=0) {
				readA=Alignments.readA-1;
				Factorize.readIntervalConnectionFile(str2);
				newCoordinates[readA][0]=Factorize.start;
				newCoordinates[readA][1]=Factorize.end;
			}
			str1=br1.readLine(); str2=br2.readLine();
		}
		br1.close(); br2.close();
		
		// Printing new fragments
		br1 = new BufferedReader(new FileReader(FRAGMENTS_STRINGS_FILE));
		bw1 = new BufferedWriter(new FileWriter(NEW_FRAGMENTS_FILE));
		bw2 = new BufferedWriter(new FileWriter(NEW_FRAGMENTS_LENGTHS));
		header=Math.POSITIVE_INFINITY>>1;
		for (i=0; i<Reads.nReads; i++) {
			if (aligns[i] && !inReference[i] && newCoordinates[i][0]!=-1) {
				length=newCoordinates[i][1]-newCoordinates[i][0]+1;
				IO.writeFakeHeader(header++/*Unlikely to collide with any meaningful header*/,length,null,bw1);
				str1=br1.readLine(); str1=br1.readLine();
				bw1.write(str1+"\n");
				bw2.write(length+"\n");
			}
			else { str1=br1.readLine(); str1=br1.readLine(); }
		}
		br1.close(); bw1.close(); bw2.close();
	}
	

}