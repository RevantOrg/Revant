package de.mpi_cbg.revant.factorize;

import java.util.*;
import java.io.*;
import de.mpi_cbg.revant.util.Math;

/**
 * Converts the reads-repeats PAF file produced by minimap2, into the output of LAshow.
 * The program assumes that $PAF2LAshowReadRead$ has already been executed, and it 
 * discards any read that is not listed in the conversion file produced by that program.
 *
 * Remark: this is just a simple first attempt and it has not been optimized. 
 */
public class PAF2LAshowReadRepeat {	
	/**
	 * @param args 
	 * 0: read-repeat alignment file; the program assumes that repeat IDs are $X/Y/Z$
	 * where $Y$ is the zero-based integer ID of the repeat;
	 * 1: conversion file generated by $PAF2LAshowReadRead$; every row is $X \t Y$, where 
	 * $X$ is the original string ID of a read, and $Y$ is its integer ID (zero-based).
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE_ALIGNMENTS = args[0];
		final int N_ALIGNMENTS = Integer.parseInt(args[1]);
		final String INPUT_FILE_CONVERSION = args[2];
		final String OUTPUT_FILE_ALIGNMENTS = args[3];
		final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[4]);
		final double MAX_ERROR_RATE = Double.parseDouble(args[5]);
		
		final String PAF_SEPARATOR = "\t";
		final String REPEAT_ID_SEPARATOR = "/";
		final int HASHMAP_CAPACITY = 1000000;  // Arbitrary
		final int BUFFER_SIZE_BYTES = 100000000;  // Arbitrary
		boolean orientation;
		int i, j, p, q;
		int length, minLength, lastAlignment;
		int readA, readB, startA, endA, startB, endB, diffs;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		Integer tmpInt;
		HashMap<String,Integer> readNames;
		String[] tokens;
		AlignmentRow[] alignments;
		
		System.err.println("Loading read name conversion file... ");
		readNames = new HashMap<String,Integer>(HASHMAP_CAPACITY);
		br = new BufferedReader(new FileReader(INPUT_FILE_CONVERSION),BUFFER_SIZE_BYTES);
		str=br.readLine();
		while (str!=null)  {
			tokens=str.split(PAF_SEPARATOR);
			readNames.put(tokens[0],Integer.valueOf(Integer.parseInt(tokens[1])));
			str=br.readLine();
		}
		br.close();
		
		System.err.println("Converting alignments... ");
		alignments = new AlignmentRow[N_ALIGNMENTS];
		br = new BufferedReader(new FileReader(INPUT_FILE_ALIGNMENTS),BUFFER_SIZE_BYTES);
		str=br.readLine(); lastAlignment=-1; minLength=Math.POSITIVE_INFINITY;
		while (str!=null)  {
			tokens=str.split(PAF_SEPARATOR);
			tmpInt=readNames.get(tokens[0]);
			if (tmpInt==null) {
				str=br.readLine();
				continue;
			}
			readA=tmpInt.intValue();
			p=tokens[5].indexOf(REPEAT_ID_SEPARATOR);
			q=tokens[5].indexOf(REPEAT_ID_SEPARATOR,p+1);			
			readB=Integer.parseInt(tokens[5].substring(p+1,q))+1;  // The original value is zero-based
			orientation=tokens[4].charAt(0)=='+';
			startA=Integer.parseInt(tokens[2]);
			endA=Integer.parseInt(tokens[3]);
			startB=Integer.parseInt(tokens[7]);
			endB=Integer.parseInt(tokens[8]);
			length=(endB-startB+endA-startA)>>1;
			diffs=Integer.parseInt(tokens[10])-Integer.parseInt(tokens[9]);
			if (length>=MIN_ALIGNMENT_LENGTH && ((double)(diffs<<1))/(endA-startA+1+endB-startB+1)<=MAX_ERROR_RATE) {
				minLength=Math.min(minLength,length);
				if (orientation) alignments[++lastAlignment] = new AlignmentRow(readA,startA>=0?startA:0,endA,readB,startB>=0?startB:0,endB,true,diffs);
				else alignments[++lastAlignment] = new AlignmentRow(readA,startA>=0?startA:0,endA,readB,endB-1,startB-1>=0?startB-1:0,false,diffs);
			}
			str=br.readLine();
		}
		br.close();
		
		// Sorting and removing duplicates (which might occur in the original PAF file).
		Arrays.sort(alignments,0,lastAlignment+1);
		j=0;
		for (i=1; i<=lastAlignment; i++) {
			if (!alignments[i].equals(alignments[j])) {
				j++;
				if (j!=i) alignments[j].copyFrom(alignments[i]);
			}
		}
		lastAlignment=j;
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_ALIGNMENTS),BUFFER_SIZE_BYTES);
		bw.write("\n"); bw.write("\n");  // Fake header
		for (i=0; i<=lastAlignment; i++) bw.write(alignments[i].readA+"  "+alignments[i].readB+"  "+(alignments[i].orientation?'n':'c')+"  ["+alignments[i].startA+".. "+alignments[i].endA+"] x ["+alignments[i].startB+".. "+alignments[i].endB+"] ( "+alignments[i].diffs+" diffs)\n");
		bw.close();
		System.err.println(" DONE  minAlignmentLength="+minLength);		
	}
	
	
	private static class AlignmentRow implements Comparable {
		public int readA, startA, endA, readB, startB, endB, diffs;
		public boolean orientation;
		
		public AlignmentRow(int ra, int sa, int ea, int rb, int sb, int eb, boolean o, int d) {
			this.readA=ra; this.startA=sa; this.endA=ea;
			this.readB=rb; this.startB=sb; this.endB=eb;
			this.orientation=o; this.diffs=d;
		}
		
		public int compareTo(Object other) {
			AlignmentRow otherAlignment = (AlignmentRow)other;
			if (readA<otherAlignment.readA) return -1;
			else if (readA>otherAlignment.readA) return 1;
			if (readB<otherAlignment.readB) return -1;
			else if (readB>otherAlignment.readB) return 1;
			if (orientation && !otherAlignment.orientation) return -1;
			else if (!orientation && otherAlignment.orientation) return 1;
			if (startA<otherAlignment.startA) return -1;
			else if (startA>otherAlignment.startA) return 1;
			if (startB<otherAlignment.startB) return -1;
			else if (startB>otherAlignment.startB) return 1;
			if (endA<otherAlignment.endA) return -1;
			else if (endA>otherAlignment.endA) return 1;
			if (endB<otherAlignment.endB) return -1;
			else if (endB>otherAlignment.endB) return 1;
			return 0;
		}
		
		public boolean equals(Object other) {
			AlignmentRow otherAlignment = (AlignmentRow)other;
			return readA==otherAlignment.readA && startA==otherAlignment.startA && endA==otherAlignment.endA && 
				   readB==otherAlignment.readB && startB==otherAlignment.startB && endB==otherAlignment.endB && 
				   orientation==otherAlignment.orientation;
		}
		
		public void copyFrom(AlignmentRow otherAlignment) {
			readA=otherAlignment.readA; startA=otherAlignment.startA; endA=otherAlignment.endA;
			readB=otherAlignment.readB; startB=otherAlignment.startB; endB=otherAlignment.endB;
			orientation=otherAlignment.orientation; 
			// A PAF file might contain identical alignments with different diffs.
			diffs=Math.min(diffs,otherAlignment.diffs);
		}
		
		public String toString() {
			return readA+"["+startA+".."+endA+"] x "+readB+"["+startB+".."+endB+"] "+orientation+","+diffs;
		}
	}

}