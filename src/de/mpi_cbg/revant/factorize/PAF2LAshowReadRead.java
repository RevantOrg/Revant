package de.mpi_cbg.revant.factorize;

import java.util.*;
import java.io.*;
import de.mpi_cbg.revant.util.Math;

/**
 * Converts the reads-reads PAF file produced by minimap2, into the output of LAshow.
 *
 * Remark: reads with no alignment are not listed in the read lengths file printed in 
 * output and are never mentioned in the alignments file. Any downstream step should use 
 * the number of reads with an alignment as the effective number of reads.
 *
 * Remark: this is just a simple first attempt and it has not been optimized. 
 */
public class PAF2LAshowReadRead {
	/**
	 * @param args 4=keep only alignments of this length or longer.
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE_ALIGNMENTS = args[0];
		final String INPUT_FILE_LENGTHS = args[1];
		final String OUTPUT_FILE_ALIGNMENTS = args[2];
		final String OUTPUT_FILE_LENGTHS = args[3];
		final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[4]);
		final double MAX_ERROR_RATE = Double.parseDouble(args[5]);
		
		final String PAF_SEPARATOR = "\t";
		final int HASHMAP_CAPACITY = 1000000;  // Arbitrary
		final int BUFFER_SIZE_BYTES = 100000000;  // Arbitrary
		boolean orientation;
		int i, j;
		int tmpInt, idGenerator, length, minLength, lastAlignment;
		int readA, readB, startA, endA, startB, endB, diffs;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		HashMap<String,Integer> sequenceNames;
		String[] tokens;
		Pair[] pairs;
		AlignmentRow[] alignments;
		
		System.err.print("Translating sequence names into numbers... ");
		idGenerator=1;
		sequenceNames = new HashMap<String,Integer>(HASHMAP_CAPACITY);
		br = new BufferedReader(new FileReader(INPUT_FILE_ALIGNMENTS),BUFFER_SIZE_BYTES);
		str=br.readLine(); lastAlignment=-1;
		while (str!=null)  {
			tokens=str.split(PAF_SEPARATOR);
			if (!sequenceNames.containsKey(tokens[0])) sequenceNames.put(tokens[0],Integer.valueOf(idGenerator++));
			if (!sequenceNames.containsKey(tokens[5])) sequenceNames.put(tokens[5],Integer.valueOf(idGenerator++));
			lastAlignment++;
			str=br.readLine();
		}
		br.close(); 
		System.err.println(" DONE");
		
		System.err.print("Converting alignments... ");
		alignments = new AlignmentRow[(lastAlignment+1)<<1];
		br = new BufferedReader(new FileReader(INPUT_FILE_ALIGNMENTS),BUFFER_SIZE_BYTES);
		str=br.readLine(); lastAlignment=-1; minLength=Math.POSITIVE_INFINITY;
		while (str!=null)  {
			tokens=str.split(PAF_SEPARATOR);
			if (tokens[0].equalsIgnoreCase(tokens[5])) {
				// Removing self-alignments
				str=br.readLine();
				continue;
			}
			readA=sequenceNames.get(tokens[0]).intValue();
			readB=sequenceNames.get(tokens[5]).intValue();
			orientation=tokens[4].charAt(0)=='+';
			startA=Integer.parseInt(tokens[2]);
			endA=Integer.parseInt(tokens[3]);
			startB=Integer.parseInt(tokens[7]);
			endB=Integer.parseInt(tokens[8]);
			length=(endB-startB+endA-startA)>>1;
			diffs=Integer.parseInt(tokens[10])-Integer.parseInt(tokens[9]);
			if (length>=MIN_ALIGNMENT_LENGTH && ((double)(diffs<<1))/(endA-startA+1+endB-startB+1)<=MAX_ERROR_RATE) {
				minLength=Math.min(minLength,length);
				if (orientation) {
					alignments[++lastAlignment] = new AlignmentRow(readA,startA>=0?startA:0,endA,readB,startB>=0?startB:0,endB,true,diffs);
					alignments[++lastAlignment] = new AlignmentRow(readB,startB>=0?startB:0,endB,readA,startA>=0?startA:0,endA,true,diffs);
				}
				else {
					alignments[++lastAlignment] = new AlignmentRow(readA,startA>=0?startA:0,endA,readB,endB-1,startB-1>=0?startB-1:0,false,diffs);
					alignments[++lastAlignment] = new AlignmentRow(readB,startB>=0?startB:0,endB,readA,endA-1,startA-1>=0?startA-1:0,false,diffs);
				}
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
		
		System.err.print("Converting read lengths... ");
		pairs = new Pair[sequenceNames.size()];
		br = new BufferedReader(new FileReader(INPUT_FILE_LENGTHS),BUFFER_SIZE_BYTES);
		j=0; str=br.readLine();
		while (str!=null)  {
			tokens=str.split(PAF_SEPARATOR);
			if (sequenceNames.containsKey(tokens[0])) pairs[j++] = new Pair(sequenceNames.get(tokens[0]).intValue(),Integer.parseInt(tokens[1]));
			str=br.readLine();
		}
		br.close();
		System.err.println(" DONE");
		Arrays.sort(pairs,0,j);
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_LENGTHS),BUFFER_SIZE_BYTES);
		for (i=0; i<j; i++) bw.write(pairs[i].length+"\n");
		bw.close();
		
		System.err.println("Printing mapping oldID -> newID:");
		Set<Map.Entry<String,Integer>> set = sequenceNames.entrySet();
		Iterator<Map.Entry<String,Integer>> iterator = set.iterator();
		while (iterator.hasNext()) {
			Map.Entry<String,Integer> entry = iterator.next();
			System.out.println(entry.getKey()+PAF_SEPARATOR+entry.getValue().intValue());
		}
	}
	
	
	private static class Pair implements Comparable {
		public int read, length;
		
		public Pair(int r, int l) {
			this.read=r; this.length=l;
		}
		
		public int compareTo(Object other) {
			Pair otherPair = (Pair)other;
			if (read<otherPair.read) return -1;
			else if (read>otherPair.read) return 1;
			else return 0;
		}
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
	}

}