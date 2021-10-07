package de.mpi_cbg.revant.biology;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;


/**
 * Basic implementation of a string graph. Every read is mapped to two nodes that 
 * represent its endpoints. Only suffix-prefix alignments are used for building edges.
 */
public class StringGraph {

	public static void main(String[] args) throws IOException {
		final String ALIGNMENTS_FILE = args[0];
		final int N_READS = Integer.parseInt(args[1]);
		final double MAX_ERROR = Double.parseDouble(args[2]);
		final int DISTANCE_THRESHOLD = Integer.parseInt(args[3]);
		final int READ_LENGTH = Integer.parseInt(args[4]);
		final String DOT_FILE = args[5];
		
		int i, j, k;
		int from, to;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		int[] lastNeighbor;
		int[][] neighbors;
		
		// Allocating memory
		Reads.readLengths = new int[N_READS];
		Math.set(Reads.readLengths,N_READS-1,READ_LENGTH);
		neighbors = new int[(N_READS)<<1][2];
		lastNeighbor = new int[(N_READS)<<1];
		for (i=0; i<N_READS; i++) {
			neighbors[i<<1][0]=(i<<1)+1;
			lastNeighbor[i<<1]=0;
			neighbors[(i<<1)+1][0]=i<<1;
			lastNeighbor[(i<<1)+1]=0;
		}
		
		// Building the graph
		br = new BufferedReader(new FileReader(ALIGNMENTS_FILE));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			if ( (2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2)>MAX_ERROR ||
				 (Alignments.startA>DISTANCE_THRESHOLD && Alignments.endA<READ_LENGTH-DISTANCE_THRESHOLD) ||
				 (Alignments.startB>DISTANCE_THRESHOLD && Alignments.endB<READ_LENGTH-DISTANCE_THRESHOLD) ||
				 (Alignments.startA<=DISTANCE_THRESHOLD && Alignments.endA>=READ_LENGTH-DISTANCE_THRESHOLD) ||
				 (Alignments.startB<=DISTANCE_THRESHOLD && Alignments.endB>=READ_LENGTH-DISTANCE_THRESHOLD)
			   ) {
				str=br.readLine();
				continue;
			}
			from=-1; to=-1;
			if (Alignments.startA<=DISTANCE_THRESHOLD) from=(Alignments.readA-1)<<1;
			else if (Alignments.endA>=READ_LENGTH-DISTANCE_THRESHOLD) from=((Alignments.readA-1)<<1)+1;
			if (Alignments.startB<=DISTANCE_THRESHOLD) to=(Alignments.readB-1)<<1;
			else if (Alignments.endB>=READ_LENGTH-DISTANCE_THRESHOLD) to=((Alignments.readB-1)<<1)+1;
			lastNeighbor[from]++;
			if (lastNeighbor[from]==neighbors[from].length) {
				int[] newArray = new int[neighbors[from].length<<1];
				System.arraycopy(neighbors[from],0,newArray,0,neighbors[from].length);
				neighbors[from]=newArray;
			}
			neighbors[from][lastNeighbor[from]]=to;
			lastNeighbor[to]++;
			if (lastNeighbor[to]==neighbors[to].length) {
				int[] newArray = new int[neighbors[to].length<<1];
				System.arraycopy(neighbors[to],0,newArray,0,neighbors[to].length);
				neighbors[to]=newArray;
			}
			neighbors[to][lastNeighbor[to]]=from;
			str=br.readLine();
		}
		br.close();
		
		// Removing duplicated edges
		for (i=0; i<(N_READS)<<1; i++) {
			if (lastNeighbor[i]<=0) continue;
			Arrays.sort(neighbors[i],0,lastNeighbor[i]+1);
			k=0;
			for (j=1; j<=lastNeighbor[i]; j++) {
				if (neighbors[i][j]!=neighbors[i][k]) neighbors[i][++k]=neighbors[i][j];
			}
			lastNeighbor[i]=k;
		}
		
		// Printing the graph
		bw = new BufferedWriter(new FileWriter(DOT_FILE));
		bw.write("graph G {");
		for (i=0; i<(N_READS)<<1; i++) {
			for (j=0; j<lastNeighbor[i]; j++) {
				if (neighbors[i][j]>i) bw.write(i+" -- "+neighbors[i][j]+";\n");
			}
		}
		bw.write("}");
		bw.close();
	}

}