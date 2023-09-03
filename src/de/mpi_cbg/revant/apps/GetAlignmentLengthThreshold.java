package de.mpi_cbg.revant.apps;

import java.io.*;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;

/**
 * Prints to STDOUT the max length of an alignment of a given type. This is used for later
 * enumerating only k-mers that are short enough to be contained in a suffix-prefix
 * alignment.
 *
 * Remark: one could return instead the min L such that, say, 90% of all alignments have
 * length <=L. However, using such a threshold to avoid considering unique k-mers does
 * make the assembly graph less connected in practice.
 *
 * Remark: k-mer enumeration should stop at max read-read alignment length, instead of at
 * a specific value of k. We currently cap it at a max k just for speed.
 */
public class GetAlignmentLengthThreshold {
	/**
	 * @param args 
	 * 1: alignment type: 0=suffix-prefix overlap; 1=local substring; 2=full containment
     * or full identity;
     * 4: histogram of all alignment lengths.
	 */
	public static void main(String[] args) throws IOException {
		final String ALIGNMENTS_FILE = args[0];
        final int TYPE = Integer.parseInt(args[1]);
        final int AVG_READ_LENGTH = Integer.parseInt(args[2]);
        final String OUTPUT_HISTOGRAM = args[3];
        
        final int QUANTUM = IO.quantum;
        final int N_CELLS = AVG_READ_LENGTH/QUANTUM;
        final int IDENTITY_THRESHOLD = QUANTUM;
        
		int i;
        int length, lengthA, lengthB, max;
        long nAlignments;
		String str;
        BufferedReader br;
        BufferedWriter bw;
        long[] histogram;
        
        histogram = new long[N_CELLS];
		br = new BufferedReader(new FileReader(ALIGNMENTS_FILE));
		str=br.readLine(); str=br.readLine();  // Skipping header
        str=br.readLine();
        max=0;
		while (str!=null)  {
			Alignments.readAlignmentFile(str);
            if (Alignments.readAlignmentFile_getType(IDENTITY_THRESHOLD,str)==TYPE) {
                lengthA=Alignments.endA-Alignments.startA+1;
                lengthB=Alignments.endB-Alignments.startB+1;
                length=Math.max(lengthA,lengthB);
                max=Math.max(max,length);
                histogram[Math.min(length/QUANTUM,N_CELLS-1)]++;
            }
            str=br.readLine();
        }
        br.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_HISTOGRAM));
        for (i=0; i<N_CELLS; i++) bw.write((i*QUANTUM)+","+histogram[i]+"\n");
        bw.write(max+",-1\n");
        bw.close();
        System.out.println(max+"");
	}

}