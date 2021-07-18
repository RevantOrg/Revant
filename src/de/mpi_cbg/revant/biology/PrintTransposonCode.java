package de.mpi_cbg.revant.biology;

import java.io.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Alignments;

/**
 * 
 */
public class PrintTransposonCode {
	
	/**
	 * 
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(INPUT_FILE));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			System.out.println("if (read=="+Alignments.readB+" && Intervals.jaccardSimilarity(start,end,"+Alignments.startB+","+Alignments.endB+")>=JACCARD_THRESHOLD) System.err.println(\"Transposon "+(Alignments.readA-1)+" intersects an interval of type=\"+IntervalGraph.nodesArray[x].type+\" in read \"+read+\" -> component \"+IntervalGraph.nodesArray[x].components[0]);");
			str=br.readLine();
		}
		br.close();
	}
	
}
