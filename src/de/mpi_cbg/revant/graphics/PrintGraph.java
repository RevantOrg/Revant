package de.mpi_cbg.revant.graphics;

import java.io.*;
import java.awt.image.*;
import java.awt.*;
import javax.imageio.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.intervalgraph.IntervalGraph;


/**
 * 
 */
public class PrintGraph {
	
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		final String OUTPUT_FILE = args[1];
		final int COVERAGE = Integer.parseInt(args[2]);
		final int QUANTUM = 10;
		int i, j;
		int maxDegree;
		int[] originalIDs, histogram;
		int[] degreeHistogram = new int[70];
		int[] clusteringCoefficientHistogram = new int[70];
		int[] tmpArray1, tmpArray2;
		
		IO.initialize();
		IO.coverage=COVERAGE;
		System.err.print("Loading interval graph "+INPUT_FILE+"... ");
		originalIDs=IntervalGraph.deserialize(INPUT_FILE,true,true);
		maxDegree=IntervalGraph.turnOffEdges(-2,false);  // Printing all edges
		//IntervalGraph.turnOnEdges(); maxDegree=0;
		System.err.println("done. Max degree: "+maxDegree);
		IntervalGraph.toDot(OUTPUT_FILE,null,null,Math.POSITIVE_INFINITY);
		
		
		
		
/*		
for (i=0; i<IntervalGraph.nNodes; i++) {
	if (IntervalGraph.nodesArray[i].read==6083 && IntervalGraph.nodesArray[i].start==26722 && IntervalGraph.nodesArray[i].end==28277) {
		System.err.println("neighbors of node 6083[26722..28277]:");
		for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
			System.err.println("EDGE: "+IntervalGraph.neighbors[i][j]);
			System.err.println("NODE: "+IntervalGraph.nodesArray[IntervalGraph.neighbors[i][j].getTo(i)]);
		}
	}
}
*/	
	
	
	
	
		
		
		// Stats, for temporary use only...
		int nContainment = 0;
		int nOverlap = 0;
		int nInsertion = 0;
		int nSharedSubstring = 0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
				if (IntervalGraph.neighbors[i][j].containment!=-1) nContainment++;
				if (IntervalGraph.neighbors[i][j].overlap!=-1) nOverlap++;
				if (IntervalGraph.neighbors[i][j].insertion!=-1) nInsertion++;
				if (IntervalGraph.neighbors[i][j].sharedSubstring!=-1) nSharedSubstring++;
			}
		}
		System.err.println("containment="+nContainment+" overlap="+nOverlap+" insertion="+nInsertion+" sharedSubstring="+nSharedSubstring);
		
		
		// Drawing edges
/*		Reads.nReads=7189;
		Reads.maxReadLength=54186;
		histogram = new int[Reads.maxReadLength/QUANTUM+1];
		Reads.loadReadLengths("/Users/ramseysnow/Desktop/histone/fourthExperiment/readLengths.txt");
		IntervalGraph.getDegreeHistogram(degreeHistogram,maxDegree,true);
		tmpArray1 = new int[maxDegree];
		tmpArray2 = new int[maxDegree];
		IntervalGraph.getClusteringCoefficientHistogram(clusteringCoefficientHistogram,true,tmpArray1,tmpArray2);
		PrintReadAlignmentsGene.drawEdges(6206,Factorize.INTERVAL_ALIGNMENT,6354,8587,"/Users/ramseysnow/Desktop/histone/fourthExperiment/allReads","edges-x","png",QUANTUM,histogram,degreeHistogram,clusteringCoefficientHistogram,tmpArray1,tmpArray2,false);
*/



	}
	
	
}