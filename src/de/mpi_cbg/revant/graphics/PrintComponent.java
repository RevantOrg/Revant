package de.mpi_cbg.revant.graphics;

import java.awt.image.*;
import java.io.*;
import javax.imageio.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.intervalgraph.IntervalGraph;


public class PrintComponent {

	public static void main(String[] args) throws Exception {
		final String GRAPH_FILE = args[0];
		final String INPUT_DIR = args[1];
		final String OUTPUT_DIR = args[2];
		final int X_CONTEXT = 2000;
		final int Y_CONTEXT = 400;
		int i;
		int read;
		int[] tmpArray = new int[3];
		BufferedImage snapshot;
		
		IO.initialize();
		System.err.print("Loading interval graph "+GRAPH_FILE+"... ");
		IntervalGraph.deserialize(GRAPH_FILE,true,true);
		Reads.nReads=7189;
		Reads.maxReadLength=54186;
		Reads.loadReadLengths("/Users/ramseysnow/Desktop/histone/fourthExperiment/readLengths.txt");
		System.err.println("done");
		
		for (i=0; i<IntervalGraph.nNodes; i++) {
			read=IntervalGraph.nodesArray[i].read;			
			snapshot=PrintReadAlignments.getIntervalContext(INPUT_DIR+"/read"+(read+1)+".png",read,IntervalGraph.nodesArray[i].start,IntervalGraph.nodesArray[i].end,IntervalGraph.nodesArray[i].type,X_CONTEXT,Y_CONTEXT,tmpArray);
			if (snapshot!=null) ImageIO.write(snapshot,"png",new File(OUTPUT_DIR+"/node"+i+".png"));
		}
	}

}