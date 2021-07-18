package de.mpi_cbg.revant.intervalgraph;

import java.io.*;

/**
 * 
 */
public class PrintAllAssemblyPaths {
	
	/**
	 * 
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_GRAPH_DIR = args[0];
		final String INPUT_GRAPH_ID = args[1];
		final String INPUT_KERNEL_ID = args[2];
		
		System.err.print("Loading bidirected graph... ");
		BidirectedGraph.deserialize(INPUT_GRAPH_DIR+"/"+INPUT_GRAPH_ID+"-kernel"+INPUT_KERNEL_ID+".bdgraph");
		System.err.println("done.");
		BidirectedGraph.print();
		
		System.err.print("Printing all paths... ");
		BidirectedGraph.printAllPaths(INPUT_GRAPH_DIR+"/"+INPUT_GRAPH_ID+"-kernel"+INPUT_KERNEL_ID+"-labels.txt",INPUT_GRAPH_DIR+"/"+INPUT_GRAPH_ID+"-kernel"+INPUT_KERNEL_ID+"-paths.txt");
		System.err.println("done.");
	}
	
	
}