package de.mpi_cbg.revant.graphics;

import java.awt.image.*;
import java.io.*;
import javax.imageio.*;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Reads;


public class PrintEdge {

	public static void main(String[] args) throws Exception {
		
		String inputFile = "/Users/ramseysnow/Desktop/histone/fourthExperiment/allReads/read50.png";
		int read = 49;
		int intervalFirst = 8630;
		int intervalLast = 10791;
		int intervalType = Constants.INTERVAL_PERIODIC;
		int xContext = 2000;
		int yContext = 400;
		String outputFile = "snapshot.png";
		
		Reads.nReads=7189;
		Reads.maxReadLength=54186;
		Reads.loadReadLengths("/Users/ramseysnow/Desktop/histone/fourthExperiment/readLengths.txt");
		
		int[] testArray = new int[3];
		BufferedImage snapshot = PrintReadAlignments.getIntervalContext(inputFile,read,intervalFirst,intervalLast,intervalType,xContext,yContext,testArray);
		if (snapshot!=null) ImageIO.write(snapshot,"png",new File(outputFile));
	}




	





}