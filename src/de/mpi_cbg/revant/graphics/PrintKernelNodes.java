package de.mpi_cbg.revant.graphics;

import java.awt.image.*;
import java.io.*;
import javax.imageio.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Reads;


public class PrintKernelNodes {
	
	private static final String[] SUFFIXES = new String[] {"-"+IO.LABEL_COORDINATES_LABEL+".txt","-"+IO.LABEL_COORDINATES_LABEL+"-"+IO.CYCLIC_LABEL+".txt","-"+IO.ONE_NODE_ONLY_LABEL+"-"+IO.LABEL_COORDINATES_LABEL+".txt","-"+IO.TOO_MANY_PATHS_LABEL+"-"+IO.LABEL_COORDINATES_LABEL+".txt"};
	

	/**
	 * java PrintKernelNodes /Users/ramseysnow/Desktop/histone/fourthExperiment/step1/10-clusters 0	33 /Users/ramseysnow/Desktop/histone/fourthExperiment/allReads  /Users/ramseysnow/Desktop/histone/fourthExperiment/step1/10-clusters/allKernels
	 */
	public static void main(String[] args) throws Exception {
		final String KERNELS_DIR = args[0];
		final int GRAPH_ID = Integer.parseInt(args[1]);
		final int LAST_KERNEL = Integer.parseInt(args[2]);
		final String PICTURES_DIR = args[3];
		final String OUTPUT_DIR = args[4];
		final int X_CONTEXT = 10;
		final int Y_CONTEXT = -1;
		int i, j, p, q;
		int row, read, start, end;
		String str, fileName;
		FileReader fr;
		BufferedReader br;
		BufferedImage snapshot;
		int[] tmpArray = new int[3];
			
		Reads.nReads=7189;
		Reads.maxReadLength=54186;
		Reads.loadReadLengths("/Users/ramseysnow/Desktop/histone/fourthExperiment/readLengths.txt");
		for (i=0; i<=LAST_KERNEL; i++) {
			fileName=GRAPH_ID+"-kernel"+i;
			fr=null;
			for (j=0; j<SUFFIXES.length; j++) {
				try { 
					fr = new FileReader(KERNELS_DIR+"/"+fileName+SUFFIXES[j]);
					break;
				}
				catch (FileNotFoundException e) { }
			}
			if (fr==null) {
				System.err.println("WARNING: kernel <"+fileName+"> not found.");
				continue;
			}
			br = new BufferedReader(fr);
			str=br.readLine(); row=0;
			while (str!=null) {
				p=str.indexOf(",");
				read=Integer.parseInt(str.substring(0,p));
				p++; q=str.indexOf(",",p);
				start=Integer.parseInt(str.substring(p,q));
				p=q+1;
				end=Integer.parseInt(str.substring(p));
				snapshot=PrintReadAlignments.getIntervalContext(PICTURES_DIR+"/read"+(read+1)+".png",read,start,end,-1,X_CONTEXT,Y_CONTEXT,tmpArray);
				if (snapshot!=null) {
					ImageIO.write(snapshot,"png",new File(OUTPUT_DIR+"/"+fileName+"-"+row+".png"));
					System.out.println("Printing "+OUTPUT_DIR+"/"+fileName+"-"+row+".png");
				}
				str=br.readLine(); row++;
			}
		}
	}

}