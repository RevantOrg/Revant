package de.mpi_cbg.revant.consensus;

import java.io.*;
import java.util.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Colors;


/**
 * Used in the experiment in which we mapped my lungfish modules to 1x of the reads.
 * Second step, after MappingExcercise.java.
 * Collects all collision files.
 * Shows statistics on "collisions", defined in Factorize.assemblyStatistics_impl().
 */
public class MappingExcerciseCollisions {
		
	/**
	 * java MappingExcerciseCollisions ROOT_DIR
	 */
	public static void main(String[] args) throws Exception {
		String ROOT_DIR = args[0];
		int N_BLOCKS = 120;
		int N_MODULES = 2020;
		boolean found;
		int i, j, k;
		int tag, count, edgeSize;
		String str;
		BufferedReader br;
		String[] tokens;
		int[] lastCollision = new int[N_MODULES];
		int[][] collisions = new int[N_MODULES][8];
		boolean[] printLabel = new boolean[N_MODULES];
		
		// Cumulating all counts
		collisions = new int[N_MODULES][8];
		lastCollision = new int[N_MODULES];
		Math.set(lastCollision,N_MODULES-1,-1);
		for (i=0; i<N_BLOCKS; i++) {
			br = new BufferedReader(new FileReader(ROOT_DIR+"/stats-"+(i+1)+".txt"));	
			for (i=0; i<N_MODULES; i++) {
				str=br.readLine();
				tokens=str.split(",");
				for (j=1; j<tokens.length; j+=2) {
					tag=Integer.parseInt(tokens[j]);
					count=Integer.parseInt(tokens[j+1]);
					found=false;
					for (k=0; k<=lastCollision[i]; k+=2) {
						if (collisions[i][k]==tag) {
							found=true;
							collisions[i][k+1]++;
							break;
						}
					}
					if (found) continue;
					if (lastCollision[i]+2>=collisions[i].length) {
						int[] tmpArray = new int[collisions[i].length<<1];
						System.arraycopy(collisions[i],0,tmpArray,0,collisions[i].length);
						collisions[i]=tmpArray;
					}
					lastCollision[i]++;
					collisions[i][lastCollision[i]]=tag;
					lastCollision[i]++;
					collisions[i][lastCollision[i]]=count;
				}
			}
			br.close();
		}
		
		// Loading names file
		String[] moduleNames = new String[N_MODULES];
		br = new BufferedReader(new FileReader(ROOT_DIR+"/../moduleNames.txt"));	
		for (i=0; i<N_MODULES; i++) moduleNames[i]=br.readLine();
		br.close();		
		
		// Printing graph
		Math.set(printLabel,N_MODULES-1,false);
		for (i=0; i<N_MODULES; i++) {
			if (lastCollision[i]>=0) printLabel[i]=true;
			for (j=0; j<=lastCollision[i]; j+=2) printLabel[collisions[i][j]]=true;
		}
		final double PEN_SCALE = 1.0;
		final String fillColor = Colors.type2color_dot(Constants.INTERVAL_ALIGNMENT);
		System.out.println("graph G {");
		System.out.println("bgcolor=\""+Colors.COLOR_DOT_BACKGROUND+"\";\n");
		for (i=0; i<N_MODULES; i++) {
			if (printLabel[i]) System.out.println(i+" [shape=\"circle\",style=\"filled\",label=\""+moduleNames[i]+"\",color=\""+fillColor+"\",fontcolor=\""+Colors.COLOR_DOT_TEXT+"\",fontname=\"Helvetica\"];\n");
		}
		for (i=0; i<N_MODULES; i++) {
			for (j=0; j<=lastCollision[i]; j+=2) {
				edgeSize=Math.max(1,(int)(Math.sqrt(collisions[i][j+1])*PEN_SCALE));
				System.out.println(i+" -- "+collisions[i][j]+" [arrowhead=\"none\",arrowtail=\"none\",penwidth=\""+edgeSize+"\",color=\""+Colors.COLOR_DOT_EDGE+"\"];\n");
			}
		}
		System.out.println("}");
	}

}