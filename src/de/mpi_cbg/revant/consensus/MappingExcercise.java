package de.mpi_cbg.revant.consensus;

import java.io.*;
import java.util.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.ReadA;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Factorize;


/**
 * Used in the experiment in which we mapped my lungfish modules to 1x of the reads.
 * This is the first step of the analysis. Second step: MappingExcerciseCollisions.java
 * Computes collision graph and coverage histogram.
 * Does not use alignments with too high error rate.
 */
public class MappingExcercise {
		
	
	private static Factorize.StatsWindow[] statWindows;
	private static int[] statReads;
	private static int lastStatWindow;
		
		
	/**
	 * java MappingExcercise FILE_ID N_ALIGNMENTS ROOT_DIR MAX_ERROR_RATE 1> collisions.txt 2> coverageHistogram.txt 
	 */
	public static void main(String[] args) throws Exception {
		final String FILE_ID = args[0];
		final int N_ALIGNMENTS = Integer.parseInt(args[1]);
		final String ROOT_DIR = args[2];
		final double MAX_ERROR_RATE = Double.parseDouble(args[3]);
		final int N_READS = 1868764;
		final int MAX_READ_LENGTH = 500000;  // Arbitrary
		final String LENGTHS_FILE = ROOT_DIR+"/readLengths.txt";
		final int MAX_ALIGNMENTS_PER_READ = 10000;  // Arbitrary
		final String ALIGNMENTS_FILE = ROOT_DIR+"/alignments-"+FILE_ID+"-clean.txt";
		final String QUALITY_THRESHOLDS_FILE = ROOT_DIR+"/qualityThresholds.txt";
		final String MODULES_LENGTHS_FILE = ROOT_DIR+"/modulesLengths.txt";
		final int N_MODULES = 2020;
		int i, j;
		int id, nReads, nReadsWithJunctions, nCollisions, readLength;
		long surface, coveredSurface;
		BufferedReader br;
		int[] tmpArray = new int[3];
		int[] modulesLengths;
		int[] coverageHistogram = new int[101];
		int[] lastCollision = new int[N_MODULES];
		int[][] collisions = new int[N_MODULES][8];
		
		// Loading input
		Factorize.loadInput(new String[] {N_READS+"",MAX_READ_LENGTH+"",LENGTHS_FILE,N_ALIGNMENTS+"",ALIGNMENTS_FILE,QUALITY_THRESHOLDS_FILE,"null","1","3","0","1000","1000","0","null","null","null","null","null","null","null","0","null"});
		ReadA.allocateMemory(MAX_READ_LENGTH,MAX_ALIGNMENTS_PER_READ);
		modulesLengths = new int[N_MODULES];
		br = new BufferedReader(new FileReader(MODULES_LENGTHS_FILE));	
		for (i=0; i<N_MODULES; i++) modulesLengths[i]=Integer.parseInt(br.readLine());
		br.close(); br=null;
		
		// Computing statistics
		Math.set(coverageHistogram,coverageHistogram.length-1,0);
		Math.set(lastCollision,lastCollision.length-1,-1);
		i=0; statWindows=null; lastStatWindow=-1; 
		coveredSurface=0; surface=0; nReads=0; nReadsWithJunctions=0; nCollisions=0;
		while (i<Alignments.nAlignments) {
			ReadA.initialize(i,false,false);
			nReads++;
			readLength=Reads.getReadLength(ReadA.id);
			surface+=readLength;
			
			// Collecting windows
			if (statWindows==null) {
				statWindows = new Factorize.StatsWindow[ReadA.lastSortedAlignment+1];
				for (i=0; i<=ReadA.lastSortedAlignment; i++) statWindows[i] = new Factorize.StatsWindow();
			}
			if (statWindows.length<ReadA.lastSortedAlignment+1) {
				Factorize.StatsWindow[] newStatWindows = new Factorize.StatsWindow[ReadA.lastSortedAlignment+1];
				System.arraycopy(statWindows,0,newStatWindows,0,lastStatWindow+1);
				for (i=lastStatWindow+1; i<newStatWindows.length; i++) newStatWindows[i] = new Factorize.StatsWindow();
				statWindows=newStatWindows;
			}
			lastStatWindow=-1;
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				id=ReadA.sortedAlignments[i].id;
				if (Alignments.getAvgDiffs(id)>MAX_ERROR_RATE) continue;
				lastStatWindow++;
				statWindows[lastStatWindow].start=Alignments.alignments[id][3];
				statWindows[lastStatWindow].end=Alignments.alignments[id][4];
				statWindows[lastStatWindow].tag=Alignments.alignments[id][1]-1;
				statWindows[lastStatWindow].tagLength=modulesLengths[statWindows[lastStatWindow].tag];
				statWindows[lastStatWindow].otherStart=Alignments.alignments[id][5];
				statWindows[lastStatWindow].otherEnd=Alignments.alignments[id][6];
			}
			if (lastStatWindow>0) Arrays.sort(statWindows,0,lastStatWindow+1);
			
			// Computing statistics
			if (lastStatWindow>=0) {
				Factorize.assemblyStatistics_impl(statWindows,lastStatWindow,tmpArray,collisions,lastCollision);
				coverageHistogram[(int)(100*(((double)tmpArray[0])/readLength))]++;
				coveredSurface+=tmpArray[0];
				if (tmpArray[1]>0) nReadsWithJunctions++;
				nCollisions+=tmpArray[2];
			}
			i=ReadA.lastAlignment+1;
		}
		//System.out.println(surface+","+coveredSurface+","+nReads+","+nReadsWithJunctions+","+nCollisions);
		for (i=0; i<N_MODULES; i++) {
			System.out.print(i+",");
			for (j=0; j<=lastCollision[i]; j++) System.out.print(collisions[i][j]+",");
			System.out.println();
		}
		for (i=0; i<coverageHistogram.length; i++) System.err.println(i+","+coverageHistogram[i]);
	}

}