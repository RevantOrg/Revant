package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Biology;

/** 
 * Builds histograms of properties of intervals
 */
public class IntervalStats {
	
	/**
	 * Temporary space used by all procedures
	 */
	private static int[] histogramX, histogramY, histogramY2;
	private static int nTicks;
	
	
	/**
	 * maxReadLength=44766
	 */
	public static void main(String[] args) throws IOException {
		final String PERIODIC_FILE = args[0];
		final String DENSE_FILE = args[1];
		final String ALIGNMENT_FILE = args[2];
		final int MAX_LENGTH = Integer.parseInt(args[3]);
		BufferedReader periodicFile, denseFile, alignmentFile;
/*		
		// Periodic substrings: periods.
		periodicFile = new BufferedReader(new FileReader(PERIODIC_FILE));
		initPeriodHistogram();
		getHistogram(periodicFile,4);
		periodicFile.close();
		System.out.println("Periods:"); printHistogram(histogramY);
		
		// Periodic substrings: lengths.
		periodicFile = new BufferedReader(new FileReader(PERIODIC_FILE));
		initLengthHistogram(MAX_LENGTH);
		getLengthHistogram(periodicFile,3,2,5,6);
		periodicFile.close();
		System.out.println("Periodic lengths, whole:"); printHistogram(histogramY);
		System.out.println("Periodic lengths, fragments:"); printHistogram(histogramY2);
*/				
		// Dense substrings: periods.
		denseFile = new BufferedReader(new FileReader(DENSE_FILE));
		periodicFile = new BufferedReader(new FileReader(PERIODIC_FILE));
		initPeriodHistogram("");
		getDensePeriodsHistogram(denseFile,periodicFile,0,3,4,2,3,10);
		periodicFile.close(); denseFile.close();
		System.out.println("Strong periods of dense substrings:"); printHistogram(histogramY);
		System.out.println("Non-strong periods of dense substrings:"); printHistogram(histogramY2);
/*		
		// Dense substrings: lengths.
		denseFile = new BufferedReader(new FileReader(DENSE_FILE));
		initLengthHistogram(MAX_LENGTH);
		getLengthHistogram(denseFile,4,3,21,22);
		denseFile.close();
		System.out.println("Dense lengths, whole:"); printHistogram(histogramY);
		System.out.println("Dense lengths, fragments:"); printHistogram(histogramY2);
				
		// Alignment interval lengths
		alignmentFile = new BufferedReader(new FileReader(ALIGNMENT_FILE));
		initLengthHistogram(MAX_LENGTH);
		getLengthHistogram(alignmentFile,3,2,6,7);
		alignmentFile.close();
		System.out.println("Alignment interval lengths, whole:"); printHistogram(histogramY);
		System.out.println("Alignment interval lengths, fragments:"); printHistogram(histogramY2);
*/	}
	
	
	/**
	 * Stores in $histogramY$ (respectively, $histogramY2$) the distribution of strong 
	 * (respectively, non-strong) periods for dense substring intervals that are not 
	 * contained in, or approximately identical to, a periodic substring interval.
	 */
	private static final void getDensePeriodsHistogram(BufferedReader dense, BufferedReader periodic, int readCell, int startCellDense, int endCellDense, int startCellPeriodic, int endCellPeriodic, int periodCellDense) throws IOException {
		final int READ_AHEAD_LIMIT = 1000000;
		boolean contained;
		int i, read, readDense, readPeriodic, number, pos;
		int startDense, endDense, startPeriodic, endPeriodic, strongPeriodSignal;
		String strDense, strPeriodic;
		int[] hist;
		String[] tokensDense, tokensPeriodic;
		
		Math.set(histogramY,histogramY.length-1,0);
		Math.set(histogramY2,histogramY2.length-1,0);
		strDense=dense.readLine();
		readPeriodic=-1;
		while (strDense!=null) {
			tokensDense=strDense.replace('|',',').split(",");
			readDense=Integer.parseInt(tokensDense[readCell]);
			while (readPeriodic<readDense) {
				periodic.mark(READ_AHEAD_LIMIT);
				strPeriodic=periodic.readLine();
				if (strPeriodic==null) {
					tokensPeriodic=null;
					periodic.reset();
					break;
				}
				tokensPeriodic=strPeriodic.replace('|',',').split(",");
				readPeriodic=Integer.parseInt(tokensPeriodic[readCell]);		
			}
			if (readPeriodic!=readDense) {
				number=Integer.parseInt(tokensDense[periodCellDense]);
				strongPeriodSignal=Integer.parseInt(tokensDense[periodCellDense+1]);
				if (strongPeriodSignal==1) hist=histogramY;
				else hist=histogramY2;
				if (number<=0) hist[0]++;
				else {
					pos=Arrays.binarySearch(histogramX,0,nTicks,number);
					if (pos<0) {
						pos=-pos-1;
						pos=Math.min(pos,nTicks-1);
					}
					hist[pos]++;
				}
if (strongPeriodSignal==1) System.err.println("read "+readDense+" has a dense substring with strong period "+number);				
			}
			else {
				startDense=Integer.parseInt(tokensDense[startCellDense]);
				endDense=Integer.parseInt(tokensDense[endCellDense]);
				periodic.reset();
				read=readDense;
				contained=false;
				while (read==readDense) {
					strPeriodic=periodic.readLine();
					if (strPeriodic==null) {
						tokensPeriodic=null;
						break;
					}
					tokensPeriodic=strPeriodic.replace('|',',').split(",");
					read=Integer.parseInt(tokensPeriodic[readCell]);
					if (read!=readPeriodic) break;
					startPeriodic=Integer.parseInt(tokensPeriodic[startCellPeriodic]);
					endPeriodic=Integer.parseInt(tokensPeriodic[endCellPeriodic]);
					if ( Intervals.areApproximatelyIdentical(startDense,endDense,startPeriodic,endPeriodic) || 
						 Intervals.isApproximatelyContained(startDense,endDense,startPeriodic,endPeriodic) ) {
						contained=true;
						break;
					}
				}
				if (!contained) {
					number=Integer.parseInt(tokensDense[periodCellDense]);
					strongPeriodSignal=Integer.parseInt(tokensDense[periodCellDense+1]);
					if (strongPeriodSignal==1) hist=histogramY;
					else hist=histogramY2;
					if (number<=0) hist[0]++;
					else {
						pos=Arrays.binarySearch(histogramX,0,nTicks,number);
						if (pos<0) {
							pos=-pos-1;
							pos=Math.min(pos,nTicks-1);
						}
						hist[pos]++;
					}
if (strongPeriodSignal==1) System.err.println("read "+readDense+" has a dense substring with strong period "+number);					
				}
			}
			strDense=dense.readLine();
		}
	}

	
	/**
	 * Considers only cell $cell$ of each line of $br$. Updates only $histogramY$.
	 */
	private static final void getHistogram(BufferedReader br, int cell) throws IOException {
		int i, number, pos;
		String str;
		String[] tokens;
		
		Math.set(histogramY,histogramY.length-1,0);
		str=br.readLine();
		while (str!=null) {
			str=str.replace('|',',');
			tokens=str.split(",");
			number=Integer.parseInt(tokens[cell]);
			if (number<=0) histogramY[0]++;
			else {
				pos=Arrays.binarySearch(histogramX,0,nTicks,number);
				if (pos<0) {
					pos=-pos-1;
					pos=Math.min(pos,nTicks-1);
				}
				histogramY[pos]++;
			}
			str=br.readLine();
		}
	}
	
	
	/**
	 * Subtracts the value in cell $cellY$ from the value in cell $cellX$ of each line of 
	 * $br$. The result increments $histogramY$ if both $cellXPrime$ and $cellYPrime$ are
	 * one, otherwise it increment $histogramY2$.
	 */
	private static final void getLengthHistogram(BufferedReader br, int cellX, int cellY, int cellXPrime, int cellYPrime) throws IOException {
		int i, number, pos;
		String str;
		String[] tokens;
		
		Math.set(histogramY,histogramY.length-1,0);
		Math.set(histogramY2,histogramY2.length-1,0);
		str=br.readLine();
		while (str!=null) {
			str=str.replace('|',',');
			tokens=str.split(",");
			number=Integer.parseInt(tokens[cellX])-Integer.parseInt(tokens[cellY])+1;
			pos=Arrays.binarySearch(histogramX,0,nTicks,number);
			if (pos<0) {
				pos=-pos-1;
				pos=Math.min(pos,nTicks-1);
			}
			if (Integer.parseInt(tokens[cellXPrime])==1 && Integer.parseInt(tokens[cellYPrime])==1) histogramY[pos]++;
			else histogramY2[pos]++;
			str=br.readLine();
		}
	}
	
	
	private static final void initPeriodHistogram(String directory) throws IOException {
		int i;
		
		Biology.initialize(directory);
		histogramX = new int[Biology.maxPeriod+1];
		for (i=0; i<histogramX.length; i++) histogramX[i]=i;
		nTicks=Biology.maxPeriod;
		histogramY = new int[nTicks];
		histogramY2 = new int[nTicks];
	}
	
	
	private static final void initLengthHistogram(int maxLength) throws IOException {
		int i;
		
		histogramX = new int[maxLength];
		for (i=0; i<histogramX.length; i++) histogramX[i]=i+1;
		nTicks=maxLength;
		histogramY = new int[nTicks];
		histogramY2 = new int[nTicks];
	}
	
	
	/**
	 * Prints every element of $histogramX$ and $histogram$.
	 */
	private static final void printHistogram(int[] x) {
		for (int i=0; i<nTicks; i++) System.out.println(histogramX[i]+","+x[i]);
	}

	
}