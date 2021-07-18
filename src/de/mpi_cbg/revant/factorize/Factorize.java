package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Biology;
import de.mpi_cbg.revant.util.RegressionTree;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Leaves;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Histograms;

/**
 * Repeat models could be built by a variant of the following greedy algorithm, suggested 
 * by Gene Myers: 
 * (1) Heuristically select a substring W^i of a read that is covered by a large number of
 * alignments, and add W^i to a dictionary. (2) Mark as "substring of W^i" all other 
 * substrings of reads that align to W^i, and hide them from the following steps of the 
 * algorithm. (3) Proceed to finding the next W^{i+1} until some stopping conditions is 
 * satisfied (e.g. all alignments have been used). However, this approach has at least the
 * following problems:
 *
 * 1. Assume that you select as W^i a long substring of a read that contains a complex 
 * configuration of repeats (e.g. W^i = R^1 R^2 ... R^k, where R^j is an elementary 
 * repeat). Then, all other reads in the full (say, 50x) read set that contain W^i will 
 * likely be marked as containing an instance of W^i, and since the information that W^i 
 * consists of k elementary repeats is never computed nor stored in the dictionary (since 
 * we didn't factorize W^i, we just added it greedily to the dictionary), such information
 * is also not annotated to the reads that contain W^i. This could be bad for the 
 * following steps of assembly, since it prevents running heuristics based e.g. on the 
 * number of consecutive repeats. One could argue that, if we just match the full copy 
 * rather than all its substrings, the algorithm would at least create a model for every 
 * distinct repeat, but even this is not true, since if a repeat occurs only inside other 
 * repeats (i.e. it is not near-supermaximal in the read set) it might not get detected
 * as a separate unit.
 *
 * 2. It could happen that W^i is one of the prefixes of a repeat R that propagates by 
 * taking prefixes of itself, but not the longest such prefix, i.e. not the entire R. 
 * Then the algorithm would add W^i to the dictionary, it would correctly mark all the 
 * occurrences of prefixes of W^i as prefixes of W^i, but it would also mark the 
 * occurrences of W^i in all the prefixes of R of length greater than |W^i|. This would 
 * likely disconnect such longer prefixes of R from the prefixes of W^i, i.e. in the end 
 * it could happen that the prefixes of R could be marked as being substrings of many 
 * different dictionary words, rather than of just one dictionary word that equals R. 
 * In addition to being intuitively a suboptimal annotation, this would have the effect of
 * annotating, in a read, multiple consecutive repeat units rather than just one. The 
 * number of such units would depend on the heuristic, rather than on the information in 
 * the actual data. This again could prevent running heuristics based on the number of 
 * consecutive repeats in a complex region. One could argue that these problems do not 
 * happen if we greedily process alignments from longest to shortest.
 *
 * 3. Similarly, if W^i is a substring of a repeat R that propagates by taking substrings 
 * of itself, but W^i is not the entire R, then the occurrences of the whole R in reads 
 * are broken into the occurrence of W^i plus a prefix and a suffix. The fate of such 
 * prefix/suffix depends on what the algorithm exactly does, but such regions will likely
 * not be related to W^i in the final annotation of the full (say, 50x) read set, leading 
 * to the annotation of multiple artificial repeats rather than of just one, as above.
 * Again, one could argue that these problems do not occur if we greedily process 
 * alignments from longest to shortest.
 *
 * 4. (minor) Given W^i in read R, marking the occurrences of W^i in every other read S 
 * using just the alignment of R with S might not be accurate. I.e. it could be more 
 * accurate to use all the alignments of read S to estimate the first and last position of
 * the occurrence of W^i in S. 
 *
 * 5. Assume that different substrings of a repeat occur in different reads, such 
 * substrings can be assembled to get the full sequence of the repeat, but the full 
 * sequence never occurs. We would annotate such substrings as separate repeats.
 *
 * Remark: the idea of detecting repeats de novo by splitting a reference genome into  
 * short (20Kbp) overlapping substrings and computing all pairwise alignments between such 
 * substrings appeared in \cite{kaminker2002transposable}. The idea of then clustering 
 * such alignments into repeat families, and to compute the multiple alignment of each
 * cluster, appeared in \cite{flutre2011considering}. All pairwise  
 * alignments have also been used in \cite{koren2012hybrid}, but for long-read correction  
 * and just between a set of short reads and a set of long reads; in 
 * \cite{she2004structure} on the human reference, for detecting repeats near centromeres 
 * at different levels of similarity; in \cite{skaletsky2003male}, to detect tandem arrays 
 * in the Y chromosome (see e.g. Figure 5); in \cite{riethman2004mapping}, to detect 
 * segmental duplications in subtelomeric regions of human. All pairwise alignments 
 * between alpha monomers where used in \cite{alkan2007organization} to cluster them in 
 * families.
 *
 * Remark: the idea of assembling centromeres by taking advantage of the evident 
 * heterogeneity of monomeric repeats and the periodic variants that punctuate otherwise 
 * homogeneous arrays of higher-order alpha satellites, was already present in  
 * \cite{rudd2004analysis}.
 *
 * Remark: the idea of masking repeats before assembly, albeit with a manually curated 
 * database of known repeats, already appeared in \cite{myers2000whole}. A claim that
 * detecting repeats before assembly is useful, and that reads with repeats should be
 * processed with a different pipeline, appears also in \cite{chen2019new,
 * kolmogorov2018assembly}. The idea of building repeat models iteratively, by first  
 * detecting repeats from a subset of reads, then masking them and all their occurrences, 
 * and iterating using more reads, already appears in \cite{peterson2013plant} (page 280),
 * in \cite{price2005denovo}, and probably also in RECON.
 *
 * Remark: detecting tandem repeats from long reads has already been attempted in 
 * \cite{ummat2014resolving}, but using a reference.
 * 
 * Remark: our overall approach of binning all reads that belong to a repeat class and  
 * then using their sequence variation to infer subclasses is similar to 
 * \cite{horvath2000mosaic}, even though they proceed biologically as well.
 *
 * Remark: our notion of maximality of alignments is related to the notion of "penalty  
 * zone" in \cite{bandyopadhyay2010repfrag}.
 *
 * Remark: paper \cite{kolmogorov2018assembly} is also essentially factorizing repeats, by 
 * single-linkage clustering of their endpoints (although there is no difference between  
 * maximal and non-maximal endpoints).
 *
 * Remark: in a typical assembly pipeline, raw reads are "patched" before assembly, i.e.: 
 * (1) substrings with low intrinsic quality are replaced with corresponding regions of
 * high quality from other reads; (2) some chimeric reads, that are the concatenation of
 * different reads or of different parts of the same read, are cut. Our repeat detection
 * pipeline models (1) explicitly, but it does not model (2) explicitly. If the chimer
 * involves a repeat, then: in the chimeric read, factorization could create two repeat
 * intervals; in the interval graph, such intervals are likely to create insertion edges 
 * with the other intervals of the same repeat, and these should not affect the rest of 
 * the pipeline. In every other read that contains the repeat, the chimer would induce <=2
 * maximal alignments, which are unlikely to affect the factorization of those reads.
 * Patching is usually an expensive process, so repeat modeling should be done on 
 * unpatched reads if possible.
 */
public class Factorize {
	/*
	 * Output file and codes
	 */
	public static final int CODE_NO_FACTOR = -1;
	public static final int CODE_NO_FACTOR_WORTH_REPORTING = -2;
	private static BufferedWriter factorsFile, intervalsDenseFile, intervalsPeriodicFile, intervalsAlignmentFile;
	private static BufferedWriter denseSubstringsFile, periodicSubstringsFile, alignmentsFile;
	private static BufferedWriter intervalConnectionFile;
	
	/**
	 * Temporary space loaded by procedure $readConnectionFile$.
	 */	
	public static boolean isWeak, isLeftMaximal, isRightMaximal, hasLongPeriod;
	public static int type, id, start, end, implied, period, minPrefixLength, minSuffixLength, isContained;
	public static int firstMaximalStart, lastMaximalEnd;
	public static int nAssignedAlignments;
	public static double avgCoverage;
	
	/**
	 * Temporary space
	 */
	private static int[] tmpArray = new int[20];
	private static int[] tmpArrayPrime = new int[20];
	private static PeriodicSubstringInterval tmpPeriodic = new PeriodicSubstringInterval();
	private static DenseSubstring tmpDense = new DenseSubstring();
	private static AlignmentInterval tmpAlignment = new AlignmentInterval();
	private static StatsWindow[] statWindows;  // For assembly statistics
	private static int[] statReads;  // Reads reported by assembly statistics
	private static int lastStatRead;


	/**
	 * @param args
	 *  0: nReads;
	 *  1: max length of a read;
	 *  2: read lengths file;
	 *  3: read IDs file;
	 *  4: nAlignments;
	 *  5: alignments file;
	 *  6: quality thresholds file;
	 *  7: qualities file;
	 *  8: coverage of the genome. Remark: this is the number of bases that cover a base
	 *     of a chromosome. So if the organism is N-ploid, the coverage is the usual 
	 *     notion of overage divided by N.
	 *  9: minimum n. of occurrences in the genome for a substring to be considered a 
	 *     repeat; increasing this was suggested by Gene Myers and Martin Pippel, e.g. to
	 *     implement an incremental pipeline that works just on high-frequency repeats at
	 *     each stage. The idea of keeping just fragments with high coverage was already 
	 *     present in e.g. \cite{chu2016repdenovo} and references.
	 * 10: alignment file includes self-alignments? (1/0)
	 * 11: minAlignmentLength;
	 * 12: minFactorLength;
	 * 13: minLongInsertionLength;
	 * 14: path of the factors file;
	 * 15: path of the factors-intervals-dense file;
	 * 16: path of the factors-intervals-periodic file;
	 * 17: path of the factors-intervals-alignment file;
	 * 18: path of the intervals-dense file;
	 * 19: path of the intervals-periodic file;
	 * 20: path of the intervals-alignments file;
	 * 21: output all periodic factors? (1/0);
	 * 22: path of the interval-connection file.
	 */
	public static final void loadInput(String[] args) throws IOException {
		final String NULL_FILE = "null";
		boolean includesSelfAlignments;

		Reads.nReads=Integer.parseInt(args[0]);
		Reads.maxReadLength=Integer.parseInt(args[1]);
		Reads.loadReadLengths(args[2]);
		Reads.loadReadIDs(args[3],Reads.nReads);
		Alignments.nAlignments=Integer.parseInt(args[4]);
		Alignments.minAlignmentLength=Integer.parseInt(args[11]);
		if (args[5].length()>0 && !args[5].equals(NULL_FILE)) Alignments.loadAlignments(args[5]);
		if (args[6].length()>0 && !args[6].equals(NULL_FILE) && args[7].length()>0) Reads.loadQualities(args[6],args[7].equalsIgnoreCase(NULL_FILE)?null:args[7]);

		IO.coverage=Integer.parseInt(args[8]);
		IO.minCoverage=Math.max(IO.coverage-1,1);
		IO.minOccurrencesInGenome=Integer.parseInt(args[9]);		
		IO.includeSelfAlignments=Integer.parseInt(args[10])==1;
		if (IO.includeSelfAlignments) IO.minCoverage++;
		IO.minRepeatCoverage=(IO.coverage*IO.minOccurrencesInGenome)-1;
		if (IO.includeSelfAlignments) IO.minRepeatCoverage++;

		Factors.minFactorLength=Integer.parseInt(args[12]);
		Factors.maxFactorsPerRead=Reads.maxReadLength>>1;
		Reads.minLongInsertionLength=Integer.parseInt(args[13]);
		if (args[14].length()>0 && !args[13].equals(NULL_FILE)) factorsFile = new BufferedWriter(new FileWriter(args[14]),IO.BUFFER_SIZE);
		if (args[15].length()>0 && !args[15].equals(NULL_FILE)) intervalsDenseFile = new BufferedWriter(new FileWriter(args[15]),IO.BUFFER_SIZE);
		if (args[16].length()>0 && !args[16].equals(NULL_FILE)) intervalsPeriodicFile = new BufferedWriter(new FileWriter(args[16]),IO.BUFFER_SIZE);
		if (args[17].length()>0 && !args[17].equals(NULL_FILE)) intervalsAlignmentFile = new BufferedWriter(new FileWriter(args[17]),IO.BUFFER_SIZE);
		if (args[18].length()>0 && !args[18].equals(NULL_FILE)) denseSubstringsFile = new BufferedWriter(new FileWriter(args[18]),IO.BUFFER_SIZE);
		if (args[19].length()>0 && !args[19].equals(NULL_FILE)) periodicSubstringsFile = new BufferedWriter(new FileWriter(args[19]),IO.BUFFER_SIZE);
		if (args[20].length()>0 && !args[20].equals(NULL_FILE)) alignmentsFile = new BufferedWriter(new FileWriter(args[20]),IO.BUFFER_SIZE);
		if (Integer.parseInt(args[21])==1) Factors.minPeriodicCoverage=0;
		else Factors.minPeriodicCoverage=IO.minCoverage;
		if (args[22].length()>0 && !args[22].equals(NULL_FILE)) intervalConnectionFile = new BufferedWriter(new FileWriter(args[22]),IO.BUFFER_SIZE);
	}


	/**
	 * Every line in the output encodes a factor, with format $R,W$, where
	 * $R$ is the id of the read, and $W$ is the representation of a factor returned by
	 * procedure $Factor.toString$. If no factor is detected, line $R,CODE_NO_FACTOR$ is
	 * printed. If some factors are detected, but none of them deserves to be reported,
	 * line $R,CODE_NO_FACTOR_WORTH_REPORTING,p,d,r$ is printed, where $p$ is the number
	 * of periodic substrings found, $d$ is the number of dense substrings found, $r$ is
	 * the number of parentheses found.
	 */
	public static void main(String[] args) throws IOException {
		final int MAX_NEIGHBORS = 7000;
		int i, j, k;
		int nAlignments = 0;
		int nParentheses, nUnassignedAlignments;
		long time, factorizationTime;
		long[] assemblyStats = new long[20];
		Math.set(assemblyStats,assemblyStats.length-1,0);
		statWindows=null; statReads=null; lastStatRead=-1;
		System.err.print("Loading input... ");
		time=System.currentTimeMillis();
		final String BIOLOGY_DIR = args[23];
		loadInput(args);
		System.err.println("done in "+((System.currentTimeMillis()-time)/1000)+" seconds");

		// Allocating memory
		System.err.print("Allocating memory...");
		time=System.currentTimeMillis();
		Biology.initialize(BIOLOGY_DIR);
		IO.initialize();
		RegressionTree.allocateMemory(Reads.maxReadLength,Factors.minFactorLength);
		System.err.print(" 1 ");
		DensityEstimationTree.allocateMemory(Alignments.maxAlignmentsPerRead<<1);  // Since it is used also for estimating the density of events, alignment qualities, etc.
		System.err.print(" 2 ");
		Leaves.allocateMemory(Math.max(RegressionTree.leaves.length,DensityEstimationTree.leaves.length));
		System.err.print(" 3 ");
		Intervals.allocateMemory(Math.max(RegressionTree.leaves.length,DensityEstimationTree.leaves.length));
		System.err.print(" 4 ");
		ReadA.allocateMemory(Reads.maxReadLength,Alignments.maxAlignmentsPerRead);
		System.err.print(" 5 ");
		Events.allocateMemory(10/*Arbitrary, reallocated later.*/,Factors.maxFactorsPerRead<<1);
		System.err.print(" 6 ");
		PeriodicSubstrings.allocateMemory(Alignments.maxAlignmentsPerRead,Alignments.maxOccurrences>>1/* Division by two arbitrary */,Factors.maxFactorsPerRead<<1,10/*Arbitrary, reallocated later.*/);
		System.err.print(" 7 ");
		DenseSubstrings.allocateMemory(Alignments.maxAlignmentsPerRead,Math.ceil(Alignments.maxAlignmentsPerRead,DenseSubstrings.DEFAULT_MIN_PATH_LENGTH>>1),Factors.maxFactorsPerRead<<1,10/*Arbitrary, reallocated later.*/,Reads.maxReadLength);
		System.err.print(" 8 ");
		AlignmentIntervals.allocateMemory(Alignments.maxAlignmentsPerRead,10/*Arbitrary, reallocated later.*/);
		System.err.print(" 9 ");
		Factors.allocateMemory(Alignments.maxAlignmentsPerRead,Factors.maxFactorsPerRead,Alignments.maxOccurrences,Factors.minFactorLength>>2);
		System.err.print(" 10 ");
		Points.allocateMemory(Alignments.maxAlignmentsPerRead<<1);
		System.err.print(" 11 ");
		Histograms.allocateMemory(Reads.maxReadLength);
		System.err.println("done in "+((System.currentTimeMillis()-time)/1000)+" seconds  maxAlignmentsPerRead="+Alignments.maxAlignmentsPerRead);

		// Factoring
		System.err.println("Factoring reads... ");
		time=System.currentTimeMillis();
		i=0;
		int alignmentCounter = 0;
		while (i<Alignments.nAlignments) {
			// Initializing $readA$
			ReadA.initialize(i,false,true);
			System.err.println("Read "+ReadA.id+" started: L="+Reads.getReadLength(ReadA.id)+", A="+(ReadA.lastSortedAlignment+1));
			assemblyStats[11]++;



			/*
final int SELECTED_READ = 1754320;  // Starting from zero
if (ReadA.id>SELECTED_READ) break;
if (ReadA.id<SELECTED_READ) {
 	i=ReadA.lastAlignment+1;
 	continue;
}
*/


			ReadA.getCoverageHistogram(true);  // Used by $Factors.mergeFactors$
			if (!Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),0,Reads.getQualityArrayLength(ReadA.id)-1,Reads.MAX_HIGH_QUALITY_SCORE+1,false,(Alignments.minAlignmentLength>>1)/Reads.QUALITY_SPACING,true,-1)) {
				factorsFile.write(ReadA.id+",");
				factorsFile.write(CODE_NO_FACTOR+"\n");
				intervalsDenseFile.write("\n");
				intervalsPeriodicFile.write("\n");
				intervalsAlignmentFile.write("\n");
				for (j=0; j<=ReadA.lastSortedAlignment; j++) intervalConnectionFile.write("\n");
				i=ReadA.lastAlignment+1;
				continue;
			}
			alignmentCounter+=ReadA.lastSortedAlignment+1;
			ReadA.print();

			// Detecting events and splits
			Events.lastEvent=-1;
			Event.order=Event.UNSORTED;
			Factors.lastSplit=-1;
			Split.order=Split.UNSORTED;
			PeriodicSubstrings.detect();

if (IO.SHOW_STD_ERR_PRIME) { 
	IO.printErr("Final set of periodic substrings:");
	for (int x=0; x<=PeriodicSubstrings.lastSubstring; x++) {
		IO.printErr("periodic substring: "+PeriodicSubstrings.substrings[x]); 
		IO.printErr("shifts:");
		for (int y=0; y<=PeriodicSubstrings.substrings[x].lastShift; y++) IO.printErr(PeriodicSubstrings.substrings[x].shifts[y]);
	}
	IO.printErr("alignments after PeriodicSubstrings.detect():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
	System.err.println("Events PeriodicSubstrings.detect():");
	for (int x=0; x<=Events.lastEvent; x++) System.err.println(Events.events[x]);
	System.err.println("Splits PeriodicSubstrings.detect():");
	for (int x=0; x<=Factors.lastSplit; x++) System.err.println(Factors.splits[x]);
}

probe(1);


			DenseSubstrings.detect();
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignments after DenseSubstrings.detect():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
	System.err.println("Events DenseSubstrings.detect():");
	for (int x=0; x<=Events.lastEvent; x++) System.err.println(Events.events[x]);
	System.err.println("Splits DenseSubstrings.detect():");
	for (int x=0; x<=Factors.lastSplit; x++) System.err.println(Factors.splits[x]);
}			

probe(2);

			
			AlignmentIntervals.detect();
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignments after AlignmentIntervals.detect():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
	System.err.println("Events AlignmentIntervals.detect():");
	for (int x=0; x<=Events.lastEvent; x++) System.err.println(Events.events[x]);
	System.err.println("Splits AlignmentIntervals.detect():");
	for (int x=0; x<=Factors.lastSplit; x++) System.err.println(Factors.splits[x]);
}
probe(3);
				
			nParentheses=Parentheses.detect(false);
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("Events Parentheses.detect():");
	for (int x=0; x<=Events.lastEvent; x++) System.err.println(Events.events[x]);
	System.err.println("Splits Parentheses.detect():");
	for (int x=0; x<=Factors.lastSplit; x++) System.err.println(Factors.splits[x]);
}
probe(4);			
//			discardRareIntervals();

			if ( Events.lastEvent==-1 && Factors.lastSplit==-1 &&
			     AlignmentIntervals.lastInterval==-1 && DenseSubstrings.lastSubstring==-1 && PeriodicSubstrings.lastInterval==-1
			   ) {
				// No maximal event, no non-maximal interval, no B-maximal dense or
				// periodic substring.
				factorsFile.write(ReadA.id+",");
				factorsFile.write(CODE_NO_FACTOR+"\n");
				intervalsDenseFile.write("\n");
				intervalsPeriodicFile.write("\n");
				intervalsAlignmentFile.write("\n");
				for (j=0; j<=ReadA.lastSortedAlignment; j++) intervalConnectionFile.write("\n");
				PeriodicSubstrings.undo_mergeIntervals();
				i=ReadA.lastAlignment+1;
				continue;
			}
			
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("Splits before detectRandomInsertions():");
	for (int x=0; x<=Factors.lastSplit; x++) System.err.println(Factors.splits[x]);		
}

			detectRandomInsertions();
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignment intervals after detectRandomInsertions():");
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
	System.err.println("Splits after detectRandomInsertions():");
	for (int x=0; x<=Factors.lastSplit; x++) System.err.println(Factors.splits[x]);			
}
probe(5);
			// Clustering events, and adding them to $Factors.splits$.
			if (Events.lastEvent>=0) Events.cluster();


if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("SPLITS:");
	for (int x=0; x<=Factors.lastSplit; x++) IO.printErr(Factors.splits[x]); 
}
			// Estimating and printing factors			
			Factors.detect();
			if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION || PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
				PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
				PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
				if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
			}
			AlignmentIntervals.markIntervalsInPeriodic(0,AlignmentIntervals.lastInterval,false);  // Both long and short periods
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignment intervals after Factors.detect(): lastFactor="+Factors.lastFactor);
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("alignments after Factors.detect():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
}
probe(6);
			if (Factors.lastFactor==-1) {
				factorsFile.write(ReadA.id+",");
				factorsFile.write(CODE_NO_FACTOR_WORTH_REPORTING+",");
				factorsFile.write((PeriodicSubstrings.lastSubstring+1)+",");
				factorsFile.write((DenseSubstrings.lastSubstring+1)+",");
				factorsFile.write(nParentheses+"\n");
				intervalsDenseFile.write("\n");
				intervalsPeriodicFile.write("\n");
				intervalsAlignmentFile.write("\n");
				for (j=0; j<=ReadA.lastSortedAlignment; j++) intervalConnectionFile.write("\n");
				PeriodicSubstrings.undo_mergeIntervals();
				i=ReadA.lastAlignment+1;
				continue;
			}
			for (j=0; j<=Factors.lastFactor; j++) {
				factorsFile.write(ReadA.id+",");
				Factors.factors[j].writeFactorsFile(factorsFile);
				Factors.factors[j].writeIntervalsDenseFile(intervalsDenseFile);
				Factors.factors[j].writeIntervalsPeriodicFile(intervalsPeriodicFile);
				Factors.factors[j].writeIntervalsAlignmentFile(intervalsAlignmentFile);
			}


if (IO.SHOW_STD_ERR_PRIME) { 
	IO.printErr("Dense substrings before intervals2denseSubstrings():");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) {
		IO.printErr("** dense substring: "+DenseSubstrings.substrings[x]); 
	}
	IO.printErr("Alignments not implied or implied by alignment intervals:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) {
		if (!ReadA.sortedAlignments[x].isImplied() || ReadA.sortedAlignments[x].mergedToInterval!=null) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" merged2interval="+ReadA.sortedAlignments[x].mergedToInterval+" ==> "+ReadA.sortedAlignments[x].toStringPointers());
	}
	IO.printErr("alignment intervals before intervals2denseSubstrings():");
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("periodic substring intervals before intervals2denseSubstrings():");
	for (i=0; i<=PeriodicSubstrings.lastInterval; i++) IO.printErr(PeriodicSubstrings.intervals[i]);
	IO.printErr("long-period periodic substring intervals before intervals2denseSubstrings():");
	for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) IO.printErr(PeriodicSubstrings.longPeriodIntervals[i]);
	IO.printErr("Alignments before intervals2denseSubstrings():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" ==> "+ReadA.sortedAlignments[x].toStringPointers());
}
probe(7);
			AlignmentIntervals.intervals2denseSubstrings(DenseSubstrings.stackPrime,PeriodicSubstrings.stack);  // We use ${Dense,Periodic}Substrings.stack$ since we cannot use $AlignmentIntervals.stack$.
			AlignmentIntervals.cleanAlignmentIntervals_afterFixBoundaries();


if (IO.SHOW_STD_ERR_PRIME) { 
	IO.printErr("periodic substring intervals after intervals2denseSubstrings():");
	for (i=0; i<=PeriodicSubstrings.lastInterval; i++) IO.printErr(PeriodicSubstrings.intervals[i]);
	IO.printErr("Dense substrings after intervals2denseSubstrings():");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) {
		IO.printErr("** dense substring: "+DenseSubstrings.substrings[x]); 
	}
	IO.printErr("alignment intervals after intervals2denseSubstrings():");
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("Alignments after intervals2denseSubstrings():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" ==> "+ReadA.sortedAlignments[x].toStringPointers());
}
probe(8);
			
			
			AlignmentIntervals.setMaximality();
			discardRareIntervals();
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Dense substrings after discardRareIntervals():");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) IO.printErr("** dense substring: "+DenseSubstrings.substrings[x]);
	IO.printErr("alignment intervals after discardRareIntervals():");
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("periodic substring intervals after discardRareIntervals():");
	for (i=0; i<=PeriodicSubstrings.lastInterval; i++) IO.printErr(PeriodicSubstrings.intervals[i]);
	IO.printErr("long-period periodic substring intervals after discardRareIntervals():");
	for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) IO.printErr(PeriodicSubstrings.longPeriodIntervals[i]);
	IO.printErr("Factorize> alignments after discardRareIntervals():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
}
probe(9);
			// Final assignment of alignments to intervals			
			DenseSubstrings.fixBoundariesOfSubstringType(tmpArray);
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Factorize> dense substrings after fixBoundariesOfSubstringType():");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	IO.printErr("Alignments after fixBoundariesOfSubstringType():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" ==> "+ReadA.sortedAlignments[x].toStringPointers());
}	
probe(10);
			PeriodicSubstrings.fixBoundariesOfPeriodicIntervals();
			DenseSubstrings.discardPeriodicSubstrings(ReadA.lastSortedAlignment+1);
			AlignmentIntervals.cleanAlignmentIntervals_afterFixBoundaries();
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Factorize> periodic substring intervals after fixBoundariesOfPeriodicIntervals():");
	for (int x=0; x<=PeriodicSubstrings.lastInterval; x++) System.err.println(PeriodicSubstrings.intervals[x]);
	IO.printErr("long-period periodic substring intervals after fixBoundariesOfPeriodicIntervals():");
	for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) IO.printErr(PeriodicSubstrings.longPeriodIntervals[i]);
	IO.printErr("Factorize> dense substrings:");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	IO.printErr("Factorize> alignment intervals:");
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("Alignments after fixBoundariesOfPeriodicIntervals():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" ==> "+ReadA.sortedAlignments[x].toStringPointers());
}				
			markContainedIntervals();
			nUnassignedAlignments=alignments2intervals();
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Factorize> dense substrings after alignments2intervals():");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	IO.printErr("Alignments after alignments2intervals():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" ==> "+ReadA.sortedAlignments[x].toStringPointers());
	IO.printErr("Factorize> periodic substring intervals after alignments2intervals():");
	for (int x=0; x<=PeriodicSubstrings.lastInterval; x++) System.err.println(PeriodicSubstrings.intervals[x]);
	IO.printErr("long-period periodic substring intervals after alignments2intervals():");
	for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) IO.printErr(PeriodicSubstrings.longPeriodIntervals[i]);
}	
probe(11);		
			
			DenseSubstrings.fixBoundariesOfSubstringType(tmpArray);
			PeriodicSubstrings.fixBoundariesOfPeriodicIntervals();
			DenseSubstrings.discardPeriodicSubstrings(ReadA.lastSortedAlignment+1);
			AlignmentIntervals.cleanAlignmentIntervals_afterFixBoundaries();
			nUnassignedAlignments=alignments2intervals();
			/* WAS INSTEAD:
			alignments2intervals_intersection(false);
			nUnassignedAlignments=alignments2intervals_updateAlignmentPointers(true);
			*/

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Factorize> alignments after fixBoundariesOfPeriodicIntervals 2 ():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
	IO.printErr("Factorize> dense substrings after fixBoundariesOfPeriodicIntervals 2 ():");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	IO.printErr("Factorize> alignment intervals after fixBoundariesOfPeriodicIntervals 2 ():");
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("Factorize> periodic substring intervals after fixBoundariesOfPeriodicIntervals 2 ():");
	for (int x=0; x<=PeriodicSubstrings.lastInterval; x++) System.err.println(PeriodicSubstrings.intervals[x]);
	IO.printErr("long-period periodic substring intervals after fixBoundariesOfPeriodicIntervals 2 ():");
	for (int x=0; x<=PeriodicSubstrings.lastLongPeriodInterval; x++) IO.printErr(PeriodicSubstrings.longPeriodIntervals[x]);
	IO.printErr("Alignments after fixBoundariesOfPeriodicIntervals 2 ():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" ==> "+ReadA.sortedAlignments[x].toStringPointers());
}
probe(12);
			if (nUnassignedAlignments>0) generalizeReplicationTypes();
			generalizeReplicationTypes_alignmentIntervals();
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Factorize> dense substrings after generalizeReplicationTypes():");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	IO.printErr("Factorize> periodic substring intervals before discardRareIntervals Z ():");
	for (int x=0; x<=PeriodicSubstrings.lastInterval; x++) System.err.println(PeriodicSubstrings.intervals[x]);
	IO.printErr("long-period periodic substring intervals before discardRareIntervals Z ():");
	for (int x=0; x<=PeriodicSubstrings.lastLongPeriodInterval; x++) IO.printErr(PeriodicSubstrings.longPeriodIntervals[x]);
	IO.printErr("Factorize> alignment intervals before discardRareIntervals Z ():");
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("Factorize> Alignments before discardRareIntervals Z ():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" ==> "+ReadA.sortedAlignments[x].toStringPointers());	
}			
probe(13);
			discardRareIntervals();  // Executed again, since an interval might become rare after $alignments2intervals()$.	
			nUnassignedAlignments=alignments2intervals();
			/* WAS INSTEAD:
			alignments2intervals_intersection(false);
			nUnassignedAlignments=alignments2intervals_updateAlignmentPointers(true);
			*/
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Factorize> periodic substring intervals after discardRareIntervals Z ():");
	for (int x=0; x<=PeriodicSubstrings.lastInterval; x++) System.err.println(PeriodicSubstrings.intervals[x]);
	IO.printErr("long-period periodic substring intervals after discardRareIntervals Z ():");
	for (int x=0; x<=PeriodicSubstrings.lastLongPeriodInterval; x++) IO.printErr(PeriodicSubstrings.longPeriodIntervals[x]);
	IO.printErr("Factorize> alignment intervals after discardRareIntervals Z ():");
	for (i=0; i<=AlignmentIntervals.lastInterval; i++) IO.printErr(AlignmentIntervals.intervals[i]);
	IO.printErr("Factorize> Alignments after discardRareIntervals Z ():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) IO.printErr(ReadA.sortedAlignments[x].toStringBoundaries()+" ==> "+ReadA.sortedAlignments[x].toStringPointers());
}
			
probe(14);			
			// Sorting all intervals, assigning the final IDs.
			if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION) {
				PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
				PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
				if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
			}
			for (j=0; j<=PeriodicSubstrings.lastInterval; j++) PeriodicSubstrings.intervals[j].id=j;		
			PeriodicSubstrings.mergeIntervals_copyIDs();
			if (DenseSubstring.order!=DenseSubstring.STARTA) {
				DenseSubstring.order=DenseSubstring.STARTA;
				if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
			}
			for (j=0; j<=DenseSubstrings.lastSubstring; j++) DenseSubstrings.substrings[j].id=j;
			if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
				AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
				if (AlignmentIntervals.lastInterval>0) Arrays.sort(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1);
			}
			for (j=0; j<=AlignmentIntervals.lastInterval; j++) AlignmentIntervals.intervals[j].id=j;

			// Marking intervals that are approximately contained inside another
			markContainedIntervals();
			
			if (IO.CONSISTENCY_CHECKS) {
				// At this stage an alignment can still be assigned both to an alignment
				// interval and to a dense or periodic interval.
				ReadA.checkAlignmentIntervals();
			}
			assignFirstLastMaximal();
			averageCoverageAssigned();
//			DenseSubstrings.fixBoundariesOfSubstringType(tmpArray);

if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("ALIGNMENTS IMPLIED BY THE DENSE SUBSTRING OF SUBSTRING TYPE:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) {
		DenseSubstring s = ReadA.sortedAlignments[x].inDenseSubstring;
		System.err.println(ReadA.sortedAlignments[x].toStringBoundaries());
		System.err.println(ReadA.sortedAlignments[x].toStringPointers());
		System.err.println();
	}
	IO.printErr("Factorize> dense substrings after fixBoundariesOfSubstringType():");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	System.err.println("ALIGNMENTS NOT IMPLIED BY ANYTHING:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) {
		if (!ReadA.sortedAlignments[x].isImplied()) {
			System.err.println(ReadA.sortedAlignments[x].toStringBoundaries());
			System.err.println(ReadA.sortedAlignments[x].toStringPointers());
			System.err.println();
		}
	}
}
probe(15);
			// Printing
			for (j=0; j<=PeriodicSubstrings.lastInterval; j++) {
				periodicSubstringsFile.write(ReadA.id+",");
				PeriodicSubstrings.intervals[j].writePeriodicSubstringsFile(periodicSubstringsFile);
			}
			for (j=0; j<=DenseSubstrings.lastSubstring; j++) {
				denseSubstringsFile.write(ReadA.id+",");
				DenseSubstrings.substrings[j].writeDenseSubstringsFile(denseSubstringsFile);
			}
			for (j=0; j<=AlignmentIntervals.lastInterval; j++) {
				alignmentsFile.write(ReadA.id+",");
				AlignmentIntervals.intervals[j].writeAlignmentsFile(alignmentsFile);
			}			
			writeIntervalConnectionFile();
			assemblyStats[12]+=ReadA.readLength;
			assemblyStatistics(assemblyStats);
			if (i%10==0) {
				// Flushing the buffers at regular intervals. 
				// WARNING: This is just for inspection, and should be removed or made
				// less frequent in the final version.
				factorsFile.flush();
				intervalsDenseFile.flush();
				intervalsPeriodicFile.flush();
				intervalsAlignmentFile.flush();
				denseSubstringsFile.flush();
				periodicSubstringsFile.flush();
				alignmentsFile.flush();
				intervalConnectionFile.flush();
			}
			// Next iteration
			PeriodicSubstrings.undo_mergeIntervals();
			i=ReadA.lastAlignment+1;			
		}
		
		factorizationTime=System.currentTimeMillis()-time;
		System.err.println("Done in "+(factorizationTime/1000)+"s ("+Alignments.nAlignments+" alignments, "+IO.format(((double)factorizationTime)/alignmentCounter)+"ms per alignment)");
		printAssemblyStatistics(assemblyStats,true);
		deallocateStats();
		factorsFile.close();
		intervalsDenseFile.close();
		intervalsPeriodicFile.close();
		intervalsAlignmentFile.close();
		denseSubstringsFile.close();
		periodicSubstringsFile.close();
		alignmentsFile.close();
		intervalConnectionFile.close();
	}

	
	/**
	 * Removes all dense substrings and alignment intervals such that: (1) the number of 
	 * maximal alignments currently assigned to them is small; or (2) a large fraction of 
	 * the interval has low quality. If a dense substring or alignment interval has no 
	 * left- or right-context of high quality, the filter applies to the average coverage.
	 * The procedure discards periodic substring intervals as well.
	 *
	 * Remark: the average coverage of a dense substring, computed using all alignments, 
	 * might be affected by the coverage of contained repeats, and the coverage computed 
	 * using just the alignments assigned to the substring varies widely along the 
	 * substring in practice, so taking its average is questionable. Finally, the average 
	 * coverage is not a correct measure with e.g. prefix substrings, in which the 
	 * coverage histogram is a line in theory.
	 *
	 * Remark: for the reasons above, it is not useful in practice to assign an average
	 * coverage, or a coverage range, to every interval, and to use such values or ranges 
	 * to decide whether to keep or discard an edge of the interval graph. The idea of
	 * keeping just repeats with high frequency, and to assemble repeat fragments with 
	 * similar frequency, appears already in \cite{chu2016repdenovo}.
	 *
	 * Remark: maximal alignments are not a good proxy for the number of occurrences of a
	 * repeat in the genome, since e.g. the B-reads in two distinct maximal alignments 
	 * might be sampled from the same region of the genome.
	 *
	 * Remark: the average coverage check is applied also to dense substrings, since, e.g. 
	 * the suffixes of a suffix substring might be very far from one another, and the 
	 * instance of the suffix substring in the current read might be very short and 
	 * sorrounded by low-quality regions. In this case, the number of maximal alignments
	 * would be low.
	 */
	private static final void discardRareIntervals() {
		final int MIN_MAXIMAL_ALIGNMENTS = (IO.minOccurrencesInGenome-1)*IO.coverage;
		final int MIN_MAXIMAL_ALIGNMENTS_BASELINE = 4;  // Arbitrary
		final int MIN_COVERAGE = IO.minOccurrencesInGenome*IO.coverage-1;
		final int MIN_COVERAGE_WINDOW = Alignments.minAlignmentLength;
		int i, j;
		int firstJForNextI, substringsDiscarded, alignmentsDiscarded, shortPeriodDiscarded, longPeriodDiscarded, nAlignmentPointers, length, minLength;
		DenseSubstring substring, minSubstring;
		AlignmentInterval aInterval;
		PeriodicSubstringInterval pInterval, minPInterval;
	
		computeMaximalAlignments(true,true,true,ReadA.lastSortedAlignment);
	
		// Marking discarded intervals
		substringsDiscarded=0; alignmentsDiscarded=0; shortPeriodDiscarded=0; longPeriodDiscarded=0;
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			substring=DenseSubstrings.substrings[i];
			if ( ( substring.prefixReplication && !substring.suffixReplication && !substring.substringReplication && !substring.singleDeletionReplication && 
			       Reads.isRightMaximal(substring.endA,ReadA.id,true)
				 ) ||
			     ( !substring.prefixReplication && substring.suffixReplication && !substring.substringReplication && !substring.singleDeletionReplication && 
			   	   Reads.isLeftMaximal(substring.startA,ReadA.id,true)
			   	 ) ||
			     ( ((substring.prefixReplication && substring.suffixReplication) || substring.singleDeletionReplication || substring.substringReplication) && 
			   	   (Reads.isLeftMaximal(substring.startA,ReadA.id,true) || Reads.isRightMaximal(substring.endA,ReadA.id,true))
			   	 )
			   ) {
	   			if (substring.nMaximalAlignments<Math.max(MIN_MAXIMAL_ALIGNMENTS,MIN_MAXIMAL_ALIGNMENTS_BASELINE)) {
	   				substring.discarded=true;
	   				substringsDiscarded++;
					continue;
	   			}
	   			else substring.discarded=false;
			}
			else {
				if ( substring.nMaximalAlignments<Math.max(MIN_MAXIMAL_ALIGNMENTS,MIN_MAXIMAL_ALIGNMENTS_BASELINE) ||
					 ReadA.getAverageCoverage(substring.startA,substring.endA)<MIN_COVERAGE || 
					 ReadA.hasSubstringWithLowCoverage(substring.startA,substring.endA,MIN_COVERAGE,MIN_COVERAGE_WINDOW)
				   ) {
   					substring.discarded=true;
   					substringsDiscarded++;
					continue;
				}
				else substring.discarded=false;
			}
		}
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			aInterval=AlignmentIntervals.intervals[i];
			if (Reads.isLeftMaximal(aInterval.firstPosition,ReadA.id,true) || Reads.isRightMaximal(aInterval.lastPosition,ReadA.id,true)) {
				if (aInterval.nMaximalAlignmentsLeft<MIN_MAXIMAL_ALIGNMENTS && aInterval.nMaximalAlignmentsRight<MIN_MAXIMAL_ALIGNMENTS && aInterval.nMaximalAlignments()<MIN_MAXIMAL_ALIGNMENTS) {
					aInterval.discarded=true;
					alignmentsDiscarded++;
					continue;
				}
				else aInterval.discarded=false;
			}
			else {
				if ( aInterval.nMergedIntervals<MIN_COVERAGE ||
					 ReadA.getAverageCoverage(aInterval.firstPosition,aInterval.lastPosition)<MIN_COVERAGE || 
					 ReadA.hasSubstringWithLowCoverage(aInterval.firstPosition,aInterval.lastPosition,MIN_COVERAGE,MIN_COVERAGE_WINDOW)
				   ) {
   					aInterval.discarded=true;
   					alignmentsDiscarded++;
					continue;
				}
				else aInterval.discarded=false;
			}
		}		
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			pInterval=PeriodicSubstrings.intervals[i];
			if (pInterval.nImpliedAlignments<Math.max(MIN_MAXIMAL_ALIGNMENTS,MIN_MAXIMAL_ALIGNMENTS_BASELINE)) {
				pInterval.discarded=true;
				if (pInterval.hasLongPeriod) longPeriodDiscarded++;
				else shortPeriodDiscarded++;
				continue;
			}
			else pInterval.discarded=false;
		}
		for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) {
			pInterval=PeriodicSubstrings.longPeriodIntervals[i];
			if (pInterval.nImpliedAlignments<Math.max(MIN_MAXIMAL_ALIGNMENTS,MIN_MAXIMAL_ALIGNMENTS_BASELINE)) {
				pInterval.discarded=true;
				continue;
			}
			else pInterval.discarded=false;
		}
		discardRareIntervals_lowQualitySubstrings(tmpArray);
		alignmentsDiscarded+=tmpArray[0];
		substringsDiscarded+=tmpArray[1];
		shortPeriodDiscarded+=tmpArray[2];
		longPeriodDiscarded+=tmpArray[3];
	
		// Disconnecting alignments from their discarded implying intervals
		nAlignmentPointers=0;
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (substring!=null && substring.discarded) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				nAlignmentPointers++;
			}
			substring=ReadA.sortedAlignments[i].inDenseSubstring;
			if (substring!=null && substring.discarded) {
				ReadA.sortedAlignments[i].inDenseSubstring=null;
				nAlignmentPointers++;
			}
			aInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (aInterval!=null && aInterval.discarded) {
				ReadA.sortedAlignments[i].mergedToInterval=null;
				nAlignmentPointers++;
			}
			pInterval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if (pInterval!=null && pInterval.discarded) {
				ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				nAlignmentPointers++;
			}
		}
	
		// Removing discarded intervals
		if (substringsDiscarded>0) {
			j=-1;
			for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
				if (DenseSubstrings.substrings[i].discarded) continue;
				j++;
				if (j==i) continue;
				substring=DenseSubstrings.substrings[j];
				DenseSubstrings.substrings[j]=DenseSubstrings.substrings[i];
				DenseSubstrings.substrings[i]=substring;
			}
			DenseSubstrings.lastSubstring=j;
		}
		if (alignmentsDiscarded>0) {
			j=-1;
			for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
				if (AlignmentIntervals.intervals[i].discarded) continue;
				j++;
				if (j==i) continue;
				aInterval=AlignmentIntervals.intervals[j];
				AlignmentIntervals.intervals[j]=AlignmentIntervals.intervals[i];
				AlignmentIntervals.intervals[i]=aInterval;
			}
			AlignmentIntervals.lastInterval=j;
		}
		if (shortPeriodDiscarded+longPeriodDiscarded>0) {
			j=-1;
			for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
				if (PeriodicSubstrings.intervals[i].discarded) continue;
				j++;
				if (j==i) continue;
				pInterval=PeriodicSubstrings.intervals[j];
				PeriodicSubstrings.intervals[j]=PeriodicSubstrings.intervals[i];
				PeriodicSubstrings.intervals[i]=pInterval;
			}
			for (i=j+1; i<=PeriodicSubstrings.lastInterval; i++) PeriodicSubstrings.intervals[i].hasLongPeriod=false;
			PeriodicSubstrings.lastInterval=j;
		}
		if (longPeriodDiscarded>0) {
			j=-1;
			for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) {
				if (PeriodicSubstrings.longPeriodIntervals[i].discarded) continue;
				j++;
				if (j==i) continue;
				pInterval=PeriodicSubstrings.longPeriodIntervals[j];
				PeriodicSubstrings.longPeriodIntervals[j]=PeriodicSubstrings.longPeriodIntervals[i];
				PeriodicSubstrings.longPeriodIntervals[i]=pInterval;
			}
			PeriodicSubstrings.lastLongPeriodInterval=j;
		}
	
		// Removing pointers from intervals to other discared intervals
		// Alignment -> Dense substring
		if (substringsDiscarded>0) {
			i=0; j=0; firstJForNextI=-1; 
			minSubstring=null; minLength=Math.POSITIVE_INFINITY;
			while (i<=AlignmentIntervals.lastInterval) {
				substring=AlignmentIntervals.intervals[i].inDenseSubstringOfSubstringType;			
				if (substring==null || !substring.discarded) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					minSubstring=null; minLength=Math.POSITIVE_INFINITY;
					continue;
				}
				if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>=AlignmentIntervals.intervals[i].lastPosition) {
					AlignmentIntervals.intervals[i].inDenseSubstringOfSubstringType=minSubstring;
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					minSubstring=null; minLength=Math.POSITIVE_INFINITY;
					continue;
				}
				substring=DenseSubstrings.substrings[j];
				if (substring.endA<AlignmentIntervals.intervals[i].firstPosition) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && substring.endA>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
				if (substring.substringReplication && Intervals.isApproximatelyContained(AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,substring.startA,substring.endA)) {
					length=substring.endA-substring.startA+1;
					if (length<minLength) {
						minLength=length;
						minSubstring=substring;
					}
				}
				j++;
			}
		}
		// Alignment -> Periodic
		if (shortPeriodDiscarded+longPeriodDiscarded>0) {
			i=0; j=0; firstJForNextI=-1; 
			minPInterval=null; minLength=Math.POSITIVE_INFINITY;
			while (i<=AlignmentIntervals.lastInterval) {
				pInterval=AlignmentIntervals.intervals[i].inPeriodicSubstringInterval;			
				if (pInterval==null || !pInterval.discarded) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					minPInterval=null; minLength=Math.POSITIVE_INFINITY;
					continue;
				}
				if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=AlignmentIntervals.intervals[i].lastPosition) {
					AlignmentIntervals.intervals[i].inPeriodicSubstringInterval=minPInterval;
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					minPInterval=null; minLength=Math.POSITIVE_INFINITY;
					continue;
				}
				pInterval=PeriodicSubstrings.intervals[j];
				if (pInterval.lastPosition<AlignmentIntervals.intervals[i].firstPosition) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && pInterval.lastPosition>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
				if (Intervals.isApproximatelyContained(AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,pInterval.firstPosition,pInterval.lastPosition)) {
					length=pInterval.length();
					if (length<minLength) {
						minLength=length;
						minPInterval=pInterval;
					}
				}
				j++;
			}
		}
		// Dense -> Periodic
		if (shortPeriodDiscarded+longPeriodDiscarded>0) {
			i=0; j=0; firstJForNextI=-1; 
			minPInterval=null; minLength=Math.POSITIVE_INFINITY;
			while (i<=DenseSubstrings.lastSubstring) {
				pInterval=DenseSubstrings.substrings[i].periodicSubstringInterval;			
				if (pInterval==null || !pInterval.discarded) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					minPInterval=null; minLength=Math.POSITIVE_INFINITY;
					continue;
				}
				if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=DenseSubstrings.substrings[i].endA) {
					DenseSubstrings.substrings[i].periodicSubstringInterval=minPInterval;
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					minPInterval=null; minLength=Math.POSITIVE_INFINITY;
					continue;
				}
				pInterval=PeriodicSubstrings.intervals[j];
				if (pInterval.lastPosition<DenseSubstrings.substrings[i].startA) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && pInterval.lastPosition>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
				if (Intervals.isApproximatelyContained(DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,pInterval.firstPosition,pInterval.lastPosition)) {
					length=pInterval.length();
					if (length<minLength) {
						minLength=length;
						minPInterval=pInterval;
					}
				}
				j++;
			}
		}
		// Removing pointers from periodic substrings to discared intervals
		if (shortPeriodDiscarded+longPeriodDiscarded>0) {
			for (i=0; i<=PeriodicSubstrings.lastSubstring; i++) {
				pInterval=PeriodicSubstrings.substrings[i].mergedToInterval;
				if (pInterval!=null && pInterval.discarded) PeriodicSubstrings.substrings[i].mergedToInterval=null;
			}
			for (i=0; i<=PeriodicSubstrings.lastLongPeriodSubstring; i++) {
				pInterval=PeriodicSubstrings.longPeriodSubstrings[i].mergedToInterval;
				if (pInterval!=null && pInterval.discarded) PeriodicSubstrings.longPeriodSubstrings[i].mergedToInterval=null;
			}
		}
	
		if (IO.SHOW_STD_ERR_PRIME) {
			System.err.println("discardRareIntervals> Discarded "+substringsDiscarded+" rare dense substrings, "+alignmentsDiscarded+" rare alignment intervals, "+shortPeriodDiscarded+" short-period intervals, "+longPeriodDiscarded+" long-period intervals.");
			System.err.println("discardRareIntervals> Removed "+nAlignmentPointers+" pointers from alignments");			
		}
	}
	
	
	/**
	 * Marks as discarded every interval that: (1) contains a low-quality substring; (2)
	 * does not contain any other interval without a low-quality substring.
	 *
	 * Remark: the procedure does not consider intervals that are already marked as 
	 * discarded.
	 *
	 * @param out 0=number of alignment intervals discarded by the procedure; 1=number of
	 * dense substrings; 2=short-period intervals; 3=long-period intervals.
	 */
	private static final void discardRareIntervals_lowQualitySubstrings(int[] out) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int PREFIX_SUFFIX_DISTANCE = 1;  // In quality spacing units. Arbitrary.
		final double LOW_QUALITY_FRACTION = 0.25;  // Arbitrary
		final int LOW_QUALITY_LENGTH = Alignments.minAlignmentLength/Reads.QUALITY_SPACING;  // In quality spacing units. Arbitrary.
		boolean found;
		int i, j;
		int firstJForNextI, nDiscarded, nDiscardedPrime;
		DenseSubstring substring;
		AlignmentInterval aInterval;
		PeriodicSubstringInterval pInterval;
		
		// Marking intervals
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			substring=DenseSubstrings.substrings[i];
			substring.destinationAlignment=0;  // Used as temporary space
			if (substring.discarded) continue;
			Reads.readEnds2qualityEnds(ReadA.id,substring.startA,substring.endA,true,tmpArray);
			if (Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),tmpArray[0],tmpArray[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,Math.min(LOW_QUALITY_LENGTH,Math.round((substring.length()*LOW_QUALITY_FRACTION)/Reads.QUALITY_SPACING)),false,PREFIX_SUFFIX_DISTANCE)) {
				substring.destinationAlignment=1;
			}
		}
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			aInterval=AlignmentIntervals.intervals[i];
			aInterval.tmpInt1=0;  // Used as temporary space
			if (aInterval.discarded) continue;
			Reads.readEnds2qualityEnds(ReadA.id,aInterval.firstPosition,aInterval.lastPosition,true,tmpArray);
			if (Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),tmpArray[0],tmpArray[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,Math.min(LOW_QUALITY_LENGTH,Math.round((aInterval.length()*LOW_QUALITY_FRACTION)/Reads.QUALITY_SPACING)),false,PREFIX_SUFFIX_DISTANCE)) {
				aInterval.tmpInt1=1;
			}
		}
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			pInterval=PeriodicSubstrings.intervals[i];
			pInterval.surface=0;  // Used as temporary space
			if (pInterval.discarded) continue;
			Reads.readEnds2qualityEnds(ReadA.id,pInterval.firstPosition,pInterval.lastPosition,true,tmpArray);
			if (Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),tmpArray[0],tmpArray[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,Math.min(LOW_QUALITY_LENGTH,Math.round((pInterval.length()*LOW_QUALITY_FRACTION)/Reads.QUALITY_SPACING)),false,PREFIX_SUFFIX_DISTANCE)) {
				pInterval.surface=1;
			}
		}
		for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) {
			pInterval=PeriodicSubstrings.longPeriodIntervals[i];
			pInterval.surface=0;  // Used as temporary space
			if (pInterval.discarded) continue;
			Reads.readEnds2qualityEnds(ReadA.id,pInterval.firstPosition,pInterval.lastPosition,true,tmpArray);
			if (Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),tmpArray[0],tmpArray[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,Math.min(LOW_QUALITY_LENGTH,Math.round((pInterval.length()*LOW_QUALITY_FRACTION)/Reads.QUALITY_SPACING)),false,PREFIX_SUFFIX_DISTANCE)) {
				pInterval.surface=1;
			}
		}
		
		// Discarding dense substrings
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			DenseSubstrings.substrings[i].startAAdded=false;  // Used as temporary space
		}
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=DenseSubstrings.lastSubstring) {	
			if (DenseSubstrings.substrings[i].discarded || DenseSubstrings.substrings[i].destinationAlignment==0) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (j>AlignmentIntervals.lastInterval || AlignmentIntervals.intervals[j].firstPosition>=DenseSubstrings.substrings[i].endA-IDENTITY_THRESHOLD) {
				if (found) DenseSubstrings.substrings[i].startAAdded=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (AlignmentIntervals.intervals[j].lastPosition<=DenseSubstrings.substrings[i].startA) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && AlignmentIntervals.intervals[j].lastPosition>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
			if ( AlignmentIntervals.intervals[j].tmpInt1==0 && 
			     ( Intervals.isApproximatelyContained_lowQuality(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,ReadA.id) ||
				   Intervals.areApproximatelyIdentical_lowQuality(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,ReadA.id)
				 )
			   ) {
				found=true;
			}
			j++;
		}
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=DenseSubstrings.lastSubstring) {
			if (DenseSubstrings.substrings[i].discarded || DenseSubstrings.substrings[i].destinationAlignment==0) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=DenseSubstrings.substrings[i].endA-IDENTITY_THRESHOLD) {
				if (found) DenseSubstrings.substrings[i].startAAdded=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (PeriodicSubstrings.intervals[j].lastPosition<=DenseSubstrings.substrings[i].startA) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && PeriodicSubstrings.intervals[j].lastPosition>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
			if ( PeriodicSubstrings.intervals[j].surface==0 && 
			     ( Intervals.isApproximatelyContained_lowQuality(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,ReadA.id) ||
 				   Intervals.areApproximatelyIdentical_lowQuality(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,ReadA.id)
				 )
			   ) {
				found=true;
			}
			j++;
		}
		nDiscarded=0;
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			if (!DenseSubstrings.substrings[i].discarded && DenseSubstrings.substrings[i].destinationAlignment==1 && !DenseSubstrings.substrings[i].startAAdded) {
				DenseSubstrings.substrings[i].discarded=true;
				nDiscarded++;
			}
		}
		out[1]=nDiscarded;
		
		// Discarding alignment intervals
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			AlignmentIntervals.intervals[i].flag1=false;  // Used as temporary space
		}
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=AlignmentIntervals.lastInterval) {
			if (AlignmentIntervals.intervals[i].discarded || AlignmentIntervals.intervals[i].tmpInt1==0) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>=AlignmentIntervals.intervals[i].lastPosition-IDENTITY_THRESHOLD) {
				if (found) AlignmentIntervals.intervals[i].flag1=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (DenseSubstrings.substrings[j].endA<=AlignmentIntervals.intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && DenseSubstrings.substrings[j].endA>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
			if ( DenseSubstrings.substrings[j].destinationAlignment==0 && 
			     ( Intervals.isApproximatelyContained_lowQuality(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,ReadA.id) ||
				   Intervals.areApproximatelyIdentical_lowQuality(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,ReadA.id)
				 )
			   ) {
				found=true;
			}
			j++;
		}
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=AlignmentIntervals.lastInterval) {
			if (AlignmentIntervals.intervals[i].discarded || AlignmentIntervals.intervals[i].tmpInt1==0) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=AlignmentIntervals.intervals[i].lastPosition-IDENTITY_THRESHOLD) {
				if (found) AlignmentIntervals.intervals[i].flag1=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (PeriodicSubstrings.intervals[j].lastPosition<=AlignmentIntervals.intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && PeriodicSubstrings.intervals[j].lastPosition>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
			if ( PeriodicSubstrings.intervals[j].surface==0 && 
			     ( Intervals.isApproximatelyContained_lowQuality(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,ReadA.id) ||
				   Intervals.areApproximatelyIdentical_lowQuality(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,ReadA.id)
				 )
			   ) {
				found=true;
			}
			j++;
		}
		nDiscarded=0;
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			if ( !AlignmentIntervals.intervals[i].discarded && AlignmentIntervals.intervals[i].tmpInt1==1 && !AlignmentIntervals.intervals[i].flag1 && 
			     !DenseSubstrings.isPrefixOfPrefixSubstring(AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition)
			   ) {
				// Discarding an alignment interval that is a prefix of a prefix
				// substring is likely to merge all its alignments back to such prefix
				// substring, thus making the creation of the interval pointless.
				AlignmentIntervals.intervals[i].discarded=true;
				nDiscarded++;
			}
		}
		out[0]=nDiscarded;
		
		// Discarding periodic intervals: $intervals$.
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) PeriodicSubstrings.intervals[i].flag1=false;
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=PeriodicSubstrings.lastInterval) {
			if (PeriodicSubstrings.intervals[i].discarded || PeriodicSubstrings.intervals[i].surface==0) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>=PeriodicSubstrings.intervals[i].lastPosition-IDENTITY_THRESHOLD) {
				if (found) PeriodicSubstrings.intervals[i].flag1=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (DenseSubstrings.substrings[j].endA<=PeriodicSubstrings.intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<PeriodicSubstrings.lastInterval && DenseSubstrings.substrings[j].endA>=PeriodicSubstrings.intervals[i+1].firstPosition) firstJForNextI=j;
			if ( DenseSubstrings.substrings[j].destinationAlignment==0 && 
			     ( Intervals.isApproximatelyContained_lowQuality(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id) ||
				   Intervals.areApproximatelyIdentical_lowQuality(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id)
				 )
			   ) {
				found=true;
			}
			j++;
		}
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=PeriodicSubstrings.lastInterval) {
			if (PeriodicSubstrings.intervals[i].discarded || PeriodicSubstrings.intervals[i].surface==0) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (j>AlignmentIntervals.lastInterval || AlignmentIntervals.intervals[j].firstPosition>=PeriodicSubstrings.intervals[i].lastPosition-IDENTITY_THRESHOLD) {
				if (found) PeriodicSubstrings.intervals[i].flag1=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (AlignmentIntervals.intervals[j].lastPosition<=PeriodicSubstrings.intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<PeriodicSubstrings.lastInterval && AlignmentIntervals.intervals[j].lastPosition>=PeriodicSubstrings.intervals[i+1].firstPosition) firstJForNextI=j;
			if ( AlignmentIntervals.intervals[j].tmpInt1==0 && 
			     ( Intervals.isApproximatelyContained_lowQuality(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id) ||
				   Intervals.areApproximatelyIdentical_lowQuality(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id)
				 )
			   ) {
				found=true;
			}
			j++;
		}
		nDiscarded=0; nDiscardedPrime=0;
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			if (!PeriodicSubstrings.intervals[i].discarded && PeriodicSubstrings.intervals[i].surface==1 && !PeriodicSubstrings.intervals[i].flag1) {
				PeriodicSubstrings.intervals[i].discarded=true;
				if (PeriodicSubstrings.intervals[i].hasLongPeriod) nDiscardedPrime++;
				else nDiscarded++;
			}
		}
		out[2]=nDiscarded; out[3]=nDiscardedPrime;
		
		// Discarding periodic intervals: $longPeriodIntervals$.
		for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) PeriodicSubstrings.longPeriodIntervals[i].flag1=false;
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=PeriodicSubstrings.lastLongPeriodInterval) {
			if (PeriodicSubstrings.longPeriodIntervals[i].discarded || PeriodicSubstrings.longPeriodIntervals[i].surface==0) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>=PeriodicSubstrings.longPeriodIntervals[i].lastPosition-IDENTITY_THRESHOLD) {
				if (found) PeriodicSubstrings.longPeriodIntervals[i].flag1=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (DenseSubstrings.substrings[j].endA<=PeriodicSubstrings.longPeriodIntervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<PeriodicSubstrings.lastLongPeriodInterval && DenseSubstrings.substrings[j].endA>=PeriodicSubstrings.longPeriodIntervals[i+1].firstPosition) firstJForNextI=j;
			if ( DenseSubstrings.substrings[j].destinationAlignment==0 && 
			     ( Intervals.isApproximatelyContained_lowQuality(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,PeriodicSubstrings.longPeriodIntervals[i].firstPosition,PeriodicSubstrings.longPeriodIntervals[i].lastPosition,ReadA.id) ||
				   Intervals.areApproximatelyIdentical_lowQuality(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,PeriodicSubstrings.longPeriodIntervals[i].firstPosition,PeriodicSubstrings.longPeriodIntervals[i].lastPosition,ReadA.id)
				 )
			   ) {
				found=true;
			}
			j++;
		}
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=PeriodicSubstrings.lastLongPeriodInterval) {
			if (PeriodicSubstrings.longPeriodIntervals[i].discarded || PeriodicSubstrings.longPeriodIntervals[i].surface==0) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (j>AlignmentIntervals.lastInterval || AlignmentIntervals.intervals[j].firstPosition>=PeriodicSubstrings.longPeriodIntervals[i].lastPosition-IDENTITY_THRESHOLD) {
				if (found) PeriodicSubstrings.longPeriodIntervals[i].flag1=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; found=false;
				continue;
			}
			if (AlignmentIntervals.intervals[j].lastPosition<=PeriodicSubstrings.longPeriodIntervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<PeriodicSubstrings.lastLongPeriodInterval && AlignmentIntervals.intervals[j].lastPosition>=PeriodicSubstrings.longPeriodIntervals[i+1].firstPosition) firstJForNextI=j;
			if (AlignmentIntervals.intervals[j].tmpInt1==0 && Intervals.isApproximatelyContained_lowQuality(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,PeriodicSubstrings.longPeriodIntervals[i].firstPosition,PeriodicSubstrings.longPeriodIntervals[i].lastPosition,ReadA.id)) {
				found=true;
			}
			j++;
		}
		for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) {
			if (!PeriodicSubstrings.longPeriodIntervals[i].discarded && PeriodicSubstrings.longPeriodIntervals[i].surface==1 && !PeriodicSubstrings.longPeriodIntervals[i].flag1) {
				PeriodicSubstrings.longPeriodIntervals[i].discarded=true;
			}
		}
	}

	
	/**
	 * Sets the $nMaximalAlignments$ values of dense substrings and alignment intervals,
	 * and the $nImpliedAlignments$ and $longestAlignment$ values of periodic intervals.
	 */
	public static final void computeMaximalAlignments(boolean denseSubstrings, boolean alignmentIntervals, boolean periodicIntervals, int lastAlignment) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		boolean isLeftMaximal, isRightMaximal, found;
		int i, j, k;
		int length;
		DenseSubstring substring;
		AlignmentInterval aInterval;
		PeriodicSubstringInterval pInterval;
		
		if (denseSubstrings) {
			for (i=0; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].nMaximalAlignments=0;
		}
		if (alignmentIntervals) {
			for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
				AlignmentIntervals.intervals[i].nMaximalAlignmentsLeft=0;
				AlignmentIntervals.intervals[i].nMaximalAlignmentsRight=0;
				AlignmentIntervals.intervals[i].nMergedIntervals=0;
			}
		}
		if (periodicIntervals) {
			for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
				PeriodicSubstrings.intervals[i].nImpliedAlignments=0;
				PeriodicSubstrings.intervals[i].mass=0;
				PeriodicSubstrings.intervals[i].longestAlignment=0;
			}
			for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) {
				PeriodicSubstrings.longPeriodIntervals[i].nImpliedAlignments=0;
				PeriodicSubstrings.longPeriodIntervals[i].mass=0;
				PeriodicSubstrings.longPeriodIntervals[i].longestAlignment=0;
			}
		}

		// Counting
		for (i=0; i<=lastAlignment; i++) {
			isLeftMaximal=ReadA.sortedAlignments[i].isLeftMaximalB==1;
			isRightMaximal=ReadA.sortedAlignments[i].isRightMaximalB==1;
			if (denseSubstrings) {
				substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
				if (substring!=null) {
					if ((substring.prefixReplication || substring.singleDeletionReplication) && isRightMaximal) substring.nMaximalAlignments++;
					else if ((substring.suffixReplication || substring.singleDeletionReplication) && isLeftMaximal) substring.nMaximalAlignments++;
				}
				substring=ReadA.sortedAlignments[i].inDenseSubstring;
				if (substring!=null && ReadA.sortedAlignments[i].mergedToInterval==null) {
					// The alignment does not contribute to the substring if it is
					// assigned to an alignment interval.
					if (isLeftMaximal || isRightMaximal) substring.nMaximalAlignments++;
				}
			}
			if (alignmentIntervals) {
				aInterval=ReadA.sortedAlignments[i].mergedToInterval;
				if (aInterval!=null) {
					if (isLeftMaximal) aInterval.nMaximalAlignmentsLeft++;
					if (isRightMaximal) aInterval.nMaximalAlignmentsRight++;
					aInterval.nMergedIntervals++;
				}
			}
			if (periodicIntervals) {
				pInterval=ReadA.sortedAlignments[i].periodicSubstringInterval;
				if (pInterval!=null) {
					pInterval.nImpliedAlignments++;
					if (isLeftMaximal) pInterval.mass++;
					if (isRightMaximal) pInterval.mass++;
					length=ReadA.sortedAlignments[i].getALength();
					if (length<pInterval.length()-IDENTITY_THRESHOLD) pInterval.longestAlignment=Math.max(pInterval.longestAlignment,length);
				}
			}
		}
		
		// Copying counts from the real long-period intervals to long-period clones.
		if (periodicIntervals) {
			for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
				if (!PeriodicSubstrings.intervals[i].hasLongPeriod) continue;
				j=Arrays.binarySearch(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1,PeriodicSubstrings.intervals[i]);
				if (j<0) {
					System.err.println("computeMaximalAlignments> ERROR: long-period interval not found: "+PeriodicSubstrings.intervals[i].hashCode()+" :: "+PeriodicSubstrings.intervals[i]);
					System.err.println("all intervals:");
					for (int x=0; x<=PeriodicSubstrings.lastInterval; x++) System.err.println(PeriodicSubstrings.intervals[x].hashCode()+" :: "+PeriodicSubstrings.intervals[x]);
					System.err.println("long-period intervals:");
					for (int x=0; x<=PeriodicSubstrings.lastLongPeriodInterval; x++) System.err.println(PeriodicSubstrings.longPeriodIntervals[x].hashCode()+" :: "+PeriodicSubstrings.longPeriodIntervals[x]);
					System.exit(1);
				}
				found=false;
				for (k=j; k>=0; k--) {
					if (PeriodicSubstrings.longPeriodIntervals[k].equals(PeriodicSubstrings.intervals[i])) {
						PeriodicSubstrings.intervals[i].nImpliedAlignments=PeriodicSubstrings.longPeriodIntervals[k].nImpliedAlignments;
						PeriodicSubstrings.intervals[i].mass=PeriodicSubstrings.longPeriodIntervals[k].mass;
						PeriodicSubstrings.intervals[i].longestAlignment=PeriodicSubstrings.longPeriodIntervals[k].longestAlignment;
						found=true;
						break;
					}
				}
				if (found) continue;
				for (k=j+1; k<=PeriodicSubstrings.lastLongPeriodInterval; k++) {
					if (PeriodicSubstrings.longPeriodIntervals[k].equals(PeriodicSubstrings.intervals[i])) {
						PeriodicSubstrings.intervals[i].nImpliedAlignments=PeriodicSubstrings.longPeriodIntervals[k].nImpliedAlignments;
						PeriodicSubstrings.intervals[i].mass=PeriodicSubstrings.longPeriodIntervals[k].mass;
						PeriodicSubstrings.intervals[i].longestAlignment=PeriodicSubstrings.longPeriodIntervals[k].longestAlignment;
						break;
					}
				}
			}
		}
	}


	/**
	 * For every alignment that is assigned to an interval, the procedure removes the 
	 * assignment if the alignment is not contained in the interval. This could happen,
	 * since interval boundaries are the result of multiple mergings and refinements that
	 * happen throughout the whole pipeline. Performing this correction is important for
	 * having a consistent output, but it makes finding bugs upstream harder. Thus, 
	 * alignments that contain an interval and flanking sequence, or alignments that 
	 * contain part of an interval and flanking sequence, are not used. Something similar 
	 * is described in \cite{han2010mite}.
	 *
	 * For every alignment that is not already assigned to an interval (or that is not 
	 * assigned to an interval any more) the procedure assigns the alignment to a shortest
	 * interval that: (1) approximately contains the readA part of the alignment, or is 
	 * approximately identical to it; (2) does not straddle or contain another interval 
	 * that contains the readA part of the alignment. This mitigates, but does not solve, 
	 * the following problem:
	 *
	 * Let $X=X_{1}X_{2}$, $Y=Y_{1}Y_{2}$ and $Z=X_{2}Y_{1}$ be simple repeats, and let:
	 *
	 * read 1 contain an occurrence of $XY$;
	 * read 2 contain an occurrence of $X$;
	 * read 3 contain an occurrence of $Z$;
	 * read 4 contain an occurrence of $Y$.
	 *
	 * Every alignment between read 2 and read 1 involves the whole $X$, every alignment 
	 * between read 4 and read 1 involves the whole $Y$, and every alignment between read 
	 * 3 and read 1 involves the whole $Z$. If multiple reads similar to read 2 and 4 are 
	 * present in the dataset, then the occurrence of $Z$ in read 3 should actually be 
	 * split into two components $X_2$ and $Y_1$. Every alignment of read 2 with read 3 
	 * should involve $X_2$ and every alignment of read 4 with read 3 should involve 
	 * $Y_1$. However, $Z$ is not split in read 1.
	 *
	 * Assume that read 2 contains just a (not necessarily proper) suffix of $X_2$, for
	 * example because of a long low-quality region. The alignment between read 2 and read
	 * 1 could then be assigned to $Z$ in read 1. If the occurrence of $X$ in read 2 
	 * belongs to the connected component of $X$, this puts $Z$ in the same connected 
	 * component as $X$. If this happens for $Y$ as well, the resulting connected 
	 * component does not correspond to a single repeat but to a set of repeats that occur
	 * consecutively in the genome.
	 *
	 * This procedure solves the problem above by not assigning the alignment between read
	 * 2 and read 1 to any interval in read 1. However, this is not enough to solve the 
	 * problem in the general case, not even in theory. For example, if $Z$ is a dense 
	 * substring of substring type which for some reason is not split into $X_2$ and $Y_1$
	 * in read 3, then any alignment between read 2 and read 3, and between read 4 and 
	 * read 3, would recreate the problem. In practice the occurrence of $Z$ in read 3 
	 * might not get split even if $Z$ is a simple repeat.
	 *
	 * Remark: recall that, contrary to non-weak dense substrings of substring type,
	 * weak dense substrings of substring type can contain other dense substrings. The 
	 * procedure allows alignments that are contained in a weak dense substring of 
	 * substring type to be reassigned to a shorter interval.
	 *
	 * Remark: the procedure tries to reassign an alignment that is assigned to an 
	 * alignment interval, if the alignment is contained in a dense substring of substring
	 * type, but the alignment interval is not contained in the same dense substring.
	 *
	 * Remark: the procedure tries also to reassign an alignment that is assigned to 
	 * a prefix interval, and ends close to the end, but starts far from the beginning, of
	 * the interval. Such alignment might be induced by e.g. a simple repeat that ends at 
	 * the end of the prefix interval and starts inside the prefix interval. The same 
	 * holds for an alignment that starts close to the start, but ends far from the end:
	 * it might be induced e.g. by a simple repeat at the prefix, which has been detected
	 * by the pipeline but not assigned to such interval.
	 *
	 * Remark: the procedure considers also the left- and right-maximality of alignments
	 * to avoid wrong reassignments. The proper way to reassign an alignment would be to 
	 * use all the details in the code for detecting each type of interval, but it is not
	 * practical to access such code from this level of abstraction.
	 *
	 * Remark: non-assigned alignments at the end of this procedure are alignments for 
	 * which it is ambiguous which interval originates them. They should be treated
	 * accordingly in the rest of the pipeline.
	 *
	 * Remark: the procedure sorts $ReadA.sortedAlignments$ by $startA$.
	 *
	 * @return the number of unassigned alignments at the end of the procedure.
	 */
	private static final int alignments2intervals() {
		final int THRESHOLD = IO.quantum;
		final int ASSIGNED_THRESHOLD = IO.quantum;  // Arbitrary
		final int PREFIX_SUFFIX_THRESHOLD = (THRESHOLD)<<1;  // Arbitrary
		final int MIN_EXTRA_LENGTH = IO.quantum<<1;  // Arbitrary
		boolean found, previousSortByID;
		int i, j, k;
		int firstJForNextI, id, length, alignmentStart, alignmentEnd, alignmentLength;
		double intersectionFraction;
		DenseSubstring tmpSubstring;
		PeriodicSubstringInterval tmpPInterval;
		AlignmentInterval tmpAInterval;				
				
		// Ensuring the necessary order in $ReadA.sortedAlignments$
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

		// Initializing alignments
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			// Removing inconsistent pointers
			alignmentStart=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			alignmentEnd=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			alignmentLength=alignmentEnd-alignmentStart+1;
			tmpPInterval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if ( tmpPInterval!=null && 
				 ((tmpPInterval.isContained&(1<<Constants.INTERVAL_PERIODIC))==0 && (tmpPInterval.isContained&(1<<(Constants.INTERVAL_PERIODIC+1)))==0) &&
				 !Intervals.isApproximatelyContained_lowQuality(alignmentStart,alignmentEnd,tmpPInterval.firstPosition,tmpPInterval.lastPosition,ReadA.id) &&
				 !Intervals.areApproximatelyIdentical_lowQuality(alignmentStart,alignmentEnd,tmpPInterval.firstPosition,tmpPInterval.lastPosition,ReadA.id) &&
				 !Intervals.contains_lowQuality(alignmentStart,alignmentEnd,tmpPInterval.firstPosition,tmpPInterval.lastPosition,ReadA.id) &&
				 alignmentLength-Intervals.intersectionLength(alignmentStart,alignmentEnd,tmpPInterval.firstPosition,tmpPInterval.lastPosition)>=MIN_EXTRA_LENGTH
			   ) ReadA.sortedAlignments[i].periodicSubstringInterval=null;
			tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if ( tmpSubstring!=null && 
				 ( ( !Intervals.isApproximatelyContained_lowQuality(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA,ReadA.id) &&
				     !Intervals.areApproximatelyIdentical_lowQuality(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA,ReadA.id) &&
				     !Intervals.contains_lowQuality(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA,ReadA.id) &&
			         alignmentLength-Intervals.intersectionLength(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA)>=MIN_EXTRA_LENGTH
				   ) ||
				   !tmpSubstring.canBeAssigned(ReadA.sortedAlignments[i],ASSIGNED_THRESHOLD)
				 )
			   ) ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
			tmpAInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if ( tmpAInterval!=null && 
				 ( ( !Intervals.isApproximatelyContained_lowQuality(alignmentStart,alignmentEnd,tmpAInterval.firstPosition,tmpAInterval.lastPosition,ReadA.id) &&
				     !Intervals.areApproximatelyIdentical_lowQuality(alignmentStart,alignmentEnd,tmpAInterval.firstPosition,tmpAInterval.lastPosition,ReadA.id) &&
				     !Intervals.contains_lowQuality(alignmentStart,alignmentEnd,tmpAInterval.firstPosition,tmpAInterval.lastPosition,ReadA.id) &&
				     alignmentLength-Intervals.intersectionLength(alignmentStart,alignmentEnd,tmpAInterval.firstPosition,tmpAInterval.lastPosition)>=MIN_EXTRA_LENGTH
				   ) ||
				   !tmpAInterval.canBeAssigned(ReadA.sortedAlignments[i],ASSIGNED_THRESHOLD)
				 )
			   ) ReadA.sortedAlignments[i].mergedToInterval=null;
			tmpSubstring=ReadA.sortedAlignments[i].inDenseSubstring;
			if ( tmpSubstring!=null && 
				 !Intervals.isApproximatelyContained_lowQuality(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA,ReadA.id) &&
				 !Intervals.areApproximatelyIdentical_lowQuality(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA,ReadA.id) &&
				 !Intervals.contains_lowQuality(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA,ReadA.id) &&
				 alignmentLength-Intervals.intersectionLength(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA)>=MIN_EXTRA_LENGTH
			   ) ReadA.sortedAlignments[i].inDenseSubstring=null;
			
			// Initializing variables
			ReadA.sortedAlignments[i].minIntervalLength=Math.POSITIVE_INFINITY;
			ReadA.sortedAlignments[i].minIntervalID=-1;
			if (ReadA.sortedAlignments[i].periodicSubstringInterval!=null) {
				tmpPInterval=ReadA.sortedAlignments[i].periodicSubstringInterval;
				if (Intervals.areApproximatelyIdentical(alignmentStart,alignmentEnd,tmpPInterval.firstPosition,tmpPInterval.lastPosition)) ReadA.sortedAlignments[i].minIntervalType=-3;
				else ReadA.sortedAlignments[i].minIntervalType=-2;
				continue;
			}
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null) {
				tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
				if ( tmpSubstring.prefixReplication && !tmpSubstring.suffixReplication && !tmpSubstring.singleDeletionReplication && 
				     Math.abs(tmpSubstring.startA,alignmentStart)>PREFIX_SUFFIX_THRESHOLD && alignmentEnd<tmpSubstring.endA-THRESHOLD
				   ) {
					ReadA.sortedAlignments[i].minIntervalType=-2;
					continue;
				}
				if ( tmpSubstring.suffixReplication && !tmpSubstring.prefixReplication && !tmpSubstring.singleDeletionReplication &&
				     Math.abs(tmpSubstring.endA,alignmentEnd)>PREFIX_SUFFIX_THRESHOLD && alignmentStart>tmpSubstring.startA+THRESHOLD
				   ) {
					ReadA.sortedAlignments[i].minIntervalType=-2;
					continue;
				}
				if ( ((tmpSubstring.prefixReplication && tmpSubstring.suffixReplication) || tmpSubstring.singleDeletionReplication) &&
				     (Math.abs(tmpSubstring.startA,alignmentStart)>PREFIX_SUFFIX_THRESHOLD || Math.abs(tmpSubstring.endA,alignmentEnd)>PREFIX_SUFFIX_THRESHOLD)
				   ) {
					ReadA.sortedAlignments[i].minIntervalType=-2;
					continue;
				}
				ReadA.sortedAlignments[i].minIntervalType=-3;
				continue;
			}
			if (ReadA.sortedAlignments[i].mergedToInterval!=null) {
				tmpSubstring=ReadA.sortedAlignments[i].inDenseSubstring;
				if ( tmpSubstring!=null && 
					 ( Intervals.isApproximatelyContained_lowQuality(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA,ReadA.id) || 
					   Intervals.areApproximatelyIdentical_lowQuality(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA,ReadA.id)
					 ) &&
				     !Intervals.isApproximatelyContained_lowQuality(tmpAInterval.firstPosition,tmpAInterval.lastPosition,tmpSubstring.startA,tmpSubstring.endA,ReadA.id)
				   ) {
					ReadA.sortedAlignments[i].minIntervalType=-2;
					continue;   	
				}
				ReadA.sortedAlignments[i].minIntervalType=-3;
				continue;
			}
			if (ReadA.sortedAlignments[i].inDenseSubstring!=null) {
				if (ReadA.sortedAlignments[i].inDenseSubstring.isWeak) {
					tmpSubstring=ReadA.sortedAlignments[i].inDenseSubstring;
					if (Intervals.areApproximatelyIdentical(alignmentStart,alignmentEnd,tmpSubstring.startA,tmpSubstring.endA)) ReadA.sortedAlignments[i].minIntervalType=-3;
					else ReadA.sortedAlignments[i].minIntervalType=-2;
					continue;
				}
				else {
					ReadA.sortedAlignments[i].minIntervalType=-3;
					continue;
				}
			}
			ReadA.sortedAlignments[i].minIntervalType=-1;			
		}
		
		// Periodic substring intervals
		if (PeriodicSubstrings.lastInterval>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {
				if (ReadA.sortedAlignments[i].minIntervalType==-3) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				id=ReadA.sortedAlignments[i].id;
				if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=Alignments.alignments[id][4]) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (PeriodicSubstrings.intervals[j].lastPosition<=Alignments.alignments[id][3]) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && PeriodicSubstrings.intervals[j].lastPosition>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=PeriodicSubstrings.intervals[j].lastPosition-PeriodicSubstrings.intervals[j].firstPosition+1;
				if ( ( ( Intervals.isApproximatelyContained_lowQuality(Alignments.alignments[id][3],Alignments.alignments[id][4],PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.id) ||
						 Intervals.areApproximatelyIdentical_lowQuality(Alignments.alignments[id][3],Alignments.alignments[id][4],PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.id)
					   ) &&
					   ( // Case 1: alignment already assigned to a long-period interval.
					     ( ( ReadA.sortedAlignments[i].periodicSubstringInterval!=null &&
						     ReadA.sortedAlignments[i].periodicSubstringInterval.hasLongPeriod &&
						     !PeriodicSubstrings.intervals[j].hasLongPeriod &&
							 Intervals.isApproximatelyContained(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.sortedAlignments[i].periodicSubstringInterval.firstPosition,ReadA.sortedAlignments[i].periodicSubstringInterval.lastPosition) &&
							 !Intervals.areApproximatelyIdentical(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.sortedAlignments[i].periodicSubstringInterval.firstPosition,ReadA.sortedAlignments[i].periodicSubstringInterval.lastPosition) &&
							 Math.abs(Alignments.alignments[id][3],ReadA.sortedAlignments[i].periodicSubstringInterval.firstPosition)>THRESHOLD &&
						     Math.abs(Alignments.alignments[id][4],ReadA.sortedAlignments[i].periodicSubstringInterval.lastPosition)>THRESHOLD
						   ) && 
						   // We allow for ambiguity with non-periodic intervals.
						   !straddlesOrContains_periodic(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],false) 
				         ) ||
						 // Case 2: periodic substring mark only.
						 ( ( ReadA.sortedAlignments[i].periodicSubstringInterval==null && 
						     ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null
						   ) &&
						   // We allow for ambiguity with non-periodic intervals.
						   !straddlesOrContains_periodic(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],false)
						 ) ||
						 // Case 3: no periodic mark on the alignment.
						 ( ( ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null && 
						     ReadA.sortedAlignments[i].periodicSubstringInterval==null 
						   ) &&
						   // No ambiguity allowed.
						   !( straddlesOrContains_periodic(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],false) ||
					  	      straddlesOrContains_dense(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],ReadA.sortedAlignments[i],ASSIGNED_THRESHOLD) || 
					          straddlesOrContains_alignment(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],ReadA.sortedAlignments[i],ASSIGNED_THRESHOLD)
					       )		  
					     ) ||
						 // Case 4: alignment already assigned to a short-period interval.
						 // We reassign it to a shorter periodic interval, of any type,
						 // if one end coincides with the alignment end.
					     ( ( ReadA.sortedAlignments[i].periodicSubstringInterval!=null &&
						     !ReadA.sortedAlignments[i].periodicSubstringInterval.hasLongPeriod &&
							 Intervals.isApproximatelyContained(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.sortedAlignments[i].periodicSubstringInterval.firstPosition,ReadA.sortedAlignments[i].periodicSubstringInterval.lastPosition) &&
							 ( Math.abs(Alignments.alignments[id][3],PeriodicSubstrings.intervals[j].firstPosition)<=THRESHOLD ||
							   Math.abs(Alignments.alignments[id][4],PeriodicSubstrings.intervals[j].lastPosition)<=THRESHOLD
							 )
						   ) && 
						   // We allow for ambiguity with non-periodic intervals.
						   !straddlesOrContains_periodic(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],false) 
				         )
					   ) 
					 ) &&
					 // Checks on whether the alignment can be assigned are not necessary.
					 length<ReadA.sortedAlignments[i].minIntervalLength
				   ) {
					ReadA.sortedAlignments[i].minIntervalLength=length;
					ReadA.sortedAlignments[i].minIntervalType=Alignment.INTERVAL_PERIODIC;
					ReadA.sortedAlignments[i].minIntervalID=j;					
				}
				j++;
			}
		}
		
		// Dense substrings
		if (DenseSubstrings.lastSubstring>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {				
				if (ReadA.sortedAlignments[i].minIntervalType==-3) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}				
				id=ReadA.sortedAlignments[i].id;
				if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>=Alignments.alignments[id][4]) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (DenseSubstrings.substrings[j].endA<=Alignments.alignments[id][3]) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && DenseSubstrings.substrings[j].endA>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=DenseSubstrings.substrings[j].endA-DenseSubstrings.substrings[j].startA+1;
				if ( ( Intervals.areApproximatelyIdentical_lowQuality(Alignments.alignments[id][3],Alignments.alignments[id][4],DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,ReadA.id) 
					   ||
					   ( Intervals.isApproximatelyContained_lowQuality(Alignments.alignments[id][3],Alignments.alignments[id][4],DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,ReadA.id) &&
						 ReadA.sortedAlignments[i].periodicSubstringInterval==null && ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null &&
						 ( (ReadA.sortedAlignments[i].mergedToInterval!=null && ReadA.sortedAlignments[i].inDenseSubstring!=null && DenseSubstrings.substrings[j].substringReplication)
						   ||
					       !( straddlesOrContains_periodic(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,Alignments.alignments[id][3],Alignments.alignments[id][4],false) ||
	  					      straddlesOrContains_dense(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,Alignments.alignments[id][3],Alignments.alignments[id][4],ReadA.sortedAlignments[i],ASSIGNED_THRESHOLD) || 
	  				          straddlesOrContains_alignment(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,Alignments.alignments[id][3],Alignments.alignments[id][4],ReadA.sortedAlignments[i],ASSIGNED_THRESHOLD)
	  				        )
						 ) 
					   ) 
					 ) &&
					 DenseSubstrings.substrings[j].canBeAssigned(ReadA.sortedAlignments[i],THRESHOLD) &&
					 length<ReadA.sortedAlignments[i].minIntervalLength
				   ) {
					ReadA.sortedAlignments[i].minIntervalLength=length;
					ReadA.sortedAlignments[i].minIntervalType=Alignment.INTERVAL_DENSE;
					ReadA.sortedAlignments[i].minIntervalID=j;
				}
				j++;
			}
		}			
		
		// Alignment intervals
		if (AlignmentIntervals.lastInterval>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {
				if (ReadA.sortedAlignments[i].minIntervalType==-3) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				id=ReadA.sortedAlignments[i].id;
				if (j>AlignmentIntervals.lastInterval || AlignmentIntervals.intervals[j].firstPosition>=Alignments.alignments[id][4]) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (AlignmentIntervals.intervals[j].lastPosition<=Alignments.alignments[id][3]) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && AlignmentIntervals.intervals[j].lastPosition>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=AlignmentIntervals.intervals[j].lastPosition-AlignmentIntervals.intervals[j].firstPosition+1;
				if ( ( Intervals.areApproximatelyIdentical_lowQuality(Alignments.alignments[id][3],Alignments.alignments[id][4],AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,ReadA.id)
					   ||
					   ( Intervals.isApproximatelyContained_lowQuality(Alignments.alignments[id][3],Alignments.alignments[id][4],AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,ReadA.id) &&
						 ReadA.sortedAlignments[i].periodicSubstringInterval==null && ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null &&
						 ( (ReadA.sortedAlignments[i].mergedToInterval!=null && ReadA.sortedAlignments[i].inDenseSubstring!=null)
						   || 
						   !( straddlesOrContains_periodic(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],false) ||
		  				      straddlesOrContains_dense(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],ReadA.sortedAlignments[i],ASSIGNED_THRESHOLD) || 
		  				      straddlesOrContains_alignment(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,Alignments.alignments[id][3],Alignments.alignments[id][4],ReadA.sortedAlignments[i],ASSIGNED_THRESHOLD)
		  				    )
						 )
					   )
					 ) &&
					 AlignmentIntervals.intervals[j].canBeAssigned(ReadA.sortedAlignments[i],THRESHOLD) &&
					 length<ReadA.sortedAlignments[i].minIntervalLength
				   ) {
					ReadA.sortedAlignments[i].minIntervalLength=length;
					ReadA.sortedAlignments[i].minIntervalType=Alignment.INTERVAL_ALIGNMENT;
					ReadA.sortedAlignments[i].minIntervalID=j;
				}
				j++;
			}
		}
		
		// Weaker condition: assigning an unassigned alignment to its only intersection
		// with an interval, if any.
		alignments2intervals_intersection(true);
		
		// Updating pointers from alignments
		return alignments2intervals_updateAlignmentPointers(false);
	}
	
	
	/**
	 * Sets the $minIntervalType$ and $minIntervalID$ fields of an unassigned alignment to
	 * the type and ID of its only intersecting periodic substring interval, dense 
	 * substring, or alignment interval, if any.
	 *
	 * Remark: the procedure assumes all intervals and alignments to be sorted by first
	 * position.
	 *
	 * @param filter TRUE: the procedure does not initialize $minIntervalType$ and 
	 * $minIntervalID$, and it does not consider alignments with $minIntervalType==-3$ or 
	 * $minIntervalID!=-1$; FALSE: the procedure does not consider alignments that are 
	 * already implied by an interval, and it initializes $minIntervalType$ and 
	 * $minIntervalID$;
	 * @return the number of alignments assigned by the procedure.
	 */
	private static final int alignments2intervals_intersection(boolean filter) {
		final int THRESHOLD = IO.quantum;
		final double INTERSECTION_THRESHOLD = 0.7;  // Arbitrary
		int i, j;
		int firstJForNextI, id, length, out;
		
		// Storing in $minIntervalLength$ the number of intersecting intervals.
		// $minIntervalLength$ is negative if the alignment intersects with more
		// than one periodic interval.
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if (filter) {
				if (ReadA.sortedAlignments[i].minIntervalType!=-3 && ReadA.sortedAlignments[i].minIntervalID==-1) ReadA.sortedAlignments[i].minIntervalLength=0;
			}
			else {
				if (!ReadA.sortedAlignments[i].isImplied()) {
					ReadA.sortedAlignments[i].minIntervalLength=0;
					ReadA.sortedAlignments[i].minIntervalType=-1;
					ReadA.sortedAlignments[i].minIntervalID=-1;
				}
			}
		}
		
		// Counting intersections
		if (PeriodicSubstrings.lastInterval>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {
				if (filter?(ReadA.sortedAlignments[i].minIntervalType==-3||ReadA.sortedAlignments[i].minIntervalID!=-1):ReadA.sortedAlignments[i].isImplied()) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				id=ReadA.sortedAlignments[i].id;
				if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>Alignments.alignments[id][4]) {
					if (ReadA.sortedAlignments[i].minIntervalLength>1) ReadA.sortedAlignments[i].minIntervalLength=-ReadA.sortedAlignments[i].minIntervalLength; 
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (PeriodicSubstrings.intervals[j].lastPosition<Alignments.alignments[id][3]) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && PeriodicSubstrings.intervals[j].lastPosition>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=Intervals.intersectionLength(Alignments.alignments[id][3],Alignments.alignments[id][4],PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition);
				if (length>THRESHOLD) ReadA.sortedAlignments[i].minIntervalLength++;
				j++;
			}
		}
		if (DenseSubstrings.lastSubstring>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {				
				if (filter?(ReadA.sortedAlignments[i].minIntervalType==-3||ReadA.sortedAlignments[i].minIntervalID!=-1):ReadA.sortedAlignments[i].isImplied()) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}				
				id=ReadA.sortedAlignments[i].id;
				if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>Alignments.alignments[id][4]) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (DenseSubstrings.substrings[j].endA<Alignments.alignments[id][3]) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && DenseSubstrings.substrings[j].endA>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=Intervals.intersectionLength(Alignments.alignments[id][3],Alignments.alignments[id][4],DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA);
				if (length>THRESHOLD) {
					if (ReadA.sortedAlignments[i].minIntervalLength<0) ReadA.sortedAlignments[i].minIntervalLength--;
					else ReadA.sortedAlignments[i].minIntervalLength++;
				}
				j++;
			}
		}

		if (AlignmentIntervals.lastInterval>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {
				if (filter?(ReadA.sortedAlignments[i].minIntervalType==-3||ReadA.sortedAlignments[i].minIntervalID!=-1):ReadA.sortedAlignments[i].isImplied()) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				id=ReadA.sortedAlignments[i].id;
				if (j>AlignmentIntervals.lastInterval || AlignmentIntervals.intervals[j].firstPosition>Alignments.alignments[id][4]) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (AlignmentIntervals.intervals[j].lastPosition<Alignments.alignments[id][3]) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && AlignmentIntervals.intervals[j].lastPosition>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=Intervals.intersectionLength(Alignments.alignments[id][3],Alignments.alignments[id][4],AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition);
				if (length>THRESHOLD) {
					if (ReadA.sortedAlignments[i].minIntervalLength<0) ReadA.sortedAlignments[i].minIntervalLength--;
					else ReadA.sortedAlignments[i].minIntervalLength++;
				}
				j++;
			}
		}
		
		// Assigning alignments to their only intersecting interval
		out=0;
		if (PeriodicSubstrings.lastInterval>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {
				if ( (filter?(ReadA.sortedAlignments[i].minIntervalType==-3||ReadA.sortedAlignments[i].minIntervalID!=-1):ReadA.sortedAlignments[i].isImplied()) || 
				     ReadA.sortedAlignments[i].minIntervalLength<0 || (ReadA.sortedAlignments[i].minIntervalLength>1 && ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null)
				   ) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				id=ReadA.sortedAlignments[i].id;
				if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>Alignments.alignments[id][4]) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (PeriodicSubstrings.intervals[j].lastPosition<Alignments.alignments[id][3]) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && PeriodicSubstrings.intervals[j].lastPosition>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=Intervals.intersectionLength(Alignments.alignments[id][3],Alignments.alignments[id][4],PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition);
				if (length>0 && ((double)length)/ReadA.sortedAlignments[i].getALength()>=INTERSECTION_THRESHOLD) {
					ReadA.sortedAlignments[i].minIntervalType=Alignment.INTERVAL_PERIODIC;
					ReadA.sortedAlignments[i].minIntervalID=j;
					out++;
				}
				j++;
			}
		}
		if (DenseSubstrings.lastSubstring>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {				
				if ( (filter?(ReadA.sortedAlignments[i].minIntervalType==-3||ReadA.sortedAlignments[i].minIntervalID!=-1):ReadA.sortedAlignments[i].isImplied()) || 
				     ReadA.sortedAlignments[i].minIntervalLength<0 || ReadA.sortedAlignments[i].minIntervalLength>1
				   ) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}				
				id=ReadA.sortedAlignments[i].id;
				if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>Alignments.alignments[id][4]) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (DenseSubstrings.substrings[j].endA<Alignments.alignments[id][3] || !DenseSubstrings.substrings[j].substringReplication) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && DenseSubstrings.substrings[j].endA>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=Intervals.intersectionLength(Alignments.alignments[id][3],Alignments.alignments[id][4],DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA);
				if ( length>0 && ((double)length)/ReadA.sortedAlignments[i].getALength()>=INTERSECTION_THRESHOLD &&
					 DenseSubstrings.substrings[j].canBeAssigned(ReadA.sortedAlignments[i],THRESHOLD)
				   ) {
					ReadA.sortedAlignments[i].minIntervalType=Alignment.INTERVAL_DENSE;
					ReadA.sortedAlignments[i].minIntervalID=j;
					out++;
				}
				j++;
			}
		}
		if (AlignmentIntervals.lastInterval>=0) {
			i=0; j=0; firstJForNextI=-1;
			while (i<=ReadA.lastSortedAlignment) {				
				if ( (filter?(ReadA.sortedAlignments[i].minIntervalType==-3||ReadA.sortedAlignments[i].minIntervalID!=-1):ReadA.sortedAlignments[i].isImplied()) || 
				     ReadA.sortedAlignments[i].minIntervalLength<0 || ReadA.sortedAlignments[i].minIntervalLength>1
				   ) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}				
				id=ReadA.sortedAlignments[i].id;
				if (j>AlignmentIntervals.lastInterval || AlignmentIntervals.intervals[j].firstPosition>Alignments.alignments[id][4]) {
					i++;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (AlignmentIntervals.intervals[j].lastPosition<Alignments.alignments[id][3]) {
					j++;
					continue;
				}
				if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && AlignmentIntervals.intervals[j].lastPosition>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
				length=Intervals.intersectionLength(Alignments.alignments[id][3],Alignments.alignments[id][4],AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition);
				if ( length>0 && ((double)length)/ReadA.sortedAlignments[i].getALength()>=INTERSECTION_THRESHOLD &&
					 AlignmentIntervals.intervals[j].canBeAssigned(ReadA.sortedAlignments[i],THRESHOLD)
				   ) {
					ReadA.sortedAlignments[i].minIntervalType=Alignment.INTERVAL_ALIGNMENT;
					ReadA.sortedAlignments[i].minIntervalID=j;
					out++;
				}
				j++;
			}
		}
		
		return out;
	}
	
	
	/**
	 * Transforms the $minIntervalType$ and $minIntervalID$ fields of an interval into
	 * the corresponding implication pointers.
	 *
	 * @param filter TRUE: the procedure does not consider alignments that are already 
	 * implied by an interval;
	 * @return the number of alignments that are not implied by any interval at the end of
	 * the procedure.
	 */
	private static final int alignments2intervals_updateAlignmentPointers(boolean filter) {
		boolean found;
		int i, j, k;
		int out;
		PeriodicSubstringInterval tmpPInterval;
		
		out=0;
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if (filter && ReadA.sortedAlignments[i].isImplied()) continue;
			if (ReadA.sortedAlignments[i].minIntervalType==Alignment.INTERVAL_PERIODIC) {				
				tmpPInterval=PeriodicSubstrings.intervals[ReadA.sortedAlignments[i].minIntervalID];
				if (tmpPInterval.hasLongPeriod) {
					// Finding the real long-period interval of which $tmpPInterval$ is
					// the clone.
					j=Arrays.binarySearch(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1,tmpPInterval);
					if (j<0) {
						System.err.println("alignments2intervals> ERROR: long-period interval not found.");
						System.exit(1);
					}
					found=false;
					for (k=j; k>=0; k--) {
						if (PeriodicSubstrings.longPeriodIntervals[k].equals(tmpPInterval)) {
							tmpPInterval=PeriodicSubstrings.longPeriodIntervals[k];
							found=true;
							break;
						}
					}
					if (!found) {
						for (k=j+1; k<=PeriodicSubstrings.lastLongPeriodInterval; k++) {
							if (PeriodicSubstrings.longPeriodIntervals[k].equals(tmpPInterval)) {
								tmpPInterval=PeriodicSubstrings.longPeriodIntervals[k];
								break;
							}
						}
					}
				}				
				ReadA.sortedAlignments[i].periodicSubstringInterval=tmpPInterval;
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				ReadA.sortedAlignments[i].inDenseSubstring=null;
				ReadA.sortedAlignments[i].mergedToInterval=null;				
			}
			else if (ReadA.sortedAlignments[i].minIntervalType==Alignment.INTERVAL_DENSE) {
				ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				if (DenseSubstrings.substrings[ReadA.sortedAlignments[i].minIntervalID].substringReplication) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].inDenseSubstring=DenseSubstrings.substrings[ReadA.sortedAlignments[i].minIntervalID];
				}
				else {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=DenseSubstrings.substrings[ReadA.sortedAlignments[i].minIntervalID];
					ReadA.sortedAlignments[i].inDenseSubstring=null;
				}
				ReadA.sortedAlignments[i].mergedToInterval=null;
			}
			else if (ReadA.sortedAlignments[i].minIntervalType==Alignment.INTERVAL_ALIGNMENT) {
				ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				ReadA.sortedAlignments[i].inDenseSubstring=null;
				ReadA.sortedAlignments[i].mergedToInterval=AlignmentIntervals.intervals[ReadA.sortedAlignments[i].minIntervalID];
			}
			else {
				// The alignment is either already assigned to an interval, or it is not
				// assigned to any interval.
				if (filter || !ReadA.sortedAlignments[i].isImplied()) out++;
			}			
		}
		
		return out;
	}
	
	
	/**
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ to be sorted by
	 * $firstPosition$.
	 *
	 * @param orIsContained returns TRUE also if a periodic interval contains 
	 * $[intervalFirst..intervalLast]$ or is approximately identical to it.
	 */
	private static final boolean straddlesOrContains_periodic(int intervalFirst, int intervalLast, int alignmentFirst, int alignmentLast, boolean orIsContained) {
		boolean previousSortByID;
		int i, j;
		
		tmpPeriodic.firstPosition=intervalFirst;
		previousSortByID=PeriodicSubstringInterval.sortByID;
		PeriodicSubstringInterval.sortByID=false;
		i=Arrays.binarySearch(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1,tmpPeriodic);
		PeriodicSubstringInterval.sortByID=previousSortByID;
		if (i<0) {
			i=-i-1;
			i=Math.min(i,PeriodicSubstrings.lastInterval);
		}
		j=i;
		while (j>=0) {
			if ( (PeriodicSubstrings.intervals[j].firstPosition==intervalFirst && PeriodicSubstrings.intervals[j].lastPosition==intervalLast) ||
				 !( Intervals.isApproximatelyContained(alignmentFirst,alignmentLast,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition) ||
				    Intervals.areApproximatelyIdentical(alignmentFirst,alignmentLast,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition)
				  )
			   ) {
				j--;
				continue;
			}
			if ( Intervals.isApproximatelyContained(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,intervalFirst,intervalLast) ||
				 Intervals.straddles(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,intervalFirst,intervalLast) ||
				 ( orIsContained && 
				   ( Intervals.isApproximatelyContained(intervalFirst,intervalLast,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition) ||
					 Intervals.areApproximatelyIdentical(intervalFirst,intervalLast,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition)
				   )
				 )
			   ) return true;
			j--;
		}
		j=i+1;
		while (j<=PeriodicSubstrings.lastInterval) {
			if (PeriodicSubstrings.intervals[j].firstPosition>intervalLast) break;
			if ( (PeriodicSubstrings.intervals[j].firstPosition==intervalFirst && PeriodicSubstrings.intervals[j].lastPosition==intervalLast) ||
				 !( Intervals.isApproximatelyContained(alignmentFirst,alignmentLast,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition) ||
					Intervals.areApproximatelyIdentical(alignmentFirst,alignmentLast,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition)
				  )
			   ) {
				j++;
				continue;
			}
			if ( Intervals.isApproximatelyContained(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,intervalFirst,intervalLast) ||
				 Intervals.straddles(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,intervalFirst,intervalLast) ||
				 ( orIsContained && 
				   ( Intervals.isApproximatelyContained(intervalFirst,intervalLast,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition) ||
					 Intervals.areApproximatelyIdentical(intervalFirst,intervalLast,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition)
				   )
			     )
			   ) return true;
			j++;
		}
		return false;
	}
	
	
	/**
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ to be sorted by
	 * $startA$.
	 */
	private static final boolean straddlesOrContains_dense(int intervalFirst, int intervalLast, int alignmentFirst, int alignmentLast, Alignment alignment, int threshold) {
		boolean previousSortByID;
		int i, j;
		
		tmpDense.startA=intervalFirst;
		previousSortByID=DenseSubstring.sortByID;
		DenseSubstring.sortByID=false;
		i=Arrays.binarySearch(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1,tmpDense);
		DenseSubstring.sortByID=previousSortByID;
		if (i<0) {
			i=-i-1;
			i=Math.min(i,DenseSubstrings.lastSubstring);
		}
		j=i;
		while (j>=0) {
			if ( (DenseSubstrings.substrings[j].startA==intervalFirst && DenseSubstrings.substrings[j].endA==intervalLast) ||
				 !( Intervals.isApproximatelyContained(alignmentFirst,alignmentLast,DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA) ||
					Intervals.areApproximatelyIdentical(alignmentFirst,alignmentLast,DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA)
				  ) ||
				 !DenseSubstrings.substrings[j].canBeAssigned(alignment,threshold)
			   ) {
				j--;
				continue;
			}
			if ( Intervals.isApproximatelyContained(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,intervalFirst,intervalLast) ||
				 Intervals.straddles(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,intervalFirst,intervalLast)
			   ) return true;
			j--;
		}
		j=i+1;
		while (j<=DenseSubstrings.lastSubstring) {
			if (DenseSubstrings.substrings[j].startA>intervalLast) break;
			if ( (DenseSubstrings.substrings[j].startA==intervalFirst && DenseSubstrings.substrings[j].endA==intervalLast) ||
				 !( Intervals.isApproximatelyContained(alignmentFirst,alignmentLast,DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA) ||
					Intervals.areApproximatelyIdentical(alignmentFirst,alignmentLast,DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA)
				  ) ||
				 !DenseSubstrings.substrings[j].canBeAssigned(alignment,threshold)
			   ) {
				j++;
				continue;
			}
			if ( Intervals.isApproximatelyContained(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,intervalFirst,intervalLast) ||
				 Intervals.straddles(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,intervalFirst,intervalLast)
			   ) return true;
			j++;
		}
		return false;
	}
	
	
	/**
	 * Remark: the procedure assumes $AlignmentIntervals.intervals$ to be sorted by
	 * $firstPosition$.
	 */
	private static final boolean straddlesOrContains_alignment(int intervalFirst, int intervalLast, int alignmentFirst, int alignmentLast, Alignment alignment, int threshold) {
		boolean previousSortByID;
		int i, j;
		
		tmpAlignment.firstPosition=intervalFirst;
		previousSortByID=AlignmentInterval.sortByID;
		AlignmentInterval.sortByID=false;
		i=Arrays.binarySearch(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1,tmpAlignment);
		AlignmentInterval.sortByID=previousSortByID;
		if (i<0) {
			i=-i-1;
			i=Math.min(i,AlignmentIntervals.lastInterval);
		}
		j=i;
		while (j>=0) {
			if ( (AlignmentIntervals.intervals[j].firstPosition==intervalFirst && AlignmentIntervals.intervals[j].lastPosition==intervalLast) ||
				 !( Intervals.isApproximatelyContained(alignmentFirst,alignmentLast,AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition) ||
					Intervals.areApproximatelyIdentical(alignmentFirst,alignmentLast,AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition)
				  ) ||
				 !AlignmentIntervals.intervals[j].canBeAssigned(alignment,threshold)
			   ) {
				j--;
				continue;
			}
			if ( Intervals.isApproximatelyContained(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,intervalFirst,intervalLast) ||
				 Intervals.straddles(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,intervalFirst,intervalLast)
			   ) return true;
			j--;
		}
		j=i+1;
		while (j<=AlignmentIntervals.lastInterval) {
			if (AlignmentIntervals.intervals[j].firstPosition>intervalLast) break;
			if ( (AlignmentIntervals.intervals[j].firstPosition==intervalFirst && AlignmentIntervals.intervals[j].lastPosition==intervalLast) ||
				 !( Intervals.isApproximatelyContained(alignmentFirst,alignmentLast,AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition) ||
					Intervals.areApproximatelyIdentical(alignmentFirst,alignmentLast,AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition)
				  ) ||
				 !AlignmentIntervals.intervals[j].canBeAssigned(alignment,threshold)
			   ) {
				j++;
				continue;
			}
			if ( Intervals.isApproximatelyContained(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,intervalFirst,intervalLast) ||
				 Intervals.straddles(AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition,intervalFirst,intervalLast)
			   ) return true;
			j++;
		}
		return false;
	}
	
	
	/**
	 * The file has one line per alignment in LAshow, and no header. A line contains:
	 *
	 * 1. the type of interval the alignment corresponds to in readA;
	 * 2. the ID of such interval;
	 * 3. the start position of such interval;
	 * 4. the end position of such interval;
	 * 5. weak/nonweak, if dense substring interval;
	 * 6. 1 iff the assignment of the alignment to the interval was decided before calling
	 * $alignments2intervals$;
	 *
	 * 7. isLeftMaximal;
	 * 8. isRightMaximal;
	 * 9. isContained;
	 * 10. minPrefixLength (written only if the type of interval is prefix/suffix
	 * or substring);
	 * 11. minSuffixLength (written only if the type of interval is prefix/suffix);
	 *
	 * 12. period (written only if the type of interval is periodic);
	 * 13. hasLongPeriod (written only if the type of interval is periodic);
	 *
	 * 14: firstMaximalStart (written only if the type of interval is prefix/suffix);
	 * 15: lastMaximalEnd (written only if the type of interval is prefix/suffix):
	 *
	 * 16: nAssignedAlignments;
	 * 17: avgCoverage.
	 *
	 * Remark: an alignment could have a $mergedToInterval$ pointer and an 
	 * $inDenseSubstring$ pointer. Priority is given to the former.
	 *
	 * Remark: the procedure sorts $ReadA.sortedAlignments$ by $id$.
	 */
	private static final void writeIntervalConnectionFile() throws IOException {
		int i, type;
		
		// Ensuring the necessary order in $ReadA.sortedAlignments$
		if (Alignment.order!=Alignment.ID) {
			Alignment.order=Alignment.ID;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		
		// Printing
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {			
			if (ReadA.sortedAlignments[i].periodicSubstringInterval!=null) {
				intervalConnectionFile.write(Constants.INTERVAL_PERIODIC+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].periodicSubstringInterval.id+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].periodicSubstringInterval.firstPosition+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].periodicSubstringInterval.lastPosition+",");
				intervalConnectionFile.write("0,");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].minIntervalType<0?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].periodicSubstringInterval.isLeftMaximal?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].periodicSubstringInterval.isRightMaximal?"1":"0")+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].periodicSubstringInterval.isContained+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].periodicSubstringInterval.period+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].periodicSubstringInterval.hasLongPeriod?"1":"0")+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].periodicSubstringInterval.nAssignedAlignments+",");
				intervalConnectionFile.write(IO.format(ReadA.sortedAlignments[i].periodicSubstringInterval.avgCoverage)+"\n");
			}
			else if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null) {
				type=-1;
				if (ReadA.sortedAlignments[i].impliedByDenseSubstring.prefixReplication && !ReadA.sortedAlignments[i].impliedByDenseSubstring.suffixReplication) type=Constants.INTERVAL_DENSE_PREFIX;
				else if (!ReadA.sortedAlignments[i].impliedByDenseSubstring.prefixReplication && ReadA.sortedAlignments[i].impliedByDenseSubstring.suffixReplication) type=Constants.INTERVAL_DENSE_SUFFIX;
				else if (ReadA.sortedAlignments[i].impliedByDenseSubstring.prefixReplication && ReadA.sortedAlignments[i].impliedByDenseSubstring.suffixReplication && !ReadA.sortedAlignments[i].impliedByDenseSubstring.singleDeletionReplication) type=Constants.INTERVAL_DENSE_PREFIXSUFFIX;
				else if (ReadA.sortedAlignments[i].impliedByDenseSubstring.singleDeletionReplication) type=Constants.INTERVAL_DENSE_SINGLEDELETION;
				intervalConnectionFile.write(type+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.id+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.startA+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.endA+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].impliedByDenseSubstring.isWeak?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].minIntervalType<0?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].impliedByDenseSubstring.isLeftMaximal?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].impliedByDenseSubstring.isRightMaximal?"1":"0")+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.isContained+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.minPrefixLength+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.minSuffixLength+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.firstMaximalStart+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.lastMaximalEnd+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].impliedByDenseSubstring.nAssignedAlignments+",");
				intervalConnectionFile.write(IO.format(ReadA.sortedAlignments[i].impliedByDenseSubstring.avgCoverage)+"\n");
			}
			else if (ReadA.sortedAlignments[i].mergedToInterval!=null) {				
				intervalConnectionFile.write(Constants.INTERVAL_ALIGNMENT+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].mergedToInterval.id+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].mergedToInterval.firstPosition+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].mergedToInterval.lastPosition+",");
				intervalConnectionFile.write("0,");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].minIntervalType<0?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].mergedToInterval.isLeftMaximal?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].mergedToInterval.isRightMaximal?"1":"0")+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].mergedToInterval.isContained+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].mergedToInterval.nAssignedAlignments+",");
				intervalConnectionFile.write(IO.format(ReadA.sortedAlignments[i].mergedToInterval.avgCoverage)+"\n");
			}
			else if (ReadA.sortedAlignments[i].inDenseSubstring!=null) {
				intervalConnectionFile.write(Constants.INTERVAL_DENSE_SUBSTRING+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].inDenseSubstring.id+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].inDenseSubstring.startA+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].inDenseSubstring.endA+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].inDenseSubstring.isWeak?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].minIntervalType<0?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].inDenseSubstring.isLeftMaximal?"1":"0")+",");
				intervalConnectionFile.write((ReadA.sortedAlignments[i].inDenseSubstring.isRightMaximal?"1":"0")+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].inDenseSubstring.isContained+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].inDenseSubstring.minPrefixLength+",");
				intervalConnectionFile.write(ReadA.sortedAlignments[i].inDenseSubstring.nAssignedAlignments+",");
				intervalConnectionFile.write(IO.format(ReadA.sortedAlignments[i].inDenseSubstring.avgCoverage)+"\n");
			}
			else intervalConnectionFile.write("\n");
		}
	}
	
	
	/**
	 * Loads row $str$ of the connection file into corresponding global variables
	 */
	public static final void readIntervalConnectionFile(String str2) {
		int p, q;
		String tmp;
		
		p=str2.indexOf(",");
		type=Integer.parseInt(str2.substring(0,p));
		q=str2.indexOf(",",p+1);
		id=Integer.parseInt(str2.substring(p+1,q));
		p=q;
		q=str2.indexOf(",",p+1);
		start=Integer.parseInt(str2.substring(p+1,q));
		p=q;
		q=str2.indexOf(",",p+1);
		end=Integer.parseInt(str2.substring(p+1,q));  // No need to subtract 1
		p=q;
		q=str2.indexOf(",",p+1);
		tmp=str2.substring(p+1,q);
		isWeak=Integer.parseInt(tmp)==1;
		p=q+1;
		q=str2.indexOf(",",p+1);
		implied=Integer.parseInt(str2.substring(p,q));
		p=q+1;
		q=str2.indexOf(",",p+1);
		isLeftMaximal=Integer.parseInt(str2.substring(p,q))==1;
		p=q+1;
		q=str2.indexOf(",",p+1);
		isRightMaximal=Integer.parseInt(str2.substring(p,q))==1;
		p=q+1;
		period=-1; minPrefixLength=-1; minSuffixLength=-1; 
		firstMaximalStart=-1; lastMaximalEnd=-1;
		if (type==Constants.INTERVAL_ALIGNMENT) {
			q=str2.indexOf(",",p+1);
			isContained=Integer.parseInt(str2.substring(p,q));
			p=q+1; q=str2.indexOf(",",p+1);
			nAssignedAlignments=Integer.parseInt(str2.substring(p,q));
			p=q+1;
			avgCoverage=Double.parseDouble(str2.substring(p));
			return;
		}
		else if (type==Constants.INTERVAL_DENSE_SUBSTRING) {	
			q=str2.indexOf(",",p+1);
			isContained=Integer.parseInt(str2.substring(p,q));
			p=q+1; q=str2.indexOf(",",p+1);
			minPrefixLength=Integer.parseInt(str2.substring(p,q));
			p=q+1; q=str2.indexOf(",",p+1);
			nAssignedAlignments=Integer.parseInt(str2.substring(p,q));
			p=q+1;
			avgCoverage=Double.parseDouble(str2.substring(p));
			return;
		}
		q=str2.indexOf(",",p+1);
		isContained=Integer.parseInt(str2.substring(p,q));
		if (type==Constants.INTERVAL_DENSE_PREFIX || type==Constants.INTERVAL_DENSE_SUFFIX || type==Constants.INTERVAL_DENSE_PREFIXSUFFIX || type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
			p=q+1;
			q=str2.indexOf(",",p+1);
			minPrefixLength=Integer.parseInt(str2.substring(p,q));
			p=q+1; q=str2.indexOf(",",p+1);
			minSuffixLength=Integer.parseInt(str2.substring(p,q));
			p=q+1; q=str2.indexOf(",",p+1);
			firstMaximalStart=Integer.parseInt(str2.substring(p,q));
			p=q+1; q=str2.indexOf(",",p+1);
			lastMaximalEnd=Integer.parseInt(str2.substring(p,q));
			p=q+1; q=str2.indexOf(",",p+1);
			nAssignedAlignments=Integer.parseInt(str2.substring(p,q));
			p=q+1;
			avgCoverage=Double.parseDouble(str2.substring(p));
		}
		else if (type==Constants.INTERVAL_PERIODIC) {
			p=q+1;
			q=str2.indexOf(",",p+1);
			period=Integer.parseInt(str2.substring(p,q));
			p=q+1; q=str2.indexOf(",",p+1);
			hasLongPeriod=Integer.parseInt(str2.substring(p,q))==1;
			p=q+1; q=str2.indexOf(",",p+1);
			nAssignedAlignments=Integer.parseInt(str2.substring(p,q));
			p=q+1;
			avgCoverage=Double.parseDouble(str2.substring(p));
		}
	}
	
	
	/**
	 * Sets the $isContained$ value of all intervals. The procedure distinguishes just
	 * between dense substrings of substring type (bit $INTERVAL_DENSE_SUBSTRING$) and all
	 * others (bit $INTERVAL_DENSE_PREFIX$). The procedure uses bit $INTERVAL_PERIODIC$ 
	 * for short-period, and $INTERVAL_PERIODIC+1$ for long period.
	 *
	 * Remark: the procedure marks bit $INTERVAL_PERIODIC$ of an interval even when it is 
	 * contained in the union of short-period periodic intervals.
	 *
	 * Remark: the procedure assumes that all intervals and dense substrings have already
	 * been sorted by startA.
	 */
	private static final void markContainedIntervals() {
		boolean isContained, previousSortByID;
		int i, k, start, end;
		PeriodicSubstringInterval tmpPeriodic = new PeriodicSubstringInterval();
		DenseSubstring tmpSubstring = new DenseSubstring();
		AlignmentInterval tmpAlignment = new AlignmentInterval();
		
		i=(PeriodicSubstrings.lastInterval+1)<<1;
		if (tmpArray.length<i) tmpArray = new int[i];
		
		// Periodic
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			PeriodicSubstrings.intervals[i].isContained=0;
			start=PeriodicSubstrings.intervals[i].firstPosition;
			end=PeriodicSubstrings.intervals[i].lastPosition;
			// Periodic
			isContained_periodic(start,end,i,false,tmpArray);
			if (tmpArray[0]==1) PeriodicSubstrings.intervals[i].isContained|=1<<Constants.INTERVAL_PERIODIC;
			if (tmpArray[1]==1) PeriodicSubstrings.intervals[i].isContained|=1<<(Constants.INTERVAL_PERIODIC+1);
			isContained=isContainedInUnion_periodic(start,end,i,false,tmpArray);
			if (isContained) PeriodicSubstrings.intervals[i].isContained|=1<<Constants.INTERVAL_PERIODIC;
			// Dense
			if (DenseSubstrings.lastSubstring>=0) {
				tmpSubstring.startA=start;
				previousSortByID=DenseSubstring.sortByID;
				DenseSubstring.sortByID=false;
				k=Arrays.binarySearch(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1,tmpSubstring);
				DenseSubstring.sortByID=previousSortByID;
				if (k<0) {
					k=-k-1;
					k=Math.min(k,DenseSubstrings.lastSubstring);
				}
				isContained_dense(start,end,k,true,tmpArray);
				if (tmpArray[0]==1) PeriodicSubstrings.intervals[i].isContained|=1<<Constants.INTERVAL_DENSE_PREFIX;  // We do not distinguish among dense substring types.
				if (tmpArray[1]==1) PeriodicSubstrings.intervals[i].isContained|=1<<Constants.INTERVAL_DENSE_SUBSTRING;
				else {
					isContained=isContainedInUnion_dense(start,end,k,true,tmpArray);
					if (isContained) PeriodicSubstrings.intervals[i].isContained|=1<<Constants.INTERVAL_DENSE_SUBSTRING;
				}
			}
			// Alignment
			if (AlignmentIntervals.lastInterval>=0) {
				tmpAlignment.firstPosition=start;
				previousSortByID=AlignmentInterval.sortByID;
				AlignmentInterval.sortByID=false;
				k=Arrays.binarySearch(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1,tmpAlignment);
				AlignmentInterval.sortByID=previousSortByID;
				if (k<0) {
					k=-k-1;
					k=Math.min(k,AlignmentIntervals.lastInterval);
				}
				isContained=isContained_alignment(start,end,k,true);
				if (isContained) PeriodicSubstrings.intervals[i].isContained|=1<<Constants.INTERVAL_ALIGNMENT;
			}
		}
		
		// Dense
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			DenseSubstrings.substrings[i].isContained=0;
			start=DenseSubstrings.substrings[i].startA;
			end=DenseSubstrings.substrings[i].endA;
			// Dense
			isContained_dense(start,end,i,false,tmpArray);
			if (tmpArray[0]==1) DenseSubstrings.substrings[i].isContained|=1<<Constants.INTERVAL_DENSE_PREFIX; 
			if (tmpArray[1]==1) DenseSubstrings.substrings[i].isContained|=1<<Constants.INTERVAL_DENSE_SUBSTRING;
			else {
				isContained=isContainedInUnion_dense(start,end,i,false,tmpArray);
				if (isContained) DenseSubstrings.substrings[i].isContained|=1<<Constants.INTERVAL_DENSE_SUBSTRING;
			}
			// Periodic
			if (PeriodicSubstrings.lastInterval>=0) {
				tmpPeriodic.firstPosition=start;
				previousSortByID=PeriodicSubstringInterval.sortByID;
				PeriodicSubstringInterval.sortByID=false;
				k=Arrays.binarySearch(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1,tmpPeriodic);
				PeriodicSubstringInterval.sortByID=previousSortByID;
				if (k<0) {
					k=-k-1;
					k=Math.min(k,PeriodicSubstrings.lastInterval);
				}
				isContained_periodic(start,end,k,true,tmpArray);
				if (tmpArray[0]==1) DenseSubstrings.substrings[i].isContained|=1<<Constants.INTERVAL_PERIODIC;
				if (tmpArray[1]==1) DenseSubstrings.substrings[i].isContained|=1<<(Constants.INTERVAL_PERIODIC+1);
				isContained=isContainedInUnion_periodic(start,end,k,true,tmpArray);
				if (isContained) DenseSubstrings.substrings[i].isContained|=1<<Constants.INTERVAL_PERIODIC;
			}
			// Alignment
			if (AlignmentIntervals.lastInterval>=0) {
				tmpAlignment.firstPosition=start;
				previousSortByID=AlignmentInterval.sortByID;
				AlignmentInterval.sortByID=false;
				k=Arrays.binarySearch(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1,tmpAlignment);
				AlignmentInterval.sortByID=previousSortByID;
				if (k<0) {
					k=-k-1;
					k=Math.min(k,AlignmentIntervals.lastInterval);
				}
				isContained=isContained_alignment(start,end,k,true);
				if (isContained) DenseSubstrings.substrings[i].isContained|=1<<Constants.INTERVAL_ALIGNMENT;
			}
		}

		// Alignment
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			AlignmentIntervals.intervals[i].isContained=0;
			start=AlignmentIntervals.intervals[i].firstPosition;
			end=AlignmentIntervals.intervals[i].lastPosition;
			// Alignment
			isContained=isContained_alignment(start,end,i,false);
			if (isContained) AlignmentIntervals.intervals[i].isContained|=1<<Constants.INTERVAL_ALIGNMENT;
			// Periodic
			if (PeriodicSubstrings.lastInterval>=0) {
				tmpPeriodic.firstPosition=start;
				previousSortByID=PeriodicSubstringInterval.sortByID;
				PeriodicSubstringInterval.sortByID=false;
				k=Arrays.binarySearch(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1,tmpPeriodic);
				PeriodicSubstringInterval.sortByID=previousSortByID;
				if (k<0) {
					k=-k-1;
					k=Math.min(k,PeriodicSubstrings.lastInterval);
				}
				isContained_periodic(start,end,k,true,tmpArray);
				if (tmpArray[0]==1) AlignmentIntervals.intervals[i].isContained|=1<<Constants.INTERVAL_PERIODIC;
				if (tmpArray[1]==1) AlignmentIntervals.intervals[i].isContained|=1<<(Constants.INTERVAL_PERIODIC+1);
				isContained=isContainedInUnion_periodic(start,end,k,true,tmpArray);
				if (isContained) DenseSubstrings.substrings[i].isContained|=1<<Constants.INTERVAL_PERIODIC;
			}
			// Dense
			if (DenseSubstrings.lastSubstring>=0) {
				tmpSubstring.startA=start;
				previousSortByID=DenseSubstring.sortByID;
				DenseSubstring.sortByID=false;
				k=Arrays.binarySearch(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1,tmpSubstring);
				DenseSubstring.sortByID=previousSortByID;
				if (k<0) {
					k=-k-1;
					k=Math.min(k,DenseSubstrings.lastSubstring);
				}
				isContained_dense(start,end,k,true,tmpArray);
				if (tmpArray[0]==1) AlignmentIntervals.intervals[i].isContained|=1<<Constants.INTERVAL_DENSE_PREFIX;
				if (tmpArray[1]==1) AlignmentIntervals.intervals[i].isContained|=1<<Constants.INTERVAL_DENSE_SUBSTRING;
				else {
					isContained=isContainedInUnion_dense(start,end,k,true,tmpArray);
					if (isContained) AlignmentIntervals.intervals[i].isContained|=1<<Constants.INTERVAL_DENSE_SUBSTRING;
				}
			}
		}
	}
	
	
	/**
	 * @return out[0]: contained in short-period interval; out[1]: contained in 
	 * long-period interval (1/0).
	 */
	private static final void isContained_periodic(int start, int end, int k, boolean includeK, int[] out) {
		int i;
		
		out[0]=0; out[1]=0;
		for (i=k-1; i>=0; i--) {
			if ( Intervals.isApproximatelyContained(start,end,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition) &&
				 !Intervals.areApproximatelyIdentical(start,end,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition)
			   ) {
				out[PeriodicSubstrings.intervals[i].hasLongPeriod?1:0]=1;
				if (out[0]==1 && out[1]==1) return;
			}
		}
		for (i=includeK?k:k+1; i>=0 && i<=PeriodicSubstrings.lastInterval; i++) {
			if (PeriodicSubstrings.intervals[i].firstPosition>=end) break;
			if ( Intervals.isApproximatelyContained(start,end,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition) &&
				 !Intervals.areApproximatelyIdentical(start,end,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition)
			   ) {
   				out[PeriodicSubstrings.intervals[i].hasLongPeriod?1:0]=1;
   				if (out[0]==1 && out[1]==1) return;
			}
		}
	}
	
	
	/**
	 * @return true iff $[start..end]$ is approximately contained in a union of 
	 * short-period periodic substring intervals without large gaps.
	 */
	private static final boolean isContainedInUnion_periodic(int start, int end, int k, boolean includeK, int[] tmpArray) {
		final int GAP_THRESHOLD = IO.quantum;
		int i, j, p;
		int tmp1, tmp2, first, last, lastInterval;
		
		// Collecting all periodic substring intervals that intersect with $[start..end]$.
		j=-1; p=-1;
		for (i=k-1; i>=0; i--) {
			if (PeriodicSubstrings.intervals[i].hasLongPeriod) continue;
			if (PeriodicSubstrings.intervals[i].lastPosition<start) continue;
			if (PeriodicSubstrings.intervals[i].lastPosition>=end) return true;
			tmpArray[++j]=PeriodicSubstrings.intervals[i].firstPosition;
			tmpArray[++j]=PeriodicSubstrings.intervals[i].lastPosition;
		}
		p=j-1;
		for (i=includeK?k:k+1; i<=PeriodicSubstrings.lastInterval; i++) {
			if (PeriodicSubstrings.intervals[i].firstPosition>=end) break;
			if (PeriodicSubstrings.intervals[i].hasLongPeriod) continue;
			tmpArray[++j]=PeriodicSubstrings.intervals[i].firstPosition;
			tmpArray[++j]=PeriodicSubstrings.intervals[i].lastPosition;
		}
		lastInterval=j-1;
		
		// Reversing $tmpArray[0..p+1]$, so all intervals are sorted by first position.
		for (i=p; i>=(p+2)>>1; i-=2) {
			tmp1=tmpArray[p-i]; tmp2=tmpArray[p-i+1];
			tmpArray[p-i]=tmpArray[i]; tmpArray[p-i+1]=tmpArray[i+1];
			tmpArray[i]=tmp1; tmpArray[i+1]=tmp2;
		}
		
		// Finding large gaps
		first=tmpArray[0]; last=tmpArray[1];
		for (i=2; i<=lastInterval; i+=2) {
			if (tmpArray[i]>last+GAP_THRESHOLD) return false;
			if (tmpArray[i+1]>last) last=tmpArray[i+1];
		}
		return Intervals.isApproximatelyContained(start,end,first,last);
	}
	
	
	/**
	 * @return out[0]: contained in a dense substring that is not of substring type (1/0);
	 * out[1]: contained in a dense substring of substring type.
	 */
	private static final void isContained_dense(int start, int end, int k, boolean includeK, int[] out) {
		int i;
		
		for (i=k-1; i>=0; i--) {
			if ( Intervals.isApproximatelyContained(start,end,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA) &&
				 !Intervals.areApproximatelyIdentical(start,end,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA)
			   ) {
				out[DenseSubstrings.substrings[i].substringReplication?1:0]=1;   	
				if (out[0]==1 && out[1]==1) return;
			}
		}
		for (i=includeK?k:k+1; i<=DenseSubstrings.lastSubstring; i++) {
			if (DenseSubstrings.substrings[i].startA>=end) break;
			if ( Intervals.isApproximatelyContained(start,end,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA) &&
				 !Intervals.areApproximatelyIdentical(start,end,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA)
			   ) {
   				out[DenseSubstrings.substrings[i].substringReplication?1:0]=1;   	
   				if (out[0]==1 && out[1]==1) return;
			}
		}
	}
	
	
	/**
	 * @return true iff $[start..end]$ is approximately contained in a union of 
	 * substring-type intervals without large gaps.
	 */
	private static final boolean isContainedInUnion_dense(int start, int end, int k, boolean includeK, int[] tmpArray) {
		final int GAP_THRESHOLD = IO.quantum;
		int i, j, p;
		int tmp1, tmp2, first, last, lastInterval;
		
		// Collecting all substring intervals that intersect with $[start..end]$.
		j=-1; p=-1;
		for (i=k-1; i>=0; i--) {
			if (!DenseSubstrings.substrings[i].substringReplication) continue;
			if (DenseSubstrings.substrings[i].endA<start) continue;
			if (DenseSubstrings.substrings[i].startA>=end) return true;
			tmpArray[++j]=DenseSubstrings.substrings[i].startA;
			tmpArray[++j]=DenseSubstrings.substrings[i].endA;
		}
		p=j-1;
		for (i=includeK?k:k+1; i<=DenseSubstrings.lastSubstring; i++) {
			if (DenseSubstrings.substrings[i].startA>=end) break;
			if (!DenseSubstrings.substrings[i].substringReplication) continue;
			tmpArray[++j]=DenseSubstrings.substrings[i].startA;
			tmpArray[++j]=DenseSubstrings.substrings[i].endA;
		}
		lastInterval=j-1;
		
		// Reversing $tmpArray[0..p+1]$, so all intervals are sorted by first position.
		for (i=p; i>=(p+2)>>1; i-=2) {
			tmp1=tmpArray[p-i]; tmp2=tmpArray[p-i+1];
			tmpArray[p-i]=tmpArray[i]; tmpArray[p-i+1]=tmpArray[i+1];
			tmpArray[i]=tmp1; tmpArray[i+1]=tmp2;
		}
		
		// Finding large gaps
		first=tmpArray[0]; last=tmpArray[1];
		for (i=2; i<=lastInterval; i+=2) {
			if (tmpArray[i]>last+GAP_THRESHOLD) return false;
			if (tmpArray[i+1]>last) last=tmpArray[i+1];
		}
		return Intervals.isApproximatelyContained(start,end,first,last);
	}
	
		
	private static final boolean isContained_alignment(int start, int end, int k, boolean includeK) {
		int i;
		
		for (i=k-1; i>=0; i--) {
			if ( Intervals.isApproximatelyContained(start,end,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition) &&
				 !Intervals.areApproximatelyIdentical(start,end,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition)
			   ) return true;
		}
		for (i=includeK?k:k+1; i<=AlignmentIntervals.lastInterval; i++) {
			if (AlignmentIntervals.intervals[i].firstPosition>=end) break;
			if ( Intervals.isApproximatelyContained(start,end,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition) &&
				 !Intervals.areApproximatelyIdentical(start,end,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition)
			   ) return true;
		}
		return false;
	}
	
	
	/**
	 * Given the files produced by the factorization pipeline, the procedure computes the
	 * number of periodic intervals ($out[0]$), dense substring intervals ($out[1]$), and 
	 * alignment intervals ($out[2]$) that have not been assigned to any alignment.
	 * Not being assigned to any interval does not necessarily imply an error in the 
	 * interval or in the assignment of alignments, as described below. However, an 
	 * unassigned interval will not be included in the interval graph by the following 
	 * steps of the pipeline.
	 *
	 * Remark: the interval of a dense substring of substring type might not have any 
	 * alignment assigned, since e.g. all its alignments might have been assigned to 
	 * alignment intervals that align with peaks in the delta function inside the dense 
	 * substring. In practice this happens frequently with weak dense substrings, since 
	 * they often have biases in the first/last positions of their alignments.
	 *
	 * Remark: all the alignments assigned to a weak dense substring of substring type are
	 * also potentially reassigned in $alignments2intervals()$.
	 *
	 * Remark: a weak dense substring of substring type might not have any alignment
	 * assigned, since alignments might be assigned instead to a shorter weak dense
	 * substring of substring type that has high Jaccard similarity with it but is not
	 * merged with it.
	 *
	 * Remark: a prefix/suffix interval might not have any alignment assigned, since e.g.
	 * $DenseSubstrings.getSubstrings()$ at the end reassigns an alignment that was 
	 * implied both by a prefix and a suffix substring interval to exactly one interval, 
	 * even randomly in ambiguous cases. This is the price to pay for the decision of 
	 * assigning an alignment to just one dense substring.
	 */
	public static final void intervalsNotAssigned(String aligmentsFilePath, String connectionFilePath, String periodicSubstringsPath, String denseSubstringsPath, String alignmentsPath, int[] out) throws IOException {
		final int MAX_INTERVALS_PER_READ = 30000;
		final int READ_AHEAD_LIMIT = MAX_INTERVALS_PER_READ*100;  // Characters
		int i;
		int id, type, read, previousRead, lastIntervalPeriodic, lastIntervalDense, lastIntervalAlignment;
		int notAssignedPeriodic, notAssignedDense, notAssignedAlignment;
		String str1, str2;
		BufferedReader alignmentsFile, connectionFile, periodicBuffer, denseBuffer, alignmentsBuffer;
		boolean[] assignedPeriodic, assignedDense, assignedAlignment;
		int[] tmpArray = new int[30];
		int[][] periodicIntervals, denseIntervals, alignmentIntervals;
		
		alignmentsFile = new BufferedReader(new FileReader(aligmentsFilePath),IO.BUFFER_SIZE);
		connectionFile = new BufferedReader(new FileReader(connectionFilePath),IO.BUFFER_SIZE);
		periodicBuffer = new BufferedReader(new FileReader(periodicSubstringsPath),IO.BUFFER_SIZE);
		periodicBuffer.mark(READ_AHEAD_LIMIT);
		denseBuffer = new BufferedReader(new FileReader(denseSubstringsPath),IO.BUFFER_SIZE);
		denseBuffer.mark(READ_AHEAD_LIMIT);
		alignmentsBuffer = new BufferedReader(new FileReader(alignmentsPath),IO.BUFFER_SIZE);
		alignmentsBuffer.mark(READ_AHEAD_LIMIT);
		periodicIntervals = new int[MAX_INTERVALS_PER_READ][11];
		denseIntervals = new int[MAX_INTERVALS_PER_READ][14];
		alignmentIntervals = new int[MAX_INTERVALS_PER_READ][8];
		assignedPeriodic = new boolean[MAX_INTERVALS_PER_READ];
		assignedDense = new boolean[MAX_INTERVALS_PER_READ];
		assignedAlignment = new boolean[MAX_INTERVALS_PER_READ];
		
		str1=alignmentsFile.readLine();
		str1=alignmentsFile.readLine();  // Skipping the first two lines
		str1=alignmentsFile.readLine();
		str2=connectionFile.readLine();
		previousRead=-1;
		notAssignedPeriodic=0; notAssignedDense=0; notAssignedAlignment=0;
		lastIntervalPeriodic=-1; lastIntervalDense=-1; lastIntervalAlignment=-1;
		while (str1!=null) {
			Alignments.readAlignmentFile(str1);
			read=Alignments.readA-1;  // Starting from zero
			if (read!=previousRead) {
				// Printing unassigned intervals
				if (previousRead!=-1) {
					for (i=0; i<=lastIntervalPeriodic; i++) {
						if (assignedPeriodic[i]) continue;
						notAssignedPeriodic++;
						System.out.println(previousRead+","+Constants.INTERVAL_PERIODIC+","+periodicIntervals[i][0]+","+periodicIntervals[i][1]+","+periodicIntervals[i][2]);
					}
					Math.set(assignedPeriodic,lastIntervalPeriodic,false);
					for (i=0; i<=lastIntervalDense; i++) {
						if (assignedDense[i]) continue;
						notAssignedDense++;
						System.out.println(previousRead+","+denseIntervals[i][1]+","+denseIntervals[i][0]+","+denseIntervals[i][2]+","+denseIntervals[i][3]+", weak="+denseIntervals[i][6]);
					}
					Math.set(assignedDense,lastIntervalDense,false);
					for (i=0; i<=lastIntervalAlignment; i++) {
						if (assignedAlignment[i]) continue;
						notAssignedAlignment++;
						System.out.println(previousRead+","+Constants.INTERVAL_ALIGNMENT+","+alignmentIntervals[i][0]+","+alignmentIntervals[i][1]+","+alignmentIntervals[i][2]);
					}
					Math.set(assignedAlignment,lastIntervalAlignment,false);
				}
				// Next iteration
				lastIntervalPeriodic=PeriodicSubstringInterval.loadPeriodicIntervals(read,periodicBuffer,READ_AHEAD_LIMIT,periodicIntervals,tmpArray);
				Math.set(assignedPeriodic,lastIntervalPeriodic,false);
				lastIntervalDense=DenseSubstring.loadDenseIntervals(read,denseBuffer,READ_AHEAD_LIMIT,denseIntervals,tmpArray);
				Math.set(assignedDense,lastIntervalDense,false);
				lastIntervalAlignment=AlignmentInterval.loadAlignmentIntervals(read,alignmentsBuffer,READ_AHEAD_LIMIT,alignmentIntervals,tmpArray);
				Math.set(assignedDense,lastIntervalAlignment,false);
				previousRead=read;
			}
			if (str2.length()==0) {
				str1=alignmentsFile.readLine();
				str2=connectionFile.readLine();
				continue;
			}
			Factorize.readIntervalConnectionFile(str2);
			id=Factorize.id; type=Factorize.type;
			if (type==Constants.INTERVAL_PERIODIC) {
				for (i=0; i<=lastIntervalPeriodic; i++) {
					if (id==periodicIntervals[i][0]) {
						assignedPeriodic[i]=true;
						break;
					}
				}
			}
			else if (type>=Constants.INTERVAL_DENSE_PREFIX && Factorize.type<=Constants.INTERVAL_DENSE_SINGLEDELETION) {
				for (i=0; i<=lastIntervalDense; i++) {
					if (id==denseIntervals[i][0]) {
						assignedDense[i]=true;
						break;
					}
				}
			}
			else if (type==Constants.INTERVAL_ALIGNMENT) {
				for (i=0; i<=lastIntervalAlignment; i++) {
					if (id==alignmentIntervals[i][0]) {
						assignedAlignment[i]=true;
						break;
					}
				}
			}
			else {
				System.err.println("ERROR: invalid interval type:");
				System.err.println(str2);
				System.exit(1);
			}
			str1=alignmentsFile.readLine();
			str2=connectionFile.readLine();
		}
		for (i=0; i<=lastIntervalPeriodic; i++) {
			if (assignedPeriodic[i]) continue;
			notAssignedPeriodic++;
			System.out.println(previousRead+","+Constants.INTERVAL_PERIODIC+","+periodicIntervals[i][0]+","+periodicIntervals[i][1]+","+periodicIntervals[i][2]);
		}
		periodicBuffer.close();
		for (i=0; i<=lastIntervalDense; i++) {
			if (assignedDense[i]) continue;
			notAssignedDense++;
			System.out.println(previousRead+","+denseIntervals[i][1]+","+denseIntervals[i][0]+","+denseIntervals[i][2]+","+denseIntervals[i][3]);
		}
		denseBuffer.close();
		for (i=0; i<=lastIntervalAlignment; i++) {
			if (assignedAlignment[i]) continue;
			notAssignedAlignment++;
			System.out.println(previousRead+","+Constants.INTERVAL_ALIGNMENT+","+alignmentIntervals[i][0]+","+alignmentIntervals[i][1]+","+alignmentIntervals[i][2]);
		}
		alignmentsBuffer.close();
		alignmentsFile.close(); connectionFile.close();
		out[0]=notAssignedPeriodic; out[1]=notAssignedDense; out[2]=notAssignedAlignment;
	}
	
	
	/**
	 * Sets the $minPrefixLength,minSuffixLength,firstMaximalStart,lastMaximalEnd$ fields
	 * of substrings of prefix, suffix, prefix+suffix, and substring type.
	 *
	 * Remark: the latter fields are the last position of a B-right-maximal alignment that 
	 * is not approximately equal to the end of the substring -- respectively, the first 
	 * position of a B-left-maximal alignment that is not approximately equal to the start 
	 * of the substring. These values are absolute, i.e. not relative to the first 
	 * position of a substring. If the substring does not contain a spurious end/start, 
	 * those positions are just outside the equality threshold from the end/start of the 
	 * substring. 
	 */
	private static final void assignFirstLastMaximal() {
		final int THRESHOLD = IO.quantum;
		int i;
		int id, start, end, length;
		DenseSubstring substring;
		
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			substring=DenseSubstrings.substrings[i];
			// Using temporary fields $leftAlignment$ and $rightAlignment$ to cumulate
			// min start position and max end position.
			substring.leftAlignment=Math.POSITIVE_INFINITY;
			substring.rightAlignment=-1;
			// Using temporary fields $previousSumStartA$ and $previousSumEndA$ to
			// cumulate min prefix length and min suffix length.
			substring.previousSumStartA=Math.POSITIVE_INFINITY;
			substring.previousSumEndA=Math.POSITIVE_INFINITY;
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			id=ReadA.sortedAlignments[i].id;
			start=Alignments.alignments[id][3];
			end=Alignments.alignments[id][4];
			length=end-start+1;
			if (substring!=null) {
				if ((substring.prefixReplication||substring.singleDeletionReplication) && ReadA.sortedAlignments[i].isRightMaximalB==1) {
					if (Math.abs(start,substring.startA)<=THRESHOLD && length<substring.previousSumStartA) substring.previousSumStartA=length;
					if (end<substring.endA-THRESHOLD && end>substring.rightAlignment) substring.rightAlignment=end;
				}
				if ((substring.suffixReplication||substring.singleDeletionReplication) && ReadA.sortedAlignments[i].isLeftMaximalB==1) {
					if (Math.abs(end,substring.endA)<=THRESHOLD && length<substring.previousSumEndA) substring.previousSumEndA=length;
					if (start>substring.startA+THRESHOLD && start<substring.leftAlignment) substring.leftAlignment=start;
				}
			}
			substring=ReadA.sortedAlignments[i].inDenseSubstring;
			if (substring!=null) {
				if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && ReadA.sortedAlignments[i].isRightMaximalB==1) {
					if (length<substring.previousSumStartA) substring.previousSumStartA=length;
				}
			}
		}
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			substring=DenseSubstrings.substrings[i];
			if (DenseSubstrings.substrings[i].leftAlignment!=Math.POSITIVE_INFINITY) DenseSubstrings.substrings[i].firstMaximalStart=DenseSubstrings.substrings[i].leftAlignment;
			else DenseSubstrings.substrings[i].firstMaximalStart=-1;
			if (DenseSubstrings.substrings[i].rightAlignment!=-1) DenseSubstrings.substrings[i].lastMaximalEnd=DenseSubstrings.substrings[i].rightAlignment;
			else DenseSubstrings.substrings[i].lastMaximalEnd=-1;
			if (DenseSubstrings.substrings[i].previousSumStartA!=Math.POSITIVE_INFINITY) DenseSubstrings.substrings[i].minPrefixLength=DenseSubstrings.substrings[i].previousSumStartA;
			else DenseSubstrings.substrings[i].minPrefixLength=-1;
			if (DenseSubstrings.substrings[i].previousSumEndA!=Math.POSITIVE_INFINITY) DenseSubstrings.substrings[i].minSuffixLength=DenseSubstrings.substrings[i].previousSumEndA;
			else DenseSubstrings.substrings[i].minSuffixLength=-1;
		}
	}
	
	
	/**
	 * Computes the number of alignments assigned to each interval, and the avg. coverage
	 * of the interval by its assigned alignments.
	 */
	private static final void averageCoverageAssigned() {
		boolean previousSortByID;
		int i, j;
		PeriodicSubstringInterval periodicInterval;
		DenseSubstring denseSubstring;
		AlignmentInterval alignmentInterval;
		Alignment alignment;
		
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			periodicInterval=PeriodicSubstrings.intervals[i];
			periodicInterval.avgCoverage=0.0;
			periodicInterval.nAssignedAlignments=0;
		}
		for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) {
			periodicInterval=PeriodicSubstrings.longPeriodIntervals[i];
			periodicInterval.avgCoverage=0.0;
			periodicInterval.nAssignedAlignments=0;
		}
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			denseSubstring=DenseSubstrings.substrings[i];
			denseSubstring.avgCoverage=0.0;
			denseSubstring.nAssignedAlignments=0;
		}
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			alignmentInterval=AlignmentIntervals.intervals[i];
			alignmentInterval.avgCoverage=0.0;
			alignmentInterval.nAssignedAlignments=0;
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			alignment=ReadA.sortedAlignments[i];
			periodicInterval=alignment.periodicSubstringInterval;
			if (periodicInterval!=null) {
				periodicInterval.nAssignedAlignments++;
				periodicInterval.avgCoverage+=Intervals.intersectionLength(alignment.startA(),alignment.endA(),periodicInterval.firstPosition,periodicInterval.lastPosition);
			}
			denseSubstring=alignment.impliedByDenseSubstring;
			if (denseSubstring!=null) {
				denseSubstring.nAssignedAlignments++;
				denseSubstring.avgCoverage+=Intervals.intersectionLength(alignment.startA(),alignment.endA(),denseSubstring.startA,denseSubstring.endA);
			}
			denseSubstring=alignment.inDenseSubstring;
			if (denseSubstring!=null) {
				denseSubstring.nAssignedAlignments++;
				denseSubstring.avgCoverage+=Intervals.intersectionLength(alignment.startA(),alignment.endA(),denseSubstring.startA,denseSubstring.endA);
			}
			alignmentInterval=alignment.mergedToInterval;
			if (alignmentInterval!=null) {
				alignmentInterval.nAssignedAlignments++;
				alignmentInterval.avgCoverage+=Intervals.intersectionLength(alignment.startA(),alignment.endA(),alignmentInterval.firstPosition,alignmentInterval.lastPosition);
			}
		}
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			periodicInterval=PeriodicSubstrings.intervals[i];
			periodicInterval.avgCoverage/=periodicInterval.length();
		}
		for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) {
			periodicInterval=PeriodicSubstrings.longPeriodIntervals[i];
			periodicInterval.avgCoverage/=periodicInterval.length();
		}
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			denseSubstring=DenseSubstrings.substrings[i];
			denseSubstring.avgCoverage/=denseSubstring.length();
		}
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			alignmentInterval=AlignmentIntervals.intervals[i];
			alignmentInterval.avgCoverage/=alignmentInterval.length();
		}
		
		// Copying values from $PeriodicSubstrings.longPeriodIntervals$ to
		// $PeriodicSubstrings.intervals$.
		previousSortByID=PeriodicSubstringInterval.sortByID;
		PeriodicSubstringInterval.sortByID=false;
		for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) {
			periodicInterval=PeriodicSubstrings.longPeriodIntervals[i];
			j=Arrays.binarySearch(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1,periodicInterval);
			if (j<0) {
				System.err.println("averageCoverageAssigned> ERROR: long-period interval not found?");
				System.err.println(periodicInterval);
				System.exit(1);
			}
			PeriodicSubstrings.intervals[j].nAssignedAlignments=periodicInterval.nAssignedAlignments;
			PeriodicSubstrings.intervals[j].avgCoverage=periodicInterval.avgCoverage;
		}
		PeriodicSubstringInterval.sortByID=previousSortByID;
	}
	
	
	/**
	 * [THIS PROCEDURE IS CURRENTLY UNUSED]
	 * Builds the coverage histogram of every interval, using just the alignments 
	 * assigned to it after $alignments2intervals()$. In practice such histograms are very
	 * wide and do not conform to the theory (e.g. prefix substrings do not have a line
	 * with slope as a histogram, and alignments don't have a horizontal line as a
	 * histogram). This prevents using such histograms to decide whether to keep an edge
	 * of the interval graph. Using the raw coverage histogram wouldn't be correct either,
	 * since an interval can contain arbitrary configurations of other repeats.
	 *
	 * The procedure is supposed to work on the following temporary space:
	 * private static int[] coverageEvents, coverageDeltas;
	 * private static double[] coverageHistogram;
	 *
	 * that should have been allocated globally as follows, and that might be enlarged
	 * by the procedure:
	 * coverageEvents = new int[10000];  // Arbitrary, will be resized if necessary.
	 * coverageDeltas = new int[Reads.maxReadLength];
	 * coverageHistogram = new double[500];  // Arbitrary, will be resized if necessary.
	 */
	private static final void assignCoverage(int[] coverageEvents, int[] coverageDeltas, double[] coverageHistogram) {
		boolean currentType, newType;  // TRUE=dense, FALSE=alignment.
		int i, j, k;
		int lastEvent, first, last, pos, surface, length;
		int previousEvent, currentCoverage, maxCoverage, currentEvent;
		DenseSubstring currentDense, newDense;
		AlignmentInterval currentAlignment, newAlignment;
		
		Alignment.order=Alignment.IMPLYINGINTERVAL_STARTA;
		if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment);
		lastEvent=-1;
		first=Alignments.alignments[ReadA.sortedAlignments[0].id][3];
		last=Alignments.alignments[ReadA.sortedAlignments[0].id][4];
		currentType=false; currentDense=null; currentAlignment=null;
		if (ReadA.sortedAlignments[0].impliedByDenseSubstring!=null) {
			currentType=true;
			currentDense=ReadA.sortedAlignments[0].impliedByDenseSubstring;
			length=currentDense.endA-currentDense.startA+1;
			pos=Math.max(first-currentDense.startA,0);
			coverageEvents[++lastEvent]=pos; 
			coverageDeltas[pos]++;
			pos=Math.min(last-currentDense.startA,length-1);
			coverageEvents[++lastEvent]=pos; 
			coverageDeltas[pos]--;
		}
		else if (ReadA.sortedAlignments[0].inDenseSubstring!=null) {
			currentType=true;
			currentDense=ReadA.sortedAlignments[0].inDenseSubstring;
			length=currentDense.endA-currentDense.startA+1;
			pos=Math.max(first-currentDense.startA,0);
			coverageEvents[++lastEvent]=pos; 
			coverageDeltas[pos]++;
			pos=Math.min(last-currentDense.startA,length-1);
			coverageEvents[++lastEvent]=pos; 
			coverageDeltas[pos]--;
		}
		else if (ReadA.sortedAlignments[0].mergedToInterval!=null) {
			currentType=false;
			currentAlignment=ReadA.sortedAlignments[0].mergedToInterval;
			length=currentAlignment.lastPosition-currentAlignment.firstPosition+1;
			pos=Math.max(first-currentAlignment.firstPosition,0);
			coverageEvents[++lastEvent]=pos; 
			coverageDeltas[pos]++;
			pos=Math.min(last-currentAlignment.firstPosition,length-1);
			coverageEvents[++lastEvent]=pos; 
			coverageDeltas[pos]--;
		}
		else return;
		if (coverageEvents.length<((ReadA.lastSortedAlignment+1)<<1)) coverageEvents = new int[(ReadA.lastSortedAlignment+1)<<1];
		if (coverageHistogram.length<ReadA.lastSortedAlignment+1) coverageHistogram = new double[(ReadA.lastSortedAlignment+1)];
		Math.set(coverageDeltas,ReadA.readLength-1,0);
		for (i=1; i<=ReadA.lastSortedAlignment; i++) {
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null && ReadA.sortedAlignments[i].inDenseSubstring==null && ReadA.sortedAlignments[i].mergedToInterval==null) break;
			newType=false; newDense=null; newAlignment=null;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null) {
				newType=true;
				newDense=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			}
			else if (ReadA.sortedAlignments[i].inDenseSubstring!=null) {
				newType=true;
				newDense=ReadA.sortedAlignments[i].inDenseSubstring;
			}
			else if (ReadA.sortedAlignments[i].mergedToInterval!=null) {
				newType=false;
				newAlignment=ReadA.sortedAlignments[i].mergedToInterval;
			}
			if ( newType!=currentType || 
				 (newType && newDense.id!=currentDense.id) ||
				 (!newType && newAlignment.id!=currentAlignment.id)
		       ) {
				// Printing coverage histogram of the current interval
				if (lastEvent>0) Arrays.sort(coverageEvents,0,lastEvent+1);
				k=0; previousEvent=coverageEvents[0];
				for (j=1; j<=lastEvent; j++) {			
					if (coverageEvents[j]==previousEvent) continue;
					k++;
					if (k!=j) coverageEvents[k]=coverageEvents[j];
					previousEvent=coverageEvents[j];
				}
				lastEvent=k;
				previousEvent=0; currentEvent=Math.POSITIVE_INFINITY;
				currentCoverage=0; maxCoverage=0; surface=0;
				for (j=0; j<=lastEvent; j++) {
					currentEvent=coverageEvents[j];
					coverageHistogram[currentCoverage]+=currentEvent-previousEvent;
					surface+=currentEvent-previousEvent;
					currentCoverage+=coverageDeltas[currentEvent];
					if (currentCoverage>maxCoverage) maxCoverage=currentCoverage;
					previousEvent=currentEvent;
				}
				last=currentType?currentDense.endA-currentDense.startA:currentAlignment.lastPosition-currentAlignment.firstPosition;
				coverageHistogram[currentCoverage]+=last-currentEvent+1;
				surface+=last-currentEvent+1;
				if (IO.CONSISTENCY_CHECKS) {
					if (surface!=length) {
						System.err.println("assignCoverage> ERROR: wrong surface: "+surface+" != "+length);
						System.exit(1);
					}
					if (currentCoverage<0) {
						System.err.println("assignCoverage> ERROR: wrong coverage at the end: "+currentCoverage);
						System.exit(1);
					}
				}
				for (j=0; j<=maxCoverage; j++) coverageHistogram[j]/=surface;
				// The histogram is ready at this point
				
				// Next interval
				for (j=0; j<=lastEvent; j++) coverageDeltas[coverageEvents[j]]=0;
				Math.set(coverageHistogram,maxCoverage,0.0);
				currentType=newType; lastEvent=-1;
				if (newType) {
					currentDense=newDense;
					length=currentDense.endA-currentDense.startA+1;
				}
				else {
					currentAlignment=newAlignment;
					length=currentAlignment.lastPosition-currentAlignment.firstPosition+1;
				}
			}
			first=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			last=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (currentType) {
				pos=Math.max(first-currentDense.startA,0);
				coverageEvents[++lastEvent]=pos; 
				coverageDeltas[pos]++;
				pos=Math.min(last-currentDense.startA,length-1);
				coverageEvents[++lastEvent]=pos;
				coverageDeltas[pos]--;
			}
			else {
				pos=Math.max(first-currentAlignment.firstPosition,0);
				coverageEvents[++lastEvent]=pos; 
				coverageDeltas[pos]++;
				pos=Math.min(last-currentAlignment.firstPosition,length-1);
				coverageEvents[++lastEvent]=pos;
				coverageDeltas[pos]--;
			}
		}
		// Handling the final interval
		if (lastEvent>0) Arrays.sort(coverageEvents,0,lastEvent+1);
		k=0; previousEvent=coverageEvents[0];
		for (j=1; j<=lastEvent; j++) {			
			if (coverageEvents[j]==previousEvent) continue;
			k++;
			if (k!=j) coverageEvents[k]=coverageEvents[j];
			previousEvent=coverageEvents[j];
		}
		lastEvent=k;
		previousEvent=0; currentEvent=Math.POSITIVE_INFINITY;
		currentCoverage=0; maxCoverage=0; surface=0;
		for (j=0; j<=lastEvent; j++) {
			currentEvent=coverageEvents[j];
			coverageHistogram[currentCoverage]+=currentEvent-previousEvent;
			surface+=currentEvent-previousEvent;
			currentCoverage+=coverageDeltas[currentEvent];
			if (currentCoverage>maxCoverage) maxCoverage=currentCoverage;
			previousEvent=currentEvent;
		}
		last=currentType?currentDense.endA-currentDense.startA:currentAlignment.lastPosition-currentAlignment.firstPosition;
		coverageHistogram[currentCoverage]+=last-currentEvent+1;
		surface+=last-currentEvent+1;
		if (IO.CONSISTENCY_CHECKS) {
			if (surface!=length) {
				System.err.println("assignCoverage> ERROR: wrong surface: "+surface+" != "+length);
				System.exit(1);
			}
			if (currentCoverage<0) {
				System.err.println("assignCoverage> ERROR: wrong coverage at the end: "+currentCoverage);
				System.exit(1);
			}
		}
		for (j=0; j<=maxCoverage; j++) coverageHistogram[j]/=surface;
		// The histogram is ready at this point		
	}
	
	
	/**
	 * Detects the first/last position of all long random insertions, by fitting a
	 * regression tree on the quality histogram of a read (which contains just one value 
	 * every $Reads.QUALITY_SPACING$ positions).
	 *
	 * Remark: the procedure does not use positions that are strictly inside dense or 
	 * periodic intervals, since they are more likely to be small drops in quality than
	 * random insertions.
	 */
	public static final void detectRandomInsertions() {
		final int BIN_LENGTH = 1;
		final int MIN_INTERVAL_LENGTH = 2;  // Units of size $Reads.QUALITY_SPACING$
		final int MAX_DIFFERENCE = Reads.MAX_HIGH_QUALITY_SCORE;  // Arbitrary
		final int SHORT_LEAF_LENGTH = (Reads.QUALITY_SPACING)<<1;  // Arbitrary
		final int DISTANCE_THRESHOLD = IO.quantum;  // Arbitrary
		int i, j;
		int first, last, previousMinIntervalLength, nLocalMaxima;
		int nonMaximumLeaf, minimumLeaf, localMinimum, pos, length, to;
		double largestValue, sum;
		PeriodicSubstringInterval tmpPInterval = new PeriodicSubstringInterval();
		DenseSubstring tmpSubstring = new DenseSubstring();
		
		if (Reads.virtualQualities || Histograms.isConstant(Reads.getQualityArray(ReadA.id),0,Reads.getQualityArrayLength(ReadA.id)-1,BIN_LENGTH)) return;

		// Adding splits to $Factors$
		// Remark: all values in $Reads.qualities$ are nonzero.
		previousMinIntervalLength=RegressionTree.minIntervalLength;
		RegressionTree.minIntervalLength=MIN_INTERVAL_LENGTH;
		largestValue=Histograms.largestValue(Reads.getQualityArray(ReadA.id),0,Reads.getQualityArrayLength(ReadA.id)-1,MAX_DIFFERENCE);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("detectRandomInsertions> 1");
		nLocalMaxima=RegressionTree.buildRegressionTree(Reads.getQualityArray(ReadA.id),0,Reads.getQualityArrayLength(ReadA.id)-1,Reads.MAX_QUALITY_STD,1,largestValue,true,Alignments.minAlignmentLength);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("detectRandomInsertions> 2 nLocalMaxima="+nLocalMaxima);
		RegressionTree.minIntervalLength=previousMinIntervalLength;		
		if (nLocalMaxima>0) {
			for (i=0; i<=RegressionTree.lastLeaf; i++) {
				first=RegressionTree.leaves[i].firstPoint*Reads.QUALITY_SPACING;
				last=(RegressionTree.leaves[i].lastPoint+1)*Reads.QUALITY_SPACING;				
if (IO.SHOW_STD_ERR_PRIME) System.err.println("detectRandomInsertions> 2.5 for local maximum ["+first+".."+last+"]: "+(last-first+1<=SHORT_LEAF_LENGTH)+", "+Reads.isRandomInsertion(ReadA.id,first,last,true));
				if (RegressionTree.leaves[i].isLocalMaximum && (last-first+1<=SHORT_LEAF_LENGTH || Reads.isRandomInsertion(ReadA.id,first,last,true))) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("detectRandomInsertions> 3 considering local maximum ["+first+".."+last+"]");
					nonMaximumLeaf=-1; minimumLeaf=-1;
					for (j=i-1; j>=0; j--) {
						if (RegressionTree.leaves[j].isLocalMaximum) break;
						nonMaximumLeaf=j;
						if (RegressionTree.leaves[j].isLocalMinimum) {
							minimumLeaf=j;
							break;
						}
					}
					if (minimumLeaf!=-1) localMinimum=(int)(Histograms.getCenterOfMass(Reads.getQualityArray(ReadA.id),RegressionTree.leaves[minimumLeaf].lastPoint+1,RegressionTree.leaves[i].firstPoint-1)*Reads.QUALITY_SPACING);
					else if (nonMaximumLeaf!=-1) localMinimum=(int)(Histograms.getCenterOfMass(Reads.getQualityArray(ReadA.id),RegressionTree.leaves[nonMaximumLeaf].firstPoint,RegressionTree.leaves[i].firstPoint-1)*Reads.QUALITY_SPACING);
					else localMinimum=first;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("detectRandomInsertions> 4 localMinimum="+localMinimum);
					if ( localMinimum>IO.quantum && 
					     !PeriodicSubstrings.inPeriodicSubstringInterval(localMinimum,DISTANCE_THRESHOLD,tmpPInterval) &&
						 !DenseSubstrings.inDenseSubstring(localMinimum,DISTANCE_THRESHOLD,tmpSubstring)
					   ) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("detectRandomInsertions> 5");
						Factors.lastSplit++;
						Factors.splits[Factors.lastSplit].clear();
						Factors.splits[Factors.lastSplit].position=localMinimum;
						Factors.splits[Factors.lastSplit].nUnknown=1;
					}
					nonMaximumLeaf=-1; minimumLeaf=-1;
					for (j=i+1; j<=RegressionTree.lastLeaf; j++) {
						if (RegressionTree.leaves[j].isLocalMaximum) break;
						nonMaximumLeaf=j;
						if (RegressionTree.leaves[j].isLocalMinimum) {
							minimumLeaf=j;
							break;
						}
					}
					if (minimumLeaf!=-1) localMinimum=(int)(Histograms.getCenterOfMass(Reads.getQualityArray(ReadA.id),RegressionTree.leaves[i].lastPoint+1,RegressionTree.leaves[minimumLeaf].firstPoint-1)*Reads.QUALITY_SPACING);
					else if (nonMaximumLeaf!=-1) localMinimum=(int)(Histograms.getCenterOfMass(Reads.getQualityArray(ReadA.id),RegressionTree.leaves[i].lastPoint+1,RegressionTree.leaves[nonMaximumLeaf].lastPoint-1)*Reads.QUALITY_SPACING);
					else localMinimum=last;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("detectRandomInsertions> 6 localMinimum="+localMinimum);
					if ( localMinimum<ReadA.readLength-IO.quantum &&
						 !PeriodicSubstrings.inPeriodicSubstringInterval(localMinimum,DISTANCE_THRESHOLD,tmpPInterval) &&
						 !DenseSubstrings.inDenseSubstring(localMinimum,DISTANCE_THRESHOLD,tmpSubstring)
					   ) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("detectRandomInsertions> 7");
						Factors.lastSplit++;
						Factors.splits[Factors.lastSplit].clear();
						Factors.splits[Factors.lastSplit].position=localMinimum;
						Factors.splits[Factors.lastSplit].nUnknown=1;
					}
				}
			}
			
			// Searching for maximal low-quality substrings inside non-maximum leaves.
			for (i=0; i<=RegressionTree.lastLeaf; i++) {
				if (RegressionTree.leaves[i].isLocalMaximum) continue;
				detectRandomInsertions_simpleScan(RegressionTree.leaves[i].firstPoint,RegressionTree.leaves[i].lastPoint,MIN_INTERVAL_LENGTH,DISTANCE_THRESHOLD,tmpPInterval,tmpSubstring);
			}
		}
		else {
			// No local maximum, or no split. Searching for maximal low-quality 
			// substrings with a simple scan.
			detectRandomInsertions_simpleScan(0,Reads.getQualityArrayLength(ReadA.id)-1,MIN_INTERVAL_LENGTH,DISTANCE_THRESHOLD,tmpPInterval,tmpSubstring);
		}
	}
	
	
	/**
	 * @param first,last positions in $Reads.qualities[ReadA.id]$;
	 * @param window minimum length of a random insertion to be detected;
	 * @param tmpPInterval,tmpSubstring temporary space.
	 */
	private static final void detectRandomInsertions_simpleScan(int first, int last, int window, int distanceThreshold, PeriodicSubstringInterval tmpPInterval, DenseSubstring tmpSubstring) {
		int i, j;
		int to, pos;
		double sum;
		final double[] array = Reads.getQualityArray(ReadA.id);
		
		i=first; sum=0.0; to=Math.min(i+window-1,last);
		for (j=i; j<=to; j++) sum+=array[j];
		while (i<=last) {
			if (sum/(to-i+1)>=Reads.MIN_RANDOM_QUALITY_SCORE) {
				// Low-quality window
				pos=i*Reads.QUALITY_SPACING;
				if ( pos>IO.quantum && 
				     !PeriodicSubstrings.inPeriodicSubstringInterval(pos,distanceThreshold,tmpPInterval) &&
					 !DenseSubstrings.inDenseSubstring(pos,distanceThreshold,tmpSubstring)
				   ) {
					Factors.lastSplit++;
					Factors.splits[Factors.lastSplit].clear();
					Factors.splits[Factors.lastSplit].position=pos;
					Factors.splits[Factors.lastSplit].nUnknown=1;
				}
				j=i+window-1;
				while (j+1<=last && (sum+array[j+1])/(j-i+2)>=Reads.MIN_RANDOM_QUALITY_SCORE) {
					j++;
					sum+=array[j];
				}
				pos=(j+1)*Reads.QUALITY_SPACING-1;
				if ( pos<ReadA.readLength-IO.quantum && 
				     !PeriodicSubstrings.inPeriodicSubstringInterval(pos,distanceThreshold,tmpPInterval) &&
					 !DenseSubstrings.inDenseSubstring(pos,distanceThreshold,tmpSubstring)
				   ) {
					Factors.lastSplit++;
					Factors.splits[Factors.lastSplit].clear();
					Factors.splits[Factors.lastSplit].position=pos;
					Factors.splits[Factors.lastSplit].nUnknown=1;
				}
				i=j+2;
				sum=0.0;
				to=Math.min(i+window-1,last);
				for (j=i; j<=to; j++) sum+=array[j];
			}
			else {
				sum-=array[i];
				i++; to=Math.min(i+window-1,last);
				if (to-i+1==window) sum+=array[to];
			}
		}
	}
	
	
	/**
	 * Tries to use alignments that have not been assigned to any interval, in order to 
	 * transform an alignment interval or a dense substring to a more general replication 
	 * type. See $generalizeReplicationTypes_discardIntervals()$ for an explanation of 
	 * which intervals we consider for generalization.
	 *
	 * Remark: this procedure is intended to be used at the very end of factorization.
	 *
	 * Remark: at the end of the procedure, all intervals are sorted by first position.
	 */
	private static final void generalizeReplicationTypes() {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int ALIGNMENTS_THRESHOLD = 3*IO.coverage;  // Arbitrary
		final double SURFACE_THRESHOLD = 0.8;  // Arbitrary
		int i, j;
		int firstJForNextI, lastAlignment, nDiscardedDense, nDiscardedAlignment;
		int alignmentStartA, alignmentEndA, substringStartA, substringEndA, intervalStartA, intervalEndA;
		int nTransformedAlignmentIntervals;
		DenseSubstring tmpSubstring = new DenseSubstring();
		AlignmentInterval tmpInterval;
		int[] alignments;
		
		// Ensuring the necessary order in alignments and intervals
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		if (DenseSubstring.order!=DenseSubstring.STARTA) {
			DenseSubstring.order=DenseSubstring.STARTA;
			if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
		}
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (AlignmentIntervals.lastInterval>0) Arrays.sort(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1);
		}
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		}
		
		// Discarding dense substrings and alignment intervals
		generalizeReplicationTypes_discardIntervals(DISTANCE_THRESHOLD,tmpArray);
		nDiscardedDense=tmpArray[0];
		nDiscardedAlignment=tmpArray[1];
		
		// Generalizing dense substrings
		if (nDiscardedDense<DenseSubstrings.lastSubstring+1) {
			alignments=DenseSubstrings.stackPrime;
			i=0; j=0; firstJForNextI=-1; lastAlignment=-1;
			substringStartA=DenseSubstrings.substrings[i].startA;
			substringEndA=DenseSubstrings.substrings[i].endA;
			alignmentStartA=ReadA.sortedAlignments[j].startA();
			alignmentEndA=ReadA.sortedAlignments[j].endA();
			while (i<=DenseSubstrings.lastSubstring) {
				if (DenseSubstrings.substrings[i].discarded) {
					i++;
					substringStartA=DenseSubstrings.substrings[i].startA;
					substringEndA=DenseSubstrings.substrings[i].endA;
					continue;
				}
				if (j>ReadA.lastSortedAlignment || alignmentStartA>=substringEndA) {
					if (lastAlignment>=0) {
						generalizeReplicationTypes_getAlignmentStats(alignments,lastAlignment,substringStartA,substringEndA,DISTANCE_THRESHOLD,tmpArray);
						if (DenseSubstrings.substrings[i].generalizeType(tmpArray,(lastAlignment+1)/5,DISTANCE_THRESHOLD,ALIGNMENTS_THRESHOLD,SURFACE_THRESHOLD,tmpArrayPrime)) {
							generalize(DenseSubstrings.substrings[i],alignments,lastAlignment);
						}
					}
					i++;
					substringStartA=DenseSubstrings.substrings[i].startA;
					substringEndA=DenseSubstrings.substrings[i].endA;
					if (firstJForNextI!=-1) {
						j=firstJForNextI;
						alignmentStartA=ReadA.sortedAlignments[j].startA();
						alignmentEndA=ReadA.sortedAlignments[j].endA();
					}
					firstJForNextI=-1;
					lastAlignment=-1;
					continue;
				}
				if (alignmentEndA<=substringStartA || ReadA.sortedAlignments[j].isImplied()) {
					j++;
					if (j<=ReadA.lastSortedAlignment) {
						alignmentStartA=ReadA.sortedAlignments[j].startA();
						alignmentEndA=ReadA.sortedAlignments[j].endA();
					}
					continue;
				}
				if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && alignmentEndA>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
				if ( ( Intervals.isApproximatelyContained_lowQuality(alignmentStartA,alignmentEndA,substringStartA,substringEndA,ReadA.id) ||
					   Intervals.areApproximatelyIdentical_lowQuality(alignmentStartA,alignmentEndA,substringStartA,substringEndA,ReadA.id)
					 )  &&
					 ( (alignmentStartA>substringStartA+DISTANCE_THRESHOLD && ReadA.sortedAlignments[j].isLeftMaximalB==1) ||
					   (alignmentEndA<substringEndA-DISTANCE_THRESHOLD && ReadA.sortedAlignments[j].isRightMaximalB==1)
					 )
				   ) {
					alignments[++lastAlignment]=Math.max(alignmentStartA,substringStartA);
					alignments[++lastAlignment]=Math.min(alignmentEndA,substringEndA);
					alignments[++lastAlignment]=alignmentStartA<substringStartA-DISTANCE_THRESHOLD?0:ReadA.sortedAlignments[j].isLeftMaximalB;
					alignments[++lastAlignment]=alignmentEndA>substringEndA+DISTANCE_THRESHOLD?0:ReadA.sortedAlignments[j].isRightMaximalB;
					alignments[++lastAlignment]=j;
				}
				j++;
				if (j<=ReadA.lastSortedAlignment) {
					alignmentStartA=ReadA.sortedAlignments[j].startA();
					alignmentEndA=ReadA.sortedAlignments[j].endA();
				}
			}
		}
		
		// Generalizing alignment intervals
		nTransformedAlignmentIntervals=0;
		if (nDiscardedAlignment<AlignmentIntervals.lastInterval+1) {
			alignments=AlignmentIntervals.stack;
			i=0; j=0; firstJForNextI=-1; lastAlignment=-1;
			intervalStartA=AlignmentIntervals.intervals[i].firstPosition;
			intervalEndA=AlignmentIntervals.intervals[i].lastPosition;
			alignmentStartA=ReadA.sortedAlignments[j].startA();
			alignmentEndA=ReadA.sortedAlignments[j].endA();
			while (i<=AlignmentIntervals.lastInterval) {
				if (AlignmentIntervals.intervals[i].discarded) {
					i++;
					intervalStartA=AlignmentIntervals.intervals[i].firstPosition;
					intervalEndA=AlignmentIntervals.intervals[i].lastPosition;
					continue;
				}
				if (j>ReadA.lastSortedAlignment || alignmentStartA>=intervalEndA) {
					if (lastAlignment>=0) {
						generalizeReplicationTypes_getAlignmentStats(alignments,lastAlignment,intervalStartA,intervalEndA,DISTANCE_THRESHOLD,tmpArray);
						if (AlignmentIntervals.intervals[i].generalizeType(tmpArray,(lastAlignment+1)/5,DISTANCE_THRESHOLD,ALIGNMENTS_THRESHOLD,SURFACE_THRESHOLD,tmpSubstring,tmpArrayPrime)) {
							generalize(AlignmentIntervals.intervals[i],alignments,lastAlignment,tmpSubstring);
							nTransformedAlignmentIntervals++;
							AlignmentIntervals.intervals[i].discarded=false;
						}
						else AlignmentIntervals.intervals[i].discarded=true;
					}
					else AlignmentIntervals.intervals[i].discarded=true;
					i++;
					intervalStartA=AlignmentIntervals.intervals[i].firstPosition;
					intervalEndA=AlignmentIntervals.intervals[i].lastPosition;
					if (firstJForNextI!=-1) {
						j=firstJForNextI;
						alignmentStartA=ReadA.sortedAlignments[j].startA();
						alignmentEndA=ReadA.sortedAlignments[j].endA();
					}
					firstJForNextI=-1;
					lastAlignment=-1;
					continue;
				}
				if (alignmentEndA<=intervalStartA || ReadA.sortedAlignments[j].isImplied()) {
					j++;
					if (j<=ReadA.lastSortedAlignment) {
						alignmentStartA=ReadA.sortedAlignments[j].startA();
						alignmentEndA=ReadA.sortedAlignments[j].endA();
					}
					continue;
				}
				if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && alignmentEndA>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
				if ( ( Intervals.isApproximatelyContained_lowQuality(alignmentStartA,alignmentEndA,intervalStartA,intervalEndA,ReadA.id) ||
					   Intervals.areApproximatelyIdentical_lowQuality(alignmentStartA,alignmentEndA,intervalStartA,intervalEndA,ReadA.id)
					 ) &&
					 ( (alignmentStartA>intervalStartA+DISTANCE_THRESHOLD && ReadA.sortedAlignments[j].isLeftMaximalB==1) ||
					   (alignmentEndA<intervalEndA-DISTANCE_THRESHOLD && ReadA.sortedAlignments[j].isRightMaximalB==1)
					 )
				   ) {
					alignments[++lastAlignment]=Math.max(alignmentStartA,intervalStartA);
					alignments[++lastAlignment]=Math.min(alignmentEndA,intervalEndA);
					alignments[++lastAlignment]=alignmentStartA<intervalStartA-DISTANCE_THRESHOLD?0:ReadA.sortedAlignments[j].isLeftMaximalB;
					alignments[++lastAlignment]=alignmentEndA>intervalEndA+DISTANCE_THRESHOLD?0:ReadA.sortedAlignments[j].isRightMaximalB;
					alignments[++lastAlignment]=j;
				}
				j++;
				if (j<=ReadA.lastSortedAlignment) {
					alignmentStartA=ReadA.sortedAlignments[j].startA();
					alignmentEndA=ReadA.sortedAlignments[j].endA();
				}
			}
		}
		generalizeReplicationTypes_updatePointersFromAlignments();
		
		// Removing transformed intervals from $AlignmentIntervals$, and sorting
		// $DenseSubstrings$.
		if (nTransformedAlignmentIntervals>0) {
			j=-1;
			for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
				if (!AlignmentIntervals.intervals[i].discarded) continue;
				j++;
				if (j!=i) {
					tmpInterval=AlignmentIntervals.intervals[j];
					AlignmentIntervals.intervals[j]=AlignmentIntervals.intervals[i];
					AlignmentIntervals.intervals[i]=tmpInterval;
				}
			}
			AlignmentIntervals.lastInterval=j;
			DenseSubstring.order=DenseSubstring.STARTA;
			if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
		}
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].discarded=false;
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) AlignmentIntervals.intervals[i].discarded=false;
	}


	/** 
	 * Marks dense substrings and alignment intervals whose type should not be 
	 * generalized. For simplicity, we choose not to generalize any interval that 
	 * straddles, contains, or is contained in another interval of any type. Otherwise, 
	 * changing the type of an interval might require concatenations, or checking 
	 * containment and merging, whose logic belongs to previous steps of the pipeline and 
	 * would be cumbersome to replicate here.
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings$, $AlignmentIntervals.
	 * intervals$ and $PeriodicSubstrings.intervals$ to be sorted by their first position.
	 *
	 * @param tmpArray output values: 0=n. discarded dense substrings; 1=n. discarded 
	 * intervals.
	 */
	private static final void generalizeReplicationTypes_discardIntervals(int distanceThreshold, int[] tmpArray) {
		int i, j;
		int firstJForNextI, nDiscarded, currentStartA, currentEndA, otherStartA, otherEndA;
		
		nDiscarded=0;
		if (DenseSubstrings.lastSubstring>=0) {
			// Dense substrings -> Dense substrings
			for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
				if (DenseSubstrings.substrings[i].substringReplication) {
					DenseSubstrings.substrings[i].discarded=true;
					nDiscarded++;
					continue;
				}
				if (i<DenseSubstrings.lastSubstring && DenseSubstrings.substrings[i+1].startA<DenseSubstrings.substrings[i].endA-distanceThreshold) {
					DenseSubstrings.substrings[i].discarded=true;
					nDiscarded++;
					continue;
				}
				DenseSubstrings.substrings[i].discarded=false;
				for (j=i-1; j>=0; j--) {
					if (DenseSubstrings.substrings[j].endA<=DenseSubstrings.substrings[i].startA+distanceThreshold) continue;
					DenseSubstrings.substrings[i].discarded=true;
					nDiscarded++;
					break;
				}
			}
			// Dense substrings -> Alignment intervals
			if (nDiscarded<DenseSubstrings.lastSubstring+1 && AlignmentIntervals.lastInterval>=0) {
				i=0; j=0; firstJForNextI=-1;
				currentStartA=DenseSubstrings.substrings[i].startA;
				currentEndA=DenseSubstrings.substrings[i].endA;
				otherStartA=AlignmentIntervals.intervals[j].firstPosition;
				otherEndA=AlignmentIntervals.intervals[j].lastPosition;
				while (i<=DenseSubstrings.lastSubstring) {
					if (DenseSubstrings.substrings[i].discarded || j>AlignmentIntervals.lastInterval || otherStartA>=currentEndA) {
						i++;
						currentStartA=DenseSubstrings.substrings[i].startA;
						currentEndA=DenseSubstrings.substrings[i].endA;
						if (firstJForNextI!=-1) {
							j=firstJForNextI;
							otherStartA=AlignmentIntervals.intervals[j].firstPosition;
							otherEndA=AlignmentIntervals.intervals[j].lastPosition;
						}
						firstJForNextI=-1;
						continue;
					}
					if (otherEndA<currentStartA) {
						j++;
						otherStartA=AlignmentIntervals.intervals[j].firstPosition;
						otherEndA=AlignmentIntervals.intervals[j].lastPosition;
						continue;
					}
					if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && otherEndA>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
					if ( (otherEndA>currentStartA+distanceThreshold && otherStartA<=currentEndA) || 
						 (otherStartA<currentEndA-distanceThreshold && otherEndA>=currentStartA)
					   ) {
						DenseSubstrings.substrings[i].discarded=true;
						nDiscarded++;
					}
					j++;
					otherStartA=AlignmentIntervals.intervals[j].firstPosition;
					otherEndA=AlignmentIntervals.intervals[j].lastPosition;
				}
			}
			// Dense substrings -> Periodic substring intervals
			if (nDiscarded<DenseSubstrings.lastSubstring+1 && PeriodicSubstrings.lastInterval>=0) {
				i=0; j=0; firstJForNextI=-1;
				currentStartA=DenseSubstrings.substrings[i].startA;
				currentEndA=DenseSubstrings.substrings[i].endA;
				otherStartA=PeriodicSubstrings.intervals[j].firstPosition;
				otherEndA=PeriodicSubstrings.intervals[j].lastPosition;
				while (i<=DenseSubstrings.lastSubstring) {
					if (DenseSubstrings.substrings[i].discarded || j>PeriodicSubstrings.lastInterval || otherStartA>=currentEndA) {
						i++;
						currentStartA=DenseSubstrings.substrings[i].startA;
						currentEndA=DenseSubstrings.substrings[i].endA;
						if (firstJForNextI!=-1) {
							j=firstJForNextI;
							otherStartA=PeriodicSubstrings.intervals[j].firstPosition;
							otherEndA=PeriodicSubstrings.intervals[j].lastPosition;
						}
						firstJForNextI=-1;
						continue;
					}
					if (otherEndA<currentStartA) {
						j++;
						otherStartA=PeriodicSubstrings.intervals[j].firstPosition;
						otherEndA=PeriodicSubstrings.intervals[j].lastPosition;
						continue;
					}
					if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && otherEndA>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
					if ( (otherEndA>currentStartA+distanceThreshold && otherStartA<=currentEndA) || 
						 (otherStartA<currentEndA-distanceThreshold && otherEndA>=currentStartA)
					   ) {
						DenseSubstrings.substrings[i].discarded=true;
						nDiscarded++;
					}
					j++;
					otherStartA=PeriodicSubstrings.intervals[j].firstPosition;
					otherEndA=PeriodicSubstrings.intervals[j].lastPosition;
				}
			}
		}
		tmpArray[0]=nDiscarded;
		
		nDiscarded=0;
		if (AlignmentIntervals.lastInterval>=0) {
			// Alignment intervals -> Alignment intervals
			for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
				if (i<AlignmentIntervals.lastInterval && AlignmentIntervals.intervals[i+1].firstPosition<AlignmentIntervals.intervals[i].lastPosition-distanceThreshold) {
					AlignmentIntervals.intervals[i].discarded=true;
					nDiscarded++;
					continue;
				}
				AlignmentIntervals.intervals[i].discarded=false;
				for (j=i-1; j>=0; j--) {
					if (AlignmentIntervals.intervals[j].lastPosition<=AlignmentIntervals.intervals[i].firstPosition+distanceThreshold) continue;
					AlignmentIntervals.intervals[i].discarded=true;
					nDiscarded++;
					break;
				}
			}
			// Alignment intervals -> Dense substrings
			if (nDiscarded<AlignmentIntervals.lastInterval+1 && DenseSubstrings.lastSubstring>=0) {
				i=0; j=0; firstJForNextI=-1;
				currentStartA=AlignmentIntervals.intervals[i].firstPosition;
				currentEndA=AlignmentIntervals.intervals[i].lastPosition;
				otherStartA=DenseSubstrings.substrings[j].startA;
				otherEndA=DenseSubstrings.substrings[j].endA;
				while (i<=AlignmentIntervals.lastInterval) {
					if (AlignmentIntervals.intervals[i].discarded || j>DenseSubstrings.lastSubstring || otherStartA>=currentEndA) {
						i++;
						currentStartA=AlignmentIntervals.intervals[i].firstPosition;
						currentEndA=AlignmentIntervals.intervals[i].lastPosition;
						if (firstJForNextI!=-1) {
							j=firstJForNextI;
							otherStartA=DenseSubstrings.substrings[j].startA;
							otherEndA=DenseSubstrings.substrings[j].endA;
						}
						firstJForNextI=-1;
						continue;
					}
					if (otherEndA<currentStartA) {
						j++;
						otherStartA=DenseSubstrings.substrings[j].startA;
						otherEndA=DenseSubstrings.substrings[j].endA;
						continue;
					}
					if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && otherEndA>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
					if ( (otherEndA>currentStartA+distanceThreshold && otherStartA<=currentEndA) || 
						 (otherStartA<currentEndA-distanceThreshold && otherEndA>=currentStartA)
					   ) {
						AlignmentIntervals.intervals[i].discarded=true;
						nDiscarded++;
					}
					j++;
					otherStartA=DenseSubstrings.substrings[j].startA;
					otherEndA=DenseSubstrings.substrings[j].endA;
				}
			}
			// Alignment intervals -> Periodic substring intervals
			if (nDiscarded<AlignmentIntervals.lastInterval+1 && PeriodicSubstrings.lastInterval>=0) {
				i=0; j=0; firstJForNextI=-1;
				currentStartA=AlignmentIntervals.intervals[i].firstPosition;
				currentEndA=AlignmentIntervals.intervals[i].lastPosition;
				otherStartA=PeriodicSubstrings.intervals[j].firstPosition;
				otherEndA=PeriodicSubstrings.intervals[j].lastPosition;
				while (i<=AlignmentIntervals.lastInterval) {
					if (AlignmentIntervals.intervals[i].discarded || j>PeriodicSubstrings.lastInterval || otherStartA>=currentEndA) {
						i++;
						currentStartA=AlignmentIntervals.intervals[i].firstPosition;
						currentEndA=AlignmentIntervals.intervals[i].lastPosition;
						if (firstJForNextI!=-1) {
							j=firstJForNextI;
							otherStartA=PeriodicSubstrings.intervals[j].firstPosition;
							otherEndA=PeriodicSubstrings.intervals[j].lastPosition;
						}
						firstJForNextI=-1;
						continue;
					}
					if (otherEndA<currentStartA) {
						j++;
						otherStartA=PeriodicSubstrings.intervals[j].firstPosition;
						otherEndA=PeriodicSubstrings.intervals[j].lastPosition;
						continue;
					}
					if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && otherEndA>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
					if ( (otherEndA>currentStartA+distanceThreshold && otherStartA<=currentEndA) || 
						 (otherStartA<currentEndA-distanceThreshold && otherEndA>=currentStartA)
					   ) {
						AlignmentIntervals.intervals[i].discarded=true;
						nDiscarded++;
					}
					j++;
					otherStartA=PeriodicSubstrings.intervals[j].firstPosition;
					otherEndA=PeriodicSubstrings.intervals[j].lastPosition;
				}
			}
		}
		tmpArray[1]=nDiscarded;
	}

	
	/**
	 * Stores in $out$ the following statistics on $alignments[0..lastAlignment]$:
	 *
	 * 0: n. left-maximal alignments;
	 * 1: n. right-maximal alignments;
	 * 2: n. alignments that are both left- and right-maximal;
	 * 3: 1 iff no left-maximal alignment is right-maximal and far from $endA$;
	 * 4: 1 iff no right-maximal alignment is left-maximal and far from $startA$;
	 * 5: min length of a prefix alignment;
	 * 6: min length of a suffix alignment;
	 * 7: min length of an alignment that is both left- and right-maximal;
	 * 8: leftmost start of a left-maximal alignment;
	 * 9: rightmost end of a right-maximal alignment;
	 * 10: surface of the union of all alignments.
	 *
	 * @param $alignments[0..lastAlignment]$ consists of blocks of 5 elements: 
	 * $alignmentStart, alignmentEnd, isLeftMaximalB, isRightMaximalB, alignmentID$;
	 * all alignments are assumed to be strictly contained inside an interval;
	 * @param startA,endA first and last positions of the interval;
	 * @param distanceThreshold threshold on distance identity.
	 */
	private static final void generalizeReplicationTypes_getAlignmentStats(int[] alignments, int lastAlignment, int startA, int endA, int distanceThreshold, int[] out) {
		boolean leftMaximalAligned, rightMaximalAligned;
		int i;
		int nLeftMaximal, nRightMaximal, nLeftRightMaximal;
		int minPrefLength, minSufLength, minSubstringLength;
		int firstMaxStart, lastMaxEnd, blockStart, blockEnd, surface, length;

		nLeftMaximal=0; nRightMaximal=0; nLeftRightMaximal=0;
		leftMaximalAligned=true; rightMaximalAligned=true;
		minPrefLength=Math.POSITIVE_INFINITY;
		minSufLength=Math.POSITIVE_INFINITY;
		minSubstringLength=Math.POSITIVE_INFINITY;
		firstMaxStart=Math.POSITIVE_INFINITY; lastMaxEnd=0;
		surface=0; blockStart=-1; blockEnd=-1;
		for (i=0; i<=lastAlignment; i+=5) {
			if (alignments[i+2]==1 && alignments[i]>startA+distanceThreshold && alignments[i]<=endA) {
				nLeftMaximal++;
				if (alignments[i+3]==1 && alignments[i+1]<endA-distanceThreshold) leftMaximalAligned=false;
				else {
					length=endA-alignments[i]+1;
					minSufLength=Math.min(minSufLength,length);
				}
				firstMaxStart=Math.min(alignments[i],firstMaxStart);
			}
			if (alignments[i+3]==1 && alignments[i+1]<endA-distanceThreshold && alignments[i+1]>=startA) {
				nRightMaximal++;
				if (alignments[i+2]==1 && alignments[i]>startA+distanceThreshold) rightMaximalAligned=false;
				else {
					length=alignments[i+1]-startA+1;
					minPrefLength=Math.min(minPrefLength,length);
				}
				lastMaxEnd=Math.max(alignments[i+1],lastMaxEnd);
			}
			if (alignments[i+2]==1 && alignments[i]>startA+distanceThreshold && alignments[i+3]==1 && alignments[i+1]<endA-distanceThreshold) {
				nLeftRightMaximal++;
				length=alignments[i+1]-alignments[i]+1;
				minSubstringLength=Math.min(length,minSubstringLength);
			}
			if (alignments[i]<blockEnd-distanceThreshold) blockEnd=Math.max(blockEnd,alignments[i+1]);
			else {
				if (blockStart!=-1 && blockEnd!=-1) surface+=blockEnd-blockStart+1;
				blockStart=alignments[i]; blockEnd=alignments[i+1];
			}
		}
		if (blockStart!=-1 && blockEnd!=-1) surface+=blockEnd-blockStart+1;  // Last block
		
		// Outputting
		out[0]=nLeftMaximal;
		out[1]=nRightMaximal;
		out[2]=nLeftRightMaximal;
		out[3]=leftMaximalAligned?1:0;
		out[4]=rightMaximalAligned?1:0;
		out[5]=minPrefLength;
		out[6]=minSufLength;
		out[7]=minSubstringLength;
		out[8]=firstMaxStart;
		out[9]=lastMaxEnd;
		out[10]=surface;
	}
	
	
	private static final void generalize(DenseSubstring substring, int[] alignments, int lastAlignment) {
		int i;
		Alignment alignment;
		
		for (i=0; i<=lastAlignment; i+=5) {
			alignment=ReadA.sortedAlignments[alignments[i+4]];
			if (substring.substringReplication) {
				alignment.inDenseSubstring=substring;
				alignment.impliedByDenseSubstring=null;
				alignment.periodicSubstringInterval=null;
				alignment.mergedToInterval=null;
			}
			else {
				alignment.impliedByDenseSubstring=substring;
				alignment.inDenseSubstring=null;
				alignment.periodicSubstringInterval=null;
				alignment.mergedToInterval=null;
			}
		}
	}
	
	
	private static final void generalize(AlignmentInterval interval, int[] alignments, int lastAlignment, DenseSubstring tmpSubstring) {
		int i;
		Alignment alignment;
		
		for (i=0; i<=lastAlignment; i+=5) {
			alignment=ReadA.sortedAlignments[alignments[i+4]];
			alignment.inDenseSubstring=null;
			alignment.impliedByDenseSubstring=null;
			alignment.periodicSubstringInterval=null;
			alignment.mergedToInterval=interval;
		}
		DenseSubstrings.lastSubstring++;
		DenseSubstrings.substrings[DenseSubstrings.lastSubstring].clone(tmpSubstring);
		interval.denseSubstring=DenseSubstrings.substrings[DenseSubstrings.lastSubstring];
	}
	
	
	private static final void generalizeReplicationTypes_updatePointersFromAlignments() {
		int i;
		DenseSubstring substring;

		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null && ReadA.sortedAlignments[i].impliedByDenseSubstring.substringReplication) {
				ReadA.sortedAlignments[i].inDenseSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
			}
			if (ReadA.sortedAlignments[i].mergedToInterval!=null && !ReadA.sortedAlignments[i].mergedToInterval.discarded) {
				substring=ReadA.sortedAlignments[i].mergedToInterval.denseSubstring;
				if (substring.substringReplication) {
					ReadA.sortedAlignments[i].inDenseSubstring=substring;
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				}
				else {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=substring;
					ReadA.sortedAlignments[i].inDenseSubstring=null;
				}
				ReadA.sortedAlignments[i].mergedToInterval=null;
			}
		}
	}	

	
	/**
	 * It can happen that an alignment interval contains exactly one dense substring, 
	 * intersects no other interval, and the boundaries of the interval and of the dense
	 * substring are similar. The procedure discards such alignment interval, and assigns
	 * its alignments to the dense substring.
	 *
	 * Remark: this procedure is intended to be used at the very end of factorization.
	 *
	 * Remark: the procedure assumes that all interval lists and $ReadA.sortedAlignments$ 
	 * are sorted by first position.
	 */
	private static final void generalizeReplicationTypes_alignmentIntervals() {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int DISTANCE_THRESHOLD = IO.quantum<<1;  // Arbitrary
		boolean found;
		int i, j;
		int firstJForNextI, firstPosition, lastPosition;
		AlignmentInterval tmpInterval;
		DenseSubstring substring;
		if (AlignmentIntervals.lastInterval<0) return;
		
		// Counting intersections
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			AlignmentIntervals.intervals[i].tmpInt1=0;
			AlignmentIntervals.intervals[i].flag1=false;
			AlignmentIntervals.intervals[i].denseSubstring=null;
		}
		if (PeriodicSubstrings.lastInterval>=0) {
			i=0; j=0; firstJForNextI=-1;
			firstPosition=AlignmentIntervals.intervals[0].firstPosition;
			lastPosition=AlignmentIntervals.intervals[0].lastPosition;
			while (i<=AlignmentIntervals.lastInterval) {
				if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>lastPosition) {
					i++;
					firstPosition=AlignmentIntervals.intervals[i].firstPosition;
					lastPosition=AlignmentIntervals.intervals[i].lastPosition;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (PeriodicSubstrings.intervals[j].lastPosition<firstPosition) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && PeriodicSubstrings.intervals[j].lastPosition>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
				if (Intervals.intersect(firstPosition,lastPosition,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition)) AlignmentIntervals.intervals[i].tmpInt1++;
				j++;
			}
		}
		if (DenseSubstrings.lastSubstring>=0) {
			i=0; j=0; firstJForNextI=-1;
			firstPosition=AlignmentIntervals.intervals[0].firstPosition;
			lastPosition=AlignmentIntervals.intervals[0].lastPosition;
			while (i<=AlignmentIntervals.lastInterval) {
				if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>lastPosition) {
					i++;
					firstPosition=AlignmentIntervals.intervals[i].firstPosition;
					lastPosition=AlignmentIntervals.intervals[i].lastPosition;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (DenseSubstrings.substrings[j].endA<firstPosition) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && DenseSubstrings.substrings[j].endA>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
				if (Intervals.intersect(firstPosition,lastPosition,DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA)) {
					AlignmentIntervals.intervals[i].tmpInt1++;
					AlignmentIntervals.intervals[i].flag1=true;
				}
				j++;
			}
		}
		if (AlignmentIntervals.lastInterval>=0) {
			i=0; j=0; firstJForNextI=-1;
			firstPosition=AlignmentIntervals.intervals[0].firstPosition;
			lastPosition=AlignmentIntervals.intervals[0].lastPosition;
			while (i<=AlignmentIntervals.lastInterval) {
				if (j>AlignmentIntervals.lastInterval || AlignmentIntervals.intervals[j].firstPosition>lastPosition) {
					i++;
					firstPosition=AlignmentIntervals.intervals[i].firstPosition;
					lastPosition=AlignmentIntervals.intervals[i].lastPosition;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (j==i || AlignmentIntervals.intervals[j].lastPosition<firstPosition) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && AlignmentIntervals.intervals[j].lastPosition>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
				if (Intervals.intersect(firstPosition,lastPosition,AlignmentIntervals.intervals[j].firstPosition,AlignmentIntervals.intervals[j].lastPosition)) AlignmentIntervals.intervals[i].tmpInt1++;
				j++;
			}
		}
		
		// Assigning alignment intervals to their only intersecting dense substring
		found=false;
		if (DenseSubstrings.lastSubstring>=0) {
			i=0; j=0; firstJForNextI=-1;
			firstPosition=AlignmentIntervals.intervals[0].firstPosition;
			lastPosition=AlignmentIntervals.intervals[0].lastPosition;
			while (i<=AlignmentIntervals.lastInterval) {				
				if (AlignmentIntervals.intervals[i].tmpInt1!=1 || !AlignmentIntervals.intervals[i].flag1) {
					i++;
					firstPosition=AlignmentIntervals.intervals[i].firstPosition;
					lastPosition=AlignmentIntervals.intervals[i].lastPosition;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>lastPosition) {
					i++;
					firstPosition=AlignmentIntervals.intervals[i].firstPosition;
					lastPosition=AlignmentIntervals.intervals[i].lastPosition;
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (DenseSubstrings.substrings[j].endA<firstPosition) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<AlignmentIntervals.lastInterval && DenseSubstrings.substrings[j].endA>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextI=j;
				if (!Intervals.intersect(firstPosition,lastPosition,DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA)) {
					j++;
					continue;
				}
				if ( ( DenseSubstrings.substrings[j].prefixReplication && !DenseSubstrings.substrings[j].suffixReplication && !DenseSubstrings.substrings[j].substringReplication && !DenseSubstrings.substrings[j].singleDeletionReplication &&
					   Math.abs(DenseSubstrings.substrings[j].startA,firstPosition)<=IDENTITY_THRESHOLD && 
					   DenseSubstrings.substrings[j].endA<=lastPosition && DenseSubstrings.substrings[j].endA>=lastPosition-DISTANCE_THRESHOLD
					 ) ||
				     ( !DenseSubstrings.substrings[j].prefixReplication && DenseSubstrings.substrings[j].suffixReplication && !DenseSubstrings.substrings[j].substringReplication && !DenseSubstrings.substrings[j].singleDeletionReplication &&
					   Math.abs(DenseSubstrings.substrings[j].endA,lastPosition)<=IDENTITY_THRESHOLD && 
					   DenseSubstrings.substrings[j].startA>=firstPosition && DenseSubstrings.substrings[j].startA<=firstPosition+DISTANCE_THRESHOLD
					 ) ||
				     ( (DenseSubstrings.substrings[j].substringReplication || DenseSubstrings.substrings[j].singleDeletionReplication) &&
					   (DenseSubstrings.substrings[j].startA>=firstPosition-IDENTITY_THRESHOLD && DenseSubstrings.substrings[j].startA<=firstPosition+DISTANCE_THRESHOLD) && 
					   (DenseSubstrings.substrings[j].endA<=lastPosition+IDENTITY_THRESHOLD && DenseSubstrings.substrings[j].endA>=lastPosition-DISTANCE_THRESHOLD)
				     )
				   ) {
					   AlignmentIntervals.intervals[i].denseSubstring=DenseSubstrings.substrings[j];
					   found=true;
				   }
				j++;
			}
		}
		if (!found) return;
		
		// Compacting intervals
		j=-1;
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			if (AlignmentIntervals.intervals[i].denseSubstring!=null) continue;
			j++;
			if (j!=i) {
				tmpInterval=AlignmentIntervals.intervals[j];
				AlignmentIntervals.intervals[j]=AlignmentIntervals.intervals[i];
				AlignmentIntervals.intervals[i]=tmpInterval;
			}
		}
		AlignmentIntervals.lastInterval=j;
		
		// Resetting pointers from alignments
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			tmpInterval=ReadA.sortedAlignments[i].mergedToInterval;
			substring=tmpInterval==null?null:tmpInterval.denseSubstring;
			if (substring!=null) {
				if (substring.substringReplication) {
					ReadA.sortedAlignments[i].inDenseSubstring=substring;
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				}
				else {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=substring;
					ReadA.sortedAlignments[i].inDenseSubstring=null;
				}
				ReadA.sortedAlignments[i].mergedToInterval=null;
			}
		}
		
		// Updating boundaries of dense substrings
		DenseSubstrings.fixBoundariesOfSubstringType(tmpArray);
	}
	
	
	
	
	
	
	
	
	// --------------------------- STATISTICS FOR ASSEMBLY -------------------------------
	
	/**
	 * Cumulates to $out$ some statistics that might be useful for assembly, computed on 
	 * the current read. $out$ has the following format:
	 *
	 * 0: number of factorized reads;
	 * 1: sum of lengths of all high-quality substrings of factorized reads;
	 * 2: sum of read positions that belong to a short-period interval;
 	 * 3: number of factorized reads whose majority of high-quality positions belong to
     * short-period intervals;
	 * 4: sum of read positions that belong to a long-period interval;
 	 * 5: number of factorized reads whose majority of high-quality positions belong to
     * long-period intervals;
	 * 6: sum of read positions that belong to any interval (might be greater than cell 1 
	 * since intervals might contain low-quality regions);
	 * 7: number of factorized reads whose majority of high-quality positions belong to
	 * any interval;
	 * 8: number of reads in cell 7 that might contain distinct repeat modules;
	 * 9: number of collisions;
	 * 10: number of factorized reads whose prefix or suffix belong to an interval;
	 * 11: not altered
	 * 12: not altered
	 * 13: number of times a periodic interval is flanked by a high-quality substring that
	 * is not covered by any repeat;
	 * 14: number of times we see configuration XYZ, where X,Z are parts of a read covered
	 * by periodic intervals, and Y is fully covered by intervals that are not periodic.
	 *
	 * The term "majority" in cells 3,5,7 means all high-quality positions, except for a
	 * small constant value. In cell 8, distinct modules are assumed to exist in a read
	 * whenever there is a position $i$ such that intervals either end at $i$ or start at 
	 * $i+1$, i.e. no interval contains both positions.
	 */
	private static final void assemblyStatistics(long[] out) {
		final int MAX_HOLE = (Alignments.minAlignmentLength<<1)/3;  // Arbitrary
		final int GROWTH_RATE = 100;  // Arbitrary
		int i;
		int sum, start, end, lastWindow, highQualitySurface, nIntervals, tag;
		PeriodicSubstringInterval pInterval;
		
		out[0]++;
		highQualitySurface=(int)(ReadA.readLength*Histograms.getFraction(Reads.getQualityArray(ReadA.id),0,Reads.getQualityArrayLength(ReadA.id)-1,Reads.MAX_HIGH_QUALITY_SCORE+1,false));
		out[1]+=highQualitySurface;
		
		// Surface covered by short-period periodic
		sum=0; start=-1; end=-1;
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			pInterval=PeriodicSubstrings.intervals[i];
			if (pInterval.hasLongPeriod) continue;
			if (start==-1) {
				start=pInterval.firstPosition;
				end=pInterval.lastPosition;
			}
			else {
				if (pInterval.firstPosition>end) {
					sum+=end-start+1;
					start=pInterval.firstPosition;
					end=pInterval.lastPosition;
				}
				else end=Math.max(end,pInterval.lastPosition);
			}
		}
		if (start!=-1) sum+=end-start+1;
		out[2]+=sum;
		if (highQualitySurface-sum<=MAX_HOLE) out[3]++;
		
		// Surface covered by long-period periodic
		sum=0; start=-1; end=-1;
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			pInterval=PeriodicSubstrings.intervals[i];
			if (!pInterval.hasLongPeriod) continue;
			if (start==-1) {
				start=pInterval.firstPosition;
				end=pInterval.lastPosition;
			}
			else {
				if (pInterval.firstPosition>end) {
					sum+=end-start+1;
					start=pInterval.firstPosition;
					end=pInterval.lastPosition;
				}
				else end=Math.max(end,pInterval.lastPosition);
			}
		}
		if (start!=-1) sum+=end-start+1;
		out[4]+=sum;
		if (highQualitySurface-sum<=MAX_HOLE) out[5]++;
		
		// Surface covered by any interval
		nIntervals=PeriodicSubstrings.lastInterval+DenseSubstrings.lastSubstring+AlignmentIntervals.lastInterval+3;
		if (statWindows==null || statWindows.length<nIntervals) statWindows = new StatsWindow[nIntervals];
		lastWindow=-1; tag=-1;
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			lastWindow++;
			if (statWindows[lastWindow]==null) statWindows[lastWindow] = new StatsWindow(PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,++tag,0,0,0,Constants.INTERVAL_PERIODIC);
			else statWindows[lastWindow].set(PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,++tag,0,0,0,Constants.INTERVAL_PERIODIC);
		}
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			lastWindow++;
			if (statWindows[lastWindow]==null) statWindows[lastWindow] = new StatsWindow(DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,++tag,0,0,0,Constants.INTERVAL_DENSE_SUBSTRING);
			else statWindows[lastWindow].set(DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,++tag,0,0,0,Constants.INTERVAL_DENSE_SUBSTRING);
		}
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			lastWindow++;
			if (statWindows[lastWindow]==null) statWindows[lastWindow] = new StatsWindow(AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,++tag,0,0,0,Constants.INTERVAL_ALIGNMENT);
			else statWindows[lastWindow].set(AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,++tag,0,0,0,Constants.INTERVAL_ALIGNMENT);
		}
		if (lastWindow>0) Arrays.sort(statWindows,0,lastWindow+1);
		assemblyStatistics_impl(statWindows,lastWindow,tmpArray,null,null);
		out[6]+=tmpArray[0]; out[9]+=tmpArray[2]; out[10]+=tmpArray[3];
		if (highQualitySurface-tmpArray[0]<=MAX_HOLE) out[7]++;
		if ((highQualitySurface-tmpArray[0]<=MAX_HOLE || tmpArray[3]==1) && tmpArray[1]>0) {
			out[8]++;
			lastStatRead++;
			if (statReads==null) statReads = new int[GROWTH_RATE];
			else if (lastStatRead==statReads.length) {
				int[] newReads = new int[lastStatRead+GROWTH_RATE];
				System.arraycopy(statReads,0,newReads,0,statReads.length);
				statReads=newReads;
			}
			statReads[lastStatRead]=ReadA.id;
		}
		if (PeriodicSubstrings.lastInterval!=-1) {
			assemblyStatistics_impl_periodic(statWindows,lastWindow,tmpArray);
			out[13]+=tmpArray[0]; out[14]+=tmpArray[1];
		}
	}
	
	
	/**
	 * Analyzes $statWindows$, which is assumed to be sorted.
	 *
	 * @param out output array: 0=surface covered by any window; 1=number of positions 
	 * $i$ such that windows either end at $i$ or start at $i+1$, i.e. no window 
	 * contains both positions; 2=number of approximately identical pairs of windows with 
	 * different tag, and such that each window is approximately of the same length as its
	 * tag (we call such events "collisions"); 3=one iff the read starts or ends with a
	 * window;
	 * @param collisions output array: if not null, the procedure stores in row $i$ all 
	 * distinct tags $>i$ that have a collision with tag $i$ (not necessarily sorted);
	 * each tag is followed by the number of collisions found.
	 */
	public static final void assemblyStatistics_impl(StatsWindow[] statWindows, int lastWindow, int[] out, int[][] collisions, int[] lastCollision) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int IDENTITY_THRESHOLD_ENDS = IO.quantum>>1;
		final int DISTANCE_THRESHOLD = Alignments.minAlignmentLength;
		boolean found, startOrEnd;
		int i, j, k;
		int sum, start, end, nJunctions, nCollisions;
		int currentStart, currentEnd, currentTag, minTag, maxTag;
		
		sum=0; start=-1; end=-1; nJunctions=0; nCollisions=0; startOrEnd=false;
		for (i=0; i<=lastWindow; i++) {
			// Collisions
			if (collisions!=null && statWindows[i].isFullCopy(IDENTITY_THRESHOLD_ENDS)) {
				currentStart=statWindows[i].start; currentEnd=statWindows[i].end; 
				currentTag=statWindows[i].tag;
				for (j=i+1; j<=lastWindow; j++) {
					if (statWindows[j].start>currentStart+IDENTITY_THRESHOLD_ENDS) break;
					if (Math.abs(statWindows[j].end,currentEnd)>IDENTITY_THRESHOLD_ENDS) continue;
					if (statWindows[j].tag!=currentTag && statWindows[j].isFullCopy(IDENTITY_THRESHOLD_ENDS)) {
						nCollisions++;		
						minTag=Math.min(currentTag,statWindows[j].tag);
						maxTag=Math.max(currentTag,statWindows[j].tag);
						found=false;
						for (k=0; k<=lastCollision[minTag]; k+=2) {
							if (collisions[minTag][k]==maxTag) {
								collisions[minTag][k+1]++;
								found=true;
								break;
							}
						}
						if (found) continue;
						if (lastCollision[minTag]+2>=collisions[minTag].length) {
							int[] tmpArray = new int[collisions[minTag].length<<1];
							System.arraycopy(collisions[minTag],0,tmpArray,0,collisions[minTag].length);
							collisions[minTag]=tmpArray;
						}
						lastCollision[minTag]++;
						collisions[minTag][lastCollision[minTag]]=maxTag;
						lastCollision[minTag]++;
						collisions[minTag][lastCollision[minTag]]=1;
					}
				}
			}
			// Surface
			if (start==-1) {
				start=statWindows[i].start;
				end=statWindows[i].end;
				continue;
			}
			if (statWindows[i].start>end) {
				sum+=end-start+1;
				if (statWindows[i].start<=end+IDENTITY_THRESHOLD) nJunctions++;
				start=statWindows[i].start;
				end=statWindows[i].end;
			}
			else end=Math.max(end,statWindows[i].end);
			// Start/end
			if (statWindows[i].start<=IDENTITY_THRESHOLD || statWindows[i].end>=ReadA.readLength-IDENTITY_THRESHOLD) startOrEnd=true;
		}
		if (start!=-1) sum+=end-start+1;
		out[0]=sum; out[1]=nJunctions; out[2]=nCollisions; out[3]=startOrEnd?1:0;
	}
	
	
	/**
	 * Remark: $statWindows$ is altered by the procedure.
	 * 
	 * @param statWindows assumed to be sorted and containing at least one periodic 
	 * interval;
	 * @param out stores in cell 0 the number of times a periodic interval is flanked by a
	 * high-quality substring that is not covered by any repeat, and in cell 1 the number 
	 * of times we see configuration XYZ, where X,Z are parts of a read covered by
	 * periodic intervals, and Y is fully covered by intervals that are not periodic.
	 */
	public static final void assemblyStatistics_impl_periodic(StatsWindow[] statWindows, int lastWindow, int[] out) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int FLANKING_THRESHOLD = Alignments.minAlignmentLength>>2;  // Arbitrary
		final int MIN_EMBEDDING_SIZE = Alignments.minAlignmentLength;
		final int EMBEDDING_THRESHOLD = IO.quantum<<1;  // Arbitrary
		boolean found, foundLeft, foundRight;
		int i, j, k;
		int window, start, end, nFlanks, nEmbeddings;
		
		// Merging all periodic windows
		window=-1; end=-1;
		for (i=0; i<=lastWindow; i++) {
			if (statWindows[i].type!=Constants.INTERVAL_PERIODIC) continue;
			if (end==-1) {
				window=i; end=statWindows[i].end;
				continue;
			}
			if (statWindows[i].start<=end+IDENTITY_THRESHOLD) {
				end=Math.max(end,statWindows[i].end);
				statWindows[i].type=-1;
			}
			else {
				statWindows[window].end=end;
				window=i; end=statWindows[i].end;
			}
		}
		statWindows[window].end=end;		
		
		// Flanking non-repetitive regions
		nFlanks=0;
		// First periodic window
		for (i=0; i<=lastWindow; i++) {
			if (statWindows[i].type!=Constants.INTERVAL_PERIODIC) continue;
			if (statWindows[i].start<FLANKING_THRESHOLD) break;
			foundRight=false;
			for (j=i-1; j>=0; j--) {
				if (statWindows[j].start<statWindows[i].start-IDENTITY_THRESHOLD && statWindows[j].end>=statWindows[i].start-FLANKING_THRESHOLD) {
					foundRight=true;
					break;
				}
			}
			if (!foundRight && !Reads.hasLowQuality(ReadA.id,statWindows[i].start-FLANKING_THRESHOLD,statWindows[i].start-1,true)) nFlanks++;
			break;
		}
		// Other periodic windows
		while (i<=lastWindow) {
			if (statWindows[i].type!=Constants.INTERVAL_PERIODIC) {
				i++;
				continue;
			}
			found=false;
			for (j=i+1; j<=lastWindow; j++) {
				if (statWindows[j].type!=Constants.INTERVAL_PERIODIC) continue;
				found=true;
				if (statWindows[j].start<statWindows[i].end+FLANKING_THRESHOLD) {
					i=j; break;
				}
				foundLeft=false; foundRight=false;
				for (k=i+1; k<j; k++) {
					if (statWindows[k].type==-1) continue;
					if ( (statWindows[k].start<statWindows[i].end-IDENTITY_THRESHOLD && statWindows[k].end>statWindows[i].end+IDENTITY_THRESHOLD) || 
						 (statWindows[k].start>=statWindows[i].end-IDENTITY_THRESHOLD && statWindows[k].start<=statWindows[i].end+FLANKING_THRESHOLD)
					   ) foundLeft=true;
					if (statWindows[k].start<statWindows[j].start-IDENTITY_THRESHOLD && statWindows[k].end>statWindows[j].start-FLANKING_THRESHOLD) foundRight=true;
				}
				if (!foundLeft && !Reads.hasLowQuality(ReadA.id,statWindows[i].end+1,statWindows[i].end+FLANKING_THRESHOLD,true)) nFlanks++;
				if (!foundRight && !Reads.hasLowQuality(ReadA.id,statWindows[j].start-FLANKING_THRESHOLD,statWindows[j].start-1,true)) nFlanks++;
				i=j; break;
			}
			if (!found) {  // Last periodic window
				if (statWindows[i].end<Reads.getReadLength(ReadA.id)-FLANKING_THRESHOLD) {
					foundLeft=false;
					for (j=i+1; j<=lastWindow; j++) {
						if (statWindows[j].type==-1) continue;
						if ( (statWindows[j].start<statWindows[i].end-IDENTITY_THRESHOLD && statWindows[j].end>statWindows[i].end+IDENTITY_THRESHOLD) || 
						     (statWindows[j].start>=statWindows[i].end-IDENTITY_THRESHOLD && statWindows[j].start<=statWindows[i].end+FLANKING_THRESHOLD)
						   ) {
							foundLeft=true;
							break;
						}
					}
					if (!foundLeft && !Reads.hasLowQuality(ReadA.id,statWindows[i].end+1,statWindows[i].end+FLANKING_THRESHOLD,true)) nFlanks++;
				}
				break;
			}
		}
		
		// Repeats surrounded by periodic regions
		nEmbeddings=0; i=0;
		while (i<=lastWindow) {
			if (statWindows[i].type!=Constants.INTERVAL_PERIODIC) {
				i++;
				continue;
			}
			found=false;
			for (j=i+1; j<=lastWindow; j++) {
				if (statWindows[j].type!=Constants.INTERVAL_PERIODIC) continue;
				found=true;
				if (statWindows[j].start<statWindows[i].end+MIN_EMBEDDING_SIZE-IDENTITY_THRESHOLD) {
					i=j;
					break;
				}
				start=-1; end=-1;
				for (k=i+1; k<j; k++) {
					if (statWindows[k].type==-1) continue;
					if (end!=-1 && statWindows[k].start>end+IDENTITY_THRESHOLD) {
						if (statWindows[k].start>statWindows[i].end+IDENTITY_THRESHOLD) {
							start=Math.POSITIVE_INFINITY; break;
						}
						else {
							start=statWindows[k].start; end=statWindows[k].end;
							continue;
						}
					}
					if (start==-1) {
						start=statWindows[k].start; end=statWindows[k].end;
					}
					else end=Math.max(end,statWindows[k].end);
				}
				if (start<=statWindows[i].end+EMBEDDING_THRESHOLD && end>=statWindows[j].start-EMBEDDING_THRESHOLD) nEmbeddings++;
				i=j;
				break;
			}
			if (!found) break;  // Last periodic window
		}
		
		out[0]=nFlanks; out[1]=nEmbeddings;
	}
	
	
	public static class StatsWindow implements Comparable {
		public int start, end, tag, tagLength, type;
		public int otherStart, otherEnd;
		
		public StatsWindow() { 
			start=-1; end=-1; tag=-1; type=-1;
			otherStart=-1; otherEnd=-1;
		}
		
		public StatsWindow(int s, int e, int t, int tl, int os, int oe, int tp) {
			start=s; end=e; tag=t; tagLength=tl; type=tp;
			otherStart=os; otherEnd=oe;
		}
		
		public void set(int s, int e, int t, int tl, int os, int oe, int tp) {
			start=s; end=e; tag=t; tagLength=tl; type=tp;
			otherStart=os; otherEnd=oe;
		}
		
		public int compareTo(Object other) {
			StatsWindow otherWindow = (StatsWindow)other;
			if (start<otherWindow.start) return -1;
			else if (start>otherWindow.start) return 1;
			return 0;
		}
		
		public boolean isFullCopy(int threshold) {
			return Math.abs(otherStart)<=threshold && Math.abs(otherEnd,tagLength)<=threshold;
		}
	}
	
	
	/**
	 * @param out cell 11 is assumed to contain the total number of reads that have been 
	 * seen (but not necessarily factorized); cell 12 is assumed to contain the sum of 
	 * lengths of all reads that have been factorized.
	 */
	private static final void printAssemblyStatistics(long[] out, boolean humanReadable) {
		System.out.println("Assembly statistics: ");
		if (!humanReadable) {
			for (int i=0; i<=14; i++) System.out.print(out[i]+",");
			System.out.println();
		}
		else {
			System.out.println("High-quality read surface: "+IO.getPercent(out[1],out[12])+"%");
			System.out.println("Reads with at least one interval: "+IO.getPercent(out[0],out[11])+"% of all reads");
			System.out.println("Total short-period surface: "+IO.getPercent(out[2],out[1])+"% of high-quality");
			System.out.println("Reads fully short-period: "+IO.getPercent(out[3],out[0])+"% of reads with an interval");
			System.out.println("Total long-period surface: "+IO.getPercent(out[4],out[1])+"% of high-quality");
			System.out.println("Reads fully long-period: "+IO.getPercent(out[5],out[0])+"% of reads with an interval");
			System.out.println("Total surface in intervals: "+IO.getPercent(out[6],out[1])+"% of high-quality");
			System.out.println("Reads fully covered by intervals: "+IO.getPercent(out[7],out[0])+"% of reads with an interval");
			System.out.println("Reads fully covered by intervals, likely with distinct repeats: "+IO.getPercent(out[8],out[0])+"% of reads with an interval");
			System.out.println("Reads with start or end in an interval: "+IO.getPercent(out[10],out[0])+"% of reads with an interval");
			System.out.println("Number of colliding pairs of intervals: "+out[9]);
			System.out.println("Unique/periodic boundaries: "+out[13]);
			System.out.println("Insertions of repeats into periodic: "+out[14]);
		}
		System.out.println("Reads fully covered by intervals, or with intervals at start/end, and likely with distinct repeats: ");
		for (int i=0; i<=lastStatRead; i++) System.out.print(statReads[i]+",");
		System.out.println();
	}
	
	
	private static final void deallocateStats() {
		for (int i=0; i<statWindows.length; i++) statWindows[i]=null;
		statWindows=null; statReads=null;
	}
	
	
	
	
	
	






	private static final void probe(int id) {
/*		final int FROM1 = 7209;
		final int TO1 = 9863;
//		final int FROM2 = 0;
//		final int TO2 = 0;
//		final int FROM_B = 4450;
//		final int TO_B = 7760;
		final int THRESHOLD = IO.quantum;
		int startA, endA, startB, endB;
		
		for (int i=0; i<=ReadA.lastSortedAlignment; i++) {
			startA=ReadA.sortedAlignments[i].startA();
			endA=ReadA.sortedAlignments[i].endA();
			//startB=ReadA.sortedAlignments[i].startB();
			//endB=ReadA.sortedAlignments[i].endB();
			if (startA==FROM1 && endA==TO1  Math.abs(startA,FROM1)<=THRESHOLD && Math.abs(endA,TO1)<=THRESHOLD) {
				System.err.println("PROBE "+id+"> alignment found: "+ReadA.sortedAlignments[i].toStringBoundaries()+" => "+ReadA.sortedAlignments[i].toStringPointers());
				System.err.println("impliedByPrefixSubstringPrime="+ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime);
				System.err.println("impliedBySuffixSubstringPrime="+ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime);
			}
			
			
			
			// if (Math.abs(endA,TO)<=THRESHOLD) {
			// 	if (ReadA.sortedAlignments[i].periodicSubstringInterval!=null && ReadA.sortedAlignments[i].periodicSubstringInterval.lastPosition>TO+500) {
			// 		System.err.println("PROBE "+id+"> alignment assigned to wrong short-period interval? "+ReadA.sortedAlignments[i].toStringBoundaries()+" => "+ReadA.sortedAlignments[i].toStringPointers());
			// 	else System.err.println("PROBE "+id+"> selected alignment correct? "+ReadA.sortedAlignments[i].toStringBoundaries()+" => "+ReadA.sortedAlignments[i].toStringPointers());
			// }
		}
*/
	}


}