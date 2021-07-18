package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import java.io.IOException;
import java.io.BufferedWriter;

import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Leaf;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.DAG;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.RegressionTree;
import de.mpi_cbg.revant.util.Histograms;
import de.mpi_cbg.revant.util.Leaves;


public class DenseSubstrings {
	/**
	 * Parameters of the pipeline
	 */
	private static int MAX_NEIGHBORS;
	private static final int identityThreshold = IO.quantum;  // Two positions are considered identical if they differ by at most this threshold
	private static final int maxDistance = IO.quantum*3;  // Only pairs of alignments whose starting positions in $readA$ are at most $maxDistance$ apart are considered when searching for dense substrings of substring type.
	private static final int repeatThreshold = identityThreshold<<1;  // A dense substring of substring type is considered a simple repeat (and thus discarded) if the first (or, respectively, the last) position of its last alignment is within this distance from the first (respectively, last) position of its first alignment.
	private static final int maxDistanceEnd = Alignments.minAlignmentLength-identityThreshold;
	public static final int maxDistanceEndWeak = Alignments.minAlignmentLength>>1;
	private static final int maxDistanceEndWeakPrime = Alignments.minAlignmentLength;
	private static final int distanceThreshold = IO.quantum<<1;  // Used by procedure $splitSubstrings$
	public static final int minDeletionLength = IO.quantum;  // Minimum length for a deletion to be detected
	private static final int eventThreshold = IO.quantum<<1;
	private static final double mergeJaccardThreshold = 0.7;  // Used instead of $Intervals.jaccardThreshold$ for detecting identical dense substrings
	private static final double mergePrefixSuffixThreshold = 0.75;
	private static final int prefixSuffixThreshold = IO.quantum*3;  // A dense substring of prefix (respectively, suffix) type must have at least this fraction of its interval in readA covered by end positions of prefixes (start positions of suffixes).
	public static final int DEFAULT_MIN_PATH_LENGTH = 5;  // Arbitrary
	public static final int DEFAULT_MAX_PATH_LENGTH = 20;
	private static final int MAX_PATHLENGTH_DIFFERENCE = 4;  // If there are "few" path length values, and if such values differ by at most this quantity, we assume they are uniformly distributed and we don't try to split them.
	private static int MIN_DELETION_OCCURRENCES;  // Minimum number of occurrences to decide a single deletion
	private static final int MAX_DELETION_DIFFERENCE = IO.quantum;  // If there are "few" observed single-deletion lengths, this is the maximum difference between such lengths to decide that they are uniformly distributed.
	private static final double MIN_LENGTH_BOUNDARY_DENSITY_RATIO = 2.0;  // Only length bias candidates whose density is at least this number times the density of their corresponding non-candidate region are considered as length biases.
	private static final int MIN_SHIFT_DIFFERENCE = IO.quantum;  // Minimum value of a shift in prefix/suffix substrings to be detected
	private static final int MIN_ALIGNMENTS_FOR_LENGTH_BIAS = 30;  // Minimum number of alignments, implied by a prefix/suffix substring, to estimate biases in their prefix/suffix lengths.
	public static final int LONG_PERIOD_THRESHOLD = IO.coverage*10;  // Arbitrary

	/**
	 * Parameters of the pipeline, estimated from the data.
	 */
	private static int minPathLength;  // Minimum length of a path to represent a dense substring

	/**
	 * Output data structures
	 */
	public static DenseSubstring[] substrings;
	public static int lastSubstring;
	public static Interval[] lengthBiasIntervals;
	public static int lastLengthBiasInterval;

	/**
	 * DAG data structures
	 */
	private static int[][] inNeighbors, outNeighbors;  // DAG of alignments
	private static int[] nInNeighbors, nInNeighborsPrime, nOutNeighbors;
	private static int[] minimalVertices;  // Minimal vertices of the DAG
	private static int[] components;  // Connected component of each vertex of the DAG
	private static int[] componentLongestPath;  // Length of a longest path in each connected component of the DAG
	private static int[] nVerticesPerComponent;  // Number of minimal vertices per connected component of the DAG
	private static int[] sorted2original;  // Maps a rank in the topological order of the DAG to the ID of the corresponding vertex
	private static int[] original2sorted;  // Reverses $sorted2original$
	private static int[] distances;  // Length of a best path from a source alignment to all other alignments
	private static int[] predecessors;  // Predecessor of each vertex in the best path from a source alignment
	private static int[] sumOfLengths;  // Sum of the lengths of all alignments in the best path from a source alignment
	private static int[] tmp;

	/**
	 * Data structures for computing splits
	 */
	public static int[] leftSplits, rightSplits, tmpSplits;
	public static double[] deltaLeft, deltaRight;
	
	/**
	 * Data structures for estimating periods
	 */
	private static Point[] shifts;

	/**
	 * Data structures for connected components
	 */
	public static int[] stack, stackPrime;

	/**
	 * Temporary space
	 */
	private static Point[] singleDeletions, lengthPoints, weakLeft, weakRight;
	private static double[] deltaPrefix, deltaSuffix;
	private static int[] countsPrefix, countsSuffix;
	private static Point[] deltaPointsPrefix, deltaPointsSuffix;
	private static int lastSingleDeletion;
	private static int[] singleDeletionTmp;  // Cell 0: estimated number of distinct deletion lengths (-1 if no deletion; 0 if lengths could not be estimated). Cell 1: shortest length. Cell 2: longest length.
	private static int[] lengthBoundariesPrefix, lengthBoundariesSuffix;  // First, center of mass, and last value, of all intervals of alignment lengths, in prefix/suffix substrings, that have a peak number of alignments.
	private static boolean[] lengthBoundariesFlags;
	private static int lastLengthBoundaryPrefix, lastLengthBoundarySuffix;
	private static int lastWeakLeft, lastWeakRight;
	private static Leaf[] tmpLeaves;
	private static Interval tmpInterval;
	private static PeriodicSubstringInterval tmpPInterval;
	private static Alignment[] tmpAlignments;
	private static int lastTmpAlignment;
	private static int[] weakTmp, shiftTmp;
	private static final int GROWTH_RATE = 10;
	private static Event[] deltaPoints_backupLeft, deltaPoints_backupRight;
	private static int lastDeltaPoint_backupLeft, lastDeltaPoint_backupRight;
	private static DenseSubstring tmpDenseSubstring;


	public static final void allocateMemory(int maxAlignments, int maxSubstrings, int maxSplits, int maxNeighbors, int maxReadLength) {
		int i;		
		
		MIN_DELETION_OCCURRENCES=(IO.minOccurrencesInGenome-1)*IO.coverage;
		MAX_NEIGHBORS=maxNeighbors;

		// Output data structures
		substrings = new DenseSubstring[maxSubstrings];  // This is not $maxFactors$ since, in the current implementation, we store a large number of substrings in this array before merging them.
		for (i=0; i<maxSubstrings; i++) substrings[i] = new DenseSubstring();
		lengthBiasIntervals = new Interval[maxAlignments];
		for (i=0; i<maxAlignments; i++) lengthBiasIntervals[i] = new Interval();

		// DAG data structures
		inNeighbors = new int[maxAlignments][MAX_NEIGHBORS];
		nInNeighbors = new int[maxAlignments];
		nInNeighborsPrime = new int[maxAlignments];
		outNeighbors = new int[maxAlignments][MAX_NEIGHBORS];
		nOutNeighbors = new int[maxAlignments];
		minimalVertices = new int[maxAlignments];
		components = new int[maxAlignments];
		componentLongestPath = new int[maxAlignments];
		nVerticesPerComponent = new int[maxAlignments];
		sorted2original = new int[maxAlignments];
		original2sorted = new int[maxAlignments];
		distances = new int[maxAlignments];
		predecessors = new int[maxAlignments];
		sumOfLengths = new int[maxAlignments];
		tmp = new int[7];

		// Data structures for computing splits
		leftSplits = new int[maxSplits*3];
		rightSplits = new int[maxSplits*3];
		tmpSplits = new int[maxSplits*3];
		deltaLeft = new double[maxReadLength];
		deltaRight = new double[maxReadLength];
		
		// Data structures for estimating periods
		shifts = new Point[maxAlignments<<2];

		// Data structures for connected components
		stack = new int[maxSubstrings];
		stackPrime = new int[maxAlignments];

		// Temporary space
		singleDeletions = new Point[maxAlignments<<1];
		for (i=0; i<singleDeletions.length; i++) singleDeletions[i] = new Point();
		singleDeletionTmp = new int[3];
		lengthPoints = new Point[maxAlignments>>1];
		for (i=0; i<lengthPoints.length; i++) lengthPoints[i] = new Point();
		deltaPrefix = new double[maxReadLength];
		deltaSuffix = new double[maxReadLength];
		countsPrefix = new int[maxReadLength];
		countsSuffix = new int[maxReadLength];
		lengthBoundariesPrefix = new int[maxAlignments];
		lengthBoundariesSuffix = new int[maxAlignments];
		lengthBoundariesFlags = new boolean[maxAlignments];
		tmpLeaves = new Leaf[DensityEstimationTree.leaves.length];
		for (i=0; i<tmpLeaves.length; i++) tmpLeaves[i] = new Leaf();
		tmpInterval = new Interval();
		tmpPInterval = new PeriodicSubstringInterval();
		tmpAlignments = new Alignment[maxAlignments];
		for (i=0; i<tmpAlignments.length; i++) tmpAlignments[i] = new Alignment();
		deltaPointsPrefix = new Point[maxAlignments<<1];
		for (i=0; i<deltaPointsPrefix.length; i++) deltaPointsPrefix[i] = new Point();
		deltaPointsSuffix = new Point[maxAlignments<<1];
		for (i=0; i<deltaPointsSuffix.length; i++) deltaPointsSuffix[i] = new Point();
		weakLeft = new Point[maxAlignments];
		for (i=0; i<weakLeft.length; i++) weakLeft[i] = new Point();
		weakRight = new Point[maxAlignments];
		for (i=0; i<weakRight.length; i++) weakRight[i] = new Point();
		weakTmp = new int[3];
		shiftTmp = new int[2];
		deltaPoints_backupLeft = new Event[ReadA.deltaPoints.length];
		for (i=0; i<deltaPoints_backupLeft.length; i++) deltaPoints_backupLeft[i] = new Event();
		deltaPoints_backupRight = new Event[ReadA.deltaPoints.length];
		for (i=0; i<deltaPoints_backupRight.length; i++) deltaPoints_backupRight[i] = new Event();
		tmpDenseSubstring = new DenseSubstring();
	}


	private static final void clearNeighborMatrices(int nAlignments) {
		int i, j;

		for (i=0; i<nAlignments; i++) {
			Math.set(inNeighbors[i],inNeighbors[i].length-1,-1);
		}
		Math.set(nInNeighbors,nAlignments-1,0);
		for (i=0; i<nAlignments; i++) {
			Math.set(outNeighbors[i],outNeighbors[i].length-1,-1);
		}
		Math.set(nOutNeighbors,nAlignments-1,0);
	}


	public static final void detect() {
		boolean mergedToDense;
		int i, nAlignmentsLow, nAlignmentsHigh;
		DenseSubstring denseSubstring;

		nAlignmentsHigh=getSubstrings();

		
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER GETSUBSTRINGS ONLY, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"], leftMax="+ReadA.sortedAlignments[i].isLeftMaximal+"("+ReadA.sortedAlignments[i].isLeftMaximalB+") rightMax="+ReadA.sortedAlignments[i].isRightMaximal+"("+ReadA.sortedAlignments[i].isRightMaximal+") || "+(ReadA.sortedAlignments[i].impliedByDenseSubstring!=null?ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode():"null")+","+(ReadA.sortedAlignments[i].inDenseSubstring!=null?ReadA.sortedAlignments[i].inDenseSubstring.hashCode():"null")+","+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null?ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode():"null")+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()) );
	}
}	
		
		if (lastSubstring>=0) {
			nAlignmentsLow=mergeSubstrings(nAlignmentsHigh-1);
			
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER GETSUBSTRINGS DENSE AND MERGE, read "+ReadA.id+", nAlignmentsLow="+nAlignmentsLow+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+ (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null?ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode():"null")+","+(ReadA.sortedAlignments[i].inDenseSubstring!=null?ReadA.sortedAlignments[i].inDenseSubstring.hashCode():"null")+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()) );
	}
}	

			mergedToDense=discardPeriodicSubstrings(nAlignmentsHigh);
			
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER discardPeriodicSubstrings, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null?"null":ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}			
			
			discardPeriodicSubstrings_prefixSuffix(nAlignmentsHigh);
			
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER discardPeriodicSubstrings_prefixSuffix, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null?"null":ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}				
			
			// Sorting again, because of $discardPeriodicSubstrings$.
			Alignment.order=Alignment.IMPLIEDBYDENSE_INDENSE_PERIODIC_STARTA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment);
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) {
					nAlignmentsLow=i;
					break;
				}
			}
			for (i=nAlignmentsLow; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) {
					nAlignmentsHigh=i;
					break;
				}
			}
			
			discardPrefixSuffix_weakSubstring(nAlignmentsHigh);
			
			// Sorting again, because of $discardPrefixSuffix_weakSubstring$.
			Alignment.order=Alignment.IMPLIEDBYDENSE_INDENSE_PERIODIC_STARTA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment);
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) {
					nAlignmentsLow=i;
					break;
				}
			}
			for (i=nAlignmentsLow; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) {
					nAlignmentsHigh=i;
					break;
				}
			}			
			if (nAlignmentsHigh-nAlignmentsLow>0) {
if (IO.SHOW_INTERACTIVE) IO.printErr("START SPLITTING DENSE SUBSTRINGS OF read "+ReadA.id+":");
				Alignment.order=Alignment.STARTA_ENDA;
				if (nAlignmentsHigh>nAlignmentsLow+1) Arrays.sort(ReadA.sortedAlignments,nAlignmentsLow,nAlignmentsHigh);
				splitSubstrings(nAlignmentsLow,nAlignmentsHigh-1);
				
				// Sorting again, since $splitSubstrings$ might have changed the 
				// start/end positions of dense substrings.
				DenseSubstring.order=DenseSubstring.STARTA;
				if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
				
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER GETSUBSTRINGS DENSE AND MERGE AND SPLIT, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh+" nAlignmentsLow="+nAlignmentsLow);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}
					
				// Merging again, since $splitSubstrings$ might have changed the
				// start/end positions of dense substrings.
				Alignment.order=Alignment.STARTA_ENDA;
				if (nAlignmentsHigh>1) Arrays.sort(ReadA.sortedAlignments,0,nAlignmentsHigh);
				nAlignmentsLow=mergeSubstrings(nAlignmentsHigh-1);
				
				
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER MERGESUBSTRINGS, read "+ReadA.id+":");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh+" nAlignmentsLow="+nAlignmentsLow);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}				
				
				if (nAlignmentsHigh-nAlignmentsLow>0) {
					Alignment.order=Alignment.STARTA_ENDA;
					if (nAlignmentsHigh>nAlignmentsLow+1) Arrays.sort(ReadA.sortedAlignments,nAlignmentsLow,nAlignmentsHigh);
					estimateSubstringPeriods(nAlignmentsLow,nAlignmentsHigh-1);
				}
			}

if (IO.SHOW_INTERACTIVE) {
	IO.printErr("BEFORE discardPeriodicSubstrings, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh+" nAlignmentsLow="+nAlignmentsLow);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null?"null":ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}

			
			getSingleDeletions(nAlignmentsLow-1);
			
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER getSingleDeletions, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null?"null":ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}				
			
			Alignment.order=Alignment.STARTA_ENDA;
			if (nAlignmentsLow>1) Arrays.sort(ReadA.sortedAlignments,0,nAlignmentsLow);
			addEvents(nAlignmentsLow-1);
			
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER addEvents, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null?"null":ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}				
			
			markStartEnd(nAlignmentsHigh-1);
			
			
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER markStartEnd, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null?"null":ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}
			
			setMaximality();
			
			
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("AFTER setMaximality, read "+ReadA.id+":");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments: nAlignmentsHigh="+nAlignmentsHigh);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"]"+(ReadA.sortedAlignments[i].isLeftMaximal)+(ReadA.sortedAlignments[i].isRightMaximal)+","+(ReadA.sortedAlignments[i].isLeftMaximalB)+(ReadA.sortedAlignments[i].isRightMaximalB)+","+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"null":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null?"null":ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode())+"|"+(ReadA.sortedAlignments[i].periodicSubstringInterval==null?"null":ReadA.sortedAlignments[i].periodicSubstringInterval.hashCode()));
	}
}

			if (mergedToDense) {
				Alignment.order=Alignment.UNSORTED;
				PeriodicSubstrings.splitLongPeriodIntervals();
			}
			
			// At this point it could still happen that an alignment assigned to a dense
			// substring is not compatible with it. We disconnect all such alignments, and
			// update substring boundaries if necessary.
			for (i=0; i<=lastSubstring; i++) {
				substrings[i].previousSumStartA=Math.POSITIVE_INFINITY;
				substrings[i].previousSumEndA=0;
			}
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				denseSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
				if (denseSubstring==null) continue;
				if (!denseSubstring.canBeAssigned(ReadA.sortedAlignments[i],identityThreshold)) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					continue;
				}
				if (denseSubstring.prefixReplication || denseSubstring.singleDeletionReplication) denseSubstring.previousSumEndA=Math.max(denseSubstring.previousSumEndA,ReadA.sortedAlignments[i].endA());
				if (denseSubstring.suffixReplication || denseSubstring.singleDeletionReplication) denseSubstring.previousSumStartA=Math.min(denseSubstring.previousSumStartA,ReadA.sortedAlignments[i].startA());
			}
			for (i=0; i<=lastSubstring; i++) {
				if (substrings[i].previousSumStartA!=Math.POSITIVE_INFINITY && substrings[i].previousSumStartA>substrings[i].startA+identityThreshold) substrings[i].startA=substrings[i].previousSumStartA;
				if (substrings[i].previousSumEndA!=0 && substrings[i].previousSumEndA<substrings[i].endA-identityThreshold) substrings[i].endA=substrings[i].previousSumEndA;				
			}			
		}
	}


	/**
	 * Let $R[i..j]$ be an occurrence of the whole maximal repeat $A$ in read $R$. If $A$
	 * replicates in the genome by copying arbitrary substrings of itself, we should see
	 * a set of alignments $[i'..j']$ with $i \leq i' \leq j' \leq j$ and no bias in the
	 * position of $i'$ and $j'$. If $A$ replicates by taking arbitrary prefixes up to a
	 * minimum length, we should see a set of alignments $[i'..j']$ with $i' \approx i$
	 * and $j' \in [i..j]$, with no bias in the position of $j'$ (symmetrically for
	 * suffixes). Substring $[i..j]$ thus becomes \emph{dense} with maximal events.
	 * The total number of such alignments is approximately proportional to the
	 * frequency of $A$ in the read set. Some of the alignments in a dense region might
	 * not be maximal, e.g. because the dense region is adjacent to two long random
	 * insertions in $R$. Thus, in practice, we require maximality only on readB.
	 *
	 * The procedure tries to find all maximally-long dense regions, using a greedy
	 * strategy based on computing longest paths in a DAG whose vertices are all the
	 * alignments of $R$, with any $R'$, and whose arcs connect a source alignment to a
	 * destination alignment iff they overlap in $R$, and if the distance between the
	 * starting position in $R$ of the destination alignment and the starting position in 
	 * $R$ of the source alignment is at most $maxDistance$.
	 *
	 * Remark: longest paths are just a greedy heuristic, and they have a number of
	 * issues. For example, they might merge a prefix substring that is nested inside a
	 * suffix substring, with the rightmost piece of the suffix substring, making their
	 * right boundaries coincide even though they do not.
	 *
	 * Remark: the same repeat could use multiple modes of replication. The procedure
	 * tries to detect all such modes. However, a dense substring that replicates by
	 * taking substrings of itself cannot have other modes of replication.
	 *
	 * Remark: assume that a maximal repeat that occurs entirely in read $R$ occurs also
	 * multiple times in other reads, but always broken by long random insertions. Assume
	 * also that the alignment algorithm outputs maximal alignments. This does not
	 * create a substring that is dense of maximal events in $R$, since maximal events
	 * occur at most around the beginning and the end of the repeat in $R$.
	 *
	 * Remark: a longest path found by this procedure could cross a number of adjacent,
	 * distinct maximal repeats, if such repeats also replicate together in the genome,
	 * i.e. if their concatenation replicates by taking prefixes, suffixes or substrings.
	 * It is left to the following stages of the pipeline to separate such repeats.
	 *
	 * Remark: the procedure assumes that periodic substrings have already been detected.
	 * In theory, a dense substring could be contained in a short-period substring: this 
	 * captures the case of, e.g., a period of a periodic substring that diverges and 
	 * replicates autonomously, by taking, say, prefixes of itself. However, in practice, 
	 * short-period intervals contain many alignments that are not implied by a periodic 
	 * substring, e.g. because the aligner does not return all alignments, or because
	 * periodic substring paths are not long enough: such alignments are likely to trace 
	 * dense substrings just by chance. Thus, the procedure does not allow dense 
	 * substrings to use alignments implied by a periodic substring of any type, or 
	 * alignments contained in a maximal range of straddling short-period intervals 
	 * (alignment intervals inside short-period intervals will instead be allowed, if they
	 * align to peaks). However, a dense substring can contain a periodic substring of any
	 * type, or even straddle one. This captures the case of, e.g., a repeat that 
	 * replicates by taking prefixes of itself, and that contains a periodic substring. 
	 * Even a dense substring that replicates by taking substrings of itself can contain a
	 * periodic substring.
	 *
	 * Remark: the procedure does not assume anything about how dense substrings intersect
	 * with other dense substrings. E.g. dense substrings are allowed to nest into one 
	 * another, and to straddle each other. A dense substring that replicates by taking 
	 * substrings of itself cannot contain other dense substrings, but it can straddle 
	 * them, and it can also straddle other dense substrings of substring type.
	 *
	 * Remark: the DAGs built by this procedure typically have a number of connected
	 * components, some of which are induced by noise and expected to be small. The
	 * procedure estimates the distribution of the length of the longest path of all
	 * connected components of size bigger than one, and it sets a threshold that removes
	 * the leftmost peak, which is assumed to be induced by noise. Such threshold (the
	 * minimum length of a longest path of a connected component that is not induced by
	 * noise) is used as an estimate of the minimum length of a longest path (in a
	 * connected component that is not induced by noise) that should be considered a dense
	 * substring.
	 *
	 * Remark: it is important to use the length of longest paths rather than the size of
	 * a connected component. Consider e.g. a connected component of prefix type that
	 * consists of $k$ prefixes of increasing length, each appearing $h$ times in the read
	 * set. The size of the connected component, $hk$, can be significantly larger than
	 * $k$.
	 *
	 * Remark: removing the first peak as described above can remove true dense
	 * substrings. E.g. assume that there are multiple suffix substrings whose suffix
	 * regions are of different lengths, but that all have the same number of suffixes in
	 * the genome. Since there is a lower bound on the distance between the starting
	 * positions of two consecutive suffixes, the shorter suffix region corresponds to a
	 * connected component with shorter longest path, and it could thus be removed. This
	 * does not change if e.g. we divide the number of suffixes by the length of the
	 * suffix region.
	 *
	 * Remark: the procedure stores all biases in prefix/suffix lengths in array
	 * $lengthBiasIntervals$, sorted by starting position but not necessarily by ending
	 * position.
	 *
	 * Remark: after detecting dense substrings of substring, prefix, and suffix type, the
	 * procedure detects also "weak" prefix/suffix substrings. A weak prefix substring
	 * consists of a sequence of alignments with consecutive end positions, which are not 
	 * implied by any of the previous substring types, and which are B-right-maximal but 
	 * not necessarily B-left-maximal. Such a sequence is classified as a prefix substring 
	 * iff all the B-left-maximal alignments that belong to the sequence have similar 
	 * starting positions, and no alignment in the sequence that is not B-left-maximal
	 * starts before them. Otherwise, the sequence is classified as a dense substring of
	 * substring type.
	 *
	 * Remark: weak dense substrings are a way to exploit some of the alignments that
	 * are not used by the normal type of dense substrings, e.g. to trace a suffix 
	 * substring even though its B-left- and B-right-maximal alignments are too far from
	 * each other, or too few, to be detected. This also helps in detecting dense 
	 * substrings of the genome that are longer than the max. length of a read: e.g. a 
	 * read that spans the left boundary of a dense substring of suffix type in the genome
	 * contains a weak dense substring of suffix type which consists of a number of 
	 * B-left-maximal alignments, only some of which are also B-right-maximal, and such 
	 * that all B-right-maximal alignments end at the right end of the read.	
	 * However, there is no easy way to tell whether e.g. a weak dense substring of suffix
	 * type comes from a long dense substring of suffix type in the genome or not.
	 * 
	 * Remark: weak dense substrings of substring type are not treated like normal dense 
	 * substrings in the rest of the pipeline. See $cleanDenseSubstringsOfSubstringType$.
	 *
	 * Remark: weak dense substrings could pick up non-weak dense substrings that were 
	 * discarded because e.g. not long enough. This is not very likely, since too short 
	 * weak dense substrings are themselves discarded.
	 *
	 * @return the number of alignments in $ReadA.sortedAlignments$ that are not implied
	 * by a periodic substring (yes, "periodic" here is correct, although it might seem
	 * that it should be replaced by "dense").
	 */
	 private static final int getSubstrings() {
		final int MIN_LENGTH_DIFFERENCE = identityThreshold<<1;
	 	final double MIN_INTERVAL_LENGTH = 3.0;
		final double MIN_LOCAL_MAX_DISTANCE = 5.0;
		final int MIN_TRIMMED_LENGTH = Alignments.minAlignmentLength;  // Arbitrary
		final int IDENTITY_THRESHOLD_AMBIGUOUS = identityThreshold<<1;
		boolean closeToPrefix, closeToSuffix;
		int i, j;
		int startA1, startA2, endA1, endA2;
		int destination, pathLength, nAlignments, maximalStartA, maximalEndA;
		int nComponents, longestPathInComponent, nMarkedAlignments;
		int minBLeftMaximal, maxBLeftMaximal, minBRightMaximal, maxBRightMaximal;
		int leftAlignment, rightAlignment;
		int minBLeft, maxBRight, lengthPrefix, lengthSuffix;
		int tmpOrder;
		DenseSubstring prefixSubstring, suffixSubstring;
		Point[] tmpPoints;

		// Ensuring the necessary order in $ReadA.sortedAlignments$
		if (Alignment.order!=Alignment.IMPL_PERIODICSUBSTRING_STARTA_ENDA) {
			Alignment.order=Alignment.IMPL_PERIODICSUBSTRING_STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		lastLengthBiasInterval=-1;
		
		// Building the substring DAG of all alignments of $readA$ with any $readB$.
		i=0;
		while (i<=ReadA.lastSortedAlignment) {
			if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) break;
			i++;
		}
		nAlignments=i;
		if (nAlignments==0) {
			lastSubstring=-1;
			return 0;
		}
		clearNeighborMatrices(nAlignments);
		for (j=0; j<nAlignments; j++) {
			if (ReadA.sortedAlignments[j].isLeftMaximalB!=1 || ReadA.sortedAlignments[j].isRightMaximalB!=1) continue;
			if (ReadA.sortedAlignments[j].lowQualityStart) continue;
			if (ReadA.sortedAlignments[j].inPeriodicSubstring) continue;
			startA2=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
			endA2=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
			for (i=j-1; i>=0; i--) {
				startA1=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				if (startA2>startA1+maxDistance) break;
				if (startA2<=startA1+identityThreshold) continue;
				endA1=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
				if ( startA2>=endA1-identityThreshold ||
					 endA2<=endA1+identityThreshold ||
				     endA2>endA1+maxDistanceEnd
				   ) continue;
				if (ReadA.sortedAlignments[i].isLeftMaximalB!=1 || ReadA.sortedAlignments[i].isRightMaximalB!=1) continue;
				if (ReadA.sortedAlignments[i].lowQualityEnd) continue;
				if (ReadA.sortedAlignments[i].inPeriodicSubstring) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(nAlignments,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stackPrime);
		if (nComponents<=0 || nComponents>nAlignments) {
			IO.printCriticalErr("Error while detecting dense substrings: wrong number of connected components. nComponents="+nComponents+" nAlignments="+nAlignments);
			System.exit(1);
		}
		System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nAlignments);
		i=DAG.topologicalSort(nAlignments,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
		if (i!=0) {
			IO.printCriticalErr("Error while detecting dense substrings: the prefix DAG of all alignments of read "+ReadA.id+" contains a cycle that involves node "+(i-1));
			System.exit(1);
		}
		if (nComponents>lengthPoints.length) {
			tmpPoints = new Point[nComponents];
			System.arraycopy(lengthPoints,0,tmpPoints,0,lengthPoints.length);
			for (i=lengthPoints.length; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
			lengthPoints=tmpPoints;
		}
		longestPathInComponent=DAG.longestPathInComponent(nAlignments,sorted2original,nInNeighbors,inNeighbors,components,lengthPoints,distances);
		if (longestPathInComponent>0) {
			for (i=0; i<nComponents; i++) componentLongestPath[i]=(int)lengthPoints[i].position;
			if (nComponents==1) minPathLength=DEFAULT_MIN_PATH_LENGTH;
			else {
				minPathLength=estimateFromPathLengthPoints(nComponents,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE);
				if (minPathLength<DEFAULT_MIN_PATH_LENGTH) minPathLength=DEFAULT_MIN_PATH_LENGTH;
				if (minPathLength>DEFAULT_MAX_PATH_LENGTH) minPathLength=DEFAULT_MAX_PATH_LENGTH;
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstrings> minPathLength="+minPathLength+" for substring type");
			}
			if (IO.SHOW_STD_ERR) IO.printErr("nComponents in substring DAG: "+nComponents+" longestPathInComponent="+longestPathInComponent);

			// Detecting all dense substrings that replicate by taking substrings of
			// themselves. We don't mark implied alignments, since every alignment with
			// this type of dense substring would be marked.
			lastSubstring=-1;
			i=0;
			while (i<nAlignments) {
				if ( ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
					 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
					 componentLongestPath[components[i]]<minPathLength
				   ) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstring FAILED from alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"] "+componentLongestPath[components[i]]+" :: "+minPathLength);						 
					i++;
					continue;
				}
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("calling getSubstring from alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"]");
				getSubstring(i,nAlignments,minPathLength);
				destination=tmp[0];
				pathLength=tmp[1];
				maximalStartA=tmp[2];
				maximalEndA=tmp[3];
				if (destination==-1) {
					i++;
					continue;
				}
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("success! from alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"] maximalEndA="+maximalEndA);
				lastSubstring++;
				substrings[lastSubstring].isWeak=false;
				substrings[lastSubstring].id=0;
				substrings[lastSubstring].startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				substrings[lastSubstring].maximalStartA=maximalStartA;
				substrings[lastSubstring].sumStartA=substrings[lastSubstring].startA;
				substrings[lastSubstring].nStartA=1;
				substrings[lastSubstring].minPrefixLength=-1;
				substrings[lastSubstring].endA=Alignments.alignments[ReadA.sortedAlignments[destination].id][4];
				substrings[lastSubstring].maximalEndA=maximalEndA;
				substrings[lastSubstring].sumEndA=substrings[lastSubstring].endA;
				substrings[lastSubstring].nEndA=1;
				substrings[lastSubstring].minSuffixLength=-1;
				substrings[lastSubstring].prefixReplication=false;
				substrings[lastSubstring].suffixReplication=false;
				substrings[lastSubstring].substringReplication=true;
				substrings[lastSubstring].singleDeletionReplication=false;
				substrings[lastSubstring].leftAlignment=i;
				substrings[lastSubstring].rightAlignment=destination;
				substrings[lastSubstring].pathLength=pathLength;
				substrings[lastSubstring].clearImpliedAlignmentsStats();
				substrings[lastSubstring].clearPrefixSuffixBias();
				substrings[lastSubstring].trim(tmp,MIN_TRIMMED_LENGTH);
				// We don't enter a dense substring of this type
				i++;
				while (i<nAlignments && Alignments.alignments[ReadA.sortedAlignments[i].id][3]<=substrings[lastSubstring].endA) i++;
			}
			if (lastSubstring>=0) {
				DenseSubstring.order=DenseSubstring.STARTA;
				if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
			}

			// Marking all alignments whose interval in readA falls inside a dense
			// substring that replicates by taking substrings of itself.
			markAlignments(nAlignments-1);
		}

		// Building the prefix DAG of all alignments of $readA$ with any $readB$, and
		// sorting it topologically.
		clearNeighborMatrices(nAlignments);
		for (i=0; i<nAlignments; i++) {
			if (ReadA.sortedAlignments[i].inDenseSubstring!=null) continue;
			if (ReadA.sortedAlignments[i].isLeftMaximalB!=1 || ReadA.sortedAlignments[i].isRightMaximalB!=1) continue;
			if (ReadA.sortedAlignments[i].lowQualityEnd) continue;
			if (ReadA.sortedAlignments[i].inPeriodicSubstring) continue;
			startA1=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			endA1=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			for (j=i-1; j>=0; j--) {
				startA2=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
				if (startA2<startA1-identityThreshold) break;
				if (ReadA.sortedAlignments[j].inDenseSubstring!=null) continue;
				if (ReadA.sortedAlignments[j].isLeftMaximalB!=1 || ReadA.sortedAlignments[j].isRightMaximalB!=1) continue;
				if (ReadA.sortedAlignments[j].inPeriodicSubstring) continue;
				endA2=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
				if (endA2<startA1) continue;  // Only alignments that overlap in $readA$ can be connected
				if (endA2>endA1 && endA2-endA1>identityThreshold && endA2-endA1<=maxDistanceEnd) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
					inNeighbors[j][nInNeighbors[j]++]=i;
				}
			}
			for (j=i+1; j<nAlignments; j++) {
				startA2=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
				if (startA2>startA1+identityThreshold) break;
				if (ReadA.sortedAlignments[j].inDenseSubstring!=null) continue;
				if (ReadA.sortedAlignments[j].isLeftMaximalB!=1 || ReadA.sortedAlignments[j].isRightMaximalB!=1) continue;
				if (ReadA.sortedAlignments[j].inPeriodicSubstring) continue;
				if (startA2>endA1) break;  // Only alignments that overlap in $readA$ can be connected
				endA2=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
				if (endA2>endA1 && endA2-endA1>identityThreshold && endA2-endA1<=maxDistanceEnd) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
					inNeighbors[j][nInNeighbors[j]++]=i;
				}
			}
		}
		nComponents=DAG.getConnectedComponents(nAlignments,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stackPrime);
		if (nComponents<=0 || nComponents>nAlignments) {
			IO.printCriticalErr("Error while detecting dense substrings: wrong number of connected components.");
			System.exit(1);
		}
		System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nAlignments);
		i=DAG.topologicalSort(nAlignments,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
		if (i!=0) {
			IO.printCriticalErr("Error while detecting dense substrings: the prefix DAG of all alignments of read "+ReadA.id+" contains a cycle that involves node "+(i-1));
			System.exit(1);
		}
		if (nComponents>lengthPoints.length) {
			tmpPoints = new Point[nComponents];
			System.arraycopy(lengthPoints,0,tmpPoints,0,lengthPoints.length);
			for (i=lengthPoints.length; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
			lengthPoints=tmpPoints;
		}
		longestPathInComponent=DAG.longestPathInComponent(nAlignments,sorted2original,nInNeighbors,inNeighbors,components,lengthPoints,distances);
		if (longestPathInComponent>0) {
			for (i=0; i<nComponents; i++) componentLongestPath[i]=(int)lengthPoints[i].position;
			if (nComponents==1) minPathLength=DEFAULT_MIN_PATH_LENGTH;
			else {
				minPathLength=estimateFromPathLengthPoints(nComponents,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE);
				if (minPathLength<DEFAULT_MIN_PATH_LENGTH) minPathLength=DEFAULT_MIN_PATH_LENGTH;
				if (minPathLength>DEFAULT_MAX_PATH_LENGTH) minPathLength=DEFAULT_MAX_PATH_LENGTH;
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstrings> minPathLength="+minPathLength+" for prefix type");
			}
			if (IO.SHOW_STD_ERR_PRIME) IO.printErr("nComponents in prefix DAG: "+nComponents+" longestPathInComponent="+longestPathInComponent);

			// Detecting all dense substrings of prefix type
			i=0;
			while (i<nAlignments) {
				if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
					 ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
					 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
					 ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
					 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
					 componentLongestPath[components[i]]<minPathLength
				   ) {
					i++;
					continue;
				}
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstring> is there a dense substring of prefix type from alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"]?  (connected component="+components[i]+")");
				getPrefixSubstring(i,nAlignments,minPathLength);
				destination=tmp[0];
				pathLength=tmp[1];
				maximalEndA=tmp[2];
				if (destination==-1) {
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr(" NO!");
					i++;
					continue;
				}
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr(" YES, pathLength="+pathLength);
				lastSubstring++;
				substrings[lastSubstring].isWeak=false;
				substrings[lastSubstring].id=0;
				substrings[lastSubstring].startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				substrings[lastSubstring].maximalStartA=-1;
				substrings[lastSubstring].sumStartA=substrings[lastSubstring].startA;
				substrings[lastSubstring].nStartA=1;
				substrings[lastSubstring].minPrefixLength=Alignments.alignments[ReadA.sortedAlignments[i].id][4]-Alignments.alignments[ReadA.sortedAlignments[i].id][3]+1;
				substrings[lastSubstring].endA=Alignments.alignments[ReadA.sortedAlignments[destination].id][4];
				substrings[lastSubstring].maximalEndA=maximalEndA;
				substrings[lastSubstring].sumEndA=substrings[lastSubstring].endA;
				substrings[lastSubstring].nEndA=1;
				substrings[lastSubstring].minSuffixLength=-1;
				substrings[lastSubstring].prefixReplication=true;
				substrings[lastSubstring].suffixReplication=false;
				substrings[lastSubstring].substringReplication=false;
				substrings[lastSubstring].singleDeletionReplication=false;
				substrings[lastSubstring].leftAlignment=i;
				substrings[lastSubstring].rightAlignment=destination;
				substrings[lastSubstring].pathLength=pathLength;
				substrings[lastSubstring].clearImpliedAlignmentsStats();
				substrings[lastSubstring].clearPrefixSuffixBias();
				nMarkedAlignments=markImpliedAlignmentsPath(substrings[lastSubstring],nAlignments-1,minPathLength);
				if (nMarkedAlignments<minPathLength) {
					// Discarding the substring if not enough alignments are assigned 
					// to it.
					unmarkImpliedAlignments(substrings[lastSubstring],nAlignments-1);
					lastSubstring--;
				}
				else substrings[lastSubstring].trim(tmp,MIN_TRIMMED_LENGTH);
				i++;
			}
		}
		
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("IN GETSUBSTRINGS, BEFORE NON-WEAK SUFFIX>");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) IO.printErr(i+","+ReadA.sortedAlignments[i].toStringBoundaries()+" impliedByPrefixSubstring="+ReadA.sortedAlignments[i].impliedByPrefixSubstring+" impliedBySuffixSubstring="+ReadA.sortedAlignments[i].impliedBySuffixSubstring+" impliedByPrefixSubstringPrime="+ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime+" impliedBySuffixSubstringPrime="+ReadA.sortedAlignments[i].impliedBySuffixSubstring+" inDenseSubstring="+ReadA.sortedAlignments[i].inDenseSubstring);
}		
		
		
		// Building the suffix DAG of all alignments of $readA$ with any $readB$, and
		// sorting it topologically. The DAG includes also alignments that are implied by
		// prefix substrings, and alignments whose length is at a peak prefix length,
		// since an alignment can belong simultaneously to a prefix and a suffix
		// substring.
		clearNeighborMatrices(nAlignments);
		for (i=0; i<nAlignments; i++) {
			if (ReadA.sortedAlignments[i].inDenseSubstring!=null) continue;
			if (ReadA.sortedAlignments[i].isLeftMaximalB!=1 || ReadA.sortedAlignments[i].isRightMaximalB!=1) continue;
			if (ReadA.sortedAlignments[i].inPeriodicSubstring) continue;
			startA1=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			endA1=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			for (j=i+1; j<nAlignments; j++) {
				startA2=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
				if (startA2>startA1+maxDistanceEnd) break;
				if (startA2<=startA1+identityThreshold) continue;
				if (ReadA.sortedAlignments[j].inDenseSubstring!=null) continue;
				if (ReadA.sortedAlignments[j].isLeftMaximalB!=1 || ReadA.sortedAlignments[j].isRightMaximalB!=1) continue;
				if (ReadA.sortedAlignments[j].lowQualityStart) continue;
				if (ReadA.sortedAlignments[j].inPeriodicSubstring) continue;
				endA2=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
				if (Intervals.intersectionLength(startA1,endA1,startA2,endA2)==0) continue;  // Only alignments that overlap in $readA$ can be connected
				if (Math.abs(endA2,endA1)<=identityThreshold) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
					inNeighbors[j][nInNeighbors[j]++]=i;
				}
			}
		}
		nComponents=DAG.getConnectedComponents(nAlignments,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stackPrime);
		if (nComponents<=0 || nComponents>nAlignments) {
			IO.printCriticalErr("Error while detecting dense substrings: wrong number of connected components.");
			System.exit(1);
		}
		System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nAlignments);
		i=DAG.topologicalSort(nAlignments,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
		if (i!=0) {
			IO.printCriticalErr("Error while detecting dense substrings: the suffix DAG of all alignments of read "+ReadA.id+" contains a cycle that involves node "+(i-1));
			System.exit(1);
		}
		if (nComponents>lengthPoints.length) {
			tmpPoints = new Point[nComponents];
			System.arraycopy(lengthPoints,0,tmpPoints,0,lengthPoints.length);
			for (i=lengthPoints.length; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
			lengthPoints=tmpPoints;
		}
		longestPathInComponent=DAG.longestPathInComponent(nAlignments,sorted2original,nInNeighbors,inNeighbors,components,lengthPoints,distances);
		if (longestPathInComponent>0) {
			for (i=0; i<nComponents; i++) {
				componentLongestPath[i]=(int)lengthPoints[i].position;
if (IO.SHOW_STD_ERR) IO.printErr("suffix component "+i+": longestPath="+componentLongestPath[i]);
			}
			if (nComponents==1) minPathLength=DEFAULT_MIN_PATH_LENGTH;
			else {
				minPathLength=estimateFromPathLengthPoints(nComponents,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE);
				if (minPathLength<DEFAULT_MIN_PATH_LENGTH) minPathLength=DEFAULT_MIN_PATH_LENGTH;
				if (minPathLength>DEFAULT_MAX_PATH_LENGTH) minPathLength=DEFAULT_MAX_PATH_LENGTH;
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstrings> minPathLength="+minPathLength+" for suffix type");
			}
			if (IO.SHOW_STD_ERR) IO.printErr("nComponents in suffix DAG: "+nComponents+" longestPathInComponent="+longestPathInComponent);

			// Detecting all dense substrings of suffix type
			i=0;
			while (i<nAlignments) {
if (componentLongestPath[components[i]]<minPathLength) {
	if (IO.SHOW_STD_ERR_PRIME) {
		IO.printErr("getSubstring> suffix connected component too small of alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"], "+componentLongestPath[components[i]]+"<"+minPathLength);
		IO.printErr("getSubstring> inDenseSubstring? "+ReadA.sortedAlignments[i].inDenseSubstring);
		IO.printErr("getSubstring> isLeftMaximalB? "+ReadA.sortedAlignments[i].isLeftMaximalB+" isRightMaximalB? "+ReadA.sortedAlignments[i].isRightMaximalB);
	}
}
				if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
					 (nInNeighbors[i]>0 && ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null) ||  // Always using the minimal vertices of the DAG: otherwise, $minPathLength$ might not be achievable.
					 (nInNeighbors[i]>0 && ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null) ||
					 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
					 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
					 ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
					 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
					 componentLongestPath[components[i]]<minPathLength ) {
					i++;
					continue;
				}
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstring> is there a dense substring of suffix type from alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"], component "+components[i]+"?");
				getSuffixSubstring(i,nAlignments,minPathLength);
				destination=tmp[0];
				pathLength=tmp[1];
				maximalStartA=tmp[2];
				if (destination==-1) {
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr(" NO!");
					i++;
					continue;
				}
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr(" YES");
				lastSubstring++;
				substrings[lastSubstring].isWeak=false;
				substrings[lastSubstring].id=0;
				substrings[lastSubstring].startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				substrings[lastSubstring].maximalStartA=maximalStartA;
				substrings[lastSubstring].sumStartA=substrings[lastSubstring].startA;
				substrings[lastSubstring].nStartA=1;
				substrings[lastSubstring].minPrefixLength=-1;
				substrings[lastSubstring].endA=Alignments.alignments[ReadA.sortedAlignments[destination].id][4];
				substrings[lastSubstring].maximalEndA=-1;
				substrings[lastSubstring].sumEndA=substrings[lastSubstring].endA;
				substrings[lastSubstring].nEndA=1;
				substrings[lastSubstring].minSuffixLength=Alignments.alignments[ReadA.sortedAlignments[destination].id][4]-Alignments.alignments[ReadA.sortedAlignments[destination].id][3]+1;
				substrings[lastSubstring].prefixReplication=false;
				substrings[lastSubstring].suffixReplication=true;
				substrings[lastSubstring].substringReplication=false;
				substrings[lastSubstring].singleDeletionReplication=false;
				substrings[lastSubstring].leftAlignment=i;
				substrings[lastSubstring].rightAlignment=destination;
				substrings[lastSubstring].pathLength=pathLength;
				substrings[lastSubstring].clearImpliedAlignmentsStats();
				substrings[lastSubstring].clearPrefixSuffixBias();
				nMarkedAlignments=markImpliedAlignmentsPath(substrings[lastSubstring],nAlignments-1,minPathLength);
				if (nMarkedAlignments<minPathLength) {
					// Discarding the substring if not enough alignments are assigned 
					// to it.
					unmarkImpliedAlignments(substrings[lastSubstring],nAlignments-1);
					lastSubstring--;
				}
				else substrings[lastSubstring].trim(tmp,MIN_TRIMMED_LENGTH);
				i++;
			}
		}
		
		
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("IN GETSUBSTRINGS, BEFORE WEAK PREFIX AND BEFORE EXTEND>");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) IO.printErr(i+","+ReadA.sortedAlignments[i].toStringBoundaries()+" impliedByPrefixSubstring="+ReadA.sortedAlignments[i].impliedByPrefixSubstring+" inDense="+ReadA.sortedAlignments[i].inDenseSubstring);
}		
		
		
		// Extending to the left and to the right all dense substrings of substring type
		nAlignments=extendSubstrings(nAlignments);





if (IO.SHOW_INTERACTIVE) {
	IO.printErr("IN GETSUBSTRINGS, BEFORE WEAK PREFIX>");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"], leftMax="+ReadA.sortedAlignments[i].isLeftMaximal+" rightMax="+ReadA.sortedAlignments[i].isRightMaximal+" || "+(ReadA.sortedAlignments[i].impliedByDenseSubstring!=null?ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode():"null")+","+(ReadA.sortedAlignments[i].inDenseSubstring!=null?ReadA.sortedAlignments[i].inDenseSubstring.hashCode():"null")+" impliedByPrefixSubstring="+ReadA.sortedAlignments[i].impliedByPrefixSubstring+" impliedBySuffixSubstring="+ReadA.sortedAlignments[i].impliedBySuffixSubstring+" PERIODIC="+ReadA.sortedAlignments[i].impliedByPeriodicSubstring);
	}
}





		// Building the weak prefix DAG of all alignments of $readA$, with any $readB$,
		// and sorting it topologically. The DAG includes also alignments that are
		// implied by suffix substrings, and alignments whose length is at a peak suffix
		// length, since an alignment can belong simultaneously to a prefix and a suffix
		// substring.
		tmpOrder=Alignment.order;
		Alignment.order=Alignment.ENDA;
		if (nAlignments>1) Arrays.sort(ReadA.sortedAlignments,0,nAlignments);
		clearNeighborMatrices(nAlignments);
		for (i=0; i<nAlignments; i++) {
			if ( ReadA.sortedAlignments[i].inDenseSubstring!=null || 
				 ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null || 
				 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null || 
				 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
				 ReadA.sortedAlignments[i].lowQualityEnd ||
				 ReadA.sortedAlignments[i].inPeriodicSubstring
			   ) continue;
			startA1=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			endA1=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			for (j=i+1; j<nAlignments; j++) {
				startA2=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
				endA2=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
				if (endA2>endA1+maxDistanceEndWeak) break;
				if (endA2<=endA1+identityThreshold || Math.abs(startA2,startA1)>maxDistanceEndWeakPrime) continue;
				if ( ReadA.sortedAlignments[j].inDenseSubstring!=null || 
					 ReadA.sortedAlignments[j].impliedByPrefixSubstring!=null || 
					 ReadA.sortedAlignments[j].impliedByPrefixSubstringPrime!=null || 
					 ReadA.sortedAlignments[j].isRightMaximalB!=1 ||
					 ReadA.sortedAlignments[j].inPeriodicSubstring
				   ) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(nAlignments,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stackPrime);
		if (nComponents<=0 || nComponents>nAlignments) {
			IO.printCriticalErr("Error while detecting weak dense substrings: wrong number of connected components.");
			System.exit(1);
		}		
		System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nAlignments);
		i=DAG.topologicalSort(nAlignments,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
		if (i!=0) {
			IO.printCriticalErr("Error while detecting dense substrings: the weak prefix DAG of all alignments of read "+ReadA.id+" contains a cycle that involves node "+(i-1));
			System.exit(1);
		}
		if (nComponents>lengthPoints.length) {
			tmpPoints = new Point[nComponents];
			System.arraycopy(lengthPoints,0,tmpPoints,0,lengthPoints.length);
			for (i=lengthPoints.length; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
			lengthPoints=tmpPoints;
		}
		longestPathInComponent=DAG.longestPathInComponent(nAlignments,sorted2original,nInNeighbors,inNeighbors,components,lengthPoints,distances);
		if (longestPathInComponent>0) {
			for (i=0; i<nComponents; i++) componentLongestPath[i]=(int)lengthPoints[i].position;
			if (nComponents==1) minPathLength=DEFAULT_MIN_PATH_LENGTH;
			else {
				minPathLength=estimateFromPathLengthPoints(nComponents,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE);
				if (minPathLength<DEFAULT_MIN_PATH_LENGTH) minPathLength=DEFAULT_MIN_PATH_LENGTH;
				if (minPathLength>DEFAULT_MAX_PATH_LENGTH) minPathLength=DEFAULT_MAX_PATH_LENGTH;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstrings> minPathLength="+minPathLength+" for weak prefix type");
			}
			if (IO.SHOW_STD_ERR) IO.printErr("nComponents in weak prefix DAG: "+nComponents+" longestPathInComponent="+longestPathInComponent);

			// Detecting all weak dense substrings of prefix type
			i=0;
			while (i<nAlignments) {
				if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
					 ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
					 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
					 (nInNeighbors[i]>0 && ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null) ||  // Always using the minimal vertices of the DAG: otherwise, $minPathLength$ might not be achievable.
					 (nInNeighbors[i]>0 && ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null) ||
					 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
					 componentLongestPath[components[i]]<minPathLength ) {
					i++;
					continue;
				}
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstring> is there a dense substring of weak prefix type from alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"]?  (connected component="+components[i]+")");
				getPrefixSubstring_weak(i,nAlignments,minPathLength);
				destination=tmp[0];
				if (destination==-1) {
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr(" NO!");
					i++;
					continue;
				}
				pathLength=tmp[1];
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr(" YES, pathLength="+pathLength);
				maximalEndA=tmp[2];
				minBLeftMaximal=tmp[3];
				maxBLeftMaximal=tmp[4];
				leftAlignment=tmp[5];
				minBLeft=tmp[6];
				lastSubstring++;
				substrings[lastSubstring].isWeak=true;
				substrings[lastSubstring].id=0;
				substrings[lastSubstring].destinationAlignment=destination;
				substrings[lastSubstring].startA=minBLeft;
				substrings[lastSubstring].maximalStartA=-1;
				substrings[lastSubstring].sumStartA=substrings[lastSubstring].startA;
				substrings[lastSubstring].nStartA=1;
				substrings[lastSubstring].endA=Alignments.alignments[ReadA.sortedAlignments[destination].id][4];
				substrings[lastSubstring].maximalEndA=maximalEndA;
				substrings[lastSubstring].sumEndA=substrings[lastSubstring].endA;
				substrings[lastSubstring].nEndA=1;
				if ( minBLeftMaximal==-1 || 
					 (Math.abs(minBLeftMaximal,minBLeft)<=identityThreshold && maxBLeftMaximal-minBLeftMaximal<=identityThreshold<<1) ) {
					substrings[lastSubstring].prefixReplication=true;
					substrings[lastSubstring].minPrefixLength=Alignments.alignments[ReadA.sortedAlignments[i].id][4]-substrings[lastSubstring].startA+1;
					substrings[lastSubstring].substringReplication=false;
				}
				else {
					substrings[lastSubstring].prefixReplication=false;
					substrings[lastSubstring].minPrefixLength=-1;
					substrings[lastSubstring].substringReplication=true;
				}
				substrings[lastSubstring].suffixReplication=false;
				substrings[lastSubstring].minSuffixLength=-1;
				substrings[lastSubstring].singleDeletionReplication=false;
				substrings[lastSubstring].leftAlignment=leftAlignment;
				substrings[lastSubstring].rightAlignment=destination;
				substrings[lastSubstring].pathLength=pathLength;
				substrings[lastSubstring].clearImpliedAlignmentsStats();
				substrings[lastSubstring].clearPrefixSuffixBias();
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("detected weak prefix substring ["+substrings[lastSubstring].startA+".."+substrings[lastSubstring].endA+"] prefixReplication="+substrings[lastSubstring].prefixReplication+" substringReplication="+substrings[lastSubstring].substringReplication);
				nMarkedAlignments=markImpliedAlignmentsPath_weak(substrings[lastSubstring],nAlignments-1,true);
				if (nMarkedAlignments<minPathLength) {
					// Discarding the substring if not enough alignments are assigned 
					// to it.
					unmarkImpliedAlignments_weak_prefix(substrings[lastSubstring],nAlignments-1);
					lastSubstring--;
				}
				else substrings[lastSubstring].trim(tmp,MIN_TRIMMED_LENGTH);
				i++;
			}
		}
		Alignment.order=tmpOrder;
		if (nAlignments>1) Arrays.sort(ReadA.sortedAlignments,0,nAlignments);




if (IO.SHOW_INTERACTIVE) {
	IO.printErr("IN GETSUBSTRINGS, BEFORE WEAK SUFFIX>");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"], leftMax="+ReadA.sortedAlignments[i].isLeftMaximal+" rightMax="+ReadA.sortedAlignments[i].isRightMaximal+" || "+(ReadA.sortedAlignments[i].impliedByDenseSubstring!=null?ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode():"null")+","+(ReadA.sortedAlignments[i].inDenseSubstring!=null?ReadA.sortedAlignments[i].inDenseSubstring.hashCode():"null")+"!!"+(ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null?ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hashCode():"null") );
	}
}


		// Building the weak suffix DAG of all alignments of $readA$, with any $readB$,
		// and sorting it topologically. The DAG includes also alignments that are
		// implied by prefix substrings, and alignments whose length is at a peak prefix
		// length, since an alignment can belong simultaneously to a prefix and a suffix
		// substring.
		clearNeighborMatrices(nAlignments);
		for (i=0; i<nAlignments; i++) {
			if ( ReadA.sortedAlignments[i].inDenseSubstring!=null || 
				 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
				 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
				 ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
				 ReadA.sortedAlignments[i].inPeriodicSubstring
			   ) continue;
			startA1=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			endA1=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			for (j=i+1; j<nAlignments; j++) {
				startA2=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
				endA2=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
				if (startA2>startA1+maxDistanceEndWeak) break;
				if (startA2<=startA1+identityThreshold || Math.abs(endA2,endA1)>maxDistanceEndWeakPrime) continue;
				if (ReadA.sortedAlignments[j].lowQualityStart) continue;
				if ( ReadA.sortedAlignments[j].inDenseSubstring!=null || 
					 ReadA.sortedAlignments[j].impliedBySuffixSubstring!=null ||
					 ReadA.sortedAlignments[j].impliedBySuffixSubstringPrime!=null ||
					 ReadA.sortedAlignments[j].isLeftMaximalB!=1 ||
					 ReadA.sortedAlignments[j].inPeriodicSubstring
				   ) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(nAlignments,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stackPrime);
		if (nComponents<=0 || nComponents>nAlignments) {
			IO.printCriticalErr("Error while detecting weak dense substrings: wrong number of connected components.");
			System.exit(1);
		}		
		System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nAlignments);
		i=DAG.topologicalSort(nAlignments,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
		if (i!=0) {
			IO.printCriticalErr("Error while detecting weak dense substrings: the weak suffix DAG of all alignments of read "+ReadA.id+" contains a cycle that involves node "+(i-1));
			System.exit(1);
		}
		if (nComponents>lengthPoints.length) {
			tmpPoints = new Point[nComponents];
			System.arraycopy(lengthPoints,0,tmpPoints,0,lengthPoints.length);
			for (i=lengthPoints.length; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
			lengthPoints=tmpPoints;
		}
		longestPathInComponent=DAG.longestPathInComponent(nAlignments,sorted2original,nInNeighbors,inNeighbors,components,lengthPoints,distances);
		if (longestPathInComponent>0) {
			for (i=0; i<nComponents; i++) componentLongestPath[i]=(int)lengthPoints[i].position;
			if (nComponents==1) minPathLength=DEFAULT_MIN_PATH_LENGTH;
			else {
				minPathLength=estimateFromPathLengthPoints(nComponents,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE);
				if (minPathLength<DEFAULT_MIN_PATH_LENGTH) minPathLength=DEFAULT_MIN_PATH_LENGTH;
				if (minPathLength>DEFAULT_MAX_PATH_LENGTH) minPathLength=DEFAULT_MAX_PATH_LENGTH;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstrings> minPathLength="+minPathLength+" for weak suffix type");
			}
			if (IO.SHOW_STD_ERR_PRIME) IO.printErr("nComponents in weak suffix DAG: "+nComponents+" longestPathInComponent="+longestPathInComponent);

			// Detecting all weak dense substrings of suffix type
			i=0;
			while (i<nAlignments) {
				if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
					 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
					 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
					 (nInNeighbors[i]>0 && ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null) ||  // Always using the minimal vertices of the DAG: otherwise, $minPathLength$ might not be achievable.
					 (nInNeighbors[i]>0 && ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null) ||
					 ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
					 componentLongestPath[components[i]]<minPathLength ) {
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstring> weak suffix type not possible from alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"] "+ReadA.sortedAlignments[i].inDenseSubstring+" "+(ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null)+" "+(ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null)+" "+(ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null)+" "+(ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null)+" "+(ReadA.sortedAlignments[i].isLeftMaximalB!=1)+" componentLongestPath="+componentLongestPath[components[i]]+", minPathLength="+minPathLength); 
					i++;
					continue;
				}
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getSubstring> is there a dense substring of weak suffix type from alignment ["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"]?  (connected component="+components[i]+")");
				getSuffixSubstring_weak(i,nAlignments,minPathLength);
				destination=tmp[0];
				if (destination==-1) {
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr(" NO!");
					i++;
					continue;
				}
				pathLength=tmp[1];
				maximalStartA=tmp[2];
				minBRightMaximal=tmp[3];
				maxBRightMaximal=tmp[4];
				rightAlignment=tmp[5];
				maxBRight=tmp[6];
if (IO.SHOW_STD_ERR_PRIME) IO.printErr(" YES, pathLength="+pathLength+" maximalStartA="+maximalStartA+" minBRightMaximal="+minBRightMaximal+" maxBRightMaximal="+maxBRightMaximal+" rightAlignment="+rightAlignment+" maxBRight="+maxBRight);
				lastSubstring++;
				substrings[lastSubstring].isWeak=true;
				substrings[lastSubstring].id=0;
				substrings[lastSubstring].destinationAlignment=destination;
				substrings[lastSubstring].endA=maxBRight;
				substrings[lastSubstring].maximalEndA=-1;
				substrings[lastSubstring].sumEndA=substrings[lastSubstring].endA;
				substrings[lastSubstring].nEndA=1;
				substrings[lastSubstring].startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				substrings[lastSubstring].maximalStartA=maximalStartA;
				substrings[lastSubstring].sumStartA=substrings[lastSubstring].startA;
				substrings[lastSubstring].nStartA=1;
				if ( maxBRightMaximal==-1 || 
					 (Math.abs(maxBRightMaximal,maxBRight)<=identityThreshold && maxBRightMaximal-minBRightMaximal<=identityThreshold<<1) ) {
					substrings[lastSubstring].suffixReplication=true;
					substrings[lastSubstring].minSuffixLength=substrings[lastSubstring].endA-Alignments.alignments[ReadA.sortedAlignments[destination].id][3]+1;
					substrings[lastSubstring].substringReplication=false;
				}
				else {
					substrings[lastSubstring].suffixReplication=false;
					substrings[lastSubstring].minSuffixLength=-1;
					substrings[lastSubstring].substringReplication=true;
				}
				substrings[lastSubstring].prefixReplication=false;
				substrings[lastSubstring].minPrefixLength=-1;
				substrings[lastSubstring].singleDeletionReplication=false;
				substrings[lastSubstring].rightAlignment=rightAlignment;
				substrings[lastSubstring].leftAlignment=i;
				substrings[lastSubstring].pathLength=pathLength;
				substrings[lastSubstring].clearImpliedAlignmentsStats();
				substrings[lastSubstring].clearPrefixSuffixBias();
				nMarkedAlignments=markImpliedAlignmentsPath_weak(substrings[lastSubstring],nAlignments-1,false);
				if (nMarkedAlignments<minPathLength) {
					// Discarding the substring if not enough alignments are assigned 
					// to it.
					unmarkImpliedAlignments_weak_suffix(substrings[lastSubstring],nAlignments-1);
					lastSubstring--;
				}
				else substrings[lastSubstring].trim(tmp,MIN_TRIMMED_LENGTH);
				i++;
			}
		}
		
		
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("IN GETSUBSTRINGS, AFTER WEAK SUFFIX>");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"], leftMax="+ReadA.sortedAlignments[i].isLeftMaximal+" rightMax="+ReadA.sortedAlignments[i].isRightMaximal+" || "+(ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null?ReadA.sortedAlignments[i].impliedByPrefixSubstring.hashCode():"null")+","+(ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null?ReadA.sortedAlignments[i].impliedBySuffixSubstring.hashCode():"null") );
	}
}		
		
		
		// Discarding alignments after trimming (executed only once, since the removed
	 	// alignments likely contain low-quality regions, and thus are not likely to be
		// useful for building other substrings).
		discardAlignmentsAfterTrimming(0,nAlignments-1);
		

if (IO.SHOW_INTERACTIVE) {
	IO.printErr("IN GETSUBSTRINGS, AFTER discardAlignmentsAfterTrimming>");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"], leftMax="+ReadA.sortedAlignments[i].isLeftMaximal+" rightMax="+ReadA.sortedAlignments[i].isRightMaximal+" || "+(ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null?ReadA.sortedAlignments[i].impliedByPrefixSubstring.hashCode():"null")+","+(ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null?ReadA.sortedAlignments[i].impliedBySuffixSubstring.hashCode():"null") );
	}
}

		
		
		
		// Enforcing straddling/containment with dense substrings of substring type
		cleanDenseSubstringsOfSubstringType(nAlignments-1);



if (IO.SHOW_INTERACTIVE) {
	IO.printErr("IN GETSUBSTRINGS, AFTER CLEAN>");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+"\n");
	IO.printErr("alignments:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"], leftMax="+ReadA.sortedAlignments[i].isLeftMaximal+" rightMax="+ReadA.sortedAlignments[i].isRightMaximal+" || "+(ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null?ReadA.sortedAlignments[i].impliedByPrefixSubstring.hashCode():"null")+","+(ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null?ReadA.sortedAlignments[i].impliedBySuffixSubstring.hashCode():"null") );
	}
}



		// Sorting dense substrings by startA
		DenseSubstring.order=DenseSubstring.STARTA;
		if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);

		// Setting field $impliedByDenseSubstring$ of every alignment. If an alignment is
		// implied both by a prefix and by a suffix substring, it is assigned to the
		// substring for which the alignment is the only one enforcing a boundary.
		// If this cannot be established, the alignment is assigned to the shortest
		// substring, or, if they have similar lengths, to the one that might be part of
		// a single-deletion substring (to be detected later). In case it is not possible
		// to decide, the alignment is assigned to one substring or the other at random.
		//
		// Remark: after this phase completes, fields $impliedBy{Prefix,Suffix}Substring$
	 	// and $impliedBy{Prefix,Suffix}SubstringPrime$ of an alignment are not
		// meaningful any more, and are all reset to NULL. The only field that is kept in
		// the following stages of the pipeline is $impliedByDenseSubstring$.
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].nextSuffixSubstring=null;
			for (j=i+1; j<=lastSubstring; j++) {
				if (substrings[j].suffixReplication) {
					substrings[i].nextSuffixSubstring=substrings[j];
					break;
				}
			}
		}
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].previousPrefixSubstring=null;
			for (j=i-1; j>=0; j--) {
				if (substrings[j].prefixReplication) {
					substrings[i].previousPrefixSubstring=substrings[j];
					break;
				}
			}
		}
		for (i=0; i<nAlignments; i++) {
			ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
			prefixSubstring=ReadA.sortedAlignments[i].impliedByPrefixSubstring;
			suffixSubstring=ReadA.sortedAlignments[i].impliedBySuffixSubstring;
			if (prefixSubstring==null && suffixSubstring==null) continue;
			if (prefixSubstring!=null && suffixSubstring==null) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=prefixSubstring;
				continue;
			}
			else if (prefixSubstring==null && suffixSubstring!=null) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=suffixSubstring;
				continue;
			}
			closeToPrefix=Math.abs(ReadA.sortedAlignments[i].endA(),prefixSubstring.endA)<=identityThreshold;
			closeToSuffix=Math.abs(ReadA.sortedAlignments[i].startA(),suffixSubstring.startA)<=identityThreshold;
			if (closeToPrefix && !closeToSuffix) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=prefixSubstring;
				continue;
			}
			else if (!closeToPrefix && closeToSuffix) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=suffixSubstring;
				continue;
			}
			else if (closeToPrefix && closeToSuffix) {
				j=preserveBoundary(i,prefixSubstring,suffixSubstring);
				if (j==0) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=prefixSubstring;
					continue;
				}
				else if (j==1) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=suffixSubstring;
					continue;
				}
			}
			closeToPrefix=Math.abs(ReadA.sortedAlignments[i].startA(),prefixSubstring.startA)<=IDENTITY_THRESHOLD_AMBIGUOUS;
			closeToSuffix=Math.abs(ReadA.sortedAlignments[i].endA(),suffixSubstring.endA)<=IDENTITY_THRESHOLD_AMBIGUOUS;
			if (closeToPrefix && !closeToSuffix) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=prefixSubstring;
				continue;
			}
			else if (!closeToPrefix && closeToSuffix) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=suffixSubstring;
				continue;
			}
			lengthPrefix=prefixSubstring.endA-prefixSubstring.startA+1;
			lengthSuffix=suffixSubstring.endA-suffixSubstring.startA+1;
			if (Math.abs(lengthPrefix,lengthSuffix)>=MIN_LENGTH_DIFFERENCE) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=lengthPrefix<lengthSuffix?prefixSubstring:suffixSubstring;
				continue;
			}
			if ( ( prefixSubstring.nextSuffixSubstring==null ||
			       (prefixSubstring.nextSuffixSubstring==suffixSubstring && suffixSubstring.nextSuffixSubstring==null)
			     ) &&
			     ( suffixSubstring.previousPrefixSubstring!=null &&
			       (suffixSubstring.previousPrefixSubstring!=prefixSubstring || prefixSubstring.previousPrefixSubstring!=null)
			     )
			   ) ReadA.sortedAlignments[i].impliedByDenseSubstring=suffixSubstring;
			else if ( ( suffixSubstring.previousPrefixSubstring==null ||
					    (suffixSubstring.previousPrefixSubstring==prefixSubstring && prefixSubstring.previousPrefixSubstring==null)
					  ) &&
					  ( prefixSubstring.nextSuffixSubstring!=null &&
					    (prefixSubstring.nextSuffixSubstring!=suffixSubstring || suffixSubstring.nextSuffixSubstring!=null)
					  )
				   ) ReadA.sortedAlignments[i].impliedByDenseSubstring=prefixSubstring;
			else {  // Ambiguity
				if (Math.random.nextBoolean()) ReadA.sortedAlignments[i].impliedByDenseSubstring=prefixSubstring;
				else ReadA.sortedAlignments[i].impliedByDenseSubstring=suffixSubstring;
			}
		}
		for (i=0; i<nAlignments; i++) {
			ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
			ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
			ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;
			ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
		}

		// Sorting $lengthBiasIntervals$
		if (lastLengthBiasInterval>0) {
			sortLengthBiasIntervals();
			tmpInterval.lastPosition=-1;
			tmpInterval.readB=lengthBiasIntervals[0].readB;
			tmpInterval.orientation=lengthBiasIntervals[0].orientation;
		}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("lengthBiasIntervals:");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastLengthBiasInterval; i++) IO.printErr("["+lengthBiasIntervals[i].firstPosition+".."+lengthBiasIntervals[i].lastPosition+","+lengthBiasIntervals[i].center+"], "); }
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("");

		return nAlignments;
	}


	/**
	 * Assume that the $alignmentID$-th alignment in $ReadA.sortedAlignments$ is implied
	 * by both $prefixSubstring$ and $suffixSubstring$, and that its first position is
	 * close to the first position of $suffixSubstring$, and its last position is close to
	 * the last position of $prefixSubstring$. The procedure looks for other alignments,
	 * implied by those substrings, that end close to the corresponding positions.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments$ to be sorted by $startA$.
	 *
	 * @return 0: $alignmentID$ should be assigned to the prefix substring; 1: to the 
	 * suffix substring; 2: no substring boundary would be preserved after removing 
	 * $alignmentID$; 3: both substring boundaries would be preserved.
	 */
	private static final int preserveBoundary(int alignmentID, DenseSubstring prefixSubstring, DenseSubstring suffixSubstring) {
		boolean preservedPrefix, preservedSuffix;
		int i, j;
		int startA, endA;
		
		// Prefix
		startA=prefixSubstring.startA; endA=prefixSubstring.endA; preservedPrefix=false;
		for (i=alignmentID-1; i>=0; i--) {
			if (ReadA.sortedAlignments[i].startA()<startA-identityThreshold) break;
			if (ReadA.sortedAlignments[i].impliedByPrefixSubstring!=prefixSubstring && !isAlignmentImpliedByPrefix(i,prefixSubstring,false)) continue;
			if (ReadA.sortedAlignments[i].isRightMaximalB==1 && Math.abs(ReadA.sortedAlignments[i].endA(),endA)<=identityThreshold) {
				preservedPrefix=true;
				break;
			}
		}
		if (!preservedPrefix) {
			for (i=alignmentID+1; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].startA()>=endA) break;
				if (ReadA.sortedAlignments[i].impliedByPrefixSubstring!=prefixSubstring && !isAlignmentImpliedByPrefix(i,prefixSubstring,false)) continue;
				if (ReadA.sortedAlignments[i].isRightMaximalB==1 && Math.abs(ReadA.sortedAlignments[i].endA(),endA)<=identityThreshold) {
					preservedPrefix=true;
					break;
				}
			}
		}
		
		// Suffix
		startA=suffixSubstring.startA; endA=suffixSubstring.endA; preservedSuffix=false;
		for (i=alignmentID-1; i>=0; i--) {
			if (ReadA.sortedAlignments[i].startA()<startA-identityThreshold) break;
			if (ReadA.sortedAlignments[i].impliedBySuffixSubstring!=suffixSubstring && !isAlignmentImpliedBySuffix(i,suffixSubstring,false)) continue;
			if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && Math.abs(ReadA.sortedAlignments[i].startA(),startA)<=identityThreshold) {
				preservedSuffix=true;
				break;
			}
		}
		if (!preservedSuffix) {
			for (i=alignmentID+1; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].startA()>startA+identityThreshold) break;
				if (ReadA.sortedAlignments[i].impliedBySuffixSubstring!=suffixSubstring && !isAlignmentImpliedBySuffix(i,suffixSubstring,false)) continue;
				if (ReadA.sortedAlignments[i].isLeftMaximalB==1) { 
					preservedSuffix=true;
					break;
				}
			}
		}
		
		if (!preservedPrefix && preservedSuffix) return 0;
		else if (preservedPrefix && !preservedSuffix) return 1;
		else if (!preservedPrefix && !preservedSuffix) return 2;
		else return 3;
	}
	

	/**
	 * Fits a density estimation tree on $lengthPoints$, and returns a separation point
	 * between the first and the second local maximum. The first local maximum is assumed
	 * to correspond to the length of longest paths of connected components that arise by
	 * chance.
	 *
	 * Remark: components with longest path zero are not used.
	 */
	private static final int estimateFromPathLengthPoints(int nComponents, double minIntervalLength, double minLocalMaxDistance) {
		final int MIN_PATH_LENGTH = 3;
		int from, tmp, nLocalMaximumLeaves, lastComponent;

		lastComponent=Points.sortAndCompact(lengthPoints,nComponents-1);
		if (lengthPoints[lastComponent].position<=DEFAULT_MIN_PATH_LENGTH) return DEFAULT_MIN_PATH_LENGTH;
		from=lengthPoints[0].position==0.0?1:0;

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("estimateFromPathLengthPoints:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=lastComponent; x++) IO.printErr(lengthPoints[x].position+","+lengthPoints[x].getMass()); }

		if (from>=lastComponent) return DEFAULT_MIN_PATH_LENGTH;
		if (lastComponent+1<=Points.FEW_POINTS) {
			if (lengthPoints[lastComponent].position-lengthPoints[from].position<=MAX_PATHLENGTH_DIFFERENCE) return DEFAULT_MIN_PATH_LENGTH;
			tmp=Points.getRoughThreshold(lengthPoints,lastComponent,true,MIN_PATH_LENGTH,false);
			return tmp==-1?DEFAULT_MIN_PATH_LENGTH:(int)lengthPoints[tmp].position;
		}
		if (Points.areUniformlyDistributed(lengthPoints,from,lastComponent,true,(lengthPoints[lastComponent].position-lengthPoints[from].position)/Points.DEFAULT_NBINS)) return DEFAULT_MIN_PATH_LENGTH;
		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(lengthPoints,from,lastComponent,minIntervalLength,minLocalMaxDistance,true,-1,-1,-1,false,true);
		if (nLocalMaximumLeaves==0) return DEFAULT_MIN_PATH_LENGTH;
		DensityEstimationTree.markRunsOfLocalMaximumLeaves(lengthPoints);


if (IO.SHOW_STD_ERR) IO.printErr("DET:");
if (IO.SHOW_STD_ERR) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr(DensityEstimationTree.leaves[x]); }


		return (int)DensityEstimationTree.separateRun(0,lengthPoints,from,lastComponent,false);  // The last parameter is false, since it doesn't make sense to make the threshold larger by putting it in the middle of a gap rather than at its first position.
	}


	/**
	 * Searches for a maximal dense substring, of prefix type, that starts from alignment
	 * $ReadA.sortedAlignments[source]$. This is done by computing the longest path, of
	 * length at least $minPathLength$, between $ReadA.sortedAlignments[source]$
	 * and every other vertex in the prefix DAG, such that all alignments in the path are
	 * both B-left- and B-right-maximal, and such that the distance between the last
	 * position of the destination alignment and the last position of the source alignment
	 * is greater than $identityThreshold$.
	 *
	 * Remark: the procedure assumes that function $setMaximality$ has already been called
	 * on all elements of $ReadA.sortedAlignments$.
	 *
	 * @return $tmp[0]$: the position in $ReadA.sortedAlignments[0..nVertices-1]$ of an
	 * interval that maximizes the distance from $source$ in the prefix DAG, and that ends
	 * the farthest from $source$, or -1 if no valid path from $source$ is found;
	 * $tmp[1]$: the length of a longest path from $source$ in the prefix DAG;
 	 * $tmp[2]$: the rightmost end of an alignment that is reachable from $source$ in the 
 	 * DAG, that ends at most at the end of $tmp[0]$, and that is B-right-maximal.
	 */
	private static final void getPrefixSubstring(int source, int nVertices, int minPathLength) {
		boolean found;
		int i, j, k, s;
		int state, sourceStartA, sourceEndA, sourceComponent, vertex, rightmostAlignment, distance, longestDistance;
		int start, end, rightmostEnd, rightmostRightMaximalEnd, length;
		int lastSortedVertex;

		// Computing the longest distance from $source$ to every other vertex
		Math.set(distances,nVertices-1,-1);
		Math.set(predecessors,nVertices-1,-1);
		sourceStartA=Alignments.alignments[ReadA.sortedAlignments[source].id][3];
		sourceComponent=components[source];
		distances[source]=0;
		s=original2sorted[source];
 		lastSortedVertex=s;
		for (j=s+1; j<nVertices; j++) {
			vertex=sorted2original[j];
			if (components[vertex]!=sourceComponent) {
				lastSortedVertex=j-1;
				break;
			}
			if (ReadA.sortedAlignments[vertex].isLeftMaximalB!=1 || ReadA.sortedAlignments[vertex].isRightMaximalB!=1) continue;
			if (Math.abs(Alignments.alignments[ReadA.sortedAlignments[vertex].id][3],sourceStartA)>identityThreshold) continue;
			for (i=0; i<nInNeighbors[vertex]; i++) {
				k=inNeighbors[vertex][i];
				if (ReadA.sortedAlignments[k].isLeftMaximalB!=1 || ReadA.sortedAlignments[k].isRightMaximalB!=1) continue;
				if (distances[k]!=-1) distance=distances[k];
				else continue;
				distance++;
				if (distance>distances[vertex]) {
					distances[vertex]=distance;
					predecessors[vertex]=k;
				}
			}
		}
		if (j==nVertices) lastSortedVertex=j-1;

		// Computing the rightmost alignment whose distance from $source$ in the DAG is
		// maximum.
		longestDistance=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]>longestDistance) longestDistance=distances[vertex];
		}
		if (longestDistance+1<minPathLength) {  // No path from $source$ with the required structure
if (IO.SHOW_STD_ERR) IO.printErr("getPrefixSubstring> failure 1, minPathLength="+minPathLength+" longestDistance+1="+(longestDistance+1));
			Math.set(tmp,2,-1);
			return;
		}
if (IO.SHOW_STD_ERR) IO.printErr("longestDistance="+longestDistance);
		rightmostEnd=-1; rightmostAlignment=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]<longestDistance) continue;
			end=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if (end>rightmostEnd) {
				rightmostAlignment=vertex;
				rightmostEnd=end;
			}
		}
if (IO.SHOW_STD_ERR) IO.printErr("rightmostEnd="+rightmostEnd);
		length=rightmostEnd-Alignments.alignments[ReadA.sortedAlignments[source].id][4];
		
		// Computing the end of the rightmost B-right-maximal alignment
		rightmostRightMaximalEnd=-1;
		for (i=s; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]==-1) continue;
			end=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if (end>rightmostEnd) continue;
			if (ReadA.sortedAlignments[vertex].isRightMaximalB==1 && !ReadA.sortedAlignments[vertex].lowQualityEnd && end>rightmostRightMaximalEnd) rightmostRightMaximalEnd=end;
		}
		
		// No path from $source$ with the required structure
		if (length<prefixSuffixThreshold+1) {
			if (length>identityThreshold) {
				// Disabling all maximal alignment ends in $[j..rightmostRightMaximalEnd-
				// identityThreshold]$, where $j$ is the last position of $source$: such
				// events are likely to create a second peak close to
				// $rightmostRightMaximalEnd$.
				j=Alignments.alignments[ReadA.sortedAlignments[source].id][4];		
				if (rightmostRightMaximalEnd-identityThreshold>j) {
if (IO.SHOW_STD_ERR) IO.printErr("DISABLING MAXIMAL EVENTS IN ["+j+".."+(rightmostRightMaximalEnd-identityThreshold)+"] (prefix)");
					for (i=source; i<nVertices; i++) {
						if (Alignments.alignments[ReadA.sortedAlignments[i].id][3]>rightmostRightMaximalEnd) break;
						start=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
						end=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
						if ( ReadA.sortedAlignments[i].isRightMaximal==1 &&
						     Math.abs(sourceStartA,start)<=identityThreshold &&
						     end>=j && end<rightmostRightMaximalEnd-identityThreshold
						   ) ReadA.sortedAlignments[i].endAdded=true;
					}
				}
			}
			Math.set(tmp,2,-1);
			return;
		}

		// Checking whether a prefix end is significantly different from the end of
		// $source$
		if (length>2*identityThreshold) {  // 2 is used, since $source$ has the leftmost end in the path.
			found=false;
			sourceEndA=Alignments.alignments[ReadA.sortedAlignments[source].id][4];
			vertex=predecessors[rightmostAlignment];
			while (vertex!=source) {
				end=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
				if (end>=sourceEndA+2*identityThreshold && ReadA.sortedAlignments[vertex].isRightMaximalB==1) {
					found=true;
					break;
				}
				vertex=predecessors[vertex];
			}
			if (!found) {
				Math.set(tmp,2,-1);
				return;
			}
		}

if (IO.SHOW_STD_ERR) IO.printErr("rightmostRightMaximalEnd="+rightmostRightMaximalEnd);
		tmp[0]=rightmostAlignment;
		tmp[1]=longestDistance;
		tmp[2]=rightmostRightMaximalEnd;
	}


	/**
	 * Searches for a maximal dense substring, of suffix type, that starts from alignment
	 * $ReadA.sortedAlignments[source]$. This is done by computing the longest path, of
	 * length at least $minPathLength$, between $ReadA.sortedAlignments[source]$
	 * and every other vertex in the suffix DAG, such that all alignments in the path are
	 * both B-left- and B-right-maximal, and such that the distance between the first
	 * position of an alignment with longest distance, and the first position of $source$
	 * is greater than $identityThreshold$.
	 *
	 * Remark: the procedure assumes that function $setMaximality$ has already been called
	 * on all elements of $ReadA.sortedAlignments$.
	 *
	 * @return $tmp[0]$: the position in $ReadA.sortedAlignments[0..nVertices-1]$ of an
	 * interval that maximizes the distance from $source$ in the suffix DAG, and that ends
	 * the farthest from $source$, or -1 if no valid path from $source$ is found;
	 * $tmp[1]$: the length of a longest path from $source$ in the suffix DAG;
	 * $tmp[2]$: the leftmost start of an alignment that is reachable from $source$ in the 
 	 * DAG, that starts at least at the start of $source$, and that is B-left-maximal.
	 */
	private static final void getSuffixSubstring(int source, int nVertices, int minPathLength) {
		boolean found;
		int i, j, k, s;
		int state, sourceStartA, sourceEndA, sourceComponent, vertex, rightmostAlignment, distance, longestDistance;
		int start, end, rightmostStart, rightmostEnd, length, leftmostLeftMaximalStart;
		int lastSortedVertex;

		// Computing the longest distance from $source$ to every other vertex
		Math.set(distances,nVertices-1,-1);
		Math.set(predecessors,nVertices-1,-1);
		sourceEndA=Alignments.alignments[ReadA.sortedAlignments[source].id][4];
		sourceComponent=components[source];
		distances[source]=0;
		s=original2sorted[source];
		lastSortedVertex=s;
		for (j=s+1; j<nVertices; j++) {
			vertex=sorted2original[j];
			if (components[vertex]!=sourceComponent) {
				lastSortedVertex=j-1;
				break;
			}
			if (ReadA.sortedAlignments[vertex].isLeftMaximalB!=1 || ReadA.sortedAlignments[vertex].isRightMaximalB!=1) continue;
			if (Math.abs(Alignments.alignments[ReadA.sortedAlignments[vertex].id][4],sourceEndA)>identityThreshold) continue;
			for (i=0; i<nInNeighbors[vertex]; i++) {
				k=inNeighbors[vertex][i];
				if (ReadA.sortedAlignments[k].isLeftMaximalB!=1 || ReadA.sortedAlignments[k].isRightMaximalB!=1) continue;
				if (distances[k]!=-1) distance=distances[k];
				else continue;
				distance++;
				if (distance>distances[vertex]) {
					distances[vertex]=distance;
					predecessors[vertex]=k;
				}
			}
		}
		if (j==nVertices) lastSortedVertex=j-1;

		// Computing the rightmost alignment whose distance from $source$ in the DAG is
		// maximum.
		longestDistance=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]>longestDistance) longestDistance=distances[vertex];
		}
		if (longestDistance+1<minPathLength) {  // No path from $source$ with the required structure
if (IO.SHOW_STD_ERR) IO.printErr("getSuffixSubstring> failure 1, minPathLength="+minPathLength+" longestDistance+1="+(longestDistance+1));
			Math.set(tmp,2,-1);
			return;
		}
if (IO.SHOW_STD_ERR) IO.printErr("getSuffixSubstring> longestDistance="+longestDistance);
		rightmostEnd=-1; rightmostStart=-1; rightmostAlignment=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]<longestDistance) continue;
if (IO.SHOW_STD_ERR) IO.printErr("alignment ["+Alignments.alignments[ReadA.sortedAlignments[vertex].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[vertex].id][4]+"] has the longest distance");
			end=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if (end>rightmostEnd) {
				rightmostAlignment=vertex;
				rightmostEnd=end;
			}
			start=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			if (start>rightmostStart) rightmostStart=start;
		}
		length=rightmostStart-Alignments.alignments[ReadA.sortedAlignments[source].id][3];
		sourceStartA=Alignments.alignments[ReadA.sortedAlignments[source].id][3];

		// Computing the start of the leftmost B-left-maximal alignment
		leftmostLeftMaximalStart=Math.POSITIVE_INFINITY;
		for (i=s; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]==-1) continue;
			start=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			if (start>rightmostEnd) continue;
			if (ReadA.sortedAlignments[vertex].isLeftMaximalB==1 && !ReadA.sortedAlignments[vertex].lowQualityStart && start<leftmostLeftMaximalStart) leftmostLeftMaximalStart=start;
		}

		// No path from $source$ with the required structure
		if (length<prefixSuffixThreshold) {
			if (length>identityThreshold && rightmostStart>=leftmostLeftMaximalStart+identityThreshold) {
				// Disabling all maximal alignment starts in $[leftmostLeftMaximalStart+
				// identityThreshold..rightmostStart]$: such events are likely to create a
				// second peak close to $leftmostLeftMaximalStart$.
if (IO.SHOW_STD_ERR) IO.printErr("DISABLING MAXIMAL EVENTS IN ["+(leftmostLeftMaximalStart+identityThreshold)+".."+rightmostStart+"] (suffix)");
				for (i=source+1; i<nVertices; i++) {
					start=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					end=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					if (start>rightmostStart) break;
					if ( ReadA.sortedAlignments[i].isLeftMaximal==1 &&
						 Math.abs(sourceEndA,end)<=identityThreshold &&
					     start>=leftmostLeftMaximalStart+identityThreshold && start<=rightmostStart
					   ) ReadA.sortedAlignments[i].startAdded=true;
				}
			}
if (IO.SHOW_STD_ERR) IO.printErr("getSuffixSubstring> failure 2. rightmostStart="+rightmostStart+", rightmostEnd="+rightmostEnd+", first position of source = "+Alignments.alignments[ReadA.sortedAlignments[source].id][3]);
			Math.set(tmp,2,-1);
			return;
		}

		// Checking whether a suffix start is significantly different from
		// $rightmostStart$
		if (length>2*identityThreshold) {  // 2 is used, since $rightmostStart$ is, indeed, the rightmost.
			found=false;
			vertex=predecessors[rightmostAlignment];
			while (vertex!=source) {
				start=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
				if (start<=rightmostStart-2*identityThreshold && ReadA.sortedAlignments[vertex].isLeftMaximalB==1) {
					found=true;
					break;
				}
				vertex=predecessors[vertex];
			}
			if (!found) {
if (IO.SHOW_STD_ERR) IO.printErr("getSuffixSubstring> failure 3");
				Math.set(tmp,2,-1);
				return;
			}
		}
		tmp[0]=rightmostAlignment;
		tmp[1]=longestDistance;
		tmp[2]=leftmostLeftMaximalStart;
	}


	/**
	 * Assume that $ReadA.sortedAlignments[0..nVertices-1]$ contains all the alignments of
	 * readA, with any readB, in any orientation, sorted by $startA$. Such order induces
	 * a DAG whose vertices are the elements of $ReadA.sortedAlignments$, and in which
	 * there is an arc $(v,w)$ iff $ReadA.sortedAlignments[w]$ starts in readA after the
	 * start of $ReadA.sortedAlignments[v]$, before the end of $ReadA.sortedAlignments[v]$,
	 * and not farther than $maxDistance$ from the start of $ReadA.sortedAlignments[v]$.
	 * Moreover, $ReadA.sortedAlignments[w]$ must end at most $maxDistanceEnd$ positions
	 * after the end of $ReadA.sortedAlignments[v]$.
	 * This procedure searches for a maximal dense substring of $readA$ that starts from
	 * alignment $ReadA.sortedAlignments[source]$. This is done by computing the longest
	 * path, of length at least $minPathLength$, between $ReadA.sortedAlignments
	 * [source]$ and every other vertex in the DAG, such that all alignments are both
	 * B-left- and B-right-maximal.
	 *
	 * This captures the case of a repeat that replicates in the genome by taking a
	 * substring of itself, or by having more than one internal deletion. The procedure
	 * does not try to distinguish between such cases. Having multiple internal deletions
	 * is very common for transposable elements: see e.g. \cite{quesneville2003detection}.
	 * An instance of a repeat with internal deletions might even look like a sequence of 
	 * adjacent dense substrings of substring type.
	 *
	 * Remark: let $v$ be the vertex associated with $source$, and let $w$ be a
	 * destination vertex. The procedure discards substring $W=[v.start..w.end]$ if $w$'s
	 * start or end positions are too close to those of $v$, since this would be a simple
	 * repeat rather than a dense substring of substring type.
	 *
	 * Remark: the procedure assumes that function $setMaximality$ has already been called
	 * on all elements of $ReadA.sortedAlignments$.
	 *
	 * @return $tmp[0]$: the position in $ReadA.sortedAlignments$ of an alignment that
	 * maximizes the distance from $ReadA.sortedAlignments[source]$ in the DAG, and that
	 * ends the farthest from $ReadA.sortedAlignments[source]$, or -1 if no valid path
	 * from $ReadA.sortedAlignments[source]$ is found; $tmp[1]$: the length of a longest
	 * path from $source$ in the DAG.
	 */
	private static final void getSubstring(int source, int nVertices, int minPathLength) {
		int i, j, k, s, vertex;
		int endA, startA1, endA1, startA2, endA2;
		int rightmostAlignment, distance, longestDistance, rightmostEndA;
		int alignmentLength, leftmostLeftMaximalStartA, rightmostRightMaximalEndA;
		int sourceComponent, lastSortedVertex;

		// Computing the longest distance from $ReadA.sortedAlignments[source]$ to every
		// other alignment.
		Math.set(distances,ReadA.lastSortedAlignment,-1);
		Math.set(predecessors,ReadA.lastSortedAlignment,-1);
		Math.set(sumOfLengths,ReadA.lastSortedAlignment,0);
		sourceComponent=components[source];
		distances[source]=0;
		sumOfLengths[source]=Alignments.alignments[ReadA.sortedAlignments[source].id][4]-Alignments.alignments[ReadA.sortedAlignments[source].id][3]+1;
		s=original2sorted[source];
 		lastSortedVertex=s;
 		longestDistance=0;
		for (j=s+1; j<nVertices; j++) {
// if (IO.SHOW_STD_ERR) IO.printErr("getSubstring> considering alignment ["+Alignments.alignments[ReadA.sortedAlignments[j].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[j].id][4]+"]");
			vertex=sorted2original[j];
			if (components[vertex]!=sourceComponent) {
				lastSortedVertex=j-1;
				break;
			}
			if (ReadA.sortedAlignments[vertex].isLeftMaximalB!=1 || ReadA.sortedAlignments[vertex].isRightMaximalB!=1) continue;
			alignmentLength=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4]-Alignments.alignments[ReadA.sortedAlignments[vertex].id][3]+1;
			for (i=0; i<nInNeighbors[vertex]; i++) {
				k=inNeighbors[vertex][i];
				if (distances[k]==-1) continue;
				if (ReadA.sortedAlignments[k].isLeftMaximalB!=1 || ReadA.sortedAlignments[k].isRightMaximalB!=1) continue;
				distance=distances[k]+1;
				if (distance>distances[vertex]) {
					distances[vertex]=distance;
					predecessors[vertex]=k;
					sumOfLengths[vertex]=sumOfLengths[k]+alignmentLength;
				}
			}
			if (distances[vertex]>longestDistance) longestDistance=distances[vertex];
if (IO.SHOW_STD_ERR) IO.printErr("getSubstring> distance="+distances[j]);
		}
		if (j==nVertices) lastSortedVertex=j-1;
		if (longestDistance+1<minPathLength) {
if (IO.SHOW_STD_ERR) IO.printErr("failure 1: longestDistance+1="+(longestDistance+1)+" minPathLength="+minPathLength);
			Math.set(tmp,1,-1);
			return;
		}

		// Computing the alignment with rightmost end, such that: (1) its starting
		// position in $readA$ and the starting position of
		// $ReadA.sortedAlignments[source]$ are distinct;
		// (2) its ending position in $readA$ and the ending position of
		// $ReadA.sortedAlignments[source]$ are distinct;
		// (3) its distance from $ReadA.sortedAlignments[source]$ in the DAG is maximum.
		startA1=Alignments.alignments[ReadA.sortedAlignments[source].id][3];
		endA1=Alignments.alignments[ReadA.sortedAlignments[source].id][4];
		longestDistance=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]==-1) continue;
			if ( Math.abs(Alignments.alignments[ReadA.sortedAlignments[vertex].id][3],startA1)<=repeatThreshold ||
			     Math.abs(Alignments.alignments[ReadA.sortedAlignments[vertex].id][4],endA1)<=repeatThreshold ) continue;
			if (distances[vertex]>longestDistance) longestDistance=distances[vertex];
		}
if (IO.SHOW_STD_ERR) IO.printErr("getSubstring> longestDistance="+longestDistance);
		if (longestDistance+1<minPathLength) {  // No path from $source$ with the required structure
if (IO.SHOW_STD_ERR) IO.printErr("failure 2, minPathLength="+minPathLength);
			Math.set(tmp,1,-1);
			return;
		}
		tmp[1]=longestDistance;
		rightmostEndA=-1; rightmostAlignment=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]==-1) continue;
			endA=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if ( distances[vertex]!=longestDistance ||
				 Math.abs(Alignments.alignments[ReadA.sortedAlignments[vertex].id][3],startA1)<=repeatThreshold ||
			     Math.abs(endA,endA1)<=repeatThreshold ) continue;
			if (endA>rightmostEndA) {
				rightmostAlignment=vertex;
				rightmostEndA=endA;
			}
		}
		tmp[0]=rightmostAlignment;

		// Computing the start of the leftmost B-left-maximal alignment
		leftmostLeftMaximalStartA=-1;
		for (i=s; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			startA1=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			if (startA1>rightmostEndA) break;
			if (ReadA.sortedAlignments[vertex].isLeftMaximalB==1 && !ReadA.sortedAlignments[vertex].lowQualityStart) {
				leftmostLeftMaximalStartA=startA1;
				break;
			}
		}
		tmp[2]=leftmostLeftMaximalStartA;

		// Computing the last position of the rightmost B-right-maximal alignment
		rightmostRightMaximalEndA=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			startA1=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			endA1=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if (ReadA.sortedAlignments[vertex].isRightMaximalB==1 && !ReadA.sortedAlignments[vertex].lowQualityEnd && endA1<=rightmostEndA && endA1>rightmostRightMaximalEndA) rightmostRightMaximalEndA=endA1;
		}
		tmp[3]=rightmostRightMaximalEndA;
	}


	/**
	 * Marks all alignments whose interval in readA has a significant overlap with a dense
	 * substring that replicates by taking substrings of itself.
	 *
	 * Remark: the procedure assumes that $substrings$ contains all and only the dense
	 * substrings that replicate by taking substrings of themselves, sorted by startA.
	 * The procedure also assumes $ReadA.sortedAlignments$ to be sorted by startA,
	 * and it assumes all alignments implied by periodic substrings to be located in
	 * $ReadA.sortedAlignments[lastAlignment+1..]$.
	 *
	 * Remark: the procedure does not use marked alignments to update the $startA$,
	 * $endA$, $sum*A$ and $n*A$ values of the corresponding dense substring. This is
	 * because, inside a dense substring of substring type, starting/ending positions of
	 * alignments are assumed to be uniformly distributed. Thus, averaging maximal events
	 * that are close to the beginning/end of the dense substring would just move the
	 * boundaries towards the interior.
	 */
	private static final void markAlignments(int lastAlignment) {
		int i, j;
		int alignmentStart, alignmentEnd;
		int firstJForNextSubstring;
		int length;

		// Marking all alignments that have a significant intersection with a dense
		// substring
		for (i=0; i<=lastSubstring; i++) substrings[i].minPrefixLength=Math.POSITIVE_INFINITY;
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			ReadA.sortedAlignments[i].inDenseSubstring=null;
			ReadA.sortedAlignments[i].startInDenseSubstring=false;
			ReadA.sortedAlignments[i].endInDenseSubstring=false;
		}
		i=0; j=0; firstJForNextSubstring=-1;
		if (lastSubstring==-1) return;
		while (i<=lastSubstring) {
			if (j>lastAlignment) {
				i++;
				if (firstJForNextSubstring!=-1) j=firstJForNextSubstring;
				firstJForNextSubstring=-1;
				continue;
			}
			if (ReadA.sortedAlignments[j].inPeriodicSubstring) {
				j++;
				continue;
			}
			alignmentStart=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
			if (alignmentStart>substrings[i].endA) {
				i++;
				if (firstJForNextSubstring!=-1) j=firstJForNextSubstring;
				firstJForNextSubstring=-1;
				continue;
			}
			alignmentEnd=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
			if (alignmentEnd<substrings[i].startA) {
				j++;
				continue;
			}
			if (i<lastSubstring && firstJForNextSubstring==-1 && alignmentEnd>=substrings[i+1].startA) firstJForNextSubstring=j;
			if (alignmentStart>=substrings[i].startA && alignmentStart<=substrings[i].endA) ReadA.sortedAlignments[j].startInDenseSubstring=true;
			if (alignmentEnd>=substrings[i].startA && alignmentEnd<=substrings[i].endA) ReadA.sortedAlignments[j].endInDenseSubstring=true;
			if (alignmentStart<=substrings[i].startA && alignmentEnd>=substrings[i].endA) {  // Discarding alignments that contain the dense substring
				j++;
				continue;
			}
			if (Intervals.isApproximatelyContained(alignmentStart,alignmentEnd,substrings[i].startA,substrings[i].endA)) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markAlignments> alignment ["+alignmentStart+".."+alignmentEnd+"] is approximately contained in substring ["+substrings[i].startA+".."+substrings[i].endA+"]");
				ReadA.sortedAlignments[j].inDenseSubstring=substrings[i];
				if (ReadA.sortedAlignments[j].isLeftMaximalB==1 && ReadA.sortedAlignments[j].isRightMaximalB==1) {
					length=ReadA.sortedAlignments[j].getALength();
					if (length<substrings[i].minPrefixLength) substrings[i].minPrefixLength=length;
				}
			}
			j++;
		}
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].minPrefixLength==Math.POSITIVE_INFINITY) substrings[i].minPrefixLength=-1;
		}
	}


	/**
	 * Sets $impliedBy{Prefix,Suffix}Substring$ for all alignments in the longest path
	 * from $destination$ to $source$, using $predecessors$. Such alignments are used just
	 * for deciding the boundaries and replication type of a dense substring, and will be
	 * discarded in the following steps of the pipeline. Then, the procedure calls
	 * $markImpliedAlignmentsAggressive$.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments$ to be sorted by startA, and
	 * all alignments implied by periodic substrings to be located in
	 * $ReadA.sortedAlignments[lastAlignment+1..]$.
	 *
	 * @return the number of alignments marked by the procedure.
	 */
	private static final int markImpliedAlignmentsPath(DenseSubstring substring, int lastAlignment, int minPathLength) {
		int vertex, out;
		
		// Marking alignments in the longest path from source to destination
		out=0;
		vertex=substring.rightAlignment;  // $rightAlignment$ is always the destination of the path
		while (vertex!=-1) {
			if (substring.prefixReplication && ReadA.sortedAlignments[vertex].impliedByPrefixSubstring==null) {
				substring.nImpliedAlignmentsPrefix++;
				ReadA.sortedAlignments[vertex].impliedByPrefixSubstring=substring;
				out++;
			}
			else if (substring.suffixReplication && ReadA.sortedAlignments[vertex].impliedBySuffixSubstring==null) {
				substring.nImpliedAlignmentsSuffix++;
				ReadA.sortedAlignments[vertex].impliedBySuffixSubstring=substring;
				out++;
			}
			vertex=predecessors[vertex];
		}
		if (out<minPathLength) return out;

		// Marking other alignments
		out+=markImpliedAlignmentsAggressive(substring,lastAlignment);
		return out;
	}


	/**
	 * Marks as implied all the alignments that behave like a prefix (respectively, a
	 * suffix) of a dense substring that replicates by taking prefixes (respectively,
	 * suffixes) of itself. An implied alignment can start (end) far to the right (left)
	 * from the start (end) of $substring$ if it is not B-left- (B-right-) maximal.
	 *
	 * The procedure tries to estimate peaks in the $\delta$ function computed on prefix/
	 * suffix alignment endpoints. Prefix/suffix alignments whose endpoints belongs to 
	 * such peaks are not marked as implied, since they are more likely to be induced by 
	 * repeats than by the prefix/suffix substring. Peaks are detected on function 
	 * $\delta$ rather than on the number of maximal events, since such number if affected
	 * by the frequency of repeats.
	 *
	 * Remark: $minPrefixLength$ is not likely to be a local maximum of function $\delta$,
	 * assuming uniform probability of choosing a prefix end inside a prefix substring. 
	 * Actually it is likely to have a lower value of function $\delta$ than the other 
	 * positions inside the substring, since the $minPrefixLength$ position is either 
	 * contained inside a prefix or it is chosen as an end, whereas other positions $x$
	 * (excluding the last one) are either contained, or end, or not contained (i.e.
	 * prefixes can end before $x$).
	 *
	 * Remark: the procedure uses implied alignments to update the boundaries, the
	 * $sum*A$ and $n*A$ quantities, and $min{Prefix,Suffix}Length$, of $substring$.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments$ to be sorted by startA, and
	 * all alignments implied by periodic substrings to be located in
	 * $ReadA.sortedAlignments[lastAlignment+1..]$.
	 *
	 * @return the number of alignments marked by the procedure.
	 */
	private static final int markImpliedAlignmentsAggressive(DenseSubstring substring, int lastAlignment) {
		final int LENGTH_THRESHOLD = identityThreshold<<1;  // Arbitrary
		final int MIN_COUNT = 3;  // Min. n. of events to be used in a $\delta$ histogram
		final int MIN_BANDWIDTH = (identityThreshold<<1)/3;
		final int MAX_BANDWIDTH = identityThreshold;
		int i, j;
		int id, pos, posPrime, length, nLengths, maxPrefixLength, maxSuffixLength, from, to;
		final int substringLength, first, last, prefixOrigin, suffixOrigin;
		int newStartA, newSumStartA, newNStartA, newEndA, newSumEndA, newNEndA;
		int bandwidth, lastDeltaPointPrefix, lastDeltaPointSuffix;
		int minPrefixLength_boundary, minSuffixLength_boundary, lastPoint;
		int out;

		// Updating $minPrefixLength$ and $minSuffixLength$. This is necessary,
		// since an alignment might not be in the longest path but still have the shortest
		// prefix (suffix) length (e.g. in a prefix substring, an alignment that does not
		// start the leftmost, but that ends the leftmost).
		substringLength=substring.endA-substring.startA+1;
		first=substring.getFromAlignment(identityThreshold);
		last=substring.getToAlignment(identityThreshold,lastAlignment);
		maxPrefixLength=-1; maxSuffixLength=-1;
		for (i=first; i<=last; i++) {	
			if (ReadA.sortedAlignments[i].isLeftMaximalB!=1 || ReadA.sortedAlignments[i].isRightMaximalB!=1) continue;
			if (ReadA.sortedAlignments[i].inPeriodicSubstring) continue;
			id=ReadA.sortedAlignments[i].id;
			length=ReadA.sortedAlignments[i].getALength();
			if ( substring.prefixReplication &&
				 Math.abs(Alignments.alignments[id][3],substring.sumStartA/substring.nStartA)<=identityThreshold &&
				 Alignments.alignments[id][4]<=substring.endA+identityThreshold ) {
				if (length<substring.minPrefixLength) substring.minPrefixLength=length;
				if (length>maxPrefixLength) maxPrefixLength=length;
			}
			if ( substring.suffixReplication &&
				 Math.abs(Alignments.alignments[id][4],substring.sumEndA/substring.nEndA)<=identityThreshold &&
				 Alignments.alignments[id][3]>=substring.startA-identityThreshold ) {
				if (length<substring.minSuffixLength) substring.minSuffixLength=length;
				if (length>maxSuffixLength) maxSuffixLength=length;
			}
		}
		maxPrefixLength+=IO.quantum; maxSuffixLength+=IO.quantum;
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("markImpliedAlignmentsAggressive> substring ["+substring.startA+".."+substring.endA+"]");
	IO.printErr("minPrefixLength="+substring.minPrefixLength+", minSuffixLength="+substring.minSuffixLength);
}
		
		// Filling $deltaPointsPrefix$ and $deltaPointsSuffix$
		pos=substring.sumStartA/substring.nStartA;
		posPrime=substring.sumEndA/substring.nEndA;
		prefixOrigin=substring.minPrefixLength-LENGTH_THRESHOLD; 
		if (substring.prefixReplication) Math.set(countsPrefix,maxPrefixLength-prefixOrigin,0);
		suffixOrigin=substring.minSuffixLength-LENGTH_THRESHOLD;		
		if (substring.suffixReplication) Math.set(countsSuffix,maxSuffixLength-suffixOrigin,0);				
		lastDeltaPointPrefix=-1; lastDeltaPointSuffix=-1;
		substring.hasLengthBoundariesPrefix=true; substring.hasLengthBoundariesSuffix=true;
		for (i=first; i<=last; i++) {
			if (ReadA.sortedAlignments[i].inDenseSubstring!=null) continue;
			if (ReadA.sortedAlignments[i].inPeriodicSubstring) continue;
			id=ReadA.sortedAlignments[i].id;
			if (substring.prefixReplication && isAlignmentImpliedByPrefix(i,substring,false)) {
				length=Alignments.alignments[id][4]-pos;
				if (length<substringLength-LENGTH_THRESHOLD) {
					lastDeltaPointPrefix++;
					deltaPointsPrefix[lastDeltaPointPrefix].position=length-prefixOrigin;
					deltaPointsPrefix[lastDeltaPointPrefix].mass=1;
					for (j=0; j<=length-prefixOrigin; j++) countsPrefix[j]++;
				}
			}
			if (substring.suffixReplication && isAlignmentImpliedBySuffix(i,substring,false)) {
				length=posPrime-Alignments.alignments[id][3];
				if (length<substringLength-LENGTH_THRESHOLD) {
					lastDeltaPointSuffix++;
					deltaPointsSuffix[lastDeltaPointSuffix].position=length-suffixOrigin;
					deltaPointsSuffix[lastDeltaPointSuffix].mass=1;
					for (j=0; j<=length-suffixOrigin; j++) countsSuffix[j]++;
				}
			}
		}
		
		// Computing $delta{Prefix,Suffix}$ and finding peaks
		substring.period=-1;
		lastLengthBoundaryPrefix=-1;
		if (lastDeltaPointPrefix+1>=MIN_ALIGNMENTS_FOR_LENGTH_BIAS) {
			lastDeltaPointPrefix=Points.sortAndCompact(deltaPointsPrefix,lastDeltaPointPrefix);
			bandwidth=Points.estimateDensity(deltaPointsPrefix,lastDeltaPointPrefix,deltaPrefix,MIN_BANDWIDTH,MAX_BANDWIDTH);
			for (i=0; i<=lastDeltaPointPrefix; i++) {
				from=(int)deltaPointsPrefix[i].position+1;
				to=Math.min((int)deltaPointsPrefix[i].position+bandwidth,countsPrefix.length-1);
				for (j=from; j<=to; j++) countsPrefix[j]+=deltaPointsPrefix[i].getMass();
			}
			for (i=0; i<=maxPrefixLength-prefixOrigin; i++) {
				if (countsPrefix[i]>=MIN_COUNT) deltaPrefix[i]/=countsPrefix[i];
				else deltaPrefix[i]=0;
			}
			estimateFromDelta(deltaPrefix,maxPrefixLength,substring,true,first,last,prefixOrigin,!inPeriodicInterval(substring,tmpPInterval),LENGTH_THRESHOLD);
		}
		substring.hasLengthBoundariesPrefix=lastLengthBoundaryPrefix>=0;
		lastLengthBoundarySuffix=-1;
		if (lastDeltaPointSuffix+1>=MIN_ALIGNMENTS_FOR_LENGTH_BIAS) {
			lastDeltaPointSuffix=Points.sortAndCompact(deltaPointsSuffix,lastDeltaPointSuffix);	
			bandwidth=Points.estimateDensity(deltaPointsSuffix,lastDeltaPointSuffix,deltaSuffix,MIN_BANDWIDTH,MAX_BANDWIDTH);	
			for (i=0; i<=lastDeltaPointSuffix; i++) {
				from=(int)deltaPointsSuffix[i].position+1;
				to=Math.min((int)deltaPointsSuffix[i].position+bandwidth,countsSuffix.length-1);
				for (j=from; j<=to; j++) countsSuffix[j]+=deltaPointsSuffix[i].getMass();
			}
			for (i=0; i<=maxSuffixLength-suffixOrigin; i++) {
				if (countsSuffix[i]>=MIN_COUNT) deltaSuffix[i]/=countsSuffix[i];
				else deltaSuffix[i]=0;
			}
			estimateFromDelta(deltaSuffix,maxSuffixLength,substring,false,first,last,suffixOrigin,!inPeriodicInterval(substring,tmpPInterval),LENGTH_THRESHOLD);
		}
		substring.hasLengthBoundariesSuffix=lastLengthBoundarySuffix>=0;
		

if (lastLengthBoundaryPrefix>=0) {
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("lengthBoundariesPrefix of dense substring id="+substring.id+"=["+substring.startA+".."+substring.endA+"] lastLengthBoundaryPrefix="+lastLengthBoundaryPrefix+" prefixOrigin="+prefixOrigin+":");
	if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastLengthBoundaryPrefix; i++) IO.printErr((/*prefixOrigin+*/lengthBoundariesPrefix[i])+","); }
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("");
}
if (lastLengthBoundarySuffix>=0) {
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("lengthBoundariesSuffix of dense substring id="+substring.id+"=["+substring.startA+".."+substring.endA+"] lastLengthBoundarySuffix="+lastLengthBoundarySuffix+"  suffixOrigin="+suffixOrigin+":");
	if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastLengthBoundarySuffix; i++) IO.printErr((/*suffixOrigin+*/lengthBoundariesSuffix[i])+","); }
	if (IO.SHOW_STD_ERR_PRIME) IO.printErr("");
}


		// Marking
		out=0;
		newStartA=substring.startA;
		newEndA=substring.endA;
		newSumStartA=substring.sumStartA;
		newNStartA=substring.nStartA;
		newSumEndA=substring.sumEndA;
		newNEndA=substring.nEndA;
		for (i=first; i<=last; i++) {
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("markImpliedAlignmentsAggressive> substring "+(substring.prefixReplication?"prefix":"-")+"+"+(substring.suffixReplication?"suffix":"-")+" ["+substring.startA+".."+substring.endA+"] considering alignment id="+ReadA.sortedAlignments[i].id+"=["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"]...");
	IO.printErr("markImpliedAlignmentsAggressive> substring sumStartA/nStartA="+(substring.sumStartA/substring.nStartA)+" sumEndA/nEndA="+(substring.sumEndA/substring.nEndA));
}


			if (ReadA.sortedAlignments[i].inDenseSubstring!=null) continue;
			if (ReadA.sortedAlignments[i].inPeriodicSubstring) continue;
			id=ReadA.sortedAlignments[i].id;
			if (substring.prefixReplication) {
				pos=substring.sumStartA/substring.nStartA;
				length=(Alignments.alignments[id][4]-pos)-prefixOrigin;
				j=lastLengthBoundaryPrefix>=0?Arrays.binarySearch(lengthBoundariesPrefix,0,lastLengthBoundaryPrefix+1,length):-1;
				if ( lastLengthBoundaryPrefix>=0 && 
				     ( j>=0 || 
				       (-j-1)%3!=0 || 
				       ( -j-1<=lastLengthBoundaryPrefix && 
						 ( Math.abs(length,lengthBoundariesPrefix[-j])<=identityThreshold || 
				           (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=identityThreshold)
						 )
					   )
					 )
				   ) {
					if (ReadA.sortedAlignments[i].impliedByPrefixSubstring==substring) {
			   			if (ReadA.sortedAlignments[i].inDenseSubstring==null) ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=ReadA.sortedAlignments[i].impliedByPrefixSubstring;
						ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
			   		}
			   		else if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring==null && 
					          ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime==null && 
							  ReadA.sortedAlignments[i].inDenseSubstring==null &&
					          isAlignmentImpliedByPrefix(i,substring,false)
					        ) {
						ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=substring;
					}
				}
				else if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring==null && 
				          ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime==null && 
						  ReadA.sortedAlignments[i].inDenseSubstring==null &&
				          isAlignmentImpliedByPrefix(i,substring,false)
			            ) {
					ReadA.sortedAlignments[i].impliedByPrefixSubstring=substring;
					if (Math.abs(Alignments.alignments[id][3],pos)<=identityThreshold && ReadA.sortedAlignments[i].isLeftMaximalB==1) {
						substring.nImpliedAlignmentsPrefix++;
						newSumStartA+=Alignments.alignments[id][3];
						newNStartA++;
					}
					if (Alignments.alignments[id][3]<newStartA) newStartA=Alignments.alignments[id][3];
					if (ReadA.sortedAlignments[i].isRightMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityEnd && Alignments.alignments[id][4]>substring.maximalEndA) substring.maximalEndA=Alignments.alignments[id][4];
					out++;
				}
			}
			if (substring.suffixReplication) {				
				pos=substring.sumEndA/substring.nEndA;
				length=(pos-Alignments.alignments[id][3])-suffixOrigin;
				j=lastLengthBoundarySuffix>=0?Arrays.binarySearch(lengthBoundariesSuffix,0,lastLengthBoundarySuffix+1,length):-1;
				if ( lastLengthBoundarySuffix>=0 && 
				     ( j>=0 || 
				       (-j-1)%3!=0 || 
					   ( -j-1<=lastLengthBoundarySuffix &&
				         ( Math.abs(length,lengthBoundariesSuffix[-j])<=identityThreshold || 
				           (-j-1==0?false:Math.abs(length,lengthBoundariesSuffix[-j-3])<=identityThreshold)
						 )
					   )
					 )
				   ) {
					if (ReadA.sortedAlignments[i].impliedBySuffixSubstring==substring) {
						if (ReadA.sortedAlignments[i].inDenseSubstring==null) ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=ReadA.sortedAlignments[i].impliedBySuffixSubstring;
			   			ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
			   		}
			   		else if ( ReadA.sortedAlignments[i].impliedBySuffixSubstring==null && 
					          ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime==null && 
							  ReadA.sortedAlignments[i].inDenseSubstring==null &&
					          isAlignmentImpliedBySuffix(i,substring,false)
					        ) {
						ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=substring;
					}
				}
				else if ( ReadA.sortedAlignments[i].impliedBySuffixSubstring==null && 
				          ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime==null && 
						  ReadA.sortedAlignments[i].inDenseSubstring==null &&
				          isAlignmentImpliedBySuffix(i,substring,false)
				        ) {
					ReadA.sortedAlignments[i].impliedBySuffixSubstring=substring;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("MARKED!");
					if (Math.abs(Alignments.alignments[id][4],pos)<=identityThreshold && ReadA.sortedAlignments[i].isRightMaximalB==1) {
						substring.nImpliedAlignmentsSuffix++;
						newSumEndA+=Alignments.alignments[id][4];
						newNEndA++;
					}
					if (Alignments.alignments[id][4]>newEndA) newEndA=Alignments.alignments[id][4];
					if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityStart && Alignments.alignments[id][3]>=substring.startA && Alignments.alignments[id][3]<substring.maximalStartA) substring.maximalStartA=Alignments.alignments[id][3];
					out++;
				}
			}
		}
		substring.startA=newStartA;
		substring.endA=newEndA;
		substring.sumStartA=newSumStartA;
		substring.nStartA=newNStartA;
		substring.sumEndA=newSumEndA;
		substring.nEndA=newNEndA;

		// Updating $minPrefixLength$ and $minSuffixLength$
		minPrefixLength_boundary=-1;
		if (substring.prefixReplication && lastLengthBoundaryPrefix>=0) {
			length=substring.minPrefixLength-prefixOrigin;
			j=Arrays.binarySearch(lengthBoundariesPrefix,0,lastLengthBoundaryPrefix+1,length);
			if ( j>=0 || 
			     (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundaryPrefix &&
			       ( Math.abs(length,lengthBoundariesPrefix[-j])<=identityThreshold || 
			         (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=identityThreshold)
				   )
				 )
		       ) {
				minPrefixLength_boundary=j;
				substring.minPrefixLength=Math.POSITIVE_INFINITY;
				for (i=first; i<=last; i++) {
					if (ReadA.sortedAlignments[i].impliedByPrefixSubstring!=substring || ReadA.sortedAlignments[i].isLeftMaximalB!=1 || ReadA.sortedAlignments[i].isRightMaximalB!=1) continue;
					if (ReadA.sortedAlignments[i].inPeriodicSubstring) continue;
					id=ReadA.sortedAlignments[i].id;
					if (Math.abs(Alignments.alignments[id][3],substring.startA)>identityThreshold) continue;
					length=ReadA.sortedAlignments[i].getALength();
				 	if (length<substring.minPrefixLength) substring.minPrefixLength=length;
				}
				if (substring.minPrefixLength==Math.POSITIVE_INFINITY) substring.minPrefixLength=-1;
			}
		}
		minSuffixLength_boundary=-1;
		if (substring.suffixReplication && lastLengthBoundarySuffix>=0) {
			length=substring.minSuffixLength-suffixOrigin;
			j=Arrays.binarySearch(lengthBoundariesSuffix,0,lastLengthBoundarySuffix+1,length);
			if ( j>=0 || 
			     (-j-1)%3!=0 || 
			     ( -j-1<=lastLengthBoundarySuffix &&
				   ( Math.abs(length,lengthBoundariesPrefix[-j])<=identityThreshold || 
			         (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=identityThreshold)
				   )
				 )
			   ) {
				minSuffixLength_boundary=j;
				substring.minSuffixLength=Math.POSITIVE_INFINITY;
				for (i=first; i<=last; i++) {
					if (ReadA.sortedAlignments[i].impliedBySuffixSubstring!=substring || ReadA.sortedAlignments[i].isLeftMaximalB!=1 || ReadA.sortedAlignments[i].isRightMaximalB!=1) continue;
					if (ReadA.sortedAlignments[i].inPeriodicSubstring) continue;
					id=ReadA.sortedAlignments[i].id;
					if (Math.abs(Alignments.alignments[id][4],substring.endA)>identityThreshold) continue;
					length=ReadA.sortedAlignments[i].getALength();
				 	if (length<substring.minSuffixLength) substring.minSuffixLength=length;
				}
				if (substring.minSuffixLength==Math.POSITIVE_INFINITY) substring.minSuffixLength=-1;
			}
		}

		// Updating $lengthBiasIntervals$
		if (lastLengthBoundaryPrefix>=0) {
			pos=substring.sumStartA/substring.nStartA;
			for (i=minPrefixLength_boundary==-1?0:(minPrefixLength_boundary/3+1)*3; i<lastLengthBoundaryPrefix; i+=3) {
				lastLengthBiasInterval++;
				lengthBiasIntervals[lastLengthBiasInterval].firstPosition=pos+prefixOrigin+lengthBoundariesPrefix[i];
				lengthBiasIntervals[lastLengthBiasInterval].center=pos+prefixOrigin+lengthBoundariesPrefix[i+1];
				lengthBiasIntervals[lastLengthBiasInterval].lastPosition=pos+prefixOrigin+lengthBoundariesPrefix[i+2];
				lengthBiasIntervals[lastLengthBiasInterval].readB=ReadA.id;
				lengthBiasIntervals[lastLengthBiasInterval].orientation=true;
			}
		}
		if (lastLengthBoundarySuffix>=0) {
			pos=substring.sumEndA/substring.nEndA;
			for (i=minSuffixLength_boundary==-1?0:(minSuffixLength_boundary/3+1)*3; i<lastLengthBoundarySuffix; i+=3) {
				lastLengthBiasInterval++;
				lengthBiasIntervals[lastLengthBiasInterval].firstPosition=pos-suffixOrigin-lengthBoundariesSuffix[i+2];
				lengthBiasIntervals[lastLengthBiasInterval].center=pos-suffixOrigin-lengthBoundariesSuffix[i+1];
				lengthBiasIntervals[lastLengthBiasInterval].lastPosition=pos-suffixOrigin-lengthBoundariesSuffix[i];
				lengthBiasIntervals[lastLengthBiasInterval].readB=ReadA.id;
				lengthBiasIntervals[lastLengthBiasInterval].orientation=true;
			}
		}
		
		// Setting $period$ if it is not already set.
		if (substring.period==-1) {
			lastPoint=getShiftsPrefixSuffix(substring,first,last,true);
			if (lastPoint>=0 && Points.mass(lengthPoints,0,lastPoint)>=MIN_ALIGNMENTS_FOR_LENGTH_BIAS) {
				substring.period=PeriodicSubstrings.estimatePeriod(lengthPoints,lastPoint,PeriodicSubstrings.MIN_PERIOD,PeriodicSubstrings.minIntervalLength,PeriodicSubstrings.minLocalMaxDistance);
				substring.strongPeriodSignal=substring.period>0&&PeriodicSubstrings.periodTmp[4]>1;
			}
			else {
				substring.period=0;
				substring.strongPeriodSignal=false;
			}
		}
		
		return out;
	}
	
	
	/**
	 * Unmarks all alignments that were implied by a non-weak prefix/suffix substring.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments[0..lastAlignment]$ to be
	 * sorted by $startA$.
	 */
	private static final void unmarkImpliedAlignments(DenseSubstring substring, int lastAlignment) {
		int i;
		final int first = substring.getFromAlignment(identityThreshold);
		final int last = substring.getToAlignment(identityThreshold,lastAlignment);
		for (i=first; i<=last; i++) {
			if (ReadA.sortedAlignments[i].impliedByPrefixSubstring==substring) ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
			if (ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime==substring) ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;
			if (ReadA.sortedAlignments[i].impliedBySuffixSubstring==substring) ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
			if (ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime==substring) ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
		}
	}
	
	
	/**
	 * Symmetrical to $markImpliedAlignmentsAggressive_weak_prefix()$.
	 */
	private static final void unmarkImpliedAlignments_weak_prefix(DenseSubstring substring, int lastAlignment) {
		final int THRESHOLD = Alignments.minAlignmentLength;  // Arbitrary
		final int last = substring.getFromAlignmentPrime(identityThreshold,lastAlignment);
		final int pos = substring.startA;
		int i;
		
		if (substring.prefixReplication) {
			for (i=last; i>=0; i--) {
				if (ReadA.sortedAlignments[i].endA()<pos-THRESHOLD) break;
				if (ReadA.sortedAlignments[i].impliedByPrefixSubstring==substring) ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
				if (ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime==substring) ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;
			}
		}
		else if (substring.substringReplication) {
			for (i=last; i>=0; i--) {					
				if (ReadA.sortedAlignments[i].endA()<pos-THRESHOLD) break;
				if (ReadA.sortedAlignments[i].inDenseSubstring==substring) ReadA.sortedAlignments[i].inDenseSubstring=null;
			}
			for (i=last+1; i<=lastAlignment; i++) {
				if (ReadA.sortedAlignments[i].inDenseSubstring==substring) ReadA.sortedAlignments[i].inDenseSubstring=null;
			}
		}
	}
	
	
	/**
	 * Symmetrical to $markImpliedAlignmentsAggressive_weak_suffix()$.
	 */
	private static final void unmarkImpliedAlignments_weak_suffix(DenseSubstring substring, int lastAlignment) {
		final int THRESHOLD = Alignments.minAlignmentLength;  // Arbitrary
		final int first = substring.getFromAlignment(identityThreshold);
		final int pos = substring.endA;
		int i;
		
		if (substring.suffixReplication) {
			for (i=first; i<=lastAlignment; i++) {
				if (ReadA.sortedAlignments[i].startA()>pos+THRESHOLD) break;
				if (ReadA.sortedAlignments[i].impliedBySuffixSubstring==substring) ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
				if (ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime==substring) ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
			}
		}
		else if (substring.substringReplication) {
			for (i=first; i<=lastAlignment; i++) {
				if (ReadA.sortedAlignments[i].startA()>pos+THRESHOLD) break;
				if (ReadA.sortedAlignments[i].inDenseSubstring==substring) ReadA.sortedAlignments[i].inDenseSubstring=null;
			}
			for (i=first-1; i>=0; i--) {
				if (ReadA.sortedAlignments[i].inDenseSubstring==substring) ReadA.sortedAlignments[i].inDenseSubstring=null;
			}			
		}
	}
	
	
	/**
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ to be sorted by 
	 * $firstPosition$.
	 *
	 * @param tmpInterval temporary space.
	 */
	private static final boolean inPeriodicInterval(DenseSubstring substring, PeriodicSubstringInterval tmpInterval) {
		boolean previousSortByID;
		int i, j;
		
		tmpInterval.firstPosition=substring.startA;
		previousSortByID=PeriodicSubstringInterval.sortByID;
		PeriodicSubstringInterval.sortByID=false;
		i=Arrays.binarySearch(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1,tmpInterval);
		PeriodicSubstringInterval.sortByID=previousSortByID;
		if (i<0) i=-i-1;
		for (j=i-1; j>=0; j--) {
			if (PeriodicSubstrings.intervals[j].lastPosition<=substring.startA) continue;
			if ( Intervals.isApproximatelyContained(substring.startA,substring.endA,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition) ||
				 Intervals.areApproximatelyIdentical(substring.startA,substring.endA,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition)
			   ) return true;
		}
		for (j=i; j<=PeriodicSubstrings.lastInterval; j++) {
			if (PeriodicSubstrings.intervals[j].firstPosition>=substring.endA) break;
			if ( Intervals.isApproximatelyContained(substring.startA,substring.endA,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition) ||
				 Intervals.areApproximatelyIdentical(substring.startA,substring.endA,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition)
			   ) return true;
		}
		return false;
	}
	
	
	/**
	 * @param alignment position of the alignment in $ReadA.sortedAlignments$.
	 */
	private static final boolean isAlignmentImpliedByPrefix(int alignment, DenseSubstring substring, boolean strict) {
		final int id = ReadA.sortedAlignments[alignment].id;
		final int pos = substring.sumStartA/substring.nStartA;
		
		if (Alignments.alignments[id][4]>substring.endA+identityThreshold) return false;
		if ( Math.abs(Alignments.alignments[id][3],pos)<=identityThreshold &&
			 ( !substring.hasLengthBoundariesPrefix ||
			   ( ReadA.sortedAlignments[alignment].isRightMaximalB==1 &&
				 Alignments.alignments[id][4]>=pos+substring.minPrefixLength-identityThreshold
			   )
			 )
		   ) return true;
		if ( !strict &&
			 Alignments.alignments[id][3]>pos+identityThreshold &&
			 ReadA.sortedAlignments[alignment].isLeftMaximalB!=1 &&
			 ReadA.sortedAlignments[alignment].isRightMaximalB==1 &&
			 Alignments.alignments[id][4]>=pos+substring.minPrefixLength-identityThreshold
		   ) return true;
		return false;
	}
	
	
	/**
	 * @param alignment position of the alignment in $ReadA.sortedAlignments$.
	 */
	private static final boolean isAlignmentImpliedBySuffix(int alignment, DenseSubstring substring, boolean strict) {	
		final int id = ReadA.sortedAlignments[alignment].id;
		final int pos = substring.sumEndA/substring.nEndA;
		
		if (Alignments.alignments[id][3]<substring.startA-identityThreshold) return false;
		if ( Math.abs(Alignments.alignments[id][4],pos)<=identityThreshold &&
			 ( !substring.hasLengthBoundariesSuffix ||
			   ( ReadA.sortedAlignments[alignment].isLeftMaximalB==1 &&
				 Alignments.alignments[id][3]<=pos-substring.minSuffixLength+identityThreshold
			   )
			 )
		   ) return true;
		if ( !strict &&
			 Alignments.alignments[id][4]<pos-identityThreshold &&
			 ReadA.sortedAlignments[alignment].isRightMaximalB!=1 &&
			 ReadA.sortedAlignments[alignment].isLeftMaximalB==1 &&
			 Alignments.alignments[id][3]<=pos-substring.minSuffixLength+identityThreshold
		   ) return true;
		return false;
	}


	/**
	 * Fits a regression tree on a $\delta$ histogram of lengths, detects high local 
	 * maxima using function $Leaves.setHighLocalMax$, and stores in array 
	 * $lengthBoundaries{Prefix,Suffix}$ the pairs of boundaries of all runs of local-
	 * maximum leaves that contain a high local maximum, in length order.
	 *
	 * Remark: $lengthBoundaries{Prefix,Suffix}$ is relative to $delta$.
	 *
	 * @param removePeriodBoundaries if TRUE, the procedure tries to detect whether the 
	 * prefixes (respectively, suffixes) of the dense substring are taken at periodic 
	 * distances. If this happens, only boundaries that are not too close to a multiple of
	 * the period are reported. This is because, if e.g. prefixes are taken at periodic 
	 * offsets, it is very likely for local maxima in prefix length to arise at such 
	 * periods.
	 * @param offset amount subtracted from every element in $lengthBoundaries$ when 
	 * comparing it to the period.
	 */
	private static final void estimateFromDelta(double[] delta, int maxLength, DenseSubstring substring, boolean prefix, int firstAlignment, int lastAlignment, int origin, boolean removePeriodBoundaries, int offset) {
		final double CONSTANT_THRESHOLD = 0.1;
		final double MIN_HIGH = 2.5;
		final double ZERO_RATIO = 1000.0;  // Arbitrary
		final int MAX_WIDTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		boolean justFirstLocalMax, isHigh;
		int i, j;
		int nLocalMaximumLeaves, firstLocalMaximum, lastLeaf, last, lastPrime;
		double center, mass;
		int[] vector;
		
		// Building regression tree and marking high local maxima
		if (prefix) lastLengthBoundaryPrefix=-1;
		else lastLengthBoundarySuffix=-1;
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("deltaLength: dense substring: "+substring);
	IO.printErr("maxLength="+maxLength+" origin="+origin);
	for (i=0; i<=maxLength-origin; i++) IO.printErr(delta[i]); 
}
	
		final int from = Histograms.firstNonzero(delta,0,maxLength-origin);
		final int to = Histograms.lastNonzero(delta,0,maxLength-origin);
		if (from==-1 || Histograms.isConstant(delta,from,to,CONSTANT_THRESHOLD)) {
			// Computing prefix/suffix bias statistics to be shown in output
			removeBoundariesAtPeriods(null,-1,substring,prefix,firstAlignment,lastAlignment,0);
			return;
		}
		final double MAX_STD = (Histograms.max(delta,from,to)-Histograms.min(delta,from,to))/20;
		nLocalMaximumLeaves=RegressionTree.buildRegressionTree(delta,from,to,MAX_STD,1,0,true,Alignments.minAlignmentLength);
		if (nLocalMaximumLeaves==0) {
			// Computing prefix/suffix bias statistics to be shown in output
			removeBoundariesAtPeriods(null,-1,substring,prefix,firstAlignment,lastAlignment,0);
			return;
		}
		justFirstLocalMax=Histograms.oneMaxAndConstant(delta,from,to,CONSTANT_THRESHOLD,MAX_WIDTH)||Histograms.allMaxAndConstant(delta,from,to,CONSTANT_THRESHOLD,MAX_WIDTH);
		if (!justFirstLocalMax) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Regression tree on deltaLength (before markRuns): MAX_STD="+MAX_STD);
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=RegressionTree.lastLeaf; i++) IO.printErr(RegressionTree.leaves[i]); }
			Leaves.setHighLocalMax(RegressionTree.leaves,RegressionTree.lastLeaf,to,true,true,null,0,false,MAX_WIDTH);
			// Enforcing all high-enough local maxima to be marked as high, even if
			// $setHighLocalMax()$ didn't mark them.
			for (i=0; i<=RegressionTree.lastLeaf; i++) {
				if (RegressionTree.leaves[i].isLocalMaximum && RegressionTree.leaves[i].lastPoint-RegressionTree.leaves[i].firstPoint+1<=MAX_WIDTH && Leaves.getLocalMaxRatio(i,RegressionTree.leaves,RegressionTree.lastLeaf,to,true,true,ZERO_RATIO)>=MIN_HIGH) RegressionTree.leaves[i].isHigh=true;
			}
		}
		RegressionTree.markRunsOfLocalMaximumLeaves();

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Regression tree on deltaLength (after setHighLocalMax and markRuns): lastLocalMaximum="+RegressionTree.lastLocalMaximum);
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=RegressionTree.lastLeaf; i++) IO.printErr(RegressionTree.leaves[i]); }


		// Building array $lengthBoundaries{Prefix,Suffix}$
		firstLocalMaximum=0; isHigh=false;
		for (i=0; i<=RegressionTree.lastLocalMaximum; i++) {
			if (RegressionTree.leaves[i].isHigh) isHigh=true;
			if (!RegressionTree.leaves[i].marked) continue;
			if (!isHigh) {
				firstLocalMaximum=i+1;
				isHigh=false;
				continue;
			}
			center=0.0; mass=0.0;
			for (j=RegressionTree.leaves[firstLocalMaximum].firstPoint; j<=RegressionTree.leaves[i].lastPoint; j++) {
				center+=j*delta[j];
				mass+=delta[j];
			}
			center/=mass;
			if (prefix) {
				lengthBoundariesPrefix[++lastLengthBoundaryPrefix]=RegressionTree.leaves[firstLocalMaximum].firstPoint;
				lengthBoundariesPrefix[++lastLengthBoundaryPrefix]=Math.round(center);
				lengthBoundariesPrefix[++lastLengthBoundaryPrefix]=RegressionTree.leaves[i].lastPoint;
			}
			else {
				lengthBoundariesSuffix[++lastLengthBoundarySuffix]=RegressionTree.leaves[firstLocalMaximum].firstPoint;
				lengthBoundariesSuffix[++lastLengthBoundarySuffix]=Math.round(center);
				lengthBoundariesSuffix[++lastLengthBoundarySuffix]=RegressionTree.leaves[i].lastPoint;
			}
			firstLocalMaximum=i+1;
			isHigh=false;
		}
		if (prefix) {
			vector=lengthBoundariesPrefix;
			last=lastLengthBoundaryPrefix;
		}
		else {
			vector=lengthBoundariesSuffix;
			last=lastLengthBoundarySuffix;
		}
		if (last==-1) {
			// Computing prefix/suffix bias statistics to be shown in output
			removeBoundariesAtPeriods(null,-1,substring,prefix,firstAlignment,lastAlignment,0);
			return;
		}

		// Filtering $lengthBoundaries{Prefix,Suffix}$
		lastPrime=removeBoundariesAtPeriods(removePeriodBoundaries?vector:null,removePeriodBoundaries?last:-1,substring,prefix,firstAlignment,lastAlignment,offset);
		if (removePeriodBoundaries) last=lastPrime;
		
		// Not needed any more: 
		// last=removeBoundariesWithLowDensity(vector,last,delta,maxLength,origin);

		// Enforcing all positions in $lengthBoundaries{Prefix,Suffix}$ to be distinct.
		for (i=1; i<=last; i++) {
			if (vector[i]<=vector[i-1]) vector[i]=vector[i-1]+1;
		}

		if (prefix) lastLengthBoundaryPrefix=last;
		else lastLengthBoundarySuffix=last;
	}


	/**
	 * Marks an alignment that starts (ends) inside the suffix (prefix) region of a dense
	 * substring that replicates by taking suffixes (prefixes) of itself, and that ends
	 * (starts) far to the right (left) of the end (start) of the dense substring.
	 *
	 * Remark: the procedure assumes $substrings$ to be sorted by startA, and
	 * $ReadA.sortedAlignments[0..lastAlignment]$ to contain all the alignments that are
	 * not implied by a periodic substring, sorted by startA.
	 */
	public static final void markStartEnd(int lastAlignment) {
		int i, j;
		int firstAlignmentForNextSubstring;

		// Initializing
		for (i=0; i<=lastAlignment; i++) {
			ReadA.sortedAlignments[i].startInDenseSuffix=false;
			ReadA.sortedAlignments[i].endInDensePrefix=false;
		}

		// Marking
		i=0; j=0; firstAlignmentForNextSubstring=-1;
		while (i<=lastSubstring) {

if (IO.SHOW_STD_ERR) IO.printErr("markStartEnd> comparing alignment ("+j+"|"+ReadA.sortedAlignments[j].id+")["+Alignments.alignments[ReadA.sortedAlignments[j].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[j].id][4]+"] to dense substring "+substrings[i]);

			if (j>lastAlignment || Alignments.alignments[ReadA.sortedAlignments[j].id][3]>substrings[i].endA) {
				i++;
				if (firstAlignmentForNextSubstring!=-1) {
					j=firstAlignmentForNextSubstring;
					firstAlignmentForNextSubstring=-1;
				}
				continue;
			}
			if (Alignments.alignments[ReadA.sortedAlignments[j].id][4]<substrings[i].startA) {
				j++;
				continue;
			}
			if (i<lastSubstring && firstAlignmentForNextSubstring==-1 && Alignments.alignments[ReadA.sortedAlignments[j].id][4]>=substrings[i+1].startA) firstAlignmentForNextSubstring=j;
			if (substrings[i].substringReplication) {
				j++;
				continue;
			}
			if ( substrings[i].prefixReplication &&
			     Alignments.alignments[ReadA.sortedAlignments[j].id][4]>=substrings[i].sumStartA/substrings[i].nStartA+substrings[i].minPrefixLength-identityThreshold &&
				 Alignments.alignments[ReadA.sortedAlignments[j].id][3]<substrings[i].startA-identityThreshold ) ReadA.sortedAlignments[j].endInDensePrefix=true;
			if ( substrings[i].suffixReplication &&
			     Alignments.alignments[ReadA.sortedAlignments[j].id][3]<=substrings[i].sumEndA/substrings[i].nEndA-substrings[i].minSuffixLength+identityThreshold &&
				 Alignments.alignments[ReadA.sortedAlignments[j].id][4]>substrings[i].endA+identityThreshold ) ReadA.sortedAlignments[j].startInDenseSuffix=true;
			j++;
		}
	}


	/**
	 * Merges substrings with high Jaccard similarity, setting the replication modes of
	 * the resulting substring to the union of the replication modes of all merged
	 * substrings. Then, concatenates compatible prefix/suffix substrings (see procedure 
	 * $getConnectedComponent_concatenation()$).
	 *
	 * Remark: non-weak dense substrings that replicate by taking substrings of themselves
	 * can be merged with other non-weak dense substrings only if the latter slightly 
	 * contain them.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments[0..lastAlignment]$ to be
	 * sorted by startA, endA, and it assumes all alignments implied by periodic
	 * substrings to be located in $ReadA.sortedAlignments[lastAlignment+1..]$.
	 * The procedure assumes $substrings$ to be sorted by startA.
	 *
	 * @return the number of alignments in $ReadA.sortedAlignments$ that are implied by a
	 * dense substring.
	 */
	private static final int mergeSubstrings(int lastAlignment) {
		boolean atLeastOneConcatenation;
		int i, j;
		int previousLastSubstring;
		DenseSubstring tmpSubstring;
		Alignment tmpAlignment;



for (i=0; i<=lastAlignment; i++) {
	tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
	if (tmpSubstring==null) continue;
	boolean found = false;
	for (j=0; j<=lastSubstring; j++) {
		if (tmpSubstring==substrings[j]) {
			found=true;
			break;
		}
	}
	if (!found) {
		System.err.println("mergeSubstrings> 0  alignment has unknown impliedByDenseSubstring: ");
		System.err.println(ReadA.sortedAlignments[i].toString());
		System.err.println(ReadA.sortedAlignments[i].toStringPointers());
	}
}



		// Merging substrings with identical intervals
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].isMerged=false;
			substrings[i].representative=substrings[i];
		}
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].isMerged) continue;
			getConnectedComponent_merge(i);
		}
		j=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].representative!=substrings[i]) continue;
			j++;
			tmpSubstring=substrings[j];
			substrings[j]=substrings[i];
			substrings[i]=tmpSubstring;
		}
		if (j<lastSubstring) {
			lastSubstring=j;			
			for (i=0; i<=lastAlignment; i++) {  // Updating alignment pointers
				if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null) {
					tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring.representative;
					if (tmpSubstring.substringReplication) {
						ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
						ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring;				
					}
					else ReadA.sortedAlignments[i].impliedByDenseSubstring=tmpSubstring;
				}
				else if (ReadA.sortedAlignments[i].inDenseSubstring!=null) {
					ReadA.sortedAlignments[i].inDenseSubstring=ReadA.sortedAlignments[i].inDenseSubstring.representative;
				}
			}
		}
		

if (IO.SHOW_INTERACTIVE) {
	IO.printErr("DenseSubstrings> after mergeIdenticalSubstrings():");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+" =====> representative="+substrings[i].representative);
	IO.printErr("DenseSubstrings> alignments:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		int x = ReadA.sortedAlignments[i].id;
		IO.printErr("["+Alignments.alignments[x][3]+".."+Alignments.alignments[x][4]+"], "+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+", "+(ReadA.sortedAlignments[i].inDenseSubstring==null?"":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
	}
	IO.printErr("DenseSubstrings> alignments implied by dense:");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		int x = ReadA.sortedAlignments[i].id;
		if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null) IO.printErr("["+Alignments.alignments[x][3]+".."+Alignments.alignments[x][4]+"], "+ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode());
	}
}


		// Concatenating compatible prefix/suffix substrings
		countAlignmentsAtBoundaries(lastAlignment);
		for (i=0; i<=lastSubstring; i++) substrings[i].isConcatenated=false;
		atLeastOneConcatenation=false;
		for (i=0; i<lastSubstring; i++) {
			if (substrings[i].representative!=substrings[i] || (substrings[i].prefixReplication && substrings[i].suffixReplication)) continue;
			j=getConnectedComponent_concatenation(i);
			if (j>0) atLeastOneConcatenation=true;
		}
		
		
if (IO.SHOW_INTERACTIVE) {
	IO.printErr("DenseSubstrings> after getConnectedComponent_concatenation():");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]+" =====> representative="+substrings[i].representative);
}		
		
		if (atLeastOneConcatenation) {
			updateConcatenatedSumN(0,lastSubstring);			
			updateConcatenatedReplicationModes();
			j=-1;
			for (i=0; i<=lastSubstring; i++) {
				if (substrings[i].representative!=substrings[i]) continue;
				j++;
				if (j!=i) {
					tmpSubstring=substrings[j];
					substrings[j]=substrings[i];
					substrings[i]=tmpSubstring;
				}
			}
			previousLastSubstring=lastSubstring; lastSubstring=j;
			for (i=0; i<=lastAlignment; i++) {  // Updating alignment pointers
				if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null) {
					tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring.representative;
					if (tmpSubstring.substringReplication) {
						ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
						ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring;
					}
					else ReadA.sortedAlignments[i].impliedByDenseSubstring=tmpSubstring;
				}
				else if (ReadA.sortedAlignments[i].inDenseSubstring!=null) {
					ReadA.sortedAlignments[i].inDenseSubstring=ReadA.sortedAlignments[i].inDenseSubstring.representative;
				}
			}
			findPeaksAfterConcatenation(lastAlignment,previousLastSubstring,stack);
		}
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("DenseSubstrings> after closing implied pointers in alignments:");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]);
	IO.printErr("DenseSubstrings> alignments:");
	for (i=0; i<=lastAlignment; i++) {
		int x = ReadA.sortedAlignments[i].id;
		IO.printErr("["+Alignments.alignments[x][3]+".."+Alignments.alignments[x][4]+"], "+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
	}
}
		
		// Enforcing straddling/containment with dense substrings of substring type
		cleanDenseSubstringsOfSubstringType(lastAlignment);
		DenseSubstring.order=DenseSubstring.STARTA;
		if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("DenseSubstrings> after cleanDenseSubstringsOfSubstringType():");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]);
}
		
		// Concatenating weak dense substrings
		concatenateWeakSubstrings(lastAlignment,stack);
		DenseSubstring.order=DenseSubstring.STARTA;
		if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);



if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("DenseSubstrings> after concatenateWeakSubstrings():");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]);
	IO.printErr("DenseSubstrings> alignments:");
	for (i=0; i<=lastAlignment; i++) {
		int x = ReadA.sortedAlignments[i].id;
		IO.printErr("["+Alignments.alignments[x][3]+".."+Alignments.alignments[x][4]+"], "+(ReadA.sortedAlignments[i].impliedByDenseSubstring==null?"":ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode())+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
	}
}

		// Collecting at the beginning of $ReadA.sortedAlignments$ all alignments that are
		// implied by a dense substring, sorted by startA, endA.
		Alignment.order=Alignment.IMPLIEDBYDENSE_INDENSE_STARTA;
		if (lastAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,lastAlignment+1);
		j=0;
		for (i=0; i<=lastAlignment; i++) {
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) {
				j=i;
				break;
			}
		}
		return j;
	}
	
	
	/**
	 * Remark: if a dense substring is contained in a long-period periodic interval, it 
	 * is discarded iff it is longer than the period, or if it is weak and has few maximal
	 * alignments. This is because, in practice, long-period intervals contain many non-
	 * implied alignments: we aggressively keep only those with a strong signal. Non-
	 * implied alignments with maximal events inside a long-period interval may be due, 
	 * e.g., to fragments of other occurrences of the periodic substring, with a different 
	 * phase. We could alternatively keep a dense substring contained in a long-period 
	 * interval, iff it has enough maximal alignments that are close to its first/last 
	 * positions; we haven't tried this because it seems less general than what we use.
	 * 
	 * Remark: the procedure assumes $substrings$ to be sorted by startA, and 
	 * $PeriodicSubstrings.{intervals,longPeriodIntervals}$ to be sorted by 
	 * $firstPosition$.
	 *
	 * Remark: dense substrings are not deleted but just pushed out of visibility. The
	 * corresponding implied alignments still point to hidden objects.
	 *
	 * Remark: the procedure assumes that field $mergedToInterval$ of all periodic
	 * substrings has already been set.
	 *
	 * Remark: the procedure uses array $lengthPoints$ as temporary space.
	 *
	 * @return TRUE iff a dense substring of substring type was merged to a long-period 
	 * interval.
	 */
	public static final boolean discardPeriodicSubstrings(int nAlignmentsHigh) {
		final int N_CONSECUTIVE_POINTS = 2;
		final double QUANTILE = 0.75;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		final double LONG_PERIOD_TOLERANCE = 0.8;  // Arbitrary
		final int DENSITY_RATIO = 20;  // Arbitrary
		final int BOUNDARY_ALIGNMENTS_THRESHOLD = IO.coverage*50;  // Arbitrary
		boolean found, atLeastOneFound;
		boolean previousSortByID, out;
		int i, j, k;
		int startJ, endJ, length, minLength;
		int firstJForNextI, newLastSubstring, rangeFrom, rangeFirst, rangeLast;
		DenseSubstring tmpSubstring;
		PeriodicSubstring tmpPSubstring;
		PeriodicSubstringInterval minInterval, tmpInterval;
		final PeriodicSubstringInterval STRADDLING_INTERVAL = new PeriodicSubstringInterval();
		
		out=false;
		Factorize.computeMaximalAlignments(true,false,true,ReadA.lastSortedAlignment);
		for (i=0; i<=lastSubstring; i++) substrings[i].periodicSubstringInterval=null;
		if (PeriodicSubstrings.lastInterval==-1) return out;
		
		// Discarding dense substrings: contained.
		countAlignmentsAtBoundaries(nAlignmentsHigh-1);
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) PeriodicSubstrings.intervals[i].mergedToDense=false;
		for (i=0; i<=PeriodicSubstrings.lastLongPeriodInterval; i++) PeriodicSubstrings.longPeriodIntervals[i].mergedToDense=false;
		i=0; j=0; firstJForNextI=-1; k=-1; atLeastOneFound=false;
		found=false; minInterval=null; minLength=Math.POSITIVE_INFINITY;
		while (i<=lastSubstring) {
			if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=substrings[i].endA) {
				if (!found) {
					k++;
					if (k!=i) {
						tmpSubstring=substrings[k];
						substrings[k]=substrings[i];
						substrings[i]=tmpSubstring;
					}
				}
				else {
					substrings[i].periodicSubstringInterval=minInterval;
					if (substrings[i].substringReplication) {
						minInterval.mergedToDense=true;
						if (minInterval.hasLongPeriod) out=true;
					}
					atLeastOneFound=true;
					found=false; minInterval=null; minLength=Math.POSITIVE_INFINITY;
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (PeriodicSubstrings.intervals[j].lastPosition<=substrings[i].startA) {
				j++;
				continue;
			}
			if (i<lastSubstring && firstJForNextI==-1 && PeriodicSubstrings.intervals[j].lastPosition>=substrings[i+1].startA) firstJForNextI=j;
			startJ=PeriodicSubstrings.intervals[j].firstPosition;
			endJ=PeriodicSubstrings.intervals[j].lastPosition;	
			if ( Intervals.areApproximatelyIdentical_lowQuality(substrings[i].startA,substrings[i].endA,startJ,endJ,ReadA.id) ||
				 ( Intervals.isApproximatelyContained_lowQuality(substrings[i].startA,substrings[i].endA,startJ,endJ,ReadA.id) &&
				   substrings[i].previousSumStartA<BOUNDARY_ALIGNMENTS_THRESHOLD &&
				   ( !PeriodicSubstrings.intervals[j].hasLongPeriod || 
					 ( (PeriodicSubstrings.intervals[j].period>0 && substrings[i].length()>=PeriodicSubstrings.intervals[j].period*LONG_PERIOD_TOLERANCE) ||
					   (PeriodicSubstrings.intervals[j].period<=0 && !PeriodicSubstrings.intervals[j].equalsOnePeriod && substrings[i].length()>PeriodicSubstrings.intervals[j].longestAlignment+identityThreshold) ||
					   (substrings[i].isWeak && ((double)substrings[i].nMaximalAlignments)/substrings[i].length()<=DENSITY_RATIO*((double)PeriodicSubstrings.intervals[j].mass)/PeriodicSubstrings.intervals[j].length())
					 )
				   )
				   /* &&  // Commented because this level of caution is not needed in practice.
				   !( (substrings[i].startA>startJ+identityThreshold && !Reads.isLeftMaximal(substrings[i].startA,ReadA.id,true)) &&
				      (substrings[i].endA<endJ-identityThreshold && !Reads.isRightMaximal(substrings[i].endA,ReadA.id,true))
				   )*/
				 )
			   ) {
				length=endJ-startJ;
				if (length<minLength) {
					minLength=length;
					minInterval=PeriodicSubstrings.intervals[j];  // Might be long-period
					found=true;
				}
			}
			j++;
		}
		newLastSubstring=k;
		
		// Discarding dense substrings: straddling.
		rangeFrom=0;
		rangeFirst=PeriodicSubstrings.intervals[0].firstPosition;
		rangeLast=PeriodicSubstrings.intervals[0].lastPosition;
		j=0; firstJForNextI=-1;
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			if (PeriodicSubstrings.intervals[i].firstPosition<=rangeLast+identityThreshold) {
				rangeLast=Math.max(rangeLast,PeriodicSubstrings.intervals[i].lastPosition);
				continue;
			}
			if (firstJForNextI!=-1) j=firstJForNextI;
			discardPeriodicSubstrings_straddling(rangeFirst,rangeLast,rangeFrom,i-1,j,STRADDLING_INTERVAL,tmp);
			firstJForNextI=tmp[0];
			atLeastOneFound|=tmp[1]>0;
			rangeFrom=i;
			rangeFirst=PeriodicSubstrings.intervals[i].firstPosition;
			rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
		}
		if (firstJForNextI!=-1) j=firstJForNextI;
		discardPeriodicSubstrings_straddling(rangeFirst,rangeLast,rangeFrom,i-1,j,STRADDLING_INTERVAL,tmp);
		atLeastOneFound|=tmp[1]>0;
		// Compacting dense substrings
		j=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].periodicSubstringInterval!=null) continue;
			j++;
			if (j!=i) {
				tmpSubstring=substrings[j];
				substrings[j]=substrings[i];
				substrings[i]=tmpSubstring;
			}
		}
		newLastSubstring=j;
		
		// Resetting pointers from long-period clones to the real long-period intervals
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].periodicSubstringInterval==null || substrings[i].periodicSubstringInterval==STRADDLING_INTERVAL || !substrings[i].periodicSubstringInterval.hasLongPeriod) continue;
			previousSortByID=PeriodicSubstringInterval.sortByID;
			PeriodicSubstringInterval.sortByID=false;
			j=Arrays.binarySearch(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1,substrings[i].periodicSubstringInterval);
			PeriodicSubstringInterval.sortByID=previousSortByID;
			if (j<0) {
				System.err.println("discardPeriodicSubstrings> ERROR: long-period interval not found.");
				System.exit(1);
			}
			found=false;
			for (k=j; k>=0; k--) {
				if (PeriodicSubstrings.longPeriodIntervals[k].equals(substrings[i].periodicSubstringInterval)) {
					PeriodicSubstrings.longPeriodIntervals[k].mergedToDense=substrings[i].periodicSubstringInterval.mergedToDense;
					substrings[i].periodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[k];
					found=true;
					break;
				}
			}
			if (found) continue;
			for (k=j+1; k<=PeriodicSubstrings.lastLongPeriodInterval; k++) {
				if (PeriodicSubstrings.longPeriodIntervals[k].equals(substrings[i].periodicSubstringInterval)) {
					PeriodicSubstrings.longPeriodIntervals[k].mergedToDense=substrings[i].periodicSubstringInterval.mergedToDense;
					substrings[i].periodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[k];
					break;
				}
			}
		}
		lastSubstring=newLastSubstring;
		
		// Updating pointers from alignments
		if (!atLeastOneFound) return out;
		for (i=0; i<nAlignmentsHigh; i++) {	
			if ( ReadA.sortedAlignments[i].impliedByDenseSubstring!=null && 
			     ReadA.sortedAlignments[i].impliedByDenseSubstring.periodicSubstringInterval!=null 
			   ) {
				tmpInterval=ReadA.sortedAlignments[i].impliedByDenseSubstring.periodicSubstringInterval;
				if (tmpInterval==STRADDLING_INTERVAL) ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				else {
					ReadA.sortedAlignments[i].periodicSubstringInterval=tmpInterval;
					tmpPSubstring=PeriodicSubstrings.findPeriodicSubstring(ReadA.sortedAlignments[i].startA(),ReadA.sortedAlignments[i].endA(),tmpInterval);
					if (tmpPSubstring!=null) ReadA.sortedAlignments[i].impliedByPeriodicSubstring=tmpPSubstring;
					else {
						tmpPSubstring=PeriodicSubstrings.getArtificialSubstring(tmpInterval);
						if (tmpPSubstring==null) tmpPSubstring=PeriodicSubstrings.addArtificialSubstring(tmpInterval);
						ReadA.sortedAlignments[i].impliedByPeriodicSubstring=tmpPSubstring;
					}
				}
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
			}
			if ( ReadA.sortedAlignments[i].inDenseSubstring!=null && 
			     ReadA.sortedAlignments[i].inDenseSubstring.periodicSubstringInterval!=null 
		       ) {
				tmpInterval=ReadA.sortedAlignments[i].inDenseSubstring.periodicSubstringInterval;
				if (tmpInterval==STRADDLING_INTERVAL) ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				else {
					ReadA.sortedAlignments[i].periodicSubstringInterval=tmpInterval;
					tmpPSubstring=PeriodicSubstrings.findPeriodicSubstring(ReadA.sortedAlignments[i].startA(),ReadA.sortedAlignments[i].endA(),tmpInterval);
					if (tmpPSubstring!=null) ReadA.sortedAlignments[i].impliedByPeriodicSubstring=tmpPSubstring;
					else {
						tmpPSubstring=PeriodicSubstrings.getArtificialSubstring(tmpInterval);
						if (tmpPSubstring==null) tmpPSubstring=PeriodicSubstrings.addArtificialSubstring(tmpInterval);
						ReadA.sortedAlignments[i].impliedByPeriodicSubstring=tmpPSubstring;
					}
				}
				ReadA.sortedAlignments[i].inDenseSubstring=null;
			}
		}
		
		return out;
	}
	
	
	/**
	 * Let $[rangeFirst..rangeLast]$ be a maximal range of contained or straddling 
	 * periodic intervals, and let $[rangeFrom..rangeTo]$ be the corresponding compact 
	 * range in $PeriodicSubstrings.intervals$. The procedure discards every dense 
	 * substring such that: (1) it is identical to or contained in $[rangeFirst..
	 * rangeLast]$; (2) it straddles at least two intervals in $[rangeFrom..rangeTo]$, 
	 * such intervals are adjacent, and at least one of them has short period.
	 *
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ and $DenseSubstrings.
	 * substrings$ to be sorted by first position.
	 *
	 * @param fromSubstring the dense substring from which the search should start;
	 * @param artificialInterval periodic interval to be used to mark discarded dense
	 * substrings;
	 * @param out[0]: the value of $fromSubstring$ to be used for the next range, or -1 if 
	 * no such value can be computed; out[1]: the number of substrings that have been 
	 * discarded by the procedure (possibly zero).
	 */
	private static final void discardPeriodicSubstrings_straddling(int rangeFirst, int rangeLast, int rangeFrom, int rangeTo, int fromSubstring, PeriodicSubstringInterval artificialInterval, int[] out) {
		final int INTERNAL_THRESHOLD = (identityThreshold)<<1;
		boolean found;
		int j, k, h;
		int firstJForNextI, firstPosition, lastPosition, substringStartA, substringEndA;
		
		out[0]=-1; out[1]=0;
		j=fromSubstring; firstJForNextI=-1;
		while (j<=lastSubstring) {
			substringStartA=substrings[j].startA;
			substringEndA=substrings[j].endA;
			if (substringStartA>=rangeLast) break;
			if (substringEndA<=rangeFirst || substrings[j].periodicSubstringInterval!=null) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && substringEndA>rangeLast) firstJForNextI=j;
			if ( !Intervals.areApproximatelyIdentical_lowQuality(substringStartA,substringEndA,rangeFirst,rangeLast,ReadA.id) &&
				 !( Intervals.isApproximatelyContained_lowQuality(substringStartA,substringEndA,rangeFirst,rangeLast,ReadA.id) &&
  				    !( (substringStartA>rangeFirst+identityThreshold && !Reads.isLeftMaximal(substringStartA,ReadA.id,true)) &&
  				       (substringEndA<rangeLast-identityThreshold && !Reads.isRightMaximal(substringEndA,ReadA.id,true))
  				    ) 
				 )
			   ) {
				j++;
				continue;
		    }
	   		for (k=rangeFrom; k<=rangeTo; k++) {
				firstPosition=PeriodicSubstrings.intervals[k].firstPosition;
				if (firstPosition>=substringEndA) break;
				lastPosition=PeriodicSubstrings.intervals[k].lastPosition;
				if (lastPosition<=substringStartA) continue;
				if (PeriodicSubstrings.intervals[k].hasLongPeriod) continue;
				found=false;
				if (lastPosition>=substringStartA+INTERNAL_THRESHOLD && lastPosition<=substringEndA-INTERNAL_THRESHOLD) {
					for (h=k+1; h<=rangeTo; h++) {
						if (PeriodicSubstrings.intervals[h].firstPosition>substringEndA-INTERNAL_THRESHOLD) break;
						if (PeriodicSubstrings.intervals[h].firstPosition<lastPosition-identityThreshold) continue;
						if ( PeriodicSubstrings.intervals[h].firstPosition<=lastPosition+identityThreshold || 
							 Reads.isRandomInsertion(ReadA.id,lastPosition+1,PeriodicSubstrings.intervals[h].firstPosition-1,true)
						   ) {
						    found=true;
						    break;
						}
					}
				}
				if (!found) {
					if (firstPosition>=substringStartA+INTERNAL_THRESHOLD && firstPosition<=substringEndA-INTERNAL_THRESHOLD) {
						for (h=k-1; h>=rangeFrom; h--) {
							if ( Math.abs(PeriodicSubstrings.intervals[h].lastPosition,firstPosition)<=identityThreshold ||
								 ( PeriodicSubstrings.intervals[h].lastPosition<firstPosition-identityThreshold && 
								   PeriodicSubstrings.intervals[h].lastPosition>=substringStartA+INTERNAL_THRESHOLD &&
								   Reads.isRandomInsertion(ReadA.id,PeriodicSubstrings.intervals[h].lastPosition+1,firstPosition-1,true)
								 )
							   ) {
								found=true;
								break;
							}
						}
					}
				}
				if (found) {
					substrings[j].periodicSubstringInterval=artificialInterval;
					out[1]++;
					break;
				}
			}
			j++;
		}
		out[0]=firstJForNextI==-1?j:firstJForNextI;
	}
	
	
	/*	
	 * Checks whether every substring of prefix and suffix type is a repeat that 
	 * replicates by a single internal deletion (rather than by taking just a prefix or 
	 * just a suffix).
	 *
	 * Then, tries to merge every substring $S$ of exclusively suffix type with a 
	 * substring $P$ of exclusively prefix type that starts the rightmost to the left of 
	 * $S$ in readA ($P$ may or may not overlap $S$). This pair can correspond both to a 
	 * repeat that replicates by a single internal deletion, or by taking 
	 * prefixes/suffixes that never come close to its last/first position. The procedure
	 * merges the substrings only in the first case.
	 *
	 * Remark: when a repeat replicates with a single internal deletion, such deletion can
	 * be so long that the corresponding prefix substring and suffix substring do not
	 * overlap in readA. The procedure handles this case.
	 *
	 * Remark: if another dense or periodic substring would be contained in a candidate 
	 * single-deletion substring, such substring is not created. This is to avoid merging
	 * distinct repeats that happen to be truncated and adjacent to one another multiple
	 * times in the genome.
	 *
	 * Remark: the procedure does not use dense substrings that are both contained in 
	 * (possibly different) short-period periodic substrings, to build a single-deletion
	 * substring. If the two containing periodic substrings P and Q are induced by the 
	 * insertion of a repeat inside a single periodic substring P', and if P' occurs in 
	 * other reads without an insertion, then there are likely alignments that are not 
	 * marked as implied by P or Q, but that are still induced by them. Such alignments 
	 * might form prefix or suffix substrings in the read, whose pairs of alignments might 
	 * be adjacent in the other reads that contain P' without the insertion. If P and Q 
	 * have different periods, the pattern of alignments used to detect a single-deletion 
	 * substring might be induced just by the fact that P and Q are adjacent in other 
	 * reads, without a repeat between them.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments[0..lastAlignment]$ to contain
	 * all the alignments implied by dense substrings; it assumes $substrings$ to be 
	 * sorted by startA, and $PeriodicSubstrings.substrings$ to be sorted by $minStartA$.
	 */
	private static final void getSingleDeletions(int lastAlignment) {
		final int THRESHOLD = IO.quantum;
		boolean found, found1, found2;
		int i, j, k;
		DenseSubstring tmpSubstring;
		PeriodicSubstring tmpPeriodic = new PeriodicSubstring();
		
		// Deciding whether to proceed or not with the following, costly steps of the
		// pipeline.
		found1=false;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].prefixReplication && substrings[i].suffixReplication) {
				found1=true;
				break;
			}
		}
		found2=false;
		for (i=1; i<=lastSubstring; i++) {
			if (substrings[i].suffixReplication && !substrings[i].prefixReplication && !substrings[i].substringReplication) {
				for (j=i-1; j>=0; j--) {
					if (substrings[j].prefixReplication && !substrings[j].suffixReplication && !substrings[j].substringReplication) {
						found2=true;
						break;
					}
				}
			}
		}
		if (!found1 && !found2) return;

		// Ensuring the necessary order in $ReadA.sortedAlignments[0..lastAlignment]$
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			// Assigning an ID to each dense substring is necessary for the following
			// sort.
			DenseSubstrings.substrings[i].id=i;
		}
		Alignment.order=Alignment.IMPLIEDBYDENSE_READB_ORIENTATION_STARTB;
		if (lastAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,lastAlignment+1);
		Alignment.order=Alignment.UNSORTED;  // Since we sorted just the interval $[0..lastAlignment]$
		for (i=0; i<=lastSubstring; i++) substrings[i].firstAlignment=-1;
		if (ReadA.sortedAlignments[0].impliedByDenseSubstring!=null) {
			ReadA.sortedAlignments[0].impliedByDenseSubstring.firstAlignment=0;
			for (i=1; i<=lastAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=ReadA.sortedAlignments[i-1].impliedByDenseSubstring) ReadA.sortedAlignments[i].impliedByDenseSubstring.firstAlignment=i;
			}
		}

		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Read "+ReadA.id+" getSingleDeletions> substrings:   (lastAlignment="+lastAlignment+")");
	for (i=0; i<=lastSubstring; i++) IO.printErr(substrings[i]);		
}

		// Detecting candidate single-deletion substrings, by checking substrings that
		// use both prefix replication and suffix replication.
		if (found1) {
			for (i=0; i<=lastSubstring; i++) {
				if (!substrings[i].prefixReplication || !substrings[i].suffixReplication || inPeriodicSubstring(substrings[i].startA,substrings[i].endA,tmpPeriodic)) continue;
				isSingleDeletion(substrings[i],substrings[i],lastAlignment,tmpPeriodic);
				if (singleDeletionTmp[0]!=-1) {
					substrings[i].singleDeletionReplication=true;
					// $prefixReplication$ and $suffixReplication$ are kept to true, since
					// there likely are prefix/suffix alignments that cannot be combined
					// with a corresponding suffix/prefix alignment.
					substrings[i].nSingleDeletions=singleDeletionTmp[0]==0?-1:singleDeletionTmp[0];
					substrings[i].shortestSingleDeletion=singleDeletionTmp[1];
					substrings[i].longestSingleDeletion=singleDeletionTmp[2];
				}
			}
		}

		// Detecting candidate single-deletion substrings, by merging a strictly-prefix
		// and a strictly-suffix substring.
		if (found2) {
			for (i=0; i<=lastSubstring; i++) substrings[i].representative=substrings[i];
			k=0;
			for (i=1; i<=lastSubstring; i++) {
				found=false;
				if (substrings[i].suffixReplication && !substrings[i].prefixReplication && !substrings[i].substringReplication && !inPeriodicSubstring(substrings[i].startA,substrings[i].endA,tmpPeriodic)) {
					for (j=k; j>=0; j--) {
						if (substrings[j].endA>substrings[i].endA+THRESHOLD) continue;
						if (!substrings[j].prefixReplication || substrings[j].suffixReplication || substrings[j].substringReplication || substrings[j].singleDeletionReplication) break;
						if (inPeriodicSubstring(substrings[j].startA,substrings[j].endA,tmpPeriodic)) break;
						isSingleDeletion(substrings[j],substrings[i],lastAlignment,tmpPeriodic);
						if (singleDeletionTmp[0]==-1) break;
						if (containsPeriodicSubstring(substrings[j].startA,substrings[i].endA,tmpPeriodic)) break;
						found=true;
						substrings[j].prefixReplication=false;
						substrings[j].singleDeletionReplication=true;
						substrings[j].endA=substrings[i].endA;
						substrings[j].sumEndA=substrings[i].sumEndA;
						substrings[j].nEndA=substrings[i].nEndA;
						substrings[j].nSingleDeletions=singleDeletionTmp[0]==0?-1:singleDeletionTmp[0];
						substrings[j].shortestSingleDeletion=singleDeletionTmp[1];
						substrings[j].longestSingleDeletion=singleDeletionTmp[2];
						substrings[j].nImpliedAlignmentsSuffix=substrings[i].nImpliedAlignmentsSuffix;
						substrings[j].period=Math.min(substrings[i].period,substrings[j].period);
						substrings[j].strongPeriodSignal=substrings[j].period>0&&(substrings[j].strongPeriodSignal||substrings[i].strongPeriodSignal);
						substrings[j].hasLengthBoundariesSuffix=substrings[i].hasLengthBoundariesSuffix;
						substrings[i].representative=substrings[j];
						break;
					}
				}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Read "+ReadA.id+" getSingleDeletions> found="+found+" for dense substring "+i);				
				if (found) continue;
				k++;
				if (k!=i) {
					tmpSubstring=substrings[k];
					substrings[k]=substrings[i];
					substrings[i]=tmpSubstring;
				}
			}
			lastSubstring=k;
			
			// Updating pointers from alignments
			for (i=0; i<=lastAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) continue;
				ReadA.sortedAlignments[i].impliedByDenseSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring.representative;
			}
		}
	}


	/**
	 * Computes the connected component of relation $Intervals.areApproximatelyIdentical$
	 * to which $substrings[fromSubstring]$ belongs. The procedure also sets the 
	 * $representative$ field of all substrings to point to the result of the merge.
	 *
	 * Remark: if a connected component contains a non-weak substring, the representative
	 * of the component is marked as non-weak.
	 *
	 * Remark: the procedure assumes the representative field of every substring to be
	 * initialized to the substring itself.
	 *
	 * Remark: the period of the representative substring is set to the minimum of the 
	 * periods of the substrings in its component, without any check for compatibility.
	 */
	private static final void getConnectedComponent_merge(int fromSubstring) {
		boolean newStrongPeriodSignal;
		int i, j;
		int top, newStartA, newEndA, newSumStartA, newNStartA, newSumEndA, newNEndA, newPeriod;
		double lengthI, lengthJ;

		top=0; stack[top]=fromSubstring;
		newStartA=substrings[fromSubstring].startA;
		newEndA=substrings[fromSubstring].endA;
		newSumStartA=substrings[fromSubstring].sumStartA;
		newNStartA=substrings[fromSubstring].nStartA;
		newSumEndA=substrings[fromSubstring].sumEndA;
		newNEndA=substrings[fromSubstring].nEndA;
		newPeriod=substrings[fromSubstring].period;
		newStrongPeriodSignal=substrings[fromSubstring].strongPeriodSignal;
		while (top>=0) {
			i=stack[top--];
			lengthI=(substrings[i].nEndA!=0?substrings[i].sumEndA/substrings[i].nEndA:substrings[i].endA)-(substrings[i].nStartA!=0?substrings[i].sumStartA/substrings[i].nStartA:substrings[i].startA)+1;
			for (j=fromSubstring+1; j<=lastSubstring; j++) {
				if (j==i || substrings[j].isMerged) continue;
				lengthJ=(substrings[j].nEndA!=0?substrings[j].sumEndA/substrings[j].nEndA:substrings[j].endA)-(substrings[j].nStartA!=0?substrings[j].sumStartA/substrings[j].nStartA:substrings[j].startA)+1;
				if ( Intervals.areApproximatelyIdentical(substrings[i].startA,substrings[i].endA,substrings[j].startA,substrings[j].endA) &&
				     !( Math.abs(substrings[i].startA,substrings[j].startA)<=identityThreshold &&
				        Math.min(lengthI,lengthJ)<=mergePrefixSuffixThreshold*Math.max(lengthI,lengthJ) &&
				        ( substrings[i].differentReplication(substrings[j]) ||
				          ( substrings[i].suffixReplication && substrings[j].suffixReplication &&
				            Math.abs((substrings[i].nEndA!=0?substrings[i].sumEndA/substrings[i].nEndA:substrings[i].endA)-substrings[i].minSuffixLength,(substrings[j].nEndA!=0?substrings[j].sumEndA/substrings[j].nEndA:substrings[j].endA)-substrings[j].minSuffixLength)>identityThreshold
				          )
				        )
				      ) &&
				     !( Math.min(lengthI,lengthJ)<=mergePrefixSuffixThreshold*Math.max(lengthI,lengthJ) &&
				        Math.abs(substrings[i].endA,substrings[j].endA)<=identityThreshold &&
				        ( substrings[i].differentReplication(substrings[j]) ||
				          ( substrings[i].prefixReplication && substrings[j].prefixReplication &&
				            Math.abs((substrings[i].nStartA!=0?substrings[i].sumStartA/substrings[i].nStartA:substrings[i].startA)+substrings[i].minPrefixLength,(substrings[j].nStartA!=0?substrings[j].sumStartA/substrings[j].nStartA:substrings[j].startA)+substrings[j].minPrefixLength)>identityThreshold
				          )
				        )
				      )
				   ) {
   					substrings[j].isMerged=true;
					substrings[j].representative=substrings[fromSubstring];
					if (!substrings[j].isWeak) substrings[fromSubstring].isWeak=false;
					stack[++top]=j;
					substrings[fromSubstring].isMerged=true;
   					if (substrings[j].pathLength>substrings[fromSubstring].pathLength) substrings[fromSubstring].pathLength=substrings[j].pathLength;
   					if (substrings[j].startA<newStartA) {
   						newStartA=substrings[j].startA;
   						substrings[fromSubstring].leftAlignment=substrings[j].leftAlignment;
   					}
					if (substrings[j].endA>newEndA) {
   						newEndA=substrings[j].endA;
   						substrings[fromSubstring].rightAlignment=substrings[j].rightAlignment;
   					}
   					newSumStartA+=substrings[j].sumStartA;
   					newNStartA+=substrings[j].nStartA;
   					newSumEndA+=substrings[j].sumEndA;
   					newNEndA+=substrings[j].nEndA;
					if (substrings[j].period<newPeriod) newPeriod=substrings[j].period;
					if (substrings[j].strongPeriodSignal) newStrongPeriodSignal=true;
   					if (substrings[j].substringReplication) {
						if (substrings[fromSubstring].substringReplication) {
							if (substrings[j].minPrefixLength<substrings[fromSubstring].minPrefixLength) substrings[fromSubstring].minPrefixLength=substrings[j].minPrefixLength;
						}
						else {
	   						substrings[fromSubstring].substringReplication=true;
	   						substrings[fromSubstring].prefixReplication=false; 
							substrings[fromSubstring].minPrefixLength=Math.min(substrings[j].minPrefixLength,substrings[fromSubstring].minPrefixLength>0?substrings[fromSubstring].minPrefixLength:Math.POSITIVE_INFINITY);
							substrings[fromSubstring].minPrefixLength=Math.min(substrings[fromSubstring].minPrefixLength,substrings[fromSubstring].minSuffixLength>0?substrings[fromSubstring].minSuffixLength:Math.POSITIVE_INFINITY);
							substrings[fromSubstring].suffixReplication=false; 
							substrings[fromSubstring].minSuffixLength=-1;
						}
   					}
   					else {
						if (!substrings[fromSubstring].substringReplication) {
	   						if (substrings[j].prefixReplication) {
	   							if (!substrings[fromSubstring].prefixReplication || substrings[fromSubstring].minPrefixLength<substrings[j].minPrefixLength) substrings[fromSubstring].minPrefixLength=substrings[j].minPrefixLength;
	   							substrings[fromSubstring].prefixReplication=true;
	   							if (!substrings[j].prefixesAreUniformlyDistributed) substrings[fromSubstring].prefixesAreUniformlyDistributed=false;
	   							if (substrings[fromSubstring].minPrefixDistanceBias==-1 || substrings[j].minPrefixDistanceBias<substrings[fromSubstring].minPrefixDistanceBias) substrings[fromSubstring].minPrefixDistanceBias=substrings[j].minPrefixDistanceBias;
	   							if (substrings[fromSubstring].maxPrefixDistanceBias==-1 || substrings[j].maxPrefixDistanceBias>substrings[fromSubstring].maxPrefixDistanceBias) substrings[fromSubstring].maxPrefixDistanceBias=substrings[j].maxPrefixDistanceBias;					
								substrings[fromSubstring].hasLengthBoundariesPrefix|=substrings[j].hasLengthBoundariesPrefix;
	   						}
	   						if (substrings[j].suffixReplication) {
	   							if (!substrings[fromSubstring].suffixReplication || substrings[j].minSuffixLength<substrings[fromSubstring].minSuffixLength) substrings[fromSubstring].minSuffixLength=substrings[j].minSuffixLength;
	   							substrings[fromSubstring].suffixReplication=true;
	   							if (!substrings[j].suffixesAreUniformlyDistributed) substrings[fromSubstring].suffixesAreUniformlyDistributed=false;
	   							if (substrings[fromSubstring].minSuffixDistanceBias==-1 || substrings[j].minSuffixDistanceBias<substrings[fromSubstring].minSuffixDistanceBias) substrings[fromSubstring].minSuffixDistanceBias=substrings[j].minSuffixDistanceBias;
	   							if (substrings[fromSubstring].maxSuffixDistanceBias==-1 || substrings[j].maxSuffixDistanceBias>substrings[fromSubstring].maxSuffixDistanceBias) substrings[fromSubstring].maxSuffixDistanceBias=substrings[j].maxSuffixDistanceBias;					
								substrings[fromSubstring].hasLengthBoundariesSuffix|=substrings[j].hasLengthBoundariesSuffix;
	   						}
	   						substrings[fromSubstring].nImpliedAlignmentsPrefix+=substrings[j].nImpliedAlignmentsPrefix;
	   						substrings[fromSubstring].nImpliedAlignmentsSuffix+=substrings[j].nImpliedAlignmentsSuffix;
						}
						else {
							if (substrings[j].minPrefixLength>0 && substrings[j].minPrefixLength<substrings[fromSubstring].minPrefixLength) substrings[fromSubstring].minPrefixLength=substrings[j].minPrefixLength;
							if (substrings[j].minSuffixLength>0 && substrings[j].minSuffixLength<substrings[fromSubstring].minPrefixLength) substrings[fromSubstring].minPrefixLength=substrings[j].minSuffixLength;
						}
   					}
   					// $substrings[fromSubstring].singleDeletionReplication$ is not 
					// set yet
				}
			}
		}
		substrings[fromSubstring].startA=newStartA;
		substrings[fromSubstring].endA=newEndA;
		substrings[fromSubstring].sumStartA=newSumStartA;
		substrings[fromSubstring].nStartA=newNStartA;
		substrings[fromSubstring].sumEndA=newSumEndA;
		substrings[fromSubstring].nEndA=newNEndA;
		substrings[fromSubstring].period=newPeriod;
		substrings[fromSubstring].strongPeriodSignal=newPeriod>0?newStrongPeriodSignal:false;
	}


	/**
	 * Let $(i,j,k)$ be a dense substring in $substrings$ that replicates by taking
	 * prefixes of itself, where $[i..k]$ is the span of the substring and $[j..k]$ is the
	 * range of the observed ending positions of its prefixes. Let $(i',j',k')$ be another
	 * substring of prefix type in $substrings$. We say that the two dense substrings are
	 * equivalent iff the distance between $i$ and $i'$ is small, and if $(j,k)$
	 * intersects $(j',k')$ or is adjacent to it. A similar definition holds for dense
	 * substrings that replicate by taking suffixes of themselves.
	 *
	 * Let $G$ be the graph whose vertices are the elements of $substrings$, and whose
	 * edges are equivalence relations. The procedure finds the connected component of
	 * $fromSubstring$ in $G$. Every dense substring in such component, except
	 * $fromSubstring$, is assigned the same value of $representative$ (which equals
	 * $substrings[fromSubstring]$). Elements of $Events.events$ that are close to
	 * concatenation junctions are marked as to be removed. The boundaries of the
	 * representative substring of the connected component are suitably updated.
	 *
	 * Remark: the procedure assumes that all substrings in $substrings[0..fromSubstring
	 * -1]$ have already been assigned a connected component different from the one of
	 * $fromSubstring$, and that the $representative$ field of every substring equals
	 * either the substring itself or another substring.
	 *
	 * Remark: the procedure handles substrings that replicate by taking both prefixes and
	 * suffixes of themselves.
	 *
	 * Remark: if a connected component contains a non-weak substring, the representative
	 * of the component is marked as non-weak.
	 *
	 * Remark: the period of the representative substring is set to the minimum of the 
	 * periods of the substrings in its component, without any check for compatibility.
	 *
	 * @return the number of substrings in the connected component of $fromSubstring$,
	 * excluding $fromSubstring$.
 	 */
	private static final int getConnectedComponent_concatenation(int fromSubstring) {
		final int BOUNDARY_ALIGNMENTS_THRESHOLD = IO.coverage*50;  // Arbitrary
		boolean inConnectedComponent, newStrongPeriodSignal;
		int i, j, top, count;
		int xiPrefix, yiPrefix, xjPrefix, yjPrefix, xiSuffix, yiSuffix, xjSuffix, yjSuffix;
		int newStartA, newEndA, newPeriod;
		DenseSubstring representative;

		count=0;
		top=0; stack[0]=fromSubstring;
		representative=substrings[fromSubstring].representative;
		representative.previousSumStartA=representative.sumStartA;
		representative.previousNStartA=representative.nStartA;
		representative.previousSumEndA=representative.sumEndA;
		representative.previousNEndA=representative.nEndA;
		newStartA=representative.startA;
		newEndA=representative.endA;
		newPeriod=representative.period;
		newStrongPeriodSignal=representative.strongPeriodSignal;
		while (top>=0) {
			i=stack[top];
			top--;
			xiPrefix=-1; yiPrefix=-1; xiSuffix=-1; yiSuffix=-1;
			if (substrings[i].prefixReplication) {
				xiPrefix=substrings[i].sumStartA/substrings[i].nStartA+substrings[i].minPrefixLength;
				yiPrefix=substrings[i].endA;
			}
			else if (substrings[i].suffixReplication) {
				xiSuffix=substrings[i].startA;
				yiSuffix=substrings[i].sumEndA/substrings[i].nEndA-substrings[i].minSuffixLength;
			}
			for (j=fromSubstring+1; j<=lastSubstring; j++) {
				if (j==i || substrings[j].representative!=substrings[j]) continue;
				inConnectedComponent=false;
				xjPrefix=-1; yjPrefix=-1; xjSuffix=-1; yjSuffix=-1;
				if (substrings[j].prefixReplication) {
					xjPrefix=substrings[j].sumStartA/substrings[j].nStartA+substrings[j].minPrefixLength;
					yjPrefix=substrings[j].endA;
				}
				else if (substrings[j].suffixReplication) {
					xjSuffix=substrings[j].startA;
					yjSuffix=substrings[j].sumEndA/substrings[j].nEndA-substrings[j].minSuffixLength;
				}
				if ( substrings[i].prefixReplication && substrings[j].prefixReplication && !substrings[j].suffixReplication &&
				     Math.abs(substrings[i].sumStartA/substrings[i].nStartA,substrings[j].sumStartA/substrings[j].nStartA)<=identityThreshold &&
					 ( (yiPrefix<yjPrefix-identityThreshold && substrings[i].previousSumStartA<BOUNDARY_ALIGNMENTS_THRESHOLD) ||
					   (yjPrefix<yiPrefix-identityThreshold && substrings[j].previousSumStartA<BOUNDARY_ALIGNMENTS_THRESHOLD)
					 ) &&
				     (Intervals.intersectionLength(xiPrefix,yiPrefix,xjPrefix,yjPrefix)>0 || xjPrefix-yiPrefix<=identityThreshold || xiPrefix-yjPrefix<=identityThreshold) 
				   ) {
					inConnectedComponent=true;
					substrings[j].representative=representative;
					representative.isConcatenated=true;
					substrings[i].isConcatenated=true;
					substrings[j].isConcatenated=true;
					if (!substrings[j].isWeak) representative.isWeak=false;
					if (substrings[j].minPrefixLength<representative.minPrefixLength) representative.minPrefixLength=substrings[j].minPrefixLength;
					representative.nImpliedAlignmentsPrefix+=substrings[j].nImpliedAlignmentsPrefix;
					if (!substrings[j].prefixesAreUniformlyDistributed) representative.prefixesAreUniformlyDistributed=false;
					if (representative.minPrefixDistanceBias==-1 || substrings[j].minPrefixDistanceBias<representative.minPrefixDistanceBias) representative.minPrefixDistanceBias=substrings[j].minPrefixDistanceBias;
					if (representative.maxPrefixDistanceBias==-1 || substrings[j].maxPrefixDistanceBias>representative.maxPrefixDistanceBias) representative.maxPrefixDistanceBias=substrings[j].maxPrefixDistanceBias;
					representative.hasLengthBoundariesPrefix|=substrings[j].hasLengthBoundariesPrefix;
					if (substrings[j].startA<newStartA) newStartA=substrings[j].startA;
					if (substrings[j].endA>newEndA) newEndA=substrings[j].endA;
					if (substrings[j].period<newPeriod) newPeriod=substrings[j].period;
					if (substrings[j].strongPeriodSignal) newStrongPeriodSignal=true;
				}
				if ( substrings[i].suffixReplication && substrings[j].suffixReplication && !substrings[j].prefixReplication &&
				     Math.abs(substrings[i].endA,substrings[j].endA)<=identityThreshold &&
					 ( (xiSuffix>xjSuffix+identityThreshold && substrings[i].previousSumStartA<BOUNDARY_ALIGNMENTS_THRESHOLD) ||
					   (xjSuffix>xiSuffix+identityThreshold && substrings[j].previousSumStartA<BOUNDARY_ALIGNMENTS_THRESHOLD)
					 ) &&
					 (Intervals.intersectionLength(xiSuffix,yiSuffix,xjSuffix,yjSuffix)>0 || xjSuffix-yiSuffix<=identityThreshold || xiSuffix-yjSuffix<=identityThreshold) 
				   ) {
					inConnectedComponent=true;
					substrings[j].representative=representative;
					representative.isConcatenated=true;
					substrings[i].isConcatenated=true;
					substrings[j].isConcatenated=true;
					if (!substrings[j].isWeak) representative.isWeak=false;
					if (substrings[j].minSuffixLength<representative.minSuffixLength) representative.minSuffixLength=substrings[j].minSuffixLength;
					representative.nImpliedAlignmentsSuffix+=substrings[j].nImpliedAlignmentsSuffix;
					if (!substrings[j].suffixesAreUniformlyDistributed) representative.suffixesAreUniformlyDistributed=false;
					if (representative.minSuffixDistanceBias==-1 || substrings[j].minSuffixDistanceBias<representative.minSuffixDistanceBias) representative.minSuffixDistanceBias=substrings[j].minSuffixDistanceBias;
					if (representative.maxSuffixDistanceBias==-1 || substrings[j].maxSuffixDistanceBias>representative.maxSuffixDistanceBias) representative.maxSuffixDistanceBias=substrings[j].maxSuffixDistanceBias;
					representative.hasLengthBoundariesSuffix|=substrings[j].hasLengthBoundariesSuffix;
					if (substrings[j].startA<newStartA) newStartA=substrings[j].startA;
					if (substrings[j].endA>newEndA) newEndA=substrings[j].endA;
					if (substrings[j].period<newPeriod) newPeriod=substrings[j].period;
					if (substrings[j].strongPeriodSignal) newStrongPeriodSignal=true;
				}
				if (!inConnectedComponent) continue;
				count++;
				stack[++top]=j;
			}
		}
		representative.startA=newStartA;
		representative.endA=newEndA;
		representative.period=newPeriod;
		representative.strongPeriodSignal=newPeriod>0?newStrongPeriodSignal:false;
		return count;
	}


	/**
	 * Updates the $sum*A$ and $n*A$ values of concatenated representative substrings.
	 * Only $substrings[from..to]$ is considered.
	 *
	 * Remark: this is done in a single pass, after the connected components have been
	 * built, since building a connected component moves both the left and the right
	 * boundary of a representative string, and only points that are close to the final
	 * boundaries should contribute to its $sum*A$ and $n*A$.
	 */
	private static final void updateConcatenatedSumN(int from, int to) {
		final int THRESHOLD = identityThreshold;
		int i;
		
		// Deleting $sumStartA/nStartA$ values of representatives if they are too far from
		// their new $startA$.
		for (i=from; i<=to; i++) {
			if (substrings[i].representative!=substrings[i]) continue;
			if (Math.abs(substrings[i].startA,substrings[i].sumStartA/substrings[i].nStartA)>THRESHOLD) {
				substrings[i].sumStartA=0;
				substrings[i].nStartA=0;
			}
			if (Math.abs(substrings[i].endA,substrings[i].sumEndA/substrings[i].nEndA)>THRESHOLD) {
				substrings[i].sumEndA=0;
				substrings[i].nEndA=0;
			}
		}
		
		// Updating $sumStartA,nStartA$.
		for (i=from; i<=to; i++) {
			if (substrings[i].representative==substrings[i]) continue;
			if (Math.abs(substrings[i].startA,substrings[i].representative.startA)<=THRESHOLD) {
				substrings[i].representative.sumStartA+=substrings[i].sumStartA;
				substrings[i].representative.nStartA+=substrings[i].nStartA;
			}
			if (Math.abs(substrings[i].endA,substrings[i].representative.endA)<=THRESHOLD) {
				substrings[i].representative.sumEndA+=substrings[i].sumEndA;
				substrings[i].representative.nEndA+=substrings[i].nEndA;
			}
			if (substrings[i].maximalStartA>=substrings[i].representative.startA && substrings[i].maximalStartA<substrings[i].representative.maximalStartA) substrings[i].representative.maximalStartA=substrings[i].maximalStartA;
			if (substrings[i].maximalEndA<=substrings[i].representative.endA && substrings[i].maximalEndA>substrings[i].representative.maximalEndA) substrings[i].representative.maximalEndA=substrings[i].maximalEndA;
		}
		
		// Ensuring nonzero $sumStartA,nStartA$.
		for (i=from; i<=to; i++) {
			if (substrings[i].representative!=substrings[i]) continue;
			if (substrings[i].sumStartA==0 || substrings[i].nStartA==0) {
				substrings[i].sumStartA=substrings[i].startA;
				substrings[i].nStartA=1;
			}
			if (substrings[i].sumEndA==0 || substrings[i].nEndA==0) {
				substrings[i].sumEndA=substrings[i].endA;
				substrings[i].nEndA=1;
			}
		}		
	}


	/**
	 * Updates the replication mode of concatenated representative substrings, using all
	 * other substrings in their equivalence class.
	 * A substring with prefix replication can make its representative of prefix and
	 * suffix replication, or of substring replication, depending on their overlap.
	 * A symmetrical argument holds for substrings with suffix replication.
	 *
	 * Remark: this procedure could create dense substrings of substring type with a 
	 * period. Such period will be overwritten by procedure $estimateSubstringPeriods$.
	 */
	private static final void updateConcatenatedReplicationModes() {
		for (int i=0; i<=lastSubstring; i++) {
			if (substrings[i].representative==substrings[i]) continue;
			if (substrings[i].substringReplication) {
				if (!substrings[i].representative.substringReplication) {
					substrings[i].representative.substringReplication=true;
					substrings[i].representative.prefixReplication=false;
					substrings[i].representative.minPrefixLength=Math.min(substrings[i].minPrefixLength,substrings[i].representative.minPrefixLength>0?substrings[i].representative.minPrefixLength:Math.POSITIVE_INFINITY);
					substrings[i].representative.minPrefixLength=Math.min(substrings[i].representative.minPrefixLength,substrings[i].representative.minSuffixLength>0?substrings[i].representative.minSuffixLength:Math.POSITIVE_INFINITY);
					substrings[i].representative.suffixReplication=false;
					substrings[i].representative.minSuffixLength=-1;
				}
				else {
					if (substrings[i].minPrefixLength<substrings[i].representative.minPrefixLength) substrings[i].representative.minPrefixLength=substrings[i].minPrefixLength;
				}
			}
			else if ( substrings[i].prefixReplication &&
					  !substrings[i].representative.prefixReplication &&
					  substrings[i].representative.suffixReplication &&
					  !substrings[i].representative.substringReplication ) {
				if (Intervals.jaccardSimilarity(substrings[i].startA,substrings[i].endA,substrings[i].representative.startA,substrings[i].representative.endA)>=Intervals.jaccardThreshold) {
					substrings[i].representative.prefixReplication=true;
					substrings[i].representative.minPrefixLength=substrings[i].minPrefixLength;
				}
				else {
					substrings[i].representative.substringReplication=true;
					substrings[i].representative.prefixReplication=false;
					substrings[i].representative.minPrefixLength=Math.min(substrings[i].minPrefixLength,substrings[i].representative.minSuffixLength);						
					substrings[i].representative.minPrefixLength=Math.min(substrings[i].representative.minPrefixLength,substrings[i].minSuffixLength>0?substrings[i].minSuffixLength:Math.POSITIVE_INFINITY);
					substrings[i].representative.suffixReplication=false;
					substrings[i].representative.minSuffixLength=-1;
				}
			}
			else if ( substrings[i].suffixReplication &&
				      substrings[i].representative.prefixReplication &&
				      !substrings[i].representative.suffixReplication &&
				      !substrings[i].representative.substringReplication ) {
				if (Intervals.jaccardSimilarity(substrings[i].startA,substrings[i].endA,substrings[i].representative.startA,substrings[i].representative.endA)>=Intervals.jaccardThreshold) {
					substrings[i].representative.suffixReplication=true;
					substrings[i].representative.minSuffixLength=substrings[i].minSuffixLength;
				}
				else {
					substrings[i].representative.substringReplication=true;
					substrings[i].representative.prefixReplication=false;
					substrings[i].representative.minPrefixLength=Math.min(substrings[i].minSuffixLength,substrings[i].representative.minPrefixLength);
					substrings[i].representative.minPrefixLength=Math.min(substrings[i].representative.minPrefixLength,substrings[i].minPrefixLength>0?substrings[i].minPrefixLength:Math.POSITIVE_INFINITY);
					substrings[i].representative.suffixReplication=false;
					substrings[i].representative.minSuffixLength=-1;
				}
			}
		}
	}


	/**
	 * Checks whether the pair $(prefixSubstring,suffixSubstring)$ is a repeat that
	 * replicates by a single internal deletion. This is done by scanning the set of all
	 * alignments implied by $prefixSubstring$ and $suffixSubstring$, checking whether the
	 * Jaccard similarity of the two sets is large enough. An alignment in a set belongs
	 * to the "intersection" iff there is a compatible alignment in the other set.
	 *
	 * A repeat that replicates by a single internal deletion could use multiple distinct
	 * deletion lengths. The procedure tries to estimate the number of such lengths by
	 * fitting a density estimation tree on the set of all observed deletion lengths,
	 * assuming that such lengths are spaced apart by at least $minDeletionLength$.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments[0..lastAlignment]$ to be
	 * sorted by the implying dense substring, readB, orientation, startB, and it assumes
	 * $DenseSubstring.firstAlignment$ to have been already computed and to be nonnegative
	 * (indeed, every dense substring has at least $minImpliedAlignments$ implied
	 * alignments).
	 *
	 * Remark: a dense substring could have no implied alignment, since alignments that 
	 * are implied both by a prefix and by a suffix substring are reassigned at the end of
	 * $getSubstrings$.
	 *
	 * @param tmpPeriodic temporary space;
	 * @return $singleDeletionTmp[0]$: the number of distinct deletion lengths; or -1 if
	 * the pair $(prefixSubstring,suffixSubstring)$ does not replicate by a single
	 * internal deletion; or 0 if it does, but we couldn't estimate the number of distinct
	 * deletion lengths. $singleDeletionTmp[1]$: shortest deletion (center of mass of the
	 * leftmost peak in the DET). $singleDeletionTmp[2]$: longest deletion (center of
	 * mass of the rightmost peak in the DET).
	 */
	private static final void isSingleDeletion(DenseSubstring prefixSubstring, DenseSubstring suffixSubstring, int lastAlignment, PeriodicSubstring tmpPeriodic) {
		final int MIN_INTERVAL_LENGTH = minDeletionLength;
		final int MIN_LOCAL_MAX_DISTANCE = minDeletionLength;
		final int MIN_HOLE_LENGTH = MIN_LOCAL_MAX_DISTANCE;
		final int BIN_LENGTH = minDeletionLength<<1;
		final int MAX_DISTANCE = IO.quantum;  // Used to estimate single deletions
		final int COVERAGE_THRESHOLD = IO.minRepeatCoverage;
		boolean straddling, eventsInStraddlingPart;
		int i, j, p;
		int firstJForNextI, delta, deletion;
		int readBPrefix, readBSuffix, orientationPrefix, orientationSuffix;
		int startBPrefix, startBSuffix, endBPrefix, endBSuffix;
		int nOccurrences, nOccurrencesPrime, nLocalMaximumLeaves;
		int straddlingPartStart, straddlingPartEnd;


if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("isSingleDeletion> prefixSubstring="+prefixSubstring.id+"="+prefixSubstring);
	IO.printErr("isSingleDeletion> alignments:");
	for (i=0; i<=lastAlignment; i++) IO.printErr("["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"] impliedByDenseSubstring="+ReadA.sortedAlignments[i].impliedByDenseSubstring.id+" readB="+Alignments.alignments[ReadA.sortedAlignments[i].id][1]+" orientation="+Alignments.alignments[ReadA.sortedAlignments[i].id][2]+" startB="+Alignments.alignments[ReadA.sortedAlignments[i].id][5]);
	IO.printErr("prefixSubstring.firstAlignment="+prefixSubstring.firstAlignment+" suffixSubstring.firstAlignment="+suffixSubstring.firstAlignment);
}
		
		
		Math.set(singleDeletionTmp,2,-1);
		straddling=prefixSubstring.endA>=suffixSubstring.startA;
		straddlingPartStart=straddling?suffixSubstring.startA:-1;
		straddlingPartEnd=straddling?prefixSubstring.endA:-1;
		if (!straddling && ReadA.getAverageCoverage(prefixSubstring.endA+1,suffixSubstring.startA-1)<COVERAGE_THRESHOLD) return;

		// Initializing markings
		i=prefixSubstring.firstAlignment;
		if (i==-1) return;
		while (i<=lastAlignment && ReadA.sortedAlignments[i].impliedByDenseSubstring==prefixSubstring) {
			ReadA.sortedAlignments[i].markedPrefix=false;
			i++;
		}
		j=suffixSubstring.firstAlignment;
		if (j==-1) return;
		while (j<=lastAlignment && ReadA.sortedAlignments[j].impliedByDenseSubstring==suffixSubstring) {
			ReadA.sortedAlignments[j].markedSuffix=false;
			j++;
		}

		// Collecting deletion lengths
		i=prefixSubstring.firstAlignment;
		j=suffixSubstring.firstAlignment;
		firstJForNextI=-1;
		lastSingleDeletion=-1;
		eventsInStraddlingPart=false;
		while (i<=lastAlignment && ReadA.sortedAlignments[i].impliedByDenseSubstring==prefixSubstring) {
			p=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (p>straddlingPartStart+identityThreshold && p<straddlingPartEnd-identityThreshold) eventsInStraddlingPart=true;
			if (j>lastAlignment || ReadA.sortedAlignments[j].impliedByDenseSubstring!=suffixSubstring) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			p=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
			if (p>straddlingPartStart+identityThreshold && p<straddlingPartEnd-identityThreshold) eventsInStraddlingPart=true;
			readBPrefix=Alignments.alignments[ReadA.sortedAlignments[i].id][1];
			readBSuffix=Alignments.alignments[ReadA.sortedAlignments[j].id][1];
			orientationPrefix=Alignments.alignments[ReadA.sortedAlignments[i].id][2];
			orientationSuffix=Alignments.alignments[ReadA.sortedAlignments[j].id][2];
			startBPrefix=Alignments.alignments[ReadA.sortedAlignments[i].id][5];
			startBSuffix=Alignments.alignments[ReadA.sortedAlignments[j].id][5];
			endBPrefix=Alignments.alignments[ReadA.sortedAlignments[i].id][6];
			endBSuffix=Alignments.alignments[ReadA.sortedAlignments[j].id][6];
			if (readBSuffix>readBPrefix || (orientationPrefix==1&&orientationSuffix==0) || orientationPrefix==1?startBSuffix>endBPrefix+MAX_DISTANCE:startBSuffix>=startBPrefix) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (readBSuffix<readBPrefix || (orientationPrefix==0&&orientationSuffix==1) || (orientationPrefix==1?startBSuffix<startBPrefix:endBSuffix<startBPrefix-MAX_DISTANCE)) {
				j++;
				continue;
			}
			if ( i<lastAlignment && firstJForNextI==-1 &&
				 ReadA.sortedAlignments[i+1].impliedByDenseSubstring==prefixSubstring &&
			     Alignments.alignments[ReadA.sortedAlignments[i+1].id][1]==readBPrefix &&
			     Alignments.alignments[ReadA.sortedAlignments[i+1].id][2]==orientationPrefix &&
			     (orientationPrefix==1?startBSuffix>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][5]:endBSuffix>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][5]-MAX_DISTANCE) ) firstJForNextI=j;
			delta=orientationPrefix==1?Math.abs(startBSuffix,endBPrefix):Math.abs(startBPrefix,endBSuffix);
			deletion=Alignments.alignments[ReadA.sortedAlignments[j].id][3]-Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (delta<=MAX_DISTANCE && deletion>=minDeletionLength) {
				ReadA.sortedAlignments[i].markedPrefix=true;
				ReadA.sortedAlignments[j].markedSuffix=true;
				lastSingleDeletion++;
				singleDeletions[lastSingleDeletion].position=deletion;
				singleDeletions[lastSingleDeletion].mass=1;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
			}
			else j++;
		}
		if (lastSingleDeletion==-1) return;
		if (straddling && !eventsInStraddlingPart && inPeriodicSubstring(prefixSubstring.startA,straddlingPartStart-1,tmpPeriodic) && inPeriodicSubstring(straddlingPartEnd+1,suffixSubstring.endA,tmpPeriodic)) return;

		// Deciding whether the pair of substrings represents a single deletion
		nOccurrences=0;
		i=prefixSubstring.firstAlignment;
		while (i<=lastAlignment && ReadA.sortedAlignments[i].impliedByDenseSubstring==prefixSubstring) {
			if (ReadA.sortedAlignments[i].markedPrefix) nOccurrences++;
			i++;
		}
		nOccurrencesPrime=0;
		j=suffixSubstring.firstAlignment;
		while (j<=lastAlignment && ReadA.sortedAlignments[j].impliedByDenseSubstring==suffixSubstring) {
			if (ReadA.sortedAlignments[j].markedSuffix) nOccurrencesPrime++;
			j++;
		}
		if (nOccurrences<MIN_DELETION_OCCURRENCES && nOccurrencesPrime<MIN_DELETION_OCCURRENCES) return;

		// Building a density estimation tree on the set of all deletion lengths
		if (lastSingleDeletion==0) {
			singleDeletionTmp[0]=1;
			singleDeletionTmp[1]=(int)singleDeletions[0].position;
			singleDeletionTmp[2]=singleDeletionTmp[1];
			return;
		}
		lastSingleDeletion=Points.sortAndCompact(singleDeletions,lastSingleDeletion);
		if (lastSingleDeletion==0) {
			singleDeletionTmp[0]=1;
			singleDeletionTmp[1]=(int)singleDeletions[0].position;
			singleDeletionTmp[2]=singleDeletionTmp[1];
			return;
		}


if (IO.SHOW_STD_ERR) IO.printErr("isSingleDeletion> lengths:");
if (IO.SHOW_STD_ERR) { for (i=0; i<=lastSingleDeletion; i++) IO.printErr(singleDeletions[i].position+","+singleDeletions[i].mass); }


		if (lastSingleDeletion+1<=Points.FEW_POINTS) {
			if (singleDeletions[lastSingleDeletion].position-singleDeletions[0].position<=MAX_DELETION_DIFFERENCE) {
				singleDeletionTmp[0]=1;
				singleDeletionTmp[1]=Math.round(Points.getCenterOfMass(singleDeletions,0,lastSingleDeletion,false,-1,lastSingleDeletion));
				singleDeletionTmp[2]=singleDeletionTmp[1];
				return;
			}
			singleDeletionTmp[0]=0;
			singleDeletionTmp[1]=(int)singleDeletions[0].position;
			singleDeletionTmp[2]=(int)singleDeletions[lastSingleDeletion].position;
			return;
		}
		if (Points.areUniformlyDistributed(singleDeletions,0,lastSingleDeletion,true,BIN_LENGTH)) {
			singleDeletionTmp[0]=1;
			singleDeletionTmp[1]=Math.round(Points.getCenterOfMass(singleDeletions,0,lastSingleDeletion,false,-1,lastSingleDeletion));
			singleDeletionTmp[2]=singleDeletionTmp[1];
			return;
		}
  		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(singleDeletions,0,lastSingleDeletion,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE,true,-1,MIN_HOLE_LENGTH,-1,false,true);
		if (nLocalMaximumLeaves==0) {
			singleDeletionTmp[0]=0;
			singleDeletionTmp[1]=(int)singleDeletions[0].position;
			singleDeletionTmp[2]=(int)singleDeletions[lastSingleDeletion].position;
			return;
		}
		DensityEstimationTree.markRunsOfLocalMaximumLeaves(singleDeletions);

		// Computing the center of mass of the first and last run of local-maximum leaves
		singleDeletionTmp[0]=DensityEstimationTree.nRuns;
		singleDeletionTmp[1]=Math.round(DensityEstimationTree.getCenterOfMassOfRun(0,singleDeletions,lastSingleDeletion,-1));
		singleDeletionTmp[2]=Math.round(DensityEstimationTree.getCenterOfMassOfRun(DensityEstimationTree.nRuns-1,singleDeletions,lastSingleDeletion,-1));
	}


	/**
	 * Adds to $Events.events$ all maximal first/last positions of alignments implied by
	 * dense substrings, that are close enough to the first/last position of their
	 * implying dense substring, and that are not contained in a maximal range of short-
	 * period intervals, or in a dense substring of substring type.
	 *
	 * Remark: the procedure does not affect dense substrings that replicate by taking
	 * substrings of themselves.
	 *
	 * Remark: the procedure assumes that $ReadA.sortedAlignments[0..lastAlignment]$
	 * contains all and only the alignments implied by dense substrings, sorted by first
	 * position.
	 */
	private static final void addEvents(int lastAlignment) {
		int i, j;
		int pos, first, last, rangeFirst, rangeLast;
		DenseSubstring substring;
		PeriodicSubstringInterval tmpInterval = new PeriodicSubstringInterval();
		
		// Ensuring the necessary orders
		if (DenseSubstring.order!=DenseSubstring.STARTA) {
			DenseSubstring.order=DenseSubstring.STARTA;
			if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		}
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		}

		// Initializing substrings
		for (i=0; i<=lastAlignment; i++) {
			substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (substring!=null) {
				substring.startAAdded=false;
				substring.endAAdded=false;
			}
			ReadA.sortedAlignments[i].markedPrefix=false;  // Used as a flag
			ReadA.sortedAlignments[i].markedSuffix=false;  // Used as a flag
		}
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].leftAlignment=0;  // Used as a flag
			substrings[i].rightAlignment=0;  // Used as a flag
		}

		// Marking alignments that start/end inside a short-period range.
		if (PeriodicSubstrings.lastInterval>=0) {
			rangeFirst=PeriodicSubstrings.intervals[0].firstPosition;
			rangeLast=PeriodicSubstrings.intervals[0].lastPosition;
			j=0;
			for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
				if (PeriodicSubstrings.intervals[i].firstPosition<=rangeLast+identityThreshold) {
					rangeLast=Math.max(rangeLast,PeriodicSubstrings.intervals[i].lastPosition);
					continue;
				}
				j=addEvents_markAlignments(j,lastAlignment,rangeFirst,rangeLast);
				rangeFirst=PeriodicSubstrings.intervals[i].firstPosition;
				rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
			}
			addEvents_markAlignments(j,lastAlignment,rangeFirst,rangeLast);
		}
		
		// Marking alignments that start/end inside a substring type.
		j=0;
		for (i=0; i<=lastSubstring; i++) {
			if (!substrings[i].substringReplication || substrings[i].isWeak) continue;
			j=addEvents_markAlignments(j,lastAlignment,substrings[i].startA,substrings[i].endA);
		}

		// Adding alignment boundaries to $Events$
		for (i=0; i<=lastAlignment; i++) {
			substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (substring==null) continue;
			first=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			last=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if ( ReadA.sortedAlignments[i].isLeftMaximal==1 &&
			     Math.abs(first,substring.sumStartA/substring.nStartA)<=eventThreshold &&
			     (substring.prefixReplication || substring.singleDeletionReplication) &&
				 !ReadA.sortedAlignments[i].markedPrefix
			   ) {
				Events.lastEvent++;
				Events.ensureSpace_events(Events.lastEvent+1);
				Events.events[Events.lastEvent].clear();
				Events.events[Events.lastEvent].position=first;
				Events.events[Events.lastEvent].nOpen=1;
				ReadA.sortedAlignments[i].startAdded=true;
				substring.startAAdded=true;
if (IO.SHOW_STD_ERR) IO.printErr("addEvents()> added event "+first);
			}
			if ( ReadA.sortedAlignments[i].isRightMaximal==1 &&
			     Math.abs(last,substring.sumEndA/substring.nEndA)<=eventThreshold &&
			     (substring.suffixReplication || substring.singleDeletionReplication) &&
				 !ReadA.sortedAlignments[i].markedSuffix
			   ) {
				Events.lastEvent++;
				Events.ensureSpace_events(Events.lastEvent+1);
				Events.events[Events.lastEvent].clear();
				Events.events[Events.lastEvent].position=last;
				Events.events[Events.lastEvent].nClosed=1;
				ReadA.sortedAlignments[i].endAdded=true;
				substring.endAAdded=true;
if (IO.SHOW_STD_ERR) IO.printErr("addEvents()> added event "+last);
			}
		}
		
		// Marking substrings that start/end inside a short-period range.
		if (PeriodicSubstrings.lastInterval>=0) {
			rangeFirst=PeriodicSubstrings.intervals[0].firstPosition;
			rangeLast=PeriodicSubstrings.intervals[0].lastPosition;
			j=0;
			for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
				if (PeriodicSubstrings.intervals[i].firstPosition<=rangeLast+identityThreshold) {
					rangeLast=Math.max(rangeLast,PeriodicSubstrings.intervals[i].lastPosition);
					continue;
				}
				j=addEvents_markSubstrings(j,rangeFirst,rangeLast,false);
				rangeFirst=PeriodicSubstrings.intervals[i].firstPosition;
				rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
			}
			addEvents_markSubstrings(j,rangeFirst,rangeLast,false);
		}
		
		// Marking substrings that start/end inside a substring type.
		j=0;
		for (i=0; i<=lastSubstring; i++) {
			if (!substrings[i].substringReplication || substrings[i].isWeak) continue;
			j=addEvents_markSubstrings(j,substrings[i].startA,substrings[i].endA,false);
		}

		// Adding missing substring boundaries to $Factors.splits$
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].substringReplication) continue;
			if (!substrings[i].startAAdded && substrings[i].leftAlignment==0) {
				Factors.lastSplit++;
				Factors.splits[Factors.lastSplit].clear();
				if ( substrings[i].startA<=identityThreshold ||
				     Reads.isRandomInsertion(ReadA.id,substrings[i].startA-Reads.minLongInsertionLength,substrings[i].startA-1,true) ||
				     substrings[i].maximalStartA==-1 
				   ) pos=substrings[i].startA;
				else pos=substrings[i].maximalStartA;
				Factors.splits[Factors.lastSplit].position=pos;
				Factors.splits[Factors.lastSplit].nOpen=1;
				substrings[i].startAAdded=true;
if (IO.SHOW_STD_ERR) IO.printErr("addEvents()>1 added SPLIT "+Factors.splits[Factors.lastSplit].position);
			}
			if (!substrings[i].endAAdded && substrings[i].rightAlignment==0) {
				Factors.lastSplit++;
				Factors.splits[Factors.lastSplit].clear();
				if ( Reads.getReadLength(ReadA.id)-substrings[i].endA<=identityThreshold ||
				     Reads.isRandomInsertion(ReadA.id,substrings[i].endA+1,substrings[i].endA+Reads.minLongInsertionLength,true) ||
				     substrings[i].maximalEndA==-1 
				   ) pos=substrings[i].endA;
				else pos=substrings[i].maximalEndA;
				Factors.splits[Factors.lastSplit].position=pos;
				Factors.splits[Factors.lastSplit].nClosed=1;
				substrings[i].endAAdded=true;
if (IO.SHOW_STD_ERR) IO.printErr("addEvents()>2 added SPLIT "+Factors.splits[Factors.lastSplit].position);
			}
		}
	}
	
	
	/**
	 * Remark: the procedure uses the temporary flags $marked{Prefix,Suffix}$ of an 
	 * alignment.
	 *
	 * @param alignmentFrom first element of $ReadA.sortedAlignments$ to be used; 
	 * @return first element of $ReadA.sortedAlignments$ to be used for the next range.
	 */
	private static final int addEvents_markAlignments(int alignmentFrom, int alignmentTo, int rangeFirst, int rangeLast) {
		int j, firstJForNext;
		int firstPosition, lastPosition;
		
		j=alignmentFrom; firstJForNext=-1;
		while (j<=alignmentTo) {
			firstPosition=ReadA.sortedAlignments[j].startA();
			lastPosition=ReadA.sortedAlignments[j].endA();
			if (lastPosition<rangeFirst) {
				j++;
				continue;
			}
			if (firstPosition>rangeLast) break;
			if (firstJForNext==-1 && lastPosition>rangeLast) firstJForNext=j;
			if (firstPosition>=rangeFirst && firstPosition<=rangeLast) ReadA.sortedAlignments[j].markedPrefix=true;
			if (lastPosition>=rangeFirst && lastPosition<=rangeLast) ReadA.sortedAlignments[j].markedSuffix=true;
			j++;
		}
		return firstJForNext==-1?j:firstJForNext;
	}
	
	
	/**
	 * Remark: the procedure uses the temporary fields ${left,right}Alignment$ of a dense
	 * substring as flags.
	 *
	 * @param substringFrom first element of $substrings$ to be used; 
	 * @param substringType marks only substrings that are not (FALSE) or that are (TRUE)
	 * of substring type;
	 * @return first element of $substrings$ to be used for the next range.
	 */
	private static final int addEvents_markSubstrings(int substringFrom, int rangeFirst, int rangeLast, boolean substringType) {
		int j, firstJForNext;
		
		j=substringFrom; firstJForNext=-1;
		while (j<=lastSubstring) {
			if (substrings[j].endA<rangeFirst || substrings[j].substringReplication!=substringType) {
				j++;
				continue;
			}
			if (substrings[j].startA>rangeLast) break;
			if (firstJForNext==-1 && substrings[j].endA>rangeLast) firstJForNext=j;
			if (substrings[j].startA>=rangeFirst && substrings[j].startA<=rangeLast) substrings[j].leftAlignment=1;
			if (substrings[j].endA>=rangeFirst && substrings[j].endA<=rangeLast) substrings[j].rightAlignment=1;
			j++;
		}
		return firstJForNext==-1?j:firstJForNext;
	}


	/**
	 * It could happen that a dense substring contains more than one occurrence of a
	 * maximal repeat, for example because a concatenation of maximal repeats replicates
	 * by taking prefixes of the entire concatenation. The boundaries of all such repeats,
	 * in all dense substrings, will be detected in the following steps of the pipeline,
	 * by estimating the density of maximal events that come from alignments that are
	 * not implied by dense or periodic substrings. However, dense substrings that
	 * replicate by taking substrings of themselves need to be split by finding peaks in
	 * function $\delta$: see procedure $ReadA.getDeltaPeaks$.
	 *
	 * Remark: the procedure computes two distinct histograms, one for left-maximal and
	 * one for right-maximal events. This makes it easier to detect e.g. the start of a 
	 * simple repeat and the end of another simple repeat that are close to each other.
	 *
	 * Remark: this procedure is applied also to dense substrings of substring type that
	 * are created by procedure $mergeSubstrings$, and to weak dense substrings of
	 * substring type.
	 *
	 * Remark: the alignments, inside the dense substring, that align with the splitpoints
	 * will likely generate alignment intervals, but the alignment that do not align with 
	 * the splitpoints will likely still be assigned to the dense substring. This might
	 * make the dense substring an overlap fork.
	 *
	 * Remark: the procedure assumes $substrings$ to be sorted by startA, and
	 * $ReadA.sortedAlignments[firstAlignment..lastAlignment]$ to contain all alignments
	 * that are neither implied by a dense nor by a periodic substring, sorted by startA.
	 */
	public static final void splitSubstrings(int firstAlignment, int lastAlignment) {
		final double CONSTANT_THRESHOLD = 0.1;
		final int SPLIT_TOLERANCE = IO.quantum<<1;
		final double MIN_HIGH = 0.4;
		final int MAX_INTERIOR_DISTANCE = IO.quantum<<2;
		final int MAX_PEAK_RADIUS = Alignments.minAlignmentLength>>1;  // Arbitrary
		final int MAX_WIDTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		final int MIN_SPACING_FOR_HISTOGRAM = IO.quantum<<2;  // Arbitrary
		final int MIN_MASS_HIGH = 50*IO.coverage;  // Arbitrary
		int i, j;
		int from, to, lastLeftSplit, lastRightSplit, firstAlignmentForNextSubstring;
		int startA, endA, interiorStart, interiorEnd, interiorStartPrime, interiorEndPrime;
		int rangeFirst, rangeLast;
		DenseSubstring tmpSubstring;
		
		// Marking substrings that start/end inside a short-period range.
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].leftAlignment=0;  // Used as a flag
			substrings[i].rightAlignment=0;  // Used as a flag
		}
		if (PeriodicSubstrings.lastInterval>=0) {
			i=0;
			while (PeriodicSubstrings.intervals[i].hasLongPeriod) i++;
			if (i<=PeriodicSubstrings.lastInterval) {
				rangeFirst=PeriodicSubstrings.intervals[i].firstPosition;
				rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
				j=0;
				for (i=i+1; i<=PeriodicSubstrings.lastInterval; i++) {
					if (PeriodicSubstrings.intervals[i].hasLongPeriod) continue;
					if (PeriodicSubstrings.intervals[i].firstPosition<=rangeLast+identityThreshold) {
						rangeLast=Math.max(rangeLast,PeriodicSubstrings.intervals[i].lastPosition);
						continue;
					}
					j=addEvents_markSubstrings(j,rangeFirst,rangeLast,true);
					rangeFirst=PeriodicSubstrings.intervals[i].firstPosition;
					rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
				}
				addEvents_markSubstrings(j,rangeFirst,rangeLast,true);
			}
		}

		// Splitting
		firstAlignmentForNextSubstring=firstAlignment;
		for (i=0; i<=lastSubstring; i++) {
			if (!substrings[i].substringReplication) continue;
			startA=substrings[i].startA; endA=substrings[i].endA;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> substring ["+startA+".."+endA+"]:");			
			// Deciding $interiorStart$ and $interiorEnd$
			// Left histogram
			interiorStart=startA;
			lastLeftSplit=-1;
			ReadA.getDeltaHistogram(startA,endA,startA,endA,firstAlignmentForNextSubstring,lastAlignment,i<lastSubstring?substrings[i+1].startA:Math.POSITIVE_INFINITY,true,deltaLeft,false);
			cloneDeltaPoints(true);
			if (ReadA.deltaTmp[0]>=(endA-MIN_SPACING_FOR_HISTOGRAM-startA+1)/MIN_SPACING_FOR_HISTOGRAM) {
				from=Histograms.firstNonzero(deltaLeft,0,endA-startA);
				if (from!=-1) {
					to=Histograms.lastNonzero(deltaLeft,0,endA-startA);
					lastLeftSplit=Histograms.getPeaks(deltaLeft,from,to,leftSplits,0,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,tmpSplits,-1,MAX_WIDTH,deltaPoints_backupLeft,lastDeltaPoint_backupLeft,MIN_MASS_HIGH,Alignments.minAlignmentLength);
					if ( lastLeftSplit>=0 && 
						 leftSplits[1]-startA<=SPLIT_TOLERANCE && 
					     containsAlignment(leftSplits[1],substrings[i].endA,firstAlignment,lastAlignment) &&
						 Reads.isLeftMaximal(substrings[i].startA,ReadA.id,true)
				       ) {
						//substrings[i].startA=leftSplits[1];  // Commented since it is not useful in practice
						if (leftSplits[2]<endA-SPLIT_TOLERANCE && leftSplits[2]<=leftSplits[1]+MAX_INTERIOR_DISTANCE) interiorStart=leftSplits[2];
						else interiorStart=leftSplits[1];	
					}
				}
			}
			// Right histogram
			interiorEnd=endA;
			lastRightSplit=-1;
			ReadA.getDeltaHistogram(startA,endA,startA,endA,firstAlignmentForNextSubstring,lastAlignment,i<lastSubstring?substrings[i+1].startA:Math.POSITIVE_INFINITY,false,deltaRight,false);
			cloneDeltaPoints(false);
			firstAlignmentForNextSubstring=ReadA.deltaTmp[1];
			if (ReadA.deltaTmp[0]>=(endA-MIN_SPACING_FOR_HISTOGRAM-startA+1)/MIN_SPACING_FOR_HISTOGRAM) {
				from=Histograms.firstNonzero(deltaRight,0,endA-startA);
				if (from!=-1) {
					to=Histograms.lastNonzero(deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> calling Histograms.getPeaks() with first="+from+" last="+to+" distanceThreshold="+distanceThreshold);
					lastRightSplit=Histograms.getPeaks(deltaRight,from,to,rightSplits,0,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,tmpSplits,-1,MAX_WIDTH,deltaPoints_backupRight,lastDeltaPoint_backupRight,MIN_MASS_HIGH,Alignments.minAlignmentLength);
					if ( lastRightSplit>=0 && 
					     endA-rightSplits[lastRightSplit-1]<=SPLIT_TOLERANCE && 
					     containsAlignment(substrings[i].startA,rightSplits[lastRightSplit-1],firstAlignment,lastAlignment) &&
						 Reads.isRightMaximal(substrings[i].endA,ReadA.id,true)
				       ) {
						//substrings[i].endA=rightSplits[lastRightSplit-1];  // Commented since it is not useful in practice
						if (rightSplits[lastRightSplit-2]>startA+SPLIT_TOLERANCE && rightSplits[lastRightSplit-2]>=rightSplits[lastRightSplit-1]-MAX_INTERIOR_DISTANCE) interiorEnd=rightSplits[lastRightSplit-2];
						else interiorEnd=rightSplits[lastRightSplit-1];
					}
				}
			}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> interiorStart="+interiorStart+" interiorEnd="+interiorEnd);

			// Splitting the interior
			if (interiorStart==startA && interiorEnd==endA) {
				for (j=1; j<=lastLeftSplit-1; j+=3) {
					Factors.lastSplit++;
					Factors.splits[Factors.lastSplit].clear();
					Factors.splits[Factors.lastSplit].position=leftSplits[j];
					Factors.splits[Factors.lastSplit].nOpen=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> (1) added left split "+leftSplits[j]);
				}
				for (j=1; j<=lastRightSplit-1; j+=3) {
					Factors.lastSplit++;
					Factors.splits[Factors.lastSplit].clear();
					Factors.splits[Factors.lastSplit].position=rightSplits[j];
					Factors.splits[Factors.lastSplit].nClosed=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> (1) added right split "+rightSplits[j]);
				}
			}
			else {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("interior: ["+interiorStart+".."+interiorEnd+"]");				
				// Recomputing the peaks of the left histogram in the interior
				interiorStartPrime=interiorStart;
				from=Histograms.firstNonzero(deltaLeft,interiorStart-startA,interiorEnd-startA);
				if (from!=-1) {
					to=Histograms.lastNonzero(deltaLeft,interiorStart-startA,interiorEnd-startA);
					lastLeftSplit=Histograms.getPeaks(deltaLeft,from,to,leftSplits,0,startA,IO.maxDeltaStd,interiorStart==startA?distanceThreshold:0,interiorEnd==endA?distanceThreshold:0,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,tmpSplits,-1,MAX_WIDTH,deltaPoints_backupLeft,lastDeltaPoint_backupLeft,MIN_MASS_HIGH,Alignments.minAlignmentLength);
					if ( interiorStart==startA && 
					     lastLeftSplit>=0 && 
					     leftSplits[1]-startA<=SPLIT_TOLERANCE && 
					     containsAlignment(leftSplits[1],substrings[i].endA,firstAlignment,lastAlignment) &&
						 Reads.isLeftMaximal(substrings[i].startA,ReadA.id,true)
				       ) {
						//substrings[i].startA=leftSplits[1];  // Commented since it is not useful in practice
						if (leftSplits[2]<interiorEnd-SPLIT_TOLERANCE && leftSplits[2]<=leftSplits[1]+MAX_INTERIOR_DISTANCE) interiorStartPrime=leftSplits[2];
						else interiorStartPrime=leftSplits[1];
					}
				}
				// Recomputing the peaks of the right histogram in the interior
				interiorEndPrime=interiorEnd;
				from=Histograms.firstNonzero(deltaRight,interiorStart-startA,interiorEnd-startA);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> recomputing peaks of right histogram in the interior:");
				if (from!=-1) {
					to=Histograms.lastNonzero(deltaRight,interiorStart-startA,interiorEnd-startA);
					lastRightSplit=Histograms.getPeaks(deltaRight,from,to,rightSplits,0,startA,IO.maxDeltaStd,interiorStart==startA?distanceThreshold:0,interiorEnd==endA?distanceThreshold:0,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,tmpSplits,-1,MAX_WIDTH,deltaPoints_backupRight,lastDeltaPoint_backupRight,MIN_MASS_HIGH,Alignments.minAlignmentLength);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> from="+from+" to="+to);
					if ( interiorEnd==endA && 
					     lastRightSplit>=0 && 
					     endA-rightSplits[lastRightSplit-1]<=SPLIT_TOLERANCE && 
					     containsAlignment(substrings[i].startA,rightSplits[lastRightSplit-1],firstAlignment,lastAlignment) &&
						 Reads.isRightMaximal(substrings[i].endA,ReadA.id,true)
					   ) {
						//substrings[i].endA=rightSplits[lastRightSplit-1];  // Commented since it is not useful in practice
						if (rightSplits[lastRightSplit-2]>interiorStart+SPLIT_TOLERANCE && rightSplits[lastRightSplit-2]>=rightSplits[lastRightSplit-1]-MAX_INTERIOR_DISTANCE) interiorEndPrime=rightSplits[lastRightSplit-2];
						else interiorEndPrime=rightSplits[lastRightSplit-1];
					}
				}
				if (interiorStartPrime==interiorStart && interiorEndPrime==interiorEnd) {
					for (j=1; j<=lastLeftSplit-1; j+=3) {
						Factors.lastSplit++;
						Factors.splits[Factors.lastSplit].clear();
						Factors.splits[Factors.lastSplit].position=leftSplits[j];
						Factors.splits[Factors.lastSplit].nOpen=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> (2) added left split "+leftSplits[j]);						
					}
					for (j=1; j<=lastRightSplit-1; j+=3) {
						Factors.lastSplit++;
						Factors.splits[Factors.lastSplit].clear();
						Factors.splits[Factors.lastSplit].position=rightSplits[j];
						Factors.splits[Factors.lastSplit].nClosed=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> (2) added right split "+rightSplits[j]);						
					}
				}
				else {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("interiorPrime: ["+interiorStartPrime+".."+interiorEndPrime+"]");
					// Recomputing the peaks of the left histogram in the new interior
					from=Histograms.firstNonzero(deltaLeft,interiorStartPrime-startA,interiorEndPrime-startA);
					if (from!=-1) {
						to=Histograms.lastNonzero(deltaLeft,interiorStartPrime-startA,interiorEndPrime-startA);
						lastLeftSplit=Histograms.getPeaks(deltaLeft,from,to,leftSplits,0,startA,IO.maxDeltaStd,interiorStartPrime==startA?distanceThreshold:0,interiorEndPrime==endA?distanceThreshold:0,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,tmpSplits,-1,MAX_WIDTH,deltaPoints_backupLeft,lastDeltaPoint_backupLeft,MIN_MASS_HIGH,Alignments.minAlignmentLength);
						if (lastLeftSplit>=0) {
							for (j=1; j<=lastLeftSplit-1; j+=3) {
								if (leftSplits[j]-startA<=SPLIT_TOLERANCE || endA-leftSplits[j]<=SPLIT_TOLERANCE) continue;
								Factors.lastSplit++;
								Factors.splits[Factors.lastSplit].clear();
								Factors.splits[Factors.lastSplit].position=leftSplits[j];
								Factors.splits[Factors.lastSplit].nOpen=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> (3) added left split "+leftSplits[j]);
							}
						}
					}
					// Recomputing the peaks of the right histogram in the new interior
					from=Histograms.firstNonzero(deltaRight,interiorStartPrime-startA,interiorEndPrime-startA);
					if (from!=-1) {
						to=Histograms.lastNonzero(deltaRight,interiorStartPrime-startA,interiorEndPrime-startA);
						lastRightSplit=Histograms.getPeaks(deltaRight,from,to,rightSplits,0,startA,IO.maxDeltaStd,interiorStartPrime==startA?distanceThreshold:0,interiorEndPrime==endA?distanceThreshold:0,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,tmpSplits,-1,MAX_WIDTH,deltaPoints_backupRight,lastDeltaPoint_backupRight,MIN_MASS_HIGH,Alignments.minAlignmentLength);
						if (lastRightSplit>=0) {
							for (j=1; j<=lastRightSplit-1; j+=3) {
								if (rightSplits[j]-startA<=SPLIT_TOLERANCE || endA-rightSplits[j]<=SPLIT_TOLERANCE) continue;
								Factors.lastSplit++;
								Factors.splits[Factors.lastSplit].clear();
								Factors.splits[Factors.lastSplit].position=rightSplits[j];
								Factors.splits[Factors.lastSplit].nClosed=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitSubstrings> (3) added right split "+rightSplits[j]);
							}
						}
					}
				}
			}
			if (substrings[i].leftAlignment==0) {
				Factors.lastSplit++;
				Factors.splits[Factors.lastSplit].clear();
				Factors.splits[Factors.lastSplit].position=substrings[i].startA;
				Factors.splits[Factors.lastSplit].nOpen=1;
			}
			if (substrings[i].rightAlignment==0) {
				Factors.lastSplit++;
				Factors.splits[Factors.lastSplit].clear();
				Factors.splits[Factors.lastSplit].position=substrings[i].endA;
				Factors.splits[Factors.lastSplit].nClosed=1;
			}
		}
	}
	
	
	/**
	 * @return TRUE iff at least one alignment in $[firstAlignment..lastAlignment]$ is
	 * approximately contained in $[start..end]$.
	 */
	private static final boolean containsAlignment(int start, int end, int firstAlignment, int lastAlignment) {
		boolean found;
		int i, startA, endA;
		
		found=false;
		for (i=firstAlignment; i<=lastAlignment; i++) {
			startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (Intervals.isApproximatelyContained(startA,endA,start,end)) {
				found=true;
				break;
			}
		}
		return found;
	}


	/**
	 * Returns true if $position$ belongs to one of the intervals in
	 * $lengthBiasIntervals$, or if $position$ is at distance at most $maxDistance$ from
	 * the center of an interval.
	 *
	 * Remark: the procedure assumes $lengthBiasIntervals$ to be sorted.
	 */
	public static final boolean inLengthBiasInterval(int position, int maxDistance) {
		int i, j;

		if (lastLengthBiasInterval==-1) return false;
		if (lastLengthBiasInterval==0) {
			if (lengthBiasIntervals[0].firstPosition<=position && lengthBiasIntervals[0].lastPosition>=position) return true;
			if (Math.abs(position-lengthBiasIntervals[0].center)<=maxDistance) return true;
			return false;
		}
		tmpInterval.firstPosition=position;
		i=Arrays.binarySearch(lengthBiasIntervals,0,lastLengthBiasInterval+1,tmpInterval);
		if (i>=0) return true;
		i=-i-1;  // Leftmost interval whose $firstPosition$ is greater than $position$
		if (i<=lastLengthBiasInterval && lengthBiasIntervals[i].center-position<=maxDistance) return true;
		j=i-1;
		while (j>=0) {
			if (lengthBiasIntervals[j].lastPosition>=position) return true;
			if (position-lengthBiasIntervals[j].center<=maxDistance) return true;
			j--;
		}
		j=i+1;
		while (j<=lastLengthBiasInterval) {
			if (lengthBiasIntervals[j].firstPosition-position>maxDistance) break;
			if (lengthBiasIntervals[j].center-position<=maxDistance) return true;
			j++;
		}
		return false;
	}


	/**
	 * Tries to decide whether the prefixes of prefix substring (respectively, the
	 * suffixes of suffix substring) $substring$ are equally spaced by a period. If so,
	 * the procedure removes from $lengthBoundaries$ all candidate boundaries that are too
	 * close to a multiple of the period.
	 *
	 * Remark: the procedure also sets values $period$, $*AreUniformlyDistributed$,
	 * $min*DistanceBias$ and $max*DistanceBias$ of $substring$.
	 *
	 * @param offset amount subtracted from every element in $lengthBoundaries$ when 
	 * comparing it to the period;
	 * @return the last element in $lengthBoundaries$ after removal.
	 */
	private static final int removeBoundariesAtPeriods(int[] lengthBoundaries, int lastBoundary, DenseSubstring substring, boolean prefix, int firstAlignment, int lastAlignment, int offset) {
		final int PERIOD_THRESHOLD = MIN_SHIFT_DIFFERENCE>>2;
		boolean found;
		int i, j, k;
		int from, to;
		final int lastPoint, period;

		lastPoint=getShiftsPrefixSuffix(substring,firstAlignment,lastAlignment,prefix);
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("removeBoundariesAtPeriods> shifts of dense substring ["+substring.startA+".."+substring.endA+"]");
	for (i=0; i<=lastPoint; i++) IO.printErr(lengthPoints[i].position+","+lengthPoints[i].mass); 
}

		period=PeriodicSubstrings.estimatePeriod(lengthPoints,lastPoint,PeriodicSubstrings.MIN_PERIOD,PeriodicSubstrings.minIntervalLength,PeriodicSubstrings.minLocalMaxDistance);  // Initializes array $PeriodicSubstrings.periodTmp$
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("removeBoundariesAtPeriods> PERIOD="+period);		
		substring.period=period;
		substring.strongPeriodSignal=period>0&&PeriodicSubstrings.periodTmp[4]>1;
		if (lastBoundary>=0 && period>0) {
			k=-1;
			for (i=1; i<=lastBoundary; i+=3) {
				found=false;
				from=Math.max(lengthBoundaries[i]-PERIOD_THRESHOLD,lengthBoundaries[i-1]);
				to=Math.min(lengthBoundaries[i]+PERIOD_THRESHOLD,lengthBoundaries[i+1]);
				for (j=from; j<=to; j++) {
					if ((j-offset)%period==0) {
						found=true;
						break;
					}
				}
				if (found) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("removeBoundariesAtPeriods> discarded boundary ["+lengthBoundaries[i-1]+".."+lengthBoundaries[i]+".."+lengthBoundaries[i+1]+"]");		
					continue;
				}
				lengthBoundaries[++k]=lengthBoundaries[i-1];
				lengthBoundaries[++k]=lengthBoundaries[i];
				lengthBoundaries[++k]=lengthBoundaries[i+1];
			}
		}
		else k=lastBoundary;
		if (prefix) {
			if (PeriodicSubstrings.periodTmp[0]>1) {
				substring.prefixesAreUniformlyDistributed=false;
				substring.minPrefixDistanceBias=PeriodicSubstrings.periodTmp[1];
				substring.maxPrefixDistanceBias=PeriodicSubstrings.periodTmp[2];
			}
			else {
				substring.prefixesAreUniformlyDistributed=true;
				substring.minPrefixDistanceBias=-1;
				substring.maxPrefixDistanceBias=-1;
			}
		}
		else {
			if (PeriodicSubstrings.periodTmp[0]>1) {
				substring.suffixesAreUniformlyDistributed=false;
				substring.minSuffixDistanceBias=PeriodicSubstrings.periodTmp[1];
				substring.maxSuffixDistanceBias=PeriodicSubstrings.periodTmp[2];
			}
			else {
				substring.suffixesAreUniformlyDistributed=true;
				substring.minSuffixDistanceBias=-1;
				substring.maxSuffixDistanceBias=-1;
			}
		}
		return k;
	}

	
	/**
	 * Stores in $lengthPoints$ the sorted and compacted list of all differences between 
	 * consecutive endpoints of all alignments implied by $substring$, of prefix (if 
	 * $prefix=true$) or suffix (if $prefix=false$) type, sorted by endpoint.
	 *
	 * @param first, last an interval of $ReadA.sortedAlignments$ that contains all the
	 * alignments implied by $substring$;
	 * @return the last point in $lengthPoints$.
	 */
	private static final int getShiftsPrefixSuffix(DenseSubstring substring, int first, int last, boolean prefix) {
		final int THRESHOLD = IO.quantum;
		int i;
		int tmp, lastPoint;
		Alignment alignment;
		
		lastTmpAlignment=-1;
		for (i=first; i<=last; i++) {
			alignment=ReadA.sortedAlignments[i];
			if (alignment.inPeriodicSubstring) continue;
			if ( ( prefix && isAlignmentImpliedByPrefix(i,substring,false) && 
				   ( alignment.isRightMaximal==1 ||
					 ( alignment.isRightMaximalB==1 &&
					   Math.abs(Alignments.alignments[alignment.id][4],substring.endA)<=THRESHOLD
					 )
				   )
				 ) ||
				 ( !prefix && isAlignmentImpliedBySuffix(i,substring,false) &&
  				   ( alignment.isLeftMaximal==1 || 
					 ( alignment.isLeftMaximalB==1 &&
  					   Math.abs(Alignments.alignments[alignment.id][3],substring.startA)<=THRESHOLD
					 )
  				   )
				 ) 
			   ) tmpAlignments[++lastTmpAlignment]=alignment;
		}
		tmp=Alignment.order;
		Alignment.order=prefix?Alignment.ENDA:Alignment.STARTA_ENDA;
		if (lastTmpAlignment>0) Arrays.sort(tmpAlignments,0,lastTmpAlignment+1);
		Alignment.order=tmp;
		lastPoint=-1;
		for (i=1; i<=lastTmpAlignment; i++) {
			lastPoint++;
			lengthPoints[lastPoint].position=Alignments.alignments[tmpAlignments[i].id][prefix?4:3]-Alignments.alignments[tmpAlignments[i-1].id][prefix?4:3];
			lengthPoints[lastPoint].mass=1;
		}
		lastPoint=Points.sortAndCompact(lengthPoints,lastPoint);
		
if (IO.SHOW_STD_ERR) {
	IO.printErr("SHIFTS FOR DENSE SUBSTRING OF PREFIX/SUFFIX TYPE ["+substring.startA+".."+substring.endA+"]:");		
	for (int x=0; x<=lastPoint; x++) IO.printErr(lengthPoints[x]);
}


		return lastPoint;
	}
	
	
	/**
	 * Stores in $lengthPoints$ the sorted and compacted set of shifts between all pairs 
	 * of distinct B-maximal events inside the interval $[start..end]$ of a dense 
	 * substring of substring type. The procedure does not use alignments implied by 
	 * another dense substring (since such substring is necessarily straddling) or by a 
	 * periodic substring.
	 *
	 * @return the last point in $lengthPoints$. $shiftTmp[0]$: number of alignments used; 
	 * $shiftTmp[1]$: the first alignment from $firstAlignment$ that either covers $next$,
	 * or starts to the right of $end$.
	 */
	private static final int getShiftsSubstring(int start, int end, int firstAlignment, int lastAlignment, int next) {
		boolean isLeftMaximalB, isRightMaximalB;
		int i, j, k, startA, endA, nUsedAlignments, firstAlignmentForNext;
		int lastPoint;
		Point[] tmpPoints;

		// Filling $lengthPoints$ with positions
		nUsedAlignments=0; firstAlignmentForNext=-1; lastPoint=-1;
		for (i=firstAlignment; i<=lastAlignment; i++) {
			endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (firstAlignmentForNext==-1 && endA>=next) firstAlignmentForNext=i;
			if (endA<start || endA>end) continue;
			if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null || ReadA.sortedAlignments[i].impliedByDenseSubstring!=null) continue;
			startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			if (startA>end) break;
			if (startA<start) continue;
			nUsedAlignments++;
			if (ReadA.sortedAlignments[i].isLeftMaximalB==1) {
				lastPoint++;
				if (lastPoint==lengthPoints.length) {
					tmpPoints = new Point[lengthPoints.length+GROWTH_RATE];
					System.arraycopy(lengthPoints,0,tmpPoints,0,lengthPoints.length);
					for (k=lengthPoints.length; k<tmpPoints.length; k++) tmpPoints[k] = new Point();
					lengthPoints=tmpPoints;
				}
				lengthPoints[lastPoint].position=startA-start;
				lengthPoints[lastPoint].mass=1;
			}
			if (ReadA.sortedAlignments[i].isRightMaximalB==1) {
				lastPoint++;
				if (lastPoint==lengthPoints.length) {
					tmpPoints = new Point[lengthPoints.length+GROWTH_RATE];
					System.arraycopy(lengthPoints,0,tmpPoints,0,lengthPoints.length);
					for (k=lengthPoints.length; k<tmpPoints.length; k++) tmpPoints[k] = new Point();
					lengthPoints=tmpPoints;
				}
				lengthPoints[lastPoint].position=endA-start;
				lengthPoints[lastPoint].mass=1;
			}
		}
		shiftTmp[0]=nUsedAlignments;
		shiftTmp[1]=firstAlignmentForNext!=-1?firstAlignmentForNext:i;
		
		// Filling $lengthPoints$ with shifts
		lastPoint=Points.sortAndCompact(lengthPoints,lastPoint);
		for (i=1; i<=lastPoint; i++) {
			lengthPoints[i-1].position=lengthPoints[i].position-lengthPoints[i-1].position;
			lengthPoints[i-1].mass=1;
		}
		lastPoint=Points.sortAndCompact(lengthPoints,lastPoint-1);
		return lastPoint;
	}
	
	
	/**
	 * Stores in $lengthPoints$ the sorted and compacted set of lengths of all alignments 
	 * that occur fully inside the dense substring of substring type $[start..end]$. 
	 * The procedure does not use alignments implied by another dense substring (since 
	 * such substring is necessarily straddling) or by a periodic substring, nor 
	 * alignments that are not B-left-maximal or not B-right-maximal.
	 *
	 * @return the last point in $lengthPoints$. $shiftTmp[0]$: number of alignments used; 
	 * $shiftTmp[1]$: the first alignment from $firstAlignment$ that either covers $next$,
	 * or starts to the right of $end$.
	 */
	private static final int getLengthsSubstring(int start, int end, int firstAlignment, int lastAlignment, int next) {
		boolean isLeftMaximalB, isRightMaximalB;
		int i, j, startA, endA, nUsedAlignments, firstAlignmentForNext;
		int lastPoint;

		// Filling $lengthPoints$ with lengths
		nUsedAlignments=0; firstAlignmentForNext=-1; lastPoint=-1;
		for (i=firstAlignment; i<=lastAlignment; i++) {
			endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (firstAlignmentForNext==-1 && endA>=next) firstAlignmentForNext=i;
			if (endA<start || endA>end) continue;
			if ( ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null || 
				 ReadA.sortedAlignments[i].impliedByDenseSubstring!=null ||
				 ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
				 ReadA.sortedAlignments[i].isRightMaximalB!=1
			   ) continue;
			startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			if (startA>end) break;
			if (startA<start) continue;
			nUsedAlignments++;
			lastPoint++;
			lengthPoints[lastPoint].position=endA-startA+1;
			lengthPoints[lastPoint].mass=1;
		}
		shiftTmp[0]=nUsedAlignments;
		shiftTmp[1]=firstAlignmentForNext!=-1?firstAlignmentForNext:i;
		lastPoint=Points.sortAndCompact(lengthPoints,lastPoint);		
		return lastPoint;
	}
	
	


	/**
	 * (Not called anywhere)
	 * (Superseded by using ratios rather than values in the definition of high maxima)
	 *
	 * Removes from $lengthBoundaries$ all candidate boundaries whose density is too small
	 * compared to the density of the non-candidate region.
	 *
	 * @return the last element in $lengthBoundaries$ after removal.
	 */
	private static final int removeBoundariesWithLowDensity(int[] lengthBoundaries, int lastBoundary, double[] delta, int maxLength, int origin) {
		int i, j, k;
		double densityNoncandidate, densityCandidate, lengthNoncandidate, previous;

		densityNoncandidate=0.0; lengthNoncandidate=0.0;
		j=0; previous=-1;
		for (i=0; i<=maxLength-origin; i++) {
			if (j>lastBoundary || i<lengthBoundaries[j]) densityNoncandidate+=delta[i];
			else if (i==lengthBoundaries[j]) lengthNoncandidate+=i-previous-1;
			else if (i==lengthBoundaries[j+2]) {
				previous=i;
				j+=3;
			}
		}
		lengthNoncandidate+=maxLength-origin-previous;
		if (lengthNoncandidate==0.0) return lastBoundary;
		densityNoncandidate/=lengthNoncandidate;		
		k=-1;
		for (i=0; i<=lastBoundary; i+=3) {
			densityCandidate=0.0;
			for (j=lengthBoundaries[i]; j<=lengthBoundaries[i+2]; j++) densityCandidate+=delta[j];
			densityCandidate/=lengthBoundaries[i+2]-lengthBoundaries[i]+1.0;
			if (densityCandidate>=densityNoncandidate*MIN_LENGTH_BOUNDARY_DENSITY_RATIO) {
				lengthBoundaries[++k]=lengthBoundaries[i];
				lengthBoundaries[++k]=lengthBoundaries[i+1];
				lengthBoundaries[++k]=lengthBoundaries[i+2];
			}
		}
		return k;
	}
	
	
	/**
	 * Searches for a maximal dense substring, of weak prefix type, that starts from 
	 * alignment $ReadA.sortedAlignments[source]$. This is analogous to procedure
	 * $getPrefixSubstring$.
	 *		
	 * @return $tmp[0]$: the position in $ReadA.sortedAlignments[0..nVertices-1]$ of an
	 * interval that maximizes the distance from $source$ in the prefix DAG, and that ends
	 * the farthest from $source$, or -1 if no valid path from $source$ is found;
	 * $tmp[1]$: the length of a longest path from $source$ in the prefix DAG;
	 * $tmp[2]$: the rightmost end of a B-right-maximal alignment in the chosen longest 
	 * path;
	 * $tmp[3],tmp[4]$: the smallest (largest) start of a B-left-maximal alignment in the
	 * chosen longest path, or -1 if no such event can be found;
	 * $tmp[5]$: the position in $ReadA.sortedAlignments[0..nVertices-1]$ of a 
	 * B-left-maximal alignment, in the chosen path, that starts at $tmp[3]$ in readA;
	 * $tmp[6]$: the smallest start of an alignment in the chosen longest path.
	 */
	private static final void getPrefixSubstring_weak(int source, int nVertices, int minPathLength) {
		int i, j, k, s;
		int sourceComponent, vertex, rightmostAlignment, distance, longestDistance;
		int start, end, rightmostEnd, rightmostRightMaximalEnd, length;
		int lastSortedVertex;
		int minBLeft, minBLeftMaximal, maxBLeftMaximal, minBLeftMaximalAlignment;

		// Computing the longest distance from $source$ to every other vertex
		Math.set(distances,nVertices-1,-1);
		Math.set(predecessors,nVertices-1,-1);
		sourceComponent=components[source];
		distances[source]=0;
		s=original2sorted[source];
 		lastSortedVertex=s;
		for (j=s+1; j<nVertices; j++) {
			vertex=sorted2original[j];
			if (components[vertex]!=sourceComponent) {
				lastSortedVertex=j-1;
				break;
			}
			for (i=0; i<nInNeighbors[vertex]; i++) {
				k=inNeighbors[vertex][i];
				if (distances[k]!=-1) distance=distances[k];
				else continue;
				distance++;
				if (distance>distances[vertex]) {
					distances[vertex]=distance;
					predecessors[vertex]=k;
				}
			}
		}
		if (j==nVertices) lastSortedVertex=j-1;

		// Computing the rightmost alignment whose distance from $source$ in the DAG is
		// maximum.
		longestDistance=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]>longestDistance) longestDistance=distances[vertex];
		}
		if (longestDistance+1<minPathLength) {  // No path from $source$ with the required structure
			Math.set(tmp,6,-1);
			return;
		}
		rightmostEnd=-1; rightmostAlignment=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]<longestDistance) continue;
			end=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if (end>rightmostEnd) {
				rightmostAlignment=vertex;
				rightmostEnd=end;
			}
		}
		length=rightmostEnd-Alignments.alignments[ReadA.sortedAlignments[source].id][4];
		if (length<=prefixSuffixThreshold) {
			Math.set(tmp,6,-1);
			return;
		}
		
		// Computing statistics on the alignments in the longest path from $source$
		rightmostRightMaximalEnd=-1;
		minBLeft=Math.POSITIVE_INFINITY;
		minBLeftMaximal=Math.POSITIVE_INFINITY; maxBLeftMaximal=-1; 
		minBLeftMaximalAlignment=-1;
		vertex=rightmostAlignment;
		while (vertex!=-1) {
			end=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if (ReadA.sortedAlignments[vertex].isRightMaximalB==1 && !ReadA.sortedAlignments[vertex].lowQualityEnd && end>rightmostRightMaximalEnd) rightmostRightMaximalEnd=end;
			start=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			if (start<minBLeft) minBLeft=start;
			if (ReadA.sortedAlignments[vertex].isLeftMaximalB==1 && !ReadA.sortedAlignments[vertex].lowQualityStart) {
				if (start<minBLeftMaximal) {
					minBLeftMaximal=start;
					minBLeftMaximalAlignment=vertex;
				}
				if (start>maxBLeftMaximal) maxBLeftMaximal=start;
			}
			vertex=predecessors[vertex];
		}
		tmp[0]=rightmostAlignment;
		tmp[1]=longestDistance;
		tmp[2]=rightmostRightMaximalEnd;
		if (maxBLeftMaximal==-1) {
			tmp[3]=-1;
			tmp[4]=-1;
			tmp[5]=-1;
		}
		else {
			tmp[3]=minBLeftMaximal;
			tmp[4]=maxBLeftMaximal;
			tmp[5]=minBLeftMaximalAlignment;
		}
		tmp[6]=minBLeft;
	}
	
	
	/**
	 * Searches for a maximal dense substring, of weak suffix type, that starts from 
	 * alignment $ReadA.sortedAlignments[source]$. This is analogous to procedure
	 * $getSuffixSubstring$.
	 *		
	 * @return $tmp[0]$: the position in $ReadA.sortedAlignments[0..nVertices-1]$ of an
	 * interval that maximizes the distance from $source$ in the prefix DAG, and that ends
	 * the farthest from $source$, or -1 if no valid path from $source$ is found;
	 * $tmp[1]$: the length of a longest path from $source$ in the prefix DAG;
	 * $tmp[2]$: the leftmost end of a B-left-maximal alignment in the chosen longest 
	 * path;
	 * $tmp[3],tmp[4]$: the smallest (largest) end of a B-right-maximal alignment in the
	 * chosen longest path, or -1 if no such event can be found;
	 * $tmp[5]$: the position in $ReadA.sortedAlignments[0..nVertices-1]$ of a 
	 * B-right-maximal alignment, in the chosen path, that ends at $tmp[4]$ in readA;
	 * $tmp[6]$: the largest end of an alignment in the chosen longest path.
	 */
	private static final void getSuffixSubstring_weak(int source, int nVertices, int minPathLength) {
		int i, j, k, s;
		int sourceComponent, vertex, rightmostAlignment, distance, longestDistance;
		int start, end, rightmostStart, leftmostLeftMaximalStart, length;
		int lastSortedVertex;
		int maxBRight, minBRightMaximal, maxBRightMaximal, maxBRightMaximalAlignment;

		// Computing the longest distance from $source$ to every other vertex
		Math.set(distances,nVertices-1,-1);
		Math.set(predecessors,nVertices-1,-1);
		sourceComponent=components[source];
		distances[source]=0;
		s=original2sorted[source];
 		lastSortedVertex=s;
		for (j=s+1; j<nVertices; j++) {
			vertex=sorted2original[j];
			if (components[vertex]!=sourceComponent) {
				lastSortedVertex=j-1;
				break;
			}
			for (i=0; i<nInNeighbors[vertex]; i++) {
				k=inNeighbors[vertex][i];
				if (distances[k]!=-1) distance=distances[k];
				else continue;
				distance++;
				if (distance>distances[vertex]) {
					distances[vertex]=distance;
					predecessors[vertex]=k;
				}
			}
		}
		if (j==nVertices) lastSortedVertex=j-1;

		// Computing the rightmost alignment whose distance from $source$ in the DAG is
		// maximum.
		longestDistance=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]>longestDistance) longestDistance=distances[vertex];
		}
		if (longestDistance+1<minPathLength) {  // No path from $source$ with the required structure
			Math.set(tmp,6,-1);
if (IO.SHOW_STD_ERR) IO.printErr("path length not enough = "+(longestDistance+1)+" < "+minPathLength);			
			return;
		}
		rightmostStart=-1; rightmostAlignment=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]<longestDistance) continue;
			start=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			if (start>rightmostStart) {
				rightmostAlignment=vertex;
				rightmostStart=start;
			}
		}
		length=rightmostStart-Alignments.alignments[ReadA.sortedAlignments[source].id][3];
		if (length<=prefixSuffixThreshold) {
			Math.set(tmp,6,-1);
			return;
		}
		
		// Computing statistics on the alignments in the longest path from $source$
		leftmostLeftMaximalStart=Math.POSITIVE_INFINITY;
		maxBRight=-1;
		minBRightMaximal=Math.POSITIVE_INFINITY; maxBRightMaximal=-1; 
		maxBRightMaximalAlignment=-1;
		vertex=rightmostAlignment;
		while (vertex!=-1) {
			start=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			if (ReadA.sortedAlignments[vertex].isLeftMaximalB==1 && !ReadA.sortedAlignments[vertex].lowQualityStart && start<leftmostLeftMaximalStart) leftmostLeftMaximalStart=start;
			end=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if (end>maxBRight) maxBRight=end;
			if (ReadA.sortedAlignments[vertex].isRightMaximalB==1 && !ReadA.sortedAlignments[vertex].lowQualityEnd) {
				if (end<minBRightMaximal) minBRightMaximal=end;
				if (end>maxBRightMaximal) {
					maxBRightMaximal=end;
					maxBRightMaximalAlignment=vertex;
				}
			}
			vertex=predecessors[vertex];
		}
		tmp[0]=rightmostAlignment;
		tmp[1]=longestDistance;
		tmp[2]=leftmostLeftMaximalStart;
		if (maxBRightMaximal==-1) {
			tmp[3]=-1;
			tmp[4]=-1;
			tmp[5]=-1;
		}
		else {
			tmp[3]=minBRightMaximal;
			tmp[4]=maxBRightMaximal;
			tmp[5]=maxBRightMaximalAlignment;
		}
		tmp[6]=maxBRight;
	}
	
	
	/**
	 * Remark: the procedure updates $sumStartA,nStartA$ of prefix substrings and 
	 * $sumEndA,nEndA$ of suffix substrings, and it refines the boundaries of weak 
	 * substrings that replicate by taking substrings of themselves: for details see 
	 * procedure $markImpliedAlignmentsAggressive_weak$.
	 *
	 * @param prefixDirection TRUE iff $ReadA.sortedAlignments$ is sorted as while finding
	 * weak prefix substrings;
	 * @return the number of alignments marked by the procedure.
	 */
	private static final int markImpliedAlignmentsPath_weak(DenseSubstring substring, int lastAlignment, boolean prefixDirection) {
		int vertex, length, out;

		// Marking alignments in the longest path from source to destination
		out=0;
		if (substring.substringReplication) {
			lastWeakLeft=-1;
			lastWeakRight=-1;
			substring.minPrefixLength=Math.POSITIVE_INFINITY;
		}
		vertex=substring.destinationAlignment;
		while (vertex!=-1) {
			if ( substring.prefixReplication && 
			     ReadA.sortedAlignments[vertex].inDenseSubstring==null && 
			     ReadA.sortedAlignments[vertex].impliedByPrefixSubstring==null && 
			     ReadA.sortedAlignments[vertex].impliedByPrefixSubstringPrime==null
		       ) {
				substring.nImpliedAlignmentsPrefix++;
				ReadA.sortedAlignments[vertex].impliedByPrefixSubstring=substring;
				if (ReadA.sortedAlignments[vertex].isLeftMaximalB==1) {
					substring.sumStartA+=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
					substring.nStartA++;
				}
				out++;
			}
			else if ( substring.suffixReplication && 
			          ReadA.sortedAlignments[vertex].inDenseSubstring==null && 
			          ReadA.sortedAlignments[vertex].impliedBySuffixSubstring==null && 
			          ReadA.sortedAlignments[vertex].impliedBySuffixSubstringPrime==null
		            ) {
				substring.nImpliedAlignmentsSuffix++;
				ReadA.sortedAlignments[vertex].impliedBySuffixSubstring=substring;
				if (ReadA.sortedAlignments[vertex].isRightMaximalB==1) {
					substring.sumEndA+=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
					substring.nEndA++;
				}
				out++;
			}
			else if (substring.substringReplication) {
				// Collecting B-maximal events inside $substring$, rather than marking.
				if (ReadA.sortedAlignments[vertex].isLeftMaximalB==1) {
					lastWeakLeft++;
					weakLeft[lastWeakLeft].position=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
					weakLeft[lastWeakLeft].mass=1;
				}
				if (ReadA.sortedAlignments[vertex].isRightMaximalB==1) {
					lastWeakRight++;
					weakRight[lastWeakRight].position=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
					weakRight[lastWeakRight].mass=1;
				}
				if (ReadA.sortedAlignments[vertex].isLeftMaximalB==1 && ReadA.sortedAlignments[vertex].isRightMaximalB==1) {
					length=ReadA.sortedAlignments[vertex].getALength();
					if (length<substring.minPrefixLength) substring.minPrefixLength=length;
				}
			}
			vertex=predecessors[vertex];
		}

		// Marking other alignments
		out+=markImpliedAlignmentsAggressive_weak(substring,lastAlignment,prefixDirection);
		if (substring.minPrefixLength==Math.POSITIVE_INFINITY) substring.minPrefixLength=-1;
		if (substring.minSuffixLength==Math.POSITIVE_INFINITY) substring.minSuffixLength=-1;
		return out;
	}


	/**
	 * (1) Marks as implied all the alignments that behave like those in the longest path 
	 * of a weak dense substring of prefix/suffix type. (2) Detects local maxima in the 
	 * distribution of non-implied B-left- and B-right-maximal events inside a weak dense 
	 * substring of substring type, and resets the boundaries of the substring to the 
	 * leftmost (respectively, rightmost) B-left-maximal (respectively, B-right-maximal) 
	 * run of local-maximum leaves in the DET. (3) Marks alignments inside such boundaries
	 * as belonging to a dense substring of substring type.
	 *
	 * Point (2) is useful, since a dense substring of substring type is not highly 
	 * constrained, and could include alignments that belong to nearby modules that should 
	 * be reported separately.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments$ to be sorted by endA for 
	 * prefix substrings, and $ReadA.sortedAlignments$ to be sorted by startA for suffix 
	 * substrings. It also assumes all alignments implied by periodic substrings to be 
	 * located in $ReadA.sortedAlignments[lastAlignment+1..]$.
	 *
	 * Remark: the procedure updates $sumStartA,nStartA,startA$ of a prefix substring, and
	 * $sumEndA,nEndA,endA$ of a suffix substring.
	 *
	 * Remark: like $markImpliedAlignmentsAggressive()$, the procedure tries to 
	 * estimate peaks in the density of maximal events inside a prefix/suffix region.
	 * The procedure sets $period$ as well, for dense substrings of prefix 
	 * (respectively, suffix) type only, using just the last (respectively, first) 
	 * position of all alignments implied by the substring.
	 *
	 * @param prefixDirection TRUE iff $ReadA.sortedAlignments$ is sorted as if we are 
	 * finding weak prefix substrings;
	 * @return the number of alignments marked by the procedure.
	 */
	private static final int markImpliedAlignmentsAggressive_weak(DenseSubstring substring, int lastAlignment, boolean prefixDirection) {
		final int BIN_LENGTH = Math.max((substring.endA-substring.startA+1)/10,identityThreshold<<1);
		final int GAP_THRESHOLD = (3*identityThreshold)>>1;  // Arbitrary
		final int minIntervalLength = identityThreshold;
		final int minLocalMaxDistance = minIntervalLength;
		int i;
		int from, startA, endA, first, last, lastPoint, length, out;
		int leftmostStartA, rightmostEndA;


if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markImpliedAlignmentsAggressive_weak> considering dense substring id="+substring.id+"=["+substring.startA+".."+substring.endA+"]");
	
		out=0;
		if (substring.prefixReplication) {
			substring.period=-1;
			first=markImpliedAlignmentsAggressive_weak_prefix(substring,lastAlignment);
			out+=weakTmp[2];
			last=substring.getFromAlignmentPrime(identityThreshold,lastAlignment);
			// Setting $period$ if it is not already set.
			if (substring.period==-1) {
				lastPoint=getShiftsPrefixSuffix(substring,first,last,true);
				if (lastPoint>=0 && Points.mass(lengthPoints,0,lastPoint)>=MIN_ALIGNMENTS_FOR_LENGTH_BIAS) {
					substring.period=PeriodicSubstrings.estimatePeriod(lengthPoints,lastPoint,PeriodicSubstrings.MIN_PERIOD,PeriodicSubstrings.minIntervalLength,PeriodicSubstrings.minLocalMaxDistance);
					substring.strongPeriodSignal=substring.period>0&&PeriodicSubstrings.periodTmp[4]>1;
				}
				else {
					substring.period=0;
					substring.strongPeriodSignal=false;
				}
			}
		}
		else if (substring.suffixReplication) {
			substring.period=-1;
			first=substring.getFromAlignment(identityThreshold);
			last=markImpliedAlignmentsAggressive_weak_suffix(substring,lastAlignment);
			out+=weakTmp[2];
			// Setting $period$ if it is not already set.
			if (substring.period==-1) {
				lastPoint=getShiftsPrefixSuffix(substring,first,last,false);
				if (lastPoint>=0 && Points.mass(lengthPoints,0,lastPoint)>=MIN_ALIGNMENTS_FOR_LENGTH_BIAS) {
					substring.period=PeriodicSubstrings.estimatePeriod(lengthPoints,lastPoint,PeriodicSubstrings.MIN_PERIOD,PeriodicSubstrings.minIntervalLength,PeriodicSubstrings.minLocalMaxDistance);
					substring.strongPeriodSignal=substring.period>0&&PeriodicSubstrings.periodTmp[4]>1;
				} 
				else {
					substring.period=0;
					substring.strongPeriodSignal=false;
				}
			}
		}
		else {  
			
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markImpliedAlignmentsAggressive_weak> SUBSTRING TYPE - ORIGINAL BOUNDARIES: ["+substring.startA+".."+substring.endA+"]");			
			
			// Collecting B-maximal events inside $substring$. Remark: $weakLeft$ and
			// $weakRight$ can already contain events from 
			// $markImpliedAlignmentsPath_weak$.
			if (prefixDirection) {
				from=substring.getFromAlignmentPrime(identityThreshold,lastAlignment);
				for (i=from; i>=0; i--) {
					endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					if (endA<substring.startA) break;
					if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].inPeriodicSubstring
					   ) continue;
					if (endA<=substring.endA && ReadA.sortedAlignments[i].isRightMaximalB==1) {
						lastWeakRight++;
						weakRight[lastWeakRight].position=endA;
						weakRight[lastWeakRight].mass=1;
					}
					startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					if (startA>=substring.startA && startA<=substring.endA && ReadA.sortedAlignments[i].isLeftMaximalB==1) {
						lastWeakLeft++;
						weakLeft[lastWeakLeft].position=startA;
						weakLeft[lastWeakLeft].mass=1;
					}
				}
				for (i=from+1; i<=lastAlignment; i++) {
					if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].inPeriodicSubstring
					   ) continue;
					startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					if (startA>=substring.startA && startA<=substring.endA && ReadA.sortedAlignments[i].isLeftMaximalB==1) {
						lastWeakLeft++;
						weakLeft[lastWeakLeft].position=startA;
						weakLeft[lastWeakLeft].mass=1;
					}
				}
			}
			else {
				from=substring.getFromAlignment(identityThreshold);
				for (i=from; i<=lastAlignment; i++) {
					startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					if (startA>substring.endA) break;
					if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].inPeriodicSubstring
					   ) continue;
					if (startA>=substring.startA && ReadA.sortedAlignments[i].isLeftMaximalB==1) {
						lastWeakLeft++;
						weakLeft[lastWeakLeft].position=startA;
						weakLeft[lastWeakLeft].mass=1;
					}
					endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					if (endA>=substring.startA && endA<=substring.endA && ReadA.sortedAlignments[i].isRightMaximalB==1) {
						lastWeakRight++;
						weakRight[lastWeakRight].position=endA;
						weakRight[lastWeakRight].mass=1;
					}
				}
				for (i=from-1; i>=0; i--) {
					if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
					     ReadA.sortedAlignments[i].inPeriodicSubstring
					   ) continue;
					endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					if (endA>=substring.startA && endA<=substring.endA && ReadA.sortedAlignments[i].isRightMaximalB==1) {
						lastWeakRight++;
						weakRight[lastWeakRight].position=endA;
						weakRight[lastWeakRight].mass=1;
					}
				}
			}
			
			// Trying to refine the boundaries of $substring$. Refinement is performed
			// iff a leftmost (rightmost) alignment inside the refined boundaries is not
			// too far from the left (right) refined boundary.
			Intervals.refineBoundaries(substring.startA,substring.endA,substring.startA,substring.endA,weakLeft,lastWeakLeft,weakRight,lastWeakRight,BIN_LENGTH,minIntervalLength,minLocalMaxDistance,weakTmp);
			if (weakTmp[0]!=substring.startA || weakTmp[1]!=substring.endA) {
				leftmostStartA=Math.POSITIVE_INFINITY; rightmostEndA=-1;
				if (prefixDirection) {
					from=substring.getFromAlignmentPrime(identityThreshold,lastAlignment);
					for (i=from; i>=0; i--) {
						endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
						if (endA<substring.startA) break;
						if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
							 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
							 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
							 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
							 ReadA.sortedAlignments[i].inDenseSubstring!=null ||
							 ReadA.sortedAlignments[i].inPeriodicSubstring
							) continue;
						startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
						if ( Intervals.isApproximatelyContained(startA,endA,weakTmp[0],weakTmp[1]) ||
							 Intervals.areApproximatelyIdentical(startA,endA,weakTmp[0],weakTmp[1])
						   ) {
							if (startA<leftmostStartA) leftmostStartA=startA;
							if (endA>rightmostEndA) rightmostEndA=endA;
						}
					}
					for (i=from+1; i<=lastAlignment; i++) {
						if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
							 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
							 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
							 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
							 ReadA.sortedAlignments[i].inDenseSubstring!=null ||
							 ReadA.sortedAlignments[i].inPeriodicSubstring
						   ) continue;
						startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
						endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
						if ( Intervals.isApproximatelyContained(startA,endA,weakTmp[0],weakTmp[1]) ||
							 Intervals.areApproximatelyIdentical(startA,endA,weakTmp[0],weakTmp[1])
						   ) {
   							if (startA<leftmostStartA) leftmostStartA=startA;
   							if (endA>rightmostEndA) rightmostEndA=endA;
						}
					}
				}
				else {
					from=substring.getFromAlignment(identityThreshold);
					for (i=from; i<=lastAlignment; i++) {
						startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
						if (startA>substring.endA) break;
						if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
							 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
							 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
							 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
							 ReadA.sortedAlignments[i].inDenseSubstring!=null ||
							 ReadA.sortedAlignments[i].inPeriodicSubstring
						   ) continue;
						endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
						if ( Intervals.isApproximatelyContained(startA,endA,weakTmp[0],weakTmp[1]) ||
							 Intervals.areApproximatelyIdentical(startA,endA,weakTmp[0],weakTmp[1])
						   ) {
  							if (startA<leftmostStartA) leftmostStartA=startA;
  							if (endA>rightmostEndA) rightmostEndA=endA;
						}
					}
					for (i=from-1; i>=0; i--) {
						if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||
							 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
							 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
							 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
						     ReadA.sortedAlignments[i].inDenseSubstring!=null ||
							 ReadA.sortedAlignments[i].inPeriodicSubstring
						   ) continue;
						startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
						endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
						if ( Intervals.isApproximatelyContained(startA,endA,weakTmp[0],weakTmp[1]) ||
							 Intervals.areApproximatelyIdentical(startA,endA,weakTmp[0],weakTmp[1])
						   ) {
 							if (startA<leftmostStartA) leftmostStartA=startA;
 							if (endA>rightmostEndA) rightmostEndA=endA;
						}
					}
				}
				if (leftmostStartA<=weakTmp[0]+GAP_THRESHOLD && rightmostEndA>=weakTmp[1]-GAP_THRESHOLD) {
					substring.startA=weakTmp[0];
					substring.endA=weakTmp[1];
					if (substring.maximalStartA<substring.startA) substring.maximalStartA=Math.POSITIVE_INFINITY;
					if (substring.maximalEndA>substring.endA) substring.maximalEndA=-1;
				}
			}
   			substring.sumStartA=substring.startA;
   			substring.nStartA=1;
			substring.sumEndA=substring.endA;
			substring.nEndA=1;
			
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markImpliedAlignmentsAggressive_weak> REFINED BOUNDARIES: ["+substring.startA+".."+substring.endA+"]");
			
			
			// Marking alignments
			if (prefixDirection) {
				from=substring.getFromAlignmentPrime(identityThreshold,lastAlignment);
				for (i=from; i>=0; i--) {					
					endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					if (endA<substring.startA) break;
					if (endA<=substring.endA) ReadA.sortedAlignments[i].endInDenseSubstring=true;
					startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					if (startA>=substring.startA && startA<=substring.endA) ReadA.sortedAlignments[i].startInDenseSubstring=true;
					if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null || 
						 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null || 
						 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].inPeriodicSubstring
					   ) continue;
					if ( ( Intervals.isApproximatelyContained(startA,endA,substring.startA,substring.endA) || 
					       Intervals.areApproximatelyIdentical(startA,endA,substring.startA,substring.endA) 
						 ) &&
						 ( ReadA.sortedAlignments[i].inDenseSubstring==null ||
						   substring.endA-substring.startA<ReadA.sortedAlignments[i].inDenseSubstring.endA-ReadA.sortedAlignments[i].inDenseSubstring.startA
						 )
					   ) {
						 ReadA.sortedAlignments[i].inDenseSubstring=substring;
						 if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && ReadA.sortedAlignments[i].isRightMaximalB==1) {
							 length=ReadA.sortedAlignments[i].getALength();
							 if (length<substring.minPrefixLength) substring.minPrefixLength=length;
						 }
						 if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityStart && startA>=substring.startA && startA<substring.maximalStartA) substring.maximalStartA=startA;
						 if (ReadA.sortedAlignments[i].isRightMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityEnd && endA<=substring.endA && endA>substring.maximalEndA) substring.maximalEndA=endA;
						 out++;
					}
				}
				for (i=from+1; i<=lastAlignment; i++) {
					startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					if (startA>=substring.startA && startA<=substring.endA) ReadA.sortedAlignments[i].startInDenseSubstring=true;
					if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null || 
						 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null || 
						 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].inPeriodicSubstring
					   ) continue;
					endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					if ( ( Intervals.isApproximatelyContained(startA,endA,substring.startA,substring.endA) || 
					       Intervals.areApproximatelyIdentical(startA,endA,substring.startA,substring.endA)
  						 ) &&
  						 ( ReadA.sortedAlignments[i].inDenseSubstring==null ||
  						   substring.endA-substring.startA<ReadA.sortedAlignments[i].inDenseSubstring.endA-ReadA.sortedAlignments[i].inDenseSubstring.startA
  						 )
					   ) {
						 ReadA.sortedAlignments[i].inDenseSubstring=substring;
						 if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && ReadA.sortedAlignments[i].isRightMaximalB==1) {
							 length=ReadA.sortedAlignments[i].getALength();
							 if (length<substring.minPrefixLength) substring.minPrefixLength=length;
						 }
						 if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityStart && startA>=substring.startA && startA<substring.maximalStartA) substring.maximalStartA=startA;
						 if (ReadA.sortedAlignments[i].isRightMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityEnd && endA<=substring.endA && endA>substring.maximalEndA) substring.maximalEndA=endA;
						 out++;
					}
				}
			}
			else {				
				from=substring.getFromAlignment(identityThreshold);
				for (i=from; i<=lastAlignment; i++) {
					startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					if (startA>substring.endA) break;
					if (startA>=substring.startA) ReadA.sortedAlignments[i].startInDenseSubstring=true;
					endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					if (endA>=substring.startA && endA<=substring.endA) ReadA.sortedAlignments[i].endInDenseSubstring=true;
					if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null || 
					     ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null || 
						 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].inPeriodicSubstring
					   ) continue;
					if ( ( Intervals.isApproximatelyContained(startA,endA,substring.startA,substring.endA) || 
					       Intervals.areApproximatelyIdentical(startA,endA,substring.startA,substring.endA)
  						 ) &&
  						 ( ReadA.sortedAlignments[i].inDenseSubstring==null ||
  						   substring.endA-substring.startA<ReadA.sortedAlignments[i].inDenseSubstring.endA-ReadA.sortedAlignments[i].inDenseSubstring.startA
  						 )
					   ) {
						 ReadA.sortedAlignments[i].inDenseSubstring=substring;
						 if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && ReadA.sortedAlignments[i].isRightMaximalB==1) {
							 length=ReadA.sortedAlignments[i].getALength();
							 if (length<substring.minPrefixLength) substring.minPrefixLength=length;
						 }
						 if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityStart && startA>=substring.startA && startA<substring.maximalStartA) substring.maximalStartA=startA;
						 if (ReadA.sortedAlignments[i].isRightMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityEnd && endA<=substring.endA && endA>substring.maximalEndA) substring.maximalEndA=endA;
						 out++;
					}
				}
				for (i=from-1; i>=0; i--) {
					endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					if (endA>=substring.startA && endA<=substring.endA) ReadA.sortedAlignments[i].endInDenseSubstring=true;
					if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null || 
						 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null || 
						 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
						 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
						 ReadA.sortedAlignments[i].inPeriodicSubstring
					   ) continue;
					startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					if ( ( Intervals.isApproximatelyContained(startA,endA,substring.startA,substring.endA) || 
					       Intervals.areApproximatelyIdentical(startA,endA,substring.startA,substring.endA)
  						 ) &&
  						 ( ReadA.sortedAlignments[i].inDenseSubstring==null ||
  						   substring.endA-substring.startA<ReadA.sortedAlignments[i].inDenseSubstring.endA-ReadA.sortedAlignments[i].inDenseSubstring.startA
  						 )
				       ) {
						 ReadA.sortedAlignments[i].inDenseSubstring=substring;
						 if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && ReadA.sortedAlignments[i].isRightMaximalB==1) {
							 length=ReadA.sortedAlignments[i].getALength();
							 if (length<substring.minPrefixLength) substring.minPrefixLength=length;
						 }
						 if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityStart && startA>=substring.startA && startA<substring.maximalStartA) substring.maximalStartA=startA;
						 if (ReadA.sortedAlignments[i].isRightMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityEnd && endA<=substring.endA && endA>substring.maximalEndA) substring.maximalEndA=endA;
						 out++;
					}
				}
			}
			if (substring.maximalStartA==Math.POSITIVE_INFINITY) substring.maximalStartA=-1;
			
			// $period$ will be set later, for all dense substrings of substring type in
			// a single step: see procedure $detect$.
		}
		return out;
	}
	
	
	/**
	 * Conceptually identical to $markImpliedAlignmentsAggressive()$.
	 *
	 * @return the first alignment that intersects $substring$, in the prefix order of 
	 * $ReadA.sortedAlignments$. The number of alignments marked by the procedure is 
	 * stored instead in $weakTmp[2]$.
	 */
	private static final int markImpliedAlignmentsAggressive_weak_prefix(DenseSubstring substring, int lastAlignment) {
		final int LENGTH_THRESHOLD = identityThreshold<<1;
		final int MIN_BANDWIDTH = (identityThreshold<<1)/3;
		final int MAX_BANDWIDTH = identityThreshold;
		final int MIN_COUNT = 3;  // Min. n. of events to be used in a $\delta$ histogram
		final int last = substring.getFromAlignmentPrime(identityThreshold,lastAlignment);
		final int first, prefixOrigin;
		int i, j;
		int pos, id, startA, endA, length, from, to;
		int newMinPrefixLength, maxPrefixLength, substringLength, bandwidth, nMarkedAlignments;
		int newStartA, newSumStartA, newNStartA, minPrefixLength_boundary, lastDeltaPointPrefix;
		
		
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markImpliedAlignmentsAggressive_weak_prefix> considering dense substring id="+substring.id+"=["+substring.startA+".."+substring.endA+"] lastAlignment="+lastAlignment);		
		
		// Updating $minPrefixLength$
		pos=substring.sumStartA/substring.nStartA;
		newMinPrefixLength=substring.minPrefixLength;
		maxPrefixLength=-1;
		for (i=last; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			endA=Alignments.alignments[id][4];
			if (endA<=pos) break;
			if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
				 (ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null && ReadA.sortedAlignments[i].impliedByPrefixSubstring!=substring) ||
				 (ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null && ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=substring) ||
				 ReadA.sortedAlignments[i].inPeriodicSubstring
			   ) continue;
			if (isAlignmentImpliedByPrefix(i,substring,false)) {
				length=endA-pos;
				if (length<newMinPrefixLength) newMinPrefixLength=length;
				if (length>maxPrefixLength) maxPrefixLength=length;
			}
		}
		first=i+1;
		if (maxPrefixLength==-1) return first;
		substring.minPrefixLength=newMinPrefixLength;
		maxPrefixLength+=IO.quantum;

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markImpliedAlignmentsAggressive_weak_prefix> 1");
		
		// Filling $deltaPointsPrefix$
		substringLength=substring.endA-substring.startA+1;
		prefixOrigin=substring.minPrefixLength-LENGTH_THRESHOLD; 
		Math.set(countsPrefix,maxPrefixLength-prefixOrigin,0);
		substring.hasLengthBoundariesPrefix=true;
		lastDeltaPointPrefix=-1;
		for (i=last; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			endA=Alignments.alignments[id][4];
			if (endA<=substring.startA) break;
			if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
				 (ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null && ReadA.sortedAlignments[i].impliedByPrefixSubstring!=substring) ||
				 (ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null && ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=substring) ||
				 ReadA.sortedAlignments[i].inPeriodicSubstring
			   ) continue;
			if (isAlignmentImpliedByPrefix(i,substring,false)) {
				length=endA-pos;
				if (length<substringLength-LENGTH_THRESHOLD) {
					lastDeltaPointPrefix++;
					deltaPointsPrefix[lastDeltaPointPrefix].position=length-prefixOrigin;
					deltaPointsPrefix[lastDeltaPointPrefix].mass=1;
					for (j=0; j<=length-prefixOrigin; j++) countsPrefix[j]++;
				}
			}
		}
	
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("markImpliedAlignmentsAggressive_weak_prefix> 2 deltaPointsPrefix:");
	for (int x=0; x<=lastDeltaPointPrefix; x++) IO.printErr(deltaPointsPrefix[x]);
}
		
		// Computing $deltaPrefix$ and finding peaks
		lastLengthBoundaryPrefix=-1;
		if (lastDeltaPointPrefix+1>=MIN_ALIGNMENTS_FOR_LENGTH_BIAS) {
			lastDeltaPointPrefix=Points.sortAndCompact(deltaPointsPrefix,lastDeltaPointPrefix);	
			bandwidth=Points.estimateDensity(deltaPointsPrefix,lastDeltaPointPrefix,deltaPrefix,MIN_BANDWIDTH,MAX_BANDWIDTH);
			for (i=0; i<=lastDeltaPointPrefix; i++) {
				from=(int)deltaPointsPrefix[i].position+1;
				to=Math.min((int)deltaPointsPrefix[i].position+bandwidth,countsPrefix.length-1);
				for (j=from; j<=to; j++) countsPrefix[j]+=deltaPointsPrefix[i].getMass();
			}
			for (i=0; i<=maxPrefixLength-prefixOrigin; i++) {
				if (countsPrefix[i]>=MIN_COUNT) deltaPrefix[i]/=countsPrefix[i];
				else deltaPrefix[i]=0;
			}
			estimateFromDelta(deltaPrefix,maxPrefixLength,substring,true,first,last,prefixOrigin,!inPeriodicInterval(substring,tmpPInterval),LENGTH_THRESHOLD);
		}
		substring.hasLengthBoundariesPrefix=lastLengthBoundaryPrefix>=0;
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("markImpliedAlignmentsAggressive_weak_prefix> 3 lengthBoundariesPrefix:");
	for (int x=0; x<=lastLengthBoundaryPrefix; x++) IO.printErr(lengthBoundariesPrefix[x]);
}
		
		// Marking
		nMarkedAlignments=0;
		newStartA=substring.startA;
		newSumStartA=substring.sumStartA;
		newNStartA=substring.nStartA;
		for (i=last; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			endA=Alignments.alignments[id][4];
			if (endA<pos) break;
			startA=Alignments.alignments[id][3];
			length=(endA-pos)-prefixOrigin;
			j=lastLengthBoundaryPrefix>=0?Arrays.binarySearch(lengthBoundariesPrefix,0,lastLengthBoundaryPrefix+1,length):-1;
			if ( lastLengthBoundaryPrefix>=0 && 
			     ( j>=0 || 
			       (-j-1)%3!=0 || 
			       ( -j-1<=lastLengthBoundaryPrefix && 
					 ( Math.abs(length,lengthBoundariesPrefix[-j])<=identityThreshold || 
			           (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=identityThreshold)
					 )
				   )
			     )
			   ) {
				if (ReadA.sortedAlignments[i].impliedByPrefixSubstring==substring) {
					if (ReadA.sortedAlignments[i].inDenseSubstring==null) ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=ReadA.sortedAlignments[i].impliedByPrefixSubstring;
		   			ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
		   		}
		   		else if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring==null && 
				          ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime==null && 
				          ReadA.sortedAlignments[i].inDenseSubstring==null && 
						  !ReadA.sortedAlignments[i].inPeriodicSubstring &&
				          isAlignmentImpliedByPrefix(i,substring,false) 
				        ) ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=substring;
			}
			else if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring==null && 
			          ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime==null && 
					  ReadA.sortedAlignments[i].inDenseSubstring==null && 
					  !ReadA.sortedAlignments[i].inPeriodicSubstring &&
			          isAlignmentImpliedByPrefix(i,substring,false)
			        ) {
				ReadA.sortedAlignments[i].impliedByPrefixSubstring=substring;
				if (Math.abs(startA,pos)<=identityThreshold && ReadA.sortedAlignments[i].isLeftMaximalB==1) {
					substring.nImpliedAlignmentsPrefix++;
					newSumStartA+=startA;
					newNStartA++;
				}
				if (startA<newStartA) newStartA=startA;
				if (ReadA.sortedAlignments[i].isRightMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityEnd && endA<=substring.endA && endA>substring.maximalEndA) substring.maximalEndA=endA;
				nMarkedAlignments++;
			}
		}
		substring.sumStartA=newSumStartA;
		substring.nStartA=newNStartA;
		substring.startA=newStartA;
		
		// Updating $minPrefixLength$
		pos=substring.sumStartA/substring.nStartA;
		minPrefixLength_boundary=-1;
		length=substring.minPrefixLength-prefixOrigin;
		j=lastLengthBoundaryPrefix>=0?Arrays.binarySearch(lengthBoundariesPrefix,0,lastLengthBoundaryPrefix+1,length):-1;
		if ( lastLengthBoundaryPrefix>=0 && 
		     ( j>=0 || 
		       (-j-1)%3!=0 || 
		       ( -j-1<=lastLengthBoundaryPrefix &&
				 ( Math.abs(length,lengthBoundariesPrefix[-j])<=identityThreshold || 
		           (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=identityThreshold)
				 )
			   )
		     )
		   ) {
			minPrefixLength_boundary=j;
			substring.minPrefixLength=Math.POSITIVE_INFINITY;
			for (i=last; i>=0; i--) {
				endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
				if (endA<=pos) break;
				startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=substring || 
					 ReadA.sortedAlignments[i].isLeftMaximalB!=1 || 
					 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
					 ReadA.sortedAlignments[i].inPeriodicSubstring
				   ) continue;
				if (Math.abs(startA,substring.startA)>identityThreshold) continue;
				length=ReadA.sortedAlignments[i].getALength();
			 	if (length<substring.minPrefixLength) substring.minPrefixLength=length;
			}
			if (substring.minPrefixLength==Math.POSITIVE_INFINITY) substring.minPrefixLength=-1;
		}
		
		// Updating $lengthBiasIntervals$
		if (lastLengthBoundaryPrefix>=0) {
			for (i=minPrefixLength_boundary==-1?0:(minPrefixLength_boundary/3+1)*3; i<lastLengthBoundaryPrefix; i+=3) {
				lastLengthBiasInterval++;
				lengthBiasIntervals[lastLengthBiasInterval].firstPosition=pos+prefixOrigin+lengthBoundariesPrefix[i];
				lengthBiasIntervals[lastLengthBiasInterval].center=pos+prefixOrigin+lengthBoundariesPrefix[i+1];
				lengthBiasIntervals[lastLengthBiasInterval].lastPosition=pos+prefixOrigin+lengthBoundariesPrefix[i+2];
				lengthBiasIntervals[lastLengthBiasInterval].readB=ReadA.id;
				lengthBiasIntervals[lastLengthBiasInterval].orientation=true;
			}
		}
		
		weakTmp[2]=nMarkedAlignments;
		return first;
	}
	
	
	/**
	 * Conceptually identical to $markImpliedAlignmentsAggressive()$.
	 *
	 * @return the last alignment that intersects $substring$, in the suffix order of 
	 * $ReadA.sortedAlignments$. The number of alignments marked by the procedure is 
	 * stored instead in $weakTmp[2]$.
	 */
	private static final int markImpliedAlignmentsAggressive_weak_suffix(DenseSubstring substring, int lastAlignment) {
		final int LENGTH_THRESHOLD = identityThreshold<<1;
		final int MIN_BANDWIDTH = (identityThreshold<<1)/3;
		final int MAX_BANDWIDTH = identityThreshold;
		final int MIN_COUNT = 3;  // Min. n. of events to be used in a $\delta$ histogram
		final int first = substring.getFromAlignment(identityThreshold);
		final int last, suffixOrigin;
		int i, j;
		int pos, id, startA, endA, length, from, to, nMarkedAlignments;
		int newMinSuffixLength, maxSuffixLength, substringLength, bandwidth;
		int newEndA, newSumEndA, newNEndA, minSuffixLength_boundary, lastDeltaPointSuffix;
		
		
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markImpliedAlignmentsAggressive_weak_suffix> considering dense substring id="+substring.id+"=["+substring.startA+".."+substring.endA+"] first="+first+" lastAlignment="+lastAlignment);
		
		
		// Updating $minSuffixLength$
		pos=substring.sumEndA/substring.nEndA;
		newMinSuffixLength=substring.minSuffixLength;
		maxSuffixLength=-1;
		for (i=first; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA>=pos) break;
			if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
				 (ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null && ReadA.sortedAlignments[i].impliedBySuffixSubstring!=substring) ||
				 (ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null && ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=substring) ||
				 ReadA.sortedAlignments[i].inPeriodicSubstring
			   ) continue;
			if (isAlignmentImpliedBySuffix(i,substring,false)) {
				length=pos-startA;			
				if (length<newMinSuffixLength) newMinSuffixLength=length;
				if (length>maxSuffixLength) maxSuffixLength=length;
			}
		}
		last=i-1;
		if (maxSuffixLength==-1) return last;
		substring.minSuffixLength=newMinSuffixLength;
		maxSuffixLength+=IO.quantum;		
		
		// Filling $deltaPointsSuffix$
		substringLength=substring.endA-substring.startA+1;		
		suffixOrigin=substring.minSuffixLength-LENGTH_THRESHOLD; 
		Math.set(countsSuffix,maxSuffixLength-suffixOrigin,0);
		substring.hasLengthBoundariesSuffix=true;
		lastDeltaPointSuffix=-1;
		for (i=first; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA>=substring.endA) break;
			if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
				 (ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null && ReadA.sortedAlignments[i].impliedBySuffixSubstring!=substring) ||
				 (ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null && ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=substring) ||
				 ReadA.sortedAlignments[i].inPeriodicSubstring
			   ) continue;
			if (isAlignmentImpliedBySuffix(i,substring,false)) {
				length=pos-startA;
				if (length<substringLength-LENGTH_THRESHOLD) {
					lastDeltaPointSuffix++;
					deltaPointsSuffix[lastDeltaPointSuffix].position=length-suffixOrigin;
					deltaPointsSuffix[lastDeltaPointSuffix].mass=1;
					for (j=0; j<=length-suffixOrigin; j++) countsSuffix[j]++;
				}
			}
		}
		
		// Computing $deltaSuffix$ and finding peaks
		lastLengthBoundarySuffix=-1;
		if (lastDeltaPointSuffix+1>=MIN_ALIGNMENTS_FOR_LENGTH_BIAS) {
			lastDeltaPointSuffix=Points.sortAndCompact(deltaPointsSuffix,lastDeltaPointSuffix);			
			bandwidth=Points.estimateDensity(deltaPointsSuffix,lastDeltaPointSuffix,deltaSuffix,MIN_BANDWIDTH,MAX_BANDWIDTH);
			for (i=0; i<=lastDeltaPointSuffix; i++) {
				from=(int)deltaPointsSuffix[i].position+1;
				to=Math.min((int)deltaPointsSuffix[i].position+bandwidth,countsSuffix.length-1);
				for (j=from; j<=to; j++) countsSuffix[j]+=deltaPointsSuffix[i].getMass();
			}
			for (i=0; i<=maxSuffixLength-suffixOrigin; i++) {
				if (countsSuffix[i]>=MIN_COUNT) deltaSuffix[i]/=countsSuffix[i];
				else deltaSuffix[i]=0;
			}
			estimateFromDelta(deltaSuffix,maxSuffixLength,substring,false,first,last,suffixOrigin,!inPeriodicInterval(substring,tmpPInterval),LENGTH_THRESHOLD);
		}
		substring.hasLengthBoundariesSuffix=lastLengthBoundarySuffix>=0;
		
		// Marking
		nMarkedAlignments=0;
		newEndA=substring.endA;
		newSumEndA=substring.sumEndA;
		newNEndA=substring.nEndA;
		for (i=first; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA>=pos) break;
			endA=Alignments.alignments[id][4];
			length=(pos-startA)-suffixOrigin;
			j=lastLengthBoundarySuffix>=0?Arrays.binarySearch(lengthBoundariesSuffix,0,lastLengthBoundarySuffix+1,length):-1;
			if ( lastLengthBoundarySuffix>=0 && 
			     ( j>=0 || 
			       (-j-1)%3!=0 || 
				   ( -j-1<=lastLengthBoundarySuffix &&
			         ( Math.abs(length,lengthBoundariesSuffix[-j])<=identityThreshold || 
			           (-j-1==0?false:Math.abs(length,lengthBoundariesSuffix[-j-3])<=identityThreshold)
					 )
				   )
			     )
			   ) {
				if (ReadA.sortedAlignments[i].impliedBySuffixSubstring==substring) {
					if (ReadA.sortedAlignments[i].inDenseSubstring==null) ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=ReadA.sortedAlignments[i].impliedBySuffixSubstring;
		   			ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
		   		}
		   		else if ( ReadA.sortedAlignments[i].impliedBySuffixSubstring==null && 
				          ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime==null && 
						  ReadA.sortedAlignments[i].inDenseSubstring==null &&
						  !ReadA.sortedAlignments[i].inPeriodicSubstring &&
				          isAlignmentImpliedBySuffix(i,substring,false)
			            ) {
					ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=substring;
				}
			}
			else if ( ReadA.sortedAlignments[i].impliedBySuffixSubstring==null && 
			          ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime==null && 
					  ReadA.sortedAlignments[i].inDenseSubstring==null &&
					  !ReadA.sortedAlignments[i].inPeriodicSubstring &&
			          isAlignmentImpliedBySuffix(i,substring,false)
			        ) {
				ReadA.sortedAlignments[i].impliedBySuffixSubstring=substring;
				ReadA.sortedAlignments[i].inDenseSubstring=null;
				if (Math.abs(endA,pos)<=identityThreshold && ReadA.sortedAlignments[i].isRightMaximalB==1) {
					substring.nImpliedAlignmentsSuffix++;
					newSumEndA+=endA;
					newNEndA++;
				}
				if (endA>newEndA) newEndA=endA;
				if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityStart && startA>=substring.startA && startA<substring.maximalStartA) substring.maximalStartA=startA;
				nMarkedAlignments++;
			}
		}
		substring.sumEndA=newSumEndA;
		substring.nEndA=newNEndA;
		substring.endA=newEndA;
		
		// Updating $minSuffixLength$
		pos=substring.sumEndA/substring.nEndA;
		minSuffixLength_boundary=-1;
		length=substring.minSuffixLength-suffixOrigin;
		j=lastLengthBoundarySuffix>=0?Arrays.binarySearch(lengthBoundariesSuffix,0,lastLengthBoundarySuffix+1,length):-1;
		if ( lastLengthBoundarySuffix>=0 && 
		     ( j>=0 || 
		       (-j-1)%3!=0 || 
		       ( -j-1<=lastLengthBoundarySuffix && 
				 ( Math.abs(length,lengthBoundariesSuffix[-j])<=identityThreshold || 
		           (-j-1==0?false:Math.abs(length,lengthBoundariesSuffix[-j-3])<=identityThreshold)
				 )
			   )
		     )
		   ) {
			minSuffixLength_boundary=j;
			substring.minSuffixLength=Math.POSITIVE_INFINITY;
			for (i=first; i<=lastAlignment; i++) {
				startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				if (startA>=pos) break;
				endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
				if ( ReadA.sortedAlignments[i].impliedBySuffixSubstring!=substring || 
					 ReadA.sortedAlignments[i].isLeftMaximalB!=1 || 
					 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
				     ReadA.sortedAlignments[i].inPeriodicSubstring
				   ) continue;
				if (Math.abs(endA,substring.endA)>identityThreshold) continue;
				length=ReadA.sortedAlignments[i].getALength();
			 	if (length<substring.minSuffixLength) substring.minSuffixLength=length;
			}
			if (substring.minSuffixLength==Math.POSITIVE_INFINITY) substring.minSuffixLength=-1;
		}
		
		// Updating $lengthBiasIntervals$
		if (lastLengthBoundarySuffix>=0) {
			for (i=minSuffixLength_boundary==-1?0:(minSuffixLength_boundary/3+1)*3; i<lastLengthBoundarySuffix; i+=3) {
				lastLengthBiasInterval++;
				lengthBiasIntervals[lastLengthBiasInterval].firstPosition=pos-suffixOrigin-lengthBoundariesSuffix[i+2];
				lengthBiasIntervals[lastLengthBiasInterval].center=pos-suffixOrigin-lengthBoundariesSuffix[i+1];
				lengthBiasIntervals[lastLengthBiasInterval].lastPosition=pos-suffixOrigin-lengthBoundariesSuffix[i];
				lengthBiasIntervals[lastLengthBiasInterval].readB=ReadA.id;
				lengthBiasIntervals[lastLengthBiasInterval].orientation=true;
			}
		}
		
		weakTmp[2]=nMarkedAlignments;
		return last;
	}
	
	
	/**
	 * Merges overlapping non-weak dense substrings of substring type into a single 
	 * substring, and ensures that no dense substring is contained in a non-weak dense 
	 * substring of substring type.
	 *
	 * Remark: adjacent non-weak dense substrings of substring type are not merged, since 
	 * they could be distinct substring repeats.
	 *
	 * Remark: weak dense substrings of substring type: (1) are not merged with each 
	 * other or with other substrings, unless their intervals are approximately identical;
	 * (2) are merged with a (weak or non-weak) dense substring of substring type iff they
	 * are contained in the latter; (3) they do not delete dense substrings contained in 
	 * them; (4) they do not unmark as implied the alignments of dense substrings 
	 * contained in them.
	 *
	 * Remark: the procedure sorts $substrings$ by $DenseSubstring.PRSR_STARTA$.
	 */
	private static final void cleanDenseSubstringsOfSubstringType(int lastAlignment) {
		boolean atLeastOneConcatenation;
		int i, j, k;
		int found, firstSubstringReplication, firstJForNextI;
		DenseSubstring tmpSubstring;		
		
		// Ensuring the necessary order in $substrings$
		if (DenseSubstring.order!=DenseSubstring.PRSR_STARTA) {
			DenseSubstring.order=DenseSubstring.PRSR_STARTA;
			if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		}
		firstSubstringReplication=lastSubstring+1;
		for (i=lastSubstring; i>=0; i--) {
			if (!substrings[i].substringReplication) {
				firstSubstringReplication=i+1;
				break;
			}
		}
		if (i==-1) firstSubstringReplication=0;
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].discarded=false;
			substrings[i].representative=substrings[i];
		}
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("cleanDenseSubstringsOfSubstringType> dense substrings before merging:");
	for (int x=0; x<=lastSubstring; x++) IO.printErr(substrings[x]);
}		
		
		// Merging non-weak dense substrings of substring type with large intersection
		atLeastOneConcatenation=false;
		for (i=firstSubstringReplication; i<lastSubstring; i++) {
			if (substrings[i].discarded || substrings[i].isWeak) continue;
			j=getConnectedComponent_substringType(i);
			if (j>0) atLeastOneConcatenation=true;
		}
		if (atLeastOneConcatenation) updateConcatenatedSumN(firstSubstringReplication,lastSubstring);

		// Discarding prefix/suffix substrings contained in a non-weak dense substring of
		// substring type
		i=firstSubstringReplication;
		j=0; firstJForNextI=-1;
		while (i<=lastSubstring) {
			if (substrings[i].isWeak || j==firstSubstringReplication || substrings[j].startA>=substrings[i].endA) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (i<lastSubstring && firstJForNextI==-1 && substrings[j].startA>=substrings[i+1].startA) firstJForNextI=j;
			if (substrings[j].discarded) {
				j++;
				continue;
			}
			if ( Intervals.isApproximatelyContained_lowQuality(substrings[j].startA,substrings[j].endA,substrings[i].startA,substrings[i].endA,ReadA.id) ||
				 Intervals.areApproximatelyIdentical_lowQuality(substrings[j].startA,substrings[j].endA,substrings[i].startA,substrings[i].endA,ReadA.id)
			   ) {
				substrings[j].discarded=true;
				substrings[j].representative=substrings[i];
			}
			j++;
		}

		// Discarding dense substrings of substring type contained in a non-weak dense
		// substring of substring type.
		for (i=firstSubstringReplication; i<=lastSubstring; i++) {
			if (substrings[i].discarded) continue;
			found=-1;
			for (j=i-1; j>=firstSubstringReplication; j--) {
				if (substrings[j].isWeak && !substrings[i].isWeak) continue;
				if ( Intervals.isApproximatelyContained_lowQuality(substrings[i].startA,substrings[i].endA,substrings[j].startA,substrings[j].endA,ReadA.id) ||
				     ( Intervals.areApproximatelyIdentical_lowQuality(substrings[i].startA,substrings[i].endA,substrings[j].startA,substrings[j].endA,ReadA.id) &&
					   substrings[i].length()<substrings[j].length()
					 )
				   ) {
					found=j;
					break;
				}
			}
			if (found!=-1) {
				substrings[i].discarded=true;
				substrings[i].representative=substrings[found];
				continue;
			}
			for (j=i+1; j<=lastSubstring; j++) {
				if (substrings[j].startA>=substrings[i].endA) break;
				if (substrings[j].isWeak && !substrings[i].isWeak) continue;
				if ( Intervals.isApproximatelyContained_lowQuality(substrings[i].startA,substrings[i].endA,substrings[j].startA,substrings[j].endA,ReadA.id) ||
					 ( Intervals.areApproximatelyIdentical_lowQuality(substrings[i].startA,substrings[i].endA,substrings[j].startA,substrings[j].endA,ReadA.id) &&
					   substrings[i].length()<substrings[j].length()
					 )
				   ) {
					found=j;
					break;
				}
			}
			if (found!=-1) {
				substrings[i].discarded=true;
				substrings[i].representative=substrings[found];
			}
		}
		
		// Removing discarded substrings and assigning them the closure of their
		// $representative$ pointer. The $representative$ pointer of all non-discarded
		// substrings is equal to the substring itself.
		k=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].discarded) {
				substrings[i].representative=substrings[i].closure();
				if (substrings[i].substringReplication && substrings[i].minPrefixLength<substrings[i].representative.minPrefixLength) substrings[i].representative.minPrefixLength=substrings[i].minPrefixLength;
				continue;
			}
			k++;
			if (k!=i) {
				tmpSubstring=substrings[k];
				substrings[k]=substrings[i];
				substrings[i]=tmpSubstring;
			}
		}
		lastSubstring=k;
		
		// Updating pointers from alignments to implied substrings
		for (i=0; i<=lastAlignment; i++) {
			if (ReadA.sortedAlignments[i].inDenseSubstring!=null && ReadA.sortedAlignments[i].inDenseSubstring.discarded) {
				ReadA.sortedAlignments[i].inDenseSubstring=ReadA.sortedAlignments[i].inDenseSubstring.representative;
				continue;
			}
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null && ReadA.sortedAlignments[i].impliedByDenseSubstring.discarded) {
				tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring.representative;
				if (tmpSubstring.substringReplication) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
					ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring;
				}
				else ReadA.sortedAlignments[i].impliedByDenseSubstring=tmpSubstring;
			}
			if (ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null && ReadA.sortedAlignments[i].impliedByPrefixSubstring.discarded) {
				tmpSubstring=ReadA.sortedAlignments[i].impliedByPrefixSubstring.representative;
				if (tmpSubstring.substringReplication) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
					ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring;
				}
				else ReadA.sortedAlignments[i].impliedByPrefixSubstring=tmpSubstring;
			}
			if (ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null && ReadA.sortedAlignments[i].impliedBySuffixSubstring.discarded) {
				tmpSubstring=ReadA.sortedAlignments[i].impliedBySuffixSubstring.representative;
				if (tmpSubstring.substringReplication) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
					ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring;
				}
				else ReadA.sortedAlignments[i].impliedBySuffixSubstring=tmpSubstring;
			}
			if (ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null && ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime.discarded) {
				tmpSubstring=ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime.representative;
				if (tmpSubstring.substringReplication) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
					ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring;
				}
				else ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=tmpSubstring;
			}
			if (ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null && ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime.discarded) {
				tmpSubstring=ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime.representative;
				if (tmpSubstring.substringReplication) {
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;
					ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;
					ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
					ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring;
				}
				else ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=tmpSubstring;
			}
		}
	}
	
	
	/**
	 * Merges non-weak dense substrings with large intersection. Straddling substrings 
	 * could originate from the merge of prefix/suffix substrings done by 
	 * $getConnectedComponent_concatenation$, but not from the initial detection of dense 
	 * substrings of substring type, which explicitly forbids straddling and containment.
	 *
	 * Remark: the procedure assumes the $representative$ field of every substring to be
	 * initialized to the substring itself.
	 *
	 * Remark: there is no need to handle periods at this stage, since periods have yet to
	 * be assigned to dense substrings of substring type.
	 *
	 * @return the number of substrings in the connected component of $fromSubstring$,
	 * excluding $fromSubstring$.
	 */
	private static final int getConnectedComponent_substringType(int fromSubstring) {
		final int INTERSECTION_THRESHOLD = Alignments.minAlignmentLength;
		int i, j;
		int top, count, newStartA, newEndA;

		newStartA=substrings[fromSubstring].startA;
		newEndA=substrings[fromSubstring].endA;
		count=0; top=0; stack[0]=fromSubstring;
		while (top>=0) {
			i=stack[top--];
			for (j=fromSubstring+1; j<=lastSubstring; j++) {
				if (j==i || substrings[j].discarded || !substrings[j].substringReplication || substrings[j].isWeak) continue;
				if (Intervals.intersectionLength(substrings[i].startA,substrings[i].endA,substrings[j].startA,substrings[j].endA)>=INTERSECTION_THRESHOLD) { 
					substrings[j].discarded=true;
					substrings[j].representative=substrings[fromSubstring];
					if (substrings[j].startA<newStartA) newStartA=substrings[j].startA;
					if (substrings[j].endA>newEndA) newEndA=substrings[j].endA;
					stack[++top]=j;
					count++;
				}
			}
		}
		substrings[fromSubstring].startA=newStartA;
		substrings[fromSubstring].endA=newEndA;
		return count;
	}
	
	
	/**
	 * Detects whether a dense substring is adjacent to a high-quality substring to
	 * the left and to the right of its readA.
	 */
	private static final void setMaximality() {
		for (int i=0; i<=lastSubstring; i++) {
			substrings[i].isLeftMaximal=Reads.isLeftMaximal(substrings[i].startA,ReadA.id,true);
			substrings[i].isRightMaximal=Reads.isRightMaximal(substrings[i].endA,ReadA.id,true);
		}
	}
	
	
	/**
	 * Estimates the period of all dense substrings of substring type, including weak
	 * substrings.
	 *
	 * Remark: the period of other dense substrings is estimated by 
	 * $removeBoundariesAtPeriods()$.
	 */
	private static final void estimateSubstringPeriods(int firstAlignment, int lastAlignment) {
		int i;
		int lastShift, firstAlignmentForNextSubstring;
		int firstAlignmentForNextSubstringPrime;

		firstAlignmentForNextSubstring=firstAlignment;
		firstAlignmentForNextSubstringPrime=firstAlignment;
		for (i=0; i<=lastSubstring; i++) {
			if (!substrings[i].substringReplication) continue;			
			lastShift=getShiftsSubstring(substrings[i].startA,substrings[i].endA,firstAlignmentForNextSubstring,lastAlignment,i<lastSubstring?substrings[i+1].startA:Math.POSITIVE_INFINITY);
			firstAlignmentForNextSubstring=shiftTmp[1];
			if (lastShift>=0) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("ESTIMATING PERIOD FOR DENSE SUBSTRING OF SUBSTRING TYPE ["+substrings[i].startA+".."+substrings[i].endA+"], read "+ReadA.id+":");
				substrings[i].period=PeriodicSubstrings.estimatePeriod(lengthPoints,lastShift,PeriodicSubstrings.MIN_PERIOD,PeriodicSubstrings.minIntervalLength,PeriodicSubstrings.minLocalMaxDistance);
				substrings[i].strongPeriodSignal=substrings[i].period>0&&PeriodicSubstrings.periodTmp[4]>1;
			}
		}
	}

	
	/**
	 * Remark: the procedure assumes $PeriodicSubstrings.substrings$ to be sorted by 
	 * $minStartA$.
	 *
	 * @param tmpPeriodic temporary space;
	 * @return TRUE iff $[start..end]$ is approximately identical to, or contained in, a 
	 * short-period periodic substring.
	 */
	private static final boolean inPeriodicSubstring(int start, int end, PeriodicSubstring tmpPeriodic) {
		int i, j;
		PeriodicSubstring periodicSubstring;
		
		tmpPeriodic.minStartA=start;
		i=Arrays.binarySearch(PeriodicSubstrings.substrings,0,PeriodicSubstrings.lastSubstring+1,tmpPeriodic);
		if (i<0) {
			i=-i-1;
			i=Math.min(i,PeriodicSubstrings.lastSubstring);
		}
		for (j=i; j>=0; j--) {
			periodicSubstring=PeriodicSubstrings.substrings[j];
			if ( Intervals.areApproximatelyIdentical(start,end,periodicSubstring.minStartA,periodicSubstring.maxEndA) ||
				 Intervals.isApproximatelyContained(start,end,periodicSubstring.minStartA,periodicSubstring.maxEndA)
			   ) return true;
		}
		for (j=i+1; j<=PeriodicSubstrings.lastSubstring; j++) {
			periodicSubstring=PeriodicSubstrings.substrings[j];
			if (periodicSubstring.minStartA>=end) break;
			if ( Intervals.areApproximatelyIdentical(start,end,periodicSubstring.minStartA,periodicSubstring.maxEndA) ||
				 Intervals.isApproximatelyContained(start,end,periodicSubstring.minStartA,periodicSubstring.maxEndA)
			   ) return true;
		}		
		return false;
	}
	
	
	/**
	 * Remark: the procedure assumes $PeriodicSubstrings.substrings$ to be sorted by 
	 * $minStartA$.
	 *
	 * @param tmpPeriodic temporary space;
	 * @return TRUE iff $[start..end]$ is approximately identical to, or contains, a 
	 * periodic substring.
	 */
	private static final boolean containsPeriodicSubstring(int start, int end, PeriodicSubstring tmpPeriodic) {
		int i, j;
		PeriodicSubstring periodicSubstring;
		
		tmpPeriodic.minStartA=start;
		i=Arrays.binarySearch(PeriodicSubstrings.substrings,0,PeriodicSubstrings.lastSubstring+1,tmpPeriodic);
		if (i<0) {
			i=-i-1;
			i=Math.min(i,PeriodicSubstrings.lastSubstring);
		}
		for (j=i; j>=0; j--) {
			periodicSubstring=PeriodicSubstrings.substrings[j];
			if ( Intervals.areApproximatelyIdentical(start,end,periodicSubstring.minStartA,periodicSubstring.maxEndA) ||
				 Intervals.isApproximatelyContained(periodicSubstring.minStartA,periodicSubstring.maxEndA,start,end)
			   ) return true;
		}
		for (j=i+1; j<=PeriodicSubstrings.lastSubstring; j++) {
			periodicSubstring=PeriodicSubstrings.substrings[j];
			if (periodicSubstring.minStartA>=end) break;
			if ( Intervals.areApproximatelyIdentical(start,end,periodicSubstring.minStartA,periodicSubstring.maxEndA) ||
				 Intervals.isApproximatelyContained(periodicSubstring.minStartA,periodicSubstring.maxEndA,start,end)
			   ) return true;
		}
		return false;
	}
	
	
	/**
	 * Tries to extend to the left and to the right every dense substring of substring
	 * type, by greedily adding maximal, non-implied alignments, that are at distance 
	 * at most $identityThreshold$ from the current endpoint, and that are compatible with
	 * the dense substring. This is useful, since the substring DAG built by 
	 * $getSubstrings()$ does not connect two intervals if their starting or ending 
	 * positions are at distance at most $identityThreshold$. Extending a dense substring
	 * helps factorization by e.g. putting the first position of the substring farther 
	 * away from a peak of left-maximal events, which could go undetected if too close to 
	 * the first position of the substring.
	 *
	 * The procedure tries also to extend prefix or suffix substrings using alignments
	 * that are marked with $inPeriodicSubstring=true$. This is useful, since a prefix or
	 * suffix substring that overlaps with a periodic substring might have some of its 
	 * alignments marked as $inPeriodicSubstring$: such alignments are not used by 
	 * $getSubstrings()$, so the resulting substring might become erroneously shorter.
	 * See procedure $extendSubstrings_impl()$ for details. Due to the latter procedure, 
	 * $nAlignments$ might increase after this procedure completes.
	 *
	 * Remark: the procedure might increase $nAlignments$, and it might alter the order of 
	 * $ReadA.sortedAlignments[0..nAlignments-1]$, but it restores it at the end. The
	 * procedure alters also the order of $substrings[0..lastSubstring]$, but it does not
	 * restore it at the end.
	 *
	 * Remark: the procedure just sets fields $impliedBy{Prefix,Suffix}Substring$ and
	 * $inDenseSubstring$ of an alignment, but it does not set field 
	 * $impliedByDenseSubstring$.
	 *
	 * @param nAlignments the first $nAlignments$ of $ReadA.sortedAlignments$ are assumed
	 * not to be implied by periodic substrings. $ReadA.sortedAlignments[nAlignments..]$
	 * is assumed to be sorted by startA: this is preserved at the end of the procedure.
	 * @return the new (possibly bigger) value of $nAlignments$.
	 */
	private static final int extendSubstrings(int nAlignments) {
		boolean periodicAlignmentChanged;
		int i, j, k;
		int firstJForNextI, firstAlignment, maxStart, minStart, minEnd, maxEnd;
		int alignmentID, alignmentIDPrime, alignmentStart, alignmentStartPrime, alignmentEnd, alignmentEndPrime;
		int previousAlignmentOrder;
		Alignment tmpAlignment;
		
		// Initializing temporary space
		previousAlignmentOrder=Alignment.order;
		for (i=0; i<=lastSubstring; i++) {
			if (!substrings[i].substringReplication) {
				i++;
				continue;
			}
			substrings[i].leftAlignment=-1;
			substrings[i].rightAlignment=-1;
		}
		periodicAlignmentChanged=false;
		
		// Adding alignments to the left of dense substrings of substring type
		Alignment.order=Alignment.STARTA_ENDA;
		if (nAlignments>1) Arrays.sort(ReadA.sortedAlignments,0,nAlignments);
		DenseSubstring.order=DenseSubstring.STARTA;
		if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		i=0; j=0; firstJForNextI=-1; firstAlignment=-1; maxStart=-1;
		while (i<=lastSubstring) {
			if (!substrings[i].substringReplication) {
				i++;
				continue;
			}
			if (j==nAlignments || Alignments.alignments[ReadA.sortedAlignments[j].id][3]>substrings[i].endA) {
				if (firstAlignment!=-1) {
					minStart=substrings[i].startA;
					for (k=firstAlignment; k>=0; k--) {
						alignmentIDPrime=ReadA.sortedAlignments[k].id;
						alignmentStartPrime=Alignments.alignments[alignmentIDPrime][3];
						alignmentEndPrime=Alignments.alignments[alignmentIDPrime][4];
						if (alignmentStartPrime<minStart-identityThreshold) break;
						if (alignmentEndPrime<=substrings[i].startA+identityThreshold || alignmentEndPrime>substrings[i].endA+identityThreshold) continue;
						if (ReadA.sortedAlignments[k].inDenseSubstring==substrings[i]) {
							minStart=alignmentStartPrime;
							continue;
						}
						if ( ReadA.sortedAlignments[k].inDenseSubstring!=null ||
							 ReadA.sortedAlignments[k].impliedByPrefixSubstring!=null || 
							 ReadA.sortedAlignments[k].impliedBySuffixSubstring!=null || 
							 ReadA.sortedAlignments[k].impliedByPrefixSubstringPrime!=null || 
							 ReadA.sortedAlignments[k].impliedBySuffixSubstringPrime!=null ||
							 ReadA.sortedAlignments[k].lowQualityStart || ReadA.sortedAlignments[k].lowQualityEnd ||
							 (ReadA.sortedAlignments[k].isLeftMaximalB!=1 && ReadA.sortedAlignments[k].isRightMaximalB!=1) ||
							 ReadA.sortedAlignments[k].inPeriodicSubstring
						   ) continue;
						minStart=alignmentStartPrime;
						substrings[i].sumStartA+=alignmentStartPrime; 
						substrings[i].nStartA++;
						ReadA.sortedAlignments[k].inDenseSubstring=substrings[i];
					}
					substrings[i].leftAlignment=minStart;  // Using an unrelated variable as temporary space
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; firstAlignment=-1; maxStart=-1;
				continue;
			}
			alignmentID=ReadA.sortedAlignments[j].id;
			alignmentStart=Alignments.alignments[alignmentID][3];
			alignmentEnd=Alignments.alignments[alignmentID][4];
			if (alignmentEnd<substrings[i].startA) {
				j++;
				continue;
			}
			if (i<lastSubstring && firstJForNextI==-1 && alignmentEnd>=substrings[i+1].startA) firstJForNextI=j;
			if (alignmentStart>=substrings[i].startA-identityThreshold && alignmentStart<=substrings[i].startA && alignmentStart>maxStart) {
				maxStart=alignmentStart;
				firstAlignment=j;
			}
			j++;
		}
		
		// Adding alignments to the right of prefix substrings
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastSubstring) {
			if (!substrings[i].prefixReplication || substrings[i].suffixReplication || substrings[i].singleDeletionReplication || substrings[i].substringReplication) {
				i++;
				firstAlignment=-1;
				continue;
			}
			if (j==nAlignments || ReadA.sortedAlignments[j].startA()>substrings[i].startA+identityThreshold) {
				if (firstAlignment!=-1) {
					maxEnd=extendSubstrings_impl(firstAlignment,j-1,substrings[i],nAlignments);
					if (maxEnd<0) {
						substrings[i].leftAlignment=-1-maxEnd;
						periodicAlignmentChanged=true;
					}
					else substrings[i].leftAlignment=maxEnd;
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; firstAlignment=-1;
				continue;
			}
			if (i<lastSubstring && firstJForNextI==-1 && ReadA.sortedAlignments[j].endA()>=substrings[i+1].startA) firstJForNextI=j;
			if (ReadA.sortedAlignments[j].startA()<substrings[i].startA-identityThreshold) {
				j++;
				continue;
			}
			if (firstAlignment==-1) firstAlignment=j;
			j++;
		}
		
		// Adding alignments to the right of dense substrings of substring type
		Alignment.order=Alignment.ENDA_REVERSE;
		if (nAlignments>1) Arrays.sort(ReadA.sortedAlignments,0,nAlignments);
		DenseSubstring.order=DenseSubstring.ENDA_REVERSE;
		if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		i=0; j=0; firstJForNextI=-1; firstAlignment=-1; minEnd=Math.POSITIVE_INFINITY;
		while (i<=lastSubstring) {
			if (!substrings[i].substringReplication) {
				i++;
				continue;
			}
			if (j==nAlignments || Alignments.alignments[ReadA.sortedAlignments[j].id][4]<substrings[i].startA) {
				if (firstAlignment!=-1) {
					maxEnd=substrings[i].endA;
					for (k=firstAlignment; k>=0; k--) {
						alignmentIDPrime=ReadA.sortedAlignments[k].id;
						alignmentStartPrime=Alignments.alignments[alignmentIDPrime][3];
						alignmentEndPrime=Alignments.alignments[alignmentIDPrime][4];
						if (alignmentEndPrime>maxEnd+identityThreshold) break;
						if (alignmentStartPrime>=substrings[i].endA-identityThreshold || alignmentStartPrime<substrings[i].startA-identityThreshold) continue;
						if (ReadA.sortedAlignments[k].inDenseSubstring==substrings[i]) {
							maxEnd=alignmentEndPrime;
							continue;
						}
						if ( ReadA.sortedAlignments[k].inDenseSubstring!=null ||
							 ReadA.sortedAlignments[k].impliedByPrefixSubstring!=null || 
							 ReadA.sortedAlignments[k].impliedBySuffixSubstring!=null || 
							 ReadA.sortedAlignments[k].impliedByPrefixSubstringPrime!=null || 
							 ReadA.sortedAlignments[k].impliedBySuffixSubstringPrime!=null ||
							 ReadA.sortedAlignments[k].lowQualityStart || ReadA.sortedAlignments[k].lowQualityEnd ||
							 (ReadA.sortedAlignments[k].isRightMaximalB!=1 && ReadA.sortedAlignments[k].isLeftMaximalB!=1) ||
							 ReadA.sortedAlignments[k].inPeriodicSubstring
						   ) continue;
						maxEnd=alignmentEndPrime;
						substrings[i].sumEndA+=alignmentEndPrime;
						substrings[i].nEndA++;
						ReadA.sortedAlignments[k].inDenseSubstring=substrings[i];
					}
					substrings[i].rightAlignment=maxEnd;  // Using as temporary space an unrelated variable
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; firstAlignment=-1; minEnd=Math.POSITIVE_INFINITY;
				continue;
			}
			alignmentID=ReadA.sortedAlignments[j].id;
			alignmentStart=Alignments.alignments[alignmentID][3];
			alignmentEnd=Alignments.alignments[alignmentID][4];
			if (alignmentStart>substrings[i].endA) {
				j++;
				continue;
			}
			if (i<lastSubstring && firstJForNextI==-1 && alignmentStart<=substrings[i+1].endA) firstJForNextI=j;
			if (alignmentEnd>=substrings[i].endA && alignmentEnd<=substrings[i].endA+identityThreshold && alignmentEnd<minEnd) {
				minEnd=alignmentEnd;
				firstAlignment=j;
			}
			j++;
		}
		
		// Adding alignments to the left of suffix substrings
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastSubstring) {
			if (!substrings[i].suffixReplication || substrings[i].prefixReplication || substrings[i].singleDeletionReplication || substrings[i].substringReplication) {
				i++;
				firstAlignment=-1;
				continue;
			}
			if (j==nAlignments || ReadA.sortedAlignments[j].endA()<substrings[i].endA-identityThreshold) {
				if (firstAlignment!=-1) {
					minStart=extendSubstrings_impl(firstAlignment,j-1,substrings[i],nAlignments);
					if (minStart<0) {
						substrings[i].leftAlignment=-1-minStart;
						periodicAlignmentChanged=true;
					}
					else substrings[i].leftAlignment=minStart;
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1; firstAlignment=-1;
				continue;
			}
			if (i<lastSubstring && firstJForNextI==-1 && ReadA.sortedAlignments[j].startA()<=substrings[i+1].endA) firstJForNextI=j;
			if (ReadA.sortedAlignments[j].endA()>substrings[i].endA+identityThreshold) {
				j++;
				continue;
			}
			if (firstAlignment==-1) firstAlignment=j;
			j++;
		}
		
		// Assigning new values to $startA$ and $endA$.
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].singleDeletionReplication) continue;
			if (substrings[i].prefixReplication && !substrings[i].suffixReplication && !substrings[i].substringReplication && substrings[i].leftAlignment!=-1) {
				substrings[i].endA=substrings[i].leftAlignment;
				substrings[i].maximalEndA=substrings[i].leftAlignment;
				substrings[i].sumEndA=substrings[i].endA;
				substrings[i].nEndA=1;
			}
			else if (substrings[i].suffixReplication && !substrings[i].prefixReplication && !substrings[i].substringReplication && substrings[i].leftAlignment!=-1) {
				substrings[i].startA=substrings[i].leftAlignment;
				substrings[i].maximalStartA=substrings[i].leftAlignment;
				substrings[i].sumStartA=substrings[i].startA;
				substrings[i].nStartA=1;
			}
			else if (substrings[i].substringReplication) {
				if (substrings[i].leftAlignment!=-1) {
					substrings[i].startA=substrings[i].leftAlignment;
					substrings[i].maximalStartA=substrings[i].leftAlignment;
					substrings[i].sumStartA=substrings[i].startA;
					substrings[i].nStartA=1;
				}
				if (substrings[i].rightAlignment!=-1) {
					substrings[i].endA=substrings[i].rightAlignment;
					substrings[i].maximalEndA=substrings[i].rightAlignment;
					substrings[i].sumEndA=substrings[i].endA;
					substrings[i].nEndA=1;
				}
			}
		}
		
		// Making sure that the range $[0..nAlignments]$ contains all and only the
		// alignments that are not implied by a periodic substring (this might have been
		// violated by $extendSubstrings_impl()$), and restoring its original order.
		j=nAlignments-1;
		if (periodicAlignmentChanged) {
			for (i=nAlignments; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) continue;
				j++;
				if (j!=i) {
					tmpAlignment=ReadA.sortedAlignments[j];
					ReadA.sortedAlignments[j]=ReadA.sortedAlignments[i];
					ReadA.sortedAlignments[i]=tmpAlignment;
				}
			}
		}
		Alignment.order=previousAlignmentOrder;
		if (j+1!=nAlignments) {
			if (j>0) Arrays.sort(ReadA.sortedAlignments,0,j+1);
			if (j+1<ReadA.lastSortedAlignment) Arrays.sort(ReadA.sortedAlignments,nAlignments,ReadA.lastSortedAlignment+1);
		}
		nAlignments=j+1;
		
		return nAlignments;
	}
	
	
	/**
	 * Let $substring$ be a prefix substring, let $ReadA.sortedAlignments[0..nAlignments
	 * -1]$ be sorted by startA, and let $[firstAlignment..lastAlignment] \subseteq 
	 * [0..nAlignments-1]$ be the maximal interval of $ReadA.sortedAlignments$ whose 
	 * starting positions are all close to the starting position of $substring$. 
	 * The procedure extends $substring$ by greedily merging maximal alignments whose last 
	 * position is iteratively reachable from the last position of $substring$, and such
	 * that they are not implied by a dense or periodic substring, and their 
	 * $inPeriodicSubstring$ field is set to true. Then, the procedure marks all 
	 * alignments that are implied by the extended substring, possibly reassigning to the
	 * extended substring some alignments that were previously assigned to a periodic 
	 * substring.
	 *
	 * A symmetrical argument holds if $substring$ is of suffix type: in this case, 
	 * $ReadA.sortedAlignments[firstAlignment..lastAlignment]$ is assumed to be sorted by 
	 * decreasing endA.
	 *
	 * Remark: the procedure uses the same distance threshold as $getSubstrings()$ for 
	 * extension.
	 *
	 * Remark: the procedure just sets fields $impliedBy{Prefix,Suffix}Substring$ of an 
	 * alignment, but it does not set field $impliedByDenseSubstring$.
	 *
	 * @param nAlignments the first $nAlignments$ of $ReadA.sortedAlignments$ are assumed
	 * not to be implied by periodic substrings. $ReadA.sortedAlignments[nAlignments..]$
	 * is assumed to be sorted by startA. After the procedure completes, some alignments
	 * in this range might have been reassigned from a periodic substring to $substring$.
	 * @return let X be the farthest coordinate of a newly-marked alignment; if an 
	 * alignment in $[nAlignments..]$ has been reassigned from a periodic substring to 
	 * $substring$, the procedure returns -1-X; otherwise, the procedure returns X.
	 */
	private static final int extendSubstrings_impl(int firstAlignment, int lastAlignment, DenseSubstring substring, int nAlignments) {
		final boolean prefixOrSuffix = substring.prefixReplication;
		boolean periodicAlignmentChanged;
		int i;
		int top, alignmentID, alignmentEnd, out;
		
		// Extending substring
		top=0; stack[0]=-1; out=prefixOrSuffix?substring.endA:substring.startA;
		while (top>=0) {
			alignmentID=stack[top--];
			alignmentEnd=alignmentID==-1?(prefixOrSuffix?substring.endA:substring.startA):(prefixOrSuffix?ReadA.sortedAlignments[alignmentID].endA():ReadA.sortedAlignments[alignmentID].startA());
			for (i=firstAlignment; i<=lastAlignment; i++) {
				if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
					 ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null || 
					 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null || 
					 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null || 
					 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
					 !ReadA.sortedAlignments[i].inPeriodicSubstring ||
					 ( prefixOrSuffix && 
					   ( ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
						 ReadA.sortedAlignments[i].endA()<alignmentEnd-maxDistanceEnd ||
					     ReadA.sortedAlignments[i].endA()>alignmentEnd+maxDistanceEnd
					   )
					 ) ||
					 ( !prefixOrSuffix && 
					   ( ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
						 ReadA.sortedAlignments[i].startA()>alignmentEnd+maxDistanceEnd ||
					     ReadA.sortedAlignments[i].startA()<alignmentEnd-maxDistanceEnd
					   )
					 )
				   ) continue;
				if (prefixOrSuffix) {
					ReadA.sortedAlignments[i].impliedByPrefixSubstring=substring;
					out=Math.max(out,ReadA.sortedAlignments[i].endA());
				}
				else {
					ReadA.sortedAlignments[i].impliedBySuffixSubstring=substring;
					out=Math.min(out,ReadA.sortedAlignments[i].startA());
				}
				stack[++top]=i;
			}
		}
		
		// Marking alignments implied by the extended substring and not already implied
		// by a periodic substring.
		i=firstAlignment;
		while (i<nAlignments) {
			if ( (prefixOrSuffix && ReadA.sortedAlignments[i].startA()>=substring.endA) ||
				 (!prefixOrSuffix && ReadA.sortedAlignments[i].endA()<=substring.startA)
			   ) break;
			if ( ReadA.sortedAlignments[i].inDenseSubstring!=null ||
				 ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null || 
				 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null || 
				 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null || 
				 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
				 !ReadA.sortedAlignments[i].inPeriodicSubstring ||
				 ( prefixOrSuffix &&
				   ( Math.abs(ReadA.sortedAlignments[i].startA(),substring.startA)>identityThreshold ||
					 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
					 ReadA.sortedAlignments[i].endA()<substring.startA+substring.minPrefixLength-identityThreshold ||
				     ReadA.sortedAlignments[i].endA()>out
				   )
				 ) ||
				 ( !prefixOrSuffix && 
				   ( Math.abs(ReadA.sortedAlignments[i].endA(),substring.endA)>identityThreshold ||
					 ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
					 ReadA.sortedAlignments[i].startA()>substring.endA-substring.minSuffixLength+identityThreshold ||
					 ReadA.sortedAlignments[i].startA()<out
				   )
				 )
			   ) {
				 i++;
				 continue;
			}
			if (prefixOrSuffix) ReadA.sortedAlignments[i].impliedByPrefixSubstring=substring;
			else ReadA.sortedAlignments[i].impliedBySuffixSubstring=substring;
			i++;
		}
		
		// Marking alignments implied by the extended substring and already implied by a
		// periodic substring.
		i=nAlignments; periodicAlignmentChanged=false;
		while (i<=ReadA.lastSortedAlignment) {
			if ( ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null ||   
				 ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null ||
				 ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null ||
				 ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null ||
				 ReadA.sortedAlignments[i].impliedByDenseSubstring!=null ||
				 ReadA.sortedAlignments[i].inDenseSubstring!=null
			   ) {
				i++;
			 	continue;
			}
			if ( ReadA.sortedAlignments[i].periodicSubstringInterval!=null && 
				 Intervals.areApproximatelyIdentical(ReadA.sortedAlignments[i].startA(),ReadA.sortedAlignments[i].endA(),ReadA.sortedAlignments[i].periodicSubstringInterval.firstPosition,ReadA.sortedAlignments[i].periodicSubstringInterval.lastPosition)
			   ) {
				i++;
				continue;
			}
			if (prefixOrSuffix) {
				if (ReadA.sortedAlignments[i].startA()>substring.startA+identityThreshold) break;
			   	if ( ReadA.sortedAlignments[i].startA()<substring.startA-identityThreshold ||
					 ReadA.sortedAlignments[i].isRightMaximalB!=1 ||
				     ReadA.sortedAlignments[i].endA()<substring.startA+substring.minPrefixLength-identityThreshold ||
			         ReadA.sortedAlignments[i].endA()>out
			       ) {
					i++;
					continue;
			    }
				ReadA.sortedAlignments[i].impliedByPrefixSubstring=substring;
				ReadA.sortedAlignments[i].impliedByPeriodicSubstring=null;
				ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				periodicAlignmentChanged=true;
				i++;
			}
			else {
				if ( Math.abs(ReadA.sortedAlignments[i].endA(),substring.endA)>identityThreshold ||
					 ReadA.sortedAlignments[i].isLeftMaximalB!=1 ||
					 ReadA.sortedAlignments[i].startA()>substring.endA-substring.minSuffixLength+identityThreshold ||
					 ReadA.sortedAlignments[i].startA()<out
				   ) {
					i++;
					continue;
				}
				ReadA.sortedAlignments[i].impliedBySuffixSubstring=substring;
				ReadA.sortedAlignments[i].impliedByPeriodicSubstring=null;
				ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				periodicAlignmentChanged=true;
				i++;
		    }
		}
		
		return periodicAlignmentChanged?-1-out:out;
	}
	
	
	/**
	 * Builds maximal intervals of pairwise-straddling weak dense substrings, of any type.
	 * For every such interval that contains at least one weak dense substring of 
	 * substring type, the procedure keeps just one of the dense substrings of substring 
	 * type in the interval, it resets its boundaries to those of the interval, and it
	 * discards all other weak dense substrings. This is useful since, in practice, 
	 * regions with multiple straddling weak dense substrings do not show clear boundaries
	 * between the substrings.
	 *
	 * Remark: the procedure uses also non-weak dense sustrings of substring type, if 
	 * their density of maximal events is compatible with other weak susbtrings. 
	 * Overlapping non-weak dense substrings of substring type have already been merged by
	 * $cleanDenseSubstringsOfSubstringType()$.
	 *
	 * Remark: the procedure merges also a set of substrings for which there is no 
	 * straddling, i.e. that are just contained in one another, iff: (1) the longest 
	 * substring of the set if not of substring type; (2) there is a long other substring
	 * of substring type.
	 *
	 * Remark: substrings with very different average coverage (computed using only their 
	 * implied alignments) are not concatenated. Substrings whose boundary is approx.
	 * identical to the boundary of a periodic interval, or of a non-weak substring of
	 * prefix/suffix type, are not concatenated: existing strong boundaries are used as a 
	 * proxy for peaks of maximal events, and we don't concatenate substrings if their 
	 * boundaries are peaks of maximal events.
	 *
	 * Remark: the procedure assumes $substrings$ to be sorted by $startA$.
	 *
	 * @param tmp temporary space, of size at least equal to $3*(lastSubstring+1)*2$;
	 * @return the last element of $tmp$ used, if a concatenation was performed; -1
	 * otherwise.
	 */
	public static final int concatenateWeakSubstrings(int lastAlignment, int[] tmp) {
		final int MIN_INTERSECTION_LENGTH = IO.quantum<<1;  // Arbitrary
		final int DENSITY_THRESHOLD = 5;  // Arbitrary
		final double WINDOWS_DENSITY_THRESHOLD = 2.0;  // Arbitrary
		final double CONTAINMENT_RATIO = 0.7;  // Arbitrary
		final double COVERAGE_RATIO = 2.0;  // Arbitrary
		final int COVERAGE_WINDOW = IO.quantum;  // Arbitrary
		final int MIN_COVERAGE = IO.coverage<<2;  // Arbitrary
		final int LARGE_WINDOW = IO.coverage*10;  // Arbitrary
		boolean atLeastOneSubstringReplication, atLeastOneStraddling;
		int i, j, k, p, q;
		int startA1, endA1, startA2, endA2, nSubstrings, nMaximalAlignments, maxWindow, newMark;
		int substringReplication, lastInterval, representative, firstKForNextI, previousOrder;
		int length, longest, longestSubstring, longestLength, longestSubstringLength;
		int surface, coverage, rangeCoverage;
		double densityI, densityJ;
		DenseSubstring tmpSubstring;
		PeriodicSubstringInterval tmpPeriodic = new PeriodicSubstringInterval();

		Factorize.computeMaximalAlignments(true,false,false,lastAlignment);
		computeDensityWindows(lastAlignment);
		
		// Computing maximal intervals of straddling weak substrings
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].previousSumStartA=0;  // Used just as a flag
		}
		k=-3; i=0; atLeastOneSubstringReplication=false;
		while (i<=lastSubstring) {
			if (substrings[i].previousSumStartA==2 || (!substrings[i].isWeak && !substrings[i].substringReplication) || substrings[i].rightWindow>=LARGE_WINDOW) {
				i++;
				continue;
			}
			k+=3;
			tmp[k]=substrings[i].maximalStartA>=0?substrings[i].maximalStartA:substrings[i].startA;
			tmp[k+1]=substrings[i].maximalEndA>0?substrings[i].maximalEndA:substrings[i].endA;
			substringReplication=substrings[i].substringReplication?i:-1;
			nMaximalAlignments=substrings[i].nMaximalAlignments;
			maxWindow=substrings[i].maxWindow;
			longest=i; longestLength=substrings[i].length();
			longestSubstring=substrings[i].substringReplication?i:-1;
			longestSubstringLength=substrings[i].substringReplication?substrings[i].length():0;
			surface=substrings[i].surface;
			substrings[i].previousSumStartA=1;
			j=i+1;
			while (j<=lastSubstring) {
				if (substrings[j].previousSumStartA==2 || (!substrings[j].isWeak && !substrings[j].substringReplication) || substrings[i].leftWindow>=LARGE_WINDOW) {
					j++;
					continue;
				}
				startA2=substrings[j].maximalStartA>=0?substrings[j].maximalStartA:substrings[j].startA;
				endA2=substrings[j].maximalEndA>0?substrings[j].maximalEndA:substrings[j].endA;
				if (startA2>tmp[k+1]-MIN_INTERSECTION_LENGTH) break;
				densityI=((double)nMaximalAlignments)/Math.max(tmp[k+1]-tmp[k]+1-Alignments.minAlignmentLength,1);
				densityJ=((double)substrings[j].nMaximalAlignments)/Math.max(substrings[j].length()-Alignments.minAlignmentLength,1);
				if ( !substrings[j].isWeak &&
				     ( Math.max(densityI,densityJ)>Math.min(densityI,densityJ)*DENSITY_THRESHOLD || 
					   (startA2>tmp[k]+identityThreshold && substrings[j].leftWindow>=maxWindow*WINDOWS_DENSITY_THRESHOLD) ||
					   (endA2<tmp[k+1]-identityThreshold && substrings[j].rightWindow>=maxWindow*WINDOWS_DENSITY_THRESHOLD)
					 )
				   ) {
					j++;
					continue;
				}
				coverage=substrings[j].surface/substrings[j].length();
				rangeCoverage=surface/(tmp[k+1]-tmp[k]+1);
				if (Math.max(coverage,rangeCoverage)>=Math.min(coverage,rangeCoverage)*COVERAGE_RATIO) {
					j++;
					continue;
				}
				if ( isBoundaryOfPeriodicInterval(startA2,true,tmpPeriodic) || isBoundaryOfSubstring(startA2,true,j) ||
					 isBoundaryOfPeriodicInterval(tmp[k+1],false,tmpPeriodic) || isBoundaryOfSubstring(tmp[k+1],false,j)
				   ) {
					j++;
					continue;
				}
				if (ReadA.getAverageCoverage(tmp[k+1]+1,Math.min(tmp[k+1]+COVERAGE_WINDOW-1,ReadA.readLength-1))<MIN_COVERAGE) {
					j++;
					continue;
				}
				if (ReadA.getAverageCoverage(Math.max(startA2-COVERAGE_WINDOW+1,0),startA2-1)<MIN_COVERAGE) {
					j++;
					continue;
				}
				if (ReadA.getAverageCoverage(startA2,tmp[k+1])<MIN_COVERAGE) {
					j++;
					continue;
				}
				tmp[k+1]=Math.max(tmp[k+1],endA2);
				if (substringReplication==-1 && substrings[j].substringReplication) substringReplication=j;
				nMaximalAlignments+=substrings[j].nMaximalAlignments;
				maxWindow=Math.max(maxWindow,substrings[j].maxWindow);
				length=substrings[j].length();
				if (length>longestLength) {
					longest=j;
					longestLength=length;
				}
				if (substrings[j].substringReplication && length>longestSubstringLength) {
					longestSubstring=j;
					longestSubstringLength=length;
				}
				surface+=substrings[j].surface;
				substrings[j].previousSumStartA=1;
				j++;
			}
			// Discarding the interval if it contains two adjacent substrings.
			if (substringReplication!=-1) {
				for (p=i; p<j; p++) {
					if (substrings[p].previousSumStartA!=1) continue;
					endA1=substrings[p].maximalEndA>0?substrings[p].maximalEndA:substrings[p].endA;
					for (q=p+1; q<j; q++) {
						if (substrings[q].startA>=substrings[p].endA) break;
						if (substrings[q].previousSumStartA!=1) continue;
						startA2=substrings[q].maximalStartA>=0?substrings[q].maximalStartA:substrings[q].startA;
						if (Math.abs(startA2,endA1)<=identityThreshold) {
							substringReplication=-1;
							break;
						}
					}
					if (substringReplication==-1) break;
				}
			}
			// Discarding the interval if there is no straddling.
			if (substringReplication!=-1) {
				atLeastOneStraddling=false;
				for (p=i; p<j; p++) {
					if (substrings[p].previousSumStartA!=1) continue;
					startA1=substrings[p].startA; endA1=substrings[p].endA;
					for (q=p+1; q<j; q++) {
						if (substrings[q].startA>=endA1) break;
						if (substrings[q].previousSumStartA!=1) continue;
						if (Intervals.straddles(startA1,endA1,substrings[q].startA,substrings[q].endA)) {
							atLeastOneStraddling=true;
							break;
						}
					}
					if (atLeastOneStraddling) break;
				}
				if ( !atLeastOneStraddling &&
					 !(!substrings[longest].substringReplication && longestSubstring>=0 && longestSubstringLength>=longestLength*CONTAINMENT_RATIO)
				   ) substringReplication=-1;
			}
			tmp[k+2]=substringReplication;
			if (substringReplication!=-1) atLeastOneSubstringReplication=true;
			// Resetting $previousSumStartA$ marks.
			newMark=substringReplication>=0?2:0;
			j=i;
			while (j<=lastSubstring) {
				if (substrings[j].previousSumStartA==2 || (!substrings[j].isWeak && !substrings[j].substringReplication)) {
					j++;
					continue;
				}
				if (substrings[j].startA>tmp[k+1]-MIN_INTERSECTION_LENGTH) break;
				if (substrings[j].previousSumStartA==1) substrings[j].previousSumStartA=newMark;
				j++;
			}
			if (j==i) break;
			i=j;
		}
		lastInterval=k;
		if (lastInterval==-3 || !atLeastOneSubstringReplication) return -1;
		
		// Updating weak substrings
		for (i=0; i<=lastSubstring; i++) {
			if (!substrings[i].isWeak && !substrings[i].substringReplication) continue;
			substrings[i].discarded=false;
			substrings[i].representative=null;
		}
		Math.set(tmp,lastInterval+3,(lastInterval+3)*2-1,0);
		for (k=0; k<=lastInterval; k+=3) {
			i=tmp[k+2];
			if (i==-1) continue;
			tmp[lastInterval+3+k]=substrings[i].startA;
			substrings[i].startA=tmp[k];
			tmp[lastInterval+3+k+1]=substrings[i].endA;
			substrings[i].endA=tmp[k+1];
			substrings[i].sumStartA=tmp[k];
			substrings[i].nStartA=1;
			substrings[i].sumEndA=tmp[k+1];
			substrings[i].nEndA=1;
			substrings[i].isWeak=true;
		}
		k=0; i=0; firstKForNextI=-1; representative=-1;
		while (i<=lastSubstring) {
			if (k>lastInterval || tmp[k]>substrings[i].endA || (!substrings[i].isWeak && !substrings[i].substringReplication) || substrings[i].previousSumStartA!=2) {
				if (representative!=-1) {
					substrings[i].discarded=true;
					substrings[i].representative=substrings[representative];
					tmpSubstring=substrings[representative];
					if (substrings[i].maximalStartA>=0) tmpSubstring.maximalStartA=tmpSubstring.maximalStartA>=0?Math.min(tmpSubstring.maximalStartA,substrings[i].maximalStartA):substrings[i].maximalStartA;
					if (substrings[i].maximalEndA>0) tmpSubstring.maximalEndA=tmpSubstring.maximalEndA>0?Math.max(tmpSubstring.maximalEndA,substrings[i].maximalEndA):substrings[i].maximalEndA;
				}
				i++;
				if (firstKForNextI>=0) k=firstKForNextI;
				firstKForNextI=-1;
				representative=-1;
				continue;
			}
			if (tmp[k+1]<substrings[i].startA || tmp[k+2]==-1) {
				k+=3;
				continue;
			}
			if (firstKForNextI==-1 && i<lastSubstring && tmp[k+1]>=substrings[i+1].startA) firstKForNextI=k;
			if (i!=tmp[k+2] && Intervals.isContained(substrings[i].maximalStartA>=0?substrings[i].maximalStartA:substrings[i].startA,substrings[i].maximalEndA>0?substrings[i].maximalEndA:substrings[i].endA,tmp[k],tmp[k+1])) representative=tmp[k+2];
			k+=3;
		}
		
		for (k=0; k<=lastInterval; k+=3) {
			i=tmp[k+2];
			if (i==-1) continue;
			substrings[i].startA=Math.min(substrings[i].startA,tmp[lastInterval+3+k]);
			substrings[i].endA=Math.max(substrings[i].endA,tmp[lastInterval+3+k+1]);
		}
		for (i=0; i<=lastSubstring; i++) {
			if ((substrings[i].isWeak || substrings[i].substringReplication) && substrings[i].discarded) {
				substrings[i].representative.startA=Math.min(substrings[i].representative.startA,substrings[i].startA);
				substrings[i].representative.endA=Math.max(substrings[i].representative.endA,substrings[i].endA);
			}
		}
		j=-1;
		for (i=0; i<=lastSubstring; i++) {
			if ((substrings[i].isWeak || substrings[i].substringReplication) && substrings[i].discarded) continue;
			j++;
			if (j!=i) {
				tmpSubstring=substrings[j];
				substrings[j]=substrings[i];
				substrings[i]=tmpSubstring;
			}
		}
		lastSubstring=j;
		
		// Updating pointers from alignments
		for (i=0; i<=lastAlignment; i++) {
			tmpSubstring=ReadA.sortedAlignments[i].inDenseSubstring;
			if (tmpSubstring!=null && (tmpSubstring.isWeak || tmpSubstring.substringReplication) && tmpSubstring.discarded) {
				ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring.representative;
			}
			tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (tmpSubstring!=null && (tmpSubstring.isWeak || tmpSubstring.substringReplication) && tmpSubstring.discarded) {
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring.representative;
			}
		}
		
		return 3*(lastInterval+1)-1;
	}
	
	
	/**
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ to be sorted by 
	 * firstPosition, and to contain both short- and long-period intervals.
	 *
	 * @param tmpInterval temporary space;
	 * @return TRUE iff $pos$ is approximately identical to the first (if 
	 * $firstOrLast=true$) or last (if $firstOrLast=false$) position of a periodic 
	 * interval.
	 */
	private static final boolean isBoundaryOfPeriodicInterval(int pos, boolean firstOrLast, PeriodicSubstringInterval tmpInterval) {
		int i, j;
		if (PeriodicSubstrings.lastInterval==-1) return false;
		
		tmpInterval.firstPosition=pos;
		i=Arrays.binarySearch(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1,tmpInterval);
		if (firstOrLast) {
			if (i>=0) return true;
			i=-i-1;
			if ( (i<=PeriodicSubstrings.lastInterval && PeriodicSubstrings.intervals[i].firstPosition<=pos+identityThreshold) ||
				 (i-1>=0 && PeriodicSubstrings.intervals[i-1].firstPosition>=pos-identityThreshold)
			   ) return true;
		}
		else {
			if (i<0) i=-i-1;
			for (j=i-1; j>=0; j--) {
				if (Math.abs(PeriodicSubstrings.intervals[j].lastPosition,pos)<=identityThreshold) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Remark: the procedure assumes $substrings$ to be sorted by first position.
	 *
	 * @param index position in $substrings$ from which the search starts;
	 * @return TRUE iff $pos$ is approximately identical to the first (if 
	 * $firstOrLast=true$) or last (if $firstOrLast=false$) position of a non-weak 
	 * substring of prefix (if $firstOrLast=true$), suffix (if $firstOrLast=false$),
	 * prefix+suffix, or single deletion type.
	 */
	private static final boolean isBoundaryOfSubstring(int pos, boolean firstOrLast, int index) {
		int i, j;
		if (lastSubstring==-1) return false;
		
		if (firstOrLast) {
			for (j=index-1; j>=0; j--) {
				if (substrings[j].startA<pos-identityThreshold) break;
				if (substrings[j].startA>pos+identityThreshold) continue;
				if (!substrings[j].isWeak && (substrings[j].prefixReplication || substrings[j].singleDeletionReplication)) return true;
			}
			for (j=index; j<=lastSubstring; j++) {
				if (substrings[j].startA>pos+identityThreshold) break;
				if (substrings[j].startA<pos-identityThreshold) continue;
				if (!substrings[j].isWeak && (substrings[j].prefixReplication || substrings[j].singleDeletionReplication)) return true;
			}
		}
		else {
			for (j=index-1; j>=0; j--) {
				if (substrings[j].startA>=pos) continue;
				if (!substrings[j].isWeak && (substrings[j].suffixReplication || substrings[j].singleDeletionReplication) && Math.abs(substrings[j].endA,pos)<=identityThreshold) return true;
			}
			for (j=index; j<=lastSubstring; j++) {
				if (substrings[j].startA>=pos) break;
				if (!substrings[j].isWeak && (substrings[j].suffixReplication || substrings[j].singleDeletionReplication) && Math.abs(substrings[j].endA,pos)<=identityThreshold) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * It could happen that a peak in the alignment endpoints of a prefix substring is not 
	 * detected because it is too close to the end of the substring, and that the prefix 
	 * substring is then concatenated to another one, forming a longer prefix substring
	 * in which the peak is far from the endpoint and can be detected. The procedure tries
	 * to find peaks in prefix substrings after concatenation, using just the alignments 
	 * that are already assigned to the prefix substring, and removing the assignment of
	 * those that end at peaks.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments[0..lastAlignment]$ to be
	 * sorted by startA, and it assumes all alignments implied by periodic
	 * substrings to be located in $ReadA.sortedAlignments[lastAlignment+1..]$.
	 * The procedure assumes also that all and only the substrings removed by 
	 * concatenation are located in $substrings[lastSubstring+1..previousLastSubstring]$.
	 *
	 * Remark: the procedure uses fields $previousSum*A,previousN*A$ of representative 
	 * dense substrings, which are assumed to have been set by procedure 
	 * $getConnectedComponent_concatenation()$, to handle the fact that a representative 
	 * was not necessarily the longest substring in its set before concatenation.
	 */
	private static final void findPeaksAfterConcatenation(int lastAlignment, int previousLastSubstring, int[] tmpArray) {
		final int THRESHOLD = identityThreshold;
		int i, j;
		int candidate, lastCandidate, startA, endA;
		
		for (i=0; i<=lastSubstring; i++) {
			if (!substrings[i].isConcatenated) continue;
			
			// Collecting candidate positions for new peaks
			startA=substrings[i].startA;
			endA=substrings[i].endA;
			lastCandidate=-1;
			candidate=substrings[i].prefixReplication?substrings[i].previousSumEndA/substrings[i].previousNEndA:substrings[i].previousSumStartA/substrings[i].previousNStartA;
			if (candidate>startA+THRESHOLD && candidate<endA-THRESHOLD) tmpArray[++lastCandidate]=candidate;
			for (j=lastSubstring+1; j<=previousLastSubstring; j++) {
				if (substrings[j].isConcatenated && substrings[j].representative==substrings[i]) {
					candidate=substrings[j].prefixReplication?substrings[j].sumEndA/substrings[j].nEndA:substrings[j].sumStartA/substrings[j].nStartA;
					if (candidate>startA+THRESHOLD && candidate<endA-THRESHOLD) tmpArray[++lastCandidate]=candidate;
				}
			}
			if (lastCandidate==-1) continue;
			if (lastCandidate>0) Arrays.sort(tmpArray,0,lastCandidate+1);
			
			// Detecting peaks
			if (substrings[i].prefixReplication) findPeaksAfterConcatenation_prefix(substrings[i],lastAlignment,tmpArray,lastCandidate);
			if (substrings[i].suffixReplication) findPeaksAfterConcatenation_suffix(substrings[i],lastAlignment,tmpArray,lastCandidate);
		}
		if (lastLengthBiasInterval>0) sortLengthBiasIntervals();
	}
	
	
	/**
	 * A simplified version of what is done in $markImpliedAlignmentsAggressive*$ 
	 * procedures. It only looks for peaks close to $candidates[0..lastCandidate]$.
	 *
	 * @param substring a prefix substring that has changed after concatenation.
	 */
	private static final void findPeaksAfterConcatenation_prefix(DenseSubstring substring, int lastAlignment, int[] candidates, int lastCandidate) {
		final int MIN_BANDWIDTH = (identityThreshold<<1)/3;
		final int MAX_BANDWIDTH = identityThreshold;
		final int MIN_COUNT = 3;  // Min. n. of events to be used in a $\delta$ histogram
		final int THRESHOLD = identityThreshold;
		final int NEW_BOUNDARIES_THRESHOLD = (THRESHOLD)<<1;  // Arbitrary
		final double UNASSIGNED_THRESHOLD = 0.8;  // Arbitrary
		final int substringLength = substring.endA-substring.startA+1;
		final int prefixOrigin = substring.startA;
		int i, j, k;
		int from, to, pos, id, length, startA, endA, lastDeltaPointPrefix, bandwidth;
		int newSumStartA, newNStartA, newEndA, newMinPrefixLength, first, last;
		int nUnassigned, nAlignments;
		
		pos=ReadA.binarySearch(substring.startA,0,lastAlignment);
		if (pos<0) pos=-pos-1;

		// Filling $deltaPointsPrefix$
		Math.set(countsPrefix,substringLength-1,0);
		lastDeltaPointPrefix=-1; first=-1; last=-1;
		for (i=pos; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA<=substring.startA-THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			first=i;
			if (last==-1) last=i;
			lastDeltaPointPrefix++;
			endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			deltaPointsPrefix[lastDeltaPointPrefix].position=endA-prefixOrigin;
			deltaPointsPrefix[lastDeltaPointPrefix].mass=1;
			for (j=0; j<=endA-prefixOrigin; j++) countsPrefix[j]++;
		}
		for (i=pos+1; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA>=substring.endA) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			last=i;
			if (first==-1) first=i;
			lastDeltaPointPrefix++;
			endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			deltaPointsPrefix[lastDeltaPointPrefix].position=endA-prefixOrigin;
			deltaPointsPrefix[lastDeltaPointPrefix].mass=1;
			for (j=0; j<=endA-prefixOrigin; j++) countsPrefix[j]++;
		}
		if (lastDeltaPointPrefix+1<MIN_ALIGNMENTS_FOR_LENGTH_BIAS) return;
	
		// Computing $deltaPrefix$ and finding peaks
		lastDeltaPointPrefix=Points.sortAndCompact(deltaPointsPrefix,lastDeltaPointPrefix);	
		bandwidth=Points.estimateDensity(deltaPointsPrefix,lastDeltaPointPrefix,deltaPrefix,MIN_BANDWIDTH,MAX_BANDWIDTH);
		for (i=0; i<=lastDeltaPointPrefix; i++) {
			from=(int)deltaPointsPrefix[i].position+1;
			to=Math.min((int)deltaPointsPrefix[i].position+bandwidth,countsPrefix.length-1);
			for (j=from; j<=to; j++) countsPrefix[j]+=deltaPointsPrefix[i].getMass();
		}
		for (i=0; i<substringLength; i++) {
			if (countsPrefix[i]>=MIN_COUNT) deltaPrefix[i]/=countsPrefix[i];
			else deltaPrefix[i]=0;
		}
		estimateFromDelta(deltaPrefix,substringLength,substring,true,first,last,0,false,-1);
		if (lastLengthBoundaryPrefix<0) return;
		lastLengthBoundaryPrefix=cleanBoundaries(lengthBoundariesPrefix,lastLengthBoundaryPrefix,candidates,lastCandidate,THRESHOLD);
		if (lastLengthBoundaryPrefix<0) return;

		// Stopping if most alignments would get unmarked
		nUnassigned=0; nAlignments=0;
		for (i=pos; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (startA<=substring.startA-THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			nAlignments++;
			length=endA-prefixOrigin;
			j=Arrays.binarySearch(lengthBoundariesPrefix,0,lastLengthBoundaryPrefix+1,length);
			if ( j>=0 || 
				 (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundaryPrefix &&
				   ( Math.abs(length,lengthBoundariesPrefix[-j])<=THRESHOLD || 
				     (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=THRESHOLD)
				   )
				 )		 
			   ) nUnassigned++;
		}
		for (i=pos+1; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (startA>=substring.endA) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			nAlignments++;
			length=endA-prefixOrigin;
			j=Arrays.binarySearch(lengthBoundariesPrefix,0,lastLengthBoundaryPrefix+1,length);
			if ( j>=0 || 
				 (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundaryPrefix &&
				   ( Math.abs(length,lengthBoundariesPrefix[-j])<=THRESHOLD || 
				     (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=THRESHOLD)
				   )
				 )
			   ) nUnassigned++;
		}
		if (nUnassigned>=nAlignments*UNASSIGNED_THRESHOLD) return;
	
		// Unmarking alignments
		for (i=pos; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (startA<=substring.startA-THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			length=endA-prefixOrigin;
			j=Arrays.binarySearch(lengthBoundariesPrefix,0,lastLengthBoundaryPrefix+1,length);
			if ( j>=0 || 
				 (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundaryPrefix &&
				   ( Math.abs(length,lengthBoundariesPrefix[-j])<=THRESHOLD || 
				     (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=THRESHOLD)
				   )
				 )
			   ) ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
		}
		for (i=pos+1; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (startA>=substring.endA) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			length=endA-prefixOrigin;
			j=Arrays.binarySearch(lengthBoundariesPrefix,0,lastLengthBoundaryPrefix+1,length);
			if ( j>=0 || 
				 (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundaryPrefix &&
				   ( Math.abs(length,lengthBoundariesPrefix[-j])<=THRESHOLD || 
				     (-j-1==0?false:Math.abs(length,lengthBoundariesPrefix[-j-3])<=THRESHOLD)
				   )
				 )
			   ) ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
		}
		
		// Updating properties of $substring$
		newSumStartA=0; newNStartA=0; newMinPrefixLength=substringLength;
		for (i=pos; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA<substring.startA-NEW_BOUNDARIES_THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring || startA>substring.startA+NEW_BOUNDARIES_THRESHOLD) continue;
			newSumStartA+=startA; newNStartA++;
			endA=Alignments.alignments[id][4];
			length=endA-prefixOrigin+1;
			if (length<newMinPrefixLength) newMinPrefixLength=length;
		}
		for (i=pos+1; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA>substring.startA+NEW_BOUNDARIES_THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			newSumStartA+=startA; newNStartA++;
			endA=Alignments.alignments[id][4];
			length=endA-prefixOrigin+1;
			if (length<newMinPrefixLength) newMinPrefixLength=length;
		}
		substring.sumStartA=newSumStartA;
		substring.nStartA=newNStartA;
		if (newNStartA>0) substring.startA=newSumStartA/newNStartA;
		substring.minPrefixLength=newMinPrefixLength;
		substring.hasLengthBoundariesPrefix=true;
	
		// Updating $lengthBiasIntervals$
		for (i=0; i<lastLengthBoundaryPrefix; i+=3) {
			lastLengthBiasInterval++;
			lengthBiasIntervals[lastLengthBiasInterval].firstPosition=prefixOrigin+lengthBoundariesPrefix[i];
			lengthBiasIntervals[lastLengthBiasInterval].center=prefixOrigin+lengthBoundariesPrefix[i+1];
			lengthBiasIntervals[lastLengthBiasInterval].lastPosition=prefixOrigin+lengthBoundariesPrefix[i+2];
			lengthBiasIntervals[lastLengthBiasInterval].readB=ReadA.id;
			lengthBiasIntervals[lastLengthBiasInterval].orientation=true;
		}
	}
	
	
	/**
	 * Symmetrical to $findPeaksAfterConcatenation_prefix()$.
	 *
	 * @param substring a suffix substring that has changed after concatenation.
	 */
	private static final void findPeaksAfterConcatenation_suffix(DenseSubstring substring, int lastAlignment, int[] candidates, int lastCandidate) {
		final int MIN_BANDWIDTH = (identityThreshold<<1)/3;
		final int MAX_BANDWIDTH = identityThreshold;
		final int MIN_COUNT = 3;  // Min. n. of events to be used in a $\delta$ histogram
		final int THRESHOLD = identityThreshold;
		final int NEW_BOUNDARIES_THRESHOLD = (THRESHOLD)<<1;  // Arbitrary
		final double UNASSIGNED_THRESHOLD = 0.8;  // Arbitrary
		final int substringLength = substring.endA-substring.startA+1;
		final int suffixOrigin = substring.endA;
		int i, j;
		int from, to, pos, id, length, startA, endA, lastDeltaPointSuffix, bandwidth;
		int newSumEndA, newNEndA, newStartA, newMinSuffixLength, first, last;
		int nUnassigned, nAlignments;
		
		pos=ReadA.binarySearch(substring.startA,0,lastAlignment);
		if (pos<0) pos=-pos-1;

		// Filling $deltaPointsSuffix$
		Math.set(countsSuffix,substringLength-1,0);
		lastDeltaPointSuffix=-1; first=-1; last=-1;
		for (i=pos; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA<=substring.startA-THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			first=i;
			if (last==-1) last=i;
			lastDeltaPointSuffix++;
			deltaPointsSuffix[lastDeltaPointSuffix].position=suffixOrigin-startA;
			deltaPointsSuffix[lastDeltaPointSuffix].mass=1;
			for (j=0; j<=suffixOrigin-startA; j++) countsSuffix[j]++;
		}
		for (i=pos+1; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA>=substring.endA) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			last=i;
			if (first==-1) first=i;
			lastDeltaPointSuffix++;
			deltaPointsSuffix[lastDeltaPointSuffix].position=suffixOrigin-startA;
			deltaPointsSuffix[lastDeltaPointSuffix].mass=1;
			for (j=0; j<=suffixOrigin-startA; j++) countsSuffix[j]++;
		}
		if (lastDeltaPointSuffix+1<MIN_ALIGNMENTS_FOR_LENGTH_BIAS) return;
	
		// Computing $deltaSuffix$ and finding peaks
		lastDeltaPointSuffix=Points.sortAndCompact(deltaPointsSuffix,lastDeltaPointSuffix);	
		bandwidth=Points.estimateDensity(deltaPointsSuffix,lastDeltaPointSuffix,deltaSuffix,MIN_BANDWIDTH,MAX_BANDWIDTH);
		for (i=0; i<=lastDeltaPointSuffix; i++) {
			from=(int)deltaPointsSuffix[i].position+1;
			to=Math.min((int)deltaPointsSuffix[i].position+bandwidth,countsSuffix.length-1);
			for (j=from; j<=to; j++) countsSuffix[j]+=deltaPointsSuffix[i].getMass();
		}
		for (i=0; i<substringLength; i++) {
			if (countsSuffix[i]>=MIN_COUNT) deltaSuffix[i]/=countsSuffix[i];
			else deltaSuffix[i]=0;
		}
		estimateFromDelta(deltaSuffix,substringLength,substring,true,first,last,0,false,-1);
		if (lastLengthBoundarySuffix<0) return;
		lastLengthBoundarySuffix=cleanBoundaries(lengthBoundariesSuffix,lastLengthBoundarySuffix,candidates,lastCandidate,THRESHOLD);
		if (lastLengthBoundarySuffix<0) return;
		
		// Stopping if most alignments would get unmarked
		nUnassigned=0; nAlignments=0;
		for (i=pos; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (startA<substring.startA-THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			nAlignments++;
			length=suffixOrigin-startA;
			j=Arrays.binarySearch(lengthBoundariesSuffix,0,lastLengthBoundarySuffix+1,length);
			if ( j>=0 || 
				 (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundarySuffix &&
				   ( Math.abs(length,lengthBoundariesSuffix[-j])<=THRESHOLD || 
				     (-j-1==0?false:Math.abs(length,lengthBoundariesSuffix[-j-3])<=THRESHOLD)
				   )
				 )
			   ) nUnassigned++;
		}
		for (i=pos+1; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (startA>=substring.endA) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			nAlignments++;
			length=suffixOrigin-startA;
			j=Arrays.binarySearch(lengthBoundariesSuffix,0,lastLengthBoundarySuffix+1,length);
			if ( j>=0 || 
				 (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundarySuffix &&
				   ( Math.abs(length,lengthBoundariesSuffix[-j])<=THRESHOLD || 
				     (-j-1==0?false:Math.abs(length,lengthBoundariesSuffix[-j-3])<=THRESHOLD)
				   )
				 )
			   ) nUnassigned++;
		}
		if (nUnassigned>=nAlignments*UNASSIGNED_THRESHOLD) return;
	
		// Unmarking alignments
		for (i=pos; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (startA<substring.startA-THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			length=suffixOrigin-startA;
			j=Arrays.binarySearch(lengthBoundariesSuffix,0,lastLengthBoundarySuffix+1,length);
			if ( j>=0 || 
				 (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundarySuffix &&
				   ( Math.abs(length,lengthBoundariesSuffix[-j])<=THRESHOLD || 
				     (-j-1==0?false:Math.abs(length,lengthBoundariesSuffix[-j-3])<=THRESHOLD)
				   )
				 )
			   ) ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
		}
		for (i=pos+1; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (startA>=substring.endA) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			length=suffixOrigin-startA;
			j=Arrays.binarySearch(lengthBoundariesSuffix,0,lastLengthBoundarySuffix+1,length);
			if ( j>=0 || 
				 (-j-1)%3!=0 || 
				 ( -j-1<=lastLengthBoundarySuffix &&
				   ( Math.abs(length,lengthBoundariesSuffix[-j])<=THRESHOLD || 
				     (-j-1==0?false:Math.abs(length,lengthBoundariesSuffix[-j-3])<=THRESHOLD)
				   )
				 )
			   ) ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
		}
		
		// Updating properties of $substring$
		newSumEndA=0; newNEndA=0; newMinSuffixLength=substringLength;
		for (i=pos; i>=0; i--) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA<substring.startA-THRESHOLD) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			endA=Alignments.alignments[id][4];
			if (Math.abs(endA,substring.endA)>NEW_BOUNDARIES_THRESHOLD) continue;
			newSumEndA+=endA; newNEndA++;
			length=suffixOrigin-startA+1;
			if (length<newMinSuffixLength) newMinSuffixLength=length;
		}
		for (i=pos+1; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[i].id;
			startA=Alignments.alignments[id][3];
			if (startA>=substring.endA) break;
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=substring) continue;
			endA=Alignments.alignments[id][4];
			if (Math.abs(endA,substring.endA)>NEW_BOUNDARIES_THRESHOLD) continue;
			newSumEndA+=endA; newNEndA++;
			length=suffixOrigin-startA+1;
			if (length<newMinSuffixLength) newMinSuffixLength=length;
		}
		substring.sumEndA=newSumEndA;
		substring.nEndA=newNEndA;
		if (newNEndA>0) substring.endA=newSumEndA/newNEndA;
		substring.minSuffixLength=newMinSuffixLength;
		substring.hasLengthBoundariesSuffix=true;
		
		// Updating $lengthBiasIntervals$
		for (i=0; i<lastLengthBoundaryPrefix; i+=3) {
			lastLengthBiasInterval++;
			lengthBiasIntervals[lastLengthBiasInterval].firstPosition=suffixOrigin-lengthBoundariesSuffix[i];
			lengthBiasIntervals[lastLengthBiasInterval].center=suffixOrigin-lengthBoundariesSuffix[i+1];
			lengthBiasIntervals[lastLengthBiasInterval].lastPosition=suffixOrigin-lengthBoundariesSuffix[i+2];
			lengthBiasIntervals[lastLengthBiasInterval].readB=ReadA.id;
			lengthBiasIntervals[lastLengthBiasInterval].orientation=true;
		}
	}
	
	
	/**
	 * Keeps the intervals in $boundaries$ that contain an element of $candidates$,
	 * and the intervals whose center is close enough to an element of $candidates$.
	 *
	 * @return the last element of $boundaries$ after cleaning, or -1 if $boundaries$ gets
	 * empty.
	 */
	private static final int cleanBoundaries(int[] boundaries, int lastBoundary, int[] candidates, int lastCandidate, int threshold) {
		int i, j, k;
		
		Math.set(lengthBoundariesFlags,lastBoundary,false);
		for (i=0; i<=lastCandidate; i++) {
			j=Arrays.binarySearch(boundaries,0,lastBoundary+1,candidates[i]);
			if ( j>=0 || 
			     (-j-1)%3!=0 || 
				 ( -j-1<=lastBoundary &&
			       ( Math.abs(candidates[i],boundaries[-j])<=threshold || 
			         (-j-1==0?false:Math.abs(candidates[i],boundaries[-j-3])<=threshold)
				   )
				 )
		       ) {
				k=j>=0?j/3:(-j-1)/3;
				lengthBoundariesFlags[k]=true;
				lengthBoundariesFlags[k+1]=true;
				lengthBoundariesFlags[k+2]=true;
			}
		}
		j=-1;
		for (i=0; i<=lastBoundary; i++) {
			if (lengthBoundariesFlags[i]) {
				j++;
				if (j!=i) boundaries[j]=boundaries[i];
			}
		}
		return j;
	}
	
	
	/**
	 * Sorts $lengthBiasIntervals$ and removes contained intervals, if any.
	 */
	private static final void sortLengthBiasIntervals() {
		int i, j;
		Interval tmp;
		
		Interval.order=Interval.READB_ORIENTATION_FIRSTPOSITION;
		if (lastLengthBiasInterval>0) Arrays.sort(lengthBiasIntervals,0,lastLengthBiasInterval+1);
		for (i=0; i<=lastLengthBiasInterval; i++) {
			// Using $avgDiffs$ as a flag
			lengthBiasIntervals[i].avgDiffs=-1.0;
		}
		for (i=0; i<lastLengthBiasInterval; i++) {
			if (lengthBiasIntervals[i].avgDiffs==1.0) continue;
			for (j=i+1; j<=lastLengthBiasInterval; j++) {
				if (lengthBiasIntervals[j].firstPosition>=lengthBiasIntervals[i].lastPosition) break;
				if (lengthBiasIntervals[j].lastPosition<=lengthBiasIntervals[i].lastPosition) lengthBiasIntervals[j].avgDiffs=1.0;
			}
		}
		j=-1;
		for (i=0; i<=lastLengthBiasInterval; i++) {
			if (lengthBiasIntervals[i].avgDiffs==1.0) continue;
			j++;
			if (j!=i) {
				tmp=lengthBiasIntervals[j];
				lengthBiasIntervals[j]=lengthBiasIntervals[i];
				lengthBiasIntervals[i]=tmp;
			}
		}
		lastLengthBiasInterval=j;
	}
	

	/**
	 * Sets to null the $impliedBy*Substring$, $impliedBy*SubstringPrime$, and
	 * $inDenseSubstring$ pointers of an alignment, iff the readA interval of the 
	 * alignment straddles the dense substring pointed to, but does not have a large 
	 * enough intersection with it.
	 */
	private static final void discardAlignmentsAfterTrimming(int firstAlignment, int lastAlignment) {
		final double INTERSECTION_RATIO = 0.75;  // Arbitrary
		int i;
		int startA, endA, substringStartA, substringEndA, length;
		
		for (i=firstAlignment; i<=lastAlignment; i++) {
			startA=ReadA.sortedAlignments[i].startA();
			endA=ReadA.sortedAlignments[i].endA();
			length=endA-startA+1;
			if (ReadA.sortedAlignments[i].impliedByPrefixSubstring!=null) {
				substringStartA=ReadA.sortedAlignments[i].impliedByPrefixSubstring.startA;
				substringEndA=ReadA.sortedAlignments[i].impliedByPrefixSubstring.endA;
				if (Intervals.straddles(startA,endA,substringStartA,substringEndA) && Intervals.intersectionLength(startA,endA,substringStartA,substringEndA)<length*INTERSECTION_RATIO) ReadA.sortedAlignments[i].impliedByPrefixSubstring=null;				
			}
			if (ReadA.sortedAlignments[i].impliedBySuffixSubstring!=null) {
				substringStartA=ReadA.sortedAlignments[i].impliedBySuffixSubstring.startA;
				substringEndA=ReadA.sortedAlignments[i].impliedBySuffixSubstring.endA;
				if (Intervals.straddles(startA,endA,substringStartA,substringEndA) && Intervals.intersectionLength(startA,endA,substringStartA,substringEndA)<length*INTERSECTION_RATIO) ReadA.sortedAlignments[i].impliedBySuffixSubstring=null;				
			}
			if (ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime!=null) {
				substringStartA=ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime.startA;
				substringEndA=ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime.endA;
				if (Intervals.straddles(startA,endA,substringStartA,substringEndA) && Intervals.intersectionLength(startA,endA,substringStartA,substringEndA)<length*INTERSECTION_RATIO) ReadA.sortedAlignments[i].impliedByPrefixSubstringPrime=null;				
			}
			if (ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime!=null) {
				substringStartA=ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime.startA;
				substringEndA=ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime.endA;
				if (Intervals.straddles(startA,endA,substringStartA,substringEndA) && Intervals.intersectionLength(startA,endA,substringStartA,substringEndA)<length*INTERSECTION_RATIO) ReadA.sortedAlignments[i].impliedBySuffixSubstringPrime=null;
			}
			if (ReadA.sortedAlignments[i].inDenseSubstring!=null) {
				substringStartA=ReadA.sortedAlignments[i].inDenseSubstring.startA;
				substringEndA=ReadA.sortedAlignments[i].inDenseSubstring.endA;
				if (Intervals.straddles(startA,endA,substringStartA,substringEndA) && Intervals.intersectionLength(startA,endA,substringStartA,substringEndA)<length*INTERSECTION_RATIO) ReadA.sortedAlignments[i].inDenseSubstring=null;				
			}
		}
	}
	
	
	/**
	 * Ensures that dense substrings of substring type are no larger than their implied
	 * alignments. This might happen, since e.g. the leftmost alignments in a dense 
	 * substring might have been reassigned to an alignment interval.
	 *
	 * Remark: only B-maximal alignments are used for updating the boundaries; only
	 * alignment ends that are not too far from the original ends of the substring are 
	 * used.
	 *
	 * @param tmpArray temporary space, with at least two cells.
	 */
	public static final void fixBoundariesOfSubstringType(int[] tmpArray) {
		final int DISTANCE_THRESHOLD = IO.quantum<<2;  // Arbitrary
		final int WINDOW_LARGE = 3;  // Arbitrary
		final int WINDOW_SMALL = 1;
		final double MIN_RATIO = 0.8;  // Arbitrary
		final int MIN_ALIGNMENT_LENGTH = Alignments.minAlignmentLength-identityThreshold;
		int i;
		int id, start, end;
		DenseSubstring substring;
		AlignmentInterval interval;
		PeriodicSubstringInterval periodicInterval;
		
		for (i=0; i<=lastSubstring; i++) {
			substring=substrings[i];
			if (substring.substringReplication) {
				// Using temporary fields $leftAlignment$ and $rightAlignment$ to cumulate
				// min start position and max end position.
				substring.leftAlignment=Math.POSITIVE_INFINITY;
				substring.rightAlignment=-1;
				// Using $destinationAlignment$ and $firstAlignment$ to cumulate the 
				// number of assigned alignments inside the old (respectively, new)
				// boundaries.
				substring.destinationAlignment=0;
				substring.firstAlignment=0;
			}
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (substring!=null) continue;
			interval=ReadA.sortedAlignments[i].mergedToInterval;
			if (interval!=null) continue;
			periodicInterval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if (periodicInterval!=null) continue;
			substring=ReadA.sortedAlignments[i].inDenseSubstring;
			if (substring==null) continue;
			substring.destinationAlignment++;
			id=ReadA.sortedAlignments[i].id;
			start=Alignments.alignments[id][3];
			end=Alignments.alignments[id][4];
			tmpArray[0]=start; tmpArray[1]=end;
			Reads.trim(tmpArray,WINDOW_LARGE,true,true);
			if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<MIN_ALIGNMENT_LENGTH) continue;
			Reads.trim(tmpArray,WINDOW_SMALL,true,true);
			if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<MIN_ALIGNMENT_LENGTH) continue;
			start=tmpArray[0]; end=tmpArray[1];
			if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && start>=substring.startA-DISTANCE_THRESHOLD && start<substring.leftAlignment) substring.leftAlignment=start;
			if (ReadA.sortedAlignments[i].isRightMaximalB==1 && end<=substring.endA+DISTANCE_THRESHOLD && end>substring.rightAlignment) substring.rightAlignment=end;
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (substring!=null) continue;
			interval=ReadA.sortedAlignments[i].mergedToInterval;
			if (interval!=null) continue;
			periodicInterval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if (periodicInterval!=null) continue;
			substring=ReadA.sortedAlignments[i].inDenseSubstring;
			if (substring==null) continue;
			id=ReadA.sortedAlignments[i].id;
			start=Alignments.alignments[id][3];
			end=Alignments.alignments[id][4];
			tmpArray[0]=start; tmpArray[1]=end;
			Reads.trim(tmpArray,WINDOW_LARGE,true,true);
			if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<MIN_ALIGNMENT_LENGTH) continue;
			Reads.trim(tmpArray,WINDOW_SMALL,true,true);
			if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<MIN_ALIGNMENT_LENGTH) continue;
			start=tmpArray[0]; end=tmpArray[1];
			if (Intervals.isApproximatelyContained(start,end,substring.leftAlignment,substring.rightAlignment)) substring.firstAlignment++;
		}
		for (i=0; i<=lastSubstring; i++) {
			substring=substrings[i];
			if (substring.substringReplication && substring.firstAlignment>=substring.destinationAlignment*MIN_RATIO) {
				if (substring.leftAlignment!=Math.POSITIVE_INFINITY) substring.startA=substring.leftAlignment;
				if (substring.rightAlignment!=-1) substring.endA=substring.rightAlignment;
			}
		}
	}
	
	
	/**
	 * Remark: the procedure assumes $substrings$ to be sorted by $startA$.
	 * 
	 * @param tmpSubstring temporary space;
	 * @return TRUE iff $position$ belongs to a substring in $substrings$, and at distance 
	 * at least $distanceThreshold$ from each end.
	 */
	public static final boolean inDenseSubstring(int position, int distanceThreshold, DenseSubstring tmpSubstring) {
		int i;
		if (lastSubstring==-1) return false;
		
		tmpSubstring.startA=position;
		i=Arrays.binarySearch(substrings,0,lastSubstring+1,tmpSubstring);
		if (i<0) i=-i-1;
		i--;
		while (i>=0) {
			if (position>substrings[i].startA+distanceThreshold && position<substrings[i].endA-distanceThreshold) return true;
			i--;
		}
		return false;
	}
	
	
	/**
	 * Sets the $leftWindow,rightWindow,maxWindow$ fields of all dense substrings.
	 *
	 * Remark: the procedure sorts $ReadA.sortedAlignments$ by implied-by-dense and in-
	 * dense substring pointers.
	 */
	public static final void computeDensityWindows(int lastAlignment) {
		final int WINDOW_SIZE = IO.quantum<<1;  // Arbitrary
		int i, j;
		DenseSubstring tmpSubstring;
		
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].leftWindow=0;
			substrings[i].rightWindow=0;
			substrings[i].maxWindow=0;
			substrings[i].surface=0;
			// Assigning an ID to each dense substring is necessary for the following
			// sort.
			substrings[i].id=i;
		}
		Alignment.order=Alignment.IMPLIEDBYDENSE_INDENSE_STARTA;
		if (lastAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,lastAlignment+1);
		
		// Implied by dense substring
		i=0; tmpSubstring=ReadA.sortedAlignments[0].impliedByDenseSubstring;
		while (i<=lastAlignment && tmpSubstring!=null) {
			j=i+1;
			while (j<=lastAlignment && ReadA.sortedAlignments[j].impliedByDenseSubstring!=null && ReadA.sortedAlignments[j].impliedByDenseSubstring==tmpSubstring) j++;
			computeDensityWindows_impl(tmpSubstring,WINDOW_SIZE,i,j-1);
			i=j; tmpSubstring=i<=lastAlignment?ReadA.sortedAlignments[i].impliedByDenseSubstring:null;
		}
		
		// In dense substring of substring type
		tmpSubstring=i<=lastAlignment?ReadA.sortedAlignments[i].inDenseSubstring:null;
		while (i<=lastAlignment && tmpSubstring!=null) {
			j=i+1;
			while (j<=lastAlignment && ReadA.sortedAlignments[j].inDenseSubstring!=null && ReadA.sortedAlignments[j].inDenseSubstring==tmpSubstring) j++;
			computeDensityWindows_impl(tmpSubstring,WINDOW_SIZE,i,j-1);
			i=j; tmpSubstring=i<=lastAlignment?ReadA.sortedAlignments[i].inDenseSubstring:null;
		}		
	}


	/**
	 * Sets the $leftWindow,rightWindow,maxWindow$ fields of $susbstring$.
	 *
	 * @param $ReadA.sortedAlignments[firstAlignment..lastAlignment]$ is assumed to 
	 * contain all and only the alignments assigned to $substring$ (sorted in any order).
	 */
	private static final void computeDensityWindows_impl(DenseSubstring substring, int windowSize, int firstAlignment, int lastAlignment) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i, j;
		int lastPoint, sum, maxSum;
		Point[] tmpPoints = singleDeletions;  // Just used as temporary space
		
		// Collecting maximal events
		lastPoint=-1;
		for (i=firstAlignment; i<=lastAlignment; i++) {
			if (ReadA.sortedAlignments[i].isLeftMaximal==1) {
				lastPoint++;
				tmpPoints[lastPoint].position=ReadA.sortedAlignments[i].startA();
				tmpPoints[lastPoint].mass=1;
			}
			if (ReadA.sortedAlignments[i].isRightMaximal==1) {
				lastPoint++;
				tmpPoints[lastPoint].position=ReadA.sortedAlignments[i].endA();
				tmpPoints[lastPoint].mass=1;
			}
			substring.surface+=Intervals.intersectionLength(ReadA.sortedAlignments[i].startA(),ReadA.sortedAlignments[i].endA(),substring.startA,substring.endA);
		}
		lastPoint=Points.sortAndCompact(tmpPoints,lastPoint);
		if (lastPoint==-1) return;
		
		// Computing windows
		sum=0;
		for (i=0; i<=lastPoint; i++) {
			if (tmpPoints[i].position<substring.startA-IDENTITY_THRESHOLD) continue;
			if (tmpPoints[i].position>substring.startA+IDENTITY_THRESHOLD) break;
			sum+=tmpPoints[i].mass;
		}
		substring.leftWindow=sum;
		sum=0;
		for (i=lastPoint; i>=0; i--) {
			if (tmpPoints[i].position>substring.endA+IDENTITY_THRESHOLD) continue;
			if (tmpPoints[i].position<substring.endA-IDENTITY_THRESHOLD) break;
			sum+=tmpPoints[i].mass;
		}
		substring.rightWindow=sum;
		i=0;
		while (i<=lastPoint && tmpPoints[i].position<substring.startA) i++;
		j=i+1; sum=tmpPoints[i].mass; maxSum=0;
		while (i<=lastPoint && tmpPoints[i].position<=substring.endA-windowSize+1) {
			if (j>lastPoint || tmpPoints[j].position>=tmpPoints[i].position+windowSize) {
				maxSum=Math.max(maxSum,sum);
				i++; j=i+1; sum=i<=lastPoint?tmpPoints[i].mass:0;
				continue;
			}
			sum+=tmpPoints[j].mass;
			j++;
		}
		substring.maxWindow=maxSum;
	}
	
	
	/**
	 * Discards every pair (X,Y) of substrings such that: (1) X is strictly-prefix, Y is
	 * strictly-suffix; (2) X and Y overlap; (3) there is a long-period interval whose 
	 * first position is similar to the first position of X, and whose last position is 
	 * similar to the last position of Y. Such substrings are more likely to belong to 
	 * the periodic interval than to be independent repeats.
	 * A similar criterion is applied also to non-overlapping pairs of weak substrings.
	 *
	 * Remark: the procedure assumes $substrings$ to be sorted by startA and 
	 * $PeriodicSubstrings.longPeriodIntervals$ to be sorted by firstPosition.
	 */
	private static final void discardPeriodicSubstrings_prefixSuffix(int nAlignmentsHigh) {
		final double SURFACE_RATIO = 0.9;  // Arbitrary
		boolean reassigned;
		int i, j, k;
		int firstJForNextI, found, distance;
		DenseSubstring tmpSubstring;
		PeriodicSubstringInterval tmpInterval;
		
		// Discarding substrings: first criterion.
		for (i=0; i<=lastSubstring; i++) substrings[i].periodicSubstringInterval=null;
		for (i=0; i<nAlignmentsHigh; i++) {
			tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (tmpSubstring!=null) tmpSubstring.periodicSubstringInterval=null;
		}
		reassigned=false;
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastSubstring) {
			if ( !substrings[i].prefixReplication || substrings[i].suffixReplication || substrings[i].substringReplication || substrings[i].singleDeletionReplication ||
			     j>PeriodicSubstrings.lastLongPeriodInterval || PeriodicSubstrings.longPeriodIntervals[j].firstPosition>substrings[i].startA+identityThreshold
			   ) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (PeriodicSubstrings.longPeriodIntervals[j].firstPosition<substrings[i].startA-identityThreshold) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastSubstring && PeriodicSubstrings.longPeriodIntervals[j].firstPosition>=substrings[i+1].startA-identityThreshold) firstJForNextI=j;
			if (substrings[i].endA<=PeriodicSubstrings.longPeriodIntervals[j].lastPosition+identityThreshold) {
				distance=(int)((1.0-SURFACE_RATIO)*PeriodicSubstrings.longPeriodIntervals[j].length());
				found=-1;
				for (k=i+1; k<=lastSubstring; k++) {
					if (substrings[k].startA>substrings[i].endA+distance) break;
					if ( !substrings[k].suffixReplication || substrings[k].prefixReplication || substrings[k].substringReplication || substrings[k].singleDeletionReplication ||
					     Math.abs(substrings[k].endA,PeriodicSubstrings.longPeriodIntervals[j].lastPosition)>identityThreshold
					   ) continue;
					found=k;
					break;
				}
				if (found!=-1) {
					substrings[i].periodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[j];
					substrings[found].periodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[j];
					reassigned=true;
				}
			}
			j++;
		}
		
		// Discarding substrings: second criterion.
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastSubstring) {
			if ( substrings[i].periodicSubstringInterval!=null || !substrings[i].isWeak ||
				 (!substrings[i].prefixReplication && !substrings[i].singleDeletionReplication && !substrings[i].substringReplication) ||
			     j>PeriodicSubstrings.lastLongPeriodInterval || PeriodicSubstrings.longPeriodIntervals[j].firstPosition>substrings[i].startA+identityThreshold
			   ) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (substrings[j].periodicSubstringInterval!=null || PeriodicSubstrings.longPeriodIntervals[j].firstPosition<substrings[i].startA-identityThreshold) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastSubstring && PeriodicSubstrings.longPeriodIntervals[j].firstPosition>=substrings[i+1].startA-identityThreshold) firstJForNextI=j;
			if (substrings[i].endA<=PeriodicSubstrings.longPeriodIntervals[j].lastPosition+identityThreshold) {
				distance=(int)((1.0-SURFACE_RATIO)*PeriodicSubstrings.longPeriodIntervals[j].length());
				found=-1;
				for (k=i+1; k<=lastSubstring; k++) {
					if (substrings[k].startA>substrings[i].endA+distance) break;
					if ( !substrings[k].isWeak || 
						 (!substrings[k].suffixReplication && !substrings[k].substringReplication && !substrings[k].singleDeletionReplication) ||
					     Math.abs(substrings[k].endA,PeriodicSubstrings.longPeriodIntervals[j].lastPosition)>identityThreshold
					   ) continue;
					found=k;
					break;
				}
				if (found!=-1) {
					substrings[i].periodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[j];
					substrings[found].periodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[j];
					reassigned=true;
				}
			}
			j++;
		}
		if (!reassigned) return;
		
		// Compacting substrings
		j=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].periodicSubstringInterval!=null) continue;
			j++;
			if (j!=i) {
				tmpSubstring=substrings[j];
				substrings[j]=substrings[i];
				substrings[i]=tmpSubstring;
			}
		}
		lastSubstring=j;
		
		// Updating pointers from alignments
		for (i=0; i<nAlignmentsHigh; i++) {
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) continue;
			tmpInterval=ReadA.sortedAlignments[i].impliedByDenseSubstring.periodicSubstringInterval;
			if (tmpInterval!=null) {
				ReadA.sortedAlignments[i].periodicSubstringInterval=tmpInterval;
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
			}
		}
	}
	
	
	/**
	 * Discards every pair (X,Y) of substrings such that: (1) X is strictly-prefix, Y is
	 * strictly-suffix; (2) X and Y are separated by a small space; (3) there is a weak 
	 * dense susbtring of substring type whose first position is similar to the first 
	 * position of X, and whose last position is similar to the last position of Y. 
	 * This procedure is essentially the same as $discardPeriodicSubstrings_
	 * prefixSuffix()$.
	 *
	 * Remark: the procedure assumes $substrings$ to be sorted by startA.
	 */
	public static final void discardPrefixSuffix_weakSubstring(int nAlignmentsHigh) {
		final double SURFACE_RATIO = 0.9;  // Arbitrary
		boolean found;
		int i, j, k;
		int firstJForNextI, startA, endA, length, distance;
		DenseSubstring substring;
		
		// Discarding substrings
		for (i=0; i<=lastSubstring; i++) substrings[i].representative=null;
		for (i=0; i<nAlignmentsHigh; i++) {
			substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (substring!=null) substring.representative=null;
		}
		found=false;
		for (i=0; i<=lastSubstring; i++) {
			if (!substrings[i].substringReplication || !substrings[i].isWeak) continue;
		 	startA=substrings[i].startA; endA=substrings[i].endA;
			length=endA-startA+1; distance=(int)((1.0-SURFACE_RATIO)*length);
			for (j=i-1; j>=0; j--) {
				if (substrings[j].startA<startA-identityThreshold) break;
				if ( substrings[j].endA>endA+identityThreshold ||
					 !substrings[j].prefixReplication || substrings[j].suffixReplication || substrings[j].substringReplication || substrings[j].singleDeletionReplication
				   ) continue;
				for (k=j+1; k<=lastSubstring; k++) {
					if (substrings[k].startA>=endA-identityThreshold) break;
					if ( k==i ||
						 substrings[k].prefixReplication || !substrings[k].suffixReplication || substrings[k].substringReplication || substrings[k].singleDeletionReplication ||
						 Math.abs(substrings[k].endA,endA)>identityThreshold
					   ) continue;
					if (!(substrings[j].isWeak && substrings[k].isWeak) && substrings[k].startA>substrings[j].endA+distance) continue;
					if (substrings[j].representative==null || length<substrings[j].representative.length()) substrings[j].representative=substrings[i];
					if (substrings[k].representative==null || length<substrings[k].representative.length()) substrings[k].representative=substrings[i];
					found=true;
				}
			}
			for (j=i+1; j<=lastSubstring; j++) {
				if (substrings[j].startA>startA+identityThreshold) break;
				if ( substrings[j].endA>endA+identityThreshold ||
					 !substrings[j].prefixReplication || substrings[j].suffixReplication || substrings[j].substringReplication || substrings[j].singleDeletionReplication
				   ) continue;
				for (k=j+1; k<=lastSubstring; k++) {
					if (substrings[k].startA>=endA-identityThreshold) break;
					if ( substrings[k].prefixReplication || !substrings[k].suffixReplication || substrings[k].substringReplication || substrings[k].singleDeletionReplication ||
						 Math.abs(substrings[k].endA,endA)>identityThreshold
					   ) continue;
					if (!(substrings[j].isWeak && substrings[k].isWeak) && substrings[k].startA>substrings[j].endA+distance) continue;
					if (substrings[j].representative==null || length<substrings[j].representative.length()) substrings[j].representative=substrings[i];
					if (substrings[k].representative==null || length<substrings[k].representative.length()) substrings[k].representative=substrings[i];
					found=true;
				}
			}
		}
		if (!found) return;
		
		// Compacting substrings
		j=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].representative!=null) continue;
			j++;
			if (j!=i) {
				substring=substrings[j];
				substrings[j]=substrings[i];
				substrings[i]=substring;
			}
		}
		lastSubstring=j;
		
		// Updating pointers from alignments
		for (i=0; i<nAlignmentsHigh; i++) {
			substring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
			if (substring==null || substring.representative==null) continue;
			ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
			ReadA.sortedAlignments[i].inDenseSubstring=substring.representative;
		}
	}
	
	
	/**
	 * Copies the content of $ReadA.lastDeltaPoint$ into one of 
	 * $deltaPoints_backup{Left,Right}$.
	 */
	private static final void cloneDeltaPoints(boolean leftOrRight) {
		int i;
		
		if (leftOrRight) {
			lastDeltaPoint_backupLeft=ReadA.lastDeltaPoint;
			for (i=0; i<=lastDeltaPoint_backupLeft; i++) deltaPoints_backupLeft[i].copyFrom(ReadA.deltaPoints[i]);
		}
		else {
			lastDeltaPoint_backupRight=ReadA.lastDeltaPoint;
			for (i=0; i<=lastDeltaPoint_backupRight; i++) deltaPoints_backupRight[i].copyFrom(ReadA.deltaPoints[i]);
		}
	}
	
	
	/**
	 * Stores in the temporary variable $previousSumStartA$ of every dense substring, the 
	 * number of its implied alignments that are approximately identical to the entire 
	 * substring.
	 */
	private static final void countAlignmentsAtBoundaries(int lastAlignment) {
		int i;
		Alignment alignment;
		DenseSubstring substring;
		
		for (i=0; i<=lastSubstring; i++) substrings[i].previousSumStartA=0;
		for (i=0; i<=lastAlignment; i++) {
			alignment=ReadA.sortedAlignments[i];
			substring=alignment.impliedByDenseSubstring;
			if (substring!=null && Intervals.areApproximatelyIdentical(alignment.startA(),alignment.endA(),substring.startA,substring.endA)) substring.previousSumStartA++;
			substring=alignment.inDenseSubstring;
			if (substring!=null && Intervals.areApproximatelyIdentical(alignment.startA(),alignment.endA(),substring.startA,substring.endA)) substring.previousSumStartA++;
		}
	}
	
	
	/**
	 * Remark: the procedure assumes $substrings$ to be sorted by $startA$.
	 *
	 * @return TRUE iff $[from..to]$ is the prefix of a prefix substring, or a suffix of a
	 * suffix substring.
	 */
	public static final boolean isPrefixOfPrefixSubstring(int from, int to) {
		int i, j;
		
		tmpDenseSubstring.startA=from;
		j=Arrays.binarySearch(substrings,0,lastSubstring+1,tmpDenseSubstring);
		if (j<0) j=-j-1;
		for (i=j; i<=lastSubstring; i++) {
			if (substrings[i].startA>from+identityThreshold) break;
			if ( Intervals.isApproximatelyContained(from,to,substrings[i].startA,substrings[i].endA) && 
				 (substrings[i].prefixReplication || substrings[i].singleDeletionReplication)
			   ) return true;
		}
		for (i=j-1; i>=0; i--) {
			if (substrings[i].endA<to-identityThreshold) continue;
			if ( Intervals.isApproximatelyContained(from,to,substrings[i].startA,substrings[i].endA) && 
				 ( ((substrings[i].prefixReplication || substrings[i].singleDeletionReplication) && Math.abs(from,substrings[i].startA)<=identityThreshold) || 
				   ((substrings[i].suffixReplication || substrings[i].singleDeletionReplication) && Math.abs(to,substrings[i].endA)<=identityThreshold)	 
				 )
			   ) return true;
		}
		return false;
	}
	
	
	/**
	 * Remark: the procedure assumes $substrings$ to be sorted by $startA$.
	 *
	 * @return TRUE iff $pos$ is strictly inside a dense substring of substring type.
	 */
	public static final boolean inSubstringType(int pos) {
		int i, j;
		
		tmpDenseSubstring.startA=pos;
		j=Arrays.binarySearch(substrings,0,lastSubstring+1,tmpDenseSubstring);
		if (j<0) j=-j-1;
		for (i=j-1; i>=0; i--) {
			if (substrings[i].substringReplication && substrings[i].startA<pos-identityThreshold && substrings[i].endA>pos+identityThreshold) return true;
		}
		return false;
	}
	

}