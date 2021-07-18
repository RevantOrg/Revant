package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import de.mpi_cbg.revant.util.Biology;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Leaf;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.DAG;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Leaves;
import de.mpi_cbg.revant.util.Histograms;


public class PeriodicSubstrings {
	/**
	 * Parameters of the pipeline
	 */
	public static int MAX_NEIGHBORS, MAX_SUBSTRINGS;
	public static final int MAX_DISTINCT_SHIFTS = 500;  // Maximum number of distinct observed shifts
	public static final int MIN_DISTINCT_SHIFTS = 10;
	public static final int MIN_PERIOD = IO.quantum;  // Length of the shortest period to be considered. The procedure can detect arbitrarily small exact periods by setting down this parameter. However, setting it too small can lead to underestimating the period.
	public static final int MIN_PERIOD_LONG = Alignments.minAlignmentLength;  // For long-period substrings. Arbitrary.
	public static final int IDENTITY_THRESHOLD_PERIODIC = MIN_PERIOD>>1;  // Two positions are considered identical if they differ by at most this threshold
	public static final int IDENTITY_THRESHOLD_PERIODIC_LONG = IO.quantum;
	public static double maxLengthRatio = 1.5;  // An alignment $w$ of type 1 or 2 is connected to a predecessor alignment $v$ if the ratio between the length of $w$ and the length of the shortest (respectively, longest) alignment in every path to $v$ (with a specific state) is at most $maxLengthRatio$.
	public static int minIntervalLength = MIN_PERIOD>>2;  // Intervals of a density estimation tree are not split if they are this length or smaller
	public static int minLocalMaxDistance = MIN_PERIOD;  // Minimum distance between two local maxima of a density estimation tree for them to be considered distinct clusters
	private static final int DEFAULT_MIN_PATH_LENGTH = 4;
	private static final int MAX_PATHLENGTH_DIFFERENCE = 4;  // If there are "few" path length values, and if such values differ by at most this quantity, we assume they are uniformly distributed and we don't try to split them.

	/**
	 * Parameters of the pipeline, estimated from the data.
	 */
	public static int minPathLength;  // Minimum length of a path to represent a periodic substring

	/**
	 * Output data structures
	 */
	public static PeriodicSubstring[] substrings;
	public static PeriodicSubstringInterval[] intervals;
	public static int lastSubstring, lastSubstringForFactoring, lastInterval;
	public static int lastSubstring_artificial, lastLongPeriodSubstring_artificial;
	private static Leaf[] tmpLeaves;
	public static PeriodicSubstring[] longPeriodSubstrings;
	public static int lastLongPeriodSubstring;
	public static PeriodicSubstringInterval[] longPeriodIntervals;
	public static int lastLongPeriodInterval;

	/**
	 * Number of states in the periodic substring automaton
	 */
	private static final int N_STATES = 9;

	/**
	 * DAG data structures
	 */
	private static int[][] inNeighbors, inNeighborStates, outNeighbors;  // The DAG of alignments
	private static int[] nInNeighbors, nInNeighborsPrime, nOutNeighbors;
	private static int[] minimalVertices;  // The minimal vertices of the DAG
	private static int[] components;  // Connected component of each vertex of the DAG
	private static int[] nVerticesPerComponent;  // Number of minimal vertices per connected component of the DAG
	private static int[] sorted2original;  // Maps a rank in the topological order of the DAG to the ID of the corresponding vertex
	private static int[] original2sorted;  // Reverses $sorted2original$
	private static int[][] distances;  // The length of a longest path from a source alignment to all other alignments, for every possible state.
	private static int[][] predecessors;  // The predecessor of each vertex (columns), for every possible state (rows), in the longest path from a source alignment.
	private static int[][] predecessorStates;  // The state of the predecessor of every vertex (columns), for every possible state (rows), in the longest path from a source alignment.
	private static int[][] minIntervalLengthA, maxIntervalLengthA;  // Min/max length of an interval in $readA$ over all paths from the source alignment to every alignment, for every possible state.
	private static int[][] minIntervalLengthB, maxIntervalLengthB;  // Min/max length of an interval in $readB$ over all paths from the source alignment to every alignment, for every possible state.
	private static int[] tmp;

	/**
	 * Data structures for estimating periods
	 */
	private static int[] alignments;  // Stores the position of every alignment in the longest path returned by procedure $getCandidateSubstring$, relative to $firstAlignment$.
	private static int[] alignmentStates;  // The state of every alignment in $alignments$
	private static int lastAlignment;  // Last element in $alignments$ and $alignmentStates$
	private static Point[] shifts, shiftsPrime;
	private static int lastShift;

	/**
	 * Data structures for computing splits
	 */
	public static int[] leftSplits, rightSplits;
	public static boolean[] leftSplitFlags;  // 1=beginning of a range.
	public static boolean[] rightSplitFlags;  // 1=end of a range.
	public static boolean[] leftSplitLowQuality, rightSplitLowQuality;

	/**
	 * Data structures for connected components
	 */
	public static int[] stack, stackPrime;
	private static Point[] lengthPoints, lengthPointsPrime;
	private static int lastLengthPointPrime;
	private static int[] distancesPrime;
	
	/**
	 * Temporary space
	 */
	public static int[] periodTmp;  // See procedure $estimatePeriod$ for details
	private static int previousLongPeriodOrder;
	private static boolean splitsUpdated;  // Used by $splitShortPeriodIntervals()$.
	public static int lastLeftSplit, lastRightSplit;  // Used by $splitShortPeriodIntervals()$.
	private static int lastLeftSplit_long, lastRightSplit_long;  // Used by $splitLongPeriodIntervals()$.
	private static PeriodicSubstringInterval tmpPeriodicInterval;


	public static final void allocateMemory(int maxAlignments, int maxSubstrings, int maxSplits, int maxNeighbors) {
		int i;
		
		MAX_NEIGHBORS=maxNeighbors;
		MAX_SUBSTRINGS=maxSubstrings;

		// Output data structures
		substrings = new PeriodicSubstring[MAX_SUBSTRINGS];
		for (i=0; i<substrings.length; i++) {
			substrings[i] = new PeriodicSubstring();
			substrings[i].allocateMemory(MIN_DISTINCT_SHIFTS);
			substrings[i].hasLongPeriod=false;
		}
		intervals = new PeriodicSubstringInterval[MAX_SUBSTRINGS];
		for (i=0; i<maxSubstrings; i++) {
			intervals[i] = new PeriodicSubstringInterval();
			intervals[i].allocateMemory(MIN_DISTINCT_SHIFTS);
			intervals[i].hasLongPeriod=false;
		}
		longPeriodSubstrings = new PeriodicSubstring[MAX_SUBSTRINGS];
		for (i=0; i<substrings.length; i++) {
			longPeriodSubstrings[i] = new PeriodicSubstring();
			longPeriodSubstrings[i].allocateMemory(MIN_DISTINCT_SHIFTS);
			longPeriodSubstrings[i].hasLongPeriod=true;
		}
		longPeriodIntervals = new PeriodicSubstringInterval[MAX_SUBSTRINGS];
		for (i=0; i<maxSubstrings; i++) {
			longPeriodIntervals[i] = new PeriodicSubstringInterval();
			longPeriodIntervals[i].allocateMemory(MIN_DISTINCT_SHIFTS);
			longPeriodIntervals[i].hasLongPeriod=true;
		}

		// DAG data structures
		inNeighbors = new int[maxAlignments][MAX_NEIGHBORS];
		inNeighborStates = new int[maxAlignments][MAX_NEIGHBORS];
		nInNeighbors = new int[maxAlignments];
		nInNeighborsPrime = new int[maxAlignments];
		outNeighbors = new int[maxAlignments][MAX_NEIGHBORS];
		nOutNeighbors = new int[maxAlignments];
		minimalVertices = new int[maxAlignments];
		components = new int[maxAlignments];
		nVerticesPerComponent = new int[maxAlignments];
		sorted2original = new int[maxAlignments];
		original2sorted = new int[maxAlignments];
		distances = new int[maxAlignments][N_STATES];
		predecessors = new int[maxAlignments][N_STATES];
		predecessorStates = new int[maxAlignments][N_STATES];
		minIntervalLengthA = new int[maxAlignments][N_STATES];
		maxIntervalLengthA = new int[maxAlignments][N_STATES];
		minIntervalLengthB = new int[maxAlignments][N_STATES];
		maxIntervalLengthB = new int[maxAlignments][N_STATES];
		tmp = new int[2];

		// Data structures for estimating periods
		alignments = new int[maxAlignments];
		alignmentStates = new int[maxAlignments];
		shifts = new Point[maxAlignments<<2];
		for (i=0; i<shifts.length; i++) shifts[i] = new Point();
		shiftsPrime = new Point[maxAlignments<<2];
		for (i=0; i<shiftsPrime.length; i++) shiftsPrime[i] = new Point();

		// Data structures for computing splits
		leftSplits = new int[maxSplits*3];
		rightSplits = new int[maxSplits*3];
		leftSplitFlags = new boolean[maxSplits*3];
		rightSplitFlags = new boolean[maxSplits*3];
		leftSplitLowQuality = new boolean[maxSplits*3];
		rightSplitLowQuality = new boolean[maxSplits*3];

		// Data structures for connected components
		stack = new int[maxSubstrings];
		stackPrime = new int[maxAlignments];
		lengthPoints = new Point[maxAlignments];
		for (i=0; i<lengthPoints.length; i++) lengthPoints[i] = new Point();
		lengthPointsPrime = new Point[maxAlignments];
		for (i=0; i<lengthPointsPrime.length; i++) lengthPointsPrime[i] = new Point();
		distancesPrime = new int[maxAlignments];
		
		// Temporary space
		periodTmp = new int[5];
		tmpLeaves = new Leaf[DensityEstimationTree.leaves.length];
		for (i=0; i<tmpLeaves.length; i++) tmpLeaves[i] = new Leaf();
		tmpPeriodicInterval = new PeriodicSubstringInterval();
	}


	private static final void clearNeighborMatrices(int nAlignments) {
		int i, j;

		for (i=0; i<nAlignments; i++) {
			Math.set(inNeighbors[i],inNeighbors[i].length-1,-1);
		}
		for (i=0; i<nAlignments; i++) {
			Math.set(inNeighborStates[i],inNeighborStates[i].length-1,-1);
		}
		Math.set(nInNeighbors,nAlignments-1,0);
		for (i=0; i<nAlignments; i++) {
			Math.set(outNeighbors[i],outNeighbors[i].length-1,-1);
		}
		Math.set(nOutNeighbors,nAlignments-1,0);
	}


	/**
	 * Remark: it could happen that a substring of a repeated periodic substring occurs by
	 * itself multiple times in the genome (for example, a nonperiodic degenerate version
	 * of the period, that aligns to just one instance of the period of the repeated
	 * periodic substring). The boundaries of all such repeats, in all periodic substrings,
	 * will be detected in the following steps of the pipeline, by estimating the density
	 * of maximal events that come from alignments that are not implied by dense or
	 * periodic substrings.
	 *
	 * Remark: a periodic substring can contain a dense substring that replicates by 
	 * taking substrings of itself: this is already split using function $\delta$ by 
	 * procedure $DenseSubstrings.splitSubstrings$.
	 *
	 * Remark: a substring of a periodic substring could occur by itself multiple times in
	 * the genome, and still behave as a (different) periodic substring. The boundaries of
	 * such nested periodic substrings are detected by this procedure.
	 *
	 * For all the reasons above, there is no need to split a periodic substring using 
	 * function $\delta$.
	 *
	 * Remark: at the end of the procedure, some periodic substrings might not point to 
	 * any interval.
	 */
	public static final void detect() {
		boolean newShortPeriodInterval;
		int i, j;
		int period, nElements;
		double sum;
		
		// Initializing reused space
		lastSubstring=-1; lastInterval=-1;
		lastLeftSplit=-1; lastRightSplit=-1;
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			ReadA.sortedAlignments[i].impliedByPeriodicSubstring=null;
			ReadA.sortedAlignments[i].isPrefixOfPeriodicSubstring=false;
			ReadA.sortedAlignments[i].isSuffixOfPeriodicSubstring=false;
			ReadA.sortedAlignments[i].inPeriodicSubstring=false;
			ReadA.sortedAlignments[i].periodicSubstringInterval=null;
		}
		
		// Detecting substrings and intervals
		getSubstrings();
		if (lastSubstring!=-1) {
			for (i=0; i<=lastSubstring; i++) substrings[i].mergedToInterval=null;
			correctCoverage(substrings,lastSubstring);
			removeSubstringsForFactoring();
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println();
	System.err.println("SHORT-PERIOD PERIODIC SUBSTRINGS after removeSubstringsForFactoring():   only up to lastSubstringForFactoring="+lastSubstringForFactoring);
	for (int x=0; x<=lastSubstringForFactoring; x++) System.err.println(substrings[x].hashCode()+" :: "+substrings[x]);
	System.err.println();
}
		
			// Assigning an estimate of the period to each periodic substring is useful
			// when the estimate of the period of a periodic interval fails (see below).
			// Such estimate is also used by $atPeriodicSubstringEnd()$.
			for (i=0; i<=lastSubstring; i++) substrings[i].period=estimatePeriod(substrings[i].shifts,substrings[i].lastShift,MIN_PERIOD,minIntervalLength,minLocalMaxDistance);
			shortPeriod2longPeriod();
			initArtificialSubstrings();			
			getIntervals();
			markImpliedAlignments_shortPeriodIntervals(true);
			setIntervalMaximality(intervals,lastInterval);
			estimatePeriodOfShortPeriodIntervals();
			PeriodicSubstring.order=PeriodicSubstring.MINSTARTA_MAXENDA;
			if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
			lastSubstringForFactoring=-1;  // Invalidated by sorting, and never to be used after this point.
			markPrefixSuffix(substrings,lastSubstring,intervals,lastInterval);			
			setIntervalPointers();
		}
		else lastLongPeriodSubstring=-1;
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println();
	System.err.println("SHORT-PERIOD PERIODIC SUBSTRINGS DETECTED:");
	for (int x=0; x<=lastSubstring; x++) System.err.println(substrings[x].hashCode()+" :: "+substrings[x]);
	System.err.println();
	System.err.println("alignments:");
	for (int x=0; x<ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
	System.err.println();
	System.err.println("short-period intervals:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x]);
}
if (IO.CONSISTENCY_CHECKS) {
	for (int x=0; x<=lastSubstring; x++) {
		for (int y=1; y<=substrings[x].lastShift; y++) {
			if (substrings[x].shifts[y].position<substrings[x].shifts[y-1].position) {
				System.err.println("ERROR 1: shifts array not sorted.");
				System.exit(1);
			}
		}
	}
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) {
		if (ReadA.sortedAlignments[x].impliedByPeriodicSubstring!=null && !ReadA.sortedAlignments[x].impliedByPeriodicSubstring.hasLongPeriod && lastSubstring==-1) {
			System.err.println("The following alignment is assigned to a short-period periodic substring, but no such substring exists:");
			System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
			System.exit(1);
		}
		if (ReadA.sortedAlignments[x].impliedByPeriodicSubstring!=null && !ReadA.sortedAlignments[x].impliedByPeriodicSubstring.hasLongPeriod && ReadA.sortedAlignments[x].periodicSubstringInterval==null && ReadA.sortedAlignments[x].impliedByPeriodicSubstring.mergedToInterval!=null) {
			System.err.println("The following alignment is assigned to a short-period periodic substring but to no interval:");
			System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
			System.exit(1);
		}
	}
}

		
		// Detecting substrings with long periods
		splitsUpdated=false;
		newShortPeriodInterval=getSubstrings_longPeriods();
		if (newShortPeriodInterval) {
			if (!splitsUpdated) splitShortPeriodIntervals();
			if (lastInterval>0) {
				lastInterval=filterIntervalsWithPeaks(intervals,lastInterval,substrings,lastSubstring,false,true,true);
				lastInterval=filterIntervalsWithPeaks(intervals,lastInterval,substrings,lastSubstring,false,true,true);
			}
			// $filterIntervalsWithPeaks()$ and $filterIntervalsWithPeaks_smallSurface()$
			// might create new short-period intervals that are similar to existing long-
			// period intervals: such long-period intervals should be removed.
			discardLongPeriodIntervals_simple();
			// Updating short periods
			for (i=0; i<=lastInterval; i++) {
				period=estimatePeriod(intervals[i].shifts,intervals[i].lastShift,MIN_PERIOD,minIntervalLength,minLocalMaxDistance);
				if (period>0) {
					intervals[i].period=period;
					continue;
				}
				if (intervals[i].period>0) continue;
				sum=0.0; nElements=0;
				for (j=0; j<=lastSubstring; j++) {
					if (substrings[j].mergedToInterval==intervals[i] && substrings[j].period>0) {
						sum+=substrings[j].period;
						nElements++;
					}
				}
				period=sum==0.0?0:(int)(sum/nElements);
				intervals[i].period=period;
			}
			markImpliedAlignments_shortPeriodIntervals(true);
		}
		else {
			for (i=0; i<=lastInterval; i++) intervals[i].isContainedInShort=false;
		}
		setIntervalPointers();

if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("BEFORE mergeIntervals:");
	System.err.println("short-period intervals:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("long-period intervals:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
}


		// Copying long-period intervals into the set of all periodic intervals.
		//
		// Remark: this must be undone before factoring the next read.
		//
		mergeIntervals();
		markImpliedAlignments_peaks();
		addEvents(false);
	

if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("AFTER mergeIntervals:");
	System.err.println("short-period intervals:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("long-period intervals:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
Alignment.order=Alignment.STARTA_ENDA;
Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
	System.err.println((lastLongPeriodInterval+1)+" long-period periodic intervals detected");
}
		
		
		
if (IO.CONSISTENCY_CHECKS) {
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) {
		PeriodicSubstringInterval pInterval = ReadA.sortedAlignments[x].periodicSubstringInterval;
		if (pInterval==null) continue;
		if (pInterval.lastPosition>=Reads.getReadLength(ReadA.id)) {
			System.err.println("ERROR: the following alignment is assigned to a periodic substring that exceeds the read length: "+ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
			System.exit(1);
		}
	}		
}	
		
		
	}
	
	
	/**
	 * Remark: the procedure assumes $substrings$, $longPeriodSubstrings$ and $intervals$
	 * to be sorted by first position.
	 */
	private static final void estimatePeriodOfShortPeriodIntervals() {
		int i, j;
		int firstJForNextI, period;
		
		// Estimating periods using shifts
		for (i=0; i<=lastInterval; i++) {
			period=estimatePeriod(intervals[i].shifts,intervals[i].lastShift,MIN_PERIOD,minIntervalLength,minLocalMaxDistance);
			if (period>0) intervals[i].period=period;
			else {
				intervals[i].period=-1;
				intervals[i].tmpInt1=0;  // Used as temporary space
				intervals[i].tmpInt2=0;  // Used as temporary space
			}
		}
		
		// Estimating period using substring periods
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastInterval) {
			if (intervals[i].period>=0 || j>lastSubstring || substrings[j].minStartA>=intervals[i].lastPosition) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (substrings[j].maxEndA<=intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastInterval && substrings[j].maxEndA>=intervals[i+1].firstPosition) firstJForNextI=j;
			if (substrings[j].period>0 && substrings[j].mergedToInterval==intervals[i]) {
				intervals[i].tmpInt1+=substrings[j].period;
				intervals[i].tmpInt2++;
			}
			j++;
		}
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastInterval) {
			if (intervals[i].period>=0 || j>lastLongPeriodSubstring || longPeriodSubstrings[j].minStartA>=intervals[i].lastPosition) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (longPeriodSubstrings[j].maxEndA<=intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastInterval && longPeriodSubstrings[j].maxEndA>=intervals[i+1].firstPosition) firstJForNextI=j;
			if (longPeriodSubstrings[j].period>0 && longPeriodSubstrings[j].mergedToInterval==intervals[i]) {
				intervals[i].tmpInt1+=longPeriodSubstrings[j].period;
				intervals[i].tmpInt2++;
			}
			j++;
		}
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].period>=0) continue;
			intervals[i].period=intervals[i].tmpInt1==0.0?0:(intervals[i].tmpInt1/intervals[i].tmpInt2);
		}
	}


	/**
	 * Tries to detect all the occurrences of repeated periodic substrings of the genome
	 * in a read. A repeated periodic substring of the genome is a maximal substring $W$
	 * with period $P$ smaller than the minimum alignment length, and such that either:
	 * (1) the genome contains at least one other maximal substring with the same period
	 * $P$ (but possibly with different length), or (2) $|W|$ is greater than the read
	 * length (e.g. in telomeres).
	 *
	 * Assume that $R[i..j]$ is an occurrence of a repeated periodic substring that also
	 * occurs in read $R'[i'..j']$, and assume that the alignment algorithm finds maximal
	 * alignments. Then, we expect to see a sequence of overlapping maximal alignments
	 * between $R$ and $R'$ that span prefixes of $[i..j]$ of increasing length in $R$,
	 * and suffixes of $[i'..j']$ of increasing length in $R'$. We also expect to see a
	 * sequence of maximal overlapping alignments between $R$ and $R'$ that span suffixes
	 * of $[i..j]$ of decreasing length in $R$, and prefixes of $[i'..j']$ of decreasing
	 * length in $R'$. Between these two states (denoted throughout by state 0 and state
	 * 3, respectively), we expect to see a sequence of overlapping maximal alignments of
	 * size approximately $j'-i'+1$ if $j'-i'$ is smaller than $j-i$ (state 1), or a 
	 * sequence of overlapping alignments of size $j-i+1$ if $j'-i'$ is bigger than $j-i$ 
	 * (state 2).
	 *
	 * This procedure finds maximal sequences of overlapping alignments that traverse such
	 * states (see procedure $getSubstring$). We use also alignments that are not maximal,
	 * since they could trace a repeated periodic substring as well (e.g. if both the A
	 * and the B read coincide with repeated periodic substrings).
	 *
	 * The procedure does not return repeated periodic substrings directly, but a set of
	 * \emph{candidate substrings}, to be refined downstream, each induced by a given B
	 * read in a given orientation, each with its own set of possible starting and ending
	 * positions in the A read and in the B read (i.e. not necessarily with just one
	 * starting and one ending position), and each with its own estimate of the period.
	 *
	 * Array $substrings$ stores a description of all such candidate substrings.
     * The set of starting and ending positions of all candidate substrings is stored in
     * array $Events.events$.
     *
     * Remark: the procedure handles the case in which $R[i..j]$ is not in the same phase
     * as $R'[i'..j']$. This creates additional states.
     *
     * Remark: the intervals in readA of distinct candidate periodic substrings can be
     * nested in each other, or straddle one another. This captures the case, e.g., of a
     * periodic substring in which an instance of the period has mutated, while still
     * being periodic, and occurs by itself in another read without the containing
     * periodic substring.
     *
     * Remark: each DAG built by this procedure typically has a number of connected
	 * components, some of which are induced by noise and expected to be small. The
	 * procedure estimates the distribution of the length of the longest path of all the 
	 * connected components of size bigger than one, in all DAGs, and it sets a threshold 
	 * that removes the leftmost peak, which is assumed to be induced by noise. Such 
	 * threshold (the shortest length of a longest path of a connected component that is 
	 * not induced by noise) is used as an estimate of the shortest length of a longest 
	 * path (inside a connected component that is not induced by noise) that should be 
	 * considered a periodic substring. All DAGs use the same threshold.
	 */
	private static final void getSubstrings() {
		final double MIN_INTERVAL_LENGTH = 3.0;
		final double MIN_LOCAL_MAX_DISTANCE = 5.0;
		final int MIN_PERIOD = 1;
		final int MAX_PERIOD = Alignments.minAlignmentLength>>1;
		final int MIN_TRIMMED_LENGTH = Alignments.minAlignmentLength;  // Arbitrary
		int i, j, id;
		int firstAlignment, nAlignments;
		int startA1, endA1, startA2, endA2, startB1, endB1, startB2, endB2;
		int minStartA, maxStartA, minEndA, maxEndA, minStartB, maxEndB;
		int sumStartA, nStartA, sumEndA, nEndA;
		int readB, source, destination, pathLength, orientation;
		int nComponents, longestPathInComponent;

		// Ensuring the necessary order in $ReadA.sortedAlignments$
		if (Alignment.order!=Alignment.READB_ORIENTATION_STARTA_ENDA) {
			Alignment.order=Alignment.READB_ORIENTATION_STARTA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

		// Estimating $minPathLength$ from the connected components of all DAGs
		firstAlignment=0; lastLengthPointPrime=-1;
		do {
			readB=Alignments.alignments[ReadA.sortedAlignments[firstAlignment].id][1]-1;
			orientation=Alignments.alignments[ReadA.sortedAlignments[firstAlignment].id][2];
			i=firstAlignment+1;
			while ( i<=ReadA.lastSortedAlignment &&
			        Alignments.alignments[ReadA.sortedAlignments[i].id][1]-1==readB &&
			        Alignments.alignments[ReadA.sortedAlignments[i].id][2]==orientation ) i++;
			nAlignments=i-firstAlignment;
			buildDAG(firstAlignment,nAlignments,MIN_PERIOD,MAX_PERIOD,IDENTITY_THRESHOLD_PERIODIC);
			nComponents=DAG.getConnectedComponents(nAlignments,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stackPrime);
			if (nComponents<=0 || nComponents>nAlignments) {
				IO.printCriticalErr("Error while detecting periodic substrings: wrong number of connected components.");
				System.exit(1);
			}
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nAlignments);
			i=DAG.topologicalSort(nAlignments,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting periodic substrings: the DAG of all alignments of read "+ReadA.id+" with read "+readB+" contains a cycle that involves node "+(i-1));
				System.exit(1);
			}
			longestPathInComponent=DAG.longestPathInComponent(nAlignments,sorted2original,nInNeighbors,inNeighbors,components,lengthPoints,distancesPrime);
			if (longestPathInComponent>0) {
				for (i=0; i<nComponents; i++) {
					lastLengthPointPrime++;
					lengthPointsPrime[lastLengthPointPrime].position=lengthPoints[i].position;
					lengthPointsPrime[lastLengthPointPrime].mass=1;
				}
			}
			firstAlignment+=nAlignments;
		}
		while (firstAlignment<=ReadA.lastSortedAlignment);
		if (lastLengthPointPrime<=0) minPathLength=DEFAULT_MIN_PATH_LENGTH;
		else {
			minPathLength=estimateFromLengthPointsPrime(MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE);
			if (minPathLength<DEFAULT_MIN_PATH_LENGTH) minPathLength=DEFAULT_MIN_PATH_LENGTH;
		}		

		// Detecting periodic substrings
		firstAlignment=0; lastSubstring=-1;
		do {
			// Building the DAG of all alignments of readA with a specific readB, in a
			// specific orientation, and sorting it topologically.
			readB=Alignments.alignments[ReadA.sortedAlignments[firstAlignment].id][1]-1;
			orientation=Alignments.alignments[ReadA.sortedAlignments[firstAlignment].id][2];
			i=firstAlignment+1;
			while ( i<=ReadA.lastSortedAlignment &&
			        Alignments.alignments[ReadA.sortedAlignments[i].id][1]-1==readB &&
			        Alignments.alignments[ReadA.sortedAlignments[i].id][2]==orientation ) i++;
			nAlignments=i-firstAlignment;
			buildDAG(firstAlignment,nAlignments,MIN_PERIOD,MAX_PERIOD,IDENTITY_THRESHOLD_PERIODIC);
			nComponents=DAG.getConnectedComponents(nAlignments,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stackPrime);
			if (nComponents<=0 || nComponents>nAlignments) {
				IO.printCriticalErr("Error while detecting periodic substrings: wrong number of connected components.");
				System.exit(1);
			}
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nAlignments);
			i=DAG.topologicalSort(nAlignments,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting periodic substrings: the directed graph of all alignments between read "+ReadA.id+" and read "+readB+" contains a cycle that involves node "+(i-1));
				System.exit(1);
			}
			longestPathInComponent=DAG.longestPathInComponent(nAlignments,sorted2original,nInNeighbors,inNeighbors,components,lengthPoints,distancesPrime);
			if (longestPathInComponent<=0) {
				firstAlignment+=nAlignments;
				continue;
			}

			// Computing maximal candidate repeated periodic substrings in readA induced
			// by readB, considering every alignment of readA with readB as a possible
			// source alignment. The intervals of such substrings in readA can intersect.
			for (source=0; source<nAlignments; source++) {
				id=ReadA.sortedAlignments[firstAlignment+source].id;
				if ( ReadA.sortedAlignments[firstAlignment+source].impliedByPeriodicSubstring!=null ||
				     lengthPoints[components[source]].position<minPathLength ) continue;
				getSubstring(firstAlignment,source,nAlignments,minPathLength,IDENTITY_THRESHOLD_PERIODIC);
				destination=tmp[0]; pathLength=tmp[1];
				if (destination==-1) continue;
				getAlignmentsInLongestPath(source,destination);
				lastSubstring++;
				minStartA=Alignments.alignments[id][3];
				maxStartA=minStartA; nStartA=1; sumStartA=minStartA;
				minEndA=Alignments.alignments[id][4];
				maxEndA=minEndA; nEndA=0; sumEndA=0;
				minStartB=Alignments.alignments[id][5];
				maxEndB=Alignments.alignments[id][6];
				for (i=0; i<pathLength+1; i++) {
					id=ReadA.sortedAlignments[firstAlignment+alignments[i]].id;
					startA1=Alignments.alignments[id][3];
					endA1=Alignments.alignments[id][4];
					startB1=Alignments.alignments[id][5];
					endB1=Alignments.alignments[id][6];
					if (startA1<minStartA) minStartA=startA1;
					if (endA1>maxEndA) maxEndA=endA1;
					if (startB1<minStartB) minStartB=startB1;
					if (endB1>maxEndB) maxEndB=endB1;
					if (alignmentStates[i]==0 || alignmentStates[i]==2 || alignmentStates[i]==6) {
						if (startA1>maxStartA) maxStartA=startA1;
						nStartA++; sumStartA+=startA1;
					}
					if (alignmentStates[i]==2 || alignmentStates[i]==3 || alignmentStates[i]==5 || alignmentStates[i]==6 || alignmentStates[i]==7 || alignmentStates[i]==8 || alignments[i]==destination) {
						if (endA1<minEndA) minEndA=endA1;
						nEndA++; sumEndA+=endA1;
					}
				}
				substrings[lastSubstring].hasLongPeriod=false;
				substrings[lastSubstring].minStartA=minStartA;
				substrings[lastSubstring].maxStartA=maxStartA;
				substrings[lastSubstring].nStartA=nStartA;
				substrings[lastSubstring].sumStartA=sumStartA;
				substrings[lastSubstring].minEndA=minEndA;
				substrings[lastSubstring].maxEndA=maxEndA;
				substrings[lastSubstring].nEndA=nEndA;
				substrings[lastSubstring].sumEndA=sumEndA;
				substrings[lastSubstring].minStartBForward=minStartB;
				substrings[lastSubstring].maxEndBForward=maxEndB;
				substrings[lastSubstring].readB=readB;
				substrings[lastSubstring].pathLength=pathLength;
				if (orientation==1) {
					substrings[lastSubstring].orientation=true;
					substrings[lastSubstring].minStartB=substrings[lastSubstring].minStartBForward;
					substrings[lastSubstring].maxEndB=substrings[lastSubstring].maxEndBForward;
				}
				else {
					substrings[lastSubstring].orientation=false;
					substrings[lastSubstring].minStartB=Reads.getReadLength(readB)-substrings[lastSubstring].maxEndBForward-1;
					substrings[lastSubstring].maxEndB=Reads.getReadLength(readB)-substrings[lastSubstring].minStartBForward-1;
				}
				collectShifts(firstAlignment,substrings[lastSubstring]);
				markImpliedAlignments(firstAlignment,source,destination,nAlignments,lastSubstring,IDENTITY_THRESHOLD_PERIODIC);
				substrings[lastSubstring].trim(tmp,MIN_TRIMMED_LENGTH);
			}
			firstAlignment+=nAlignments;
		}
		while (firstAlignment<=ReadA.lastSortedAlignment);
		if (lastSubstring>=0) {
			discardAlignmentsAfterTrimming(0,ReadA.lastSortedAlignment,false);
			PeriodicSubstring.order=PeriodicSubstring.READB_ORIENTATION_MINSTARTA;
			if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		}
		

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("OUTPUT OF GETSUBSTRINGS PERIODIC: (lastSubstring="+lastSubstring+")");
	for (int x=0; x<=lastSubstring; x++) IO.printErr(substrings[x].hashCode()+" :: "+substrings[x]); 
	IO.printErr("alignments:");
Alignment.order=Alignment.STARTA_ENDA;
Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);	
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
}

	}


	/**
	 * Builds the DAG of periodic substrings using all alignments in
	 * $ReadA.sortedAlignments[firstAlignment..firstAlignment+nAlignments-1]$, which are 
	 * supposed to use the same readB in the same orientation.
	 *
	 * @param minPeriod,maxPeriod alignments whose starting positions in readA or readB 
	 * are not at a distance in $[minPeriod..maxPeriod]$ are not connected;
	 * @param identityThresholdPeriodic to decide identity of two positions.
	 */
	private static final void buildDAG(int firstAlignment, int nAlignments, int minPeriod, int maxPeriod, int identityThresholdPeriodic) {
		boolean orientation;
		int i, j;
		int id, startA1, endA1, startB1, endB1, startA2, endA2, startB2, endB2;
		int longestA, shortestA, longestB, shortestB;

		clearNeighborMatrices(nAlignments);
		for (i=0; i<nAlignments; i++) {
			id=ReadA.sortedAlignments[firstAlignment+i].id;
			startA1=Alignments.alignments[id][3];
			endA1=Alignments.alignments[id][4];
			startB1=Alignments.alignments[id][5];
			endB1=Alignments.alignments[id][6];
			orientation=Alignments.alignments[id][2]==1;
			for (j=0; j<nAlignments; j++) {
				if (j==i) continue;
				id=ReadA.sortedAlignments[firstAlignment+j].id;
				startA2=Alignments.alignments[id][3];
				endA2=Alignments.alignments[id][4];
				startB2=Alignments.alignments[id][5];
				endB2=Alignments.alignments[id][6];
				if ( Intervals.intersectionLength(startA1,endA1,startA2,endA2)<=IO.quantum ||
					 Intervals.intersectionLength(startB1,endB1,startB2,endB2)<=IO.quantum ) {
					// Only alignments that intersect both in readA and in readB can be
					// connected.
					continue;  
				}				
				if ( (Math.abs(startA2,startA1)>identityThresholdPeriodic && Math.abs(startA2,startA1)<minPeriod) || 
					 Math.abs(startA2,startA1)>maxPeriod || 
					 (Math.abs(startB2,startB1)>identityThresholdPeriodic && Math.abs(startB2,startB1)<minPeriod) || 
					 Math.abs(startB2,startB1)>maxPeriod ) continue;
				longestA=Math.max(endA1-startA1,endA2-startA2)+1;
				shortestA=Math.min(endA1-startA1,endA2-startA2)+1;
				longestB=Math.max(endB1-startB1,endB2-startB2)+1;
				shortestB=Math.min(endB1-startB1,endB2-startB2)+1;
				// Original states
				if ( Math.abs(startA2,startA1)<=identityThresholdPeriodic &&
					 endA2-endA1>identityThresholdPeriodic &&
					 (orientation?startB1-startB2>identityThresholdPeriodic:endB2-endB1>identityThresholdPeriodic) &&
					 (orientation?Math.abs(endB2,endB1)<=identityThresholdPeriodic:Math.abs(startB2,startB1)<=identityThresholdPeriodic) ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=0;  // prefix to prefix and prefix to interior.
					nInNeighbors[j]++;
				}
				else if ( startA2-startA1>identityThresholdPeriodic &&
						  endA2-endA1>identityThresholdPeriodic &&
						  Math.abs(endB2,endB1)<=identityThresholdPeriodic &&
						  Math.abs(startB2,startB1)<=identityThresholdPeriodic &&
						  longestA<=maxLengthRatio*shortestA ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=1;  // $j'-i'<j-i$, interior to interior.
					nInNeighbors[j]++;
				}
				else if ( Math.abs(startA2,startA1)<=identityThresholdPeriodic &&
						  Math.abs(endA2,endA1)<=identityThresholdPeriodic &&
						  (orientation?endB1-endB2>identityThresholdPeriodic:startB2-startB1>identityThresholdPeriodic) &&
						  (orientation?startB1-startB2>identityThresholdPeriodic:endB2-endB1>identityThresholdPeriodic) &&
						  longestB<=maxLengthRatio*shortestB ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=2;  // $j'-i'>j-i$, whole to whole.
					nInNeighbors[j]++;
				}
				else if ( startA2-startA1>identityThresholdPeriodic &&
						  Math.abs(endA2,endA1)<=identityThresholdPeriodic &&
						  (orientation?Math.abs(startB2,startB1)<=identityThresholdPeriodic:Math.abs(endB2,endB1)<=identityThresholdPeriodic) &&
						  (orientation?endB1-endB2>identityThresholdPeriodic:startB2-startB1>identityThresholdPeriodic) ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=3;  // interior to suffix and suffix to suffix.
					nInNeighbors[j]++;
				}
				// States added to handle phase difference
				else if ( startA2-startA1>identityThresholdPeriodic &&
						  endA2-endA1>identityThresholdPeriodic &&
						  (orientation?startB1-startB2>identityThresholdPeriodic:endB2-endB1>identityThresholdPeriodic) &&
						  (orientation?Math.abs(endB2,endB1)<=identityThresholdPeriodic:Math.abs(startB2,startB1)<=identityThresholdPeriodic) ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=4;  // $j'-i'<j-i$, prefix to interior.
					nInNeighbors[j]++;
				}
				else if ( startA2-startA1>identityThresholdPeriodic &&
						  endA2-endA1>identityThresholdPeriodic &&
						  (orientation?Math.abs(startB2,startB1)<=identityThresholdPeriodic:Math.abs(endB2,endB1)<=identityThresholdPeriodic) &&
						  (orientation?endB1-endB2>identityThresholdPeriodic:startB2-startB1>identityThresholdPeriodic) ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=5;  // $j'-i'<j-i$, interior to suffix.
					nInNeighbors[j]++;
				}
				else if ( Math.abs(startA2,startA1)<=identityThresholdPeriodic &&
						  endA2-endA1>identityThresholdPeriodic &&
						  (orientation?startB1-startB2>identityThresholdPeriodic:endB2-endB1>identityThresholdPeriodic) &&
						  (orientation?endB1-endB2>identityThresholdPeriodic:startB2-startB1>identityThresholdPeriodic) ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=6;  // $j'-i'>j-i$, prefix to whole.
					nInNeighbors[j]++;
				}
				else if ( startA2-startA1>identityThresholdPeriodic &&
						  Math.abs(endA2,endA1)<=identityThresholdPeriodic &&
						  (orientation?startB1-startB2>identityThresholdPeriodic:endB2-endB1>identityThresholdPeriodic) &&
						  (orientation?endB1-endB2>identityThresholdPeriodic:startB2-startB1>identityThresholdPeriodic) ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=7;  // $j'-i'>j-i$, whole to suffix.
					nInNeighbors[j]++;
				}
				else if ( startA2-startA1>identityThresholdPeriodic &&
						  endA2-endA1>identityThresholdPeriodic &&
						  (orientation?endB1-endB2>identityThresholdPeriodic:startB2-startB1>identityThresholdPeriodic) &&
						  (orientation?startB1-startB2>identityThresholdPeriodic:endB2-endB1>identityThresholdPeriodic) &&
						  longestA<=maxLengthRatio*shortestA ) {
					Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
					outNeighbors[i][nOutNeighbors[i]++]=j;
					Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+1);
					inNeighbors[j][nInNeighbors[j]]=i;
					Math.ensureSpace(inNeighborStates,j,nInNeighbors[j]+1);
					inNeighborStates[j][nInNeighbors[j]]=8;  // $j'-i'=j-i$, prefix to suffix.
					nInNeighbors[j]++;
				}
			}
		}
	}


	/**
	 * Fits a density estimation tree on $lengthPointsPrime$, and returns a separation 
	 * point between the first local maximum and the rest (the first local maximum is 
	 * assumed to correspond to connected components that arise by chance). Since 
	 * $lengthPointsPrime$ contains longest-path lengths from all connected components of 
	 * all reads, and since the longest-path length is proportional to the sum of the 
	 * lengths of a periodic substring in readA and readB, we expect $lengthPointsPrime$ 
	 * to be a continuous function with just one true peak around one, and no obvious 
	 * separation point between such peak and the rest.
	 *
	 * Remark: setting the threshold is a very delicate process. Missing a periodic 
	 * substring element affects the coverage estimate, and might introduce spurious 
	 * events, since the corresponding alignments might not be marked as implied.
	 */
	private static final int estimateFromLengthPointsPrime(double minIntervalLength, double minLocalMaxDistance) {
		final int MAX_PATH_LENGTH = 10;  // Arbitrary
		int tmp, nLocalMaximumLeaves, out;

		lastLengthPointPrime=Points.sortAndCompact(lengthPointsPrime,lastLengthPointPrime);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("estimateFromLengthPointsPrime:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=lastLengthPointPrime; x++) IO.printErr(lengthPointsPrime[x].position+","+lengthPointsPrime[x].getMass()); }

		if (lastLengthPointPrime+1<=Points.FEW_POINTS) {
			if (lengthPointsPrime[lastLengthPointPrime].position-lengthPointsPrime[0].position<=MAX_PATHLENGTH_DIFFERENCE) return DEFAULT_MIN_PATH_LENGTH;
			tmp=Points.getRoughThreshold(lengthPointsPrime,lastLengthPointPrime,true,DEFAULT_MIN_PATH_LENGTH,true);
			return (tmp==-1||lengthPointsPrime[tmp].position>MAX_PATH_LENGTH)?DEFAULT_MIN_PATH_LENGTH:(int)lengthPointsPrime[tmp].position;
		}
		if (Points.areUniformlyDistributed(lengthPointsPrime,0,lastLengthPointPrime,true,(lengthPointsPrime[lastLengthPointPrime].position-lengthPointsPrime[0].position)/Points.DEFAULT_NBINS)) return DEFAULT_MIN_PATH_LENGTH;
		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(lengthPointsPrime,0,lastLengthPointPrime,minIntervalLength,minLocalMaxDistance,true,-1,-1,-1,false,true);
		if (nLocalMaximumLeaves==0) return DEFAULT_MIN_PATH_LENGTH;
		DensityEstimationTree.markRunsOfLocalMaximumLeaves(lengthPointsPrime);


if (IO.SHOW_STD_ERR_PRIME) IO.printErr("DET:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr(DensityEstimationTree.leaves[x]); }

		out=(int)lengthPointsPrime[DensityEstimationTree.leaves[Leaves.lastLeafOfRun(DensityEstimationTree.leaves,DensityEstimationTree.lastLeaf,0)].lastPoint].position;
		return (out<DEFAULT_MIN_PATH_LENGTH||out>MAX_PATH_LENGTH)?DEFAULT_MIN_PATH_LENGTH:out;
	}


	/**
	 * The procedure sets $impliedByPeriodicSubstring$ to $substrings[substring]$ for all
	 * alignments in the longest path from $source$ to $destination$, using $alignments$.
	 * Then, it scans all the remaining alignments that lie between $source$ and 
	 * $destination$ in the order of startA: an alignment that intersects the readA and 
	 * readB intervals of $substrings[substring]$ with the same orientation is assigned 
	 * $impliedByPeriodicSubstring$. However, an alignment that intersects just the readA 
	 * interval of $substrings[substring]$, and matches it to a different interval of 
	 * readB, or to the same interval of readB in the opposite orientation, is not marked.
     *
	 * The alignments that are marked by this procedure will be discarded by the following
	 * steps of the pipeline.
	 *
	 * Remark: the procedure assumes that there is no peak in the distribution of lengths
	 * or positions of alignments, since this is not likely to happen (peaks in alignment
	 * quality are more likely).
	 *
	 * Remark: implied alignments are used to update the boundaries of
	 * $substrings[substring]$.
	 */
	private static final void markImpliedAlignments(int firstAlignment, int source, int destination, int nAlignments, int substring, int identityThresholdPeriodic) {
		final int MAXIMALITY_THRESHOLD = IO.quantum;
		int i, vertex;
		int id, first, last, newMinStartA, newMaxEndA, newMinStartB, newMaxEndB;
		int startA, endA, startB, endB, readB;
		int state, maxDistance, max;
		int sourceID, sourceReadB, sourceOrientation;

		// Marking all alignments in the longest path from $source$ to $destination$
if (IO.SHOW_STD_ERR) IO.printErr("Longest-path alignments:");
		for (i=0; i<=lastAlignment; i++) {
			if (ReadA.sortedAlignments[firstAlignment+alignments[i]].impliedByPeriodicSubstring==null) {
				ReadA.sortedAlignments[firstAlignment+alignments[i]].impliedByPeriodicSubstring=substrings[substring];
if (IO.SHOW_STD_ERR) IO.printErr(ReadA.sortedAlignments[firstAlignment+alignments[i]].toStringBoundaries());
			}
		}

		// Deciding the first and the last candidate alignment for further markings
		max=Math.max(IO.quantum,identityThresholdPeriodic);
		sourceID=ReadA.sortedAlignments[firstAlignment+source].id;
		sourceReadB=Alignments.alignments[sourceID][1];
		sourceOrientation=Alignments.alignments[sourceID][2];
		for (i=source-1; i>=0; i--) {
			id=ReadA.sortedAlignments[firstAlignment+i].id;
			if (Alignments.alignments[id][1]!=sourceReadB || Alignments.alignments[id][2]!=sourceOrientation || Alignments.alignments[id][3]<substrings[substring].minStartA-max) break;
		}
		first=i+1;
		for (i=destination+1; i<nAlignments; i++) {
			id=ReadA.sortedAlignments[firstAlignment+i].id;
			if (Alignments.alignments[id][1]!=sourceReadB || Alignments.alignments[id][2]!=sourceOrientation || Alignments.alignments[id][3]>substrings[substring].maxEndA) break;
		}
		last=i-1;

		// Marking the remaining alignments implied by $substrings[substring]$
		newMinStartA=substrings[substring].minStartA;
		newMaxEndA=substrings[substring].maxEndA;
		newMinStartB=substrings[substring].minStartBForward;
		newMaxEndB=substrings[substring].maxEndBForward;
		for (i=first; i<=last; i++) {
			id=ReadA.sortedAlignments[firstAlignment+i].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			startB=Alignments.alignments[id][5];
			endB=Alignments.alignments[id][6];
			readB=Alignments.alignments[id][1];
			if ( startA>=substrings[substring].minStartA-identityThresholdPeriodic &&
				 endA<=substrings[substring].maxEndA+identityThresholdPeriodic &&
				 startB>=substrings[substring].minStartBForward-identityThresholdPeriodic &&
				 endB<=substrings[substring].maxEndBForward+identityThresholdPeriodic &&
				 ( Intervals.isApproximatelyContained(startA,endA,substrings[substring].minStartA,substrings[substring].maxEndA) &&
				   !( (startA>substrings[substring].minStartA+identityThresholdPeriodic && (!Reads.isLeftMaximal(startA,ReadA.id,true) || ReadA.sortedAlignments[firstAlignment+i].lowQualityStart)) &&
					  (endA<substrings[substring].maxEndA-identityThresholdPeriodic && (!Reads.isRightMaximal(endA,ReadA.id,true) || ReadA.sortedAlignments[firstAlignment+i].lowQualityEnd))
				   )
				 ) &&
				 ( Intervals.isApproximatelyContained(startB,endB,substrings[substring].minStartBForward,substrings[substring].maxEndBForward) &&
   				   !( (startB>substrings[substring].minStartBForward+identityThresholdPeriodic && !Reads.isLeftMaximal(startB,readB-1,true)) &&
   					  (endB<substrings[substring].maxEndBForward-identityThresholdPeriodic && !Reads.isRightMaximal(endB,readB-1,true))
   				   )	  
				 )
			   ) {
				if (ReadA.sortedAlignments[firstAlignment+i].impliedByPeriodicSubstring==null) ReadA.sortedAlignments[firstAlignment+i].impliedByPeriodicSubstring=substrings[substring];
				// Using an implied alignment to update the bounds of the containing
				// periodic substring
				if (startA<newMinStartA) newMinStartA=startA;
				if (endA>newMaxEndA) newMaxEndA=endA;
				if (startB<newMinStartB) newMinStartB=startB;
				if (endB>newMaxEndB) newMaxEndB=endB;
			}
		}

		// Updating $substrings[substring]$
		substrings[substring].minStartA=newMinStartA;
		substrings[substring].maxEndA=newMaxEndA;
		substrings[substring].minStartBForward=newMinStartB;
		substrings[substring].maxEndBForward=newMaxEndB;
		if (substrings[substring].orientation) {
			substrings[substring].minStartB=substrings[substring].minStartBForward;
			substrings[substring].maxEndB=substrings[substring].maxEndBForward;
		}
		else {
			substrings[substring].minStartB=Reads.getReadLength(substrings[substring].readB)-substrings[substring].maxEndBForward-1;
			substrings[substring].maxEndB=Reads.getReadLength(substrings[substring].readB)-substrings[substring].minStartBForward-1;
		}
	}
	
	
	/**
	 * Detects whether a periodic substring interval is adjacent to a high-quality 
	 * substring to the left and to the right in its readA.
	 */
	private static final void setIntervalMaximality(PeriodicSubstringInterval[] intervals, int lastInterval) {
		for (int i=0; i<=lastInterval; i++) {
			intervals[i].isLeftMaximal=Reads.isLeftMaximal(intervals[i].firstPosition,ReadA.id,true);
			intervals[i].isRightMaximal=Reads.isRightMaximal(intervals[i].lastPosition,ReadA.id,true);
		}
	}


	/**
	 * Searches for a maximal repeated periodic substring of $readA$ that starts at a
	 * specific alignment $v$ with a $readB$. This is done by computing the longest path,
	 * of length at least $minPathLength$, with states 0^i->3^j, 0^i->1^j->3^k, or
	 * 0^i->2^j->3^k between vertex $v$ and every other vertex in the DAG built by
	 * procedure $getSubstrings$.
	 *
	 * We assume that any subset of the set of all alignments of a repeated periodic
	 * substring of $readA$ with its corresponding substring of $readB$ can be missing.
	 * Thus, we don't impose the longest path to reach a vertex with state 3, but we
	 * consider every possible state.
	 *
	 * Remark: we handle the case in which $R[i..j]$ is not in the same phase as
	 * $R'[i'..j']$, using the additional states introduced by procedure $getSubstrings$.
	 * See the <a href="figures/automaton.png">diagram of valid state transitions</a>.
	 *
	 * Remark: we could add specifications on whether two alignments should be connected,
	 * based on their left- and right-maximality in readA and in readB, considered
	 * separately. This is likely to make the diagram of valid transitions and the code
	 * more complex, without a significant benefit in accuracy.
	 *
	 * Remark: we don't expect the path returned by this procedure to cross two adjacent
	 * repeated periodic substrings with distinct periods.
	 *
	 * Remark: this is just a heuristic, which maximizes just the surface of a periodic
	 * substring in readA, and which does not necessarily maximize the surface in readB.
	 * Rather than longest paths, we should find \emph{maximal} paths in the DAG, whose
	 * labels conform to the <a href="figures/automaton.png">diagram of valid state
	 * transitions</a>.
	 *
	 * @param firstAlignment first alignment of $readA$ with a specific $readB$ in
	 * $ReadA.sortedAlignments$;
	 * @param source relative position of alignment $v$ in the list of all the alignments
	 * of $readA$ with a specific $readB$;
	 * @param nVertices number of alignments of $readA$ with a specific $readB$;
	 * @return $tmp[0]$: the index (relative to $firstAlignment$) of an interval $w$ that
	 * maximizes the distance from $source$ in the DAG, and that ends the farthest from
	 * $source$, or -1 if no valid path from $source$ is found; $tmp[1]$: the length
	 * of a longest path from $source$ in the DAG.
	 */
	private static final void getSubstring(int firstAlignment, int source, int nVertices, int minPathLength, int identityThresholdPeriodic) {
		boolean orientation;
		int i, j, k, s;
		int id, vertex, rightmostAlignment, state, predecessorState, distance, longestDistance, end, rightmostEnd;
		int sourceStartA, sourceEndA, sourceStartB, sourceEndB, sourceComponent, intervalLengthA, intervalLengthB;
		int lastSortedVertex;

		// Computing the longest distance from $source$ to every other vertex
		for (i=0; i<nVertices; i++) {
			Math.set(distances[i],N_STATES-1,-1);
			Math.set(predecessors[i],N_STATES-1,-1);
			Math.set(predecessorStates[i],N_STATES-1,-1);
			Math.set(maxIntervalLengthA[i],N_STATES-1,-1);
			Math.set(minIntervalLengthA[i],N_STATES-1,-1);
			Math.set(maxIntervalLengthB[i],N_STATES-1,-1);
			Math.set(minIntervalLengthB[i],N_STATES-1,-1);
		}
		id=ReadA.sortedAlignments[firstAlignment+source].id;
		sourceStartA=Alignments.alignments[id][3];
		sourceEndA=Alignments.alignments[id][4];
		sourceStartB=Alignments.alignments[id][5];
		sourceEndB=Alignments.alignments[id][6];
		orientation=Alignments.alignments[id][2]==1;
		sourceComponent=components[source];
		maxIntervalLengthA[source][0]=sourceEndA-sourceStartA+1;
		minIntervalLengthA[source][0]=maxIntervalLengthA[source][0];
		distances[source][0]=0;
		s=original2sorted[source];
		lastSortedVertex=s;
if (IO.SHOW_STD_ERR) IO.printErr("source="+(firstAlignment+source)+" => s="+s+" nVertices="+nVertices);
		for (j=s+1; j<nVertices; j++) {
			vertex=sorted2original[j];
			if (components[vertex]!=sourceComponent) {
				lastSortedVertex=j-1;
				break;
			}
if (IO.SHOW_STD_ERR) IO.printErr("--> considering alignment "+(firstAlignment+vertex));
			id=ReadA.sortedAlignments[firstAlignment+vertex].id;
			intervalLengthA=Alignments.alignments[id][4]-Alignments.alignments[id][3]+1;
			intervalLengthB=Alignments.alignments[id][6]-Alignments.alignments[id][5]+1;
			for (i=0; i<nInNeighbors[vertex]; i++) {
				state=inNeighborStates[vertex][i];
				k=inNeighbors[vertex][i];
				distance=-1;
				// Original states
				if (state==0) {
					if ( Math.abs(Alignments.alignments[id][3],sourceStartA)>identityThresholdPeriodic ||
					     (orientation?Math.abs(Alignments.alignments[id][6],sourceEndB)>identityThresholdPeriodic:Math.abs(Alignments.alignments[id][5],sourceStartB)>identityThresholdPeriodic) ) continue;
					if (distances[k][0]!=-1) distance=distances[k][0];
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][0]) {
						distances[vertex][0]=distance;
						predecessors[vertex][0]=k;
						predecessorStates[vertex][0]=0;
						minIntervalLengthA[vertex][0]=intervalLengthA;
						maxIntervalLengthA[vertex][0]=intervalLengthA;
						minIntervalLengthB[vertex][0]=intervalLengthB;
						maxIntervalLengthB[vertex][0]=intervalLengthB;
					}
				}
				else if (state==1) {
					if (orientation?Math.abs(Alignments.alignments[id][6],sourceEndB)>identityThresholdPeriodic:Math.abs(Alignments.alignments[id][5],sourceStartB)>identityThresholdPeriodic) continue;
					predecessorState=-1;
					if ( distances[k][0]!=-1 &&
					     intervalLengthA<=maxLengthRatio*minIntervalLengthA[k][0] &&
					     maxIntervalLengthA[k][0]<=maxLengthRatio*intervalLengthA &&
					     intervalLengthB<=maxLengthRatio*minIntervalLengthB[k][0] &&
					     maxIntervalLengthB[k][0]<=maxLengthRatio*intervalLengthB ) {
						distance=distances[k][0];
						predecessorState=0;
					}
					if ( distances[k][1]!=-1 &&
					     distances[k][1]>distance &&
					     intervalLengthA<=maxLengthRatio*minIntervalLengthA[k][1] &&
					     maxIntervalLengthA[k][1]<=maxLengthRatio*intervalLengthA &&
					     intervalLengthB<=maxLengthRatio*minIntervalLengthB[k][1] &&
					     maxIntervalLengthB[k][1]<=maxLengthRatio*intervalLengthB ) {
						distance=distances[k][1];
						predecessorState=1;
					}
					if ( distances[k][4]!=-1 &&
					     distances[k][4]>distance &&
					     intervalLengthA<=maxLengthRatio*minIntervalLengthA[k][4] &&
					     maxIntervalLengthA[k][4]<=maxLengthRatio*intervalLengthA &&
					     intervalLengthB<=maxLengthRatio*minIntervalLengthB[k][4] &&
					     maxIntervalLengthB[k][4]<=maxLengthRatio*intervalLengthB ) {
						distance=distances[k][4];
						predecessorState=4;
					}
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][1]) {
						distances[vertex][1]=distance;
						predecessors[vertex][1]=k;
						predecessorStates[vertex][1]=predecessorState;
						minIntervalLengthA[vertex][1]=Math.min(intervalLengthA,minIntervalLengthA[k][predecessorState]);
						maxIntervalLengthA[vertex][1]=Math.max(intervalLengthA,maxIntervalLengthA[k][predecessorState]);
						minIntervalLengthB[vertex][1]=Math.min(intervalLengthB,minIntervalLengthB[k][predecessorState]);
						maxIntervalLengthB[vertex][1]=Math.max(intervalLengthB,maxIntervalLengthB[k][predecessorState]);
					}
				}
				else if (state==2) {
					if (Math.abs(Alignments.alignments[id][3],sourceStartA)>identityThresholdPeriodic) continue;
					predecessorState=-1;
					if ( distances[k][0]!=-1 &&
					     intervalLengthA<=maxLengthRatio*minIntervalLengthA[k][0] &&
					     maxIntervalLengthA[k][0]<=maxLengthRatio*intervalLengthA &&
					     intervalLengthB<=maxLengthRatio*minIntervalLengthB[k][0] &&
					     maxIntervalLengthB[k][0]<=maxLengthRatio*intervalLengthB ) {
						distance=distances[k][0];
						predecessorState=0;
					}
					if ( distances[k][2]!=-1 &&
					     distances[k][2]>distance &&
					     intervalLengthA<=maxLengthRatio*minIntervalLengthA[k][2] &&
					     maxIntervalLengthA[k][2]<=maxLengthRatio*intervalLengthA &&
					     intervalLengthB<=maxLengthRatio*minIntervalLengthB[k][2] &&
					     maxIntervalLengthB[k][2]<=maxLengthRatio*intervalLengthB ) {
						distance=distances[k][2];
						predecessorState=2;
					}
					if ( distances[k][6]!=-1 &&
					     distances[k][6]>distance &&
					     intervalLengthA<=maxLengthRatio*minIntervalLengthA[k][6] &&
					     maxIntervalLengthA[k][6]<=maxLengthRatio*intervalLengthA &&
					     intervalLengthB<=maxLengthRatio*minIntervalLengthB[k][6] &&
					     maxIntervalLengthB[k][6]<=maxLengthRatio*intervalLengthB ) {
						distance=distances[k][6];
						predecessorState=6;
					}
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][2]) {
						distances[vertex][2]=distance;
						predecessors[vertex][2]=k;
						predecessorStates[vertex][2]=predecessorState;
						minIntervalLengthA[vertex][2]=Math.min(intervalLengthA,minIntervalLengthA[k][predecessorState]);
						maxIntervalLengthA[vertex][2]=Math.max(intervalLengthA,maxIntervalLengthA[k][predecessorState]);
						minIntervalLengthB[vertex][2]=Math.min(intervalLengthB,minIntervalLengthB[k][predecessorState]);
						maxIntervalLengthB[vertex][2]=Math.max(intervalLengthB,maxIntervalLengthB[k][predecessorState]);
					}
				}
				else if (state==3) {
					predecessorState=-1;
					if (distances[k][0]!=-1) {
						distance=distances[k][0];
						predecessorState=0;
					}
					if (distances[k][1]!=-1 && distances[k][1]>distance) {
						distance=distances[k][1];
						predecessorState=1;
					}
					if (distances[k][2]!=-1 && distances[k][2]>distance) {
						distance=distances[k][2];
						predecessorState=2;
					}
					if (distances[k][3]!=-1 && distances[k][3]>distance) {
						distance=distances[k][3];
						predecessorState=3;
					}
					if (distances[k][5]!=-1 && distances[k][5]>distance) {
						distance=distances[k][5];
						predecessorState=5;
					}
					if (distances[k][7]!=-1 && distances[k][7]>distance) {
						distance=distances[k][7];
						predecessorState=7;
					}
					if (distances[k][8]!=-1 && distances[k][8]>distance) {
						distance=distances[k][8];
						predecessorState=8;
					}
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][3]) {
						distances[vertex][3]=distance;
						predecessors[vertex][3]=k;
						predecessorStates[vertex][3]=predecessorState;
					}
				}
				// States added to handle phase difference
				else if (state==4) {
					if (orientation?Math.abs(Alignments.alignments[id][6],sourceEndB)>identityThresholdPeriodic:Math.abs(Alignments.alignments[id][5],sourceStartB)>identityThresholdPeriodic) continue;
					if (distances[k][0]!=-1) distance=distances[k][0];
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][4]) {
						distances[vertex][4]=distance;
						predecessors[vertex][4]=k;
						predecessorStates[vertex][4]=0;
						minIntervalLengthA[vertex][4]=intervalLengthA;
						maxIntervalLengthA[vertex][4]=intervalLengthA;
						minIntervalLengthB[vertex][4]=intervalLengthB;
						maxIntervalLengthB[vertex][4]=intervalLengthB;
					}
				}
				else if (state==5) {
					predecessorState=-1;
					if (distances[k][1]!=-1 && distances[k][1]>distance) {
						distance=distances[k][1];
						predecessorState=1;
					}
					if (distances[k][4]!=-1 && distances[k][4]>distance) {
						distance=distances[k][4];
						predecessorState=4;
					}
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][5]) {
						distances[vertex][5]=distance;
						predecessors[vertex][5]=k;
						predecessorStates[vertex][5]=predecessorState;
					}
				}
				else if (state==6) {
					if (Math.abs(Alignments.alignments[id][3],sourceStartA)>identityThresholdPeriodic) continue;
					if (distances[k][0]!=-1) distance=distances[k][0];
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][6]) {
						distances[vertex][6]=distance;
						predecessors[vertex][6]=k;
						predecessorStates[vertex][6]=0;
						minIntervalLengthA[vertex][6]=intervalLengthA;
						maxIntervalLengthA[vertex][6]=intervalLengthA;
						minIntervalLengthB[vertex][6]=intervalLengthB;
						maxIntervalLengthB[vertex][6]=intervalLengthB;
					}
				}
				else if (state==7) {
					predecessorState=-1;
					if (distances[k][2]!=-1) {
						distance=distances[k][2];
						predecessorState=2;
					}
					if (distances[k][6]!=-1 && distances[k][6]>distance) {
						distance=distances[k][6];
						predecessorState=6;
					}
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][7]) {
						distances[vertex][7]=distance;
						predecessors[vertex][7]=k;
						predecessorStates[vertex][7]=predecessorState;
					}
				}
				else if (state==8) {
					if (distances[k][0]!=-1) distance=distances[k][0];
					if (distance==-1) continue;
					distance++;
					if (distance>distances[vertex][8]) {
						distances[vertex][8]=distance;
						predecessors[vertex][8]=k;
						predecessorStates[vertex][8]=0;
					}
				}
			}
		}
		if (j==nVertices) lastSortedVertex=j-1;

		// Computing the rightmost alignment, of any state, whose distance from $source$
		// in the graph is maximum.
		longestDistance=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			for (j=0; j<N_STATES; j++) {
				if (distances[vertex][j]>longestDistance) longestDistance=distances[vertex][j];
			}
		}
		if (longestDistance+1<minPathLength) {  // No path from $source$ with the required structure
			Math.set(tmp,1,-1);
			return;
		}
		tmp[1]=longestDistance;
		rightmostEnd=-1; rightmostAlignment=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			id=ReadA.sortedAlignments[firstAlignment+vertex].id;
			for (j=0; j<N_STATES; j++) {
				if (distances[vertex][j]<longestDistance) continue;
				end=Alignments.alignments[id][4];
				if (end>rightmostEnd) {
					rightmostAlignment=vertex;
					rightmostEnd=end;
				}
			}
		}
		tmp[0]=rightmostAlignment;
	}


	/**
	 * Stores in $alignments$ and $alignmentStates$ the set of all alignments and
	 * corresponding states in the longest path from $source$ to $destination$ encoded by
	 * $predecessors$. Alignments are sorted from $destination$ to $source$.
	 *
	 * Remark: the values stored in $alignments$ are relative to $firstAlignment$.
	 */
	private static final void getAlignmentsInLongestPath(int source, int destination) {
		int i;
		int alignment, nextAlignment, state, maxDistance;


		alignment=destination;
		state=-1;
		maxDistance=-1;
		for (i=0; i<N_STATES; i++) {
			if (distances[destination][i]>maxDistance) {
				maxDistance=distances[destination][i];
				state=i;
			}
		}

		lastAlignment=-1;
		do {
			lastAlignment++;
			alignments[lastAlignment]=alignment;
			alignmentStates[lastAlignment]=state;
			nextAlignment=predecessors[alignment][state];
			state=predecessorStates[alignment][state];
			alignment=nextAlignment;
		} while (alignment!=-1);

		if (alignments[lastAlignment]!=source) {
			IO.printCriticalErr("Error in procedure $getSubstring$: the longest path from the source vertex does not start at the source vertex.");
			System.exit(1);
		}
		if (lastAlignment!=maxDistance) {
			IO.printCriticalErr("Error in procedure $getSubstring$: the longest path from the source vertex to the destination vertex contains "+(lastAlignment+1)+" alignments, rather than the claimed "+(maxDistance+1)+".");
			System.exit(1);
		}
	}


	/**
	 * Stores in $substring.shifts$ all distinct shifts between consecutive alignments in
	 * $alignments$. To save space, the procedure collects all (not necessarily 
	 * distinct) observed shifts in temporary array $shifts$ first.
	 *
	 * Remark: the procedure assumes that all ${min,max}{Start,End}{A,B}$ properties of 
	 * $substring$ have already been set.
	 *
	 * @param firstAlignment the first alignment between $readA$ and $readB$ in
	 * $ReadA.sortedAlignments$.
	 */
	private static final void collectShifts(int firstAlignment, PeriodicSubstring substring) {
		final int THRESHOLD = IO.quantum;
		boolean orientation;
		int i, j;
		int id, shift;
		int startA, endA, startB, endB;
		int previousStartA, previousEndA, previousStartB, previousEndB, isLeftMaximal, isRightMaximal;
		Alignment alignment;
		
		alignment=ReadA.sortedAlignments[firstAlignment+alignments[0]];
		id=alignment.id;
		orientation=Alignments.alignments[id][2]==1;
		isLeftMaximal=alignment.isLeftMaximal;
		isRightMaximal=alignment.isRightMaximal;
		startA=Alignments.alignments[id][3];
		endA=Alignments.alignments[id][4];
		startB=Alignments.alignments[id][5];
		endB=Alignments.alignments[id][6];
		if ( isLeftMaximal==1 || 
			 Math.abs(startA,substring.minStartA)<=THRESHOLD || 
			 (orientation?Math.abs(startB,substring.minStartB)<=THRESHOLD:Math.abs(endB,substring.maxEndB)<=THRESHOLD)
		   ) previousStartA=startA;
		else previousStartA=-1;
		if ( isRightMaximal==1 ||
		     Math.abs(endA,substring.maxEndA)<=THRESHOLD ||
			 (orientation?Math.abs(endB,substring.maxEndB)<=THRESHOLD:Math.abs(startB,substring.minStartB)<=THRESHOLD)
		   ) previousEndA=endA;
		else previousEndA=-1;
		if ( (orientation?isLeftMaximal:isRightMaximal)==1 || 
		     Math.abs(startB,substring.minStartB)<=THRESHOLD ||
			 (orientation?Math.abs(startA,substring.minStartA)<=THRESHOLD:Math.abs(endA,substring.maxEndA)<=THRESHOLD)
		   ) previousStartB=startB;
		else previousStartB=-1;
		if ( (orientation?isRightMaximal:isLeftMaximal)==1 || 
		     Math.abs(endB,substring.maxEndB)<=THRESHOLD ||
			 (orientation?Math.abs(endA,substring.maxEndA)<=THRESHOLD:Math.abs(startA,substring.minStartA)<=THRESHOLD)
	       ) previousEndB=endB;
		else previousEndB=-1;
		lastShift=-1;
		for (i=1; i<=lastAlignment; i++) {
			alignment=ReadA.sortedAlignments[firstAlignment+alignments[i]];
			id=alignment.id;
			orientation=Alignments.alignments[id][2]==1;
			isLeftMaximal=alignment.isLeftMaximal;
			startA=Alignments.alignments[id][3];
			isRightMaximal=alignment.isRightMaximal;
			endA=Alignments.alignments[id][4];
			startB=Alignments.alignments[id][5];
			endB=Alignments.alignments[id][6];
			if ( isLeftMaximal==1 || 
				 Math.abs(startA,substring.minStartA)<=THRESHOLD || 
				 (orientation?Math.abs(startB,substring.minStartB)<=THRESHOLD:Math.abs(endB,substring.maxEndB)<=THRESHOLD)
			   ) {
				if (previousStartA!=-1) {
					shift=Math.abs(startA,previousStartA);
					lastShift++;
					ensureSpace_shifts(lastShift+1);
					shifts[lastShift].position=shift;
					shifts[lastShift].mass=1;
				}
				previousStartA=startA;
			}
			if ( isRightMaximal==1 || 
			     Math.abs(endA,substring.maxEndA)<=THRESHOLD ||
				 (orientation?Math.abs(endB,substring.maxEndB)<=THRESHOLD:Math.abs(startB,substring.minStartB)<=THRESHOLD)
			   ) {
				if (previousEndA!=-1) {
					shift=Math.abs(endA,previousEndA);
					lastShift++;
					ensureSpace_shifts(lastShift+1);
					shifts[lastShift].position=shift;
					shifts[lastShift].mass=1;
				}
				previousEndA=endA;
			}
			if ( (orientation?isLeftMaximal:isRightMaximal)==1 || 
			     Math.abs(startB,substring.minStartB)<=THRESHOLD ||
				 (orientation?Math.abs(startA,substring.minStartA)<=THRESHOLD:Math.abs(endA,substring.maxEndA)<=THRESHOLD)
			   ) {
				if (previousStartB!=-1) {
					shift=Math.abs(startB,previousStartB);
					lastShift++;
					ensureSpace_shifts(lastShift+1);
					shifts[lastShift].position=shift;
					shifts[lastShift].mass=1;
				}
				previousStartB=startB;
			}
			if ( (orientation?isRightMaximal:isLeftMaximal)==1 || 
			     Math.abs(endB,substring.maxEndB)<=THRESHOLD ||
				 (orientation?Math.abs(endA,substring.maxEndA)<=THRESHOLD:Math.abs(startA,substring.minStartA)<=THRESHOLD)
		       ) {
				if (previousEndB!=-1) {
					shift=Math.abs(endB,previousEndB);
					lastShift++;
					ensureSpace_shifts(lastShift+1);
					shifts[lastShift].position=shift;
					shifts[lastShift].mass=1;
				}
				previousEndB=endB;
			}
		}
		lastShift=Points.sortAndCompact(shifts,lastShift);
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("collectShifts> substring: "+substring);
	IO.printErr("collectShifts> collected shifts: (lastShift="+lastShift+")");
	for (int x=0; x<=lastShift; x++) System.err.println(shifts[x]);
	for (int x=1; x<=substring.lastShift; x++) {
		if (substring.shifts[x].position<substring.shifts[x-1].position) {
			System.err.println("ERROR: collectShifts() produced a non-sorted file.");
			System.exit(1);
		}
	}
}

		substring.cloneShifts(shifts,0,lastShift,MAX_DISTINCT_SHIFTS,Points.tmpPoints);
	}

	
	/**
	 * @return TRUE iff the readA part of a periodic substring corresponds to one unit of
	 * its period; the procedure analyzes just the longest-path alignments in 
	 * $alignments$.
	 */
	private static final boolean equalsOnePeriod(int firstAlignment, int minStartA, int maxEndA, int minStartB, int maxEndB, int threshold) {
		final int INTERSECTION_THRESHOLD = Alignments.minAlignmentLength>>1;  // Arbitrary
		int i, j;
		int id, startA, endA, startB, endB, nPrefix, nSuffix, nFull;
		int idPrime, startAPrime, endAPrime, startBPrime, endBPrime;
		
		// Checking alignment types
		nFull=0; nPrefix=0; nSuffix=0;
		for (i=0; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[firstAlignment+alignments[i]].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			startB=Alignments.alignments[id][5];
			endB=Alignments.alignments[id][6];
			if (Math.abs(startA,minStartA)<=threshold && Math.abs(endA,maxEndA)<=threshold) nFull++;
			else if (Math.abs(startA,minStartA)<=threshold && endA<maxEndA-threshold && Math.abs(endB,maxEndB)<=threshold) {
				nPrefix++;
				if (nPrefix>1) return false;
			}
			else if (startA>minStartA+threshold && Math.abs(endA,maxEndA)<=threshold && Math.abs(startB,minStartB)<=threshold) {
				nSuffix++;
				if (nSuffix>1) return false;
			}
			else if ( (startA>minStartA+threshold && ReadA.sortedAlignments[firstAlignment+alignments[i]].isLeftMaximal==1) || 
			          (endA<maxEndA-threshold && ReadA.sortedAlignments[firstAlignment+alignments[i]].isRightMaximal==1)
			        ) {
				return false;
			}
		}
		if (nFull==0) return false;
		
		// Checking whether full alignments overlap in readB.
		for (i=0; i<=lastAlignment; i++) {
			id=ReadA.sortedAlignments[firstAlignment+alignments[i]].id;
			startA=Alignments.alignments[id][3];
			endA=Alignments.alignments[id][4];
			if (Math.abs(startA,minStartA)>threshold || Math.abs(endA,maxEndA)>threshold) continue;
			startB=Alignments.alignments[id][5];
			endB=Alignments.alignments[id][6];
			for (j=i+1; j<=lastAlignment; j++) {
				idPrime=ReadA.sortedAlignments[firstAlignment+alignments[j]].id;
				startAPrime=Alignments.alignments[idPrime][3];
				endAPrime=Alignments.alignments[idPrime][4];
				if (Math.abs(startAPrime,minStartA)>threshold || Math.abs(endAPrime,maxEndA)>threshold) continue;
				startBPrime=Alignments.alignments[idPrime][5];
				endBPrime=Alignments.alignments[idPrime][6];
				if (Intervals.intersectionLength(startB,endB,startBPrime,endBPrime)>=INTERSECTION_THRESHOLD) return false;
			}
		}
		return true;
	}	


	/**
	 * Tries to estimate the length of a period from a collection of observed shifts in 
	 * $readA$ and $readB$, by estimating the density distribution of shifts using a 
	 * density estimation tree, and by inferring the period from the peaks in such
	 * distribution as follows. The distribution of all shifts should have a peak around
	 * zero, a peak around the shortest period, and a peak around every integer multiple
	 * of the shortest period. The procedure assumes that the distribution consists of an
	 * arbitrary subset of such peaks. It finds the second peak (since the first peak
	 * might be the zero peak), and it divides it by increasing integers, until all other
	 * peaks (including possibly the first) are a multiple of the result of the division,
	 * which is reported as the period. If there are less than two peaks in the 
	 * distribution of shifts, the procedure does not estimate the period.
	 * For robustness, only high peaks are considered.
	 * See procedure $estimatePeriodFromPeak$ for details.
	 *
	 * The procedure sets $periodTmp$ as follows:
	 *
     * 0: number of local maxima found (not necessarily high, including the one around 
	 * zero), or zero if the procedure cannot estimate the number of local maxima;
	 * 1: center of mass of the first local-maximum run;
     * 2: center of mass of the last local-maximum run;
	 * 3: position of the largest point in the first local-maximum run;
	 * 4: number of high local maxima.
	 *
	 * Remark: if there is just one peak in the distribution of shifts, the procedure uses 
	 * it to compute the period only if it is bigger than a minimum, and if it is a 
	 * biologically known period (see $Biology.isKnownPeriod$).
	 *
	 * Remark: for a long-period substring, one might think of estimating the period by 
	 * building a DET on maximal events inside the substring, finding peaks, and computing
	 * the avg. distance between consecutive peaks. This is not correct in general, since 
	 * there might be multiple periodic regions in the genome, with the same period but 
	 * different phase. When aligned to the current interval of readA, such regions could 
	 * create a possibly uniform density of maximal events, so even the existence of peaks
	 * is not obvious in general.
	 *
	 * Remark: $shifts$ is assumed to be sorted and compacted.
	 *
	 * @return an estimate of the length of a period (not necessarily the shortest
	 * period), or zero if such estimate fails.
	 */
	public static final int estimatePeriod(Point[] shifts, int lastShift, int minPeriod, int minIntervalLength, int minLocalMaxDistance) {
		final int BIN_LENGTH = minPeriod>>1;
		final double MASS_THRESHOLD = 0.1;  // Arbitrary
		final double MASS_RATIO_THRESHOLD = 5;  // Arbitrary
		final int MAX_WIDTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		int MIN_MASS_FOR_PERIOD = 10;  // Arbitrary
		boolean isHigh;
		int i, j;
		int shift, period, first, last, nLocalMaximumLeaves, largestMass, nHigh, lowerBound;
		int lastLeaf;
		double center, massRatio;
		double[] out = new double[2];
		Leaf[] lvs;
		if (lastShift==-1) return 0;
		if (lastShift==0) return (int)(shifts[0].position);


if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("estimatePeriod> OBSERVED SHIFTS:");
	for (int x=0; x<=lastShift; x++) IO.printErr(shifts[x].position+","+shifts[x].getMass()); 
}


		if (Points.areUniformlyDistributed(shifts,0,lastShift,true,BIN_LENGTH)) {
			if (shifts[lastShift].position-shifts[0].position<=minPeriod && Points.mass(shifts,0,lastShift)>=MIN_MASS_FOR_PERIOD) {
				// Considering this case as a single peak, even if it might be smaller
				// than $minPeriod$.
				periodTmp[0]=1;
				periodTmp[1]=Math.round(Points.getCenterOfMass(shifts,0,lastShift,false,-1,lastShift));
				periodTmp[2]=periodTmp[1];
				periodTmp[3]=(int)shifts[lastShift].position;
				periodTmp[4]=0;
				return (minPeriod==MIN_PERIOD_LONG||Biology.isKnownPeriod(periodTmp[1]))?periodTmp[1]:0;
			}
			periodTmp[0]=0;
			periodTmp[1]=-1;
			periodTmp[2]=-1;
			periodTmp[3]=-1;
			periodTmp[4]=-1;
			return 0;
		}
		if (Points.mass(shifts,0,lastShift)<=Points.FEW_POINTS) {
			Points.getRightmostInterval(shifts,lastShift,minPeriod,MIN_MASS_FOR_PERIOD,out);
			if (out[0]>=0) {
				// Considering this case as a single peak, even if it might be smaller
				// than $minPeriod$.
				periodTmp[0]=1;
				periodTmp[1]=(int)out[0];
				periodTmp[2]=periodTmp[1];
				periodTmp[3]=(int)out[1];
				periodTmp[4]=0;
				return (minPeriod==MIN_PERIOD_LONG||Biology.isKnownPeriod(periodTmp[1]))?periodTmp[1]:0;
			}
			periodTmp[0]=0;
			periodTmp[1]=-1;
			periodTmp[2]=-1;
			periodTmp[3]=-1;
			periodTmp[4]=-1;
			return 0;
		}

		largestMass=Points.largestMass(shifts,lastShift,-1);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("largestMass="+largestMass);
		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(shifts,0,lastShift,minIntervalLength,minLocalMaxDistance,true,largestMass,minPeriod,minPeriod>>1,true,true);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr(">>>> LEAVES:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr("firstPoint="+DensityEstimationTree.leaves[x].firstPoint+" lastPoint="+DensityEstimationTree.leaves[x].lastPoint+" isLocalMaximum="+DensityEstimationTree.leaves[x].isLocalMaximum+" marked="+DensityEstimationTree.leaves[x].marked+" isHigh="+DensityEstimationTree.leaves[x].isHigh); }
		
		if (nLocalMaximumLeaves==0) {
			periodTmp[0]=0;
			periodTmp[1]=-1;
			periodTmp[2]=-1;
			periodTmp[3]=-1;
			periodTmp[4]=-1;
			return 0;
		}
		else {
			lvs=DensityEstimationTree.leaves;
			DensityEstimationTree.leaves=tmpLeaves;
			tmpLeaves=lvs;
			lastLeaf=DensityEstimationTree.lastLeaf;
			Leaves.setHighLocalMax(tmpLeaves,lastLeaf,lastShift,true,false,shifts,Points.mass(shifts,0,lastShift),true,MAX_WIDTH);
			lvs=tmpLeaves;
			tmpLeaves=DensityEstimationTree.leaves;
			DensityEstimationTree.leaves=lvs;
			DensityEstimationTree.lastLeaf=lastLeaf;
			DensityEstimationTree.markRunsOfLocalMaximumLeaves(shifts);
			
if (IO.SHOW_STD_ERR_PRIME) IO.printErr(">>>> LEAVES (after markRuns):");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr("firstPoint="+DensityEstimationTree.leaves[x].firstPoint+" lastPoint="+DensityEstimationTree.leaves[x].lastPoint+" isLocalMaximum="+DensityEstimationTree.leaves[x].isLocalMaximum+" marked="+DensityEstimationTree.leaves[x].marked+" isHigh="+DensityEstimationTree.leaves[x].isHigh); }					
			
			// Estimating the period using high runs of local-maximum leaves
			periodTmp[3]=(int)DensityEstimationTree.separateRun(0,shifts,0,lastShift,false);
			if (DensityEstimationTree.nRuns==1) {
				// Using the only run if its center of mass is at least $minPeriod$
				lowerBound=Math.max((int)(Points.mass(shifts,0,lastShift)*MASS_THRESHOLD),2);
				MIN_MASS_FOR_PERIOD=Math.min(MIN_MASS_FOR_PERIOD,lowerBound);
				periodTmp[0]=1;
				periodTmp[1]=Math.round(DensityEstimationTree.getCenterOfMassOfRun(0,shifts,lastShift,-1));
				periodTmp[2]=periodTmp[1];
				periodTmp[4]=0;
				return periodTmp[1]>=minPeriod && 
					   DensityEstimationTree.getMassOfRun(0,shifts)>=MIN_MASS_FOR_PERIOD && 
					   (minPeriod==MIN_PERIOD_LONG||Biology.isKnownPeriod(periodTmp[1])) ? periodTmp[1] : 0;
			}
			periodTmp[0]=DensityEstimationTree.nRuns;
			periodTmp[1]=Math.round(DensityEstimationTree.getCenterOfMassOfRun(0,shifts,lastShift,-1));
			periodTmp[2]=Math.round(DensityEstimationTree.getCenterOfMassOfRun(DensityEstimationTree.nRuns-1,shifts,lastShift,-1));
			nHigh=Leaves.getNumberOHighRuns(DensityEstimationTree.leaves,DensityEstimationTree.lastLocalMaximum);
			periodTmp[4]=nHigh;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("nHigh="+nHigh);			
			if (nHigh==0) return 0;
			else if (nHigh==1) {
				// Using the only high run if its center of mass is at least $minPeriod$.
				lowerBound=Math.max((int)(Points.mass(shifts,0,lastShift)*MASS_THRESHOLD),2);
				MIN_MASS_FOR_PERIOD=Math.min(MIN_MASS_FOR_PERIOD,lowerBound);
				DensityEstimationTree.getCenterOfMassOfHighRun(0,true,shifts,out);
				period=(int)out[0];
				return period>=minPeriod && 
					   DensityEstimationTree.getMassOfHighRun(0,shifts)>=MIN_MASS_FOR_PERIOD && 
					   (minPeriod==MIN_PERIOD_LONG||Biology.isKnownPeriod(period)) ? period : 0;
			}
			else if (nHigh==2) {
				DensityEstimationTree.getCenterOfMassOfHighRun(0,true,shifts,out);
				period=(int)out[0];
				if (period<minPeriod) {
					// The first high run is around zero
					lowerBound=Math.max((int)(Points.mass(shifts,Leaves.getLastPointOfHighRun(DensityEstimationTree.leaves,DensityEstimationTree.lastLocalMaximum,0)+1,lastShift)*MASS_THRESHOLD),2);
					MIN_MASS_FOR_PERIOD=Math.min(MIN_MASS_FOR_PERIOD,lowerBound);
					DensityEstimationTree.getCenterOfMassOfHighRun(shifts[lastShift].position,false,shifts,out);
					period=(int)out[0];
					return period>=minPeriod && DensityEstimationTree.getMassOfHighRun(1,shifts)>=MIN_MASS_FOR_PERIOD ? period : 0;
				}
				else {
					// The first high run is not around zero
					lowerBound=Math.max((int)(Points.mass(shifts,0,lastShift)*MASS_THRESHOLD),2);
					MIN_MASS_FOR_PERIOD=Math.min(MIN_MASS_FOR_PERIOD,lowerBound);
					DensityEstimationTree.findHighRun(minPeriod,shifts,MIN_MASS_FOR_PERIOD,true);
					if (DensityEstimationTree.highRunTmp[0]==-1) return 0;
					return (int)Points.getCenterOfMass(shifts,DensityEstimationTree.leaves[DensityEstimationTree.highRunTmp[0]].firstPoint,DensityEstimationTree.leaves[DensityEstimationTree.highRunTmp[1]].lastPoint,false,-1,lastShift);
				}
			}
			else {
				DensityEstimationTree.getCenterOfMassOfHighRun(0,true,shifts,out);
				period=(int)out[0];
				if (period<minPeriod) {
					// The first high run is around zero
					lowerBound=Math.max((int)(Points.mass(shifts,Leaves.getLastPointOfHighRun(DensityEstimationTree.leaves,DensityEstimationTree.lastLocalMaximum,0)+1,lastShift)*MASS_THRESHOLD),2);
					MIN_MASS_FOR_PERIOD=Math.min(MIN_MASS_FOR_PERIOD,lowerBound);
					massRatio=DensityEstimationTree.findHighRun(minPeriod,shifts,MIN_MASS_FOR_PERIOD,false);
				}
				else {
					// The first high run is not around zero
					lowerBound=Math.max((int)(Points.mass(shifts,0,lastShift)*MASS_THRESHOLD),2);
					MIN_MASS_FOR_PERIOD=Math.min(MIN_MASS_FOR_PERIOD,lowerBound);
					massRatio=DensityEstimationTree.findHighRun(minPeriod,shifts,MIN_MASS_FOR_PERIOD,true);
				}	
				if (DensityEstimationTree.highRunTmp[0]==-1) return 0;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("(3) massRatio="+massRatio+" MASS_RATIO_THRESHOLD="+MASS_RATIO_THRESHOLD);
				if (massRatio>=MASS_RATIO_THRESHOLD) return (int)Points.getCenterOfMass(shifts,DensityEstimationTree.leaves[DensityEstimationTree.highRunTmp[0]].firstPoint,DensityEstimationTree.leaves[DensityEstimationTree.highRunTmp[1]].lastPoint,false,-1,lastShift);
				else {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("(3) estimating period with high run ["+DensityEstimationTree.highRunTmp[0]+".."+DensityEstimationTree.highRunTmp[1]+"] MIN_MASS_FOR_PERIOD="+MIN_MASS_FOR_PERIOD);
					return estimatePeriodFromPeak(DensityEstimationTree.highRunTmp[0],DensityEstimationTree.highRunTmp[1],shifts,lastShift,nHigh,minPeriod);
				}
			}
		}
	}


	/**
	 * Assume that a density estimation tree has been built on a set of observed shifts,
	 * and let $[firstLeaf..lastLeaf]$ be a high run of local-maximum leaves in the tree. 
	 * For every integer $x \in [firstLeaf.firstPoint.position..
	 * lastLeaf.lastPoint.position]$, we assign a score to integer $x'=Math.round(x/d)$, 
	 * for increasing integers $d$, as follows. For every other high run of local-maximum 
	 * leaves in the tree, we increment the score of $x'$ by one if the high run contains 
	 * an integer that is a multiple of $x'$. The procedure returns the largest 
	 * $x'>=minPeriod$ among those with maximum score.
	 *
	 * Remark: only periods that are close enough to an observed shift are considered.
	 * Only their multiples that are close enough to an observed shift contribute to their
	 * score.
	 *
	 * Remark: the procedure returns zero if the maximum score of an $x'$ is smaller than
	 * a fraction of the total number of high runs $nHigh$. In this case, the set of 
	 * shifts is not likely to be induced by an underlying period.
	 *
	 * Remark: the procedure assumes that procedure
	 * $DensityEstimationTree.markRunsOfLocalMaximumLeaves$ has already been executed.
	 */
	private static final int estimatePeriodFromPeak(int firstLeaf, int lastLeaf, Point[] shifts, int lastShift, int nHigh, int minPeriod) {
		final int DISTANCE_THRESHOLD = IO.quantum;  // Distance from an observed shift
		final int MAX_DIVISOR = 6;  // Arbitrary
		final double THRESHOLD = 0.5;
		boolean isHigh, found;
		int i, j, k, h, m, p;
		int firstLocalMaximum, divisor, period, bestPeriod, score, bestScore;
		int firstPoint, lastPoint, firstPosition, lastPosition;
		int firstPointPrime, lastPointPrime, firstPositionPrime, lastPositionPrime;
		int minPosition, maxPosition, minPositionPrime, maxPositionPrime;
		
		for (i=0; i<=lastShift; i++) {
			if (shifts[i].position!=0) break;
		}
		minPeriod=Math.max(minPeriod,(int)(shifts[i].position));
		bestScore=0; bestPeriod=-1;
		firstPoint=DensityEstimationTree.leaves[firstLeaf].firstPoint;
		lastPoint=DensityEstimationTree.leaves[lastLeaf].lastPoint;
		firstPosition=(int)shifts[firstPoint].position;
		lastPosition=(int)shifts[lastPoint].position;
		for (i=lastPoint; i>=firstPoint; i--) {
			p=(int)shifts[i].position;
			if (p<minPeriod) break;
			minPosition=Math.max(p-DISTANCE_THRESHOLD,firstPosition);
			minPosition=Math.max(minPosition,minPeriod);
			maxPosition=Math.min(p+DISTANCE_THRESHOLD,lastPosition);
			for (j=maxPosition; j>=minPosition; j--) {
				period=j; divisor=1;
				while (period>=minPeriod && divisor<=MAX_DIVISOR) {
					score=0; firstLocalMaximum=0; isHigh=false;
					for (k=0; k<=DensityEstimationTree.lastLocalMaximum; k++) {
						if (DensityEstimationTree.leaves[k].isHigh) isHigh=true;
						if (!DensityEstimationTree.leaves[k].marked) continue;
						if (!isHigh || (firstLocalMaximum==firstLeaf&&k==lastLeaf)) {
							firstLocalMaximum=k+1;
							isHigh=false;
							continue;
						}
						firstPointPrime=DensityEstimationTree.leaves[firstLocalMaximum].firstPoint;
						lastPointPrime=DensityEstimationTree.leaves[k].lastPoint;
						firstPositionPrime=(int)shifts[firstPointPrime].position;
						lastPositionPrime=(int)shifts[lastPointPrime].position;
						found=false;
						for (h=firstPointPrime; h<=lastPointPrime; h++) {
							p=(int)shifts[h].position;
							minPositionPrime=Math.max(p-DISTANCE_THRESHOLD,firstPositionPrime);
							maxPositionPrime=Math.min(p+DISTANCE_THRESHOLD,lastPositionPrime);
							for (m=minPositionPrime; m<=maxPositionPrime; m++) {
								if (m!=0 && m%period==0) {
									found=true;
									break;
								}
							}
							if (found) break;
						}
						if (found) score++;
						firstLocalMaximum=k+1;
						isHigh=false;
					}
					if (score>bestScore) {
						bestScore=score;
						bestPeriod=period;
					}
					else if (score==bestScore && period>bestPeriod) bestPeriod=period;
					period=Math.round(j,++divisor);
				}
			}
		}
		return ((nHigh<=3&&bestScore>=1) || (nHigh>3&&bestScore>=THRESHOLD*nHigh))?bestPeriod:0;
	}


	/**
	 * Detects periodic substrings that are almost equivalent to one another (see
	 * procedure $getConnectedComponent$), and for each such set elects a substring as
	 * representative, resets its values of $minStartA$, $maxEndA$, $minStartB$, $maxEndB$
	 * and $shifts$ to the union of all the substrings in the set, and marks all the 
	 * substrings in the set (excluding the representative) as not to be used in computing
	 * corrected coverage (see procedure $Factors.computeCoverages$).
     *
     * Then, the procedure marks some representative periodic substrings as not to be
     * used in the computation of corrected coverage. Such objects are those that arise
     * from a periodic substring in $readA$ being equal to its reverse-complement, and
     * from a periodic substring in a $readB$ being broken down into multiple substrings
     * by low-quality insertions.
	 *
	 * Remark: the procedure updates the set of shifts of each substring. 
	 * If $strings==longPeriodSubstrings$, it also updates the period.
	 */
	private static final void correctCoverage(PeriodicSubstring[] strings, int lastString) {
		boolean found, orientation;
		int i, j, firstJForNextI;
		int readB, readBLength, nMerges, period;
		double intersectionLength;
		PeriodicSubstring tmp, from, to, destination;

		// Ensuring the necessary order in $strings$
		if (strings==substrings && PeriodicSubstring.order!=PeriodicSubstring.READB_ORIENTATION_MINSTARTA) {
			PeriodicSubstring.order=PeriodicSubstring.READB_ORIENTATION_MINSTARTA;
			if (lastString>0) Arrays.sort(strings,0,lastString+1);
		}
		if (strings==longPeriodSubstrings) {
			if (PeriodicSubstring.order_longPeriod!=PeriodicSubstring.READB_ORIENTATION_MINSTARTA) {
				PeriodicSubstring.order_longPeriod=PeriodicSubstring.READB_ORIENTATION_MINSTARTA;
				if (lastString>0) Arrays.sort(strings,0,lastString+1);
			}
		}
		for (i=0; i<=lastString; i++) {
			strings[i].shouldBeCountedForCoverage=true;
			strings[i].threadStart=strings[i];
			strings[i].equalsReverseComplement=false;
		}

		// Detecting almost-equivalent periodic substrings
		for (i=0; i<=lastString; i++) strings[i].representative=null;
		for (i=0; i<=lastString; i++) {
			if (strings[i].representative!=null) continue;
			strings[i].representative=strings[i];
			nMerges=getConnectedComponent(i,strings,lastString,strings==substrings);
			if (nMerges>0 && strings==longPeriodSubstrings) {
				period=estimateLongPeriod(strings[i].shifts,strings[i].lastShift,null);
				strings[i].period=Math.max(strings[i].period,period);
			}
		}		

		// Discarding non-representative substrings, and updating
		// $impliedByPeriodicSubstring$ pointers in $Alignment$ objects.
		j=-1;
		for (i=0; i<=lastString; i++) {
			if (strings[i].representative!=strings[i]) continue;
			j++;
			if (j!=i) {
				tmp=strings[j];
				strings[j]=strings[i];
				strings[i]=tmp;
			}
		}
		lastString=j;
		if (strings==substrings) {
			lastSubstring=lastString;
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null || ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hasLongPeriod) continue;
				ReadA.sortedAlignments[i].impliedByPeriodicSubstring=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.representative;
			}
		}
		else if (strings==longPeriodSubstrings) {
			lastLongPeriodSubstring=lastString;
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null || !ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hasLongPeriod) continue;
				ReadA.sortedAlignments[i].impliedByPeriodicSubstring=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.representative;
			}
		}

		// Detecting duplicate candidate substrings that arise when a periodic substring
		// equals its reverse complement.
		for (i=0; i<lastString; i++) {
			if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Using periodic substring "+i+" to kill the others");
			readBLength=Reads.getReadLength(strings[i].readB);
			if (IO.SHOW_STD_ERR_PRIME) IO.printErr("readBLength="+readBLength);
			for (j=i+1; j<=lastString; j++) {
				if (IO.SHOW_STD_ERR_PRIME) IO.printErr("   trying to kill periodic substring "+j);
				if (strings[j].readB!=strings[i].readB) break;
				if (strings[j].orientation==strings[i].orientation) continue;
				if (strings[j].minStartA>strings[i].maxEndA) break;
				if (!strings[j].shouldBeCountedForCoverage) continue;
				if (Intervals.jaccardSimilarity(strings[i].minStartA,strings[i].maxEndA,strings[j].minStartA,strings[j].maxEndA)<Intervals.jaccardThreshold) continue;
				if (IO.SHOW_STD_ERR_PRIME) IO.printErr("   high similarity of readA intervals");
				intersectionLength=Intervals.intersectionLength(strings[i].minStartB,strings[i].maxEndB,readBLength-strings[j].maxEndB-1,readBLength-strings[j].minStartB-1);  // It is always the case that $strings[i].orientation=true$ and $strings[j].orientation=false$.
				if (intersectionLength==0) continue;
				if (Intervals.jaccardSimilarity(strings[i].minStartB,strings[i].maxEndB,readBLength-strings[j].maxEndB-1,readBLength-strings[j].minStartB-1)<Intervals.jaccardThreshold) continue;  // It is always the case that $strings[i].orientation=true$ and $strings[j].orientation=false$.
				if (IO.SHOW_STD_ERR_PRIME) IO.printErr("   high similarity of readB intervals --> killed");
				strings[j].shouldBeCountedForCoverage=false;
				strings[i].equalsReverseComplement=true;
				strings[j].equalsReverseComplement=true;
				strings[i].equalsOnePeriod&=strings[j].equalsOnePeriod;
			}
		}

		// Detecting duplicate candidate substrings that arise when a single periodic
		// substring in $readB$ is broken into multiple periodic substrings by a set of
		// random insertions.
		readB=-1; orientation=false; firstJForNextI=-1;
		for (i=0; i<=lastString; i++) {  // We use also substrings that should not be counted for coverage
			if (strings[i].readB!=readB) {
				readB=strings[i].readB;
				orientation=strings[i].orientation;
				firstJForNextI=-1;
			}
			else if (strings[i].orientation!=orientation) {
				orientation=strings[i].orientation;
				firstJForNextI=-1;
			}
			if (firstJForNextI==-1) {
				j=i+1;
				if (i<lastString && strings[i+1].readB==strings[i].readB && strings[i+1].orientation==strings[i].orientation && strings[i].maxEndA>=strings[i+1].minStartA) firstJForNextI=i;
				else firstJForNextI=-1;
			}
			else {
				j=firstJForNextI;
				firstJForNextI=-1;
			}
			while (j<=lastString) {
				if (strings[j].readB!=strings[i].readB || strings[j].orientation!=strings[i].orientation) break;
				if (strings[j].maxEndA<strings[i].minStartA) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<lastString && strings[i+1].readB==strings[i].readB && strings[i+1].orientation==strings[i].orientation && j<i+1 && strings[j].maxEndA>=strings[i+1].minStartA) firstJForNextI=j;
				if (j==i) {
					j++;
					continue;
				}
				intersectionLength=Intervals.intersectionLength(strings[i].minStartA,strings[i].maxEndA,strings[j].minStartA,strings[j].maxEndA);
				if (intersectionLength==0) break;
				if (Intervals.jaccardSimilarity(strings[i].minStartA,strings[i].maxEndA,strings[j].minStartA,strings[j].maxEndA)<Intervals.jaccardThreshold) {
					j++;
					continue;
				}
				if (strings[j].minStartB<=strings[i].maxEndB) {
					j++;
					continue;
				}
				if (Reads.isRandomInsertion(strings[i].readB,strings[i].maxEndB+1,strings[j].minStartB-1,strings[i].orientation)) {
					strings[j].shouldBeCountedForCoverage=false;
					strings[j].threadStart=strings[i];
					if (strings[j].equalsReverseComplement) strings[i].equalsReverseComplement=true;
					strings[i].equalsOnePeriod&=strings[j].equalsOnePeriod;
				}
				j++;
			}
		}

		// Computing the closure of $threadStart$ links
		for (i=0; i<=lastString; i++) {
			from=strings[i];
			to=strings[i].threadStart;
			while (to!=from) {
				to=to.threadStart;
				from=to;
			}
			destination=to;
			from=strings[i];
			to=strings[i].threadStart;
			while (to!=from) {
				from.threadStart=destination;
				from=to;
				to=to.threadStart;
			}
		}
	}


	/**
	 * Let $([i..j],[x..y])$ and $([i'..j'],[x'..y'])$ be two periodic substrings in
	 * $substrings$, where the intervals indicate the minimum and the maximum position of
	 * the periodic substring in $readA$ and in a specific $readB$ and orientation. 
	 * If $equivalenceType=true$, we say that the two periodic substrings are equivalent 
	 * if either: (1) $[i..j]$ is almost entirely in $[i'..j']$ (or vice versa), and 
	 * $[x..y]$ is almost entirely in $[x'..y']$ (or viceversa); (2) the Jaccard 
	 * similarity of $[i..j]$ and $[i'..j']$ is high, and the Jaccard similarity of 
	 * $[x..y]$ and $[x'..y']$ is high (but the order of the intervals migh differ in 
	 * readA and readB). If $equivalenceType=false$, equivalence requires the intervals in
	 * both reads to be approximately identical.
	 *
	 * Let $G$ be the graph whose vertices are the elements of $substrings$, and whose
	 * edges are equivalence relations. The procedure finds the connected component of
	 * $fromSubstring$ in $G$. Every periodic substring in such component, except
	 * $fromSubstring$, is assigned the same value of $representative$ (which equals
	 * $substrings[fromSubstring]$), and it is marked as not to be counted for computing
	 * the periodic coverage of factors. The intervals of $fromSubstring$ in $readA$ and
	 * $readB$ are set to the union of all the intervals in its closure, and the $shifts$
	 * arrays of all substrings in the same component are merged into the one of
	 * $substrings[fromSubstring]$.
	 *
	 * Remark: the procedure assumes $substrings$ to be sorted by $readB$, $orientation$,
	 * $minStartA$.
	 *
	 * Remark: the procedure assumes that all substrings in $substrings[0..fromSubstring
	 * -1]$ have already been assigned a connected component different from the one of
	 * $fromSubstring$.
	 *
	 * @return the number of substrings in the connected component of $fromSubstring$,
	 * excluding $fromSubstring$.
 	 */
	private static final int getConnectedComponent(int fromSubstring, PeriodicSubstring[] substrings, int lastSubstring, boolean equivalenceType) {
		final double JACCARD_THRESHOLD = 0.25;  // Arbitrary
		int i, j, top, count;
		Point[] tmp;

		count=0;
		top=0; stack[0]=fromSubstring;
		lastShift=Points.simpleClone(substrings[fromSubstring].shifts,substrings[fromSubstring].lastShift,shifts);
		while (top>=0) {
			i=stack[top];
			top--;
			for (j=fromSubstring+1; j<=lastSubstring; j++) {
				if (j==i || substrings[j].representative!=null) continue;
				if ( substrings[j].readB!=substrings[i].readB ||
					 (substrings[j].orientation^substrings[i].orientation) ||
					 substrings[j].minStartA>substrings[i].maxEndA 
				   ) break;
				if ( ( equivalenceType && 
					   !( ( ( Intervals.isApproximatelyContained(substrings[i].minStartA,substrings[i].maxEndA,substrings[j].minStartA,substrings[j].maxEndA) ||
					          Intervals.isApproximatelyContained(substrings[j].minStartA,substrings[j].maxEndA,substrings[i].minStartA,substrings[i].maxEndA)
  				            ) &&
  					        ( Intervals.isApproximatelyContained(substrings[i].minStartB,substrings[i].maxEndB,substrings[j].minStartB,substrings[j].maxEndB) ||
					          Intervals.isApproximatelyContained(substrings[j].minStartB,substrings[j].maxEndB,substrings[i].minStartB,substrings[i].maxEndB)
					        )
						  ) ||
						  ( Intervals.jaccardSimilarity(substrings[i].minStartA,substrings[i].maxEndA,substrings[j].minStartA,substrings[j].maxEndA)>=JACCARD_THRESHOLD &&
							Intervals.jaccardSimilarity(substrings[i].minStartB,substrings[i].maxEndB,substrings[j].minStartB,substrings[j].maxEndB)>=JACCARD_THRESHOLD
						  )
					    )
					 ) ||
				     ( !equivalenceType && 
  					   !( Intervals.areApproximatelyIdentical(substrings[i].minStartA,substrings[i].maxEndA,substrings[j].minStartA,substrings[j].maxEndA) &&
    					  Intervals.areApproximatelyIdentical(substrings[i].minStartB,substrings[i].maxEndB,substrings[j].minStartB,substrings[j].maxEndB)
  					    )
  					 )
				   ) continue;
				count++;
				substrings[j].representative=substrings[i].representative;
				substrings[j].representative.equalsOnePeriod&=substrings[j].equalsOnePeriod;
				substrings[j].shouldBeCountedForCoverage=false;
				stack[++top]=j;
				if (substrings[j].minStartA<substrings[fromSubstring].minStartA) substrings[fromSubstring].minStartA=substrings[j].minStartA;
				if (substrings[j].maxEndA>substrings[fromSubstring].maxEndA) substrings[fromSubstring].maxEndA=substrings[j].maxEndA;
				if (substrings[j].minStartB<substrings[fromSubstring].minStartB) substrings[fromSubstring].minStartB=substrings[j].minStartB;
				if (substrings[j].maxEndB>substrings[fromSubstring].maxEndB) substrings[fromSubstring].maxEndB=substrings[j].maxEndB;
				if (substrings[fromSubstring].orientation) {
					substrings[fromSubstring].minStartBForward=substrings[fromSubstring].minStartB;
					substrings[fromSubstring].maxEndBForward=substrings[fromSubstring].maxEndB;
				}
				else {
					substrings[fromSubstring].minStartBForward=Reads.getReadLength(substrings[fromSubstring].readB)-substrings[fromSubstring].maxEndB-1;
					substrings[fromSubstring].maxEndBForward=Reads.getReadLength(substrings[fromSubstring].readB)-substrings[fromSubstring].minStartB-1;
				}
				lastShift=Points.merge(shifts,lastShift,substrings[j].shifts,substrings[j].lastShift,shiftsPrime);
				tmp=shifts;
				shifts=shiftsPrime;
				shiftsPrime=tmp;	
			}
		}
		substrings[fromSubstring].cloneShifts(shifts,0,lastShift,MAX_DISTINCT_SHIFTS,Points.tmpPoints);
		return count;
	}


	/**
	 * Merges elements of $substrings[0..lastSubstringForFactoring]$ (which is assumed not
	 * to be empty) with similar readA intervals; then, discards intervals that are 
	 * contained in a short-period range but that do not align to its peaks of maximal 
	 * events: see $filterIntervalsWithPeaks()$ for details.
	 *
	 * Remark: the resulting intervals can be nested or straddling.
	 *
	 * Remark: the procedure does not alter $substrings$.
	 */
	private static final void getIntervals() {
		int i;

		// Detecting intervals
		for (i=0; i<=lastSubstring; i++) substrings[i].mergedToInterval=null;
		lastInterval=-1;
		for (i=0; i<=lastSubstringForFactoring; i++) {
			if (substrings[i].mergedToInterval!=null) continue;
			lastInterval=getInterval(substrings,i,lastSubstringForFactoring,intervals,lastInterval);
		}
		if (lastInterval==-1) return;
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
		if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		estimatePeriodOfShortPeriodIntervals();
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("getIntervals> Periodic substring intervals before filtering with peaks:");
	{ for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x].hashCode()+" :: "+intervals[x]); }
}		
		
		// Filtering
		splitShortPeriodIntervals();
		if (lastInterval>0) {
			lastInterval=filterIntervalsWithPeaks(intervals,lastInterval,substrings,lastSubstringForFactoring,false,true,true);
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("getIntervals> Periodic substring intervals after filtering with peaks 1:");
	{ for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x].hashCode()+" :: "+intervals[x]); }
}
			
			lastInterval=filterIntervalsWithPeaks(intervals,lastInterval,substrings,lastSubstringForFactoring,false,true,true);
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("getIntervals> Periodic substring intervals after filtering with peaks 2:");
	{ for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x].hashCode()+" :: "+intervals[x]); }
}
			
		}
	}

	
	/**
	 * Computes the connected component of relation $Intervals.areApproximatelyIdentical$
	 * to which the readA part of $substrings[fromSubstring]$ belongs, and adds a new 
	 * interval to $intervals$ corresponding to the merge of all substrings in the 
	 * connected component.
	 *
	 * Remark: the procedure incurs the risk that the boundaries of the connected
	 * component become very different from the boundaries of $substrings[fromSubstring]$.
	 *
	 * Remark: when $intervals==longPeriodIntervals$: (1) the procedure might create both 
	 * long-period and short-period intervals inside $intervals$; (2) at the end of the 
	 * procedure, a long-period substring might be assigned to a short-period interval.
	 *
	 * Remark: using $sum*A$ and $n*A$ for merging does not give good results in practice,
	 * even when these quantities are correctly updated during the whole pipeline.
	 *
	 * @param lastSubstring last substring to be used;
	 * @return the updated value of $lastInterval$.
	 */
	private static final int getInterval(PeriodicSubstring[] substrings, int fromSubstring, int lastSubstring, PeriodicSubstringInterval[] intervals, int lastInterval) {
		int i, j;
		int top, min, max, sumFirst, sumLast, nMergedIntervals, longestAlignment;
		Point[] tmp;

		lastInterval++;	
		if (intervals[lastInterval]==null) {
			intervals[lastInterval] = new PeriodicSubstringInterval();
			intervals[lastInterval].allocateMemory(MAX_DISTINCT_SHIFTS);
		}
		substrings[fromSubstring].mergedToInterval=intervals[lastInterval];
		intervals[lastInterval].equalsOnePeriod=substrings[fromSubstring].equalsOnePeriod;
		intervals[lastInterval].hasLongPeriod=substrings[fromSubstring].hasLongPeriod;
		top=0; stack[top]=fromSubstring;
		min=substrings[fromSubstring].minStartA;
		max=substrings[fromSubstring].maxEndA;
		longestAlignment=substrings[fromSubstring].longestAlignment;
		sumFirst=min; sumLast=max;
		nMergedIntervals=1;
		lastShift=Points.simpleClone(substrings[fromSubstring].shifts,substrings[fromSubstring].lastShift,shifts);
		while (top>=0) {
			i=stack[top--];
			for (j=fromSubstring+1; j<=lastSubstring; j++) {
				if (substrings[j].mergedToInterval!=null) continue;
				if (Intervals.areApproximatelyIdentical(substrings[j].minStartA,substrings[j].maxEndA,substrings[i].minStartA,substrings[i].maxEndA)) {
					if (substrings[j].minStartA<min) min=substrings[j].minStartA;
					if (substrings[j].maxEndA>max) max=substrings[j].maxEndA;
					sumFirst+=substrings[j].minStartA;
					sumLast+=substrings[j].maxEndA;
					nMergedIntervals++;
					substrings[j].mergedToInterval=intervals[lastInterval];
					intervals[lastInterval].equalsOnePeriod&=substrings[j].equalsOnePeriod;
					intervals[lastInterval].hasLongPeriod&=substrings[j].hasLongPeriod;
					longestAlignment=Math.max(longestAlignment,substrings[j].longestAlignment);
					lastShift=Points.merge(shifts,lastShift,substrings[j].shifts,substrings[j].lastShift,shiftsPrime);	
					tmp=shifts;
					shifts=shiftsPrime;
					shiftsPrime=tmp;
					stack[++top]=j;
				}
			}
		}
		intervals[lastInterval].id=0;
		intervals[lastInterval].firstPosition=sumFirst/nMergedIntervals;
		intervals[lastInterval].lastPosition=sumLast/nMergedIntervals;
		intervals[lastInterval].nMergedIntervals=nMergedIntervals;
		intervals[lastInterval].cloneShifts(shifts,0,lastShift,MAX_DISTINCT_SHIFTS,Points.tmpPoints);
		intervals[lastInterval].discarded=false;
		intervals[lastInterval].isContained=0;
		intervals[lastInterval].longestAlignment=longestAlignment;
		return lastInterval;
	}

	
	/**
	 * Removes periodic substrings that should contribute neither to $Events.events$ nor 
	 * to $intervals$, i.e. which are not a proof for splitting readA. The procedure 
	 * removes all periodic substrings that straddle two adjacent periodic substrings 
	 * $XY$, such that neither $X$ nor $Y$ themselves straddle two adjacent periodic 
	 * substrings: a straddling periodic substring just proves that the genome contains 
	 * substring $Y'X'$ where $X'$ (respectively, $Y'$) is a version of $X$ (respectively,
	 * $Y$) with possibly different length.
	 *
	 * Periodic substrings are not actually removed, but just positioned after 
	 * $lastSubstringForFactoring$ in $substrings$. The procedure sorts $substrings[0..
	 * lastSubstringForFactoring]$ by $minStartA$ and $maxEndA$. However, interval
	 * $[lastSubstringForFactoring+1..lastSubstring]$ is not sorted.
	 *
	 * Remark: a substring that straddles $XY$ proves that the periods of $X$ and $Y$ are 
	 * similar. A periodic substring cannot straddle just one other periodic substring, 
	 * e.g. it cannot contain a periodic suffix and a nonperiodic prefix.
	 * 
	 * Remark: just removing all straddling periodic substrings might lose information.
	 * For example, let $ABC$ be the concatenation of 3 strings, and let $A$, $AB$, $BC$,
	 * $C$ be periodic substrings. If we remove $AB$ and $BC$, the information that $B$ is 
	 * periodic is lost. If $B$ is short enough, there might not be other alignments that 
	 * prove $B$ to be periodic.
	 */
	private static final void removeSubstringsForFactoring() {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int EVENT_THRESHOLD = IO.quantum;
		final int ADJACENCY_THRESHOLD = IO.quantum;
		final double INTERNAL_STRADDLING_THRESHOLD = IO.quantum<<1;  // Arbitrary
		final double EXTERNAL_STRADDLING_THRESHOLD = IO.quantum<<1;  // Arbitrary
		int i, j, k;
		int first;
		PeriodicSubstring tmp;
		
		lastSubstringForFactoring=lastSubstring;
		if (lastSubstring==0) return;
		
		// Ensuring the necessary order in $substrings$
		if (PeriodicSubstring.order!=PeriodicSubstring.MINSTARTA_MAXENDA) {
			PeriodicSubstring.order=PeriodicSubstring.MINSTARTA_MAXENDA;
			if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		}

		// Marking periodic substrings that straddle two other periodic susbtrings
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].isStraddling=false;
			first=i;
			while (first<=lastSubstring && substrings[first].minStartA<=substrings[i].minStartA+IDENTITY_THRESHOLD) first++;
			first--;
			if (first==i) first=i-1;
			for (j=first; j>=0; j--) {
				if (j==i) continue;				
				if ( Math.abs(substrings[j].minStartA,substrings[i].minStartA)<=IDENTITY_THRESHOLD && 
				 	 substrings[j].maxEndA<substrings[i].maxEndA-IDENTITY_THRESHOLD ) {
					for (k=j+1; k<=lastSubstring; k++) {
						if (k==i) continue;
						if ( substrings[k].minStartA>substrings[i].minStartA+IDENTITY_THRESHOLD &&
							 substrings[k].minStartA<=substrings[i].maxEndA-INTERNAL_STRADDLING_THRESHOLD &&
							 Math.abs(substrings[k].minStartA,substrings[j].maxEndA)<=ADJACENCY_THRESHOLD &&
							 substrings[k].maxEndA-substrings[i].maxEndA>=EXTERNAL_STRADDLING_THRESHOLD
						   ) {
							substrings[i].isStraddling=true;
							break;
						}
					}
				}
				else if ( substrings[j].minStartA<=substrings[i].minStartA-EXTERNAL_STRADDLING_THRESHOLD &&
					      substrings[j].maxEndA-substrings[i].minStartA>=INTERNAL_STRADDLING_THRESHOLD &&
						  substrings[j].maxEndA<substrings[i].maxEndA-IDENTITY_THRESHOLD
				        ) {
					for (k=j+1; k<=lastSubstring; k++) {
						if (k==i) continue;
						if ( substrings[k].minStartA>substrings[i].minStartA+IDENTITY_THRESHOLD &&
							 Math.abs(substrings[k].minStartA,substrings[j].maxEndA)<=ADJACENCY_THRESHOLD &&
							 ( Math.abs(substrings[k].maxEndA,substrings[i].maxEndA)<=IDENTITY_THRESHOLD ||
							   ( substrings[k].minStartA<=substrings[i].maxEndA-INTERNAL_STRADDLING_THRESHOLD &&
								 substrings[k].maxEndA-substrings[i].maxEndA>=EXTERNAL_STRADDLING_THRESHOLD
							   )
							 )
						   ) {
							substrings[i].isStraddling=true;
							break;
						}
					}
				}
				if (substrings[i].isStraddling) break;
			}
		}
		
		// Discarding just periodic substrings that straddle two nonstraddling periodic
		// susbtrings
		for (i=0; i<=lastSubstring; i++) {
			substrings[i].shouldBeUsedForFactoring=true;
			if (!substrings[i].isStraddling) continue;
			first=i;
			while (first<=lastSubstring && substrings[first].minStartA<=substrings[i].minStartA+IDENTITY_THRESHOLD) first++;
			first--;
			if (first==i) first=i-1;
			for (j=first; j>=0; j--) {
				if (j==i || substrings[j].isStraddling) continue;
				if ( Math.abs(substrings[j].minStartA,substrings[i].minStartA)<=IDENTITY_THRESHOLD && 
				 	 substrings[j].maxEndA<substrings[i].maxEndA-IDENTITY_THRESHOLD ) {
					for (k=j+1; k<=lastSubstring; k++) {
						if (k==i || substrings[k].isStraddling) continue;
						if ( substrings[k].minStartA>substrings[i].minStartA+IDENTITY_THRESHOLD &&
							 substrings[k].minStartA<=substrings[i].maxEndA-INTERNAL_STRADDLING_THRESHOLD &&
							 Math.abs(substrings[k].minStartA,substrings[j].maxEndA)<=ADJACENCY_THRESHOLD &&
							 substrings[k].maxEndA-substrings[i].maxEndA>=EXTERNAL_STRADDLING_THRESHOLD
						   ) {
							substrings[i].shouldBeUsedForFactoring=false;
							break;
						}
					}
				}
				else if ( substrings[j].minStartA<=substrings[i].minStartA-EXTERNAL_STRADDLING_THRESHOLD &&
					      substrings[j].maxEndA-substrings[i].minStartA>=INTERNAL_STRADDLING_THRESHOLD &&
						  substrings[j].maxEndA<substrings[i].maxEndA-IDENTITY_THRESHOLD
				        ) {
					for (k=j+1; k<=lastSubstring; k++) {
						if (k==i || substrings[k].isStraddling) continue;
						if ( substrings[k].minStartA>substrings[i].minStartA+IDENTITY_THRESHOLD &&
							 Math.abs(substrings[k].minStartA,substrings[j].maxEndA)<=ADJACENCY_THRESHOLD &&
							 ( Math.abs(substrings[k].maxEndA,substrings[i].maxEndA)<=IDENTITY_THRESHOLD ||
							   ( substrings[k].minStartA<=substrings[i].maxEndA-INTERNAL_STRADDLING_THRESHOLD &&
								 substrings[k].maxEndA-substrings[i].maxEndA>=EXTERNAL_STRADDLING_THRESHOLD
							   )
							 )
						   ) {
							substrings[i].shouldBeUsedForFactoring=false;
							break;
						}
					}
				}
				if (!substrings[i].shouldBeUsedForFactoring) break;
			}
		}		
		k=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].shouldBeUsedForFactoring) {
				k++;
				if (k!=i) {
					tmp=substrings[k];
					substrings[k]=substrings[i];
					substrings[i]=tmp;
				}
			}
		}
		lastSubstringForFactoring=k;
	}
	
	
	/**
	 * Creates events from maximal alignments that are close enough to either a peak of 
	 * maximal events inside a maximal short-period range, or the beginning/end of a long-
	 * period interval (the beginning/end of short-period ranges should be contained in
	 * the peaks).
	 *
	 * @param computePeaks recomputes the peaks of maximal events inside maximal short-
	 * period ranges (otherwise, the peaks are assumed to have already been computed).
	 */
	private static final void addEvents(boolean computePeaks) {
		final int EVENT_THRESHOLD = IO.quantum;
		boolean found, previousSortByID;
		int i, j, k;
		int startA, endA;
		Alignment alignment;
		PeriodicSubstringInterval tmpInterval = new PeriodicSubstringInterval();
		
		if (computePeaks && lastInterval>=0) splitShortPeriodIntervals();
		
		// Creating events using left-maximal alignments
		if (PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (lastLongPeriodInterval>0) Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			alignment=ReadA.sortedAlignments[i];
			if (alignment.isLeftMaximal!=1) continue;
			startA=alignment.startA();
			found=false;
			if (lastLeftSplit>=0) {
				j=position2split(startA,leftSplits,lastLeftSplit,EVENT_THRESHOLD);
				found=(j>=0)&&(Math.abs(startA,leftSplits[3*j+1])<=EVENT_THRESHOLD);
			}
			if (!found && lastLongPeriodInterval>=0) {
				tmpInterval.firstPosition=startA;
				previousSortByID=PeriodicSubstringInterval.sortByID;
				PeriodicSubstringInterval.sortByID=false;
				j=Arrays.binarySearch(longPeriodIntervals,0,lastLongPeriodInterval+1,tmpInterval);
				PeriodicSubstringInterval.sortByID=previousSortByID;
				if (j>=0) found=true;
				else {
					j=-j-1;
					found = (j<=lastLongPeriodInterval && longPeriodIntervals[j].firstPosition<=startA+EVENT_THRESHOLD) ||
						    (j>0 && longPeriodIntervals[j-1].firstPosition>=startA-EVENT_THRESHOLD);
				}
			}
			if (!found) continue;
			Events.lastEvent++;
			Events.ensureSpace_events(Events.lastEvent+1);
			Events.events[Events.lastEvent].clear();
			Events.events[Events.lastEvent].position=startA;
			Events.events[Events.lastEvent].nOpen=1;
			ReadA.sortedAlignments[i].startAdded=true;
		}
		
		// Creating events using right-maximal alignments
		if (PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.LASTPOSITION) {
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.LASTPOSITION;
			if (lastLongPeriodInterval>0) Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			alignment=ReadA.sortedAlignments[i];
			if (alignment.isRightMaximal!=1) continue;
			endA=alignment.endA();
			found=false;
			if (lastRightSplit>=0) {
				j=position2split(endA,rightSplits,lastRightSplit,EVENT_THRESHOLD);
				found=(j>=0)&&(Math.abs(endA,rightSplits[3*j+1])<=EVENT_THRESHOLD);
			}
			if (!found && lastLongPeriodInterval>=0) {
				tmpInterval.lastPosition=endA;
				previousSortByID=PeriodicSubstringInterval.sortByID;
				PeriodicSubstringInterval.sortByID=false;
				j=Arrays.binarySearch(longPeriodIntervals,0,lastLongPeriodInterval+1,tmpInterval);
				PeriodicSubstringInterval.sortByID=previousSortByID;
				if (j>=0) found=true;
				else {
					j=-j-1;
					found = (j<=lastLongPeriodInterval && longPeriodIntervals[j].lastPosition<=endA+EVENT_THRESHOLD) ||
						    (j>0 && longPeriodIntervals[j-1].lastPosition>=endA-EVENT_THRESHOLD);
				}
			}
			if (!found) continue;
			Events.lastEvent++;
			Events.ensureSpace_events(Events.lastEvent+1);
			Events.events[Events.lastEvent].clear();
			Events.events[Events.lastEvent].position=endA;
			Events.events[Events.lastEvent].nClosed=1;
			ReadA.sortedAlignments[i].endAdded=true;
		}
		
		// Resetting the order of $longPeriodIntervals$.
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
		if (lastLongPeriodInterval>0) Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
	}


	/**
	 * Sets $isPrefixOfPeriodicSubstring$, $isSuffixOfPeriodicSubstring$, 
	 * $inPeriodicSubstring$ for every alignment in $ReadA.sortedAlignments$, using only 
	 * substrings in $substrings$. The $inPeriodicSubstring$ field is set using maximal
	 * ranges of straddling or adjacent short-period intervals, and it is assigned also
	 * to alignments that are not fully contained in a range, but that have a large 
	 * intersection with a range.
	 *
	 * Remark: the procedure assumes $substrings$ to be sorted by $minStartA$, and 
	 * $intervals$ to be sorted by $firstPosition$.
	 */
	private static final void markPrefixSuffix(PeriodicSubstring[] substrings, int lastSubstring, PeriodicSubstringInterval[] intervals, int lastInterval) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final double MIN_INTERSECTION_RATIO = 0.8;  // Arbitrary
		int i, j;
		int firstJForNextI, startA, endA, rangeFirst, rangeLast, firstAlignment, length;
		
		// Ensuring the necessary order in $ReadA.sortedAlignments$.
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		
		// Setting $is*OfPeriodicSubstring$
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			ReadA.sortedAlignments[i].isPrefixOfPeriodicSubstring=false;
			ReadA.sortedAlignments[i].isSuffixOfPeriodicSubstring=false;
			ReadA.sortedAlignments[i].inPeriodicSubstring=false;
		}
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastSubstring) {
			if (j>ReadA.lastSortedAlignment || Alignments.alignments[ReadA.sortedAlignments[j].id][3]>substrings[i].maxEndA) {
				i++;
				if (firstJForNextI>=0) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			endA=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
			if (firstJForNextI==-1 && i<lastSubstring && endA>=substrings[i+1].minStartA) firstJForNextI=j;
			startA=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
			if (startA<substrings[i].minStartA-DISTANCE_THRESHOLD) {
				j++;
				continue;
			}
			if (startA<=substrings[i].minStartA+DISTANCE_THRESHOLD && endA<substrings[i].maxEndA) ReadA.sortedAlignments[j].isPrefixOfPeriodicSubstring=true;
			if (startA>substrings[i].minStartA && Math.abs(endA,substrings[i].maxEndA)<=DISTANCE_THRESHOLD) ReadA.sortedAlignments[j].isSuffixOfPeriodicSubstring=true;
			j++;
		}
		
		// Setting $inPeriodicSubstring$
		rangeFirst=-1; rangeLast=-1; firstAlignment=0;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].hasLongPeriod) continue;
			if (rangeFirst==-1 || rangeLast==-1) {
				rangeFirst=intervals[i].firstPosition;
				rangeLast=intervals[i].lastPosition;
				continue;
			}
			else if (intervals[i].firstPosition<=rangeLast+DISTANCE_THRESHOLD) {
				if (intervals[i].lastPosition>rangeLast) rangeLast=intervals[i].lastPosition;
				continue;
			}
			// Current range
			firstJForNextI=-1;
			for (j=firstAlignment; j<=ReadA.lastSortedAlignment; j++) {
				startA=ReadA.sortedAlignments[j].startA();
				if (startA>=rangeLast) break;
				endA=ReadA.sortedAlignments[j].endA();
				if (endA<=rangeFirst) continue;
				if (firstJForNextI==-1 && endA>=intervals[i].firstPosition) firstJForNextI=j;
				if ( Intervals.isApproximatelyContained(startA,endA,rangeFirst,rangeLast) ||  
					 Intervals.areApproximatelyIdentical(startA,endA,rangeFirst,rangeLast) ||
					 Intervals.intersectionLength(startA,endA,rangeFirst,rangeLast)>=(endA-startA+1)*MIN_INTERSECTION_RATIO
				   ) ReadA.sortedAlignments[j].inPeriodicSubstring=true;
			}
			// Next range
			rangeFirst=intervals[i].firstPosition;
			rangeLast=intervals[i].lastPosition;
			if (firstJForNextI!=-1) firstAlignment=firstJForNextI;
		}
		// Last range
		for (j=firstAlignment; j<=ReadA.lastSortedAlignment; j++) {
			startA=ReadA.sortedAlignments[j].startA();
			if (startA>=rangeLast) break;
			endA=ReadA.sortedAlignments[j].endA();
			if (endA<=rangeFirst) continue;
			if ( Intervals.isApproximatelyContained(startA,endA,rangeFirst,rangeLast) ||  
				 Intervals.areApproximatelyIdentical(startA,endA,rangeFirst,rangeLast) ||
				 Intervals.intersectionLength(startA,endA,rangeFirst,rangeLast)>=(endA-startA+1)*MIN_INTERSECTION_RATIO
			   ) ReadA.sortedAlignments[j].inPeriodicSubstring=true;
		}
	}
	
	
	/**
	 * Remark: the procedure checks also the first multiple of the period from the left
	 * (respectively, right), in addition to the start (respectively, end) of the periodic
	 * substring. This is because such multiples are likely to have the second-largest 
	 * value of delta (built on nonimplied maximal alignments) among all multiples, thus
	 * creating a factor boundary, and such boundary is likely to have a large number of 
	 * noisy events because it is close to the beginning/end of the periodic substring.
	 *
	 * Remark: the procedure assumes that field $period$ of every periodic substring has
	 * already been computed.
	 *
	 * Remark: the procedure could work on periodic substring intervals rather than on
	 * periodic substrings.
	 *
	 * Remark: the procedure does not assume any order in $substrings$.
	 *
	 * @return TRUE iff the distance between $x$ and the start (if $left=true$) of a 
	 * periodic substring is at most $threshold$.
	 */
	public static final boolean atPeriodicSubstringEnd(int x, boolean left, int threshold) {
		int i;
		
		if (left) {
			for (i=0; i<=lastSubstring; i++) {
				if ( (substrings[i].nStartA>0 && Math.abs(x,substrings[i].sumStartA/substrings[i].nStartA)<=threshold) ||
					 Math.abs(x,substrings[i].minStartA)<=threshold ||
					 (substrings[i].nStartA>0 && Math.abs(x,substrings[i].sumStartA/substrings[i].nStartA+substrings[i].period)<=threshold) ||
					 Math.abs(x,substrings[i].minStartA+substrings[i].period)<=threshold ) return true;
			}
		}
		else {
			for (i=0; i<=lastSubstring; i++) {
				if ( (substrings[i].nEndA>0 && Math.abs(x,substrings[i].sumEndA/substrings[i].nEndA)<=threshold) ||
					 Math.abs(x,substrings[i].minEndA)<=threshold ||
					 (substrings[i].nEndA>0 && Math.abs(x,substrings[i].sumEndA/substrings[i].nEndA-substrings[i].period)<=threshold) ||
					 Math.abs(x,substrings[i].minEndA-substrings[i].period)<=threshold ) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Sets the $periodicSubstringInterval$ pointer of all alignments that are implied by
	 * a short-period periodic substring. If such periodic substring points to a short-
	 * period interval, and if the pointed interval is not discarded, the alignment points
	 * to the same interval. Otherwise, the alignment points to the shortest short-period 
	 * interval that approximately contains it.
	 *
	 * Remark: the procedure assumes $lastInterval>=0$.
	 */
	private static final void setIntervalPointers() {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final double INTERSECTION_THRESHOLD = 0.9;  // Arbitrary
		int i, j;
		int firstJForNextI, startA, endA, length;
		int minLength, minInterval, minLengthPrime, minIntervalPrime;
		PeriodicSubstring substring;
		
		// Ensuring the necessary order in $ReadA.sortedAlignments$ and in $intervals$
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		
		// Updating pointers from alignments to intervals
		i=0; j=0; firstJForNextI=-1;
		minLength=Math.POSITIVE_INFINITY; minInterval=-1;
		minLengthPrime=Math.POSITIVE_INFINITY; minIntervalPrime=-1;
		while (i<=ReadA.lastSortedAlignment) {
			substring=ReadA.sortedAlignments[i].impliedByPeriodicSubstring;
			if (substring==null) {
				i++;
				if (firstJForNextI>=0) j=firstJForNextI;
				firstJForNextI=-1;
				minLength=Math.POSITIVE_INFINITY; minInterval=-1;
				minLengthPrime=Math.POSITIVE_INFINITY; minIntervalPrime=-1;
				continue;
			}
			if ( (!substring.hasLongPeriod && substring.mergedToInterval!=null && !substring.mergedToInterval.discarded) ||
				 (substring.hasLongPeriod && substring.mergedToInterval!=null && !substring.mergedToInterval.hasLongPeriod && !substring.mergedToInterval.discarded)
			   ) {
				ReadA.sortedAlignments[i].periodicSubstringInterval=substring.mergedToInterval;
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				ReadA.sortedAlignments[i].inDenseSubstring=null;
				ReadA.sortedAlignments[i].mergedToInterval=null;
				i++;
				if (firstJForNextI>=0) j=firstJForNextI;
				firstJForNextI=-1;
				minLength=Math.POSITIVE_INFINITY; minInterval=-1;
				minLengthPrime=Math.POSITIVE_INFINITY; minIntervalPrime=-1;
				continue;
			}
			if (substring.hasLongPeriod && substring.mergedToInterval!=null && substring.mergedToInterval.hasLongPeriod) {
				i++;
				if (firstJForNextI>=0) j=firstJForNextI;
				firstJForNextI=-1;
				minLength=Math.POSITIVE_INFINITY; minInterval=-1;
				minLengthPrime=Math.POSITIVE_INFINITY; minIntervalPrime=-1;
				continue;
			}
			if (substring.mergedToInterval==null && (ReadA.sortedAlignments[i].periodicSubstringInterval==null || ReadA.sortedAlignments[i].periodicSubstringInterval.hasLongPeriod)) {
				i++;
				if (firstJForNextI>=0) j=firstJForNextI;
				firstJForNextI=-1;
				minLength=Math.POSITIVE_INFINITY; minInterval=-1;
				minLengthPrime=Math.POSITIVE_INFINITY; minIntervalPrime=-1;
				continue;
			}
			endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (j>lastInterval || intervals[j].firstPosition>=endA) {
				if (minInterval!=-1) {
					ReadA.sortedAlignments[i].periodicSubstringInterval=intervals[minInterval];
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].inDenseSubstring=null;
					ReadA.sortedAlignments[i].mergedToInterval=null;
				}
				else if (minIntervalPrime!=-1) {
					ReadA.sortedAlignments[i].periodicSubstringInterval=intervals[minIntervalPrime];
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].inDenseSubstring=null;
					ReadA.sortedAlignments[i].mergedToInterval=null;
				}
				else ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				i++;
				if (firstJForNextI>=0) j=firstJForNextI;
				firstJForNextI=-1;
				minLength=Math.POSITIVE_INFINITY; minInterval=-1;
				minLengthPrime=Math.POSITIVE_INFINITY; minIntervalPrime=-1;
				continue;
			}
			startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			if (firstJForNextI==-1 && i<ReadA.lastSortedAlignment && intervals[j].lastPosition>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextI=j;
			length=intervals[j].lastPosition-intervals[j].firstPosition+1;
			if ( ( Intervals.isApproximatelyContained_lowQuality(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id) ||
				   Intervals.areApproximatelyIdentical_lowQuality(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id)
				 ) && 
				 !( (startA>intervals[j].firstPosition+DISTANCE_THRESHOLD && (!Reads.isLeftMaximal(startA,ReadA.id,true) || ReadA.sortedAlignments[i].lowQualityStart)) &&
				    (endA<intervals[j].lastPosition-DISTANCE_THRESHOLD && (!Reads.isRightMaximal(endA,ReadA.id,true) || ReadA.sortedAlignments[i].lowQualityEnd))
				 ) &&
			     length<minLength
			   ) {
				minLength=length;
				minInterval=j;
			}
			if (Intervals.intersectionLength(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition)>=INTERSECTION_THRESHOLD*(endA-startA+1) && length<minLengthPrime) {
				minLengthPrime=length;
				minIntervalPrime=j;
			}
			j++;
		}
	}
	
	
	/**
	 * Marks as implied by a short-period interval all alignments that are approximately
	 * identical to it, or that are strictly contained in it and either cover a large 
	 * fraction of it, or they are not maximal (even though the alignments do not form a 
	 * periodic pattern).
	 *
	 * Remark: this criterion makes sense only for alignments that are long WRT their 
	 * containing short-period interval. Short alignments that do not form a periodic 
	 * pattern should be factorized separately.
	 *
	 * Remark: the procedure sets the $impliedByPeriodicSubstring$ field of alignments to 
	 * an artificial periodic substring that points to the short-period interval. We do 
	 * this, rather than setting the $periodicSubstringInterval$ field of alignments, for 
	 * compatibility with the following steps of the pipeline. Artificial periodic 
	 * substrings are added to the lists of substrings in the range 
	 * $[lastSubstring+1..lastSubstring_artificial]$.
	 *
	 * @param intervalsComposition TRUE: $intervals$ contains only short-period intervals;
	 * FALSE: $intervals$ contains both short- and long-period intervals.
	 */
	private static final void markImpliedAlignments_shortPeriodIntervals(boolean intervalsComposition) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final double LENGTH_THRESHOLD = 0.7;  // Arbitrary
		int i, j;
		int firstJForNextI, minInterval, minLength, intervalLength;
		int intervalStartA, intervalEndA, alignmentStartA, alignmentEndA;
		PeriodicSubstring artificialSubstring;
		
		// Ensuring the necessary order in $ReadA.sortedAlignments$ and in $intervals$.
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION || (!intervalsComposition && PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION)) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			if (!intervalsComposition) PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		
		// Marking
		i=0; j=0; firstJForNextI=-1; 
		minInterval=-1; minLength=Math.POSITIVE_INFINITY;
		while (i<=ReadA.lastSortedAlignment) {
			if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) {
				i++;
				if (firstJForNextI>=0) j=firstJForNextI;
				firstJForNextI=-1;
				minInterval=-1; minLength=Math.POSITIVE_INFINITY;				
				continue;
			}
			intervalStartA=intervals[j].firstPosition;
			alignmentEndA=ReadA.sortedAlignments[i].endA();
			if (j>lastInterval || intervalStartA>=alignmentEndA) {
				if (minInterval!=-1) {
					artificialSubstring=getArtificialSubstring(intervals[minInterval]);
					if (artificialSubstring==null) artificialSubstring=addArtificialSubstring(intervals[minInterval]);
					ReadA.sortedAlignments[i].impliedByPeriodicSubstring=artificialSubstring;
				}
				i++;
				if (firstJForNextI>=0) j=firstJForNextI;
				firstJForNextI=-1;
				minInterval=-1; minLength=Math.POSITIVE_INFINITY;
				continue;
			}
			if (intervals[j].hasLongPeriod) {
				j++;
				continue;
			}
			intervalEndA=intervals[j].lastPosition;
			if (firstJForNextI==-1 && i<ReadA.lastSortedAlignment && intervalEndA>=ReadA.sortedAlignments[i+1].startA()) firstJForNextI=j;
			intervalLength=intervals[j].length();
			alignmentStartA=ReadA.sortedAlignments[i].startA();
			if ( Intervals.areApproximatelyIdentical_lowQuality(alignmentStartA,alignmentEndA,intervalStartA,intervalEndA,ReadA.id) ||
				 ( Intervals.isApproximatelyContained_lowQuality(alignmentStartA,alignmentEndA,intervalStartA,intervalEndA,ReadA.id) && 
			       ( ReadA.sortedAlignments[i].getALength()>=intervalLength*LENGTH_THRESHOLD ||
					 (ReadA.sortedAlignments[i].isLeftMaximal!=1 && ReadA.sortedAlignments[i].isRightMaximal!=1) ||
					 (Math.abs(alignmentStartA,intervalStartA)<=DISTANCE_THRESHOLD && ReadA.sortedAlignments[i].isRightMaximal!=1) ||
					 (Math.abs(alignmentEndA,intervalEndA)<=DISTANCE_THRESHOLD && ReadA.sortedAlignments[i].isLeftMaximal!=1)
				   ) &&
				   !( (alignmentStartA>intervalStartA+DISTANCE_THRESHOLD && (!Reads.isLeftMaximal(alignmentStartA,ReadA.id,true) || ReadA.sortedAlignments[i].lowQualityStart) ) &&
					  (alignmentEndA<intervalEndA-DISTANCE_THRESHOLD && (!Reads.isRightMaximal(alignmentEndA,ReadA.id,true) || ReadA.sortedAlignments[i].lowQualityEnd))
				   ) &&
				   intervalLength<minLength
				 )
		       ) {
				minInterval=j;
				minLength=intervalLength;
			}
			j++;
		}
	}
	
	
	public static final void initArtificialSubstrings() {
		lastSubstring_artificial=lastSubstring;
	}
	
	
	public static final PeriodicSubstring getArtificialSubstring(PeriodicSubstringInterval interval) {
		for (int i=lastSubstring+1; i<=lastSubstring_artificial; i++) {
			if (substrings[i].mergedToInterval==interval) return substrings[i];
		}
		return null;
	}
	
	
	public static final PeriodicSubstring addArtificialSubstring(PeriodicSubstringInterval interval) {
		final int GROWTH_RATE = 10;
		int i;
		
		lastSubstring_artificial++;
		if (lastSubstring_artificial==substrings.length) {
			PeriodicSubstring[] array = new PeriodicSubstring[substrings.length+GROWTH_RATE];
			System.arraycopy(substrings,0,array,0,lastSubstring+1);
			for (i=lastSubstring+1; i<array.length; i++) array[i] = new PeriodicSubstring();
			substrings=array;
		}
		substrings[lastSubstring_artificial].setFromInterval(interval);
		return substrings[lastSubstring_artificial];
	}
	
	
	/**
	 * Sets the $mergedToInterval$ field to $mergedToInterval.representative$, if the 
	 * latter is not null.
	 */
	private static final void updateArtificialSubstringsPointers(boolean skipNonDiscardedIntervals) {
		int i;
		int leftSplit, rightSplit;
		PeriodicSubstring substring;
		PeriodicSubstringInterval tmpInterval, representative;
		
		for (i=lastSubstring+1; i<=lastSubstring_artificial; i++) {
			substring=substrings[i];
			tmpInterval=substring.mergedToInterval;
			if (tmpInterval==null || (skipNonDiscardedIntervals&&!tmpInterval.discarded)) continue;
			representative=tmpInterval.representative;
			substring.mergedToInterval=representative;
		}
	}
	
	
	private static final void initArtificialSubstrings_longPeriod() {
		lastLongPeriodSubstring_artificial=lastLongPeriodSubstring;
	}
	
	
	private static final PeriodicSubstring getArtificialSubstring_longPeriod(PeriodicSubstringInterval interval) {
		for (int i=lastLongPeriodSubstring+1; i<=lastLongPeriodSubstring_artificial; i++) {
			if (longPeriodSubstrings[i].mergedToInterval==interval) return longPeriodSubstrings[i];
		}
		return null;
	}
	
	
	private static final PeriodicSubstring addArtificialSubstring_longPeriod(PeriodicSubstringInterval interval) {
		final int GROWTH_RATE = 10;
		int i;
		
		lastLongPeriodSubstring_artificial++;
		if (lastLongPeriodSubstring_artificial==longPeriodSubstrings.length) {
			PeriodicSubstring[] array = new PeriodicSubstring[longPeriodSubstrings.length+GROWTH_RATE];
			System.arraycopy(longPeriodSubstrings,0,array,0,lastLongPeriodSubstring+1);
			for (i=lastLongPeriodSubstring+1; i<array.length; i++) array[i] = new PeriodicSubstring();
			longPeriodSubstrings=array;
		}
		longPeriodSubstrings[lastLongPeriodSubstring_artificial].setFromInterval(interval);
		return longPeriodSubstrings[lastLongPeriodSubstring_artificial];
	}
	
	
	
	
	
	
	
	
	// ---------------------------------- LONG PERIODS -----------------------------------
	
	/**
	 * Simplified version of $detect()$ for periodic substrings with long periods.
	 *
	 * Remark: the procedure assumes that $lastLongPeriodSubstring$ has already been
	 * initialized.
	 *
	 * Remark: like short-period periodic substrings/intervals, long-period periodic 
	 * substrings/intervals can straddle in readA.
	 *
	 * Remark: the procedure creates new $PeriodicSubstring$ and 
	 * $PeriodicSubstringInterval$ objects, and stores them in separate lists, so it 
	 * should not interfere with procedure $detect()$. Only at the end, and only 
	 * intervals, are merged with those of regular periodic substrings.
	 *
	 * @return TRUE iff new short-period intervals have been created by the procedure
	 * (yes, "short-period" is correct, even though it might seem that "long-period" would
	 * be correct).
	 */
	private static final boolean getSubstrings_longPeriods() {
		final int MAX_PERIOD = Alignments.minAlignmentLength*10;  // Arbitrary
		final int EVENT_THRESHOLD = IO.quantum;
		final int MIN_TRIMMED_LENGTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		boolean newShortPeriodInterval, newShortPeriodIntervalPrime;
		int i, j;
		int id, firstAlignment, nAlignments, nImpliedAlignments;
		int startA1, endA1, startA2, endA2, startB1, endB1, startB2, endB2;
		int minStartA, maxStartA, minEndA, maxEndA, minStartB, maxEndB;
		int sumStartA, nStartA, sumEndA, nEndA;
		int readB, source, destination, pathLength, orientation;
		int nComponents, longestPathInComponent, longestAlignment, lastLongPeriodSubstringPrime;
		int newShortPeriodSubstrings;
		PeriodicSubstring tmpSubstring;

if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("getSubstrings_longPeriods> SHORT-PERIOD INTERVALS at the very beginning:   MAX_PERIOD="+MAX_PERIOD);
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
}


		// Ensuring the necessary order in $ReadA.sortedAlignments$
		if (Alignment.order!=Alignment.READB_ORIENTATION_STARTA_ENDA) {
			Alignment.order=Alignment.READB_ORIENTATION_STARTA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

		// Detecting periodic substrings
		firstAlignment=0;
		do {
			// Building the DAG of all alignments of readA with a specific readB, in a
			// specific orientation, and sorting it topologically.
			readB=Alignments.alignments[ReadA.sortedAlignments[firstAlignment].id][1]-1;			
			orientation=Alignments.alignments[ReadA.sortedAlignments[firstAlignment].id][2];
			i=firstAlignment+1;
			while ( i<=ReadA.lastSortedAlignment &&
			        Alignments.alignments[ReadA.sortedAlignments[i].id][1]-1==readB &&
			        Alignments.alignments[ReadA.sortedAlignments[i].id][2]==orientation ) i++;
			nAlignments=i-firstAlignment;
			buildDAG(firstAlignment,nAlignments,MIN_PERIOD_LONG,MAX_PERIOD,IDENTITY_THRESHOLD_PERIODIC_LONG);
			nComponents=DAG.getConnectedComponents(nAlignments,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stackPrime);
			if (nComponents<=0 || nComponents>nAlignments) {
				IO.printCriticalErr("Error while detecting periodic substrings: wrong number of connected components.");
				System.exit(1);
			}			
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nAlignments);
			i=DAG.topologicalSort(nAlignments,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting periodic substrings: the directed graph of all alignments between read "+ReadA.id+" and read "+readB+" contains a cycle that involves node "+(i-1));
				System.exit(1);
			}

			// Computing maximal candidate repeated periodic substrings in readA induced
			// by readB, considering every alignment of readA with readB as a possible
			// source alignment. The intervals of such substrings in readA can intersect.
			for (source=0; source<nAlignments; source++) {
				id=ReadA.sortedAlignments[firstAlignment+source].id;
				if (ReadA.sortedAlignments[firstAlignment+source].impliedByPeriodicSubstring!=null) continue;
				getSubstring(firstAlignment,source,nAlignments,DEFAULT_MIN_PATH_LENGTH,IDENTITY_THRESHOLD_PERIODIC_LONG);
				destination=tmp[0]; pathLength=tmp[1];						
				if (destination==-1) continue;
				getAlignmentsInLongestPath(source,destination);
				lastLongPeriodSubstring++;
				minStartA=Alignments.alignments[id][3];
				maxStartA=minStartA; nStartA=1; sumStartA=minStartA;
				minEndA=Alignments.alignments[id][4];
				maxEndA=minEndA; nEndA=0; sumEndA=0;
				minStartB=Alignments.alignments[id][5];
				maxEndB=Alignments.alignments[id][6];
				for (i=0; i<pathLength+1; i++) {
					if ( ReadA.sortedAlignments[firstAlignment+alignments[i]].impliedByPeriodicSubstring!=null &&
						 !ReadA.sortedAlignments[firstAlignment+alignments[i]].impliedByPeriodicSubstring.hasLongPeriod
					   ) continue;
					id=ReadA.sortedAlignments[firstAlignment+alignments[i]].id;
					startA1=Alignments.alignments[id][3];
					endA1=Alignments.alignments[id][4];
					startB1=Alignments.alignments[id][5];
					endB1=Alignments.alignments[id][6];
					if (startA1<minStartA) minStartA=startA1;
					if (endA1>maxEndA) maxEndA=endA1;
					if (startB1<minStartB) minStartB=startB1;
					if (endB1>maxEndB) maxEndB=endB1;
					if (alignmentStates[i]==0 || alignmentStates[i]==2 || alignmentStates[i]==6) {
						if (startA1>maxStartA) maxStartA=startA1;
						nStartA++; sumStartA+=startA1;
					}
					if (alignmentStates[i]==2 || alignmentStates[i]==3 || alignmentStates[i]==5 || alignmentStates[i]==6 || alignmentStates[i]==7 || alignmentStates[i]==8 || alignments[i]==destination) {
						if (endA1<minEndA) minEndA=endA1;
						nEndA++; sumEndA+=endA1;
					}
				}
				longPeriodSubstrings[lastLongPeriodSubstring].hasLongPeriod=true;
				longPeriodSubstrings[lastLongPeriodSubstring].equalsOnePeriod=equalsOnePeriod(firstAlignment,minStartA,maxEndA,minStartB,maxEndB,IDENTITY_THRESHOLD_PERIODIC_LONG);
				longPeriodSubstrings[lastLongPeriodSubstring].minStartA=minStartA;
				longPeriodSubstrings[lastLongPeriodSubstring].maxStartA=maxStartA;
				longPeriodSubstrings[lastLongPeriodSubstring].nStartA=nStartA;
				longPeriodSubstrings[lastLongPeriodSubstring].sumStartA=sumStartA;
				longPeriodSubstrings[lastLongPeriodSubstring].minEndA=minEndA;
				longPeriodSubstrings[lastLongPeriodSubstring].maxEndA=maxEndA;
				longPeriodSubstrings[lastLongPeriodSubstring].nEndA=nEndA;
				longPeriodSubstrings[lastLongPeriodSubstring].sumEndA=sumEndA;
				longPeriodSubstrings[lastLongPeriodSubstring].minStartBForward=minStartB;
				longPeriodSubstrings[lastLongPeriodSubstring].maxEndBForward=maxEndB;
				longPeriodSubstrings[lastLongPeriodSubstring].readB=readB;
				longPeriodSubstrings[lastLongPeriodSubstring].pathLength=pathLength;
				if (orientation==1) {
					longPeriodSubstrings[lastLongPeriodSubstring].orientation=true;
					longPeriodSubstrings[lastLongPeriodSubstring].minStartB=longPeriodSubstrings[lastLongPeriodSubstring].minStartBForward;
					longPeriodSubstrings[lastLongPeriodSubstring].maxEndB=longPeriodSubstrings[lastLongPeriodSubstring].maxEndBForward;
				}
				else {
					longPeriodSubstrings[lastLongPeriodSubstring].orientation=false;
					longPeriodSubstrings[lastLongPeriodSubstring].minStartB=Reads.getReadLength(readB)-longPeriodSubstrings[lastLongPeriodSubstring].maxEndBForward-1;
					longPeriodSubstrings[lastLongPeriodSubstring].maxEndB=Reads.getReadLength(readB)-longPeriodSubstrings[lastLongPeriodSubstring].minStartBForward-1;
				}
				collectShifts(firstAlignment,longPeriodSubstrings[lastLongPeriodSubstring]);
				longPeriodSubstrings[lastLongPeriodSubstring].period=estimateLongPeriod(longPeriodSubstrings[lastLongPeriodSubstring].shifts,longPeriodSubstrings[lastLongPeriodSubstring].lastShift,null);
				// Marking as implied just alignments in the longest path from $source$ to
				// $destination$. We don't call $markImpliedAlignments()$ because we want
				// to assign alignments induced by repeats contained in the long-period
				// substring, to such contained repeats rather than to the container.
				nImpliedAlignments=0; longestAlignment=0;
				for (i=0; i<=lastAlignment; i++) {
					if (ReadA.sortedAlignments[firstAlignment+alignments[i]].impliedByPeriodicSubstring==null) {
						ReadA.sortedAlignments[firstAlignment+alignments[i]].impliedByPeriodicSubstring=longPeriodSubstrings[lastLongPeriodSubstring];
						longestAlignment=Math.max(longestAlignment,ReadA.sortedAlignments[firstAlignment+alignments[i]].getALength());
						nImpliedAlignments++;
					}
				}
				if (nImpliedAlignments<DEFAULT_MIN_PATH_LENGTH) {
					for (i=0; i<=lastAlignment; i++) {
						if (ReadA.sortedAlignments[firstAlignment+alignments[i]].impliedByPeriodicSubstring==longPeriodSubstrings[lastLongPeriodSubstring]) ReadA.sortedAlignments[firstAlignment+alignments[i]].impliedByPeriodicSubstring=null;
					}
					lastLongPeriodSubstring--;
				}
				else {
					longPeriodSubstrings[lastLongPeriodSubstring].longestAlignment=longestAlignment;
					longPeriodSubstrings[lastLongPeriodSubstring].trim(tmp,MIN_TRIMMED_LENGTH);
				}
			}
			firstAlignment+=nAlignments;
		}
		while (firstAlignment<=ReadA.lastSortedAlignment);		
		if (lastLongPeriodSubstring==-1) {
			lastLongPeriodInterval=-1;
			return false;
		}
		discardAlignmentsAfterTrimming(0,ReadA.lastSortedAlignment,true);
		Alignment.order=Alignment.STARTA_ENDA;
		if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		for (i=0; i<=lastLongPeriodSubstring; i++) longPeriodSubstrings[i].mergedToInterval=null;
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("LONG-PERIOD SUBSTRINGS IMMEDIATELY AFTER DETECTION: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
	IO.printErr("");
	System.err.println("alignments immediately after detection:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
}
		
		// Filtering and merging periodic substrings
		PeriodicSubstring.order_longPeriod=PeriodicSubstring.READB_ORIENTATION_MINSTARTA;
		if (lastLongPeriodSubstring>0) Arrays.sort(longPeriodSubstrings,0,lastLongPeriodSubstring+1);
		PeriodicSubstring.order=PeriodicSubstring.READB_ORIENTATION_MINSTARTA;
		if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		discardLongPeriodSubstrings();

		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("LONG-PERIOD SUBSTRINGS AFTER DISCARD 1: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
	IO.printErr("");
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
}
		
		
		newShortPeriodSubstrings=longPeriod2shortPeriod_substrings();
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("LONG-PERIOD SUBSTRINGS AFTER DISCARD 1.5: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")  newShortPeriodSubstrings="+newShortPeriodSubstrings);
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
	IO.printErr("");
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
}		
		
		if (newShortPeriodSubstrings>0) {
			PeriodicSubstring.order_longPeriod=PeriodicSubstring.LONGPERIOD_SHORTPERIOD;
			PeriodicSubstring.order=PeriodicSubstring.LONGPERIOD_SHORTPERIOD;
			if (lastLongPeriodSubstring>0) Arrays.sort(longPeriodSubstrings,0,lastLongPeriodSubstring+1);
			lastLongPeriodSubstringPrime=lastLongPeriodSubstring-newShortPeriodSubstrings;
		}
		else lastLongPeriodSubstringPrime=lastLongPeriodSubstring;
		correctCoverage(longPeriodSubstrings,lastLongPeriodSubstringPrime);
		j=lastLongPeriodSubstring;
		for (i=0; i<newShortPeriodSubstrings; i++) {
			j++;
			if (j==lastLongPeriodSubstringPrime+1+i) continue;
			tmpSubstring=longPeriodSubstrings[j];
		 	longPeriodSubstrings[j]=longPeriodSubstrings[lastLongPeriodSubstringPrime+1+i];
			longPeriodSubstrings[lastLongPeriodSubstringPrime+1+i]=tmpSubstring;
		}
		lastLongPeriodSubstringPrime=lastLongPeriodSubstring;
		lastLongPeriodSubstring=j;
		for (i=lastLongPeriodSubstringPrime+1; i<=lastLongPeriodSubstring; i++) {
			longPeriodSubstrings[i].shouldBeCountedForCoverage=true;
			longPeriodSubstrings[i].threadStart=longPeriodSubstrings[i];
			longPeriodSubstrings[i].equalsReverseComplement=false;
		}
		PeriodicSubstring.order_longPeriod=PeriodicSubstring.MINSTARTA_MAXENDA;
		if (newShortPeriodSubstrings>0) PeriodicSubstring.order=PeriodicSubstring.MINSTARTA_MAXENDA;
		if (lastLongPeriodSubstring>0) Arrays.sort(longPeriodSubstrings,0,lastLongPeriodSubstring+1);
		initArtificialSubstrings_longPeriod();
		if (IO.SHOW_STD_ERR) System.err.println((lastLongPeriodSubstring+1)+" long-period periodic substrings detected");


if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("LONG-PERIOD SUBSTRINGS AFTER DISCARD 2: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
	IO.printErr("");
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
}

		// Building intervals		
		for (i=0; i<=lastLongPeriodSubstring; i++) longPeriodSubstrings[i].mergedToInterval=null;
		lastLongPeriodInterval=-1;
		for (i=0; i<=lastLongPeriodSubstring; i++) {
			if (longPeriodSubstrings[i].mergedToInterval!=null) continue;		
			lastLongPeriodInterval=getInterval(longPeriodSubstrings,i,lastLongPeriodSubstring,longPeriodIntervals,lastLongPeriodInterval);
		}		
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
		if (lastLongPeriodInterval>0) Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("LONG-PERIOD INTERVALS --> 0:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
	System.err.println("SHORT-PERIOD INTERVALS --> 0:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
	IO.printErr("LONG-PERIOD SUBSTRINGS: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
}			
		
		
		newShortPeriodInterval=longPeriod2shortPeriod_moveIntervals();
		mergeLongPeriodIntervals();
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("LONG-PERIOD INTERVALS --> 1:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
	System.err.println("SHORT-PERIOD INTERVALS --> 1:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
	IO.printErr("LONG-PERIOD SUBSTRINGS: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
}		

		
		// Transforming long-period intervals into short-period intervals.
		discardLongPeriodIntervals(true,false);
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("LONG-PERIOD INTERVALS --> 1.5:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
	System.err.println("SHORT-PERIOD INTERVALS --> 1.5:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	IO.printErr("LONG-PERIOD SUBSTRINGS: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
}			
		
		newShortPeriodIntervalPrime=longPeriod2shortPeriod();
		if (newShortPeriodIntervalPrime) {
			newShortPeriodInterval|=newShortPeriodIntervalPrime;
			newShortPeriodInterval|=longPeriod2shortPeriod();
			// Executed again, since the procedure changes the set of short-period
			// intervals.
		}
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("LONG-PERIOD INTERVALS --> 2:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
	System.err.println("SHORT-PERIOD INTERVALS --> 2:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
	IO.printErr("LONG-PERIOD SUBSTRINGS: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
}		
		
		if (lastLongPeriodInterval==-1) {
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].periodicSubstringInterval==null && ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) {
					ReadA.sortedAlignments[i].periodicSubstringInterval=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.mergedToInterval;
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].inDenseSubstring=null;
					ReadA.sortedAlignments[i].mergedToInterval=null;
				}
			}
			return newShortPeriodInterval;
		}
		discardLongPeriodIntervals(false,newShortPeriodInterval);

		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("LONG-PERIOD INTERVALS --> 3:  newShortPeriodInterval="+newShortPeriodInterval);
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
	System.err.println("SHORT-PERIOD INTERVALS --> 3:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
	IO.printErr("LONG-PERIOD SUBSTRINGS: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
}					
		
		if (lastLongPeriodInterval==-1) {
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				if (ReadA.sortedAlignments[i].periodicSubstringInterval==null && ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) {
					ReadA.sortedAlignments[i].periodicSubstringInterval=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.mergedToInterval;
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].inDenseSubstring=null;
					ReadA.sortedAlignments[i].mergedToInterval=null;
				}
			}
			return newShortPeriodInterval;
		}
		setIntervalMaximality(longPeriodIntervals,lastLongPeriodInterval);		
		for (i=0; i<=lastLongPeriodInterval; i++) {
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("getSubstrings_longPeriods> estimating period of interval: "+longPeriodIntervals[i]);
	IO.printErr("OBSERVED SHIFTS (sorted and compacted):");
	for (int x=0; x<=longPeriodIntervals[i].lastShift; x++) IO.printErr(longPeriodIntervals[i].shifts[x].position+","+longPeriodIntervals[i].shifts[x].getMass()); 
}
			longPeriodIntervals[i].period=estimateLongPeriod(longPeriodIntervals[i].shifts,longPeriodIntervals[i].lastShift,longPeriodIntervals[i]);
		}
		markImpliedAlignments_longPeriodIntervals();



if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("LONG-PERIOD INTERVALS --> 3.5:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
	System.err.println("SHORT-PERIOD INTERVALS --> 3.5:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
	IO.printErr("LONG-PERIOD SUBSTRINGS: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
}	

		// The following is not needed for periodic substrings with long periods:
		// PeriodicSubstring.order_longPeriod=PeriodicSubstring.MINSTARTA_MAXENDA;
		// Arrays.sort(longPeriodSubstrings,0,lastLongPeriodSubstring+1);
		// markPrefixSuffix(longPeriodSubstrings,lastLongPeriodSubstring);
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if (ReadA.sortedAlignments[i].periodicSubstringInterval==null && ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) {
				ReadA.sortedAlignments[i].periodicSubstringInterval=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.mergedToInterval;  // Could be a short-period interval for a long-period substring, because of e.g. $discardLongPeriodIntervals()$.
				ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				ReadA.sortedAlignments[i].inDenseSubstring=null;
				ReadA.sortedAlignments[i].mergedToInterval=null;
			}
		}
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("LONG-PERIOD INTERVALS --> 3.6:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x].hashCode()+" :: "+longPeriodIntervals[x]);
	System.err.println("SHORT-PERIOD INTERVALS --> 3.6:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x].hashCode()+" :: "+intervals[x]);
	System.err.println("alignments:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
	IO.printErr("LONG-PERIOD SUBSTRINGS: (lastLongPeriodSubstring="+lastLongPeriodSubstring+")");
	for (int x=0; x<=lastLongPeriodSubstring; x++) IO.printErr(longPeriodSubstrings[x].hashCode()+" :: "+longPeriodSubstrings[x]);
}
	
	
if (IO.CONSISTENCY_CHECKS) {
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) {
		if (ReadA.sortedAlignments[x].impliedByPeriodicSubstring!=null) {
			if (!ReadA.sortedAlignments[x].impliedByPeriodicSubstring.hasLongPeriod && lastSubstring==-1) {
				System.err.println("The following alignment is assigned to a short-period periodic substring, but no such substring exists:");
				System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
				System.exit(1);
			}
			else if (ReadA.sortedAlignments[x].impliedByPeriodicSubstring.hasLongPeriod && lastLongPeriodSubstring==-1) {
				System.err.println("The following alignment is assigned to a long-period periodic substring, but no such substring exists:");
				System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
				System.exit(1);
			}
		}
	}
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) {
		if (ReadA.sortedAlignments[x].impliedByPeriodicSubstring!=null && !ReadA.sortedAlignments[x].impliedByPeriodicSubstring.hasLongPeriod && ReadA.sortedAlignments[x].periodicSubstringInterval==null && ReadA.sortedAlignments[x].impliedByPeriodicSubstring.mergedToInterval!=null) {
			System.err.println("(8) The following alignment is assigned to a short-period periodic substring but to no interval:");
			System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" -> "+ReadA.sortedAlignments[x].toStringPointers());
			System.exit(1);
		}
	}
	for (int x=0; x<=lastLongPeriodSubstring; x++) {
		for (int y=1; y<=longPeriodSubstrings[x].lastShift; y++) {
			if (longPeriodSubstrings[x].shifts[y].position<longPeriodSubstrings[x].shifts[y-1].position) {
				System.err.println("ERROR 2: shifts array not sorted.");
				System.exit(1);
			}
		}
	}
	for (int x=0; x<=lastSubstring; x++) {
		for (int y=1; y<=substrings[x].lastShift; y++) {
			if (substrings[x].shifts[y].position<substrings[x].shifts[y-1].position) {
				System.err.println("ERROR 3: shifts array not sorted.");
				System.exit(1);
			}
		}
	}
}
		
		return newShortPeriodInterval;		
	}
	
	
	/**
	 * Consider two long-period intervals as equivalent iff they straddle by at least S
	 * positions without being contained in one another. The procedure computes the 
	 * transitive closure of such equivalence relation, and merges all intervals in the 
	 * closure (as well as their shifts).
	 *
	 * Remark: the procedure merges only long-period intervals that are not contained in, 
	 * or identical to, a max range of straddling short-period intervals.
	 *
	 * Remark: the procedure assumes $longPeriodIntervals$ to be sorted by first position.
	 */
	private static final void mergeLongPeriodIntervals() {
		final int MIN_INTERSECTION_LENGTH = Alignments.minAlignmentLength;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i, j, k;
		int firstJForNextI, rangeFirst, rangeLast, lastShift;
		int top, firstPosition, lastPosition, longestAlignment;
		PeriodicSubstringInterval tmpInterval;
		Point[] tmpShifts;
		if (lastLongPeriodInterval<=0) return;
		
		// Marking long-period intervals contained in the union of short-period intervals.
		for (i=0; i<=lastLongPeriodInterval; i++) longPeriodIntervals[i].flag1=false;
		if (lastInterval!=-1) {
			rangeFirst=intervals[0].firstPosition; rangeLast=intervals[0].lastPosition;
			j=0; firstJForNextI=-1;
			for (i=1; i<=lastInterval; i++) {
				if (intervals[i].firstPosition<=rangeLast+IDENTITY_THRESHOLD) {
					rangeLast=Math.max(rangeLast,intervals[i].lastPosition);
					continue;
				}
				// Current range
				while (j<=lastLongPeriodInterval) {
					if (longPeriodIntervals[j].lastPosition<=rangeFirst) {
						j++;
						continue;
					}
					if (longPeriodIntervals[j].firstPosition>=rangeLast) break;
					if (firstJForNextI==-1 && longPeriodIntervals[j].lastPosition>rangeLast) firstJForNextI=j;
					if ( Intervals.areApproximatelyIdentical_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast,ReadA.id) ||
					     Intervals.areApproximatelyIdentical_lowCoverage(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast) ||
						 Intervals.isApproximatelyContained_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast,ReadA.id)
					   ) longPeriodIntervals[j].flag1=true;
					j++;
				}
				// Next range
				rangeFirst=intervals[i].firstPosition;
				rangeLast=intervals[i].lastPosition;
				if (firstJForNextI!=-1) j=firstJForNextI;
			}
			// Last range
			while (j<=lastLongPeriodInterval) {
				if (longPeriodIntervals[j].lastPosition<=rangeFirst) {
					j++;
					continue;
				}
				if (longPeriodIntervals[j].firstPosition>=rangeLast) break;
				if ( Intervals.areApproximatelyIdentical_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast,ReadA.id) ||
				     Intervals.areApproximatelyIdentical_lowCoverage(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast) ||
					 Intervals.isApproximatelyContained_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast,ReadA.id)
				   ) longPeriodIntervals[j].flag1=true;
				j++;
			}
		}
		
		// Merging intervals with large intersection
		for (i=0; i<=lastLongPeriodInterval; i++) longPeriodIntervals[i].representative=null;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].representative!=null || longPeriodIntervals[i].flag1) continue;
			rangeFirst=longPeriodIntervals[i].firstPosition;
			rangeLast=longPeriodIntervals[i].lastPosition;
			lastShift=Points.simpleClone(longPeriodIntervals[i].shifts,longPeriodIntervals[i].lastShift,shifts);
			longestAlignment=longPeriodIntervals[i].longestAlignment;
			stack[0]=i; top=0;
			while (top>=0) {
				j=stack[top--];
				firstPosition=longPeriodIntervals[j].firstPosition;
				lastPosition=longPeriodIntervals[j].lastPosition;
				for (k=j+1; k<=lastLongPeriodInterval; k++) {
					if (longPeriodIntervals[k].firstPosition>lastPosition-MIN_INTERSECTION_LENGTH) break;
					if (longPeriodIntervals[k].representative!=null || longPeriodIntervals[k].flag1 || longPeriodIntervals[k].firstPosition<=firstPosition+IDENTITY_THRESHOLD || longPeriodIntervals[k].lastPosition<=lastPosition+IDENTITY_THRESHOLD) continue;
					longPeriodIntervals[k].representative=longPeriodIntervals[i];
					longPeriodIntervals[i].simpleMerge(longPeriodIntervals[k]);
					rangeLast=Math.max(rangeLast,longPeriodIntervals[k].lastPosition);
					longestAlignment=Math.max(longestAlignment,longPeriodIntervals[k].longestAlignment);
					lastShift=Points.merge(shifts,lastShift,longPeriodIntervals[k].shifts,longPeriodIntervals[k].lastShift,shiftsPrime);
					tmpShifts=shifts;
					shifts=shiftsPrime;
					shiftsPrime=tmpShifts;
					stack[++top]=k;
				}
			}
			longPeriodIntervals[i].firstPosition=rangeFirst;
			longPeriodIntervals[i].lastPosition=rangeLast;
			longPeriodIntervals[i].longestAlignment=longestAlignment;
			longPeriodIntervals[i].cloneShifts(shifts,0,lastShift,MAX_DISTINCT_SHIFTS,Points.tmpPoints);
		}
		
		// Compacting intervals
		j=-1;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].representative!=null) continue;
			j++;
			if (j!=i) {
				tmpInterval=longPeriodIntervals[j];
				longPeriodIntervals[j]=longPeriodIntervals[i];
				longPeriodIntervals[i]=tmpInterval;
			}
		}
		lastLongPeriodInterval=j;
		
		// Merging intervals that are approximately identical
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].representative!=null) continue;
			rangeFirst=longPeriodIntervals[i].firstPosition;
			rangeLast=longPeriodIntervals[i].lastPosition;
			lastShift=Points.simpleClone(longPeriodIntervals[i].shifts,longPeriodIntervals[i].lastShift,shifts);
			longestAlignment=longPeriodIntervals[i].longestAlignment;
			stack[0]=i; top=0;
			while (top>=0) {
				j=stack[top--];
				firstPosition=longPeriodIntervals[j].firstPosition;
				lastPosition=longPeriodIntervals[j].lastPosition;
				for (k=j+1; k<=lastLongPeriodInterval; k++) {
					if (longPeriodIntervals[k].firstPosition>=lastPosition) break;
					if (!Intervals.areApproximatelyIdentical(longPeriodIntervals[k].firstPosition,longPeriodIntervals[k].lastPosition,firstPosition,lastPosition)) continue;
					longPeriodIntervals[k].representative=longPeriodIntervals[i];
					longPeriodIntervals[i].simpleMerge(longPeriodIntervals[k]);
					rangeLast=Math.max(rangeLast,longPeriodIntervals[k].lastPosition);
					longestAlignment=Math.max(longestAlignment,longPeriodIntervals[k].longestAlignment);
					lastShift=Points.merge(shifts,lastShift,longPeriodIntervals[k].shifts,longPeriodIntervals[k].lastShift,shiftsPrime);
					tmpShifts=shifts;
					shifts=shiftsPrime;
					shiftsPrime=tmpShifts;
					stack[++top]=k;
				}
			}
			longPeriodIntervals[i].lastPosition=rangeLast;
			longPeriodIntervals[i].longestAlignment=longestAlignment;
			longPeriodIntervals[i].cloneShifts(shifts,0,lastShift,MAX_DISTINCT_SHIFTS,Points.tmpPoints);
		}
		
		// Compacting intervals
		j=-1;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].representative!=null) continue;
			j++;
			if (j!=i) {
				tmpInterval=longPeriodIntervals[j];
				longPeriodIntervals[j]=longPeriodIntervals[i];
				longPeriodIntervals[i]=tmpInterval;
			}
		}
		lastLongPeriodInterval=j;
		
		// Updating pointers from periodic substrings
		for (i=0; i<=lastLongPeriodSubstring; i++) {
			if (!longPeriodSubstrings[i].hasLongPeriod) continue;
			tmpInterval=longPeriodSubstrings[i].mergedToInterval;
			if (tmpInterval==null) continue;
			while (tmpInterval.representative!=null) {
				longPeriodSubstrings[i].mergedToInterval=tmpInterval.representative;
				tmpInterval=tmpInterval.representative;
			}
		}
	}
	
	
	/**
	 * If $interval==null$, estimates the period of a long-period periodic substring from 
	 * its set of shifts. If $interval!=null$, estimates the period of a long-period 
	 * periodic substring interval, using the union of the shifts of all its periodic 
	 * substrings, as well as all the periodic substrings assigned to $interval$.
	 */
	private static final int estimateLongPeriod(Point[] shifts, int lastShift, PeriodicSubstringInterval interval) {
		final int MIN_INTERVAL_LENGTH = MIN_PERIOD_LONG;
		final int MIN_LOCAL_MAX_DISTANCE = MIN_PERIOD_LONG;
		int i, out;
		
		out=estimatePeriod(shifts,lastShift,MIN_PERIOD_LONG,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE);
		if (interval==null) return out;
		else {
			for (i=0; i<=lastLongPeriodSubstring; i++) {
				if (longPeriodSubstrings[i].mergedToInterval==interval) out=Math.max(out,longPeriodSubstrings[i].period);
			}
			return out;
		}
	}
	
	
	/**
	 * If a long-period interval L contains a short-period interval S, and if S is long 
	 * compared to the period of L, the procedure assigns to L the average of the 
	 * periods of all short-period intervals that are assigned to L, and it moves L to 
	 * $intervals$. I.e. the procedure transforms L into a short-period interval. The 
	 * shifts of S are also merged to those of L. If L instead just straddles short-period
	 * intervals, it is just marked as short-period and moved to $intervals$.
	 *
	 * Remark: transforming L into short-period is not allowed if L overlaps two short-
	 * period intervals that are adjacent, or if L contains two short-period intervals 
	 * separated by a high-quality non-periodic sequence. The presence of non-periodic 
	 * sequence implies that L does not have short period everywhere. The presence of two 
	 * adjacent intervals implies that L contains short-period substrings with distinct 
	 * periods (otherwise their intervals would have been merged). L is also not 
	 * transformed if it is contained in a maximal range of straddling or adjacent short-
	 * period intervals, or if it equals one copy of its period (a long-period interval 
	 * that contains a relatively long short-period interval should not be collapsed with 
	 * the latter).
	 *
	 * Remark: L is allowed to be transformed into short-period if it connects two short-
	 * period intervals such that at least one of them extends beyond L. This is because
	 * it might happen that a short-period interval is broken down into sub-intervals, 
	 * some of which are long-period.
	 *
	 * Remark: for simplicity, the boundaries of L are not updated with those of S.
	 *
	 * Remark: at the end of this procedure, a long-period substring might be assigned to
	 * a short-period interval.
	 *
	 * Remark: the procedure assumes that both lists of intervals are sorted by
	 * $firstPosition$.
	 *
	 * @return TRUE iff at least one long-period interval was transformed into
	 * short-period by the procedure.
	 */
	private static final boolean longPeriod2shortPeriod() {
		final double PERIOD_RATIO = 1.5;  // Arbitrary
		final int MIN_N_SHIFTS = 20*IO.coverage;  // Long-period substrings with fewer shifts than this are assumed not to have a reliable period estimate. Arbitrary value.
		final double LENGTH_RATIO = 0.3;  // Arbitrary
		final int DISTANCE_THRESHOLD = IO.quantum;  // Arbitrary
		final double INTERSECTION_RATIO = 0.5;  // Arbitrary
		boolean atLeastOneTransform, reliablePeriod;
		int i, j, k;
		int firstJForNextI, period, firstPosition, lastPosition, nAssigned, denominator, length, intersectionLength;
		int nTransformable, out, rangeFirst, rangeLast, surface, rangeIntervalFirst, intersectionFirst, intersectionLast;
		double avgPeriod;
		PeriodicSubstring substring;
		PeriodicSubstringInterval tmpInterval, interval, representative;
		Point[] tmpPoints;
		
		if (lastInterval==-1) return false;
		
		// Marking long-period intervals that can't be transformed into short-period (1/2)
		nTransformable=0;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].equalsOnePeriod) {
				longPeriodIntervals[i].canBeTransformed=false;
			}
			else {
				longPeriodIntervals[i].canBeTransformed=true;
				nTransformable++;
			}
		}
		if (nTransformable==0) return false;
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastLongPeriodInterval) {
			if (j>lastInterval || intervals[j].firstPosition>=longPeriodIntervals[i].lastPosition) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (intervals[j].lastPosition<=longPeriodIntervals[i].firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastLongPeriodInterval && intervals[j].lastPosition>=longPeriodIntervals[i+1].firstPosition) firstJForNextI=j;
			if ( Intervals.isApproximatelyContained(longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition,intervals[j].firstPosition,intervals[j].lastPosition) &&
		 	     !( (longPeriodIntervals[i].firstPosition>intervals[j].firstPosition+DISTANCE_THRESHOLD && !Reads.isLeftMaximal(longPeriodIntervals[i].firstPosition,ReadA.id,true)) &&
		 	        (longPeriodIntervals[i].lastPosition<intervals[j].lastPosition-DISTANCE_THRESHOLD && !Reads.isRightMaximal(longPeriodIntervals[i].lastPosition,ReadA.id,true))
		 	     )
			   ) {
   				longPeriodIntervals[i].canBeTransformed=false;
   				nTransformable--;
   				i++;
   				if (firstJForNextI!=-1) j=firstJForNextI;
   				firstJForNextI=-1;
   				continue;
			}
			if (Intervals.isApproximatelyContained(intervals[j].firstPosition,intervals[j].lastPosition,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition)) {
				if (intervals[j].lastPosition>=longPeriodIntervals[i].lastPosition-DISTANCE_THRESHOLD || !Reads.isRightMaximal(intervals[j].lastPosition,ReadA.id,true)) {
					j++;
					continue;
				}
				for (k=j+1; k<=lastInterval; k++) {
					if (intervals[k].firstPosition>longPeriodIntervals[i].lastPosition-DISTANCE_THRESHOLD) break;
					if (intervals[k].firstPosition<intervals[j].lastPosition-DISTANCE_THRESHOLD || !Reads.isLeftMaximal(intervals[k].firstPosition,ReadA.id,true)) continue;
					if ( intervals[k].firstPosition<=intervals[j].lastPosition+DISTANCE_THRESHOLD || 
					     ( Intervals.isApproximatelyContained(intervals[k].firstPosition,intervals[k].lastPosition,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition) &&
						   !Reads.isRandomInsertion(ReadA.id,intervals[j].lastPosition+1,intervals[k].firstPosition-1,true)
						 )
					   ) {
		   				longPeriodIntervals[i].canBeTransformed=false;
		   				nTransformable--;
		   				i++;
		   				if (firstJForNextI!=-1) j=firstJForNextI;
		   				firstJForNextI=-1;
		   				continue;
					}
				}
			}
			else if (intervals[j].lastPosition>longPeriodIntervals[i].firstPosition+DISTANCE_THRESHOLD && intervals[j].lastPosition<longPeriodIntervals[i].lastPosition-DISTANCE_THRESHOLD) {
				if (!Reads.isRightMaximal(intervals[j].lastPosition,ReadA.id,true)) {
					j++;
					continue;
				}
				for (k=j+1; k<=lastInterval; k++) {
					if (intervals[k].firstPosition>intervals[j].lastPosition+DISTANCE_THRESHOLD) break;
					if (intervals[k].firstPosition<intervals[j].lastPosition-DISTANCE_THRESHOLD || !Reads.isLeftMaximal(intervals[k].firstPosition,ReadA.id,true)) continue;
	   				longPeriodIntervals[i].canBeTransformed=false;
	   				nTransformable--;
	   				i++;
	   				if (firstJForNextI!=-1) j=firstJForNextI;
	   				firstJForNextI=-1;
	   				continue;
				}
			}
			j++;
		}
		if (nTransformable==0) return false;
		
		// Marking long-period intervals that can't be transformed into short-period (2/2)
		rangeIntervalFirst=0;
		rangeFirst=intervals[0].firstPosition;
		rangeLast=intervals[0].lastPosition;
		j=0;
		for (i=1; i<=lastInterval; i++) {
			if (intervals[i].firstPosition<=rangeLast+IO.quantum) {
				if (intervals[i].lastPosition>rangeLast) rangeLast=intervals[i].lastPosition;
				continue;
			}
			tmp[0]=j;
			longPeriod2shortPeriod_impl(rangeFirst,rangeLast,DISTANCE_THRESHOLD,tmp);
			j=tmp[0]; nTransformable-=tmp[1];
			// Next range
			rangeIntervalFirst=i;
			rangeFirst=intervals[i].firstPosition;
			rangeLast=intervals[i].lastPosition;
		}
		tmp[0]=j;
		longPeriod2shortPeriod_impl(rangeFirst,rangeLast,DISTANCE_THRESHOLD,tmp);
		nTransformable-=tmp[1];
		if (nTransformable==0) return false;	
		
		// Assigning short-period intervals to long-period intervals
		atLeastOneTransform=false;
		for (i=0; i<=lastInterval; i++) intervals[i].representative=null;
		for (i=0; i<=lastLongPeriodInterval; i++) longPeriodIntervals[i].representative=longPeriodIntervals[i];
		i=0; j=0; firstJForNextI=-1; avgPeriod=0.0; nAssigned=0; denominator=0;
		lastShift=Points.simpleClone(longPeriodIntervals[0].shifts,longPeriodIntervals[0].lastShift,shifts);
		firstPosition=longPeriodIntervals[0].firstPosition;
		lastPosition=longPeriodIntervals[0].lastPosition;
		length=longPeriodIntervals[0].length();
		reliablePeriod = longPeriodIntervals[0].period>0 && Points.mass(longPeriodIntervals[0].shifts,0,longPeriodIntervals[0].lastShift)>=MIN_N_SHIFTS;
		rangeFirst=-1; rangeLast=-1; surface=0;
		while (i<=lastLongPeriodInterval) {
			if (!longPeriodIntervals[i].canBeTransformed) {
				i++;
				avgPeriod=0.0; nAssigned=0; denominator=0;
				lastShift=Points.simpleClone(longPeriodIntervals[i].shifts,longPeriodIntervals[i].lastShift,shifts);
				firstPosition=longPeriodIntervals[i].firstPosition;
				lastPosition=longPeriodIntervals[i].lastPosition;
				length=longPeriodIntervals[i].length();
				reliablePeriod = longPeriodIntervals[i].period>0 && Points.mass(longPeriodIntervals[i].shifts,0,longPeriodIntervals[i].lastShift)>=MIN_N_SHIFTS;
				rangeFirst=-1; rangeLast=-1; surface=0;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (j>lastInterval || intervals[j].firstPosition>=lastPosition) {
				if (rangeLast!=-1) surface+=rangeLast-rangeFirst+1;
				if (nAssigned!=0) {
					longPeriodIntervals[i].hasLongPeriod=false;
					atLeastOneTransform=true;
					longPeriodIntervals[i].cloneShifts(shifts,0,lastShift,MAX_DISTINCT_SHIFTS,Points.tmpPoints);
					if (nAssigned==1 && avgPeriod!=0) longPeriodIntervals[i].period=(int)avgPeriod;
					else {
						period=estimatePeriod(longPeriodIntervals[i].shifts,longPeriodIntervals[i].lastShift,MIN_PERIOD,minIntervalLength,minLocalMaxDistance);
						if (period!=0) longPeriodIntervals[i].period=period;
						else if (avgPeriod!=0) longPeriodIntervals[i].period=(int)(avgPeriod/denominator);
					}
				}
				else if ( surface>0 &&
					      ( (reliablePeriod && surface>=longPeriodIntervals[i].period*PERIOD_RATIO) ||
      				        surface>=longPeriodIntervals[i].longestAlignment ||
   				            surface>=length*LENGTH_RATIO
						  )
				        ) {
					longPeriodIntervals[i].hasLongPeriod=false;
					atLeastOneTransform=true;
				}
				i++;
				avgPeriod=0.0; nAssigned=0; denominator=0;
				lastShift=Points.simpleClone(longPeriodIntervals[i].shifts,longPeriodIntervals[i].lastShift,shifts);
				firstPosition=longPeriodIntervals[i].firstPosition;
				lastPosition=longPeriodIntervals[i].lastPosition;
				length=longPeriodIntervals[i].length();
				reliablePeriod = longPeriodIntervals[i].period>0 && Points.mass(longPeriodIntervals[i].shifts,0,longPeriodIntervals[i].lastShift)>=MIN_N_SHIFTS;
				rangeFirst=-1; rangeLast=-1; surface=0;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (intervals[j].lastPosition<=firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastLongPeriodInterval && intervals[j].lastPosition>=longPeriodIntervals[i+1].firstPosition) firstJForNextI=j;
			if ( Intervals.isApproximatelyContained_lowQuality(intervals[j].firstPosition,intervals[j].lastPosition,firstPosition,lastPosition,ReadA.id) &&
 				 ( (reliablePeriod && intervals[j].length()>=longPeriodIntervals[i].period*PERIOD_RATIO) ||
 				   intervals[j].length()>=longPeriodIntervals[i].longestAlignment ||
 				   intervals[j].length()>=length*LENGTH_RATIO
				 )
			   ) {
				intervals[j].representative=longPeriodIntervals[i];
				longPeriodIntervals[i].simpleMerge(intervals[j]);
				nAssigned++;
				if (intervals[j].period!=0) {
					denominator++;
					avgPeriod+=intervals[j].period;
				}
				lastShift=Points.merge(shifts,lastShift,intervals[j].shifts,intervals[j].lastShift,shiftsPrime);
				tmpPoints=shifts;
				shifts=shiftsPrime;
				shiftsPrime=tmpPoints;
			}
			intersectionLength=Intervals.intersectionLength(intervals[j].firstPosition,intervals[j].lastPosition,firstPosition,lastPosition);
			if (intersectionLength>0) {
				intersectionFirst=Math.max(intervals[j].firstPosition,firstPosition);
				intersectionLast=Math.min(intervals[j].lastPosition,lastPosition);
				if (rangeLast==-1 || intersectionFirst>rangeLast+DISTANCE_THRESHOLD) {
					if (rangeLast!=-1) surface+=rangeLast-rangeFirst+1;
					rangeFirst=intersectionFirst;
					rangeLast=intersectionLast;
				}
				else rangeLast=Math.max(rangeLast,intersectionLast);
			}
			j++;
		}
		if (!atLeastOneTransform) {
			for (i=0; i<=lastLongPeriodInterval; i++) longPeriodIntervals[i].representative=null;
			return false;
		}
		
		for (i=0; i<=lastLongPeriodInterval; i++) longPeriodIntervals[i].representative=null;
		longPeriod2shortPeriod_moveIntervals();
		return true;
	}
	
	
	/**
	 * @param tmpArray input/output values: 0=first element of $longPeriodIntervals$ to 
	 * be used; at the end of the procedure, this cell contains the corresponding value 
	 * for the next invocation of the procedure;  1=number of long-period intervals for 
	 * which the procedure has set $canBeTransformed=false$.
	 */
	private static final void longPeriod2shortPeriod_impl(int rangeFirst, int rangeLast, int distanceThreshold, int[] tmpArray) {
		int j, firstJForNext, out;
		
		j=tmpArray[0]; firstJForNext=-1; out=0;
		while (j<=lastLongPeriodInterval) {
			if (longPeriodIntervals[j].firstPosition>=rangeLast) break;
			if (longPeriodIntervals[j].lastPosition<=rangeFirst) {
				j++;
				continue;
			}
			if (firstJForNext==-1 && longPeriodIntervals[j].lastPosition>rangeLast) firstJForNext=j;
			if (!longPeriodIntervals[j].canBeTransformed) {
				j++;
				continue;
			}
			if ( Intervals.areApproximatelyIdentical_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast,ReadA.id) ||
				 Intervals.areApproximatelyIdentical_lowCoverage(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast)
			   ) {
				longPeriodIntervals[j].canBeTransformed=false;
				out++;
			}
			else if ( Intervals.isApproximatelyContained(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast) &&
			          !( (longPeriodIntervals[j].firstPosition>rangeFirst+distanceThreshold && !Reads.isLeftMaximal(longPeriodIntervals[j].firstPosition,ReadA.id,true)) &&
			             (longPeriodIntervals[j].lastPosition<rangeLast-distanceThreshold && !Reads.isRightMaximal(longPeriodIntervals[j].lastPosition,ReadA.id,true))
			          )
				    ) {
				longPeriodIntervals[j].canBeTransformed=false;
				out++;
			}
			j++;
		}
		tmpArray[0]=firstJForNext==-1?j:firstJForNext;
		tmpArray[1]=out;
	}
	
	
	/**
	 * It might happen that a short-period periodic substring captures a long-period
	 * periodic substring. The procedure moves such substrings to $longPeriodSubstrings$
	 * and initializes $lastLongPeriodSubstring$.
	 *
	 * Remark: even after this, it might happen that a short-period periodic interval is 
	 * assigned a long period: see e.g. $discardLongPeriodSubstrings()$ and 
	 * $discardLongPeriodIntervals()$. E.g. the estimation of the period of a short-period
	 * interval might have failed, a similar long-period interval might have succeeded in 
	 * estimating the period, and the two intervals might have been merged. We choose to 
	 * output a short-period periodic interval with a long period, since we don't want to 
	 * lose the information that the interval might contain a short period. Interpreting 
	 * such annotation is left to the user.
	 */
	private static final void shortPeriod2longPeriod() {
		int i, j;
		PeriodicSubstring tmpSubstring;

		j=-1; lastLongPeriodSubstring=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].period>=MIN_PERIOD_LONG) {
				lastLongPeriodSubstring++;
				tmpSubstring=substrings[i];
				substrings[i]=longPeriodSubstrings[lastLongPeriodSubstring];
				substrings[i].hasLongPeriod=false;
				longPeriodSubstrings[lastLongPeriodSubstring]=tmpSubstring;
				longPeriodSubstrings[lastLongPeriodSubstring].hasLongPeriod=true;
			}
			else {
				j++;
				if (j!=i) {
					tmpSubstring=substrings[j];
					substrings[j]=substrings[i];
					substrings[i]=tmpSubstring;
				}
			}
		}
		lastSubstring=j;
		lastSubstringForFactoring=lastSubstring==-1?-1:Math.min(lastSubstringForFactoring,lastSubstring);
	}
	
	
	/**
	 * Discards long-period periodic substrings that are almost equivalent to a short-
	 * period periodic substring, updating the global variables $longPeriodSubstrings$ and
	 * $lastLongPeriodSubstring$.
	 *
	 * Remark: if a short-period periodic substring implies several long-period periodic
	 * substrings, and if it has no assigned period, the procedure assigns to it the max 
	 * period of all implied long-period substrings. However, the procedure does not merge
	 * the shifts of a long-period substring to those of the implying short-period
	 * substring.
	 *
	 * Remark: the procedure assumes $substrings$ and $longPeriodSubstrings$ to be sorted
	 * by $readB$, $orientation$, $minStartA$.
	 *
	 * @return the last long-period substring after the procedure completes.
	 */
	private static final void discardLongPeriodSubstrings() {
		final int DISTANCE_THRESHOLD = IO.quantum;
		boolean orientation;
		int i, j;
		int readB, minStartA, maxEndA, firstJForNextI;
		PeriodicSubstring tmp;

		for (i=0; i<=lastLongPeriodSubstring; i++) longPeriodSubstrings[i].representative=longPeriodSubstrings[i];
		for (i=0; i<=lastSubstring; i++) substrings[i].tmpInt=0;
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastLongPeriodSubstring) {
			readB=longPeriodSubstrings[i].readB;
			orientation=longPeriodSubstrings[i].orientation;
			minStartA=longPeriodSubstrings[i].minStartA;
			maxEndA=longPeriodSubstrings[i].maxEndA;
			if (j>lastSubstring || substrings[j].readB>readB || (substrings[j].orientation==false&&orientation==true) || substrings[j].minStartA>=maxEndA) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (substrings[j].readB<readB || (substrings[j].orientation==true&&orientation==false) || substrings[j].maxEndA<=minStartA) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastLongPeriodSubstring && longPeriodSubstrings[i+1].readB==readB && longPeriodSubstrings[i+1].orientation==orientation && substrings[j].maxEndA>=longPeriodSubstrings[i+1].minStartA) firstJForNextI=j;
			if ( Intervals.areApproximatelyIdentical_lowQuality(minStartA,maxEndA,substrings[j].minStartA,substrings[j].maxEndA,ReadA.id) &&
				 Intervals.areApproximatelyIdentical_lowQuality(longPeriodSubstrings[i].minStartB,longPeriodSubstrings[i].maxEndB,substrings[j].minStartB,substrings[j].maxEndB,readB)
			   ) {
				longPeriodSubstrings[i].representative=substrings[j];
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			j++;
		}
		j=-1;
		for (i=0; i<=lastLongPeriodSubstring; i++) {
			if (longPeriodSubstrings[i].representative!=longPeriodSubstrings[i]) {
				if (longPeriodSubstrings[i].representative.period==0) longPeriodSubstrings[i].representative.tmpInt=Math.max(longPeriodSubstrings[i].representative.tmpInt,longPeriodSubstrings[i].period);
				continue;
			}
			j++;
			if (j!=i) {
				tmp=longPeriodSubstrings[j];
				longPeriodSubstrings[j]=longPeriodSubstrings[i];
				longPeriodSubstrings[i]=tmp;
			}
		}
		lastLongPeriodSubstring=j;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].period==0 && substrings[i].tmpInt!=0) substrings[i].period=substrings[i].tmpInt;
		}
		
		// Updating pointers from alignments
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if ( ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null || 
				 !ReadA.sortedAlignments[i].impliedByPeriodicSubstring.hasLongPeriod ||
				 ReadA.sortedAlignments[i].impliedByPeriodicSubstring.representative==ReadA.sortedAlignments[i].impliedByPeriodicSubstring
			   ) continue;
			// Remark: to keep the code simple, intervals that get reassigned from a
			// long-period to a short-period substring are not used to change the
			// properties of the short-period substring or interval.
			ReadA.sortedAlignments[i].impliedByPeriodicSubstring=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.representative;
			ReadA.sortedAlignments[i].periodicSubstringInterval=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.mergedToInterval;
			ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
			ReadA.sortedAlignments[i].inDenseSubstring=null;
			ReadA.sortedAlignments[i].mergedToInterval=null;
		}
	}
	
	
	/**
	 * Discards long-period periodic intervals that are approx. equivalent to a short-
	 * period interval, or that are approx. equivalent to or contained in a maximal range
	 * of straddling or adjacent short-period intervals, updating the period of the short-
	 * period periodic intervals and the global variables $longPeriodIntervals$ and 
	 * $lastLongPeriodInterval$. At the end of the procedure, some long-period periodic 
	 * substrings might have their $mergedToInterval$ field set to a short-period 
	 * interval. After that, the procedure discards also long-period intervals that are 
	 * contained in another long-period interval, updating the set of shifts of the 
	 * container.
	 *
	 * Remark: the procedure does not merge the shifts of a long-period interval to those
	 * of the implying short-period interval, but it might update the period of the 
	 * short-period interval.
	 *
	 * Remark: the procedure assumes that $longPeriod2shortPeriod()$ has already been 
	 * executed, i.e. that no long-period interval has a short period.
	 *
	 * Remark: the procedure assumes $intervals$ and $longPeriodIntervals$ to be sorted
	 * by $firstPosition$.
	 *
	 * @param justEqualToShort if TRUE, the procedure removes just intervals that are 
	 * approx. identical to a short-period one;
	 * @param computePeaks recomputes the peaks of maximal events inside maximal short-
	 * period ranges (otherwise, the peaks are assumed to have already been computed); not
	 * used if $justEqualToShort=true$;
	 * @param containedInLong_intersectsShort TRUE: when trying to remove a long-period 
	 * interval X that is contained in another long-period interval Y, don't remove X if X
	 * or Y intersect a short-period interval.
	 */
	private static final void discardLongPeriodIntervals(boolean justEqualToShort, boolean computePeaks) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int INTERSECTION_THRESHOLD = IO.quantum<<2;  // Arbitrary
		final double JACCARD_THRESHOLD = 0.9;  // Arbitrary
		final int MAX_DISTANCE = Alignments.minAlignmentLength>>1;  // Arbitrary
		final int MIN_INTERSECTION_LENGTH = Alignments.minAlignmentLength;  // Arbitrary
		int i, j, k;
		int firstPosition, lastPosition, firstJForNextI, longestAlignment, top;
		int rangeFirst, rangeLast;
		PeriodicSubstringInterval tmpInterval;
		Point[] tmpShifts;

		if (!justEqualToShort && lastInterval>=0 && computePeaks) splitShortPeriodIntervals();

		// Short-period -> Long period
		for (i=0; i<=lastLongPeriodInterval; i++) {
			longPeriodIntervals[i].representative=null;
			longPeriodIntervals[i].intersectsShort=false;
		}
		for (i=0; i<=lastInterval; i++) {
			intervals[i].surface=0;  // Used as temporary space
			intervals[i].tmpInt1=intervals[i].firstPosition;  // Used as temporary space
			intervals[i].tmpInt2=intervals[i].lastPosition;  // Used as temporary space
		}
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastLongPeriodInterval) {
			firstPosition=longPeriodIntervals[i].firstPosition;
			lastPosition=longPeriodIntervals[i].lastPosition;
			if (j>lastInterval || intervals[j].firstPosition>=lastPosition) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (intervals[j].lastPosition<=firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastLongPeriodInterval && intervals[j].lastPosition>=longPeriodIntervals[i+1].firstPosition) firstJForNextI=j;
			if ( Intervals.areApproximatelyIdentical_lowQuality(firstPosition,lastPosition,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id) ||
				 ( ( Intervals.areApproximatelyIdentical_lowCoverage(firstPosition,lastPosition,intervals[j].firstPosition,intervals[j].lastPosition) ||
			         ( Intervals.jaccardSimilarity(firstPosition,lastPosition,intervals[j].firstPosition,intervals[j].lastPosition)>=JACCARD_THRESHOLD &&
				       Math.abs(firstPosition,intervals[j].firstPosition)<=MAX_DISTANCE && Math.abs(lastPosition,intervals[j].lastPosition)<=MAX_DISTANCE 
				     )
				   ) &&
				   position2split(firstPosition,leftSplits,lastLeftSplit,IDENTITY_THRESHOLD)==-1 && 
				   position2split(lastPosition,rightSplits,lastRightSplit,IDENTITY_THRESHOLD)==-1
				 )  // Passing the interval to the next steps if it aligns to peaks
			   ) {
				longPeriodIntervals[i].representative=intervals[j];
				intervals[j].simpleMerge(longPeriodIntervals[i]);
				intervals[j].tmpInt1=Math.min(intervals[j].tmpInt1,longPeriodIntervals[i].firstPosition);
				intervals[j].tmpInt2=Math.max(intervals[j].tmpInt2,longPeriodIntervals[i].lastPosition);
				longPeriodIntervals[i].intersectsShort=true;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (Intervals.intersectionLength(firstPosition,lastPosition,intervals[j].firstPosition,intervals[j].lastPosition)>=INTERSECTION_THRESHOLD) longPeriodIntervals[i].intersectsShort=true;
			j++;
		}
		j=-1;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].representative!=null) {
				if (longPeriodIntervals[i].representative.period==0) longPeriodIntervals[i].representative.surface=Math.max(longPeriodIntervals[i].representative.surface,longPeriodIntervals[i].period);
				continue;
			}
			j++;
			if (j!=i) {
				tmpInterval=longPeriodIntervals[j];
				longPeriodIntervals[j]=longPeriodIntervals[i];
				longPeriodIntervals[i]=tmpInterval;
			}
		}
		lastLongPeriodInterval=j;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].period==0 && intervals[i].surface!=0) intervals[i].period=intervals[i].surface;
			intervals[i].firstPosition=intervals[i].tmpInt1;
			intervals[i].lastPosition=intervals[i].tmpInt2;
		}

		// Updating pointers from periodic substrings
        for (i=0; i<=lastLongPeriodSubstring; i++) {
			if (!longPeriodSubstrings[i].hasLongPeriod) continue;
            tmpInterval=longPeriodSubstrings[i].mergedToInterval;
            if (tmpInterval==null) continue;
			while (tmpInterval.representative!=null) {
				longPeriodSubstrings[i].mergedToInterval=tmpInterval.representative;
				tmpInterval=tmpInterval.representative;
			}
        }
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("discardLongPeriodIntervals> 7  surviving long-period intervals:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x]);
	System.err.println("discardLongPeriodIntervals> 7  all intervals:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x]);
}

		
		// Short-period peaks -> Long period
		if (!justEqualToShort && lastInterval>=0) {
			lastLongPeriodInterval=filterIntervalsWithPeaks(longPeriodIntervals,lastLongPeriodInterval,longPeriodSubstrings,lastLongPeriodSubstring,true,false,true);
		}
		else {
			for (i=0; i<=lastLongPeriodInterval; i++) longPeriodIntervals[i].isContainedInShort=false;
		}
if (IO.SHOW_STD_ERR_PRIME) { 
	System.err.println("discardLongPeriodIntervals> 8  surviving long-period intervals:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x]);
	System.err.println("discardLongPeriodIntervals> 8  all intervals:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x]);
}
		
		
		// Marking long-period intervals contained in a max range of short-period
		// intervals.
		for (i=0; i<=lastLongPeriodInterval; i++) longPeriodIntervals[i].intersectsShort=false;
		i=0;
		rangeFirst=intervals[0].firstPosition;
		rangeLast=intervals[0].lastPosition;
		j=0;
		for (i=1; i<=lastInterval; i++) {
			if (intervals[i].firstPosition<=rangeLast+IDENTITY_THRESHOLD) {
				rangeLast=Math.max(rangeLast,intervals[i].lastPosition);
				continue;
			}
			while (j<=lastLongPeriodInterval) {
				if (longPeriodIntervals[j].firstPosition>=rangeLast) break;
				if (longPeriodIntervals[j].lastPosition<=rangeFirst) {
					j++;
					continue;
				}
				if ( Intervals.areApproximatelyIdentical(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast) ||
					 Intervals.isApproximatelyContained(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast)
				   ) longPeriodIntervals[j].intersectsShort=true;
				j++;
			}
			// Next range
			rangeFirst=intervals[i].firstPosition;
			rangeLast=intervals[i].lastPosition;
		}
		while (j<=lastLongPeriodInterval) {
			if (longPeriodIntervals[j].firstPosition>=rangeLast) break;
			if (longPeriodIntervals[j].lastPosition<=rangeFirst) {
				j++;
				continue;
			}
			if ( Intervals.areApproximatelyIdentical(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast) ||
				 Intervals.isApproximatelyContained(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,rangeFirst,rangeLast)
			   ) longPeriodIntervals[j].intersectsShort=true;
			j++;
		}

		// Long-period -> Long period
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].intersectsShort || longPeriodIntervals[i].representative!=null) continue;
			firstPosition=longPeriodIntervals[i].firstPosition;
			lastPosition=longPeriodIntervals[i].lastPosition;
			longestAlignment=longPeriodIntervals[i].longestAlignment;
			lastShift=Points.simpleClone(longPeriodIntervals[i].shifts,longPeriodIntervals[i].lastShift,shifts);
			for (j=i-1; j>=0; j--) {
				if (longPeriodIntervals[j].intersectsShort || longPeriodIntervals[j].representative!=null) continue;
				if ( ( Intervals.isApproximatelyContained_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,firstPosition,lastPosition,ReadA.id) ||
					   Intervals.areApproximatelyIdentical_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,firstPosition,lastPosition,ReadA.id)
					 ) /* &&  // Commented because this level of caution is not needed in practice.
			 	     !( (longPeriodIntervals[j].firstPosition>firstPosition+IDENTITY_THRESHOLD && !Reads.isLeftMaximal(longPeriodIntervals[j].firstPosition,ReadA.id,true)) &&
			 	        (longPeriodIntervals[j].lastPosition<lastPosition-IDENTITY_THRESHOLD && !Reads.isRightMaximal(longPeriodIntervals[j].lastPosition,ReadA.id,true))
			 	     )*/
				   ) {
					longPeriodIntervals[j].representative=longPeriodIntervals[i];
					longPeriodIntervals[i].simpleMerge(longPeriodIntervals[j]);
					lastShift=Points.merge(shifts,lastShift,longPeriodIntervals[j].shifts,longPeriodIntervals[j].lastShift,shiftsPrime);
					tmpShifts=shifts;
					shifts=shiftsPrime;
					shiftsPrime=tmpShifts;
				}
			}
			for (j=i+1; j<=lastLongPeriodInterval; j++) {
				if (longPeriodIntervals[j].firstPosition>=lastPosition) break;
				if (longPeriodIntervals[j].intersectsShort || longPeriodIntervals[j].representative!=null) continue;
				if ( ( Intervals.isApproximatelyContained_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,firstPosition,lastPosition,ReadA.id) ||
					   Intervals.areApproximatelyIdentical_lowQuality(longPeriodIntervals[j].firstPosition,longPeriodIntervals[j].lastPosition,firstPosition,lastPosition,ReadA.id)
					 ) /* &&  // Commented because this level of caution is not needed in practice.
		 	         !( (longPeriodIntervals[j].firstPosition>firstPosition+IDENTITY_THRESHOLD && !Reads.isLeftMaximal(longPeriodIntervals[j].firstPosition,ReadA.id,true)) &&
		 	            (longPeriodIntervals[j].lastPosition<lastPosition-IDENTITY_THRESHOLD && !Reads.isRightMaximal(longPeriodIntervals[j].lastPosition,ReadA.id,true))
		 	         )*/
				   ) {
					longPeriodIntervals[j].representative=longPeriodIntervals[i];
					longPeriodIntervals[i].simpleMerge(longPeriodIntervals[j]);
					lastShift=Points.merge(shifts,lastShift,longPeriodIntervals[j].shifts,longPeriodIntervals[j].lastShift,shiftsPrime);
					tmpShifts=shifts;
					shifts=shiftsPrime;
					shiftsPrime=tmpShifts;
				}
			}
			longPeriodIntervals[i].cloneShifts(shifts,0,lastShift,MAX_DISTINCT_SHIFTS,Points.tmpPoints);
		}
		j=-1;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (!longPeriodIntervals[i].intersectsShort && longPeriodIntervals[i].representative!=null) continue;
			j++;
			if (j!=i) {
				tmpInterval=longPeriodIntervals[j];
				longPeriodIntervals[j]=longPeriodIntervals[i];
				longPeriodIntervals[i]=tmpInterval;
			}
		}
		lastLongPeriodInterval=j;
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("discardLongPeriodIntervals> 10  surviving long-period intervals:");
	for (int x=0; x<=lastLongPeriodInterval; x++) System.err.println(longPeriodIntervals[x]);
	System.err.println("discardLongPeriodIntervals> 10  all intervals:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x]);
}
	
		// Updating pointers from periodic substrings
		for (i=0; i<=lastLongPeriodSubstring; i++) {
			if (!longPeriodSubstrings[i].hasLongPeriod) continue;
			tmpInterval=longPeriodSubstrings[i].mergedToInterval;
			if (tmpInterval==null || tmpInterval.isContainedInShort) continue;
			while (tmpInterval.representative!=null) {
				longPeriodSubstrings[i].mergedToInterval=tmpInterval.representative;
				tmpInterval=tmpInterval.representative;
			}
		}
	}
		
		
	/**
	 * Clones all elements of $longPeriodIntervals$ into $intervals$.
	 *
	 * Remark: this function must be undone before processing the next read, by calling
	 * $undo_mergeIntervals()$.
	 *
	 * Remark: long-period periodic substrings are not merged to short-period periodic
	 * substrings, so the rest of the pipeline that uses them does not need to be changed.
	 *
	 * Remark: the procedure does not merge-sort the two lists of intervals, but it just
	 * copies and sorts, to avoid reallocating $intervals$.
	 */
	private static final void mergeIntervals() {
		if (lastLongPeriodInterval==-1) return;
		final int NEW_LENGTH = lastInterval+lastLongPeriodInterval+2;
		int i;
		
		if (intervals.length<NEW_LENGTH) {
			PeriodicSubstringInterval[] tmp = new PeriodicSubstringInterval[NEW_LENGTH];
			System.arraycopy(intervals,0,tmp,0,lastInterval+1);
			intervals=tmp;
		}
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (intervals[lastInterval+1+i]==null) intervals[lastInterval+1+i] = new PeriodicSubstringInterval();
			intervals[lastInterval+1+i].clone(longPeriodIntervals[i]);
		}
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
		if (NEW_LENGTH>1) Arrays.sort(intervals,0,NEW_LENGTH);
		lastInterval=NEW_LENGTH-1;
	}
	
	
	/**
	 * Resets all $hasLongPeriod$ marks in $intervals$. This should be done before 
	 * factorizing the next read.
	 */
	public static final void undo_mergeIntervals() {
		for (int i=0; i<=lastInterval; i++) intervals[i].hasLongPeriod=false;
	}
	
	
	/**
	 * Copies back, to the objects in $longPeriodIntervals$, the ID fields of their clones 
	 * in $intervals$.
	 *
	 * Remark: the procedure assumes $intervals$ to be a superset of 
	 * $longPeriodIntervals$, and it assumes both arrays to be sorted by $firstPosition$.
	 */
	public static final void mergeIntervals_copyIDs() {
		int i, j, firstJForNextI;

		i=0; j=0; firstJForNextI=-1;
		while (i<=lastLongPeriodInterval) {
			if (j>lastInterval) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (firstJForNextI==-1 && i<lastLongPeriodInterval && intervals[j].firstPosition==longPeriodIntervals[i+1].firstPosition) firstJForNextI=j;
			if (intervals[j].firstPosition>longPeriodIntervals[i].firstPosition) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (!intervals[j].hasLongPeriod || intervals[j].firstPosition<longPeriodIntervals[i].firstPosition) {
				j++;
				continue;
			}
			if (intervals[j].lastPosition==longPeriodIntervals[i].lastPosition) {
				longPeriodIntervals[i].id=intervals[j].id;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			j++;
		}
	}

	
	/**
	 * Assigns to a long-period periodic substring interval, all alignments that 
	 * correspond to: (1) a short substring of the long period, that occurs in a readB 
	 * sorrounded by low-quality regions; (2) a short prefix or suffix of the long period,
	 * that occurs in a readB truncated on the other side by a low-quality region. Such 
	 * substrings should have multiple alignments with the same long-period periodic 
	 * interval of readA, but they should not display the pattern of alignments used to 
	 * detect periodic substrings. The procedure also assigns to a long-period interval
	 * all alignments whose readA part is approx. identical to the interval.
	 *
	 * Remark: the procedure applies also to long-period intervals contained in a short-
	 * period interval.
	 *
	 * Remark: it could happen that another instance of a long-period periodic substring
	 * occurs in a readB with a different phase. If this happens to a substring of the 
	 * long period, then this is not detected by the procedure, since it would create a 
	 * maximal event inside the periodic interval of readA.
	 */
	private static final void markImpliedAlignments_longPeriodIntervals() {
		final int COUNT_THRESHOLD = 2;  // Arbitrary
		int i;
		final int previousOrder = Alignment.order;
		Alignment alignment;
		PeriodicSubstring artificialSubstring;
		PeriodicSubstringInterval interval;
		PeriodicSubstringInterval tmpInterval = new PeriodicSubstringInterval();
		
		// Ensuring the necessary order in $ReadA.sortedAlignments$ and in
		// $longPeriodIntervals$.
		if (Alignment.order!=Alignment.READB_STARTB_ORIENTATION) {
			Alignment.order=Alignment.READB_STARTB_ORIENTATION;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		if (PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (lastLongPeriodInterval>0) Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
		}
		
		// Assigning alignments to intervals
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			alignment=ReadA.sortedAlignments[i];
			if (alignment.impliedByPeriodicSubstring!=null || alignment.periodicSubstringInterval!=null) continue;
			// First criterion
			interval=getLongPeriodInterval(alignment,tmpInterval);
			if (interval!=null) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("markImpliedAlignments_longPeriodIntervals> candidate source alignment: "+alignment.toStringBoundaries());
				if (hasCompatibleAlignments(i,interval,COUNT_THRESHOLD,tmpInterval)) assignCompatibleAlignments(i,interval,tmpInterval);
			}
			// Second criterion. 
			// We mark using an artificial periodic substring, rather than an interval:
			// see $assignCompatibleAlignments()$ for details.
			interval=getLongPeriodInterval_identical(alignment,tmpInterval);
			if (interval!=null) {
				artificialSubstring=getArtificialSubstring_longPeriod(interval);
				if (artificialSubstring==null) artificialSubstring=addArtificialSubstring_longPeriod(interval);
				ReadA.sortedAlignments[i].impliedByPeriodicSubstring=artificialSubstring;
				interval.nImpliedAlignments++;
				interval.longestAlignment=Math.max(interval.longestAlignment,alignment.getALength());
			}
		}
		
		// Restoring the original order in $ReadA.sortedAlignments$
		if (Alignment.order!=previousOrder) {
			Alignment.order=previousOrder;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
	}
	
	
	/**
	 * Computes the long-period interval with the following properties: (1) it contains 
	 * the readA part of $alignment$; (2) the readA part of $alignment$ is shorter than 
	 * its long period, if any; (3) the alignment does not have a maximality event inside
	 * the interval. The procedure returns NULL if $alignment$ does not belong to any 
	 * such long-period interval.
	 * 
	 * Remark: the procedure assumes $longPeriodIntervals$ to be sorted by 
	 * $firstPosition$, it assumes periods to have already been assigned to long-period 
	 * intervals, and it assumes that an alignment belongs to at most one long-period 
	 * interval in readA.
	 *
	 * @param tmpInterval temporary space.
	 */
	private static final PeriodicSubstringInterval getLongPeriodInterval(Alignment alignment, PeriodicSubstringInterval tmpInterval) {
		final double SHORTER_THAN_PERIOD_THRESHOLD = 0.8;  // Arbitrary
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int startA = Alignments.alignments[alignment.id][3];
		final int endA = Alignments.alignments[alignment.id][4];
		final int lengthA = endA-startA+1;
		final int orientation = Alignments.alignments[alignment.id][2];
		boolean found, previousSortByID;
		int i, p;
		
		if (alignment.isLeftMaximalB==1 && alignment.isRightMaximalB==1) return null;
		tmpInterval.firstPosition=startA;
		previousSortByID=PeriodicSubstringInterval.sortByID;
		PeriodicSubstringInterval.sortByID=false;
		p=Arrays.binarySearch(longPeriodIntervals,0,lastLongPeriodInterval+1,tmpInterval);
		PeriodicSubstringInterval.sortByID=previousSortByID;
		if (p<0) p=-p-1;
		for (i=p-1; i>=0; i--) {
if (IO.SHOW_STD_ERR_PRIME) System.err.print("getLongPeriodInterval> alignment=["+startA+".."+endA+"] :: interval=["+longPeriodIntervals[i].firstPosition+".."+longPeriodIntervals[i].lastPosition+"] ??? ");
			if (longPeriodIntervals[i].lastPosition<=startA) continue;
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.print(Intervals.isApproximatelyContained(startA,endA,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition)?"1,":"0,");
	System.err.print(Intervals.areApproximatelyIdentical(startA,endA,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition)?"1,":"0,");
	System.err.print((longPeriodIntervals[i].period==0 || lengthA<=longPeriodIntervals[i].period*SHORTER_THAN_PERIOD_THRESHOLD)?"1,":"0,");
	System.err.print((startA>longPeriodIntervals[i].firstPosition+DISTANCE_THRESHOLD && (orientation==1?alignment.isLeftMaximal==1:alignment.isRightMaximal==1))?"1,":"0,");
	System.err.print((endA<longPeriodIntervals[i].lastPosition-DISTANCE_THRESHOLD && (orientation==1?alignment.isRightMaximal==1:alignment.isLeftMaximal==1))?"1,":"0,");
	System.err.println();
}			
			if ( Intervals.isApproximatelyContained(startA,endA,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition) &&
				 !Intervals.areApproximatelyIdentical(startA,endA,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition) &&
				 (longPeriodIntervals[i].period==0 || lengthA<=longPeriodIntervals[i].period*SHORTER_THAN_PERIOD_THRESHOLD) &&
				 !( ( startA>longPeriodIntervals[i].firstPosition+DISTANCE_THRESHOLD && 
					  (alignment.isLeftMaximal==1 || !Reads.isLeftMaximal(startA,ReadA.id,true) || alignment.lowQualityStart)
					) ||
				    ( endA<longPeriodIntervals[i].lastPosition-DISTANCE_THRESHOLD && 
					  (alignment.isRightMaximal==1 || !Reads.isRightMaximal(endA,ReadA.id,true) || alignment.lowQualityEnd)
				    )
				  )
			   ) return longPeriodIntervals[i];
		}
		for (i=p; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].firstPosition>=endA) break;
			if ( Intervals.isApproximatelyContained(startA,endA,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition) &&
				 !Intervals.areApproximatelyIdentical(startA,endA,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition) &&
				 (longPeriodIntervals[i].period==0 || lengthA<=longPeriodIntervals[i].period*SHORTER_THAN_PERIOD_THRESHOLD) &&
				 !( ( startA>longPeriodIntervals[i].firstPosition+DISTANCE_THRESHOLD && 
					  (alignment.isLeftMaximal==1 || !Reads.isLeftMaximal(startA,ReadA.id,true) || alignment.lowQualityStart)
					) ||
				    ( endA<longPeriodIntervals[i].lastPosition-DISTANCE_THRESHOLD && 
					  (alignment.isRightMaximal==1 || !Reads.isRightMaximal(endA,ReadA.id,true) || alignment.lowQualityEnd)
				    )
				  )
			   ) return longPeriodIntervals[i];
		}
		return null;
	}
	
	
	/**
	 * Returns a shortest long-period interval that is approx. identical to the readA part
	 * of $alignment$, or NULL if no such interval can be found.
	 * 
	 * Remark: the procedure assumes $longPeriodIntervals$ to be sorted by 
	 * $firstPosition$.
	 *
	 * @param tmpInterval temporary space.
	 */
	private static final PeriodicSubstringInterval getLongPeriodInterval_identical(Alignment alignment, PeriodicSubstringInterval tmpInterval) {
		final int startA = Alignments.alignments[alignment.id][3];
		final int endA = Alignments.alignments[alignment.id][4];
		boolean previousSortByID;
		int i, p;
		int minInterval, minLength, length;
		
		tmpInterval.firstPosition=startA;
		previousSortByID=PeriodicSubstringInterval.sortByID;
		PeriodicSubstringInterval.sortByID=false;
		p=Arrays.binarySearch(longPeriodIntervals,0,lastLongPeriodInterval+1,tmpInterval);
		PeriodicSubstringInterval.sortByID=previousSortByID;
		if (p<0) p=-p-1;
		minInterval=-1; minLength=Math.POSITIVE_INFINITY;
		for (i=p-1; i>=0; i--) {
			if (longPeriodIntervals[i].lastPosition<=startA) continue;
			if (Intervals.areApproximatelyIdentical(startA,endA,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition)) {
				length=longPeriodIntervals[i].length();
				if (length<minLength) {
					minLength=length;
					minInterval=i;
				}
			}
		}
		for (i=p; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].firstPosition>=endA) break;
			if (Intervals.areApproximatelyIdentical(startA,endA,longPeriodIntervals[i].firstPosition,longPeriodIntervals[i].lastPosition)) {
				length=longPeriodIntervals[i].length();
				if (length<minLength) {
					minLength=length;
					minInterval=i;
				}
			}
		}
		return minInterval!=-1?longPeriodIntervals[minInterval]:null;
	}
	
	
	/**
	 * Decides if there are at least $countThreshold$ other alignments with the following
	 * properties: (1) they are approximately identical to $fromAlignment$ in readB; 
	 * (2) they belong to $fromInterval$ in readA, with the same orientation as 
	 * $fromAlignment$ and without creating maximal events inside $fromInterval$; 
	 * (3) they are not already implied by a periodic substring or assigned to a periodic 
	 * interval.
	 *
	 * @param countThreshold lower bound; the procedure can lower this value based on the 
	 * period of $fromInterval$; if the final threshold is zero, the procedure returns 
	 * true optimistically;
	 * @param fromAlignment ID in $ReadA.sortedAlignments$;
	 * @param fromInterval the long-period interval to which $fromAlignment$ belongs in
	 * readA;
	 * @param tmpInterval temporary space.
	 */
	private static final boolean hasCompatibleAlignments(int fromAlignment, PeriodicSubstringInterval fromInterval, int countThreshold, PeriodicSubstringInterval tmpInterval) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int SEARCH_THRESHOLD = IO.quantum<<1;
		int i;
		int id, readB, orientation, fromStartB, fromEndB, startA, endA, startB, endB, count;
		final int fromAlignmentID = ReadA.sortedAlignments[fromAlignment].id;
		Alignment otherAlignment;
		
		if (fromInterval.period!=0) countThreshold=Math.min(countThreshold+1,fromInterval.length()/fromInterval.period)-1;
		if (countThreshold<=0) return true;
		readB=Alignments.alignments[fromAlignmentID][1];
		orientation=Alignments.alignments[fromAlignmentID][2];
		fromStartB=Alignments.alignments[fromAlignmentID][5];
		fromEndB=Alignments.alignments[fromAlignmentID][6];
		count=0;
		for (i=fromAlignment-1; i>=0; i--) {
			otherAlignment=ReadA.sortedAlignments[i];
			id=otherAlignment.id;
			startB=Alignments.alignments[id][5]; endB=Alignments.alignments[id][6];
			if (Alignments.alignments[id][1]!=readB || startB<fromStartB-SEARCH_THRESHOLD) break;
			if (Alignments.alignments[id][2]!=orientation) continue;
			if (otherAlignment.isLeftMaximalB==1 && otherAlignment.isRightMaximalB==1) continue;
			if (otherAlignment.impliedByPeriodicSubstring!=null || otherAlignment.periodicSubstringInterval!=null) continue;
			startA=Alignments.alignments[id][3]; endA=Alignments.alignments[id][4];
			if ( Intervals.areApproximatelyIdentical(startB,endB,fromStartB,fromEndB) &&
				 Intervals.isApproximatelyContained(startA,endA,fromInterval.firstPosition,fromInterval.lastPosition) &&
				 !( ( startA>fromInterval.firstPosition+DISTANCE_THRESHOLD && 
					  (otherAlignment.isLeftMaximal==1 || !Reads.isLeftMaximal(startA,ReadA.id,true) || otherAlignment.lowQualityStart)
					) ||
				    ( endA<fromInterval.lastPosition-DISTANCE_THRESHOLD && 
					  (otherAlignment.isRightMaximal==1 || !Reads.isRightMaximal(endA,ReadA.id,true) || otherAlignment.lowQualityEnd)
				    )
				  )
			   ) {
			   count++;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("hasCompatibleAlignments> compatible alignment found: "+otherAlignment.toStringBoundaries());
			   if (count>=countThreshold) return true;
			}
		}
		for (i=fromAlignment+1; i<=ReadA.lastSortedAlignment; i++) {
			otherAlignment=ReadA.sortedAlignments[i];
			id=otherAlignment.id;
			startB=Alignments.alignments[id][5]; endB=Alignments.alignments[id][6];
			if (Alignments.alignments[id][1]!=readB || startB>fromStartB+SEARCH_THRESHOLD) break;
			if (Alignments.alignments[id][2]!=orientation) continue;
			if (otherAlignment.isLeftMaximalB==1 && otherAlignment.isRightMaximalB==1) continue;
			if (otherAlignment.impliedByPeriodicSubstring!=null || otherAlignment.periodicSubstringInterval!=null) continue;
			startA=Alignments.alignments[id][3]; endA=Alignments.alignments[id][4];
			if ( Intervals.areApproximatelyIdentical(startB,endB,fromStartB,fromEndB) &&
				 Intervals.isApproximatelyContained(startA,endA,fromInterval.firstPosition,fromInterval.lastPosition) &&
				 !( ( startA>fromInterval.firstPosition+DISTANCE_THRESHOLD && 
					  (otherAlignment.isLeftMaximal==1 || !Reads.isLeftMaximal(startA,ReadA.id,true) || otherAlignment.lowQualityStart)
					) ||
				    ( endA<fromInterval.lastPosition-DISTANCE_THRESHOLD && 
					  (otherAlignment.isRightMaximal==1 || !Reads.isRightMaximal(endA,ReadA.id,true) || otherAlignment.lowQualityEnd)
				    )
				  )
			   ) {
			   count++;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("hasCompatibleAlignments> compatible alignment found: "+otherAlignment.toStringBoundaries());
			   if (count>=countThreshold) return true;
		    }
		}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("hasCompatibleAlignments> final count: "+count+" countThreshold="+countThreshold);
		return false;
	}
	
	
	/**
	 * Same as $hasCompatibleAlignments()$, but sets the $impliedByPeriodicSubstring$ 
	 * field of $fromAlignment$, and of all other marked alignments, to an artificial 
	 * periodic substring that points to $fromInterval$. We do this, rather than setting 
	 * to $fromInterval$ the $periodicSubstringInterval$ field of alignments, for 
	 * compatibility with the other steps of the pipeline. Artificial substrings are 
	 * added to $longPeriodSubstrings$ in the interval 
	 * $[lastLongPeriodSubstring+1..lastLongPeriodSubstring_artificial]$
	 */
	private static final void assignCompatibleAlignments(int fromAlignment, PeriodicSubstringInterval fromInterval, PeriodicSubstringInterval tmpInterval) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int SEARCH_THRESHOLD = IO.quantum<<1;
		int i;
		int id, readB, orientation, fromStartB, fromEndB, startA, endA, startB, endB;
		final int fromAlignmentID = ReadA.sortedAlignments[fromAlignment].id;
		Alignment otherAlignment;
		PeriodicSubstring artificialSubstring;
		
		artificialSubstring=getArtificialSubstring_longPeriod(fromInterval);
		if (artificialSubstring==null) artificialSubstring=addArtificialSubstring_longPeriod(fromInterval);
		ReadA.sortedAlignments[fromAlignment].impliedByPeriodicSubstring=artificialSubstring;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("assignCompatibleAlignments> compatible alignment assigned: "+ReadA.sortedAlignments[fromAlignment].toStringBoundaries());
		readB=Alignments.alignments[fromAlignmentID][1];
		orientation=Alignments.alignments[fromAlignmentID][2];
		fromStartB=Alignments.alignments[fromAlignmentID][5];
		fromEndB=Alignments.alignments[fromAlignmentID][6];
		for (i=fromAlignment-1; i>=0; i--) {
			otherAlignment=ReadA.sortedAlignments[i];
			id=otherAlignment.id;
			startB=Alignments.alignments[id][5]; endB=Alignments.alignments[id][6];
			if (Alignments.alignments[id][1]!=readB || startB<fromStartB-SEARCH_THRESHOLD) break;
			if (Alignments.alignments[id][2]!=orientation) continue;
			if (otherAlignment.isLeftMaximalB==1 && otherAlignment.isRightMaximalB==1) continue;
			if (otherAlignment.impliedByPeriodicSubstring!=null || otherAlignment.periodicSubstringInterval!=null) continue;
			startA=Alignments.alignments[id][3]; endA=Alignments.alignments[id][4];
			if ( Intervals.areApproximatelyIdentical(startB,endB,fromStartB,fromEndB) &&
				 Intervals.isApproximatelyContained(startA,endA,fromInterval.firstPosition,fromInterval.lastPosition) &&
				 !( ( startA>fromInterval.firstPosition+DISTANCE_THRESHOLD && 
					  (otherAlignment.isLeftMaximal==1 || !Reads.isLeftMaximal(startA,ReadA.id,true) || otherAlignment.lowQualityStart)
					) ||
				    ( endA<fromInterval.lastPosition-DISTANCE_THRESHOLD && 
					  (otherAlignment.isRightMaximal==1 || !Reads.isRightMaximal(endA,ReadA.id,true) || otherAlignment.lowQualityEnd)
				    )
				  )
			   ) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("assignCompatibleAlignments> compatible alignment assigned: "+otherAlignment.toStringBoundaries());
				   otherAlignment.impliedByPeriodicSubstring=artificialSubstring;
				   fromInterval.nImpliedAlignments++;
				   fromInterval.longestAlignment=Math.max(fromInterval.longestAlignment,otherAlignment.getALength());
			}
		}
		for (i=fromAlignment+1; i<=ReadA.lastSortedAlignment; i++) {
			otherAlignment=ReadA.sortedAlignments[i];
			id=otherAlignment.id;
			startB=Alignments.alignments[id][5];
			endB=Alignments.alignments[id][6];
			if (Alignments.alignments[id][1]!=readB || startB>fromStartB+SEARCH_THRESHOLD) break;
			if (Alignments.alignments[id][2]!=orientation) continue;
			if (otherAlignment.isLeftMaximalB==1 && otherAlignment.isRightMaximalB==1) continue;
			if (otherAlignment.impliedByPeriodicSubstring!=null || otherAlignment.periodicSubstringInterval!=null) continue;
			startA=Alignments.alignments[id][3]; endA=Alignments.alignments[id][4];
			if ( Intervals.areApproximatelyIdentical(startB,endB,fromStartB,fromEndB) &&
				 Intervals.isApproximatelyContained(startA,endA,fromInterval.firstPosition,fromInterval.lastPosition) &&
				 !( ( startA>fromInterval.firstPosition+DISTANCE_THRESHOLD && 
					  (otherAlignment.isLeftMaximal==1 || !Reads.isLeftMaximal(startA,ReadA.id,true) || otherAlignment.lowQualityStart)
					) ||
				    ( endA<fromInterval.lastPosition-DISTANCE_THRESHOLD && 
					  (otherAlignment.isRightMaximal==1 || !Reads.isRightMaximal(endA,ReadA.id,true) || otherAlignment.lowQualityEnd)
				    )
				  )
			   ) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("assignCompatibleAlignments> compatible alignment assigned: "+otherAlignment.toStringBoundaries());
				   otherAlignment.impliedByPeriodicSubstring=artificialSubstring;
				   fromInterval.nImpliedAlignments++;
				   fromInterval.longestAlignment=Math.max(fromInterval.longestAlignment,otherAlignment.getALength());
			}
		}
	}
	
	
	/**
	 * Clones all elements in $longPeriodSubstrings$ to the end of $substrings$.
	 */
	public static final void mergePeriodicSubstrings() {
		if (lastLongPeriodSubstring==-1) return;
		final int NEW_LENGTH = lastSubstring+lastLongPeriodSubstring+2;
		int i;
		
		if (substrings.length<NEW_LENGTH) {
			PeriodicSubstring[] tmp = new PeriodicSubstring[NEW_LENGTH];
			System.arraycopy(substrings,0,tmp,0,lastSubstring+1);
			substrings=tmp;
		}
		for (i=0; i<=lastLongPeriodSubstring; i++) {
			if (substrings[lastSubstring+1+i]==null) {
				substrings[lastSubstring+1+i] = new PeriodicSubstring();
				substrings[lastSubstring+1+i].allocateMemory(MAX_DISTINCT_SHIFTS);
			}
			substrings[lastSubstring+1+i].clone(longPeriodSubstrings[i]);
		}
		lastSubstring=NEW_LENGTH-1;
		previousLongPeriodOrder=PeriodicSubstring.order_longPeriod;
	}
	
	
	/**
	 * Removes all long-period substrings from $substrings$, but keeps the new objects 
	 * allocated by $mergePeriodicSubstrings()$ in $substrings$. 
	 * Restores $PeriodicSubstring.order_longPeriod$ to the value it had before
	 * $mergePeriodicSubstrings()$ was called.
	 */
	public static final void undo_mergePeriodicSubstrings() {
		int i, j;
		PeriodicSubstring tmp;
		
		j=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i].hasLongPeriod) continue;
			j++;
			if (j!=i) {
				tmp=substrings[i];
				substrings[i]=substrings[j];
				substrings[j]=tmp;
			}
		}
		lastSubstring=j;
		for (i=lastSubstring+1; i<substrings.length; i++) {
			if (substrings[i]!=null) substrings[i].hasLongPeriod=false;
		}
		PeriodicSubstring.order_longPeriod=previousLongPeriodOrder;
	}
	
	
	/**
	 * Sets to NULL the $impliedByPeriodicSubstring$ pointer of an alignment, iff the 
	 * readA interval of the alignment straddles the periodic substring pointed to, but
	 * does not have a large-enough intersection with it.
	 */
	private static final void discardAlignmentsAfterTrimming(int firstAlignment, int lastAlignment, boolean updateLongestAlignment) {
		final double INTERSECTION_RATIO = 0.75;  // Arbitrary
		int i;
		int startA, endA, substringStartA, substringEndA, length;
		
		for (i=firstAlignment; i<=lastAlignment; i++) {
			startA=ReadA.sortedAlignments[i].startA();
			endA=ReadA.sortedAlignments[i].endA();
			length=endA-startA+1;
			if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) {
				substringStartA=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.minStartA;
				substringEndA=ReadA.sortedAlignments[i].impliedByPeriodicSubstring.maxEndA;
				if (updateLongestAlignment) ReadA.sortedAlignments[i].impliedByPeriodicSubstring.longestAlignment=0;
				if (Intervals.straddles(startA,endA,substringStartA,substringEndA) && Intervals.intersectionLength(startA,endA,substringStartA,substringEndA)<length*INTERSECTION_RATIO) ReadA.sortedAlignments[i].impliedByPeriodicSubstring=null;
			}
		}
		if (!updateLongestAlignment) return;
		for (i=firstAlignment; i<=lastAlignment; i++) {
			startA=ReadA.sortedAlignments[i].startA();
			endA=ReadA.sortedAlignments[i].endA();
			length=endA-startA+1;
			if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) ReadA.sortedAlignments[i].impliedByPeriodicSubstring.longestAlignment=Math.max(length,ReadA.sortedAlignments[i].impliedByPeriodicSubstring.longestAlignment);
		}
	}
	
	
	/**
	 * @return a shortest periodic substring that is merged to $interval$ and contains 
	 * $[startA..endA]$; otherwise, the only periodic substring that is merged to 
	 * $interval$; otherwise, NULL.
	 */
	public static final PeriodicSubstring findPeriodicSubstring(int startA, int endA, PeriodicSubstringInterval interval) {
		boolean shortOrLong;
		int i;
		int last, length, minLength, substringID;
		PeriodicSubstring minSubstring;
		PeriodicSubstring[] strings;
		
		// Searching the main list
		minSubstring=null; minLength=Math.POSITIVE_INFINITY;
		strings=interval.hasLongPeriod?longPeriodSubstrings:substrings;
		last=interval.hasLongPeriod?lastLongPeriodSubstring:lastSubstring;
		substringID=-1; shortOrLong=true;
		for (i=0; i<=last; i++) {
			if (strings[i].mergedToInterval!=interval) continue;
			if (substringID==-1) substringID=i;
			else substringID=-2;
			if ( Intervals.isApproximatelyContained_lowQuality(startA,endA,strings[i].minStartA,strings[i].maxEndA,ReadA.id) ||
				 Intervals.areApproximatelyIdentical_lowQuality(startA,endA,strings[i].minStartA,strings[i].maxEndA,ReadA.id)
			   ) {
			   length=strings[i].length();
			   if (length<minLength) {
				   minLength=length;
				   minSubstring=strings[i];
			   }
			}
		}
		
		// Searching the other list (since a long-period interval might originate from a
		// short-period substring, and a short-period interval might originate from a
		// long-period substring).
		strings=interval.hasLongPeriod?substrings:longPeriodSubstrings;
		last=interval.hasLongPeriod?lastSubstring:lastLongPeriodSubstring;
		for (i=0; i<=last; i++) {
			if (strings[i].mergedToInterval!=interval) continue;
			if (substringID==-1) {
				substringID=i;
				shortOrLong=false;
			}
			else substringID=-2;
			if ( Intervals.isApproximatelyContained_lowQuality(startA,endA,strings[i].minStartA,strings[i].maxEndA,ReadA.id) ||
				 Intervals.areApproximatelyIdentical_lowQuality(startA,endA,strings[i].minStartA,strings[i].maxEndA,ReadA.id)
			   ) {
			   length=strings[i].length();
			   if (length<minLength) {
				   minLength=length;
				   minSubstring=strings[i];
			   }
			}
		}
		if (minSubstring!=null) return minSubstring;
		
		// Returning the only substring that is merged to $interval$.
		if (substringID>=0) {
			if (shortOrLong) {
				return interval.hasLongPeriod?longPeriodSubstrings[substringID]:substrings[substringID];
			}
			else return interval.hasLongPeriod?substrings[substringID]:longPeriodSubstrings[substringID];
		}
		return null;
	}
	
	
	/**
	 * Sets the boundaries of a periodic interval to the average of all boundaries of 
	 * alignments assigned to it, that are not too far from the boundaries of the 
	 * interval. This is useful, since alignments from many interval types might have been 
	 * assigned to a periodic substring interval throughout the pipeline.
	 *
	 * Remark: the procedure does not update the boundaries if their new value is too 
	 * close to the old or new value of another periodic interval.
	 *
	 * Remark: only B-maximal alignments are used for updating the boundaries.
	 *
	 * Remark: the procedure assumes $longPeriodIntervals$ to be sorted by 
	 * $firstPosition$. At the end of the procedure, $intervals$ and $longPeriodIntervals$
	 * are sorted by $firstPosition$.
	 */
	public static final void fixBoundariesOfPeriodicIntervals() {
		final int DISTANCE_THRESHOLD = IO.quantum<<1;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int WINDOW_LARGE = 3;  // Arbitrary
		final int WINDOW_SMALL = 1;
		final int MIN_COVERAGE = 3*IO.quantum;  // Arbitrary
		final int QUALITY_WINDOW = (Alignments.minAlignmentLength>>1)/Reads.QUALITY_SPACING;  // Arbitrary
		boolean found, foundLeft, foundRight, atLeastOneChanged, atLeastOneChanged_longPeriod, previousSortByID;
		int i, j, k;
		int startA, endA, newValue;
		PeriodicSubstringInterval interval;
		
		// Trimming interval boundaries
		for (i=0; i<=lastInterval; i++) {
			periodTmp[0]=intervals[i].firstPosition; 
			periodTmp[1]=intervals[i].lastPosition;
			Reads.trim(periodTmp,WINDOW_LARGE,true,true);
			if (periodTmp[0]==-1 || periodTmp[1]==-1 || periodTmp[1]-periodTmp[0]+1<Alignments.minAlignmentLength-IDENTITY_THRESHOLD) continue;
			Reads.trim(periodTmp,WINDOW_SMALL,true,true);
			if (periodTmp[0]==-1 || periodTmp[1]==-1 || periodTmp[1]-periodTmp[0]+1<Alignments.minAlignmentLength-IDENTITY_THRESHOLD) continue;
			intervals[i].firstPosition=periodTmp[0];
			intervals[i].lastPosition=periodTmp[1];
		}
		atLeastOneChanged_longPeriod=false;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			periodTmp[0]=longPeriodIntervals[i].firstPosition; 
			periodTmp[1]=longPeriodIntervals[i].lastPosition;
			Reads.trim(periodTmp,WINDOW_LARGE,true,true);
			if (periodTmp[0]==-1 || periodTmp[1]==-1 || periodTmp[1]-periodTmp[0]+1<Alignments.minAlignmentLength-IDENTITY_THRESHOLD) continue;
			Reads.trim(periodTmp,WINDOW_SMALL,true,true);
			if (periodTmp[0]==-1 || periodTmp[1]==-1 || periodTmp[1]-periodTmp[0]+1<Alignments.minAlignmentLength-IDENTITY_THRESHOLD) continue;
			if (longPeriodIntervals[i].firstPosition!=periodTmp[0]) atLeastOneChanged_longPeriod=true;
			longPeriodIntervals[i].firstPosition=periodTmp[0];
			longPeriodIntervals[i].lastPosition=periodTmp[1];
		}
		
		// Ensuring that no interval is significantly larger/smaller than its implied
		// alignments.
		for (i=0; i<=lastInterval; i++) {
			intervals[i].tmpInt2=Math.POSITIVE_INFINITY;  // Used to cumulate a min
			intervals[i].tmpInt4=0;  // Used to cumulate a max
			intervals[i].tmpInt1=0;  // Used to cumulate a left coverage
			intervals[i].tmpInt3=0;  // Used to cumulate a right coverage
		}
		for (i=0; i<=lastLongPeriodInterval; i++) {
			longPeriodIntervals[i].tmpInt2=Math.POSITIVE_INFINITY;
			longPeriodIntervals[i].tmpInt4=0;
			intervals[i].tmpInt1=0;
			intervals[i].tmpInt3=0;
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			interval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if (interval==null) continue;
			startA=ReadA.sortedAlignments[i].startA();
			endA=ReadA.sortedAlignments[i].endA();
			periodTmp[0]=startA; periodTmp[1]=endA;
			Reads.trim(periodTmp,WINDOW_LARGE,true,true);
			if (periodTmp[0]==-1 || periodTmp[1]==-1 || periodTmp[1]-periodTmp[0]+1<Alignments.minAlignmentLength-IDENTITY_THRESHOLD) continue;
			Reads.trim(periodTmp,WINDOW_SMALL,true,true);
			if (periodTmp[0]==-1 || periodTmp[1]==-1 || periodTmp[1]-periodTmp[0]+1<Alignments.minAlignmentLength-IDENTITY_THRESHOLD) continue;
			startA=periodTmp[0]; endA=periodTmp[1];
			interval.tmpInt2=Math.min(interval.tmpInt2,startA);
			interval.tmpInt4=Math.max(interval.tmpInt4,endA);
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			interval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if (interval==null) continue;
			interval.tmpInt1+=Intervals.intersectionLength(interval.tmpInt2,interval.firstPosition-1,ReadA.sortedAlignments[i].startA(),ReadA.sortedAlignments[i].endA());
			interval.tmpInt3+=Intervals.intersectionLength(interval.lastPosition+1,interval.tmpInt4,ReadA.sortedAlignments[i].startA(),ReadA.sortedAlignments[i].endA());
		}
		
		// Copying values from $longPeriodIntervals$ to $intervals$.
		if (atLeastOneChanged_longPeriod && lastLongPeriodInterval>0) {
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			previousSortByID=PeriodicSubstringInterval.sortByID;
			PeriodicSubstringInterval.sortByID=false;
			Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
			PeriodicSubstringInterval.sortByID=previousSortByID;
		}
		for (i=0; i<=lastInterval; i++) {
			if (!intervals[i].hasLongPeriod) continue;
			previousSortByID=PeriodicSubstringInterval.sortByID;
			PeriodicSubstringInterval.sortByID=false;
			j=Arrays.binarySearch(longPeriodIntervals,0,lastLongPeriodInterval+1,intervals[i]);
			PeriodicSubstringInterval.sortByID=previousSortByID;
			if (j<0) {
				System.err.println("fixBoundariesOfPeriodicIntervals> (0) ERROR: long-period interval not found: "+intervals[i]+" order_longPeriod="+PeriodicSubstringInterval.order_longPeriod+" lastLongPeriodInterval="+lastLongPeriodInterval+" order="+PeriodicSubstringInterval.order);
				System.exit(1);
			}
			found=false;
			for (k=j; k>=0; k--) {
				if (longPeriodIntervals[k].equals(intervals[i])) {
					intervals[i].tmpInt2=longPeriodIntervals[k].tmpInt2;
					intervals[i].tmpInt4=longPeriodIntervals[k].tmpInt4;
					intervals[i].tmpInt1=longPeriodIntervals[k].tmpInt1;
					intervals[i].tmpInt3=longPeriodIntervals[k].tmpInt3;
					found=true;
					break;
				}
			}
			if (found) continue;
			for (k=j+1; k<=lastLongPeriodInterval; k++) {
				if (longPeriodIntervals[k].equals(intervals[i])) {
					intervals[i].tmpInt2=longPeriodIntervals[k].tmpInt2;
					intervals[i].tmpInt4=longPeriodIntervals[k].tmpInt4;
					intervals[i].tmpInt1=longPeriodIntervals[k].tmpInt1;
					intervals[i].tmpInt3=longPeriodIntervals[k].tmpInt3;
					found=true;
					break;
				}
			}
			if (!found) {
				System.err.println("fixBoundariesOfPeriodicIntervals> (0.5) ERROR: long-period interval not found: "+intervals[i]+" order_longPeriod="+PeriodicSubstringInterval.order_longPeriod+" lastLongPeriodInterval="+lastLongPeriodInterval);
				System.err.println("intervals:");
				for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x]);
				System.exit(1);
			}
		}
		// Updating boundaries and sorting
		atLeastOneChanged=false;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].tmpInt2!=Math.POSITIVE_INFINITY) {
				Reads.readEnds2qualityEnds(ReadA.id,intervals[i].tmpInt2,intervals[i].firstPosition-1,true,stack);
				if ( intervals[i].firstPosition<intervals[i].tmpInt2-DISTANCE_THRESHOLD ||
					 ( intervals[i].firstPosition>intervals[i].tmpInt2+DISTANCE_THRESHOLD &&
					   intervals[i].tmpInt1/(intervals[i].firstPosition-intervals[i].tmpInt2)>=MIN_COVERAGE &&
					   !Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),stack[0],stack[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,QUALITY_WINDOW,true,-1)
					 )
				   ) {
					intervals[i].firstPosition=intervals[i].tmpInt2;
					atLeastOneChanged=true;
				}
			}
			if (intervals[i].tmpInt4!=0) {
				Reads.readEnds2qualityEnds(ReadA.id,intervals[i].lastPosition,intervals[i].tmpInt4-1,true,stack);
				if ( intervals[i].lastPosition>intervals[i].tmpInt4+DISTANCE_THRESHOLD ||
					 ( intervals[i].lastPosition<intervals[i].tmpInt4-DISTANCE_THRESHOLD &&
					   intervals[i].tmpInt3/(intervals[i].tmpInt4-intervals[i].lastPosition)>=MIN_COVERAGE &&
					   !Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),stack[0],stack[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,QUALITY_WINDOW,true,-1)
					 )
				   ) intervals[i].lastPosition=intervals[i].tmpInt4;
			}
		}
		atLeastOneChanged_longPeriod=false;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].tmpInt2!=Math.POSITIVE_INFINITY) {
				Reads.readEnds2qualityEnds(ReadA.id,longPeriodIntervals[i].tmpInt2,longPeriodIntervals[i].firstPosition-1,true,stack);
				if ( longPeriodIntervals[i].firstPosition<longPeriodIntervals[i].tmpInt2-DISTANCE_THRESHOLD ||
					 ( longPeriodIntervals[i].firstPosition>longPeriodIntervals[i].tmpInt2+DISTANCE_THRESHOLD &&
					   longPeriodIntervals[i].tmpInt1/(longPeriodIntervals[i].firstPosition-longPeriodIntervals[i].tmpInt2)>=MIN_COVERAGE &&
					   !Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),stack[0],stack[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,QUALITY_WINDOW,true,-1)
					 )
				   ) {
					longPeriodIntervals[i].firstPosition=longPeriodIntervals[i].tmpInt2;
					atLeastOneChanged_longPeriod=true;
				}
			}
			if (longPeriodIntervals[i].tmpInt4!=0) {
				Reads.readEnds2qualityEnds(ReadA.id,longPeriodIntervals[i].lastPosition,longPeriodIntervals[i].tmpInt4-1,true,stack);
				if ( longPeriodIntervals[i].lastPosition>longPeriodIntervals[i].tmpInt4+DISTANCE_THRESHOLD ||
					 ( longPeriodIntervals[i].lastPosition<longPeriodIntervals[i].tmpInt4-DISTANCE_THRESHOLD &&
					   longPeriodIntervals[i].tmpInt3/(longPeriodIntervals[i].tmpInt4-longPeriodIntervals[i].lastPosition)>=MIN_COVERAGE &&
					   !Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),stack[0],stack[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,QUALITY_WINDOW,true,-1)
					 )
				   ) longPeriodIntervals[i].lastPosition=longPeriodIntervals[i].tmpInt4;
			}
		}
		if (atLeastOneChanged && lastInterval>0) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			Arrays.sort(intervals,0,lastInterval+1);
		}
		if (atLeastOneChanged_longPeriod && lastLongPeriodInterval>0) {
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
		}
		
		// Computing new boundaries
		for (i=0; i<=lastInterval; i++) {
			intervals[i].nStartA=0; intervals[i].sumStartA=0;
			intervals[i].nEndA=0; intervals[i].sumEndA=0;
		}
		for (i=0; i<=lastLongPeriodInterval; i++) {
			longPeriodIntervals[i].nStartA=0; longPeriodIntervals[i].sumStartA=0;
			longPeriodIntervals[i].nEndA=0; longPeriodIntervals[i].sumEndA=0;
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			interval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if (interval==null) continue;
			startA=ReadA.sortedAlignments[i].startA();
			if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && Math.abs(startA,interval.firstPosition)<=DISTANCE_THRESHOLD) {
				interval.sumStartA+=startA;
				interval.nStartA++;
			}
			endA=ReadA.sortedAlignments[i].endA();
			if (ReadA.sortedAlignments[i].isRightMaximalB==1 && Math.abs(endA,interval.lastPosition)<=DISTANCE_THRESHOLD) {
				interval.sumEndA+=endA;
				interval.nEndA++;
			}
		}
		// Copying the new boundaries from $longPeriodIntervals$ to $intervals$.
		for (i=0; i<=lastInterval; i++) {
			if (!intervals[i].hasLongPeriod) continue;
			previousSortByID=PeriodicSubstringInterval.sortByID;
			PeriodicSubstringInterval.sortByID=false;
			j=Arrays.binarySearch(longPeriodIntervals,0,lastLongPeriodInterval+1,intervals[i]);
			PeriodicSubstringInterval.sortByID=previousSortByID;
			if (j<0) {
				System.err.println("fixBoundariesOfPeriodicIntervals> (1) ERROR: long-period interval not found: "+intervals[i]+" order_longPeriod="+PeriodicSubstringInterval.order_longPeriod+" lastLongPeriodInterval="+lastLongPeriodInterval);
				System.exit(1);
			}
			found=false;
			for (k=j; k>=0; k--) {
				if (longPeriodIntervals[k].equals(intervals[i])) {
					intervals[i].sumStartA=longPeriodIntervals[k].sumStartA;
					intervals[i].nStartA=longPeriodIntervals[k].nStartA;
					intervals[i].sumEndA=longPeriodIntervals[k].sumEndA;
					intervals[i].nEndA=longPeriodIntervals[k].nEndA;
					found=true;
					break;
				}
			}
			if (found) continue;
			for (k=j+1; k<=lastLongPeriodInterval; k++) {
				if (longPeriodIntervals[k].equals(intervals[i])) {
					intervals[i].sumStartA=longPeriodIntervals[k].sumStartA;
					intervals[i].nStartA=longPeriodIntervals[k].nStartA;
					intervals[i].sumEndA=longPeriodIntervals[k].sumEndA;
					intervals[i].nEndA=longPeriodIntervals[k].nEndA;
					found=true;
					break;
				}
			}
			if (!found) {
				System.err.println("fixBoundariesOfPeriodicIntervals> (2) ERROR: long-period interval not found: "+intervals[i]+" order_longPeriod="+PeriodicSubstringInterval.order_longPeriod+" lastLongPeriodInterval="+lastLongPeriodInterval);
				System.err.println("intervals:");
				for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x]);
				System.exit(1);
			}
		}
		// Discarding new boundaries if they are too similar to those of other intervals
		for (i=0; i<=lastInterval; i++) {
			intervals[i].flag1=true;
			intervals[i].flag2=true;
		}
		for (i=0; i<=lastLongPeriodInterval; i++) {
			longPeriodIntervals[i].flag1=true;
			longPeriodIntervals[i].flag2=true;
		}
		for (i=0; i<=lastInterval; i++) {
			startA=intervals[i].firstPosition;
			endA=intervals[i].lastPosition;
			// Left boundary
			foundLeft=false;
			if (intervals[i].nStartA>0) {
				newValue=intervals[i].sumStartA/intervals[i].nStartA;
				for (j=i-1; j>=0; j--) {
					if (intervals[j].lastPosition<=startA) break;
					if ( (intervals[j].nStartA>0 && Math.abs(newValue,intervals[j].sumStartA/intervals[j].nStartA)<=IDENTITY_THRESHOLD) || 
					     Math.abs(newValue,intervals[j].firstPosition)<=IDENTITY_THRESHOLD
					   ) {
						foundLeft=true;
						break;
					}
				}
				if (!foundLeft) {
					for (j=i+1; j<=lastInterval; j++) {
						if (intervals[j].firstPosition>=endA) break;
						if ( (intervals[j].nStartA>0 && Math.abs(newValue,intervals[j].sumStartA/intervals[j].nStartA)<=IDENTITY_THRESHOLD) || 
						     Math.abs(newValue,intervals[j].firstPosition)<=IDENTITY_THRESHOLD
						   ) {
							foundLeft=true;
							break;
						}
					}
				}
				if (foundLeft) intervals[i].flag1=false;
			}
			// Right boundary
			foundRight=false;
			if (intervals[i].nEndA>0) {
				newValue=intervals[i].sumEndA/intervals[i].nEndA;
				for (j=i+1; j<=lastInterval; j++) {
					if (intervals[j].firstPosition>=endA) break;
					if ( (intervals[j].nEndA>0 && Math.abs(newValue,intervals[j].sumEndA/intervals[j].nEndA)<=IDENTITY_THRESHOLD) || 
					     Math.abs(newValue,intervals[j].lastPosition)<=IDENTITY_THRESHOLD
					   ) {
						foundRight=true;
						break;
					}
				}
				if (!foundRight) {
					for (j=i-1; j>=0; j--) {
						if (intervals[j].lastPosition<=startA) continue;
						if ( (intervals[j].nEndA>0 && Math.abs(newValue,intervals[j].sumEndA/intervals[j].nEndA)<=IDENTITY_THRESHOLD) || 
						     Math.abs(newValue,intervals[j].lastPosition)<=IDENTITY_THRESHOLD
						   ) {
							foundRight=true;
							break;
						}
					}
				}
				if (foundRight) intervals[i].flag2=false;
			}
			if (!intervals[i].hasLongPeriod || (!foundLeft && !foundRight)) continue;
			// Updating $longPeriodIntervals$.
			previousSortByID=PeriodicSubstringInterval.sortByID;
			PeriodicSubstringInterval.sortByID=false;
			j=Arrays.binarySearch(longPeriodIntervals,0,lastLongPeriodInterval+1,intervals[i]);
			PeriodicSubstringInterval.sortByID=previousSortByID;
			if (j<0) {
				System.err.println("fixBoundariesOfPeriodicIntervals> (3) ERROR: long-period interval not found: "+intervals[i]+" order_longPeriod="+PeriodicSubstringInterval.order_longPeriod+" lastLongPeriodInterval="+lastLongPeriodInterval);
				System.exit(1);
			}
			found=false;
			for (k=j; k>=0; k--) {
				if (longPeriodIntervals[k].equals(intervals[i])) {
					if (foundLeft) longPeriodIntervals[k].flag1=false;
					if (foundRight) longPeriodIntervals[k].flag2=false;
					found=true;
					break;
				}
			}
			if (found) continue;
			for (k=j+1; k<=lastLongPeriodInterval; k++) {
				if (longPeriodIntervals[k].equals(intervals[i])) {
					if (foundLeft) longPeriodIntervals[k].flag1=false;
					if (foundRight) longPeriodIntervals[k].flag2=false;
					found=true;
					break;
				}
			}
			if (!found) {
				System.err.println("fixBoundariesOfPeriodicIntervals> (4) ERROR: long-period interval not found: "+intervals[i]+" order_longPeriod="+PeriodicSubstringInterval.order_longPeriod+" lastLongPeriodInterval="+lastLongPeriodInterval);
				System.err.println("intervals:");
				for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x]);
				System.exit(1);
			}
		}
		// Updating boundaries and sorting
		atLeastOneChanged=false;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].nStartA!=0 && (intervals[i].flag1 || intervals[i].flag2)) {
				newValue=intervals[i].sumStartA/intervals[i].nStartA;
				if (newValue!=intervals[i].firstPosition) {
					intervals[i].firstPosition=newValue;
					atLeastOneChanged=true;
				}
			}
			if (intervals[i].nEndA!=0 && (intervals[i].flag2 || intervals[i].flag1)) intervals[i].lastPosition=intervals[i].sumEndA/intervals[i].nEndA;
		}
		atLeastOneChanged_longPeriod=false;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].nStartA!=0 && (longPeriodIntervals[i].flag1 || longPeriodIntervals[i].flag2)) {
				newValue=longPeriodIntervals[i].sumStartA/longPeriodIntervals[i].nStartA;
				if (newValue!=longPeriodIntervals[i].firstPosition) {
					longPeriodIntervals[i].firstPosition=newValue;
					atLeastOneChanged_longPeriod=true;
				}
			}
			if (longPeriodIntervals[i].nEndA!=0 && (longPeriodIntervals[i].flag2 || longPeriodIntervals[i].flag1)) longPeriodIntervals[i].lastPosition=longPeriodIntervals[i].sumEndA/longPeriodIntervals[i].nEndA;
		}
		if (atLeastOneChanged && lastInterval>0) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			Arrays.sort(intervals,0,lastInterval+1);
		}
		if (atLeastOneChanged_longPeriod && lastLongPeriodInterval>0) {
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
		}
	}
	
	
	/**
	 * Stores in arrays ${left,right}Splits$ the sorted peaks induced by all {left,right}-
	 * maximal alignments inside maximum ranges of straddling or adjacent short-period 
	 * intervals. All alignments inside such ranges are used, not just those implied by 
	 * periodic substrings. Short-period periodic intervals with uniform period should 
	 * have approx. uniform distribution of maximal events and uniform coverage (except 
	 * near the boundaries); however, in practice, peaks of maximal events and 
	 * discontinuities in coverage can be seen. 
	 * The procedure stores the last element of $*Splits$ in global variable $last*Split$.
	 *
	 * Remark: the procedure flags the peaks that correspond to the beginning/end of an 
	 * interval in arrays ${left,right}SplitFlags$, and the peaks that coincide with 
	 * low-quality regions in arrays ${left,right}SplitLowQuality$.
	 *
	 * Remark: the procedure assumes that a period estimate has already been assigned to
	 * all intervals.
	 *
	 * Remark: the procedure uses also arrays $DenseSubstrings.delta{Left,Right}$ and
	 * $DenseSubstrings.tmpSplits$.
	 *
	 * Remark: at the end of the procedure, global variable $splitsUpdated$ is set to 
	 * true.
	 */
	private static final void splitShortPeriodIntervals() {
		final double CONSTANT_THRESHOLD = 0.13;  // Arbitrary
		final double MIN_HIGH = 0.3;  // Arbitrary leaf density
		final int MAX_PEAK_RADIUS = Alignments.minAlignmentLength>>2;  // Arbitrary
		final int MAX_WIDTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		final int LOW_QUALITY_WINDOW = IO.quantum;  // Arbitrary
		final int MIN_PERIOD_FOR_HISTOGRAM = MIN_PERIOD_LONG>>1;
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_MASS_HIGH = 50*IO.coverage;  // Arbitrary
		boolean isConcatenation;
		int i, j;
		int firstAlignment, startA, endA, from, to;
		int previousLastLeftSplit, previousLastRightSplit, firstIntervalInRange;
		int rangeFirst, rangeLast, periodNumerator, periodDenominator, readLength;
		int distanceThreshold = IO.quantum<<1;
		int minAlignmentsForHistogram;
		final double[] array = Reads.getQualityArray(ReadA.id);
		
		// Ensuring the necessary order in intervals and alignments
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		
		// Iterating over maximal ranges of straddling or adjacent short-period intervals.
		lastLeftSplit=-1; lastRightSplit=-1; previousLastLeftSplit=-1;
		rangeFirst=intervals[0].firstPosition;
		rangeLast=intervals[0].lastPosition;
		periodNumerator=intervals[0].period; periodDenominator=1;
		firstAlignment=0; firstIntervalInRange=0; isConcatenation=false;
		for (i=1; i<=lastInterval; i++) {
			Reads.readEnds2qualityEnds(ReadA.id,rangeLast+1,intervals[i].firstPosition-1,true,periodTmp);
			if (intervals[i].firstPosition<=rangeLast+IDENTITY_THRESHOLD) {
				if ( !( intervals[i].firstPosition>=rangeLast-IDENTITY_THRESHOLD && 
					    (intervals[i].firstPosition<=rangeLast || !Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),periodTmp[0],periodTmp[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,1,true,-1))
					 )
				   ) {
					rangeLast=Math.max(rangeLast,intervals[i].lastPosition);
					if (intervals[i].period>0) {
						periodNumerator+=intervals[i].period;
						periodDenominator++;
					}
					isConcatenation|=splitShortPeriodIntervals_isConcatenation(i,firstIntervalInRange,IDENTITY_THRESHOLD);
					continue;
				}
			}
			startA=rangeFirst; endA=rangeLast;
			// Left-maximal
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> LEFT HISTOGRAM OF ["+startA+".."+endA+"]:");
			ReadA.getDeltaHistogram_periodic(startA,endA,firstAlignment,endA+1,true,DenseSubstrings.deltaLeft);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> ReadA.deltaTmp[0]="+ReadA.deltaTmp[0]);
			if (ReadA.deltaTmp[0]>=(endA-MIN_PERIOD_FOR_HISTOGRAM-startA+1)/MIN_PERIOD_FOR_HISTOGRAM) {
				from=Histograms.firstNonzero(DenseSubstrings.deltaLeft,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> firstNonzero outputs "+from);
				if (from!=-1) {
					to=Histograms.lastNonzero(DenseSubstrings.deltaLeft,0,endA-startA);	
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> lastNonzero outputs "+to);
					previousLastLeftSplit=lastLeftSplit;
					lastLeftSplit=Histograms.getPeaks(DenseSubstrings.deltaLeft,from,to,leftSplits,lastLeftSplit+1,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,DenseSubstrings.tmpSplits,isConcatenation?-1:periodNumerator/periodDenominator,MAX_WIDTH,ReadA.deltaPoints,ReadA.lastDeltaPoint,MIN_MASS_HIGH,Alignments.minAlignmentLength);
					for (j=previousLastLeftSplit+1; j<=lastLeftSplit; j++) leftSplitFlags[j]=false;
					if (lastLeftSplit>previousLastLeftSplit && startA>=leftSplits[previousLastLeftSplit+2]-MAX_PEAK_RADIUS) {
						leftSplits[previousLastLeftSplit+1]=startA;
						leftSplitFlags[previousLastLeftSplit+1]=true;
						leftSplitFlags[previousLastLeftSplit+2]=true;
						leftSplitFlags[previousLastLeftSplit+3]=true;
					}
					else if (previousLastLeftSplit>0 && startA>=leftSplits[previousLastLeftSplit-1]-MAX_PEAK_RADIUS) {
						leftSplits[previousLastLeftSplit-2]=startA;
						leftSplitFlags[previousLastLeftSplit-2]=true;
						leftSplitFlags[previousLastLeftSplit-1]=true;
						leftSplitFlags[previousLastLeftSplit]=true;
					}
				}
			}
			// Right-maximal
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> RIGHT HISTOGRAM OF ["+startA+".."+endA+"]:");
			ReadA.getDeltaHistogram_periodic(startA,endA,firstAlignment,endA+1,false,DenseSubstrings.deltaRight);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> ReadA.deltaTmp[0]="+ReadA.deltaTmp[0]);
			if (ReadA.deltaTmp[0]>=(endA-MIN_PERIOD_FOR_HISTOGRAM-startA+1)/MIN_PERIOD_FOR_HISTOGRAM) {
				from=Histograms.firstNonzero(DenseSubstrings.deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> firstNonzero outputs "+from);
				if (from!=-1) {
					to=Histograms.lastNonzero(DenseSubstrings.deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> lastNonzero outputs "+to);
					previousLastRightSplit=lastRightSplit;
					lastRightSplit=Histograms.getPeaks(DenseSubstrings.deltaRight,from,to,rightSplits,lastRightSplit+1,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,DenseSubstrings.tmpSplits,isConcatenation?-1:periodNumerator/periodDenominator,MAX_WIDTH,ReadA.deltaPoints,ReadA.lastDeltaPoint,MIN_MASS_HIGH,Alignments.minAlignmentLength);
					for (j=previousLastRightSplit+1; j<=lastRightSplit; j++) rightSplitFlags[j]=false;
					if (lastRightSplit>previousLastRightSplit && endA<=rightSplits[lastRightSplit-1]+MAX_PEAK_RADIUS) {
						rightSplits[lastRightSplit]=endA;
						rightSplitFlags[lastRightSplit]=true;
						rightSplitFlags[lastRightSplit-1]=true;
						rightSplitFlags[lastRightSplit-2]=true;
					}
					else if (previousLastRightSplit>0 && endA<=rightSplits[previousLastRightSplit-1]+MAX_PEAK_RADIUS) {
						rightSplits[previousLastRightSplit]=endA;
						rightSplitFlags[previousLastRightSplit]=true;
						rightSplitFlags[previousLastRightSplit-1]=true;
						rightSplitFlags[previousLastRightSplit-2]=true;
					}
				}
			}	
			// Next range
			rangeFirst=intervals[i].firstPosition;
			rangeLast=intervals[i].lastPosition;
			periodNumerator=intervals[i].period; periodDenominator=1;
			firstAlignment=ReadA.deltaTmp[1];
			isConcatenation=false; firstIntervalInRange=i;
		}
		startA=rangeFirst; endA=rangeLast;
		// Last left-maximal
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> LEFT HISTOGRAM OF ["+startA+".."+endA+"]:");
		ReadA.getDeltaHistogram_periodic(startA,endA,firstAlignment,endA+1,true,DenseSubstrings.deltaLeft);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> ReadA.deltaTmp[0]="+ReadA.deltaTmp[0]);
		if (ReadA.deltaTmp[0]>=(endA-MIN_PERIOD_FOR_HISTOGRAM-startA+1)/MIN_PERIOD_FOR_HISTOGRAM) {
			from=Histograms.firstNonzero(DenseSubstrings.deltaLeft,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> firstNonzero outputs "+from);
			if (from!=-1) {
				to=Histograms.lastNonzero(DenseSubstrings.deltaLeft,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> lastNonzero outputs "+to);
				previousLastLeftSplit=lastLeftSplit;
				lastLeftSplit=Histograms.getPeaks(DenseSubstrings.deltaLeft,from,to,leftSplits,lastLeftSplit+1,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,DenseSubstrings.tmpSplits,isConcatenation?-1:periodNumerator/periodDenominator,MAX_WIDTH,ReadA.deltaPoints,ReadA.lastDeltaPoint,MIN_MASS_HIGH,Alignments.minAlignmentLength);
				for (j=previousLastLeftSplit+1; j<=lastLeftSplit; j++) leftSplitFlags[j]=false;
				if (lastLeftSplit>previousLastLeftSplit && startA>=leftSplits[previousLastLeftSplit+2]-MAX_PEAK_RADIUS) {
					leftSplits[previousLastLeftSplit+1]=startA;
					leftSplitFlags[previousLastLeftSplit+1]=true;
					leftSplitFlags[previousLastLeftSplit+2]=true;
					leftSplitFlags[previousLastLeftSplit+3]=true;
				}
				else if (previousLastLeftSplit>0 && startA>=leftSplits[previousLastLeftSplit-1]-MAX_PEAK_RADIUS) {
					leftSplits[previousLastLeftSplit-2]=startA;
					leftSplitFlags[previousLastLeftSplit-2]=true;
					leftSplitFlags[previousLastLeftSplit-1]=true;
					leftSplitFlags[previousLastLeftSplit]=true;
				}
			}
		}
		// Last right-maximal
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> RIGHT HISTOGRAM OF ["+startA+".."+endA+"]:");
		ReadA.getDeltaHistogram_periodic(startA,endA,firstAlignment,endA+1,false,DenseSubstrings.deltaRight);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> ReadA.deltaTmp[0]="+ReadA.deltaTmp[0]);
		if (ReadA.deltaTmp[0]>=(endA-MIN_PERIOD_FOR_HISTOGRAM-startA+1)/MIN_PERIOD_FOR_HISTOGRAM) {
			from=Histograms.firstNonzero(DenseSubstrings.deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> firstNonzero outputs "+from);
			if (from!=-1) {
				to=Histograms.lastNonzero(DenseSubstrings.deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitShortPeriodIntervals> lastNonzero outputs "+to);
				previousLastRightSplit=lastRightSplit;
				lastRightSplit=Histograms.getPeaks(DenseSubstrings.deltaRight,from,to,rightSplits,lastRightSplit+1,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,DenseSubstrings.tmpSplits,isConcatenation?-1:periodNumerator/periodDenominator,MAX_WIDTH,ReadA.deltaPoints,ReadA.lastDeltaPoint,MIN_MASS_HIGH,Alignments.minAlignmentLength);
				for (j=previousLastRightSplit+1; j<=lastRightSplit; j++) rightSplitFlags[j]=false;
				if (lastRightSplit>=0 && endA<=rightSplits[lastRightSplit-1]+MAX_PEAK_RADIUS) {
					rightSplits[lastRightSplit]=endA;
					rightSplitFlags[lastRightSplit]=true;
					rightSplitFlags[lastRightSplit-1]=true;
					rightSplitFlags[lastRightSplit-2]=true;
				}
				else if (previousLastRightSplit>0 && endA<=rightSplits[previousLastRightSplit-1]+MAX_PEAK_RADIUS) {
					rightSplits[previousLastRightSplit]=endA;
					rightSplitFlags[previousLastRightSplit]=true;
					rightSplitFlags[previousLastRightSplit-1]=true;
					rightSplitFlags[previousLastRightSplit-2]=true;
				}
			}
		}
		splitsUpdated=true;
		
		// Assigning qualities to splits
		readLength=Reads.getReadLength(ReadA.id);
		Math.set(leftSplitLowQuality,lastLeftSplit,false);
		for (i=1; i<=lastLeftSplit; i+=3) {
			from=Math.max(0,Math.min(leftSplits[i-1],leftSplits[i]-LOW_QUALITY_WINDOW));
			to=Math.min(readLength-1,Math.max(leftSplits[i+1],leftSplits[i]+LOW_QUALITY_WINDOW));
			Reads.readEnds2qualityEnds(ReadA.id,from,to,true,periodTmp);
			j=leftSplits[i]/Reads.QUALITY_SPACING;
			if ( array[j]>=Reads.MIN_RANDOM_QUALITY_SCORE || 
				 array[Math.max(j-1,0)]>=Reads.MIN_RANDOM_QUALITY_SCORE || 
			     array[Math.min(j+1,Reads.getQualityArrayLength(ReadA.id)-1)]>=Reads.MIN_RANDOM_QUALITY_SCORE || 
				 Histograms.hasSubstringWithAverage(array,periodTmp[0],periodTmp[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,1,true,-1) ||
			     Reads.hasLowQuality(ReadA.id,from,to,true)
			   ) {
				leftSplitLowQuality[i-1]=true;
				leftSplitLowQuality[i]=true;
				leftSplitLowQuality[i+1]=true;
			}
		}
		Math.set(rightSplitLowQuality,lastRightSplit,false);
		for (i=1; i<=lastRightSplit; i+=3) {
			from=Math.max(0,Math.min(rightSplits[i-1],rightSplits[i]-LOW_QUALITY_WINDOW));
			to=Math.min(readLength-1,Math.max(rightSplits[i+1],rightSplits[i]+LOW_QUALITY_WINDOW));
			Reads.readEnds2qualityEnds(ReadA.id,from,to,true,periodTmp);
			j=rightSplits[i]/Reads.QUALITY_SPACING;
			if ( array[j]>=Reads.MIN_RANDOM_QUALITY_SCORE || 
				 array[Math.max(j-1,0)]>=Reads.MIN_RANDOM_QUALITY_SCORE || 
			     array[Math.min(j+1,Reads.getQualityArrayLength(ReadA.id)-1)]>=Reads.MIN_RANDOM_QUALITY_SCORE || 
				 Histograms.hasSubstringWithAverage(array,periodTmp[0],periodTmp[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,1,true,-1) ||
			     Reads.hasLowQuality(ReadA.id,from,to,true)
			   ) {
				rightSplitLowQuality[i-1]=true;
				rightSplitLowQuality[i]=true;
				rightSplitLowQuality[i+1]=true;
			}
		}		
		
if (IO.SHOW_STD_ERR_PRIME) {		
	System.err.print("splitShortPeriodIntervals> at the end, deltaLeft: ");
	for (int x=0; x<=lastLeftSplit; x++) System.err.print(leftSplits[x]+",");
	System.err.print(" flags: ");
	for (int x=0; x<=lastLeftSplit; x++) System.err.print(leftSplitFlags[x]+",");
	System.err.print(" low quality? ");
	for (int x=0; x<=lastLeftSplit; x++) System.err.print(leftSplitLowQuality[x]+",");
	System.err.println();
	System.err.print("splitShortPeriodIntervals> at the end, deltaRight: ");
	for (int x=0; x<=lastRightSplit; x++) System.err.print(rightSplits[x]+",");
	System.err.print(" flags: ");
	for (int x=0; x<=lastRightSplit; x++) System.err.print(rightSplitFlags[x]+",");
	System.err.print(" low quality? ");
	for (int x=0; x<=lastRightSplit; x++) System.err.print(rightSplitLowQuality[x]+",");
	System.err.println();
}	
	}
	
	
	/**
	 * Remark: the procedure assumes $intervals$ to be sorted by first position.
	 *
	 * @return TRUE iff there is an interval in $intervals[firstIntervalInRange..$ that 
	 * ends where $intervals[interval]$ starts.
	 */
	private static final boolean splitShortPeriodIntervals_isConcatenation(int interval, int firstIntervalInRange, int identityThreshold) {
		final int first = intervals[interval].firstPosition;

		for (int i=interval-1; i>=firstIntervalInRange; i--) {
			if (Math.abs(intervals[i].lastPosition,first)<=identityThreshold) return true;			
		}
		return false;
	}
	
	
	/**
	 * Discards elements of $inputIntervals$ that are contained in a maximal range of 
	 * straddling or adjacent short-period intervals, and such that their first position
	 * is not inside a peak of maximal events in $leftSplits$, or their last position is 
	 * not inside a peak in $rightSplits$. For every discarded interval, the procedure 
	 * merges its shifts to those of a non-discarded shortest short-period interval that 
	 * approximately contains it, if any. Finally, the procedure marks as 
	 * $isContainedInShort$ all input intervals that are approx. contained in a maximal 
	 * range of straddling or adjacent short-period intervals.
	 *
	 * Remark: an interval is discarded also if it is aligned just to peaks of low 
	 * quality. This is essentially equivalent to the merging of nonmaximal contained 
	 * substrings done in $removeSubstringsForFactoring()$, but applied to maximal ranges.
	 *
	 * Remark: other approaches for filtering contained intervals are not effective in 
	 * practice. For example, we tried to keep, for every interval, the number of 
	 * substrings that have been merged into it (a proxy for its frequency in the genome,
	 * althout they are not necessarily equal). After the merge, an interval that is 
	 * contained in another interval can be discarded iff its number of substrings is too 
	 * small (the threshold being decided by fitting a DET on the numbers of substrings of 
	 * all contained intervals). Note also that B-maximality for contained periodic
	 * substrings is not necessarily a signal for splitting.
	 *
	 * Remark: the procedure assumes $intervals$ to contain just short-period intervals.
	 *
	 * Remark: the $*Splits$ arrays are assumed to have been already filled by
	 * $splitShortPeriodIntervals()$.
	 *
	 * Remark: the procedure works correctly even when $inputIntervals==intervals$.
	 * 
	 * Remark: the procedure assumes $intervals$ and $inputIntervals$ to be sorted by 
	 * first position.
	 *
	 * @param removeIdentical considers for removal also elements of $inputIntervals$ that
	 * are approximately identical to a maximal range;
	 * @return the last element in $inputIntervals$ after filtering.
	 */
	private static final int filterIntervalsWithPeaks(PeriodicSubstringInterval[] inputIntervals, int lastInputInterval, PeriodicSubstring[] inputSubstrings, int lastInputSubstring, boolean removeIdentical, boolean updateContainerShifts, boolean updateContainerPeriod) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int LOW_COUNT = 1;
		final int MIN_LENGTH_LOWCOUNT = IO.quantum<<1;
		final double RANGE_RATIO = 0.75;  // Arbitrary
		boolean isConcatenation, adjacentCompartments;
		int i, j, k;
		int pos, rangeFirst, rangeLast, firstKForNextI, length, minLength, nDiscarded;
		int rangeKeptFirst, rangeKeptLast, rangeSurface, firstIntervalInRange;
		PeriodicSubstringInterval artificialIntervalContained = new PeriodicSubstringInterval();
		PeriodicSubstringInterval artificialIntervalIdentical = new PeriodicSubstringInterval();
		PeriodicSubstringInterval tmpInterval, minInterval;
		
		// Marking input intervals contained in the union of short-period intervals.
		for (i=0; i<=lastInputInterval; i++) {
			inputIntervals[i].representative=null;
			inputIntervals[i].isContainedInShort=true;
		}
		rangeFirst=intervals[0].firstPosition;
		rangeLast=intervals[0].lastPosition;
		isConcatenation=false; firstIntervalInRange=0;
		j=0;
		for (i=1; i<=lastInterval; i++) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> "+(inputIntervals==intervals)+" adding interval ["+intervals[i].firstPosition+".."+intervals[i].lastPosition+"] to interval ["+rangeFirst+".."+rangeLast+"]");
			Reads.readEnds2qualityEnds(ReadA.id,rangeLast+1,intervals[i].firstPosition-1,true,periodTmp);
			if (intervals[i].firstPosition<=rangeLast+IDENTITY_THRESHOLD) {
				if ( !( inputIntervals==intervals && 
					    intervals[i].firstPosition>=rangeLast-IDENTITY_THRESHOLD && 
					    (intervals[i].firstPosition<=rangeLast || !Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),periodTmp[0],periodTmp[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,1,true,-1))
					 )
				   ) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> interval ["+rangeFirst+".."+Math.max(rangeLast,intervals[i].lastPosition)+"] built");
					rangeLast=Math.max(rangeLast,intervals[i].lastPosition);
					isConcatenation|=splitShortPeriodIntervals_isConcatenation(i,firstIntervalInRange,IDENTITY_THRESHOLD);
					continue;
				}
			}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> marking. current range: ["+rangeFirst+".."+rangeLast+"] j="+j+" isConcatenation="+isConcatenation);
			j=filterIntervalsWithPeaks_impl(j,inputIntervals,lastInputInterval,rangeFirst,rangeLast,removeIdentical,isConcatenation,DISTANCE_THRESHOLD,artificialIntervalContained,artificialIntervalIdentical);
			// Next range
			rangeFirst=intervals[i].firstPosition;
			rangeLast=intervals[i].lastPosition;
			isConcatenation=false; firstIntervalInRange=i;
		}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> marking. current range: ["+rangeFirst+".."+rangeLast+"] j="+j+" isConcatenation="+isConcatenation);
		filterIntervalsWithPeaks_impl(j,inputIntervals,lastInputInterval,rangeFirst,rangeLast,removeIdentical,isConcatenation,DISTANCE_THRESHOLD,artificialIntervalContained,artificialIntervalIdentical);


if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 1 non-marked:");
	for (int x=0; x<=lastInputInterval; x++) {
		if (inputIntervals[x].representative!=artificialIntervalContained && inputIntervals[x].representative!=artificialIntervalIdentical) System.err.println(inputIntervals[x]);
	}
}	


if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> input intervals at the very beginning:");
	System.err.println("artificialIntervalIdentical="+artificialIntervalIdentical);
	System.err.println("artificialIntervalContained="+artificialIntervalContained);
	for (int x=0; x<=lastInputInterval; x++) System.err.println(inputIntervals[x]+" ==REPRESENTATIVE==> "+inputIntervals[x].representative);
}

		// Discarding input intervals
		for (i=0; i<=lastInputInterval; i++) inputIntervals[i].isContainedInShort=false;  // Resetting $isContainedInShort$.
		nDiscarded=0;
		rangeFirst=inputIntervals[0].firstPosition;
		rangeLast=inputIntervals[0].lastPosition;
		rangeKeptFirst=-1; rangeKeptLast=-1; rangeSurface=0;
		firstIntervalInRange=0;
		for (i=0; i<=lastInputInterval; i++) {
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> considering for discard interval "+inputIntervals[i]);
	System.err.println("        current intervals:");
	for (int x=0; x<=lastInputInterval; x++) System.err.println("        "+inputIntervals[x]);
}
			Reads.readEnds2qualityEnds(ReadA.id,rangeLast+1,intervals[i].firstPosition-1,true,periodTmp);
			if (inputIntervals==intervals) {
				if ( inputIntervals[i].firstPosition<=rangeLast+IO.quantum && 
					 intervals[i].firstPosition>=rangeLast-IO.quantum && 
					 (intervals[i].firstPosition<=rangeLast || !Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),periodTmp[0],periodTmp[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,1,true,-1))
				   ) adjacentCompartments=true;
				else adjacentCompartments=false;
				if (inputIntervals[i].firstPosition<=rangeLast+IO.quantum && !adjacentCompartments) {
					rangeLast=Math.max(rangeLast,inputIntervals[i].lastPosition);
				}
				else {
					if (rangeKeptFirst>=0 && rangeKeptLast>=0) rangeSurface+=rangeKeptLast-rangeKeptFirst+1;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 3 rangeSurface="+rangeSurface);
					if (rangeSurface<(rangeLast-rangeFirst+1)*RANGE_RATIO) nDiscarded+=filterIntervalsWithPeaks_smallSurface(firstIntervalInRange,i-1,rangeFirst,rangeLast,DISTANCE_THRESHOLD);
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 4");
	System.err.println("        current intervals:");
	for (int x=0; x<=lastInputInterval; x++) System.err.println("        "+inputIntervals[x]);
}
					rangeSurface=0; firstIntervalInRange=i;
					rangeFirst=inputIntervals[i].firstPosition;
					rangeLast=inputIntervals[i].lastPosition;
					rangeKeptFirst=-1; rangeKeptLast=-1;
				}
			}
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 5  interval currently: "+intervals[i]);
	System.err.println("        current intervals:");
	for (int x=0; x<=lastInputInterval; x++) System.err.println("        "+inputIntervals[x]);
}
			
			// Handling current interval
			inputIntervals[i].discarded=false;
			if (inputIntervals[i].representative==null) {
				if (rangeKeptFirst==-1 && rangeKeptLast==-1) {
					rangeKeptFirst=inputIntervals[i].firstPosition;
					rangeKeptLast=inputIntervals[i].lastPosition;
				}
				else {
					if (inputIntervals[i].firstPosition<=rangeKeptLast+IO.quantum) {
						rangeKeptLast=Math.max(rangeKeptLast,inputIntervals[i].lastPosition);
					}
					else {
						rangeSurface+=rangeKeptLast-rangeKeptFirst+1;
						rangeKeptFirst=inputIntervals[i].firstPosition;
						rangeKeptLast=inputIntervals[i].lastPosition;
					}
				}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 5.1");
				continue;
			}
			inputIntervals[i].isContainedInShort=true;
			if (removeIdentical && inputIntervals[i].representative==artificialIntervalIdentical) {
				// Discarding intervals that coincide with a short-period range.
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 5.2");
				inputIntervals[i].discarded=true;
				nDiscarded++;
				inputIntervals[i].representative=null;
				continue;
			}
			pos=inputIntervals[i].firstPosition;
			j=position2split(pos,leftSplits,lastLeftSplit,DISTANCE_THRESHOLD);
			if (j==-1) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 6  interval currently: "+intervals[i]);
				inputIntervals[i].discarded=true;
				nDiscarded++;
			}
			else {
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 7 j="+j+" pos="+pos);
	System.err.println("        current intervals:");
	for (int x=0; x<=lastInputInterval; x++) System.err.println("        "+inputIntervals[x]);
}
				pos=inputIntervals[i].lastPosition;
				k=position2split(pos,rightSplits,lastRightSplit,DISTANCE_THRESHOLD);
				if (k==-1) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 8");
					inputIntervals[i].discarded=true;
					nDiscarded++;
				}
				else {
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 9 j="+j+" flag="+leftSplitFlags[3*j]+" k="+k+" flag="+rightSplitFlags[3*k]);
	System.err.println("filterIntervalsWithPeaks> 9 "+(leftSplitLowQuality[3*j] || leftSplitFlags[3*j]));
	System.err.println("filterIntervalsWithPeaks> 9 "+(rightSplitLowQuality[3*k] || rightSplitFlags[3*k]));
	System.err.println("filterIntervalsWithPeaks> 9 "+(leftSplitFlags[3*j] && rightSplitFlags[3*k]));
	System.err.println("filterIntervalsWithPeaks> 9 "+(inputIntervals==intervals?filterIntervalsWithPeaks_isCovered(i,firstIntervalInRange,DISTANCE_THRESHOLD):true));
}
					if ( removeIdentical && leftSplitFlags[3*j] && rightSplitFlags[3*k] &&
						 (inputIntervals==intervals?inputIntervals[i].isContainedInShort:true)
					   ) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 9.1");
						// Discarding intervals that coincide with a short-period range.
						inputIntervals[i].discarded=true;
						nDiscarded++;
					}
					else if ( (leftSplitLowQuality[3*j] || leftSplitFlags[3*j]) &&
						      (rightSplitLowQuality[3*k] || rightSplitFlags[3*k]) &&
							  !(leftSplitFlags[3*j] && rightSplitFlags[3*k]) &&
							  (inputIntervals==intervals?filterIntervalsWithPeaks_isCovered(i,firstIntervalInRange,DISTANCE_THRESHOLD):true)
					        ) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 9.2");
						// Discarding intervals with both endpoints on a low-quality peak.
						inputIntervals[i].discarded=true;
						nDiscarded++;
					}
					else {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 9.3");
						if (rangeKeptFirst==-1 && rangeKeptLast==-1) {
	if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 10");
							rangeKeptFirst=inputIntervals[i].firstPosition;
							rangeKeptLast=inputIntervals[i].lastPosition;
						}
						else {
	if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 11");
							if (inputIntervals[i].firstPosition<=rangeKeptLast+IO.quantum) {
								rangeKeptLast=Math.max(rangeKeptLast,inputIntervals[i].lastPosition);
							}
							else {
								rangeSurface+=rangeKeptLast-rangeKeptFirst+1;
								rangeKeptFirst=inputIntervals[i].firstPosition;
								rangeKeptLast=inputIntervals[i].lastPosition;
							}
						}
					}
				}
			}
			inputIntervals[i].representative=null;
		}
		// Last range
		if (rangeKeptFirst>=0 && rangeKeptLast>=0) rangeSurface+=rangeKeptLast-rangeKeptFirst+1;
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 12 rangeSurface="+rangeSurface+" rangeFirst="+rangeFirst+" rangeLast="+rangeLast);
	for (int x=0; x<=lastInputInterval; x++) System.err.println(inputIntervals[x]);
}
		if (inputIntervals==intervals && rangeSurface<(rangeLast-rangeFirst+1)*RANGE_RATIO) {
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 13");
	for (int x=0; x<=lastInputInterval; x++) System.err.println(inputIntervals[x]);
}
			nDiscarded+=filterIntervalsWithPeaks_smallSurface(firstIntervalInRange,lastInputInterval,rangeFirst,rangeLast,DISTANCE_THRESHOLD);
		}
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 14 nDiscarded="+nDiscarded+". input intervals:");
	for (int x=0; x<=lastInputInterval; x++) System.err.println(inputIntervals[x]);
}

		
		// Enforcing at most one input interval per pair of splits
		nDiscarded+=discardIntervalsOnSameSplits(inputIntervals,lastInputInterval,DISTANCE_THRESHOLD);
		if (nDiscarded==0) return lastInputInterval;
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> 15 nDiscarded="+nDiscarded+". input intervals:");
	for (int x=0; x<=lastInputInterval; x++) System.err.println(inputIntervals[x]);
}
		
		// Merging the shifts of every discarded input interval to those of its shortest
		// short-period container, if any, and updating the container's period.
		for (i=0; i<=lastInterval; i++) intervals[i].surface=0;  // Using $surface$ as temporary space.
		firstKForNextI=0; k=0;
		for (i=0; i<=lastInputInterval; i++) {
			if (!inputIntervals[i].discarded) continue;
			minInterval=null; minLength=Math.POSITIVE_INFINITY;
			if (firstKForNextI!=-1) k=firstKForNextI;
			firstKForNextI=-1;
			while (k<=lastInterval) {
				if (intervals[k].discarded) {
					k++;
					continue;
				}
				if (intervals[k].firstPosition>=inputIntervals[i].lastPosition) break;
				if (firstKForNextI==-1 && i<lastInputInterval && intervals[k].lastPosition>=inputIntervals[i+1].firstPosition) firstKForNextI=k;
				if ( ( Intervals.isApproximatelyContained(inputIntervals[i].firstPosition,inputIntervals[i].lastPosition,intervals[k].firstPosition,intervals[k].lastPosition) ||
					   Intervals.areApproximatelyIdentical_lowQuality(inputIntervals[i].firstPosition,inputIntervals[i].lastPosition,intervals[k].firstPosition,intervals[k].lastPosition,ReadA.id) ||
					   Intervals.areApproximatelyIdentical_lowCoverage(inputIntervals[i].firstPosition,inputIntervals[i].lastPosition,intervals[k].firstPosition,intervals[k].lastPosition)
					 ) &&
					 !( (inputIntervals[i].firstPosition>intervals[k].firstPosition+DISTANCE_THRESHOLD && !Reads.isLeftMaximal(inputIntervals[i].firstPosition,ReadA.id,true)) &&
					    (inputIntervals[i].lastPosition<intervals[k].lastPosition-DISTANCE_THRESHOLD && !Reads.isRightMaximal(inputIntervals[i].lastPosition,ReadA.id,true))
					 )
				   ) {
					length=intervals[k].length();
					if (length<minLength) {
					   minLength=length;
					   minInterval=intervals[k];
					}
				}
				k++;
			}
			if (minInterval!=null) {
				if (updateContainerShifts) {
					lastShift=Points.merge(inputIntervals[i].shifts,inputIntervals[i].lastShift,minInterval.shifts,minInterval.lastShift,shifts);
					minInterval.cloneShifts(shifts,0,lastShift,MAX_DISTINCT_SHIFTS,Points.tmpPoints);
				}
				if (updateContainerPeriod && minInterval.period==0) minInterval.surface=Math.max(minInterval.surface,inputIntervals[i].period);
				minInterval.simpleMerge(inputIntervals[i]);
			}
			inputIntervals[i].representative=minInterval;
		}
		if (updateContainerPeriod) {
			for (i=0; i<=lastInterval; i++) {
				if (intervals[i].period==0 && intervals[i].surface!=0) intervals[i].period=intervals[i].surface;
			}
		}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 16");
		
		// Updating pointers from input substrings
		for (i=0; i<=lastInputSubstring; i++) {
			tmpInterval=inputSubstrings[i].mergedToInterval;
            if (tmpInterval==null || !tmpInterval.discarded) continue;
            inputSubstrings[i].mergedToInterval=tmpInterval.representative;   
		}
		if (inputSubstrings==substrings) {
			// Updating also artificial substrings
			updateArtificialSubstringsPointers(true);
			// Updating also short-period substrings in $longPeriodSubstrings$.
			for (i=0; i<=lastLongPeriodSubstring; i++) {
				tmpInterval=longPeriodSubstrings[i].mergedToInterval;
	            if (tmpInterval==null || tmpInterval.hasLongPeriod || !tmpInterval.discarded) continue;
	            longPeriodSubstrings[i].mergedToInterval=tmpInterval.representative;   
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("filterIntervalsWithPeaks> UPDATED POINTER FROM LONG-PERIOD SUBSTRING: "+longPeriodSubstrings[i].hashCode());
	System.err.println("filterIntervalsWithPeaks> OLD INTERVAL POINTER: "+tmpInterval.hashCode()+" :: "+tmpInterval);
	System.err.println("filterIntervalsWithPeaks> NEW INTERVAL POINTER: "+(longPeriodSubstrings[i].mergedToInterval==null?"null":longPeriodSubstrings[i].mergedToInterval.hashCode())+" :: "+longPeriodSubstrings[i].mergedToInterval);
}
			}
		}
		for (i=0; i<=lastInputInterval; i++) inputIntervals[i].representative=null;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("filterIntervalsWithPeaks> 17");
		
		// Compacting input intervals
		j=-1;
		for (i=0; i<=lastInputInterval; i++) {
			if (inputIntervals[i].discarded) continue;
			j++;
			if (j!=i) {
				tmpInterval=inputIntervals[i];
				inputIntervals[i]=inputIntervals[j];
				inputIntervals[j]=tmpInterval;
			}
		}
		return j;
	}
	
	
	/**
	 * The procedure is called when the non-discarded short-period intervals in a maximal
	 * range $[rangeFirst..rangeLast]$ have too small a total surface. The procedure 
	 * finds an interval that is approximately identical to the max range, and makes sure
	 * it is not discarded (if no such interval is found, the procedure resets an existing
	 * interval to it). All non-discarded intervals are kept.
	 *
	 * @param {from,to}Interval compact range in $intervals$;
	 * @return quantity to be added to $nDiscarded$ after the procedure completes.
	 */
	private static final int filterIntervalsWithPeaks_smallSurface(int fromInterval, int toInterval, int rangeFirst, int rangeLast, int distanceThreshold) {
		boolean isLeftMaximal, isRightMaximal, hasLongPeriod, equalsOnePeriod, mergedToDense, firstPositionChanged;
		boolean closeToStart, closeToEnd;
		int i;
		int found, out;
		
		found=-1; isLeftMaximal=false; isRightMaximal=false; out=0;
		hasLongPeriod=true; equalsOnePeriod=true; mergedToDense=false; firstPositionChanged=false;
		for (i=fromInterval; i<=toInterval; i++) {
			closeToStart=Math.abs(intervals[i].firstPosition,rangeFirst)<=distanceThreshold;
			closeToEnd=Math.abs(intervals[i].lastPosition,rangeLast)<=distanceThreshold;
			if (closeToStart && closeToEnd) {
				if (found==-1) {
					found=i;
					if (intervals[i].discarded) {
						intervals[i].discarded=false;
						out--;
					}
				}
				else {
					if (!intervals[i].discarded) {
						intervals[i].discarded=true;
						out++;
					} 
				}
				if (!intervals[i].hasLongPeriod) hasLongPeriod=false;
				if (!intervals[i].equalsOnePeriod) equalsOnePeriod=false;
				mergedToDense|=intervals[i].mergedToDense;
				isLeftMaximal|=intervals[i].isLeftMaximal;
				isRightMaximal|=intervals[i].isRightMaximal;
			}
			else {
				if (intervals[i].discarded) {
					if (!intervals[i].hasLongPeriod) hasLongPeriod=false;
					if (!intervals[i].equalsOnePeriod) equalsOnePeriod=false;
					mergedToDense|=intervals[i].mergedToDense;
					if (closeToStart) isLeftMaximal|=intervals[i].isLeftMaximal;
					if (closeToEnd) isRightMaximal|=intervals[i].isRightMaximal;	
				}
			}
		}
		if (found==-1) {
			for (i=fromInterval; i<=toInterval; i++) {
				if (intervals[i].discarded) {
					// This must happen at least once.
					found=i;
					break;
				}
			}
			intervals[found].discarded=false; out--;
			if (intervals[found].firstPosition!=rangeFirst) firstPositionChanged=true;
			intervals[found].firstPosition=rangeFirst;
			intervals[found].sumStartA=rangeFirst;
			intervals[found].nStartA=1;
			intervals[found].lastPosition=rangeLast;
			intervals[found].sumEndA=rangeLast;
			intervals[found].nEndA=1;
			intervals[found].isLeftMaximal=isLeftMaximal;
			intervals[found].isRightMaximal=isRightMaximal;
		}
		intervals[found].hasLongPeriod=hasLongPeriod;
		intervals[found].equalsOnePeriod=equalsOnePeriod;
		intervals[found].mergedToDense=mergedToDense;
		if (firstPositionChanged) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			Arrays.sort(intervals,fromInterval,toInterval+1);
		}
		return out;
	}
	
	
	/**
	 * Remark: the procedure assumes $intervals$ to be sorted by first position.
	 * 
	 * @return FALSE iff there is a position of $intervals[interval]$ that is: 
	 * (1) not covered by any other interval in $intervals[fromInterval..]$ which is not 
	 * approx. contained in or identical to $intervals[interval]$; (2) more than 
	 * $identityThreshold$ away from the beginning/end of $intervals[interval]$.
	 */
	private static final boolean filterIntervalsWithPeaks_isCovered(int interval, int fromInterval, int identityThreshold) {
		int i;
		int hole, boundary;
		final int intervalFirst = intervals[interval].firstPosition;
		final int intervalLast = intervals[interval].lastPosition;
		
		hole=intervalFirst+identityThreshold; boundary=-1; i=fromInterval;
		while (hole<=intervalLast-identityThreshold) {
			if (i>lastInterval || intervals[i].firstPosition>intervalLast-identityThreshold) return boundary>=intervalLast-identityThreshold;
			if ( i==interval || 
				 intervals[i].lastPosition<hole ||
				 Intervals.isApproximatelyContained(intervals[i].firstPosition,intervals[i].lastPosition,intervalFirst,intervalLast) ||
			     Intervals.areApproximatelyIdentical(intervals[i].firstPosition,intervals[i].lastPosition,intervalFirst,intervalLast)
			   ) {
				i++;
				continue;
			}
			if (intervals[i].firstPosition>hole) {
				if (boundary==-1) return false;
				hole=boundary+1; boundary=-1;
				continue;
			}
			boundary=Math.max(boundary,intervals[i].lastPosition);
			i++;
		}
		return true;
	}
	
	
	/**
	 * Out of all the non-discarded intervals in $inputIntervals$ that share the same left 
	 * and right split in $*Splits$, the procedure keeps just a longest one, discards the 
	 * rest, and resets the boundaries of the longest one to those of the union.
	 * If $inputIntervals!=intervals$, the procedure discards also elements of 
	 * $inputIntervals$ that use the same splits as an element of $intervals$.
	 *
	 * Remark: the procedure assumes $inputIntervals$ and $intervals$ to be sorted by 
	 * first position.
	 *
	 * Remark: the procedure uses the $surface$ field of an interval as temporary space.
	 *
	 * @param distanceThreshold used to assign splits with $position2split()$;
	 * @return the number of discarded intervals.
	 */
	private static final int discardIntervalsOnSameSplits(PeriodicSubstringInterval[] inputIntervals, int lastInputInterval, int distanceThreshold) {
		boolean firstPositionChanged;
		int i, j;
		int lastJ, leftSplit, rightSplit, lastPosition, representative, out;
		int maxLength, minStart, maxEnd, firstJForNextI;
		if (lastLeftSplit<0 && lastRightSplit<0) return 0;
		
		// $inputIntervals$ -> $inputIntervals$
		for (i=0; i<=lastInputInterval; i++) {
			inputIntervals[i].surface=-1;  // Used as temporary space
			inputIntervals[i].tmpInt1=-1;
			inputIntervals[i].tmpInt2=-1;
		}
		out=0;
		for (i=0; i<=lastInputInterval; i++) {
			if (inputIntervals[i].discarded || inputIntervals[i].surface!=-1) continue;
			leftSplit=position2split(inputIntervals[i].firstPosition,leftSplits,lastLeftSplit,distanceThreshold);
			rightSplit=position2split(inputIntervals[i].lastPosition,rightSplits,lastRightSplit,distanceThreshold);
			if (leftSplit==-1 || rightSplit==-1) continue;
			inputIntervals[i].surface=i; 
			maxLength=inputIntervals[i].length(); 
			minStart=inputIntervals[i].firstPosition;
			maxEnd=inputIntervals[i].lastPosition;
			lastJ=i; lastPosition=inputIntervals[i].lastPosition;
			for (j=i+1; j<=lastInputInterval; j++) {
				if (inputIntervals[j].firstPosition>=lastPosition) break;
				if (inputIntervals[j].discarded) continue;
				if ( position2split(inputIntervals[j].firstPosition,leftSplits,lastLeftSplit,distanceThreshold)==leftSplit &&
					 position2split(inputIntervals[j].lastPosition,rightSplits,lastRightSplit,distanceThreshold)==rightSplit
				   ) {
				    inputIntervals[j].surface=i;
					maxLength=Math.max(maxLength,inputIntervals[j].length());
					minStart=Math.min(minStart,inputIntervals[j].firstPosition);
					maxEnd=Math.max(maxEnd,inputIntervals[j].lastPosition);
					lastJ=j;
				}
			}
			representative=-1;
			for (j=i; j<=lastJ; j++) {
				if (inputIntervals[j].surface==i && inputIntervals[j].length()==maxLength) {
					representative=j;
					break;
				}
			}
			inputIntervals[representative].tmpInt1=minStart;
			inputIntervals[representative].tmpInt2=maxEnd;
			for (j=i; j<=lastJ; j++) {
				if (inputIntervals[j].surface==i && j!=representative) {
					inputIntervals[j].discarded=true;
					out++;
				}
			}
		}
		firstPositionChanged=false;
		for (i=0; i<=lastInputInterval; i++) {
			if (inputIntervals[i].discarded || inputIntervals[i].tmpInt1==-1 || inputIntervals[i].tmpInt2==-1) continue;
			if (inputIntervals[i].tmpInt1!=inputIntervals[i].firstPosition) firstPositionChanged=true;
			inputIntervals[i].firstPosition=inputIntervals[i].tmpInt1;
			inputIntervals[i].lastPosition=inputIntervals[i].tmpInt2;
		}
		if (firstPositionChanged) {
			if (inputIntervals==longPeriodIntervals) PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			else PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			Arrays.sort(inputIntervals,0,lastInputInterval+1);
		}
		
		// $inputIntervals$ -> $intervals$
		if (inputIntervals!=intervals && lastInterval>=0) {
			for (i=0; i<=lastInputInterval; i++) inputIntervals[i].surface=0;
			i=0; 
			leftSplit=position2split(inputIntervals[i].firstPosition,leftSplits,lastLeftSplit,distanceThreshold);
			rightSplit=position2split(inputIntervals[i].lastPosition,rightSplits,lastRightSplit,distanceThreshold);
			j=0; firstJForNextI=-1;
			while (i<=lastInputInterval) {
				if ( inputIntervals[i].discarded || inputIntervals[i].surface==1 || leftSplit==-1 || rightSplit==-1 ||
				     j>lastInterval || intervals[j].firstPosition>=inputIntervals[i].lastPosition
				   ) {
					if (inputIntervals[i].surface==1) {
						inputIntervals[i].discarded=true;
						out++;
					}
					i++;
					if (i<=lastInputInterval) {
						leftSplit=position2split(inputIntervals[i].firstPosition,leftSplits,lastLeftSplit,distanceThreshold);
						rightSplit=position2split(inputIntervals[i].lastPosition,rightSplits,lastRightSplit,distanceThreshold);
					}
					if (firstJForNextI!=-1) j=firstJForNextI;
					firstJForNextI=-1;
					continue;
				}
				if (intervals[j].lastPosition<=inputIntervals[i].firstPosition) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && i<lastInputInterval && intervals[j].lastPosition>=inputIntervals[i+1].firstPosition) firstJForNextI=j;
				if ( position2split(intervals[j].firstPosition,leftSplits,lastLeftSplit,distanceThreshold)==leftSplit && 
					 position2split(intervals[j].lastPosition,rightSplits,lastRightSplit,distanceThreshold)==rightSplit
				   ) inputIntervals[i].surface=1;
				j++;
			}
		}
		
		return out;
	}
	
	
	/**
	 * [THIS PROCEDURE IS NOT CURRENTLY USED]
	 * Sets $isContainedInShort=true$ for all elements of $intervals$ that might not be 
	 * essential for preserving the surface of a maximal range. Given a maximal 
	 * range of short-period intervals, the procedure finds every sub-range of length at 
	 * least $minLength$ and covered by at most $lowCount>=1$ short-period intervals of 
	 * the range, and it sets $isContainedInShort=false$ for all the intervals covering 
	 * the range.
	 *
	 * Remark: the procedure assumes $intervals$ to contain just short-period intervals.
	 *
	 * @param tmpStart,tmpEnd temporary space, assumed to be at least as long as the max.
	 * number of intervals in a short-period range.
	 */
	private static final void markEssentialIntervals(int lowCount, int minLength, int[] tmpStart, int[] tmpEnd, PeriodicSubstringInterval tmpInterval) {
		int i, j;
		int rangeFirst, rangeLast, rangeFirstInterval, lowRangeFirst;
		int lastStart, lastEnd, nextStart, nextEnd, count;
		
		for (i=0; i<=lastInterval; i++) intervals[i].isContainedInShort=true;
		rangeFirstInterval=0;
		rangeFirst=intervals[0].firstPosition;
		rangeLast=intervals[0].lastPosition;
		for (i=1; i<=lastInterval; i++) {
			if (intervals[i].firstPosition<=rangeLast+IO.quantum) {
				if (intervals[i].lastPosition>rangeLast) rangeLast=intervals[i].lastPosition;
				continue;
			}
			
			// Processing current range
			lastStart=-1;
			for (j=rangeFirstInterval; j<i; j++) tmpStart[++lastStart]=intervals[j].firstPosition;
			lastEnd=-1;
			for (j=rangeFirstInterval; j<i; j++) tmpEnd[++lastEnd]=intervals[j].lastPosition;
			if (lastEnd>0) Arrays.sort(tmpEnd,0,lastEnd+1);
			nextStart=1; count=1;
			while (nextStart<=lastStart && tmpStart[nextStart]==tmpStart[0]) {
				nextStart++;
				count++;
			}
			nextEnd=0; lowRangeFirst=count<=lowCount?tmpStart[0]:-1;
			while (nextEnd<=lastEnd) {
				if (nextStart>lastStart) {
					count--;
					if (count==lowCount) lowRangeFirst=tmpEnd[nextEnd]+1;
					nextEnd++;
					continue;
				}
				if (tmpStart[nextStart]==tmpEnd[nextEnd]) {
					nextStart++;
					nextEnd++;
					continue;
				}
				if (tmpEnd[nextEnd]<tmpStart[nextStart]) {
					count--;
					if (count==lowCount) lowRangeFirst=tmpEnd[nextEnd]+1;
					nextEnd++;
					continue;
				}
				else {
					count++;
					if (count==lowCount+1 && tmpStart[nextStart]>=lowRangeFirst+minLength) markEssentialIntervals_impl(lowRangeFirst,tmpStart[nextStart]-1,tmpInterval);
					nextStart++;
				}
			}
			if (count<=lowCount && tmpEnd[lastEnd]>=lowRangeFirst+minLength-1) markEssentialIntervals_impl(lowRangeFirst,tmpEnd[lastEnd],tmpInterval);
			
			// Next range
			rangeFirstInterval=i;
			rangeFirst=intervals[i].firstPosition;
			rangeLast=intervals[i].lastPosition;
		}
		
		// Last range
		lastStart=-1;
		for (j=rangeFirstInterval; j<i; j++) tmpStart[++lastStart]=intervals[j].firstPosition;
		lastEnd=-1;
		for (j=rangeFirstInterval; j<i; j++) tmpEnd[++lastEnd]=intervals[j].lastPosition;
		if (lastEnd>0) Arrays.sort(tmpEnd,0,lastEnd+1);
		nextStart=1; count=1;
		while (nextStart<=lastStart && tmpStart[nextStart]==tmpStart[0]) {
			nextStart++;
			count++;
		}
		nextEnd=0; lowRangeFirst=count<=lowCount?tmpStart[0]:-1;
		while (nextEnd<=lastEnd) {
			if (nextStart>lastStart) {
				count--;
				if (count==lowCount) lowRangeFirst=tmpEnd[nextEnd]+1;
				nextEnd++;
				continue;
			}
			if (tmpStart[nextStart]==tmpEnd[nextEnd]) {
				nextStart++;
				nextEnd++;
				continue;
			}
			if (tmpEnd[nextEnd]<tmpStart[nextStart]) {
				count--;
				if (count==lowCount) lowRangeFirst=tmpEnd[nextEnd]+1;
				nextEnd++;
				continue;
			}
			else {
				count++;
				if (count==lowCount+1 && tmpStart[nextStart]>=lowRangeFirst+minLength) markEssentialIntervals_impl(lowRangeFirst,tmpStart[nextStart]-1,tmpInterval);
				nextStart++;
			}
		}
		if (count<=lowCount && tmpEnd[lastEnd]>=lowRangeFirst+minLength-1) markEssentialIntervals_impl(lowRangeFirst,tmpEnd[lastEnd],tmpInterval);
	}
	
	
	/**
	 * Sets field $isContainedInShort$ to FALSE for all elements of $intervals$ that 
	 * strictly contain $[lowRangeFirst..lowRangeLast]$.
	 */
	private static final void markEssentialIntervals_impl(int lowRangeFirst, int lowRangeLast, PeriodicSubstringInterval tmpInterval) {
		boolean previousSortByID;
		int i, j;
		int first, last;
		
		tmpInterval.firstPosition=lowRangeFirst;
		previousSortByID=PeriodicSubstringInterval.sortByID;
		PeriodicSubstringInterval.sortByID=false;
		i=Arrays.binarySearch(intervals,0,lastInterval+1,tmpInterval);
		PeriodicSubstringInterval.sortByID=previousSortByID;
		if (i<0) i=-i-1;
		j=i-1;
		while (j>=0) {
			last=intervals[j].lastPosition;
			if (last<=lowRangeFirst) {
				j--;
				continue;
			}
			first=intervals[j].firstPosition;
			if (Intervals.isContained(lowRangeFirst,lowRangeLast,first,last)) intervals[j].isContainedInShort=false;
			j--;
		}
		j=i;
		while (j<=lastInterval) {
			first=intervals[j].firstPosition;
			if (first>=lowRangeLast) break;
			last=intervals[j].lastPosition;
			if (Intervals.isContained(lowRangeFirst,lowRangeLast,first,last)) intervals[j].isContainedInShort=false;
			j++;
		}
	}
	
	
	/**
	 * Remark: if $inputIntervals==intervals$, the procedure marks only input intervals 
	 * for which $isContainedInShort=true$.
	 *
	 * @param artificialInterval{Contained,Identical} used to distinguish between 
	 * intervals that are identical to the range or just contained in it;
	 * @param isConcatenation TRUE iff $[rangeFirst..rangeLast]$ is the concatenation of
	 * adjacent but distinct short-period intervals; a long-period interval that is
	 * identical to a maximal short-period range is marked only if $isConcatenation$ is
	 * false, since otherwise the concatenation could be itself the unit of a long-period 
	 * periodic substring;
	 * @return the first value of $fromInputInterval$ to be used in the following range.
	 */
	private static final int filterIntervalsWithPeaks_impl(int fromInputInterval, PeriodicSubstringInterval[] inputIntervals, int lastInputInterval, int rangeFirst, int rangeLast, boolean removeIdentical, boolean isConcatenation, int distanceThreshold, PeriodicSubstringInterval artificialIntervalContained, PeriodicSubstringInterval artificialIntervalIdentical) {
		final double INTERSECTION_THRESHOLD = 0.9;
		final int RANGE_LENGTH = rangeLast-rangeFirst+1;
		final int PERIOD_THRESHOLD = Alignments.minAlignmentLength<<1;
		final int IDENTITY_THRESHOLD = IO.quantum;
		int j, firstJForNext, leftSplit, rightSplit;
		
		j=fromInputInterval; firstJForNext=-1;
		while (j<=lastInputInterval) {
			if (inputIntervals[j].lastPosition<=rangeFirst) {
				j++;
				continue;
			}
			if (inputIntervals[j].firstPosition>=rangeLast) break;
			if (firstJForNext==-1 && inputIntervals[j].lastPosition>rangeLast) firstJForNext=j;
			leftSplit=position2split(inputIntervals[j].firstPosition,leftSplits,lastLeftSplit,IDENTITY_THRESHOLD);
			rightSplit=position2split(inputIntervals[j].lastPosition,rightSplits,lastRightSplit,IDENTITY_THRESHOLD);
			if ( Intervals.areApproximatelyIdentical_lowQuality(inputIntervals[j].firstPosition,inputIntervals[j].lastPosition,rangeFirst,rangeLast,ReadA.id) ||
			     Intervals.areApproximatelyIdentical_lowCoverage(inputIntervals[j].firstPosition,inputIntervals[j].lastPosition,rangeFirst,rangeLast) ||
				 ( Intervals.isApproximatelyContained(inputIntervals[j].firstPosition,inputIntervals[j].lastPosition,rangeFirst,rangeLast) &&
				   Intervals.intersectionLength(inputIntervals[j].firstPosition,inputIntervals[j].lastPosition,rangeFirst,rangeLast)>=RANGE_LENGTH*INTERSECTION_THRESHOLD &&
				   (leftSplit==-1 || rightSplit==-1)
				 )
			   ) {
				// Identity
				if ( removeIdentical && inputIntervals==intervals?
				     inputIntervals[j].isContainedInShort:
				     ( !isConcatenation || 
					   ( inputIntervals[j].period>0 && inputIntervals[j].period<PERIOD_THRESHOLD )
					 )
				   ) {
					inputIntervals[j].representative=artificialIntervalIdentical;
				}
			}
			if (Intervals.isApproximatelyContained(inputIntervals[j].firstPosition,inputIntervals[j].lastPosition,rangeFirst,rangeLast)) {
				// Containment
				if ( inputIntervals==intervals?
				     inputIntervals[j].isContainedInShort:
				     ( !spansRangeBoundary(leftSplit,rightSplit) || 
					   ( inputIntervals[j].period>0 && inputIntervals[j].period<PERIOD_THRESHOLD )
					 )
				   ) {
					inputIntervals[j].representative=artificialIntervalContained;
				}
			}
			j++;
		}
		return firstJForNext==-1?j:firstJForNext;
	}
	
	
	/**
	 * @return TRUE iff $leftSplit$ and $rightSplit$ are valid splits, and if there is a 
	 * range boundary between the two.
	 */
	private static final boolean spansRangeBoundary(int leftSplit, int rightSplit) {
		if (leftSplit==-1 || rightSplit==-1) return false;
		int i;
		final int FIRST_RIGHT = rightSplits[3*rightSplit];
		final int LAST_LEFT = leftSplits[3*leftSplit+2];
		
		for (i=3*(leftSplit+1); i<=lastLeftSplit; i+=3) {
			if (leftSplits[i+2]>=FIRST_RIGHT) break;
			if (leftSplitFlags[i]) return true;
		}
		for (i=3*(rightSplit-1); i>=0; i-=3) {
			if (rightSplits[i]<=LAST_LEFT) break;
			if (rightSplitFlags[i]) return true;
		}
		return false;
	}
	
	
	/**
	 * Remark: the procedure assumes $intervals$ to be sorted by $firstPosition$.
	 * 
	 * @param tmpInterval temporary space;
	 * @return TRUE iff $position$ is inside an interval in $intervals$, and at distance 
	 * at least $distanceThreshold$ from each end.
	 */
	public static final boolean inPeriodicSubstringInterval(int position, int distanceThreshold, PeriodicSubstringInterval tmpInterval) {
		int i;
		if (lastInterval==-1) return false;
		
		tmpInterval.firstPosition=position;
		i=Arrays.binarySearch(intervals,0,lastInterval+1,tmpInterval);
		if (i<0) i=-i-1;
		i--;
		while (i>=0) {
			if (position>intervals[i].firstPosition+distanceThreshold && position<intervals[i].lastPosition-distanceThreshold) return true;
			i--;
		}
		return false;
	}
	
	
	/**
	 * Assume that an alignment is contained in a maximal short-period range, its 
	 * boundaries fall inside two peaks in such range, and there is just one periodic
	 * interval P (with long or short period) whose boundaries fall inside the same peaks.
	 * The procedure marks such alignment as implied by P, and by a (possibly artificial) 
	 * periodic substring that points to P (i.e. the procedure sets both the 
	 * $impliedByPeriodicSubstring$ and the $periodicSubstringInterval$ field of the
	 * alignment).
	 *
	 * Remark: this procedure is useful, since it avoids the later creation of alignment 
	 * intervals that are very similar to periodic intervals, but that fail the approx.
	 * identical test.
	 *
	 * Remark: the $*Splits$ arrays are assumed to have been already filled by
	 * $splitShortPeriodIntervals()$.
	 */
	private static final void markImpliedAlignments_peaks() {
		final int THRESHOLD_SMALL = IO.quantum;
		final int THRESHOLD_LARGE = (IO.quantum)<<1;
		
		int i, j;
		int firstJForNextI, pos, interval, rangeFirst, rangeLast, rangeFirstInterval;
		int leftSplit, rightSplit;
		PeriodicSubstring pSubstring;
		
		// Ensuring the necessary order in $ReadA.sortedAlignments$ and $intervals$.
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
		if (lastLongPeriodInterval>0) Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
		if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		
		// Marking
		for (i=0; i<=lastInterval; i++) {
			if (!intervals[i].hasLongPeriod) break;
		}
		if (i>lastInterval) return;
		rangeFirstInterval=i;
		rangeFirst=intervals[i].firstPosition;
		rangeLast=intervals[i].lastPosition;
		j=0; firstJForNextI=-1;
		for (++i; i<=lastInterval; i++) {
			if (intervals[i].hasLongPeriod) continue;
			if (intervals[i].firstPosition<=rangeLast+IO.quantum) {
				rangeLast=Math.max(rangeLast,intervals[i].lastPosition);
				continue;
			}
			// Current range
			if (firstJForNextI!=-1) j=firstJForNextI;
			firstJForNextI=-1;
			while (j<=ReadA.lastSortedAlignment) {
				pos=ReadA.sortedAlignments[j].startA();
				if (pos>=rangeLast-THRESHOLD_LARGE && firstJForNextI==-1) firstJForNextI=j;
				if (pos>=rangeLast) break;
				leftSplit=position2split(pos,leftSplits,lastLeftSplit,THRESHOLD_LARGE);
				if (leftSplit==-1 || (leftSplitFlags[3*leftSplit] && pos<leftSplits[3*leftSplit+1]-THRESHOLD_SMALL)) {
					j++;
					continue;
				}
				pos=ReadA.sortedAlignments[j].endA();
				if (pos<=rangeFirst) {
					j++;
					continue;
				}
				if (pos>rangeLast+THRESHOLD_LARGE) {
					j++;
					continue;
				}
				rightSplit=position2split(pos,rightSplits,lastRightSplit,THRESHOLD_LARGE);
				if (rightSplit==-1 || (rightSplitFlags[3*rightSplit] && pos>rightSplits[3*rightSplit+1]+THRESHOLD_SMALL)) {
					j++;
					continue;
				}
				interval=markImpliedAlignments_peaks_getInterval(rangeFirstInterval,i-1,rangeFirst,rangeLast,leftSplit,rightSplit,THRESHOLD_LARGE);
				if (interval>=0) {
					pSubstring=markImpliedAlignments_peaks_interval2substring(intervals[interval],leftSplit,rightSplit);
					ReadA.sortedAlignments[j].impliedByPeriodicSubstring=pSubstring;
					ReadA.sortedAlignments[j].periodicSubstringInterval=pSubstring.mergedToInterval;
					if (pSubstring.mergedToInterval!=null) {
						pSubstring.mergedToInterval.nImpliedAlignments++;
						pSubstring.mergedToInterval.longestAlignment=Math.max(pSubstring.mergedToInterval.longestAlignment,ReadA.sortedAlignments[j].getALength());
					}
				}
				j++;
			}
			// Next range
			rangeFirstInterval=i;
			rangeFirst=intervals[i].firstPosition;
			rangeLast=intervals[i].lastPosition;
		}
		// Last range
		if (firstJForNextI!=-1) j=firstJForNextI;
		while (j<=ReadA.lastSortedAlignment) {
			pos=ReadA.sortedAlignments[j].startA();
			if (pos>=rangeLast) break;
			leftSplit=position2split(pos,leftSplits,lastLeftSplit,THRESHOLD_LARGE);
			if (leftSplit==-1 || (leftSplitFlags[3*leftSplit] && pos<leftSplits[3*leftSplit+1]-THRESHOLD_SMALL)) {
				j++;
				continue;
			}
			pos=ReadA.sortedAlignments[j].endA();
			if (pos<=rangeFirst || pos>rangeLast+THRESHOLD_LARGE) {
				j++;
				continue;
			}
			rightSplit=position2split(pos,rightSplits,lastRightSplit,THRESHOLD_LARGE);
			if (rightSplit==-1 || (rightSplitFlags[3*rightSplit] && pos>rightSplits[3*rightSplit+1]+THRESHOLD_SMALL)) {
				j++;
				continue;
			}
			interval=markImpliedAlignments_peaks_getInterval(rangeFirstInterval,i-1,rangeFirst,rangeLast,leftSplit,rightSplit,THRESHOLD_LARGE);
			if (interval>=0) {
				pSubstring=markImpliedAlignments_peaks_interval2substring(intervals[interval],leftSplit,rightSplit);
				ReadA.sortedAlignments[j].impliedByPeriodicSubstring=pSubstring;
				ReadA.sortedAlignments[j].periodicSubstringInterval=pSubstring.mergedToInterval;
				if (pSubstring.mergedToInterval!=null) {
					pSubstring.mergedToInterval.nImpliedAlignments++;
					pSubstring.mergedToInterval.longestAlignment=Math.max(pSubstring.mergedToInterval.longestAlignment,ReadA.sortedAlignments[j].getALength());
				}
			}
			j++;
		}
	}
	
	
	/**
	 * @return the ID of a peak in $splits$ (i.e. $i$ for the $i$-th block of 3 cells) 
	 * that contains position $pos$, or whose center is at distance at most $threshold$
	 * from $pos$; the procedure returns -1 if no such split can be found.
	 */
	public static final int position2split(int pos, int[] splits, int lastSplit, int threshold) {
		int k, kPrime;
		int distanceLeft, distanceRight;
		if (lastSplit<0) return -1;
		
		k=Arrays.binarySearch(splits,0,lastSplit+1,pos);
		if (k<0) {
			kPrime=-k-1;
			if (kPrime%3==0) {
				distanceRight=(kPrime<=lastSplit&&splits[kPrime+1]-pos<=threshold)?splits[kPrime+1]-pos:-1;
				distanceLeft=(kPrime>0&&pos-splits[kPrime-2]<=threshold)?pos-splits[kPrime-2]:-1;
				if (distanceLeft==-1 && distanceRight==-1) return -1;
				if (distanceRight==-1) k=kPrime-2;
				else if (distanceLeft==-1) k=kPrime+1;
				else if (distanceLeft<=distanceRight) k=kPrime-2;
				else k=kPrime+1;
			}
			else k=kPrime;
		}
		return k/3;
	}
	
	
	/**
	 * @return the position in $intervals$ of the only periodic interval (short or long) 
	 * such that: (1) is contained in the maximal short-period range $[rangeFirst..
	 * rangeLast]$ (which corresponds to the compact range $[fromInterval..toInterval]$ in 
	 * $intervals$); (2) its boundaries fall into the $leftPeak$-th peak in $leftSplits$, 
	 * and into the $rightPeak$-th peak in $rightSplits$; the procedure returns -1 if no 
	 * such interval can be found.
	 */
	private static final int markImpliedAlignments_peaks_getInterval(int fromInterval, int toInterval, int rangeFirst, int rangeLast, int leftPeak, int rightPeak, int threshold) {
		int i;
		int pos, count, interval, leftSplit, rightSplit;
		
		count=0; interval=-1;
		for (i=fromInterval; i<=toInterval; i++) {
			pos=intervals[i].firstPosition;
			leftSplit=position2split(pos,leftSplits,lastLeftSplit,threshold);
			if (leftSplit==-1 || leftSplit!=leftPeak) continue;
			pos=intervals[i].lastPosition;
			if (pos>rangeLast+threshold) continue;
			rightSplit=position2split(pos,rightSplits,lastRightSplit,threshold);
			if (rightSplit==-1 || rightSplit!=rightPeak) continue;
			count++;
			if (count>1) return -1;
			interval=i;
		}
		for (i=toInterval+1; i<=lastInterval; i++) {
			pos=intervals[i].firstPosition;
			if (pos>=rangeLast) break;
			leftSplit=position2split(pos,leftSplits,lastLeftSplit,threshold);
			if (leftSplit==-1 || leftSplit!=leftPeak) continue;
			pos=intervals[i].lastPosition;
			if (pos>rangeLast+threshold) continue;
			rightSplit=position2split(pos,rightSplits,lastRightSplit,threshold);
			if (rightSplit==-1 || rightSplit!=rightPeak) continue;
			count++;
			if (count>1) return -1;
			interval=i;
		}
		for (i=fromInterval-1; i>=0; i--) {
			pos=intervals[i].firstPosition;
			if (pos<rangeFirst-threshold) break;
			leftSplit=position2split(pos,leftSplits,lastLeftSplit,threshold);
			if (leftSplit==-1 || leftSplit!=leftPeak) continue;
			pos=intervals[i].lastPosition;
			if (pos<=rangeFirst || pos>rangeLast+threshold) continue;
			rightSplit=position2split(pos,rightSplits,lastRightSplit,threshold);
			if (rightSplit==-1 || rightSplit!=rightPeak) continue;
			count++;
			if (count>1) return -1;
			interval=i;
		}
		return interval;
	}
	
	
	/**
	 * Remark: the procedure assumes $longPeriodIntervals$ to be sorted by first position.
	 *
	 * @param interval an interval in $intervals$; the procedure converts it into an 
	 * interval in $longPeriodIntervals$ if necessary.
	 */
	private static final PeriodicSubstring markImpliedAlignments_peaks_interval2substring(PeriodicSubstringInterval interval, int leftSplit, int rightSplit) {
		boolean previousSortByID;
		int k;
		PeriodicSubstringInterval pInterval;
		PeriodicSubstring out;
		
		pInterval=interval; out=null;
		if (pInterval.hasLongPeriod) {
			previousSortByID=PeriodicSubstringInterval.sortByID;
			PeriodicSubstringInterval.sortByID=false;
			k=Arrays.binarySearch(longPeriodIntervals,0,lastLongPeriodInterval+1,pInterval);
			PeriodicSubstringInterval.sortByID=previousSortByID;
			if (k<0) {
				System.err.println("markImpliedAlignments_peaks_interval2substring> ERROR: long-period interval not found.");
				System.exit(1);
			}
			pInterval=longPeriodIntervals[k];			
		}
		out=findPeriodicSubstring(leftSplits[3*leftSplit+1],rightSplits[3*rightSplit+1],pInterval);
		if (out==null) {
			out=getArtificialSubstring(pInterval);
			if (out==null) out=addArtificialSubstring(pInterval);
		}
		return out;
	}
	
	
	/**
	 * Adds to $Factors.splits$ the peaks induced by all {left,right}-maximal alignments 
	 * inside maximum ranges of straddling or adjacent long-period intervals that have
	 * been merged to a dense substring of substring type. All alignments inside such 
	 * ranges are used, not just those implied by periodic substrings. This is useful,
	 * since otherwise, once a dense substring of substring type is merged to a long-
	 * period interval, it is not split, so its internal biases would go unnoticed.
	 *
	 * Remark: the procedure uses also peaks spaced by the estimated long period. This is
	 * correct, since alignment intervals created in this way will be filtered out later.
	 *
	 * Remark: the procedure uses arrays $DenseSubstrings.{left,right}Splits$,
	 * $DenseSubstrings.delta{Left,Right}$ and $DenseSubstrings.tmpSplits$.
	 * The procedure stores the last element of $DenseSubstrings.*Splits$ in global 
	 * variable $last*Split_long$.
	 *
	 * Remark: the procedure assumes that at least one long-period interval was merged to
	 * a dense substring of substring type.
	 */
	public static final void splitLongPeriodIntervals() {
		final double CONSTANT_THRESHOLD = 0.15;
		final double MIN_HIGH = 0.3;  // Arbitrary leaf density
		final int MAX_PEAK_RADIUS = Alignments.minAlignmentLength>>2;  // Arbitrary
		final int MAX_WIDTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int MIN_PERIOD_FOR_HISTOGRAM = (MIN_PERIOD_LONG)<<1;  // Arbitrary
		final int MIN_MASS_HIGH = 50*IO.coverage;  // Arbitrary
		int i, j;
		int firstAlignment, startA, endA, from, to, previousLastLeftSplit, previousLastRightSplit;
		int rangeFirst, rangeLast, periodNumerator, periodDenominator;
		int distanceThreshold = IO.quantum<<1;
		PeriodicSubstring tmpSubstring;
		PeriodicSubstringInterval tmpInterval;
		
		// Ensuring the necessary order in intervals and alignments
		if (PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (lastLongPeriodInterval>0) Arrays.sort(longPeriodIntervals,0,lastLongPeriodInterval+1);
		}
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		
		// Iterating over maximal ranges of straddling or adjacent long-period intervals.
		lastLeftSplit_long=-1; lastRightSplit_long=-1;
		firstAlignment=0;
		i=0;
		while (i<=lastLongPeriodInterval && !longPeriodIntervals[i].mergedToDense) i++;
		if (i>lastLongPeriodInterval) return;
		rangeFirst=longPeriodIntervals[i].firstPosition;
		rangeLast=longPeriodIntervals[i].lastPosition;
		periodNumerator=longPeriodIntervals[i].period; periodDenominator=1;
		for (i=i+1; i<=lastLongPeriodInterval; i++) {
			if (!longPeriodIntervals[i].mergedToDense) continue;
			if (longPeriodIntervals[i].firstPosition<=rangeLast+IO.quantum) {
				if (longPeriodIntervals[i].lastPosition>rangeLast) rangeLast=longPeriodIntervals[i].lastPosition;
				if (longPeriodIntervals[i].period>0) {
					periodNumerator+=longPeriodIntervals[i].period;
					periodDenominator++;
				}
				continue;
			}
			startA=rangeFirst; endA=rangeLast;
			// Left-maximal
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> LEFT HISTOGRAM OF ["+startA+".."+endA+"]:  firstAlignment="+firstAlignment);
			ReadA.getDeltaHistogram_periodic(startA,endA,firstAlignment,endA+1,true,DenseSubstrings.deltaLeft);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> ReadA.deltaTmp[0]="+ReadA.deltaTmp[0]);
			if (ReadA.deltaTmp[0]>=(endA-MIN_PERIOD_FOR_HISTOGRAM-startA+1)/MIN_PERIOD_FOR_HISTOGRAM) {
				from=Histograms.firstNonzero(DenseSubstrings.deltaLeft,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> firstNonzero outputs "+from);
				if (from!=-1) {
					to=Histograms.lastNonzero(DenseSubstrings.deltaLeft,0,endA-startA);	
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> lastNonzero outputs "+to);
					previousLastLeftSplit=lastLeftSplit_long;
					lastLeftSplit_long=Histograms.getPeaks(DenseSubstrings.deltaLeft,from,to,DenseSubstrings.leftSplits,lastLeftSplit_long+1,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,DenseSubstrings.tmpSplits,-1/*periodNumerator/periodDenominator*/,MAX_WIDTH,ReadA.deltaPoints,ReadA.lastDeltaPoint,MIN_MASS_HIGH,Alignments.minAlignmentLength);
					for (j=previousLastLeftSplit+1; j<=lastLeftSplit_long; j+=3) {
						Factors.lastSplit++;
						Factors.splits[Factors.lastSplit].clear();
						Factors.splits[Factors.lastSplit].position=DenseSubstrings.leftSplits[j];
						Factors.splits[Factors.lastSplit].nOpen=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitLongPeriodIntervals> added left split "+DenseSubstrings.leftSplits[j]);
					}
				}
			}
			// Right-maximal
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> RIGHT HISTOGRAM OF ["+startA+".."+endA+"]:  firstAlignment="+firstAlignment);
			ReadA.getDeltaHistogram_periodic(startA,endA,firstAlignment,endA+1,false,DenseSubstrings.deltaRight);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> ReadA.deltaTmp[0]="+ReadA.deltaTmp[0]);
			if (ReadA.deltaTmp[0]>=(endA-MIN_PERIOD_FOR_HISTOGRAM-startA+1)/MIN_PERIOD_FOR_HISTOGRAM) {
				from=Histograms.firstNonzero(DenseSubstrings.deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> firstNonzero outputs "+from);
				if (from!=-1) {
					to=Histograms.lastNonzero(DenseSubstrings.deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> lastNonzero outputs "+to);
					previousLastRightSplit=lastRightSplit_long;
					lastRightSplit_long=Histograms.getPeaks(DenseSubstrings.deltaRight,from,to,DenseSubstrings.rightSplits,lastRightSplit_long+1,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,DenseSubstrings.tmpSplits,-1/*periodNumerator/periodDenominator*/,MAX_WIDTH,ReadA.deltaPoints,ReadA.lastDeltaPoint,MIN_MASS_HIGH,Alignments.minAlignmentLength);
					for (j=previousLastRightSplit+1; j<=lastRightSplit_long; j+=3) {
						Factors.lastSplit++;
						Factors.splits[Factors.lastSplit].clear();
						Factors.splits[Factors.lastSplit].position=DenseSubstrings.rightSplits[j];
						Factors.splits[Factors.lastSplit].nClosed=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitLongPeriodIntervals> added right split "+DenseSubstrings.rightSplits[j]);
					}
				}
			}
			// Next range
			rangeFirst=longPeriodIntervals[i].firstPosition;
			rangeLast=longPeriodIntervals[i].lastPosition;
			periodNumerator=intervals[i].period; periodDenominator=1;
			firstAlignment=ReadA.deltaTmp[1];
		}
		startA=rangeFirst; endA=rangeLast;
		// Last left-maximal
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> LEFT HISTOGRAM OF ["+startA+".."+endA+"]:  firstAlignment="+firstAlignment);
		ReadA.getDeltaHistogram_periodic(startA,endA,firstAlignment,endA+1,true,DenseSubstrings.deltaLeft);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> ReadA.deltaTmp[0]="+ReadA.deltaTmp[0]);
		if (ReadA.deltaTmp[0]>=(endA-MIN_PERIOD_FOR_HISTOGRAM-startA+1)/MIN_PERIOD_FOR_HISTOGRAM) {
			from=Histograms.firstNonzero(DenseSubstrings.deltaLeft,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> firstNonzero outputs "+from);
			if (from!=-1) {
				to=Histograms.lastNonzero(DenseSubstrings.deltaLeft,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> lastNonzero outputs "+to);
				previousLastLeftSplit=lastLeftSplit_long;
				lastLeftSplit_long=Histograms.getPeaks(DenseSubstrings.deltaLeft,from,to,DenseSubstrings.leftSplits,lastLeftSplit_long+1,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,DenseSubstrings.tmpSplits,-1/*periodNumerator/periodDenominator*/,MAX_WIDTH,ReadA.deltaPoints,ReadA.lastDeltaPoint,MIN_MASS_HIGH,Alignments.minAlignmentLength);
				for (j=previousLastLeftSplit+1; j<=lastLeftSplit_long; j+=3) {
					Factors.lastSplit++;
					Factors.splits[Factors.lastSplit].clear();
					Factors.splits[Factors.lastSplit].position=DenseSubstrings.leftSplits[j];
					Factors.splits[Factors.lastSplit].nOpen=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitLongPeriodIntervals> added left split "+DenseSubstrings.leftSplits[j]);
				}
			}
		}
		// Last right-maximal
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> RIGHT HISTOGRAM OF ["+startA+".."+endA+"]:  firstAlignment="+firstAlignment);
		ReadA.getDeltaHistogram_periodic(startA,endA,firstAlignment,endA+1,false,DenseSubstrings.deltaRight);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> ReadA.deltaTmp[0]="+ReadA.deltaTmp[0]);
		if (ReadA.deltaTmp[0]>=(endA-MIN_PERIOD_FOR_HISTOGRAM-startA+1)/MIN_PERIOD_FOR_HISTOGRAM) {
			from=Histograms.firstNonzero(DenseSubstrings.deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> firstNonzero outputs "+from);
			if (from!=-1) {
				to=Histograms.lastNonzero(DenseSubstrings.deltaRight,0,endA-startA);
if (IO.SHOW_STD_ERR_PRIME) System.err.println("splitLongPeriodIntervals> lastNonzero outputs "+to);
				previousLastRightSplit=lastRightSplit_long;
				lastRightSplit_long=Histograms.getPeaks(DenseSubstrings.deltaRight,from,to,DenseSubstrings.rightSplits,lastRightSplit_long+1,startA,IO.maxDeltaStd,distanceThreshold,distanceThreshold,true,CONSTANT_THRESHOLD,MIN_HIGH,MAX_PEAK_RADIUS,DenseSubstrings.tmpSplits,-1/*periodNumerator/periodDenominator*/,MAX_WIDTH,ReadA.deltaPoints,ReadA.lastDeltaPoint,MIN_MASS_HIGH,Alignments.minAlignmentLength);
				for (j=previousLastRightSplit+1; j<=lastRightSplit_long; j+=3) {
					Factors.lastSplit++;
					Factors.splits[Factors.lastSplit].clear();
					Factors.splits[Factors.lastSplit].position=DenseSubstrings.rightSplits[j];
					Factors.splits[Factors.lastSplit].nClosed=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitLongPeriodIntervals> added right split "+DenseSubstrings.rightSplits[j]);
				}
			}
		}
		
		// Unmarking maximal alignments between peaks
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if (ReadA.sortedAlignments[i].isLeftMaximal!=1 || ReadA.sortedAlignments[i].isRightMaximal!=1) continue;
			tmpInterval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if (tmpInterval==null || !tmpInterval.hasLongPeriod || !tmpInterval.mergedToDense) continue;
			startA=ReadA.sortedAlignments[i].startA();
			j=Arrays.binarySearch(DenseSubstrings.leftSplits,0,lastLeftSplit_long+1,startA);
			if ( !( j>=0 || 
				    ( -j-1>0 && -j-1<=lastLeftSplit_long &&
				      ( (-j-1)%3!=0 ||
				        Math.abs(DenseSubstrings.leftSplits[-j],startA)<=DISTANCE_THRESHOLD ||
					    Math.abs(DenseSubstrings.leftSplits[-j-2],startA)<=DISTANCE_THRESHOLD
				      )
				    )
			      )
			   ) continue;
			endA=ReadA.sortedAlignments[i].endA();
			j=Arrays.binarySearch(DenseSubstrings.rightSplits,0,lastRightSplit_long+1,endA);
			if ( !( j>=0 || 
				    ( -j-1>0 && -j-1<=lastRightSplit_long &&
				      ( (-j-1)%3!=0 ||
				        Math.abs(DenseSubstrings.rightSplits[-j],endA)<=DISTANCE_THRESHOLD ||
					    Math.abs(DenseSubstrings.rightSplits[-j-2],endA)<=DISTANCE_THRESHOLD
				      )
				    )
			      )
			   ) continue;
			ReadA.sortedAlignments[i].periodicSubstringInterval=null;
			ReadA.sortedAlignments[i].impliedByPeriodicSubstring=null;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("splitLongPeriodIntervals> UNMARKING alignment "+ReadA.sortedAlignments[i].toStringBoundaries());
		}
	}
	
	
	/**
	 * Sets $hasLongPeriod=false$ for all long-period substrings whose readB interval
	 * is approximately identical to the readB interval of a short-period substring.
	 * Such substrings are not moved from $longPeriodSubstrings$ to $substrings$, so they
	 * will create alignment intervals in the following steps of the pipeline.
	 *
	 * @return the number of long-period substring that have been marked as short-period 
	 * by the procedure.
	 */
	private static final int longPeriod2shortPeriod_substrings() {
		int i, j;
		int firstJForNextI, previousOrder, previousOrderLong, out;
		
		// Ensuring the necessary order in $substrings$ and $longPeriodSubstrings$.
		previousOrder=PeriodicSubstring.order;
		if (PeriodicSubstring.order!=PeriodicSubstring.READB_ORIENTATION_MINSTARTB) {
			PeriodicSubstring.order=PeriodicSubstring.READB_ORIENTATION_MINSTARTB;
			if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		}
		previousOrderLong=PeriodicSubstring.order_longPeriod;
		if (PeriodicSubstring.order_longPeriod!=PeriodicSubstring.READB_ORIENTATION_MINSTARTB) {
			PeriodicSubstring.order_longPeriod=PeriodicSubstring.READB_ORIENTATION_MINSTARTB;
			if (lastLongPeriodSubstring>0) Arrays.sort(longPeriodSubstrings,0,lastLongPeriodSubstring+1);
		}
		
		// Marking
		i=0; j=0; firstJForNextI=-1; out=0;
		while (i<=lastLongPeriodSubstring) {
			if (j>lastSubstring || substrings[j].readB!=longPeriodSubstrings[i].readB || substrings[j].orientation!=longPeriodSubstrings[i].orientation) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (substrings[j].maxEndB<=longPeriodSubstrings[i].minStartB) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastLongPeriodSubstring && longPeriodSubstrings[i+1].readB==longPeriodSubstrings[i].readB && longPeriodSubstrings[i+1].orientation==longPeriodSubstrings[i].orientation && substrings[j].maxEndB>=longPeriodSubstrings[i+1].minStartB) firstJForNextI=j;
			if (Intervals.areApproximatelyIdentical_lowQuality(substrings[j].minStartB,substrings[j].maxEndB,longPeriodSubstrings[i].minStartB,longPeriodSubstrings[i].maxEndB,longPeriodSubstrings[i].readB)) {
				longPeriodSubstrings[i].hasLongPeriod=false;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				out++;
				continue;
			}
			j++;
		}
		
		// Restoring the original order
		if (PeriodicSubstring.order!=previousOrder) {
			PeriodicSubstring.order=previousOrder;
			if (lastSubstring>0) Arrays.sort(substrings,0,lastSubstring+1);
		}
		if (PeriodicSubstring.order_longPeriod!=previousOrderLong) {
			PeriodicSubstring.order_longPeriod=previousOrderLong;
			if (lastLongPeriodSubstring>0) Arrays.sort(longPeriodSubstrings,0,lastLongPeriodSubstring+1);
		}
		
		return out;
	}
	
	
	/**
	 * Moves from $longPeriodIntervals$ to $intervals$ all intervals with 
	 * $hasLongPeriod=false$, and sets their $representative$ field to NULL.
	 *
	 * @return true iff at least one interval was moved.
	 */
	private static final boolean longPeriod2shortPeriod_moveIntervals() {
		final int GROWTH_RATE = 10;  // Arbitrary
		int i, j, k;
		int previousLastInterval;
		PeriodicSubstringInterval tmpInterval;
		PeriodicSubstringInterval[] newIntervals;
		
		j=-1; previousLastInterval=lastInterval;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (!longPeriodIntervals[i].hasLongPeriod) {
				longPeriodIntervals[i].representative=null;
				lastInterval++;
				if (lastInterval==intervals.length) {
					newIntervals = new PeriodicSubstringInterval[intervals.length+GROWTH_RATE];
					System.arraycopy(intervals,0,newIntervals,0,intervals.length);
					for (k=intervals.length; k<newIntervals.length; k++) {
						newIntervals[k] = new PeriodicSubstringInterval();
						newIntervals[k].hasLongPeriod=false;
					}
					intervals=newIntervals;
				}
				tmpInterval=intervals[lastInterval];
				intervals[lastInterval]=longPeriodIntervals[i];
				longPeriodIntervals[i]=tmpInterval;
				longPeriodIntervals[i].hasLongPeriod=true;
			}
			else {
				j++;
				if (j==i) continue;
				tmpInterval=longPeriodIntervals[j];
				longPeriodIntervals[j]=longPeriodIntervals[i];
				longPeriodIntervals[i]=tmpInterval;
			}
		}
		lastLongPeriodInterval=j;
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
		if (lastInterval>previousLastInterval) {
			Arrays.sort(intervals,0,lastInterval+1);
			estimatePeriodOfShortPeriodIntervals();
			return true;
		}
		return false;
	}
	
	
	/**
	 * Remark: the procedure assumes $intervals$ to be sorted by $firstPosition$.
	 *
	 * @return TRUE iff $pos$ is strictly inside a short-period interval.
	 */
	public static final boolean inShortPeriodInterval(int pos) {
		final int THRESHOLD = IO.quantum;
		int i, j;
		
		tmpPeriodicInterval.firstPosition=pos;
		j=Arrays.binarySearch(intervals,0,lastInterval+1,tmpPeriodicInterval);
		if (j<0) j=-j-1;
		for (i=j-1; i>=0; i--) {
			if (!intervals[i].hasLongPeriod && intervals[i].firstPosition<pos-THRESHOLD && intervals[i].lastPosition>pos+THRESHOLD) return true;
		}
		return false;
	}
	
	
	private static final void ensureSpace_shifts(int nShifts) {
		final double GROWTH_RATE = 1.5;  // Arbitrary
		if (shifts.length>=nShifts) return;
		
		Point[] newShifts = new Point[(int)(nShifts*GROWTH_RATE)];
		System.arraycopy(shifts,0,newShifts,0,shifts.length);
		for (int i=shifts.length; i<newShifts.length; i++) newShifts[i] = new Point();
		shifts=newShifts;
	}
	
	
	/**
	 * Simple variant of $discardLongPeriodIntervals()$ that discards just long-period
	 * intervals that are approx. identical to short-period intervals, and that updates
	 * pointers in periodic substrings and alignments.
	 */
	private static final void discardLongPeriodIntervals_simple() {
		int i, j;
		int firstPosition, lastPosition, firstJForNextI;
		PeriodicSubstringInterval tmpInterval;
		
		// Discarding intervals
		for (i=0; i<=lastInterval; i++) {
			intervals[i].tmpInt1=intervals[i].firstPosition;  // Used as temporary space
			intervals[i].tmpInt2=intervals[i].lastPosition;  // Used as temporary space
		}		
		for (i=0; i<=lastLongPeriodInterval; i++) longPeriodIntervals[i].representative=null;
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastLongPeriodInterval) {
			firstPosition=longPeriodIntervals[i].firstPosition;
			lastPosition=longPeriodIntervals[i].lastPosition;
			if (j>lastInterval || intervals[j].firstPosition>=lastPosition) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (intervals[j].lastPosition<=firstPosition) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastLongPeriodInterval && intervals[j].lastPosition>=longPeriodIntervals[i+1].firstPosition) firstJForNextI=j;
			if (Intervals.areApproximatelyIdentical(firstPosition,lastPosition,intervals[j].firstPosition,intervals[j].lastPosition)) {
				longPeriodIntervals[i].representative=intervals[j];
				intervals[j].simpleMerge(longPeriodIntervals[i]);
				intervals[j].tmpInt1=Math.min(intervals[j].tmpInt1,longPeriodIntervals[i].firstPosition);
				intervals[j].tmpInt2=Math.max(intervals[j].tmpInt2,longPeriodIntervals[i].lastPosition);
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			j++;
		}
		j=-1;
		for (i=0; i<=lastLongPeriodInterval; i++) {
			if (longPeriodIntervals[i].representative!=null) continue;
			j++;
			if (j!=i) {
				tmpInterval=longPeriodIntervals[j];
				longPeriodIntervals[j]=longPeriodIntervals[i];
				longPeriodIntervals[i]=tmpInterval;
			}
		}
		lastLongPeriodInterval=j;
		for (i=0; i<=lastInterval; i++) {
			intervals[i].firstPosition=intervals[i].tmpInt1;
			intervals[i].lastPosition=intervals[i].tmpInt2;
		}

		// Updating pointers from periodic substrings and from alignments
        for (i=0; i<=lastLongPeriodSubstring; i++) {
			if (!longPeriodSubstrings[i].hasLongPeriod) continue;
            tmpInterval=longPeriodSubstrings[i].mergedToInterval;
            if (tmpInterval==null) continue;
			if (tmpInterval.representative!=null) longPeriodSubstrings[i].mergedToInterval=tmpInterval.representative;
        }
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if ( ReadA.sortedAlignments[i].periodicSubstringInterval!=null && 
				 ReadA.sortedAlignments[i].periodicSubstringInterval.hasLongPeriod &&
				 ReadA.sortedAlignments[i].periodicSubstringInterval.representative!=null
			   ) ReadA.sortedAlignments[i].periodicSubstringInterval=ReadA.sortedAlignments[i].periodicSubstringInterval.representative;
		}
	}
	

}