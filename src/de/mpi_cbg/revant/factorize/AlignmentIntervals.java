package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.DAG;
import de.mpi_cbg.revant.util.Histograms;
import de.mpi_cbg.revant.util.DensityEstimationTree;


public class AlignmentIntervals {
	/**
	 * Parameters of the pipeline
	 */
	private static int MAX_NEIGHBORS;
	public static final int MAX_DISTINCT_SHIFTS = 100;  // Maximum number of distinct observed shifts
	private static final int MIN_MERGED_THRESHOLD = (IO.coverage<<1)+1;  // Arbitrary
	private static final int MIN_INTERVALS_IN_NONMAXIMAL = (IO.minOccurrencesInGenome*IO.coverage-1)<<1;  // Arbitrary
	private static final int MIN_INTERVALS_IN_MAXIMAL = 
		(MIN_INTERVALS_IN_NONMAXIMAL)<<1;  // Arbitrary
	private static final int MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX = (MIN_INTERVALS_IN_MAXIMAL)<<1;  // Arbitrary
	private static final int MIN_ALIGNMENTS_LARGE = IO.coverage*50;  // Arbitrary
	private static final int MIN_PATH_LENGTH = Math.max( 
		4,  // Arbitrary
		(IO.minOccurrencesInGenome-1)*IO.coverage  // To comply with $Factorize.discardRareIntervals()$
	);
	
	/**
	 * Output data structures
	 */
	public static AlignmentInterval[] intervals;
	public static int lastInterval;
	
	/**
	 * DAG data structures
	 */
	private static int[][] inNeighbors, outNeighbors;  // DAG of alignments
	private static int[] nInNeighbors, nInNeighborsPrime, nOutNeighbors;
	private static int[] minimalVertices;  // Minimal vertices of the DAG
	private static int[] components;  // Connected component of each vertex of the DAG
	private static int[] nVerticesPerComponent;  // Number of minimal vertices per connected component of the DAG
	private static int[] sorted2original;  // Maps a rank in the topological order of the DAG to the ID of the corresponding vertex
	private static int[] original2sorted;  // Reverses $sorted2original$
	private static int[] distances;  // Length of a best path from a source alignment to all other alignments
	private static int[] predecessors;  // Predecessor of each vertex in the best path from a source alignment
	
	/**
	 * Temporary space
	 */
	public static int[] stack;
	public static Point[] nMergedIntervalsPoints;
	private static Alignment tmpAlignment;
	private static DenseSubstring tmpSubstring;
	private static Point[] left, right;
	private static int lastLeft, lastRight;
	private static final int[] tmp = new int[2];
	private static int[][] positions;
	
	
	public static final void allocateMemory(int maxAlignments, int maxNeighbors) {
		int i;
		
		MAX_NEIGHBORS=maxNeighbors;
		
		// Output data structures
		intervals = new AlignmentInterval[maxAlignments];
		for (i=0; i<maxAlignments; i++) intervals[i] = new AlignmentInterval();
		
		// DAG data structures
		inNeighbors = new int[maxAlignments][MAX_NEIGHBORS];
		nInNeighbors = new int[maxAlignments];
		nInNeighborsPrime = new int[maxAlignments];
		outNeighbors = new int[maxAlignments][MAX_NEIGHBORS];
		nOutNeighbors = new int[maxAlignments];
		minimalVertices = new int[maxAlignments];
		components = new int[maxAlignments];
		nVerticesPerComponent = new int[maxAlignments];
		sorted2original = new int[maxAlignments];
		original2sorted = new int[maxAlignments];
		distances = new int[maxAlignments];
		predecessors = new int[maxAlignments];
		
		// Temporary space
		stack = new int[maxAlignments];
		nMergedIntervalsPoints = new Point[maxAlignments];
		for (i=0; i<nMergedIntervalsPoints.length; i++) nMergedIntervalsPoints[i] = new Point();
		tmpAlignment = new Alignment();
		tmpSubstring = new DenseSubstring();
		left = new Point[maxAlignments];
		for (i=0; i<left.length; i++) left[i] = new Point();
		right = new Point[maxAlignments];
		for (i=0; i<right.length; i++) right[i] = new Point();
		positions = new int[maxAlignments][4];
	}
	
	
	/**
	 * The procedure assumes that periodic substrings and dense substrings have already
	 * been computed.
	 *
	 * Remark: alignment intervals could detect simple repeats, or even modules inside 
	 * TEs. See e.g. \cite{tempel2010moduleorganizer}.
	 */
	public static final void detect() {
		boolean isSorted;
		int i, firstImplied;
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("AlignmentIntervals.detect> alignments before getAlignmentIntervals():");
Alignment.order=Alignment.IMPL_STARTA;
Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
	for (i=0; i<=ReadA.lastSortedAlignment; i++) System.err.println(ReadA.sortedAlignments[i].toStringBoundaries()+" :: "+ReadA.sortedAlignments[i].toStringPointers());
	System.err.println();
	IO.printErr("AlignmentIntervals.detect> periodic substring intervals at the very beginning:");
	for (i=0; i<=PeriodicSubstrings.lastInterval; i++) IO.printErr(PeriodicSubstrings.intervals[i]);
}		
		
		
		
		
		firstImplied=getAlignmentIntervals();
		if (firstImplied==0) return;
		

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignment intervals after getAlignmentIntervals():");
	for (i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments after getAlignmentIntervals():");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],null,"+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		else IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode()+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		IO.printErr("   isLeftMaximal="+ReadA.sortedAlignments[i].isLeftMaximal+" isRightMaximal="+ReadA.sortedAlignments[i].isRightMaximal);
		IO.printErr("   isLeftMaximalB="+ReadA.sortedAlignments[i].isLeftMaximalB+" isRightMaximalB="+ReadA.sortedAlignments[i].isRightMaximalB);
		IO.printErr("   mergedToInterval="+ReadA.sortedAlignments[i].mergedToInterval);
		IO.printErr(ReadA.sortedAlignments[i].toStringPointers());
	}
	System.err.println("Consistency checks 1");
	checkConsistency();
}
		
		
		mergeSimilarIntervals(firstImplied);
				
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignment intervals after mergeSimilarIntervals():");
	for (i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments after mergeSimilarIntervals():");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],null,"+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		else IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode()+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		IO.printErr("   isLeftMaximal="+ReadA.sortedAlignments[i].isLeftMaximal+" isRightMaximal="+ReadA.sortedAlignments[i].isRightMaximal);
		IO.printErr("   mergedToInterval="+ReadA.sortedAlignments[i].mergedToInterval);
		IO.printErr(ReadA.sortedAlignments[i].toStringPointers());
	}
	System.err.println("Consistency checks 2");
	checkConsistency();
}			
		
		
/*		isSorted=filterIntervalsInPeriodicSubstrings(firstImplied);
		if (!isSorted) {
			Alignment.order=Alignment.IMPL_STARTA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignment intervals after filterIntervalsInPeriodicSubstrings():");
	for (i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments after filterIntervalsInPeriodicSubstrings():");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],null,"+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		else IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode()+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		IO.printErr("   isLeftMaximal="+ReadA.sortedAlignments[i].isLeftMaximal+" isRightMaximal="+ReadA.sortedAlignments[i].isRightMaximal);
		IO.printErr("   mergedToInterval="+ReadA.sortedAlignments[i].mergedToInterval);
		IO.printErr(ReadA.sortedAlignments[i].toStringPointers());
	}
	System.err.println("Consistency checks 3");
	checkConsistency();
}
*/

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignment intervals immediately before cleanAlignmentIntervals():");
	for (i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments immediately before cleanAlignmentIntervals():");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],null,"+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		else IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode()+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		IO.printErr("   isLeftMaximal="+ReadA.sortedAlignments[i].isLeftMaximal+" isRightMaximal="+ReadA.sortedAlignments[i].isRightMaximal);
		IO.printErr("   mergedToInterval="+ReadA.sortedAlignments[i].mergedToInterval);
		IO.printErr(ReadA.sortedAlignments[i].toStringPointers());
	}
}	


		
		
		cleanAlignmentIntervals(firstImplied);
		// $ReadA.sortedAlignments$ might not be sorted at this point.
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("periodic substring intervals immediately after cleanAlignmentIntervals():");
	for (i=0; i<=PeriodicSubstrings.lastInterval; i++) IO.printErr(PeriodicSubstrings.intervals[i]);
	
	IO.printErr("alignment intervals immediately after cleanAlignmentIntervals():");
	for (i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments after cleanAlignmentIntervals():");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],null,"+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		else IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode()+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		IO.printErr("   isLeftMaximal="+ReadA.sortedAlignments[i].isLeftMaximal+" isRightMaximal="+ReadA.sortedAlignments[i].isRightMaximal);
		IO.printErr("   mergedToInterval="+ReadA.sortedAlignments[i].mergedToInterval);
		IO.printErr(ReadA.sortedAlignments[i].toStringPointers());
	}
	System.err.println("Consistency checks 4");
	checkConsistency();
}		
		
		
		markIntervalsInDenseSubstrings();
		// Now $intervals$ is already sorted by $firstPosition$
		

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("detect> ALIGNMENT INTERVALS:");
	for (i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("");
	IO.printErr("alignments after setMaximality():");
	for (i=0; i<=ReadA.lastSortedAlignment; i++) {
		if (ReadA.sortedAlignments[i].impliedByDenseSubstring==null) IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],null,"+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		else IO.printErr(i+",["+Alignments.alignments[ReadA.sortedAlignments[i].id][3]+".."+Alignments.alignments[ReadA.sortedAlignments[i].id][4]+"],"+ReadA.sortedAlignments[i].impliedByDenseSubstring.hashCode()+","+(ReadA.sortedAlignments[i].inDenseSubstring==null?"null":ReadA.sortedAlignments[i].inDenseSubstring.hashCode()));
		IO.printErr("   isLeftMaximal="+ReadA.sortedAlignments[i].isLeftMaximal+" isRightMaximal="+ReadA.sortedAlignments[i].isRightMaximal);
		IO.printErr("   isLeftMaximalB="+ReadA.sortedAlignments[i].isLeftMaximalB+" isRightMaximalB="+ReadA.sortedAlignments[i].isRightMaximalB);
		IO.printErr("   mergedToInterval="+ReadA.sortedAlignments[i].mergedToInterval);
		IO.printErr(ReadA.sortedAlignments[i].toStringPointers());
	}
}
		
		// At this point, an alignment can have both $mergedToInterval!=null$ and
	 	// $inDenseSubstring!=null$.
		if (IO.CONSISTENCY_CHECKS) checkConsistency();
	}
	

	/**
	 * Merges elements of $ReadA.sortedAlignments$ with compatible readA intervals. See
	 * also $getAlignmentIntervals_weak()$ for details. The procedure considers only 
	 * alignments that are neither implied by a periodic substring nor by a dense
	 * substring.
	 *
	 * Remark: if two repeats AB and BC overlap at B, but there is no left- and right-
	 * maximal alignment that spans just B, the procedure does not report B as a separate 
	 * interval, even though there is an alignment that includes B, starts at the first 
	 * position of B, and is left-maximal, and an alignment that includes B, ends at the 
	 * last position of B, and is right-maximal. This is because B is always either 
	 * preceded by A or followed by C, i.e. every occurrence of B in the genome can be 
	 * inferred from the occurrences of AB and BC.
	 *
	 * Remark: assume that readA is contained in a repeat that is longer than the longest 
	 * read (for example, the mitochondrial genome). At the end of this procedure, such 
	 * read contains either a single alignment interval (which is neither left- nor
	 * right-maximal), or multiple alignment intervals (that are neither left- nor right-
	 * maximal), contained in one another, generated by all alignments of other reads with
	 * prefixes and suffixes of readA.
	 *
	 * Remark: alignments that stop at a short, low-quality region, are not merged to 
	 * the (usually fewer) alignments that cross the region. If they were assigned to the 
	 * same interval, such interval might be discarded in later stages of the pipeline, 
	 * because the signal at its endpoints might be too weak, and all its alignments might
	 * get reassigned to a container interval. Not merging is not likely to be a problem 
	 * in later stages of the pipeline, since it could be resolved e.g. by removing
	 * non-supermaximal repeats.
	 *
	 * Remark: the procedure sets the $nMaximalAlignments*$ fields of every interval.
	 *
	 * Remark: the procedure changes the order of $ReadA.sortedAlignments$.
	 * At the end of the procedure, $intervals$ is not sorted.
	 *
	 * Remark: the procedure uses the $is{Left,Right}Maximal$ fields of 
	 * $AlignmentInterval$ as temporary space. The correct value of such fields has to be
	 * set by $setMaximality()$.
	 *
	 * @return the procedure puts all alignments in $ReadA.sortedAlignments$ that are 
	 * implied by a dense or periodic substring, in a suffix block $[X..]$, where $X$ is 
	 * returned in output; alignments in each of the two blocks are sorted by startA.
	 */
	private static final int getAlignmentIntervals() {
		final int THRESHOLD = (3*IO.quantum)>>1;
		final int MAX_DISTANCE_THRESHOLD = DenseSubstrings.maxDistanceEndWeak;
		final int MIN_INTERVALS_FOR_LONG_REPEAT = 2;
		final int MIN_TRIMMED_LENGTH = Alignments.minAlignmentLength;  // Arbitrary
		final int SPLITS_THRESHOLD = IO.quantum;
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int N_ALIGNMENTS_THRESHOLD = 50*IO.coverage;  // Arbitrary
		boolean isMaximal, trimLeft, trimRight;
		int i, j, k;
		int firstImplied, min, max, minLeft, maxLeft, minJ, maxJ, start, end, nAlignments;
		int idI, idJ, sumFirst, sumLast, startA, endA, nMergedIntervals, nMaximalAlignmentsLeft, nMaximalAlignmentsRight;
		int isLeftMaximal, isRightMaximal, discarded;
		AlignmentInterval tmpInterval;
		PeriodicSubstring tmpSubstring;
		lastInterval=-1;
		if (ReadA.lastSortedAlignment==-1) return 0;
		
		// Disconnecting alignments that are implied by a periodic substring, but that 
		// form nonetheless a short alignment interval strictly inside the periodic
		// substring. This happens rarely.
		if (PeriodicSubstrings.lastInterval>=0) {
			if (Alignment.order!=Alignment.IMPL_PERIODICSUBSTRING_STARTA_ENDA) {
				Alignment.order=Alignment.IMPL_PERIODICSUBSTRING_STARTA_ENDA;
				if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
			}
			i=0;
			while (i<=ReadA.lastSortedAlignment && ReadA.sortedAlignments[i].impliedByPeriodicSubstring==null) i++;
			if (i<=ReadA.lastSortedAlignment) {
				for (j=i; j<=ReadA.lastSortedAlignment; j++) {
					ReadA.sortedAlignments[j].markedPrefix=false;  // Temporary flag
				}
				while (i<=ReadA.lastSortedAlignment) {
					startA=ReadA.sortedAlignments[i].startA();
					endA=ReadA.sortedAlignments[i].endA();
					tmpSubstring=ReadA.sortedAlignments[i].impliedByPeriodicSubstring;
					if ( Math.abs(startA,tmpSubstring.minStartA)<=IDENTITY_THRESHOLD || Math.abs(endA,tmpSubstring.maxEndA)<=IDENTITY_THRESHOLD ||
						 (tmpSubstring.nStartA>0 && Math.abs(startA,tmpSubstring.sumStartA/tmpSubstring.nStartA)<=IDENTITY_THRESHOLD) ||
						 (tmpSubstring.nEndA>0 && Math.abs(endA,tmpSubstring.sumEndA/tmpSubstring.nEndA)<=IDENTITY_THRESHOLD)
					   ) {
						i++;
						continue;
					}
					nAlignments=1; j=i+1;
					while (j<=ReadA.lastSortedAlignment) {
						if (ReadA.sortedAlignments[j].startA()>startA+IDENTITY_THRESHOLD) break;
						tmpSubstring=ReadA.sortedAlignments[j].impliedByPeriodicSubstring;
						if (tmpSubstring!=ReadA.sortedAlignments[i].impliedByPeriodicSubstring) break;
						if ( Math.abs(ReadA.sortedAlignments[j].startA(),tmpSubstring.minStartA)<=IDENTITY_THRESHOLD || Math.abs(ReadA.sortedAlignments[j].endA(),tmpSubstring.maxEndA)<=IDENTITY_THRESHOLD ||
							 (tmpSubstring.nStartA>0 && Math.abs(ReadA.sortedAlignments[j].startA(),tmpSubstring.sumStartA/tmpSubstring.nStartA)<=IDENTITY_THRESHOLD) ||
							 (tmpSubstring.nEndA>0 && Math.abs(ReadA.sortedAlignments[j].endA(),tmpSubstring.sumEndA/tmpSubstring.nEndA)<=IDENTITY_THRESHOLD)
						   ) {
							j++;
							continue;
						}
						if (Intervals.areApproximatelyIdentical(ReadA.sortedAlignments[j].startA(),ReadA.sortedAlignments[j].endA(),startA,endA)) {
							nAlignments++;
							ReadA.sortedAlignments[j].markedPrefix=true;
						}
						j++;
					}
					if (nAlignments>=N_ALIGNMENTS_THRESHOLD) {	
						k=j-1;
						while (k>i) {
							if (ReadA.sortedAlignments[k].markedPrefix) {
								ReadA.sortedAlignments[k].impliedByPeriodicSubstring=null;
								ReadA.sortedAlignments[k].periodicSubstringInterval=null;
							}
							k--;
						}
						ReadA.sortedAlignments[i].impliedByPeriodicSubstring=null;
						ReadA.sortedAlignments[i].periodicSubstringInterval=null;
					}
					i=j;
				}
			}
		}
		
		// Ensuring the necessary order in $sortedAlignments$
		if (Alignment.order!=Alignment.IMPL_STARTA) {
			Alignment.order=Alignment.IMPL_STARTA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		firstImplied=ReadA.lastSortedAlignment+1;
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if (ReadA.sortedAlignments[i].impliedByDenseSubstring!=null || ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null || ReadA.sortedAlignments[i].equalToSubstringType()) {
				firstImplied=i;
				break;
			}
		}
		if (firstImplied==0) return 0;

		// Detecting intervals from alignments that are both left- and right-maximal.
		//
		// Remark: here we do a greedy left-to-right merging in which, for every interval
		// that is not already merged, we merge all following intervals in startA order
		// that are approximately identical to it. Everywhere else in the code we compute
		// instead the connected component of the $areApproximatelyIdentical$ relation.
		// Building connected components here is very likely to result in inaccurate
		// boundaries due to drifting (e.g. when alignments are inside a dense substring).
		for (i=0; i<firstImplied; i++) ReadA.sortedAlignments[i].mergedToInterval=null;
		for (i=0; i<firstImplied; i++) {
			if (ReadA.sortedAlignments[i].mergedToInterval!=null || ReadA.sortedAlignments[i].isLeftMaximal!=1 || ReadA.sortedAlignments[i].isRightMaximal!=1) continue;
			idI=ReadA.sortedAlignments[i].id;
			sumFirst=Alignments.alignments[idI][3]; 
			sumLast=Alignments.alignments[idI][4];
			min=Alignments.alignments[idI][3];
			max=Alignments.alignments[idI][4];
			lastInterval++;
			ReadA.sortedAlignments[i].mergedToInterval=intervals[lastInterval];
			nMergedIntervals=1; 
			nMaximalAlignmentsLeft=ReadA.sortedAlignments[i].isLeftMaximal==1?1:0;
			nMaximalAlignmentsRight=ReadA.sortedAlignments[i].isRightMaximal==1?1:0;
			for (j=i+1; j<firstImplied; j++) {
				idJ=ReadA.sortedAlignments[j].id;
				if (Alignments.alignments[idJ][3]>Alignments.alignments[idI][4]) break;
				if (ReadA.sortedAlignments[j].mergedToInterval!=null || ReadA.sortedAlignments[j].isLeftMaximal!=1 || ReadA.sortedAlignments[j].isRightMaximal!=1) continue;
				if (Intervals.areApproximatelyIdentical(Alignments.alignments[idJ][3],Alignments.alignments[idJ][4],Alignments.alignments[idI][3],Alignments.alignments[idI][4])) {
					sumFirst+=Alignments.alignments[idJ][3];
					sumLast+=Alignments.alignments[idJ][4];
					min=Math.min(min,Alignments.alignments[idJ][3]);
					max=Math.max(max,Alignments.alignments[idJ][4]);
					nMergedIntervals++; 
					nMaximalAlignmentsLeft+=ReadA.sortedAlignments[j].isLeftMaximal==1?1:0;
					nMaximalAlignmentsRight+=ReadA.sortedAlignments[j].isRightMaximal==1?1:0;
					ReadA.sortedAlignments[j].mergedToInterval=intervals[lastInterval];	
				}
			}
			j=sumFirst/nMergedIntervals;
			intervals[lastInterval].firstPosition=min<j-THRESHOLD?min:j;
			j=sumLast/nMergedIntervals;
			intervals[lastInterval].lastPosition=max>j+THRESHOLD?max:j;
			intervals[lastInterval].nMergedIntervals=nMergedIntervals;
			intervals[lastInterval].nMaximalAlignmentsLeft=nMaximalAlignmentsLeft;
			intervals[lastInterval].nMaximalAlignmentsRight=nMaximalAlignmentsRight;
			intervals[lastInterval].isLongRepeat=false;
			intervals[lastInterval].isLeftMaximal=true;
			intervals[lastInterval].isRightMaximal=true;
			intervals[lastInterval].id=0;
		}
		
		// Marking all remaining maximal alignments that could be explained by the
		// detected intervals
		for (i=0; i<=lastInterval; i++) {
			intervals[i].id=i;
			positions[i][0]=intervals[i].nMergedIntervals;
			positions[i][1]=intervals[i].nMergedIntervals;
			positions[i][2]=intervals[i].firstPosition*positions[i][0];
			positions[i][3]=intervals[i].lastPosition*positions[i][1];
		}
		for (i=0; i<firstImplied; i++) {
			if (ReadA.sortedAlignments[i].mergedToInterval!=null) continue;
			if (IO.CONSISTENCY_CHECKS && ReadA.sortedAlignments[i].isLeftMaximal==1 && ReadA.sortedAlignments[i].isRightMaximal==1) {
				System.err.println("getAlignmentIntervals> ERROR 1");
				System.exit(1);
			}
			if (ReadA.sortedAlignments[i].isLeftMaximal!=1 && ReadA.sortedAlignments[i].isRightMaximal!=1) {
				// Avoiding nonmaximal alignments at this stage
				continue;
			}
			idI=ReadA.sortedAlignments[i].id;
			start=Alignments.alignments[idI][3];
			end=Alignments.alignments[idI][4];
			tmpInterval=keepMaximalAlignment(ReadA.sortedAlignments[i].isRightMaximal==1,ReadA.sortedAlignments[i],THRESHOLD);
			if (tmpInterval!=null) {
				ReadA.sortedAlignments[i].mergedToInterval=tmpInterval;
				tmpInterval.nMergedIntervals++;
				if (ReadA.sortedAlignments[i].isLeftMaximal==1) tmpInterval.nMaximalAlignmentsLeft++;
				if (ReadA.sortedAlignments[i].isRightMaximal==1) tmpInterval.nMaximalAlignmentsRight++;
				idJ=tmpInterval.id;
				if (Math.abs(start,tmpInterval.firstPosition)<=THRESHOLD) {
					positions[idJ][0]++;
					positions[idJ][2]+=start;
				}
				if (Math.abs(end,tmpInterval.lastPosition)<=THRESHOLD) {
					positions[idJ][1]++;
					positions[idJ][3]+=end;
				}
			}
		}
		for (i=0; i<=lastInterval; i++) {
			intervals[i].firstPosition=positions[i][2]/positions[i][0];
			intervals[i].lastPosition=positions[i][3]/positions[i][1];
		}		
		
		// Detecting intervals from alignments that are just left- or right-maximal.
		getAlignmentIntervals_weak(firstImplied,THRESHOLD,MAX_DISTANCE_THRESHOLD);
		
		// Marking all remaining alignments that could be explained by the detected
		// intervals. Remark: there might be both non-maximal and maximal alignments
		// remaining.
		for (i=0; i<=lastInterval; i++) {
			intervals[i].id=i;
			positions[i][0]=intervals[i].nMergedIntervals;
			positions[i][1]=intervals[i].nMergedIntervals;
			positions[i][2]=intervals[i].firstPosition*positions[i][0];
			positions[i][3]=intervals[i].lastPosition*positions[i][1];
		}
		for (i=0; i<firstImplied; i++) {
			if (ReadA.sortedAlignments[i].mergedToInterval!=null) continue;
			idI=ReadA.sortedAlignments[i].id;
			start=Alignments.alignments[idI][3];
			end=Alignments.alignments[idI][4];
			isMaximal=ReadA.sortedAlignments[i].isLeftMaximal==1 || ReadA.sortedAlignments[i].isRightMaximal==1;
			if (isMaximal) tmpInterval=keepMaximalAlignment(ReadA.sortedAlignments[i].isRightMaximal==1,ReadA.sortedAlignments[i],THRESHOLD);
			else tmpInterval=keepNonmaximalAlignment(ReadA.sortedAlignments[i],THRESHOLD);
			if (tmpInterval!=null) {
				ReadA.sortedAlignments[i].mergedToInterval=tmpInterval;
				tmpInterval.nMergedIntervals++;
				if (ReadA.sortedAlignments[i].isLeftMaximal==1) tmpInterval.nMaximalAlignmentsLeft++;
				if (ReadA.sortedAlignments[i].isRightMaximal==1) tmpInterval.nMaximalAlignmentsRight++;
				idJ=tmpInterval.id;
				if (Math.abs(start,tmpInterval.firstPosition)<=THRESHOLD) {
					positions[idJ][0]++;
					positions[idJ][2]+=start;
				}
				if (Math.abs(end,tmpInterval.lastPosition)<=THRESHOLD) {
					positions[idJ][1]++;
					positions[idJ][3]+=end;
				}
			}
		}
		for (i=0; i<=lastInterval; i++) {
			intervals[i].firstPosition=positions[i][2]/positions[i][0];
			intervals[i].lastPosition=positions[i][3]/positions[i][1];
		}	
		
		// Detecting intervals from the remaining alignments.
		// Remark: there might be both non-maximal and maximal alignments remaining.
		// Remark: we do a greedy merging without connected components as mentioned above.
		Alignment.order=Alignment.STARTA_ENDA;
		if (firstImplied>0) Arrays.sort(ReadA.sortedAlignments,0,firstImplied);
		for (i=0; i<firstImplied; i++) {
			if (ReadA.sortedAlignments[i].mergedToInterval!=null) continue;
			idI=ReadA.sortedAlignments[i].id;
			sumFirst=Alignments.alignments[idI][3]; 
			sumLast=Alignments.alignments[idI][4];
			min=Alignments.alignments[idI][3];
			max=Alignments.alignments[idI][4];
			isLeftMaximal=ReadA.sortedAlignments[i].isLeftMaximal;
			isRightMaximal=ReadA.sortedAlignments[i].isRightMaximal;
			isMaximal=(isLeftMaximal==1)&&(isRightMaximal==1);
			lastInterval++;
			intervals[lastInterval].id=0;
			intervals[lastInterval].isLeftMaximal=isLeftMaximal==1;
			intervals[lastInterval].isRightMaximal=isRightMaximal==1;
			intervals[lastInterval].isLongRepeat=false;
			ReadA.sortedAlignments[i].mergedToInterval=intervals[lastInterval];
			nMergedIntervals=1;
			nMaximalAlignmentsLeft=isLeftMaximal==1?1:0;
			nMaximalAlignmentsRight=isRightMaximal==1?1:0;
			for (j=i+1; j<firstImplied; j++) {
				idJ=ReadA.sortedAlignments[j].id;
				if (Alignments.alignments[idJ][3]>Alignments.alignments[idI][4]) break;
				if (ReadA.sortedAlignments[j].mergedToInterval!=null) continue;
				if ( ReadA.sortedAlignments[j].isLeftMaximal==isLeftMaximal && 
				     ReadA.sortedAlignments[j].isRightMaximal==isRightMaximal && 
				     Intervals.areApproximatelyIdentical(Alignments.alignments[idJ][3],Alignments.alignments[idJ][4],Alignments.alignments[idI][3],Alignments.alignments[idI][4])
				   ) {
					sumFirst+=Alignments.alignments[idJ][3];
					sumLast+=Alignments.alignments[idJ][4];
					min=Math.min(min,Alignments.alignments[idJ][3]);
					max=Math.max(max,Alignments.alignments[idJ][4]);
					nMergedIntervals++;
					if (isLeftMaximal==1) nMaximalAlignmentsLeft++;
					if (isRightMaximal==1) nMaximalAlignmentsRight++;
					ReadA.sortedAlignments[j].mergedToInterval=intervals[lastInterval];
				}
			}
			j=sumFirst/nMergedIntervals;
			intervals[lastInterval].firstPosition=min<j-THRESHOLD?min:j;
			j=sumLast/nMergedIntervals;
			intervals[lastInterval].lastPosition=max>j+THRESHOLD?max:j;
			intervals[lastInterval].nMergedIntervals=nMergedIntervals;
			intervals[lastInterval].nMaximalAlignmentsLeft=nMaximalAlignmentsLeft;
			intervals[lastInterval].nMaximalAlignmentsRight=nMaximalAlignmentsRight;
		}
		
		// Trimming only intervals whose boundaries are not close to a split (otherwise
		// such intervals might be discarded by $Factors.assignAlignmentIntervals()$).
		if (Split.order!=Split.POSITION) {
			// Ensuring the necessary order on $splits$.
			Split.order=Split.POSITION;
			if (Factors.lastSplit>0) Arrays.sort(Factors.splits,0,Factors.lastSplit+1);
		}
		for (i=0; i<=lastInterval; i++) {
			j=intervals[i].firstPosition;
			trimLeft=(!DenseSubstrings.inSubstringType(j)&&!PeriodicSubstrings.inShortPeriodInterval(j))||!Factors.isCloseToSplit(j,SPLITS_THRESHOLD);
			j=intervals[i].lastPosition;
			trimRight=(!DenseSubstrings.inSubstringType(j)&&!PeriodicSubstrings.inShortPeriodInterval(j))||!Factors.isCloseToSplit(j,SPLITS_THRESHOLD);
			intervals[i].trim(tmp,MIN_TRIMMED_LENGTH,trimLeft,trimRight);
		}
		Split.order=Split.UNSORTED;  // For downstream compatibility
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<firstImplied; i++) {
				if (ReadA.sortedAlignments[i].mergedToInterval==null) {
					System.err.println("getAlignmentIntervals> ERROR 3");
					System.exit(1);
				}
			}
		}
		discardAlignmentsAfterTrimming(0,firstImplied-1);
		AlignmentInterval.order=AlignmentInterval.UNSORTED;
		return firstImplied;
	}
	
	
	/**
	 * Builds intervals from longest paths of right-maximal but non-left-maximal 
	 * alignments that share the same end position, and such that the start positions of 
	 * every adjacent pair of alignments are at most $maxDistanceThreshold$ from each 
	 * other (and possibly identical). The same holds for sets of left-maximal but 
	 * non-right-maximal alignments.
	 * 
	 * Remark: the procedure does not merge two right-maximal alignments such that the 
	 * shorter if a suffix of the longer, and readA has low quality before or after the 
	 * beginning of the shorter.
	 *
	 * Remark: the procedure assumes that all alignments that are both left- and 
	 * right-maximal have been already merged to an interval.
	 *
	 * Remark: the procedure works just on $ReadA.sortedAlignments[0..firstImplied)$.
	 */
	private static final void getAlignmentIntervals_weak(int firstImplied, int identityThreshold, int maxDistanceThreshold) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_PATH_LENGTH = 3;  // Arbitrary
		final int MERGED_INTERVALS_THRESHOLD = 30*IO.coverage;  // Arbitrary
		int i, j;
		int startA1, endA1, startA2, endA2, nComponents, largestComponent, destination;
		AlignmentInterval tmpInterval;
		int[] componentSize;
		
		// Building the DAG of all right-maximal and non-left-maximal alignments, and
		// sorting it topologically.
		Alignment.order=Alignment.STARTA_ENDA;
		if (firstImplied>0) Arrays.sort(ReadA.sortedAlignments,0,firstImplied);
		clearNeighborMatrices(firstImplied);
		for (i=0; i<firstImplied; i++) {
			startA1=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			if (ReadA.sortedAlignments[i].mergedToInterval!=null || ReadA.sortedAlignments[i].isRightMaximal!=1) continue;
			endA1=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			for (j=i+1; j<firstImplied; j++) {
				startA2=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
				if (startA2>startA1+maxDistanceThreshold) break;
				if ( ReadA.sortedAlignments[j].mergedToInterval!=null || 
					 ReadA.sortedAlignments[j].isRightMaximal!=1 ||
					 (startA2>startA1+IDENTITY_THRESHOLD && (!Reads.isLeftMaximal(startA2,ReadA.id,true)||ReadA.sortedAlignments[j].lowQualityStart))
				   ) continue;
				endA2=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
				if (Math.abs(endA1,endA2)>identityThreshold) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(firstImplied,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stack);
		if (nComponents<=0 || nComponents>firstImplied) {
			IO.printCriticalErr("Error while detecting weak alignment intervals: wrong number of connected components.");
			System.exit(1);
		}
		componentSize = new int[nComponents];
		largestComponent=DAG.largestComponent(components,firstImplied,componentSize,nComponents);
		if (largestComponent>=MIN_PATH_LENGTH) {
			// Detecting all right-maximal alignment intervals
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,firstImplied);
			i=DAG.topologicalSort(firstImplied,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting weak alignment intervals: the weak suffix DAG of all alignments of read "+ReadA.id+" contains a cycle that involves node "+(i-1));
				System.exit(1);
			}
			i=0;
			while (i<firstImplied) {
				if (ReadA.sortedAlignments[i].mergedToInterval!=null || ReadA.sortedAlignments[i].isRightMaximal!=1 || componentSize[components[i]]<MIN_PATH_LENGTH) {
					i++;
					continue;
				}
				destination=getAlignmentInterval_weak(true,i,firstImplied,MIN_PATH_LENGTH,IDENTITY_THRESHOLD);
				if (destination==-1) {
					i++;
					continue;
				}
				lastInterval++;
				intervals[lastInterval].id=0;
				intervals[lastInterval].isLeftMaximal=false;
				intervals[lastInterval].isRightMaximal=true;
				intervals[lastInterval].firstPosition=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				markImpliedAlignments_weak(true,destination,intervals[lastInterval],IDENTITY_THRESHOLD);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("detected right-maximal alignment interval ["+intervals[lastInterval].firstPosition+".."+intervals[lastInterval].lastPosition+"]");
				i++;
			}
		}
		
		// Building the DAG of all left-maximal and non-right-maximal alignments, and
		// sorting it topologically.
		Alignment.order=Alignment.ENDA_REVERSE;
		if (firstImplied>0) Arrays.sort(ReadA.sortedAlignments,0,firstImplied);
		clearNeighborMatrices(firstImplied);
		for (i=0; i<firstImplied; i++) {
			if (ReadA.sortedAlignments[i].mergedToInterval!=null || ReadA.sortedAlignments[i].isLeftMaximal!=1) continue;
			startA1=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			endA1=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			for (j=i+1; j<firstImplied; j++) {
				endA2=Alignments.alignments[ReadA.sortedAlignments[j].id][4];
				if (endA2<endA1-maxDistanceThreshold) break;
				startA2=Alignments.alignments[ReadA.sortedAlignments[j].id][3];
				if (Math.abs(startA1,startA2)>identityThreshold) continue;
				if ( ReadA.sortedAlignments[j].mergedToInterval!=null || 
					 ReadA.sortedAlignments[j].isLeftMaximal!=1 ||
					 (endA2<endA1-IDENTITY_THRESHOLD && (!Reads.isRightMaximal(endA2,ReadA.id,true)||ReadA.sortedAlignments[j].lowQualityEnd))
				   ) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(firstImplied,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stack);
		if (nComponents<=0 || nComponents>firstImplied) {
			IO.printCriticalErr("Error while detecting weak alignment intervals: wrong number of connected components.");
			System.exit(1);
		}
		if (nComponents>componentSize.length) componentSize = new int[nComponents];
		largestComponent=DAG.largestComponent(components,firstImplied,componentSize,nComponents);
		if (largestComponent>=MIN_PATH_LENGTH) {
			// Detecting all left-maximal alignment intervals
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,firstImplied);
			i=DAG.topologicalSort(firstImplied,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting weak alignment intervals: the weak prefix DAG of all alignments of read "+ReadA.id+" contains a cycle that involves node "+(i-1));
				System.exit(1);
			}
			i=0;
			while (i<firstImplied) {
				if (ReadA.sortedAlignments[i].mergedToInterval!=null || ReadA.sortedAlignments[i].isLeftMaximal!=1 || componentSize[components[i]]<MIN_PATH_LENGTH) {
					i++;
					continue;
				}
				destination=getAlignmentInterval_weak(false,i,firstImplied,MIN_PATH_LENGTH,IDENTITY_THRESHOLD);
				if (destination==-1) {
					i++;
					continue;
				}
				lastInterval++;
				intervals[lastInterval].id=0;
				intervals[lastInterval].isLeftMaximal=true;
				intervals[lastInterval].isRightMaximal=false;
				intervals[lastInterval].lastPosition=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
				markImpliedAlignments_weak(false,destination,intervals[lastInterval],IDENTITY_THRESHOLD);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("detected left-maximal alignment interval ["+intervals[lastInterval].firstPosition+".."+intervals[lastInterval].lastPosition+"] MIN_PATH_LENGTH="+MIN_PATH_LENGTH);
				i++;
			}
		}
		
		// Removing intervals that are implied by other intervals
		for (i=0; i<=lastInterval; i++) intervals[i].representative=intervals[i];
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].representative!=intervals[i]) continue;
			startA1=intervals[i].firstPosition; endA1=intervals[i].lastPosition;
			for (j=i-1; j>=0; j--) {
				if (intervals[j].representative!=intervals[j]) continue;
				startA2=intervals[j].firstPosition; endA2=intervals[j].lastPosition;
				if ( intervals[j].nMergedIntervals<MERGED_INTERVALS_THRESHOLD &&
					 ( (startA2>startA1+IDENTITY_THRESHOLD && Math.abs(endA2,endA1)<=IDENTITY_THRESHOLD && !intervals[j].isLeftMaximal) ||
					   (Math.abs(startA2,startA1)<=IDENTITY_THRESHOLD && endA2<endA1-IDENTITY_THRESHOLD && !intervals[j].isRightMaximal)
					 )
				   ) {
					   intervals[j].representative=intervals[i];
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("discarded interval "+intervals[j]+" because implied by "+intervals[i]);
				   }
			}
			for (j=i+1; j<=lastInterval; j++) {
				if (intervals[j].representative!=intervals[j]) continue;
				startA2=intervals[j].firstPosition; endA2=intervals[j].lastPosition;
				if ( intervals[j].nMergedIntervals<MERGED_INTERVALS_THRESHOLD &&
					 ( (startA2>startA1+IDENTITY_THRESHOLD && Math.abs(endA2,endA1)<=IDENTITY_THRESHOLD && !intervals[j].isLeftMaximal) ||
					   (Math.abs(startA2,startA1)<=IDENTITY_THRESHOLD && endA2<endA1-IDENTITY_THRESHOLD && !intervals[j].isRightMaximal)
					 )
				   ) {
					   intervals[j].representative=intervals[i];
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("discarded interval "+intervals[j]+" because implied by "+intervals[i]);
				   }
			}
		}
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].representative!=intervals[i]) continue;
			j++;
			if (j!=i) {
				tmpInterval=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=tmpInterval;
			}
		}
		lastInterval=j;
		
		// Updating pointers from alignments
		for (i=0; i<firstImplied; i++) {
			tmpInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (tmpInterval!=null) ReadA.sortedAlignments[i].mergedToInterval=tmpInterval.closure();
		}
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
	
	
	/**
	 * Detects a right-maximal but non-left-maximal alignment interval (if 
	 * $rightOrLeft=TRUE$), or a left-maximal but non-right-maximal interval (if 
	 * $rightOrLeft=FALSE$) that starts from alignment $ReadA.sortedAlignments[source]$ in
	 * the topologically sorted DAG.
	 *
	 * @return the position in $ReadA.sortedAlignments[0..nVertices-1]$ of an interval 
	 * that maximizes the distance from $source$ in the DAG; -1 if no path from $source$ 
	 * contains at least $minPathLength$ vertices.
	 */
	private static final int getAlignmentInterval_weak(boolean rightOrLeft, int source, int nVertices, int minPathLength, int distanceThreshold) {
		int i, j, k, s;
		int sourceComponent, vertex, selectedAlignment, distance, longestDistance;
		int start, end, leftmostStart, rightmostEnd;
		int lastSortedVertex;
		final int sourceStartA = ReadA.sortedAlignments[source].startA();
		final int sourceEndA = ReadA.sortedAlignments[source].endA();

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
			if ( (rightOrLeft && Math.abs(ReadA.sortedAlignments[vertex].endA(),sourceEndA)>distanceThreshold) ||
				 (!rightOrLeft && Math.abs(ReadA.sortedAlignments[vertex].startA(),sourceStartA)>distanceThreshold)
			   ) continue;
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

		// Computing an alignment whose distance from $source$ in the DAG is maximum, and
		// whose ending position is the rightmost (respectively, whose starting position
		// is the leftmost).
		longestDistance=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]>longestDistance) longestDistance=distances[vertex];
		}
		if (longestDistance+1<minPathLength) return -1;
		selectedAlignment=-1; leftmostStart=Math.POSITIVE_INFINITY; rightmostEnd=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]<longestDistance) continue;
			start=Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			if (start<leftmostStart) {
				leftmostStart=start;
				if (!rightOrLeft) selectedAlignment=vertex;
			}
			end=Alignments.alignments[ReadA.sortedAlignments[vertex].id][4];
			if (end>rightmostEnd) {
				rightmostEnd=end;
				if (rightOrLeft) selectedAlignment=vertex;
			}
		}
		return selectedAlignment;
	}
	
	
	/**
	 * Sets $mergedToInterval=interval$ for all alignments in a longest path to 
	 * $destination$, using $predecessors$. Sets the last position (if right-maximal) or 
	 * the first position (if left-maximal) of $interval$ to the average over all 
	 * intervals in the longest path. Sets also the fields $nMergedIntervals$ and 
	 * $isLongRepeat$ of $interval$.
	 *
	 * @param rightOrLeft TRUE=right-maximal, FALSE=left-maximal;
	 */
	private static final void markImpliedAlignments_weak(boolean rightOrLeft, int destination, AlignmentInterval interval, int identityThreshold) {
		int vertex, nAlignments;
		int p, boundary;
		double sum;

		interval.isLongRepeat=false;
		sum=0.0; nAlignments=0; vertex=destination;
		boundary=rightOrLeft?-1:Math.POSITIVE_INFINITY;
		while (vertex!=-1) {
			ReadA.sortedAlignments[vertex].mergedToInterval=interval;
			p=rightOrLeft?Alignments.alignments[ReadA.sortedAlignments[vertex].id][4]:Alignments.alignments[ReadA.sortedAlignments[vertex].id][3];
			sum+=p;
			if (rightOrLeft && p>boundary) boundary=p;
			if (!rightOrLeft && p<boundary) boundary=p;
			nAlignments++;
			if ( (rightOrLeft && Alignments.alignments[ReadA.sortedAlignments[vertex].id][3]>interval.firstPosition+identityThreshold) ||
				 (!rightOrLeft && Alignments.alignments[ReadA.sortedAlignments[vertex].id][4]<interval.lastPosition-identityThreshold)
			   ) interval.isLongRepeat=true;
			vertex=predecessors[vertex];
		}
		p=(int)(sum/nAlignments);
		if (rightOrLeft) interval.lastPosition=Math.abs(boundary,p)<=identityThreshold?p:boundary;
		else interval.firstPosition=Math.abs(boundary,p)<=identityThreshold?p:boundary;
		interval.nMergedIntervals=nAlignments;
		if (rightOrLeft) interval.nMaximalAlignmentsRight=nAlignments;
		else interval.nMaximalAlignmentsLeft=nAlignments;
	}
	
	
	/**
	 * The right-maximal but non-left-maximal alignment $alignment$ should be kept iff no 
	 * existing right-maximal interval ends close to its end, and approximately contains 
	 * it. A symmetrical argument holds for a left-maximal but non-right-maximal 
	 * alignment.
	 *
	 * Remark: the procedure does not assign a right-maximal alignment to an interval, if 
	 * $leftEnd$ is far from the beginning of the interval and readA has low quality 
	 * before or after $leftEnd$.
	 *
	 * Remark: the procedure scans all current alignment intervals, since it does not 
	 * assume they are sorted.
	 *
	 * @patram rightOrLeft TRUE=right-maximal, FALSE=left-maximal;
	 * @return NULL if the alignment should be kept; otherwise, a shortest right-maximal 
	 * (respectively, left-maximal) alignment interval that contains the given alignment.
	 */
	private static final AlignmentInterval keepMaximalAlignment(boolean rightOrLeft, Alignment alignment, int identityThreshold) {
		boolean isMaximal;
		int i;
		int length, minLength, leftEnd, rightEnd;
		AlignmentInterval minInterval;
		
		leftEnd=alignment.startA(); rightEnd=alignment.endA();
		minInterval=null; minLength=Math.POSITIVE_INFINITY;
		if (rightOrLeft) {
			isMaximal=Reads.isLeftMaximal(leftEnd,ReadA.id,true);
			for (i=0; i<=lastInterval; i++) {
				if ( !intervals[i].isRightMaximal || 
				     Math.abs(intervals[i].lastPosition,rightEnd)>identityThreshold || 
					 !( Intervals.isApproximatelyContained(leftEnd,rightEnd,intervals[i].firstPosition,intervals[i].lastPosition) ||
						Intervals.areApproximatelyIdentical(leftEnd,rightEnd,intervals[i].firstPosition,intervals[i].lastPosition)
					 ) ||
					 (leftEnd>intervals[i].firstPosition+identityThreshold && (!isMaximal||alignment.lowQualityStart))
				   ) continue;
				length=intervals[i].lastPosition-intervals[i].firstPosition+1;
				if (length<minLength) {
					minLength=length;
					minInterval=intervals[i];
				}
			}
		}
		else {
			isMaximal=Reads.isRightMaximal(rightEnd,ReadA.id,true);
			for (i=0; i<=lastInterval; i++) {
				if ( !intervals[i].isLeftMaximal || 
				     Math.abs(intervals[i].firstPosition,leftEnd)>identityThreshold || 
					 !( Intervals.isApproximatelyContained(leftEnd,rightEnd,intervals[i].firstPosition,intervals[i].lastPosition) ||
						Intervals.areApproximatelyIdentical(leftEnd,rightEnd,intervals[i].firstPosition,intervals[i].lastPosition)
					 ) ||
					 (rightEnd<intervals[i].lastPosition-identityThreshold && (!isMaximal||alignment.lowQualityEnd))
				   ) continue;
				length=intervals[i].lastPosition-intervals[i].firstPosition+1;
				if (length<minLength) {
					minLength=length;
					minInterval=intervals[i];
				}
			}
		}
		return minInterval;
	}
	
	
	/**
	 * The alignment $alignment$, which is neither left- nor right-maximal, should be kept
	 * iff no interval that is left-maximal, right-maximal, or both, approximately 
	 * contains it.
	 *
	 * Remark: the procedure does not assign an alignment to an interval, if $leftEnd$ is 
	 * far from the beginning of the interval and readA has low quality before or after 
	 * $leftEnd$ (a symmetric argument holds for $rightEnd$).
	 *
	 * Remark: the procedure assumes $intervals[0..lastInterval]$ to contain just 
	 * intervals that are left-maximal, right-maximal, or both.
	 *
	 * Remark: the procedure scans all current alignment intervals, since it does not 
	 * assume they are sorted.
	 *
	 * @return NULL if the alignment should be kept; otherwise, a shortest left-maximal,
	 * right-maximal, or both, alignment interval, that contains the given alignment.
	 */
	private static final AlignmentInterval keepNonmaximalAlignment(Alignment alignment, int identityThreshold) {
		int i;
		int length, minLength, leftEnd, rightEnd;
		AlignmentInterval minInterval;
		
		leftEnd=alignment.startA(); rightEnd=alignment.endA();
		minInterval=null; minLength=Math.POSITIVE_INFINITY;
		for (i=0; i<=lastInterval; i++) {
			if ( !Intervals.isApproximatelyContained(leftEnd,rightEnd,intervals[i].firstPosition,intervals[i].lastPosition) && 
				 !Intervals.areApproximatelyIdentical(leftEnd,rightEnd,intervals[i].firstPosition,intervals[i].lastPosition)
			   ) continue;
			if ( (leftEnd>intervals[i].firstPosition+identityThreshold && (!Reads.isLeftMaximal(leftEnd,ReadA.id,true)||alignment.lowQualityStart)) || 
				 (rightEnd<intervals[i].lastPosition-identityThreshold && (!Reads.isRightMaximal(rightEnd,ReadA.id,true)||alignment.lowQualityEnd))
			   ) continue;
			length=intervals[i].lastPosition-intervals[i].firstPosition+1;
			if (length<minLength) {
				minLength=length;
				minInterval=intervals[i];
			}
		}
		return minInterval;
	}
	
	
	/**
	 * Removes alignment intervals that are: (1) too similar to a dense substring or to a
	 * periodic substring interval; (2) just right-maximal, and such that there is a dense
	 * or periodic substring that contains the interval and ends close to the end of the
	 * interval, and the interval is the merge of too few alignments; (3) neither left- 
	 * nor right-maximal, and such that a dense or periodic substring contains the 
	 * interval, or the interval is the merge of too few alignments; (4) right-maximal and
	 * the prefix of a prefix substring, but they are not the merge of enough alignments 
	 * (a symmetrical criterion holds for suffix substrings); (5) straddling periodic 
	 * intervals.
	 * The last condition in point (3) was suggested by Gene Myers, and it is useful since
	 * such non-contained intervals might be strings that just happen to occur multiple 
	 * times in the dataset because of randomness in how the subset of all reads was 
	 * sampled, i.e. they might not be repeats of the genome. Point (4) is useful, since
	 * such intervals might be generated by prefix/suffix biases detected by the dense 
	 * substrings pipeline but too weak in practice.
	 *
	 * Remark: note that using average coverage in the latter condition above is wrong, 
	 * since a non-maximal interval might contain any configuration of repeats.
	 *
	 * Remark: an alignment interval might be identical to e.g. a long-period periodic
	 * interval that coincides with one occurrence of the period, and such that this 
	 * single occurrence appears also in another read.
	 *
	 * Remark: the procedure does not assign an alignment interval to a dense or periodic
	 * substring interval, if its ends are adjacent to low-quality regions in readA.
	 *
	 * Remark: the procedure takes into account that some periodic substrings might not
	 * point to any periodic substring interval.
	 *
	 * Remark: the procedure has different meanings depending on when it is called, since 
	 * $AlignmentInterval.is*Maximal$ has different meanings at different steps of 
	 * $detect()$.
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ to be sorted by $startA$
	 * and $PeriodicSubstrings.intervals$ to be sorted by $firstPosition$. The procedure 
	 * also assumes that $PeriodicSubstrings.{left,right}Splits$ contains an updated 
	 * version of all peaks inside maximal ranges of overlapping short-period intervals.
	 *
	 * @param firstImplied the procedure assumes $ReadA.sortedAlignments$ to be sorted by 
	 * $IMPL_STARTA$, and all and only the alignments implied by a dense or periodic 
	 * substring to be located starting from $firstImplied$; this might not be true any 
	 * more after the procedure completes;
	 * @return TRUE iff $ReadA.sortedAlignments$ is still sorted by $IMPL_STARTA$ at the 
	 * end of the procedure.
	 */
	private static final boolean cleanAlignmentIntervals(int firstImplied) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int SPLIT_THRESHOLD = IO.quantum<<1;  // Arbitrary
		boolean denseOrPeriodic, foundDense, foundPeriodic, isSorted;
		int i, j, k, firstJForNext;
		int start1, start2, start3, end1, end2, end3, length, minLength;
		int previousLastInterval, rangeFirst, rangeLast, rangeFrom, rangeTo;
		int leftSplit, rightSplit;
		AlignmentInterval tmp;
		PeriodicSubstringInterval periodicSubstringInterval;
		
		// Ensuring the necessary order in $intervals$
		if (lastInterval==-1) return true;
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		for (i=0; i<=lastInterval; i++) {
			intervals[i].discarded=false;
			intervals[i].denseSubstring=null;
			intervals[i].inPeriodicSubstringInterval=null;
			intervals[i].flag2=true;  // TRUE = If $intervals[i]$ gets discarded, we should try to reassign, to another alignment interval, the alignments assigned to $intervals[i]$.
		}

		// Comparing intervals to dense substrings
		foundDense=cleanAlignmentIntervals_dense(DISTANCE_THRESHOLD);

		// Comparing intervals to periodic substring intervals
		foundPeriodic=false;
		if (PeriodicSubstrings.lastInterval>=0) {
			i=0; j=0; firstJForNext=-1; 
			periodicSubstringInterval=null; minLength=Math.POSITIVE_INFINITY;
			while (i<=lastInterval) {
				if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=intervals[i].lastPosition) {
					if (periodicSubstringInterval!=null) {
						intervals[i].inPeriodicSubstringInterval=periodicSubstringInterval;
						foundPeriodic=true;
						minLength=Math.POSITIVE_INFINITY; periodicSubstringInterval=null;
					}
					i++;
					if (firstJForNext!=-1) j=firstJForNext;
					firstJForNext=-1;
					continue;
				}
				if (PeriodicSubstrings.intervals[j].lastPosition<=intervals[i].firstPosition) {
					j++;
					continue;
				}
				if (i<lastInterval && firstJForNext==-1 && PeriodicSubstrings.intervals[j].lastPosition>=intervals[i+1].firstPosition) firstJForNext=j;
				length=PeriodicSubstrings.intervals[j].lastPosition-PeriodicSubstrings.intervals[j].firstPosition;
				leftSplit=PeriodicSubstrings.position2split(intervals[i].firstPosition,PeriodicSubstrings.leftSplits,PeriodicSubstrings.lastLeftSplit,SPLIT_THRESHOLD);
				rightSplit=PeriodicSubstrings.position2split(intervals[i].lastPosition,PeriodicSubstrings.rightSplits,PeriodicSubstrings.lastRightSplit,SPLIT_THRESHOLD);
				if ( Intervals.areApproximatelyIdentical_lowQuality(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,intervals[i].firstPosition,intervals[i].lastPosition,ReadA.id) ||
					 ( (leftSplit!=-1 && PeriodicSubstrings.position2split(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.leftSplits,PeriodicSubstrings.lastLeftSplit,SPLIT_THRESHOLD)==leftSplit) &&
					   (rightSplit!=-1 && PeriodicSubstrings.position2split(PeriodicSubstrings.intervals[j].lastPosition,PeriodicSubstrings.rightSplits,PeriodicSubstrings.lastRightSplit,SPLIT_THRESHOLD)==rightSplit)
					 )
				   ) {
					if (length<minLength) {
						minLength=length;
						periodicSubstringInterval=PeriodicSubstrings.intervals[j];
					}
				}
				else if ( Intervals.isApproximatelyContained_lowQuality(intervals[i].firstPosition,intervals[i].lastPosition,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.id) &&
					      !( (intervals[i].firstPosition>PeriodicSubstrings.intervals[j].firstPosition+DISTANCE_THRESHOLD && !Reads.isLeftMaximal(intervals[i].firstPosition,ReadA.id,true)) &&
							 (intervals[i].lastPosition<PeriodicSubstrings.intervals[j].lastPosition-DISTANCE_THRESHOLD && !Reads.isRightMaximal(intervals[i].lastPosition,ReadA.id,true)) &&
							 PeriodicSubstrings.intervals[j].nMergedIntervals>=MIN_INTERVALS_IN_MAXIMAL
						  )
				        ) {
					if ( (!intervals[i].isLeftMaximal || (leftSplit!=-1 && (PeriodicSubstrings.leftSplitLowQuality[3*leftSplit]||PeriodicSubstrings.leftSplitFlags[3*leftSplit]))) && 
					     (!intervals[i].isRightMaximal || (rightSplit!=-1 && (PeriodicSubstrings.rightSplitLowQuality[3*rightSplit]||PeriodicSubstrings.rightSplitFlags[3*rightSplit])))
					   ) {
						if (length<minLength) {
							minLength=length;
							periodicSubstringInterval=PeriodicSubstrings.intervals[j];
						}
					}
					else if ( (!intervals[i].isLeftMaximal || (leftSplit!=-1 && PeriodicSubstrings.leftSplitLowQuality[3*leftSplit])) && 
						      (intervals[i].isRightMaximal && (rightSplit==-1 || !PeriodicSubstrings.rightSplitLowQuality[3*rightSplit])) && 
					          ( Math.abs(intervals[i].lastPosition,PeriodicSubstrings.intervals[j].lastPosition)<=DISTANCE_THRESHOLD ||
								( !PeriodicSubstrings.intervals[j].hasLongPeriod &&
								  intervals[i].nMaximalAlignmentsRight<MIN_INTERVALS_IN_MAXIMAL
							    ) ||
								( PeriodicSubstrings.intervals[j].hasLongPeriod &&
								  cleanAlignmentIntervals_discardWithLongPeriod(intervals[i],PeriodicSubstrings.intervals[j],false,true,MIN_INTERVALS_IN_MAXIMAL,MIN_INTERVALS_IN_NONMAXIMAL)
								)
							  )
					        ) {
						if (length<minLength) {
							minLength=length;
							periodicSubstringInterval=PeriodicSubstrings.intervals[j];
						}
					}
					else if ( (intervals[i].isLeftMaximal && (leftSplit==-1 || !PeriodicSubstrings.leftSplitLowQuality[3*leftSplit])) && 
						      (!intervals[i].isRightMaximal || (rightSplit!=-1 && PeriodicSubstrings.rightSplitLowQuality[3*rightSplit])) &&
					          ( Math.abs(intervals[i].firstPosition,PeriodicSubstrings.intervals[j].firstPosition)<=DISTANCE_THRESHOLD ||
								( !PeriodicSubstrings.intervals[j].hasLongPeriod &&
								  intervals[i].nMaximalAlignmentsLeft<MIN_INTERVALS_IN_MAXIMAL
							    ) ||
								( PeriodicSubstrings.intervals[j].hasLongPeriod &&
								  cleanAlignmentIntervals_discardWithLongPeriod(intervals[i],PeriodicSubstrings.intervals[j],true,false,MIN_INTERVALS_IN_MAXIMAL,MIN_INTERVALS_IN_NONMAXIMAL)
								)
							  )
					        ) {
						if (length<minLength) {
							minLength=length;
							periodicSubstringInterval=PeriodicSubstrings.intervals[j];
						}
					}
					else if ( (intervals[i].isLeftMaximal && (leftSplit==-1 || !PeriodicSubstrings.leftSplitLowQuality[3*leftSplit])) && 
						      (intervals[i].isRightMaximal && (rightSplit==-1 || !PeriodicSubstrings.rightSplitLowQuality[3*rightSplit])) &&
					  		  ( ( !PeriodicSubstrings.intervals[j].hasLongPeriod &&
								  intervals[i].nMaximalAlignmentsLeft<MIN_INTERVALS_IN_MAXIMAL && intervals[i].nMaximalAlignmentsRight<MIN_INTERVALS_IN_MAXIMAL
							    ) || 
								( PeriodicSubstrings.intervals[j].hasLongPeriod &&
								  cleanAlignmentIntervals_discardWithLongPeriod(intervals[i],PeriodicSubstrings.intervals[j],true,true,MIN_INTERVALS_IN_MAXIMAL,MIN_INTERVALS_IN_NONMAXIMAL)
								)
					          )
					        ) {
						if (length<minLength) {
							minLength=length;
							periodicSubstringInterval=PeriodicSubstrings.intervals[j];
						}
					}
				}
				j++;
			}
			
			// Discarding alignment intervals with the same peaks as a periodic interval
			i=0;
			leftSplit=PeriodicSubstrings.position2split(intervals[i].firstPosition,PeriodicSubstrings.leftSplits,PeriodicSubstrings.lastLeftSplit,SPLIT_THRESHOLD);
			rightSplit=PeriodicSubstrings.position2split(intervals[i].lastPosition,PeriodicSubstrings.rightSplits,PeriodicSubstrings.lastRightSplit,SPLIT_THRESHOLD);
			j=0; firstJForNext=-1; periodicSubstringInterval=null; minLength=Math.POSITIVE_INFINITY;
			while (i<=lastInterval) {
				if ( intervals[i].inPeriodicSubstringInterval!=null || leftSplit==-1 || rightSplit==-1 ||
				     j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=intervals[i].lastPosition-IDENTITY_THRESHOLD
				   ) {
					if (periodicSubstringInterval!=null) {
						intervals[i].inPeriodicSubstringInterval=periodicSubstringInterval;
						foundPeriodic=true;
						periodicSubstringInterval=null; minLength=Math.POSITIVE_INFINITY;
					}
					i++;
					if (i<=lastInterval) {
						leftSplit=PeriodicSubstrings.position2split(intervals[i].firstPosition,PeriodicSubstrings.leftSplits,PeriodicSubstrings.lastLeftSplit,SPLIT_THRESHOLD);
						rightSplit=PeriodicSubstrings.position2split(intervals[i].lastPosition,PeriodicSubstrings.rightSplits,PeriodicSubstrings.lastRightSplit,SPLIT_THRESHOLD);
					}
					if (firstJForNext!=-1) j=firstJForNext;
					firstJForNext=-1;
					continue;
				}
				if (PeriodicSubstrings.intervals[j].lastPosition<=intervals[i].firstPosition) {
					j++;
					continue;
				}
				if (firstJForNext==-1 && i<lastInterval && PeriodicSubstrings.intervals[j].lastPosition>=intervals[i+1].firstPosition) firstJForNext=j;
				if ( PeriodicSubstrings.position2split(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.leftSplits,PeriodicSubstrings.lastLeftSplit,SPLIT_THRESHOLD)==leftSplit && 
					 PeriodicSubstrings.position2split(PeriodicSubstrings.intervals[j].lastPosition,PeriodicSubstrings.rightSplits,PeriodicSubstrings.lastRightSplit,SPLIT_THRESHOLD)==rightSplit
				   ) {
					length=PeriodicSubstrings.intervals[j].length();
					if (length<minLength) {
						minLength=length;
						periodicSubstringInterval=PeriodicSubstrings.intervals[j];
					}
				}
				j++;
			}
			
			// Resetting pointers from long-period clones to the real long-period
			// intervals.
			for (i=0; i<=lastInterval; i++) {
				if (intervals[i].inPeriodicSubstringInterval==null || !intervals[i].inPeriodicSubstringInterval.hasLongPeriod) continue;
				j=Arrays.binarySearch(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1,intervals[i].inPeriodicSubstringInterval);
				if (j<0) {
					System.err.println("cleanAlignmentIntervals> ERROR: long-period interval not found.");
					System.exit(1);
				}
				foundPeriodic=false;
				for (k=j; k>=0; k--) {
					if (PeriodicSubstrings.longPeriodIntervals[k].equals(intervals[i].inPeriodicSubstringInterval)) {
						intervals[i].inPeriodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[k];
						foundPeriodic=true;
						break;
					}
				}
				if (foundPeriodic) continue;
				for (k=j+1; k<=PeriodicSubstrings.lastLongPeriodInterval; k++) {
					if (PeriodicSubstrings.longPeriodIntervals[k].equals(intervals[i].inPeriodicSubstringInterval)) {
						intervals[i].inPeriodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[k];
						break;
					}
				}
			}
		}
		
		// Compacting intervals
		previousLastInterval=lastInterval;
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].denseSubstring!=null) {
				intervals[i].discarded=true;
				if (intervals[i].denseSubstring.substringReplication) intervals[i].flag2=false;
				continue;
			}
			if (intervals[i].inPeriodicSubstringInterval!=null) {
				intervals[i].discarded=true;
				intervals[i].flag2=false;
				continue;
			}
			leftSplit=PeriodicSubstrings.position2split(intervals[i].firstPosition,PeriodicSubstrings.leftSplits,PeriodicSubstrings.lastLeftSplit,SPLIT_THRESHOLD);
			rightSplit=PeriodicSubstrings.position2split(intervals[i].lastPosition,PeriodicSubstrings.rightSplits,PeriodicSubstrings.lastRightSplit,SPLIT_THRESHOLD);
			if ( (!intervals[i].isLeftMaximal || (leftSplit!=-1 && PeriodicSubstrings.leftSplitLowQuality[3*leftSplit])) && 
				 (!intervals[i].isRightMaximal || (rightSplit!=-1 && PeriodicSubstrings.rightSplitLowQuality[3*rightSplit])) && 
			     intervals[i].nMergedIntervals<MIN_INTERVALS_IN_NONMAXIMAL
		       ) {
				intervals[i].discarded=true;
				continue;
			}
			j++;
			if (j!=i) {
				tmp=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=tmp;
			}
		}
		lastInterval=j;
		if (lastInterval==previousLastInterval) return true;
		
		// Updating pointers from alignments
		isSorted=discardIntervals_reassignAlignments(firstImplied,3,true,true);
		if (!isSorted) Alignment.order=Alignment.UNSORTED;		
		
		// Discarding alignment intervals that straddle periodic intervals
		if (PeriodicSubstrings.lastInterval>=0) {
			rangeFrom=0;
			rangeFirst=PeriodicSubstrings.intervals[0].firstPosition;
			rangeLast=PeriodicSubstrings.intervals[0].lastPosition;
			j=0; firstJForNext=-1;
			for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
				if (PeriodicSubstrings.intervals[i].firstPosition<=rangeLast+DISTANCE_THRESHOLD) {
					rangeLast=Math.max(rangeLast,PeriodicSubstrings.intervals[i].lastPosition);
					continue;
				}
				if (firstJForNext!=-1) j=firstJForNext;
				firstJForNext=cleanAlignmentIntervals_straddling(rangeFirst,rangeLast,rangeFrom,i-1,j);
				rangeFrom=i;
				rangeFirst=PeriodicSubstrings.intervals[i].firstPosition;
				rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
			}
			if (firstJForNext!=-1) j=firstJForNext;
			cleanAlignmentIntervals_straddling(rangeFirst,rangeLast,rangeFrom,i-1,j);
		}
		
		// Discarding alignment intervals that form prefix/suffix paths inside long- and
		// short-period intervals.
		cleanAlignmentIntervals_longPeriodEnds(Math.POSITIVE_INFINITY/*Aggressive*/);
		cleanAlignmentIntervals_shortPeriodEnds(Math.POSITIVE_INFINITY/*Aggressive*/);
		
		// Compacting intervals
		previousLastInterval=lastInterval;
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].discarded) continue;
			j++;
			if (j!=i) {
				tmp=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=tmp;
			}
		}
		lastInterval=j;
		
		// Updating pointers from alignments
		if (lastInterval<previousLastInterval) discardIntervals_reassignAlignments(firstImplied,0,true,true);

		return isSorted;
	}
	
	
	/**
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ and $intervals$ to be
	 * sorted by first position.
	 */
	public static final boolean cleanAlignmentIntervals_dense(int distanceThreshold) {
		final double MERGED_INTERVALS_THRESHOLD = 2.0;  // Arbitrary
		boolean foundDense;
		int i, j;
		int firstJForNext, length, minLength;
		DenseSubstring denseSubstring;
		if (DenseSubstrings.lastSubstring<0) return false;
		
		foundDense=false;
		i=0; j=0; firstJForNext=-1; 
		minLength=Math.POSITIVE_INFINITY; denseSubstring=null;
		while (i<=lastInterval) {			
			if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>=intervals[i].lastPosition) {
				if (denseSubstring!=null) {
					intervals[i].denseSubstring=denseSubstring;
					foundDense=true;
					minLength=Math.POSITIVE_INFINITY; denseSubstring=null;						
				}
				i++;
				if (firstJForNext!=-1) j=firstJForNext;
				firstJForNext=-1;
				continue;
			}
			if (DenseSubstrings.substrings[j].endA<=intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (i<lastInterval && firstJForNext==-1 && DenseSubstrings.substrings[j].endA>=intervals[i+1].firstPosition) firstJForNext=j;
			length=DenseSubstrings.substrings[j].endA-DenseSubstrings.substrings[j].startA;
			if ( Intervals.areApproximatelyIdentical(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,intervals[i].firstPosition,intervals[i].lastPosition) ||
				 ( Intervals.areApproximatelyIdentical_lowQuality(DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,intervals[i].firstPosition,intervals[i].lastPosition,ReadA.id) &&
				   intervals[i].nMergedIntervals<=DenseSubstrings.substrings[j].nMaximalAlignments*MERGED_INTERVALS_THRESHOLD
				 )
			   ) {
				if (length<minLength) {
					minLength=length;
					denseSubstring=DenseSubstrings.substrings[j];
				}
			}
			else if ( Intervals.isApproximatelyContained_lowQuality(intervals[i].firstPosition,intervals[i].lastPosition,DenseSubstrings.substrings[j].startA,DenseSubstrings.substrings[j].endA,ReadA.id) &&
				      !( (intervals[i].firstPosition>DenseSubstrings.substrings[j].startA+distanceThreshold && !Reads.isLeftMaximal(intervals[i].firstPosition,ReadA.id,true)) &&
						 (intervals[i].lastPosition<DenseSubstrings.substrings[j].endA-distanceThreshold && !Reads.isRightMaximal(intervals[i].lastPosition,ReadA.id,true))
					  )
			        ) {
				if (!intervals[i].isLeftMaximal && !intervals[i].isRightMaximal && intervals[i].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX) {
					if (length<minLength) {
						minLength=length;
						denseSubstring=DenseSubstrings.substrings[j];
					}
				}
				else if ( !intervals[i].isLeftMaximal && intervals[i].isRightMaximal && 
				          ( ( Math.abs(intervals[i].lastPosition,DenseSubstrings.substrings[j].endA)<=distanceThreshold &&
							  intervals[i].nMaximalAlignmentsLeft<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX && 
							  intervals[i].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
							) ||
							( (DenseSubstrings.substrings[j].prefixReplication || DenseSubstrings.substrings[j].singleDeletionReplication || DenseSubstrings.substrings[j].substringReplication) &&
							  intervals[i].nMaximalAlignmentsRight<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX &&
							  intervals[i].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
						    )
						  )
				        ) {
					if (length<minLength) {
						minLength=length;
						denseSubstring=DenseSubstrings.substrings[j];
					}
				}
				else if ( intervals[i].isLeftMaximal && !intervals[i].isRightMaximal &&
				          ( ( Math.abs(intervals[i].firstPosition,DenseSubstrings.substrings[j].startA)<=distanceThreshold &&
							  intervals[i].nMaximalAlignmentsRight<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX && 
							  intervals[i].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
							) ||
							( (DenseSubstrings.substrings[j].suffixReplication || DenseSubstrings.substrings[j].singleDeletionReplication || DenseSubstrings.substrings[j].substringReplication) &&
							  intervals[i].nMaximalAlignmentsLeft<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX &&
							  intervals[i].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
							)
						  )
				        ) {
					if (length<minLength) {
						minLength=length;
						denseSubstring=DenseSubstrings.substrings[j];
					}
				}
				else if (DenseSubstrings.substrings[j].substringReplication && !intervals[i].checkAlignments(MIN_INTERVALS_IN_MAXIMAL,MIN_ALIGNMENTS_LARGE)) {
					if (length<minLength) {
						minLength=length;
						denseSubstring=DenseSubstrings.substrings[j];
					}
				}
				else if ( ( (DenseSubstrings.substrings[j].prefixReplication || DenseSubstrings.substrings[j].singleDeletionReplication) && 
					        intervals[i].isRightMaximal && 
						    Math.abs(intervals[i].firstPosition,DenseSubstrings.substrings[j].startA)<=distanceThreshold &&
						    intervals[i].nMaximalAlignmentsRight<MIN_INTERVALS_IN_MAXIMAL &&
							intervals[i].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
						  ) ||
						  ( (DenseSubstrings.substrings[j].suffixReplication || DenseSubstrings.substrings[j].singleDeletionReplication) && 
					        intervals[i].isLeftMaximal && 
						    Math.abs(intervals[i].lastPosition,DenseSubstrings.substrings[j].endA)<=distanceThreshold &&
						    intervals[i].nMaximalAlignmentsLeft<MIN_INTERVALS_IN_MAXIMAL &&
							intervals[i].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
						  )		
				        ) {
					if (length<minLength) {
						minLength=length;
						denseSubstring=DenseSubstrings.substrings[j];
					}
				}
			}
			j++;
		}
		return foundDense;
	}
	
	
	/**
	 * Essentially the same as $discardPeriodicSubstrings_straddling()$ in 
	 * $DenseSubstrings$, but for alignment intervals.
	 *
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ and $intervals$ to be 
	 * sorted by first position.
	 *
	 * @param fromInterval the alignment interval from which the search should start;
	 * @return the value of $fromInterval$ to be used for the next range, or -1 if 
	 * no such value can be computed.
	 */
	private static final int cleanAlignmentIntervals_straddling(int rangeFirst, int rangeLast, int rangeFrom, int rangeTo, int fromInterval) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int INTERNAL_THRESHOLD = (3*IDENTITY_THRESHOLD)>>1;
		boolean found;
		int j, k, h;
		int firstJForNext, firstPosition, lastPosition, intervalStartA, intervalEndA;
	
		j=fromInterval; firstJForNext=-1;
		while (j<=lastInterval) {			
			intervalStartA=intervals[j].firstPosition;
			intervalEndA=intervals[j].lastPosition;
			if (intervalStartA>=rangeLast) break;
			if (intervalEndA<=rangeFirst || intervals[j].discarded) {
				j++;
				continue;
			}
			if (firstJForNext==-1 && intervalEndA>rangeLast) firstJForNext=j;
			if ( !Intervals.areApproximatelyIdentical_lowQuality(intervalStartA,intervalEndA,rangeFirst,rangeLast,ReadA.id) &&
				 ( !Intervals.isApproximatelyContained_lowQuality(intervalStartA,intervalEndA,rangeFirst,rangeLast,ReadA.id) ||
				   (intervalStartA>rangeFirst+IDENTITY_THRESHOLD && !Reads.isLeftMaximal(intervalStartA,ReadA.id,true)) ||
				   (intervalEndA<rangeLast-IDENTITY_THRESHOLD && !Reads.isRightMaximal(intervalEndA,ReadA.id,true))
				 )
			   ) {
				j++;
				continue;
		    }
	   		for (k=rangeFrom; k<=rangeTo; k++) {
				firstPosition=PeriodicSubstrings.intervals[k].firstPosition;
				if (firstPosition>=intervalEndA) break;
				lastPosition=PeriodicSubstrings.intervals[k].lastPosition;
				if (lastPosition<=intervalStartA) continue;
				if (PeriodicSubstrings.intervals[k].hasLongPeriod) continue;
				found=false;
				if (lastPosition>=intervalStartA+INTERNAL_THRESHOLD && lastPosition<=intervalEndA-INTERNAL_THRESHOLD) {
					for (h=k+1; h<=rangeTo; h++) {
						if (PeriodicSubstrings.intervals[h].firstPosition>intervalEndA-INTERNAL_THRESHOLD) break;
						if (PeriodicSubstrings.intervals[h].firstPosition<lastPosition-IDENTITY_THRESHOLD) continue;
						if ( PeriodicSubstrings.intervals[h].firstPosition<=lastPosition+IDENTITY_THRESHOLD || 
							 Reads.isRandomInsertion(ReadA.id,lastPosition+1,PeriodicSubstrings.intervals[h].firstPosition-1,true)
						   ) {
						    found=true;
						    break;
					   }
					}
				}
				if (!found) {
					if (firstPosition>=intervalStartA+INTERNAL_THRESHOLD && firstPosition<=intervalEndA-INTERNAL_THRESHOLD) {
						for (h=k-1; h>=rangeFrom; h--) {
							if ( Math.abs(PeriodicSubstrings.intervals[h].lastPosition,firstPosition)<=IDENTITY_THRESHOLD ||
								 ( PeriodicSubstrings.intervals[h].lastPosition<firstPosition-IDENTITY_THRESHOLD && 
								   PeriodicSubstrings.intervals[h].lastPosition>=intervalStartA+INTERNAL_THRESHOLD &&
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
					intervals[j].discarded=true;
					intervals[j].flag2=false;
					break;
				}
			}
			j++;
		}
		return firstJForNext;
	}
	
	
	/**
	 * @param periodicInterval assumed to have long period and to contain 
	 * $alignmentInterval$;
	 * @return TRUE iff $alignmentInterval$ should be discarded using $periodicInterval$.
	 */
	private static final boolean cleanAlignmentIntervals_discardWithLongPeriod(AlignmentInterval alignmentInterval, PeriodicSubstringInterval periodicInterval, boolean left, boolean right, int minIntervalsLong, int minIntervalsShort) {
		final double LONG_PERIOD_TOLERANCE = 0.8;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		
		return ( periodicInterval.period>0 && 
			     ( ( alignmentInterval.length()>=periodicInterval.period*LONG_PERIOD_TOLERANCE &&
			         (left?alignmentInterval.nMaximalAlignmentsLeft<minIntervalsLong:true) && 
					 (right?alignmentInterval.nMaximalAlignmentsRight<minIntervalsLong:true)
			       ) ||
			       ( alignmentInterval.length()<=periodicInterval.period*LONG_PERIOD_TOLERANCE &&
			         (left?alignmentInterval.nMaximalAlignmentsLeft<minIntervalsShort:true) && 
					 (right?alignmentInterval.nMaximalAlignmentsRight<minIntervalsShort:true)
			       )
			     ) 
			   ) ||
		       ( periodicInterval.period<=0 && 
			     ( ( !periodicInterval.equalsOnePeriod &&
					 alignmentInterval.length()>periodicInterval.longestAlignment+IDENTITY_THRESHOLD &&
				     (left?alignmentInterval.nMaximalAlignmentsLeft<minIntervalsLong:true) && 
					 (right?alignmentInterval.nMaximalAlignmentsRight<minIntervalsLong:true)
			       ) ||
			       ( alignmentInterval.length()<=periodicInterval.longestAlignment &&
				     (left?alignmentInterval.nMaximalAlignmentsLeft<minIntervalsShort:true) && 
					 (right?alignmentInterval.nMaximalAlignmentsRight<minIntervalsShort:true)
			       )
			     )
			   );
	}


	/**
	 * If enough alignment intervals end close to the last position of a long-period 
	 * interval, and are contained in the long-period interval, they are all marked as 
	 * discarded. A similar criterion holds for the first position of the long-period
	 * interval.
	 *
	 * Remark: the procedure does not use alignment intervals for which $discarded=true$.
	 *
	 * Remark: the procedure assumes $intervals$ and $PeriodicsSubstrings.
	 * longPeriodIntervals$ to be sorted by first position. The procedure sorts such 
	 * arrays in a different order, and restores the original order at the end. However, 
	 * the order of intervals that start at the same position might not be preserved.
	 *
	 * @param minMaximalAlignments alignment intervals with at least this number of 
	 * maximal alignments inside a periodic interval are not discarded.
	 */
	private static final void cleanAlignmentIntervals_longPeriodEnds(int minMaximalAlignments) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_PATH_LENGTH = 4;  // Arbitrary
		int i, j, k;
		int firstJForNextI, nIntervals, firstPosition, lastPosition;
		if (PeriodicSubstrings.lastLongPeriodInterval==-1) return;
		
		// First position
		for (i=0; i<=lastInterval; i++) intervals[i].flag1=false;
		i=0; j=0; firstJForNextI=-1; nIntervals=0;
		firstPosition=PeriodicSubstrings.longPeriodIntervals[0].firstPosition;
		lastPosition=PeriodicSubstrings.longPeriodIntervals[0].lastPosition;
		while (i<=PeriodicSubstrings.lastLongPeriodInterval) {
			if (j>lastInterval || intervals[j].firstPosition>firstPosition+IDENTITY_THRESHOLD) {
				if (nIntervals>=MIN_PATH_LENGTH) {
					for (k=j; k>=0; k--) {
						if (intervals[k].firstPosition<firstPosition-IDENTITY_THRESHOLD) break;
						if (intervals[k].flag1) {
							intervals[k].discarded=true;
							intervals[k].flag2=false;
						}
					}
				}
				i++; nIntervals=0;
				firstPosition=PeriodicSubstrings.longPeriodIntervals[i].firstPosition;
				lastPosition=PeriodicSubstrings.longPeriodIntervals[i].lastPosition;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (intervals[j].discarded || intervals[j].firstPosition<firstPosition-IDENTITY_THRESHOLD || intervals[j].nMaximalAlignmentsRight>=minMaximalAlignments) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<PeriodicSubstrings.lastLongPeriodInterval && intervals[j].firstPosition>=PeriodicSubstrings.longPeriodIntervals[i+1].firstPosition-IDENTITY_THRESHOLD) firstJForNextI=j;
			if (Intervals.isApproximatelyContained(intervals[j].firstPosition,intervals[j].lastPosition,firstPosition,lastPosition)) {
				intervals[j].flag1=true;
				nIntervals++;
			}
			j++;
		}		
		
		// Last position
		AlignmentInterval.order=AlignmentInterval.LASTPOSITION;
		if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.LASTPOSITION;
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.LASTPOSITION;
		if (PeriodicSubstrings.lastLongPeriodInterval>0) Arrays.sort(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1);
		for (i=0; i<=lastInterval; i++) intervals[i].flag1=false;
		i=0; j=0; firstJForNextI=-1; nIntervals=0;
		firstPosition=PeriodicSubstrings.longPeriodIntervals[0].firstPosition;
		lastPosition=PeriodicSubstrings.longPeriodIntervals[0].lastPosition;
		while (i<=PeriodicSubstrings.lastLongPeriodInterval) {
			if (j>lastInterval || intervals[j].lastPosition>lastPosition+IDENTITY_THRESHOLD) {
				if (nIntervals>=MIN_PATH_LENGTH) {
					for (k=Math.min(j,lastInterval); k>=0; k--) {
						if (intervals[k].lastPosition<lastPosition-IDENTITY_THRESHOLD) break;
						if (intervals[k].flag1) {
							intervals[k].discarded=true;
							intervals[k].flag2=false;
						}
					}
				}
				i++; nIntervals=0;
				firstPosition=PeriodicSubstrings.longPeriodIntervals[i].firstPosition;
				lastPosition=PeriodicSubstrings.longPeriodIntervals[i].lastPosition;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (intervals[j].discarded || intervals[j].lastPosition<lastPosition-IDENTITY_THRESHOLD || intervals[j].nMaximalAlignmentsLeft>=minMaximalAlignments) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<PeriodicSubstrings.lastLongPeriodInterval && intervals[j].lastPosition>=PeriodicSubstrings.longPeriodIntervals[i+1].lastPosition-IDENTITY_THRESHOLD) firstJForNextI=j;
			if (Intervals.isApproximatelyContained(intervals[j].firstPosition,intervals[j].lastPosition,firstPosition,lastPosition)) {
				intervals[j].flag1=true;
				nIntervals++;
			}
			j++;
		}
		
		// Restoring the original order in $intervals$ and $PeriodicSubstrings.
		// longPeriodIntervals$.
		AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
		if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
		if (PeriodicSubstrings.lastLongPeriodInterval>0) Arrays.sort(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1);
	}
	
	
	/**
	 * Essentially identical to $cleanAlignmentIntervals_longPeriodEnds()$.
	 */
	private static final void cleanAlignmentIntervals_shortPeriodEnds(int minMaximalAlignments) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_PATH_LENGTH = 4;  // Arbitrary
		int i, j, k;
		int firstJForNextI, nIntervals, firstPosition, lastPosition;
		if (PeriodicSubstrings.lastInterval==-1) return;
		
		// First position
		for (i=0; i<=lastInterval; i++) intervals[i].flag1=false;
		i=0; 
		while (PeriodicSubstrings.intervals[i].hasLongPeriod) i++;
		if (i>PeriodicSubstrings.lastInterval) return;
		j=0; firstJForNextI=-1; nIntervals=0;
		firstPosition=PeriodicSubstrings.intervals[i].firstPosition;
		lastPosition=PeriodicSubstrings.intervals[i].lastPosition;
		while (i<=PeriodicSubstrings.lastInterval) {
			if (PeriodicSubstrings.intervals[i].hasLongPeriod || j>lastInterval || intervals[j].firstPosition>firstPosition+IDENTITY_THRESHOLD) {
				if (nIntervals>=MIN_PATH_LENGTH) {
					for (k=j; k>=0; k--) {
						if (intervals[k].firstPosition<firstPosition-IDENTITY_THRESHOLD) break;
						if (intervals[k].flag1) {
							intervals[k].discarded=true;
							intervals[k].flag2=false;
						}
					}
				}
				i++; nIntervals=0;
				while (PeriodicSubstrings.intervals[i].hasLongPeriod) i++;
				if (i>PeriodicSubstrings.lastInterval) break;
				firstPosition=PeriodicSubstrings.intervals[i].firstPosition;
				lastPosition=PeriodicSubstrings.intervals[i].lastPosition;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (intervals[j].discarded || intervals[j].firstPosition<firstPosition-IDENTITY_THRESHOLD || intervals[j].nMaximalAlignmentsRight>=minMaximalAlignments) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<PeriodicSubstrings.lastInterval && intervals[j].firstPosition>=PeriodicSubstrings.intervals[i+1].firstPosition-IDENTITY_THRESHOLD) firstJForNextI=j;
			if ( Intervals.isApproximatelyContained(intervals[j].firstPosition,intervals[j].lastPosition,firstPosition,lastPosition) &&
				 PeriodicSubstrings.position2split(intervals[j].lastPosition,PeriodicSubstrings.rightSplits,PeriodicSubstrings.lastRightSplit,IDENTITY_THRESHOLD)==-1
			   ) {
				intervals[j].flag1=true;
				nIntervals++;
			}
			j++;
		}		
		
		// Last position
		AlignmentInterval.order=AlignmentInterval.LASTPOSITION;
		if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.LASTPOSITION;
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.LASTPOSITION;
		if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		for (i=0; i<=lastInterval; i++) intervals[i].flag1=false;
		i=0; 
		while (PeriodicSubstrings.intervals[i].hasLongPeriod) i++;
		if (i>PeriodicSubstrings.lastInterval) return;
		j=0; firstJForNextI=-1; nIntervals=0;
		firstPosition=PeriodicSubstrings.intervals[i].firstPosition;
		lastPosition=PeriodicSubstrings.intervals[i].lastPosition;
		while (i<=PeriodicSubstrings.lastInterval) {
			if (PeriodicSubstrings.intervals[i].hasLongPeriod || j>lastInterval || intervals[j].lastPosition>lastPosition+IDENTITY_THRESHOLD) {
				if (nIntervals>=MIN_PATH_LENGTH) {
					for (k=Math.min(j,lastInterval); k>=0; k--) {
						if (intervals[k].lastPosition<lastPosition-IDENTITY_THRESHOLD) break;
						if (intervals[k].flag1) {
							intervals[k].discarded=true;
							intervals[k].flag2=false;
						}
					}
				}
				i++; nIntervals=0;
				while (PeriodicSubstrings.intervals[i].hasLongPeriod) i++;
				if (i>PeriodicSubstrings.lastInterval) break;
				firstPosition=PeriodicSubstrings.intervals[i].firstPosition;
				lastPosition=PeriodicSubstrings.intervals[i].lastPosition;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (intervals[j].discarded || intervals[j].lastPosition<lastPosition-IDENTITY_THRESHOLD || intervals[j].nMaximalAlignmentsLeft>=minMaximalAlignments) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<PeriodicSubstrings.lastInterval && intervals[j].lastPosition>=PeriodicSubstrings.intervals[i+1].lastPosition-IDENTITY_THRESHOLD) firstJForNextI=j;
			if ( Intervals.isApproximatelyContained(intervals[j].firstPosition,intervals[j].lastPosition,firstPosition,lastPosition) &&
				 PeriodicSubstrings.position2split(intervals[j].firstPosition,PeriodicSubstrings.leftSplits,PeriodicSubstrings.lastLeftSplit,IDENTITY_THRESHOLD)==-1
			   ) {
				intervals[j].flag1=true;
				nIntervals++;
			}
			j++;
		}
		
		// Restoring the original order in $intervals$ and $PeriodicSubstrings.intervals$.
		AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
		if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
		if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
	}

	
	/**
	 * Procedures $fixBoundariesOf*()$ in $DenseSubstrings$ and $PeriodicSubstrings$ might
	 * extend or contract periodic intervals and substring-type intervals. After this, an
	 * alignment interval might become approximately identical to such an interval: the 
	 * procedure discards all these alignment intervals, and reassigns their implied 
	 * alignments. The procedure is a simplified version of $cleanAlignmentIntervals()$.
	 *
	 * Remark: contraction should have been already handled for alignment intervals inside
	 * short-period or substring-type, since they must be aligned to peaks; contraction is
	 * not handled for long-period intervals. Extension is not handled for any type.
	 *
	 * Remark: the procedure assumes all intervals to be sorted by first position.
	 */
	public static final void cleanAlignmentIntervals_afterFixBoundaries() {
		final int DISTANCE_THRESHOLD = IO.quantum;
		boolean foundDense, foundPeriodic;
		int i, j, k;
		int firstJForNext, length, minLength, previousLastInterval;
		AlignmentInterval tmp;
		DenseSubstring denseSubstring;
		PeriodicSubstringInterval periodicSubstringInterval;
		
		for (i=0; i<=lastInterval; i++) {
			intervals[i].discarded=false;
			intervals[i].denseSubstring=null;
			intervals[i].inPeriodicSubstringInterval=null;
			intervals[i].flag2=true;
		}

		// Comparing intervals to dense substrings
		foundDense=cleanAlignmentIntervals_dense(DISTANCE_THRESHOLD);

		// Comparing intervals to periodic substring intervals
		foundPeriodic=false;
		if (PeriodicSubstrings.lastInterval>=0) {
			i=0; j=0; firstJForNext=-1; 
			periodicSubstringInterval=null; minLength=Math.POSITIVE_INFINITY;
			while (i<=lastInterval) {
				if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=intervals[i].lastPosition) {
					if (periodicSubstringInterval!=null) {
						intervals[i].inPeriodicSubstringInterval=periodicSubstringInterval;
						foundPeriodic=true;
						minLength=Math.POSITIVE_INFINITY; periodicSubstringInterval=null;
					}
					i++;
					if (firstJForNext!=-1) j=firstJForNext;
					firstJForNext=-1;
					continue;
				}
				if (PeriodicSubstrings.intervals[j].lastPosition<=intervals[i].firstPosition) {
					j++;
					continue;
				}				
				if (i<lastInterval && firstJForNext==-1 && PeriodicSubstrings.intervals[j].lastPosition>=intervals[i+1].firstPosition) firstJForNext=j;
				length=PeriodicSubstrings.intervals[j].lastPosition-PeriodicSubstrings.intervals[j].firstPosition;
				if (Intervals.areApproximatelyIdentical_lowQuality(PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,intervals[i].firstPosition,intervals[i].lastPosition,ReadA.id)) {
					if (length<minLength) {
						minLength=length;
						periodicSubstringInterval=PeriodicSubstrings.intervals[j];
					}
				}
				else if (Intervals.isApproximatelyContained_lowQuality(intervals[i].firstPosition,intervals[i].lastPosition,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.id)) {
					if ( ( Math.abs(intervals[i].firstPosition,PeriodicSubstrings.intervals[j].firstPosition)<=DISTANCE_THRESHOLD && 
						   (!intervals[i].isRightMaximal || intervals[i].nMaximalAlignmentsRight<MIN_INTERVALS_IN_MAXIMAL)
						 ) ||
						 ( Math.abs(intervals[i].lastPosition,PeriodicSubstrings.intervals[j].lastPosition)<=DISTANCE_THRESHOLD && 
						   (!intervals[i].isLeftMaximal || intervals[i].nMaximalAlignmentsLeft<MIN_INTERVALS_IN_MAXIMAL)
						 )
					   ) {
						if (length<minLength) {
							minLength=length;
							periodicSubstringInterval=PeriodicSubstrings.intervals[j];
						}
					}
				}
				j++;
			}
			// Resetting pointers from long-period clones to the real long-period
			// intervals.
			for (i=0; i<=lastInterval; i++) {
				if (intervals[i].inPeriodicSubstringInterval==null || !intervals[i].inPeriodicSubstringInterval.hasLongPeriod) continue;
				j=Arrays.binarySearch(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1,intervals[i].inPeriodicSubstringInterval);
				if (j<0) {
					System.err.println("cleanAlignmentIntervals_afterFixBoundaries> ERROR: long-period interval not found.");
					System.exit(1);
				}
				foundPeriodic=false;
				for (k=j; k>=0; k--) {
					if (PeriodicSubstrings.longPeriodIntervals[k].equals(intervals[i].inPeriodicSubstringInterval)) {
						intervals[i].inPeriodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[k];
						foundPeriodic=true;
						break;
					}
				}
				if (foundPeriodic) continue;
				for (k=j+1; k<=PeriodicSubstrings.lastLongPeriodInterval; k++) {
					if (PeriodicSubstrings.longPeriodIntervals[k].equals(intervals[i].inPeriodicSubstringInterval)) {
						intervals[i].inPeriodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[k];
						break;
					}
				}
			}
		}
		
		// Compacting intervals
		previousLastInterval=lastInterval;
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].denseSubstring!=null) {
				intervals[i].discarded=true;
				if (intervals[i].denseSubstring.substringReplication) intervals[i].flag2=false;
				continue;
			}
			if (intervals[i].inPeriodicSubstringInterval!=null) {
				intervals[i].discarded=true;
				intervals[i].flag2=false;
				continue;
			}
			j++;
			if (j!=i) {
				tmp=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=tmp;
			}
		}
		lastInterval=j;
		if (lastInterval==previousLastInterval) return;
		
		// Updating pointers from alignments
		discardIntervals_reassignAlignments(ReadA.lastSortedAlignment+1,3,true,true);
	}
	
	
	/**
	 * To every interval that is contained in a dense substring of substring type, the 
	 * procedure assigns the shortest such containing substring.
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ to be sorted by 
	 * $startA$.
	 */
	public static final void markIntervalsInDenseSubstrings() {
		final int DISTANCE_THRESHOLD = IO.quantum;
		int i, j;
		int firstJForNextI, lastPoint, length, minLength;
		DenseSubstring tmpSubstring, minSubstring;
		
		// Ensuring the necessary order in $intervals$
		if (lastInterval==-1) return;
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		
		// Marking
		for (i=0; i<=lastInterval; i++) intervals[i].inDenseSubstringOfSubstringType=null;
		i=0; j=0; firstJForNextI=-1; lastPoint=-1;
		minSubstring=null; minLength=Math.POSITIVE_INFINITY;
		while (i<=lastInterval) {
			if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>intervals[i].lastPosition) {
				if (minSubstring!=null) {
					intervals[i].inDenseSubstringOfSubstringType=minSubstring;
					minSubstring=null;
					minLength=Math.POSITIVE_INFINITY;
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			tmpSubstring=DenseSubstrings.substrings[j];
			if (tmpSubstring.endA<intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (i<lastInterval && firstJForNextI==-1 && tmpSubstring.endA>=intervals[i+1].firstPosition) firstJForNextI=j;
			if ( tmpSubstring.substringReplication && 
			     Intervals.isApproximatelyContained(intervals[i].firstPosition,intervals[i].lastPosition,tmpSubstring.startA,tmpSubstring.endA) &&
				 !( (intervals[i].firstPosition>tmpSubstring.startA+DISTANCE_THRESHOLD && !Reads.isLeftMaximal(intervals[i].firstPosition,ReadA.id,true)) &&
					(intervals[i].lastPosition<tmpSubstring.endA-DISTANCE_THRESHOLD && !Reads.isRightMaximal(intervals[i].lastPosition,ReadA.id,true))
				 )
			   ) {
				length=tmpSubstring.endA-tmpSubstring.startA;
				if (length<minLength) {
					minLength=length;
					minSubstring=tmpSubstring;
				}
			}
			j++;
		}
	}
	
	
	/**
	 * To every interval that is contained in a periodic interval, the procedure assigns 
	 * the shortest such interval to field $inPeriodicSubstringInterval$ (this variable is
	 * set to null otherwise). The procedure sets $inPeriodicRange=true$ for every 
	 * interval that is contained in a maximal range of straddling periodic intervals.
	 *
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ to be sorted by 
	 * $firstPosition$.
	 *
	 * @param first,last only this range of $intervals$ is considered; the procedure 
	 * assumes this range to be already sorted by $firstPosition$;
	 * @param justShort only short-period intervals are considered.
	 */
	public static final void markIntervalsInPeriodic(int first, int last, boolean justShort) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		int i, j;
		int firstJForNextI, lastPoint, length, minLength, rangeFirst, rangeLast;
		PeriodicSubstringInterval tmpInterval, minInterval;
		
		// Setting $inPeriodicSubstringInterval$.
		for (i=first; i<=last; i++) intervals[i].inPeriodicSubstringInterval=null;
		i=first; j=0; firstJForNextI=-1; lastPoint=-1;
		minInterval=null; minLength=Math.POSITIVE_INFINITY;
		while (i<=last) {
			if (j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>intervals[i].lastPosition) {
				if (minInterval!=null) {
					intervals[i].inPeriodicSubstringInterval=minInterval;
					minInterval=null;
					minLength=Math.POSITIVE_INFINITY;
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			tmpInterval=PeriodicSubstrings.intervals[j];
			if (tmpInterval.lastPosition<intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (i<last && firstJForNextI==-1 && tmpInterval.lastPosition>=intervals[i+1].firstPosition) firstJForNextI=j;
			if ( (justShort?!tmpInterval.hasLongPeriod:true) && 
				 ( Intervals.isApproximatelyContained(intervals[i].firstPosition,intervals[i].lastPosition,tmpInterval.firstPosition,tmpInterval.lastPosition) ||
				   Intervals.areApproximatelyIdentical(intervals[i].firstPosition,intervals[i].lastPosition,tmpInterval.firstPosition,tmpInterval.lastPosition)
				 )
			   ) {
				length=tmpInterval.lastPosition-tmpInterval.firstPosition;
				if (length<minLength) {
					minLength=length;
					minInterval=tmpInterval;
				}
			}
			j++;
		}
		
		// Setting $inPeriodicRange$.
		for (i=0; i<=lastInterval; i++) intervals[i].inPeriodicRange=false;
		i=0;
		if (justShort) {
			while (i<=PeriodicSubstrings.lastInterval && PeriodicSubstrings.intervals[i].hasLongPeriod) i++;
		}
		if (i>PeriodicSubstrings.lastInterval) return;
		rangeFirst=PeriodicSubstrings.intervals[i].firstPosition;
		rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
		j=first; firstJForNextI=-1;
		for (++i; i<=PeriodicSubstrings.lastInterval; i++) {
			if (justShort && PeriodicSubstrings.intervals[i].hasLongPeriod) continue;
			if (PeriodicSubstrings.intervals[i].firstPosition<=rangeLast+IO.quantum) {
				if (PeriodicSubstrings.intervals[i].lastPosition>rangeLast) rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
				continue;
			}
			if (firstJForNextI!=-1) j=firstJForNextI;
			firstJForNextI=-1;
			while (j<=last) {
				if (intervals[j].firstPosition>=rangeLast) break;
				if (intervals[j].lastPosition<rangeFirst) {
					j++;
					continue;
				}
				if (firstJForNextI==-1 && intervals[j].lastPosition>rangeLast) firstJForNextI=j;
				if ( Intervals.isApproximatelyContained(intervals[j].firstPosition,intervals[j].lastPosition,rangeFirst,rangeLast) ||
					 Intervals.areApproximatelyIdentical(intervals[j].firstPosition,intervals[j].lastPosition,rangeFirst,rangeLast)
				   ) intervals[j].inPeriodicRange=true;
				j++;
			}
			// Next range
			rangeFirst=PeriodicSubstrings.intervals[i].firstPosition;
			rangeLast=PeriodicSubstrings.intervals[i].lastPosition;
		}
		// Last range
		if (firstJForNextI!=-1) j=firstJForNextI;
		while (j<=last) {
			if (intervals[j].firstPosition>=rangeLast) break;
			if (intervals[j].lastPosition<rangeFirst) {
				j++;
				continue;
			}
			if ( Intervals.isApproximatelyContained(intervals[j].firstPosition,intervals[j].lastPosition,rangeFirst,rangeLast) ||
				 Intervals.areApproximatelyIdentical(intervals[j].firstPosition,intervals[j].lastPosition,rangeFirst,rangeLast)
			   ) intervals[j].inPeriodicRange=true;
			j++;
		}
	}
	
	
	/**
	 * Detects whether an alignment interval is adjacent to a high-quality substring to
	 * the left and to the right of its readA.
	 */
	public static final void setMaximality() {
		for (int i=0; i<=lastInterval; i++) {
			intervals[i].isLeftMaximal=Reads.isLeftMaximal(intervals[i].firstPosition,ReadA.id,true);
			intervals[i].isRightMaximal=Reads.isRightMaximal(intervals[i].lastPosition,ReadA.id,true);
		}
	}

	
	/**
	 * Assume that some alignment intervals have been marked as discarded and have been
	 * pushed out of $intervals[0..lastInterval]$. The procedure tries to reassign every
	 * alignment in $ReadA.sortedAlignments[0..firstImplied)$ whose $mergedToInterval$ 
	 * pointer is a discarded interval, to a shortest surviving interval that is 
	 * compatible with the alignment in terms of maximality (possibly to no interval).
	 *
	 * Remark: the procedure assumes $intervals[0..lastInterval]$ to be sorted by 
	 * $firstPosition$.
	 *
	 * @param assignmentMode if no surviving interval can be assigned to an alignment, the
	 * procedure: (0) does not do anything; (1) sets the $inDenseSubstring$ field of the 
	 * alignment to the $inDenseSubstringOfSubstringType$ field of its discarded alignment
	 * interval, if any; (2) sets the $periodicSubstringInterval$ field of the alignment 
	 * to the $inPeriodicSubstringInterval$ field of its discarded alignment interval, if 
	 * any; (3) sets one of the $impliedByDenseSubstring$, $inDenseSubstring$, and 
	 * $periodicSubstringInterval$ fields of the alignment to the $denseSubstring$ or 
	 * $inPeriodicSubstringInterval$ fields of its discarded alignment interval, if any.
	 * @param incrementCounts increments fields $nMergedIntervals,nMaximalAlignments*$ of
	 * an alignment interval;
	 * @param onlyMarked tries to reassign to an alignment interval, only alignments whose
	 * discarded alignment interval has $flag2=true$;
	 * @return TRUE iff the procedure has assigned the $impliedByDenseSubstring$ or the
	 * $periodicSubstringInterval$ pointer of an element of $ReadA.sortedAlignments$.
	 */
	public static final boolean discardIntervals_reassignAlignments(int firstImplied, int assignmentMode, boolean incrementCounts, boolean onlyMarked) {
		final int THRESHOLD = IO.quantum;
		boolean isRightMaximal, isLeftMaximal, denseOrPeriodic, out, previousSortByID;
		int i, j, p;
		int startA, endA, start1, end1, start2, end2, start3, end3, length, minLength;
		DenseSubstring substring;
		AlignmentInterval oldInterval, minInterval;
		AlignmentInterval tmpInterval = new AlignmentInterval();
		PeriodicSubstringInterval periodicSubstringInterval;
		
		out=false;
		for (i=0; i<firstImplied; i++) {
			oldInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (oldInterval==null || !oldInterval.discarded) continue;
			isRightMaximal=ReadA.sortedAlignments[i].isRightMaximal==1;
			isLeftMaximal=ReadA.sortedAlignments[i].isLeftMaximal==1;
			startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			minInterval=null; minLength=Math.POSITIVE_INFINITY;
			if (!onlyMarked || oldInterval.flag2) {
				tmpInterval.firstPosition=startA;
				tmpInterval.lastPosition=endA;
				previousSortByID=AlignmentInterval.sortByID;
				AlignmentInterval.sortByID=false;
				p=Arrays.binarySearch(intervals,0,lastInterval+1,tmpInterval);
				AlignmentInterval.sortByID=previousSortByID;
				if (p<0) p=-p-1;
				p=Math.min(p,lastInterval);		
				for (j=p; j>=0; j--) {
					if ( intervals[j].canBeAssigned(ReadA.sortedAlignments[i],THRESHOLD) &&
					     ( Intervals.areApproximatelyIdentical_lowQuality(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id) ||
						   Intervals.isApproximatelyContained_lowQuality(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id) ||
						   Intervals.contains_lowQuality(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id)
						 )
					   ) {
						length=intervals[j].length();
						if (length<minLength) {
							minLength=length;
							minInterval=intervals[j];
						}
					}
				}
				for (j=p+1; j<=lastInterval; j++) {
					if (intervals[j].firstPosition>=endA-THRESHOLD) break;
					if ( intervals[j].canBeAssigned(ReadA.sortedAlignments[i],THRESHOLD) &&
						 ( Intervals.areApproximatelyIdentical_lowQuality(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id) ||
	                       Intervals.isApproximatelyContained_lowQuality(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id) ||
						   Intervals.contains_lowQuality(startA,endA,intervals[j].firstPosition,intervals[j].lastPosition,ReadA.id)
						 )
					   ) {
						length=intervals[j].length();
						if (length<minLength) {
							minLength=length;
							minInterval=intervals[j];
						}
					}
				}
			}
			ReadA.sortedAlignments[i].mergedToInterval=null;
			
			// Assignment modes
			if (assignmentMode==0) {
				if (minInterval!=null) {
					ReadA.sortedAlignments[i].mergedToInterval=minInterval;
					if (incrementCounts) {
						minInterval.nMergedIntervals++;
						if (isLeftMaximal) minInterval.nMaximalAlignmentsLeft++;
						if (isRightMaximal) minInterval.nMaximalAlignmentsRight++;
					}
				}
			}
			else if (assignmentMode==1) {
				substring=oldInterval.inDenseSubstringOfSubstringType;
				if ( substring!=null && 
				     (ReadA.sortedAlignments[i].inDenseSubstring==null || substring.length()<ReadA.sortedAlignments[i].inDenseSubstring.length()) &&
					 (minInterval==null || substring.length()<minLength || !Intervals.isApproximatelyContained(minInterval.firstPosition,minInterval.lastPosition,substring.startA,substring.endA))
				   ) {
					ReadA.sortedAlignments[i].inDenseSubstring=substring;
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
				}
				else if ( minInterval!=null && 
				          ( substring==null || 
							(minLength<substring.length() && Intervals.isApproximatelyContained(minInterval.firstPosition,minInterval.lastPosition,substring.startA,substring.endA))
						  ) &&
						  (ReadA.sortedAlignments[i].inDenseSubstring==null || minLength<ReadA.sortedAlignments[i].inDenseSubstring.length())
				        ) {
					ReadA.sortedAlignments[i].mergedToInterval=minInterval;
					if (incrementCounts) {
						minInterval.nMergedIntervals++;
						if (isLeftMaximal) minInterval.nMaximalAlignmentsLeft++;
						if (isRightMaximal) minInterval.nMaximalAlignmentsRight++;
					}
				}
			}
			else if (assignmentMode==2) {
				periodicSubstringInterval=oldInterval.inPeriodicSubstringInterval;
				if ( periodicSubstringInterval!=null && 
				     (ReadA.sortedAlignments[i].periodicSubstringInterval==null || periodicSubstringInterval.length()<ReadA.sortedAlignments[i].periodicSubstringInterval.length()) &&
					 (minInterval==null || periodicSubstringInterval.length()<minLength || !Intervals.isApproximatelyContained(minInterval.firstPosition,minInterval.lastPosition,periodicSubstringInterval.firstPosition,periodicSubstringInterval.lastPosition))
				   ) {
					ReadA.sortedAlignments[i].periodicSubstringInterval=periodicSubstringInterval;
					ReadA.sortedAlignments[i].inDenseSubstring=null;
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					out=true;
				}
				else if ( minInterval!=null && 
				          ( periodicSubstringInterval==null || 
							(minLength<periodicSubstringInterval.length() && Intervals.isApproximatelyContained(minInterval.firstPosition,minInterval.lastPosition,periodicSubstringInterval.firstPosition,periodicSubstringInterval.lastPosition))
						  ) &&
						  (ReadA.sortedAlignments[i].periodicSubstringInterval==null || minLength<ReadA.sortedAlignments[i].periodicSubstringInterval.length())
				        ) {
					ReadA.sortedAlignments[i].mergedToInterval=minInterval;
					if (incrementCounts) {
						minInterval.nMergedIntervals++;
						if (isLeftMaximal) minInterval.nMaximalAlignmentsLeft++;
						if (isRightMaximal) minInterval.nMaximalAlignmentsRight++;
					}
					ReadA.sortedAlignments[i].periodicSubstringInterval=null;
				}
			}
			else if (assignmentMode==3) {	
				substring=oldInterval.denseSubstring;
				periodicSubstringInterval=oldInterval.inPeriodicSubstringInterval;		
				if (substring==null && periodicSubstringInterval==null) {
					if (minInterval!=null) {
						ReadA.sortedAlignments[i].mergedToInterval=minInterval;
						if (incrementCounts) {
							minInterval.nMergedIntervals++;
							if (isLeftMaximal) minInterval.nMaximalAlignmentsLeft++;
							if (isRightMaximal) minInterval.nMaximalAlignmentsRight++;
						}
					}
					continue;
				}
				if ((substring!=null && substring.canBeAssigned(ReadA.sortedAlignments[i],THRESHOLD)) && periodicSubstringInterval!=null) {
					start1=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
					end1=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
					start2=substring.startA;
					end2=substring.endA;
					start3=periodicSubstringInterval.firstPosition;
					end3=periodicSubstringInterval.lastPosition;
					if (Intervals.jaccardSimilarity(start1,end1,start2,end2)>Intervals.jaccardSimilarity(start1,end1,start3,end3)) denseOrPeriodic=true;
					else denseOrPeriodic=false;
				}
				else if ((substring!=null && substring.canBeAssigned(ReadA.sortedAlignments[i],THRESHOLD)) && periodicSubstringInterval==null) denseOrPeriodic=true;
				else denseOrPeriodic=false;
				if (denseOrPeriodic) {
					if (substring.substringReplication) {
						if ( (minInterval==null || substring.length()<minLength || !Intervals.isApproximatelyContained(minInterval.firstPosition,minInterval.lastPosition,substring.startA,substring.endA)) &&
							 (ReadA.sortedAlignments[i].inDenseSubstring==null || substring.length()<ReadA.sortedAlignments[i].inDenseSubstring.length())
						   ) {
							ReadA.sortedAlignments[i].inDenseSubstring=substring;
							ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
						}
						else if ( minInterval!=null && 
							      (minLength<substring.length() && Intervals.isApproximatelyContained(minInterval.firstPosition,minInterval.lastPosition,substring.startA,substring.endA)) && 
								  (ReadA.sortedAlignments[i].inDenseSubstring==null || minLength<ReadA.sortedAlignments[i].inDenseSubstring.length())
						        ) {
							ReadA.sortedAlignments[i].mergedToInterval=minInterval;
							if (incrementCounts) {
								minInterval.nMergedIntervals++;
								if (isLeftMaximal) minInterval.nMaximalAlignmentsLeft++;
								if (isRightMaximal) minInterval.nMaximalAlignmentsRight++;
							}
						}
					}
					else {
						if ( (minInterval==null || substring.length()<minLength) &&
						     (ReadA.sortedAlignments[i].impliedByDenseSubstring==null || substring.length()<ReadA.sortedAlignments[i].impliedByDenseSubstring.length())
						   ) {
							ReadA.sortedAlignments[i].impliedByDenseSubstring=substring;
							ReadA.sortedAlignments[i].inDenseSubstring=null;
							out=true;
						}
						else if ( minInterval!=null && minLength<substring.length() && 
								  (ReadA.sortedAlignments[i].impliedByDenseSubstring==null || minLength<ReadA.sortedAlignments[i].impliedByDenseSubstring.length())
						        ) {
							ReadA.sortedAlignments[i].mergedToInterval=minInterval;
							if (incrementCounts) {
								minInterval.nMergedIntervals++;
								if (isLeftMaximal) minInterval.nMaximalAlignmentsLeft++;
								if (isRightMaximal) minInterval.nMaximalAlignmentsRight++;
							}
							ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
						}
					}
				}
				else {
					if ( periodicSubstringInterval!=null &&
						 (ReadA.sortedAlignments[i].periodicSubstringInterval==null || periodicSubstringInterval.length()<ReadA.sortedAlignments[i].periodicSubstringInterval.length()) &&
						 (minInterval==null || periodicSubstringInterval.length()<minLength || !Intervals.isApproximatelyContained(minInterval.firstPosition,minInterval.lastPosition,periodicSubstringInterval.firstPosition,periodicSubstringInterval.lastPosition))
					   ) {
						ReadA.sortedAlignments[i].periodicSubstringInterval=periodicSubstringInterval;
						ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
						ReadA.sortedAlignments[i].inDenseSubstring=null;
						out=true;
					}
					else if ( minInterval!=null && 
				          	  ( periodicSubstringInterval==null || 
								(minLength<periodicSubstringInterval.length() && Intervals.isApproximatelyContained(minInterval.firstPosition,minInterval.lastPosition,periodicSubstringInterval.firstPosition,periodicSubstringInterval.lastPosition))
							  ) &&
						      (ReadA.sortedAlignments[i].periodicSubstringInterval==null || minLength<ReadA.sortedAlignments[i].periodicSubstringInterval.length())
				            ) {
						ReadA.sortedAlignments[i].mergedToInterval=minInterval;
						if (incrementCounts) {
							minInterval.nMergedIntervals++;
							if (isLeftMaximal) minInterval.nMaximalAlignmentsLeft++;
							if (isRightMaximal) minInterval.nMaximalAlignmentsRight++;
						}
						ReadA.sortedAlignments[i].periodicSubstringInterval=null;
					}
				}
			}
		}
		return out;
	}
	
	
	/**
	 * Tries to ensure that no two alignment intervals are approximately identical. This 
	 * is done by greedily merging intervals from left to right, so the desired outcome is
	 * not guaranteed. The procedure merges also intervals that are contained in others,
	 * iff their $nMergedIntervals$ field is not too large WRT the container.
	 *
	 * Remark: the procedure updates the following properties of the result of a merge: 
	 * $firstPosition,lastPosition,nMergedIntervals,nMaximalAlignments*,isLongRepeat$.
	 */
	private static final void mergeSimilarIntervals(int firstImplied) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final double MAX_MERGED_INTERVALS_RATIO = 2.0;  // Arbitrary
		boolean isLeftMaximal, isRightMaximal, highQualityLeft, highQualityRight;
		int i, j;
		int type, sumStart, sumEnd, firstPosition, lastPosition;
		int nIntervals, nIntervalsLeft, nIntervalsRight;
		int nMergedIntervals, nMaximalAlignmentsLeft, nMaximalAlignmentsRight, min, max, minInterval, minLength;
		AlignmentInterval tmpInterval;
		
		// Ensuring the necessary order in $intervals$
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		
		// Merging identical intervals and updating boundaries
		for (i=0; i<=lastInterval; i++) intervals[i].discarded=false;
		for (i=0; i<lastInterval; i++) {
			if (intervals[i].discarded) continue;
			firstPosition=intervals[i].firstPosition;
			lastPosition=intervals[i].lastPosition;
			nMergedIntervals=intervals[i].nMergedIntervals;
			sumStart=firstPosition*nMergedIntervals; 
			sumEnd=lastPosition*nMergedIntervals;
			min=firstPosition; max=lastPosition;
			isLeftMaximal=intervals[i].isLeftMaximal;
			isRightMaximal=intervals[i].isRightMaximal;
			nIntervals=nMergedIntervals; nIntervalsLeft=nMergedIntervals; nIntervalsRight=nMergedIntervals;
			nMaximalAlignmentsLeft=intervals[i].nMaximalAlignmentsLeft;
			nMaximalAlignmentsRight=intervals[i].nMaximalAlignmentsRight;
			for (j=i+1; j<=lastInterval; j++) {
				if (intervals[j].firstPosition>=lastPosition) break;
				if (intervals[j].discarded) continue;
				if (Intervals.areApproximatelyIdentical(intervals[j].firstPosition,intervals[j].lastPosition,sumStart/nIntervalsLeft,sumEnd/nIntervalsRight)) {
					nMergedIntervals=intervals[j].nMergedIntervals;
					sumStart+=intervals[j].firstPosition*nMergedIntervals;
					sumEnd+=intervals[j].lastPosition*nMergedIntervals;
					min=Math.min(min,intervals[j].firstPosition);
					max=Math.max(max,intervals[j].lastPosition);
					nIntervals+=nMergedIntervals; nIntervalsLeft+=nMergedIntervals; nIntervalsRight+=nMergedIntervals;
					nMaximalAlignmentsLeft+=intervals[j].nMaximalAlignmentsLeft;
					nMaximalAlignmentsRight+=intervals[j].nMaximalAlignmentsRight;
					intervals[i].isLongRepeat|=intervals[j].isLongRepeat;
					isLeftMaximal|=intervals[j].isLeftMaximal;
					isRightMaximal|=intervals[j].isRightMaximal;
					intervals[j].discarded=true;
				}
			}
			j=sumStart/nIntervalsLeft;
			intervals[i].firstPosition=min<j-DISTANCE_THRESHOLD?min:j;
			j=sumEnd/nIntervalsRight;
			intervals[i].lastPosition=max>j+DISTANCE_THRESHOLD?max:j;
			intervals[i].nMergedIntervals=nIntervals;
			intervals[i].nMaximalAlignmentsLeft=nMaximalAlignmentsLeft;
			intervals[i].nMaximalAlignmentsRight=nMaximalAlignmentsRight;
			intervals[i].isLeftMaximal=isLeftMaximal;
			intervals[i].isRightMaximal=isRightMaximal;
		}
		
		// Merging contained intervals and possibly updating boundaries
		for (i=0; i<lastInterval; i++) {
			if (intervals[i].discarded) continue;			
			isLeftMaximal=intervals[i].isLeftMaximal;
			isRightMaximal=intervals[i].isRightMaximal;
			if (isLeftMaximal && isRightMaximal) continue;
			firstPosition=intervals[i].firstPosition;
			lastPosition=intervals[i].lastPosition;
			highQualityLeft=Reads.isLeftMaximal(firstPosition,ReadA.id,true);
			highQualityRight=Reads.isRightMaximal(lastPosition,ReadA.id,true);
			nMergedIntervals=intervals[i].nMergedIntervals;
			minInterval=-1; minLength=Math.POSITIVE_INFINITY;
			for (j=i+1; j<=lastInterval; j++) {
				if (intervals[j].firstPosition>=lastPosition) break;
				if (intervals[j].discarded || nMergedIntervals>intervals[j].nMergedIntervals*MAX_MERGED_INTERVALS_RATIO) continue;
				if ( Intervals.isApproximatelyContained(firstPosition,lastPosition,intervals[j].firstPosition,intervals[j].lastPosition) &&
					 !( firstPosition>intervals[j].firstPosition+DISTANCE_THRESHOLD && 
						(isLeftMaximal || !highQualityLeft)
					 ) &&
					 !( lastPosition<intervals[j].lastPosition-DISTANCE_THRESHOLD && 
					    (isRightMaximal || !highQualityRight)
					 ) &&
					 intervals[j].length()<minLength
				   ) {
					minInterval=j;
					minLength=intervals[j].length();
				}
			}
			for (j=i-1; j>=0; j--) {
				if (intervals[j].lastPosition<firstPosition) continue;
				if (intervals[j].discarded || nMergedIntervals>intervals[j].nMergedIntervals*MAX_MERGED_INTERVALS_RATIO) continue;
				if ( Intervals.isApproximatelyContained(firstPosition,lastPosition,intervals[j].firstPosition,intervals[j].lastPosition) &&
					 !( firstPosition>intervals[j].firstPosition+DISTANCE_THRESHOLD && 
						(isLeftMaximal || !highQualityLeft)
					 ) &&
					 !( lastPosition<intervals[j].lastPosition-DISTANCE_THRESHOLD && 
					    (isRightMaximal || !highQualityRight)
					 ) &&
					 intervals[j].length()<minLength
				   ) {
					minInterval=j;
					minLength=intervals[j].length();
				}
			}
			if (minInterval==-1) continue;
			intervals[i].discarded=true;
			if (firstPosition<=intervals[minInterval].firstPosition+DISTANCE_THRESHOLD) {
				j=(intervals[minInterval].firstPosition*intervals[minInterval].nMergedIntervals+firstPosition*nMergedIntervals)/(intervals[minInterval].nMergedIntervals+nMergedIntervals);
				min=Math.min(firstPosition,intervals[minInterval].firstPosition);
				intervals[minInterval].firstPosition=min<j-DISTANCE_THRESHOLD?min:j;
				intervals[minInterval].isLeftMaximal|=isLeftMaximal;
			}
			if (lastPosition>=intervals[minInterval].lastPosition-DISTANCE_THRESHOLD) {
				j=(intervals[minInterval].lastPosition*intervals[minInterval].nMergedIntervals+lastPosition*nMergedIntervals)/(intervals[minInterval].nMergedIntervals+nMergedIntervals);
				max=Math.max(lastPosition,intervals[minInterval].lastPosition);
				intervals[minInterval].lastPosition=max>j+DISTANCE_THRESHOLD?max:j;
				intervals[minInterval].isRightMaximal|=isRightMaximal;
			}
			intervals[minInterval].nMergedIntervals+=nMergedIntervals;
		}
		
		// Compacting intervals
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].discarded) continue;
			j++;
			if (j!=i) {
				tmpInterval=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=tmpInterval;
			}
		}
		lastInterval=j;
		
		// Updating pointers from alignments
		AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
		if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);  // Since the first position might have changed
		discardIntervals_reassignAlignments(firstImplied,0,false,false);
	}
	
	
	/**
	 * Checks whether all alignments point to a valid alignment interval.
	 */
	public static final void checkConsistency() {
		final int THRESHOLD = Alignments.minAlignmentLength;
		final int QUALITY_THRESHOLD = IO.quantum*3;
		boolean found;
		int i, j;
		int start, end;
		AlignmentInterval tmpInterval;
		
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			tmpInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (tmpInterval==null) continue;
			found=false;
			for (j=0; j<=lastInterval; j++) {
				if (intervals[j]==tmpInterval) {
					found=true;
					break;
				}
			}
			if (!found) {
				System.err.println("ERROR: the following alignment points to an absent interval:");
				System.err.println("ALIGNMENT: "+ReadA.sortedAlignments[i].toStringBoundaries());
				System.err.println("ABSENT INTERVAL: "+tmpInterval);
				System.exit(1);
			}
			start=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			end=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if ( ( start<tmpInterval.firstPosition-THRESHOLD && 
				   !ReadA.sortedAlignments[i].lowQualityStart &&
				   !Reads.hasLowQuality(ReadA.id,start,tmpInterval.firstPosition-1,true) &&
				   !Reads.hasLowQuality(ReadA.id,start,start+QUALITY_THRESHOLD-1,true)
				 ) || 
			     ( end>tmpInterval.lastPosition+THRESHOLD && 
				   !ReadA.sortedAlignments[i].lowQualityEnd && 
				   !Reads.hasLowQuality(ReadA.id,tmpInterval.lastPosition+1,end,true) &&
				   !Reads.hasLowQuality(ReadA.id,end-QUALITY_THRESHOLD+1,end,true)
				 )
			   ) {
				System.err.println("ERROR: the following alignment points to a shorter interval?!");
				System.err.println("ALIGNMENT: "+ReadA.sortedAlignments[i].toStringBoundaries());
				System.err.println("SHORTER INTERVAL: "+tmpInterval);
				System.exit(1);
			}
		}
	}








	// ----------------------- INTERVALS TO SUBSTRINGS PROCEDURES ------------------------

	/**
	 * @param tmp1,tmp2 temporary space, of length at least equal to the number of 
	 * alignment intervals. This must not be a pointer to $AlignmentIntervals.stack$.
	 */
	public static final void intervals2denseSubstrings(int[] tmp1, int[] tmp2) {
		int from;
		if (lastInterval<MIN_PATH_LENGTH) return;


if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings at the very beginning:");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);	
	IO.printErr("alignment intervals at the very beginning:");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments at the very beginning:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());	
}
		
		
		markLowQualityIntervals();
		
		// Transforming intervals into weak substrings of substring type
		from=DenseSubstrings.lastSubstring+1;
		if (DenseSubstring.order!=DenseSubstring.STARTA) {
			DenseSubstring.order=DenseSubstring.STARTA;
			if (from-1>0) Arrays.sort(DenseSubstrings.substrings,0,from);
		}
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		intervals2substring(tmp1,tmp2);
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after intervals2substring():   from="+from+" DenseSubstrings.lastSubstring="+DenseSubstrings.lastSubstring);
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);		
	IO.printErr("alignments after intervals2substring():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());	
	IO.printErr("alignment intervals after intervals2substring():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
}
	
		if (DenseSubstrings.lastSubstring!=from-1) {
			DenseSubstrings.fixBoundariesOfSubstringType(tmp);
			markAlignmentsInNewDenseSubstrings(from);
			DenseSubstrings.fixBoundariesOfSubstringType(tmp);
			DenseSubstrings.computeDensityWindows(ReadA.lastSortedAlignment);
			from=cleanNewDenseSubstrings(from);
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after cleanNewDenseSubstrings():");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);		
	IO.printErr("alignments after cleanNewDenseSubstrings():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());	
	IO.printErr("alignment intervals after cleanNewDenseSubstrings():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
}
			
			
			cleanNewDenseSubstrings_periodic(from);
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after cleanNewDenseSubstrings_periodic():");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);	
	IO.printErr("alignments after cleanNewDenseSubstrings_periodic():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());	
	IO.printErr("alignment intervals after cleanNewDenseSubstrings_periodic():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
}			
			
			DenseSubstring.order=DenseSubstring.STARTA;
			if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
			markIntervalsInDenseSubstrings();
		}


if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignment intervals after AlignmentIntervals.intervals2denseSubstrings():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments after AlignmentIntervals.intervals2denseSubstrings():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
}

		
		// Discarding prefix/suffix substrings inside substrings of substring type.
		intervals2prefixSuffix_inSubstring(tmp1);


if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after intervals2prefixSuffix_inSubstring():");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);		
	IO.printErr("alignment intervals after AlignmentIntervals.intervals2prefixSuffix_inSubstring():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments after AlignmentIntervals.intervals2prefixSuffix_inSubstring():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
}
		
		// Transforming intervals into weak prefix/suffix substrings
		from=DenseSubstrings.lastSubstring+1;
		intervals2prefixSuffix(tmp1);
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after intervals2prefixSuffix():  from="+from);
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);
	IO.printErr("intervals2denseSubstrings> alignment intervals after intervals2prefixSuffix():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);	
}		
		
		if (DenseSubstrings.lastSubstring!=from-1) {
			DenseSubstrings.computeDensityWindows(ReadA.lastSortedAlignment);
			DenseSubstring.order=DenseSubstring.STARTA;
			if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,from,DenseSubstrings.lastSubstring+1);
			from=cleanNewDenseSubstrings(from);
			

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after cleanNewDenseSubstrings():  DenseSubstrings.lastSubstring="+DenseSubstrings.lastSubstring+" from="+from);
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);	
	IO.printErr("intervals2denseSubstrings> alignment intervals after cleanNewDenseSubstrings():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);		
}			
			
			
			cleanNewDenseSubstrings_periodic(from);
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after cleanNewDenseSubstrings_periodic():");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);		
	IO.printErr("intervals2denseSubstrings> alignment intervals after cleanNewDenseSubstrings_periodic():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);	
}			
			
			
			DenseSubstring.order=DenseSubstring.STARTA;
			if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
		}
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("alignment intervals after AlignmentIntervals.intervals2prefixSuffix():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
	IO.printErr("alignments after AlignmentIntervals.intervals2prefixSuffix():");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());
}

		
		// Concatenating weak dense substrings
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings before concatenateWeakSubstrings():");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);	
	IO.printErr("intervals2denseSubstrings> alignment intervals before concatenateWeakSubstrings():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);		
}		
		from=concatenateWeakSubstrings();
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after concatenateWeakSubstrings():   from="+from);
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);		
	IO.printErr("intervals2denseSubstrings> alignment intervals after concatenateWeakSubstrings():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);	
}

		if (from!=-1) {
			DenseSubstrings.computeDensityWindows(ReadA.lastSortedAlignment);
			from=cleanNewDenseSubstrings(from);
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after cleanNewDenseSubstrings():");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);		
	IO.printErr("intervals2denseSubstrings> alignment intervals after cleanNewDenseSubstrings():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
}
			
			cleanNewDenseSubstrings_periodic(from);
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after cleanNewDenseSubstrings_periodic():");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);		
	IO.printErr("intervals2denseSubstrings> alignment intervals after cleanNewDenseSubstrings_periodic():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
}
			
			cleanNewDenseSubstrings_alignmentIntervals(from);
			
			
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2denseSubstrings> dense substrings after cleanNewDenseSubstrings_alignmentIntervals():");
	for (int i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr(DenseSubstrings.substrings[i]);		
	IO.printErr("intervals2denseSubstrings> alignment intervals after cleanNewDenseSubstrings_alignmentIntervals():");
	for (int i=0; i<=lastInterval; i++) IO.printErr(intervals[i]);
}			
			
			
			DenseSubstring.order=DenseSubstring.STARTA;
			if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
			DenseSubstrings.discardPrefixSuffix_weakSubstring(ReadA.lastSortedAlignment+1);
		}
		else {
			DenseSubstrings.discardPrefixSuffix_weakSubstring(ReadA.lastSortedAlignment+1);
		}
	}


	/**
	 * Sets fields $hasSubstringWithLowCoverage$ and $hasSubstringWithLowQuality$ of every
	 * alignment interval.
	 */
	private static final void markLowQualityIntervals() {
		final int MIN_COVERAGE = IO.minOccurrencesInGenome*IO.coverage-1;
		final int WINDOW = Alignments.minAlignmentLength>>1;  // Arbitrary
		final int QUALITY_WINDOW = WINDOW/Reads.QUALITY_SPACING;
		
		for (int i=0; i<=lastInterval; i++) {
			intervals[i].hasSubstringWithLowCoverage=ReadA.hasSubstringWithLowCoverage(intervals[i].firstPosition,intervals[i].lastPosition,MIN_COVERAGE,WINDOW);
			Reads.readEnds2qualityEnds(ReadA.id,intervals[i].firstPosition,intervals[i].lastPosition,true,tmp);
			intervals[i].hasSubstringWithLowQuality=Histograms.hasSubstringWithAverage(Reads.getQualityArray(ReadA.id),tmp[0],tmp[1],Reads.MIN_RANDOM_QUALITY_SCORE,true,QUALITY_WINDOW,true,-1);
		}
	}
	

	/**
	 * Computes the closure of $intervals2substring_intervalsAreCompatible()$ from each
	 * interval, and creates a new weak dense substring of substring type iff the closure 
	 * contains at least $minDistinctMaximalEvents$ distinct maximal events.
	 * The procedure is a greedy, simplified version of 
	 * $DenseSubstrings.get*Substring_weak()$, and it is useful since: (1) the detection 
	 * of a dense substring of substring type might have failed for several reasons at 
	 * previous stages of the pipeline; (2) all procedures in $DenseSubstrings$ require
	 * alignments to be close to one another to create a dense substring, whereas this
	 * procedure has very weak distance bounds, and can thus detect regions that behave 
	 * like dense substrings but have maximal alignments whose endpoints are very distant 
	 * from one another (such regions occur in practice).
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ to be sorted by 
	 * $startA$, and $intervals$ to be sorted by $firstPosition$. At the end of the 
	 * procedure, the two blocks $DenseSubstrings.substrings[0..X]$ and $DenseSubstrings.
	 * substrings[X+1..lastSubstring]$ are sorted by $startA$, where $X$ is the value of 
	 * $DenseSubstrings.lastSubstring$ before the procedure starts. $intervals$ is also
	 * still sorted by $firstPosition$ when the procedure ends.
	 *
	 * @param tmp1,tmp2 temporary space, of length at least equal to the number of 
	 * alignment intervals. This must not be a pointer to $AlignmentIntervals.stack$.
	 */
	private static final void intervals2substring(int[] tmp1, int[] tmp2) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int DISTANCE_THRESHOLD = Alignments.minAlignmentLength<<1;  // Arbitrary
		final int MIN_DISTINCT_MAXIMAL_EVENTS = IO.coverage*4;  // Arbitrary
		final int MIN_INTERSECTION_LENGTH = (Alignments.minAlignmentLength)>>2;  // Arbitrary
		final int MIN_OVERLAP_LENGTH = IO.quantum<<1;  // Arbitrary
		int i, j;
		int previousLastSubstring, firstJForNext, startA, endA;
		AlignmentInterval oldInterval;
		DenseSubstring firstInRange, tmpSubstring;
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2substring> alignment intervals at the very beginning:");
	for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x]);
	IO.printErr("alignments at the very beginning:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());	
}
		
		// Transforming intervals into dense substrings
		for (i=0; i<=lastInterval; i++) {
			intervals[i].discarded=false;
			intervals[i].denseSubstring=null;
		}
		previousLastSubstring=DenseSubstrings.lastSubstring;
		for (i=0; i<=lastInterval; i++) {
			if (!intervals2substring_useInterval(i)) continue;
			intervals2substring_getConnectedComponent(i,IDENTITY_THRESHOLD,MIN_INTERSECTION_LENGTH,DISTANCE_THRESHOLD,MIN_DISTINCT_MAXIMAL_EVENTS,stack,tmp1,tmp2);
		}
		if (DenseSubstrings.lastSubstring==previousLastSubstring) return;
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("intervals2substring> Z1 new dense substrings created:");
	for (int x=previousLastSubstring+1; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	IO.printErr("intervals2substring> alignment intervals:");
	for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x]+" == denseSubstring ==> "+intervals[x].denseSubstring+" discarded="+intervals[x].discarded);
}			
		
		// Merging straddling dense substrings
		for (i=previousLastSubstring+1; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].representative=null;
		firstInRange=DenseSubstrings.substrings[previousLastSubstring+1];
		for (i=previousLastSubstring+2; i<=DenseSubstrings.lastSubstring; i++) {
			if ((DenseSubstrings.substrings[i].maximalStartA>=0?DenseSubstrings.substrings[i].maximalStartA:DenseSubstrings.substrings[i].startA)>(firstInRange.maximalEndA>0?firstInRange.maximalEndA:firstInRange.endA)-MIN_OVERLAP_LENGTH) {
				firstInRange=DenseSubstrings.substrings[i];
				continue;
			}
			firstInRange.concatenateSubstringType(DenseSubstrings.substrings[i]);
			DenseSubstrings.substrings[i].representative=firstInRange;
		}
		j=previousLastSubstring;
		for (i=previousLastSubstring+1; i<=DenseSubstrings.lastSubstring; i++) {
			if (DenseSubstrings.substrings[i].representative==null) {
				j++;
				tmpSubstring=DenseSubstrings.substrings[j];
				DenseSubstrings.substrings[j]=DenseSubstrings.substrings[i];
				DenseSubstrings.substrings[i]=tmpSubstring;
			}
		}
		DenseSubstrings.lastSubstring=j;
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("intervals2substring> Z1.5 new dense substrings after concatenation:");
	for (int x=previousLastSubstring+1; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	IO.printErr("intervals2substring> Z1.5 alignment intervals after concatenation:");
	for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x]+" == denseSubstring ==> "+intervals[x].denseSubstring+" discarded="+intervals[x].discarded);
}		
		
		// Discarding intervals that are approximately identical to or contained in the 
	 	// result of the merge.
		i=previousLastSubstring+1; j=0; firstJForNext=-1;
		while (i<=DenseSubstrings.lastSubstring) {
			endA=DenseSubstrings.substrings[i].endA;
			if (j>lastInterval || intervals[j].firstPosition>=endA+IDENTITY_THRESHOLD) {
				i++;
				if (firstJForNext!=-1) j=firstJForNext;
				firstJForNext=-1;
				continue;
			}
			startA=DenseSubstrings.substrings[i].startA;
			if ( intervals[j].discarded || intervals[j].lastPosition<=startA-IDENTITY_THRESHOLD || 
			     ( !intervals2substring_markInterval(j) && 
				   !( Math.abs(intervals[j].firstPosition,startA)<=IDENTITY_THRESHOLD && 
					  intervals[j].lastPosition<endA-IDENTITY_THRESHOLD &&
					  intervals[j].nMaximalAlignmentsRight<MIN_INTERVALS_IN_MAXIMAL &&
					  intervals[j].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
				   ) &&
   				   !( Math.abs(intervals[j].lastPosition,endA)<=IDENTITY_THRESHOLD && 
   					  intervals[j].firstPosition>startA+IDENTITY_THRESHOLD &&
   					  intervals[j].nMaximalAlignmentsLeft<MIN_INTERVALS_IN_MAXIMAL &&
					  intervals[j].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
   				   )
				 )
			   ) {
				j++;
				continue;
			}
			if (firstJForNext==-1 && i<DenseSubstrings.lastSubstring && intervals[j].lastPosition>=DenseSubstrings.substrings[i+1].startA) firstJForNext=j;
			if ( intervals2substring_discardInterval(j,startA,endA,IDENTITY_THRESHOLD) ||
			     ( Math.abs(intervals[j].firstPosition,startA)<=IDENTITY_THRESHOLD && 
				   intervals[j].lastPosition<endA-IDENTITY_THRESHOLD &&
				   intervals[j].nMaximalAlignmentsRight<MIN_INTERVALS_IN_MAXIMAL &&
				   intervals[j].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
			     ) ||
				 ( Math.abs(intervals[j].lastPosition,endA)<=IDENTITY_THRESHOLD && 
				   intervals[j].firstPosition>startA+IDENTITY_THRESHOLD &&
				   intervals[j].nMaximalAlignmentsLeft<MIN_INTERVALS_IN_MAXIMAL &&
				   intervals[j].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX
				 )
			   ) {
				intervals[j].discarded=true;
				intervals[j].denseSubstring=DenseSubstrings.substrings[i];
			}
			j++;
		}
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("intervals2substring> Z2 new dense substrings created:");
	for (int x=previousLastSubstring+1; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
	IO.printErr("intervals2substring> alignment intervals after creation:");
	for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x]+" == denseSubstring ==> "+intervals[x].denseSubstring+" discarded="+intervals[x].discarded);
	IO.printErr("alignments after creation:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());	
}		
		
		// Compacting intervals
		j=-1;
		for (i=0; i<=lastInterval; i++) {		
			if (intervals[i].discarded) continue;
			j++;
			if (j!=i) {
				oldInterval=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=oldInterval;
			}
		}
		lastInterval=j;
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2substring> alignment intervals at Z3:");
	for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x]+" == denseSubstring ==> "+intervals[x].denseSubstring+" discarded="+intervals[x].discarded);
}
		
		// Updating pointers from alignments
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			oldInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (oldInterval==null || !oldInterval.discarded) continue;
			ReadA.sortedAlignments[i].inDenseSubstring=oldInterval.denseSubstring.representative==null?oldInterval.denseSubstring:oldInterval.denseSubstring.representative;
			ReadA.sortedAlignments[i].mergedToInterval=null;
		}
		for (i=0; i<=lastInterval; i++) intervals[i].denseSubstring=null;
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("intervals2substring> alignment intervals at the very end:");
	for (int x=0; x<=lastInterval; x++) IO.printErr(intervals[x]);
	IO.printErr("intervals2substring> alignments at the very end:");
	for (int x=0; x<=ReadA.lastSortedAlignment; x++) System.err.println(ReadA.sortedAlignments[x].toStringBoundaries()+" :: "+ReadA.sortedAlignments[x].toStringPointers());	
}		
		
	}
	
	
	/**
	 * @param fromInterval assumed to satisfy $intervals2substring_useInterval()$;
	 * @param minIntersectionLength min. length of an intersection between two intervals;
	 * @param maxOverhangLength max. length of a straddling overhang;
	 * @param stack stack used for computing the closure;
	 * @param usedIntervals temporary space;
	 * @param maximalEvents temporary space.
	 */
	private static final void intervals2substring_getConnectedComponent(int fromInterval, int identityThreshold, int minIntersectionLength, int maxOverhangLength, int minDistinctMaximalEvents, int[] stack, int[] usedIntervals, int[] maximalEvents) {
		final int MIN_INTERVALS = 2;  // Arbitrary
		int i, j;
		int top, lastMaximal, lastUsedInterval, firstPosition, lastPosition, startA, endA;
		int substringStartA, substringEndA, maximalStartA, maximalEndA, nDistinctMaximal;
		DenseSubstring newSubstring;

if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> CONSIDERING FROM ALIGN. INTERVAL "+intervals[fromInterval]);
		// Computing the closure of $intervals2substring_intervalsAreCompatible()$.
		substringStartA=intervals[fromInterval].firstPosition; 
		substringEndA=intervals[fromInterval].lastPosition;
		lastUsedInterval=0; usedIntervals[lastUsedInterval]=fromInterval;
		top=0; stack[top]=fromInterval;
		intervals[fromInterval].discarded=true;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> marked as discarded interval ["+intervals[fromInterval].firstPosition+".."+intervals[fromInterval].lastPosition+"]");
		while (top>=0) {
			i=stack[top--];
			firstPosition=intervals[i].firstPosition; lastPosition=intervals[i].lastPosition;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent>    JUST POPPED FROM THE STACK ["+firstPosition+".."+lastPosition+"]");
			for (j=0; j<=lastInterval; j++) {
				startA=intervals[j].firstPosition;
				if (startA>=lastPosition-minIntersectionLength) break;
				if (j==i || !intervals2substring_useInterval(j)) continue;
				endA=intervals[j].lastPosition;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent>       considering other interval from the one popped from the stack ["+startA+".."+endA+"]");
				if (intervals2substring_intervalsAreCompatible(startA,endA,firstPosition,lastPosition,minIntersectionLength,maxOverhangLength)) {
   					intervals[j].discarded=true;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> marked as discarded interval ["+intervals[j].firstPosition+".."+intervals[j].lastPosition+"]");
					usedIntervals[++lastUsedInterval]=j;
					substringStartA=Math.min(substringStartA,startA);
					substringEndA=Math.max(substringEndA,endA);
					stack[++top]=j;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> (1) discarded interval "+intervals[j]);
				}
			}
		}
		if (lastUsedInterval+1<MIN_INTERVALS) {
			for (i=0; i<=lastUsedInterval; i++) {
				j=usedIntervals[i];
				intervals[j].discarded=false;
			}
			return;
		}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> substringStartA="+substringStartA+" substringEndA="+substringEndA);
		
		// Adding intervals that are approximately identical to or contained in the 
	 	// candidate new dense substring.
		for (i=0; i<=lastInterval; i++) {
			firstPosition=intervals[i].firstPosition; lastPosition=intervals[i].lastPosition;
			if (firstPosition>=substringEndA+identityThreshold) break;
			if (lastPosition<=substringStartA-identityThreshold) continue;
			if ( ( Intervals.areApproximatelyIdentical(firstPosition,lastPosition,substringStartA,substringEndA) ||
				   Intervals.isApproximatelyContained(firstPosition,lastPosition,substringStartA,substringEndA) 
				 ) &&
				 intervals2substring_useInterval(i)
			   ) {
				intervals[i].discarded=true;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> marked as discarded interval ["+intervals[i].firstPosition+".."+intervals[i].lastPosition+"]");
				usedIntervals[++lastUsedInterval]=i;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> (2) discarded interval "+intervals[i]);
			}
		}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> substringStartA="+substringStartA+" substringEndA="+substringEndA);
		
		// Computing the number of maximal events
		lastMaximal=-1;
		substringStartA=Math.POSITIVE_INFINITY; substringEndA=-1;
		maximalStartA=Math.POSITIVE_INFINITY; maximalEndA=-1;
		for (i=0; i<=lastUsedInterval; i++) {
			j=usedIntervals[i];
			startA=intervals[j].firstPosition;
			substringStartA=Math.min(substringStartA,startA);
			if (intervals[j].isLeftMaximal) {
				maximalStartA=Math.min(maximalStartA,startA);
				maximalEvents[++lastMaximal]=startA;
			}
			endA=intervals[j].lastPosition;
			substringEndA=Math.max(substringEndA,endA);
			if (intervals[j].isRightMaximal) {
				maximalEndA=Math.max(maximalEndA,endA);
				maximalEvents[++lastMaximal]=endA;
			}
		}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> lastMaximal="+lastMaximal);
		if (lastMaximal>0) Arrays.sort(maximalEvents,0,lastMaximal+1);
		nDistinctMaximal=1;
		for (j=1; j<=lastMaximal; j++) {
			if (maximalEvents[j]>maximalEvents[j-1]+identityThreshold) nDistinctMaximal++;
		}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> nDistinctMaximal="+nDistinctMaximal+" minDistinctMaximalEvents="+minDistinctMaximalEvents);
		if (nDistinctMaximal<minDistinctMaximalEvents) {
			for (i=0; i<=lastUsedInterval; i++) {
				j=usedIntervals[i];
				intervals[j].discarded=false;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> marked as NOT DISCARDED interval ["+intervals[j].firstPosition+".."+intervals[j].lastPosition+"]");
			}
			return;
		}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> 1");
		
		// Building the new dense substring
		DenseSubstrings.lastSubstring++;
		newSubstring=DenseSubstrings.substrings[DenseSubstrings.lastSubstring];
		newSubstring.isWeak=true;
		newSubstring.id=0;
		newSubstring.startA=substringStartA;
		newSubstring.endA=substringEndA;
		newSubstring.nStartA=1;
		newSubstring.sumStartA=substringStartA;
		newSubstring.nEndA=1;
		newSubstring.sumEndA=substringEndA;
		newSubstring.prefixReplication=false;
		newSubstring.minPrefixLength=-1;
		newSubstring.suffixReplication=false;
		newSubstring.minSuffixLength=-1;
		newSubstring.substringReplication=true;
		newSubstring.singleDeletionReplication=false;
		newSubstring.maximalStartA=maximalStartA;
		newSubstring.maximalEndA=maximalEndA;
		newSubstring.clearImpliedAlignmentsStats();
		newSubstring.clearPrefixSuffixBias();
		newSubstring.period=0;
		newSubstring.strongPeriodSignal=false;
		
		// Marking used intervals, and discarding other intervals that are approximately
		// identical to or contained in the new dense substring.
		for (i=0; i<=lastUsedInterval; i++) {
			j=usedIntervals[i];
			intervals[j].discarded=true;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> marked as discarded interval ["+intervals[j].firstPosition+".."+intervals[j].lastPosition+"]");
			intervals[j].denseSubstring=newSubstring;
		}
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].firstPosition>=substringEndA+identityThreshold) break;
			if (intervals[i].lastPosition<=substringStartA-identityThreshold || !intervals2substring_markInterval(i)) continue;
			if (intervals2substring_discardInterval(i,substringStartA,substringEndA,identityThreshold)) {
				intervals[i].discarded=true;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("intervals2substring_getConnectedComponent> marked as discarded interval ["+intervals[i].firstPosition+".."+intervals[i].lastPosition+"]");
				intervals[i].denseSubstring=newSubstring;
			}
		}
	}
	
	
	private static final boolean intervals2substring_useInterval(int i) {
		return intervals2substring_markInterval(i) && 
			   (intervals[i].isLeftMaximal || intervals[i].isRightMaximal);
	}
	
	
	private static final boolean intervals2substring_markInterval(int i) {
		return intervals2substring_markIntervalPrime(i) && !intervals[i].checkAlignments(MIN_INTERVALS_IN_MAXIMAL,MIN_ALIGNMENTS_LARGE) && intervals[i].nMergedIntervals<MIN_INTERVALS_IN_MAXIMAL_PREFIXSUFFIX;
	}
	
	
	private static final boolean intervals2substring_markIntervalPrime(int i) {
		return !( intervals[i].discarded ||
			      intervals[i].inDenseSubstringOfSubstringType!=null ||
				  intervals[i].hasSubstringWithLowCoverage || 
				  intervals[i].hasSubstringWithLowQuality
			   );
	}

	
	private static final boolean intervals2substring_intervalsAreCompatible(int firstI, int lastI, int firstJ, int lastJ, int minIntersectionLength, int maxOverhangLength) {
		final int COVERAGE_WINDOW = IO.quantum;  // Arbitrary
		final int MIN_COVERAGE = IO.coverage<<2;  // Arbitrary
		final int left = Math.max(firstI,firstJ);
		final int right = Math.min(lastI,lastJ);
		
		return Intervals.straddles(firstI,lastI,firstJ,lastJ) && 
			   (Math.min(lastI,lastJ)-Math.max(firstI,firstJ)+1)>=minIntersectionLength &&
			   Math.abs(firstI,firstJ)<=maxOverhangLength && Math.abs(lastI,lastJ)<=maxOverhangLength &&
			   (left<COVERAGE_WINDOW || ReadA.getAverageCoverage(left-COVERAGE_WINDOW,left-1)>=MIN_COVERAGE) &&
			   (right>=ReadA.readLength-COVERAGE_WINDOW || ReadA.getAverageCoverage(right+1,right+COVERAGE_WINDOW)>=MIN_COVERAGE) &&
			   ReadA.getAverageCoverage(left,right)>=MIN_COVERAGE;
	}
	
	
	/**
	 * Remark: containment rules mirror those in $cleanAlignmentIntervals()$.
	 *
	 * @param startA,endA of the dense substring of substring type;
	 * @param threshold distance threshold;
	 * @return TRUE iff $interval$ should be discarded.
	 */
	private static final boolean intervals2substring_discardInterval(int i, int startA, int endA, int threshold) {
		final int first = intervals[i].firstPosition;
		final int last = intervals[i].lastPosition;

		if (!intervals2substring_markIntervalPrime(i)) return false;
		return Intervals.areApproximatelyIdentical(first,last,startA,endA) ||
			   ( Intervals.isApproximatelyContained(first,last,startA,endA) && 
			     ( (!intervals[i].isLeftMaximal && !intervals[i].isRightMaximal) ||
				   !intervals[i].checkAlignments(MIN_INTERVALS_IN_MAXIMAL,MIN_ALIGNMENTS_LARGE) ||
			       ( !intervals[i].isLeftMaximal && intervals[i].isRightMaximal && 
			         Math.abs(last,endA)<=threshold
			       ) || 
			       ( intervals[i].isLeftMaximal && !intervals[i].isRightMaximal && 
			         Math.abs(first,startA)<=threshold
			       )
			     )
			   );
	}


	/**
	 * The procedure discards alignment intervals that are contained in a dense substring
	 * of substring type, and that trace a long prefix or suffix substring. The procedure 
	 * is a simplified version of $getSubstrings()$, which works on alignment intervals 
	 * rather than on alignments.
	 *	
	 * Remark: the procedure assumes that $markIntervalsInDenseSubstrings()$ has already
	 * been executed. At the end of the procedure, $intervals$ is sorted by 
	 * $lastPosition$.
	 * 
	 * @param componentSize temporary space, of length at least equal to the number of 
	 * alignment intervals.
	 */
	private static final void intervals2prefixSuffix_inSubstring(int[] componentSize) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MAX_DISTANCE_THRESHOLD = DenseSubstrings.maxDistanceEndWeak;
		final int MERGED_INTERVALS_THRESHOLD = IO.coverage*20;  // Arbitrary
		final int COVERAGE_WINDOW = IO.quantum;  // Arbitrary
		final int MIN_COVERAGE = IO.coverage<<1;  // Arbitrary
		boolean atLeastOneDiscarded;
		int i, j;
		int nIntervals, nComponents, largestComponent, destination;
		int start1, end1, start2, end2;
		AlignmentInterval oldInterval;
		
		if (lastInterval==-1) return;
		nIntervals=lastInterval+1;
		for (i=0; i<=lastInterval; i++) intervals[i].discarded=false;
		atLeastOneDiscarded=false;
		
		// Building the prefix DAG, and sorting it topologically.
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (nIntervals>1) Arrays.sort(intervals,0,nIntervals);
		}
		clearNeighborMatrices(nIntervals);
		for (i=0; i<=lastInterval; i++) {
			if ( intervals[i].inDenseSubstringOfSubstringType==null || 
				 intervals[i].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) ||
				 intervals[i].hasSubstringWithLowCoverage || 
				 intervals[i].hasSubstringWithLowQuality
			   ) continue;
			start1=intervals[i].firstPosition;
			end1=intervals[i].lastPosition;
			for (j=i-1; j>=0; j--) {
				if ( intervals[j].inDenseSubstringOfSubstringType==null || 
					 intervals[j].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) ||
					 intervals[j].hasSubstringWithLowCoverage || 
					 intervals[j].hasSubstringWithLowQuality
				   ) continue;
				start2=intervals[j].firstPosition;
				if (start2<start1-IDENTITY_THRESHOLD) break;
				end2=intervals[j].lastPosition;
				if (end2<=end1 || end2>end1+MAX_DISTANCE_THRESHOLD) continue;				
				if (ReadA.getAverageCoverage(end1,Math.min(end1+COVERAGE_WINDOW-1,ReadA.readLength-1))<MIN_COVERAGE) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
			for (j=i+1; j<=lastInterval; j++) {
				if ( intervals[j].inDenseSubstringOfSubstringType==null || 
					 intervals[j].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) ||
					 intervals[j].hasSubstringWithLowCoverage || 
					 intervals[j].hasSubstringWithLowQuality
				   ) continue;
				start2=intervals[j].firstPosition;
				if (start2>start1+IDENTITY_THRESHOLD) break;
				end2=intervals[j].lastPosition;
				if (end2<=end1 || end2>end1+MAX_DISTANCE_THRESHOLD) continue;
				if (ReadA.getAverageCoverage(end1,Math.min(end1+COVERAGE_WINDOW-1,ReadA.readLength-1))<MIN_COVERAGE) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(nIntervals,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stack);
		if (nComponents<=0 || nComponents>nIntervals) {
			IO.printCriticalErr("Error while detecting substrings from alignment intervals: wrong number of connected components.");
			System.exit(1);
		}
		largestComponent=DAG.largestComponent(components,nIntervals,componentSize,nComponents);
		if (largestComponent>=MIN_PATH_LENGTH) {
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nIntervals);
			i=DAG.topologicalSort(nIntervals,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting substrings from alignment intervals: the prefix DAG of intervals contains a cycle that involves node "+(i-1));
				System.exit(1);
			}
			i=0;
			while (i<=lastInterval) {
				if (intervals[i].discarded || componentSize[components[i]]<MIN_PATH_LENGTH) {
					i++;
					continue;
				}
				atLeastOneDiscarded|=intervals2prefixSuffix_getSubstring(true,i,nIntervals,MIN_PATH_LENGTH,false)!=-1;
				i++;
			}
		}
		
		// Building the suffix DAG, and sorting it topologically.
		if (AlignmentInterval.order!=AlignmentInterval.LASTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.LASTPOSITION;
			if (nIntervals>1) Arrays.sort(intervals,0,nIntervals);
		}
		clearNeighborMatrices(nIntervals);
		for (i=0; i<=lastInterval; i++) {
			if ( intervals[i].inDenseSubstringOfSubstringType==null || 
				 intervals[i].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) ||
				 intervals[i].hasSubstringWithLowCoverage || 
				 intervals[i].hasSubstringWithLowQuality
			   ) continue;
			start1=intervals[i].firstPosition;
			end1=intervals[i].lastPosition;
			for (j=i-1; j>=0; j--) {
				if ( intervals[j].inDenseSubstringOfSubstringType==null || 
					 intervals[j].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) ||
					 intervals[j].hasSubstringWithLowCoverage || 
					 intervals[j].hasSubstringWithLowQuality
				   ) continue;
				end2=intervals[j].lastPosition;
				if (end2<end1-IDENTITY_THRESHOLD) break;
				start2=intervals[j].firstPosition;
				if (start2<=start1 || start2>start1+MAX_DISTANCE_THRESHOLD) continue;
				if (ReadA.getAverageCoverage(Math.max(start2-COVERAGE_WINDOW+1,0),start2-1)<MIN_COVERAGE) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
			for (j=i+1; j<=lastInterval; j++) {
				if ( intervals[j].inDenseSubstringOfSubstringType==null || 
					 intervals[j].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) ||
					 intervals[j].hasSubstringWithLowCoverage || 
					 intervals[j].hasSubstringWithLowQuality
				   ) continue;
				end2=intervals[j].lastPosition;
				if (end2>end1+IDENTITY_THRESHOLD) break;
				start2=intervals[j].firstPosition;
				if (start2<=start1 || start2>start1+MAX_DISTANCE_THRESHOLD) continue;
				if (ReadA.getAverageCoverage(Math.max(start2-COVERAGE_WINDOW+1,0),start2-1)<MIN_COVERAGE) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(nIntervals,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stack);
		if (nComponents<=0 || nComponents>nIntervals) {
			IO.printCriticalErr("Error while detecting substrings from alignment intervals: wrong number of connected components.");
			System.exit(1);
		}
		largestComponent=DAG.largestComponent(components,nIntervals,componentSize,nComponents);
		if (largestComponent>=MIN_PATH_LENGTH) {
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nIntervals);
			i=DAG.topologicalSort(nIntervals,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting substrings from alignment intervals: the suffix DAG of intervals contains a cycle that involves node "+(i-1));
				System.exit(1);
			}
			i=0;
			while (i<=lastInterval) {
				if (intervals[i].discarded || componentSize[components[i]]<MIN_PATH_LENGTH) {
					i++;
					continue;
				}
				atLeastOneDiscarded|=intervals2prefixSuffix_getSubstring(false,i,nIntervals,MIN_PATH_LENGTH,false)!=-1;
				i++;
			}
		}
		if (!atLeastOneDiscarded) return;
		
		// Compacting intervals
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].discarded) continue;
			j++;
			if (j!=i) {
				oldInterval=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=oldInterval;
			}
		}
		lastInterval=j;
		
		// Updating pointers from alignments
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			oldInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (oldInterval==null || !oldInterval.discarded) continue;
			ReadA.sortedAlignments[i].inDenseSubstring=oldInterval.inDenseSubstringOfSubstringType;
			ReadA.sortedAlignments[i].mergedToInterval=null;
		}
	}
	
	
	/**
	 * Same as $intervals2prefixSuffix_inSubstring()$, but for intervals that are not 
	 * contained in a dense substring of substring type. The procedure requires them to 
	 * have few assigned alignments, and it creates corresponding weak dense substrings.
	 *
	 * Remark: at the end of the procedure, block $DenseSubstrings.substrings[0..X]$ is
	 * still sorted in the same order as at the beginning, but $DenseSubstrings.
	 * substrings[X+1..lastSubstring]$ is not sorted in any order (X is the value of 
	 * $DenseSubstrings.lastSubstring$ before the procedure is called).
	 */
	private static final void intervals2prefixSuffix(int[] componentSize) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MAX_DISTANCE_THRESHOLD = DenseSubstrings.maxDistanceEndWeak;
		final int MERGED_INTERVALS_THRESHOLD = MIN_INTERVALS_IN_NONMAXIMAL;  // Arbitrary
		final int COVERAGE_WINDOW = IO.quantum;  // Arbitrary
		final int MIN_COVERAGE = IO.coverage<<1;  // Arbitrary
		int i, j;
		int nIntervals, nComponents, largestComponent, destination;
		int start1, end1, start2, end2, previousLastSubstring;
		AlignmentInterval oldInterval;
		
		if (lastInterval==-1) return;
		nIntervals=lastInterval+1;
		for (i=0; i<=lastInterval; i++) {
			intervals[i].discarded=false;
			intervals[i].denseSubstring=null;
		}
		previousLastSubstring=DenseSubstrings.lastSubstring;
		
		// Building the prefix DAG, and sorting it topologically.
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (nIntervals>1) Arrays.sort(intervals,0,nIntervals);
		}
		clearNeighborMatrices(nIntervals);
		for (i=0; i<=lastInterval; i++) {
			if ( intervals[i].inDenseSubstringOfSubstringType!=null || 
				 intervals[i].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) || 
				 intervals[i].hasSubstringWithLowCoverage || 
				 intervals[i].hasSubstringWithLowQuality ||
				 !intervals[i].isRightMaximal	 
			   ) continue;
			start1=intervals[i].firstPosition;
			end1=intervals[i].lastPosition;
			for (j=i-1; j>=0; j--) {
				if ( intervals[j].inDenseSubstringOfSubstringType!=null || 
					 intervals[j].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) || 
					 intervals[j].hasSubstringWithLowCoverage || 
					 intervals[j].hasSubstringWithLowQuality ||
					 !intervals[j].isRightMaximal
				   ) continue;
				start2=intervals[j].firstPosition;
				if (start2<start1-IDENTITY_THRESHOLD) break;
				end2=intervals[j].lastPosition;
				if (end2<end1 || end2>end1+MAX_DISTANCE_THRESHOLD) continue;
				if (ReadA.getAverageCoverage(end1,Math.min(end1+COVERAGE_WINDOW-1,ReadA.readLength-1))<MIN_COVERAGE) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
			for (j=i+1; j<=lastInterval; j++) {
				if ( intervals[j].inDenseSubstringOfSubstringType!=null || 
					 intervals[j].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) || 
					 intervals[j].hasSubstringWithLowCoverage || 
					 intervals[j].hasSubstringWithLowQuality ||
					 !intervals[j].isRightMaximal
				   ) continue;
				start2=intervals[j].firstPosition;
				if (start2>start1+IDENTITY_THRESHOLD) break;
				end2=intervals[j].lastPosition;
				if (end2<end1 || end2>end1+MAX_DISTANCE_THRESHOLD) continue;
				if (ReadA.getAverageCoverage(end1,Math.min(end1+COVERAGE_WINDOW-1,ReadA.readLength-1))<MIN_COVERAGE) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(nIntervals,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stack);
		if (nComponents<=0 || nComponents>nIntervals) {
			IO.printCriticalErr("Error while detecting substrings from alignment intervals: wrong number of connected components.");
			System.exit(1);
		}
		largestComponent=DAG.largestComponent(components,nIntervals,componentSize,nComponents);
		if (largestComponent>=MIN_PATH_LENGTH) {
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nIntervals);
			i=DAG.topologicalSort(nIntervals,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting substrings from alignment intervals: the prefix DAG of intervals contains a cycle that involves node "+(i-1));
				System.exit(1);
			}
			i=0;
			while (i<=lastInterval) {
				if (intervals[i].discarded || componentSize[components[i]]<MIN_PATH_LENGTH) {
					i++;
					continue;
				}
				intervals2prefixSuffix_getSubstring(true,i,nIntervals,MIN_PATH_LENGTH,true);
				i++;
			}
		}
		
		// Building the suffix DAG, and sorting it topologically.
		if (AlignmentInterval.order!=AlignmentInterval.LASTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.LASTPOSITION;
			if (nIntervals>1) Arrays.sort(intervals,0,nIntervals);
		}
		clearNeighborMatrices(nIntervals);
		for (i=0; i<=lastInterval; i++) {
			if ( intervals[i].inDenseSubstringOfSubstringType!=null || 
				 intervals[i].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) || 
				 intervals[i].hasSubstringWithLowCoverage || 
				 intervals[i].hasSubstringWithLowQuality ||
				 !intervals[i].isLeftMaximal
			   ) continue;
			start1=intervals[i].firstPosition;
			end1=intervals[i].lastPosition;
			for (j=i-1; j>=0; j--) {
				if ( intervals[j].inDenseSubstringOfSubstringType!=null || 
					 intervals[j].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) || 
					 intervals[j].hasSubstringWithLowCoverage || 
					 intervals[j].hasSubstringWithLowQuality ||
					 !intervals[j].isLeftMaximal
				   ) continue;
				end2=intervals[j].lastPosition;
				if (end2<end1-IDENTITY_THRESHOLD) break;
				start2=intervals[j].firstPosition;
				if (start2<start1 || start2>start1+MAX_DISTANCE_THRESHOLD) continue;
				if (ReadA.getAverageCoverage(Math.max(start2-COVERAGE_WINDOW+1,0),start2-1)<MIN_COVERAGE) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
			for (j=i+1; j<=lastInterval; j++) {
				if ( intervals[j].inDenseSubstringOfSubstringType!=null || 
					 intervals[j].checkAlignments(MERGED_INTERVALS_THRESHOLD,MIN_ALIGNMENTS_LARGE) || 
					 intervals[j].hasSubstringWithLowCoverage || 
					 intervals[j].hasSubstringWithLowQuality ||
					 !intervals[j].isLeftMaximal
				   ) continue;
				end2=intervals[j].lastPosition;
				if (end2>end1+IDENTITY_THRESHOLD) break;
				start2=intervals[j].firstPosition;
				if (start2<start1 || start2>start1+MAX_DISTANCE_THRESHOLD) continue;
				if (ReadA.getAverageCoverage(Math.max(start2-COVERAGE_WINDOW+1,0),start2-1)<MIN_COVERAGE) continue;
				Math.ensureSpace(outNeighbors,i,nOutNeighbors[i]+2);
				outNeighbors[i][nOutNeighbors[i]++]=j;
				Math.ensureSpace(inNeighbors,j,nInNeighbors[j]+2);
				inNeighbors[j][nInNeighbors[j]++]=i;
			}
		}
		nComponents=DAG.getConnectedComponents(nIntervals,inNeighbors,nInNeighbors,outNeighbors,nOutNeighbors,components,stack);
		if (nComponents<=0 || nComponents>nIntervals) {
			IO.printCriticalErr("Error while detecting substrings from alignment intervals: wrong number of connected components.");
			System.exit(1);
		}
		largestComponent=DAG.largestComponent(components,nIntervals,componentSize,nComponents);
		if (largestComponent>=MIN_PATH_LENGTH) {
			System.arraycopy(nInNeighbors,0,nInNeighborsPrime,0,nIntervals);
			i=DAG.topologicalSort(nIntervals,inNeighbors,nInNeighborsPrime,outNeighbors,nOutNeighbors,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (i!=0) {
				IO.printCriticalErr("Error while detecting substrings from alignment intervals: the suffix DAG of intervals contains a cycle that involves node "+(i-1));
				System.exit(1);
			}
			i=0;
			while (i<=lastInterval) {
				if (intervals[i].discarded || componentSize[components[i]]<MIN_PATH_LENGTH) {
					i++;
					continue;
				}
				intervals2prefixSuffix_getSubstring(false,i,nIntervals,MIN_PATH_LENGTH,true);
				i++;
			}
		}
		if (DenseSubstrings.lastSubstring==previousLastSubstring) return;
		
		// Compacting intervals
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].discarded) continue;
			j++;
			if (j!=i) {
				oldInterval=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=oldInterval;
			}
		}
		lastInterval=j;
		
		// Updating pointers from alignments
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			oldInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (oldInterval==null || !oldInterval.discarded) continue;
			if (ReadA.sortedAlignments[i].inDenseSubstring==null) ReadA.sortedAlignments[i].impliedByDenseSubstring=oldInterval.denseSubstring;
			ReadA.sortedAlignments[i].mergedToInterval=null;
		}
		for (i=0; i<=lastInterval; i++) intervals[i].denseSubstring=null;
	}
	
	
	/**
	 * Remark: the procedure is a simplified version of $DenseSubstrings.get*Substring()$,
	 * which runs on alignment intervals rather than on alignments.
	 *
	 * @return the position in $intervals$ of an interval that maximizes the distance from
	 * $source$ in the DAG, and whose ending position is the rightmost (if 
	 * $prefixOrSuffix=TRUE$), or whose starting position is the leftmost (if 
	 * $prefixOrSuffix=FALSE$). The procedure marks as discarded all intervals in a 
	 * longest path from $source$ to the selected interval, and it returns -1 if no path 
	 * from $source$ contains at least $minPathLength$ vertices.
	 * @param createDense TRUE: creates a new weak dense substring corresponding to the
	 * intervals in a longest path, and it assigns it to field $denseSubstring$ of such
	 * alignment intervals.
	 */
	private static final int intervals2prefixSuffix_getSubstring(boolean prefixOrSuffix, int source, int nVertices, int minPathLength, boolean createDense) {
		int i, j, k, s;
		int sourceComponent, vertex, selectedInterval, distance, longestDistance;
		int start, end, leftmostStart, rightmostEnd, firstPosition, lastPosition;
		int lastSortedVertex, startA, endA, maximalStartA, maximalEndA;
		DenseSubstring newSubstring;

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

		// Computing an interval whose distance from $source$ in the DAG is maximum, and
		// whose ending position is the rightmost (respectively, whose starting position
		// is the leftmost).
		longestDistance=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]>longestDistance) longestDistance=distances[vertex];
		}
		if (longestDistance+1<minPathLength) return -1;
		selectedInterval=-1; leftmostStart=Math.POSITIVE_INFINITY; rightmostEnd=-1;
		for (i=s+1; i<=lastSortedVertex; i++) {
			vertex=sorted2original[i];
			if (distances[vertex]<longestDistance) continue;
			start=intervals[vertex].firstPosition;
			if (start<leftmostStart) {
				leftmostStart=start;
				if (!prefixOrSuffix) selectedInterval=vertex;
			}
			end=intervals[vertex].lastPosition;
			if (end>rightmostEnd) {
				rightmostEnd=end;
				if (prefixOrSuffix) selectedInterval=vertex;
			}
		}
		
		// Discarding all intervals in a longest path to $selectedInterval$, and creating
		// a new dense substring if needed.
		if (createDense) {
			DenseSubstrings.lastSubstring++;
			newSubstring=DenseSubstrings.substrings[DenseSubstrings.lastSubstring];
			newSubstring.isWeak=true;
			newSubstring.id=0;
			newSubstring.prefixReplication=prefixOrSuffix;
			newSubstring.minPrefixLength=-1;
			newSubstring.suffixReplication=!prefixOrSuffix;
			newSubstring.minSuffixLength=-1;
			newSubstring.substringReplication=false;
			newSubstring.singleDeletionReplication=false;
			newSubstring.clearImpliedAlignmentsStats();
			newSubstring.clearPrefixSuffixBias();
			newSubstring.period=0;
			newSubstring.strongPeriodSignal=false;
			newSubstring.hasLengthBoundariesPrefix=false;
			newSubstring.hasLengthBoundariesSuffix=false;
		}
		else newSubstring=null;
		startA=Math.POSITIVE_INFINITY; endA=0;
		maximalStartA=Math.POSITIVE_INFINITY; maximalEndA=-1;
		vertex=selectedInterval;
		while (vertex!=-1) {
			intervals[vertex].discarded=true;
			intervals[vertex].denseSubstring=newSubstring;
			firstPosition=intervals[vertex].firstPosition;
			lastPosition=intervals[vertex].lastPosition;
			if (firstPosition<startA) startA=firstPosition;
			if (lastPosition>endA) endA=lastPosition;
			if (intervals[vertex].isLeftMaximal && firstPosition<maximalStartA) maximalStartA=firstPosition;
			if (intervals[vertex].isRightMaximal && lastPosition>maximalEndA) maximalEndA=lastPosition;
			vertex=predecessors[vertex];
		}
		if (createDense) {
			newSubstring.startA=startA;
			newSubstring.endA=endA;
			newSubstring.nStartA=1;
			newSubstring.sumStartA=startA;
			newSubstring.nEndA=1;
			newSubstring.sumEndA=endA;
			newSubstring.maximalStartA=maximalStartA;
			newSubstring.maximalEndA=maximalEndA;
		}
		return selectedInterval;
	}
	
	
	/**
	 * Makes sure that newly-created dense substrings are consistent with existing dense 
	 * substrings.
	 *
	 * Remark: this procedure enforces just basic consistency criteria. In particular, it 
	 * does not concatenate weak dense substrings, as in $DenseSubstrings.
	 * concatenateWeakSubstrings()$.
	 *
	 * Remark: if an old weak dense substring of substring type is contained in a new one,
	 * it is merged to the new one iff the density of maximal events is not too large
	 * compared to the new one. Thus, after this procedure completes, a weak dense 
	 * substring of substring type might be contained in another weak dense substring of
	 * substring type. This is contrary to the logic enforced in all steps of the pipeline
	 * before this one, but it is useful in practice.
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings[0..from-1]$ and
	 * $DenseSubstrings.substrings[from..]$ to be sorted by $startA$.
	 *
	 * @param from $DenseSubstrings.substrings[from..]$ is assumed to contain all and only
	 * the new weak substrings;
	 * @return position X, where block $DenseSubstrings.substrings[X..]$ contains all and 
	 * only the surviving new weak dense substrings; both blocks $DenseSubstrings.
	 * substrings[0..X-1]$ and $DenseSubstrings.substrings[X..]$ are internally sorted by 
	 * $startA$.
	 */
	private static final int cleanNewDenseSubstrings(int from) {
		final int DENSITY_THRESHOLD = 5;  // Arbitrary
		boolean atLeastOneDiscarded;
		int i, j;
		int startI, endI, startJ, endJ, firstJForNextI, representative;
		int previousLastSubstring, newFrom;
		double densityI, densityJ;
		DenseSubstring tmpSubstring;
		
		Factorize.computeMaximalAlignments(true,false,false,ReadA.lastSortedAlignment);
		previousLastSubstring=DenseSubstrings.lastSubstring;
		
		// Discarding new weak substrings that are contained in or identical to old (weak
		// and non-weak) substrings.
		for (i=from; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].representative=null;
		i=from; j=0; firstJForNextI=-1; atLeastOneDiscarded=false;
		startI=DenseSubstrings.substrings[i].startA;
		endI=DenseSubstrings.substrings[i].endA;
		startJ=DenseSubstrings.substrings[j].startA;
		endJ=DenseSubstrings.substrings[j].endA;
		representative=-1;
		while (i<=DenseSubstrings.lastSubstring) {
			if (j==from || startJ>=endI) {
				if (representative!=-1) {
					DenseSubstrings.substrings[i].representative=DenseSubstrings.substrings[representative];
					atLeastOneDiscarded=true;
				}
				i++;
				if (i<=DenseSubstrings.lastSubstring) {
					startI=DenseSubstrings.substrings[i].startA;
					endI=DenseSubstrings.substrings[i].endA;
				}
				representative=-1;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					if (j<from) {
						startJ=DenseSubstrings.substrings[j].startA;
						endJ=DenseSubstrings.substrings[j].endA;
					}
				}
				firstJForNextI=-1;
				continue;
			}
			if (endJ<startI) {
				j++;
				if (j<from) {
					startJ=DenseSubstrings.substrings[j].startA;
					endJ=DenseSubstrings.substrings[j].endA;
				}
				continue;
			}
			if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && endJ>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
			if (isDiscarded(i,j,startI,endI,-1,startJ,endJ,-1,-1,true)) representative=j;
			j++;
			if (j<from) {
				startJ=DenseSubstrings.substrings[j].startA;
				endJ=DenseSubstrings.substrings[j].endA;
			}
		}
		if (atLeastOneDiscarded) {
			j=from-1;
			for (i=from; i<=DenseSubstrings.lastSubstring; i++) {
				if (DenseSubstrings.substrings[i].representative!=null) continue;
				j++;
				if (j!=i) {
					tmpSubstring=DenseSubstrings.substrings[j];
					DenseSubstrings.substrings[j]=DenseSubstrings.substrings[i];
					DenseSubstrings.substrings[i]=tmpSubstring;
				}	
			}
			DenseSubstrings.lastSubstring=j;
		}
		
		// Discarding old weak substrings that are contained in or identical to new ones,
		// and old non-weak substrings that are identical to new ones.
		for (i=0; i<from; i++) DenseSubstrings.substrings[i].representative=null;
		i=0; j=from; firstJForNextI=-1; atLeastOneDiscarded=false;
		startI=DenseSubstrings.substrings[i].startA;
		endI=DenseSubstrings.substrings[i].endA;
		densityI=((double)DenseSubstrings.substrings[i].nMaximalAlignments)/Math.max(DenseSubstrings.substrings[i].length()-Alignments.minAlignmentLength,1);
		startJ=DenseSubstrings.substrings[j].startA;
		endJ=DenseSubstrings.substrings[j].endA;
		densityJ=((double)DenseSubstrings.substrings[j].nMaximalAlignments)/Math.max(DenseSubstrings.substrings[j].length()-Alignments.minAlignmentLength,1);
		representative=-1;
		while (i<from) {
			if (j>DenseSubstrings.lastSubstring || startJ>=endI) {
				if (representative!=-1) {
					DenseSubstrings.substrings[i].representative=DenseSubstrings.substrings[representative];
					DenseSubstrings.substrings[representative].hasLengthBoundariesPrefix|=DenseSubstrings.substrings[i].hasLengthBoundariesPrefix;
					DenseSubstrings.substrings[representative].hasLengthBoundariesSuffix|=DenseSubstrings.substrings[i].hasLengthBoundariesSuffix;
					atLeastOneDiscarded=true;
				}
				i++;
				if (i<from) {
					startI=DenseSubstrings.substrings[i].startA;
					endI=DenseSubstrings.substrings[i].endA;
					densityI=((double)DenseSubstrings.substrings[i].nMaximalAlignments)/Math.max(DenseSubstrings.substrings[i].length()-Alignments.minAlignmentLength,1);
				}
				representative=-1;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					startJ=DenseSubstrings.substrings[j].startA;
					endJ=DenseSubstrings.substrings[j].endA;
					densityJ=((double)DenseSubstrings.substrings[j].nMaximalAlignments)/Math.max(DenseSubstrings.substrings[j].length()-Alignments.minAlignmentLength,1);
				}
				firstJForNextI=-1;
				continue;
			}
			if (endJ<startI) {
				j++;
				if (j<=DenseSubstrings.lastSubstring) {
					startJ=DenseSubstrings.substrings[j].startA;
					endJ=DenseSubstrings.substrings[j].endA;
					densityJ=((double)DenseSubstrings.substrings[j].nMaximalAlignments)/Math.max(DenseSubstrings.substrings[j].length()-Alignments.minAlignmentLength,1);
				}
				continue;
			}
			if (firstJForNextI==-1 && i<from-1 && endJ>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
			if (isDiscarded(i,j,startI,endI,densityI,startJ,endJ,densityJ,DENSITY_THRESHOLD,DenseSubstrings.substrings[i].isWeak||DenseSubstrings.substrings[i].substringReplication)) {
				representative=j;
			}
			j++;
			if (j<=DenseSubstrings.lastSubstring) {
				startJ=DenseSubstrings.substrings[j].startA;
				endJ=DenseSubstrings.substrings[j].endA;
				densityJ=((double)DenseSubstrings.substrings[j].nMaximalAlignments)/Math.max(DenseSubstrings.substrings[j].length()-Alignments.minAlignmentLength,1);
			}
		}
		newFrom=from;
		if (atLeastOneDiscarded) {
			j=-1;
			for (i=0; i<from; i++) {
				if (DenseSubstrings.substrings[i].representative!=null) continue;
				j++;
				if (j!=i) {
					tmpSubstring=DenseSubstrings.substrings[j];
					DenseSubstrings.substrings[j]=DenseSubstrings.substrings[i];
					DenseSubstrings.substrings[i]=tmpSubstring;
				}
			}
			newFrom=j+1;
			for (i=from; i<=DenseSubstrings.lastSubstring; i++) {
				j++;
				tmpSubstring=DenseSubstrings.substrings[j];
				DenseSubstrings.substrings[j]=DenseSubstrings.substrings[i];
				DenseSubstrings.substrings[i]=tmpSubstring;
			}
			DenseSubstrings.lastSubstring=j;
		}
		
		// Updating pointers from alignments
		if (DenseSubstrings.lastSubstring!=previousLastSubstring) {
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
				if (tmpSubstring!=null && tmpSubstring.representative!=null) {
					tmpSubstring=tmpSubstring.representative;
					if (tmpSubstring.substringReplication) {
						ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
						ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring;
					}
					else ReadA.sortedAlignments[i].impliedByDenseSubstring=tmpSubstring;
				}
				tmpSubstring=ReadA.sortedAlignments[i].inDenseSubstring;
				if (tmpSubstring!=null && tmpSubstring.representative!=null) ReadA.sortedAlignments[i].inDenseSubstring=tmpSubstring.representative;
			}
			for (i=0; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].representative=null;
		}
		
		return newFrom;
	}
	
	
	/**
	 * Discards substring $i$ if it is identical to or contained in substring $j$, and if
	 * its average density of maximal events, or its number of maximal events at the 
	 * boundaries, is not large compared to substring $j$.
	 *
	 * @param densityThreshold if negative, the procedure does not use densities and 
	 * maximal events at the boundaries to discard substring $i$;
	 * @param useContainment TRUE iff containment should also be checked, in addition to
	 * identity;
	 * @return TRUE iff $DenseSubstrings.substrings[i]$ is discarded by 
	 * $DenseSubstrings.substrings[j]$.
	 */
	private static final boolean isDiscarded(int i, int j, int startI, int endI, double densityI, int startJ, int endJ, double densityJ, double densityThreshold, boolean useContainment) {
		final double WINDOWS_DENSITY_THRESHOLD = 2.0;  // Arbitrary
		
		if ( DenseSubstrings.substrings[i].substringReplication &&
			 DenseSubstrings.substrings[j].substringReplication &&
			 ( Intervals.areApproximatelyIdentical_lowQuality(startI,endI,startJ,endJ,ReadA.id) || 
			   ( useContainment && 
				 Intervals.isApproximatelyContained_lowQuality(startI,endI,startJ,endJ,ReadA.id) && 
				 ( densityThreshold<0?true:(densityI<=densityJ*densityThreshold && isDiscarded_windows(i,j,WINDOWS_DENSITY_THRESHOLD)) )
			   )
		     )
		   ) return true;
		else if ( ((DenseSubstrings.substrings[i].prefixReplication && DenseSubstrings.substrings[i].suffixReplication) || DenseSubstrings.substrings[i].singleDeletionReplication) &&
			      ( ( (DenseSubstrings.substrings[j].prefixReplication && DenseSubstrings.substrings[j].suffixReplication) || 
					   DenseSubstrings.substrings[j].singleDeletionReplication || 
					   DenseSubstrings.substrings[j].substringReplication
					) &&
			 	    ( Intervals.areApproximatelyIdentical_lowQuality(startI,endI,startJ,endJ,ReadA.id) || 
			          ( useContainment &&
					    Intervals.isApproximatelyContained_lowQuality(startI,endI,startJ,endJ,ReadA.id) && 
					    ( densityThreshold<0?true:(densityI<=densityJ*densityThreshold && isDiscarded_windows(i,j,WINDOWS_DENSITY_THRESHOLD)) )
					  )
				    )
			      )
		        ) return true;
		else if ( DenseSubstrings.substrings[i].prefixReplication &&
			      ( ( (DenseSubstrings.substrings[j].prefixReplication || DenseSubstrings.substrings[j].singleDeletionReplication) &&
			 	      ( Intervals.areApproximatelyIdentical_lowQuality(startI,endI,startJ,endJ,ReadA.id) || 
			            ( useContainment &&
						  Intervals.isApproximatelyContained_lowQuality(startI,endI,startJ,endJ,ReadA.id) && 
						  !DenseSubstrings.substrings[i].isLeftMaximal && 
						  ( densityThreshold<0?true:(densityI<=densityJ*densityThreshold && isDiscarded_windows(i,j,WINDOWS_DENSITY_THRESHOLD)) )
						)
				      )
			        ) ||
			        ( DenseSubstrings.substrings[j].substringReplication &&
				      ( Intervals.areApproximatelyIdentical_lowQuality(startI,endI,startJ,endJ,ReadA.id) ||
				        ( useContainment &&
						  Intervals.isApproximatelyContained_lowQuality(startI,endI,startJ,endJ,ReadA.id) && 
						  ( densityThreshold<0?true:(densityI<=densityJ*densityThreshold && isDiscarded_windows(i,j,WINDOWS_DENSITY_THRESHOLD)) )
						)
				      )
			        )
			      )
		        ) return true;
		else if ( DenseSubstrings.substrings[i].suffixReplication &&
			      ( ( (DenseSubstrings.substrings[j].suffixReplication || DenseSubstrings.substrings[j].singleDeletionReplication) &&
			 	      ( Intervals.areApproximatelyIdentical_lowQuality(startI,endI,startJ,endJ,ReadA.id) || 
			            ( useContainment &&
						  Intervals.isApproximatelyContained_lowQuality(startI,endI,startJ,endJ,ReadA.id) && 
						  !DenseSubstrings.substrings[i].isRightMaximal && 
						  ( densityThreshold<0?true:(densityI<=densityJ*densityThreshold && isDiscarded_windows(i,j,WINDOWS_DENSITY_THRESHOLD)) )
						)
				      )
			        ) ||
			        ( DenseSubstrings.substrings[j].substringReplication &&
				      ( Intervals.areApproximatelyIdentical_lowQuality(startI,endI,startJ,endJ,ReadA.id) ||
				        ( useContainment &&
						  Intervals.isApproximatelyContained_lowQuality(startI,endI,startJ,endJ,ReadA.id) && 
						  ( densityThreshold<0?true:(densityI<=densityJ*densityThreshold && isDiscarded_windows(i,j,WINDOWS_DENSITY_THRESHOLD)) )
						)
				      )
			        )
			      )
		        ) return true;
		return false;
	}
	
	
	private static final boolean isDiscarded_windows(int i, int j, double threshold) {
		final int IDENTITY_THRESHOLD = IO.quantum;

		return !(DenseSubstrings.substrings[i].startA>DenseSubstrings.substrings[j].startA+IDENTITY_THRESHOLD && DenseSubstrings.substrings[i].leftWindow>DenseSubstrings.substrings[j].maxWindow*threshold) &&
			   !(DenseSubstrings.substrings[i].endA<DenseSubstrings.substrings[j].endA-IDENTITY_THRESHOLD && DenseSubstrings.substrings[i].rightWindow>DenseSubstrings.substrings[j].maxWindow*threshold);
	}


	/**
	 * Discards new weak substrings that are contained in, or identical to, a periodic
	 * substring interval of any type. This is essentially the same as 
	 * $DenseSubstrings.discardPeriodicSubstrings()$.
	 *
	 * Remark: the procedure sets the $periodicSubstringInterval$ field of alignments, but
	 * it does not set their $impliedByPeriodicSubstring$ pointer. This is because the
	 * procedure is intended to be called at the very end of the pipeline, when no step
	 * depends on such pointers any more. If such pointers are needed, we could create
	 * artificial periodic substrings from the corresponding intervals, as done in 
	 * $PeriodicSubstrings.assignCompatibleAlignments()$.
	 *
	 * Remark: the procedure uses array $nMergedIntervalsPoints$ as temporary space.
	 *
	 * @param from $DenseSubstrings.substrings[from..]$ is assumed to contain all and only
	 * the new weak substrings, sorted by $startA$.
	 */
	private static final void cleanNewDenseSubstrings_periodic(int from) {
		final int N_CONSECUTIVE_POINTS = 2;
		final double QUANTILE = 0.75;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		boolean atLeastOneDiscarded, hasLongPeriodJ, found;
		int i, j, k;
		int startI, endI, startJ, endJ, periodJ, longestAlignmentJ, minInterval, length, minLength, firstJForNextI;
		int previousLastSubstring, newLastSubstring, threshold, longPeriodThreshold, lastPoint;
		int minIntervalLength, minLocalMaxDistance, nLocalMaximumLeaves;
		DenseSubstring tmpSubstring;
		
		for (i=from; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].periodicSubstringInterval=null;
		if (PeriodicSubstrings.lastInterval==-1) return;
		previousLastSubstring=DenseSubstrings.lastSubstring;
		

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("cleanNewDenseSubstrings_periodic> periodic intervals at the very beginning:");
	for (int x=0; x<=PeriodicSubstrings.lastInterval; x++) IO.printErr(PeriodicSubstrings.intervals[x].hashCode()+" :: "+PeriodicSubstrings.intervals[x]);
}
		
		// Estimating $longPeriodThreshold$ using dense substrings that are strictly
		// contained in a long-period interval.
		Factorize.computeMaximalAlignments(true,false,false,ReadA.lastSortedAlignment);
		longPeriodThreshold=DenseSubstrings.LONG_PERIOD_THRESHOLD;
		i=from; j=0; firstJForNextI=-1; lastPoint=-1;
		minInterval=-1; minLength=Math.POSITIVE_INFINITY;
		while (i<=DenseSubstrings.lastSubstring) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("cleanNewDenseSubstrings_periodic> dense="+DenseSubstrings.substrings[i].hashCode()+" periodic="+PeriodicSubstrings.intervals[j].hashCode());
			if (!DenseSubstrings.substrings[i].isWeak || j>PeriodicSubstrings.lastInterval || PeriodicSubstrings.intervals[j].firstPosition>=DenseSubstrings.substrings[i].endA) {
				if (minInterval!=-1) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("cleanNewDenseSubstrings_periodic> new weak dense="+DenseSubstrings.substrings[i].hashCode()+" contained in periodic="+PeriodicSubstrings.intervals[minInterval].hashCode()+"::"+PeriodicSubstrings.intervals[minInterval]+" => nMaximalAlignments="+DenseSubstrings.substrings[i].nMaximalAlignments);
					lastPoint++;
					nMergedIntervalsPoints[lastPoint].position=DenseSubstrings.substrings[i].nMaximalAlignments;
					nMergedIntervalsPoints[lastPoint].mass=1;
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				minInterval=-1; minLength=Math.POSITIVE_INFINITY;
				continue;
			}
			if (PeriodicSubstrings.intervals[j].lastPosition<=DenseSubstrings.substrings[i].startA) {
				j++;
				continue;
			}
			if (i<DenseSubstrings.lastSubstring && firstJForNextI==-1 && PeriodicSubstrings.intervals[j].lastPosition>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
			if ( PeriodicSubstrings.intervals[j].hasLongPeriod && 
			     ( Intervals.isApproximatelyContained_lowQuality(DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.id) &&
				   !Intervals.areApproximatelyIdentical_lowQuality(DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,PeriodicSubstrings.intervals[j].firstPosition,PeriodicSubstrings.intervals[j].lastPosition,ReadA.id)
				 )	 
			   ) {
				length=PeriodicSubstrings.intervals[j].lastPosition-PeriodicSubstrings.intervals[j].firstPosition;
				if (length<minLength) {
					minLength=length;
					minInterval=j;
				}
			}
			j++;
		}
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("cleanNewDenseSubstrings_periodic> array of counts:");
	for (int x=0; x<=lastPoint; x++) IO.printErr(nMergedIntervalsPoints[x]+"");
}
		if (lastPoint>0) {
			lastPoint=Points.sortAndCompact(nMergedIntervalsPoints,lastPoint);				
			if (lastPoint+1<=Points.FEW_POINTS) {
				i=Points.getRoughThreshold(nMergedIntervalsPoints,lastPoint,true,DenseSubstrings.LONG_PERIOD_THRESHOLD,false);
				if (i!=-1) {
					threshold=(int)nMergedIntervalsPoints[i].position;
					longPeriodThreshold=Math.min(longPeriodThreshold,threshold);
				}
			}
			else {
				minIntervalLength=Math.max( (int)Points.distanceQuantile(nMergedIntervalsPoints,lastPoint,N_CONSECUTIVE_POINTS,true,QUANTILE),
											(int)((nMergedIntervalsPoints[lastPoint].position-nMergedIntervalsPoints[0].position)/MIN_INTERVAL_LENGTH_FACTOR)
										  );
				minLocalMaxDistance=minIntervalLength<<1;
				if (!Points.areUniformlyDistributed(nMergedIntervalsPoints,0,lastPoint,true,(nMergedIntervalsPoints[lastPoint].position-nMergedIntervalsPoints[0].position)/Points.DEFAULT_NBINS)) {
					nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(nMergedIntervalsPoints,0,lastPoint,minIntervalLength,minLocalMaxDistance,true,-1,-1,-1,false,true);
					if (nLocalMaximumLeaves>0) {
						DensityEstimationTree.markRunsOfLocalMaximumLeaves(nMergedIntervalsPoints);
						threshold=(int)DensityEstimationTree.separateRun(0,nMergedIntervalsPoints,0,lastPoint,true);
						longPeriodThreshold=Math.min(longPeriodThreshold,threshold);
					}
				}
			}
		}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("cleanNewDenseSubstrings_periodic> longPeriodThreshold="+longPeriodThreshold);
	
		// Discarding substrings
		atLeastOneDiscarded=false;
		i=from; j=0; firstJForNextI=-1;
		startI=DenseSubstrings.substrings[i].startA;
		endI=DenseSubstrings.substrings[i].endA;
		startJ=PeriodicSubstrings.intervals[j].firstPosition;
		endJ=PeriodicSubstrings.intervals[j].lastPosition;
		hasLongPeriodJ=PeriodicSubstrings.intervals[j].hasLongPeriod;
		periodJ=PeriodicSubstrings.intervals[j].period;
		longestAlignmentJ=PeriodicSubstrings.intervals[j].longestAlignment;
		minInterval=-1; minLength=Math.POSITIVE_INFINITY;
		while (i<=DenseSubstrings.lastSubstring) {			
			if (j==PeriodicSubstrings.lastInterval+1 || startJ>=endI) {
				if (minInterval!=-1) {
					DenseSubstrings.substrings[i].periodicSubstringInterval=PeriodicSubstrings.intervals[minInterval];
					atLeastOneDiscarded=true;
				}
				i++; minInterval=-1; minLength=Math.POSITIVE_INFINITY;
				if (i<=DenseSubstrings.lastSubstring) {
					startI=DenseSubstrings.substrings[i].startA;
					endI=DenseSubstrings.substrings[i].endA;
				}
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					if (j<=PeriodicSubstrings.lastInterval) {
						startJ=PeriodicSubstrings.intervals[j].firstPosition;
						endJ=PeriodicSubstrings.intervals[j].lastPosition;
						hasLongPeriodJ=PeriodicSubstrings.intervals[j].hasLongPeriod;
						periodJ=PeriodicSubstrings.intervals[j].period;
						longestAlignmentJ=PeriodicSubstrings.intervals[j].longestAlignment;
					}
				}
				firstJForNextI=-1;
				continue;
			}
			if (endJ<=startI) {
				j++;
				if (j<=PeriodicSubstrings.lastInterval) {
					startJ=PeriodicSubstrings.intervals[j].firstPosition;
					endJ=PeriodicSubstrings.intervals[j].lastPosition;
					hasLongPeriodJ=PeriodicSubstrings.intervals[j].hasLongPeriod;
					periodJ=PeriodicSubstrings.intervals[j].period;
					longestAlignmentJ=PeriodicSubstrings.intervals[j].longestAlignment;
				}
				continue;
			}
			if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && endJ>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;			
			if ( Intervals.areApproximatelyIdentical_lowQuality(startI,endI,startJ,endJ,ReadA.id) || 
			     ( Intervals.isApproximatelyContained_lowQuality(startI,endI,startJ,endJ,ReadA.id) &&
				   ( !hasLongPeriodJ || 
					 ( (periodJ>0 && endI-startI+1>periodJ) || 
					   (periodJ<=0 && endI-startI+1>longestAlignmentJ) ||
					   DenseSubstrings.substrings[i].nMaximalAlignments<longPeriodThreshold
				     )
				   )
				 )
			   ) {
				length=endJ-startJ;
				if (length<minLength) {
					minLength=length;
					minInterval=j;
				}
			}
			j++;
			if (j<=PeriodicSubstrings.lastInterval) {
				startJ=PeriodicSubstrings.intervals[j].firstPosition;
				endJ=PeriodicSubstrings.intervals[j].lastPosition;
				hasLongPeriodJ=PeriodicSubstrings.intervals[j].hasLongPeriod;
				periodJ=PeriodicSubstrings.intervals[j].period;
				longestAlignmentJ=PeriodicSubstrings.intervals[j].longestAlignment;
			}
		}
		if (atLeastOneDiscarded) {
			j=from-1;
			for (i=from; i<=DenseSubstrings.lastSubstring; i++) {
				if (DenseSubstrings.substrings[i].periodicSubstringInterval!=null) continue;
				j++;
				if (j!=i) {
					tmpSubstring=DenseSubstrings.substrings[j];
					DenseSubstrings.substrings[j]=DenseSubstrings.substrings[i];
					DenseSubstrings.substrings[i]=tmpSubstring;
				}	
			}
			newLastSubstring=j;
			// Resetting pointers from long-period clones to the real long-period
			// intervals.
			for (i=from; i<=DenseSubstrings.lastSubstring; i++) {
				if (DenseSubstrings.substrings[i].periodicSubstringInterval==null || !DenseSubstrings.substrings[i].periodicSubstringInterval.hasLongPeriod) continue;
				j=Arrays.binarySearch(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1,DenseSubstrings.substrings[i].periodicSubstringInterval);
				if (j<0) {
					System.err.println("cleanNewDenseSubstrings_periodic> ERROR: long-period interval not found.");
					System.exit(1);
				}
				found=false;
				for (k=j; k>=0; k--) {
					if (PeriodicSubstrings.longPeriodIntervals[k].equals(DenseSubstrings.substrings[i].periodicSubstringInterval)) {
						DenseSubstrings.substrings[i].periodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[k];
						found=true;
						break;
					}
				}
				if (found) continue;
				for (k=j+1; k<=PeriodicSubstrings.lastLongPeriodInterval; k++) {
					if (PeriodicSubstrings.longPeriodIntervals[k].equals(DenseSubstrings.substrings[i].periodicSubstringInterval)) {
						DenseSubstrings.substrings[i].periodicSubstringInterval=PeriodicSubstrings.longPeriodIntervals[k];
						break;
					}
				}
			}
			DenseSubstrings.lastSubstring=newLastSubstring;
		}
	
		// Updating pointers from alignments
		if (DenseSubstrings.lastSubstring!=previousLastSubstring) {
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				tmpSubstring=ReadA.sortedAlignments[i].impliedByDenseSubstring;
				if (tmpSubstring!=null && tmpSubstring.periodicSubstringInterval!=null) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("cleanNewDenseSubstrings_periodic> (1) discarding alignment "+ReadA.sortedAlignments[i].toStringBoundaries());
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].periodicSubstringInterval=tmpSubstring.periodicSubstringInterval;
				}
				tmpSubstring=ReadA.sortedAlignments[i].inDenseSubstring;
				if (tmpSubstring!=null && tmpSubstring.periodicSubstringInterval!=null) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("cleanNewDenseSubstrings_periodic> (2) discarding alignment "+ReadA.sortedAlignments[i].toStringBoundaries());
					ReadA.sortedAlignments[i].inDenseSubstring=null;
					ReadA.sortedAlignments[i].periodicSubstringInterval=tmpSubstring.periodicSubstringInterval;
				}
			}
			for (i=from; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].periodicSubstringInterval=null;
		}
	}


	/**
	 * Discards alignment intervals that are approximately identical to a dense substring 
	 * produced by a concatenation, and reassigns their alignments to the substring.
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings[from..]$ to contain all 
	 * and only the new substrings resulting from a concatenation, sorted by $startA$.
	 * All such substrings are assumed to be of substring type.
	 */
	private static final void cleanNewDenseSubstrings_alignmentIntervals(int from) {
		boolean atLeastOneDiscarded;
		int i, j;
		int startI, endI, startJ, endJ, firstJForNextI;
		AlignmentInterval tmpInterval;
		DenseSubstring tmpSubstring;
		
		// Ensuring the necessary order in $intervals$
		if (lastInterval==-1) return;
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}

		// Discarding intervals
		for (i=0; i<=lastInterval; i++) intervals[i].denseSubstring=null;
		i=0; j=from; firstJForNextI=-1; atLeastOneDiscarded=false;
		startI=intervals[i].firstPosition;
		endI=intervals[i].lastPosition;
		startJ=DenseSubstrings.substrings[j].startA;
		endJ=DenseSubstrings.substrings[j].endA;
		while (i<=lastInterval) {
			if (j==DenseSubstrings.lastSubstring+1 || startJ>=endI) {
				i++;
				if (i<=lastInterval) {
					startI=intervals[i].firstPosition;
					endI=intervals[i].lastPosition;
				}
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					startJ=DenseSubstrings.substrings[j].startA;
					endJ=DenseSubstrings.substrings[j].endA;
				}
				firstJForNextI=-1;
				continue;
			}
			if (endJ<startI) {
				j++;
				if (j<=DenseSubstrings.lastSubstring) {
					startJ=DenseSubstrings.substrings[j].startA;
					endJ=DenseSubstrings.substrings[j].endA;
				}
				continue;
			}
			if (firstJForNextI==-1 && i<lastInterval && endJ>=intervals[i+1].firstPosition) firstJForNextI=j;
			if (Intervals.areApproximatelyIdentical(startI,endI,startJ,endJ)) {
				intervals[i].denseSubstring=DenseSubstrings.substrings[j];
				atLeastOneDiscarded=true;
			}
			j++;
			if (j<=DenseSubstrings.lastSubstring) {
				startJ=DenseSubstrings.substrings[j].startA;
				endJ=DenseSubstrings.substrings[j].endA;
			}
		}
		if (atLeastOneDiscarded) {
			j=-1;
			for (i=0; i<=lastInterval; i++) {
				if (intervals[i].denseSubstring!=null) continue;
				j++;
				if (j!=i) {
					tmpInterval=intervals[j];
					intervals[j]=intervals[i];
					intervals[i]=tmpInterval;
				}	
			}
			lastInterval=j;
		}
		
		// Updating pointers from alignments
		if (atLeastOneDiscarded) {
			for (i=0; i<=ReadA.lastSortedAlignment; i++) {
				tmpInterval=ReadA.sortedAlignments[i].mergedToInterval;
				if (tmpInterval!=null && tmpInterval.denseSubstring!=null) {
					ReadA.sortedAlignments[i].mergedToInterval=null;
					ReadA.sortedAlignments[i].impliedByDenseSubstring=null;
					ReadA.sortedAlignments[i].inDenseSubstring=tmpInterval.denseSubstring;
				}
			}
			for (i=0; i<=lastInterval; i++) intervals[i].denseSubstring=null;
		}
	}


	/**
	 * Sets the $inDenseSubstring$ field of every alignment that is contained inside a 
	 * newly-created dense substring of substring type, and that is not already assigned 
	 * to a dense substring or periodic substring interval. If the alignment is already
	 * assigned to an alignment interval, and if such alignment interval is approximately 
	 * identical to a new dense susbtring, the alignment is assigned to the new dense 
	 * substring.
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings[from..]$ to contain all 
	 * and only the newly-created dense substrings of substring type, sorted by $startA$.
	 */
	private static final void markAlignmentsInNewDenseSubstrings(int from) {
		boolean alignmentIntervalDiscarded;
		int i, j;
		int firstJForNextI, length, minLength, startA, endA;
		DenseSubstring tmpSubstring, minSubstring;
		AlignmentInterval aInterval;
		
		// Enforcing the necessary order in $ReadA.sortedAlignments$
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}
		
		// Marking
		alignmentIntervalDiscarded=false;
		for (i=0; i<=lastInterval; i++) intervals[i].discarded=false;
		i=0; j=from; firstJForNextI=-1; minSubstring=null; minLength=Math.POSITIVE_INFINITY;
		while (i<=ReadA.lastSortedAlignment) {
			if ( ReadA.sortedAlignments[i].periodicSubstringInterval!=null || 
				 ReadA.sortedAlignments[i].impliedByDenseSubstring!=null || 
				 ReadA.sortedAlignments[i].inDenseSubstring!=null
			   ) {
			   i++;
			   continue;
			}
			aInterval=ReadA.sortedAlignments[i].mergedToInterval;
			startA=ReadA.sortedAlignments[i].startA();
			endA=ReadA.sortedAlignments[i].endA();
			if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>endA) {
				if (minSubstring!=null) {
					ReadA.sortedAlignments[i].inDenseSubstring=minSubstring;
					ReadA.sortedAlignments[i].mergedToInterval=null;
	   				if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && ReadA.sortedAlignments[i].isRightMaximalB==1) {
						length=ReadA.sortedAlignments[i].getALength();
						if (length<minSubstring.minPrefixLength) minSubstring.minPrefixLength=length;
					}
					if (ReadA.sortedAlignments[i].isLeftMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityStart && startA>=minSubstring.startA && (minSubstring.maximalStartA>=0?(startA<minSubstring.maximalStartA):true)) minSubstring.maximalStartA=startA;
					if (ReadA.sortedAlignments[i].isRightMaximalB==1 && !ReadA.sortedAlignments[i].lowQualityEnd && endA<=minSubstring.endA && endA>minSubstring.maximalEndA) minSubstring.maximalEndA=endA;
					minSubstring=null;
					minLength=Math.POSITIVE_INFINITY;
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			tmpSubstring=DenseSubstrings.substrings[j];
			if (tmpSubstring.endA<startA) {
				j++;
				continue;
			}
			if (i<ReadA.lastSortedAlignment && firstJForNextI==-1 && tmpSubstring.endA>=ReadA.sortedAlignments[i+1].startA()) firstJForNextI=j;
			if (tmpSubstring.substringReplication) {
				if ( aInterval==null &&
			         ( Intervals.isApproximatelyContained(startA,endA,tmpSubstring.startA,tmpSubstring.endA) ||
				       Intervals.areApproximatelyIdentical(startA,endA,tmpSubstring.startA,tmpSubstring.endA)
				     )
				   ) {
					length=tmpSubstring.endA-tmpSubstring.startA;
					if (length<minLength) {
						minLength=length;
						minSubstring=tmpSubstring;
					}
				}
				if ( aInterval!=null &&
				     Intervals.areApproximatelyIdentical(aInterval.firstPosition,aInterval.lastPosition,tmpSubstring.startA,tmpSubstring.endA)
				   ) {
					alignmentIntervalDiscarded=true;
					aInterval.discarded=true;
					length=tmpSubstring.endA-tmpSubstring.startA;
					if (length<minLength) {
						minLength=length;
						minSubstring=tmpSubstring;
					}
				}
			}
			j++;
		}
		
		// Removing alignment intervals
		if (!alignmentIntervalDiscarded) return;
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].discarded) continue;
			j++;
			if (j!=i) {
				aInterval=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=aInterval;
			}
		}
		lastInterval=j;
	}


	/**
	 * Sets to null the $mergedToInterval$ pointer of an alignment, iff the readA interval
	 * of the alignment straddles the alignment interval pointed to, but does not have a 
	 * large enough intersection with it.
	 */
	private static final void discardAlignmentsAfterTrimming(int firstAlignment, int lastAlignment) {
		final double INTERSECTION_RATIO = 0.75;  // Arbitrary
		int i;
		int startA, endA, substringStartA, substringEndA, length;
		
		for (i=firstAlignment; i<=lastAlignment; i++) {
			startA=ReadA.sortedAlignments[i].startA();
			endA=ReadA.sortedAlignments[i].endA();
			length=endA-startA+1;
			if (ReadA.sortedAlignments[i].mergedToInterval!=null) {
				substringStartA=ReadA.sortedAlignments[i].mergedToInterval.firstPosition;
				substringEndA=ReadA.sortedAlignments[i].mergedToInterval.lastPosition;
				if (Intervals.straddles(startA,endA,substringStartA,substringEndA) && Intervals.intersectionLength(startA,endA,substringStartA,substringEndA)<length*INTERSECTION_RATIO) ReadA.sortedAlignments[i].mergedToInterval=null;
			}
		}
	}
	
	
	/**
	 * Just a wrapper of $DenseSubstrings.concatenateWeakSubstrings()$.
	 *
	 * Remark: concatenating weak dense substrings of substring type, and even merging a
	 * non-weak dense substring of substring type contained in a weak dense substring of
	 * substring type, makes sense at this stage, since the result of such merges is still
	 * marked as weak, and all intervals that would be produced by splitting the result of
	 * the merge have very likely already been detected.
	 *
	 * Remark: the procedure assumes all dense substrings to be sorted by $startA$.
	 *
	 * @return new substrings that result from the concatenation are stored in $[from..
	 * lastSubstring]$, sorted by startA, where $from$ is returned in output; range 
	 * $[0..from-1]$ is kept sorted by startA as well; the output is -1 if no 
	 * concatenation was performed.
	 */
	private static final int concatenateWeakSubstrings() {
		boolean found;
		int i, j;
		int last, firstJForNextI, out;
		DenseSubstring tmpSubstring;
		int[] tmpArray = DenseSubstrings.stack;
		
		// Concatenating
		concatenateWeakSubstrings_markStraddlingIntervals();
		last=DenseSubstrings.concatenateWeakSubstrings(ReadA.lastSortedAlignment,tmpArray);
		if (last==-1) return -1;
		concatenateWeakSubstrings_discardStraddlingIntervals();
		
		// Moving new substrings to the end of the list
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].discarded=false;  // Using $discarded$ just as a temporary flag.
		i=0; j=0; firstJForNextI=-1; found=false;
		while (i<=DenseSubstrings.lastSubstring) {
			if (j>last || tmpArray[j]>DenseSubstrings.substrings[i].startA) {
				DenseSubstrings.substrings[i].discarded=found;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				found=false;
				continue;
			}
			if (tmpArray[j+2]==-1 || tmpArray[j]<DenseSubstrings.substrings[i].startA) {
				j+=3;
				continue;
			}
			if (firstJForNextI==-1 && i<DenseSubstrings.lastSubstring && tmpArray[j]>=DenseSubstrings.substrings[i+1].startA) firstJForNextI=j;
			if ( (tmpArray[j]==DenseSubstrings.substrings[i].startA || tmpArray[j]==DenseSubstrings.substrings[i].maximalStartA) && 
			     (tmpArray[j+1]==DenseSubstrings.substrings[i].endA || tmpArray[j+1]==DenseSubstrings.substrings[i].maximalEndA)
			   ) found=true;
			j+=3;
		}
		j=-1;
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			if (DenseSubstrings.substrings[i].discarded) continue;
			j++;
			if (j!=i) {
				tmpSubstring=DenseSubstrings.substrings[j];
				DenseSubstrings.substrings[j]=DenseSubstrings.substrings[i];
				DenseSubstrings.substrings[i]=tmpSubstring;
			}
		}
		out=j+1;
		DenseSubstring.order=DenseSubstring.STARTA;
		if (DenseSubstrings.lastSubstring>out) Arrays.sort(DenseSubstrings.substrings,out,DenseSubstrings.lastSubstring+1);
		for (i=out; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].discarded=false;
		
		return out;
	}
	
	
	/**
	 * Sets field $tmpInt1=1$ for every alignment interval that straddles two dense 
	 * substrings X and Y, such that neither X nor Y is contained in the interval (X and Y
	 * can straddle each other). $tmpInt1$ is set to zero otherwise.
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ to be sorted by first 
	 * position.
	 */
	private static final void concatenateWeakSubstrings_markStraddlingIntervals() {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i, j, k;
		int firstJForNextI, intervalStart, intervalEnd;
		int substring1Start, substring1End, substring2Start, substring2End;
		
		// Ensuring the necessary order in $intervals$
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("concatenateWeakSubstrings_markStraddlingIntervals> alignment intervals at the very beginning:");
	for (int x=0; x<=lastInterval; x++) System.err.println(intervals[x]);
	System.err.println("concatenateWeakSubstrings_markStraddlingIntervals> dense substrings at the very beginning:");
	for (int x=0; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
}	
		
		
		for (i=0; i<=lastInterval; i++) intervals[i].tmpInt1=0;
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastInterval) {
			intervalStart=intervals[i].firstPosition;
			if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>=intervalStart-IDENTITY_THRESHOLD) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			substring1End=DenseSubstrings.substrings[j].endA;
			if (firstJForNextI==-1 && i<lastInterval && substring1End>=intervals[i+1].firstPosition) firstJForNextI=j;
			intervalEnd=intervals[i].lastPosition;
			if (substring1End<=intervalStart+IDENTITY_THRESHOLD || substring1End>=intervalEnd-IDENTITY_THRESHOLD) {
				j++;
				continue;
			}
			for (k=j+1; k<=DenseSubstrings.lastSubstring; k++) {
				substring2Start=DenseSubstrings.substrings[k].startA;
				if (substring2Start>=intervalEnd-IDENTITY_THRESHOLD) break;
				if (substring2Start<=intervalStart+IDENTITY_THRESHOLD) continue;
				substring2End=DenseSubstrings.substrings[k].endA;
				if (substring2End<=intervalEnd+IDENTITY_THRESHOLD) continue;
				intervals[i].tmpInt1=1;
if (IO.SHOW_STD_ERR_PRIME) System.err.println("concatenateWeakSubstrings_markStraddlingIntervals> marked interval "+intervals[i]);				
				break;
			}
			j++;
		}
	}
	
	
	/**
	 * Discards alignment intervals that: (1) straddled two dense substrings before 
	 * $concatenateWeakSubstrings()$; (2) do not straddle two dense subsbtrings after
	 * $concatenateWeakSubstrings()$; (3) have few alignments. The alignments that pointed
	 * to such intervals are not reassigned to other intervals or dense substrings.
	 *
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ to be sorted by first 
	 * position.
	 */
	private static final void concatenateWeakSubstrings_discardStraddlingIntervals() {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_MAXIMAL_ALIGNMENTS = MIN_INTERVALS_IN_MAXIMAL;  // Arbitrary
		boolean straddling;
		int i, j, k;
		int firstJForNextI, intervalStart, intervalEnd;
		int substring1Start, substring1End, substring2Start, substring2End;
		AlignmentInterval tmpInterval;
		
		// Ensuring the necessary order in $intervals$
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
		}
		
		// Detecting intervals to be discarded
		for (i=0; i<=lastInterval; i++) intervals[i].discarded=false;
		i=0; j=0; firstJForNextI=-1; straddling=false;
		while (i<=lastInterval) {
			intervalStart=intervals[i].firstPosition;
			if (intervals[i].tmpInt1!=1 || intervals[i].nMaximalAlignmentsLeft>=MIN_MAXIMAL_ALIGNMENTS || intervals[i].nMaximalAlignmentsRight>=MIN_MAXIMAL_ALIGNMENTS) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				straddling=false;
				continue;
			}
			if (j>DenseSubstrings.lastSubstring || DenseSubstrings.substrings[j].startA>=intervalStart-IDENTITY_THRESHOLD) {
				intervals[i].discarded=!straddling;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				straddling=false;
				continue;
			}
			substring1End=DenseSubstrings.substrings[j].endA;
			if (firstJForNextI==-1 && i<lastInterval && substring1End>=intervals[i+1].firstPosition) firstJForNextI=j;
			intervalEnd=intervals[i].lastPosition;
			if (substring1End<=intervalStart+IDENTITY_THRESHOLD || substring1End>=intervalEnd-IDENTITY_THRESHOLD) {
				j++;
				continue;
			}
			for (k=j+1; k<=DenseSubstrings.lastSubstring; k++) {
				substring2Start=DenseSubstrings.substrings[k].startA;
				if (substring2Start>=intervalEnd-IDENTITY_THRESHOLD) break;
				if (substring2Start<=intervalStart+IDENTITY_THRESHOLD) continue;
				substring2End=DenseSubstrings.substrings[k].endA;
				if (substring2End<=intervalEnd+IDENTITY_THRESHOLD) continue;
				straddling=true;
				break;
			}
			j++;
		}
		
		// Compacting intervals
		j=-1;
		for (i=0; i<=lastInterval; i++) {
			if (intervals[i].discarded) continue;
			j++;
			if (j!=i) {
				tmpInterval=intervals[j];
				intervals[j]=intervals[i];
				intervals[i]=tmpInterval;
			}
		}
		lastInterval=j;
		
		// Updating pointers from alignments
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			tmpInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (tmpInterval==null) continue;
			if (tmpInterval.discarded) ReadA.sortedAlignments[i].mergedToInterval=null;
		}
	}
	
	
	
	

}