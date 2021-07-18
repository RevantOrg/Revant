package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Histograms;


public class ReadA {
	public static int id;  // Starts from zero
	public static int readLength;
	public static int firstAlignment;
	public static int lastAlignment;
	public static int nAlignments;

	/**
	 * Parameters of the pipeline
	 */
	public static double maxDiffsAlignmentFraction = 0.5;  // Used by procedure $estimateMaxDiffs$

	/**
	 * Parameters of the pipeline, estimated from the data.
	 */
	public static double maxDiffs;  // Maximum value of $diffs/(\ell_a+\ell_b)$ for an alignment to involve two instances of the same repeat, rather than two variants of the repeat. $\ell_i$ is the length of the interval of an alignment in read $i$.
	public static int nVariants;  // Lower bound on the number of variants
	public static double maxDiffsRead;  // A version of $maxDiffs$ computed on the whole readA

	/**
	 * Data structures for histograms
	 */
	public static double[] coverage;

	/**
	 * Data structures for sorting alignments
	 */
	public static Alignment[] sortedAlignments;
	public static int lastSortedAlignment;

	/**
	 * Data structures for alignment qualities with a given substring of $readA$
	 */
	private static Point[] alignmentQualities;
	public static int lastQuality;

	/**
	 * Temporary space
	 */
	public static int[] deltaTmp;
	private static int[] counts;  // Temporary space used by $getDeltaHistogram$
	public static Event[] deltaPoints;
	public static int lastDeltaPoint;


	public static final void allocateMemory(int maxReadLength, int maxAlignments) {
		int i;

		coverage = new double[maxReadLength];
		counts = new int[maxReadLength];
		sortedAlignments = new Alignment[maxAlignments];
		for (i=0; i<maxAlignments; i++) sortedAlignments[i] = new Alignment();
		alignmentQualities = new Point[maxAlignments];
		for (i=0; i<alignmentQualities.length; i++) alignmentQualities[i] = new Point();
		deltaTmp = new int[2];
		deltaPoints = new Event[maxAlignments<<1];
		for (i=0; i<deltaPoints.length; i++) deltaPoints[i] = new Event();
	}


	/**
	 * Sets $id$, $readLength$, $firstAlignment$, $lastAlignment$, $sortedAlignments$.
	 *
	 * @param fa first alignment.
	 */
	public static final void initialize(int fa, boolean removeContainedAlignments, boolean estimateMaxDiffs) {
		int i;
		double backup;

		firstAlignment=fa;
		id=Alignments.alignments[firstAlignment][0]-1;
		readLength=Reads.getReadLength(id);
		lastAlignment=firstAlignment+1;
		while (lastAlignment<Alignments.nAlignments && Alignments.alignments[lastAlignment][0]==Alignments.alignments[firstAlignment][0]) lastAlignment++;
		lastAlignment--;
		nAlignments=lastAlignment-firstAlignment+1;
		lastSortedAlignment=-1; lastQuality=-1;
		for (i=firstAlignment; i<=lastAlignment; i++) {
			lastSortedAlignment++;
			sortedAlignments[lastSortedAlignment].initialize(i);
		}
		Alignment.order=Alignment.READB_ORIENTATION_STARTA;

		// Setting maximality
		if (removeContainedAlignments) removeContainedAlignments();
		for (i=0; i<=lastSortedAlignment; i++) sortedAlignments[i].setMaximality();

		// Estimating $maxDiffsRead$
		if (estimateMaxDiffs) {
			backup=maxDiffsAlignmentFraction;
			maxDiffsAlignmentFraction=0.0;
			estimateMaxDiffs(0,readLength-1,0,lastSortedAlignment,false);
			maxDiffsAlignmentFraction=backup;
			maxDiffsRead=maxDiffs;
		}
		else maxDiffsRead=-1;
	}


	/**
	 * Removes some alignments in $sortedAlignments$ that are contained in other
	 * alignments. We assume that a contained alignment provides more information than all
	 * its containing alignments, iff the minimum of the difference between the average
	 * $diffs$ of a containing alignment and the average $diffs$ of the alignment is
	 * higher than a positive threshold. Such a threshold is computed by estimating the
	 * density of all such minimum differences in readA, and by keeping the rightmost 
	 * positive local maximum.
	 *
	 * Alternatively, we could use ratios rather than differences.
	 */
	private static final void removeContainedAlignments() {
		final int N_CONSECUTIVE_POINTS = 2;
		final double QUANTILE = 0.75;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		final double MAX_DELTA_DIFFERENCE = Alignments.MAX_ALIGNMENT_ERROR/5;
		int i, j;
		int readBI, readBJ, orientationI, orientationJ, startAI, endAI, startBI, endBI, startAJ, endAJ, startBJ, endBJ;
		int nLocalMaximumLeaves, tmp;
		double diffsI, diffsJ, delta, deltaThreshold;
		double minIntervalLength, minLocalMaxDistance;
		Alignment tmpAlignment;

		// Ensuring the necessary order in $sortedAlignments$
		if (Alignment.order!=Alignment.READB_ORIENTATION_STARTA_ENDA) {
			Alignment.order=Alignment.READB_ORIENTATION_STARTA_ENDA;
			if (lastSortedAlignment>0) Arrays.sort(sortedAlignments,0,lastSortedAlignment+1);
		}

		// Marking contained alignments and computing $minDeltaDiff$
		for (i=0; i<=lastSortedAlignment; i++) {
			sortedAlignments[i].isContained=false;
			sortedAlignments[i].minDeltaDiff=Math.POSITIVE_INFINITY;
		}
		i=0; j=i+1;
		readBI=Alignments.alignments[sortedAlignments[i].id][1];
		orientationI=Alignments.alignments[sortedAlignments[i].id][2];
		startAI=Alignments.alignments[sortedAlignments[i].id][3];
		endAI=Alignments.alignments[sortedAlignments[i].id][4];
		startBI=Alignments.alignments[sortedAlignments[i].id][5];
		endBI=Alignments.alignments[sortedAlignments[i].id][6];
		diffsI=Alignments.getAvgDiffs(i);
		while (i<=lastSortedAlignment) {
			if (j>lastSortedAlignment) {
				i++;
				if (i>lastSortedAlignment) break;
				readBI=Alignments.alignments[sortedAlignments[i].id][1];
				orientationI=Alignments.alignments[sortedAlignments[i].id][2];
				startAI=Alignments.alignments[sortedAlignments[i].id][3];
				endAI=Alignments.alignments[sortedAlignments[i].id][4];
				startBI=Alignments.alignments[sortedAlignments[i].id][5];
				endBI=Alignments.alignments[sortedAlignments[i].id][6];
				diffsI=Alignments.getAvgDiffs(i);
				j=i+1;
				continue;
			}
			readBJ=Alignments.alignments[sortedAlignments[j].id][1];
			orientationJ=Alignments.alignments[sortedAlignments[j].id][2];
			startAJ=Alignments.alignments[sortedAlignments[j].id][3];
			if (readBJ!=readBI || orientationJ!=orientationI || startAJ>endAI) {
				i++;
				if (i>lastSortedAlignment) break;
				readBI=Alignments.alignments[sortedAlignments[i].id][1];
				orientationI=Alignments.alignments[sortedAlignments[i].id][2];
				startAI=Alignments.alignments[sortedAlignments[i].id][3];
				endAI=Alignments.alignments[sortedAlignments[i].id][4];
				startBI=Alignments.alignments[sortedAlignments[i].id][5];
				endBI=Alignments.alignments[sortedAlignments[i].id][6];
				diffsI=Alignments.getAvgDiffs(i);
				j=i+1;
				continue;
			}
			endAJ=Alignments.alignments[sortedAlignments[j].id][4];
			if (endAJ>endAI) {
				j++;
				continue;
			}
			startBJ=Alignments.alignments[sortedAlignments[j].id][5];
			endBJ=Alignments.alignments[sortedAlignments[j].id][6];
			if (startBJ>=startBI && endBJ<=endBI) {
				sortedAlignments[j].isContained=true;
				diffsJ=Alignments.getAvgDiffs(sortedAlignments[j].id);
				delta=diffsI-diffsJ;
				if (delta>0 && delta<sortedAlignments[j].minDeltaDiff) sortedAlignments[j].minDeltaDiff=delta;
			}
			j++;
		}
		
		// Fitting density estimation tree on $minDeltaDiff$ values
		lastQuality=-1;
		for (i=0; i<=lastSortedAlignment; i++) {
			if (!sortedAlignments[i].isContained || sortedAlignments[i].minDeltaDiff==Math.POSITIVE_INFINITY) continue;
			lastQuality++;
			alignmentQualities[lastQuality].position=sortedAlignments[i].minDeltaDiff;
			alignmentQualities[lastQuality].mass=1;
		}
		if (lastQuality==-1) deltaThreshold=Math.POSITIVE_INFINITY;
		else {
			lastQuality=Points.sortAndCompact(alignmentQualities,lastQuality);

if (IO.SHOW_STD_OUT) IO.printOut("MINDELATDIFF OF CONTAINED ALIGNMENTS, READ "+id);
if (IO.SHOW_STD_OUT) { for (int x=0; x<=lastQuality; x++) IO.printOut(alignmentQualities[x]); }

			if (lastQuality==0) deltaThreshold=Math.POSITIVE_INFINITY;
			else if (lastQuality+1<=Points.FEW_POINTS) {
				if (alignmentQualities[lastQuality].position-alignmentQualities[0].position<=MAX_DELTA_DIFFERENCE) deltaThreshold=Math.POSITIVE_INFINITY;
				else {
					tmp=Points.getRoughThreshold(alignmentQualities,lastQuality,false,0,true);
					deltaThreshold=tmp==-1?Math.POSITIVE_INFINITY:alignmentQualities[tmp].position;
				}
			}
			else {
				minIntervalLength=Math.max( Points.distanceQuantile(alignmentQualities,lastQuality,N_CONSECUTIVE_POINTS,false,QUANTILE),
											(alignmentQualities[lastQuality].position-alignmentQualities[0].position)/MIN_INTERVAL_LENGTH_FACTOR
										  );
				minLocalMaxDistance=minIntervalLength*2.0;
				if (Points.areUniformlyDistributed(alignmentQualities,0,lastQuality,false,(alignmentQualities[lastQuality].position-alignmentQualities[0].position)/Points.DEFAULT_NBINS)) deltaThreshold=Math.POSITIVE_INFINITY;
				else {
					nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(alignmentQualities,0,lastQuality,minIntervalLength,minLocalMaxDistance,false,-1,-1,-1,false,true);
					if (nLocalMaximumLeaves==0) deltaThreshold=Math.POSITIVE_INFINITY;
					else {
						DensityEstimationTree.markRunsOfLocalMaximumLeaves(alignmentQualities);
						deltaThreshold=DensityEstimationTree.separateRun(DensityEstimationTree.nRuns==1?-1:DensityEstimationTree.nRuns-2,alignmentQualities,0,lastQuality,true);
					}
				}
if (IO.SHOW_STD_OUT) IO.printOut("DELTA THRESHOLD = "+deltaThreshold+" READ "+id);
			}
		}

		// Removing contained alignments whose $minDeltaDiff$ is too small
		j=0; i=-1;
		while (j<=lastSortedAlignment) {
			if (!sortedAlignments[j].isContained || sortedAlignments[j].minDeltaDiff>=deltaThreshold) {
				i++;
				if (j!=i) {
					tmpAlignment=sortedAlignments[i];
					sortedAlignments[i]=sortedAlignments[j];
					sortedAlignments[j]=tmpAlignment;
				}
			}
else if (IO.SHOW_STD_ERR) IO.printErr("removed contained alignment ["+Alignments.alignments[sortedAlignments[j].id][3]+".."+Alignments.alignments[sortedAlignments[j].id][4]+"]");
			j++;
		}
		lastSortedAlignment=i;
	}


	public static final void print() {
		if (IO.SHOW_STD_ERR) IO.printErr("Read "+id+". Length: "+readLength+". Alignments: ["+firstAlignment+".."+lastAlignment+"].");
// 		if (IO.SHOW_STD_ERR) IO.printErr("Uncorrected coverage histogram (averages of consecutive "+IO.quantum+"bp bins):");
// 		Histograms.printHistogram(coverage,0,readLength-1,IO.quantum,0);
	}


	/**
	 * Stores in array $coverage$ the number of alignments that cover each position of
	 * the read.
	 *
	 * @param correct if TRUE, we do not count alignments for which another alignment
	 * exists that involves approximately the same substrings of $readA$ and of $readB$,
	 * but in opposite orientation. The resulting coverage of a position $i$ is an
	 * approximation of the number of distinct substrings of other reads that cover $i$.
	 *
	 * Remark: it is left to other procedures downstream to refine this rought estimate.
	 */
	public static final void getCoverageHistogram(boolean correct) {
		boolean found;
		int i, j;
		int readB, previousReadB, readBLength, startA, startB, endA, endB, firstInReverse;
		
		Math.set(coverage,readLength-1,0.0);
		previousReadB=-1; firstInReverse=-1;
		for (i=firstAlignment; i<=lastAlignment; i++) {
			readB=Alignments.alignments[i][1]-1;
			startA=Alignments.alignments[i][3];
			endA=Alignments.alignments[i][4];
			found=false;
			if (correct && Alignments.alignments[i][2]==1) {
				if (readB!=previousReadB) {
					previousReadB=readB;
					firstInReverse=-1;
					j=i+1;
					while (j<=lastAlignment) {
						if (Alignments.alignments[j][1]!=readB) break;
						if (Alignments.alignments[j][2]!=Alignments.alignments[i][2]) {
							firstInReverse=j;
							break;
						}
						j++;
					}
				}
				if (firstInReverse!=-1) {
					readBLength=Reads.getReadLength(readB);
					startB=Alignments.alignments[i][5];
					endB=Alignments.alignments[i][6];
					for (j=firstInReverse; j<=lastAlignment; j++) {
						if (Alignments.alignments[j][1]-1!=readB) break;
						if (Intervals.jaccardSimilarity(startA,endA,Alignments.alignments[j][3],Alignments.alignments[j][4])<Intervals.jaccardThreshold) continue;
						if (Intervals.jaccardSimilarity(startB,endB,Alignments.alignments[j][5],Alignments.alignments[j][6])<Intervals.jaccardThreshold) continue;
						found=true;
						break;
					}
				}
			}
			if (!found) {
				for (j=startA; j<=endA; j++) coverage[j]+=1.0;
			}
		}
	}


	public static final double getAverageCoverage(int start, int end) {
		return Histograms.getAverage(coverage,start,end);
	}


	/**
	 * Wraps $Histograms.getFraction()$ to count coverages that are at least 
	 * $minCoverage$.
	 */
	public static final double getFractionAtCoverage(int start, int end, double minCoverage) {
		return Histograms.getFraction(coverage,start,end,minCoverage,true);
	}


	/**
	 * Wraps $Histograms.hasSubstringWithAverage()$ to detect windows with average 
	 * coverage less than $minCoverage$.
	 */
	public static final boolean hasSubstringWithLowCoverage(int start, int end, double minCoverage, int window) {		
		return Histograms.hasSubstringWithAverage(coverage,start,end,minCoverage,false,window,true,-1);
	}
	

	/**
	 * Estimates $maxDiffs$ and $nVariants$ by fitting a density estimation tree on
	 * the qualities of all alignments that intersect at least $maxDiffsAlignmentFraction$
	 * fraction of $readA[start..end]$. $maxDiffs$ is set to a separation point between
	 * the first and the second local maximum of the density. $nVariants$ (a lower
	 * bound on the number of variants) is set to the number of local-maximum runs minus
	 * one.
	 *
	 * Remark: the procedure uses just $sortedAlignments[firstAlignment..lastAlignment]$,
	 * and if $sortedByStartA=true$ it assumes that such interval is sorted by startA.
	 *
	 * @return if $sortedByStartA=true$, the first alignment from $firstAlignment$ whose
	 * interval in $readA$ ends at least at $end$. If no such alignment exists, the
	 * procedure returns $firstAlignment$.
	 */
	public static final int estimateMaxDiffs(int start, int end, int firstAlignment, int lastAlignment, boolean sortedByStartA) {
		final int N_CONSECUTIVE_POINTS = 2;
		final int MAX_NPEAKS = 5;  // Maximum number of peaks to expect inside the range of alignment qualities
		final int MIN_INTERSECTION_LENGTH = Math.min((int)(maxDiffsAlignmentFraction*(end-start+1)),Alignments.minAlignmentLength);
		final double MAX_DIFFERENCE = Alignments.MAX_ALIGNMENT_ERROR/5;
		int i, j, out;
		int nLocalMaximumLeaves, startA, endA;
		double minIntervalLength, minLocalMaxDistance;

if (IO.SHOW_STD_ERR_PRIME) IO.printOut("estimateMaxDiffs> ["+start+".."+end+"] firstAlignment="+firstAlignment+" lastAlignment="+lastAlignment);

		// Collecting alignment qualities
		lastQuality=-1;
		out=-1;
		for (i=firstAlignment; i<=lastAlignment; i++) {
			startA=Alignments.alignments[sortedAlignments[i].id][3];
			if (sortedByStartA && startA>end) break;
			endA=Alignments.alignments[sortedAlignments[i].id][4];
			if (sortedByStartA && endA>=end && out==-1) out=i;
			if (Intervals.intersectionLength(startA,endA,start,end)<MIN_INTERSECTION_LENGTH) continue;
			lastQuality++;
			alignmentQualities[lastQuality].position=Alignments.getAvgDiffs(sortedAlignments[i].id);
			alignmentQualities[lastQuality].mass=1;
		}
		if (sortedByStartA && out==-1) out=firstAlignment;
		lastQuality=Points.sortAndCompact(alignmentQualities,lastQuality);
		
	
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("estimateMaxDiffs> qualities: (lastQuality="+lastQuality+")");
	for (i=0; i<=lastQuality; i++) IO.printErr(alignmentQualities[i].position+","+alignmentQualities[i].mass);
}
		

		// Fitting density estimation tree
		if (lastQuality==-1) {
			nVariants=0;
			maxDiffs=0.0;
		}
		else if (lastQuality==0) {
			nVariants=0;
			maxDiffs=alignmentQualities[0].position;
		}
		else if (lastQuality+1<=Points.FEW_POINTS) {
			if (alignmentQualities[lastQuality].position-alignmentQualities[0].position<=MAX_DIFFERENCE) {
				nVariants=0;
				maxDiffs=alignmentQualities[lastQuality].position;
			}
			else {
				i=Points.getRoughThreshold(alignmentQualities,lastQuality,true,0,true);
				if (i==-1) {
					nVariants=0;
					maxDiffs=alignmentQualities[lastQuality].position;
				}
				else {
					if (i>0) i--;
					j=Points.binarySearch(alignmentQualities,lastQuality,Alignments.MIN_ALIGNMENT_ERROR);
					if (j<0) j=-j-1;
					if (i>=j-1) {
						// The first peak intersects the typical error region: setting $maxDiffs$ to the maximum.
						nVariants=0;
						maxDiffs=alignmentQualities[lastQuality].position;
					}
					else {
						nVariants=j-i-1;
						maxDiffs=alignmentQualities[i].position;
					}
				}
			}
		}
		else {
			minIntervalLength=(alignmentQualities[lastQuality].position-alignmentQualities[0].position)/MAX_NPEAKS;
			minLocalMaxDistance=minIntervalLength;
			if (Points.areUniformlyDistributed(alignmentQualities,0,lastQuality,false,(alignmentQualities[lastQuality].position-alignmentQualities[0].position)/Points.DEFAULT_NBINS)) {
				// No local maximum: setting $maxDiffs$ to the maximum.
				nVariants=0;
				maxDiffs=alignmentQualities[lastQuality].position;
			}
			else {
				nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(alignmentQualities,0,lastQuality,minIntervalLength,minLocalMaxDistance,false,-1,-1,-1,false,true);
				if (nLocalMaximumLeaves==0) {
					// No local-maximum leaf: setting $maxDiffs$ to the maximum.
					nVariants=0;
					maxDiffs=alignmentQualities[lastQuality].position;
				}
				else {
					DensityEstimationTree.markRunsOfLocalMaximumLeaves(alignmentQualities);
					if (alignmentQualities[DensityEstimationTree.leaves[0].lastPoint].position>=Alignments.MIN_ALIGNMENT_ERROR) {
						// The first peak intersects the typical error region: setting $maxDiffs$ to the maximum.
						nVariants=0;
						maxDiffs=alignmentQualities[lastQuality].position;
					}
					else {
	
	
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Leaves of the DET: (minIntervalLength="+minIntervalLength+" minLocalMaxDistance="+minLocalMaxDistance+")");
	for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr(DensityEstimationTree.leaves[x]+"");
}						
						
						
						nVariants=DensityEstimationTree.nRuns-1;
						maxDiffs=DensityEstimationTree.separateRun(0,alignmentQualities,0,lastQuality,true);
					}
				}
			}
		}

if (IO.SHOW_STD_ERR_PRIME) IO.printOut("diffs ["+start+".."+end+"]:  (maxDiffs="+maxDiffs+", nVariants="+nVariants+")");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=lastQuality; x++) IO.printOut(alignmentQualities[x]); }

		return out;
	}


	/**
	 * Computes, for every position $i$ of $readA[interiorStart..interiorEnd]$, the 
	 * fraction of times $i$ is close to the first position of a left-maximal alignment or
	 * to the last position of a right-maximal alignment that covers $i$. Such values are 
	 * stored in array $delta[0..interiorEnd-interiorStart]$. If $start$ is smaller than 
	 * $interiorStart$, a left-maximal alignment with starting position in $[start..
	 * interiorStart-1]$ is not used to compute $delta$, but its left-maximal event is 
	 * sent to $Events.events$, i.e. it is used to estimate the first position of the 
	 * substring rather than to split the substring. A symmetrical argument holds if 
	 * $end>interiorEnd$. This behavior can be switched off by setting 
	 * $start=interiorStart$ and $end=interiorEnd$.
	 *
	 * Remark: the procedure performs a kernel estimation of the density of maximal 
	 * events, using a bandwidth derived from the points themselves (see procedure 
	 * $Points.estimateDensity$). Then, it divides such estimate by the number of times 
	 * each position is covered by an alignment.
	 * 
	 * Remark: the procedure actually uses B-maximality, rather than maximality, to
	 * compute function $\delta$.
	 *
	 * Remark: the procedure does not use alignments that have been marked as implied
	 * by a periodic (respectively, dense) substring. Specifically, the procedure assumes
	 * that $sortedAlignments[firstAlignment..lastAlignment]$ contains all and only such
	 * non-implied alignments of readA, sorted by startA. It does not alter such order.
	 *
	 * Remark: the procedure uses also alignments that are significantly outside
	 * $[start..end]$.
	 *
	 * Remark: the procedure takes the conservative approach of not using any alignment
	 * that is the prefix/suffix of a periodic substring. Since such alignments are not
	 * implied, they are likely noise caused by the fact that the aligner does not report 
	 * all alignments in periodic substrings. All such noisy alignments likely form the 
	 * characteristic delta histogram of a periodic substring. The peaks of such histogram
	 * are likely higher than those of e.g. a simple repeat inside the dense substring 
	 * that straddles the periodic substring (i.e. a signal we want to detect), since the 
	 * set of all such noisy alignments is potentially very large (because it depends both
	 * on the number of occurrences of the periodic substring in other reads, and on their
	 * lengths), inducing a large denominator. This conservative approach makes it 
	 * impossible to detect simple repeats that are the prefix/suffix of a periodic 
	 * substring.
	 *
	 * Remark: assume that $readA[start..end]$ is a periodic substring. An alignment whose
	 * intersection with readA is approximately the entire $[start..end]$ interval, but
	 * that is not implied by a periodic substring, suggests the existence of substrings
	 * of $readA[start..end]$ that are distinct repeats (or it suggests that $end-start+1
	 * \approx Alignments.minAlignmentLength$).
	 *
	 * Remark: we could compute a fast approximation of the delta function by using a
	 * random subset of all non-implied alignments, or just the $k$ alignments with
	 * largest alignment score, or with largest average quality in the aligned region of
	 * readB. However, the time bottleneck is very likely the construction of a regression
	 * tree on the histogram, rather than the construction of the histogram itself.
	 *
	 * @param left use just left- (TRUE) or right- (FALSE) maximal alignments;
	 * @param next a position >= end;
	 * @param strict TRUE: if an alignment is significantly inside (respectively, outside)
	 * $[start..end]$, and if its normalized number of differences is bigger than
	 * $maxDiffs$ (respectively, of $maxDiffsRead$), it is not used (the procedure assumes
	 * that $maxDiffs$ has already been computed). The reason for this is that such 
	 * alignments might come from variants of the substring rather than from its copies, 
	 * thus they might not provide reliable information for splitting. In practice this is
	 * rarely useful, and we keep this functionality just for historical reasons.
	 * @return $deltaTmp[0]$: number of alignments used to compute $delta$; $deltaTmp[1]$:
	 * the first alignment from $firstAlignment$ that either covers $next$, or starts to
	 * the right of $end$.
	 */
	public static final void getDeltaHistogram(int interiorStart, int interiorEnd, int start, int end, int firstAlignment, int lastAlignment, int next, boolean left, double[] delta, boolean strict) {
		final int MIN_BANDWIDTH = (IO.quantum<<1)/3;
		final int MAX_BANDWIDTH = IO.quantum;
		boolean isLeftMaximal, isRightMaximal, isLeftMaximalB, isRightMaximalB, startAOutside, endAOutside;
		int i, j, from, to, startA, endA, nUsedAlignments, firstAlignmentForNext;
		int bandwidth;
		double alignmentDiffs;

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getDeltaHistogram> firstAlignment="+firstAlignment+" lastAlignment="+lastAlignment+" LEFT="+left);

		// Filling $deltaPoints$
		Math.set(counts,counts.length-1,0);
		nUsedAlignments=0; firstAlignmentForNext=-1; lastDeltaPoint=-1;	
		for (i=firstAlignment; i<=lastAlignment; i++) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getDeltaHistogram> considering alignment "+i+"=["+Alignments.alignments[sortedAlignments[i].id][3]+".."+Alignments.alignments[sortedAlignments[i].id][4]+"], end="+end);
			endA=Alignments.alignments[sortedAlignments[i].id][4];
			if (firstAlignmentForNext==-1 && endA>=next) firstAlignmentForNext=i;
			if (sortedAlignments[i].isPrefixOfPeriodicSubstring || sortedAlignments[i].isSuffixOfPeriodicSubstring || endA<start) continue;
			startA=Alignments.alignments[sortedAlignments[i].id][3];
			if (startA>end) break;
			isLeftMaximal=sortedAlignments[i].isLeftMaximal==1;
			isRightMaximal=sortedAlignments[i].isRightMaximal==1;
			isLeftMaximalB=sortedAlignments[i].isLeftMaximalB==1;
			isRightMaximalB=sortedAlignments[i].isRightMaximalB==1;
			if ((left&&(!isLeftMaximalB||startA<start)) || (!left&&(!isRightMaximalB||endA>end))) continue;
			if (strict) {
				alignmentDiffs=Alignments.getAvgDiffs(sortedAlignments[i].id);
				if ( (Intervals.isApproximatelyContained(startA,endA,start,end) && alignmentDiffs>maxDiffs) ||
				     (!Intervals.isApproximatelyContained(startA,endA,start,end) && alignmentDiffs>maxDiffsRead) 
				   ) continue;
			}
			if (left && isLeftMaximal && startA>=start && startA<interiorStart) {
				Events.lastEvent++;
				Events.ensureSpace_events(Events.lastEvent+1);
				Events.events[Events.lastEvent].clear();
				Events.events[Events.lastEvent].position=startA;
				Events.events[Events.lastEvent].nOpen=1;
				sortedAlignments[i].startAdded=true;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getDeltaHistogram> added to Events "+startA);
			}
			if (!left && isRightMaximal && endA>interiorEnd && endA<=end) {
				Events.lastEvent++;
				Events.ensureSpace_events(Events.lastEvent+1);
				Events.events[Events.lastEvent].clear();
				Events.events[Events.lastEvent].position=endA;
				Events.events[Events.lastEvent].nClosed=1;
				sortedAlignments[i].endAdded=true;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getDeltaHistogram> added to Events "+endA);
			}
			startAOutside=(startA<interiorStart)||(startA>interiorEnd);
			endAOutside=(endA<interiorStart)||(endA>interiorEnd);
			if ((left&&(!isLeftMaximalB||startAOutside)) || (!left&&(!isRightMaximalB||endAOutside))) continue;

if (IO.SHOW_STD_ERR_PRIME) IO.printOut("getDeltaHistogram> to compute delta, i use alignment ("+startA+","+endA+")");

			nUsedAlignments++;
			from=Math.max(interiorStart,startA);
			to=Math.min(interiorEnd,endA);
			for (j=from; j<=to; j++) counts[j-interiorStart]++;
			if (left && isLeftMaximalB && !startAOutside) {
				lastDeltaPoint++;
				deltaPoints[lastDeltaPoint].position=startA-interiorStart;
				deltaPoints[lastDeltaPoint].mass=1;
				deltaPoints[lastDeltaPoint].nOpen=1;
				deltaPoints[lastDeltaPoint].nClosed=0;
			}
			if (!left && isRightMaximalB && !endAOutside) {
				lastDeltaPoint++;
				deltaPoints[lastDeltaPoint].position=endA-interiorStart;
				deltaPoints[lastDeltaPoint].mass=1;
				deltaPoints[lastDeltaPoint].nClosed=1;
				deltaPoints[lastDeltaPoint].nOpen=0;
			}
		}
		deltaTmp[0]=nUsedAlignments;
		deltaTmp[1]=firstAlignmentForNext!=-1?firstAlignmentForNext:i;
		if (lastDeltaPoint==-1) {
			Math.set(delta,interiorEnd-interiorStart,0.0);
			return;
		}
		
		// Estimating the density of $deltaPoints$ and updating $counts$
		lastDeltaPoint=Points.sortAndCompact(deltaPoints,lastDeltaPoint);
		if (Points.areUniformlyDistributed(deltaPoints,0,lastDeltaPoint,true,(deltaPoints[lastDeltaPoint].position-deltaPoints[0].position+1)/Points.DEFAULT_NBINS)) {
			deltaTmp[0]=0;
			Math.set(delta,interiorEnd-interiorStart,0.0);
			return;
		}
		bandwidth=Points.estimateDensity(deltaPoints,lastDeltaPoint,delta,MIN_BANDWIDTH,MAX_BANDWIDTH);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getDeltaHistogram> bandwidthSubstring="+bandwidth);
		for (i=0; i<=lastDeltaPoint; i++) {
			if (deltaPoints[i].nOpen>0) {
				from=Math.max((int)deltaPoints[i].position-bandwidth,0);
				to=(int)deltaPoints[i].position-1;
				for (j=from; j<=to; j++) counts[j]+=deltaPoints[i].nOpen;
			}
			if (deltaPoints[i].nClosed>0) {
				from=(int)deltaPoints[i].position+1;
				to=Math.min((int)deltaPoints[i].position+bandwidth,counts.length-1);
				for (j=from; j<=to; j++) counts[j]+=deltaPoints[i].nClosed;
			}
		}
		
		// Computing $delta$
		for (i=interiorStart; i<=interiorEnd; i++) {
			if (counts[i-interiorStart]>0) delta[i-interiorStart]/=counts[i-interiorStart];
		}

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("getDeltaHistogram> deltaTmp[1]="+deltaTmp[1]+" lastAlignment="+lastAlignment);
	IO.printErr("getDeltaHistogram> DELTA HISTOGRAM "+(left?"LEFT":"RIGHT")+" OF ["+start+","+interiorStart+".."+interiorEnd+","+end+"]:");
	for (int x=interiorStart; x<=interiorEnd; x++) IO.printErr(delta[x-interiorStart]);
}

	}
	
	
	/**
	 * Simplified version of $getDeltaHistogram()$ designed for short-period periodic 
	 * substrings. The procedure uses all alignments in 
	 * $sortedAlignments[firstAlignment..]$, including those that are not implied by a 
	 * periodic substring.
     *
	 * Remark: if $left=true$ (respectively, if $left=false$), an alignment is used to 
	 * compute the histogram even if it is not left- (respectively, right-) maximal. This 
	 * is because, otherwise, the histogram of left-maximal and of right-maximal events 
	 * might become very different. Assume e.g. that there is just a single other readB
	 * containing the periodic interval, and that in such readB the periodic interval is 
	 * left-maximal but not right-maximal. This induces several left-B-maximal but not 
	 * right-B-maximal alignments in readA, thus the left-histogram and the right-
	 * histogram of readA will be normalized by very different counts.
	 *
	 * Remark: the procedure assumes that $sortedAlignments[firstAlignment..]$ is sorted 
	 * by startA.
	 *
	 * @return FALSE iff the distribution of maximal events is approx uniform and no delta
	 * histogram has been computed.
	 */
	public static final void getDeltaHistogram_periodic(int start, int end, int firstAlignment, int next, boolean left, double[] delta) {
		final int MIN_BANDWIDTH = (IO.quantum<<1)/3;
		final int MAX_BANDWIDTH = IO.quantum;
		final int IDENTITY_THRESHOLD = IO.quantum;
		boolean isLeftMaximalB, isRightMaximalB;
		int i, j;
		int from, to, startA, endA, bandwidth;
		int firstAlignmentForNext, nUsedAlignments;

		// Filling $deltaPoints$
		Math.set(counts,counts.length-1,0);
		nUsedAlignments=0; firstAlignmentForNext=-1; lastDeltaPoint=-1;	
		for (i=firstAlignment; i<=lastSortedAlignment; i++) {
			endA=Alignments.alignments[sortedAlignments[i].id][4];
			if (firstAlignmentForNext==-1 && endA>=next) firstAlignmentForNext=i;
			if (endA<start) continue;
			startA=Alignments.alignments[sortedAlignments[i].id][3];
			if (startA>end) break;
			isLeftMaximalB=sortedAlignments[i].isLeftMaximalB==1;
			isRightMaximalB=sortedAlignments[i].isRightMaximalB==1;
			if ((left&&(startA<start||startA>end)) || (!left&&(endA>end||endA<start))) continue;
			nUsedAlignments++;
			from=Math.max(startA,start); to=Math.min(endA,end);
			for (j=from; j<=to; j++) counts[j-start]++;
			if (left && (isLeftMaximalB || startA<=IDENTITY_THRESHOLD)) {
				lastDeltaPoint++;
				deltaPoints[lastDeltaPoint].position=startA-start;
				deltaPoints[lastDeltaPoint].mass=1;
				deltaPoints[lastDeltaPoint].nOpen=1;
				deltaPoints[lastDeltaPoint].nClosed=0;
			}
			if (!left && (isRightMaximalB || endA>=Reads.getReadLength(id)-IDENTITY_THRESHOLD)) {
				lastDeltaPoint++;
				deltaPoints[lastDeltaPoint].position=endA-start;
				deltaPoints[lastDeltaPoint].mass=1;
				deltaPoints[lastDeltaPoint].nClosed=1;
				deltaPoints[lastDeltaPoint].nOpen=0;
			}
		}
		deltaTmp[0]=nUsedAlignments;
		deltaTmp[1]=firstAlignmentForNext!=-1?firstAlignmentForNext:i;
		if (lastDeltaPoint==-1) {
			Math.set(delta,end-start,0.0);
			return;
		}
		
		// Estimating the density of $deltaPoints$ and updating $counts$
		lastDeltaPoint=Points.sortAndCompact(deltaPoints,lastDeltaPoint);
		if (Points.areUniformlyDistributed(deltaPoints,0,lastDeltaPoint,true,(deltaPoints[lastDeltaPoint].position-deltaPoints[0].position+1)/Points.DEFAULT_NBINS)) {
			deltaTmp[0]=0;
			Math.set(delta,end-start,0.0);
			return;
		}
		bandwidth=Points.estimateDensity(deltaPoints,lastDeltaPoint,delta,MIN_BANDWIDTH,MAX_BANDWIDTH);
		for (i=0; i<=lastDeltaPoint; i++) {
			if (deltaPoints[i].nOpen>0) {
				from=Math.max((int)deltaPoints[i].position-bandwidth,0);
				to=(int)deltaPoints[i].position-1;
				for (j=from; j<=to; j++) counts[j]+=deltaPoints[i].nOpen;
			}
			if (deltaPoints[i].nClosed>0) {
				from=(int)deltaPoints[i].position+1;
				to=Math.min((int)deltaPoints[i].position+bandwidth,counts.length-1);
				for (j=from; j<=to; j++) counts[j]+=deltaPoints[i].nClosed;
			}
		}
		
		// Computing $delta$
		for (i=start; i<=end; i++) {
			if (counts[i-start]>0) delta[i-start]/=counts[i-start];
		}

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("getDeltaHistogram_periodic> deltaTmp[1]="+deltaTmp[1]+" lastDeltaPoint="+lastDeltaPoint);
	IO.printErr("getDeltaHistogram_periodic> DELTA HISTOGRAM "+(left?"LEFT":"RIGHT")+" OF ["+start+".."+end+"]:  bandwidth="+bandwidth);
	for (int x=start; x<=end; x++) IO.printErr(delta[x-start]);
}

	}
	
	
	/**
	 * Checks whether every alignment in $sortedAlignments$ points to at most one output
	 * interval.
	 */
	public static final void checkAlignmentIntervals() {
		final int THRESHOLD = IO.quantum*3;  // Arbitrary
		int i;
		int nPointers, start, end;
		
		for (i=0; i<=lastSortedAlignment; i++) {
			nPointers=0;
			if (sortedAlignments[i].periodicSubstringInterval!=null) nPointers++;
			if (sortedAlignments[i].impliedByDenseSubstring!=null) nPointers++;
			if (sortedAlignments[i].inDenseSubstring!=null) nPointers++;
			if (sortedAlignments[i].mergedToInterval!=null) nPointers++;
/*			if (nPointers>1 && sortedAlignments[i].mergedToInterval==null) {
				System.err.println("Error: alignment "+i+" points to more than one interval: ");
UNCOMMENT THIS				System.err.println("ALIGNMENT: "+sortedAlignments[i].toStringBoundaries());
				System.err.println("periodic interval: "+sortedAlignments[i].periodicSubstringInterval);
				System.err.println("impliedByDense interval: "+sortedAlignments[i].impliedByDenseSubstring);
				System.err.println("inDense interval: "+sortedAlignments[i].inDenseSubstring);
				System.err.println("alignment interval: "+sortedAlignments[i].mergedToInterval);
				System.exit(1);
			}
*/			if (sortedAlignments[i].inDenseSubstring!=null && !sortedAlignments[i].inDenseSubstring.substringReplication) {
				System.err.println("Error: alignment "+i+" has inDenseSubstring that does not point to a substring type");
				System.err.println("alignment: "+sortedAlignments[i].toStringBoundaries()+" -> "+sortedAlignments[i].toStringPointers());
				System.err.println("pointed substring: "+sortedAlignments[i].inDenseSubstring);
				System.err.println();
				System.err.println("All dense substrings:");
				for (int x=0; x<=DenseSubstrings.lastSubstring; x++) System.err.println(DenseSubstrings.substrings[x]);
				System.exit(1);
			}
/*			start=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			end=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if ( sortedAlignments[i].mergedToInterval!=null && (start<sortedAlignments[i].mergedToInterval.firstPosition-THRESHOLD || end>sortedAlignments[i].mergedToInterval.lastPosition+THRESHOLD) &&
				 !Intervals.areApproximatelyIdentical_lowQuality(start,end,sortedAlignments[i].mergedToInterval.firstPosition,sortedAlignments[i].mergedToInterval.lastPosition,id) &&
UNCOMMENT THIS				 !Intervals.isApproximatelyContained_lowQuality(start,end,sortedAlignments[i].mergedToInterval.firstPosition,sortedAlignments[i].mergedToInterval.lastPosition,id) &&
			     !Intervals.contains_lowQuality(start,end,sortedAlignments[i].mergedToInterval.firstPosition,sortedAlignments[i].mergedToInterval.lastPosition,id)
			   ) {
				System.err.println("ERROR: the following alignment points to an alignment interval with different boundaries?!");
				System.err.println("ALIGNMENT: "+ReadA.sortedAlignments[i].toStringBoundaries());
				System.err.println("SHORTER INTERVAL: "+sortedAlignments[i].mergedToInterval);
				System.exit(1);
			}
			if ( sortedAlignments[i].impliedByDenseSubstring!=null && (start<sortedAlignments[i].impliedByDenseSubstring.startA-THRESHOLD || end>sortedAlignments[i].impliedByDenseSubstring.endA+THRESHOLD) &&
			 	 !Intervals.areApproximatelyIdentical_lowQuality(start,end,sortedAlignments[i].impliedByDenseSubstring.startA,sortedAlignments[i].impliedByDenseSubstring.endA,id) &&
			 	 !Intervals.isApproximatelyContained_lowQuality(start,end,sortedAlignments[i].impliedByDenseSubstring.startA,sortedAlignments[i].impliedByDenseSubstring.endA,id) &&
				 !Intervals.contains_lowQuality(start,end,sortedAlignments[i].impliedByDenseSubstring.startA,sortedAlignments[i].impliedByDenseSubstring.endA,id)
			   ) {
				System.err.println("ERROR: the following alignment points to a dense substring with different boundaries?!");
				System.err.println("ALIGNMENT: "+ReadA.sortedAlignments[i].toStringBoundaries());
				System.err.println("SHORTER DENSE SUBSTRING: "+sortedAlignments[i].impliedByDenseSubstring);
				System.exit(1);
			}
*/			/* Commented because the boundaries of a dense substring of substring type can 
			   change a lot during detection.
			if (sortedAlignments[i].inDenseSubstring!=null && (start<sortedAlignments[i].inDenseSubstring.startA-THRESHOLD || end>sortedAlignments[i].inDenseSubstring.endA+THRESHOLD)) {
				System.err.println("ERROR: the following alignment points to a dense substring with different boundaries?!");
				System.err.println("ALIGNMENT: "+ReadA.sortedAlignments[i].toStringBoundaries());
				System.err.println("SHORTER DENSE SUBSTRING: "+sortedAlignments[i].inDenseSubstring);
				System.exit(1);
			}*/
			/* Commented because an alignment assigned to a substring interval could be 
			   reassigned to a periodic interval.
			if (sortedAlignments[i].periodicSubstringInterval!=null && (start<sortedAlignments[i].periodicSubstringInterval.firstPosition-THRESHOLD || end>sortedAlignments[i].periodicSubstringInterval.lastPosition+THRESHOLD)) {
				System.err.println("ERROR: the following alignment points to a periodic substring interval with different boundaries?!");
				System.err.println("ALIGNMENT: "+ReadA.sortedAlignments[i].toStringBoundaries());
				System.err.println("SHORTER INTERVAL: "+sortedAlignments[i].periodicSubstringInterval);
				System.exit(1);
			}*/
		}
	}
	
	
	/**
	 * Runs binary search over $sortedAlignments[firstAlignment..lastAlignment]$ using 
	 * startA as key, assuming that the array interval is sorted by startA.
	 *
	 * @return if $startA$ is not found, the procedure returns $-1-x$, where $x$ is the 
	 * first occurrence of an alignment whose startA is bigger than $startA$.
	 */
	public static final int binarySearch(int startA, int firstAlignment, int lastAlignment) {
		int from, to, middle, middleValue;
		
		from=firstAlignment; to=lastAlignment;
		while (from<=to) {
			middle=(from+to)>>1;
			middleValue=Alignments.alignments[sortedAlignments[middle].id][3];
			if (middleValue==startA) return middle;
			else if (middleValue>startA) to=middle-1;
			else if (middleValue<startA) from=middle+1;
		}
		return -1-(to+1);
	}


	public static final void printUnassignedAlignments() {
		System.err.println("Alignments not assigned to any interval:");
		for (int i=0; i<=lastSortedAlignment; i++) {	
			if (sortedAlignments[i].impliedByDenseSubstring!=null || ReadA.sortedAlignments[i].inDenseSubstring!=null || ReadA.sortedAlignments[i].mergedToInterval!=null || ReadA.sortedAlignments[i].periodicSubstringInterval!=null) continue;
			System.err.println(ReadA.sortedAlignments[i].toStringBoundaries());
		}
	}

}