package de.mpi_cbg.revant.factorize;

import java.util.Arrays;

import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Histograms;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Leaves;


public class Factors {
	/**
	 * Parameters of the pipeline
	 */
	public static int minFactorLength;  // Minimum length, smaller than $minAlignmentLength$, of a factor we aim at reconstructing.
	public static int maxFactorsPerRead;
	public static int minPeriodicCoverage;  // Minimum coverage for a periodic factor to be output. Can be smaller than $minRepeatCoverage$.
	public static int minSplitDistance = 80;  // Minimum distance between two splits for them to be considered distinct. Arbitrary.
	public static double constantCoverageThreshold;  // If the coverage histogram varies by at most this quantity in an interval, it is considered constant.
	public static int cleanWindow = 1;  // Used by $cleanHighQualityFactors$. In units of $Reads.QUALITY_SPACING$.
	public static int maxLengthBoundaryDistance = IO.quantum<<1;
	private static final double DEFAULT_MIN_COVERAGE_RATIO = 2.0;
	private static final double DEFAULT_BIAS_RATIO = 2.0;

	/**
	 * Parameters of the pipeline, estimated from the data.
	 */
	public static double minCoverageRatio;
	public static double biasRatio01, biasRatio02, biasRatio03, biasRatio31, biasRatio32, biasRatioRightLeft, biasRatioLeftRight;  // Decide whether the alignments of a factor tend to continue to its left, to its right, or in both directions.

	/**
	 * Input data structures
	 */
	public static Split[] splits;
	public static int lastSplit;

	/**
	 * Output data structures
	 */
	public static Factor[] factors;
	public static int lastFactor;

	/**
	 * Temporary space
	 */
	private static Point[] ratios, shifts, tmpShifts;
	private static int lastRatio, lastShift;
	private static Interval[] intervals1, intervals2;
	private static int lastInterval1, lastInterval2;
	private static int firstAlignmentForPreviousFactor;
	private static int firstAlignmentForCurrentFactor;
	private static int firstAlignmentForNextFactor;
	private static int lastIdenticalPeriodFirstFactor;
	private static boolean lastIdenticalPeriodOutput;
	private static double[] tmpQualities;
	private static Split tmpSplit;


	public static final void allocateMemory(int maxAlignments, int maxFactors, int maxOccurrences, int boundaryRefLength) {
		int i;
		
		constantCoverageThreshold=Math.max(1,IO.coverage-1);
		Boundaries.boundaryRefinementLength=boundaryRefLength;
		
		splits = new Split[maxFactors<<1];
		for (i=0; i<splits.length; i++) splits[i] = new Split(-1);
		factors = new Factor[maxFactors];
		for (i=0; i<maxFactors; i++) factors[i] = new Factor();
		ratios = new Point[maxFactors<<1];
		for (i=0; i<ratios.length; i++) ratios[i] = new Point();
		intervals1 = new Interval[maxAlignments];
		for (i=0; i<maxAlignments; i++) intervals1[i] = new Interval();
		lastInterval1=-1;
		intervals2 = new Interval[maxAlignments];
		for (i=0; i<maxAlignments; i++) intervals2[i] = new Interval();
		lastInterval2=-1;
		shifts = new Point[maxAlignments<<2];
		for (i=0; i<shifts.length; i++) shifts[i] = new Point();
		tmpShifts = new Point[shifts.length];
		for (i=0; i<tmpShifts.length; i++) tmpShifts[i] = new Point();
		OccurrenceTypes.allocateMemory();
		Occurrences.allocateMemory(maxAlignments,maxOccurrences);
		tmpQualities = new double[Reads.maxReadLength/Reads.QUALITY_SPACING+1];
		tmpSplit = new Split(0);
	}


	/**
	 * Stores in $factors$ all the factors of a read $R$, i.e. the subset of all
	 * substrings obtained by cutting $R$ at splits, which satisfy the following
	 * properties: (1) most positions in the substring have quality at least
	 * $Reads.qualityThreshold$; (2) the corrected coverage of the substring is at least
	 * $minRepeatCoverage$. Substrings of $R$ that do not have this corrected coverage
	 * come either from unique regions of the genome (if they have high quality),
	 * or from random substrings inserted during sequencing.
	 *
	 * Remark: a factor is the occurrence of a left-maximal or of a right-maximal repeat
	 * of the read set, but it is not necessarily the occurrence of a maximal repeat of
	 * the read set. For example, assume that string $A$ is left-maximal but not right-
	 * maximal, since it is always followed by maximal repeat $B$ in the read set: then,
	 * if $AB$ occurs in a read, $A$ is a factor of the read.
	 *
	 * Remark: a factor is not necessarily an occurrence of a submaximal repeat of the
	 * read set, and it is not necessarily an occurrence of a maximal repeat of the
	 * genome. However, every occurrence in a read $R$ of a maximal repeat of the genome
	 * is a concatenation of one or more factors, not necessarily adjacent in $R$.
	 *
	 * Remark: any parenthesization of a read is potentially correct, so we don't enforce
	 * balanced parentheses in the output.
	 */
	public static final void detect() {
		int i, nDiscardedFactors;
		
		// Backing up the original quality track
		System.arraycopy(Reads.getQualityArray(ReadA.id),0,tmpQualities,0,Reads.getQualityArrayLength(ReadA.id));		

		// Deriving factor boundaries from splitpoints
		splits2factors();

		// Refining factor boundaries and discarding low-quality factors
		markLowQualityFactors();

if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AFTER 1:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}

		cleanLowQualityFactors();
		

if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AFTER 2:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}
		
		
		cleanHighQualityFactors();
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AFTER 3:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}		
		
		
		estimateNVariants();  // Assigns $maxDiffs$ to each factor
		Occurrences.estimateCorrectedCoverage();
		if (lastFactor<=1) minCoverageRatio=DEFAULT_MIN_COVERAGE_RATIO;
		else minCoverageRatio=estimateMinCoverageRatio();
		refineBoundaries();
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AFTER 4:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}		
		
		nDiscardedFactors=discardLowQualityFactors();
		cleanFactors(nDiscardedFactors);
		if (lastFactor<0) {
			resetIDs();
			return;
		}


if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AFTER 5:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}	


		// Merging high-quality factors
		estimateNVariants();
		Occurrences.estimateCorrectedCoverage();
		if (lastFactor<=1) minCoverageRatio=DEFAULT_MIN_COVERAGE_RATIO;
		else minCoverageRatio=estimateMinCoverageRatio();
		estimatePeriods();
		computeAlignmentStatistics();
		computePeriodicSubstringStatistics();
		assignDenseSubstrings(true);  // Only substrings that contain factors
		assignPeriodicSubstringIntervals(true);  // Only intervals that contain factors
		estimateMinNAlignments();
		markBracketedFactors();
		assignAlignmentIntervals(false);
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS JUST BEFORE MERGEFACTORS:");
	for (int x=0; x<=lastFactor; x++) {
		System.err.println(factors[x]);
		System.err.println(factors[x].toStringStatistics());
		System.err.println();
	}
}		
		
		
		nDiscardedFactors=mergeFactors();
		cleanFactors(nDiscardedFactors);
		if (lastFactor<0) {
			resetIDs();
			return;
		}


if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AFTER MERGEFACTORS:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}	

		// Discarding low-coverage factors
		estimateNVariants();
		Occurrences.estimateCorrectedCoverage();		
		nDiscardedFactors=discardLowCoverageFactors();
		cleanFactors(nDiscardedFactors);
		if (lastFactor<0) {
			resetIDs();
			return;
		}


if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AFTER MERGEFACTORS PRIME:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}

		// Merging the remaining factors again
		estimateNVariants();
		Occurrences.estimateCorrectedCoverage();
		estimatePeriods();
		computeAlignmentStatistics();
		computePeriodicSubstringStatistics();
		for (i=0; i<=lastFactor; i++) {
			factors[i].keepLeftBoundary=false;
			factors[i].keepRightBoundary=false;
			factors[i].shouldBeOutput=true;
		}
		assignDenseSubstrings(false);  // Both containing and contained in a factor
		assignPeriodicSubstringIntervals(false);  // Both containing and contained
		estimateMinNAlignments();
		markBracketedFactors();
		assignAlignmentIntervals(false);
		nDiscardedFactors=mergeFactors2();
		cleanFactors(nDiscardedFactors);
		if (lastFactor<0) {
			resetIDs();
			return;
		}
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AFTER MERGEFACTORS2:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}			
		
		
		// Discarding intervals that do not comply with the final factors
		estimateNVariants();
		Occurrences.estimateCorrectedCoverage();
		estimatePeriods();
		computeAlignmentStatistics();
		computePeriodicSubstringStatistics();
		for (i=0; i<=lastFactor; i++) {
			factors[i].keepLeftBoundary=false;
			factors[i].keepRightBoundary=false;
			factors[i].shouldBeOutput=true;
		}
		assignDenseSubstrings(false);
		assignPeriodicSubstringIntervals(false);
		estimateMinNAlignments();
		markBracketedFactors();
		assignAlignmentIntervals(true);
		
		// Resetting the original quality track
		System.arraycopy(tmpQualities,0,Reads.getQualityArray(ReadA.id),0,Reads.getQualityArrayLength(ReadA.id));
		resetIDs();
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("FACTORS:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);		
}	
		
	}


	/**
	 * Creates factors by coalescing all splits that belong to the same run which spans at
	 * most $minSplitDistance$ consecutive positions. A split event is assigned to the
	 * center of mass of each run.
	 */
	private static final void splits2factors() {
		int i;
		int position, firstPositionInRun, nPositionsInRun;
		int splitPoint, previousSplitPoint;
		int nOpen, nClosed, nUnknown, previousNOpen, previousNClosed, previousNUnknown;
		double sum;

		// Ensuring the necessary order on $splits$.
		if (Split.order!=Split.POSITION) {
			Split.order=Split.POSITION;
			if (lastSplit>0) Arrays.sort(splits,0,lastSplit+1);
		}
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> SPLITS:");
	for (int x=0; x<=lastSplit; x++) System.err.print(splits[x]+" :: ");
	System.err.println();
}		

		// Coalescing splits
		lastFactor=-1;
		previousSplitPoint=0; previousNOpen=0; previousNClosed=0; previousNUnknown=0;
		firstPositionInRun=splits[0].position; sum=firstPositionInRun; nPositionsInRun=1;
		nOpen=splits[0].nOpen; nClosed=splits[0].nClosed; nUnknown=splits[0].nUnknown;
		i=1;
		while (i<=lastSplit+1) {
			if (i==lastSplit+1) position=-1;
			else {
				position=splits[i].position;
				if (position-firstPositionInRun<=minSplitDistance) {
					sum+=position; nPositionsInRun++;
					nOpen+=splits[i].nOpen; nClosed+=splits[i].nClosed; nUnknown+=splits[i].nUnknown;
					i++;
					continue;
				}
			}
			splitPoint=(int)(sum/nPositionsInRun);	
			if (splitPoint>minSplitDistance && splitPoint-1!=previousSplitPoint) {
				lastFactor++;
				factors[lastFactor].initialize(previousSplitPoint,splitPoint-1);
				factors[lastFactor].nOpenLeft=previousNOpen;
				factors[lastFactor].nClosedLeft=previousNClosed;
				factors[lastFactor].nUnknownLeft=previousNUnknown;
				factors[lastFactor].nOpenRight=nOpen;
				factors[lastFactor].nClosedRight=nClosed;
				factors[lastFactor].nUnknownRight=nUnknown;
			}
			// Next iteration
			firstPositionInRun=position; sum=firstPositionInRun; nPositionsInRun=1;
			if (splitPoint>minSplitDistance) {
				previousSplitPoint=splitPoint; 
				previousNOpen=nOpen; previousNClosed=nClosed; previousNUnknown=nUnknown;
			}
			else {
				previousNOpen+=nOpen; previousNClosed+=nClosed; previousNUnknown+=nUnknown;
			}
			if (i<=lastSplit) {
				nOpen=splits[i].nOpen; nClosed=splits[i].nClosed; nUnknown=splits[i].nUnknown;
			}
			i++;
		}
		if (previousSplitPoint!=ReadA.readLength-1) {
			lastFactor++;
			factors[lastFactor].initialize(previousSplitPoint,ReadA.readLength-1);
			factors[lastFactor].nOpenLeft=previousNOpen;
			factors[lastFactor].nClosedLeft=previousNClosed;
			factors[lastFactor].nUnknownLeft=previousNUnknown;
		}

		// At this point, factors are automatically sorted by $firstPosition$ (and also
		// by $lastPosition$).
		Factor.order=Factor.FIRSTPOSITION;
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("--> FACTORS AT THE VERY BEGINNING:");
	for (int x=0; x<=lastFactor; x++) System.err.println(factors[x]);
}		
		
	}


	/**
	 * Marks factors with not enough positions at quality at least $Reads.qualityThreshold$.
	 *
	 * @return the number of low-quality factors found.
	 */
	private static final int markLowQualityFactors() {
		int i, out;

		out=0;
		for (i=0; i<=lastFactor; i++) {
			if (Reads.isRandomInsertion(ReadA.id,factors[i].firstPosition,factors[i].lastPosition,true)) {
				factors[i].hasHighQuality=false;
				out++;
			}
		}
		return out;
	}


	/**
	 * For every low-quality factor, the procedure sets the quality of all its positions
	 * to $Reads.MAX_QUALITY$, except for its maximal prefix and its maximal suffix with
	 * quality score smaller than $Reads.MIN_RANDOM_QUALITY_SCORE$, as determined by 
	 * procedures $Histograms.maximalPrefix$ and $Histograms.maximalSuffix$.
	 *
	 * High-quality substrings fully inside low-quality factors are likely to be short,
	 * they will not be detected by the following steps of the pipeline, and they can
	 * perturb the second derivative computations of procedure $refineBoundaries$.
	 *
	 * Remark: the procedure alters $Reads.qualities$.
	 */
	private static final void cleanLowQualityFactors() {
		final int WINDOW = 4;  // Multiple of $Reads.QUALITY_SPACING$. Arbitrary.
		int i, j;
		int from, to;
		double[] array;

		for (i=0; i<=lastFactor; i++) {
			if (factors[i].hasHighQuality) continue;
			from=Histograms.maximalPrefix(Reads.getQualityArray(ReadA.id),factors[i].firstPosition/Reads.QUALITY_SPACING,factors[i].lastPosition/Reads.QUALITY_SPACING,false,Reads.MIN_RANDOM_QUALITY_SCORE,Math.min(WINDOW,Math.max(factors[i].length>>2,1)));
			from++;
			to=Histograms.maximalSuffix(Reads.getQualityArray(ReadA.id),factors[i].firstPosition/Reads.QUALITY_SPACING,factors[i].lastPosition/Reads.QUALITY_SPACING,false,Reads.MIN_RANDOM_QUALITY_SCORE,Math.min(WINDOW,Math.max(factors[i].length>>2,1)));
			to--;
			array=Reads.getQualityArray(ReadA.id);
			for (j=from; j<=to; j++) array[j]=Reads.MAX_QUALITY_SCORE;
		}
	}


	/**
	 * Let $[i..j]$ be a maximal substring such that no position in $[i..j]$ has high
	 * quality, and such that the average quality of a window of size $cleanWindow$ before
	 * $i$ and after $j$ is high. The procedure removes all such substrings $[i..j]$ from
	 * all high-quality factors. Such substrings are likely to be short, they will not
	 * be detected by the following steps of the pipeline, and they can perturb the second
	 * derivative computations of $refineBoundaries$.
	 *
	 * Remark: the procedure alters $Reads.qualities$.
	 */
	private static final void cleanHighQualityFactors() {
		int i, j, k;
		int firstPosition, newQuality;
		double avgQualityLeft, avgQualityRight;
		final double[] array = Reads.getQualityArray(ReadA.id);

		for (i=0; i<=lastFactor; i++) {
			if (!factors[i].hasHighQuality) continue;
			firstPosition=-1;
			for (j=factors[i].firstPosition/Reads.QUALITY_SPACING+cleanWindow; j<=factors[i].lastPosition/Reads.QUALITY_SPACING-cleanWindow; j++) {
				if (array[j]<=Reads.MAX_HIGH_QUALITY_SCORE) continue;
				avgQualityLeft=Histograms.getAverage(array,j-cleanWindow,j-1);
				if (firstPosition==-1 && avgQualityLeft<=Reads.MAX_HIGH_QUALITY_SCORE) firstPosition=j;
				avgQualityRight=Histograms.getAverage(array,j+1,j+cleanWindow);
				if (firstPosition!=-1 && avgQualityRight<=Reads.MAX_HIGH_QUALITY_SCORE) {
					newQuality=(int)((avgQualityLeft+avgQualityRight)/2);
					for (k=firstPosition; k<=j; k++) array[k]=newQuality;
					firstPosition=-1;
				}
			}
		}
	}


	/**
	 * Uses the second derivative of the coverage histogram and of the quality histogram
	 * to adjust the boundaries between two consecutive factors, as well as the left
	 * boundary of the first factor and the right boundary of the last factor. See
	 * procedure $Boundaries.refineBoundary()$.
	 */
	private static final void refineBoundaries() {
		int i, start, end;
		double fraction;
		Factor previousFactor, nextFactor;

		// Degenerate case: exactly one factor.
		if (lastFactor==0) {
			previousFactor=null; nextFactor=null;
			Boundaries.refineBoundary(factors[0],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,Reads.getQualityArray(ReadA.id),1,false,true);
			Boundaries.refineBoundary(factors[0],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,ReadA.coverage,1,true,false);
			Boundaries.refineBoundary(factors[0],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,false,Reads.getQualityArray(ReadA.id),-1,false,true);
			Boundaries.refineBoundary(factors[0],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,false,ReadA.coverage,-1,true,false);
			return;
		}

		// First factor
		previousFactor=null; nextFactor=factors[1];
		if (factors[0].hasHighCoverage) {
			Boundaries.refineBoundary(factors[0],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,Reads.getQualityArray(ReadA.id),1,false,true);
			Boundaries.refineBoundary(factors[0],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,ReadA.coverage,1,true,false);
		}
		else if (lastFactor>=1 && factors[1].hasHighCoverage) {
			start=factors[0].lastPosition-(Boundaries.boundaryRefinementLength<factors[0].length?Boundaries.boundaryRefinementLength:factors[0].length>>1);
			fraction=ReadA.getFractionAtCoverage(start,factors[0].lastPosition,factors[1].getCorrectedCoverage());
			if (fraction>0) Boundaries.refineBoundary(factors[0],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,false,ReadA.coverage,-1,false,false);
			if (!factors[0].hasHighQuality) {
				fraction=Histograms.getFraction(Reads.getQualityArray(ReadA.id),start/Reads.QUALITY_SPACING,factors[0].lastPosition/Reads.QUALITY_SPACING,Reads.MIN_RANDOM_QUALITY_SCORE,false);
				if (fraction>0) Boundaries.refineBoundary(factors[0],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,false,Reads.getQualityArray(ReadA.id),-1,true,true);
			}
		}

		// All intermediate low-coverage factors
		for (i=1; i<lastFactor; i++) {
			if (factors[i].hasHighCoverage) continue;
			previousFactor=factors[i-1]; nextFactor=factors[i+1];
			if (Boundaries.boundaryRefinementLength<factors[i].length) {
				start=factors[i].lastPosition-Boundaries.boundaryRefinementLength;
				end=factors[i].firstPosition+Boundaries.boundaryRefinementLength;
			}
			else {
				start=factors[i].lastPosition-(factors[i].length>>1);
				end=factors[i].firstPosition+(factors[i].length>>1);
			}
			if (factors[i-1].hasHighCoverage) {
				fraction=ReadA.getFractionAtCoverage(factors[i].firstPosition,end,factors[i-1].getCorrectedCoverage());
				if (fraction>0) Boundaries.refineBoundary(factors[i],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,ReadA.coverage,1,false,false);
			}
			if (factors[i+1].hasHighCoverage) {
				fraction=ReadA.getFractionAtCoverage(start,factors[i].lastPosition,factors[i+1].getCorrectedCoverage());
				if (fraction>0) Boundaries.refineBoundary(factors[i],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,false,ReadA.coverage,-1,false,false);
			}
			if (factors[i].hasHighQuality) continue;
			if (factors[i-1].hasHighQuality) {
				fraction=Histograms.getFraction(Reads.getQualityArray(ReadA.id),factors[i].firstPosition/Reads.QUALITY_SPACING,end/Reads.QUALITY_SPACING,Reads.MIN_RANDOM_QUALITY_SCORE,false);
				if (fraction>0) Boundaries.refineBoundary(factors[i],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,Reads.getQualityArray(ReadA.id),1,true,true);
			}
			if (factors[i+1].hasHighQuality) {
				fraction=Histograms.getFraction(Reads.getQualityArray(ReadA.id),start/Reads.QUALITY_SPACING,factors[i].lastPosition/Reads.QUALITY_SPACING,Reads.MIN_RANDOM_QUALITY_SCORE,false);
				if (fraction>0) Boundaries.refineBoundary(factors[i],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,false,Reads.getQualityArray(ReadA.id),-1,true,true);
			}
		}

		// Last factor
		previousFactor=factors[lastFactor-1]; nextFactor=null;
		if (factors[lastFactor].hasHighCoverage) {
			Boundaries.refineBoundary(factors[lastFactor],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,false,Reads.getQualityArray(ReadA.id),-1,false,true);
			Boundaries.refineBoundary(factors[lastFactor],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,false,ReadA.coverage,-1,true,false);
		}
		else if (lastFactor>=1 && factors[lastFactor-1].hasHighCoverage) {
			end=factors[lastFactor].firstPosition+(Boundaries.boundaryRefinementLength<factors[lastFactor].length?Boundaries.boundaryRefinementLength:factors[lastFactor].length>>1);
			fraction=ReadA.getFractionAtCoverage(factors[lastFactor].firstPosition,end,factors[lastFactor-1].getCorrectedCoverage());
			if (fraction>0) Boundaries.refineBoundary(factors[lastFactor],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,ReadA.coverage,1,false,false);
			if (!factors[lastFactor].hasHighQuality) {
				fraction=Histograms.getFraction(Reads.getQualityArray(ReadA.id),factors[lastFactor].firstPosition/Reads.QUALITY_SPACING,end/Reads.QUALITY_SPACING,Reads.MIN_RANDOM_QUALITY_SCORE,false);
				if (fraction>0) Boundaries.refineBoundary(factors[lastFactor],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,Reads.getQualityArray(ReadA.id),1,true,true);
			}
		}

		// All intermediate pairs of high-coverage factors
		for (i=1; i<=lastFactor; i++) {
			if (!factors[i].hasHighCoverage || !factors[i-1].hasHighCoverage) continue;
			previousFactor=factors[i-1]; nextFactor=null;
			if ( factors[i].getCorrectedCoverage()<factors[i-1].getCorrectedCoverage()*minCoverageRatio &&
			     factors[i-1].getCorrectedCoverage()<factors[i].getCorrectedCoverage()*minCoverageRatio ) continue;
			Boundaries.refineBoundary(factors[i],previousFactor,nextFactor,factors[0].firstPosition,factors[lastFactor].lastPosition,true,ReadA.coverage,0,false,false);
		}
	}


	/**
	 * Computing alignment statistics for each factor ($nAlignments$,
	 * $nMaximalAlignments$, $leftScore$, $rightScore$).
	 *
	 * Remark: the procedure does not assume any order of $ReadA.sortedAlignments$.
	 */
	private static final void computeAlignmentStatistics() {
		boolean implied;
		int i, j, k;
		int alignmentStart, alignmentEnd;
		int firstJForCurrentAlignment, firstJForNextAlignment;

		// Ensuring the necessary order in $ReadA.sortedAlignments$
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

		// Computing statistics
		for (i=0; i<=lastFactor; i++) factors[i].clearNAlignments();
		i=0; j=0; firstJForCurrentAlignment=0; firstJForNextAlignment=-1;
		while (i<=ReadA.lastSortedAlignment) {
			alignmentStart=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			if (alignmentStart>factors[lastFactor].lastPosition) break;
			alignmentEnd=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (alignmentEnd<factors[0].firstPosition) {
				i++;
				continue;
			}
			if (j>lastFactor || factors[j].firstPosition>alignmentEnd) {
				implied=ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null||ReadA.sortedAlignments[i].impliedByDenseSubstring!=null;
				for (k=firstJForCurrentAlignment; k<=j-1; k++) factors[k].updateAlignmentStatistics(alignmentStart,alignmentEnd,ReadA.sortedAlignments[i].isLeftMaximal==1,ReadA.sortedAlignments[i].isRightMaximal==1,implied,k==0?null:factors[k-1],k==lastFactor?null:factors[k+1]);
				i++;
				if (firstJForNextAlignment!=-1) {
					j=firstJForNextAlignment;
					firstJForNextAlignment=-1;
					firstJForCurrentAlignment=j;
				}
				continue;
			}
			if (factors[j].lastPosition<alignmentStart) {
				j++;
				firstJForCurrentAlignment=j;
				continue;
			}
			if (i<ReadA.lastSortedAlignment && firstJForNextAlignment==-1 && factors[j].lastPosition>=Alignments.alignments[ReadA.sortedAlignments[i+1].id][3]) firstJForNextAlignment=j;
			j++;
		}
	}


	/**
	 * Computing, for each factor, the counts based on periodic substrings that are used
	 * to decide oversplits.
	 *
	 * Remark: the procedure does not assume any order of $PeriodicSubstrings.substrings$.
	 */
	private static final void computePeriodicSubstringStatistics() {
		int i, j, k;
		int substringStart, substringEnd;
		int firstJForCurrentSubstring, firstJForNextSubstring;

		// Ensuring the necessary order in $PeriodicSubstrings.substrings$
		if (PeriodicSubstring.order!=PeriodicSubstring.MINSTARTA && PeriodicSubstring.order!=PeriodicSubstring.MINSTARTA_MAXENDA) {
			PeriodicSubstring.order=PeriodicSubstring.MINSTARTA;
			if (PeriodicSubstrings.lastSubstring>0) Arrays.sort(PeriodicSubstrings.substrings,0,PeriodicSubstrings.lastSubstring+1);
		}

		// Computing statistics
		for (i=0; i<=lastFactor; i++) {
			factors[i].periodicSubstringLeft=false;
			factors[i].periodicSubstringRight=false;
			factors[i].periodicSubstringSpanning=false;
		}
		i=0; j=0; firstJForCurrentSubstring=j; firstJForNextSubstring=-1;
		while (i<=PeriodicSubstrings.lastSubstring) {
			substringStart=PeriodicSubstrings.substrings[i].minStartA;
			if (substringStart>factors[lastFactor].lastPosition) break;
			substringEnd=PeriodicSubstrings.substrings[i].maxEndA;
			if (substringEnd<factors[0].firstPosition) {
				i++;
				continue;
			}
			if (j>lastFactor || factors[j].firstPosition>substringEnd) {
				for (k=firstJForCurrentSubstring; k<=j-1; k++) factors[k].updatePeriodicSubstringStatistics(substringStart,substringEnd);
				i++;
				if (firstJForNextSubstring!=-1) {
					j=firstJForNextSubstring;
					firstJForNextSubstring=-1;
					firstJForCurrentSubstring=j;
				}
				continue;
			}
			if (factors[j].lastPosition<substringStart) {
				j++;
				firstJForCurrentSubstring=j;
				continue;
			}
			if (firstJForNextSubstring==-1 && i<PeriodicSubstrings.lastSubstring && factors[j].lastPosition>=PeriodicSubstrings.substrings[i+1].minStartA) firstJForNextSubstring=j;
			j++;
		}
	}


	/**
	 * Marks a factor that belongs to a dense substring of substring type as left-
	 * (respectively, right-) bracketed iff it has a large enough number of left-maximal 
	 * (right-maximal) alignments.
	 *
	 * Remark: the procedure assumes that $minNAlignmentsLeft$ and $minNAlignmentsRight$ 
	 * have already been computed for all dense substrings of substring type.
	 */
	private static final void markBracketedFactors() {
		int i, observed;

		for (i=0; i<=lastFactor; i++) {
			factors[i].isLeftBracketed=false;
			factors[i].isRightBracketed=false;
			if (!factors[i].inDenseSubstringOfSubstringType) continue;

			// Strict notion
			observed=factors[i].nMaximalAlignments[0]+factors[i].nMaximalAlignments[2]+factors[i].nMaximalAlignments[3];
			if (observed>=Factor.minMaximalAlignments && factors[i].hasLargeNAlignments(true,observed)) factors[i].isLeftBracketed=true;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markBracketedFactors> factor "+i+"=["+factors[i].firstPosition+".."+factors[i].lastPosition+"] isLeftBracketed="+factors[i].isLeftBracketed+" observed="+observed);
			observed=factors[i].nMaximalAlignments[1]+factors[i].nMaximalAlignments[2]+factors[i].nMaximalAlignments[4];
			if (observed>=Factor.minMaximalAlignments && factors[i].hasLargeNAlignments(false,observed)) factors[i].isRightBracketed=true;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markBracketedFactors> factor "+i+"=["+factors[i].firstPosition+".."+factors[i].lastPosition+"] isRightBracketed="+factors[i].isRightBracketed+" observed="+observed);
			
			// Relaxed notion
			if (!factors[i].isLeftBracketed && factors[i].isRightBracketed) {
				observed=factors[i].nMaximalAlignments[0]+factors[i].nMaximalAlignments[1]+factors[i].nMaximalAlignments[2]+factors[i].nMaximalAlignments[3];
				if (factors[i].hasLargeNAlignments(true,observed)) factors[i].isLeftBracketed=true;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markBracketedFactors> factor "+i+"=["+factors[i].firstPosition+".."+factors[i].lastPosition+"] isLeftBracketed="+factors[i].isLeftBracketed+" observed="+observed);
			}
			else if (!factors[i].isRightBracketed && factors[i].isLeftBracketed) {
				observed=factors[i].nMaximalAlignments[0]+factors[i].nMaximalAlignments[1]+factors[i].nMaximalAlignments[2]+factors[i].nMaximalAlignments[4];
				if (factors[i].hasLargeNAlignments(false,observed)) factors[i].isRightBracketed=true;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("markBracketedFactors> factor "+i+"=["+factors[i].firstPosition+".."+factors[i].lastPosition+"] isRightBracketed="+factors[i].isRightBracketed+" observed="+observed);
			}
		}
	}


	/**
	 * Removes splits, between consecutive factors, that are likely to be
	 * wrong. The splits that form the left and the right boundary of a factor are
	 * discarded using maximality, type of dense substring, length biases in dense
	 * substrings, maximal alignment statistics, alignment statistics, coverage, period,
	 * and internal symmetries, of the factor and of its neighbors. The procedure assumes
	 * that the previous steps of the pipeline can only oversplit a factor, and never
	 * undersplit a factor.
	 *
	 * Remark: the procedure assumes that a number of other procedures have been already
	 * executed. See $detect$ for details.
	 *
	 * @return the number of factors that are discarded as a result of merging.
	 */
	private static final int mergeFactors() {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i, out;
		int keepLeftSplit;  // -1: don't keep the left boundary of the current factor. 1: keep it. 0: can't say.
		int firstFactorInRun;

		// Ensuring the necessary order in $ReadA.sortedAlignments$.
		// This is required just by procedure $haveIdenticalPeriod$.
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

		// Estimating $biasRatio*$
		biasRatio01=estimateBiasRatio(0);
		biasRatio02=estimateBiasRatio(1);
		biasRatio03=estimateBiasRatio(2);
		biasRatio31=estimateBiasRatio(3);
		biasRatio32=estimateBiasRatio(4);
		biasRatioRightLeft=estimateBiasRatio(5);
		biasRatioLeftRight=estimateBiasRatio(6);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("biasRatio01="+biasRatio01+" biasRatio02="+biasRatio02+" biasRatio03="+biasRatio03+" biasRatio31="+biasRatio31+" biasRatio32="+biasRatio32+" biasRatioRightLeft="+biasRatioRightLeft+" biasRatioLeftRight="+biasRatioLeftRight);

		// Merging
		out=0;
		firstFactorInRun=0; keepLeftSplit=1;
		lastIdenticalPeriodFirstFactor=Math.POSITIVE_INFINITY;
		firstAlignmentForPreviousFactor=0;
		firstAlignmentForCurrentFactor=0;
		firstAlignmentForNextFactor=0;
		for (i=0; i<=lastFactor; i++) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("mergeHighCoverageFactors> factor "+i+" isLeftMaximal="+factors[i].isLeftMaximal()+" isRightMaximal="+factors[i].isRightMaximal()+" inDenseSubstringOfSubstringType="+factors[i].inDenseSubstringOfSubstringType);
			firstAlignmentForPreviousFactor=firstAlignmentForCurrentFactor;
			firstAlignmentForCurrentFactor=firstAlignmentForNextFactor;
			if (keepLeftBoundary(factors[i],i>0?factors[i-1]:null,i,keepLeftSplit)) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("mergeHighCoverageFactors> factor "+i+": left boundary kept");
				firstFactorInRun=i;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("mergeHighCoverageFactors> firstFactorInRun="+firstFactorInRun);
			}
			else {
				factors[i].shouldBeOutput=false;
				out++;
				factors[firstFactorInRun].lastPosition=factors[i].lastPosition;
				factors[firstFactorInRun].length=factors[firstFactorInRun].lastPosition-factors[firstFactorInRun].firstPosition+1;
				factors[firstFactorInRun].isRightBracketed=factors[i].isRightBracketed;
				factors[firstFactorInRun].nOpenRight=factors[i].nOpenRight;
				factors[firstFactorInRun].nClosedRight=factors[i].nClosedRight;
				factors[firstFactorInRun].nUnknownRight=factors[i].nUnknownRight;
				if (factors[i].firstPosition<=factors[firstFactorInRun].firstPosition+IDENTITY_THRESHOLD) {
					factors[firstFactorInRun].nOpenLeft+=factors[i].nOpenLeft;
					factors[firstFactorInRun].nClosedLeft+=factors[i].nClosedLeft;
					factors[firstFactorInRun].nUnknownLeft+=factors[i].nUnknownLeft;
				}
				
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("mergeHighCoverageFactors> factor "+i+": left boundary not kept, firstFactorInRun="+firstFactorInRun);
			}
			keepLeftSplit=keepRightBoundary(factors[i],i<lastFactor?factors[i+1]:null,i,lastFactor);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("mergeHighCoverageFactors> factor "+i+": right boundary = "+keepLeftSplit);
		}
		return out;
	}

	
 	/**
 	 * Estimates the minimum number, of maximal alignments of a factor that are not 
 	 * implied by a dense or periodic substring, that signals a boundary inside a dense 
 	 * substring that replicates by taking substrings of itself. See procedure 
 	 * $estimateFromRatios$ for details.
 	 *
 	 * Remark: factors that are too close to the beginning (end) of the dense substring
 	 * are not used to estimate the minimum number of alignments of left (right) 
 	 * boundaries.
 	 *
 	 * Remark: the procedure takes the conservative approach of not using factor start/end
 	 * positions that are too close to the start/end of a periodic substring. It is true 
 	 * that variable $nMaximalAlignments$ of a factor counts just alignments that are not 
 	 * implied, however such alignments inside a periodic substring can be noise caused by 
 	 * the fact that the aligner does not report all alignments in periodic substrings. 
 	 * This is similar to what is done in $ReadA.getDeltaHistogram()$. See 
 	 * $PeriodicSubstrings.atPeriodicSubstringEnd()$. Note that the same cannot be said 
 	 * for other types of substring, since they are not known to induce noise alignments.
	 */
	private static final void estimateMinNAlignments() {
		for (int i=0; i<=DenseSubstrings.lastSubstring; i++) {
			if (!DenseSubstrings.substrings[i].substringReplication) {
				DenseSubstrings.substrings[i].minNAlignmentsLeft=-1;
				DenseSubstrings.substrings[i].minNAlignmentsRight=-1;
				continue;
			}
			DenseSubstrings.substrings[i].minNAlignmentsLeft=estimateMinNAlignments_impl(DenseSubstrings.substrings[i].id,true);
			DenseSubstrings.substrings[i].minNAlignmentsRight=estimateMinNAlignments_impl(DenseSubstrings.substrings[i].id,false);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("estimateMinRatios> dense substring "+DenseSubstrings.substrings[i].id+"=["+DenseSubstrings.substrings[i].startA+".."+DenseSubstrings.substrings[i].endA+"] minNAlignmentsLeft="+DenseSubstrings.substrings[i].minNAlignmentsLeft+" minNAlignmentsRight="+DenseSubstrings.substrings[i].minNAlignmentsRight);
		}
	}
	

	/*
	 * @return -3 if the dense substring contains at most one factor; -2 if the dense 
	 * substring contains no maximal factor; -1 if no useful ratio exists.
	 */
	private static final double estimateMinNAlignments_impl(int denseSubstringID, boolean left) {
		final int N_CONSECUTIVE_POINTS = 2;
		final double QUANTILE = 0.75;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		final double DEFAULT_MIN_N_ALIGNMENTS_MIN = IO.coverage<<1;  // Arbitrary
		final double DEFAULT_MIN_N_ALIGNMENTS_MAX = IO.coverage*20;  // Arbitrary
		final int MAX_DIFFERENCE = IO.coverage<<1;
		final int THRESHOLD = IO.quantum<<1;
		final int minIntervalLength, minLocalMaxDistance;
		int i, j;
		int nFactors, nMaximalFactors;
		double out;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("beginning of estimateMinNAlignments("+denseSubstringID+")");
		lastRatio=-1; nFactors=0; nMaximalFactors=0;
		for (i=0; i<=lastFactor; i++) {
			j=factors[i].inDenseSubstring(denseSubstringID);
			if ( j==-1 ||
			     !( Intervals.isApproximatelyContained(factors[i].firstPosition,factors[i].lastPosition,factors[i].denseSubstrings[j].startA,factors[i].denseSubstrings[j].endA) &&
			        !Intervals.areApproximatelyIdentical(factors[i].firstPosition,factors[i].lastPosition,factors[i].denseSubstrings[j].startA,factors[i].denseSubstrings[j].endA)
			      )
			   ) continue;
			nFactors++;
			if ((left && factors[i].isLeftMaximal()) || (!left && factors[i].isRightMaximal())) nMaximalFactors++;
			// Left boundary
			if ( left && factors[i].isLeftMaximal() && 
			     Math.abs(factors[i].firstPosition,factors[i].denseSubstrings[j].sumStartA/factors[i].denseSubstrings[j].nStartA)>THRESHOLD &&
				 Math.abs(factors[i].firstPosition,factors[i].denseSubstrings[j].sumEndA/factors[i].denseSubstrings[j].nEndA)>=Alignments.minAlignmentLength &&
			     !PeriodicSubstrings.atPeriodicSubstringEnd(factors[i].firstPosition,true,THRESHOLD)
			   ) {				   
				lastRatio++;
				ratios[lastRatio].position=factors[i].nMaximalAlignments[0]+factors[i].nMaximalAlignments[2]+factors[i].nMaximalAlignments[3];
				ratios[lastRatio].mass=1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("added left nAlignments = "+ratios[lastRatio].position+" by factor ["+factors[i].firstPosition+".."+factors[i].lastPosition+"] for dense substring id="+denseSubstringID);				
			}
			// Right boundary
			else if ( !left && factors[i].isRightMaximal() && 
			          Math.abs(factors[i].lastPosition,factors[i].denseSubstrings[j].sumEndA/factors[i].denseSubstrings[j].nEndA)>THRESHOLD &&
					  Math.abs(factors[i].lastPosition,factors[i].denseSubstrings[j].sumStartA/factors[i].denseSubstrings[j].nStartA)>=Alignments.minAlignmentLength &&
					  !PeriodicSubstrings.atPeriodicSubstringEnd(factors[i].lastPosition,false,THRESHOLD)
			        ) {
				lastRatio++;
				ratios[lastRatio].position=factors[i].nMaximalAlignments[1]+factors[i].nMaximalAlignments[2]+factors[i].nMaximalAlignments[4];
				ratios[lastRatio].mass=1;
			}
		}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("nFactors="+nFactors+" lastRatio="+lastRatio);
		if (nFactors<=1) return -3;
		if (nMaximalFactors==0) return -2;
		if (lastRatio==-1) return -1;
		if (lastRatio==0) return DEFAULT_MIN_N_ALIGNMENTS_MIN;
		lastRatio=Points.sortAndCompact(ratios,lastRatio);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("NMAXIMALALIGNMENTS "+(left?"LEFT":"RIGHT")+" for dense substring "+denseSubstringID);
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=lastRatio; x++) IO.printErr(ratios[x]); }
		if (lastRatio==0) return DEFAULT_MIN_N_ALIGNMENTS_MIN;
		minIntervalLength=Math.max( (int)Points.distanceQuantile(ratios,lastRatio,N_CONSECUTIVE_POINTS,true,QUANTILE),
									(int)((ratios[lastRatio].position-ratios[0].position)/MIN_INTERVAL_LENGTH_FACTOR)
								  );
		minLocalMaxDistance=minIntervalLength<<1;
		out=estimateFromRatios(minIntervalLength,minLocalMaxDistance,DEFAULT_MIN_N_ALIGNMENTS_MIN,MAX_DIFFERENCE,true);
		if (out>DEFAULT_MIN_N_ALIGNMENTS_MAX) out=DEFAULT_MIN_N_ALIGNMENTS_MAX;
		return out;
	}


	/**
	 * Fits a density estimation tree on $ratios$, and returns a separation point between
	 * the first and the second local maximum. The first local maximum is assumed to
	 * correspond to random ratios.
	 *
	 * @param maxDifference if there are "few" values, and if such values differ by at 
	 * most this quantity, the procedure assumes they are uniformly distributed and does
	 * not try to split them.
	 */
	private static final double estimateFromRatios(double minIntervalLength, double minLocalMaxDistance, double defaultRatio, double maxDifference, boolean integerPositions) {
		int tmp, nLocalMaximumLeaves;

		if (lastRatio+1<=Points.FEW_POINTS) {
			if (ratios[lastRatio].position-ratios[0].position<=maxDifference) return defaultRatio;
			tmp=Points.getRoughThreshold(ratios,lastRatio,true,0,false);
 if (IO.SHOW_STD_ERR_PRIME) IO.printErr("---> 4 tmp="+tmp+" ratios[tmp]="+ratios[tmp].position);
			return tmp==-1?defaultRatio:ratios[tmp].position;
		}
		if (Points.areUniformlyDistributed(ratios,0,lastRatio,false,(ratios[lastRatio].position-ratios[0].position)/Points.DEFAULT_NBINS)) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("estimateFromRatios> case 0");
			return defaultRatio;
		}
		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(ratios,0,lastRatio,minIntervalLength,minLocalMaxDistance,integerPositions,-1,-1,-1,false,true);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("estimateFromRatios> TREE OF RATIOS:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr(DensityEstimationTree.leaves[x]); }
		if (nLocalMaximumLeaves==0) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("estimateFromRatios> case 1");
			return defaultRatio;
		}
		DensityEstimationTree.markRunsOfLocalMaximumLeaves(ratios);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("estimateFromRatios> AFTER MARKRUNS:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr(DensityEstimationTree.leaves[x]); }
		if (DensityEstimationTree.nRuns==1 && DensityEstimationTree.leaves[Leaves.lastLeafOfRun(DensityEstimationTree.leaves,DensityEstimationTree.lastLeaf,0)].lastPoint==lastRatio) {
			// Separating the only run from the points \emph{before it}.
			return DensityEstimationTree.separateRun(-1,ratios,0,lastRatio,true);
		}
		return DensityEstimationTree.separateRun(0,ratios,0,lastRatio,true);
	}


	/**
	 * Estimates the minimum ratio, between the corrected coverage of two adjacent
	 * factors, that signals a boundary between such factors. See procedure
	 * $estimateFromRatios$ for details.
	 *
	 * @return -1 if no useful ratio exists.
	 */
	private static final double estimateMinCoverageRatio() {
		final double MIN_INTERVAL_LENGTH = 1.0;
		final double MIN_LOCAL_MAX_DISTANCE = 2.0;
		final double MAX_DIFFERENCE = 2.0;
		int i;
		double t, min, max;

		lastRatio=-1;
		for (i=1; i<=lastFactor; i++) {
			min=factors[i-1].getCorrectedCoverage();
			max=factors[i].getCorrectedCoverage();
			if (min<=0 || max<=0) continue;
			if (max<min) {
				t=min;
				min=max;
				max=t;
			}
			lastRatio++;
			ratios[lastRatio].position=max/min;
			ratios[lastRatio].mass=1;
		}
		if (lastRatio==-1) return -1;
		if (lastRatio==0) return DEFAULT_MIN_COVERAGE_RATIO;
		lastRatio=Points.sortAndCompact(ratios,lastRatio);
if (IO.SHOW_STD_ERR) IO.printErr("COVERAGE RATIOS:  lastRatio="+lastRatio);
if (IO.SHOW_STD_ERR) { for (int x=0; x<=lastRatio; x++) IO.printErr(ratios[x]); }
		if (lastRatio==0) return DEFAULT_MIN_COVERAGE_RATIO;
		return estimateFromRatios(MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE,DEFAULT_MIN_COVERAGE_RATIO,MAX_DIFFERENCE,false);
	}


	/**
	 * Estimates the smallest ratio >1, between two elements of array $Factor.nAlignments$,
	 * that makes a factor belong to a given type (inside factor, spanning factor, right
	 * factor, left factor), as used in procedures $keep{Left,Right}Boundary$.
	 *
	 * @param type: 0=0/1, 1=0/2, 2=0/3, 3=3/1, 4=3/2, 5=right/left, 6=left/right.
	 */
	private static final double estimateBiasRatio(int type) {
		final double MIN_INTERVAL_LENGTH = 1.0;
		final double MIN_LOCAL_MAX_DISTANCE = 2.0;
		final double MAX_DIFFERENCE = 2.0;
		int i;
		double n=0, d=0;

		lastRatio=-1;
		for (i=0; i<=lastFactor; i++) {
			if (type==0) {
				n=factors[i].nAlignments[0];
				d=factors[i].nAlignments[1];
			}
			else if (type==1) {
				n=factors[i].nAlignments[0];
				d=factors[i].nAlignments[2];
			}
			else if (type==2) {
				n=factors[i].nAlignments[0];
				d=factors[i].nAlignments[3];
			}
			else if (type==3) {
				n=factors[i].nAlignments[3];
				d=factors[i].nAlignments[1];
			}
			else if (type==4) {
				n=factors[i].nAlignments[3];
				d=factors[i].nAlignments[2];
			}
			else if (type==5) {
				n=factors[i].rightScore;
				d=factors[i].leftScore;
			}
			else if (type==6) {
				n=factors[i].leftScore;
				d=factors[i].rightScore;
			}
			if (n<d || n==0) continue;
			if (d==0) d=1;
			lastRatio++;
			ratios[lastRatio].position=n/d;
			ratios[lastRatio].mass=1;
		}
		if (lastRatio<=0) return DEFAULT_BIAS_RATIO;
		lastRatio=Points.sortAndCompact(ratios,lastRatio);
if (IO.SHOW_STD_ERR) IO.printErr("BIAS RATIOS, type "+type+":");
if (IO.SHOW_STD_ERR) { for (int x=0; x<=lastRatio; x++) IO.printErr(ratios[x]); }
		if (lastRatio==0) return DEFAULT_BIAS_RATIO;
		return estimateFromRatios(MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE,DEFAULT_BIAS_RATIO,MAX_DIFFERENCE,false);
	}


	/**
	 * Decides whether to keep the left boundary of a factor, using maximality,
	 * type of dense substrings, biases in dense substring lengths, maximal alignment
	 * statistics, alignment statistics, coverage, period, and internal symmetries, of the
	 * factor and of the factor to its left.
	 *
	 * Remark: we could always merge factors whose length is smaller than a minimum factor
	 * length. For full generality, the current procedure does not do that.
	 *
	 * Remark: a factor without a strong bias in its covering alignments could also come
	 * from undersplitting. This procedure assumes that just oversplitting happened
	 * upstream.
	 */
	private static final boolean keepLeftBoundary(Factor factor, Factor previousFactor, int factorID, int keepLeftSplit) {
		boolean out;
		int i;
		int observed, expected;
		int tmp1, tmp2;

		if (factorID==0) return true;
		if (keepLeftSplit==1) return true;
		if (factor.firstPosition>previousFactor.lastPosition+1) return true;
if (IO.SHOW_STD_ERR) IO.printErr(factor.firstPosition+" inLengthBiasInterval? "+DenseSubstrings.inLengthBiasInterval(factor.firstPosition,maxLengthBoundaryDistance));
		if (DenseSubstrings.inLengthBiasInterval(factor.firstPosition,maxLengthBoundaryDistance)) return true;
	 	if (factor.isLeftMaximal()) {
	 		if (!factor.inDenseSubstringOfSuffixType && !factor.inDenseSubstringOfSubstringType && factor.periodicCoverage==0) {
				// Left-maximality is enough to keep the left boundary, regardless of the
				// vote of the previous factor.
if (IO.SHOW_STD_ERR) IO.printErr("      keepLeftBoundary> case 1");
				return true;
			}
			else if (factor.inDenseSubstringOfSubstringType) {
				// Left-maximality is not enough. If the number of left-maximal alignments
				// at the boundary is high enough, we keep the left split, regardless of
				// the vote of the previous factor.
				if (factor.isLeftBracketed) {
if (IO.SHOW_STD_ERR) IO.printErr("      keepLeftBoundary> case 2");
					return true;
				}
			}
			else {
				// In all other cases we fall back to the standard merging below, since
				// we can separate skeleton alignments from all other alignments.
			}
		}
		// The factor is not left-maximal, or left-maximality is not enough to decide the
		// boundary.
		if ( factor.nAlignments[0]>factor.nAlignments[1]*biasRatio01 &&
			 factor.nAlignments[0]>factor.nAlignments[2]*biasRatio02 &&
			 factor.nAlignments[0]>factor.nAlignments[3]*biasRatio03 ) {
			// Factor with most alignments inside
if (IO.SHOW_STD_ERR) IO.printErr("      keepLeftBoundary> case 3");
			return true;
		}
		else if ( factor.nAlignments[3]>factor.nAlignments[1]*biasRatio31 &&
				  factor.nAlignments[3]>factor.nAlignments[2]*biasRatio32 ) {
			// Spanning factor
			tmp1=firstAlignmentForCurrentFactor;
			tmp2=firstAlignmentForNextFactor;
			firstAlignmentForCurrentFactor=firstAlignmentForPreviousFactor;
			out=keepLeftBoundaryCore(factor,previousFactor,factorID,keepLeftSplit,minCoverageRatio);
			firstAlignmentForCurrentFactor=tmp1;
			firstAlignmentForNextFactor=tmp2;
if (IO.SHOW_STD_ERR) IO.printErr("      keepLeftBoundary> case 4");
			return out;
		}
		else if (factor.rightScore>factor.leftScore*biasRatioRightLeft) {
			// Right factor
if (IO.SHOW_STD_ERR) IO.printErr("      keepLeftBoundary> case 5");
			return keepLeftSplit!=-1;
		}
		else if (factor.leftScore>factor.rightScore*biasRatioLeftRight) {
			// Left factor
			tmp1=firstAlignmentForCurrentFactor;
			tmp2=firstAlignmentForNextFactor;
			firstAlignmentForCurrentFactor=firstAlignmentForPreviousFactor;
			out=keepLeftBoundaryCore(factor,previousFactor,factorID,keepLeftSplit,minCoverageRatio);
			firstAlignmentForCurrentFactor=tmp1;
			firstAlignmentForNextFactor=tmp2;
if (IO.SHOW_STD_ERR) IO.printErr("      keepLeftBoundary> case 6");
			return out;
		}
		else {
			// Factor without any strong bias
			if (factor.periodicCoverage>0 || factor.lastDenseSubstring>0 || factor.inDenseSubstringOfSubstringType) {
				// Factor with no bias, inside a complex region. Merged with its
				// neighbors, unless something prevents it.
				tmp1=firstAlignmentForCurrentFactor;
				tmp2=firstAlignmentForNextFactor;
				firstAlignmentForCurrentFactor=firstAlignmentForPreviousFactor;
				out=keepLeftBoundaryCore(factor,previousFactor,factorID,keepLeftSplit,minCoverageRatio);
				firstAlignmentForCurrentFactor=tmp1;
				firstAlignmentForNextFactor=tmp2;
if (IO.SHOW_STD_ERR) IO.printErr("      keepLeftBoundary> case 7");
				return out;
			}
			else  {
				// Factor with no bias, outside a complex region. Votes zero for the left
				// boundary.
if (IO.SHOW_STD_ERR) IO.printErr("      keepLeftBoundary> case 8");
				return keepLeftSplit!=-1;
			}
		}
	}


	/**
	 * Decides whether to keep the left boundary of a factor, using the following:
	 * - coverage;
	 * - dense substrings;
	 * - period;
	 * - periodic substring intervals;
	 * - alignment intervals.
	 */
	private static final boolean keepLeftBoundaryCore(Factor factor, Factor previousFactor, int factorID, int keepLeftSplit, double minCoverageRatio) {
		return ( keepLeftSplit==0 &&
		         ( previousFactor.getCorrectedCoverage()*minCoverageRatio<factor.getCorrectedCoverage() ||
				   previousFactor.getCorrectedCoverage()>factor.getCorrectedCoverage()*minCoverageRatio ||
				   !previousFactor.sameDenseSubstrings(factor) ||
				   (previousFactor.periodicCoverage>0)^(factor.periodicCoverage>0) ||
				   (factor.periodicCoverage>0 && !previousFactor.samePeriodicSubstringIntervals(factor)) ||
				   !previousFactor.sameAlignmentIntervals(factor)
				 )
			   );
	}


	/**
	 * Decides whether to keep the right boundary of a factor, using maximality,
	 * type of dense substrings, maximal alignment statistics, alignment statistics,
	 * coverage, period, and internal symmetries, of the factor and of the factor to its
	 * right.
	 *
	 * Remark: we could merge factors whose length is smaller than a minimum factor
	 * length. For full generality, the current procedure does not do that.
	 *
	 * Remark: a factor without a strong bias in its covering alignments could also come
	 * from undersplitting. This procedure assumes that just oversplitting happened
	 * upstream.
	 *
	 * @param minRatio at least one;
	 * @return a vote: -1=remove the boundary; 1=keep the boundary; 0=no strong signal.
	 */
	private static final int keepRightBoundary(Factor factor, Factor nextFactor, int factorID, int lastFactorID) {
		int i;
		int observed, expected;

		if (factorID==lastFactorID) return 1;
		if (factor.lastPosition<nextFactor.firstPosition-1) return 1;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr(factor.lastPosition+" inLengthBiasInterval? "+DenseSubstrings.inLengthBiasInterval(factor.lastPosition,maxLengthBoundaryDistance));
		if (DenseSubstrings.inLengthBiasInterval(factor.lastPosition,maxLengthBoundaryDistance)) return 1;
	 	if (factor.isRightMaximal()) {
	 		// Right-maximal factor
			if (!factor.inDenseSubstringOfPrefixType && !factor.inDenseSubstringOfSubstringType && factor.periodicCoverage==0) {
				// Right-maximality is enough to keep the right boundary
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("      keepRightBoundary> case 1");
				return 1;
			}
			else if (factor.inDenseSubstringOfSubstringType) {
				// Right-maximality is not enough: checking also the number of
				// right-maximal alignments at the boundary.
				if (factor.isRightBracketed) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("      keepRightBoundary> case 2");
					return 1;
				}
			}
			else {
				// In all other cases we fall back to the standard merging below, since
				// we can separate skeleton alignments from all other alignments.
			}
		}
		// The factor is not right-maximal, or right-maximality is not enough to decide
		// the boundary.
		if ( factor.nAlignments[0]>factor.nAlignments[1]*biasRatio01 &&
			 factor.nAlignments[0]>factor.nAlignments[2]*biasRatio02 &&
			 factor.nAlignments[0]>factor.nAlignments[3]*biasRatio03 ) {
			// Factor with most alignments inside
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("      keepRightBoundary> case 3");
			return 1;
		}
		else if ( factor.nAlignments[3]>factor.nAlignments[1]*biasRatio31 &&
				  factor.nAlignments[3]>factor.nAlignments[2]*biasRatio32 ) {
			// Spanning factor
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("      keepRightBoundary> case 4");
			return keepRightBoundaryCore(factor,nextFactor,factorID,minCoverageRatio)?1:-1;
		}
		else if (factor.rightScore>factor.leftScore*biasRatioRightLeft) {
			// Right factor
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("      keepRightBoundary> case 5");
			return keepRightBoundaryCore(factor,nextFactor,factorID,minCoverageRatio)?1:-1;
		}
		else if (factor.leftScore>factor.rightScore*biasRatioLeftRight) {
			// Left factor
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("      keepRightBoundary> case 6");
			return 0;
		}
		else {
			// Factor without any strong bias
			if (factor.periodicCoverage>0 || factor.lastDenseSubstring>0 || factor.inDenseSubstringOfSubstringType) {
				// Factor with no bias, inside a complex region. Merged with its
				// neighbors, unless something prevents it.
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("      keepRightBoundary> case 7");
				return keepRightBoundaryCore(factor,nextFactor,factorID,minCoverageRatio)?1:-1;
			}
			else  {
				// Factors with no bias, outside complex regions. Votes zero for the right
				// boundary.
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("      keepRightBoundary> case 8");
				return 0;
			}
		}
	}


	/**
	 * Decides whether to keep the right boundary of a factor, using the following:
	 * - coverage;
	 * - dense substrings;
	 * - period;
	 * - periodic substring intervals;
	 * - alignment intervals.
	 */
	private static final boolean keepRightBoundaryCore(Factor factor, Factor nextFactor, int factorID, double minCoverageRatio) {
		return ( nextFactor.getCorrectedCoverage()*minCoverageRatio<factor.getCorrectedCoverage() ||
				 nextFactor.getCorrectedCoverage()>factor.getCorrectedCoverage()*minCoverageRatio ||
				 !factor.sameDenseSubstrings(nextFactor) ||
				 (nextFactor.periodicCoverage>0)^(factor.periodicCoverage>0) ||
				 (factor.periodicCoverage>0 && !factor.samePeriodicSubstringIntervals(nextFactor)) ||
				 !factor.sameAlignmentIntervals(nextFactor)
			   );
	}


	/**
	 * @return TRUE iff $factors[factor]$ is the leftmost in at least one of the dense
	 * substrings that cover it.
	 */
	private static final boolean leftmostInDenseSubstring(int factor) {
		if (factor==0) return true;
		for (int i=0; i<=factors[factor].lastDenseSubstring; i++) {
			if (factors[factor-1].inDenseSubstring(factors[factor].denseSubstrings[i].id)==-1) return true;
		}
		return false;
	}


	/**
	 * @return TRUE iff $factors[factor]$ is the rightmost in at least one of the dense
	 * substrings that cover it.
	 */
	private static final boolean rightmostInDenseSubstring(int factor) {
		if (factor==lastFactor) return true;
		for (int i=0; i<=factors[factor].lastDenseSubstring; i++) {
			if (factors[factor+1].inDenseSubstring(factors[factor].denseSubstrings[i].id)==-1) return true;
		}
		return false;
	}


	/**
	 * Merges two adjacent factors $X$ and $Y$ iff both $X.keepRightBoundary$ and
	 * $Y.keepLeftBoundary$ are false, i.e. iff no interval induced by dense substrings,
	 * periodic substrings, or alignments, starts or ends at the boundary between $X$ and
	 * $Y$.
	 */
	private static final int mergeFactors2() {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i, out;
		int firstFactorInRun;

		for (i=1; i<=lastFactor; i++) factors[i].shouldBeOutput=true;
		out=0;
		firstFactorInRun=0;
		for (i=1; i<=lastFactor; i++) {

if (IO.SHOW_STD_ERR) IO.printErr("mergeFactors2> factor "+i+"=["+factors[i].firstPosition+".."+factors[i].lastPosition+"]");
if (IO.SHOW_STD_ERR) IO.printErr("mergeFactors2> factor["+i+"].keepLeftBoundary="+factors[i].keepLeftBoundary+" factor["+(i-1)+"].keepRightBoundary="+(i>0 && factors[i-1].keepRightBoundary));

			if (factors[i].keepLeftBoundary || factors[i-1].keepRightBoundary) firstFactorInRun=i;
			else {
				factors[i].shouldBeOutput=false;
				out++;
				factors[firstFactorInRun].lastPosition=factors[i].lastPosition;
				factors[firstFactorInRun].length=factors[firstFactorInRun].lastPosition-factors[firstFactorInRun].firstPosition+1;
				factors[firstFactorInRun].isRightBracketed=factors[i].isRightBracketed;
				factors[firstFactorInRun].nOpenRight=factors[i].nOpenRight;
				factors[firstFactorInRun].nClosedRight=factors[i].nClosedRight;
				factors[firstFactorInRun].nUnknownRight=factors[i].nUnknownRight;
				if (factors[i].firstPosition<=factors[firstFactorInRun].firstPosition+IDENTITY_THRESHOLD) {
					factors[firstFactorInRun].nOpenLeft+=factors[i].nOpenLeft;
					factors[firstFactorInRun].nClosedLeft+=factors[i].nClosedLeft;
					factors[firstFactorInRun].nUnknownLeft+=factors[i].nUnknownLeft;
				}
			}
		}
		return out;
	}


	/**
	 * Removes from $factors$ all elements for which $shouldBeOutput=false$.
	 */
	private static final void cleanFactors(int nDiscardedFactors) {
		Factor.order=Factor.SHOULDBEOUTPUT_FIRSTPOSITION;
		if (lastFactor>0) Arrays.sort(factors,0,lastFactor+1);
		lastFactor=lastFactor-nDiscardedFactors;
	}


	/**
	 * Assigns to each factor the set of all dense substrings, of any type, that contain
	 * it (if $containing=true$), or that either contain it or are contained in it (if
	 * $containing=false$).
	 *
	 * Remark: before factor merging, a factor can only be contained inside dense
	 * substrings, and it cannot contain any dense substring. After merging, it can
	 * contain one or more dense substrings
	 *
	 * Remark: the procedure assumes $factors$ to be sorted by $firstPosition$.
	 *
	 * Remark: assume that a factor coincides with a dense substring that replicates by
	 * taking prefixes of itself. Recall that we decide whether an occurrence is full or
	 * partial using a threshold. It could thus happen that the factor has no prefix
	 * coverage and it has just full coverage, so the fact that it replicates by taking
	 * prefixes of itself could be lost. This is the reason why we don't discard dense
	 * substrings that intersect just one factor.
	 */
	private static final void assignDenseSubstrings(boolean containing) {
		int i, j;
		int firstJForNextSubstring, distance;

		// Ensuring the necessary order in $DenseSubstrings.substrings$ and in
		// $ReadA.sortedAlignments$
		if (DenseSubstring.order!=DenseSubstring.STARTA) {
			DenseSubstring.order=DenseSubstring.STARTA;
			if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
		}
		// Assigning IDs to dense substrings. These are necessary for the following
		// steps, and might differ from the final assignment of IDs reported in output.
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].id=i;
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("assignDenseSubstrings> list of dense substrings: ("+containing+")");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=DenseSubstrings.lastSubstring; i++) IO.printErr("id="+DenseSubstrings.substrings[i].id+" ["+DenseSubstrings.substrings[i].startA+".."+DenseSubstrings.substrings[i].endA+"]"); }


		// Assigning dense substrings to factors
		for (i=0; i<=lastFactor; i++) {
			factors[i].inDenseSubstring=false;
			factors[i].inDenseSubstringOfPrefixType=false;
			factors[i].inDenseSubstringOfSuffixType=false;
			factors[i].inDenseSubstringOfSubstringType=false;
			factors[i].lastDenseSubstring=-1;
		}
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) {
			DenseSubstrings.substrings[i].firstFactor=-1;
			DenseSubstrings.substrings[i].lastFactor=-1;
		}
		i=0; j=0; firstJForNextSubstring=-1;
		while (i<=DenseSubstrings.lastSubstring) {
			if (j>lastFactor || factors[j].firstPosition>=DenseSubstrings.substrings[i].endA) {
				i++;
				if (firstJForNextSubstring!=-1) j=firstJForNextSubstring;
				firstJForNextSubstring=-1;
				continue;
			}
			if (factors[j].lastPosition<=DenseSubstrings.substrings[i].startA) {
				j++;
				continue;
			}
			if (i<DenseSubstrings.lastSubstring && firstJForNextSubstring==-1 && factors[j].lastPosition>=DenseSubstrings.substrings[i+1].startA) firstJForNextSubstring=j;
			if ( Intervals.isApproximatelyContained_lowQuality(factors[j].firstPosition,factors[j].lastPosition,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,ReadA.id) ||
				 Intervals.areApproximatelyIdentical_lowQuality(factors[j].firstPosition,factors[j].lastPosition,DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,ReadA.id)
			   ) {
				factors[j].inDenseSubstring=true;
				if (DenseSubstrings.substrings[i].prefixReplication) factors[j].inDenseSubstringOfPrefixType=true;
				if (DenseSubstrings.substrings[i].suffixReplication) factors[j].inDenseSubstringOfSuffixType=true;
				if (DenseSubstrings.substrings[i].substringReplication) factors[j].inDenseSubstringOfSubstringType=true;
				factors[j].lastDenseSubstring++;
				factors[j].ensureSpace_denseSubstrings(factors[j].lastDenseSubstring+1);
				factors[j].denseSubstrings[factors[j].lastDenseSubstring]=DenseSubstrings.substrings[i];
				distance=DenseSubstrings.substrings[i].startA-factors[j].firstPosition;
				if ( (distance>=0 && distance<=Factor.ANCHORING_DISTANCE && distance<=(factors[j].length>>1)) ||
				     (distance<0 && -distance<=Factor.ANCHORING_DISTANCE && (j==0 || -distance<=(factors[j-1].length>>1)))
				   ) {
					factors[j].keepLeftBoundary=true;
					if (DenseSubstrings.substrings[i].firstFactor==-1 || j<DenseSubstrings.substrings[i].firstFactor) DenseSubstrings.substrings[i].firstFactor=j;
				}
				distance=DenseSubstrings.substrings[i].endA-factors[j].lastPosition;
				if ( (distance<0 && -distance<=Factor.ANCHORING_DISTANCE && -distance<=(factors[j].length>>1)) ||
				     (distance>=0 && distance<=Factor.ANCHORING_DISTANCE && (j==lastFactor || distance<=(factors[j+1].length>>1)))
				   ) {
					factors[j].keepRightBoundary=true;
					if (DenseSubstrings.substrings[i].lastFactor==-1 || j>DenseSubstrings.substrings[i].lastFactor) DenseSubstrings.substrings[i].lastFactor=j;
				}
			}
			else if (!containing && Intervals.isApproximatelyContained_lowQuality(DenseSubstrings.substrings[i].startA,DenseSubstrings.substrings[i].endA,factors[j].firstPosition,factors[j].lastPosition,ReadA.id)) {
				factors[j].lastDenseSubstring++;
				factors[j].ensureSpace_denseSubstrings(factors[j].lastDenseSubstring+1);
				factors[j].denseSubstrings[factors[j].lastDenseSubstring]=DenseSubstrings.substrings[i];
				distance=DenseSubstrings.substrings[i].startA-factors[j].firstPosition;
				if (distance<0) distance=-distance;
				if (distance<=Factor.ANCHORING_DISTANCE) {
					factors[j].keepLeftBoundary=true;
					if (DenseSubstrings.substrings[i].firstFactor==-1) DenseSubstrings.substrings[i].firstFactor=j;
				}
				distance=DenseSubstrings.substrings[i].endA-factors[j].lastPosition;
				if (distance<0) distance=-distance;
				if (distance<=Factor.ANCHORING_DISTANCE) {
					factors[j].keepRightBoundary=true;
					if (DenseSubstrings.substrings[i].lastFactor==-1) DenseSubstrings.substrings[i].lastFactor=j;
				}
			}
			j++;
		}

		// Sorting by ID the dense substrings of each factor
		DenseSubstring.order=DenseSubstring.ID;
		for (i=0; i<=lastFactor; i++) {
			if (factors[i].lastDenseSubstring>0) Arrays.sort(factors[i].denseSubstrings,0,factors[i].lastDenseSubstring+1);
		}
	}


	/**
	 * Assigns to each factor the set of all periodic substring intervals that contain
	 * it (if $containing=true$), or that either contain it or are contained in it (if
	 * $containing=false$). In the second case, intervals that are contained in a factor 
	 * are grouped in a compact block at the end of array $PeriodicSubstrings.intervals$, 
	 * and they are marked as discarded.
	 *
	 * Remark: before factor merging, a factor can only be contained inside periodic
	 * substring intervals, and it cannot contain any periodic substring interval. After
	 * merging, it can contain one or more periodic substring intervals.
	 *
	 * Remark: the procedure assumes $factors$ to be sorted by $firstPosition$.
	 */
	private static final void assignPeriodicSubstringIntervals(boolean containing) {
		int i, j, k, firstJForNextInterval;
		int distance;
		PeriodicSubstringInterval tmp;

		// Ensuring the necessary order in $PeriodicSubstrings.intervals$
		if (PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (PeriodicSubstrings.lastLongPeriodInterval>0) Arrays.sort(PeriodicSubstrings.longPeriodIntervals,0,PeriodicSubstrings.lastLongPeriodInterval+1);
		}
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		}
		// Assigning IDs to periodic substring intervals. These are necessary for the
		// following steps, and might differ from the final assignment of IDs reported in
		// output.
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) PeriodicSubstrings.intervals[i].id=i;

		// Assigning intervals to factors
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			PeriodicSubstrings.intervals[i].discarded=false;
			PeriodicSubstrings.intervals[i].nFactors=0;
			PeriodicSubstrings.intervals[i].firstFactor=0;
			PeriodicSubstrings.intervals[i].lastFactor=0;
		}
		for (i=0; i<=lastFactor; i++) {
			factors[i].lastPeriodicSubstringInterval=-1;
			factors[i].inPeriodicSubstringInterval=false;
		}
		i=0; j=0; firstJForNextInterval=-1;
		while (i<=PeriodicSubstrings.lastInterval) {
			if (j>lastFactor || factors[j].firstPosition>=PeriodicSubstrings.intervals[i].lastPosition) {
				i++;
				if (firstJForNextInterval!=-1) j=firstJForNextInterval;
				firstJForNextInterval=-1;
				continue;
			}
			if (factors[j].lastPosition<=PeriodicSubstrings.intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (i<PeriodicSubstrings.lastInterval && firstJForNextInterval==-1 && factors[j].lastPosition>=PeriodicSubstrings.intervals[i+1].firstPosition) firstJForNextInterval=j;
			if ( !PeriodicSubstrings.intervals[i].hasLongPeriod &&
				 ( Intervals.isApproximatelyContained_lowQuality(factors[j].firstPosition,factors[j].lastPosition,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id) ||
				   Intervals.areApproximatelyIdentical_lowQuality(factors[j].firstPosition,factors[j].lastPosition,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id)
				 )
			   ) factors[j].inPeriodicSubstringInterval=true;
			if ( ( containing &&
			       ( Intervals.isApproximatelyContained_lowQuality(factors[j].firstPosition,factors[j].lastPosition,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id) ||
			     	 Intervals.areApproximatelyIdentical_lowQuality(factors[j].firstPosition,factors[j].lastPosition,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id)
			       )
			     ) ||
			     ( !containing &&
			       ( Intervals.isApproximatelyContained_lowQuality(factors[j].firstPosition,factors[j].lastPosition,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id) ||
			     	 Intervals.areApproximatelyIdentical_lowQuality(factors[j].firstPosition,factors[j].lastPosition,PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,ReadA.id) ||
			         Intervals.isApproximatelyContained_lowQuality(PeriodicSubstrings.intervals[i].firstPosition,PeriodicSubstrings.intervals[i].lastPosition,factors[j].firstPosition,factors[j].lastPosition,ReadA.id)
			       )
			     )
			   ) {
				factors[j].lastPeriodicSubstringInterval++;
				factors[j].ensureSpace_periodicSubstringIntervals(factors[j].lastPeriodicSubstringInterval+1);
				factors[j].periodicSubstringIntervals[factors[j].lastPeriodicSubstringInterval]=PeriodicSubstrings.intervals[i];
				PeriodicSubstrings.intervals[i].nFactors++;
				distance=PeriodicSubstrings.intervals[i].firstPosition-factors[j].firstPosition;
				if ( (distance>=0 && distance<=Factor.ANCHORING_DISTANCE && distance<=(factors[j].length>>1)) ||
				     (distance<0 && -distance<=Factor.ANCHORING_DISTANCE && (j==0 || -distance<=(factors[j-1].length>>1)))
				   ) {
					factors[j].keepLeftBoundary=true;
					if (PeriodicSubstrings.intervals[i].firstFactor==-1 || j<PeriodicSubstrings.intervals[i].firstFactor) PeriodicSubstrings.intervals[i].firstFactor=j;
				}
				distance=PeriodicSubstrings.intervals[i].lastPosition-factors[j].lastPosition;
				if ( (distance<0 && -distance<=Factor.ANCHORING_DISTANCE && -distance<=(factors[j].length>>1)) ||
				     (distance>=0 && distance<=Factor.ANCHORING_DISTANCE && (j==lastFactor || distance<=(factors[j+1].length>>1)))
				   ) {
					factors[j].keepRightBoundary=true;
					if (PeriodicSubstrings.intervals[i].lastFactor==-1 || j>PeriodicSubstrings.intervals[i].lastFactor) PeriodicSubstrings.intervals[i].lastFactor=j;
				}
				j++;
				continue;
			}
			j++;
		}

		// Sorting the intervals of each factor by their ID
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.ID;
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.ID;
		for (i=0; i<=lastFactor; i++) {
			if (factors[i].lastPeriodicSubstringInterval>0) Arrays.sort(factors[i].periodicSubstringIntervals,0,factors[i].lastPeriodicSubstringInterval+1);
		}
		PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
		PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;

		// Compacting intervals that are contained in a factor
		if (!containing) {
			for (i=0; i<=lastFactor; i++) {
				k=-1;
				for (j=0; j<=factors[i].lastPeriodicSubstringInterval; j++) {
					if ( factors[i].periodicSubstringIntervals[j].nFactors==1 &&
						 Intervals.isApproximatelyContained_lowQuality(factors[i].periodicSubstringIntervals[j].firstPosition,factors[i].periodicSubstringIntervals[j].lastPosition,factors[i].firstPosition,factors[i].lastPosition,ReadA.id) &&
						 !Intervals.areApproximatelyIdentical_lowQuality(factors[i].periodicSubstringIntervals[j].firstPosition,factors[i].periodicSubstringIntervals[j].lastPosition,factors[i].firstPosition,factors[i].lastPosition,ReadA.id)
					   ) {
						factors[i].periodicSubstringIntervals[j].discarded=true;
					   	continue;
					}
					k++;
					if (k!=j) {
						tmp=factors[i].periodicSubstringIntervals[k];
						factors[i].periodicSubstringIntervals[k]=factors[i].periodicSubstringIntervals[j];
						factors[i].periodicSubstringIntervals[j]=tmp;
					}
				}
				factors[i].lastPeriodicSubstringInterval=k;
			}
			j=-1;
			for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
				if (PeriodicSubstrings.intervals[i].discarded) continue;
				j++;
				if (j!=i) {
					tmp=PeriodicSubstrings.intervals[j];
					PeriodicSubstrings.intervals[j]=PeriodicSubstrings.intervals[i];
					PeriodicSubstrings.intervals[i]=tmp;
				}
			}
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.UNSORTED;
		}
	}


	/**
	 * Assigns to each factor all alignment intervals that approximately contain it. Then,
	 * it marks as discarded every interval such that:
	 *
	 * 1. it contains no factor;
	 *
	 * 2. it is approximately contained in a dense substring of substring type, and:
	 *
	 * 2.1 its first (respectively, last) position is not close to the first 
	 * (respectively, last) position of a factor; or
	 * 2.2 the number of alignments that start (respectively, end) at the first 
	 * (respectively, last) position of its first (respectively, last) factor is not 
	 * significantly larger than expected, in any dense substring that contains the 
	 * factor;
	 *
	 * 3. it is approximately contained in a short-period interval, and:
	 *
	 * 3.1 its first (respectively, last) position is not close to the first 
	 * (respectively, last) position of a factor.
	 *
	 * Alignment intervals that are marked as discarded are grouped into a compact 
	 * interval at the end of array $ReadA.intervals$, and they are removed from the array
	 * $alignmentIntervals$ of each factor. Iff $discard=true$, they are also pushed out 
	 * of visibility in array $AlignmentIntervals.intervals$ and removed from alignment
	 * pointers.
	 * 
	 * Remark: the procedure does not discard a maximal alignment interval inside a short-
	 * period interval, if it has high error rate and enough alignments, even though it 
	 * does not contain any factor. This is because the corresponding peaks of maximal 
	 * events might have been too faint to be detected (if the number of alignments in the
	 * interval was not large), yet the error rate signal might be strong enough to detect
	 * a different repeat.
	 *
	 * Remark: the procedure keeps alignment intervals that cover exactly one factor.
	 *
	 * Remark: the procedure assumes $factors$ to be sorted by $firstPosition$.
	 */
	private static final void assignAlignmentIntervals(boolean discard) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int SHORT_FACTOR_LENGTH = IO.quantum<<1;  // Arbitrary
		final int MIN_OCCURRENCES_NONMAXIMAL = IO.coverage*40;  // Arbitrary
		final double FACTORS_SURFACE_THRESHOLD = 0.5;  // Arbitrary
		final int MIN_ALIGNMENTS_IN_PERIODIC = IO.coverage*10;  // Arbitrary
		final double ERROR_RATE_RATIO = 2.0;  // Arbitrary
		boolean isLeftMaximal, isRightMaximal;
		int i, j, k;
		int id, firstJForNextInterval, distance, type, lastInterval, previousThreshold;
		Factor factorLeft, factorRight;
		AlignmentInterval interval;
		DenseSubstring substring;
		int[] tmpFactors;

		// Ensuring the necessary order in $AlignmentIntervals.intervals$
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (AlignmentIntervals.lastInterval>0) Arrays.sort(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1);
		}

		// Assigning IDs to alignment intervals. These are necessary for the following
		// steps, and might differ from the final assignment of IDs reported in output.
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) AlignmentIntervals.intervals[i].id=i;

		// Assigning alignment intervals to factors, and vice versa.
		previousThreshold=Intervals.absoluteThreshold;
		Intervals.absoluteThreshold=Factor.ANCHORING_DISTANCE;
		for (i=0; i<=lastFactor; i++) factors[i].lastAlignmentInterval=-1;
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			AlignmentIntervals.intervals[i].cleanFactors();
			AlignmentIntervals.intervals[i].isAnchoredLeft=null;
			AlignmentIntervals.intervals[i].isAnchoredRight=null;
			AlignmentIntervals.intervals[i].containerFactor=null;
		}
		i=0; j=0; firstJForNextInterval=-1;
		while (i<=AlignmentIntervals.lastInterval) {
			if (j>lastFactor || factors[j].firstPosition>=AlignmentIntervals.intervals[i].lastPosition) {
				i++;
				if (firstJForNextInterval!=-1) j=firstJForNextInterval;
				firstJForNextInterval=-1;
				continue;
			}
			if (factors[j].lastPosition<=AlignmentIntervals.intervals[i].firstPosition) {
				j++;
				continue;
			}
			if (i<AlignmentIntervals.lastInterval && firstJForNextInterval==-1 && factors[j].lastPosition>=AlignmentIntervals.intervals[i+1].firstPosition) firstJForNextInterval=j;
			if ( Intervals.isApproximatelyContained_lowQuality(factors[j].firstPosition,factors[j].lastPosition,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,ReadA.id) ||
			     Intervals.areApproximatelyIdentical_lowQuality(factors[j].firstPosition,factors[j].lastPosition,AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,ReadA.id)
			   ) {
				factors[j].lastAlignmentInterval++;
				factors[j].ensureSpace_alignmentIntervals(factors[j].lastAlignmentInterval+1);
				factors[j].alignmentIntervals[factors[j].lastAlignmentInterval]=AlignmentIntervals.intervals[i];
				AlignmentIntervals.intervals[i].lastFactor++;
				if (AlignmentIntervals.intervals[i].factors.length<=AlignmentIntervals.intervals[i].lastFactor) {
					tmpFactors = new int[AlignmentIntervals.intervals[i].factors.length+AlignmentInterval.FACTORS_GROWTH_RATE];
					System.arraycopy(AlignmentIntervals.intervals[i].factors,0,tmpFactors,0,AlignmentIntervals.intervals[i].factors.length);
					AlignmentIntervals.intervals[i].factors=tmpFactors;
				}
				AlignmentIntervals.intervals[i].factors[AlignmentIntervals.intervals[i].lastFactor]=j;
				if (AlignmentIntervals.intervals[i].isAnchoredLeft==null) {
					distance=Math.abs(AlignmentIntervals.intervals[i].firstPosition,factors[j].firstPosition);
					if (distance<=Factor.ANCHORING_DISTANCE) AlignmentIntervals.intervals[i].isAnchoredLeft=factors[j];
				}
				distance=Math.abs(AlignmentIntervals.intervals[i].lastPosition,factors[j].lastPosition);
				if (distance<=Factor.ANCHORING_DISTANCE) AlignmentIntervals.intervals[i].isAnchoredRight=factors[j];
			}
			if ( Intervals.isApproximatelyContained_lowQuality(AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,factors[j].firstPosition,factors[j].lastPosition,ReadA.id) &&
				 !Intervals.areApproximatelyIdentical_lowQuality(AlignmentIntervals.intervals[i].firstPosition,AlignmentIntervals.intervals[i].lastPosition,factors[j].firstPosition,factors[j].lastPosition,ReadA.id)
			   ) {
				if ( AlignmentIntervals.intervals[i].isLeftMaximal && !AlignmentIntervals.intervals[i].isRightMaximal && 
					 Math.abs(AlignmentIntervals.intervals[i].firstPosition,factors[j].firstPosition)<=Factor.ANCHORING_DISTANCE
				   ) AlignmentIntervals.intervals[i].containerFactor=factors[j];
				else if ( AlignmentIntervals.intervals[i].isRightMaximal && !AlignmentIntervals.intervals[i].isLeftMaximal && 
					      Math.abs(AlignmentIntervals.intervals[i].lastPosition,factors[j].lastPosition)<=Factor.ANCHORING_DISTANCE
				        ) AlignmentIntervals.intervals[i].containerFactor=factors[j];
				else if ( !AlignmentIntervals.intervals[i].isLeftMaximal && !AlignmentIntervals.intervals[i].isRightMaximal && 
					      AlignmentIntervals.intervals[i].nMergedIntervals>=MIN_OCCURRENCES_NONMAXIMAL
				        ) AlignmentIntervals.intervals[i].containerFactor=factors[j];
			}
			j++;
		}
		Intervals.absoluteThreshold=previousThreshold;

		// Ensuring the necessary order in $AlignmentIntervals.intervals$
		AlignmentInterval.order=AlignmentInterval.FACTORS;
		if (AlignmentIntervals.lastInterval>0) Arrays.sort(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1);



if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("ALIGNMENT INTERVALS TO FACTORS:");
	for (int x=0; x<=AlignmentIntervals.lastInterval; x++) {
		System.err.println("alignment interval: "+AlignmentIntervals.intervals[x]);
		System.err.println("    isAnchoredLeft: "+AlignmentIntervals.intervals[x].isAnchoredLeft);
		System.err.println("   isAnchoredRight: "+AlignmentIntervals.intervals[x].isAnchoredRight);
		System.err.println("   containerFactor: "+AlignmentIntervals.intervals[x].containerFactor);
		System.err.print("           factors: ");
		for (int y=0; y<=AlignmentIntervals.intervals[x].lastFactor; y++) System.err.print(factors[AlignmentIntervals.intervals[x].factors[y]]+", ");
		System.err.println();
		System.err.println("inDenseSubstringOfSubstringType: "+(AlignmentIntervals.intervals[x].inDenseSubstringOfSubstringType!=null?AlignmentIntervals.intervals[x].inDenseSubstringOfSubstringType.hashCode():"null"));
		System.err.println("inPeriodicSubstringInterval: "+(AlignmentIntervals.intervals[x].inPeriodicSubstringInterval!=null?AlignmentIntervals.intervals[x].inPeriodicSubstringInterval.hashCode():"null"));
		System.err.println(); System.err.println();
	}
}


		// Assigning to each alignment interval its shortest containing short-period
		// interval, if any.
		if (AlignmentInterval.order!=AlignmentInterval.FIRSTPOSITION) {
			AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
			if (AlignmentIntervals.lastInterval>0) Arrays.sort(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1);
		}
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION || PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		}
		AlignmentIntervals.markIntervalsInPeriodic(0,AlignmentIntervals.lastInterval,true);
		assignAlignmentIntervals_getErrorRate();

		// Discarding alignment intervals
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			AlignmentIntervals.intervals[i].discarded=false;
			AlignmentIntervals.intervals[i].representative=null;
		}
		j=-1;
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			if (AlignmentIntervals.intervals[i].isAnchoredLeft!=null) {
				AlignmentIntervals.intervals[i].isAnchoredLeft.keepLeftBoundary=true;
				if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setting keepLeftBoundary=true for factor ["+AlignmentIntervals.intervals[i].isAnchoredLeft.firstPosition+".."+AlignmentIntervals.intervals[i].isAnchoredLeft.lastPosition+"] = "+AlignmentIntervals.intervals[i].isAnchoredLeft);
			}
			if (AlignmentIntervals.intervals[i].isAnchoredRight!=null) {
				AlignmentIntervals.intervals[i].isAnchoredRight.keepRightBoundary=true;
				if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setting keepRightBoundary=true for factor ["+AlignmentIntervals.intervals[i].isAnchoredRight.firstPosition+".."+AlignmentIntervals.intervals[i].isAnchoredRight.lastPosition+"] = "+AlignmentIntervals.intervals[i].isAnchoredRight);
			}
			if (AlignmentIntervals.intervals[i].containerFactor!=null) {
				if (AlignmentIntervals.intervals[i].isLeftMaximal) {
					AlignmentIntervals.intervals[i].containerFactor.keepLeftBoundary=true;
					if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setting keepLeftBoundary=true for factor ["+AlignmentIntervals.intervals[i].containerFactor.firstPosition+".."+AlignmentIntervals.intervals[i].containerFactor.lastPosition+"] = "+AlignmentIntervals.intervals[i].containerFactor);
				}
				if (AlignmentIntervals.intervals[i].isRightMaximal) {
					AlignmentIntervals.intervals[i].containerFactor.keepRightBoundary=true;
					if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setting keepRightBoundary=true for factor ["+AlignmentIntervals.intervals[i].containerFactor.firstPosition+".."+AlignmentIntervals.intervals[i].containerFactor.lastPosition+"] = "+AlignmentIntervals.intervals[i].containerFactor);
				}
			}			
			if ( ( AlignmentIntervals.intervals[i].lastFactor<0 &&
				   ( AlignmentIntervals.intervals[i].inDenseSubstringOfSubstringType!=null ||
				     ( AlignmentIntervals.intervals[i].inPeriodicSubstringInterval!=null &&
					   ( (!AlignmentIntervals.intervals[i].isLeftMaximal && !AlignmentIntervals.intervals[i].isRightMaximal) ||
						 AlignmentIntervals.intervals[i].errorRate/AlignmentIntervals.intervals[i].nMergedIntervals<ERROR_RATE_RATIO*(AlignmentIntervals.intervals[i].inPeriodicSubstringInterval.errorRate/AlignmentIntervals.intervals[i].inPeriodicSubstringInterval.nImpliedAlignments) ||
					     AlignmentIntervals.intervals[i].nMergedIntervals<MIN_ALIGNMENTS_IN_PERIODIC
					   )
					 ) ||
					 ( AlignmentIntervals.intervals[i].inPeriodicSubstringInterval==null &&
					   AlignmentIntervals.intervals[i].containerFactor==null		 
					 )
				   )
				 ) ||
				 ( AlignmentIntervals.intervals[i].lastFactor>=0 &&
				   getFactorsSurface(AlignmentIntervals.intervals[i])<AlignmentIntervals.intervals[i].length()*FACTORS_SURFACE_THRESHOLD &&
				   !factorBoundaryOrLowQuality(AlignmentIntervals.intervals[i])
				 )
			   ) {
				AlignmentIntervals.intervals[i].discarded=true;
				if (IO.SHOW_STD_ERR_PRIME) {
					IO.printErr("(0)discarded interval "+AlignmentIntervals.intervals[i]);
					IO.printErr("container factor: "+AlignmentIntervals.intervals[i].containerFactor);
				}
				continue;
			}
			if (AlignmentIntervals.intervals[i].inDenseSubstringOfSubstringType!=null) {
				if (IO.SHOW_STD_ERR_PRIME) IO.printErr("in dense substring of substring type interval "+AlignmentIntervals.intervals[i]);
				factorLeft=AlignmentIntervals.intervals[i].isAnchoredLeft;
				factorRight=AlignmentIntervals.intervals[i].isAnchoredRight;
				if ( factorLeft==null ||
				     ( factorLeft.inDenseSubstringOfSubstringType && 
					   ( ( Reads.isLeftRightMaximal(factorLeft.firstPosition,ReadA.id,true) && 
					       !factorLeft.isLeftMaximal() && !isNextFactorMaximal(factorLeft,true,SHORT_FACTOR_LENGTH)
					     ) ||
						 (!factorLeft.checkSplits(true) && !checkSplitsOfNextFactor(factorLeft,true,SHORT_FACTOR_LENGTH))
					   )
				     ) ||
				     factorRight==null ||
				     ( factorRight.inDenseSubstringOfSubstringType && 
					   ( ( Reads.isLeftRightMaximal(factorRight.lastPosition,ReadA.id,true) && 
					       !factorRight.isRightMaximal() && !isNextFactorMaximal(factorRight,false,SHORT_FACTOR_LENGTH)
					     ) ||
						 (!factorRight.checkSplits(false) && !checkSplitsOfNextFactor(factorRight,false,SHORT_FACTOR_LENGTH))
					   )
					 )
				   ) {
					AlignmentIntervals.intervals[i].discarded=true;
					if (IO.SHOW_STD_ERR_PRIME) {
						IO.printErr("(1)discarded interval "+AlignmentIntervals.intervals[i]+", isAnchoredLeft="+(AlignmentIntervals.intervals[i].isAnchoredLeft!=null)+" isAnchoredRight="+(AlignmentIntervals.intervals[i].isAnchoredRight!=null));
						if (factorLeft==null) IO.printErr("  factorLeft: null");
						else IO.printErr("  factorLeft: ["+factorLeft.firstPosition+".."+factorLeft.lastPosition+"] leftMaximal="+Reads.isLeftMaximal(factorLeft.firstPosition,ReadA.id,true)+" isLeftMaximal="+factorLeft.isLeftMaximal()+" checkSplits="+factorLeft.checkSplits(true));
						if (factorRight==null) IO.printErr("  factorRight: null");
						else IO.printErr("  factorRight: ["+factorRight.firstPosition+".."+factorRight.lastPosition+"] rightMaximal="+Reads.isRightMaximal(factorRight.lastPosition,ReadA.id,true)+" isRightMaximal="+factorRight.isRightMaximal()+" checkSplits="+factorRight.checkSplits(false));
					}
					continue;
				}
				if ( factorLeft.inDenseSubstringOfSubstringType && !factorLeft.isFirstFactor(true) && 
				     !isBracketed(factorLeft,AlignmentIntervals.intervals[i].firstPosition,true) &&
					 factorRight.inDenseSubstringOfSubstringType && !factorRight.isFirstFactor(false) &&
					 !isBracketed(factorRight,AlignmentIntervals.intervals[i].lastPosition,false)
				   ) {
					AlignmentIntervals.intervals[i].discarded=true;
					if (IO.SHOW_STD_ERR_PRIME) IO.printErr("(2)discarded interval "+AlignmentIntervals.intervals[i]+", left factor: ["+factorLeft.firstPosition+".."+factorLeft.lastPosition+"]");
					continue;
				}
			}
			else if (AlignmentIntervals.intervals[i].inPeriodicRange && AlignmentIntervals.intervals[i].lastFactor>=0) {
				if (IO.SHOW_STD_ERR_PRIME) IO.printErr("in short-period interval "+AlignmentIntervals.intervals[i]);
				factorLeft=AlignmentIntervals.intervals[i].isAnchoredLeft;
				factorRight=AlignmentIntervals.intervals[i].isAnchoredRight;
				if ( factorLeft==null ||
				     ( factorLeft.inPeriodicSubstringInterval && 
					   ( ( Reads.isLeftRightMaximal(factorLeft.firstPosition,ReadA.id,true) && 
					       !factorLeft.isLeftMaximal() && !isNextFactorMaximal(factorLeft,true,SHORT_FACTOR_LENGTH)
						 ) ||
						 (!factorLeft.checkSplits(true) && !checkSplitsOfNextFactor(factorLeft,true,SHORT_FACTOR_LENGTH))
					   )
				     ) ||
				     factorRight==null ||
				     ( factorRight.inPeriodicSubstringInterval && 
					   ( ( Reads.isLeftRightMaximal(factorRight.lastPosition,ReadA.id,true) && 
					       !factorRight.isRightMaximal() && !isNextFactorMaximal(factorRight,true,SHORT_FACTOR_LENGTH)
					     ) ||
						 (!factorRight.checkSplits(false) && !checkSplitsOfNextFactor(factorRight,false,SHORT_FACTOR_LENGTH))
					   )
					 )
				   ) {
					AlignmentIntervals.intervals[i].discarded=true;
					if (IO.SHOW_STD_ERR_PRIME) {
						IO.printErr("(3)discarded interval "+AlignmentIntervals.intervals[i]+", isAnchoredLeft="+(AlignmentIntervals.intervals[i].isAnchoredLeft!=null)+" isAnchoredRight="+(AlignmentIntervals.intervals[i].isAnchoredRight!=null));
						if (factorLeft==null) IO.printErr("  factorLeft: null");
						else IO.printErr("  factorLeft: ["+factorLeft.firstPosition+".."+factorLeft.lastPosition+"] leftMaximal="+Reads.isLeftMaximal(factorLeft.firstPosition,ReadA.id,true)+" isLeftMaximal="+factorLeft.isLeftMaximal());
						if (factorRight==null) IO.printErr("  factorRight: null");
						else IO.printErr("  factorRight: ["+factorRight.firstPosition+".."+factorRight.lastPosition+"] rightMaximal="+Reads.isRightMaximal(factorRight.lastPosition,ReadA.id,true)+" isRightMaximal="+factorRight.isRightMaximal());
					}
					continue;
				}
			}
			j++;
			if (j!=i) {
				interval=AlignmentIntervals.intervals[j];
				AlignmentIntervals.intervals[j]=AlignmentIntervals.intervals[i];
				AlignmentIntervals.intervals[i]=interval;
			}
		}
		lastInterval=j;
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) AlignmentIntervals.intervals[i].id=i;  // Compacting IDs
		
		// Assigning to each alignment interval its shortest containing periodic interval
		// of any type, if any.
		AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
		if (lastInterval>0) Arrays.sort(AlignmentIntervals.intervals,0,lastInterval+1);
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION || PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		}
		AlignmentIntervals.markIntervalsInPeriodic(0,lastInterval,false);
		if (AlignmentIntervals.lastInterval>lastInterval) Arrays.sort(AlignmentIntervals.intervals,lastInterval+1,AlignmentIntervals.lastInterval+1);
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION || PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		}
		AlignmentIntervals.markIntervalsInPeriodic(lastInterval+1,AlignmentIntervals.lastInterval,false);
		AlignmentIntervals.markIntervalsInDenseSubstrings();
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) AlignmentIntervals.intervals[i].denseSubstring=AlignmentIntervals.intervals[i].inDenseSubstringOfSubstringType;
		
		// Removing discarded intervals from factors
		for (i=0; i<=lastFactor; i++) {
			k=-1;
			for (j=0; j<=factors[i].lastAlignmentInterval; j++) {
				if (factors[i].alignmentIntervals[j].discarded) continue;
				k++;
				if (k!=j) {
					interval=factors[i].alignmentIntervals[k];
					factors[i].alignmentIntervals[k]=factors[i].alignmentIntervals[j];
					factors[i].alignmentIntervals[j]=interval;
				}
			}
			factors[i].lastAlignmentInterval=k;
		}
		if (!discard) {
			AlignmentInterval.order=AlignmentInterval.UNSORTED;
			return;
		}
		
		// Updating pointers from alignments: trying to assign alignments to surviving
		// intervals and, if this fails, to a containing dense substring of substring
		// type, or periodic interval.
		AlignmentIntervals.lastInterval=lastInterval;
		AlignmentInterval.order=AlignmentInterval.FIRSTPOSITION;
		if (AlignmentIntervals.lastInterval>0) Arrays.sort(AlignmentIntervals.intervals,0,AlignmentIntervals.lastInterval+1);
		AlignmentIntervals.discardIntervals_reassignAlignments(ReadA.lastSortedAlignment+1,3,true,false);
	}

	
	/**
	 * Computes the average error rate for short-period intervals, and for alignment 
	 * intervals contained in a short-period range.
	 *
	 * Remark: the procedure assumes that $inPeriodicRange$ is correctly set for all
	 * alignment intervals contained in a short-period range.
	 */
	private static final void assignAlignmentIntervals_getErrorRate() {
		int i;
		int lastRate;
		AlignmentInterval alignmentInterval;
		PeriodicSubstringInterval periodicInterval;
		
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) {
			periodicInterval=PeriodicSubstrings.intervals[i];
			if (periodicInterval.hasLongPeriod) continue;
			periodicInterval.errorRate=0;
			periodicInterval.nImpliedAlignments=0;
		}
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) {
			alignmentInterval=AlignmentIntervals.intervals[i];
			if (!alignmentInterval.inPeriodicRange) continue;
			alignmentInterval.errorRate=0;
			alignmentInterval.nMergedIntervals=0;
		}
		
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			alignmentInterval=ReadA.sortedAlignments[i].mergedToInterval;
			if (alignmentInterval!=null && alignmentInterval.inPeriodicRange) {
				alignmentInterval.nMergedIntervals++;
				alignmentInterval.errorRate+=ReadA.sortedAlignments[i].errorRate();
			}
			periodicInterval=ReadA.sortedAlignments[i].periodicSubstringInterval;
			if (periodicInterval!=null) {
				periodicInterval.nImpliedAlignments++;
				periodicInterval.errorRate+=ReadA.sortedAlignments[i].errorRate();
			}
		}
	}

	
	/**
	 * @return the number of positions in $interval$ that also belong to one of its
	 * factors.
	 */
	private static final int getFactorsSurface(AlignmentInterval interval) {
		int i;
		int surface, from, to;
	
		surface=0;
		for (i=0; i<=interval.lastFactor; i++) {
			from=Math.max(interval.firstPosition,factors[interval.factors[i]].firstPosition);
			to=Math.min(interval.lastPosition,factors[interval.factors[i]].lastPosition);
			surface+=to-from+1;
		}
		return surface;
	}
	
	
	/**
	 * @return TRUE iff both ends of $interval$ are close to the corresponding end of a 
	 * factor assigned to the interval, or if one end is close and the other is adjacent
	 * to a low-quality region of readA.
	 */
	private static final boolean factorBoundaryOrLowQuality(AlignmentInterval interval) {
		final int THRESHOLD = IO.quantum;
		boolean foundFirst, foundLast;
		int i;
		
		foundFirst=false; foundLast=false;
		for (i=0; i<=interval.lastFactor; i++) {
			if (Math.abs(interval.firstPosition,factors[interval.factors[i]].firstPosition)<=THRESHOLD) foundFirst=true;
			if (Math.abs(interval.lastPosition,factors[interval.factors[i]].lastPosition)<=THRESHOLD) foundLast=true;
		}
		return (foundFirst && foundLast) || 
			   (foundFirst && !Reads.isRightMaximal(interval.lastPosition,ReadA.id,true)) ||
			   (foundLast && !Reads.isLeftMaximal(interval.firstPosition,ReadA.id,true));
	}
	
	
	/**
	 * Remark: the procedure assumes $factors$ to be sorted by $firstPosition$.
	 *
	 * @param leftOrRight TRUE=left, FALSE=right.
	 */
	private static final boolean isBracketed(Factor factor, int position, boolean leftOrRight) {
		final int THRESHOLD = IO.quantum;
		int i;
		
		if ((leftOrRight&&factor.isLeftBracketed) || (!leftOrRight&&factor.isRightBracketed)) return true;
		i=Arrays.binarySearch(factors,0,lastFactor+1,factor);
		if (i<0) {
			System.err.println("Factors.isBracketed> ERROR: factor not found?!");
			System.exit(1);
		}
		if ( (leftOrRight && i<lastFactor && Math.abs(position,factors[i+1].firstPosition)<=THRESHOLD && factors[i+1].isLeftBracketed) ||
			 (!leftOrRight && i>0 && Math.abs(position,factors[i-1].lastPosition)<=THRESHOLD && factors[i-1].isRightBracketed)
		   ) return true;
		return false;
	}


	/**
	 * Field $shouldBeOutput$ of a factor is set to false if the factor is not periodic
	 * and its corrected coverage is smaller than $IO.minRepeatCoverage$, or if the
	 * factor is periodic and its corrected coverage is smaller than
	 * $minPeriodicCoverage$.
	 *
	 * @return the number of low-coverage factors discarded.
	 */
	private static final int discardLowCoverageFactors() {
		int i, out, coverage;

		out=0;
		for (i=0; i<=lastFactor; i++) {
			coverage=factors[i].getCorrectedCoverage();
			if ( (factors[i].periodicCoverage>0 && coverage<minPeriodicCoverage) ||
			     (factors[i].periodicCoverage==0 && coverage<IO.minRepeatCoverage) ) {
				factors[i].shouldBeOutput=false;
				out++;
			}
		}
		return out;
	}


	/**
	 * Field $shouldBeOutput$ of a factor is set to false if the factor has low quality.
	 *
	 * @return the number of low-quality factors discarded.
	 */
	private static final int discardLowQualityFactors() {
		final int MIN_COVERAGE = (IO.minOccurrencesInGenome<<1)*IO.coverage-1;
		int i, out;

		out=0;
		for (i=0; i<=lastFactor; i++) {
			if (!factors[i].hasHighQuality && ReadA.getAverageCoverage(factors[i].firstPosition,factors[i].lastPosition)<MIN_COVERAGE) {
				factors[i].shouldBeOutput=false;
				out++;
			}
		}
		return out;
	}


	/**
	 * Estimates the period of every periodic factor, using just the short-period periodic
	 * substrings that contain the factor.
	 */
	private static final void estimatePeriods() {
		int i, j, firstJForNextFactor;
		Point[] tmp;

		// Ensuring the necessary order in $PeriodicSubstrings.substrings$
		if (PeriodicSubstring.order!=PeriodicSubstring.MINSTARTA && PeriodicSubstring.order!=PeriodicSubstring.MINSTARTA_MAXENDA) {
			PeriodicSubstring.order=PeriodicSubstring.MINSTARTA;
			if (PeriodicSubstrings.lastSubstring>0) Arrays.sort(PeriodicSubstrings.substrings,0,PeriodicSubstrings.lastSubstring+1);
		}

		// Estimating periods
		for (i=0; i<=lastFactor; i++) factors[i].period=-1;
		i=0; j=0; firstJForNextFactor=-1; lastShift=-1;
		while (i<=lastFactor) {
			if (factors[i].periodicCoverage==0) {
				i++;
				firstJForNextFactor=-1;
				lastShift=-1;
				continue;
			}
			if (j>PeriodicSubstrings.lastSubstring || PeriodicSubstrings.substrings[j].minStartA>factors[i].lastPosition) {
				if (lastShift!=-1) factors[i].period=PeriodicSubstrings.estimatePeriod(shifts,lastShift,PeriodicSubstrings.MIN_PERIOD,PeriodicSubstrings.minIntervalLength,PeriodicSubstrings.minLocalMaxDistance);
				i++;
				if (firstJForNextFactor!=-1) j=firstJForNextFactor;
				firstJForNextFactor=-1;
				lastShift=-1;
				continue;
			}
			if (PeriodicSubstrings.substrings[j].maxEndA<factors[i].firstPosition) {
				j++;
				continue;
			}
			if (i<lastFactor && firstJForNextFactor==-1 && PeriodicSubstrings.substrings[j].maxEndA>=factors[i+1].firstPosition) firstJForNextFactor=j;
			if ( Intervals.isApproximatelyContained(factors[i].firstPosition,factors[i].lastPosition,PeriodicSubstrings.substrings[j].minStartA,PeriodicSubstrings.substrings[j].maxEndA) ||
				 Intervals.areApproximatelyIdentical(factors[i].firstPosition,factors[i].lastPosition,PeriodicSubstrings.substrings[j].minStartA,PeriodicSubstrings.substrings[j].maxEndA)
			   ) {
				if (lastShift==-1) {
					Points.simpleClone(PeriodicSubstrings.substrings[j].shifts,PeriodicSubstrings.substrings[j].lastShift,shifts);
					lastShift=PeriodicSubstrings.substrings[j].lastShift;
				}
				else {
					lastShift=Points.merge(shifts,lastShift,PeriodicSubstrings.substrings[j].shifts,PeriodicSubstrings.substrings[j].lastShift,tmpShifts);
					tmp=shifts;
					shifts=tmpShifts;
					tmpShifts=tmp;
				}
			}
			j++;
		}
	}


	/**
	 * Estimates whether periodic factors $firstFactor$ and $firstFactor+1$ in array
	 * $Factors$ have identical periods. This is done by finding alignments $A1$ and $A2$
	 * such that: (1) they have the same readB and orientation; (2) $A1$ intersects
	 * $firstFactor$, $A2$ intersects $firstFactor+1$; (3) the projections on readB of
	 * their intersections with the factors have a high Jaccard similarity, or are
	 * contained in one another; (4) the average diffs of $A1$ (respectively, $A2$) is at
	 * most $factors[firstFactor].maxDiffs$ (respectively, $factors[firstFactor+1].maxDiffs$).
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments$ to be sorted by $startA$.
	 *
	 * Remark: the procedure considers only alignments whose position in
	 * $ReadA.sortedAlignments$ is at least $firstAlignmentForCurrentFactor$. At the end
	 * of the procedure, $firstAlignmentForNextFactor$ is set to the first alignment in
	 * $ReadA.sortedAlignments$ that intersects $firstFactor+1$.
	 *
	 * Remark: the procedure stores the last value of $firstFactor$ on which it has been
	 * called, and the corresponding output. If it is called again on the same value of
	 * $firstFactor$, it returns the cached output and it quits, without updating
	 * $firstAlignmentForNextFactor$. If it is called on the last value of $firstFactor$
	 * plus one, it reuses the intervals already stored in $intervals2$.
	 */
	private static final boolean haveIdenticalPeriod(int firstFactor) {
		int i, j, firstJ, firstJForNextI;
		int alignmentID, alignmentStart, alignmentEnd, alignmentLength;
		int interval1Length, interval2Length;
		double intersection;
		Interval[] tmpIntervals;

		// Cached result
		if (firstFactor==lastIdenticalPeriodFirstFactor) {
			if (IO.SHOW_STD_ERR) IO.printErr("haveIdenticalPeriod("+firstFactor+")="+lastIdenticalPeriodOutput);
			return lastIdenticalPeriodOutput;
		}

		// Collecting intervals in $readB$ induced by the intersection of $firstFactor$
		// and $lastFactor$ with alignments
		if (firstFactor==lastIdenticalPeriodFirstFactor+1) {
			tmpIntervals=intervals1;
			intervals1=intervals2;
			intervals2=tmpIntervals;
			lastInterval1=lastInterval2;
		}
		else lastInterval1=-1;
		lastInterval2=-1; firstAlignmentForNextFactor=-1;
		for (i=firstAlignmentForCurrentFactor; i<=ReadA.lastSortedAlignment; i++) {
			alignmentID=ReadA.sortedAlignments[i].id;
			alignmentStart=Alignments.alignments[alignmentID][3];
			if (alignmentStart>factors[firstFactor+1].lastPosition) break;
			alignmentEnd=Alignments.alignments[alignmentID][4];
			if (alignmentEnd<factors[firstFactor].firstPosition) continue;
			if (firstAlignmentForNextFactor==-1 && alignmentEnd>=factors[firstFactor+1].firstPosition) firstAlignmentForNextFactor=i;
			alignmentLength=alignmentEnd-alignmentStart+1;
			if (firstFactor!=lastIdenticalPeriodFirstFactor+1) {
				if ( Intervals.isApproximatelyContained(alignmentStart,alignmentEnd,factors[firstFactor].firstPosition,factors[firstFactor].lastPosition) ||
				     Intervals.isApproximatelyContained(factors[firstFactor].firstPosition,factors[firstFactor].lastPosition,alignmentStart,alignmentEnd) ) {
					lastInterval1++;
					intervals1[lastInterval1].firstPosition=Alignments.getF(alignmentID,factors[firstFactor].firstPosition);
					intervals1[lastInterval1].lastPosition=Alignments.getG(alignmentID,factors[firstFactor].lastPosition);
					intervals1[lastInterval1].readB=Alignments.alignments[alignmentID][1]-1;
					intervals1[lastInterval1].orientation=Alignments.alignments[alignmentID][2]==1;
					intervals1[lastInterval1].alignmentID=alignmentID;
					intervals1[lastInterval1].avgDiffs=Alignments.getAvgDiffs(alignmentID);
				}
			}
			if ( Intervals.isApproximatelyContained(alignmentStart,alignmentEnd,factors[firstFactor+1].firstPosition,factors[firstFactor+1].lastPosition) ||
				 Intervals.isApproximatelyContained(factors[firstFactor+1].firstPosition,factors[firstFactor+1].lastPosition,alignmentStart,alignmentEnd) ) {
			    lastInterval2++;
			    intervals2[lastInterval2].firstPosition=Alignments.getF(alignmentID,factors[firstFactor+1].firstPosition);
				intervals2[lastInterval2].lastPosition=Alignments.getG(alignmentID,factors[firstFactor+1].lastPosition);
				intervals2[lastInterval2].readB=Alignments.alignments[alignmentID][1]-1;
				intervals2[lastInterval2].orientation=Alignments.alignments[alignmentID][2]==1;
				intervals2[lastInterval2].alignmentID=alignmentID;
				intervals2[lastInterval2].avgDiffs=Alignments.getAvgDiffs(alignmentID);
			}
		}

		// Finding two compatible intervals from $firstFactor$ and $firstFactor+1$
		Interval.order=Interval.READB_ORIENTATION_FIRSTPOSITION;
		if (firstFactor!=lastIdenticalPeriodFirstFactor+1) {
			if (lastInterval1>0) Arrays.sort(intervals1,0,lastInterval1+1);
		}
		if (lastInterval2>0) Arrays.sort(intervals2,0,lastInterval2+1);



// if (IO.SHOW_STD_ERR) IO.printErr("intervals1 (length="+(lastInterval1+1)+"): ");
// for (i=0; i<lastInterval1; i++) if (IO.SHOW_STD_ERR) IO.printErr("(i="+i+","+intervals1[i]+") ");
// if (IO.SHOW_STD_ERR) IO.printErr("intervals2 (length="+(lastInterval2+1)+"): ");
// for (i=0; i<lastInterval2; i++) if (IO.SHOW_STD_ERR) IO.printErr("(i="+i+","+intervals2[i]+") ");



double max=0.0;


		firstJ=0;
		for (i=0; i<=lastInterval1; i++) {
			interval1Length=intervals1[i].lastPosition-intervals1[i].firstPosition+1;
			firstJForNextI=-1;
			for (j=firstJ; j<=lastInterval2; j++) {
				if ( intervals2[j].readB>intervals1[i].readB ||
				     (!intervals2[j].orientation && intervals1[i].orientation) ||
				     intervals2[j].firstPosition>intervals1[i].lastPosition
				   ) {
					if (firstJForNextI!=-1) firstJ=firstJForNextI;
					else firstJ=j;
					break;
				}
				if ( intervals2[j].readB<intervals1[i].readB ||
				     (intervals2[j].orientation && !intervals1[i].orientation) ||
				     intervals2[j].lastPosition<intervals1[i].firstPosition
				   ) continue;
				if ( firstJForNextI==-1 && i<lastInterval1 &&
				     intervals1[i+1].readB==intervals1[i].readB &&
				     intervals1[i+1].orientation==intervals1[i].orientation &&
				     intervals2[j].lastPosition>=intervals1[i+1].firstPosition ) firstJForNextI=j;
				interval2Length=intervals2[j].lastPosition-intervals2[j].firstPosition+1;
				intersection=Intervals.intersectionLength(intervals1[i].firstPosition,intervals1[i].lastPosition,intervals2[j].firstPosition,intervals2[j].lastPosition);

if (intersection/interval1Length>max) max=intersection/interval1Length;
if (intersection/interval2Length>max) max=intersection/interval2Length;
if (intersection>interval1Length) if (IO.SHOW_STD_ERR) IO.printErr("ERROR: intersection="+intersection+" interval1=["+intervals1[i].firstPosition+".."+intervals1[i].lastPosition+"] interval2=["+intervals2[j].firstPosition+".."+intervals2[j].lastPosition+"] length1="+interval1Length+" length2="+interval2Length);
if (intersection>interval2Length) if (IO.SHOW_STD_ERR) IO.printErr("ERROR: intersection="+intersection+" interval2=["+intervals2[j].firstPosition+".."+intervals2[j].lastPosition+"] interval1=["+intervals1[i].firstPosition+".."+intervals1[i].lastPosition+"] length2="+interval2Length+" length1="+interval1Length);


if (IO.SHOW_STD_ERR) IO.printErr("HIP> "+intervals1[i].avgDiffs+"::"+factors[firstFactor].maxDiffs+", "+intervals2[j].avgDiffs+"::"+factors[firstFactor+1].maxDiffs);

				if ( ( Intervals.isApproximatelyContained(intervals1[i].firstPosition,intervals1[i].lastPosition,intervals2[j].firstPosition,intervals2[j].lastPosition) ||
				       Intervals.isApproximatelyContained(intervals2[j].firstPosition,intervals2[j].lastPosition,intervals1[i].firstPosition,intervals1[i].lastPosition)
				     ) &&
				     (intervals1[i].avgDiffs<=factors[firstFactor].maxDiffs && intervals2[j].avgDiffs<=factors[firstFactor+1].maxDiffs) ) {
				   	lastIdenticalPeriodFirstFactor=firstFactor;
					lastIdenticalPeriodOutput=true;
					if (IO.SHOW_STD_ERR) IO.printErr("haveIdenticalPeriod("+firstFactor+")="+lastIdenticalPeriodOutput+" because of alignments "+intervals1[i].alignmentID+" and "+intervals2[j].alignmentID+" in LAshow, max="+max);
					return lastIdenticalPeriodOutput;
				}
			}
		}
		lastIdenticalPeriodFirstFactor=firstFactor;
		lastIdenticalPeriodOutput=false;
		if (IO.SHOW_STD_ERR) IO.printErr("haveIdenticalPeriod("+firstFactor+")="+lastIdenticalPeriodOutput+", max="+max);
		return lastIdenticalPeriodOutput;
	}


	/**
	 * Assigns to every factor an estimated lower bound on the number of its variants.
	 */
	private static final void estimateNVariants() {
		int i, firstAlignment;

		// Ensuring the necessary order in $ReadA.sortedAlignments$
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

		// Estimating number of variants
		firstAlignment=0;
		for (i=0; i<=lastFactor; i++) {
if (IO.SHOW_STD_ERR) IO.printErr("estimateNVariants> factor ["+factors[i].firstPosition+".."+factors[i].lastPosition+"]"); 
			firstAlignment=ReadA.estimateMaxDiffs(factors[i].firstPosition,factors[i].lastPosition,firstAlignment,ReadA.lastSortedAlignment,true);
			factors[i].maxDiffs=ReadA.maxDiffs;
			factors[i].nVariants=ReadA.nVariants;
		}
	}
	
	
	/**
	 * @return TRUE iff the factor that follows $factor$ in sequence order exists and is 
	 * left-maximal.
	 */
	private static final boolean isNextFactorMaximal(Factor factor, boolean leftOrRight, int threshold) {
		int i, j;
		
		if (Factor.order==Factor.FIRSTPOSITION) i=Arrays.binarySearch(factors,0,lastFactor+1,factor);
		else {
			i=-1;
			for (j=0; j<=lastFactor; j++) {
				if (factors[j].firstPosition==factor.firstPosition && factors[j].lastPosition==factor.lastPosition) {
					i=j;
					break;
				}
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			if (i<0) {
				System.err.println("isNextFactorLeftMaximal> ERROR: factor not found in the list of factors.");
				System.exit(1);
			}
			if (factors[i].lastPosition!=factor.lastPosition) {
				System.err.println("isNextFactorLeftMaximal> ERROR: end positions of factors do not match?!");
				System.exit(1);
			}
		}
		
		if (i>0) {
			if (leftOrRight) {
				if (factors[i-1].firstPosition>=factors[i].firstPosition-threshold && factors[i-1].isLeftMaximal()) return true;
			}
			if (!leftOrRight) {
				if (factors[i-1].lastPosition>=factors[i].lastPosition-threshold && factors[i-1].isRightMaximal()) return true;
			}
		}
		if (i<lastFactor) {
			if (leftOrRight) {
				if (factors[i+1].firstPosition<=factors[i].firstPosition+threshold && factors[i+1].isLeftMaximal()) return true;
			}
			if (!leftOrRight) {
				if (factors[i+1].lastPosition<=factors[i].lastPosition+threshold && factors[i+1].isRightMaximal()) return true;
			}
		}
		return false;
	}
	
	
	private static final boolean checkSplitsOfNextFactor(Factor factor, boolean leftOrRight, int threshold) {
		int i, j;
		
		if (Factor.order==Factor.FIRSTPOSITION) i=Arrays.binarySearch(factors,0,lastFactor+1,factor);
		else {
			i=-1;
			for (j=0; j<=lastFactor; j++) {
				if (factors[j].firstPosition==factor.firstPosition && factors[j].lastPosition==factor.lastPosition) {
					i=j;
					break;
				}
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			if (i<0) {
				System.err.println("isNextFactorLeftMaximal> ERROR: factor not found in the list of factors.");
				System.exit(1);
			}
			if (factors[i].lastPosition!=factor.lastPosition) {
				System.err.println("isNextFactorLeftMaximal> ERROR: end positions of factors do not match?!");
				System.exit(1);
			}
		}
		
		if (i>0) {
			if (leftOrRight) {
				if (factors[i-1].firstPosition>=factors[i].firstPosition-threshold && factors[i-1].checkSplits(leftOrRight)) return true;
			}
			if (!leftOrRight) {
				if (factors[i-1].lastPosition>=factors[i].lastPosition-threshold && factors[i-1].checkSplits(leftOrRight)) return true;
			}
		}
		if (i<lastFactor) {
			if (leftOrRight) {
				if (factors[i+1].firstPosition<=factors[i].firstPosition+threshold && factors[i+1].checkSplits(leftOrRight)) return true;
			}
			if (!leftOrRight) {
				if (factors[i+1].lastPosition<=factors[i].lastPosition+threshold && factors[i+1].checkSplits(leftOrRight)) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Resets to zero all the IDs assigned while detecting factors.
	 */
	private static final void resetIDs() {
		int i;
		
		for (i=0; i<=DenseSubstrings.lastSubstring; i++) DenseSubstrings.substrings[i].id=0;
		for (i=0; i<=PeriodicSubstrings.lastInterval; i++) PeriodicSubstrings.intervals[i].id=0;
		for (i=0; i<=AlignmentIntervals.lastInterval; i++) AlignmentIntervals.intervals[i].id=0;
	}
	
	
	/**
	 * Remark: the procedure assumes $splits$ to be sorted by position.
	 *
	 * @return TRUE iff $pos$ is at distance at most $threshold$ from an element of 
	 * $splits$.
	 */
	public static final boolean isCloseToSplit(int pos, int threshold) {
		int j;
		
		tmpSplit.position=pos;
		j=Arrays.binarySearch(splits,0,lastSplit+1,tmpSplit);
		if (j>=0) return true;
		j=-j-1;
		if (j<=lastSplit && splits[j].position<=pos+threshold) return true;
		if (j-1>=0 && splits[j-1].position>=pos-threshold) return true;
		return false;
	}
	

}