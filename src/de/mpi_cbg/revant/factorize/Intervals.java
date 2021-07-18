package de.mpi_cbg.revant.factorize;

import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Leaf;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Leaves;


public class Intervals {
	/**
	 * Parameters of the pipeline
	 */
	public static double jaccardThreshold = 0.8;  // Lower bound on the Jaccard similarity between two intervals for them to be considered approximately the same.
	public static double relativeThreshold = 0.25;
	public static double relativeThresholdPrime = 0.5;
	public static int absoluteThreshold = IO.quantum;
	public static int absoluteThresholdPrime = IO.quantum<<1;
	
	/**
	 * Temporary space
	 */
	private static double[] tmp = new double[2];
	private static Leaf[] tmpLeaves;
	
	
	public static final void allocateMemory(int maxLeaves) {
		tmpLeaves = new Leaf[maxLeaves];
		for (int i=0; i<maxLeaves; i++) tmpLeaves[i] = new Leaf();
	}


	public static final boolean intersect(int start1, int end1, int start2, int end2) {
		if (start2>end1 || start1>end2) return false;
		return true;
	}


	/**
	 * @return the length of the intersection between the intervals $[start1..start2]$
	 * and $[end1..end2]$.
	 */
	public static final int intersectionLength(int start1, int end1, int start2, int end2) {
		if (start2>end1 || start1>end2) return 0;
		return Math.min(end2,end1)-Math.max(start2,start1)+1;
	}


	/**
	 * @return the length of the union of the intervals $[start1..start2]$
	 * and $[end1..end2]$.
	 */
	public static final int unionLength(int start1, int end1, int start2, int end2) {
		if (start2>end1 || start1>end2) return end1-start1+1 + end2-start2+1;
		return Math.max(end2,end1)-Math.min(start2,start1)+1;
	}


	/**
	 * @return the Jaccard similarity between intervals $[start1..end1]$ and
	 * $[start2..end2]$.
	 */
	public static final double jaccardSimilarity(int start1, int end1, int start2, int end2) {
		double intersection = intersectionLength(start1,end1,start2,end2);
		return intersection/unionLength(start1,end1,start2,end2);
	}
	
	
	/**
	 * @return the size of the difference between union and intersection of the intervals.
	 */
	public static final double difference(int start1, int end1, int start2, int end2) {
		return unionLength(start1,end1,start2,end2)-intersectionLength(start1,end1,start2,end2);
	}


	/**
	 * @return TRUE iff $[start1..end1] \subseteq [start2..end2]$.
	 */
	public static final boolean isContained(int start1, int end1, int start2, int end2) {
		return start1>=start2 && end1<=end2;
	}


	/**
	 * @return TRUE iff $[start1..end1]$ strictly contains $[start2..end2]$, if 
	 * $[start1..end1] \setminus [start2..end2]$ has low quality, and if $[start1..end1]$ 
	 * is not too big compared to $[start2..end2]$.
	 */
	public static final boolean contains_lowQuality(int start1, int end1, int start2, int end2, int readID) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		final int LENGTH_THRESHOLD = Math.min(length1,length2);
		if (length1<=length2 || length1>length2+LENGTH_THRESHOLD) return false;
		
		return (start1<start2-absoluteThreshold && Reads.hasLowQuality(readID,start1,start2-1,true)) && 
			   (end1>end2+absoluteThreshold && Reads.hasLowQuality(readID,end2+1,end1,true));
	}


	/**
	 * @return TRUE iff $[start1..end1] \subseteq [start2..end2]$ or
	 * $[start2..end2] \subseteq [start1..end1]$.
	 */
	public static final boolean eitherIsContained(int start1, int end1, int start2, int end2) {
		return (start2>=start1 && end2<=end1) || (start1>=start2 && end1<=end2);
	}


	/**
	 * @return TRUE iff $[start1..end1]$ is approximately contained in $[start2..end2]$.
	 */
	public static final boolean isApproximatelyContained(int start1, int end1, int start2, int end2) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		if (length1>=length2) return false;
		final int DISTANCE_THRESHOLD = absoluteThreshold;
		final int SMALL_DISTANCE_THRESHOLD = 10;
			
		return (start1>start2-SMALL_DISTANCE_THRESHOLD && end1<=end2+DISTANCE_THRESHOLD) || (end1<end2+SMALL_DISTANCE_THRESHOLD && start1>=start2-DISTANCE_THRESHOLD);
	}

	
	/**
	 * @return TRUE iff $[start1..end1]$ is approximately contained in $[start2..end2]$, 
	 * or if $[start1..end1] \setminus [start2..end2]$ has low quality, and if the vice 
	 * versa is not true.
	 */
	public static final boolean isApproximatelyContained_lowQuality(int start1, int end1, int start2, int end2, int readID) {
		return isApproximatelyContained_lowQuality_impl(start1,end1,start2,end2,readID) && 
			   !isApproximatelyContained_lowQuality_impl(start2,end2,start1,end1,readID);
	}
	
	
	private static final boolean isApproximatelyContained_lowQuality_impl(int start1, int end1, int start2, int end2, int readID) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		final int LENGTH_THRESHOLD = (int)(Math.min(length1,length2)*0.5);
		if (length1>length2+LENGTH_THRESHOLD) return false;
		final int DISTANCE_THRESHOLD = absoluteThreshold;
		final int SMALL_DISTANCE_THRESHOLD = 10;  // Arbitrary
		final double INTERSECTION_THRESHOLD = 0.5;  // Arbitrary
		
		if (intersectionLength(start1,end1,start2,end2)<length1*INTERSECTION_THRESHOLD) return false;
		return ( start1>start2-SMALL_DISTANCE_THRESHOLD && 
			     ( end1<=end2+DISTANCE_THRESHOLD || 
				   ( start1<end2-DISTANCE_THRESHOLD && Reads.hasLowQuality(readID,end2+1,end1,true) ) 
				 ) 
			   ) ||
			   ( end1<end2+SMALL_DISTANCE_THRESHOLD && 
				 ( start1>=start2-DISTANCE_THRESHOLD || 
				   ( end1>start2+DISTANCE_THRESHOLD && Reads.hasLowQuality(readID,start1,start2-1,true) )
				 ) 
			   );
	}


	/**
	 * @return TRUE iff $[start1..end1]$ and $[start2..end2]$ are approximately identical.
	 */
	public static final boolean areApproximatelyIdentical(int start1, int end1, int start2, int end2) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		final int DISTANCE_THRESHOLD = Math.min(absoluteThreshold,(int)(Math.min(length1,length2)*0.25));
		
		return Math.abs(start1,start2)<=DISTANCE_THRESHOLD && Math.abs(end1,end2)<=DISTANCE_THRESHOLD;
	}
	
	
	/**
	 * @return TRUE iff at least one of the ends of $[start1..end1]$ and $[start2..end2]$ 
	 * coincide, and if the difference on the other end has low quality.
	 */
	public static final boolean areApproximatelyIdentical_lowQuality(int start1, int end1, int start2, int end2, int readID) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		final int LENGTH_THRESHOLD = (int)(Math.min(length1,length2)*0.5);
		if (Math.abs(length1,length2)>LENGTH_THRESHOLD) return false;
		final int DISTANCE_THRESHOLD = Math.min(absoluteThreshold,(int)(Math.min(length1,length2)*0.25));
		
		return ( Math.abs(start1,start2)<=DISTANCE_THRESHOLD && 
			     ( Math.abs(end1,end2)<=DISTANCE_THRESHOLD || Reads.hasLowQuality(readID,Math.min(end1,end2)+1,Math.max(end1,end2),true) )
			   ) ||
			   ( Math.abs(end1,end2)<=DISTANCE_THRESHOLD && 
				 ( Math.abs(start1,start2)<=DISTANCE_THRESHOLD || Reads.hasLowQuality(readID,Math.min(start1,start2),Math.max(start1,start2)-1,true) )
			   );
	}
	
	
	/**
	 * @return TRUE if intervals $[start1..end1]$ and $[start2..end2]$ in readA coincide,
	 * or if the coverage in their symmetric difference is low compared to the average 
	 * coverage in the intersection of the intervals.
	 */
	public static final boolean areApproximatelyIdentical_lowCoverage(int start1, int end1, int start2, int end2) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		final int LENGTH_THRESHOLD = (int)(Math.min(length1,length2)*0.5);
		if (Math.abs(length1,length2)>LENGTH_THRESHOLD) return false;
		final int DISTANCE_THRESHOLD = Math.min(absoluteThreshold,(int)(Math.min(length1,length2)*0.25));
		
		final double COVERAGE_RATIO = 0.25;  // Arbitrary
		final int left = Math.max(start1,start2);
		final int right = Math.min(end1,end2);
		final double coverage = ReadA.getAverageCoverage(left,right);
		return ( Math.abs(start1,start2)<=DISTANCE_THRESHOLD || ReadA.getAverageCoverage(Math.min(start1,start2),left-1)<=coverage*COVERAGE_RATIO ) &&
			   ( Math.abs(end1,end2)<=DISTANCE_THRESHOLD || ReadA.getAverageCoverage(right+1,Math.max(end1,end2))<=coverage*COVERAGE_RATIO );
	}
	
	
	/**
	 * @return TRUE iff $[start1..end1]$ straddles $[start2..end2]$ in any direction.
	 */
	public static final boolean straddles(int start1, int end1, int start2, int end2) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		final int DISTANCE_THRESHOLD = Math.min(absoluteThreshold,(int)(Math.min(length1,length2)*0.25));
		
		return (end1<end2-DISTANCE_THRESHOLD && end1>start2+DISTANCE_THRESHOLD && start1<start2-DISTANCE_THRESHOLD) || 
			   (start1>start2+DISTANCE_THRESHOLD && start1<end2-DISTANCE_THRESHOLD && end1>end2+DISTANCE_THRESHOLD);
	}
	
	
	/**
	 * Let $[startAPrime..endAPrime] \subseteq [startA..endA]$ be intervals, and let 
	 * $left$ (respectively, $right$) collect a set of B-left- (respectively B-right-) 
	 * maximal events inside $[startA..endA]$ (not necessarily sorted and compacted). The 
	 * procedure resets $startA$ (respectively, $endA$) to the high peak of $left$
	 * (respectively, $right$) that is closest to $startAPrime$ (respectively,
	 * $endAPrime$).
	 *
	 * @return $out[0]$ and $out[1]$ store the new boundaries.
	 */
	public static final void refineBoundaries(int startA, int endA, int startAPrime, int endAPrime, Point[] left, int lastLeft, Point[] right, int lastRight, int binLength, int minIntervalLength, int minLocalMaxDistance, int[] out) {
		final int MIN_LENGTH = Alignments.minAlignmentLength;
		final int MAX_REFINEMENT_DISTANCE = Alignments.minAlignmentLength;
		final int MIN_HOLE_LENGTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		final int SMALL_LEAF_DISTANCE = IO.quantum>>2;  // Arbitrary
		final int MAX_WIDTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		int newStartA, newEndA;
		int largestMass, nLocalMaximumLeaves, nHigh, lastLeaf;
		double newStartADensity, newEndADensity;
		Leaf[] lvs;
		
		newStartA=-1; newEndA=-1;
		newStartADensity=-1; newEndADensity=-1;
		lastLeft=Points.sortAndCompact(left,lastLeft);
		if (!Points.areUniformlyDistributed(left,0,lastLeft,true,binLength)) {
			largestMass=Points.largestMass(left,lastLeft,-1);
			nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(left,0,lastLeft,minIntervalLength,minLocalMaxDistance,true,largestMass,MIN_HOLE_LENGTH,SMALL_LEAF_DISTANCE,true,true);		


if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("refineBoundaries> DET left:");
	for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) System.err.println(DensityEstimationTree.leaves[x]);
	System.err.println();
	System.err.println("refineBoundaries> points left:");
	for (int x=0; x<=lastLeft; x++) System.err.println(left[x]);
	System.err.println();
}

			if (nLocalMaximumLeaves>0) {
				lvs=DensityEstimationTree.leaves;
				DensityEstimationTree.leaves=tmpLeaves;
				tmpLeaves=lvs;
				lastLeaf=DensityEstimationTree.lastLeaf;
				nHigh=Leaves.setHighLocalMax(tmpLeaves,lastLeaf,lastLeft,true,false,left,Points.mass(left,0,lastLeft),true,MAX_WIDTH);
				if (nHigh>0) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("refineBoundaries> nHigh="+nHigh);					
					lvs=tmpLeaves;
					tmpLeaves=DensityEstimationTree.leaves;
					DensityEstimationTree.leaves=lvs;
					DensityEstimationTree.lastLeaf=lastLeaf;
					DensityEstimationTree.markRunsOfLocalMaximumLeaves(left);
					
if (IO.SHOW_STD_ERR_PRIME) {				
	System.err.println("refineBoundaries> DET left prime:");
	for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) System.err.println(DensityEstimationTree.leaves[x]);
	System.err.println();					
}

					DensityEstimationTree.getCenterOfMassOfHighRun(startAPrime,true,left,tmp);
					newStartA=(int)tmp[0];
					newStartADensity=tmp[1];
					
if (IO.SHOW_STD_ERR_PRIME) System.err.println("refineBoundaries> newStartA="+newStartA+" newStartADensity="+newStartADensity);					
					
				}
			}
		}
		lastRight=Points.sortAndCompact(right,lastRight);
		if (!Points.areUniformlyDistributed(right,0,lastRight,true,binLength)) {
			largestMass=Points.largestMass(right,lastRight,-1);
			nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(right,0,lastRight,minIntervalLength,minLocalMaxDistance,true,largestMass,-1,-1,true,true);			
			
if (IO.SHOW_STD_ERR_PRIME) {
	System.err.println("refineBoundaries> DET right:");
	for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) System.err.println(DensityEstimationTree.leaves[x]);
	System.err.println();
	System.err.println("refineBoundaries> points right:");
	for (int x=0; x<=lastRight; x++) System.err.println(right[x]);
	System.err.println();			
}

			if (nLocalMaximumLeaves>0) {
				lvs=DensityEstimationTree.leaves;
				DensityEstimationTree.leaves=tmpLeaves;
				tmpLeaves=lvs;
				lastLeaf=DensityEstimationTree.lastLeaf;
				nHigh=Leaves.setHighLocalMax(tmpLeaves,lastLeaf,lastRight,true,false,right,Points.mass(right,0,lastRight),true,MAX_WIDTH);				
				if (nHigh>0) {
					lvs=tmpLeaves;
					tmpLeaves=DensityEstimationTree.leaves;
					DensityEstimationTree.leaves=lvs;
					DensityEstimationTree.lastLeaf=lastLeaf;
					DensityEstimationTree.markRunsOfLocalMaximumLeaves(right);
					
if (IO.SHOW_STD_ERR_PRIME) {					
	System.err.println("refineBoundaries> DET right prime:");
	for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) System.err.println(DensityEstimationTree.leaves[x]);
	System.err.println();					
}

					DensityEstimationTree.getCenterOfMassOfHighRun(endAPrime,false,right,tmp);
					newEndA=(int)tmp[0];
					newEndADensity=tmp[1];
					

if (IO.SHOW_STD_ERR_PRIME) System.err.println("refineBoundaries> newEndA="+newEndA+" newEndADensity="+newEndADensity);										
					
				}
			}
		}
		if ( newStartA!=-1 && 
			 Math.abs(newStartA,startA)<=MAX_REFINEMENT_DISTANCE && (
			 ( newEndA==-1 && endA-newStartA+1>=MIN_LENGTH ) ||
			 ( newEndA!=-1 && (
			   	newEndA-newStartA+1>=MIN_LENGTH || 
			   	(endA-newStartA+1>=MIN_LENGTH && newStartADensity>=newEndADensity)
			   )
			 )
		   )
		) out[0]=newStartA;
		else out[0]=startA;
		if ( newEndA!=-1 && 
			 Math.abs(newEndA,endA)<=MAX_REFINEMENT_DISTANCE && (
			 ( newStartA==-1 && newEndA-startA+1>=MIN_LENGTH ) ||
			 ( newStartA!=-1 && (
			    newEndA-newStartA+1>=MIN_LENGTH || 
				(newEndA-startA+1>=MIN_LENGTH && newEndADensity>newStartADensity)
			   )
			 )
		   )
		) out[1]=newEndA;
		else out[1]=endA;
		
if (IO.SHOW_STD_ERR_PRIME) System.err.println("refineBoundaries> new boundaries: "+out[0]+".."+out[1]);		
		
	}
	
	
	/**
	 * @param prefixOrSuffix TRUE=prefix, FALSE=suffix;
	 * @return TRUE iff $[start1..end1]$ is approximately a proper prefix (respectively, 
	 * suffix) of $[start2..end2]$.
	 */
	public static final boolean isApproximatePrefixOrSuffix(int start1, int end1, int start2, int end2, boolean prefixOrSuffix) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		if (length1>=length2) return false;
		final int DISTANCE_THRESHOLD = absoluteThreshold;
		
		return prefixOrSuffix?(Math.abs(start1,start2)<=DISTANCE_THRESHOLD&&end1<end2-DISTANCE_THRESHOLD):(Math.abs(end1,end2)<=DISTANCE_THRESHOLD&&start1>start2+DISTANCE_THRESHOLD);
	}
	
	
	
	
	// ------------------------- PROCEDURES BASED ON ALIGNMENTS --------------------------
	
	/**
	 * @return TRUE iff alignment $[alignmentStart1..alignmentEnd1] x
	 * [alignmentStart2..alignmentEnd2]$ maps the whole interval $[start1..end1]$ onto 
	 * the whole interval $[start2..end2]$, possibly in opposite orientation. Intervals 
	 * are assumed to either coincide with alignment intervals, or to be contained in 
	 * alignment intervals.
	 */
	public static final boolean areApproximatelyIdentical(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean orientation) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		final int DISTANCE_THRESHOLD = absoluteThreshold;  // Arbitrary
		final double RATIO = ((double)(alignmentEnd2-alignmentStart2+1))/(alignmentEnd1-alignmentStart1+1);
		boolean correct;
		int delta1, delta2;
		
		if ( areApproximatelyIdentical(alignmentStart1,alignmentEnd1,start1,end1) &&  
			 areApproximatelyIdentical(alignmentStart2,alignmentEnd2,start2,end2)
		   ) return true;
		if ( !isApproximatelyContained(start1,end1,alignmentStart1,alignmentEnd1) ||
			 !isApproximatelyContained(start2,end2,alignmentStart2,alignmentEnd2)
		   ) return false;
		correct=false;
		delta1=start1-alignmentStart1; delta2=orientation?start2-alignmentStart2:alignmentEnd2-end2;
		if ( (delta1>0 && delta2>0 && Math.abs((int)(delta1*RATIO),delta2)<=DISTANCE_THRESHOLD) ||
			 (delta1<0 && delta2<=DISTANCE_THRESHOLD) ||
			 (delta2<0 && delta1<=DISTANCE_THRESHOLD)
		   ) correct=true;
		if (!correct) return false;
		correct=false;
		delta1=alignmentEnd1-end1; delta2=orientation?alignmentEnd2-end2:start2-alignmentStart2;
		if ( (delta1>0 && delta2>0 && Math.abs((int)(delta1*RATIO),delta2)<=DISTANCE_THRESHOLD) ||
			 (delta1<0 && delta2<=DISTANCE_THRESHOLD) ||
			 (delta2<0 && delta1<=DISTANCE_THRESHOLD)
		   ) correct=true;
		return correct;
	}
	
	
	/**
	 * @param out output array: 0/1=first/last absolute position of the projection of 
	 * $[start1..end1]$ inside $[start2..end2]$;
	 * @return TRUE iff alignment $[alignmentStart1..alignmentEnd1] x
	 * [alignmentStart2..alignmentEnd2]$ maps the whole interval $[start1..end1]$ onto 
	 * a proper subinterval of the whole interval $[start2..end2]$, possibly in opposite 
	 * orientation. $[start2..end2]$ might not be contained in $[alignmentStart2..
	 * alignmentEnd2]$.
	 */
	public static final boolean isApproximatelyContained(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean orientation, int[] out) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		if (length1>=length2) return false;
		final int DISTANCE_THRESHOLD = absoluteThreshold;  // Arbitrary
		final double RATIO = ((double)(alignmentEnd2-alignmentStart2+1))/(alignmentEnd1-alignmentStart1+1);
		final int delta1 = start1-alignmentStart1;
		final int delta2 = alignmentEnd1-end1;
		boolean leftIdentical, leftAfter, rightIdentical, rightBefore;
		int outFirst, outLast;
		
		if (areApproximatelyIdentical(start1,end1,alignmentStart1,alignmentEnd1)) {
			if ( !areApproximatelyIdentical(start2,end2,alignmentStart2,alignmentEnd2) &&
		 	   	 isApproximatelyContained(start2,end2,alignmentStart2,alignmentEnd2)
	   		   ) return false;
		}
		else if (!isApproximatelyContained(start1,end1,alignmentStart1,alignmentEnd1)) return false;
		if (orientation) {
			outFirst=alignmentStart2+(int)(delta1*RATIO);
			outLast=alignmentEnd2-(int)(delta2*RATIO);
			leftIdentical=(delta1>=0||delta1>=-DISTANCE_THRESHOLD) && Math.abs(outFirst,start2)<=DISTANCE_THRESHOLD;
			leftAfter=(delta1>=0||delta1>=-DISTANCE_THRESHOLD) && outFirst>start2+DISTANCE_THRESHOLD;
			rightIdentical=(delta2>=0||delta2>=-DISTANCE_THRESHOLD) && Math.abs(outLast,end2)<=DISTANCE_THRESHOLD;
			rightBefore=(delta2>=0||delta2>=-DISTANCE_THRESHOLD) && outLast<end2-DISTANCE_THRESHOLD;
		}
		else {
			outFirst=alignmentStart2+(int)(delta2*RATIO);
			outLast=alignmentEnd2-(int)(delta1*RATIO);
			leftIdentical=(delta1>=0||delta1>=-DISTANCE_THRESHOLD) && Math.abs(outLast,end2)<=DISTANCE_THRESHOLD;
			leftAfter=(delta1>=0||delta1>=-DISTANCE_THRESHOLD) && outLast<end2-DISTANCE_THRESHOLD;
			rightIdentical=(delta2>=0||delta2>=-DISTANCE_THRESHOLD) && Math.abs(outFirst,start2)<=DISTANCE_THRESHOLD;
			rightBefore=(delta2>=0||delta2>=-DISTANCE_THRESHOLD) && outFirst>start2+DISTANCE_THRESHOLD;
		}
		if ((leftIdentical && rightBefore) || (leftAfter && rightIdentical) || (leftAfter && rightBefore)) {
			out[0]=outFirst; out[1]=outLast;
			return true;
		}
		return false;
	}

	
	/**
	 * @param prefixOrSuffix TRUE=prefix, FALSE=suffix;
	 * @return TRUE iff alignment $[alignmentStart1..alignmentEnd1] x
	 * [alignmentStart2..alignmentEnd2]$ maps the whole interval $[start1..end1]$ to a 
	 * proper prefix (respectively, suffix) of interval $[start2..end2]$, possibily in
	 * opposite orientation.
	 */
	public static final boolean isApproximatePrefixOrSuffix(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean orientation, boolean prefixOrSuffix) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		final int DISTANCE_THRESHOLD = absoluteThreshold;  // Arbitrary
		final double RATIO = ((double)(alignmentEnd2-alignmentStart2+1))/(alignmentEnd1-alignmentStart1+1);
		boolean correct;
		int delta1, delta2;
		
		if (areApproximatelyIdentical(alignmentStart1,alignmentEnd1,start1,end1)) return isApproximatePrefixOrSuffix(alignmentStart2,alignmentEnd2,start2,end2,prefixOrSuffix);
		if (!isApproximatelyContained(start1,end1,alignmentStart1,alignmentEnd1)) return false;
		correct=false;
		if (prefixOrSuffix) {
			delta1=orientation?start1-alignmentStart1:alignmentEnd1-end1;
			delta2=start2-alignmentStart2;
		}
		else {
			delta1=orientation?alignmentEnd1-end1:start1-alignmentStart1;
			delta2=alignmentEnd2-end2;
		}
		if ( (delta1>0 && delta2>0 && Math.abs((int)(delta1*RATIO),delta2)<=DISTANCE_THRESHOLD) ||
			 (delta1<0 && delta2<=DISTANCE_THRESHOLD) ||
			 (delta2<0 && delta1<=DISTANCE_THRESHOLD)
		   ) correct=true;
		if (!correct) return false;
		if (prefixOrSuffix) return end2-(alignmentStart2+(orientation?end1-alignmentStart1:alignmentEnd1-start1)*RATIO)>DISTANCE_THRESHOLD;
		else return (alignmentStart2+(orientation?start1-alignmentStart1:alignmentEnd1-end1)*RATIO)-start2>DISTANCE_THRESHOLD;
	}
	
	
	/**
	 * Let $[startX..endX]$ be two intervals in readA, and let $[alignmentXStartA..
	 * alignmentXEndA] x [alignmentXStartB..alignmentXEndB]$ be two alignments. The 
	 * procedure stores in $out[0,1]$ (respectively, in $out[2,3]$) an estimate of the 
	 * first and last position, in readA, of the portion of $[start1..end1]$ that belongs 
	 * to the intersection of the projections of $[start1..end1]$ by $alignment1$, and of 
	 * $[start2..end2]$ by $alignment2$, in readB. $out[4]$ stores the orientation 
	 * (1=normal, 0=RC) of the resulting artificial, same-read, alignment.
	 *
	 * Remark: the $diffs$ field of the artificial alignment is not estimated.
	 *
	 * Remark: the procedure modifies $out$ even when it terminates unsuccessfully.
	 *
	 * @return false iff any of the projections involved is empty.
	 */
	public static final boolean project(int start1, int end1, int start2, int end2, int alignment1StartA, int alignment1EndA, int alignment1StartB, int alignment1EndB, boolean alignment1Orientation, int alignment2StartA, int alignment2EndA, int alignment2StartB, int alignment2EndB, boolean alignment2Orientation, int[] out) {
		boolean success;
		final int intersectionStart, intersectionEnd;
		int projectionStart1, projectionEnd1, projectionStart2, projectionEnd2;
		
		success=project(start1,end1,alignment1StartA,alignment1EndA,alignment1StartB,alignment1EndB,alignment1Orientation,out,0);
		if (!success) return false;
		projectionStart1=out[0]; projectionEnd1=out[1];
		success=project(start2,end2,alignment2StartA,alignment2EndA,alignment2StartB,alignment2EndB,alignment2Orientation,out,0);
		if (!success) return false;
		projectionStart2=out[0]; projectionEnd2=out[1];
		if (intersectionLength(projectionStart1,projectionEnd1,projectionStart2,projectionEnd2)==0) return false;
		intersectionStart=Math.max(projectionStart1,projectionStart2);
		intersectionEnd=Math.min(projectionEnd1,projectionEnd2);
		success=project(intersectionStart,intersectionEnd,alignment1StartB,alignment1EndB,alignment1StartA,alignment1EndA,alignment1Orientation,out,0);
		if (!success) return false;
		success=project(intersectionStart,intersectionEnd,alignment2StartB,alignment2EndB,alignment2StartA,alignment2EndA,alignment2Orientation,out,2);
		if (!success) return false;
		out[4]=alignment1Orientation==alignment2Orientation?1:0;
		return true;
	}
	
	
	/**
	 * Let $[start..end]$ be an interval in readA, and let $[alignmentStartA..
	 * alignmentEndA] x [alignmentStartB..alignmentEndB]$ be an alignment that might or 
	 * might not intersect $[start..end]$ in readA. The procedure stores in $out[p..]$ an 
	 * estimate of the first and last position in readB of the portion of $[start..end]$ 
	 * that intersects the alignment in readA.
	 *
	 * @return false iff the projection is empty (e.g. because the alignment does not
	 * intersect the interval in readA). In this case $out$ is not modified.
	 */
	public static final boolean project(int start, int end, int alignmentStartA, int alignmentEndA, int alignmentStartB, int alignmentEndB, boolean orientation, int[] out, int p) {
		final double RATIO = ((double)(alignmentEndB-alignmentStartB+1))/(alignmentEndA-alignmentStartA+1);
		int projectionStart, projectionEnd;
		
		projectionStart=-1; projectionEnd=-1;
		if (intersectionLength(start,end,alignmentStartA,alignmentEndA)==0) return false;
		if (orientation) {
			if (start<alignmentStartA) projectionStart=alignmentStartB;
			else projectionStart=(int)(alignmentStartB+(start-alignmentStartA)*RATIO);
			if (end>alignmentEndA) projectionEnd=alignmentEndB;
			else projectionEnd=(int)(alignmentEndB-(alignmentEndA-end)*RATIO);
		}
		else {
			if (end>alignmentEndA) projectionStart=alignmentStartB;
			else projectionStart=(int)(alignmentStartB+(alignmentEndA-end)*RATIO);
			if (start<alignmentStartA) projectionEnd=alignmentEndB;
			else projectionEnd=(int)(alignmentEndB-(start-alignmentStartA)*RATIO);
		}
		if (projectionStart>=projectionEnd) return false;
		out[p]=projectionStart; out[p+1]=projectionEnd;
		return true;
	}
	
	
	/**
	 * @param intervalsX a sorted sequence of disjoint intervals; every interval is a pair
	 * $[start..end]$; $lastX$ is the last position of an $end$ element;
	 * @param minIntervalLength only intervals of length $>=minIntervalLength$ contribute
	 * to the output;
	 * @return if $mode=true$, the sum of all intersections, and if $mode=false$, the 
	 * length of a longest intersection, between an interval in $intervalsI$ and an 
	 * interval in $intervalsJ$.
	 */
	public static final int intersectionLength(int[] intervalsI, int lastI, int[] intervalsJ, int lastJ, int minIntervalLength, boolean mode) {
		int i, j;
		int firstJForNextI, intersection, out;
		
		if (IO.CONSISTENCY_CHECKS) {
			for (i=2; i<lastI; i+=2) {
				if (intervalsI[i]<=intervalsI[i-1]) {
					System.err.println("Intervals.intersectionLength> ERROR 1");
					System.err.print("intervalsI: ");
					for (int x=0; x<lastI; x+=2) System.err.print("["+intervalsI[x]+".."+intervalsI[x+1]+"], ");
					System.err.println();
					System.err.print("intervalsJ: ");
					for (int x=0; x<lastJ; x+=2) System.err.print("["+intervalsJ[x]+".."+intervalsJ[x+1]+"], ");
					System.err.println();
					System.exit(1);
				}
			}
			for (i=2; i<lastJ; i+=2) {
				if (intervalsJ[i]<=intervalsJ[i-1]) {
					System.err.println("Intervals.intersectionLength> ERROR 2");
					System.err.print("intervalsI: ");
					for (int x=0; x<lastI; x+=2) System.err.print("["+intervalsI[x]+".."+intervalsI[x+1]+"], ");
					System.err.println();
					System.err.print("intervalsJ: ");
					for (int x=0; x<lastJ; x+=2) System.err.print("["+intervalsJ[x]+".."+intervalsJ[x+1]+"], ");
					System.err.println();
					System.exit(1);
				}
			}
		}
		
		out=0; i=0; j=0; firstJForNextI=-1;
		while (i<lastI && j<lastJ) {
			if (intervalsJ[j+1]<intervalsI[i]) {
				j+=2;
				continue;
			}
			if (firstJForNextI==-1 && i<lastI-1 && intervalsJ[j+1]>=intervalsI[i+2]) firstJForNextI=j;
			if (intervalsJ[j]>intervalsI[i+1]) {
				i+=2;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			intersection=intersectionLength(intervalsI[i],intervalsI[i+1],intervalsJ[j],intervalsJ[j+1]);
			if (intersection>=minIntervalLength) {
				if (mode) out+=intersection;
				else out=Math.max(out,intersection);
			}
			j+=2;
			if (j>lastJ) {
				i+=2;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
			}
		}
		return out;
	}


	/**
	 * @param intervalsX a (possibly empty) sorted sequence of disjoint intervals; every 
	 * interval is a pair $[start..end]$; $lastX$ is the last position of an $end$ 
	 * element;
	 * @param id constant ID added to $out$ (see below);
	 * @param out output array (assumed to be large enough): starting at $from$ 
	 * (inclusive), contains the sorted set of all intervals of length at least 
	 * $minIntervalLength$ in $intervalsI \setminus intervalsJ$, where every interval is a
	 * triplet $id,start,end$;
	 * @return the last element in $out$ after the difference ($from-1$ if the difference
	 * is empty).
	 */
	public static final int difference(int[] intervalsI, int lastI, int[] intervalsJ, int lastJ, int[] out, int from, int id, int minIntervalLength) {
		int i, j;
		int firstJForNextI, gapStart, lastOut;
		
		if (IO.CONSISTENCY_CHECKS) {
			for (i=2; i<lastI; i+=2) {
				if (intervalsI[i]<=intervalsI[i-1]) {
					System.err.println("Intervals.difference> ERROR 1");
					System.exit(1);
				}
			}
			for (i=2; i<lastJ; i+=2) {
				if (intervalsJ[i]<=intervalsJ[i-1]) {
					System.err.println("Intervals.difference> ERROR 2");
					System.exit(1);
				}
			}
		}
		lastOut=from-1;
		if (lastI==-1) return lastOut;
		if (lastJ==-1) {
			for (i=0; i<lastI; i+=2) {
				out[++lastOut]=id;
				out[++lastOut]=intervalsI[i];
				out[++lastOut]=intervalsI[i+1];
			}
			return lastOut;
		}
		i=0; j=0; firstJForNextI=-1;
		gapStart=lastI>0?intervalsI[0]:-1;
		while (i<lastI) {
			if (j>lastJ) {
				if (gapStart<=intervalsI[i+1] && intervalsI[i+1]-gapStart+1>=minIntervalLength) {
					out[++lastOut]=id;
					out[++lastOut]=gapStart;
					out[++lastOut]=intervalsI[i+1];
				}
				i+=2;
				if (i<lastI) gapStart=intervalsI[i];
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			if (intervalsJ[j+1]<intervalsI[i]) {
				j+=2;
				continue;
			}
			if (firstJForNextI==-1 && i<lastI-1 && intervalsJ[j+1]>=intervalsI[i+2]) firstJForNextI=j;
			if (intervalsJ[j]>intervalsI[i+1]) {
				if (gapStart<=intervalsI[i+1] && intervalsI[i+1]-gapStart+1>=minIntervalLength) {
					out[++lastOut]=id;
					out[++lastOut]=gapStart;
					out[++lastOut]=intervalsI[i+1];
				}
				i+=2; 
				if (i<lastI) gapStart=intervalsI[i];
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			if (intervalsJ[j]<=intervalsI[i]) gapStart=intervalsJ[j+1]+1;
			else if (intervalsJ[j]-gapStart>=minIntervalLength) {
				out[++lastOut]=id;
				out[++lastOut]=gapStart;
				out[++lastOut]=intervalsJ[j]-1;
				gapStart=intervalsJ[j+1]+1;
			}
			else gapStart=intervalsJ[j+1]+1;
			j+=2;
		}
		return lastOut;
	}

	
	/**
	 * @param intervalsX a (possibly empty) sorted sequence of disjoint intervals; every 
	 * interval is a pair $[start..end]$; $lastX$ is the last position of an $end$ 
	 * element;
	 * @param out the union of $intervalsI$ and $intervalsJ$; the array is assumed to be
	 * large enough;
	 * @return the last element in $out$.
	 */
	public static final int union(int[] intervalsI, int lastI, int[] intervalsJ, int lastJ, int[] out) {
		int i, j, k;
		int last, currentStart, currentEnd;
		if (lastI==-1 && lastJ==-1) return -1;
		
		if (IO.CONSISTENCY_CHECKS) {
			for (i=2; i<lastI; i+=2) {
				if (intervalsI[i]<=intervalsI[i-1]) {
					System.err.println("Intervals.union> ERROR 1");
					System.err.print("intervalsI: ");
					for (int x=0; x<lastI; x+=2) System.err.print("["+intervalsI[x]+".."+intervalsI[x+1]+"], ");
					System.err.println();
					System.err.print("intervalsJ: ");
					for (int x=0; x<lastJ; x+=2) System.err.print("["+intervalsJ[x]+".."+intervalsJ[x+1]+"], ");
					System.err.println();
					System.exit(1);
				}
			}
			for (i=2; i<lastJ; i+=2) {
				if (intervalsJ[i]<=intervalsJ[i-1]) {
					System.err.println("Intervals.union> ERROR 2");
					System.err.print("intervalsI: ");
					for (int x=0; x<lastI; x+=2) System.err.print("["+intervalsI[x]+".."+intervalsI[x+1]+"], ");
					System.err.println();
					System.err.print("intervalsJ: ");
					for (int x=0; x<lastJ; x+=2) System.err.print("["+intervalsJ[x]+".."+intervalsJ[x+1]+"], ");
					System.err.println();
					System.exit(1);
				}
			}
		}
		
		// Merge-sorting the intervals by start position
		i=0; j=0; k=-1;
		while (i<lastI && j<lastJ) {
			if (intervalsI[i]<intervalsJ[j]) {
				out[++k]=intervalsI[i++]; out[++k]=intervalsI[i++];
			}
			else if (intervalsJ[j]<intervalsI[i]) {
				out[++k]=intervalsJ[j++]; out[++k]=intervalsJ[j++];
			}
			else {
				out[++k]=intervalsI[i++]; out[++k]=intervalsI[i++];
				out[++k]=intervalsJ[j++]; out[++k]=intervalsJ[j++];
			}
		}
		while (i<lastI) { out[++k]=intervalsI[i++]; out[++k]=intervalsI[i++]; }
		while (j<lastJ) { out[++k]=intervalsJ[j++]; out[++k]=intervalsJ[j++]; }
		last=k;
		if (last==-1) return -1;
		
		// Merging overlapping intervals
		k=-1; currentStart=out[0]; currentEnd=out[1];
		for (i=2; i<last; i+=2) {
			if (out[i]>currentEnd) { 
				out[++k]=currentStart; out[++k]=currentEnd; 
				currentStart=out[i]; currentEnd=out[i+1];
			}
			else currentEnd=Math.max(currentEnd,out[i+1]);
		}
		out[++k]=currentStart; out[++k]=currentEnd;
		
		return k;
	}
	
	
	/**
	 * @param intervals a (possibly empty) sorted sequence of disjoint intervals; every 
	 * interval is a pair $[start..end]$; $last$ is the last position of an $end$ element;
	 * @param length length of the segment on which intervals are assumed to lie;
	 * @param out output array, storing the complement of $intervals$ in the segment; the 
	 * array is assumed to be large enough;
	 * @return the last element in $out$.
	 */
	public static final int complement(int[] intervals, int last, int length, int[] out) {
		int i, j;
		
		if (last==-1) {
			out[0]=0; out[1]=length-1;
			return 1;
		}
		j=-1;
		if (intervals[0]>0) {
			out[++j]=0; out[++j]=intervals[0]-1;
		}
		for (i=2; i<last; i+=2) {
			out[++j]=intervals[i-1]+1;
			out[++j]=intervals[i]-1;
		}
		if (intervals[last]<length-1) {
			out[++j]=intervals[last]+1; out[++j]=length-1;
		}
		return j;
	}
	
	
	

}