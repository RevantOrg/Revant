package de.mpi_cbg.revant.util;

import java.util.Arrays;
import de.mpi_cbg.revant.util.Math;


public class Points {

	public static final int FEW_POINTS = 10;  // Some procedures behave differently when points are few.
	public static final int DEFAULT_NBINS = 10;  // Default number of bins to be used in $areUniformlyDistributed()$ when no prior knowledge of bin length exists.
	private static double[] distances;
	private static int lastDistance;
	private static Point tmpPoint;
	public static Point[] tmpPoints;
	private static int lastTmpPoint;


	public static final void allocateMemory(int maxPoints) {
		if (distances==null || distances.length<maxPoints) distances = new double[maxPoints];
		if (tmpPoint==null) tmpPoint = new Point();
		allocateTmpPoints(maxPoints);
	}
	
	
	public static final void allocateTmpPoints(int maxPoints) {
		if (tmpPoints==null || tmpPoints.length<maxPoints) {
			tmpPoints = new Point[maxPoints];
			for (int i=0; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
		}
	}
	
	
	public static final void deallocateMemory() {
		distances=null;
		tmpPoint=null;
		for (int i=0; i<tmpPoints.length; i++) tmpPoints[i]=null;
		tmpPoints=null;
	}


	/**
	 * Sorts $points$ by $position$, and transforms runs of points with the same position
	 * into a single point, whose counts are the sum of the counts of all points in the
	 * run.
	 *
	 * @return the value of $lastPoint$ after compaction.
	 */
	public static final int sortAndCompact(Point[] points, int lastPoint) {
		int i, j;
		
		if (lastPoint<=0) return lastPoint;
		Point.order=Point.POSITION;
		if (lastPoint>0) Arrays.sort(points,0,lastPoint+1);
		i=0; j=i+1;
		while (j<=lastPoint) {
			if (points[j].position==points[i].position) points[i].sum(points[j]);
			else points[++i].copyFrom(points[j]);
			j++;
		}
		return i;
	}
	
	
	/**
	 * Merges two sorted and compacted arrays of points into a sorted and compacted 
	 * destination array of points.
	 *
	 * @return the last element of $destination$.
	 */
	public static final int merge(Point[] points1, int lastPoint1, Point[] points2, int lastPoint2, Point[] destination) {
		int i1, i2, j;
		
		
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("merging the following lists: ");
	for (int x=0; x<=lastPoint1; x++) IO.printErr(points1[x].position+" ");	
	IO.printErr("\n");
	for (int x=0; x<=lastPoint2; x++) IO.printErr(points2[x].position+" ");	
	IO.printErr("\n");
}		
		
		
		i1=0; i2=0; j=-1;
		while (i1<=lastPoint1 && i2<=lastPoint2) {
			if (points1[i1].position<points2[i2].position) {
				j++;
				destination[j].position=points1[i1].position;
				destination[j].mass=points1[i1].mass;
				i1++;
			}
			else if (points2[i2].position<points1[i1].position) {
				j++;
				destination[j].position=points2[i2].position;
				destination[j].mass=points2[i2].mass;
				i2++;
			}
			else {
				j++;
				destination[j].position=points1[i1].position;
				destination[j].mass=points1[i1].mass+points2[i2].mass;
				if (points1[i1].mass>0 && points2[i2].mass>0 && destination[j].mass<0) destination[j].mass=Integer.MAX_VALUE;  // Overflow
				i1++; i2++;
			}
		}
		while (i1<=lastPoint1) {
			j++;
			destination[j].position=points1[i1].position;
			destination[j].mass=points1[i1].mass;
			i1++;
		}
		while (i2<=lastPoint2) {
			j++;
			destination[j].position=points2[i2].position;
			destination[j].mass=points2[i2].mass;
			i2++;
		}
		

if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("result of the merge: ");
	for (int x=0; x<=j; x++) IO.printErr(destination[x].position+" ");		
}
		
		return j;
	}
	
	
	/**
	 * Copies the values of $source[0..last]$ into $destination[0..last]$, assuming that 
	 * the latter is long enough.
	 *
	 * @return the input value of $last$.
	 */
	public static final int simpleClone(Point[] source, int last, Point[] destination) {
		for (int i=0; i<=last; i++) destination[i].copyFrom(source[i]);
		return last;
	}
	
	
	/**
	 * Copies the values of the points in $source[from..to]$ (which are assumed to be 
	 * sorted and distinct) into the points at the beginning of $destination$. Since 
	 * $destination$ can be significantly smaller than $source$, the procedure greedily 
	 * collapses into a single point a maximal run of consecutive points that are at 
	 * distance at most $d$ from the first point of the run. The procedure uses the 
	 * smallest value of $d$ (possibly zero) that makes $source$ fit into $destination$.
	 *
	 * @return the last element in $destination$.
	 */
	public static final int clone(Point[] source, int from, int to, Point[] destination) {
		int i, j, d;
		
		d=-1;
		do {
			d++;
			tmpPoints[0].position=source[from].position;
			tmpPoints[0].mass=source[from].mass;
			j=0;
			for (i=from+1; i<=to; i++) {
				if (source[i].position<=tmpPoints[j].position+d) tmpPoints[j].mass+=source[i].mass;
				else {
					j++;
					tmpPoints[j].position=source[i].position;
					tmpPoints[j].mass=source[i].mass;
				}
			}
		} while (j>=destination.length);
		for (i=0; i<=j; i++) {
			destination[i].position=tmpPoints[i].position;
			destination[i].mass=tmpPoints[i].mass;
		}
		return j;
	}


	/**
	 * Returns the smallest distance $d$ between the first and the last of $nPoints$
	 * consecutive points in $points$, such that a $quantile$ fraction of such distances
	 * are smaller than or equal to $d$, where $quantile$ is in $(0..1]$.
	 *
	 * @param integerPositions TRUE=$points$ stores integer positions;
	 * @return -1 if $lastPoint+1$ is smaller than $nPoints$.
	 */
	public static final double distanceQuantile(Point[] points, int lastPoint, int nPoints, boolean integerPositions, double quantile) {		
		if (lastPoint+1<nPoints) return -1;
		if (distances==null || distances.length<lastPoint+1) distances = new double[lastPoint+1];
		lastDistance=-1;
		for (int i=nPoints-1; i<=lastPoint; i++) distances[++lastDistance]=points[i].position-points[i-nPoints+1].position+(integerPositions?1:0);
		if (lastDistance==0) return distances[0];
		if (lastDistance>0) Arrays.sort(distances,0,lastDistance+1);
		if (((double)lastDistance)/(lastDistance+1)<quantile) return distances[lastDistance-1];
		return distances[(int)Math.ceil(quantile*(lastDistance+1))-1];
	}


	/**
	 * Estimates a rough splitting of $points$, discarding points at position smaller than
	 * $smallestPosition$. If the maximum mass of a point is large enough WRT the minimum
	 * mass of a point, the procedure returns the point immediately to the right (if
	 * $side=true$) or to the left (if $side=false$) of the leftmost (if $direction=true$)
	 * or of the rightmost (if $direction=false$) point with maximum mass. Otherwise, the 
	 * procedure returns a point $i$ with large distance between $points[i]$ and 
	 * $points[i-1]$, where "large" is defined WRT the smallest distance between 
	 * consecutive points: the selected point is the leftmost if $direction=true$, the 
	 * rightmost otherwise; the procedure returns $i$ if $side=true$, $i-1$ otherwise.
	 *
	 * Remark: the procedure assumes $points$ to be sorted and compacted.
	 *
	 * @return -1 if $lastPoint=0$.
	 */
	public static final int getRoughThreshold(Point[] points, int lastPoint, boolean direction, int smallestPosition, boolean side) {
		final double MIN_MASS_RATIO = 8.0;
		final double MIN_DISTANCE_RATIO = 5.0;
		final double QUANTILE = 0.5;
		final double MAX_DISTANCE_THRESHOLD = 0.75;
		int i, max, min, firstPoint;
		int largestDistancePoint, smallestDistancePoint;

		firstPoint=0;
		while (firstPoint<=lastPoint && points[firstPoint].position<smallestPosition) firstPoint++;
		if (lastPoint-firstPoint<=0) return -1;
		if (lastPoint-firstPoint==1) {
			if (firstPoint==0) return lastPoint;
			if (points[lastPoint].position-points[lastPoint-1].position>points[firstPoint].position-points[firstPoint-1].position) return side?lastPoint:firstPoint;
			else return side?firstPoint:firstPoint-1;
		}
		min=Math.POSITIVE_INFINITY; max=0;
		for (i=firstPoint; i<=lastPoint; i++) {
			if (points[i].getMass()<min) min=points[i].getMass();
			if (points[i].getMass()>max) max=points[i].getMass();
		}
		if (max>=MIN_MASS_RATIO*min) {
			if (direction) {
				for (i=firstPoint; i<=lastPoint; i++) {
					if (points[i].getMass()==max) {
						if (side) return i==lastPoint?i:i+1;
						else return i==0?i:i-1;
					}
				}
			}
			else {
				for (i=lastPoint; i>=firstPoint; i--) {
					if (points[i].getMass()==max) {
						if (side) return i==lastPoint?i:i+1;
						else return i==0?i:i-1;
					}
				}
			}
		}
		else {
			largestDistancePoint=largestDistance(points,firstPoint,lastPoint,direction);
			smallestDistancePoint=smallestDistance(points,firstPoint,lastPoint,direction);
if (smallestDistancePoint!=0) if (IO.SHOW_STD_ERR_PRIME) IO.printErr("smallestDistance="+(points[smallestDistancePoint].position-points[smallestDistancePoint-1].position));
if (largestDistancePoint!=0) if (IO.SHOW_STD_ERR_PRIME) IO.printErr("largestDistance="+(points[largestDistancePoint].position-points[largestDistancePoint-1].position));
			if ( largestDistancePoint!=0 && smallestDistancePoint!=0 &&
			     points[largestDistancePoint].position-points[largestDistancePoint-1].position>=MIN_DISTANCE_RATIO*(points[smallestDistancePoint].position-points[smallestDistancePoint-1].position)
			   ) {
 if (IO.SHOW_STD_ERR_PRIME) IO.printErr("---> 1 smallestPosition="+smallestPosition);
				   return atDistance(points,firstPoint,lastPoint,(points[largestDistancePoint].position-points[largestDistancePoint-1].position)*MAX_DISTANCE_THRESHOLD,direction,side);
			   }
		}
 if (IO.SHOW_STD_ERR_PRIME) IO.printErr("---> 2 smallestPosition="+smallestPosition);
		return (int)quantile(points,firstPoint,lastPoint,direction,QUANTILE,true);
	}
	
	
	/**
	 * Returns a point $i$ such that the distance between $points[i]$ and $points[i-1]$ is
	 * at least $minDistance$, or $lastPoint$ if $lastPoint=firstPoint$. Such point is the
	 * leftmost if $direction=true$, the rightmost otherwise; the procedure returns $i$ if
	 * $side=true$, $i-1$ otherwise.
	 *
	 * Remark: the procedure assumes $points$ to be sorted.
	 */
	public static final int atDistance(Point[] points, int firstPoint, int lastPoint, double minDistance, boolean direction, boolean side) {
		int i, from, largestPoint;
		double distance, largestDistance;

		if (lastPoint-firstPoint==0) return lastPoint;
		from=firstPoint==0?1:firstPoint;
		largestDistance=0; largestPoint=0;
		if (direction) {
			for (i=from; i<=lastPoint; i++) {
				distance=points[i].position-points[i-1].position;
				if (distance>=minDistance) {
					largestPoint=i;
					break;
				}
			}
		}
		else {
			for (i=lastPoint; i>=from; i--) {
				distance=points[i].position-points[i-1].position;
				if (distance>=minDistance) {
					largestPoint=i;
					break;
				}
			}
		}
 if (IO.SHOW_STD_ERR_PRIME) IO.printErr("---> 3 largestPoint="+largestPoint);		
	 	return side?largestPoint:largestPoint-1;
	}


	/**
	 * Returns a point $i$ that maximizes the distance between $points[i]$ and
	 * $points[i-1]$, or $lastPoint$ if $lastPoint=firstPoint$. Such point is the leftmost
	 * if $direction=true$, the rightmost otherwise.
	 *
	 * Remark: the procedure assumes $points$ to be sorted.
	 */
	public static final int largestDistance(Point[] points, int firstPoint, int lastPoint, boolean direction) {
		int i, from, largestPoint;
		double distance, largestDistance;

		if (lastPoint-firstPoint==0) return lastPoint;
		from=firstPoint==0?1:firstPoint;
		largestDistance=0; largestPoint=0;
		if (direction) {
			for (i=from; i<=lastPoint; i++) {
				distance=points[i].position-points[i-1].position;
				if (distance>largestDistance) {
					largestDistance=distance;
					largestPoint=i;
				}
			}
		}
		else {
			for (i=lastPoint; i>=from; i--) {
				distance=points[i].position-points[i-1].position;
				if (distance>largestDistance) {
					largestDistance=distance;
					largestPoint=i;
				}
			}
		}
		return largestPoint;
	}


	/**
	 * Returns a point $i$ that minimizes the distance between $points[i]$ and
	 * $points[i-1]$, or zero if $lastPoint=0$. Such point is the leftmost if
	 * $direction=true$, the rightmost otherwise.
	 *
	 * Remark: the procedure assumes $points$ to be sorted.
	 */
	public static final int smallestDistance(Point[] points, int firstPoint, int lastPoint, boolean direction) {
		int i, smallestPoint;
		double distance, smallestDistance;

		if (lastPoint-firstPoint==0) return 0;
		smallestDistance=Math.POSITIVE_INFINITY; smallestPoint=0;
		if (direction) {
			for (i=firstPoint+1; i<=lastPoint; i++) {
				distance=points[i].position-points[i-1].position;
				if (distance<smallestDistance) {
					smallestDistance=distance;
					smallestPoint=i;
				}
			}
		}
		else {
			for (i=lastPoint; i>firstPoint; i--) {
				distance=points[i].position-points[i-1].position;
				if (distance<smallestDistance) {
					smallestDistance=distance;
					smallestPoint=i;
				}
			}
		}
		return smallestPoint;
	}


	/**
	 * @param threshold if positive, uses also points at distance at most $threshold$ from
	 * from $first$ and $last$;
	 * @param lastPoint last element of $points$; used only when $threshold>0$;
	 * @return the center of mass of $points[first..last]$. If $observed=true$, returns
	 * the position field of the point in $points[first..last]$ that is closest to the 
	 * center of mass, if such closest point is close enough to the center of mass.
	 */
	public static final double getCenterOfMass(Point[] points, int first, int last, boolean observed, double threshold, int lastPoint) {
		final double DISTANCE_THRESHOLD = 0.25;
		int i, mass, minPoint;
		double centerOfMass, distance, minDistance, p;

		// Computing the center of mass
		centerOfMass=0.0; mass=0;
		for (i=first; i<=last; i++) {
			centerOfMass+=points[i].position*points[i].getMass();
			mass+=points[i].getMass();
		}
		if (threshold>0) {
			p=points[first].position;
			i=first-1;
			while (i>=0 && points[i].position>=p-threshold) {
				centerOfMass+=points[i].position*points[i].getMass();
				mass+=points[i].getMass();
				i--;
			}
			p=points[last].position;
			i=last+1;
			while (i<=lastPoint && points[i].position<=p+threshold) {
				centerOfMass+=points[i].position*points[i].getMass();
				mass+=points[i].getMass();
				i++;
			}
		}
		centerOfMass/=mass;
		if (!observed) return centerOfMass;

		// Finding the observed point that is closest to the center of mass
		minDistance=Math.POSITIVE_INFINITY; minPoint=-1;
		for (i=first; i<=last; i++) {
			distance=points[i].position-centerOfMass;
			if (distance<0) distance=-distance;
			if (distance<minDistance) {
				minDistance=distance;
				minPoint=i;
			}
		}
		return minDistance>(points[last].position-points[first].position)*DISTANCE_THRESHOLD?centerOfMass:points[minPoint].position;
	}


	/**
	 * Remark: the procedure assumes that $sortAndCompact$ has already been run on
	 * $points$.
	 */
	public static final int binarySearch(Point[] points, int lastPoint, double position) {
		tmpPoint.position=position;
		return Arrays.binarySearch(points,0,lastPoint+1,tmpPoint);
	}


	/**
	 * Tries to detect points with an anomalously high mass. This is done by fitting a
	 * DET on the set of all masses of $points$, and by returning the first value of a
	 * selected run of local-maximum leaves, if there are at least two runs of local-
	 * maximum leaves. To be considered "high", the mass of a point must be at least equal
	 * to a given absolute threshold.
	 *
	 * Remark: the procedure assumes that $sortAndCompact$ has already been run on
	 * $points$.
	 *
	 * @param maxDifference if there are "few" distinct masses, and if such masses differ
	 * by at most this quantity, the procedure assumes that they are uniformly
	 * distributed and it does not try to split them;
	 * @return -1 if no point with anomalously high mass can be detected.
	 */
	public static final int largestMass(Point[] points, int lastPoint, double maxDifference) {
		final int MIN_INTERVAL_LENGTH = 3;
		final int MIN_LOCAL_MAX_DISTANCE = 10;
		final int MIN_MASS = 3;
		final int MIN_MASS_FOR_HIGH = 5;  // Minimum mass for a point to be considered high
		final double SECOND_SELECTED_RUN_RATIO = 2.0;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		int i;
		int nLocalMaximumLeaves, selectedRun, currentRun;

		// Collecting masses
		allocateTmpPoints(lastPoint+1);
		lastTmpPoint=-1;
		for (i=0; i<=lastPoint; i++) {
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=points[i].getMass();
			tmpPoints[lastTmpPoint].mass=1;
		}
		if (lastTmpPoint<=0) return -1;
		lastTmpPoint=sortAndCompact(tmpPoints,lastTmpPoint);
		if (lastTmpPoint==0 || tmpPoints[lastTmpPoint].position<MIN_MASS_FOR_HIGH) return -1;

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Points.largestMass> points after sortAndCompact:");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastTmpPoint; i++) IO.printErr(tmpPoints[i]); }

		if (lastTmpPoint+1<=FEW_POINTS) {
			if (tmpPoints[lastTmpPoint].position-tmpPoints[0].position<=maxDifference) return -1;
			i=getRoughThreshold(tmpPoints,lastTmpPoint,false,MIN_MASS,true);
			if (i==-1) return -1;
			return (int)Math.max(tmpPoints[i].position,MIN_MASS_FOR_HIGH);
		}

		// Finding the selected run of local-maximum leaves
		if (areUniformlyDistributed(tmpPoints,0,lastTmpPoint,true,(tmpPoints[lastTmpPoint].position-tmpPoints[0].position)/DEFAULT_NBINS)) return -1;
		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(tmpPoints,0,lastTmpPoint,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE,true,-1,-1,-1,false,true);


if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Points.largestMass> DET: nLocalMaximumLeaves="+nLocalMaximumLeaves);
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=DensityEstimationTree.lastLeaf; i++) IO.printErr(DensityEstimationTree.leaves[i]); }


		if (nLocalMaximumLeaves==0) return -1;
		else if (nLocalMaximumLeaves==1) {
			DensityEstimationTree.markRunsOfLocalMaximumLeaves(tmpPoints);
			return (int)DensityEstimationTree.separateRun(0,tmpPoints,0,lastTmpPoint,true);
		}
		DensityEstimationTree.markRunsOfLocalMaximumLeaves(tmpPoints);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Points.largestMass> nRuns="+DensityEstimationTree.nRuns);

		if (DensityEstimationTree.nRuns==1) return (int)DensityEstimationTree.separateRun(0,tmpPoints,0,lastTmpPoint,true);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Points.largestMass> DET on masses:");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=DensityEstimationTree.lastLeaf; i++) IO.printErr(DensityEstimationTree.leaves[i]); }

		selectedRun=DensityEstimationTree.nRuns==2?1:DensityEstimationTree.getRunWithLargestDistance(tmpPoints,false);
		if (selectedRun==1) {
			i=Leaves.lastLeafOfRun(DensityEstimationTree.leaves,DensityEstimationTree.lastLeaf,0);
			if (tmpPoints[DensityEstimationTree.leaves[i+1].firstPoint].position<SECOND_SELECTED_RUN_RATIO*tmpPoints[DensityEstimationTree.leaves[i].lastPoint].position) {
				if (DensityEstimationTree.nRuns==2) return -1;
				selectedRun++;
			}
		}
		i=nLocalMaximumLeaves-1; currentRun=DensityEstimationTree.nRuns;
		while (i>=0) {
			if (DensityEstimationTree.leaves[i].marked) {
				currentRun--;
				if (currentRun==selectedRun) break;
			}
			i--;
		}
		i--;
		while (i>=0 && !DensityEstimationTree.leaves[i].marked) i--;
		i++;
		return (int)Math.max(tmpPoints[DensityEstimationTree.leaves[i].firstPoint].position,MIN_MASS_FOR_HIGH);
	}


	/**
	 * Uses Pearson's chi-squared test to decide whether $points[first..last]$, with
	 * corresponding masses, is distributed uniformly on a set of equal bins of length
	 * $binLength$ that divide the length $points[last].position-points[first].position$.
	 *
	 * Remark: $binLength$ is reset if the resulting number of bins is bigger than
	 * $Math.CHI_SQUARE.length+1$.
	 *
	 * Remark: the procedure assumes that procedure $sortAndCompact$ has already been run
	 * on $points$.
	 */
	public static final boolean areUniformlyDistributed(Point[] points, int first, int last, boolean integerPositions, double binLength) {
		int j;
		int nPoints, nBins, nDegreesOfFreedom;
		double i, length, expectation, sum, chiSquare;

		nPoints=0;
		for (j=first; j<=last; j++) nPoints+=points[j].getMass();
		length=points[last].position-points[first].position+(integerPositions?1:0);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("areUniformlyDistributed> length="+length+" binLength="+binLength);

		nBins=(int)Math.ceil(length/binLength);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("areUniformlyDistributed> nBins="+nBins);

		nDegreesOfFreedom=nBins-1;
		if (nDegreesOfFreedom>Math.CHI_SQUARE.length) {
			nBins=Math.CHI_SQUARE.length;
			nDegreesOfFreedom=nBins-1;
			binLength=integerPositions?Math.ceil((int)length,nBins):length/nBins;
			if (binLength==0.0) return false;
		}
		else if (nBins<DEFAULT_NBINS) {
			nBins=DEFAULT_NBINS;
			nDegreesOfFreedom=nBins-1;
			binLength=integerPositions?Math.ceil((int)length,nBins):length/nBins;
			if (binLength==0.0) return false;
		}
		expectation=((double)nPoints)/nBins;
		chiSquare=0.0; sum=0;
		i=points[first].position+binLength; j=first;
		while (j<=last) {
			if (points[j].position<i) {
				sum+=points[j].getMass();
				j++;
			}
			else {
				chiSquare+=(sum-expectation)*(sum-expectation)/expectation;
				i+=binLength;
				sum=0;
			}
		}
		chiSquare+=(sum-expectation)*(sum-expectation)/expectation;
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("areUniformlyDistributed> chiSquare="+chiSquare+" 99.9%="+Math.CHI_SQUARE[nDegreesOfFreedom-1][4]+" nBins="+nBins+" binLength="+binLength);
		return chiSquare<=Math.CHI_SQUARE[nDegreesOfFreedom-1][4];
	}
	
	
	/**
	 * Uses the Kolmogorov-Smirnov test to decide whether $points[first..last]$, with
	 * corresponding masses, comes from a continuous uniform distribution, with 
	 * significance level $significance$.
	 *
	 * Remark: many data points might be required in practice to reject uniformity.
	 *
	 * Remark: this procedure is currently not used anywhere in the code.
	 *
	 * Remark: the procedure assumes that procedure $sortAndCompact$ has already been run
	 * on $points$.
	 *
	 * @param significance significance id: 0=0.1, 1=0.05, 2=0.02, 3=0.01.
	 */
	public static final boolean areUniformlyDistributed_KS(Point[] points, int first, int last, boolean integerPositions, int significance) {
		int i, totalMass;
		double length, currentMass, difference, maxDifference;
		
		// Computing KS statistic
		totalMass=0;
		for (i=first; i<=last; i++) totalMass+=points[i].mass;
		length=points[last].position-points[first].position+(integerPositions?1:0);
		maxDifference=0.0; currentMass=0.0;
		for (i=first; i<=last; i++) {
			currentMass+=points[i].mass;
			difference=Math.abs(currentMass/totalMass-(points[i].position-points[first].position+(integerPositions?1:0))/length);
			if (difference>maxDifference) maxDifference=difference;
		}
		
		// Deciding confidence
		if (totalMass<=Math.KOLMOGOROV_SMIRNOV_ABSOLUTE.length) return maxDifference<=Math.KOLMOGOROV_SMIRNOV_ABSOLUTE[totalMass-1][significance];
		else return maxDifference<=Math.KOLMOGOROV_SMIRNOV_ABSOLUTE_LARGE[significance]/Math.sqrt(totalMass);
	}


	/**
	 * Checks whether there is an interval $[x..y]$ in $points$ such that the distance
	 * between any two consecutive points inside the interval is at most $threshold$,
	 * and such that the distance between $x$ and $y$ is at least $runLength$.
	 */
	public static final boolean hasRun(Point[] points, int last, double threshold, double runLength, boolean integerPositions) {
		int i;
		int firstPointInRun;

		firstPointInRun=0;
		for (i=1; i<=last; i++) {
			if (points[i].position-points[i-1].position>threshold) {
				if (points[i-1].position-points[firstPointInRun].position+(integerPositions?1:0)>=runLength) {
if (IO.SHOW_STD_ERR) IO.printErr("hasRun> ["+points[firstPointInRun].position+".."+points[i-1].position+"]");
					return true;
				}
				firstPointInRun=i;
			}
		}
		return false;
	}


	/**
	 * Stores in $histogram$ an estimate of the density of $points$ using rectangular 
	 * kernels. The bandwidth is computed by procedure $kernelBandwidth$, and it is 
	 * returned in output.
	 *
	 * Remark: $points$ is assumed to have integer positions, and to be sorted and 
	 * compacted.
	 *
	 * @param minimumBandwidth the procedure uses this if $kernelBandwidth()$ returns a
	 * smaller value.
	 */
	public static final int estimateDensity(Point[] points, int lastPoint, double[] histogram, int minimumBandwidth, int maximumBandwidth) {
		int i, j, from, to, bandwidth;
		double interquartileRange, mass;
		
		Math.set(histogram,histogram.length-1,0.0);
		if (lastPoint==0) {
			if (points[0].position<histogram.length) histogram[(int)points[0].position]=points[0].getMass();
			return 0;
		}
		interquartileRange=quantile(points,0,lastPoint,false,0.25,false)-quantile(points,0,lastPoint,true,0.25,false);
		bandwidth=(int)kernelBandwidth(lastPoint+1,Math.sqrt(variance(points,lastPoint)),interquartileRange);
		if (bandwidth<minimumBandwidth && (points[lastPoint].position-points[0].position+1)>minimumBandwidth*3) bandwidth=minimumBandwidth;
		if (bandwidth>maximumBandwidth) bandwidth=maximumBandwidth;
		for (i=0; i<=lastPoint; i++) {
			from=Math.max((int)points[i].position-bandwidth,0);
			to=Math.min((int)points[i].position+bandwidth,histogram.length-1);
			mass=points[i].getMass();
			for (j=from; j<=to; j++) histogram[j]+=mass;
		}
		return bandwidth;
	}
	
	
	/**
	 * @param points must be sorted and compacted;
	 * @param indexOrPosition TRUE=index, FALSE=position.
	 */
	public static final double quantile(Point[] points, int firstPoint, int lastPoint, boolean leftToRight, double quantile, boolean indexOrPosition) {
		int i, m;
		double mass;
		
		mass=mass(points,firstPoint,lastPoint)*quantile;
		m=0;
		if (leftToRight) { for (i=firstPoint; i<=lastPoint && m<mass; i++) m+=points[i].getMass(); }
		else { for (i=lastPoint; i>=firstPoint && m<mass; i--) m+=points[i].getMass(); }
		if (i<firstPoint) i=firstPoint;
		if (i>lastPoint) i=lastPoint-1;
		return indexOrPosition?i:points[i].position;
	}
	
	
	/**
	 * @param points does not need to be sorted or compacted.
	 */
	public static final int mass(Point[] points, int firstPoint, int lastPoint) {
		int i, mass;
		
		mass=0;
		for (i=firstPoint; i<=lastPoint; i++) mass+=points[i].getMass();
		return mass;
	}
	
	
	/**
	 * @param from,to positions (rather than indexes in $points$);
	 * @param points[0..lastPoint] assumed to be sorted and compacted.
	 */
	public static final int mass(Point[] points, int lastPoint, double from, double to) {
		int i, mass;
		
		tmpPoint.position=from;
		i=Arrays.binarySearch(points,0,lastPoint+1,tmpPoint);
		if (i<0) i=-i-1;
		mass=0;
		while (i<=lastPoint && points[i].position<=to) {
			mass+=points[i].getMass();
			i++;
		}
		return mass;
	}
	
	
	/**
	 * Estimates the bandwidth of a kernel density estimator, using Silverman's "rule of 
	 * thumb" for the Gaussian kernel. This is the 
	 * <a href="http://stat.ethz.ch/R-manual/R-patched/library/stats/html/density.html">
	 * default bandwidth used by R</a>. See also 
	 * <a href="http://stat.ethz.ch/R-manual/R-patched/library/stats/html/bandwidth.html">
	 * bw</a>. The procedure returns a fraction of this bandwidth.
	 *
	 * @param interquartileRange difference between the third quartile and the first quartile.
	 */
	private static final double kernelBandwidth(int nPoints, double standardDeviation, double interquartileRange) {
		final double ADJUST = 1.0/3;
		final double x = standardDeviation;
		final double y = interquartileRange/1.34;	
		return ADJUST*0.9*(y>0.0&&y<x?y:x)/Math.pow(nPoints,1.0/5);
	}
	
	
	/**
	 * @param points does not need to be sorted or compacted.
	 */
	public static final double variance(Point[] points, int lastPoint) {
		int i, mass;
		double delta, expectation, variance;
		
		expectation=expectation(points,lastPoint);
		variance=0.0; mass=0;
		for (i=0; i<=lastPoint; i++) {
			delta=points[i].position-expectation;
			variance+=delta*delta*points[i].getMass();
			mass+=points[i].getMass();
		}
		return variance/mass;
	}
	
	
	/**
	 * @param points does not need to be sorted or compacted.
	 */
	public static final double expectation(Point[] points, int lastPoint) {
		int i, mass;
		double expectation;
		
		expectation=0.0; mass=0;
		for (i=0; i<=lastPoint; i++) {
			expectation+=points[i].position*points[i].getMass();
			mass+=points[i].getMass();
		}
		return expectation/mass;
	}
	
	
	/**
	 * @param points sorted and compacted;
	 * @param minIntervalLength ignored if non-positive;
	 * @return a value that separates the leftmost peak in $points$, or $defaultThreshold$
	 * if such value cannot be found.
	 */
	public static final double getThreshold(Point[] points, int lastPoint, double defaultThreshold, boolean integerPositions, double minIntervalLength) {
		final int N_CONSECUTIVE_POINTS = 2;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		final double QUANTILE = 0.8;
		int i, nLocalMaximumLeaves;
		double threshold, minLocalMaxDistance;
	
		threshold=defaultThreshold;
		if (lastPoint>0) {
			if (lastPoint+1<=FEW_POINTS) {
				i=getRoughThreshold(points,lastPoint,true,0,true);
				if (i!=-1) threshold=points[i].position;
			}
			else {
				if (minIntervalLength<=0) minIntervalLength=Math.max( distanceQuantile(points,lastPoint,N_CONSECUTIVE_POINTS,integerPositions,QUANTILE),
																	  (points[lastPoint].position-points[0].position)/MIN_INTERVAL_LENGTH_FACTOR
										  						    );
				minLocalMaxDistance=minIntervalLength*2.0;
				if (!areUniformlyDistributed(points,0,lastPoint,false,(points[lastPoint].position-points[0].position)/DEFAULT_NBINS)) {
					DensityEstimationTree.allocateMemory(lastPoint+1);
					nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(points,0,lastPoint,minIntervalLength,minLocalMaxDistance,integerPositions,-1,-1,-1,false,true);
					if (nLocalMaximumLeaves>0) {
						DensityEstimationTree.markRunsOfLocalMaximumLeaves(points);
						threshold=DensityEstimationTree.separateRun(0,points,0,lastPoint,false);
					}
					DensityEstimationTree.deallocateMemory();
				}
			}
		}
		return threshold;
	}
	
	
	/**
	 * Conceptually identical to $getThreshold()$, but just returns the number of peaks
	 * in $points[0..lastPoint]$ (possibly zero).
	 */
	public static final int hasPeaks(Point[] points, int lastPoint, boolean integerPositions, double minIntervalLength) {
		final int N_CONSECUTIVE_POINTS = 2;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		final double QUANTILE = 0.8;
		int out, nLocalMaximumLeaves;
		double minLocalMaxDistance;
		
		if (areUniformlyDistributed(points,0,lastPoint,integerPositions,(points[lastPoint].position-points[0].position)/DEFAULT_NBINS)) return 0;
		if (minIntervalLength<=0) minIntervalLength=Math.max( distanceQuantile(points,lastPoint,N_CONSECUTIVE_POINTS,integerPositions,QUANTILE),
															  (points[lastPoint].position-points[0].position)/MIN_INTERVAL_LENGTH_FACTOR
								  						    );
		minLocalMaxDistance=minIntervalLength*2.0;
		DensityEstimationTree.allocateMemory(lastPoint+1);
		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(points,0,lastPoint,minIntervalLength,minLocalMaxDistance,integerPositions,-1,-1,-1,false,true);
		if (nLocalMaximumLeaves==0) return 0;
		DensityEstimationTree.markRunsOfLocalMaximumLeaves(points);
		out=DensityEstimationTree.getNumberOfRuns();
		DensityEstimationTree.deallocateMemory();
		return out;
	}
	
	
	/**
	 * @param points sorted and compacted;
	 * @return the position of the leftmost local minimum point, starting from $pos$ 
	 * included, or -1 if no such point can be found.
	 */
	public static final int firstLocalMin(int pos, Point[] points, int lastPoint) {
		for (int i=pos==0?1:pos; i<lastPoint; i++) {
			if (points[i].mass<points[i-1].mass && points[i].mass<points[i+1].mass) return i;
		}
		return -1;
	}
	
	
	/**
	 * @param points sorted and compacted;
	 * @return the position of the leftmost local maximum point, starting from $pos$ 
	 * included, or -1 if no such point can be found.
	 */
	public static final int firstLocalMax(int pos, Point[] points, int lastPoint) {
		int i;
		
		if (lastPoint==0) return pos==0?0:-1;
		if (pos==lastPoint) return points[lastPoint].mass>points[lastPoint-1].mass?lastPoint:-1;
		if (pos==0 && points[0].mass>points[1].mass) return 0;
		for (i=pos==0?1:pos; i<lastPoint; i++) {
			if (points[i].mass>points[i-1].mass && points[i].mass>points[i+1].mass) return i;
		}
		if (points[lastPoint].mass>points[lastPoint-1].mass) return lastPoint;
		return -1;
	}
	
	
	/**
	 * @param points sorted and compacted;
	 * @return the number of local-maximum points.
	 */
	public static final int nLocalMax(Point[] points, int lastPoint) {
		int i, out;
		
		if (lastPoint==0) return 1;
		out=0;
		if (points[0].mass>points[1].mass) out++;
		for (i=1; i<lastPoint; i++) {
			if (points[i].mass>points[i-1].mass && points[i].mass>points[i+1].mass) out++;
		}
		if (points[lastPoint].mass>points[lastPoint-1].mass) out++;
		return out;
	}
	
	
	/**
	 * Computes the center of mass of the rightmost interval of $points$ that has length 
	 * at most $lengthThreshold$ and mass at least $massThreshold$.
	 *
	 * @param points assumed to be sorted an compacted;
	 * @param out output array: 0: center of mass; 1: position of the rightmost point of
	 * the interval; $out$ is set to all -1 values if no valid interval can be found.
	 */
	public static final void getRightmostInterval(Point[] points, int lastPoint, int lengthThreshold, int massThreshold, double[] out) {
		int i, j;
		int nElements;
		double sum;
		
		i=lastPoint; sum=0.0; nElements=0;
		while (i>=0) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("getRightmostInterval> i="+i);
			for (j=i; j>=0; j--) {
if (IO.SHOW_STD_ERR_PRIME) System.err.println("getRightmostInterval>    j="+j);
				if (points[i].position-points[j].position>lengthThreshold) break;
				sum+=points[j].mass*points[j].position;
				nElements+=points[j].mass;
			}
if (IO.SHOW_STD_ERR_PRIME) System.err.println("getRightmostInterval> range=["+points[j+1].position+".."+points[i].position+"] mass="+nElements+" theshold="+massThreshold);
			if (nElements>=massThreshold) {
				out[0]=sum/nElements;
				out[1]=points[i].position;
				return;
			}
			i=j; sum=0.0; nElements=0;
		}
		out[0]=-1; 
		out[1]=-1;
		return;
	}
	

}