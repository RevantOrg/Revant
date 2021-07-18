package de.mpi_cbg.revant.util;

import de.mpi_cbg.revant.util.Math;

/**
 * Procedures, on arrays of $Leaf$ objects, that do not require the underlying array of
 * $Point$ objects.
 */
public class Leaves {
	
	public static final double ZERO = 5E-5;  // Arbitrary
	
	/**
	 * Temporary data structures
	 */
	private static Point[] points;
	private static int lastPoint;


	public static final void allocateMemory(int maxLeaves) {
		if (points==null || points.length<maxLeaves) {
			points = new Point[maxLeaves];
			for (int i=0; i<maxLeaves; i++) points[i] = new Point();
		}
	}


	/**
	 * Returns the position of the last leaf of the $run$-th run of local-maximum leaves
	 * in $leaves$.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed, and that $run \in [0..nRuns-1]$.
	 */
	public static final int lastLeafOfRun(Leaf[] leaves, int lastLeaf, int run) {
		int i, currentRun, out;

		currentRun=-1; out=-1;
		for (i=0; i<=lastLeaf; i++) {
			if (!leaves[i].isLocalMaximum) break;
			if (leaves[i].marked) {
				currentRun++;
				if (currentRun==run) {
					out=i;
					break;
				}
			}
		}
		return out;
	}


	/**
	 * Marks as "high" all local-maximum leaves whose value of $getLocalMaxRatio()$ is 
	 * bigger than a threshold. The threshold is estimated by fitting a density estimation
	 * tree on the set of values of all local maxima, and by separating the first local-
	 * maximum from the second local-maximum. The procedure assumes there is no high leaf 
	 * if the ratio between the largest and the smallest value of a leaf is smaller than a
	 * given constant.
	 *
	 * Remark: since the procedure builds a DET, it must be used with care when $leaves$
	 * is itself the array of leaves of a DET. I.e. the array of the DET must be replaced
	 * with a temporary one before calling the procedure, and when the procedure 
     * completes the original array might have to be put back in the DET. See e.g. 
	 * $Intervals.refineBoundaries()$ for an example.
	 *
	 * Remark: the procedure sets $isHigh$ to true for all leaves that it estimates to be 
	 * high, but it does not change the previous $isHigh$ value of other leaves.
	 *
	 * Remark: the DET requires parameter $minIntervalLength$ in input. This is estimated
	 * as a quantile of all distances between the first and the last of a number of
	 * consecutive points. Parameter $minLocalMaxDistance$ is set to twice
	 * $minIntervalLength$.
	 *
	 * @param lastLeafPoint largest point used in a leaf;
	 * @param mode TRUE: $leaves$ is in tree order; FALSE: all local maxima are grouped at
	 * the beginning of $leaves$;
	 * @param inputPoints (ignored if null) if $leaves$ comes from a DET, $inputPoints$ is
	 * used to impose a minimum mass on a local maximum to be marked as high;
	 * @param inputPoints_totalMass total mass in $inputPoints$;
	 * @param maxWidth local maxima wider than this are not marked as high;
	 * @return the total number of high local-maximum leaves, including those that were
	 * marked as high before calling the procedure and that are not marked as high by the 
	 * procedure.
	 */
	public static final int setHighLocalMax(Leaf[] leaves, int lastLeaf, int lastLeafPoint, boolean mode, boolean skipAdjacentMaxima, Point[] inputPoints, int inputPoints_totalMass, boolean detOrRegressionTree, double maxWidth) {
		final int N_CONSECUTIVE_POINTS = 2;
		final double QUANTILE = 0.5;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		final double MIN_HIGH_RATIO = 2.0;
		final double LOW_THRESHOLD = 5.0;  // Arbitrary
		final double HIGH_THRESHOLD = 100.0;  // Arbitrary
		final int MIN_HIGH_MASS_LOWERBOUND = 5;  // Arbitrary
		final double MIN_HIGH_MASS_RATIO = 0.3;  // Arbitrary
		final int MIN_HIGH_MASS = Math.min(MIN_HIGH_MASS_LOWERBOUND,(int)(inputPoints_totalMass*MIN_HIGH_MASS_RATIO));  // Arbitrary
		final double ZERO_RATIO = 2000.0;  // Arbitrary
		int i, tmp;
		int nLocalMaximumLeaves, out;
		double tmpDouble, distance, threshold, defaultThreshold;
		double minIntervalLength, minLocalMaxDistance;
		
		// Computing the number of high leaves before the procedure is called
		out=0;
		for (i=0; i<=lastLeaf; i++) {
			if (leaves[i].isHigh) out++;
		}

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setHighLocalMax> MIN_HIGH_MASS="+MIN_HIGH_MASS);

		// Estimating threshold
		lastPoint=-1;
		for (i=0; i<=lastLeaf; i++) {
			if ( !leaves[i].isLocalMaximum || 
				 ((detOrRegressionTree&&inputPoints!=null)?Points.mass(inputPoints,leaves[i].firstPoint,leaves[i].lastPoint)<MIN_HIGH_MASS:false)
			   ) continue;
			tmpDouble=getLocalMaxRatio(i,leaves,lastLeaf,lastLeafPoint,mode,skipAdjacentMaxima,ZERO_RATIO);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("getLocalMaxRatio("+leaves[i].firstPoint+".."+leaves[i].lastPoint+")="+tmpDouble+", mode="+mode);
			if (tmpDouble!=-1.0 && tmpDouble!=ZERO_RATIO) {
				lastPoint++;
				points[lastPoint].position=tmpDouble;
				points[lastPoint].mass=1;
			}
		}
		if (lastPoint==-1) {
			// No local maximum has a valid ratio. We assume this is because $leaves$ 
			// contains just local maxima and zero-density leaves, so we mark all local
			// maxima as high.
			out=0;
			for (i=0; i<=lastLeaf; i++) {
				if ( leaves[i].isLocalMaximum && 
				     ((detOrRegressionTree&&inputPoints!=null)?Points.mass(inputPoints,leaves[i].firstPoint,leaves[i].lastPoint)>=MIN_HIGH_MASS:true) && 
				     (detOrRegressionTree?(inputPoints!=null?inputPoints[leaves[i].lastPoint].position-inputPoints[leaves[i].firstPoint].position<=maxWidth:true):leaves[i].lastPoint-leaves[i].firstPoint+1<=maxWidth)
				   ) {
					leaves[i].isHigh=true;
					out++;
				}
			}
			return out;
		}
		if (lastPoint>0) lastPoint=Points.sortAndCompact(points,lastPoint);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setHighLocalMax> points:");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastPoint; i++) IO.printErr(points[i].position+","+points[i].getMass()); }


		defaultThreshold=points[0].position>LOW_THRESHOLD?points[0].position:LOW_THRESHOLD;
		if (lastPoint==0) threshold=points[lastPoint].position;  // A single local maximum is always considered high
		else if (points[lastPoint].position<points[0].position*MIN_HIGH_RATIO) threshold=defaultThreshold;
		else if (lastPoint+1<=Points.FEW_POINTS) {
			// It is typically hard to set meaningful thresholds on the values of leaves,
			// thus we don't try to decide whether few points are uniformly distributed,
			// and we always split.
			tmp=Points.getRoughThreshold(points,lastPoint,true,0,false);
			threshold=tmp==-1?defaultThreshold:points[tmp].position;
		}
		else {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setHighLocalMax> X2");
			minIntervalLength=Math.max( Points.distanceQuantile(points,lastPoint,N_CONSECUTIVE_POINTS,false,QUANTILE),
										(points[lastPoint].position-points[0].position)/MIN_INTERVAL_LENGTH_FACTOR
			 						  );
			minLocalMaxDistance=minIntervalLength*2.0;
			if (Points.areUniformlyDistributed(points,0,lastPoint,false,(points[lastPoint].position-points[0].position)/Points.DEFAULT_NBINS)) threshold=defaultThreshold;
			else {
				nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(points,0,lastPoint,minIntervalLength,minLocalMaxDistance,false,-1,-1,-1,false,true);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setHighLocalMax> nLocalMaximumLeaves="+nLocalMaximumLeaves);
				if (nLocalMaximumLeaves==0) threshold=defaultThreshold;
				else {
					DensityEstimationTree.markRunsOfLocalMaximumLeaves(points);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("DET OF LOCAL MAXIMA:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr("["+points[DensityEstimationTree.leaves[x].firstPoint].position+".."+points[DensityEstimationTree.leaves[x].lastPoint].position+"] isLocalMaximum="+DensityEstimationTree.leaves[x].isLocalMaximum+" isLocalMinimum="+DensityEstimationTree.leaves[x].isLocalMinimum); }


					threshold=DensityEstimationTree.separateRun(0,points,0,lastPoint,true);
				}
			}
		}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setHighLocalMax> threshold before: "+threshold);
		if (threshold>LOW_THRESHOLD) {
			if (threshold<HIGH_THRESHOLD) threshold=LOW_THRESHOLD;
			else threshold=HIGH_THRESHOLD;
		}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("setHighLocalMax> threshold after: "+threshold);

		// Marking high local maxima (local maxima with no valid ratio are marked as high
		// as well).
		for (i=0; i<=lastLeaf; i++) {
			if ( leaves[i].isLocalMaximum && 
			     ((detOrRegressionTree&&inputPoints!=null)?Points.mass(inputPoints,leaves[i].firstPoint,leaves[i].lastPoint)>=MIN_HIGH_MASS:true) && 
				 (detOrRegressionTree?(inputPoints!=null?inputPoints[leaves[i].lastPoint].position-inputPoints[leaves[i].firstPoint].position<=maxWidth:true):leaves[i].lastPoint-leaves[i].firstPoint+1<=maxWidth)
			   ) {
				tmpDouble=getLocalMaxRatio(i,leaves,lastLeaf,lastLeafPoint,mode,skipAdjacentMaxima,ZERO_RATIO);
				if (tmpDouble==-1 || tmpDouble>=threshold) {
					leaves[i].isHigh=true;
				}
			}
		}
		out=0;
		for (i=0; i<=lastLeaf; i++) {
			if (leaves[i].isHigh) out++;
		}
		return out;
	}


	/**
	 * Returns the ID of the first local-minimum leaf. The procedure assumes that
	 * $markRunsOfLocalMaximumLeaves$ has already been executed.
	 *
	 * @return $lastLeaf+1$ if no local-minimum leaf exists.
	 */
	public static final int firstLocalMinimum(Leaf[] leaves, int lastLeaf) {
		int i = 0;
		while (i<=lastLeaf && !leaves[i].isLocalMinimum) i++;
		return i;
	}


	/**
	 * Returns the ID of the first leaf that is neither a local maximum nor a local
	 * minimum. The procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 *
	 * @return $lastLeaf+1$ if no normal leaf exists.
	 */
	public static final int firstNormal(Leaf[] leaves, int lastLeaf) {
		int i = 0;
		while (i<=lastLeaf && (leaves[i].isLocalMaximum||leaves[i].isLocalMinimum)) i++;
		return i;
	}


	/**
	 * Returns the ID of the leaf immediately to the right (if $direction=true$) or to the
	 * left (if $direction=false$) of $leaf$, which is not a local maximum, and which
	 * fully occurs before (respectively, after) $boundary$. If $mode=true$, priority is 
	 * given to local-minimum leaves.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 *
	 * @param firstLocalMinimum ID of the first local-minimum leaf;
	 * @param firstNormal ID of the first leaf that is neither a local-maximum nor a
	 * local-minimum;
	 * @return -1 if no non-maximum leaf can be found in the given direction.
	 */
	public static final int nextNonMaximum(Leaf[] leaves, int lastLeaf, int leaf, boolean direction, int boundary, int firstLocalMinimum, int firstNormal, boolean mode) {
		int i;
		int nextLocalMinimum, nextNormal;

		nextLocalMinimum=-1;
		nextNormal=-1;
		if (direction) {
			for (i=firstLocalMinimum; i<firstNormal; i++) {
				if (leaves[i].firstPoint>leaves[leaf].lastPoint && leaves[i].lastPoint<boundary) {
					nextLocalMinimum=i;
					break;
				}
			}
			if (mode && nextLocalMinimum>=0) return nextLocalMinimum;
			for (i=firstNormal; i<=lastLeaf; i++) {
				if (leaves[i].firstPoint>leaves[leaf].lastPoint && leaves[i].lastPoint<boundary) {
					nextNormal=i;
					break;
				}
			}
			if (nextNormal>=0) {
				if (nextLocalMinimum>=0) return leaves[nextLocalMinimum].firstPoint<leaves[nextNormal].firstPoint?nextLocalMinimum:nextNormal;
				else return nextNormal;
			}
		}
		else {
			for (i=firstLocalMinimum; i<firstNormal && leaves[i].lastPoint<leaves[leaf].firstPoint; i++) nextLocalMinimum=i;
			if (mode && nextLocalMinimum>=0 && leaves[nextLocalMinimum].firstPoint>boundary) return nextLocalMinimum;
			for (i=firstNormal; i<=lastLeaf && leaves[i].lastPoint<leaves[leaf].firstPoint; i++) nextNormal=i;
			if (nextNormal>=0 && leaves[nextNormal].firstPoint>boundary) {
				if (nextLocalMinimum>=0) return leaves[nextLocalMinimum].lastPoint>leaves[nextNormal].lastPoint?nextLocalMinimum:nextNormal;
				else return nextNormal;
			}
		}
		return -1;
	}


	/**
	 * Returns the ID of the leftmost local-minimum leaf whose interval of points is
	 * contained in $[firstPoint..lastPoint]$. The procedure assumes that
	 * $markRunsOfLocalMaximumLeaves$ has already been executed.
	 *
	 * @return -1 if no such local-minimum leaf can be found.
	 */
	public static final int includedLocalMinimum(Leaf[] leaves, int lastLeaf, int firstPoint, int lastPoint, int firstLocalMinimum, int firstNormal) {
		for (int i=firstLocalMinimum; i<firstNormal; i++) {
			if (leaves[i].firstPoint<firstPoint || leaves[i].lastPoint>lastPoint) continue;
			return i;
		}
		return -1;
	}
	
	
	/**
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	public static final int getNumberOHighRuns(Leaf[] leaves, int lastLocalMaximum) {
		boolean isHigh;
		int i, first, out;
		
		out=0; first=0;
		while (first<=lastLocalMaximum) {
			i=first; isHigh=false;
			while (!leaves[i].marked) {
				if (leaves[i].isHigh) isHigh=true;
				i++;
			}
			if (leaves[i].isHigh) isHigh=true;
			if (isHigh) out++;
			first=i+1;
		}
		return out;
	}
	
	
	/**
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 *
	 * @return the last point of the $run$-th high run, or -1 if no such high run can be
	 * found.
	 */
	public static final int getLastPointOfHighRun(Leaf[] leaves, int lastLocalMaximum, int run) {
		boolean isHigh;
		int i, first, currentHighRun;
		
		currentHighRun=-1; first=0;
		while (first<=lastLocalMaximum) {
			i=first; isHigh=false;
			while (!leaves[i].marked) {
				if (leaves[i].isHigh) isHigh=true;
				i++;
			}
			if (leaves[i].isHigh) isHigh=true;
			if (isHigh) {
				currentHighRun++;
				if (currentHighRun==run) return leaves[i].lastPoint;
			}
			first=i+1;
		}
		return -1;
	}
	
	
	/**
	 * @param mode TRUE: $leaves$ is in tree order; FALSE: all local maxima are grouped at
	 * the beginning of $leaves$ and sorted in tree order;
	 * @param skipAdjacentMaxima starts the search after a sequence of local-maximum 
	 * leaves that are consecutive to each other and to $leaf$ in the specified direction;
	 * @param zeroRatio the procedure assigns this default ratio when the denominator of 
	 * the ratio would otherwise be zero;
	 * @return the ratio between the value of $leaf$ and the average value of its adjacent
	 * leaves in tree order; -1.0 if neither a valid left leaf nor a valid right leaf can 
	 * be found.
	 */
	public static final double getLocalMaxRatio(int leaf, Leaf[] leaves, int lastLeaf, int lastPoint, boolean mode, boolean skipAdjacentMaxima, double zeroRatio) {
		double leftValue, rightValue, denominator;
		
		if (mode) {
			rightValue=nextLeafValue_treeOrder(leaf,leaves,lastLeaf,true,skipAdjacentMaxima);
			leftValue=nextLeafValue_treeOrder(leaf,leaves,lastLeaf,false,skipAdjacentMaxima);
		}
		else {
			rightValue=nextLeafValue_localMaxOrder(leaf,leaves,lastLeaf,lastPoint,true,skipAdjacentMaxima);
			leftValue=nextLeafValue_localMaxOrder(leaf,leaves,lastLeaf,lastPoint,false,skipAdjacentMaxima);
		}
		if (leftValue==-1.0) denominator=rightValue;
		else if (rightValue==-1.0) denominator=leftValue;
		else denominator=Math.min(leftValue,rightValue);
		if (denominator==-1.0) return -1.0;
		else if (denominator==0.0) return zeroRatio;
		else return Math.min(leaves[leaf].value/denominator,zeroRatio);
	}
	
	
	/**
	 * @param leaves assumed to be sorted in tree order;
	 * @param skipAdjacentMaxima starts the search after a sequence of local-maximum 
	 * leaves that are consecutive to each other and to $leaf$ in the specified direction;
	 * @return the value of the leftmost leaf to the right of $leaf$ (if $rightLeaf=
	 * true$), or of the rightmost leaf to the left of $leaf$, that is a local minimum; if
	 * such a leaf cannot be found, the procedure returns the value of a closest leaf that
	 * is not a local maximum; otherwise, the procedure returns -1.0. The procedure 
	 * assumes $leaves[leaf]$ to be a local maximum.
	 */
	private static final double nextLeafValue_treeOrder(int leaf, Leaf[] leaves, int lastLeaf, boolean rightLeaf, boolean skipAdjacentMaxima) {
		int i, j;
		double normalValue;
	
		normalValue=-1.0;
		if (rightLeaf) {
			if (leaf==lastLeaf) return -1.0;
			j=leaf+1;
			if (skipAdjacentMaxima) {
				while (j<=lastLeaf && leaves[j].isLocalMaximum && leaves[j].firstPoint==leaves[j-1].lastPoint+1) j++;
			}
			for (i=j; i<=lastLeaf; i++) {
				if (leaves[i].isLocalMaximum) break;
				if (leaves[i].isLocalMinimum) {
					if (leaves[i].value>ZERO) return leaves[i].value;
					if (i>0 && !leaves[i-1].isLocalMaximum && leaves[i-1].value>ZERO) return leaves[i-1].value;
					break;
				}
				if (normalValue==-1.0 && leaves[i].value>ZERO) normalValue=leaves[i].value;
			}
		}
		else {
			if (leaf==0) return -1.0;
			j=leaf-1;
			if (skipAdjacentMaxima) {
				while (j>=0 && leaves[j].isLocalMaximum && leaves[j].lastPoint==leaves[j+1].firstPoint-1) j--;
			}	
			for (i=j; i>=0; i--) {
				if (leaves[i].isLocalMaximum) break;
				if (leaves[i].isLocalMinimum) {
					if (leaves[i].value>ZERO) return leaves[i].value;
					if (i<lastLeaf && !leaves[i+1].isLocalMaximum && leaves[i+1].value>ZERO) return leaves[i+1].value;
					break;
				}
				if (normalValue==-1.0 && leaves[i].value>ZERO) normalValue=leaves[i].value;
			}
		}
		return normalValue;
	}
	
	
	/**
	 * Like $nextLeafValue_treeOrder()$, but assumes that all local maxima are grouped at 
	 * the beginning of $leaves$ and sorted by tree order, and that $leaf$ is a local 
	 * maximum.
	 */
	private static final double nextLeafValue_localMaxOrder(int leaf, Leaf[] leaves, int lastLeaf, int lastPoint, boolean rightLeaf, boolean skipAdjacentMaxima) {
		int i, j, p;
		int point, nextMaximumPoint, localMinFirstPoint, localMinLastPoint, normalPoint;
		int selectedLeaf;
		double localMinValue, normalValue;
		
		localMinValue=-1.0; normalValue=-1.0; 
		if (rightLeaf) {
			if (leaf==lastLeaf) return -1.0;
			point=leaves[leaf].lastPoint+1;
			if (point>lastPoint) return -1.0;
			j=leaf+1;
			if (skipAdjacentMaxima) {
				while (j<=lastLeaf && leaves[j].isLocalMaximum && leaves[j].firstPoint==leaves[j-1].lastPoint+1) j++;
			}
			nextMaximumPoint=leaves[j].isLocalMaximum?leaves[j].firstPoint:Math.POSITIVE_INFINITY;
			localMinFirstPoint=Math.POSITIVE_INFINITY; localMinLastPoint=-1; normalPoint=Math.POSITIVE_INFINITY;
			for (i=j; i<=lastLeaf; i++) {
				if (leaves[i].isLocalMaximum) continue;
				p=leaves[i].firstPoint;
				if (p>=nextMaximumPoint) continue;
				if (p>=point) {
					if (leaves[i].isLocalMinimum && p<localMinFirstPoint) {
						localMinFirstPoint=p;
						localMinLastPoint=leaves[i].lastPoint;
						localMinValue=leaves[i].value;
					}
					if (!leaves[i].isLocalMinimum && p<normalPoint && leaves[i].value>ZERO) {
						normalPoint=p;
						normalValue=leaves[i].value;
					}
				}
			}
		}
		else {
			point=leaves[leaf].firstPoint-1;
			if (point<0) return -1.0;
			j=leaf-1;
			if (skipAdjacentMaxima) {
				while (j>=0 && leaves[j].lastPoint==leaves[j+1].firstPoint-1) j--;
			}
			nextMaximumPoint=j==-1?Math.NEGATIVE_INFINITY:leaves[j].lastPoint;
			localMinFirstPoint=-1; localMinLastPoint=-1; normalPoint=-1;
			for (i=leaf+1; i<=lastLeaf; i++) {
				if (leaves[i].isLocalMaximum) continue;
				p=leaves[i].lastPoint;
				if (p<=nextMaximumPoint) continue;
				if (p<=point) {
					if (leaves[i].isLocalMinimum && p>localMinLastPoint) {
						localMinLastPoint=p;
						localMinFirstPoint=leaves[i].firstPoint;
						localMinValue=leaves[i].value;
					}
					if (!leaves[i].isLocalMinimum && p>normalPoint && leaves[i].value>ZERO) {
						normalPoint=p;
						normalValue=leaves[i].value;
					}
				}
			}
		}
		if (localMinValue!=-1) {
			if (localMinValue>ZERO) return localMinValue;
			selectedLeaf=-1;
			if (rightLeaf) {
				for (i=lastLeaf; i>=0; i--) {
					if (leaves[i].isLocalMaximum) break;
					if (leaves[i].isLocalMinimum || leaves[i].value<=ZERO) continue;
					if (leaves[i].lastPoint<=localMinFirstPoint-1 && (selectedLeaf==-1 || leaves[i].lastPoint>leaves[selectedLeaf].lastPoint)) selectedLeaf=i;
				}
				if (selectedLeaf!=-1) return leaves[selectedLeaf].value; 
			}
			else {
				for (i=lastLeaf; i>=0; i--) {
					if (leaves[i].isLocalMaximum) break;
					if (leaves[i].isLocalMinimum || leaves[i].value<=ZERO) continue;
					if (leaves[i].firstPoint>=localMinLastPoint+1 && (selectedLeaf==-1 || leaves[i].firstPoint<leaves[selectedLeaf].firstPoint)) selectedLeaf=i;
				}
				if (selectedLeaf!=-1) return leaves[selectedLeaf].value;
			}
		}
		return normalValue;
	}

	
	/**
	 * @param points if null, the procedure assumes that $leaves$ comes from a histogram;
	 * @param integerPositions used only if $points!=null$;
	 * @return the surface of all local-max leaves.
	 */
	public static final double getLocalMaxSurface(Leaf[] leaves, int lastLeaf, Point[] points, boolean integerPositions) {
		int i;
		double surface;
		
		surface=0.0;
		if (points==null) {
			for (i=0; i<=lastLeaf; i++) {
				if (leaves[i].isLocalMaximum) surface+=leaves[i].lastPoint-leaves[i].firstPoint+1;
			}
		}
		else {
			for (i=0; i<=lastLeaf; i++) {
				if (leaves[i].isLocalMaximum) surface+=points[leaves[i].lastPoint].position-points[leaves[i].firstPoint].position+(integerPositions?1:0);
			}
		}
		return surface;
	}
	
	
	/**
	 * @param points if null, the procedure assumes that $leaves$ comes from a histogram;
	 * @param integerPositions used only if $points!=null$;
	 */
	public static final double getNonzeroSurface(Leaf[] leaves, int lastLeaf, Point[] points, boolean integerPositions) {
		int i;
		double surface;
		
		surface=0;
		if (points==null) {
			for (i=0; i<=lastLeaf; i++) {
				if (leaves[i].value>ZERO) surface+=leaves[i].lastPoint-leaves[i].firstPoint+1;
			}
		}
		else {
			for (i=0; i<=lastLeaf; i++) {
				if (leaves[i].value>ZERO) surface+=points[leaves[i].lastPoint].position-points[leaves[i].firstPoint].position+(integerPositions?1:0);
			}
		}
		return surface;
	}
	
	
	/**
	 * @param minDensity,maxWidth only local-maxima with density and width compatible with
	 * such thresholds are considered;
	 * @return TRUE iff there are at least $minLocalMaxima$ local maxima in $leaves$, and 
	 * if at least $peaksRatio$ of them are at distance $k*period$ from another peak (up 
	 * to $distanceThreshold$ tolerance).
	 */
	public static final boolean peaksSpacedByPeriod(Leaf[] leaves, int lastLeaf, double[] histogram, int period, int minLocalMaxima, int distanceThreshold, double peaksRatio, double minDensity, int maxWidth) {
		final int MAX_MULTIPLE = 5;  // Arbitrary
		boolean found;
		int i, j, p;
		int center, numerator, denominator;
		
		for (p=1; p<=MAX_MULTIPLE; p++) {
			numerator=0; denominator=0;
			for (i=0; i<=lastLeaf; i++) {
				if (!leaves[i].isLocalMaximum || leaves[i].value<minDensity || leaves[i].lastPoint-leaves[i].firstPoint+1>maxWidth) continue;
				denominator++;
				center=(int)Histograms.getCenterOfMass(histogram,leaves[i].firstPoint,leaves[i].lastPoint);
				found=false;
				for (j=i-1; j>=0; j--) {
					if (!leaves[j].isLocalMaximum || leaves[j].value<minDensity || leaves[j].lastPoint-leaves[j].firstPoint+1>maxWidth) continue;
					if (Math.abs(center-(int)Histograms.getCenterOfMass(histogram,leaves[j].firstPoint,leaves[j].lastPoint),p*period)<=distanceThreshold) {
						found=true;
						break;
					}
				}
				if (!found) {
					for (j=i+1; j<=lastLeaf; j++) {
						if (!leaves[j].isLocalMaximum || leaves[j].value<minDensity || leaves[j].lastPoint-leaves[j].firstPoint+1>maxWidth) continue;
						if (Math.abs((int)Histograms.getCenterOfMass(histogram,leaves[j].firstPoint,leaves[j].lastPoint)-center,p*period)<=distanceThreshold) {
							found=true;
							break;
						}
					}
				}
				if (found) numerator++;
			}
			if (denominator>=minLocalMaxima && numerator>=denominator*peaksRatio) return true;
		}
		return false;
	}


}