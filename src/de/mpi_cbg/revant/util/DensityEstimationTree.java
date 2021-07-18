package de.mpi_cbg.revant.util;

import java.util.Arrays;
import de.mpi_cbg.revant.util.Math;

/**
 * A version of [1] modified for working in one dimension and without training.
 *
 * [1] Ram, Parikshit, and Alexander G. Gray. "Density estimation trees." Proceedings of 
 * the 17th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining.
 * ACM, 2011.
 */
public class DensityEstimationTree {
	
	/**
	 * Parameters of the pipeline
	 */
	private static final int MIN_LEAF_MASS = 3;  // Leaves with mass smaller than this are discarded
	public static int minLocalMaxMass = 2;  // Minimum mass of a local-maximum leaf for it not to be discarded
	public static double minIntervalLength;
	public static double minLocalMaxDistance;
	public static boolean integerPositions;
	public static double splitErrorRatio = 2.0;  // An interval can be split into intervals of length smaller than $minIntervalLength$, iff the new error (negative) is at most this multiple of the old error (negative).

	/**
	 * Output data structures
	 */
	public static Leaf[] leaves;  // The sorted list of leaves
	public static int lastLeaf;  // Position of the last leaf in $leaves$
	public static int nRuns;  // Number of runs of local-maximum leaves, computed by $markRunsOfLocalMaximumLeaves$.
	public static int lastLocalMaximum;  // Position of the last local-maximum leaf in $leaves$ after $markRunsOfLocalMaximumLeaves$
	public static int lastLocalMinimum;  // Position of the last local-minimum leaf in $leaves$ after $markRunsOfLocalMaximumLeaves$

	/**
	 * Temporary data structures
	 */
	private static double[] stack;  // Temporary space used to traverse the tree
	private static double[] tmp;
	public static int[] highRunTmp;


	public static final void allocateMemory(int maxLeaves) {
		if (leaves==null || leaves.length<maxLeaves) {
			leaves = new Leaf[maxLeaves];
			for (int i=0; i<maxLeaves; i++) leaves[i] = new Leaf();
		}
		if (stack==null || stack.length<maxLeaves*3) stack = new double[maxLeaves*3];
		if (tmp==null) tmp = new double[6];
		if (highRunTmp==null) highRunTmp = new int[2];
	}
	
	
	public static final void deallocateMemory() {
		for (int i=0; i<leaves.length; i++) leaves[i]=null;
		leaves=null;
		stack=null;
		tmp=null;
		highRunTmp=null;
	}


	/**
	 * Builds the intervals of a density estimation tree on the points in
	 * $points[firstPoint..lastPoint]$. Leaves are stored in array $leaves[0..lastLeaf]$.
	 * $points$ is assumed to be sorted by $position$.
	 *
	 * @param minIntervalLength intervals are not split if they are this length or
	 * smaller;
	 * @param minLocalMaxDistance minimum distance between two local maxima for them to be
	 * considered distinct clusters;
	 * @param integerPositions TRUE iff points have integer positions;
	 * @param largestMass -1 if procedure $addPointsWithLargestMass$ should not be called;
	 * @param minHoleLength -1 if procedure $removeHoles$ should not be called;
	 * @param shortEnds allows intervals shorter than $minIntervalLength$ at the beginning
	 * or end of $points[firstPoint..lastPoint]$;
	 * @return the number of local-maximum leaves.
	 */
	public static final int buildDensityEstimationTree(Point[] points, int firstPoint, int lastPoint, double minIntervalLength, double minLocalMaxDistance, boolean integerPositions, int largestMass, double minHoleLength, double smallLeafDistance, boolean shortEnds, boolean removeSmallLeaves) {
		boolean split;
		int i, top, topPrime;
		int intervalFirst, intervalLast;
		double intervalLength, splitPoint, intervalError, totalMass, intervalMass, intervalDensity;
		DensityEstimationTree.minIntervalLength=minIntervalLength;
		DensityEstimationTree.minLocalMaxDistance=minLocalMaxDistance;
		DensityEstimationTree.integerPositions=integerPositions;

		// Pushing the root onto $stack$
		totalMass=0.0;  // Not normalized
		for (i=firstPoint; i<=lastPoint; i++) totalMass+=points[i].getMass();
		intervalMass=totalMass;
		intervalLength=points[lastPoint].position-points[firstPoint].position+(integerPositions?1:0);
		intervalDensity=intervalMass/intervalLength;  // Not normalized
		intervalError=-1.0/intervalLength;  // Normalized
		stack[0]=firstPoint; stack[1]=lastPoint; stack[2]=intervalError;
		top=0; lastLeaf=-1;

		// Building the tree in preorder
		while (top>=0) {
			intervalFirst=(int)stack[top];
			intervalLast=(int)stack[top+1];
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("buildDensityEstimationTree> considering interval ["+points[intervalFirst].position+".."+points[intervalLast].position+"]  intervalFirst="+intervalFirst+" intervalLast="+intervalLast);
			intervalError=stack[top+2];  // Normalized
			intervalLength=points[intervalLast].position-points[intervalFirst].position+(integerPositions?1:0);
			intervalMass=Math.sqrt(-intervalError*intervalLength);  // Normalized
			intervalDensity=intervalMass/intervalLength;  // Normalized
			top-=3;
			split=true;
			if (intervalFirst==intervalLast || intervalLength<minIntervalLength) split=false;
			else {
				split(points,intervalFirst,intervalLast,totalMass);
				if ( tmp[0]==-1 ||
				     tmp[1]+tmp[2]>=intervalError ||
				     ( tmp[5]==0 && 
					   !( tmp[1]+tmp[2]<=intervalError*splitErrorRatio || (shortEnds&&(intervalFirst==firstPoint||intervalLast==lastPoint)) )
					 )
				   ) split=false;
			}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("buildDensityEstimationTree> split="+split+" tmp[1]+tmp[2]="+(tmp[1]+tmp[2])+" intervalError="+intervalError+" intervalMass="+intervalMass+" intervalLength="+intervalLength+" minIntervalLength="+minIntervalLength);
			if (!split) {
				lastLeaf++;
				leaves[lastLeaf].firstPoint=intervalFirst;
				leaves[lastLeaf].lastPoint=intervalLast;
				leaves[lastLeaf].value=intervalDensity*totalMass;  // Not normalized
			}
			else {
				splitPoint=tmp[0];
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("buildDensityEstimationTree> splitting at point "+points[(int)splitPoint].position);
				if (tmp[5]==1) {
					// Recurring on the right child
					topPrime=top+3;
					stack[topPrime]=splitPoint+1;
					stack[topPrime+1]=intervalLast;
					stack[topPrime+2]=tmp[2];
					// Recurring on the left child
					topPrime+=3;
					stack[topPrime]=intervalFirst;
					stack[topPrime+1]=splitPoint;
					stack[topPrime+2]=tmp[1];
					top=topPrime;
				}
				else {
					// Saving the left child
					lastLeaf++;
					leaves[lastLeaf].firstPoint=intervalFirst;
					leaves[lastLeaf].lastPoint=(int)splitPoint;
					leaves[lastLeaf].value=tmp[3];
					// Saving the right child
					lastLeaf++;
					leaves[lastLeaf].firstPoint=(int)splitPoint+1;
					leaves[lastLeaf].lastPoint=intervalLast;
					leaves[lastLeaf].value=tmp[4];
				}
			}
		}

		// Removing holes
		if (minHoleLength!=-1) removeHoles(minHoleLength,points,smallLeafDistance);
		setLocalMaxMin(points);

		// Adding new leaves for points with high mass
		for (i=0; i<=lastLeaf; i++) leaves[i].isPointWithLargestMass=false;
		if (integerPositions && largestMass>0) addPointsWithLargestMass(points,firstPoint,lastPoint,largestMass);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("largestMass="+largestMass);
		for (i=0; i<=lastLeaf; i++) {
			if (leaves[i].isPointWithLargestMass) leaves[i].isHigh=true;
			else leaves[i].isHigh=false;
		}		
		setLocalMaxMin(points);
		if (removeSmallLeaves) discardSmallLeaves(points);
		computeDistances(points);
		return setLocalMaxMin(points);
	}


	/**
	 * Tries to split interval $points[firstPoint..lastPoint]$ so that both of the
	 * resulting intervals have length at least $minIntervalLength$. If this is not
	 * possible, the procedure splits the interval at the global point that minimizes the
	 * error.
	 *
	 * Remark: $points$ is assumed to be sorted by $position$.
	 *
	 * @param totalMass total mass of the dataset (not normalized);
	 * @return $tmp[0]$: the point $i$ that minimizes the sum of the errors of the two
	 * partitions $points[firstPoint..i],points[i+1..lastPoint]$. The intervals induced
	 * by splitting at this point are both of length at least $minIntervalLength$ if
	 * $tmp[5]=1$, they are not if $tmp[5]=0$. $tmp[0]=-1$ if no splitting point can be
	 * found.
	 * $tmp[1]$: the error of the left partition (normalized);
	 * $tmp[2]$: the error of the right partition (normalized);
	 * $tmp[3]$: the density of the left partition (not normalized);
	 * $tmp[4]$: the density of the right partition (not normalized).
	 * The return value of the procedure is the mass of $points[firstPoint..lastPoint]$
	 * (not normalized).
	 */
	private static final double split(Point[] points, int firstPoint, int lastPoint, double totalMass) {
		final double LOWER_OVERFLOW_THRESHOLD = 1e-10;
		final double UPPER_OVERFLOW_THRESHOLD = 1e10;
		int i, j, from, to;
		int minPointAllowed, minPointGlobal;
		double lengthLeft, lengthRight, mass, massLeft, massRight, intervalMass, errorLeft, errorRight;
		double minErrorAllowed, minErrorGlobal;
		double minErrorAllowedLeft, minErrorGlobalLeft, minErrorAllowedRight, minErrorGlobalRight;
		double minDensityAllowedLeft, minDensityGlobalLeft, minDensityAllowedRight, minDensityGlobalRight;

		from=firstPoint+1;
		to=lastPoint-2;
		if (to<from) {
			tmp[0]=-1;
			tmp[1]=1; tmp[2]=1;  // Error is always negative
			tmp[3]=-1; tmp[4]=-1;
			tmp[5]=-1;
		}
		massLeft=points[firstPoint].getMass();
		massRight=0.0;
		for (i=from; i<=lastPoint; i++) massRight+=points[i].getMass();
		intervalMass=massLeft+massRight;
		minPointAllowed=-1; minPointGlobal=-1;
		minErrorAllowed=1; minErrorGlobal=1;  // Error is always negative
		minErrorAllowedLeft=1; minErrorAllowedRight=1;  // Error is always negative
		minErrorGlobalLeft=1; minErrorGlobalRight=1;  // Error is always negative
		minDensityAllowedLeft=-1; minDensityAllowedRight=-1;
		minDensityGlobalLeft=-1; minDensityGlobalRight=-1;
		for (i=from; i<=to; i++) {
			mass=points[i].getMass();
			massLeft+=mass; 
			if (massLeft<0 || massLeft<=LOWER_OVERFLOW_THRESHOLD || massLeft>=UPPER_OVERFLOW_THRESHOLD) {
				break;  // Overflow
			}
			massRight-=mass;
			if (massRight<0 || massRight<=LOWER_OVERFLOW_THRESHOLD || massRight>=UPPER_OVERFLOW_THRESHOLD) {
				break;  // Overflow
			}
			lengthLeft=points[i].position-points[firstPoint].position+(integerPositions?1:0);
			lengthRight=points[lastPoint].position-points[i+1].position+(integerPositions?1:0);
			errorLeft=-(massLeft/totalMass)*(massLeft/totalMass)/lengthLeft;
			errorRight=-(massRight/totalMass)*(massRight/totalMass)/lengthRight;
			if (lengthLeft>=minIntervalLength && lengthRight>=minIntervalLength && errorLeft+errorRight<minErrorAllowed) {
				minErrorAllowed=errorLeft+errorRight;
				minErrorAllowedLeft=errorLeft;
				minErrorAllowedRight=errorRight;
				minDensityAllowedLeft=massLeft/lengthLeft;
				minDensityAllowedRight=massRight/lengthRight;
				minPointAllowed=i;
			}
			if (errorLeft+errorRight<minErrorGlobal) {
				minErrorGlobal=errorLeft+errorRight;
				minErrorGlobalLeft=errorLeft;
				minErrorGlobalRight=errorRight;
				minDensityGlobalLeft=massLeft/lengthLeft;
				minDensityGlobalRight=massRight/lengthRight;
				minPointGlobal=i;
			}
		}

		// Outputting
		if (minPointAllowed!=-1) {
			tmp[0]=minPointAllowed;
			tmp[1]=minErrorAllowedLeft;
			tmp[2]=minErrorAllowedRight;
			tmp[3]=minDensityAllowedLeft;
			tmp[4]=minDensityAllowedRight;
			tmp[5]=1;
		}
		else {
			tmp[0]=minPointGlobal;
			tmp[1]=minErrorGlobalLeft;
			tmp[2]=minErrorGlobalRight;
			tmp[3]=minDensityGlobalLeft;
			tmp[4]=minDensityGlobalRight;
			tmp[5]=0;
		}
		return intervalMass;
	}


	/**
	 * Sets fields $isLocalMaximum$ and $isLocalMinimum$ of each leaf in $leaves$.
	 *
	 * Let a run be a sequence of consecutive leaves with identical value, and such that
	 * the distance between two consecutive leaves is smaller than $minLocalMaxDistance$.
	 * We say that a run is left (respectively, right) maximum if the value of the leaves
	 * in the run is bigger than the value of the leaves in the previous (respectively,
	 * next) adjacent run, or if such adjacent run is at least $minLocalMaxDistance$
	 * positions away. We say that a run is left (respectively, right) minimum if the
	 * value of the leaves in the run is smaller than the value of the leaves in the
	 * previous (respectively, next) adjacent run, and if such adjacent run is less than
	 * $minLocalMaxDistance$ positions away. A run is a local maximum (minimum) if it is
	 * both left and right maximum (minimum). The procedure marks as local maximum just
	 * the leaf in the run that contains the center of mass of all points in the run. The
	 * procedure marks as local minimum all leaves in a local-minimum run.
	 *
	 * @return the number of local maxima found.
	 */
	private static final int setLocalMaxMin(Point[] points) {
		final double THRESHOLD = 0.1;  // Arbitrary
		final double EPSILON_THRESHOLD = 0.0001;  // Arbitrary
		boolean inLeaf, found;
		int i, j, k;
		int mass, out;
		double centerOfMass, epsilon;

		for (i=0; i<=lastLeaf; i++) {
			leaves[i].isLocalMaximum=false;
			leaves[i].isLocalMinimum=false;
		}
		if (lastLeaf==0) {
			if (points[leaves[0].lastPoint].position-points[leaves[0].firstPoint].position>=minLocalMaxDistance) return 0;
			leaves[0].isLocalMaximum=true;
			return 1;
		}
		out=0;
		i=0;
		while (i<=lastLeaf) {
			centerOfMass=0.0;
			mass=0;
			for (k=leaves[i].firstPoint; k<=leaves[i].lastPoint; k++) {
				mass+=points[k].getMass();
				centerOfMass+=points[k].position*points[k].getMass();
			}
			j=i+1;
			while ( j<=lastLeaf &&
			        points[leaves[j].firstPoint].position-points[leaves[j-1].lastPoint].position<minLocalMaxDistance &&
					Math.abs(leaves[j].value-leaves[i].value)<=Math.max(leaves[j].value,leaves[i].value)*THRESHOLD
			      ) {
				for (k=leaves[j].firstPoint; k<=leaves[j].lastPoint; k++) {
					mass+=points[k].getMass();
					centerOfMass+=points[k].position*points[k].getMass();
				}
				j++;
			}
			if ( ( i==0 ||  // Local maximum test
			       points[leaves[i].firstPoint].position-points[leaves[i-1].lastPoint].position>=minLocalMaxDistance ||
			       leaves[i].value>leaves[i-1].value 
				 ) &&
			     ( j>lastLeaf ||
			       points[leaves[j].firstPoint].position-points[leaves[j-1].lastPoint].position>=minLocalMaxDistance ||
			       leaves[i].value>leaves[j].value 
				 )
			   ) {
			   	out++;
				centerOfMass=integerPositions?Math.round(centerOfMass/mass):centerOfMass/mass;
				inLeaf=false;
				if (j==i+1) {
					leaves[i].isLocalMaximum=true;
					inLeaf=true;
				}
				else {
					for (k=i; k<j; k++) {
						epsilon=EPSILON_THRESHOLD*(points[leaves[k].lastPoint].position-points[leaves[k].firstPoint].position);
						if (points[leaves[k].firstPoint].position-(k==i?epsilon:0)<=centerOfMass && centerOfMass<=points[leaves[k].lastPoint].position+(k==j-1?epsilon:0)) {
							leaves[k].isLocalMaximum=true;
							inLeaf=true;
							break;
						}
					}
				}
				if (!inLeaf) {  // The center of mass might fall between two consecutive leaves
					found=false;
					for (k=i; k<j-1; k++) {
						if (points[leaves[k].lastPoint].position<=centerOfMass && centerOfMass<=points[leaves[k+1].firstPoint].position) {
							leaves[k].isLocalMaximum=true;
							found=true;
							break;
						}
					}
					if (IO.CONSISTENCY_CHECKS && !found) {
						System.err.println("setLocalMaxMin> ERROR: center of mass "+centerOfMass+" not found among the leaves of a local-max run: ["+i+".."+(j-1)+"]");
						System.err.println("setLocalMaxMin> leaves:");
						for (int x=0; x<=lastLeaf; x++) System.err.println(leaves[x]);
						System.exit(1);
					}
				}
			}
			else if ( ( i>0 &&   // Local minimum test
						( points[leaves[i].firstPoint].position-points[leaves[i-1].lastPoint].position<minLocalMaxDistance &&
					      leaves[i].value<leaves[i-1].value
					    )
					  ) &&
					  ( j<=lastLeaf &&
					  	( points[leaves[j].firstPoint].position-points[leaves[j-1].lastPoint].position<minLocalMaxDistance &&
					      leaves[i].value<leaves[j].value
					    )
					  )
					) {
				for (k=i; k<j; k++) leaves[k].isLocalMinimum=true;
			}
			i=j;
		}
		return out;
	}


	/**
	 * Removes leaves whose mass is smaller than $MIN_LEAF_MASS$, if they are not local
	 * maxima and if they have not been added by procedure $addPointsWithLargestMass$.
	 *
	 * Remark: the procedure assumes that $setLocalMaxMin$ has already been executed.
	 */
	public static final void discardSmallLeaves(Point[] points) {
		int i, j, k;
		int mass;
		Leaf tmp;

		i=-1;
		for (j=0; j<=lastLeaf; j++) {
			mass=0;
			for (k=leaves[j].firstPoint; k<=leaves[j].lastPoint; k++) mass+=points[k].getMass();
			if (mass>=MIN_LEAF_MASS || (leaves[j].isLocalMaximum && mass>=minLocalMaxMass) || leaves[j].isPointWithLargestMass) {
				i++;
				if (i!=j) {
					tmp=leaves[i];
					leaves[i]=leaves[j];
					leaves[j]=tmp;
				}
			}
			else {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("discardSmallLeaves> discarded leaf "+leaves[j]);
			}
		}
		lastLeaf=i;
	}


	/**
	 * Assigns $distanceFromPrevious$ to every leaf.
	 * The procedure assumes that leaves are sorted by their first point.
	 */
	public static final void computeDistances(Point[] points) {
		leaves[0].distanceFromPrevious=0;
		for (int i=1; i<=lastLeaf; i++) leaves[i].distanceFromPrevious=points[leaves[i].firstPoint].position-points[leaves[i-1].lastPoint].position;
	}


	/**
	 * Moves all local-maximum leaves to the beginning of array $leaves$, followed by
	 * all local-minimum leaves (see $Leaf.compareTo$). Then, marks the last element of
	 * every run of local-maximum leaves that are too close to one another (as decided by
	 * threshold $minLocalMaxDistance$), and such that they are all high or all not high.
	 * Local maxima that are of the same high/not-high type and that are too close to one
	 * another are likely the result of oversplitting a cluster.
	 *
	 * Remark: the procedure assumes that there is at least one local-maximum leaf, and
	 * that local-maximum leaves have already been marked as high/not-high.
	 *
	 * Remark: the procedure stores in $nRuns$ the number of runs of local-maximum leaves,
	 * and in $lastLocalMaximum$ and $lastLocalMinimum$ the position in $leaves$ of the
	 * last local maximum/minimum after sorting.
	 */
	public static final void markRunsOfLocalMaximumLeaves(Point[] points) {
		boolean isHigh;
		int i, firstLeaf;

		if (lastLeaf>0) Arrays.sort(leaves,0,lastLeaf+1);
		for (i=0; i<=lastLeaf; i++) leaves[i].marked=false;
		firstLeaf=0; nRuns=0; isHigh=leaves[0].isHigh;
		for (i=1; i<=lastLeaf; i++) {
			if (!leaves[i].isLocalMaximum) break;
			if ( leaves[i].isHigh==isHigh &&
				 points[leaves[i].lastPoint].position-points[leaves[firstLeaf].firstPoint].position<minLocalMaxDistance
			   ) continue;
			leaves[i-1].marked=true;
			nRuns++;
			firstLeaf=i;
			isHigh=leaves[i].isHigh;
		}
		leaves[i-1].marked=true;
		nRuns++;
		lastLocalMaximum=i-1;
		while (i<=lastLeaf) {
			if (!leaves[i].isLocalMinimum) break;
			i++;
		}
		lastLocalMinimum=i-1;
	}


	/**
	 * @return the mass of the heaviest run of local-maximum leaves. This is the sum of
	 * the masses of all events whose position lies between the first event of the first
	 * local-maximum peak in the run, and the last event in the last local-maximum peak in
	 * the run.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	public static final int getMassOfHeaviestRun(Point[] points) {
		int i, k;
		int mass, maxMass, firstLocalMaximum, firstPoint, lastPoint;

		maxMass=0;
		firstLocalMaximum=0;
		k=0;  // Pointer to events
		for (i=0; i<=lastLocalMaximum; i++) {
			if (!leaves[i].marked) continue;
			firstPoint=leaves[firstLocalMaximum].firstPoint;
			lastPoint=leaves[i].lastPoint;
			while (k<firstPoint) k++;
			mass=0;
			while (k<=lastPoint) mass+=points[k++].getMass();
			if (mass>maxMass) maxMass=mass;
			firstLocalMaximum=i+1;
		}
		return maxMass;
	}

	
	/**
	 * @return the mass of the $run$-th run of local-maximum leaves. This is the sum of
	 * the masses of all events whose position lies between the first event of the first
	 * local-maximum peak in the run, and the last event in the last local-maximum peak in
	 * the run.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	public static final int getMassOfRun(int run, Point[] points) {
		int i, currentRun, out;

		currentRun=-1; out=0;
		for (i=0; i<=lastLeaf; i++) {
			if (!leaves[i].isLocalMaximum) break;
			out+=Points.mass(points,leaves[i].firstPoint,leaves[i].lastPoint);
			if (leaves[i].marked) {
				currentRun++;
				if (currentRun==run) return out;
				out=0;
			}
		}
		return out;
	}
	
	
	/**
	 * Identical to $getMassOfRun()$, but returns the mass of the $run$-th high run.
	 */
	public static final int getMassOfHighRun(int run, Point[] points) {
		boolean isHigh;
		int i, currentRun, out;

		currentRun=-1; out=0; isHigh=false;
		for (i=0; i<=lastLeaf; i++) {
			if (!leaves[i].isLocalMaximum) break;
			isHigh|=leaves[i].isHigh;
			out+=Points.mass(points,leaves[i].firstPoint,leaves[i].lastPoint);
			if (leaves[i].marked) {
				if (isHigh) {
					currentRun++;
					if (currentRun==run) return out;
				}
				out=0; isHigh=false;
			}
		}
		return out;
	}


	/**
	 * @param lastPoint used only if $threshold>0$;
	 * @return the center of mass of the run $run$ of local-maximum leaves, where $run$ is
	 * indexed from zero and assumed to be valid.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	public static final double getCenterOfMassOfRun(int run, Point[] points, int lastPoint, double threshold) {
		int i;
		int currentRun, firstLocalMaximum;
		double out;

		out=-1;
		currentRun=-1;
		firstLocalMaximum=0;
		i=0;
		while (currentRun<run) {
			if (!leaves[i].marked) {
				i++;
				continue;
			}
			currentRun++;
			if (currentRun!=run) {
				i++;
				firstLocalMaximum=i;
				continue;
			}
			out=Points.getCenterOfMass(points,leaves[firstLocalMaximum].firstPoint,leaves[i].lastPoint,false,threshold,lastPoint);
			break;
		}
		return out;
	}
	
	
	/**
	 * @return $out[0]$: the closest center of mass of a high run of local-maximum leaves, 
	 * with respect to $position$. The search proceeds from the left (if $direction=true$)
	 * or right (if $direction=false$); -1 if no high run can be found.
	 * $out[1]$: maximum density of a high local-maximum leaf inside the run.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	public static final void getCenterOfMassOfHighRun(double position, boolean direction, Point[] points, double[] out) {
		boolean isHigh;
		int i;
		int firstLocalMaximum;
		double centerOfMass, maxDensity, distance, minDistance, minCenterOfMass, minDensity;

		centerOfMass=-1; maxDensity=-1; isHigh=false; minDensity=-1;
		minDistance=Math.POSITIVE_INFINITY; minCenterOfMass=-1;
		if (direction) {
			firstLocalMaximum=0;
			i=0;
			while (i<=lastLocalMaximum) {
				if (leaves[i].isHigh) {
					isHigh=true;
					if (leaves[i].value>maxDensity) maxDensity=leaves[i].value;
				}
				if (!leaves[i].marked) {
					i++;
					continue;
				}
				if (!isHigh) {
					i++;
					firstLocalMaximum=i;
					isHigh=false; maxDensity=-1;
					continue;
				}
				centerOfMass=Points.getCenterOfMass(points,leaves[firstLocalMaximum].firstPoint,leaves[i].lastPoint,false,-1,-1);
				distance=Math.abs(centerOfMass-position);
				if (distance<minDistance) {
					minDistance=distance;
					minCenterOfMass=centerOfMass;
					minDensity=maxDensity;
				}
				else break;
				i++;
				firstLocalMaximum=i;
				isHigh=false; maxDensity=-1;
			}
		}
		else {
			firstLocalMaximum=lastLocalMaximum;
			isHigh=leaves[firstLocalMaximum].isHigh;
			maxDensity=isHigh?leaves[firstLocalMaximum].value:-1;
			i=lastLocalMaximum-1;
			while (i>=-1) {
				if (i==-1 || leaves[i].marked) {
					if (isHigh) {
						centerOfMass=Points.getCenterOfMass(points,leaves[i+1].firstPoint,leaves[firstLocalMaximum].lastPoint,false,-1,-1);						
						distance=Math.abs(centerOfMass-position);
						if (distance<minDistance) {
							minDistance=distance;
							minCenterOfMass=centerOfMass;
							minDensity=maxDensity;
						}
						else break;
						if (i>=0) {
							firstLocalMaximum=i;
							isHigh=leaves[firstLocalMaximum].isHigh;
							maxDensity=isHigh?leaves[firstLocalMaximum].value:-1;
						}
						i--;
						continue;
					}
					if (i>=0) {
						firstLocalMaximum=i;
						isHigh=leaves[firstLocalMaximum].isHigh;
						maxDensity=isHigh?leaves[firstLocalMaximum].value:-1;
					}
					i--;
					continue;
				}
				if (leaves[i].isHigh) {
					isHigh=true;
					if (leaves[i].value>maxDensity) maxDensity=leaves[i].value;
				}
				i--;
			}
		}
		out[0]=minCenterOfMass;
		out[1]=minDensity;
	}


	/**
	 * Says not to cluster all points in $[firstPoint..lastPoint]$ that fall inside those
	 * local-minimum leaves that occur between two runs of local-maximum leaves. We can't
	 * discard the points inside all local-minimum leaves for clustering, since a local-
	 * minimum leaf could be located inside a cluster that has been oversplit.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed, and that all points in the range have $isUsedForClustering=true$.
	 */
	public static final void discardLocalMinima(Point[] points, int firstPoint, int lastPoint) {
		int i, j, k;

		j=0;  // Pointer to local-minimum leaves
		j=lastLocalMaximum+1;
		if (j>lastLeaf) return;
		k=firstPoint;  // Pointer to events
		for (i=1; i<=lastLeaf; i++) {
			if (!leaves[i].isLocalMaximum) break;
			if (leaves[i].marked || !leaves[i-1].marked) continue;
			while (j<=lastLocalMinimum && leaves[j].firstPoint<leaves[i-1].lastPoint) j++;
			while (j<=lastLocalMinimum && leaves[j].firstPoint<leaves[i].firstPoint) {
				while (k<leaves[j].firstPoint) k++;
				while (k<=leaves[j].lastPoint) {
					points[k].isUsedForClustering=false;
					k++;
				}
				j++;
			}
		}
	}


	/**
	 * Decides a point of separation between the $run$-th run of local-maximum leaves and
	 * the $run+1$-th run, or between the $run$-th run and $lastPoint$ (if
	 * $run=nRuns-1$), or between $firstPoint$ and the first run (if $run=-1$).
	 * The point of separation is set to the center (if $gapPosition=TRUE$) or to the
	 * first position (if $gapPosition=FALSE$) of the leftmost gap of size at least
	 * $minLocalMaxDistance$ between the two maxima, or to the center of mass (or to the
	 * first position) of the leftmost local-minimum leaf between the two maxima, or to 
	 * the center of mass (or to the first position) of all points between the two maxima.
	 *
	 * If $run=-1$, or if $run$ equals the rightmost run, and if there is no gap and no
	 * local-minimum leaf in the corresponding intervals, the position of the farthest
	 * point inside a leaf adjacent to the run is used as boundary.
	 *
	 * If $run+1$ is the rightmost run, if there is no local-minimum leaf between the
	 * runs, and if the leftmost gap between the runs is adjacent to $run+1$, then the
	 * position of the rightmost point inside a leaf adjacent to $run$ is used as
	 * boundary. This is because the rightmost run can be very far from the rest of the
	 * distribution, and using the gap would discard all points between the two runs.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 *
	 * @param $run \in [-1..nRuns-1]$.
	 */
	public static final double separateRun(int run, Point[] points, int firstPoint, int lastPoint, boolean gapPosition) {
		int i, x, y, lastLeafOfRun, minLeaf, maxLeaf;
		double j, from, to;

		if (run==-1) {
			lastLeafOfRun=-1;
			from=points[firstPoint].position;
			to=points[leaves[0].firstPoint].position;
		}
		else {
			lastLeafOfRun=Leaves.lastLeafOfRun(leaves,lastLeaf,run);
			if (leaves[lastLeafOfRun].lastPoint==lastPoint) return points[lastPoint].position;
			from=points[leaves[lastLeafOfRun].lastPoint].position;
			if (run<nRuns-1) to=points[leaves[lastLeafOfRun+1].firstPoint].position;
			else to=points[lastPoint].position;
		}
		i=findLocalMinimumLeaf(points,from,to);

if (IO.SHOW_STD_ERR) {
	IO.printErr("separateRun> run="+run+" firstPoint="+firstPoint+" lastPoint="+lastPoint+" i="+i+" number of runs: "+nRuns+"  from="+from+" to="+to+" lastLeafOfRun="+lastLeafOfRun);
	if (i>=0) System.err.println("local minimum leaf: ["+points[leaves[i].firstPoint].position+".."+points[leaves[i].lastPoint].position+"]");
}

		if (i>=0) return gapPosition?Points.getCenterOfMass(points,leaves[i].firstPoint,leaves[i].lastPoint,true,-1,-1):points[leaves[i].firstPoint].position;
		i=findGap(points,from,to);
		if (i>=0 && (run<nRuns-2 || i!=lastLeafOfRun+1)) {
			j=points[leaves[i].firstPoint].position;
			return j-(gapPosition?leaves[i].distanceFromPrevious/2:leaves[i].distanceFromPrevious);
		}

		// Neither a local-minimum leaf nor a valid gap can be found
		if (run==-1) {
			maxLeaf=-1;
			for (i=1; i<=lastLeaf; i++) {
				if (leaves[i].lastPoint==leaves[0].firstPoint-1) {
					maxLeaf=i;
					break;
				}
			}
			if (maxLeaf==-1 || leaves[maxLeaf].firstPoint==firstPoint) return points[leaves[0].firstPoint].position;
			return points[leaves[maxLeaf].firstPoint].position;
		}
		else {
			x=leaves[lastLeafOfRun].lastPoint+1;
			if (run==nRuns-1) {
				minLeaf=-1;
				for (i=0; i<=lastLeaf; i++) {
					if (leaves[i].firstPoint==x) {
						minLeaf=i;
						break;
					}
				}
				if (minLeaf==-1 || leaves[minLeaf].lastPoint==lastPoint) return points[leaves[lastLeafOfRun].lastPoint].position;
				return gapPosition?points[leaves[minLeaf].lastPoint].position:points[leaves[lastLeafOfRun].lastPoint].position;
			}
			else {
				minLeaf=-1;
				for (i=0; i<=lastLeaf; i++) {
					if (leaves[i].firstPoint==x) {
						minLeaf=i;
						break;
					}
				}
				if (minLeaf==-1 || minLeaf==lastLeafOfRun+1) return gapPosition?(points[leaves[lastLeafOfRun].lastPoint].position+points[leaves[lastLeafOfRun+1].firstPoint].position)/2:points[leaves[lastLeafOfRun].lastPoint].position;
				return gapPosition?points[leaves[minLeaf].lastPoint].position:points[leaves[lastLeafOfRun].lastPoint].position;
			}
		}
	}


	/**
	 * Returns the position in $leaves$ of the leftmost local-minimum leaf whose points
	 * are fully contained in $[from..to]$, or -1 if no such leaf can be found.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	private static final int findLocalMinimumLeaf(Point[] points, double from, double to) {
		int i;
		double j;

if (IO.SHOW_STD_ERR) IO.printErr("findLocalMinimumLeaf> lastLocalMaximum="+lastLocalMaximum+" lastLocalMinimum="+lastLocalMinimum+" from="+from+" to="+to);

		i=lastLocalMaximum+1;
		while (i<=lastLocalMinimum) {
			j=points[leaves[i].firstPoint].position;
			if (j<from || j>to) {
				i++;
				continue;
			}
			j=points[leaves[i].lastPoint].position;
			if (j<from || j>to) {
				i++;
				continue;
			}
			return i;
		}
		return -1;
	}


	/**
	 * Returns the position of the leftmost leaf in $leaves$ in string order that is at
	 * distance at least $minLocalMaxDistance$ from the previous leaf in string order, and
	 * such that the gap between the leaf and the previous leaf is fully contained inside
	 * $[from..to]$. The procedure returns -1 if no such leaf can be found.
	 *
	 * Remark: local-minimum leaves are not considered.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	private static final int findGap(Point[] points, double from, double to) {
		int i, minLeaf;
		double j, minStart;

		if (to-from<minLocalMaxDistance) return -1;
		minLeaf=-1; minStart=Math.POSITIVE_INFINITY;
		for (i=0; i<=lastLocalMaximum; i++) {
			if (leaves[i].distanceFromPrevious<minLocalMaxDistance) continue;
			j=points[leaves[i].firstPoint].position;
			if (j<from || j>to) continue;
			if (j-leaves[i].distanceFromPrevious<from) continue;
			if (j<minStart) {
				minStart=j;
				minLeaf=i;
			}
		}
		for (i=lastLocalMinimum+1; i<=lastLeaf; i++) {
			if (leaves[i].distanceFromPrevious<minLocalMaxDistance) continue;
			j=points[leaves[i].firstPoint].position;
			if (j<from || j>to) continue;
			if (j-leaves[i].distanceFromPrevious<from) continue;
			if (j<minStart) {
				minStart=j;
				minLeaf=i;
			}
		}
		return minLeaf;
	}


	/**
	 * Assume that $points$ have integer positions, and that some points have very large
	 * mass. The DET could assign such a point to a long leaf with overall smaller
	 * density, and the density of such leaf could be smaller than the density of other
	 * leaves (which e.g. happen to be shorter). Thus, points with very large mass could
	 * be eventually missed, both as local maxima and as high local maxima.
	 * The procedure assigns every point whose mass is at least $largestMass$ to a new
	 * leaf of length one, splitting in two the original leaf that contains it. This
	 * happens only if the leaf that originally contains the point is not already short 
	 * and not already a local maximum.
	 *
	 * Remark: the new leaf, or the containing small leaf, is assigned
	 * $isPointWithLargestMass=true$.
	 *
	 * Remark: there could be runs of multiple consecutive points with mass at least
	 * $largestMass$, however such points could belong to distinct leaves. The
	 * procedure does not handle such runs in an ad hoc way.
	 *
	 * Remark: the procedure assumes that points are located at integer positions, and
	 * that $leaves$ and $points$ are distinct and sorted by position.
	 *
	 * @param largestMass typically the output of procedure $Points.largestMass()$;
	 * @return TRUE iff new points were added.
	 */
	private static final boolean addPointsWithLargestMass(Point[] points, int firstPoint, int lastPoint, int largestMass) {
		final int MIN_LEAF_LENGTH = 5;
		int i, j, k, p;
		double oldLength, newLength;
		double oldMass, rightMass, pointMass;
		Leaf tmp;

		// Adding single-point leaves
		j=lastLeaf; k=lastLeaf;
		for (i=lastPoint; i>=firstPoint; i--) {
			pointMass=points[i].getMass();
			if (pointMass<largestMass) continue;
			while (j>=0 && k<leaves.length-1) {
				if (i>=leaves[j].firstPoint && i<=leaves[j].lastPoint) {
					oldLength=points[leaves[j].lastPoint].position-points[leaves[j].firstPoint].position+1;
					if (oldLength<MIN_LEAF_LENGTH) {
						leaves[j].isPointWithLargestMass=true;
						break;
					}
					if (leaves[j].isLocalMaximum) break;
					oldMass=leaves[j].value*oldLength;
					if (i==leaves[j].firstPoint) {
						newLength=points[leaves[j].lastPoint].position-points[leaves[j].firstPoint+1].position+1;
						leaves[j].value=(oldMass-pointMass)/newLength;
						leaves[j].firstPoint++;
					}
					else if (i==leaves[j].lastPoint) {
						newLength=points[leaves[j].lastPoint-1].position-points[leaves[j].firstPoint].position+1;
						leaves[j].value=(oldMass-pointMass)/newLength;
						leaves[j].lastPoint--;
					}
					else {
						newLength=points[i-1].position-points[leaves[j].firstPoint].position+1;
						rightMass=0.0;
						for (p=i+1; p<=leaves[j].lastPoint; p++) rightMass+=points[p].getMass();
						leaves[j].value=(oldMass-pointMass-rightMass)/newLength;
						newLength=points[leaves[j].lastPoint].position-points[i+1].position+1;
						k++;
						leaves[k].firstPoint=i+1;
						leaves[k].lastPoint=leaves[j].lastPoint;
						leaves[k].value=rightMass/newLength;
						leaves[j].lastPoint=i-1;
					}
					k++;
					leaves[k].firstPoint=i;
					leaves[k].lastPoint=i;
					leaves[k].value=pointMass;
					leaves[k].isPointWithLargestMass=true;
					break;
				}
				j--;
			}
		}
		if (k==lastLeaf) return false;

		// Sorting
		for (i=lastLeaf+1; i<=k; i++) {
			for (j=i-1; j>=0; j--) {
				if (leaves[j].firstPoint<leaves[j+1].firstPoint) break;
				tmp=leaves[j];
				leaves[j]=leaves[j+1];
				leaves[j+1]=tmp;
			}
		}
		lastLeaf=k;
		return true;
	}


	/**
	 * @return the ID of the leftmost (if $direction=true$) or rightmost (if $direction=
	 * false$) run of local-maximum leaves whose distance from the previous run is
	 * maximum.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	public static final int getRunWithLargestDistance(Point[] points, boolean direction) {
		int i;
		int currentRun, maxRunLeft, maxRunRight, firstLocalMaximum, firstPoint;
		double distance, maxDistance;

		maxDistance=0.0; maxRunLeft=-1; maxRunRight=-1;
		currentRun=-1; firstLocalMaximum=0;
		for (i=0; i<=lastLocalMaximum; i++) {
			if (!leaves[i].marked) continue;
			currentRun++;
			firstPoint=leaves[firstLocalMaximum].firstPoint;
			distance=currentRun==0?0.0:points[leaves[firstLocalMaximum].firstPoint].position-points[leaves[firstLocalMaximum-1].lastPoint].position;
			if (distance>maxDistance) {
				maxDistance=distance;
				maxRunLeft=currentRun;
				maxRunRight=currentRun;
			}
			else if (distance==maxDistance) maxRunRight=currentRun;
			firstLocalMaximum=i+1;
		}
		return direction?maxRunLeft:maxRunRight;
	}


	/**
	 * Partitions a leaf into distinct leaves if it contains two consecutive points whose
	 * distance is at least $minHoleLength$. This is useful, since a leaf of the DET might
	 * contain a peak as well as parts of an adjacent peak, because splitting might have
	 * stopped too early.
	 *
	 * Remark: the procedure assumes $leaves$ to be the set of leaves produced by the
	 * sequence of recursive splits.
	 *
	 * @param smallLeafDistance short leaves created by hole removal are merged to 
	 * neighboring leaves if they are at distance at most $smallLeafDistance$ (this is not
	 * done if $smallLeafDistance < 0$).
	 */
	private static final void removeHoles(double minHoleLength, Point[] points, double smallLeafDistance) {
		int i, j, k;
		int currentLeaf, originalFirstPoint, originalLastPoint;
		int previousLeafMass, currentLeafMass, nextLeafMass;
		double currentMass, currentLength, distanceLeft, distanceRight;
		Leaf tmp;

		// Splitting leaves at holes
		for (i=0; i<=lastLeaf; i++) leaves[i].fromHole=0;
		k=lastLeaf;
		for (i=0; i<=lastLeaf; i++) {
			currentLeaf=i;
			originalFirstPoint=leaves[i].firstPoint;
			originalLastPoint=leaves[i].lastPoint;
			if (originalLastPoint==originalFirstPoint) continue;
			currentMass=points[leaves[i].firstPoint].getMass();
			j=originalFirstPoint+1;
			while (j<=originalLastPoint) {
				if (points[j].position-points[j-1].position>=minHoleLength) {
					leaves[currentLeaf].lastPoint=j-1;
					currentLength=points[j-1].position-points[leaves[currentLeaf].firstPoint].position+(integerPositions?1.0:0.0);
					leaves[currentLeaf].value=currentMass/currentLength;
					leaves[currentLeaf].fromHole=1;
					k++;
					leaves[k].firstPoint=j;
					leaves[k].fromHole=1;
					currentLeaf=k;
					currentMass=points[j].getMass();
				}
				else currentMass+=points[j].getMass();
				j++;
			}
			if (currentLeaf!=i) {
				leaves[currentLeaf].lastPoint=originalLastPoint;
				currentLength=points[leaves[currentLeaf].lastPoint].position-points[leaves[currentLeaf].firstPoint].position+(integerPositions?1.0:0.0);
				leaves[currentLeaf].value=currentMass/currentLength;
			}
		}
		if (k==lastLeaf) return;


if (IO.SHOW_STD_ERR_PRIME) IO.printErr("removeHoles> added leaves:");
if (IO.SHOW_STD_ERR_PRIME) { for (i=lastLeaf+1; i<=k; i++) IO.printErr(leaves[i]); }


		// Sorting
		for (i=lastLeaf+1; i<=k; i++) {
			for (j=i-1; j>=0; j--) {
				if (leaves[j].firstPoint<leaves[j+1].firstPoint) break;
				tmp=leaves[j];
				leaves[j]=leaves[j+1];
				leaves[j+1]=tmp;
			}
		}
		lastLeaf=k;
		
		// Trying to merge short leaves created by hole removal
		if (lastLeaf==0 || smallLeafDistance<0) return;
		previousLeafMass=0; 
		currentLeafMass=leaves[0].getMass(points); 
		nextLeafMass=leaves[1].getMass(points);
		for (i=0; i<=lastLeaf; i++) {
			if (leaves[i].fromHole==0) {
				previousLeafMass=currentLeafMass;
				currentLeafMass=nextLeafMass;
				nextLeafMass=i+2<=lastLeaf?leaves[i+2].getMass(points):0;
				continue;
			}
			distanceRight=i+1<=lastLeaf?points[leaves[i+1].firstPoint].position-points[leaves[i].firstPoint].position:Math.POSITIVE_INFINITY;
			distanceLeft=i>0?points[leaves[i].lastPoint].position-points[leaves[i-1].lastPoint].position:Math.POSITIVE_INFINITY;
			if (distanceRight<=smallLeafDistance && distanceLeft>smallLeafDistance && leaves[i+1].fromHole==0) {
				leaves[i].fromHole=-1;
				leaves[i+1].firstPoint=leaves[i].firstPoint;
				currentLength=points[leaves[i+1].lastPoint].position-points[leaves[i].firstPoint].position+(integerPositions?1.0:0.0);
				leaves[i+1].value=(currentLeafMass+nextLeafMass)/currentLength;
			}
			else if (distanceLeft<=smallLeafDistance && distanceRight>smallLeafDistance && leaves[i-1].fromHole==0) {
				leaves[i].fromHole=-1;
				leaves[i-1].lastPoint=leaves[i].lastPoint;
				currentLength=points[leaves[i].lastPoint].position-points[leaves[i-1].firstPoint].position+(integerPositions?1.0:0.0);
				leaves[i-1].value=(previousLeafMass+currentLeafMass)/currentLength;
			}
			previousLeafMass=currentLeafMass;
			currentLeafMass=nextLeafMass;
			nextLeafMass=i+2<=lastLeaf?leaves[i+2].getMass(points):0;
		}
		
		// Discarding merged leaves
		j=-1;
		for (i=0; i<=lastLeaf; i++) {
			if (leaves[i].fromHole==-1) continue;
			j++;
			tmp=leaves[j];
			leaves[j]=leaves[i];
			leaves[i]=tmp;
		}
		lastLeaf=j;
	}
	
	
	/**
	 * Sets $highRunTmp[0],highRunTmp[1]$ to the first and last leaf of a high run such
	 * that: (1) its rightmost point is at least $min$; (2) its mass is at least 
	 * $minMass$; (3) it has the largest mass; (4) it can (if $includeFirstRun=true$) or
	 * cannot be (if $includeFirstRun=false$) the first high run.
	 * The elements of $highRunTmp$ are set to -1 if the desired high run does not exist.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 *
	 * @return the ratio between the largest mass and the second-largest mass of a high 
	 * run with the desired properties, or zero if at most one such run can be found.
	 */
	public static final double findHighRun(double min, Point[] points, int minMass, boolean includeFirstRun) {
		boolean isHigh;
		int i;
		int first, mass, largestMass, secondLargestMass, largestRunFirst, largestRunLast;
		
		// Finding the first high run
		findFirstHighRun(points);
		first=includeFirstRun?highRunTmp[0]:highRunTmp[1]+1;
		highRunTmp[0]=-1;
		highRunTmp[1]=-1;
		
		// Finding the desired high run starting from $first$.
		largestMass=0; secondLargestMass=0; largestRunFirst=-1; largestRunLast=-1;
		i=-1;
		while (first<=lastLocalMaximum) {
			i=first; isHigh=false; mass=0;
			while (!leaves[i].marked) {
				isHigh|=leaves[i].isHigh;
				mass+=leaves[i].getMass(points);
				i++;
			}
			isHigh|=leaves[i].isHigh;
			mass+=leaves[i].getMass(points);
			if (isHigh && points[leaves[i].lastPoint].position>=min && mass>=minMass && mass>largestMass) {
				secondLargestMass=largestMass; largestMass=mass;
				largestRunFirst=first;
				largestRunLast=i;
			}
			first=i+1;
		}
		if (largestRunFirst>0 && largestRunLast>0) {
			highRunTmp[0]=largestRunFirst;
			highRunTmp[1]=largestRunLast;
			return secondLargestMass==0?0:((double)largestMass)/secondLargestMass;
		}
		else return 0;
	}
	
	
	/**
	 * Sets $highRunTmp[0],highRunTmp[1]$ to the first and last leaf of the first high 
	 * run. The elements of $highRunTmp$ are set to -1 if there is no high run.
	 *
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	public static final void findFirstHighRun(Point[] points) {
		boolean isHigh;
		int i, first, mass;
		
		highRunTmp[0]=-1;
		highRunTmp[1]=-1;
		first=0; i=-1;
		while (first<=lastLocalMaximum) {
			i=first; isHigh=false;
			while (!leaves[i].marked) {
				isHigh|=leaves[i].isHigh;
				i++;
			}
			isHigh|=leaves[i].isHigh;
			if (isHigh) break;
			first=i+1;
		}
		if (first>lastLocalMaximum) return;
		highRunTmp[0]=first;
		highRunTmp[1]=i;
	}

	
	/**
	 * Remark: the procedure assumes that $markRunsOfLocalMaximumLeaves$ has already been
	 * executed.
	 */
	public static final int getNumberOfRuns() {
		int i, count;
		
		count=0;
		for (i=0; i<=lastLocalMaximum; i++) {
			if (leaves[i].marked) count++;
		}
		return count;
	}

}