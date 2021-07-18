package de.mpi_cbg.revant.util;

import java.util.Arrays;
import de.mpi_cbg.revant.util.Math;


/**
 * See e.g.:
 * Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. The elements of statistical
 * learning. Springer, Berlin: Springer series in statistics, 2001.
 */
public class RegressionTree {
	/**
	 * Parameters of the pipeline
	 */
	public static int minIntervalLength = IO.quantum;  // Intervals of length smaller than or equal to this are not split
	public static int minLocalMaxDistance;  // Minimum distance between two local maxima for them to be considered distinct peaks
	public static double splitErrorRatio = 3.0;  // An interval can be split into intervals of length smaller than $minIntervalLength$, iff the new variance is at most the old variance divided by this number.

	/**
	 * Output data structures
	 */
	public static Leaf[] leaves;  // The sorted list of leaves
	public static int lastLeaf;  // Position of the last leaf in $leaves$
	public static int nRuns;  // Number of runs of local-maximum leaves
	public static int lastLocalMaximum;  // Position of the last local-maximum leaf in $leaves$ after $markRunsOfLocalMaximumLeaves$.
	public static int lastLocalMinimum;  // Position of the last local-minimum leaf in $leaves$ after $markRunsOfLocalMaximumLeaves$.

	/**
	 * Temporary data structures
	 */
	private static double[] stack;  // Temporary space used to traverse the tree
	private static double[] tmp;


	public static final void allocateMemory(int maxLeaves, int mlmd) {
		int i;
		
		minLocalMaxDistance=mlmd;
		leaves = new Leaf[maxLeaves];
		for (i=0; i<maxLeaves; i++) leaves[i] = new Leaf();
		stack = new double[maxLeaves<<2];
		tmp = new double[6];
	}


	/**
	 * Builds the intervals of a regression tree on values $x[firstPoint..lastPoint]$,
	 * storing its leaves in array $leaves[0..lastLeaf]$. Intervals with standard
	 * deviation at most $maxStd$ are not split.
	 *
	 * Remark: the procedure assumes that $[firstPoint..lastPoint]$ does not have a long 
	 * prefix or a long suffix of zeros.
	 *
	 * @param step only multiples of $step$ are used as candidate split points of a node.
	 * Setting $step$ to values bigger than one is useful for speed;
	 * @param shortEnds allows intervals shorter than $minIntervalLength$ at the beginning
	 * or end of $x[firstPoint..lastPoint]$;
	 * @return the number of local-maximum leaves, if any.
	 */
	public static final int buildRegressionTree(double[] x, int firstPoint, int lastPoint, double maxStd, int step, double largestValue, boolean shortEnds, int overflowOperations) {
		int i, top, topPrime, length, out;
		int intervalFirst, intervalLast;
		double sum, expectation, variance, splitPoint;
		double intervalExpectation, intervalVariance;

		// Pushing the root onto $stack$
		length=lastPoint-firstPoint+1;
		sum=0.0;
		for (i=firstPoint; i<=lastPoint; i++) sum+=x[i];
		expectation=sum/length;
		sum=0.0;
		for (i=firstPoint; i<=lastPoint; i++) sum+=(x[i]-expectation)*(x[i]-expectation);
		variance=sum/length;
		stack[0]=firstPoint; stack[1]=lastPoint; stack[2]=expectation; stack[3]=variance;
		top=0; lastLeaf=-1;

		// Building the tree in preorder
		while (top>=0) {
			intervalFirst=(int)stack[top];
			intervalLast=(int)stack[top+1];
			intervalExpectation=stack[top+2];
			intervalVariance=stack[top+3];
			top-=4;

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("RegressionTree> considering interval ["+intervalFirst+".."+intervalLast+"] intervalExpectation="+intervalExpectation+" intervalVariance="+intervalVariance+" maxStd="+maxStd+" minIntervalLength="+minIntervalLength);

			if (!shouldSplit(x,intervalFirst,intervalLast,intervalExpectation,intervalVariance,maxStd)) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("RegressionTree> should NOT split");
				lastLeaf++;
				leaves[lastLeaf].firstPoint=intervalFirst;
				leaves[lastLeaf].lastPoint=intervalLast;
				leaves[lastLeaf].value=intervalExpectation;
				continue;
			}
			split(x,intervalFirst,intervalLast,step,overflowOperations);			
			splitPoint=tmp[0];
			if ( splitPoint==-1 ||
				 (tmp[2]>intervalVariance && tmp[4]>intervalVariance) ||					 
			     ( tmp[5]==0 && 
				   !( (tmp[2]*splitErrorRatio<=intervalVariance && tmp[4]*splitErrorRatio<=intervalVariance) || (shortEnds&&(intervalFirst==firstPoint||intervalLast==lastPoint)) )
				 )
			   ) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("RegressionTree> splitting failed");
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("reason: splitPoint="+splitPoint+" tmp2="+tmp[2]+" "+" tmp4="+tmp[4]+" tmp5="+tmp[5]+" intervalVariance="+intervalVariance+" intervalFirst="+intervalFirst+" intervalLast="+intervalLast+" step="+step);
				lastLeaf++;
				leaves[lastLeaf].firstPoint=intervalFirst;
				leaves[lastLeaf].lastPoint=intervalLast;
				leaves[lastLeaf].value=intervalExpectation;
				continue;
			}
			if (tmp[5]==1) {
				// Recurring on the right child
				topPrime=top+4;
				stack[topPrime]=splitPoint+1;
				stack[topPrime+1]=intervalLast;
				stack[topPrime+2]=tmp[3];
				stack[topPrime+3]=tmp[4];
				// Recurring on the left child
				topPrime+=4;
				stack[topPrime]=intervalFirst;
				stack[topPrime+1]=splitPoint;
				stack[topPrime+2]=tmp[1];
				stack[topPrime+3]=tmp[2];
				top=topPrime;
			}
			else {
				// Saving the left child
				lastLeaf++;
				leaves[lastLeaf].firstPoint=intervalFirst;
				leaves[lastLeaf].lastPoint=(int)splitPoint;
				leaves[lastLeaf].value=tmp[1];
				// Saving the right child
				lastLeaf++;
				leaves[lastLeaf].firstPoint=(int)splitPoint+1;
				leaves[lastLeaf].lastPoint=intervalLast;
				leaves[lastLeaf].value=tmp[3];
			}
		}

		// Initializing leaves
		if (lastLeaf==0) {
			leaves[0].isLocalMaximum=false;
			leaves[0].isLocalMinimum=false;
			leaves[0].marked=false;
			leaves[0].isHigh=false;
			return 0;
		}
		out=setLocalMaxMin(x);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("RegressionTree> leaves 1: ["+firstPoint+".."+lastPoint+"]");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastLeaf; i++) IO.printErr(leaves[i]); }	

		// Adding new leaves for points with high value
		for (i=0; i<=lastLeaf; i++) leaves[i].isPointWithLargestMass=false;
		if (largestValue>0) addPointsWithLargestValue(x,firstPoint,lastPoint,largestValue);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("largestValue="+largestValue+" lastLeaf="+lastLeaf);
		for (i=0; i<=lastLeaf; i++) {
			if (leaves[i].isPointWithLargestMass) leaves[i].isHigh=true;
			else leaves[i].isHigh=false;
		}
		
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("RegressionTree> leaves 2:");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastLeaf; i++) IO.printErr(leaves[i]); }		
		
		
		out=setLocalMaxMin(x);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("RegressionTree> leaves 3:");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastLeaf; i++) IO.printErr(leaves[i]); }

		return out;
	}


	/**
	 * Tries to split interval $[start..end]$ of a regression tree on values $x$, so that
	 * both of the resulting intervals have length at least $minIntervalLength$. If this
	 * is not possible, the procedure splits the interval at the global point that
	 * minimizes the variance.
	 *
	 * Remark: the procedure assumes that $[start..end]$ does not have a long prefix or a
	 * long suffix of zeros.
	 *
	 * Remark: for speed, the procedure considers just increments of $step$ as candidate
	 * split points.
	 *
	 * @return $tmp[0]$: the point $i$ that minimizes the sum of the variances of the two
	 * partitions $x[start..i],x[i+1..end]$. The intervals induced by splitting at this
	 * point are both of length at least $minIntervalLength$ if $tmp[5]=1$, they are not
	 * if $tmp[5]=0$.
	 * $tmp[1]$: the expectation of the left partition;
	 * $tmp[2]$: the variance of the left partition;
	 * $tmp[3]$: the expectation of the right partition;
	 * $tmp[4]$: the variance of the right partition.
	 */
	private static final void split(double[] x, int start, int end, int step, int overflowOperations) {
		final double LOWER_OVERFLOW_THRESHOLD = 1e-10;
		final double UPPER_OVERFLOW_THRESHOLD = 1e10;
		int i, j;
		int minPointAllowed, minPointGlobal;
		double sumLeft, sumRight;
		double expectationLeft, expectationRight;
		double minExpectationLeftAllowed, minExpectationLeftGlobal;
		double minExpectationRightAllowed, minExpectationRightGlobal;
		double variance, varianceLeft, varianceRight;
		double minVarianceAllowed, minVarianceGlobal;
		double minVarianceLeftAllowed, minVarianceLeftGlobal;
		double minVarianceRightAllowed, minVarianceRightGlobal;

		sumRight=0.0;
		for (i=start; i<=end; i++) sumRight+=x[i];
		sumLeft=0.0;
		minPointAllowed=-1; minPointGlobal=-1;
		minVarianceAllowed=Math.POSITIVE_INFINITY; minVarianceGlobal=Math.POSITIVE_INFINITY;
		minExpectationLeftAllowed=-1; minExpectationLeftGlobal=-1;
		minVarianceLeftAllowed=-1; minVarianceLeftGlobal=-1;
		minExpectationRightAllowed=-1; minExpectationRightGlobal=-1;
		minVarianceRightAllowed=-1; minVarianceRightGlobal=-1;
		for (i=start; i<end; i++) {
			sumLeft+=x[i];
			if (sumLeft<0 || sumLeft>=UPPER_OVERFLOW_THRESHOLD || (i-start+1>=overflowOperations && sumLeft!=0 && sumLeft<=LOWER_OVERFLOW_THRESHOLD)) {
				break;  // Overflow
			}
			sumRight-=x[i];
			if (sumRight<0 || sumRight>=UPPER_OVERFLOW_THRESHOLD || (i-start+1>=overflowOperations && sumRight!=0 && sumRight<=LOWER_OVERFLOW_THRESHOLD)) {
				break;  // Overflow
			}
			if ((i-start)%step!=0) continue;
			expectationLeft=sumLeft/(i-start+1);
			expectationRight=sumRight/(end-i);
			varianceLeft=0.0;
			for (j=start; j<=i; j++) varianceLeft+=(x[j]-expectationLeft)*(x[j]-expectationLeft);
			if (varianceLeft<0 || varianceLeft>=UPPER_OVERFLOW_THRESHOLD || (i-start+1>=overflowOperations && varianceLeft!=0 && varianceLeft<=LOWER_OVERFLOW_THRESHOLD)) {
				break;  // Overflow, or left interval almost constant (unlikely for so many operations).
			}
			varianceRight=0.0;
			for (j=i+1; j<=end; j++) varianceRight+=(x[j]-expectationRight)*(x[j]-expectationRight);
			if (varianceRight<0 || varianceRight>=UPPER_OVERFLOW_THRESHOLD || (i-start+1>=overflowOperations && varianceRight!=0 && varianceRight<=LOWER_OVERFLOW_THRESHOLD)) {
				break;  // Overflow, or right interval almost constant.
			}
			variance=varianceLeft+varianceRight;
			if (i-start+1>=minIntervalLength && end-i>=minIntervalLength && variance<minVarianceAllowed) {
				minVarianceAllowed=variance;
				minExpectationLeftAllowed=expectationLeft;
				minExpectationRightAllowed=expectationRight;
				minVarianceLeftAllowed=varianceLeft/(i-start+1);
				minVarianceRightAllowed=varianceRight/(end-i);
				minPointAllowed=i;
			}
			if (variance<minVarianceGlobal) {
				minVarianceGlobal=variance;
				minExpectationLeftGlobal=expectationLeft;
				minExpectationRightGlobal=expectationRight;
				minVarianceLeftGlobal=varianceLeft/(i-start+1);
				minVarianceRightGlobal=varianceRight/(end-i);
				minPointGlobal=i;
			}
		}

		// Outputting
		if (minPointAllowed!=-1) {
			tmp[0]=minPointAllowed;
			tmp[1]=minExpectationLeftAllowed;
			tmp[2]=minVarianceLeftAllowed;
			tmp[3]=minExpectationRightAllowed;
			tmp[4]=minVarianceRightAllowed;
			tmp[5]=1;
		}
		else {
			tmp[0]=minPointGlobal;
			tmp[1]=minExpectationLeftGlobal;
			tmp[2]=minVarianceLeftGlobal;
			tmp[3]=minExpectationRightGlobal;
			tmp[4]=minVarianceRightGlobal;
			tmp[5]=0;
		}
	}


	/**
	 * Decides whether the interval $[start..end]$ of a regression tree on values $x$,
	 * with average $expectation$ and with variance $variance$ inside the interval, should
	 * be split or not.
	 *
	 * Remark: we don't use a chi-squared test since the histogram has real values.
	 *
	 * @return TRUE if the size of the interval is bigger than $minIntervalLength$, and if
	 * the standard deviation is bigger than $maxStd$.
	 */
	private static final boolean shouldSplit(double[] x, int start, int end, double expectation, double variance, double maxStd) {
		final int WINDOW = minIntervalLength;  // Arbitrary
		final double RATIO = 3.0;  // Arbitrary
		
		if (end-start+1<=minIntervalLength) return false;
		if (Math.sqrt(variance)<=maxStd) return hasHighWindowRatio(x,start,end,WINDOW,RATIO);
		return true;
	}
	
	
	/**
	 * Checks whether a window of length $window$ of vector $x[start..end]$ has
	 * average bigger than $ratio$ times the average of an adjacent window of the same 
	 * length.
	 */
	private static final boolean hasHighWindowRatio(double[] x, int start, int end, int window, double ratio) {
		final int ONE_WINDOW = window;
		final int TWO_WINDOWS = (ONE_WINDOW)<<1;
		int i;
		double sumLeft, sumRight;

		if (end-start+1<TWO_WINDOWS) return false;
		sumLeft=0.0;
		for (i=start; i<start+ONE_WINDOW; i++) sumLeft+=x[i];
		sumRight=0.0;
		for (i=start+ONE_WINDOW; i<start+TWO_WINDOWS; i++) sumRight+=x[i];
		if (sumLeft>sumRight*ratio || sumRight>sumLeft*ratio) return true;
		
		i=start;
		while ((i+TWO_WINDOWS)<=end) {
			sumLeft+=x[i+ONE_WINDOW]-x[i];
			sumRight+=x[i+TWO_WINDOWS]-x[i+ONE_WINDOW];	
			if (sumLeft>sumRight*ratio || sumRight>sumLeft*ratio) return true;
			i++;
		}
		return false;
	}


	/**
	 * Sets fields $isLocalMaximum$ and $isLocalMinimum$ of each leaf in $leaves$.
	 *
	 * Let a run be a sequence of consecutive leaves with similar value. We say that a
	 * run is left (respectively, right) maximum if the value of the leaves in the run is
	 * bigger than the value of the leaves in the previous (respectively, next) adjacent
	 * run. We say that a run is left (respectively, right) minimum if the value of the
	 * leaves in the run is smaller than the value of the leaves in the previous
	 * (respectively, next) adjacent run. A run is a local maximum (minimum) if it is
	 * both left and right maximum (minimum). The procedure marks as local maximum 
	 * (minimum) all leaves in a local-maximum (-minimum) run.
	 *
	 * @return the number of local maxima found.
	 */
	private static final int setLocalMaxMin(double[] x) {
		final double THRESHOLD = 0.1;  // Arbitrary
		final double OUTLIER_RATIO = 3.0;  // Arbitrary
		final int LEAF_RATIO = 2;  // Arbitrary
		int i, j, k;
		int windowLeft_start, windowLeft_end, windowRight_start, windowRight_end, out;
		double lastValue;

		for (i=0; i<=lastLeaf; i++) {
			leaves[i].isLocalMaximum=false;
			leaves[i].isLocalMinimum=false;
		}
		out=0;
		if (lastLeaf==0) return out;
		i=0;
		while (i<=lastLeaf) {
			j=i+1; lastValue=leaves[i].value;
			while ( j<=lastLeaf &&
			        Math.abs(leaves[j].value-leaves[i].value)<=Math.max(leaves[j].value,leaves[i].value)*THRESHOLD
			      ) {
				lastValue=leaves[j].value;
				j++;
			}
			windowLeft_start=i>0?leaves[i].firstPoint-Math.min(minIntervalLength,(leaves[i-1].lastPoint-leaves[i-1].firstPoint+1)/LEAF_RATIO):-1;
			windowLeft_end=leaves[i].firstPoint+Math.min(minIntervalLength,(leaves[i].lastPoint-leaves[i].firstPoint+1)/LEAF_RATIO);
			if (j<=lastLeaf) {
				windowRight_start=leaves[j].firstPoint-Math.min(minIntervalLength,(leaves[j-1].lastPoint-leaves[j-1].firstPoint+1)/LEAF_RATIO);
				windowRight_end=leaves[j].firstPoint+Math.min(minIntervalLength,(leaves[j].lastPoint-leaves[j].firstPoint+1)/LEAF_RATIO);
			}
			else {
				windowRight_start=-1;
				windowRight_end=-1;
			}
			if ( ( i==0 ||   // Local maximum test
			       leaves[i].value>leaves[i-1].value ||
				   Histograms.findOutlier(x,windowLeft_start,windowLeft_end,OUTLIER_RATIO,true)
			     ) &&
			     ( j>lastLeaf ||
			       lastValue>leaves[j].value ||
				   Histograms.findOutlier(x,windowRight_start,windowRight_end,OUTLIER_RATIO,true)
				 )
			   ) {
			   	out++;
				for (k=i; k<j; k++) leaves[k].isLocalMaximum=true;
			}
			else if ( ( i==0 ||   // Local minimum test
					    leaves[i].value<leaves[i-1].value ||
					    Histograms.findOutlier(x,windowLeft_start,windowLeft_end,OUTLIER_RATIO,false)
					  ) &&
					  ( j>lastLeaf ||
					    lastValue<leaves[j].value ||
					    Histograms.findOutlier(x,windowRight_start,windowRight_end,OUTLIER_RATIO,false)
					  )
					) {
				for (k=i; k<j; k++) leaves[k].isLocalMinimum=true;
			}
			i=j;
		}
		return out;
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
	 */
	public static final void markRunsOfLocalMaximumLeaves() {
		boolean isHigh;
		int i, firstLeaf;

		if (lastLeaf>0) Arrays.sort(leaves,0,lastLeaf+1);
		for (i=0; i<lastLeaf; i++) leaves[i].marked=false;
		firstLeaf=0; nRuns=0; isHigh=leaves[0].isHigh;
		for (i=1; i<=lastLeaf; i++) {
			if (!leaves[i].isLocalMaximum) break;
			if ( leaves[i].isHigh==isHigh &&
				 leaves[i].lastPoint-leaves[firstLeaf].firstPoint<minLocalMaxDistance
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
	 * Assume that some positions in $x[firstPoint..lastPoint]$ have very large value.
	 * The tree could assign such a position to a long leaf with overall smaller
	 * value, and the value of such leaf could be smaller than the value of other
	 * leaves (which e.g. happen to be shorter). Thus, positions with very large value
	 * could be eventually missed, both as local maxima and as high local maxima.
	 * The procedure assigns every position whose value is at least $largestValue$ to a
	 * new leaf of length one, splitting in two the original leaf that contains it. This
	 * happens only if the leaf that originally contains the position is not already
	 * short and not already a local maximum. $largestValue$ is typically the output of 
	 * procedure $Histograms.largestValue()$.
	 *
	 * Remark: there could be runs of multiple consecutive positions with value at least
	 * $largestValue$, however such positions could belong to distinct leaves. The
	 * procedure does not handle such runs in an ad hoc way.
	 *
	 * Remark: the procedure assumes that $leaves$ are distinct and sorted by position.
	 *
	 * @return TRUE iff new points were added.
	 */
	private static final boolean addPointsWithLargestValue(double[] x, int firstPoint, int lastPoint, double largestValue) {
		final int MIN_LEAF_LENGTH = 5;
		int i, j, k, p;
		double oldLength, newLength;
		double oldSum, rightSum, pointValue;
		Leaf tmp;

		// Adding single-point leaves
		j=lastLeaf; k=lastLeaf;
		for (i=lastPoint; i>=firstPoint && k<leaves.length-1; i--) {
			pointValue=x[i];
			if (pointValue<largestValue) continue;
			while (j>=0) {
				if (i>=leaves[j].firstPoint && i<=leaves[j].lastPoint) {
					oldLength=leaves[j].lastPoint-leaves[j].firstPoint+1;
					if (oldLength<MIN_LEAF_LENGTH) {
						leaves[j].isPointWithLargestMass=true;
						break;
					}
					if (leaves[j].isLocalMaximum) break;
					oldSum=leaves[j].getMass(x);  // Slow, but reduces numerical instability compared to $leaves[j].value*oldLength$
					if (i==leaves[j].firstPoint) {
						newLength=oldLength-1;
						leaves[j].value=(oldSum-pointValue)/newLength;
						leaves[j].firstPoint++;
					}
					else if (i==leaves[j].lastPoint) {
						newLength=oldLength-1;
						leaves[j].value=(oldSum-pointValue)/newLength;
						leaves[j].lastPoint--;
					}
					else {
						newLength=i-leaves[j].firstPoint;
						rightSum=0.0;
						for (p=i+1; p<=leaves[j].lastPoint; p++) rightSum+=x[p];
						leaves[j].value=(oldSum-pointValue-rightSum)/newLength;
						newLength=leaves[j].lastPoint-i;
						k++;
						leaves[k].firstPoint=i+1;
						leaves[k].lastPoint=leaves[j].lastPoint;
						leaves[k].value=rightSum/newLength;
						leaves[j].lastPoint=i-1;
					}
					k++;
					leaves[k].firstPoint=i;
					leaves[k].lastPoint=i;
					leaves[k].value=pointValue;
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

}