package de.mpi_cbg.revant.util;

import de.mpi_cbg.revant.util.Math;


public class Histograms {
	
	private static Point[] valuePoints;
	private static int lastValuePoint;


	public static final void allocateMemory(int maxPoints) {
		valuePoints = new Point[maxPoints];
		for (int i=0; i<valuePoints.length; i++) valuePoints[i] = new Point();
	}


	/**
	 * @return the average of $histogram[start..end]$.
	 */
	public static final double getAverage(double[] histogram, int start, int end) {
		int i;
		double sum;

		// Checking boundaries
		if (start<0) start=0;
		else if (start>=histogram.length) start=histogram.length-1;
		if (end<start) end=start;
		else if (end>=histogram.length) end=histogram.length-1;

		// Computing
		sum=0.0;
		for (i=start; i<=end; i++) sum+=histogram[i];
		return sum/(end-start+1);
	}


	/**
	 * Remark: if all values in $histogram[start..end]$ are zero, the procedure returns
	 * zero.
	 *
	 * @return the center of mass of $histogram[start..end]$.
	 */
	public static final double getCenterOfMass(double[] histogram, int start, int end) {
		int i;
		double mass, centerOfMass;

		// Checking boundaries
		if (start<0) start=0;
		else if (start>=histogram.length) start=histogram.length-1;
		if (end<start) end=start;
		else if (end>=histogram.length) end=histogram.length-1;

		// Computing
		mass=0.0; centerOfMass=0.0;
		for (i=start; i<=end; i++) {
			mass+=histogram[i];
			centerOfMass+=i*histogram[i];
		}
		return mass==0.0?0.0:centerOfMass/mass;
	}


	/**
	 * Given a window $[start..start+window-1]$ in $histogram$, returns the rightmost
	 * $i \in [0..start-1]$ such that $histogram[i..i+window-1]$ has average greater than
	 * (if $sign=TRUE$) or at most (if $sign=FALSE$) $threshold$.
	 *
	 * @param length number of elements in $histogram$;
	 * @param orientation orientation in which $histogram$ has to be read;
	 * @return -1 if no window with the required average can be found; otherwise, the
	 * value of $i$ in the corresponding orientation.
	 */
	public static final int skipLeft(double[] histogram, int length, boolean orientation, int start, int window, double threshold, boolean sign) {
		int i, end;
		double sum;

		if (orientation) {
			// Checking boundaries
			if (start<0) start=0;
			if (start>=length) start=length-1;
			end=start+window-1;
			if (end>=length) end=length-1;
			sum=0.0;
			for (i=start; i<=end; i++) sum+=histogram[i];
			while (start>length-window) {
				start--;
				sum+=histogram[start];
			}

			// Computing
			for (i=start-1; i>=0; i--) {
				sum+=histogram[i];
				sum-=histogram[i+window];
				if ( (sign && sum>threshold*window) ||
					 ((!sign) && sum<=threshold*window) ) return i;
			}
		}
		else {
			// Checking boundaries
			start=length-start-1;
			if (start<0) start=0;
			if (start>=length) start=length-1;
			end=start-window+1;
			if (end<0) end=0;
			sum=0.0;
			for (i=start; i>=end; i--) sum+=histogram[i];
			while (start<window-1) {
				start++;
				sum+=histogram[start];
			}

			// Computing
			for (i=start+1; i<length; i++) {
				sum+=histogram[i];
				sum-=histogram[i-window];
				if ( (sign && sum>threshold*window) ||
					 ((!sign) && sum<=threshold*window) ) return length-i-1;
			}
		}
		return -1;
	}


	/**
	 * Given a window $[start..start+window-1]$ in $histogram$, returns the leftmost
	 * $i \in [start+1..length-window]$ such that $histogram[i..i+window-1]$ has average
	 * greater than (if $sign=TRUE$) or at most (if $sign=FALSE$) $threshold$.
	 *
	 * @param length number of elements in $histogram$;
	 * @param orientation orientation in which $histogram$ must be read;
	 * @return -1 if no window with the required average can be found; otherwise, the
	 * value of $i$ in the corresponding orientation.
	 */
	public static final int skipRight(double[] histogram, int length, boolean orientation, int start, int window, double threshold, boolean sign) {
		int i, end;
		double sum;

		if (orientation) {
			// Checking boundaries
			if (start<0) start=0;
			else if (start>=length) start=length-1;
			end=start+window-1;
			if (end>=length) end=length-1;
			sum=0.0;
			for (i=start; i<=end; i++) sum+=histogram[i];

			// Computing
			for (i=start+window; i<length; i++) {
				sum+=histogram[i];
				sum-=histogram[i-window];
				if ( (sign && sum>threshold*window) ||
					 ((!sign) && sum<=threshold*window) ) return i-window+1;
			}
		}
		else {
			// Checking boundaries
			start=length-start-1;
			if (start<0) start=0;
			else if (start>=length) start=length-1;
			end=start-window+1;
			if (end<0) end=0;
			sum=0.0;
			for (i=start; i>=end; i--) sum+=histogram[i];

			// Computing
			for (i=start-window; i>=0; i--) {
				sum+=histogram[i];
				sum-=histogram[i+window];
				if ( (sign && sum>threshold*window) ||
					 ((!sign) && sum<=threshold*window) ) return length-(i+window-1)-1;
			}
		}
		return -1;
	}


	/**
	 * @return the fraction of $histogram[start..end]$ that has value at least (if 
	 * $sign=true$) or smaller than (if $sign=false$) $max$.
	 */
	public static final double getFraction(double[] histogram, int start, int end, double max, boolean sign) {
		int i;
		double mass;

		// Checking boundaries
		if (start<0) start=0;
		else if (start>=histogram.length) start=histogram.length-1;
		if (end<start) end=start;
		else if (end>=histogram.length) end=histogram.length-1;

		// Computing
		mass=0.0;
		if (sign) {
			for (i=start; i<=end; i++) {
				if (histogram[i]>=max) mass+=1.0;
			}
		}
		else {
			for (i=start; i<=end; i++) {
				if (histogram[i]<max) mass+=1.0;
			}
		}
		return mass/(end-start+1);
	}


	/**
	 * Let $[start..X]$ (respectively, $[Y..end]$) be the longest prefix (respectively,
	 * suffix) of $[start..end]$ such that all its substrings of length $window$
	 * have average less than $average$. The procedure returns one if
	 * $histogram[X+1..Y-1]$ contains a substring of length $window$ whose average is less
	 * than $average$. It returns -1 if $Y-X-1$ is smaller than $window$. It returns zero
	 * otherwise.
	 *
	 * Remark: the prefix/suffix are excluded, since this function is called only by
	 * $Reads.isRandomInsertion$, and the query interval of that function might include a
	 * prefix/suffix that belongs to high-quality substrings of the read.
	 *
	 * @param length number of elements in $histogram$.
	 */
	public static final int hasSubstringWithAverage(double[] histogram, int length, int start, int end, boolean orientation, double average, int window) {
		int i, first, last, prefixLast, suffixFirst;
		int newStart, newEnd;
		double sum;

		// Checking boundaries
		if (orientation) {
			newStart=start;
			newEnd=end;
		}
		else {
			newStart=length-end-1;
			newEnd=length-start-1;
		}
		if (newStart<0) newStart=0;
		else if (newStart>=length) newStart=length-1;
		if (newEnd<newStart) newEnd=newStart;
		else if (newEnd>=length) newEnd=length-1;

		// Finding the longest prefix/suffix such that all its substrings of length
		// $substringLength$ have average less than $average$
		sum=0.0;
		last=newStart+window-1;
		if (last>newEnd) return -1;
		for (i=newStart; i<=last; i++) sum+=histogram[i];
		if (sum>=average*window) prefixLast=newStart-1;
		else {
			for (i=newStart+window; i<=newEnd; i++) {
				sum+=histogram[i]-histogram[i-window];
				if (sum>=average*window) break;
			}
			prefixLast=i-1;
		}
		sum=0.0;
		first=newEnd-window+1;
		if (first<newStart) return -1;
		for (i=first; i<=newEnd; i++) sum+=histogram[i];
		if (sum>=average*window) suffixFirst=newEnd+1;
		else {
			for (i=newEnd-window; i>=newStart; i--) {
				sum+=histogram[i]-histogram[i+window];
				if (sum>=average*window) break;
			}
			suffixFirst=i+1;
		}
		if (suffixFirst-prefixLast-1<window) return -1;

		// Searching for windows
		sum=0.0;
		last=prefixLast+window;
		for (i=prefixLast+1; i<=last; i++) sum+=histogram[i];
		if (sum<average*window) return 1;
		for (i=prefixLast+window+1; i<=suffixFirst-1; i++) {
			sum+=histogram[i]-histogram[i-window];
			if (sum<average*window) return 1;
		}
		return 0;
	}
	
	
	/**
	 * Simplified version of the previous procedure: it just searches for a substring with
	 * average at least (if $sign=TRUE$) or less than (if $sign=FALSE$) $average$.
	 *
	 * @param prefixSuffix TRUE=the substring can be a prefix or suffix of $[start..end]$; 
	 * FALSE=the substring cannot be a prefix or suffix;
	 * @param prefixSuffixDistance max distance from a boundary; used only if 
	 * $prefixSuffix=false$.
	 */
	public static final boolean hasSubstringWithAverage(double[] histogram, int start, int end, double average, boolean sign, int window, boolean prefixSuffix, int prefixSuffixDistance) {
		final double THRESHOLD = average*Math.min(window,end-start+1);
		int i;
		double sum;
		if (end-start+1<window) return false;
		
		sum=0.0;
		for (i=start; i<Math.min(start+window,histogram.length); i++) sum+=histogram[i];
		if (((sign && sum>=THRESHOLD) || (!sign && sum<THRESHOLD)) && prefixSuffix) return true;
		for (i=start+window; i<=end; i++) {
			sum+=histogram[i]-histogram[i-window];
			if ( ((sign && sum>=THRESHOLD) || (!sign && sum<THRESHOLD)) && 
				 (prefixSuffix || (i<end-prefixSuffixDistance && i-window+1>start+prefixSuffixDistance))
			   ) return true;
		}
		return false;
	}


	/**
	 * @return the largest $i \leq end$ such that, for all $j \in [start..i]$,
	 * $histogram[min\{j-window+1,start\}..j]$ has average at least (if $sign=true$) or
	 * smaller than (if $sign=false$) $threshold$. The procedure returns $start-1$ if no
	 * such $i$ can be found.
	 */
	public static final int maximalPrefix(double[] histogram, int start, int end, boolean sign, double threshold, int window) {
		int i;
		double sum;

		// Checking boundaries
		if (start<0) start=0;
		else if (start>=histogram.length) start=histogram.length-1;
		if (end<start) end=start;
		else if (end>=histogram.length) end=histogram.length-1;

		// Computing
		sum=0.0;
		for (i=start; i<=end && i<=start+window-1; i++) {
			sum+=histogram[i];
			if ( (sign && sum<threshold*(i-start+1)) ||
				 (!sign && sum>=threshold*(i-start+1)) ) return i-1;
		}
		while (i<=end) {
			sum+=histogram[i]-histogram[i-window];
			if ( (sign && sum<threshold*window) ||
			     (!sign && sum>=threshold*window) ) return i-1;
			i++;
		}
		return end;
	}


	/**
	 * @return the smallest $i \geq start$ such that, for all $j \in [i..end]$,
	 * $histogram[j..\min\{end,j+window-1\}]$ has average at least (if $sign=true$) or
	 * smaller than (if $sign=false$) $threshold$. The procedure returns $end+1$ if no
	 * such $i$ can be found.
	 */
	public static final int maximalSuffix(double[] histogram, int start, int end, boolean sign, double threshold, int window) {
		int i;
		double sum;

		// Checking boundaries
		if (start<0) start=0;
		else if (start>=histogram.length) start=histogram.length-1;
		if (end<start) end=start;
		else if (end>=histogram.length) end=histogram.length-1;

		// Computing
		sum=0.0;
		for (i=end; i>=start && i>=end-window+1; i--) {
			sum+=histogram[i];
			if ( (sign && sum<threshold*(end-i+1)) ||
			     (!sign && sum>=threshold*(end-i+1)) ) return i+1;
		}
		while (i>=start) {
			sum+=histogram[i]-histogram[i+window];
			if ( (sign && sum<threshold*window) ||
			     (!sign && sum>=threshold*window) ) return i+1;
			i--;
		}
		return start;
	}


	/**
	 * @param q prints the average of $histogram[start..end]$ every $q$ positions,
	 * adding to each position value $offset$.
	 */
	public static final void printHistogram(double[] histogram, int start, int end, int q, int offset) {
		int i, j;

		// Checking boundaries
		if (start<0) start=0;
		else if (start>=histogram.length) start=histogram.length-1;
		if (end<start) end=start;
		else if (end>=histogram.length) end=histogram.length-1;

		// Computing
		for (i=start; i<=end; i+=q) {
			j=i+q-1;
			if (j>end) j=end;
			if (IO.SHOW_STD_OUT) IO.printOut((offset+i)+": "+getAverage(histogram,i,j)+"");
		}
	}


	/**
	 * Estimates the peaks of a piecewise-linear $histogram[first..last]$ as follows: (1)
	 * computes a regression tree on $histogram$ and finds intervals whose average is a
	 * local maximum with respect to their adjacent intervals; (2) marks some local maxima
	 * as "high", using the distribution of the values of all local maxima; (3) computes
	 * runs of local-maximum intervals that are too close to each other; (4) assigns a
	 * peak to the center of mass of each run that contains a high local maximum.
	 *
	 * Array $peaks[pos..x]$ contains the sorted sequence $f_0,p_0,l_0,f_1,p_1,l_1,...$,
	 * where $p_i$ is the position of a peak plus offset $start$, $f_i$ (respectively,
	 * $l_i$) is the first (respectively, last) position of the rightmost (respectively,
	 * leftmost) leaf that is not a local maximum and that occurs before (respectively,
	 * after) the run of local-maximum leaves associated with $p_i$, plus offset $start$.
	 * Value $x$ is returned in output. The procedure tries to combine the peaks in
	 * $pos,pos+1,pos+2$ and $pos-3,pos-2,pos-1$ if they are approx. identical.
	 *
	 * If $isDelta=TRUE$, the procedure behaves as follows:
	 *
	 * First, the procedure removes the high mark from all high leaves with low density.
	 *
	 * Then, if a run contains exactly two local maxima, the cluster center is either 
	 * assigned to the midpoint of the leftmost local-minimum leaf between such local 
	 * maxima, or to the midpoint between the last position of the first local maximum and
	 * the first position of the second local maximum, if no local-minimum leaf is found 
	 * between the local maxima.
	 *
	 * The special treatment for runs with exactly two local maxima is designed for
	 * function $\delta$, and it is motivated by the fact that two adjacent repeats could
	 * produce two, very close peaks at their interface, separated by a very steep local
	 * minimum (e.g. if the two repeats never occur together elsewhere). Assigning the
	 * peak based on the center of mass would move the peak closer to the repeat with
	 * higher value of $\delta$, which is not necessarily correct.
	 *
	 * @param maxStd intervals of the regression tree with standard deviation at most this
	 * much are not split;
	 * @param distanceThreshold* all runs, except possibly the first and the last, induce 
	 * a peak at their center of mass. The first run induces a peak at its first position 
	 * if such first position is at distance at most $distanceThresholdFirst$ from 
	 * $first$; the last run induces a peak at its last position if such last position is 
	 * at distance at most $distanceThresholdLast$ from $last$.
	 * @param minHigh enforces that all high-enough local maxima are marked as high, even 
	 * if $setHighLocalMax()$ doesn't mark them; non-positive values disable such 
	 * behavior;
	 * @param tmp temporary array, of length at least equal to $peaks$;
	 * @param period used only if $isDelta=true$ and $period>-1$; if: (1) the set of all 
	 * local maxima covers a large fraction of $[first..last]$, or (2) $period$ is large,
	 * and most local maxima have another local maximum at distance equal to $period$, the 
	 * procedure marks as high only the local maxima that are close to $first$ or $last$; 
	 * this is because short-period intervals (condition 1) and long-period intervals
	 * (condition 2) can have a large number of strong, uniformly spaced peaks, and only 
	 * some of them might end up being marked as high by the normal procedure;
	 * @param maxWidth maximum width of a local maximum to be considered high;
	 * @param points if not null, it is used to compute the total number of events inside 
	 * a local maximum: if such number if at least $minMassHigh$, the local maximum is
	 * marked as high, even though it would not have been marked as high otherwise;
	 * @return $pos-1$ if $histogram$ is never split by the regression tree.
	 */
	public static final int getPeaks(double[] histogram, int first, int last, int[] peaks, int pos, int start, double maxStd, int distanceThresholdFirst, int distanceThresholdLast, boolean isDelta, double constantThreshold, double minHigh, int maxPeakRadius, int[] tmp, int period, int maxWidth, Point[] points, int lastPoint, int minMassHigh, int minAlignmentLength) {
		double MIN_PEAK_DENSITY_DELTA = 0;  // Set to a quantile of the density of all leaves of the regression tree on $histogram[first..last]$. A constant would not be correct, since leaf densities depend e.g. on the number of alignments that cover the leaf.
		double MIN_PEAK_DENSITY_DELTA_QUANTILE = 0.5;  // Arbitrary
		double MIN_PEAK_DENSITY_DELTA_RATIO = 50.0;  // Arbitrary
		final double SURFACE_THRESHOLD = 0.75;  // Arbitrary
		final double PEAKS_RATIO = 0.5;  // Arbitrary
		final int MIN_LOCAL_MAXIMA = 4;  // Arbitrary
		final int MIN_PEAK_DENSITY_WINDOW = IO.quantum>>1;  // Arbitrary
		final int MIN_LONG_PERIOD = minAlignmentLength>>1;
		boolean isHigh, defaultHighMarking, justFirstLocalMax;
		int i, j, peak, next, lastLeaf, firstLocalMaximum;
		int firstLocalMinimum, firstNormal, nLocalMaximumLeaves, center, boundary;
		final int initialPos = pos;
		double largestValue, surface, nonzeroSurface;
		Leaf[] leaves;

		if (isConstant(histogram,first,last,constantThreshold)) return pos-1;
		largestValue=(isDelta&&period>-1)?0:largestValue(histogram,first,last,-1);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Histogram.getPeaks> largestValue="+largestValue+" start="+start);

		nLocalMaximumLeaves=RegressionTree.buildRegressionTree(histogram,first,last,maxStd,1,largestValue,true,minAlignmentLength);
		if (RegressionTree.lastLeaf==0 || nLocalMaximumLeaves==0) return pos-1;
		justFirstLocalMax=oneMaxAndConstant(histogram,first,last,constantThreshold,maxWidth)||allMaxAndConstant(histogram,first,last,constantThreshold,maxWidth);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Histogram.getPeaks()> oneMaxAndConstant="+justFirstLocalMax);
		if (!justFirstLocalMax) {
			if (isDelta) {
				lastLeaf=-1;
				for (i=0; i<=RegressionTree.lastLeaf; i++) {
					// We do not use leaves of size one, since they represent single
					// points.
					if (RegressionTree.leaves[i].lastPoint==RegressionTree.leaves[i].firstPoint) continue;
					lastLeaf++;
					valuePoints[lastLeaf].position=RegressionTree.leaves[i].value;
					valuePoints[lastLeaf].mass=1;
				}
				if (lastLeaf<=0) MIN_PEAK_DENSITY_DELTA=constantThreshold;
				else {
					lastLeaf=Points.sortAndCompact(valuePoints,lastLeaf);
					if (lastLeaf<=0) MIN_PEAK_DENSITY_DELTA=constantThreshold;
					else {
						MIN_PEAK_DENSITY_DELTA=-1;
						for (i=lastLeaf; i>0; i--) {
							if (valuePoints[i-1].position>Leaves.ZERO && valuePoints[i].position>=valuePoints[i-1].position*MIN_PEAK_DENSITY_DELTA_RATIO) {
								MIN_PEAK_DENSITY_DELTA=valuePoints[i].position;
								break;
							}
						}
						if (MIN_PEAK_DENSITY_DELTA==-1) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Histogram.getPeaks> computing quantile");
							MIN_PEAK_DENSITY_DELTA=Points.quantile(valuePoints,0,lastLeaf,true,MIN_PEAK_DENSITY_DELTA_QUANTILE,false);
						}
						MIN_PEAK_DENSITY_DELTA=Math.max(MIN_PEAK_DENSITY_DELTA,constantThreshold);
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr("Histogram.getPeaks> MIN_PEAK_DENSITY_DELTA="+MIN_PEAK_DENSITY_DELTA);
	IO.printErr("Histogram.getPeaks> values: ");
	for (int x=0; x<=lastLeaf; x++) IO.printErr(valuePoints[x].position+","+valuePoints[x].mass);
	IO.printErr("");
}
					}
				}
			}
			defaultHighMarking=true;
			if (isDelta && period>-1) {
				// Special case: delta histogram of a periodic region.
				surface=Leaves.getLocalMaxSurface(RegressionTree.leaves,RegressionTree.lastLeaf,null,false);
				nonzeroSurface=Leaves.getNonzeroSurface(RegressionTree.leaves,RegressionTree.lastLeaf,null,false);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Histogram.getPeaks> surface ratio = "+(surface/nonzeroSurface)+" PERIOD="+period+" for range ["+first+".."+last+"] peaksSpacedByPeriod="+Leaves.peaksSpacedByPeriod(RegressionTree.leaves,RegressionTree.lastLeaf,histogram,period,MIN_LOCAL_MAXIMA,Math.max(distanceThresholdFirst,distanceThresholdLast),PEAKS_RATIO,MIN_PEAK_DENSITY_DELTA,maxWidth)+" distance="+Math.max(distanceThresholdFirst,distanceThresholdLast));
				if (surface>=nonzeroSurface*SURFACE_THRESHOLD || (period>=MIN_LONG_PERIOD && Leaves.peaksSpacedByPeriod(RegressionTree.leaves,RegressionTree.lastLeaf,histogram,period,MIN_LOCAL_MAXIMA,Math.max(distanceThresholdFirst,distanceThresholdLast),PEAKS_RATIO,MIN_PEAK_DENSITY_DELTA,maxWidth))) {
					for (i=0; i<=RegressionTree.lastLeaf; i++) {
						if (!RegressionTree.leaves[i].isLocalMaximum || RegressionTree.leaves[i].lastPoint-RegressionTree.leaves[i].firstPoint+1>maxWidth) continue;
						RegressionTree.leaves[i].isHigh=RegressionTree.leaves[i].lastPoint<=first+distanceThresholdFirst||RegressionTree.leaves[i].firstPoint>=last-distanceThresholdLast;
					}
					defaultHighMarking=false;
				}
			}
			if (defaultHighMarking) {
				// Not a delta histogram, or delta histogram of a non-periodic region, or
				// delta histogram of a periodic region with not enough local max surface.
				Leaves.setHighLocalMax(RegressionTree.leaves,RegressionTree.lastLeaf,last,true,true,null,0,false,maxWidth);
				if (minHigh>0) {
					for (i=0; i<=RegressionTree.lastLeaf; i++) {
						if ( RegressionTree.leaves[i].isLocalMaximum && RegressionTree.leaves[i].lastPoint-RegressionTree.leaves[i].firstPoint+1<=maxWidth && 
							 (RegressionTree.leaves[i].value>=minHigh || (points!=null && Points.mass(points,lastPoint,RegressionTree.leaves[i].firstPoint,RegressionTree.leaves[i].lastPoint)>=minMassHigh))
						   ) RegressionTree.leaves[i].isHigh=true;
					}
				}
				if (isDelta) {
					for (i=0; i<=RegressionTree.lastLeaf; i++) {
						if ( RegressionTree.leaves[i].isLocalMaximum && RegressionTree.leaves[i].isHigh && 
							 (points==null || Points.mass(points,lastPoint,RegressionTree.leaves[i].firstPoint,RegressionTree.leaves[i].lastPoint)<minMassHigh) &&
							 RegressionTree.leaves[i].value<MIN_PEAK_DENSITY_DELTA && !hasSubstringWithAverage(histogram,RegressionTree.leaves[i].firstPoint,RegressionTree.leaves[i].lastPoint,MIN_PEAK_DENSITY_DELTA,true,MIN_PEAK_DENSITY_WINDOW,true,-1)
						   ) RegressionTree.leaves[i].isHigh=false;
					}
				}
			}
		}
		RegressionTree.markRunsOfLocalMaximumLeaves();
		leaves=RegressionTree.leaves;
		lastLeaf=RegressionTree.lastLeaf;
if (IO.SHOW_STD_ERR_PRIME) {
	for (i=0; i<=lastLeaf; i++) {
		IO.printErr("Histogram.getPeaks> ["+leaves[i].firstPoint+".."+leaves[i].lastPoint+"] isLocalMaximum="+leaves[i].isLocalMaximum+" isLocalMinimum="+leaves[i].isLocalMinimum+" marked="+leaves[i].marked+" value="+leaves[i].value+" high="+leaves[i].isHigh);
	}
}


		// Computing peaks
		firstLocalMinimum=RegressionTree.lastLocalMaximum+1;
		firstNormal=RegressionTree.lastLocalMinimum+1;
		firstLocalMaximum=0; isHigh=false;
		for (i=0; i<=RegressionTree.lastLocalMaximum; i++) {
			if (leaves[i].isHigh) isHigh=true;
			if (!leaves[i].marked) continue;
			if (!isHigh) {
				firstLocalMaximum=i+1;
				isHigh=false;
				continue;
			}
			peak=start;
			if (firstLocalMaximum==0 && leaves[firstLocalMaximum].lastPoint-first<=distanceThresholdFirst) peak+=leaves[firstLocalMaximum].firstPoint;
			else if ( (i==lastLeaf||!leaves[i+1].isLocalMaximum) && (last-leaves[i].firstPoint<=distanceThresholdLast) ) peak+=leaves[i].lastPoint;
			else {
				if (i==firstLocalMaximum+1 && isDelta) {
					next=-1;
					if (firstLocalMinimum<=lastLeaf) next=Leaves.includedLocalMinimum(leaves,lastLeaf,leaves[firstLocalMaximum].lastPoint+1,leaves[i].firstPoint-1,firstLocalMinimum,firstNormal);
					if (next>=0) {
						center=(int)Math.ceil(getCenterOfMass(histogram,leaves[next].firstPoint,leaves[next].lastPoint));
						if (center<=leaves[next].firstPoint || center>=leaves[next].lastPoint) center=(leaves[next].firstPoint+leaves[next].lastPoint)>>1;
						peak+=center;
					}
					else peak+=(leaves[firstLocalMaximum].lastPoint+leaves[i].firstPoint)>>1;
				}
				else {
					center=(int)Math.ceil(getCenterOfMass(histogram,leaves[firstLocalMaximum].firstPoint,leaves[i].lastPoint));
					if (center<=leaves[firstLocalMaximum].firstPoint || center>=leaves[i].lastPoint) center=(leaves[firstLocalMaximum].firstPoint+leaves[i].lastPoint)>>1;
					peak+=center;
				}
			}
			next=-1;
			boundary=start+(firstLocalMaximum==0?first-1:leaves[firstLocalMaximum-1].lastPoint);
			if (boundary<peak-maxPeakRadius) boundary=peak-maxPeakRadius;
			if (leaves[firstLocalMaximum].firstPoint>first) next=Leaves.nextNonMaximum(leaves,lastLeaf,firstLocalMaximum,false,boundary-start,firstLocalMinimum,firstNormal,true);
			tmp[pos]=next;
			if (next==-1) {
				peaks[pos]=start+leaves[firstLocalMaximum].firstPoint;
				if (peaks[pos]<peak-maxPeakRadius) peaks[pos]=peak-maxPeakRadius;
			}
			else {
				if (pos>0 && tmp[pos-1]==next) {
					center=(int)getCenterOfMass(histogram,leaves[next].firstPoint,leaves[next].lastPoint);
					if (center<=leaves[next].firstPoint || center>=leaves[next].lastPoint) center=(leaves[next].firstPoint+leaves[next].lastPoint)>>1;
					peaks[pos]=start+center+1;
					peaks[pos-1]=start+center;
					if (peaks[pos]<peak-maxPeakRadius) peaks[pos]=peak-maxPeakRadius;
					if (peaks[pos-1]>peaks[pos-2]+maxPeakRadius) peaks[pos-1]=peaks[pos-2]+maxPeakRadius;
				}
				else {
					if (leaves[next].isLocalMinimum) {
						center=(int)getCenterOfMass(histogram,leaves[next].firstPoint,leaves[next].lastPoint);
						peaks[pos]=start+(center==0?(leaves[next].firstPoint+leaves[next].lastPoint)>>1:center);
					}
					else peaks[pos]=start+leaves[next].firstPoint;
					if (peaks[pos]<peak-maxPeakRadius) peaks[pos]=peak-maxPeakRadius;
				}
			}
			pos++;
			peaks[pos]=peak;
			pos++;
			next=-1;
			boundary=start+((i==firstLocalMinimum-1||i==firstNormal-1||i==lastLeaf)?last+1:leaves[i+1].firstPoint);
			if (boundary>peak+maxPeakRadius) boundary=peak+maxPeakRadius;
			if (leaves[i].lastPoint<last) next=Leaves.nextNonMaximum(leaves,lastLeaf,i,true,boundary-start,firstLocalMinimum,firstNormal,true);
			tmp[pos]=next;
			if (next==-1) peaks[pos]=start+leaves[i].lastPoint;
			else {
				if (leaves[next].isLocalMinimum) {
					center=(int)getCenterOfMass(histogram,leaves[next].firstPoint,leaves[next].lastPoint);
					peaks[pos]=start+(center==0?(leaves[next].firstPoint+leaves[next].lastPoint)>>1:center);
				}
				else peaks[pos]=start+leaves[next].lastPoint;
			}
			if (peaks[pos]>peak+maxPeakRadius) peaks[pos]=peak+maxPeakRadius;
			pos++;
			firstLocalMaximum=i+1;
			isHigh=false;
if (IO.SHOW_STD_ERR_PRIME) {
	IO.printErr(">> assigned peak "+peak+" >> peaks: "+peaks[pos-3]+","+peaks[pos-2]+","+peaks[pos-1]);
	IO.printErr(">> peaks positions: "+(pos-3)+","+(pos-2)+","+(pos-1));
	IO.printErr(">> values in tmp: "+tmp[pos-3]+",--,"+tmp[pos-1]);
	System.err.println(":::: LEAVES OF THE REGRESSION TREE:");
	for (int x=0; x<=lastLeaf; x++) System.err.println(leaves[x]);
}
		}
		
		// Merging the first new peak to the last old peak if they are close enough.
		if (initialPos>0 && pos>initialPos && Math.abs(peaks[initialPos+1],peaks[initialPos-2])<=Math.min(distanceThresholdFirst,distanceThresholdLast)) {
			boundary=Math.min(peaks[initialPos-3],peaks[initialPos]);
			peaks[initialPos-3]=initialPos-3>0?Math.max(boundary,peaks[initialPos-4]+1):boundary;
			peaks[initialPos-2]=(peaks[initialPos-2]+peaks[initialPos+1])>>1;
			boundary=Math.max(peaks[initialPos-1],peaks[initialPos+2]);
			peaks[initialPos-1]=initialPos+3<pos?Math.min(boundary,peaks[initialPos+3]-1):boundary;
			for (i=initialPos; i<pos-3; i++) peaks[i]=peaks[i+3];
			pos-=3;
		}
		
		return pos-1;
	}

	
	/**
	 * If $RegressionTree$ contains a run of local-maximum leaves such that the histogram 
	 * outside the run is constant, the procedure marks such run as high, and it marks all
	 * other local-maximum leaves as not high.
	 *
	 * Remark: the run can be anywhere in the histogram. The closest local minima on each
	 * side are used to delimit the regions to be tested by 
	 * $isConstant(constantThreshold)$.
	 *
	 * Remark: the procedure assumes the leaves in $RegressionTree$ to be sorted by 
	 * position.
	 * 
	 * @param histogram[first..last] the histogram on which the regression tree was built;
	 * @param maxWidth maximum width of a run of local-maximum leaves to be considered by
	 * the procedure;
	 * @return TRUE iff $RegressionTree$ contains such a run of local-maximum leaves.
	 */
	public static final boolean oneMaxAndConstant(double[] histogram, int first, int last, double constantThreshold, int maxWidth) {
		final double MAX_LEFTRIGHT_RATIO = 2.0;  // Arbitrary
		boolean isConstantLeft, isConstantRight;
		int i, j;
		int min, minLeft, minRight, normal, to, previousLocalMax, from;
		double avgLeft, avgRight;
		
		i=0; previousLocalMax=-1;
		while (i<=RegressionTree.lastLeaf) {
			if (!RegressionTree.leaves[i].isLocalMaximum) {
				i++;
				continue;
			}
			if (RegressionTree.leaves[i].lastPoint-RegressionTree.leaves[i].firstPoint+1>maxWidth) {
				for (j=i+1; j<=RegressionTree.lastLeaf; j++) {
					if (!RegressionTree.leaves[j].isLocalMaximum) break;
					else i=j;
				}
				previousLocalMax=i;
				i++;
				continue;
			}
			// Left
			min=-1; normal=-1; to=Math.max(0,previousLocalMax+1);
			for (j=i-1; j>=to; j--) {
				if (RegressionTree.leaves[j].isLocalMinimum) {
					min=j;
					break;
				}
				else if (normal==-1) normal=j;
			}
			minLeft=min>=0?min:(normal>=0?normal:i);
			// Right
			from=RegressionTree.leaves[i].firstPoint;
			for (j=i+1; j<=RegressionTree.lastLeaf; j++) {
				if (!RegressionTree.leaves[j].isLocalMaximum) break;
				else i=j;
			}
			if (RegressionTree.leaves[i].lastPoint-from+1>maxWidth || getAverage(histogram,from,RegressionTree.leaves[i].lastPoint)<=constantThreshold) {
				previousLocalMax=i;
				i++;
				continue;
			}
			min=-1; normal=-1;
			for (j=i+1; j<=RegressionTree.lastLeaf; j++) {
				if (RegressionTree.leaves[j].isLocalMaximum) break;
				if (RegressionTree.leaves[j].isLocalMinimum) {
					min=j;
					break;
				}
				else if (normal==-1) normal=j;
			}
			minRight=min>=0?min:(normal>=0?normal:i);
			// Constancy test
			avgLeft=RegressionTree.leaves[minLeft].firstPoint==first?0:getAverage(histogram,first,RegressionTree.leaves[minLeft].firstPoint-1);
			avgRight=RegressionTree.leaves[minRight].lastPoint==last?0:getAverage(histogram,RegressionTree.leaves[minRight].lastPoint+1,last);
			isConstantLeft=RegressionTree.leaves[minLeft].firstPoint==first?true:isConstant(histogram,first,RegressionTree.leaves[minLeft].firstPoint-1,constantThreshold);
			isConstantRight=RegressionTree.leaves[minRight].lastPoint==last?true:isConstant(histogram,RegressionTree.leaves[minRight].lastPoint+1,last,constantThreshold);
			if ( isConstantLeft && isConstantRight &&
				 ( (avgLeft>Leaves.ZERO && avgRight>Leaves.ZERO && Math.max(avgLeft,avgRight)<=Math.min(avgLeft,avgRight)*MAX_LEFTRIGHT_RATIO) ||
				   (avgLeft<=Leaves.ZERO && avgRight<=constantThreshold) ||
				   (avgRight<=Leaves.ZERO && avgLeft<=constantThreshold)
				 )
			   ) {
				j=i;
				while (j>=0 && RegressionTree.leaves[j].isLocalMaximum) {
					RegressionTree.leaves[j].isHigh=true;
					j--;
				}
				while (j>=0) {
					RegressionTree.leaves[j].isHigh=false;
					j--;
				}
				for (j=i+1; j<=RegressionTree.lastLeaf; j++) RegressionTree.leaves[j].isHigh=false;
				return true;
			}
			previousLocalMax=i;
			i++;
		}
		return false;
	}
	
	
	/**
	 * Same as $oneMaxAndConstant()$, but checks the case in which the histogram is 
	 * constant everywhere, except in all runs of local-maximum leaves of total width at 
	 * most $maxWidth$ each.
	 *
	 * @return TRUE iff $RegressionTree$ contains such runs of local-maximum leaves.
	 */
	public static final boolean allMaxAndConstant(double[] histogram, int first, int last, double constantThreshold, int maxWidth) {
		final double MAX_LEFTRIGHT_RATIO = 2.0;  // Arbitrary
		boolean isConstantLeft, isConstantRight;
		int i, j;
		int min, minLeft, minRight, normal, to, previousLocalMax, previousRangeStart, lastLocalMax;
		double avgLeft, previousAvgLeft, avgRight, minAvg, maxAvg;
		
		// Checking constancy
		i=0; previousLocalMax=-1; previousRangeStart=first; previousAvgLeft=Leaves.ZERO;
		minAvg=Math.POSITIVE_INFINITY; maxAvg=0.0;
		while (i<=RegressionTree.lastLeaf) {
			if (!RegressionTree.leaves[i].isLocalMaximum) {
				i++;
				continue;
			}
			if (RegressionTree.leaves[i].lastPoint-RegressionTree.leaves[i].firstPoint+1>maxWidth) {
				for (j=i+1; j<=RegressionTree.lastLeaf; j++) {
					if (!RegressionTree.leaves[j].isLocalMaximum) break;
					else i=j;
				}
				previousLocalMax=i;
				i++;
				continue;
			}
			// Right
			lastLocalMax=i;
			for (j=i+1; j<=RegressionTree.lastLeaf; j++) {
				if (!RegressionTree.leaves[j].isLocalMaximum) break;
				else lastLocalMax=j;
			}
			if (RegressionTree.leaves[lastLocalMax].lastPoint-RegressionTree.leaves[i].firstPoint+1>maxWidth || getAverage(histogram,RegressionTree.leaves[i].firstPoint,RegressionTree.leaves[lastLocalMax].lastPoint)<=constantThreshold) {
				previousLocalMax=lastLocalMax;
				i=lastLocalMax+1;
				continue;
			}
			// Left
			min=-1; normal=-1; to=Math.max(0,previousLocalMax+1);
			for (j=i-1; j>=to; j--) {
				if (RegressionTree.leaves[j].isLocalMinimum) {
					min=j;
					break;
				}
				else if (normal==-1) normal=j;
			}
			minLeft=min>=0?min:(normal>=0?normal:i);
			// Constancy test (left).
			avgLeft=RegressionTree.leaves[minLeft].firstPoint==previousRangeStart?0:getAverage(histogram,previousRangeStart,RegressionTree.leaves[minLeft].firstPoint-1);
			minAvg=Math.max(minAvg,avgLeft); maxAvg=Math.max(maxAvg,avgLeft);
			isConstantLeft=RegressionTree.leaves[minLeft].firstPoint==previousRangeStart?true:isConstant(histogram,previousRangeStart,RegressionTree.leaves[minLeft].firstPoint-1,constantThreshold);
			if ( !isConstantLeft || 
				 (previousAvgLeft>Leaves.ZERO && avgLeft>Leaves.ZERO && Math.max(previousAvgLeft,avgLeft)>Math.min(previousAvgLeft,avgLeft)*MAX_LEFTRIGHT_RATIO) ||
  				 (previousAvgLeft<=Leaves.ZERO && avgLeft>constantThreshold) ||
  				 (avgLeft<=Leaves.ZERO && previousAvgLeft>constantThreshold)
			   ) return false;
			previousAvgLeft=avgLeft;
			min=-1; normal=-1;
			for (j=lastLocalMax+1; j<=RegressionTree.lastLeaf; j++) {
				if (RegressionTree.leaves[j].isLocalMaximum) break;
				if (RegressionTree.leaves[j].isLocalMinimum) {
					min=j;
					break;
				}
				else if (normal==-1) normal=j;
			}
			minRight=min>=0?min:(normal>=0?normal:lastLocalMax);
			previousRangeStart=RegressionTree.leaves[minRight].lastPoint+1;
			previousLocalMax=lastLocalMax;
			i=lastLocalMax+1;
		}
		// Last range
		avgRight=previousRangeStart>=last?0:getAverage(histogram,previousRangeStart,last);
		minAvg=Math.max(minAvg,avgRight); maxAvg=Math.max(maxAvg,avgRight);
		isConstantRight=previousRangeStart>=last?true:isConstant(histogram,previousRangeStart,last,constantThreshold);
		if ( !isConstantRight || 
			 (previousAvgLeft>Leaves.ZERO && avgRight>Leaves.ZERO && Math.max(previousAvgLeft,avgRight)>Math.min(previousAvgLeft,avgRight)*MAX_LEFTRIGHT_RATIO) || 
			 (previousAvgLeft<=Leaves.ZERO && avgRight>constantThreshold) ||
			 (avgRight<=Leaves.ZERO && previousAvgLeft>constantThreshold) ||
			 maxAvg>minAvg*MAX_LEFTRIGHT_RATIO
		   ) return false;
		
		// Marking high maxima
		for (i=0; i<=RegressionTree.lastLeaf; i++) RegressionTree.leaves[i].isHigh=false;
		i=0;
		while (i<=RegressionTree.lastLeaf) {
			if (!RegressionTree.leaves[i].isLocalMaximum) {
				i++;
				continue;
			}
			if (RegressionTree.leaves[i].lastPoint-RegressionTree.leaves[i].firstPoint+1>maxWidth) {
				for (j=i+1; j<=RegressionTree.lastLeaf; j++) {
					if (!RegressionTree.leaves[j].isLocalMaximum) break;
					else i=j;
				}
				i++;
				continue;
			}
			lastLocalMax=i;
			for (j=i+1; j<=RegressionTree.lastLeaf; j++) {
				if (!RegressionTree.leaves[j].isLocalMaximum) break;
				else lastLocalMax=j;
			}
			if (RegressionTree.leaves[lastLocalMax].lastPoint-RegressionTree.leaves[i].firstPoint+1>maxWidth || getAverage(histogram,RegressionTree.leaves[i].firstPoint,RegressionTree.leaves[lastLocalMax].lastPoint)<=constantThreshold) {
				i=lastLocalMax+1;
				continue;
			}
			for (j=i; j<=lastLocalMax; j++) RegressionTree.leaves[j].isHigh=true;
			i=lastLocalMax+1;
		}
		return true;
	}


	/**
	 * @return the position in $[firstPosition..lastPosition]$ with maximum absolute
	 * second derivative of function $f$. The second derivative at position $x$ is
	 * estimated as the difference between the average of all position in $[x-step/2..
	 * x+step/2]$, and the average of all position in $(x+step/2..x+3*step/2]$ and in
	 * $[x-3*step/2..x-step/2)$. Positions outside $f$ are assumed to have value zero.
	 *
	 * @param fLength number of elements in array $f$;
	 * @param direction returns the leftmost (-1) or the rightmost (1) of all positions
	 * with maximum second derivative. If $direction=0$, the procedure returns the center
	 * of mass of all positions with maximum second derivative.
	 */
	public static final int getMaxSecondDerivative(int firstPosition, int lastPosition, double[] f, int fLength, int step, int direction) {
		int i, j, span, first, last, leftmostMax, rightmostMax, nMax, centerOfMass;
		double value, valueLeft, valueRight;
		double derivative, maxDerivative;

		// Checking boundaries
		if (firstPosition<0) firstPosition=0;
		else if (firstPosition>=fLength) firstPosition=fLength-1;
		if (lastPosition<firstPosition) lastPosition=firstPosition;
		else if (lastPosition>=fLength) lastPosition=fLength-1;

		// Computing
		maxDerivative=-1; leftmostMax=-1; rightmostMax=-1; nMax=0; centerOfMass=0;
		span=step+(step>>1);
		for (i=firstPosition; i<=lastPosition; i++) {
			// Estimating the value of the current position
			first=i-(step>>1);
			if (first<0) first=0;
			last=i+(step>>1);
			if (last>=fLength) last=fLength-1;
			value=0.0;
			for (j=first; j<=last; j++) value+=f[j];
			value/=1+((step>>1)<<1);

			// Estimating the value of the previous position
			first=i-span;
			if (first<0) first=0;
			last=i-(step>>1)-1;
			if (last>=fLength) last=fLength-1;
			valueLeft=0.0;
			for (j=first; j<=last; j++) valueLeft+=f[j];
			valueLeft/=step;

			// Estimating the value of the next position
			first=i+1+(step>>1);
			if (first<0) first=0;
			last=i+span;
			if (last>=fLength) last=fLength-1;
			valueRight=0.0;
			for (j=first; j<=last; j++) valueRight+=f[j];
			valueRight/=step;

			// Computing derivatives
			derivative=(valueRight-value)-(value-valueLeft);
			if (derivative<0) derivative=-derivative;
			if (derivative==maxDerivative) {
				nMax++;
				centerOfMass+=i;
				rightmostMax=i;
			}
			else if (derivative>maxDerivative) {
				maxDerivative=derivative;
				nMax=1;
				centerOfMass=i;
				leftmostMax=i;
				rightmostMax=i;
			}
		}

		// Output
		switch (direction) {
			case -1: return leftmostMax;
			case 0: return centerOfMass/nMax;
			case 1: return rightmostMax;
		}
		return -1;
	}
	
	
	/**
	 * @return true iff the difference between the largest and the smallest value in
	 * $histogram[firstPosition..lastPosition]$ is at most $threshold$.
	 */
	public static final boolean isConstant(double[] histogram, int firstPosition, int lastPosition, double threshold) {
		int i;
		double min, max;

		// Checking boundaries
		if (firstPosition<0) firstPosition=0;
		else if (firstPosition>=histogram.length) firstPosition=histogram.length-1;
		if (lastPosition<firstPosition) lastPosition=firstPosition;
		else if (lastPosition>=histogram.length) lastPosition=histogram.length-1;
		if (lastPosition-firstPosition<1) return true;

		// Computing
		min=histogram[firstPosition]; max=min;
		for (i=firstPosition+1; i<=lastPosition; i++) {
			if (histogram[i]<min) {
				min=histogram[i];
				if (max-min>threshold) return false;
			}
			if (histogram[i]>max) {
				max=histogram[i];
				if (max-min>threshold) return false;
			}
		}
		return true;
	}
	
	
	/**
	 * Uses Pearson's chi-squared test to decide whether $histogram[firstPosition..lastPosition]$, 
	 * \emph{with integer values}, is distributed uniformly on a set of equal bins of length 
	 * $binLength$ that divide the length $lastPosition-firstPosition+1$.
	 *
	 * Remark: $binLength$ is reset if the resulting number of bins is bigger than
	 * $Math.CHI_SQUARE.length+1$.
	 */
	public static final boolean isConstant(double[] histogram, int firstPosition, int lastPosition, int binLength) {
		final int MIN_LENGTH = 10;
		int i, j;
		int nBins, nDegreesOfFreedom, length;
		double mass, expectation, sum, chiSquare;
		if (lastPosition-firstPosition<1) return false;
		
		mass=0;
		for (j=firstPosition; j<=lastPosition; j++) mass+=histogram[j];
		length=lastPosition-firstPosition+1;

if (IO.SHOW_STD_ERR) IO.printErr("areUniformlyDistributed> length="+length+" binLength="+binLength);

		nBins=(int)Math.ceil(length/binLength);

if (IO.SHOW_STD_ERR) IO.printErr("areUniformlyDistributed> nBins="+nBins);

		nDegreesOfFreedom=nBins-1;
		if (nDegreesOfFreedom>Math.CHI_SQUARE.length) {
			nBins=Math.CHI_SQUARE.length;
			nDegreesOfFreedom=nBins-1;
			binLength=(int)Math.ceil(length/nBins);
		}
		expectation=mass/nBins;
		chiSquare=0.0; sum=0;
		i=binLength; j=firstPosition;
		while (j<=lastPosition) {
			if (j<i) {
				sum+=histogram[j];
				j++;
			}
			else {
				chiSquare+=(sum-expectation)*(sum-expectation)/expectation;
				i+=binLength;
				sum=0;
			}
		}
		chiSquare+=(sum-expectation)*(sum-expectation)/expectation;
if (IO.SHOW_STD_ERR) IO.printErr("areUniformlyDistributed> chiSquare="+chiSquare+" 99.9%="+Math.CHI_SQUARE[nDegreesOfFreedom-1][4]+" nBins="+nBins+" binLength="+binLength);
		return chiSquare<=Math.CHI_SQUARE[nDegreesOfFreedom-1][4];
	}


	/**
	 * Tries to detect positions with an anomalously high value. This is done by fitting a
	 * DET on the set of all nonzero values of $histogram[start..end]$, and by returning 
	 * the first value of a selected run of local-maximum leaves, if there are at least 
	 * two runs of local-maximum leaves.
	 *
	 * Remark: the DET requires parameter $minIntervalLength$ in input. This is estimated
	 * as a quantile of all distances between the first and the last of a number of
	 * consecutive points. Parameter $minLocalMaxDistance$ is set to twice
	 * $minIntervalLength$.
	 *
	 * @param maxDifference if there are "few" distinct values, and if such values differ
	 * by at most this quantity, the procedure assumes that they are uniformly
	 * distributed and it does not try to split them;
	 * @return -1 if no position with anomalously high value can be detected.
	 */
	public static final double largestValue(double[] histogram, int start, int end, double maxDifference) {
		final double QUANTILE = 0.8;
		final int N_CONSECUTIVE_POINTS = 2;
		final double SECOND_SELECTED_RUN_RATIO = 2.0;
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;
		int i;
		int nLocalMaximumLeaves, selectedRun, currentRun;
		double minIntervalLength, minLocalMaxDistance;

		// Collecting nonzero values
		lastValuePoint=-1;
		for (i=start; i<=end; i++) {
			if (histogram[i]==0) continue;
			lastValuePoint++;
			valuePoints[lastValuePoint].position=histogram[i];
			valuePoints[lastValuePoint].mass=1;
		}
		if (lastValuePoint<=0) return -1;
		lastValuePoint=Points.sortAndCompact(valuePoints,lastValuePoint);
		
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Histograms.largestValue> values:");
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=lastValuePoint; i++) IO.printErr(valuePoints[i]); }		
		
		
		if (lastValuePoint==0) return -1;
		if (lastValuePoint==1) return 1;
		if (lastValuePoint+1<=Points.FEW_POINTS) {
			if (valuePoints[lastValuePoint].position-valuePoints[0].position<=maxDifference) return -1;
			i=Points.getRoughThreshold(valuePoints,lastValuePoint,false,0,true);
			if (i==-1) return -1;
			return valuePoints[i].position;
		}

		// Finding the selected run of local-maximum leaves
		minIntervalLength=Math.max( (valuePoints[lastValuePoint].position-valuePoints[0].position)/MIN_INTERVAL_LENGTH_FACTOR,
									Points.distanceQuantile(valuePoints,lastValuePoint,N_CONSECUTIVE_POINTS,false,QUANTILE)
								  );
		minLocalMaxDistance=minIntervalLength*2.0;
		if (Points.areUniformlyDistributed(valuePoints,0,lastValuePoint,false,(valuePoints[lastValuePoint].position-valuePoints[0].position)/Points.DEFAULT_NBINS)) return -1;
		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(valuePoints,0,lastValuePoint,minIntervalLength,minLocalMaxDistance,false,-1,-1,-1,false,true);
		
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Histograms.largestValue> DET on valuePoints: minLocalMaxDistance="+minLocalMaxDistance);
if (IO.SHOW_STD_ERR_PRIME) { for (i=0; i<=DensityEstimationTree.lastLeaf; i++) IO.printErr(DensityEstimationTree.leaves[i]); }


		if (nLocalMaximumLeaves<2) return -1;
		DensityEstimationTree.markRunsOfLocalMaximumLeaves(valuePoints);
		if (DensityEstimationTree.nRuns<2) return -1;
		selectedRun=DensityEstimationTree.nRuns==2?1:DensityEstimationTree.getRunWithLargestDistance(valuePoints,false);
		if (selectedRun==1) {
			i=Leaves.lastLeafOfRun(DensityEstimationTree.leaves,DensityEstimationTree.lastLeaf,0);
			if (valuePoints[DensityEstimationTree.leaves[i+1].firstPoint].position<SECOND_SELECTED_RUN_RATIO*valuePoints[DensityEstimationTree.leaves[i].lastPoint].position) {
				if (DensityEstimationTree.nRuns==2) return -1;
				selectedRun++;
			}
		}
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Histograms.largestValue> selectedRun="+selectedRun+" nRuns="+DensityEstimationTree.nRuns);
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
		return valuePoints[DensityEstimationTree.leaves[i].firstPoint].position;
	}
	
	
	/**
	 * Let $histogram$ be an array of nonnegative integers. The procedure computes the 
	 * array that stores at position $i$ the number of times integer $i$ occurs in 
	 * $histogram$, and it sums such array to $out$. Values that exceed $out.length-1$ are
	 * added to the last position of $out$.
	 */
	public static final void valueDistribution(double[] histogram, int first, int last, int[] out) {
		final int length = out.length;
		for (int i=first; i<=last; i++) out[histogram[i]>=length?length-1:(int)histogram[i]]++;
	}
	
	
	/**
	 * @return the largest value in $histogram[first..last]$.
	 */
	public static final double max(double[] histogram, int first, int last) {
		int i;
		double max;
		
		max=Math.NEGATIVE_INFINITY;
		for (i=first; i<=last; i++) {
			if (histogram[i]>max) max=histogram[i];
		}
		return max;
	}
	
	
	/**
	 * @return the smallest value in $histogram[first..last]$.
	 */
	public static final double min(double[] histogram, int first, int last) {
		int i;
		double min;
		
		min=Math.POSITIVE_INFINITY;
		for (i=first; i<=last; i++) {
			if (histogram[i]<min) min=histogram[i];
		}
		return min;
	}


	/**
	 * @return the position of the first nonzero element in $histogram[first..last]$; 
	 * -1 if $histogram[first..last]$ contains only zeros.
	 */
	public static final int firstNonzero(double[] histogram, int first, int last) {
		for (int i=first; i<=last; i++) {
			if (histogram[i]!=0) return i;
		}
		return -1;
	}
	
	
	/**
	 * @return the position of the last nonzero element in $histogram[first..last]$; 
	 * -1 if $histogram[first..last]$ contains only zeros.
	 */
	public static final int lastNonzero(double[] histogram, int first, int last) {
		for (int i=last; i>=first; i--) {
			if (histogram[i]!=0) return i;
		}
		return -1;
	}


	/**
	 * @return TRUE iff there are points $i < j < k \in [start..end]$ such that 
	 * $x[j] \leq x[i]/ratio$ and $x[j] \leq x[k]/ratio$ and $x[j]>0$ 
	 * (if $direction=true$), or
	 * $x[i] \leq x[j]/ratio$ and $x[k] \leq x[j]/ratio$ and $x[i]>0$ and $x[k]>0$ 
	 * (if $direction=false$).
	 */
	public static final boolean findOutlier(double[] x, int start, int end, double ratio, boolean direction) {
		boolean found;
		int i, j;
		
		if (direction) {
			for (i=start+1; i<end; i++) {
				if (x[i]<=0) continue;
				found=false;
				for (j=i-1; j>=start; j--) {
					if (x[i]<=x[j]/ratio) {
						found=true;
						break;
					}
				}
				if (!found) continue;
				found=false;
				for (j=i+1; j<=end; j++) {
					if (x[i]<=x[j]/ratio) {
						found=true;
						break;
					}
				}
				if (found) return true;
			}
		}
		else {
			for (i=start+1; i<end; i++) {
				if (x[i]<=0) continue;
				found=false;
				for (j=i-1; j>=start; j--) {
					if (x[j]>0 && x[j]<=x[i]/ratio) {
						found=true;
						break;
					}
				}
				if (!found) continue;
				found=false;
				for (j=i+1; j<=end; j++) {
					if (x[j]>0 && x[j]<=x[i]/ratio) {
						found=true;
						break;
					}
				}
				if (found) return true;
			}
		}
		return false;
	}
	
	
	private static final double getNonzeroSurface(double[] histogram, int start, int end) {
		int i;
		double surface;
		
		surface=0.0;
		for (i=start; i<=end; i++) {
			if (histogram[i]>Leaves.ZERO) surface++;
		}
		return surface;
	}
	
	
	/**
	 * Searches for the leftmost substring of size $window$ with average at least (if 
	 * $sign=TRUE$) or less than (if $sign=FALSE$) $average$, inside the interval 
	 * $[start..end]$.
	 *
	 * @return the first position of a substring with those properties, or -1 if no such 
	 * substring can be found.
	 */
	public static final int firstSubstringWithAverage(double[] histogram, int start, int end, double average, boolean sign, int window) {
		final double THRESHOLD = average*window;
		int i;
		double sum;
		if (end-start+1<window) return -1;
		
		sum=0.0;
		for (i=start; i<Math.min(start+window,end+1); i++) sum+=histogram[i];
		if ((sign && sum>=THRESHOLD) || (!sign && sum<THRESHOLD)) return start;
		for (i=start+window; i<=end; i++) {
			sum+=histogram[i]-histogram[i-window];
			if ((sign && sum>=THRESHOLD) || (!sign && sum<THRESHOLD)) return i-window+1;
		}
		return -1;
	}
	
	
	/**
	 * @return the rightmost $j<=end$ such that for all $i \in [start+window-1..j]$, 
	 * substring $[i-window+1..i]$ has average at least (if $sign=TRUE$) or less than (if 
	 * $sign=FALSE$) $average$. The procedure returns -1 if no such substring can be 
	 * found.
	 */
	public static final int lastSubstringWithAverage(double[] histogram, int start, int end, double average, boolean sign, int window) {
		final double THRESHOLD = average*window;
		int i;
		double sum;
		if (end-start+1<window) return -1;
		
		sum=0.0;
		for (i=start; i<start+window; i++) sum+=histogram[i];
		if ((sign && sum<THRESHOLD) || (!sign && sum>=THRESHOLD)) return -1;
		for (i=start+window; i<=end; i++) {
			sum+=histogram[i]-histogram[i-window];
			if ((sign && sum<THRESHOLD) || (!sign && sum>=THRESHOLD)) return i-1;
		}
		return end;
	}


}