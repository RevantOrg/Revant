package de.mpi_cbg.revant.factorize;

import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;


public class PeriodicSubstring implements Comparable {
	/**
	 * Sorting status of all periodic substrings
	 */
	public static final int UNSORTED = -1;
	public static final int MINSTARTA = 0;
	public static final int MINSTARTA_MAXENDA = 1;
	public static final int READB_ORIENTATION_MINSTARTA = 2;
	public static final int READB_MINSTARTB = 3;
	public static final int SHOULDBECOUNTED_READB_MINSTARTB = 4;
	public static final int READB_MINSTARTB_FORWARD = 5;
	public static final int READB_ORIENTATION_MINSTARTB = 6;
	public static final int LONGPERIOD_SHORTPERIOD = 7;
	public static int order = UNSORTED;
	public static int order_longPeriod = UNSORTED;

	/**
	 * Properties of a periodic substring
	 */
	public int minStartA, maxStartA;  // Minimum/maximum starting position in $readA$
	public int nStartA, sumStartA;
	public int minEndA, maxEndA;  // Minimum/maximum ending position in $readA$
	public int nEndA, sumEndA;
	public int nAlignments;
	public int minStartB;  // Minimum starting position in $readB$ considered in orientation $orientation$
	public int maxEndB;  // Maximum ending position in $readB$ considered in orientation $orientation$
	public boolean orientation;  // Orientation of $readB$: forward (TRUE) or reverse (FALSE).
	public int readB;  // ID of $readB$. The ID of $readA$ is implicit from the context.
    public boolean shouldBeCountedForCoverage;  // Decides whether this substring is counted in the corrected coverage
    public PeriodicSubstring threadStart;  // Pointer to a periodic substring that represents the thread of the current periodic substring. A \emph{thread} is a maximal sequence of periodic substrings that highly overlap in readA, and that are separated in readB by just random insertions.
    public boolean equalsReverseComplement;
    public int minStartBForward, maxEndBForward;  // $minStartB$ and $maxEndB$ in the forward orientation
    public int pathLength;
	public boolean isLeftMaximal, isRightMaximal;  // Marks the presence of nonrandom sequence to the left/right of the periodic substring in readB
	public Point[] shifts;
	public int lastShift;
	public boolean isStraddling;
	public PeriodicSubstringInterval mergedToInterval;
	public int period;
	public boolean hasLongPeriod;
	public boolean equalsOnePeriod;
	public int longestAlignment;

    /**
     * Connected component variables
     */
    public PeriodicSubstring representative;  // Pointer to the representative of the connected component of this substring. Can be equal to $this$.

	/**
	 * Flags
	 */
	public boolean containsExactOccurrence;  // Marks whether the periodic substring contains and exact occurrence of a factor
	public boolean shouldBeUsedForFactoring;
	
	/**
	 * Temporary space
	 */
	public int tmpInt;


	public PeriodicSubstring() { 
		// NOP
	}
	
	
	/**
	 * Creates an artificial periodic substring, copying the values from a periodic 
	 * substring interval.
	 */
	public PeriodicSubstring(PeriodicSubstringInterval interval) {
		setFromInterval(interval);
	}
	
	
	public final void setFromInterval(PeriodicSubstringInterval interval) {
		minStartA=interval.firstPosition;
		maxStartA=minStartA;
		nStartA=1;
		sumStartA=minStartA;
		minEndA=interval.lastPosition;
		maxEndA=minEndA;
		nEndA=1;
		sumEndA=minEndA;
		minStartB=-1;
		maxEndB=-1;
		orientation=true;
		readB=-1;
		shouldBeCountedForCoverage=false;
		threadStart=null;
		equalsReverseComplement=false;
		minStartBForward=-1;
		maxEndBForward=-1;
		pathLength=0;
		isLeftMaximal=interval.isLeftMaximal;
		isRightMaximal=interval.isRightMaximal;
		lastShift=-1;
		isStraddling=false;
		mergedToInterval=interval;
		period=interval.period;
		hasLongPeriod=interval.hasLongPeriod;
		representative=this;
		containsExactOccurrence=false;
		shouldBeUsedForFactoring=false;
		equalsOnePeriod=interval.equalsOnePeriod;
	}


	/**
	 * @param maxShifts upper bound on the number of shifts to estimate the period.
	 */
	public final void allocateMemory(int maxShifts) {
		shifts = new Point[maxShifts];
		for (int i=0; i<shifts.length; i++) shifts[i] = new Point();
		lastShift=-1;
	}
	
	
	/**
	 * Copies all fields of $from$ onto this substring.
	 */
	public final void clone(PeriodicSubstring from) {
		minStartA=from.minStartA;
		maxStartA=from.maxStartA;
		nStartA=from.nStartA;
		sumStartA=from.sumStartA;
		minEndA=from.minEndA;
		maxEndA=from.maxEndA;
		nEndA=from.nEndA;
		sumEndA=from.sumEndA;
		minStartB=from.minStartB;
		maxEndB=from.maxEndB;
		orientation=from.orientation;
		readB=from.readB;
		shouldBeCountedForCoverage=from.shouldBeCountedForCoverage;
		threadStart=from.threadStart;
		equalsReverseComplement=from.equalsReverseComplement;
		minStartBForward=from.minStartBForward;
		maxEndBForward=from.maxEndBForward;
		pathLength=from.pathLength;
		isLeftMaximal=from.isLeftMaximal;
		isRightMaximal=from.isRightMaximal;
		if (shifts==null || shifts.length<from.lastShift+1) {
			shifts = new Point[from.lastShift+1];
			for (int i=0; i<shifts.length; i++) shifts[i] = new Point();
		}
		Points.simpleClone(from.shifts,from.lastShift,shifts);
		lastShift=from.lastShift;
		isStraddling=from.isStraddling;
		mergedToInterval=from.mergedToInterval;
		period=from.period;
		hasLongPeriod=from.hasLongPeriod;
		representative=from.representative;
		containsExactOccurrence=from.containsExactOccurrence;
		shouldBeUsedForFactoring=from.shouldBeUsedForFactoring;
		equalsOnePeriod=from.equalsOnePeriod;
	}


    public final int compareTo(Object other) {
		PeriodicSubstring otherSubstring = (PeriodicSubstring)other;
    	int ord;
		
		if (!hasLongPeriod) ord=order;
		else ord=order_longPeriod;
		
    	if (ord==MINSTARTA) {
			if (minStartA<otherSubstring.minStartA) return -1;
			else if (minStartA>otherSubstring.minStartA) return 1;
		}
		else if (ord==MINSTARTA_MAXENDA) {
			if (minStartA<otherSubstring.minStartA) return -1;
			else if (minStartA>otherSubstring.minStartA) return 1;
			if (maxEndA<otherSubstring.maxEndA) return -1;
			else if (maxEndA>otherSubstring.maxEndA) return 1;
		}
		else if (ord==READB_ORIENTATION_MINSTARTA) {
			if (readB<otherSubstring.readB) return -1;
			else if (readB>otherSubstring.readB) return 1;
			if (orientation && (!otherSubstring.orientation)) return -1;
			else if ((!orientation) && otherSubstring.orientation) return 1;
			if (minStartA<otherSubstring.minStartA) return -1;
			else if (minStartA>otherSubstring.minStartA) return 1;
		}
		else if (ord==READB_MINSTARTB) {
			if (readB<otherSubstring.readB) return -1;
			else if (readB>otherSubstring.readB) return 1;
			if (minStartB<otherSubstring.minStartB) return -1;
			else if (minStartB>otherSubstring.minStartB) return 1;
		}
		else if (ord==SHOULDBECOUNTED_READB_MINSTARTB) {
			if (shouldBeCountedForCoverage && !otherSubstring.shouldBeCountedForCoverage) return -1;
			else if (!shouldBeCountedForCoverage && otherSubstring.shouldBeCountedForCoverage) return 1;
			if (readB<otherSubstring.readB) return -1;
			else if (readB>otherSubstring.readB) return 1;
			if (minStartB<otherSubstring.minStartB) return -1;
			else if (minStartB>otherSubstring.minStartB) return 1;
		}
		else if (ord==READB_MINSTARTB_FORWARD) {
			if (readB<otherSubstring.readB) return -1;
			else if (readB>otherSubstring.readB) return 1;
			if (minStartBForward<otherSubstring.minStartBForward) return -1;
			else if (minStartBForward>otherSubstring.minStartBForward) return 1;
		}
		else if (ord==READB_ORIENTATION_MINSTARTB) {
			if (readB<otherSubstring.readB) return -1;
			else if (readB>otherSubstring.readB) return 1;
			if (orientation && (!otherSubstring.orientation)) return -1;
			else if ((!orientation) && otherSubstring.orientation) return 1;
			if (minStartB<otherSubstring.minStartB) return -1;
			else if (minStartB>otherSubstring.minStartB) return 1;
		}
		else if (ord==LONGPERIOD_SHORTPERIOD) {
			if (hasLongPeriod && !otherSubstring.hasLongPeriod) return -1;
			else if (!hasLongPeriod && otherSubstring.hasLongPeriod) return 1;
		}
		return 0;
    }


    public String toString() {
    	return "("+minStartA+","+maxStartA+")..("+minEndA+","+maxEndA+") x readB="+readB+"["+minStartBForward+".."+maxEndBForward+"]"+orientation+", shouldBeCountedForCoverage="+shouldBeCountedForCoverage+", isRepresentative="+(representative==this)+", containsExactOccurrence="+containsExactOccurrence+" PATHLENGTH="+pathLength+" period="+period+" hasLongPeriod="+hasLongPeriod+" isLeftMaximal="+isLeftMaximal+" isRightMaximal="+isRightMaximal+" shouldBeUsedForFactoring="+shouldBeUsedForFactoring+" equalsOnePeriod="+equalsOnePeriod+" ---> "+mergedToInterval;
		
		/* New format, just for studying periodic prefix/suffix substrings:
		return ReadA.id+","+minStartA+","+maxEndA+","+readB+","+(orientation?minStartB:Reads.getReadLength(readB)-maxEndB-1)+","+(orientation?maxEndB:Reads.getReadLength(readB)-minStartB-1);
		*/
    }


	/**
	 * Tries to remove low-quality regions from the start/end of the substring.
	 * Trimming is not performed if it would create an interval of length smaller than
	 * $minLength$.
	 *
	 * Remark: $is{Left,Right}Maximal$ and $isStraddling$ might get invalidated by the 
	 * procedure.
	 */
	public final void trim(int[] tmpArray, int minLength) {
		final int WINDOW_LARGE = 3;  // Arbitrary
		final int WINDOW_SMALL = 1;
		
		tmpArray[0]=minStartA; tmpArray[1]=maxEndA;
		Reads.trim(tmpArray,WINDOW_LARGE,true,true);
		if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<minLength) return;
		Reads.trim(tmpArray,WINDOW_SMALL,true,true);
		if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<minLength) return;
		if (tmpArray[0]!=minStartA) {
			minStartA=tmpArray[0];
			sumStartA=nStartA*minStartA;
			if (minStartA>maxStartA) maxStartA=minStartA;
			isLeftMaximal=false;
		}
		if (tmpArray[1]!=maxEndA) {
			maxEndA=tmpArray[1];
			sumEndA=nEndA*maxEndA;
			if (maxEndA<minEndA) minEndA=maxEndA;
			isRightMaximal=false;
		}
	}
	
	
	public final int length() {
		return maxEndA-minStartA+1;
	}
	
	
	/**
	 * Copies the values of the points in $source[from..to]$ (which are assumed to be 
	 * sorted and distinct) into $shifts$, allowing $shifts$ to be at most 
	 * $maxDestinationLength$ long. Since this amount can be significantly smaller than 
	 * $source$, the procedure greedily collapses into a single point a maximal run of 
	 * consecutive points that are at distance at most $d$ from the first point of the 
	 * run. The procedure uses the smallest value of $d$ (possibly zero) that makes 
	 * $source$ fit into $maxDestinationLength$ units.
	 * 
	 * @param tmpPoints temporary space, assumed to be large enough to handle any size of
	 * $source$.
	 */
	public final void cloneShifts(Point[] source, int from, int to, int maxDestinationLength, Point[] tmpPoints) {
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
		} while (j>=maxDestinationLength);
		if (shifts==null || shifts.length<j+1) {
			shifts = new Point[j+1];
			for (i=0; i<shifts.length; i++) shifts[i] = new Point();
		}
		for (i=0; i<=j; i++) {
			shifts[i].position=tmpPoints[i].position;
			shifts[i].mass=tmpPoints[i].mass;
		}
		lastShift=j;
	}

}