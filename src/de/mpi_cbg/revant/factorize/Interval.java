package de.mpi_cbg.revant.factorize;


public class Interval implements Comparable {
	/**
	 * Sorting status of all intervals
	 */
	public static final int UNSORTED = -1;
	public static final int READB_ORIENTATION_FIRSTPOSITION = 0;
	public static final int FIRSTPOSITION = 1;
	public static int order = UNSORTED;

	/**
	 * Properties of an interval
	 */
	public int firstPosition, lastPosition, readB;
	public boolean orientation;
	public int alignmentID;
	public double avgDiffs;
	public double center;


	public final int compareTo(Object other) {
		Interval otherInterval = (Interval)other;
		if (order==READB_ORIENTATION_FIRSTPOSITION) {
			if (readB<otherInterval.readB) return -1;
			else if (readB>otherInterval.readB) return 1;
			if (orientation && !otherInterval.orientation) return -1;
			else if (!orientation && otherInterval.orientation) return 1;
			if (firstPosition<otherInterval.firstPosition) return -1;
			else if (firstPosition>otherInterval.firstPosition) return 1;
		}
		return 0;
	}


	public String toString() {
		return "readB="+readB+",orientation="+orientation+",first="+firstPosition+",last="+lastPosition;
	}

}