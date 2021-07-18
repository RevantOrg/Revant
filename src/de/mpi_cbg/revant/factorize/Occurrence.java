package de.mpi_cbg.revant.factorize;


public class Occurrence implements Comparable {
	/**
	 * Sorting status of all occurrences
	 */
	public static final int UNSORTED = -1;
	public static final int FIRSTPOSITION = 0;
	public static int order = UNSORTED;

	/**
	 * Properties of an occurrence
	 */
	public boolean orientation;
	public int firstPositionInForwardOrientation, lastPositionInForwardOrientation;
	public boolean isExact;
	public int type;  // One of the nonperiodic coverage types of $Factor$


	public final int compareTo(Object other) {
		Occurrence otherOccurrence = (Occurrence)other;
		if (order==FIRSTPOSITION) {
			if (firstPositionInForwardOrientation<otherOccurrence.firstPositionInForwardOrientation) return -1;
			else if (firstPositionInForwardOrientation>otherOccurrence.firstPositionInForwardOrientation) return 1;
		}
		return 0;
	}


	public String toString() {
		return (orientation?"1":"0")+","+firstPositionInForwardOrientation+","+lastPositionInForwardOrientation;
	}

}