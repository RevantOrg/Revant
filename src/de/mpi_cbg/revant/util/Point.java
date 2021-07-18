package de.mpi_cbg.revant.util;


public class Point implements Comparable {
	/**
	 * Sorting status of all events
	 */
	public static final int UNSORTED = -1;
	public static final int POSITION = 0;
	public static int order = UNSORTED;

	/**
	 * Properties of a point
	 */
	public double position;
	public int mass;

	/**
	 * Temporary space used for clustering
	 */
	public int cluster;  // ID of the cluster assigned to this parenthesis
	public boolean isUsedForClustering;


	public void clear() {
		position=-1;
		mass=0;
		cluster=-1;
		isUsedForClustering=false;
	}


	public int getMass() {
		return mass;
	}


	public void sum(Point otherPoint) {
		mass+=otherPoint.getMass();
	}


	public void copyFrom(Point otherPoint) {
		position=otherPoint.position;
		mass=otherPoint.getMass();
		cluster=otherPoint.cluster;
		isUsedForClustering=otherPoint.isUsedForClustering;
	}


	public int compareTo(Object other) {
		Point otherPoint = (Point)other;
		if (order==POSITION) {
			if (position<otherPoint.position) return -1;
			else if (position>otherPoint.position) return 1;
		}
		return 0;
	}


	public String toString() {
		return position+","+getMass();
	}

}