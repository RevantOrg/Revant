package de.mpi_cbg.revant.util;


/**
 * The leaf of a density estimation tree and of a regression tree.
 */
public class Leaf implements Comparable {
	/**
	 * Key properties of a leaf
	 */
	public int firstPoint, lastPoint;  // First and last position of the interval
	public double value;  // Generic value associated with the interval
	public double distanceFromPrevious;  // Distance between the first point of the leaf and the last point of the previous leaf in string order

	/**
	 * Flags
	 */
	public boolean isLocalMinimum, isLocalMaximum, isHigh, isPointWithLargestMass, marked;
	public int fromHole;


	public final void clone(Leaf toLeaf) {
		// Key properties of a leaf
		toLeaf.firstPoint=firstPoint;
		toLeaf.lastPoint=lastPoint;
		toLeaf.value=value;
		toLeaf.distanceFromPrevious=distanceFromPrevious;

		// Flags
		toLeaf.isLocalMinimum=isLocalMinimum;
		toLeaf.isLocalMaximum=isLocalMaximum;
		toLeaf.isHigh=isHigh;
		toLeaf.marked=marked;
		toLeaf.isPointWithLargestMass=isPointWithLargestMass;
	}


	/**
	 * Sorting leaves in three blocks: (1) $isLocalMaximum=true$;
	 * (2) $isLocalMinimum=true$; (3) others. Inside each block, leaves are sorted by
	 * $firstPoint$.
	 */
	public final int compareTo(Object other) {
		Leaf otherLeaf = (Leaf)other;
		if (isLocalMaximum) {
			if (!otherLeaf.isLocalMaximum) return -1;
		}
		else if (isLocalMinimum) {
			if (otherLeaf.isLocalMaximum) return 1;
			if (!otherLeaf.isLocalMinimum) return -1;
		}
		else {
			if (otherLeaf.isLocalMaximum) return 1;
			if (otherLeaf.isLocalMinimum) return 1;
		}
		if (firstPoint<otherLeaf.firstPoint) return -1;
		else if (otherLeaf.firstPoint<firstPoint) return 1;
		return 0;
	}


	public String toString() {
		return "["+firstPoint+".."+lastPoint+"] isLocalMaximum="+isLocalMaximum+" isLocalMinimum="+isLocalMinimum+" value="+value+" marked="+marked+" isHigh="+isHigh;
	}


	public int getMass(Point[] points) {
		int out = 0;
		for (int i=firstPoint; i<=lastPoint; i++) out+=points[i].getMass();
		return out;
	}
	
	
	public double getMass(double[] x) {
		double out = 0;
		for (int i=firstPoint; i<=lastPoint; i++) out+=x[i];
		return out;
	}
	
}