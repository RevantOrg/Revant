package de.mpi_cbg.revant.util;


public class Constants {
	/**
	 * Types of interval
	 */
	public static final int INTERVAL_ALIGNMENT = 0;
	public static final int INTERVAL_DENSE_PREFIX = 1;
	public static final int INTERVAL_DENSE_SUFFIX = 2;
	public static final int INTERVAL_DENSE_PREFIXSUFFIX = 3;
	public static final int INTERVAL_DENSE_SUBSTRING = 4;
	public static final int INTERVAL_DENSE_SINGLEDELETION = 5;
	public static final int INTERVAL_PERIODIC = 6;
	public static final int INTERVAL_UNKNOWN = 7;
	
	/**
	 * Overlap IDs, also used as bit masks.
	 */
	public static final int OVERLAP_PREFIX_PREFIX = 0x00000001;
	public static final int OVERLAP_PREFIX_SUFFIX = (OVERLAP_PREFIX_PREFIX)<<1;
	public static final int OVERLAP_SUFFIX_PREFIX = (OVERLAP_PREFIX_PREFIX)<<2;
	public static final int OVERLAP_SUFFIX_SUFFIX = (OVERLAP_PREFIX_PREFIX)<<3;
	public static final int OVERLAP_MASK = OVERLAP_PREFIX_PREFIX|OVERLAP_PREFIX_SUFFIX|OVERLAP_SUFFIX_PREFIX|OVERLAP_SUFFIX_SUFFIX;
	
	/**
	 * Types of interval graph edge, also used as bit masks.
	 */
	public static final int CONTAINMENT_IDENTICAL = (OVERLAP_SUFFIX_SUFFIX)<<1;  // =16
	public static final int CONTAINMENT_ONE_IN_TWO = (CONTAINMENT_IDENTICAL)<<1;
	public static final int CONTAINMENT_TWO_IN_ONE = (CONTAINMENT_IDENTICAL)<<2;
	public static final int INSERTION_ONE_IN_TWO = (CONTAINMENT_IDENTICAL)<<3;
	public static final int INSERTION_TWO_IN_ONE = (CONTAINMENT_IDENTICAL)<<4;
	public static final int INSERTION_BOTH = (CONTAINMENT_IDENTICAL)<<5;
	public static final int SHARED_SUBSTRING = (CONTAINMENT_IDENTICAL)<<6;
	public static final int CONTAINMENT_MASK = CONTAINMENT_IDENTICAL|CONTAINMENT_ONE_IN_TWO|CONTAINMENT_TWO_IN_ONE;
	public static final int INSERTION_MASK = INSERTION_ONE_IN_TWO|INSERTION_TWO_IN_ONE|INSERTION_BOTH;
	
	/**
	 * Interval graph containment subtypes
	 */
	public static final int CONTAINMENT_PREFIX = 0;
	public static final int CONTAINMENT_SUFFIX = 1;
	public static final int CONTAINMENT_SUBSTRING = 2;
	
	
	/**
	 * Remark: the procedure handles also periodic types.
	 *
	 * @param hasLongPeriod* ignored if $type*$ is not periodic;
	 * @return true iff $type1$ can be the type of a substring of an interval of type 
	 * $type2$ that occurs in another read.
	 */
	public static final boolean compareTypes(int type1, int type2, boolean hasLongPeriod1, boolean hasLongPeriod2, int orientation) {
		if (type2==INTERVAL_ALIGNMENT) {
			return type1==INTERVAL_ALIGNMENT;
		}
		else if (type2==INTERVAL_DENSE_PREFIX) {
			if (orientation==0) return type1==INTERVAL_DENSE_PREFIX || type1==INTERVAL_ALIGNMENT;
			else if (orientation==1) return type1==INTERVAL_DENSE_SUFFIX || type1==INTERVAL_ALIGNMENT;
			else return type1==INTERVAL_DENSE_PREFIX || type1==INTERVAL_DENSE_SUFFIX || type1==Constants.INTERVAL_ALIGNMENT;
		}
		else if (type2==INTERVAL_DENSE_SUFFIX) {
			if (orientation==0) return type1==INTERVAL_DENSE_SUFFIX || type1==INTERVAL_ALIGNMENT;
			else if (orientation==1) return type1==INTERVAL_DENSE_PREFIX || type1==INTERVAL_ALIGNMENT;
			else return type1==INTERVAL_DENSE_PREFIX || type1==INTERVAL_DENSE_SUFFIX || type1==Constants.INTERVAL_ALIGNMENT;
		}
		else if (type2==INTERVAL_DENSE_PREFIXSUFFIX) {
			return type1==INTERVAL_DENSE_PREFIXSUFFIX || type1==INTERVAL_DENSE_PREFIX || type1==INTERVAL_DENSE_SUFFIX || type1==INTERVAL_ALIGNMENT;
		}
		else if (type2==INTERVAL_DENSE_SINGLEDELETION) {
			return type1==INTERVAL_DENSE_SINGLEDELETION || type1==INTERVAL_DENSE_PREFIXSUFFIX || type1==INTERVAL_DENSE_PREFIX || type1==INTERVAL_DENSE_SUFFIX || type1==INTERVAL_ALIGNMENT;
		}
		else if (type2==INTERVAL_DENSE_SUBSTRING) {
			return type1==INTERVAL_DENSE_SUBSTRING || type1==INTERVAL_DENSE_PREFIXSUFFIX || type1==INTERVAL_DENSE_PREFIX || type1==INTERVAL_DENSE_SUFFIX || type1==INTERVAL_ALIGNMENT;
		}
		else if (type2==INTERVAL_PERIODIC) {
			if (hasLongPeriod2) {
				return type1==INTERVAL_ALIGNMENT || (type1==INTERVAL_PERIODIC && hasLongPeriod1);
			}
			else {
				return type1==INTERVAL_PERIODIC && !hasLongPeriod1;
			}
		}
		return false;
	}
	
	
	/**
	 * @return a version of interval-type $type$ after putting the interval in the forward
	 * (if $orientation=0$) or reverse (if $orientation=1$) orientation; $orientation=2$ 
	 * is treated like $orientation=0$.
	 */
	public static final int typeInOrientation(int type, int orientation) {
		if (type==INTERVAL_DENSE_PREFIX && orientation==1) return INTERVAL_DENSE_SUFFIX;
		else if (type==INTERVAL_DENSE_SUFFIX && orientation==1) return INTERVAL_DENSE_PREFIX;
		else return type;
	}
	
	
	/**
	 * @return the type that comes from merging an interval of $type1$ and an interval of 
	 * $type2$, or -2 if merging the two interval types is not valid.
	 */
	public static final int mergeTypes(int type1, int type2) {
		if (type1==type2) return type1;
		final int typeA = Math.min(type1,type2);
		final int typeB = Math.max(type1,type2);
		
		if (typeA==INTERVAL_ALIGNMENT) return -2;
		else if (typeA==INTERVAL_DENSE_PREFIX || typeA==INTERVAL_DENSE_SUFFIX) {
			switch (typeB) {
				case INTERVAL_DENSE_PREFIXSUFFIX: return INTERVAL_DENSE_PREFIXSUFFIX;
				case INTERVAL_DENSE_SUBSTRING: return INTERVAL_DENSE_SUBSTRING;
				case INTERVAL_DENSE_SINGLEDELETION: return INTERVAL_DENSE_SINGLEDELETION;
				default: return -2;
			}
		}
		else if (typeA==INTERVAL_DENSE_PREFIXSUFFIX) {
			switch (typeB) {
				case INTERVAL_DENSE_SUBSTRING: return INTERVAL_DENSE_SUBSTRING;
				case INTERVAL_DENSE_SINGLEDELETION: return INTERVAL_DENSE_SINGLEDELETION;
				default: return -2;
			}
		}
		else if (typeA==INTERVAL_DENSE_SUBSTRING) {
			if (typeB==INTERVAL_DENSE_SINGLEDELETION) return INTERVAL_DENSE_SUBSTRING;
			else return -2;
		}
		else return -2;
	}
	
	
	/**
	 * @return the type of overlap seen when going in the opposite direction through 
	 * $type$.
	 */
	public static final int oppositeDirectionOverlap(int type) {
		switch (type) {
			case OVERLAP_PREFIX_PREFIX: return OVERLAP_PREFIX_PREFIX;
			case OVERLAP_PREFIX_SUFFIX: return OVERLAP_SUFFIX_PREFIX;
			case OVERLAP_SUFFIX_PREFIX: return OVERLAP_PREFIX_SUFFIX;
			case OVERLAP_SUFFIX_SUFFIX: return OVERLAP_SUFFIX_SUFFIX;
		}
		return -1;
	}
	
	
	/**
	 * @return the overlap obtained from $type$ when interval 1 is reverse-complemented
	 * ($rc=0$), when interval 2 is reverse-complemented ($rc=1$), or when both are
	 * reverse-complemented ($rc=2$).
	 */
	public static final int reverseComplementOverlap(int type, int rc) {
		if (type==OVERLAP_PREFIX_PREFIX) {
			switch (rc) {
				case 0: return OVERLAP_SUFFIX_PREFIX;
				case 1: return OVERLAP_PREFIX_SUFFIX;
				case 2: return OVERLAP_SUFFIX_SUFFIX;
			}
		}
		else if (type==OVERLAP_PREFIX_SUFFIX) {
			switch (rc) {
				case 0: return OVERLAP_SUFFIX_SUFFIX;
				case 1: return OVERLAP_PREFIX_PREFIX;
				case 2: return OVERLAP_SUFFIX_PREFIX;
			}
		}
		else if (type==OVERLAP_SUFFIX_PREFIX) {
			switch (rc) {
				case 0: return OVERLAP_PREFIX_PREFIX;
				case 1: return OVERLAP_SUFFIX_SUFFIX;
				case 2: return OVERLAP_PREFIX_SUFFIX;
			}
		}
		else if (type==OVERLAP_SUFFIX_SUFFIX) {
			switch (rc) {
				case 0: return OVERLAP_PREFIX_SUFFIX;
				case 1: return OVERLAP_SUFFIX_PREFIX;
				case 2: return OVERLAP_PREFIX_PREFIX;
			}
		}
		return -1;
	}
	
	
	public static final int overlap2periodicOverlap(int overlap) {
		return overlap<<4;
	}
	
	
	/**
	 * @param id in {0,1,2,3}.
	 */
	public static final int id2overlap(int id) {
		return (OVERLAP_PREFIX_PREFIX)<<id;
	}
	
	
	/**
	 * @return a number in {0,1,2,3}.
	 */
	public static final int overlap2id(int overlap) {
		return Math.numberOfTrailingZeros(overlap);
	}
	
	
	/**
	 * @return 0=same orientation; 1=opposite; -1=unknown.
	 */
	public static final int overlap2orientation(int overlap) {
		switch (overlap) {
			case OVERLAP_PREFIX_PREFIX: return 1;
			case OVERLAP_PREFIX_SUFFIX: return 0;
			case OVERLAP_SUFFIX_PREFIX: return 0;
			case OVERLAP_SUFFIX_SUFFIX: return 1;
		}
		return -1;
	}
	
	
	/**
	 * @param sideX TRUE=the overlap involves a prefix of $nodeX$;
	 * FALSE=the overlap involves a suffix of $nodeX$.
	 */
	public static final int sidesToOverlap(boolean side1, boolean side2) {
		if (side1) {
			if (side2) return OVERLAP_PREFIX_PREFIX;
			else return OVERLAP_PREFIX_SUFFIX;
		}
		else {
			if (side2) return OVERLAP_SUFFIX_PREFIX;
			else return OVERLAP_SUFFIX_SUFFIX;
		}
	}
	
	
	/**
	 * @return the type of containment seen when going in the opposite direction through 
	 * $type$.
	 */
	public static final int oppositeDirectionContainment(int type) {
		switch (type) {
			case CONTAINMENT_ONE_IN_TWO: return CONTAINMENT_TWO_IN_ONE;
			case CONTAINMENT_TWO_IN_ONE: return CONTAINMENT_ONE_IN_TWO;
			case CONTAINMENT_IDENTICAL: return CONTAINMENT_IDENTICAL;
		}
		return -1;
	}
	
	
	/**
	 * @return the type of insertion seen when going in the opposite direction through 
	 * $type$.
	 */
	public static final int oppositeDirectionInsertion(int type) {
		switch (type) {
			case INSERTION_ONE_IN_TWO: return INSERTION_TWO_IN_ONE;
			case INSERTION_TWO_IN_ONE: return INSERTION_ONE_IN_TWO;
			case INSERTION_BOTH: return INSERTION_BOTH;
		}
		return -1;
	}
	
	
}