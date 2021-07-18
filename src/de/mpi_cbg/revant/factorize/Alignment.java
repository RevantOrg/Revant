package de.mpi_cbg.revant.factorize;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;


public class Alignment implements Comparable {
	/**
	 * Sorting status of all alignments
	 */
	public static final int UNSORTED = -1;
	public static final int STARTA_ENDA = 0;
	public static final int READB_ORIENTATION_STARTA = 1;
	public static final int READB_ORIENTATION_STARTA_ENDA = 2;  // Decreasing $endA$, for testing containment.
	public static final int SHOULDBEUSED_READB_ORIENTATION_F = 3;
	public static final int IMPL_PERIODICSUBSTRING_STARTA_ENDA = 4;
	public static final int DISTANCE_STARTA = 5;
	public static final int IMPLIEDBYDENSE_READB_ORIENTATION_STARTB = 6;
	public static final int IMPL_STARTA = 7;
	public static final int ENDA = 8;
	public static final int IMPLIEDBYDENSE_INDENSE_PERIODIC_STARTA = 9;
	public static final int ID = 10;
	public static final int ENDA_REVERSE = 11;
	public static final int IMPLYINGINTERVAL_STARTA = 12;
	public static final int READB_STARTB_ORIENTATION = 13;
	public static final int IMPLIEDBYDENSE_INDENSE_STARTA = 14;
	public static int order = UNSORTED;
	public static int factorStart;  // Used when $order$ is based on $f$
	public static boolean sortByID = true;  // Enables final sorting by ID in $compareTo$.

	/**
	 * Properties of an alignment
	 */
	public int id;  // Absolute position in $Alignments.alignments$
	public PeriodicSubstring impliedByPeriodicSubstring;  // Pointer to a representative periodic substring that tells the same story as the alignment
	public DenseSubstring impliedByPrefixSubstring, impliedBySuffixSubstring;
	public DenseSubstring impliedByPrefixSubstringPrime, impliedBySuffixSubstringPrime;
	public boolean startInDenseSubstring, endInDenseSubstring;  // The interval of the alignment in readA starts/ends inside a dense substring that replicates by taking substrings of itself
	public boolean startInDenseSuffix, endInDensePrefix;  // The interval of the alignment in readA starts (ends) inside a dense substring that replicates by taking suffixes (prefixes) of itself
	public boolean shouldBeUsed;
	public boolean isContained;
	public double minDeltaDiff;  // Minimum difference between the average $diffs$ of this alignment and the average $diffs$ of a containing alignment
	public int isLeftMaximal, isRightMaximal;
	public int isLeftMaximalB, isRightMaximalB;
	public boolean startAdded, endAdded;
	public int intervalType;
	public boolean isPrefixOfPeriodicSubstring, isSuffixOfPeriodicSubstring, inPeriodicSubstring;  // Short-period only
	public boolean lowQualityStart, lowQualityEnd;
	
	/**
	 * Output intervals
	 */
	public PeriodicSubstringInterval periodicSubstringInterval;
	public DenseSubstring impliedByDenseSubstring;  // Pointer to a dense substring whose backbone the interval of the alignment in readA is part of.
	public DenseSubstring inDenseSubstring;  // The interval of the alignment in readA has a significant overlap with a dense substring that replicates by taking substrings of itself
	public AlignmentInterval mergedToInterval;
	public int minIntervalLength, minIntervalType, minIntervalID;
	public static final int INTERVAL_PERIODIC = 0;
	public static final int INTERVAL_DENSE = 1;
	public static final int INTERVAL_ALIGNMENT = 2;

	/**
	 * Temporary space
	 */
	public int distance;
	public boolean markedPrefix, markedSuffix;


	public String toString() {
		return "id="+id+"|"	+(impliedByPeriodicSubstring!=null)+","+inDenseSubstring+","+startInDenseSubstring+","+endInDenseSubstring+","+shouldBeUsed+","+isContained+","+isLeftMaximal+","+isRightMaximal;
	}
	
	
	public String toStringBoundaries() {
		return "["+startA()+".."+endA()+"] x readB="+readB()+"["+startB()+".."+endB()+"]fwd (orientation="+orientation()+") isLeftMaximal="+isLeftMaximal+" isRightMaximal="+isRightMaximal+" isLeftMaximalB="+isLeftMaximalB+" isRightMaximalB="+isRightMaximalB;
	}

	
	public String toStringPointers() {
		return "mergedToInterval="+(mergedToInterval==null?"null":mergedToInterval.hashCode())+"\n"+
			   "impliedByDenseSubstring="+(impliedByDenseSubstring==null?"null":impliedByDenseSubstring.hashCode())+"\n"+
			   "periodicSubstringInterval="+(periodicSubstringInterval==null?"null":periodicSubstringInterval.hashCode())+"\n"+
			   "impliedByPeriodicSubstring="+(impliedByPeriodicSubstring==null?"null":impliedByPeriodicSubstring.hashCode())+"\n"+
			   (impliedByPeriodicSubstring!=null?"~~~> "+(impliedByPeriodicSubstring.mergedToInterval==null?"null":impliedByPeriodicSubstring.mergedToInterval.hashCode())+" :: "+impliedByPeriodicSubstring.mergedToInterval+"\n":"")+
			   "inDenseSubstring="+(inDenseSubstring==null?"null":inDenseSubstring.hashCode());
	}


	public final void initialize(int identifier) {
		id=identifier;
		impliedByPrefixSubstring=null;
		impliedBySuffixSubstring=null;
		impliedByPrefixSubstringPrime=null;
		impliedBySuffixSubstringPrime=null;
		impliedByPeriodicSubstring=null;
		startInDenseSubstring=false;
		endInDenseSubstring=false;
		shouldBeUsed=true;
		isContained=false;
		isLeftMaximal=-1; isRightMaximal=-1;
		isLeftMaximalB=-1; isRightMaximalB=-1;
		startAdded=false; endAdded=false;
		isPrefixOfPeriodicSubstring=false;
		isSuffixOfPeriodicSubstring=false;
		
		// Intervals
		periodicSubstringInterval=null;
		impliedByDenseSubstring=null;
		inDenseSubstring=null;
		mergedToInterval=null;
	}

	
	/**
	 * Sets $is*Maximal$, $is*MaximalB$, $lowQualityStart$, $lowQualityEnd$.
	 */
	public final void setMaximality() {
		final int LOW_QUALITY_THRESHOLD = IO.quantum;
		
		Alignments.isLeftMaximal(id);
		isLeftMaximal=Alignments.tmp[0];
		isLeftMaximalB=Alignments.tmp[1];
		Alignments.isRightMaximal(id);
		isRightMaximal=Alignments.tmp[0];
		isRightMaximalB=Alignments.tmp[1];
		lowQualityStart=hasLowQualityStart(LOW_QUALITY_THRESHOLD);
		lowQualityEnd=hasLowQualityEnd(LOW_QUALITY_THRESHOLD);
	}


	private final boolean hasLowQualityStart(int distanceThreshold) {
		final int startA = Alignments.alignments[id][3];
		final int endA = Math.min(startA+distanceThreshold,Alignments.alignments[id][4]);
		return Reads.hasLowQuality(Alignments.alignments[id][0]-1,startA,endA,true);
	}
	
	
	private final boolean hasLowQualityEnd(int distanceThreshold) {
		final int endA = Alignments.alignments[id][4];
		final int startA = Math.max(endA-distanceThreshold,Alignments.alignments[id][3]);
		return Reads.hasLowQuality(Alignments.alignments[id][0]-1,startA,endA,true);
	}


	public final int getALength() {
		return Alignments.alignments[id][4]-Alignments.alignments[id][3]+1;
	}


	public final int getBLength() {
		return Alignments.alignments[id][6]-Alignments.alignments[id][5]+1;
	}


	/**
	 * @return true iff the interval of the alignment in readA intersects a dense
	 * substring that replicates by taking substrings of itself.
	 */
	public final boolean intersectsDenseSubstring() {
		return startInDenseSubstring || endInDenseSubstring || inDenseSubstring!=null;
	}


	public final int compareTo(Object other) {
		Alignment otherAlignment = (Alignment)other;
		boolean equalToSubstringType, equalToSubstringType_other;
		DenseSubstring tmpSubstring, otherSubstring;
		
		if (order==STARTA_ENDA) {
			if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
			else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
			if (Alignments.alignments[id][4]<Alignments.alignments[otherAlignment.id][4]) return -1;
			else if (Alignments.alignments[otherAlignment.id][4]<Alignments.alignments[id][4]) return 1;
		}
		else if (order==ENDA) {
			if (Alignments.alignments[id][4]<Alignments.alignments[otherAlignment.id][4]) return -1;
			else if (Alignments.alignments[otherAlignment.id][4]<Alignments.alignments[id][4]) return 1;
		}
		else if (order==ENDA_REVERSE) {
			if (Alignments.alignments[id][4]<Alignments.alignments[otherAlignment.id][4]) return 1;
			else if (Alignments.alignments[otherAlignment.id][4]<Alignments.alignments[id][4]) return -1;
		}
		else if (order==IMPL_PERIODICSUBSTRING_STARTA_ENDA) {
			if (impliedByPeriodicSubstring==null && otherAlignment.impliedByPeriodicSubstring!=null) return -1;
			else if (impliedByPeriodicSubstring!=null && otherAlignment.impliedByPeriodicSubstring==null) return 1;
			if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
			else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
			if (Alignments.alignments[id][4]<Alignments.alignments[otherAlignment.id][4]) return -1;
			else if (Alignments.alignments[otherAlignment.id][4]<Alignments.alignments[id][4]) return 1;
		}
		else if (order==READB_ORIENTATION_STARTA) {
			if (Alignments.alignments[id][1]<Alignments.alignments[otherAlignment.id][1]) return -1;
			else if (Alignments.alignments[otherAlignment.id][1]<Alignments.alignments[id][1]) return 1;
			if (Alignments.alignments[id][2]==1 && Alignments.alignments[otherAlignment.id][2]==0) return -1;
			else if (Alignments.alignments[id][2]==0 && Alignments.alignments[otherAlignment.id][2]==1) return 1;
			if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
			else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
		}
		else if (order==READB_ORIENTATION_STARTA_ENDA) {
			if (Alignments.alignments[id][1]<Alignments.alignments[otherAlignment.id][1]) return -1;
			else if (Alignments.alignments[otherAlignment.id][1]<Alignments.alignments[id][1]) return 1;
			if (Alignments.alignments[id][2]==1 && Alignments.alignments[otherAlignment.id][2]==0) return -1;
			else if (Alignments.alignments[id][2]==0 && Alignments.alignments[otherAlignment.id][2]==1) return 1;
			if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
			else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
			if (Alignments.alignments[id][4]>Alignments.alignments[otherAlignment.id][4]) return -1;
			else if (Alignments.alignments[otherAlignment.id][4]>Alignments.alignments[id][4]) return 1;
		}
		else if (order==SHOULDBEUSED_READB_ORIENTATION_F) {
			if (shouldBeUsed && (!otherAlignment.shouldBeUsed)) return -1;
			else if ((!shouldBeUsed) && otherAlignment.shouldBeUsed) return 1;
			if (Alignments.alignments[id][1]<Alignments.alignments[otherAlignment.id][1]) return -1;
			else if (Alignments.alignments[otherAlignment.id][1]<Alignments.alignments[id][1]) return 1;
			if (Alignments.alignments[id][2]==1 && Alignments.alignments[otherAlignment.id][2]==0) return -1;
			else if (Alignments.alignments[id][2]==0 && Alignments.alignments[otherAlignment.id][2]==1) return 1;
			int f = Alignments.getF(id,factorStart);
			int otherF = Alignments.getF(otherAlignment.id,factorStart);
			if (f<otherF) return -1;
			else if (f>otherF) return 1;
		}
		else if (order==DISTANCE_STARTA) {
			if (distance<otherAlignment.distance) return -1;
			else if (distance>otherAlignment.distance) return 1;
			if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
			else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
		}
		else if (order==IMPLIEDBYDENSE_READB_ORIENTATION_STARTB) {
			if (impliedByDenseSubstring!=null && otherAlignment.impliedByDenseSubstring==null) return -1;
			else if (impliedByDenseSubstring==null && otherAlignment.impliedByDenseSubstring!=null) return 1;
			else if (impliedByDenseSubstring!=null && otherAlignment.impliedByDenseSubstring!=null) {
				if (impliedByDenseSubstring.id<otherAlignment.impliedByDenseSubstring.id) return -1;
				else if (otherAlignment.impliedByDenseSubstring.id<impliedByDenseSubstring.id) return 1;
			}
			if (Alignments.alignments[id][1]<Alignments.alignments[otherAlignment.id][1]) return -1;
			else if (Alignments.alignments[otherAlignment.id][1]<Alignments.alignments[id][1]) return 1;
			if (Alignments.alignments[id][2]==1 && Alignments.alignments[otherAlignment.id][2]==0) return -1;
			else if (Alignments.alignments[id][2]==0 && Alignments.alignments[otherAlignment.id][2]==1) return 1;
			if (Alignments.alignments[id][5]<Alignments.alignments[otherAlignment.id][5]) return -1;
			else if (Alignments.alignments[otherAlignment.id][5]<Alignments.alignments[id][5]) return 1;
		}
		else if (order==IMPL_STARTA) {
			// Remark: alignments that are approx. identical to a dense substring of
			// substring type are put after those that are not.
			if ( (impliedByDenseSubstring==null && impliedByPeriodicSubstring==null && !equalToSubstringType()) && 
				 (otherAlignment.impliedByDenseSubstring!=null || otherAlignment.impliedByPeriodicSubstring!=null || otherAlignment.equalToSubstringType()) 
			   ) return -1;
			else if ( (impliedByDenseSubstring!=null || impliedByPeriodicSubstring!=null || equalToSubstringType()) && 
				      (otherAlignment.impliedByDenseSubstring==null && otherAlignment.impliedByPeriodicSubstring==null && !otherAlignment.equalToSubstringType())
					) return 1;
			if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
			else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
		}
		else if (order==IMPLIEDBYDENSE_INDENSE_PERIODIC_STARTA) {
			if (impliedByPeriodicSubstring!=null && otherAlignment.impliedByPeriodicSubstring==null) return 1;
			if (impliedByPeriodicSubstring==null && otherAlignment.impliedByPeriodicSubstring!=null) return -1;
			if (inDenseSubstring==null && otherAlignment.inDenseSubstring!=null) return -1;
			if (inDenseSubstring!=null && otherAlignment.inDenseSubstring==null) return 1;
			if (impliedByDenseSubstring!=null && otherAlignment.impliedByDenseSubstring==null) return -1;
			if (impliedByDenseSubstring==null && otherAlignment.impliedByDenseSubstring!=null) return 1;
			if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
			else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
		}
		else if (order==IMPLIEDBYDENSE_INDENSE_STARTA) {
			if (impliedByDenseSubstring!=null && otherAlignment.impliedByDenseSubstring==null) return -1;
			if (impliedByDenseSubstring==null && otherAlignment.impliedByDenseSubstring!=null) return 1;
			if (impliedByDenseSubstring!=null && otherAlignment.impliedByDenseSubstring!=null) {
				if (impliedByDenseSubstring.id<otherAlignment.impliedByDenseSubstring.id) return -1;
				else if (impliedByDenseSubstring.id>otherAlignment.impliedByDenseSubstring.id) return 1;
			}
			if (inDenseSubstring!=null && otherAlignment.inDenseSubstring==null) return -1;
			if (inDenseSubstring==null && otherAlignment.inDenseSubstring!=null) return 1;
			if (inDenseSubstring!=null && otherAlignment.inDenseSubstring!=null) {
				if (inDenseSubstring.id<otherAlignment.inDenseSubstring.id) return -1;
				else if (inDenseSubstring.id>otherAlignment.inDenseSubstring.id) return 1;
			}
			if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
			else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
			if (Alignments.alignments[id][4]<Alignments.alignments[otherAlignment.id][4]) return -1;
			else if (Alignments.alignments[otherAlignment.id][4]<Alignments.alignments[id][4]) return 1;
		}
		else if (order==IMPLYINGINTERVAL_STARTA) {
			// Order 1: dense, alignmentInterval, periodic, not implied.
			// Order 2: alignment startA.
			if (impliedByDenseSubstring!=null || inDenseSubstring!=null) {
				if (otherAlignment.impliedByDenseSubstring==null && otherAlignment.inDenseSubstring==null) return -1;
				tmpSubstring=impliedByDenseSubstring!=null?impliedByDenseSubstring:inDenseSubstring;
				otherSubstring=otherAlignment.impliedByDenseSubstring!=null?otherAlignment.impliedByDenseSubstring:otherAlignment.inDenseSubstring;
				if (tmpSubstring.id<otherSubstring.id) return -1;
				else if (tmpSubstring.id>otherSubstring.id) return 1;
				if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
				else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
			}
			else if (mergedToInterval!=null) {
				if (otherAlignment.impliedByDenseSubstring!=null || otherAlignment.inDenseSubstring!=null) return 1;
				else if (otherAlignment.mergedToInterval==null) return -1;
				if (mergedToInterval.id<otherAlignment.mergedToInterval.id) return -1;
				else if (mergedToInterval.id>otherAlignment.mergedToInterval.id) return 1;
				if (Alignments.alignments[id][3]<Alignments.alignments[otherAlignment.id][3]) return -1;
				else if (Alignments.alignments[otherAlignment.id][3]<Alignments.alignments[id][3]) return 1;
			}
			else if (impliedByPeriodicSubstring!=null) {
				if (otherAlignment.impliedByDenseSubstring!=null || otherAlignment.inDenseSubstring!=null || otherAlignment.mergedToInterval!=null) return 1;
				else if (otherAlignment.impliedByPeriodicSubstring==null) return -1;
			}
			else {
				// The alignment is implied by no interval
				if (otherAlignment.impliedByDenseSubstring!=null || otherAlignment.inDenseSubstring!=null || otherAlignment.mergedToInterval!=null || otherAlignment.impliedByPeriodicSubstring!=null) return 1;
			}
		}
		else if (order==READB_STARTB_ORIENTATION) {
			if (Alignments.alignments[id][1]<Alignments.alignments[otherAlignment.id][1]) return -1;
			else if (Alignments.alignments[otherAlignment.id][1]<Alignments.alignments[id][1]) return 1;
			if (Alignments.alignments[id][5]<Alignments.alignments[otherAlignment.id][5]) return -1;
			else if (Alignments.alignments[otherAlignment.id][5]<Alignments.alignments[id][5]) return 1;
			if (Alignments.alignments[id][2]==1 && Alignments.alignments[otherAlignment.id][2]==0) return -1;
			else if (Alignments.alignments[id][2]==0 && Alignments.alignments[otherAlignment.id][2]==1) return 1;
		}
		
		// Final sorting by $id$, to make sure we can exactly reproduce an order after
		// having sorted by a different criterion.
		if (sortByID) {
			if (id<otherAlignment.id) return -1;
			else if (id>otherAlignment.id) return 1;
		}
		return 0;
	}


	public final int readB() {
		return Alignments.alignments[id][1];
	}
	
	
	public final boolean orientation() {
		return Alignments.alignments[id][2]==1;
	}
	

	public final int startA() {
		return Alignments.alignments[id][3];
	}
	
	
	public final int endA() {
		return Alignments.alignments[id][4];
	}
	
	
	public final int startB() {
		return Alignments.alignments[id][5];
	}
	
	
	public final int endB() {
		return Alignments.alignments[id][6];
	}
	
	
	public boolean isImplied() {
		return mergedToInterval!=null || impliedByDenseSubstring!=null || periodicSubstringInterval!=null || inDenseSubstring!=null;
	}
	
	
	public double errorRate() {
		return ((double)(Alignments.alignments[id][7]<<1))/(getALength()+getBLength());
	}
	
	
	public boolean equalToSubstringType() {
		return inDenseSubstring!=null && Intervals.areApproximatelyIdentical(Alignments.alignments[id][3],Alignments.alignments[id][4],inDenseSubstring.startA,inDenseSubstring.endA);
	}
	
}