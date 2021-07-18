package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import de.mpi_cbg.revant.util.IO;


public class Parentheses {

	/**
	 * Given an alignment $(R[i..j],R'[i'..j'])$, we say that the pair $(i,i')$ is an
	 * open parenthesis if we can prove that the alignment is left-maximal (see procedure
	 * $IO.isLeftMaximal$). We say that the pair $(j,j')$ is a closed parenthesis if we
	 * can prove that the alignment is right-maximal (see procedure $IO.isRightMaximal$).
	 * The procedure stores in $Events.events$ all parentheses that do not belong to
	 * alignments implied by periodic or dense substrings, and that are not located inside
	 * a dense substring that replicates by taking substrings of itself, or inside a 
	 * short-period periodic interval.
	 * Parentheses are one of the ways of splitting readA. The first type of parenthesis
	 * belongs to alignments that are the backbone of the corresponding substrings, so
	 * it cannot be used for splitting such substrings. The second type of parenthesis
	 * falls inside the only type of dense substring that has to be split using function
	 * $\delta$ rather than by clustering events.
	 *
	 * We say that $R'[i'..j']$ is a \emph{witness} of $R[i..j]$ iff $(i,i')$ is an open
	 * parenthesis and $(j,j')$ is a closed parenthesis. In other words, $R$ and $R'$
	 * prove the existence of a maximal repeat in the genome. Not all occurrences of
	 * maximal repeats in $R$ have a witness: for example, consider a read $R$ in which
	 * maximal repeats $A$, $B$ and $C$ of the genome occur consecutively as $ABC$, and
	 * assume that in the genome, whenever $B$ occurs, it is either preceded by $A$ or it
	 * is followed by $C$. Then, no witness exists for the occurrence of $B$ in $R$.
	 * However, if just a substring of a maximal repeat of the genome occurs in a read
	 * $R$ (for example because $R$ sampled just a prefix of the repeat, or because of a
	 * long random insertion), then such substring cannot have a witness.
	 *
	 * Remark: a read can contain both open and closed parentheses, and such parentheses
	 * can occur in any order. A read could also contain just open or just closed
	 * parentheses.
	 *
	 * Remark: the lists of all dense and periodic substrings must have already
	 * been computed.
	 *
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ to be sorted by
	 * $firstPosition$.
	 *
	 * @param useUnknownEvents uses also the start/end of alignments whose left/right
	 * maximality is unknown;
	 * @return the number of parentheses detected.
	 */
	public static final int detect(boolean useUnknownEvents) {
		final int MAXIMAL_ALIGNMENTS_THRESHOLD = 20*IO.coverage;  // Arbitrary
		boolean previousSortByID_periodic, previousSortByID_dense;
		int i;
		int firstEvent, position, periodicFlags, substringFlags;
		int previousShortPeriodOrder, previousLongPeriodOrder, previousDenseOrder;
		PeriodicSubstringInterval tmpInterval = new PeriodicSubstringInterval();
		DenseSubstring tmpSubstring = new DenseSubstring();
		AlignmentInterval mergedToInterval;
		int[] tmpArray = new int[2];

		// Marking useful positions
		previousSortByID_periodic=PeriodicSubstringInterval.sortByID;
		previousShortPeriodOrder=PeriodicSubstringInterval.order;
		previousLongPeriodOrder=PeriodicSubstringInterval.order_longPeriod;
		PeriodicSubstringInterval.sortByID=false;
		if (PeriodicSubstringInterval.order!=PeriodicSubstringInterval.FIRSTPOSITION || PeriodicSubstringInterval.order_longPeriod!=PeriodicSubstringInterval.FIRSTPOSITION) {
			PeriodicSubstringInterval.order=PeriodicSubstringInterval.FIRSTPOSITION;
			PeriodicSubstringInterval.order_longPeriod=PeriodicSubstringInterval.FIRSTPOSITION;
			if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		}
		previousSortByID_dense=DenseSubstring.sortByID;
		previousDenseOrder=DenseSubstring.order;
		DenseSubstring.sortByID=false;
		if (DenseSubstring.order!=DenseSubstring.STARTA) {
			DenseSubstring.order=DenseSubstring.STARTA;		
			if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
		}
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null || ReadA.sortedAlignments[i].impliedByDenseSubstring!=null) {
				ReadA.sortedAlignments[i].markedPrefix=false;
				ReadA.sortedAlignments[i].markedSuffix=false;
			}
			else {
				periodicFlags=usePositions_periodic(ReadA.sortedAlignments[i],useUnknownEvents,tmpInterval,tmpArray);
				substringFlags=usePositions_substring(ReadA.sortedAlignments[i],useUnknownEvents,tmpSubstring);
				mergedToInterval=ReadA.sortedAlignments[i].mergedToInterval;
				ReadA.sortedAlignments[i].markedPrefix=(((periodicFlags&1)!=0)&&((substringFlags&1)!=0))||(mergedToInterval!=null&&mergedToInterval.nMaximalAlignmentsLeft>=MAXIMAL_ALIGNMENTS_THRESHOLD);
				ReadA.sortedAlignments[i].markedSuffix=(((periodicFlags&2)!=0)&&((substringFlags&2)!=0))||(mergedToInterval!=null&&mergedToInterval.nMaximalAlignmentsRight>=MAXIMAL_ALIGNMENTS_THRESHOLD);
			}
		}
		
		// Adding useful positions to $Events$.
		firstEvent=Events.lastEvent+1;
		for (i=0; i<=ReadA.lastSortedAlignment; i++) {
			if (ReadA.sortedAlignments[i].markedPrefix) {
				position=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				Events.lastEvent++;
				Events.ensureSpace_events(Events.lastEvent+1);
				Events.events[Events.lastEvent].clear();
				Events.events[Events.lastEvent].position=position;
				Events.events[Events.lastEvent].nOpen=1;
				ReadA.sortedAlignments[i].startAdded=true;
			}
			else ReadA.sortedAlignments[i].startAdded=false;
			if (ReadA.sortedAlignments[i].markedSuffix) {
				position=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
				Events.lastEvent++;
				Events.ensureSpace_events(Events.lastEvent+1);
				Events.events[Events.lastEvent].clear();
				Events.events[Events.lastEvent].position=position;
				Events.events[Events.lastEvent].nClosed=1;
				ReadA.sortedAlignments[i].endAdded=true;
			}
			else ReadA.sortedAlignments[i].endAdded=false;
		}
		
		// Restoring original orders
		PeriodicSubstringInterval.sortByID=previousSortByID_periodic;
		PeriodicSubstringInterval.order=previousShortPeriodOrder;
		PeriodicSubstringInterval.order_longPeriod=previousLongPeriodOrder;
		if (PeriodicSubstrings.lastInterval>0) Arrays.sort(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1);
		DenseSubstring.sortByID=previousSortByID_dense;
		DenseSubstring.order=previousDenseOrder;
		if (DenseSubstrings.lastSubstring>0) Arrays.sort(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1);
		
		return Events.lastEvent-firstEvent+1;
	}
	
	
	/**
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ to be sorted by
	 * $firstPosition$.
	 *
	 * @param tmpInterval temporary space;
	 * @return mask of flags telling whether the first (bit 0 from LSB) and the last (bit 
	 * 1 from LSB) position in readA of $alignment$ should be used to create an event, 
	 * based on their containment in periodic intervals.
	 */
	private static final int usePositions_periodic(Alignment alignment, boolean useUnknownEvents, PeriodicSubstringInterval tmpInterval, int[] tmpArray) {
		final int MIN_MAXIMAL_ALIGNMENTS = IO.coverage*10;  // Arbitrary
		final int WINDOW_LARGE = 3;  // Arbitrary
		final int WINDOW_SMALL = 1;
		
		tmpArray[0]=alignment.startA(); tmpArray[1]=alignment.endA();
		Reads.trim(tmpArray,WINDOW_LARGE,true,true);
		if (tmpArray[0]==-1 || tmpArray[1]==-1) return 0;
		Reads.trim(tmpArray,WINDOW_SMALL,true,true);
		if (tmpArray[0]==-1 || tmpArray[1]==-1) return 0;
		final int inPeriodicStart = inPeriodic(tmpArray[0],tmpInterval);
		final int inPeriodicEnd = inPeriodic(tmpArray[1],tmpInterval);
		final boolean useStart = (useUnknownEvents||alignment.isLeftMaximal==1) &&
			                       ( ((inPeriodicStart&1)==0 && (inPeriodicStart&2)==0) ||
			                         ((inPeriodicStart&1)!=0 && (inPeriodicEnd&1)==0 && (inPeriodicEnd&2)==0 && alignment.mergedToInterval!=null && alignment.mergedToInterval.nMaximalAlignmentsLeft>=MIN_MAXIMAL_ALIGNMENTS) ||
			                         ((inPeriodicStart&1)==0 && (inPeriodicStart&2)!=0 && (inPeriodicStart&4)!=0)
								   );
		final boolean useEnd = (useUnknownEvents||alignment.isRightMaximal==1) &&
			                   ( ((inPeriodicEnd&1)==0 && (inPeriodicEnd&2)==0) ||
			                     ((inPeriodicEnd&1)!=0 && (inPeriodicStart&1)==0 && (inPeriodicStart&2)==0 && alignment.mergedToInterval!=null && alignment.mergedToInterval.nMaximalAlignmentsRight>=MIN_MAXIMAL_ALIGNMENTS) ||
			                     ((inPeriodicEnd&1)==0 && (inPeriodicEnd&2)!=0 && (inPeriodicEnd&8)!=0)
							   );
		return (useStart?1:0)|(useEnd?2:0);
	}
	
	
	/**
	 * Remark: the procedure assumes $PeriodicSubstrings.intervals$ to be sorted by
	 * $firstPosition$. 
	 *
	 * @return flags telling whether $position$ is contained in a periodic interval. 
	 * Bits (from LSB): 0=short-period interval; 1=long-period interval; 2=far from the 
	 * first position of all containing long-period intervals; 3=far from the last 
	 * position of all containing long-period intervals.
	 */
	private static final int inPeriodic(int position, PeriodicSubstringInterval tmpInterval) {
		final int IDENTITY_THRESHOLD = (3*IO.quantum)>>1;  // Arbitrary
		boolean inShort, inLong, farFromFirst, farFromLast;
		int i, j;
		int out;
		PeriodicSubstringInterval interval;
		if (PeriodicSubstrings.lastInterval==-1) return 0;
		
		tmpInterval.firstPosition=position;
		i=Arrays.binarySearch(PeriodicSubstrings.intervals,0,PeriodicSubstrings.lastInterval+1,tmpInterval);
		if (i>=0) {
			i++;
			while (i<=PeriodicSubstrings.lastInterval && PeriodicSubstrings.intervals[i].firstPosition==position) i++;
		}
		else i=-i-1;
		j=i;
		inShort=false; inLong=false; farFromFirst=true; farFromLast=true;
		i=j-1;
		while (i>=0) {
			interval=PeriodicSubstrings.intervals[i];
			if (interval.lastPosition<position-IDENTITY_THRESHOLD) {
				i--;
				continue;
			}
			if (!interval.hasLongPeriod) inShort=true;
			else {
				inLong=true;
				if (position<=PeriodicSubstrings.intervals[i].firstPosition+IDENTITY_THRESHOLD) farFromFirst=false;
				if (position>=PeriodicSubstrings.intervals[i].lastPosition-IDENTITY_THRESHOLD) farFromLast=false;
			}
			i--;
		}
		i=j;
		while (i<=PeriodicSubstrings.lastInterval) {
			interval=PeriodicSubstrings.intervals[i];
			if (interval.firstPosition>position+IDENTITY_THRESHOLD) break;
			if (!interval.hasLongPeriod) inShort=true;
			else {
				inLong=true;
				farFromFirst=false;
				if (position>=PeriodicSubstrings.intervals[i].lastPosition-IDENTITY_THRESHOLD) farFromLast=false;
			}
			i++;
		}
		return (inShort?1:0)|(inLong?2:0)|(farFromFirst?4:0)|(farFromLast?8:0);
	}
	
	
	/**
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ to be sorted by
	 * $startA$.
	 *
	 * @param tmpSubstring temporary space;
	 * @return mask of flags telling whether the first (bit 0 from LSB) and the last (bit 
	 * 1 from LSB) position in readA of $alignment$ should be used to create an event, 
	 * based on their containment in dense substrings of substring type.
	 */
	private static final int usePositions_substring(Alignment alignment, boolean useUnknownEvents, DenseSubstring tmpSubstring) {
		final int MIN_MAXIMAL_ALIGNMENTS = IO.coverage*10;  // Arbitrary
		
		final boolean inSubstringStart = inSubstring(alignment.startA(),tmpSubstring);
		final boolean inSubstringEnd = inSubstring(alignment.endA(),tmpSubstring);
		final boolean useStart = (useUnknownEvents||alignment.isLeftMaximal==1) &&
			                     ( !inSubstringStart ||
								   (!inSubstringEnd && alignment.mergedToInterval!=null && alignment.mergedToInterval.nMaximalAlignmentsLeft>=MIN_MAXIMAL_ALIGNMENTS)
								 );
		final boolean useEnd = (useUnknownEvents||alignment.isRightMaximal==1) &&
					           ( !inSubstringEnd ||
							     (!inSubstringStart && alignment.mergedToInterval!=null && alignment.mergedToInterval.nMaximalAlignmentsRight>=MIN_MAXIMAL_ALIGNMENTS)
							   );
		return (useStart?1:0)|(useEnd?2:0);
	}
	
	
	/**
	 * Remark: the procedure assumes $DenseSubstrings.substrings$ to be sorted by
	 * $startA$.
	 *
	 * @param tmpSubstring temporary space;
	 * @return TRUE iff $position$ is contained in a dense substring of substring type, 
	 * and it is far enough from its boundaries.
	 */
	private static final boolean inSubstring(int position, DenseSubstring tmpSubstring) {
		final int BOUNDARY_DISTANCE = IO.quantum;
		int i, j;
		DenseSubstring substring;
		if (DenseSubstrings.lastSubstring==-1) return false;
		
		tmpSubstring.startA=position;
		i=Arrays.binarySearch(DenseSubstrings.substrings,0,DenseSubstrings.lastSubstring+1,tmpSubstring);
		if (i>=0) {
			i++;
			while (i<=DenseSubstrings.lastSubstring && DenseSubstrings.substrings[i].startA==position) i++;
		}
		else i=-i-1;
		j=i;
		i=j-1;
		while (i>=0) {
			substring=DenseSubstrings.substrings[i];
			if (substring.endA<position-BOUNDARY_DISTANCE) {
				i--;
				continue;
			}
			if (substring.substringReplication) return true;
			i--;
		}
		i=j;
		while (i<=DenseSubstrings.lastSubstring) {
			substring=DenseSubstrings.substrings[i];
			if (substring.startA>position+BOUNDARY_DISTANCE) break;
			if (substring.substringReplication) return true;
			i++;
		}
		return false;
	}
	

}