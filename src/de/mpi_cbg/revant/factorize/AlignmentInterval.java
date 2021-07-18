package de.mpi_cbg.revant.factorize;

import java.io.IOException;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class AlignmentInterval implements Comparable {
	/**
	 * Parameters of the pipeline
	 */
	public static final int MAX_FACTORS_PER_INTERVAL = 200;
	public static final int FACTORS_GROWTH_RATE = 10;

	/**
	 * Sorting status of all intervals
	 */
	public static final int UNSORTED = -1;
	public static final int FIRSTPOSITION = 0;
	public static final int FACTORS = 1;
	public static final int LASTPOSITION = 2;
	public static int order = UNSORTED;
	public static boolean sortByID = true;  // Enables final sorting by ID in $compareTo$.

	/**
	 * Properties of an interval
	 */
	public int id;
	public int firstPosition, lastPosition;
	public Factor isAnchoredLeft, isAnchoredRight, containerFactor;
	public int[] factors;  // IDs of all factors that have a significant intersection with the interval
	public int lastFactor;
	public boolean discarded;
	public int nMergedIntervals;  // Number of alignments merged to this interval
	public int nMaximalAlignmentsLeft, nMaximalAlignmentsRight;  // Number of maximal alignments assigned to the interval
	public DenseSubstring inDenseSubstringOfSubstringType;
	public boolean isLongRepeat;  // TRUE iff the interval is likely to belong to a simple repeat that is longer than the max. length of a read.
	public boolean isLeftMaximal, isRightMaximal;
	public boolean hasSubstringWithLowCoverage, hasSubstringWithLowQuality;
	public double errorRate;
	public int nAssignedAlignments;
	public double avgCoverage;

	/**
	 * Intervals
	 */
	public PeriodicSubstringInterval inPeriodicSubstringInterval;  // Used by $markIntervalsInPeriodic()$ and $discardIntervals_reassignAlignments()$ in $AlignmentIntervals$, and by $assignAlignmentIntervals()$ in $Factors$.
	public PeriodicSubstringInterval periodicSubstringInterval;  // Used by $cleanAlignmentIntervals()$ and $discardIntervals_reassignAlignments()$ in $AlignmentIntervals$.
	public DenseSubstring denseSubstring;
	public int isContained;  // Bitmask, where bit $Factorize.INTERVAL_*$ marks whether the current interval is contained in a different interval of that type.
	public AlignmentInterval representative;
	public boolean inPeriodicRange;
	
	/**
	 * Temporary space
	 */
	public int tmpInt1;
	public boolean flag1, flag2;

	
	public AlignmentInterval() {
		factors = new int[MAX_FACTORS_PER_INTERVAL];
		cleanFactors();
	}


	public final void cleanFactors() {
		Math.set(factors,factors.length-1,-1);
		lastFactor=-1;
	}


	public final int compareTo(Object other) {
		int sameFactors;
		AlignmentInterval otherInterval = (AlignmentInterval)other;

		if (order==FIRSTPOSITION) {
			if (firstPosition<otherInterval.firstPosition) return -1;
			else if (firstPosition>otherInterval.firstPosition) return 1;
		}
		else if (order==LASTPOSITION) {
			if (lastPosition<otherInterval.lastPosition) return -1;
			else if (lastPosition>otherInterval.lastPosition) return 1;
		}
		else if (order==FACTORS) {
			sameFactors=sameFactorsAs(otherInterval);
			if (sameFactors!=0) return sameFactors;
		}
		// Final sorting by $id$, to make sure we can exactly reproduce an order after
		// having sorted by a different criterion.
		if (sortByID) {
			if (id<otherInterval.id) return -1;
			else if (id>otherInterval.id) return 1;
		}
		return 0;
	}
	
	
	/**
	 * @return same convention as $compareTo$. 0=same factors; 1: $this$ is lex. bigger 
	 * than $otherInterval$; -1: $this$ is lex. smaller than $otherInterval$.
	 */
	public final int sameFactorsAs(AlignmentInterval otherInterval) {
		int i, max;
		
		max=lastFactor;
		if (otherInterval.lastFactor<max) max=otherInterval.lastFactor;
		for (i=0; i<=max; i++) {
			if (factors[i]<otherInterval.factors[i]) return -1;
			if (factors[i]>otherInterval.factors[i]) return 1;
		}
		if (lastFactor<otherInterval.lastFactor) return -1;
		else if (lastFactor>otherInterval.lastFactor) return 1;
		return 0;
	}


	public final void writeAlignmentsFile(BufferedWriter file) throws IOException {
		file.write(id+",");
		file.write(firstPosition+",");
		file.write(lastPosition+",");
		file.write(isLongRepeat?"1,":"0,");
		file.write(isLeftMaximal?"1,":"0,");
		file.write(isRightMaximal?"1,":"0,");
		file.write(isContained+"|");
		file.write(nAssignedAlignments+",");
		file.write(IO.format(avgCoverage)+"\n");
	}
	
	
	/**
	 * Loads the numbers in $str[start..]$ into $out[from..]$, which must have at least 10
	 * cells.
	 */
    public static final void readAlignmentsFile(String str, int start, int[] out, int from) {
		int i, p, q;
		
		i=from;
		p=str.indexOf(",",start);
		out[i++]=Integer.parseInt(str.substring(start,p));
		q=str.indexOf(",",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf("|",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q;
		out[i++]=(int)Double.parseDouble(str.substring(p+1));
	}
	
	
	/**
	 * Loads in $out[0..X]$ the set of all intervals in file $br$ that belong to $read$, 
	 * where $X$ is returned in output. Every line in the file is assumed to be prefixed 
	 * by a read ID. The procedure starts reading from the current position of $br$.
	 *
	 * @param out columns: 0=id, 1=start, 2=end; 3=isLeftMaximal; 4=isRightMaximal;
	 * 5=isContained; 6=nAssignedAlignments; 7=avgCoverage;
	 * @param tokens temporary space, of size at least 10.
	 */
	public static final int loadAlignmentIntervals(int read, BufferedReader br, int readAheadLimit, int[][] out, int[] tokens) throws IOException {
		int i, r, lastInterval;
		String str;
		
		lastInterval=-1;
		str=br.readLine();
		while (str!=null) {
			i=str.indexOf(",");
			r=Integer.parseInt(str.substring(0,i));
			if (r<read) {
				str=br.readLine();
				continue;
			}
			else if (r>read) {
				br.reset();
				return lastInterval;
			}
			br.mark(readAheadLimit);
			readAlignmentsFile(str,i+1,tokens,0);
			lastInterval++;
			out[lastInterval][0]=tokens[0];
			out[lastInterval][1]=tokens[1];
			out[lastInterval][2]=tokens[2];
			out[lastInterval][3]=tokens[5];
			out[lastInterval][4]=tokens[6];
			out[lastInterval][5]=tokens[7];
			out[lastInterval][6]=tokens[8];
			out[lastInterval][7]=tokens[9];
			str=br.readLine();
		}
		return lastInterval;
	}
	
	
	/**
	 * Loads the intervals of a set of reads. The content of $out$ is identical to the
	 * other $loadAlignmentIntervals$ procedure. Read IDs are stored in $outReads$ instead.
	 *
	 * @param reads[0..lastRead] sorted list of distinct read IDs;
	 * @param tmp temporary space, with a number of rows at least equal to the maximum 
	 * number of intervals in a read, and with at least 8 columns;
	 * @param tokens temporary space, of size at least 10;
	 * @return the last interval loaded.
	 */
	public static final int loadAlignmentIntervals(int[] reads, int lastRead, BufferedReader br, int readAheadLimit, int[] outRead, int[][] out, int[][] tmp, int[] tokens) throws IOException {
		int i, j, k;
		int read, lastInterval;
		
		k=-1;
		for (i=0; i<=lastRead; i++) {
			read=reads[i];
			lastInterval=loadAlignmentIntervals(read,br,readAheadLimit,tmp,tokens);
			for (j=0; j<=lastInterval; j++) {
				k++;
				outRead[k]=read;
				System.arraycopy(tmp[j],0,out[k],0,6);
			}
		}
		return k;
	}
	
	
	/**
	 * A version of $loadAlignmentIntervals$ that just counts.
	 */
	public static final int countAlignmentIntervals(int[] reads, int lastRead, BufferedReader br, int readAheadLimit, int[][] tmp, int[] tokens) throws IOException {
		int i, out;
		
		out=0;
		for (i=0; i<=lastRead; i++) out+=1+loadAlignmentIntervals(reads[i],br,readAheadLimit,tmp,tokens);
		return out;
	}


	public String toString() {
		return hashCode()+" id="+id+" ["+firstPosition+".."+lastPosition+"] isLeftMaximal="+isLeftMaximal+" isRightMaximal="+isRightMaximal+" nMergedIntervals="+nMergedIntervals+" nMaximalAlignmentsLeft="+nMaximalAlignmentsLeft+" nMaximalAlignmentsRight="+nMaximalAlignmentsRight+" periodicSubstringInterval="+(periodicSubstringInterval==null?"null":periodicSubstringInterval.hashCode())+" inPeriodicSubstringInterval="+(inPeriodicSubstringInterval==null?"null":inPeriodicSubstringInterval.hashCode())+" inDenseSubstringOfSubstringType="+inDenseSubstringOfSubstringType+" isLongRepeat="+isLongRepeat+" discarded="+discarded;
	}
	
	
	/**
	 * Remark: the procedure does not assign an alignment $[x..y]$ to an interval, if $x$ 
	 * is far from the beginning of the interval and readA has low quality before $x$
	 * (a symmetrical argument applies to $y$).
	 *
	 * @return TRUE iff $alignment$ can be assigned to this alignment interval.
	 */
	public final boolean canBeAssigned(Alignment alignment, int threshold) {
		final int alignmentID = alignment.id;
		final int startA = Alignments.alignments[alignmentID][3];
		final int endA = Alignments.alignments[alignmentID][4];
		
		if ( (alignment.isLeftMaximal==1 && startA>firstPosition+threshold) ||
			 (alignment.isRightMaximal==1 && endA<lastPosition-threshold) ||
			 (isLeftMaximal && firstPosition>startA+threshold) ||
			 (isRightMaximal && lastPosition<endA-threshold) ||
			 (startA>firstPosition+threshold && !Reads.isLeftMaximal(startA,ReadA.id,true)) ||
			 (endA<lastPosition-threshold && !Reads.isRightMaximal(endA,ReadA.id,true))
		   ) return false;
		return true;
	}
	
	
	/**
	 * Tries to remove low-quality regions from the start/end of the interval.
	 * Trimming is not performed if it would create an interval of length smaller than
	 * $minLength$.
	 *
	 * Remark: $isAnchored{Left,Right}$, $nMaximalAlignments{Left,Right}$, $in*$, 
	 * $is{Left,Right}Maximal$, $denseSubstring$, $isContained$, and possibly other 
	 * fields, might get invalidated by the procedure.
	 */
	public final void trim(int[] tmpArray, int minLength, boolean trimLeft, boolean trimRight) {
		final int WINDOW_LARGE = 3;  // Arbitrary
		final int WINDOW_SMALL = 1;
		
		tmpArray[0]=firstPosition; tmpArray[1]=lastPosition;
		Reads.trim(tmpArray,WINDOW_LARGE,trimLeft,trimRight);
		if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<minLength) return;
		Reads.trim(tmpArray,WINDOW_SMALL,trimLeft,trimRight);
		if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<minLength) return;
		if (tmpArray[0]!=firstPosition) {
			firstPosition=tmpArray[0];
			isLeftMaximal=false;
		}
		if (tmpArray[1]!=lastPosition) {
			lastPosition=tmpArray[1];
			isRightMaximal=false;
		}
	}
	
	
	public final int length() {
		return lastPosition-firstPosition+1;
	}
	
	
	public final int nMaximalAlignments() {
		return (nMaximalAlignmentsLeft+nMaximalAlignmentsRight)>>1;
	}
	
	
	/**
	 * Tries to transform the interval into a dense substring, using the B-maximal 
	 * alignments in $alignments$.
	 *
	 * Remark: the resulting substring is marked as weak.
	 *
	 * @param alignmentStats output of $Factorize.generalizeReplicationTypes_
	 * getAlignmentStats()$;
	 * @param distanceThreshold threshold on distance identity;
	 * @param alignmentsThreshold minimum number of B-maximal alignments to transform the 
	 * interval;
	 * @param surfaceThreshold the alignments in $alignments$ must cover at least this
	 * fraction of the length of the interval for it to be transformed;
	 * @param tmpArray temporary space, with at least 3 cells;
	 * @return TRUE iff the interval was transformed into a dense substring, whose 
	 * representation is stored in $tmpSubstring$.
	 */
	public final boolean generalizeType(int[] alignmentStats, int nAlignments, int distanceThreshold, int alignmentsThreshold, double surfaceThreshold, DenseSubstring tmpSubstring, int[] tmpArray) {
		boolean transformed;
		if (nAlignments<alignmentsThreshold) return false;
		final int nLeftMaximal = alignmentStats[0];
		final int nRightMaximal = alignmentStats[1];
		final int nLeftRightMaximal = alignmentStats[2];
		final boolean leftMaximalAligned = alignmentStats[3]==1;
		final boolean rightMaximalAligned = alignmentStats[4]==1;
		final int minPrefLength = alignmentStats[5];
		final int minSufLength = alignmentStats[6];
		final int minSubstringLength = alignmentStats[7];
		final int firstMaxStart = alignmentStats[8];
		final int lastMaxEnd = alignmentStats[9];
		final int surface = alignmentStats[10];
		int length;
		if (surface<(lastPosition-firstPosition+1)*surfaceThreshold) return false;

		// Generalizing replication type
		transformed=false;
		if ( nLeftRightMaximal>=alignmentsThreshold ||
			 (nLeftMaximal+nRightMaximal>=alignmentsThreshold) ||
			 (nLeftMaximal>=alignmentsThreshold && !leftMaximalAligned) ||
		     (nRightMaximal>=alignmentsThreshold && !rightMaximalAligned)
		   ) {
			tmpSubstring.prefixReplication=false;
			tmpSubstring.suffixReplication=false;
			tmpSubstring.singleDeletionReplication=false;
			tmpSubstring.substringReplication=true;
			tmpArray[0]=minPrefLength;
			tmpArray[1]=minSufLength;
			tmpArray[2]=minSubstringLength;
			length=Math.minPositive(tmpArray,2);
			tmpSubstring.minPrefixLength=length!=-1?length:0;
			transformed=true;
		}
		else if (nLeftMaximal>=alignmentsThreshold && nRightMaximal>=alignmentsThreshold && leftMaximalAligned && rightMaximalAligned) {
			tmpSubstring.prefixReplication=true;
			tmpSubstring.suffixReplication=true;
			tmpSubstring.singleDeletionReplication=false;  // We don't try to estimate it.
			tmpSubstring.substringReplication=false;
			tmpSubstring.minPrefixLength=minPrefLength;
			tmpSubstring.minSuffixLength=minSufLength;
			transformed=true;
		}
		else if (nLeftMaximal>=alignmentsThreshold && leftMaximalAligned) {
			tmpSubstring.prefixReplication=false;
			tmpSubstring.suffixReplication=true;
			tmpSubstring.singleDeletionReplication=false;
			tmpSubstring.substringReplication=false;
			tmpSubstring.minSuffixLength=minSufLength;
			transformed=true;
		}
		else if (nRightMaximal>=alignmentsThreshold && rightMaximalAligned) {
			tmpSubstring.prefixReplication=true;
			tmpSubstring.suffixReplication=false;
			tmpSubstring.singleDeletionReplication=false;
			tmpSubstring.substringReplication=false;
			tmpSubstring.minPrefixLength=minPrefLength;
			transformed=true;
		}
		if (!transformed) return false;
		tmpSubstring.startA=firstPosition;
		tmpSubstring.endA=lastPosition;
		tmpSubstring.nStartA=1;
		tmpSubstring.sumStartA=tmpSubstring.startA;
		tmpSubstring.nEndA=1;
		tmpSubstring.sumEndA=tmpSubstring.endA;
		tmpSubstring.isLeftMaximal=isLeftMaximal;
		tmpSubstring.isRightMaximal=isRightMaximal;
		tmpSubstring.firstMaximalStart=firstMaxStart;
		tmpSubstring.lastMaximalEnd=lastMaxEnd;
		tmpSubstring.nMaximalAlignments+=nAlignments;
		tmpSubstring.isContained=isContained;
		tmpSubstring.isWeak=true;
		return true;
	}
	
	
	/**
	 * @return TRUE iff the interval is maximal with at least $minMaxAlignments$ maximal 
	 * alignments on one side, or if the interval contains at least $minAlignments$ in 
	 * total.
	 */
	public final boolean checkAlignments(int minMaxAlignments, int minAlignments) {
		boolean out = (isLeftMaximal && nMaximalAlignmentsLeft>=minMaxAlignments) ||
		              (isRightMaximal && nMaximalAlignmentsRight>=minMaxAlignments) ||
		              nMergedIntervals>=minAlignments;
		return out;
	}
	
	
	/**
	 * Returns the last element in the chain of pointers that starts from
	 * $representative$. No chain of pointers should end with null, since the
	 * $representative$ field of every interval should have been initialized to the
	 * interval itself.
	 */
	public final AlignmentInterval closure() {
		AlignmentInterval out = this;
		while (out!=null && out.representative!=out) out=out.representative;
		if (out==null) {
			IO.printCriticalErr("ERROR: a chain of $representative$ pointers ends with null.");
			System.exit(1);
		}
		return out;
	}
	

}