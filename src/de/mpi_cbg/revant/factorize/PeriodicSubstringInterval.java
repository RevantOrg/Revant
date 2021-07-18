package de.mpi_cbg.revant.factorize;

import java.io.IOException;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class PeriodicSubstringInterval implements Comparable {
	/**
	 * Sorting status of all intervals
	 */
	public static final int UNSORTED = -1;
	public static final int FIRSTPOSITION = 0;
	public static final int ID = 1;
	public static final int LASTPOSITION = 2;
	public static int order = UNSORTED;
	public static int order_longPeriod = UNSORTED;
	public static boolean sortByID = true;  // Enables final sorting by ID in $compareTo$.

	/**
	 * Properties of an interval
	 */
	public int id;
	public int firstPosition, lastPosition;
	public int nFactors;  // Number of factors that have a significant intersection with the interval
	public int firstFactor, lastFactor;
	public boolean discarded;
	public int nMergedIntervals;  // Number of intervals of which this interval is the merge
	public int mass;  // Number of maximal events in this interval
	public int period;
	public Point[] shifts;
	public int lastShift;
	public boolean isLeftMaximal, isRightMaximal;
	public boolean hasLongPeriod;
	public int nImpliedAlignments;
	public boolean equalsOnePeriod;
	public boolean mergedToDense;
	public int longestAlignment;
	public double errorRate;
	public int nAssignedAlignments;
	public double avgCoverage;
	
	/**
	 * Intervals
	 */
	public int isContained;  // Bitmask, where bit $Factorize.INTERVAL_*$ marks whether the current interval is contained in a different interval of that type.
	
	/**
	 * Temporary space
	 */
	public PeriodicSubstringInterval representative;
	public int surface;
	public boolean canBeTransformed, isContainedInShort, intersectsShort, flag1, flag2;
	public int nStartA, sumStartA, nEndA, sumEndA;
	public int tmpInt1, tmpInt2, tmpInt3, tmpInt4;
	
	
	/**
	 * @param maxShifts upper bound on the number of shifts to estimate the period.
	 */
	public final void allocateMemory(int maxShifts) {
		shifts = new Point[maxShifts];
		for (int i=0; i<shifts.length; i++) shifts[i] = new Point();
		lastShift=-1;
		representative=null;
	}
	
	
	public final void deallocate() {
		for (int i=0; i<shifts.length; i++) shifts[i]=null;
		shifts=null;
	}
	
	
	public final void clone(PeriodicSubstringInterval other) {
		id=other.id;
		firstPosition=other.firstPosition;
		lastPosition=other.lastPosition;
		nFactors=other.nFactors;
		firstFactor=other.firstFactor;
		lastFactor=other.lastFactor;
		discarded=other.discarded;
		nMergedIntervals=other.nMergedIntervals;
		mass=other.mass;
		period=other.period;
		if (shifts==null || shifts.length<other.lastShift+1) {
			shifts = new Point[other.lastShift+1];
			for (int i=0; i<shifts.length; i++) shifts[i] = new Point();
		}
		Points.simpleClone(other.shifts,other.lastShift,shifts);
		lastShift=other.lastShift;
		isLeftMaximal=other.isLeftMaximal;
		isRightMaximal=other.isRightMaximal;
		hasLongPeriod=other.hasLongPeriod;
		nImpliedAlignments=other.nImpliedAlignments;
		isContained=other.isContained;
		representative=other.representative;  // Copying the pointer
		surface=other.surface;
		canBeTransformed=other.canBeTransformed;
		equalsOnePeriod=other.equalsOnePeriod;
		mergedToDense=other.mergedToDense;
		longestAlignment=other.longestAlignment;
		errorRate=other.errorRate;
		nAssignedAlignments=other.nAssignedAlignments;
		avgCoverage=other.avgCoverage;
	}

	
	/**
	 * Adds to $this$ just a few properties of $other$.
	 */
	public final void simpleMerge(PeriodicSubstringInterval other) {
		nMergedIntervals+=other.nMergedIntervals;
		mass+=other.mass;
		nImpliedAlignments+=other.nImpliedAlignments;
		equalsOnePeriod&=other.equalsOnePeriod;
		mergedToDense|=other.mergedToDense;
		longestAlignment=Math.max(longestAlignment,other.longestAlignment);
		errorRate+=other.errorRate;
		nAssignedAlignments+=other.nAssignedAlignments;
		avgCoverage=(avgCoverage*length()+other.avgCoverage*other.length())/((length()+other.length())/2);
	}
	

	public final int compareTo(Object other) {
		PeriodicSubstringInterval otherInterval = (PeriodicSubstringInterval)other;
		int ord;

		if (!hasLongPeriod) ord=order;
		else ord=order_longPeriod;
		
		if (ord==FIRSTPOSITION) {
			if (firstPosition<otherInterval.firstPosition) return -1;
			else if (firstPosition>otherInterval.firstPosition) return 1;
		}
		else if (ord==ID) {
			if (id<otherInterval.id) return -1;
			else if (id>otherInterval.id) return 1;
			return 0;
		}
		else if (ord==LASTPOSITION) {
			if (lastPosition<otherInterval.lastPosition) return -1;
			else if (lastPosition>otherInterval.lastPosition) return 1;
		}
		// Final sorting by $id$, to make sure we can exactly reproduce an order after
		// having sorted by a different criterion.
		if (sortByID) {
			if (id<otherInterval.id) return -1;
			else if (id>otherInterval.id) return 1;
		}
		return 0;
	}
	
	
	public boolean equals(Object other) {
		PeriodicSubstringInterval otherInterval = (PeriodicSubstringInterval)other;
		return firstPosition==otherInterval.firstPosition && 
			   lastPosition==otherInterval.lastPosition;
	}


	public final void writePeriodicSubstringsFile(BufferedWriter file) throws IOException {
		file.write(id+"|");
		file.write(firstPosition+",");
		file.write(lastPosition+",");
		file.write(period+"|");
		file.write(isLeftMaximal?"1,":"0,");
		file.write(isRightMaximal?"1,":"0,");
		file.write(isContained+",");
		file.write(hasLongPeriod?"1,":"0,");
		file.write(equalsOnePeriod?"1":"0");
		file.write("|");
		file.write(nAssignedAlignments+",");
		file.write(IO.format(avgCoverage)+"\n");
	}
	
	
	/**
	 * Loads the numbers in $str[start..]$ into $out[from..]$, which must have at least 11
	 * cells.
	 */
    public static final void readPeriodicSubstringsFile(String str, int start, int[] out, int from) {
		int i, p, q;		

		i=from;
		p=str.indexOf("|",start);
		out[i++]=Integer.parseInt(str.substring(start,p));
		q=str.indexOf(",",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		out[i++]=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf("|",p+1);
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
	 * 5=isContained; 6=period; 7=hasLongPeriod; 8=equalsOnePeriod; 9=nAssignedAlignments;
	 * 10=avgCoverage.
	 * @param tokens temporary space, of size at least 11.
	 */
	public static final int loadPeriodicIntervals(int read, BufferedReader br, int readAheadLimit, int[][] out, int[] tokens) throws IOException {
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
			readPeriodicSubstringsFile(str,i+1,tokens,0);
			lastInterval++;
			out[lastInterval][0]=tokens[0];
			out[lastInterval][1]=tokens[1];
			out[lastInterval][2]=tokens[2];
			out[lastInterval][3]=tokens[4];
			out[lastInterval][4]=tokens[5];
			out[lastInterval][5]=tokens[6];
			out[lastInterval][6]=tokens[3];
			out[lastInterval][7]=tokens[7];
			out[lastInterval][8]=tokens[8];
			out[lastInterval][9]=tokens[9];
			out[lastInterval][10]=tokens[10];
			str=br.readLine();
		}
		return lastInterval;
	}
	
	
	/**
	 * Loads the intervals of a set of reads. The content of $out$ is identical to the
	 * other $loadPeriodicIntervals$ procedure. Read IDs are stored in $outReads$ instead.
	 *
	 * @param reads[0..lastRead] sorted list of distinct read IDs;
	 * @param tmp temporary space, with a number of rows at least equal to the maximum 
	 * number of intervals in a read, and with at least 11 columns;
	 * @param tokens temporary space, of size at least 11;
	 * @return the last interval loaded.
	 */
	public static final int loadPeriodicIntervals(int[] reads, int lastRead, BufferedReader br, int readAheadLimit, int[] outRead, int[][] out, int[][] tmp, int[] tokens) throws IOException {
		int i, j, k;
		int read, lastInterval;
		
		k=-1;
		for (i=0; i<=lastRead; i++) {
			read=reads[i];
			lastInterval=loadPeriodicIntervals(read,br,readAheadLimit,tmp,tokens);
			for (j=0; j<=lastInterval; j++) {
				k++;
				outRead[k]=read;
				System.arraycopy(tmp[j],0,out[k],0,11);
			}
		}
		return k;
	}
	
	
	/**
	 * A version of $loadPeriodicIntervals$ that just counts.
	 */
	public static final int countPeriodicIntervals(int[] reads, int lastRead, BufferedReader br, int readAheadLimit, int[][] tmp, int[] tokens) throws IOException {
		int i, out;
		
		out=0;
		for (i=0; i<=lastRead; i++) out+=1+loadPeriodicIntervals(reads[i],br,readAheadLimit,tmp,tokens);
		return out;
	}
	
	
	public String toString() {
		return hashCode()+" :: interval "+id+"=["+firstPosition+".."+lastPosition+"] nFactors="+nFactors+" discarded="+discarded+" mass="+mass+" period="+period+" hasLongPeriod="+hasLongPeriod+" isContained="+isContained+" nImpliedAlignments="+nImpliedAlignments+" --REPR--> "+(representative==null?"null":representative.hashCode())+" longestAlignment="+longestAlignment;
	}
	
	
	public int length() {
		return lastPosition-firstPosition+1;
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