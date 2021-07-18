package de.mpi_cbg.revant.factorize;

import java.io.IOException;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Constants;


public class DenseSubstring implements Comparable {
	/**
	 * Sorting status of all dense substrings
	 */
	public static final int UNSORTED = -1;
	public static final int STARTA = 0;
	public static final int PRSR_STARTA = 1;
	public static final int LENGTH = 2;
	public static final int ID = 3;
	public static final int ENDA = 4;
	public static final int ENDA_REVERSE = 5;
	public static int order = UNSORTED;
	public static boolean sortByID = true;  // Enables final sorting by ID in $compareTo$.

	/**
	 * Properties of a dense substring
	 */
	public int id;
	public int startA;  // Lower bound on the first position in $readA$
	public int endA;  // Upper bound on the last position in $readA$
	public int nStartA, sumStartA;
	public int nEndA, sumEndA;
	public boolean prefixReplication;
	public boolean suffixReplication;
	public boolean substringReplication;  // A repeat that replicates either by taking substrings of itself, or by applying more than one internal deletion.
	public boolean singleDeletionReplication;  // A repeat that replicates by applying one internal deletion, or by taking prefixes and suffixes that never come close to the first/last position.
// 	public int nFactors;  // Number of factors that have a significant intersection with the dense substring
	public DenseSubstring representative;
	public int nSingleDeletions, shortestSingleDeletion, longestSingleDeletion;  // -1 if $singleDeletionReplication=true$ but the values could not be estimated.
	public int minPrefixLength, minSuffixLength;  // Estimate of the minimum length of a prefix/suffix occurrence of the dense substring in a readB. If the dense substring is of substring type, $minPrefixLength$ contains the shortest length of a left- and right-maximal substring in a readB.
	public int pathLength;
	public int maximalStartA, maximalEndA;
	public int nImpliedAlignmentsPrefix, nImpliedAlignmentsSuffix;  // Only implied alignments that behave like prefixes/suffixes of a prefix/suffix substring are counted.
	public double minNAlignmentsLeft, minNAlignmentsRight;  // Threshold in the number of maximal alignments, not implied by dense or periodic substrings, used to decide whether to keep the boundary of a factor or not. Defined only when $substringReplication=true$.
	public int firstFactor, lastFactor;
	public DenseSubstring nextSuffixSubstring, previousPrefixSubstring;
	public boolean prefixesAreUniformlyDistributed, suffixesAreUniformlyDistributed;
	public int minPrefixDistanceBias, maxPrefixDistanceBias;
	public int minSuffixDistanceBias, maxSuffixDistanceBias;
	public boolean isLeftMaximal, isRightMaximal;
	public int period;
	public boolean strongPeriodSignal;
	public int firstMaximalStart;  // Absolute position in the read, excluding the beginning of the substring.
	public int lastMaximalEnd;  // Absolute position in the read, excluding the end of the substring.
	public int nMaximalAlignments;  // Number of maximal alignments that suggest distinct occurrence of this dense substring in the genome.
	public boolean hasLengthBoundariesPrefix, hasLengthBoundariesSuffix;
	public int nAssignedAlignments;
	public double avgCoverage;
	
	/**
	 * Intervals
	 */
	public PeriodicSubstringInterval periodicSubstringInterval;
	public int isContained;  // Bitmask, where bit $Constants.INTERVAL_*$ marks whether the current interval is contained in a different interval of that type.
	
	/**
	 * Density windows
	 */
	public int leftWindow, rightWindow, maxWindow;
	public int surface;
	
	/**
	 * Temporary space
	 */
	public int firstAlignment;
	public boolean isMerged, isConcatenated;
	public boolean startAAdded, endAAdded;
	public int destinationAlignment, leftAlignment, rightAlignment;  // Temporary variables used by $DenseSubstrings.markImpliedAlignments*$ procedures
	public boolean isWeak;
	public boolean discarded;
	public int previousSumStartA, previousSumEndA, previousNStartA, previousNEndA;


	/**
	 * Copies all fields of $from$ onto this substring.
	 */
	public final void clone(DenseSubstring from) {
		id=from.id;
		startA=from.startA;
		endA=from.endA;
		nStartA=from.nStartA;
		sumStartA=from.sumStartA;
		nEndA=from.nEndA;
		sumEndA=from.sumEndA;
		prefixReplication=from.prefixReplication;
		suffixReplication=from.suffixReplication;
		substringReplication=from.substringReplication;
		singleDeletionReplication=from.singleDeletionReplication;
		representative=from.representative;
		nSingleDeletions=from.nSingleDeletions;
		shortestSingleDeletion=from.shortestSingleDeletion;
		longestSingleDeletion=from.longestSingleDeletion;
		minPrefixLength=from.minPrefixLength;
		minSuffixLength=from.minSuffixLength;
		pathLength=from.pathLength;
		maximalStartA=from.maximalStartA;
		maximalEndA=from.maximalEndA;
		nImpliedAlignmentsPrefix=from.nImpliedAlignmentsPrefix;
		nImpliedAlignmentsSuffix=from.nImpliedAlignmentsSuffix;
		minNAlignmentsLeft=from.minNAlignmentsLeft;
		minNAlignmentsRight=from.minNAlignmentsRight;
		firstFactor=from.firstFactor;
		lastFactor=from.lastFactor;
		nextSuffixSubstring=from.nextSuffixSubstring;
		previousPrefixSubstring=from.previousPrefixSubstring;
		prefixesAreUniformlyDistributed=from.prefixesAreUniformlyDistributed;
		suffixesAreUniformlyDistributed=from.suffixesAreUniformlyDistributed;
		minPrefixDistanceBias=from.minPrefixDistanceBias;
		maxPrefixDistanceBias=from.maxPrefixDistanceBias;
		minSuffixDistanceBias=from.minSuffixDistanceBias;
		maxSuffixDistanceBias=from.maxSuffixDistanceBias;
		isLeftMaximal=from.isLeftMaximal;
		isRightMaximal=from.isRightMaximal;
		period=from.period;
		strongPeriodSignal=from.strongPeriodSignal;
		firstMaximalStart=from.firstMaximalStart;
		lastMaximalEnd=from.lastMaximalEnd;
		nMaximalAlignments=from.nMaximalAlignments;
		periodicSubstringInterval=from.periodicSubstringInterval;
		isContained=from.isContained;
		leftWindow=from.leftWindow;
		rightWindow=from.rightWindow;
		maxWindow=from.maxWindow;
		hasLengthBoundariesPrefix=from.hasLengthBoundariesPrefix;
		hasLengthBoundariesSuffix=from.hasLengthBoundariesSuffix;
		nAssignedAlignments=from.nAssignedAlignments;
		avgCoverage=from.avgCoverage;
	}


	/**
	 * Returns the last element in the chain of pointers that starts from
	 * $representative$. No chain of pointers should end with null, since the
	 * $representative$ field of every substring should have been initialized to the
	 * substring itself.
	 */
	public final DenseSubstring closure() {
		DenseSubstring out = this;
		while (out!=null && out.representative!=out) out=out.representative;
		if (out==null) {
			IO.printCriticalErr("ERROR: a chain of $representative$ pointers ends with null.");
			System.exit(1);
		}
		return out;
	}


	/**
	 * @return TRUE iff the substring differs in at least one replication mode from
	 * substring $other$.
	 */
	public final boolean differentReplication(DenseSubstring other) {
		return prefixReplication!=other.prefixReplication ||
			   suffixReplication!=other.suffixReplication ||
			   substringReplication!=other.substringReplication ||
			   singleDeletionReplication!=other.singleDeletionReplication;
	}


	public final void clearImpliedAlignmentsStats() {
		nImpliedAlignmentsPrefix=0;
		nImpliedAlignmentsSuffix=0;
	}
	
	
	public final void clearPrefixSuffixBias() {
		prefixesAreUniformlyDistributed=true;
		minPrefixDistanceBias=-1;
		maxPrefixDistanceBias=-1;
		hasLengthBoundariesPrefix=false;
		suffixesAreUniformlyDistributed=true;
		minSuffixDistanceBias=-1;
		maxSuffixDistanceBias=-1;
		hasLengthBoundariesSuffix=false;
	}


    public final int compareTo(Object other) {
		DenseSubstring otherSubstring = (DenseSubstring)other;
    	int length, otherLength;
    	if (order==STARTA) {
			if (startA<otherSubstring.startA) return -1;
			else if (startA>otherSubstring.startA) return 1;
		}
		else if (order==PRSR_STARTA) {
			if ( (prefixReplication||suffixReplication) && !otherSubstring.prefixReplication && !otherSubstring.suffixReplication ) return -1;
			else if ( (otherSubstring.prefixReplication||otherSubstring.suffixReplication) && !prefixReplication && !suffixReplication ) return 1;
			if (startA<otherSubstring.startA) return -1;
			else if (startA>otherSubstring.startA) return 1;
		}
		else if (order==LENGTH) {
			length=endA-startA+1;
			otherLength=otherSubstring.endA-otherSubstring.startA+1;
			if (length<otherLength) return -1;
			else if (length>otherLength) return 1;
		}
		else if (order==ID) {
			if (id<otherSubstring.id) return -1;
			else if (id>otherSubstring.id) return 1;
		}
		else if (order==ENDA) {
			if (endA<otherSubstring.endA) return -1;
			else if (endA>otherSubstring.endA) return 1;
		}
		else if (order==ENDA_REVERSE) {
			if (endA<otherSubstring.endA) return 1;
			else if (endA>otherSubstring.endA) return -1;
		}
		// Final sorting by $id$, to make sure we can exactly reproduce an order after
		// having sorted by a different criterion.
		if (sortByID) {
			if (id<otherSubstring.id) return -1;
			else if (id>otherSubstring.id) return 1;
		}
		return 0;
    }
	

    public final void writeDenseSubstringsFile(BufferedWriter file) throws IOException {
		file.write(id+"|");
		file.write(isWeak?"1":"0");
		file.write("|");
		file.write(startA+",");
		file.write(endA+"|");
		if (prefixReplication) {
			file.write("1,");
			file.write(minPrefixLength+",");
		}
		else if (substringReplication) {
			file.write("0,");
			file.write(minPrefixLength+",");
		}
		else file.write("0,0,");
		if (suffixReplication) {
			file.write("1,");
			file.write(minSuffixLength+",");
		}
		else file.write("0,0,");
		file.write(substringReplication?"1,":"0,");
		file.write(period==-1?"0,":period+",");
		file.write(strongPeriodSignal?"1":"0");
		file.write("|");
		file.write(prefixesAreUniformlyDistributed?"0,":"1,");
		file.write(minPrefixDistanceBias+","+maxPrefixDistanceBias+",");
		file.write(suffixesAreUniformlyDistributed?"0,":"1,");
		file.write(minSuffixDistanceBias+","+maxSuffixDistanceBias);
		file.write("|");
		file.write(singleDeletionReplication?"1,":"0,");
		file.write(nSingleDeletions+",");
		file.write(shortestSingleDeletion+",");
		file.write(longestSingleDeletion+"|");
		file.write(isLeftMaximal?"1,":"0,");
		file.write(isRightMaximal?"1,":"0,");
		file.write(isContained+",");
		file.write(minPrefixLength+",");
		file.write(minSuffixLength+"|");
		file.write(firstMaximalStart+",");
		file.write(lastMaximalEnd+"|");
		file.write(nAssignedAlignments+",");
		file.write(IO.format(avgCoverage)+"\n");
	}
	
	
	/**
	 * Loads the numbers in $str[start..]$ into $out[from..]$, which must have at least 30
	 * cells.
	 */
    public static final void readDenseSubstringsFile(String str, int start, int[] out, int from) {
		int i, p, q;
		
		i=from;
		p=str.indexOf("|",start);
		out[i++]=Integer.parseInt(str.substring(start,p));
		q=str.indexOf("|",p+1);
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
		p=q; q=str.indexOf(",",p+1);
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
	 * @param out columns: 0=id, 1=type, 2=start, 3=end; 4=isLeftMaximal; 
	 * 5=isRightMaximal; 6=isWeak; 7=isContained; 8=minPrefixLength; 9=minSuffixLength;
	 * 10=firstMaximalStart; 11=lastMaximalEnd; 12=nAssignedAlignments; 13=avgCoverage;
	 * @param tokens temporary space, of size at least 30.
	 */
	public static final int loadDenseIntervals(int read, BufferedReader br, int readAheadLimit, int[][] out, int[] tokens) throws IOException {
		boolean prefixReplication, suffixReplication, substringReplication, singleDeletionReplication;
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
			readDenseSubstringsFile(str,i+1,tokens,0);
			lastInterval++;
			out[lastInterval][0]=tokens[0];
			out[lastInterval][1]=getType(tokens);
			out[lastInterval][2]=tokens[2];
			out[lastInterval][3]=tokens[3];
			out[lastInterval][4]=tokens[21];
			out[lastInterval][5]=tokens[22];
			out[lastInterval][6]=tokens[1];
			out[lastInterval][7]=tokens[23];
			out[lastInterval][8]=tokens[24];
			out[lastInterval][9]=tokens[25];
			out[lastInterval][10]=tokens[26];
			out[lastInterval][11]=tokens[27];
			out[lastInterval][12]=tokens[28];
			out[lastInterval][13]=tokens[29];
			str=br.readLine();
		}
		return lastInterval;
	}
	
	
	/**
	 * Loads the intervals of a set of reads. The content of $out$ is identical to the
	 * other $loadDenseIntervals$ procedure. Read IDs are stored in $outReads$ instead.
	 *
	 * @param reads[0..lastRead] sorted list of distinct read IDs;
	 * @param tmp temporary space, with a number of rows at least equal to the maximum 
	 * number of intervals in a read, and with at least 14 columns;
	 * @param tokens temporary space, of size at least 30;
	 * @return the last interval loaded.
	 */
	public static final int loadDenseIntervals(int[] reads, int lastRead, BufferedReader br, int readAheadLimit, int[] outRead, int[][] out, int[][] tmp, int[] tokens) throws IOException {
		int i, j, k;
		int read, lastInterval;
		
		k=-1;
		for (i=0; i<=lastRead; i++) {
			read=reads[i];
			lastInterval=loadDenseIntervals(read,br,readAheadLimit,tmp,tokens);
			for (j=0; j<=lastInterval; j++) {
				k++;
				outRead[k]=read;
				System.arraycopy(tmp[j],0,out[k],0,10);
			}
		}
		return k;
	}
	
	
	/**
	 * A version of $loadDenseIntervals$ that just counts.
	 */
	public static final int countDenseIntervals(int[] reads, int lastRead, BufferedReader br, int readAheadLimit, int[][] tmp, int[] tokens) throws IOException {
		int i, out;
		
		out=0;
		for (i=0; i<=lastRead; i++) out+=1+loadDenseIntervals(reads[i],br,readAheadLimit,tmp,tokens);
		return out;
	}
	
	
	/**
	 * @param tokens output of $readDenseSubstringsFile$.
	 */
	public static final int getType(int[] tokens) {
		boolean prefixReplication = tokens[4]==1;
		boolean suffixReplication = tokens[6]==1;
		boolean substringReplication = tokens[8]==1;
		boolean singleDeletionReplication = tokens[17]==1;
		if (prefixReplication && !suffixReplication && !substringReplication && !singleDeletionReplication) return Constants.INTERVAL_DENSE_PREFIX;
		else if (!prefixReplication && suffixReplication && !substringReplication && !singleDeletionReplication) return Constants.INTERVAL_DENSE_SUFFIX;
		else if (prefixReplication && suffixReplication && !substringReplication && !singleDeletionReplication) return Constants.INTERVAL_DENSE_PREFIXSUFFIX;
		else if (!prefixReplication && !suffixReplication && substringReplication && !singleDeletionReplication) return Constants.INTERVAL_DENSE_SUBSTRING;
		else if (singleDeletionReplication) return Constants.INTERVAL_DENSE_SINGLEDELETION;
		System.err.println("ERROR: interval of unidentified type");
		System.exit(1);
		return -1;
	}


	public String toString() {
		return hashCode()+"=["+startA+","+maximalStartA+".."+endA+","+maximalEndA+"] avgEndA="+(((double)sumEndA)/nEndA)+" id="+id+" pathLength="+pathLength+" isWeak="+isWeak+" representativeID="+(representative==null?"null":representative.id)+" prefixReplication="+prefixReplication+(prefixReplication?", minPrefixLength="+minPrefixLength:"")+", suffixReplication="+suffixReplication+(suffixReplication?", minSuffixLength="+minSuffixLength:"")+", substringReplication="+substringReplication+", singleDeletionReplication="+singleDeletionReplication+" ("+sumStartA+"/"+nStartA+"), ("+sumEndA+"/"+nEndA+"), prefixesAreUniformlyDistributed="+prefixesAreUniformlyDistributed+" suffixesAreUniformlyDistributed="+suffixesAreUniformlyDistributed+" discarded="+discarded+" representative.substringReplication="+(representative==null?"-":representative.substringReplication)+" hasLengthBoundariesPrefix="+hasLengthBoundariesPrefix+" hasLengthBoundariesSuffix="+hasLengthBoundariesSuffix+" nMaximalAlignments="+nMaximalAlignments+" leftAlignment="+leftAlignment+" rightAlignment="+rightAlignment+"\n";
	}


	/**
	 * Returns a lower bound on the smallest index of an alignment in
	 * $ReadA.sortedAlignments$ that could be implied by this dense substring.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments$ to be sorted by startA.
	 */
	public final int getFromAlignment(int identityThreshold) {
		int i;
		for (i=leftAlignment-1; i>=0; i--) {
			if (Alignments.alignments[ReadA.sortedAlignments[i].id][3]<startA-identityThreshold) break;
		}
		return i+1;
	}
	
	
	/**
	 * The procedure assumes $ReadA.sortedAlignments$ to be sorted by endA.
	 */
	public final int getFromAlignmentPrime(int identityThreshold, int lastAlignment) {
		int i;
		for (i=rightAlignment+1; i<=lastAlignment; i++) {
			if (Alignments.alignments[ReadA.sortedAlignments[i].id][4]>endA+identityThreshold) break;
		}
		return i-1;
	}


	/**
	 * Returns an upper bound on the largest index of an alignment in
	 * $ReadA.sortedAlignments$ that could be implied by this dense substring.
	 *
	 * Remark: the procedure assumes $ReadA.sortedAlignments$ to be sorted by startA.
	 */
	public final int getToAlignment(int identityThreshold, int lastAlignment) {
		int i;
		for (i=rightAlignment+1; i<=lastAlignment; i++) {
			if (Alignments.alignments[ReadA.sortedAlignments[i].id][3]>endA) break;
		}
		return i-1;
	}

	
	/**
	 * @return TRUE iff $alignment$ can be assigned to this dense substring.
	 */
	public final boolean canBeAssigned(Alignment alignment, int threshold) {
		if (substringReplication) return true;
		final int alignmentID = alignment.id;
		final int alignmentStart = Alignments.alignments[alignmentID][3];
		final int alignmentEnd = Alignments.alignments[alignmentID][4];
		if (prefixReplication && !suffixReplication && !singleDeletionReplication) {
			if (alignment.isLeftMaximal==1 && alignmentStart>startA+threshold) return false;
			// Stronger condition not used any more: || (alignment.isRightMaximal==1 && alignmentEnd<startA+minPrefixLength-threshold)
		}
		if (suffixReplication && !prefixReplication && !singleDeletionReplication) {
			if (alignment.isRightMaximal==1 && alignmentEnd<endA-threshold) return false;
			// Stronger condition not used any more: || (alignment.isLeftMaximal==1 && alignmentStart>endA-minSuffixLength+threshold)
		}
		if ((prefixReplication && suffixReplication) || singleDeletionReplication) {
			if (alignment.isLeftMaximal==1 && alignmentStart>startA+threshold && alignment.isRightMaximal==1 && alignmentEnd<endA-threshold) return false;
			// Stronger conditions not used, as above.
		}
		return true;
	}

	
	public int length() {
		return endA-startA+1;
	}
	
	
	/**
	 * Tries to remove low-quality regions from the start/end of the substring.
	 * Trimming is not performed if it would create an interval of length smaller than
	 * $minLength$.
	 *
	 * Remark: $maximal{Start,End}A$, $firstMaximalStart,lastMaximalEnd$, 
	 * $nMaximalAlignments$, $is{Left,Right}Maximal$ and $isContained$ might get 
	 * invalidated by the procedure.
	 */
	public final void trim(int[] tmpArray, int minLength) {
		final int WINDOW_LARGE = 3;  // Arbitrary
		final int WINDOW_SMALL = 1;
		int delta;
		
		tmpArray[0]=startA; tmpArray[1]=endA;
		Reads.trim(tmpArray,WINDOW_LARGE,true,true);
		if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<minLength) return;
		Reads.trim(tmpArray,WINDOW_SMALL,true,true);
		if (tmpArray[0]==-1 || tmpArray[1]==-1 || tmpArray[1]-tmpArray[0]+1<minLength) return;
		if (tmpArray[0]!=startA) {
			delta=tmpArray[0]-startA;
			startA=tmpArray[0];
			sumStartA=nStartA*startA;
			if ((prefixReplication||singleDeletionReplication) && minPrefixLength>0) minPrefixLength-=delta;
			if (maximalStartA<startA) maximalStartA=-1;
			if (firstMaximalStart<startA) firstMaximalStart=-1;
			isLeftMaximal=false;
		}
		if (tmpArray[1]!=endA) {
			delta=endA-tmpArray[1];
			endA=tmpArray[1];
			sumEndA=nEndA*endA;
			if ((suffixReplication||singleDeletionReplication) && minSuffixLength>0) minSuffixLength-=delta;
			if (maximalEndA>endA) maximalEndA=-1;
			if (lastMaximalEnd>endA) lastMaximalEnd=-1;
			isRightMaximal=false;
		}
	}
	
	
	/**
	 * Tries to transform the replication type of the dense substring (which is assumed 
	 * not to be of substring type) into a more general type, using the B-maximal 
	 * alignments in $alignments$.
	 *
	 * Remark: the procedure does not change the weak type of the dense substring. This
	 * conforms with $DenseSubstrings.getConnectedComponent_merge()$, where merging a non-
	 * weak and a weak substring gives a non-weak substring.
	 *
	 * @param alignmentStats output of $Constants.generalizeReplicationTypes_
	 * getAlignmentStats()$;
	 * @param distanceThreshold threshold on distance identity;
	 * @param alignmentsThreshold minimum number of B-maximal alignments to generalize the 
	 * type;
	 * @param surfaceThreshold the alignments in $alignments$ must cover at least this
	 * fraction of the length of the substring for the type to be generalized;
	 * @param tmpArray temporary space, with at least 5 cells;
	 * @return TRUE iff the type of the substring was generalized.
	 */
	public final boolean generalizeType(int[] alignmentStats, int nAlignments, int distanceThreshold, int alignmentsThreshold, double surfaceThreshold, int[] tmpArray) {
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
		if (surface<(endA-startA+1)*surfaceThreshold) return false;
		
		// Generalizing replication type
		if (prefixReplication && !suffixReplication && !singleDeletionReplication) {
			if (nLeftMaximal<alignmentsThreshold) return false;
			if (leftMaximalAligned) {
				suffixReplication=true;
				minSuffixLength=minSufLength;
				suffixesAreUniformlyDistributed=true;  // We don't try to estimate this.
			}
			else {
				prefixReplication=false;
				substringReplication=true;
				tmpArray[0]=minPrefixLength;
				tmpArray[1]=minPrefLength;
				tmpArray[2]=minSufLength;
				tmpArray[3]=minSubstringLength;
				length=Math.minPositive(tmpArray,3);
				minPrefixLength=length!=-1?length:0;
			}
			firstMaximalStart=Math.min(firstMaximalStart,firstMaxStart);
			nMaximalAlignments+=nAlignments;
			return true;
		}
		else if (suffixReplication && !prefixReplication && !singleDeletionReplication) {
			if (nRightMaximal<alignmentsThreshold) return false;
			if (rightMaximalAligned) {
				prefixReplication=true;
				minPrefixLength=minPrefLength;
				prefixesAreUniformlyDistributed=true;  // We don't try to estimate this.
			}
			else {
				suffixReplication=false;
				substringReplication=true;
				tmpArray[0]=minSuffixLength;
				tmpArray[1]=minPrefLength;
				tmpArray[2]=minSufLength;
				tmpArray[3]=minSubstringLength;
				length=Math.minPositive(tmpArray,3);
				minPrefixLength=length!=-1?length:0;
			}
			lastMaximalEnd=Math.max(lastMaximalEnd,lastMaxEnd);
			nMaximalAlignments+=nAlignments;
			return true;
		}
		else if ((prefixReplication && suffixReplication) || singleDeletionReplication) {
			if (nLeftRightMaximal<alignmentsThreshold) return false;
			prefixReplication=false;
			suffixReplication=false;
			singleDeletionReplication=false;
			substringReplication=true;
			tmpArray[0]=minPrefixLength;
			tmpArray[1]=minSuffixLength;
			tmpArray[2]=minPrefLength;
			tmpArray[3]=minSufLength;
			tmpArray[4]=minSubstringLength;
			length=Math.minPositive(tmpArray,4);
			minPrefixLength=length!=-1?length:0;
			firstMaximalStart=Math.min(firstMaximalStart,firstMaxStart);
			lastMaximalEnd=Math.max(lastMaximalEnd,lastMaxEnd);
			nMaximalAlignments+=nAlignments;
			return true;
		}
		return false;
	}
	
	
	/**
	 * Concatenates $rightSubstring$ (assumed to be of substring type) to the right of the
	 * current substring, which is assumed to be of substring type as well.
	 */
	public final void concatenateSubstringType(DenseSubstring rightSubstring) {
		if (rightSubstring.endA>endA) {
			endA=rightSubstring.endA;
			sumEndA=rightSubstring.sumEndA;
			nEndA=rightSubstring.nEndA;
			isRightMaximal=rightSubstring.isRightMaximal;
		}
		maximalEndA=Math.max(maximalEndA,rightSubstring.maximalEndA);
		lastMaximalEnd=Math.max(lastMaximalEnd,rightSubstring.lastMaximalEnd);
		nMaximalAlignments+=rightSubstring.nMaximalAlignments;
	}
	

}