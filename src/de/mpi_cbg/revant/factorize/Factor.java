package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import java.io.IOException;
import java.io.BufferedWriter;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class Factor implements Comparable {
	/**
	 * Parameters of the pipeline
	 */
	public static int maxBracketingDistance = 200;
	public static double leftRightThreshold = 0.1;
	public static int minMaximalAlignments = 3;  // Minimum number of alignments for declaring maximality
	public static double factorIntersectionFraction = 0.5;  // Used by $compute*Statistics$
	public static final int ANCHORING_DISTANCE = (IO.quantum*3)>>1;

	/**
	 * Sorting status of all factors
	 */
	public static final int UNSORTED = -1;
	public static final int FIRSTPOSITION = 0;
	public static final int LASTPOSITION = 1;
	public static final int SHOULDBEOUTPUT_FIRSTPOSITION = 2;
	public static int order = UNSORTED;

	/**
	 * Properties of a factor
	 */
	public static final int N_NONPERIODIC_COVERAGE_TYPES = 6;  // 0=full; 1=prefix; 2=suffix; 3=one substring; 4=more than one substring (i.e. at least one deletion); 5=partial, but type cannot be decided. -1 for periodic factors.
	public static final int MIN_INTERVALS_PER_FACTOR = 10;
	public int firstPosition, lastPosition, length;
	public int periodicCoverage, periodicCoverageIdentity;
	public int[] nonPeriodicCoverage;  // Periodic factors use just cell 0
	public int[] nonPeriodicCoverageIdentity;  // Periodic factors use just cell 0
	public int period;  // -1: the factor is not periodic; 0: the factor is periodic, but it is not possible to estimate a period; >0: estimate of a period.
	public boolean hasHighQuality;
	public boolean hasHighCoverage;  // Refers to the raw, uncorrected coverage.
	public boolean shouldBeOutput;
	public int nVariants;
	public double maxDiffs;
	public boolean inDenseSubstring;   // The factor is contained in a dense substring
	public boolean inDenseSubstringOfPrefixType;  // The factor is contained in a dense substring that replicates by taking prefixes of itself
	public boolean inDenseSubstringOfSuffixType;  // The factor is contained in a dense substring that replicates by taking suffixes of itself
	public boolean inDenseSubstringOfSubstringType;  // The factor is contained in a dense substring that replicates by taking substrings of itself.
	public boolean inPeriodicSubstringInterval;  // The factor is contained in a short-period interval.
	
	/**
	 * Copied from the splits that generated the factor
	 */
	public int nOpenLeft, nClosedLeft, nUnknownLeft, nOpenRight, nClosedRight, nUnknownRight;
	

	/**
	 * Cells:
	 *
	 * 0: There is an overlap between two occurrences of type 0 (full) in opposite
	 * orientations. 1=both occurrences are exact; 2=one occurrence is exact, the other is
	 * not.
	 *
	 * 1: There is an overlap between an occurrence of type 1 (prefix) and an occurrence
	 * of type 2 (suffix), in opposite orientations. Both occurrences are exact.
	 *
	 * 2: There is an overlap between an occurrence of type 1 (prefix) and an occurrence
	 * of type 2 (suffix), in opposite orientations. One occurrence is exact, the other is
	 * not.
	 *
	 * 3: There is an overlap between two occurrences of type 3 (one substring), in
	 * opposite orientations. Both occurrences are exact.
	 *
	 * 4: There is an overlap between two occurrences of type 3 (one substring), in
	 * opposite orientations. One occurrence is exact, the other is not.
	 *
	 * Periodic factors use just cell 0, with values 0/1.
	 */
	public int[] equalsReverseComplement;
	public boolean hasSelfSimilarity;  // True iff there is an instance of the factor in a readB, for which there is evidence that it can be built with at least two distinct replication modes. Only for nonperiodic factors.

	/**
	 * Number of alignments:
	 * 0: mostly inside the factor;
	 * 1: left alignment;
	 * 2: right alignment;
	 * 3: spanning alignment.
	 */
	public int[] nAlignments;

	/**
	 * Number of alignments that are not implied by a dense or periodic substring, and:
	 *
	 * Whose readA interval has high Jaccard similarity with the factor, and that are:
	 * 0: only left-maximal;
	 * 1: only right-maximal;
	 * 2: both left- and right-maximal.
	 *
	 * 3: left-bracketing;
	 * 4: right-bracketing.
	 */
	public int[] nMaximalAlignments;
	public double avgLeftBracketing, avgRightBracketing;
	public double stdLeftBracketing, stdRightBracketing;
	public boolean isLeftBracketed, isRightBracketed;
	public int leftScore, rightScore;

	/**
	 * Periodic substring statistics
	 */
	public boolean periodicSubstringLeft, periodicSubstringRight, periodicSubstringSpanning;

	/**
	 * Intervals
	 */
	public DenseSubstring[] denseSubstrings;  // All dense substring intervals that have a significant intersection with the current factor and with at least one other factor
	public PeriodicSubstringInterval[] periodicSubstringIntervals;  // All periodic substring intervals that have a significant intersection with the current factor and with at least one other factor
	public AlignmentInterval[] alignmentIntervals;  // All alignment intervals that have a significant intersection with the current factor and with at least one other factor
	public int lastDenseSubstring, lastPeriodicSubstringInterval, lastAlignmentInterval;
	public boolean keepLeftBoundary, keepRightBoundary;
	public static final int BRACKETING_BIN = Alignments.minAlignmentLength>>2;
	public static final int N_BRACKETING_BINS = 40;
	private static final DenseSubstring tmpSubstring = new DenseSubstring();  // Temporary space used for search



	public final void initialize(int firstPosition, int lastPosition) {
		this.firstPosition=firstPosition;
		this.lastPosition=lastPosition;
		length=lastPosition-firstPosition+1;
		hasHighQuality=true;
		hasHighCoverage=true;
		shouldBeOutput=true;

		nonPeriodicCoverage = new int[N_NONPERIODIC_COVERAGE_TYPES];
		nonPeriodicCoverageIdentity = new int[N_NONPERIODIC_COVERAGE_TYPES];
		periodicCoverage=0;
		periodicCoverageIdentity=0;
		period=-1;
		leftScore=0;
		rightScore=0;
		nVariants=0;
		maxDiffs=0.0;
		equalsReverseComplement = new int[5];

		denseSubstrings = new DenseSubstring[MIN_INTERVALS_PER_FACTOR];
		periodicSubstringIntervals = new PeriodicSubstringInterval[MIN_INTERVALS_PER_FACTOR];
		lastDenseSubstring=-1;
		lastPeriodicSubstringInterval=-1;
		inDenseSubstring=false;
		inDenseSubstringOfSubstringType=false;

		nAlignments = new int[4];
		nMaximalAlignments = new int[5];

		alignmentIntervals = new AlignmentInterval[MIN_INTERVALS_PER_FACTOR];
		lastAlignmentInterval=-1;
		
		nOpenLeft=0; nClosedLeft=0; nUnknownLeft=0; 
		nOpenRight=0; nClosedRight=0; nUnknownRight=0;
	}


	public final void clearNAlignments() {
		int i;
		Math.set(nAlignments,3,0);
		Math.set(nMaximalAlignments,4,0);
		leftScore=0; rightScore=0;
		avgLeftBracketing=0.0;
		avgRightBracketing=0.0;
		stdLeftBracketing=0.0;
		stdRightBracketing=0.0;
	}







	// ------------------------------ nMaximalAlignments ---------------------------------


	public final boolean isLeftBracketing(int alignmentStart, int alignmentEnd, boolean isLeftMaximal, Factor previousFactor) {
		if (!isLeftMaximal || alignmentEnd<=lastPosition) return false;
		int distance = alignmentStart-firstPosition;
		if (distance<0) distance=-distance;
		if (distance>maxBracketingDistance) return false;
		return (alignmentStart>=firstPosition && distance<=(length>>1)) ||
		       (alignmentStart<firstPosition && (previousFactor==null||distance<=(previousFactor.length>>1)));
	}


	public final boolean isRightBracketing(int alignmentStart, int alignmentEnd, boolean isRightMaximal, Factor nextFactor) {
		if (!isRightMaximal || alignmentStart>=firstPosition) return false;
		int distance = alignmentEnd-lastPosition;
		if (distance<0) distance=-distance;
		if (distance>maxBracketingDistance) return false;
		return (alignmentEnd<=lastPosition && distance<=(length>>1)) ||
		       (alignmentEnd>lastPosition && (nextFactor==null||distance<=(nextFactor.length>>1)));
	}


	public final int getTotalMaximalAlignments() {
		int out = 0;
		for (int i=0; i<=4; i++) out+=nMaximalAlignments[i];
		return out;
	}


	public final boolean isLeftMaximal() {
		return nMaximalAlignments[0]+nMaximalAlignments[2]+nMaximalAlignments[3]>=minMaximalAlignments;
	}


	public final boolean isRightMaximal() {
		return nMaximalAlignments[1]+nMaximalAlignments[2]+nMaximalAlignments[4]>=minMaximalAlignments;
	}


	public final boolean isMaximalRepeat() {
		return hasHighQuality && isLeftMaximal() && isRightMaximal();
	}




	// -------------------------------- nAlignments --------------------------------------

	public final boolean isRightAlignment(int alignmentStart, int alignmentEnd) {
		return alignmentEnd>lastPosition &&
			   ( alignmentStart>=firstPosition ||
			     ( alignmentStart<firstPosition &&
				   firstPosition-alignmentStart<=(alignmentEnd-alignmentStart+1)*leftRightThreshold &&
			       alignmentStart>IO.quantum &&
			       alignmentEnd-lastPosition>firstPosition-alignmentStart
			     )
			   );
	}


	public final boolean isLeftAlignment(int alignmentStart, int alignmentEnd, int readLength) {
		return alignmentStart<firstPosition &&
			   ( alignmentEnd<=lastPosition ||
				 ( alignmentEnd>lastPosition &&
				   alignmentEnd-lastPosition<=(alignmentEnd-alignmentStart+1)*leftRightThreshold &&
				   alignmentEnd<readLength-IO.quantum &&
				   firstPosition-alignmentStart>alignmentEnd-lastPosition
				 )
			   );
	}


	public final boolean isSpanningAlignment(int alignmentStart, int alignmentEnd, int readLength) {
		double small = (alignmentEnd-alignmentStart+1)*leftRightThreshold;
		return alignmentStart<firstPosition &&
			   alignmentEnd>lastPosition &&
			   alignmentEnd-lastPosition>=small &&
			   firstPosition-alignmentStart>=small;
	}




	/**
	 * @param npc deltas (possibly negative) to $nonPeriodicCoverage$;
	 * @param npci deltas (possibly negative) to $nonPeriodicCoverageIdentity$;
	 * @param pc delta (possibly negative) to $periodicCoverage$;
	 * @param pci delta (possibly negative) to $periodicCoverageIdentity$.
	 */
	public final void incrementCorrectedCoverage(int[] npc, int[] npci, int pc, int pci) {
		int i;

		for (i=0; i<npc.length; i++) nonPeriodicCoverage[i]+=npc[i];
		for (i=0; i<npci.length; i++) nonPeriodicCoverageIdentity[i]+=npci[i];
		periodicCoverage+=pc;
		periodicCoverageIdentity+=pci;
	}


	/**
	 * We count all occurrences, of any replication type, of any length.
	 */
	public final int getNonPeriodicCoverage() {
		int out = 0;
		for (int i=0; i<N_NONPERIODIC_COVERAGE_TYPES; i++) out+=nonPeriodicCoverage[i];
		return out;
	}


	/**
	 * We count all occurrences, of any replication type, of any length.
	 */
	public final int getNonPeriodicCoverageIdentity() {
		int out = 0;
		for (int i=0; i<N_NONPERIODIC_COVERAGE_TYPES; i++) out+=nonPeriodicCoverageIdentity[i];
		return out;
	}


	public final int getCorrectedCoverage() {
		return periodicCoverage+getNonPeriodicCoverage();
	}


	/**
	 * Remark: the procedure assumes that all dense substrings have been assigned an ID, 
	 * and that $denseSubstrings$ has been sorted by ID.
	 *
	 * @return the position in array $denseSubstrings$ at which the dense substring with
	 * ID $id$ is located, or -1 if no such position can be found.
	 */
	public int inDenseSubstring(int id) {
		int out, tmpCriterion;

		tmpCriterion=DenseSubstring.order;  // $DenseSubstring.order=DenseSubstring.ID$ is not necessarily true.
		DenseSubstring.order=DenseSubstring.ID;
		tmpSubstring.id=id;
		out=Arrays.binarySearch(denseSubstrings,0,lastDenseSubstring+1,tmpSubstring);
		if (out<0) out=-1;
		DenseSubstring.order=tmpCriterion;
		return out;
	}


	public final int compareTo(Object other) {
		Factor otherFactor = (Factor)other;
		if (order==FIRSTPOSITION) {
			if (firstPosition<otherFactor.firstPosition) return -1;
			else if (firstPosition>otherFactor.firstPosition) return 1;
		}
		else if (order==LASTPOSITION) {
			if (lastPosition<otherFactor.lastPosition) return -1;
			else if (lastPosition>otherFactor.lastPosition) return 1;
		}
		else if (order==SHOULDBEOUTPUT_FIRSTPOSITION) {
			if (shouldBeOutput && !otherFactor.shouldBeOutput) return -1;
			else if (!shouldBeOutput && otherFactor.shouldBeOutput) return 1;
			if (firstPosition<otherFactor.firstPosition) return -1;
			else if (firstPosition>otherFactor.firstPosition) return 1;
		}
		return 0;
	}


	public String toString() {
		return "["+firstPosition+".."+lastPosition+"] nOL="+nOpenLeft+" nCL="+nClosedLeft+" nOR="+nOpenRight+" nCR="+nClosedRight;
	}


	public String toStringStatistics() {
		int i;
		String out = "";
		for (i=0; i<3; i++) out+=nAlignments[i]+",";
		out+=nAlignments[3]+"|";
		for (i=0; i<4; i++) out+=nMaximalAlignments[i]+",";
		out+=nMaximalAlignments[4];
		return out+", leftScore="+leftScore+", rightScore="+rightScore+", periodicSubstringLeft="+periodicSubstringLeft+", periodicSubstringRight="+periodicSubstringRight+", periodicSubstringSpanning="+periodicSubstringSpanning+", length="+length+", correctedCoverage="+getCorrectedCoverage()+" lastDenseSubstring="+lastDenseSubstring+" inDenseSubstringOfSubstringType="+inDenseSubstringOfSubstringType+" isLeftMaximal="+isLeftMaximal()+" isRightMaximal="+isRightMaximal();
	}


	/**
	 * Remark: $period$ has value -1 if the factor is not periodic; it has value 0 if
	 * the factor is periodic and it has not been possible to estimate its period.
	 * $periodicCoverage$ has value 0 if the factor is not periodic.
	 * The maximum value of $nVariants$ is decided by constant $MAX_NPEAKS$ inside procedure
	 * $ReadA.estimateMaxDiffs()$.
	 */
	public final void writeFactorsFile(BufferedWriter file) throws IOException {
		int i;

		file.write(firstPosition+",");
		file.write(lastPosition+"|");
		file.write(isMaximalRepeat()?"1|":"0|");
		file.write(period+"|");
		for (i=0; i<N_NONPERIODIC_COVERAGE_TYPES-1; i++) {
			file.write(nonPeriodicCoverage[i]+",");
			file.write(nonPeriodicCoverageIdentity[i]+",");
		}
		file.write(nonPeriodicCoverage[N_NONPERIODIC_COVERAGE_TYPES-1]+",");
		file.write(nonPeriodicCoverageIdentity[N_NONPERIODIC_COVERAGE_TYPES-1]+"|");
		file.write(periodicCoverage+",");
		file.write(periodicCoverageIdentity+"|");
		for (i=0; i<equalsReverseComplement.length-1; i++) file.write(equalsReverseComplement[i]+",");
		file.write(equalsReverseComplement[equalsReverseComplement.length-1]+"|");
		file.write(hasSelfSimilarity?"1|":"0|");
		file.write(nVariants+"\n");
	}


	public final void writeIntervalsDenseFile(BufferedWriter file) throws IOException {
		for (int i=0; i<=lastDenseSubstring; i++) file.write(denseSubstrings[i].id+",");
		file.write("\n");
	}


	public final void writeIntervalsPeriodicFile(BufferedWriter file) throws IOException {
		for (int i=0; i<=lastPeriodicSubstringInterval; i++) file.write(periodicSubstringIntervals[i].id+",");
		file.write("\n");
	}


	public final void writeIntervalsAlignmentFile(BufferedWriter file) throws IOException {
		for (int i=0; i<=lastAlignmentInterval; i++) file.write(alignmentIntervals[i].id+",");
		file.write("\n");
	}


	/**
	 * Updates $nAlignments$, $nMaximalAlignments$, $leftScore$, $rightScore$, using an
	 * alignment.
	 *
	 * @param implied TRUE iff the alignment is implied by a dense or periodic substring.
	 */
	public final void updateAlignmentStatistics(int alignmentStart, int alignmentEnd, boolean isLeftMaximal, boolean isRightMaximal, boolean implied, Factor previousFactor, Factor nextFactor) {
		final boolean contained = Intervals.isApproximatelyContained(alignmentStart,alignmentEnd,firstPosition,lastPosition);
		
		// Updating $nAlignments$, $leftScore$, $rightScore$.
		if (contained) nAlignments[0]++;
		else if (isLeftAlignment(alignmentStart,alignmentEnd,ReadA.readLength)) {
			nAlignments[1]++;
			leftScore+=Intervals.intersectionLength(alignmentStart,alignmentEnd,firstPosition,lastPosition);
		}
		else if (isRightAlignment(alignmentStart,alignmentEnd)) {
			nAlignments[2]++;
			rightScore+=Intervals.intersectionLength(alignmentStart,alignmentEnd,firstPosition,lastPosition);
		}
		else if (isSpanningAlignment(alignmentStart,alignmentEnd,ReadA.readLength)) nAlignments[3]++;

		// Updating $nMaximalAlignments$
		if (implied || (!isLeftMaximal && !isRightMaximal)) return;	
		if (contained || Intervals.jaccardSimilarity(alignmentStart,alignmentEnd,firstPosition,lastPosition)>=Intervals.jaccardThreshold) {
			if (isLeftMaximal && !isRightMaximal) nMaximalAlignments[0]++;
			else if (!isLeftMaximal && isRightMaximal) nMaximalAlignments[1]++;
			else if (isLeftMaximal && isRightMaximal) nMaximalAlignments[2]++;
		}
		else {
			if (isLeftBracketing(alignmentStart,alignmentEnd,isLeftMaximal,previousFactor)) {
				nMaximalAlignments[3]++;
				avgLeftBracketing+=alignmentStart;
			}
			else if (isRightBracketing(alignmentStart,alignmentEnd,isRightMaximal,nextFactor)) {
				nMaximalAlignments[4]++;
				avgRightBracketing+=alignmentEnd;
			}
		}
	}


	/**
	 * Updates $periodicSubstringLeft$, $periodicSubstringRight$,
	 * $periodicSubstringSpanning$, using a given periodic substring.
	 */
	public final void updatePeriodicSubstringStatistics(int substringStart, int substringEnd) {
		int substringLength;

		substringLength=substringEnd-substringStart+1;
		if ( Intervals.isApproximatelyContained(substringStart,substringEnd,firstPosition,lastPosition) ||
		     Intervals.intersectionLength(substringStart,substringEnd,firstPosition,lastPosition)<factorIntersectionFraction*length) return;
		if ( (substringStart>=firstPosition && substringEnd>=lastPosition) ||
			 ( substringStart<firstPosition && substringEnd>lastPosition &&
			   firstPosition-substringStart<=substringLength*leftRightThreshold &&
			   substringStart>IO.quantum &&
			   substringEnd-lastPosition>firstPosition-substringStart ) ) periodicSubstringRight=true;
		else if ( (substringStart<=firstPosition && substringEnd<=lastPosition) ||
				  ( substringStart<firstPosition && substringEnd>lastPosition &&
					substringEnd-lastPosition<=substringLength*leftRightThreshold &&
					substringEnd<ReadA.readLength-IO.quantum &&
					firstPosition-substringStart>substringEnd-lastPosition ) ) periodicSubstringLeft=true;
		else if (substringStart<firstPosition && substringEnd>lastPosition) periodicSubstringSpanning=true;
	}


	/**
	 * Returns TRUE iff $denseSubstrings$ and $otherFactor.denseSubstrings$ are identical.
	 *
	 * Remark: the procedure assumes $denseSubstrings$ and $otherFactor.denseSubstrings$
	 * to have been sorted by the same criterion, and it assumes dense substrings to have
	 * been assigned an ID.
	 */
	public final boolean sameDenseSubstrings(Factor otherFactor) {
		if (lastDenseSubstring!=otherFactor.lastDenseSubstring) return false;
		for (int i=0; i<=lastDenseSubstring; i++) {
			if (denseSubstrings[i].id!=otherFactor.denseSubstrings[i].id) return false;
		}
		return true;
	}


	/**
	 * Returns TRUE iff $periodicSubstringIntervals$ and
	 * $otherFactor.periodicSubstringIntervals$ are identical.
	 *
	 * Remark: the procedure assumes $periodicSubstringIntervals$ and
	 * $otherFactor.periodicSubstringIntervals$ to have been sorted by the same criterion,
	 * and it assumes periodic substring intervals to have been assigned an ID.
	 */
	public final boolean samePeriodicSubstringIntervals(Factor otherFactor) {
		if (lastPeriodicSubstringInterval!=otherFactor.lastPeriodicSubstringInterval) return false;
		for (int i=0; i<=lastPeriodicSubstringInterval; i++) {
			if (periodicSubstringIntervals[i].id!=otherFactor.periodicSubstringIntervals[i].id) return false;
		}
		return true;
	}


	/**
	 * Returns TRUE iff $alignmentIntervals$ and $otherFactor.alignmentIntervals$ are
	 * identical.
	 *
	 * Remark: the procedure assumes $alignmentIntervals$ and
	 * $otherFactor.alignmentIntervals$ to have been sorted by the same criterion, and it 
	 * assumes alignment intervals to have been assigned an ID.
	 */
	public final boolean sameAlignmentIntervals(Factor otherFactor) {
		if (lastAlignmentInterval!=otherFactor.lastAlignmentInterval) return false;
		for (int i=0; i<=lastAlignmentInterval; i++) {
			if (alignmentIntervals[i].id!=otherFactor.alignmentIntervals[i].id) return false;
		}
		return true;
	}
	
	
	/**
	 * @return TRUE iff there is at least one dense substring of substring type in
	 * $denseSubstrings$ with $observed>=minNAlignments{Left,Right}$.
	 */
	public final boolean hasLargeNAlignments(boolean left, int observed) {
		int i;
		
		if (left) {
			for (i=0; i<=lastDenseSubstring; i++) {
				if (denseSubstrings[i].substringReplication && denseSubstrings[i].minNAlignmentsLeft!=-1 && observed>=denseSubstrings[i].minNAlignmentsLeft) return true;
			}
		}
		else {
			for (i=0; i<=lastDenseSubstring; i++) {
				if (denseSubstrings[i].substringReplication && denseSubstrings[i].minNAlignmentsRight!=-1 && observed>=denseSubstrings[i].minNAlignmentsRight) return true;
			}
		}
		return false;
	}
		
		
	/**
	 * @return TRUE iff there is at least one dense substring of substring type in
	 * $denseSubstrings$ with the current factor as the first factor from the left/right.
	 */
	public final boolean isFirstFactor(boolean left) {
		int i;
	
		if (left) {
			for (i=0; i<=lastDenseSubstring; i++) {
				if (denseSubstrings[i].substringReplication && denseSubstrings[i].firstFactor!=-1 && Factors.factors[denseSubstrings[i].firstFactor]==this) return true;
			}
		}
		else {
			for (i=0; i<=lastDenseSubstring; i++) {
				if (denseSubstrings[i].substringReplication && denseSubstrings[i].lastFactor!=-1 && Factors.factors[denseSubstrings[i].lastFactor]==this) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * @param leftOrRight TRUE=beginning of the factor; FALSE=end of the factor;
	 * @return TRUE iff the side of the factor conforms to the $n{Open,Closed}$ 
	 * information from the splits that generated the factor.
	 */
	public final boolean checkSplits(boolean leftOrRight) {
		return leftOrRight?(nOpenLeft>0||nClosedLeft==0):(nClosedRight>0||nOpenRight==0);
	}


	public final void ensureSpace_denseSubstrings(int n) {
		final double GROWTH_RATE = 1.5;  // Arbitrary
		if (denseSubstrings.length>=n) return;
		
		DenseSubstring[] newDenseSubstrings = new DenseSubstring[(int)(n*GROWTH_RATE)];
		System.arraycopy(denseSubstrings,0,newDenseSubstrings,0,denseSubstrings.length);
		for (int i=denseSubstrings.length; i<newDenseSubstrings.length; i++) newDenseSubstrings[i] = new DenseSubstring();
		denseSubstrings=newDenseSubstrings;
	}
	
	
	public final void ensureSpace_periodicSubstringIntervals(int n) {
		final double GROWTH_RATE = 1.5;  // Arbitrary
		if (periodicSubstringIntervals.length>=n) return;
		
		PeriodicSubstringInterval[] newPeriodicSubstringIntervals = new PeriodicSubstringInterval[(int)(n*GROWTH_RATE)];
		System.arraycopy(periodicSubstringIntervals,0,newPeriodicSubstringIntervals,0,periodicSubstringIntervals.length);
		for (int i=periodicSubstringIntervals.length; i<newPeriodicSubstringIntervals.length; i++) newPeriodicSubstringIntervals[i] = new PeriodicSubstringInterval();
		periodicSubstringIntervals=newPeriodicSubstringIntervals;
	}
	
	
	public final void ensureSpace_alignmentIntervals(int n) {
		final double GROWTH_RATE = 1.5;  // Arbitrary
		if (alignmentIntervals.length>=n) return;
		
		AlignmentInterval[] newAlignmentIntervals = new AlignmentInterval[(int)(n*GROWTH_RATE)];
		System.arraycopy(alignmentIntervals,0,newAlignmentIntervals,0,alignmentIntervals.length);
		for (int i=alignmentIntervals.length; i<newAlignmentIntervals.length; i++) newAlignmentIntervals[i] = new AlignmentInterval();
		alignmentIntervals=newAlignmentIntervals;
	}


}