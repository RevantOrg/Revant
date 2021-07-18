package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class Occurrences {
	/**
	 * Parameters of the pipeline
	 */
	public static double overlapThreshold = 0.5;  // Used by $getCost$ to detect highly-overlapping alignments in $readB$. Arbitrary.
	public static double surfaceThreshold = Intervals.jaccardThreshold;  // We say that a factor of $readA$ occurs in a given $readB$ if at least $surfaceThreshold$ fraction of the factor is covered by alignments to a substring of $readB$, possibly broken down by long random insertions.
	public static double exactSurfaceThreshold = 0.5;
	public static double factorUnitFraction = 0.2;  // Unit of length, as a fraction of the length of a factor. Arbitrary.
	
	/**
	 * Temporary space
	 */
	private static int[] deltaNonperiodicCoverage, deltaNonperiodicCoverageIdentity;
	private static int deltaPeriodicCoverage, deltaPeriodicCoverageIdentity;
	private static Occurrence[] occurrences;
	private static int lastOccurrence;
	private static boolean[] hasZeroCost;
	private static int[] tmp = new int[10];

	
	public static final void allocateMemory(int maxAlignments, int maxOccurrences) {
		int i;
		
		occurrences = new Occurrence[maxOccurrences];
		lastOccurrence=-1;
		for (i=0; i<maxOccurrences; i++) occurrences[i] = new Occurrence();
		deltaNonperiodicCoverage = new int[Factor.N_NONPERIODIC_COVERAGE_TYPES];
		deltaNonperiodicCoverageIdentity = new int[Factor.N_NONPERIODIC_COVERAGE_TYPES];
		hasZeroCost = new boolean[maxAlignments];
	}


	/**
	 * Computes a corrected version of the coverage of a factor. This integer number
	 * supersedes the (rational) average coverage of the factor computed from the
	 * $ReadA.coverage$ histogram.
	 *
	 * 1. If a factor $W$ is not periodic, its corrected coverage is the number of
	 * distinct subsequences $V$, in the read set, such that: (1) the set of all
	 * alignments between $V$ and $W$ covers at least a fraction $surfaceThreshold$ of
	 * $W$; (2) either $V$ has no gaps (i.e. it is a substring of a read), or all gaps
	 * between the solid characters of $V$ are long low-quality insertions; (3) the
	 * intervals in $readB$ of the alignments between $V$ and $W$ have approximately the
	 * same order as the corresponding intervals in $readA$ (see procedure
	 * $estimateOccurrences$).
	 *
	 * 2. If a factor is short-period periodic with period $p$, its corrected coverage is 
	 * the number of occurrences, in the read set, of maximal substrings with period $p$,
	 * possibly fragmented by long, low-quality insertions, and possibly longer or shorter
	 * than the factor.
	 *
	 * 3. A factor can be part of a periodic substring of $readA$, and at the same time
	 * occur outside periodic substrings in other reads (e.g. the factor could be a
	 * specific instance of the period of a periodic substring of $readA$, which occurs by
	 * itself in other reads). In this case, the corrected coverage of the factor is the
	 * sum of (1) and (2).
	 *
	 * Replacing the average coverage in $ReadA.coverage$ with the corrected coverage is
	 * necessary for the following reasons:
	 *
 	 * 1. The coverage of a position of a factor that comes from a dense substring is a
	 * function of the number of other reads in which the corresponding repeat occurs, and
	 * of the probability that an alignment covers that position (which is itself a
	 * function of the length of the factor and of the average length of an alignment).
	 * This value could be smaller than the number of occurrences of the repeat in the
	 * read set.
	 *
	 * 2. The coverage of a position of a factor that comes from a periodic substring
	 * equals the number of alignments that cover the position. This number is
	 * proportional to the length of the substring in a $readB$, divided by the period,
	 * which might be larger than the number of occurrences of the repeated periodic
	 * substring in the read set.
	 *
	 * A $readB$ contributes $k_f+k_r$ to the corrected coverage of a factor, where $k_f$
	 * (respectively, $k_r$) is an estimate of the number of occurrences of the entire
	 * factor, possibly broken down in fragments (some of which possibly shorter than
	 * $IO.minAlignmentLength$), in the forward (respectively, reverse complement)
	 * direction of $readB$. Such estimates are obtained by computing the length of a
	 * shortest path in a DAG whose vertices are alignments between $readA$ and $readB$,
	 * and whose arcs have cost one iff they connect pairs of alignments that induce
	 * specific shifts in $readA$ and $readB$: see procedure $estimateOccurrences$.
	 *
	 * Corrected coverage is used by this procedure to decide whether to mark a factor as
	 * low-coverage or not. Field $hasHighCoverage$ of a factor is set to false iff its
	 * corrected coverage is smaller than $Factors.minRepeatCoverage$ and if the factor is
	 * not periodic.
	 *
	 * The procedure tries to correct for factors that are the reverse-complement of
	 * themselves, by making two occurrences that overlap in $readB$ and are in opposite
	 * orientations count as one: see procedure $estimateNOccurrences$.
	 *
	 * Remark: the procedure assumes that $markLowQualityFactors$ has already been
	 * executed.
	 *
	 * Remark: the procedure tries also to estimate the number of exact occurrences of a
	 * factor (as opposed to mutated variants). A periodic substring that intersects a
	 * factor is assumed to contain an exact occurrence of the factor if there is an
	 * alignment implied by the periodic substring, mostly inside the factor, and with
	 * relative $diffs$ at most the $maxDiffs$ of the factor. This is because we assume
	 * that the periods of a periodic factor are approximately uniform.
	 * A nonperiodic factor is assumed to have an exact occurrence, iff the fraction of
	 * the surface of the factor that is covered by alignments with the occurrence, with
	 * relative $diffs$ at most the $maxDiffs$ of the factor, is high enough.
	 *
	 * Remark: the procedure tries also to break down the corrected coverage of a
	 * nonperiodic factor by \emph{replication type}, where a replication type is a way
	 * of copying the factor by applying a set of deletions. See procedures
	 * $estimateOccurrences$ and $estimateNOccurrences$.
	 */
	public static final void estimateCorrectedCoverage() {
		final int MIN_COVERAGE = IO.minOccurrencesInGenome*IO.coverage-1;
		int i, j, factor;
		int firstAlignment, firstAlignmentPrime, firstAlignmentForNextFactor, lastAlignment;
		int firstPeriodic, firstPeriodicPrime, firstPeriodicForNextFactor, lastPeriodic;
		int readB, previousReadB, previousOrientation;
		int blockStart, nOccurrences, nUsedAlignments;
		int startA, endA;
		double diffs;
		PeriodicSubstring ps;

		// Temporarily cloning long-period substrings onto $substrings$.
		PeriodicSubstrings.mergePeriodicSubstrings();

		// Ensuring the necessary order in $PeriodicSubstrings.substrings$ and
		// $ReadA.sortedAlignments$.
		PeriodicSubstring.order=PeriodicSubstring.MINSTARTA_MAXENDA;
		PeriodicSubstring.order_longPeriod=PeriodicSubstring.MINSTARTA_MAXENDA;
		if (PeriodicSubstrings.lastSubstring>0) Arrays.sort(PeriodicSubstrings.substrings,0,PeriodicSubstrings.lastSubstring+1);
		if (Alignment.order!=Alignment.STARTA_ENDA) {
			Alignment.order=Alignment.STARTA_ENDA;
			if (ReadA.lastSortedAlignment>0) Arrays.sort(ReadA.sortedAlignments,0,ReadA.lastSortedAlignment+1);
		}

		// Computing corrected coverages
		firstPeriodic=0; firstAlignment=0;
		for (factor=0; factor<=Factors.lastFactor; factor++) {
			Factors.factors[factor].periodicCoverage=0;
			Factors.factors[factor].periodicCoverageIdentity=0;
			Math.set(Factors.factors[factor].nonPeriodicCoverage,Factor.N_NONPERIODIC_COVERAGE_TYPES-1,0);
			Math.set(Factors.factors[factor].nonPeriodicCoverageIdentity,Factor.N_NONPERIODIC_COVERAGE_TYPES-1,0);
			Factors.factors[factor].hasHighCoverage=false;
			Math.set(Factors.factors[factor].equalsReverseComplement,Factors.factors[factor].equalsReverseComplement.length-1,0);
			Factors.factors[factor].hasSelfSimilarity=false;
			if (!Factors.factors[factor].hasHighQuality && ReadA.getAverageCoverage(Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition)<MIN_COVERAGE) continue;

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("ECC> checking whether factor ["+Factors.factors[factor].firstPosition+".."+Factors.factors[factor].lastPosition+"] intersects periodic substrings...");

			// Marking all representative periodic substrings $W$ such that there is at
			// least one alignment implied by $W$, between the factor and $W$, mostly
			// intersecting the factor, and with $diffs$ at most the $maxDiffs$ of the
			// factor.
			estimateCorrectedCoverage_markPeriodic(factor,firstAlignment,PeriodicSubstrings.substrings,PeriodicSubstrings.lastSubstring);

			// Collecting intersections with periodic substrings
			firstPeriodicPrime=firstPeriodic;
			estimateCorrectedCoverage_periodicCoverage(factor,PeriodicSubstrings.substrings,PeriodicSubstrings.lastSubstring,firstPeriodic,tmp);
			firstPeriodic=tmp[0]; lastPeriodic=tmp[1];

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("ECC> considering factor ["+Factors.factors[factor].firstPosition+".."+Factors.factors[factor].lastPosition+"]...");

			// Marking suitable alignments in $ReadA.sortedAlignments$. An alignment is
			// marked iff it is not implied by a periodic substring, and its
			// intersection with $factor$ is large compared to the length of the alignment
			// (i.e. if the alignment lies inside $factor$, or it approximately coincides
			// with $factor$), or to the length of $factor$ (i.e. if $factor$ lies inside
			// the alignment).
			//
			// Remark: an alignment that is implied by a dense substring is used.
			// Remark: an alignment that is contained in another alignment is used, since
			// it can be useful for estimating the exact surface of an occurrence (see
			// $ReadA.removeContainedAlignments$).
			nUsedAlignments=0;
			firstAlignmentPrime=firstAlignment; firstAlignmentForNextFactor=-1;
			for (i=firstAlignment; i<=ReadA.lastSortedAlignment; i++) {
				startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
				if (startA>Factors.factors[factor].lastPosition) {
					if (firstAlignmentForNextFactor!=-1) firstAlignment=firstAlignmentForNextFactor;
					else firstAlignment=i;
					break;
				}
				endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
				if (endA<Factors.factors[factor].firstPosition) {
					ReadA.sortedAlignments[i].shouldBeUsed=false;
					continue;
				}
				if (firstAlignmentForNextFactor==-1 && factor<Factors.lastFactor && endA>=Factors.factors[factor+1].firstPosition) firstAlignmentForNextFactor=i;
				if (ReadA.sortedAlignments[i].impliedByPeriodicSubstring!=null) {
					ReadA.sortedAlignments[i].shouldBeUsed=false;
					continue;
				}
				if ( !Intervals.isApproximatelyContained(startA,endA,Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition) &&
					 !Intervals.areApproximatelyIdentical(startA,endA,Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition) &&
					 !Intervals.isApproximatelyContained(Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition,startA,endA) 
				   ) {
					ReadA.sortedAlignments[i].shouldBeUsed=false;
					continue;
				}
				ReadA.sortedAlignments[i].shouldBeUsed=true;
				nUsedAlignments++;
			}
			if (i>ReadA.lastSortedAlignment) {
				if (firstAlignmentForNextFactor!=-1) firstAlignment=firstAlignmentForNextFactor;
				else firstAlignment=i;
			}

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("ECC> nUsedAlignments="+nUsedAlignments+" firstAlignment="+firstAlignmentPrime+" lastAlignment="+(i-1));

			if (nUsedAlignments==0) {
				Factors.factors[factor].hasHighCoverage=Factors.factors[factor].periodicCoverage>Factors.minPeriodicCoverage;
				continue;
			}
			lastAlignment=i-1;

			// Sorting the alignments $X$ of $factor$ by $shouldBeUsed$, $readB$,
			// $orientation$, $f_X$ (see procedure $estimateOccurrences$).
			Alignment.order=Alignment.SHOULDBEUSED_READB_ORIENTATION_F;
			Alignment.factorStart=Factors.factors[factor].firstPosition;
			if (lastAlignment+1-firstAlignmentPrime>1) Arrays.sort(ReadA.sortedAlignments,firstAlignmentPrime,lastAlignment+1);

			// Sorting $PeriodicSubstrings.substrings[firstPeriodicPrime..lastPeriodic]$
			// by $readB$ and $minStartBForward$ (see procedure $estimateNOccurrences$).
			PeriodicSubstring.order=PeriodicSubstring.READB_MINSTARTB_FORWARD;
			PeriodicSubstring.order_longPeriod=PeriodicSubstring.READB_MINSTARTB_FORWARD;
			if (lastPeriodic+1-firstPeriodicPrime>1) Arrays.sort(PeriodicSubstrings.substrings,firstPeriodicPrime,lastPeriodic+1);

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("ECC> up to now, factor.periodicCoverage="+Factors.factors[factor].periodicCoverage+" periodicCoverageIdentity="+Factors.factors[factor].periodicCoverageIdentity+" nonPeriodicCoverage="+Factors.factors[factor].nonPeriodicCoverage+" nonPeriodicCoverageIdentity="+Factors.factors[factor].nonPeriodicCoverageIdentity);

			// Estimating the number of occurrences of the factor in each $readB$
			blockStart=firstAlignmentPrime;
			readB=Alignments.alignments[ReadA.sortedAlignments[firstAlignmentPrime].id][1]-1;
			previousReadB=readB;
			previousOrientation=Alignments.alignments[ReadA.sortedAlignments[firstAlignmentPrime].id][2];
			lastOccurrence=-1;
			j=firstPeriodicPrime;  // Pointer to $PeriodicSubstrings.substrings$
			for (i=firstAlignmentPrime+1; i<=lastAlignment; i++) {
				if (!ReadA.sortedAlignments[i].shouldBeUsed) break;
				readB=Alignments.alignments[ReadA.sortedAlignments[i].id][1]-1;
				if (readB!=previousReadB) {
					estimateOccurrences(blockStart,i-1,Factors.factors[factor]);
					j=readB2firstPeriodic(previousReadB,j,lastPeriodic);
					estimateNOccurrences(Factors.factors[factor],previousReadB,j);
					Factors.factors[factor].incrementCorrectedCoverage(deltaNonperiodicCoverage,deltaNonperiodicCoverageIdentity,deltaPeriodicCoverage,deltaPeriodicCoverageIdentity);
					if (IO.SHOW_STD_ERR) IO.printErr("ECC> read="+previousReadB+", in both orientations, contributes the following deltaCoverage: "+deltaNonperiodicCoverage+","+deltaNonperiodicCoverageIdentity+"|"+deltaPeriodicCoverage+","+deltaPeriodicCoverageIdentity);
					lastOccurrence=-1;
					blockStart=i;
					previousReadB=Alignments.alignments[ReadA.sortedAlignments[i].id][1]-1;
					previousOrientation=Alignments.alignments[ReadA.sortedAlignments[i].id][2];
				}
				else if (Alignments.alignments[ReadA.sortedAlignments[i].id][2]!=previousOrientation) {
					estimateOccurrences(blockStart,i-1,Factors.factors[factor]);
					blockStart=i;
					previousOrientation=Alignments.alignments[ReadA.sortedAlignments[i].id][2];
				}
			}
			estimateOccurrences(blockStart,lastAlignment,Factors.factors[factor]);
			j=readB2firstPeriodic(previousReadB,j,lastPeriodic);
			estimateNOccurrences(Factors.factors[factor],previousReadB,j);
			Factors.factors[factor].incrementCorrectedCoverage(deltaNonperiodicCoverage,deltaNonperiodicCoverageIdentity,deltaPeriodicCoverage,deltaPeriodicCoverageIdentity);
			if (IO.SHOW_STD_ERR) IO.printErr("ECC> read="+previousReadB+", in both orientations, contributes the following deltaCoverage: "+deltaNonperiodicCoverage+","+deltaNonperiodicCoverageIdentity+"|"+deltaPeriodicCoverage+","+deltaPeriodicCoverageIdentity);
			Factors.factors[factor].hasHighCoverage=(Factors.factors[factor].getCorrectedCoverage()>=IO.minRepeatCoverage)||(Factors.factors[factor].periodicCoverage>Factors.minPeriodicCoverage);

			// Sorting back (not necessarily in a stable way) interval
			// $PeriodicSubstrings.substrings[firstPeriodicPrime..lastPeriodic]$ by
			// $minStartA$ and $maxEndA$, for the next iteration.
			PeriodicSubstring.order=PeriodicSubstring.MINSTARTA_MAXENDA;
			PeriodicSubstring.order_longPeriod=PeriodicSubstring.MINSTARTA_MAXENDA;
			if (lastPeriodic+1-firstPeriodicPrime>1) Arrays.sort(PeriodicSubstrings.substrings,firstPeriodicPrime,lastPeriodic+1);

			// Sorting back (not necessarily in a stable way) all the alignments of
			// $factor$ by $startA$, then by $endA$, for the next iteration.
			Alignment.order=Alignment.STARTA_ENDA;
			if (lastAlignment+1-firstAlignmentPrime>1) Arrays.sort(ReadA.sortedAlignments,firstAlignmentPrime,lastAlignment+1);
		}
		
		// Removing long-period substrings from $substrings$.
		PeriodicSubstrings.undo_mergePeriodicSubstrings();
	}
	
	
	/**
	 * Marks all representative periodic substrings $W$ such that there is at least one 
	 * alignment implied by $W$, between $factor$ and $W$, mostly intersecting $factor$, 
	 * and with $diffs$ at most the $maxDiffs$ of the factor.
	 */
	private static final void estimateCorrectedCoverage_markPeriodic(int factor, int firstAlignment, PeriodicSubstring[] periodicSubstrings, int lastSubstring) {
		int i;
		int startA, endA;
		double diffs;
		PeriodicSubstring ps;
		
		for (i=0; i<=lastSubstring; i++) periodicSubstrings[i].containsExactOccurrence=false;
		for (i=firstAlignment; i<=ReadA.lastSortedAlignment; i++) {
			startA=Alignments.alignments[ReadA.sortedAlignments[i].id][3];
			if (startA>Factors.factors[factor].lastPosition) break;
			endA=Alignments.alignments[ReadA.sortedAlignments[i].id][4];
			if (endA<Factors.factors[factor].firstPosition) continue;
			ps=ReadA.sortedAlignments[i].impliedByPeriodicSubstring;
			if (ps!=null && ps.representative!=null && ps.representative.threadStart!=null) {
				if ( !Intervals.isApproximatelyContained(startA,endA,Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition) &&
					 !Intervals.areApproximatelyIdentical(startA,endA,Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition) &&
					 !Intervals.isApproximatelyContained(Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition,startA,endA)
				   ) continue;
				diffs=Alignments.getAvgDiffs(ReadA.sortedAlignments[i].id);
				if (diffs>Factors.factors[factor].maxDiffs) continue;
				ps.representative.threadStart.containsExactOccurrence=true;  // $ps$ might not be representative
			}
		}
	}
	
	
	/**
	 * Collects intersections with periodic substrings.
	 *
	 * @param firstPeriodic the first element to consider in $periodicSubstrings$;
	 * @return out 0: the value of $firstPeriodic$ to be used with the next factor; 
	 * 1: the last element of $periodicSubstrings$ considered by this procedure.
	 */
	private static final void estimateCorrectedCoverage_periodicCoverage(int factor, PeriodicSubstring[] periodicSubstrings, int lastSubstring, int firstPeriodic, int[] out) {
		int i;
		int firstPeriodicForNextFactor, firstPeriodicPrime;
		
		firstPeriodicForNextFactor=-1;
		for (i=firstPeriodic; i<=lastSubstring; i++) {
			if (periodicSubstrings[i].minStartA>=Factors.factors[factor].lastPosition) {
				if (firstPeriodicForNextFactor!=-1) firstPeriodic=firstPeriodicForNextFactor;
				else firstPeriodic=i;
				break;
			}
			if (periodicSubstrings[i].maxEndA<=Factors.factors[factor].firstPosition) continue;
			if (factor<Factors.lastFactor && firstPeriodicForNextFactor==-1 && periodicSubstrings[i].maxEndA>=Factors.factors[factor+1].firstPosition) firstPeriodicForNextFactor=i;
			if (!periodicSubstrings[i].shouldBeCountedForCoverage) continue;
			if ( Intervals.isApproximatelyContained_lowQuality(Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition,periodicSubstrings[i].minStartA,periodicSubstrings[i].maxEndA,ReadA.id) ||
				 Intervals.areApproximatelyIdentical_lowQuality(Factors.factors[factor].firstPosition,Factors.factors[factor].lastPosition,periodicSubstrings[i].minStartA,periodicSubstrings[i].maxEndA,ReadA.id)
			   ) {			
				Factors.factors[factor].periodicCoverage++;
				if (periodicSubstrings[i].containsExactOccurrence) Factors.factors[factor].periodicCoverageIdentity++;
				if (periodicSubstrings[i].equalsReverseComplement) Factors.factors[factor].equalsReverseComplement[0]=1;
			}
		}
		if (i>lastSubstring) {
			if (firstPeriodicForNextFactor!=-1) firstPeriodic=firstPeriodicForNextFactor;
			else firstPeriodic=i;
		}
		out[0]=firstPeriodic; out[1]=i-1;
	}

	
	/**
	 * Let $factor$ be a factor of $readA$. Given an alignment $X$ of $readA$, with a
	 * given $readB$, that intersects $factor$ in $readA$, call $f_X$ (respectively,
	 * $g_X$) the first position (respectively, the last position) of the alignment in
	 * $readB$ that aligns with the intersection between the alignment and $factor$ in
	 * $readA$. The procedure builds the DAG of all alignments in
	 * $ReadA.sortedAlignments[blockStart..blockEnd]$, where $ReadA.sortedAlignments$ is
	 * assumed to be sorted by $shouldBeUsed$, $readB$, and orientation, then by $f_X$,
	 * and where the interval $[blockStart..blockEnd]$ is assumed to correspond to the set
	 * of all alignments between $factor$ and $readB$, in a specific orientation.
	 *
	 * Given two alignments $X$ and $Y$, we say that $X$ precedes $Y$ iff $f_X$ is smaller
	 * than $f_Y$. Such a pair $(X,Y)$ receives a 0/1 cost as computed by procedure 
	 * $getCost$. Let $X(Y)$ be the leftmost alignment that precedes $Y$ such that $(X,Y)$
	 * has cost zero. The procedure builds the DAG in which every $Y$ is connected to all
	 * alignments $X$ such that $f_X \in [f_{X(Y)}..f_Y)$, assigning to each arc its
	 * corresponding cost. If $Y$ has no $X(Y)$, then it is connected just to the
	 * rightmost preceding $X$ such that $(X,Y)$ has cost one.
	 *
	 * Consider the shortest distance, in the DAG, between the alignment at position
	 * $blockStart$ and every other alignment in $[blockStart..blockEnd]$ up to the last
	 * usable alignment. Such distances are monotonically nondecreasing. Let a run be a
 	 * maximal interval of alignments with the same distance. The surface of a run is the
	 * number of positions in $factor$ that are covered by the alignments in the run. The
	 * procedure stores all runs in array $occurrences$.
	 *
	 * The procedure also computes the surface of a run that is covered by identity
	 * alignments, i.e. by alignments that are mostly inside $factor$, and that have a
	 * relative $diffs$ at most equal to $factor.maxDiffs$. If such surface is at least
	 * $exactSurfaceThreshold$ times the length of $factor$, the occurrence associated
	 * with the run is assumed to be an exact copy of $factor$. If $factor$ is periodic,
	 * just one identity alignment is enough to consider the occurrence an exact copy of
	 * $factor$: this is because we assume $factor$ to contain approximately homogeneous
	 * periods.
	 *
	 * Remark: the procedure alters the order of $ReadA.sortedAlignments[blockStart..
	 * blockEnd]$.
	 *
	 * Remark: a specific instance of a factor in readB could be built from the full
	 * factor in more than one way, where every such way corresponds to a distinct
	 * subsequence of the factor. Moreover, a \emph{substring} of an instance could be
	 * built using yet other subsequences of the factor. The aligner could report the
	 * alignments that correspond to all such possibilities. The shortest path method used
	 * by this procedure could assign the same shortest distance to a number of such
	 * distinct subsequences, effectively collapsing them into a single one. However, it
	 * could also assign different shortest distances to different subsequences, thus
	 * creating a number of objects in $occurrences$, possibly of different replication
	 * type, that all overlap in readB. This implies that the value of the shortest
	 * distance is not necessarily monotonic on the alignments sorted by $f$. Example:
	 *
	 * $readB[a..b]=factor[0..c]=factor[d..e]*factor[f..g]=factor[h..i]*factor[j..k]$.
	 * Let $x'$ be the projection of $x$ on readB. Let $f' < j'$, $d' > a$, $h' > d'$. 
	 * Then $readB[d'..e']*readB[f'..g']$ is an occurrence at distance 1 from the first
	 * alignment, $readB[h'..i']*readB[j'..k']$ is an occurrence at distance 2 from the
	 * first alignment, and distance is not monotonic on the order $[a..b],[d'..e'],
	 * [h'..i'],[f'..g'],[j'..k']$.
	 */
	private static final void estimateOccurrences(int blockStart, int blockEnd, Factor factor) {
		boolean orientation, containsExactAlignment;
		boolean hasDeletion, hasLeftMaximalAlignment, hasRightMaximalAlignment;
		int i, j, firstUsable, lastUsable, lastZero;
		int nVertices, readB, shiftA, shiftB, firstA, lastA, maxLastA, distance;
		int runFirst, runLast;
		int quantum;
		int f, firstF, g, lastG, begin, end;
		int surface, exactSurface;
		double dilation, maxDilation, diffs;

		readB=Alignments.alignments[ReadA.sortedAlignments[blockStart].id][1]-1;
		orientation=Alignments.alignments[ReadA.sortedAlignments[blockStart].id][2]==1;
		quantum=Math.min(IO.quantum,(int)(factor.length*factorUnitFraction));  // Unit of length in readA. Since factors can be short, a fixed constant might not work well.
		nVertices=blockEnd-blockStart+1;
		firstUsable=0;
		lastUsable=-1;
		for (i=blockStart; i<=blockEnd; i++) {
			if (!ReadA.sortedAlignments[i].shouldBeUsed) break;
			lastUsable=i-blockStart;
		}
		if (lastUsable==-1) return;
if (IO.SHOW_STD_ERR) IO.printErr("ENO>>> number of vertices in the DAG: "+(lastUsable-firstUsable+1)+" lastUsable="+ReadA.sortedAlignments[blockStart+lastUsable].id+" (in LAshow)");

		// Finding the maximum dilation of an alignment in the block (can be <1).
		maxDilation=0;
		for (i=blockStart; i<=blockEnd; i++) {
			dilation=((double)(Alignments.alignments[ReadA.sortedAlignments[i].id][6]-Alignments.alignments[ReadA.sortedAlignments[i].id][5])) / (Alignments.alignments[ReadA.sortedAlignments[i].id][4]-Alignments.alignments[ReadA.sortedAlignments[i].id][3]);
			if (dilation>maxDilation) maxDilation=dilation;
		}

		// Computing the length of a shortest path from $firstUsable$ to all other
		// vertices up to $lastUsable$.
		for (i=0; i<nVertices; i++) ReadA.sortedAlignments[blockStart+i].distance=Math.POSITIVE_INFINITY;
		ReadA.sortedAlignments[blockStart+firstUsable].distance=0;
		for (j=firstUsable+1; j<=lastUsable; j++) {
			lastZero=-1;
			for (i=j-1; i>=firstUsable; i--) {
				if (getCost(ReadA.sortedAlignments[blockStart+i].id,ReadA.sortedAlignments[blockStart+j].id,readB,orientation,factor,maxDilation,quantum)==0) {
					hasZeroCost[i]=true;
					lastZero=i;
				}
				else hasZeroCost[i]=false;
			}
			if (lastZero==-1) {
				i=j-1;
				if (ReadA.sortedAlignments[blockStart+i].distance!=Math.POSITIVE_INFINITY) {
					distance=ReadA.sortedAlignments[blockStart+i].distance+1;
					if (distance<ReadA.sortedAlignments[blockStart+j].distance) ReadA.sortedAlignments[blockStart+j].distance=distance;
				}
			}
			else {
				for (i=j-1; i>=lastZero; i--) {
					if (ReadA.sortedAlignments[blockStart+i].distance==Math.POSITIVE_INFINITY) continue;
					distance=ReadA.sortedAlignments[blockStart+i].distance+(hasZeroCost[i]?0:1);
					if (distance<ReadA.sortedAlignments[blockStart+j].distance) ReadA.sortedAlignments[blockStart+j].distance=distance;
				}
			}
		}
		if (ReadA.sortedAlignments[blockStart+lastUsable].distance==Math.POSITIVE_INFINITY) {
			IO.printCriticalErr("Error: no path from the first usable alignment to the last usable alignment.");
			System.exit(1);
		}

		// Sorting alignments between $firstUsable$ and $lastUsable$ by distance, then by
		// $startA$. This is necessary, since the same instance of $factor$ in readB could
		// have been built in distinct ways from $factor$, and these ways can make
		// distance non-monotonic in the original order by $f$.
		// The same distance could be assigned to distinct ways of building the instance
		// from $factor$: in this case, they will be collapsed into a single element of
		// $occurrences$.
		if (Alignment.order!=Alignment.DISTANCE_STARTA) {
			Alignment.order=Alignment.DISTANCE_STARTA;
			if (lastUsable+1-firstUsable>1) Arrays.sort(ReadA.sortedAlignments,blockStart+firstUsable,blockStart+lastUsable+1);
		}

		// Splitting equal-distance runs every time two consecutive alignments have cost
		// one according to $getCost$: this happens only when distinct occurrences have
		// been assigned the same distance. The result, however, could oversplit the true
		// occurrences.
		// Computing the fraction of $factor$ covered by the alignments in each run, and
		// the fraction of such covered surface that is covered by the identity alignments
		// in each run.
		runFirst=firstUsable; runLast=runFirst;
		for (i=firstUsable+1; i<=lastUsable+1; i++) {
			if ( i<=lastUsable &&
			     ReadA.sortedAlignments[blockStart+i].distance==ReadA.sortedAlignments[blockStart+i-1].distance &&
			     getCost(ReadA.sortedAlignments[blockStart+i-1].id,ReadA.sortedAlignments[blockStart+i].id,readB,orientation,factor,maxDilation,quantum)==0
			   ) {
				runLast=i;
				continue;
			}
			firstA=Alignments.alignments[ReadA.sortedAlignments[blockStart+runFirst].id][3];
			if (firstA<factor.firstPosition) firstA=factor.firstPosition;
			begin=firstA;
			maxLastA=Alignments.alignments[ReadA.sortedAlignments[blockStart+runFirst].id][4];
			if (maxLastA>factor.lastPosition) maxLastA=factor.lastPosition;
			surface=maxLastA-firstA+1;
			firstF=Alignments.getF(ReadA.sortedAlignments[blockStart+runFirst].id,factor.firstPosition);
			lastG=Alignments.getG(ReadA.sortedAlignments[blockStart+runFirst].id,factor.lastPosition);
			hasDeletion=false;
			for (j=runFirst+1; j<=runLast; j++) {
				lastA=Alignments.alignments[ReadA.sortedAlignments[blockStart+j].id][4];
				if (lastA>factor.lastPosition) lastA=factor.lastPosition;
				if (lastA>maxLastA) {
					firstA=Alignments.alignments[ReadA.sortedAlignments[blockStart+j].id][3];
					if (firstA<factor.firstPosition) firstA=factor.firstPosition;
					if (firstA<=maxLastA) surface+=lastA-maxLastA;
					else {
						surface+=lastA-firstA+1;
						if (firstA-maxLastA-1>=DenseSubstrings.minDeletionLength) hasDeletion=true;
					}
					maxLastA=lastA;
				}
				f=Alignments.getF(ReadA.sortedAlignments[blockStart+j].id,factor.firstPosition);
				if (f<firstF) firstF=f;
				g=Alignments.getG(ReadA.sortedAlignments[blockStart+j].id,factor.lastPosition);
				if (g>lastG) lastG=g;
			}
			end=maxLastA;
			hasLeftMaximalAlignment=false; hasRightMaximalAlignment=false;
			for (j=runFirst; j<=runLast; j++) {
				if ( ReadA.sortedAlignments[blockStart+j].isLeftMaximal==1 &&
					 Math.abs(Alignments.getF(ReadA.sortedAlignments[blockStart+j].id,factor.firstPosition),firstF)<=quantum*maxDilation
				   ) hasLeftMaximalAlignment=true;
				if ( ReadA.sortedAlignments[blockStart+j].isRightMaximal==1 &&
					 Math.abs(Alignments.getG(ReadA.sortedAlignments[blockStart+j].id,factor.lastPosition),lastG)<=quantum*maxDilation
				   ) hasRightMaximalAlignment=true;
			}
			exactSurface=0; maxLastA=-1;
			for (j=runFirst; j<=runLast; j++) {
				firstA=Alignments.alignments[ReadA.sortedAlignments[blockStart+j].id][3];
				lastA=Alignments.alignments[ReadA.sortedAlignments[blockStart+j].id][4];
				if (!Intervals.isApproximatelyContained(firstA,lastA,factor.firstPosition,factor.lastPosition)) continue;
				diffs=Alignments.getAvgDiffs(ReadA.sortedAlignments[blockStart+j].id);
				if (diffs>factor.maxDiffs) continue;
				if (lastA>factor.lastPosition) lastA=factor.lastPosition;
				if (lastA>maxLastA) {
					if (firstA<factor.firstPosition) firstA=factor.firstPosition;
					if (firstA<=maxLastA) exactSurface+=lastA-maxLastA;
					else exactSurface+=lastA-firstA+1;
					maxLastA=lastA;
				}
			}
			lastOccurrence++;
			occurrences[lastOccurrence].orientation=orientation;
			occurrences[lastOccurrence].firstPositionInForwardOrientation=orientation?firstF:Reads.getReadLength(readB)-lastG-1;
			occurrences[lastOccurrence].lastPositionInForwardOrientation=orientation?lastG:Reads.getReadLength(readB)-firstF-1;
			if (factor.periodicCoverage>0) {
				occurrences[lastOccurrence].isExact=exactSurface>0;
				occurrences[lastOccurrence].type=-1;
			}
			else {
				occurrences[lastOccurrence].isExact=exactSurface>=surface*exactSurfaceThreshold;
				// More than one way of building the instance of $factor$ in readB could
				// be collapsed into the same element of $occurrences$. Assigning the
				// simplest explanation. If the instance in readB is both a prefix and a
				// suffix, prefix is arbitrarily preferred.
				if (surface>=surfaceThreshold*(factor.lastPosition-factor.firstPosition+1)) occurrences[lastOccurrence].type=0;
				else if (begin-factor.firstPosition<=quantum && factor.lastPosition-end>quantum && hasRightMaximalAlignment) occurrences[lastOccurrence].type=1;
				else if (begin-factor.firstPosition>quantum && factor.lastPosition-end<=quantum && hasLeftMaximalAlignment) occurrences[lastOccurrence].type=2;
				else if (begin-factor.firstPosition>quantum && hasLeftMaximalAlignment && factor.lastPosition-end>quantum && hasRightMaximalAlignment) occurrences[lastOccurrence].type=3;
				else if (hasDeletion) occurrences[lastOccurrence].type=4;
				else occurrences[lastOccurrence].type=5;
			}
			// We don't need to sort back, since the alignments in a run will not be used
			// again during the processing of $[blockStart..blockEnd]$, and since
			// procedure $estimateCorrectedCoverage$ already restores the correct order
			// before moving to the next factor.
			runFirst=i; runLast=runFirst;
		}
	}
	
	
	/**
	 * Computes the cost of arc $(i,j)$ in the DAG built by procedure
	 * $estimateOccurrences$, where $i$ precedes $j$ according to that procedure.
	 *
	 * 0. If intervals $(f_i,g_i)$ and $(f_j,g_j)$ in readB have a Jaccard similarity of
	 * at least $overlapThreshold$, or if one is contained in the other, and if intervals
	 * $i$ and $j$ overlap in readA, pair $(i,j)$ has cost zero.
	 * 1. Otherwise, if the order between $f_i$ and $f_j$ violates the order of the
	 * corresponding points in readA, pair $(i,j)$ has cost one.
	 * 2. Otherwise, if the offset between $f_i$ and $f_j$ is similar to the offset
	 * between the first position of $i$ and the first position of $j$, and if $i$ and $j$
	 * overlap, the cost is zero.
	 * 3. Otherwise, if the offset between the last position of $i$ and the first position
	 * of $j$ in readA is small, and if the substring of readB between $g_i$ and $f_j$ is
	 * a random insertion, the cost is zero.
	 * 4. Otherwise, if the offset between $g_i$ and $f_j$ is small, and the offset
	 * between the last position of $i$ and the first position of $j$ is large (internal
	 * deletion), the cost is zero.
	 * 5. Otherwise, the cost is one.
	 *
	 * Remark: condition (0) handles the case of a periodic substring in readA that
	 * occurs also in readB, possibly with a different length, but such that the pair of
	 * intervals in readA and readB is not detected as part of a repeated periodic
	 * substring. Because of the periodicity, alignments that start consecutively in readB
	 * could violate the order in readA, thus leading to an over-estimation of the number
	 * of occurrences of the factor in readB.
	 *
	 * Remark: condition (0) handles also the case of an alignment that is contained in
	 * another one.
	 *
	 * Remark: alignments that fully contain the factor are handled correctly (both pairs
	 * of such alignments, and pairs in which just one alignment is fully containing).
	 *
	 * @param $idi$ (respectively, $idj$) is the id of alignment $i$ (respectively, $j$)
	 * in $Alignments.alignments$;
	 * @param maxDilation maximum dilation observed in an alignment between $readA$ and
	 * $readB$ in the block of $ReadA.sortedAlignments$ that starts at $blockStart$ (can
	 * be smaller than 1). The distance between $f_i$ and $f_j$ is considered "too large" 
	 * if the ratio between this distance, and the corresponding distance in $readA$, 
	 * exceeds the maximum dilation ratio observed in all alignments of $readA$ with 
	 * $readB$.
	 * @param quantum distance threshold to consider the intersections of the two
	 * alignments with $factor$ in readA as overlapping;
	 * @return 0 or 1.
	 */
	private static final int getCost(int idi, int idj, int readB, boolean orientation, Factor factor, double maxDilation, int quantum) {
		boolean overlapA;
		int shiftA, shiftB;
		int startAI, endAI, startAJ, endAJ;
		int fi, gi, fj, gj;

		startAI=Alignments.alignments[idi][3];
		endAI=Alignments.alignments[idi][4];
		startAJ=Alignments.alignments[idj][3];
		endAJ=Alignments.alignments[idj][4];
		fi=Alignments.getF(idi,factor.firstPosition);
		gi=Alignments.getG(idi,factor.lastPosition);
		fj=Alignments.getF(idj,factor.firstPosition);
		gj=Alignments.getG(idj,factor.lastPosition);

		// Case 0
		if ( Intervals.intersectionLength(startAI,endAI,startAJ,endAJ)>0 &&
		     (Intervals.jaccardSimilarity(fi,gi,fj,gj)>=overlapThreshold || Intervals.eitherIsContained(fi,gi,fj,gj))
		   ) return 0;

		// Case 1
		shiftA=Alignments.getAShift(idi,idj,factor.firstPosition);
		if (shiftA<-quantum) return 1;

		// Case 2
		shiftB=Alignments.getBShift(idi,idj,factor.firstPosition);  // Never negative
		if ( (shiftA>0 && shiftB<=shiftA*maxDilation && Intervals.intersectionLength(startAI,endAI,startAJ,endAJ)>0) ||
			 (shiftA==0 && shiftB<=maxDilation*quantum) 
		   ) return 0;

		// Cases 3 and 4
		if (startAJ<factor.firstPosition) startAJ=factor.firstPosition;
		if (endAI>factor.lastPosition) endAI=factor.lastPosition;
		if ( Math.abs(startAJ,endAI)<=quantum &&
			 fj-gi-1>=Reads.minLongInsertionLength &&
			 Reads.isRandomInsertion(readB,gi+1,fj-1,orientation) ) return 0;
		if ( startAJ-endAI>=DenseSubstrings.minDeletionLength &&
		     fj-gi-1<=quantum*maxDilation ) return 0;

		return 1;
	}
	
	
	/**
	 * Estimates the total number of occurrences of $factor$ in $readB$, in both the
	 * forward and in the reverse-complement orientation, using array $occurrences$.
	 *
	 * Two elements in $occurrences$ can overlap significantly and have the same
	 * orientation (for example, when $factor$ is contained inside two alignments, and
	 * slightly different estimates of the position of $factor$ in $readB$ are computed
	 * from each alignment). Two elements in $occurrences$ can overlap significantly and
	 * have opposite orientation if $factor$ equals its reverse complement: in this case,
	 * the procedure sets the cells of $factor.equalsReverseComplement$, depending on
	 * whether it detects a full+full overlap, a prefix+suffix overlap, or a
	 * substring+substring overlap. The procedure scans the sequence of elements in
	 * $occurrences$, growing a run as long as a new element overlaps significantly with
	 * the current run, and creating a new run otherwise. The number of runs is reported
	 * in output.
	 *
	 * Remark: two occurrences could overlap significantly because the factor has a long
	 * border, i.e. a short period. If the occurrences are at least $minAlignmentLength$,
	 * this implies that the readB instance of the factor should belong to a periodic
	 * substring (but it doesn't). If the occurrences are shorter than
	 * $minAlignmentLength$, then it could be that the instance is indeed a periodic
	 * substring, and that it has not been detected as such before. For example, let $UVW$
	 * in readB be such that $V$ is periodic, $|V|$ is smaller than $minAlignmentLength$, 
	 * $|UV| \geq minAlignmentLength$, $|VW| \geq minAlignmentLength$. Let readA contain 
	 * string $UV'W$, where $|V'|<|V|$ and $V'$ has the same period as $V$. Then, the 
	 * aligner could have aligned the $UV'$ substring and the $V'W$ substring in the two 
	 * reads, and such aligned $V'$ substrings could overlap in readB. We don't try to 
	 * detect such events here, both for simplicity, and because the boundaries of an 
	 * occurrence in readB are possibly inaccurate estimates.
	 *
	 * If $factor$ is periodic, the procedure concatenates runs that have, at their left
	 * or at their right, another run, or a suitable other element of
	 * $PeriodicSubstrings.substrings$, such that the intervening space in $readB$ is
	 * smaller than $IO.quantum$, or it is a random insertion. The result of a maximal
	 * sequence of concatenations contributes just one to the coverage. See $Rope$ for a
	 * description of how threads and occurrences are concatenated. See procedure
	 * $PeriodicSubstrings.correctCoverage()$ for a similar process.
	 *
	 * Remark: the concatenation procedure described above works also if an occurrence in
	 * readB is itself fragmented by random insertions.
	 *
	 * Remark: the procedure assumes $PeriodicSubstrings.substrings$ to be sorted by
	 * $readB$ and $minStartBForward$.
	 *
	 * Remark: the procedure writes the global variables $deltaNonperiodicCoverage$,
	 * $deltaNonperiodicCoverageIdentity$, $deltaPeriodicCoverage$,
	 * $deltaPeriodicCoverageIdentity$, but it doesn't update the corresponding values of
	 * $factor$.
	 *
	 * Remark: at this point of the pipeline, a factor is unlikely to \emph{contain} a
	 * maximal substring that equals the reverse-complement of itself. To be detected by
	 * an alignment, such substring must either be fully inside the factor and at least
	 * $IO.minAlignmentLength$ positions long, or it must be a prefix or a suffix of the
	 * factor. In both cases, it is easy to prove that such substring is a maximal repeat
	 * of the read set, so it should have been detected as a distinct factor by the
	 * previous steps of the pipeline.
	 *
	 * Remark: two elements in $occurrences$ can overlap significantly and have the same
	 * orientation, also if they describe alternative ways of obtaining the same instance
	 * of $factor$ in readB (or a substring of such instance) by deletions. All such
	 * overlapping elements of $occurrences$ are collapsed into a single instance, whose
	 * type is assigned according to the simplest explanation (see $OccurrenceTypes$).
	 *
	 * @param firstPeriodicSubstring the first element in $PeriodicSubstrings.substrings$
	 * related to $readB$.
	 */
	private static final void estimateNOccurrences(Factor factor, int readB, int firstPeriodicSubstring) {
		boolean isExact, previousIsExact, mergeLeft;
		int i, j, k;
		int firstI, previousFirstI, firstJ, firstJForNextOccurrence;
		int firstPosition, lastPosition, previousFirstPosition, previousLastPosition, readBLength;
		int firstStart, lastEnd, periodicSubstringStart, periodicSubstringEnd;
		int type;
		PeriodicSubstring firstSubstring, lastSubstring;
		Rope previousRope, currentRope;
		Rope tmpRope1 = new Rope();
		Rope tmpRope2 = new Rope();
		Math.set(deltaNonperiodicCoverage,deltaNonperiodicCoverage.length-1,0);
		Math.set(deltaNonperiodicCoverageIdentity,deltaNonperiodicCoverageIdentity.length-1,0);
		deltaPeriodicCoverage=0; deltaPeriodicCoverageIdentity=0;
		if (lastOccurrence==-1) return;

		// Ensuring the necessary order in $occurrences$
		if (lastOccurrence>0) {
			Occurrence.order=Occurrence.FIRSTPOSITION;
			if (lastOccurrence>0) Arrays.sort(occurrences,0,lastOccurrence+1);
		}

		// Nonperiodic factor
		if (factor.periodicCoverage==0) {
			firstPosition=occurrences[0].firstPositionInForwardOrientation;
			lastPosition=occurrences[0].lastPositionInForwardOrientation;
			OccurrenceTypes.clear();
			firstI=0;
			for (i=1; i<=lastOccurrence+1; i++) {
				if ( i==lastOccurrence+1 ||
				     occurrences[i].firstPositionInForwardOrientation>lastPosition ||
					 ( !Intervals.eitherIsContained(occurrences[i].firstPositionInForwardOrientation,occurrences[i].lastPositionInForwardOrientation,firstPosition,lastPosition) &&
					   Intervals.jaccardSimilarity(occurrences[i].firstPositionInForwardOrientation,occurrences[i].lastPositionInForwardOrientation,firstPosition,lastPosition)<Intervals.jaccardThreshold
					 )
				   ) {
				   	// Using just occurrences that span approximately the whole interval
				   	// $[firstPosition..lastPosition]$ in readB to fill matrix $types$.
				   	// This is necessary, since there could be occurrences that are just
				   	// substrings of $[firstPosition..lastPosition]$ in readB, e.g. if
				   	// the factor contains repeated substrings.
					for (j=firstI; j<i; j++) {
						if (Intervals.jaccardSimilarity(occurrences[j].firstPositionInForwardOrientation,occurrences[j].lastPositionInForwardOrientation,firstPosition,lastPosition)>=Intervals.jaccardThreshold) OccurrenceTypes.add(occurrences[j]);
					}
				   	type=OccurrenceTypes.getType();
				   	if (type>0) {
				   		deltaNonperiodicCoverage[type-1]++;
				   		deltaNonperiodicCoverageIdentity[type-1]++;
				   		OccurrenceTypes.setEqualsReverseComplement(type-1,factor);
				   	}
				   	else {
				   		deltaNonperiodicCoverage[-type-1]++;
				   		OccurrenceTypes.setEqualsReverseComplement(-type-1,factor);
				   	}
				   	OccurrenceTypes.setHasSelfSimilarity(factor);
				   	OccurrenceTypes.clear();
					if (i<=lastOccurrence) {
						firstPosition=occurrences[i].firstPositionInForwardOrientation;
						lastPosition=occurrences[i].lastPositionInForwardOrientation;
						firstI=i;
					}
					continue;
				}
				if (occurrences[i].lastPositionInForwardOrientation>lastPosition) lastPosition=occurrences[i].lastPositionInForwardOrientation;
			}
			return;
		}

		// Periodic factor

if (IO.SHOW_STD_ERR) IO.printErr("ENO> n. occurrences = "+(lastOccurrence+1));
if (IO.SHOW_STD_ERR) IO.printErr("ENO> firstPeriodicSubstring="+PeriodicSubstrings.substrings[firstPeriodicSubstring]);


		readBLength=Reads.getReadLength(readB);
		firstPosition=occurrences[0].firstPositionInForwardOrientation;
		lastPosition=occurrences[0].lastPositionInForwardOrientation;
		previousIsExact=false; isExact=occurrences[0].isExact;
		previousLastPosition=-1; previousFirstPosition=-1;
		firstI=0; previousFirstI=-1;
		firstJ=firstPeriodicSubstring; firstJForNextOccurrence=-1;
		previousRope=null; currentRope=null;
		for (i=1; i<=lastOccurrence+1; i++) {
			if ( i==lastOccurrence+1 ||
			     occurrences[i].firstPositionInForwardOrientation>lastPosition ||
			     ( !Intervals.eitherIsContained(occurrences[i].firstPositionInForwardOrientation,occurrences[i].lastPositionInForwardOrientation,firstPosition,lastPosition) &&
			       Intervals.jaccardSimilarity(occurrences[i].firstPositionInForwardOrientation,occurrences[i].lastPositionInForwardOrientation,firstPosition,lastPosition)<Intervals.jaccardThreshold
			     )
			   ) {

if (IO.SHOW_STD_ERR) IO.printErr("ENO> considering run of occurrences ["+firstPosition+".."+lastPosition+"] firstJ="+firstJ+" firstPeriodicSubstring="+firstPeriodicSubstring+" previousFirstPosition="+previousFirstPosition+" previousLastPosition="+previousLastPosition+" previousFirstI="+previousFirstI);

				lastEnd=-1; firstStart=-1;
				firstSubstring=null; lastSubstring=null;
				mergeLeft=false;
				for (j=firstJ; j<=PeriodicSubstrings.lastSubstring; j++) {
					if (PeriodicSubstrings.substrings[j].readB!=readB) break;

if (IO.SHOW_STD_ERR) IO.printErr("ENO> considering periodic substring "+j+"...");

					if ( Intervals.jaccardSimilarity(PeriodicSubstrings.substrings[j].minStartA,PeriodicSubstrings.substrings[j].maxEndA,factor.firstPosition,factor.lastPosition)<Intervals.jaccardThreshold &&
					     !Intervals.eitherIsContained(PeriodicSubstrings.substrings[j].minStartA,PeriodicSubstrings.substrings[j].maxEndA,factor.firstPosition,factor.lastPosition)
					   ) continue;
					periodicSubstringStart=PeriodicSubstrings.substrings[j].minStartBForward;
					periodicSubstringEnd=PeriodicSubstrings.substrings[j].maxEndBForward;

if (IO.SHOW_STD_ERR) IO.printErr("ENO> ok, using it: ["+periodicSubstringStart+".."+periodicSubstringEnd+"]");

					if (periodicSubstringEnd>=firstPosition) {
						if (periodicSubstringStart<=firstPosition) {
							lastEnd=periodicSubstringEnd;
							lastSubstring=PeriodicSubstrings.substrings[j];
							if (previousLastPosition!=-1 && firstStart==-1) {
								firstStart=periodicSubstringStart;
								firstSubstring=PeriodicSubstrings.substrings[j];
							}
						}
						if ( ( Intervals.jaccardSimilarity(periodicSubstringStart,periodicSubstringEnd,firstPosition,lastPosition)>=Intervals.jaccardThreshold ||
						       Intervals.eitherIsContained(periodicSubstringStart,periodicSubstringEnd,firstPosition,lastPosition)
						     ) &&
						     PeriodicSubstrings.substrings[j].orientation!=occurrences[firstI].orientation 
						   ) factor.equalsReverseComplement[0]=1;
						firstJForNextOccurrence=j;
						break;
					}
					if (periodicSubstringStart<=previousLastPosition) {
						if (periodicSubstringEnd>=previousLastPosition) {
							if (firstStart==-1) {
								firstStart=periodicSubstringStart;
								firstSubstring=PeriodicSubstrings.substrings[j];
							}
							if (periodicSubstringEnd>lastEnd) {
								lastEnd=periodicSubstringEnd;
								lastSubstring=PeriodicSubstrings.substrings[j];
							}
						}
						if ( ( Intervals.jaccardSimilarity(periodicSubstringStart,periodicSubstringEnd,previousFirstPosition,previousLastPosition)>=Intervals.jaccardThreshold ||
						       Intervals.eitherIsContained(periodicSubstringStart,periodicSubstringEnd,previousFirstPosition,previousLastPosition)
						     ) &&
						     (previousFirstI>=0 && PeriodicSubstrings.substrings[j].orientation!=occurrences[previousFirstI].orientation)
						   ) factor.equalsReverseComplement[0]=1;
						continue;
					}
					if (previousLastPosition!=-1 && firstStart==-1) {
						firstStart=periodicSubstringStart;
						firstSubstring=PeriodicSubstrings.substrings[j];
					}
					if (periodicSubstringEnd>lastEnd) {
						lastEnd=periodicSubstringEnd;
						lastSubstring=PeriodicSubstrings.substrings[j];
					}
				}


if (IO.SHOW_STD_ERR) IO.printErr("ENO> previousRope!=null?"+(previousRope!=null)+" firstStart="+firstStart+" previousLastPosition="+previousLastPosition+" lastEnd="+lastEnd);
if (previousLastPosition!=-1) if (IO.SHOW_STD_ERR) IO.printErr("ENO> isRandomInsertion? "+Reads.isRandomInsertion(readB,previousLastPosition+1,firstStart-1,true));



				if ( firstStart!=-1 &&
					 ( firstStart-previousLastPosition-1<=IO.quantum ||
					   Reads.isRandomInsertion(readB,previousLastPosition+1,firstStart-1,true)
					 )
				   ) {
				   	mergeLeft=true;
					if (previousRope==null) {
						previousRope=tmpRope1;
						previousRope.initialize(previousIsExact,firstSubstring);
						updateDeltas();
					}
					else {
						previousRope.concatenate(firstSubstring);
						updateDeltas();
					}
				}
				if (lastEnd==-1) {
					if ( previousLastPosition==-1 ||
					     (firstPosition>previousLastPosition && !Reads.isRandomInsertion(readB,previousLastPosition+1,firstPosition-1,true))
					   ) {
						currentRope=null;
						deltaNonperiodicCoverage[0]++;
						if (isExact) deltaNonperiodicCoverageIdentity[0]++;
					}
					else {
						if (previousRope!=null) {
							previousRope.concatenateOccurrence(isExact);
							currentRope=tmpRope2;
							currentRope.clone(previousRope);
							updateDeltas();
						}
						else {
							currentRope=tmpRope2;
							currentRope.initialize(previousIsExact,isExact);
							updateDeltas();
						}
					}
				}
				else {
					if ( firstPosition-lastEnd-1<=IO.quantum ||
						 Reads.isRandomInsertion(readB,lastEnd+1,firstPosition-1,true) ) {
if (IO.SHOW_STD_ERR) IO.printErr("-->1: firstPosition="+firstPosition+" lastEnd="+lastEnd);
						if (mergeLeft && lastSubstring.threadStart==firstSubstring.threadStart) {
							previousRope.concatenateOccurrence(isExact);
							currentRope=tmpRope2;
							currentRope.clone(previousRope);
						}
						else {
							currentRope=tmpRope2;
							currentRope.initialize(lastSubstring,isExact);
						}
						updateDeltas();
					}
					else {
if (IO.SHOW_STD_ERR) IO.printErr("-->2");
						currentRope=null;
						deltaNonperiodicCoverage[0]++;
						if (isExact) deltaNonperiodicCoverageIdentity[0]++;
					}
				}

				// Updating variables for next iteration
				previousLastPosition=lastPosition;
				previousFirstPosition=firstPosition;
				previousFirstI=firstI;
				previousIsExact=isExact;
				if (currentRope==null) previousRope=null;
				else {
					if (previousRope==null) previousRope=tmpRope1;
					previousRope.clone(currentRope);
				}
				firstJ=firstJForNextOccurrence!=-1?firstJForNextOccurrence:j;
				firstJForNextOccurrence=-1;
				if (i<=lastOccurrence) {
					firstPosition=occurrences[i].firstPositionInForwardOrientation;
					lastPosition=occurrences[i].lastPositionInForwardOrientation;
					isExact=occurrences[i].isExact;
					firstI=i;
				}
				continue;
			}
			if (occurrences[i].lastPositionInForwardOrientation>lastPosition) lastPosition=occurrences[i].lastPositionInForwardOrientation;
			if (occurrences[i].orientation!=occurrences[firstI].orientation) factor.equalsReverseComplement[0]=1;
			if (occurrences[i].isExact) isExact=true;
if (IO.SHOW_STD_ERR) IO.printErr("ENO> i incorporated, in current run, occurrence "+occurrences[i]);
		}


if (IO.SHOW_STD_ERR) IO.printErr("(before) deltaNonperiodicCoverage="+deltaNonperiodicCoverage+" deltaPeriodicCoverage="+deltaPeriodicCoverage);


		// Handling periodic substrings to the right of the last occurrence
		firstStart=-1; firstSubstring=null;
		for (j=firstJ; j<=PeriodicSubstrings.lastSubstring; j++) {
			if (PeriodicSubstrings.substrings[j].readB!=readB) break;
			if ( Intervals.jaccardSimilarity(PeriodicSubstrings.substrings[j].minStartA,PeriodicSubstrings.substrings[j].maxEndA,factor.firstPosition,factor.lastPosition)<Intervals.jaccardThreshold &&
				 !Intervals.eitherIsContained(PeriodicSubstrings.substrings[j].minStartA,PeriodicSubstrings.substrings[j].maxEndA,factor.firstPosition,factor.lastPosition)
			   ) continue;
			periodicSubstringStart=PeriodicSubstrings.substrings[j].minStartBForward;
			periodicSubstringEnd=PeriodicSubstrings.substrings[j].maxEndBForward;


if (IO.SHOW_STD_ERR) IO.printErr("ENO> considering periodic substring "+j+": ["+periodicSubstringStart+".."+periodicSubstringEnd+"]");


			if (periodicSubstringStart<=previousLastPosition) {
				if (periodicSubstringEnd>=previousLastPosition && firstStart==-1) {
					firstStart=periodicSubstringStart;
					firstSubstring=PeriodicSubstrings.substrings[j];
				}
				if ( ( Intervals.jaccardSimilarity(periodicSubstringStart,periodicSubstringEnd,previousFirstPosition,previousLastPosition)>=Intervals.jaccardThreshold ||
					   Intervals.eitherIsContained(periodicSubstringStart,periodicSubstringEnd,previousFirstPosition,previousLastPosition)
					 ) &&
					 PeriodicSubstrings.substrings[j].orientation!=occurrences[previousFirstI].orientation ) factor.equalsReverseComplement[0]=1;
				continue;
			}
			if (firstStart==-1) {
				firstStart=periodicSubstringStart;
				firstSubstring=PeriodicSubstrings.substrings[j];
				break;
			}
		}


if (IO.SHOW_STD_ERR) IO.printErr("ENO> firstStart="+firstStart+" previousLastPosition="+previousLastPosition);
if (previousLastPosition!=-1) if (IO.SHOW_STD_ERR) IO.printErr("ENO> isRandomInsertion? "+Reads.isRandomInsertion(readB,previousLastPosition+1,firstStart-1,true));


		if ( firstStart!=-1 &&
			 ( firstStart-previousLastPosition-1<=IO.quantum ||
			   Reads.isRandomInsertion(readB,previousLastPosition+1,firstStart-1,true)
			 )
		   ) {
			if (previousRope==null) {
				previousRope=tmpRope1;
				previousRope.initialize(previousIsExact,firstSubstring);
				updateDeltas();
			}
			else {
				previousRope.concatenate(firstSubstring);
				updateDeltas();
			}
		}


if (IO.SHOW_STD_ERR) IO.printErr("(after) deltaNonperiodicCoverage="+deltaNonperiodicCoverage+" deltaPeriodicCoverage="+deltaPeriodicCoverage);
	}
		
		
	/**
	 * Updates the global variables $delta*$ using the corresponding variables in $Rope$.
	 */
	private static final void updateDeltas() {
		deltaNonperiodicCoverage[0]+=Rope.deltaNonperiodicCoverage;
		deltaNonperiodicCoverageIdentity[0]+=Rope.deltaNonperiodicCoverageIdentity;
		deltaPeriodicCoverage+=Rope.deltaPeriodicCoverage;
		deltaPeriodicCoverageIdentity+=Rope.deltaPeriodicCoverageIdentity;
	}
	
	
	/**
	 * Assume that $PeriodicSubstrings.substrings[from..to]$ is sorted by $readB$. The
	 * procedure returns the smallest $i \in [from..to]$ such that
	 * $PeriodicSubstrings.substrings[i].readB=readB$. It returns $from$ even if no such
	 * $i$ can be found.
	 */
	private static final int readB2firstPeriodic(int readB, int from, int to) {
		if (PeriodicSubstrings.substrings[from].readB==readB) return from;
		for (int i=from+1; i<=to; i++) {
			if (PeriodicSubstrings.substrings[i].readB==readB) return i;
			if (PeriodicSubstrings.substrings[i].readB>readB) break;
		}
		return from;
	}


}