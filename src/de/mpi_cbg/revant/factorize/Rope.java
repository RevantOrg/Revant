package de.mpi_cbg.revant.factorize;

/**
 * Let $[x..y]$ be a periodic factor of readA. Given a specific readB, a \emph{rope} is a
 * maximal sequence:
 *
 * $P_{1}^{n_1} O_{1}^{m_1} P_{2}^{n_2} O_{2}^{m_2} ... P_{k}^{n_k} O_{k}^{m_k} 
 * P_{k+1}^{n_{k+1}}$
 *
 * such that:
 *
 * 1. $P_{i}^{n_i}$ is a thread of $n_i$ intervals in readB of periodic substrings of
 * readA that have a significant overlap with $readA[x..y]$, as defined by procedure
 * $PeriodicSubstrings.correctCoverage$.
 *
 * 2. $O_{i}^{m_i}$ is a sequence of $m_i$ occurrences of $readA[x..y]$ in $readB$, as
 * defined by procedure $Factors.estimateNOccurrences$, such that every two consecutive
 * occurrences are either "too close", or they are separated by a long random insertion.
 *
 * 3. The last interval of $P_{i}^{n_i}$ and the first interval of $O_{i}^{m_i}$ are
 * either "too close", or they are separated by a long random insertion. The same holds
 * for the last interval of $O_{i}^{m_i}$ and the first interval of $P_{i+1}^{n_{i+1}}$.
 *
 * 4. The sequence contains at least one $O_{i}^{m_i}$, and at least two intervals.
 * The number of threads of periodic substrings can be zero.
 *
 * Maximality is in the sense of concatenating the rope with additional occurrences, or
 * with threads of periodic substrings, while still satisfying (1-4).
 *
 * A rope contributes one to the coverage of the factor. Specifically, a rope contributes
 * one to the \emph{periodic} coverage of the factor. Every occurrence that does not
 * belong to a rope contributes one to the nonperiodic coverage of the factor. Thus, a
 * rope can be seen as a generalization of a thread of periodic substrings, in which
 * distinct threads are concatenated using occurrences.
 *
 * Remark: a rope might not be extended to the right or to the left, if the intervening
 * substring of readB is not a random insertion. This could happen when such substring is
 * a sequence of random insertions and of periodic substrings shorter than
 * $Alignments.minAlignmentLength$. Thus, the number of ropes is an upper bound on the
 * number of periodic occurrences of a periodic factor.
 */
public class Rope {
	/**
	 * Properties of a rope
	 */
	private boolean containsExactOccurrence;
	private PeriodicSubstring lastThreadStart;  // Representative of the rightmost thread of periodic substrings that has been concatenated to the rope

	/**
	 * Deltas, reset by each procedure.
	 */
	public static int deltaNonperiodicCoverage, deltaNonperiodicCoverageIdentity;
	public static int deltaPeriodicCoverage, deltaPeriodicCoverageIdentity;


	public final void clone(Rope otherRope) {
		containsExactOccurrence=otherRope.containsExactOccurrence;
		lastThreadStart=otherRope.lastThreadStart;
	}


	/**
	 * Creates a new rope by concatenating a periodic substring to the right of an
	 * occurrence
	 */
	public final void initialize(boolean occurrenceIsExact, PeriodicSubstring substring) {
		deltaNonperiodicCoverage=0; deltaNonperiodicCoverageIdentity=0;
		deltaPeriodicCoverage=0; deltaPeriodicCoverageIdentity=0;

		// Computing deltas
		deltaNonperiodicCoverage--;
		if (occurrenceIsExact) deltaNonperiodicCoverageIdentity--;
		deltaPeriodicCoverage++;
		if (occurrenceIsExact||substring.threadStart.containsExactOccurrence) deltaPeriodicCoverageIdentity++;
		if (substring.threadStart.shouldBeCountedForCoverage) {
			deltaPeriodicCoverage--;
			if (substring.threadStart.containsExactOccurrence) deltaPeriodicCoverageIdentity--;
		}

		// Updating the rope
		containsExactOccurrence=occurrenceIsExact||substring.threadStart.containsExactOccurrence;
		lastThreadStart=substring.threadStart;
	}


	/**
	 * Creates a new rope by concatenating an occurrence to the right of a periodic
	 * substring
	 */
	public final void initialize(PeriodicSubstring substring, boolean occurrenceIsExact) {
		deltaNonperiodicCoverage=0; deltaNonperiodicCoverageIdentity=0;
		deltaPeriodicCoverage=0; deltaPeriodicCoverageIdentity=0;

		// Computing deltas
		if (substring.threadStart.shouldBeCountedForCoverage) {
			deltaPeriodicCoverage--;
			if (substring.threadStart.containsExactOccurrence) deltaPeriodicCoverageIdentity--;
		}
		deltaPeriodicCoverage++;
		if (occurrenceIsExact||substring.threadStart.containsExactOccurrence) deltaPeriodicCoverageIdentity++;

		// Updating the rope
		containsExactOccurrence=occurrenceIsExact||substring.threadStart.containsExactOccurrence;
		lastThreadStart=substring.threadStart;
	}


	/**
	 * Creates a new rope by concatenating an occurrence to the right of another occurrence
	 */
	public final void initialize(boolean occurrenceLeftIsExact, boolean occurrenceRightIsExact) {
		deltaNonperiodicCoverage=0; deltaNonperiodicCoverageIdentity=0;
		deltaPeriodicCoverage=0; deltaPeriodicCoverageIdentity=0;

		// Computing deltas
		deltaNonperiodicCoverage--;
		if (occurrenceLeftIsExact) deltaNonperiodicCoverageIdentity--;
		deltaPeriodicCoverage++;
		if (occurrenceLeftIsExact || occurrenceRightIsExact) deltaPeriodicCoverageIdentity++;

		// Updating the rope
		containsExactOccurrence=occurrenceLeftIsExact||occurrenceRightIsExact;
	}


	/**
	 * Concatenates a periodic substring to the right of the rope
	 */
	public final void concatenate(PeriodicSubstring substring) {
		deltaNonperiodicCoverage=0; deltaNonperiodicCoverageIdentity=0;
		deltaPeriodicCoverage=0; deltaPeriodicCoverageIdentity=0;

		// Computing deltas
		if (!containsExactOccurrence && substring.threadStart.containsExactOccurrence) deltaPeriodicCoverageIdentity++;
		if (substring.threadStart!=lastThreadStart && substring.threadStart.shouldBeCountedForCoverage) {
			deltaPeriodicCoverage--;
			if (substring.threadStart.containsExactOccurrence) deltaPeriodicCoverageIdentity--;
		}

		// Updating the rope
		if (substring.threadStart.containsExactOccurrence) containsExactOccurrence=true;
		lastThreadStart=substring.threadStart;
	}


	/**
	 * Concatenates an occurrence to the right of the rope
	 */
	public final void concatenateOccurrence(boolean occurrenceIsExact) {
		deltaNonperiodicCoverage=0; deltaNonperiodicCoverageIdentity=0;
		deltaPeriodicCoverage=0; deltaPeriodicCoverageIdentity=0;

		// Computing deltas
		if (!containsExactOccurrence && occurrenceIsExact) deltaPeriodicCoverageIdentity++;

		// Updating the rope
		if (occurrenceIsExact) containsExactOccurrence=true;
	}

}