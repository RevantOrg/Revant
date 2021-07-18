package de.mpi_cbg.revant.factorize;


/**
 * The list of all the observed replication types of occurrences that are mapped to the
 * same instance of a factor in readB, broken down by orientation and quality.
 */
public class OccurrenceTypes {
	/**
	 * Rows: replication types (see $Factor$).
	 * Column 0: forward orientation.
	 * Column 1: forward orientation, exact.
	 * Column 2: reverse complement orientation.
	 * Column 3: reverse complement orientation, exact.
	 */
	private static boolean[][] types;


	public static final void allocateMemory() {
		types = new boolean[Factor.N_NONPERIODIC_COVERAGE_TYPES][4];
	}


	public static final void clear() {
		for (int i=0; i<types.length; i++) {
			for (int j=0; j<types[i].length; j++) types[i][j]=false;
		}
	}


	public static final void add(Occurrence occurrence) {
		if (occurrence.orientation) types[occurrence.type][occurrence.isExact?1:0]=true;
		else types[occurrence.type][occurrence.isExact?3:2]=true;
	}


	/**
	 * Decides the replication type of an instance of a factor in readB, using the types
	 * of all elements of $occurrences$ that are collapsed into the instance. The
	 * procedure returns the simplest explanation with identity-level quality. If no
	 * occurrence with identity-level quality exists, the procedure returns the simplest
	 * explanation. Prefix is arbitrarily preferred over suffix.
	 *
	 * @return X, where $|X|$=replication type+1, and $X$ is positive iff the instance is
	 * exact.
	 */
	public static final int getType() {
		int i;

		// Full
		if (types[0][1] || types[0][3]) return 1;
		else if (types[0][0] || types[0][2]) return -1;

		// Other types
		for (i=1; i<=5; i++) {
			if (types[i][1] || types[i][3]) return i+1;
		}
		for (i=1; i<=5; i++) {
			if (types[i][0] || types[i][2]) return -(i+1);
		}

		// All cells of $types$ are false
		return -6;
	}


	/**
	 * Sets to one some cells of array $factor.equalsReverseComplement$, using matrix
	 * $types$.
	 *
	 * @param type the replication type that has been assigned to the instance of $factor$
	 * in readB onto which all the occurrences that contributed to $types$ have been
	 * collapsed.
	 */
	public static final void setEqualsReverseComplement(int type, Factor factor) {
		if (type==0) {
			if (types[0][1] && types[0][3]) factor.equalsReverseComplement[0]=1;
			else if ( (types[0][0]&&types[0][3]) || (types[0][1]&&types[0][2]) ) factor.equalsReverseComplement[0]=2;
		}
		else if (type==1 || type==2) {
			if ( (types[1][1]&&types[2][3]) || (types[1][3]&&types[2][1]) ) factor.equalsReverseComplement[1]=1;
			else if ( (types[1][0]&&types[2][3]) || (types[1][1]&&types[2][2]) ) factor.equalsReverseComplement[2]=1;
		}
		else if (type==3) {
			if (types[3][1] && types[3][3]) factor.equalsReverseComplement[3]=1;
			else if ( (types[3][0]&&types[3][3]) || (types[3][1]&&types[3][2]) ) factor.equalsReverseComplement[4]=1;
		}
	}


	/**
	 * Sets $factor.hasSelfSimilarity$ using matrix $types$
	 */
	public static final void setHasSelfSimilarity(Factor factor) {
		boolean found;
		int i, j;

		found=false;
		for (i=0; i<=4; i++) {
			for (j=0; j<types[i].length; j++) {
				if (!types[i][j]) continue;
				if (!found) {
					found=true;
					break;
				}
				else {
					factor.hasSelfSimilarity=true;
					return;
				}
			}
		}
	}

}