package de.mpi_cbg.revant.factorize;

import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;


/**
 * Compares the period of every factor in $factorization.txt$ to its true period, based
 * on $trueFactorization.txt$. The boundaries of the factor are not checked.
 */
public class CheckPeriods {
	/**
	 * Parameters of the pipeline
	 */
	private static final double intersectionThreshold = 0.2;  // Fraction of a factor length to decide an intersection
	private static final double identityThreshold = 0.8;  // Fraction of the length of a true factor to decide its inclusion in $substrings$

	/**
	 * Reused variables
	 */
	private static int[][] repeats;
	private static int[][] factors;
	private static int[][] trueFactors;
	private static final int MAX_SUBSTRINGS_PER_FACTOR = 1000;
	private static final int MAX_FACTORS_PER_READ = 100;
	private static final double ERROR_BIN = 0.1;
	private static final int N_ERROR_BINS = 200;
	private static int[] substrings;  // Stores just $repeatID$
	private static int lastSubstring;


	/**
	 * java CheckPeriods ../simulations/forPeriodicSubstrings/genome9/factorization.txt ../simulations/forPeriodicSubstrings/genome9/trueFactorization.txt ../simulations/forPeriodicSubstrings/genome9/repeats.txt 195154 98 ../simulations/forPeriodicSubstrings/genome9/LAshow.txt ../simulations/forPeriodicSubstrings/genome9/readLengths.txt 1> ../simulations/forPeriodicSubstrings/genome9/periodHistogram.txt 2> ../simulations/forPeriodicSubstrings/genome9/periodLog.txt
	 */
	public static void main(String[] args) throws IOException {
		int i, j, p;
		int read, previousRead, length, factorLength, truePeriod;
		double intersection, error;
		int[] histogram;

		loadFactors(args[0]);
		loadTrueFactors(args[1]);
		loadRepeats(args[2]);

		substrings = new int[MAX_SUBSTRINGS_PER_FACTOR];
		Alignments.nAlignments=Integer.parseInt(args[3]);
		Reads.nReads=Integer.parseInt(args[4]);
		Alignments.loadAlignments(args[5]);
		Reads.loadReadLengths(args[6]);
		histogram = new int[N_ERROR_BINS<<1];

		p=0;  // First position of the block of a read in $trueFactors$
		previousRead=-1;
		for (i=0; i<factors.length; i++) {
			read=factors[i][0];
			if (read!=previousRead) {
				while (p<trueFactors.length && trueFactors[p][0]<read) p++;
				previousRead=read;
			}
			factorLength=factors[i][2]-factors[i][1]+1;

			// Building the set of substrings that corresponds to the factor. Only true
			// factors with a large enough intersection with the factor create a
			// substring.
			lastSubstring=-1;
			for (j=p; j<trueFactors.length; j++) {
				if (trueFactors[j][0]!=read || trueFactors[j][1]>factors[i][2]) break;
				if (trueFactors[j][2]<factors[i][1]) continue;
				intersection=Intervals.intersectionLength(factors[i][1],factors[i][2],trueFactors[j][1],trueFactors[j][2]);
				length=trueFactors[j][2]-trueFactors[j][1]+1;
				if (intersection/factorLength<intersectionThreshold && intersection/length<identityThreshold) continue;
				substrings[++lastSubstring]=trueFactors[j][4];
			}

			// Analyzing the period
			truePeriod=getPeriod(factorLength);
			if (truePeriod==-1) {
				if (factors[i][5]!=-1) System.err.println("Error: factor "+i+" is not periodic, but it is reported as periodic.");
				continue;
			}
			else if (factors[i][5]==-1) {
				System.err.println("Error: factor "+i+" is periodic, but it is reported as nonperiodic.");
				continue;
			}
			if (factors[i][5]==0) {
				System.err.println("Warning: factor "+i+" is periodic, but its period was not estimated.");
				continue;
			}
			error=((double)(truePeriod-factors[i][5]))/factors[i][5];
			histogram[(int)(N_ERROR_BINS+error/ERROR_BIN)]++;
		}
		printHistogram(histogram);
	}


	/**
	 * @return the length of a repeat that forms a longest run in $substrings$, and which
	 * has largest surface among all of the repeats with longest run. -1 if
	 * the length of such longest run is smaller than two, or if the surface of such
	 * longest run is too small compared to $factorLength$.
	 */
	private static final int getPeriod(int factorLength) {
		int i;
		int lastRepeat, maxRepeat, frequency, maxFrequency, surface, maxSurface;

		lastRepeat=-1; maxRepeat=-1;
		frequency=-1; maxFrequency=-1;
		surface=-1; maxSurface=-1;
		for (i=0; i<=lastSubstring; i++) {
			if (substrings[i]!=lastRepeat) {
				if (frequency>maxFrequency || (frequency==maxFrequency&&surface>maxSurface)) {
					maxFrequency=frequency;
					maxRepeat=lastRepeat;
					maxSurface=surface;
				}
				lastRepeat=substrings[i];
				frequency=1;
				surface=repeats[substrings[i]][0];
			}
			else {
				frequency++;
				surface+=repeats[substrings[i]][0];
			}
		}
		if (frequency>maxFrequency || (frequency==maxFrequency&&surface>maxSurface)) {
			maxFrequency=frequency;
			maxRepeat=lastRepeat;
			maxSurface=surface;
		}
		if (maxFrequency<2 || ((double)maxSurface)/factorLength<identityThreshold) return -1;
		return repeats[maxRepeat][0];
	}


	private static final void printHistogram(int[] histogram) {
		for (int i=0; i<histogram.length; i++) System.out.println((i>=N_ERROR_BINS?(i-N_ERROR_BINS+1)*ERROR_BIN:(i-N_ERROR_BINS)*ERROR_BIN)+","+histogram[i]);
	}


	private static final void loadFactors(String path) throws IOException {
		int i, j;
		String str;
		String[] tokens;
		BufferedReader br;

		// Determining the number of factors
		br = new BufferedReader(new FileReader(path));
		i=0;
		str=br.readLine();
		while (str!=null && str.length()>0) {
			i++;
			str=br.readLine();
		}
		br.close();
		factors = new int[i][9];

		// Loading factors
		br = new BufferedReader(new FileReader(path));
		i=0;
		str=br.readLine();
		while (str!=null && str.length()>0) {
			tokens=str.split(",");
			for (j=0; j<tokens.length; j++) factors[i][j]=Integer.parseInt(tokens[j]);
			i++;
			str=br.readLine();
		}
		br.close();
	}


	private static final void loadTrueFactors(String path) throws IOException {
		int i, j;
		String str;
		String[] tokens;
		BufferedReader br;

		// Determining the number of factors
		br = new BufferedReader(new FileReader(path));
		i=0;
		str=br.readLine();
		while (str!=null && str.length()>0) {
			i++;
			str=br.readLine();
		}
		br.close();
		trueFactors = new int[i][8];

		// Loading factors
		br = new BufferedReader(new FileReader(path));
		i=0;
		str=br.readLine();
		while (str!=null && str.length()>0) {
			tokens=str.split(",");
			for (j=0; j<tokens.length; j++) trueFactors[i][j]=Integer.parseInt(tokens[j]);
			i++;
			str=br.readLine();
		}
		br.close();
	}


	private static final void loadRepeats(String path) throws IOException {
		int i, j;
		String str;
		String[] tokens;
		BufferedReader br;

		// Determining the number of repeats
		br = new BufferedReader(new FileReader(path));
		i=0;
		str=br.readLine();
		while (str!=null && str.length()>0) {
			i++;
			str=br.readLine();
		}
		br.close();
		repeats = new int[i][2];

		// Loading factors
		br = new BufferedReader(new FileReader(path));
		i=0;
		str=br.readLine();
		while (str!=null && str.length()>0) {
			tokens=str.split(",");
			for (j=0; j<tokens.length; j++) repeats[i][j]=Integer.parseInt(tokens[j]);
			i++;
			str=br.readLine();
		}
		br.close();
	}

}