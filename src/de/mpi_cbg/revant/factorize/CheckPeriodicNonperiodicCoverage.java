package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;


/**
 * Plots the distribution of differences between full and periodic coverage for all
 * periodic estimated factors.
 */
public class CheckPeriodicNonperiodicCoverage {

	private static int[][] factors;
	private static final double DIFFERENCE_BIN = 0.1;


	/**
	 * java CheckPeriodicNonperiodicCoverage ../simulations/forPeriodicSubstrings/genome9/factorization.txt > ../simulations/forPeriodicSubstrings/genome9/periodicNonperiodicCoverage.txt
	 */
	public static void main(String[] args) throws IOException {
		int i, bin;
		double difference;
		int[] histogram2 = new int[250];

		loadFactors(args[0]);
		for (i=0; i<factors.length; i++) {
			if (factors[i][5]==-1) continue;
			difference=factors[i][6]-factors[i][7];
			difference/=factors[i][7];

if (difference>5) System.err.println("factor "+i+", relative difference: "+difference);

			bin=(int)(difference/DIFFERENCE_BIN);
			if (bin>=histogram2.length) bin=histogram2.length-1;
			histogram2[bin]++;
		}
		printHistogram(histogram2);
	}


	private static final void printHistogram(int[] histogram) {
		for (int i=0; i<histogram.length; i++) System.out.println(i*DIFFERENCE_BIN+","+histogram[i]);
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

}