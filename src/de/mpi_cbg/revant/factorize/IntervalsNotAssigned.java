package de.mpi_cbg.revant.factorize;

import java.io.*;


public class IntervalsNotAssigned {

	/**
	 * java IntervalsNotAssigned ../simulations/drosophila/LAshow-all.txt ../simulations/drosophila/connection-all.txt ../simulations/drosophila/intervals-periodic-all.txt ../simulations/drosophila/intervals-dense-all.txt ../simulations/drosophila/intervals-alignments-all.txt
	 *
	 * See $Factorize.intervalsNotAssigned$ for details.
	 */
	public static void main(String[] args) throws IOException {
		final String aligmentsFilePath = args[0];
		final String connectionFilePath = args[1];
		final String periodicSubstringsPath = args[2];
		final String denseSubstringsPath = args[3];
		final String alignmentsPath = args[4];
		int[] out = new int[3];
		Factorize.intervalsNotAssigned(aligmentsFilePath,connectionFilePath,periodicSubstringsPath,denseSubstringsPath,alignmentsPath,out);
		System.out.println("Intervals not assigned: "+out[0]+" periodic, "+out[1]+" dense, "+out[2]+" alignment.");
	}

}