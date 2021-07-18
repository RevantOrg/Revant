package de.mpi_cbg.revant.util;

import java.util.Random;
import java.io.*;

/**
 * Takes a random subset of the reads in a FASTA file.
 */
public class RandomSubset {

	public static void main(String[] args) throws IOException {
		final int MIN_BASES = Integer.parseInt(args[0]);
		final int N_READS = Integer.parseInt(args[1]);
		final String LENGTHS_FILE = args[2];
		
		int i, sum;
		Random random = new Random();
		BufferedReader br;
		int[] lengths = new int[N_READS];
		boolean[] selected = new boolean[N_READS];
		
		br = new BufferedReader(new FileReader(LENGTHS_FILE));
		for (i=0; i<N_READS; i++) lengths[i]=Integer.parseInt(br.readLine());
		br.close();
		for (i=0; i<N_READS; i++) selected[i]=false;
		sum=0;
		while (sum<MIN_BASES) {
			i=random.nextInt(N_READS);
			if (selected[i]) continue;
			selected[i]=true; sum+=lengths[i];
			System.out.println(i+"");
		}
	}

}