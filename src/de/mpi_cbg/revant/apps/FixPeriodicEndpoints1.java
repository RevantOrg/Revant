package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.factorize.Reads;

/**
 * Collects short non-periodic spacers between periodic endpoints, and assigns arbitrary
 * breakpoints to them.
 *
 * Remark: this program processes all reads and all alignments sequentially. This is done
 * just for simplicity. Collecting spacers could be done in parallel on different chunks 
 * of reads, but assigning breakpoints needs all spacers in memory. 
 */
public class FixPeriodicEndpoints1 {
	
	public static void main(String[] args) throws IOException {
		final int MAX_SPACER_LENGTH = Integer.parseInt(args[0]);
		final int MIN_ALIGNMENT_LENGTH_READ_REPEAT = Integer.parseInt(args[1]);
		final int N_READS = Integer.parseInt(args[2]);
		final String READ_IDS_FILE = args[3];
		final String READ_LENGTHS_FILE = args[4];
		final String ALPHABET_FILE = args[5];
		final String TRANSLATED_READS_CHARACTERS_FILE = args[6];
		final String TRANSLATED_READS_BOUNDARIES_FILE = args[7];
		final String FULLY_UNIQUE_FILE = args[8];
		final int N_FULLY_UNIQUE = Integer.parseInt(args[9]);
		final String FULLY_CONTAINED_FILE = args[10];
		final int N_FULLY_CONTAINED = Integer.parseInt(args[11]);
		final String READ_READ_ALIGNMENTS_FILE = args[12];
		final int N_BLOCKS = Integer.parseInt(args[13]);  // Of the alignments file
		final String LAST_READ_FILE = args[14];  // Last block included
		final String OUTPUT_PREFIX = args[15];
		
		int i, j;
		int maxBlockLength;
		String str1, str2;
		BufferedReader br1, br2;
		BufferedWriter bw1, bw2;
		int[] lastRead;
		int[] tmpArray = new int[2];
		
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		maxBlockLength=RepeatAlphabet.loadAllBoundaries(TRANSLATED_READS_CHARACTERS_FILE,true,true,TRANSLATED_READS_BOUNDARIES_FILE);
		RepeatAlphabet.loadReadsFully(FULLY_UNIQUE_FILE,N_FULLY_UNIQUE,FULLY_CONTAINED_FILE,N_FULLY_CONTAINED);
		RepeatAlphabet.loadSpacers(MAX_SPACER_LENGTH,maxBlockLength);
		RepeatAlphabet.loadSpacerNeighbors(READ_READ_ALIGNMENTS_FILE,MIN_ALIGNMENT_LENGTH_READ_REPEAT,tmpArray);		
		RepeatAlphabet.assignBreakpoints();
		
RepeatAlphabet.printSpacerNeighbors("/Users/ramseysnow/Downloads/SIMULATED-REPBASE/spacers.dot");		
		
		lastRead = new int[N_BLOCKS];
		br1 = new BufferedReader(new FileReader(LAST_READ_FILE));
		for (i=0; i<N_BLOCKS; i++) lastRead[i]=Integer.parseInt(br1.readLine());
		br1.close();
		RepeatAlphabet.serializeSpacers(OUTPUT_PREFIX,lastRead);
		j=0;
		bw1 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"ids-0.txt"));
		bw2 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"lengths-0.txt"));
		br1 = new BufferedReader(new FileReader(READ_IDS_FILE));
		br2 = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null) {
			if (Integer.parseInt(str1)>lastRead[j]) {
				bw1.close(); bw2.close();
				j++;
				bw1 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"ids-"+j+".txt"));
				bw2 = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"lengths-"+j+".txt"));
			}
			bw1.write(str1); bw1.newLine();
			bw2.write(str2); bw2.newLine();
			str1=br1.readLine(); str2=br2.readLine();
		}
		bw1.close(); bw2.close(); br1.close(); br2.close();
	}

}