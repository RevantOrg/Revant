package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Reads;

/**
 * Collects new unique characters induced by removing rare characters.
 * Prints to STDOUT the new value of $maxOpenLength_unique$.
 */
public class CleanTranslatedReads1 {
	/**
	 * @param args 8: keep periodic characters untouched (1/0).
	 */
	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String ALPHABET_COUNTS_FILE = args[1];
		final int N_READS = Integer.parseInt(args[2]);
		final String READ_IDS_FILE = args[3];
		final String READ_LENGTHS_FILE = args[4];
        final int AVG_READ_LENGTH = Integer.parseInt(args[5]);
        final int SPANNING_BPS = Integer.parseInt(args[6]);
        final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[7]);  // Read-repeat
        final long GENOME_LENGTH = Long.parseLong(args[8]);  // Of one haplotype
        final int N_HAPLOTYPES = Integer.parseInt(args[9]);
		final String TRANSLATED_READS_CHARACTERS_FILE = args[10];
		final String TRANSLATED_READS_BOUNDARIES_FILE = args[11];
		final boolean KEEP_PERIODIC = Integer.parseInt(args[12])==1;
		final String OUTPUT_FILE = args[13];
		
        final double SIGNIFICANCE_LEVEL = 0.05;  // Conventional
        final int MIN_MISSING_LENGTH = IO.quantum;  // Arbitrary
        
		int i;
		String str1, str2;
		RepeatAlphabet.Character tmpChar = new RepeatAlphabet.Character();
		BufferedReader br1, br2;
		BufferedWriter bw;
		
		Reads.nReads=N_READS;
		Reads.loadReadIDs(READ_IDS_FILE,Reads.nReads);
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		RepeatAlphabet.loadAlphabetCount(ALPHABET_COUNTS_FILE,RepeatAlphabet.lastAlphabet+1);		
		br1 = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
		br2 = new BufferedReader(new FileReader(TRANSLATED_READS_BOUNDARIES_FILE));
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		i=0; str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null) {
			RepeatAlphabet.cleanTranslatedRead_collectCharacterInstances(str1,str2,Reads.readLengths[i],KEEP_PERIODIC,IO.quantum,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,MIN_ALIGNMENT_LENGTH,MIN_MISSING_LENGTH,SIGNIFICANCE_LEVEL,bw,tmpChar);
			i++; str1=br1.readLine(); str2=br2.readLine();
		}
		br1.close(); br2.close(); bw.close();
		System.out.println(RepeatAlphabet.maxOpenLength_unique+"");
	}

}