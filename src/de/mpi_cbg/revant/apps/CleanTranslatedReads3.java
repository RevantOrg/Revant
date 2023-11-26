package de.mpi_cbg.revant.apps;

import java.io.*;

import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.Reads;

/**
 * Updates the translation of every read to use the alphabet built by
 * $CleanTranslatedReads2$.
 */
public class CleanTranslatedReads3 {
	
	public static void main(String[] args) throws IOException {
		final int N_READS = Integer.parseInt(args[0]);
		final String READ_IDS_FILE = args[1];
		final String READ_LENGTHS_FILE = args[2];
		final String ALPHABET_FILE_OLD = args[3];
		final String ALPHABET_COUNTS_FILE_OLD = args[4];
		final String TRANSLATED_FILE_CHARACTERS_OLD = args[5];
		final String TRANSLATED_FILE_BOUNDARIES_OLD = args[6];
        final int AVG_READ_LENGTH = Integer.parseInt(args[7]);
        final int SPANNING_BPS = Integer.parseInt(args[8]);
        final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[9]);  // Read-repeat
        final long GENOME_LENGTH = Long.parseLong(args[10]);  // Of one haplotype
        final int N_HAPLOTYPES = Integer.parseInt(args[11]);
		final boolean KEEP_PERIODIC = Integer.parseInt(args[12])==1;
		final String ALPHABET_FILE_NEW = args[13];
		final String OLD2NEW_FILE = args[14];
		final String TRANSLATED_FILE_CHARACTERS_NEW = args[15];
		final String TRANSLATED_FILE_BOUNDARIES_NEW = args[16];
		final String HISTOGRAM_FILE = args[17];
		final String FULLY_UNIQUE_FILE_NEW = args[18];
		final int LAST_TRANSLATED_READ = Integer.parseInt(args[19]);  // Same as in $TranslateReads.java$. -1 if no read has been updated yet.
		
        final double SIGNIFICANCE_LEVEL = 0.05;  // Conventional
        final int MIN_MISSING_LENGTH = IO.quantum;  // Arbitrary
        final int MAX_HISTOGRAM_LENGTH = 1000;  // Arbitrary
		
		int i;
		int lastUnique_old, lastPeriodic_old, lastAlphabet_old;
		int lastUnique_new, lastPeriodic_new, lastAlphabet_new;
		int nBlocks;
		String str1, str2;
		RepeatAlphabet.Character tmpChar = new RepeatAlphabet.Character();
		BufferedReader br1, br2;
		BufferedWriter bw1, bw2, bw3;
		int[] old2new;
		long[] blocksHistogram;
		RepeatAlphabet.Character[] oldAlphabet, newAlphabet;
		
		// Loading old alphabet
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_OLD,2);
		RepeatAlphabet.loadAlphabetCount(ALPHABET_COUNTS_FILE_OLD,RepeatAlphabet.lastAlphabet+1);
		oldAlphabet=RepeatAlphabet.alphabet;
		lastUnique_old=RepeatAlphabet.lastUnique;
		lastPeriodic_old=RepeatAlphabet.lastPeriodic;
		lastAlphabet_old=RepeatAlphabet.lastAlphabet;
		br1 = new BufferedReader(new FileReader(OLD2NEW_FILE));
		old2new = new int[lastAlphabet_old-lastUnique_old];
		for (i=0; i<old2new.length; i++) old2new[i]=Integer.parseInt(br1.readLine());
		br1.close();
		
		// Updating the translation
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE_NEW,2);
		newAlphabet=RepeatAlphabet.alphabet;
		lastUnique_new=RepeatAlphabet.lastUnique;
		lastPeriodic_new=RepeatAlphabet.lastPeriodic;
		lastAlphabet_new=RepeatAlphabet.lastAlphabet;
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,Reads.nReads);
		blocksHistogram = new long[MAX_HISTOGRAM_LENGTH+1];
		Math.set(blocksHistogram,MAX_HISTOGRAM_LENGTH,0L);
		br1 = new BufferedReader(new FileReader(TRANSLATED_FILE_CHARACTERS_OLD));
		br2 = new BufferedReader(new FileReader(TRANSLATED_FILE_BOUNDARIES_OLD));
		bw1 = new BufferedWriter(new FileWriter(TRANSLATED_FILE_CHARACTERS_NEW));
		bw2 = new BufferedWriter(new FileWriter(TRANSLATED_FILE_BOUNDARIES_NEW));
		bw3 = new BufferedWriter(new FileWriter(FULLY_UNIQUE_FILE_NEW));
		i=LAST_TRANSLATED_READ+1; str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null && i<N_READS) {            
            nBlocks=RepeatAlphabet.cleanTranslatedRead_updateTranslation(str1,str2,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,old2new,Reads.getReadLength(i),KEEP_PERIODIC,IO.quantum,N_READS,AVG_READ_LENGTH,SPANNING_BPS,GENOME_LENGTH,N_HAPLOTYPES,MIN_ALIGNMENT_LENGTH,MIN_MISSING_LENGTH,SIGNIFICANCE_LEVEL,bw1,bw2,tmpChar);
			blocksHistogram[nBlocks]++;
			bw1.newLine(); bw2.newLine();
			if (nBlocks==0) bw3.write(Reads.readIDs[i]+"\n");
			i++; str1=br1.readLine(); str2=br2.readLine();
		}
		br1.close(); br2.close(); bw1.close(); bw2.close(); bw3.close();
		bw1 = new BufferedWriter(new FileWriter(HISTOGRAM_FILE));
		for (i=0; i<=MAX_HISTOGRAM_LENGTH; i++) bw1.write(blocksHistogram[i]+"\n");
		bw1.close();
	}

}