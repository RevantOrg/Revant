package de.mpi_cbg.revant.apps;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import java.io.*;

/**
 * Expands the alphabet by wobbling existing periodic characters, as well as characters 
 * that are adjacent to periodic characters in some translation. Wobbling means creating a
 * character that is similar to an existing one but has slightly different length.
 */
public class WobbleCreateAlphabet {

	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final String TRANSLATED_READS_CHARACTERS_FILE = args[1];
		final int WOBBLE_LENGTH = Integer.parseInt(args[2]);
		final String OUTPUT_FILE_ALPHABET = args[3];
		final String OUTPUT_FILE_OLD2NEW = args[4];
		
		int i;
		int nBlocks, nFlags, lastUnique_old, lastPeriodic_old, lastAlphabet_old, lastUnique_new, lastPeriodic_new, lastAlphabet_new;
		String str;
		BufferedReader br;
		boolean[] flags;
		int[] tmpArray;
		int[] out = new int[3];
		RepeatAlphabet.Character[] alphabet_old, alphabet_new;
		
		System.err.println("Finding non-periodic characters to be wobbled... ");
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		alphabet_old=RepeatAlphabet.alphabet; lastUnique_old=RepeatAlphabet.lastUnique; lastPeriodic_old=RepeatAlphabet.lastPeriodic; lastAlphabet_old=RepeatAlphabet.lastAlphabet;
		flags = new boolean[RepeatAlphabet.lastAlphabet+1];
		Math.set(flags,RepeatAlphabet.lastAlphabet,false);
		tmpArray = new int[100];  // Arbitrary
		br = new BufferedReader(new FileReader(TRANSLATED_READS_CHARACTERS_FILE));
		str=br.readLine();
		while (str!=null) {
			nBlocks=1+((str.length()+1)>>1);  // Loose upper bound
			if (tmpArray.length<nBlocks) tmpArray = new int[nBlocks];
			RepeatAlphabet.wobble_markAlphabet(str,flags,tmpArray);
			str=br.readLine();
		}
		br.close();
		nFlags=0;
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) nFlags+=flags[i]?1:0;
		
		System.err.println("Wobbling "+nFlags+" characters...");
		alphabet_new=RepeatAlphabet.wobble_extendAlphabet(flags,nFlags,WOBBLE_LENGTH,IO.quantum,out);
		lastUnique_new=out[0]; lastPeriodic_new=out[1]; lastAlphabet_new=out[2];
		RepeatAlphabet.alphabet=alphabet_new;
		RepeatAlphabet.lastUnique=lastUnique_new; RepeatAlphabet.lastPeriodic=lastPeriodic_new; RepeatAlphabet.lastAlphabet=lastAlphabet_new;
		RepeatAlphabet.compactInstances();
		RepeatAlphabet.closeAlphabetByRC();
		// $compactInstances()$ might be needed again: see $closeAlphabetByRC()$.
		RepeatAlphabet.compactInstances();
		System.err.println("Wobbling completed: "+(lastUnique_new-lastUnique_old)+" new unique characters; "+(lastPeriodic_new-lastUnique_new-lastPeriodic_old+lastUnique_old)+" new periodic characters; "+(lastAlphabet_new-lastPeriodic_new-lastAlphabet_old+lastPeriodic_old)+" new non-periodic characters.");
		RepeatAlphabet.serializeAlphabet(OUTPUT_FILE_ALPHABET);
		RepeatAlphabet.wobble_buildOld2New(alphabet_old,lastAlphabet_old,alphabet_new,lastAlphabet_new,OUTPUT_FILE_OLD2NEW);
	}

}