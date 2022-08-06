package de.mpi_cbg.revant.apps;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import java.io.*;

/**
 * Expands the alphabet by wobbling existing periodic characters, as well as all the non-
 * periodic characters flagged in some input bitvectors. Wobbling means creating a
 * character that is similar to an existing one but has slightly different length.
 */
public class WobbleCreateAlphabet2 {

	public static void main(String[] args) throws IOException {
		final String ALPHABET_FILE = args[0];
		final int WOBBLE_LENGTH = Integer.parseInt(args[1]);
		final String FLAGS_FILE_PREFIX = args[2];
		final int LAST_FLAG_FILE = Integer.parseInt(args[3]);
		final String OUTPUT_FILE_ALPHABET = args[4];
		final String OUTPUT_FILE_OLD2NEW = args[5];
		
		int i, j;
		int nFlags, lastUnique_old, lastPeriodic_old, lastAlphabet_old, lastUnique_new, lastPeriodic_new, lastAlphabet_new;
		String str;
		BufferedReader br;
		boolean[] flags;
		int[] out = new int[3];
		RepeatAlphabet.Character[] alphabet_old, alphabet_new;
		
		System.err.println("Merging bitvectors...");
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		flags = new boolean[RepeatAlphabet.lastAlphabet+1];
		Math.set(flags,RepeatAlphabet.lastAlphabet,false);
		for (i=0; i<=LAST_FLAG_FILE; i++) {
			br = new BufferedReader(new FileReader(FLAGS_FILE_PREFIX+"-"+i+".txt"));
			for (j=0; j<=RepeatAlphabet.lastAlphabet; j++) flags[j]|=Integer.parseInt(br.readLine())==1;
			br.close();
		}
		nFlags=0;
		for (i=0; i<=RepeatAlphabet.lastAlphabet; i++) nFlags+=flags[i]?1:0;
		
		System.err.println("Wobbling "+nFlags+" characters...");
		alphabet_old=RepeatAlphabet.alphabet; lastUnique_old=RepeatAlphabet.lastUnique; lastPeriodic_old=RepeatAlphabet.lastPeriodic; lastAlphabet_old=RepeatAlphabet.lastAlphabet;
		alphabet_new=RepeatAlphabet.wobble_extendAlphabet(flags,nFlags,WOBBLE_LENGTH,IO.quantum,out);
		lastUnique_new=out[0]; lastPeriodic_new=out[1]; lastAlphabet_new=out[2];
		RepeatAlphabet.alphabet=alphabet_new;
		RepeatAlphabet.lastUnique=lastUnique_new; RepeatAlphabet.lastPeriodic=lastPeriodic_new; RepeatAlphabet.lastAlphabet=lastAlphabet_new;
		RepeatAlphabet.compactInstances_weak();
		RepeatAlphabet.closeAlphabetByRC();
		System.err.println("Wobbled alphabet created: "+(lastUnique_new-lastUnique_old)+" new unique characters; "+(lastPeriodic_new-lastUnique_new-lastPeriodic_old+lastUnique_old)+" new periodic characters; "+(lastAlphabet_new-lastPeriodic_new-lastAlphabet_old+lastPeriodic_old)+" new non-periodic characters.");
		RepeatAlphabet.serializeAlphabet(OUTPUT_FILE_ALPHABET);
		RepeatAlphabet.wobble_buildOld2New(alphabet_old,lastAlphabet_old,alphabet_new,lastAlphabet_new,OUTPUT_FILE_OLD2NEW);
	}

}