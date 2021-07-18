package de.mpi_cbg.revant.consensus;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.IO;


/**
 * Assume that we are given a set of blocks $B_x..B_y$, and a set of repeat tracks on 
 * those blocks to be used to soft-mask alignments (let's call such repeat tracks "input 
 * repeat tracks"). Assume also that we computed all (soft-masked) pairwise alignments 
 * inside $B_x..B_y$, and that we ran our repeat inference up to the consensi. Finally, 
 * assume that we mapped the consensi back to $B_x..B_y$: let's call "output repeat track"
 * the merge of all alignments of such consensi to the reads in $B_x..B_y$. The program 
 * compares the output repeat track to the merge of all intervals in $B_x..B_y$ 
 * produced by $Factorize.java$, and to the merge of all intervals in $B_x..B_y$ that 
 * belong to a basin descriptor file.
 */
public class PippelCheckRepeatTrackOfDatabase {

	/**
	 * @param args
	 * 2: input repeat track ("null" discards it);
	 * 3: input format (see $FilterAlignments.loadRepeatTrack()$);
	 * 4: output repeat track;
	 * 5: output format (see $FilterAlignments.loadRepeatTrack()$).
	 */
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		final int N_READS = Integer.parseInt(args[1]);
		final String REPEAT_TRACK_INPUT = args[2];
		final byte REPEAT_TRACK_INPUT_FORMAT = Byte.parseByte(args[3]);
		final String REPEAT_TRACK_OUTPUT = args[4];
		final byte REPEAT_TRACK_OUTPUT_FORMAT = Byte.parseByte(args[5]);
		final String QUALITY_THRESHOLDS_FILE = args[6];
		final String QUALITIES_FILE = args[7].equalsIgnoreCase("null")?null:args[7];
		
		final String FACTORIZE_OUTPUT_DIR = STEP1_DIR+"/..";
		final String READ_LENGTHS_FILE = FACTORIZE_OUTPUT_DIR+"/"+IO.READS_LENGTHS;
		final String READ_IDS_FILE = FACTORIZE_OUTPUT_DIR+"/"+IO.READS_IDS;
		ConsensusStep1.checkConsensiRepeatTrack(STEP1_DIR,REPEAT_TRACK_INPUT,REPEAT_TRACK_INPUT_FORMAT,REPEAT_TRACK_OUTPUT,REPEAT_TRACK_OUTPUT_FORMAT,N_READS,READ_LENGTHS_FILE,READ_IDS_FILE,QUALITY_THRESHOLDS_FILE,QUALITIES_FILE);
	}


}