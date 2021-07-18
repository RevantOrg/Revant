package de.mpi_cbg.revant.consensus;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.AlignmentInterval;
import de.mpi_cbg.revant.factorize.DenseSubstring;
import de.mpi_cbg.revant.factorize.PeriodicSubstringInterval;
import de.mpi_cbg.revant.factorize.Intervals;
import de.mpi_cbg.revant.factorize.FilterAlignments;
import de.mpi_cbg.revant.intervalgraph.IntervalGraph;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;
import de.mpi_cbg.revant.fragments.FragmentsStep1;
import de.mpi_cbg.revant.fragments.FragmentsStep3;


/**
 * Similar to $FragmentsStep3$, but for the case in which the reference is the output of a
 * consensus algorithm.
 *
 * Remark: the program builds the histogram of error rates for all fragment-consensus 
 * alignments. For short-period repeats such histogram might be distorted by the aligner:
 * see $compareConsensusDominatorAlignments()$ for details.
 *
 * Remark: "structural" errors in a reference (e.g. the reference is the concatenation of 
 * several repeats) do not seem to affect the distribution of fragment-consensus alignment
 * quality.
 */
public class ConsensusStep1 {
	
	/**
	 * Temporary space
	 */
	private static int lastPairAllIntervals, lastPairAllBasins;
	private static Pair[] pairsAllIntervals, pairsAllBasins, pairs1, pairs2, pairs3;
	private static Pair[][] pairsAllAlignments;
	private static int[] lastPairAllAlignments;
	private static int[] tmpArray, tmpArray2;
	private static int[][] track;
	private static int[] lastTrack;
	private static long[] tmpLong;
	private static int[] frags2cons, frags2dom;
	private static double[][] minError;
	private static long[][] surface;
	private static FragConsAlignment[] consensusWindows, dominatorWindows;
	
	
	/**
	 * @param args 
	 *  0: Step1 directory;
	 *  1: use only fragment-consensus alignments with error rate <= this;
	 *  2: min length of an alignment to use for computing the coverage plot of a new 
	 *     reference (for debugging);
	 *  3: suffix of consensus strings (excluding ".fasta"; e.g. ".daccord.consensus");
	 *  4: total number of reads in the block used for factorization;
	 *  5: min alignment length used in repeat inference;
	 *  6: (0/1) analyze how consensi map to the blocks used for repeat inference;
	 *  7: LAshow file used to compute the consensus-consensus alignment graph;
	 *  8: format of fragment-dominator alignments (0/1, see 
	 *     $compareConsensusDominatorAlignments()$);
	 *  9: format of fragment-consensus alignments;
	 * 10: repeat track built by mapping all consensi to the DB blocks they were inferred
	 *     from, and by merging the intervals (use "null" to discard it);
	 * 11: format of the consensi track above (see $FilterAlignments.loadRepeatTrack()$).
	 */
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		final double FRAGMENT_CONSENSUS_ERROR = Double.parseDouble(args[1]);
		final int MIN_ALIGNMENT_LENGTH_FOR_COVERAGE = Integer.parseInt(args[2]);
		final String CONSENSUS_SUFFIX = args[3];
		final int N_READS = Integer.parseInt(args[4]);
		final int MIN_ALIGNMENT_LENGTH_IN_INFERENCE = Integer.parseInt(args[5]);
		final boolean BLOCK_CONSENSUS_ANALYSIS = Integer.parseInt(args[6])==1;
		final String CONSENSUS_CONSENSUS_ALIGNMENTS = args[7];
		final boolean FD_ALIGNMENTS_FORMAT = Integer.parseInt(args[8])==1;
		final boolean FC_ALIGNMENTS_FORMAT = Integer.parseInt(args[9])==1;
		final String CONSENSI_REPEAT_TRACK = args[10].equalsIgnoreCase("null")?null:args[10];
		final byte CONSENSI_REPEAT_TRACK_FORMAT = Byte.parseByte(args[11]);
		
		final String STEP4_DIR = STEP1_DIR+"/finalOutput/step4";
		final String STEP5_DIR = STEP4_DIR+"/step5";
		final String FRAGMENTS_COORDINATES_DIR = STEP5_DIR+"/fragments";
		final String FRAGMENTS_ALIGNMENTS_DIR = STEP5_DIR+"/fragments-strings-alignments";
		final String DOMINATOR_ALIGNMENTS_DIR = STEP5_DIR+"/fragments-strings-alignments/fragments-strings-alignments-new";
		final String CONSENSI_DIR = FRAGMENTS_ALIGNMENTS_DIR+"/consensi";
		final String STATS_DIR = CONSENSI_DIR+"/stats";
		final String REFERENCE_LENGTHS_DIR = CONSENSI_DIR+"/consensi-strings";
		final String OLD_REFERENCE_LENGTHS_DIR = FRAGMENTS_ALIGNMENTS_DIR+"/references-strings-new";
		final String FRAGMENTS_LENGTHS_DIR = FRAGMENTS_ALIGNMENTS_DIR+"/fragments-strings-new";
		final String ERROR_RATE_HISTOGRAM_FILE = STATS_DIR+"/errorRateHistogram.txt";
		final String ERROR_RATE_HISTOGRAM_SHORTPERIOD_FILE = STATS_DIR+"/errorRateHistogram-shortPeriod.txt";
		final String ERROR_RATE_HISTOGRAM_LONGPERIOD_FILE = STATS_DIR+"/errorRateHistogram-longPeriod.txt";
		final String ERROR_RATE_HISTOGRAM_UNKNOWNPERIOD_FILE = STATS_DIR+"/errorRateHistogram-unknownPeriod.txt";
		final String ERROR_RATE_HISTOGRAM_CYCLIC_FILE = STATS_DIR+"/errorRateHistogram-cyclic.txt";
		final String ERROR_RATE_MATRIX_FILE = STATS_DIR+"/errorRateMatrix.txt";
		final String LENGTH_HISTOGRAM_FILE = STATS_DIR+"/lengthHistogram.txt";
		final String LENGTH_HISTOGRAM_LONGPERIOD_FILE = STATS_DIR+"/lengthHistogram-longPeriod.txt";		
		final String LENGTH_HISTOGRAM_UNKNOWNPERIOD_FILE = STATS_DIR+"/lengthHistogram-unknownPeriod.txt";
		final String FRAGMENT_REFERENCE_HISTOGRAM_FILE = STATS_DIR+"/fragmentReferenceHistogram.txt";
		final String BLOCK_REFERENCE_HISTOGRAM_FILE = STATS_DIR+"/blockReferenceHistogram.txt";
		final String BLOCK_REFERENCE_HISTOGRAM_IDS_FILE = STATS_DIR+"/blockReferenceHistogram-ids.txt";
		final String BLOCK_REFERENCE_HISTOGRAM_PERIODIC_FILE = STATS_DIR+"/blockReferenceHistogram-periodic.txt";
		final String BLOCK_REFERENCE_QUALITIES_FILE = STATS_DIR+"/blockReferenceQualities.txt";
		final String UTILITIES_FILE_C = STATS_DIR+"/utilities-consensus.txt";
		final String UTILITIES_FILE_R = STATS_DIR+"/utilities-repeat.txt";
		final String INTERVAL_QUALITIES_FILE = STATS_DIR+"/intervalQualities.txt";
		final String READ_LENGTHS_FILE = STEP1_DIR+"/../"+IO.READS_LENGTHS;
		final String READ_IDS_FILE = STEP1_DIR+"/../"+IO.READS_IDS;
		final String QUALITY_THRESHOLDS_FILE = STEP1_DIR+"/../qualityThresholds.txt";
		final String QUALITIES_FILE = STEP1_DIR+"/../"+IO.READS_PHRED_PREFIX+IO.getPhredSuffix(STEP1_DIR+"/../"+IO.READS_PHRED_PREFIX);
		final String TRACK_FILE = null;
		final String TRACK_INTERSECTION_FILE = null;
		//final String TRACK_FILE = STEP1_DIR+"/../reads-dust.txt";
		//final String TRACK_INTERSECTION_FILE = STATS_DIR+"/trackIntersection-dust.txt";
		//final String TRACK_FILE = STEP1_DIR+"/../reads-tandem.txt";
		//final String TRACK_INTERSECTION_FILE = STATS_DIR+"/trackIntersection-tandem.txt";
		final String REPEAT_TRACK_INPUT_FILE = STEP1_DIR+"/../input/repmod-union.txt";
		final byte REPEAT_TRACK_INPUT_FORMAT = 0;
		final boolean CONSENSUS_DOMINATOR_VERBOSE = true;

		final int CAPACITY = 100;  // Arbitrary
		final int LENGTH_THRESHOLD = IO.quantum<<2;  // Arbitrary
		final double MIN_ERROR = 0.0;
		final double MAX_ERROR = 0.4;
		final int N_HISTOGRAM_BINS = 100;  // Arbitrary
		final double HISTOGRAM_QUANTUM = (MAX_ERROR-MIN_ERROR)/N_HISTOGRAM_BINS;
		final int MIN_PERIOD_LONG = 1000;  // Arbitrary
		final double SURFACE_THRESHOLD = 0.8;  // Arbitrary
		final double[] ERROR_RATES = new double[] {0.3,0.23,0.2,0.17};
		final byte MIN_LOW_QUALITY = 45;  //51;
		final int IDENTITY_THRESHOLD = IO.quantum;
		
		boolean isShortPeriod;
		int i, j, k, p;
		int sum, length, maxReferenceLength, nTooLong, currentRead, lastTmpWindow, nAlignments;
		int nFragments, nReferences, totalNFragments, totalNFragmentsNonShortPeriod, totalNReferences, nDescriptors, nReferencesWithLowCoverage;
		int PATH_ID, BASIN_TYPE, PERIOD;
		double errorRate;
		String str, id, idPrime;
		String BASINS_FILE, FRAGMENTS_LENGTHS_FILE, REFERENCE_LENGTHS_FILE, CONNECTION_FILE, FR_ALIGNMENTS_FILE;
		File file;
		BufferedReader br;
		BufferedWriter bw, lengthHistogram, lengthHistogramLongPeriod, lengthHistogramUnknownPeriod, fragmentReferenceHistogram;
		BufferedWriter blockReferenceHistogram, blockReferenceHistogramIDs, blockReferenceHistogramPeriodic, blockReferenceQualities;
		BufferedWriter intervalQualitiesFile, trackIntersectionFile;
		long[] tmpIO = new long[36];
		long[] out = new long[10];
		long[] consensusDominatorStats = new long[10];
		int[] fragmentStats = new int[7];
		int[] totalFragmentStats = new int[7];
		int[] referenceLengths = new int[CAPACITY];
		int[] errorRateHistogram = new int[N_HISTOGRAM_BINS];
		int[] errorRateHistogramShortPeriod = new int[N_HISTOGRAM_BINS];
		int[] errorRateHistogramLongPeriod = new int[N_HISTOGRAM_BINS];
		int[] errorRateHistogramUnknownPeriod = new int[N_HISTOGRAM_BINS];
		int[] errorRateHistogramCyclic = new int[N_HISTOGRAM_BINS];
		int[] errorRateHistogramInteractive = new int[N_HISTOGRAM_BINS];
		int[][] errorRateMatrix, qualityHistogram;
		int[] erHistogram, readIDsBackup, readLengthsBackup, qualityHistogramIntervals;
		String[] list;
		FragmentsStep1.PermutationWindow[] tmpWindows;
		int[][] referenceCoverage = new int[CAPACITY][0];
		int[][] fragmentReferenceStats = new int[CAPACITY][2];
		int[][] coverageFromBlock = new int[CAPACITY][1000];
		long[][] trackIntersection;
		double[][] consensusDominatorDetails;
		
		if (TRACK_FILE!=null) trackIntersection = new long[ERROR_RATES.length][3];
		else trackIntersection=null;
		IO.initialize();
		pairsAllIntervals = new Pair[1000];
		for (i=0; i<pairsAllIntervals.length; i++) pairsAllIntervals[i] = new Pair();
		pairsAllBasins = new Pair[1000];
		for (i=0; i<pairsAllBasins.length; i++) pairsAllBasins[i] = new Pair();
		pairsAllAlignments = new Pair[ERROR_RATES.length][1000];
		for (i=0; i<pairsAllAlignments.length; i++) {
			for (j=0; j<pairsAllAlignments[i].length; j++) pairsAllAlignments[i][j] = new Pair();
		}
		lastPairAllAlignments = new int[ERROR_RATES.length];
		pairs1 = new Pair[CAPACITY];
		for (i=0; i<pairs1.length; i++) pairs1[i] = new Pair();
		pairs2 = new Pair[CAPACITY];
		for (i=0; i<pairs2.length; i++) pairs2[i] = new Pair();
		pairs3 = new Pair[CAPACITY];
		for (i=0; i<pairs3.length; i++) pairs3[i] = new Pair();
		qualityHistogram = new int[ERROR_RATES.length][52];
		qualityHistogramIntervals = new int[52];
		Reads.nReads=N_READS;
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,Reads.nReads);
		readIDsBackup=Reads.readIDs;
		tmpArray = new int[100];
		tmpLong = new long[100];
		tmpArray2 = new int[Reads.maxReadLength];
		readLengthsBackup=Reads.readLengths;
		Reads.loadQualities(QUALITY_THRESHOLDS_FILE,QUALITIES_FILE);
		totalNFragments=0; totalNFragmentsNonShortPeriod=0; totalNReferences=0; nReferencesWithLowCoverage=0;
		Math.set(totalFragmentStats,totalFragmentStats.length-1,0);
		Math.set(errorRateHistogram,N_HISTOGRAM_BINS-1,0);
		Math.set(errorRateHistogramShortPeriod,N_HISTOGRAM_BINS-1,0);
		Math.set(errorRateHistogramLongPeriod,N_HISTOGRAM_BINS-1,0);
		Math.set(errorRateHistogramUnknownPeriod,N_HISTOGRAM_BINS-1,0);
		Math.set(errorRateHistogramCyclic,N_HISTOGRAM_BINS-1,0);
		tmpWindows = new FragmentsStep1.PermutationWindow[CAPACITY];
		for (i=0; i<CAPACITY; i++) tmpWindows[i] = new FragmentsStep1.PermutationWindow();
		lengthHistogram = new BufferedWriter(new FileWriter(LENGTH_HISTOGRAM_FILE));
		lengthHistogramLongPeriod = new BufferedWriter(new FileWriter(LENGTH_HISTOGRAM_LONGPERIOD_FILE));
		lengthHistogramUnknownPeriod = new BufferedWriter(new FileWriter(LENGTH_HISTOGRAM_UNKNOWNPERIOD_FILE));
		fragmentReferenceHistogram = new BufferedWriter(new FileWriter(FRAGMENT_REFERENCE_HISTOGRAM_FILE));
		blockReferenceHistogram = new BufferedWriter(new FileWriter(BLOCK_REFERENCE_HISTOGRAM_FILE));
		blockReferenceHistogramIDs = new BufferedWriter(new FileWriter(BLOCK_REFERENCE_HISTOGRAM_IDS_FILE));
		blockReferenceHistogramPeriodic = new BufferedWriter(new FileWriter(BLOCK_REFERENCE_HISTOGRAM_PERIODIC_FILE));
		blockReferenceQualities = new BufferedWriter(new FileWriter(BLOCK_REFERENCE_QUALITIES_FILE));
		intervalQualitiesFile = new BufferedWriter(new FileWriter(INTERVAL_QUALITIES_FILE));
		if (TRACK_FILE!=null) trackIntersectionFile = new BufferedWriter(new FileWriter(TRACK_INTERSECTION_FILE));
		else trackIntersectionFile=null;
		frags2cons = new int[CAPACITY];
		frags2dom = new int[CAPACITY];
		minError = new double[2][CAPACITY];
		surface = new long[2][CAPACITY];
		consensusWindows = new FragConsAlignment[CAPACITY];
		for (i=0; i<consensusWindows.length; i++) consensusWindows[i] = new FragConsAlignment();
		dominatorWindows = new FragConsAlignment[CAPACITY];
		for (i=0; i<dominatorWindows.length; i++) dominatorWindows[i] = new FragConsAlignment();
		Math.set(consensusDominatorStats,consensusDominatorStats.length-1,0);
		
		file = new File(CONSENSI_DIR);
		list=file.list();
		System.err.println("Descriptors:");
		nDescriptors=0;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			System.err.println(i+": "+list[i]);
			nDescriptors++;
		}
		


//consensusAlignmentsGraph(true,REFERENCE_LENGTHS_DIR+"/LAshow-consensi-consensi.txt",MIN_ALIGNMENT_LENGTH_FOR_COVERAGE,/* int nConsensi*/316,REFERENCE_LENGTHS_DIR,REFERENCE_LENGTHS_DIR+"/list-consensi.txt",REFERENCE_LENGTHS_DIR+"/connection.txt",REFERENCE_LENGTHS_DIR+"/consensus-consensus.dot");
		
		
		
		errorRateMatrix = new int[nDescriptors][N_HISTOGRAM_BINS];
		Math.set(errorRateMatrix,0);
		consensusDominatorDetails = new double[nDescriptors][6];
		Math.set(consensusDominatorDetails,0.0);
		lastPairAllBasins=loadAllBasinIntervals(list,STEP4_DIR);
		lastPairAllIntervals=loadAllIntervals(STEP1_DIR,qualityHistogramIntervals);
		if (TRACK_FILE!=null) {
			track = new int[N_READS][10]; lastTrack = new int[N_READS];
			loadTrack(TRACK_FILE,MIN_ALIGNMENT_LENGTH_IN_INFERENCE,track,lastTrack);
		}
		for (i=0; i<qualityHistogramIntervals.length; i++) intervalQualitiesFile.write(qualityHistogramIntervals[i]+",");
		intervalQualitiesFile.write("\n");
		intervalQualitiesFile.close();
		p=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			p++;
			
			// Initializing data structures for the current repeat
			id=list[i].substring(IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH);
			System.err.print("Loading data structures for repeat "+id+"... ");
			
			// Basin descriptor
			BASINS_FILE=STEP4_DIR+"/"+IO.BASIN_PREFIX+id+".txt";
			br = new BufferedReader(new FileReader(BASINS_FILE));
			str=br.readLine(); br.close();
			IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
			PATH_ID=tmpArray[2]; BASIN_TYPE=tmpArray[3]; PERIOD=tmpArray[4];
			
			// Loading fragment lengths
			FRAGMENTS_LENGTHS_FILE=FRAGMENTS_LENGTHS_DIR+"/fragments-"+id+"-lengths.txt";
			br = new BufferedReader(new FileReader(FRAGMENTS_LENGTHS_FILE));
			nFragments=0; str=br.readLine();
			while (str!=null) {
				nFragments++;
				str=br.readLine();
			}
			br.close();
			Reads.nReads=nFragments;
			Reads.loadReadLengths(FRAGMENTS_LENGTHS_FILE);
			Reads.loadReadIDs(0,nFragments);
			totalNFragments+=nFragments;
			if (BASIN_TYPE!=Constants.INTERVAL_PERIODIC || PERIOD>=MIN_PERIOD_LONG) totalNFragmentsNonShortPeriod+=nFragments;
			
			// Loading reference lengths
			idPrime=id.replace('-','_');
			REFERENCE_LENGTHS_FILE=REFERENCE_LENGTHS_DIR+"/"+idPrime+CONSENSUS_SUFFIX+"-lengths.txt";
			br = new BufferedReader(new FileReader(REFERENCE_LENGTHS_FILE));
			nReferences=0; maxReferenceLength=0; str=br.readLine();
			while (str!=null) {
				length=Integer.parseInt(str);
				referenceLengths[nReferences++]=length;
				maxReferenceLength=Math.max(maxReferenceLength,length);
				str=br.readLine();
			}
			br.close();
			for (j=0; j<nReferences; j++) {
				if (referenceCoverage[j].length<referenceLengths[j]) referenceCoverage[j] = new int[referenceLengths[j]];
				Math.set(referenceCoverage[j],referenceLengths[j]-1,0);
			}
			totalNReferences+=nReferences;
			System.err.println("done");
			
			// Loading old reference lengths
			REFERENCE_LENGTHS_FILE=OLD_REFERENCE_LENGTHS_DIR+"/reference-"+id+"-lengths.txt";
			br = new BufferedReader(new FileReader(REFERENCE_LENGTHS_FILE));
			for (j=0; j<nReferences; j++) {
				str=br.readLine();
				if (BASIN_TYPE!=Constants.INTERVAL_PERIODIC) lengthHistogram.write(Integer.parseInt(str)+","+referenceLengths[j]+"\n");
				else if (PERIOD>=MIN_PERIOD_LONG) {
					lengthHistogramLongPeriod.write(Integer.parseInt(str)+","+referenceLengths[j]+"\n");
//if ( (((double)Integer.parseInt(str))-referenceLengths[j])/referenceLengths[j]>=0.4 ) System.err.println("VITTU> id="+id+"  anomaly="+((((double)Integer.parseInt(str))-referenceLengths[j])/referenceLengths[j]));
				}
				else if (PERIOD==0) lengthHistogramUnknownPeriod.write(Integer.parseInt(str)+","+referenceLengths[j]+"\n");
			}
			br.close();
			
			// Checking whether fragments are too long
			nTooLong=0;
			if (BASIN_TYPE!=Constants.INTERVAL_PERIODIC && PATH_ID!=-1) {
				for (j=0; j<nFragments; j++) {
					length=Reads.getReadLength(j);
					if (length>maxReferenceLength+LENGTH_THRESHOLD) {
						nTooLong++;
						System.err.println("WARNING: Fragment "+j+" has length "+length+" > maxReferenceLength="+maxReferenceLength);
					}
				}
			}			
			
			// Checking fragment-reference alignments
			Math.set(fragmentStats,fragmentStats.length-1,0);
			Math.set(fragmentReferenceStats,0);
			fragmentStats[6]=nTooLong;
			if (BASIN_TYPE==Constants.INTERVAL_PERIODIC) {
				if (PERIOD==0) erHistogram=errorRateHistogramUnknownPeriod;
				else if (PERIOD < MIN_PERIOD_LONG) erHistogram=errorRateHistogramShortPeriod;
				else erHistogram=errorRateHistogramLongPeriod;
				blockReferenceHistogramPeriodic.write(PERIOD+"\n");
			}
			else {
				if (PATH_ID==-1) {
					erHistogram=errorRateHistogramCyclic;
					blockReferenceHistogramPeriodic.write("-1\n");
				}
				else {
					erHistogram=errorRateHistogram;
					blockReferenceHistogramPeriodic.write("-2\n");
				}
			}			
			FR_ALIGNMENTS_FILE=CONSENSI_DIR+"/"+list[i]+"/LAshow-fragments-consensus.txt";
			br = new BufferedReader(new FileReader(FR_ALIGNMENTS_FILE),IO.BUFFER_SIZE);
			currentRead=-1; lastTmpWindow=-1;
			str=br.readLine(); str=br.readLine();  // Skipping header
			str=br.readLine(); nAlignments=0;
			while (str!=null) {
				nAlignments++;
				Alignments.readAlignmentFile(str);
				errorRate=((double)(Alignments.diffs<<1))/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2);
				if (errorRate>FRAGMENT_CONSENSUS_ERROR) {
					str=br.readLine();
					continue;
				}
				if (Alignments.readA-1!=currentRead) {
					if (currentRead!=-1) {
						FragmentsStep3.alignmentStats(currentRead,tmpWindows,lastTmpWindow,fragmentStats,referenceLengths,referenceCoverage,MIN_ALIGNMENT_LENGTH_FOR_COVERAGE);
						fragmentsPerReference(tmpWindows,lastTmpWindow,fragmentReferenceStats);
						if (BASIN_TYPE!=Constants.INTERVAL_PERIODIC) {
							if (PATH_ID!=-1) incrementErrorRateHistogram_nonPeriodic(currentRead,SURFACE_THRESHOLD,tmpWindows,lastTmpWindow,errorRateMatrix[p],MIN_ERROR,MAX_ERROR,N_HISTOGRAM_BINS,HISTOGRAM_QUANTUM);
							else incrementErrorRateHistogram_trivial(tmpWindows,lastTmpWindow,errorRateMatrix[p],MIN_ERROR,MAX_ERROR,N_HISTOGRAM_BINS,HISTOGRAM_QUANTUM);
						}
						else {
							if (PERIOD!=0) incrementErrorRateHistogram_periodic(currentRead,SURFACE_THRESHOLD,tmpWindows,lastTmpWindow,errorRateMatrix[p],MIN_ERROR,MAX_ERROR,N_HISTOGRAM_BINS,HISTOGRAM_QUANTUM);
							else incrementErrorRateHistogram_trivial(tmpWindows,lastTmpWindow,errorRateMatrix[p],MIN_ERROR,MAX_ERROR,N_HISTOGRAM_BINS,HISTOGRAM_QUANTUM);
						}
					}
					currentRead=Alignments.readA-1; 
					if (Alignments.diffs<=FragmentsStep1.TOO_FEW_DIFFS) {
						if (Alignments.startA<=IDENTITY_THRESHOLD && Alignments.endA>=Reads.getReadLength(Alignments.readA-1)-IDENTITY_THRESHOLD) {
							System.err.println("ERROR: Fragment "+(Alignments.readA-1)+" is too similar to reference "+(Alignments.readB-1)+" ("+Alignments.diffs+" diffs):");
							System.err.println("   "+str);
							System.exit(1);
						}
						// A fragment might still contain the reference as a substring. We
						// ignore such an alignment.
						str=br.readLine();
						continue;
					}
					tmpWindows[0].set(Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.diffs);
					lastTmpWindow=0;
					str=br.readLine();
					continue;
				}
				if (Alignments.diffs<=FragmentsStep1.TOO_FEW_DIFFS) {
					if (Alignments.startA<=IDENTITY_THRESHOLD && Alignments.endA>=Reads.getReadLength(Alignments.readA-1)-IDENTITY_THRESHOLD) {
						System.err.println("ERROR: Fragment "+(Alignments.readA-1)+" is too similar to reference "+(Alignments.readB-1)+" ("+Alignments.diffs+" diffs):");
						System.err.println("   "+str);
						System.exit(1);
					}
					// A fragment might still contain the reference as a substring. We
					// ignore such an alignment.
					str=br.readLine();
					continue;
				}
				lastTmpWindow++;
				if (lastTmpWindow==tmpWindows.length) {
					FragmentsStep1.PermutationWindow[] newArray = new FragmentsStep1.PermutationWindow[tmpWindows.length+CAPACITY];
					System.arraycopy(tmpWindows,0,newArray,0,tmpWindows.length);
					for (j=tmpWindows.length; j<newArray.length; j++) newArray[j] = new FragmentsStep1.PermutationWindow();
					tmpWindows=newArray;
				}
				tmpWindows[lastTmpWindow].set(Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.diffs);
				str=br.readLine();
			}
			br.close();
			if (currentRead!=-1) {
				FragmentsStep3.alignmentStats(currentRead,tmpWindows,lastTmpWindow,fragmentStats,referenceLengths,referenceCoverage,MIN_ALIGNMENT_LENGTH_FOR_COVERAGE);
				fragmentsPerReference(tmpWindows,lastTmpWindow,fragmentReferenceStats);
				if (BASIN_TYPE!=Constants.INTERVAL_PERIODIC) {
					if (PATH_ID!=-1) incrementErrorRateHistogram_nonPeriodic(currentRead,SURFACE_THRESHOLD,tmpWindows,lastTmpWindow,errorRateMatrix[p],MIN_ERROR,MAX_ERROR,N_HISTOGRAM_BINS,HISTOGRAM_QUANTUM);
					else incrementErrorRateHistogram_trivial(tmpWindows,lastTmpWindow,errorRateMatrix[p],MIN_ERROR,MAX_ERROR,N_HISTOGRAM_BINS,HISTOGRAM_QUANTUM);
				}
				else {
					if (PERIOD!=0) incrementErrorRateHistogram_periodic(currentRead,SURFACE_THRESHOLD,tmpWindows,lastTmpWindow,errorRateMatrix[p],MIN_ERROR,MAX_ERROR,N_HISTOGRAM_BINS,HISTOGRAM_QUANTUM);
					else incrementErrorRateHistogram_trivial(tmpWindows,lastTmpWindow,errorRateMatrix[p],MIN_ERROR,MAX_ERROR,N_HISTOGRAM_BINS,HISTOGRAM_QUANTUM);
				}
			}
			for (j=0; j<fragmentStats.length; j++) totalFragmentStats[j]+=fragmentStats[j];
			for (j=0; j<nReferences; j++) {
				sum=fragmentReferenceStats[j][0]+fragmentReferenceStats[j][1];
				fragmentReferenceHistogram.write(((double)fragmentReferenceStats[j][0])/sum+","+((double)fragmentReferenceStats[j][1])/sum+"\n");
			}
			for (j=0; j<N_HISTOGRAM_BINS; j++) erHistogram[j]+=errorRateMatrix[p][j];
			nReferencesWithLowCoverage+=checkConsensusCoverage(id,nReferences,referenceLengths,referenceCoverage);
			if ( BASIN_TYPE==Constants.INTERVAL_PERIODIC && PERIOD>0 && (PERIOD)<MIN_PERIOD_LONG && 
			     nAlignments<nFragments*FragmentsStep1.N_FRAGMENTS_REFERENCE_ALIGNMENTS_THRESHOLD
			   )  System.err.println("WARNING: Too few fragment-reference alignments in short-period repeat "+id+": k-mer capping by the aligner?");
			
			// Comparing consensus alignments to dominator alignments
			isShortPeriod=(BASIN_TYPE==Constants.INTERVAL_PERIODIC)&&((PERIOD)<MIN_PERIOD_LONG);
			compareConsensusDominatorAlignments(isShortPeriod,FR_ALIGNMENTS_FILE,FC_ALIGNMENTS_FORMAT,DOMINATOR_ALIGNMENTS_DIR+"/"+list[i]+"/LAshow-fragments-reference.txt",FD_ALIGNMENTS_FORMAT,list[i].substring(IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH),nFragments,nReferences,MIN_ALIGNMENT_LENGTH_IN_INFERENCE,out,CONSENSUS_DOMINATOR_VERBOSE);
			consensusDominatorDetails[p][0]=((double)out[0])/out[1];
			consensusDominatorDetails[p][1]=isShortPeriod?-1:((double)out[7])/out[8];
			consensusDominatorDetails[p][2]=out[4];
			consensusDominatorDetails[p][3]=isShortPeriod?-1:out[5];
			consensusDominatorDetails[p][4]=isShortPeriod?-1:out[6];
			consensusDominatorDetails[p][5]=nFragments;
			if (!isShortPeriod) {
				for (j=0; j<consensusDominatorStats.length; j++) consensusDominatorStats[j]+=out[j];
			}
			
			// Checking whole block-consensus alignments
			if (BLOCK_CONSENSUS_ANALYSIS) {
				Reads.nReads=N_READS;
				Reads.readLengths=readLengthsBackup;
				Reads.readIDs=readIDsBackup;
				Reads.firstRead=readIDsBackup[0]; Reads.lastRead=readIDsBackup[N_READS-1]; 
				analyzeBlockAlignments(CONSENSI_DIR+"/"+list[i]+"/LAshow-block-consensus.txt",BASINS_FILE,1500,ERROR_RATES,tmpIO,qualityHistogram,trackIntersection,MIN_LOW_QUALITY);
				for (j=0; j<ERROR_RATES.length; j++) {
					for (k=0; k<9; k++) blockReferenceHistogram.write(tmpIO[9*j+k]+",");
				}
				blockReferenceHistogram.write("\n");
				blockReferenceHistogramIDs.write(id+"\n");
				for (j=0; j<qualityHistogram.length; j++) {
					for (k=0; k<qualityHistogram[j].length; k++) blockReferenceQualities.write(qualityHistogram[j][k]+",");
				}
				blockReferenceQualities.write("\n");
				if (TRACK_FILE!=null) {
					for (j=0; j<ERROR_RATES.length; j++) trackIntersectionFile.write(trackIntersection[j][0]+","+trackIntersection[j][1]+","+trackIntersection[j][2]+","+tmpIO[9*j+6]+",");
					trackIntersectionFile.write("\n");
				}
				blockCoverage(CONSENSI_DIR+"/"+list[i]+"/LAshow-block-consensus.txt",MIN_ALIGNMENT_LENGTH_FOR_COVERAGE,nReferences,referenceLengths,coverageFromBlock);
				for (j=0; j<nReferences; j++) {
					bw = new BufferedWriter(new FileWriter(STATS_DIR+"/blockCoverage-"+id+"-"+j+".txt"));
					for (k=0; k<referenceLengths[j]; k++) bw.write(coverageFromBlock[j][k]+"\n");
					bw.close();
				}
			}

			// Printing statistics
			bw = new BufferedWriter(new FileWriter(STATS_DIR+"/fragmentStats-"+id+".txt"));
			for (j=0; j<fragmentStats.length; j++) bw.write(fragmentStats[j]+",");
			bw.close();
			for (j=0; j<nReferences; j++) {
				bw = new BufferedWriter(new FileWriter(STATS_DIR+"/referenceCoverage-"+id+"-"+j+".txt"));
				for (k=0; k<referenceLengths[j]; k++) bw.write(referenceCoverage[j][k]+"\n");
				bw.close();
			}
		}
		lengthHistogram.close(); lengthHistogramLongPeriod.close(); lengthHistogramUnknownPeriod.close();
		fragmentReferenceHistogram.close(); 
		blockReferenceHistogram.close();
		blockReferenceHistogramIDs.close();
		blockReferenceHistogramPeriodic.close();
		blockReferenceQualities.close();
		if (TRACK_FILE!=null) trackIntersectionFile.close();
		
		// Computing utility
		if (BLOCK_CONSENSUS_ANALYSIS) getUtility(CONSENSI_DIR,list,ERROR_RATES,UTILITIES_FILE_C,UTILITIES_FILE_R);
		
		// Displaying global statistics
		System.out.println();
		System.out.println("Overall statistics:  ("+totalNFragments+" fragments, "+totalNReferences+" references)");
		System.out.println("Consensi with low-coverage substrings: "+nReferencesWithLowCoverage+" ("+IO.format(100*((double)nReferencesWithLowCoverage)/totalNReferences)+"%)");
		System.out.println("Too long fragments: "+totalFragmentStats[6]+" ("+IO.format(100*((double)totalFragmentStats[6])/totalNFragments)+"%)");
		nFragments=totalNFragments-totalFragmentStats[0]-totalFragmentStats[1]-totalFragmentStats[2]-totalFragmentStats[3]-totalFragmentStats[4]-totalFragmentStats[5];
		System.out.println("Fragments with no alignment: "+nFragments+" ("+IO.format(100*((double)nFragments)/totalNFragments)+"%)");
		System.out.println("Fragments with holes: "+totalFragmentStats[5]+" ("+IO.format(100*((double)totalFragmentStats[5])/totalNFragments)+"%)");
		System.out.println("Type 4 fragments: "+totalFragmentStats[4]+" ("+IO.format(100*((double)totalFragmentStats[4])/totalNFragments)+"%)");
		System.out.println("Type 3 fragments: "+totalFragmentStats[3]+" ("+IO.format(100*((double)totalFragmentStats[3])/totalNFragments)+"%)");
		System.out.println("Type 2 fragments: "+totalFragmentStats[2]+" ("+IO.format(100*((double)totalFragmentStats[2])/totalNFragments)+"%)");
		System.out.println("Type 1 fragments: "+totalFragmentStats[1]+" ("+IO.format(100*((double)totalFragmentStats[1])/totalNFragments)+"%)");
		System.out.println("Type 0 fragments: "+totalFragmentStats[0]+" ("+IO.format(100*((double)totalFragmentStats[0])/totalNFragments)+"%)");
		System.out.println();
		System.out.println("Comparison of consensus/dominator alignments (excluding short-period):");
		System.out.println("Alignments to dominators: "+consensusDominatorStats[0]);
		System.out.println("Alignments to consensi: "+consensusDominatorStats[1]+" ("+IO.format(((double)consensusDominatorStats[1])/consensusDominatorStats[0])+" times those to dominators)");
		System.out.println("Fragments that align to a dominator: "+consensusDominatorStats[2]+" ("+IO.getPercent(consensusDominatorStats[2],totalNFragmentsNonShortPeriod)+"% of all fragments)");
		System.out.println("Fragments that align to a consensus: "+consensusDominatorStats[3]+" ("+IO.getPercent(consensusDominatorStats[3],totalNFragmentsNonShortPeriod)+"% of all fragments)");
		System.out.println("Fragments that align to a dominator but do not align to a consensus: "+consensusDominatorStats[4]+" ("+IO.getPercent(consensusDominatorStats[4],consensusDominatorStats[2])+"% of those to dominators)");
		System.out.println("Fragments with consensus error rate worse than dominator error rate: "+consensusDominatorStats[5]+" ("+IO.getPercent(consensusDominatorStats[5],consensusDominatorStats[2]-consensusDominatorStats[4])+"% of those to both)");
		System.out.println("Fragments with consensus surface smaller than dominator surface: "+consensusDominatorStats[6]+" ("+IO.getPercent(consensusDominatorStats[6],consensusDominatorStats[2]-consensusDominatorStats[4])+"%)");
		System.out.println("Total dominator surface: "+consensusDominatorStats[7]);
		System.out.println("Total consensus surface: "+consensusDominatorStats[8]+" ("+IO.format(((double)(consensusDominatorStats[8]))/consensusDominatorStats[7])+" times)");
		System.err.println();
		System.err.println("Consensus-dominator details (including short-period):");
		System.err.println("IDS:");
		j=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			j++;
			System.err.println(j+": "+list[i].substring(IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH));
		}
		System.err.println();
		System.err.println("nConsensusAlignments/nDominatorAlignments:");
		j=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			j++;
			System.err.println(j+", "+IO.format(1.0/consensusDominatorDetails[j][0]));
		}
		System.err.println();
		System.err.println("consensusSurface/dominatorSurface:");
		j=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			j++;
			System.err.println(j+", "+IO.format(1.0/consensusDominatorDetails[j][1]));
		}
		System.err.println();
		System.err.println("fragments that do not align to consensus:");
		j=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			j++;
			System.err.println(j+", "+IO.format(consensusDominatorDetails[j][2])+", "+consensusDominatorDetails[j][5]);
		}
		System.err.println();
		System.err.println("fragments with worse error rate in consensus:");
		j=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			j++;
			System.err.println(j+", "+IO.format(consensusDominatorDetails[j][3])+", "+consensusDominatorDetails[j][5]);
		}
		System.err.println();
		System.err.println("fragments with smaller surface in consensus:");
		j=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			j++;
			System.err.println(j+", "+IO.format(consensusDominatorDetails[j][4])+", "+consensusDominatorDetails[j][5]);
		}
		System.err.println();
		
		// Printing error rate histograms
		bw = new BufferedWriter(new FileWriter(ERROR_RATE_HISTOGRAM_FILE));
		for (i=0; i<N_HISTOGRAM_BINS; i++) bw.write((MIN_ERROR+i*HISTOGRAM_QUANTUM)+","+errorRateHistogram[i]+"\n");
		bw.close();
		bw = new BufferedWriter(new FileWriter(ERROR_RATE_HISTOGRAM_SHORTPERIOD_FILE));
		for (i=0; i<N_HISTOGRAM_BINS; i++) bw.write((MIN_ERROR+i*HISTOGRAM_QUANTUM)+","+errorRateHistogramShortPeriod[i]+"\n");
		bw.close();
		bw = new BufferedWriter(new FileWriter(ERROR_RATE_HISTOGRAM_LONGPERIOD_FILE));
		for (i=0; i<N_HISTOGRAM_BINS; i++) bw.write((MIN_ERROR+i*HISTOGRAM_QUANTUM)+","+errorRateHistogramLongPeriod[i]+"\n");
		bw.close();
		bw = new BufferedWriter(new FileWriter(ERROR_RATE_HISTOGRAM_UNKNOWNPERIOD_FILE));
		for (i=0; i<N_HISTOGRAM_BINS; i++) bw.write((MIN_ERROR+i*HISTOGRAM_QUANTUM)+","+errorRateHistogramUnknownPeriod[i]+"\n");
		bw.close();
		bw = new BufferedWriter(new FileWriter(ERROR_RATE_HISTOGRAM_CYCLIC_FILE));
		for (i=0; i<N_HISTOGRAM_BINS; i++) bw.write((MIN_ERROR+i*HISTOGRAM_QUANTUM)+","+errorRateHistogramCyclic[i]+"\n");
		bw.close();
		bw = new BufferedWriter(new FileWriter(ERROR_RATE_MATRIX_FILE));
		for (i=0; i<nDescriptors; i++) {
			for (j=0; j<N_HISTOGRAM_BINS; j++) bw.write(errorRateMatrix[i][j]+",");
			bw.write("\n");
		}
		bw.close();
		
		// Printing consensus-consensus alignment graph
		consensusAlignmentsGraph(true,CONSENSUS_CONSENSUS_ALIGNMENTS,totalNReferences,REFERENCE_LENGTHS_DIR,REFERENCE_LENGTHS_DIR+"/list-consensi.txt",MIN_ALIGNMENT_LENGTH_FOR_COVERAGE,REFERENCE_LENGTHS_DIR+"/connection.txt",REFERENCE_LENGTHS_DIR+"/consensus-consensus.dot");
		
		// Checking repeat track built from mapping the consensi.
		if (CONSENSI_REPEAT_TRACK!=null) checkConsensiRepeatTrack(STEP1_DIR,REPEAT_TRACK_INPUT_FILE,REPEAT_TRACK_INPUT_FORMAT,CONSENSI_REPEAT_TRACK,CONSENSI_REPEAT_TRACK_FORMAT,N_READS,READ_LENGTHS_FILE,READ_IDS_FILE,QUALITY_THRESHOLDS_FILE,QUALITIES_FILE);
	}
	
	
	/**
	 * Checks the coverage of every consensus of basin $id$, using a sliding window.
	 *
	 * @return number of consensi with a window of low coverage.
	 */
	private static final int checkConsensusCoverage(String id, int nReferences, int[] referenceLengths, int[][] referenceCoverage) {
		final int WINDOW_LENGTH = 1000;  // Arbitrary
		final int MIN_CONSENSUS_COVERAGE = 3;  // Arbitrary
		final int MIN_COVERAGE = MIN_CONSENSUS_COVERAGE*WINDOW_LENGTH;
		int i, j;
		int length, sum, minSum, minEnd, out;
		
		out=0;
		for (i=0; i<nReferences; i++) {
			sum=0;
			for (j=0; j<WINDOW_LENGTH; j++) sum+=referenceCoverage[i][j];
			minSum=sum; minEnd=WINDOW_LENGTH-1; length=referenceLengths[i];
			for (j=WINDOW_LENGTH; j<length; j++) {
				sum+=referenceCoverage[i][j]-referenceCoverage[i][j-WINDOW_LENGTH];
				if (sum<minSum) { minSum=sum; minEnd=j; }
			}
			if (minSum<MIN_COVERAGE) {
				System.err.println("ERROR: Basin "+id+", consensus "+i+"/"+nReferences+", window end "+minEnd+": avg. coverage very low ("+(((double)minSum)/WINDOW_LENGTH)+").");
				out++;
			}
		}
		return out;
	}
	
	
	
	
	// -------------------------------- ERROR RATE ---------------------------------------

	/**	
	 * Assume that a consensus covers at least the $surfaceThreshold$ fraction of fragment 
	 * $fragmentID$ with alignments in the same orientation. The procedure finds the 
	 * consensus with min avg error rate of its alignments, and it uses just the error 
	 * rates of those alignments to increment $histogram$. If no such consensus can be
	 * found, the procedure does not increment $histogram$ and returns FALSE.
	 *
	 * Remark: the procedure uses alignments of any length.
	 */
	private static final boolean incrementErrorRateHistogram_nonPeriodic(int fragmentID, double surfaceThreshold, FragmentsStep1.PermutationWindow[] windows, int lastWindow, int[] histogram, double minHistogramError, double maxHistogramError, int nHistogramBins, double histogramQuantum) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_SURFACE = (int)(Reads.getReadLength(fragmentID)*surfaceThreshold);
		boolean orientation, currentOrientation, minOrientation;
		int i;
		int readB, startA, endA, currentReadB, currentStartA, currentEndA;
		int currentSurface, currentDenominator, minReadB;
		double avgError, minError, currentNumerator;

		// Finding the best consensus string and orientation
		FragmentsStep1.PermutationWindow.order=FragmentsStep1.PermutationWindow.ORDER_READB_ORIENTATION_STARTA_STARTB;
		if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
		currentReadB=-1; currentOrientation=false; currentStartA=-1; currentEndA=-1;
		currentSurface=0; currentNumerator=0.0; currentDenominator=0;
		minError=Math.POSITIVE_INFINITY; minReadB=-1; minOrientation=false;
		for (i=0; i<=lastWindow; i++) {
			readB=windows[i].readB; orientation=windows[i].orientation;
			startA=windows[i].startA; endA=windows[i].endA;
			if (readB!=currentReadB || orientation!=currentOrientation) {
				currentSurface+=currentEndA-currentStartA+1;
				if (currentReadB!=-1 && currentSurface>=MIN_SURFACE) {
					avgError=currentNumerator/currentDenominator;
					if (avgError<minError) {
						minError=avgError;
						minReadB=currentReadB; minOrientation=currentOrientation;
					}
				}
				currentReadB=readB; currentOrientation=orientation;
				currentStartA=startA; currentEndA=endA; currentSurface=0;
				currentNumerator=((double)(windows[i].diffs<<1))/(endA-startA+windows[i].endB-windows[i].startB+2);
				currentDenominator=1;
				continue;
			}
			if (startA>currentEndA+IDENTITY_THRESHOLD) {
				currentSurface+=currentEndA-currentStartA+1;
				currentStartA=startA; currentEndA=endA;
			}
			else currentEndA=Math.max(currentEndA,endA);
			currentNumerator+=((double)(windows[i].diffs<<1))/(endA-startA+windows[i].endB-windows[i].startB+2);
			currentDenominator++;
		}
		currentSurface+=currentEndA-currentStartA+1;
		if (currentReadB!=-1 && currentSurface>=MIN_SURFACE) {
			avgError=currentNumerator/currentDenominator;
			if (avgError<minError) {
				minError=avgError;
				minReadB=currentReadB; minOrientation=currentOrientation;
			}
		}
		
		// Incrementing histogram
		if (minReadB==-1) return false;
		for (i=0; i<=lastWindow; i++)  {
			if (windows[i].readB<minReadB || (windows[i].readB==minReadB && windows[i].orientation && !minOrientation)) continue;
			if (windows[i].readB>minReadB) break;
			avgError=((double)(windows[i].diffs<<1))/(windows[i].endA-windows[i].startA+windows[i].endB-windows[i].startB+2);
			if (avgError<minHistogramError) histogram[0]++;
			else if (avgError>=maxHistogramError) histogram[nHistogramBins-1]++;
			else histogram[(int)((avgError-minHistogramError)/histogramQuantum)]++;
		}
		return true;
	}
	
	
	/**
	 * Let's call "red" an alignment A (to any consensus) such that its interval in the 
	 * fragment is not covered by a set of other alignments (to any consensus, possibly to
	 * distinct consensi) with lower error rate than A. The procedure adds to $histogram$ 
	 * only red alignments. This is done because the period of a long periodic fragment 
	 * might slightly change over the string, so disjoint substrings of the fragment might
	 * align better to different consensi.
	 *
	 * Remark: this applies both to short and long periods. The consensus of a long-period
	 * kernel is likely as long as the period itself, and distinct consensi should 
	 * represent variants of the period. The consensus of a short-period kernel is much 
	 * longer than the period, and a single consensus might contain several distinct 
	 * variants of the period.
	 *
	 * @return FALSE iff less that $surfaceThreshold$ fraction of $fragmentID$'s surface 
	 * is covered by alignments. We take this as an indication that the fragment contains 
	 * variants of the period that might align with high error to the consensi, so we do
	 * not update $histogram$ in this case.
	 */
	private static final boolean incrementErrorRateHistogram_periodic(int fragmentID, double surfaceThreshold, FragmentsStep1.PermutationWindow[] windows, int lastWindow, int[] histogram, double minHistogramError, double maxHistogramError, int nHistogramBins, double histogramQuantum) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MAX_GAP = IO.quantum<<1;
		final int MIN_SURFACE = (int)(Reads.getReadLength(fragmentID)*surfaceThreshold);
		boolean found;
		int i, j;
		int startA, endA, currentStartA, currentEndA, surface;
		double errorI, errorJ;
		
		FragmentsStep1.PermutationWindow.order=FragmentsStep1.PermutationWindow.ORDER_STARTA;
		if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
		
		// Checking if enough surface is covered
		currentStartA=-1; currentEndA=-1; surface=0;
		for (i=0; i<=lastWindow; i++) {
			startA=windows[i].startA; endA=windows[i].endA;
			if (currentStartA==-1) {
				currentStartA=startA; currentEndA=endA;
				continue;
			}
			else if (startA>currentEndA+IDENTITY_THRESHOLD) {
				surface+=currentEndA-currentStartA+1;
				currentStartA=startA; currentEndA=endA;
				continue;
			}
			else currentEndA=Math.max(currentEndA,endA);
		}
		surface+=currentEndA-currentStartA+1;
		if (surface<MIN_SURFACE) return false;
		
		// Incrementing histogram
		for (i=0; i<=lastWindow; i++) {
			startA=windows[i].startA; endA=windows[i].endA;
			errorI=((double)(windows[i].diffs<<1))/(endA-startA+windows[i].endB-windows[i].startB+2);
			found=true; currentStartA=-1; currentEndA=-1;
			for (j=0; j<=lastWindow; j++) {
				if (windows[j].startA>=endA-IDENTITY_THRESHOLD) break;
				if (windows[j].endA<=startA+IDENTITY_THRESHOLD) continue;
				errorJ=((double)(windows[j].diffs<<1))/(windows[j].endA-windows[j].startA+windows[j].endB-windows[j].startB+2);
				if (errorJ>=errorI) continue;
				if (currentStartA==-1) {
					currentStartA=windows[j].startA; currentEndA=windows[j].endA;
					if (currentStartA>startA+MAX_GAP) {
						found=false;
						break;
					}
				}
				else if (windows[j].startA>currentEndA+MAX_GAP) {
					found=false;
					break;
				}
				else currentEndA=Math.max(currentEndA,windows[j].endA);
			}
			if (found && (currentStartA==-1 || currentEndA<endA-MAX_GAP)) found=false;
			if (found) continue;
			if (errorI<minHistogramError) histogram[0]++;
			else if (errorI>=maxHistogramError) histogram[nHistogramBins-1]++;
			else histogram[(int)((errorI-minHistogramError)/histogramQuantum)]++;
		}
		return true;
	}

	
	/**
	 * Uses all windows to increment $histogram$.
	 */
	private static final void incrementErrorRateHistogram_trivial(FragmentsStep1.PermutationWindow[] windows, int lastWindow, int[] histogram, double minHistogramError, double maxHistogramError, int nHistogramBins, double histogramQuantum) {
		int i;
		double error;
		
		for (i=0; i<=lastWindow; i++) {
			error=((double)(windows[i].diffs<<1))/(windows[i].endA-windows[i].startA+windows[i].endB-windows[i].startB+2);
			if (error<minHistogramError) histogram[0]++;
			else if (error>=maxHistogramError) histogram[nHistogramBins-1]++;
			else histogram[(int)((error-minHistogramError)/histogramQuantum)]++;
		}
	}
	
	
	





	// ------------------------------ FALSE POSITIVES ------------------------------------
	
	/**
	 * Loads in global variable $pairsAllBasins[0..X]$ the sorted and compacted set of all 
	 * intervals in all basins, where $X$ is returned in output.
	 */
	private static final int loadAllBasinIntervals(String[] list, String step4dir) throws IOException {
		int i, j;
		int lastPair;
		String str;
		BufferedReader basinDescriptor;
		
		lastPair=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			basinDescriptor = new BufferedReader(new FileReader(step4dir+"/"+IO.BASIN_PREFIX+list[i].substring(IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH)+".txt"));
			str=basinDescriptor.readLine();  // Skipping header
			str=basinDescriptor.readLine(); 
			while (str!=null) {
				IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
				lastPair++;
				if (lastPair==pairsAllBasins.length) {
					Pair[] newPairs = new Pair[pairsAllBasins.length<<1];
					System.arraycopy(pairsAllBasins,0,newPairs,0,pairsAllBasins.length);
					for (j=pairsAllBasins.length; j<newPairs.length; j++) newPairs[j] = new Pair();
					pairsAllBasins=newPairs;
				}
				pairsAllBasins[lastPair].read=tmpArray[0];
				pairsAllBasins[lastPair].start=tmpArray[1];
				pairsAllBasins[lastPair].end=tmpArray[2];
				str=basinDescriptor.readLine();
			}
			basinDescriptor.close();
		}
		lastPair=compactPairs(pairsAllBasins,lastPair,false);
		return lastPair;
	}
	
	
	/**
	 * Loads in global variable $pairsAllIntervals[0..X]$ the sorted and compacted set of 
	 * all intervals detected by factorization (independent of whether they belong to a
	 * basin selected in Step 5 or not, or to a basin or not).
	 *
	 * @param qualityHistogram output array: histogram of quality values of all positions 
	 * inside an interval;
	 * @return the last element X of $pairsAllIntervals$.
	 */
	private static final int loadAllIntervals(String step1dir, int[] qualityHistogram) throws IOException {
		int i, j;
		int lastPair, from, to, value;
		String str;
		BufferedReader br;
		double[] qualities;
		
		lastPair=-1;
		
		// Alignment intervals
		br = new BufferedReader(new FileReader(step1dir+"/../intervals-alignments-all.txt"));
		str=br.readLine();
		while (str!=null) {
			lastPair++;
			if (lastPair==pairsAllIntervals.length) {
				Pair[] newPairs = new Pair[pairsAllIntervals.length<<1];
				System.arraycopy(pairsAllIntervals,0,newPairs,0,pairsAllIntervals.length);
				for (j=pairsAllIntervals.length; j<newPairs.length; j++) newPairs[j] = new Pair();
				pairsAllIntervals=newPairs;
			}
			i=str.indexOf(",");
			pairsAllIntervals[lastPair].read=Integer.parseInt(str.substring(0,i));
			AlignmentInterval.readAlignmentsFile(str,i+1,tmpArray,0);
			pairsAllIntervals[lastPair].start=tmpArray[1];
			pairsAllIntervals[lastPair].end=tmpArray[2];
			str=br.readLine();
		}
		br.close();
		
		// Dense intervals
		br = new BufferedReader(new FileReader(step1dir+"/../intervals-dense-all.txt"));
		str=br.readLine();
		while (str!=null) {
			lastPair++;
			if (lastPair==pairsAllIntervals.length) {
				Pair[] newPairs = new Pair[pairsAllIntervals.length<<1];
				System.arraycopy(pairsAllIntervals,0,newPairs,0,pairsAllIntervals.length);
				for (j=pairsAllIntervals.length; j<newPairs.length; j++) newPairs[j] = new Pair();
				pairsAllIntervals=newPairs;
			}
			i=str.indexOf(",");
			pairsAllIntervals[lastPair].read=Integer.parseInt(str.substring(0,i));
			DenseSubstring.readDenseSubstringsFile(str,i+1,tmpArray,0);
			pairsAllIntervals[lastPair].start=tmpArray[2];
			pairsAllIntervals[lastPair].end=tmpArray[3];
			str=br.readLine();
		}
		br.close();
		
		// Periodic intervals
		br = new BufferedReader(new FileReader(step1dir+"/../intervals-periodic-all.txt"));
		str=br.readLine();
		while (str!=null) {
			lastPair++;
			if (lastPair==pairsAllIntervals.length) {
				Pair[] newPairs = new Pair[pairsAllIntervals.length<<1];
				System.arraycopy(pairsAllIntervals,0,newPairs,0,pairsAllIntervals.length);
				for (j=pairsAllIntervals.length; j<newPairs.length; j++) newPairs[j] = new Pair();
				pairsAllIntervals=newPairs;
			}
			i=str.indexOf(",");
			pairsAllIntervals[lastPair].read=Integer.parseInt(str.substring(0,i));
			PeriodicSubstringInterval.readPeriodicSubstringsFile(str,i+1,tmpArray,0);
			pairsAllIntervals[lastPair].start=tmpArray[1];
			pairsAllIntervals[lastPair].end=tmpArray[2];
			str=br.readLine();
		}
		br.close();
		
		// Merging all intervals and computing quality histogram
		lastPair=compactPairs(pairsAllIntervals,lastPair,false);
		Math.set(qualityHistogram,qualityHistogram.length-1,0);
		for (i=0; i<=lastPair; i++) {
			qualities=Reads.getQualityArray(pairsAllIntervals[i].read);
			Reads.readEnds2qualityEnds(pairsAllIntervals[i].read,pairsAllIntervals[i].start,pairsAllIntervals[i].end,true,tmpArray);
			from=Math.min(tmpArray[0],qualities.length-1)+1;
			to=Math.min(tmpArray[1],qualities.length-1)-1;
			// Conservative from/to values, to completely exclude boundary
			// tracepoints from quality.
			for (j=from; j<=to; j++) {
				value=(int)qualities[j];
				qualityHistogram[value>=qualityHistogram.length?qualityHistogram.length-1:value]++;
			}
		}
		
		return lastPair;
	}
	
	
	/**
	 * Remark: the procedure uses the global variables $pairs1,pairs2,tmpArray$.
	 *
	 * @param alignmentsFile contains alignments between the consensi of a kernel and
	 * the entire set of reads S used for repeat inference; let A be the total surface of 
	 * S covered by such alignments (regardless of which consensus they use);
	 * @param basinDescriptor the basin descriptor file of a kernel; let B be the total 
	 * surface of S covered by intervals in $basinDescriptor$;
	 * @param errorRates for every $i$, use only alignments with error rate at most 
	 * $errorRates[i]$;
	 * @param out output array; for every cell $i$ of $errorRates$:
	 * out[9*i+0] = A \setminus B 
	 * out[9*i+1] = A \intersect B 
	 * out[9*i+2] = B \setminus A
	 * out[9*i+3] = A \setminus B' 
	 * out[9*i+4] = A \intersect B' 
	 * out[9*i+5] = B' \setminus A
	 * out[9*i+6] = A \setminus B''
	 * out[9*i+7] = A \intersect B''
	 * out[9*i+8] = B'' \setminus A
	 * where $B'$ is the total surface of S covered by intervals in \emph{any} basin 
	 * descriptor (which are assumed to be stored in $pairsAllBasins[0..
	 * lastPairAllBasins]$), and $B''$ is the total surface of S covered by any interval
	 * (which are assumed to be stored in $pairsAllIntervals[0..lastPairAllIntervals]$);
	 * @param qualityHistogram output array; cell $(i,j)$ contains the number of quality
	 * units of $A \setminus B''$ that have quality $j$ and belong to the $A$ computed 
	 * with the $i$-th error rate;
	 * @param trackIntersection if not NULL, stores at row $i$ the following statistics
	 * for the $i$-th error rate: (0) the sum of all intersections between $A 
	 * \setminus B''$ and $track$; (1) the sum of all intersections between $A \setminus 
	 * B''$ and regions with quality $>=minLowQuality$; (2) the intersection between (0) 
	 * and (1).
	 */
	private static final void analyzeBlockAlignments(String alignmentsFile, String basinDescriptorFile, int minAlignmentLength, double[] errorRates, long[] out, int[][] qualityHistogram, long[][] trackIntersection, byte minLowQuality) throws IOException {
		int i, j, k, h, p, r;
		int surface, intersection, lastPair1, lastPair2, lastPair3, lastLowQuality;
		long basinSurface, alignmentMinusBasin, alignmentIntersectBasin, basinMinusAlignment;
		int from, to, value, read;
		double error;
		String str;
		BufferedReader alignments, basinDescriptor;
		double[] qualities;
		
		// Loading basin descriptor intervals
		basinDescriptor = new BufferedReader(new FileReader(basinDescriptorFile));
		lastPair2=-1; basinSurface=0;
		str=basinDescriptor.readLine();  // Skipping header
		str=basinDescriptor.readLine(); 
		while (str!=null) {
			IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
			lastPair2++;
			if (lastPair2==pairs2.length) {
				Pair[] newPairs = new Pair[pairs2.length<<1];
				System.arraycopy(pairs2,0,newPairs,0,pairs2.length);
				for (i=pairs2.length; i<newPairs.length; i++) newPairs[i] = new Pair();
				pairs2=newPairs;
			}
			pairs2[lastPair2].read=tmpArray[0];
			pairs2[lastPair2].start=tmpArray[1];
			pairs2[lastPair2].end=tmpArray[2];
			str=basinDescriptor.readLine();
		}
		basinDescriptor.close();
		lastPair2=compactPairs(pairs2,lastPair2,false);
		for (i=0; i<=lastPair2; i++) basinSurface+=pairs2[i].end-pairs2[i].start+1;
		
		// Loading all alignment intervals
		alignments = new BufferedReader(new FileReader(alignmentsFile));
		lastPair3=-1; 
		str=alignments.readLine(); str=alignments.readLine();  // Skipping header
		str=alignments.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			if (Alignments.endA-Alignments.startA+1<minAlignmentLength) {
				str=alignments.readLine();
				continue;
			}
			lastPair3++;
			if (lastPair3==pairs3.length) {
				Pair[] newPairs = new Pair[pairs3.length<<1];
				System.arraycopy(pairs3,0,newPairs,0,pairs3.length);
				for (i=pairs3.length; i<newPairs.length; i++) newPairs[i] = new Pair();
				pairs3=newPairs;
			}
			pairs3[lastPair3].read=Alignments.readA-1;
			pairs3[lastPair3].start=Alignments.startA;
			pairs3[lastPair3].end=Alignments.endA;
			pairs3[lastPair3].error=((double)(Alignments.diffs<<1))/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2);
			pairs3[lastPair3].readB=Alignments.readB-1;
			str=alignments.readLine();
		}
		alignments.close();
		
		// Merging aligment intervals and basin intervals, at every error rate.
		Math.set(qualityHistogram,0);
		if (trackIntersection!=null) Math.set(trackIntersection,0L);
		for (r=0; r<errorRates.length; r++) {			
			error=errorRates[r];			
			lastPair1=-1;
			for (i=0; i<=lastPair3; i++) {
				if (pairs3[i].error>error) continue;
				lastPair1++;
				if (lastPair1==pairs1.length) {
					Pair[] newPairs = new Pair[pairs1.length<<1];
					System.arraycopy(pairs1,0,newPairs,0,pairs1.length);
					for (j=pairs1.length; j<newPairs.length; j++) newPairs[j] = new Pair();
					pairs1=newPairs;
				}
				pairs1[lastPair1].copyFrom(pairs3[i]);
			}
			lastPair1=compactPairs(pairs1,lastPair1,false);
			// Merging with the intervals of the current basin
			j=0; alignmentMinusBasin=0; alignmentIntersectBasin=0;
			for (i=0; i<=lastPair1; i++) {
				while (j<=lastPair2 && pairs2[j].read<pairs1[i].read) {
					j++;
					continue;
				}
				surface=pairs1[i].end-pairs1[i].start+1;
				if (j>lastPair2 || pairs2[j].read>pairs1[i].read) {
					alignmentMinusBasin+=surface;
					continue;
				}
				for (k=j; k<=lastPair2; k++) {
					if (pairs2[k].read!=pairs1[i].read) break;
					intersection=Intervals.intersectionLength(pairs2[k].start,pairs2[k].end,pairs1[i].start,pairs1[i].end);
					alignmentIntersectBasin+=intersection;
					surface-=intersection;
				}
				alignmentMinusBasin+=surface;
			}
			j=0; basinMinusAlignment=0;
			for (i=0; i<=lastPair2; i++) {
				while (j<=lastPair1 && pairs1[j].read<pairs2[i].read) {
					j++;
					continue;
				}
				surface=pairs2[i].end-pairs2[i].start+1;
				if (j>lastPair1 || pairs1[j].read>pairs2[i].read) {
					basinMinusAlignment+=surface;
					continue;
				}
				for (k=j; k<=lastPair1; k++) {
					if (pairs1[k].read!=pairs2[i].read) break;
					surface-=Intervals.intersectionLength(pairs1[k].start,pairs1[k].end,pairs2[i].start,pairs2[i].end);
				}
				basinMinusAlignment+=surface;
			}
			out[9*r+0]=alignmentMinusBasin; out[9*r+1]=alignmentIntersectBasin; out[9*r+2]=basinMinusAlignment;
			if (IO.CONSISTENCY_CHECKS && out[9*r+1]+out[9*r+2]!=basinSurface) {
				System.err.println("ERROR: (1)  out["+(9*r+1)+"]+out["+(9*r+2)+"]="+(out[9*r+1]+out[9*r+2])+" != basinSurface="+basinSurface+" :: "+basinDescriptorFile);
				System.exit(1);
			}
			// Merging with the intervals of every basin
			j=0; alignmentMinusBasin=0; alignmentIntersectBasin=0;
			for (i=0; i<=lastPair1; i++) {
				while (j<=lastPairAllBasins && pairsAllBasins[j].read<pairs1[i].read) {
					j++;
					continue;
				}
				surface=pairs1[i].end-pairs1[i].start+1;
				if (j>lastPairAllBasins || pairsAllBasins[j].read>pairs1[i].read) {
					alignmentMinusBasin+=surface;
					continue;
				}
				for (k=j; k<=lastPairAllBasins; k++) {
					if (pairsAllBasins[k].read!=pairs1[i].read) break;
					intersection=Intervals.intersectionLength(pairsAllBasins[k].start,pairsAllBasins[k].end,pairs1[i].start,pairs1[i].end);
					alignmentIntersectBasin+=intersection;
					surface-=intersection;
				}
				alignmentMinusBasin+=surface;
			}
			j=0; basinMinusAlignment=0;
			for (i=0; i<=lastPairAllBasins; i++) {
				while (j<=lastPair1 && pairs1[j].read<pairsAllBasins[i].read) {
					j++;
					continue;
				}
				surface=pairsAllBasins[i].end-pairsAllBasins[i].start+1;
				if (j>lastPair1 || pairs1[j].read>pairsAllBasins[i].read) {
					basinMinusAlignment+=surface;
					continue;
				}
				for (k=j; k<=lastPair1; k++) {
					if (pairs1[k].read!=pairsAllBasins[i].read) break;
					surface-=Intervals.intersectionLength(pairs1[k].start,pairs1[k].end,pairsAllBasins[i].start,pairsAllBasins[i].end);
				}
				basinMinusAlignment+=surface;
			}
			out[9*r+3]=alignmentMinusBasin; out[9*r+4]=alignmentIntersectBasin; out[9*r+5]=basinMinusAlignment;
			// Merging with the intervals of every repeat
			j=0; alignmentMinusBasin=0; alignmentIntersectBasin=0;
			for (i=0; i<=lastPair1; i++) {
				read=pairs1[i].read;
				while (j<=lastPairAllIntervals && pairsAllIntervals[j].read<read) {
					j++;
					continue;
				}
				surface=pairs1[i].end-pairs1[i].start+1;
				if (j>lastPairAllIntervals || pairsAllIntervals[j].read>read) {
					alignmentMinusBasin+=surface;
					qualities=Reads.getQualityArray(read);
					Reads.readEnds2qualityEnds(read,pairs1[i].start,pairs1[i].end,true,tmpArray);
					from=Math.min(tmpArray[0],qualities.length-1)+1;
					to=Math.min(tmpArray[1],qualities.length-1)-1;
					// Conservative from/to values, to completely exclude boundary
					// tracepoints from quality.
					for (h=from; h<=to; h++) {
						value=(int)qualities[h];
						qualityHistogram[r][value>=qualityHistogram[r].length?qualityHistogram[r].length-1:value]++;
					}
					if (trackIntersection!=null) {
						tmpArray[0]=pairs1[i].start; tmpArray[1]=pairs1[i].end;
						trackIntersection[r][0]+=Intervals.intersectionLength(tmpArray,1,track[read],lastTrack[read],1,true);
						lastLowQuality=Reads.getLowQualityIntervals(read,minLowQuality,tmpArray[0],tmpArray[1],tmpArray2);
						if (lastLowQuality>0) {
							for (h=0; h<lastLowQuality; h+=2) trackIntersection[r][1]+=tmpArray2[h+1]-tmpArray2[h]+1;
							trackIntersection[r][2]+=Intervals.intersectionLength(track[read],lastTrack[read],tmpArray2,lastLowQuality,1,true);
						}
					}
					continue;
				}
				p=pairs1[i].start; intersection=0;
				for (k=j; k<=lastPairAllIntervals; k++) {
					if (pairsAllIntervals[k].read!=read) break;
					intersection=Intervals.intersectionLength(pairsAllIntervals[k].start,pairsAllIntervals[k].end,pairs1[i].start,pairs1[i].end);
					if (intersection==0) continue;
					alignmentIntersectBasin+=intersection;
					surface-=intersection;
					if (pairsAllIntervals[k].start>p) {
						qualities=Reads.getQualityArray(read);
						Reads.readEnds2qualityEnds(read,p,pairsAllIntervals[k].start-1,true,tmpArray);
						from=Math.min(tmpArray[0],qualities.length-1)+1;
						to=Math.min(tmpArray[1],qualities.length-1)-1;
						// Conservative from/to values, to completely exclude boundary
						// tracepoints from quality.
						for (h=from; h<=to; h++) {
							value=(int)qualities[h];
							qualityHistogram[r][value>=qualityHistogram[r].length?qualityHistogram[r].length-1:value]++;
						}
						if (trackIntersection!=null) {
							tmpArray[0]=p; tmpArray[1]=pairsAllIntervals[k].start-1;
							trackIntersection[r][0]+=Intervals.intersectionLength(tmpArray,1,track[read],lastTrack[read],1,true);
							lastLowQuality=Reads.getLowQualityIntervals(read,minLowQuality,tmpArray[0],tmpArray[1],tmpArray2);
							if (lastLowQuality>0) {
								for (h=0; h<lastLowQuality; h+=2) trackIntersection[r][1]+=tmpArray2[h+1]-tmpArray2[h]+1;
								trackIntersection[r][2]+=Intervals.intersectionLength(track[read],lastTrack[read],tmpArray2,lastLowQuality,1,true);
							}
						}
					}
					p=pairsAllIntervals[k].end+1;
				}
				alignmentMinusBasin+=surface;
				if (intersection>0 && p<=pairs1[i].end) {
					qualities=Reads.getQualityArray(read);
					Reads.readEnds2qualityEnds(read,p,pairs1[i].end,true,tmpArray);
					from=Math.min(tmpArray[0],qualities.length-1)+1;
					to=Math.min(tmpArray[1],qualities.length-1)-1;
					// Conservative from/to values, to completely exclude boundary
					// tracepoints from quality.
					for (h=from; h<=to; h++) {
						value=(int)qualities[h];
						qualityHistogram[r][value>=qualityHistogram[r].length?qualityHistogram[r].length-1:value]++;
					}
					if (trackIntersection!=null) {
						tmpArray[0]=p; tmpArray[1]=pairs1[i].end;
						trackIntersection[r][0]+=Intervals.intersectionLength(tmpArray,1,track[read],lastTrack[read],1,true);
						lastLowQuality=Reads.getLowQualityIntervals(read,minLowQuality,tmpArray[0],tmpArray[1],tmpArray2);
						if (lastLowQuality>0) {
							for (h=0; h<lastLowQuality; h+=2) trackIntersection[r][1]+=tmpArray2[h+1]-tmpArray2[h]+1;
							trackIntersection[r][2]+=Intervals.intersectionLength(track[read],lastTrack[read],tmpArray2,lastLowQuality,1,true);
						}
					}
				}
			}
			j=0; basinMinusAlignment=0;
			for (i=0; i<=lastPairAllIntervals; i++) {
				while (j<=lastPair1 && pairs1[j].read<pairsAllIntervals[i].read) {
					j++;
					continue;
				}
				surface=pairsAllIntervals[i].end-pairsAllIntervals[i].start+1;
				if (j>lastPair1 || pairs1[j].read>pairsAllIntervals[i].read) {
					basinMinusAlignment+=surface;
					continue;
				}
				for (k=j; k<=lastPair1; k++) {
					if (pairs1[k].read!=pairsAllIntervals[i].read) break;
					surface-=Intervals.intersectionLength(pairs1[k].start,pairs1[k].end,pairsAllIntervals[i].start,pairsAllIntervals[i].end);
				}
				basinMinusAlignment+=surface;
			}
			out[9*r+6]=alignmentMinusBasin; out[9*r+7]=alignmentIntersectBasin; out[9*r+8]=basinMinusAlignment;	
		}
	}
	
	
	/**
	 * Sorts and merges all overlapping or adjacent elements of $pairs$.
	 * 
	 * @param readBmode TRUE=merges only pairs with the same readB, and uses readB as the
	 * first sort key; FALSE=merges all pairs and does not sort by readB;
	 * @return the new value of $lastPair$.
	 */
	private static final int compactPairs(Pair[] pairs, int lastPair, boolean readBmode) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i, j;
		int currentRead, currentReadB, currentStart, currentEnd;
		
		if (lastPair>0) {
			Pair.order=readBmode?Pair.ORDER_READB_READ_START:Pair.ORDER_READ_START;
			Arrays.sort(pairs,0,lastPair+1);
		}
		j=-1;
		currentRead=pairs[0].read; currentReadB=pairs[0].readB;
		currentStart=pairs[0].start; currentEnd=pairs[0].end;
		if (readBmode) {
			for (i=1; i<=lastPair; i++) {
				if (pairs[i].readB!=currentReadB || pairs[i].read!=currentRead || pairs[i].start>currentEnd+IDENTITY_THRESHOLD) {
					j++;
					pairs[j].read=currentRead; pairs[j].readB=currentReadB;
					pairs[j].start=currentStart; pairs[j].end=currentEnd;
					currentRead=pairs[i].read; currentReadB=pairs[i].readB;
					currentStart=pairs[i].start; currentEnd=pairs[i].end;
					continue;
				}
				currentEnd=Math.max(currentEnd,pairs[i].end);
			}
			j++;
			pairs[j].read=currentRead; pairs[j].readB=currentReadB;
			pairs[j].start=currentStart; pairs[j].end=currentEnd;
		}
		else {
			for (i=1; i<=lastPair; i++) {
				if (pairs[i].read!=currentRead || pairs[i].start>currentEnd+IDENTITY_THRESHOLD) {
					j++;
					pairs[j].read=currentRead; pairs[j].start=currentStart; pairs[j].end=currentEnd;
					currentRead=pairs[i].read; currentStart=pairs[i].start; currentEnd=pairs[i].end;
					continue;
				}
				currentEnd=Math.max(currentEnd,pairs[i].end);
			}
			j++;
			pairs[j].read=currentRead; pairs[j].start=currentStart; pairs[j].end=currentEnd;
		}
		return j;
	}
	
	
	public static class Pair implements Comparable {
		public static final int ORDER_READ_START = 0;
		public static final int ORDER_READB_READ_START = 1;
		public static int order;
		
		public int read, start, end;
		private double error;
		public int readB;
		public int id1, id2, id3;
		
		public Pair() {
			read=-1; start=-1; end=-1; error=0.0; readB=-1;
		}
		
		public Pair(int r, int s, int e, double er, int rb) {
			this.read=r; this.start=s; this.end=e; this.error=er; this.readB=rb;
		}
		
		public void copyFrom(Pair other) {
			this.read=other.read; this.start=other.start; this.end=other.end;
			this.error=other.error; this.readB=other.readB;
			this.id1=other.id1; this.id2=other.id2; this.id3=other.id3;
		}
		
		public int compareTo(Object other) {
			Pair otherPair = (Pair)other;
			if (order==ORDER_READB_READ_START) {
				if (readB<otherPair.readB) return -1;
				else if (readB>otherPair.readB) return 1;
			}
			if (read<otherPair.read) return -1;
			else if (read>otherPair.read) return 1;
			if (start<otherPair.start) return -1;
			else if (start>otherPair.start) return 1;
			return 0;
		}
		
		public String toString() {
			return "("+read+","+start+","+end+","+error+","+readB+")";
		}
	}
	
	
	/**
	 * Stores in $coverageFromBlock[i]$ the coverage histogram of the $i$-th consensus by 
	 * the alignments to the block stored in $alignmentsFile$.
	 */
	private static final void blockCoverage(String alignmentsFile, int minAlignmentLength, int nReferences, int[] referenceLengths, int[][] coverageFromBlock) throws IOException {
		int i;
		int from, to;
		String str;
		BufferedReader alignments;
		
		for (i=0; i<nReferences; i++) {
			if (coverageFromBlock[i].length<referenceLengths[i]) coverageFromBlock[i] = new int[referenceLengths[i]];
			Math.set(coverageFromBlock[i],referenceLengths[i]-1,0);
		}
		alignments = new BufferedReader(new FileReader(alignmentsFile));
		str=alignments.readLine(); str=alignments.readLine();  // Skipping header
		str=alignments.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			if (Alignments.endA-Alignments.startA+1<minAlignmentLength) {
				str=alignments.readLine();
				continue;
			}
			from=Math.min(Alignments.startB,referenceLengths[Alignments.readB-1]-1);
			to=Math.min(Alignments.endB,referenceLengths[Alignments.readB-1]-1);
			for (i=from; i<=to; i++) coverageFromBlock[Alignments.readB-1][i]++;
			str=alignments.readLine();
		}
		alignments.close();
	}
	
	
	/**
	 * Loads in $track,lastTrack$ the track in file $trackFile$, merging pairs of 
	 * consecutive intervals that are at distance $<=distanceThreshold$ from each other.
	 */
	public static final void loadTrack(String trackFile, int distanceThreshold, int[][] track, int[] lastTrack) throws IOException {
		final int MAX_READ_LENGTH = Reads.maxReadLength;
		final int INCREMENT = 100;  // Arbitrary, even.
		int i;
		int start, end, currentStart, currentEnd, lastInterval, readID;
		String str;
		BufferedReader br;
		String[] tokens;
		
		br = new BufferedReader(new FileReader(trackFile));
		readID=-1;
		str=br.readLine();
		while (str!=null) {
			readID++;
			if (str.length()==0) {
				lastTrack[readID]=-1;
				str=br.readLine();
				continue;
			}
			lastInterval=-1; currentStart=-1; currentEnd=-1;
			tokens=str.split(" ");
			for (i=0; i<tokens.length; i+=2) {
				start=Integer.parseInt(tokens[i]);
				end=Integer.parseInt(tokens[i+1]);
				if (currentStart==-1) { currentStart=start; currentEnd=end; }
				else {
					if (start<=currentEnd+distanceThreshold) currentEnd=end;
					else {
						if (lastInterval+2>=track[readID].length) {
							int[] newArray = new int[track[readID].length+INCREMENT];
							System.arraycopy(track[readID],0,newArray,0,track[readID].length);
							track[readID]=newArray;
						}
						track[readID][++lastInterval]=currentStart;
						track[readID][++lastInterval]=currentEnd;
						currentStart=start; currentEnd=end;
					}
				}
			}
			if (currentStart!=-1) {
				if (lastInterval+2>=track[readID].length) {
					int[] newArray = new int[track[readID].length+INCREMENT];
					System.arraycopy(track[readID],0,newArray,0,track[readID].length);
					track[readID]=newArray;
				}
				track[readID][++lastInterval]=currentStart;
				track[readID][++lastInterval]=currentEnd;
			}
			lastTrack[readID]=lastInterval;
			str=br.readLine();
		}
	}
	
	
	
	
	
	
	
	
	// ---------------------------------- UTILITY ----------------------------------------
	
	/**
	 * Computes the "utility" of a consensus, i.e. the surface of the whole block covered 
	 * by alignments with the consensus, and not covered by any other alignment with a 
	 * different consensus (possibly from the same repeat).
	 *
	 * @param list all files and subdirectories under $consensiDir$;
	 * @param outputFileConsensus contains a comma-separated row for every cell of 
	 * $errorRates$, storing the sequence of per-consesus utilities in the order induced 
	 * by $list$ and by the consensus IDs of each repeat;
	 * @param outputFileRepeat same as above, but for per-repeat utilities.
	 */
	private static final void getUtility(String consensiDir, String[] list, double[] errorRates, String outputFileConsensus, String outputFileRepeat) throws IOException {
		int i, j, r;
		int from, id1, id2, id3, nConsensi;
		BufferedWriter bwC, bwR;
		
		bwC = new BufferedWriter(new FileWriter(outputFileConsensus));
		bwR = new BufferedWriter(new FileWriter(outputFileRepeat));
		getUtility_loadAllAlignmentIntervals(consensiDir,list,errorRates);
		for (r=0; r<errorRates.length; r++) {
			System.err.println("getUtility> computing errorRate="+errorRates[r]+"...");
			from=0;
			id1=pairsAllAlignments[r][0].id1; 
			id2=pairsAllAlignments[r][0].id2;
			id3=pairsAllAlignments[r][0].id3;
			for (i=1; i<=lastPairAllAlignments[r]; i++) {
				if (pairsAllAlignments[r][i].id1!=id1 || pairsAllAlignments[r][i].id2!=id2 || pairsAllAlignments[r][i].id3!=id3) {
					nConsensi=getUtility_impl(r,from,i-1,tmpLong);
					for (j=0; j<nConsensi; j++) bwC.write(tmpLong[j]+",");
					bwR.write(tmpLong[nConsensi]+",");
					from=i;
					id1=pairsAllAlignments[r][i].id1; 
					id2=pairsAllAlignments[r][i].id2;
					id3=pairsAllAlignments[r][i].id3;
					continue;
				}
			}
			nConsensi=getUtility_impl(r,from,lastPairAllAlignments[r],tmpLong);
			for (j=0; j<nConsensi; j++) bwC.write(tmpLong[j]+",");
			bwR.write(tmpLong[nConsensi]+",");
			bwC.write("\n"); bwR.write("\n");
		}
		bwC.close(); bwR.close();
	}
		
		
	/**
	 * Remark: the procedure uses global variables $pairs1,pairs2$.
	 *
	 * @param r error rate ID;
	 * @param from,to all the alignment intervals in $pairsAllAlignments[r]$ coming from
	 * the same repeat (but possibly from multiple consensi);
	 * @param out output array: out[i]=difference between the surface covered by the i-th 
	 * consensus and the surface covered by any other consensus; the $nConsensi$-th cell 
	 * stores the same difference, but merging all consensi of the same repeat together;
	 * @return the number of consensi of the repeat.
	 */
	private static final int getUtility_impl(int r, int from, int to, long[] out) {
		int i, j, k, c;
		int nConsensi, lastPair1, lastPair2, surface, intersection, startI, endI;
		Pair pair;
		
		nConsensi=-1;
		for (i=from; i<=to; i++) nConsensi=Math.max(nConsensi,pairsAllAlignments[r][i].readB);
		nConsensi++;
		Math.set(out,nConsensi,0L);
		
		// 1. Per-consensus analysis
		for (c=0; c<nConsensi; c++) {
			// Compacting the intervals of the consensus
			lastPair1=-1;
			for (i=from; i<=to; i++) {
				pair=pairsAllAlignments[r][i];
				if (pair.readB!=c) continue;
				lastPair1++;
				if (lastPair1==pairs1.length) {
					Pair[] newPairs = new Pair[pairs1.length<<1];
					System.arraycopy(pairs1,0,newPairs,0,pairs1.length);
					for (j=pairs1.length; j<newPairs.length; j++) newPairs[j] = new Pair();
					pairs1=newPairs;
				}
				pairs1[lastPair1].copyFrom(pair);
			}
			lastPair1=compactPairs(pairs1,lastPair1,false);
			// Compacting the intervals that are not of the consensus
			lastPair2=-1;
			for (i=0; i<from; i++) {
				pair=pairsAllAlignments[r][i];
				lastPair2++;
				if (lastPair2==pairs2.length) {
					Pair[] newPairs = new Pair[pairs2.length<<1];
					System.arraycopy(pairs2,0,newPairs,0,pairs2.length);
					for (j=pairs2.length; j<newPairs.length; j++) newPairs[j] = new Pair();
					pairs2=newPairs;
				}
				pairs2[lastPair2].copyFrom(pair);
			}
			for (i=from; i<=to; i++) {
				pair=pairsAllAlignments[r][i];
				if (pair.readB==c) continue;
				lastPair2++;
				if (lastPair2==pairs2.length) {
					Pair[] newPairs = new Pair[pairs2.length<<1];
					System.arraycopy(pairs2,0,newPairs,0,pairs2.length);
					for (j=pairs2.length; j<newPairs.length; j++) newPairs[j] = new Pair();
					pairs2=newPairs;
				}
				pairs2[lastPair2].copyFrom(pair);
			}
			for (i=to+1; i<=lastPairAllAlignments[r]; i++) {
				pair=pairsAllAlignments[r][i];
				lastPair2++;
				if (lastPair2==pairs2.length) {
					Pair[] newPairs = new Pair[pairs2.length<<1];
					System.arraycopy(pairs2,0,newPairs,0,pairs2.length);
					for (j=pairs2.length; j<newPairs.length; j++) newPairs[j] = new Pair();
					pairs2=newPairs;
				}
				pairs2[lastPair2].copyFrom(pair);
			}
			lastPair2=compactPairs(pairs2,lastPair2,false);
			// Subtracting the two sets of intervals
			j=0;
			for (i=0; i<=lastPair1; i++) {
				while (j<=lastPair2 && pairs2[j].read<pairs1[i].read) {
					j++;
					continue;
				}
				surface=pairs1[i].end-pairs1[i].start+1;
				if (j>lastPair2 || pairs2[j].read>pairs1[i].read) {
					out[c]+=surface;
					continue;
				}
				startI=pairs1[i].start; endI=pairs1[i].end;
				for (k=j; k<=lastPair2; k++) {
					if (pairs2[k].read!=pairs1[i].read) break;
					intersection=Intervals.intersectionLength(pairs2[k].start,pairs2[k].end,startI,endI);
					surface-=intersection;
				}
				out[c]+=surface;
			}
		}
		
		// 2. Per-repeat analysis
		// Compacting all intervals of the repeat
		lastPair1=-1;
		for (i=from; i<=to; i++) {
			pair=pairsAllAlignments[r][i];
			lastPair1++;
			if (lastPair1==pairs1.length) {
				Pair[] newPairs = new Pair[pairs1.length<<1];
				System.arraycopy(pairs1,0,newPairs,0,pairs1.length);
				for (j=pairs1.length; j<newPairs.length; j++) newPairs[j] = new Pair();
				pairs1=newPairs;
			}
			pairs1[lastPair1].copyFrom(pair);
		}
		lastPair1=compactPairs(pairs1,lastPair1,false);
		// Compacting the intervals that are not of the repeat
		lastPair2=-1;
		for (i=0; i<from; i++) {
			pair=pairsAllAlignments[r][i];
			lastPair2++;
			if (lastPair2==pairs2.length) {
				Pair[] newPairs = new Pair[pairs2.length<<1];
				System.arraycopy(pairs2,0,newPairs,0,pairs2.length);
				for (j=pairs2.length; j<newPairs.length; j++) newPairs[j] = new Pair();
				pairs2=newPairs;
			}
			pairs2[lastPair2].copyFrom(pair);
		}
		for (i=to+1; i<=lastPairAllAlignments[r]; i++) {
			pair=pairsAllAlignments[r][i];
			lastPair2++;
			if (lastPair2==pairs2.length) {
				Pair[] newPairs = new Pair[pairs2.length<<1];
				System.arraycopy(pairs2,0,newPairs,0,pairs2.length);
				for (j=pairs2.length; j<newPairs.length; j++) newPairs[j] = new Pair();
				pairs2=newPairs;
			}
			pairs2[lastPair2].copyFrom(pair);
		}
		lastPair2=compactPairs(pairs2,lastPair2,false);
		// Subtracting the two sets of intervals
		j=0;
		for (i=0; i<=lastPair1; i++) {
			while (j<=lastPair2 && pairs2[j].read<pairs1[i].read) {
				j++;
				continue;
			}
			surface=pairs1[i].end-pairs1[i].start+1;
			if (j>lastPair2 || pairs2[j].read>pairs1[i].read) {
				out[nConsensi]+=surface;
				continue;
			}
			startI=pairs1[i].start; endI=pairs1[i].end;
			for (k=j; k<=lastPair2; k++) {
				if (pairs2[k].read!=pairs1[i].read) break;
				intersection=Intervals.intersectionLength(pairs2[k].start,pairs2[k].end,startI,endI);
				surface-=intersection;
			}
			out[nConsensi]+=surface;
		}
		return nConsensi;
	}
	
	
	/**
	 * Loads in global variable $pairsAllAlignments[r]$ the sorted and compacted set of 
	 * all alignment intervals in the entire block, with error rate at most 
	 * $errorRates[r]$. Intervals from the same file in $list$ are stored in a compact 
	 * block, but intervals are otherwise not guaranteed to be sorted and are not 
	 * compacted (they are just the original intervals from the files).
	 *
	 * @param list all files and subdirectories under $consensiDir$;
	 * @param errorRates decreasing values.
	 */
	private static final void getUtility_loadAllAlignmentIntervals(String consensiDir, String[] list, double[] errorRates) throws IOException {
		int i, j, k, p, q;
		int lastPair;
		double error;
		String str;
		BufferedReader alignmentsFile;
		Pair currentPair;
		
		// Loading pairs at the highest error
		lastPair=-1;
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			alignmentsFile = new BufferedReader(new FileReader(consensiDir+"/"+list[i]+"/LAshow-block-consensus.txt"));
			str=alignmentsFile.readLine();   // Skipping header
			str=alignmentsFile.readLine();
			str=alignmentsFile.readLine(); 
			while (str!=null) {
				Alignments.readAlignmentFile(str);
				error=((double)(Alignments.diffs<<1))/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2);
				if (error>errorRates[0]) {
					str=alignmentsFile.readLine();
					continue;
				}
				lastPair++;
				if (lastPair==pairsAllAlignments[0].length) {
					Pair[] newPairs = new Pair[pairsAllAlignments[0].length<<1];
					System.arraycopy(pairsAllAlignments[0],0,newPairs,0,pairsAllAlignments[0].length);
					for (j=pairsAllAlignments[0].length; j<newPairs.length; j++) newPairs[j] = new Pair();
					pairsAllAlignments[0]=newPairs;
				}
				pairsAllAlignments[0][lastPair].read=Alignments.readA-1;
				pairsAllAlignments[0][lastPair].start=Alignments.startA;
				pairsAllAlignments[0][lastPair].end=Alignments.endA;
				pairsAllAlignments[0][lastPair].error=error;
				pairsAllAlignments[0][lastPair].readB=Alignments.readB-1;
				p=IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH; q=list[i].indexOf('-',p+1);
				pairsAllAlignments[0][lastPair].id1=Integer.parseInt(list[i].substring(p,q));
				p=q+1; q=list[i].indexOf('-',p+1);
				pairsAllAlignments[0][lastPair].id2=Integer.parseInt(list[i].substring(p,q));
				pairsAllAlignments[0][lastPair].id3=Integer.parseInt(list[i].substring(q+1));
				str=alignmentsFile.readLine();
			}
			alignmentsFile.close();
		}
		
		// Creating all other pairs
		Math.set(lastPairAllAlignments,errorRates.length-1,-1);
		lastPairAllAlignments[0]=lastPair;
		for (i=0; i<=lastPair; i++) {
			currentPair=pairsAllAlignments[0][i];
			for (j=1; j<errorRates.length; j++) {
				if (currentPair.error<=errorRates[j]) {
					lastPairAllAlignments[j]++;
					if (lastPairAllAlignments[j]==pairsAllAlignments[j].length) {
						Pair[] newPairs = new Pair[pairsAllAlignments[j].length<<1];
						System.arraycopy(pairsAllAlignments[j],0,newPairs,0,pairsAllAlignments[j].length);
						for (k=pairsAllAlignments[j].length; k<newPairs.length; k++) newPairs[k] = new Pair();
						pairsAllAlignments[j]=newPairs;
					}
					pairsAllAlignments[j][lastPairAllAlignments[j]].copyFrom(currentPair);
				}
			}
		}
	}








	// ------------------------ CONSENSUS-CONSENSUS ALIGNMENTS ---------------------------
	
	/**
	 * Draws an undirected graph in which every node is a consensus string and every edge
	 * is an alignment in $alignmentsFile$. If $mode=TRUE$, only the following edges are 
	 * plotted: (1) any edge between two consensi of the same repeat; (2) identity edges 
	 * between two consensi of different repeats. This is because our pipeline explicitly 
	 * allows distinct repeats to share substrings, to overlap, and to be contained in one 
	 * another, so the only problematic case is when two repeats are approx. identical 
	 * end to end. If $mode=FALSE$, all edges are plotted.
	 *
	 * @param alignmentsFile all consensus-consensus alignments; it should have been 
	 * passed through $FilterAlignments$ to enforce a specific error rate and 
	 * minAlignmentLength;
	 * @param nConsensi total number of consensus strings over all repeats;
	 * @param consensusLengthsDir directory containing the length files of each repeat, 
	 * i.e. a file that stores the string length of every consensus of a repeat;	
	 * @param repeatIDsFile list of distinct repeat IDs;
	 * @param connectionFile output file: artificial connection file of the graph;
	 * @param dotFile output file: DOT representation of the graph.
	 */
	private static final void consensusAlignmentsGraph(boolean mode, String alignmentsFile, int nConsensi, String consensusLengthsDir, String repeatIDsFile, int minAlignmentLength, String connectionFile, String dotFile) throws IOException {
		final char CONSENSUS_SEPARATOR = '_';
		final double PEN_SCALE = 100.0;  // Edge width is % error rate.
		boolean sameRepeat;
		int i, j, k, p, q;
		int lastLength, readA, readI, previousReadI, lengthI, readJ, repeatID, id0, id1, id2, edgeSize, to;
		double errorRate;
		String str1, str2;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		String color;
		IntervalGraph.Node node;
		IntervalGraph.Edge edge;
		int[] consensusLengths;
		int[][] repeatIDs;
		
		IO.initialize();
		
		// Loading the length of every consensus
		consensusLengths = new int[nConsensi];
		lastLength=-1;
		br1 = new BufferedReader(new FileReader(repeatIDsFile));
		str1=br1.readLine();
		while (str1!=null) {
			p=str1.lastIndexOf(".");
			br2 = new BufferedReader(new FileReader(consensusLengthsDir+"/"+str1.substring(0,p)+"-lengths.txt"));
			str2=br2.readLine();
			while (str2!=null) {
				lastLength++;
				consensusLengths[lastLength]=Integer.parseInt(str2);
				str2=br2.readLine();
			}
			br2.close();
			str1=br1.readLine();
		}
		br1.close();
		
		// Loading the class of every consensus
		repeatIDs = new int[nConsensi][3];
		br1 = new BufferedReader(new FileReader(repeatIDsFile));
		str1=br1.readLine(); repeatID=-1;
		while (str1!=null) {
			p=str1.indexOf(CONSENSUS_SEPARATOR);
			id0=Integer.parseInt(str1.substring(0,p));
			q=str1.indexOf(CONSENSUS_SEPARATOR,p+1);
			id1=Integer.parseInt(str1.substring(p+1,q));
			p=q; q=str1.indexOf(".",p+1);
			id2=Integer.parseInt(str1.substring(p+1,q));
			q=str1.lastIndexOf(".");
			br2 = new BufferedReader(new FileReader(consensusLengthsDir+"/"+str1.substring(0,q)+"-lengths.txt"));
			str2=br2.readLine();
			while (str2!=null) {
				repeatID++;
				repeatIDs[repeatID][0]=id0; 
				repeatIDs[repeatID][1]=id1; 
				repeatIDs[repeatID][2]=id2;
				str2=br2.readLine();
			}
			br2.close();
			str1=br1.readLine();
		}
		br1.close();
		
		// Building artificial connection file
		br1 = new BufferedReader(new FileReader(alignmentsFile),IO.BUFFER_SIZE);
		bw = new BufferedWriter(new FileWriter(connectionFile),IO.BUFFER_SIZE);
		str1=br1.readLine(); str1=br1.readLine();  // Skipping header
		str1=br1.readLine();
		while (str1!=null) {
			Alignments.readAlignmentFile(str1);
			readA=Alignments.readA-1;
			bw.write(Constants.INTERVAL_ALIGNMENT+","+readA+",0,"+(consensusLengths[readA]-1)+",0,0,0,0,0,1,1\n");
			str1=br1.readLine();
		}
		br1.close(); bw.close();
		
		// Building the graph
		Alignments.minAlignmentLength=minAlignmentLength;
		Reads.loadReadIDs(0,nConsensi-1);
		Reads.readLengths=consensusLengths;
		IntervalGraph.maxAlignmentsPerRead=IntervalGraph.buildIntervalGraph(alignmentsFile,connectionFile,null,null,null,false,false,false,false,false,false,false,null,null,false);
		IntervalGraph.sortNodesArray();
		IntervalGraph.setAvgDiffs();
		bw = new BufferedWriter(new FileWriter(dotFile));
		bw.write("graph G {\n");
		previousReadI=-1;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			readI=IntervalGraph.nodesArray[i].read;
			while (previousReadI<readI-1) {
				previousReadI++;
				bw.write(previousReadI+" [class1=\""+repeatIDs[previousReadI][0]+"\",class2=\""+repeatIDs[previousReadI][0]+CONSENSUS_SEPARATOR+repeatIDs[previousReadI][1]+"\",class3=\""+repeatIDs[previousReadI][0]+CONSENSUS_SEPARATOR+repeatIDs[previousReadI][1]+CONSENSUS_SEPARATOR+repeatIDs[previousReadI][2]+"\"];\n");
			}
			bw.write(readI+" [class1=\""+repeatIDs[readI][0]+"\",class2=\""+repeatIDs[readI][0]+CONSENSUS_SEPARATOR+repeatIDs[readI][1]+"\",class3=\""+repeatIDs[readI][0]+CONSENSUS_SEPARATOR+repeatIDs[readI][1]+CONSENSUS_SEPARATOR+repeatIDs[readI][2]+"\"];\n");
			previousReadI=readI;
		}
		for (i=0; i<IntervalGraph.nNodes; i++) {
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) IntervalGraph.neighbors[i][j].printed=-1;
		}
		for (i=0; i<IntervalGraph.nNodes; i++) {
			readI=IntervalGraph.nodesArray[i].read;
			lengthI=IntervalGraph.nodesArray[i].length();
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
				edge=IntervalGraph.neighbors[i][j];
				to=edge.getTo(i);
				errorRate=(2*edge.avgDiffs)/(lengthI+IntervalGraph.nodesArray[to].length());
				if (edge.printed!=-1) continue;				
				readJ=IntervalGraph.nodesArray[to].read;
				sameRepeat=true;
				for (k=0; k<repeatIDs[readI].length; k++) {
					if (repeatIDs[readI][k]!=repeatIDs[readJ][k]) {
						sameRepeat=false;
						break;
					}
				}
				if (!mode || sameRepeat || edge.containment==Constants.CONTAINMENT_IDENTICAL) {
					if (edge.containment!=-1) color="green";
					else if (edge.overlap!=-1) color="blue";
					else if (edge.insertion!=-1) color="red";
					else if (edge.sharedSubstring!=-1) color="yellow";
					else color="gray";
					edgeSize=Math.max(1,(int)(edge.avgDiffs*PEN_SCALE));
					bw.write(readI+" -- "+readJ+" [color=\""+color+"\",weight=\""+edgeSize+"\",error=\""+IO.format(errorRate)+"\"];\n");
					edge.printed=1;
				}
			}
		}
		bw.write("}\n");
		bw.close();
	}





	
	
	
	// ------------------------------- SPECIFICITY ---------------------------------------

	/**
	 * @param windows $0..lastWindow$ is assumed to contain all windows of the same
	 * fragment;
	 * @param out for every reference (row): number of fragments that align only to that 
	 * reference (column 0), or that align also to other references (1).
	 */
	private static final void fragmentsPerReference(FragmentsStep1.PermutationWindow[] windows, int lastWindow, int[][] out) {
		int i, lastReference;
		
		lastReference=windows[0].readB;
		for (i=1; i<=lastWindow; i++) {
			if (windows[i].readB!=lastReference) {
				lastReference=-1;
				break;
			}
		}
		if (lastReference!=-1) {
			out[lastReference][0]++;
			return;
		}
		lastReference=windows[0].readB;
		for (i=1; i<=lastWindow; i++) {
			if (windows[i].readB!=lastReference) {
				out[lastReference][1]++;
				lastReference=windows[i].readB;
			}
		}
		out[lastReference][1]++;
	}
	
	
	
	
	
	
	
	
	// ------------------------- FRAG-CONSENSUS VS FRAG-DOMINATOR ------------------------
	
	/**
	 * Performs a detailed comparison between the alignments of fragments to their
	 * dominator sequences, and the alignments of the same fragments to their consensi.
	 *
	 * Remark: if the repeat is short-period, the n. of fragment-consensus alignments 
	 * might be much smaller than the n. of fragment-dominator alignments, since the 
	 * aligner might discard k-mers that are too frequent, and the consensus is more
	 * likely to contain frequent k-mers. For this reason, the error rate histogram for 
	 * short-period repeats might be distorted as well.
	 *	
	 * Remark: the procedure uses global variables $frags2cons,frags2dom,minError,
	 * consensusWindows,dominatorWindows,tmpArray$, which are assumed not to be null.
	 *
	 * @param nReferences number of reference sequences in the dominator and in the 
	 * consensus files;
   	 * @param *AlignmentsFormat format of LAshow files: TRUE=fragments and references are 
	 * in different DBs, and the IDs of both start from 1; FALSE=references are in block 1
	 * and fragments are in block 2 of the same DB, so fragment IDs start from 
	 * $nReferences$;
	 * @param minAlignmentLength the value used in repeat inference;
	 * @param out output array containing the number of:
	 * 0: dominator alignments;
	 * 1: consensus alignments; 
	 * 2: fragments that align to a dominator; 
	 * 3: fragments that align to a consensus;
	 * <-- (the procedure stops here if the repeat is short-period)
	 * 4: fragments that aligned to a dominator but do not align to a consensus;
	 * 5: fragments with consensus error rate worse than dominator error rate;
	 * 6: fragments with consensus surface smaller than dominator surface;
	 * 7: total surface that aligns to some dominator, over all fragments;
	 * 8: total surface that aligns to some consensus, over all fragments.
	 */
	private static final void compareConsensusDominatorAlignments(boolean isShortPeriod, String consensusAlignments, boolean consensusAlignmentsFormat, String dominatorAlignments, boolean dominatorAlignmentsFormat, String moduleID, int nFragments, int nReferences, int minAlignmentLength, long[] out, boolean verbose) throws IOException {
		final int SURFACE_THRESHOLD = minAlignmentLength>>2;  // Arbitrary
		final int CAPACITY = 1000;  // Arbitrary
		int i, j, n, p;
		int fragment, currentReadA, currentReadB, currentFragment, currentAlignments;
		int currentStart, currentEnd, currentSurface;
		int lastConsensusWindow, lastDominatorWindow, lastFrag2cons, lastFrag2dom;
		long sum;
		double error, currentError;
		String str;
		BufferedReader br;
		Math.set(out,out.length-1,0);
		
		// Loading frag-cons alignments
		br = new BufferedReader(new FileReader(consensusAlignments),IO.BUFFER_SIZE);
		lastConsensusWindow=-1;
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		n=0;
		while (str!=null) {
			n++;
			Alignments.readAlignmentFile(str);
			lastConsensusWindow++;
			if (lastConsensusWindow==consensusWindows.length) {
				FragConsAlignment[] newArray = new FragConsAlignment[consensusWindows.length<<1];
				System.arraycopy(consensusWindows,0,newArray,0,consensusWindows.length);
				for (j=consensusWindows.length; j<newArray.length; j++) newArray[j] = new FragConsAlignment();
				consensusWindows=newArray;
			}
			consensusWindows[lastConsensusWindow].set(Alignments.readA-1,Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.diffs);
			str=br.readLine();
		}
		br.close();
		out[1]=n;
		j=1; currentFragment=consensusWindows[0].readA;
		for (i=1; i<=lastConsensusWindow; i++) {
			fragment=consensusWindows[i].readA;
			if (fragment!=currentFragment) {
				j++;
				currentFragment=fragment;
			}
		}
		if (frags2cons.length<j) frags2cons = new int[j];
		currentFragment=consensusWindows[0].readA;
		j=0; frags2cons[0]=consensusAlignmentsFormat?currentFragment:currentFragment-nReferences;
		for (i=1; i<=lastConsensusWindow; i++) {
			fragment=consensusWindows[i].readA;
			if (fragment!=currentFragment) {
				frags2cons[++j]=consensusAlignmentsFormat?fragment:fragment-nReferences;
				currentFragment=fragment;
			}
		}
		lastFrag2cons=j; out[3]=lastFrag2cons+1;
		
		// Loading frag-dom alignments
		br = new BufferedReader(new FileReader(dominatorAlignments),IO.BUFFER_SIZE);
		lastDominatorWindow=-1;
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		n=0;
		while (str!=null) {
			n++;
			Alignments.readAlignmentFile(str);
			lastDominatorWindow++;
			if (lastDominatorWindow==dominatorWindows.length) {
				FragConsAlignment[] newArray = new FragConsAlignment[dominatorWindows.length<<1];
				System.arraycopy(dominatorWindows,0,newArray,0,dominatorWindows.length);
				for (j=dominatorWindows.length; j<newArray.length; j++) newArray[j] = new FragConsAlignment();
				dominatorWindows=newArray;
			}
			dominatorWindows[lastDominatorWindow].set(Alignments.readA-1,Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.diffs);
			str=br.readLine();
		}
		br.close();
		out[0]=n;
		j=1; currentFragment=dominatorWindows[0].readA;
		for (i=1; i<=lastDominatorWindow; i++) {
			fragment=dominatorWindows[i].readA;
			if (fragment!=currentFragment) {
				j++;
				currentFragment=fragment;
			}
		}
		if (frags2dom.length<j) frags2dom = new int[j];
		currentFragment=dominatorWindows[0].readA;
		j=0; frags2dom[0]=dominatorAlignmentsFormat?currentFragment:currentFragment-nReferences;
		for (i=1; i<=lastDominatorWindow; i++) {
			fragment=dominatorWindows[i].readA;
			if (fragment!=currentFragment) {
				frags2dom[++j]=dominatorAlignmentsFormat?fragment:fragment-nReferences;
				currentFragment=fragment;
			}
		}
		lastFrag2dom=j; out[2]=lastFrag2dom+1;
		if (isShortPeriod) return;
		
		// Checking existence of alignments
		if (tmpArray.length<lastFrag2dom+1) tmpArray = new int[lastFrag2dom+1];
		i=Math.setMinus(frags2dom,lastFrag2dom,frags2cons,lastFrag2cons,tmpArray);
		if (i!=-1) {
			out[4]+=i+1;
			if (verbose) {
				System.err.print("ERROR: module "+moduleID+": the following fragments (zero-based) aligned to some dominator, but do not align to any consensus: ");
				for (j=0; j<=i; j++) System.err.print(tmpArray[j]+", ");
				System.err.println();
			}
		}
		
		// Comparing avg. error rates
		if (minError[0].length<lastFrag2dom+1) minError = new double[2][lastFrag2dom+1];
		Math.set(minError[0],lastFrag2dom,Math.POSITIVE_INFINITY);
		Math.set(minError[1],lastFrag2dom,Math.POSITIVE_INFINITY);
		// Computing best avg. error of every fragment to a dominator
		j=0;
		currentReadA=dominatorWindows[0].readA; currentReadB=dominatorWindows[0].readB;
		currentError=dominatorWindows[0].errorRate(); currentAlignments=1;
		for (i=0; i<=lastDominatorWindow; i++) {
			if (dominatorWindows[i].readA!=currentReadA || dominatorWindows[i].readB!=currentReadB) {
				error=currentError/currentAlignments;
				minError[0][j]=Math.min(minError[0][j],error);
				if (dominatorWindows[i].readA!=currentReadA) {
					j++;
					currentReadA=dominatorWindows[i].readA;
					currentError=dominatorWindows[i].errorRate(); currentAlignments=1;
				}
				currentReadB=dominatorWindows[i].readB;
			}
			else {
				currentError+=dominatorWindows[i].errorRate();
				currentAlignments++;
			}
		}
		error=currentError/currentAlignments;
		minError[0][j]=Math.min(minError[0][j],error);
		// Computing best avg. error of every fragment to a consensus
		j=0;
		currentReadA=consensusWindows[0].readA; currentReadB=consensusWindows[0].readB;
		currentError=consensusWindows[0].errorRate(); currentAlignments=1;
		for (i=0; i<=lastConsensusWindow; i++) {
			if (consensusWindows[i].readA!=currentReadA || consensusWindows[i].readB!=currentReadB) {
				error=currentError/currentAlignments;
				p=consensusAlignmentsFormat?currentReadA:currentReadA-nReferences;
				while (j<=lastFrag2dom) {
					if (frags2dom[j]<p) j++;
					else if (frags2dom[j]>=p) break;
				}
				if (j<=lastFrag2dom && frags2dom[j]==p) minError[1][j]=Math.min(minError[1][j],error);
				if (consensusWindows[i].readA!=currentReadA) {
					j++;
					currentReadA=consensusWindows[i].readA;
					currentError=consensusWindows[i].errorRate(); currentAlignments=1;
				}
				currentReadB=consensusWindows[i].readB;
			}
			else {
				currentError+=consensusWindows[i].errorRate();
				currentAlignments++;
			}
		}
		error=currentError/currentAlignments;
		p=consensusAlignmentsFormat?currentReadA:currentReadA-nReferences;
		while (j<=lastFrag2dom) {
			if (frags2dom[j]<p) j++;
			else if (frags2dom[j]>=p) break;
		}
		if (j<=lastFrag2dom && frags2dom[j]==p) minError[1][j]=Math.min(minError[1][j],error);
		// Comparing
		n=0;
		for (i=0; i<=lastFrag2dom; i++) {
			if (minError[1][i]!=Math.POSITIVE_INFINITY && minError[1][i]>minError[0][i]) {
				if (verbose) System.err.println("ERROR: fragment "+frags2dom[i]+" (0-based from first fragment) aligned to its closest dominator with avg. error "+minError[0][i]+", but it aligns to its closest consensus with avg. error "+minError[1][i]);
				n++;
			}
		}
		out[5]=n;
		
		// Comparing surfaces
		if (surface[0].length<nFragments) surface = new long[2][nFragments];
		Math.set(surface[0],nFragments-1,0);
		Math.set(surface[1],nFragments-1,0);
		Arrays.sort(dominatorWindows,0,lastDominatorWindow+1);
		currentReadA=dominatorWindows[0].readA; currentSurface=0;
		currentStart=dominatorWindows[0].startA; currentEnd=dominatorWindows[0].endA;
		for (i=1; i<=lastDominatorWindow; i++) {
			if (dominatorWindows[i].readA!=currentReadA) {
				currentSurface+=currentEnd-currentStart+1;
				surface[0][dominatorAlignmentsFormat?currentReadA:currentReadA-nReferences]=currentSurface;
				currentReadA=dominatorWindows[i].readA; currentSurface=0;
				currentStart=dominatorWindows[i].startA; currentEnd=dominatorWindows[i].endA;
			}
			else if (dominatorWindows[i].startA>currentEnd) {
				currentSurface+=currentEnd-currentStart+1;
				currentStart=dominatorWindows[i].startA; currentEnd=dominatorWindows[i].endA;
			}
			else currentEnd=Math.max(currentEnd,dominatorWindows[i].endA);
		}
		currentSurface+=currentEnd-currentStart+1;
		surface[0][dominatorAlignmentsFormat?currentReadA:currentReadA-nReferences]=currentSurface;
		Arrays.sort(consensusWindows,0,lastConsensusWindow+1);
		currentReadA=consensusWindows[0].readA; currentSurface=0;
		currentStart=consensusWindows[0].startA; currentEnd=consensusWindows[0].endA;
		for (i=1; i<=lastConsensusWindow; i++) {
			if (consensusWindows[i].readA!=currentReadA) {
				currentSurface+=currentEnd-currentStart+1;
				surface[1][consensusAlignmentsFormat?currentReadA:currentReadA-nReferences]=currentSurface;
				currentReadA=consensusWindows[i].readA; currentSurface=0;
				currentStart=consensusWindows[i].startA; currentEnd=consensusWindows[i].endA;
			}
			else if (consensusWindows[i].startA>currentEnd) {
				currentSurface+=currentEnd-currentStart+1;
				currentStart=consensusWindows[i].startA; currentEnd=consensusWindows[i].endA;
			}
			else currentEnd=Math.max(currentEnd,consensusWindows[i].endA);
		}
		currentSurface+=currentEnd-currentStart+1;
		surface[1][consensusAlignmentsFormat?currentReadA:currentReadA-nReferences]=currentSurface;
		n=0;
		for (i=0; i<nFragments; i++) {
			if (surface[1][i]!=0 && surface[1][i]<surface[0][i]-SURFACE_THRESHOLD) {
				if (verbose) System.err.println("ERROR: fragment "+i+" has "+surface[0][i]+" bp covered by alignments to dominators, but just "+surface[1][i]+" bp covered by alignments to consensi.");
				n++;
			}
		}
		out[6]=n;
		sum=0;
		for (i=0; i<nFragments; i++) sum+=surface[0][i];
		out[7]=sum;
		sum=0;
		for (i=0; i<nFragments; i++) sum+=surface[1][i];
		out[8]=sum;
	}
	
	
	private static class FragConsAlignment implements Comparable {
		public boolean orientation;
		public int readA, readB, startA, endA, startB, endB, diffs;
		
		public FragConsAlignment() {
			this.readA=-1;
			this.readB=-1;
			this.orientation=false;
			this.startA=-1;
			this.endA=-1;
			this.startB=-1;
			this.endB=-1;
			this.diffs=0;
		}
		
		public FragConsAlignment(int ra, int rb, boolean or, int s, int e, int sb, int eb, int d) {
			set(ra,rb,or,s,e,sb,eb,d);
		}
		
		public void set(int ra, int rb, boolean or, int s, int e, int sb, int eb, int d) {
			readA=ra;
			readB=rb;
			orientation=or;
			startA=s;
			endA=e;
			startB=sb;
			endB=eb;
			diffs=d;
		}
		
		public double errorRate() {
			return ((double)(diffs<<1))/(endB-startB+endA-startA+2);
		}
		
		/**
		 * Sorts by readA,startA.
		 */
		public int compareTo(Object other) {
			FragConsAlignment otherAlignment = (FragConsAlignment)other;
			if (readA<otherAlignment.readA) return -1;
			else if (readA>otherAlignment.readA) return 1;
			if (startA<otherAlignment.startA) return -1;
			else if (startA>otherAlignment.startA) return 1;
			return 0;
		}
	}








	// ------------------------------- REPEAT TRACKS -------------------------------------

	/**
	 * Assume that we are given a set of blocks $B_x..B_y$, and a set of repeat tracks on 
	 * those blocks to be used to soft-mask alignments (let's call such repeat tracks 
	 * "input repeat tracks"). Assume also that we computed all (soft-masked) pairwise 
	 * alignments inside $B_x..B_y$, and that we ran our repeat inference up to the 
	 * consensi step. Finally, assume that we mapped the consensi back to $B_x..B_y$: 
	 * let's call "output repeat track" the merge of all alignments of such consensi to 
	 * the reads in $B_x..B_y$. The procedure compares the output repeat track to the 
	 * merge of all intervals in $B_x..B_y$ produced by $Factorize.java$, and to the merge
	 * of all intervals in $B_x..B_y$ that belong to a basin descriptor file.
	 *
	 * @param repeatTrackInput input repeat track ("null" discards it);
	 * @param repeatTrack*Format see $FilterAlignments.loadRepeatTrack()$.
	 */
	public static void checkConsensiRepeatTrack(String step1Dir, String repeatTrackInputPath, byte repeatTrackInputFormat, String repeatTrackOutputPath, byte repeatTrackOutputFormat, int nReads, String readLengthsFile, String readIDsFile, String qualityThresholdsFile, String qualitiesFile) throws IOException {
		final String FACTORIZE_OUTPUT_DIR = step1Dir+"/..";
		final String BASINS_DIR = step1Dir+"/"+IO.TAGS_DIR+"/"+IO.STEP4_DIR;
		final String BASINS_LIST = BASINS_DIR+"/"+IO.STEP5_DIR+"/"+IO.STEP5_DESCRIPTORS_FILE;
		final int BASIN_MINUS_FACTOR_THRESHOLD = 2000;  // Arbitrary
		final int BASIN_MINUS_REPEAT_TRACK_THRESHOLD = 2000;  // Arbitrary
		final int FACTOR_MINUS_QUALITY_THRESHOLD = 2000;  // Arbitrary
		
		int i, j;
		int last, sum;
		long repeatTrackInputSurface, repeatTrackOutputSurface, factorTrackSurface, basinTrackSurface, qualityTrackSurface;
		long factorBasinIntersection, factorQualityIntersection, basinOutputIntersection;
		long[][] intersections;  // (F,B,Ri,Ro,Q) x (F,B,Ri,Ro,Q)
		String str;
		BufferedReader br;
		int[] tmpArray, lastRepeatTrackInput, lastRepeatTrackOutput, lastFactorTrack, lastBasinTrack, lastQualityTrack;
		int[][] repeatTrackInput, repeatTrackOutput, factorTrack, basinTrack, qualityTrack;	
		
		Reads.nReads=nReads;
		Reads.maxReadLength=Reads.loadReadLengths(readLengthsFile);
		Reads.loadReadIDs(readIDsFile,nReads);
		Reads.loadQualities(qualityThresholdsFile,qualitiesFile);
		final int FIRST_READ = Reads.firstRead;
		final int LAST_READ = Reads.lastRead;
		tmpArray = new int[Reads.maxReadLength];
		IO.initialize();
		
		// Loading the complement of the input repeat track
		repeatTrackInput = new int[nReads][0];
		lastRepeatTrackInput = new int[nReads];
		if (repeatTrackInputPath.equalsIgnoreCase("null")) Math.set(lastRepeatTrackInput,nReads-1,-1);
		else FilterAlignments.loadRepeatTrack(repeatTrackInputPath,repeatTrackInputFormat,FIRST_READ,LAST_READ,repeatTrackInput,lastRepeatTrackInput);
		for (i=0; i<nReads; i++) {
			last=Intervals.complement(repeatTrackInput[i],lastRepeatTrackInput[i],Reads.getReadLength(FIRST_READ+i),tmpArray);
			if (last+1>repeatTrackInput[i].length) repeatTrackInput[i] = new int[last+1];
			System.arraycopy(tmpArray,0,repeatTrackInput[i],0,last+1);
			lastRepeatTrackInput[i]=last;
		}
		repeatTrackInputSurface=0;
		for (i=0; i<nReads; i++) {
			last=lastRepeatTrackInput[i];
			for (j=0; j<last; j+=2) repeatTrackInputSurface+=repeatTrackInput[i][j+1]-repeatTrackInput[i][j]+1;
		}
		
		// Loading the output repeat track
		repeatTrackOutput = new int[nReads][0];
		lastRepeatTrackOutput = new int[nReads];
		FilterAlignments.loadRepeatTrack(repeatTrackOutputPath,repeatTrackOutputFormat,FIRST_READ,LAST_READ,repeatTrackOutput,lastRepeatTrackOutput);
		repeatTrackOutputSurface=0;
		for (i=0; i<nReads; i++) {
			last=lastRepeatTrackOutput[i];
			for (j=0; j<last; j+=2) repeatTrackOutputSurface+=repeatTrackOutput[i][j+1]-repeatTrackOutput[i][j]+1;
		}
		
		// Loading factorization intervals
		factorTrack = new int[nReads][0];
		lastFactorTrack = new int[nReads];
		Math.set(lastFactorTrack,nReads-1,-1);
		loadFactorTrack(FACTORIZE_OUTPUT_DIR+"/intervals-alignments-all.txt",0,FIRST_READ,LAST_READ,factorTrack,lastFactorTrack,tmpArray);
		loadFactorTrack(FACTORIZE_OUTPUT_DIR+"/intervals-dense-all.txt",1,FIRST_READ,LAST_READ,factorTrack,lastFactorTrack,tmpArray);
		loadFactorTrack(FACTORIZE_OUTPUT_DIR+"/intervals-periodic-all.txt",2,FIRST_READ,LAST_READ,factorTrack,lastFactorTrack,tmpArray);
		compactTracks(nReads,FIRST_READ,factorTrack,lastFactorTrack,false);
		factorTrackSurface=0;
		for (i=0; i<nReads; i++) {
			last=lastFactorTrack[i];
			for (j=0; j<last; j+=2) factorTrackSurface+=factorTrack[i][j+1]-factorTrack[i][j]+1;
		}
		
		// Loading basin descriptors
		basinTrack = new int[nReads][0];
		lastBasinTrack = new int[nReads];
		Math.set(lastBasinTrack,nReads-1,-1);
		br = new BufferedReader(new FileReader(BASINS_LIST));
		str=br.readLine();
		while (str!=null) {
			loadFactorTrack(BASINS_DIR+"/"+str.trim(),3,FIRST_READ,LAST_READ,basinTrack,lastBasinTrack,tmpArray);
			str=br.readLine();
		}
		br.close();
		compactTracks(nReads,FIRST_READ,basinTrack,lastBasinTrack,false);
		basinTrackSurface=0;
		for (i=0; i<nReads; i++) {
			last=lastBasinTrack[i];
			for (j=0; j<last; j+=2) basinTrackSurface+=basinTrack[i][j+1]-basinTrack[i][j]+1;
		}
		
		// Loading quality
		qualityTrack = new int[nReads][0];
		lastQualityTrack = new int[nReads];
		Math.set(lastQualityTrack,nReads-1,-1);		
		for (i=0; i<nReads; i++) {
			last=Reads.getQualityAsTrack(FIRST_READ+i,Reads.MAX_HIGH_QUALITY_SCORE,tmpArray);
			if (last+1>qualityTrack[i].length) qualityTrack[i] = new int[last+1];
			System.arraycopy(tmpArray,0,qualityTrack[i],0,last+1);
			lastQualityTrack[i]=last;
		}
		qualityTrackSurface=0;
		for (i=0; i<nReads; i++) {
			last=lastQualityTrack[i];
			for (j=0; j<last; j+=2) qualityTrackSurface+=qualityTrack[i][j+1]-qualityTrack[i][j]+1;
		}
		
		// Computing intersections
		intersections = new long[5][5];  // (F,B,Ri,Ro,Q) x (F,B,Ri,Ro,Q)
		Math.set(intersections,0);
		for (i=0; i<nReads; i++) {
			factorBasinIntersection=Intervals.intersectionLength(factorTrack[i],lastFactorTrack[i],basinTrack[i],lastBasinTrack[i],1,true);
			intersections[0][1]+=factorBasinIntersection;
			intersections[0][2]+=Intervals.intersectionLength(factorTrack[i],lastFactorTrack[i],repeatTrackInput[i],lastRepeatTrackInput[i],1,true);
			intersections[0][3]+=Intervals.intersectionLength(factorTrack[i],lastFactorTrack[i],repeatTrackOutput[i],lastRepeatTrackOutput[i],1,true);
			factorQualityIntersection=Intervals.intersectionLength(factorTrack[i],lastFactorTrack[i],qualityTrack[i],lastQualityTrack[i],1,true);
			intersections[0][4]+=factorQualityIntersection;
			intersections[1][2]+=Intervals.intersectionLength(basinTrack[i],lastBasinTrack[i],repeatTrackInput[i],lastRepeatTrackInput[i],1,true);
			basinOutputIntersection=Intervals.intersectionLength(basinTrack[i],lastBasinTrack[i],repeatTrackOutput[i],lastRepeatTrackOutput[i],1,true);
			intersections[1][3]+=basinOutputIntersection;
			intersections[1][4]+=Intervals.intersectionLength(basinTrack[i],lastBasinTrack[i],qualityTrack[i],lastQualityTrack[i],1,true);
			intersections[2][3]+=Intervals.intersectionLength(repeatTrackInput[i],lastRepeatTrackInput[i],repeatTrackOutput[i],lastRepeatTrackOutput[i],1,true);
			intersections[2][4]+=Intervals.intersectionLength(repeatTrackInput[i],lastRepeatTrackInput[i],qualityTrack[i],lastQualityTrack[i],1,true);
			intersections[3][4]+=Intervals.intersectionLength(repeatTrackOutput[i],lastRepeatTrackOutput[i],qualityTrack[i],lastQualityTrack[i],1,true);
			// Consistency check: factors should be high-quality.
			sum=0;
			for (j=0; j<lastFactorTrack[i]; j+=2) sum+=factorTrack[i][j+1]-factorTrack[i][j]+1;
			if ( sum-factorQualityIntersection>=FACTOR_MINUS_QUALITY_THRESHOLD && 
				 Intervals.difference(factorTrack[i],lastFactorTrack[i],qualityTrack[i],lastQualityTrack[i],tmpArray,0,-1,FACTOR_MINUS_QUALITY_THRESHOLD)>=0
			   ) {
				System.err.println("ERROR: read "+(FIRST_READ+i)+" has "+(sum-factorQualityIntersection)+" bps in factor intervals that are not high-quality?!");
				System.err.print("factor track: ");
				for (int x=0; x<lastFactorTrack[i]; x+=2) System.err.print("("+factorTrack[i][x]+","+factorTrack[i][x+1]+") ");
				System.err.println();
				System.err.print("quality track: ");
				for (int x=0; x<lastQualityTrack[i]; x+=2) System.err.print("("+qualityTrack[i][x]+","+qualityTrack[i][x+1]+") ");
				System.err.println();
				System.exit(1);
			}
			// Consistency check: the basin track should be contained in the factors
			// track.
			sum=0;
			for (j=0; j<lastBasinTrack[i]; j+=2) sum+=basinTrack[i][j+1]-basinTrack[i][j]+1;
			if (sum-factorBasinIntersection>=BASIN_MINUS_FACTOR_THRESHOLD) {
				System.err.println("ERROR: read "+(FIRST_READ+i)+" has "+(sum-factorBasinIntersection)+" bps in basin intervals that are not factor intervals?!");
				System.err.print("basin track: ");
				for (int x=0; x<lastBasinTrack[i]; x+=2) System.err.print("("+basinTrack[i][x]+","+basinTrack[i][x+1]+") ");
				System.err.println();
				System.exit(1);
			}
			// Consistency check: the basin track should be contained in the output repeat
			// track.
			if (sum-basinOutputIntersection>=BASIN_MINUS_REPEAT_TRACK_THRESHOLD) {
				System.err.println("WARNING: read "+(FIRST_READ+i)+" has "+(sum-basinOutputIntersection)+" bps in basin intervals that are not in the output repeat track?!");
				System.err.print("basin track: ");
				for (int x=0; x<lastBasinTrack[i]; x+=2) System.err.print("("+basinTrack[i][x]+","+basinTrack[i][x+1]+") ");
				System.err.println();
				System.err.print("output repeat track: ");
				for (int x=0; x<lastRepeatTrackOutput[i]; x+=2) System.err.print("("+repeatTrackOutput[i][x]+","+repeatTrackOutput[i][x+1]+") ");
				System.err.println();
			}
		}
		
		// Outputting
		System.err.println("Tracks:");
		System.err.println("F = union of all factorization intervals");
		System.err.println("B = union of all basin intervals");
		System.err.println("Ri = input repeat track");
		System.err.println("Ro = output repeat track");
		System.err.println("Q = high quality track");
		System.err.println();
		System.err.println("T1, T2, |T1 \\setminus T2|/|T1|, |T1 \\intersect T2|/|T1 \\union T2|, |T2 \\setminus T1|/|T2|");
		System.err.println();
		System.err.println("    F,  B, "+IO.getPercent(factorTrackSurface-intersections[0][1],factorTrackSurface)+"%, "+IO.getPercent(intersections[0][1],factorTrackSurface+basinTrackSurface-intersections[0][1])+"%, "+IO.getPercent(basinTrackSurface-intersections[0][1],basinTrackSurface)+"%");
		System.err.println("*   F,-Ri, "+IO.getPercent(factorTrackSurface-intersections[0][2],factorTrackSurface)+"%, "+IO.getPercent(intersections[0][2],factorTrackSurface+repeatTrackInputSurface-intersections[0][2])+"%, "+IO.getPercent(repeatTrackInputSurface-intersections[0][2],repeatTrackInputSurface)+"%");
		System.err.println("*   F, Ro, "+IO.getPercent(factorTrackSurface-intersections[0][3],factorTrackSurface)+"%, "+IO.getPercent(intersections[0][3],factorTrackSurface+repeatTrackOutputSurface-intersections[0][3])+"%, "+IO.getPercent(repeatTrackOutputSurface-intersections[0][3],repeatTrackOutputSurface)+"%");
		System.err.println("    F,  Q, "+IO.getPercent(factorTrackSurface-intersections[0][4],factorTrackSurface)+"%, "+IO.getPercent(intersections[0][4],factorTrackSurface+qualityTrackSurface-intersections[0][4])+"%, "+IO.getPercent(qualityTrackSurface-intersections[0][4],qualityTrackSurface)+"%");
		System.err.println("*   B,-Ri, "+IO.getPercent(basinTrackSurface-intersections[1][2],basinTrackSurface)+"%, "+IO.getPercent(intersections[1][2],basinTrackSurface+repeatTrackInputSurface-intersections[1][2])+"%, "+IO.getPercent(repeatTrackInputSurface-intersections[1][2],repeatTrackInputSurface)+"%");
		System.err.println("*   B, Ro, "+IO.getPercent(basinTrackSurface-intersections[1][3],basinTrackSurface)+"%, "+IO.getPercent(intersections[1][3],basinTrackSurface+repeatTrackOutputSurface-intersections[1][3])+"%, "+IO.getPercent(repeatTrackOutputSurface-intersections[1][3],repeatTrackOutputSurface)+"%");
		System.err.println("    B,  Q, "+IO.getPercent(basinTrackSurface-intersections[1][4],basinTrackSurface)+"%, "+IO.getPercent(intersections[1][4],basinTrackSurface+qualityTrackSurface-intersections[1][4])+"%, "+IO.getPercent(qualityTrackSurface-intersections[1][4],qualityTrackSurface)+"%");
		System.err.println("* -Ri, Ro, "+IO.getPercent(repeatTrackInputSurface-intersections[2][3],repeatTrackInputSurface)+"%, "+IO.getPercent(intersections[2][3],repeatTrackInputSurface+repeatTrackOutputSurface-intersections[2][3])+"%, "+IO.getPercent(repeatTrackOutputSurface-intersections[2][3],repeatTrackOutputSurface)+"%");
		System.err.println("  -Ri,  Q, "+IO.getPercent(repeatTrackInputSurface-intersections[2][4],repeatTrackInputSurface)+"%, "+IO.getPercent(intersections[2][4],repeatTrackInputSurface+qualityTrackSurface-intersections[2][4])+"%, "+IO.getPercent(qualityTrackSurface-intersections[2][4],qualityTrackSurface)+"%");
		System.err.println("   Ro,  Q, "+IO.getPercent(repeatTrackOutputSurface-intersections[3][4],repeatTrackOutputSurface)+"%, "+IO.getPercent(intersections[3][4],repeatTrackOutputSurface+qualityTrackSurface-intersections[3][4])+"%, "+IO.getPercent(qualityTrackSurface-intersections[3][4],qualityTrackSurface)+"%");
	}
	
	
	/**
	 * Appends the contents of $file$ to $factorTrack,lastFactorTrack$.
	 *
	 * @param file list of factorization intervals;
	 * @param fileType 0=alignment intervals; 1=dense; 2=periodic; 3=basin descriptor;
	 * @param tokens temporary space, assumed to be large enough.
	 */
	private static final void loadFactorTrack(String file, int fileType, int firstRead, int lastRead, int[][] factorTrack, int[] lastFactorTrack, int[] tokens) throws IOException {
		int i;
		int last, length, currentRead, currentStart, currentEnd, read, start, end;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(file));
		if (fileType==3) str=br.readLine();  // Skipping header
		str=br.readLine();
		currentRead=-1; currentStart=-1; currentEnd=-1;
		while (str!=null) {
			switch (fileType) {
				case 0: i=str.indexOf(","); read=Integer.parseInt(str.substring(0,i))-firstRead; AlignmentInterval.readAlignmentsFile(str,i+1,tokens,0); start=tokens[1]; end=tokens[2]; break;
				case 1: i=str.indexOf(","); read=Integer.parseInt(str.substring(0,i))-firstRead; DenseSubstring.readDenseSubstringsFile(str,i+1,tokens,0); start=tokens[2]; end=tokens[3]; break;
				case 2: i=str.indexOf(","); read=Integer.parseInt(str.substring(0,i))-firstRead; PeriodicSubstringInterval.readPeriodicSubstringsFile(str,i+1,tokens,0); start=tokens[1]; end=tokens[2]; break;
				case 3: IntervalGraphStep3.readBasinDescriptorRecord(str,tokens); read=tokens[0]-firstRead; start=tokens[1]; end=tokens[2]; break;
				default: read=-1; start=-1; end=-1;
			} 
			if (read!=currentRead) {
				if (currentRead!=-1) {
					last=lastFactorTrack[currentRead];
					length=factorTrack[currentRead].length;
					if (last+2>=length) {
						int[] newTrack = new int[Math.max(length<<1,last+3)];
						System.arraycopy(factorTrack[currentRead],0,newTrack,0,length);
						factorTrack[currentRead]=newTrack;
					}
					factorTrack[currentRead][last+1]=currentStart;
					factorTrack[currentRead][last+2]=currentEnd;
					lastFactorTrack[currentRead]+=2;
				}
				currentRead=read; currentStart=start; currentEnd=end;
			}
			else if (start>currentEnd) {
				last=lastFactorTrack[currentRead];
				length=factorTrack[currentRead].length;
				if (last+2>=length) {
					int[] newTrack = new int[Math.max(length<<1,last+3)];
					System.arraycopy(factorTrack[currentRead],0,newTrack,0,length);
					factorTrack[currentRead]=newTrack;
				}
				factorTrack[currentRead][last+1]=currentStart;
				factorTrack[currentRead][last+2]=currentEnd;
				lastFactorTrack[currentRead]+=2;
				currentStart=start; currentEnd=end;
			}
			else currentEnd=Math.max(currentEnd,end);
			str=br.readLine();
		}
		br.close();
		if (currentRead!=-1) {
			last=lastFactorTrack[currentRead];
			length=factorTrack[currentRead].length;
			if (last+2>=length) {
				int[] newTrack = new int[Math.max(length<<1,last+3)];
				System.arraycopy(factorTrack[currentRead],0,newTrack,0,length);
				factorTrack[currentRead]=newTrack;
			}
			factorTrack[currentRead][last+1]=currentStart;
			factorTrack[currentRead][last+2]=currentEnd;
			lastFactorTrack[currentRead]+=2;
		}	
	}
	
	
	/**
	 * Merges all overlapping intervals in $factorTrack$.
	 */
	public static final void compactTracks(int nReads, int firstRead, int[][] factorTrack, int[] lastFactorTrack, boolean verbose) {
		int i, j, k, h;
		int max, last, currentStart, currentEnd;
		FragmentsStep1.Pair[] pairs;
		
		// Allocating shared data structures
		max=0;
		for (i=0; i<nReads; i++) max=Math.max(max,lastFactorTrack[i]+1);
		max>>=1;
		pairs = new FragmentsStep1.Pair[max];
		for (i=0; i<max; i++) pairs[i] = new FragmentsStep1.Pair();
		
		// Compacting
		for (i=0; i<nReads; i++) {
			last=lastFactorTrack[i];
			if (last<=1) continue;
			k=-1;
			for (j=0; j<last; j+=2) {
				k++;
				pairs[k].start=factorTrack[i][j];
				pairs[k].end=factorTrack[i][j+1];
			}
			Arrays.sort(pairs,0,k+1);			
			if (verbose) {
				System.err.print("SORTED FACTOR TRACK OF READ "+i+"-th = "+(firstRead+i)+": ");			
				for (int x=0; x<=k; x++) System.err.print("("+pairs[x].start+","+pairs[x].end+")");
				System.err.println();
			}
			h=-1; currentStart=pairs[0].start; currentEnd=pairs[0].end;
			for (j=1; j<=k; j++) {
				if (pairs[j].start>currentEnd) {
					h++;
					pairs[h].start=currentStart; pairs[h].end=currentEnd;
					currentStart=pairs[j].start; currentEnd=pairs[j].end;
				}
				else currentEnd=Math.max(currentEnd,pairs[j].end);
			}
			h++;
			pairs[h].start=currentStart; pairs[h].end=currentEnd;
			k=-1;
			for (j=0; j<=h; j++) {
				factorTrack[i][++k]=pairs[j].start;
				factorTrack[i][++k]=pairs[j].end;
			}
			lastFactorTrack[i]=k;
			if (verbose) {
				System.err.print("COMPACTED FACTOR TRACK OF READ "+i+"-th = "+(firstRead+i)+": ");				
				for (int x=0; x<lastFactorTrack[i]; x+=2) System.err.print("("+factorTrack[i][x]+","+factorTrack[i][x+1]+")");
				System.err.println();
			}
		}
	}


}