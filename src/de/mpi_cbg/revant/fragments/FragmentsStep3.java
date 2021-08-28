package de.mpi_cbg.revant.fragments;

import java.util.Arrays;
import java.util.HashSet;
import java.io.*;
import java.awt.*;
import java.awt.image.*;
import javax.imageio.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Stream;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Colors;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Intervals;
import de.mpi_cbg.revant.intervalgraph.IntervalGraph;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;


/**
 * Checks whether the sets of reference strings and fragments built by $FragmentsStep2$   
 * are consistent, and prints statistics for each module (counts of fragment types,  
 * coverage histogram of each reference). Then, selects a high-similarity set of  
 * non-conflicting alignments for every (fragment,reference) pair: this is useful, since
 * a fragment might have many alignments with the reference, and we want every character 
 * of a fragment to vote at most once for the consensus, and every character of the 
 * reference to have at most one vote from a specific fragment.
 *
 * Remark: the fragment-reference alignments used at this stage might be very different
 * from those used during repeat inference. E.g. now we might be using a smaller  
 * minAlignmentLength, or a bigger error rate; moreover, the aligner is now working on a
 * much smaller problem, so its heuristics might behave differently.
 */
public class FragmentsStep3 {
	/**
	 * Encoding of a conflict graph
	 */
	public static IntervalGraph.Edge[] edgePool;  // Pool of reused conflict edges
	private static int lastInPool;
	private static int[][] neighbors;
	private static int[] lastNeighbor;
	
	/**
	 * Number of conflict graphs with a given max claw size.
	 */
	private static int[] clawHistogram, clawHistogramPeriodic;
	
	/**
	 * Temporary space
	 */
	public static FragmentsStep1.FFAlignment[] alignments;
	private static boolean[] isNeighbor, currentSolution, tmpBoolean;
	private static int[] tmpArray1, tmpArray2, tmpArray3;
	public static int[] stack;
	private static long[] tmpArrayLong = new long[10];
	private static Stream[] cliques;

	
	/**
	 * Remark: in every LAshow file, readA is a fragment ID and readB is a reference ID.
	 *
	 * @param args 
	 * 0: Step1 directory;
	 * 1: min length of an alignment to use for computing the coverage plot of a new 
	 *    reference (for debugging, it should be the same as in repeat detection);
	 * 2: 1=filter alignments; 0=compute statistics on conflict graphs (slow);
	 * 3: 1=when filtering alignments, prints conflict graphs as rectangles (slow);
	 * 4: format of LAshow files: 1=fragments and references are in different DBs, and the
	 *    IDs of both start from 1; 0=references are in block 1 and fragments are in block
	 *    2 of the same DB, so fragment IDs start from $nReferences$.
	 */
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		final int MIN_ALIGNMENT_LENGTH_FOR_COVERAGE = Integer.parseInt(args[1]);
		final boolean FILTER_MODE = Integer.parseInt(args[2])==1;
		final boolean PRINT_CONFLICT_GRAPHS = Integer.parseInt(args[3])==1;
		final boolean LASHOW_FORMAT = Integer.parseInt(args[4])==1;
		
		final String STEP4_DIR = STEP1_DIR+"/"+IO.TAGS_DIR+"/"+IO.STEP4_DIR;
		final String STEP5_DIR = STEP4_DIR+"/"+IO.STEP5_DIR;
		final String FRAGMENTS_COORDINATES_DIR = STEP5_DIR+"/"+IO.FRAGMENTS_LABEL;
		final String FRAGMENTS_ALIGNMENTS_DIR = STEP5_DIR+"/fragments-strings-alignments";
		final String FRAGMENTS_ALIGNMENTS_DIR_NEW = FRAGMENTS_ALIGNMENTS_DIR+"/fragments-strings-alignments-new";
		final String STATS_DIR = FRAGMENTS_ALIGNMENTS_DIR_NEW+"/stats";
		final String ERROR_RATE_HISTOGRAM_FILE = STATS_DIR+"/errorRateHistogram.txt";
		final String FRAGMENTS_LENGTHS_DIR = FRAGMENTS_ALIGNMENTS_DIR+"/fragments-strings-new";
		final String REFERENCE_LENGTHS_DIR = FRAGMENTS_ALIGNMENTS_DIR+"/references-strings-new";

		final int CAPACITY = 100;  // Arbitrary
		final int POOL_CAPACITY = 1000;  // Arbitrary
		final int LENGTH_THRESHOLD = IO.quantum<<2;  // Arbitrary
		final double MIN_ERROR = 0.0;
		final double MAX_ERROR = 0.4;
		final int N_HISTOGRAM_BINS = 100;  // Arbitrary
		final double HISTOGRAM_QUANTUM = (MAX_ERROR-MIN_ERROR)/N_HISTOGRAM_BINS;
		final int MAX_CLAW_SIZE = 100;   // Arbitrary
		final int MIN_CONFLICT_INTERSECTION = IO.quantum<<1;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		
		int i, j, k;
		int length, maxReferenceLength, nTooLong, currentRead, lastTmpWindow;
		int nFragments, nReferences, totalNFragments, totalNReferences, clawSize, maxClawSize;
		int PATH_ID, BASIN_TYPE, PERIOD;
		double errorRate;
		String str, id;
		String BASINS_FILE, FRAGMENTS_LENGTHS_FILE, REFERENCE_LENGTHS_FILE, CONNECTION_FILE, FR_ALIGNMENTS_FILE;
		File file;
		BufferedReader br;
		BufferedWriter bw;
		int[] fragmentStats = new int[7];
		int[] totalFragmentStats = new int[7];
		int[] referenceLengths = new int[CAPACITY];
		int[] errorRateHistogram = new int[N_HISTOGRAM_BINS];
		int[] tmpArray = new int[10];
		String[] list;
		FragmentsStep1.PermutationWindow[] tmpWindows;
		int[][] referenceCoverage = new int[CAPACITY][0];
		
		alignments = new FragmentsStep1.FFAlignment[CAPACITY];
		for (i=0; i<alignments.length; i++) alignments[i] = new FragmentsStep1.FFAlignment();
		clawHistogram = new int[MAX_CLAW_SIZE+1];
		Math.set(clawHistogram,MAX_CLAW_SIZE,0);
		clawHistogramPeriodic = new int[MAX_CLAW_SIZE+1];
		Math.set(clawHistogramPeriodic,MAX_CLAW_SIZE,0);
		IO.initialize();
		totalNFragments=0; totalNReferences=0; maxClawSize=0;
		Math.set(totalFragmentStats,totalFragmentStats.length-1,0);
		Math.set(errorRateHistogram,N_HISTOGRAM_BINS-1,0);
		tmpWindows = new FragmentsStep1.PermutationWindow[CAPACITY];
		for (i=0; i<CAPACITY; i++) tmpWindows[i] = new FragmentsStep1.PermutationWindow();
		edgePool = new IntervalGraph.Edge[POOL_CAPACITY];
		for (i=0; i<edgePool.length; i++) edgePool[i] = new IntervalGraph.Edge();
		stack = new int[POOL_CAPACITY];
		file = new File(FRAGMENTS_ALIGNMENTS_DIR_NEW);
		list=file.list();
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			
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
			
			// Loading reference lengths
			REFERENCE_LENGTHS_FILE=REFERENCE_LENGTHS_DIR+"/reference-"+id+"-lengths.txt";
			br = new BufferedReader(new FileReader(REFERENCE_LENGTHS_FILE));
			nReferences=0; maxReferenceLength=0; str=br.readLine();
			while (str!=null) {
				length=Integer.parseInt(str);
				nReferences++;
				if (nReferences==referenceLengths.length) {
					int[] newLengths = new int[referenceLengths.length<<1];
					System.arraycopy(referenceLengths,0,newLengths,0,referenceLengths.length);
					referenceLengths=newLengths;
				}
				referenceLengths[nReferences]=length;
				maxReferenceLength=Math.max(maxReferenceLength,length);
				str=br.readLine();
			}
			br.close();
			if (nReferences>referenceCoverage.length) {
				int[][] newCoverage = new int[nReferences][0];
				for (j=0; j<referenceCoverage.length; j++) {
					newCoverage[j]=referenceCoverage[j];
					referenceCoverage[j]=null;
				}
				referenceCoverage=newCoverage;
			}
			for (j=0; j<nReferences; j++) {
				if (referenceCoverage[j].length<referenceLengths[j]) referenceCoverage[j] = new int[referenceLengths[j]];
				Math.set(referenceCoverage[j],referenceLengths[j]-1,0);
			}
			totalNReferences+=nReferences;
			System.err.println("done");
			
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
			
			// Checking fragments-references alignments
			Math.set(fragmentStats,fragmentStats.length-1,0);
			fragmentStats[6]=nTooLong;
			FR_ALIGNMENTS_FILE=FRAGMENTS_ALIGNMENTS_DIR_NEW+"/"+list[i]+"/"+IO.FRAGMENTS_REFERENCE_FILE+".txt";
			br = new BufferedReader(new FileReader(FR_ALIGNMENTS_FILE),IO.BUFFER_SIZE);
			currentRead=-1; lastTmpWindow=-1;
			str=br.readLine(); str=br.readLine();  // Skipping header
			str=br.readLine();
			while (str!=null) {
				Alignments.readAlignmentFile(str);
				if (!LASHOW_FORMAT) Alignments.readA-=nReferences;
				errorRate=((double)(Alignments.diffs<<1))/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2);
				if (errorRate<MIN_ERROR) errorRateHistogram[0]++;
				else if (errorRate>=MAX_ERROR) errorRateHistogram[N_HISTOGRAM_BINS-1]++;
				else errorRateHistogram[(int)((errorRate-MIN_ERROR)/HISTOGRAM_QUANTUM)]++;
				if (Alignments.readA-1!=currentRead) {
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
					if (currentRead!=-1) alignmentStats(currentRead,tmpWindows,lastTmpWindow,fragmentStats,referenceLengths,referenceCoverage,MIN_ALIGNMENT_LENGTH_FOR_COVERAGE);
					currentRead=Alignments.readA-1;
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
			if (currentRead!=-1) alignmentStats(currentRead,tmpWindows,lastTmpWindow,fragmentStats,referenceLengths,referenceCoverage,MIN_ALIGNMENT_LENGTH_FOR_COVERAGE);
			for (j=0; j<fragmentStats.length; j++) totalFragmentStats[j]+=fragmentStats[j];
			
			// Printing statistics
			bw = new BufferedWriter(new FileWriter(STATS_DIR+"/fragmentStats-"+id+".txt"));
			for (j=0; j<fragmentStats.length; j++) bw.write(fragmentStats[j]+",");
			bw.close();
			for (j=0; j<nReferences; j++) {
				bw = new BufferedWriter(new FileWriter(STATS_DIR+"/referenceCoverage-"+id+"-"+j+".txt"));
				for (k=0; k<referenceLengths[j]; k++) bw.write(referenceCoverage[j][k]+"\n");
				bw.close();
			}
			
			// Building conflict graphs
			clawSize=filterAlignments(BASIN_TYPE,FR_ALIGNMENTS_FILE,MIN_CONFLICT_INTERSECTION,referenceLengths,FILTER_MODE,BASIN_TYPE==Constants.INTERVAL_PERIODIC?clawHistogramPeriodic:clawHistogram,PRINT_CONFLICT_GRAPHS);
			maxClawSize=Math.max(maxClawSize,clawSize);
		}
		
		// Symmetrizing bitvectors
		buildSymmetricBitvectors(FRAGMENTS_ALIGNMENTS_DIR_NEW);
		
		// Displaying global statistics
		System.out.println();
		System.out.println("Overall statistics:  ("+totalNFragments+" fragments, "+totalNReferences+" references)");
		System.out.println("Too long fragments: "+totalFragmentStats[6]+" ("+IO.format(100*((double)totalFragmentStats[6])/totalNFragments)+"%)");
		System.out.println("Fragments with holes: "+totalFragmentStats[5]+" ("+IO.format(100*((double)totalFragmentStats[5])/totalNFragments)+"%)");
		System.out.println("Type 4 fragments: "+totalFragmentStats[4]+" ("+IO.format(100*((double)totalFragmentStats[4])/totalNFragments)+"%)");
		System.out.println("Type 3 fragments: "+totalFragmentStats[3]+" ("+IO.format(100*((double)totalFragmentStats[3])/totalNFragments)+"%)");
		System.out.println("Type 2 fragments: "+totalFragmentStats[2]+" ("+IO.format(100*((double)totalFragmentStats[2])/totalNFragments)+"%)");
		System.out.println("Type 1 fragments: "+totalFragmentStats[1]+" ("+IO.format(100*((double)totalFragmentStats[1])/totalNFragments)+"%)");
		System.out.println("Type 0 fragments: "+totalFragmentStats[0]+" ("+IO.format(100*((double)totalFragmentStats[0])/totalNFragments)+"%)");
		
		// Printing histograms
		bw = new BufferedWriter(new FileWriter(ERROR_RATE_HISTOGRAM_FILE));
		for (i=0; i<N_HISTOGRAM_BINS; i++) bw.write((MIN_ERROR+i*HISTOGRAM_QUANTUM)+","+errorRateHistogram[i]+"\n");
		bw.close();
		
		
		if (!FILTER_MODE) {
			System.err.println("clawHistogram:");
			for (i=0; i<=MAX_CLAW_SIZE; i++) System.err.println(i+","+clawHistogram[i]);
			System.err.println("clawHistogramPeriodic:");
			for (i=0; i<=MAX_CLAW_SIZE; i++) System.err.println(i+","+clawHistogramPeriodic[i]);
		}
		
	}


	/**
	 * Increments $fragmentStats$ and $referenceCoverage$.
	 *
	 * Remark: ideally, the coverage of a reference by its fragments of length at least 
	 * equal to the minAlignmentLength used during repeat detection, should be uniform 
	 * (except for boundary effects), the intervals contained in a reference during 
	 * factorization should have been modeled as separate repeats, and their fragments 
	 * should not have been assigned to this reference. However, in practice there can be 
	 * coverage peaks, for example because:
	 *
	 * 1. a contained interval might be assigned to the same kernel as its container;
	 * 1.1 e.g. a dense substring might be contained in another dense substring, and the 
	 * two are not forcibly separated by $IntervalGraphStep2.containment2insertion()$;
	 * even alignment intervals might not be separated by such procedure, if maximality
	 * and orientation are not compatible;
	 * 1.2 the kernel graph built in Step 3 can merge distinct path kernels, e.g. because
	 * the contained does not pass a frequency threshold; all its occurrences are then 
	 * marked with the tags of its containers (contrary to Step 5, which does not reassign
	 * the intervals of redundant repeats); see $IntervalGraphStep3.{markRepeats(),
	 * getKernelFrequencyThreshold()}$.
	 * 2. the alignments assigned by factorization to a dense substring of substring type 
	 * (or to a set of dense substrings that are eventually assigned the same kernel tag) 
	 * might already have peaks (this happens very frequently). This might happen because 
	 * we did not factorize correctly, and because peaks computed just with the alignments 
	 * assigned to the dense substring might not be visible when considering all 
	 * alignments during factorization.
	 * 3. what looks like an alignment interval in the factorization of one read, might 
	 * map to completely different types of intervals in other reads, and thereby its
	 * kernel might be discarded e.g. because of not enough frequency in Step 3;
	 * 4. there can be regions of low coverage, for example because there could be small
	 * substrings of not very high quality in the reference.
	 *
	 * Weaker reasons for not comparing the shape of the coverage of the reference that 
	 * we see here, to the shape of the coverage computed with the alignments used during 
	 * factorization of the same read, are that the two sets of alignments use different 
	 * settings (e.g. different minAlignmentLength) and thus the output of the aligner 
	 * might be different, and that the fragments we use here can be much shorter than
	 * the original intervals.
	 *
	 * @param fragmentStats number of fragments that: 
	 * 0: have a single containment alignment to at least one of the references;
	 * 1: they are not (0), but they are fully covered by alignments to one specific 
	 * reference, in one specific orientation;
	 * 2: they are not (1), but they are fully covered by alignments to one or more 
	 * references, in any orientation; 
	 * 3: like (2), but limited to a prefix or suffix;
	 * 4: like (2), but limited to a substring that is not a prefix/suffix;
	 * 5: none of the above;
	 * @param referenceCoverage number of alignments that cover each position of each 
	 * reference;
	 * @param minAlignmentLength min length of an alignment used for $referenceCoverage$;
	 * to detect errors, it should be the same as in the repeat detection pipeline.
	 */
	public static final void alignmentStats(int fragmentID, FragmentsStep1.PermutationWindow[] windows, int lastWindow, int[] fragmentStats, int[] referenceLengths, int[][] referenceCoverage, int minAlignmentLength) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int FRAGMENT_LENGTH_THRESHOLD = Reads.getReadLength(fragmentID)-IDENTITY_THRESHOLD;
		boolean orientation, currentOrientation, found, hasHole;
		boolean type0, type1, type2, type3, type4;
		int i, j;
		int readB, startA, endA, currentReadB, currentStartA, currentEndA;
		int from, to, reference;
		
		// Type 0
		type0=false;
		for (i=0; i<=lastWindow; i++) {
			if (windows[i].startA<=IDENTITY_THRESHOLD && windows[i].endA>=FRAGMENT_LENGTH_THRESHOLD) {
				type0=true;
				break;
			}
		}
		
		// Type 1
		if (type0) type1=true;
		else {
			type1=false;
			FragmentsStep1.PermutationWindow.order=FragmentsStep1.PermutationWindow.ORDER_READB_ORIENTATION_STARTA_STARTB;
			if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
			currentReadB=-1; currentOrientation=false; currentStartA=-1; currentEndA=-1;
			found=true;
			for (i=0; i<=lastWindow; i++) {
				readB=windows[i].readB; orientation=windows[i].orientation;
				startA=windows[i].startA; endA=windows[i].endA;
				if (readB!=currentReadB || orientation!=currentOrientation) {
					if (currentReadB!=-1) {
						if (currentEndA<FRAGMENT_LENGTH_THRESHOLD) found=false;
						if (found) {
							type1=true;
							break;
						}
					}
					currentReadB=readB; currentOrientation=orientation;
					currentStartA=startA; currentEndA=endA;
					found=currentStartA<=IDENTITY_THRESHOLD;
					continue;
				}
				if (startA>currentEndA+IDENTITY_THRESHOLD) {
					found=false;
					currentStartA=startA; currentEndA=endA;
				}
				else currentEndA=Math.max(currentEndA,endA);
			}
			if (currentReadB!=-1) {
				if (currentEndA<FRAGMENT_LENGTH_THRESHOLD) found=false;
				if (found) type1=true;
			}
		}	
		
		// Type 2
		if (type1) type2=true;
		else {
			type2=false;
			FragmentsStep1.PermutationWindow.order=FragmentsStep1.PermutationWindow.ORDER_STARTA;
			if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
			currentStartA=windows[0].startA; 
			if (currentStartA<=IDENTITY_THRESHOLD) {
				currentEndA=windows[0].endA;
				for (i=1; i<=lastWindow; i++) {
					startA=windows[i].startA; endA=windows[i].endA;
					if (startA>currentEndA+IDENTITY_THRESHOLD) {
						type2=false;
						break;
					}
					else currentEndA=Math.max(currentEndA,endA);
				}
				if (currentEndA<FRAGMENT_LENGTH_THRESHOLD) type2=false;
			}
		}
		
		// Type 3,4 = prefix/suffix, substring fully covered
		if (type2) { type3=false; type4=false; }
		else {
			if (FragmentsStep1.PermutationWindow.order!=FragmentsStep1.PermutationWindow.ORDER_STARTA) {
				FragmentsStep1.PermutationWindow.order=FragmentsStep1.PermutationWindow.ORDER_STARTA;
				if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
			}
			hasHole=false;
			currentStartA=windows[0].startA; currentEndA=windows[0].endA;
			for (i=1; i<=lastWindow; i++) {
				startA=windows[i].startA; endA=windows[i].endA;
				if (startA>currentEndA+IDENTITY_THRESHOLD) {
					hasHole=true;
					break;
				}
				else currentEndA=Math.max(currentEndA,endA);
			}
			if (hasHole) { type3=false; type4=false; }
			else {
				if (currentStartA<=IDENTITY_THRESHOLD || currentEndA>=FRAGMENT_LENGTH_THRESHOLD) { type3=true; type4=false; }
				else { type3=false; type4=true; }
			}
		}
		
		if (type0) fragmentStats[0]++;
		else if (type1) fragmentStats[1]++;
		else if (type2) fragmentStats[2]++;
		else if (type3) fragmentStats[3]++;
		else if (type4) fragmentStats[4]++;
		else fragmentStats[5]++;
		
		// References
		FragmentsStep1.PermutationWindow.order=FragmentsStep1.PermutationWindow.ORDER_READB_STARTB;
		if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
		for (i=0; i<=lastWindow; i++) {
			reference=windows[i].readB;
			from=Math.max(windows[i].startB,0);
			to=Math.min(windows[i].endB,referenceLengths[reference]-1);
			if (to-from+1>=minAlignmentLength) {
				for (j=from; j<=to; j++) referenceCoverage[reference][j]++;
			}
		}
	}
	
	
	
	
	
	
	
	
	// ------------------------------ FILTERING ALIGNMENTS -------------------------------
	
	/**
	 * For every fragment-reference pair, keeps only a subset of alignments such that:
	 * (1) every pair of alignments share $< minIntersection$ bps in the fragment and in 
	 * the reference; (2) the set of alignments maximizes the total number of matching
	 * characters in the fragment. This is done by computing a max-weight independent set 
	 * in a "conflict graph" in which every alignment is a node, and two alignments are 
	 * adjacent iff they share $>=minIntersection$ bps in some sequence. 
	 *
	 * This problem was first defined and proven to be NP-hard in 
	 * \cite{bafna1996nonoverlapping}. In that paper, however, they assume 
	 * $minIntersection=1$ and that no two alignments are contained in one another in 
	 * either read.
	 *
	 * Remark: this allows a fragment to be a permutation of a subsequence (i.e. of a 
	 * non-adjacent set of substrings) of the reference. This is useful in microsats, a
	 * fragment of which might contain internal variants of the period that have not been
	 * assigned to distinct repeats. However, this might be too permissive for non-
	 * periodic repeats.
	 *
	 * @param alignmentsFile readA=fragment ID; readB=reference ID;
	 * @param mode TRUE=filters alignments; FALSE=does not filter alignments, and instead
	 * collects statistics on max claw size (see $getMaxClawSize()$ for definitions);
	 * @param histogram if $mode=FALSE$, the procedure adds to $histogram$ the number of 
	 * talons of each size (this array is just incremented, never reset to zero);
	 * @param printGraphs TRUE: draws a figure for every conflict graph, representing it
	 * as a set of rectangles;
	 * @return if $mode=TRUE$, the total number of alignments kept; otherwise, the max 
	 * number of talons in a claw of the conflict graph.
	 */
	public static final int filterAlignments(int basinType, String alignmentsFile, int minIntersection, int[] referenceLengths, boolean mode, int[] histogram, boolean printGraphs) throws IOException {
		final int p = alignmentsFile.lastIndexOf(".");
		int i, j;
		int currentReadA, currentReadB, maxNeighbors, maxClawSize, clawSize, totalSize;
		int firstAlignmentID;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(alignmentsFile.substring(0,p)+"-keep.txt"));
		br = new BufferedReader(new FileReader(alignmentsFile));
		br.readLine(); br.readLine();  // Skipping header
		str=br.readLine();
		i=-1; currentReadA=-1; currentReadB=-1; maxClawSize=0; totalSize=0; firstAlignmentID=0;
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			if (Alignments.readA!=currentReadA || Alignments.readB!=currentReadB) {
				if (currentReadA!=-1 && currentReadB!=-1) {
					maxNeighbors=buildConflictGraph(basinType,alignments,i+1,minIntersection);
					if (tmpArray1==null || tmpArray1.length<maxNeighbors) tmpArray1 = new int[maxNeighbors];
					if (tmpArray2==null || tmpArray2.length<maxNeighbors) tmpArray2 = new int[maxNeighbors];
					if (tmpArray3==null || tmpArray3.length<maxNeighbors) tmpArray3 = new int[maxNeighbors];
					buildComplementGraph(i+1);
					if (mode) {
						if (isNeighbor==null || isNeighbor.length<i+1) isNeighbor = new boolean[i+1];
						if (currentSolution==null || currentSolution.length<i+1) currentSolution = new boolean[i+1];
						if (tmpBoolean==null || tmpBoolean.length<i+1) tmpBoolean = new boolean[i+1];
						if (cliques==null) cliques = new Stream[i+1];
						else if (cliques.length<i+1) {
							Stream[] newCliques = new Stream[i+1];
							System.arraycopy(cliques,0,newCliques,0,cliques.length);
							for (j=0; j<cliques.length; j++) cliques[j]=null;
							cliques=newCliques;
						}
						totalSize+=maxScoreIndependentSet(basinType,i+1,minIntersection);
						if (printGraphs && i>0) drawConflictGraph(alignments,i+1,referenceLengths,alignmentsFile,currentSolution);
						for (j=0; j<=i; j++) bw.write(currentSolution[j]?"1\n":"0\n");
					}
					else {
						clawSize=getMaxClawSize(i+1);
						histogram[clawSize]++;
						maxClawSize=Math.max(maxClawSize,clawSize);
						System.err.println("buildConflictGraphs> current clawHistogram:");
						int x = clawHistogram.length-1;
						while (x>=0 && clawHistogram[x]==0) x--;
						for (int y=0; y<=x; y++) System.err.println(y+","+clawHistogram[y]);
						System.err.println();
						System.err.println("buildConflictGraphs> current clawHistogramPeriodic:");
						x=clawHistogramPeriodic.length-1;
						while (x>=0 && clawHistogramPeriodic[x]==0) x--;
						for (int y=0; y<=x; y++) System.err.println(y+","+clawHistogramPeriodic[y]);
						System.err.println();
					}
				}
				currentReadA=Alignments.readA; currentReadB=Alignments.readB;
				firstAlignmentID+=i+1;
				i=0; alignments[0].set(Alignments.readA-1,Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
				alignments[0].id=0;
				alignments[0].score=Alignments.endA-Alignments.startA+1-Alignments.diffs;
				str=br.readLine();
				continue;
			}
			i++;
			if (i==alignments.length) {
				FragmentsStep1.FFAlignment[] newAlignments = new FragmentsStep1.FFAlignment[alignments.length<<1];
				System.arraycopy(alignments,0,newAlignments,0,alignments.length);
				alignments=newAlignments;
			}
			if (alignments[i]==null) alignments[i] = new FragmentsStep1.FFAlignment(Alignments.readA-1,Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
			else alignments[i].set(Alignments.readA-1,Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
			alignments[i].id=i;
			alignments[i].score=Alignments.endA-Alignments.startA+1-Alignments.diffs;
			str=br.readLine();
		}
		br.close();
		if (currentReadA!=-1 && currentReadB!=-1) {
			maxNeighbors=buildConflictGraph(basinType,alignments,i+1,minIntersection);
			if (tmpArray1==null || tmpArray1.length<maxNeighbors) tmpArray1 = new int[maxNeighbors];
			if (tmpArray2==null || tmpArray2.length<maxNeighbors) tmpArray2 = new int[maxNeighbors];
			if (tmpArray3==null || tmpArray3.length<maxNeighbors) tmpArray3 = new int[maxNeighbors];
			buildComplementGraph(i+1);
			if (mode) {
				if (isNeighbor==null || isNeighbor.length<i+1) isNeighbor = new boolean[i+1];
				if (currentSolution==null || currentSolution.length<i+1) currentSolution = new boolean[i+1];
				if (tmpBoolean==null || tmpBoolean.length<i+1) tmpBoolean = new boolean[i+1];
				if (cliques==null) cliques = new Stream[i+1];
				else if (cliques.length<i+1) {
					Stream[] newCliques = new Stream[i+1];
					System.arraycopy(cliques,0,newCliques,0,cliques.length);
					for (j=0; j<cliques.length; j++) cliques[j]=null;
					cliques=newCliques;
				}
				totalSize+=maxScoreIndependentSet(basinType,i+1,minIntersection);
				if (printGraphs && i>0) drawConflictGraph(alignments,i+1,referenceLengths,alignmentsFile,currentSolution);
				for (j=0; j<=i; j++) bw.write(currentSolution[j]?"1\n":"0\n");
			}
			else {
				clawSize=getMaxClawSize(i+1);
				histogram[clawSize]++;
				maxClawSize=Math.max(maxClawSize,clawSize);
				System.err.println("buildConflictGraphs> current clawHistogram:");
				int x = clawHistogram.length-1;
				while (x>=0 && clawHistogram[x]==0) x--;
				for (int y=0; y<=x; y++) System.err.println(y+","+clawHistogram[y]);
				System.err.println();
				System.err.println("buildConflictGraphs> current clawHistogramPeriodic:");
				x=clawHistogramPeriodic.length-1;
				while (x>=0 && clawHistogramPeriodic[x]==0) x--;
				for (int y=0; y<=x; y++) System.err.println(y+","+clawHistogramPeriodic[y]);
				System.err.println();
			}
		}
		bw.close();
		return mode?totalSize:maxClawSize;
	}
	
	
	/**
	 * Stores in global variables $neighbors,lastNeighbor$ the conflict graph of all the 
	 * alignments in $alignments$. A conflict is either a containment, or an overlap of 
	 * length $>=minIntersection$, in either read, regardless of orientation. Node IDs
	 * equal alignment IDs.
	 *
	 * If the repeat is periodic, only containments/overlaps in the fragment are 
	 * considered conflicts. This is because a fragment that is longer than the reference 
	 * contains several copies of the same period, and every such copy might align best to
	 * the same or to overlapping substrings of the reference: all such copies should vote
	 * for the consensus of the substrings of the reference they align best to. This is 
	 * not correct for non-periodic repeats, since different substrings of the fragment
	 * should not vote for the same or for overlapping substrings of the reference.
	 *
	 * @param alignments readA=fragment ID; readB=reference ID;
	 * @param minIntersection assumed to be much smaller than $minAlignmentLength$;
	 * @return the max number of neighbors of a node in the conflict graph.
	 */
	private static final int buildConflictGraph(int basinType, FragmentsStep1.FFAlignment[] alignments, int nAlignments, int minIntersection) {
		final int CAPACITY = 10;  // Arbitrary
		final boolean IS_PERIODIC = basinType==Constants.INTERVAL_PERIODIC;
		int i, j, k;
		int last, startI, endI, idI, idJ, previous, maxNeighbors;
		
		if (neighbors==null || neighbors.length<nAlignments) neighbors = new int[nAlignments][CAPACITY];
		if (lastNeighbor==null || lastNeighbor.length<nAlignments) lastNeighbor = new int[nAlignments];
		Math.set(lastNeighbor,nAlignments-1,-1);
		if (nAlignments==1) return 0;
		
		// Detecting readA conflicts
		FragmentsStep1.FFAlignment.order=FragmentsStep1.FFAlignment.ORDER_READA_STARTA;
		Arrays.sort(alignments,0,nAlignments);
		for (i=0; i<nAlignments; i++) {
			startI=alignments[i].startA; endI=alignments[i].endA; idI=alignments[i].id;
			for (j=i+1; j<nAlignments; j++) {
				if (alignments[j].startA>endI-minIntersection+1) break;
				idJ=alignments[j].id;
				lastNeighbor[idI]++;
				if (lastNeighbor[idI]==neighbors[idI].length) {
					int[] newNeighbors = new int[neighbors[idI].length<<1];
					System.arraycopy(neighbors[idI],0,newNeighbors,0,neighbors[idI].length);
					neighbors[idI]=newNeighbors;
				}
				neighbors[idI][lastNeighbor[idI]]=idJ;
				lastNeighbor[idJ]++;
				if (lastNeighbor[idJ]==neighbors[idJ].length) {
					int[] newNeighbors = new int[neighbors[idJ].length<<1];
					System.arraycopy(neighbors[idJ],0,newNeighbors,0,neighbors[idJ].length);
					neighbors[idJ]=newNeighbors;
				}
				neighbors[idJ][lastNeighbor[idJ]]=idI;
			}
		}
		
		// Detecting readB conflicts
		if (!IS_PERIODIC) {
			FragmentsStep1.FFAlignment.order=FragmentsStep1.FFAlignment.ORDER_READB_STARTB;
			Arrays.sort(alignments,0,nAlignments);
			for (i=0; i<nAlignments; i++) {
				startI=alignments[i].startB; endI=alignments[i].endB; idI=alignments[i].id;
				for (j=i+1; j<nAlignments; j++) {
					if (alignments[j].startB>endI-minIntersection+1) break;
					idJ=alignments[j].id;
					lastNeighbor[idI]++;
					if (lastNeighbor[idI]==neighbors[idI].length) {
						int[] newNeighbors = new int[neighbors[idI].length<<1];
						System.arraycopy(neighbors[idI],0,newNeighbors,0,neighbors[idI].length);
						neighbors[idI]=newNeighbors;
					}
					neighbors[idI][lastNeighbor[idI]]=idJ;
					lastNeighbor[idJ]++;
					if (lastNeighbor[idJ]==neighbors[idJ].length) {
						int[] newNeighbors = new int[neighbors[idJ].length<<1];
						System.arraycopy(neighbors[idJ],0,newNeighbors,0,neighbors[idJ].length);
						neighbors[idJ]=newNeighbors;
					}
					neighbors[idJ][lastNeighbor[idJ]]=idI;
				}
			}
		}	
		
		// Sorting edges and removing duplicates
		maxNeighbors=0;
		for (i=0; i<nAlignments; i++) {
			if (lastNeighbor[i]<=0) continue; 
			Arrays.sort(neighbors[i],0,lastNeighbor[i]+1);
			previous=neighbors[i][0];
			k=0; last=lastNeighbor[i];
			for (j=1; j<=last; j++) {
				if (neighbors[i][j]==previous) continue;
				k++;
				neighbors[i][k]=neighbors[i][j];
				previous=neighbors[i][j];
			}
			lastNeighbor[i]=k;
			maxNeighbors=Math.max(maxNeighbors,k+1);
		}
		
		if (IO.CONSISTENCY_CHECKS) {
			int nEdges = 0;
			for (i=0; i<nAlignments; i++) nEdges+=lastNeighbor[i]+1;
			if (nEdges%2!=0) {
				System.err.println("buildConflictGraph> ERROR 1: asymmetric edges.");
				System.err.println("Conflict graph:");
				for (int x=0; x<nAlignments; x++) {
					System.err.print(x+": ");
					for (int y=0; y<=lastNeighbor[x]; y++) System.err.print(neighbors[x][y]+",");
					System.err.println();
				}
				System.exit(1);
			}
			for (i=0; i<nAlignments; i++) {
				for (j=0; j<=lastNeighbor[i]; j++) {
					int neighbor = neighbors[i][j];
					if (Arrays.binarySearch(neighbors[neighbor],0,lastNeighbor[neighbor]+1,i)<0) {
						System.err.println("buildConflictGraph> ERROR 2: "+neighbor+" is a neighbor of node "+i+", but node "+i+" is not a neighbor of node "+neighbor+"?!");
						for (int x=0; x<=lastNeighbor[neighbor]; x++) System.err.print(neighbors[neighbor][x]+",");
						System.err.println();
						System.exit(1);
					}
				}
			}
		}
		
		return maxNeighbors;
	}
	
	
	/**
	 * Stores in $IntervalGraph$ the complement of the conflict graph encoded by the 
	 * global variables $neighbors,lastNeighbor$ (i.e. a graph with the same set of nodes 
	 * as the original, but such that there is an edge in the complement iff there is no 
	 * edge in the original). This is used for finding claws, since an independent set in 
	 * the original graph is a clique in the complement, and there is already an optimized
	 * clique enumeration procedure in $IntervalGraph.getMaximalCliques()$. Node IDs equal
	 * alignment IDs.
	 *
	 * Remark: the interval graph built by the procedure is not a fully-fledged interval
	 * graph. It is just enough to make $IntervalGraph.getMaximalCliques()$ work.
	 *
	 * Remark: the procedure uses the global variables $tmpArray{1,2},edgePool$.
	 */
	private static final void buildComplementGraph(int nNodes) {
		final int INCREMENT = 1000000;  // Arbitrary
		int i, j, k, h;
		int last2, neighbor, length;
		
		if (tmpArray1==null || tmpArray1.length<nNodes) tmpArray1 = new int[nNodes];
		if (tmpArray2==null || tmpArray2.length<nNodes) tmpArray2 = new int[nNodes];
		for (i=0; i<nNodes; i++) tmpArray1[i]=i;
		IntervalGraph.nNodes=nNodes;
		if (IntervalGraph.neighbors==null) IntervalGraph.neighbors = new IntervalGraph.Edge[nNodes][0];
		else if (IntervalGraph.neighbors.length<nNodes) {
			IntervalGraph.Edge[][] newNeighbors = new IntervalGraph.Edge[nNodes][0];
			length=IntervalGraph.neighbors.length;
			for (i=0; i<length; i++) {
				newNeighbors[i]=IntervalGraph.neighbors[i];
				IntervalGraph.neighbors[i]=null;
			}
			IntervalGraph.neighbors=newNeighbors;
		}
		for (i=0; i<nNodes; i++) {
			length=IntervalGraph.neighbors[i].length;
			for (j=0; j<length; j++) IntervalGraph.neighbors[i][j]=null;
		}
		if (IntervalGraph.nNeighbors==null || IntervalGraph.nNeighbors.length<nNodes) IntervalGraph.nNeighbors = new int[nNodes];
		Math.set(IntervalGraph.nNeighbors,nNodes-1,-1);
		lastInPool=-1;
		for (i=0; i<nNodes; i++) {
			last2=Math.setMinus(tmpArray1,nNodes-1,neighbors[i],lastNeighbor[i],tmpArray2);
			if (IntervalGraph.neighbors[i]==null) IntervalGraph.neighbors[i] = new IntervalGraph.Edge[last2+1];
			else if (IntervalGraph.neighbors[i].length<last2+1) {
				IntervalGraph.Edge[] newNeighbors = new IntervalGraph.Edge[last2+1];
				for (j=0; j<IntervalGraph.neighbors[i].length; j++) IntervalGraph.neighbors[i][j]=null;
				IntervalGraph.neighbors[i]=newNeighbors;
			}
			k=-1;
			for (j=0; j<=last2; j++) {
				neighbor=tmpArray2[j];
				if (neighbor==i) continue;
				lastInPool++;
				if (lastInPool==edgePool.length) {
					IntervalGraph.Edge[] newPool = new IntervalGraph.Edge[edgePool.length+INCREMENT];
					System.arraycopy(edgePool,0,newPool,0,edgePool.length);
					for (h=edgePool.length; h<newPool.length; h++) newPool[h] = new IntervalGraph.Edge();
					for (h=0; h<edgePool.length; h++) edgePool[h]=null;
					edgePool=newPool;
				}
				k++;
				IntervalGraph.neighbors[i][k]=edgePool[lastInPool];
				IntervalGraph.neighbors[i][k].nodeID1=Math.min(i,neighbor);
				IntervalGraph.neighbors[i][k].nodeID2=Math.max(i,neighbor);
				IntervalGraph.neighbors[i][k].on=true;
			}
			IntervalGraph.nNeighbors[i]=k+1;
		}
	}
	
	
	/**
	 * A set of nodes $V$ is a \emph{claw} iff its induced subgraph is a star. The nodes
	 * that are not the center of the star are called \emph{talons}, and the size of a
	 * claw is its number of talons. For every node $i$ in the conflict graph, the 
	 * procedure computes the size of a largest claw whose center is $i$, and it returns
	 * the max overall.
	 *
	 * Most conflict graphs of non-periodic and periodic repeats have claws of size <=5
	 * (see $clawSize.png$). However, we observed graphs of periodic repeats with claws of
	 * size up to 15.
	 *
	 * Remark: if no two alignments are contained in one another in any dimension, and if
	 * the overlap threshold is one, then the conflict graph is 5-claw-free 
	 * \cite{bafna1996nonoverlapping}. It is easy to see that this holds even with a 
	 * threshold >1. However, allowing containment breaks this property.
	 *
	 * Remark: the procedure assumes that $buildComplementGraph()$ has already been
	 * called, since it works by enumerating all the maximal cliques in the subgraph of
	 * the complement graph induced by the neighbors of $i$.
	 *
	 * Remark: the procedure can be very slow in periodic repeats, since both the conflict
	 * graph and its complement can be dense.
	 * 
	 * @return 0 if there is just one node in the graph.
	 */
	private static final int getMaxClawSize(int nNodes) {
		int i;
		int last, last2, neighbor, maxCliqueSize, out;
		
		if (nNodes==1) return 0;
		out=1;
		for (i=0; i<nNodes; i++) {
			last=lastNeighbor[i];
			if (last==0) continue;
			IntervalGraph.getMaximalCliques(neighbors[i],0,lastNeighbor[i],null,-1,-1,Math.POSITIVE_INFINITY,false);
			maxCliqueSize=IntervalGraph.maxCliqueSize;
			System.err.println("getMaxClawSize> node "+i+"/"+nNodes+" nNeighbors="+(lastNeighbor[i]+1)+" maxCliqueSize="+maxCliqueSize);
			if (IO.CONSISTENCY_CHECKS && maxCliqueSize>lastNeighbor[i]+1) {
				System.err.println("getMaxClawSize> ERROR: maxCliqueSize="+maxCliqueSize+" > nNeighbors="+(lastNeighbor[i]+1));
				System.exit(1);
			}
			out=Math.max(out,maxCliqueSize);
		}
		return out;
	}
	
	
	/**
	 * Finds a high-scoring independent set of the conflict graph, by first greedily 
	 * adding a node that maximizes the score. The greedy solution has score >=S/(d-1) in 
	 * graphs without d-claws, where S is the optimum. Then, the procedure improves the 
	 * greedy solution by iteratively swapping the talons of a claw with their neighbors 
	 * in the solution, if this improves the score. We use the algorithm in 
	 * \cite{berman2000d}, which applies any improving claw, since enumerating all claws 
	 * to find the one that maximizes some function (e.g. as done in 
	 * \cite{chandra2001greedy}) is too slow in practice. This should give a solution with
	 * score >=S/(d/2).
	 *
	 * Remark: alternatively, one could iteratively look for an independent set of a fixed 
	 * size $t$ outside the current solution, and then swap it with its neighbors in the 
	 * current solution if this improves the score \cite{bafna1996nonoverlapping}. This is
	 * probably faster in practice, but it might achieve a worse approximation ratio. Of 
	 * course one could also express the problem as an ILP and feed it to a solver, or 
	 * use a heuristic on top of an LP relaxation \cite{chan2012approximation}.
	 *
	 * Remark: better techniques might exist if we take into account that our rectangles
	 * are approximately squares (e.g. \cite{erlebach2005polynomial}).
	 *
	 * Remark: the procedure tries to consider small claws first, since enumerating large 
	 * claws is expensive.
	 *
	 * Remark: the procedure uses global variables $tmpArray*,tmpBoolean$.
	 *
	 * @return the number of nodes in the solution.
	 */
	private static final int maxScoreIndependentSet(int basinType, int nNodes, int minIntersection) {
		final double SAMPLING_RATE = 0.05;  // To limit the time of $IntervalGraph.getMaximalCliques()$ in dense graphs. Arbitrary.
		final int MAX_N_CLIQUES = 10000;  // To limit the time of $IntervalGraph.getMaximalCliques()$ in dense graphs. Arbitrary.
		final boolean IS_PERIODIC = basinType==Constants.INTERVAL_PERIODIC;
		int i, j;
		int last, lastClique, currentSize, alignmentID;
		int neighbor, nNeighbors, maxNode;
		long greedyScore, improvedScore, currentScore, deltaScore, maxScore, nTraversedCliques;
		
		if (nNodes==1) {
			currentSolution[0]=true;
			return 1;
		}
		
		// First greedy pass
		System.err.println("maxScoreIndependentSet> Computing the greedy solution... ");
		FragmentsStep1.FFAlignment.order=FragmentsStep1.FFAlignment.ORDER_ID;
		Arrays.sort(alignments,0,nNodes);
		Math.set(isNeighbor,nNodes-1,false);
		Math.set(currentSolution,nNodes-1,false);
		currentSize=0; currentScore=0; nNeighbors=0;
		while (currentSize+nNeighbors<nNodes) {
			maxScore=0; maxNode=-1;
			for (i=0; i<nNodes; i++) {
				if (currentSolution[i] || isNeighbor[i]) continue;
				if (alignments[i].score>maxScore) {
					maxScore=alignments[i].score;
					maxNode=i;
				}
			}
			currentSolution[maxNode]=true;
			currentSize++;
			currentScore+=maxScore*maxScore;
			for (i=0; i<=lastNeighbor[maxNode]; i++) {
				neighbor=neighbors[maxNode][i];
				if (!isNeighbor[neighbor]) {
					nNeighbors++;
					isNeighbor[neighbor]=true;
				}
			}
		}
		System.err.println("maxScoreIndependentSet> The greedy solution has score "+currentScore+" and "+currentSize+" nodes ("+IO.getPercent(currentSize,nNodes)+"% of all nodes).");		
		if (IO.CONSISTENCY_CHECKS) checkSolution(basinType,currentSolution,nNodes,minIntersection,"greedy");
		greedyScore=currentScore;
		
		// Improvements
		System.err.println("maxScoreIndependentSet> Improving the greedy solution... ");
		for (i=0; i<nNodes; i++) alignments[i].nNeighbors=IntervalGraph.nNeighbors[i];
		FragmentsStep1.FFAlignment.order=FragmentsStep1.FFAlignment.ORDER_NNEIGHBORS;
		Arrays.sort(alignments,0,nNodes);
		if (tmpArray1==null || tmpArray1.length<nNodes) tmpArray1 = new int[nNodes];
		Math.set(tmpArray1,nNodes-1,-1);
		for (i=0; i<nNodes; i++) tmpArray1[alignments[i].id]=i;
		Math.set(tmpBoolean,nNodes-1,false);
		cliques = new Stream[nNodes]; lastClique=-1; nTraversedCliques=0;
		while (true) {
			deltaScore=0;
			for (i=0; i<nNodes; i++) {			
				alignmentID=alignments[i].id;
				last=lastNeighbor[alignmentID];
				if (last==0) deltaScore=maxScoreIndependentSet_impl_oneTalon(neighbors[alignmentID][0],currentSolution,currentScore,tmpArray1,tmpArrayLong);
				else {
					if (i>lastClique) {
						cliques[i] = new Stream(256);
						if (IS_PERIODIC) IntervalGraph.getMaximalCliques(neighbors[alignmentID],0,lastNeighbor[alignmentID],cliques[i],-1,SAMPLING_RATE,MAX_N_CLIQUES,false);
						else IntervalGraph.getMaximalCliques(neighbors[alignmentID],0,lastNeighbor[alignmentID],cliques[i],-1,-1,Math.POSITIVE_INFINITY,false);
						lastClique=i;
					}
					deltaScore=maxScoreIndependentSet_impl(cliques[i],currentSolution,currentScore,tmpBoolean,tmpArray1,tmpArrayLong);
					nTraversedCliques+=tmpArrayLong[0];
				}
				if (deltaScore>0) {
					currentScore+=deltaScore;
					currentSize+=(int)(tmpArrayLong[1]-tmpArrayLong[2]);
					System.err.println("maxScoreIndependentSet> A claw of size "+tmpArrayLong[1]+" with center "+alignmentID+" ("+i+"/"+nNodes+") increases the score by "+deltaScore);
					break;
				}
				if (nTraversedCliques>=MAX_N_CLIQUES) break;
			}
			if (deltaScore==0) break;
		}
		System.err.println("maxScoreIndependentSet> The refined solution has score "+currentScore+" and "+currentSize+" nodes ("+IO.getPercent(currentSize,nNodes)+"% of all nodes).");
		if (currentScore==greedyScore) return currentSize;
		if (IO.CONSISTENCY_CHECKS) checkSolution(basinType,currentSolution,nNodes,minIntersection,"refined");
		improvedScore=currentScore;
		
		// Second greedy pass. This is useful, since improvements do not guarantee the
		// solution to be maximal, i.e. a vertex cover (e.g. a node $v$ might be adjacent
		// only to a node $w$ in the greedy solution; then $w$ might get swapped out of
		// the refined solution, while $v$ might still not be in the refined solution).
		System.err.println("maxScoreIndependentSet> Ensuring that the refined solution is maximal... ");
		FragmentsStep1.FFAlignment.order=FragmentsStep1.FFAlignment.ORDER_ID;
		Arrays.sort(alignments,0,nNodes);
		Math.set(isNeighbor,nNodes-1,false); nNeighbors=0;
		for (i=0; i<nNodes; i++) {
			if (!currentSolution[i]) continue;
			last=lastNeighbor[i];
			for (j=0; j<=last; j++) {
				neighbor=neighbors[i][j];
				if (!isNeighbor[neighbor]) {
					nNeighbors++;
					isNeighbor[neighbor]=true;
				}
			}
		}
		while (currentSize+nNeighbors<nNodes) {
			maxScore=0; maxNode=-1;
			for (i=0; i<nNodes; i++) {
				if (currentSolution[i] || isNeighbor[i]) continue;
				if (alignments[i].score>maxScore) {
					maxScore=alignments[i].score;
					maxNode=i;
				}
			}
			currentSolution[maxNode]=true;
			currentSize++;
			currentScore+=maxScore*maxScore;
			for (i=0; i<=lastNeighbor[maxNode]; i++) {
				neighbor=neighbors[maxNode][i];
				if (!isNeighbor[neighbor]) {
					nNeighbors++;
					isNeighbor[neighbor]=true;
				}
			}
		}
		System.err.println("maxScoreIndependentSet> The final solution has score "+currentScore+" and "+currentSize+" nodes ("+IO.getPercent(currentSize,nNodes)+"% of all nodes).");		
		if (IO.CONSISTENCY_CHECKS && currentScore>improvedScore) checkSolution(basinType,currentSolution,nNodes,minIntersection,"final");
		
		return currentSize;
	}

	
	/**
	 * A variant of $IntervalGraph.markMaximalCliques()$ that stops at the first clique
	 * that, if applied to the current solution, improves the sum of squared scores
	 * \cite{berman2000d}; the procedure applies such a clique to the current solution.
	 *
	 * Remark: the procedure uses the global variable $stack$.
	 *
	 * @param currentSolution represented as a bitvector flagging the nodes;
	 * @param tmpBoolean temporary array, of size at least equal to the number of nodes;
	 * the procedure assumes it to be set to all FALSE, and resets it at the end;
	 * @param old2new $alignments$ might be sorted by a different criterion than ID: this
	 * array stores the position of the alignment with ID $i$ in the current version of
	 * $alignments$;
	 * @param out output array: 0=n. of cliques examined; 1=size of the clique used to 
	 * improve the current solution; 2=size of its neighborhood in the solution;
	 * @return the increase in score (zero if no clique increases the score).
	 */
	private static final long maxScoreIndependentSet_impl(Stream in, boolean[] currentSolution, long currentScore, boolean[] tmpBoolean, int[] old2new, long[] out) {
		boolean fullyContained;
		int i, j, k;
		int top, element, node, neighbor, score, nNeighbors;
		long nCliques, deltaScore;
		final int nElements = in.nElements();
		
		top=-1; nCliques=0;
		for (i=0; i<nElements; i++) {
			element=in.getElementAt(i);
			if (element==IntervalGraph.MAXIMAL_CLIQUE_FOUND) {
				nCliques++;
				// Discarding the clique if it is fully contained in the current solution
				fullyContained=true;
				for (j=0; j<=top; j++) {
					if (!currentSolution[stack[j]]) {
						fullyContained=false;
						break;
					}
				}
				if (fullyContained) {
					top--;
					continue;
				}				
				// Measuring the change in score
				deltaScore=0; nNeighbors=0;
				for (j=0; j<=top; j++) {
					node=stack[j];
					if (!currentSolution[node]) {
						score=alignments[old2new[node]].score;
						deltaScore+=score*score;
					}
					for (k=0; k<=lastNeighbor[node]; k++) {
						neighbor=neighbors[node][k];
						if (currentSolution[neighbor] && !tmpBoolean[neighbor]) {
							nNeighbors++;
							score=alignments[old2new[neighbor]].score;
							deltaScore-=score*score;
							tmpBoolean[neighbor]=true;
						}
					}
				}
				if (deltaScore<=0) {
					for (j=0; j<=top; j++) {
						node=stack[j];
						for (k=0; k<=lastNeighbor[node]; k++) tmpBoolean[neighbors[node][k]]=false;
					}
					top--;
					continue;
				}
				// Updating the current solution
				for (j=0; j<=top; j++) {
					node=stack[j];
					for (k=0; k<=lastNeighbor[node]; k++) currentSolution[neighbors[node][k]]=false;
				}
				for (j=0; j<=top; j++) currentSolution[stack[j]]=true;
				out[0]=nCliques; out[1]=top+1; out[2]=nNeighbors;
				// Cleaning $tmpBoolean$.
				for (j=0; j<=top; j++) {
					node=stack[j];
					for (k=0; k<=lastNeighbor[node]; k++) tmpBoolean[neighbors[node][k]]=false;
				}
				return deltaScore;
			}
			else if (element==IntervalGraph.MAXIMAL_CLIQUE_UP) {
				top--;
				continue;
			}
			else if (element==IntervalGraph.MAXIMAL_CLIQUE_NOT_FOUND) {
				top--;
				continue;
			}
			else {
				top++;
				if (top==stack.length) {
					int[] newStack = new int[stack.length<<1];
					System.arraycopy(stack,0,newStack,0,stack.length);
					stack=newStack;
				}
				stack[top]=element;
				/*if (IO.CONSISTENCY_CHECKS) {   // Commented out since it is working
					for (int x=0; x<top; x++) {
						if (stack[x]==element) {
							System.err.println("maxScoreIndependentSet_impl> ERROR: duplicated element in the stack?!");
							for (int y=0; y<=top; y++) System.err.print(stack[y]+", ");
							System.err.println();
							System.exit(1);
						}
					}
				}*/
			}
		}
		return 0;
	}
	
	
	/**
	 * Variant of $maxScoreIndependentSet_impl()$ for claws with one talon ($node$ is such
	 * a talon).
	 */
	private static final long maxScoreIndependentSet_impl_oneTalon(int node, boolean[] currentSolution, long currentScore, int[] old2new, long[] out) {
		int i;
		int neighbor, nNeighbors, score;
		final int last = lastNeighbor[node];
		long deltaScore;
		
		// Discarding the clique if it is fully contained in the current solution
		if (currentSolution[node]) return 0;
		
		// Measuring the change in score
		deltaScore=0; nNeighbors=0;
		score=alignments[old2new[node]].score;
		deltaScore+=score*score;
		for (i=0; i<=last; i++) {
			neighbor=neighbors[node][i];
			if (currentSolution[neighbor]) {
				nNeighbors++;
				score=alignments[old2new[neighbor]].score;
				deltaScore-=score*score;
			}
		}
		if (deltaScore<=0) return 0;
		
		// Updating the current solution
		for (i=0; i<=last; i++) currentSolution[neighbors[node][i]]=false;
		currentSolution[node]=true;
		out[0]=0; out[1]=1; out[2]=nNeighbors;
		return deltaScore;
	}
	
	
	/**
	 * Checks if the solution has conflicts.
	 * 
	 * Remark: the procedure uses global variable $tmpArray3$.
	 */
	private static final void checkSolution(int basinType, boolean[] currentSolution, int nNodes, int minIntersection, String solutionID) {
		final boolean IS_PERIODIC = basinType==Constants.INTERVAL_PERIODIC;
		int i, j;
		int last, nodeI, nodeJ, startA, endA, startB, endB;
		
		if (nNodes>1) {
			FragmentsStep1.FFAlignment.order=FragmentsStep1.FFAlignment.ORDER_ID;
			Arrays.sort(alignments,0,nNodes);
		}
		if (tmpArray3==null || tmpArray3.length<nNodes) tmpArray3 = new int[nNodes];
		last=-1;
		for (i=0; i<nNodes; i++) {
			if (!currentSolution[i]) continue;
			tmpArray3[++last]=i;
		}
		if (last==0) return;
		else if (last==-1) {
			System.err.println("ERROR: the "+solutionID+" solution has zero nodes?!");
			System.exit(1);
		}
		for (i=0; i<=last; i++) {
			nodeI=tmpArray3[i];
			startA=alignments[nodeI].startA; endA=alignments[nodeI].endA;
			startB=alignments[nodeI].startB; endB=alignments[nodeI].endB;
			for (j=i+1; j<=last; j++) {
				nodeJ=tmpArray3[j];
				if (Intervals.intersectionLength(startA,endA,alignments[nodeJ].startA,alignments[nodeJ].endA)>=minIntersection) {
					System.err.println("ERROR: the "+solutionID+" solution contains overlapping alignments in the fragment:");
					System.err.println(alignments[nodeI]);
					System.err.println(alignments[nodeJ]);
					System.exit(1);
				}
				if (!IS_PERIODIC && Intervals.intersectionLength(startB,endB,alignments[nodeJ].startB,alignments[nodeJ].endB)>=minIntersection) {
					System.err.println("ERROR: the "+solutionID+" solution contains overlapping alignments in the reference:");
					System.err.println(alignments[nodeI]);
					System.err.println(alignments[nodeJ]);
					System.exit(1);
				}
			}
		}
	}
	
	
	private static final Stroke EDGE_STROKE = new BasicStroke(1);
	private static final Stroke DEFAULT_STROKE = new BasicStroke(2);
	private static final Stroke SELECTED_STROKE = new BasicStroke(8);
	private static final Color DEFAULT_COLOR = new Color(Colors.GRADIENT_TO);
	private static final Color SELECTED_COLOR = new Color(Colors.GRADIENT_FROM);
	private static final Color BACKGROUND_COLOR = new Color(Colors.COLOR_BACKGROUND);
	private static final Color EDGE_COLOR = new Color(Colors.COLOR_BACKGROUND_LIGHT);
	private static final int CIRCLE_RADIUS_PIXELS = 10;
	private static final int CIRCLE_DIAMETER_PIXELS = (CIRCLE_RADIUS_PIXELS)<<1;
	private static final int SCALE_FACTOR = 10;
	private static final double SCORE_FACTOR = 0.1;
	
	
	/**
	 * Draws a conflict graph as a set of rectangles. The dimensions of the image are the
	 * length of the fragment (derived from $Reads.getReadLength()$) and the length of the
	 * reference (derived from $referenceLengths$). This idea comes from 
	 * \cite{bafna1996nonoverlapping}.
	 *
	 * @param selected flags the alignments selected by $maxScoreIndependentSet()$;
	 * @param prefix prefix of the output file.
	 */
	private static final void drawConflictGraph(FragmentsStep1.FFAlignment[] alignments, int nAlignments, int[] referenceLengths, String prefix, boolean[] selected) throws IOException {
		final int readA = alignments[0].readA;
		final int readB = alignments[0].readB;
		int i, j;
		int last, neighbor, centerA, centerB, centerAPrime, centerBPrime;
		int circleRadius, circleDiameter;
		BufferedImage image;
		Graphics2D graphics;
		
		image = new BufferedImage(Math.ceil(Reads.getReadLength(readA),SCALE_FACTOR),Math.ceil(referenceLengths[readB],SCALE_FACTOR),BufferedImage.TYPE_INT_RGB);
		graphics=image.createGraphics();
		graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		graphics.setBackground(BACKGROUND_COLOR);
		if (nAlignments>0) {
			FragmentsStep1.FFAlignment.order=FragmentsStep1.FFAlignment.ORDER_ID;
			Arrays.sort(alignments,0,nAlignments);
		}
		// Drawing edges
		graphics.setStroke(EDGE_STROKE); graphics.setColor(EDGE_COLOR);
		for (i=0; i<nAlignments; i++) {
			centerA=(alignments[i].startA+alignments[i].endA+1)>>1;
			centerB=(alignments[i].startB+alignments[i].endB+1)>>1;
			last=lastNeighbor[i];
			for (j=0; j<=last; j++) {
				neighbor=neighbors[i][j];
				centerAPrime=(alignments[neighbor].startA+alignments[neighbor].endA+1)>>1;
				centerBPrime=(alignments[neighbor].startB+alignments[neighbor].endB+1)>>1;
				graphics.drawLine(centerA/SCALE_FACTOR,centerB/SCALE_FACTOR,centerAPrime/SCALE_FACTOR,centerBPrime/SCALE_FACTOR);
			}
		}
		// Drawing rectangles
		for (i=0; i<nAlignments; i++) {
			if (selected[i]) {
				graphics.setStroke(SELECTED_STROKE);
				graphics.setColor(SELECTED_COLOR);
			}
			else {
				graphics.setStroke(DEFAULT_STROKE);
				graphics.setColor(DEFAULT_COLOR);
			}
			graphics.drawRect(alignments[i].startA/SCALE_FACTOR,alignments[i].startB/SCALE_FACTOR,(alignments[i].endA-alignments[i].startA+1)/SCALE_FACTOR,(alignments[i].endB-alignments[i].startB+1)/SCALE_FACTOR);
			centerA=(alignments[i].startA+alignments[i].endA+1)>>1;
			centerB=(alignments[i].startB+alignments[i].endB+1)>>1;
			circleRadius=(int)((alignments[i].score/SCALE_FACTOR)*SCORE_FACTOR);
			circleDiameter=circleRadius<<1;
			graphics.fillOval(centerA/SCALE_FACTOR-circleRadius,centerB/SCALE_FACTOR-circleRadius,circleDiameter,circleDiameter);
		}
		ImageIO.write(image,"png",new File(prefix+"-"+readA+"-"+readB+".png"));
	}
	
	
	/**
	 * Prints the conflict graph and its complement to STDERR. Node IDs equal alignment 
	 * IDs.
	 */
	private static final void printConflictGraph(int nNodes) {
		int i, j;
		int last;
		
		System.err.println("CONFLICT GRAPH:");
		System.err.println("graph G {");
		for (i=0; i<nNodes; i++) {
			last=lastNeighbor[i];
			for (j=0; j<=last; j++) System.err.println(i+" -- "+neighbors[i][j]+";");
		}
		System.err.println("}");
		System.err.println();
		System.err.println("COMPLEMENT GRAPH:");		
		System.err.println("graph G {");
		for (i=0; i<IntervalGraph.nNodes; i++) {
			last=IntervalGraph.nNeighbors[i]-1;
			for (j=0; j<=last; j++) System.err.println(i+" -- "+IntervalGraph.neighbors[i][j].getTo(i)+";");
		}
		System.err.println("}");
	}
	
	
	/**
	 * Uses the bitvectors built for the fragment-reference alignments file, to build the 
	 * bitvectors of every corresponding reference-fragment alignment file.
	 *
	 * @param path root directory of fragments-reference alignments.
	 */
	private static final void buildSymmetricBitvectors(String path) throws IOException {
		int i;
		int nOnes;
		String str1, str2;
		BufferedReader br1, br2;
		BufferedWriter bw;
		File file;
		FragmentsStep1.FFAlignment tmpAlignment, newAlignment;
		HashSet<FragmentsStep1.FFAlignment> set;
		String[] list;
		
		tmpAlignment = new FragmentsStep1.FFAlignment();
		set = new HashSet<FragmentsStep1.FFAlignment>(10000);
		file = new File(path);
		list=file.list();
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			file = new File(path+"/"+list[i]+"/"+IO.REFERENCE_FRAGMENTS_FILE+".txt");
			if (!file.exists()) {
				System.err.println("WARNING: file not found: "+path+"/"+list[i]+"/"+IO.REFERENCE_FRAGMENTS_FILE+".txt");
				continue;
			}
			// Loading fragment-reference alignments
			set.clear();
			br1 = new BufferedReader(new FileReader(path+"/"+list[i]+"/"+IO.FRAGMENTS_REFERENCE_FILE+".txt"));
			br2 = new BufferedReader(new FileReader(path+"/"+list[i]+"/"+IO.FRAGMENTS_REFERENCE_FILE+"-keep.txt"));
			str1=br1.readLine(); str1=br1.readLine();  // Skipping header
			str1=br1.readLine();
			str2=br2.readLine();
			while (str2!=null) {
				if (str2.charAt(0)=='1') {
					Alignments.readAlignmentFile(str1);
					newAlignment = new FragmentsStep1.FFAlignment(Alignments.readA,Alignments.readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
					newAlignment.inCanonicalForm();
					set.add(newAlignment);
				}
				str1=br1.readLine(); str2=br2.readLine();
			}
			br1.close(); br2.close();
			
			// Filtering reference-fragment alignments
			br1 = new BufferedReader(new FileReader(path+"/"+list[i]+"/"+IO.REFERENCE_FRAGMENTS_FILE+".txt"));
			bw = new BufferedWriter(new FileWriter(path+"/"+list[i]+"/"+IO.REFERENCE_FRAGMENTS_FILE+"-keep.txt"));
			str1=br1.readLine(); str1=br1.readLine();  // Skipping header
			str1=br1.readLine(); nOnes=0;
			while (str1!=null) {
				Alignments.readAlignmentFile(str1);
				tmpAlignment.set(Alignments.readA,Alignments.readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
				tmpAlignment.inCanonicalForm();
				if (set.contains(tmpAlignment)) {
					nOnes++;
					bw.write("1\n");
				}
				else bw.write("0\n");
				str1=br1.readLine();
			}
			br1.close(); bw.close();
			if (nOnes!=set.size()) {
				System.err.println("ERROR: written "+nOnes+" ones instead of "+set.size()+", module "+list[i]);
				System.exit(1);
			}
		}
		set.clear();
	}
	

}