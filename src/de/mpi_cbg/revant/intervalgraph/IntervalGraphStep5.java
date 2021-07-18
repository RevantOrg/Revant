package de.mpi_cbg.revant.intervalgraph;

import java.io.*;
import java.util.Arrays;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Leaf;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Leaves;
import de.mpi_cbg.revant.util.DAG;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Colors;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Intervals;
import de.mpi_cbg.revant.factorize.Factorize;


/**
 * 
 *  
 * 
 * Remark: for simplicity, the program does not create new kernel tags files, but it just 
 * filters the existing kernel descriptor files.
 */
public class IntervalGraphStep5 {
	
	/** 
	 * Global data structures
	 */
	private static int nNodes, nPathKernels;
	private static byte[] pathKernelPeriodic;
	private static int[] pathKernelLengths;
	private static IntervalGraph.Node[] nodesArray;
	private static int[][] kernel2descriptor;
	private static byte[] labels;
	
	/**
	 * Temporary space for the $markExposed()$ procedure.
	 */
	private static int lastKernelWithTrack, lastSelectedKernel;
	private static int[][] counts, tracks, cTrackStart, cTrackEnd;
	private static int[] lastTrack, lastCTrack;
	private static int[] kernelsWithTrack, selectedKernels;
	private static long[] kernelSurface;
	private static IntervalGraph.Node[] tmpWindows;
	
	/**
	 * Temporary space
	 */
	private static byte[] tmpByte1 = new byte[50];
	private static int[] tmpArray1 = new int[50];
	private static int[] tmpArray2 = new int[50];
	private static int[] tmpArray3 = new int[50];
	private static int[] tmpIO = new int[10];
	private static long[] tmpLong = new long[2];
	private static boolean[] tmpBoolean = new boolean[50];
	private static Point[] tmpPoints;
	private static Leaf[] tmpLeaves;
	
	
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		IO.coverage=Integer.parseInt(args[1]);
		final String QUALITY_THRESHOLDS_FILE = args[2];
		final int N_READS = Integer.parseInt(args[3]);
		Alignments.minAlignmentLength=Integer.parseInt(args[4]);
		final String QUALITIES_FILE = args[5].equalsIgnoreCase("null")?null:args[5];
		final int MEDIAN_READ_LENGTH = Integer.parseInt(args[6]);
		
		final double MIN_SURFACE_FRACTION = 0.1;  // Arbitrary
		IO.minOccurrencesInGenome=3;  // Arbitrary
		IO.minRepeatCoverage=(IO.minOccurrencesInGenome*IO.coverage)-1;
		final String READ_LENGTHS_FILE = STEP1_DIR+"/../"+IO.READS_LENGTHS;
		final String READS_IDS_FILE = STEP1_DIR+"/../"+IO.READS_IDS;
		final String OUTPUT_DIR = STEP1_DIR+"/"+IO.TAGS_DIR+"/"+IO.STEP4_DIR+"/"+IO.STEP5_DIR;
		final int MIN_EXPOSED_SURFACE = Math.min(IO.quantum<<2,Alignments.minAlignmentLength>>1);  // Arbitrary
		
		int i;
		int nRemoved, nMarked, nSelected, nSelectedShortPeriod, nSelectedLongPeriod;
		BufferedWriter bw;
		
		// Loading kernels
		IO.initialize();
		Reads.nReads=N_READS;
		Reads.loadReadIDs(READS_IDS_FILE,N_READS);
		Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadQualities(QUALITY_THRESHOLDS_FILE,QUALITIES_FILE);
		i=loadNodes(STEP1_DIR);
		System.err.println("Step5> "+nNodes+" total nodes and " +nPathKernels+" path kernels loaded, i="+i);
		if (i==-1) return;
		tmpPoints = new Point[nPathKernels];
		for (i=0; i<nPathKernels; i++) tmpPoints[i] = new Point();
		tmpLeaves = new Leaf[nPathKernels];
		for (i=0; i<nPathKernels; i++) tmpLeaves[i] = new Leaf();
		
		// Filtering kernel descriptors
		labels = new byte[nPathKernels];
		nMarked=markExposed(MEDIAN_READ_LENGTH,MIN_EXPOSED_SURFACE);
		System.err.println("STEP5> "+nMarked+" kernels kept because exposed ("+IO.getPercent(nMarked,nPathKernels)+"%)");
		if (IO.SHOW_INTERACTIVE) {
			System.err.println("labels:");
			for (int x=0; x<nPathKernels; x++) System.err.println(x+","+labels[x]);
		}		
		keepNonexposed(true);
		nRemoved=markSameSurface(MIN_SURFACE_FRACTION,true,true);
		System.err.println("STEP5> Removed "+nRemoved+" short-period path kernels with similar surface ("+IO.getPercent(nRemoved,nMarked)+"% of surviving kernels)");		
		nRemoved=markSameSurface(MIN_SURFACE_FRACTION,false,true);
		System.err.println("STEP5> Removed "+nRemoved+" long-period path kernels with similar surface ("+IO.getPercent(nRemoved,nMarked)+"% of surviving kernels)");
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		if (tmpByte1==null || tmpByte1.length<nPathKernels) tmpByte1 = new byte[nPathKernels];
		for (i=0; i<IntervalGraph.nNodes; i++) removeKernels(IntervalGraph.nodesArray[i]);
		bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+IO.STEP5_DESCRIPTORS_FILE));
		nSelected=0; nSelectedShortPeriod=0; nSelectedLongPeriod=0;
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]!=0) {
				bw.write(IO.BASIN_DESCRIPTOR_PREFIX+"-"+kernel2descriptor[i][0]+"-"+kernel2descriptor[i][1]+"-"+kernel2descriptor[i][2]+".txt \n");
				nSelected++;
				if (pathKernelPeriodic[i]==1) nSelectedShortPeriod++;
				else if (pathKernelPeriodic[i]==2) nSelectedLongPeriod++;
			}
		}
		bw.close();
		System.err.println("STEP5> Selected "+nSelected+" path kernels: "+nSelectedShortPeriod+" short-period, "+nSelectedLongPeriod+" long-period.");
		
		// Writing adjacency information
		markAdjacent(OUTPUT_DIR+"/"+IO.STEP5_ADJACENCY_FILE);
	}
	
	
	/**
	 * Wraps $IntervalGraphStep4.loadNodes()$ to load all nodes in all connected 
	 * components, and the basins of all path kernel descriptors.
	 *
	 * Remark: the procedure sorts $nodesArray$ by $read,start$.
	 *
	 * Remark: only intervals that belong to a kernel descriptor receive a kernel tag, 
	 * even though there might have been other tags in the original graph (e.g. a 
	 * connected component of the Step1 graph might have been too infrequent to be kept, 
	 * but it was still assigned a tag). So e.g. exposed windows are not decided based on 
	 * all available tags.
	 *
	 * Remark: as in Step4, kernels with unknown period are classified as short-period.
	 *
	 * @return -1 if the descriptors directory contains less than two descriptors, zero
	 * otherwise.
	 */
	private static final int loadNodes(String step1dir) throws IOException {
		final String DESCRIPTORS_DIR = step1dir+"/"+IO.TAGS_DIR+"/"+IO.STEP4_DIR;
		final int PREFIX_LENGTH = IO.BASIN_DESCRIPTOR_PREFIX.length();
		final int INCREMENT = 1000;  // Arbitrary
		int i, j, k, p, q, r;
		int lastComponentID, component;
		String str;
		IntervalGraph.Node node;
		BufferedReader br;
		File directory, file;
		int[] componentIDs;
		String[] files;
		
		// Allocating global data structures
		directory = new File(DESCRIPTORS_DIR);
		files=directory.list(); 
		nNodes=0; nPathKernels=0;
		for (i=0; i<files.length; i++) {
			if (files[i].length()<=PREFIX_LENGTH || !files[i].substring(0,PREFIX_LENGTH).equalsIgnoreCase(IO.BASIN_DESCRIPTOR_PREFIX)) continue;
			p=files[i].indexOf("-"); q=files[i].indexOf("-",p+1); r=files[i].indexOf("-",q+1);
			if (p>=0 && q>=0 && r>=0) nPathKernels++;
			br = new BufferedReader(new FileReader(DESCRIPTORS_DIR+"/"+files[i]));
			str=br.readLine();  // Header
			str=br.readLine();
			while (str!=null) {
				nNodes++;
				str=br.readLine();
			}
			br.close();
		}
		if (nPathKernels<2) return -1;
		pathKernelPeriodic = new byte[nPathKernels];
		pathKernelLengths = new int[nPathKernels];
		nodesArray = new IntervalGraph.Node[nNodes];
		kernel2descriptor = new int[nPathKernels][3];
		componentIDs = new int[nPathKernels];
		lastComponentID=-1;
		for (i=0; i<files.length; i++) {
			if (files[i].length()<=PREFIX_LENGTH || !files[i].substring(0,PREFIX_LENGTH).equalsIgnoreCase(IO.BASIN_DESCRIPTOR_PREFIX)) continue;
			p=files[i].indexOf("-"); q=files[i].indexOf("-",p+1);
			componentIDs[++lastComponentID]=Integer.parseInt(files[i].substring(p+1,q));
		}
		Arrays.sort(componentIDs);
		lastComponentID=0;
		for (i=1; i<nPathKernels; i++) {
			if (componentIDs[i]!=componentIDs[lastComponentID]) {
				lastComponentID++;
				componentIDs[lastComponentID]=componentIDs[i];
			}
		}
		System.err.println("loadNodes> "+nPathKernels+" total path kernels, "+(lastComponentID+1)+" components.");
		
		// Loading every Step2 component
		nNodes=0; nPathKernels=0;
		for (i=0; i<=lastComponentID; i++) {
			component=componentIDs[i];
			j=IntervalGraphStep4.loadNodes(DESCRIPTORS_DIR,IO.BASIN_DESCRIPTOR_PREFIX+"-"+component+"-",step1dir+"/"+component+".graph",false,false,true);
			if (j==-2) continue;
			if (nPathKernels!=0) {
				for (j=0; j<IntervalGraph.nNodes; j++) {
					node=IntervalGraph.nodesArray[j];
					for (k=0; k<=node.lastKernel; k++) node.kernels[k]+=nPathKernels;
					for (k=0; k<=node.lastPathWithStart; k++) node.pathsWithStart[k]+=node.pathsWithStart[k]>=0?nPathKernels:-nPathKernels;
					for (k=0; k<=node.lastPathWithEnd; k++) node.pathsWithEnd[k]+=node.pathsWithEnd[k]>=0?nPathKernels:-nPathKernels;
				}
			}
			if (nNodes+IntervalGraph.nNodes>nodesArray.length) {
				IntervalGraph.Node[] newArray = new IntervalGraph.Node[nNodes+IntervalGraph.nNodes+INCREMENT];
				System.arraycopy(nodesArray,0,newArray,0,nNodes);
				nodesArray=newArray;
			}
			System.arraycopy(IntervalGraph.nodesArray,0,nodesArray,nNodes,IntervalGraph.nNodes);
			if (nPathKernels+IntervalGraphStep3.nPathKernels>pathKernelLengths.length) {
				int[] newArray = new int[nPathKernels+IntervalGraphStep3.nPathKernels+INCREMENT];
				System.arraycopy(pathKernelLengths,0,newArray,0,nPathKernels);
				pathKernelLengths=newArray;
			}
			System.arraycopy(IntervalGraphStep3.pathKernelLengths,0,pathKernelLengths,nPathKernels,IntervalGraphStep3.nPathKernels);
			if (nPathKernels+IntervalGraphStep3.nPathKernels>pathKernelPeriodic.length) {
				byte[] newArray = new byte[nPathKernels+IntervalGraphStep3.nPathKernels+INCREMENT];
				System.arraycopy(pathKernelPeriodic,0,newArray,0,nPathKernels);
				pathKernelPeriodic=newArray;
			}
			System.arraycopy(IntervalGraphStep3.pathKernelPeriodic,0,pathKernelPeriodic,nPathKernels,IntervalGraphStep3.nPathKernels);
			if (nPathKernels+IntervalGraphStep3.nPathKernels>kernel2descriptor.length) {
				int[][] newArray = new int[nPathKernels+IntervalGraphStep3.nPathKernels+INCREMENT][0];
				System.arraycopy(kernel2descriptor,0,newArray,0,nPathKernels);
				kernel2descriptor=newArray;
			}
			for (j=0; j<IntervalGraphStep3.nPathKernels; j++) {
				kernel2descriptor[nPathKernels+j][0]=IntervalGraphStep4.kernel2descriptor[j][0];
				kernel2descriptor[nPathKernels+j][1]=IntervalGraphStep4.kernel2descriptor[j][1];
				kernel2descriptor[nPathKernels+j][2]=IntervalGraphStep4.kernel2descriptor[j][2];
			}
			nNodes+=IntervalGraph.nNodes; nPathKernels+=IntervalGraphStep3.nPathKernels;
		}
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		Arrays.sort(nodesArray,0,nNodes);		
		return 0;
	}
	
	
	/**
	 * Sets $labels[i]=1$ iff path kernel $i$ has enough exposed windows, where a window 
	 * is the union of a maximal set of interval graph nodes that can be merged according 
	 * to $IntervalGraphStep3.printKernelTags_shouldBeMerged()$, and a window is exposed 
	 * if it has just one kernel, and if it contains a long enough interval that does not
	 * belong to any other kernel.
	 *
	 * A repeat with few exposed windows occurs rarely by itself in the genome, and is
	 * considered just as a substring shared by several other repeats, or a substring that
	 * occurs multiple times inside the same repeat, or a substring induced by two or more
	 * repeats that are frequently adjacent. Note that every such substring that is long 
	 * enough creates an interval in the factorization, and likely a distinct connected 
	 * component and kernel in the interval graph.
	 * 
	 * Every surviving configuration of non-periodic kernels can be legitimate, and there 
	 * is not enough evidence to remove any of the involved repeats. In particular, high 
	 * Jaccard similarity between the surfaces of two non-periodic kernels can still be 
	 * possible and is not necessarily redundant.
	 *	
	 * Remark: an exposed window can be just a substring of the full path kernel sequence.
	 *
	 * Remark: the procedure does not use short-period kernels to remove exposed surface
	 * from other short-period kernels. I.e. an occurrence of a short-period kernel can be
	 * covered only by occurrences of non-periodic or long-period kernels. Similarly, a
	 * long-period kernel can be covered only by occurrences of non-periodic kernels.
	 * Short-period kernels intersecting one another, and long-period kernels intersecting
	 * one another, are handled separately downstream.
	 *
	 * Remark: alternatively, one could remove every repeat that is contained at least 
	 * once in another repeat, since it may be argued that it is not useful when mapping 
	 * the modules to the reads. This is better done later in the pipeline, at the level 
	 * of consensus sequences.
	 *
	 * Remark: some intervals of prefix/suffix/substring type produced by factorization 
	 * might straddle multiple periodic intervals, and a path kernel containing one such
	 * interval might have to be discarded. The procedure does not do this explicitly.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$.
	 *
	 * @param longKernelLength nonperiodic kernels of length at least this much are 
	 * considered extremely long;
	 * @return the number of marks in $labels$.
	 */
	private static final int markExposed(int longKernelLength, int minExposedSurface) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int CAPACITY = 6;  // Arbitrary
		int minExposedWindows = IO.minOccurrencesInGenome*IO.coverage;  // Arbitrary, >1.
		boolean merged;
		int i, j;
		int read, lastWindow, nExposedWindows, out;
		IntervalGraph.Node node;
		IntervalGraph.Node[] windows;
		
		// Allocating global data structures
		tracks = new int[nPathKernels][CAPACITY];
		lastTrack = new int[nPathKernels];
		Math.set(lastTrack,nPathKernels-1,-1);
		kernelsWithTrack = new int[nPathKernels];
		cTrackStart = new int[nPathKernels][CAPACITY];
		cTrackEnd = new int[nPathKernels][CAPACITY];
		lastCTrack = new int[nPathKernels];
		Math.set(lastCTrack,nPathKernels-1,-1);
		selectedKernels = new int[nPathKernels];
		counts = new int[nPathKernels][3];
		Math.set(counts,0);
		if (tmpBoolean==null || tmpBoolean.length<nPathKernels) tmpBoolean = new boolean[nPathKernels];
		Math.set(tmpBoolean,nPathKernels-1,false);
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		
		// Finding exposed windows
		nExposedWindows=0;
		windows = new IntervalGraph.Node[CAPACITY];
		tmpWindows = new IntervalGraph.Node[CAPACITY];
		for (i=0; i<tmpWindows.length; i++) tmpWindows[i] = new IntervalGraph.Node();
		node=nodesArray[0]; read=node.read;
		lastWindow=0; 
		IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
		tmpWindows[0] = new IntervalGraph.Node();
		for (i=1; i<nNodes; i++) {
			node=nodesArray[i];
			if (node.read!=read) {
				nExposedWindows+=markExposed_impl(windows,lastWindow,minExposedSurface,longKernelLength);
				read=node.read;
				lastWindow=0; 
				IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
				continue;
			}
			merged=false;
			for (j=lastWindow; j>=0; j--) {
				if (IntervalGraphStep3.printKernelTags_shouldBeMerged(node,windows[j],IDENTITY_THRESHOLD,pathKernelLengths,pathKernelPeriodic,false)) {
					IntervalGraphStep3.printKernelTags_merge(node,windows[j],false,IDENTITY_THRESHOLD);
					merged=true;
					// Merging with multiple windows is allowed
				}
			}
			if (!merged) {
				lastWindow++;
				windows=IntervalGraphStep3.printKernelTags_ensureArray(windows,lastWindow+1,CAPACITY);
				IntervalGraphStep3.printKernelTags_setArray(windows,lastWindow,node);
			}
		}
		nExposedWindows+=markExposed_impl(windows,lastWindow,minExposedSurface,longKernelLength);
		Math.set(labels,nPathKernels-1,(byte)0);
		if (nExposedWindows==0) return 0;
		if (IO.SHOW_INTERACTIVE) {
			System.err.println("nUnique,nHighQuality,nExposed");
			for (int x=0; x<nPathKernels; x++) System.err.println(x+": "+kernel2descriptor[x][0]+"-"+kernel2descriptor[x][1]+"-"+kernel2descriptor[x][2]+" :: "+counts[x][0]+","+counts[x][1]+","+counts[x][2]);
		}		
		minExposedWindows=markExposed_getThreshold(counts,minExposedWindows);
		System.err.println("markExposed> minExposedWindows="+minExposedWindows);
		out=0;
		for (i=0; i<nPathKernels; i++) {
			if (counts[i][0]>=minExposedWindows && (counts[i][1]<minExposedWindows || counts[i][2]>=minExposedWindows)) {
				labels[i]=1;
				out++;
			}
		}
		return out;
	}
	
	
	/**
     * Increments $counts[i][0]$ by the number of windows with just kernel $i$; 
	 * $counts[i][1]$ by the subset of previous windows that have high-quality sequence on
	 * each side; $counts[i][2]$ by the subset of previous windows that have a substring 
	 * of length $>=minExposedSurface$ that is not covered by other kernels.
	 *
	 * Remark: high-quality checks are waived if the kernel is periodic, since: (1) it is
	 * often the case that long tandems (both short- and long-period) are longer than a 
	 * read; (2) periodic intervals can still be filtered downstream by set cover.
	 * The same applies to very long kernels ($\geq longKernelLength$).
	 * 
	 * Remark: the procedure uses the global array $tmpBoolean$, assumes it is initialized
	 * to all zeros, and resets it before exiting. The procedure uses also global array
	 * $tmpWindows$, possibly resizing it.
	 *
	 * @param windows assumed to be sorted by $start$;
	 * @return the total number of exposed windows found.
	 */
	private static final int markExposed_impl(IntervalGraph.Node[] windows, int lastWindow, int minExposedSurface, int longKernelLength) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int QUALITY_WINDOW = IO.quantum<<2;  // Arbitrary. Works well in practice.
		boolean highQualityLeft, highQualityRight, hasNonperiodic, hasShortPeriod, hasLongPeriod;
		int i, j, p;
		int kernel, lastKernel, kernelID, surface, out;
		IntervalGraph.Node window;
		
		// Building complement tracks
		lastKernelWithTrack=buildTracks(windows,lastWindow);
		lastSelectedKernel=-1;
		for (i=0; i<=lastWindow; i++) {
			hasNonperiodic=false; hasShortPeriod=false; hasLongPeriod=false;
			lastKernel=windows[i].lastKernel;
			for (j=0; j<=lastKernel; j++) {
				switch (pathKernelPeriodic[windows[i].kernels[j]]) {
					case 0: hasNonperiodic=true; break;
					case 1: hasShortPeriod=true; break;
					case 2: hasLongPeriod=true; break;
				}
			}
			if (lastKernel>0 && hasNonperiodic) continue;			
			for (j=0; j<=lastKernel; j++) {
				kernel=windows[i].kernels[j];
				if (pathKernelPeriodic[kernel]==1 && hasLongPeriod) continue;				
				if (!tmpBoolean[kernel]) {
					tmpBoolean[kernel]=true;
					selectedKernels[++lastSelectedKernel]=kernel;
				}
			}
		}
		for (i=0; i<=lastWindow; i++) {  // Cleaning $tmpBoolean$.
			window=windows[i];
			lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) tmpBoolean[window.kernels[j]]=false;
		}		
		buildComplementTracks();
		
		// Counting windows with enough exposed surface
		out=0;
		for (i=0; i<=lastWindow; i++) {
			hasNonperiodic=false; hasShortPeriod=false; hasLongPeriod=false;
			lastKernel=windows[i].lastKernel;
			for (j=0; j<=lastKernel; j++) {
				switch (pathKernelPeriodic[windows[i].kernels[j]]) {
					case 0: hasNonperiodic=true; break;
					case 1: hasShortPeriod=true; break;
					case 2: hasLongPeriod=true; break;
				}
			}
			if (lastKernel>0 && hasNonperiodic) continue;		
			highQualityLeft = windows[i].start>IDENTITY_THRESHOLD && !Reads.hasLowQuality(windows[i].read,windows[i].start-QUALITY_WINDOW,windows[i].start-1,true);
			highQualityRight = windows[i].end<Reads.getReadLength(windows[i].read)-IDENTITY_THRESHOLD && !Reads.hasLowQuality(windows[i].read,windows[i].end+1,windows[i].end+QUALITY_WINDOW,true);
			for (j=0; j<=lastKernel; j++) {
				kernel=windows[i].kernels[j];
				if (pathKernelPeriodic[kernel]==1 && hasLongPeriod) continue;
				counts[kernel][0]++;
				if (pathKernelPeriodic[kernel]!=0 || pathKernelLengths[kernel]>=longKernelLength || (highQualityLeft && highQualityRight)) {
					counts[kernel][1]++;
					kernelID=Math.linearSearch_unsorted(selectedKernels,0,lastSelectedKernel+1,kernel,false);
					surface=exposedSurface(windows[i].start,windows[i].end,cTrackStart[kernelID],cTrackEnd[kernelID],lastCTrack[kernelID],minExposedSurface);
					if (surface!=0) {
						counts[kernel][2]++;
						out++;
					}
				}
			}
		}
		cleanTracks(windows,lastWindow);  // Cleaning tracks for the next iteration
		return out;
	}
	
	
	/**
	 * For every kernel ID $k$, merges all compatible elements of $windows$ tagged with 
	 * $k$, and stores them in $tracks[k]$. At the end of the procedure, 
	 * $kernelsWithTrack[0..X]$ contains the (not necessarily sorted) list of distinct 
	 * kernels with a window in $windows$, where $X$ is returned in output.
	 *
	 * Remark: adjacent windows with the same kernel and orientation are not merged.
	 * Straddling windows with the same kernel but different orientation are not merged.
	 *
	 * Remark: the procedure uses the global array $tmpBoolean$, assumes it is initialized
	 * to all zeros, and resets it before exiting. The procedure assumes global array 
	 * $lastTrack$ to be initialized to all -1.
	 */
	private static final int buildTracks(IntervalGraph.Node[] windows, int lastWindow) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int INCREMENT = 15;  // Arbitrary, multiple of 3.
		byte orientation;
		int i, j;
		int last, out, windowStart, windowEnd, trackStart, trackEnd, trackOrientation;
		int kernel, lastKernel;
		IntervalGraph.Node window;
		
		out=-1;
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; windowStart=window.start; windowEnd=window.end; 
			lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) {
				kernel=window.kernels[j];
				if (!tmpBoolean[kernel]) {
					tmpBoolean[kernel]=true;
					kernelsWithTrack[++out]=kernel;
				}
				orientation=window.kernelOrientations[j];
				if (lastTrack[kernel]==-1) {
					tracks[kernel][0]=windowStart;
					tracks[kernel][1]=windowEnd;
					tracks[kernel][2]=orientation;
					lastTrack[kernel]=2;
					continue;
				}
				last=lastTrack[kernel];
				trackOrientation=tracks[kernel][last];
				trackEnd=tracks[kernel][last-1];
				trackStart=tracks[kernel][last-2];
				if ( Intervals.areApproximatelyIdentical(windowStart,windowEnd,trackStart,trackEnd) ||
					 Intervals.isApproximatelyContained(windowStart,windowEnd,trackStart,trackEnd) ||
				     Intervals.isApproximatelyContained(trackStart,trackEnd,windowStart,windowEnd) ||
				     (windowStart<trackEnd-IDENTITY_THRESHOLD && orientation==trackOrientation)
				   ) tracks[kernel][last-1]=Math.max(trackEnd,windowEnd);
				else {
					if (last+3>=tracks[kernel].length) {
						int[] newTracks = new int[tracks[kernel].length+INCREMENT];
						System.arraycopy(tracks[kernel],0,newTracks,0,tracks[kernel].length);
						tracks[kernel]=newTracks;
					}
					tracks[kernel][last+1]=windowStart;
					tracks[kernel][last+2]=windowEnd;
					tracks[kernel][last+3]=orientation;
					lastTrack[kernel]+=3;
				}
			}
		}
		
		// Cleaning $tmpBoolean$.
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) tmpBoolean[window.kernels[j]]=false;
		}
		
		return out;
	}
	
	
	private static final void cleanTracks(IntervalGraph.Node[] windows, int lastWindow) {
		int i, j;
		int lastKernel;
		IntervalGraph.Node window;
		
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) lastTrack[window.kernels[j]]=-1;
		}
	}
	
	
	/**
	 * If K is the $i$-th element in $selectedKernels[0..lastSelectedKernel]$, the 
	 * procedure sets $cTrackStart[i],cTrackEnd[i]$ to the union of all the rows of 
	 * $tracks$ that correspond to an element of $kernelsWithTrack[0..
	 * lastKernelWithTrack]$ that is not K.
	 *
	 * Remark: if K is short-period, other short-period kernels are not used for building 
	 * its complement track. If K is long-period, other long-period or short-period 
	 * kernels are not used for building its complement track.
	 *
	 * Remark: to support any filtering criterion downstream, contained/overlapping 
	 * intervals in a complement track are not merged; just duplicated intervals are 
	 * collapsed.
	 *
	 * Remark: the procedure uses global array $tmpArray1$.
	 *
	 * @param tracks a window of the input tracks is represented by 3 cells.
	 */
	private static final void buildComplementTracks() {
		final int INCREMENT = 32;  // Arbitrary, multiple of 2.
		boolean shortPeriodI, longPeriodI;
		int i, j, k, p, q;
		int last, fromKernel, toKernel, length, start, end;
		
		// Merging arrays of intervals
		Math.set(lastCTrack,lastSelectedKernel,-1);
		for (i=0; i<=lastKernelWithTrack; i++) {
			fromKernel=kernelsWithTrack[i];
			shortPeriodI=pathKernelPeriodic[fromKernel]==1;
			longPeriodI=pathKernelPeriodic[fromKernel]==2;
			for (j=0; j<=lastSelectedKernel; j++) {
				toKernel=selectedKernels[j];
				if ( toKernel==fromKernel || 
					 (shortPeriodI && pathKernelPeriodic[toKernel]==1) || 
					 (longPeriodI && (pathKernelPeriodic[toKernel]==1 || pathKernelPeriodic[toKernel]==2))
				   ) continue;
				p=0; q=0; last=-1;
				while (p<=lastTrack[fromKernel] && q<=lastCTrack[j]) {
					if (tracks[fromKernel][p]>cTrackStart[j][q]) {
						if (last+2>=tmpArray1.length) {
							int[] newArray = new int[tmpArray1.length+INCREMENT];
							System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
							tmpArray1=newArray;
						}
						tmpArray1[++last]=cTrackStart[j][q];
						tmpArray1[++last]=cTrackEnd[j][q];
						q++;
						continue;
					}
					if (tracks[fromKernel][p]<cTrackStart[j][q]) {
						if (last+2>=tmpArray1.length) {
							int[] newArray = new int[tmpArray1.length+INCREMENT];
							System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
							tmpArray1=newArray;
						}
						tmpArray1[++last]=tracks[fromKernel][p++];
						tmpArray1[++last]=tracks[fromKernel][p++];
						p++;
						continue;
					}
					if (last+4>=tmpArray1.length) {
						int[] newArray = new int[tmpArray1.length+INCREMENT];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					tmpArray1[++last]=tracks[fromKernel][p++];
					tmpArray1[++last]=tracks[fromKernel][p++];
					tmpArray1[++last]=cTrackStart[j][q];
					tmpArray1[++last]=cTrackEnd[j][q];
					p++; q++;
				}
				if (p<=lastTrack[fromKernel]) {
					length=((lastTrack[fromKernel]-p+1)/3)<<1;
					if (last+length>=tmpArray1.length) {
						int[] newArray = new int[last+1+length];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					while (p<=lastTrack[fromKernel]) {
						tmpArray1[++last]=tracks[fromKernel][p++];
						tmpArray1[++last]=tracks[fromKernel][p++];
						p++;
					}
				}
				else if (q<=lastCTrack[j]) {
					length=(lastCTrack[j]-q+1)<<1;
					if (last+length>=tmpArray1.length) {
						int[] newArray = new int[last+1+length];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					while (q<=lastCTrack[j]) {
						tmpArray1[++last]=cTrackStart[j][q];
						tmpArray1[++last]=cTrackEnd[j][q];
						q++;
					}
				}
				last=((last+1)>>1)-1;
				if (last>=cTrackStart[j].length) {
					cTrackStart[j] = new int[last+INCREMENT];
					cTrackEnd[j] = new int[last+INCREMENT];
				}
				for (k=0; k<=last; k++) {
					cTrackStart[j][k]=tmpArray1[k<<1];
					cTrackEnd[j][k]=tmpArray1[(k<<1)+1];
				}
				lastCTrack[j]=last;				
			}
		}
		
		// Removing duplicated intervals
		for (i=0; i<=lastSelectedKernel; i++) {
			last=lastCTrack[i];
			for (j=0; j<=last; j++) {
				start=cTrackStart[i][j];
				if (start==-1) continue;
				end=cTrackEnd[i][j];
				for (k=j+1; k<=last; k++) {
					if (cTrackStart[i][k]!=start) break;
					if (cTrackEnd[i][k]==end) cTrackStart[i][k]=-1;
				}
			}
			k=-1;
			for (j=0; j<=last; j++) {
				if (cTrackStart[i][j]==-1) continue;
				k++;
				cTrackStart[i][k]=cTrackStart[i][j];
				cTrackEnd[i][k]=cTrackEnd[i][j];
			}
			lastCTrack[i]=k;
		}
	}
	
	
	/**
	 * Remark: the procedure uses the global array $tmpWindows$, possibly resizing it.
	 *
	 * @param trackStart,trackEnd intervals in the track can be contained/straddling;
	 * @return the sum of lengths of all sub-intervals of $[start..end]$ that are of 
	 * length at least $minIntervalLength$ and that do not belong to any interval in the
	 * track.
	 */
	private static final int exposedSurface(int start, int end, int[] trackStart, int[] trackEnd, int lastTrack, int minIntervalLength) {
		int i, j, k, h;
		int startJ, endJ, previousOrder, currentStart, currentEnd, length, surface;
		
		i=Arrays.binarySearch(trackStart,0,lastTrack+1,start);
		if (i<0) i=-i-1;
		k=-1;
		for (j=i; j<=lastTrack; j++) {
			startJ=trackStart[j]; 
			if (startJ>end) break;
			k++;
			if (k==tmpWindows.length) {
				IntervalGraph.Node[] newArray = new IntervalGraph.Node[tmpWindows.length<<1];
				System.arraycopy(tmpWindows,0,newArray,0,tmpWindows.length);
				tmpWindows=newArray;
			}
			if (tmpWindows[k]==null) tmpWindows[k] = new IntervalGraph.Node();
			tmpWindows[k].start=startJ;
			tmpWindows[k].end=Math.min(end,trackEnd[j]);
		}
		for (j=i-1; j>=0; j--) {
			endJ=trackEnd[j]; 
			if (endJ<start) continue;
			k++;
			if (k==tmpWindows.length) {
				IntervalGraph.Node[] newArray = new IntervalGraph.Node[tmpWindows.length<<1];
				System.arraycopy(tmpWindows,0,newArray,0,tmpWindows.length);
				tmpWindows=newArray;
			}
			if (tmpWindows[k]==null) tmpWindows[k] = new IntervalGraph.Node();
			tmpWindows[k].start=start;
			tmpWindows[k].end=Math.min(end,endJ);
		}
		if (k==-1) return end-start+1;
		if (k>0) {
			previousOrder=IntervalGraph.Node.order;
			IntervalGraph.Node.order=IntervalGraph.Node.START;
			Arrays.sort(tmpWindows,0,k+1);
			IntervalGraph.Node.order=previousOrder;
		}
		currentStart=tmpWindows[0].start; currentEnd=tmpWindows[0].end;
		length=currentStart-start;
		surface=length>=minIntervalLength?length:0;
		for (i=1; i<=k; i++) {
			if (tmpWindows[i].start>currentEnd) {
				length=tmpWindows[i].start-currentEnd-1;
				surface+=length>=minIntervalLength?length:0;
				currentStart=tmpWindows[i].start;
				currentEnd=tmpWindows[i].end;
			}
			else currentEnd=Math.max(currentEnd,tmpWindows[i].end);
		}
		length=end-currentEnd;
		surface+=length>=minIntervalLength?length:0;
		return surface;
	}
	
	
	/**
	 * Fits a density estimation tree on the number of exposed windows of every kernel, 
	 * and returns a value that separates the first local maximum from the rest.
	 *
	 * Remark: the procedure assumes $tmpPoints,tmpLeaves$ to be of size at least equal to
	 * the number of kernels.
	 *
	 * @param minThreshold the returned value is not smaller than this.
	 */
	private static final int markExposed_getThreshold(int[][] counts, int minThreshold) {
		final int N_CONSECUTIVE_POINTS = 2;  // Arbitrary
		final double QUANTILE = 0.75;  // Arbitrary
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;  // Arbitrary
		final int MAX_DIFFERENCE = 10;  // Arbitrary
		boolean isPeriodic;
		int i, j;
		int last, nLocalMaximumLeaves, nRuns, boundary, firstLocalMinimum, firstNormal;
		int leaf, nextLeaf, lastLeaf, lastPoint, nHigh;
		int minIntervalLength, minLocalMaxDistance, threshold;
		Leaf[] lvs;
		
		last=-1;
		for (i=0; i<nPathKernels; i++) {
			last++;
			tmpPoints[last].position=counts[i][2];
			tmpPoints[last].mass=1;
		}
		last=Points.sortAndCompact(tmpPoints,last);
		if (IO.SHOW_INTERACTIVE) {
			System.err.println("markExposed_getThreshold> points:");
			for (int x=0; x<=last; x++) System.err.println(tmpPoints[x].position+","+tmpPoints[x].mass);
		}
		threshold=minThreshold;
		if (last+1<=Points.FEW_POINTS) {
			if (tmpPoints[last].position-tmpPoints[0].position>MAX_DIFFERENCE) {
				i=Points.getRoughThreshold(tmpPoints,last,true,0,true);
				if (i!=-1) threshold=Math.max(threshold,(int)tmpPoints[i].position);
			}
		}
		else {	
			minIntervalLength=minThreshold;
			minLocalMaxDistance=2*minIntervalLength;
			if (!Points.areUniformlyDistributed(tmpPoints,0,last,true,(tmpPoints[last].position-tmpPoints[0].position)/Points.DEFAULT_NBINS)) {
				DensityEstimationTree.allocateMemory(last+1);
				nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(tmpPoints,0,last,minIntervalLength,minLocalMaxDistance,true,-1,-1,-1,false,false);
				if (nLocalMaximumLeaves>0) {
					DensityEstimationTree.markRunsOfLocalMaximumLeaves(tmpPoints);
					if (IO.SHOW_INTERACTIVE) {
						System.err.println("markExposed_getThreshold> leaves:");
						for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) System.err.println(DensityEstimationTree.leaves[x]);
					}
					leaf=Leaves.lastLeafOfRun(DensityEstimationTree.leaves,DensityEstimationTree.lastLeaf,0);
					nRuns=DensityEstimationTree.getNumberOfRuns();
					boundary=nRuns==1?(int)tmpPoints[last].position:(int)tmpPoints[DensityEstimationTree.leaves[leaf+1].firstPoint].position;
					firstLocalMinimum=Leaves.firstLocalMinimum(DensityEstimationTree.leaves,DensityEstimationTree.lastLeaf);
					firstNormal=Leaves.firstNormal(DensityEstimationTree.leaves,DensityEstimationTree.lastLeaf);
					nextLeaf=Leaves.nextNonMaximum(DensityEstimationTree.leaves,DensityEstimationTree.lastLeaf,leaf,true,boundary,firstLocalMinimum,firstNormal,false);
					if (IO.SHOW_INTERACTIVE) System.err.println("markExposed_getThreshold> firstLocalMinimum="+firstLocalMinimum+" firstNormal="+firstNormal+" leaf="+leaf+" nextLeaf="+nextLeaf);
					threshold=(int)tmpPoints[DensityEstimationTree.leaves[nextLeaf==-1?leaf:nextLeaf].lastPoint].position;
				}
				DensityEstimationTree.deallocateMemory();
			}
		}		
		return threshold;
	}
	
	
	
	
	
	
	
	
	/*
	 * Checks whether the beginning/end of a path kernel is always followed by the 
	 * beginning/end of another (possibly straddling) path kernel. This might happen, 
	 * since a cluster in the interval graph might have been erroneously oversplit into 
	 * multiple clusters, and every such cluster might produce a path kernel that is just 
	 * part of the full repeat (since inter-cluster edges are not used for assembling the 
	 * sequence of a repeat). If the interval graph is a line, we might also have 
	 * oversplit it into multiple clusters that are instead parts of the same repeat.
	 *
	 * Remark: this implicit form of assembly is useful just for non-periodic path kernels
	 * that moreover do not have a clear boundary on at least one end. The procedure 
	 * analyzes all non-periodic path kernels, for simplicity.
	 *
	 * Remark: the procedure builds a Markov chain of order one just between beginning and
	 * end of path kernels. This is not the full Markov chain of path kernels, since in
	 * general fragments not at the beginning/end might follow other fragments.
	 *
	 * Remark: the procedure just detects correlations. This is because no correlation has
	 * been observed in practice yet.
	 *
	 * Remark: path kernels marked with a zero in the global array $labels$ are not 
	 * considered.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$.
	 *
	 * @return the number of correlations found.
	 */
	private static final int markAdjacent(String outputPath) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int CAPACITY = 4;  // Arbitrary
		final int INCREMENT = 15;  // Arbitrary. Multiple of 3.
		final int MIN_CORRELATION_OCCURRENCES = IO.minOccurrencesInGenome*IO.coverage;  // Arbitrary
		final double CORRELATION_THRESHOLD = 0.8;  // Arbitrary
		final int N_IDS = nPathKernels<<1;
		boolean merged;
		int i, j, k;
		int read, lastWindow, neighbor, nFull, nPartial;
		IntervalGraph.Node node;
		BufferedWriter bw;
		IntervalGraph.Node[] windows;
		int[] lastNeighbor, highQualityCounts;
		int[][] neighbors;
		
		// Building the adjacency graph
		if (tmpBoolean==null || tmpBoolean.length<N_IDS) tmpBoolean = new boolean[N_IDS];
		Math.set(tmpBoolean,N_IDS-1,false);
		lastNeighbor = new int[N_IDS];
		Math.set(lastNeighbor,N_IDS-1,-1);
		highQualityCounts = new int[N_IDS];
		Math.set(highQualityCounts,N_IDS-1,0);
		neighbors = new int[N_IDS][INCREMENT];
		windows = new IntervalGraph.Node[CAPACITY];
		node=nodesArray[0];
		read=node.read;
		lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
		for (i=1; i<nNodes; i++) {
			node=nodesArray[i];
			if (node.read!=read) {
				markAdjacent_impl(windows,lastWindow,neighbors,lastNeighbor,highQualityCounts,tmpBoolean);
				read=node.read;
				lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
				continue;
			}
			merged=false;
			for (j=lastWindow; j>=0; j--) {  // Merging with multiple windows allowed
				if (IntervalGraphStep3.printKernelTags_shouldBeMerged(node,windows[j],IDENTITY_THRESHOLD,pathKernelLengths,pathKernelPeriodic,false)) {
					IntervalGraphStep3.printKernelTags_merge(node,windows[j],false,IDENTITY_THRESHOLD);
					merged=true;
				}
			}
			if (!merged) {
				lastWindow++;
				windows=IntervalGraphStep3.printKernelTags_ensureArray(windows,lastWindow+1,CAPACITY);
				IntervalGraphStep3.printKernelTags_setArray(windows,lastWindow,node);
			}
		}
		markAdjacent_impl(windows,lastWindow,neighbors,lastNeighbor,highQualityCounts,tmpBoolean);
		
		// Detecting correlations
		bw = new BufferedWriter(new FileWriter(outputPath));
		nFull=0; nPartial=0;
		for (i=0; i<N_IDS; i++) {
			k=i>>1;
			if (labels[k]==0 || pathKernelPeriodic[k]!=0 || pathKernelLengths[k]==0) continue;
			for (j=0; j<=lastNeighbor[i]; j+=3) {
				if (neighbors[i][j+1]>highQualityCounts[i]) {
					System.err.println("markAdjacent> ERROR in counts: i="+i+", neighbor="+neighbors[i][j]+", "+neighbors[i][j+1]+">"+highQualityCounts[i]);
					System.exit(1);
				}
				neighbor=neighbors[i][j];
				k=neighbor>>1;
				if (labels[k]==0 || pathKernelPeriodic[k]!=0 || pathKernelLengths[k]==0) continue;
				if (highQualityCounts[i]<MIN_CORRELATION_OCCURRENCES || highQualityCounts[neighbor]<MIN_CORRELATION_OCCURRENCES) continue;
				if (neighbors[i][j+1]>=highQualityCounts[i]*CORRELATION_THRESHOLD) {
					nPartial++;
					bw.write(kernel2descriptor[i>>1][0]+"-"+kernel2descriptor[i>>1][1]+"-"+kernel2descriptor[i>>1][2]+(i%2==0?"-start":"-end")+ " -> "+kernel2descriptor[neighbor>>1][0]+"-"+kernel2descriptor[neighbor>>1][1]+"-"+kernel2descriptor[neighbor>>1][2]+(neighbor%2==0?"-start":"-end")+"\n");
					for (k=0; k<lastNeighbor[neighbor]; k+=3) {
						if (neighbors[neighbor][k]==i && neighbors[neighbor][k+1]>=highQualityCounts[neighbor]*CORRELATION_THRESHOLD) {
							bw.write(kernel2descriptor[i>>1][0]+"-"+kernel2descriptor[i>>1][1]+"-"+kernel2descriptor[i>>1][2]+(i%2==0?"-start":"-end")+ " <-> "+kernel2descriptor[neighbor>>1][0]+"-"+kernel2descriptor[neighbor>>1][1]+"-"+kernel2descriptor[neighbor>>1][2]+(neighbor%2==0?"-start":"-end")+"\n");
							nFull++;
						}						
					}
				}
			}
		}
		bw.close();
		if (nFull+nPartial>0) {
			if (IO.SHOW_INTERACTIVE) markAdjacent_printGraph(neighbors,lastNeighbor,highQualityCounts,N_IDS,MIN_CORRELATION_OCCURRENCES,0.1/*minProbability for display*/);
		}
		else System.err.println("markAdjacent> No correlation found");
		return nFull>>1;
	}

	
	/**
	 * @param flags temporary space of size at least $2*nPathKernels$; the procedure 
	 * assumes that they are initialized to all zeros, and resets them to all zeros before
	 * returning;
	 * @return the number of new elements added to $neighbors$ by the procedure.
	 */
	private static final int markAdjacent_impl(IntervalGraph.Node[] windows, int lastWindow, int[][] neighbors, int[] lastNeighbor, int[] highQualityCounts, boolean[] flags) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		final int QUALITY_WINDOW = IO.quantum<<2;  // Arbitrary. Works well in practice.
		int i, j;
		int out, startI, endI, startJ, endJ, read;
		
		out=0;
		for (i=0; i<=lastWindow; i++) {
			startI=windows[i].start; endI=windows[i].end; read=windows[i].read;
			if (windows[i].lastPathWithStart!=-1 && startI>DISTANCE_THRESHOLD && !Reads.hasLowQuality(read,startI-QUALITY_WINDOW,startI-1,true)) {
				markAdjacent_incrementHighQualityCounts(windows[i],true,highQualityCounts);
				for (j=i-1; j>=0; j--) {
					startJ=windows[j].start; endJ=windows[j].end;
					if (windows[j].lastPathWithEnd==-1 || endJ<startI-DISTANCE_THRESHOLD || endJ>=endI-DISTANCE_THRESHOLD || startJ>=startI-DISTANCE_THRESHOLD) continue;
					out+=markAdjacent_addNeighbors(windows[i].pathsWithStart,windows[i].lastPathWithStart,windows[j].pathsWithEnd,windows[j].lastPathWithEnd,startI-endJ,neighbors,lastNeighbor,flags);
				}
				for (j=i-1; j>=0; j--) {
					startJ=windows[j].start; endJ=windows[j].end;
					if (windows[j].lastPathWithEnd==-1 || endJ<startI-DISTANCE_THRESHOLD || endJ>=endI-DISTANCE_THRESHOLD || startJ>=startI-DISTANCE_THRESHOLD) continue;
					markAdjacent_cleanFlags(windows[i].pathsWithStart,windows[i].lastPathWithStart,windows[j].pathsWithEnd,windows[j].lastPathWithEnd,flags);
				}
			}
			if (windows[i].lastPathWithEnd!=-1 && endI<Reads.getReadLength(read)-DISTANCE_THRESHOLD && !Reads.hasLowQuality(read,endI+1,endI+QUALITY_WINDOW,true)) {
				markAdjacent_incrementHighQualityCounts(windows[i],false,highQualityCounts);
				for (j=i+1; j<=lastWindow; j++) {
					startJ=windows[j].start; endJ=windows[j].end;
					if (startJ>endI+DISTANCE_THRESHOLD) break;
					if (windows[j].lastPathWithStart==-1 || startJ<=startI+DISTANCE_THRESHOLD || endJ<=endI+DISTANCE_THRESHOLD) continue;
					out+=markAdjacent_addNeighbors(windows[i].pathsWithEnd,windows[i].lastPathWithEnd,windows[j].pathsWithStart,windows[j].lastPathWithStart,startJ-endI,neighbors,lastNeighbor,flags);
				}
				for (j=i+1; j<=lastWindow; j++) {
					startJ=windows[j].start; endJ=windows[j].end;
					if (startJ>endI+DISTANCE_THRESHOLD) break;
					if (windows[j].lastPathWithStart==-1 || startJ<=startI+DISTANCE_THRESHOLD || endJ<=endI+DISTANCE_THRESHOLD) continue;
					markAdjacent_cleanFlags(windows[i].pathsWithEnd,windows[i].lastPathWithEnd,windows[j].pathsWithStart,windows[j].lastPathWithStart,flags);
				}
			}
		}
		return out;
	}
	
	
	/**
	 * @param startOrEnd TRUE=start.
	 */
	private static final void markAdjacent_incrementHighQualityCounts(IntervalGraph.Node window, boolean startOrEnd, int[] highQualityCounts) {
		int i;
		int kernel, kernelPrime;
		
		if (startOrEnd) {
			for (i=0; i<=window.lastPathWithStart; i++) {
				kernel=window.pathsWithStart[i];
				kernelPrime=kernel>=0?kernel:-1-kernel;
				if (labels[kernelPrime]==0 || pathKernelPeriodic[kernelPrime]!=0 || pathKernelLengths[kernelPrime]==0) continue;
				highQualityCounts[kernel>=0?(kernel<<1):1+((-1-kernel)<<1)]++;
			}
		}
		else {
			for (i=0; i<=window.lastPathWithEnd; i++) {
				kernel=window.pathsWithEnd[i];
				kernelPrime=kernel>=0?kernel:-1-kernel;
				if (labels[kernelPrime]==0 || pathKernelPeriodic[kernelPrime]!=0 || pathKernelLengths[kernelPrime]==0) continue;
				highQualityCounts[kernel>=0?(kernel<<1):1+((-1-kernel)<<1)]++;
			}
		}
	}
	
	
	/**
	 * Adds to $neighbors$ all tuples $(x,y,delta)$, where $x \in array1$, $y \in array2$,
	 * and $x,y$ are the ends of distinct kernel IDs. $delta$ is an offset: positive means
	 * that the kernels are adjacent on those ends, negative means that they overlap.
	 *
	 * Remark: the procedure uses global arrays $tmpArray{1,2,3}$.
	 *
	 * @param array2Flags of size at least $2*nPathKernels$; the procedure does not add 
	 * $y$ values that are marked in this bitvector;
	 * @return the number of new elements added to $neighbors$ by the procedure.
	 */
	private static final int markAdjacent_addNeighbors(int[] array1, int last1, int[] array2, int last2, int delta, int[][] neighbors, int[] lastNeighbor, boolean[] array2Flags) {
		int i, j;
		int l1, l2, l3, max, out, kernelI, kernelJ, idI, idJ;
		
		max=Math.max(last1+1,last2+1);
		if (max>tmpArray1.length) tmpArray1 = new int[max];
		if (max>tmpArray2.length) tmpArray2 = new int[max];
		if (max>tmpArray3.length) tmpArray3 = new int[max];
		System.arraycopy(array1,0,tmpArray1,0,last1+1);
		l1=Math.makePositive(tmpArray1,0,last1);
		System.arraycopy(array2,0,tmpArray2,0,last2+1);
		l2=Math.makePositive(tmpArray2,0,last2);
		l3=Math.setIntersection(tmpArray1,0,l1,tmpArray2,0,l2,tmpArray3,0);
		out=0;
		for (i=0; i<=last1; i++) {
			kernelI=array1[i];
			if (kernelI>=0) idI=kernelI<<1;
			else {
				kernelI=-1-kernelI;
				idI=1+(kernelI<<1);
			}
			if (labels[kernelI]==0 || pathKernelPeriodic[kernelI]!=0 || pathKernelLengths[kernelI]==0 || Arrays.binarySearch(tmpArray3,0,l3+1,kernelI)>=0) continue;
			for (j=0; j<=last2; j++) {
				kernelJ=array2[j];
				if (kernelJ>=0) idJ=kernelJ<<1;
				else {
					kernelJ=-1-kernelJ;
					idJ=1+(kernelJ<<1);
				}
				if (labels[kernelJ]==0 || pathKernelPeriodic[kernelJ]!=0 || pathKernelLengths[kernelJ]==0 || Arrays.binarySearch(tmpArray3,0,l3+1,kernelJ)>=0) continue;
				if (array2Flags[idJ]) continue;
				out+=markAdjacent_addNeighbor(idI,idJ,delta,neighbors,lastNeighbor)?1:0;
				array2Flags[idJ]=true;
			}
		}
		return out;
	}
	
	
	/**
	 * Variant of $markAdjacent_addNeighbors()$ that just resets $array2Flags$ to all 
	 * zeros.
	 */
	private static final void markAdjacent_cleanFlags(int[] array1, int last1, int[] array2, int last2, boolean[] array2Flags) {
		int i, j;
		int l1, l2, l3, kernelI, kernelJ, idI, idJ;
		
		System.arraycopy(array1,0,tmpArray1,0,last1+1);
		l1=Math.makePositive(tmpArray1,0,last1);
		System.arraycopy(array2,0,tmpArray2,0,last2+1);
		l2=Math.makePositive(tmpArray2,0,last2);
		l3=Math.setIntersection(tmpArray1,0,l1,tmpArray2,0,l2,tmpArray3,0);
		for (i=0; i<=last1; i++) {
			kernelI=array1[i];
			if (kernelI>=0) idI=kernelI<<1;
			else {
				kernelI=-1-kernelI;
				idI=1+(kernelI<<1);
			}
			if (labels[kernelI]==0 || pathKernelPeriodic[kernelI]!=0 || pathKernelLengths[kernelI]==0 || Arrays.binarySearch(tmpArray3,0,l3+1,kernelI)>=0) continue;
			for (j=0; j<=last2; j++) {
				kernelJ=array2[j];
				if (kernelJ>=0) idJ=kernelJ<<1;
				else {
					kernelJ=-1-kernelJ;
					idJ=1+(kernelJ<<1);
				}
				if (labels[kernelJ]==0 || pathKernelPeriodic[kernelJ]!=0 || pathKernelLengths[kernelJ]==0 || Arrays.binarySearch(tmpArray3,0,l3+1,kernelJ)>=0) continue;
				array2Flags[idJ]=false;
			}
		}
	}
	
	
	/**
	 * Adds $id2$ to $neighbors[id1]$.
	 * 
	 * @return TRUE if the procedure added a new element to $neighbors[id1]$.
	 */
	private static final boolean markAdjacent_addNeighbor(int id1, int id2, int delta, int[][] neighbors, int[] lastNeighbor) {
		final int INCREMENT = 30;  // Arbitrary. Multiple of 3.
		int i, last;
		
		last=lastNeighbor[id1];
		for (i=0; i<=last; i+=3) {
			if (neighbors[id1][i]==id2) {
				neighbors[id1][i+1]++;
				neighbors[id1][i+2]+=delta;
				return false;
			}
		}
		if (last+3>=neighbors[id1].length) {
			int[] newArray = new int[neighbors[id1].length+INCREMENT];
			System.arraycopy(neighbors[id1],0,newArray,0,neighbors[id1].length);
			neighbors[id1]=newArray;
		}
		neighbors[id1][last+1]=id2;
		neighbors[id1][last+2]=1;
		neighbors[id1][last+3]=delta;
		lastNeighbor[id1]+=3;
		return true;
	}
	
	
	/**
	 * Plots just nodes with at least $minCorrelationOccurrences$ high-quality occurrences
	 * and arcs with $minProbability$ probability ([0..1]). Without such constrains, every
	 * pair of nodes should be connected by two directed arcs.
	 */
	private static final void markAdjacent_printGraph(int[][] neighbors, int[] lastNeighbor, int[] highQualityCounts, int nIDs, int minCorrelationOccurrences, double minProbability) {
		final double PEN_SCALE = 10.0;
		final double MIN_PROBABILITY = 0.1;  // Arbitrary
		boolean found;
		int i, j;
		int neighbor, length, freqI, edgeSize;
		double probability;
		StringBuilder sb = new StringBuilder();
		
		// Marking nodes to be printed
		if (tmpBoolean==null || tmpBoolean.length<nIDs) tmpBoolean = new boolean[nIDs];
		Math.set(tmpBoolean,nIDs-1,false);
		for (i=0; i<nIDs; i++) {
			if (labels[i]==0 || pathKernelPeriodic[i]!=0 || pathKernelLengths[i]==0) continue;
			freqI=highQualityCounts[i];
			if (freqI<minCorrelationOccurrences) continue;
			found=false;
			for (j=0; j<=lastNeighbor[i]; j+=3) {
				probability=((double)neighbors[i][j+1])/highQualityCounts[i];
				if (probability>=minProbability) {
					found=true;
					break;
				}
			}
			if (!found) continue;
			tmpBoolean[i]=true;
			for (j=0; j<=lastNeighbor[i]; j+=3) {
				probability=((double)neighbors[i][j+1])/freqI;
				if (probability<minProbability) continue;
				tmpBoolean[neighbors[i][j]]=true;
			}
		}
		
		// Printing
		System.err.println("digraph G {");
		for (i=0; i<nIDs; i++) {
			if (!tmpBoolean[i]) continue;
			sb.delete(0,sb.length());
			sb.append(kernel2descriptor[i>>1][0]);
			sb.append("-");
			sb.append(kernel2descriptor[i>>1][1]);
			sb.append("-");
			sb.append(kernel2descriptor[i>>1][2]);
			sb.append(i%2==0?"-start":"-end");
			System.err.println(i+" [label=\""+sb.toString()+"\"];");
		}
		for (i=0; i<nIDs; i++) {
			if (!tmpBoolean[i]) continue;
			freqI=highQualityCounts[i];
			for (j=0; j<=lastNeighbor[i]; j+=3) {
				neighbor=neighbors[i][j];
				if (!tmpBoolean[neighbor]) continue;
				probability=((double)neighbors[i][j+1])/freqI;
				if (probability<minProbability) continue;
				edgeSize=Math.max(1,(int)(probability*PEN_SCALE));
				sb.delete(0,sb.length());
				sb.append(i+" -> ");
				sb.append(neighbor+" [prob=\"");
				sb.append(IO.format(probability));
				sb.append("\",num=\""+neighbors[i][j+1]);
				sb.append("\",denom=\""+highQualityCounts[i]);
				sb.append("\",offset=\""+(neighbors[i][j+2]/neighbors[i][j+1]));
				sb.append("\",penwidth=\""+edgeSize);
				sb.append("\"];");
				System.err.println(sb.toString());
			}
		}
		System.err.println("}");
	}
	
	
	/**
	 * Several short-period kernels might still share a large fraction of their surface.
	 * Since every short-period fragment is likely to align to several places of every
	 * other short-period fragment, short-period kernels with similar sequence are likely
	 * to produce very similar consensus sequences in the end. This procedure builds a 
	 * graph in which there is an edge between every pair of short-period kernels that 
	 * share a large-enough amount of sequence; for each connected component of this 
	 * graph, the procedure keeps just a smallest set of kernels that cover at least 
	 * $1-minSurfaceFraction$ fraction of the total sequence used by all kernels in the 
	 * component. The same is done for long-period kernels if $shortOrLong=false$.
	 *
	 * Remark: the procedure does not apply to non-periodic path kernels, since two path
	 * kernels might share a large fraction of their surface but e.g. one might have an 
	 * insertion WRT the other.
	 *
	 * Remark: the procedure might merge path kernels that belong to different Step2 
	 * clusters of a periodic interval graph, thus undoing Step2.
	 *
	 * Remark: since kernels with unknown period are classified as short-period by 
	 * $loadNodes()$, the procedure automatically includes those when $shortOrLong=true$.
	 *
	 * Remark: path kernels marked with a zero in the global array $labels$ are not 
	 * considered.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$.
	 *
	 * @param resetOrder the procedure permutes $nodesArray$: if this flag is TRUE, 
	 * $nodesArray$ is reset to its original order at the end;
	 * @return the number of short-period kernels removed by the procedure.
	 */
	private static final int markSameSurface(double minSurfaceFraction, boolean shortOrLong, boolean resetOrder) {
		final int MIN_INTERSECTION = Alignments.minAlignmentLength;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int CAPACITY = 16;  // Arbitrary
		boolean isPeriodic;
		int i, j, c;
		int read, lastWindow, nEdges, nEdgesPeriodic, nComponents, neighbor, out;
		int nPeriod, nOthers, component, length, lastPeriodComponent;
		double jaccardThreshold;
		IntervalGraph.Node node;
		int[] last, lastNeighbor, nOutNeighbors, components, periodComponents;
		IntervalGraph.Node[] windows;
		int[][] outNeighbors, component2kernels;
		long[][] neighbors;
		
		// Allocating global data structures
		nPeriod=0; nOthers=0;
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]!=0 && ((shortOrLong && pathKernelPeriodic[i]==1) || (!shortOrLong && pathKernelPeriodic[i]==2))) nPeriod++;
			else nOthers++;
		}
		if (nPeriod==0) return 0;
		if (tracks==null) tracks = new int[nPathKernels][CAPACITY];
		if (lastTrack==null) lastTrack = new int[nPathKernels];
		Math.set(lastTrack,nPathKernels-1,-1);
		if (kernelsWithTrack==null) kernelsWithTrack = new int[nPathKernels];
		if (tmpBoolean==null || tmpBoolean.length<nPathKernels) tmpBoolean = new boolean[nPathKernels];
		Math.set(tmpBoolean,nPathKernels-1,false);
		neighbors = new long[nPathKernels][CAPACITY];
		lastNeighbor = new int[nPathKernels];
		Math.set(lastNeighbor,nPathKernels-1,-1);
		if (kernelSurface==null || kernelSurface.length<nPathKernels) kernelSurface = new long[nPathKernels];
		Math.set(kernelSurface,nPathKernels-1,0);
		
		// Building the intersection graph
		windows = new IntervalGraph.Node[CAPACITY];
		node=nodesArray[0]; read=node.read;
		lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
		nEdges=0; nEdgesPeriodic=0;
		for (i=1; i<nNodes; i++) {
			node=nodesArray[i];
			if (node.read!=read) {
				nEdges+=markSameSurface_impl(windows,lastWindow,shortOrLong,neighbors,lastNeighbor,kernelSurface);
				read=node.read;
				lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
				continue;
			}
			lastWindow++;
			windows=IntervalGraphStep3.printKernelTags_ensureArray(windows,lastWindow+1,CAPACITY);
			IntervalGraphStep3.printKernelTags_setArray(windows,lastWindow,node);
		}
		nEdges+=markSameSurface_impl(windows,lastWindow,shortOrLong,neighbors,lastNeighbor,kernelSurface);
		if (nEdges==0) return 0;
		
		// Detecting connected components
		outNeighbors = new int[nPathKernels][CAPACITY];
		nOutNeighbors = new int[nPathKernels];
		Math.set(nOutNeighbors,nPathKernels-1,0);
		nEdges=0;
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]==0 || (shortOrLong && pathKernelPeriodic[i]!=1) || (!shortOrLong && pathKernelPeriodic[i]!=2)) continue;
			for (j=0; j<=lastNeighbor[i]; j+=4) {
				neighbor=(int)neighbors[i][j];
				if (neighbors[i][j+3]>=MIN_INTERSECTION) {
					if (nOutNeighbors[i]==outNeighbors[i].length) {
						int[] newArray = new int[outNeighbors[i].length+CAPACITY];
						System.arraycopy(outNeighbors[i],0,newArray,0,outNeighbors[i].length);
						outNeighbors[i]=newArray;
					}
					outNeighbors[i][nOutNeighbors[i]++]=neighbor;
					if (nOutNeighbors[neighbor]==outNeighbors[neighbor].length) {
						int[] newArray = new int[outNeighbors[neighbor].length+CAPACITY];
						System.arraycopy(outNeighbors[neighbor],0,newArray,0,outNeighbors[neighbor].length);
						outNeighbors[neighbor]=newArray;
					}
					outNeighbors[neighbor][nOutNeighbors[neighbor]++]=i;
					nEdges++;
				}
			}
		}
		if (nEdges==0) return 0;
		if (IO.SHOW_INTERACTIVE) markSameSurface_printConnectedComponents(outNeighbors,nOutNeighbors,neighbors,lastNeighbor);
		if (tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		if (tmpArray2.length<nPathKernels) tmpArray2 = new int[nPathKernels];
		if (tmpArray3.length<nPathKernels) tmpArray3 = new int[nPathKernels];
		Math.set(tmpArray1,nPathKernels-1,0);
		components=tmpArray2;
		nComponents=DAG.getConnectedComponents(nPathKernels,null,tmpArray1,outNeighbors,nOutNeighbors,components,tmpArray3);
		System.err.println("markSameSurface> Found "+(nComponents-nOthers)+" connected components of "+nPeriod+" surviving "+(shortOrLong?"short":"long")+"-period kernels");
		
		// Building set cover
		periodComponents = new int[nPeriod];
		Math.set(tmpBoolean,nPathKernels-1,false);
		j=-1;
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]==0 || (shortOrLong && pathKernelPeriodic[i]!=1) || (!shortOrLong && pathKernelPeriodic[i]!=2)) continue;
			component=components[i];
			if (tmpBoolean[component]) continue;
			tmpBoolean[component]=true;
			periodComponents[++j]=component;
		}
		lastPeriodComponent=j;
		if (lastPeriodComponent>0) Arrays.sort(periodComponents,0,lastPeriodComponent+1);
		component2kernels = new int[lastPeriodComponent+1][CAPACITY];
		last = new int[lastPeriodComponent+1];
		Math.set(last,lastPeriodComponent,-1);
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]==0 || (shortOrLong && pathKernelPeriodic[i]!=1) || (!shortOrLong && pathKernelPeriodic[i]!=2)) continue;
			c=Arrays.binarySearch(periodComponents,0,lastPeriodComponent+1,components[i]);
			last[c]++;
			if (last[c]+1>component2kernels[c].length) {
				int[] newArray = new int[component2kernels[c].length+CAPACITY];
				System.arraycopy(component2kernels[c],0,newArray,0,component2kernels[c].length);
				component2kernels[c]=newArray;
			}			
			component2kernels[c][last[c]]=i;
		}
		j=-1;
		for (i=0; i<=lastPeriodComponent; i++) {
			if (last[i]==0) continue;
			j++;
			periodComponents[j]=periodComponents[i];
			if (last[i]+1<component2kernels[i].length) {
				int[] newArray = new int[last[i]+1];
				System.arraycopy(component2kernels[i],0,newArray,0,last[i]+1);
				component2kernels[j]=newArray;
			}
			else component2kernels[j]=component2kernels[i];
		}
		if (j==-1) return 0;
		lastPeriodComponent=j;
		Math.set(last,lastPeriodComponent,-1);
		setCover(periodComponents,lastPeriodComponent,component2kernels,last,components,minSurfaceFraction,resetOrder);
		
		// Removing kernels that are not in the set cover
		System.err.println("markSameSurface> Kernels selected in each component:");
		for (i=0; i<=lastPeriodComponent; i++) {
			System.err.print(i+": ");
			for (j=0; j<=last[i]; j++) System.err.print(component2kernels[i][j]+",");
			System.err.println();
		}
		out=0;
		for (i=0; i<=lastPeriodComponent; i++) {
			length=component2kernels[i].length;
			for (j=last[i]+1; j<length; j++) labels[component2kernels[i][j]]=0;
			out+=length-last[i]-1;
		}
		return out;
	}
	
	
	/**
	 * Remark: the procedure uses $tmpArray1$.
	 *
	 * @param neighbors output matrix; the procedure adds edges between every distinct 
	 * pair of kernels that occur in $windows$;
	 * @param kernelSurface output array: the procedure cumulates the total surface of
	 * every kernel;
	 * @return number of new edges created by the procedure.
	 */
	private static final int markSameSurface_impl(IntervalGraph.Node[] windows, int lastWindow, boolean shortOrLong, long[][] neighbors, int[] lastNeighbor, long[] kernelSurface) {
		boolean isPeriodic, newEdge;
		int i, j;
		int kernelI, kernelJ, sumI;
		int intersectionSize, intersectionSegment, nEdges;
		
		lastKernelWithTrack=markSameSurface_buildTracks(windows,lastWindow,shortOrLong);
		if (tmpArray1==null || lastKernelWithTrack+1>tmpArray1.length) tmpArray1 = new int[lastKernelWithTrack+1];
		for (i=0; i<=lastKernelWithTrack; i++) {
			kernelI=kernelsWithTrack[i];
			if (labels[kernelI]==0 || (shortOrLong && pathKernelPeriodic[kernelI]!=1) || (!shortOrLong && pathKernelPeriodic[kernelI]!=2)) continue;
			sumI=0;
			for (j=0; j<lastTrack[kernelI]; j+=2) {
				sumI+=tracks[kernelI][j+1]-tracks[kernelI][j]+1;
			}
			tmpArray1[i]=sumI;
			kernelSurface[kernelI]+=sumI;
		}
		nEdges=0;
		for (i=0; i<=lastKernelWithTrack; i++) {
			kernelI=kernelsWithTrack[i];
			if (labels[kernelI]==0 || (shortOrLong && pathKernelPeriodic[kernelI]!=1) || (!shortOrLong && pathKernelPeriodic[kernelI]!=2)) continue;
			sumI=tmpArray1[i];
			for (j=i+1; j<=lastKernelWithTrack; j++) {
				kernelJ=kernelsWithTrack[j];
				if (labels[kernelJ]==0 || (shortOrLong && pathKernelPeriodic[kernelJ]!=1) || (!shortOrLong && pathKernelPeriodic[kernelJ]!=2)) continue;
				intersectionSize=Intervals.intersectionLength(tracks[kernelI],lastTrack[kernelI],tracks[kernelJ],lastTrack[kernelJ],1,true);
				intersectionSegment=Intervals.intersectionLength(tracks[kernelI],lastTrack[kernelI],tracks[kernelJ],lastTrack[kernelJ],1,false);
				newEdge=markSameSurface_addEdge(kernelI,kernelJ,intersectionSize,sumI+tmpArray1[j]-intersectionSize,intersectionSegment,neighbors,lastNeighbor);
				if (newEdge) nEdges++;
			}
		}
		cleanTracks(windows,lastWindow);  // For the next iteration
		return nEdges;
	}
	
	
	/**
	 * For every kernel ID $k$, merges all overlapping and adjacent elements of $windows$ 
	 * tagged with $k$, and stores them in $tracks[k]$. At the end of the procedure, 
	 * $kernelsWithTrack[0..X]$ contains the (not necessarily sorted) list of distinct 
	 * kernels with a window in $windows$, where $X$ is returned in output.
	 *
	 * Remark: the procedure uses the global array $tmpBoolean$, assumes it is initialized
	 * to all zeros, and resets it before exiting. The procedure assumes global array 
	 * $lastTrack$ to be initialized to all -1.
	 */
	private static final int markSameSurface_buildTracks(IntervalGraph.Node[] windows, int lastWindow, boolean shortOrLong) {
		final int INCREMENT = 16;  // Arbitrary, multiple of 2.
		int i, j;
		int last, out, windowStart, windowEnd, trackStart, trackEnd;
		int kernel, lastKernel;
		IntervalGraph.Node window;
		
		out=-1;
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; windowStart=window.start; windowEnd=window.end; 
			lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) {
				kernel=window.kernels[j];
				if (labels[kernel]==0 || (shortOrLong && pathKernelPeriodic[kernel]!=1) || (!shortOrLong && pathKernelPeriodic[kernel]!=2)) continue;
				if (!tmpBoolean[kernel]) {
					tmpBoolean[kernel]=true;
					kernelsWithTrack[++out]=kernel;
				}
				if (lastTrack[kernel]==-1) {
					tracks[kernel][0]=windowStart;
					tracks[kernel][1]=windowEnd;
					lastTrack[kernel]=1;
					continue;
				}
				last=lastTrack[kernel];
				trackEnd=tracks[kernel][last];
				trackStart=tracks[kernel][last-1];
				if (windowStart<=trackEnd) tracks[kernel][last]=Math.max(trackEnd,windowEnd);
				else {
					if (last+2>=tracks[kernel].length) {
						int[] newTracks = new int[tracks[kernel].length+INCREMENT];
						System.arraycopy(tracks[kernel],0,newTracks,0,tracks[kernel].length);
						tracks[kernel]=newTracks;
					}
					tracks[kernel][last+1]=windowStart;
					tracks[kernel][last+2]=windowEnd;
					lastTrack[kernel]+=2;
				}
			}
		}
		
		// Cleaning $tmpBoolean$.
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) tmpBoolean[window.kernels[j]]=false;
		}
		
		return out;
	}
	
	
	/**
	 * @return TRUE iff the procedure created a new edge.
	 */
	private static final boolean markSameSurface_addEdge(int kernelFrom, int kernelTo, int intersectionSize, int unionSize, int intersectionSegment, long[][] neighbors, int[] lastNeighbor) {
		final int INCREMENT = 16;  // Arbitrary, multiple of 4.
		boolean found;
		int i;
		final int min = Math.min(kernelFrom,kernelTo);
		final int max = Math.max(kernelFrom,kernelTo);
		int last;
		
		found=false; last=lastNeighbor[min];
		for (i=0; i<last; i+=4) {
			if (neighbors[min][i]==max) {
				found=true;
				neighbors[min][i+1]+=intersectionSize;
				neighbors[min][i+2]+=unionSize;
				neighbors[min][i+3]=Math.max((int)neighbors[min][i+3],intersectionSegment);
				break;
			}
		}
		if (found) return false;
		if (last+4>=neighbors[min].length) {
			long[] newArray = new long[neighbors[min].length+INCREMENT];
			System.arraycopy(neighbors[min],0,newArray,0,neighbors[min].length);
			neighbors[min]=newArray;
		}
		neighbors[min][last+1]=max;
		neighbors[min][last+2]=intersectionSize;
		neighbors[min][last+3]=unionSize;
		neighbors[min][last+4]=intersectionSegment;
		lastNeighbor[min]+=4;
		return true;		
	}
	
	
	/**
	 * Prints a DOT file to stderr.
	 *
	 * @param outNeighbors graph of connected components;
	 * @param neighbors graph of all overlaps from which $outNeighbors$ was built.
	 */
	private static final void markSameSurface_printConnectedComponents(int[][] outNeighbors, int[] nOutNeighbors, long[][] neighbors, int[] lastNeighbor) {
		final int PEN_SCALE = 10;
		int i, j, k;
		int neighbor, edgeSize;
		double jaccard;
		String color;
	
		System.err.println("graph G {");
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]==0 || pathKernelPeriodic[i]!=1) continue;
			switch (pathKernelPeriodic[i]) {
				case 0: color="white"; break;
				case 1: color=Colors.type2color_dot(Constants.INTERVAL_PERIODIC); break;
				case 2: color=Colors.COLOR_DOT_LONGPERIOD; break;
				default: color="";
			}
			System.err.println(i+" [label=\""+kernel2descriptor[i][0]+"-"+kernel2descriptor[i][1]+"-"+kernel2descriptor[i][2]+"\",style=\"filled\",fillcolor=\""+color+"\"];");
		}
		for (i=0; i<nPathKernels; i++) {
			if (nOutNeighbors[i]==0) continue;
			k=0;
			for (j=0; j<nOutNeighbors[i]; j++) {
				neighbor=outNeighbors[i][j];
				if (neighbor<i) continue;
				jaccard=-1;
				while (k<=lastNeighbor[i]) {
					if (neighbors[i][k]==neighbor) {
						jaccard=((double)neighbors[i][k+1])/neighbors[i][k+2];
						k+=4;
						break;
					}
					else k+=4;
				}
				edgeSize=Math.max(1,(int)(jaccard*PEN_SCALE));
				System.err.println(i+" -- "+neighbor+" [weight=\""+IO.format(jaccard)+"\",penwidth=\""+edgeSize+"\"];");
			}
		}
		System.err.println("}");
	}
	
	
	/**
	 * Let $components$ be the sorted array of all component IDs with at least two kernels
	 * each, and let $component2kernels[i]$ be the list of all distinct kernel IDs in 
	 * component $i$ (not necessarily sorted). The procedure stores in 
	 * $component2kernels[i][0..last[i]]$ the kernels chosen by the greedy set cover 
	 * algorithm that iteratively selects the kernel with largest non-covered surface, 
	 * where "surface" means the simple size of the union of all intervals of a kernel, 
	 * and "non-covered" means the surface that does not belong to a selected kernel.
	 *
	 * Remark: the procedure assumes $kernelSurface$ to have been already initialized.
	 *
	 * @param minSurfaceFraction the greedy algorithm stops when the non-covered surface
	 * of a component becomes smaller than this fraction of the total surface of the 
	 * component;
	 * @param resetOrder the procedure permutes $nodesArray$: if this flag is TRUE, 
	 * $nodesArray$ is reset to its original order at the end;
	 * @param completed temporary array, of size at least equal to $components$.
	 */
	private static final void setCover(int[] components, int lastComponent, int[][] component2kernels, int[] last, int[] kernel2component, double minSurfaceFraction, boolean resetOrder) {
		final int CAPACITY = 16;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		boolean found;
		int i, j, k;
		int tmpInt, read, component, kernel, lastKernel, lastPeriodicNode, length;
		int lastWindow, maxKernel, previousOrder, nCompleted;
		long maxSurface;
		IntervalGraph.Node node, tmpNode;
		boolean[] completed;
		long[] thresholdSurface, componentSurface, remainingSurface;
		IntervalGraph.Node[] windows;
		
		if (tmpBoolean==null || tmpBoolean.length<nPathKernels) tmpBoolean = new boolean[nPathKernels];
		Math.set(tmpBoolean,nPathKernels-1,false);
		completed = new boolean[lastComponent+1];
		Math.set(completed,lastComponent,false);		
		
		// Collecting at the beginning of $nodesArray$ all and only the nodes with a
		// kernel in $components$.
		previousOrder=IntervalGraph.Node.order;
		if (resetOrder) {
			for (i=0; i<nNodes; i++) {
				nodesArray[i].visited=nodesArray[i].nodeID;
				nodesArray[i].nodeID=i;
			}
		}
		k=-1;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			lastKernel=node.lastKernel;
			found=false;
			for (j=0; j<=lastKernel; j++) {
				component=kernel2component[node.kernels[j]];
				if (Arrays.binarySearch(components,0,lastComponent+1,component)>=0) {
					found=true;
					break;
				}
			}
			if (!found) continue;
			k++;
			if (k!=i) {
				tmpNode=nodesArray[k];
				nodesArray[k]=nodesArray[i];
				nodesArray[i]=tmpNode;
			}
		}
		lastPeriodicNode=k;
		if (lastPeriodicNode<=0) {
			Math.set(last,lastComponent,-1);
			return;
		}
		windows = new IntervalGraph.Node[CAPACITY];
		System.err.println("setCover> "+(lastPeriodicNode+1)+" interval graph nodes within a component");
		
		// Measuring the total surface of each connected component
		componentSurface = new long[lastComponent+1];
		Math.set(componentSurface,lastComponent,0);
		node=nodesArray[0]; read=node.read;
		lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
		for (i=1; i<=lastPeriodicNode; i++) {
			node=nodesArray[i];
			if (node.read!=read) {
				getComponentsSurface(windows,lastWindow,components,lastComponent,component2kernels,last,kernel2component,completed,componentSurface,(byte)0);
				read=node.read;
				lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
				continue;
			}
			lastWindow++;
			windows=IntervalGraphStep3.printKernelTags_ensureArray(windows,lastWindow+1,CAPACITY);
			IntervalGraphStep3.printKernelTags_setArray(windows,lastWindow,node);
		}
		getComponentsSurface(windows,lastWindow,components,lastComponent,component2kernels,last,kernel2component,completed,componentSurface,(byte)0);
		System.err.println("setCover> Surface of each component:");
		for (i=0; i<=lastComponent; i++) System.err.println(components[i]+": "+componentSurface[i]);		
		
		// Selecting the largest kernel from each component
		for (i=0; i<=lastComponent; i++) {
			length=component2kernels[i].length; maxSurface=0; maxKernel=-1;
			for (j=0; j<length; j++) {
				kernel=component2kernels[i][j];
				if (kernelSurface[kernel]>maxSurface) {
					maxSurface=kernelSurface[kernel];
					maxKernel=j;
				}
			}
			tmpInt=component2kernels[i][0];
			component2kernels[i][0]=component2kernels[i][maxKernel];
			component2kernels[i][maxKernel]=tmpInt;
		}
		Math.set(last,lastComponent,0);
		
		// Selecting other kernels if needed
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<=lastComponent; i++) {
				if (componentSurface[i]<kernelSurface[component2kernels[i][0]]) {
					System.err.println("setCover> ERROR in component "+i+"-th = "+components[i]+". componentSurface="+componentSurface[i]+" < kernelSurface="+kernelSurface[component2kernels[i][0]]+", kernel "+component2kernels[i][0]);
					System.exit(1);
				}				
			}
		}
		remainingSurface = new long[lastComponent+1];
		for (i=0; i<=lastComponent; i++) remainingSurface[i]=componentSurface[i]-kernelSurface[component2kernels[i][0]];
		thresholdSurface = new long[lastComponent+1];
		for (i=0; i<=lastComponent; i++) thresholdSurface[i]=(long)(componentSurface[i]*minSurfaceFraction);
		nCompleted=0;
		for (i=0; i<=lastComponent; i++) {
			if (remainingSurface[i]<=thresholdSurface[i]) {
				completed[i]=true;
				nCompleted++;
			}
		}
		if (IO.SHOW_INTERACTIVE) {
			System.err.println("setCover> Surface of a largest kernel in each component:  ("+nCompleted+" completed)");
			for (i=0; i<=lastComponent; i++) System.err.println(components[i]+": "+component2kernels[i][0]+"=>"+kernelSurface[component2kernels[i][0]]+" : "+completed[i]);
		}
		k=0;
		while (nCompleted<lastComponent+1) {
			System.err.println("setCover> Iteration "+k+" started: completed "+nCompleted+" out of "+(lastComponent+1));
			if (IO.SHOW_INTERACTIVE) {
				System.err.println("setCover> remainingSurface:");
				for (int x=0; x<=lastComponent; x++) System.err.println(components[x]+": "+remainingSurface[x]);
				System.err.println("setCover> kernels per component:");
				for (int x=0; x<=lastComponent; x++) {
					System.err.println(components[x]+": ");
					for (int y=0; y<=last[x]; y++) System.err.print(component2kernels[x][y]+",");
					System.err.print(" | ");
					for (int y=last[x]+1; y<component2kernels[x].length; y++) System.err.print(component2kernels[x][y]+",");
					System.err.println();
				}
			}
					
			// Recomputing the surface of all active kernels
			for (i=0; i<=lastComponent; i++) {
				if (completed[i]) continue;
				length=component2kernels[i].length;
				for (j=last[i]+1; j<length; j++) kernelSurface[component2kernels[i][j]]=0;
			}
			node=nodesArray[0]; read=node.read;
			lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
			for (i=1; i<=lastPeriodicNode; i++) {
				node=nodesArray[i];
				if (node.read!=read) {
					setCover_impl(windows,lastWindow,components,lastComponent,component2kernels,last,kernel2component,completed);
					read=node.read;
					lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(windows,0,node);
					continue;
				}
				lastWindow++;
				windows=IntervalGraphStep3.printKernelTags_ensureArray(windows,lastWindow+1,CAPACITY);
				IntervalGraphStep3.printKernelTags_setArray(windows,lastWindow,node);
			}
			setCover_impl(windows,lastWindow,components,lastComponent,component2kernels,last,kernel2component,completed);
			
			// Selecting the kernel with largest updated surface
			for (i=0; i<=lastComponent; i++) {
				if (completed[i]) continue;
				length=component2kernels[i].length; maxSurface=0; maxKernel=-1;
				for (j=last[i]+1; j<length; j++) {
					kernel=component2kernels[i][j];
					if (kernelSurface[kernel]>maxSurface) {
						maxSurface=kernelSurface[kernel];
						maxKernel=j;
					}
				}
				if (IO.CONSISTENCY_CHECKS && maxSurface>remainingSurface[i]) {
					System.err.println("setCover> ERROR in component "+i+"-th = "+components[i]+" length="+length+" last="+last[i]+" surface of next kernel: "+maxSurface+" > remaining surface: "+remainingSurface[i]+" next kernel="+component2kernels[i][maxKernel]);
					System.exit(1);
				}
				last[i]++;
				tmpInt=component2kernels[i][last[i]];
				component2kernels[i][last[i]]=component2kernels[i][maxKernel];
				component2kernels[i][maxKernel]=tmpInt;
				remainingSurface[i]-=maxSurface;
				if (remainingSurface[i]<=thresholdSurface[i]) {
					completed[i]=true;
					nCompleted++;
				}
			}
			
			k++;
		}
		
		// Restoring the original order in $nodesArray$.
		if (resetOrder) {
			IntervalGraph.Node.order=IntervalGraph.Node.NODE_ID;
			Arrays.sort(nodesArray,0,nNodes);
			for (i=0; i<nNodes; i++) nodesArray[i].nodeID=nodesArray[i].visited;
			IntervalGraph.Node.order=previousOrder;
		}
	}
	
	
	/**
	 * Adds to $componentSurface[i]$ the sum of lengths of the union of all tracks of
	 * kernels in $windows$ that belong to the $i$-th element of $components$.
	 *
	 * @param mode 0=all kernels; 1=only inactive kernels; 2=only active kernels;
	 * @param completed assumed to be all zeros.
	 */
	private static final void getComponentsSurface(IntervalGraph.Node[] windows, int lastWindow, int[] components, int lastComponent, int[][] component2kernels, int[] last, int[] kernel2component, boolean[] completed, long[] componentSurface, byte mode) {
		int i, j, sum;
		
		lastKernelWithTrack=setCover_buildTracks(windows,lastWindow,components,lastComponent,kernel2component,completed);
		getComponentsSurface_buildTracks(components,lastComponent,component2kernels,last,kernel2component,mode);		
		for (i=0; i<=lastComponent; i++) {
			sum=0;
			for (j=0; j<lastCTrack[i]; j+=2) sum+=cTrackStart[i][j+1]-cTrackStart[i][j]+1;
			componentSurface[i]+=sum;
		}
		cleanTracks(windows,lastWindow);  // For the next iteration
	}

	
	/**
	 * Assume that $kernelsWithTrack$ contains all kernel IDs in some component. The 
	 * procedure stores in $cTrack*[i]$ the union of all tracks of kernels in the $i$-th 
	 * element of $components$. Straddling/contained intervals from different tracks are 
	 * merged.
	 *
	 * Remark: the procedure uses $tmpArray1$.
	 *
	 * Remark: the procedure stores both start and end positions in $cTrackStart$, and it
	 * does not use $cTrackEnd$; $lastCTrack$ stores the last cell used in $cTrackStart$.
	 *
	 * @param mode 0=all kernels; 1=only inactive kernels; 2=only active kernels.
	 */
	private static final void getComponentsSurface_buildTracks(int[] components, int lastComponent, int[][] component2kernels, int[] last, int[] kernel2component, byte mode) {
		final int INCREMENT = 32;  // Arbitrary, multiple of 2.
		int i, j, k, p, q;
		int last1, fromKernel, toKernel, length, start, end, component;
		
		Math.set(lastCTrack,lastComponent,-1);
		for (i=0; i<=lastKernelWithTrack; i++) {
			fromKernel=kernelsWithTrack[i];
			component=kernel2component[fromKernel];
			j=Arrays.binarySearch(components,0,lastComponent+1,component);
			if (mode!=0) {
				k=Math.linearSearch_unsorted(component2kernels[j],0,component2kernels[j].length,fromKernel,false);
				if (k<0) {
					System.err.println("getComponentsSurface_buildTracks> ERROR: kernel "+fromKernel+" does not belong to component "+component+" j="+j);
					System.exit(1);
				}
				if ((mode==1 && k>last[j]) || (mode==2 && k<=last[j])) continue;
			}
			
			// Merging arrays
			p=0; q=0; last1=-1;
			while (p<=lastTrack[fromKernel] && q<lastCTrack[j]) {
				if (tracks[fromKernel][p]>cTrackStart[j][q]) {
					if (last1+2>=tmpArray1.length) {
						int[] newArray = new int[tmpArray1.length+INCREMENT];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					tmpArray1[++last1]=cTrackStart[j][q++];
					tmpArray1[++last1]=cTrackStart[j][q++];
					continue;
				}
				if (tracks[fromKernel][p]<cTrackStart[j][q]) {
					if (last1+2>=tmpArray1.length) {
						int[] newArray = new int[tmpArray1.length+INCREMENT];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					tmpArray1[++last1]=tracks[fromKernel][p++];
					tmpArray1[++last1]=tracks[fromKernel][p++];
					continue;
				}
				if (last1+4>=tmpArray1.length) {
					int[] newArray = new int[tmpArray1.length+INCREMENT];
					System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
					tmpArray1=newArray;
				}
				tmpArray1[++last1]=tracks[fromKernel][p++];
				tmpArray1[++last1]=tracks[fromKernel][p++];
				tmpArray1[++last1]=cTrackStart[j][q++];
				tmpArray1[++last1]=cTrackStart[j][q++];
			}
			if (p<=lastTrack[fromKernel]) {
				length=lastTrack[fromKernel]-p+1;
				if (last1+length>=tmpArray1.length) {
					int[] newArray = new int[last1+1+length];
					System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
					tmpArray1=newArray;
				}
				while (p<=lastTrack[fromKernel]) {
					tmpArray1[++last1]=tracks[fromKernel][p++];
					tmpArray1[++last1]=tracks[fromKernel][p++];
				}
			}
			else if (q<lastCTrack[j]) {
				length=lastCTrack[j]-q+1;
				if (last1+length>=tmpArray1.length) {
					int[] newArray = new int[last1+1+length];
					System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
					tmpArray1=newArray;
				}
				while (q<=lastCTrack[j]) {
					tmpArray1[++last1]=cTrackStart[j][q++];
					tmpArray1[++last1]=cTrackStart[j][q++];
				}
			}
			
			// Merging intervals
			q=-2; start=tmpArray1[0]; end=tmpArray1[1];
			for (p=2; p<last1; p+=2) {
				if (tmpArray1[p]<=end) end=Math.max(end,tmpArray1[p+1]);
				else {
					q+=2;
					tmpArray1[q]=start; tmpArray1[q+1]=end;
					start=tmpArray1[p]; end=tmpArray1[p+1];
				}
			}
			q+=2;
			tmpArray1[q]=start; tmpArray1[q+1]=end;
			last1=q+1;
			
			// Copying to $cTrackStart$.
			if (last1>=cTrackStart[j].length) cTrackStart[j] = new int[last1+INCREMENT];
			System.arraycopy(tmpArray1,0,cTrackStart[j],0,last1+1);
			lastCTrack[j]=last1;				
		}
	}
	
	
	/**
	 * Let K be an active kernel of component C. The procedure adds to $kernelSurface[K]$ 
	 * the total surface of K in $windows$ that is not covered by inactive kernels of C.
	 */
	private static final void setCover_impl(IntervalGraph.Node[] windows, int lastWindow, int[] components, int lastComponent, int[][] component2kernels, int[] last, int[] kernel2component, boolean[] completed) {
		boolean found;
		int i, j, c;
		int lastKernel, kernel, sum, intersectionLength;
		IntervalGraph.Node window;
			
		// Aborting if no window belongs to an active kernel
		found=false;
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) {
				kernel=window.kernels[j];
				c=Arrays.binarySearch(components,0,lastComponent+1,kernel2component[kernel]);
				if (c>=0 && !completed[c] && Math.linearSearch_unsorted(component2kernels[c],last[c]+1,component2kernels[c].length,kernel,false)>=0) {
					found=true;
					break;
				}
			}
		}
		if (!found) return;
		
		// Computing non-covered surface of active kernels
		lastKernelWithTrack=setCover_buildTracks(windows,lastWindow,components,lastComponent,kernel2component,completed);
		lastSelectedKernel=-1;
		for (i=0; i<=lastKernelWithTrack; i++) {
			kernel=kernelsWithTrack[i];
			c=Arrays.binarySearch(components,0,lastComponent+1,kernel2component[kernel]);
			if (Math.linearSearch_unsorted(component2kernels[c],last[c]+1,component2kernels[c].length,kernel,false)>=0) {
				selectedKernels[++lastSelectedKernel]=kernel;
				kernelsWithTrack[i]=-1-kernelsWithTrack[i];
			}
		}
		setCover_buildComplementTracks(kernel2component);
		for (i=0; i<=lastSelectedKernel; i++) {
			kernel=selectedKernels[i]; sum=0;
			for (j=0; j<lastTrack[kernel]; j+=2) sum+=tracks[kernel][j+1]-tracks[kernel][j]+1;
			intersectionLength=Intervals.intersectionLength(tracks[kernel],lastTrack[kernel],cTrackStart[i],lastCTrack[i],1,true);
			kernelSurface[kernel]+=sum-intersectionLength;			
		}
		cleanTracks(windows,lastWindow);  // For the next iteration
	}
	
	
	/**
	 * Computes the tracks in $windows$ of all kernels (both active and inactive) in 
	 * non-completed components.
	 *
	 * Remark: the procedure uses global array $tmpBoolean$, which is assumed to be all
	 * zeros at the beginning and is reset at the end.
	 *
	 * @param components NULL: tracks are computed for every kernel.
	 */
	private static final int setCover_buildTracks(IntervalGraph.Node[] windows, int lastWindow, int[] components, int lastComponent, int[] kernel2component, boolean[] completed) {
		final int INCREMENT = 16;  // Arbitrary, multiple of 2.
		int i, j, c;
		int last, out, windowStart, windowEnd, trackStart, trackEnd;
		int kernel, lastKernel, component;
		IntervalGraph.Node window;

		out=-1;
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; windowStart=window.start; windowEnd=window.end; 
			lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) {
				kernel=window.kernels[j];
				if (components!=null) {
					c=Arrays.binarySearch(components,0,lastComponent+1,kernel2component[kernel]);
					if (c<0 || completed[c]) continue;
				}
				if (!tmpBoolean[kernel]) {
					tmpBoolean[kernel]=true;
					kernelsWithTrack[++out]=kernel;
				}
				if (lastTrack[kernel]==-1) {
					tracks[kernel][0]=windowStart;
					tracks[kernel][1]=windowEnd;
					lastTrack[kernel]=1;
					continue;
				}
				last=lastTrack[kernel];
				trackEnd=tracks[kernel][last]; trackStart=tracks[kernel][last-1];
				if (windowStart<=trackEnd) tracks[kernel][last]=Math.max(trackEnd,windowEnd);
				else {
					if (last+2>=tracks[kernel].length) {
						int[] newTracks = new int[tracks[kernel].length+INCREMENT];
						System.arraycopy(tracks[kernel],0,newTracks,0,tracks[kernel].length);
						tracks[kernel]=newTracks;
					}
					tracks[kernel][last+1]=windowStart;
					tracks[kernel][last+2]=windowEnd;
					lastTrack[kernel]+=2;
				}
			}
		}
		
		// Cleaning $tmpBoolean$.
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) tmpBoolean[window.kernels[j]]=false;
		}		
		return out;
	}
	
	
	/**
	 * Assume that $kernelsWithTrack$ contains all kernels (both active and inactive) in 
	 * non-completed components, that inactive kernels are negative in $kernelsWithTrack$,
	 * and that all and only the active kernels are in $selectedKernels$. The procedure
	 * stores in $cTrackStart[i]$ the union of all tracks of inactive kernels in the
	 * component of active kernel $selectedKernels[i]$.
	 *
	 * Remark: the procedure stores both start and end positions in $cTrackStart$, and it
	 * does not use $cTrackEnd$. $lastCTrack$ is the last used cell in $cTrackStart$.
	 */
	private static final void setCover_buildComplementTracks(int[] kernel2component) {
		final int INCREMENT = 32;  // Arbitrary, multiple of 2.
		int i, j, k, p, q;
		int last, fromKernel, toKernel, length, start, end, component;
		
		Math.set(lastCTrack,lastSelectedKernel,-1);
		for (i=0; i<=lastKernelWithTrack; i++) {
			fromKernel=kernelsWithTrack[i];
			if (fromKernel<0) continue;
			component=kernel2component[fromKernel];			
			for (j=0; j<=lastSelectedKernel; j++) {
				toKernel=selectedKernels[j];
				if (kernel2component[toKernel]!=component) continue;
				
				// Merging arrays
				p=0; q=0; last=-1;
				while (p<=lastTrack[fromKernel] && q<lastCTrack[j]) {
					if (tracks[fromKernel][p]>cTrackStart[j][q]) {
						if (last+2>=tmpArray1.length) {
							int[] newArray = new int[tmpArray1.length+INCREMENT];
							System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
							tmpArray1=newArray;
						}
						tmpArray1[++last]=cTrackStart[j][q++];
						tmpArray1[++last]=cTrackStart[j][q++];
						continue;
					}
					if (tracks[fromKernel][p]<cTrackStart[j][q]) {
						if (last+2>=tmpArray1.length) {
							int[] newArray = new int[tmpArray1.length+INCREMENT];
							System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
							tmpArray1=newArray;
						}
						tmpArray1[++last]=tracks[fromKernel][p++];
						tmpArray1[++last]=tracks[fromKernel][p++];
						continue;
					}
					if (last+4>=tmpArray1.length) {
						int[] newArray = new int[tmpArray1.length+INCREMENT];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					tmpArray1[++last]=tracks[fromKernel][p++];
					tmpArray1[++last]=tracks[fromKernel][p++];
					tmpArray1[++last]=cTrackStart[j][q++];
					tmpArray1[++last]=cTrackStart[j][q++];
				}
				if (p<=lastTrack[fromKernel]) {
					length=lastTrack[fromKernel]-p+1;
					if (last+length>=tmpArray1.length) {
						int[] newArray = new int[last+1+length];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					while (p<=lastTrack[fromKernel]) {
						tmpArray1[++last]=tracks[fromKernel][p++];
						tmpArray1[++last]=tracks[fromKernel][p++];
					}
				}
				else if (q<lastCTrack[j]) {
					length=lastCTrack[j]-q+1;
					if (last+length>=tmpArray1.length) {
						int[] newArray = new int[last+1+length];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					while (q<=lastCTrack[j]) {
						tmpArray1[++last]=cTrackStart[j][q++];
						tmpArray1[++last]=cTrackStart[j][q++];
					}
				}
				
				// Merging intervals
				q=-2; start=tmpArray1[0]; end=tmpArray1[1];
				for (p=2; p<last; p+=2) {
					if (tmpArray1[p]<=end) end=Math.max(end,tmpArray1[p+1]);
					else {
						q+=2;
						tmpArray1[q]=start; tmpArray1[q+1]=end;
						start=tmpArray1[p]; end=tmpArray1[p+1];
					}
				}
				q+=2;
				tmpArray1[q]=start; tmpArray1[q+1]=end;
				last=q+1;
				
				// Copying to $cTrackStart$.
				if (last>=cTrackStart[j].length) cTrackStart[j] = new int[last+INCREMENT];
				System.arraycopy(tmpArray1,0,cTrackStart[j],0,last+1);
				lastCTrack[j]=last;				
			}
		}
	}
	
	
	/**
	 * It might happen that some kernels are marked as non-exposed by $markExposed()$
	 * because they cover one another, so deleting all of them would make some repeat 
	 * occurrences not belong to any kernel in the output. Let N and E be the union of all
	 * the occurrences of non-exposed and of exposed kernels, respectively. The procedure 
	 * iteratively picks a non-exposed kernel that covers the largest quantity of the 
	 * current $N \setminus E$ among all remaining kernels (until such covered surface 
	 * becomes too small), and marks this kernel in $labels$.
	 * 
	 * @param resetOrder the procedure permutes $nodesArray$: if this flag is TRUE, 
	 * $nodesArray$ is reset to its original order at the end.
	 */
	private static final void keepNonexposed(boolean resetOrder) {
		final int CAPACITY = 16;  // Arbitrary
		final int MIN_INTERVAL_LENGTH = Alignments.minAlignmentLength;
		final int MIN_EXTRA_SURFACE=MIN_INTERVAL_LENGTH*IO.minOccurrencesInGenome;
		boolean found;
		int i, j, k;
		int nonExposedKernels, lastSurface, previousOrder, lastKernel;
		int lastNonexposedNode, selectedKernel;
		long minNonexposedSurface, surfaceSize, previousSurfaceSize, inters;
		IntervalGraph.Node node, tmpNode;
		int[] tmpArray;
		long[] intersection;
		
		// Computing the surface of non-exposed kernels that is not covered by exposed
		// kernels.
		nonExposedKernels=0;
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]==0) nonExposedKernels++;
		}
		System.err.println("keepNonexposed> "+nonExposedKernels+" non-exposed path kernels");
		if (nonExposedKernels<2) return;
		lastSurface=nonexposedMinusExposed(MIN_INTERVAL_LENGTH);
		surfaceSize=0;
		for (i=0; i<lastSurface; i+=3) surfaceSize+=tmpArray1[i+2]-tmpArray1[i+1]+1;
		System.err.println("keepNonexposed> non-exposed surface not covered by exposed: "+surfaceSize+" intervals="+(lastSurface+1));
		if (surfaceSize<MIN_EXTRA_SURFACE) return;
		
		// Collecting at the beginning of $nodesArray$ all and only the nodes with a
		// non-exposed kernel.
		previousOrder=IntervalGraph.Node.order;
		if (resetOrder) {
			for (i=0; i<nNodes; i++) {
				nodesArray[i].visited=nodesArray[i].nodeID;
				nodesArray[i].nodeID=i;
			}			
		}
		k=-1;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			lastKernel=node.lastKernel;
			found=false;
			for (j=0; j<=lastKernel; j++) {
				if (labels[node.kernels[j]]==0) {
					found=true;
					break;
				}
			}
			if (!found) continue;
			k++;
			if (k!=i) {
				tmpNode=nodesArray[k];
				nodesArray[k]=nodesArray[i];
				nodesArray[i]=tmpNode;
			}
		}
		lastNonexposedNode=k;
		if (lastNonexposedNode==-1) return;
		System.err.println("keepNonexposed> "+(lastNonexposedNode+1)+" interval graph nodes with non-exposed kernel");
		
		// Allocating reused data structures
		if (tracks==null || tracks.length<nPathKernels) tracks = new int[nPathKernels][CAPACITY];
		if (lastTrack==null || lastTrack.length<nPathKernels) lastTrack = new int[nPathKernels];
		Math.set(lastTrack,nPathKernels-1,-1);
		if (kernelsWithTrack==null || kernelsWithTrack.length<nPathKernels) kernelsWithTrack = new int[nPathKernels];
		tmpWindows = new IntervalGraph.Node[CAPACITY];
		if (tmpArray2==null || tmpArray2.length<tmpArray1.length) tmpArray2 = new int[tmpArray1.length];
		if (tmpArray3==null || tmpArray3.length<tmpArray1.length) tmpArray3 = new int[tmpArray1.length];
		if (tmpBoolean==null || tmpBoolean.length<nPathKernels) tmpBoolean = new boolean[nPathKernels];
		Math.set(tmpBoolean,nPathKernels-1,false);
		intersection = new long[nPathKernels];
		Math.set(intersection,nPathKernels-1,0);
		
		// Greedy set cover
		while (true) {
			selectedKernel=keepNonexposed_maxIntersection(lastNonexposedNode,tmpArray1,lastSurface,MIN_INTERVAL_LENGTH,intersection);
			inters=intersection[selectedKernel];
			if (IO.CONSISTENCY_CHECKS) {
				if (inters>surfaceSize) {
					System.err.println("keepNonexposed> ERROR: the intersection is bigger than the remaining surface?! "+inters+" > "+surfaceSize);
					System.exit(1);
				}
				if (surfaceSize>0 && inters<=0) {
					System.err.println("keepNonexposed> ERROR: the intersection is zero but the remaining surface is nonzero: "+inters+" :: "+surfaceSize);
					System.exit(1);
				}
			}
			if (inters<MIN_EXTRA_SURFACE) break;
			lastSurface=keepNonexposed_subtract(selectedKernel,tmpArray1,lastSurface,lastNonexposedNode);
			tmpArray=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpArray;
			labels[selectedKernel]=1;
			previousSurfaceSize=surfaceSize;
			surfaceSize=0;
			for (i=0; i<lastSurface; i+=3) surfaceSize+=tmpArray1[i+2]-tmpArray1[i+1]+1;
			if (IO.CONSISTENCY_CHECKS && surfaceSize>=previousSurfaceSize) {
				System.err.println("keepNonexposed> ERROR: after subtracting kernel "+selectedKernel+", the surface is the same or bigger?! "+surfaceSize+" >= "+previousSurfaceSize);
				System.exit(1);
			}
			System.err.println("keepNonexposed> Selected kernel "+selectedKernel+" with intersection "+inters+" ("+MIN_EXTRA_SURFACE+"). Remaining surface: "+surfaceSize);
		}
		
		// Restoring the original order in $nodesArray$, if needed.
		if (resetOrder) {
			IntervalGraph.Node.order=IntervalGraph.Node.NODE_ID;
			Arrays.sort(nodesArray,0,nNodes);
			for (i=0; i<nNodes; i++) nodesArray[i].nodeID=nodesArray[i].visited;
			IntervalGraph.Node.order=previousOrder;
		}
	}
	
	
	/**
	 * Stores in global variable $tmpArray1[0..X]$ the sorted list of disjoint intervals, 
	 * of length $>=minIntervalLength$ each, that result from computing $NE \setminus E$, 
	 * where E (respectively, NE) is the set of disjoint intervals that result from 
	 * merging all overlapping intervals that belong to any kernel $i$ with $labels[i]=1$ 
	 * (respectively, $labels[i]=0$). $X$ is returned in output; every interval is a 
	 * triplet $(read,start,end)$; $tmpArray1$ is sorted by $read,start$.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$.
	 */
	private static final int nonexposedMinusExposed(int minIntervalLength) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int CAPACITY = 1000;  // Arbitrary
		boolean hasExposed, hasNonexposed;
		int i, j;
		int lastTrack, lastExposedTrack, lastKernel, lastOut, length;
		int currentRead, currentStart, currentExposedStart, currentEnd, currentExposedEnd;
		IntervalGraph.Node node;
		int[] track, exposedTrack;
		
		track = new int[(CAPACITY)<<1];
		exposedTrack = new int[(CAPACITY)<<1];
		if (tmpArray1==null || tmpArray1.length<3*CAPACITY) tmpArray1 = new int[3*CAPACITY];
		lastOut=-1; lastTrack=-1; currentStart=-1; currentEnd=-1; currentRead=-1;
		lastExposedTrack=-1; currentExposedStart=-1; currentExposedEnd=-1;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			lastKernel=node.lastKernel;
			if (lastKernel==-1) continue;
			hasExposed=false; hasNonexposed=false;
			for (j=0; j<=lastKernel; j++) {
				if (labels[node.kernels[j]]==1) hasExposed=true;
				else hasNonexposed=true;
				if (hasExposed && hasNonexposed) break;
			}
			if (node.read!=currentRead) {
				if (currentRead!=-1) {
					length=lastOut+1+((lastTrack+lastExposedTrack+2)>>1)*3;
					if (tmpArray1.length<length) {
						int[] newArray = new int[length+CAPACITY*3];
						System.arraycopy(tmpArray1,0,newArray,0,lastOut+1);
						tmpArray1=newArray;
					}
					lastOut=Intervals.difference(track,lastTrack,exposedTrack,lastExposedTrack,tmpArray1,lastOut+1,currentRead,minIntervalLength);
				}
				currentRead=node.read; 
				if (hasExposed) { currentExposedStart=node.start; currentExposedEnd=node.end; }
				else { currentExposedStart=-1; currentExposedEnd=-1; }
				if (hasNonexposed) { currentStart=node.start; currentEnd=node.end; }
				else { currentStart=-1; currentEnd=-1; }
				lastExposedTrack=-1; lastTrack=-1;
				continue;
			}
			if (hasExposed) {
				if (currentExposedStart==-1) { currentExposedStart=node.start; currentExposedEnd=node.end; }
				else {
					if (node.start<=currentExposedEnd+IDENTITY_THRESHOLD) currentExposedEnd=Math.max(currentExposedEnd,node.end);
					else {
						if (lastExposedTrack==exposedTrack.length) {
							int[] newArray = new int[exposedTrack.length+((CAPACITY)<<1)];
							System.arraycopy(exposedTrack,0,newArray,0,exposedTrack.length);
							exposedTrack=newArray;
						}
						exposedTrack[++lastExposedTrack]=currentExposedStart;
						exposedTrack[++lastExposedTrack]=currentExposedEnd;
						currentExposedStart=node.start; currentExposedEnd=node.end;
					}
				}
			}
			if (hasNonexposed) {
				if (currentStart==-1) { currentStart=node.start; currentEnd=node.end; }
				else {
					if (node.start<=currentEnd+IDENTITY_THRESHOLD) currentEnd=Math.max(currentEnd,node.end);
					else {
						if (lastTrack==track.length) {
							int[] newArray = new int[track.length+((CAPACITY)<<1)];
							System.arraycopy(track,0,newArray,0,track.length);
							track=newArray;
						}
						track[++lastTrack]=currentStart;
						track[++lastTrack]=currentEnd;
						currentStart=node.start; currentEnd=node.end;
					}
				}
			}
		}
		length=lastOut+1+((lastTrack+lastExposedTrack+2)>>1)*3;
		if (tmpArray1.length<length) {
			int[] newArray = new int[length];
			System.arraycopy(tmpArray1,0,newArray,0,lastOut+1);
			tmpArray1=newArray;
		}
		lastOut=Intervals.difference(track,lastTrack,exposedTrack,lastExposedTrack,tmpArray1,lastOut+1,currentRead,minIntervalLength);
		return lastOut;
	}
	
	
	/**
	 * Remark: the procedure uses $tmpArray2,tmpWindows,tmpLong$ as temporary space.
	 *
	 * @param lastNode the procedure uses only $nodesArray[0..lastNode]$, which are 
	 * assumed to contain all and only the nodes with a non-exposed kernel;
	 * @param minIntervalLength only intervals of the intersection that are 
	 * $>=minIntervalLength$ contribute to the intersection length;
	 * @param intersection temporary space, of size at least equal to $nPathKernels$;
	 * @return the ID of a non-exposed kernel ($labels[i]=0$) with largest intersection
	 * with $surface$.
	 */
	private static final int keepNonexposed_maxIntersection(int lastNode, int[] surface, int lastSurface, int minIntervalLength, long[] intersection) {
		final int CAPACITY = 16;  // Arbitrary
		int i, j;
		int read, lastWindow;
		IntervalGraph.Node node;
		
		Math.set(intersection,nPathKernels-1,0);
		tmpLong[0]=-1; tmpLong[1]=0;
		node=nodesArray[0]; read=node.read;
		lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(tmpWindows,0,node); j=0;
		for (i=1; i<=lastNode; i++) {
			node=nodesArray[i];
			if (node.read!=read) {
				while (j<lastSurface && surface[j]<read) j+=3;
				if (j<lastSurface && surface[j]==read) j=keepNonexposed_maxIntersection_impl(read,tmpWindows,lastWindow,surface,j,lastSurface,minIntervalLength,intersection,tmpArray2,tmpLong);
				read=node.read;
				lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(tmpWindows,0,node);
				continue;
			}
			lastWindow++;
			tmpWindows=IntervalGraphStep3.printKernelTags_ensureArray(tmpWindows,lastWindow+1,CAPACITY);
			IntervalGraphStep3.printKernelTags_setArray(tmpWindows,lastWindow,node);
		}
		while (j<lastSurface && surface[j]<read) j+=3;
		if (j<lastSurface && surface[j]==read) keepNonexposed_maxIntersection_impl(read,tmpWindows,lastWindow,surface,j,lastSurface,minIntervalLength,intersection,tmpArray2,tmpLong);
		return (int)tmpLong[0];
	}
	
	
	/**
	 * Remark: the procedure uses $tracks,lastTrack,kernelsWithTrack$ as temporary space.
	 *
	 * @param readID $windows[0..lastWindow]$ is assumed to be a nonempty list of 
	 * intervals from this read;
	 * @param fromSurface the first element of $surface$ with read equal to $readID$;
	 * @param minIntervalLength only intervals of the intersection that are 
	 * $>=minIntervalLength$ contribute to the intersection length;
	 * @param intersection output array (the procedure cumulates intersection lengths);
	 * @param tmpArray temporary space, assumed to be at least as long as $surface$;
	 * @param tmpLong input/output array: 0=ID of a kernel with max intersection; 1=max 
	 * intersection;
	 * @return the first element of $surface$ from $fromSurface$, with read different from
	 * $readID$.
	 */
	private static final int keepNonexposed_maxIntersection_impl(int readID, IntervalGraph.Node[] windows, int lastWindow, int[] surface, int fromSurface, int lastSurface, int minIntervalLength, long[] intersection, int[] tmpArray, long[] tmpLong) {
		final int INCREMENT = 16;  // Arbitrary, even.
		int i, j;
		int last, windowStart, windowEnd, trackEnd, kernel, lastKernel, out;
		IntervalGraph.Node window;
		
		// Building the tracks of all non-exposed kernels
		lastKernelWithTrack=-1;
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; windowStart=window.start; windowEnd=window.end;
			lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) {
				kernel=window.kernels[j];
				if (labels[kernel]!=0) continue;
				if (!tmpBoolean[kernel]) {
					tmpBoolean[kernel]=true;
					kernelsWithTrack[++lastKernelWithTrack]=kernel;
				}
				if (lastTrack[kernel]==-1) {
					tracks[kernel][0]=windowStart;
					tracks[kernel][1]=windowEnd;
					lastTrack[kernel]=1;
					continue;
				}
				last=lastTrack[kernel];
				trackEnd=tracks[kernel][last];
				if (windowStart<=trackEnd) tracks[kernel][last]=Math.max(trackEnd,windowEnd);
				else {
					if (last+2>tracks[kernel].length-1) {
						int[] newArray = new int[tracks[kernel].length+INCREMENT];
						System.arraycopy(tracks[kernel],0,newArray,0,tracks[kernel].length);
						tracks[kernel]=newArray;
					}
					tracks[kernel][++last]=windowStart; tracks[kernel][++last]=windowEnd;
					lastTrack[kernel]=last;
				}
			}
		}
		
		// Computing track intersections
		last=-1;
		for (out=fromSurface; out<lastSurface; out+=3) {
			if (surface[out]!=readID) break;
			tmpArray[++last]=surface[out+1];
			tmpArray[++last]=surface[out+2];
		}
		for (i=0; i<=lastKernelWithTrack; i++) {
			kernel=kernelsWithTrack[i];
			intersection[kernel]+=Intervals.intersectionLength(tmpArray,last,tracks[kernel],lastTrack[kernel],minIntervalLength,true);
			if (intersection[kernel]>tmpLong[1]) {
				tmpLong[0]=kernel;
				tmpLong[1]=intersection[kernel];
			}
		}
		
		// Cleaning $tmpBoolean$ and $lastTrack$.
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) {
				kernel=window.kernels[j];
				tmpBoolean[kernel]=false;
				lastTrack[kernel]=-1;
			}
		}
		
		return out;
	}
	
	
	/**
	 * Stores $surface \setminus selectedKernel$ into $tmpArray2[0..X]$, where $X$ is
	 * returned in output.
	 *
	 * Remark: the procedure uses global variables $tmpWindows,tmpArray3$ as temporary 
	 * space.
	 *
	 * @param lastNode the procedure uses only $nodesArray[0..lastNode]$, which are 
	 * assumed to contain all and only the nodes with a non-exposed kernel.
	 */
	private static final int keepNonexposed_subtract(int selectedKernel, int[] surface, int lastSurface, int lastNode) {
		final int CAPACITY = 16;  // Arbitrary
		int i, j, k;
		int read, lastWindow, lastNewSurface, length;
		IntervalGraph.Node node;
		
		lastNewSurface=-1;
		node=nodesArray[0]; read=node.read;
		lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(tmpWindows,0,node); j=0;
		for (i=1; i<=lastNode; i++) {
			node=nodesArray[i];
			if (node.read!=read) {
				k=j;
				while (k<lastSurface && surface[k]<read) k+=3;
				length=k-j;
				if (tmpArray2.length<lastNewSurface+1+length) {
					int[] newArray = new int[lastNewSurface+1+length+CAPACITY*3];
					System.arraycopy(tmpArray2,0,newArray,0,lastNewSurface+1);
					tmpArray2=newArray;
				}
				while (j<lastSurface && surface[j]<read) {
					tmpArray2[++lastNewSurface]=surface[j];
					tmpArray2[++lastNewSurface]=surface[j+1];
					tmpArray2[++lastNewSurface]=surface[j+2];
					j+=3;
				}
				if (j<lastSurface && surface[j]==read) {
					lastNewSurface=keepNonexposed_subtract_impl(selectedKernel,read,tmpWindows,lastWindow,surface,j,lastSurface,lastNewSurface);
					while (j<lastSurface && surface[j]==read) j+=3;
				}
				read=node.read;
				lastWindow=0; IntervalGraphStep3.printKernelTags_setArray(tmpWindows,0,node);
				continue;
			}
			lastWindow++;
			tmpWindows=IntervalGraphStep3.printKernelTags_ensureArray(tmpWindows,lastWindow+1,CAPACITY);
			IntervalGraphStep3.printKernelTags_setArray(tmpWindows,lastWindow,node);
		}
		k=j;
		while (k<lastSurface && surface[k]<read) k+=3;
		length=k-j;
		if (tmpArray2.length<lastNewSurface+1+length) {
			int[] newArray = new int[lastNewSurface+1+length+CAPACITY*3];
			System.arraycopy(tmpArray2,0,newArray,0,lastNewSurface+1);
			tmpArray2=newArray;
		}
		while (j<lastSurface && surface[j]<read) {
			tmpArray2[++lastNewSurface]=surface[j];
			tmpArray2[++lastNewSurface]=surface[j+1];
			tmpArray2[++lastNewSurface]=surface[j+2];
			j+=3;
		}
		if (j<lastSurface && surface[j]==read) {
			lastNewSurface=keepNonexposed_subtract_impl(selectedKernel,read,tmpWindows,lastWindow,surface,j,lastSurface,lastNewSurface);
			while (j<lastSurface && surface[j]==read) j+=3;
		}
		length=lastSurface+1-j;
		if (tmpArray2.length<lastNewSurface+1+length) {
			int[] newArray = new int[lastNewSurface+1+length+CAPACITY*3];
			System.arraycopy(tmpArray2,0,newArray,0,lastNewSurface+1);
			tmpArray2=newArray;
		}
		while (j<lastSurface) {
			tmpArray2[++lastNewSurface]=surface[j];
			tmpArray2[++lastNewSurface]=surface[j+1];
			tmpArray2[++lastNewSurface]=surface[j+2];
			j+=3;
		}
		return lastNewSurface;
	}
	
	
	/**
	 * Appends to $tmpArray2[lastNewSurface+1..X]$ the result of $surface \setminus 
	 * selectedKernel$ limited to $readID$, where X is returned in output.
	 *
	 * Remark: the procedure uses global variable $tmpArray3$ as temporary space.
	 */
	private static final int keepNonexposed_subtract_impl(int selectedKernel, int readID, IntervalGraph.Node[] windows, int lastWindow, int[] surface, int fromSurface, int lastSurface, int lastNewSurface) {
		final int INCREMENT = 16;  // Arbitrary, even.
		int i, j, k;
		int last, windowStart, windowEnd, trackEnd, kernel, lastKernel, nKernelsWithTrack, length;
		IntervalGraph.Node window;
		
		// Building the track of $selectedKernel$.
		lastTrack[selectedKernel]=-1;
		for (i=0; i<=lastWindow; i++) {
			window=windows[i]; windowStart=window.start; windowEnd=window.end;
			lastKernel=window.lastKernel;
			for (j=0; j<=lastKernel; j++) {
				kernel=window.kernels[j];
				if (kernel!=selectedKernel) continue;
				if (lastTrack[kernel]==-1) {
					tracks[kernel][0]=windowStart;
					tracks[kernel][1]=windowEnd;
					lastTrack[kernel]=1;
					break;
				}
				last=lastTrack[kernel];
				trackEnd=tracks[kernel][last];
				if (windowStart<=trackEnd) tracks[kernel][last]=Math.max(trackEnd,windowEnd);
				else {
					if (last==tracks[kernel].length) {
						int[] newArray = new int[tracks[kernel].length+INCREMENT];
						System.arraycopy(tracks[kernel],0,newArray,0,tracks[kernel].length);
						tracks[kernel]=newArray;
					}
					tracks[kernel][++last]=windowStart; tracks[kernel][++last]=windowEnd;
					lastTrack[kernel]=last;
				}
				break;
			}
		}	
		if (lastTrack[selectedKernel]==-1) {
			k=fromSurface;
			while (k<lastSurface && surface[k]==readID) k+=3;
			length=k-fromSurface;
			if (tmpArray2.length<lastNewSurface+1+length) {
				int[] newArray = new int[lastNewSurface+1+length+INCREMENT*3];
				System.arraycopy(tmpArray2,0,newArray,0,lastNewSurface+1);
				tmpArray2=newArray;
			}
			for (i=fromSurface; i<lastSurface; i+=3) {
				if (surface[i]!=readID) break;
				tmpArray2[++lastNewSurface]=surface[i];
				tmpArray2[++lastNewSurface]=surface[i+1];
				tmpArray2[++lastNewSurface]=surface[i+2];
			}
			return lastNewSurface;
		}
		
		// Subtracting tracks
		k=fromSurface;
		while (k<lastSurface && surface[k]==readID) k+=3;
		length=k-fromSurface;
		if (tmpArray3==null || tmpArray3.length<length) tmpArray3 = new int[length+INCREMENT];
		last=-1;
		for (i=fromSurface; i<lastSurface; i+=3) {
			if (surface[i]!=readID) break;
			tmpArray3[++last]=surface[i+1];
			tmpArray3[++last]=surface[i+2];
		}
		return Intervals.difference(tmpArray3,last,tracks[selectedKernel],lastTrack[selectedKernel],tmpArray2,lastNewSurface+1,readID,1);
	}
	
	
	/**
	 * Removes every kernel in $node.kernels$ that has a zero in $labels$. Some intervals 
	 * might have no kernel after this.
	 *
	 * Remark: for simplicity, the program does not reassign intervals with a removed tag 
	 * to other tags. An occurrence of a kernel with not enough exposed occurrences is 
	 * either already contained in an occurrence of another kernel, or it is exposed, and 
	 * in the latter case it is not clear which tag it should get. An occurrence of a 
	 * periodic kernel discarded by $markSameSurface()$ could be assigned to the 
	 * representative of its connected component, but it is not clear in which 
	 * orientation.
	 *
	 * Remark: the procedure uses $tmpArray1$ and $tmpByte1$, which are assumed to be of 
	 * length at least $nPathKernels$ each.
	 *
	 * @param ancestors one row per element of $kernels25$.
	 */
	private static final void removeKernels(IntervalGraph.Node node) {
		int i, j, k;
		int kernel, lastKernel, last1;
		
		// Updating kernels and orientations
		lastKernel=node.lastKernel; last1=-1;
		for (i=0; i<=lastKernel; i++) {
			kernel=node.kernels[i];
			if (labels[kernel]==0) continue;
			last1++;
			tmpArray1[last1]=kernel;
			tmpByte1[last1]=node.kernelOrientations[i];
		}
		if (last1>=0) {
			System.arraycopy(tmpArray1,0,node.kernels,0,last1+1);
			System.arraycopy(tmpByte1,0,node.kernelOrientations,0,last1+1);
		}
		node.lastKernel=last1;
		
		// Updating start and end
		if (node.lastPathWithStart!=-1) {
			last1=-1;
			for (i=0; i<=node.lastPathWithStart; i++) {
				kernel=node.pathsWithStart[i];
				k=kernel>=0?kernel:-1-kernel;
				if (labels[k]==0) continue;
				tmpArray1[++last1]=kernel;
			}
			if (last1!=-1) System.arraycopy(tmpArray1,0,node.pathsWithStart,0,last1+1);
			node.lastPathWithStart=last1;
		}
		if (node.lastPathWithEnd!=-1) {
			last1=-1;
			for (i=0; i<=node.lastPathWithEnd; i++) {
				kernel=node.pathsWithEnd[i];		
				k=kernel>=0?kernel:-1-kernel;
				if (labels[k]==0) continue;
				tmpArray1[++last1]=kernel;
			}
			if (last1!=-1) System.arraycopy(tmpArray1,0,node.pathsWithEnd,0,last1+1);
			node.lastPathWithEnd=last1;
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// ------------------------------ UNUSED PROCEDURES ----------------------------------
	
	
	private static final void markSameSurface_printAllIntersections(int[][] neighbors, int[] lastNeighbor) {
		final double MIN_JACCARD = 0.1;  // Arbitrary
		final int PEN_SCALE = 10;
		int i, j;
		int edgeSize;
		double jaccard;
		String color;
		
		System.err.println("graph G {");
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]==0) continue; 
			switch (pathKernelPeriodic[i]) {
				case 0: color="white"; break;
				case 1: color=Colors.type2color_dot(Constants.INTERVAL_PERIODIC); break;
				case 2: color=Colors.COLOR_DOT_LONGPERIOD; break;
				default: color="";
			}
			System.err.println(i+" [label=\""+kernel2descriptor[i][0]+"-"+kernel2descriptor[i][1]+"-"+kernel2descriptor[i][2]+"\",style=\"filled\",fillcolor=\""+color+"\"];");
		}
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]==0) continue;
			for (j=0; j<=lastNeighbor[i]; j+=3) {
				jaccard=((double)neighbors[i][j+1])/neighbors[i][j+2];
				if (jaccard>MIN_JACCARD) {
					edgeSize=Math.max(1,(int)(jaccard*PEN_SCALE));
					System.err.println(i+" -- "+neighbors[i][j]+" [weight=\""+IO.format(jaccard)+"\",penwidth=\""+edgeSize+"\"];");
				}
			}
		}
		System.err.println("}");
	}
	
	
	/**
	 * Fits a density estimation tree on all Jaccard similarities between periodic kernels
	 * in $neighbors$, and returns a value that separates the last local maximum from the
	 * rest.
	 *
	 * Remark: the procedure assumes $tmpPoints,tmpLeaves$ to be of size at least equal to
	 * the number of periodic edges in the matrix.
	 *
	 * @param minThreshold the returned value is not smaller than this.
	 */
	public static final double getJaccardThreshold(int[][] neighbors, int[] lastNeighbor, double minThreshold) {
		final int N_CONSECUTIVE_POINTS = 2;  // Arbitrary
		final double QUANTILE = 0.75;  // Arbitrary
		final int MIN_INTERVAL_LENGTH_FACTOR = 20;  // Arbitrary
		final double MAX_DIFFERENCE = 0.3;  // Arbitrary
		boolean isPeriodic;
		int i, j;
		int lastJaccard, nLocalMaximumLeaves, nRuns;
		int lastLeaf, lastPoint, nHigh;
		double minIntervalLength, minLocalMaxDistance, threshold;
		Leaf[] lvs;
		
		lastJaccard=-1;
		for (i=0; i<nPathKernels; i++) {
			if (labels[i]==0 || pathKernelPeriodic[i]==0) continue;
			for (j=0; j<=lastNeighbor[i]; j+=3) {
				if (neighbors[i][j+1]==0) continue;
				lastJaccard++;
				tmpPoints[lastJaccard].position=((double)neighbors[i][j+1])/neighbors[i][j+2];
				tmpPoints[lastJaccard].mass=1;
			}
		}
		lastJaccard=Points.sortAndCompact(tmpPoints,lastJaccard);
		threshold=minThreshold;
		if (lastJaccard+1<=Points.FEW_POINTS) {
			if (tmpPoints[lastJaccard].position-tmpPoints[0].position>MAX_DIFFERENCE) {
				i=Points.getRoughThreshold(tmpPoints,lastJaccard,true,0,true);
				if (i!=-1) threshold=Math.max(threshold,tmpPoints[i].position);
			}
		}
		else {	
			minIntervalLength=Math.max( Points.distanceQuantile(tmpPoints,lastJaccard,N_CONSECUTIVE_POINTS,false,QUANTILE),
										(tmpPoints[lastJaccard].position-tmpPoints[0].position)/MIN_INTERVAL_LENGTH_FACTOR
									  );
			minLocalMaxDistance=2*minIntervalLength;
			if (!Points.areUniformlyDistributed(tmpPoints,0,lastJaccard,false,(tmpPoints[lastJaccard].position-tmpPoints[0].position)/Points.DEFAULT_NBINS)) {
				DensityEstimationTree.allocateMemory(lastJaccard+1);
				nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(tmpPoints,0,lastJaccard,minIntervalLength,minLocalMaxDistance,false,-1,-1,-1,false,false);				
				if (nLocalMaximumLeaves>0) {
					lastLeaf=DensityEstimationTree.lastLeaf;
					lvs=DensityEstimationTree.leaves;
					DensityEstimationTree.leaves=tmpLeaves;
					tmpLeaves=lvs;
					nHigh=Leaves.setHighLocalMax(tmpLeaves,lastLeaf,lastJaccard,false,false,tmpPoints,Points.mass(tmpPoints,0,lastJaccard),true,Math.POSITIVE_INFINITY);
					lvs=tmpLeaves;
					tmpLeaves=DensityEstimationTree.leaves;
					DensityEstimationTree.leaves=lvs;
					DensityEstimationTree.lastLeaf=lastLeaf;
					DensityEstimationTree.markRunsOfLocalMaximumLeaves(tmpPoints);
					nRuns=DensityEstimationTree.getNumberOfRuns();
					if (nRuns==1) threshold=Math.max(threshold,DensityEstimationTree.separateRun(0,tmpPoints,0,lastJaccard,false));
					else threshold=Math.max(threshold,DensityEstimationTree.separateRun(nRuns-2,tmpPoints,0,lastJaccard,false));
				}
				DensityEstimationTree.deallocateMemory();
			}
		}	
		return threshold;
	}
	
	
	
	
	
	
	/**
	 * @param trackStart sorted;
	 * @return TRUE iff $[start..end]$ is approximately contained in a track.
	 */
	private static final boolean containedInTrack(int start, int end, int[] trackStart, int[] trackEnd, int lastTrack) {
		int i, j;
		int cStart, cEnd;
		
		i=Arrays.binarySearch(trackStart,0,lastTrack+1,start);
		if (i<0) i=-i-1;
		for (j=i; j<=lastTrack; j++) {
			cStart=trackStart[j];
			if (cStart>=end) break;
			cEnd=trackEnd[j];
			if (Intervals.isApproximatelyContained(start,end,cStart,cEnd) || Intervals.areApproximatelyIdentical(start,end,cStart,cEnd)) return true;
		}
		for (j=i-1; j>=0; j--) {
			cStart=trackStart[j]; cEnd=trackEnd[j];
			if (Intervals.isApproximatelyContained(start,end,cStart,cEnd) || Intervals.areApproximatelyIdentical(start,end,cStart,cEnd)) return true;
		}
		return false;
	}
	
	
	/**
	 * Decides if there is a sequence of tracks $(i_0,j_0),...,(i_k,j_k)$, where $k>=1$, 
	 * such that the 0-th and k-th tracks straddle $[start..end]$, and all other tracks
	 * are adjacent to one another and to the 0-th and k-th track.
	 *
	 * @param marks cells have the following meaning: 1=$[start..end]$ straddles >=2
	 * consecutive tracks; 2/3=a track covers the immediate left/right of $[start..end]$; 
	 * 4/5=a track straddles $[start..end]$ on the left/right (cell 0 is unused).
	 */
	private static final void straddlingTracks(int start, int end, int[] trackStart, int[] trackEnd, int lastTrack, int straddlingThreshold, int identityThreshold, boolean[] marks) {
		int i, j, k;
		int startJ, endJ, startK, endK;
		
		i=Arrays.binarySearch(trackStart,0,lastTrack+1,start);
		if (i<0) i=-i-1;
		for (j=i; j<=lastTrack; j++) {
			startJ=trackStart[j]; 
			if (startJ>end+identityThreshold) break;
			endJ=trackEnd[j];
			if (endJ<end+straddlingThreshold) continue;
			marks[3]=true;
			if (startJ>end-straddlingThreshold || startJ<start+straddlingThreshold) continue;
			marks[5]=true;
			for (k=i-1; k>=0; k--) {
				startK=trackStart[k]; endK=trackEnd[k];
				if (endK<start-identityThreshold || startK>start-straddlingThreshold) continue;
				marks[2]=true;
				if (endK<start+straddlingThreshold || endK>end-straddlingThreshold) continue;
				marks[4]=true;
				if (Math.abs(endK,startJ)<=identityThreshold || straddlingTracks_adjacentSequence(endK,startJ,k+1,j-1,trackStart,trackEnd,identityThreshold)) {
					marks[1]=true;
					return;
				}
			}
		}
	}
	
	
	/**
	 * Remark: the procedure uses the global array $tmpArray1$. 
	 *
	 * @param trackStart sorted;
	 * @return TRUE iff there is a sequence of adjacent tracks in $[from..to]$ (possibly
	 * just one track) such that the beginning of the first track is close to $from$ and 
	 * the end of the last track is close to $to$.
	 */
	private static final boolean straddlingTracks_adjacentSequence(int from, int to, int firstTrack, int lastTrack, int[] trackStart, int[] trackEnd, int identityThreshold) {
		final int GROWTH_RATE = 10;  // Arbitrary
		int i, j;
		int trackID, start, end, top;
		
		for (i=firstTrack; i<=lastTrack; i++) {
			if (trackStart[i]>from+identityThreshold) break;
			if (trackStart[i]<from-identityThreshold) continue;
			tmpArray1[0]=i; top=0;
			while (top>=0) {
				trackID=tmpArray1[top--];
				end=trackEnd[trackID];
				if (Math.abs(end,to)<=identityThreshold) return true;
				for (j=trackID+1; j<=lastTrack; j++) {
					start=trackStart[j];
					if (start>end+identityThreshold) break;
					if (start<end-identityThreshold) continue;
					top++;
					if (top>=tmpArray1.length) {
						int[] newArray = new int[tmpArray1.length+GROWTH_RATE];
						System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
						tmpArray1=newArray;
					}
					tmpArray1[++top]=j;
				}
			}
		}
		return false;
	}
	
	
	
	
}