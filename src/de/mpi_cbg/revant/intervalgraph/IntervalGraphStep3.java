package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Leaf;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Histograms;
import de.mpi_cbg.revant.util.Leaves;
import de.mpi_cbg.revant.util.DAG;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Colors;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Factorize;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.ReadA;
import de.mpi_cbg.revant.factorize.Intervals;
import de.mpi_cbg.revant.factorize.PeriodicSubstrings;


/**
 * Step 3: kernels inside clusters.
 *
 * Remark: some kernel descriptors produced by this program can have low frequency, for at
 * least the following reasons: (1) there is no min. frequency constraint for selecting 
 * periodic kernels (see $markRepeats_pathKernelPeriodic()$); (2) rare kernels of any type
 * are indeed discarded at the end of the program (see $discardRareKernels()$), but the 
 * threshold is very low; (3) if no path kernel is selected throughout the program, all of
 * them are output in order not to lose information.
 */
public class IntervalGraphStep3 {
	/**
	 * Parameters of the pipeline
	 */
	private static int MAX_ALIGNMENTS_PER_READ;
	private static String GRAPH_FILE;
	public static final int NEW_OVERLAPS_THRESHOLD = IO.quantum; //(3*IntervalGraph.IDENTITY_THRESHOLD)>>1;  // Arbitrary
	public static final int SOURCE_CONTAINMENT_IDENTITY = 3;  // Must be different from all $IntervalGraph.CONTAINMENT_$ constants.
	private static final int BLOCK_SIZE = 8;  // In $buildNodes2Kernel()$.
	private static final int BLOCK_SIZE_PRIME = BLOCK_SIZE-2;
	private static final int INSERTION_DISTANCE_THRESHOLD = IO.quantum<<1;  // Arbitrary
	public static int maxReadLength;  // Needs to be initialized, since it is used by $buildNodes2Kernel*(),buildKernelGraph_impl()$.
	
	/**
	 * Statistics on kernels
	 */
	public static int nKernels, nPathKernels;
	public static int[] nKernelNodes;  // Interval graph nodes that form each path kernel
	public static int[] kernelSize;  // N. of intervals in the basin of a kernel.
	public static int[][] kernelFrequency;
	private static int[][] kernelStats;
	public static int[] pathKernelLengths;  // String lengths. 0=cyclic kernel.
	public static int nPathKernelCyclic;  // Number of zeros in $pathKernelLengths$.
	public static int[] pathKernel2Kernel;
	public static byte[] pathKernelPeriodic;  // Marks path kernels that contain at least one periodic interval (1=short, 2=long).
	public static int nPathKernelPeriodic;  // Number of marks in $pathKernelPeriodic$.
	private static int[] pathKernelLongestNode;  // Length of a longest node in each path kernel
	private static int[] sorted2original, original2sorted;  // The kernel graph, topologically sorted in reverse order.
	
	/**
	 * All nodes that belong to a kernel
	 */
	public static int[] nodesInKernel;
	public static int nNodesInKernel;
	public static int[] blockKernels;
	
	/**
	 * Description of the kernel graph
	 */
	public static int[] lastInNeighbor, lastOutNeighbor;
	public static int[][] inNeighbors, outNeighbors;
	public static byte[][] inNeighborFlags;  // Each entry of $inNeighbors[i]$ corresponds to a tuple $(i_1,a_1,i_2,a_2)$, where 1=forward and 2=RC orientation; $i_x$: insertion type (see values in $addInsertions[][]$); $a_x=1$ iff all nodes of kernel $i$ are contained in the in-neighbor of kernel $i$ (0 otherwise).
	private static int[][] outNeighbors_ps, inNeighbors_ps, outNeighbors_f, inNeighbors_f;
	private static int[] lastOutNeighbor_ps, lastInNeighbor_ps, lastOutNeighbor_f, lastInNeighbor_f;
	private static int[][] descendants;  // Descendants of each node of the kernel graph
	private static int[] lastDescendant;
	
	/**
	 * Data structures for checking if an alignment or a substring is contained in a 
	 * short-period range.
	 */
	public static boolean[] shortPeriodAlignments = null; 
	public static IntervalGraph.ShortPeriodWindow[] shortPeriodWindows = null;
	public static int[] trackPointers = null;
	
	/**
	 * Temporary space
	 */
	public static Point[] tmpPoints;
	private static int lastTmpPoint;
	public static Leaf[] tmpLeaves;
	private static int[] tmpArray1 = new int[50];
	private static int[] tmpArray2 = new int[50];
	private static int[] tmpArray3 = new int[50];
	private static int[] tmpArray4 = new int[50];
	private static int[] tmpArray5 = new int[50];
	private static byte[] tmpByte1 = new byte[50];
	private static byte[] tmpByte2 = new byte[50];
	private static byte[] tmpByte3 = new byte[50];
	private static boolean[] tmpBoolean1 = new boolean[50];
	private static boolean[] tmpBoolean2 = new boolean[50];
	public static int[] tmpIO = new int[4];
	private static IntervalGraph.Node[] tmpNodes = new IntervalGraph.Node[50];
	private static Node2KernelWindow[] tmpWindows;
	private static IntervalGraph.Node tmpNode = new IntervalGraph.Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
	private static int[][] kernel2NewKernels;  // Used by $markRepeats*()$ procedures.
	private static byte[][] kernel2NewKernels_orientation;
	private static int[][] kernel2NewKernels_start, kernel2NewKernels_end;
	private static int[] lastNewKernel, kernel2NewKernels_lastStart, kernel2NewKernels_lastEnd;
	private static int insertionDistanceStart, insertionDistanceEnd;  // Written by $buildNodes2Kernel_getInsertionDistance()$.
	public static FrequencyPair[] frequencyPairs;
	
	
	public static void main(String[] args) throws IOException {
		final String INPUT_DIR = args[0];
		final String INPUT_GRAPH_ID = args[1];
		GRAPH_FILE=INPUT_DIR+"/"+INPUT_GRAPH_ID;
		IO.coverage=Integer.parseInt(args[2]);
		IO.minOccurrencesInGenome=Integer.parseInt(args[3]);
		IO.minRepeatCoverage=(IO.minOccurrencesInGenome*IO.coverage)-1;
		final int MIN_NODES_PER_KERNEL = IO.minRepeatCoverage+1;
		Alignments.minAlignmentLength=Integer.parseInt(args[4]);
		MAX_ALIGNMENTS_PER_READ=Integer.parseInt(args[5]);
		final int N_READS = Integer.parseInt(args[6]);
		final String QUALITY_THRESHOLDS_FILE = args[7];
		final int N_ALIGNMENTS = Integer.parseInt(args[8]);
		final String OUTPUT_DIR = args[9];
		final String OUTPUT_PREFIX = args[10];
		final String CHECK_LENGTHS_FILE = args[11];
		final double MIN_PERIODIC_FRACTION = 0.8;  // Arbitrary
		
		final String READ_IDS_FILE = INPUT_DIR+"/"+INPUT_GRAPH_ID+"-"+IO.READS_IDS;
		final String READ_LENGTHS_FILE = INPUT_DIR+"/"+INPUT_GRAPH_ID+"-"+IO.READS_LENGTHS;
		final String QUALITIES_FILE = INPUT_DIR+"/"+INPUT_GRAPH_ID+"-"+IO.READS_PHRED_PREFIX+IO.getPhredSuffix(INPUT_DIR+"/"+INPUT_GRAPH_ID+"-"+IO.READS_PHRED_PREFIX);
		final String SHORTPERIOD_TRACKS_FILE = INPUT_DIR+"/"+INPUT_GRAPH_ID+"-"+IO.READS_SHORTPERIOD;
		final String ALIGNMENTS_FILE = INPUT_DIR+"/"+INPUT_GRAPH_ID+"-"+IO.ALIGNMENTS_FILE;
		final String SHORTPERIOD_BITVECTOR_FILE = INPUT_DIR+"/"+INPUT_GRAPH_ID+"-"+IO.ALIGNMENTS_SHORTPERIOD;
		
		boolean overlapAdded;
		int i, j;
		int maxKernelsPerInterval, refined, lastMinimal, minFrequency, nRepeats;
		byte[] kernelLabels;
		int[] new2old, kernelFrequencyPrime, minimalVertices, counts;
		int[][] newKernelFrequency;
		IO.initialize();
		
		// Loading the graph.
		// Remark: after loading, $avgDiffs$ in the graph is a sum, not an avg., and
		// nodes are sorted by $read,start$.
		System.err.println("STEP3> Loading interval graph...");
		new2old=IntervalGraph.deserialize(GRAPH_FILE+".graph",true,true);
		if (IO.CONSISTENCY_CHECKS) {
			System.err.println("STEP3> Consistency checks started...");
			if (IntervalGraph.nNodes==1 && IntervalGraph.nodesArray[0].type!=Constants.INTERVAL_PERIODIC) {
				System.err.println("STEP3> ERROR: the graph consists of a single node, which is not periodic:");
				System.err.println(IntervalGraph.nodesArray[0]);
				System.exit(1);
			}
			IntervalGraph.checkConsistency(0,IntervalGraph.nNodes-1,new int[IntervalGraph.getMaxDegree()*10]);
			System.err.println("STEP3> Consistency checks passed");
		}
		IntervalGraph.printNew2OldArray(new2old,GRAPH_FILE+".new2old");
		new2old=null;
		if (getShortPeriodFraction()>=MIN_PERIODIC_FRACTION) {
			if (IntervalGraph.nNodes>=MIN_NODES_PER_KERNEL) {
				Reads.loadReadIDs(READ_IDS_FILE,N_READS);
				maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
				Reads.loadQualities(QUALITY_THRESHOLDS_FILE,(new File(QUALITIES_FILE)).exists()?QUALITIES_FILE:null);
				handleShortPeriodGraph(OUTPUT_DIR,OUTPUT_PREFIX,INPUT_GRAPH_ID,GRAPH_FILE);
			}
			else {
				for (i=0; i<IntervalGraph.nNodes; i++) IntervalGraph.nodesArray[i].lastKernel=-1;
				printKernelTags(OUTPUT_DIR+"/"+IO.TAGS_PREFIX+OUTPUT_PREFIX+"-"+INPUT_GRAPH_ID+".txt",OUTPUT_PREFIX+"."+INPUT_GRAPH_ID,true,false,-1);
			}
			return;
		}
		
		// Turning on all edges. Supplement edges added later will be set to ON by their
		// construction procedures.
		for (i=0; i<IntervalGraph.nNodes; i++) {
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) IntervalGraph.neighbors[i][j].on=true;
		}
		IntervalGraph.stack = new int[IntervalGraph.nNodes];
		System.err.println("STEP3> Interval graph loaded.  nNodes="+IntervalGraph.nNodes);
		
		// Loading read IDs, lengths, and qualities (used by $trimOverhang()$).
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadQualities(QUALITY_THRESHOLDS_FILE,(new File(QUALITIES_FILE)).exists()?QUALITIES_FILE:null);

		
		// -------------------------- FIXING THE TOPOLOGY --------------------------------
		IntervalGraph.setPeriodicSubstringsOfNodes_sameRead(true);
		IntervalGraph.mergePeriodicSubstringsOfNodes();		
		IntervalGraph.fixUndersplits(true);
		// Remark: resetting the $periodicSubstrings$ list of nodes after fixing
		// undersplits is not necessary, since the lists of existing nodes are still
		// correct, and new nodes are either periodic or have no periodic substring.		
		IntervalGraph.sortNodesArray();  // Needed because $IntervalGraph.nodesArray$ is not necessarily sorted at this point. Several procedures downstream assume it to be sorted, including $checkKernelsWithReads()$, $buildKernelGraph() -> buildNodes2Kernel()$, $basinStats() -> getFrequency()$.
		IntervalGraphStep2.containment2insertion(false,IntervalGraph.stack);  // Transforming containment edges into insertion edges
		IntervalGraphStep2.overlap2sharedSubstring(false,IntervalGraph.stack);
		

		// ---------------------------- BUILDING KERNELS ---------------------------------
		System.err.println("STEP3> Building kernels...");
		// Remark: we need the containing/contained marks local to this graph, rather
		// than to the graph from which this cluster may come from.
		// Remark: $fixUndersplits()$ can create new containment edges, so
		// $markContainedContainingIntervals()$ must be called after it.
		IntervalGraph.markContainedContainingIntervals(false);
		getKernels(tmpArray1);
		nKernels=tmpArray1[0];
		if (nKernels==0) {
			System.err.println("STEP3> REMARK: the graph has no kernel. This might still be correct. nNodes="+IntervalGraph.nNodes);
			System.exit(1);
		}
		System.err.println("STEP3> Adding supplementary edges...");
		IntervalGraph.addIdentityEdges(true,ALIGNMENTS_FILE);
		IntervalGraph.addContainmentEdges(true,ALIGNMENTS_FILE);
		for (i=0; i<IntervalGraph.nNodes; i++) IntervalGraph.nodesArray[i].visited=0;
		do { overlapAdded=IntervalGraph.addEdgesToOverlapNeighbors(true,ALIGNMENTS_FILE); }
		while (overlapAdded);
		IntervalGraphStep2.containment2insertion(true,IntervalGraph.stack);
		IntervalGraphStep2.overlap2sharedSubstring(true,IntervalGraph.stack);
		kernelStats=countRepeatsPerKernel(tmpIO);

if (IO.SHOW_INTERACTIVE) {
System.err.println("KERNEL STATS:");		
for (int x=0; x<nKernels; x++) {
	System.err.println(x+" :: nAssemblies="+kernelStats[x][0]+" minStringLength="+kernelStats[x][1]+" maxStringLength="+kernelStats[x][2]+" nTransitive="+kernelStats[x][3]+" oneEnd="+kernelStats[x][4]+" noNeighbor="+kernelStats[x][5]+" connectedComponents="+kernelStats[x][6]);
}
System.err.println("nodesInKernel:");
for (int x=0; x<nNodesInKernel; x++) System.err.println(IntervalGraph.nodesArray[nodesInKernel[x]]);
}
		
		nPathKernels=tmpIO[1]; nPathKernelCyclic=tmpIO[2];
		Points.allocateMemory(nPathKernels);
		if (tmpPoints==null || tmpPoints.length<nPathKernels) {
			tmpPoints = new Point[nPathKernels];
			for (i=0; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
		}
		Leaves.allocateMemory(nPathKernels);
		if (tmpLeaves==null || tmpLeaves.length<nPathKernels) {
			tmpLeaves = new Leaf[nPathKernels];
			for (i=0; i<tmpLeaves.length; i++) tmpLeaves[i] = new Leaf();
		}
		System.err.println("STEP3> "+nKernels+" kernels built:");
		printKernelHistograms();
		System.err.println();

		
		// -------------------------------- BASINS ---------------------------------------
		System.err.println("STEP3> Computing basins of "+nPathKernels+" path kernels...  ("+IntervalGraph.nNodes+" nodes)");
		getBasins(nPathKernels,true);  // Some nodes might get no kernel tag
		if (IO.CONSISTENCY_CHECKS) checkKernelLengths(0,GRAPH_FILE,CHECK_LENGTHS_FILE);
		kernelSize = new int[nPathKernels];
		kernelFrequency = new int[nPathKernels][4];
		frequencyPairs = new FrequencyPair[100];  // Arbitrary, needed by $getFrequency$.
		for (i=0; i<frequencyPairs.length; i++) frequencyPairs[i] = new FrequencyPair();
		descendants=null; lastDescendant=null;
		basinStats(tmpIO,true);
		maxKernelsPerInterval=tmpIO[0];
		System.err.println("STEP3> Basins computed");
		printBasinHistograms(tmpIO);
		kernelLabels = new byte[nPathKernels];
		if (nPathKernels==1) {
			nodesWithNoKernel(ALIGNMENTS_FILE,SHORTPERIOD_BITVECTOR_FILE,SHORTPERIOD_TRACKS_FILE,N_ALIGNMENTS);
			if (IntervalGraph.nNodes>=MIN_NODES_PER_KERNEL) {
				basinStats(tmpIO,true);
				kernelLabels[0]=1;
				printKernelTags(OUTPUT_DIR+"/"+IO.TAGS_PREFIX+OUTPUT_PREFIX+"-"+INPUT_GRAPH_ID+".txt",OUTPUT_PREFIX+"."+INPUT_GRAPH_ID,true,false,-1);
				printBasinDescriptors(kernelLabels,OUTPUT_DIR,IO.BASIN_DESCRIPTOR_PREFIX+"-"+OUTPUT_PREFIX+"-"+INPUT_GRAPH_ID,null,false);
			}
			else {
				for (i=0; i<IntervalGraph.nNodes; i++) IntervalGraph.nodesArray[i].lastKernel=-1;
				printKernelTags(OUTPUT_DIR+"/"+IO.TAGS_PREFIX+OUTPUT_PREFIX+"-"+INPUT_GRAPH_ID+".txt",OUTPUT_PREFIX+"."+INPUT_GRAPH_ID,true,false,-1);
			}
			Points.deallocateMemory();
			return;
		}		
		System.err.println("STEP3> Refining basins with sequential information...");		
		refined=checkKernelsWithReads(MAX_ALIGNMENTS_PER_READ,nPathKernels,true);
		basinStats(tmpIO,true);
		maxKernelsPerInterval=tmpIO[0];
		System.err.println("STEP3> Path kernels refined in "+refined+" nodes ("+(100*((double)refined)/IntervalGraph.nNodes)+"% of all nodes)");
		printBasinHistograms(tmpIO);		
		
		// ----------------------------- KERNEL GRAPH ------------------------------------
		System.out.println("STEP3> Building the kernel graph");
		nPathKernelPeriodic=buildKernelGraph(N_ALIGNMENTS,ALIGNMENTS_FILE,SHORTPERIOD_BITVECTOR_FILE,SHORTPERIOD_TRACKS_FILE);
		basinStats(tmpIO,true);
		updateNKernelNodes();  // Must be done after all uses of $nKernelNodes$ in the original format.
		printKernelGraph(GRAPH_FILE+".kernelGraph-before.dot",null,kernelFrequency,null);
		minFrequency=getKernelFrequencyThreshold();
		System.err.println("STEP3> minFrequency="+minFrequency+"  nodes in the interval graph: "+IntervalGraph.nNodes);
		newKernelFrequency=null; nRepeats=0;
		if (minFrequency<IntervalGraph.nNodes) {
			buildFrequencyGraph_prefixSuffix(); buildFrequencyGraph_full();
		
		
if (IO.SHOW_INTERACTIVE_2) {
System.err.println("inNeighbors_ps:");	
for (int x=0; x<nPathKernels<<1; x++) {
	System.err.print(x+": ");
	for (int y=0; y<=lastInNeighbor_ps[x]; y++) System.err.print(inNeighbors_ps[x][y]+",");
	System.err.println();
}
System.err.println("outNeighbors_ps:");	
for (int x=0; x<nPathKernels<<1; x++) {
	System.err.print(x+": ");
	for (int y=0; y<=lastOutNeighbor_ps[x]; y++) System.err.print(outNeighbors_ps[x][y]+",");
	System.err.println();
}
System.err.println("inNeighbors_f:");	
for (int x=0; x<nPathKernels; x++) {
	System.err.print(x+": ");
	for (int y=0; y<=lastInNeighbor_f[x]; y++) System.err.print(inNeighbors_f[x][y]+",");
	System.err.println();
}
System.err.println("outNeighbors_f:");	
for (int x=0; x<nPathKernels; x++) {
	System.err.print(x+": ");
	for (int y=0; y<=lastOutNeighbor_f[x]; y++) System.err.print(outNeighbors_f[x][y]+",");
	System.err.println();
}
}


			newKernelFrequency = new int[nPathKernels][4];
			nRepeats=markRepeats(minFrequency,kernelLabels,newKernelFrequency); // It might happen that no path kernel gets a positive number in $kernelLabels$, and that some interval graph nodes get no kernel tag. 
			deallocateFrequencyGraphs();
			System.err.println("STEP3> found "+nRepeats+" independent repeats (out of "+nPathKernels+" path kernels)");
		}
		else {
			// Marking all non-periodic kernel paths as non-redundant.
			// Remark: their kernel tags in all interval graph nodes are already correct.
			for (i=0; i<nPathKernels; i++) {
				if (pathKernelLengths[i]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[i]==0)) kernelLabels[i]=1;
			}
			nRepeats=nPathKernels;
		}
		if (nPathKernelPeriodic>0 || nPathKernelCyclic>0) {
			nRepeats+=markRepeats_pathKernelPeriodic(kernelLabels,false);
			printKernelGraph(GRAPH_FILE+".kernelGraph-4.dot",kernelLabels,newKernelFrequency!=null?newKernelFrequency:kernelFrequency,null);
		}
		if (nRepeats==0) {
			printKernelTags(OUTPUT_DIR+"/"+IO.TAGS_PREFIX+OUTPUT_PREFIX+"-"+INPUT_GRAPH_ID+".txt",OUTPUT_PREFIX+"."+INPUT_GRAPH_ID,true,false,-1);
			Points.deallocateMemory();
			return;
		}
		updateBidirectedGraphNodePointer();
		nodesWithNoKernel(ALIGNMENTS_FILE,SHORTPERIOD_BITVECTOR_FILE,SHORTPERIOD_TRACKS_FILE,N_ALIGNMENTS);
		shortPeriodAlignments=null; shortPeriodWindows=null; trackPointers=null;
		basinStats(tmpIO,false);
		if (IO.CONSISTENCY_CHECKS) checkKernelTags(kernelLabels);
		nRepeats=discardRareKernels(MIN_NODES_PER_KERNEL);
		if (nRepeats==0) {
			printKernelTags(OUTPUT_DIR+"/"+IO.TAGS_PREFIX+OUTPUT_PREFIX+"-"+INPUT_GRAPH_ID+".txt",OUTPUT_PREFIX+"."+INPUT_GRAPH_ID,true,false,-1);
			Points.deallocateMemory();
			return;
		}
		basinStats(tmpIO,false);
		getFrequencyPeriodic();
		printKernelTags(OUTPUT_DIR+"/"+IO.TAGS_PREFIX+OUTPUT_PREFIX+"-"+INPUT_GRAPH_ID+".txt",OUTPUT_PREFIX+"."+INPUT_GRAPH_ID,true,false,-1);
		printBasinDescriptors(kernelLabels,OUTPUT_DIR,IO.BASIN_DESCRIPTOR_PREFIX+"-"+OUTPUT_PREFIX+"-"+INPUT_GRAPH_ID,null,false);
		Points.deallocateMemory();
		if (IO.CONSISTENCY_CHECKS) checkKernelLengths(1,GRAPH_FILE,CHECK_LENGTHS_FILE);
		
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		Math.set(tmpArray1,nPathKernels-1,0);
		for (i=0; i<IntervalGraph.nNodes; i++) tmpArray1[IntervalGraph.nodesArray[i].lastKernel+1]++;
		System.err.println("Number of nodes with K kernel tags:");
		for (i=0; i<nPathKernels; i++) System.err.println(i+": "+tmpArray1[i]);
	}
	
	
	
	


	
	
	// ------------------------ KERNELS AND REPEATS PER KERNEL ---------------------------
	
	/**
	 * If two overlapping intervals are also contained in one another, they are identical,
	 * or they have an insertion, then the overlap is not trusted, since it could be 
	 * explained by the other relationships.
	 */
	private static final boolean isOverlap(IntervalGraph.Edge edge) {
		return edge.overlap!=-1 && edge.containment==-1 && edge.insertion==-1;
	}
	
	
	/**
	 * Identity relations are always trusted.
	 */
	private static final boolean isIdentity(IntervalGraph.Edge edge) {
		return edge.containment==Constants.CONTAINMENT_IDENTICAL;
	}
	
	
	/**
	 * A \emph{kernel} is a (possibly singleton) connected component of maximal intervals
	 * (in the partial order of edge containment) induced by non-supplement overlap edges 
	 * (as decided by procedure $isOverlap()$) and identity edges (as decided by 
	 * $isIdentity()$). The procedure sets the $kernels$ variable of selected nodes.
	 *
	 * Remark: the same Step 2 cluster might contain multiple kernels.
	 *
	 * Remark: we don't use insertion edges in the partial order, since they might connect
	 * e.g. the full copy of a repeat to an interval that is the concatenation of distinct
	 * repeats. Symmetrically, this implies that fragments of the same repeat might form
	 * different kernels, thus we will need to merge all kernels related to the same
	 * repeat down the line.
	 *
	 * Remark: intervals that are contained inside a union of short-period periodic 
	 * substring intervals might not be used to build kernels, since it is not clear
	 * whether they correspond to independent repeats or not. If they are not, i.e. if 
	 * they are random, then assembling them might just reconstruct the containing 
	 * periodic substring itself, and they might form overlap forks. However, the current 
	 * procedure uses them as well, delegating to the following steps of the pipeline the
	 * decision on whether they should be reported or not.
	 *
	 * Remark: assume that some intervals produced by factorization are wrong, in the 
	 * sense that they are random substrings of a repeat, rather than a full copy of the
	 * repeat. If they form a kernel, such kernel, after assembly, should approximate the 
	 * full copy of the repeat, or a long substring thereof that is not an independent
	 * repeat, and thus it will be discarded by the later stages of the pipeline.
	 *
	 * Remark: kernels are related to assembly paths in \cite{li2007novel,
	 * bandyopadhyay2010repfrag}. However, in those papers they select paths greedily 
	 * based on the weight of vertices and edges. They also proceed iteratively, removing 
	 * all substrings of reads that belong to a path and updating their graph. Kernels are
	 * also related to arcs in the repeat domain graph in \cite{zhi2006identifying}. The 
	 * fact that repeats can have shared domains is a reason for expecting kernels not to 
	 * be simple chains (paper \cite{zhi2006identifying} uses as input a set of consensus 
	 * sequences of repeats, however).
	 *
	 * Remark: the procedure does not assume $IntervalGraph.nodesArray$ to be sorted.
	 *
	 * @param out 0: number of kernels; 1(2): max. number of nodes (edges) in a kernel. 
	 */
	private static final void getKernels(int[] out) {
		final int MIN_KERNELS_PER_NODE = 5;
		int i, j;
		int top, from, lastKernel, to, nodes, maxNodes, edges, maxEdges, kernel, component;
		IntervalGraph.Node node, nodeTo;
		IntervalGraph.Edge edge;
	
		lastKernel=-1; maxNodes=0; maxEdges=0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			node.lastKernel=-1;
			node.kernels = new int[MIN_KERNELS_PER_NODE];
		}
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (!canBeKernel(node) || node.lastKernel!=-1) continue;
			kernel=++lastKernel;
			node.lastKernel=0; node.kernels[0]=kernel;
			nodes=1; edges=0;
			if (IntervalGraph.nNeighbors[i]>0) {
				IntervalGraph.stack[0]=i;
				top=0;
				while (top>=0) {
					from=IntervalGraph.stack[top--];
					for (j=0; j<IntervalGraph.nNeighbors[from]; j++) {
						edge=IntervalGraph.neighbors[from][j];
						if (!edge.on) break;
						if (!isOverlap(edge) && !isIdentity(edge)) continue;
						to=edge.getTo(from); nodeTo=IntervalGraph.nodesArray[to];
						if (!canBeKernel(nodeTo)) continue;
						edges++;
						if (nodeTo.lastKernel!=-1) continue;
						nodes++;
						nodeTo.lastKernel=0; nodeTo.kernels[0]=kernel;
						IntervalGraph.stack[++top]=to;
					}
				}
			}
			if (nodes>maxNodes) maxNodes=nodes;
			if (edges>maxEdges) maxEdges=edges;
		}
		out[0]=lastKernel+1; out[1]=maxNodes; out[2]=maxEdges;
		
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<IntervalGraph.nNodes; i++) {
				if (IntervalGraph.nodesArray[i].lastKernel==-1) continue;
				for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
					to=IntervalGraph.neighbors[i][j].getTo(i);
					if (IntervalGraph.nodesArray[to].lastKernel==-1) continue;
					if (IntervalGraph.nodesArray[i].kernels[0]!=IntervalGraph.nodesArray[to].kernels[0] && (isOverlap(IntervalGraph.neighbors[i][j]) || isIdentity(IntervalGraph.neighbors[i][j]))) {
						System.err.println("getKernels> ERROR: the following kernels are connected by an overlap or identity edge: "+IntervalGraph.nodesArray[i].kernels[0]+","+IntervalGraph.nodesArray[to].kernels[0]);
						System.exit(1);
					}
				}
			}
		}
	}
	
	
	public static final boolean canBeKernel(IntervalGraph.Node node) {
		if (node.isContained || node.lastComponent!=0 /*For discarding contained-in-short-period intervals: || (node.isContainedInterval&(1<<Constants.INTERVAL_PERIODIC))!=0*/) return false;
		return true;
	}
	
	
	/**
	 * Estimates the number of distinct repeats that have been collapsed onto the same
	 * kernel. This is done by building a bidirected graph on the kernel, and then by 
	 * building a number of DAGs from it and counting their paths. See procedures 
	 * $getKernels()$ and $BidirectedGraph.estimateAssemblies()$ for details.
	 *
	 * Remark: the resulting bidirected graph is not necessarily connected (see 
	 * $buildBidirectedGraph_mergeIdentical()$). The procedure counts paths in every 
	 * connected component (including components with just one node).
	 *
	 * Remark: the procedure considers a kernel as an assembly graph, which can contain 
	 * more that one source-sink path, i.e. it can represent more than one distinct 
	 * repeat, since distinct repeats can share substrings, and therefore they can be 
	 * collapsed onto the same kernel. E.g. if both repeats X and Y contain repeat Z,
	 * and a maximal interval that is neither left- nor right-maximal is connected by 
	 * suffix-prefix overlaps to the maximal intervals of X, Y and Z. Note, however, that 
	 * if a set of repeats share substrings, then it is not necessarily true that they are
	 * collapsed into the same kernel.
	 *
	 * Remark: at the end of the procedure, the field $bidirectedGraphNode$ of each node 
	 * in $nodesInKernel$ contains the ID of the bidirected graph node it was projected to
	 * (-1 if it was not projected to any bidirected graph node). If all nodes of a kernel
	 * have $bidirectedGraphNode=-1$, the kernel is assumed to be deleted, i.e. the 
	 * procedure can delete entire kernels.
	 *
	 * Remark: the procedure resets the $kernels$ array of every interval graph node $v$ 
	 * in the kernel, to the set of IDs of all assembly paths in the bidirected graph that 
	 * use node $v.bidirectedGraphNode$, if any.
	 *
	 * Remark: the procedure assumes $getKernels()$ to have already been executed.
	 *
	 * Remark: this procedure could be parallelized by working on each kernel 
	 * independently.
	 * 
	 * @param out output array: 0=number of deleted kernels; 1=number of distinct new 
	 * kernel tags assigned; 2=number of cyclic kernels created;
	 * @return an array of counts per old kernel; column 0=number of assemblies (-1 if the 
	 * underlying reachability DAGs contain a cycle; -2 if the kernel is deleted); 
	 * 1=min string length of an assembly; 2=max string length of an assembly; 
	 * 3=number of transitive edges; 4=number of one-end nodes; 5=number of nodes with no
	 * neighbors; 6=number of connected components of the bidirected graph.
	 */
	private static final int[][] countRepeatsPerKernel(int[] out) throws IOException {
		final int MAX_ASSEMBLIES_PER_COMPONENT = 32;  // Arbitrary
		int i, j;
		int first, last, oldNode, nNewNodes, maxNodes, maxEdges, maxDegree;
		int nComponents, nDeletedKernels, firstNewKernel, nCyclic;
		IntervalGraph.Node tmpNode;
		int[] lastOutNeighbor, nInNeighbors, component, componentLengths, componentStringLengths, nEdges, nPaths;
		IntervalGraph.Node[] tmpNodes, oldNodes_pointers, newNodes_pointers;
		int[][] kernelsSize;
		int[][] oldNodes, node2oldNode, newNodes, node2newNode, oldNode2newNode, outNeighbors, orientations, pathLengths, componentPaths;
		
		// Collecting kernel nodes
		kernelsSize = new int[nKernels][2];
		nEdges = new int[nKernels];
		getKernelsSize(kernelsSize,nEdges,tmpArray1);
		maxNodes=tmpArray1[0]; maxEdges=tmpArray1[1]; maxDegree=tmpArray1[2];
		nNodesInKernel=0;
		for (i=0; i<nKernels; i++) nNodesInKernel+=kernelsSize[i][0];
		nodesInKernel = new int[nNodesInKernel];
		sortNodesByKernel(kernelsSize);
		
		// Allocating temporary space
		kernelStats = new int[nKernels][7];
		BidirectedGraph.allocateMemory(maxNodes,maxEdges,maxDegree);
		oldNodes = new int[maxNodes][4];
		node2oldNode = new int[maxNodes][3];
		tmpNodes = new IntervalGraph.Node[maxNodes];
		newNodes = new int[maxNodes][5];
		oldNode2newNode = new int[maxNodes][2];
		node2newNode = new int[maxNodes][2];
		outNeighbors = new int[maxNodes][maxDegree];
		lastOutNeighbor = new int[maxNodes];
		orientations = new int[maxNodes][maxDegree];
		nInNeighbors = new int[maxNodes];
		component = new int[maxNodes];
		componentPaths = new int[maxNodes][10/*Arbitrary*/];
		componentLengths = new int[maxNodes];
		componentStringLengths = new int[maxNodes];
		oldNodes_pointers = new IntervalGraph.Node[maxNodes];
		newNodes_pointers = new IntervalGraph.Node[maxNodes];
		nPaths = new int[maxNodes<<1];
		pathLengths = new int[maxNodes<<1][4];
		if (tmpBoolean1==null || tmpBoolean1.length<maxNodes) tmpBoolean1 = new boolean[maxNodes];
		if (tmpByte1==null || tmpByte1.length<maxNodes) tmpByte1 = new byte[maxNodes];
		
		// Building bidirected graphs
		first=0; nDeletedKernels=0; firstNewKernel=0; nCyclic=0;
		for (i=0; i<nKernels; i++) {
			last=first+kernelsSize[i][0]-1;
			if (last==first) {
				tmpNode=IntervalGraph.nodesArray[nodesInKernel[first]];
				kernelStats[i][0]=1;
				kernelStats[i][1]=tmpNode.length();
				kernelStats[i][2]=kernelStats[i][1];
				kernelStats[i][3]=0;
				kernelStats[i][4]=0;
				kernelStats[i][5]=1;
				kernelStats[i][6]=1;
				newNodes[0][0]=tmpNode.read;
				newNodes[0][1]=tmpNode.start;
				newNodes[0][2]=tmpNode.end;
				newNodes[0][3]=tmpNode.type;
				tmpNode.bidirectedGraphNode=0;
				tmpNode.bidirectedGraphNodeOrientation=0;
				// Resetting kernel tag
				tmpNode.lastKernel=0; tmpNode.kernels[0]=firstNewKernel; 
				tmpNode.kernelOrientations = new byte[] {0};
				addPathKernelLength(tmpNode.length(),firstNewKernel,i);
				tmpNode.lastPathWithStart=0;
				tmpNode.pathsWithStart = new int[1];
				tmpNode.pathsWithStart[0]=firstNewKernel;
				tmpNode.lastPathWithEnd=0;
				tmpNode.pathsWithEnd = new int[1];
				tmpNode.pathsWithEnd[0]=-1-firstNewKernel;
				firstNewKernel++;
				// Saving bidirected graph and path labels to disk
				BidirectedGraph.clear(); BidirectedGraph.nNodes=1; BidirectedGraph.lastNeighbor[0]=-1; 
				BidirectedGraph.nodeLength[0]=tmpNode.length(); BidirectedGraph.intervalGraphPointers[0]=tmpNode;
				BidirectedGraph.printAllPaths_singleton(firstNewKernel);
				BidirectedGraph.serialize(GRAPH_FILE+"-kernel"+i+".bdgraph");
				printBidirectedGraphLabels(newNodes,1,GRAPH_FILE+"-kernel"+i+"-"+IO.ONE_NODE_ONLY_LABEL+"-"+IO.LABEL_COORDINATES_LABEL+".txt");
				BidirectedGraph.toDot(GRAPH_FILE+"-kernel"+i+".dot",newNodes);
			}
			else {
				nNewNodes=buildBidirectedGraph(first,last,oldNodes,oldNodes_pointers,node2oldNode,tmpNodes,newNodes,newNodes_pointers,oldNode2newNode,outNeighbors,lastOutNeighbor,orientations,nInNeighbors,component);				
				if (IO.CONSISTENCY_CHECKS) {
					System.err.print("IntervalGraphStep3.countRepeatsPerKernel> Consistency checks started... ");
					BidirectedGraph.checkConsistency();
					System.err.println("OK");
				}
				if (nNewNodes==0) {
					kernelStats[i][0]=-2;
					kernelStats[i][1]=0;
					kernelStats[i][2]=0;
					kernelStats[i][3]=0;
					kernelStats[i][4]=0;
					kernelStats[i][5]=0;
					kernelStats[i][6]=0;
					nDeletedKernels++;
					// Resetting kernel tags in all nodes
					for (j=first; j<=last; j++) {
						tmpNode=IntervalGraph.nodesArray[nodesInKernel[j]];
						tmpNode.lastKernel=-1;
						tmpNode.lastPathWithStart=-1;
						tmpNode.pathsWithStart=null;
						tmpNode.lastPathWithEnd=-1;
						tmpNode.pathsWithEnd=null;
					}
				}
				else {
					nComponents=BidirectedGraph.getConnectedComponents(component);
					kernelStats[i][6]=nComponents;
					BidirectedGraph.estimateAssemblies(nPaths,pathLengths,nComponents,component,componentPaths,componentLengths,componentStringLengths,kernelStats[i]);
					if (kernelStats[i][0]!=-1) {  // No cycle
						if (kernelStats[i][0]<=MAX_ASSEMBLIES_PER_COMPONENT*nComponents) {		
							// Resetting kernel tags in all nodes
							BidirectedGraph.printAllPaths(firstNewKernel,false,null);
							for (j=first; j<=last; j++) {
								tmpNode=IntervalGraph.nodesArray[nodesInKernel[j]];
								if (tmpNode.bidirectedGraphNode==-1) {
									tmpNode.lastKernel=-1;
									tmpNode.lastPathWithStart=-1;
									tmpNode.pathsWithStart=null;
									tmpNode.lastPathWithEnd=-1;
									tmpNode.pathsWithEnd=null;
								}
								else {
									tmpNode.lastKernel=BidirectedGraph.node2paths_last[tmpNode.bidirectedGraphNode];								
									if (tmpNode.kernels==null || tmpNode.kernels.length<tmpNode.lastKernel+1) tmpNode.kernels = new int[tmpNode.lastKernel+1];
									System.arraycopy(BidirectedGraph.node2paths[tmpNode.bidirectedGraphNode],0,tmpNode.kernels,0,tmpNode.lastKernel+1);
									setPathsWithEnd(first,j,firstNewKernel,node2oldNode,oldNodes,oldNode2newNode);
									setKernelOrientations(first,j,firstNewKernel);
								}
							}
							addPathKernelLength(-1,firstNewKernel,i);
							firstNewKernel+=BidirectedGraph.lastAssemblyPath+1;
							// Saving bidirected graph and path labels to disk
							BidirectedGraph.serialize(GRAPH_FILE+"-kernel"+i+".bdgraph");
							printBidirectedGraphLabels(newNodes,nNewNodes,GRAPH_FILE+"-kernel"+i+(nNewNodes==1?"-"+IO.ONE_NODE_ONLY_LABEL:"")+"-"+IO.LABEL_COORDINATES_LABEL+".txt");
							BidirectedGraph.toDot(GRAPH_FILE+"-kernel"+i+".dot",newNodes);
						}
						else {
							// Too many assembly paths make the following steps too slow.
							// Forcing each connected component of the bidirected graph
							// to have just one kernel, which is a longest assembly path.
							kernelStats[i][0]=nComponents;
							BidirectedGraph.printLongestPaths(firstNewKernel,nComponents,componentPaths,componentLengths,componentStringLengths,component);
							BidirectedGraph.propagateOrientationFromLongestPath(tmpByte1,tmpBoolean1);
							for (j=first; j<=last; j++) {
								tmpNode=IntervalGraph.nodesArray[nodesInKernel[j]];
								if (tmpNode.bidirectedGraphNode==-1) {
									tmpNode.lastKernel=-1;
									tmpNode.lastPathWithStart=-1;
									tmpNode.pathsWithStart=null;
									tmpNode.lastPathWithEnd=-1;
									tmpNode.pathsWithEnd=null;
								}
								else {
									tmpNode.lastKernel=0;
									if (tmpNode.kernels==null || tmpNode.kernels.length<1) tmpNode.kernels = new int[1];
									tmpNode.kernels[0]=firstNewKernel+component[tmpNode.bidirectedGraphNode];
									if (tmpNode.kernelOrientations==null || tmpNode.kernelOrientations.length<1) tmpNode.kernelOrientations = new byte[1];
									switch (tmpNode.bidirectedGraphNodeOrientation) {
										case 0: tmpNode.kernelOrientations[0]=tmpByte1[tmpNode.bidirectedGraphNode]; break;
										case 1: tmpNode.kernelOrientations[0]=oppositeOrientation[tmpByte1[tmpNode.bidirectedGraphNode]]; break;
										case 2: tmpNode.kernelOrientations[0]=2; break;
									}
									setPathsWithEnd(first,j,firstNewKernel,node2oldNode,oldNodes,oldNode2newNode);
								}
							}
							addPathKernelLength(-1,firstNewKernel,i);
							firstNewKernel+=BidirectedGraph.lastAssemblyPath+1;
							// Saving bidirected graph and path labels to disk
							BidirectedGraph.serialize(GRAPH_FILE+"-kernel"+i+".bdgraph");
							printBidirectedGraphLabels(newNodes,nNewNodes,GRAPH_FILE+"-kernel"+i+"-"+IO.TOO_MANY_PATHS_LABEL+"-"+IO.LABEL_COORDINATES_LABEL+".txt");
							BidirectedGraph.toDot(GRAPH_FILE+"-kernel"+i+".dot",newNodes);
						}
					}
					else {  // At least one cycle (this should be the case for kernels
						// with mostly periodic nodes). Resetting kernel tags in all nodes
						for (j=first; j<=last; j++) {
							tmpNode=IntervalGraph.nodesArray[nodesInKernel[j]];
							if (tmpNode.bidirectedGraphNode==-1) tmpNode.lastKernel=-1;
							else {
								tmpNode.lastKernel=0;
								tmpNode.kernels[0]=firstNewKernel;
								tmpNode.kernelOrientations = new byte[] {tmpNode.bidirectedGraphNodeOrientation};
							}
							tmpNode.lastPathWithStart=-1;
							tmpNode.pathsWithStart=null;
							tmpNode.lastPathWithEnd=-1;
							tmpNode.pathsWithEnd=null;
						}
						addPathKernelLength(0,firstNewKernel,i);
						firstNewKernel++; nCyclic++;				
						// Saving bidirected graph and path labels to disk
						BidirectedGraph.serialize(GRAPH_FILE+"-kernel"+i+"-cyclic.bdgraph");
						printBidirectedGraphLabels(newNodes,nNewNodes,GRAPH_FILE+"-kernel"+i+(nNewNodes==1?"-"+IO.ONE_NODE_ONLY_LABEL:"")+"-"+IO.LABEL_COORDINATES_LABEL+"-"+IO.CYCLIC_LABEL+".txt");
						BidirectedGraph.toDot(GRAPH_FILE+"-kernel"+i+"-"+IO.CYCLIC_LABEL+".dot",newNodes);
					}
				}
			}
			first=last+1;
		}
		BidirectedGraph.deallocateMemory();
		if (IO.CONSISTENCY_CHECKS && nDeletedKernels==nKernels) {
			System.err.println("countRepeatsPerKernel> ERROR: no kernel built for this cluster?!");
			System.exit(1);
		}
		
		// Building $nKernelNodes$.
		nKernelNodes = new int[firstNewKernel];
		Math.set(nKernelNodes,nKernelNodes.length-1,0);
		for (i=0; i<nNodesInKernel; i++) {
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[i]];
			for (j=0; j<=tmpNode.lastKernel; j++) nKernelNodes[tmpNode.kernels[j]]++;
		}
				
		out[0]=nDeletedKernels; out[1]=firstNewKernel; out[2]=nCyclic;
		return kernelStats;
	}
	
	
	/**
	 * Adds to $pathKernelLengths[from..]$ a single new length (if $length>=0$), or all 
	 * lengths of assembly paths in $BidirectedGraph$ (if $length < 0$).
	 * We explicitly allow $length=0$ to be stored in $pathKernelLengths$: this marks
	 * cyclic kernels.
	 * The procedure adds $kernel$ to the corresponding entries of $pathKernel2Kernel$.
	 */
	private static final void addPathKernelLength(int length, int from, int kernel) {
		final int GROWTH_RATE = 16;  // Arbitrary
		
		if (length>=0) {
			if (pathKernelLengths==null) {
				pathKernelLengths = new int[from+GROWTH_RATE];
				pathKernel2Kernel = new int[from+GROWTH_RATE];
			}
			else if (pathKernelLengths.length<=from) {
				int[] newArray = new int[from+GROWTH_RATE];
				System.arraycopy(pathKernelLengths,0,newArray,0,pathKernelLengths.length);
				pathKernelLengths=newArray;
				newArray = new int[from+GROWTH_RATE];
				System.arraycopy(pathKernel2Kernel,0,newArray,0,pathKernel2Kernel.length);
				pathKernel2Kernel=newArray;
			}
			pathKernelLengths[from]=length;
			pathKernel2Kernel[from]=kernel;
		}
		else {
			if (pathKernelLengths==null) {
				pathKernelLengths = new int[from+BidirectedGraph.lastAssemblyPath+1+GROWTH_RATE];
				pathKernel2Kernel = new int[from+BidirectedGraph.lastAssemblyPath+1+GROWTH_RATE];
			}
			else if (pathKernelLengths.length<=from+BidirectedGraph.lastAssemblyPath) {
				int[] newArray = new int[from+BidirectedGraph.lastAssemblyPath+1+GROWTH_RATE];
				System.arraycopy(pathKernelLengths,0,newArray,0,pathKernelLengths.length);
				pathKernelLengths=newArray;
				newArray = new int[from+BidirectedGraph.lastAssemblyPath+1+GROWTH_RATE];
				System.arraycopy(pathKernel2Kernel,0,newArray,0,pathKernel2Kernel.length);
				pathKernel2Kernel=newArray;
			}
			System.arraycopy(BidirectedGraph.assemblyPathStringLengths,0,pathKernelLengths,from,BidirectedGraph.lastAssemblyPath+1);
			Math.set(pathKernel2Kernel,from,from+BidirectedGraph.lastAssemblyPath,kernel);
		}
	}
	
	
	/**
	 * Remark: the procedure assumes $tmpPoints$ to be of size at least $nKernels$.
	 * Remark: the procedure works just on $kernelStats$, and it does not look at the 
	 * kernel tags of nodes.
	 */
	private static final void printKernelHistograms() {
		int i;
		int nCyclic, nDeleted;
		
		// Assemblies per kernel
		lastTmpPoint=-1; nCyclic=0; nDeleted=0;
		for (i=0; i<nKernels; i++) {
			if (kernelStats[i][0]==-1) {
				nCyclic++;
				continue;
			}
			if (kernelStats[i][0]==-2) {
				nDeleted++;
				continue;
			}
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=kernelStats[i][0];
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		System.out.println("Histogram> Assemblies per kernel:  ("+nKernels+" kernels, "+nCyclic+" cyclic, "+nDeleted+" deleted)");
		for (i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i].position+","+tmpPoints[i].mass);
		
		// Components per kernel
		lastTmpPoint=-1;
		for (i=0; i<nKernels; i++) {
			if (kernelStats[i][0]==-2) continue;
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=kernelStats[i][6];
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		System.out.println("Histogram> Connected components per kernel (including cyclic):");
		for (i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i].position+","+tmpPoints[i].mass);
		
		// Max string length
		lastTmpPoint=-1;
		for (i=0; i<nKernels; i++) {
			if (kernelStats[i][0]==-1 || kernelStats[i][0]==-2) continue;
			if (kernelStats[i][1]>kernelStats[i][2]) {
				System.err.println("ERROR in kernel "+i+": minLength="+kernelStats[i][1]+" maxLength="+kernelStats[i][2]);
				System.exit(1);
			}
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=kernelStats[i][2];
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		System.out.println("Histogram> Max assembly length per kernel:");
		for (i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i].position+","+tmpPoints[i].mass);
		
		// One-end nodes
		lastTmpPoint=-1;
		for (i=0; i<nKernels; i++) {
			if (kernelStats[i][0]==-1 || kernelStats[i][0]==-2) continue;
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=kernelStats[i][4];
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		System.out.println("Histogram> One-end nodes per kernel:");
		for (i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i].position+","+tmpPoints[i].mass);
		
		// Nodes with no neighbor
		lastTmpPoint=-1;
		for (i=0; i<nKernels; i++) {
			if (kernelStats[i][0]==-1 || kernelStats[i][0]==-2) continue;
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=kernelStats[i][5];
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		System.out.println("Histogram> Sigleton nodes per kernel:");
		for (i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i].position+","+tmpPoints[i].mass);
		
		// Transitive eges
		lastTmpPoint=-1;
		for (i=0; i<nKernels; i++) {
			if (kernelStats[i][0]==-1 || kernelStats[i][0]==-2) continue;
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=kernelStats[i][3];
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		System.out.println("Histogram> Transitive edges per kernel:");
		for (i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i].position+","+tmpPoints[i].mass);
	}
	
	
	/**
	 * Prints the coordinates for retrieving the string labels of bidirected graph nodes
	 * from the reads database.
	 *
	 * @param newNodes 0=readID, 1=start, 2=end.
	 */
	private static final void printBidirectedGraphLabels(int[][] newNodes, int nNodes, String path) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(path),IO.BUFFER_SIZE);
		for (int i=0; i<nNodes; i++) bw.write(newNodes[i][0]+","+newNodes[i][1]+","+newNodes[i][2]+"\n");
		bw.close();
	}
	
	
	/**
	 * Stores in $size[i][0]$ the number of nodes in kernel $i$ and in $nEdges[i]$ the 
	 * number of edges in kernel $i$.
	 *
	 * @param output array: 0=max number of nodes in a kernel; 1=max number of edges in a
	 * kernel; 2=max degree of a node, counted just inside its kernel.
	 */
	private static final void getKernelsSize(int[][] size, int[] nEdges, int[] out) {
		int i, j, k;
		int kernel, degree, maxNodes, maxEdges, maxDegree;
		IntervalGraph.Node node;
		IntervalGraph.Edge edge;
		
		for (i=0; i<nKernels; i++) size[i][0]=0;
		Math.set(nEdges,0,nKernels-1,0);
		maxDegree=0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			if (IntervalGraph.nodesArray[i].lastKernel==-1) continue;
			kernel=IntervalGraph.nodesArray[i].kernels[0];
			size[kernel][0]++;
			degree=0;
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
				edge=IntervalGraph.neighbors[i][j];
				if (!edge.on) break;
				node=IntervalGraph.nodesArray[edge.getTo(i)];
				for (k=0; k<=node.lastKernel; k++) {
					if (node.kernels[k]==kernel) {
						degree++;
						nEdges[kernel]++;
						break;
					}
				}
			}
			if (degree>maxDegree) maxDegree=degree;
		}
		for (i=0; i<nKernels; i++) nEdges[i]>>=1;
		
		// Building the output
		maxNodes=0; maxEdges=0;
		for (i=0; i<nKernels; i++) {
			if (size[i][0]>maxNodes) maxNodes=size[i][0];
			if (nEdges[i]>maxEdges) maxEdges=nEdges[i];
		}
		out[0]=maxNodes; out[1]=maxEdges; out[2]=maxDegree;
	}
	
	
	/**
	 * Collects in $nodesInKernel[0..nNodesInKernel-1]$ the IDs of all nodes that belong 
	 * to a kernel, sorted by kernel ID, then by nodeID.
	 *
	 * Remark: the procedure assumes $kernelSize$ to have at least two columns, only the
	 * first of which is used by the input.
	 */
	private static final void sortNodesByKernel(int[][] kernelsSize) {
		int i, kernel, lastNodeInKernel;
		
		// Sorting by kernel
		kernelsSize[0][1]=0;
		for (i=1; i<nKernels; i++) kernelsSize[i][1]=kernelsSize[i-1][1]+kernelsSize[i-1][0];
		for (i=0; i<IntervalGraph.nNodes; i++) {
			if (IntervalGraph.nodesArray[i].lastKernel==-1) continue;
			kernel=IntervalGraph.nodesArray[i].kernels[0];
			nodesInKernel[kernelsSize[kernel][1]++]=i;
		}
		lastNodeInKernel=kernelsSize[nKernels-1][1]-1;
		
		// Sorting by node ID
		kernelsSize[0][1]=0;
		for (i=1; i<nKernels; i++) kernelsSize[i][1]=kernelsSize[i-1][1]+kernelsSize[i-1][0];
		for (i=0; i<nKernels-1; i++) {
			if (kernelsSize[i+1][1]-kernelsSize[i][1]>1) Arrays.sort(nodesInKernel,kernelsSize[i][1],kernelsSize[i+1][1]);
		}
		if (lastNodeInKernel+1-kernelsSize[nKernels-1][1]>1) Arrays.sort(nodesInKernel,kernelsSize[nKernels-1][1],lastNodeInKernel+1);
	}
	
	
	/**
	 * Builds the bidirected graph of a kernel, starting from overlap edges between nodes 
	 * in the same kernel, and then simplifying the topology using supplement edges and 
	 * positions of intervals on the reads. The resulting bidirected graph can be 
	 * significantly simpler than the initial kernel graph, and it should be similar to a 
	 * simple path that encodes a full copy of a repeat. The observation that the assembly
	 * graph of a repeat should be simple was already present in 
	 * \cite{quitzau2008detecting}, which however assumed that biological repeats are not 
	 * likely to contain other repeats. We don't make such assumption, but we just assume 
	 * that repeats inside other repeats should form distinct assembly graphs, and that 
	 * $Alignments.minAlignmentLength$ is big enough to render small repeat modules 
	 * invisible (it is well-known that biological repeats have internal modules shared by
	 * other repeats).
	 *
	 * Even after such simplifications, the bidirected graph might still contain "forks", 
	 * i.e. nodes that overlap, on the same side, with many other nodes. Forks can either
	 * be caused by different repeats sharing substrings and being collapsed to the same
	 * assembly graph, or by undersplits, i.e. by the fact that a repeat contained in the 
	 * interval of the fork, and also contained in other intervals, has not been correctly
	 * factorized, nor it has been detected and corrected by the previous steps of the
	 * pipeline.
	 *
	 * Remark: nodes that are contained in others are more likely to be overlap forks,
	 * since they might be substrings that are shared by different repeats. So, choosing
	 * kernel nodes among those that are maximal by containment should reduce forks.
	 *
	 * Remark: the procedure could use the information of whether interval graph nodes are
	 * RC-palindromes to add edges to the bidirected graph that are not supported by an 
	 * alignment. The current version of the procedure does not do that.
	 *
	 * Remark: at the end of the procedure, the field $bidirectedGraphNode$ of each node 
	 * in $nodesInKernel$ contains the ID of the bidirected graph node it was projected to
	 * (-1 if it was not projected to any bidirected graph node). It can happen that all
	 * nodes in $nodesInKernel[first..last]$ have $bidirectedGraphNode=-1$.
	 *
	 * @param first,last $nodesInKernel[first..last]$ contains all and only the nodes in
	 * a kernel;
	 * @param newNodes temporary and output matrix: see $buildBidirectedGraph_
	 * mergeIdentical()$ for the required size; at the end of the procedure, this matrix
	 * stores the read, startPosition, endPosition, type, of every node in the bidirected 
	 * graph (such coordinates might not correspond to any node of the interval graph);
	 * the length of the matrix must be $>=last-first+1$; type is -2 if interval graph 
	 * nodes of incompatible type were merged in the same bidirected graph node;
	 * @param other variables: temporary space: see $buildBidirectedGraph_
	 * mergeStraddling()$ and $buildBidirectedGraph_mergeIdentical()$ for details;
	 * @return number of nodes in the bidirected graph; this can be zero, if all nodes 
	 * were simplified away.
	 */
	private static final int buildBidirectedGraph(int first, int last, int[][] oldNodes, IntervalGraph.Node[] oldNodes_pointers, int[][] node2oldNode, IntervalGraph.Node[] tmpNodes, int[][] newNodes, IntervalGraph.Node[] newNodes_pointers, int[][] oldNode2newNode, int[][] outNeighbors, int[] lastOutNeighbor, int[][] orientations, int[] nInNeighbors, int[] component) {
		boolean overlapAlignments;
		int i, j, k;
		int from, to, toPrime, node1, node2, nOldNodes, nNewNodes, result, tmp;
		int oldNode, newNode, send, orientation;
		final int kernel = IntervalGraph.nodesArray[nodesInKernel[first]].kernels[0];
		IntervalGraph.Node intervalNode;
		IntervalGraph.Edge intervalEdge;
		BidirectedGraph.BidirectedEdge bidirectedEdge;
		
		// Raw topology (built using also supplement edges).
		BidirectedGraph.clear();
		BidirectedGraph.nNodes=last-first+1;
		for (i=first; i<=last; i++) {
			BidirectedGraph.nodeLength[i-first]=IntervalGraph.nodesArray[nodesInKernel[i]].length();
			BidirectedGraph.intervalGraphPointers[i-first]=IntervalGraph.nodesArray[nodesInKernel[i]];
		}
		for (i=first; i<=last; i++) {
			from=nodesInKernel[i];
			for (j=0; j<IntervalGraph.nNeighbors[from]; j++) {
				intervalEdge=IntervalGraph.neighbors[from][j];
				intervalEdge.printed=-1;  // Using the $printed$ field of each interval graph edge to mark whether an edge has already been used by this procedure.
			}
		}
		for (i=first; i<=last; i++) {
			from=nodesInKernel[i];
			for (j=0; j<IntervalGraph.nNeighbors[from]; j++) {
				intervalEdge=IntervalGraph.neighbors[from][j];
				if (intervalEdge.printed==1) continue;
				to=intervalEdge.getTo(from);
				toPrime=Arrays.binarySearch(nodesInKernel,first,last+1,to);
				// Remark: $toPrime$ is necessary. Just comparing $IntervalGraph.
				// nodesArray[to].kernels[0]$ to $kernel$ might yield wrong nodes, since
				// the previous procedures might have reset to $kernel$ the kernel IDs of
				// some nodes that were not marked with $kernel$ before the reassignment. 
				if (toPrime<0 || !isOverlap(intervalEdge) || IntervalGraph.nodesArray[to].lastKernel==-1) continue;
				if (from==intervalEdge.nodeID1) {
					node1=i-first; 
					node2=toPrime-first;
				}
				else {
					node1=toPrime-first;
					node2=i-first;
				}
				overlapAlignments=intervalEdge.nAlignmentsOfType(1)!=Math.POSITIVE_INFINITY;
				for (k=Constants.OVERLAP_PREFIX_PREFIX; k<=Constants.OVERLAP_SUFFIX_SUFFIX; k<<=1) {
					if ((intervalEdge.overlap&k)==0) continue;
					bidirectedEdge=BidirectedGraph.getEdge();
					bidirectedEdge.node1=node1; 
					bidirectedEdge.node2=node2; 
					bidirectedEdge.type=k;
					bidirectedEdge.setOverhangs(k,intervalEdge.overhangs);
					bidirectedEdge.nAlignments=overlapAlignments?1:2;  // Assigning 2 alignments is arbitrary.
					BidirectedGraph.addEdge(bidirectedEdge);
				}
				intervalEdge.printed=1;
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			System.err.print("buildBidirectedGraph> Consistency checks (1) started... ");
			BidirectedGraph.checkConsistency();
			System.err.println("OK");
		}
		
		// Simplifications of the topology.
		// $buildBidirectedGraph_disconnectInsertions()$ should be done before
		// $buildBidirectedGraph_fixForks()$, otherwise all extensions of a fork on one
		// side might disappear.
		if (BidirectedGraph.nNodes>1) {
			result=buildBidirectedGraph_addOverlaps(first,last,NEW_OVERLAPS_THRESHOLD);
			if (IO.SHOW_STD_ERR_PRIME) System.err.println("buildBidirectedGraph> new overlaps added: "+result);
		}
		result=buildBidirectedGraph_disconnectInsertions(first,last);
		if (IO.SHOW_STD_ERR_PRIME) System.err.println("buildBidirectedGraph> nodes disconnected by insertion: "+result);
		if (BidirectedGraph.nRemoved==BidirectedGraph.nNodes) {
			for (i=first; i<=last; i++) {
				intervalNode=IntervalGraph.nodesArray[nodesInKernel[i]];
				intervalNode.bidirectedGraphNode=-1;
			}
			BidirectedGraph.resetEdgePool();  // Reuse all edges in the next call
			return 0;
		}
		else if (BidirectedGraph.nRemoved==BidirectedGraph.nNodes-1) {
			for (i=first; i<=last; i++) {
				intervalNode=IntervalGraph.nodesArray[nodesInKernel[i]];
				intervalNode.bidirectedGraphNode=-1;
			}
			BidirectedGraph.deleteRemoved();
			BidirectedGraph.resetEdgePool();
			intervalNode=BidirectedGraph.intervalGraphPointers[0];
			intervalNode.bidirectedGraphNode=0;
			intervalNode.bidirectedGraphNodeOrientation=0;
			oldNodes[0][0]=intervalNode.read;
			oldNodes[0][1]=intervalNode.start;
			oldNodes[0][2]=intervalNode.end;
			oldNodes[0][3]=intervalNode.type;
			for (i=first; i<=last; i++) {
				if (IntervalGraph.nodesArray[nodesInKernel[i]]==intervalNode) {
					node2oldNode[i-first][0]=0;
					node2oldNode[i-first][1]=0;
					node2oldNode[i-first][2]=intervalNode.length()-1;
					break;
				}
			}
			System.arraycopy(oldNodes[0],0,newNodes[0],0,4);
			oldNode2newNode[0][0]=0;
			oldNode2newNode[0][1]=0;
			return 1;
		}		
		if (BidirectedGraph.nNodes-BidirectedGraph.nRemoved>2) {
			buildBidirectedGraph_fixForks(first,last,tmpIO);
			if (IO.SHOW_STD_ERR_PRIME) System.err.println("buildBidirectedGraph> "+tmpIO[0]+" nodes and "+tmpIO[1]+" edges removed by fixing forks");
		}
		if (BidirectedGraph.nRemoved==BidirectedGraph.nNodes-1) {
			for (i=first; i<=last; i++) {
				intervalNode=IntervalGraph.nodesArray[nodesInKernel[i]];
				intervalNode.bidirectedGraphNode=-1;
			}
			BidirectedGraph.deleteRemoved();
			BidirectedGraph.resetEdgePool();
			intervalNode=BidirectedGraph.intervalGraphPointers[0];
			intervalNode.bidirectedGraphNode=0;
			intervalNode.bidirectedGraphNodeOrientation=0;
			oldNodes[0][0]=intervalNode.read;
			oldNodes[0][1]=intervalNode.start;
			oldNodes[0][2]=intervalNode.end;
			oldNodes[0][3]=intervalNode.type;
			for (i=first; i<=last; i++) {
				if (IntervalGraph.nodesArray[nodesInKernel[i]]==intervalNode) {
					node2oldNode[i-first][0]=0;
					node2oldNode[i-first][1]=0;
					node2oldNode[i-first][2]=intervalNode.length()-1;
					break;
				}
			}
			System.arraycopy(oldNodes[0],0,newNodes[0],0,4);
			oldNode2newNode[0][0]=0;
			oldNode2newNode[0][1]=0;
			return 1;
		}
		if (IO.CONSISTENCY_CHECKS) {
			System.err.print("buildBidirectedGraph> Consistency checks (2) started...");
			BidirectedGraph.checkConsistency();
			System.err.println("OK");
		}
		
		// Contractions
		nOldNodes=buildBidirectedGraph_mergeStraddling(first,last,oldNodes,oldNodes_pointers,node2oldNode,tmpNodes);
		if (nOldNodes==-1) nOldNodes=BidirectedGraph.nNodes;
		else if (nOldNodes==1) {
			for (i=first; i<=last; i++) {
				intervalNode=IntervalGraph.nodesArray[nodesInKernel[i]];
				intervalNode.bidirectedGraphNode=node2oldNode[i-first][0];
				intervalNode.bidirectedGraphNodeOrientation=0;
			}
			System.arraycopy(oldNodes[0],0,newNodes[0],0,4);
			oldNode2newNode[0][0]=0;
			oldNode2newNode[0][1]=0;
			BidirectedGraph.resetEdgePool();
			return 1;
		}		
		BidirectedGraph.sortNeighbors();
		nNewNodes=buildBidirectedGraph_mergeIdentical(first,last,oldNodes,nOldNodes,oldNodes_pointers,node2oldNode,newNodes,newNodes_pointers,oldNode2newNode,outNeighbors,lastOutNeighbor,orientations,nInNeighbors,component);
		if (nNewNodes==-1) nNewNodes=nOldNodes;
		else if (nNewNodes==1) {
			for (i=first; i<=last; i++) {
				intervalNode=IntervalGraph.nodesArray[nodesInKernel[i]];
				oldNode=node2oldNode[i-first][0];
				if (oldNode==-1) intervalNode.bidirectedGraphNode=-1;
				else {
					intervalNode.bidirectedGraphNode=oldNode2newNode[oldNode][0];
					intervalNode.bidirectedGraphNodeOrientation=(byte)oldNode2newNode[oldNode][1];
				}
			}
			BidirectedGraph.resetEdgePool();
			return 1;
		}
		BidirectedGraph.breakSimpleCycles();  // Applies also to periodic kernels
		BidirectedGraph.resetEdgePool();
		if (IO.CONSISTENCY_CHECKS) {
			System.err.print("buildBidirectedGraph> Consistency checks (3) started... ");
			BidirectedGraph.checkConsistency();
			System.err.println("OK");
		}	
		
		// Setting field $bidirectedGraphNode$ of interval graph nodes.
		for (i=first; i<=last; i++) {
			intervalNode=IntervalGraph.nodesArray[nodesInKernel[i]];
			oldNode=node2oldNode[i-first][0];
			if (oldNode==-1) {
				intervalNode.bidirectedGraphNode=-1;
				continue;
			}
			intervalNode.bidirectedGraphNode=oldNode2newNode[oldNode][0];
			intervalNode.bidirectedGraphNodeOrientation=(byte)oldNode2newNode[oldNode][1];
		}
		if (IO.CONSISTENCY_CHECKS) {
			System.err.print("buildBidirectedGraph> Consistency checks (4) started... ");
			BidirectedGraph.checkConsistency();
			System.err.println("OK");
		}
		
		return nNewNodes;
	}

	
	/**
	 * Sets the $pathsWith*,lastPathWith*$ variables of the $i$-th node of 
	 * $nodesInKernel$. All interval graph nodes in the same kernel lie in
	 * $nodesInKernel[first..j]$ with $j>=i$.
	 *
	 * Remark: the procedure assumes that $intervalNode.bidirectedGraphNode$ has already
	 * been set.
	 * 
	 * @param node2oldNode mapping generated by $buildBidirectedGraph_mergeStraddling()$;
	 * @param oldNode2newNode mapping produced by $buildBidirectedGraph_mergeIdentical()$.
	 */
	private static final void setPathsWithEnd(int first, int i, int firstTag, int[][] node2oldNode, int[][] oldNodes, int[][] oldNode2newNode) {
		final int THRESHOLD = (3*IO.quantum)>>1;  // To tolerate some drift when e.g. a chain of interval graph nodes connected by identity edges is merged.
		int oldNode, newNode, end, orientation, leftOffset, rightOffset;
		final IntervalGraph.Node intervalNode = IntervalGraph.nodesArray[nodesInKernel[i]];
				
		intervalNode.lastPathWithStart=-1;
		intervalNode.pathsWithStart=null;
		intervalNode.lastPathWithEnd=-1;
		intervalNode.pathsWithEnd=null;
		newNode=intervalNode.bidirectedGraphNode;
		if (newNode==-1) return;
		leftOffset=node2oldNode[i-first][1];
		oldNode=node2oldNode[i-first][0];
		rightOffset=oldNodes[oldNode][2]-(oldNodes[oldNode][1]+node2oldNode[i-first][2]);
		end=BidirectedGraph.isOneEnd(newNode);
		if (end==-2) return;
		orientation=oldNode2newNode[node2oldNode[i-first][0]][1];
		if (end==-1) {
			if (orientation==0) {
				if (leftOffset<=THRESHOLD) {
					intervalNode.lastPathWithStart=BidirectedGraph.nAssemblyPaths(newNode,firstTag,true)-1;
					if (intervalNode.lastPathWithStart>=0) {
						intervalNode.pathsWithStart = new int[intervalNode.lastPathWithStart+1];
						BidirectedGraph.getAssemblyPaths(newNode,firstTag,true,intervalNode.pathsWithStart);
					}
				}
				if (rightOffset<=THRESHOLD) {
					intervalNode.lastPathWithEnd=BidirectedGraph.nAssemblyPaths(newNode,firstTag,false)-1;
					if (intervalNode.lastPathWithEnd>=0) {
						intervalNode.pathsWithEnd = new int[intervalNode.lastPathWithEnd+1];
						BidirectedGraph.getAssemblyPaths(newNode,firstTag,false,intervalNode.pathsWithEnd);
					}
				}
			}
			else if (orientation==1) {
				if (rightOffset<=THRESHOLD) {
					intervalNode.lastPathWithStart=BidirectedGraph.nAssemblyPaths(newNode,firstTag,false)-1;
					if (intervalNode.lastPathWithStart>=0) {
						intervalNode.pathsWithStart = new int[intervalNode.lastPathWithStart+1];
						BidirectedGraph.getAssemblyPaths(newNode,firstTag,false,intervalNode.pathsWithStart);
					}
				}
				if (leftOffset<=THRESHOLD) {
					intervalNode.lastPathWithEnd=BidirectedGraph.nAssemblyPaths(newNode,firstTag,true)-1;
					if (intervalNode.lastPathWithEnd>=0) {
						intervalNode.pathsWithEnd = new int[intervalNode.lastPathWithEnd+1];
						BidirectedGraph.getAssemblyPaths(newNode,firstTag,true,intervalNode.pathsWithEnd);
					}
				}
			}
		}
		else if (end==1) {
			if (orientation==0) {
				if (leftOffset<=THRESHOLD) {
					intervalNode.lastPathWithStart=BidirectedGraph.nAssemblyPaths(newNode,firstTag,true)-1;
					if (intervalNode.lastPathWithStart>=0) {
						intervalNode.pathsWithStart = new int[intervalNode.lastPathWithStart+1];
						BidirectedGraph.getAssemblyPaths(newNode,firstTag,true,intervalNode.pathsWithStart);
					}
				}
			}
			else {
				if (leftOffset<=THRESHOLD) {
					intervalNode.lastPathWithEnd=BidirectedGraph.nAssemblyPaths(newNode,firstTag,true)-1;
					if (intervalNode.lastPathWithEnd>=0) {
						intervalNode.pathsWithEnd = new int[intervalNode.lastPathWithEnd+1];
						BidirectedGraph.getAssemblyPaths(newNode,firstTag,true,intervalNode.pathsWithEnd);
					}
				}
			}
		}
		else if (end==0) {
			if (orientation==1) {
				if (rightOffset<=THRESHOLD) {
					intervalNode.lastPathWithStart=BidirectedGraph.nAssemblyPaths(newNode,firstTag,false)-1;
					if (intervalNode.lastPathWithStart>=0) {
						intervalNode.pathsWithStart = new int[intervalNode.lastPathWithStart+1];
						BidirectedGraph.getAssemblyPaths(newNode,firstTag,false,intervalNode.pathsWithStart);
					}
				}
			}
			else {
				if (rightOffset<=THRESHOLD) {
					intervalNode.lastPathWithEnd=BidirectedGraph.nAssemblyPaths(newNode,firstTag,false)-1;
					if (intervalNode.lastPathWithEnd>=0) {
						intervalNode.pathsWithEnd = new int[intervalNode.lastPathWithEnd+1];
						BidirectedGraph.getAssemblyPaths(newNode,firstTag,false,intervalNode.pathsWithEnd);
					}
				}
			}
		}
	}
	
	
	/**
	 * Sets the $kernelOrientations$ array of the $i$-th node of $nodesInKernel$.
	 *
	 * Remark: the procedure assumes that $intervalNode.bidirectedGraphNode$ has already
	 * been set and is not -1.
	 */
	private static final void setKernelOrientations(int first, int i, int firstTag) {
		byte orientation;
		int bidirectedNode, lastKernel;
		final IntervalGraph.Node intervalNode = IntervalGraph.nodesArray[nodesInKernel[i]];
		
		intervalNode.kernelOrientations=null;
		bidirectedNode=intervalNode.bidirectedGraphNode;
		lastKernel=intervalNode.lastKernel;
		intervalNode.kernelOrientations = new byte[lastKernel+1];
		orientation=intervalNode.bidirectedGraphNodeOrientation;
		BidirectedGraph.getAssemblyPathOrientations(bidirectedNode,firstTag,intervalNode.kernelOrientations);
		if (orientation==1) {
			for (i=0; i<=lastKernel; i++) intervalNode.kernelOrientations[i]=oppositeOrientation[intervalNode.kernelOrientations[i]];
		}
		else if (orientation==2) {
			for (i=0; i<=lastKernel; i++) intervalNode.kernelOrientations[i]=2;
		}
	}

	
	/**
	 * Disconnects nodes of the kernel that insert into other nodes of the same kernel. 
	 * This is useful, since a kernel can be the union of the assembly graphs of multiple 
	 * fragments of a repeat, and insertions among nodes in the same kernel can create
	 * overlap forks (possibly on both sides of the same fork node). The procedure 
	 * considers all edges of the interval graph, including supplement edges.
	 *
	 * The procedure also removes kernel nodes that are contained into another node
	 * (which does not necessarily belong to any kernel). This is possible thanks to 
	 * supplement edges.
	 *
	 * Remark: nodes are disconnected, rather than removed, since otherwise a fragment
	 * might never belong to the output. If the kernel node that is the destination of the
	 * insertion is the concatenation of multiple repeats, and if the fragment is an 
	 * independent repeat, this would mean that the repeat would never be in the output.
	 *
	 * Remark: it is correct to assign the same kernel ID to the fragments, since this is 
	 * exactly what will be done later when detecting insertion edges across kernels.
	 *
	 * Remark: copies of the repeat, with internal deletions, are either detected as 
	 * adjacent fragments, and thus have the same effect on the kernel graph, or as a 
	 * single dense substring of substring type, so they might form a distinct kernel 
	 * graph. Having multiple internal deletions is very common for transposable elements:
	 * see e.g. \cite{quesneville2003detection}.
	 *
	 * Remark: the procedure could disconnect the remaining graph, e.g. because a 
	 * disconnected node bridges two connected components on the side in which it is not 
	 * maximal.
	 *
	 * @return the number of nodes disconnected or removed by the procedure.
	 */
	private static final int buildBidirectedGraph_disconnectInsertions(int first, int last) {
		boolean foundInsertion, foundContainment;
		int i, j;
		int from, to, kernel, nNeighbors, out;
		IntervalGraph.Node node;
		IntervalGraph.Edge intervalEdge;
		
		out=0;
		for (i=first; i<=last; i++) {
			from=nodesInKernel[i];
			node=IntervalGraph.nodesArray[from];
			if (node.type==Constants.INTERVAL_PERIODIC) continue;
			kernel=node.kernels[0];
			nNeighbors=IntervalGraph.nNeighbors[from];
			foundInsertion=false; foundContainment=false;
			for (j=0; j<nNeighbors; j++) {
				intervalEdge=IntervalGraph.neighbors[from][j];
				to=intervalEdge.getTo(from);
				node=IntervalGraph.nodesArray[to];
				if ( node.lastKernel!=-1 && node.kernels[0]==kernel && Arrays.binarySearch(nodesInKernel,first,last+1,to)>=0 &&
				     ( (intervalEdge.insertion==Constants.INSERTION_ONE_IN_TWO && intervalEdge.nodeID1==from) ||
					   (intervalEdge.insertion==Constants.INSERTION_TWO_IN_ONE && intervalEdge.nodeID2==from)
					 )
				   ) {
					foundInsertion=true;
if (IO.SHOW_INTERACTIVE) System.err.println("buildBidirectedGraph_disconnectInsertions> REMOVED NODE "+IntervalGraph.nodesArray[from]+"\n BECAUSE IT INSERTS INTO NODE "+node);
					if (foundContainment) break;
				}
				if ( (intervalEdge.containment==Constants.CONTAINMENT_ONE_IN_TWO && intervalEdge.nodeID1==from) ||
					 (intervalEdge.containment==Constants.CONTAINMENT_TWO_IN_ONE && intervalEdge.nodeID2==from)
				   ) {
   					foundContainment=true;
if (IO.SHOW_INTERACTIVE) System.err.println("buildBidirectedGraph_disconnectInsertions> REMOVED NODE "+IntervalGraph.nodesArray[from]+"\n BECAUSE IT IS CONTAINED IN NODE "+node);
   					if (foundInsertion) break;
				}
			}
			if (foundContainment) {
				BidirectedGraph.removeNode(i-first);
				out++;
			}
			else if (foundInsertion) {
				BidirectedGraph.disconnectNode(i-first);
				out++;
			}
		}		
		return out;
	}
	
	
	/**
	 * Adds overlap edges between pairs of nodes in the same kernel, using all
	 * non-supplement interval graph edges and all supplement edges (but discarding the
	 * overlap information of every edge, which should have already been used to build the
	 * kernel graph).
	 * This adds interval graph edges that are not of overlap type e.g. because of the 
	 * type or maximality of their nodes. Discarding maximality is useful for removing 
	 * overlap forks induced by repeat fragments (since a kernel graph can be the union of 
	 * the assembly graphs of many fragments of the same repeat). Adding overlap edges
	 * might make the kernel graph more complex, but it might also help simplify it later
	 * on by transitive reduction (overlap forks arise precisely because two neighbors of 
	 * the same node do not overlap).
	 *
	 * Remark: the procedure considers as candidates for a new overlap edge, all pairs of 
	 * nodes in the same kernel that are: (1) not both periodic; (2) not already connected
	 * by an overlap edge of any type.
	 *
	 * Remark: this is just a way to simplify kernel graphs. If done earlier, it would 
	 * build larger kernel graphs.
	 *
	 * Remark: it is correct to assign the same kernel ID to the fragments, since this is 
	 * exactly what will be done later when detecting insertion edges across kernels.
	 *
	 * Remark: the procedure does not assume the neighbors of a node of the bidirected 
	 * graph to be sorted.
	 *
	 * @param threshold to decide overlaps;
	 * @return the number of new overlap edges created by the procedure.
	 */
	private static final int buildBidirectedGraph_addOverlaps(int first, int last, int threshold) {
		boolean periodic1, periodic2, orientation;
		int i, j, k;
		int from, to, kernel, last1, last2, overlapType, nodeID1, nodeID2, nNeighbors, nAdded;
		int start1, end1, start2, end2, alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2; 
		IntervalGraph.Node node1, node2;
		IntervalGraph.Edge intervalEdge;
		BidirectedGraph.BidirectedEdge bidirectedEdge;
		
		nAdded=0;
		for (i=first; i<=last; i++) {
			nodeID1=i-first;
			if (BidirectedGraph.removed[nodeID1]) continue;
			from=nodesInKernel[i];
			node1=IntervalGraph.nodesArray[from]; 
			start1=node1.start; end1=node1.end; 
			periodic1=node1.type==Constants.INTERVAL_PERIODIC;
			kernel=node1.kernels[0];
			nNeighbors=IntervalGraph.nNeighbors[from];
			for (j=0; j<nNeighbors; j++) {
				intervalEdge=IntervalGraph.neighbors[from][j];
				if (intervalEdge.orientation==2 || intervalEdge.orientation==-1) continue;  // We cannot decompress the orientation of each contained alignment
				orientation=intervalEdge.orientation==0;
				to=intervalEdge.getTo(from);
				nodeID2=Arrays.binarySearch(nodesInKernel,first,last+1,to)-first;
				if (nodeID2<0 || BidirectedGraph.removed[nodeID2] || BidirectedGraph.findEdge(nodeID1,nodeID2)>=0) continue;
				node2=IntervalGraph.nodesArray[to];
				periodic2=node2.type==Constants.INTERVAL_PERIODIC;
				if (node2.lastKernel==-1 || node2.kernels[0]!=kernel || (periodic1 && periodic2)) continue;
				start2=node2.start; end2=node2.end;
				last1=intervalEdge.projectOnto(from,tmpArray1,false);
				if (last1==-1) continue;
				last2=intervalEdge.projectOnto(to,tmpArray2,false);
				if (last2==-1 || last2!=last1) continue;
				for (k=0; k<last1; k+=2) {
					// This works since $IntervalGraph.Edge.projectOnto()$ writes its
					// output in the same order for both nodes.
					alignmentStart1=start1+tmpArray1[k];
					alignmentEnd1=start1+tmpArray1[k+1];
					alignmentStart2=start2+tmpArray2[k];
					alignmentEnd2=start2+tmpArray2[k+1];		
					overlapType=IntervalGraphStep2.isOverlap_simple_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,orientation,threshold);
					if (overlapType==-1) continue;
					bidirectedEdge=BidirectedGraph.getEdge();
					bidirectedEdge.node1=nodeID1;
					bidirectedEdge.node2=nodeID2;
					bidirectedEdge.type=overlapType;
					if (overlapType==Constants.OVERLAP_PREFIX_PREFIX) {
						bidirectedEdge.overhang1=end1-alignmentEnd1;
						bidirectedEdge.overhang2=end2-alignmentEnd2;
					}
					else if (overlapType==Constants.OVERLAP_PREFIX_SUFFIX) {
						bidirectedEdge.overhang1=end1-alignmentEnd1;
						bidirectedEdge.overhang2=alignmentStart2-start2;
					}
					else if (overlapType==Constants.OVERLAP_SUFFIX_PREFIX) {
						bidirectedEdge.overhang1=alignmentStart1-start1;
						bidirectedEdge.overhang2=end2-alignmentEnd2;						
					}
					else if (overlapType==Constants.OVERLAP_SUFFIX_SUFFIX) {
						bidirectedEdge.overhang1=alignmentStart1-start1;
						bidirectedEdge.overhang2=alignmentStart2-start2;
					}
					bidirectedEdge.nAlignments=1;
					BidirectedGraph.addEdge(bidirectedEdge);
					nAdded++;
				}
			}
		}
		return nAdded;
	}
	
	
	/**
	 * Tries to fix nonperiodic forks in the bidirected graph by removing a neighbor $v$ 
	 * of a fork, whose overhang with the fork is identical to or contained in the 
	 * overhang with the fork of another neighbor $w$ of the fork (using either normal or 
	 * supplement edges), and such that: (1) $v$ and $w$ do not already overlap on the 
	 * sides of their overhangs with the fork; (2) the neighbors of $v$ on the side of its 
	 * overhang are a subset of the neighbors of $w$ on the side of its overhang
	 * (otherwise, a neighbor of $v$ could have at best a shared substring with $w$, but 
	 * they might not be similar in the remaining parts, thus the fork is justified); 
	 * (3) $v$ has only one edge on the side of its overlap with the fork.
	 * This can be seen as a transitive reduction using edges of the interval graph that 
	 * are not overlaps, as well as supplement edges. Nodes are removed just from the 
	 * bidirected graph, i.e. the purpose of the procedure is just to make kernel graphs
	 * simpler.
	 *
	 * Remark: the procedure does not disconnect the bidirected graph. In particular, if
	 * two forks share two neighbors, the procedure does not accidentally delete both 
	 * neighbors.
	 *
	 * Remark: given two neighbors of a fork $v,w$, such that the overhang of $v$ with the
	 * fork if equal to or contained in the overhang of $w$ with the fork, the procedure
	 * does not check whether a neighbor $z$ of $v$ is \emph{reachable} from $w$. This 
	 * should not be necessary, since $z$ should directly overlap $w$.
	 *
	 * Remark: the procedure does not assume the neighbors of a node of the bidirected 
	 * graph to be sorted.
	 *
	 * Remark: the procedure uses $tmpArray4$ and $tmpNodes$.
	 *
	 * @param out temporary space with at least 4 elements; at the end of the procedure, 
	 * stores the total number of bidirected graph nodes (cell 0) and edges (cell 1) 
	 * removed.
	 */
	private static final void buildBidirectedGraph_fixForks(int first, int last, int[] out) {
		final int TMP_ARRAY_UNIT = 50;  // For growing $tmpArray4$
		final int TRIM_WINDOW = 1;  // Arbitrary
		final int TRIM_DISTANCE = IO.quantum;  // Arbitrary
		boolean overhangSideFork, overhangSideForkPrime, trimmed, overhangSide1, overhangSide2;  // TRUE=the overhang is a prefix of the interval; FALSE=the overhang is a suffix.
		boolean trimmedFSide1, trimmedFSide2, trimmedOSide1, trimmedOSide2;
		int i, j, k;
		int lastNeighbor, start1, end1, start2, end2, lastRemoved, previousRemoved, id;
		int overhangStart1, overhangEnd1, overhangStart2, overhangEnd2; 
		int forkID, neighbor1ID, neighbor2ID, similarity, nDistinctRemoved, nMarked;
		final int maxDegree = BidirectedGraph.getMaxDegree();
		IntervalGraph.Node neighbor1, neighbor2, fork;
		BidirectedGraph.BidirectedEdge edge, otherEdge;
		int[] newArray;
		
		BidirectedGraph.cleanMarkedFlags(true);
		out[0]=0; out[1]=0;
		do {  // New nodes can become removable after each iteration
			out[2]=0; out[3]=0;
			for (i=first; i<=last; i++) {
				forkID=i-first;
				fork=IntervalGraph.nodesArray[nodesInKernel[i]];
				if (BidirectedGraph.removed[forkID] || fork.type==Constants.INTERVAL_PERIODIC) continue;
				lastNeighbor=BidirectedGraph.lastNeighbor[forkID];
				lastRemoved=-1;
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> CONSIDERING FORK "+IntervalGraph.nodesArray[nodesInKernel[i]]);				
				for (j=0; j<=lastNeighbor; j++) BidirectedGraph.neighbors[forkID][j].marked=false;
				nMarked=0;
				for (j=0; j<=lastNeighbor; j++) {
					edge=BidirectedGraph.neighbors[forkID][j];
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 0.1");
					if (edge.nAlignments!=1 || edge.marked) continue;
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 0.2");
					if (edge.node1==forkID) {
						neighbor1ID=edge.node2;
						neighbor1=BidirectedGraph.intervalGraphPointers[neighbor1ID];
						if (neighbor1.type==Constants.INTERVAL_PERIODIC) continue;
						start1=neighbor1.start; end1=neighbor1.end;
						overhangSideFork=edge.type==Constants.OVERLAP_SUFFIX_PREFIX||edge.type==Constants.OVERLAP_SUFFIX_SUFFIX;
						if (edge.type==Constants.OVERLAP_PREFIX_PREFIX || edge.type==Constants.OVERLAP_SUFFIX_PREFIX) {
							overhangSide1=false;
							overhangStart1=end1-edge.overhang2+1;
							overhangEnd1=end1;
						}
						else {
							overhangSide1=true;
							overhangStart1=start1;
							overhangEnd1=start1+edge.overhang2-1;
						}
					}
					else {
						neighbor1ID=edge.node1;
						neighbor1=BidirectedGraph.intervalGraphPointers[neighbor1ID];
						if (neighbor1.type==Constants.INTERVAL_PERIODIC) continue;
						start1=neighbor1.start; end1=neighbor1.end;
						overhangSideFork=edge.type==Constants.OVERLAP_PREFIX_SUFFIX||edge.type==Constants.OVERLAP_SUFFIX_SUFFIX;
						if (edge.type==Constants.OVERLAP_PREFIX_PREFIX || edge.type==Constants.OVERLAP_PREFIX_SUFFIX) {
							overhangSide1=false;
							overhangStart1=end1-edge.overhang1+1;
							overhangEnd1=end1;
						}
						else {
							overhangSide1=true;
							overhangStart1=start1;
							overhangEnd1=start1+edge.overhang1-1;
						}
					}
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> considering neighbor1: "+neighbor1);					
					if (!trimOverhang(neighbor1.read,overhangStart1,overhangEnd1,tmpArray1)) {
						edge.marked=true;
						nMarked++;
						continue;
					}
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 1  neighbor1 overhang before trimming=["+overhangStart1+".."+overhangEnd1+"]");
					trimmedFSide1=overhangSide1?tmpArray1[1]<overhangEnd1-TRIM_DISTANCE:tmpArray1[0]>overhangStart1+TRIM_DISTANCE;
					trimmedOSide1=overhangSide1?tmpArray1[0]>overhangStart1+TRIM_DISTANCE:tmpArray1[1]<overhangEnd1-TRIM_DISTANCE;
					overhangStart1=tmpArray1[0]; overhangEnd1=tmpArray1[1];
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 1  neighbor1 overhang after trimming=["+overhangStart1+".."+overhangEnd1+"]");
					for (k=j+1; k<=lastNeighbor; k++) {
						otherEdge=BidirectedGraph.neighbors[forkID][k];
						if (otherEdge.nAlignments!=1 || otherEdge.marked) continue;
						if (otherEdge.node1==forkID) {
							neighbor2ID=otherEdge.node2;
							neighbor2=BidirectedGraph.intervalGraphPointers[neighbor2ID];
							if (neighbor2.type==Constants.INTERVAL_PERIODIC || neighbor2.nodeID==neighbor1.nodeID) continue;
							start2=neighbor2.start; end2=neighbor2.end;
							overhangSideForkPrime=otherEdge.type==Constants.OVERLAP_SUFFIX_PREFIX||otherEdge.type==Constants.OVERLAP_SUFFIX_SUFFIX;
							if (overhangSideForkPrime!=overhangSideFork) continue;
							if (otherEdge.type==Constants.OVERLAP_PREFIX_PREFIX || otherEdge.type==Constants.OVERLAP_SUFFIX_PREFIX) {
								overhangSide2=false;
								overhangStart2=end2-otherEdge.overhang2+1;
								overhangEnd2=end2;
							}
							else {
								overhangSide2=true;
								overhangStart2=start2;
								overhangEnd2=start2+otherEdge.overhang2-1;
							}
						}
						else {
							neighbor2ID=otherEdge.node1;
							neighbor2=BidirectedGraph.intervalGraphPointers[neighbor2ID];
							if (neighbor2.type==Constants.INTERVAL_PERIODIC || neighbor2.nodeID==neighbor1.nodeID) continue;
							start2=neighbor2.start; end2=neighbor2.end;
							overhangSideForkPrime=otherEdge.type==Constants.OVERLAP_PREFIX_SUFFIX||otherEdge.type==Constants.OVERLAP_SUFFIX_SUFFIX;
							if (overhangSideForkPrime!=overhangSideFork) continue;
							if (otherEdge.type==Constants.OVERLAP_PREFIX_PREFIX || otherEdge.type==Constants.OVERLAP_PREFIX_SUFFIX) {
								overhangSide2=false;
								overhangStart2=end2-otherEdge.overhang1+1;
								overhangEnd2=end2;
							}
							else {
								overhangSide2=true;
								overhangStart2=start2;
								overhangEnd2=start2+otherEdge.overhang1-1;
							}
						}
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> considering neighbor2: "+neighbor2);
						if (!trimOverhang(neighbor2.read,overhangStart2,overhangEnd2,tmpArray1)) {
							otherEdge.marked=true;
							nMarked++;
							continue;
						}
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 2  neighbor2 overhang before trimming=["+overhangStart2+".."+overhangEnd2+"]");
						trimmedFSide2=overhangSide2?tmpArray1[1]<overhangEnd2-TRIM_DISTANCE:tmpArray1[0]>overhangStart2+TRIM_DISTANCE;
						trimmedOSide2=overhangSide2?tmpArray1[0]>overhangStart2+TRIM_DISTANCE:tmpArray1[1]<overhangEnd2-TRIM_DISTANCE;
						overhangStart2=tmpArray1[0]; overhangEnd2=tmpArray1[1];
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 2  neighbor2 overhang after trimming=["+overhangStart2+".."+overhangEnd2+"]");
						if (BidirectedGraph.neighborsOfForkOverlap(neighbor1ID,overhangSide1,neighbor2ID,overhangSide2)) {
							// This should be resolved later by transitive reduction
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 3");
							continue;
						}
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 4  neighbor1 overhang=["+overhangStart1+".."+overhangEnd1+"]  neighbor2 overhang=["+overhangStart2+".."+overhangEnd2+"]");
						id=-1;
						if (!trimmedFSide1) {
							if (!trimmedFSide2) {
								id=buildBidirectedGraph_fixForks_impl(0,neighbor1,neighbor1ID,overhangStart1,overhangEnd1,overhangSide1,trimmedOSide1,neighbor2,neighbor2ID,overhangStart2,overhangEnd2,overhangSide2,trimmedOSide2,fork,maxDegree);
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 4.1  id="+id);
							}
							else {
								id=buildBidirectedGraph_fixForks_impl(1,neighbor2,neighbor2ID,overhangStart2,overhangEnd2,overhangSide2,trimmedOSide2,neighbor1,neighbor1ID,overhangStart1,overhangEnd1,overhangSide1,trimmedOSide1,fork,maxDegree);
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 4.2 id="+id);
							}
						}
						else {
							if (!trimmedFSide2) {
								id=buildBidirectedGraph_fixForks_impl(1,neighbor1,neighbor1ID,overhangStart1,overhangEnd1,overhangSide1,trimmedOSide1,neighbor2,neighbor2ID,overhangStart2,overhangEnd2,overhangSide2,trimmedOSide2,fork,maxDegree);
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 4.3 id="+id);
							}
							else {
								id=buildBidirectedGraph_fixForks_impl(2,neighbor1,neighbor1ID,overhangStart1,overhangEnd1,overhangSide1,trimmedOSide1,neighbor2,neighbor2ID,overhangStart2,overhangEnd2,overhangSide2,trimmedOSide2,fork,maxDegree);
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 4.4 id="+id);
							}
						}						
						if (id>=0) {
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> 5  id="+id);
							lastRemoved++;
							if (lastRemoved==tmpArray4.length) {
								newArray = new int[tmpArray4.length+TMP_ARRAY_UNIT];
								System.arraycopy(tmpArray4,0,newArray,0,tmpArray4.length);
								tmpArray4=newArray;
							}
							tmpArray4[lastRemoved]=id;
						}
					}
				}
				// Removing the distinct nodes that have been marked for removal
				if (lastRemoved>=0) {
					if (lastRemoved>0) Arrays.sort(tmpArray4,0,lastRemoved+1);
					BidirectedGraph.removeNode(tmpArray4[0]);
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> ------ REMOVED NODE "+BidirectedGraph.intervalGraphPointers[tmpArray4[0]]);
					nDistinctRemoved=1; previousRemoved=tmpArray4[0];
					for (j=1; j<=lastRemoved; j++) {
						if (tmpArray4[j]==previousRemoved) continue;
						BidirectedGraph.removeNode(tmpArray4[j]);
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks> ------ REMOVED NODE "+BidirectedGraph.intervalGraphPointers[tmpArray4[j]]);
						previousRemoved=tmpArray4[j];
						nDistinctRemoved++;
					}
					out[2]+=nDistinctRemoved;
				}
				// Removing marked edges, and neighbors of the fork that get disconnected.
				// Remark: this removes edges from a neighbor of the fork, so subsequent 
				// calls to $neighborsAreSubset()$ might return false. We should instead
				// remove all marked edges at the very end of the process.
				if (BidirectedGraph.lastNeighbor[forkID]+1>tmpArray1.length) tmpArray1 = new int[BidirectedGraph.lastNeighbor[forkID]+1];
				out[2]+=BidirectedGraph.removeMarkedEdges(forkID,tmpArray1);
				out[3]+=nMarked;
			}
			out[0]+=out[2]; out[1]+=out[3];
		}
		while (out[2]+out[3]>0);
	}
	
	
	/**
	 * Similar to the $trim()$ function of intervals used in factorization, but tries to 
	 * trim both sides of an overlap overhang $[start..end]$, and uses a more stringent 
	 * quality threshold.
	 *
	 * @param tmpArray temporary space;
	 * @return TRUE if the overhang can be trimmed: the new interval is stored in 
	 * $tmpArray$; FALSE if the overhang is completely low-quality, or if it becomes too 
	 * short after trimming.
	 */
	private static final boolean trimOverhang(int readID, int start, int end, int[] tmpArray) {
		final byte OLD_THRESHOLD = Reads.MIN_RANDOM_QUALITY_SCORE;
		final byte NEW_THRESHOLD = (byte)(Reads.MAX_HIGH_QUALITY_SCORE+1);  // More stringent threshold than the one used in factorization.
		final int WINDOW_LARGE = 3;  // Arbitrary
		final int WINDOW_SMALL = 1;
		final int MIN_LENGTH = Alignments.minAlignmentLength>>2;  // Arbitrary
		final int OLD_READ_A = ReadA.id;
		final int length = end-start+1;
		double quality;
		
if (IO.SHOW_INTERACTIVE) System.err.println("trimOverhang> 1  Alignments.minAlignmentLength="+Alignments.minAlignmentLength);
		Reads.readEnds2qualityEnds(readID,start,end,true,tmpArray);
		quality=Histograms.getAverage(Reads.getQualityArray(readID),tmpArray[0],tmpArray[1]);
		if (quality>=NEW_THRESHOLD) return false;
		
		ReadA.id=readID;
		Reads.MIN_RANDOM_QUALITY_SCORE=NEW_THRESHOLD;
		tmpArray[0]=start; tmpArray[1]=end;
		Reads.trim(tmpArray,WINDOW_LARGE,true,true);
		if (tmpArray[0]==-1 || tmpArray[1]==-1) return false;
if (IO.SHOW_INTERACTIVE) System.err.println("trimOverhang> 3  newLength="+(tmpArray[1]-tmpArray[0]+1)+" oldLength="+length);		
		if (tmpArray[1]-tmpArray[0]+1<MIN_LENGTH) return false;
		Reads.trim(tmpArray,WINDOW_SMALL,true,true);
		if (tmpArray[0]==-1 || tmpArray[1]==-1) return false;
if (IO.SHOW_INTERACTIVE) System.err.println("trimOverhang> 4  newLength="+(tmpArray[1]-tmpArray[0]+1)+" oldLength="+length);
		if (tmpArray[1]-tmpArray[0]+1<MIN_LENGTH) return false;
if (IO.SHOW_INTERACTIVE) System.err.println("trimOverhang> 5");
		Reads.MIN_RANDOM_QUALITY_SCORE=OLD_THRESHOLD;
		ReadA.id=OLD_READ_A;
		return true;
	}
	
	
	/**
	 * Handles the following cases: ($mode=0$) neither neighbor has been trimmed on the 
	 * side of the fork; ($mode=1$) $neighbor1$ has been trimmed on the side of the fork, 
	 * $neighbor2$ might have been trimmed on the side of the fork or not; ($mode=2$) both
	 * neighbors have been trimmed on the side of the fork.
	 *
	 * Remark: the procedure assumes that no node in a kernel is contained in any other
	 * node; thus, the neighbor of a node X that is not also a neighbor of another node Y,
	 * cannot be contained in Y.
	 *
	 * @param neighbor*ID IDs in the bidirected graph;
	 * @param overhangStart*,overhangEnd* after trimming;
	 * @param maxDegree the max degree of a node in the current bidirected graph;
	 * @return $neighborXID$ if $neighborXID$ has to be removed; -1 if no neighbor should 
	 * be removed.
	 */
	private static final int buildBidirectedGraph_fixForks_impl(int mode, IntervalGraph.Node neighbor1, int neighbor1ID, int overhangStart1, int overhangEnd1, boolean overhangSide1, boolean trimmedOSide1, IntervalGraph.Node neighbor2, int neighbor2ID, int overhangStart2, int overhangEnd2, boolean overhangSide2, boolean trimmedOSide2, IntervalGraph.Node fork, int maxDegree) {
		final int similarity = similarOverhangs(mode,neighbor1,overhangStart1,overhangEnd1,overhangSide1,neighbor2,overhangStart2,overhangEnd2,overhangSide2);
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks_impl> mode="+mode+" similarity="+similarity+" overhang1=["+overhangStart1+".."+overhangEnd1+"] overhang2=["+overhangStart2+".."+overhangEnd2+"]");
		if (similarity<0) return -1;

		if (mode==0 || mode==2) {
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks_impl> similarity="+similarity);
			if (similarity==0) {
				if ( overhangEnd1-overhangStart1<overhangEnd2-overhangStart2 && 
				     BidirectedGraph.neighborsAreSubset(neighbor1ID,overhangSide1,neighbor2ID,overhangSide2) &&
					 edgesOnSide(neighbor1ID,!overhangSide1,fork,maxDegree)==1
				   ) return neighbor1ID;
				else if ( overhangEnd2-overhangStart2<overhangEnd1-overhangStart1 && 
				          BidirectedGraph.neighborsAreSubset(neighbor2ID,overhangSide2,neighbor1ID,overhangSide1) &&
						  edgesOnSide(neighbor2ID,!overhangSide2,fork,maxDegree)==1
				        ) return neighbor2ID;
			}
			else if ( similarity==1 && 
			          BidirectedGraph.neighborsAreSubset(neighbor1ID,overhangSide1,neighbor2ID,overhangSide2) &&
					  edgesOnSide(neighbor1ID,!overhangSide1,fork,maxDegree)==1
			        ) {
if (IO.SHOW_INTERACTIVE) System.err.println("fixForks_impl> 1  returning "+neighbor1ID);
						return neighbor1ID;
					}
			else if ( similarity==2 &&
			          BidirectedGraph.neighborsAreSubset(neighbor2ID,overhangSide2,neighbor1ID,overhangSide1) &&
					  edgesOnSide(neighbor2ID,!overhangSide2,fork,maxDegree)==1
		            ) return neighbor2ID;
		}
		else if (mode==1) {
			if ( similarity==1 &&
				 BidirectedGraph.neighborsAreSubset(neighbor1ID,overhangSide1,neighbor2ID,overhangSide2) &&
			     edgesOnSide(neighbor1ID,!overhangSide1,fork,maxDegree)==1
			   ) return neighbor1ID;
		}
if (IO.SHOW_INTERACTIVE) {
System.err.println("fixForks_impl> 2  returning -1");
System.err.println("fixForks_impl> neighborsAreSubset? "+BidirectedGraph.neighborsAreSubset(neighbor1ID,overhangSide1,neighbor2ID,overhangSide2));
System.err.println("fixForks_impl> edgesOnSide? "+edgesOnSide(neighbor1ID,!overhangSide1,fork,maxDegree));
System.err.println("fixForks_impl> NEIGHBORS OF neighbor1 = "+BidirectedGraph.intervalGraphPointers[neighbor1ID]);
for (int x=0; x<=BidirectedGraph.lastNeighbor[neighbor1ID]; x++) {
	System.err.println(BidirectedGraph.neighbors[neighbor1ID][x]+"\n node: "+BidirectedGraph.intervalGraphPointers[BidirectedGraph.neighbors[neighbor1ID][x].getTo(neighbor1ID)]);
}
System.err.println("fixForks_impl> NEIGHBORS OF neighbor2 = "+BidirectedGraph.intervalGraphPointers[neighbor2ID]);
for (int x=0; x<=BidirectedGraph.lastNeighbor[neighbor2ID]; x++) {
	System.err.println(BidirectedGraph.neighbors[neighbor2ID][x]+"\n node: "+BidirectedGraph.intervalGraphPointers[BidirectedGraph.neighbors[neighbor2ID][x].getTo(neighbor2ID)]);
}
}
		return -1;
	}
	
	
	/**
	 * Wraps the corresponding procedure of $BidirectedGraph$, but returns the number of 
	 * distinct nodes, reached by edges incident to $node$ on its $side$, that are either 
	 * equal to $fork$, or have no identity edge with $fork$ in the interval graph.
	 * This is useful, since using just the $BidirectedGraph$ procedure would not allow to
	 * remove any neighbor of e.g. two forks that are identical and have the same overlap
	 * neighbors.
	 *
	 * Remark: the procedure uses global array $tmpNodes$.
	 *
	 * @param maxDegree the max degree of a node in the current bidirected graph.
	 */
	private static final int edgesOnSide(int node, boolean side, IntervalGraph.Node fork, int maxDegree) {
		boolean isIdentical;
		int i, j, k;
		int nEdges, previousOrder, lastNode, nodeID, out;
		IntervalGraph.Node tmpNode, previousNode;
		IntervalGraph.Edge edge;
		
		if (tmpNodes.length<maxDegree) tmpNodes = new IntervalGraph.Node[maxDegree];
		nEdges=BidirectedGraph.edgesOnSide(node,side,tmpNodes);
		if (nEdges==1) return 1;
		
		// Compacting $tmpNodes$.
		previousOrder=IntervalGraph.Node.order;
		IntervalGraph.Node.order=IntervalGraph.Node.NODE_ID;
		Arrays.sort(tmpNodes,0,nEdges);
		IntervalGraph.Node.order=previousOrder;
		previousNode=tmpNodes[0];
		j=0;
		for (i=1; i<nEdges; i++) {
			if (tmpNodes[i]==previousNode) continue;
			j++;
			tmpNode=tmpNodes[j];
			tmpNodes[j]=tmpNodes[i];
			tmpNodes[i]=tmpNode;
		}
		if (j==0) return 1;
		lastNode=j;
		
		// Counting the distinct nodes without identity edges to $fork$.
		out=0;
		for (i=0; i<=lastNode; i++) {
			tmpNode=tmpNodes[i];
			if (tmpNode==fork) {
				out++;
				continue;
			}
			nodeID=tmpNode.nodeID;
			isIdentical=false;
			for (j=0; j<IntervalGraph.nNeighbors[nodeID]; j++) {
				edge=IntervalGraph.neighbors[nodeID][j];
				if (edge.containment==Constants.CONTAINMENT_IDENTICAL && IntervalGraph.nodesArray[edge.getTo(nodeID)]==fork) {
					isIdentical=true;
					break;
				}
			}
			if (!isIdentical) out++;
		}
		return out;
	}
	
	
	/**
	 * If $mode=0$, checks whether the two overhangs are identical, or if one is the 
	 * prefix/suffix of the other on the side of the fork. If $mode=1$, checks whether
	 * the overhang of $neighbor1$ is strictly contained in the overhang of $neighbor2$.
	 * If $mode=2$, checks whether the two overhangs are identical, or if one is contained
	 * in the other.
	 *
	 * Remark: the procedure uses $tmpArray1,tmpArray2,tmpArray3$.
	 *
	 * @param overhangStart*,overhangEnd* absolute positions;
	 * @param overhangSide* side of the neighbor that is involved in the overhang: 
	 * TRUE=prefix; FALSE=suffix;
	 * @return 0=the overhangs are approximately identical; 1=the overhang of $neighbor1$
	 * is approximately contained in the overhang of $neighbor2$; 2=vice versa; -1=none of
	 * the above; -2=contrasting signals.
	 */
	private static final int similarOverhangs(int mode, IntervalGraph.Node neighbor1, int overhangStart1, int overhangEnd1, boolean overhangSide1, IntervalGraph.Node neighbor2, int overhangStart2, int overhangEnd2, boolean overhangSide2) {
		final int neighbor1Start = neighbor1.start;
		final int neighbor2Start = neighbor2.start;
		final int neighbor1ID = neighbor1.nodeID;
		final int neighbor2ID = neighbor2.nodeID;
		final int nNeighbors = IntervalGraph.nNeighbors[neighbor1ID];
		int i, j;
		int out, last1, last2, alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		IntervalGraph.Edge edge;		
		
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 0 nNeighbors="+nNeighbors+" neighbor1ID="+neighbor1ID+" neighbor2ID="+neighbor2ID);
		out=-1;		
		for (i=0; i<nNeighbors; i++) {
			edge=IntervalGraph.neighbors[neighbor1ID][i];
			if (edge.orientation==2 || edge.orientation==-1) {
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> EXIT because orientation="+edge.orientation);
				continue; // We cannot decompress the orientation of each contained alignment
			}
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 0.1  read1="+IntervalGraph.nodesArray[edge.nodeID1].read+" read2="+IntervalGraph.nodesArray[edge.nodeID2].read);
			if ((edge.nodeID1==neighbor1ID && edge.nodeID2==neighbor2ID) || (edge.nodeID1==neighbor2ID && edge.nodeID2==neighbor1ID)) {
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> considering edge: "+edge);
				last1=edge.projectOnto(neighbor1ID,tmpArray1,true);
				if (last1==-1) continue;
				last2=edge.projectOnto(neighbor2ID,tmpArray2,true);
				if (last2==-1 || last2!=last1) continue;
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 1");
				for (j=0; j<last1; j+=2) {
					alignmentStart1=neighbor1Start+tmpArray1[j];
					alignmentEnd1=neighbor1Start+tmpArray1[j+1];
					alignmentStart2=neighbor2Start+tmpArray2[j];
					alignmentEnd2=neighbor2Start+tmpArray2[j+1];
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 2 considering alignment projections: ["+alignmentStart1+".."+alignmentEnd1+"] x ["+alignmentStart2+".."+alignmentEnd2+"]");
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 2 and overhangs: ["+overhangStart1+".."+overhangEnd1+"], ["+overhangStart2+".."+overhangEnd2+"]");
					if (mode==0) {
						if (Intervals.areApproximatelyIdentical(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,edge.orientation==0)) {
							if (out==-1) out=0;
							else if (out!=0) return -2;
						}
						else if (Intervals.isApproximatePrefixOrSuffix(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,edge.orientation==0,!overhangSide2)) {
							if (out==-1) out=1;
							else if (out!=1) return -2;
						}
						else if (Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,overhangStart1,overhangEnd1,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,edge.orientation==0,!overhangSide1)) {
							if (out==-1) out=2;
							else if (out!=2) return -2;
						}
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 3 out="+out);
					}
					else if (mode==1) {
						if (Intervals.isApproximatelyContained(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,edge.orientation==0,tmpArray3)) {
							out=1;
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 4 out="+out);
							break;
						}
					}
					else if (mode==2) {
						if (Intervals.areApproximatelyIdentical(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,edge.orientation==0)) {
							if (out==-1) out=0;
							else if (out!=0) return -2;
						}
						else if (Intervals.isApproximatelyContained(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,edge.orientation==0,tmpArray3)) {
							if (out==-1) out=1;
							else if (out!=1) return -2;
						}
						else if (Intervals.isApproximatelyContained(overhangStart2,overhangEnd2,overhangStart1,overhangEnd1,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,edge.orientation==0,tmpArray3)) {
							if (out==-1) out=2;
							else if (out!=2) return -2;
						}
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 5 out="+out);
					}
				}
			}
		}
if (IO.SHOW_INTERACTIVE) System.err.println("similarOverhangs> 6 out="+out);
		return out;
	}
	
	
	/**
	 * Merges nodes of the bidirected graph that correspond to overlapping intervals in
	 * the same read (i.e. approximately to the same substring of the same read).
	 * Merging straddling intervals in the same read is not a good idea in general: we do 
	 * so just because both intervals are believed to belong to the same repeat at this 
	 * step of the pipeline.
	 *
	 * The procedure also removes nodes of the bidirected graph that correspond to
	 * intervals that are contained in another interval of the same read. Doing so is 
	 * correct just because both intervals are believed to belong to the same repeat at 
	 * this step of the pipeline.
	 *
	 * Remark: this simplification applies also to periodic kernels.
	 *
	 * Remark: the procedure erases removed nodes from the bidirected graph completely.
	 * $node2newNode$ is set to valid values for all nodes removed by this procedure (a
	 * removed node is mapped to the new node to which a containing node in its read is 
	 * mapped). Nodes that were removed before this procedure was called have 
	 * $node2newNode[0]=-1$.
	 *
	 * Remark: the procedure does not assume the neighbors of a node of the bidirected 
	 * graph to be sorted.
	 *
	 * @param newNodes output matrix: read, startPosition, endPosition, type, of every 
	 * node in the new overlap graph; the length of the matrix must be at least 
	 * $last-first+1$; type is -2 if interval graph nodes of different type are merged in 
	 * the same bidirected graph node;
	 * @param newNodes_pointers for each row in $newNodes$, a pointer to a longest 
	 * interval graph node that was mapped to this bidirected graph node;
	 * @param node2newNode for each node $i$ in $[first..last]$: ID of the new node to 
	 * which it is mapped; relative start; relative end; this holds for nodes removed by
	 * this procedure as well;
	 * @param tmpNodes temporary space, of size at least $last-first+1$;
	 * @return the number of nodes in the new bidirected graph, or -1 if such number has 
	 * not changed; note that, even if no merge is performed, the resulting graph might 
	 * still be smaller than the original, since removed nodes are deleted; even if the 
	 * bidirected graph is not changed, $newNodes$ and $node2newNode$ are set correctly.
	 */
	private static final int buildBidirectedGraph_mergeStraddling(int first, int last, int[][] newNodes, IntervalGraph.Node[] newNodes_pointers, int[][] node2newNode, IntervalGraph.Node[] tmpNodes) {
		final double JACCARD_THRESHOLD = 0.9;  // Arbitrary
		boolean previouslyRemoved;
		int i, j;
		int nNodes, previousOrder, lastNewNode, from, top, newNode, minContainer, minLength;
		IntervalGraph.Node nodeI, nodeJ;

if (IO.SHOW_INTERACTIVE) {
	System.err.println("buildBidirectedGraph_mergeStraddling> bidir. nodes at the very beginning:");
	for (int x=0; x<BidirectedGraph.nNodes; x++) System.err.println(BidirectedGraph.removed[x]+" :: "+BidirectedGraph.intervalGraphPointers[x]);
}

		// Checking whether nodes were removed before this procedure was called.
		previouslyRemoved=false;
		for (i=first; i<=last; i++) {
			if (BidirectedGraph.removed[i-first]) {
				previouslyRemoved=true;
				break;
			}
		}

		// Removing contained nodes
		nNodes=last-first+1;
		Math.set(newNodes,-1);
		for (i=0; i<newNodes_pointers.length; i++) newNodes_pointers[i]=null;
		Math.set(node2newNode,-1);
		for (i=first; i<=last; i++) {
			if (BidirectedGraph.removed[i-first]) continue;
			nodeI=IntervalGraph.nodesArray[nodesInKernel[i]];
			for (j=first; j<=last; j++) {
				if (j==i) continue;
				nodeJ=IntervalGraph.nodesArray[nodesInKernel[j]];
				if (nodeJ.read!=nodeI.read) continue;
				if (Intervals.isContained(nodeI.start,nodeI.end,nodeJ.start,nodeJ.end)) {
					BidirectedGraph.removeNode(i-first);

if (IO.SHOW_INTERACTIVE) {
	System.err.println("buildBidirectedGraph_mergeStraddling> REMOVED node "+BidirectedGraph.intervalGraphPointers[i-first]);					
	System.err.println("buildBidirectedGraph_mergeStraddling> because contained in "+nodeJ);
}
					break;
				}
			}
		}
		for (i=first; i<=last; i++) {
			if (!BidirectedGraph.removed[i-first]) continue;
			nodeI=IntervalGraph.nodesArray[nodesInKernel[i]];
			minContainer=-1; minLength=Math.POSITIVE_INFINITY;
			for (j=first; j<=last; j++) {
				if (j==i || BidirectedGraph.removed[j-first]) continue;
				nodeJ=IntervalGraph.nodesArray[nodesInKernel[j]];
				if (nodeJ.read!=nodeI.read || !Intervals.isContained(nodeI.start,nodeI.end,nodeJ.start,nodeJ.end)) continue;
				if (nodeJ.length()<minLength) {
					minContainer=j-first;
					minLength=nodeJ.length();
				}
			}
			if (IO.CONSISTENCY_CHECKS && minContainer==-1 && !previouslyRemoved) {
				System.err.println("buildBidirectedGraph_mergeStraddling> ERROR: the following removed node has no non-removed container?!");
				System.err.println(nodeI);
				System.exit(1);
			}
			node2newNode[i-first][0]=minContainer;
		}

		// Merging non-removed straddling nodes
		j=-1;
		for (i=first; i<=last; i++) {
			if (BidirectedGraph.removed[i-first]) continue;
			j++;
			tmpNodes[j]=IntervalGraph.nodesArray[nodesInKernel[i]];
			tmpNodes[j].visited=i-first;  // Storing in $visited$ the original relative position of a node in $nodesInKernel$, i.e. its original node ID in the bidirected graph.
		}
		nNodes=j+1;
		previousOrder=IntervalGraph.Node.order;
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		if (nNodes>1) Arrays.sort(tmpNodes,0,nNodes);
		IntervalGraph.Node.order=previousOrder;
		lastNewNode=-1;
		for (i=0; i<nNodes; i++) {		
			if (node2newNode[tmpNodes[i].visited][0]!=-1) continue;
			node2newNode[tmpNodes[i].visited][0]=++lastNewNode;
			newNodes[lastNewNode][0]=tmpNodes[i].read;
			IntervalGraph.stack[0]=i; top=0;
			while (top>=0) {
				from=IntervalGraph.stack[top--];
				for (j=from+1; j<nNodes; j++) {
					if (tmpNodes[j].read!=tmpNodes[from].read || tmpNodes[j].start>=tmpNodes[from].end) break;
					if (node2newNode[tmpNodes[j].visited][0]!=-1) continue;
					if (Intervals.jaccardSimilarity(tmpNodes[from].start,tmpNodes[from].end,tmpNodes[j].start,tmpNodes[j].end)>=JACCARD_THRESHOLD) {
						node2newNode[tmpNodes[j].visited][0]=lastNewNode;
						IntervalGraph.stack[++top]=j;
					}
				}
				for (j=from-1; j>=0; j--) {
					if (tmpNodes[j].read!=tmpNodes[from].read) break;
					if (node2newNode[tmpNodes[j].visited][0]!=-1 || tmpNodes[j].end<=tmpNodes[from].start) continue;
					if (Intervals.jaccardSimilarity(tmpNodes[from].start,tmpNodes[from].end,tmpNodes[j].start,tmpNodes[j].end)>=JACCARD_THRESHOLD) {
						node2newNode[tmpNodes[j].visited][0]=lastNewNode;
						IntervalGraph.stack[++top]=j;
					}
				}
			}
		}
		for (i=0; i<nNodes; i++) {
			from=node2newNode[tmpNodes[i].visited][0];
			if (newNodes[from][1]==-1) newNodes[from][1]=tmpNodes[i].start;
			else if (tmpNodes[i].start<newNodes[from][1]) newNodes[from][1]=tmpNodes[i].start;
			if (newNodes[from][2]==-1) newNodes[from][2]=tmpNodes[i].end;
			else if (tmpNodes[i].end>newNodes[from][2]) newNodes[from][2]=tmpNodes[i].end;
			if (newNodes[from][3]==-1) newNodes[from][3]=tmpNodes[i].type;
			else if (newNodes[from][3]>=0 && newNodes[from][3]!=tmpNodes[i].type) newNodes[from][3]=-2;
			if (newNodes_pointers[from]==null || tmpNodes[i].length()>newNodes_pointers[from].length()) newNodes_pointers[from]=tmpNodes[i];
		}
		for (i=0; i<nNodes; i++) {
			from=node2newNode[tmpNodes[i].visited][0];
			node2newNode[tmpNodes[i].visited][1]=tmpNodes[i].start-newNodes[from][1];
			node2newNode[tmpNodes[i].visited][2]=tmpNodes[i].end-newNodes[from][1];
		}

		// Updating $node2newNode$ for nodes removed by this procedure.
		for (i=first; i<=last; i++) {
			if (!BidirectedGraph.removed[i-first]) continue;
			newNode=node2newNode[i-first][0];
			if (newNode==-1) continue;
			newNode=node2newNode[newNode][0];
			node2newNode[i-first][0]=newNode;
			nodeI=IntervalGraph.nodesArray[nodesInKernel[i]];
			if (nodeI.start<newNodes[newNode][1]) node2newNode[i-first][1]=0;
			else node2newNode[i-first][1]=nodeI.start-newNodes[newNode][1];
			if (nodeI.end>newNodes[newNode][2]) node2newNode[i-first][2]=newNodes[newNode][2]-newNodes[newNode][1];
			else node2newNode[i-first][2]=nodeI.end-newNodes[newNode][1];
		}
		
if (IO.SHOW_INTERACTIVE) {
	System.err.println("buildBidirectedGraph_mergeStraddling> bidir. nodes before contract:");
	for (int x=0; x<BidirectedGraph.nNodes; x++) System.err.println(BidirectedGraph.removed[x]+" :: "+BidirectedGraph.intervalGraphPointers[x]);
}
		
		if (lastNewNode+1<nNodes || BidirectedGraph.nRemoved>0) BidirectedGraph.contract(newNodes,lastNewNode,newNodes_pointers,node2newNode,true);
	
if (IO.SHOW_INTERACTIVE) {	
	System.err.println("buildBidirectedGraph_mergeStraddling> bidir. nodes AFTER contract:");
	for (int x=0; x<BidirectedGraph.nNodes; x++) System.err.println(BidirectedGraph.removed[x]+" :: "+BidirectedGraph.intervalGraphPointers[x]);		
}	
		
		return lastNewNode<last-first?lastNewNode+1:-1;
	}
	
	
	/**
	 * Merges nodes of the bidirected graph generated by $buildBidirectedGraph_
	 * mergeStraddling()$, iff: (1) they contain large sub-intervals connected by identity 
	 * edges in the interval graph (both regular edges and supplement edges); 
	 * (2) they have the same neighbors in the bidirected graph, or one of them has no
	 * neighbor.
	 *
	 * Remark: the resulting bidirected graph is not necessarily connected, since: (1) the 
	 * set of kernel nodes was not necessarily connected in the interval graph when 
	 * considering only overlap edges; (2) the procedure is not using, for merging, every 
	 * identity edge used to build the kernel. If it were, the bidirected graph would be 
	 * connected.
	 *
	 * Remark: this simplification applies also to periodic kernels.
	 *
	 * Remark: the procedure assumes the neighbors of a node of the bidirected graph to be
	 * sorted.
	 *
	 * @param oldNodes output of $buildBidirectedGraph_mergeStraddling()$;
	 * @param node2oldNode mapping generated by $buildBidirectedGraph_mergeStraddling()$;
	 * @param newNodes with at least 5 columns; row $i$ contains the read, start, end,
	 * of a longest node produced by $buildBidirectedGraph_mergeStraddling()$, that is 
	 * mapped to node $i$ by this procedure (i.e. a longest representative of new
	 * node $i$); the fourth column contains the merge of the types of all nodes that are
	 * mapped to the new node (-2 if such merge is not valid); the fifth column is used as
	 * temporary space;
	 * @param oldNode2newNode with at least two columns; mapping generated by this 
	 * procedure (column zero); column one is zero (respectively, one) if the old node is 
	 * in the same (opposite) orientation as the new node; it equals two if the old node 
	 * maps in both orientations;
	 * @param outNeighbors,orientations temporary space, of size at least $nOldNodes$
	 * times $nOldNodes$;
	 * @param lastOutNeighbor,nInNeighbors,component temporary space, of size at least 
	 * $nOldNodes$;
	 * @return the number of nodes in the new bidirected graph, or -1 if the graph was not
	 * changed; even if the graph was not changed, $oldNode2newNode$ and $newNodes$ are 
	 * set correctly.
	 */
	private static final int buildBidirectedGraph_mergeIdentical(int first, int last, int[][] oldNodes, int nOldNodes, IntervalGraph.Node[] oldNodes_pointers, int[][] node2oldNode, int[][] newNodes, IntervalGraph.Node[] newNodes_pointers, int[][] oldNode2newNode, int[][] outNeighbors, int[] lastOutNeighbor, int[][] orientations, int[] nInNeighbors, int[] component) {
		final double JACCARD_THRESHOLD = 0.9;  // Arbitrary
		final int GROWTH_RATE = 10;  // Arbitrary
		boolean atLeastOneMerge;
		int i, j, k, n, p;
		int top, from, to, oldFrom, oldTo, type;
		int lastNewNode, newNode, kernel, orientation, otherOrientation;
		IntervalGraph.Node tmpNode;
		IntervalGraph.Edge edge;
		BidirectedGraph.BidirectedEdge bidirectedEdge = BidirectedGraph.getEdge();
		int[] tmpArray;
		
		// Merging old nodes
		Math.set(oldNode2newNode,-1);
		for (i=0; i<newNodes_pointers.length; i++) newNodes_pointers[i]=null;
		lastNewNode=-1;
		kernel=IntervalGraph.nodesArray[nodesInKernel[first]].kernels[0];
		Math.set(lastOutNeighbor,nOldNodes-1,-1);
		Math.set(orientations,-1);
		atLeastOneMerge=false;
		for (i=first; i<=last; i++) {
			from=nodesInKernel[i];
			n=IntervalGraph.nNeighbors[from];
			for (j=0; j<n; j++) {
				edge=IntervalGraph.neighbors[from][j];
				edge.printed=-1;  // Using the $printed$ field of each interval graph edge to mark whether an edge was already used by this procedure.
			}
		}
		for (i=first; i<=last; i++) {
			from=nodesInKernel[i];
			oldFrom=node2oldNode[i-first][0];
			if (oldFrom==-1 || IntervalGraph.nodesArray[from].length()<(oldNodes[oldFrom][2]-oldNodes[oldFrom][1]+1)*JACCARD_THRESHOLD) continue;
			n=IntervalGraph.nNeighbors[from];
			for (j=0; j<n; j++) {
				edge=IntervalGraph.neighbors[from][j];
				if (edge.printed==1 || !isIdentity(edge)) continue;
				to=IntervalGraph.neighbors[from][j].getTo(from);
				tmpNode=IntervalGraph.nodesArray[to];
				if (tmpNode.lastKernel==-1 || tmpNode.kernels[0]!=kernel) continue;
				p=Arrays.binarySearch(nodesInKernel,first,last+1,to);
				if (p<0) continue;
				oldTo=node2oldNode[p-first][0];
				if (oldTo==-1 || IntervalGraph.nodesArray[to].length()<(oldNodes[oldTo][2]-oldNodes[oldTo][1]+1)*JACCARD_THRESHOLD) continue;
				orientation=IntervalGraph.neighbors[from][j].orientation;
				if (!BidirectedGraph.sameNeighbors(oldFrom,oldTo,orientation,bidirectedEdge) && !(BidirectedGraph.lastNeighbor[oldFrom]==-1 || BidirectedGraph.lastNeighbor[oldTo]==-1)) continue;
				k=Math.linearSearch_unsorted(outNeighbors[oldFrom],0,lastOutNeighbor[oldFrom]+1,oldTo,false);
				if (k>=0) {
					if (orientations[oldFrom][k]!=2) {
						if (orientation==2 || orientation!=orientations[oldFrom][k]) orientations[oldFrom][k]=2;
					}
				}
				else {
					lastOutNeighbor[oldFrom]++;
					if (lastOutNeighbor[oldFrom]>=outNeighbors[oldFrom].length) {
						tmpArray = new int[outNeighbors[oldFrom].length+GROWTH_RATE];
						System.arraycopy(outNeighbors[oldFrom],0,tmpArray,0,outNeighbors[oldFrom].length);
						outNeighbors[oldFrom]=tmpArray;
						tmpArray = new int[orientations[oldFrom].length+GROWTH_RATE];
						System.arraycopy(orientations[oldFrom],0,tmpArray,0,orientations[oldFrom].length);
						orientations[oldFrom]=tmpArray;
					}
					outNeighbors[oldFrom][lastOutNeighbor[oldFrom]]=oldTo;
					orientations[oldFrom][lastOutNeighbor[oldFrom]]=orientation;
				}
				k=Math.linearSearch_unsorted(outNeighbors[oldTo],0,lastOutNeighbor[oldTo]+1,oldFrom,false);
				if (k>=0) {
					if (orientations[oldTo][k]!=2) {
						if (orientation==2 || orientation!=orientations[oldTo][k]) orientations[oldTo][k]=2;
					}
				}
				else {
					lastOutNeighbor[oldTo]++;
					if (lastOutNeighbor[oldTo]>=outNeighbors[oldTo].length) {
						tmpArray = new int[outNeighbors[oldTo].length+GROWTH_RATE];
						System.arraycopy(outNeighbors[oldTo],0,tmpArray,0,outNeighbors[oldTo].length);
						outNeighbors[oldTo]=tmpArray;
						tmpArray = new int[orientations[oldTo].length+GROWTH_RATE];
						System.arraycopy(orientations[oldTo],0,tmpArray,0,orientations[oldTo].length);
						orientations[oldTo]=tmpArray;
					}
					outNeighbors[oldTo][lastOutNeighbor[oldTo]]=oldFrom;
					orientations[oldTo][lastOutNeighbor[oldTo]]=orientation;
				}
				atLeastOneMerge=true;
				IntervalGraph.neighbors[from][j].printed=1;
			}
		}
		if (!atLeastOneMerge) {
			for (i=0; i<nOldNodes; i++) {
				oldNode2newNode[i][0]=i;
				oldNode2newNode[i][1]=0;
			}
			for (i=0; i<nOldNodes; i++) {
				newNodes[i][0]=oldNodes[i][0];
				newNodes[i][1]=oldNodes[i][1];
				newNodes[i][2]=oldNodes[i][2];
				newNodes[i][3]=oldNodes[i][3];
				newNodes_pointers[i]=oldNodes_pointers[i];
			}
			return -1;
		}
		for (i=0; i<nOldNodes; i++) lastOutNeighbor[i]++;
		lastNewNode=DAG.getConnectedComponents(nOldNodes,null,nInNeighbors,outNeighbors,lastOutNeighbor,component,IntervalGraph.stack)-1;
		for (i=0; i<nOldNodes; i++) lastOutNeighbor[i]--;
		Math.set(newNodes,-1);
		for (i=0; i<nOldNodes; i++) {
			newNode=component[i];
			oldNode2newNode[i][0]=newNode;
			if (newNodes[newNode][0]==-1 || oldNodes[i][2]-oldNodes[i][1]>newNodes[newNode][2]-newNodes[newNode][1]) {
				newNodes[newNode][0]=oldNodes[i][0];
				newNodes[newNode][1]=oldNodes[i][1];
				newNodes[newNode][2]=oldNodes[i][2];
				newNodes[newNode][4]=i;
				newNodes_pointers[newNode]=oldNodes_pointers[i];
			}
		}
		
		// Orienting old nodes with respect to the representative new node
		for (i=0; i<=lastNewNode; i++) {
			from=newNodes[i][4];
			oldNode2newNode[from][1]=0;
			IntervalGraph.stack[0]=from;
			top=0;
			while (top>=0) {
				from=IntervalGraph.stack[top--];
				orientation=oldNode2newNode[from][1];
				for (j=0; j<=lastOutNeighbor[from]; j++) {
					to=outNeighbors[from][j];
					if (orientations[from][j]==0) otherOrientation=orientation;
					else if (orientations[from][j]==1) {
						if (orientation==0) otherOrientation=1;
						else if (orientation==1) otherOrientation=0;
						else otherOrientation=2;
					}
					else otherOrientation=2;
					if (oldNode2newNode[to][1]==-1) {
						oldNode2newNode[to][1]=otherOrientation;
						top++;
						if (top==IntervalGraph.stack.length) {
							int[] newStack = new int[IntervalGraph.stack.length<<1];
							System.arraycopy(IntervalGraph.stack,0,newStack,0,top);
							IntervalGraph.stack=newStack;
						}
						IntervalGraph.stack[top]=to;
						continue;
					}
					else if (oldNode2newNode[to][1]==otherOrientation || oldNode2newNode[to][1]==2) continue;
					else {
						oldNode2newNode[to][1]=2;
						top++;
						if (top==IntervalGraph.stack.length) {
							int[] newStack = new int[IntervalGraph.stack.length<<1];
							System.arraycopy(IntervalGraph.stack,0,newStack,0,top);
							IntervalGraph.stack=newStack;
						}
						IntervalGraph.stack[top]=to;
					}
				}
			}
		}
		
		// Deciding the type of $newNodes$.
		for (i=0; i<nOldNodes; i++) {
			newNode=oldNode2newNode[i][0];
			type=Constants.typeInOrientation(oldNodes[i][3],oldNode2newNode[i][1]);
			if (newNodes[newNode][3]==-1) newNodes[newNode][3]=type;
			else newNodes[newNode][3]=Constants.mergeTypes(newNodes[newNode][3],type);
		}
		
		// Contracting the bidirected graph
		BidirectedGraph.contract(newNodes,lastNewNode,newNodes_pointers,oldNode2newNode,false);
		return lastNewNode+1;
	}
	

	
	
	
	
	
	
	// ------------------------------------- BASINS --------------------------------------
	
	/**
	 * The \emph{basin} of a kernel is the set of all intervals that are reachable from 
	 * kernel nodes by strict containment and identity arcs. Rather than the basin of a 
	 * kernel, the procedure builds the basin of every \emph{assembly path} of a kernel, 
	 * i.e. it assumes that every kernel node $v$ is annotated with the IDs of all 
	 * assembly paths of its bidirected graph that use node $v.bidirectedGraphNode$, and 
	 * it propagates such tags to other nodes in the interval graph. This gives distinct
	 * tags to e.g. distinct connected components of the bidirected graph, and to assembly
	 * paths that use different one-end nodes.
	 * 
	 * Remark: the same interval graph node can belong to multiple basins. An interval 
	 * graph node may also belong to no basin, e.g. if it belongs only to the basin of 
	 * kernels that have been simplified away by $buildBidirectedGraph()$.
	 *
	 * Remark: a kernel can have a basin of size one, if e.g. it consists of a single 
	 * maximal node whose all incident containment edges are of type 
	 * $CONTAINMENT_IDENTICAL$ and lead to nodes that are not maximal.
	 *
	 * Remark: the kernel tags in each kernel node are assumed to be sorted. The
	 * kernel tags in every node of a basin will be sorted as well.
	 *
	 * Remark: the procedure assumes that $nodesInKernel[0..nNodesInKernel-1]$ contains 
	 * the IDs of all kernel nodes (in no particular order).
	 *
	 * Remark: this procedure could be parallelized by working on each kernel 
	 * independently.
	 *
	 * @param resetBasins if TRUE the procedure preventively sets $lastKernel$ to -1 for
	 * all nodes of all basins (excluding kernel nodes).
	 */
	private static final void getBasins(int maxKernelsPerInterval, boolean resetBasins) {
		final int GROWTH_RATE = 30;  // Arbitrary 
		int i;
		int top, node, to;
		int[] tmp = new int[maxKernelsPerInterval];
		byte[] tmpBytes = new byte[maxKernelsPerInterval];
		
		// Setting $lastKernel=-1$ in all basins of all kernels.
		if (resetBasins) {
			for (i=0; i<nNodesInKernel; i++) {
				node=nodesInKernel[i];
				if (IntervalGraph.nNeighbors[node]==0) continue;
				top=getBasins_clean(node,false/*From a kernel node: only strict containment*/,IntervalGraph.stack,-1,-1);
				while (top>=0) top=getBasins_clean(IntervalGraph.stack[top],true,IntervalGraph.stack,top-1,-1);
			}
		}
		
		// Propagating the $kernels$ array of every node in a kernel.
		for (i=0; i<nNodesInKernel; i++) {
			node=nodesInKernel[i];
			if (IntervalGraph.nodesArray[node].lastKernel==-1 || IntervalGraph.nNeighbors[node]==0) continue;
			if (IntervalGraph.stack.length<3*IntervalGraph.nNeighbors[node]) {
				int[] newStack = new int[3*IntervalGraph.nNeighbors[node]+GROWTH_RATE];
				System.arraycopy(IntervalGraph.stack,0,newStack,0,IntervalGraph.stack.length);
				IntervalGraph.stack=newStack;
			}
			top=getBasins_impl(node,node,SOURCE_CONTAINMENT_IDENTITY,0,false/*From a kernel node: only strict containment*/,IntervalGraph.stack,-1,tmp,tmpBytes);
			while (top>=0) {
				to=IntervalGraph.stack[top-2];
				if (IntervalGraph.stack.length<(top+1)+3*IntervalGraph.nNeighbors[to]) {
					int[] newStack = new int[(top+1)+3*IntervalGraph.nNeighbors[to]+GROWTH_RATE];
					System.arraycopy(IntervalGraph.stack,0,newStack,0,IntervalGraph.stack.length);
					IntervalGraph.stack=newStack;
				}
				top=getBasins_impl(node,to,IntervalGraph.stack[top-1],IntervalGraph.stack[top],true,IntervalGraph.stack,top-3,tmp,tmpBytes);
			}
			if (i%100==0) System.err.println("getBasins> "+IO.getPercent(i,nNodesInKernel)+"% done");
		}
		
		// Ensures that only nodes with compatible length are assigned to a kernel tag
		getBasins_ensureKernelLengths();
	}
	
	
	/**
	 * @param to node to be processed;
	 * @param useIdenticalEdges uses also ON edges of type $CONTAINMENT_IDENTICAL$;
	 * @param top the current top element in $stack$;
	 * @param cleaningMode -1: deletes all kernel tags; $i>=0$: deletes just kernel tag 
	 * $i$;
	 * @return the new value of $top$ in $stack$.
	 */
	private static final int getBasins_clean(int to, boolean useIdenticalEdges, int[] stack, int top, int cleaningMode) {
		int i, j, p;
		int out, from;
		IntervalGraph.Node node;
		IntervalGraph.Edge edge;
		
		out=top;
		for (j=0; j<IntervalGraph.nNeighbors[to]; j++) {		
			edge=IntervalGraph.neighbors[to][j];
			if (!edge.on) break;
			from=edge.getTo(to);
			node=IntervalGraph.nodesArray[from];
			if ( node.isContained && 
				 ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==from && edge.nodeID2==to) ||
				   (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID2==from && edge.nodeID1==to) ||
				   (useIdenticalEdges && edge.containment==Constants.CONTAINMENT_IDENTICAL)
				 )
			   ) {
			   if (node.visited!=-1) {
				   node.visited=-1;
				   if (cleaningMode>=0) getBasins_clean_removeKernel(node,cleaningMode);
				   else {
					   node.lastKernel=-1;
					   node.lastPathWithStart=-1;
					   node.lastPathWithEnd=-1;
					   node.kernelOrientations=null;
				   }
				   IntervalGraph.stack[++out]=from;
			   }
			}
		}
		return out;
	}
	
	
	/**
	 * Removes just $kernel>=0$ from the kernel tags of $node$.
	 */
	private static final void getBasins_clean_removeKernel(IntervalGraph.Node node, int kernel) {
		int i, p;
		
		p=Arrays.binarySearch(node.kernels,0,node.lastKernel+1,kernel);
		if (p>=0) {
		   for (i=p; i<node.lastKernel; i++) node.kernels[i]=node.kernels[i+1];
		   if (node.kernelOrientations!=null) {
			   for (i=p; i<node.lastKernel; i++) node.kernelOrientations[i]=node.kernelOrientations[i+1];
		   }
		   node.lastKernel--;
		   if (node.lastPathWithStart>=0) {
			   p=Arrays.binarySearch(node.pathsWithStart,0,node.lastPathWithStart+1,kernel);
			   if (p>=0) {
				   for (i=p; i<node.lastPathWithStart; i++) node.pathsWithStart[i]=node.pathsWithStart[i+1];
				   node.lastPathWithStart--;
			   }
			   p=Arrays.binarySearch(node.pathsWithStart,0,node.lastPathWithStart+1,-1-kernel);
			   if (p>=0) {
				   for (i=p; i<node.lastPathWithStart; i++) node.pathsWithStart[i]=node.pathsWithStart[i+1];
				   node.lastPathWithStart--;
			   }
		   }
		   if (node.lastPathWithEnd>=0) {
			   p=Arrays.binarySearch(node.pathsWithEnd,0,node.lastPathWithEnd+1,kernel);
			   if (p>=0) {
				   for (i=p; i<node.lastPathWithEnd; i++) node.pathsWithEnd[i]=node.pathsWithEnd[i+1];
				   node.lastPathWithEnd--;
			   }
			   p=Arrays.binarySearch(node.pathsWithEnd,0,node.lastPathWithEnd+1,-1-kernel);
			   if (p>=0) {
				   for (i=p; i<node.lastPathWithEnd; i++) node.pathsWithEnd[i]=node.pathsWithEnd[i+1];
				   node.lastPathWithEnd--;
			   }
		   }
		}
	}
	
	
	/**
	 * Remark: the procedure propagates kernel tags also along supplement edges.
	 *
	 * @param to node to be processed;
	 * @param source node in a kernel;
	 * @param sourceContainmentTo containment type of $to$ WRT $source$;
	 * @param sourceOrientationTo orientation of $to$ WRT $source$;
	 * @param useIdenticalEdges uses also ON edges of type $CONTAINMENT_IDENTICAL$;
	 * @param top the current top element in $stack$;
	 * @param tmp,tmpBytes temporary space, of size at least equal to the maximum number 
	 * of kernels per interval;
	 * @return the new value of $top$ in $stack$.
	 */
	private static final int getBasins_impl(int source, int to, int sourceContainmentTo, int sourceOrientationTo, boolean useIdenticalEdges, int[] stack, int top, int[] tmp, byte[] tmpBytes) {
		final int GROWTH_RATE = 10;  // Arbitrary
		int j;
		int out, from, lastKernel, sourceContainment, sourceOrientation, containmentSubtype;
		IntervalGraph.Node sourceNode, fromNode;
		IntervalGraph.Edge edge;
		
		sourceNode=IntervalGraph.nodesArray[source];
//if (sourceNode.kernels[0]==7/*sourceNode.read==12359 && Math.abs(sourceNode.start,4672)<=100 && Math.abs(sourceNode.end,9895)<=100*/) System.err.println("VITTU> getBasins_impl() called on node "+to+" "+IntervalGraph.nodesArray[to]);
		out=top;
		for (j=0; j<IntervalGraph.nNeighbors[to]; j++) {
			edge=IntervalGraph.neighbors[to][j];
			if (!edge.on) break;
			from=edge.getTo(to);
			fromNode=IntervalGraph.nodesArray[from];
			if ( fromNode.isContained && 
				 ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==from && edge.nodeID2==to) ||
				   (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID2==from && edge.nodeID1==to) ||
				   (useIdenticalEdges && edge.containment==Constants.CONTAINMENT_IDENTICAL) 
				 )
			   ) {
				if (fromNode.visited==source) continue;
//if (sourceNode.kernels[0]==43/*sourceNode.read==12359 && Math.abs(sourceNode.start,4672)<=100 && Math.abs(sourceNode.end,9895)<=100*/) {
//	System.err.println("   explored node "+fromNode.read+"["+fromNode.start+".."+fromNode.end+"]");				
//	System.err.println("   EDGE: "+edge);
//}
//if (sourceNode.kernels[0]==7) System.err.println(to+" -> "+from+" [source=\""+source+"\",color=\""+(edge.containment==Constants.CONTAINMENT_IDENTICAL?"yellow":"blue")+"\"];");

			   	fromNode.visited=source;				
				sourceOrientation=getBasins_getSourceOrientation[sourceOrientationTo+1][edge.orientation+1];  // Can be -1
				if (edge.containment==Constants.CONTAINMENT_IDENTICAL) containmentSubtype=3;
				else if (edge.containmentSubtype==-1) containmentSubtype=2;
				else containmentSubtype=edge.containmentSubtype;
				sourceContainment=getBasins_getSourceContainment[sourceContainmentTo==Constants.CONTAINMENT_IDENTICAL?3:sourceContainmentTo][sourceOrientationTo+1][containmentSubtype];
				if (sourceContainment==3) sourceContainment=Constants.CONTAINMENT_IDENTICAL;
				getBasins_impl_lastPaths(fromNode,sourceNode,sourceContainment,sourceOrientation,tmp,GROWTH_RATE);  // Must be done before the next $continue$ for full generality.
				lastKernel=Math.setUnion(IntervalGraph.nodesArray[source].kernels,IntervalGraph.nodesArray[source].lastKernel,fromNode.kernels,fromNode.lastKernel,tmp);
				getBasins_impl_kernelOrientations(fromNode,sourceNode,sourceOrientation,tmp,lastKernel,tmpBytes);
			   	if (fromNode.kernelOrientations==null || fromNode.kernelOrientations.length<lastKernel+1) fromNode.kernelOrientations = new byte[lastKernel+1+GROWTH_RATE];
			   	System.arraycopy(tmpBytes,0,fromNode.kernelOrientations,0,lastKernel+1);
			   	if (lastKernel==fromNode.lastKernel) continue;
			   	fromNode.lastKernel=lastKernel;
			   	if (fromNode.kernels.length<lastKernel+1) fromNode.kernels = new int[lastKernel+1+GROWTH_RATE];
			   	System.arraycopy(tmp,0,fromNode.kernels,0,lastKernel+1);
			   	stack[++out]=from;
				stack[++out]=sourceContainment;
				stack[++out]=sourceOrientation;
			}
		}
		return out;
	}
	
	
	/**
	 * Decides how the neighbor of node $to$ along $edge$ is oriented WRT the source node
	 * of the propagation.
	 *
	 * Dimensions:
	 * 0: $sourceOrientationTo$: 0=unknown, 1=forward, 2=reverse, 3=both.
	 * 1: $edge.orientation$: 0=unknown, 1=forward, 2=reverse, 3=both.
	 *
	 * Value in a cell: -1=unknown, 0=forward, 1=reverse, 2=both.
	 */
	private static final int[][] getBasins_getSourceOrientation = new int[][] {
		{-1,-1,-1,-1}, {-1,0,1,2}, {-1,1,0,2}, {2,2,2,2}
	};
	
	
	/**
	 * Decides how the neighbor of node $to$ along $edge$ is contained in the source node
	 * of the propagation.
	 *
	 * Dimensions: 
	 * 0: $sourceContainmentTo$: 0=prefix, 1=suffix, 2=substring, 3=identical.
	 * 1: $sourceOrientationTo$: 0=unknown, 1=forward, 2=reverse, 3=both.
	 * 2: $edge.containmentSubtype$: 0=prefix, 1=suffix, 2=substring, 3=identical.
	 *
	 * Value in a cell: 0=prefix, 1=suffix, 2=substring, 3=identical.
	 */
	private static final int[][][] getBasins_getSourceContainment = new int[][][] {
		{ {2,2,2,2}, {0,2,2,0}, {2,0,2,0}, {2,2,2,0} },
	    { {2,2,2,2}, {2,1,2,1}, {1,2,2,1}, {2,2,2,1} },
		{ {2,2,2,2}, {2,2,2,2}, {2,2,2,2}, {2,2,2,2} },
		{ {2,2,2,2}, {0,1,2,3}, {1,0,2,3}, {2,2,2,3} }
	};
	
	
	/**
	 * Propagates the $pathsWith*,lastPathWith*$ fields of $sourceNode$ to $fromNode$.
	 *
	 * @param sourceContainment containment type of $fromNode$ WRT $sourceNode$;
	 * @param sourceOrientation orientation of $fromNode$ WRT $sourceNode$; the procedure
	 * does not do anything if $sourceOrientation \notin \{0,1,2\}$.
     */
	private static final void getBasins_impl_lastPaths(IntervalGraph.Node fromNode, IntervalGraph.Node sourceNode, int sourceContainment, int sourceOrientation, int[] tmp, int growthRate) {
		int last;

		if (sourceOrientation==0 || sourceOrientation==2) {
			if ( (sourceContainment==Constants.CONTAINMENT_PREFIX || sourceContainment==SOURCE_CONTAINMENT_IDENTITY) && 
			     sourceNode.lastPathWithStart>=0
			   ) {
				if (fromNode.lastPathWithStart==-1) {
					fromNode.lastPathWithStart=sourceNode.lastPathWithStart;
					if (fromNode.pathsWithStart==null || fromNode.pathsWithStart.length<=fromNode.lastPathWithStart) fromNode.pathsWithStart = new int[fromNode.lastPathWithStart+1+growthRate];
					System.arraycopy(sourceNode.pathsWithStart,0,fromNode.pathsWithStart,0,fromNode.lastPathWithStart+1);
				}
			   	else {
				   	last=Math.setUnion(sourceNode.pathsWithStart,sourceNode.lastPathWithStart,fromNode.pathsWithStart,fromNode.lastPathWithStart,tmp);
				   	if (last>fromNode.lastPathWithStart) {
					   	fromNode.lastPathWithStart=last;
					   	if (fromNode.pathsWithStart.length<last+1) fromNode.pathsWithStart = new int[last+1+growthRate];
					   	System.arraycopy(tmp,0,fromNode.pathsWithStart,0,last+1);
				   	}
			   	}
		   	}
		   	else if ( (sourceContainment==Constants.CONTAINMENT_SUFFIX || sourceContainment==SOURCE_CONTAINMENT_IDENTITY) && 
			          sourceNode.lastPathWithEnd>=0
			        ) {
			   	if (fromNode.lastPathWithEnd==-1) {
				   	fromNode.lastPathWithEnd=sourceNode.lastPathWithEnd;
				   	if (fromNode.pathsWithEnd==null || fromNode.pathsWithEnd.length<=fromNode.lastPathWithEnd) fromNode.pathsWithEnd = new int[fromNode.lastPathWithEnd+1+growthRate];
				   	System.arraycopy(sourceNode.pathsWithEnd,0,fromNode.pathsWithEnd,0,fromNode.lastPathWithEnd+1);
			   	}
			   	else {
				   	last=Math.setUnion(sourceNode.pathsWithEnd,sourceNode.lastPathWithEnd,fromNode.pathsWithEnd,fromNode.lastPathWithEnd,tmp);
				   	if (last>fromNode.lastPathWithEnd) {
					   	fromNode.lastPathWithEnd=last;
					   	if (fromNode.pathsWithEnd.length<last+1) fromNode.pathsWithEnd = new int[last+1+growthRate];
						System.arraycopy(tmp,0,fromNode.pathsWithEnd,0,last+1);
				   	}
			   	}
		   	}
	   	}
		if (sourceOrientation==1 || sourceOrientation==2) {
			if ( (sourceContainment==Constants.CONTAINMENT_PREFIX || sourceContainment==SOURCE_CONTAINMENT_IDENTITY) && 
			     sourceNode.lastPathWithStart>=0
			   ) {
				if (fromNode.lastPathWithEnd==-1) {
					fromNode.lastPathWithEnd=sourceNode.lastPathWithStart;
					if (fromNode.pathsWithEnd==null || fromNode.pathsWithEnd.length<sourceNode.lastPathWithStart+1) fromNode.pathsWithEnd = new int[sourceNode.lastPathWithStart+1+growthRate];
					System.arraycopy(sourceNode.pathsWithStart,0,fromNode.pathsWithEnd,0,sourceNode.lastPathWithStart+1);
				}
			   	else {
				   	last=Math.setUnion(sourceNode.pathsWithStart,sourceNode.lastPathWithStart,fromNode.pathsWithEnd,fromNode.lastPathWithEnd,tmp);
				   	if (last>fromNode.lastPathWithEnd) {
					   	fromNode.lastPathWithEnd=last;
					   	if (fromNode.pathsWithEnd.length<last+1) fromNode.pathsWithEnd = new int[last+1+growthRate];
					   	System.arraycopy(tmp,0,fromNode.pathsWithEnd,0,last+1);
				   	}
			   	}
		   	}
		   	else if ( (sourceContainment==Constants.CONTAINMENT_SUFFIX || sourceContainment==SOURCE_CONTAINMENT_IDENTITY) && 
			          sourceNode.lastPathWithEnd>=0
			        ) {
			   	if (fromNode.lastPathWithStart==-1) {
				   	fromNode.lastPathWithStart=sourceNode.lastPathWithEnd;
				   	if (fromNode.pathsWithStart==null || fromNode.pathsWithStart.length<sourceNode.lastPathWithEnd+1) fromNode.pathsWithStart = new int[sourceNode.lastPathWithEnd+1+growthRate];
				   	System.arraycopy(sourceNode.pathsWithEnd,0,fromNode.pathsWithStart,0,sourceNode.lastPathWithEnd+1);
			   	}
			   	else {
				   	last=Math.setUnion(sourceNode.pathsWithEnd,sourceNode.lastPathWithEnd,fromNode.pathsWithStart,fromNode.lastPathWithStart,tmp);
				   	if (last>fromNode.lastPathWithStart) {
					   	fromNode.lastPathWithStart=last;
					   	if (fromNode.pathsWithStart.length<last+1) fromNode.pathsWithStart = new int[last+1+growthRate];
					   	System.arraycopy(tmp,0,fromNode.pathsWithStart,0,last+1);
				   	}
			   	}
		   	}
	   	}
	}
	
	
	/**
	 * Let $union[0..lastUnion]$ be the union of the $kernels$ field of $fromNode$ and 
	 * $sourceNode$. The procedure stores in $out[0..lastUnion]$ the $kernelOrientations$ 
	 * values of all elements in $union$.
	 *
	 * @param sourceOrientation orientation of $fromNode$ WRT $sourceNode$; if $\notin
	 * $\{0,1,2\}$, the corresponding positions of $out$ are set to -1.
     */
	private static final void getBasins_impl_kernelOrientations(IntervalGraph.Node fromNode, IntervalGraph.Node sourceNode, int sourceOrientation, int[] union, int lastUnion, byte[] out) {
		int i, j;
		
		Math.set(out,lastUnion,(byte)(-1));
		
		// Copying the old orientations of $fromNode$.
		if (fromNode.kernelOrientations!=null) {
			i=0; j=0;
			while (i<=lastUnion && j<=fromNode.lastKernel) {
				if (union[i]<fromNode.kernels[j]) {
					i++;
					continue;
				}
				else if (union[i]>fromNode.kernels[j]) {
					j++;
					continue;
				}
				out[i]=fromNode.kernelOrientations[j];
				i++; j++;
			}
		}
		
		// Adding the orientations of $sourceNode$.
		if (sourceOrientation==0) {
			i=0; j=0;
			while (i<=sourceNode.lastKernel && j<=lastUnion) {
				if (sourceNode.kernels[i]<union[j]) {
					i++;
					continue;
				}
				else if (sourceNode.kernels[i]>union[j]) {
					j++;
					continue;
				}
				if (out[j]==-1) out[j]=sourceNode.kernelOrientations[i];
				else out[j]=addOrientation[out[j]][sourceNode.kernelOrientations[i]];
				i++; j++;
			}
		}
		else if (sourceOrientation==1) {
			i=0; j=0;
			while (i<=sourceNode.lastKernel && j<=lastUnion) {
				if (sourceNode.kernels[i]<union[j]) {
					i++;
					continue;
				}
				else if (sourceNode.kernels[i]>union[j]) {
					j++;
					continue;
				}
				if (out[j]==-1) out[j]=oppositeOrientation[sourceNode.kernelOrientations[i]];
				else out[j]=addOrientation[out[j]][oppositeOrientation[sourceNode.kernelOrientations[i]]];
				i++; j++;
			}
		}
		else if (sourceOrientation==2) {
			i=0; j=0;
			while (i<=sourceNode.lastKernel && j<=lastUnion) {
				if (sourceNode.kernels[i]<union[j]) {
					i++;
					continue;
				}
				else if (sourceNode.kernels[i]>union[j]) {
					j++;
					continue;
				}
				out[j]=2;
				i++; j++;
			}
		}
	}
	
	
	/**
	 * Cell $i$ contains orientation $i$ taken in the opposite direction.
	 * Orientation values: 0=forward, 1=RC, 2=both.
	 */
	private static final byte[] oppositeOrientation = {1,0,2};
	
	
	/**
	 * Cell $(i,j)$ contains the result of adding orientation $i$ and orientation $j$.
	 * Orientation values: 0=forward, 1=RC, 2=both.
	 */
	public static final byte[][] addOrientation = new byte[][] {
		{0,2,2}, {2,1,2}, {2,2,2}	
	};
	
	
	/**
	 * Removes a path kernel tag from a node iff the string length of the node is greater
	 * than the string length of the path kernel. This might happen because e.g. the 
	 * identity edges that are followed by $getBasins()$ use a distance threshold, and
	 * identity paths from a maximal node to a node in the basin might be hundreds of
	 * nodes long (this has been observed in practice), thus distance thresholds might
	 * add up arbitrarily. The whole $IntervalGraphStep3$ assumes that the kernel paths 
	 * built by assembling maximal nodes by containment are close in length and sequence 
	 * to the underlying repeat modules, and that all nodes in a basin are smaller pieces
	 * of such kernel paths.
	 *
	 * Remark: the length criterion should be used after the propagation of $getBasins()$ 
	 * has completed, rather than to stop the propagation; otherwise, some nodes with 
	 * consistent length might not get the kernel tag (since length is not monotonic along
	 * identity edges).
	 *
	 * Remark: the procedure might make a node have a kernel tag that none of its 
	 * ancestors by containment has. Fixing this inconsistency takes too much effort:* we 
	 * don't do it for simplicity. However, it might be the reason why a node is reported 
	 * by $FragmentsStep1$.
	 *
	 * Remark: there might be a repeat that replicates \emph{by expansion}, i.e. every 
	 * copy is a slightly longer version of a previous copy. This might create long
	 * paths of identity edges in a basin, and the procedure would remove the tag after
	 * a few steps. We don't model such a replication mode anywhere in the pipeline.
	 *
	 * @return the number of nodes that have been stripped of all kernel tags by the 
	 * procedure. This numer is negligible in practice.
	 */
	private static final int getBasins_ensureKernelLengths() {
		final int LENGTH_THRESHOLD = IO.quantum<<2;  // Arbitrary
		final int nNodes = IntervalGraph.nNodes;
		int i, j, k;
		int last, nodeLength, kernel, kernelPrime, nNoKernel;
		IntervalGraph.Node node;
		
		nNoKernel=0;
		for (i=0; i<nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			last=node.lastKernel;
			if (last==-1) continue;
			nodeLength=node.length();
			k=-1;
			for (j=0; j<=last; j++) {
				kernel=node.kernels[j];
				if ( pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0) ||
				     nodeLength<=pathKernelLengths[kernel]+LENGTH_THRESHOLD
				   ) {
				   k++;
				   node.kernels[k]=kernel;
				   node.kernelOrientations[k]=node.kernelOrientations[j];
			   }
			}
			node.lastKernel=k;
			if (k==-1) {
				node.lastPathWithStart=-1; node.lastPathWithEnd=-1;
				nNoKernel++;
				continue;
			}
			last=node.lastPathWithStart;
			k=-1;
			for (j=0; j<=last; j++) {
				kernel=node.pathsWithStart[j];
				kernelPrime=kernel>=0?kernel:-1-kernel;
				if ( pathKernelLengths[kernelPrime]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernelPrime]!=0) ||
				     nodeLength<=pathKernelLengths[kernelPrime]+LENGTH_THRESHOLD
				   ) node.pathsWithStart[++k]=kernel;
			}
			node.lastPathWithStart=k;
			last=node.lastPathWithEnd;
			k=-1;
			for (j=0; j<=last; j++) {
				kernel=node.pathsWithEnd[j];
				kernelPrime=kernel>=0?kernel:-1-kernel;
				if ( pathKernelLengths[kernelPrime]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernelPrime]!=0) ||
				     nodeLength<=pathKernelLengths[kernelPrime]+LENGTH_THRESHOLD
				   ) node.pathsWithEnd[++k]=kernel;
			}
			node.lastPathWithEnd=k;
		}
		if (nNoKernel!=0) System.err.println("getBasins_ensureKernelLengths> The procedure removed all kernel tags from "+nNoKernel+" nodes ("+IO.format(100*((double)nNoKernel)/nNodes)+"%).");
		return nNoKernel;
	}
	
	
	/**
	 * Computes $kernelSize$ and $kernelFrequency$ by iterating over interval graph nodes.
	 * See $getFrequency()$ for more details.
	 *
	 * Remark: if the global variable $descendants!=NULL$, the frequency of kernel paths
	 * is corrected using that information (see $getFrequency()$ for details).
	 *
	 * @param checkKernelsWithZeroNodes used only for consistency checks;
	 * @param out output array: 0=max number of path kernels to which an interval graph 
	 * node belongs; 1=number of nodes with one kernel; 2=number of nodes with no kernel.
	 */
	public static final void basinStats(int[] out, boolean checkKernelsWithZeroNodes) {
		int i, j;
		int from, read, lastKernel, max, nNodesWithOneKernel, nNodesWithNoKernel;
		IntervalGraph.Node tmpNode;
		boolean[] used = null;
		
		Math.set(kernelSize,nPathKernels-1,0);
		Math.set(kernelFrequency,0);
		Math.set(out,2,0);
		if (IntervalGraph.nNodes==0) return;
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		if (tmpArray2==null || tmpArray2.length<nPathKernels) tmpArray2 = new int[nPathKernels];
		if (tmpArray3==null || tmpArray3.length<nPathKernels) tmpArray3 = new int[nPathKernels];
		if (tmpArray4==null || tmpArray4.length<nPathKernels) tmpArray4 = new int[nPathKernels];
		max=0; nNodesWithOneKernel=0; nNodesWithNoKernel=0;
		from=0; read=IntervalGraph.nodesArray[0].read;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			lastKernel=tmpNode.lastKernel;
			if (lastKernel==-1) nNodesWithNoKernel++;
			else if (lastKernel==0) nNodesWithOneKernel++;
			if (lastKernel>max) max=lastKernel;
			for (j=0; j<=lastKernel; j++) kernelSize[tmpNode.kernels[j]]++;
			if (tmpNode.read!=read) {
				if (used==null || used.length<i-from) used = new boolean[(i-from)<<1];
				getFrequency(read,from,i-1,used);
				read=tmpNode.read;
				from=i;
			}
		}
		if (used==null || used.length<IntervalGraph.nNodes-from) used = new boolean[IntervalGraph.nNodes-from];
		getFrequency(read,from,IntervalGraph.nNodes-1,used);
		out[0]=max+1; out[1]=nNodesWithOneKernel; out[2]=nNodesWithNoKernel;
		
		// Consistency check
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<nPathKernels; i++) {
				if (checkKernelsWithZeroNodes && kernelSize[i]==0) {
					System.err.println("ERROR: path kernel "+i+" has size zero");
					System.exit(1);
				}
				if (pathKernelLengths[i]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[i]==0) && kernelSize[i]==1 && kernelFrequency[i][0]==0) {
					System.err.println("ERROR: path kernel "+i+" consists of just one interval graph node but has no full occurrence?!");
					System.err.println("Nodes in path kernel "+i+":  pathKernelLength="+pathKernelLengths[i]);
					for (int x=0; x<nNodesInKernel; x++) {
						tmpNode=IntervalGraph.nodesArray[nodesInKernel[x]];
						if (tmpNode.lastKernel>=0 && Arrays.binarySearch(tmpNode.kernels,0,tmpNode.lastKernel+1,i)>=0) {
							System.err.println("node: "+tmpNode);
							System.err.print("pathsWithStart: ");
							for (int y=0; y<=tmpNode.lastPathWithStart; y++) System.err.print(tmpNode.pathsWithStart[y]+",");
							System.err.println();
							System.err.print("pathsWithEnd: ");
							for (int y=0; y<=tmpNode.lastPathWithEnd; y++) System.err.print(tmpNode.pathsWithEnd[y]+",");
							System.err.println();
						}
					}
					System.exit(1);
				}
			}
		}
	}
	
	
	/**
	 * Stores in $kernelFrequency[i]$ an estimate of the number of times the following 
	 * parts of non-periodic kernel $i$ are observed in a read: 0=full copy of the kernel;
	 * 1=just a prefix of the kernel; 2=just a suffix of the kernel; 3=fragments that are 
	 * neither a prefix nor a suffix.
	 * A kernel $K$ is assumed to occur in full iff we observe a sequence of interval 
	 * graph nodes $N_1,N_2,...,N_n$ such that the beginning of $N_1$ maps to one end of 
	 * $K$, the end of $N_n$ maps to the other end of $K$, no two intervals $N_i,N_j$
	 * overlap, the substrings of the read between such intervals (if any) are of low 
	 * quality, and the distance between $N_1$ and $N_n$ matches the string length of $K$.
	 * Another meaning of "full copy" might be that the repeat occurs in full inside a
	 * read, but with insertions of other repeats: this is not supported here.
	 *
	 * The counts of prefix/suffix occurrences are incremented by just looking at the
	 * $pathsWith*$ arrays of single nodes. We could be more careful and store, in each
	 * interval graph node, for each one of its kernels, a pointer to a corresponding
	 * bidirected graph node. Then, we could try to extend every kernel in, e.g., the
	 * $pathsWithStart$ array of a node, until we reach a sequence of bidirected graph 
	 * nodes that is compatible with the kernel path, and uniquely identifies it, and we
	 * could increment just the count of such kernel path (in case of ambiguity, we could 
	 * increment the count of no kernel path). This would also require storing the full 
	 * sequence of bidirected graph nodes for each kernel path. We omit such approach for 
	 * simplicity.
	 * 		
	 * Remark: the procedure handles the case in which distinct kernel paths have a full 
	 * copy that starts from the same interval graph node. In this case, the full copy is
	 * not assigned to any kernel path. This is possible but not likely, since different 
	 * bidirected graph nodes are not expected to share high-quality substrings, thus 
	 * interval graph nodes that map to them should not overlap in a read.
     *
	 * Remark: the procedure assumes low-quality regions to be replacements, not 
	 * insertions.
	 *
	 * Remark: if the global variable $descendants!=NULL$, non-full frequencies are 
	 * incremented only for the kernels $X$ of a node such that no other kernel in the 
	 * same node belongs to $descendants[X]$.
	 *
	 * Remark: the procedure tries to handle the case in which an interval graph node 
	 * occurs at the beginning, end, and in the middle of a kernel path at the same time, 
	 * possibly reverse-complemented.
	 * 
	 * Remark: nodes without $pathsWith*$ arrays (e.g. those of cyclic kernels) contribute
	 * only to the substring count (i.e. they do not contribute to the full/prefix/suffix 
	 * counts).
	 *
	 * Remark: the procedure assumes arrays $pathsWith*$ of every interval graph node to 
	 * be sorted.
	 *
	 * @param from,to all interval graph nodes in the same read are in 
	 * $IntervalGraph.nodesArray[from..to]$, sorted by $start$;
	 * @param used temporary array, of size at least $to-from+1$.
	 */
	private static final void getFrequency(int read, int from, int to, boolean[] used) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int OVERLAP_THRESHOLD = IO.quantum;
		final int LENGTH_THRESHOLD = IO.quantum<<1;  // Arbitrary
		boolean isFull, found;
		int i, j, k, h;
		int kernel, nFull, fullKernel, lastNodeInFull, end, otherEnd, previousEnd;
		int last, last1, last5, lastPair, position;
		IntervalGraph.Node firstNode, lastNode, middleNode, previousNode;
		
		if (tmpArray5==null || tmpArray5.length<nPathKernels) tmpArray5 = new int[nPathKernels];
		Math.set(used,to-from,false);
		Math.set(tmpArray1,nPathKernels-1,0);		
		
		// Full occurrences
		for (i=from; i<=to; i++) {	
			if (used[i-from]) continue;
			firstNode=IntervalGraph.nodesArray[i];			
			nFull=0; fullKernel=-1; lastNodeInFull=-1;
			for (j=0; j<=firstNode.lastPathWithStart; j++) {
				end=firstNode.pathsWithStart[j];
				otherEnd=-1-firstNode.pathsWithStart[j];
				kernel=end>=0?end:-1-end;
				if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
				if (tmpArray1[kernel]==1) {
					continue;  // The kernel occurs in full in the other direction
				}
				for (k=i; k<=to; k++) {
					lastNode=IntervalGraph.nodesArray[k];
					if (lastNode.start>=firstNode.start+pathKernelLengths[kernel]) break;  
					if (used[k-from]) continue;
					if (lastNode.lastPathWithEnd==-1 || Arrays.binarySearch(lastNode.pathsWithEnd,0,lastNode.lastPathWithEnd+1,otherEnd)<0) continue;
					if (Math.abs(lastNode.end-firstNode.start+1,pathKernelLengths[kernel])>LENGTH_THRESHOLD) {
						continue;  // This is correct also when $k==i$.
					}					
					isFull=true; previousEnd=firstNode.end;
					for (h=i+1; h<=k; h++) {
						middleNode=IntervalGraph.nodesArray[h];
						if (middleNode.lastKernel==-1 || Arrays.binarySearch(middleNode.kernels,0,middleNode.lastKernel+1,kernel)<0) continue;
						if (middleNode.start>previousEnd+OVERLAP_THRESHOLD && !Reads.hasLowQuality(read,previousEnd+1,middleNode.start-1,true)) {
							isFull=false;
							break;
						}
						previousEnd=Math.max(previousEnd,middleNode.end);
					}
					if (isFull) {
						nFull++;
						fullKernel=kernel;
						lastNodeInFull=k;
						tmpArray1[kernel]=1;
						break;
					}
				}
				if (nFull>1) break;
			}
			if (nFull==1) {
				kernelFrequency[fullKernel][0]++;
				used[i-from]=true;
				for (j=i+1; j<=lastNodeInFull; j++) {
					if (used[j-from]) continue;
					middleNode=IntervalGraph.nodesArray[j];
					if (Arrays.binarySearch(middleNode.kernels,0,middleNode.lastKernel+1,fullKernel)>=0) used[j-from]=true;
				}
			}
			// Cleaning up $tmpArray1$ for the next iteration.
			for (j=0; j<=firstNode.lastPathWithStart; j++) {
				end=firstNode.pathsWithStart[j];
				kernel=end>=0?end:-1-end;
				tmpArray1[kernel]=0;
			}
		}
		
		// Prefixes, suffixes, substrings. If an interval graph node occurs both at one
		// end and in the middle of a kernel, it contributes to the count of that end even
		// if the occurrence we are considering is in reality the one in the middle.
		lastPair=-1;
		Math.set(tmpArray4,nPathKernels-1,0);
		for (i=from; i<=to; i++) {
			firstNode=IntervalGraph.nodesArray[i];
			if (used[i-from]) continue;
			if (firstNode.lastKernel==-1) continue;
			last5=startUnionEnd(firstNode,tmpArray1,tmpArray2,tmpArray5);
			last=startIntersectEnd(firstNode,tmpArray1,tmpArray2,tmpArray3);
			if (last>=0) {
				for (j=0; j<=last; j++) {
					kernel=tmpArray3[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
					if ( descendants!=null && 
						 ( Math.nonemptyIntersection(tmpArray5,0,last5,descendants[kernel],0,lastDescendant[kernel],kernel) ||
						   getFrequency_inDescendant(i,firstNode,kernel,true,tmpArray1) ||
						   getFrequency_inDescendant(i,firstNode,kernel,false,tmpArray1)
						 )
					   ) continue;
					if (tmpArray2[j]==0) {
						if (!getFrequency_inKernel(firstNode.start,read,kernel,i)) lastPair=getFrequency_addFrequencyPair(firstNode.start,-1-kernel/*Incrementing just the end of the kernel. Arbitrary.*/,lastPair);
						else if (!getFrequency_inKernel(firstNode.end,read,kernel,i)) lastPair=getFrequency_addFrequencyPair(firstNode.end,-1-kernel/*Incrementing just the end of the kernel. Arbitrary.*/,lastPair);
					}
					else if (tmpArray2[j]==1) {
						if (!getFrequency_inKernel(firstNode.start,read,kernel,i)) lastPair=getFrequency_addFrequencyPair(firstNode.start,kernel,lastPair);
						else if (!getFrequency_inKernel(firstNode.end,read,kernel,i)) lastPair=getFrequency_addFrequencyPair(firstNode.end,kernel,lastPair);
					}
					else if (tmpArray2[j]==2) {
						if (!getFrequency_inKernel(firstNode.start,read,kernel,i)) lastPair=getFrequency_addFrequencyPair(firstNode.start,-1-kernel,lastPair);
						else if (!getFrequency_inKernel(firstNode.end,read,kernel,i)) lastPair=getFrequency_addFrequencyPair(firstNode.end,-1-kernel,lastPair);
					}
				}
				System.arraycopy(firstNode.pathsWithStart,0,tmpArray1,0,firstNode.lastPathWithStart+1);
				last1=Math.makePositive(tmpArray1,0,firstNode.lastPathWithStart);
				h=firstNode.lastPathWithStart;  // Last negative element
				for (j=0; j<=firstNode.lastPathWithStart; j++) {
					if (firstNode.pathsWithStart[j]>=0) {
						h=j-1;
						break;
					}
				}
				j=h; k=0;
				while (j>=0 && k<=last) {
					kernel=-1-firstNode.pathsWithStart[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
						j--;
						continue;
					}
					if (kernel<tmpArray3[k]) {
						if ( ( descendants==null || 
						       ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
							     !getFrequency_inDescendant(i,firstNode,kernel,true,tmpArray2)
							   )
							 ) &&
							 !getFrequency_inKernel(firstNode.start,read,kernel,i)
						   ) {
							lastPair=getFrequency_addFrequencyPair(firstNode.start,-1-kernel,lastPair);
							tmpArray4[kernel]=1;
						}
						j--;
						continue;
					}
					else if (kernel>tmpArray3[k]) {
						k++;
						continue;
					}
					j--; k++;
				}
				while (j>=0) {
					kernel=-1-firstNode.pathsWithStart[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
						j--;
						continue;
					}
					if ( ( descendants==null || 
					       ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
						     !getFrequency_inDescendant(i,firstNode,kernel,true,tmpArray2)
						   )
						 ) &&
						 !getFrequency_inKernel(firstNode.start,read,kernel,i)
					   ) {
						lastPair=getFrequency_addFrequencyPair(firstNode.start,-1-kernel,lastPair);
						tmpArray4[kernel]=1;
					}
					j--;
				}
				j=h+1; k=0;
				while (j<=firstNode.lastPathWithStart && k<=last) {
					kernel=firstNode.pathsWithStart[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
						j++;
						continue;
					}
					if (kernel<tmpArray3[k]) {
						if ( tmpArray4[kernel]==0 && 
				             ( descendants==null || 
							   ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
								 !getFrequency_inDescendant(i,firstNode,kernel,true,tmpArray2)
							   )
							 ) &&
							 !getFrequency_inKernel(firstNode.start,read,kernel,i)
					       ) lastPair=getFrequency_addFrequencyPair(firstNode.start,kernel,lastPair);
						j++;
						continue;
					}
					else if (kernel>tmpArray3[k]) {
						k++;
						continue;
					}
					j++; k++;
				}
				while (j<=firstNode.lastPathWithStart) {
					kernel=firstNode.pathsWithStart[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
						j++;
						continue;
					}
					if ( tmpArray4[kernel]==0 && 
						 ( descendants==null || 
						   ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
							 !getFrequency_inDescendant(i,firstNode,kernel,true,tmpArray2)
						   )
						 ) &&
						 !getFrequency_inKernel(firstNode.start,read,kernel,i)
					   ) lastPair=getFrequency_addFrequencyPair(firstNode.start,kernel,lastPair);
					j++;
				}
				for (j=0; j<=firstNode.lastPathWithStart; j++) {
					kernel=firstNode.pathsWithStart[j];
					if (kernel<0) kernel=-1-kernel;
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
					tmpArray4[kernel]=0;
				}
				System.arraycopy(firstNode.pathsWithEnd,0,tmpArray1,0,firstNode.lastPathWithEnd+1);
				last1=Math.makePositive(tmpArray1,0,firstNode.lastPathWithEnd);
				h=firstNode.lastPathWithEnd;  // Last negative element
				for (j=0; j<=firstNode.lastPathWithEnd; j++) {
					if (firstNode.pathsWithEnd[j]>=0) {
						h=j-1;
						break;
					}
				}
				j=h; k=0;
				while (j>=0 && k<=last) {
					kernel=-1-firstNode.pathsWithEnd[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
						j--;
						continue;
					}
					if (kernel<tmpArray3[k]) {
						if ( ( descendants==null || 
						       ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
							     !getFrequency_inDescendant(i,firstNode,kernel,false,tmpArray2)
							   )
							 ) &&
							 !getFrequency_inKernel(firstNode.end,read,kernel,i)
						   ) {
							lastPair=getFrequency_addFrequencyPair(firstNode.end,-1-kernel,lastPair);
							tmpArray4[kernel]=1;
						}
						j--;
						continue;
					}
					else if (kernel>tmpArray3[k]) {
						k++;
						continue;
					}
					j--; k++;
				}
				while (j>=0) {
					kernel=-1-firstNode.pathsWithEnd[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
						j--;
						continue;
					}
					if ( ( descendants==null || 
					       ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
						     !getFrequency_inDescendant(i,firstNode,kernel,false,tmpArray2)
						   )
						 ) &&
					     !getFrequency_inKernel(firstNode.end,read,kernel,i)
					   ) {
						lastPair=getFrequency_addFrequencyPair(firstNode.end,-1-kernel,lastPair);
						tmpArray4[kernel]=1;
					}
					j--;
				}
				j=h+1; k=0;
				while (j<=firstNode.lastPathWithEnd && k<=last) {
					kernel=firstNode.pathsWithEnd[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
						j++;
						continue;
					}
					if (kernel<tmpArray3[k]) {
						if ( tmpArray4[kernel]==0 && 
							 ( descendants==null || 
							   ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
								 !getFrequency_inDescendant(i,firstNode,kernel,false,tmpArray2)
							   )
							 ) &&
							 !getFrequency_inKernel(firstNode.end,read,kernel,i)
						   ) lastPair=getFrequency_addFrequencyPair(firstNode.end,kernel,lastPair);
						j++;
						continue;
					}
					else if (kernel>tmpArray3[k]) {
						k++;
						continue;
					}
					j++; k++;
				}
				while (j<=firstNode.lastPathWithEnd) {
					kernel=firstNode.pathsWithEnd[j];
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
						j++;
						continue;
					}
					if ( tmpArray4[kernel]==0 &&
						 ( descendants==null || 
						   ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
							 !getFrequency_inDescendant(i,firstNode,kernel,false,tmpArray2)
						   )
						 ) &&
						 !getFrequency_inKernel(firstNode.end,read,kernel,i)
					   ) lastPair=getFrequency_addFrequencyPair(firstNode.end,kernel,lastPair);
					j++;
				}
				for (j=0; j<=firstNode.lastPathWithEnd; j++) {
					kernel=firstNode.pathsWithEnd[j];
					if (kernel<0) kernel=-1-kernel;
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
					tmpArray4[kernel]=0;
				}
			}
			else {
				if (firstNode.lastPathWithStart>=0) {
					System.arraycopy(firstNode.pathsWithStart,0,tmpArray1,0,firstNode.lastPathWithStart+1);
					last1=Math.makePositive(tmpArray1,0,firstNode.lastPathWithStart);
					for (j=0; j<=firstNode.lastPathWithStart; j++) {
						end=firstNode.pathsWithStart[j];
						kernel=end>=0?end:-1-end;
						if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
						if ( ( descendants==null || 
						       ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
							     !getFrequency_inDescendant(i,firstNode,kernel,true,tmpArray2)
							   )
							 ) &&
							 !getFrequency_inKernel(firstNode.start,read,kernel,i)
						   ) {
							if (end>=0) lastPair=getFrequency_addFrequencyPair(firstNode.start,kernel,lastPair);
							else lastPair=getFrequency_addFrequencyPair(firstNode.start,-1-kernel,lastPair);
						}
					}
				}
				if (firstNode.lastPathWithEnd>=0) {
					System.arraycopy(firstNode.pathsWithEnd,0,tmpArray1,0,firstNode.lastPathWithEnd+1);
					last1=Math.makePositive(tmpArray1,0,firstNode.lastPathWithEnd);
					for (j=0; j<=firstNode.lastPathWithEnd; j++) {
						end=firstNode.pathsWithEnd[j];
						kernel=end>=0?end:-1-end;
						if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
						if ( ( descendants==null || 
						       ( !Math.nonemptyIntersection(tmpArray1,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel) &&
							     !getFrequency_inDescendant(i,firstNode,kernel,false,tmpArray2)
							   )
							 ) &&
							 !getFrequency_inKernel(firstNode.end,read,kernel,i)
						   ) {
							if (end>=0) lastPair=getFrequency_addFrequencyPair(firstNode.end,kernel,lastPair);
							else lastPair=getFrequency_addFrequencyPair(firstNode.end,-1-kernel,lastPair);
						}
					}
				}
			}
			last=neitherStartNorEnd(firstNode,tmpArray1,tmpArray2,tmpArray3);
			for (j=0; j<=last; j++) {
				kernel=tmpArray2[j];
				if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
				if (descendants==null || !Math.nonemptyIntersection(firstNode.kernels,0,firstNode.lastKernel,descendants[kernel],0,lastDescendant[kernel],kernel)) {
					kernelFrequency[kernel][3]++;
				}
			}
		}
		
		// Updating prefix/suffix frequencies with distinct start/end positions.
		if (lastPair==-1) return;
		if (lastPair>0) Arrays.sort(frequencyPairs,0,lastPair+1);
		kernel=frequencyPairs[0].kernel; position=frequencyPairs[0].position;	
		for (i=1; i<=lastPair; i++) {
			if (frequencyPairs[i].kernel!=kernel) {
				if (kernel>=0) kernelFrequency[kernel][1]++;
				else kernelFrequency[-1-kernel][2]++;
				kernel=frequencyPairs[i].kernel; position=frequencyPairs[i].position;
			}
			else if (frequencyPairs[i].position>position+IDENTITY_THRESHOLD) {
				if (kernel>=0) kernelFrequency[kernel][1]++;
				else kernelFrequency[-1-kernel][2]++;
				position=frequencyPairs[i].position;
			}
		}
		if (kernel>=0) kernelFrequency[kernel][1]++;
		else kernelFrequency[-1-kernel][2]++;		
	}
	
	
	/**
	 * Appends a frequency pair to the global variable $frequencyPairs$.
	 *
	 * @return the new value of $lastPair$.
	 */
	private static final int getFrequency_addFrequencyPair(int position, int kernel, int lastPair) {
		lastPair++;
		if (lastPair==frequencyPairs.length) {
			FrequencyPair[] newPairs = new FrequencyPair[frequencyPairs.length<<1];
			System.arraycopy(frequencyPairs,0,newPairs,0,frequencyPairs.length);
			for (int i=frequencyPairs.length; i<newPairs.length; i++) newPairs[i] = new FrequencyPair();
			frequencyPairs=newPairs;
		}
		frequencyPairs[lastPair].position=position;
		frequencyPairs[lastPair].kernel=kernel;
		return lastPair;	
	}
	
	
	public static class FrequencyPair implements Comparable {
		public int position;
		public int kernel;  // $k$: start; $-1-k$: end.
		
		public int compareTo(Object other) {
			FrequencyPair otherPair = (FrequencyPair)other;
			if (kernel<otherPair.kernel) return -1;
			else if (kernel>otherPair.kernel) return 1;
			if (position<otherPair.position) return -1;
			else if (position>otherPair.position) return 1;
			return 0;
		}
		
		public String toString() { return "("+position+","+kernel+")"; }
	}
	
	
	/**
	 * @param startOrEnd the procedure works on a position $p$ that is the start (TRUE) or
	 * the end (FALSE) of $node$;
	 * @param tmpArray1 temporary space;
	 * @return TRUE iff there is an interval graph node $v$ different from $node$ such 
	 * that: (1) the start of $v$ is close to the start of $node$, or the end of $v$ is 
	 * close to the end of $node$; (2) an end of a descendant of $kernel$ occurs at the 
	 * end of $v$ that is close to $p$.
	 */
	private static final boolean getFrequency_inDescendant(int nodeIndex, IntervalGraph.Node node, int kernel, boolean startOrEnd, int[] tmpArray) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i;
		int last1;
		final int nNodes = IntervalGraph.nNodes;
		final int nodeRead = node.read;
		final int nodeStart = node.start;
		final int nodeEnd = node.end;
		IntervalGraph.Node tmpNode;
		
		if (startOrEnd) {
			for (i=nodeIndex-1; i>=0; i--) {
				tmpNode=IntervalGraph.nodesArray[i];
				if (tmpNode.read!=nodeRead || tmpNode.start<nodeStart-IDENTITY_THRESHOLD) break;
				if (tmpNode.lastPathWithStart==-1) continue;
				System.arraycopy(tmpNode.pathsWithStart,0,tmpArray,0,tmpNode.lastPathWithStart+1);
				last1=Math.makePositive(tmpArray,0,tmpNode.lastPathWithStart);
				if (Math.nonemptyIntersection(tmpArray,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel)) return true;
			}
			for (i=nodeIndex+1; i<nNodes; i++) {
				tmpNode=IntervalGraph.nodesArray[i];
				if (tmpNode.read!=nodeRead || tmpNode.start>nodeStart+IDENTITY_THRESHOLD) break;
				if (tmpNode.lastPathWithStart==-1) continue;
				System.arraycopy(tmpNode.pathsWithStart,0,tmpArray,0,tmpNode.lastPathWithStart+1);
				last1=Math.makePositive(tmpArray,0,tmpNode.lastPathWithStart);
				if (Math.nonemptyIntersection(tmpArray,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel)) return true;
			}
		}
		else {
			for (i=nodeIndex-1; i>=0; i--) {
				tmpNode=IntervalGraph.nodesArray[i];
				if (tmpNode.read!=nodeRead) break;
				if (Math.abs(tmpNode.end,nodeEnd)>IDENTITY_THRESHOLD || tmpNode.lastPathWithEnd==-1) continue;
				System.arraycopy(tmpNode.pathsWithEnd,0,tmpArray,0,tmpNode.lastPathWithEnd+1);
				last1=Math.makePositive(tmpArray,0,tmpNode.lastPathWithEnd);
				if (Math.nonemptyIntersection(tmpArray,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel)) return true;
			}
			for (i=nodeIndex+1; i<nNodes; i++) {
				tmpNode=IntervalGraph.nodesArray[i];
				if (tmpNode.read!=nodeRead || tmpNode.start>=nodeEnd) break;
				if (Math.abs(tmpNode.end,nodeEnd)>IDENTITY_THRESHOLD || tmpNode.lastPathWithEnd==-1) continue;
				System.arraycopy(tmpNode.pathsWithEnd,0,tmpArray,0,tmpNode.lastPathWithEnd+1);
				last1=Math.makePositive(tmpArray,0,tmpNode.lastPathWithEnd);
				if (Math.nonemptyIntersection(tmpArray,0,last1,descendants[kernel],0,lastDescendant[kernel],kernel)) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Remark: this is used for removing $position$ if it is an end of kernel $k$ that is 
	 * also contained in kernel $k$. This might make no occurrence of $k$ count, in a
	 * maximal set of straddling/contained intervals tagged with $k$.
	 *
	 * @param nodeIndex position >=0 in $IntervalGraph.nodesArray$ of a node in $read$ 
	 * whose start is close to $position$;
	 * @return TRUE iff an interval graph node tagged with $kernel$ strictly contains 
	 * $position$ in $read$.
	 */
	private static final boolean getFrequency_inKernel(int position, int read, int kernel, int nodeIndex) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i, j;
		int last;
		final int nNodes = IntervalGraph.nNodes;
		IntervalGraph.Node tmpNode;

		j=nodeIndex; tmpNode=IntervalGraph.nodesArray[j];
		while (j>0 && tmpNode.read==read && tmpNode.start>=position-IDENTITY_THRESHOLD) {
			j--; tmpNode=IntervalGraph.nodesArray[j];
		}
		for (i=j; i>=0; i--) {
			tmpNode=IntervalGraph.nodesArray[i];
			if (tmpNode.read!=read) break;
			last=tmpNode.lastKernel;
			if (last==-1) continue;
			if (tmpNode.end>position+IDENTITY_THRESHOLD && Arrays.binarySearch(tmpNode.kernels,0,last+1,kernel)>=0) return true;
		}
		for (i=j+1; i<nNodes; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			if (tmpNode.read!=read || tmpNode.start>=position-IDENTITY_THRESHOLD) break;
			last=tmpNode.lastKernel;
			if (last==-1) continue;
			if (tmpNode.end>position+IDENTITY_THRESHOLD && Arrays.binarySearch(tmpNode.kernels,0,last+1,kernel)>=0) return true;
		}
		return false;
	}
	
	
	/**
	 * Stores in the global variable $descendants[X][0..lastDescendant[X]]$ the sorted 
	 * list of kernels $Y$ such that $Y$ descends from $X$ in the kernel graph.
	 *
	 * Remark: the procedure assumes that all arcs in the kernel graph are such that all 
	 * nodes of the contained kernel are contained in the container kernel.
	 *
	 * Remark: the procedure assumes all elements of $inNeighbors$ and $outNeighbors$ to 
	 * be non-negative, and it assumes every $inNeighbors[i]$ to be sorted.
	 *
	 * @param includePeriodic FALSE=does not include cyclic kernels or nodes marked in 
	 * $pathKernelPeriodic$, and does not compute the descendants list of such nodes;
	 * @param stack,pushed temporary space, of size at least equal to the number of path
	 * kernels.
	 */
	private static final void getFrequency_buildDescendantsList(boolean includePeriodic, int[] stack, int[] pushed) {
		final int GROWTH_RATE = 8;  // Arbitrary
		int i, j, k;
		int neighbor, kernel, top, tmpInt;
		IntervalGraph.Node tmpNode;
		int[] tmpArray;
		
		if (descendants==null) descendants = new int[nPathKernels][0];
		if (lastDescendant==null) lastDescendant = new int[nPathKernels];
		Math.set(lastDescendant,nPathKernels-1,-1);
		Math.set(pushed,nPathKernels-1,-1);
		for (i=0; i<nPathKernels; i++) {
			if (!includePeriodic && (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0))) continue;
			if (descendants[i]==null || descendants[i].length<lastOutNeighbor[i]+1) descendants[i] = new int[lastOutNeighbor[i]+1];
			top=0; stack[top]=i; pushed[i]=i;
			while (top>=0) {
				kernel=stack[top--];
				for (j=0; j<=lastOutNeighbor[kernel]; j++) {
					neighbor=outNeighbors[kernel][j];
					if (!includePeriodic && (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0))) continue;
					if (Arrays.binarySearch(descendants[i],0,lastDescendant[i]+1,neighbor)<0) {
						lastDescendant[i]++;
						if (lastDescendant[i]==descendants[i].length) {
							tmpArray = new int[lastDescendant[i]+GROWTH_RATE];
							System.arraycopy(descendants[i],0,tmpArray,0,lastDescendant[i]);
							descendants[i]=tmpArray;
						}
						k=lastDescendant[i];
						descendants[i][k]=neighbor;
						while (k>0 && descendants[i][k]<descendants[i][k-1]) {
							tmpInt=descendants[i][k-1];
							descendants[i][k-1]=descendants[i][k];
							descendants[i][k]=tmpInt;
							k--;
						}
					}
					if (pushed[neighbor]!=i) {
						stack[++top]=neighbor;
						pushed[neighbor]=i;
					}
				}
			}
		}
	}
	
	
	/**
	 * Stores in $tmpArray3[0..X]$ the result of $node.pathsWithStart \intersect 
	 * node.pathsWithEnd$, where $X$ is returned in output, and intersection works on
	 * IDs after they have been internally transformed to non-negative. The procedure 
	 * marks in $tmpArray2[i]$ the ends of $tmpArray3[i]$ that appear at each end of 
	 * $node$: 0=different kernel ends; 1=same ends, beginning of the kernel path; 2=same 
	 * ends, end of the kernel path.
	 *
	 * @param tmpArray* temporary space, of size at least $node.lastKernel+1$.
	 */
	private static final int startIntersectEnd(IntervalGraph.Node node, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) {
		int i, j, k;
		int tmpInt, previous, last, last1, last2, last3;
		
		// Trivial cases
		if (node.lastPathWithStart==-1 || node.lastPathWithEnd==-1) return -1;
		if (node.lastPathWithStart==0 && node.lastPathWithEnd==0) {
			i=node.pathsWithStart[0];
			if (i<0) i=-1-i;
			j=node.pathsWithEnd[0];
			if (j<0) j=-1-j;
			if (i!=j) return -1;
			tmpArray3[0]=i;
			if (node.pathsWithStart[0]>=0 && node.pathsWithEnd[0]>=0) tmpArray2[0]=1;
			else if (node.pathsWithStart[0]<0 && node.pathsWithEnd[0]<0) tmpArray2[0]=2;
			else tmpArray2[0]=0;
			return 0;
		}
		
		// Finding the kernels that occur both at the start and at the end of $node$.
		last=node.lastPathWithStart;
		System.arraycopy(node.pathsWithStart,0,tmpArray1,0,last+1);
		last1=Math.makePositive(tmpArray1,0,last);
		last=node.lastPathWithEnd;
		System.arraycopy(node.pathsWithEnd,0,tmpArray2,0,last+1);
		last2=Math.makePositive(tmpArray2,0,last);
		last3=Math.setIntersection(tmpArray1,0,last1,tmpArray2,0,last2,tmpArray3,0);
		if (last3==-1) return -1;
		
		// Computing the ends of each kernel
		Math.set(tmpArray2,last3,-1);
		k=-1;
		for (i=node.lastPathWithStart; i>=0; i--) {
			if (node.pathsWithStart[i]<0) {
				k=i;
				break;
			}
		}
		if (k>=0) {
			i=k; j=0;
			while (i>=0 && j<=last3) {
				tmpInt=-1-node.pathsWithStart[i];
				if (tmpInt<tmpArray3[j]) {
					i--;
					continue;
				}
				if (tmpArray3[j]<tmpInt) {
					j++;
					continue;
				}
				tmpArray2[j]=2;
				i--; j++;
			}
		}
		if (k+1<=node.lastPathWithStart) {
			i=k+1; j=0;
			while (i<=node.lastPathWithStart && j<=last3) {
				tmpInt=node.pathsWithStart[i];
				if (tmpInt<tmpArray3[j]) {
					i++;
					continue;
				}
				if (tmpArray3[j]<tmpInt) {
					j++;
					continue;
				}
				if (tmpArray2[j]==-1) tmpArray2[j]=1;
				else if (tmpArray2[j]==0 || tmpArray2[j]==1) { /* NOP */ } 
				else tmpArray2[j]=0;
				i++; j++;
			}
		}
		k=-1;
		for (i=node.lastPathWithEnd; i>=0; i--) {
			if (node.pathsWithEnd[i]<0) {
				k=i;
				break;
			}
		}
		if (k>=0) {
			i=k; j=0;
			while (i>=0 && j<=last3) {
				tmpInt=-1-node.pathsWithEnd[i];
				if (tmpInt<tmpArray3[j]) {
					i--;
					continue;
				}
				if (tmpArray3[j]<tmpInt) {
					j++;
					continue;
				}
				if (tmpArray2[j]==-1) tmpArray2[j]=2;
				else if (tmpArray2[j]==0 || tmpArray2[j]==2) { /* NOP */ } 
			 	else tmpArray2[j]=0;
				i--; j++;
			}
		}
		if (k+1<=node.lastPathWithEnd) {
			i=k+1; j=0;
			while (i<=node.lastPathWithEnd && j<=last3) {
				tmpInt=node.pathsWithEnd[i];
				if (tmpInt<tmpArray3[j]) {
					i++;
					continue;
				}
				if (tmpArray3[j]<tmpInt) {
					j++;
					continue;
				}
				if (tmpArray2[j]==-1) tmpArray2[j]=1;
				else if (tmpArray2[j]==0 || tmpArray2[j]==1) { /* NOP */ } 
				else tmpArray2[j]=0;
				i++; j++;
			}
		}
		return last3;
	}
	
	
	/**
	 * Stores in $tmpArray3[0..X]$ the result of $node.pathsWithStart \cup 
	 * node.pathsWithEnd$, where $X$ is returned in output, and union works on IDs after 
	 * they have been internally transformed to non-negative.
	 *
	 * @param tmpArray* temporary space, of size at least $node.lastKernel+1$.
	 */
	private static final int startUnionEnd(IntervalGraph.Node node, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) {
		int i;
		int last, last1, last2;
		
		// Trivial cases
		if (node.lastPathWithStart==-1 && node.lastPathWithEnd==-1) return -1;
		else if (node.lastPathWithStart==0 && node.lastPathWithEnd==-1) {
			i=node.pathsWithStart[0];
			if (i<0) i=-1-i;
			tmpArray3[0]=i;
			return 0;
		}
		else if (node.lastPathWithStart==-1 && node.lastPathWithEnd==0) {
			i=node.pathsWithEnd[0];
			if (i<0) i=-1-i;
			tmpArray3[0]=i;
			return 0;
		}
		
		// Computing the union
		last=node.lastPathWithStart;
		if (last>=0) {
			System.arraycopy(node.pathsWithStart,0,tmpArray1,0,last+1);
			last1=Math.makePositive(tmpArray1,0,last);
		}
		else last1=-1;
		last=node.lastPathWithEnd;
		if (last>=0) {
			System.arraycopy(node.pathsWithEnd,0,tmpArray2,0,last+1);
			last2=Math.makePositive(tmpArray2,0,last);
		}
		else last2=-1;
		return Math.setUnion(tmpArray1,last1,tmpArray2,last2,tmpArray3);
	}
	
	
	/**
	 * Stores in $tmpArray1[0..X]$ the set $node.kernels \setminus pathsWithStart$, where 
	 * $X$ is returned in output and $pathsWithStart$ is internally transformed to 
	 * non-negative.
	 *
	 * @param tmpArray2 temporary space.
	 */
	private static final int notStart(IntervalGraph.Node node, int[] tmpArray1, int[] tmpArray2) {
		int i, j;
		int previous, last, last1, last2, tmpInt;
		
		last1=-1; last=node.lastPathWithStart;
		if (last>=0) {
			System.arraycopy(node.pathsWithStart,0,tmpArray2,0,last+1);
			last2=Math.makePositive(tmpArray2,0,last);
			last1=Math.setMinus(node.kernels,node.lastKernel,tmpArray2,last2,tmpArray1);
		}
		else {
			last1=node.lastKernel;
			System.arraycopy(node.kernels,0,tmpArray1,0,last1+1);
		}
		return last1;
	}
	
	
	/**
	 * Stores in $tmpArray1[0..X]$ the set $node.kernels \setminus pathsWithEnd$, where 
	 * $X$ is returned in output and $pathsWithEnd$ is internally transformed to 
	 * non-negative.
	 *
	 * @param tmpArray2 temporary space.
	 */
	private static final int notEnd(IntervalGraph.Node node, int[] tmpArray1, int[] tmpArray2) {
		int i, j;
		int previous, last, last1, last2, tmpInt;
		
		last1=-1; last=node.lastPathWithEnd;
		if (last>=0) {
			System.arraycopy(node.pathsWithEnd,0,tmpArray2,0,last+1);
			last2=Math.makePositive(tmpArray2,0,last);
			last1=Math.setMinus(node.kernels,node.lastKernel,tmpArray2,last2,tmpArray1);
		}
		else {
			last1=node.lastKernel;
			System.arraycopy(node.kernels,0,tmpArray1,0,last1+1);
		}
		return last1;
	}
	
	
	/**
	 * Stores in $tmpArray2[0..X]$ the result of $node.kernels \setminus 
	 * node.pathsWithStart \setminus node.pathsWithEnd$, where $X$ is returned in output 
	 * and each $pathsWith*$ is internally transformed to non-negative.
	 *
	 * @param tmpArray* temporary space, of size at least $node.lastKernel+1$.
	 */
	private static final int neitherStartNorEnd(IntervalGraph.Node node, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) {
		int i, j;
		int previous, last, last1, last2, last3, tmpInt;
		
		// Subtracting $pathsWithStart$.
		last1=notStart(node,tmpArray1,tmpArray3);
		
		// Subtracting $pathsWithEnd$.
		last=node.lastPathWithEnd; last2=-1;
		if (last>=0) {
			System.arraycopy(node.pathsWithEnd,0,tmpArray3,0,last+1);
			last3=Math.makePositive(tmpArray3,0,last);
			last2=Math.setMinus(tmpArray1,last1,tmpArray3,last3,tmpArray2);
		}
		else {
			last2=last1;
			System.arraycopy(tmpArray1,0,tmpArray2,0,last1+1);
		}
		
		return last2;
	}
	
	
	/**
	 * Remark: the procedure assumes $tmpPoints$ to be of size at least $nPathKernels$.
	 * Remark: the procedure reads just $kernelSize$ and $kernelFrequency$, and it does 
	 * not scan interval graph nodes.
	 *
	 * @param out output array of $basinStats()$.
	 */
	private static final void printBasinHistograms(int[] out) {
		int i;
		
		System.out.println("Path kernels: "+nPathKernels+". Nodes: "+IntervalGraph.nNodes+". Nodes with one path kernel: "+out[1]+". Nodes with no path kernel: "+out[2]);
		
		// Kernel size
		lastTmpPoint=-1;
		for (i=0; i<nPathKernels; i++) {
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=kernelSize[i];
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		System.out.println("Histogram of path kernel basin sizes:");
		for (i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i].position+","+tmpPoints[i].mass);
		
		// Kernel frequency
		lastTmpPoint=-1;
		for (i=0; i<nPathKernels; i++) {
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=kernelFrequency[i][0]+(kernelFrequency[i][1]+kernelFrequency[i][2])>>1;
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		System.out.println("Histogram of estimated path kernel frequencies:");
		for (i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i].position+","+tmpPoints[i].mass);
	}

	
	/**
	 * Given an assignment of kernels to nodes of the interval graph, the procedure tries 
	 * to verify the kernels of all intervals that are not left-maximal or not 
	 * right-maximal, using the reads in which they occur. The procedure looks to the left
	 * and to the right of every interval $I$ inside its read, and keeps a kernel-path tag
	 * $K$ iff the same kernel occurs in a nearby interval $J$ whose distance from $I$ is 
	 * within the estimated string length of path $K$ (assuming long low-quality regions 
	 * to be substitutions).
	 *
	 * Remark: the procedure does not change the kernel tags of a node if no tag occurs 
	 * within suitable distance.
	 *
	 * Remark: the procedure might make a node have a kernel tag that none of its 
	 * ancestors by containment has. Fixing this inconsistency takes too much effort:* we 
	 * don't do it for simplicity. However, it might be the reason why a node is reported 
	 * by $FragmentsStep1$.
	 *
	 * Remark: the procedure assumes that $IntervalGraph.nodesArray$ is sorted by read, 
	 * and the $kernels$ array of every node to be sorted in ascending order.
	 *
	 * @param kernelsSorted if false, the procedure sorts the $kernels$ array of every 
	 * node in ascending order;
	 * @return the number of nodes that have been corrected by the procedure.
	 */
	private static final int checkKernelsWithReads(int maxAlignmentsPerRead, int maxKernelsPerInterval, boolean kernelsSorted) {
		int i, j;
		int id, readA, previousReadA, lastNodeInBatch, tmpOrder, lastCorrectedNode;
		IntervalGraph.Node node;
		int[] notFoundLeft = new int[maxKernelsPerInterval];
		int[] notFoundRight = new int[maxKernelsPerInterval];
		if (tmpArray1==null || tmpArray1.length<maxKernelsPerInterval<<1) tmpArray1 = new int[maxKernelsPerInterval<<1];
		if (tmpArray2==null || tmpArray2.length<maxKernelsPerInterval<<1) tmpArray2 = new int[maxKernelsPerInterval<<1];
		if (tmpArray3==null || tmpArray3.length<maxKernelsPerInterval<<1) tmpArray3 = new int[maxKernelsPerInterval<<1];
		if (tmpArray4==null || tmpArray4.length<nNodesInKernel) tmpArray4 = new int[nNodesInKernel];
		System.arraycopy(nodesInKernel,0,tmpArray4,0,nNodesInKernel);
		if (nNodesInKernel>1) Arrays.sort(tmpArray4,0,nNodesInKernel);
		IntervalGraph.Node[] batch = new IntervalGraph.Node[maxAlignmentsPerRead];
		
		if (!kernelsSorted) {
			for (i=0; i<IntervalGraph.nNodes; i++) {
				if (IntervalGraph.nodesArray[i].lastKernel>0) Arrays.sort(IntervalGraph.nodesArray[i].kernels,0,IntervalGraph.nodesArray[i].lastKernel+1);
			}
		}
		
		// Loading in $batch$ the set of all distinct intervals (excluding periodic
		// intervals) in the same read.
		previousReadA=-1; lastNodeInBatch=-1; lastCorrectedNode=-1; j=0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			readA=IntervalGraph.nodesArray[i].read;
			if (readA!=previousReadA) {
				if (lastNodeInBatch>0) {
					tmpOrder=IntervalGraph.Node.order;
					IntervalGraph.Node.order=IntervalGraph.Node.START;
					if (lastNodeInBatch>0) Arrays.sort(batch,0,lastNodeInBatch+1);
					IntervalGraph.Node.order=tmpOrder;
					lastCorrectedNode=checkKernelsWithReads_impl(batch,lastNodeInBatch,notFoundLeft,notFoundRight,lastCorrectedNode+1,tmpArray1,tmpArray2,tmpArray3);					
				}
				previousReadA=readA;
				lastNodeInBatch=-1;
			}
			while (j<nNodesInKernel && tmpArray4[j]<i) j++;		
			if ((j==nNodesInKernel || tmpArray4[j]!=i) && (IntervalGraph.nodesArray[i].type!=Constants.INTERVAL_PERIODIC || IntervalGraph.nodesArray[i].hasLongPeriod)) batch[++lastNodeInBatch]=IntervalGraph.nodesArray[i];
		}
		if (lastNodeInBatch>0) {
			tmpOrder=IntervalGraph.Node.order;
			IntervalGraph.Node.order=IntervalGraph.Node.START;
			if (lastNodeInBatch>0) Arrays.sort(batch,0,lastNodeInBatch+1);
			IntervalGraph.Node.order=tmpOrder;
			lastCorrectedNode=checkKernelsWithReads_impl(batch,lastNodeInBatch,notFoundLeft,notFoundRight,lastCorrectedNode+1,tmpArray1,tmpArray2,tmpArray3);
		}
		for (i=0; i<batch.length; i++) batch[i]=null;
		batch=null;
		
		return lastCorrectedNode+1;
	}
	
	
	/**
	 * Remark: the procedure assumes the $kernels$ array of every node to be sorted in
	 * ascending order.
	 *
	 * @param batch set of distinct node pointers that do not point to periodic substring
	 * intervals, sorted by node start position;
	 * @param notFoundLeft,notFoundRight,tmp1,tmp2,tmp3 temporary space;
	 * @param correctedNodes nodes that have been corrected by the procedure;
	 * @return the procedure stores in global variable $tmpArray4[from...]$ the nodes it 
	 * corrected, and it returns the index of the last such node.
	 */
	private static final int checkKernelsWithReads_impl(IntervalGraph.Node[] batch, int lastNodeInBatch, int[] notFoundLeft, int[] notFoundRight, int from, int[] tmp1, int[] tmp2, int[] tmp3) {
		final int READ = batch[0].read;
		final int READ_LENGTH = Reads.getReadLength(READ);
		final int THRESHOLD = IO.quantum;  // Arbitrary
		int i, j, k, i1, i2;
		int length, nUnreachable, lastNotFoundLeft, lastNotFoundRight, lastNotFound, lastCorrectedNode;
		IntervalGraph.Node node;
		int[] swap, notFound;
		
		lastCorrectedNode=from-1;
		for (i=0; i<=lastNodeInBatch; i++) {
			node=batch[i];
			if (node.lastKernel==0 || (node.isLeftMaximal && node.isRightMaximal)) continue;
			
			// Checking on the left side
			if (!node.isLeftMaximal && node.start>=Alignments.minAlignmentLength) {
				lastNotFoundLeft=notStart(node,notFoundLeft,tmp3);
				for (j=i-1; j>=0; j--) {
					if (batch[j].isRightMaximal || batch[j].end>=node.start) continue;
					length=node.end-batch[j].start+1;
					// Deciding whether there is hope to remove elements from $notFound$.
					nUnreachable=0;
					for (k=0; k<=lastNotFoundLeft; k++) {
						if (length>pathKernelLengths[notFoundLeft[k]]+THRESHOLD) nUnreachable++;
					}
					if (nUnreachable==lastNotFoundLeft+1) break;
					// Updating $notFoundLeft$.
					i1=0; i2=0; k=-1;
					while (i1<=lastNotFoundLeft && i2<=batch[j].lastKernel) {
						if (notFoundLeft[i1]<batch[j].kernels[i2]) tmp1[++k]=notFoundLeft[i1++];
						else if (notFoundLeft[i1]>batch[j].kernels[i2]) i2++;
						else {
							if (length>pathKernelLengths[notFoundLeft[i1]]+THRESHOLD) tmp1[++k]=notFoundLeft[i1];
							i1++; i2++;
						}
					}
					while (i1<=lastNotFoundLeft) tmp1[++k]=notFoundLeft[i1++];
					lastNotFoundLeft=k;
					swap=notFoundLeft; notFoundLeft=tmp1; tmp1=swap;
				}
			}
			else lastNotFoundLeft=-1;
			
			// Checking on the right side
			if (!node.isRightMaximal && node.end<READ_LENGTH-Alignments.minAlignmentLength) {
				lastNotFoundRight=notEnd(node,notFoundRight,tmp3);
				for (j=i+1; j<=lastNodeInBatch; j++) {
					if (batch[j].isLeftMaximal || batch[j].start<=node.end) continue;
					length=batch[j].end-node.start+1;
					// Deciding whether there is hope to remove elements from $notFound$.
					nUnreachable=0;
					for (k=0; k<=lastNotFoundRight; k++) {
						if (batch[j].start-node.start>pathKernelLengths[notFoundRight[k]]) nUnreachable++;
					}
					if (nUnreachable==lastNotFoundRight+1) break;
					// Updating $notFoundRight$.
					i1=0; i2=0; k=-1;
					while (i1<=lastNotFoundRight && i2<=batch[j].lastKernel) {
						if (notFoundRight[i1]<batch[j].kernels[i2]) tmp1[++k]=notFoundRight[i1++];
						else if (notFoundRight[i1]>batch[j].kernels[i2]) i2++;
						else {
							if (length>pathKernelLengths[notFoundRight[i1]]+THRESHOLD) tmp1[++k]=notFoundRight[i1];
							i1++; i2++;
						}
					}
					while (i1<=lastNotFoundRight) tmp1[++k]=notFoundRight[i1++];
					lastNotFoundRight=k;
					swap=notFoundRight; notFoundRight=tmp1; tmp1=swap;
				}
			}
			else lastNotFoundRight=-1;
			
			// Subtracting the intersection of $notFoundLeft$ and $notFoundRight$ from
			// $node.kernels$.
			if (lastNotFoundLeft==-1) {
				if (lastNotFoundRight==-1) continue;
				notFound=notFoundRight;
				lastNotFound=lastNotFoundRight;
			}
			else {
				if (lastNotFoundRight==-1) {
					notFound=notFoundLeft;
					lastNotFound=lastNotFoundLeft;
				}
				else {
					lastNotFound=Math.setIntersection(notFoundLeft,0,lastNotFoundLeft,notFoundRight,0,lastNotFoundRight,tmp1,0);
					notFound=tmp1;
				}
			}
			if (lastNotFound==node.lastKernel) {
				// Keeping all kernels if none of them was found to the left or right
				continue;
			}
			// Removing cyclic kernels from $notFound$, if any.
			i2=-1;
			for (i1=0; i1<=lastNotFound; i1++) {
				k=notFound[i1];
				if (pathKernelLengths[k]>0) notFound[++i2]=k;
			}
			lastNotFound=i2;
			// Subtracting
			i1=0; i2=0; k=-1;
			while (i1<=node.lastKernel && i2<=lastNotFound) {
				if (node.kernels[i1]<notFound[i2]) {
					k++;
					tmp2[k]=node.kernels[i1];
					tmp3[k]=node.kernelOrientations[i1];
					i1++;
				}
				else if (node.kernels[i1]>notFound[i2]) i2++;
				else { i1++; i2++; }
			}
			while (i1<=node.lastKernel) {
				k++;
				tmp2[k]=node.kernels[i1];
				tmp3[k]=node.kernelOrientations[i1];
				i1++;
			}
			node.lastKernel=k;
			System.arraycopy(tmp2,0,node.kernels,0,k+1);
			for (j=0; j<=k; j++) node.kernelOrientations[j]=(byte)tmp3[j];
			// It's correct to update the $kernels$ array inside the loop (instead of
			// updating all arrays after the loop).
			//
			// Updating $pathsWith*$.
			k=-1;
			for (i1=0; i1<=node.lastPathWithStart; i1++) {
				i2=node.pathsWithStart[i1];
				if (i2<0) i2=-1-i2;
				if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,i2)>=0) node.pathsWithStart[++k]=node.pathsWithStart[i1];
			}
			node.lastPathWithStart=k;
			k=-1;
			for (i1=0; i1<=node.lastPathWithEnd; i1++) {
				i2=node.pathsWithEnd[i1];
				if (i2<0) i2=-1-i2;
				if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,i2)>=0) node.pathsWithEnd[++k]=node.pathsWithEnd[i1];
			}
			node.lastPathWithEnd=k;
			lastCorrectedNode++;
			if (lastCorrectedNode==tmpArray4.length) {
				int[] newArray = new int[tmpArray4.length<<1];
				System.arraycopy(tmpArray4,0,newArray,0,tmpArray4.length);
				tmpArray4=newArray;
			}
			tmpArray4[lastCorrectedNode]=node.nodeID;
		}
		return lastCorrectedNode;
	}

	
	
	
	
	
	
	
	// ---------------------------------- KERNEL GRAPH -----------------------------------
	
	/**
	 * Builds a bitvector $S$ such that $S[i]=1$ iff the $i$-th alignment in 
	 * $alignmentsPath$ is contained in the union of a set of overlapping short-period 
	 * intervals in either readA or readB.
	 *
	 * @param shortPeriodPath file containing a similar bitvector that marks alignments
	 * contained in short-period intervals in readA only.
	 */
	public static final boolean[] loadShortPeriodBitvector(int nAlignments, String alignmentsPath, String shortPeriodPath) throws IOException {
		final int HASHSET_CAPACITY = nAlignments>>2;  // Arbitrary
		int i;
		String str1, str2;
		AlignmentPair newPair;
		AlignmentPair tmpPair = new AlignmentPair();
		BufferedReader alignmentsFile, shortPeriodFile;
		HashSet<AlignmentPair> shortPeriodAlignments = new HashSet<AlignmentPair>(HASHSET_CAPACITY,1);
		boolean[] out = new boolean[nAlignments];
	
		// Marking alignments contained in a short-period interval in readA
		alignmentsFile = new BufferedReader(new FileReader(alignmentsPath),IO.BUFFER_SIZE);
		shortPeriodFile = new BufferedReader(new FileReader(shortPeriodPath),IO.BUFFER_SIZE);
		alignmentsFile.readLine();
		alignmentsFile.readLine();  // Skipping the first two lines
		str1=alignmentsFile.readLine(); str2=shortPeriodFile.readLine();
		i=-1;
		while (str1!=null) {
			i++;
			if (str2.charAt(0)=='1') {
				out[i]=true;
				Alignments.readAlignmentFile(str1);
				tmpPair.set(Alignments.readA-1,Alignments.startA,Alignments.endA,Alignments.readB-1,Alignments.startB,Alignments.endB,Alignments.orientation);
				if (!shortPeriodAlignments.contains(tmpPair)) {
					newPair = new AlignmentPair();
					tmpPair.clone(newPair);
					shortPeriodAlignments.add(newPair);
				}
			}
			else out[i]=false;
			str1=alignmentsFile.readLine(); str2=shortPeriodFile.readLine();
		}
		alignmentsFile.close(); alignmentsFile=null;
		shortPeriodFile.close(); shortPeriodFile=null;
	
		// Marking alignments contained in a short-period interval in readB
		alignmentsFile = new BufferedReader(new FileReader(alignmentsPath),IO.BUFFER_SIZE);
		alignmentsFile.readLine();
		alignmentsFile.readLine();  // Skipping the first two lines
		str1=alignmentsFile.readLine();
		i=-1;
		while (str1!=null) {
			i++;
			if (!out[i]) {
				Alignments.readAlignmentFile(str1);
				tmpPair.set(Alignments.readA-1,Alignments.startA,Alignments.endA,Alignments.readB-1,Alignments.startB,Alignments.endB,Alignments.orientation);
				if (shortPeriodAlignments.contains(tmpPair)) out[i]=true;
			}
			str1=alignmentsFile.readLine();
		}
		alignmentsFile.close(); alignmentsFile=null;
		shortPeriodAlignments.clear(); shortPeriodAlignments=null;
	
		return out;
	}
	
	
	/**
	 * All procedures either impose or assume that objects are in canonical form.
	 */
	private static class AlignmentPair implements Comparable {
		public int read1, start1, end1;
		public int read2, start2, end2;
		public boolean orientation;
		
		public AlignmentPair() { 
			orientation=true;
			read1=-1; start1=-1; end1=-1;
			read2=-1; start2=-1; end2=-1;
		}
		
		public AlignmentPair(int r1, int s1, int e1, int r2, int s2, int e2, boolean o) {
			set(r1,s1,e1,r2,s2,e2,o);
		}
		
		public void set(int r1, int s1, int e1, int r2, int s2, int e2, boolean o) {
			orientation=o;
			if (r1<r2) {
				read1=r1; start1=s1; end1=e1;
				read2=r2; start2=s2; end2=e2;
			}
			else {
				read1=r2; start1=s2; end1=e2;
				read2=r1; start2=s1; end2=e1;
			}
		}
		
		/**
		 * Sets $to$ to be a copy of the current object.
		 */
		public void clone(AlignmentPair to) {
			to.read1=read1; to.start1=start1; to.end1=end1;
			to.read2=read2; to.start2=start2; to.end2=end2;
			to.orientation=orientation;
		}
		
		public int compareTo(Object other) {
			AlignmentPair otherPair = (AlignmentPair)other;
			if (read1<otherPair.read1) return -1;
			else if (read1>otherPair.read1) return 1;
			if (read2<otherPair.read2) return -1;
			else if (read2>otherPair.read2) return 1;
			if (start1<otherPair.start1) return -1;
			else if (start1>otherPair.start1) return 1;
			if (start2<otherPair.start2) return -1;
			else if (start2>otherPair.start2) return 1;
			return 0;
		}
		
		public boolean equals(Object other) {
			AlignmentPair otherPair = (AlignmentPair)other;
			return orientation==otherPair.orientation &&
				   read1==otherPair.read1 && start1==otherPair.start1 && end1==otherPair.end1 &&
				   read2==otherPair.read2 && start2==otherPair.start2 && end2==otherPair.end2;
		}
		
		public int hashCode() {
			return (read1+"-"+start1+"-"+end1+"-"+read2+"-"+start2+"-"+end2+"-"+(orientation?"1":"0")).hashCode();
		}
	}
	
	
	/**
	 * Transforms $nodesInKernel$ into a possibly larger array, in which all nodes marked 
	 * with the same kernel tag are located in a compact block, and in which a node occurs
	 * in the blocks of all kernels to which it belongs. The procedure writes in global 
	 * array $blockKernels$ the kernel ID of every element in the new $nodesInKernel$, and
	 * sets $nNodesInKernel$ to the size of the new $nodesInKernel$. Every block is sorted
	 * by $bidirectedGraphNode$.
	 *
	 * Remark: the procedure discards nodes with no $bidirectedGraphNode$ pointer.
	 */
	private static final void expandNodesInKernel() {
		int i, j, k;
		int count, last, kernel, previousOrder, previousFirst, previousKernel, maxLength;
		IntervalGraph.Node tmpNode;
		int[] blockStarts, newNodesInKernel;
		if (nNodesInKernel==0 || nPathKernels==0) return;		
		
		// Computing blocks
		blockStarts = new int[nPathKernels];
		Math.set(blockStarts,nPathKernels-1,0);
		count=0;
		for (i=0; i<nNodesInKernel; i++) {
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[i]];
			if (tmpNode.bidirectedGraphNode==-1) continue;
			last=tmpNode.lastKernel;
			count+=last+1;
			for (j=0; j<=last; j++) blockStarts[tmpNode.kernels[j]]++;
		}
		newNodesInKernel = new int[count];
		j=blockStarts[0]; blockStarts[0]=0;
		for (i=1; i<nPathKernels; i++) {
			k=blockStarts[i];
			blockStarts[i]=j;
			j+=k;
		}
		blockKernels = new int[count];
		for (i=1; i<nPathKernels; i++) {
			for (j=blockStarts[i-1]; j<blockStarts[i]; j++) blockKernels[j]=i-1;
		}
		for (j=blockStarts[nPathKernels-1]; j<count; j++) blockKernels[j]=nPathKernels-1;
		for (i=0; i<nNodesInKernel; i++) {
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[i]];
			if (tmpNode.bidirectedGraphNode==-1) continue;
			last=tmpNode.lastKernel;
			for (j=0; j<=last; j++) {
				kernel=tmpNode.kernels[j];
				newNodesInKernel[blockStarts[kernel]++]=nodesInKernel[i];
			}
		}
		nodesInKernel=newNodesInKernel;
		nNodesInKernel=count;
		
		// Sorting blocks
		previousOrder=IntervalGraph.Node.order;
		IntervalGraph.Node.order=IntervalGraph.Node.BIDIRECTED;
		previousFirst=0; previousKernel=blockKernels[0];
		for (i=1; i<nNodesInKernel; i++) {
			if (blockKernels[i]==previousKernel) continue;
			if (i>previousFirst+1) Arrays.sort(nodesInKernel,previousFirst,i);
			previousFirst=i; previousKernel=blockKernels[i];
		}
		if (nNodesInKernel>previousFirst+1) Arrays.sort(nodesInKernel,previousFirst,nNodesInKernel);
		IntervalGraph.Node.order=previousOrder;
	}
	
	
	public static final void buildPathKernelLongestNode() {
		int i, j;
		int length, kernel, lastKernel;
		IntervalGraph.Node tmpNode;
		
		pathKernelLongestNode = new int[nPathKernels];
		Math.set(pathKernelLongestNode,nPathKernels-1,0);
		for (i=0; i<nNodesInKernel; i++) {
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[i]];
			length=tmpNode.length();
			lastKernel=tmpNode.lastKernel;
			for (j=0; j<=lastKernel; j++) {
				kernel=tmpNode.kernels[j];
				pathKernelLongestNode[kernel]=Math.max(pathKernelLongestNode[kernel],length);
			}
		}
	}
	
	
	/**
	 * Builds the DAG in which an arc $(v,w)$ means that some nodes in kernel $w$ are 
	 * contained in, or insert into, nodes in the basin of kernel $v$. Containment/
	 * insertion are decided using all alignments in $alignmentsFilePath$, without taking 
	 * into account node types or the assignment of alignments to nodes. A node in $w$
	 * might be a permutation of a substring of $v$: the procedure does not detect this 
	 * and assumes that no kernel is a permutation of another kernel.
	 *
	 * Remark: alignments that have been assigned to a short-period interval during 
	 * factorization are not used to build the kernel graph. Assume e.g. that kernel 1 has
	 * sequence $U_1 V_1$ and kernel 2 has sequence $U_2 V_2$, where $U_1 > U_2$,
	 * $V_1$ is a proper prefix of $V_2$, $U_i$ have short period, and $|U_2 V_2| >= 
	 * |U_1 V_1|$. Kernel 1 might be considered contained in or identical to kernel 2
	 * because its sequence might be fully covered by alignments with (fregments of) 
	 * kernel 2. Note that the short-period intervals discarded alignments have been 
	 * assigned to might not belong to the current interval graph.
	 *
	 * Remark: rather than building the kernel graph at this step, we could instead just 
	 * report all kernel path strings in output, compute all their pairwise alignments,
	 * and do something equivalent to the kernel graph later. These two approaches could
	 * be combined. However, if a repeat has many fragments in the genome, the number of 
	 * path kernels might be large, and computing al pairwise alignments of such fragments
	 * might be slower than just scanning the list of alignments we already have.
	 *
	 * Remark: to enforce that the graph is acyclic, we add an arc only if the string 
	 * lengths of its adjacent kernels are compatible: see $buildKernelGraph_impl_
	 * filter()$ for details.
	 *
	 * Remark: if a kernel is such that all nodes that map to its bidirected graph nodes 
	 * are of substring type, then it is likely that its fragments already belong to the 
	 * kernel's basin. However, because of differential fragmentation (defined in 
	 * $buildNodes2Kernel()$), this is not guaranteed. Symmetrically, because of 
	 * differential fragmentation the kernel might itself be a fragment of another kernel,
	 * or be identical to it. Thus, even kernels with only substring-type nodes should 
	 * be included in the kernel graph.
	 *
	 * Remark: cyclic kernels (i.e. those with $pathKernelLengths[i]=0$) are included in 
	 * the kernel graph. See $expandPathKernelPeriodic()$ for motivations.
	 *
	 * Remark: the procedure assumes kernel basins, $kernelStats$, $kernelSize$, 
	 * $kernelFrequency$ to have already been computed; $IntervalGraph.nodesArray$ to be 
	 * sorted by $read,start$, and the kernel tags of every node to be sorted.
	 *
	 * @return the number of path kernels marked in $pathKernelPeriodic$.
	 */
	private static final int buildKernelGraph(int nAlignments, String alignmentsFilePath, String shortPeriodBitvectorPath, String shortPeriodTrackPath) throws IOException {
		int i, j;
		int nNodes, out;
		final int previousOrder = IntervalGraph.Node.order;
		// Data structures for the $nodes2kernel$ matrix:
		IntervalGraph.Node[] nodes;
		int[][] nodes2kernel;
		int[] lastNode2kernel;
		// Temporary space:
		IntervalGraph.Node tmpNode;
		
		expandNodesInKernel();
		out=buildPathKernelPeriodic();
		System.err.println("buildKernelGraph> "+out+" initial path kernels marked as periodic");
		buildPathKernelLongestNode();
		
		// Building $nodes2kernel$.
		nodes = new IntervalGraph.Node[nNodesInKernel];
		for (i=0; i<nNodesInKernel; i++) nodes[i]=IntervalGraph.nodesArray[nodesInKernel[i]];
		nNodes=nNodesInKernel;
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		if (nNodes>1) Arrays.sort(nodes,0,nNodes);
		IntervalGraph.Node.order=previousOrder;
		j=0; tmpNode=nodes[0];
		for (i=1; i<nNodes; i++) {  // Removing duplicates
			if (nodes[i]==tmpNode) continue;
			tmpNode=nodes[i];
			j++;
			if (j!=i) nodes[j]=nodes[i];
		}
		nNodes=j+1;
		if (shortPeriodAlignments==null) shortPeriodAlignments=loadShortPeriodBitvector(nAlignments,alignmentsFilePath,shortPeriodBitvectorPath);
		if (shortPeriodWindows==null) {
			shortPeriodWindows=IntervalGraph.loadShortPeriodTrack(shortPeriodTrackPath);
			trackPointers=IntervalGraph.buildTrackPointers(shortPeriodWindows);
		}
		nodes2kernel = new int[nNodes][0];
		lastNode2kernel = new int[nNodes];
		buildNodes2Kernel(nodes,nNodes,alignmentsFilePath,nodes2kernel,lastNode2kernel,nAlignments);
		out=buildKernelGraph_main(nAlignments,alignmentsFilePath,nodes,nNodes,nodes2kernel,lastNode2kernel,out);
		// $shortPeriodAlignments,shortPeriodWindows,trackPointers$ should not be
		// deallocated now, since they might be useful for other procedures downstream.
		return out;
	}


	/**
	 * Remark: the procedure uses the global variables $shortPeriodAlignments,
	 * shortPeriodWindows,trackPointers$, which should have already been initialized.
	 *
	 * @param out number of path kernels marked in $pathKernelPeriodic$ before this 
	 * procedure is called;
	 * @return number of path kernels marked in $pathKernelPeriodic$ after this procedure 
	 * is called.
	 */
	public static final int buildKernelGraph_main(int nAlignments, String alignmentsFilePath, IntervalGraph.Node[] nodes, int nNodes, int[][] nodes2kernel, int[] lastNode2kernel, int out) throws IOException {
		final int IDENTITY_THRESHOLD = 100;
		int i;
		int from, arcs, tmpArcs;
		int[] lastPermutationWindow;
		int[][] permutationWindows;
		
		
int fabio = 22;
if (IO.SHOW_INTERACTIVE) {	
System.err.println("nodes2kernel for kernel path "+fabio);		
for (int x=0; x<nNodes; x++) {
	if (Arrays.binarySearch(nodes[x].kernels,0,nodes[x].lastKernel+1,fabio)>=0) {
		System.err.println("this node has kernel "+fabio+": "+nodes[x]);
		System.err.print("nodes2kernel: ");
		for (int y=0; y<=lastNode2kernel[x]; y++) System.err.print(nodes2kernel[x][y]+((y+1)%BLOCK_SIZE_PRIME==0?",  ":","));
		System.err.println();
	}
}
}
		
		// Building the kernel graph from $nodes2kernel$.
		lastInNeighbor = new int[nPathKernels];
		Math.set(lastInNeighbor,nPathKernels-1,-1);
		lastOutNeighbor = new int[nPathKernels];
		Math.set(lastOutNeighbor,nPathKernels-1,-1);
		inNeighbors = new int[nPathKernels][0];
		outNeighbors = new int[nPathKernels][0];
		inNeighborFlags = new byte[nPathKernels][0];
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		if (tmpArray2==null || tmpArray2.length<nPathKernels) tmpArray2 = new int[nPathKernels];
		if (tmpArray3==null || tmpArray3.length<nPathKernels) tmpArray3 = new int[nPathKernels];
		if (tmpArray4==null || tmpArray4.length<nPathKernels) tmpArray4 = new int[nPathKernels];
		from=0; arcs=0;
		for (i=1; i<nNodesInKernel; i++) {
			if (blockKernels[i]==blockKernels[i-1]) continue;
			buildKernelGraph_impl(nodes2kernel,lastNode2kernel,nodes,nNodes,from,i-1,blockKernels[from]);
			arcs+=lastInNeighbor[blockKernels[from]]+1;
			from=i;
		}
		buildKernelGraph_impl(nodes2kernel,lastNode2kernel,nodes,nNodes,from,nNodesInKernel-1,blockKernels[from]);
		arcs+=lastInNeighbor[blockKernels[from]]+1;
		System.err.println("buildKernelGraph> kernel graph built ("+arcs+" total arcs)");
		
if (IO.SHOW_INTERACTIVE) {
System.err.print("in-neighbors (1) of kernel path "+fabio+": ");	
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighbors[fabio][x]+",");
System.err.println();
System.err.print("inNeighborFlags: ");
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighborFlags[fabio][4*x+0]+","+inNeighborFlags[fabio][4*x+1]+","+inNeighborFlags[fabio][4*x+2]+","+inNeighborFlags[fabio][4*x+3]+",  ");
System.err.println();
}

		
		// Simplifying the graph
		tmpArcs=buildKernelGraph_removeArcs();  // Arcs with $a_X=0$ for some $X$ (but not for both) might still be present.
		System.err.println("buildKernelGraph> removed "+tmpArcs+" arcs");

		
if (IO.SHOW_INTERACTIVE) {
System.err.print("in-neighbors (2) of kernel path "+fabio+": ");	
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighbors[fabio][x]+",");
System.err.println();
System.err.print("inNeighborFlags: ");
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighborFlags[fabio][4*x+0]+","+inNeighborFlags[fabio][4*x+1]+","+inNeighborFlags[fabio][4*x+2]+","+inNeighborFlags[fabio][4*x+3]+",  ");
System.err.println();
}
		
		tmpArcs=buildKernelGraph_removeSameReadArcs();
		System.err.println("buildKernelGraph> removed "+tmpArcs+" same-read permutation arcs");


if (IO.SHOW_INTERACTIVE) {
System.err.print("in-neighbors (3) of kernel path "+fabio+": ");	
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighbors[fabio][x]+",");
System.err.println();
System.err.print("inNeighborFlags: ");
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighborFlags[fabio][4*x+0]+","+inNeighborFlags[fabio][4*x+1]+","+inNeighborFlags[fabio][4*x+2]+","+inNeighborFlags[fabio][4*x+3]+",  ");
System.err.println();
}

		permutationWindows = new int[nNodes][0];
		lastPermutationWindow = new int[nNodes];
		buildKernelGraph_buildPermutationWindows(permutationWindows,lastPermutationWindow,alignmentsFilePath,nAlignments,nodes,nNodes);
		tmpArcs=buildKernelGraph_removePermutations(nodes,nNodes,permutationWindows,lastPermutationWindow,IDENTITY_THRESHOLD);
		System.err.println("buildKernelGraph> removed "+tmpArcs+" permutation arcs");
		
if (IO.SHOW_INTERACTIVE) {
System.err.print("in-neighbors (3.1) of kernel path "+fabio+": ");	
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighbors[fabio][x]+",");
System.err.println();
System.err.print("inNeighborFlags: ");
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighborFlags[fabio][4*x+0]+","+inNeighborFlags[fabio][4*x+1]+","+inNeighborFlags[fabio][4*x+2]+","+inNeighborFlags[fabio][4*x+3]+",  ");
System.err.println();
}
		
		tmpArcs=buildKernelGraph_removePermutations_longPeriods();
		System.err.println("buildKernelGraph> removed "+tmpArcs+" long-period substring arcs");
		
if (IO.SHOW_INTERACTIVE) {
System.err.print("in-neighbors (3.2) of kernel path "+fabio+": ");	
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighbors[fabio][x]+",");
System.err.println();
System.err.print("inNeighborFlags: ");
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighborFlags[fabio][4*x+0]+","+inNeighborFlags[fabio][4*x+1]+","+inNeighborFlags[fabio][4*x+2]+","+inNeighborFlags[fabio][4*x+3]+",  ");
System.err.println();
}		
		
		tmpArcs=buildKernelGraph_removeDeletions(nodes,nNodes,permutationWindows,lastPermutationWindow,IDENTITY_THRESHOLD);
		System.err.println("buildKernelGraph> removed "+tmpArcs+" deletion arcs");
		
if (IO.SHOW_INTERACTIVE) {
System.err.print("in-neighbors (3.3) of kernel path "+fabio+": ");	
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighbors[fabio][x]+",");
System.err.println();
System.err.print("inNeighborFlags: ");
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighborFlags[fabio][4*x+0]+","+inNeighborFlags[fabio][4*x+1]+","+inNeighborFlags[fabio][4*x+2]+","+inNeighborFlags[fabio][4*x+3]+",  ");
System.err.println();
}
		
		getFrequency_buildDescendantsList(true,tmpArray1,tmpArray2);
		tmpArcs=buildKernelGraph_addSameReadArcs_prime();
		System.err.println("buildKernelGraph> added "+tmpArcs+" same-read arcs (2)");
		
if (IO.SHOW_INTERACTIVE) {
System.err.print("in-neighbors (4) of kernel path "+fabio+": ");	
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighbors[fabio][x]+",");
System.err.println();
System.err.print("inNeighborFlags: ");
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighborFlags[fabio][4*x+0]+","+inNeighborFlags[fabio][4*x+1]+","+inNeighborFlags[fabio][4*x+2]+","+inNeighborFlags[fabio][4*x+3]+",  ");
System.err.println();
}
		
		tmpArcs=buildKernelGraph_removeInsertions();
		System.err.println("buildKernelGraph> removed insertions in "+tmpArcs+" arcs");
		
if (IO.SHOW_INTERACTIVE) {
System.err.print("in-neighbors (5) of kernel path "+fabio+": ");	
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighbors[fabio][x]+",");
System.err.println();
System.err.print("inNeighborFlags: ");
for (int x=0; x<=lastInNeighbor[fabio]; x++) System.err.print(inNeighborFlags[fabio][4*x+0]+","+inNeighborFlags[fabio][4*x+1]+","+inNeighborFlags[fabio][4*x+2]+","+inNeighborFlags[fabio][4*x+3]+",  ");
System.err.println();
}
		
		// Initializing descendants for the following steps of the pipeline
		getFrequency_buildDescendantsList(false,tmpArray1,tmpArray2);
		
		// Expanding $pathKernelPeriodic$ marks using the kernel graph.
		if (tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		if (tmpArray2.length<nPathKernels) tmpArray2 = new int[nPathKernels];
		if (tmpArray3.length<nPathKernels) tmpArray3 = new int[nPathKernels];
		if (out>=0) out+=expandPathKernelPeriodic(tmpArray1,tmpArray2,tmpArray3);
		System.err.println("buildKernelGraph> "+out+" path kernels marked as periodic in total");
		
		permutationWindows=null; lastPermutationWindow=null;
		return out;
	}
	
	
	/**
	 * Builds matrix $nodes2kernel$ which stores, for every node $n$ in $nodes$, the set
	 * of tuples $(a,b,i,o,s,e)$ where $a$ is a kernel of node $n$, and $b$ is a kernel of 
	 * one or more other interval graph nodes (not necessarily in $nodes$), such that the 
	 * full sequence of node $n$ is contained in the full sequence of kernel $b$, or it 
	 * inserts into it (insertion/containment depend on $a$). $i$ stores which end of $a$ 
	 * is seen to insert into $b$, if any (1=none, 2=start, 3=end, 4=both), and $o$ stores
	 * the orientation of $a$ with respect to $b$ (0=same, 1=RC). $s$ ($e$) stores an 
	 * estimate of the max distance of an insertion of the start (end) of $a$ from the 
	 * following end of $b$ (zero if there is no insertion). Tuples are sorted by $a,b,o$.
	 *
	 * Assume that the full sequence of kernel $a$ is indeed contained in the full 
	 * sequence of another kernel $b$. It could happen that the nodes that form the path 
	 * of $a$ in its assembly graph, fragment $a$ in a completely different way from how 
	 * the substring of $b$ that corresponds to $a$ is fragmented by the nodes in its 
	 * assembly graph, or even by the nodes in its basin. We call this \emph{differential
	 * fragmentation}. The procedure tries to handle differential fragmentation by 
	 * collecting all substrings of node $i$ that are mapped by alignments to a node of
	 * kernel $b$, by merging such substrings, and by marking $i$ as contained/inserted 
	 * into $b$ iff most of $i$ is covered by such substrings. 
	 * An example of differential fragentation could be a one-node kernel that is a 
	 * fragment of a multi-node kernel, but which doesn't insert in any node in the basin 
	 * of the container kernel, because of how the container repeat is cut into pieces.
	 * This is actually more likely to happen for maximal nodes by containment, on which 
	 * kernels are built. Because of differential fragmentation, no interval in the basin 
	 * of the contained kernel might belong to the basin of the container kernel.
	 * All alignments are used, without taking into account interval types or the 
	 * assignment of alignments to nodes, since the fact that kernels that are contained 
	 * in one another form different assembly graphs is due to lack of containment edges 
	 * among their nodes, which is most likely caused by inverval type mismatches or wrong 
	 * assignments of alignments to intervals. Without such issues, differential 
	 * fragmentation would not put the sequence of $a$ and $b$ in different kernel graphs.
	 * For this reason, the procedure can be seen as a variant of 
	 * $IntervalGraph.addContainmentEdges()$.
	 *
	 * Remark: the full sequence of node $n$ might be covered by windows with kernel $b$,
	 * but such windows might occur in a different order in the full sequence of $b$, i.e.
	 * $n$ might be a rotation or a permutation of a substring of $b$. The procedure does 
	 * not detect this and assumes that no kernel is a permutation of another kernel.
	 *
	 * Remark: we could have used insertion and shared-substring edges incident to nodes
	 * in $nodes$, rather than all alignments. However, this would not have solved cases 
	 * in which alignments are assigned to wrong intervals.
	 *
	 * Remark: the procedure assumes $IntervalGraph.nodesArray$ to be sorted by $read,
	 * start$.
	 *
	 * Remark: the procedure uses the global variables $shortPeriodAlignments,
	 * shortPeriodWindows,trackPointers$, which should have already been initialized.
	 *
	 * @param nodes[0..nNodes-1] all nodes in a kernel, sorted by $read,start$;
	 * @param alignmentsFilePath file containing all alignments.
	 */
	public static final void buildNodes2Kernel(IntervalGraph.Node[] nodes, int nNodes, String alignmentsFilePath, int[][] nodes2kernel, int[] lastNode2kernel, int totalNAlignments) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int COMPACTION_THRESHOLD = BLOCK_SIZE << 10;  // Arbitrary. Should not be too small, otherwise the procedure becomes very slow.
		final int MAX_GAP = IO.quantum<<1;  // Arbitrary
		final int OUTPUT_EVERY_ALIGNMENTS = 1000;  // Arbitrary
		boolean isConcatenationFrom, isConcatenationTo;
		int i, j, k, p, q;
		int nAlignments, readA, readB, nodeStart, nodeEnd, maxGap, end;
		int kernelFrom, kernelTo, isInsertion, distanceStart, distanceEnd;
		int currentKernelFrom, currentKernelTo, orientation, currentOrientation;
		int insertionDS, insertionDE, operation;
		final int minRead = nodes[0].read;
		final int maxRead = nodes[nNodes-1].read;
		final int minStart = nodes[0].start;
		int maxEnd;
		final int previousOrder = IntervalGraph.Node.order;
		String str;
		BufferedReader alignmentsFile;
		IntervalGraph.Node fromNode, toNode;
		int[] state = new int[7];
		
		// Initializing temporary space
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		if (tmpArray2==null || tmpArray2.length<nPathKernels) tmpArray2 = new int[nPathKernels];
		if (tmpArray3==null || tmpArray3.length<nPathKernels) tmpArray3 = new int[nPathKernels];
		if (tmpWindows==null || tmpWindows.length<nPathKernels) {
			tmpWindows = new Node2KernelWindow[nPathKernels];
			for (i=0; i<nPathKernels; i++) tmpWindows[i] = new Node2KernelWindow(-1,-1,-1,-1,-1,-1,0,0);
		}
		maxEnd=0; i=nNodes-1;
		while (i>=0 && nodes[i].read==maxRead) maxEnd=Math.max(maxEnd,nodes[i--].end);
		
		// Building windows
		Math.set(lastNode2kernel,nNodes-1,-1);
		alignmentsFile = new BufferedReader(new FileReader(alignmentsFilePath),IO.BUFFER_SIZE);
		// Skipping the first two lines
		alignmentsFile.readLine();
		alignmentsFile.readLine();
		str=alignmentsFile.readLine();
		nAlignments=0;
		state[0]=-1; state[1]=-1; state[2]=-1; state[3]=-1; state[4]=-1; state[5]=-1; state[6]=-1;
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		while (str!=null) {
			operation=buildNodes2Kernel_skipAlignment(str,minRead,maxRead,minStart,maxEnd,nodes,nNodes,IntervalGraph.nodesArray,IntervalGraph.nNodes,tmpNode,state);
			if (operation==0) break;
			else if (operation==1) {
				nAlignments++;
				if (nAlignments%OUTPUT_EVERY_ALIGNMENTS==0) System.err.println("buildNodes2Kernel> "+IO.getPercent(nAlignments,totalNAlignments)+"%");
				str=alignmentsFile.readLine();
				continue;
			}
			p=state[0]; q=state[1];
			// Adding windows
			readA=Alignments.readA-1; readB=Alignments.readB-1;
			for (i=p-1; i>=0; i--) {
				fromNode=nodes[i];
				if (fromNode.read!=readA) break;
				if (fromNode.end<=Alignments.startA) continue;
				isConcatenationFrom=(fromNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readA,fromNode.start,fromNode.end,shortPeriodWindows,trackPointers);
				for (j=q-1; j>=0; j--) {
					toNode=IntervalGraph.nodesArray[j];
					if (toNode.read!=readB) break;
					if (toNode.end<=Alignments.startB || toNode.lastKernel==-1) continue;
					if (shortPeriodAlignments[nAlignments] && (fromNode.type!=Constants.INTERVAL_PERIODIC || fromNode.hasLongPeriod || toNode.type!=Constants.INTERVAL_PERIODIC || toNode.hasLongPeriod)) continue;
					isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readB,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
					buildNodes2Kernel_addWindows(fromNode,i,isConcatenationFrom,toNode,isConcatenationTo,nodes2kernel,lastNode2kernel,IDENTITY_THRESHOLD);
					if (lastNode2kernel[i]+1>COMPACTION_THRESHOLD) buildNodes2Kernel_compactWindows(nodes2kernel,lastNode2kernel,i,IDENTITY_THRESHOLD);
				}
				for (j=q; j<IntervalGraph.nNodes; j++) {
					toNode=IntervalGraph.nodesArray[j];
					if (toNode.read!=readB || toNode.start>=Alignments.endB) break;
					if (toNode.lastKernel==-1) continue;
					if (shortPeriodAlignments[nAlignments] && (fromNode.type!=Constants.INTERVAL_PERIODIC || fromNode.hasLongPeriod || toNode.type!=Constants.INTERVAL_PERIODIC || toNode.hasLongPeriod)) continue;
					isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readB,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
					buildNodes2Kernel_addWindows(fromNode,i,isConcatenationFrom,toNode,isConcatenationTo,nodes2kernel,lastNode2kernel,IDENTITY_THRESHOLD);
					if (lastNode2kernel[i]+1>COMPACTION_THRESHOLD) buildNodes2Kernel_compactWindows(nodes2kernel,lastNode2kernel,i,IDENTITY_THRESHOLD);
				}
			}
			for (i=p; i<nNodes; i++) {
				fromNode=nodes[i];
				if (fromNode.read!=readA || fromNode.start>=Alignments.endA) break;
				isConcatenationFrom=(fromNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readA,fromNode.start,fromNode.end,shortPeriodWindows,trackPointers);
				for (j=q-1; j>=0; j--) {
					toNode=IntervalGraph.nodesArray[j];
					if (toNode.read!=readB) break;
					if (toNode.end<=Alignments.startB || toNode.lastKernel==-1) continue;
					if (shortPeriodAlignments[nAlignments] && (fromNode.type!=Constants.INTERVAL_PERIODIC || fromNode.hasLongPeriod || toNode.type!=Constants.INTERVAL_PERIODIC || toNode.hasLongPeriod)) continue;
					isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readB,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
					buildNodes2Kernel_addWindows(fromNode,i,isConcatenationFrom,toNode,isConcatenationTo,nodes2kernel,lastNode2kernel,IDENTITY_THRESHOLD);
					if (lastNode2kernel[i]+1>COMPACTION_THRESHOLD) buildNodes2Kernel_compactWindows(nodes2kernel,lastNode2kernel,i,IDENTITY_THRESHOLD);
				}
				for (j=q; j<IntervalGraph.nNodes; j++) {
					toNode=IntervalGraph.nodesArray[j];
					if (toNode.read!=readB || toNode.start>=Alignments.endB) break;
					if (toNode.lastKernel==-1) continue;
					if (shortPeriodAlignments[nAlignments] && (fromNode.type!=Constants.INTERVAL_PERIODIC || fromNode.hasLongPeriod || toNode.type!=Constants.INTERVAL_PERIODIC || toNode.hasLongPeriod)) continue;
					isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readB,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
					buildNodes2Kernel_addWindows(fromNode,i,isConcatenationFrom,toNode,isConcatenationTo,nodes2kernel,lastNode2kernel,IDENTITY_THRESHOLD);
					if (lastNode2kernel[i]+1>COMPACTION_THRESHOLD) buildNodes2Kernel_compactWindows(nodes2kernel,lastNode2kernel,i,IDENTITY_THRESHOLD);
				}
			}
			nAlignments++;
			if (nAlignments%OUTPUT_EVERY_ALIGNMENTS==0) System.err.println("buildNodes2Kernel> "+IO.getPercent(nAlignments,totalNAlignments)+"%");
			str=alignmentsFile.readLine();
		}
		alignmentsFile.close();
		IntervalGraph.Node.order=previousOrder;

int fabio = 22;
if (IO.SHOW_INTERACTIVE) {		
for (int x=0; x<nNodes; x++) {
	if (Arrays.binarySearch(nodes[x].kernels,0,nodes[x].lastKernel+1,fabio)>=0) {
		System.err.println("buildNodes2Kernel> windows of node "+x+":  (maxReadLength="+maxReadLength+")   node="+nodes[x]);
		for (int y=0; y<=lastNode2kernel[x]; y++) System.err.print(nodes2kernel[x][y]+((y+1)%BLOCK_SIZE==0?",  ":","));
		System.err.println();
	}
}
}
		
		
		// Transforming windows into kernel IDs
		for (i=0; i<IntervalGraph.nNodes; i++) {  // Backing up $nodeID$.
			fromNode=IntervalGraph.nodesArray[i];
			fromNode.visited=fromNode.nodeID; fromNode.nodeID=i;
		}
		for (i=0; i<nNodes; i++) {
			if (nodes[i].lastKernel>0) buildNodes2Kernel_addWindows_sameNode(nodes[i],i,nodes2kernel,lastNode2kernel);
			buildNodes2Kernel_addWindows_sameRead(nodes[i],i,nodes2kernel,lastNode2kernel,IDENTITY_THRESHOLD);
		}
		for (i=0; i<IntervalGraph.nNodes; i++) {  // Restoring $nodeID$.
			fromNode=IntervalGraph.nodesArray[i]; fromNode.nodeID=fromNode.visited;
		}
		for (i=0; i<nNodes; i++) {
			if (lastNode2kernel[i]==-1) continue;
			if (lastNode2kernel[i]>0) buildNodes2Kernel_compactWindows(nodes2kernel,lastNode2kernel,i,IDENTITY_THRESHOLD);
			
if (IO.SHOW_INTERACTIVE) {
if (Arrays.binarySearch(nodes[i].kernels,0,nodes[i].lastKernel+1,fabio)>=0) {
	System.err.println("buildNodes2Kernel> compacted windows of node "+nodes[i]+":  (maxReadLength="+maxReadLength+")");
	for (int y=0; y<=lastNode2kernel[i]; y++) System.err.print(nodes2kernel[i][y]+((y+1)%BLOCK_SIZE==0?",  ":","));
	System.err.println();
}
}
			
			
			nodeStart=nodes[i].start; nodeEnd=nodes[i].end;
			k=-1;
			kernelFrom=nodes2kernel[i][0];
			kernelTo=nodes2kernel[i][1];
			distanceStart=nodes2kernel[i][2]-nodeStart;
			maxGap=distanceStart;
			end=nodes2kernel[i][3];
			isInsertion=nodes2kernel[i][4];
			orientation=nodes2kernel[i][5];
			insertionDS=nodes2kernel[i][6];
			insertionDE=nodes2kernel[i][7];
			for (j=BLOCK_SIZE; j<=lastNode2kernel[i]; j+=BLOCK_SIZE) {
				currentKernelFrom=nodes2kernel[i][j];
				currentKernelTo=nodes2kernel[i][j+1];
				currentOrientation=nodes2kernel[i][j+5];
				if (currentKernelFrom!=kernelFrom || currentKernelTo!=kernelTo || currentOrientation!=orientation) {
					distanceEnd=nodeEnd-end;
					maxGap=Math.max(maxGap,distanceEnd);
					if (maxGap<=MAX_GAP) {
						nodes2kernel[i][++k]=kernelFrom;
						nodes2kernel[i][++k]=kernelTo;
						// If the leftmost/rightmost windows are too far from the ends of
						// the node, we force them to be insertions. The correct way to 
						// do it would be to add a new "unknown" state, and to handle it
						// downstream. However, such unknown state would be treated as if
						// it were an insertion by all procedures downstream.
						if (distanceStart>IDENTITY_THRESHOLD) {
							isInsertion=addInsertions[isInsertion][2];
							insertionDS=maxReadLength;
						}
						if (distanceEnd>IDENTITY_THRESHOLD) {
							isInsertion=addInsertions[isInsertion][3];
							insertionDE=maxReadLength;
						}
						nodes2kernel[i][++k]=isInsertion;
						nodes2kernel[i][++k]=orientation;
						nodes2kernel[i][++k]=insertionDS;
						nodes2kernel[i][++k]=insertionDE;
					}
					kernelFrom=currentKernelFrom;
					kernelTo=currentKernelTo;
					distanceStart=nodes2kernel[i][j+2]-nodeStart;
					maxGap=distanceStart;
					end=nodes2kernel[i][j+3];
					isInsertion=nodes2kernel[i][j+4];
					orientation=currentOrientation;
					insertionDS=nodes2kernel[i][j+6];
					insertionDE=nodes2kernel[i][j+7];
					continue;
				}
				maxGap=Math.max(maxGap,nodes2kernel[i][j+2]-end);  // Windows do not overlap
				end=nodes2kernel[i][j+3];  // Windows do not overlap
				isInsertion=addInsertions[isInsertion][nodes2kernel[i][j+4]];
				insertionDS=Math.max(insertionDS,nodes2kernel[i][j+6]);
				insertionDE=Math.max(insertionDE,nodes2kernel[i][j+7]);
			}
			distanceEnd=nodeEnd-end;
			maxGap=Math.max(maxGap,nodeEnd-end);
			if (maxGap<=MAX_GAP) {
				nodes2kernel[i][++k]=kernelFrom;
				nodes2kernel[i][++k]=kernelTo;
				// See comment above on leftmost/rightmost windows that are too far from
				// the ends of the node.
				if (distanceStart>IDENTITY_THRESHOLD) {
					isInsertion=addInsertions[isInsertion][2];
					insertionDS=maxReadLength;
				}
				if (distanceEnd>IDENTITY_THRESHOLD) {
					isInsertion=addInsertions[isInsertion][3];
					insertionDE=maxReadLength;
				}
				nodes2kernel[i][++k]=isInsertion;
				nodes2kernel[i][++k]=orientation;
				nodes2kernel[i][++k]=insertionDS;
				nodes2kernel[i][++k]=insertionDE;
			}
			lastNode2kernel[i]=k;
		}
	}
	
	
	/**
	 * Tries to skip an alignment, and to replace global binary searches with linear 
	 * scans of few consecutive elements. This is useful e.g. in periodic components with 
	 * a large number of alignments per interval graph node.
	 *
	 * @param state 0=p; 1=q; 2=blockReadA; 3=blockStartA; 4=blockReadB; 5=blockStartB; 
	 * 6=previousReadA;
	 * @param tmpNode temporary space;
	 * @return 0=skip the alignment and break the loop in the caller;
	 * 1=skip the alignment and continue the loop in the caller;
	 * 2=do not skip the alignment and use the returned values of $p,q$, which are assumed
	 * to be the output of $Arrays.binarySearch()$ when searching for each range of the 
	 * alignment in the corresponding list of nodes.
     */
	private static final int buildNodes2Kernel_skipAlignment(String str, int minRead, int maxRead, int minStart, int maxEnd, IntervalGraph.Node[] nodesA, int nNodesA, IntervalGraph.Node[] nodesB, int nNodesB, IntervalGraph.Node tmpNode, int[] state) {
		int blockReadA = state[2];
		int blockStartA = state[3];
		int blockReadB = state[4];
		int blockStartB = state[5];
		int previousReadA = state[6];
		int i;
		int p, q, readA, readB;
		
		i=Alignments.readAlignmentFile_readA(str);
		readA=Alignments.readA-1;
		if (readA>maxRead) return 0;
		if (readA<minRead) return 1;
		Alignments.readAlignmentFile_rest(str,i);
		if ((readA==minRead && Alignments.endA<=minStart) || (readA==maxRead && Alignments.startA>=maxEnd)) return 1;
		if (readA==blockReadA) {
			p=blockStartA;
			while (p<nNodesA) {
				if (nodesA[p].read!=readA || nodesA[p].start>=Alignments.startA) break;
				p++;
			}
		}
		else {
			p=blockStartA+1;
			while (p<nNodesA && nodesA[p].read<readA) p++;
			if (p==nNodesA || nodesA[p].read!=readA) return 1;
			blockReadA=readA; blockStartA=p;
			state[2]=blockReadA; state[3]=blockStartA;
			while (p<nNodesA) {
				if (nodesA[p].read!=readA || nodesA[p].start>=Alignments.startA) break;
				p++;
			}
		}
		readB=Alignments.readB-1;
		if (readB==blockReadB) {
			q=blockStartB;
			while (q<nNodesB) {
				if (nodesB[q].read!=readB || nodesB[q].start>=Alignments.startB) break;
				q++;
			}
		}
		else {
			if (readA==previousReadA) {
				q=blockStartB+1;
				while (q<nNodesB && nodesB[q].read<readB) q++;
				if (q==nNodesB || nodesB[q].read!=readB) return 1;
				blockReadB=readB; blockStartB=q;
				state[4]=blockReadB; state[5]=blockStartB;
				while (q<nNodesB) {
					if (nodesB[q].read!=readB || nodesB[q].start>=Alignments.startB) break;
					q++;
				}
			}
			else {
				tmpNode.read=readB; tmpNode.start=Alignments.startB; tmpNode.id=-1;
				q=Arrays.binarySearch(nodesB,0,nNodesB,tmpNode);
				if (q<0) q=Math.min(-1-q,nNodesB-1);
				if (nodesB[q].read!=readB && (q>0?(nodesB[q-1].read!=readB):true)) return 1;
				previousReadA=readA; blockReadB=readB; 
				state[6]=previousReadA; state[4]=blockReadB;
				blockStartB=q;
				while (blockStartB>0 && nodesB[blockStartB].read>=readB) blockStartB--;
				blockStartB++;
				state[5]=blockStartB;
			}
		}
		state[0]=p; state[1]=q;
		return 2;
	}
	
	
	/**
	 * Let $fromNode$ be a kernel node, assume that $Alignments$ stores an alignment with 
	 * nonempty intersection $[x..y]$ with $fromNode$ and nonempty intersection with 
	 * $toNode$. The procedure appends to $nodes2kernel[fromNode]$ a tuple 
	 * $(a,b,x,y,i,o,s,e)$ for every pair $(a,b)$ where $a \in fromNode.kernels$ and 
	 * $b \in toNode.kernels \setminus fromNode.kernels$. Depending on $a,b$, such copy 
	 * can represent an insertion or a containment: $i$ stores which end of $a$ is seen to
	 * insert into $b$, if any (1=none, 2=start, 3=end, 4=both); $o$ stores the 
	 * orientation of $a$ with respect to $b$ implied by the window (0=forward, 1=RC, 
	 * 2=both); $s$ ($e$) stores the distance of the start (end) of $a$ from a following 
	 * end of $b$, if this is an insertion (0 otherwise).
	 *
	 * Remark: the procedure does not use kernel end information in $pathsWith*$ arrays,
	 * i.e. all the elements of such arrays are internally transformed to non-negative.
	 *
	 * Remark: the procedure sets $o=0$ for every tuple with $a$ or $b$ cyclic. This is 
	 * useful, since the real orientations of the nodes of a cyclic kernel WRT a reference
	 * sequence of the kernel is unknown. Having the same orientation helps windows to 
	 * cover the full surface of every node, and to get a nonempty intersection when 
	 * comparing the lists of kernels of different bidirected graph nodes in
	 * $buildKernelGraph_impl_allFlags()$.
	 *
	 * Remark: the procedure just appends windows, without compacting them.
	 *
	 * Remark: the procedure assumes that global variables $tmpArray{1,2,3}$ have length 
	 * at least equal to $fromNode.lastKernel+1$.
	 *
	 * @param fromNodeIndex position of $fromNode$ in $nodes2kernel$;
	 * @param *NodeIsConcatenation tells whether $*Node$ is a concatenation of several
	 * distinct short-period intervals;
	 * @param threshold identity threshold.
	 */
	private static final void buildNodes2Kernel_addWindows(IntervalGraph.Node fromNode, int fromNodeIndex, boolean fromNodeIsConcatenation, IntervalGraph.Node toNode, boolean toNodeIsConcatenation, int[][] nodes2kernel, int[] lastNode2kernel, int threshold) {
		final int MIN_WINDOW_SIZE = Alignments.minAlignmentLength>>2;  // Arbitrary
		boolean projectionStartInside, projectionEndInside;
		int windowStart, windowEnd, windowStartPrime, windowEndPrime;
		final int readFrom = fromNode.read;
		final int readTo = toNode.read;
		
		windowStart=Math.max(Alignments.startA,fromNode.start);
		windowEnd=Math.min(Alignments.endA,fromNode.end);
		Intervals.project(toNode.start,toNode.end,Alignments.startB,Alignments.endB,Alignments.startA,Alignments.endA,Alignments.orientation,tmpIO,0);
		windowStart=Math.max(windowStart,tmpIO[0]);
		windowEnd=Math.min(windowEnd,tmpIO[1]);
		if (windowEnd-windowStart+1<MIN_WINDOW_SIZE) {
			// Remark: it can be $windowStart>windowEnd$ since e.g. $[tmpIO[0]..tmpIO[1]]$
			// might have no intersection with $fromNode$.
			return;
		}
		// Discarding the window if it belongs to a short-period interval
		if ((fromNode.type!=Constants.INTERVAL_PERIODIC || fromNode.hasLongPeriod) && !fromNodeIsConcatenation && IntervalGraph.inShortPeriodTrack(readFrom,windowStart,windowEnd,shortPeriodWindows,trackPointers)) return;
		if ((toNode.type!=Constants.INTERVAL_PERIODIC || toNode.hasLongPeriod) && !toNodeIsConcatenation) {
			windowStartPrime=Math.max(Alignments.startB,toNode.start);
			windowEndPrime=Math.min(Alignments.endB,toNode.end);
			Intervals.project(fromNode.start,fromNode.end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,tmpIO,0);
			windowStartPrime=Math.max(windowStartPrime,tmpIO[0]);
			windowEndPrime=Math.min(windowEndPrime,tmpIO[1]);
			if (IntervalGraph.inShortPeriodTrack(readTo,windowStartPrime,windowEndPrime,shortPeriodWindows,trackPointers)) return;
		}
		Intervals.project(windowStart,windowEnd,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,tmpIO,0);
		if (Alignments.orientation) {
			// No need to use the output of $project()$ in some cases.
			if (Math.abs(windowStart,Alignments.startA)<=threshold && Math.abs(Alignments.startB,toNode.start)<=threshold) projectionStartInside=false;
			else projectionStartInside=tmpIO[0]>toNode.start+threshold;
			if (Math.abs(windowEnd,Alignments.endA)<=threshold && Math.abs(Alignments.endB,toNode.end)<=threshold) projectionEndInside=false;
			else projectionEndInside=tmpIO[1]<toNode.end-threshold;
		}
		else {
			// No need to use the output of $project()$ in some cases.
			if (Math.abs(windowEnd,Alignments.endA)<=threshold && Math.abs(Alignments.startB,toNode.start)<=threshold) projectionStartInside=false;
			else projectionStartInside=tmpIO[0]>toNode.start+threshold;
			if (Math.abs(windowStart,Alignments.startA)<=threshold && Math.abs(Alignments.endB,toNode.end)<=threshold) projectionEndInside=false;
			else projectionEndInside=tmpIO[1]<toNode.end-threshold;
		}		
		buildNodes2Kernel_addWindows_impl(windowStart,windowEnd,Alignments.orientation,projectionStartInside,projectionEndInside,fromNode,fromNodeIndex,toNode,nodes2kernel,lastNode2kernel,threshold);
	}

	
	/**
	 * Variant of $buildNodes2Kernel_addWindows()$ for the case in which $fromNode$ and
	 * $toNode$ belong to the same read.
	 */
	private static final void buildNodes2Kernel_addWindows_sameRead_impl(IntervalGraph.Node fromNode, int fromNodeIndex, boolean fromNodeIsConcatenation, IntervalGraph.Node toNode, boolean toNodeIsConcatenation, int[][] nodes2kernel, int[] lastNode2kernel, int threshold) {
		final int MIN_WINDOW_SIZE = Alignments.minAlignmentLength>>2;  // Arbitrary
		boolean projectionStartInside, projectionEndInside;
		int windowStart, windowEnd;
		final int readFrom = fromNode.read;
		
		windowStart=Math.max(fromNode.start,toNode.start);
		windowEnd=Math.min(fromNode.end,toNode.end);
		if (windowEnd-windowStart+1<MIN_WINDOW_SIZE) return;
		// Discarding the window if it belongs to a short-period interval
		if ( ( ((fromNode.type!=Constants.INTERVAL_PERIODIC || fromNode.hasLongPeriod) && !fromNodeIsConcatenation) ||
			   ((toNode.type!=Constants.INTERVAL_PERIODIC || toNode.hasLongPeriod) && !toNodeIsConcatenation)
			 ) && 
			 IntervalGraph.inShortPeriodTrack(readFrom,windowStart,windowEnd,shortPeriodWindows,trackPointers)
		   ) return;
		projectionStartInside=windowStart>toNode.start+threshold;
		projectionEndInside=windowEnd<toNode.end-threshold;		
		buildNodes2Kernel_addWindows_impl(windowStart,windowEnd,true,projectionStartInside,projectionEndInside,fromNode,fromNodeIndex,toNode,nodes2kernel,lastNode2kernel,threshold);
	}
	
	
	private static final void buildNodes2Kernel_addWindows_impl(int windowStart, int windowEnd, boolean alignmentOrientation, boolean projectionStartInside, boolean projectionEndInside, IntervalGraph.Node fromNode, int fromNodeIndex, IntervalGraph.Node toNode, int[][] nodes2kernel, int[] lastNode2kernel, int threshold) {
		boolean toNodeStartHasKernelTo, toNodeEndHasKernelTo;
		byte kernelFromInStart, kernelFromInEnd, isInsertion, kernelsOrientation, orientationFrom, orientationTo;
		int i, j;
		int lastKernel, kernelFrom, kernelTo;
		final int lastKernelTo = Math.setMinus(toNode.kernels,toNode.lastKernel,fromNode.kernels,fromNode.lastKernel,tmpArray3);
		
		if (windowStart<=fromNode.start+threshold && windowEnd>=fromNode.end-threshold) {
			lastKernel=Math.setUnion(fromNode.pathsWithStart,fromNode.lastPathWithStart,fromNode.pathsWithEnd,fromNode.lastPathWithEnd,tmpArray1);
			lastKernel=Math.makePositive(tmpArray1,0,lastKernel);
			for (i=0; i<=lastKernel; i++) {
				kernelFrom=tmpArray1[i];
				orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
				if (orientationFrom==-1) continue;
				kernelFromInStart=1;  // 1=no; 2=start of kernel; 3=end of kernel; 4=both.
				if (fromNode.lastPathWithStart>=0) {
					if (Arrays.binarySearch(fromNode.pathsWithStart,0,fromNode.lastPathWithStart+1,kernelFrom)>=0) kernelFromInStart=2;
					if (Arrays.binarySearch(fromNode.pathsWithStart,0,fromNode.lastPathWithStart+1,-1-kernelFrom)>=0) kernelFromInStart=addInsertions[kernelFromInStart][3];
				}
				kernelFromInEnd=1;  // 1=no; 2=start of kernel; 3=end of kernel; 4=both.
				if (fromNode.lastPathWithEnd>=0) {
					if (Arrays.binarySearch(fromNode.pathsWithEnd,0,fromNode.lastPathWithEnd+1,kernelFrom)>=0) kernelFromInEnd=2;
					if (Arrays.binarySearch(fromNode.pathsWithEnd,0,fromNode.lastPathWithEnd+1,-1-kernelFrom)>=0) kernelFromInEnd=addInsertions[kernelFromInEnd][3];
				}
				for (j=0; j<=lastKernelTo; j++) {
					kernelTo=tmpArray3[j];
					orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
					if (orientationTo==-1) continue;
					toNodeStartHasKernelTo=toNode.lastPathWithStart>=0 && (Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,kernelTo)>=0 || Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,-1-kernelTo)>=0);
					toNodeEndHasKernelTo=toNode.lastPathWithEnd>=0 && (Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,kernelTo)>=0 || Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,-1-kernelTo)>=0);
					isInsertion=1;  // 1=no; 2=start of kernel; 3=end of kernel; 4=both.					
					if ( (alignmentOrientation && (projectionStartInside || !toNodeStartHasKernelTo)) ||
						 (!alignmentOrientation && (projectionEndInside || !toNodeEndHasKernelTo))
					   ) isInsertion=kernelFromInStart;
					if ( (alignmentOrientation && (projectionEndInside || !toNodeEndHasKernelTo)) ||
					     (!alignmentOrientation && (projectionStartInside || !toNodeStartHasKernelTo))
					   ) isInsertion=addInsertions[isInsertion][kernelFromInEnd];
					if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
					else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];
					buildNodes2Kernel_getInsertionDistance(alignmentOrientation,kernelFromInStart,kernelFromInEnd,toNode.start,toNode.end,toNodeStartHasKernelTo,toNodeEndHasKernelTo,tmpIO[0],tmpIO[1],projectionStartInside,projectionEndInside);
					
					
					
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 1  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"] orient="+alignmentOrientation+"  isInsertion="+isInsertion);
	System.err.println("FUCKKKKKK> 1  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 1  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 1  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 1  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 1  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 1  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.println("FUCKKKKKK> 1  kernelFromInStart="+kernelFromInStart);
	System.err.println("FUCKKKKKK> 1  kernelFromInEnd="+kernelFromInEnd);
	System.err.print("FUCKKKKKK> 1  fromNode.pathsWithStart=");
	for (int x=0; x<=fromNode.lastPathWithStart; x++) System.err.print(fromNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 1  fromNode.pathsWithEnd=");
	for (int x=0; x<=fromNode.lastPathWithEnd; x++) System.err.print(fromNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 1  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 1  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 1  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
/*	System.err.println("FUCKKKKKK> 1  containment neighbors of toNode: ");
	for (int x=0; x<IntervalGraph.nNeighbors[toNode.nodeID]; x++) {
		IntervalGraph.Edge tmpEdge = IntervalGraph.neighbors[toNode.nodeID][x];
		if (tmpEdge.isType_containment()) System.err.println(tmpEdge);  //   IntervalGraph.nodesArray[tmpEdge.getTo(toNode.nodeID)]);
	}*/
	System.err.println();
}
					
					
					buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,isInsertion,kernelsOrientation,insertionDistanceStart,insertionDistanceEnd);
				}
			}
			lastKernel=Math.setMinus(fromNode.kernels,fromNode.lastKernel,tmpArray1,lastKernel,tmpArray2);
			for (i=0; i<=lastKernel; i++) {
				kernelFrom=tmpArray2[i];
				orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
				if (orientationFrom==-1) continue;
				for (j=0; j<=lastKernelTo; j++) {
					kernelTo=tmpArray3[j];
					orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
					if (orientationTo==-1) continue;
					if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
					else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];
					
					
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 1.1  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"] orient="+alignmentOrientation);
	System.err.println("FUCKKKKKK> 1.1  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 1.1  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 1.1  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 1.1  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 1.1  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 1.1  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.print("FUCKKKKKK> 1.1  fromNode.pathsWithStart=");
	for (int x=0; x<=fromNode.lastPathWithStart; x++) System.err.print(fromNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 1.1  fromNode.pathsWithEnd=");
	for (int x=0; x<=fromNode.lastPathWithEnd; x++) System.err.print(fromNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 1.1  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 1.1  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 1.1  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
/*	System.err.println("FUCKKKKKK> 1  containment neighbors of toNode: ");
	for (int x=0; x<IntervalGraph.nNeighbors[toNode.nodeID]; x++) {
		IntervalGraph.Edge tmpEdge = IntervalGraph.neighbors[toNode.nodeID][x];
		if (tmpEdge.isType_containment()) System.err.println(tmpEdge);  //   IntervalGraph.nodesArray[tmpEdge.getTo(toNode.nodeID)]);
	}*/
	System.err.println();
}					
					
					
					
					buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,1,kernelsOrientation,0,0);
				}
			}
		}
		else if (windowStart<=fromNode.start+threshold) {
			if (fromNode.lastPathWithStart>=0) {
				System.arraycopy(fromNode.pathsWithStart,0,tmpArray1,0,fromNode.lastPathWithStart+1);
				lastKernel=Math.makePositive(tmpArray1,0,fromNode.lastPathWithStart);
				for (i=0; i<=lastKernel; i++) {
					kernelFrom=tmpArray1[i];
					orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
					if (orientationFrom==-1) continue;
					kernelFromInStart=1;
					if (Arrays.binarySearch(fromNode.pathsWithStart,0,fromNode.lastPathWithStart+1,kernelFrom)>=0) kernelFromInStart=2;
					if (Arrays.binarySearch(fromNode.pathsWithStart,0,fromNode.lastPathWithStart+1,-1-kernelFrom)>=0) kernelFromInStart=addInsertions[kernelFromInStart][3];
					for (j=0; j<=lastKernelTo; j++) {
						kernelTo=tmpArray3[j];
						if (!buildNodes2kernel_useWindow(fromNode,1,toNode,projectionStartInside,projectionEndInside,kernelFrom,kernelTo)) continue;
						orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
						if (orientationTo==-1) continue;
						toNodeStartHasKernelTo=toNode.lastPathWithStart>=0 && (Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,kernelTo)>=0 || Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,-1-kernelTo)>=0);
						toNodeEndHasKernelTo=toNode.lastPathWithEnd>=0 && (Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,kernelTo)>=0 || Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,-1-kernelTo)>=0);
						isInsertion=1;
						if ( (alignmentOrientation && (projectionStartInside || !toNodeStartHasKernelTo)) ||
							 (!alignmentOrientation && (projectionEndInside || !toNodeEndHasKernelTo))
						   ) isInsertion=kernelFromInStart;
						if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
						else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];						
						buildNodes2Kernel_getInsertionDistance(alignmentOrientation,kernelFromInStart,1,toNode.start,toNode.end,toNodeStartHasKernelTo,toNodeEndHasKernelTo,tmpIO[0],tmpIO[1],projectionStartInside,projectionEndInside);
						
						
							
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 3  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"]  isInsertion="+isInsertion);
	System.err.println("FUCKKKKKK> 3  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 3  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 3  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 3  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 3  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 3  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.println("FUCKKKKKK> 3  kernelFromInStart="+kernelFromInStart);
	System.err.print("FUCKKKKKK> 3  fromNode.pathsWithStart=");
	for (int x=0; x<=fromNode.lastPathWithStart; x++) System.err.print(fromNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
	System.err.println("FUCKKKKKK> 3  projectionStartInside="+projectionStartInside+" projectionEndInside="+projectionEndInside);
	System.err.println();
}
						
						buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,isInsertion,kernelsOrientation,insertionDistanceStart,insertionDistanceEnd);
					}
				}
				lastKernel=Math.setMinus(fromNode.kernels,fromNode.lastKernel,tmpArray1,lastKernel,tmpArray2);
				for (i=0; i<=lastKernel; i++) {
					kernelFrom=tmpArray2[i];
					orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
					if (orientationFrom==-1) continue;
					for (j=0; j<=lastKernelTo; j++) {
						kernelTo=tmpArray3[j];
						if (!buildNodes2kernel_useWindow(fromNode,1,toNode,projectionStartInside,projectionEndInside,kernelFrom,kernelTo)) continue;
						orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
						if (orientationTo==-1) continue;
						if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
						else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];
						
						
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 3.1  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"] orient="+alignmentOrientation);
	System.err.println("FUCKKKKKK> 3.1  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 3.1  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 3.1  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 3.1  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 3.1  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 3.1  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.print("FUCKKKKKK> 3.1  fromNode.pathsWithStart=");
	for (int x=0; x<=fromNode.lastPathWithStart; x++) System.err.print(fromNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3.1  fromNode.pathsWithEnd=");
	for (int x=0; x<=fromNode.lastPathWithEnd; x++) System.err.print(fromNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3.1  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3.1  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3.1  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
/*	System.err.println("FUCKKKKKK> 1  containment neighbors of toNode: ");
	for (int x=0; x<IntervalGraph.nNeighbors[toNode.nodeID]; x++) {
		IntervalGraph.Edge tmpEdge = IntervalGraph.neighbors[toNode.nodeID][x];
		if (tmpEdge.isType_containment()) System.err.println(tmpEdge);  //   IntervalGraph.nodesArray[tmpEdge.getTo(toNode.nodeID)]);
	}*/
	System.err.println();
}						
						
						
						buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,1,kernelsOrientation,0,0);
					}
				}
			}
			else {
				for (i=0; i<=fromNode.lastKernel; i++) {
					kernelFrom=fromNode.kernels[i];
					orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
					if (orientationFrom==-1) continue;
					for (j=0; j<=lastKernelTo; j++) {
						kernelTo=tmpArray3[j];
						if (!buildNodes2kernel_useWindow(fromNode,1,toNode,projectionStartInside,projectionEndInside,kernelFrom,kernelTo)) continue;
						orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
						if (orientationTo==-1) continue;
						if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
						else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];
						
						
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 3.2  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"] orient="+alignmentOrientation);
	System.err.println("FUCKKKKKK> 3.2  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 3.2  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 3.2  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 3.2  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 3.2  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 3.2  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.print("FUCKKKKKK> 3.2 fromNode.pathsWithStart=");
	for (int x=0; x<=fromNode.lastPathWithStart; x++) System.err.print(fromNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3.2  fromNode.pathsWithEnd=");
	for (int x=0; x<=fromNode.lastPathWithEnd; x++) System.err.print(fromNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3.2  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3.2  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 3.2  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
/*	System.err.println("FUCKKKKKK> 1  containment neighbors of toNode: ");
	for (int x=0; x<IntervalGraph.nNeighbors[toNode.nodeID]; x++) {
		IntervalGraph.Edge tmpEdge = IntervalGraph.neighbors[toNode.nodeID][x];
		if (tmpEdge.isType_containment()) System.err.println(tmpEdge);  //   IntervalGraph.nodesArray[tmpEdge.getTo(toNode.nodeID)]);
	}*/
	System.err.println();
}
						
						
						
						buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,1,kernelsOrientation,0,0);
					}
				}
			}
		}
		else if (windowEnd>=fromNode.end-threshold) {
			if (fromNode.lastPathWithEnd>=0) {
				System.arraycopy(fromNode.pathsWithEnd,0,tmpArray1,0,fromNode.lastPathWithEnd+1);
				lastKernel=Math.makePositive(tmpArray1,0,fromNode.lastPathWithEnd);
				for (i=0; i<=lastKernel; i++) {
					kernelFrom=tmpArray1[i];
					orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
					if (orientationFrom==-1) continue;
					kernelFromInEnd=1;
					if (Arrays.binarySearch(fromNode.pathsWithEnd,0,fromNode.lastPathWithEnd+1,kernelFrom)>=0) kernelFromInEnd=2;
					if (Arrays.binarySearch(fromNode.pathsWithEnd,0,fromNode.lastPathWithEnd+1,-1-kernelFrom)>=0) kernelFromInEnd=addInsertions[kernelFromInEnd][3];
					for (j=0; j<=lastKernelTo; j++) {
						kernelTo=tmpArray3[j];
						if (!buildNodes2kernel_useWindow(fromNode,0,toNode,projectionStartInside,projectionEndInside,kernelFrom,kernelTo)) continue;
						orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
						if (orientationTo==-1) continue;
						toNodeStartHasKernelTo=toNode.lastPathWithStart>=0 && (Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,kernelTo)>=0 || Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,-1-kernelTo)>=0);
						toNodeEndHasKernelTo=toNode.lastPathWithEnd>=0 && (Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,kernelTo)>=0 || Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,-1-kernelTo)>=0);
						isInsertion=1;
						if ( (alignmentOrientation && (projectionEndInside || toNode.lastPathWithEnd<0 || (Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,kernelTo)<0 && Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,-1-kernelTo)<0))) ||
						     (!alignmentOrientation && (projectionStartInside || toNode.lastPathWithStart<0 || (Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,kernelTo)<0 && Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,-1-kernelTo)<0)))
						   ) isInsertion=kernelFromInEnd;
						if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
						else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];
						buildNodes2Kernel_getInsertionDistance(alignmentOrientation,1,kernelFromInEnd,toNode.start,toNode.end,toNodeStartHasKernelTo,toNodeEndHasKernelTo,tmpIO[0],tmpIO[1],projectionStartInside,projectionEndInside);
						
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 6  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"] orient="+alignmentOrientation+"  isInsertion="+isInsertion);
	System.err.println("FUCKKKKKK> 6  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 6  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 6  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 6  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 6  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 6  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.println("FUCKKKKKK> 6  kernelFromInEnd="+kernelFromInEnd);
	System.err.print("FUCKKKKKK> 6  fromNode.pathsWithEnd=");
	for (int x=0; x<=fromNode.lastPathWithEnd; x++) System.err.print(fromNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
	System.err.println("FUCKKKKKK> 6  projectionEndInside="+projectionEndInside+" projectionStartInside="+projectionStartInside);
	System.err.println();
}
						
						buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,isInsertion,kernelsOrientation,insertionDistanceStart,insertionDistanceEnd);
					}
				}
				lastKernel=Math.setMinus(fromNode.kernels,fromNode.lastKernel,tmpArray1,lastKernel,tmpArray2);
				for (i=0; i<=lastKernel; i++) {
					kernelFrom=tmpArray2[i];
					orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
					if (orientationFrom==-1) continue;
					for (j=0; j<=lastKernelTo; j++) {
						kernelTo=tmpArray3[j];
						if (!buildNodes2kernel_useWindow(fromNode,0,toNode,projectionStartInside,projectionEndInside,kernelFrom,kernelTo)) continue;
						orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
						if (orientationTo==-1) continue;
						if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
						else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];
						
						
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 6.1  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"] orient="+alignmentOrientation);
	System.err.println("FUCKKKKKK> 6.1  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 6.1  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 6.1  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 6.1  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 6.1  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 6.1  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.print("FUCKKKKKK> 6.1  fromNode.pathsWithEnd=");
	for (int x=0; x<=fromNode.lastPathWithEnd; x++) System.err.print(fromNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.1  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.1  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.1  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
	System.err.println("FUCKKKKKK> 6.1  projectionEndInside="+projectionEndInside+" projectionStartInside="+projectionStartInside);
	System.err.println();
}
						
						
						
						buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,1,kernelsOrientation,0,0);
					}
				}
			}
			else {
				for (i=0; i<=fromNode.lastKernel; i++) {
					kernelFrom=fromNode.kernels[i];
					orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
					if (orientationFrom==-1) continue;
					for (j=0; j<=lastKernelTo; j++) {
						kernelTo=tmpArray3[j];
						if (!buildNodes2kernel_useWindow(fromNode,0,toNode,projectionStartInside,projectionEndInside,kernelFrom,kernelTo)) continue;
						orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
						if (orientationTo==-1) continue;
						if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
						else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];
						
						
						
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 6.2  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"] orient="+alignmentOrientation);
	System.err.println("FUCKKKKKK> 6.2  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 6.2  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 6.2  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 6.2  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 6.2  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 6.2  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.print("FUCKKKKKK> 6.2  fromNode.pathsWithEnd=");
	for (int x=0; x<=fromNode.lastPathWithEnd; x++) System.err.print(fromNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.2  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.2  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.2  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
	System.err.println("FUCKKKKKK> 6.2  projectionEndInside="+projectionEndInside+" projectionStartInside="+projectionStartInside);
	System.err.println();
}
						
						
						
						
						buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,1,kernelsOrientation,0,0);
					}
				}
			}
		}
		else {
			for (i=0; i<=fromNode.lastKernel; i++) {
				kernelFrom=fromNode.kernels[i];				
				orientationFrom=fromNode.kernelOrientations[Arrays.binarySearch(fromNode.kernels,0,fromNode.lastKernel+1,kernelFrom)];
				if (orientationFrom==-1) continue;
				for (j=0; j<=lastKernelTo; j++) {
					kernelTo=tmpArray3[j];
					if (!buildNodes2kernel_useWindow(fromNode,2,toNode,projectionStartInside,projectionEndInside,kernelFrom,kernelTo)) continue;
					orientationTo=toNode.kernelOrientations[Arrays.binarySearch(toNode.kernels,0,toNode.lastKernel+1,kernelTo)];
					if (orientationTo==-1) continue;
					if (pathKernelLengths[kernelFrom]==0 || pathKernelLengths[kernelTo]==0) kernelsOrientation=0;
					else kernelsOrientation=getKernelsOrientation[orientationFrom][orientationTo][alignmentOrientation?0:1];
					
					
if (IO.SHOW_INTERACTIVE && /*isInsertion!=1 &&*/ kernelFrom==22 && kernelTo==7) {
	System.err.println("FUCKKKKKK> 6.3  "+Alignments.readA+"["+Alignments.startA+".."+Alignments.endA+"] x "+Alignments.readB+"["+Alignments.startB+".."+Alignments.endB+"] orient="+alignmentOrientation);
	System.err.println("FUCKKKKKK> 6.3  insertionDistanceStart="+insertionDistanceStart+" insertionDistanceEnd="+insertionDistanceEnd);
	System.err.println("FUCKKKKKK> 6.3  fromNode: "+fromNode);
	System.err.println("FUCKKKKKK> 6.3  orientationFrom: "+orientationFrom);
	System.err.println("FUCKKKKKK> 6.3  toNode: "+toNode);
	System.err.println("FUCKKKKKK> 6.3  orientationTo: "+orientationTo);
	System.err.println("FUCKKKKKK> 6.3  windowStart="+windowStart+" windowEnd="+windowEnd+" kernelsOrientation="+kernelsOrientation);
	System.err.print("FUCKKKKKK> 6.3  fromNode.pathsWithEnd=");
	for (int x=0; x<=fromNode.lastPathWithEnd; x++) System.err.print(fromNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.3  toNode.pathsWithStart=");
	for (int x=0; x<=toNode.lastPathWithStart; x++) System.err.print(toNode.pathsWithStart[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.3  toNode.pathsWithEnd=");
	for (int x=0; x<=toNode.lastPathWithEnd; x++) System.err.print(toNode.pathsWithEnd[x]+",");
	System.err.println();
	System.err.print("FUCKKKKKK> 6.3  toNode.kernels=");
	for (int x=0; x<=toNode.lastKernel; x++) System.err.print(toNode.kernels[x]+",");
	System.err.println();
	System.err.println("FUCKKKKKK> 6.3  projectionEndInside="+projectionEndInside+" projectionStartInside="+projectionStartInside);
	System.err.println();
}
					
					
					
					buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelFrom,kernelTo,windowStart,windowEnd,1,kernelsOrientation,0,0);
				}
			}
		}
	}
	
	
	/**
	 * Sets the global variables $insertionDistance{Start,End}$, which store the distance
	 * between the start (resp. end) of $kernelFrom$, and the following start or end of 
	 * $kernelTo$, as witnessed by an alignment between interval graph nodes $fromNode$
	 * and $toNode$. Insertion distance is set to zero if there is no insertion, and to
	 * $maxReadLength$ if the alignment maps an end of $kernelFrom$ to a position of 
	 * $toNode$ whose distance from the end of $kernelTo$ cannot be estimated.
	 *
	 * See $buildNodes2Kernel_addWindows()$ for the meaning of input arguments.
	 */
	private static final void buildNodes2Kernel_getInsertionDistance(boolean alignmentOrientation, int kernelFromInStart, int kernelFromInEnd, int toNodeStart, int toNodeEnd, boolean toNodeStartHasKernelTo, boolean toNodeEndHasKernelTo, int projectionStart, int projectionEnd, boolean projectionStartInside, boolean projectionEndInside) {
		int distance;
		
		insertionDistanceStart=0; insertionDistanceEnd=0;
		if (alignmentOrientation) {
			if (projectionStartInside) {
				if (toNodeStartHasKernelTo) {
					distance=projectionStart-toNodeStart;
					if (kernelFromInStart==2 || kernelFromInStart==4) insertionDistanceStart=Math.max(insertionDistanceStart,distance);
					if (kernelFromInStart==3 || kernelFromInStart==4) insertionDistanceEnd=Math.max(insertionDistanceEnd,distance);
				}
				else {
					if (kernelFromInStart==2 || kernelFromInStart==4) insertionDistanceStart=maxReadLength;
					if (kernelFromInStart==3 || kernelFromInStart==4) insertionDistanceEnd=maxReadLength;
				}
			}
			else if (!toNodeStartHasKernelTo) {
				if (kernelFromInStart==2 || kernelFromInStart==4) insertionDistanceStart=maxReadLength;
				if (kernelFromInStart==3 || kernelFromInStart==4) insertionDistanceEnd=maxReadLength;
			}
			if (projectionEndInside) {
				if (toNodeEndHasKernelTo) {
					distance=toNodeEnd-projectionEnd;
					if (kernelFromInEnd==2 || kernelFromInEnd==4) insertionDistanceStart=Math.max(insertionDistanceStart,distance);
					if (kernelFromInEnd==3 || kernelFromInEnd==4) insertionDistanceEnd=Math.max(insertionDistanceEnd,distance);
				}
				else {
					if (kernelFromInEnd==2 || kernelFromInEnd==4) insertionDistanceStart=maxReadLength;
					if (kernelFromInEnd==3 || kernelFromInEnd==4) insertionDistanceEnd=maxReadLength;
				}
			}
			else if (!toNodeEndHasKernelTo) {
				if (kernelFromInEnd==2 || kernelFromInEnd==4) insertionDistanceStart=maxReadLength;
				if (kernelFromInEnd==3 || kernelFromInEnd==4) insertionDistanceEnd=maxReadLength;
			}
		}
		else {
			if (projectionStartInside) {
				if (toNodeStartHasKernelTo) {
					distance=projectionStart-toNodeStart;
					if (kernelFromInEnd==2 || kernelFromInEnd==4) insertionDistanceStart=Math.max(insertionDistanceStart,distance);
					if (kernelFromInEnd==3 || kernelFromInEnd==4) insertionDistanceEnd=Math.max(insertionDistanceEnd,distance);
				}
				else {
					if (kernelFromInEnd==2 || kernelFromInEnd==4) insertionDistanceStart=maxReadLength;
					if (kernelFromInEnd==3 || kernelFromInEnd==4) insertionDistanceEnd=maxReadLength;
				}
			}
			else if (!toNodeStartHasKernelTo) {
				if (kernelFromInEnd==2 || kernelFromInEnd==4) insertionDistanceStart=maxReadLength;
				if (kernelFromInEnd==3 || kernelFromInEnd==4) insertionDistanceEnd=maxReadLength;
			}
			if (projectionEndInside) {
				if (toNodeEndHasKernelTo) {
					distance=toNodeEnd-projectionEnd;
					if (kernelFromInStart==2 || kernelFromInStart==4) insertionDistanceStart=Math.max(insertionDistanceStart,distance);
					if (kernelFromInStart==3 || kernelFromInStart==4) insertionDistanceEnd=Math.max(insertionDistanceEnd,distance);
				}
				else {
					if (kernelFromInStart==2 || kernelFromInStart==4) insertionDistanceStart=maxReadLength;
					if (kernelFromInStart==3 || kernelFromInStart==4) insertionDistanceEnd=maxReadLength;
				}
			}
			else if (!toNodeEndHasKernelTo) {
				if (kernelFromInStart==2 || kernelFromInStart==4) insertionDistanceStart=maxReadLength;
				if (kernelFromInStart==3 || kernelFromInStart==4) insertionDistanceEnd=maxReadLength;
			}
		}
	}
	
	
	/**
	 * Adds artificial windows between every pair of kernels in $fromNode.kernels$.
	 */
	private static final void buildNodes2Kernel_addWindows_sameNode(IntervalGraph.Node fromNode, int fromNodeIndex, int[][] nodes2kernel, int[] lastNode2kernel) {
		int i, j;
		int kernelI, orientationI;
		final int windowStart = fromNode.start;
		final int windowEnd = fromNode.end;
		final int nWindows = (fromNode.lastKernel+1)*(fromNode.lastKernel+2);
		
		if (lastNode2kernel[fromNodeIndex]+BLOCK_SIZE*nWindows>=nodes2kernel[fromNodeIndex].length) {
			int[] tmpArray = new int[lastNode2kernel[fromNodeIndex]+1+BLOCK_SIZE*nWindows];
			System.arraycopy(nodes2kernel[fromNodeIndex],0,tmpArray,0,lastNode2kernel[fromNodeIndex]+1);
			nodes2kernel[fromNodeIndex]=tmpArray;
		}
		for (i=0; i<=fromNode.lastKernel; i++) {
			kernelI=fromNode.kernels[i]; orientationI=fromNode.kernelOrientations[i];
			for (j=0; j<=fromNode.lastKernel; j++) {
				if (j!=i) buildNodes2Kernel_addWindow(nodes2kernel,lastNode2kernel,fromNodeIndex,kernelI,fromNode.kernels[j],windowStart,windowEnd,1,getKernelsOrientation[orientationI][fromNode.kernelOrientations[j]][0],0,0);
			}
		}
	}
	
	
	/**
	 * Adds artificial windows between $fromNode$ and every other intersecting node in the
	 * same read.
	 *
	 * Remark: the procedure assumes that the $nodeID$ field of $fromNode$ is a pointer to
	 * its position in $IntervalGraph.nodesArray$, and that the latter is sorted by read,
	 * start.
	 *
	 * Remark: the procedure uses the global variables $shortPeriodWindows,trackPointers$,
	 * which should have already been initialized.
	 */
	private static final void buildNodes2Kernel_addWindows_sameRead(IntervalGraph.Node fromNode, int fromNodeIndex, int[][] nodes2kernel, int[] lastNode2kernel, int threshold) {
		final int fromRead = fromNode.read;
		final int fromStart = fromNode.start;
		final int fromEnd = fromNode.end;
		final boolean isConcatenationFrom = (fromNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(fromRead,fromStart,fromEnd,shortPeriodWindows,trackPointers);
		final int n = IntervalGraph.nNodes;
		boolean isConcatenationTo;
		int i;
		IntervalGraph.Node toNode;
		
		for (i=fromNode.nodeID-1; i>=0; i--) {
			toNode=IntervalGraph.nodesArray[i];
			if (toNode.read!=fromRead) break;
			if (toNode.end<=fromStart+threshold) continue;
			isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(fromRead,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
			buildNodes2Kernel_addWindows_sameRead_impl(fromNode,fromNodeIndex,isConcatenationFrom,toNode,isConcatenationTo,nodes2kernel,lastNode2kernel,threshold);
		}
		for (i=fromNode.nodeID+1; i<n; i++) {
			toNode=IntervalGraph.nodesArray[i];
			if (toNode.read!=fromRead || toNode.start>=fromEnd-threshold) break;
			isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(fromRead,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
			buildNodes2Kernel_addWindows_sameRead_impl(fromNode,fromNodeIndex,isConcatenationFrom,toNode,isConcatenationTo,nodes2kernel,lastNode2kernel,threshold);
		}
	}
	
	
	/**
	 * Remark: this captures the case of e.g. Long Terminal Repeat retrotransposons, whose
	 * suffix and prefix have the same sequence.
	 *
	 * Remark: the procedure assumes that $pathKernelPeriodic[]$ has already been built.
	 *
	 * @param leftOrRight the procedure assumes that: (0) $windowStart$ is far from 
	 * $fromNode.start$; (1) $windowEnd$ is far from $fromNode.end$; (2) both are far;
	 * @return FALSE iff the window projects an end of a kernel of $toNode$, to the middle
	 * of $fromNode$: this means that $fromNode$ contains a prefix or a suffix of a kernel
	 * of $toNode$, but it is incompatible with any kernel of $fromNode$ being contained 
	 * or inserting in a kernel of $toNode$.
	 */
	private static final boolean buildNodes2kernel_useWindow(IntervalGraph.Node fromNode, int leftOrRight, IntervalGraph.Node toNode, boolean projectionStartInside, boolean projectionEndInside, int kernelFrom, int kernelTo) {
		if (pathKernelLengths[kernelFrom]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernelFrom]!=0) || pathKernelLengths[kernelTo]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernelTo]!=0)) return true;
		if (leftOrRight==0 || leftOrRight==2) {
			if (Alignments.orientation) {
				if ( ( !projectionStartInside && toNode.lastPathWithStart>=0 && 
					   ( Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,kernelTo)>=0 || 
					     Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,-1-kernelTo)>=0
					   )
				     )
				   ) return false;
			}
			else {
				if ( ( !projectionEndInside && toNode.lastPathWithEnd>=0 && 
					   ( Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,kernelTo)>=0 || 
					     Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,-1-kernelTo)>=0
					   )
				     )
				   ) return false;
			}
		}
		else if (leftOrRight==1 || leftOrRight==2) {
			if (Alignments.orientation) {
				if ( ( !projectionEndInside && toNode.lastPathWithEnd>=0 && 
	 				   ( Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,kernelTo)>=0 || 
	 				     Arrays.binarySearch(toNode.pathsWithEnd,0,toNode.lastPathWithEnd+1,-1-kernelTo)>=0
	 				   )
	 			     )
				   ) return false;
			}
			else {
				if ( ( !projectionStartInside && toNode.lastPathWithStart>=0 && 
	 				   ( Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,kernelTo)>=0 || 
	 				     Arrays.binarySearch(toNode.pathsWithStart,0,toNode.lastPathWithStart+1,-1-kernelTo)>=0
	 				   )
	 			     )
				   ) return false;
			}
		}
		return true;
	}
	
	
	/*
	 * 0=unused; 1=no insertion; 2=start of kernel; 3=end of kernel; 4=both.
	 */
	private static final byte[][] addInsertions = {
		{-1,-1,-1,-1,-1},
		{-1,1,2,3,4},
		{-1,2,2,4,4},
		{-1,3,4,3,4},
		{-1,4,4,4,4}
	};
	
	
	/*
	 * First dimension: insertion from which to remove.
	 * 0=unused; 1=no insertion; 2=start of kernel; 3=end of kernel; 4=both.
	 * Second dimension: insertion side to remove.
	 * 0=start; 1=end; 3=both.
	 */
	private static final byte[][] removeInsertion = { 
		{-1,-1,-1},
		{1,1,1},
		{1,2,1},
		{3,1,1},
		{3,2,1}
	};
	
	
	/*
	 * First/second dimension: orientation of an interval graph node WRT a kernel.
	 * 0=forward, 1=RC, 2=both.
	 * Third dimension: orientation of an alignment between two interval graph nodes.
	 * 0=forward, 1=RC.
	 * Value of a cell: orientation between the kernels induced by the alignment.
	 * 0=forward, 1=RC, 2=both.
	 */
	private static final byte[][][] getKernelsOrientation = {
		{ {0,1},{1,0},{2,2} },
		{ {1,0},{0,1},{2,2} },
	    { {2,2},{2,2},{2,2} }
	};
	
	
	/**
	 * Appends a window to $nodes2kernel[fromNode]$, without performing any compaction.
	 *
	 * @param kernelFrom,kernelTo non-negative;
	 * @param insertionDistanceStart distance between the start of $kernelFrom$ and the
	 * closest start/end of $kernelTo$ with respect to which it is an insertion.
	 */
	private static final void buildNodes2Kernel_addWindow(int[][] nodes2kernel, int[] lastNode2kernel, int fromNode, int kernelFrom, int kernelTo, int windowStart, int windowEnd, int isInsertion, int kernelsOrientation, int insertionDistanceStart, int insertionDistanceEnd) {
		final int MIN_LENGTH = BLOCK_SIZE << 2;  // Arbitrary
		final int GROWTH_RATE = BLOCK_SIZE << 8;  // Arbitrary. Should not be too small, otherwise the procedure becomes very slow.
		
		if (nodes2kernel[fromNode]==null || nodes2kernel[fromNode].length==0) {
			nodes2kernel[fromNode] = new int[MIN_LENGTH];
			lastNode2kernel[fromNode]=-1;
		}
		if (lastNode2kernel[fromNode]+BLOCK_SIZE>=nodes2kernel[fromNode].length) {
			int[] tmpArray = new int[(lastNode2kernel[fromNode]+1)+GROWTH_RATE];
			System.arraycopy(nodes2kernel[fromNode],0,tmpArray,0,lastNode2kernel[fromNode]+1);
			nodes2kernel[fromNode]=tmpArray;
		}
		nodes2kernel[fromNode][++lastNode2kernel[fromNode]]=kernelFrom;
		nodes2kernel[fromNode][++lastNode2kernel[fromNode]]=kernelTo;
		nodes2kernel[fromNode][++lastNode2kernel[fromNode]]=windowStart;
		nodes2kernel[fromNode][++lastNode2kernel[fromNode]]=windowEnd;
		nodes2kernel[fromNode][++lastNode2kernel[fromNode]]=isInsertion;
		nodes2kernel[fromNode][++lastNode2kernel[fromNode]]=kernelsOrientation;
		nodes2kernel[fromNode][++lastNode2kernel[fromNode]]=isInsertion==2||isInsertion==4?insertionDistanceStart:0;
		nodes2kernel[fromNode][++lastNode2kernel[fromNode]]=isInsertion==3||isInsertion==4?insertionDistanceEnd:0;
	}

	
	/**
	 * Merges all overlapping windows in $nodes2kernel[fromNode]$ that are related to the
	 * same pair of kernels in the same orientation (every window with two orientations is
	 * transformed into two windows with opposite orientations). Sorts 
	 * $nodes2kernel[fromNode]$ by $kernel1,kernel2,orientation,start$.
	 */
	private static final void buildNodes2Kernel_compactWindows(int[][] nodes2kernel, int[] lastNode2kernel, int fromNode, int threshold) {
		int i, j;
		int lastWindow, currentKernelFrom, currentKernelTo, currentStart, currentEnd;
		int currentInsertion, currentOrientation, currentInsertionDS, currentInsertionDE;
		
		// Sorting
		if (tmpWindows.length<((lastNode2kernel[fromNode]+1)/BLOCK_SIZE)<<1) {
			Node2KernelWindow[] newTmpWindows = new Node2KernelWindow[((lastNode2kernel[fromNode]+1)/BLOCK_SIZE)<<1];
			System.arraycopy(tmpWindows,0,newTmpWindows,0,tmpWindows.length);
			for (i=tmpWindows.length; i<newTmpWindows.length; i++) newTmpWindows[i] = new Node2KernelWindow();
			tmpWindows=newTmpWindows;
		}
		j=-1;
		for (i=0; i<=lastNode2kernel[fromNode]; i+=BLOCK_SIZE) {
			currentOrientation=nodes2kernel[fromNode][i+5];
			if (currentOrientation==2) {
				tmpWindows[++j].set(nodes2kernel[fromNode][i],nodes2kernel[fromNode][i+1],nodes2kernel[fromNode][i+2],nodes2kernel[fromNode][i+3],nodes2kernel[fromNode][i+4],0,nodes2kernel[fromNode][i+6],nodes2kernel[fromNode][i+7]);
				tmpWindows[++j].set(nodes2kernel[fromNode][i],nodes2kernel[fromNode][i+1],nodes2kernel[fromNode][i+2],nodes2kernel[fromNode][i+3],nodes2kernel[fromNode][i+4],1,nodes2kernel[fromNode][i+6],nodes2kernel[fromNode][i+7]);
			}
			else tmpWindows[++j].set(nodes2kernel[fromNode][i],nodes2kernel[fromNode][i+1],nodes2kernel[fromNode][i+2],nodes2kernel[fromNode][i+3],nodes2kernel[fromNode][i+4],currentOrientation,nodes2kernel[fromNode][i+6],nodes2kernel[fromNode][i+7]);
		}
		lastWindow=j;
		if (lastWindow>0) Arrays.sort(tmpWindows,0,lastWindow+1);
		
		// Compacting
		currentKernelFrom=tmpWindows[0].kernelFrom;
		currentKernelTo=tmpWindows[0].kernelTo;
		currentStart=tmpWindows[0].start;
		currentEnd=tmpWindows[0].end;
		currentInsertion=tmpWindows[0].isInsertion;
		currentOrientation=tmpWindows[0].orientation;
		currentInsertionDS=tmpWindows[0].insertionDistanceStart;
		currentInsertionDE=tmpWindows[0].insertionDistanceEnd;
		j=-1;
		for (i=1; i<=lastWindow; i++) {
			if (tmpWindows[i].kernelFrom!=currentKernelFrom || tmpWindows[i].kernelTo!=currentKernelTo || tmpWindows[i].orientation!=currentOrientation || tmpWindows[i].start>currentEnd+threshold) {
				tmpWindows[++j].set(currentKernelFrom,currentKernelTo,currentStart,currentEnd,currentInsertion,currentOrientation,currentInsertionDS,currentInsertionDE);
				currentKernelFrom=tmpWindows[i].kernelFrom;
				currentKernelTo=tmpWindows[i].kernelTo;
				currentStart=tmpWindows[i].start;
				currentEnd=tmpWindows[i].end;
				currentInsertion=tmpWindows[i].isInsertion;
				currentOrientation=tmpWindows[i].orientation;
				currentInsertionDS=tmpWindows[i].insertionDistanceStart;
				currentInsertionDE=tmpWindows[i].insertionDistanceEnd;
				continue;
			}
			currentEnd=Math.max(currentEnd,tmpWindows[i].end);
			currentInsertion=addInsertions[currentInsertion][tmpWindows[i].isInsertion];
			currentInsertionDS=Math.max(currentInsertionDS,tmpWindows[i].insertionDistanceStart);
			currentInsertionDE=Math.max(currentInsertionDE,tmpWindows[i].insertionDistanceEnd);
		}
		tmpWindows[++j].set(currentKernelFrom,currentKernelTo,currentStart,currentEnd,currentInsertion,currentOrientation,currentInsertionDS,currentInsertionDE);
		lastWindow=j;
		
		// Updating $nodes2kernel[fromNode]$.
		if (nodes2kernel[fromNode].length<(lastWindow+1)*BLOCK_SIZE) nodes2kernel[fromNode] = new int[(lastWindow+1)*BLOCK_SIZE];
		j=0;
		for (i=0; i<=lastWindow; i++) {
			nodes2kernel[fromNode][j++]=tmpWindows[i].kernelFrom;
			nodes2kernel[fromNode][j++]=tmpWindows[i].kernelTo;
			nodes2kernel[fromNode][j++]=tmpWindows[i].start;
			nodes2kernel[fromNode][j++]=tmpWindows[i].end;
			nodes2kernel[fromNode][j++]=tmpWindows[i].isInsertion;
			nodes2kernel[fromNode][j++]=tmpWindows[i].orientation;
			nodes2kernel[fromNode][j++]=tmpWindows[i].insertionDistanceStart;
			nodes2kernel[fromNode][j++]=tmpWindows[i].insertionDistanceEnd;
		}
		lastNode2kernel[fromNode]=j-1;
	}
	
	
	private static class Node2KernelWindow implements Comparable {
		public int kernelFrom, kernelTo, start, end, isInsertion, orientation;
		public int insertionDistanceStart, insertionDistanceEnd;
		
		public Node2KernelWindow() {
			this.kernelFrom=-1;
			this.kernelTo=-1;
			this.start=-1;
			this.end=-1;
			isInsertion=-1;
			orientation=-1;
			insertionDistanceStart=0;
			insertionDistanceEnd=0;
		}
		
		public Node2KernelWindow(int kf, int kt, int s, int e, int i, int o, int si, int ei) {
			set(kf,kt,s,e,i,o,si,ei);
		}
		
		public void set(int kf, int kt, int s, int e, int i, int o, int si, int ei) {
			kernelFrom=kf;
			kernelTo=kt;
			start=s;
			end=e;
			isInsertion=i;
			orientation=o;
			insertionDistanceStart=si;
			insertionDistanceEnd=ei;
		}
		
		public int compareTo(Object other) {
			Node2KernelWindow otherWindow = (Node2KernelWindow)other;
			if (kernelFrom<otherWindow.kernelFrom) return -1;
			else if (kernelFrom>otherWindow.kernelFrom) return 1;
			if (kernelTo<otherWindow.kernelTo) return -1;
			else if (kernelTo>otherWindow.kernelTo) return 1;
			if (orientation<otherWindow.orientation) return -1;
			else if (orientation>otherWindow.orientation) return 1;
			if (start<otherWindow.start) return -1;
			else if (start>otherWindow.start) return 1;
			return 0;
		}	
	}
	
	
	/**
	 * Builds the kernel graph from the $nodes2kernel$ matrix. $inNeighbors[i]$ is sorted,
	 * but $outNeighbors[i]$ is not necessarily sorted. The procedure computes also 
	 * $inNeighborFlags[i]$ (see its declaration).
	 *
	 * Remark: the procedure assumes $nodesInKernel[from..to]$ to be sorted by 
	 * $bidirectedGraphNode$.
	 *
	 * @param from,to interval in $nodesInKernel$.
	 */
	private static final void buildKernelGraph_impl(int[][] nodes2kernel, int[] lastNode2kernel, IntervalGraph.Node[] nodes, int nNodes, int from, int to, int kernel) {
		final int GROWTH_RATE = 10;  // Arbitrary
		int i, j, k, p, q, r;
		int last1, last2, last4, value, source, bidirected, previous, previousBidirected, kernelTo, previousKernelTo;
		int oldInsertion, newInsertion, orientation, out;
		final int previousOrder = IntervalGraph.Node.order;
		int[] tmpArray;
		IntervalGraph.Node tmpNode;

		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		
		// Building $inNeighbors$ and $outNeighbors$.
		last2=-1;
		for (i=from; i<=to; i++) {
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[i]];
			j=Arrays.binarySearch(nodes,0,nNodes,tmpNode);
			if (lastNode2kernel[j]==-1) continue;
			last1=-1; previousKernelTo=-1;
			for (p=0; p<=lastNode2kernel[j]; p+=BLOCK_SIZE_PRIME) {
				if (nodes2kernel[j][p]>kernel) break;
				if (nodes2kernel[j][p]<kernel) continue;
				kernelTo=nodes2kernel[j][p+1];
				if (kernelTo!=previousKernelTo) {
					tmpArray1[++last1]=kernelTo;
					previousKernelTo=kernelTo;
				}
			}
			if (last2==-1) {
				tmpArray=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpArray;
				last2=last1;
			}
			else {
				last2=Math.setUnion(tmpArray1,last1,tmpArray2,last2,tmpArray3);
				tmpArray=tmpArray2; tmpArray2=tmpArray3; tmpArray3=tmpArray;
			}
		}
		lastInNeighbor[kernel]=last2;
		if (last2==-1) return;
		last2=buildKernelGraph_impl_filter(kernel,tmpArray2,last2);
		lastInNeighbor[kernel]=last2;
		if (last2==-1) return;
		if (inNeighbors[kernel]==null || inNeighbors[kernel].length<last2+1) inNeighbors[kernel] = new int[last2+1];
		System.arraycopy(tmpArray2,0,inNeighbors[kernel],0,last2+1);
		for (i=0; i<=last2; i++) {
			source=tmpArray2[i];
			lastOutNeighbor[source]++;
			if (outNeighbors[source]==null || outNeighbors[source].length==0) outNeighbors[source] = new int[GROWTH_RATE];
			else if (lastOutNeighbor[source]==outNeighbors[source].length) {
				tmpArray = new int[outNeighbors[source].length+GROWTH_RATE];
				System.arraycopy(outNeighbors[source],0,tmpArray,0,outNeighbors[source].length);
				outNeighbors[source]=tmpArray;
			}
			outNeighbors[source][lastOutNeighbor[source]]=kernel;
		}		
		
		// Initializing $inNeighborFlags$, for every arc and for every orientation (even
	 	// if an orientation is never observed).
		if (inNeighborFlags[kernel]==null || inNeighborFlags[kernel].length<(lastInNeighbor[kernel]+1)<<2) inNeighborFlags[kernel] = new byte[(lastInNeighbor[kernel]+1)<<2];
		for (i=0; i<=lastInNeighbor[kernel]; i++) {
			j=i<<2;
			inNeighborFlags[kernel][j+0]=1;
			inNeighborFlags[kernel][j+1]=0;
			inNeighborFlags[kernel][j+2]=1;
			inNeighborFlags[kernel][j+3]=0;
		}
		
		// Building $inNeighborFlags$: insertion flags.
		if (tmpArray1==null || tmpArray1.length<(lastInNeighbor[kernel]+1)<<1) tmpArray1 = new int[(lastInNeighbor[kernel]+1)<<1];
		if (tmpArray2==null || tmpArray2.length<(lastInNeighbor[kernel]+1)<<1) tmpArray2 = new int[(lastInNeighbor[kernel]+1)<<1];
		Math.set(tmpArray1,((lastInNeighbor[kernel]+1)<<1)-1,0);  // Insertion distance start for each neighbor and orientation
		Math.set(tmpArray2,((lastInNeighbor[kernel]+1)<<1)-1,0);  // Insertion distance end for each neighbor and orientation
		for (i=from; i<=to; i++) {
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[i]];
			j=Arrays.binarySearch(nodes,0,nNodes,tmpNode);
			previousKernelTo=-1; q=-1;
			for (p=0; p<=lastNode2kernel[j]; p+=BLOCK_SIZE_PRIME) {
				if (nodes2kernel[j][p]>kernel) break;
				if (nodes2kernel[j][p]<kernel) continue;
				kernelTo=nodes2kernel[j][p+1];
				if (kernelTo!=previousKernelTo) q=Arrays.binarySearch(inNeighbors[kernel],0,lastInNeighbor[kernel]+1,kernelTo);
				else { /* $q$ holds already the correct value */ }
				if (q<0) {
					// Could have been removed by $buildKernelGraph_impl_filter()$.
					continue;
				}
				orientation=nodes2kernel[j][p+3];
				r=(q<<1)+(orientation==0?0:1);
				if (tmpArray1[r]!=maxReadLength) tmpArray1[r]=Math.max(tmpArray1[r],nodes2kernel[j][p+4]);
				if (tmpArray2[r]!=maxReadLength) tmpArray2[r]=Math.max(tmpArray2[r],nodes2kernel[j][p+5]);
			}
		}
		for (i=0; i<=lastInNeighbor[kernel]; i++) {
			if (tmpArray1[(i<<1)+0]>=INSERTION_DISTANCE_THRESHOLD) inNeighborFlags[kernel][(i<<2)+0]=addInsertions[inNeighborFlags[kernel][(i<<2)+0]][2];
			if (tmpArray1[(i<<1)+1]>=INSERTION_DISTANCE_THRESHOLD) inNeighborFlags[kernel][(i<<2)+2]=addInsertions[inNeighborFlags[kernel][(i<<2)+2]][2];
		}
		for (i=0; i<=lastInNeighbor[kernel]; i++) {
			if (tmpArray2[(i<<1)+0]>=INSERTION_DISTANCE_THRESHOLD) inNeighborFlags[kernel][(i<<2)+0]=addInsertions[inNeighborFlags[kernel][(i<<2)+0]][3];
			if (tmpArray2[(i<<1)+1]>=INSERTION_DISTANCE_THRESHOLD) inNeighborFlags[kernel][(i<<2)+2]=addInsertions[inNeighborFlags[kernel][(i<<2)+2]][3];
		}
		IntervalGraph.Node.order=previousOrder;	
		
		// Building $inNeighborFlags$: all-kernel-nodes flags.
		// If orientation $X$ is never observed for an arc, its $a_X$ value is zero, and
		// the following steps will have to avoid orientations of arcs with $a_X=0$.
		buildKernelGraph_impl_allFlags(0,nodes2kernel,lastNode2kernel,nodes,nNodes,from,to,kernel);
		buildKernelGraph_impl_allFlags(1,nodes2kernel,lastNode2kernel,nodes,nNodes,from,to,kernel);
	}
	
	
	/**
	 * Sets the $a_X$ flags of $inNeighborFlags[kernel]$, where $X=orientation$. $a_X$ is 
	 * set to 1 for all and only the kernels in $U_0 \intersect U_1 \intersect ... U_k$, 
	 * where $U_i$ is the union of the $nodes2kernel$ lists of all nodes in 
	 * $nodesInKernel[from..to]$ that are mapped to the same node $i$ of the bidirected 
	 * graph of which the kernel is a path.
	 *
	 * Remark: the procedure assumes $nodesInKernel[from..to]$ to be sorted by 
	 * $bidirectedGraphNode$.
	 */
	private static final void buildKernelGraph_impl_allFlags(int orientation, int[][] nodes2kernel, int[] lastNode2kernel, IntervalGraph.Node[] nodes, int nNodes, int from, int to, int kernel) {
		int i, j, p;
		int last1, last2, last4, bidirected, previousBidirected, kernelTo;
		final int previousOrder = IntervalGraph.Node.order;
		int[] tmpArray;
		IntervalGraph.Node tmpNode;
		
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		
int fabio = 271;		
if (IO.SHOW_INTERACTIVE && kernel==fabio) System.err.println("buildKernelGraph_impl_allFlags> kernel="+kernel+" orientation="+orientation);
		
		// Computing intersection of unions
		previousBidirected=-1; last2=-1; last4=-2;
		for (i=from; i<=to; i++) {
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[i]];			
			bidirected=tmpNode.bidirectedGraphNode;
if (IO.SHOW_INTERACTIVE && kernel==fabio) System.err.println("buildKernelGraph_impl_allFlags> bidirected="+bidirected+" previousBidirected="+previousBidirected+"  tmpNode="+tmpNode);
			if (bidirected!=previousBidirected) {
				if (previousBidirected!=-1) {
					
if (IO.SHOW_INTERACTIVE && kernel==fabio) {
	System.err.print("buildKernelGraph_impl_allFlags> union: ");
	for (int x=0; x<=last2; x++) System.err.print(tmpArray2[x]+",");
	System.err.println();
	System.err.print("buildKernelGraph_impl_allFlags> intersection: ");
	for (int x=0; x<=last4; x++) System.err.print(tmpArray4[x]+",");
	System.err.println();
}
					
					
					if (last2==-1) {
						last4=-1;
						break;
					}
					else {
						if (last4==-2) {
							System.arraycopy(tmpArray2,0,tmpArray4,0,last2+1);
							last4=last2;
						}
						else {
							last4=Math.setIntersection(tmpArray2,0,last2,tmpArray4,0,last4,tmpArray3,0);
							if (last4==-1) break;
							tmpArray=tmpArray4; tmpArray4=tmpArray3; tmpArray3=tmpArray;
						}
					}
				}
				previousBidirected=bidirected;
				last2=-1;
				j=Arrays.binarySearch(nodes,0,nNodes,tmpNode);
				if (lastNode2kernel[j]==-1) continue;
				for (p=0; p<=lastNode2kernel[j]; p+=BLOCK_SIZE_PRIME) {
					if (nodes2kernel[j][p]>kernel) break;
					if (nodes2kernel[j][p]<kernel || nodes2kernel[j][p+3]!=orientation) continue;
					tmpArray2[++last2]=nodes2kernel[j][p+1];
				}
			}
			else {
				last1=-1;
				j=Arrays.binarySearch(nodes,0,nNodes,tmpNode);
				if (lastNode2kernel[j]==-1) continue;
				for (p=0; p<=lastNode2kernel[j]; p+=BLOCK_SIZE_PRIME) {
					if (nodes2kernel[j][p]>kernel) break;
					if (nodes2kernel[j][p]<kernel || nodes2kernel[j][p+3]!=orientation) continue;
					tmpArray1[++last1]=nodes2kernel[j][p+1];
				}
				if (last2==-1) {
					tmpArray=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpArray;
					last2=last1;
				}
				else {
					last2=Math.setUnion(tmpArray1,last1,tmpArray2,last2,tmpArray3);
					tmpArray=tmpArray2; tmpArray2=tmpArray3; tmpArray3=tmpArray;
				}
			}
		}
		if (last4!=-1) {
			
if (IO.SHOW_INTERACTIVE && kernel==fabio) {
	System.err.print("buildKernelGraph_impl_allFlags> union: ");
	for (int x=0; x<=last2; x++) System.err.print(tmpArray2[x]+",");
	System.err.println();
	System.err.print("buildKernelGraph_impl_allFlags> intersection: ");
	for (int x=0; x<=last4; x++) System.err.print(tmpArray4[x]+",");
	System.err.println();
}
			
			
			if (last2==-1) last4=-1;
			else {
				if (last4==-2) {
					System.arraycopy(tmpArray2,0,tmpArray4,0,last2+1);
					last4=last2;
				}
				else {
					last4=Math.setIntersection(tmpArray2,0,last2,tmpArray4,0,last4,tmpArray3,0);
					if (last4!=-1) {
						tmpArray=tmpArray4; tmpArray4=tmpArray3; tmpArray3=tmpArray;
					}
				}
			}			
		}
if (IO.SHOW_INTERACTIVE && kernel==fabio) System.err.println("buildKernelGraph_impl_allFlags> last4="+last4+" lastInNeighbor="+lastInNeighbor[kernel]);	
if (IO.SHOW_INTERACTIVE && kernel==fabio) {
	System.err.print("buildKernelGraph_impl_allFlags> final result: ");
	for (int x=0; x<=last4; x++) System.err.print(tmpArray4[x]+",");
	System.err.println();
}
		
		// Updating $inNeighborFlags$.
		// No need to call $buildKernelGraph_impl_filter()$, since we are intersecting
		// with $inNeighbors[kernel]$ which is already filtered.
		i=0; j=0;
		while (i<=last4 && j<=lastInNeighbor[kernel]) {
			if (tmpArray4[i]<inNeighbors[kernel][j]) i++;
			else if (tmpArray4[i]>inNeighbors[kernel][j]) j++;
			else {
				kernelTo=tmpArray4[i];
				p=Arrays.binarySearch(inNeighbors[kernel],0,lastInNeighbor[kernel]+1,kernelTo);
				inNeighborFlags[kernel][(p<<2)+(orientation==0?1:3)]=1;
				i++; j++;
			}
		}
		
		IntervalGraph.Node.order=previousOrder;
	}
	
	
	/**
	 * Removes from $fromKernels$ all kernels whose string length does not make it safe to 
	 * assume that they can contain $toKernel$.
	 *
	 * Remark: using a threshold rather than strict equality/inequality would not 
	 * guarantee the kernel graph to be acyclic.
	 *
	 * @return the last element of $fromKernels$ after removal.
	 */
	private static final int buildKernelGraph_impl_filter(int toKernel, int[] fromKernels, int lastFromKernel) {
		int i, j;
		int fromKernel, fromLength;
		final int toLength = pathKernelLengths[toKernel]==0?pathKernelLongestNode[toKernel]:pathKernelLengths[toKernel];
		
		j=-1;
		for (i=0; i<=lastFromKernel; i++) {
			fromKernel=fromKernels[i];
			fromLength=pathKernelLengths[fromKernel]==0?pathKernelLongestNode[fromKernel]:pathKernelLengths[fromKernel];
			if (fromLength>toLength || (fromLength==toLength && fromKernel>toKernel)) {
				j++;
				if (j!=i) fromKernels[j]=fromKernel;
			}
		}
		return j;
	}
	
	
	/**
	 * Adds to $oldFlags[from..from+3]$ the $inNeighborFlags$ array in $newFlags[0..3]$.
	 */
	private static final void addInNeighborFlags(byte[] oldFlags, int from, byte[] newFlags) {
		if (newFlags[1]==1) {
			if (oldFlags[from+1]==0) {
				oldFlags[from+1]=1;
				oldFlags[from]=newFlags[0];
			}
			else oldFlags[from]=addInsertions[oldFlags[from]][newFlags[0]];
		}
		if (newFlags[3]==1) {
			if (oldFlags[from+3]==0) {
				oldFlags[from+3]=1;
				oldFlags[from+2]=newFlags[2];
			}
			else oldFlags[from+2]=addInsertions[oldFlags[from+2]][newFlags[2]];
		}
	}
	
	
	/**
	 * 0=atStart, 1=atEnd.
	 */
	private static final byte[][] getFlag = new byte[][] {
		{4,2}, {3,1}
	};
	
	
	/**
	 * First dimension: flag \in {1,2,3,4} (0 is unused).
	 * Second dimension: orientation \in {1,2,3} (0 is unused).
	 * Cell: array of 4 elements.
	 */
	private static final byte[][][] getInNeighborFlags = new byte[][][] {
		{ {1,0,1,0}, {1,0,1,0}, {1,0,1,0}, {1,0,1,0} },
		{ {1,0,1,0}, {1,1,1,0}, {1,0,1,1}, {1,1,1,1} },
		{ {1,0,1,0}, {2,1,1,0}, {1,0,2,1}, {2,1,2,1} },
		{ {1,0,1,0}, {3,1,1,0}, {1,0,3,1}, {3,1,3,1} },
		{ {1,0,1,0}, {4,1,1,0}, {1,0,4,1}, {4,1,4,1} }
	};

	
	/**
	 * Multiple full copies of the same repeat might occur in the same read: since 
	 * factorization does not use alignments of a read to itself, there is no alignment
	 * between those copies, and they might thus become the nodes of two distinct path 
	 * kernels. If those kernels have a large enough basin, they might still be connected 
	 * by an arc in the kernel graph. The procedure tries to handle the case in which the 
	 * two kernels are not reachable in the kernel graph. More generally, two distinct
	 * path kernels that correspond to the same repeat might not have an arc in the kernel
	 * graph, because at least one pair of their nodes lie in the same read. The procedure
	 * creates identity arcs between every pair of path kernels that: (1) are not already 
	 * reachable in the kernel graph; (2) have at least one pair of kernel nodes that 
	 * occur in the same read without overlapping; (3) have similar string lengths; 
	 * (4) their basins have large Jaccard similarity.
	 *
	 * Remark: the procedure assumes that global variables $descendants,lastDescendant$
	 * have already been initialized (see $getFrequency_buildDescendantsList()$).
	 *
	 * Remark: the procedure assumes that every row of $inNeighbors$ is sorted, and it is 
	 * left sorted by the procedure. A row of $outNeighbors$ is not necessarily sorted.
	 *
	 * @return number of arcs created by the procedure.
	 */
	private static final int buildKernelGraph_addSameReadArcs_prime() {
		final int LENGTH_THRESHOLD = IO.quantum<<1;  // Arbitrary
		final int OVERLAP_THRESHOLD = IO.quantum;  // Arbitrary
		final int GROWTH_RATE = 20;  // Arbitrary (multiple of 5).
		final double JACCARD_THRESHOLD = 0.9;  // Arbitrary
		boolean kernelIPeriodic, kernelJPeriodic;
		int i, j, k, p, q;
		int previous, last, kernel, lastKernelPair, lengthI, lengthJ, kernelI, kernelJ;
		int kernelFrom, kernelTo;
		int start1, end1, start2, end2, out;
		IntervalGraph.Node tmpNode;
		int[] newArray, lastKernelRead;
		int[][] kernelReads;
		SameReadPair[] pairs;
		
		// Building, for each kernel, the set of reads used by its nodes in
		// $nodesInKernel$.
		kernelReads = new int[nPathKernels][0];
		lastKernelRead = new int[nPathKernels];
		Math.set(lastKernelRead,nPathKernels-1,-1);
		kernel=blockKernels[0]; 
		tmpArray1[0]=IntervalGraph.nodesArray[nodesInKernel[0]].read;
		last=0;
		for (i=1; i<nNodesInKernel; i++) {
			if (blockKernels[i]!=kernel) {
				if (last>0) Arrays.sort(tmpArray1,0,last+1);
				j=0; previous=tmpArray1[0];
				for (k=1; k<=last; k++) {
					if (tmpArray1[k]==previous) continue;
					tmpArray1[++j]=tmpArray1[k];
					previous=tmpArray1[k];
				}
				lastKernelRead[kernel]=j;
				kernelReads[kernel]=Arrays.copyOf(tmpArray1,j+1);
				kernel=blockKernels[i];
				tmpArray1[0]=IntervalGraph.nodesArray[nodesInKernel[i]].read;
				last=0;
			}
			else {
				last++;
				if (last==tmpArray1.length) {
					newArray = new int[tmpArray1.length+GROWTH_RATE];
					System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
					tmpArray1=newArray;
				}
				tmpArray1[last]=IntervalGraph.nodesArray[nodesInKernel[i]].read;
			}
		}
		if (last>0) Arrays.sort(tmpArray1,0,last+1);
		j=0; previous=tmpArray1[0];
		for (k=1; k<=last; k++) {
			if (tmpArray1[k]==previous) continue;
			tmpArray1[++j]=tmpArray1[k];
			previous=tmpArray1[k];
		}
		lastKernelRead[kernel]=j;
		kernelReads[kernel]=Arrays.copyOf(tmpArray1,j+1);
		
		// Computing the pairs of kernels for which we should evaluate the Jaccard
		// similarity of their basins.
		pairs = new SameReadPair[nPathKernels];
		last=-1;
		for (i=0; i<nPathKernels; i++) pairs[++last] = new SameReadPair(i,pathKernelLengths[i]);
		if (last>0) Arrays.sort(pairs,0,last+1);
		lastKernelPair=-1;
		for (i=0; i<=last; i++) {
			lengthI=pairs[i].length;
			kernelI=pairs[i].kernelID;
			kernelIPeriodic=pathKernelLengths[kernelI]==0||(pathKernelPeriodic!=null&&pathKernelPeriodic[kernelI]!=0);
			for (j=i+1; j<=last; j++) {
				lengthJ=pairs[j].length;
				kernelJ=pairs[j].kernelID;
				kernelJPeriodic=pathKernelLengths[kernelJ]==0||(pathKernelPeriodic!=null&&pathKernelPeriodic[kernelJ]!=0);
				if (kernelIPeriodic!=kernelJPeriodic) continue;
				if (!kernelJPeriodic && lengthJ>lengthI+LENGTH_THRESHOLD) break;
				if ( !Math.nonemptyIntersection(kernelReads[kernelI],0,lastKernelRead[kernelI],kernelReads[kernelJ],0,lastKernelRead[kernelJ]) ||
					 Arrays.binarySearch(descendants[kernelI],0,lastDescendant[kernelI]+1,kernelJ)>=0 || 
					 Arrays.binarySearch(descendants[kernelJ],0,lastDescendant[kernelJ]+1,kernelI)>=0 ||
					 !haveNonoverlappingPair(kernelI,kernelJ,OVERLAP_THRESHOLD)
				   ) continue;
				if (lastKernelPair+5>=tmpArray1.length) {
					newArray = new int[tmpArray1.length+GROWTH_RATE];
					System.arraycopy(tmpArray1,0,newArray,0,tmpArray1.length);
					tmpArray1=newArray;
				}
				tmpArray1[++lastKernelPair]=kernelI;
				tmpArray1[++lastKernelPair]=kernelJ;
				tmpArray1[++lastKernelPair]=0;  // Union
				tmpArray1[++lastKernelPair]=0;  // Intersection
				tmpArray1[++lastKernelPair]=0;  // Orientation
			}
		}
		pairs=null;
		System.err.println("buildKernelGraph_addSameReadArcs_prime> "+((lastKernelPair+1)/5)+" candidate kernel pairs");
		
		// Computing the Jaccard similarity of basins
		for (i=0; i<IntervalGraph.nNodes; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			if (tmpNode.lastKernel==-1) continue;
			for (j=0; j<=lastKernelPair; j+=5) {
				p=Arrays.binarySearch(tmpNode.kernels,0,tmpNode.lastKernel+1,tmpArray1[j]);
				q=Arrays.binarySearch(tmpNode.kernels,0,tmpNode.lastKernel+1,tmpArray1[j+1]);
				if (p>=0 || q>=0) tmpArray1[j+2]++;
				if (p>=0 && q>=0) tmpArray1[j+3]++;
				if (tmpNode.lastPathWithStart>=0) {
					start1=Arrays.binarySearch(tmpNode.pathsWithStart,0,tmpNode.lastPathWithStart+1,tmpArray1[j]);
					end1=Arrays.binarySearch(tmpNode.pathsWithStart,0,tmpNode.lastPathWithStart+1,-1-tmpArray1[j]);
					if (start1>=0 || end1>=0) {
						start2=Arrays.binarySearch(tmpNode.pathsWithStart,0,tmpNode.lastPathWithStart+1,tmpArray1[j+1]);
						end2=Arrays.binarySearch(tmpNode.pathsWithStart,0,tmpNode.lastPathWithStart+1,-1-tmpArray1[j+1]);
						if ((start1>=0 && start2>=0) || (end1>=0 && end2>=0)) tmpArray1[j+4]=addOrientation_prime[tmpArray1[j+4]][1];
						if ((start1>=0 && end2>=0) || (end1>=0 && start2>=0)) tmpArray1[j+4]=addOrientation_prime[tmpArray1[j+4]][2];
					}
				}
				if (tmpNode.lastPathWithEnd>=0) {
					start1=Arrays.binarySearch(tmpNode.pathsWithEnd,0,tmpNode.lastPathWithEnd+1,tmpArray1[j]);
					end1=Arrays.binarySearch(tmpNode.pathsWithEnd,0,tmpNode.lastPathWithEnd+1,-1-tmpArray1[j]);
					if (start1>=0 || end1>=0) {
						start2=Arrays.binarySearch(tmpNode.pathsWithEnd,0,tmpNode.lastPathWithEnd+1,tmpArray1[j+1]);
						end2=Arrays.binarySearch(tmpNode.pathsWithEnd,0,tmpNode.lastPathWithEnd+1,-1-tmpArray1[j+1]);
						if ((start1>=0 && start2>=0) || (end1>=0 && end2>=0)) tmpArray1[j+4]=addOrientation_prime[tmpArray1[j+4]][1];
						if ((start1>=0 && end2>=0) || (end1>=0 && start2>=0)) tmpArray1[j+4]=addOrientation_prime[tmpArray1[j+4]][2];
					}
				}
			}
		}
		
		// Adding arcs
		out=0;
		for (i=0; i<=lastKernelPair; i+=5) {
			if (IO.SHOW_INTERACTIVE) System.err.println("buildKernelGraph_addSameReadArcs_prime> candidate pair "+tmpArray1[i]+","+tmpArray1[i+1]+" jaccard="+(((double)tmpArray1[i+3])/tmpArray1[i+2])+" union="+tmpArray1[i+2]+" intersection="+tmpArray1[i+3]+" orientation="+tmpArray1[i+4]);
			if (tmpArray1[i+3]<tmpArray1[i+2]*JACCARD_THRESHOLD) continue;
			kernelI=tmpArray1[i];
			lengthI=pathKernelLengths[kernelI]!=0?pathKernelLengths[kernelI]:pathKernelLongestNode[kernelI];
			kernelJ=tmpArray1[i+1];
			lengthJ=pathKernelLengths[kernelJ]!=0?pathKernelLengths[kernelJ]:pathKernelLongestNode[kernelJ];
			if (lengthI>lengthJ) { kernelFrom=kernelI; kernelTo=kernelJ; }
			else if (lengthI==lengthJ) {
				if (kernelI>kernelJ) { kernelFrom=kernelI; kernelTo=kernelJ; }
				else { kernelFrom=kernelJ; kernelTo=kernelI; }
			}
			else { kernelFrom=kernelJ; kernelTo=kernelI; }
			buildKernelGraph_addSameReadArcs_impl(kernelFrom,0/*Fake start*/,0/*Fake end*/,kernelTo,0/*Fake start*/,0/*Fake end*/,tmpArray1[i+4]);
			out++;
		}
		return out;
	}
	
	
	/**
	 * 0=unknown, 1=forward, 2=reverse, 3=both.
	 */
	private static final byte[][] addOrientation_prime = new byte[][] {   
		{0,1,2,3}, {1,1,3,3}, {2,3,2,3}, {3,3,3,3}
	};
	
	
	private static class SameReadPair implements Comparable {
		
		public int kernelID, length;
		
		public SameReadPair(int k, int l) {
			kernelID=k;
			length=l;
		}
		
		public int compareTo(Object other) {
			SameReadPair otherPair = (SameReadPair)other;
			if (length<otherPair.length) return -1;
			else if (length>otherPair.length) return 1;
			else return 0;
		}
	}
	
	
	/**
	 * @param orientation 1=forward, 2=reverse, 3=both;
	 * @return TRUE iff the procedure added an arc.
	 */
	private static final boolean buildKernelGraph_addSameReadArcs_impl(int kernelFrom, int kernelFromStart, int kernelFromEnd, int kernelTo, int kernelToStart, int kernelToEnd, int orientation) {
		final int GROWTH_RATE = 10;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		final byte flag = getFlag[(Math.abs(kernelFromStart,kernelToStart)<=IDENTITY_THRESHOLD)?1:0][(Math.abs(kernelFromEnd,kernelToEnd)<=IDENTITY_THRESHOLD)?1:0];
		final byte[] newFlags = getInNeighborFlags[flag][orientation];
		byte tmpB;
		int i, j, p, q;
		int tmpI, last;
		byte[] tmpByte;
		int[] tmpInt;
	
		// The arc might already exist
		i=Arrays.binarySearch(inNeighbors[kernelTo],0,lastInNeighbor[kernelTo]+1,kernelFrom);
		if (i>=0) {
			addInNeighborFlags(inNeighborFlags[kernelTo],i<<2,newFlags);
			return false;
		}
	
		// Adding to $inNeighbors$.
		last=++lastInNeighbor[kernelTo];
		if (inNeighbors[kernelTo]==null || inNeighbors[kernelTo].length==last) {
			tmpInt = new int[last+GROWTH_RATE];
			System.arraycopy(inNeighbors[kernelTo],0,tmpInt,0,last);
			inNeighbors[kernelTo]=tmpInt;
		}
		if (inNeighborFlags[kernelTo]==null || inNeighborFlags[kernelTo].length<(last+1)<<2) {
			tmpByte = new byte[(last+GROWTH_RATE)<<2];
			System.arraycopy(inNeighborFlags[kernelTo],0,tmpByte,0,last<<2);
			inNeighborFlags[kernelTo]=tmpByte;
		}
		inNeighbors[kernelTo][last]=kernelFrom;
		System.arraycopy(newFlags,0,inNeighborFlags[kernelTo],last<<2,4);
		i=last;
		while (i>0 && inNeighbors[kernelTo][i]<inNeighbors[kernelTo][i-1]) {
			tmpI=inNeighbors[kernelTo][i-1];
			inNeighbors[kernelTo][i-1]=inNeighbors[kernelTo][i];
			inNeighbors[kernelTo][i]=tmpI;
			p=(i-1)<<2; q=i<<2;
			for (j=0; j<=3; j++) {
				tmpB=inNeighborFlags[kernelTo][p+j];
				inNeighborFlags[kernelTo][p+j]=inNeighborFlags[kernelTo][q+j];
				inNeighborFlags[kernelTo][q+j]=tmpB;
			}
			i--;
		}
	
		// Adding to $outNeighbors$.
		last=++lastOutNeighbor[kernelFrom];
		if (outNeighbors[kernelFrom]==null || outNeighbors[kernelFrom].length==last) {
			tmpInt = new int[lastOutNeighbor[kernelFrom]+GROWTH_RATE];
			System.arraycopy(outNeighbors[kernelFrom],0,tmpInt,0,last);
			outNeighbors[kernelFrom]=tmpInt;
		}
		outNeighbors[kernelFrom][last]=kernelTo;
	
		return true;
	}
	
	
	/**
	 * Consider two path kernels whose references consist of just one interval graph node 
	 * each, and both nodes straddle in the same read (but are not contained nor they
	 * coincide); assume also that there is an interval graph arc between the two path 
	 * kernels. In practice it is more likely that the two path kernels are a permutation 
	 * of one another, rather than being contained in one another, so the procedure 
	 * conservatively deletes the arc.
	 *
	 * Remark: $buildKernelGraph_removePermutations()$ cannot detect if the two path 
	 * kernels are a permutation of one another, since there are no alignments between 
	 * their reference nodes.
	 *
	 * @return number of arcs of the kernel graph removed by the procedure.
	 */
	private static final int buildKernelGraph_removeSameReadArcs() {
		boolean found;
		int i, j, k;
		int firstI, firstJ, kernelI, kernelJ, readI, readJ;
		int startI, startJ, endI, endJ, out;
		int from, fromStart, fromEnd, to, toStart, toEnd;
		IntervalGraph.Node tmpNode;
		
		out=0; firstI=0; kernelI=blockKernels[0];
		for (i=1; i<nNodesInKernel; i++) {
			if (blockKernels[i]==kernelI) continue;
			if (i-firstI>1 || pathKernelLengths[kernelI]==0) {
				firstI=i; kernelI=blockKernels[i];
				continue;
			}
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[firstI]];
			readI=tmpNode.read; startI=tmpNode.start; endI=tmpNode.end;
			firstJ=i; kernelJ=blockKernels[i];
			for (j=i+1; j<nNodesInKernel; j++) {
				if (blockKernels[j]==kernelJ) continue;
				if (j-firstJ>1 || pathKernelLengths[kernelJ]==0) {
					firstJ=j; kernelJ=blockKernels[j];
					continue;
				}
				tmpNode=IntervalGraph.nodesArray[nodesInKernel[firstJ]];
				if (tmpNode.read==readI && Intervals.straddles(startI,endI,tmpNode.start,tmpNode.end)) {
					found=false;
					for (k=0; k<=lastInNeighbor[kernelI]; k++) {
						if (inNeighbors[kernelI][k]==kernelJ) {
							removeArc(kernelI,k);
							out++; found=true;
							break;
						}
					}
					if (!found) {
						for (k=0; k<=lastInNeighbor[kernelJ]; k++) {
							if (inNeighbors[kernelJ][k]==kernelI) {
								removeArc(kernelJ,k);
								out++;
								break;
							}
						}
					}
				}
				firstJ=j; kernelJ=blockKernels[j];
			}
			if (firstJ==nNodesInKernel-1 && pathKernelLengths[kernelJ]!=0) {
				tmpNode=IntervalGraph.nodesArray[nodesInKernel[firstJ]];
				if (tmpNode.read==readI && Intervals.straddles(startI,endI,tmpNode.start,tmpNode.end)) {
					found=false;
					for (k=0; k<=lastInNeighbor[kernelI]; k++) {
						if (inNeighbors[kernelI][k]==kernelJ) {
							removeArc(kernelI,k);
							out++; found=true;
							break;
						}
					}
					if (!found) {
						for (k=0; k<=lastInNeighbor[kernelJ]; k++) {
							if (inNeighbors[kernelJ][k]==kernelI) {
								removeArc(kernelJ,k);
								out++;
								break;
							}
						}
					}
				}
			}
			firstI=i; kernelI=blockKernels[i];
		}
		return out;
	}
	
	
	/**
	 * Stores in $nKernelNodes[i]$ the first position of the block of path kernel $i$ in 
	 * $nodesInKernel$.
	 */
	public static final void updateNKernelNodes() {
		int i;
		int kernel;
		
		Math.set(nKernelNodes,nPathKernels-1,-1);
		kernel=blockKernels[0]; nKernelNodes[kernel]=0;
		for (i=1; i<nNodesInKernel; i++) {
			if (blockKernels[i]==kernel) continue;
			kernel=blockKernels[i];
			nKernelNodes[kernel]=i;
		}
	}
	
	
	/**
	 * @param threshold minimum length of an overlap;
	 * @return TRUE iff $kernel1$ and $kernel2$ have a pair of non-overlapping nodes in 
	 * the same read.
	 */
	private static final boolean haveNonoverlappingPair(int kernel1, int kernel2, int threshold) {
		int i, j, p, q;
		int read, start, end;
		IntervalGraph.Node tmpNode;
		
		// Remark: we search rather than using $nKernelNodes$, since we don't assume that
		// $updateNKernelNodes()$ has been executed yet.
		p=Arrays.binarySearch(blockKernels,0,nNodesInKernel,kernel1);
		while (p>=0 && blockKernels[p]==kernel1) p--;
		p++;
		q=Arrays.binarySearch(blockKernels,0,nNodesInKernel,kernel2);
		while (q>=0 && blockKernels[q]==kernel2) q--;
		q++;
		
		for (i=p; i<nNodesInKernel; i++) {
			if (blockKernels[i]!=kernel1) break;
			tmpNode=IntervalGraph.nodesArray[nodesInKernel[i]];
			read=tmpNode.read; start=tmpNode.start; end=tmpNode.end;
			for (j=q; j<nNodesInKernel; j++) {
				if (blockKernels[j]!=kernel2) break;
				tmpNode=IntervalGraph.nodesArray[nodesInKernel[j]];
				if (tmpNode.read==read && Intervals.intersectionLength(tmpNode.start,tmpNode.end,start,end)<threshold) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Removes all arcs $(u,v)$ of the kernel graph such that in no orientation the 
	 * following holds: all kernel nodes of $v$ that are mapped to a bidirected graph node
	 * are fully contained in the basin of $u$.
	 *
	 * @return the number of arcs removed by the procedure.
	 */
	private static final int buildKernelGraph_removeArcs() {
		byte tmpB;
		int i, j, k, h, p, q;
		int neighbor, last, out, tmpI;
		
		out=0;
		for (i=0; i<nPathKernels; i++) {
			last=lastInNeighbor[i];
			for (j=0; j<=last; j++) {
				neighbor=inNeighbors[i][j];
				if (inNeighborFlags[i][(j<<2)+1]==0 && inNeighborFlags[i][(j<<2)+3]==0) {
					out++;
					inNeighbors[i][j]=-1-neighbor;
					for (k=0; k<=lastOutNeighbor[neighbor]; k++) {
						if (outNeighbors[neighbor][k]==i) {
							outNeighbors[neighbor][k]=-1-outNeighbors[neighbor][k];
							break;
						}
					}
				}
			}
		}
		for (i=0; i<nPathKernels; i++) {
			k=-1;
			for (j=0; j<=lastInNeighbor[i]; j++) {
				if (inNeighbors[i][j]<0) continue;
				k++;
				if (k!=j) {
					tmpI=inNeighbors[i][k];
					inNeighbors[i][k]=inNeighbors[i][j];
					inNeighbors[i][j]=tmpI;
					p=k<<2; q=j<<2;
					for (h=0; h<=3; h++) {
						tmpB=inNeighborFlags[i][p+h];
						inNeighborFlags[i][p+h]=inNeighborFlags[i][q+h];
						inNeighborFlags[i][q+h]=tmpB;
					}
				}
			}
			lastInNeighbor[i]=k;
		}
		for (i=0; i<nPathKernels; i++) {
			k=-1;
			for (j=0; j<=lastOutNeighbor[i]; j++) {
				if (outNeighbors[i][j]<0) continue;
				k++;
				if (k!=j) {
					tmpI=outNeighbors[i][k];
					outNeighbors[i][k]=outNeighbors[i][j];
					outNeighbors[i][j]=tmpI;
				}
			}
			lastOutNeighbor[i]=k;
		}
		return out;
	}
	
	
	/**
	 * Removes insertions from $inNeighborFlags$ if two non-periodic kernels have similar 
	 * lengths. Kernels with similar lengths can have an insertion in the kernel graph 
	 * because of the distance threshold used for deciding identity between positions 
	 * throughout the code. E.g. the kernel nodes from which the kernel tag of an in-
	 * neighbor of a kernel was propagated, might be a bit larger than the real repeat 
	 * they capture, and a kernel in their basin might be a bit shorter, thus it might be 
	 * marked as a substring of the kernel nodes even if it coincides with an end of the 
	 * true repeat, thus the kernel tag might not be copied to its $pathsWith*$ arrays. 
	 * This can cause a real containment to be mistaken for an insertion.
	 *
	 * It might happen that two kernels with similar lengths are such that the longer 
	 * contains a short substring that is not in the shorter (e.g. when the longer is
	 * $UV$ and the shorter is $V$). We try to detect this by computing the difference
	 * between the set of in-neighbors of the shorter and of the longer.
	 *
	 * @return number of arcs for which insertions have been removed.
	 */
	private static final int buildKernelGraph_removeInsertions() {
		final int LENGTH_THRESHOLD = IO.quantum<<1;  // Arbitrary
		final double IN_NEIGHBORS_THRESHOLD = 0.5;  // Arbitrary
		int i, j, p;
		int length, neighbor, out;
		
		out=0;
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) continue;
			length=pathKernelLengths[i];
			for (j=0; j<=lastInNeighbor[i]; j++) {
				neighbor=inNeighbors[i][j];
				if (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0)) continue;
				if (pathKernelLengths[neighbor]-length<=LENGTH_THRESHOLD && lastInNeighbor[neighbor]>=0 && missingInNeighbors(i,j)<(lastInNeighbor[i]+1)*IN_NEIGHBORS_THRESHOLD) {
					p=j<<2;
					inNeighborFlags[i][p+0]=1;
					inNeighborFlags[i][p+2]=1;
					out++;
				}
			}
		}
		return out;
	}
	
	
	/**
	 * @return the number of in-neighbors of $shorterKernel$ that are not also in-
	 * neighbors of the kernel at $inNeighbors[shorterKernel][p]$.
	 */
	private static final int missingInNeighbors(int shorterKernel, int p) {
		int i1, i2, j;
		int p1, p2, out;
		final int longerKernel = inNeighbors[shorterKernel][p];
		final int last1 = lastInNeighbor[shorterKernel];
		final int last2 = lastInNeighbor[longerKernel];
		
		if (last2==-1) return last1+1;
		i1=0; i2=0; out=0;
		while (i1<=last1 && i2<=last2) {
			if (i1==p) i1++;
			else if (inNeighbors[shorterKernel][i1]<inNeighbors[longerKernel][i2]) { out++; i1++; }
			else if (inNeighbors[shorterKernel][i1]>inNeighbors[longerKernel][i2]) i2++;
			else {
				p1=i1<<2; p2=i2<<2;
				if (inNeighborFlags[shorterKernel][p1+1]!=inNeighborFlags[longerKernel][p2+1] || inNeighborFlags[shorterKernel][p1+3]!=inNeighborFlags[longerKernel][p2+3]) out++;
				i1++; i2++;
			}
		}
		out+=last1+1-i1;
		return out;
	}
	
	
	/**
	 * Stores in $permutationWindows[i]$ all sorted alignments between $nodes[i]$ and any
	 * other element of $nodes$. See $buildPermutationWindows_addWindows()$ for details.
	 *
	 * Remark: the procedure uses the global variables $shortPeriodAlignments,
	 * shortPeriodWindows,trackPointers$, which should have already been initialized.
	 *
	 * @param permutationWindows,lastPermutationWindow with at least $nNodes$ rows each.
	 */
	private static final void buildKernelGraph_buildPermutationWindows(int[][] permutationWindows, int[] lastPermutationWindow, String alignmentsFilePath, int totalNAlignments, IntervalGraph.Node[] nodes, int nNodes) throws IOException {
		final int BLOCK_SIZE = 6;
		final int OUTPUT_EVERY_ALIGNMENTS = 1000;  // Arbitrary
		final int CAPACITY = 10*BLOCK_SIZE;  // Arbitrary
		final int minRead = nodes[0].read;
		final int maxRead = nodes[nNodes-1].read;
		final int minStart = nodes[0].start;
		final int previousOrder = IntervalGraph.Node.order;
		boolean isConcatenationFrom, isConcatenationTo;
		int i, j, p, q;
		int max, maxEnd, nAlignments, operation, readA, readB, nWindows;
		String str;
		BufferedReader alignmentsFile;
		IntervalGraph.Node fromNode, toNode;
		int[] state = new int[7];
		PermutationWindow[] tmpWindows;
		
		// Building $permutationWindows$.
		maxEnd=0; i=nNodes-1;
		while (i>=0 && nodes[i].read==maxRead) maxEnd=Math.max(maxEnd,nodes[i--].end);
		for (i=0; i<nNodes; i++) permutationWindows[i] = new int[CAPACITY];
		Math.set(lastPermutationWindow,nNodes-1,-1);
		alignmentsFile = new BufferedReader(new FileReader(alignmentsFilePath),IO.BUFFER_SIZE);
		// Skipping the first two lines
		alignmentsFile.readLine();
		alignmentsFile.readLine();
		str=alignmentsFile.readLine();
		nAlignments=0;
		state[0]=-1; state[1]=-1; state[2]=-1; state[3]=-1; state[4]=-1; state[5]=-1; state[6]=-1;
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		while (str!=null) {
			operation=buildNodes2Kernel_skipAlignment(str,minRead,maxRead,minStart,maxEnd,nodes,nNodes,nodes,nNodes,tmpNode,state);
			if (operation==0) break;
			else if (operation==1 || shortPeriodAlignments[nAlignments]) {
				nAlignments++;
				if (nAlignments%OUTPUT_EVERY_ALIGNMENTS==0) System.err.println("buildPermutationWindows> "+IO.getPercent(nAlignments,totalNAlignments)+"%");
				str=alignmentsFile.readLine();
				continue;
			}
			p=state[0]; q=state[1];
			readA=Alignments.readA-1; readB=Alignments.readB-1;
			for (i=p-1; i>=0; i--) {
				fromNode=nodes[i];
				if (fromNode.read!=readA) break;
				if ((fromNode.type==Constants.INTERVAL_PERIODIC && fromNode.period>0 && !fromNode.hasLongPeriod) || fromNode.lastKernel==-1 || fromNode.end<=Alignments.startA) continue;
				isConcatenationFrom=(fromNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readA,fromNode.start,fromNode.end,shortPeriodWindows,trackPointers);
				for (j=q-1; j>=0; j--) {
					toNode=nodes[j];
					if (toNode.read!=readB) break;
					if ((toNode.type==Constants.INTERVAL_PERIODIC && toNode.period>0 && !toNode.hasLongPeriod) || toNode.lastKernel==-1 || toNode.end<=Alignments.startB) continue;
					isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readB,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
					buildPermutationWindows_addWindows(fromNode,i,isConcatenationFrom,toNode,j,isConcatenationTo,permutationWindows,lastPermutationWindow);
				}
				for (j=q; j<nNodes; j++) {
					toNode=nodes[j];
					if (toNode.read!=readB || toNode.start>=Alignments.endB) break;
					if ((toNode.type==Constants.INTERVAL_PERIODIC && toNode.period>0 && !toNode.hasLongPeriod) || toNode.lastKernel==-1) continue;
					isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readB,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
					buildPermutationWindows_addWindows(fromNode,i,isConcatenationFrom,toNode,j,isConcatenationTo,permutationWindows,lastPermutationWindow);
				}
			}
			for (i=p; i<nNodes; i++) {
				fromNode=nodes[i];
				if (fromNode.read!=readA || fromNode.start>=Alignments.endA) break;
				if ((fromNode.type==Constants.INTERVAL_PERIODIC && fromNode.period>0 && !fromNode.hasLongPeriod) || fromNode.lastKernel==-1 || fromNode.end<=Alignments.startA) continue;
				isConcatenationFrom=(fromNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readA,fromNode.start,fromNode.end,shortPeriodWindows,trackPointers);
				for (j=q-1; j>=0; j--) {
					toNode=nodes[j];
					if (toNode.read!=readB) break;
					if ((toNode.type==Constants.INTERVAL_PERIODIC && toNode.period>0 && !toNode.hasLongPeriod) || toNode.lastKernel==-1 || toNode.end<=Alignments.startB) continue;
					isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readB,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
					buildPermutationWindows_addWindows(fromNode,i,isConcatenationFrom,toNode,j,isConcatenationTo,permutationWindows,lastPermutationWindow);
				}
				for (j=q; j<nNodes; j++) {
					toNode=nodes[j];
					if (toNode.read!=readB || toNode.start>=Alignments.endB) break;
					if ((toNode.type==Constants.INTERVAL_PERIODIC && toNode.period>0 && !toNode.hasLongPeriod) || toNode.lastKernel==-1) continue;
					isConcatenationTo=(toNode.type!=Constants.INTERVAL_PERIODIC)&&IntervalGraph.inShortPeriodTrack(readB,toNode.start,toNode.end,shortPeriodWindows,trackPointers);
					buildPermutationWindows_addWindows(fromNode,i,isConcatenationFrom,toNode,j,isConcatenationTo,permutationWindows,lastPermutationWindow);
				}
			}
			nAlignments++;
			if (nAlignments%OUTPUT_EVERY_ALIGNMENTS==0) System.err.println("buildPermutationWindows> "+IO.getPercent(nAlignments,totalNAlignments)+"%");
			str=alignmentsFile.readLine();
		}
		alignmentsFile.close();
		IntervalGraph.Node.order=previousOrder;
		
		// Sorting
		max=0;
		for (i=0; i<nNodes; i++) {
			if (lastPermutationWindow[i]>max) max=lastPermutationWindow[i];
		}
		tmpWindows = new PermutationWindow[(max+1)/BLOCK_SIZE];
		for (i=0; i<tmpWindows.length; i++) tmpWindows[i] = new PermutationWindow();
		for (i=0; i<nNodes; i++) {
			if (lastPermutationWindow[i]==-1 || lastPermutationWindow[i]==BLOCK_SIZE-1) continue;
			nWindows=0;
			for (j=0; j<lastPermutationWindow[i]; j+=BLOCK_SIZE) tmpWindows[nWindows++].set(permutationWindows[i][j+0],permutationWindows[i][j+1],permutationWindows[i][j+2],permutationWindows[i][j+3],permutationWindows[i][j+4],permutationWindows[i][j+5]);
			Arrays.sort(tmpWindows,0,nWindows);
			tmpWindows[0].writeTo(permutationWindows[i],0);
			for (j=1; j<nWindows; j++) tmpWindows[j].writeTo(permutationWindows[i],j*BLOCK_SIZE);
		}
	}
	
	
	/**
	 * Adds the alignment string loaded in $Alignments$ to row $fromNodeIndex$ of 
	 * $windows$. The alignment is stored in the following block: 
	 * 0=ID of the other node in the alignment; 1=orientation of the alignment (0=fwd,
	 * 1=rc); 2=start of the common window inside $fromNode$; 3=end inside $fromNode$; 
	 * 4=start in $toNode$ (in the forward orientation); 5=end in $toNode$ (in the forward
	 * orientation).
	 *
	 * Remark: since the same alignment occurs twice (symmetrized) in the alignments file,
	 * the same block is added twice, to different rows.
	 *
	 * @param *Index position of the node in the nodes array;
	 * @param *NodeIsConcatenation tells whether $*Node$ is a concatenation of several
	 * distinct short-period intervals.
	 */
	private static final void buildPermutationWindows_addWindows(IntervalGraph.Node fromNode, int fromNodeIndex, boolean fromNodeIsConcatenation, IntervalGraph.Node toNode, int toNodeIndex, boolean toNodeIsConcatenation, int[][] windows, int[] lastWindow) throws IOException {
		final int MIN_WINDOW_SIZE = Alignments.minAlignmentLength>>2;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int BLOCK_SIZE = 6;
		final int INCREMENT = 10;  // In number of blocks. Arbitrary.
		int i;
		int last, windowStart, windowEnd, windowStartPrime, windowEndPrime;
		final int readFrom = fromNode.read;
		final int readTo = toNode.read;
		
		windowStart=Math.max(Alignments.startA,fromNode.start);
		windowEnd=Math.min(Alignments.endA,fromNode.end);
		Intervals.project(toNode.start,toNode.end,Alignments.startB,Alignments.endB,Alignments.startA,Alignments.endA,Alignments.orientation,tmpIO,0);
		windowStart=Math.max(windowStart,tmpIO[0]);
		windowEnd=Math.min(windowEnd,tmpIO[1]);
		if (windowEnd-windowStart+1<MIN_WINDOW_SIZE) {
			// Remark: it can be $windowStart>windowEnd$ since e.g. $[tmpIO[0]..tmpIO[1]]$
			// might have no intersection with $fromNode$.
			return;
		}
		if (!fromNodeIsConcatenation && IntervalGraph.inShortPeriodTrack(readFrom,windowStart,windowEnd,shortPeriodWindows,trackPointers)) return;
		windowStartPrime=Math.max(Alignments.startB,toNode.start);
		windowEndPrime=Math.min(Alignments.endB,toNode.end);
		Intervals.project(fromNode.start,fromNode.end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,tmpIO,0);
		windowStartPrime=Math.max(windowStartPrime,tmpIO[0]);
		windowEndPrime=Math.min(windowEndPrime,tmpIO[1]);
		if (!toNodeIsConcatenation && IntervalGraph.inShortPeriodTrack(readTo,windowStartPrime,windowEndPrime,shortPeriodWindows,trackPointers)) return;
		
		// Adding
		last=lastWindow[fromNodeIndex];
		lastWindow[fromNodeIndex]+=BLOCK_SIZE;
		if (lastWindow[fromNodeIndex]>=windows[fromNodeIndex].length) {
			int[] newWindows = new int[windows[fromNodeIndex].length+INCREMENT*BLOCK_SIZE];
			System.arraycopy(windows[fromNodeIndex],0,newWindows,0,last+1);
			windows[fromNodeIndex]=newWindows;
		}
		windows[fromNodeIndex][last+1]=toNodeIndex;
		windows[fromNodeIndex][last+2]=Alignments.orientation?0:1;
		windows[fromNodeIndex][last+3]=windowStart;
		windows[fromNodeIndex][last+4]=windowEnd;
		windows[fromNodeIndex][last+5]=windowStartPrime;
		windows[fromNodeIndex][last+6]=windowEndPrime;
	}
	
	
	private static class PermutationWindow implements Comparable {
		public boolean orientation;
		public int otherNode, start, end, otherStart, otherEnd;
		
		public PermutationWindow() {
			this.otherNode=-1;
			this.orientation=false;
			this.start=-1;
			this.end=-1;
			this.otherStart=-1;
			this.otherEnd=-1;
		}
		
		public PermutationWindow(int on, int or, int s, int e, int os, int oe) {
			set(on,or,s,e,os,oe);
		}
		
		public void set(int on, int or, int s, int e, int os, int oe) {
			otherNode=on;
			orientation=or==0;
			start=s;
			end=e;
			otherStart=os;
			otherEnd=oe;
		}
		
		public int compareTo(Object other) {
			PermutationWindow otherWindow = (PermutationWindow)other;
			if (otherNode<otherWindow.otherNode) return -1;
			else if (otherNode>otherWindow.otherNode) return 1;
			if (orientation && !otherWindow.orientation) return -1;
			else if (!orientation && otherWindow.orientation) return 1;
			if (start<otherWindow.start) return -1;
			else if (start>otherWindow.start) return 1;
			if (otherStart<otherWindow.otherStart) return -1;
			else if (otherStart>otherWindow.otherStart) return 1;
			return 0;
		}
		
		public boolean equals(Object other) {
			PermutationWindow otherWindow = (PermutationWindow)other;
			return orientation==otherWindow.orientation && otherNode==otherWindow.otherNode && start==otherWindow.start && end==otherWindow.end && otherStart==otherWindow.otherStart && otherEnd==otherWindow.otherEnd;
		}
		
		/**
		 * @param from first position of $array$.
		 */
		public void writeTo(int[] array, int from) {
			array[from+0]=otherNode;
			array[from+1]=orientation?0:1;
			array[from+2]=start;
			array[from+3]=end;
			array[from+4]=otherStart;
			array[from+5]=otherEnd;
		}
	}
	
	
	/**
	 * Some pairs of path kernels might be connected by an arc in the kernel graph, not 
	 * because one is contained in the other, but because one is covered by substrings of
	 * the other (e.g. because both consist of a set of smaller repeats, and one is a 
	 * permutation of the repeats of the other). The procedure tries to detect this by
	 * finding pairs of interval graph nodes $A,B$ that belong to the sequence of some 
	 * path kernels, and such that: (1) there are two, not necessarily disjoint substrings 
	 * $A[i..j],A[i'..j']$ of $A$, with $i'>i$, that are aligned to substrings $B[p..q],
	 * B[p'..q']$ of $B$, with $p' <= p$; (2) $A[i..j],A[i'..j']$ are never aligned to any
	 * pair of substrings $B[x..y],B[x'..y']$ with $x'>x$.
	 *
	 * Remark: the procedure handles the case in which $B[p..q]$ and $B[p'..q']$ coincide.
	 *
	 * Remark: this simple filter does not take into account the full sequence of a path
	 * kernel (which might consist of multiple nodes), nor the basins of path kernels.
	 * This would make the code significantly more complex.
	 *
	 * Remark: the long period of a path kernel might be a permutation of the long period 
	 * of another one, and this is detected by the procedure. The case in which the long
	 * period of a path kernel is a substring of the long period of another one is handled
	 * by $buildKernelGraph_removePermutations_longPeriods()$.
	 *
	 * @param permutationWindows output of $buildKernelGraph_buildPermutationWindows()$;
	 * @param nodes assumed to be sorted by $read,start$;
	 * @return number of arcs of the kernel graph removed by the procedure.
	 */
	private static final int buildKernelGraph_removePermutations(IntervalGraph.Node[] nodes, int nNodes, int[][] permutationWindows, int[] lastPermutationWindow, int threshold) throws IOException {
		final int BLOCK_SIZE = 6;
		final int DISTANCE_THRESHOLD = threshold<<1;  // Arbitrary
		final int MIN_PERMUTATION_UNIT = Alignments.minAlignmentLength;  // Arbitrary
		boolean hasPermutation, found, hasLongPeriodI, hasLongPeriodJ;
		int i, j, k;
		int jPrime, kPrime, periodI, periodJ, last, distance, nRemovedArcs;
		int startJ, endJ, orientationJ, otherStartJ, otherEndJ, otherNodeJ, projectedStartJ, projectedEndJ;
		int startK, endK, orientationK, otherStartK, otherEndK, otherNodeK, projectedStartK, projectedEndK;
		

if (IO.SHOW_INTERACTIVE) {
System.err.println("removePermutations> nodes:");
for (int x=0; x<nNodes; x++) System.err.println(x+": "+nodes[x]);
System.err.println();
System.err.println("removePermutations> permutationWindows:");
for (int x=0; x<nNodes; x++) {
	System.err.print(x+":   ");
	for (int y=0; y<=lastPermutationWindow[x]; y+=BLOCK_SIZE) {
		System.err.print(permutationWindows[x][y+0]+",");
		System.err.print(permutationWindows[x][y+1]+",");
		System.err.print(permutationWindows[x][y+2]+",");
		System.err.print(permutationWindows[x][y+3]+",");
		System.err.print(permutationWindows[x][y+4]+",");
		System.err.print(permutationWindows[x][y+5]+", ");
	}
	System.err.println();
}		
}		
		
		
		nRemovedArcs=0;
		for (i=0; i<nNodes; i++) {
			last=lastPermutationWindow[i];
			if (last==BLOCK_SIZE-1) continue;
			periodI=nodes[i].period;
			hasLongPeriodI=false;
			if (nodes[i].type==Constants.INTERVAL_PERIODIC) {
				if (nodes[i].hasLongPeriod) hasLongPeriodI=true;
				else if (periodI>0) continue;
			}
			j=0;
			while (j<=last) {
				otherNodeJ=permutationWindows[i][j+0];
				periodJ=nodes[otherNodeJ].period;
				hasLongPeriodJ=false;
				if (nodes[otherNodeJ].type==Constants.INTERVAL_PERIODIC) {
					if (nodes[otherNodeJ].hasLongPeriod) hasLongPeriodJ=true;
					else if (periodJ>0) {
						j+=BLOCK_SIZE;
						continue;
					}
				}
				orientationJ=permutationWindows[i][j+1];
				startJ=permutationWindows[i][j+2]; 
				endJ=permutationWindows[i][j+3]; 
				otherStartJ=permutationWindows[i][j+4];
				otherEndJ=permutationWindows[i][j+5];
				hasPermutation=false;
				for (k=j+BLOCK_SIZE; k<=last; k+=BLOCK_SIZE) {
					otherNodeK=permutationWindows[i][k+0];
					orientationK=permutationWindows[i][k+1];
					startK=permutationWindows[i][k+2]; 
					endK=permutationWindows[i][k+3]; 
					if (Intervals.isApproximatelyContained(startK,endK,startJ,endJ)) continue;
					otherStartK=permutationWindows[i][k+4];
					otherEndK=permutationWindows[i][k+5];
					if (otherNodeK!=otherNodeJ || orientationK!=orientationJ || (hasLongPeriodI && startK>startJ+periodI-MIN_PERMUTATION_UNIT)) break;
					distance=startK-startJ;
					if ( distance>threshold && 
						 ( (orientationJ==0 && otherStartK<=otherStartJ+threshold && (hasLongPeriodJ?Math.abs(otherStartK,otherStartJ)<=periodJ-MIN_PERMUTATION_UNIT:true)) ||
						   (orientationJ==1 && otherEndK>=otherEndJ-threshold && (hasLongPeriodJ?Math.abs(otherEndK,otherEndJ)<=periodJ-MIN_PERMUTATION_UNIT:true))
						 ) &&
						 (hasLongPeriodI||hasLongPeriodJ?Intervals.straddles(otherStartK,otherEndK,otherStartJ,otherEndJ):true)
					   ) {
						jPrime=j+BLOCK_SIZE;
						while (jPrime<=last && permutationWindows[i][jPrime+0]==otherNodeJ && permutationWindows[i][jPrime+1]==orientationJ && permutationWindows[i][jPrime+2]<=startJ+threshold) jPrime+=BLOCK_SIZE;
						jPrime-=BLOCK_SIZE;
						found=false;
						while (jPrime>=0) {
							if (permutationWindows[i][jPrime+0]!=otherNodeJ || permutationWindows[i][jPrime+1]!=orientationJ) break;
							if (permutationWindows[i][jPrime+3]<endJ-threshold) {
								jPrime-=BLOCK_SIZE;
								continue;
							}
							Intervals.project(startJ,endJ,permutationWindows[i][jPrime+2],permutationWindows[i][jPrime+3],permutationWindows[i][jPrime+4],permutationWindows[i][jPrime+5],permutationWindows[i][jPrime+1]==0,tmpIO,0);
							projectedStartJ=tmpIO[0]; projectedEndJ=tmpIO[1];
							kPrime=k+BLOCK_SIZE;
							while (kPrime<=last && permutationWindows[i][kPrime+0]==otherNodeJ && permutationWindows[i][kPrime+1]==orientationJ && permutationWindows[i][kPrime+2]<=startK+threshold) kPrime+=BLOCK_SIZE;
							kPrime-=BLOCK_SIZE;
							while (kPrime>=0) {
								if (permutationWindows[i][kPrime+0]!=otherNodeJ || permutationWindows[i][kPrime+1]!=orientationJ) break;
								if (permutationWindows[i][kPrime+3]<endK-threshold) {
									kPrime-=BLOCK_SIZE;
									continue;
								}
								Intervals.project(startK,endK,permutationWindows[i][kPrime+2],permutationWindows[i][kPrime+3],permutationWindows[i][kPrime+4],permutationWindows[i][kPrime+5],permutationWindows[i][kPrime+1]==0,tmpIO,0);
								projectedStartK=tmpIO[0]; projectedEndK=tmpIO[1];
								if ( (orientationJ==0 && projectedStartK>projectedStartJ+threshold && Math.abs(projectedStartK-projectedStartJ,distance)<=DISTANCE_THRESHOLD) ||
									 (orientationJ==1 && projectedEndK<projectedEndJ-threshold && Math.abs(projectedEndJ-projectedEndK,distance)<=DISTANCE_THRESHOLD)
								   ) {
									found=true;
									break;
								}
								kPrime-=BLOCK_SIZE;
							}
							if (found) break;
							jPrime-=BLOCK_SIZE;
						}
						if (!found) {
							hasPermutation=true;
							break;
						}
					}
				}
				if (hasPermutation) {
					nRemovedArcs+=buildKernelGraph_removePermutations_impl(nodes[i],nodes[otherNodeJ],orientationJ);
					j+=BLOCK_SIZE;
					while (j<=last && permutationWindows[i][j+0]==otherNodeJ && permutationWindows[i][j+1]==orientationJ) j+=BLOCK_SIZE;
				}
				else j+=BLOCK_SIZE;
			}
		}
		return nRemovedArcs;
	}

	
	/**
	 * @param alignmentOrientation orientation between $fromNode$ and $toNode$ in which a
	 * permutation was detected;
	 * @return the number of arcs in the kernel graph removed by the procedure.
	 */
	private static final int buildKernelGraph_removePermutations_impl(IntervalGraph.Node fromNode, IntervalGraph.Node toNode, int alignmentOrientation) {
		boolean found;
		byte orientationI, orientationJ, kernelsOrientation;
		int i, j, k;
		int last, kernelI, kernelJ, out;
		
		out=0;
		for (i=0; i<=fromNode.lastKernel; i++) {
			kernelI=fromNode.kernels[i];
			if (pathKernelLengths[kernelI]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernelI]==1)) continue;
			orientationI=fromNode.kernelOrientations[i];
			for (j=0; j<=toNode.lastKernel; j++) {
				kernelJ=toNode.kernels[j];
				if (kernelJ==kernelI || pathKernelLengths[kernelJ]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernelJ]==1)) continue;
				orientationJ=toNode.kernelOrientations[j];
				kernelsOrientation=getKernelsOrientation[orientationI][orientationJ][alignmentOrientation];
				found=false; last=lastInNeighbor[kernelI];
				for (k=0; k<=last; k++) {
					if (inNeighbors[kernelI][k]!=kernelJ) continue;
					found=true;
					if (kernelsOrientation==0) {
						if (inNeighborFlags[kernelI][(k<<2)+1]==1) {
							if (inNeighborFlags[kernelI][(k<<2)+3]==1) inNeighborFlags[kernelI][(k<<2)+1]=0;
							else {
								removeArc(kernelI,k);
								out++;
							}
						}
					}
					else if (kernelsOrientation==1) {
						if (inNeighborFlags[kernelI][(k<<2)+3]==1) {
							if (inNeighborFlags[kernelI][(k<<2)+1]==1) inNeighborFlags[kernelI][(k<<2)+3]=0;
							else {
								removeArc(kernelI,k);
								out++;
							}
						}
					}
					else {
						removeArc(kernelI,k);
						out++;
					}
					break;
				}
				if (found) continue;
				last=lastInNeighbor[kernelJ];
				for (k=0; k<=last; k++) {
					if (inNeighbors[kernelJ][k]!=kernelI) continue;
					if (kernelsOrientation==0) {
						if (inNeighborFlags[kernelJ][(k<<2)+1]==1) {
							if (inNeighborFlags[kernelJ][(k<<2)+3]==1) inNeighborFlags[kernelJ][(k<<2)+1]=0;
							else {
								removeArc(kernelJ,k);
								out++;
							}
						}
					}
					else if (kernelsOrientation==1) {
						if (inNeighborFlags[kernelJ][(k<<2)+3]==1) {
							if (inNeighborFlags[kernelJ][(k<<2)+1]==1) inNeighborFlags[kernelJ][(k<<2)+3]=0;
							else {
								removeArc(kernelJ,k);
								out++;
							}
						}
					}
					else {
						removeArc(kernelJ,k);
						out++;
					}
					break;
				}
			}
		}
		return out;
	}
	
	
	/**
	 * Removes the arc between $kernel$ and its $p$-th in-neighbor.
	 */
	private static final void removeArc(int kernel, int p) {
		int i, j;
		int iPrime, last;
		final int otherKernel = inNeighbors[kernel][p];
		
		last=lastInNeighbor[kernel];
		for (i=p; i<last; i++) inNeighbors[kernel][i]=inNeighbors[kernel][i+1];
		for (i=p; i<last; i++) {
			iPrime=i<<2;
			inNeighborFlags[kernel][iPrime+0]=inNeighborFlags[kernel][iPrime+0+4];
			inNeighborFlags[kernel][iPrime+1]=inNeighborFlags[kernel][iPrime+1+4];
			inNeighborFlags[kernel][iPrime+2]=inNeighborFlags[kernel][iPrime+2+4];
			inNeighborFlags[kernel][iPrime+3]=inNeighborFlags[kernel][iPrime+3+4];
		}
		lastInNeighbor[kernel]--;
		last=lastOutNeighbor[otherKernel];
		for (i=0; i<=last; i++) {
			if (outNeighbors[otherKernel][i]!=kernel) continue;
			for (j=i; j<last; j++) outNeighbors[otherKernel][j]=outNeighbors[otherKernel][j+1];
			lastOutNeighbor[otherKernel]--;
			break;
		}
	}
	
	
	/**
	 * Some pairs of long-period path kernels might be connected by an arc in the kernel 
	 * graph because the period of one is a substring of the period of the other. The 
	 * procedure tries to remove all such arcs based on average period estimates.
	 *
	 * @return number of arcs removed by the procedure.
	 */
	private static final int buildKernelGraph_removePermutations_longPeriods() {
		final double FRACTION = 0.8;  // Arbitrary
		byte tmpB;
		int i, j, k, h, p, q;
		int neighbor, out, tmpI, nPeriodic;
		int containedPeriod, containerPeriod, containedLength, containerLength;
		int[] selected;
		int[][] avgPeriods;
		
		// Computing avg. periods
		if (pathKernelPeriodic==null) return 0;
		nPeriodic=0;
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelPeriodic[i]!=0) nPeriodic++;
		}
		if (nPeriodic<=1) return 0;
		selected = new int[nPeriodic];
		nPeriodic=-1;
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelPeriodic[i]!=0) selected[++nPeriodic]=i;
		}
		nPeriodic++;
		avgPeriods = new int[nPeriodic][2];
		getAvgPeriods(selected,nPeriodic,avgPeriods);
		
		// Removing arcs
		out=0; h=0;
		for (i=0; i<nPathKernels; i++) {
			containedLength=pathKernelLengths[i];
			if (containedLength==0 || pathKernelPeriodic[i]==0) continue;
			while (h<nPeriodic && selected[h]<i) h++;
			containedPeriod=avgPeriods[h][0];
			if (containedPeriod==0) containedPeriod=containedLength;  // Unknown period
			else if (containedPeriod<PeriodicSubstrings.MIN_PERIOD_LONG) {
				continue;  // Short period
			}
			for (j=0; j<=lastInNeighbor[i]; j++) {
				neighbor=inNeighbors[i][j];
				containerLength=pathKernelLengths[neighbor];
				if (containerLength==0 || pathKernelPeriodic[neighbor]!=2) continue;
				containerPeriod=avgPeriods[Arrays.binarySearch(selected,0,nPeriodic,neighbor)][0];
				if (containerPeriod>0 && containedPeriod<containerPeriod*FRACTION) {
					out++;
					inNeighbors[i][j]=-1-neighbor;
					for (k=0; k<=lastOutNeighbor[neighbor]; k++) {
						if (outNeighbors[neighbor][k]==i) {
							outNeighbors[neighbor][k]=-1-i;
							break;
						}
					}	
				}
			}
		}
		for (i=0; i<nPathKernels; i++) {
			k=-1;
			for (j=0; j<=lastInNeighbor[i]; j++) {
				if (inNeighbors[i][j]<0) continue;
				k++;
				if (k!=j) {
					tmpI=inNeighbors[i][k];
					inNeighbors[i][k]=inNeighbors[i][j];
					inNeighbors[i][j]=tmpI;
					p=k<<2; q=j<<2;
					for (h=0; h<=3; h++) {
						tmpB=inNeighborFlags[i][p+h];
						inNeighborFlags[i][p+h]=inNeighborFlags[i][q+h];
						inNeighborFlags[i][q+h]=tmpB;
					}
				}
			}
			lastInNeighbor[i]=k;
		}
		for (i=0; i<nPathKernels; i++) {
			k=-1;
			for (j=0; j<=lastOutNeighbor[i]; j++) {
				if (outNeighbors[i][j]<0) continue;
				k++;
				if (k!=j) {
					tmpI=outNeighbors[i][k];
					outNeighbors[i][k]=outNeighbors[i][j];
					outNeighbors[i][j]=tmpI;
				}
			}
			lastOutNeighbor[i]=k;
		}
		return out;
	}
	
	
	/**
	 * Remark: the procedure uses $tmpArray1$, and it scans all interval graph nodes.
	 *
	 * @param selected kernel tags, assumed to be of periodic path kernels;
	 * @param avgPeriods output array: 0=avg. period, 1=number of interval graph nodes
	 * used to compute the average.
	 */
	private static final void getAvgPeriods(int[] selected, int nSelected, int[][] avgPeriods) {
		int i, j, k;
		int last1;
		IntervalGraph.Node node;
		
		Math.set(avgPeriods,nSelected-1,0);
		if (tmpArray1==null || tmpArray1.length<nSelected) tmpArray1 = new int[nSelected];
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastKernel==-1 || node.type!=Constants.INTERVAL_PERIODIC || node.period<=0) continue;
			last1=Math.setIntersection(node.kernels,0,node.lastKernel,selected,0,nSelected-1,tmpArray1,0);
			if (last1==-1) continue;
			k=0;
			for (j=0; j<=last1; j++) {
				while (k<nSelected && selected[k]<tmpArray1[j]) k++;
				avgPeriods[k][0]+=node.period;
				avgPeriods[k][1]++;
			}
		}
		for (i=0; i<nSelected; i++) {
			if (avgPeriods[i][1]!=0) avgPeriods[i][0]/=avgPeriods[i][1];
		}
	}
	
	
	/**
	 * Some pairs of path kernels might be connected by an arc in the kernel graph, 
	 * because one equals a substring of the other with some parts deleted. The procedure 
	 * tries to detect this by finding pairs of interval graph nodes $A,B$ that belong to 
	 * the sequence of some path kernels, and such that $A$ is fully covered by alignments
	 * with $B$, but the union of such alignments in $B$ is not a compact interval.
	 *
	 * Remark: this simple filter does not take into account the full sequence of a path
	 * kernel (which might consist of multiple nodes).
	 *
	 * Remark: the long period of a path kernel might be equal to the long period of 
	 * another one, with some deletions, and this is detected by the procedure. 
	 *
	 * Remark: deletions of short-period intervals are not detected. If an interval is not
	 * periodic, but it contains a short-period part and that part has a deletion, this is
	 * not detected.
	 *
	 * @param permutationWindows output of $buildKernelGraph_buildPermutationWindows()$;
	 * @param nodes assumed to be sorted by $read,start$;
	 * @return number of arcs of the kernel graph removed by the procedure.
	 */
	private static final int buildKernelGraph_removeDeletions(IntervalGraph.Node[] nodes, int nNodes, int[][] permutationWindows, int[] lastPermutationWindow, int threshold) {
		final int BLOCK_SIZE = 6;
		final int DISTANCE_THRESHOLD = threshold<<1;  // Arbitrary
		boolean coversNodeI, coversNodeIPrime, hasLongPeriod;
		int i, j, k, p, q;
		int max, last, lastPair, nRemovedArcs, firstPair;
		int nodeStart, nodeEnd, nodePeriod, otherNode, orientation, start, end;
		int otherNodeJ, orientationJ, startJ, endJ, otherStartJ, otherEndJ;
		Pair[] pairs;
		
		
if (IO.SHOW_INTERACTIVE) {
System.err.println("removeDeletions> nodes:");
for (int x=0; x<nNodes; x++) System.err.println(x+": "+nodes[x]);
System.err.println();
System.err.println("removeDeletions> permutationWindows:");
for (int x=0; x<nNodes; x++) {
	System.err.print(x+":   ");
	for (int y=0; y<=lastPermutationWindow[x]; y+=BLOCK_SIZE) {
		System.err.print(permutationWindows[x][y+0]+",");
		System.err.print(permutationWindows[x][y+1]+",");
		System.err.print(permutationWindows[x][y+2]+",");
		System.err.print(permutationWindows[x][y+3]+",");
		System.err.print(permutationWindows[x][y+4]+",");
		System.err.print(permutationWindows[x][y+5]+", ");
	}
	System.err.println();
}
}
		
		max=0;
		for (i=0; i<nNodes; i++) {
			if (lastPermutationWindow[i]>max) max=lastPermutationWindow[i];
		}
		pairs = new Pair[(max+1)/BLOCK_SIZE];
		for (i=0; i<pairs.length; i++) pairs[i] = new Pair();
		nRemovedArcs=0;
		for (i=0; i<nNodes; i++) {
			last=lastPermutationWindow[i];
			if (last==BLOCK_SIZE-1) continue;
			nodeStart=nodes[i].start; nodeEnd=nodes[i].end; nodePeriod=nodes[i].period;
			hasLongPeriod=false;
			if (nodes[i].type==Constants.INTERVAL_PERIODIC) {
				if (nodes[i].hasLongPeriod) hasLongPeriod=true;
				else if (nodePeriod>0) continue;
			}
			j=0; otherNode=-1; orientation=-1; start=-1; end=-1; startJ=-1; endJ=-1;
			coversNodeI=false; lastPair=-1;
			while (j<=last) {
				otherNodeJ=permutationWindows[i][j+0];
				if (nodes[otherNodeJ].type==Constants.INTERVAL_PERIODIC && nodes[otherNodeJ].period>0 && !nodes[otherNodeJ].hasLongPeriod) {
					j+=BLOCK_SIZE;
					continue;
				}
				orientationJ=permutationWindows[i][j+1];
				startJ=permutationWindows[i][j+2];
				endJ=permutationWindows[i][j+3]; 
				otherStartJ=permutationWindows[i][j+4];
				otherEndJ=permutationWindows[i][j+5];
				if (otherNode==-1) {
					otherNode=otherNodeJ; orientation=orientationJ; start=startJ;
					if (start>nodeStart+threshold) {
						coversNodeI=false;
						j+=BLOCK_SIZE;
						while (j<=last && permutationWindows[i][j+0]==otherNodeJ && permutationWindows[i][j+1]==orientationJ) j+=BLOCK_SIZE;
						continue;
					}
					end=endJ; coversNodeI=true;
					lastPair=0; pairs[0].start=otherStartJ; pairs[0].end=otherEndJ;
					pairs[0].startJ=startJ; pairs[0].endJ=endJ;
					j+=BLOCK_SIZE;
					continue;
				}
				else if (otherNodeJ!=otherNode || orientationJ!=orientation) {
					if (coversNodeI && end>=nodeEnd-threshold) {
						Pair.order=Pair.ORDER_START;
						if (lastPair>0) Arrays.sort(pairs,0,lastPair+1);
						coversNodeIPrime=false;
						p=pairs[0].start; q=pairs[0].end; firstPair=0;
						for (k=1; k<=lastPair; k++) {
							if (pairs[k].start>q+threshold) {
								if (removeDeletions_covers(pairs,firstPair,k-1,startJ,endJ,threshold)) coversNodeIPrime=true;
								p=pairs[k].start; q=pairs[k].end; firstPair=k;								
							}
							else q=Math.max(q,pairs[k].end);
						}
						if (removeDeletions_covers(pairs,firstPair,lastPair,startJ,endJ,threshold)) coversNodeIPrime=true;
						if (!coversNodeIPrime) nRemovedArcs+=buildKernelGraph_removePermutations_impl(nodes[i],nodes[otherNode],orientation);
					}
					// Next iteration
					otherNode=otherNodeJ; orientation=orientationJ; start=startJ;
					if (start>nodeStart+threshold) {
						coversNodeI=false;
						j+=BLOCK_SIZE;
						while (j<=last && permutationWindows[i][j+0]==otherNodeJ && permutationWindows[i][j+1]==orientationJ) j+=BLOCK_SIZE;
						continue;
					}
					end=endJ; coversNodeI=true;
					lastPair=0; pairs[0].start=otherStartJ; pairs[0].end=otherEndJ;
					pairs[0].startJ=startJ; pairs[0].endJ=endJ;
					j+=BLOCK_SIZE;
					continue;
				}
				if (startJ>end+threshold) {
					coversNodeI=false;
					j+=BLOCK_SIZE;
					while (j<=last && permutationWindows[i][j+0]==otherNodeJ && permutationWindows[i][j+1]==orientationJ) j+=BLOCK_SIZE;
					continue;
				}
				end=Math.max(end,endJ);
				lastPair++; pairs[lastPair].start=otherStartJ; pairs[lastPair].end=otherEndJ;
				pairs[lastPair].startJ=startJ; pairs[lastPair].endJ=endJ;
				j+=BLOCK_SIZE;
			}
			// Last iteration
			if (coversNodeI && end>=nodeEnd-threshold) {
				Pair.order=Pair.ORDER_START;
				if (lastPair>0) Arrays.sort(pairs,0,lastPair+1);
				coversNodeIPrime=false;
				p=pairs[0].start; q=pairs[0].end; firstPair=0;
				for (k=1; k<=lastPair; k++) {
					if (pairs[k].start>q+threshold) {
						if (removeDeletions_covers(pairs,firstPair,k-1,startJ,endJ,threshold)) coversNodeIPrime=true;
						p=pairs[k].start; q=pairs[k].end; firstPair=k;								
					}
					else q=Math.max(q,pairs[k].end);
				}
				if (removeDeletions_covers(pairs,firstPair,lastPair,startJ,endJ,threshold)) coversNodeIPrime=true;
				if (!coversNodeIPrime) nRemovedArcs+=buildKernelGraph_removePermutations_impl(nodes[i],nodes[otherNode],orientation);
			}
		}
		return nRemovedArcs;
	}
	
	
	/**
	 * Remark: the procedure sorts $pairs[from..to]$ by $startJ$.
	 *
	 * @return TRUE iff the union of the $[startJ..endJ]$ intervals of all elements in
	 * $pairs[from..to]$ covers $[nodeStart..nodeEnd]$.
	 */
	private static final boolean removeDeletions_covers(Pair[] pairs, int from, int to, int nodeStart, int nodeEnd, int threshold) {
		int i;
		int endJ;
		
		Pair.order=Pair.ORDER_STARTJ;
		if (to>from) Arrays.sort(pairs,from,to+1);
		if (pairs[from].startJ>nodeStart+threshold) return false;
		endJ=pairs[from].endJ;
		for (i=from+1; i<=to; i++) {
			if (pairs[i].startJ>endJ+threshold) return false;
			endJ=Math.max(endJ,pairs[i].endJ);
		}
		return endJ>=nodeEnd-threshold;
	}
	
	
	private static class Pair implements Comparable {
		public static final int ORDER_START = 0;
		public static final int ORDER_STARTJ = 1;
		public static int order;
		
		int start, end, startJ, endJ;
		
		public int compareTo(Object other) {
			Pair otherPair = (Pair)other;
			if (order==ORDER_START) {
				if (start<otherPair.start) return -1;
				else if (start>otherPair.start) return 1;
			}
			else {
				if (startJ<otherPair.startJ) return -1;
				else if (startJ>otherPair.startJ) return 1;
			}
			return 0;
		}
	}
	
	
	/**
	 * Fits a density estimation tree on $kernelFrequency[*][0..2]$ and returns a value 
	 * that separates the first local maximum from the rest. The first local maximum is 
	 * assumed to correspond to concatenations of repeats that have been wrongly
	 * undersplit.
	 *
	 * Remark: the procedure assumes $tmpPoints$ to be of size at least $nPathKernels$.
	 */
	public static final int getKernelFrequencyThreshold() {
		final int DEFAULT_THRESHOLD = (IO.minOccurrencesInGenome*IO.coverage)<<2;  // Multiplication by 4 rather than by 2 arbitrary.
		final int N_CONSECUTIVE_POINTS = 2;  // Arbitrary
		final double QUANTILE = 0.75;  // Arbitrary
		final int MIN_INTERVAL_LENGTH_FACTOR = 100;  // Arbitrary
		final int MAX_DIFFERENCE = DEFAULT_THRESHOLD*3;  // Arbitrary
		int i;
		int lastFrequency, minIntervalLength, minLocalMaxDistance, nLocalMaximumLeaves;
		int lastLeaf, lastPoint, nHigh;
		double threshold;
		Leaf[] lvs;
		
		threshold=DEFAULT_THRESHOLD;
		if (nPathKernels==1) return (int)threshold;
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) continue;
			tmpPoints[i].position=kernelFrequency[i][0]+(kernelFrequency[i][1]+kernelFrequency[i][2])>>1;
			tmpPoints[i].mass=1;
		}
		lastFrequency=Points.sortAndCompact(tmpPoints,nPathKernels-1);
		if (lastFrequency+1<=Points.FEW_POINTS) {
			if (tmpPoints[lastFrequency].position>MAX_DIFFERENCE) {
				i=Points.getRoughThreshold(tmpPoints,lastFrequency,true,0,true);
				if (i!=-1) threshold=Math.max(threshold,tmpPoints[i].position);
			}
		}
		else {	
			minIntervalLength=Math.max( (int)Points.distanceQuantile(tmpPoints,lastFrequency,N_CONSECUTIVE_POINTS,true,QUANTILE),
										(int)((tmpPoints[lastFrequency].position-tmpPoints[0].position)/MIN_INTERVAL_LENGTH_FACTOR)
									  );
			minLocalMaxDistance=minIntervalLength<<1;
			if (!Points.areUniformlyDistributed(tmpPoints,0,lastFrequency,true,(tmpPoints[lastFrequency].position-tmpPoints[0].position)/Points.DEFAULT_NBINS)) {
				DensityEstimationTree.allocateMemory(lastFrequency+1);
				nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(tmpPoints,0,lastFrequency,minIntervalLength,minLocalMaxDistance,true,-1,-1,-1,false,false);				
				if (nLocalMaximumLeaves>0) {
					lastLeaf=DensityEstimationTree.lastLeaf;
					lvs=DensityEstimationTree.leaves;
					DensityEstimationTree.leaves=tmpLeaves;
					tmpLeaves=lvs;
					nHigh=Leaves.setHighLocalMax(tmpLeaves,lastLeaf,lastFrequency,false,false,tmpPoints,Points.mass(tmpPoints,0,lastFrequency),true,Math.POSITIVE_INFINITY);
					lvs=tmpLeaves;
					tmpLeaves=DensityEstimationTree.leaves;
					DensityEstimationTree.leaves=lvs;
					DensityEstimationTree.lastLeaf=lastLeaf;
					DensityEstimationTree.markRunsOfLocalMaximumLeaves(tmpPoints);
					threshold=Math.max(threshold,(int)DensityEstimationTree.separateRun(0,tmpPoints,0,lastFrequency,true));
				}
				DensityEstimationTree.deallocateMemory();
			}
		}	
		return (int)threshold;
	}
	
	
	/**
	 * Builds the transitively-reduced DAG that tells how to propagate the counts of the 
	 * prefix/suffix of each kernel, to its containing kernels in the kernel graph. The 
	 * start (respectively, end) position of every kernel $i$ is a node $2i$ 
	 * (respectively, $2i+1$), and there is an arc between end $p$ of kernel $i$ and end 
	 * $q$ of kernel $j$ iff: (1) there is an arc $(j,i)$ in the kernel graph; (2) there 
	 * is no insertion of $p$ into $j$; (3) $q$ is the end of $j$ that corresponds to $p$ 
	 * according to the orientation of arc $(j,i)$.
	 *
	 * Remark: transitive reduction is performed just to simplify the DAG, but it is not 
	 * strictly necessary (to avoid propagating the same counts multiple times to the same
	 * kernel, transitive reduction is not enough).
	 *
	 * Remark: for simplicity, the procedure assumes that insertions of both the start and
	 * the end of a kernel $i$ into a kernel $j$ (in the same orientation) mean that $i$
	 * is a substring of $j$, and thus no count is propagated. However, $i$ might also 
	 * occur twice in $j$ (in the same orientation), once as a prefix and once as a 
	 * suffix, and in this case both counts should be propagated.
	 *
	 * Remark: the DAG is represented by matrices $outNeighbors_ps$ (every row is sorted) 
	 * and $inNeighbors_ps$ (rows are not necessarily sorted).
	 */
	public static final void buildFrequencyGraph_prefixSuffix() {
		int i, j, k, p;
		int from, to, neighbor, nComponents, isDAG;
		boolean[] isActive;
		int[] lastInNeighborPrime, components, minimalVertices, nVerticesPerComponent, longestPath;
		int[] sorted2original, original2sorted;
		
		// Building the graph
		inNeighbors_ps = new int[nPathKernels<<1][0];
		lastInNeighbor_ps = new int[nPathKernels<<1];
		Math.set(lastInNeighbor_ps,(nPathKernels<<1)-1,-1);
		outNeighbors_ps = new int[nPathKernels<<1][0];
		lastOutNeighbor_ps = new int[nPathKernels<<1];
		Math.set(lastOutNeighbor_ps,(nPathKernels<<1)-1,-1);
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) continue;
			for (j=0; j<=lastInNeighbor[i]; j++) {
				neighbor=inNeighbors[i][j];
				if (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0)) continue;
				p=j<<2;
				from=i<<1;  // Start of kernel $i$.
				if (inNeighborFlags[i][p+1]==1 && (inNeighborFlags[i][p+0]==1 || inNeighborFlags[i][p+0]==3)) {
					to=neighbor<<1;
					addArc(from,to,inNeighbors_ps,lastInNeighbor_ps,outNeighbors_ps,lastOutNeighbor_ps);
				}
				if (inNeighborFlags[i][p+3]==1 && (inNeighborFlags[i][p+2]==1 || inNeighborFlags[i][p+2]==3)) {
   					to=(neighbor<<1)+1;
   					addArc(from,to,inNeighbors_ps,lastInNeighbor_ps,outNeighbors_ps,lastOutNeighbor_ps);
				}
				from=(i<<1)+1;  // End of kernel $i$.
				if (inNeighborFlags[i][p+1]==1 && (inNeighborFlags[i][p+0]==1 || inNeighborFlags[i][p+0]==2)) {
					to=(neighbor<<1)+1;
					addArc(from,to,inNeighbors_ps,lastInNeighbor_ps,outNeighbors_ps,lastOutNeighbor_ps);
				}
				if (inNeighborFlags[i][p+3]==1 && (inNeighborFlags[i][p+2]==1 || inNeighborFlags[i][p+2]==2)) {
   					to=neighbor<<1;
   					addArc(from,to,inNeighbors_ps,lastInNeighbor_ps,outNeighbors_ps,lastOutNeighbor_ps);
				}
			}
		}
		
		// Marking redundant arcs
		if (tmpArray1==null || tmpArray1.length<nPathKernels<<1) tmpArray1 = new int[nPathKernels<<1];
		if (tmpArray2==null || tmpArray2.length<nPathKernels<<1) tmpArray2 = new int[nPathKernels<<1];
		if (tmpArray3==null || tmpArray3.length<nPathKernels<<1) tmpArray3 = new int[nPathKernels<<1];
		if (tmpArray4==null || tmpArray4.length<nPathKernels<<1) tmpArray4 = new int[nPathKernels<<1];
		sorted2original = new int[nPathKernels<<1]; 
		original2sorted = new int[nPathKernels<<1];
		for (i=0; i<nPathKernels<<1; i++) lastInNeighbor_ps[i]++;
		for (i=0; i<nPathKernels<<1; i++) lastOutNeighbor_ps[i]++;
		components=tmpArray1;
		nComponents=DAG.getConnectedComponents(nPathKernels<<1,inNeighbors_ps,lastInNeighbor_ps,outNeighbors_ps,lastOutNeighbor_ps,components,IntervalGraph.stack);
		if (nComponents<=0 || nComponents>nPathKernels<<1) {
			IO.printCriticalErr("buildFrequencyGraph_prefixSuffix> ERROR: wrong number of connected components.");
			System.exit(1);
		}
		lastInNeighborPrime=tmpArray2;
		System.arraycopy(lastInNeighbor_ps,0,lastInNeighborPrime,0,nPathKernels<<1);  // Since the following procedure modifies its $lastInNeighbor$ argument.
		minimalVertices=tmpArray3; nVerticesPerComponent=tmpArray4;
		isDAG=DAG.topologicalSort(nPathKernels<<1,inNeighbors_ps,lastInNeighborPrime,outNeighbors_ps,lastOutNeighbor_ps,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
		if (isDAG!=0) {
			System.err.println("buildFrequencyGraph_prefixSuffix> ERROR: the graph contains a cycle that involves node "+(isDAG-1));
			System.exit(1);
		}
		longestPath=tmpArray2; isActive = new boolean[nPathKernels<<1];
		DAG.transitiveReduction(nPathKernels<<1,inNeighbors_ps,lastInNeighbor_ps,outNeighbors_ps,lastOutNeighbor_ps,sorted2original,original2sorted,components,isActive,longestPath);
		for (i=0; i<nPathKernels<<1; i++) lastInNeighbor_ps[i]--;
		for (i=0; i<nPathKernels<<1; i++) lastOutNeighbor_ps[i]--;
		
		// Removing redundant arcs
		for (i=0; i<nPathKernels<<1; i++) {
			if (pathKernelLengths[i>>1]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i>>1]!=0)) continue;
			k=-1;
			for (j=0; j<=lastInNeighbor_ps[i]; j++) {
				if (inNeighbors_ps[i][j]>=0) continue;
				k++;
				inNeighbors_ps[i][k]=-1-inNeighbors_ps[i][j];
			}
			lastInNeighbor_ps[i]=k;
		}
		for (i=0; i<nPathKernels<<1; i++) {
			if (pathKernelLengths[i>>1]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i>>1]!=0)) continue;
			k=-1;
			for (j=0; j<=lastOutNeighbor_ps[i]; j++) {
				if (outNeighbors_ps[i][j]>=0) continue;
				k++;
				outNeighbors_ps[i][k]=-1-outNeighbors_ps[i][j];
			}
			lastOutNeighbor_ps[i]=k;
		}
		
		// Sorting $outNeighbors_ps$.
		for (i=0; i<nPathKernels<<1; i++) {
			if (pathKernelLengths[i>>1]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[i>>1]==0) && lastOutNeighbor_ps[i]>0) Arrays.sort(outNeighbors_ps[i],0,lastOutNeighbor_ps[i]+1);
		}
	}
	
	
	/**
	 * Builds the transitively-reduced DAG that tells how to propagate the counts of the 
	 * full occurrences of a kernel, to its containing kernels in the kernel graph. The 
	 * DAG before transitive reduction coincides with the subgraph of the kernel graph
	 * induced by identity arcs.
	 *
	 * Remark: transitive reduction is done just to simplify the DAG, but it not strictly 
	 * necessary (it is not enough to avoid propagating the same counts multiple times to
	 * the same kernel).
	 *
	 * Remark: the DAG is represented by matrices $outNeighbors_f,inNeighbors_f$ (their
	 * rows are not necessarily sorted).
	 */
	public static final void buildFrequencyGraph_full() {
		int i, j, k, p;
		int from, to, neighbor, nComponents, isDAG;
		boolean[] isActive;
		int[] lastInNeighborPrime, components, minimalVertices, nVerticesPerComponent, longestPath;
		int[] sorted2original, original2sorted;
		
		// Building the graph
		inNeighbors_f = new int[nPathKernels][0];
		lastInNeighbor_f = new int[nPathKernels];
		Math.set(lastInNeighbor_f,nPathKernels-1,-1);
		outNeighbors_f = new int[nPathKernels][0];
		lastOutNeighbor_f = new int[nPathKernels];
		Math.set(lastOutNeighbor_f,nPathKernels-1,-1);
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) continue;
			for (j=0; j<=lastInNeighbor[i]; j++) {
				neighbor=inNeighbors[i][j];
				if (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0)) continue;
				p=j<<2; from=i; to=neighbor;
				if ( (inNeighborFlags[i][p+1]==1 && inNeighborFlags[i][p+0]==1) ||
					 (inNeighborFlags[i][p+3]==1 && inNeighborFlags[i][p+2]==1)
				   ) addArc(from,to,inNeighbors_f,lastInNeighbor_f,outNeighbors_f,lastOutNeighbor_f);
			}
		}
		
		// Marking redundant arcs
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		if (tmpArray2==null || tmpArray2.length<nPathKernels) tmpArray2 = new int[nPathKernels];
		if (tmpArray3==null || tmpArray3.length<nPathKernels) tmpArray3 = new int[nPathKernels];
		if (tmpArray4==null || tmpArray4.length<nPathKernels) tmpArray4 = new int[nPathKernels];
		sorted2original = new int[nPathKernels]; 
		original2sorted = new int[nPathKernels];
		for (i=0; i<nPathKernels; i++) lastInNeighbor_f[i]++;
		for (i=0; i<nPathKernels; i++) lastOutNeighbor_f[i]++;
		components=tmpArray1;
		nComponents=DAG.getConnectedComponents(nPathKernels,inNeighbors_f,lastInNeighbor_f,outNeighbors_f,lastOutNeighbor_f,components,IntervalGraph.stack);
		if (nComponents<=0 || nComponents>nPathKernels) {
			IO.printCriticalErr("buildFrequencyGraph_full> ERROR: wrong number of connected components.");
			System.exit(1);
		}
		lastInNeighborPrime=tmpArray2;
		System.arraycopy(lastInNeighbor_f,0,lastInNeighborPrime,0,nPathKernels);  // Since the following procedure modifies its $lastInNeighbor$ argument.
		minimalVertices=tmpArray3; nVerticesPerComponent=tmpArray4;
		isDAG=DAG.topologicalSort(nPathKernels,inNeighbors_f,lastInNeighborPrime,outNeighbors_f,lastOutNeighbor_f,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
		if (isDAG!=0) {
			System.err.println("buildFrequencyGraph_full> ERROR: the graph contains a cycle that involves node "+(isDAG-1));
			System.exit(1);
		}
		longestPath=tmpArray2; isActive = new boolean[nPathKernels];
		DAG.transitiveReduction(nPathKernels,inNeighbors_f,lastInNeighbor_f,outNeighbors_f,lastOutNeighbor_f,sorted2original,original2sorted,components,isActive,longestPath);
		for (i=0; i<nPathKernels; i++) lastInNeighbor_f[i]--;
		for (i=0; i<nPathKernels; i++) lastOutNeighbor_f[i]--;
		
		// Removing redundant arcs
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) continue;
			k=-1;
			for (j=0; j<=lastInNeighbor_f[i]; j++) {
				if (inNeighbors_f[i][j]>=0) continue;
				k++;
				inNeighbors_f[i][k]=-1-inNeighbors_f[i][j];
			}
			lastInNeighbor_f[i]=k;
		}
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) continue;
			k=-1;
			for (j=0; j<=lastOutNeighbor_f[i]; j++) {
				if (outNeighbors_f[i][j]>=0) continue;
				k++;
				outNeighbors_f[i][k]=-1-outNeighbors_f[i][j];
			}
			lastOutNeighbor_f[i]=k;
		}
	}
	
	
	private static final void addArc(int from, int to, int[][] in, int[] lastIn, int[][] out, int[] lastOut) {
		final int GROWTH_RATE = 10;  // Arbitrary
		
		if (out[from]==null || out[from].length==0) out[from] = new int[lastOut[from]+1+GROWTH_RATE];
		else if (out[from].length<=lastOut[from]+1) {
			int[] newArray = new int[lastOut[from]+1+GROWTH_RATE];
			System.arraycopy(out[from],0,newArray,0,lastOut[from]+1);
			out[from]=newArray;
		}
		out[from][++lastOut[from]]=to;
		if (in[to]==null || in[to].length==0) in[to] = new int[lastIn[to]+1+GROWTH_RATE];
		else if (in[to].length<=lastIn[to]+1) {
			int[] newArray = new int[lastIn[to]+GROWTH_RATE];
			System.arraycopy(in[to],0,newArray,0,lastIn[to]+1);
			in[to]=newArray;
		}
		in[to][++lastIn[to]]=from;
	}
	
	
	public static final void deallocateFrequencyGraphs() {
		int i;
		
		for (i=0; i<nPathKernels; i++) inNeighbors_ps[i]=null;
		inNeighbors_ps=null; lastInNeighbor_ps=null;
		for (i=0; i<nPathKernels; i++) outNeighbors_ps[i]=null;
		outNeighbors_ps=null; lastOutNeighbor_ps=null;
		for (i=0; i<nPathKernels; i++) inNeighbors_f[i]=null;
		inNeighbors_f=null; lastInNeighbor_f=null;
		for (i=0; i<nPathKernels; i++) outNeighbors_f[i]=null;
		outNeighbors_f=null; lastOutNeighbor_f=null;
	}
	
	
	/**
	 * @param fromVertex a vertex in the $*_ps$ frequency DAG;
	 * @param visited temporary space, assumed to be all FALSE; the procedure resets it to
	 * this state at the end;
	 * @return the sum of the $kernelFrequency$ entries of all reachable ancestors of 
	 * $fromVertex$ in the DAG ($fromVertex$ included); traversal stops at vertices that
	 * correspond to a kernel $k$ with $labels[k]!=0$.
	 */
	private static final int getNewFrequency_ps(int fromVertex, byte[] labels, int[] stack, boolean[] visited) {
		int i;
		int vertex, neighbor, sum, top;
		
		sum=kernelFrequency[fromVertex>>1][fromVertex%2==0?1:2];
		stack[0]=fromVertex; top=0;
		while (top>=0) {
			vertex=stack[top--];
			for (i=0; i<=lastInNeighbor_ps[vertex]; i++) {
				neighbor=inNeighbors_ps[vertex][i];
				if (visited[neighbor] || labels[neighbor>>1]!=0) continue;				
				sum+=kernelFrequency[neighbor>>1][neighbor%2==0?1:2];
				visited[neighbor]=true;
				stack[++top]=neighbor;
			}
		}
		
		// Cleaning up for next call
		stack[0]=fromVertex; top=0;
		while (top>=0) {
			vertex=stack[top--]; 
			for (i=0; i<=lastInNeighbor_ps[vertex]; i++) {
				neighbor=inNeighbors_ps[vertex][i];
				if (!visited[neighbor]) continue;
				visited[neighbor]=false;
				stack[++top]=neighbor;
			}
		}
		return sum;
	}
	
	
	/**
	 * Remark: substring counts are not updated. Note that e.g. the prefix counts of a 
	 * kernel might have to be added to the substring counts of a containing kernel, and 
	 * to do so correctly we should create another DAG.
	 * 
	 * @param fromVertex a vertex in the $*_f$ frequency DAG;
	 * @param visited temporary space, assumed to be all FALSE; the procedure resets it to
	 * this state at the end;
	 * @return the sum of the $kernelFrequency$ entries of all reachable descendants of 
	 * $fromVertex$ in the DAG ($fromVertex$ included); traversal stops at vertices that
	 * correspond to a kernel $k$ with $labels[k]!=0$.
	 */
	private static final int getNewFrequency_f(int fromVertex, byte[] labels, int[] stack, boolean[] visited) {
		int i;
		int vertex, neighbor, sum, top;
		
		sum=kernelFrequency[fromVertex][0];
		stack[0]=fromVertex; top=0;
		while (top>=0) {
			vertex=stack[top--];
			for (i=0; i<=lastInNeighbor_f[vertex]; i++) {
				neighbor=inNeighbors_f[vertex][i];
				if (visited[neighbor] || labels[neighbor>>1]!=0) continue;
				sum+=kernelFrequency[neighbor][0];
				visited[neighbor]=true;
				stack[++top]=neighbor;
			}
		}
		
		// Cleaning up for next call
		stack[0]=fromVertex; top=0;
		while (top>=0) {
			vertex=stack[top--];
			for (i=0; i<=lastInNeighbor_f[vertex]; i++) {
				neighbor=inNeighbors_f[vertex][i];
				if (!visited[neighbor]) continue;
				visited[neighbor]=false;
				stack[++top]=neighbor;
			}
		}
		return sum;
	}
	
	
	/**
	 * @param startOrEnd the start (TRUE) or end (FALSE) of $kernel$ that has high 
	 * frequency (the procedure assumes that there is just one end with high frequency);
	 * @return TRUE iff $kernel$ is contained in at least one of its in-neighbors, so that 
	 * the end of $kernel$ marked by $startOrEnd$ is aligned with an end of the in-
	 * neighbor. This might be an overly-aggressive criterion to conclude that $kernel$ is
	 * a fragment of another repeat.
	 */
	private static final boolean markRepeats_isPrefixOrSuffix(int kernel, boolean startOrEnd) {
		int i, p;
		int neighbor;
		
		if (startOrEnd) {
			for (i=0; i<=lastInNeighbor[kernel]; i++) {
				neighbor=inNeighbors[kernel][i];
				if (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0)) continue;
				p=i<<2;
				if ( (inNeighborFlags[kernel][p+1]==1 && (inNeighborFlags[kernel][p+0]==1 || inNeighborFlags[kernel][p+0]==3)) ||
					 (inNeighborFlags[kernel][p+3]==1 && (inNeighborFlags[kernel][p+2]==1 || inNeighborFlags[kernel][p+2]==3))
				   ) return true;
			}
		}
		else {
			for (i=0; i<=lastInNeighbor[kernel]; i++) {
				neighbor=inNeighbors[kernel][i];
				if (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0)) continue;
				p=i<<2;
				if ( (inNeighborFlags[kernel][p+1]==1 && (inNeighborFlags[kernel][p+0]==1 || inNeighborFlags[kernel][p+0]==2)) ||
					 (inNeighborFlags[kernel][p+3]==1 && (inNeighborFlags[kernel][p+2]==1 || inNeighborFlags[kernel][p+2]==2))
				   ) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * @return TRUE iff $kernel$ is contained without insertion in at least one of its 
	 * in-neighbors.
	 */
	private static final boolean markRepeats_isIdentical(int kernel) {
		int i, p;
		int neighbor;
		
		for (i=0; i<=lastInNeighbor[kernel]; i++) {
			neighbor=inNeighbors[kernel][i];
			if (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0)) continue;
			p=i<<2;
			if ( (inNeighborFlags[kernel][p+1]==1 && inNeighborFlags[kernel][p+3]==0 && inNeighborFlags[kernel][p+0]==1) || 
				 (inNeighborFlags[kernel][p+1]==0 && inNeighborFlags[kernel][p+3]==1 && inNeighborFlags[kernel][p+2]==1) ||
				 (inNeighborFlags[kernel][p+1]==1 && inNeighborFlags[kernel][p+3]==1 && inNeighborFlags[kernel][p+0]==1 && inNeighborFlags[kernel][p+2]==1)
			   ) return true;
		}
		return false;
	}
	
	
	/**
	 * Decides which kernel paths represent the full copy of a repeat, and which are just
	 * fragments of a repeat, by traversing the kernel graph bottom-up and inspecting the 
	 * number of prefix, suffix, and full occurrences of each kernel path. The procedure 
	 * assumes that fragmentation occurs at random positions inside a repeat, so if a 
	 * kernel path is a fragment of a repeat, it should have few full occurrences in the 
	 * read set, and few occurrences of the end(s) of the kernel path that corresponds to 
	 * the random event (prefix/suffix). If the kernel is a substring of a repeat, both 
	 * ends should be rare. If the kernel is e.g. the prefix of a repeat, only its suffix 
	 * should be rare. If a kernel is marked as a fragment, its frequencies are propagated
	 * to its in-neighbors.
	 *
	 * A kernel path might also represent an erroneous concatenation of multiple repeats.
	 * This can happen because of issues in factorization, and because the element of a 
	 * cluster in the interval graph are detected using domain-independent graph 
	 * clustering: many objects at the boundary of clusters might be wrong concatenations 
	 * of repeats, and some of them might be assigned to a specific cluster (rather than
	 * being removed from all clusters) in Step 2.
	 *
	 * Remark: it might happen that a kernel path is a substring of a rotation of a 
	 * repeat. This could happen e.g. if two full copies of the repeat are adjacent in the
	 * same read, and an interval around their boundary was produced by mistake during 
	 * factorization. A random rotation is not likely contained in any other kernel, and 
	 * it is likely to have low frequency at both ends. If it has high frequency at one 
	 * end, it is handled like a wrong concatenation (see below).
	 *
	 * Remark: if a kernel path is exactly the concatenation of several repeats, it should
	 * be handled correctly, since only intervals that are not tagged with its component 
	 * repeats contribute to its frequencies (including its full frequency: see 
	 * $getFrequency()$ for details). I.e. if it is not itself a repeat, then it should be
	 * handled like any other wrong concatenation, since the prefix/suffix frequencies of 
	 * its first/last repeat should not be propagated to it.
	 *
	 * Call "red" a kernel path with both high frequency start and end, and "blue" a 
	 * kernel with high frequency in just its start or end. Blue/red kernel paths might 
	 * descend from red kernel paths in the kernel graph: this might happen because of the
	 * choice of $minFrequency$, or because they are real repeats contained in other 
	 * repeats. The procedure removes red kernels with a red ancestor iff they always 
	 * occur in the read set as contained in one of their ancestors. Blue kernels are 
	 * removed if either: (1) they have a red ancestor, and they always occur contained in
	 * one of their ancestors; or (2) they have no red ancestor, and their high-frequency 
	 * end aligns with an end of a red descendant. I.e. the procedure assumes that a blue 
	 * kernel whose end aligns with the end of a red descendant is a wrong concatenation 
	 * of the descendant.
	 *
	 * Remark: a blue kernel might have its high-frequency end that aligns to a red, but
	 * this alignment might not be detected and it might instead be mistaken for an 
	 * insertion, due to the distance thresholds used to decide insertions. E.g. an 
	 * interval of the blue kernel whose end aligns to an interval of the red kernel, 
	 * might not have stored the tags of the blue kernel in its $pathsWith*$ array, 
	 * because in the interval graph it is connected by an identity edge to an interval 
	 * that is a substring of the blue kernel; combining the relative offsets of such 
	 * intervals WRT the kernel one would get a small enough offset, but this is 
	 * impossible to do in practice, since identity edges in the interval graph do not 
	 * store overhangs). The same problem might stop the propagation of end counts across
	 * the kernel graph.
	 *
	 * A removed kernel is assigned to all its surviving ancestors in the kernel graph. 
	 * The procedure resets the $kernels$ array of every interval graph node to contain 
	 * only the surviving kernels. See $markRepeats_removeRedundantKernels()$ for details.
	 *
	 * Remark: the procedure might not mark any path kernel as an independent repeat.
	 *
	 * Remark: the kernel graph might consist of multiple connected components, and some 
	 * of them might have no path kernel marked as an independent repeat.
	 *
	 * Remark: the procedure initializes $sorted2original,original2sorted$.
	 *
	 * @param minFrequency counts are considered rare iff they are smaller than this;
	 * @param labels output array with a number in [-1..5] assigned to every path kernel
	 * (see the code for details);
	 * @param newFrequencies temporary space, of size at least equal to the size of 
	 * $kernelFrequency$; at the end of the procedure, the columns of full, prefix and 
	 * suffix counts contain the updated values after bottom-up traversal (substring 
	 * counts are not updated);
	 * @return the number of repeats with a positive value in $labels$ (might be zero).
	 */
	public static final int markRepeats(int minFrequency, byte[] labels, int[][] newFrequencies) {
		final int REDUNDANT_THRESHOLD = IO.coverage<<1;  // Arbitrary
		final int GROWTH_RATE = 8;  // Arbitrary
		boolean highFrequencyStart, highFrequencyEnd;
		int i, j, p, q;
		int kernel, vertex, last, neighbor, out, last4;
		int nKernels25, from, read, nonRedundant;
		IntervalGraph.Node tmpNode;
		boolean[] hasAncestor, endAligns, selectedLabels, oneBidirectedNode;
		int[] kernels25, lastAncestor, isRedundant, ancestors_lastStart, ancestors_lastEnd;
		int[][] ancestors, ancestors_start, ancestors_end;
		byte[][] ancestorOrientations;
		
		// Initializing arrays
		if (tmpArray1==null || tmpArray1.length<(nPathKernels<<1)) tmpArray1 = new int[nPathKernels<<1];
		if (tmpArray2==null || tmpArray2.length<(nPathKernels<<1)) tmpArray2 = new int[nPathKernels<<1];
		if (tmpArray3==null || tmpArray3.length<(nPathKernels<<1)) tmpArray3 = new int[nPathKernels<<1];
		if (tmpArray4==null || tmpArray4.length<(nPathKernels<<1)) tmpArray4 = new int[nPathKernels<<1];
		if (tmpArray5==null || tmpArray5.length<(nPathKernels<<1)) tmpArray5 = new int[nPathKernels<<1];
		if (tmpByte1==null || tmpByte1.length<(nPathKernels<<1)) tmpByte1 = new byte[nPathKernels<<1];
		if (tmpByte2==null || tmpByte2.length<(nPathKernels<<1)) tmpByte2 = new byte[nPathKernels<<1];
		if (tmpByte3==null || tmpByte3.length<(nPathKernels<<1)) tmpByte3 = new byte[nPathKernels<<1];
		if (tmpBoolean1==null || tmpBoolean1.length<(nPathKernels<<1)) tmpBoolean1 = new boolean[nPathKernels<<1];
		if (tmpBoolean2==null || tmpBoolean2.length<(nPathKernels<<1)) tmpBoolean2 = new boolean[nPathKernels<<1];
		selectedLabels = new boolean[6];
		Math.set(labels,nPathKernels-1,(byte)(-1));
		
		sortKernelGraph();
		
		// Assigning the following tags to $labels$:
		// -1: cyclic, or non-cyclic and periodic;
		//  0: neither high frequency start nor end;
		//  1: both high frequency start and end;
		//  2: high-frequency in just start or end.
		Math.set(tmpBoolean1,tmpBoolean1.length-1,false);
		out=0;
		for (i=0; i<nPathKernels; i++) {
			kernel=sorted2original[i];
			if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
			newFrequencies[kernel][0]=getNewFrequency_f(kernel,labels,tmpArray1,tmpBoolean1);
			newFrequencies[kernel][1]=getNewFrequency_ps(kernel<<1,labels,tmpArray1,tmpBoolean1);
			newFrequencies[kernel][2]=getNewFrequency_ps((kernel<<1)+1,labels,tmpArray1,tmpBoolean1);
			highFrequencyStart=newFrequencies[kernel][0]+newFrequencies[kernel][1]>=minFrequency;
			highFrequencyEnd=newFrequencies[kernel][0]+newFrequencies[kernel][2]>=minFrequency;
			if (highFrequencyStart && highFrequencyEnd) {
				if (markRepeats_isIdentical(kernel)) labels[kernel]=0;
				else {
					labels[kernel]=1;
					out++;
				}
			}
			else if (highFrequencyStart) {
				if (markRepeats_isPrefixOrSuffix(kernel,true)) labels[kernel]=0;
				else {
					labels[kernel]=2;
					out++;
				}
			}
			else if (highFrequencyEnd) {
				if (markRepeats_isPrefixOrSuffix(kernel,false)) labels[kernel]=0;
				else {
					labels[kernel]=2;
					out++;
				}
			}
			else labels[kernel]=0;
		}
		
		
if (IO.CONSISTENCY_CHECKS) {
	for (int x=0; x<nPathKernels; x++) {
		if (labels[x]<0 && !(pathKernelLengths[x]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[x]!=0))) {
			System.err.println("markrepeats> 0  ERROR: inconsistent labels value:  "+x+" -> "+labels[x]+", "+pathKernelLengths[x]+", "+(pathKernelPeriodic==null?"null":pathKernelPeriodic[x]));
			for (int y=0; y<nPathKernels; y++) System.err.println(y+" :: "+labels[y]+", "+pathKernelLengths[y]+", "+(pathKernelPeriodic==null?"null":pathKernelPeriodic[y]));
			System.exit(1);
		}
	}
}
		
		System.err.println("markRepeats> Kernels with a high-frequency end: "+out);
		if (GRAPH_FILE!=null) try { printKernelGraph(GRAPH_FILE+".kernelGraph-1.dot",labels,newFrequencies,null); } catch(IOException e) { }
		if (out==0) {
			removeAllKernelTags();
			return 0;
		}
		
if (IO.CONSISTENCY_CHECKS) {
	for (int x=0; x<nPathKernels; x++) {
		if (labels[x]<0 && !(pathKernelLengths[x]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[x]!=0))) {
			System.err.println("markrepeats> 1  ERROR: inconsistent labels value:  "+x+" -> "+labels[x]+", "+pathKernelLengths[x]+", "+(pathKernelPeriodic==null?"null":pathKernelPeriodic[x]));
			for (int y=0; y<nPathKernels; y++) System.err.println(y+" :: "+labels[y]+", "+pathKernelLengths[y]+", "+(pathKernelPeriodic==null?"null":pathKernelPeriodic[y]));
			System.exit(1);
		}
	}
}
		
		// Assigning the following tags to $labels$. We call "red" a node with both high
		// frequency start and end, and "blue" a node with high frequency in just its
		// start or end.
		// 
		// -1: cyclic, or non-cyclic and periodic;
	 	//  0: neither red nor blue;
		//  1: red, without any red ancestor; 
		//  2: red, with a red ancestor;
		//  3: blue, without any red ancestor, and its high-frequency end does not align
		//  to the end of a red;
		//  4: blue, without any red ancestor, and its high-frequency end aligns to the
		//  end of a red;
		//  5: blue, with a red ancestor.
		//
		hasAncestor=tmpBoolean1;
		Math.set(hasAncestor,nPathKernels-1,false);
		for (i=nPathKernels-1; i>=0; i--) {
			kernel=sorted2original[i];
			if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0) || !(hasAncestor[kernel] || labels[kernel]==1)) continue;
			for (j=0; j<=lastOutNeighbor[kernel]; j++) hasAncestor[outNeighbors[kernel][j]]=true;
		}
		endAligns=tmpBoolean2;
		Math.set(endAligns,(nPathKernels<<1)-1,false);
		for (i=0; i<nPathKernels; i++) {
			kernel=sorted2original[i];
			if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;			
			vertex=kernel<<1;
			if (endAligns[vertex] || labels[kernel]==1) {
				for (j=0; j<=lastOutNeighbor_ps[vertex]; j++) endAligns[outNeighbors_ps[vertex][j]]=true;
			}
			vertex+=1;
			if (endAligns[vertex] || labels[kernel]==1) {
				for (j=0; j<=lastOutNeighbor_ps[vertex]; j++) endAligns[outNeighbors_ps[vertex][j]]=true;
			}
		}
		nKernels25=0;
		for (i=0; i<nPathKernels; i++) {
			kernel=sorted2original[i];
			if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
			if (labels[kernel]==1) {
				if (hasAncestor[kernel]) {
					labels[kernel]=2;
					nKernels25++;
				}
			}
			else if (labels[kernel]==2) {
				if (hasAncestor[kernel]) {
					labels[kernel]=5;
					nKernels25++;
				}
				else {
					highFrequencyStart=newFrequencies[kernel][0]+newFrequencies[kernel][1]>=minFrequency;
					highFrequencyEnd=newFrequencies[kernel][0]+newFrequencies[kernel][2]>=minFrequency;
					if ((highFrequencyStart && !endAligns[kernel<<1]) || (highFrequencyEnd && !endAligns[(kernel<<1)+1])) labels[kernel]=3;
					else labels[kernel]=4;
				}
			}
		}
		kernels25 = new int[nKernels25];
		if (nKernels25>0) {
			nKernels25=0;
			for (i=0; i<nPathKernels; i++) {
				kernel=sorted2original[i];
				if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
				if (labels[kernel]==2 || labels[kernel]==5) kernels25[nKernels25++]=kernel;
			}
			if (nKernels25>1) Arrays.sort(kernels25,0,nKernels25);
		}
		System.err.println("markRepeats> Kernels of type 2 or 5: "+nKernels25);
		if (GRAPH_FILE!=null) try { printKernelGraph(GRAPH_FILE+".kernelGraph-2.dot",labels,newFrequencies,null); } catch(IOException e) { }

		
if (IO.CONSISTENCY_CHECKS) {
	for (int x=0; x<nPathKernels; x++) {
		if (labels[x]<0 && !(pathKernelLengths[x]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[x]!=0))) {
			System.err.println("markrepeats> 2  ERROR: inconsistent labels value:  "+x+" -> "+labels[x]+", "+pathKernelLengths[x]+", "+(pathKernelPeriodic==null?"null":pathKernelPeriodic[x]));
			for (int y=0; y<nPathKernels; y++) System.err.println(y+" :: "+labels[y]+", "+pathKernelLengths[y]+", "+(pathKernelPeriodic==null?"null":pathKernelPeriodic[y]));
			System.exit(1);
		}
	}
}
		
		
		// Propagating the tag of a type-{1,2,3,5} kernel to all type-{1,2,3,5} kernels
		// that descend from it in the kernel graph (=building the $ancestors$ matrix).
		kernel2NewKernels = new int[nPathKernels][GROWTH_RATE];
		lastNewKernel = new int[nPathKernels];
		kernel2NewKernels_orientation = new byte[nPathKernels][GROWTH_RATE];
		kernel2NewKernels_start = new int[nPathKernels][0];
		kernel2NewKernels_lastStart = new int[nPathKernels];
		kernel2NewKernels_end = new int[nPathKernels][0];
		kernel2NewKernels_lastEnd = new int[nPathKernels];
		if (nKernels25>0) {
			selectedLabels[0]=false; selectedLabels[1]=true; selectedLabels[2]=true; 
			selectedLabels[3]=true; selectedLabels[4]=false; selectedLabels[5]=true;
			markRepeats_propagateKernels(selectedLabels,true,sorted2original,labels,kernel2NewKernels,lastNewKernel,kernel2NewKernels_orientation,true,kernel2NewKernels_start,kernel2NewKernels_lastStart,kernel2NewKernels_end,kernel2NewKernels_lastEnd,GROWTH_RATE);
			ancestors = new int[nKernels25][GROWTH_RATE];
			ancestorOrientations = new byte[nKernels25][GROWTH_RATE];
			lastAncestor = new int[nKernels25];
			Math.set(lastAncestor,nKernels25-1,-1);
			ancestors_start = new int[nKernels25][GROWTH_RATE];
			ancestors_lastStart = new int[nKernels25];
			Math.set(ancestors_lastStart,nKernels25-1,-1);
			ancestors_end = new int[nKernels25][GROWTH_RATE];
			ancestors_lastEnd = new int[nKernels25];
			Math.set(ancestors_lastEnd,nKernels25-1,-1);
			j=0;
			for (i=0; i<nPathKernels; i++) {
				if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) continue;
				if (labels[i]!=2 && labels[i]!=5) continue;
				// Removing $i$ from the ancestors of $i$.
				p=Arrays.binarySearch(kernel2NewKernels[i],0,lastNewKernel[i]+1,i);
				if (p>=0) { 
					p++;
					while (p<=lastNewKernel[i]) {
						kernel2NewKernels[i][p-1]=kernel2NewKernels[i][p];
						kernel2NewKernels_orientation[i][p-1]=kernel2NewKernels_orientation[i][p];
						p++;
					}
					lastNewKernel[i]--;
				}
				q=-1;
				for (p=0; p<=kernel2NewKernels_lastStart[i]; p++) {
					if (kernel2NewKernels_start[i][p]!=i && kernel2NewKernels_start[i][p]!=-1-i) kernel2NewKernels_start[i][++q]=kernel2NewKernels_start[i][p];
				}
				kernel2NewKernels_lastStart[i]=q;
				q=-1;
				for (p=0; p<=kernel2NewKernels_lastEnd[i]; p++) {
					if (kernel2NewKernels_end[i][p]!=i && kernel2NewKernels_end[i][p]!=-1-i) kernel2NewKernels_end[i][++q]=kernel2NewKernels_end[i][p];
				}
				kernel2NewKernels_lastEnd[i]=q;
				// Copying to ancestor arrays
				if (lastNewKernel[i]>=ancestors[j].length) ancestors[j] = new int[lastNewKernel[i]+1];
				System.arraycopy(kernel2NewKernels[i],0,ancestors[j],0,lastNewKernel[i]+1);
				if (lastNewKernel[i]>=ancestorOrientations[j].length) ancestorOrientations[j] = new byte[lastNewKernel[i]+1];
				System.arraycopy(kernel2NewKernels_orientation[i],0,ancestorOrientations[j],0,lastNewKernel[i]+1);
				lastAncestor[j]=lastNewKernel[i];
				if (kernel2NewKernels_lastStart[i]>=ancestors_start[j].length) ancestors_start[j] = new int[kernel2NewKernels_lastStart[i]+1];
				System.arraycopy(kernel2NewKernels_start[i],0,ancestors_start[j],0,kernel2NewKernels_lastStart[i]+1);
				ancestors_lastStart[j]=kernel2NewKernels_lastStart[i];
				if (kernel2NewKernels_lastEnd[i]>=ancestors_end[j].length) ancestors_end[j] = new int[kernel2NewKernels_lastEnd[i]+1];
				System.arraycopy(kernel2NewKernels_end[i],0,ancestors_end[j],0,kernel2NewKernels_lastEnd[i]+1);
				ancestors_lastEnd[j]=kernel2NewKernels_lastEnd[i];
				j++;
			}
		}
		else {
			ancestors=null; ancestorOrientations=null; lastAncestor=null;
			ancestors_start=null; ancestors_lastStart=null;
			ancestors_end=null; ancestors_lastEnd=null;
		}
		
if (IO.CONSISTENCY_CHECKS) {
	for (int x=0; x<nPathKernels; x++) {
		if (labels[x]<0 && !(pathKernelLengths[x]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[x]!=0))) {
			System.err.println("markrepeats> 3  ERROR: inconsistent labels value:  "+x+" -> "+labels[x]+", "+pathKernelLengths[x]+", "+(pathKernelPeriodic==null?"null":pathKernelPeriodic[x]));
			for (int y=0; y<nPathKernels; y++) System.err.println(y+" :: "+labels[y]+", "+pathKernelLengths[y]+", "+(pathKernelPeriodic==null?"null":pathKernelPeriodic[y]));
			System.exit(1);
		}
	}
}

		
		// Propagating the tag of a type-{1,2,3,5} kernel, to all kernels that descend
		// from it in the kernel graph and that are not themselves of type {1,2,3,5}, and
		// updating the $kernels$ field of interval graph nodes.
		// Remark: interval graph nodes assigned just to kernels with no type-{1,2,3,5}
		// ancestor in the kernel graph have no kernel tag after this.
		selectedLabels[0]=false; selectedLabels[1]=true; selectedLabels[2]=true; 
		selectedLabels[3]=true; selectedLabels[4]=false; selectedLabels[5]=true;
		markRepeats_propagateKernels(selectedLabels,false,sorted2original,labels,kernel2NewKernels,lastNewKernel,kernel2NewKernels_orientation,true,kernel2NewKernels_start,kernel2NewKernels_lastStart,kernel2NewKernels_end,kernel2NewKernels_lastEnd,GROWTH_RATE);
		for (i=0; i<IntervalGraph.nNodes; i++) markRepeats_updateKernels(IntervalGraph.nodesArray[i],true,kernel2NewKernels,lastNewKernel,kernel2NewKernels_orientation,kernel2NewKernels_start,kernel2NewKernels_lastStart,kernel2NewKernels_end,kernel2NewKernels_lastEnd,-1);
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) continue;
			if (labels[i]==4) {
				labels[i]=0;
				out--;
			}
		}
		
		// Detecting redundant type-{2,5} kernels, setting to zero their $labels$ field,
		// and removing redundant kernels from the $kernels$ field of interval graph
		// nodes.
		if (nKernels25==0) {
			if (IO.CONSISTENCY_CHECKS) markRepeats_checkKernelsWithZeroLabel(labels);
			return out;
		}
		Math.set(tmpArray2,nKernels25-1,0);
		isRedundant=tmpArray2;
		oneBidirectedNode=tmpBoolean1;
		markKernelsWithOneBidirectedNode(oneBidirectedNode);
		nonRedundant=0; from=0; read=IntervalGraph.nodesArray[0].read;
		for (i=1; i<IntervalGraph.nNodes; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			if (tmpNode.read!=read) {
				nonRedundant+=markRepeats_markRedundantKernels(from,i-1,kernels25,nKernels25-1,ancestors,lastAncestor,isRedundant,REDUNDANT_THRESHOLD,oneBidirectedNode);
				if (nonRedundant==nKernels25) break;
				from=i; read=tmpNode.read;
			}
		}
		if (nonRedundant<nKernels25) nonRedundant+=markRepeats_markRedundantKernels(from,IntervalGraph.nNodes-1,kernels25,nKernels25-1,ancestors,lastAncestor,isRedundant,REDUNDANT_THRESHOLD,oneBidirectedNode);
		if (nonRedundant<nKernels25) {
			for (i=0; i<nKernels25; i++) {
				if (isRedundant[i]<REDUNDANT_THRESHOLD) {
					labels[kernels25[i]]=0;
					out--;
				}
			}
		}
		System.err.println("markRepeats> Non-redundant type-25 kernels: "+nonRedundant+" out of "+nKernels25);
		if (GRAPH_FILE!=null) try { printKernelGraph(GRAPH_FILE+".kernelGraph-3.dot",labels,newFrequencies,null); } catch(IOException e) { }
		markRepeats_removeRedundantAncestors(kernels25,nKernels25-1,labels,ancestors,lastAncestor,ancestorOrientations,ancestors_start,ancestors_lastStart,ancestors_end,ancestors_lastEnd);
		if (nonRedundant<nKernels25) {
			for (i=0; i<IntervalGraph.nNodes; i++) markRepeats_removeRedundantKernels(IntervalGraph.nodesArray[i],kernels25,nKernels25-1,labels,ancestors,lastAncestor,ancestorOrientations,ancestors_start,ancestors_lastStart,ancestors_end,ancestors_lastEnd,-1);
		}
		if (IO.CONSISTENCY_CHECKS) markRepeats_checkKernelsWithZeroLabel(labels);
		return out;
	}
	
	
	/**
	 * Removes all non-cyclic and non-periodic kernel tags from every interval graph node.
	 */
	private static final void removeAllKernelTags() {
		int i, j, k;
		int kernel;
		IntervalGraph.Node tmpNode;
		
		for (i=0; i<IntervalGraph.nNodes; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			if (tmpNode.lastKernel==-1) continue;
			k=-1;
			for (j=0; j<=tmpNode.lastKernel; j++) {
				kernel=tmpNode.kernels[j];
				if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
					k++;
					tmpNode.kernels[k]=kernel;
					tmpNode.kernelOrientations[k]=tmpNode.kernelOrientations[j];
				}
			}
			tmpNode.lastKernel=k;
			if (tmpNode.lastPathWithStart>=0) {
				k=-1;
				for (j=0; j<=tmpNode.lastPathWithStart; j++) {
					kernel=tmpNode.pathsWithStart[j];
					if (kernel<0) kernel=-1-kernel;
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) tmpNode.pathsWithStart[++k]=tmpNode.pathsWithStart[j];
				}
				tmpNode.lastPathWithStart=k;
			}
			if (tmpNode.lastPathWithEnd>=0) {
				k=-1;
				for (j=0; j<=tmpNode.lastPathWithEnd; j++) {
					kernel=tmpNode.pathsWithEnd[j];
					if (kernel<0) kernel=-1-kernel;
					if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) tmpNode.pathsWithEnd[++k]=tmpNode.pathsWithEnd[j];
				}
				tmpNode.lastPathWithEnd=k;
			}
		}
	}
	
	
	/**
	 * Sorts the kernel graph topologically, in reverse order.
	 *
	 * Remark: the procedure sets global variables $sorted2original,original2sorted$, and
	 * it uses temporary arrays $tmpArray{1,2,3,4}$.
	 */
	private static final void sortKernelGraph() {
		int i;
		int nComponents, isDAG;
		int[] components, lastOutNeighborPrime, minimalVertices, nVerticesPerComponent;
		
		if (sorted2original==null || sorted2original.length<nPathKernels) sorted2original = new int[nPathKernels];
		if (original2sorted==null || original2sorted.length<nPathKernels) original2sorted = new int[nPathKernels];
		for (i=0; i<nPathKernels; i++) lastInNeighbor[i]++;
		for (i=0; i<nPathKernels; i++) lastOutNeighbor[i]++;
		components=tmpArray1;
		nComponents=DAG.getConnectedComponents(nPathKernels,inNeighbors,lastInNeighbor,outNeighbors,lastOutNeighbor,components,IntervalGraph.stack);
		if (nComponents<=0 || nComponents>nPathKernels) {
			IO.printCriticalErr("sortKernelGraph> ERROR: wrong number of connected components.");
			System.exit(1);
		}
		lastOutNeighborPrime=tmpArray2;
		System.arraycopy(lastOutNeighbor,0,lastOutNeighborPrime,0,nPathKernels);  // Since the following procedure modifies its $lastInNeighbor$ argument.
		minimalVertices=tmpArray3; nVerticesPerComponent=tmpArray4;
		isDAG=DAG.topologicalSort(nPathKernels,/*Inverting the direction of all arcs*/outNeighbors,lastOutNeighborPrime,/*Inverting the direction of all arcs*/inNeighbors,lastInNeighbor,components,nComponents,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
		if (isDAG!=0) {
			System.err.println("sortKernelGraph> ERROR: the kernel graph contains a cycle that involves node "+(isDAG-1));
			System.exit(1);
		}
		for (i=0; i<nPathKernels; i++) lastInNeighbor[i]--;
		for (i=0; i<nPathKernels; i++) lastOutNeighbor[i]--;
	}
	
	
	/**
	 * Keeps only the kernels whose label is set to true in $selectedLabels$, propagates
	 * such kernel tags to all their descendants in the kernel graph (if $mode=TRUE$), or
	 * to all their descendants that do not have themselves a selected label (if 
	 * $mode=FALSE$), and stores the result in $kernel2NewKernels$.
	 *
	 * @param computeAll TRUE: the procedure computes also 
	 * $kernel2NewKernels_{start,end}$;
	 * @param sorted2original topological sort of the kernel graph (in reverse order);
	 * @param kernel2NewKernels_orientation output matrix: cell $(i,j)$ contains the 
	 * orientation of kernel $kernel2NewKernels[i][j]$ WRT kernel $i$.
	 * @param kernel2NewKernels_{start,end} output matrix: the $i$-th row contains the IDs
	 * of the ends of all kernels in $kernel2NewKernels[i]$ that are aligned to the start
	 * (resp. end) of kernel $i$.
	 */
	private static final void markRepeats_propagateKernels(boolean[] selectedLabels, boolean mode, int[] sorted2original, byte[] labels, int[][] kernel2NewKernels, int[] lastNewKernel, byte[][] kernel2NewKernels_orientation, boolean computeAll, int[][] kernel2NewKernels_start, int[] kernel2NewKernels_lastStart, int[][] kernel2NewKernels_end, int[] kernel2NewKernels_lastEnd, int growthRate) {
		byte arcOrientation, insertionFwd, insertionRev;
		int i, j, k;
		int last, kernel, neighbor;
		
		Math.set(lastNewKernel,nPathKernels-1,-1);
		if (kernel2NewKernels_lastStart!=null) Math.set(kernel2NewKernels_lastStart,nPathKernels-1,-1);
		if (kernel2NewKernels_lastEnd!=null) Math.set(kernel2NewKernels_lastEnd,nPathKernels-1,-1);
		for (i=nPathKernels-1; i>=0; i--) {
			kernel=sorted2original[i];
			if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
			if (selectedLabels[labels[kernel]]) {
				if (mode) {
					// Updating kernels and orientations
					tmpArray1[0]=kernel;
					last=Math.setUnion(kernel2NewKernels[kernel],lastNewKernel[kernel],tmpArray1,0,tmpArray2);
					updateOrientations_singleton(kernel,(byte)0,kernel2NewKernels[kernel],lastNewKernel[kernel],tmpArray2,last,kernel2NewKernels_orientation[kernel],tmpByte1);
					if (last>=kernel2NewKernels_orientation[kernel].length) kernel2NewKernels_orientation[kernel] = new byte[last+growthRate];
					System.arraycopy(tmpByte1,0,kernel2NewKernels_orientation[kernel],0,last+1);
					if (last>=kernel2NewKernels[kernel].length) kernel2NewKernels[kernel] = new int[last+growthRate];
					System.arraycopy(tmpArray2,0,kernel2NewKernels[kernel],0,last+1);
					lastNewKernel[kernel]=last;
					if (computeAll) {
						// Updating start
						tmpArray1[0]=kernel;
						last=Math.setUnion(kernel2NewKernels_start[kernel],kernel2NewKernels_lastStart[kernel],tmpArray1,0,tmpArray2);
						if (last>=kernel2NewKernels_start[kernel].length) kernel2NewKernels_start[kernel] = new int[last+growthRate];
						System.arraycopy(tmpArray2,0,kernel2NewKernels_start[kernel],0,last+1);
						kernel2NewKernels_lastStart[kernel]=last;
						// Updating end
						tmpArray1[0]=-1-kernel;
						last=Math.setUnion(kernel2NewKernels_end[kernel],kernel2NewKernels_lastEnd[kernel],tmpArray1,0,tmpArray2);
						if (last>=kernel2NewKernels_end[kernel].length) kernel2NewKernels_end[kernel] = new int[last+growthRate];
						System.arraycopy(tmpArray2,0,kernel2NewKernels_end[kernel],0,last+1);
						kernel2NewKernels_lastEnd[kernel]=last;
					}
				}
				else {
					kernel2NewKernels[kernel][0]=kernel;
					lastNewKernel[kernel]=0;
					if (kernel2NewKernels_orientation[kernel]==null || kernel2NewKernels_orientation[kernel].length==0) kernel2NewKernels_orientation[kernel] = new byte[1];
					kernel2NewKernels_orientation[kernel][0]=0;
					if (computeAll) {
						if (kernel2NewKernels_start[kernel]==null || kernel2NewKernels_start[kernel].length==0) kernel2NewKernels_start[kernel] = new int[1];
						kernel2NewKernels_start[kernel][0]=kernel;
						kernel2NewKernels_lastStart[kernel]=0;
						if (kernel2NewKernels_end[kernel]==null || kernel2NewKernels_end[kernel].length==0) kernel2NewKernels_end[kernel] = new int[1];
						kernel2NewKernels_end[kernel][0]=-1-kernel;
						kernel2NewKernels_lastEnd[kernel]=0;
					}
				}
			}
			for (j=0; j<=lastOutNeighbor[kernel]; j++) {
				neighbor=outNeighbors[kernel][j];
				if (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0)) continue;
				if (selectedLabels[labels[neighbor]]) {
					if (!mode) continue;
				}
				k=Arrays.binarySearch(inNeighbors[neighbor],0,lastInNeighbor[neighbor]+1,kernel);
				arcOrientation=inNeighborFlags2orientation[inNeighborFlags[neighbor][(k<<2)+1]][inNeighborFlags[neighbor][(k<<2)+3]];
				insertionFwd=inNeighborFlags[neighbor][(k<<2)+0];
				insertionRev=inNeighborFlags[neighbor][(k<<2)+2];
				// Updating kernels and orientations
				if (lastNewKernel[neighbor]==-1) {
					last=lastNewKernel[kernel];
					if (last>=kernel2NewKernels[neighbor].length) kernel2NewKernels[neighbor] = new int[last+growthRate];
					System.arraycopy(kernel2NewKernels[kernel],0,kernel2NewKernels[neighbor],0,last+1);
					lastNewKernel[neighbor]=last;
					if (last>=kernel2NewKernels_orientation[neighbor].length) kernel2NewKernels_orientation[neighbor] = new byte[last+growthRate];
					System.arraycopy(kernel2NewKernels_orientation[kernel],0,kernel2NewKernels_orientation[neighbor],0,last+1);
					if (arcOrientation==1) {
						for (k=0; k<=last; k++) {
							if (kernel2NewKernels_orientation[neighbor][k]!=-1) kernel2NewKernels_orientation[neighbor][k]=oppositeOrientation[kernel2NewKernels_orientation[neighbor][k]];
						}
					}
					else if (arcOrientation==2) {
						for (k=0; k<=last; k++) kernel2NewKernels_orientation[neighbor][k]=2;
					}
				}
				else {
					last=Math.setUnion(kernel2NewKernels[kernel],lastNewKernel[kernel],kernel2NewKernels[neighbor],lastNewKernel[neighbor],tmpArray1);					
					updateOrientations(kernel2NewKernels[kernel],lastNewKernel[kernel],kernel2NewKernels_orientation[kernel],arcOrientation,kernel2NewKernels[neighbor],lastNewKernel[neighbor],kernel2NewKernels_orientation[neighbor],tmpArray1,last,tmpByte1);
					if (last>=kernel2NewKernels_orientation[neighbor].length) kernel2NewKernels_orientation[neighbor] = new byte[last+growthRate];
					System.arraycopy(tmpByte1,0,kernel2NewKernels_orientation[neighbor],0,last+1);
					if (last>=kernel2NewKernels[neighbor].length) kernel2NewKernels[neighbor] = new int[last+growthRate];
					System.arraycopy(tmpArray1,0,kernel2NewKernels[neighbor],0,last+1);
					lastNewKernel[neighbor]=last;
				}
				// Updating start/end
				if (computeAll) {
					updateStart(kernel,neighbor,arcOrientation,insertionFwd,insertionRev,kernel2NewKernels_start,kernel2NewKernels_lastStart,kernel2NewKernels_end,kernel2NewKernels_lastEnd,growthRate,tmpArray1);
					updateEnd(kernel,neighbor,arcOrientation,insertionFwd,insertionRev,kernel2NewKernels_start,kernel2NewKernels_lastStart,kernel2NewKernels_end,kernel2NewKernels_lastEnd,growthRate,tmpArray1);
				}
			}
		}
	}
	
	
	/**
	 * Stores in $newOrientations$ the list of orientations for all kernels in 
	 * $newKernels$, where $newKernels = oldKernels1 U oldKernels2$ and $oldKernels1$ is 
	 * taken in $orientation1$.
	 */
	private static final void updateOrientations(int[] oldKernels1, int lastOldKernel1, byte[] oldOrientations1, byte orientation1, int[] oldKernels2, int lastOldKernel2, byte[] oldOrientations2, int[] newKernels, int lastNewKernel, byte[] newOrientations) {
		byte o;
		int i, j;
		
		Math.set(newOrientations,lastNewKernel,(byte)(-2));
		i=0; j=0;
		while (j<=lastOldKernel2) {
			if (newKernels[i]<oldKernels2[j]) {
				i++;
				continue;
			}
			else if (newKernels[i]>oldKernels2[j]) {
				j++;
				continue;
			}
			newOrientations[i]=oldOrientations2[j];
			i++; j++;
		}
		i=0; j=0;
		while (j<=lastOldKernel1) {
			if (newKernels[i]<oldKernels1[j]) {
				i++;
				continue;
			}
			else if (newKernels[i]>oldKernels1[j]) {
				j++;
				continue;
			}
			if (oldOrientations1[j]==-1 || orientation1==-1) newOrientations[i]=-1;
			else if (newOrientations[i]!=-1) {
				o=transformOrientation[oldOrientations1[j]][orientation1];
				if (newOrientations[i]==-2) newOrientations[i]=o;
				else newOrientations[i]=addOrientation[newOrientations[i]][o];
			}
			i++; j++;
		}
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<=lastNewKernel; i++) {
				if (newOrientations[i]==-2) {
					System.err.println("updateOrientations> ERROR: an orientation equals -2");
					System.exit(1);
				}
			}
		}
	}
	
	
	/**
	 * Stores in $newOrientations$ the list of orientations for all kernels in 
	 * $newKernels$, where $newKernels = oldKernels \setunion {kernel}$, and $kernel$
	 * is in orientation $kernelOrientation$.
	 */
	private static final void updateOrientations_singleton(int kernel, byte kernelOrientation, int[] oldKernels, int lastOldKernel, int[] newKernels, int lastNewKernel, byte[] oldOrientations, byte[] newOrientations) {
		int i, j, k;
		
		Math.set(newOrientations,lastNewKernel,(byte)(-1));
		i=0; j=0; k=-1;
		while (j<=lastOldKernel) {
			if (newKernels[i]==kernel) k=i;
			if (newKernels[i]<oldKernels[j]) {
				i++;
				continue;
			}
			else if (newKernels[i]>oldKernels[j]) {
				j++;
				continue;
			}
			newOrientations[i]=oldOrientations[j];
			i++; j++;
		}
		if (k==-1) newOrientations[lastNewKernel]=kernelOrientation;
		else {
			if (newOrientations[k]==-1 || kernelOrientation==-1) newOrientations[k]=kernelOrientation;
			else newOrientations[k]=addOrientation[kernelOrientation][newOrientations[k]];
		}
	}
	
	
	/**
	 * Cell $(a_1,a_2)$ tells which orientation a block $(i_1,a_1,i_2,a_2)$ of 
	 * $inNeighborFlags$ encodes. Orientation values: -1=unknown, 0=forward, 1=RC, 2=both.
	 */
	private static final byte[][] inNeighborFlags2orientation = new byte[][] {
		{-1,1}, {0,2}
	};
	
	
	/**
	 * Cell $(i,j)$ tells how orientation $i$ appears when seen from orientation $j$.
	 * Orientation values: 0=forward, 1=RC, 2=both.
	 */
	private static final byte[][] transformOrientation = new byte[][] {
		{0,1,2}, {1,0,2}, {2,2,2}
	};


	/**
	 * Adds to $kernel2NewKernels_start[toKernel]$ the entries of 
	 * $kernel2NewKernels_start[fromKernel]$ or $kernel2NewKernels_end[fromKernel]$,
	 * depending on the orientation of the $(fromKernel,toKernel)$ arc of the kernel graph
	 * that connects them.
	 *
	 * @param tmpArray temporary space, of size at least $2*nPathKernels$.
	 */
	private static final void updateStart(int fromKernel, int toKernel, byte arcOrientation, byte insertionFwd, byte insertionRev, int[][] kernel2NewKernels_start, int[] kernel2NewKernels_lastStart, int[][] kernel2NewKernels_end, int[] kernel2NewKernels_lastEnd, int growthRate, int[] tmpArray) {
		int last;
		
		if ((arcOrientation==0 || arcOrientation==2) && (insertionFwd==1 || insertionFwd==3)) {
			last=kernel2NewKernels_lastStart[fromKernel];
			if (last>=0) {
				if (kernel2NewKernels_lastStart[toKernel]==-1) {
					if (last>=kernel2NewKernels_start[toKernel].length) kernel2NewKernels_start[toKernel] = new int[last+growthRate];
					System.arraycopy(kernel2NewKernels_start[fromKernel],0,kernel2NewKernels_start[toKernel],0,last+1);
					kernel2NewKernels_lastStart[toKernel]=last;
				}
				else {
					last=Math.setUnion(kernel2NewKernels_start[fromKernel],last,kernel2NewKernels_start[toKernel],kernel2NewKernels_lastStart[toKernel],tmpArray);
					if (last>=kernel2NewKernels_start[toKernel].length) kernel2NewKernels_start[toKernel] = new int[last+growthRate];
					System.arraycopy(tmpArray,0,kernel2NewKernels_start[toKernel],0,last+1);
					kernel2NewKernels_lastStart[toKernel]=last;
				}
			}
		}
		if ((arcOrientation==1 || arcOrientation==2) && (insertionRev==1 || insertionRev==3)) {
			last=kernel2NewKernels_lastEnd[fromKernel];
			if (last>=0) {
				if (kernel2NewKernels_lastStart[toKernel]==-1) {
					if (last>=kernel2NewKernels_start[toKernel].length) kernel2NewKernels_start[toKernel] = new int[last+growthRate];
					System.arraycopy(kernel2NewKernels_end[fromKernel],0,kernel2NewKernels_start[toKernel],0,last+1);
					kernel2NewKernels_lastStart[toKernel]=last;
				}
				else {
					last=Math.setUnion(kernel2NewKernels_end[fromKernel],last,kernel2NewKernels_start[toKernel],kernel2NewKernels_lastStart[toKernel],tmpArray);
					if (last>=kernel2NewKernels_start[toKernel].length) kernel2NewKernels_start[toKernel] = new int[last+growthRate];
					System.arraycopy(tmpArray,0,kernel2NewKernels_start[toKernel],0,last+1);
					kernel2NewKernels_lastStart[toKernel]=last;
				}
			}
		}
	}
	

	/**
	 * Symmetrical to $updateStart()$.
	 */
	private static final void updateEnd(int fromKernel, int toKernel, byte arcOrientation, byte insertionFwd, byte insertionRev, int[][] kernel2NewKernels_start, int[] kernel2NewKernels_lastStart, int[][] kernel2NewKernels_end, int[] kernel2NewKernels_lastEnd, int growthRate, int[] tmpArray) {
		int last;

		if ((arcOrientation==0 || arcOrientation==2) && (insertionFwd==1 || insertionFwd==2)) {
			last=kernel2NewKernels_lastEnd[fromKernel];
			if (last>=0) {
				if (kernel2NewKernels_lastEnd[toKernel]==-1) {
					if (last>=kernel2NewKernels_end[toKernel].length) kernel2NewKernels_end[toKernel] = new int[last+growthRate];
					System.arraycopy(kernel2NewKernels_end[fromKernel],0,kernel2NewKernels_end[toKernel],0,last+1);
					kernel2NewKernels_lastEnd[toKernel]=last;
				}
				else {
					last=Math.setUnion(kernel2NewKernels_end[fromKernel],last,kernel2NewKernels_end[toKernel],kernel2NewKernels_lastEnd[toKernel],tmpArray);
					if (last>=kernel2NewKernels_end[toKernel].length) kernel2NewKernels_end[toKernel] = new int[last+growthRate];
					System.arraycopy(tmpArray,0,kernel2NewKernels_end[toKernel],0,last+1);
					kernel2NewKernels_lastEnd[toKernel]=last;
				}
			}
		}
		if ((arcOrientation==1 || arcOrientation==2) && (insertionRev==1 || insertionRev==2)) {
			last=kernel2NewKernels_lastStart[fromKernel];
			if (last>=0) {
				if (kernel2NewKernels_lastEnd[toKernel]==-1) {
					if (last>=kernel2NewKernels_end[toKernel].length) kernel2NewKernels_end[toKernel] = new int[last+growthRate];
					System.arraycopy(kernel2NewKernels_start[fromKernel],0,kernel2NewKernels_end[toKernel],0,last+1);
					kernel2NewKernels_lastEnd[toKernel]=last;
				}
				else {
					last=Math.setUnion(kernel2NewKernels_start[fromKernel],last,kernel2NewKernels_end[toKernel],kernel2NewKernels_lastEnd[toKernel],tmpArray);
					if (last>=kernel2NewKernels_end[toKernel].length) kernel2NewKernels_end[toKernel] = new int[last+growthRate];
					System.arraycopy(tmpArray,0,kernel2NewKernels_end[toKernel],0,last+1);
					kernel2NewKernels_lastEnd[toKernel]=last;
				}
			}
		}
	}

	
	/**
	 * Resets $node.kernels$ to the union of the corresponding rows of 
	 * $kernel2NewKernels$.
	 *
	 * Remark: after this, $node$ might have no kernel tag, since rows of 
	 * $kernel2NewKernels$ might be empty.
	 * Remark: the procedure updates also arrays $kernelOrientations$ and $pathsWith*$.
	 * Remark: the procedure uses $tmpArray{1,2,3}$, which are assumed to be of size at 
	 * least $2*nPathKernels$ each.
	 *
	 * @param last4 if not -1: keeps only the new kernels that appear in $tmpArray4[0..
	 * last4]$, if any (this option is currently unused);
	 * @param mode TRUE=only non-periodic path kernels; FALSE=only periodic path kernels;
	 * in the latter case, every periodic or cyclic tag is deleted from the $pathsWith*$ 
	 * arrays.
	 */
	private static final void markRepeats_updateKernels(IntervalGraph.Node node, boolean mode, int[][] kernel2NewKernels, int[] lastNewKernel, byte[][] kernel2NewKernels_orientation, int[][] kernel2NewKernels_start, int[] kernel2NewKernels_lastStart, int[][] kernel2NewKernels_end, int[] kernel2NewKernels_lastEnd, int last4) {
		byte kernelOrientation;
		int i, j;
		int kernel, lastKernel, last1, last2, inContainer;
		byte[] tmpBytes;
		int[] tmpInts;
		
		// Mapping $node.kernels$ to $kernel2NewKernels$.
		lastKernel=node.lastKernel; last1=-1;
		for (i=0; i<=lastKernel; i++) {
			kernel=node.kernels[i];
			kernelOrientation=node.kernelOrientations[i];
			if ( (mode && (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0))) || 
			     (!mode && pathKernelLengths[kernel]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[kernel]==0))
			   ) {
				// These tags are not altered
				if (last1==-1) {
					tmpArray1[0]=kernel;
					tmpByte1[0]=kernelOrientation;
					last1=0;
				}
				else {
					tmpArray3[0]=kernel;
					last2=Math.setUnion(tmpArray3,0,tmpArray1,last1,tmpArray2);
					updateOrientations_singleton(kernel,kernelOrientation,tmpArray1,last1,tmpArray2,last2,tmpByte1,tmpByte2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
					tmpBytes=tmpByte1; tmpByte1=tmpByte2; tmpByte2=tmpBytes;
					last1=last2;
				}
			}
			else if (lastNewKernel[kernel]>=0) {
				last2=Math.setUnion(kernel2NewKernels[kernel],lastNewKernel[kernel],tmpArray1,last1,tmpArray2);
				updateOrientations(kernel2NewKernels[kernel],lastNewKernel[kernel],kernel2NewKernels_orientation[kernel],kernelOrientation,tmpArray1,last1,tmpByte1,tmpArray2,last2,tmpByte2);
				tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				tmpBytes=tmpByte1; tmpByte1=tmpByte2; tmpByte2=tmpBytes;
				last1=last2;
			}
			else {
				// The kernel tag is removed if $kernel2NewKernels$ ha no elements for it.
			}
		}
		
		// Keeping only the new kernels that appear in $tmpArray4$, if any.
		if (last4!=-1) {
			last2=-1; j=0; inContainer=0;
			for (i=0; i<=last1; i++) {
				kernel=tmpArray1[i];
				kernelOrientation=tmpByte1[i];
				if ( (mode && (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0))) || 
				     (!mode && pathKernelLengths[kernel]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[kernel]==0))
				   ) {
					// These tags are not altered
					last2++;
					tmpArray2[last2]=kernel;
					tmpByte2[last2]=kernelOrientation;
				}
				else {
					while (j<=last4 && tmpArray4[j]<kernel) j++;
					if (j<=last4 && tmpArray4[j]==kernel) {
						last2++;
						tmpArray2[last2]=kernel;
						tmpByte2[last2]=kernelOrientation;
						inContainer++;
					}
				}
			}
			if (IO.CONSISTENCY_CHECKS) {
				if (last2>last1 || Math.setIntersectionSize(tmpArray2,0,last2,tmpArray1,0,last1,Math.POSITIVE_INFINITY)!=last2+1) {
					System.err.println("markRepeats_updateKernels> ERROR: not subset");
					System.exit(1);
				}
			}
			if (inContainer>0) {
				tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				tmpBytes=tmpByte1; tmpByte1=tmpByte2; tmpByte2=tmpBytes;
				last1=last2;
			}
		}
		
		// Copying $tmpArray1,tmpByte1$ to $node$.
		if (last1>=0) {
			if (last1>=node.kernels.length) node.kernels = new int[last1+1];
			System.arraycopy(tmpArray1,0,node.kernels,0,last1+1);
			if (last1>=node.kernelOrientations.length) node.kernelOrientations = new byte[last1+1];
			System.arraycopy(tmpByte1,0,node.kernelOrientations,0,last1+1);
		}
		node.lastKernel=last1;
		
		// Updating start and end
		if (!mode) {
			j=-1;
			for (i=0; i<=node.lastPathWithStart; i++) {
				kernel=node.pathsWithStart[i];
				if (kernel<0) kernel=-1-kernel;
				if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
				node.pathsWithStart[++j]=node.pathsWithStart[i];
			}
			node.lastPathWithStart=j;
			j=-1;
			for (i=0; i<=node.lastPathWithEnd; i++) {
				kernel=node.pathsWithEnd[i];
				if (kernel<0) kernel=-1-kernel;
				if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) continue;
				node.pathsWithEnd[++j]=node.pathsWithEnd[i];
			}
			node.lastPathWithEnd=j;
			return;
		}
		if (node.lastPathWithStart!=-1) {
			last1=-1;
			for (i=0; i<=node.lastPathWithStart; i++) {
				kernel=node.pathsWithStart[i];
				if (kernel>=0 && kernel2NewKernels_lastStart[kernel]!=-1) {
					last1=Math.setUnion(tmpArray1,last1,kernel2NewKernels_start[kernel],kernel2NewKernels_lastStart[kernel],tmpArray2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				}
				else if (kernel<0 && kernel2NewKernels_lastEnd[-1-kernel]!=-1) {
					last1=Math.setUnion(tmpArray1,last1,kernel2NewKernels_end[-1-kernel],kernel2NewKernels_lastEnd[-1-kernel],tmpArray2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				}
				else { /* The kernel tag is removed */ }
			}
			j=-1;
			for (i=0; i<=last1; i++) {
				kernel=tmpArray1[i];
				if (kernel<0) kernel=-1-tmpArray1[i];
				if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,kernel)>=0) tmpArray1[++j]=tmpArray1[i];
			}
			last1=j;
			if (last1!=-1) {
				if (node.pathsWithStart.length<last1+1) node.pathsWithStart = new int[last1+1];
				System.arraycopy(tmpArray1,0,node.pathsWithStart,0,last1+1);
			}
			node.lastPathWithStart=last1;
		}
		if (node.lastPathWithEnd!=-1) {
			last1=-1;
			for (i=0; i<=node.lastPathWithEnd; i++) {
				kernel=node.pathsWithEnd[i];
				if (kernel>=0 && kernel2NewKernels_lastStart[kernel]!=-1) {
					last1=Math.setUnion(tmpArray1,last1,kernel2NewKernels_start[kernel],kernel2NewKernels_lastStart[kernel],tmpArray2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				}
				else if (kernel<0 && kernel2NewKernels_lastEnd[-1-kernel]!=-1) {
					last1=Math.setUnion(tmpArray1,last1,kernel2NewKernels_end[-1-kernel],kernel2NewKernels_lastEnd[-1-kernel],tmpArray2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				}
				else { /* The kernel tag is removed */ }
			}
			j=-1;
			for (i=0; i<=last1; i++) {
				kernel=tmpArray1[i];
				if (kernel<0) kernel=-1-tmpArray1[i];
				if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,kernel)>=0) tmpArray1[++j]=tmpArray1[i];
			}
			last1=j;
			if (last1!=-1) {
				if (node.pathsWithEnd.length<last1+1) node.pathsWithEnd = new int[last1+1];
				System.arraycopy(tmpArray1,0,node.pathsWithEnd,0,last1+1);
			}
			node.lastPathWithEnd=last1;
		}
	}
	
	
	/**
	 * The procedure stores in $isRedundant[i]$ the number of interval graph nodes that 
	 * are marked with the $i$-th kernel in $kernels25$ and that are not contained
	 * in the interval of any kernel in $ancestors[i]$ (we call "exposed" such intervals).
	 * This number is not incremented after it reaches $threshold$. 
	 * This can be seen as a test of near-supermaximality, limited to ancestors in the 
	 * kernel graph.
	 *
	 * Remark: the procedure does not check if the intervals that were used to build a 
	 * kernel K are exposed or not. Assume that K was built from multiple intervals. Then,
	 * the interval I that contains the prefix of K has low quality on the end in which it 
	 * is interrupted, but K might be the proper prefix of another kernel K', so I might 
	 * not be contained in any interval of K'. I cannot belong to the basin of K', either,
	 * since it is by construction maximal by containment. If I is an internal substring 
	 * of K, then it is interrupted on both sides, so no interval of K' is likely to 
	 * contain it. If K consists of a single interval I, and if I is adjacent to at least
	 * one low-quality region, the same issues arise. Otherwise, an interval of K' 
	 * containing I might indeed exist.
	 *
	 * Remark: K might consist of a single interval, and many instances of K that are 
	 * connected by identity edges might belong to the bidirected graph. We use such 
	 * instances as candidate exposed, as well, iff both their ends are adjacent to high-
	 * quality parts of a read; i.e. we assume that no ancestor of K with similar length 
	 * as K (if any) occurs at those instances.
	 *
	 * Remark: for simplicity, the procedure compares just nodes of the interval graph.
	 * It could be updated to use the more accurate approach of kernel tracks as done by
	 * $markExposed()$ in Step5.
	 *
	 * Remark: the procedure assumes that $updateNKernelNodes()$ has already been run.
	 *
	 * Remark: the procedure uses $tmpArray1$, which is assumed to be of size at least 
	 * $nPathKernels$.
	 *
	 * @param from,to all and only the interval graph nodes in the same read are assumed 
	 * to belong to this interval of $IntervalGraph.nodesArray$, sorted by $start$;
	 * @param oneBidirectedNode output of $markKernelsWithOneBidirectedNode()$;
	 * @return the number of elements in $kernels25$ whose $isRedundant$ count has been
	 * pushed to $threshold$ by the procedure.
	 */
	private static final int markRepeats_markRedundantKernels(int from, int to, int[] kernels25, int lastKernel25, int[][] ancestors, int[] lastAncestor, int[] isRedundant, int threshold, boolean[] oneBidirectedNode) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int QUALITY_WINDOW = IO.quantum<<1;  // Arbitrary
		boolean found;
		int i, j, k, p;
		int last1, kernel, read, start, end, out, bidirectedGraphNode;
		IntervalGraph.Node node, otherNode;
		
		out=0;
		for (i=from; i<=to; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastKernel==-1) continue;
			last1=Math.setIntersection(node.kernels,0,node.lastKernel,kernels25,0,lastKernel25,tmpArray1,0);
			if (last1==-1) continue;
			read=node.read; start=node.start; end=node.end;
			if (!node.isLeftMaximal && !node.isRightMaximal) continue;
			p=0;
			for (j=0; j<=last1; j++) {
				kernel=tmpArray1[j];
				while (p<=lastKernel25 && kernels25[p]<kernel) p++;
				if (isRedundant[p]==threshold) continue;
				if (!oneBidirectedNode[kernel]) {
					found=false;
					for (k=nKernelNodes[kernel]; k<blockKernels.length; k++) {
						if (blockKernels[k]!=kernel) break;
						if (IntervalGraph.nodesArray[nodesInKernel[k]]==node) {
							found=true;
							break;
						}
					}
					if (found) continue; 
				}
				else {
					if (!node.isLeftMaximal || !node.isRightMaximal) continue;
				}
				if (Math.nonemptyIntersection(node.kernels,0,node.lastKernel,ancestors[p],0,lastAncestor[p])) continue;
				found=false;
				for (k=i-1; k>=from; k--) {
					otherNode=IntervalGraph.nodesArray[k];
					if (otherNode.end<=node.start) continue;
					if (!Intervals.isApproximatelyContained(node.start,node.end,otherNode.start,otherNode.end)) continue;
					if (Math.nonemptyIntersection(otherNode.kernels,0,otherNode.lastKernel,ancestors[p],0,lastAncestor[p])) {
						found=true;
						break;
					}
				}
				if (!found) {
					for (k=i+1; k<=to; k++) {
						otherNode=IntervalGraph.nodesArray[k];
						if (otherNode.start>=node.end) break;
						if (!Intervals.isApproximatelyContained(node.start,node.end,otherNode.start,otherNode.end)) continue;
						if (Math.nonemptyIntersection(otherNode.kernels,0,otherNode.lastKernel,ancestors[p],0,lastAncestor[p])) {
							found=true;
							break;
						}
					}
				}
				if (!found) {
					isRedundant[p]++;
					if (isRedundant[p]==threshold) out++;
				}
			}
		}
		return out;
	}
	
	
	/**
	 * Sets $out[i]=true$ iff kernel $i$ has just one bidirected graph node.
	 * $out$ is assumed to have at least $nPathKernels$ cells.
	 */
	private static final void markKernelsWithOneBidirectedNode(boolean[] out) {
		int i, j;
		int kernel, bidirectedGraphNode;
		
		Math.set(out,nPathKernels-1,false);
		kernel=blockKernels[0]; 
		bidirectedGraphNode=IntervalGraph.nodesArray[nodesInKernel[0]].bidirectedGraphNode;
		for (i=1; i<nNodesInKernel; i++) {
			if (blockKernels[i]!=kernel) {
				if (bidirectedGraphNode>=0) out[kernel]=true;
				kernel=blockKernels[i];
				bidirectedGraphNode=IntervalGraph.nodesArray[nodesInKernel[i]].bidirectedGraphNode;
				continue;
			}
			j=IntervalGraph.nodesArray[nodesInKernel[i]].bidirectedGraphNode;
			if (j==-1) continue;
			if (bidirectedGraphNode==-1) bidirectedGraphNode=j;
			else if (bidirectedGraphNode>=0 && j!=bidirectedGraphNode) bidirectedGraphNode=-2;
		}
		if (bidirectedGraphNode>=0) out[kernel]=true;
	}
	
	
	/**
	 * Removes redundant kernels from $ancestors*$ arrays.
	 */
	private static final void markRepeats_removeRedundantAncestors(int[] kernels25, int lastKernel25, byte[] labels, int[][] ancestors, int[] lastAncestor, byte[][] ancestorOrientations, int[][] ancestors_start, int[] ancestors_lastStart, int[][] ancestors_end, int[] ancestors_lastEnd) {
		int i, j, k, p;
		int last1, kernel;
		
		for (i=0; i<=lastKernel25; i++) {
			last1=lastAncestor[i]; k=-1;
			for (j=0; j<=last1; j++) {
				kernel=ancestors[i][j];
				if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0) || labels[kernel]!=0) {
					k++;
					ancestors[i][k]=kernel;
					ancestorOrientations[i][k]=ancestorOrientations[i][j];
				}
			}
			lastAncestor[i]=k;
		}
		for (i=0; i<=lastKernel25; i++) {
			last1=ancestors_lastStart[i]; k=-1;
			for (j=0; j<=last1; j++) {
				kernel=ancestors_start[i][j];
				p=kernel>=0?kernel:-1-kernel;
				if (pathKernelLengths[p]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[p]!=0) || labels[p]!=0) ancestors_start[i][++k]=kernel;
			}
			ancestors_lastStart[i]=k;
		}
		for (i=0; i<=lastKernel25; i++) {
			last1=ancestors_lastEnd[i]; k=-1;
			for (j=0; j<=last1; j++) {
				kernel=ancestors_end[i][j];
				p=kernel>=0?kernel:-1-kernel;
				if (pathKernelLengths[p]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[p]!=0) || labels[p]!=0) ancestors_end[i][++k]=kernel;
			}
			ancestors_lastEnd[i]=k;
		}
	}
	
	
	/**
	 * Replaces every redundant kernel in $node.kernels$ with the (non-empty) set of all 
	 * its non-redundant ancestors in $ancestors$. Redundant kernels are assumed to be a 
	 * subset of $kernels25$ and to be marked with a zero in $labels$.
	 *
	 * Remark: the procedure uses temporary arrays $tmpArray{1,2,3}$ and $tmpByte{1,2,3}$,
	 * which are assumed to be of length at least $nPathKernels$ each.
	 *
	 * @param last4 if not -1: keeps only the new kernels that appear in $tmpArray4[0..
	 * last4]$, if any (this option is currently unused);
	 * @param ancestors one row per element of $kernels25$.
	 */
	private static final void markRepeats_removeRedundantKernels(IntervalGraph.Node node, int[] kernels25, int lastKernel25, byte[] labels, int[][] ancestors, int[] lastAncestor, byte[][] ancestorOrientations, int[][] ancestors_start, int[] ancestors_lastStart, int[][] ancestors_end, int[] ancestors_lastEnd, int last4) {
		byte kernelOrientation;
		int i, j, k, p;
		int kernel, lastKernel, last1, last2, last3, inContainer;
		byte[] tmpBytes;
		int[] tmpInts;
		
		// Updating kernels and orientations
		lastKernel=node.lastKernel; last1=-1; p=0;
		for (i=0; i<=lastKernel; i++) {
			kernel=node.kernels[i];
			kernelOrientation=node.kernelOrientations[i];
			if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0) || labels[kernel]!=0) {
				// These tags are not altered
				if (last1==-1) {
					tmpArray1[0]=kernel;
					tmpByte1[0]=kernelOrientation;
					last1=0;
				}
				else {
					tmpArray2[0]=kernel;
					last3=Math.setUnion(tmpArray2,0,tmpArray1,last1,tmpArray3);
					updateOrientations_singleton(kernel,kernelOrientation,tmpArray1,last1,tmpArray3,last3,tmpByte1,tmpByte2);
					tmpInts=tmpArray1; tmpArray1=tmpArray3; tmpArray3=tmpInts;
					tmpBytes=tmpByte1; tmpByte1=tmpByte2; tmpByte2=tmpBytes;
					last1=last3;
				}
			}
			else {
				while (p<=lastKernel25 && kernels25[p]<kernel) p++;
				if (p>lastKernel25 || kernels25[p]!=kernel) {
					System.err.println("removeRedundantKernels> ERROR: kernel "+kernel+" has zero label but it is not in kernels25 ?!");
					System.err.print("removeRedundantKernels> kernels25: ");
					for (int x=0; x<=lastKernel25; x++) System.err.print(kernels25[x]+",");
					System.err.println();
					System.exit(1);
				}
				last2=lastAncestor[p];
				if (last2==-1) {
					System.err.println("removeRedundantKernels> ERROR: redundant kernel "+kernel+" has no ancestor?");
					System.exit(1);
				}				
				System.arraycopy(ancestors[p],0,tmpArray2,0,last2+1);
				System.arraycopy(ancestorOrientations[p],0,tmpByte2,0,last2+1);
				k=-1;
				for (j=0; j<=last2; j++) {
					if (labels[tmpArray2[j]]!=0) {
						k++;
						tmpArray2[k]=tmpArray2[j];
						tmpByte2[k]=tmpByte2[j];
					}
				}
				last2=k;
				if (last2>=0) {
					last3=Math.setUnion(tmpArray1,last1,tmpArray2,last2,tmpArray3);
					updateOrientations(tmpArray2,last2,tmpByte2,kernelOrientation,tmpArray1,last1,tmpByte1,tmpArray3,last3,tmpByte3);
					tmpInts=tmpArray1; tmpArray1=tmpArray3; tmpArray3=tmpInts;
					tmpBytes=tmpByte1; tmpByte1=tmpByte3; tmpByte3=tmpBytes;
					last1=last3;
				}
				else {
					// The kernel tag is removed if none of its ancestors was kept.
				}
			}
		}
		
		// Keeping only the new kernels that appear in $tmpArray4$, if any.
		if (last4!=-1) {
			last2=-1; j=0; inContainer=0;
			for (i=0; i<=last1; i++) {
				kernel=tmpArray1[i];
				kernelOrientation=tmpByte1[i];
				if (pathKernelLengths[kernel]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0) || labels[kernel]!=0) {
					// These tags are not altered
					last2++;
					tmpArray2[last2]=kernel;
					tmpByte2[last2]=kernelOrientation;
				}
				else {
					while (j<=last4 && tmpArray4[j]<kernel) j++;
					if (j<=last4 && tmpArray4[j]==kernel) {
						last2++;
						tmpArray2[last2]=kernel;
						tmpByte2[last2]=kernelOrientation;
						inContainer++;
					}
				}
			}
			if (IO.CONSISTENCY_CHECKS) {
				if (last2>last1 || Math.setIntersectionSize(tmpArray2,0,last2,tmpArray1,0,last1,Math.POSITIVE_INFINITY)!=last2+1) {
					System.err.println("markRepeats_removeRedundantKernels> ERROR: not subset");
					System.exit(1);
				}
			}
			if (inContainer>0) {
				tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				tmpBytes=tmpByte1; tmpByte1=tmpByte2; tmpByte2=tmpBytes;
				last1=last2;
			}
		}
		
		// Copying $tmpArray1,tmpByte1$ to $node$.
		if (last1>=0) {
			if (last1>=node.kernels.length) node.kernels = new int[last1+1];
			if (last1>=node.kernelOrientations.length) node.kernelOrientations = new byte[last1+1];
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
				if (pathKernelLengths[k]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[k]!=0) || labels[k]!=0) {
					// These tags are not altered
					tmpArray2[0]=kernel;
					last3=Math.setUnion(tmpArray1,last1,tmpArray2,0,tmpArray3);
					tmpInts=tmpArray1; tmpArray1=tmpArray3; tmpArray3=tmpInts;
					last1=last3;
					continue;
				}
				p=Arrays.binarySearch(kernels25,0,lastKernel25+1,k);
				if (kernel>=0 && ancestors_lastStart[p]!=-1) {
					last1=Math.setUnion(tmpArray1,last1,ancestors_start[p],ancestors_lastStart[p],tmpArray2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				}
				else if (kernel<0 && ancestors_lastEnd[p]!=-1) {
					last1=Math.setUnion(tmpArray1,last1,ancestors_end[p],ancestors_lastEnd[p],tmpArray2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				}
				else { /* The kernel tag is removed */ }
			}
			j=-1;
			for (i=0; i<=last1; i++) {
				kernel=tmpArray1[i];
				if (kernel<0) kernel=-1-tmpArray1[i];
				if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,kernel)>=0) tmpArray1[++j]=tmpArray1[i];
			}
			last1=j;
			if (last1!=-1) {
				if (node.pathsWithStart.length<last1+1) node.pathsWithStart = new int[last1+1];
				System.arraycopy(tmpArray1,0,node.pathsWithStart,0,last1+1);
			}
			node.lastPathWithStart=last1;
		}
		if (node.lastPathWithEnd!=-1) {
			last1=-1;
			for (i=0; i<=node.lastPathWithEnd; i++) {
				kernel=node.pathsWithEnd[i];		
				k=kernel>=0?kernel:-1-kernel;
				if (pathKernelLengths[k]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[k]!=0) || labels[k]!=0) {
					// These tags are not altered
					tmpArray2[0]=kernel;
					last3=Math.setUnion(tmpArray1,last1,tmpArray2,0,tmpArray3);
					tmpInts=tmpArray1; tmpArray1=tmpArray3; tmpArray3=tmpInts;
					last1=last3;
					continue;
				}
				p=Arrays.binarySearch(kernels25,0,lastKernel25+1,k);
				if (kernel>=0 && ancestors_lastStart[p]!=-1) {
					last1=Math.setUnion(tmpArray1,last1,ancestors_start[p],ancestors_lastStart[p],tmpArray2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				}
				else if (kernel<0 && ancestors_lastEnd[p]!=-1) {
					last1=Math.setUnion(tmpArray1,last1,ancestors_end[p],ancestors_lastEnd[p],tmpArray2);
					tmpInts=tmpArray1; tmpArray1=tmpArray2; tmpArray2=tmpInts;
				}
				else { /* The kernel tag is removed */ }
			}
			j=-1;
			for (i=0; i<=last1; i++) {
				kernel=tmpArray1[i];
				if (kernel<0) kernel=-1-tmpArray1[i];
				if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,kernel)>=0) tmpArray1[++j]=tmpArray1[i];
			}
			last1=j;
			if (last1!=-1) {
				if (node.pathsWithEnd.length<last1+1) node.pathsWithEnd = new int[last1+1];
				System.arraycopy(tmpArray1,0,node.pathsWithEnd,0,last1+1);
			}
			node.lastPathWithEnd=last1;
		}
	}
	
	
	private static final void markRepeats_checkKernelsWithZeroLabel(byte[] labels) {
		int i, j;
		int last;
		IntervalGraph.Node tmpNode;
		
		for (i=0; i<IntervalGraph.nNodes; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			last=tmpNode.lastKernel;
			for (j=0; j<=last; j++) {
				if (labels[tmpNode.kernels[j]]==0) {
					System.err.println("checkKernels> ERROR: this interval graph node has a kernel marked as zero: "+tmpNode);
					System.exit(1);
				}
			}
		}
	}
	
	
	
	
	
	
	
	
	// -------------------- CYCLIC KERNELS AND PERIODIC PATH KERNELS ---------------------
	
	/**
	 * @return the fraction of periodic nodes with short period in the interval graph.
	 */
	private static final double getShortPeriodFraction() {
		final int nNodes = IntervalGraph.nNodes;
		int i;
		double out;
		IntervalGraph.Node node;
		
		out=0; 
		for (i=0; i<nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.type==Constants.INTERVAL_PERIODIC && !node.hasLongPeriod) out++;
		}
		return out/=nNodes;
	}
	
	
	/**
 	 * Simplified handling of an interval graph that is almost fully short-period. 
	 * Kernels and the kernel graph are not built (building the latter can be especially 
	 * slow, since there are most likely many alignments between every pair of kernels).
	 * Assembling with suffix-prefix overlaps is of questionable value in this case.
	 * Non-periodic nodes with a long substring for which there is no evidence of 
	 * periodicity are discarded.
	 */
	private static final void handleShortPeriodGraph(String outputDir, String outputPrefix, String inputGraphID, String graphFile) throws IOException {
		final int PERIODICS_SUBSTRINGS_GAP = IO.quantum<<2;  // Arbitrary
		int i;
		int maxNode, length, maxLength, maxPeriodicSubstrings, nNodesPrime;
		final int nNodes = IntervalGraph.nNodes;
		IntervalGraph.Node node;
		byte[] kernelLabels;
		int[][] newNodes;
		
		// Computing the surface covered by periodic substrings of each node
		IntervalGraph.sortNodesArray();
		IntervalGraph.setPeriodicSubstringsOfNodes_sameRead(false);
		IntervalGraph.setPeriodicSubstringsOfNodes_neighbors(true,tmpArray1);
		maxPeriodicSubstrings=IntervalGraph.mergePeriodicSubstringsOfNodes();
		maxPeriodicSubstrings<<=1;
		if (tmpArray1==null || tmpArray1.length<maxPeriodicSubstrings) tmpArray1 = new int[maxPeriodicSubstrings];
		
		// Setting as the only kernel a single short-period node of maximum length
		maxLength=0; maxNode=-1;
		for (i=0; i<nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.type==Constants.INTERVAL_PERIODIC && !node.hasLongPeriod) {
				length=node.length();
				if (length>maxLength) {
					maxLength=length;
					maxNode=i;
				}
			}
		}
		
		// Storing an artificial bidirected graph to disk
		node=IntervalGraph.nodesArray[maxNode];
		BidirectedGraph.allocateMemory(1,1,1);
		BidirectedGraph.clear(); BidirectedGraph.nNodes=1; BidirectedGraph.lastNeighbor[0]=-1; 
		BidirectedGraph.nodeLength[0]=maxLength; 
		BidirectedGraph.intervalGraphPointers[0]=node;
		BidirectedGraph.printAllPaths_singleton(0);
		BidirectedGraph.serialize(graphFile+"-kernel0.bdgraph");
		newNodes = new int[1][3];
		newNodes[0][0]=node.read; newNodes[0][1]=node.start; newNodes[0][2]=node.end; 
		printBidirectedGraphLabels(newNodes,1,graphFile+"-kernel0-"+IO.ONE_NODE_ONLY_LABEL+"-"+IO.LABEL_COORDINATES_LABEL+".txt");
		BidirectedGraph.deallocateMemory();
		
		// Assigning and storing kernel tags
		node.bidirectedGraphNode=0;
		node.lastKernel=0;
		if (node.kernels==null || node.kernels.length==0) node.kernels = new int[1];
		node.kernels[0]=0;
		if (node.kernelOrientations==null || node.kernelOrientations.length==0) node.kernelOrientations = new byte[1];
		node.kernelOrientations[0]=(byte)0;
		node.lastPathWithStart=-1; node.lastPathWithEnd=-1; nNodesPrime=1;
		for (i=0; i<maxNode; i++) {
			node=IntervalGraph.nodesArray[i];
			node.bidirectedGraphNode=-1;
			if (node.type!=Constants.INTERVAL_PERIODIC && node.gapInPeriodicSubstrings(PERIODICS_SUBSTRINGS_GAP,tmpArray1)) node.lastKernel=-1;
			else {
				nNodesPrime++;
				node.lastKernel=0;
				if (node.kernels==null || node.kernels.length==0) node.kernels = new int[1];
				node.kernels[0]=0;
				if (node.kernelOrientations==null || node.kernelOrientations.length==0) node.kernelOrientations = new byte[1];
				node.kernelOrientations[0]=(byte)0;
			}
			node.lastPathWithStart=-1; node.lastPathWithEnd=-1;
		}
		for (i=maxNode+1; i<nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			node.bidirectedGraphNode=-1;
			if (node.type!=Constants.INTERVAL_PERIODIC && node.gapInPeriodicSubstrings(PERIODICS_SUBSTRINGS_GAP,tmpArray1)) node.lastKernel=-1;
			else {
				nNodesPrime++;
				node.lastKernel=0;
				if (node.kernels==null || node.kernels.length==0) node.kernels = new int[1];
				node.kernels[0]=0;
				if (node.kernelOrientations==null || node.kernelOrientations.length==0) node.kernelOrientations = new byte[1];
				node.kernelOrientations[0]=(byte)0;
			}
			node.lastPathWithStart=-1; node.lastPathWithEnd=-1;
		}
		printKernelTags(outputDir+"/"+IO.TAGS_PREFIX+outputPrefix+"-"+inputGraphID+".txt",outputPrefix+"."+inputGraphID,true,false,-1);
		nPathKernels=1;
		kernelLabels = new byte[] {-3};
		pathKernel2Kernel = new int[] {0};
		pathKernelLengths = new int[] {maxLength};
		pathKernelPeriodic = new byte[] {1};
		kernelSize = new int[] {nNodesPrime};
		kernelFrequency = new int[1][4];
		getFrequencyPeriodic();
		printBasinDescriptors(kernelLabels,outputDir,IO.BASIN_DESCRIPTOR_PREFIX+"-"+outputPrefix+"-"+inputGraphID,null,true);
	}
	
	
	/**
	 * A path of a non-cyclic kernel can contain periodic interval graph nodes (e.g. the
	 * path might consist of a single interval graph node of periodic type). Conversely,
	 * a cyclic kernel might contain no periodic interval graph node (e.g. if the prefix 
	 * and the suffix of a repeat are identical, and if they are not factorized in
	 * distinct intervals). The procedure marks in $pathKernelPeriodic$ all the path 
	 * kernels (possibly cyclic) that contain an interval graph node of periodic type
	 * (1=short, 2=long).
	 *
	 * Remark: the procedure assumes that $expandNodesInKernel()$ has already been called.
	 *
	 * @return number of (possibly cyclic) path kernels marked by the procedure.
	 */
	public static final int buildPathKernelPeriodic() {
		boolean foundShort, foundLong;
		int i;
		int previousKernel, out;
		IntervalGraph.Node node;
		
		if (pathKernelPeriodic==null) pathKernelPeriodic = new byte[nPathKernels];
		Math.set(pathKernelPeriodic,nPathKernels-1,(byte)0);
		previousKernel=blockKernels[0];
		node=IntervalGraph.nodesArray[nodesInKernel[0]];
		if (node.type==Constants.INTERVAL_PERIODIC) {
			if (node.hasLongPeriod) { foundLong=true; foundShort=false; }
			else { foundShort=true; foundLong=false; }
		}
		else { foundShort=false; foundLong=false; }
		out=0;
		for (i=0; i<nNodesInKernel; i++) {
			node=IntervalGraph.nodesArray[nodesInKernel[i]];
			if (blockKernels[i]!=previousKernel) {
				if (foundShort) pathKernelPeriodic[previousKernel]=1;
				else if (foundLong) pathKernelPeriodic[previousKernel]=2;
				if (foundShort||foundLong) out++;
				previousKernel=blockKernels[i];
				if (node.type==Constants.INTERVAL_PERIODIC) {
					if (node.hasLongPeriod) foundLong=true;
					else foundShort=true;
				}
				else { foundShort=false; foundLong=false; }
				continue;
			}
			if (node.type==Constants.INTERVAL_PERIODIC) {
				if (node.hasLongPeriod) foundLong=true;
				else foundShort=true;
			}
		}
		if (foundShort) pathKernelPeriodic[previousKernel]=1;
		else if (foundLong) pathKernelPeriodic[previousKernel]=2;
		if (foundShort||foundLong) out++;
		return out;
	}
	
	
	/**
	 * Marks all nodes of the kernel graph that can be reached from those already marked
	 * in $pathKernelPeriodic$ by using all outgoing arcs, and all incoming identity arcs.
	 * A kernel that is reachable from long-period kernels but that is shorter than their 
	 * periods is not marked.
	 *
	 * Remark: the kernel graph built by $buildKernelGraph()$ includes periodic path 
	 * kernels and cyclic kernels. This is useful, since periodic interval graph 
	 * nodes with the same period might give rise to different kernels. This happens 
	 * because an edge in the interval graph is marked as containment using the same 
	 * criteria as for simple intervals, i.e. depending on the quality of the ends of the 
	 * intervals; but two periodic kernels joined by insertion edges should be collapsed
	 * when building a representative. The same holds for overlaps, which might be marked
	 * as shared substring edges in the interval graph. Moreover, containment/overlap is
	 * allowed only between intervals that are both long-period or both short-period, but
	 * this tag might be inaccurate. Finally, multiple occurrences of the same periodic 
	 * string might occur in the same read, and become distinct kernels for lack of 
	 * alignments.
	 *
	 * @param stack of size at least $nPathKernels$;
	 * @param tmpArray{1,2} temporary space, of size at least $nPathKernels$;
	 * @return number of new path kernels marked as periodic by the procedure.
	 */
	private static final int expandPathKernelPeriodic(int[] stack, int[] tmpArray1, int[] tmpArray2) {
		final int LONG_PERIOD_THRESHOLD = IO.quantum<<1;  // Arbitrary
		boolean valid;
		byte type;
		int i, j, p;
		int top, kernel, neighbor, period, out;
		IntervalGraph.Node node;
				
		// Estimating the length of long periods
		Math.set(tmpArray1,nPathKernels-1,0);
		Math.set(tmpArray2,nPathKernels-1,0);
		for (i=0; i<nNodesInKernel; i++) {
			node=IntervalGraph.nodesArray[nodesInKernel[i]];
			if (node.type==Constants.INTERVAL_PERIODIC && node.hasLongPeriod && node.period>0) {
				kernel=blockKernels[i];
				tmpArray1[kernel]+=node.period;
				tmpArray2[kernel]++;
			}
		}
		
		// Propagating
		out=0; top=-1;
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelPeriodic[i]!=0) stack[++top]=i;
		}
		while (top>=0) {
			kernel=stack[top--];
			type=pathKernelPeriodic[kernel];
			period=tmpArray2[kernel]==0?0:tmpArray1[kernel]/tmpArray2[kernel];
			for (j=0; j<=lastOutNeighbor[kernel]; j++) {
				neighbor=outNeighbors[kernel][j]; 
				valid=false;
				if ( pathKernelPeriodic[neighbor]==0 && 
				     (type==1 || (period>0 && pathKernelLengths[neighbor]>=period-LONG_PERIOD_THRESHOLD))
				   ) {
					valid=true; out++; top++;
					if (top==stack.length) {
						int[] newStack = new int[stack.length<<1];
						System.arraycopy(stack,0,newStack,0,stack.length);
						stack=newStack;
					}
					stack[top]=neighbor;
				}
				else if (pathKernelPeriodic[neighbor]==2 && type==1) {
					valid=true; top++;
					if (top==stack.length) {
						int[] newStack = new int[stack.length<<1];
						System.arraycopy(stack,0,newStack,0,stack.length);
						stack=newStack;
					}
					stack[top]=neighbor;
				}
				if (valid) {
					pathKernelPeriodic[neighbor]=addPeriodic[pathKernelPeriodic[neighbor]][type];
					if (pathKernelPeriodic[neighbor]==2 && period>0) {
						tmpArray1[neighbor]+=period;
						tmpArray2[neighbor]++;
					}
				}
			}
			if (pathKernelLengths[kernel]==0) continue;
			for (j=0; j<=lastInNeighbor[kernel]; j++) {
				neighbor=inNeighbors[kernel][j];
				p=j<<2;
				if ( pathKernelLengths[neighbor]!=0 &&
				     ((inNeighborFlags[kernel][p+1]==1 && inNeighborFlags[kernel][p+0]==1) || (inNeighborFlags[kernel][p+3]==1 && inNeighborFlags[kernel][p+2]==1))
				   ) {
					if (pathKernelPeriodic[neighbor]==0) {
						out++; top++;
						if (top==stack.length) {
							int[] newStack = new int[stack.length<<1];
							System.arraycopy(stack,0,newStack,0,stack.length);
							stack=newStack;
						}
						stack[top]=neighbor;
					}
					else if (pathKernelPeriodic[neighbor]==2 && type==1) {
						top++;
						if (top==stack.length) {
							int[] newStack = new int[stack.length<<1];
							System.arraycopy(stack,0,newStack,0,stack.length);
							stack=newStack;
						}
						stack[top]=neighbor;
					}
					pathKernelPeriodic[neighbor]=addPeriodic[pathKernelPeriodic[neighbor]][type];
					if (pathKernelPeriodic[neighbor]==2 && period>0) {
						tmpArray1[neighbor]+=period;
						tmpArray2[neighbor]++;
					}
				}
			}
		}
		return out;
	}
	
	
	/**
	 * Cell $(i,j)$ tells the addition of periodic type $j$ to periodic type $i$.
	 * Values: 0=not periodic, 1=short period, 2=long period.
	 */
	private static final byte[][] addPeriodic = new byte[][] {
		{0,1,2}, {1,1,1}, {2,1,2}
	};
	
	
	/**
	 * Resets $labels[i]$ for every periodic, possibly cyclic $i$: the new value is -3 if 
	 * the (path) kernel is selected as a representative, -2 if it is redundant. 
	 * Representatives are just maximal nodes in the subgraph of the kernel graph induced 
	 * by $pathKernelPeriodic$ and cyclic kernels (such nodes might not be maximal in the 
	 * whole kernel graph, since a periodic path kernel might be contained in a non-
	 * periodic path kernel), and their basins can have any size.
	 *
	 * Remark: the procedure removes any periodic or cyclic path kernel from the 
	 * $pathsWith*$ arrays.
	 *
	 * Remark: the procedure initializes $sorted2original,original2sorted$ and 
	 * $kernel2NewKernels,lastNewKernel$ if they are NULL.
	 *
	 * Remark: the procedure uses $tmpArray{1,2,3,4}$ and $tmpByte{1,2}$, which are 
	 * assumed to be of size at least $nPathKernels$ each.
	 * 
	 * @param all marks all periodic, possibly cyclic kernels as representative;
	 * @return the number of repeats marked with -3 in $labels$ (might be zero).
	 */
	public static final int markRepeats_pathKernelPeriodic(byte[] labels, boolean all) {
		final int GROWTH_RATE = 8;  // Arbitrary
		boolean found;
		byte arcOrientation;
		int i, j, k;
		int last, kernel, neighbor, nComponents, isDAG, out, last4;
		int[] components, lastOutNeighborPrime, minimalVertices, nVerticesPerComponent;
		
		if (all) {
			if (kernel2NewKernels==null) kernel2NewKernels = new int[nPathKernels][GROWTH_RATE];
			if (lastNewKernel==null) {
				lastNewKernel = new int[nPathKernels];
				Math.set(lastNewKernel,nPathKernels-1,-1);
			}
			if (kernel2NewKernels_orientation==null) kernel2NewKernels_orientation = new byte[nPathKernels][GROWTH_RATE];
			out=0;
			for (i=0; i<nPathKernels; i++) {
				if (pathKernelLengths[i]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[i]==0)) continue;
				labels[i]=(byte)(-3);
				out++;
				lastNewKernel[i]=0;
				kernel2NewKernels[i][0]=i;
				kernel2NewKernels_orientation[i][0]=0;
			}
			return out;
		}
		
		// Topologically sorting the kernel graph, if needed.
		if (sorted2original==null) sortKernelGraph();
		
		// Marking representatives
		out=0;
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[i]==0)) continue;
			found=false;
			for (j=0; j<=lastInNeighbor[i]; j++) {
				neighbor=inNeighbors[i][j];
				if (pathKernelLengths[neighbor]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[neighbor]!=0)) {
					found=true;
					break;
				}
			}
			if (found) labels[i]=(byte)(-2);
			else {
				labels[i]=(byte)(-3);
				out++;
			}
		}
		
		// Initializing arrays
		if (tmpArray1==null || tmpArray1.length<(nPathKernels<<1)) tmpArray1 = new int[nPathKernels<<1];
		if (tmpArray2==null || tmpArray2.length<(nPathKernels<<1)) tmpArray2 = new int[nPathKernels<<1];
		if (tmpArray3==null || tmpArray3.length<(nPathKernels<<1)) tmpArray3 = new int[nPathKernels<<1];
		if (tmpArray4==null || tmpArray4.length<nPathKernels) tmpArray4 = new int[nPathKernels<<1];
		if (tmpArray5==null || tmpArray5.length<nPathKernels) tmpArray5 = new int[nPathKernels<<1];
		if (tmpByte1==null || tmpByte1.length<(nPathKernels<<1)) tmpByte1 = new byte[nPathKernels<<1];
		if (tmpByte2==null || tmpByte2.length<(nPathKernels<<1)) tmpByte2 = new byte[nPathKernels<<1];
		
		// Updating $kernel2NewKernels$, and kernel tags in interval graph nodes.
		if (kernel2NewKernels==null) kernel2NewKernels = new int[nPathKernels][GROWTH_RATE];
		if (lastNewKernel==null) {
			lastNewKernel = new int[nPathKernels];
			Math.set(lastNewKernel,nPathKernels-1,-1);
		}
		if (kernel2NewKernels_orientation==null) kernel2NewKernels_orientation = new byte[nPathKernels][GROWTH_RATE];
		for (i=nPathKernels-1; i>=0; i--) {
			kernel=sorted2original[i];
			if (pathKernelLengths[kernel]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[kernel]==0)) continue;
			if (labels[kernel]==-3) {
				lastNewKernel[kernel]=0;
				kernel2NewKernels[kernel][0]=kernel;
				kernel2NewKernels_orientation[kernel][0]=0;
			}
			for (j=0; j<=lastOutNeighbor[kernel]; j++) {
				neighbor=outNeighbors[kernel][j];
				if (pathKernelLengths[neighbor]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[neighbor]==0)) continue;
				k=Arrays.binarySearch(inNeighbors[neighbor],0,lastInNeighbor[neighbor]+1,kernel);
				arcOrientation=inNeighborFlags2orientation[inNeighborFlags[neighbor][(k<<2)+1]][inNeighborFlags[neighbor][(k<<2)+3]];				
				if (lastNewKernel[neighbor]==-1) {
					last=lastNewKernel[kernel];
					if (last>=kernel2NewKernels[neighbor].length) kernel2NewKernels[neighbor] = new int[last+GROWTH_RATE];
					System.arraycopy(kernel2NewKernels[kernel],0,kernel2NewKernels[neighbor],0,last+1);
					lastNewKernel[neighbor]=last;
					if (last>=kernel2NewKernels_orientation[neighbor].length) kernel2NewKernels_orientation[neighbor] = new byte[last+GROWTH_RATE];
					System.arraycopy(kernel2NewKernels_orientation[kernel],0,kernel2NewKernels_orientation[neighbor],0,last+1);
					if (arcOrientation==1) {
						for (k=0; k<=last; k++) {
							if (kernel2NewKernels_orientation[neighbor][k]!=-1) kernel2NewKernels_orientation[neighbor][k]=oppositeOrientation[kernel2NewKernels_orientation[neighbor][k]];
						}
					}
					else if (arcOrientation==2) {
						for (k=0; k<=last; k++) kernel2NewKernels_orientation[neighbor][k]=2;
					}
				}
				else {
					last=Math.setUnion(kernel2NewKernels[kernel],lastNewKernel[kernel],kernel2NewKernels[neighbor],lastNewKernel[neighbor],tmpArray1);
					updateOrientations(kernel2NewKernels[kernel],lastNewKernel[kernel],kernel2NewKernels_orientation[kernel],arcOrientation,kernel2NewKernels[neighbor],lastNewKernel[neighbor],kernel2NewKernels_orientation[neighbor],tmpArray1,last,tmpByte1);
					if (last>=kernel2NewKernels_orientation[neighbor].length) kernel2NewKernels_orientation[neighbor] = new byte[last+GROWTH_RATE];
					System.arraycopy(tmpByte1,0,kernel2NewKernels_orientation[neighbor],0,last+1);
					if (last>=kernel2NewKernels[neighbor].length) kernel2NewKernels[neighbor] = new int[last+GROWTH_RATE];
					System.arraycopy(tmpArray1,0,kernel2NewKernels[neighbor],0,last+1);
					lastNewKernel[neighbor]=last;
				}
			}
		}
		for (i=0; i<IntervalGraph.nNodes; i++) markRepeats_updateKernels(IntervalGraph.nodesArray[i],false,kernel2NewKernels,lastNewKernel,kernel2NewKernels_orientation,null,null,null,null,-1);
		return out;
	}
	
	
	/**
	 * Overwrites $kernelFrequency[i]$ with the following estimates for every periodic or
	 * cyclic kernel $i$: 0=number of occurrences of the kernel in the genome; this is 
	 * half the number of times a maximal sequence of contained, straddling or adjacent 
	 * nodes tagged with kernel $i$ is adjacent to a high-quality substring of a read; 
	 * 1=min length and 2=max length of an occurrence; 3=number of interval graph nodes.
	 *
	 * Remark: the procedure assumes that $IntervalGraph.nodesArray$ is sorted by 
	 * $read,start$.
	 */
	public static final void getFrequencyPeriodic() throws IOException {
		int i, j;
		int lastSelectedKernel, from, previousReadA, readA, kernel, nCyclic;
		IntervalGraph.Node node;
		
		// Loading cyclic kernels
		if (tmpArray3==null || tmpArray3.length<nPathKernels) tmpArray3 = new int[nPathKernels];
		lastSelectedKernel=-1; nCyclic=0;
		for (i=0; i<nPathKernels; i++) {
			if (pathKernelLengths[i]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[i]==0)) continue;
			if (pathKernelLengths[i]==0) nCyclic++;
			tmpArray3[++lastSelectedKernel]=i;
			kernelFrequency[i][0]=0; 
			kernelFrequency[i][1]=Math.POSITIVE_INFINITY; 
			kernelFrequency[i][2]=0; kernelFrequency[i][3]=0;
		}
		if (lastSelectedKernel==-1) {
			System.err.println("No periodic or cyclic kernel found");
			return;
		}
		
		// Computing occurrences
		if (tmpArray1==null || tmpArray1.length<(lastSelectedKernel<<1)+1) tmpArray1 = new int[(lastSelectedKernel<<1)+1];
		if (tmpArray2==null || tmpArray2.length<lastSelectedKernel+1) tmpArray2 = new int[lastSelectedKernel+1];
		if (tmpArray4==null || tmpArray4.length<(lastSelectedKernel+1)*(Constants.INTERVAL_UNKNOWN+1)) {
			// One block of $Constants.INTERVAL_UNKNOWN+1$ elements per selected kernel.
			// Index $Constants.INTERVAL_UNKNOWN$ is used for long-period intervals.
			tmpArray4 = new int[(lastSelectedKernel+1)*(Constants.INTERVAL_UNKNOWN+1)];
		}
		from=0; previousReadA=IntervalGraph.nodesArray[0].read;
		for (i=1; i<IntervalGraph.nNodes; i++) {
			readA=IntervalGraph.nodesArray[i].read;
			if (readA!=previousReadA) {
				getFrequencyPeriodic_impl(previousReadA,from,i-1,tmpArray3,lastSelectedKernel);
				previousReadA=readA;
				from=i;
			}
		}
		getFrequencyPeriodic_impl(previousReadA,from,IntervalGraph.nNodes-1,tmpArray3,lastSelectedKernel);
		for (i=0; i<=lastSelectedKernel; i++) kernelFrequency[tmpArray3[i]][0]>>=1;
		
		// Printing statistics on cyclic kernels
		System.err.println(nCyclic+" cyclic kernels found:");
		for (i=0; i<=lastSelectedKernel; i++) {
			kernel=tmpArray3[i];
			if (pathKernelLengths[kernel]!=0) continue;
			System.err.print("kernelID="+pathKernel2Kernel[kernel]+" freq="+kernelFrequency[kernel][0]+" nNodes="+kernelFrequency[kernel][3]+" minLength="+kernelFrequency[kernel][1]+" maxLength="+kernelFrequency[kernel][2]+" :: Node types: ");
			System.err.print("A="+tmpArray4[i*(Constants.INTERVAL_UNKNOWN+1)+Constants.INTERVAL_ALIGNMENT]+",");
			System.err.print("P="+tmpArray4[i*(Constants.INTERVAL_UNKNOWN+1)+Constants.INTERVAL_DENSE_PREFIX]+",");
			System.err.print("S="+tmpArray4[i*(Constants.INTERVAL_UNKNOWN+1)+Constants.INTERVAL_DENSE_SUFFIX]+",");
			System.err.print("PS="+tmpArray4[i*(Constants.INTERVAL_UNKNOWN+1)+Constants.INTERVAL_DENSE_PREFIXSUFFIX]+",");
			System.err.print("SS="+tmpArray4[i*(Constants.INTERVAL_UNKNOWN+1)+Constants.INTERVAL_DENSE_SUBSTRING]+",");
			System.err.print("SD="+tmpArray4[i*(Constants.INTERVAL_UNKNOWN+1)+Constants.INTERVAL_DENSE_SINGLEDELETION]+",");
			System.err.print("PS="+tmpArray4[i*(Constants.INTERVAL_UNKNOWN+1)+Constants.INTERVAL_PERIODIC]+",");
			System.err.print("PL="+tmpArray4[i*(Constants.INTERVAL_UNKNOWN+1)+Constants.INTERVAL_UNKNOWN]);
			System.err.println();
		}		
	}
	
	
	/**
	 * Remark: the procedure uses temporary array $tmpArray1$ (respectively, $tmpArray2$),
	 * which is assumed to be at least as long as $2C$ (respectively, $C$), where $C$ is
	 * the number of cyclic kernels.
	 */
	private static final void getFrequencyPeriodic_impl(int read, int from, int to, int[] selectedKernels, int lastSelectedKernel) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int QUALITY_WINDOW = IO.quantum<<1;  // Arbitrary
		final int QUALITY_THRESHOLD = Reads.MIN_RANDOM_QUALITY_SCORE;
		boolean hasLongPeriod;
		int i, j, k;
		int lastIntersection, kernel, start, length, nodeType;
		final int readLength = Reads.getReadLength(read);
		IntervalGraph.Node node;
		
		length=(lastSelectedKernel+1)<<1;
		if (tmpArray1==null || tmpArray1.length<length) tmpArray1 = new int[length];
		Math.set(tmpArray1,length-1,-1);  // 2i=start, 2i+1=end.
		for (i=from; i<=to; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastKernel==-1) continue;
			nodeType=node.type;
			hasLongPeriod=nodeType==Constants.INTERVAL_PERIODIC&&node.hasLongPeriod;
			lastIntersection=Math.setIntersection(node.kernels,0,node.lastKernel,selectedKernels,0,lastSelectedKernel,tmpArray2,0);
			for (j=0; j<=lastIntersection; j++) {
				kernelFrequency[tmpArray2[j]][3]++;
				kernel=Arrays.binarySearch(selectedKernels,0,lastSelectedKernel+1,tmpArray2[j]);
				if (nodeType<Constants.INTERVAL_UNKNOWN) tmpArray4[kernel*(Constants.INTERVAL_UNKNOWN+1)+(hasLongPeriod?Constants.INTERVAL_UNKNOWN:nodeType)]++;
				k=kernel<<1;
				if (tmpArray1[k+1]==-1) {
					if (node.start>QUALITY_WINDOW && Histograms.getAverage(Reads.getQualityArray(read),node.start-QUALITY_WINDOW,node.start-1)<QUALITY_THRESHOLD) {
						kernelFrequency[tmpArray2[j]][0]++;
					}
					tmpArray1[k]=node.start; tmpArray1[k+1]=node.end;
				}
				else {
					if (node.start<=tmpArray1[k+1]+IDENTITY_THRESHOLD) tmpArray1[k+1]=Math.max(tmpArray1[k+1],node.end);
					else {
						length=tmpArray1[k+1]-tmpArray1[k]+1;
						kernelFrequency[tmpArray2[j]][1]=Math.min(kernelFrequency[tmpArray2[j]][1],length);
						kernelFrequency[tmpArray2[j]][2]=Math.max(kernelFrequency[tmpArray2[j]][2],length);
						tmpArray1[k]=node.start;
						if (Histograms.getAverage(Reads.getQualityArray(read),Math.max(tmpArray1[k+1]+1,node.start-QUALITY_WINDOW),node.start-1)<QUALITY_THRESHOLD) kernelFrequency[tmpArray2[j]][0]+=2;
						tmpArray1[k+1]=node.end;
					}
				}
			}
		}
		for (i=0; i<=lastSelectedKernel; i++) {
			kernel=selectedKernels[i];
			k=i<<1;
			if (tmpArray1[k+1]==-1) continue;
			if (tmpArray1[k+1]<readLength-QUALITY_WINDOW && Histograms.getAverage(Reads.getQualityArray(read),tmpArray1[k+1]+1,tmpArray1[k+1]+QUALITY_WINDOW)<QUALITY_THRESHOLD) {
				kernelFrequency[kernel][0]++;
			}
			length=tmpArray1[k+1]-tmpArray1[k]+1;
			kernelFrequency[kernel][1]=Math.min(kernelFrequency[kernel][1],length);
			kernelFrequency[kernel][2]=Math.max(kernelFrequency[kernel][2],length);
		}
	}
	
	
	
	
	
	
	
	
	// ---------------------------- FINALIZATION PROCEDURES ------------------------------
	
	/**
	 * Sets $bidirectedGraphNode=-1$ for every node in $nodesInKernel$ that does not 
	 * contain, in its $kernels$ array, any of the tags in $blockKernels[]$. This allows
	 * to use $bidirectedGraphNode$ for highlighting only the intervals in the path 
	 * kernels that survived previous filtering (see e.g. $printKernelTags_printWindow$).
	 */
	public static final void updateBidirectedGraphNodePointer() {
		int i;
		IntervalGraph.Node node;
		
		for (i=0; i<nNodesInKernel; i++) IntervalGraph.nodesArray[nodesInKernel[i]].visited=0;
		for (i=0; i<nNodesInKernel; i++) {
			node=IntervalGraph.nodesArray[nodesInKernel[i]];
			if (node.bidirectedGraphNode==-1) continue;
			if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,blockKernels[i])>=0) node.visited++;
		}
		for (i=0; i<nNodesInKernel; i++) {
			node=IntervalGraph.nodesArray[nodesInKernel[i]];
			if (node.bidirectedGraphNode==-1) continue;
			if (node.visited==0) node.bidirectedGraphNode=-1;
		}
	}
	
	
	/**
	 * The pipeline might not assign a kernel tag to some interval graph nodes, for 
	 * several reasons. E.g. an interval graph node that is used to build a kernel (and 
	 * thus belongs to $nodesInKernel$) might be removed from the bidirected graph by 
	 * $buildBidirectedGraph_{disconnectInsertions(),fixForks()}$, and then it might never
	 * be assigned to a kernel because it might not be reachable by containment.
	 * This procedure tries to assign a kernel tag to such nodes, using the matrix built
	 * by $buildNodes2Kernel()$ using the current kernels and basins.
	 *
	 * Remark: the procedure tries to avoid assigning a node to a kernel, if the node is a
	 * permutation of the kernel. However, this is done using simple and conservative 
	 * heuristics based on type and length of nodes and path kernels, rather than using 
	 * all alignments like in the simplification procedures called by 
	 * $buildKernelGraph()$. This is because such procedures cannot be easily adapted to 
	 * this case, and rewriting them would take considerable effort. So this procedure 
	 * might still assign a node to a kernel, if the node is a permutation or a deletion 
	 * of the kernel.
	 *
	 * Remark: some nodes might still have no kernel tag after this procedure completes.
	 *
	 * Remark: $buildBidirectedGraph_{mergeStraddling(),mergeIdentical()}$ remove nodes of
	 * the bidirected graph as well, but they map removed nodes to new nodes of the 
	 * bidirected graph; thus, removed nodes should get a kernel tag.
	 */
	private static final void nodesWithNoKernel(String alignmentsFilePath, String shortPeriodBitvectorPath, String shortPeriodTrackPath, int nAlignments) throws IOException {
		final int LENGTH_THRESHOLD = IO.quantum<<1;  // Arbitrary
		final double PERIOD_FRACTION = 0.1;  // Arbitrary
		final int GROWTH_RATE = 8;  // Arbitrary
		int i, j, k;
		int nNodes, nUnassigned, kernel, orientation, insertionDistance;
		int nPeriodic, kernelPeriod;
		final int previousOrder = IntervalGraph.Node.order;
		IntervalGraph.Node node;
		int[] lastNode2kernel, selected;
		IntervalGraph.Node[] nodes;
		int[][] nodes2kernel, avgPeriods;
		
		// Building $nodes2kernel$.
		nNodes=0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			if (IntervalGraph.nodesArray[i].lastKernel==-1) nNodes++;
		}
		if (nNodes==0) return;
		System.err.println("nodesWithNoKernel> Nodes without kernel at the beginning: "+nNodes+" ("+IO.getPercent(nNodes,IntervalGraph.nNodes)+"%)");
		nodes = new IntervalGraph.Node[nNodes];
		j=-1;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastKernel>=0) continue;
			nodes[++j]=node;
			if (node.kernels==null || node.kernels.length==0) node.kernels = new int[GROWTH_RATE];
			node.kernels[0]=nPathKernels;  // Artificial kernel tag
			if (node.kernelOrientations==null || node.kernelOrientations.length==0) node.kernelOrientations = new byte[GROWTH_RATE];
			node.kernelOrientations[0]=0; 
			node.lastKernel=0;
			if (node.pathsWithStart==null || node.pathsWithStart.length==0) node.pathsWithStart = new int[GROWTH_RATE];
			node.pathsWithStart[0]=node.kernels[0]; node.lastPathWithStart=0;
			if (node.pathsWithEnd==null || node.pathsWithEnd.length==0) node.pathsWithEnd = new int[GROWTH_RATE];
			node.pathsWithEnd[0]=-1-node.kernels[0]; node.lastPathWithEnd=0;
		}
		if (nNodes>1) {
			IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
			Arrays.sort(nodes,0,nNodes);
			IntervalGraph.Node.order=previousOrder;
		}
		if (pathKernelLengths==null) pathKernelLengths = new int[nPathKernels+1];
		else if (pathKernelLengths.length<nPathKernels+1) {
			int[] newLengths = new int[nPathKernels+1];
			System.arraycopy(pathKernelLengths,0,newLengths,0,nPathKernels);
			pathKernelLengths=newLengths;
		}
		pathKernelLengths[nPathKernels]=Math.POSITIVE_INFINITY;  // Artificial kernel tag
		if (pathKernelPeriodic==null) pathKernelPeriodic = new byte[nPathKernels+1];
		if (pathKernelPeriodic.length<nPathKernels+1) {
			byte[] newPeriodic = new byte[nPathKernels+1];
			System.arraycopy(pathKernelPeriodic,0,newPeriodic,0,nPathKernels);
			pathKernelPeriodic=newPeriodic;
		}
		pathKernelPeriodic[nPathKernels]=0;  // Artificial kernel tag
		if (shortPeriodAlignments==null) shortPeriodAlignments=loadShortPeriodBitvector(nAlignments,alignmentsFilePath,shortPeriodBitvectorPath);
		if (shortPeriodWindows==null) {
			shortPeriodWindows=IntervalGraph.loadShortPeriodTrack(shortPeriodTrackPath);
			trackPointers=IntervalGraph.buildTrackPointers(shortPeriodWindows);
		}
		nodes2kernel = new int[nNodes][0];
		lastNode2kernel = new int[nNodes];

if (IO.SHOW_INTERACTIVE) System.err.println("nodesWithNoKernel> building nodes2kernel...");

		buildNodes2Kernel(nodes,nNodes,alignmentsFilePath,nodes2kernel,lastNode2kernel,nAlignments);
		
if (IO.SHOW_INTERACTIVE) {	
System.err.println("nodesWithNoKernel> final nodes2kernel:");
for (int x=0; x<nNodes; x++) {
	System.err.println("node: "+nodes[x]);
	System.err.print("nodes2kernel: ");
	for (int y=0; y<=lastNode2kernel[x]; y++) System.err.print(nodes2kernel[x][y]+",");
	System.err.println();
}
}



		// Computing avg. periods
		nPeriodic=0; selected=null; avgPeriods=null;
		if (pathKernelPeriodic!=null) {
			for (i=0; i<nPathKernels; i++) {
				if (pathKernelPeriodic[i]!=0) nPeriodic++;
			}
			if (nPeriodic>=1) {
				selected = new int[nPeriodic];
				nPeriodic=-1;
				for (i=0; i<nPathKernels; i++) {
					if (pathKernelPeriodic[i]!=0) selected[++nPeriodic]=i;
				}
				nPeriodic++;
				if (nPeriodic>0) {
					avgPeriods = new int[nPeriodic][2];
					getAvgPeriods(selected,nPeriodic,avgPeriods);
				}
			}
		}
		
		// Assigning kernels using $nodes2kernel$.
		nUnassigned=0;
		for (i=0; i<nNodes; i++) {
			node=nodes[i];
			node.lastKernel=-1; node.lastPathWithStart=-1; node.lastPathWithEnd=-1;
			for (j=0; j<=lastNode2kernel[i]; j+=BLOCK_SIZE_PRIME) {
				kernel=nodes2kernel[i][j+1];
				orientation=nodes2kernel[i][j+3];
				kernelPeriod=pathKernelPeriodic[kernel]!=0?avgPeriods[Arrays.binarySearch(selected,0,nPeriodic,kernel)][0]:-1;
				if ( pathKernelLengths[kernel]==0 ||
					 ( (pathKernelPeriodic==null || pathKernelPeriodic[kernel]==0) && (node.type==Constants.INTERVAL_PERIODIC || node.length()>pathKernelLengths[kernel]+LENGTH_THRESHOLD) ) ||
					 ( pathKernelPeriodic!=null && pathKernelPeriodic[kernel]==2 && 
					   ( node.type!=Constants.INTERVAL_PERIODIC || !node.hasLongPeriod || 
						 (node.period>0 && Math.abs(node.period,kernelPeriod)>Math.max(kernelPeriod*PERIOD_FRACTION,LENGTH_THRESHOLD))
					   )
					 )
				   ) {
					// The assignment is not safe
					continue;
				}
				node.lastKernel++;
				if (node.lastKernel==node.kernels.length) {
					int[] newKernels = new int[node.kernels.length+GROWTH_RATE];
					System.arraycopy(node.kernels,0,newKernels,0,node.kernels.length);
					node.kernels=newKernels;
				}
				node.kernels[node.lastKernel]=kernel;
				if (node.lastKernel==node.kernelOrientations.length) {
					byte[] newOrientations = new byte[node.kernelOrientations.length+GROWTH_RATE];
					System.arraycopy(node.kernelOrientations,0,newOrientations,0,node.kernelOrientations.length);
					node.kernelOrientations=newOrientations;
				}
				node.kernelOrientations[node.lastKernel]=(byte)orientation;
				insertionDistance=nodes2kernel[i][j+4];
				if (insertionDistance<INSERTION_DISTANCE_THRESHOLD) {
					node.lastPathWithStart++;
					if (node.lastPathWithStart==node.pathsWithStart.length) {
						int[] newStarts = new int[node.pathsWithStart.length+GROWTH_RATE];
						System.arraycopy(node.pathsWithStart,0,newStarts,0,node.pathsWithStart.length);
						node.pathsWithStart=newStarts;
					}
					node.pathsWithStart[node.lastPathWithStart]=orientation==0?kernel:-1-kernel;
				}
				insertionDistance=nodes2kernel[i][j+5];
				if (insertionDistance<INSERTION_DISTANCE_THRESHOLD) {
					node.lastPathWithEnd++;
					if (node.lastPathWithEnd==node.pathsWithEnd.length) {
						int[] newEnds = new int[node.pathsWithEnd.length+GROWTH_RATE];
						System.arraycopy(node.pathsWithEnd,0,newEnds,0,node.pathsWithEnd.length);
						node.pathsWithEnd=newEnds;
					}
					node.pathsWithEnd[node.lastPathWithEnd]=orientation==0?-1-kernel:kernel;
				}
			}
			if (node.lastKernel==-1) {
				nUnassigned++;
				continue;
			}
			
			// Removing duplicated tags
			if (node.lastKernel>0) {
				k=0; kernel=node.kernels[0]; orientation=node.kernelOrientations[0];
				for (j=1; j<=node.lastKernel; j++) {
					if (node.kernels[j]!=kernel) {
						node.kernels[k]=kernel; node.kernelOrientations[k]=(byte)orientation;
						k++;
						kernel=node.kernels[j]; orientation=node.kernelOrientations[j];
						continue;
					}
					orientation=addOrientation[orientation][node.kernelOrientations[j]];
				}
				node.kernels[k]=kernel; node.kernelOrientations[k]=(byte)orientation;
				node.lastKernel=k;
			}
			if (node.lastPathWithStart>0) {
				Arrays.sort(node.pathsWithStart,0,node.lastPathWithStart+1);
				k=0; kernel=node.pathsWithStart[0];
				for (j=1; j<=node.lastPathWithStart; j++) {
					if (node.pathsWithStart[j]!=kernel) {
						k++; 
						node.pathsWithStart[k]=node.pathsWithStart[j];
						kernel=node.pathsWithStart[j];
					}
				}
				node.lastPathWithStart=k;
			}
			if (node.lastPathWithEnd>0) {
				Arrays.sort(node.pathsWithEnd,0,node.lastPathWithEnd+1);
				k=0; kernel=node.pathsWithEnd[0];
				for (j=1; j<=node.lastPathWithEnd; j++) {
					if (node.pathsWithEnd[j]!=kernel) {
						k++; 
						node.pathsWithEnd[k]=node.pathsWithEnd[j];
						kernel=node.pathsWithEnd[j];
					}
				}
				node.lastPathWithEnd=k;
			}
		}
		System.err.println("nodesWithNoKernel> Nodes without kernel at the end: "+nUnassigned+" ("+IO.getPercent(nUnassigned,IntervalGraph.nNodes)+"%)");
	}
	
	
	/**
	 * Removes from every interval graph node all kernel tags that are assigned to fewer 
	 * than $threshold$ nodes. Nodes that have no kernel tag after this are not assigned 
	 * to surviving kernel tags.
	 *
	 * Remark: the procedure assumes $kernelSize$ to be up to date.
	 *
	 * @return the number of surviving kernel tags.
	 */
	public static final int discardRareKernels(int threshold) {
		int i, j, k;
		int last1, last2, last3;
		IntervalGraph.Node node;
		
		// Loading rare kernel tags
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		last1=-1;
		for (i=0; i<nPathKernels; i++) {
			if (kernelSize[i]<threshold) tmpArray1[++last1]=i;
		}
		if (last1==-1) return nPathKernels;
		last2=(last1+1)<<1;
		if (tmpArray2==null || tmpArray2.length<last2) tmpArray2 = new int[last2];
		System.arraycopy(tmpArray1,0,tmpArray2,last1+1,last1+1);
		for (i=0; i<=last1; i++) tmpArray2[last1-i]=-1-tmpArray2[last1+1+i];
		if (tmpArray3==null || tmpArray3.length<nPathKernels) tmpArray3 = new int[nPathKernels];
		
		// Removing kernel tags from interval graph nodes
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastKernel>=0) {
				last3=Math.setMinus(node.kernels,node.lastKernel,tmpArray1,last1,tmpArray3);
				if (last3>=0 && last3<node.lastKernel) {
					k=0;
					for (j=0; j<=last3; j++) {
						while (node.kernels[k]<tmpArray3[j]) k++;
						node.kernelOrientations[j]=node.kernelOrientations[k];
					}
					System.arraycopy(tmpArray3,0,node.kernels,0,last3+1);
				}
				node.lastKernel=last3;
			}
			if (node.lastPathWithStart>=0) {
				last3=Math.setMinus(node.pathsWithStart,node.lastPathWithStart,tmpArray2,last2,tmpArray3);
				if (last3>=0 && last3<node.lastPathWithStart) System.arraycopy(tmpArray3,0,node.pathsWithStart,0,last3+1);
				node.lastPathWithStart=last3;
			}
			if (node.lastPathWithEnd>=0) {
				last3=Math.setMinus(node.pathsWithEnd,node.lastPathWithEnd,tmpArray2,last2,tmpArray3);
				if (last3>=0 && last3<node.lastPathWithEnd) System.arraycopy(tmpArray3,0,node.pathsWithEnd,0,last3+1);
				node.lastPathWithEnd=last3;
			}
		}
		System.out.println("discardRareKernels> Discarded "+(last1+1)+" kernels with fewer than "+threshold+" nodes ("+IO.getPercent(last1+1,nPathKernels)+"%)");
		return last1+1;
	}
	
	
	/**
	 * Checks basic consistency criteria on $kernels$, $kernelOrientations$, $pathsWith*$
	 * at the very end of the pipeline.
	 *
	 * @param kernelLabels possible values are those set by $markRepeats()$
	 * and $markRepeats_pathKernelPeriodic()$.
	 */
	private static final void checkKernelTags(byte[] kernelLabels) {
		int i, j;
		int kernel;
		IntervalGraph.Node node;
		
		for (i=0; i<IntervalGraph.nNodes; i++) {
			// Kernels and orientations
			node=IntervalGraph.nodesArray[i];
			if (node.lastKernel==-1) continue;
			if (node.kernelOrientations==null || node.kernelOrientations.length<node.lastKernel+1) {
				System.err.println("checkKernelTags> ERROR: the following node has "+(node.lastKernel+1)+" kernels but "+(node.kernelOrientations==null?"0":node.kernelOrientations.length)+" kernel orientations?!  "+node);
				System.exit(1);
			}
			for (j=0; j<=node.lastKernel; j++) {
				if (kernelLabels[node.kernels[j]]>0 || kernelLabels[node.kernels[j]]==-3) continue;
				System.err.println("checkKernelTags> ERROR: not all kernel tags of the following node are to be returned in output:");
				System.err.println("checkKernelTags> Node: "+node);
				System.err.print("checkKernelTags> Tags selected for output: ");
				for (int x=0; x<nPathKernels; x++) {
					if (kernelLabels[x]>0 || kernelLabels[x]==-3) System.err.print(x+"("+kernelLabels[x]+"), ");
				}
				System.err.println();
				System.exit(1);
			}
			for (j=1; j<=node.lastKernel; j++) {
				if (node.kernels[j]==node.kernels[j-1]) {
					System.err.print("checkKernelTags> ERROR: duplicate kernel");
					System.err.println("checkKernelTags> Node: "+node);
					System.err.print("checkKernelTags> Kernel tags:");
					for (int x=0; x<=node.lastKernel; x++) System.err.print(node.kernels[x]+",");
					System.err.println();
					System.err.print("checkKernelTags> Kernel orientations:");
					for (int x=0; x<=node.lastKernel; x++) System.err.print(node.kernelOrientations[x]+",");
					System.err.println();
					System.exit(1);
				}
			}
			
			// Start
			for (j=0; j<=node.lastPathWithStart; j++) {
				kernel=node.pathsWithStart[j];
				if (kernel<0) kernel=-1-kernel;
				if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,kernel)<0) {
					System.err.println("checkKernelTags> ERROR: the following node has kernel tag "+node.pathsWithStart[j]+" in pathsWithStart[], which does not occur in kernels[]:");
					System.err.println("checkKernelTags> Node: "+node);
					System.err.print("checkKernelTags> pathsWithStart: ");
					for (int x=0; x<=node.lastPathWithStart; x++) System.err.print(node.pathsWithStart[x]+",");
					System.err.println();
					System.err.print("checkKernelTags> kernelLabels:");
					for (int x=0; x<nPathKernels; x++) System.err.print(x+"("+kernelLabels[x]+"), ");
					System.err.println();
					System.exit(1);
				}
				if ((pathKernelLengths!=null && pathKernelLengths[kernel]==0) || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
					System.err.println("checkKernelTags> ERROR: the following node has kernel tag "+node.pathsWithStart[j]+" in pathsWithStart[], which is cyclic or periodic:");
					System.err.println("checkKernelTags> Node: "+node);
					System.err.print("checkKernelTags> pathsWithStart: ");
					for (int x=0; x<=node.lastPathWithStart; x++) System.err.print(node.pathsWithStart[x]+",");
					System.err.println();
					System.exit(1);
				}
			}
			// End
			for (j=0; j<=node.lastPathWithEnd; j++) {
				kernel=node.pathsWithEnd[j];
				if (kernel<0) kernel=-1-kernel;
				if (Arrays.binarySearch(node.kernels,0,node.lastKernel+1,kernel)<0) {
					System.err.println("checkKernelTags> ERROR: the following node has kernel tag "+node.pathsWithEnd[j]+" in pathsWithEnd[], which does not occur in kernels[]:");
					System.err.println("checkKernelTags> Node: "+node);
					System.err.print("checkKernelTags> pathsWithEnd: ");
					for (int x=0; x<=node.lastPathWithEnd; x++) System.err.print(node.pathsWithEnd[x]+",");
					System.err.println();
					System.exit(1);
				}
				if ((pathKernelLengths!=null && pathKernelLengths[kernel]==0) || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]!=0)) {
					System.err.println("checkKernelTags> ERROR: the following node has kernel tag "+node.pathsWithEnd[j]+" in pathsWithEnd[], which is cyclic or periodic:");
					System.err.println("checkKernelTags> Node: "+node);
					System.err.print("checkKernelTags> pathsWithEnd: ");
					for (int x=0; x<=node.lastPathWithEnd; x++) System.err.print(node.pathsWithEnd[x]+",");
					System.err.println();
					System.exit(1);
				}
			}
		}
	}
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	// ------------------------------ PRINTING PROCEDURES --------------------------------
	
	/**
	 * Remark: the procedure prints non-periodic, periodic, and cyclic path kernels.
	 *
	 * @param labels discarded if null; possible values are those set by $markRepeats()$
	 * and $markRepeats_pathKernelPeriodic()$;
	 * @param frequencies frequencies of every non-periodic path kernel, according to the 
	 * conventions of $markRepeats()$;
	 * @param kernelNames ignored if null.
	 */
	public static final void printKernelGraph(String path, byte[] labels, int[][] frequencies, int[][] kernelNames) throws IOException {
		final int PEN_WIDTH = 4;
		int i, j;
		int from, penWidth;
		String fillColor;
		final String forwardColor = "#509FB5";
		final String rcColor = "#2F8C6D";
		final String[] labelColors = new String[] {"#5B201D", "#A6231A","#B66B15", "#319FC7","#244254","#2F8C6D"};
		final String[] labelColors_periodic = new String[] {"#00E89B","#FEED50"};
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(path),IO.BUFFER_SIZE);
		bw.write("digraph G {\n");
		bw.write("bgcolor=\""+Colors.COLOR_DOT_BACKGROUND+"\";\n");
		for (i=0; i<nPathKernels; i++) {
			if (labels==null) fillColor=labelColors[0];
			else {
				if (labels[i]<-1) fillColor=labelColors_periodic[-labels[i]-2];
				else if (labels[i]==-1) fillColor=labelColors[0];
				else fillColor=labelColors[labels[i]];
			}
			bw.write(i+" [shape=\"plaintext\",label=\"("+(kernelNames==null?i:IO.printRow(kernelNames[i]))+",K="+(pathKernel2Kernel==null?"?":pathKernel2Kernel[i])+",L="+(pathKernelLengths!=null&&pathKernelLengths[i]==0?(pathKernelLongestNode!=null?pathKernelLongestNode[i]:"?"):(pathKernelLengths==null?"?":pathKernelLengths[i]))+") f="+(labels!=null&&labels[i]>=0?frequencies[i][0]:-1)+",p="+(labels!=null&&labels[i]>=0?frequencies[i][1]:-1)+",s="+(labels!=null&&labels[i]>=0?frequencies[i][2]:-1)+"\",style=\"filled\",color=\""+fillColor+"\",fontcolor=\""+Colors.COLOR_DOT_TEXT+"\",fontname=\"Helvetica\",fontsize=\"15\"];\n");
		}
		for (i=0; i<nPathKernels; i++) {
			for (j=0; j<=lastInNeighbor[i]; j++) {
				from=inNeighbors[i][j];
				if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0) || pathKernelLengths[from]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[from]==1)) {
					if (inNeighborFlags[i][(j<<2)+1]==1) bw.write(from+" -> "+i+" [style=\"dotted\",penwidth=\""+PEN_WIDTH+"\",color=\""+forwardColor+"\"];\n");
					if (inNeighborFlags[i][(j<<2)+3]==1) bw.write(from+" -> "+i+" [style=\"dotted\",penwidth=\""+PEN_WIDTH+"\",color=\""+rcColor+"\"];\n");
				}
				else {
					if (inNeighborFlags[i][(j<<2)+1]==1) bw.write(from+" -> "+i+" [style=\""+(inNeighborFlags[i][(j<<2)+0]>1?"bold":"dotted")+"\",penwidth=\""+PEN_WIDTH+"\",color=\""+forwardColor+"\",inNeighborFlags=\""+inNeighborFlags[i][(j<<2)+0]+"\"];\n");
					if (inNeighborFlags[i][(j<<2)+3]==1) bw.write(from+" -> "+i+" [style=\""+(inNeighborFlags[i][(j<<2)+2]>1?"bold":"dotted")+"\",penwidth=\""+PEN_WIDTH+"\",color=\""+rcColor+"\",inNeighborFlags=\""+inNeighborFlags[i][(j<<2)+2]+"\"];\n");
				}				
			}
		}
		bw.write("}\n");
		bw.close();
	}
	
	
	/**
	 * Appends to file $path$ the kernel tags of every node in $IntervalGraph.nodesArray$,
	 * using the format specified by $printKernelTags_printWindow()$. This procedure is 
	 * designed for drawing a graphical plot, thus it tries to simplify the output by 
	 * merging overlapping nodes that carry similar information, and specifically nodes 
	 * with no kernel tag and nodes with identical kernel tags: see 
	 * $printKernelTags_shouldBeMerged()$ for details on merging criteria.
	 *
	 * Remark: the procedure appends to a file, so it might be necessary to sort the file 
	 * after this procedure is called on it multiple times.
	 *
	 * Remark: some intervals printed by this procedure might have no corresponding 
	 * interval produced by factorization, since e.g. $IntervalGraph.fixUndersplits()$ can
	 * create new intervals.
	 *
	 * @param aggressive merging mode: see $printKernelTags_shouldBeMerged()$ for details;
	 * @param discardNoTag discards nodes with no kernel tag;
	 * @param artificialKernel defined in $printKernelTags_printWindow()$.
	 */
	public static final void printKernelTags(String path, String prefix, boolean aggressive, boolean discardNoTag, int artificialKernel) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int CAPACITY = 4;  // Arbitrary
		boolean merged;
		int i, j, k;
		int currentRead, inputRead, inputFrom, lastActiveWindow, nTag, nNoTag;
		int outputFrom, outputFromPrime, lastOutput;
		final int previousOrder = IntervalGraph.Node.order;
		IntervalGraph.Node node, tmpNode;
		IntervalGraph.Node[] activeWindows, output;
		BufferedWriter bw;
		
		// Building $output$.
		output = new IntervalGraph.Node[CAPACITY];
		activeWindows = new IntervalGraph.Node[CAPACITY];
		node=IntervalGraph.nodesArray[0];
		inputRead=node.read; inputFrom=0; outputFrom=0; lastOutput=-1;
		if (node.lastKernel!=-1) {  // Active windows have a kernel tag
			lastActiveWindow=0; printKernelTags_setArray(activeWindows,0,node);
			nTag=1; nNoTag=0;
		}
		else {
			lastActiveWindow=-1;
			nTag=0; nNoTag=1;
		}	
		for (i=1; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.read!=inputRead) {
				// Flushing all active windows
				if (lastActiveWindow>=0) {
					output=printKernelTags_ensureArray(output,lastOutput+1+lastActiveWindow+1,CAPACITY);
					for (j=0; j<=lastActiveWindow; j++) printKernelTags_setArray(output,++lastOutput,activeWindows[j]);
				}
				// Merging windows with no kernel tag, to some output window with a kernel
			 	// tag, and with one another.
				if (!discardNoTag && nNoTag>0) {
					outputFromPrime=lastOutput+1;
					if (nTag>0) {
						output=printKernelTags_mergeNoTag_tag(inputFrom,i-1,output,outputFrom,lastOutput,CAPACITY,tmpIO);
						lastOutput=tmpIO[0];
						if (lastOutput>outputFromPrime) lastOutput=printKernelTags_mergeNoTag_noTag(output,outputFromPrime,lastOutput,IDENTITY_THRESHOLD);
					}
					else {
						output=printKernelTags_ensureArray(output,lastOutput+1+i-inputFrom,CAPACITY);
						for (j=inputFrom; j<i; j++) printKernelTags_setArray(output,++lastOutput,IntervalGraph.nodesArray[j]);
						if (lastOutput>outputFromPrime) lastOutput=printKernelTags_mergeNoTag_noTag(output,outputFromPrime,lastOutput,IDENTITY_THRESHOLD);
					}
				}
				// Next iteration
				inputRead=node.read; inputFrom=i; outputFrom=lastOutput+1;
				if (node.lastKernel!=-1) {
					lastActiveWindow=0; printKernelTags_setArray(activeWindows,0,node);
					nTag=1; nNoTag=0;
				}
				else {
					lastActiveWindow=-1;
					nTag=0; nNoTag=1;
				}
				continue;
			}
			if (node.lastKernel==-1) {
				nNoTag++;
				continue;
			}
			nTag++;
			// Trying to merge (merging with multiple active windows is allowed).
			merged=false;
			for (j=lastActiveWindow; j>=0; j--) {
				if (printKernelTags_shouldBeMerged(node,activeWindows[j],IDENTITY_THRESHOLD,pathKernelLengths,pathKernelPeriodic,aggressive)) {
					printKernelTags_merge(node,activeWindows[j],aggressive,IDENTITY_THRESHOLD);
					merged=true;
				}
			}
			// Moving inactive windows to $output$.
			k=-1;
			for (j=0; j<=lastActiveWindow; j++) {
				if (activeWindows[j].end<node.start-IDENTITY_THRESHOLD) {
					lastOutput++;
					output=printKernelTags_ensureArray(output,lastOutput+1,CAPACITY);
					printKernelTags_setArray(output,lastOutput,activeWindows[j]);
				}
				else {
					k++;
					tmpNode=activeWindows[k]; 
					activeWindows[k]=activeWindows[j];
					activeWindows[j]=tmpNode;
				}
			}
			lastActiveWindow=k;
			// Building a new active window from $node$.
			if (!merged) {
				lastActiveWindow++;
				activeWindows=printKernelTags_ensureArray(activeWindows,lastActiveWindow+1,CAPACITY);
				printKernelTags_setArray(activeWindows,lastActiveWindow,node);
			}
		}
		// Flushing all active windows
		if (lastActiveWindow>=0) {
			output=printKernelTags_ensureArray(output,lastOutput+1+lastActiveWindow+1,CAPACITY);
			for (j=0; j<=lastActiveWindow; j++) printKernelTags_setArray(output,++lastOutput,activeWindows[j]);
		}
		// Merging windows with no kernel tag
		if (!discardNoTag && nNoTag>0) {
			outputFromPrime=lastOutput+1;
			if (nTag>0) {
				output=printKernelTags_mergeNoTag_tag(inputFrom,IntervalGraph.nNodes-1,output,outputFrom,lastOutput,CAPACITY,tmpIO);
				lastOutput=tmpIO[0];
				if (lastOutput>outputFromPrime) lastOutput=printKernelTags_mergeNoTag_noTag(output,outputFromPrime,lastOutput,IDENTITY_THRESHOLD);
			}
			else {
				output=printKernelTags_ensureArray(output,lastOutput+1+IntervalGraph.nNodes-inputFrom,CAPACITY);
				for (j=inputFrom; j<IntervalGraph.nNodes; j++) printKernelTags_setArray(output,++lastOutput,IntervalGraph.nodesArray[j]);
				if (lastOutput>outputFromPrime) lastOutput=printKernelTags_mergeNoTag_noTag(output,outputFromPrime,lastOutput,IDENTITY_THRESHOLD);
			}
		}
		
		// Sorting the windows of each read, removing duplicates, and printing $output$
		// to file.
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		j=0; currentRead=output[0].read;
		for (i=1; i<=lastOutput; i++) {
			if (output[i].read==currentRead) continue;
			if (i>j+1) {
				Arrays.sort(output,j,i);
				printKernelTags_mergeIdenticalWindows(output,j,i-1);
			}
			j=i; currentRead=output[i].read;
		}
		if (lastOutput>=j+1) {
			Arrays.sort(output,j,lastOutput+1);
			printKernelTags_mergeIdenticalWindows(output,j,lastOutput);
		}
		IntervalGraph.Node.order=previousOrder;
		bw = new BufferedWriter(new FileWriter(path,true/*Appending*/),IO.BUFFER_SIZE);
		for (i=0; i<=lastOutput; i++) {
			if (output[i]!=null) printKernelTags_printWindow(output[i],prefix,bw,artificialKernel);
		}
		bw.close();
	}

	
	public static final IntervalGraph.Node[] printKernelTags_ensureArray(IntervalGraph.Node[] array, int newLength, int growthRate) {
		final int oldLength = array.length;
		if (oldLength>=newLength) return array;
		IntervalGraph.Node[] newArray = new IntervalGraph.Node[newLength+growthRate];
		System.arraycopy(array,0,newArray,0,oldLength);
		return newArray;
	}
	
	
	public static final void printKernelTags_setArray(IntervalGraph.Node[] array, int p, IntervalGraph.Node source) {
		if (array[p]==null) array[p] = new IntervalGraph.Node();
		printKernelTags_clone(source,array[p]);
	}
	
	
	/**
	 * A version of $IntervalGraph.Node.clone()$ that copies from $node$ to $window$ just
	 * the information used by $printKernelTags()$.
	 */
	public static final void printKernelTags_clone(IntervalGraph.Node node, IntervalGraph.Node window) {
		window.read=node.read; window.start=node.start; window.end=node.end;
		window.bidirectedGraphNode=node.bidirectedGraphNode==-1?-1:0;
		window.lastKernel=node.lastKernel;
		window.lastPathWithStart=node.lastPathWithStart;
		window.lastPathWithEnd=node.lastPathWithEnd;
		if (node.lastKernel==-1) return;
		if (window.kernels==null || window.kernels.length<node.lastKernel+1) window.kernels = new int[node.lastKernel+1];
		System.arraycopy(node.kernels,0,window.kernels,0,node.lastKernel+1);
		if (window.kernelOrientations==null || window.kernelOrientations.length<node.lastKernel+1) window.kernelOrientations = new byte[node.lastKernel+1];
		System.arraycopy(node.kernelOrientations,0,window.kernelOrientations,0,node.lastKernel+1);
		if (node.lastPathWithStart!=-1) {
			if (window.pathsWithStart==null || window.pathsWithStart.length<node.lastPathWithStart+1) window.pathsWithStart = new int[node.lastPathWithStart+1];
			System.arraycopy(node.pathsWithStart,0,window.pathsWithStart,0,node.lastPathWithStart+1);
		}
		if (node.lastPathWithEnd!=-1) {
			if (window.pathsWithEnd==null || window.pathsWithEnd.length<node.lastPathWithEnd+1) window.pathsWithEnd = new int[node.lastPathWithEnd+1];
			System.arraycopy(node.pathsWithEnd,0,window.pathsWithEnd,0,node.lastPathWithEnd+1);
		}
	}
	
	
	/**
	 * Collapses all elements of $output[from..to]$ that have almost identical start/end 
	 * positions and that can be merged, setting to null all duplicates except one.
	 * Identical windows can be present at the end of $printKernelTags()$ since there can
	 * be incompatible sets of windows that give rise to the same interval when merged.
	 *
	 * @param output $[from..to]$ contains all nodes in the same read, sorted by $start$.
	 */
	private static final void printKernelTags_mergeIdenticalWindows(IntervalGraph.Node[] output, int from, int to) {
		final int IDENTITY_THRESHOLD = IO.quantum>>1;  // Arbitrary
		int i, j;
		int start, end;
		
		for (i=from; i<to; i++) {
			if (output[i]==null) continue;
			start=output[i].start; end=output[i].end;
			for (j=i+1; j<=to; j++) {
				if (output[j]==null) continue;
				if (output[j].start>start+IDENTITY_THRESHOLD) break;
				if (Math.abs(output[j].end,end)<=IDENTITY_THRESHOLD && printKernelTags_mergeIdenticalWindows_impl(output[i],output[j])) output[j]=null;
			}
		}
	}
	
	
	/**
	 * @return true if merging $window2$ into $window1$ was successful, i.e. if their
	 * $lastPathWith*$ arrays are identical.
	 */
	public static final boolean printKernelTags_mergeIdenticalWindows_impl(IntervalGraph.Node window1, IntervalGraph.Node window2) {
		int last;
		
		if ( !Math.setIdentity(window1.pathsWithStart,0,window1.lastPathWithStart,window2.pathsWithStart,0,window2.lastPathWithStart) ||
		     !Math.setIdentity(window1.pathsWithEnd,0,window1.lastPathWithEnd,window2.pathsWithEnd,0,window2.lastPathWithEnd) 
		   ) return false;
		if (window2.bidirectedGraphNode!=-1) window1.bidirectedGraphNode=0;
		last=window1.lastKernel+window2.lastKernel+1;
		if (last>=tmpArray1.length) tmpArray1 = new int[last+1];
		last=Math.setUnion(window1.kernels,window1.lastKernel,window2.kernels,window2.lastKernel,tmpArray1);
		if (last>=tmpByte1.length) tmpByte1 = new byte[last+1];
		updateOrientations(window1.kernels,window1.lastKernel,window1.kernelOrientations,(byte)0,window2.kernels,window2.lastKernel,window2.kernelOrientations,tmpArray1,last,tmpByte1);
		if (window1.kernels==null || last>=window1.kernels.length) window1.kernels = new int[last+1];
		System.arraycopy(tmpArray1,0,window1.kernels,0,last+1);		
		if (window1.kernelOrientations==null || last>=window1.kernelOrientations.length) window1.kernelOrientations = new byte[last+1];
		System.arraycopy(tmpByte1,0,window1.kernelOrientations,0,last+1);
		window1.lastKernel=last;
		window1.end=Math.max(window1.end,window2.end);
		return true;
	}
	
	
	/**
	 * Merges windows with no kernel tag, to windows with a kernel tag. Windows that 
	 * cannot be merged are appended to $output$.
	 *
	 * @param from,to all intervals in the same read in $IntervalGraph.nodesArray$;
	 * @param fromOutput all intervals in the same read in $output$ are contained in
	 * suffix $[fromOutput..lastOutput]$;
	 * @param out output array: cell 0 contains the new value of $lastOutput$ at the end
	 * of the procedure.
	 */
	private static final IntervalGraph.Node[] printKernelTags_mergeNoTag_tag(int from, int to, IntervalGraph.Node[] output, int fromOutput, int lastOutput, int capacity, int[] out) {
		boolean merged;
		int i, j;
		final int toOutput = lastOutput;
		IntervalGraph.Node tmpNode;
		
		for (i=from; i<=to; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			if (tmpNode.lastKernel!=-1) continue;
			merged=false;
			for (j=fromOutput; j<=toOutput; j++) {
				if (Intervals.isContained(tmpNode.start,tmpNode.end,output[j].start,output[j].end)) {
					merged=true;
					break;
				}
			}
			if (!merged) {
				lastOutput++;
				output=printKernelTags_ensureArray(output,lastOutput+1,capacity);
				printKernelTags_setArray(output,lastOutput,tmpNode);
			}
		}
		out[0]=lastOutput;
		return output;
	}
	
	
	/**
	 * Merges windows with no kernel tag, to other windows with no kernel tag. Both the 
	 * input windows and the result of the merge are located at the end of $output$.
	 *
	 * @param from,to suffix of $output$ that contains all windows with no tag (the 
	 * procedure assumes $to>from$);
	 * @return the last window of $output$ after the merge.
	 */ 
	private static final int printKernelTags_mergeNoTag_noTag(IntervalGraph.Node[] output, int from, int to, int identityThreshold) {
		final int previousOrder = IntervalGraph.Node.order;
		int i, j;
		IntervalGraph.Node tmpNode;
		
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		Arrays.sort(output,from,to+1);
		IntervalGraph.Node.order=previousOrder;
		j=from;
		for (i=from+1; i<=to; i++) {
			if (Intervals.isContained(output[i].start,output[i].end,output[j].start,output[j].end)) continue;
			else if (output[i].start<=output[j].end+identityThreshold) {
				output[j].end=Math.max(output[j].end,output[i].end);
				continue;
			}
			j++;
			tmpNode=output[j]; output[j]=output[i]; output[i]=tmpNode;
		}
		return j;
	}
	
	
	/**
	 * Remark: we might not want to merge adjacent short-period intervals with the same 
	 * kernel ID but very different coverage, or between which there is a drop in 
	 * coverage, since it could mean that there is a change in period offset or sequence.
	 * This is not implemented yet. Different orientations are already not merged.
	 *
	 * @param node,window assumed to have at least one kernel tag each; the procedure 
	 * assumes also that $node.start>=window.start$;
	 * @param aggressive when both $node$ and $window$ have multiple kernel tags, allows
	 * merging them even if the sets of tags are different;
	 * @return TRUE iff $node$ can be merged to $window$.
	 */
	public static final boolean printKernelTags_shouldBeMerged(IntervalGraph.Node node, IntervalGraph.Node window, int identityThreshold, int[] pathKernelLengths, byte[] pathKernelPeriodic, boolean aggressive) {
		final int IDENTITY_THRESHOLD = IO.quantum>>1;  // Arbitrary
		final boolean closeToStart = Math.abs(node.start,window.start)<=identityThreshold;
		final boolean closeToEnd = Math.abs(node.end,window.end)<=identityThreshold;
		final boolean sameOneKernel = window.lastKernel==0 && node.lastKernel==0 && window.kernels[0]==node.kernels[0] && window.kernelOrientations[0]==node.kernelOrientations[0];
		final int nodeEnd, windowStart, lastKernel;
		int i;
		int kernel, kernelLength;
		
		if (Math.abs(node.start,window.start)<=IDENTITY_THRESHOLD && Math.abs(node.end,window.end)<=IDENTITY_THRESHOLD) { 
			// Merging never disallowed
		}
		else {
			if ( (!aggressive && !window.sameKernels(node)) || 
				  (aggressive && !( window.sameKernels(node) || (node.lastKernel>0 && window.lastKernel>0) ))
			   ) return false;
		}
		if (Intervals.isContained(node.start,node.end,window.start,window.end)) {
			if (sameOneKernel) return true;  // Self-similarity of the module
			if ((!closeToStart && node.lastPathWithStart!=-1) || (!closeToEnd && node.lastPathWithEnd!=-1)) return false;
			if ( (closeToStart && Math.nonemptySetMinus(node.pathsWithStart,node.lastPathWithStart,window.pathsWithStart,window.lastPathWithStart)) ||
				 (closeToEnd && Math.nonemptySetMinus(node.pathsWithEnd,node.lastPathWithEnd,window.pathsWithEnd,window.lastPathWithEnd))
			   ) return false;
			return true;
		}
		else if (closeToStart) {
			if (sameOneKernel) return true;  // Self-similarity of the module
			if (Math.nonemptySetMinus(node.pathsWithStart,node.lastPathWithStart,window.pathsWithStart,window.lastPathWithStart)) return false;
			if ( (closeToEnd && Math.nonemptySetMinus(node.pathsWithEnd,node.lastPathWithEnd,window.pathsWithEnd,window.lastPathWithEnd)) ||
				 (!closeToEnd && window.lastPathWithEnd!=-1)
			   ) return false;
			return true;
		}
		else if (node.start<window.end-identityThreshold) {
			if (sameOneKernel) return true;  // Self-similarity of the module
			if (node.lastPathWithStart!=-1) return false;
			if ( (closeToEnd && Math.nonemptySetMinus(node.pathsWithEnd,node.lastPathWithEnd,window.pathsWithEnd,window.lastPathWithEnd)) ||
				 (!closeToEnd && window.lastPathWithEnd!=-1) 
			   ) return false;
			return true;
		}
		else if (Math.abs(window.end,node.start)<=identityThreshold) {
			if (window.lastPathWithEnd!=-1 || node.lastPathWithStart!=-1) return false;
			if (pathKernelLengths!=null) {
				nodeEnd=node.end; windowStart=window.start;
				lastKernel=node.lastKernel;
				for (i=0; i<=lastKernel; i++) {
					kernel=node.kernels[i]; kernelLength=pathKernelLengths[kernel];
					if (kernelLength==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[kernel]==1)) continue;
					if ((pathKernelPeriodic!=null && pathKernelPeriodic[kernel]==2) || Math.abs(kernelLength,nodeEnd-windowStart+1)>identityThreshold) return false;
				}
			}
			return true;
		}
		else return false;
	}
	
	
	/**
	 * Merges $node$ to $window$, which are assumed to satisfy
	 * $printKernelTags_shouldBeMerged()$.
	 */
	public static final void printKernelTags_merge(IntervalGraph.Node node, IntervalGraph.Node window, boolean aggressive, int identityThreshold) {
		int last;
		
		if (node.bidirectedGraphNode!=-1) window.bidirectedGraphNode=0;
		if (Math.abs(node.end,window.end)<=identityThreshold) {
			last=Math.setUnion(window.pathsWithEnd,window.lastPathWithEnd,node.pathsWithEnd,node.lastPathWithEnd,tmpArray1);
			if (window.pathsWithEnd==null || window.pathsWithEnd.length<last+1) window.pathsWithEnd = new int[last+1];
			System.arraycopy(tmpArray1,0,window.pathsWithEnd,0,last+1);
			window.lastPathWithEnd=last;
		}
		else if (node.end>window.end) {
			if (node.lastPathWithEnd!=-1) {
				if (window.pathsWithEnd==null || window.pathsWithEnd.length<node.lastPathWithEnd+1) window.pathsWithEnd = new int[node.lastPathWithEnd+1];
				System.arraycopy(node.pathsWithEnd,0,window.pathsWithEnd,0,node.lastPathWithEnd+1);
			}
			window.lastPathWithEnd=node.lastPathWithEnd;
		}
		window.end=Math.max(window.end,node.end);
		if (aggressive) {
			last=node.lastKernel+window.lastKernel+1;
			if (last>=tmpArray1.length) tmpArray1 = new int[last+1];
			last=Math.setUnion(node.kernels,node.lastKernel,window.kernels,window.lastKernel,tmpArray1);
			if (last>=tmpByte1.length) tmpByte1 = new byte[last+1];
			updateOrientations(window.kernels,window.lastKernel,window.kernelOrientations,(byte)0,node.kernels,node.lastKernel,node.kernelOrientations,tmpArray1,last,tmpByte1);
			if (window.kernels==null || last>=window.kernels.length) window.kernels = new int[last+1];
			System.arraycopy(tmpArray1,0,window.kernels,0,last+1);
			if (window.kernelOrientations==null || last>=window.kernelOrientations.length) window.kernelOrientations = new byte[last+1];
			System.arraycopy(tmpByte1,0,window.kernelOrientations,0,last+1);
			window.lastKernel=last;
		}
	}

	
	/**
	 * Format:
	 *
	 * read,start,end,f,x,y,z,prefix : k_0,...k_x : o_0,...,o_x : l_0,...,l_x : 
	 * s_0,...,s_y : e_0,...,e_z
	 *
	 * where: 
	 * - $f \in {0,1}$ signals if the window is used in the reference sequence of a path 
	 * kernel;
	 * - $k_i$ is the ID of path kernel $i'$ (just a number, without $prefix$);
	 * - $o_i$ is the orientation of the window WRT $k_i$;
	 * - $l_i$ is the string length of the full sequence of $k_i$;
	 * - $s_i$ is the ID of an end of path kernel $i'$ that is aligned to the beginning of
	 * the window (either $i'$ or $-1-i'$);
	 * - $e_i$ is the ID of an end of path kernel $i'$ that is aligned to the end of the
	 * window (either $i'$ or $-1-i'$).
	 *
	 * @param artificialKernel (ignored if -1) if there is just one kernel tag, and such 
	 * tag equals $artificialKernel$, the window is printed as if there were no kernel 
	 * tags.
	 */
	private static final void printKernelTags_printWindow(IntervalGraph.Node window, String prefix, BufferedWriter bw, int artificialKernel) throws IOException {
		final boolean artificial = artificialKernel>=0 && window.lastKernel==0 && window.kernels[0]==artificialKernel;
		int i;
		
		bw.write(window.read+",");
		bw.write(window.start+",");
		bw.write(window.end+",");
		bw.write((window.bidirectedGraphNode==-1?0:1)+",");
		bw.write(artificial?"-1,":window.lastKernel+",");
		bw.write(window.lastPathWithStart+",");
		bw.write(window.lastPathWithEnd+",");
		bw.write(prefix); bw.write(":");
		if (window.lastKernel==-1 || artificial) {
			bw.write("::::\n");
			return;
		}
		bw.write(window.kernels[0]+"");
		for (i=1; i<=window.lastKernel; i++) bw.write(","+window.kernels[i]);
		bw.write(":");
		bw.write(window.kernelOrientations[0]+"");
		for (i=1; i<=window.lastKernel; i++) bw.write(","+window.kernelOrientations[i]);
		bw.write(":");
		if (pathKernelLengths==null) {
			bw.write("0");
			for (i=1; i<=window.lastKernel; i++) bw.write(",0");
		}
		else {
			bw.write(pathKernelLengths[window.kernels[0]]+"");
			for (i=1; i<=window.lastKernel; i++) bw.write(","+pathKernelLengths[window.kernels[i]]);
		}
		bw.write(":");
		if (window.lastPathWithStart==-1) bw.write(":");
		else {
			bw.write(window.pathsWithStart[0]+"");
			for (i=1; i<=window.lastPathWithStart; i++) bw.write(","+window.pathsWithStart[i]);
			bw.write(":");
		}
		if (window.lastPathWithEnd!=-1) {
			bw.write(window.pathsWithEnd[0]+"");
			for (i=1; i<=window.lastPathWithEnd; i++) bw.write(","+window.pathsWithEnd[i]);
		}
		bw.write("\n");
	}
	
	
	/**
	 * Loads in $windows[0..X]$ and $prefixes[0..X]$ the set of all kernel tags in file 
	 * $br$ that belong to $read$, where $X$ is returned in output. The procedure starts 
	 * reading from the current position of $br$.
	 *
	 * @param windows output array: each row contains all the elements printed by
	 * $printKernelTags_printWindow()$, excluding $read$ and $prefix$ (which is stored in 
	 * $prefixes$).
	 */
	public static final int loadKernelTags(int read, BufferedReader br, int readAheadLimit, int[][] windows, String[] prefixes) throws IOException {
		int i, j, p, q;
		int n, start, end, f, x, y, z;
		int currentRead, lastWindow;
		String str;
		
		lastWindow=-1;
		str=br.readLine();
		while (str!=null) {
			p=str.indexOf(",");
			currentRead=Integer.parseInt(str.substring(0,p));
			if (currentRead<read) {
				str=br.readLine();
				continue;
			}
			else if (currentRead>read) {
				br.reset();
				return lastWindow;
			}
			br.mark(readAheadLimit);
			lastWindow++;
			// Loading $windows$ and $prefixes$.
			n=0;
			for (i=0; i<str.length(); i++) {
				if (str.charAt(i)==',' || str.charAt(i)==':') n++;
			}
			n++;
			n-=2;  // Removing $read$ and $prefix$.
			if (windows[lastWindow]==null || windows[lastWindow].length<n) windows[lastWindow] = new int[n];
			q=str.indexOf(",",p+1);
			start=Integer.parseInt(str.substring(p+1,q));
			p=q; q=str.indexOf(",",p+1);
			end=Integer.parseInt(str.substring(p+1,q));
			p=q; q=str.indexOf(",",p+1);
			f=Integer.parseInt(str.substring(p+1,q));
			p=q; q=str.indexOf(",",p+1);
			x=Integer.parseInt(str.substring(p+1,q));
			p=q; q=str.indexOf(",",p+1);
			y=Integer.parseInt(str.substring(p+1,q));
			p=q; q=str.indexOf(",",p+1);
			z=Integer.parseInt(str.substring(p+1,q));
			p=q; q=str.indexOf(":",p+1);
			prefixes[lastWindow]=str.substring(p+1,q);
			windows[lastWindow][0]=start; windows[lastWindow][1]=end;
			windows[lastWindow][2]=f; windows[lastWindow][3]=x; windows[lastWindow][4]=y;
			windows[lastWindow][5]=z;
			if (x==-1) {
				str=br.readLine();
				continue;
			}
			i=5;
			for (j=0; j<x; j++) {
				p=q; q=str.indexOf(",",p+1);
				windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q)); // k_i
			}
			p=q; q=str.indexOf(":",p+1);
			windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q)); // k_i
			for (j=0; j<x; j++) {
				p=q; q=str.indexOf(",",p+1);
				windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q));  // o_i
			}
			p=q; q=str.indexOf(":",p+1);
			windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q));  // o_i
			for (j=0; j<x; j++) {
				p=q; q=str.indexOf(",",p+1);
				windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q));  // l_i
			}
			p=q; q=str.indexOf(":",p+1);
			windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q));  // l_i
			if (y>=0) {  // s_i
				for (j=0; j<y; j++) {
					p=q; q=str.indexOf(",",p+1);
					windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q));  
				}
				p=q; q=str.indexOf(":",p+1);
				windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q));
			}
			else q=str.indexOf(":",q+1);
			if (z>=0) {  // e_i
				for (j=0; j<z; j++) {
					p=q; q=str.indexOf(",",p+1);
					windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1,q));
				}
				p=q;
				windows[lastWindow][++i]=Integer.parseInt(str.substring(p+1));
			}
			str=br.readLine();
		}
		return lastWindow;
	}
	
	
	/**
	 * Prints to a distinct file the description of every kernel tag selected in 
	 * $kernelLabels$. Such a file has the following format:
	 *
	 * K,L,P,T,Q,N,F_0,F_1,F_2,F_3
	 * r_i,s_i,e_i,o_i,ks_i,ke_i,f_i,t_i,p_i,ms_i,me_i,ok_i
	 * 
	 * where:
	 *
	 * the first line describes the path kernel: 
	 * - $K$ is the ID of the kernel of which the path kernel is a path; 
	 * - $L$ is the string length of the path kernel (zero if $K$ is cyclic); 
	 * - $P$ is the ID of the path kernel among all paths of kernel graph $K$ (-1 if $K$ 
	 *   is cyclic); 
	 * - $T$ is the type of repeat of the kernel tag (defined in $Factorize.java$); 
	 * - $Q$ is the length of the period if the path kernel is periodic and if a period 
	 *   can be estimated, 0 otherwise; 
	 * - $N$ is the number of lines that follow in the file; 
	 * - $F_i$ are estimates of the number of full, prefix, suffix, substring 
	 *   occurrences if the path kernel is neither periodic nor cyclic; otherwise, they 
	 *   are estimates of the number of distinct occurrences in the genome, and of the 
	 *   length of a shortest (longest) occurrence ($F_3=0$ is not used in this case).
	 *
	 * there is one $*_i$ record for each interval graph node in the basin of the path
	 * kernel, storing: 
	 * - read, start, end of the substring ($r_i,s_i,e_i$), and its orientation WRT the 
	 *   path kernel ($o_i$: 0=forward; 1=RC; 2=both); 
	 * - $ks_i$ is one (respectively, two) if the beginning (respectively, end) of the 
	 *   path kernel occurs at the start of the interval, and zero otherwise ($ke_i$ is
	 *   similarly defined for the end of the interval graph node).
	 * - $f_i$ is the ID of the node of the bidirected graph to which the interval graph 
	 *   node is mapped (-1 if it is not mapped); 
	 * - $t_i$ is the type of the interval (defined in $Factorize.java$); 
	 * - $p_i$ is the length of the period if the interval is periodic and if a period can
	 *   be estimated, 0 otherwise;
	 * - $ms_i$ (respectively, $me_i$) is one iff the start (end) of the interval is 
	 *   adjacent to a high-quality substring of a read;
	 * - $ok_i$ is one iff the interval belongs to other path kernels.
	 * 
	 * Remark: the procedure assumes $IntervalGraph.nodesArray$ to be sorted by $read,
	 * start$, and it outputs lines sorted by $r_i,s_i$. The procedure assumes that 
	 * $kernelSize$ and $kernelFrequency$ are correctly initialized.
	 *
	 * Remark: the procedure uses $tmpArray1$.
	 * 
	 * @param kernelLabels possible values are those set by $markRepeats()$ and 
	 * $markRepeats_pathKernelPeriodic()$;
	 * @param filePrefix each file has the form $filePrefix-*.txt$;
	 * @param kernel2descriptor filenames have the form $filePrefix-X.txt$, where $X$ is 
	 * the ID of a path kernel (if $kernel2descriptor=null$), or row $Y$ of 
	 * $kernel2descriptor$, where $Y$ is the ID of a path kernel (if $kernel2descriptor$ 
	 * is not null);
	 * @param forceShortPeriod makes sure that the basin decriptor header does not store a
	 * long period.
	 */
	public static final void printBasinDescriptors(byte[] kernelLabels, String outputDir, String filePrefix, int[][] kernel2descriptor, boolean forceShortPeriod) throws IOException {
		final int LENGTH_THRESHOLD = IO.quantum<<2;  // Same value as in $getBasins_ensureKernelLengths()$
		final int COUNTS_THRESHOLD = IO.coverage*IO.minOccurrencesInGenome;
		final double PERIODIC_THRESHOLD = 0.5;  // Arbitrary
		int i, j, k, h;
		int nSelected, kernel, kernelInSelected, last1, type, period, nodeLength;
		IntervalGraph.Node node;
		BufferedWriter bw;
		boolean[] hasReferenceNode;
		int[] selected, nPrintedLines;
		BufferedWriter[] files;
		int[][] counts = new int[nPathKernels][8];  
		
		if (tmpArray1==null || tmpArray1.length<nPathKernels) tmpArray1 = new int[nPathKernels];
		
		// Initializing output files
		nSelected=0;
		for (i=0; i<nPathKernels; i++) {
			if ((kernelLabels[i]>0 || kernelLabels[i]==-3) && kernelSize[i]>0) nSelected++;
		}
		if (nSelected==0) return;
		selected = new int[nSelected];
		j=-1;
		for (i=0; i<nPathKernels; i++) {
			if ((kernelLabels[i]>0 || kernelLabels[i]==-3) && kernelSize[i]>0) selected[++j]=i;
		}
		printBasinDescriptors_getCounts(selected,nSelected,counts);		
		files = new BufferedWriter[nSelected];
		for (i=0; i<nSelected; i++) files[i] = new BufferedWriter(new FileWriter(outputDir+"/"+filePrefix+"-"+(kernel2descriptor==null?selected[i]:IO.printRow(kernel2descriptor[selected[i]]))+".txt"),IO.BUFFER_SIZE);
		
		// Writing path kernel info
		k=-1;
		for (i=0; i<nPathKernels; i++) {
			if ((kernelLabels[i]>0 || kernelLabels[i]==-3) && kernelSize[i]>0) {
				k++;
				type=printBasinDescriptors_getType(i,counts,COUNTS_THRESHOLD,PERIODIC_THRESHOLD);
				period=(type==Constants.INTERVAL_PERIODIC&&((!forceShortPeriod)||counts[i][6]<PeriodicSubstrings.MIN_PERIOD_LONG))?counts[i][6]:0;
				kernel=pathKernel2Kernel[i];
				if (pathKernelLengths[i]==0) files[k].write(kernel+",0,-1,"+type+","+period+","+kernelSize[i]);
				else {
					j=i;
					while (j>=0 && pathKernel2Kernel[j]==kernel) j--;
					j++;
					files[k].write(kernel+","+pathKernelLengths[i]+","+(i-j)+","+type+","+period+","+kernelSize[i]);
				}
				for (j=0; j<=2; j++) files[k].write(","+kernelFrequency[i][j]);
				if (pathKernelLengths[i]==0 || (pathKernelPeriodic!=null && pathKernelPeriodic[i]!=0)) files[k].write(",0\n");
				else files[k].write(","+kernelFrequency[i][3]+"\n");
			}
		}
		
		// Writing basin info
		hasReferenceNode = new boolean[nSelected];
		Math.set(hasReferenceNode,hasReferenceNode.length-1,false);
		nPrintedLines = new int[nPathKernels];
		Math.set(nPrintedLines,nPathKernels-1,0);
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastKernel==-1) continue;
			nodeLength=node.length();
			last1=Math.setIntersection(node.kernels,0,node.lastKernel,selected,0,nSelected-1,tmpArray1,0);
			if (last1==-1) {
				// Only selected kernels should be present in the interval graph at this
				// stage of the pipeline.
				System.err.println("printBasinDescriptors> ERROR: the following node has some kernel tags, but none of them is to be returned in output?!");
				System.err.println("printBasinDescriptors> Node: "+node);
				System.err.print("printBasinDescriptors> Tags selected for output: ");
				for (int x=0; x<nSelected; x++) System.err.print(selected[x]+",");
				System.err.println();
				System.exit(1);
			}
			k=0;
			for (j=0; j<=last1; j++) {
				kernel=tmpArray1[j];
				if (pathKernelLengths[kernel]!=0 && (pathKernelPeriodic==null || pathKernelPeriodic[kernel]==0) && nodeLength>pathKernelLengths[kernel]+LENGTH_THRESHOLD) {
					System.err.println("printBasinDescriptors> ERROR: the following node has length "+nodeLength+" > "+pathKernelLengths[kernel]+"=pathKernelLength of path kernel "+kernel);
					System.err.println("printBasinDescriptors> "+node);
					System.exit(1);
				}
				while (k<=node.lastKernel && node.kernels[k]<kernel) k++;
				kernelInSelected=Arrays.binarySearch(selected,0,nSelected,kernel);
				bw=files[kernelInSelected];
				bw.write(node.read+","+node.start+","+node.end+","+node.kernelOrientations[k]+",");
				nPrintedLines[kernel]++;
				if (node.lastPathWithStart>=0) {
					if (Arrays.binarySearch(node.pathsWithStart,0,node.lastPathWithStart+1,kernel)>=0) bw.write("1,");
					else if (Arrays.binarySearch(node.pathsWithStart,0,node.lastPathWithStart+1,-1-kernel)>=0) bw.write("2,");
					else bw.write("0,");
				}
				else bw.write("0,");
				if (node.lastPathWithEnd>=0) {
					if (Arrays.binarySearch(node.pathsWithEnd,0,node.lastPathWithEnd+1,kernel)>=0) bw.write("1,");
					else if (Arrays.binarySearch(node.pathsWithEnd,0,node.lastPathWithEnd+1,-1-kernel)>=0) bw.write("2,");
					else bw.write("0,");
				}
				else bw.write("0,");
				bw.write(node.bidirectedGraphNode+",");
				if (node.bidirectedGraphNode!=-1) hasReferenceNode[kernelInSelected]=true;
				bw.write(node.type+",");
				bw.write((node.type==Constants.INTERVAL_PERIODIC?node.period:0)+",");
				bw.write(node.isLeftMaximal?"1,":"0,");
				bw.write(node.isRightMaximal?"1,":"0,");
				bw.write(last1>0?"1,":"0,");
				bw.write("\n");
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<nPathKernels; i++) {
				if (Arrays.binarySearch(selected,0,nSelected,i)>=0 && nPrintedLines[i]!=kernelSize[i]) {
					System.err.println("printBasinDescriptors> ERROR: number of printed lines does not match kernelSize: "+nPrintedLines[i]+"!="+kernelSize[i]+" path kernel="+i);
					System.exit(1);
				}
			}
			for (i=0; i<nSelected; i++) {
				if (!hasReferenceNode[i]) {
					System.err.println("printBasinDescriptors> ERROR: the basin descriptor file of kernel "+selected[i]+" has no interval marked as in the reference.");
					System.exit(1);
				}
			}
		}
		
		// Finalizing
		for (i=0; i<nSelected; i++) {
			files[i].close(); files[i]=null;
		}
		files=null;
	}
	
	
	/**
	 * @param selected the procedure cumulates counts just for kernel tags in this array, 
	 * which is assumed to be sorted;
	 * @param counts output array: each row stores the following information: 
	 * 0=nFull, 1=nPrefix, 2=nSuffix, 3=nSubstring, 4=nPeriodic, 5=nPositivePeriod, 
	 * 6=avgPeriod, 7=nTotal$.
	 */
	private static final void printBasinDescriptors_getCounts(int[] selected, int nSelected, int[][] counts) {
		boolean atStart, atEnd;
		int i, j, k;
		int kernel, last1;
		IntervalGraph.Node node;
		
		Math.set(counts,0);
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastKernel==-1) continue;
			last1=Math.setIntersection(node.kernels,0,node.lastKernel,selected,0,nSelected-1,tmpArray1,0);
			if (last1==-1) {
				// Only selected kernels should be present in the interval graph at this
				// stage of the pipeline.
				System.err.println("printBasinDescriptors_getCounts> ERROR: the following node has some kernel tags, but none of them is to be returned in output?!");
				System.err.println("printBasinDescriptors_getCounts> Node: "+node);
				System.err.print("printBasinDescriptors_getCounts> Tags selected for output: ");
				for (int x=0; x<nSelected; x++) System.err.print(selected[x]+",");
				System.err.println();
				System.exit(1);
			}
			k=0;
			for (j=0; j<=last1; j++) {
				kernel=tmpArray1[j];
				while (k<=node.lastKernel && node.kernels[k]<kernel) k++;
				counts[kernel][7]++;
				atStart=(node.lastPathWithStart>=0 && Arrays.binarySearch(node.pathsWithStart,0,node.lastPathWithStart+1,kernel)>=0) ||
					    (node.lastPathWithEnd>=0 && Arrays.binarySearch(node.pathsWithEnd,0,node.lastPathWithEnd+1,kernel)>=0);
				atEnd=(node.lastPathWithStart>=0 && Arrays.binarySearch(node.pathsWithStart,0,node.lastPathWithStart+1,-1-kernel)>=0) ||
					  (node.lastPathWithEnd>=0 && Arrays.binarySearch(node.pathsWithEnd,0,node.lastPathWithEnd+1,-1-kernel)>=0);
				if (atStart) {
					if (atEnd) counts[kernel][0]++;
					else counts[kernel][1]++;
				}
				else {
					if (atEnd) counts[kernel][2]++;
					else counts[kernel][3]++;
				}
				if (node.type==Constants.INTERVAL_PERIODIC) {				
					counts[kernel][4]++;
					if (node.period>0) {
						counts[kernel][5]++;
						counts[kernel][6]+=node.period;
					}
				}
			}
		}
		for (i=0; i<nSelected; i++) {
			kernel=selected[i];
			if (counts[kernel][6]!=0 && counts[kernel][5]!=0) counts[kernel][6]/=counts[kernel][5];
		}
	}
	
	
	/**
	 * @param counts output of $printBasinDescriptors_getCounts()$;
	 * @param threshold on the number of intervals stored in $counts$, to decide a type;
	 * @return a code in $Constants.INTERVAL_*$.
	 */
	private static final int printBasinDescriptors_getType(int kernelTag, int[][] counts, int threshold, double periodicThreshold) {
		if ((pathKernelPeriodic!=null && pathKernelPeriodic[kernelTag]!=0) || counts[kernelTag][4]>=counts[kernelTag][7]*periodicThreshold) return Constants.INTERVAL_PERIODIC;
		if (counts[kernelTag][0]>=threshold && counts[kernelTag][1]<threshold && counts[kernelTag][2]<threshold && counts[kernelTag][3]<threshold) return Constants.INTERVAL_ALIGNMENT;
		if (counts[kernelTag][1]>=threshold && counts[kernelTag][2]<threshold && counts[kernelTag][3]<threshold) return Constants.INTERVAL_DENSE_PREFIX;
		if (counts[kernelTag][1]<threshold && counts[kernelTag][2]>=threshold && counts[kernelTag][3]<threshold) return Constants.INTERVAL_DENSE_SUFFIX;
		if (counts[kernelTag][1]>=threshold && counts[kernelTag][2]>=threshold && counts[kernelTag][3]<threshold) return Constants.INTERVAL_DENSE_PREFIXSUFFIX;
		if (counts[kernelTag][3]>=threshold) return Constants.INTERVAL_DENSE_SUBSTRING;
		return Constants.INTERVAL_UNKNOWN;
	}
	
	
	/**
	 * @param out the values of the header described in $printBasinDescriptors()$.
	 */
	public static final void readBasinDescriptorHeader(String str, int[] out) {
		int i, p, q;
		int kernelID, pathLength, pathID, type, period, nIntervals;
		
		p=str.indexOf(",");
		kernelID=Integer.parseInt(str.substring(0,p));
		q=str.indexOf(",",p+1);
		pathLength=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		pathID=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		type=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		period=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		nIntervals=Integer.parseInt(str.substring(p+1,q));
		
		out[0]=kernelID; out[1]=pathLength; out[2]=pathID;
		out[3]=type; out[4]=period; out[5]=nIntervals;
		for (i=0; i<=2; i++) {
			p=q; q=str.indexOf(",",p+1);
			out[6+i]=Integer.parseInt(str.substring(p+1,q));
		}
		out[9]=Integer.parseInt(str.substring(q+1));
	}
	
	
	/**
	 * @param out the values of a record described in $printBasinDescriptors()$.
	 */
	public static final void readBasinDescriptorRecord(String str, int[] out) {
		int p, q;
		int read, start, end, orientation, atStart, atEnd, bidirectedGraphNode;
		int type, period, leftMaximal, rightMaximal, otherKernels;
		
		p=str.indexOf(",");
		read=Integer.parseInt(str.substring(0,p));
		q=str.indexOf(",",p+1);
		start=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		end=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		orientation=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		atStart=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		atEnd=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		bidirectedGraphNode=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		type=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		period=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		leftMaximal=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		rightMaximal=Integer.parseInt(str.substring(p+1,q));
		p=q; q=str.indexOf(",",p+1);
		if (q>=0) {
			otherKernels=Integer.parseInt(str.substring(p+1,q));
			p=q; q=str.indexOf(",",p+1);
		}
		else otherKernels=0;
		
		out[0]=read; out[1]=start; out[2]=end; out[3]=orientation; out[4]=atStart; 
		out[5]=atEnd; out[6]=bidirectedGraphNode; out[7]=type; out[8]=period; out[9]=leftMaximal;
		out[10]=rightMaximal; out[11]=otherKernels;
	}

	
	/**
	 * Checks whether the string length of every node is consistent with the string 
	 * length of every path kernel assigned to it. See $getBasins_ensureKernelLengths()$
	 * for details.
	 *
	 * @param errorID,fileID just for display purposes;
	 * @param logFile if not NULL and if there is a node with excess length, writes to the
	 * file rather than to STDERR.
	 */
	private static final void checkKernelLengths(int errorID, String fileID, String logFile) throws IOException {
		final int LENGTH_THRESHOLD = IO.quantum<<2;  // Arbitrary
		final int N_BINS = 10;  // Arbitrary
		final int QUANTUM = Alignments.minAlignmentLength;
		int i, j;
		int length, kernel, diff, maxDiff, maxDiffPrime, maxNode, maxKernel;
		IntervalGraph.Node node;
		int[] histogram = new int[N_BINS];
		
		Math.set(histogram,histogram.length-1,0);
		maxDiff=0; maxNode=-1; maxKernel=-1;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			length=node.length();
			maxDiffPrime=0;
			for (j=0; j<=node.lastKernel; j++) {
				kernel=node.kernels[j];
				if ((pathKernelPeriodic==null || pathKernelPeriodic[kernel]==0) && pathKernelLengths[kernel]!=0 && length>pathKernelLengths[kernel]+LENGTH_THRESHOLD) {
					diff=length-pathKernelLengths[kernel];
					if (diff>maxDiffPrime) maxDiffPrime=diff;
					if (diff>maxDiff) { maxDiff=diff; maxNode=i; maxKernel=kernel; }
				}
			}
			if (maxDiffPrime!=0) histogram[Math.min(maxDiffPrime/QUANTUM,N_BINS-1)]++;
		}
		
		if (maxNode==-1) System.err.println("checkKernelLengths> "+errorID+"  all nodes have consistent lengths");
		else {
			if (logFile==null) {
				System.err.println("checkKernelLengths> "+errorID+"  histogram of excess lengths:  (% of all nodes)");
				for (i=0; i<N_BINS-1; i++) System.err.println("["+(QUANTUM*i)+".."+(QUANTUM*(i+1))+"): "+(100*((double)histogram[i])/IntervalGraph.nNodes)+"%");
				System.err.println("["+(QUANTUM*(i+1))+"..: "+(100*((double)histogram[i])/IntervalGraph.nNodes)+"%");
				System.err.println("checkKernelLengths> a node with max excess "+maxDiff+" (WRT kernel "+maxKernel+", length="+pathKernelLengths[maxKernel]+"): "+IntervalGraph.nodesArray[maxNode]);
			}
			else {
				BufferedWriter bw = new BufferedWriter(new FileWriter(logFile,true/*appending*/));
				bw.write(fileID+", "+errorID+", histogram of excess lengths:  (% of all nodes)\n");
				for (i=0; i<N_BINS-1; i++) bw.write("["+(QUANTUM*i)+".."+(QUANTUM*(i+1))+"): "+(100*((double)histogram[i])/IntervalGraph.nNodes)+"%\n");
				bw.write("["+(QUANTUM*(i+1))+"..: "+(100*((double)histogram[i])/IntervalGraph.nNodes)+"%\n");
				bw.write("A node with max excess "+maxDiff+" (WRT kernel "+maxKernel+", length="+pathKernelLengths[maxKernel]+"): "+IntervalGraph.nodesArray[maxNode]+"\n");
				bw.write("\n");
				bw.close();
			}
		}
	}
	
	
	
	
	

	
	
	
	
	
	
	// --------------------------------- TYPES OF EDGE -----------------------------------
	
	/**
	 * Remark: the graph of a periodic component can contain nonperiodic intervals.
	 * E.g. it could happen that a fragment $V$ of a periodic substring is incorporated
	 * by a longer nonperiodic repeat $U$. The instance of $V$ inside $U$ should have a 
	 * pattern of alignments that marks it as a periodic substring, and such alignments 
	 * should be assigned to $V$ rather than to $U$. If this does not happen, i.e. if $V$ 
	 * is not marked as a periodic substring, but has alignments with substrings of a 
	 * periodic substring interval, and such alignments are assigned to $U$, we have a 
	 * periodic-nonperiodic edge.
	 *
	 * It could also happen that an alignment makes a nonperiodic interval completely 
	 * contained inside a periodic interval (this could happen when a degenerate substring
	 * of a periodic substring becomes an independent repeat, with not enough similarity 
	 * to other substrings or instances of the periodic substring), and that an alignment 
	 * makes a nonperiodic interval coincide with a periodic interval (could be explained 
	 * as above).
	 *
	 * Remark: in the largest periodic component of Drosophila melanogaster (the only one 
	 * with clusters), overlap edges suitable for assembly are significantly more frequent
	 * than overlap edges induced by periodic-periodic alignments, and span the entire 
	 * graph.
	 *
	 * Remark: an edge between two periodic intervals can correspond to a large number of 
	 * alignments implied by those two intervals.
	 */
	public static final int inStep3(IntervalGraph.AlignmentPair pair) {
		int i;

		i=inStep3_containment(pair);
		if (i>=0) return i;
		i=inStep3_overlap(pair);
		if (i>=0) return i;
		i=IntervalGraphStep2.isInsertion(pair);
		if (i>=0) return i;
		i=IntervalGraphStep2.isSharedSubstring(pair);
		if (i>=0) return i;
		return -1;
	}
	
	
	private static final int inStep3_containment(IntervalGraph.AlignmentPair pair) {
		// Identity
		if ( Math.abs(pair.alignmentStart1,pair.node1.start)<=IntervalGraph.IDENTITY_THRESHOLD && 
		     Math.abs(pair.alignmentEnd1-1,pair.node1.end)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(pair.alignmentStart2,pair.node2.start)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(pair.alignmentEnd2-1,pair.node2.end)<=IntervalGraph.IDENTITY_THRESHOLD
		   ) return Constants.CONTAINMENT_IDENTICAL;
		
		// Proper containment
		if (pair.node1.type!=Constants.INTERVAL_PERIODIC && pair.node2.type!=Constants.INTERVAL_PERIODIC) return IntervalGraphStep2.inStep2_containment(pair);
		return IntervalGraphStep2.isContained_alignment_alignment(pair);
	}
	
	
	private static final int inStep3_overlap(IntervalGraph.AlignmentPair pair) {
		final int type1 = pair.node1.type;
		final int type2 = pair.node2.type;
		int i;
		
		if (type1!=Constants.INTERVAL_PERIODIC && type2!=Constants.INTERVAL_PERIODIC) {
			// Assembly overlap only
			return IntervalGraphStep2.inStep2_overlap(pair);
		}
		if (type1==Constants.INTERVAL_PERIODIC && type2==Constants.INTERVAL_PERIODIC) {
			// Assembly overlap
			i=IntervalGraphStep2.isOverlap_alignment_alignment(pair);
			if (i>=0) return i;
			// Periodic overlap
			i=IntervalGraphStep2.isOverlap_simple(pair);
			if (i>=0) return Constants.overlap2periodicOverlap(i);
			else return -1;
		}
		// Assembly overlap only
		return IntervalGraphStep2.isOverlap_alignment_alignment(pair);
	}
	
	
	
	
	
	
	
	
	// ------------------------- PROCEDURES FOR EXPLORATION ONLY -------------------------
	
	/**
	 * (For exploration only)
	 * Prints the distribution of diffs values over all alignments in the graph.
	 *
	 * Remark: in the largest component of Drosophila melanogaster, the distribution 
	 * doesn't seem to have multiple peaks in a periodic component. The same happens with 
	 * $IntervalGraph.getAvgDiffsDistribution()$, and assigning min diffs rather than avg.
	 * to edges does not change it.
	 */
	private static final void getAlignmentDiffDistribution(String path, int nAlignments) throws IOException {
		if (tmpPoints==null || tmpPoints.length<nAlignments) {
			tmpPoints = new Point[nAlignments];
			for (int i=0; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
		}
		BufferedReader br = new BufferedReader(new FileReader(path));
		String str = br.readLine();
		str=br.readLine();
		int lastPoint = -1;
		while (str!=null) {
			if (str.length()==0) {
				str=br.readLine();
				continue;
			}
			Alignments.readAlignmentFile(str);
			lastPoint++;
			tmpPoints[lastPoint].position=((double)Alignments.diffs)/(Alignments.endA-Alignments.startA+1+Alignments.endB-Alignments.startB+1);
			tmpPoints[lastPoint].mass=1;
			str=br.readLine();
		}
		br.close();
		Points.sortAndCompact(tmpPoints,lastPoint);
	
		System.err.println("histogram of diffs:");
		for (int i=0; i<=lastPoint; i++) System.err.println(tmpPoints[i]);
	}
	
	
	/**
	 * (For exploration only)
	 * The procedure sets to OFF all edges with $avgDiffs>=threshold$ (where $avgDiffs$ is
	 * assumed to be an avg., not a sum). The procedure also puts all ON edges at the 
	 * beginning of each row of $IntervalGraph.neighbors$.
	 */
	private static final void turnOffWeakEdges(double threshold) {
		int i, j, from, to;
		int nEdges, nEdgesOn;
		IntervalGraph.Edge edge;
		
		nEdges=0; nEdgesOn=0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			from=i;
			if (IntervalGraph.nodesArray[from].type!=Constants.INTERVAL_PERIODIC || IntervalGraph.nodesArray[from].hasLongPeriod) {
				for (j=0; j<IntervalGraph.nNeighbors[from]; j++) IntervalGraph.neighbors[from][j].on=false;
				continue;
			}
			for (j=0; j<IntervalGraph.nNeighbors[from]; j++) {
				edge=IntervalGraph.neighbors[from][j];
				to=edge.getTo(from);
				if (IntervalGraph.nodesArray[to].type!=Constants.INTERVAL_PERIODIC || IntervalGraph.nodesArray[to].hasLongPeriod) {
					edge.on=false;
					continue;
				}
				nEdges++;
				if (edge.avgDiffs==-1 || edge.avgDiffs>=threshold) edge.on=false;
				else {
					edge.on=true;
					nEdgesOn++;
				}
			}
			IntervalGraph.Edge.order=from;
			if (IntervalGraph.nNeighbors[from]>1) Arrays.sort(IntervalGraph.neighbors[from],0,IntervalGraph.nNeighbors[from]);
		}
		System.err.println("turned off "+((nEdges-nEdgesOn)>>1)+" edges ("+((double)nEdges)/nEdgesOn+"%)");
	}
	
	
	/**
	 * (For exploration only)
	 * In the largest component of Drosophila melanogaster, the distribution doesn't have
	 * multiple peaks. Keeping only edges with large number of alignments does not
	 * disconnect the graph.
	 */
	private static final void getNAlignmentsDistribution(int nNodes, int nEdges) {
		int i, j;
		int lastN;
		
		if (IntervalGraph.tmpPoints==null || IntervalGraph.tmpPoints.length<nEdges) {
			IntervalGraph.tmpPoints = new Point[nEdges];
			for (i=0; i<IntervalGraph.tmpPoints.length; i++) IntervalGraph.tmpPoints[i] = new Point();
		}
		lastN=-1;
		for (i=0; i<nNodes; i++) {
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
				lastN++;
				IntervalGraph.tmpPoints[lastN].position=IntervalGraph.neighbors[i][j].getNAlignments();
				IntervalGraph.tmpPoints[lastN].mass=1;
			}
		}
		lastN=Points.sortAndCompact(IntervalGraph.tmpPoints,lastN);
		System.err.println("nAlignments:");
		for (i=0; i<=lastN; i++) System.err.println(IntervalGraph.tmpPoints[i]);
	}
	
	
}