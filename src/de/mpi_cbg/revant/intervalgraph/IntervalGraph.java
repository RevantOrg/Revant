package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Map;
import java.util.Set;
import java.io.*;
import java.awt.image.*;
import java.awt.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Stream;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.DAG;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Colors;
import de.mpi_cbg.revant.factorize.Intervals;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Factorize;
import de.mpi_cbg.revant.factorize.PeriodicSubstringInterval;
import de.mpi_cbg.revant.factorize.DenseSubstring;
import de.mpi_cbg.revant.factorize.AlignmentInterval;
import de.mpi_cbg.revant.factorize.PeriodicSubstrings;
import de.mpi_cbg.revant.factorize.Reads;


public class IntervalGraph {
	/**
	 * Parameters of the pipeline
	 */
	public static final int MIN_COMPONENTS = 1;
	private static final int MIN_PERIODIC_SUBSTRINGS_PER_INTERVAL = 10;
	public static final int IDENTITY_THRESHOLD = Intervals.absoluteThreshold;
	
	/**
	 * Codes for the encoding of maximal cliques
	 */
	public static final int MAXIMAL_CLIQUE_FOUND = -1;
	public static final int MAXIMAL_CLIQUE_NOT_FOUND = -2;
	public static final int MAXIMAL_CLIQUE_UP = -3;
	
	/**
	 * Encoding of the graph
	 */
	public static Node[] nodesArray;
	public static int nNodes;
	public static Edge[][] neighbors;
	public static int[] nNeighbors;
	public static HashMap<Node,Node> nodes;  // Used by $printConnectedComponents$
	
	/**
	 * Number of nodes of type:
	 * 0: alignment
	 * 1: prefix or suffix
	 * 2: prefix and suffix
	 * 3: substring
	 * 4: single deletion
	 * 5: periodic
	 */
	public static int[] nodeStats = new int[6];
	
	/**
	 * Number of edges of type:
	 *
	 * Nonperiodic-nonperiodic:
	 * 0: containment, identical
	 * 1: containment proper
	 * 2: overlap
	 * 3: insertion (one direction only)
	 * 4: shared substring
	 *
	 * Nonperiodic-periodic:
	 * 5: containment, identical
	 * 6: containment, periodic->nonperiodic
	 * 7: containment, nonperiodic->periodic
	 * 8: overlap
	 * 9: insertion, periodic->nonperiodic
	 * 10: insertion, nonperiodic->periodic
	 * 11: shared substring
	 *
	 * Periodic-periodic:
	 * 12: containment, identical
	 * 13: containment proper
	 * 14: overlap, assembly
	 * 15: overlap, periodic
	 * 16: insertion
	 * 17: shared substring
	 *
	 * 18: periodic-nonperiodic marks
	 * 
	 * 19: implied marks
	 */
	public static int[] edgeStats = new int[20];
	
	/** 
	 * Temporary space used for seed set construction.
     */
	private static boolean[] marked;
	private static int[] boundary, lastOutNeighbor, firstSeed, nInNeighbors, clique2newClique;
	public static int[] seeds, clusters;
	private static NodeCliquePair[] pairs;
	private static int[][] outNeighbors;
	
	/**
	 * Temporary space used just by this module
	 */
	public static int read1, read2, start1, end1, start2, end2;
	public static int[] stack;  // Shared stack used by all steps of the pipeline
	public static Point[] tmpPoints;
	public static int lastTmpPoint;
	private static SeedExtensionPair pair = new SeedExtensionPair(-1,-1,-1);
	private static int[] tmpArray = new int[30];
	public static int maxAlignmentsPerRead;
	private static boolean[] tmpBoolean = new boolean[(MIN_PERIODIC_SUBSTRINGS_PER_INTERVAL)<<1];
	private static Edge tmpEdge = new Edge();
	private static Edge tmpEdge2 = new Edge();
	private static Node[] reusableNodes;
	public static int maxCliqueSize;
	
	
	/**
	 * Builds the interval graph using all the alignment of the genome. The graph connects
	 * every two intervals that are assigned to the same alignment. The procedure computes
	 * also data needed in the following steps of graph partitioning.
	 *
	 * Remark: at the end of the procedure, $nodesArray$ is sorted by $READ_START_NODEID$,
	 * the $avgDiffs$ field of every edge is non-negative, and equals the \emph{sum} of 
	 * the $diffs$ values of all the alignments assigned to the edge (use procedure 
	 * $setAvgDiffs()$ to transform $avgDiffs$ from sum to avg.). Except stated otherwise,
	 * all procedures in this module assume $avgDiffs$ to be a sum rather than an avg.
	 *
	 * Remark: some alignments might be assigned to intervals in just one read, and this 
	 * could happen even in Step 2, since an alignment could be assigned to just one 
	 * interval that belongs to a Step 2 component. Some alignments might not be assigned
	 * to an interval in both reads, but this cannot happen in Step 2.
	 *
	 * Remark: using nonprimary alignments is necessary. Consider e.g. a simple repeat in 
	 * readI that aligns with a fraction of itself in readJ, sorrounded by long random
	 * insertions. The alignment would be discarded when processing readI, so it would not
	 * be primary. However, we need to connect the substring of the repeat in readJ with
	 * the full repeat in readI.
	 *
	 * Remark: the idea of building a graph on the result of an all-against-all alignment
	 * of reads, and of running an algorithm for finding communities in such graph, 
	 * already appeared in \cite{staton2015transposome,novak2013repeat,kim2006bag}.
	 *
	 * Remark: a similar graph that connects substrings of reads was described in 
	 * \cite{li2007novel,bandyopadhyay2010repfrag}, with similar criteria to connect 
	 * intervals by alignments, but those papers use a database of known repeats. Those
	 * papers also explicitly try to reconstruct nested repeats. Our notion of maximality 
	 * of alignments is also related to their notion of "penalty zone", but they don't 
	 * seem to handle it correctly sometimes. A very similar graph to ours was also used 
	 * in \cite{agarwal1994repeat}, built de novo from self-alignments, and its connected 
	 * components have been considered as families of repeats (it was also observed that
	 * well conserved repeats should have high connectivity, approaching that of cliques).
	 * For each clique, they build a minimum spanning tree to approximate evolutionary
	 * relationships. However, intervals are built by just merging alignments with Jaccard
	 * distance at least 90%. In \cite{quesneville2003detection}, an alignment is merged 
	 * to a group if one of its intervals overlaps enough with the group; groups that 
	 * share sequence locations are then merged into clusters. In \cite{bejerano2004into}
	 * they do clustering of conserved non-coding elements in multiple genomes, and they
	 * note our same problems of multi-domain sequences and wrong edges collapsing 
	 * distinct clusters to the same connected component. They use min-weight cuts and 
	 * "local articulation points", i.e. points such that the graph induced by their 
	 * neighbors has a small min-weight cut. Such points are cloned into two points and 
	 * connected to distinct subsets of neighbors.
	 * 
	 * Remark: paper \cite{kolmogorov2018assembly} is also assigning tags to repeat 
	 * occurrences, by just computing connected components of occurrences that align 
	 * pairwise.
	 *
	 * Remark: a similar "linkage graph" in which edges are local similarities, but on 
	 * multidomain proteins, is used in \cite{matsuda1999classifying}. Some proteins might
	 * contain multiple domains and thus bridge different clusters. Rather than finding 
	 * connected components (=single-linkage clustering), they find a minimal set of 
	 * covers, such that each cover is a maximal, quasi-complete subgraph, and covers can 
	 * overlap ("quasi-complete" means that every vertex must have a min. degree inside 
	 * the subgraph). They show that finding such covers is NP-hard and they describe an
	 * approximation algorithm.
	 *
	 * @param periodicSubstringsPath,denseSubstringsPath,alignmentsPath ignored if
	 * $fixDanglingEdges=FALSE$;
	 * @param discardContainedInShortPeriod discards all intervals that are contained 
	 * inside a union of short-period periodic substring intervals in their reads, and 
	 * tries to reuse the alignments assigned to them (see $IntervalGraphStep3.
	 * getKernels()$ for motivation);
	 * @param discardContainedInLongPeriod same as above, but for intervals contained in
	 * a long-period interval;
	 * @param discardPeriodicContainedInPeriodic discards only periodic intervals (of any 
	 * type) contained in another periodic interval (of any type);
	 * @param discardDenseInDense discards dense substrings contained in dense substrings
	 * of substring type;
	 * @return max number of alignments per read.
	 */
	public static final int buildIntervalGraph(String alignmentsFilePath, String connectionFilePath, String periodicSubstringsPath, String denseSubstringsPath, String alignmentsPath, boolean discardContainedInShortPeriod, boolean discardContainedInLongPeriod, boolean discardPeriodicContainedInPeriodic, boolean discardDenseInDense, boolean fixDangling, boolean fixUnassigned, boolean addSameReadEdges, String tagsDir, String tagsPrefix, boolean fullProcedure) throws IOException {
		final int MIN_ALIGNMENT_LENGTH = Alignments.minAlignmentLength-IO.quantum;  // Arbitrary
		final int HASHMAP_CAPACITY = 8000000;  // Arbitrary
		int i;
		int readA, previousReadA, readB, firstRead, lastRead, nChanged, nEdges;
		int nodeIDGenerator, alignmentsPerRead, notAssigned, doubleInsertions, disconnected;
		String str1, str2, key, tmp;
		AlignmentPair value;
		UnassignedAlignmentPair valuePrime;
		Node tmpNode, tmpNodePrime;
		Edge tmpEdge, tmpEdgePrime;
		HashMap<String,AlignmentPair> alignmentPairs = new HashMap<String,AlignmentPair>(HASHMAP_CAPACITY,1);
		HashMap<String,UnassignedAlignmentPair> unassignedAlignmentPairs = new HashMap<String,UnassignedAlignmentPair>(8000000,1);
		nodes = new HashMap<Node,Node>();
		HashMap<Edge,Edge> edges = new HashMap<Edge,Edge>();
		Iterator<Edge> edgeIterator;
		BufferedReader alignmentsFile, connectionFile;
		Node queryNode = new Node();
		Edge queryEdge = new Edge();
		
		// Loading hashmaps $nodes$ and $edges$
		alignmentsFile = new BufferedReader(new FileReader(alignmentsFilePath),IO.BUFFER_SIZE);
		connectionFile = new BufferedReader(new FileReader(connectionFilePath),IO.BUFFER_SIZE);
		nodeIDGenerator=-1; maxAlignmentsPerRead=0; previousReadA=-1; alignmentsPerRead=0;
		notAssigned=0; firstRead=-1; lastRead=-1;
		str1=alignmentsFile.readLine();
		str1=alignmentsFile.readLine();  // Skipping the first two lines
		str1=alignmentsFile.readLine();
		str2=connectionFile.readLine();
		i=0;
		while (str1!=null) {
			i++;		
			if (str2.length()==0) {
				// Updating $unassignedAlignmentPairs$.
				Alignments.readAlignmentFile(str1);
				readA=Alignments.readA-1;  // Read IDs in LAshow start from one, but in all our code they consistently start from zero.
				readB=Alignments.readB-1;
				key=IO.getCanonicalForm(readA,readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
				valuePrime=unassignedAlignmentPairs.get(key);
				if (valuePrime==null) {
					if (readA<readB) valuePrime = new UnassignedAlignmentPair(readA,readB,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,Alignments.diffs);
					else valuePrime = new UnassignedAlignmentPair(readB,readA,Alignments.startB,Alignments.endB,Alignments.startA,Alignments.endA,Alignments.orientation,Alignments.diffs);
					unassignedAlignmentPairs.put(key,valuePrime);
				}
				// Next iteration
				notAssigned++;
				str1=alignmentsFile.readLine();
				str2=connectionFile.readLine();
				continue;
			}
			Alignments.readAlignmentFile(str1);
			readA=Alignments.readA-1;  // Read IDs in LAshow start from one, but in all our code they consistently start from zero.
			if (firstRead==-1) firstRead=readA;
			lastRead=readA;
			readB=Alignments.readB-1;
			if (readA!=previousReadA) {
				if (alignmentsPerRead>maxAlignmentsPerRead) maxAlignmentsPerRead=alignmentsPerRead;
				alignmentsPerRead=1;
				previousReadA=readA;
			}
			else alignmentsPerRead++;
			Factorize.readIntervalConnectionFile(str2);
			if ( (discardContainedInShortPeriod && (Factorize.isContained&(1<<Constants.INTERVAL_PERIODIC))!=0) ||
				 (discardContainedInLongPeriod && (Factorize.isContained&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0) ||
				 ( discardPeriodicContainedInPeriodic && Factorize.type==Constants.INTERVAL_PERIODIC && 
				   ( (Factorize.isContained&(1<<Constants.INTERVAL_PERIODIC))!=0 || 
					 (Factorize.isContained&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0
				   ) 
				 ) ||
				 ( discardDenseInDense && Factorize.type>=Constants.INTERVAL_DENSE_PREFIX && Factorize.type<=Constants.INTERVAL_DENSE_SINGLEDELETION && 
				   (Factorize.isContained&(1<<Constants.INTERVAL_DENSE_SUBSTRING))!=0
				 )
			   ) {
				// An alignment that is assigned to a filtered-out interval is treated
				// like an unassigned alignment.
				key=IO.getCanonicalForm(readA,readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
				valuePrime=unassignedAlignmentPairs.get(key);
				if (valuePrime==null) {
					if (readA<readB) valuePrime = new UnassignedAlignmentPair(readA,readB,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,Alignments.diffs);
					else valuePrime = new UnassignedAlignmentPair(readB,readA,Alignments.startB,Alignments.endB,Alignments.startA,Alignments.endA,Alignments.orientation,Alignments.diffs);
					unassignedAlignmentPairs.put(key,valuePrime);
				}
				// Next iteration
				notAssigned++;
				str1=alignmentsFile.readLine();
				str2=connectionFile.readLine();
				continue;
			}
			
			// Getting the node that corresponds to the interval
			queryNode.set(readA,Factorize.type,Factorize.id,Factorize.start,Factorize.end,Factorize.isLeftMaximal,Factorize.isRightMaximal,Factorize.isWeak,Factorize.isContained,Factorize.minPrefixLength,Factorize.minSuffixLength,Factorize.period,Factorize.hasLongPeriod,Factorize.firstMaximalStart,Factorize.lastMaximalEnd,Factorize.nAssignedAlignments,Factorize.avgCoverage);
			tmpNodePrime=nodes.get(queryNode);
			if (tmpNodePrime==null) {
				tmpNode = new Node(readA,Factorize.type,Factorize.id,Factorize.start,Factorize.end,Factorize.isLeftMaximal,Factorize.isRightMaximal,Factorize.isWeak,Factorize.isContained,Factorize.minPrefixLength,Factorize.minSuffixLength,Factorize.period,Factorize.hasLongPeriod,Factorize.firstMaximalStart,Factorize.lastMaximalEnd,Factorize.nAssignedAlignments,Factorize.avgCoverage);
				nodeIDGenerator++;
				tmpNode.nodeID=nodeIDGenerator;
				nodes.put(tmpNode,tmpNode);
			}
			else tmpNode=tmpNodePrime;

			// Using the alignment to build an edge
			key=IO.getCanonicalForm(readA,readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
			value=alignmentPairs.get(key);
			if (value==null) {
				if (readA<readB) value = new AlignmentPair(tmpNode,null,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,Alignments.diffs,Factorize.implied==1);
				else value = new AlignmentPair(null,tmpNode,Alignments.startB,Alignments.endB,Alignments.startA,Alignments.endA,Alignments.orientation,Alignments.diffs,Factorize.implied==1);
				alignmentPairs.put(key,value);
			}
			else {
				alignmentPairs.remove(key);
				if (readA<readB) value.node1=tmpNode;
				else value.node2=tmpNode;
				if ( ( Intervals.isApproximatelyContained(value.alignmentStart1,value.alignmentEnd1,value.node1.start,value.node1.end) ||  
					   Intervals.areApproximatelyIdentical(value.alignmentStart1,value.alignmentEnd1,value.node1.start,value.node1.end)
					 ) &&
				     ( Intervals.isApproximatelyContained(value.alignmentStart2,value.alignmentEnd2,value.node2.start,value.node2.end) ||  
					   Intervals.areApproximatelyIdentical(value.alignmentStart2,value.alignmentEnd2,value.node2.start,value.node2.end)
					 )
				   ) {
					// Alignments that do not satisfy the conditions above can be assigned
					// the wrong Step2 type after trim().
					value.trim();
					if (value.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && value.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
						queryEdge.set(value);
						tmpEdgePrime=edges.get(queryEdge);
						if (tmpEdgePrime==null) {
							tmpEdge = new Edge(value);
							edges.put(tmpEdge,tmpEdge);
						}
						else {
							tmpEdgePrime.addTypes(queryEdge);
							tmpEdgePrime.addDiffs(queryEdge);
							tmpEdgePrime.addOrientation(queryEdge);
						}
					}
				}
			}
			
			// Updating $unassignedAlignmentPairs$.
			unassignedAlignmentPairs.remove(key);
			
			// Displaying progress
			if (i%100000==0) System.err.println("Processed alignment "+(i/100000)+"*10^5: nEdges="+edges.size()+" nNodes="+nodes.size());
			str1=alignmentsFile.readLine();
			str2=connectionFile.readLine();
		}
		if (alignmentsPerRead>maxAlignmentsPerRead) maxAlignmentsPerRead=alignmentsPerRead;
		alignmentsFile.close(); connectionFile.close();
		if (alignmentPairs.size()>0) System.err.println(alignmentPairs.size()+" alignments assigned to intervals in just one read ("+IO.getPercent(alignmentPairs.size(),i/2)+"% of all distinct alignments)");
		if (unassignedAlignmentPairs.size()>0) System.err.println(unassignedAlignmentPairs.size()+" alignments not assigned to intervals in either read ("+IO.getPercent(unassignedAlignmentPairs.size(),i/2)+"% of all distinct alignments)");
		alignmentsFile=null; connectionFile=null;
		
		// Additional edges: dangling/unassigned edges, same-read edges.
		if (fixDangling && alignmentPairs.size()>0) {
			fixDanglingEdges(alignmentPairs,nodes,edges,periodicSubstringsPath,denseSubstringsPath,alignmentsPath,maxAlignmentsPerRead,discardContainedInShortPeriod,discardContainedInLongPeriod,discardPeriodicContainedInPeriodic,discardDenseInDense);
			alignmentPairs.clear(); 
		}
		alignmentPairs=null;		
		if (fixUnassigned && unassignedAlignmentPairs.size()>0) {
			fixUnassignedAlignments(unassignedAlignmentPairs,nodes,edges,periodicSubstringsPath,denseSubstringsPath,alignmentsPath,maxAlignmentsPerRead,discardContainedInShortPeriod,discardContainedInLongPeriod,discardPeriodicContainedInPeriodic,discardDenseInDense);
			unassignedAlignmentPairs.clear(); 
		}
		unassignedAlignmentPairs=null;		
		if (addSameReadEdges) addSameReadEdges(firstRead,lastRead,periodicSubstringsPath,denseSubstringsPath,alignmentsPath,nodes,edges,maxAlignmentsPerRead);

		// Converting $nodes$ into $nodesArray$ and $edges$ into $neighbors$ and
		// $nNeighbors$.
		stack = new int[nodes.size()];
		buildNodesArray(nodes);
		nodes.clear(); nodes=null;
		buildNeighborsMatrix(edges);
		nEdges=edges.size();
		System.out.println(nNodes+" nodes and "+nEdges+" edges built ("+IO.format(((double)nEdges)/nNodes)+" edges per node on avg.)");
		edges.clear(); edges=null;
		doubleInsertions=unmarkInsertionEdges();
		if (doubleInsertions>0) System.err.println(doubleInsertions+" edges with two insertions");
		if (!fullProcedure) {
			// Remark: at this point $nodesArray$ is not necessarily sorted.
			nChanged=IntervalGraphStep2.containment2insertion(false,stack);
			System.out.println(nChanged+" edges changed from containment to insertion ("+IO.format(((double)(100*nChanged))/nEdges)+"%)");
			nChanged=IntervalGraphStep2.overlap2sharedSubstring(false,stack);
			System.out.println(nChanged+" edges changed from overlap to shared substring ("+IO.format(((double)(100*nChanged))/nEdges)+"%)");
			return maxAlignmentsPerRead;
		}		
		
		//
		// WARNING: when adding new steps below, make sure that $nodesArray$ is kept
		// sorted whenever needed.
		//
		
		// ---------------------- FIXES RELATED TO PERIODIC INTERVALS --------------------
		sortNodesArray();
		makeAlignmentIntervalsPeriodic();  // Propagating periodic tag
		fixPeriodicLoops();  // Transforming short loops into periodic
		fixMultiEdges();  // Transforming overlap multi-edges into periodic
		sortNodesArray();
		setPeriodicSubstringsOfNodes_sameRead(false);
		mergePeriodicSubstringsOfNodes();
		fixPeriodicUndersplits(true,tmpArray);  // Fixing periodic undersplits
		removePeriodicNonperiodicEdges();
		
		// -------------------- FIXES RELATED TO WEAK DENSE SUBSTRINGS -------------------
		setDenseSubstringsOfNodes_sameRead(false);		
		mergePeriodicSubstringsOfNodes();		
		fixPeriodicUndersplits(false,tmpArray);		
		redirectWeakSubstrings(true);
		redirectWeakSubstrings(false);
		
		// Fixing undersplits
/*		setPeriodicSubstringsOfNodes_neighbors(false,tmpArray);
		sortNodesArray();
		setPeriodicSubstringsOfNodes_sameRead(true,true);
		mergePeriodicSubstringsOfNodes();		
		// This is commented out since we don't currently run $fixUndersplits()$ multiple 
		// times:
		//boolean fixed;
		//do { fixed=fixUndersplits(); } while (fixed);
		fixUndersplits(true);
*/
		// Remark: resetting the $periodicSubstrings$ list of nodes after fixing
		// undersplits is not necessary, since the lists of existing nodes are still
		// correct, and new nodes are either periodic or have no periodic substring.
/*		sortNodesArray();
		fixOverlapForks();
		sortNodesArray();
		fixOverlapForks_suffixPrefix(true);
*/		
		// Removing short overlaps
/*		removeShortOverlaps();
*/				
				
		// ----------------------------- DISCONNECTING NODES -----------------------------
		disconnectLowCoverageNodes();		
		i=nNodes;
		disconnected=removeDisconnectedNodes(tagsDir,tagsPrefix);
		System.err.println(disconnected+" disconnected nodes removed ("+IO.getPercent(disconnected,i)+"%)");
		
		// Finalizing properties
		finalizeEdgeAndNodeProperties();
		
		// Adding supplement edges
/*		sortNodesArray();
		addIdentityEdges(alignmentsFilePath);
		addEdgesToOverlapNeighbors(alignmentsFilePath);
*/		
		
		// Transforming containment edges into insertion edges
		nChanged=IntervalGraphStep2.containment2insertion(false,stack);
		System.out.println(nChanged+" edges changed from containment to insertion ("+IO.format(((double)(100*nChanged))/nEdges)+"%)");
		nChanged=IntervalGraphStep2.overlap2sharedSubstring(false,stack);
		System.out.println(nChanged+" edges changed from overlap to shared substring ("+IO.format(((double)(100*nChanged))/nEdges)+"%)");
		
		// Ensuring that the final result is consistent
		if (IO.CONSISTENCY_CHECKS) {
			System.err.println("Consistency checks");
			checkConsistency(0,nNodes-1,new int[maxAlignmentsPerRead<<1]);
			System.err.println("Consistency checks passed");
		}		
		
		// Remark: at this point $nodesArray$ is sorted.
		return maxAlignmentsPerRead;
	}
	
	
	/**
	 * By $READ_START_NODEID$. Updates $neighbors$, $nNeighbors$, and pointers from 
	 * edges.
	 */
	public static final void sortNodesArray() {
		int i, j;
		int from, to, tmpNNeighbors;
		boolean[] fixed;
		int[] old2new;
		Edge[] tmpEdges;
		
		// Sorting $nodesArray$
		Node.order=Node.READ_START_NODEID;
		if (nNodes>1) Arrays.sort(nodesArray,0,nNodes);
		
		// Permuting $neighbors$ and $nNeighbors$
		fixed = new boolean[nNodes];
		for (i=0; i<nNodes; i++) {
			if (fixed[i] || nodesArray[i].nodeID==i) continue;
			tmpEdges=neighbors[i]; tmpNNeighbors=nNeighbors[i];
			to=i; from=nodesArray[to].nodeID;
			while (from!=i) {
				neighbors[to]=neighbors[from];
				nNeighbors[to]=nNeighbors[from];
				fixed[to]=true;
				to=from; from=nodesArray[to].nodeID;
			}
			neighbors[to]=tmpEdges;
			nNeighbors[to]=tmpNNeighbors;
			fixed[to]=true;
		}
		fixed=null;
		
		// Fixing pointers to nodes from edges
		old2new = new int[nNodes];
		for (i=0; i<nNodes; i++) old2new[nodesArray[i].nodeID]=i;
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				neighbors[i][j].printed=-1;  // Using $printed$ to avoid processing the same edge twice
			}
		}
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].printed!=-1) continue;
				neighbors[i][j].nodeID1=old2new[neighbors[i][j].nodeID1];
				neighbors[i][j].nodeID2=old2new[neighbors[i][j].nodeID2];
				neighbors[i][j].printed=1;
			}
		}
		old2new=null;
		
		// Reassigning $nodeID$ fields to nodes
		for (i=0; i<nNodes; i++) nodesArray[i].nodeID=i;
	}
	
	
	/**
	 * Ensures that $nodesArray$, $nNeighbors$ and $neighbors$ are large enough to 
	 * contain position $position$.
	 */
	private static final void expandTo(int position) {
		final int INCREMENT = Math.max(1000,position-nodesArray.length+1);
		int i;
		int[] tmpNNeighbors;
		Node[] tmpNodesArray;
		Edge[][] tmpNeighbors;
		
		if (position<nodesArray.length) return;
		tmpNodesArray = new Node[nodesArray.length+INCREMENT];
		System.arraycopy(nodesArray,0,tmpNodesArray,0,nodesArray.length);
		nodesArray=tmpNodesArray;
		tmpNNeighbors = new int[nNeighbors.length+INCREMENT];
		System.arraycopy(nNeighbors,0,tmpNNeighbors,0,nNeighbors.length);
		nNeighbors=tmpNNeighbors;
		tmpNeighbors = new Edge[neighbors.length+INCREMENT][0];
		for (i=0; i<neighbors.length; i++) {
			if (neighbors[i]==null) continue;
			tmpNeighbors[i] = new Edge[neighbors[i].length];
			System.arraycopy(neighbors[i],0,tmpNeighbors[i],0,neighbors[i].length);
		}
		neighbors=tmpNeighbors;
	}
	
	
	/**
	 * Unmarks edges with two types of insertion, and transforms them into identity or 
	 * shared substring if they have no other type. The criterion for identity is copied
	 * from $IntervalGraphStep2.inStep2_containment()$.
	 *
	 * @return the total number of edges with two types of insertion.
	 */
	private static final int unmarkInsertionEdges() {
		boolean longPeriod1, longPeriod2;
		int i, j;
		int type1, type2, out;
		Edge tmpEdge;
		
		out=0;
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				tmpEdge=neighbors[i][j];
				if (tmpEdge.insertion!=Constants.INSERTION_BOTH) continue;
				out++;
				tmpEdge.insertion=-1;
				if (tmpEdge.nDistinctTypes()==0) {
					type1=nodesArray[tmpEdge.nodeID1].type;
					type2=nodesArray[tmpEdge.nodeID2].type;
					longPeriod1=nodesArray[tmpEdge.nodeID1].hasLongPeriod;
					longPeriod2=nodesArray[tmpEdge.nodeID2].hasLongPeriod;
					if ( (type1!=Constants.INTERVAL_PERIODIC && type2!=Constants.INTERVAL_PERIODIC) ||
						 (type1==Constants.INTERVAL_PERIODIC && type2==Constants.INTERVAL_PERIODIC && longPeriod1==longPeriod2)
					   ) tmpEdge.containment=Constants.CONTAINMENT_IDENTICAL;
					else tmpEdge.sharedSubstring=Constants.SHARED_SUBSTRING;
				}
			}
		}
		return out;
	}
	
	
	public static final void buildNodesArray(HashMap<Node,Node> nodes) {
		Node tmpNode;
		Iterator<Node> nodeIterator;
		
		nodesArray = new Node[nodes.size()];
		nodeIterator=nodes.keySet().iterator();
		while (nodeIterator.hasNext()) {
			tmpNode=nodeIterator.next();
			nodesArray[tmpNode.nodeID]=tmpNode;
		}
		nNodes=nodesArray.length;
	}
	
	
	private static final void buildNeighborsMatrix(HashMap<Edge,Edge> edges) {
		int i, nodeID1, nodeID2;
		Edge tmpEdge;
		Iterator<Edge> edgeIterator;
		
		if (nNeighbors==null || nNeighbors.length<nNodes) nNeighbors = new int[nNodes];
		Math.set(nNeighbors,nNodes-1,0);
		edgeIterator=edges.keySet().iterator();
		while (edgeIterator.hasNext()) {
			tmpEdge=edgeIterator.next();
			nNeighbors[tmpEdge.nodeID1]++;
			nNeighbors[tmpEdge.nodeID2]++;
		}
		if (neighbors==null || neighbors.length<nNodes) neighbors = new Edge[nNodes][0];
		for (i=0; i<nNodes; i++) {
			if (nNeighbors[i]>neighbors[i].length) neighbors[i] = new Edge[nNeighbors[i]];
		}
		Math.set(nNeighbors,nNodes-1,0);
		edgeIterator=edges.keySet().iterator();
		while (edgeIterator.hasNext()) {
			tmpEdge=edgeIterator.next();
			nodeID1=tmpEdge.nodeID1; nodeID2=tmpEdge.nodeID2;
			neighbors[nodeID1][nNeighbors[nodeID1]++]=tmpEdge;
			neighbors[nodeID2][nNeighbors[nodeID2]++]=tmpEdge;
		}
	}
	
	
	/**
	 * Adds to the interval graph a shared substring edge with $avgDiffs=0$ (and one 
	 * artificial alignment, to avoid divisions by zero) for every pair of intervals in 
	 * the same read that overlap or are contained. This is useful, since the 
	 * factorization can produce a number of intervals around approximately the same
	 * region of a read.
	 *
	 * Remark: the fact that two intervals share an exact substring does not imply that
	 * their corresponding kernels are more similar than with alignments, since intervals 
	 * are just occurrences of kernels. So there is no reason for treating intra-read 
	 * edges in a different way than other shared substring edges.
	 *
	 * Remark: this procedure makes periodic intervals in the graph connected also by 
	 * shared substring edges with $avgDiffs=0$, in addition to $CONTAINMENT_IDENTICAL$
	 * edges.
	 *
	 * Remark: the procedure does not create new nodes.
	 *
	 * @param firstRead,lastRead zero-based, since read IDs in factorization files start 
	 * from zero.
	 */
	private static final void addSameReadEdges(int firstRead, int lastRead, String periodicSubstringsFile, String denseSubstringsFile, String alignmentsFile, HashMap<Node,Node> nodes, HashMap<Edge,Edge> edges, int maxIntervalsPerRead) throws IOException {
		final int READ_AHEAD_LIMIT = maxIntervalsPerRead*100;  // Characters
		final int THRESHOLD = Alignments.minAlignmentLength>>1;
		final int TYPE_COLUMN = 1;
		final int[] ID_COLUMN = new int[] {0,0,0};
		final int[] START_COLUMN = new int[] {1,2,1};
		final int[] END_COLUMN = new int[] {2,3,2};
		int i, j, k, nAdded;
		int[] lastInterval = new int[3];
		int[][] periodicIntervals = new int[maxIntervalsPerRead][11];
		int[][] denseIntervals = new int[maxIntervalsPerRead][14];
		int[][] alignmentIntervals = new int[maxIntervalsPerRead][8];
		int[][][] intervals = new int[3][0][0];
		BufferedReader periodicBuffer = new BufferedReader(new FileReader(periodicSubstringsFile),IO.BUFFER_SIZE);
		periodicBuffer.mark(READ_AHEAD_LIMIT);
		BufferedReader denseBuffer = new BufferedReader(new FileReader(denseSubstringsFile),IO.BUFFER_SIZE);
		denseBuffer.mark(READ_AHEAD_LIMIT);
		BufferedReader alignmentsBuffer = new BufferedReader(new FileReader(alignmentsFile),IO.BUFFER_SIZE);
		alignmentsBuffer.mark(READ_AHEAD_LIMIT);
		
		nAdded=0;
		for (i=firstRead; i<=lastRead; i++) { 			
			lastInterval[0]=PeriodicSubstringInterval.loadPeriodicIntervals(i,periodicBuffer,READ_AHEAD_LIMIT,periodicIntervals,tmpArray);
			lastInterval[1]=DenseSubstring.loadDenseIntervals(i,denseBuffer,READ_AHEAD_LIMIT,denseIntervals,tmpArray);
			lastInterval[2]=AlignmentInterval.loadAlignmentIntervals(i,alignmentsBuffer,READ_AHEAD_LIMIT,alignmentIntervals,tmpArray);
			intervals[0]=periodicIntervals;
			intervals[1]=denseIntervals;
			intervals[2]=alignmentIntervals;
			for (j=0; j<3; j++) {
				if (lastInterval[j]<0) continue;
				nAdded+=addSameReadEdges1(i,j,intervals[j],lastInterval[j],START_COLUMN[j],END_COLUMN[j],ID_COLUMN[j],TYPE_COLUMN,THRESHOLD,nodes,edges);
				for (k=j+1; k<3; k++) {
					if (lastInterval[k]<0) continue;
					nAdded+=addSameReadEdges2(i,j,intervals[j],lastInterval[j],START_COLUMN[j],END_COLUMN[j],ID_COLUMN[j],TYPE_COLUMN,k,intervals[k],lastInterval[k],START_COLUMN[k],END_COLUMN[k],ID_COLUMN[k],TYPE_COLUMN,THRESHOLD,nodes,edges);
				}
			}
		}
		periodicBuffer.close(); denseBuffer.close(); alignmentsBuffer.close();
		System.err.println("Added "+nAdded+" same-read edges");
	}
	
	
	/**
	 * Adds a shared substring edge for each pair of distinct overlapping or contained 
	 * intervals in $intervals[0..lastInterval]$. No new node is added to $nodes$.
	 *
	 * Remark: the procedure assumes $intervals$ to be sorted by first position.
	 *
	 * Remark: read IDs start from zero in factorization files, but start from one in
	 * the interval graph, since they are inherited from LAshow.
	 *
	 * @param read zero-based;
	 * @param type type of $intervals$: 0: periodic; 1: dense; 2: alignment;
	 * @return the number of edges created by the procedure.
	 */
	private static final int addSameReadEdges1(int read, int type, int[][] intervals, int lastInterval, int startColumn, int endColumn, int idColumn, int typeColumn, int threshold, HashMap<Node,Node> nodes, HashMap<Edge,Edge> edges) {
		final int edgeType;
		int i, j;
		int out, startI, startJ, endI, min, max;
		Node tmpNodeI = new Node(read,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Node tmpNodeJ = new Node(read,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Node tmpNodeIPrime, tmpNodeJPrime;
		Edge tmpEdge, tmpEdgePrime;
		
		if (type==0) tmpNodeI.type=Constants.INTERVAL_PERIODIC;
		else if (type==2) tmpNodeI.type=Constants.INTERVAL_ALIGNMENT;
		edgeType=Constants.SHARED_SUBSTRING;
		out=0;
		for (i=0; i<lastInterval; i++) {
			startI=intervals[i][startColumn];
			endI=intervals[i][endColumn];
			if (type==1) tmpNodeI.type=intervals[i][typeColumn];
			tmpNodeI.id=intervals[i][idColumn];
			tmpNodeIPrime=nodes.get(tmpNodeI);
			if (tmpNodeIPrime==null) continue;
			for (j=i+1; j<=lastInterval; j++) {
				startJ=intervals[j][startColumn];
				if (startJ>endI-threshold+1) break;
				if (Intervals.intersectionLength(startI,endI,startJ,intervals[j][endColumn])<threshold) continue;
				if (type==1) tmpNodeJ.type=intervals[j][typeColumn];
				tmpNodeJ.id=intervals[j][idColumn];
				tmpNodeJPrime=nodes.get(tmpNodeJ);
				if (tmpNodeJPrime==null) continue;
				if (tmpNodeIPrime.nodeID<tmpNodeJPrime.nodeID) {
					min=tmpNodeIPrime.nodeID;
					max=tmpNodeJPrime.nodeID;
				}
				else {
					min=tmpNodeJPrime.nodeID;
					max=tmpNodeIPrime.nodeID;
				}
				tmpEdge = new Edge(min,max,type==0&&tmpNodeIPrime.hasLongPeriod!=tmpNodeJPrime.hasLongPeriod?Constants.SHARED_SUBSTRING:edgeType);
				if (edges.get(tmpEdge)!=null) {
					System.err.println("ERROR: same-read edge already in $edges$?! (1) read="+read);
					System.exit(1);
				}
				edges.put(tmpEdge,tmpEdge);  // The edge is absent from $edges$, since it connects intervals in the same read.				
				out++;
			}
		}
		return out;
	}
	
	
	/**
	 * Like $addSameReadEdges1$, but for two sets of intervals.
	 *
	 * @param read zero-based.
	 */
	private static final int addSameReadEdges2( int read, int type1, int[][] intervals1, int lastInterval1, int startColumn1, int endColumn1, int idColumn1, int typeColumn1, 
	 												      int type2, int[][] intervals2, int lastInterval2, int startColumn2, int endColumn2, int idColumn2, int typeColumn2, int threshold, HashMap<Node,Node> nodes, HashMap<Edge,Edge> edges ) {
		final int edgeType;
		boolean firstJForNextIFound;
		int i, j;
		int out, start1, start2, end1, end2, firstJForNextI, min, max;
		Node tmpNode1 = new Node(read,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Node tmpNode2 = new Node(read,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Node tmpNode1Prime, tmpNode2Prime;
		Edge tmpEdge, tmpEdgePrime;
		
		if (type1==0) {
			tmpNode1.type=Constants.INTERVAL_PERIODIC;
			if (type2==0) {  
				// The procedure is never called with this combination
				tmpNode2.type=Constants.INTERVAL_PERIODIC;
			}
			else if (type2==2) tmpNode2.type=Constants.INTERVAL_ALIGNMENT;
		}
		else {
			if (type1==2) tmpNode1.type=Constants.INTERVAL_ALIGNMENT;
			if (type2==0) tmpNode2.type=Constants.INTERVAL_PERIODIC;
			else if (type2==2) tmpNode2.type=Constants.INTERVAL_ALIGNMENT;
		}
		edgeType=Constants.SHARED_SUBSTRING;
		out=0; firstJForNextI=0;
		for (i=0; i<=lastInterval1; i++) {
			start1=intervals1[i][startColumn1];
			end1=intervals1[i][endColumn1];
			if (type1==1) tmpNode1.type=intervals1[i][typeColumn1];
			tmpNode1.id=intervals1[i][idColumn1];
			tmpNode1Prime=nodes.get(tmpNode1);
			if (tmpNode1Prime==null) continue;
			firstJForNextIFound=false;
			for (j=firstJForNextI; j<=lastInterval2; j++) {
				end2=intervals2[j][endColumn2];
				if (!firstJForNextIFound && i<lastInterval1 && end2>=intervals1[i+1][startColumn1]) firstJForNextI=j;
				if (end2<start1+threshold-1) continue;
				start2=intervals2[j][startColumn2];
				if (start2>end1-threshold+1) break;
				if (Intervals.intersectionLength(start1,end1,start2,end2)<threshold) continue;
				if (type2==1) tmpNode2.type=intervals2[j][typeColumn2];
				tmpNode2.id=intervals2[j][idColumn2];
				tmpNode2Prime=nodes.get(tmpNode2);
				if (tmpNode2Prime==null) continue;
				if (tmpNode1Prime.nodeID<tmpNode2Prime.nodeID) {
					min=tmpNode1Prime.nodeID;
					max=tmpNode2Prime.nodeID;
				}
				else {
					min=tmpNode2Prime.nodeID;
					max=tmpNode1Prime.nodeID;
				}
				tmpEdge = new Edge(min,max,edgeType);
				if (edges.get(tmpEdge)!=null) {
					System.err.println("ERROR: same-read edge already in $edges$?! (2) read="+read);
					System.err.println("one node:");
					System.err.println(tmpNode1Prime);
					System.err.println("other node:");
					System.err.println(tmpNode2Prime);
					System.exit(1);
				}
				edges.put(tmpEdge,tmpEdge);  // The edge is absent from $edges$, since it connects intervals in the same read.
				out++;		
			}
		}
		return out;
	}
	
	
	/**
	 * Sets the $components[0]$ field of every node in $nodesArray$, using only edges with
	 * $on=true$ (if $onlyOn=true$), or all edges (if $onlyOn=false$), in matrix
	 * $neighbors$.
	 *
	 * @param sorted if TRUE, assumes that all edges with $on=true$ in $neighbors[i]$ are
	 * stored in a prefix of $neighbors[i]$, for every $i$;
	 * @return the number of components of size at least $minSize$.
	 */
	public static final int getConnectedComponents(int minSize, boolean sorted, boolean onlyOn) {
		int i, j;
		int top, from, lastComponent, to, size, nSmallComponents;
		
		if (stack==null || stack.length<nNodes) stack = new int[nNodes];
		lastComponent=-1; nSmallComponents=0;
		for (i=0; i<nNodes; i++) {
			nodesArray[i].lastComponent=-1;
			if (nodesArray[i].components==null) nodesArray[i].components = new int[MIN_COMPONENTS];
		}
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].lastComponent!=-1) continue;
			nodesArray[i].lastComponent=0;
			nodesArray[i].components[0]=++lastComponent;
			size=1;
			if (nNeighbors[i]==0) {
				if (minSize>1) nSmallComponents++;
				continue;
			}
			stack[0]=i;
			top=0;
			while (top>=0) {
				from=stack[top--];
				for (j=0; j<nNeighbors[from]; j++) {
					if (onlyOn && !neighbors[from][j].on) {
						if (sorted) break;
						else continue;
					}
					to=neighbors[from][j].getTo(from);
					if (nodesArray[to].lastComponent!=-1) continue;
					nodesArray[to].lastComponent=0;
					nodesArray[to].components[0]=nodesArray[from].components[0];
					size++;
					stack[++top]=to;
				}
			}
			if (size<minSize) nSmallComponents++;
		}
		return lastComponent+1-(minSize>1?nSmallComponents:0);
	}
	
	
	/**
	 * @return the total number of edges.
	 */
	public static final int getNEdges() {
		int i, nEdges;
		
		nEdges=0;
		for (i=0; i<nNodes; i++) nEdges+=nNeighbors[i];
		return nEdges>>1;
	}
	
	
	/**
	 * @return the number of ON edges.
	 */
	public static final int getNOnEdges(boolean sorted) {
		int i, j;
		int nEdges;
		
		nEdges=0;
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (!neighbors[i][j].on) {
					if (sorted) break;
					else continue;
				}
				nEdges++;
			}
		}
		return nEdges>>1;
	}
	
	
	/**
	 * @return the number of ON edges between nodes in $subgraph[first..last]$.
	 */
	public static final int getNOnEdges(int[] subgraph, int first, int last, boolean sorted) {
		int i, j, from, nEdges;
		
		nEdges=0;
		for (i=first; i<=last; i++) {
			from=subgraph[i];
			for (j=0; j<nNeighbors[from]; j++) {
				if (!neighbors[from][j].on) {
					if (sorted) break;
					else continue;
				}
				nEdges++;
			}
		}
		return nEdges>>1;
	}
	
	
	/**
	 * Sets to OFF all supplement edges, and all normal edges that are not of type $type$,
	 * and sets to ON all other edges. The procedure puts all edges with $on=true$ at the 
	 * beginning of $neighbors[i]$, for every $i$, then sorts egdes by $getTo(i)$.
	 *
	 * @param type 0=containment, 1=overlap, 2=insertion, 3=sharedSubstring, -1=turn on
	 * all normal edges of containment, overlap, insertion, or sharedSubstring type; 
	 * -2=turn on all normal edges;
	 * @param preserveContainment TRUE=does not change the current ON state of normal 
	 * edges of containment type;
	 * @return the max degree of a node, considering only ON edges.
	 */
	public static final int turnOffEdges(int type, boolean preserveContainment) {
		int i, j;
		int decision, originalOrder, out, onNeighbors;
		Edge edge;
		
		originalOrder=Edge.order;
		out=0;
		for (i=0; i<nNodes; i++) {
			onNeighbors=0;
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (edge.supplement) {
					edge.on=false;
					continue;
				}
				if (type==0) {
					if (edge.containment==-1) edge.on=false;
					else {
						if (!preserveContainment) edge.on=true;
					}
				}
				else {
					decision=-1;
					switch (type) {
						case 1: decision=edge.overlap; break;
						case 2: decision=edge.insertion; break;
						case 3: decision=edge.sharedSubstring; break;
						case -1: decision=(edge.containment!=-1||edge.overlap!=-1||edge.insertion!=-1||edge.sharedSubstring!=-1)?Math.POSITIVE_INFINITY:-1; break;
						case -2: decision=Math.POSITIVE_INFINITY; break;
					}
					edge.on=decision!=-1;
				}
				if (edge.on) onNeighbors++;
			}
			Edge.order=Edge.UNSORTED+1+i;
			if (nNeighbors[i]>1) Arrays.sort(neighbors[i],0,nNeighbors[i]);
			if (onNeighbors>out) out=onNeighbors;
		}
		Edge.order=originalOrder;
		System.err.println("nEdges in graph: "+getNEdges()+" ON="+getNOnEdges(true));
		return out;
	}
	
	
	/**
	 * Variant of $turnOffEdges$ that works only on the edges between nodes in
	 * $subgraph[first..last]$, which is assumed to be sorted. Edges between a node in the
	 * subgraph and a node outside the subgraph are set to OFF, regardless of their type.
	 */
	public static final void turnOffEdges(int[] subgraph, int first, int last, int type) {
		int i, j;
		int from, decision, originalOrder;
		Edge edge;
		
		originalOrder=Edge.order;
		for (i=first; i<=last; i++) {
			from=subgraph[i];
			for (j=0; j<nNeighbors[from]; j++) {
				edge=neighbors[from][j];
				if (Arrays.binarySearch(subgraph,first,last+1,edge.getTo(from))<0) {
					edge.on=false;
					continue;
				}
				if (edge.supplement) {
					edge.on=false;
					continue;
				}
				decision=-1;
				switch (type) {
					case 0: decision=edge.containment; break;
					case 1: decision=edge.overlap; break;
					case 2: decision=edge.insertion; break;
					case 3: decision=edge.sharedSubstring; break;
				}
				edge.on=decision!=-1;
			}
			Edge.order=Edge.UNSORTED+1+from;
			if (nNeighbors[from]>1) Arrays.sort(neighbors[from],0,nNeighbors[from]);
		}
		Edge.order=originalOrder;
	}
	
	
	/**
	 * Variant of $turnOffEdges$ that turns on all edges (including supplement).
	 */
	public static final void turnOnEdges() {
		int i, j;
		
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].on=true;
		}
	}
	
	
	/**
	 * Puts all edges with $on=true$ at the beginning of $neighbors[i]$, for every $i$ 
	 * contained in $subgraph[first..last]$, then sorts egdes by $getTo(i)$.
	 */
	public static final void sortEdges(int[] subgraph, int first, int last) {
		int i, node;
		
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			Edge.order=Edge.UNSORTED+1+node;
			if (nNeighbors[node]>1) Arrays.sort(neighbors[node],0,nNeighbors[node]);
		}
	}

	
	
	

	
	
	
	// --------------------------- RECURSIVE PEELING PROCEDURES --------------------------
	
	/**
	 * Let $subgraph[first..last]$ contain all nodes in a subgraph, which is assumed to be 
	 * connected to the rest of the graph just by edges with $on=false$. Let the 
	 * "frontier" of the subgraph be the set of all its maximal/minimal nodes by 
	 * containment, as well as the set of unary paths that start from a node that is 
	 * adjacent to just maximal/minimal nodes and to the next node in the unary path 
	 * (excluding the last node in the unary path), as well as all nodes reachable from 
	 * such nodes by $CONTAINMENT_IDENTICAL$ edges. The procedure turns off all edges in 
	 * the subgraph that are incident to at least one node in the frontier, and marks 
	 * frontier nodes.
	 *
	 * Remark: the procedure sets the maximal/minimal tag of each node of the subgraph.
	 *
	 * Remark: some periodic nodes might belong to the frontier, but the construction of 
	 * the frontier does not start from periodic nodes.
	 *
	 * @param frontier the procedure fills $frontier[first..X]$ with the IDs, in the 
	 * original graph, of all the frontier nodes of the subgraph, where $X=out[0]$.
	 * $out[1]$ (respectively, $out[2]$) contains the number of maximal (respectively, 
	 * minimal) nodes in the frontier.
	 */
	public static final void turnOffFrontier(int[] subgraph, int first, int last, int[] frontier, int[] out) {
		boolean isContained, contains;
		int i, j, k, f;
		int to, lastFrontier, node1, node2, node3, top, nMaximal, nMinimal;
		int inNeighbor, inNeighborPrime, outNeighbor, outNeighborPrime;
		Edge edge;
		
		// Counting maximal/minimal nodes
		if (stack==null || stack.length<nNodes) stack = new int[nNodes];
		nMaximal=0; nMinimal=0;
		for (i=first; i<=last; i++) {
			node1=subgraph[i];
			nodesArray[node1].isMaximal=false; nodesArray[node1].isMinimal=false;
			isContained=false; contains=false; 
			for (j=0; j<nNeighbors[node1]; j++) {
				edge=neighbors[node1][j];
				if (!edge.on) continue;
				node2=edge.getTo(node1);
				if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==node1) ||
					 (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID2==node1)
				   ) isContained=true;
				if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID2==node1) ||
					 (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID1==node1)
				   ) contains=true;
			}
			if (isContained && contains) continue;
			if (!isContained && contains) {
				nodesArray[node1].isMaximal=true;
				nMaximal++;
			}
			if (isContained && !contains) {
				nodesArray[node1].isMinimal=true;
				nMinimal++;
			}
		}
		
		// Marking minimal/maximal nodes and unary paths
		for (i=first; i<=last; i++) {
			nodesArray[subgraph[i]].visited=-1;
			nodesArray[subgraph[i]].frontier=false;
		}
		lastFrontier=first-1;
		for (i=first; i<=last; i++) {
			node1=subgraph[i];
			if (nodesArray[node1].type==Constants.INTERVAL_PERIODIC || nodesArray[node1].frontier) continue;
			isContained=false; contains=false; 
			for (j=0; j<nNeighbors[node1]; j++) {
				edge=neighbors[node1][j];
				if (!edge.on) continue;
				node2=edge.getTo(node1);
				if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==node1) ||
					 (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID2==node1)
				   ) isContained=true;
				if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID2==node1) ||
					 (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID1==node1)
				   ) contains=true;
			}
			if (isContained && contains) continue;
			frontier[++lastFrontier]=node1;
			nodesArray[node1].frontier=true;
			
			// Turning off strictly-contained edges incident to a maximal/minimal vertex,
			// and pushing the destination nodes on the stack.
			nodesArray[node1].visited=node1;
			top=-1;
			for (j=0; j<nNeighbors[node1]; j++) {
				edge=neighbors[node1][j];
				if (!edge.on) continue;
				if (edge.containment!=Constants.CONTAINMENT_ONE_IN_TWO && edge.containment!=Constants.CONTAINMENT_TWO_IN_ONE) continue;
				edge.on=false;
				node2=edge.getTo(node1);
				if (nodesArray[node2].visited!=node1 && !nodesArray[node2].frontier) {
					nodesArray[node2].visited=node1;
					stack[++top]=node2;
				}
			}
			
			// Turning off all egdes in unary paths from maximal/minimal nodes
			while (top>=0) {
				node2=stack[top--];
				inNeighbor=-1; outNeighbor=-1; inNeighborPrime=-1; outNeighborPrime=-1;
				for (j=0; j<nNeighbors[node2]; j++) {
					edge=neighbors[node2][j];
					if (!edge.on) continue;
					node3=edge.getTo(node2);
					if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==node2) ||
						 (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID2==node2)
					   ) {
						if (outNeighbor==-1) {
							outNeighbor=node3;
							outNeighborPrime=j;
						}
						else if (outNeighbor>=0 && outNeighbor!=node3) {
							outNeighbor=-2;
							outNeighborPrime=-2;
						}
					}
					else if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID2==node2) ||
						      (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID1==node2)
					        ) {
   						if (inNeighbor==-1) {
							inNeighbor=node3;
							inNeighborPrime=j;
						}
   						else if (inNeighbor>=0 && inNeighbor!=node3) {
							inNeighbor=-2;
							inNeighborPrime=-2;
						}
					}
				}
				if ( (inNeighbor>=0 && outNeighbor==-1) || 
				     (inNeighbor==-1 && outNeighbor>=0) ||
					 (inNeighbor==-1 && outNeighbor==-1)
				   ) {
					frontier[++lastFrontier]=node2;
					nodesArray[node2].frontier=true;
					if (inNeighbor>=0) {
						neighbors[node2][inNeighborPrime].on=false;
						if (nodesArray[inNeighbor].visited!=node1) {
							nodesArray[inNeighbor].visited=node1;
							stack[++top]=inNeighbor;
						}
					}
					if (outNeighbor>=0) {
						neighbors[node2][outNeighborPrime].on=false;
						if (nodesArray[outNeighbor].visited!=node1) {
							nodesArray[outNeighbor].visited=node1;
							stack[++top]=outNeighbor;
						}
					}
				}
			}
		}
		
		// Marking and disconnecting nodes reachable from the frontier by
		// $CONTAINMENT_IDENTICAL$ edges
		for (i=first; i<=last; i++) nodesArray[subgraph[i]].visited=-1;
		f=lastFrontier;
		for (i=first; i<=f; i++) {
			node1=frontier[i];
			nodesArray[node1].visited=node1;
			top=0;
			stack[top]=node1;
			while (top>=0) {
				node2=stack[top--];
				for (j=0; j<nNeighbors[node2]; j++) {
					edge=neighbors[node2][j];
					if (!edge.on) continue;
					if (edge.containment!=Constants.CONTAINMENT_IDENTICAL) continue;
					edge.on=false;
					to=edge.getTo(node2);
					if (nodesArray[to].visited==node1 || nodesArray[to].frontier) continue;
					frontier[++lastFrontier]=to;
					nodesArray[to].frontier=true;
					nodesArray[to].visited=node1;
					for (k=0; k<nNeighbors[to]; k++) {
						edge=neighbors[to][k];
						if (!edge.on) continue;
						if (edge.containment!=Constants.CONTAINMENT_ONE_IN_TWO && edge.containment!=Constants.CONTAINMENT_TWO_IN_ONE) continue;
						edge.on=false;
					}
					stack[++top]=to;
				}
			}
		}

		// Consistency checks
		if (IO.CONSISTENCY_CHECKS) {
			if (lastFrontier<first) {
				boolean onePeriodicNode = false;
				for (i=first; i<=last; i++) {
					if (nodesArray[subgraph[i]].type==Constants.INTERVAL_PERIODIC) {
						onePeriodicNode=true;
						break;
					}
				}
				if (!onePeriodicNode) {
					System.err.println("ERROR in turnOffFrontier(): empty frontier, but no periodic node.");
					System.exit(1);
				}
			}
			if (lastFrontier-first>last-first) {
				System.err.println("ERROR in turnOffFrontier(): frontier size: "+(lastFrontier-first+1)+" subgraph size: "+(last-first+1));
				System.exit(1);
			}
		}
	
		out[0]=lastFrontier; out[1]=nMaximal; out[2]=nMinimal;
	}
	
	
	/**
	 * Recursive version of $getConnectedComponents()$. Sets the $components[0]$ field of 
	 * all nodes in the subgraph induced by the nodes in $subgraph[first..last]$, using 
	 * only edges with $on=true$. The procedure assumes that the subgraph is connected to 
	 * the rest of the graph just by edges with $on=false$. Frontier nodes of the subgraph
	 * are assigned no connected component. 
	 *
	 * Remark: the procedure assumes that the nodes in the frontier of the subgraph have 
	 * already been marked.
	 *
	 * @param sorted if TRUE, assumes that all edges with $on=true$ in $neighbors[i]$ are
	 * stored in a prefix of $neighbors[i]$, for every $i$;
	 * @param firstComponent component IDs are assigned starting with $firstComponent$,
	 * inclusive;
	 * @param minSize >=1;
	 * @param out 0: the total number of connected components found in non-frontier nodes;
	 * 1: the number of components of size at least $minSize$ in non-frontier nodes.
	 */
	public static final void getConnectedComponents(int[] subgraph, int first, int last, int firstComponent, int minSize, boolean sorted, int[] out) {
		int i, j;
		int top, from, lastComponent, to, size, nSmallComponents;
		Node node;
		
		if (stack==null || stack.length<nNodes) stack = new int[nNodes];
		lastComponent=-1; nSmallComponents=0;
		for (i=first; i<=last; i++) {
			node=nodesArray[subgraph[i]];
			node.lastComponent=-1;
			if (node.components==null) node.components = new int[MIN_COMPONENTS];
		}
		for (i=first; i<=last; i++) {
			node=nodesArray[subgraph[i]];
			if (node.frontier || node.lastComponent!=-1) continue;
			node.lastComponent=0;
			lastComponent++;
			node.components[0]=firstComponent+lastComponent;			
			size=1;
			if (nNeighbors[subgraph[i]]==0) {
				if (size<minSize) nSmallComponents++;
				continue;
			}
			stack[0]=subgraph[i];
			top=0;
			while (top>=0) {
				from=stack[top--];
				for (j=0; j<nNeighbors[from]; j++) {
					if (!neighbors[from][j].on) {
						if (sorted) break;
						else continue;
					}
					to=neighbors[from][j].getTo(from);
					if (nodesArray[to].frontier || nodesArray[to].lastComponent!=-1) continue;
					nodesArray[to].lastComponent=0;
					nodesArray[to].components[0]=nodesArray[from].components[0];
					size++;
					stack[++top]=to;
				}
			}
			if (size<minSize) nSmallComponents++;
		}
		
		// Consistency checks
		if (IO.CONSISTENCY_CHECKS) {	
			if (lastComponent>last-first) {
				System.err.println("ERROR in getConnectedComponents(sorted="+sorted+"): nComponents="+(lastComponent+1)+" nNodesInSubgraph="+(last-first+1));
				System.exit(1);
			}
			for (i=first; i<=last; i++) {
				from=subgraph[i];
				if (nodesArray[from].frontier) continue;
				if (nodesArray[from].lastComponent!=0) {
					System.err.println("ERROR in getConnectedComponents(sorted="+sorted+"): non-frontier node "+from+" has "+(nodesArray[from].lastComponent+1)+" components");
					System.exit(1);
				}
				for (j=0; j<nNeighbors[from]; j++) {
					if (!neighbors[from][j].on) continue;
					if (neighbors[from][j].containment==-1) continue;
					to=neighbors[from][j].getTo(from);
					boolean found = false;
					for (int k=first; k<=last; k++) {
						if (subgraph[k]==to) {
							found=true;
							break;
						}
					}
					if (!found) {
						System.err.println("ERROR in getConnectedComponents(sorted="+sorted+"): node "+to+" not in the subgraph, but reachable from node "+from);
						System.exit(1);
					}
					if (nodesArray[to].lastComponent>0 || (nodesArray[to].lastComponent==-1 && !nodesArray[to].frontier)) {
						System.err.println("ERROR in getConnectedComponents(sorted="+sorted+"): node "+to+" has "+(nodesArray[to].lastComponent+1)+" components");
						System.exit(1);
					}
					if (!nodesArray[to].frontier && (nodesArray[to].components[0]!=nodesArray[from].components[0])) {
						System.err.println("ERROR in getConnectedComponents(sorted="+sorted+"): nodes "+from+" and "+to+" are connected but belong to different components: "+nodesArray[from].components[0]+"::"+nodesArray[to].components[0]);
						System.err.println("neighbors of "+from+" (lastComponent="+nodesArray[from].lastComponent+"):");
						for (int x=0; x<nNeighbors[from]; x++) System.err.print("("+neighbors[from][x].getTo(from)+","+neighbors[from][x].on+"),");
						System.err.println();
						System.err.println("neighbors of "+to+" (lastComponent="+nodesArray[to].lastComponent+"):");
						for (int x=0; x<nNeighbors[to]; x++) System.err.print("("+neighbors[to][x].getTo(to)+","+neighbors[to][x].on+"),");
						System.err.println();
						System.exit(1);
					}
				}
			}
		}
		
		out[0]=lastComponent+1;
		out[1]=lastComponent+1-nSmallComponents;
	}
	
	
	/**
	 * Assume that the subgraph induced by nodes in $subgraph[first..last]$: (1) is 
	 * connected to the rest of the graph just by edges with $on=false$; (2) is connected, 
	 * by following edges with $on=true$. Assume also that $turnOffFrontier()$ and 
	 * $getConnectedComponents()$ have been executed on the subgraph, in this order.
	 * The procedure propagates component tags from nodes not in the frontier of the
	 * subgraph, to nodes in the frontier of the subgraph. A node in the frontier can get 
	 * more than one tag. No new tag is created by the procedure. Then, the procedure
	 * turns on non-supplement edges that connect, to their component, frontier nodes with
	 * exactly one tag, and removes frontier tags from frontier nodes.
	 *
	 * Remark: the procedure uses every type of non-supplement containment edge for
	 * propagating tags.
	 *
	 * Remark: since the subgraph including the frontier is assumed to be connected, a
	 * simple propagation from non-frontier to frontier nodes suffices.
	 *
	 * @param tmp temporary space with at least M elements, where M is the maximum
	 * number of components per node;
	 * @param frontierStats the procedure adds to this matrix, for every node type (rows),
	 * the number of maximal and minimal nodes, in the frontier of the current subgraph, 
	 * that have more than one component; if $frontierStats$ is null the procedure ignores
	 * it; columns:
	 * 0: maximal nodes, total;
	 * 1: maximal nodes, leftmax. or rightmax. intervals only;
	 * 2: maximal nodes, leftmax. and rightmax. intervals only;
	 * 3: maximal nodes, neither left- nor right-maximal intervals;
	 * 4: minimal nodes, total;
	 * 5: minimal nodes, leftmax. or rightmax. intervals only;
	 * 6: minimal nodes, leftmax. and rightmax. intervals only;
	 * 7: minimal nodes, neither left- nor right-maximal intervals.
	 * @param subgraphComponent the component ID of $subgraph$.
	 */
	public static final void tagFrontier(int[] subgraph, int first, int last, int[] frontier, int lastFrontier, int[] tmp, int[][] frontierStats, int subgraphComponent, boolean logFrontierRatio) {
		boolean isContained, contains;
		int i, j;
		int from, to, lastComponent, previousLastComponent, top, type;
		Node node, nodeTo;
		Edge edge;
		int[] tmpPrime;
		
		if (logFrontierRatio) {
			for (i=first; i<=lastFrontier; i++) {
				nodesArray[frontier[i]].lastComponent=Math.POSITIVE_INFINITY;  // To signal more than one component to the following steps.
			}
		}
		else {
			// Propagating tags
			if (stack==null || stack.length<nNodes) stack = new int[nNodes];
			for (i=first; i<=last; i++) nodesArray[subgraph[i]].visited=-1;		
			for (i=first; i<=last; i++) {
				node=nodesArray[subgraph[i]];
				if (node.frontier) continue;
				node.visited=subgraph[i];
				top=0;
				stack[0]=subgraph[i];
				while (top>=0) {
					from=stack[top--];
					for (j=0; j<nNeighbors[from]; j++) {
						edge=neighbors[from][j];
						if (edge.supplement) continue;
						type=edge.containment;
						if (type==-1) continue;
						to=edge.getTo(from);
						nodeTo=nodesArray[to];
						if (!nodeTo.frontier || nodeTo.visited==subgraph[i]) continue;
						previousLastComponent=nodeTo.lastComponent;
						lastComponent=Math.setUnion(node.components,node.lastComponent,nodeTo.components,nodeTo.lastComponent,tmp);
						if (lastComponent==previousLastComponent) continue;
						if (nodeTo.components.length<lastComponent+1) nodeTo.components=Arrays.copyOf(tmp,lastComponent+1);
					    else System.arraycopy(tmp,0,nodeTo.components,0,lastComponent+1);
						nodeTo.lastComponent=lastComponent;
						nodeTo.visited=subgraph[i];
						stack[++top]=to;
					}
				}
			}
			// Consistency checks
			if (IO.CONSISTENCY_CHECKS) {
				for (i=first; i<=lastFrontier; i++) {
					if (nodesArray[frontier[i]].lastComponent<0) {
						System.err.println("ERROR in tagFrontier(): frontier node "+frontier[i]+" has no component");
						System.exit(1);
					}
				}
			}
		
			// Connecting back to their component all frontier nodes with just one tag.
			for (i=first; i<=last; i++) nodesArray[subgraph[i]].visited=-1;
			for (i=first; i<=last; i++) {
				node=nodesArray[subgraph[i]];
				if (node.frontier) continue;
				node.visited=subgraph[i];
				top=0;
				stack[0]=subgraph[i];
				while (top>=0) {
					from=stack[top--];
					for (j=0; j<nNeighbors[from]; j++) {
						edge=neighbors[from][j];
						if (edge.supplement) continue;
						type=edge.containment;
						if (type==-1) continue;
						to=edge.getTo(from);
						nodeTo=nodesArray[to];
						if (!nodeTo.frontier || nodeTo.lastComponent!=0 || edge.on) continue;
						edge.on=true;
						if (nodeTo.visited!=subgraph[i]) {
							nodeTo.visited=subgraph[i];
							stack[++top]=to;
						}
					}
				}
			}
			// Consistency checks
			if (IO.CONSISTENCY_CHECKS) {
				for (i=first; i<=lastFrontier; i++) {
					if (nodesArray[frontier[i]].lastComponent>0) {
						for (j=0; j<nNeighbors[frontier[i]]; j++) {
							if (neighbors[frontier[i]][j].on) {
								System.err.println("ERROR in tagFrontier(): frontier node "+frontier[i]+" with more than on component has adjacent ON edges.");
								System.exit(1);
							}
						}
					}
				}
				for (i=first; i<=last; i++) {
					if (nodesArray[subgraph[i]].lastComponent>0) {
						for (j=0; j<nNeighbors[subgraph[i]]; j++) {
							if (neighbors[subgraph[i]][j].on) {
								System.err.println("ERROR in tagFrontier(): node "+subgraph[i]+" with more than on component has adjacent ON edges.");
								System.exit(1);
							}
						}
					}
				}
			}
		}

		// Unmarking frontier nodes
		for (i=first; i<=lastFrontier; i++) nodesArray[frontier[i]].frontier=false;
		
		// Computing statistics on maximal/minimal frontier nodes with multiple components
		if (frontierStats!=null) {	
			for (i=first; i<=lastFrontier; i++) {
				from=frontier[i];
				node=nodesArray[from];
				if (node.lastComponent==0) continue;
				if (node.isMaximal) {
					frontierStats[node.type][0]++;
					if (node.isLeftMaximal) {
						if (node.isRightMaximal) frontierStats[node.type][2]++;
						else frontierStats[node.type][1]++;
					}
					else if (node.isRightMaximal) frontierStats[node.type][1]++;
					else frontierStats[node.type][3]++;					
				}
				if (node.isMinimal) {
					frontierStats[node.type][4]++;
					if (node.isLeftMaximal) {
						if (node.isRightMaximal) frontierStats[node.type][6]++;
						else frontierStats[node.type][5]++;
					}
					else if (node.isRightMaximal) frontierStats[node.type][5]++;
					else frontierStats[node.type][7]++;
				}
			}
		}
	}
	
	
	/**
	 * Assume that all nodes in $subgraph[first..last]$ have been assigned one or multiple
	 * components, and that the distinct component IDs assigned to the nodes of the 
	 * subgraph form a compact interval that starts with $firstComponent$.
	 * The procedure:
	 *
	 * - discards all nodes with more than one component, and all nodes whose component 
	 * contains less than $minComponentSize$ nodes, and stores them in no particular order
	 * in $subgraph[X+1..last]$;
	 *
	 * - sorts the remaining nodes by component, then by node ID, and stores the sorted 
	 * list in $subgraph[first..X]$;
	 *
	 * - stores in $componentSize[i][1]$ the absolute starting position, in the sorted
	 * $subgraph$, of the interval corresponding to the $i$-th large component, for $i \in 
	 * [firstComponent..firstComponent+Y]$, where $Y$ is returned in output;
	 *
	 * - stores in $componentSize[i][0]$ the number of nodes in large component $i$, for 
	 * $i \in [firstComponent..firstComponent+Y]$; thus, $X=first+S-1$, where $S$ is the 
	 * sum of all values in column 0 of $componentSize$;
	 *
	 * - keeps in $componentSize$ just large components, storing the ID of each in column
	 * 2 of $componentSize$.
	 *
	 * Remark: column 3 of $componentSize$ is used just as temporary space.
	 *
	 * @param nComponents total number of components;
	 * @param tmp temporary space of size equal to $subgraph$; only $tmp[first..last]$ is 
	 * used by the procedure;
	 * @return -1 if no component contains at least $minComponentSize$ nodes.
	 */
	public static final int sortNodesByComponent(int[] subgraph, int first, int last, int firstComponent, int nComponents, int minComponentSize, int[][] componentSize, int[] tmp) {
		int i;
		int found, pointer, component, lastComponent, nNodes;
		Node node;

		// Computing the number of nodes in each component
		nNodes=last-first+1;
		for (i=0; i<nComponents; i++) {
			componentSize[firstComponent+i][0]=0;
			componentSize[firstComponent+i][1]=-1;
			componentSize[firstComponent+i][2]=firstComponent+i;
			componentSize[firstComponent+i][3]=-1;
		}
		for (i=first; i<=last; i++) {
			node=nodesArray[subgraph[i]];
			if (node.lastComponent!=0) continue;
			componentSize[node.components[0]][0]++;
		}
		
		// Filling $tmp[first..last]$ with nodes sorted by component.
		found=-1;
		for (i=0; i<nComponents; i++) {
			if (componentSize[firstComponent+i][0]>=minComponentSize) {
				found=i;
				break;
			}
		}
		if (found==-1) return -1;
		componentSize[firstComponent+found][1]=first;
		componentSize[firstComponent+found][3]=first;
		pointer=componentSize[firstComponent+found][0];
		for (i=found+1; i<nComponents; i++) {
			if (componentSize[firstComponent+i][0]<minComponentSize) continue;
			componentSize[firstComponent+i][1]=first+pointer; 
			componentSize[firstComponent+i][3]=componentSize[firstComponent+i][1];
			pointer+=componentSize[firstComponent+i][0];
		}
		for (i=first; i<=last; i++) {
			node=nodesArray[subgraph[i]];
			if (node.lastComponent!=0) {
				tmp[first+(pointer++)]=subgraph[i];
				continue;
			}
			component=node.components[0];
			if (componentSize[component][0]<minComponentSize) {
				tmp[first+(pointer++)]=subgraph[i];
				continue;
			}
			tmp[componentSize[component][3]++]=subgraph[i];
		}
		if (first+pointer!=last+1) {
			System.err.println("ERROR in sortNodesByComponent()");
			System.exit(1);
		}
		System.arraycopy(tmp,first,subgraph,first,nNodes);
		
		// Removing small components from the current $componentSize$ interval
		lastComponent=-1;
		for (i=0; i<nComponents; i++) {
			if (componentSize[firstComponent+i][0]<minComponentSize) continue;
			lastComponent++;
			if (i!=lastComponent) System.arraycopy(componentSize[firstComponent+i],0,componentSize[firstComponent+lastComponent],0,4);
		}
		
		// Sorting nodes inside each large component
		for (i=0; i<=lastComponent; i++) {
			if (componentSize[firstComponent+i][0]>1) Arrays.sort(subgraph,componentSize[firstComponent+i][1],componentSize[firstComponent+i][1]+componentSize[firstComponent+i][0]);
		}
		return lastComponent;
	}
	
	
	
	
	
	
	
	
	// -------------------------- SEED SET EXPANSION PROCEDURES --------------------------
	
	/**
	 * Sets $degrees[i]$ for every $i$ contained in $subgraph[first..last]$ (considering 
	 * only ON edges), and stores in $tmpPoints[0..x]$ the sorted and compacted list of 
	 * distinct degrees, where $x$ is returned in output.
	 *
	 * Remark: the procedure uses the global variables $tmpPoints$ and $lastTmpPoint$.
	 *
	 * @param degrees output array, of size at least equal to the number of nodes in the
	 * interval graph;
	 * @param sorted TRUE iff a prefix of $neighbors[i]$ contains all ON edges, for all 
	 * $i$.
	 */
	public static final int getDegreeDistribution(int[] subgraph, int first, int last, boolean noMinimalMaximal, int[] degrees, boolean sorted) {
		int i, j;
		int from, nNodes;
		Node node;
		Point point;
		
		nNodes=last-first+1;
		lastTmpPoint=-1;
		for (i=first; i<=last; i++) {
			from=subgraph[i]; node=nodesArray[from]; 
			degrees[from]=getDegree(from,sorted,Math.NEGATIVE_INFINITY,null);
			if ( noMinimalMaximal && ((!node.contains && node.isContained) || (node.contains && !node.isContained)) ) continue;
			lastTmpPoint++;
			ensureTmpPoints(lastTmpPoint);
			tmpPoints[lastTmpPoint].position=degrees[from];
			tmpPoints[lastTmpPoint].mass=1;
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		return lastTmpPoint;
	}
	
	
	public static final void ensureTmpPoints(int lastTmpPoint) {
		final int GROWTH_RATE = 100;  // Arbitrary
		int i;
		
		if (tmpPoints==null) {
			tmpPoints = new Point[lastTmpPoint+GROWTH_RATE];
			for (i=0; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
		}
		else if (tmpPoints.length==lastTmpPoint) {
			Point[] newTmpPoints = new Point[tmpPoints.length+GROWTH_RATE];
			System.arraycopy(tmpPoints,0,newTmpPoints,0,tmpPoints.length);
			for (i=tmpPoints.length; i<newTmpPoints.length; i++) newTmpPoints[i] = new Point();
			tmpPoints=newTmpPoints;
		}
	}
	
	
	/**
	 * @param clusterField used iff $excludeClusters!=Math.NEGATIVE_INFINITY$;
	 * @return the number of ON edges incident to $node$, excluding those that connect to
	 * nodes whose cluster is $<=excludeClusters$.
	 */
	private static final int getDegree(int nodeID, boolean sorted, int excludeClusters, int[] clusterField) {
		int i, degree;
		Edge edge;
		
		degree=0;
		for (i=0; i<nNeighbors[nodeID]; i++) {
			edge=neighbors[nodeID][i];
			if (!edge.on) {
				if (sorted) break;
				else continue;
			}
			if (excludeClusters!=Math.NEGATIVE_INFINITY && clusterField[edge.getTo(nodeID)]<=excludeClusters) continue;
			degree++;
		}
		return degree;
	}
	
	
	/**
	 * @return the maximum degree of a node.
	 */
	public static final int getMaxDegree() {
		int i, n, max;
		
		max=0;
		for (i=0; i<nNodes; i++) {		
			n=nNeighbors[i];
			if (n>max) max=n;
		}
		return max;
	}
	
	
	/**
	 * Let $subgraph[first..last]$ be a connected subgraph, connected to the rest of the 
	 * interval graph just by OFF edges. The procedure grows cluster seeds, inside 
	 * $subgraph[first..last]$, by computing all maximal cliques, and by keeping only
	 * large cliques with small intersection. The procedure stores in the global
	 * variable $seeds[0..X]$ all the distinct nodes in the clusters, sorted by node ID,
	 * and in $clusters[0..X]$ the cluster ID of each node in $seeds$, where $X=out[0]$.
	 *
	 * Remark: we use seed set expansion to detect clusters, rather than simply removing
	 * all edges with alignment score lower than a threshold, since this simpler method
	 * does not work in the non-periodic parts of the interval graph.
	 *
	 * Remark: rather than growing clusters from seeds, we could just start with one 
	 * cluster per vertex, and iteratively move the vertex with largest modularity gain to
	 * the cluster of one of its neighbors \cite{blondel2008fast}. Alternatively, we could
	 * consider at each step the merge of two clusters that produces the largest gain in
	 * modularity \cite{newman2004fast,clauset2004finding,novak2010graph}, or we could use 
	 * the Louvain method that consists in assigning a cluster to each node, then move a
	 * node from its cluster to the one of its neighbors that maximizes the modularity 
	 * gain, then collapse all nodes in the same cluster into a single new node, and 
	 * repeat \cite{guo2017replong}. We haven't tried these approaches yet, they might be 
	 * faster. Paper \cite{estill2010taxonomy} mentions other clustering methods that seem
	 * to work well for graphs that look very similar to ours: affinity propagation 
	 * \cite{frey2007clustering}, and Markov clustering \cite{van2000graph}; they also 
	 * seem very fast. We haven't tried these approaches either. Maximal cliques have been
	 * used to find clusters in an overlap graph also in \cite{baaijens2017novo} (using 
	 * also transitive edges but not "double-transitive" edges).
	 *
	 * Remark: after the procedure completes, the subgraphs induced by all nodes in the 
     * same cluster are connected.
	 *
	 * Remark: the procedure disregards OFF edges.
	 *
	 * Remark: the procedure internally allocates a number of matrices whose size depends
	 * on the number of large maximal cliques in the subgraph. Since an upper bound to the
	 * size of such matrices is hard to compute by the caller, they are not passed as
	 * input parameters.
	 *
	 * @param minCliqueSize (>=3) minimum size of a clique to be considered;
	 * @param stream temporary space;
	 * @param maxSize (ignored if -1) rather than considering all maximal cliques as 
	 * seeds, the procedure considers all cliques of exactly this size (not necessarily 
	 * maximal);
	 * @param samplingRate if positive, uses just a random subset, of this relative size,
	 * of the nodes at the root of the enumeration tree (not building all maximal cliques
	 * might yield fewer seeds, so each of our clusters might contain more than one true 
	 * cluster);
	 * @param maxNCliques clique construction stops after it finds at least this number of
	 * maximal cliques;
	 * @param out $out[1]$ contains the number of distinct clusters in output.
	 */
	public static final void getSeeds(int[] subgraph, int first, int last, int minCliqueSize, Stream stream, int maxSize, double samplingRate, int maxNCliques, int[] out) {
		final double INTERSECTION_RATE = 0.1;  // Arbitrary
		boolean multipleComponents;
		int i, j;
		int lastSeed, nCliques;
		
		System.err.println("getSeeds> Computing maximal cliques...  (subgraph size="+(last-first+1)+" maxSize="+maxSize+" samplingRate="+samplingRate+")");
		getMaximalCliques(subgraph,first,last,stream,maxSize,samplingRate,maxNCliques,true);
		markMaximalCliques(subgraph,first,last,stream,minCliqueSize,INTERSECTION_RATE,out);
		lastSeed=out[0]; nCliques=out[1];
		if (nCliques==0) return;
		
		// Removing seeds that occur in multiple clusters.
		// Remark: assigning such nodes to clusters using greedy modularity maximization,
		// as done elsewhere in the clustering pipeline, does not work well in practice.
		if (pairs==null || pairs.length<lastSeed+1) pairs = new NodeCliquePair[lastSeed+1];
		for (i=0; i<lastSeed+1; i++) {
			if (pairs[i]==null) pairs[i] = new NodeCliquePair(seeds[i],clusters[i]);
			else {
				pairs[i].node=seeds[i];
				pairs[i].clique=clusters[i];
			}
		}
		NodeCliquePair.order=NodeCliquePair.NODE;
		if (pairs.length>1) Arrays.sort(pairs,0,pairs.length);
		i=1; j=-1; multipleComponents=false;
		while (i<=lastSeed) {
			if (pairs[i].node!=pairs[i-1].node) {
				if (!multipleComponents) {
					j++;
					seeds[j]=pairs[i-1].node;
					clusters[j]=pairs[i-1].clique;
				}
				multipleComponents=false;
			}
			else {
				if (pairs[i].clique!=pairs[i-1].clique) multipleComponents=true;
			}
			i++;
		}
		if (!multipleComponents) {
			j++;
			seeds[j]=pairs[lastSeed].node;
			clusters[j]=pairs[lastSeed].clique;
		}
		lastSeed=j;

		// Compacting component IDs after seed removal
		if (nInNeighbors==null || nInNeighbors.length<nCliques) nInNeighbors = new int[nCliques];  // We reuse array $nInNeighbors$.
		Math.set(nInNeighbors,nCliques-1,-1);
		for (i=0; i<=lastSeed; i++) nInNeighbors[clusters[i]]=1;
		j=0;
		for (i=0; i<nCliques; i++) {
			if (nInNeighbors[i]==1) nInNeighbors[i]=j++;
		}
		nCliques=j;
		for (i=0; i<=lastSeed; i++) clusters[i]=nInNeighbors[clusters[i]];
		out[0]=lastSeed; out[1]=nCliques;
		System.err.println("getSeeds> done.  lastSeed="+lastSeed+" nCliques="+nCliques);
	}
	
	
	private static class NodeCliquePair implements Comparable {
		public static final int UNSORTED = -1;
		public static final int CLIQUE = 0;  // Sorts by clique, then by node.
		public static final int NODE = 1;  // Sorts just by node
		public static int order;
		
		public int node;
		public int clique;
		
		public NodeCliquePair(int n, int c) {
			node=n;
			clique=c;
		}
		
		public int compareTo(Object other) {
			NodeCliquePair otherPair = (NodeCliquePair)other;
			if (order==CLIQUE) {
				if (clique>otherPair.clique) return 1;
				else if (clique<otherPair.clique) return -1;
				if (node>otherPair.node) return 1;
				else if (node<otherPair.node) return -1;
			}
			else if (order==NODE) {
				if (node>otherPair.node) return 1;
				else if (node<otherPair.node) return -1;
			}
			return 0;
		}
	}	
	
	
	/**
	 * Let $subgraph[first..last]$ be a subgraph connected to the rest of the interval 
	 * graph just by OFF edges. The procedure grows seeds into larger clusters, ultimately
	 * setting the cluster field of every node in $subgraph[first..last]$: see 
	 * procedures $extendSeeds_impl()$ for details. The procedure uses global arrays 
	 * $seeds[0..lastSeed]$ and $clusters[0..lastSeed]$ as input: $seeds$ is assumed to be
	 * sorted by nodeID, and $clusters[i]$ is assumed to contain the cluster ID of the 
	 * $i$-th seed node.
	 *
	 * @param nClusters number of distinct cluster IDs;
	 * @param clusterSize,deltas temporary space, of size at least $nClusters$;
	 * @param nEdges number of ON edges in the subgraph;
	 * @param queue temporary space, assumed not NULL; can grow;
	 * @param clusterField temporary space, of size at least equal to the total number of 
	 * nodes in the interval graph; at the end of the procedure, 
	 * $clusterField[subgraph[i]]$ contains the cluster assignment of node $subgraph[i]$;
	 * @param degrees degree of every node in the interval graph;
	 * @param sorted TRUE iff ON edges are at the beginning of each row of $neighbors$;
	 * @param maxAddedNodes the procedure stops after this many nodes have been added
	 * (ignored if non-positive);
	 */
	public static final void extendSeeds(int[] subgraph, int first, int last, int lastSeed, int nClusters, int[] clusterSize, int nEdges, double[] deltas, PriorityQueue<SeedExtensionPair> queue, int[] clusterField, int[] degrees, boolean sorted, int maxAddedNodes) {
		int i, j, k;
		int nNodes, nNodesInCluster, from, to;
		Edge edge;
		
		nNodes=last-first+1;
		for (i=first; i<=last; i++) clusterField[subgraph[i]]=-1;
		for (i=0; i<=lastSeed; i++) clusterField[seeds[i]]=clusters[i];
		Math.set(clusterSize,nClusters-1,0);
		for (i=0; i<=lastSeed; i++) clusterSize[clusters[i]]++;
		nNodesInCluster=lastSeed+1;
		queue.clear();
		for (i=0; i<=lastSeed; i++) {
			from=seeds[i];
			for (j=0; j<nNeighbors[from]; j++) {
				edge=neighbors[from][j];
				if (!edge.on) {
					if (sorted) break;
					else continue;
				}
				to=edge.getTo(from);
				pair.node=to;
				if (clusterField[to]!=-1 || queue.contains(pair)) continue;
				lastTmpPoint=getAdjacentClusters(to,clusterField,sorted);
				k=deltaModularity(to,nEdges,deltas,tmpPoints,lastTmpPoint,clusterField,degrees,sorted);
				queue.add(new SeedExtensionPair(to,k,deltas[k]));
			}
		}
		i=extendSeeds_impl(clusterSize,nEdges,queue,deltas,clusterField,degrees,sorted,maxAddedNodes);
		if (maxAddedNodes<=0 && i!=nNodes-nNodesInCluster) {
			System.err.println("ERROR in extendSeeds: the graph contains multiple connected components?");
			System.err.println("added: "+i+" not seeds: "+(nNodes-nNodesInCluster));
			System.exit(1);
		}
		System.err.println("extendSeeds> seeds extended. Added "+i+" nodes.");
	}
	
	
	/**
	 * Greedily extends clusters by iteratively adding the node, in the boundary of all 
	 * current clusters, that maximizes the difference in modularity score between the 
	 * current partition, and the new partition obtained by adding the node to one of the
	 * clusters. The boundary of all clusters is the set of nodes that do not belong to 
	 * any cluster, and that are adjacent to at least one node that belongs to a cluster.
	 * See procedure $deltaModularity$ for details. For an overview of seed set expansion
	 * methods, see e.g. \cite{xie2013overlapping}.
	 *
	 * Remark: we use simple greedy approaches throughout the seed extension pipeline 
	 * because finding a partition of maximal modularity is NP-hard, and because greedy
	 * works in practice.
	 *
	 * Remark: iteratively adding the node that minimizes the expansion or the 
	 * conductance, or that maximizes the coverage, or that maximizes the performance, of 
	 * a bipartition $(C_i,V \setminus C_i)$, is expected to create similar clusters in 
	 * practice. See e.g. \cite{kannan2004clusterings,brandes2007engineering} for a
	 * definition of objective functions.
	 *
	 * @param clusterSize size of each cluster before the procedure is called; the 
	 * procedure updates these counts;
	 * @param nEdges number of ON edges in the subgraph;
	 * @param queue contains all nodes in the boundary of a cluster (the boundary is grown
	 * by the procedure);
	 * @param deltas temporary space, of size at least equal to the number of clusters;
	 * @param clusterField temporary space, of size at least equal to the total number of 
	 * nodes in the interval graph;
	 * @param degrees degree of every node in the interval graph;
	 * @param maxAddedNodes the procedure stops after this many nodes have been added
	 * (ignored if non-positive);
	 * @return the number of nodes added to clusters by the procedure.
	 */
	private static final int extendSeeds_impl(int[] clusterSize, int nEdges, PriorityQueue<SeedExtensionPair> queue, double[] deltas, int[] clusterField, int[] degrees, boolean sorted, int maxAddedNodes) {
		int j, k, to, node, nAddedNodes;
		Edge edge;
		SeedExtensionPair tmpPair;
		
		nAddedNodes=0;
		while (!queue.isEmpty()) {		
			tmpPair=queue.poll();
			node=tmpPair.node;
			clusterField[node]=tmpPair.cluster;
			clusterSize[tmpPair.cluster]++;
			nAddedNodes++;
			if (nAddedNodes%1000==0) System.err.println("extendSeeds_impl> added "+nAddedNodes+" nodes");
			if (maxAddedNodes>0 && nAddedNodes>=maxAddedNodes) break;
			for (j=0; j<nNeighbors[node]; j++) {
				edge=neighbors[node][j];
				if (!edge.on) {
					if (sorted) break;
					else continue;
				}
				to=edge.getTo(node);
				if (clusterField[to]!=-1) continue;
				lastTmpPoint=getAdjacentClusters(to,clusterField,sorted);
				k=deltaModularity(to,nEdges,deltas,tmpPoints,lastTmpPoint,clusterField,degrees,sorted);
				pair.node=to;
				if (queue.contains(pair)) queue.remove(pair);
				queue.add(new SeedExtensionPair(to,k,deltas[k]));
			}
		}
		return nAddedNodes;
	}
	
	
	/**
	 * Refines clusters by: (1) improving their boundaries; (2) discarding clusters that 
	 * are too small after improving boundaries, and reassigning their nodes; (3) removing 
	 * all nodes in the boundaries (contraction); (4) discarding clusters that are too 
	 * small after removing boundary nodes, and reassigning their nodes. 
	 * See procedures $refineClusters_boundary$ and $refineClusters_contract$ for details.
	 *
	 * @param nClusters number of distinct cluster IDs;
	 * @param clusterSize size of each cluster before refinement; the procedure updates
	 * these counts;
	 * @param clusterSizePrime temporary space, of size at least equal to $nClusters$;
	 * @param queue temporary space; can grow;
	 * @param deltas temporary space, of size at least equal to $nClusters$;
	 * @param clusterTmp,outNeighbors,lastOutNeighbor,nInNeighbors,components temporary 
	 * space, of size at least equal to the number of nodes in the subgraph;
	 * @param clusterField temporary space, of size at least equal to the total number of 
	 * nodes in the interval graph;
	 * @param nEdges number of ON edges in the subgraph;
	 * @param degrees degree of every node in the interval graph;
	 * @param sorted TRUE iff ON edges are at the beginning of each row of $neighbors$;
	 * @return the number $C$ of clusters after refinement; if $C>1$, cluster IDs after 
	 * refinement form a compact interval $[0..C-1]$, otherwise compactness is not
	 * guaranteed, since the procedure exits as soon as it detects that there is just one 
	 * cluster.
	 */
	public static final int refineClusters(int[] subgraph, int first, int last, int nClusters, int[] clusterSize, int[] clusterSizePrime, int nEdges, PriorityQueue<SeedExtensionPair> queue, double[] deltas, int[] clusterTmp, int[] clusterTmpPrime, int[][] outNeighbors, int[] lastOutNeighbor, int[] nInNeighbors, int[] components, int[] clusterField, int[] degrees, boolean sorted) {
		final double CLUSTER_SIZE_RELATIVE_THRESHOLD = 0.05;  // Arbitrary
		final int CLUSTER_SIZE_ABSOLUTE_THRESHOLD = 100;  // Arbitrary
		final int CLUSTER_SIZE_THRESHOLD = Math.min(CLUSTER_SIZE_ABSOLUTE_THRESHOLD,(int)(CLUSTER_SIZE_RELATIVE_THRESHOLD*(last-first+1)));
		final double BOUNDARY_THRESHOLD = 0.25;  // Arbitrary
		int i, j, k;
		int node, nNodes, nNodesInCluster, cluster, nActiveClusters;
	
		// Refining cluster boundaries. This can move the boundary between adjacent
		// clusters, and make some clusters smaller.
		nNodes=last-first+1;
		j=0;
		do {
			queue.clear();
			for (i=first; i<=last; i++) {
				node=subgraph[i];
				cluster=clusterField[node];
				lastTmpPoint=getAdjacentClusters(node,clusterField,sorted);
				if (lastTmpPoint==0 && (int)tmpPoints[0].position==cluster) continue;
				k=deltaModularity(node,nEdges,deltas,tmpPoints,lastTmpPoint,clusterField,degrees,sorted);
				// We add also boundary nodes with $k==cluster$, since they might be
				// reassigned later.
				queue.add(new SeedExtensionPair(node,k,deltas[k]));
			}
			i=refineClusters_boundary(clusterSize,nEdges,queue,deltas,clusterField,degrees,sorted);
		} while (i>0);
		
		// Reassigning the nodes of clusters that become too small after refining
		// boundaries.
		queue.clear();
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			if (clusterSize[clusterField[node]]<CLUSTER_SIZE_THRESHOLD) clusterField[node]=-1;
		}
		nActiveClusters=0; nNodesInCluster=0;
		for (i=0; i<nClusters; i++) {
			if (clusterSize[i]<CLUSTER_SIZE_THRESHOLD) clusterSize[i]=0;
			else {
				nActiveClusters++;
				nNodesInCluster+=clusterSize[i];
			}
		}
		if (nActiveClusters<=1) {
			clusterSize[0]=nNodes;
			for (i=first; i<=last; i++) clusterField[subgraph[i]]=0;
			return 1;
		}
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			if (cluster!=-1) continue;
			lastTmpPoint=getAdjacentClusters(node,clusterField,sorted);
			if (lastTmpPoint==-1) continue;
			k=deltaModularity(node,nEdges,deltas,tmpPoints,lastTmpPoint,clusterField,degrees,sorted);
			queue.add(new SeedExtensionPair(node,k,deltas[k]));
		}
		i=extendSeeds_impl(clusterSize,nEdges,queue,deltas,clusterField,degrees,sorted,-1);
		if (i!=nNodes-nNodesInCluster) {
			System.err.println("ERROR in extendSeeds (first call by refineClusters): the graph contains multiple connected components?");
			System.err.println("added: "+i+" not seeds: "+(nNodes-nNodesInCluster));
			System.exit(1);
		}
		
		// Reassigning the nodes of clusters that would become too small after contraction
		System.arraycopy(clusterSize,0,clusterSizePrime,0,nClusters);
		for (i=first; i<=last; i++) clusterTmp[i-first]=clusterField[subgraph[i]];
		i=refineClusters_contract(subgraph,first,last,nClusters,clusterSizePrime,BOUNDARY_THRESHOLD,outNeighbors,lastOutNeighbor,nInNeighbors,components,clusterTmpPrime,clusterField,sorted);
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			if (clusterSizePrime[clusterTmp[i-first]]<CLUSTER_SIZE_THRESHOLD) clusterField[node]=-1;
			else clusterField[node]=clusterTmp[i-first];
		}
		nActiveClusters=0; nNodesInCluster=0;
		for (i=0; i<nClusters; i++) {
			if (clusterSizePrime[i]<CLUSTER_SIZE_THRESHOLD) clusterSize[i]=0;
			else {
				nActiveClusters++;
				nNodesInCluster+=clusterSize[i];
			}
		}
		if (nActiveClusters<=1) {
			clusterSize[0]=nNodes;
			for (i=first; i<=last; i++) clusterField[subgraph[i]]=0;
			return 1;
		}
		queue.clear();
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			if (cluster!=-1) continue;
			lastTmpPoint=getAdjacentClusters(node,clusterField,sorted);
			if (lastTmpPoint==-1) continue;
			k=deltaModularity(node,nEdges,deltas,tmpPoints,lastTmpPoint,clusterField,degrees,sorted);
			queue.add(new SeedExtensionPair(node,k,deltas[k]));
		}
		i=extendSeeds_impl(clusterSize,nEdges,queue,deltas,clusterField,degrees,sorted,-1);
		if (i!=nNodes-nNodesInCluster) {
			System.err.println("ERROR in extendSeeds (second call by refineClusters): the graph contains multiple connected components?");
			System.err.println("added: "+i+" not seeds: "+(nNodes-nNodesInCluster));
			System.exit(1);
		}
		
		// Contracting for real
		i=refineClusters_contract(subgraph,first,last,nClusters,clusterSize,BOUNDARY_THRESHOLD,outNeighbors,lastOutNeighbor,nInNeighbors,components,clusterTmpPrime,clusterField,sorted);
		
		// Compacting the IDs of the remaining clusters
		Math.set(clusterSizePrime,nClusters-1,-1);  // Reusing $clusterSizePrime$
		j=-1;
		for (i=0; i<nClusters; i++) {
			if (clusterSize[i]==0) continue;
			j++;
			clusterSizePrime[i]=j;
			clusterSize[j]=clusterSize[i];
		}
		nClusters=j+1;
		if (nClusters==1) return 1;
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			if (clusterField[node]==-1) continue;
			clusterField[node]=clusterSizePrime[clusterField[node]];
		}
		return nClusters;
	}
	
	
	/**
	 * Simplified version of $refineClusters$ for large subgraphs. The procedure just 
	 * removes small clusters and reassigns their nodes.
	 */
	public static final int refineClusters_large(int[] subgraph, int first, int last, int nClusters, int[] clusterSize, int[] clusterSizePrime, int nEdges, PriorityQueue<SeedExtensionPair> queue, double[] deltas, int[] clusterField, int[] degrees, boolean sorted) {
		final double CLUSTER_SIZE_RELATIVE = 0.005;  // Arbitrary
		final int CLUSTER_SIZE_ABSOLUTE = 100;  // Arbitrary
		final int MIN_CLUSTER_SIZE = Math.min(CLUSTER_SIZE_ABSOLUTE,(int)(CLUSTER_SIZE_RELATIVE*nNodes));
		final int nNodes = last-first+1;
		int i, j, k;
		int node, cluster, nActiveClusters, nNodesInCluster;
		
		// Reassigning the nodes of clusters that are too small
		queue.clear();
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			if (cluster==-1) continue;
			if (clusterSize[cluster]<MIN_CLUSTER_SIZE) clusterField[node]=-1;
		}
		nActiveClusters=0; nNodesInCluster=0;
		for (i=0; i<nClusters; i++) {
			if (clusterSize[i]<MIN_CLUSTER_SIZE) clusterSize[i]=0;
			else {
				nActiveClusters++;
				nNodesInCluster+=clusterSize[i];
			}
		}
		if (nActiveClusters<=1) {
			clusterSize[0]=nNodes;
			for (i=first; i<=last; i++) clusterField[subgraph[i]]=0;
			return 1;
		}
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			if (cluster!=-1) continue;
			lastTmpPoint=getAdjacentClusters(node,clusterField,sorted);
			if (lastTmpPoint==-1) continue;
			k=deltaModularity(node,nEdges,deltas,tmpPoints,lastTmpPoint,clusterField,degrees,sorted);
			queue.add(new SeedExtensionPair(node,k,deltas[k]));
		}
		i=extendSeeds_impl(clusterSize,nEdges,queue,deltas,clusterField,degrees,sorted,-1);
		if (i!=nNodes-nNodesInCluster) {
			System.err.println("ERROR in extendSeeds (called by refineClusters_large): the graph contains multiple connected components?");
			System.err.println("added: "+i+" not seeds: "+(nNodes-nNodesInCluster));
			System.exit(1);
		}

		// Compacting the IDs of the remaining clusters
		Math.set(clusterSizePrime,nClusters-1,-1);  // Reusing $clusterSizePrime$.
		j=-1;
		for (i=0; i<nClusters; i++) {
			if (clusterSize[i]==0) continue;
			j++;
			clusterSize[j]=clusterSize[i];
			clusterSizePrime[i]=j;
		}
		nClusters=j+1;
		if (nClusters==1) return 1;
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			if (clusterField[node]==-1) continue;
			clusterField[node]=clusterSizePrime[clusterField[node]];
		}
		return nClusters;
	}
	
	
	/**
	 * Greedily reassigns the boundary node, in the boundary of all current clusters, that 
	 * maximizes the difference in modularity score between the current partition, and the
	 * new partition obtained by moving the node to a different adjacent cluster. The 
	 * boundary of all clusters is the set of nodes that have at least one neighbor in a 
	 * different cluster than theirs.
	 *
	 * @param clusterSize the procedure updates these counts;
	 * @param nEdges total number of ON edges in the subgraph;
	 * @param queue all frontier nodes (the frontier is not changed by the procedure);
	 * @param deltas temporary space, of size at least equal to the number of clusters;
	 * @param clusterField temporary space, of size at least equal to the total number of 
	 * nodes in the interval graph;
	 * @param degrees degree of every node in the interval graph;
	 * @return the number of nodes reassigned by the procedure.
	 */
	private static final int refineClusters_boundary(int[] clusterSize, int nEdges, PriorityQueue<SeedExtensionPair> queue, double[] deltas, int[] clusterField, int[] degrees, boolean sorted) {
		int j, k;
		int to, node, nReassignedNodes, oldCluster, newCluster;
		Edge edge;
		SeedExtensionPair tmpPair;
		
		nReassignedNodes=0;	
		while (!queue.isEmpty()) {
			tmpPair=queue.poll();
			node=tmpPair.node;
			oldCluster=clusterField[node];
			newCluster=tmpPair.cluster;
			if (newCluster==oldCluster) continue;
			clusterField[node]=newCluster;
			nReassignedNodes++;
			clusterSize[oldCluster]--;
			clusterSize[newCluster]++;
			for (j=0; j<nNeighbors[node]; j++) {
				edge=neighbors[node][j];
				if (!edge.on) {
					if (sorted) break;
					else continue;
				}
				to=edge.getTo(node);
				pair.node=to;
				if (!queue.contains(pair)) continue;
				lastTmpPoint=getAdjacentClusters(to,clusterField,sorted);
				k=deltaModularity(to,nEdges,deltas,tmpPoints,lastTmpPoint,clusterField,degrees,sorted);
				queue.remove(pair);
				queue.add(new SeedExtensionPair(to,k,deltas[k]));
			}
		}
		return nReassignedNodes;
	}
	
	
	/**
	 * @return the ratio between the number of neighbors of $node$ in a different nonempty
	 * cluster than $node$, and the total number of neighbors of $node$ in a nonempty
	 * cluster; -1 if $node$ has no neighbor in a nonempty cluster.
	 */
	private static final double neighborsInDifferentCluster(int node, int[] clusterField, boolean sorted) {
		int j, to, cluster, clusterTo, denominator;
		double numerator;
		Edge edge;
	
		cluster=clusterField[node];
		numerator=0.0; denominator=0;
		for (j=0; j<nNeighbors[node]; j++) {
			edge=neighbors[node][j];
			if (!edge.on) {
				if (sorted) break;
				else continue;
			}
			to=edge.getTo(node);
			clusterTo=clusterField[to];
			if (clusterTo==-1) continue;
			denominator++;
			if (clusterTo!=cluster) numerator+=1.0;
		}
		return denominator==0?-1.0:numerator/denominator;
	}
	
	
	/**
	 * Removes all nodes in the boundary of clusters, and keeps just the largest connected
	 * component of each cluster after removal. "Boundary" in this procedure has a 
	 * stricter meaning than in e.g. $refineClusters_boundary$: a node is in the boundary 
	 * iff the ratio of its neighbors in a different cluster is greater than 
	 * $boundaryThreshold$. This is because nodes with large degree inside a cluster might
	 * be connected also to nodes in a different cluster, and we don't want to remove such
	 * high-degree nodes from the cluster.
	 *
	 * Remark: contraction is done to try and remove nodes in the "filaments" that connect 
	 * different clusters. Nodes that get unmarked after contraction are assumed not to be
	 * interesting or further processed by the caller (e.g. by creating kernels inside 
	 * them).
	 *
	 * Remark: the histogram of all ratios between the number of neighbors of a node in a 
	 * different nonempty cluster than the node, and the total number of neighbors of the
	 * node in a nonempty cluster, has a peak near zero corresponding to nodes that are 
	 * not in a cluster boundary. Except for that, the histogram has a complex shape in
	 * practice, with multiple peaks, so setting $boundaryThreshold$ to a value that
	 * separates the peak near zero from the rest might not be a good approach, because 
	 * later peaks might still correspond to boundary nodes.
	 *
	 * Remark: rather than using a threshold, we could avoid unmarking boundary nodes that
	 * are neither maximal nor minimal. This does not seem to provide any advantage in 
	 * practice, so it is not explored further.
	 *
	 * @param subgraph the interval $[first..last]$ is assumed to be sorted;
	 * @param clusterSize size of each cluster; the procedure updates these counts;
	 * @param outNeighbors,lastOutNeighbor,nInNeighbors,components,clusterTmp temporary 
	 * space, of size at least equal to the number of nodes in the subgraph;
	 * @param clusterField temporary space, of size at least equal to the total number of 
	 * nodes in the interval graph;
	 * @return the number of nodes removed by contraction.
	 */
	private static final int refineClusters_contract(int[] subgraph, int first, int last, int nClusters, int[] clusterSize, double boundaryThreshold, int[][] outNeighbors, int[] lastOutNeighbor, int[] nInNeighbors, int[] components, int[] clusterTmp, int[] clusterField, boolean sorted) {
		final int GROWTH_RATE = 10;  // Arbitrary
		int i, j;
		int to, cluster, newCluster, clusterFrom, size, node, nNodes, component, nComponents, out;
		int[] componentSize, componentCluster, largestComponent, largestComponentSize;
		
		// Removing boundary nodes
		out=0;
		nNodes=last-first+1;
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			if (neighborsInDifferentCluster(node,clusterField,sorted)<=boundaryThreshold) {
				clusterTmp[i-first]=cluster;
				continue;
			}
			lastTmpPoint=getAdjacentClusters(node,clusterField,sorted);
			if (lastTmpPoint==0 && tmpPoints[lastTmpPoint].position==cluster) clusterTmp[i-first]=cluster;
			else clusterTmp[i-first]=-1;
		}
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			newCluster=clusterTmp[i-first];
			clusterField[node]=newCluster;
			if (newCluster==-1) {
				clusterSize[cluster]--;
				out++;
			}
		}

		// Computing the connected components of each cluster after removal
		Math.set(lastOutNeighbor,nNodes-1,-1);
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			clusterFrom=clusterField[node];
			if (clusterFrom==-1) continue;
			for (j=0; j<nNeighbors[node]; j++) {
				if (!neighbors[node][j].on) {
					if (sorted) break;
					else continue;
				}
				to=neighbors[node][j].getTo(node);
				if (clusterField[to]!=clusterFrom) continue;
				lastOutNeighbor[i-first]++;
				if (lastOutNeighbor[i-first]==outNeighbors[i-first].length) {
					int[] newNeighbors = new int[outNeighbors[i-first].length+GROWTH_RATE];
					System.arraycopy(outNeighbors[i-first],0,newNeighbors,0,outNeighbors[i-first].length);
					outNeighbors[i-first]=newNeighbors;
				}
				outNeighbors[i-first][lastOutNeighbor[i-first]]=Arrays.binarySearch(subgraph,first,last+1,to)-first;
			}
		}
		for (i=0; i<nNodes; i++) lastOutNeighbor[i]++;
		Math.set(nInNeighbors,nNodes-1,0);
		nComponents=DAG.getConnectedComponents(nNodes,null,nInNeighbors,outNeighbors,lastOutNeighbor,components,stack);
		
		// Keeping just the largest component of each cluster
		componentSize = new int[nComponents];
		Math.set(componentSize,nComponents-1,0);
		componentCluster = new int[nComponents];
		Math.set(componentCluster,nComponents-1,-2);
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			if (cluster==-1) continue;
			component=components[i-first];
			componentSize[component]++;
			if (componentCluster[component]==-2) componentCluster[component]=cluster;
		}
		largestComponent = new int[nClusters];
		Math.set(largestComponent,nClusters-1,-1);
		largestComponentSize = new int[nClusters];
		Math.set(largestComponentSize,nClusters-1,0);
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			if (cluster==-1) continue;
			component=components[i-first];
			size=componentSize[component];
			if (size>largestComponentSize[cluster]) {
				largestComponentSize[cluster]=size;
				largestComponent[cluster]=component;
			}
		}
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			cluster=clusterField[node];
			if (cluster==-1) continue;
			if (components[i-first]!=largestComponent[cluster]) {
				clusterField[node]=-1;
				clusterSize[cluster]--;
				out++;
			}
		}
		return out;
	}
	
	
	/**
	 * Let $(C_0,...,C_{n-1})$ be the current partition of the subgraph into clusters, and 
	 * let $adjacentClusters[0..lastAdjacentCluster]$ be the sorted list of distinct 
	 * non-negative cluster IDs adjacent to $node$. The procedure stores in $deltas[j]$ 
	 * the difference between the modularity of the current partition, and the modularity 
	 * of a hypothetical new partition in which $node$ has been removed from its current 
	 * $C_i$ and assigned to $C_j$.
	 *
	 * For a definition of modularity, see e.g. 
	 * <a href="https://en.wikipedia.org/wiki/Modularity_(networks)">Wikipedia</a>.
	 *
	 * Remark: the procedure works also when $node$ does not belong to any $C_i$ and it 
	 * has to be added to a $C_i$.
	 *
	 * Remark: the procedure discards OFF edges, as well as edges to nodes with cluster
	 * smaller than -1.
	 *
	 * @param nEdges number of ON edges in the subgraph;
	 * @param degrees degree of every node in the interval graph;
	 * @param clusterField temporary space, of size at least equal to the total number of 
	 * nodes in the interval graph;
	 * @return if the maximum of array $deltas$ is achieved at position $i$, the procedure
	 * returns $i$; otherwise, the procedure returns the smallest $j$ with maximum 
	 * $deltas[j]$.
	 */
	private static final int deltaModularity(int node, int nEdges, double[] deltas, Point[] adjacentClusters, int lastAdjacentCluster, int[] clusterField, int[] degrees, boolean sorted) {
		int j, k;
		int to, degree, cluster, clusterTo, maxCluster;
		double d, maxDelta;
		Edge edge;
		if (lastAdjacentCluster==-1) return -1;
		
		// Computing $deltas$
		cluster=clusterField[node];
		degree=degrees[node];
		for (j=0; j<=lastAdjacentCluster; j++) deltas[(int)(adjacentClusters[j].position)]=0.0;
		for (j=0; j<nNeighbors[node]; j++) {
			edge=neighbors[node][j];
			if (!edge.on) {
				if (sorted) break;
				else continue;
			}
			to=edge.getTo(node);
			d=1.0-(degree*degrees[to])/(2.0*nEdges);
			clusterTo=clusterField[to];
			if (clusterTo<-1) continue;  // Discarding edges to clusterIDs<-1
			if (clusterTo==cluster) {
				for (k=0; k<=lastAdjacentCluster; k++) deltas[(int)(adjacentClusters[k].position)]-=d;
			}
			if (clusterTo>=0) deltas[clusterTo]+=d;
		}
		
		// Outputting
		maxDelta=Math.NEGATIVE_INFINITY; maxCluster=-1;
		for (j=0; j<=lastAdjacentCluster; j++) {
			k=(int)(adjacentClusters[j].position);
			if (deltas[k]>maxDelta) maxDelta=deltas[k];
		}
		maxCluster=-1;
		for (j=0; j<=lastAdjacentCluster; j++) {
			k=(int)(adjacentClusters[j].position);
			if (deltas[k]!=maxDelta) continue;
			if (maxCluster==-1) maxCluster=k;
			if (k==cluster) {
				maxCluster=cluster;
				break;
			}
		}
		return maxCluster;
	}
	
	
	/**
	 * Stores in $tmpPoints[0..X]$ the sorted list of non-negative distinct cluster
	 * IDs of the neighbors of $node$, where $X$ is returned in output.
	 *
	 * @return -1 if $node$ has no non-negative neighbor.
	 */
	private static final int getAdjacentClusters(int node, int[] clusterField, boolean sorted) {
		int j, to, clusterTo, lastAdjacentCluster;
		Edge edge;
		
		lastAdjacentCluster=-1;
		for (j=0; j<nNeighbors[node]; j++) {
			edge=neighbors[node][j];
			if (!edge.on) {
				if (sorted) break;
				else continue;
			}
			to=edge.getTo(node);
			clusterTo=clusterField[to];
			if (clusterTo<0) continue;
			lastAdjacentCluster++;
			ensureTmpPoints(lastAdjacentCluster);
			tmpPoints[lastAdjacentCluster].position=clusterTo;
			tmpPoints[lastAdjacentCluster].mass=1;
		}
		return lastAdjacentCluster<=0?lastAdjacentCluster:Points.sortAndCompact(tmpPoints,lastAdjacentCluster);
	}

	
	public static class SeedExtensionPair implements Comparable {
		public int node;
		public int cluster;
		public double delta;
		
		public SeedExtensionPair(int n, int c, double d) {
			node=n;
			cluster=c;
			delta=d;
		}
		
		public boolean equals(Object other) {
			SeedExtensionPair otherPair = (SeedExtensionPair)other;
			return otherPair.node==node;
		}
		
		/**
		 * Remark: pairs with larger $delta$ are considered smaller because
		 * $PriorityQueue.poll()$ returns a smallest element.
		 */
		public int compareTo(Object other) {
			SeedExtensionPair otherPair = (SeedExtensionPair)other;
			if (delta>otherPair.delta) return -1;
			else if (delta<otherPair.delta) return 1;
			return 0;
		}
	}
	
	
	/**
	 * Stores in $out$ a tree encoding of all the maximal cliques of the graph, using the
	 * $O(3^{n/3})$ algorithm and encoding described in \cite{tomita2006worst}, where $n$ 
	 * is the number of nodes.
	 *
	 * Remark: the procedure uses just ON edges, and it assumes them to be located in a 
	 * prefix of $neighbors[i]$ for all $i$, and to be sorted by the ID of their 
	 * destination node.
	 *
	 * Remark: the procedure uses the $visited$ field of nodes, and does not reset it to
	 * the value it had before the procedure was called.
	 *
	 * Remark: the procedure stores in global variable $maxCliqueSize$ the size of a 
	 * largest max clique.
	 *
	 * @param subgraph the interval $[first..last]$ is assumed to be sorted;
	 * @param out if NULL, the procedure does not write the cliques to $out$;
	 * @param maxSize if positive, the algorithm in \cite{tomita2006worst} is modified in 
	 * a similar way as described in \cite{zhang2005genome}, to enumerate all (not 
	 * necessarily maximal) cliques of size $maxSize$, and some (but not necessarily all) 
	 * maximal cliques of smaller size;
	 * @param samplingRate if positive, uses just a random subset, of this relative size,
	 * of the nodes at the root of the enumeration tree;
	 * @param maxNCliques the procedure stops after it finds at least this number of 
	 * maximal cliques;
	 * @return the number of max cliques found.
	 */
	public static final int getMaximalCliques(int[] subgraph, int first, int last, Stream out, int maxSize, double samplingRate, int maxNCliques, boolean verbose) {
		int i;
		int nNodes, boundarySize;
		
		nNodes=last-first+1;
		boundarySize=nNodes+getNOnEdges(subgraph,first,last,true);
		if (boundary==null || boundary.length<boundarySize) boundary = new int[boundarySize];
		System.arraycopy(subgraph,first,boundary,0,nNodes);
		if (marked==null || marked.length<boundarySize) marked = new boolean[boundarySize];
		lastTmpPoint=-1; 
		if (out!=null) out.clear(false);	
		maxCliqueSize=0;		
		return getMaximalCliques_impl(boundary,0,nNodes-1,nNodes,marked,out,0,false,maxSize,samplingRate,0,maxNCliques,verbose);
	}
	
	
	/**
	 * [USED FOR DEBUGGING ONLY]
	 * Resets the graph to Figure 3 in \cite{tomita2006worst}. 
	 * To be called inside $getMaximalCliques()$.
	 */
	private static final void testMaximalCliques() {
		nodesArray = new Node[10];
		for (int x=0; x<=9; x++) {
			nodesArray[x] = new Node();
			nodesArray[x].nodeID=x;
		}
		nNeighbors[0]=0;

		nNeighbors[1]=3;
		neighbors[1][0] = new Edge(); neighbors[1][0].nodeID1=1; neighbors[1][0].nodeID2=2;
		neighbors[1][1] = new Edge(); neighbors[1][1].nodeID1=1; neighbors[1][1].nodeID2=9;

		nNeighbors[2]=3;
		neighbors[2][0] = new Edge(); neighbors[2][0].nodeID1=2; neighbors[2][0].nodeID2=1;
		neighbors[2][1] = new Edge(); neighbors[2][1].nodeID1=2; neighbors[2][1].nodeID2=3;
		neighbors[2][2] = new Edge(); neighbors[2][2].nodeID1=2; neighbors[2][2].nodeID2=9;

		nNeighbors[9]=3;
		neighbors[9][0] = new Edge(); neighbors[9][0].nodeID1=9; neighbors[9][0].nodeID2=1;
		neighbors[9][1] = new Edge(); neighbors[9][1].nodeID1=9; neighbors[9][1].nodeID2=2;
		neighbors[9][2] = new Edge(); neighbors[9][2].nodeID1=9; neighbors[9][2].nodeID2=3;

		nNeighbors[3]=4;
		neighbors[3][0] = new Edge(); neighbors[3][0].nodeID1=3; neighbors[3][0].nodeID2=2;
		neighbors[3][1] = new Edge(); neighbors[3][1].nodeID1=3; neighbors[3][1].nodeID2=4;
		neighbors[3][2] = new Edge(); neighbors[3][2].nodeID1=3; neighbors[3][2].nodeID2=8;
		neighbors[3][3] = new Edge(); neighbors[3][3].nodeID1=3; neighbors[3][3].nodeID2=9;

		nNeighbors[8]=4;
		neighbors[8][0] = new Edge(); neighbors[8][0].nodeID1=8; neighbors[8][0].nodeID2=3;
		neighbors[8][1] = new Edge(); neighbors[8][1].nodeID1=8; neighbors[8][1].nodeID2=4;
		neighbors[8][2] = new Edge(); neighbors[8][2].nodeID1=8; neighbors[8][2].nodeID2=6;
		neighbors[8][3] = new Edge(); neighbors[8][3].nodeID1=8; neighbors[8][3].nodeID2=7;

		nNeighbors[4]=5;
		neighbors[4][0] = new Edge(); neighbors[4][0].nodeID1=4; neighbors[4][0].nodeID2=3;
		neighbors[4][1] = new Edge(); neighbors[4][1].nodeID1=4; neighbors[4][1].nodeID2=5;
		neighbors[4][2] = new Edge(); neighbors[4][2].nodeID1=4; neighbors[4][2].nodeID2=6;
		neighbors[4][3] = new Edge(); neighbors[4][3].nodeID1=4; neighbors[4][3].nodeID2=7;
		neighbors[4][4] = new Edge(); neighbors[4][4].nodeID1=4; neighbors[4][4].nodeID2=8;

		nNeighbors[5]=2;
		neighbors[5][0] = new Edge(); neighbors[5][0].nodeID1=5; neighbors[5][0].nodeID2=4;
		neighbors[5][1] = new Edge(); neighbors[5][1].nodeID1=5; neighbors[5][1].nodeID2=6;

		nNeighbors[6]=4;
		neighbors[6][0] = new Edge(); neighbors[6][0].nodeID1=6; neighbors[6][0].nodeID2=4;
		neighbors[6][1] = new Edge(); neighbors[6][1].nodeID1=6; neighbors[6][1].nodeID2=5;
		neighbors[6][2] = new Edge(); neighbors[6][2].nodeID1=6; neighbors[6][2].nodeID2=7;
		neighbors[6][3] = new Edge(); neighbors[6][3].nodeID1=6; neighbors[6][3].nodeID2=8;

		nNeighbors[7]=3;
		neighbors[7][0] = new Edge(); neighbors[7][0].nodeID1=7; neighbors[7][0].nodeID2=4;
		neighbors[7][1] = new Edge(); neighbors[7][1].nodeID1=7; neighbors[7][1].nodeID2=6;
		neighbors[7][2] = new Edge(); neighbors[7][2].nodeID1=7; neighbors[7][2].nodeID2=8;

		int nEdges = 0;
		for (int x=0; x<=9; x++) {
			for (int y=0; y<nNeighbors[x]; y++) neighbors[x][y].on=true;
			nEdges+=nNeighbors[x];
		}
		nEdges/=2;

		nNodes=10;
		boundary = new int[10+nEdges];
		for (int x=0; x<=9; x++) boundary[x]=x;
	}
	
	
	/**
	 * Scans the encoding of maximal cliques in $in$, and adds to the global variable 
	 * $seeds[0..X]$ all nodes in a clique of size at least $minCliqueSize$, iff at most a
	 * $intersectionRate$ fraction of the clique's nodes have already been used by other 
	 * stored cliques. Variable $clusters[0..X]$ contains the ID of the maximal clique 
	 * each node belongs to, where $X=out[0]$. $out[1]$ contains the number of distinct 
	 * cliques returned.
	 *
	 * Remark: the graph might contain clusters with very different size and clustering 
	 * coefficient, so taking e.g. the K largest cliques might not be useful, since they 
	 * might be mostly inside a small subset of all clusters.
	 *
	 * @param in the tree encoding described in \cite{tomita2006worst} produced by 
	 * procedure $getMaximalCliques()$;
	 * @return true if there is at least one clique of size $>=minCliqueSize$.
	 */
	private static final boolean markMaximalCliques(int[] subgraph, int first, int last, Stream in, int minCliqueSize, double intersectionRate, int[] out) {
		final int GROWTH_RATE = 100;  // Arbitrary
		final int ABSOLUTE_THRESHOLD = 2;  // Arbitrary
		int i, j;
		int top, lastClique, lastStoredClique, lastNode, element, nElements, size, nMarked;
		
		if (stack==null || stack.length<last-first+1) stack = new int[last-first+1];
		if (marked==null || marked.length<last-first+1) marked = new boolean[last-first+1];
		Math.set(marked,0,marked.length-1,false);
		top=-1; lastClique=-1; lastStoredClique=-1; lastNode=-1; nElements=in.nElements();
		for (i=0; i<nElements; i++) {
			element=in.getElementAt(i);
			if (element==MAXIMAL_CLIQUE_FOUND) {				
				lastClique++;
				size=top+1;
				if (size<minCliqueSize) {
					top--;
					continue;
				}
				nMarked=0;
				for (j=0; j<=top; j++) {
					if (marked[Arrays.binarySearch(subgraph,first,last+1,stack[j])]) nMarked++;
				}
				if (nMarked>Math.max(intersectionRate*size,ABSOLUTE_THRESHOLD)) {
					top--;
					continue;
				}
				lastStoredClique++;
				for (j=0; j<=top; j++) {
					marked[Arrays.binarySearch(subgraph,first,last+1,stack[j])]=true;
					lastNode++;
					if (seeds==null) seeds = new int[GROWTH_RATE];
					else if (lastNode==seeds.length) {
						int[] newSeeds = new int[seeds.length+GROWTH_RATE];
						System.arraycopy(seeds,0,newSeeds,0,seeds.length);
						seeds=newSeeds;
					}
					if (clusters==null) clusters = new int[GROWTH_RATE];
					else if (lastNode==clusters.length) {
						int[] newClusters = new int[clusters.length+GROWTH_RATE];
						System.arraycopy(clusters,0,newClusters,0,clusters.length);
						clusters=newClusters;
					}
					seeds[lastNode]=stack[j];
					clusters[lastNode]=lastStoredClique;
				}
				top--;
				continue;
			}
			else if (element==MAXIMAL_CLIQUE_UP) {
				top--;
				continue;
			}
			else if (element==MAXIMAL_CLIQUE_NOT_FOUND) {
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
				if (IO.CONSISTENCY_CHECKS) {
					for (int x=0; x<top; x++) {
						if (stack[x]==element) {
							System.err.println("ERROR: duplicated element in the stack?! 2");
							for (int y=0; y<=top; y++) System.err.print(stack[y]+", ");
							System.err.println();
							System.exit(1);
						}
					}
				}
			}
		}
		System.err.println("markMaximalCliques> "+(lastStoredClique+1)+" cliques of size >="+minCliqueSize+" stored (out of "+(lastClique+1)+" total cliques enumerated)");
		out[0]=lastNode; out[1]=lastStoredClique+1;
		return lastStoredClique>=0;
	}
	
	
	/**
	 * @return the size of a largest clique (if $mode=TRUE$), or the number of cliques (if
	 * $mode=FALSE$), stored in $in$.
	 */
	public static final int inspectMaximalCliques(Stream in, boolean mode, int first, int last) {
		int i;
		int top, element, out;
		final int nElements = in.nElements();
		
		if (stack==null || stack.length<last-first+1) stack = new int[last-first+1];
		out=0; top=-1;
		for (i=0; i<nElements; i++) {
			element=in.getElementAt(i);
			if (element==IntervalGraph.MAXIMAL_CLIQUE_FOUND) {
				if (mode) out=Math.max(out,top+1);
				else out++;
				top--;
				continue;
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
				stack[top]=element;
				if (IO.CONSISTENCY_CHECKS) {
					for (int x=0; x<top; x++) {
						if (stack[x]==element) {
							System.err.println("ERROR: duplicated element in the stack?! 2");
							for (int y=0; y<=top; y++) System.err.print(stack[y]+", ");
							System.err.println();
							System.exit(1);
						}
					}
				}
			}
		}
		return out;
	}

	
	/**
	 * Remark: the procedure stores in global variable $maxCliqueSize$ the size of a 
	 * largest max clique found.
	 * 	
	 * @param boundary sorted;
	 * @param out if NULL, the procedure does not write the cliques to $out$;
	 * @param maxSize if positive, the algorithm in \cite{tomita2006worst} is modified in 
	 * a similar way as described in \cite{zhang2005genome}, to enumerate all (not 
	 * necessarily maximal) cliques of size $maxSize$, and all maximal cliques of smaller
	 * size;
	 * @param samplingRate if positive, uses just a random subset, of this relative size,
	 * of the nodes at the root of the enumeration tree;
	 * @param maxNCliques the procedure stops when $previousNCliques$ plus the number of 
	 * cliques found by the procedure is at least $maxNCliques$;
	 * @return the number of max cliques found.
	 */
	private static final int getMaximalCliques_impl(int[] boundary, int firstBoundary, int lastBoundary, int nUnexplored, boolean[] marked, Stream out, int size, boolean logSize, int maxSize, double samplingRate, int previousNCliques, int maxNCliques, boolean verbose) {
		int i, j;
		int node, maxNode, degree, maxDegree, boundarySize;
		int firstBoundaryPrime, lastBoundaryPrime, nUnexploredPrime, result;
		
		result=0;
		boundarySize=lastBoundary-firstBoundary+1;
		if (boundarySize==0 || size==maxSize) {
			// Base case 1
			if (out!=null) out.push(MAXIMAL_CLIQUE_FOUND); 
			result++;
			maxCliqueSize=Math.max(maxCliqueSize,size);
			if (logSize) {
				lastTmpPoint++;
				ensureTmpPoints(lastTmpPoint);
				tmpPoints[lastTmpPoint].position=size;
				tmpPoints[lastTmpPoint].mass=1;
			}
			if (verbose && (previousNCliques+result)%1000000==0) System.err.println("found "+(previousNCliques+result)+" max. cliques");
			return result;
		}
		if (nUnexplored==0) {
			// Base case 2
			if (size>0 && out!=null) out.push(MAXIMAL_CLIQUE_NOT_FOUND);
			return result;
		}
	
		// Finding a node in $boundary$ with the largest number of non-explored
		// neighbors.
		if (boundarySize>1) {
			maxDegree=0; maxNode=-1;
			for (i=firstBoundary; i<=lastBoundary; i++) {
				node=boundary[i];
				if (node<0) node=-1-node;
				degree=size==0?nUnexploredNeighbors_root(node):nUnexploredNeighbors(node,boundary,firstBoundary,lastBoundary);
				if (degree>maxDegree) {
					maxDegree=degree;
					maxNode=node;
				}
			}
			if (maxDegree==nUnexplored) {
				// Base case 4
				if (size>0 && out!=null) out.push(MAXIMAL_CLIQUE_NOT_FOUND);
				return result;
			}
			if (maxNode!=-1) setIntersectionMark(boundary,firstBoundary,lastBoundary,maxNode,marked);
		}
	
		// Recursion
		for (i=firstBoundary; i<=lastBoundary && previousNCliques+result<maxNCliques; i++) {
			if (boundary[i]<0 || marked[i] || (size==0 && samplingRate>0 && Math.random.nextDouble()>samplingRate)) continue;
			node=boundary[i];
			if (out!=null) out.push(node);
			firstBoundaryPrime=lastBoundary+1;
			setIntersection(boundary,firstBoundary,lastBoundary,node,boundary,firstBoundaryPrime,tmpArray);
			lastBoundaryPrime=tmpArray[0];
			nUnexploredPrime=tmpArray[1];
			Math.set(marked,firstBoundaryPrime,lastBoundaryPrime,false);
			result+=getMaximalCliques_impl(boundary,firstBoundaryPrime,lastBoundaryPrime,nUnexploredPrime,marked,out,size+1,logSize,maxSize,samplingRate,previousNCliques+result,maxNCliques,verbose);
			boundary[i]=-1-boundary[i];
		}
		if (size>0 && out!=null) out.push(MAXIMAL_CLIQUE_UP);
		return result;
	}
	
	
	/**
	 * A variant of $Math.setIntersection$ in which: (1) $x2=neighbors[node]$, which is
	 * assumed to have all ON edges in a prefix, sorted by node ID; (2) $x1[from1..
	 * last1]$ can have negative values $x$, which are interpreted as positive by taking
	 * $-1-x$; $x1[from1..last1]$ is sorted by such interpretation; (3) negative values
	 * are copied to $y$.
	 *
	 * @param out output array: 0=last element of $y$ used; 1=number of non-negative 
	 * values copied to $y$.
	 */
	private static final void setIntersection(int[] x1, int from1, int last1, int node, int[] y, int fromY, int[] out) {
		int i1, i2, j, to, x1value, nNonNegative;
		Edge edge;
		
		i1=from1; i2=0;
		j=fromY; nNonNegative=0;
		while (i1<=last1 && i2<nNeighbors[node]) {
			edge=neighbors[node][i2];
			if (!edge.on) break;
			to=edge.getTo(node);
			x1value=x1[i1];
			if (x1value<0) x1value=-1-x1value;
			if (x1value<to) i1++;
			else if (x1value>to) i2++;
			else {
				y[j++]=x1[i1];
				if (x1[i1]>=0) nNonNegative++;
				i1++; i2++;
			}
		}
		out[0]=j-1;
		out[1]=nNonNegative;
	}
	
	
	/**
	 * The procedure assumes that: (1) all the ON edges of $neighbors[i]$ are in a prefix, 
	 * sorted by node ID; (2) all and only the currently unexplored nodes of 
	 * $boundary[firstBoundary..lastBoundary]$ are non-negative.
	 */
	private static final int nUnexploredNeighbors(int node, int[] boundary, int firstBoundary, int lastBoundary) {
		int i1, i2, size, to;
		Edge edge;
		
		i1=0; i2=firstBoundary; size=0;
		while (i1<nNeighbors[node] && i2<=lastBoundary) {
			edge=neighbors[node][i1];
			if (!edge.on) break;
			if (boundary[i2]<0) {
				i2++;
				continue;
			}
			to=edge.getTo(node);
			if (boundary[i2]<to) i2++;
			else if (boundary[i2]>to) i1++;
			else {
				size++;
				i1++; i2++;
			}
		}
		return size;
	}
	
	
	/**
	 * Variant of $nUnexploredNeighbors$ for the root node (emtpy current clique).
	 */
	private static final int nUnexploredNeighbors_root(int node) {
		if (nNeighbors[node]==0) return 0;
		if (neighbors[node][nNeighbors[node]-1].on) return nNeighbors[node];
		for (int i=0; i<nNeighbors[node]; i++) {
			if (!neighbors[node][i].on) return i;
		}
		return -1;
	}
	
	
	/**
	 * A variant of $Math.setIntersectionSize()$ in which $x2=neighbors[node]$, which is
	 * assumed to have all ON edges in a prefix, sorted by node ID.
	 */
	private static final int setIntersectionSize(int[] x1, int from1, int last1, int node) {
		int i1, i2, to, size;
		Edge edge;
		
		i1=from1; i2=0; size=0;
		while (i1<=last1 && i2<nNeighbors[node]) {
			edge=neighbors[node][i2];
			if (!edge.on) break;
			to=edge.getTo(node);
			if (x1[i1]<to) i1++;
			else if (x1[i1]>to) i2++;
			else {
				size++;
				i1++; i2++;
			}
		}
		return size;
	}
	
	
	/**
	 * Sets $marked[i]=true$ for all $i \in [from1..last1]$ such that $x1[i]$ is non-
	 * negative and belongs to the intersection of $x1[from1..last1]$ and the neighbors of
	 * $node$.
	 */
	private static final void setIntersectionMark(int[] x1, int from1, int last1, int node, boolean[] marked) {
		int i1, i2, to;
		Edge edge;
		
		i1=from1; i2=0;
		while (i1<=last1 && i2<nNeighbors[node]) {
			edge=neighbors[node][i2];
			if (!edge.on) break;
			if (x1[i1]<0) {
				i1++;
				continue;
			}
			to=edge.getTo(node);
			if (to==node) {
				i2++;
				continue;
			}
			if (x1[i1]<to) i1++;
			else if (x1[i1]>to) i2++;
			else {
				marked[i1]=true;
				i1++; i2++;
			}
		}
		return;
	}
	
	
	/**
	 * Assume that no node in $subgraph[first..last]$ is assigned a nonnegative component
	 * ID, except for the nodes $x$ such that $x=subgraph[i]$ and $target[x]>=0$. The 
	 * procedure sets the $components$ field of every node $y$ such that $y=subgraph[j]$ 
	 * and $target[y] < 0$, by propagating the $components$ field of $x$ nodes via ON
	 * edges of containment type only.
	 *
	 * @param tmp temporary space with at least M elements, where M is the maximum
	 * number of components per node.
	 */
	public static final void propagateComponentTags(int[] subgraph, int first, int last, int[] target, int[] tmp, boolean sorted) {
		int i, j;
		int source, from, to, top, lastComponent, previousLastComponent, type;
		Node node, nodeTo;
		Edge edge;
		
		if (stack==null || stack.length<nNodes) stack = new int[nNodes];
		for (i=first; i<=last; i++) nodesArray[subgraph[i]].visited=-1;
		for (i=first; i<=last; i++) {
			source=subgraph[i];
			if (target[source]<0) continue;
			node=nodesArray[source];
			node.visited=source;
			top=0;
			stack[0]=source;
			while (top>=0) {
				from=stack[top--];
				for (j=0; j<nNeighbors[from]; j++) {
					edge=neighbors[from][j];
					if (!edge.on) {
						if (sorted) break;
						else continue;
					}
					type=edge.containment;
					if (type==-1) continue;
					to=edge.getTo(from);
					if (target[to]>=0) continue;
					nodeTo=nodesArray[to];
					if (nodeTo.visited==source) continue;
					previousLastComponent=nodeTo.lastComponent;
					lastComponent=Math.setUnion(node.components,node.lastComponent,nodeTo.components,nodeTo.lastComponent,tmp);
					if (lastComponent==previousLastComponent) continue;
					if (nodeTo.components.length<lastComponent+1) nodeTo.components=Arrays.copyOf(tmp,lastComponent+1);
				    else System.arraycopy(tmp,0,nodeTo.components,0,lastComponent+1);
					nodeTo.lastComponent=lastComponent;
					nodeTo.visited=source;
					stack[++top]=to;
				}
			}
		}
	}
	
	
	/**
	 * Implements the $O(nm)$-time algorithm in \cite{brandes2001faster}.
	 * This procedure is not currently used (see below), but it is kept here for the 
	 * future.
	 *
	 * Remark: the procedure assumes that the subgraph induced by nodes in 
	 * $subgraph[first..last]$ is connected to the rest of the graph just by OFF edges. 
	 * The subgraph itself need not be connected.
	 *
	 * Remark: iteratively splitting the interval graph by removing vertices of high 
	 * betweenness is not likely to work in practice, since betweenness does not seem to 
	 * highlight vertices that connect different clusters. However, we never implemented 
	 * the full algorithm, and it might be that recursive splitting by betweenness works 
	 * even though betweenness does not seem to highlight bridges at the first iteration.
	 *
	 * Remark: another heuristic for a top-down splitting strategy could be removing
	 * nodes that connect biconnected components. We didn't try this approach, since 
	 * clusters seem to have multiple connections with one another.
	 *
	 * @param sorted TRUE iff edges with $on=true$ are at the beginning of each row of
	 * $neighbors$;
	 * @param tmpQueue,tmpOrder,lastInNeighbor,inNeighbors temporary space, of size at 
     * least equal to the total number of nodes in the interval graph. Row $i$ of
	 * $inNeighbors$ is also assumed to have at least $nNeighbors[i]$ elements;
	 * @param nShortestPaths,shortestPath,delta temporary space, of size at least equal to
	 * the total number of nodes in the interval graph;
	 * @param centrality stores the centrality of each node; must be of size at least 
	 * equal to the total number of nodes in the interval graph.
	 */
	public static final void betweennessCentrality(int[] subgraph, int first, int last, int[] tmpQueue, int[] tmpOrder, int[][] inNeighbors, int[] lastInNeighbor, int[] nShortestPaths, int[] shortestPath, double[] delta, double[] centrality, boolean sorted) {
		final int queueLength = tmpQueue.length;
		int i, j, k;
		int from, to, firstInQueue, lastInQueue, queueSize, lastInOrder;
		double min, max;
		Node nodeFrom, nodeTo;
		Edge tmpEdge;
		
		for (i=first; i<=last; i++) {
			nodeFrom=nodesArray[subgraph[i]];
			centrality[subgraph[i]]=0.0;
			nodeFrom.visited=-1;
		}
		for (i=first; i<=last; i++) {		
			// Storing in $tmpOrder$ the topologically-sorted DAG of shortest paths from
			// the current node (excluded).
			from=subgraph[i];
			nodeFrom=nodesArray[from];
			shortestPath[from]=0;
			nShortestPaths[from]=1;
			nodeFrom.visited=subgraph[i];
			tmpQueue[0]=subgraph[i]; firstInQueue=0; lastInQueue=0; queueSize=1;
			lastInOrder=-1;
			while (queueSize>0) {
				from=tmpQueue[firstInQueue]; firstInQueue=(firstInQueue+1)%queueLength; queueSize--;
				nodeFrom=nodesArray[from];
				for (j=0; j<nNeighbors[from]; j++) {
					tmpEdge=neighbors[from][j];
					if (!tmpEdge.on) {
						if (sorted) break;
						else continue;
					}
					to=tmpEdge.getTo(from);
					nodeTo=nodesArray[to];
					if (nodeTo.visited!=subgraph[i]) {
						shortestPath[to]=shortestPath[from]+1;
						nShortestPaths[to]=nShortestPaths[from];
						nodeTo.visited=subgraph[i];
						lastInNeighbor[to]=0;
						inNeighbors[to][0]=from;
						queueSize++; lastInQueue=(lastInQueue+1)%queueLength; tmpQueue[lastInQueue]=to;
						tmpOrder[++lastInOrder]=to;				
					}
					else if (shortestPath[to]==shortestPath[from]+1) {
						nShortestPaths[to]+=nShortestPaths[from];
						inNeighbors[to][++lastInNeighbor[to]]=from;
					}
				}
			}
			
			// Propagating $\delta_{s\bullet}(v)$ values backwards in the DAG
			for (j=0; j<=lastInOrder; j++) delta[tmpOrder[j]]=0.0;
			for (j=lastInOrder; j>=0; j--) {
				to=tmpOrder[j];
				for (k=0; k<=lastInNeighbor[to]; k++) {
					from=inNeighbors[to][k];
					delta[from]+=((1.0+delta[to])*nShortestPaths[from])/nShortestPaths[to];
				}
			}
			
			// Cumulating delta to centrality on all nodes (excluding the source)
			for (j=0; j<=lastInOrder; j++) centrality[tmpOrder[j]]+=delta[tmpOrder[j]];
		}
		min=Math.POSITIVE_INFINITY; max=0;
		for (i=first; i<=last; i++) {
			centrality[subgraph[i]]/=2.0;
			if (centrality[subgraph[i]]<min) min=centrality[subgraph[i]];
			if (centrality[subgraph[i]]>max) max=centrality[subgraph[i]];
		}
		for (i=first; i<=last; i++) centrality[subgraph[i]]=(centrality[subgraph[i]]-min)/(max-min);
	}
	
	
	/**
	 * Assumes that the set of all edges on which centrality was computed are ON.
	 */
	public static final void turnOffCentralNodes(double threshold, double[] centrality) {
		int i, j;
		int nEdges, nEdgesOn;
		
		nEdges=0; nEdgesOn=0;
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].on) {
					nEdges++;
					if (centrality[i]>threshold) neighbors[i][j].on=false;
					else nEdgesOn++;
				}
			}
		}
		System.err.println("nEdges="+(nEdges/2)+" nEdgesOn="+(nEdgesOn/2));
	}
	
	
	
	
	
	
	
	
	// ---------------- RESOLVING DANGLING EDGES AND UNASSIGNED INTERVALS ----------------
	
	/**
	 * Tries to add to the graph all the alignments that have been assigned to an interval
	 * in just one of readA and readB. See $fixDanglingEdges_impl$ for details.
	 *
	 * Remark: the procedure can increase both $edges$ and $nodes$, if intervals had not 
	 * been assigned to any alignment in the original factorization. Such nodes could 
	 * become components of size one in $IntervalGraphStep1$ and be subsequently 
	 * discarded.
	 *
	 * @param alignmentPairs the procedure assumes that the set of all and only dangling 
	 * edges is stored in this data structure;
	 * @param nodes,edges nodes and edges already added to the graph.
	 */
	private static final void fixDanglingEdges(HashMap<String,AlignmentPair> alignmentPairs, HashMap<Node,Node> nodes, HashMap<Edge,Edge> edges, String periodicSubstringsFile, String denseSubstringsFile, String alignmentsFile, int maxIntervalsPerRead, boolean discardContainedInShortPeriod, boolean discardContainedInLongPeriod, boolean discardPeriodicContainedInPeriodic, boolean discardDenseInDense) throws IOException {
		final int READ_AHEAD_LIMIT = maxIntervalsPerRead*100;  // Characters
		int i;
		int read, nDangling, previousRead, firstDangling, lastDangling;
		int lastPeriodicInterval, lastDenseInterval, lastAlignmentInterval;
		int fixed, fixedPeriodic, fixedDense, fixedAlignment, fixedNodes;
		Set<Map.Entry<String,AlignmentPair>> entries;
		Iterator<Map.Entry<String,AlignmentPair>> iterator;
		Map.Entry<String,AlignmentPair> entry;
		AlignmentPair pair;
		int[] out = new int[4];
		String[] tokens;
		DanglingEdge[] dangling;
		int[][] periodicIntervals, denseIntervals, alignmentIntervals;
		BufferedReader periodicBuffer = new BufferedReader(new FileReader(periodicSubstringsFile),IO.BUFFER_SIZE);
		periodicBuffer.mark(READ_AHEAD_LIMIT);
		BufferedReader denseBuffer = new BufferedReader(new FileReader(denseSubstringsFile),IO.BUFFER_SIZE);
		denseBuffer.mark(READ_AHEAD_LIMIT);
		BufferedReader alignmentsBuffer = new BufferedReader(new FileReader(alignmentsFile),IO.BUFFER_SIZE);
		alignmentsBuffer.mark(READ_AHEAD_LIMIT);
		
		// Building the sorted list of all dangling edges
		entries=alignmentPairs.entrySet();
		nDangling=entries.size();
		dangling = new DanglingEdge[nDangling];
		iterator=entries.iterator();
		i=-1;
		while (iterator.hasNext()) {
			entry=iterator.next();
			tokens=entry.getKey().split(",");
			i++;
			pair=entry.getValue();
			if (pair.node1==null) read=Integer.parseInt(tokens[0]);
			else read=Integer.parseInt(tokens[3]);
			dangling[i] = new DanglingEdge(pair,read);
		}
		if (dangling.length>1) Arrays.sort(dangling);
		
		// Fixing dangling edges
		fixedPeriodic=0; fixedDense=0; fixedAlignment=0; fixedNodes=0;
		periodicIntervals = new int[maxIntervalsPerRead][11];
		lastPeriodicInterval=-1;
		denseIntervals = new int[maxIntervalsPerRead][14];
		lastDenseInterval=-1;
		alignmentIntervals = new int[maxIntervalsPerRead][8];
		lastAlignmentInterval=-1;
		previousRead=dangling[0].read; 
		firstDangling=0; lastDangling=0; i=1;
		while (i<nDangling) {
			if (dangling[i].read!=previousRead) {
				lastPeriodicInterval=PeriodicSubstringInterval.loadPeriodicIntervals(previousRead,periodicBuffer,READ_AHEAD_LIMIT,periodicIntervals,tmpArray);
				lastDenseInterval=DenseSubstring.loadDenseIntervals(previousRead,denseBuffer,READ_AHEAD_LIMIT,denseIntervals,tmpArray);
				lastAlignmentInterval=AlignmentInterval.loadAlignmentIntervals(previousRead,alignmentsBuffer,READ_AHEAD_LIMIT,alignmentIntervals,tmpArray);
				fixDanglingEdges_impl(dangling,firstDangling,lastDangling,periodicIntervals,lastPeriodicInterval,denseIntervals,lastDenseInterval,alignmentIntervals,lastAlignmentInterval,nodes,edges,out,discardContainedInShortPeriod,discardContainedInLongPeriod,discardPeriodicContainedInPeriodic,discardDenseInDense);
				fixedPeriodic+=out[0]; fixedDense+=out[1]; fixedAlignment+=out[2]; fixedNodes+=out[3];
				previousRead=dangling[i].read;
				firstDangling=i; lastDangling=i;
			}
			else lastDangling=i; 
			i++;
		}
		lastPeriodicInterval=PeriodicSubstringInterval.loadPeriodicIntervals(previousRead,periodicBuffer,READ_AHEAD_LIMIT,periodicIntervals,tmpArray);
		lastDenseInterval=DenseSubstring.loadDenseIntervals(previousRead,denseBuffer,READ_AHEAD_LIMIT,denseIntervals,tmpArray);
		lastAlignmentInterval=AlignmentInterval.loadAlignmentIntervals(previousRead,alignmentsBuffer,READ_AHEAD_LIMIT,alignmentIntervals,tmpArray);
		fixDanglingEdges_impl(dangling,firstDangling,lastDangling,periodicIntervals,lastPeriodicInterval,denseIntervals,lastDenseInterval,alignmentIntervals,lastAlignmentInterval,nodes,edges,out,discardContainedInShortPeriod,discardContainedInLongPeriod,discardPeriodicContainedInPeriodic,discardDenseInDense);
		fixedPeriodic+=out[0]; fixedDense+=out[1]; fixedAlignment+=out[2]; fixedNodes+=out[3];
		periodicBuffer.close(); denseBuffer.close(); alignmentsBuffer.close();
		fixed=fixedPeriodic+fixedDense+fixedAlignment;
		System.err.println("Fixed "+fixed+" dangling alignments ("+IO.getPercent(fixed,nDangling)+"%): "+fixedPeriodic+" periodic ("+IO.getPercent(fixedPeriodic,fixed)+"% of fixed); "+fixedDense+" dense ("+IO.getPercent(fixedDense,fixed)+"% of fixed); "+fixedAlignment+" alignment ("+IO.getPercent(fixedAlignment,fixed)+"% of fixed); new nodes created: "+fixedNodes);
	}
	
	
	/**
	 * Let $dangling[first..last]$ contain all dangling edges in a read R. The procedure 
	 * tries to assign each dangling edge in $dangling[first..last]$ to the only interval 
	 * in R that satisfies $IntervalGraphStep2.inStep2$, adding the newly formed edge to 
	 * $edges$. If multiple such intervals exist, the procedure does not assign the 
	 * dangling edge.
	 *
	 * Remark: intervals that were not assigned to any alignment in the connection file
	 * might be added as new nodes to the graph, if they resolve a dangling edge.
	 *
	 * Remark: the procedure assumes $dangling[first..last]$ to be sorted by starting 
	 * position.
	 *
	 * @param out number of dangling edges assigned by the procedure to a periodic (0), 
	 * dense (1), or alignment (2) interval; and number of new nodes created (3).
	 */
	private static final void fixDanglingEdges_impl(DanglingEdge[] dangling, int first, int last, int[][] periodicIntervals, int lastPeriodicInterval, int[][] denseIntervals, int lastDenseInterval, int[][] alignmentIntervals, int lastAlignmentInterval, HashMap<Node,Node> nodes, HashMap<Edge,Edge> edges, int[] out, boolean discardContainedInShortPeriod, boolean discardContainedInLongPeriod, boolean discardPeriodicContainedInPeriodic, boolean discardDenseInDense) {
		final int MIN_ALIGNMENT_LENGTH = Alignments.minAlignmentLength-IO.quantum;  // Arbitrary
		boolean doPeriodic, doDense, doAlignments;
		int i, j, k, h;
		int fixedType, fixedID, edgeType;
		int firstJForNextI, firstKForNextI, firstHForNextI;
		int fixedPeriodic, fixedDense, fixedAlignment, newNodesCreated;
		int nullStart, nullEnd;
		Node tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Node tmpNodePrime;
		Edge queryEdge = new Edge();
		Edge tmpEdge;
		
		fixedPeriodic=0; fixedDense=0; fixedAlignment=0; newNodesCreated=0;
		j=0; doPeriodic=j<=lastPeriodicInterval;
		k=0; doDense=k<=lastDenseInterval;
		h=0; doAlignments=h<=lastAlignmentInterval;
		fixedType=-1; fixedID=-1;
		i=first; firstJForNextI=-1; firstKForNextI=-1; firstHForNextI=-1; 
		while (i<=last) {
			nullStart=dangling[i].getNullStart();
			nullEnd=dangling[i].getNullEnd();
			if (!doPeriodic && !doDense && !doAlignments) {
				if (fixedType>=0) {
					if (fixedType==Constants.INTERVAL_PERIODIC) fixedPeriodic++;
					else if (fixedType>=Constants.INTERVAL_DENSE_PREFIX && fixedType<=Constants.INTERVAL_DENSE_SINGLEDELETION) fixedDense++;
					else fixedAlignment++;
					tmpNodePrime=dangling[i].getNullNode();
					if (tmpNodePrime.nodeID==-1) {
						tmpNodePrime.nodeID=nodes.size();
						nodes.put(tmpNodePrime,tmpNodePrime);
						newNodesCreated++;
					}
					dangling[i].pair.trim();
					if (dangling[i].pair.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && dangling[i].pair.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
						queryEdge.set(dangling[i].pair);
						tmpEdge=edges.get(queryEdge);
						if (tmpEdge==null) {
							tmpEdge = new Edge(dangling[i].pair);
							edges.put(tmpEdge,tmpEdge);
						}
						else {
							tmpEdge.addTypes(queryEdge);
							tmpEdge.addDiffs(queryEdge);
							tmpEdge.addOrientation(queryEdge);
						}
					}
				}
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				doPeriodic=j<=lastPeriodicInterval;
				if (firstKForNextI!=-1) k=firstKForNextI; 
				firstKForNextI=-1;
				doDense=k<=lastDenseInterval;
				if (firstHForNextI!=-1) h=firstHForNextI; 
				firstHForNextI=-1;
				doAlignments=h<=lastAlignmentInterval;
				fixedType=-1; fixedID=-1;
				continue;
			}
			
			if (doPeriodic) {
				if (firstJForNextI==-1 && i<last && periodicIntervals[j][2]>=dangling[i+1].getNullStart()) firstJForNextI=j;
				if ( periodicIntervals[j][1]<nullEnd && 						 
					 ((discardContainedInShortPeriod||discardPeriodicContainedInPeriodic)?(periodicIntervals[j][5]&(1<<Constants.INTERVAL_PERIODIC))==0:true) &&
					 ((discardContainedInLongPeriod||discardPeriodicContainedInPeriodic)?(periodicIntervals[j][5]&(1<<(Constants.INTERVAL_PERIODIC+1)))==0:true)
			       ) {
					if ( Intervals.isApproximatelyContained(nullStart,nullEnd,periodicIntervals[j][1],periodicIntervals[j][2]) ||
						 Intervals.areApproximatelyIdentical(nullStart,nullEnd,periodicIntervals[j][1],periodicIntervals[j][2])
					   ) {
	   					tmpNode.read=dangling[i].read;
	   					tmpNode.type=Constants.INTERVAL_PERIODIC;
	   					tmpNode.id=periodicIntervals[j][0];
	   					tmpNodePrime=nodes.get(tmpNode);
						if (tmpNodePrime==null) {
							tmpNodePrime = new Node(tmpNode.read,tmpNode.type,tmpNode.id,periodicIntervals[j][1],periodicIntervals[j][2],periodicIntervals[j][3]==1,periodicIntervals[j][4]==1,false,periodicIntervals[j][5],-1,-1,periodicIntervals[j][6],periodicIntervals[j][7]==1,-1,-1,periodicIntervals[j][9],periodicIntervals[j][10]);
							tmpNodePrime.nodeID=-1;
						}
   						dangling[i].setNullNode(tmpNodePrime);
   						edgeType=IntervalGraphStep2.inStep2(dangling[i].pair);
   						if (edgeType>=0 && edgeType!=Constants.SHARED_SUBSTRING && edgeType!=Constants.INSERTION_BOTH) {
   							if (fixedType==-1) {
   								fixedType=Constants.INTERVAL_PERIODIC;
   								fixedID=periodicIntervals[j][0];
   							}
   							else {
   								fixedType=-2;
   								doPeriodic=false;
   								doDense=false;
   								doAlignments=false;
   							}
   						}
					}
					j++;
					if (j==lastPeriodicInterval+1) doPeriodic=false;
				}
				else doPeriodic=false;
			}
			
			if (doDense) {
				if (firstKForNextI==-1 && i<last && denseIntervals[k][3]>=dangling[i+1].getNullStart()) firstKForNextI=k;
				if ( denseIntervals[k][2]<nullEnd &&
				     (discardContainedInShortPeriod?(denseIntervals[k][7]&(1<<Constants.INTERVAL_PERIODIC))==0:true) &&
				     (discardContainedInLongPeriod?(denseIntervals[k][7]&(1<<(Constants.INTERVAL_PERIODIC+1)))==0:true) &&
					 (discardDenseInDense?(denseIntervals[k][7]&(1<<Constants.INTERVAL_DENSE_SUBSTRING))==0:true)
				   ) {
					if ( Intervals.isApproximatelyContained(nullStart,nullEnd,denseIntervals[k][2],denseIntervals[k][3]) ||
						 Intervals.areApproximatelyIdentical(nullStart,nullEnd,denseIntervals[k][2],denseIntervals[k][3])
					   ) {
						tmpNode.read=dangling[i].read;
						tmpNode.type=denseIntervals[k][1];
						tmpNode.id=denseIntervals[k][0];
						tmpNodePrime=nodes.get(tmpNode);
						if (tmpNodePrime==null) {
							tmpNodePrime = new Node(tmpNode.read,tmpNode.type,tmpNode.id,denseIntervals[k][2],denseIntervals[k][3],denseIntervals[k][4]==1,denseIntervals[k][5]==1,denseIntervals[k][6]==1,denseIntervals[k][7],denseIntervals[k][8],denseIntervals[k][9],-1,false,denseIntervals[k][10],denseIntervals[k][11],denseIntervals[k][12],denseIntervals[k][13]);
							tmpNodePrime.nodeID=-1;
						}
						dangling[i].setNullNode(tmpNodePrime);
						edgeType=IntervalGraphStep2.inStep2(dangling[i].pair);
						if (edgeType>=0 && edgeType!=Constants.SHARED_SUBSTRING && edgeType!=Constants.INSERTION_BOTH) {
							if (fixedType==-1) {
								fixedType=denseIntervals[k][1];
								fixedID=denseIntervals[k][0];
							}
							else {
								fixedType=-2;
								doPeriodic=false;
								doDense=false;
								doAlignments=false;
							}
						}
					}
					k++;
					if (k==lastDenseInterval+1) doDense=false;
				}
				else doDense=false;
			}
			
			if (doAlignments) {
				if (firstHForNextI==-1 && i<last && alignmentIntervals[h][2]>=dangling[i+1].getNullStart()) firstHForNextI=h;
				if ( alignmentIntervals[h][1]<nullEnd &&
				     (discardContainedInShortPeriod?(alignmentIntervals[h][5]&(1<<Constants.INTERVAL_PERIODIC))==0:true) &&
				     (discardContainedInLongPeriod?(alignmentIntervals[h][5]&(1<<(Constants.INTERVAL_PERIODIC+1)))==0:true)
				   ) {
					if ( Intervals.isApproximatelyContained(nullStart,nullEnd,alignmentIntervals[h][1],alignmentIntervals[h][2]) ||
						 Intervals.areApproximatelyIdentical(nullStart,nullEnd,alignmentIntervals[h][1],alignmentIntervals[h][2])
					   ) {
						tmpNode.read=dangling[i].read;
						tmpNode.type=Constants.INTERVAL_ALIGNMENT;
						tmpNode.id=alignmentIntervals[h][0];
						tmpNodePrime=nodes.get(tmpNode);
						if (tmpNodePrime==null) {
							tmpNodePrime = new Node(tmpNode.read,tmpNode.type,tmpNode.id,alignmentIntervals[h][1],alignmentIntervals[h][2],alignmentIntervals[h][3]==1,alignmentIntervals[h][4]==1,false,alignmentIntervals[h][5],-1,-1,-1,false,-1,-1,alignmentIntervals[h][6],alignmentIntervals[h][7]);
							tmpNodePrime.nodeID=-1;
						}
						dangling[i].setNullNode(tmpNodePrime);
						edgeType=IntervalGraphStep2.inStep2(dangling[i].pair);
						if (edgeType>=0 && edgeType!=Constants.SHARED_SUBSTRING && edgeType!=Constants.INSERTION_BOTH) {
							if (fixedType==-1) {
								fixedType=Constants.INTERVAL_ALIGNMENT;
								fixedID=alignmentIntervals[h][0];
							}
							else {
								fixedType=-2;
								doPeriodic=false;
								doDense=false;
								doAlignments=false;
							}
						}
					}
					h++;
					if (h==lastAlignmentInterval+1) doAlignments=false;
				}
				else doAlignments=false;
			}
		}
		out[0]=fixedPeriodic; out[1]=fixedDense; out[2]=fixedAlignment; out[3]=newNodesCreated;
	}
	
	
	private static class DanglingEdge implements Comparable {
		public AlignmentPair pair;
		public int read;  // Of the null node. Zero-based.
		public boolean nullNode;  // Of $pair$. TRUE=1, FALSE=2.
		
		public DanglingEdge(AlignmentPair p, int r) {
			pair=p; 
			read=r;
			nullNode=p.node1==null;
		}
		
		public final int getNullStart() {
			return nullNode?pair.alignmentStart1:pair.alignmentStart2;
		}
		
		public final int getNullEnd() {
			return nullNode?pair.alignmentEnd1:pair.alignmentEnd2;
		}
		
		public final void setNullNode(Node newNode) {
			if (nullNode) pair.node1=newNode;
			else pair.node2=newNode;
		}
		
		public final Node getNullNode() {
			if (nullNode) return pair.node1;
			else return pair.node2;
		}
		
		public int compareTo(Object other) {
			DanglingEdge otherEdge = (DanglingEdge)other;
			if (read<otherEdge.read) return -1;
			else if (read>otherEdge.read) return 1;
			int start = getNullStart();
			int otherStart = otherEdge.getNullStart();
			if (start<otherStart) return -1;
			else if (start>otherStart) return 1;
			return 0;
		}
	}
	
	
	/**
	 * Tries to add to the graph all the alignments that have been assigned to an interval
	 * in neither readA nor readB. See $fixUnassignedAlignments_impl()$ for details.
	 *
	 * Remark: the procedure can increase both $edges$ and $nodes$, if intervals had not 
	 * been assigned to any alignment during factorization. Such nodes could become 
	 * components of size one in $IntervalGraphStep1$ and be eventually discarded.
	 *
	 * @param unassignedAlignmentPairs the procedure assumes the set of all and only 
	 * alignments that have not been assigned in either read to be stored in this data 
	 * structure;
	 * @param nodes,edges nodes and edges already added to the graph.
	 */
	private static final void fixUnassignedAlignments(HashMap<String,UnassignedAlignmentPair> unassignedAlignmentPairs, HashMap<Node,Node> nodes, HashMap<Edge,Edge> edges, String periodicSubstringsFile, String denseSubstringsFile, String alignmentsFile, int maxIntervalsPerRead, boolean discardContainedInShortPeriod, boolean discardContainedInLongPeriod, boolean discardPeriodicContainedInPeriodic, boolean discardDenseInDense) throws IOException {
		final int READ_AHEAD_LIMIT = maxIntervalsPerRead*100;  // Characters
		boolean success;
		int i, j;
		int fixed, fixedNodes, nPairs, nSolutions, lastRead, lastPeriodicInterval, lastDenseInterval, lastAlignmentInterval, firstJForNextI;
		int tmpType, tmpID, nullStart, nullEnd;
		Set<Map.Entry<String,UnassignedAlignmentPair>> entries;
		Iterator<Map.Entry<String,UnassignedAlignmentPair>> iterator;
		BufferedReader periodicBuffer, denseBuffer, alignmentsBuffer;
		AlignmentPair alignmentPair = new AlignmentPair(null,null,-1,-1,-1,-1,false,0.0,false);
		Node tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Edge queryEdge = new Edge();
		Edge tmpEdge;
		int[] reads;
		int[] periodicIntervalsRead, denseIntervalsRead, alignmentIntervalsRead;
		int[][] fixes;
		UnassignedAlignmentPair[] pairs;
		int[][] periodicIntervals, denseIntervals, alignmentIntervals;
		int[][] tmpMatrix = new int[maxIntervalsPerRead][14];
		
		// Building the sorted list of all unassigned pairs, and of all distinct read IDs
		// that occur in an unassigned pair.
		entries=unassignedAlignmentPairs.entrySet();
		nPairs=entries.size();
		pairs = new UnassignedAlignmentPair[nPairs];
		iterator=entries.iterator();
		i=-1;
		while (iterator.hasNext()) pairs[++i]=iterator.next().getValue();
		if (pairs.length>1) Arrays.sort(pairs);
		reads = new int[nPairs<<1];
		lastRead=getReadIDs(pairs,nPairs,reads);
		
		// Loading the intervals of all reads in $reads$.
		periodicBuffer = new BufferedReader(new FileReader(periodicSubstringsFile),IO.BUFFER_SIZE);
		periodicBuffer.mark(READ_AHEAD_LIMIT);
		lastPeriodicInterval=PeriodicSubstringInterval.countPeriodicIntervals(reads,lastRead,periodicBuffer,READ_AHEAD_LIMIT,tmpMatrix,tmpArray)-1;
		periodicBuffer.close();
		periodicIntervalsRead = new int[lastPeriodicInterval+1];
		periodicIntervals = new int[lastPeriodicInterval+1][11];
		periodicBuffer = new BufferedReader(new FileReader(periodicSubstringsFile),IO.BUFFER_SIZE);
		periodicBuffer.mark(READ_AHEAD_LIMIT);
		PeriodicSubstringInterval.loadPeriodicIntervals(reads,lastRead,periodicBuffer,READ_AHEAD_LIMIT,periodicIntervalsRead,periodicIntervals,tmpMatrix,tmpArray);
		periodicBuffer.close();
		
		denseBuffer = new BufferedReader(new FileReader(denseSubstringsFile),IO.BUFFER_SIZE);
		denseBuffer.mark(READ_AHEAD_LIMIT);
		lastDenseInterval=DenseSubstring.countDenseIntervals(reads,lastRead,denseBuffer,READ_AHEAD_LIMIT,tmpMatrix,tmpArray);
		denseBuffer.close();
		denseIntervalsRead = new int[lastDenseInterval+1];
		denseIntervals = new int[lastDenseInterval+1][14];
		denseBuffer = new BufferedReader(new FileReader(denseSubstringsFile),IO.BUFFER_SIZE);
		denseBuffer.mark(READ_AHEAD_LIMIT);
		DenseSubstring.loadDenseIntervals(reads,lastRead,denseBuffer,READ_AHEAD_LIMIT,denseIntervalsRead,denseIntervals,tmpMatrix,tmpArray);
		denseBuffer.close();
		
		alignmentsBuffer = new BufferedReader(new FileReader(alignmentsFile),IO.BUFFER_SIZE);
		alignmentsBuffer.mark(READ_AHEAD_LIMIT);
		lastAlignmentInterval=AlignmentInterval.countAlignmentIntervals(reads,lastRead,alignmentsBuffer,READ_AHEAD_LIMIT,tmpMatrix,tmpArray);
		alignmentsBuffer.close();
		alignmentIntervalsRead = new int[lastAlignmentInterval+1];
		alignmentIntervals = new int[lastAlignmentInterval+1][8];
		alignmentsBuffer = new BufferedReader(new FileReader(alignmentsFile),IO.BUFFER_SIZE);
		alignmentsBuffer.mark(READ_AHEAD_LIMIT);
		AlignmentInterval.loadAlignmentIntervals(reads,lastRead,alignmentsBuffer,READ_AHEAD_LIMIT,alignmentIntervalsRead,alignmentIntervals,tmpMatrix,tmpArray);
		alignmentsBuffer.close();
		
		// Trying to fix pairs with periodic intervals.
		// Row $i$ of $fixes$ has the following format: 
		//
		// read1,type1,id1,read2,type2,id2,
		// start1,end1,isLeftMaximal1,isRightMaximal1,isWeak1,isContained1,
		// start2,end2,isLeftMaximal2,isRightMaximal2,isWeak2,isContained2,
		// period1,period2,
		// minPrefixLength1,minSuffixLength1,
		// minPrefixLength2,minSuffixLength2,
		// hasLongPeriod1,hasLongPeriod2,
		// firstMaximalStart1,lastMaximalEnd1,
		// firstMaximalStart2,lastMaximalEnd2,
	    // nAssignedAlignments1,nAssignedAlignments2,
		// avgCoverage1,avgCoverage2.
		//
		// type1=-1: pair $i$ never assigned; -2: assigned multiple times; >=0: assigned
		// exactly once.
		fixedNodes=nodes.size();
		fixes = new int[nPairs][34];
		Math.set(fixes,-1);
		for (i=0; i<nPairs; i++) {
			fixes[i][0]=pairs[i].read1;
			fixes[i][3]=pairs[i].read2;
		}
		i=0; j=0; firstJForNextI=-1;
		while (i<nPairs && j<=lastPeriodicInterval) {
			nullStart=pairs[i].alignmentStart1;
			nullEnd=pairs[i].alignmentEnd1;
			if (i<nPairs-1 && pairs[i+1].read1==pairs[i].read1 && firstJForNextI==-1 && periodicIntervals[j][2]>=pairs[i+1].alignmentStart1) firstJForNextI=j;
			if ( periodicIntervalsRead[j]>pairs[i].read1 || 
			     (periodicIntervalsRead[j]==pairs[i].read1 && periodicIntervals[j][1]>pairs[i].alignmentEnd1) 
			   ) {
				i++;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			if ( periodicIntervalsRead[j]<pairs[i].read1 || periodicIntervals[j][2]<pairs[i].alignmentStart1 ||
			     ((discardContainedInShortPeriod||discardPeriodicContainedInPeriodic) && (periodicIntervals[j][5]&(1<<Constants.INTERVAL_PERIODIC))!=0) ||
			     ((discardContainedInLongPeriod||discardPeriodicContainedInPeriodic) && (periodicIntervals[j][5]&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0)
			   ) {
				j++;
				continue;
			}
			if ( !Intervals.isApproximatelyContained(nullStart,nullEnd,periodicIntervals[j][1],periodicIntervals[j][2]) &&
				 !Intervals.areApproximatelyIdentical(nullStart,nullEnd,periodicIntervals[j][1],periodicIntervals[j][2])
			   ) {
   				j++;
   				continue;
			}
			tmpType=fixes[i][1]; tmpID=fixes[i][2];
			fixes[i][1]=Constants.INTERVAL_PERIODIC;
			fixes[i][2]=periodicIntervals[j][0];
			fixes[i][6]=periodicIntervals[j][1];
			fixes[i][7]=periodicIntervals[j][2];
			fixes[i][8]=periodicIntervals[j][3];
			fixes[i][9]=periodicIntervals[j][4];
			fixes[i][10]=0;
			fixes[i][11]=periodicIntervals[j][5];
			fixes[i][18]=periodicIntervals[j][6];
			fixes[i][20]=-1; 
			fixes[i][21]=-1;
			fixes[i][24]=periodicIntervals[j][7];
			fixes[i][26]=-1;
			fixes[i][27]=-1;
			fixes[i][30]=periodicIntervals[j][9];
			fixes[i][32]=periodicIntervals[j][10];
			nSolutions=fixUnassignedAlignments_impl(pairs[i],fixes[i],nodes,periodicIntervalsRead,periodicIntervals,lastPeriodicInterval,denseIntervalsRead,denseIntervals,lastDenseInterval,alignmentIntervalsRead,alignmentIntervals,lastAlignmentInterval,discardContainedInShortPeriod,discardContainedInLongPeriod,discardPeriodicContainedInPeriodic,discardDenseInDense);
			if (nSolutions==-2 || (nSolutions>=0 && tmpType!=-1)) {
				fixes[i][1]=-2;
				i++;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			else if (nSolutions==-1) {
				fixes[i][1]=tmpType;
				fixes[i][2]=tmpID;
			}
			j++;
		}
		
		// Trying to fix pairs with dense intervals
		i=0; j=0; firstJForNextI=-1;
		while (i<nPairs && j<=lastDenseInterval) {
			if (fixes[i][1]==-2) {
				i++;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			nullStart=pairs[i].alignmentStart1;
			nullEnd=pairs[i].alignmentEnd1;
			if (i<nPairs-1 && pairs[i+1].read1==pairs[i].read1 && firstJForNextI==-1 && denseIntervals[j][3]>=pairs[i+1].alignmentStart1) firstJForNextI=j;
			if ( denseIntervalsRead[j]>pairs[i].read1 || 
			     (denseIntervalsRead[j]==pairs[i].read1 && denseIntervals[j][2]>pairs[i].alignmentEnd1) 
			   ) {
				i++;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			if ( denseIntervalsRead[j]<pairs[i].read1 || denseIntervals[j][3]<pairs[i].alignmentStart1 ||
		         (discardContainedInShortPeriod && (denseIntervals[j][7]&(1<<Constants.INTERVAL_PERIODIC))!=0) ||
		         (discardContainedInLongPeriod && (denseIntervals[j][7]&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0) ||
			     (discardDenseInDense && (denseIntervals[j][7]&(1<<Constants.INTERVAL_DENSE_SUBSTRING))!=0)
			   ) {
				j++;
				continue;
			}
			if ( !Intervals.isApproximatelyContained(nullStart,nullEnd,denseIntervals[j][2],denseIntervals[j][3]) &&
				 !Intervals.areApproximatelyIdentical(nullStart,nullEnd,denseIntervals[j][2],denseIntervals[j][3])
			   ) {
   				j++;
   				continue;
			}
			tmpType=fixes[i][1]; tmpID=fixes[i][2];
			fixes[i][1]=denseIntervals[j][1];
			fixes[i][2]=denseIntervals[j][0];
			fixes[i][6]=denseIntervals[j][2];
			fixes[i][7]=denseIntervals[j][3];
			fixes[i][8]=denseIntervals[j][4];
			fixes[i][9]=denseIntervals[j][5];
			fixes[i][10]=denseIntervals[j][6];
			fixes[i][11]=denseIntervals[j][7];
			fixes[i][18]=-1;
			fixes[i][20]=denseIntervals[j][8];
			fixes[i][21]=denseIntervals[j][9];
			fixes[i][24]=-1;
			fixes[i][26]=denseIntervals[j][10];
			fixes[i][27]=denseIntervals[j][11];
			fixes[i][30]=denseIntervals[j][12];
			fixes[i][32]=denseIntervals[j][13];
			nSolutions=fixUnassignedAlignments_impl(pairs[i],fixes[i],nodes,periodicIntervalsRead,periodicIntervals,lastPeriodicInterval,denseIntervalsRead,denseIntervals,lastDenseInterval,alignmentIntervalsRead,alignmentIntervals,lastAlignmentInterval,discardContainedInShortPeriod,discardContainedInLongPeriod,discardPeriodicContainedInPeriodic,discardDenseInDense);
			if (nSolutions==-2 || (nSolutions>=0 && tmpType!=-1)) {
				fixes[i][1]=-2;
				i++;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			else if (nSolutions==-1) {
				fixes[i][1]=tmpType;
				fixes[i][2]=tmpID;
			}
			j++;
		}
		
		// Trying to fix pairs with alignment intervals
		i=0; j=0; firstJForNextI=-1;
		while (i<nPairs && j<=lastAlignmentInterval) {
			if (fixes[i][1]==-2) {
				i++;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			nullStart=pairs[i].alignmentStart1;
			nullEnd=pairs[i].alignmentEnd1;
			if (i<nPairs-1 && pairs[i+1].read1==pairs[i].read1 && firstJForNextI==-1 && alignmentIntervals[j][2]>=pairs[i+1].alignmentStart1) firstJForNextI=j;
			if ( alignmentIntervalsRead[j]>pairs[i].read1 || 
			     (alignmentIntervalsRead[j]==pairs[i].read1 && alignmentIntervals[j][1]>pairs[i].alignmentEnd1) 
			   ) {
				i++;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			if ( alignmentIntervalsRead[j]<pairs[i].read1 || alignmentIntervals[j][2]<pairs[i].alignmentStart1 ||
		         (discardContainedInShortPeriod && (alignmentIntervals[j][5]&(1<<Constants.INTERVAL_PERIODIC))!=0) ||
		         (discardContainedInLongPeriod && (alignmentIntervals[j][5]&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0)
			   ) {
				j++;
				continue;
			}
			if ( !Intervals.isApproximatelyContained(nullStart,nullEnd,alignmentIntervals[j][1],alignmentIntervals[j][2]) &&
				 !Intervals.areApproximatelyIdentical(nullStart,nullEnd,alignmentIntervals[j][1],alignmentIntervals[j][2])
			   ) {
   				j++;
   				continue;
			}
			tmpType=fixes[i][1]; tmpID=fixes[i][2];
			fixes[i][1]=Constants.INTERVAL_ALIGNMENT;
			fixes[i][2]=alignmentIntervals[j][0];
			fixes[i][6]=alignmentIntervals[j][1];
			fixes[i][7]=alignmentIntervals[j][2];
			fixes[i][8]=alignmentIntervals[j][3];
			fixes[i][9]=alignmentIntervals[j][4];
			fixes[i][10]=0;
			fixes[i][11]=alignmentIntervals[j][5];
			fixes[i][18]=-1; 
			fixes[i][20]=-1; 
			fixes[i][21]=-1;
			fixes[i][24]=-1;
			fixes[i][26]=-1;
			fixes[i][27]=-1;
			fixes[i][30]=alignmentIntervals[j][6];
			fixes[i][32]=alignmentIntervals[j][7];
			nSolutions=fixUnassignedAlignments_impl(pairs[i],fixes[i],nodes,periodicIntervalsRead,periodicIntervals,lastPeriodicInterval,denseIntervalsRead,denseIntervals,lastDenseInterval,alignmentIntervalsRead,alignmentIntervals,lastAlignmentInterval,discardContainedInShortPeriod,discardContainedInLongPeriod,discardPeriodicContainedInPeriodic,discardDenseInDense);
			if (nSolutions==-2 || (nSolutions>=0 && tmpType!=-1)) {
				fixes[i][1]=-2;
				i++;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			else if (nSolutions==-1) {
				fixes[i][1]=tmpType;
				fixes[i][2]=tmpID;
			}
			j++;
		}
		
		// Applying unique fixes
		fixed=0;
		for (i=0; i<nPairs; i++) {
			if (fixes[i][1]<0) continue;
			success=pairs[i].toAlignmentPair(fixes[i],nodes,alignmentPair,tmpNode); 
			if (!success) {
				// Should not fail, since it has already been invoked by
				// $fixUnassignedAlignments_impl()$.
				System.err.println("ERROR: Failed to create an alignment pair");
				System.err.println("fixes["+i+"]="+fixes[i][0]+","+fixes[i][1]+","+fixes[i][2]+",  "+fixes[i][3]+","+fixes[i][4]+","+fixes[i][5]);
				System.exit(1);
			}
			// The instruction below should not fail, since it has already been invoked by
			// $fixUnassignedAlignments_impl()$.
			queryEdge.set(alignmentPair);
			tmpEdge=edges.get(queryEdge);
			if (tmpEdge==null) {
				tmpEdge = new Edge(alignmentPair);
				edges.put(tmpEdge,tmpEdge);
			}
			else {
				tmpEdge.addTypes(queryEdge);
				tmpEdge.addDiffs(queryEdge);
				tmpEdge.addOrientation(queryEdge);
			}
			fixed++;
		}
		
		System.err.println("Fixed "+fixed+" alignments unassigned in both reads ("+IO.getPercent(fixed,nPairs)+"%). New nodes created: "+(nodes.size()-fixedNodes));
	}
	
	
	/**
	 * Stores in $out[0..X]$ the sorted list of distinct read IDs of all elements in 
	 * $pairs$, where $X$ is returned in output.
	 */
	private static final int getReadIDs(UnassignedAlignmentPair[] pairs, int nPairs, int[] out) {
		int i, j, k;
		int read, previousRead, lastRead;
		
		k=-1; previousRead=-1;
		for (i=0; i<nPairs; i++) {
			read=pairs[i].read1;
			if (read==previousRead) continue;
			out[++k]=read;
			previousRead=read;
		}
		lastRead=k;
		for (i=0; i<nPairs; i++) {
			read=pairs[i].read2;
			j=Arrays.binarySearch(out,0,lastRead+1,read);
			if (j<0) out[++k]=read;
		}
		lastRead=k;
		if (lastRead>0) Arrays.sort(out,0,lastRead+1);
		k=-1; previousRead=-1;
		for (i=0; i<=lastRead; i++) {
			read=out[i];
			if (read==previousRead) continue;
			out[++k]=read;
			previousRead=read;
		}
		return k;
	}
	
	
	/**
	 * Let $out$ be an array with the following format: 
	 *
	 * read1,type1,id1,read2,type2,id2,
	 * start1,end1,isLeftMaximal1,isRightMaximal1,isWeak1,isContained1,
	 * start2,end2,isLeftMaximal2,isRightMaximal2,isWeak2,isContained2,
	 * period1,period2,
	 * minPrefixLength1,minSuffixLength1,
	 * minPrefixLength2,minSuffixLength2,
	 * hasLongPeriod1,hasLongPeriod2,
	 * firstMaximalStart1,lastMaximalEnd1,
	 * firstMaximalStart2,lastMaximalEnd2,
	 * nAssignedAlignments1,nAssignedAlignments2,
	 * avgCoverage1,avgCoverage2.
	 *
	 * where $*1$ and $read2$ are assumed to have already been set by the caller,
	 * and where $id1,read2$ are assumed to be initialized to -1. The procedure writes in 
	 * $*2$ the only interval in $pair.read2$ that satisfies $IntervalGraphStep2.inStep2$.
	 * If multiple such intervals exist, the content of $*2$ is undefined.
	 *
	 * @return -1: $pair$ not fixed; -2: multiple possible fixes; 0: exactly one fix.
	 */
	private static final int fixUnassignedAlignments_impl(UnassignedAlignmentPair pair, int[] out, HashMap<Node,Node> nodes, int[] periodicIntervalsRead, int[][] periodicIntervals, int lastPeriodicInterval, int[] denseIntervalsRead, int[][] denseIntervals, int lastDenseInterval, int[] alignmentIntervalsRead, int[][] alignmentIntervals, int lastAlignmentInterval, boolean discardContainedInShortPeriod, boolean discardContainedInLongPeriod, boolean discardPeriodicContainedInPeriodic, boolean discardDenseInDense) {
		final int MIN_ALIGNMENT_LENGTH = Alignments.minAlignmentLength-IO.quantum;  // Arbitrary
		final int READ2 = pair.read2;
		final int START2 = pair.alignmentStart2;
		final int END2 = pair.alignmentEnd2;
		boolean fixed;
		int i;
		int first, tmpType, tmpID, edgeType;
		Node tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		AlignmentPair alignmentPair = new AlignmentPair(null,null,-1,-1,-1,-1,false,0.0,false);
		
		fixed=false;
		
		// Trying to fix $pair$ with periodic intervals
		i=Arrays.binarySearch(periodicIntervalsRead,0,lastPeriodicInterval+1,READ2);
		if (i>=0) {
			first=i;
			i--;
			while (i>=0 && periodicIntervalsRead[i]==READ2) {
				first=i;
				i--;
			}
			i=first;
			while (i<=lastPeriodicInterval && periodicIntervalsRead[i]==READ2) {
				if (periodicIntervals[i][1]>END2) break;
				if ( periodicIntervals[i][2]<START2 ||
				     ((discardContainedInShortPeriod||discardPeriodicContainedInPeriodic) && (periodicIntervals[i][5]&(1<<Constants.INTERVAL_PERIODIC))!=0) ||
				     ((discardContainedInLongPeriod||discardPeriodicContainedInPeriodic) && (periodicIntervals[i][5]&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0)
				   ) {
					i++;
					continue;
				}
				if ( !Intervals.isApproximatelyContained(START2,END2,periodicIntervals[i][1],periodicIntervals[i][2]) &&
					 !Intervals.areApproximatelyIdentical(START2,END2,periodicIntervals[i][1],periodicIntervals[i][2])
				   ) {
	   				i++;
	   				continue;
				}
				tmpType=out[4]; tmpID=out[5];
				out[4]=Constants.INTERVAL_PERIODIC;
				out[5]=periodicIntervals[i][0];
				out[12]=periodicIntervals[i][1];
				out[13]=periodicIntervals[i][2];
				out[14]=periodicIntervals[i][3];
				out[15]=periodicIntervals[i][4];
				out[16]=0;
				out[17]=periodicIntervals[i][5];
				out[19]=periodicIntervals[i][6];
				out[22]=-1; 
				out[23]=-1;
				out[25]=periodicIntervals[i][7];
				out[28]=-1;
				out[29]=-1;
				out[31]=periodicIntervals[i][9];
				out[33]=periodicIntervals[i][10];
				if (pair.toAlignmentPair(out,nodes,alignmentPair,tmpNode)) {
					alignmentPair.trim();
					edgeType=IntervalGraphStep2.inStep2(alignmentPair);
					if (edgeType>=0 && edgeType!=Constants.SHARED_SUBSTRING && edgeType!=Constants.INSERTION_BOTH && alignmentPair.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && alignmentPair.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
						if (fixed) return -2;
						fixed=true;
					}
					else {
						out[4]=tmpType;
						out[5]=tmpID;
					}
				}
				else {
					out[4]=tmpType;
					out[5]=tmpID;
				}
				i++;
			}
		}
		
		// Trying to fix $pair$ with dense intervals
		i=Arrays.binarySearch(denseIntervalsRead,0,lastDenseInterval+1,READ2);
		if (i>=0) {
			first=i;
			i--;
			while (i>=0 && denseIntervalsRead[i]==READ2) {
				first=i;
				i--;
			}
			i=first;
			while (i<=lastDenseInterval && denseIntervalsRead[i]==READ2) {
				if (denseIntervals[i][1]>END2) break;
				if ( denseIntervals[i][2]<START2 ||
			         (discardContainedInShortPeriod && (denseIntervals[i][7]&(1<<Constants.INTERVAL_PERIODIC))!=0) ||
			         (discardContainedInLongPeriod && (denseIntervals[i][7]&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0) ||
				     (discardDenseInDense && (denseIntervals[i][7]&(1<<Constants.INTERVAL_DENSE_SUBSTRING))!=0)
				   ) {
					i++;
					continue;
				}
				if ( !Intervals.isApproximatelyContained(START2,END2,denseIntervals[i][2],denseIntervals[i][3]) &&
					 !Intervals.areApproximatelyIdentical(START2,END2,denseIntervals[i][2],denseIntervals[i][3])
				   ) {
	   				i++;
	   				continue;
				}
				tmpType=out[4]; tmpID=out[5];
				out[4]=denseIntervals[i][1];
				out[5]=denseIntervals[i][0];
				out[12]=denseIntervals[i][2];
				out[13]=denseIntervals[i][3];
				out[14]=denseIntervals[i][4];
				out[15]=denseIntervals[i][5];
				out[16]=denseIntervals[i][6];
				out[17]=denseIntervals[i][7];
				out[19]=-1;
				out[22]=denseIntervals[i][8];
				out[23]=denseIntervals[i][9];
				out[25]=-1;
				out[28]=denseIntervals[i][10];
				out[29]=denseIntervals[i][11];
				out[31]=denseIntervals[i][12];
				out[33]=denseIntervals[i][13];
				if (pair.toAlignmentPair(out,nodes,alignmentPair,tmpNode)) {
					alignmentPair.trim();
					edgeType=IntervalGraphStep2.inStep2(alignmentPair);
					if (edgeType>=0 && edgeType!=Constants.SHARED_SUBSTRING && edgeType!=Constants.INSERTION_BOTH && alignmentPair.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && alignmentPair.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
						if (fixed) return -2;
						fixed=true;
					}
					else {
						out[4]=tmpType;
						out[5]=tmpID;
					}
				}
				else {
					out[4]=tmpType;
					out[5]=tmpID;
				}
				i++;
			}
		}
		
		// Trying to fix $pair$ with alignment intervals
		i=Arrays.binarySearch(alignmentIntervalsRead,0,lastAlignmentInterval+1,READ2);
		if (i>=0) {
			first=i;
			i--;
			while (i>=0 && alignmentIntervalsRead[i]==READ2) {
				first=i;
				i--;
			}
			i=first;
			while (i<=lastAlignmentInterval && alignmentIntervalsRead[i]==READ2) {
				if (alignmentIntervals[i][1]>END2) break;
				if ( alignmentIntervals[i][2]<START2 ||
   			         (discardContainedInShortPeriod && (alignmentIntervals[i][5]&(1<<Constants.INTERVAL_PERIODIC))!=0) ||
   			         (discardContainedInLongPeriod && (alignmentIntervals[i][5]&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0)
				   ) {
					i++;
					continue;
				}
				if ( !Intervals.isApproximatelyContained(START2,END2,alignmentIntervals[i][1],alignmentIntervals[i][2]) &&
					 !Intervals.areApproximatelyIdentical(START2,END2,alignmentIntervals[i][1],alignmentIntervals[i][2])
				   ) {
	   				i++;
	   				continue;
				}
				tmpType=out[4]; tmpID=out[5];
				out[4]=Constants.INTERVAL_ALIGNMENT;
				out[5]=alignmentIntervals[i][0];
				out[12]=alignmentIntervals[i][1];
				out[13]=alignmentIntervals[i][2];
				out[14]=alignmentIntervals[i][3];
				out[15]=alignmentIntervals[i][4];
				out[16]=0;
				out[17]=alignmentIntervals[i][5];
				out[19]=-1; 
				out[22]=-1; 
				out[23]=-1;
				out[25]=-1;
				out[28]=-1;
				out[29]=-1;
				out[31]=alignmentIntervals[i][6];
				out[33]=alignmentIntervals[i][7];
				if (pair.toAlignmentPair(out,nodes,alignmentPair,tmpNode)) {
					alignmentPair.trim();
					edgeType=IntervalGraphStep2.inStep2(alignmentPair);
					if (edgeType>=0 && edgeType!=Constants.SHARED_SUBSTRING && edgeType!=Constants.INSERTION_BOTH && alignmentPair.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && alignmentPair.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
						if (fixed) return -2;
						fixed=true;
					}
					else {
						out[4]=tmpType;
						out[5]=tmpID;
					}
				}
				else {
					out[4]=tmpType;
					out[5]=tmpID;
				}
				i++;
			}
		}

		return fixed?0:-1;
	}
	
	
	public static class UnassignedAlignmentPair implements Comparable {
		// read1 is the one with smaller ID
		public int read1, read2, alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;  // In the forward orientation of both reads
		public boolean orientation;
		public double diffs;  // Normalized

		
		public UnassignedAlignmentPair(int r1, int r2, int s1, int e1, int s2, int e2, boolean o, int d) {
			read1=r1;
			read2=r2;
			alignmentStart1=s1;
			alignmentEnd1=e1;
			alignmentStart2=s2;
			alignmentEnd2=e2;
			orientation=o;
			diffs=2.0*d/(alignmentEnd1-alignmentStart1+alignmentEnd2-alignmentStart2+2);
		}
		
		
		public int compareTo(Object other) {
			UnassignedAlignmentPair otherPair = (UnassignedAlignmentPair)other;
			if (read1<otherPair.read1) return -1;
			if (read1>otherPair.read1) return 1;
			if (alignmentStart1<otherPair.alignmentStart1) return -1;
			if (alignmentStart1>otherPair.alignmentStart1) return 1;
			return 0;
		}
		
		
		/**
		 * Loads in $out$ the alignment pair that correspond to the merge of this 
		 * unassigned alignment and $intervals$.
		 *
		 * Remark: the procedure might add nodes to $nodes$.
		 * 
		 * @param intervals format:
		 * read1,type1,id1,read2,type2,id2, 
		 * start1,end1,isLeftMaximal1,isRightMaximal1,isWeak1,isContained1,
		 * start2,end2,isLeftMaximal2,isRightMaximal2,isWeak2,isContained2,
		 * period1,period2;
		 * minPrefixLength1,minSuffixLength1;
		 * minPrefixLength2,minSuffixLength2;
		 * hasLongPeriod1,hasLongPeriod2,
		 * firstMaximalStart1,lastMaximalEnd1,
		 * firstMaximalStart2,lastMaximalEnd2,
		 * nAssignedAlignments1,nAssignedAlignments2,
		 * avgCoverage1,avgCoverage2.
		 * @param tmpNode temporary space;
		 * @return false iff no alignment pair can be built.
		 */
		public boolean toAlignmentPair(int[] intervals, HashMap<Node,Node> nodes, AlignmentPair out, Node tmpNode) {
			Node tmpNodePrime;
			
			out.alignmentStart1=alignmentStart1;
			out.alignmentEnd1=alignmentEnd1;
			out.alignmentStart2=alignmentStart2;
			out.alignmentEnd2=alignmentEnd2;
			out.orientation=orientation;
			out.diffs=diffs;
			
			tmpNode.read=intervals[0]; 
			tmpNode.type=intervals[1]; 
			tmpNode.id=intervals[2];
			tmpNodePrime=nodes.get(tmpNode);
			if (tmpNodePrime==null) {
				tmpNodePrime = new Node(tmpNode.read,tmpNode.type,tmpNode.id,intervals[6],intervals[7],intervals[8]==1,intervals[9]==1,intervals[10]==1,intervals[11],intervals[20],intervals[21],tmpNode.type==Constants.INTERVAL_PERIODIC?intervals[18]:-1,intervals[24]==1,intervals[26],intervals[27],intervals[30],intervals[32]);
				tmpNodePrime.nodeID=nodes.size();
				nodes.put(tmpNodePrime,tmpNodePrime);
			}
			out.node1=tmpNodePrime;
			
			tmpNode.read=intervals[3]; 
			tmpNode.type=intervals[4]; 
			tmpNode.id=intervals[5];
			tmpNodePrime=nodes.get(tmpNode);
			if (tmpNodePrime==null) {
				tmpNodePrime = new Node(tmpNode.read,tmpNode.type,tmpNode.id,intervals[12],intervals[13],intervals[14]==1,intervals[15]==1,intervals[16]==1,intervals[17],intervals[22],intervals[23],tmpNode.type==Constants.INTERVAL_PERIODIC?intervals[19]:-1,intervals[25]==1,intervals[28],intervals[29],intervals[31],intervals[33]);
				tmpNodePrime.nodeID=nodes.size();
				nodes.put(tmpNodePrime,tmpNodePrime);
			}
			out.node2=tmpNodePrime;
			
			return true;
		}
	}
	
	
	
	
	
	
	
	
	// --------------------------- GENERAL PRINTING PROCEDURES ---------------------------
	
	/**
	 * Prints a file for each kernel family, containing the set of all intervals tagged 
	 * with a kernel that belongs to that family, in the following format:
	 *
	 * readID,start,end,type,F1,...,F_k
	 *
	 * where $type$ is the type of interval, and F_1,...,F_k are all the distinct families
	 * of the node. F_1,...,F_k are printed only if $k>1$. The name of the file contains 
	 * the ID of the family.
	 *
	 * Remark: families that are not flagged in $selectFamily$ are not printed, and they 
	 * are assumed to have both the number of left-maximal and of right-maximal intervals 
	 * not bigger than the coverage. Call "frequent" a kernel with either the number of 
	 * left-maximal or of right-maximal intervals bigger than the coverage. Such families 
	 * do not contain any frequent kernel. If only frequent kernels are output by the
	 * pipeline, then an outputted family is not just a subset of the outputted kernels, 
	 * and some outputted families might contain intervals that belong to no outputted 
	 * kernel.
	 *
	 * @param nFamilies number of distinct IDs in $kernel2family$;
	 * @param step step of family construction.
	 */
	public static final void printKernelFamilies(int nFamilies, int[] kernel2family, boolean[] selectFamily, String outDir, int step) throws IOException {
		final String SUFFIX = "-family"+step+".txt";
		int i, j, f;
		int family, firstFamily, lastKernel, nOutputFiles;
		String allFamilies;
		Node tmpNode;
		boolean[] marked;
		BufferedWriter[] familyFiles;
		
		familyFiles = new BufferedWriter[nFamilies];
		marked = new boolean[nFamilies];
		firstFamily=0; nOutputFiles=IO.MAX_OPEN_FILES;
		while (firstFamily<nFamilies) {
			for (i=0; i<nOutputFiles && firstFamily+i<nFamilies; i++) {
				if (selectFamily[firstFamily+i]) familyFiles[firstFamily+i] = new BufferedWriter(new FileWriter(outDir+"/"+(firstFamily+i)+SUFFIX),IO.BUFFER_SIZE);
			}
			for (i=0; i<nNodes; i++) {
				tmpNode=nodesArray[i];
				lastKernel=tmpNode.lastKernel;
				Math.set(marked,nFamilies-1,false); allFamilies=""; f=0;
				for (j=0; j<=lastKernel; j++) {
					family=kernel2family[tmpNode.kernels[j]];
					if (family<0 || !selectFamily[family]) continue;
					if (!marked[family]) {
						marked[family]=true;
						f++;
						allFamilies+=family+",";
					}
				}
				for (j=0; j<=lastKernel; j++) {
					family=kernel2family[tmpNode.kernels[j]];
					if (family<firstFamily || family>=firstFamily+nOutputFiles) continue;
					if (marked[family]) {
						familyFiles[family].write(tmpNode.read+","+tmpNode.start+","+tmpNode.end+","+tmpNode.type+(f>1?","+allFamilies:"")+"\n");
						marked[family]=false;
					}
				}
			}
			for (i=0; i<nOutputFiles && firstFamily+i<nFamilies; i++) {
				if (familyFiles[firstFamily+i]==null) continue;
				familyFiles[firstFamily+i].close();
				familyFiles[firstFamily+i]=null;
			}
			firstFamily+=nOutputFiles;
		}
		familyFiles=null; marked=null;
	}
	
	
	/**
	 * Prints one file for each component in $components$, which is assumed to contain the
	 * sorted list of all component IDs of nodes with just one component. Every row of the
	 * file has the following format:
	 *
	 * readID,start,end,type
	 *
	 * The name of the file contains the ID of the component.
	 */
	public static final void printComponentsAsFragments(int[] components, String outDir) throws IOException {
		final String SUFFIX = "-componentFragments.txt";
		int i, component, firstComponent, nOutputFiles;
		IntervalGraph.Node node;
		BufferedWriter[] outputFiles;

		outputFiles = new BufferedWriter[components.length];
		nOutputFiles=IO.MAX_OPEN_FILES; firstComponent=0;
		while (firstComponent<components.length) {
			for (i=0; i<nOutputFiles && firstComponent+i<components.length; i++) outputFiles[i] = new BufferedWriter(new FileWriter(outDir+"/"+components[firstComponent+i]+SUFFIX),IO.BUFFER_SIZE);
			for (i=0; i<nNodes; i++) {
				node=nodesArray[i];
				if (node.lastComponent!=0) continue;
				component=Arrays.binarySearch(components,node.components[0]);
				if (component<0) {
					System.err.println("ERROR: component not found?");
					System.exit(1);
				}
				if (component<firstComponent || component>=firstComponent+nOutputFiles) continue;
				outputFiles[component].write(node.read+","+node.start+","+node.end+","+node.type+"\n");
			}
			for (i=0; i<nOutputFiles && firstComponent+i<components.length; i++) {
				if (outputFiles[firstComponent+i]==null) continue;
				outputFiles[firstComponent+i].close();
				outputFiles[firstComponent+i]=null;
			}
			firstComponent+=nOutputFiles;
		}
	}	
	
	
	/**
	 * Writes each row of the input files to potentially multiple files, associated with 
	 * all the components to which the readA interval of the row belongs. If 
	 * $strict=false$, the procedure prints also rows whose readA corresponds to a node 
	 * with more than one component: the motivation for this is that a large cluster might
	 * be mostly frontier and just few of its nodes might have a single component. 
	 * Otherwise, only rows whose readA corresponds to a node with a single component are 
	 * printed, and each row is printed to exactly one file.
	 *
	 * @param printComponent one flag per element in $componentIDs$; components not marked
	 * in this array are discarded; if the array is null, no component is discarded;
	 * @param componentIDs set of distinct connected component IDs; only components in
	 * this set are considered.
	 */
	public static final void splitInputByComponent(int[] componentIDs, boolean[] printComponent, String alignmentsFilePath, String connectionFilePath, String periodicSubstringsFile, String denseSubstringsFile, String alignmentsFile, String outputDir, boolean strict, boolean[] isComponentPeriodic) throws IOException {
		splitInputByComponent_alignmentsConnectionsFiles(componentIDs,printComponent,alignmentsFilePath,connectionFilePath,outputDir,strict,isComponentPeriodic);
		splitInputByComponent_intervalsFile(0,componentIDs,printComponent,periodicSubstringsFile,outputDir,strict,isComponentPeriodic);
		splitInputByComponent_intervalsFile(1,componentIDs,printComponent,denseSubstringsFile,outputDir,strict,isComponentPeriodic);
		splitInputByComponent_intervalsFile(2,componentIDs,printComponent,alignmentsFile,outputDir,strict,isComponentPeriodic);
	}
	
	
	/**
	 * A version of $splitInputByComponent$ specialized for the pair 
	 * (alignments,connections).
	 */
	private static final void splitInputByComponent_alignmentsConnectionsFiles(int[] componentIDs, boolean[] printComponent, String alignmentsFilePath, String connectionFilePath, String outputDir, boolean strict, boolean[] isComponentPeriodic) throws IOException {
		final String ALIGNMENTS_SUFFIX = "-alignments.txt";
		final String CONNECTIONS_SUFFIX = "-connections.txt";
		int i, j;
		int component, firstComponent, nComponents, nOutputFiles;
		String str1, str2, suffix;
		Node tmpNode;
		BufferedReader alignmentsFile, connectionFile;
		BufferedWriter[] alignmentsFiles, connectionFiles;
		
		nComponents=componentIDs.length;
		alignmentsFiles = new BufferedWriter[nComponents];
		connectionFiles = new BufferedWriter[nComponents];
		firstComponent=0; nOutputFiles=(IO.MAX_OPEN_FILES>>1)-2;
		while (firstComponent<nComponents) {
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				if (printComponent!=null && !printComponent[firstComponent+i]) continue;
				if (isComponentPeriodic!=null) suffix=isComponentPeriodic[firstComponent+i]?"-p":"-np";
				else suffix="";
				alignmentsFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+componentIDs[firstComponent+i]+suffix+ALIGNMENTS_SUFFIX),IO.BUFFER_SIZE);
				alignmentsFiles[firstComponent+i].write("\n");
				alignmentsFiles[firstComponent+i].write("\n");
				connectionFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+componentIDs[firstComponent+i]+suffix+CONNECTIONS_SUFFIX),IO.BUFFER_SIZE);
			}
			alignmentsFile = new BufferedReader(new FileReader(alignmentsFilePath),IO.BUFFER_SIZE);
			connectionFile = new BufferedReader(new FileReader(connectionFilePath),IO.BUFFER_SIZE);
			str1=alignmentsFile.readLine();
			str1=alignmentsFile.readLine();  // Skipping the first two lines
			str1=alignmentsFile.readLine();
			str2=connectionFile.readLine();
			i=0;
			while (str1!=null) {
				i++;
				if (str2.length()==0) {
					str1=alignmentsFile.readLine();
					str2=connectionFile.readLine();
					continue;
				}
				Alignments.readAlignmentFile(str1);
				Factorize.readIntervalConnectionFile(str2);
				tmpNode=nodes.get(new Node(Alignments.readA-1,Factorize.type,Factorize.id,Factorize.start,Factorize.end,Factorize.isLeftMaximal,Factorize.isRightMaximal,Factorize.isWeak,Factorize.isContained,Factorize.minPrefixLength,Factorize.minSuffixLength,Factorize.period,Factorize.hasLongPeriod,Factorize.firstMaximalStart,Factorize.lastMaximalEnd,Factorize.nAssignedAlignments,Factorize.avgCoverage));
				if (tmpNode==null) {
					System.err.println("ERROR: node not found in the hashmap");
					System.exit(1);
				}
				if ( tmpNode.lastComponent==-1 || (strict && tmpNode.lastComponent>0) ) {
					// Remark: a node can have no component, if it is a singleton
					// component with no neighbor.
					if (i%1000000==0) System.err.println("Processed alignment "+(i/1000000)+"K");
					str1=alignmentsFile.readLine();
					str2=connectionFile.readLine();
					continue;
				}
				for (i=0; i<=tmpNode.lastComponent; i++) {
					component=tmpNode.components[i];
					j=Arrays.binarySearch(componentIDs,component);
					if (j>=firstComponent && j<firstComponent+nOutputFiles && (printComponent==null||printComponent[j])) {
						alignmentsFiles[j].write(str1+"\n");
						connectionFiles[j].write(str2+"\n");
					}
				}
				// Displaying progress
				if (i%1000000==0) System.err.println("Processed alignment "+(i/1000000)+"K");
				str1=alignmentsFile.readLine();
				str2=connectionFile.readLine();
			}
			alignmentsFile.close(); connectionFile.close();
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				if (alignmentsFiles[firstComponent+i]!=null) {
					alignmentsFiles[firstComponent+i].close();
					alignmentsFiles[firstComponent+i]=null;
				}
				if (connectionFiles[firstComponent+i]!=null) {
					connectionFiles[firstComponent+i].close();
					connectionFiles[firstComponent+i]=null;
				}
			}
			firstComponent+=nOutputFiles;
		}
	}
	
	
	/**
	 * A version of $splitInputByComponent$ specialized for interval files.
	 *
	 * @type type of interval: 0=periodic, 1=dense, 2=alignment.
	 */
	private static final void splitInputByComponent_intervalsFile(int type, int[] componentIDs, boolean[] printComponent, String inputPath, String outputDir, boolean strict, boolean[] isComponentPeriodic) throws IOException {
		final String INTERVALS_PERIODIC_SUFFIX = "-intervals-periodic.txt";
		final String INTERVALS_DENSE_SUFFIX = "-intervals-dense.txt";
		final String INTERVALS_ALIGNMENT_SUFFIX = "-intervals-alignment.txt";
		int i, j;
		int component, firstComponent, nComponents, nOutputFiles;
		String str, suffix, outputSuffix;
		Node tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Node tmpNodePrime;
		BufferedReader inputFile;
		BufferedWriter[] outputFiles;
		
		switch (type) {
			case 0: outputSuffix=INTERVALS_PERIODIC_SUFFIX; break;
			case 1: outputSuffix=INTERVALS_DENSE_SUFFIX; break;
			case 2: outputSuffix=INTERVALS_ALIGNMENT_SUFFIX; break;
			default: outputSuffix="";
		}
		nComponents=componentIDs.length;
		outputFiles = new BufferedWriter[nComponents];
		firstComponent=0; nOutputFiles=IO.MAX_OPEN_FILES-1;
		while (firstComponent<nComponents) {
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				if (printComponent!=null && !printComponent[firstComponent+i]) continue;
				if (isComponentPeriodic!=null) suffix=isComponentPeriodic[firstComponent+i]?"-p":"-np";
				else suffix="";
				outputFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+componentIDs[firstComponent+i]+suffix+outputSuffix),IO.BUFFER_SIZE);
			}
			inputFile = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
			str=inputFile.readLine();
			i=0;
			while (str!=null) {
				i++;
				j=str.indexOf(",");
				tmpNode.read=Integer.parseInt(str.substring(0,j));
				if (type==0) {
					PeriodicSubstringInterval.readPeriodicSubstringsFile(str,j+1,tmpArray,0);
					tmpNode.type=Constants.INTERVAL_PERIODIC;
					tmpNode.id=tmpArray[0];
				}
				else if (type==1) {
					DenseSubstring.readDenseSubstringsFile(str,j+1,tmpArray,0);
					tmpNode.type=DenseSubstring.getType(tmpArray);
					tmpNode.id=tmpArray[0];
				}
				else {
					AlignmentInterval.readAlignmentsFile(str,j+1,tmpArray,0);
					tmpNode.type=Constants.INTERVAL_ALIGNMENT;
					tmpNode.id=tmpArray[0];
				}
				tmpNodePrime=nodes.get(tmpNode);
				if (tmpNodePrime==null) {
					// It is possible for an interval not to be assigned any alignment and
					// thus not to appear in the hashmap.
					if (i%10000==0) System.err.println("Processed interval "+(i/10000)+"K");
					str=inputFile.readLine();
					continue;
				}
				if ( tmpNodePrime.lastComponent==-1 || (strict && tmpNodePrime.lastComponent>0) ) {
					// Remark: a node can have no component, if it is a singleton
					// component with no neighbor.
					if (i%10000==0) System.err.println("Processed interval "+(i/10000)+"K");
					str=inputFile.readLine();
					continue;
				}
				for (i=0; i<=tmpNodePrime.lastComponent; i++) {
					component=tmpNodePrime.components[i];
					j=Arrays.binarySearch(componentIDs,component);
					if (j>=firstComponent && j<firstComponent+nOutputFiles && (printComponent==null||printComponent[j])) outputFiles[j].write(str+"\n");
				}
				// Displaying progress
				if (i%10000==0) System.err.println("Processed interval "+(i/10000)+"K");
				str=inputFile.readLine();
			}
			inputFile.close();
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				if (outputFiles[firstComponent+i]!=null) {
					outputFiles[firstComponent+i].close();
					outputFiles[firstComponent+i]=null;
				}
			}
			firstComponent+=nOutputFiles;
		}
	}
	
	
	/**
	 * Prints the graph induced by ON edges only.
	 *
	 * @param printedComponents if not null, the component field of each node in the DOT 
	 * file is set to -1 if it does not belong to $printedComponents$, which is assumed to
	 * be sorted;
	 * @param clusterField used only if not null.
	 */
	public static final void toDot(String path, int[] printedComponents, int[] clusterField, double edgeThreshold) throws IOException {
		int i, j;
		BufferedWriter bw = new BufferedWriter(new FileWriter(path),IO.BUFFER_SIZE);
		
		// Printing
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=-1;
		}
		bw.write("graph G {\n");
		for (i=0; i<nNodes; i++) bw.write(nodesArray[i].toDot(printedComponents,clusterField,i));
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (!neighbors[i][j].on || neighbors[i][j].printed==1) continue;
				neighbors[i][j].printed=1;
				bw.write(neighbors[i][j].toDot(edgeThreshold));
			}
		}
		bw.write("}\n");
		bw.close();
		
		// Cleaning $printed$ marks
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=-1;
		}
	}
	
	
	/**
	 * Fills the global variable $nodeStats$
	 */
	public static final void getNodeStats() {
		int i, j;
		
		Math.set(nodeStats,nodeStats.length-1,0);
		for (i=0; i<nNodes; i++) {
			switch (nodesArray[i].type) {
				case Constants.INTERVAL_ALIGNMENT: nodeStats[0]++; break;
				case Constants.INTERVAL_DENSE_PREFIX: nodeStats[1]++; break;
				case Constants.INTERVAL_DENSE_SUFFIX: nodeStats[1]++; break;
				case Constants.INTERVAL_DENSE_PREFIXSUFFIX: nodeStats[2]++; break;
				case Constants.INTERVAL_DENSE_SUBSTRING: nodeStats[3]++; break;
				case Constants.INTERVAL_DENSE_SINGLEDELETION: nodeStats[4]++; break;
				case Constants.INTERVAL_PERIODIC: nodeStats[5]++; break;
			}
		}
	}
	
	
	public static final void printNodeStats() {
		System.out.println("Node types: ");
		System.out.println(nodeStats[0]+" alignment ("+IO.format((100*(double)nodeStats[0])/nNodes)+"%)");
		System.out.println(nodeStats[1]+" prefix/suffix ("+IO.format((100*(double)nodeStats[1])/nNodes)+"%)");
		System.out.println(nodeStats[2]+" prefix+suffix ("+IO.format((100*(double)nodeStats[2])/nNodes)+"%)");
		System.out.println(nodeStats[3]+" substring ("+IO.format((100*(double)nodeStats[3])/nNodes)+"%)");
		System.out.println(nodeStats[4]+" single deletion ("+IO.format((100*(double)nodeStats[4])/nNodes)+"%)");
		System.out.println(nodeStats[5]+" periodic ("+IO.format((100*(double)nodeStats[5])/nNodes)+"%)");
	}
	
	
	/**
	 * Fills the global variable $edgeStats$.
	 *
	 * Remark: fractions do not necessarily sum up to one, since an edge can have multiple 
	 * types (especially periodic-periodic edges).
	 *
	 * @return the total number of edges in the graph.
	 */
	public static final int getEdgeStats() {
		boolean periodicI, periodicJ;
		int i, j, k, nEdges;
		
		Math.set(edgeStats,edgeStats.length-1,0);
		nEdges=0;
		for (i=0; i<nNodes; i++) {
			nEdges+=nNeighbors[i];
			periodicI=nodesArray[i].type==Constants.INTERVAL_PERIODIC;
			for (j=0; j<nNeighbors[i]; j++) {
				periodicJ=nodesArray[neighbors[i][j].getTo(i)].type==Constants.INTERVAL_PERIODIC;
				// Containment
				k=neighbors[i][j].containment;
				if (k>=0) {
					if (k==Constants.CONTAINMENT_IDENTICAL) {
						if (!periodicI) {
							if (!periodicJ) edgeStats[0]++;
							else edgeStats[5]++;
						}
						else {
							if (!periodicJ) edgeStats[5]++;
							else edgeStats[12]++;
						}
					}
					else if (k==Constants.CONTAINMENT_ONE_IN_TWO) {
						if (!periodicI && !periodicJ) edgeStats[1]++;
						else if (periodicI && periodicJ) edgeStats[13]++;
						else {
							if (nodesArray[neighbors[i][j].nodeID1].type==Constants.INTERVAL_PERIODIC) edgeStats[6]++;
							else edgeStats[7]++;
						}
					}
					else if (k==Constants.CONTAINMENT_TWO_IN_ONE) {
						if (!periodicI && !periodicJ) edgeStats[1]++;
						else if (periodicI && periodicJ) edgeStats[13]++;
						else {
							if (nodesArray[neighbors[i][j].nodeID2].type==Constants.INTERVAL_PERIODIC) edgeStats[6]++;
							else edgeStats[7]++;
						}
					}
				}
				// Overlap
				k=neighbors[i][j].overlap;
				if (k>=0) {
					if (!periodicI && !periodicJ) edgeStats[2]++;
					else if (!periodicI || !periodicJ) edgeStats[8]++;
					else {
						if (k<=Constants.OVERLAP_SUFFIX_SUFFIX) edgeStats[14]++;
						else edgeStats[15]++;
					}
				}
				// Insertion
				k=neighbors[i][j].insertion;
				if (k>=0 && k!=Constants.INSERTION_BOTH) {
					if (k==Constants.INSERTION_ONE_IN_TWO) {
						if (!periodicI && !periodicJ) edgeStats[3]++;
						else if (periodicI && periodicJ) edgeStats[16]++;
						else {
							if (nodesArray[neighbors[i][j].nodeID1].type==Constants.INTERVAL_PERIODIC) edgeStats[9]++;
							else edgeStats[10]++;
						}
					}
					else if (k==Constants.INSERTION_TWO_IN_ONE) {
						if (!periodicI && !periodicJ) edgeStats[3]++;
						else if (periodicI && periodicJ) edgeStats[16]++;
						else {
							if (nodesArray[neighbors[i][j].nodeID2].type==Constants.INTERVAL_PERIODIC) edgeStats[9]++;
							else edgeStats[10]++;
						}
					}
				}
				// Shared substring
				k=neighbors[i][j].sharedSubstring;
				if (k>=0) {
					if (!periodicI && !periodicJ) edgeStats[4]++;
					else if (periodicI && periodicJ) edgeStats[17]++;
					else edgeStats[11]++;
				}
				// Periodic-nonperiodic edges
				if (periodicI!=periodicJ) edgeStats[18]++;
				// Implied marks
				if (neighbors[i][j].implied) edgeStats[19]++;
			}
		}
		for (i=0; i<edgeStats.length; i++) edgeStats[i]>>=1;
		return nEdges>>1;
	}
	
	
	public static final void printEdgeStats(int nEdges) {
		System.out.println("Edge types (not necessarily summing to 100%): ");
		System.out.println("Nonperiodic-nonperiodic:");
		System.out.println(edgeStats[0]+" containment, identical ("+IO.format((100*(double)edgeStats[0])/nEdges)+"%)");
		System.out.println(edgeStats[1]+" containment proper ("+IO.format((100*(double)edgeStats[1])/nEdges)+"%)");
		System.out.println(edgeStats[2]+" overlap ("+IO.format((100*(double)edgeStats[2])/nEdges)+"%)");
		System.out.println(edgeStats[3]+" insertion ("+IO.format((100*(double)edgeStats[3])/nEdges)+"%)");
		System.out.println(edgeStats[4]+" shared substring ("+IO.format((100*(double)edgeStats[4])/nEdges)+"%)");
		System.out.println("Nonperiodic-periodic:");
		System.out.println(edgeStats[5]+" containment, identical ("+IO.format((100*(double)edgeStats[5])/nEdges)+"%)");
		System.out.println(edgeStats[6]+" containment, periodic->nonperiodic ("+IO.format((100*(double)edgeStats[6])/nEdges)+"%)");
		System.out.println(edgeStats[7]+" containment, nonperiodic->periodic ("+IO.format((100*(double)edgeStats[7])/nEdges)+"%)");
		System.out.println(edgeStats[8]+" overlap ("+IO.format((100*(double)edgeStats[8])/nEdges)+"%)");
		System.out.println(edgeStats[9]+" insertion, periodic->nonperiodic ("+IO.format((100*(double)edgeStats[9])/nEdges)+"%)");
		System.out.println(edgeStats[10]+" insertion, nonperiodic->periodic ("+IO.format((100*(double)edgeStats[10])/nEdges)+"%)");
		System.out.println(edgeStats[11]+" shared substring ("+IO.format((100*(double)edgeStats[11])/nEdges)+"%)");
		System.out.println("Periodic-periodic:");
		System.out.println(edgeStats[12]+" containment, identical ("+IO.format((100*(double)edgeStats[12])/nEdges)+"%)");
		System.out.println(edgeStats[13]+" containment proper ("+IO.format((100*(double)edgeStats[13])/nEdges)+"%)");
		System.out.println(edgeStats[14]+" overlap, assembly ("+IO.format((100*(double)edgeStats[14])/nEdges)+"%)");
		System.out.println(edgeStats[15]+" overlap, periodic ("+IO.format((100*(double)edgeStats[15])/nEdges)+"%)");
		System.out.println(edgeStats[16]+" insertion ("+IO.format((100*(double)edgeStats[16])/nEdges)+"%)");
		System.out.println(edgeStats[17]+" shared substring ("+IO.format((100*(double)edgeStats[17])/nEdges)+"%)");
		System.out.println();
		System.out.println(edgeStats[18]+" periodic-nonperiodic marks ("+IO.format((100*(double)edgeStats[18])/nEdges)+"%)");
		System.out.println(edgeStats[19]+" implied marks ("+IO.format((100*(double)edgeStats[19])/nEdges)+"%)");
	}
	
	
	/**
	 * Stores in $tmpPoints[0..lastTmpPoint]$ the sorted and compacted set of distinct
	 * $avgDiffs$ values of all edges in the graph (both ON and OFF edges).
	 * 
	 * Remark: the procedure assumes $avgDiffs$ to be an average, not a sum.
	 * 
	 * @param subgraph[first..last] if not null, computes statistics only for the subgraph
	 * induced by those vertices; otherwise, computes statistics for the whole graph.
	 * $subgraph[first..last]$ is assumed to be sorted.
	 */
	public static final void getAvgDiffsDistribution(int[] subgraph, int first, int last, int nEdges) {
		int i, j, from, to;
		
		lastTmpPoint=-1;
		if (subgraph==null) {
			for (i=0; i<neighbors.length; i++) {
				for (j=0; j<neighbors[i].length; j++) {
					if (neighbors[i][j].avgDiffs==-1) continue;
					lastTmpPoint++;
					ensureTmpPoints(lastTmpPoint);
					tmpPoints[lastTmpPoint].position=neighbors[i][j].avgDiffs;
					tmpPoints[lastTmpPoint].mass=1;
				}
			}
		}
		else {
			for (i=first; i<=last; i++) {
				from=subgraph[i];
				for (j=0; j<neighbors[from].length; j++) {
					to=neighbors[from][j].getTo(from);
					if (Arrays.binarySearch(subgraph,first,last+1,to)<0) continue;
					if (neighbors[from][j].avgDiffs==-1) continue;
					lastTmpPoint++;
					ensureTmpPoints(lastTmpPoint);
					tmpPoints[lastTmpPoint].position=neighbors[from][j].avgDiffs;
					tmpPoints[lastTmpPoint].mass=1;
				}
			}
		}
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
	}
	
	
	public static final void printAvgDiffsDistribution() {
		System.out.println("avgDiffs:");
		for (int i=0; i<=lastTmpPoint; i++) System.out.println(tmpPoints[i]);
	}
	
	
	
	
	
	
	
	
	// --------------------------- PERIODIC SUBSTRINGS OF NODES --------------------------
	
	/**
	 * Sets the $periodicSubstrings$ list of every nonperiodic node to the set of 
	 * substrings induced by the fact that the node straddles with or contains a periodic 
	 * substring in the same read.
	 *
	 * Remark: if a dense substring (of any type) is contained in a periodic substring, 
	 * the container is added to the $periodicSubstrings$ list, since the dense substring 
	 * might be an artefact induced by the periodic substring. This does not happen with 
	 * alignment intervals, since they might be the result of a strong signal inside the 
	 * periodic substring.
	 *
	 * Remark: the procedure assumes $nodesArray$ to have already been sorted. 
	 * The procedure does not create new intervals.
	 *
	 * @param add TRUE: adds substrings to the existing $periodicSubstrings$ list, rather 
	 * than resetting the list.
	 */
	public static final void setPeriodicSubstringsOfNodes_sameRead(boolean add) {
		boolean added;
		int i, j;
		int read, start, end, count, type;
		
		count=0;
		for (i=0; i<nNodes; i++) {
			if (!add) nodesArray[i].lastPeriodicSubstring=-1;
			type=nodesArray[i].type;
			if (type==Constants.INTERVAL_PERIODIC) continue;
			added=false;
			read=nodesArray[i].read; start=nodesArray[i].start; end=nodesArray[i].end;
			for (j=i-1; j>=0; j--) {
				if (nodesArray[j].read!=read) break;
				if (nodesArray[j].type!=Constants.INTERVAL_PERIODIC) continue;
				if ( ( !Intervals.areApproximatelyIdentical(start,end,nodesArray[j].start,nodesArray[j].end) &&
					   !Intervals.isApproximatelyContained(start,end,nodesArray[j].start,nodesArray[j].end) &&
					   ( Intervals.straddles(start,end,nodesArray[j].start,nodesArray[j].end) ||
				         Intervals.isApproximatelyContained(nodesArray[j].start,nodesArray[j].end,start,end)
					   )
					 ) ||
					 ( type>=Constants.INTERVAL_DENSE_PREFIX && type<=Constants.INTERVAL_DENSE_SINGLEDELETION &&
					   Intervals.isApproximatelyContained(start,end,nodesArray[j].start,nodesArray[j].end)
					 )
				   ) {
					nodesArray[i].addPeriodicSubstring(Math.max(start,nodesArray[j].start)-start,Math.min(end,nodesArray[j].end)-start,nodesArray[j].hasLongPeriod);
					added=true;
				}
			}
			for (j=i+1; j<nNodes; j++) {
				if (nodesArray[j].read!=read || nodesArray[j].start>=end) break;
				if (nodesArray[j].type!=Constants.INTERVAL_PERIODIC) continue;
				if ( ( !Intervals.areApproximatelyIdentical(start,end,nodesArray[j].start,nodesArray[j].end) &&
					   !Intervals.isApproximatelyContained(start,end,nodesArray[j].start,nodesArray[j].end) &&
					   ( Intervals.straddles(start,end,nodesArray[j].start,nodesArray[j].end) ||
				         Intervals.isApproximatelyContained(nodesArray[j].start,nodesArray[j].end,start,end)
					   )
					 ) ||
					 ( type>=Constants.INTERVAL_DENSE_PREFIX && type<=Constants.INTERVAL_DENSE_SINGLEDELETION &&
					   Intervals.isApproximatelyContained(start,end,nodesArray[j].start,nodesArray[j].end)
					 )
				   ) {
					nodesArray[i].addPeriodicSubstring(Math.max(start,nodesArray[j].start)-start,Math.min(end,nodesArray[j].end)-start,nodesArray[j].hasLongPeriod);
					added=true;
				}
			}
			if (added) count++;
		}
		System.err.println(count+" nonperiodic nodes have at least one same-read periodic substring ("+IO.getPercent(count,nNodes)+"% of all nodes)");
	}
	
	
	/**
	 * It can happen that a periodic interval of short length is not recognized as 
	 * periodic, but is marked instead as an alignment interval. The procedure transforms
	 * into a short-period interval, every alignment interval such that all its alignments
	 * are contained in a periodic interval in their readB, and such that some of them 
	 * involve a short-period interval. This is a basic form of tag propagation across 
	 * different reads.
	 * 
	 * Remark: the procedure keeps running tag propagation until no new interval is 
	 * transformed.
	 *
	 * Remark: the procedure does not assume $nodesArray$ to be sorted, and it does not 
	 * create new intervals.
	 */
	private static final void makeAlignmentIntervalsPeriodic() {
		final int GROWTH_RATE = 100;  // Arbitrary
		boolean firstIteration;
		int i, j;
		int to, top, topPrime;
		int[] out = new int[2];
		
		firstIteration=true; top=-1;
		do {
			out[0]=0; out[1]=0;
			if (firstIteration) {
				for (i=0; i<nNodes; i++) {
					if (makeAlignmentIntervalsPeriodic_impl(i,out)) {
						for (j=0; j<nNeighbors[i]; j++) {
							to=neighbors[i][j].getTo(i);
							if (nodesArray[to].type==Constants.INTERVAL_ALIGNMENT) {
								top++;
								if (top==stack.length) {
									int[] newStack = new int[stack.length+GROWTH_RATE];
									System.arraycopy(stack,0,newStack,0,stack.length);
									stack=newStack;
								}
								stack[top]=to;
							}
						}
					}
				}
				// Compacting the list of neighbors
				if (top>0) {
					Arrays.sort(stack,0,top+1);
					j=0;
					for (i=1; i<=top; i++) {
						if (stack[i]==stack[i-1]) continue;
						j++;
						stack[j]=stack[i];
					}
					top=j;
				}
				firstIteration=false;				
			}
			else {
				topPrime=top;
				for (i=0; i<=top; i++) {
					if (makeAlignmentIntervalsPeriodic_impl(stack[i],out)) {
						for (j=0; j<nNeighbors[stack[i]]; j++) {
							to=neighbors[stack[i]][j].getTo(stack[i]);
							if (nodesArray[to].type==Constants.INTERVAL_ALIGNMENT) {
								topPrime++;
								if (topPrime==stack.length) {
									int[] newStack = new int[stack.length+GROWTH_RATE];
									System.arraycopy(stack,0,newStack,0,stack.length);
									stack=newStack;
								}
								stack[topPrime]=to;
							}
						}
					}
				}
				// Compacting the list of neighbors
				if (topPrime>top+1) {
					Arrays.sort(stack,top+1,topPrime+1);
					j=0; stack[0]=stack[top+1];
					for (i=top+2; i<=topPrime; i++) {
						if (stack[i]==stack[i-1]) continue;
						j++;
						stack[j]=stack[i];
					}
					top=j;
				}
			} 
			System.err.println("makeAlignmentIntervalsPeriodic> nCandidates="+out[0]+" nTransformed="+out[1]);		
		} while (out[1]>0);
	}
	
	
	/**
	 * @param out 0=total number of candidates for transformation; 1=total number of nodes 
	 * transformed; these counts are incremented by the procedure;
	 * @return TRUE iff the $i$-th node has been transformed. 
	 */
	private static final boolean makeAlignmentIntervalsPeriodic_impl(int i, int[] out) {
		final int MIN_SHORT_PERIOD_ALIGNMENTS = 1;  // Arbitrary
		final double PERIODIC_ALIGNMENTS_FRACTION = 0.9;  // Arbitrary
		boolean found, containedInShortPeriod, containedInLongPeriod, firstIteration;
		int j, k, h;
		int to, periodicType, nCandidates, nTransformed, nAlignments, nPeriodicAlignments;
		int projection, top;
		int nShortPeriod, nLongPeriod, nShortPeriodFull;
		Node node, nodeTo;
		Edge edge;
		
		node=nodesArray[i];
		containedInShortPeriod=(node.isContainedInterval&(1<<Constants.INTERVAL_PERIODIC))!=0;
		containedInLongPeriod=(node.isContainedInterval&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0;
		if (node.type!=Constants.INTERVAL_ALIGNMENT) return false;
		out[0]++;
		
		// Counting periodic neighbors
		nAlignments=0; nPeriodicAlignments=0; 
		nShortPeriod=0; nLongPeriod=0; nShortPeriodFull=0;
		for (j=0; j<nNeighbors[i]; j++) {
			edge=neighbors[i][j];
			to=edge.getTo(i);
			nodeTo=nodesArray[to];
			if (nodeTo.type==Constants.INTERVAL_PERIODIC) {
				nAlignments++;
				nPeriodicAlignments++;
				if (nodeTo.hasLongPeriod) nLongPeriod++;
				else {
					nShortPeriod++;
					projection=edge.projectOnto(i,tmpArray,true);
					for (k=0; k<=projection; k+=2) {
						if (Intervals.areApproximatelyIdentical(tmpArray[k],tmpArray[k+1],0,node.length()-1)) nShortPeriodFull++;
					}
				}
				continue;
			}
			found=false;
			if (edge.containment!=-1) {
				nAlignments++;
				projection=edge.projectOnto_containment(to,tmpArray,-1);
				if (projection!=-1) {
					periodicType=nodeTo.inPeriodicSubstring(tmpArray[0],tmpArray[1]);
					if (periodicType>=0) {
						found=true;
						nPeriodicAlignments++;
						if (periodicType==1) {
							nShortPeriod++;
							projection=edge.projectOnto_containment(i,tmpArray,-1);
							for (k=0; k<=projection; k+=2) {
								if (Intervals.areApproximatelyIdentical(tmpArray[k],tmpArray[k+1],0,node.length()-1)) nShortPeriodFull++;
							}
						}
						else nLongPeriod++;
					}
				}
			}
			if (edge.insertion!=-1) {
				nAlignments++;
				projection=edge.projectOnto_insertion(to,tmpArray,-1);
				if (projection!=-1) {
					periodicType=nodeTo.inPeriodicSubstring(tmpArray[0],tmpArray[1]);
					if (periodicType>=0) {
						found=true;
						nPeriodicAlignments++;
						if (periodicType==1) {
							nShortPeriod++;
							projection=edge.projectOnto_insertion(i,tmpArray,-1);
							for (k=0; k<=projection; k+=2) {
								if (Intervals.areApproximatelyIdentical(tmpArray[k],tmpArray[k+1],0,node.length()-1)) nShortPeriodFull++;
							}
						}
						else nLongPeriod++;
					}
				}
			}
			if (edge.overlap!=-1) {
				for (k=Constants.OVERLAP_PREFIX_PREFIX; k<=Constants.OVERLAP_SUFFIX_SUFFIX; k<<=1) {
					if ((edge.overlap&k)!=0) {
						nAlignments++;
						projection=edge.projectOnto_overlap(to,tmpArray,k,-1);
						if (projection!=-1) {
							periodicType=nodeTo.inPeriodicSubstring(tmpArray[0],tmpArray[1]);
							if (periodicType>=0) {
								found=true;
								nPeriodicAlignments++;
								if (periodicType==1) {
									nShortPeriod++;
									projection=edge.projectOnto_overlap(i,tmpArray,k,-1);
									for (h=0; h<=projection; h+=2) {
										if (Intervals.areApproximatelyIdentical(tmpArray[h],tmpArray[h+1],0,node.length()-1)) nShortPeriodFull++;
									}
								}
								else nLongPeriod++;
								break;
							}
						}
					}
				}
			}
			if (edge.sharedSubstring!=-1) {
				nAlignments++;
				projection=edge.projectOnto_sharedSubstring(to,tmpArray,-1);
				if (projection!=-1) {
					periodicType=nodeTo.inPeriodicSubstring(tmpArray[0],tmpArray[1]);
					if (periodicType>=0) {
						found=true;
						nPeriodicAlignments++;
						if (periodicType==1) {
							nShortPeriod++;
							projection=edge.projectOnto_sharedSubstring(i,tmpArray,-1);
							for (k=0; k<=projection; k+=2) {
								if (Intervals.areApproximatelyIdentical(tmpArray[k],tmpArray[k+1],0,node.length()-1)) nShortPeriodFull++;
							}
						}
						else nLongPeriod++;
					}
				}
			}
		}		
		if (nPeriodicAlignments<nAlignments*PERIODIC_ALIGNMENTS_FRACTION || nShortPeriodFull<MIN_SHORT_PERIOD_ALIGNMENTS) return false;
		
		// Transforming the node into periodic
		node.type=Constants.INTERVAL_PERIODIC;
		node.period=-1;
		node.hasLongPeriod=false;
		node.lastPeriodicSubstring=-1;
		out[1]++;
		for (j=0; j<nNeighbors[i]; j++) {
			edge=neighbors[i][j];
			edge.changeType_periodic(tmpEdge,tmpEdge2);
			tmpEdge2.clone(edge);
		}
		return true;
	}
	
	
	/**
	 * It often happens that long, weak dense substrings of substring type span multiple 
	 * distinct repeats. Such repeats might actually have been detected as separate 
	 * intervals during factorization. However, it could happen that an alignment assigned 
	 * to the weak substring should instead have been assigned to a contained or 
	 * straddling interval. Due to issues in the aligner and in the input, this could 
	 * happen even if the alignment cannot be assigned to the contained interval during 
	 * factorization, based on e.g. maximality and implication.
     * The procedure takes the conservative approach of redirecting an edge incident to 
	 * such a weak substring, to a shortest non-periodic interval contained in the 
	 * substring or straddling with it, if any.
	 *
	 * Similarly, an edge between a long-period interval, and a non-periodic interval that
	 * is significantly shorter than the long period and has enough assigned alignments, 
	 * is not trustworthy, since the repeat should have been detected as a separate 
	 * interval inside the long-period interval as well, and the alignment should have 
	 * been assigned to such contained interval.
	 * The procedure takes the conservative approach of redirecting an edge incident to 
	 * such a long-period interval, to a shortest non-periodic interval contained in it, 
	 * if any; if no such contained interval can be found, the procedure removes the edge.
	 *
	 * Remark: the procedure does not create new periodic-nonperiodic edges.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$.
	 * The procedure does not create new intervals.
	 *
	 * @param weakOrLongPeriod TRUE=weak substring-type; FALSE=long period intervals.
	 */
	private static final void redirectWeakSubstrings(boolean weakOrLongPeriod) {
		final int DISTANCE = IO.quantum<<1;  // Arbitrary
		final int MIN_WEAK_LENGTH = Alignments.minAlignmentLength<<1;  // Arbitrary
		final int MIN_INTERSECTION_LENGTH = Alignments.minAlignmentLength;  // Arbitrary
		final int MIN_CONTAINED_INTERVALS = 2;  // Arbitrary
		final double MIN_JACCARD = 0.8;  // Arbitrary
		final int MAX_DIFFERENCE = Alignments.minAlignmentLength/3;  // Arbitrary
		final int GROWTH_RATE = 100;  // Arbitrary
		final int MIN_COVERAGE = IO.coverage*10;  // Arbitrary
		final double PERIOD_RATE = 0.8;  // Arbitrary
		int i, j, k, h;
		int to, read, start, end, nodeID1;
		int lastContained, lastNewEdge, lastNewEdgePrime, nSharedSubstrings, nStoredAlignments, projection, interval;
		int redirectedEdges, redirectedEdgesKept, processedNodes, nNewEdges, period;
		double avgDiffs;
		Node node, otherNode, nodeTo;
		Edge edge, newEdge;
		AlignmentPair pair = new AlignmentPair();
		int[] containedIntervals = new int[GROWTH_RATE];
		Edge[] newEdges = new Edge[GROWTH_RATE];
		
		redirectedEdges=0; redirectedEdgesKept=0; 
		processedNodes=0; nNewEdges=0;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if ( (weakOrLongPeriod && (node.type!=Constants.INTERVAL_DENSE_SUBSTRING || !node.isWeak || node.length()<MIN_WEAK_LENGTH)) ||
				 (!weakOrLongPeriod && (node.type!=Constants.INTERVAL_PERIODIC || !node.hasLongPeriod))
			   ) continue;
			read=node.read; start=node.start; end=node.end; period=node.period;
			
			// Loading non-periodic sub-intervals
			lastContained=-1;
			for (j=i-1; j>=0; j--) {
				otherNode=nodesArray[j];
				if (otherNode.read!=read) break;
				if (otherNode.end<=start || otherNode.type==Constants.INTERVAL_PERIODIC) continue;
				if ( ( weakOrLongPeriod && 
					   ( Intervals.isApproximatelyContained(otherNode.start,otherNode.end,start,end) || 
				         (!Intervals.isApproximatelyContained(start,end,otherNode.start,otherNode.end) && Intervals.intersectionLength(otherNode.start,otherNode.end,start,end)>=MIN_INTERSECTION_LENGTH)
					   )
					 ) ||
					 (!weakOrLongPeriod && Intervals.isApproximatelyContained(otherNode.start,otherNode.end,start,end))
				   ) {
					lastContained++;
					if (lastContained==containedIntervals.length) {
						int[] containedIntervalsPrime = new int[containedIntervals.length+GROWTH_RATE];
						System.arraycopy(containedIntervals,0,containedIntervalsPrime,0,containedIntervals.length);
						containedIntervals=containedIntervalsPrime;
					}
					containedIntervals[lastContained]=j;
				}
			}
			for (j=i+1; j<nNodes; j++) {
				otherNode=nodesArray[j];
				if (otherNode.read!=read || otherNode.start>=end) break;
				if (otherNode.type==Constants.INTERVAL_PERIODIC) continue;
				if ( ( weakOrLongPeriod && 
					   ( Intervals.isApproximatelyContained(otherNode.start,otherNode.end,start,end) || 
				         (!Intervals.isApproximatelyContained(start,end,otherNode.start,otherNode.end) && Intervals.intersectionLength(otherNode.start,otherNode.end,start,end)>=MIN_INTERSECTION_LENGTH)
					   )
					 ) ||
				     (!weakOrLongPeriod && Intervals.isApproximatelyContained(otherNode.start,otherNode.end,start,end))
				   ) {
					lastContained++;
					if (lastContained==containedIntervals.length) {
						int[] containedIntervalsPrime = new int[containedIntervals.length+GROWTH_RATE];
						System.arraycopy(containedIntervals,0,containedIntervalsPrime,0,containedIntervals.length);
						containedIntervals=containedIntervalsPrime;
					}
					containedIntervals[lastContained]=j;
				}
			}
			if (weakOrLongPeriod && lastContained+1<MIN_CONTAINED_INTERVALS) continue;
			
			// Redirecting edges
			lastNewEdge=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (edge==null) continue;
				to=edge.getTo(i); nodeTo=nodesArray[to];
				if ( nodeTo.read==read ||
					 (!weakOrLongPeriod && (nodeTo.type==Constants.INTERVAL_PERIODIC || (period>0 && nodeTo.length()>=period*PERIOD_RATE)))
				   ) continue;
				nSharedSubstrings=edge.sharedSubstring==-1?0:edge.getNAlignments()-(edge.containment==-1?0:1)-(edge.insertion==-1?0:1)-edge.getNOverlaps();  // Upper bound
				nStoredAlignments=edge.nStoredAlignments();
				avgDiffs=edge.avgDiffs/edge.getNAlignments();
				lastNewEdgePrime=-1;
				
				// Containment sub-edge
				if (edge.containment!=-1) {
					projection=edge.projectOnto_containment(i,tmpArray,-1);
					if (projection!=-1 && !Intervals.areApproximatelyIdentical(tmpArray[0],tmpArray[1],0,node.length()-1)) {
						interval=redirectWeakSubstrings_getContainedInterval(start+tmpArray[0],start+tmpArray[1],containedIntervals,lastContained,weakOrLongPeriod?MIN_JACCARD:Math.POSITIVE_INFINITY,weakOrLongPeriod?MAX_DIFFERENCE:-1);
						if (interval>=0) {
							pair.node1=nodesArray[interval];
							pair.alignmentStart1=start+tmpArray[0];
							pair.alignmentEnd1=start+tmpArray[1];
							pair.node2=nodesArray[to];
							projection=edge.projectOnto_containment(to,tmpArray,-1);
							if (projection!=-1) {
								pair.alignmentStart2=nodeTo.start+tmpArray[0];
								pair.alignmentEnd2=nodeTo.start+tmpArray[1];
								pair.orientation=edge.orientation==0;
								pair.diffs=avgDiffs;
								pair.inCanonicalForm(); 
								if ( ( Intervals.isApproximatelyContained(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end) ||  
									   Intervals.areApproximatelyIdentical(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end)
									 ) &&
								     ( Intervals.isApproximatelyContained(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end) ||  
									   Intervals.areApproximatelyIdentical(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end)
									 )
								   ) {
									// Alignments that do not satisfy the conditions
									// above can be assigned the wrong Step2 type after
									// trim().
									pair.trim();
									newEdge = new Edge(pair);
									if (edge.orientation==0) newEdge.orientation=0;
									else if (edge.orientation==1) newEdge.orientation=1;
									else newEdge.orientation=-1;
									lastNewEdgePrime++; lastNewEdge++;
									if (lastNewEdge==newEdges.length) {
										Edge[] newEdgesPrime = new Edge[newEdges.length+GROWTH_RATE];
										System.arraycopy(newEdges,0,newEdgesPrime,0,newEdges.length);
										newEdges=newEdgesPrime;
									}
									newEdges[lastNewEdge]=newEdge;
									edge.containment=-1;
								}
							}
						}
					}
				}
				
				// Insertion sub-edge
				if (edge.insertion!=-1) {
					projection=edge.projectOnto_insertion(i,tmpArray,-1);
					if (projection!=-1 && !Intervals.areApproximatelyIdentical(tmpArray[0],tmpArray[1],0,node.length()-1)) {
						interval=redirectWeakSubstrings_getContainedInterval(start+tmpArray[0],start+tmpArray[1],containedIntervals,lastContained,weakOrLongPeriod?MIN_JACCARD:Math.POSITIVE_INFINITY,weakOrLongPeriod?MAX_DIFFERENCE:-1);
						if (interval>=0) {
							pair.node1=nodesArray[interval];
							pair.alignmentStart1=start+tmpArray[0];
							pair.alignmentEnd1=start+tmpArray[1];
							pair.node2=nodesArray[to];
							projection=edge.projectOnto_insertion(to,tmpArray,-1);
							if (projection!=-1) {
								pair.alignmentStart2=nodeTo.start+tmpArray[0];
								pair.alignmentEnd2=nodeTo.start+tmpArray[1];
								pair.orientation=edge.orientation==0;
								pair.diffs=avgDiffs;
								pair.inCanonicalForm(); 
								if ( ( Intervals.isApproximatelyContained(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end) ||  
									   Intervals.areApproximatelyIdentical(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end)
									 ) &&
								     ( Intervals.isApproximatelyContained(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end) ||  
									   Intervals.areApproximatelyIdentical(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end)
									 )
								   ) {
									// Alignments that do not satisfy the conditions
									// above can be assigned the wrong Step2 type after
									// trim().
									pair.trim();
									newEdge = new Edge(pair);
									if (edge.orientation==0) newEdge.orientation=0;
									else if (edge.orientation==1) newEdge.orientation=1;
									else newEdge.orientation=-1;
									lastNewEdgePrime++; lastNewEdge++;
									if (lastNewEdge==newEdges.length) {
										Edge[] newEdgesPrime = new Edge[newEdges.length+GROWTH_RATE];
										System.arraycopy(newEdges,0,newEdgesPrime,0,newEdges.length);
										newEdges=newEdgesPrime;
									}
									newEdges[lastNewEdge]=newEdge;
									edge.insertion=-1;
								}
							}
						}
					}
				}
				
				// Overlap sub-edges
				if (edge.overlap!=-1) {
					for (k=Constants.OVERLAP_PREFIX_PREFIX; k<=Constants.OVERLAP_SUFFIX_SUFFIX; k<<=1) {
						if ((edge.overlap&k)!=0) {
							projection=edge.projectOnto_overlap(i,tmpArray,k,-1);
							if (projection!=-1 && !Intervals.areApproximatelyIdentical(tmpArray[0],tmpArray[1],0,node.length()-1)) {
								interval=redirectWeakSubstrings_getContainedInterval(start+tmpArray[0],start+tmpArray[1],containedIntervals,lastContained,weakOrLongPeriod?MIN_JACCARD:Math.POSITIVE_INFINITY,weakOrLongPeriod?MAX_DIFFERENCE:-1);
								if (interval>=0) {
									pair.node1=nodesArray[interval];
									pair.alignmentStart1=start+tmpArray[0];
									pair.alignmentEnd1=start+tmpArray[1];
									pair.node2=nodesArray[to];
									projection=edge.projectOnto_overlap(to,tmpArray,k,-1);
									if (projection!=-1) {
										pair.alignmentStart2=nodeTo.start+tmpArray[0];
										pair.alignmentEnd2=nodeTo.start+tmpArray[1];
										pair.orientation=edge.orientation==0;
										pair.diffs=avgDiffs;
										pair.inCanonicalForm(); 
										if ( ( Intervals.isApproximatelyContained(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end) ||  
											   Intervals.areApproximatelyIdentical(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end)
											 ) &&
										     ( Intervals.isApproximatelyContained(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end) ||  
											   Intervals.areApproximatelyIdentical(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end)
											 )
										   ) {
											// Alignments that do not satisfy the
											// conditions above can be assigned the wrong
											// Step2 type after trim().
											pair.trim();
											newEdge = new Edge(pair);
											if (edge.orientation==0) newEdge.orientation=0;
											else if (edge.orientation==1) newEdge.orientation=1;
											else newEdge.orientation=-1;
											lastNewEdgePrime++; lastNewEdge++;
											if (lastNewEdge==newEdges.length) {
												Edge[] newEdgesPrime = new Edge[newEdges.length+GROWTH_RATE];
												System.arraycopy(newEdges,0,newEdgesPrime,0,newEdges.length);
												newEdges=newEdgesPrime;
											}
											newEdges[lastNewEdge]=newEdge;
											edge.removeOverlap(k);
										}
									}
								}
							}
						}
					}
					if (edge.overlap==0) edge.overlap=-1;
				}	
				
				// Shared substring sub-edge
				if (edge.sharedSubstring!=-1) {
					projection=edge.projectOnto_sharedSubstring(i,tmpArray,-1);  // Projecting one of the possibly many shared substring intervals (only this one has been stored).
					if (projection!=-1 && !Intervals.areApproximatelyIdentical(tmpArray[0],tmpArray[1],0,node.length()-1)) {
						interval=redirectWeakSubstrings_getContainedInterval(start+tmpArray[0],start+tmpArray[1],containedIntervals,lastContained,weakOrLongPeriod?MIN_JACCARD:Math.POSITIVE_INFINITY,weakOrLongPeriod?MAX_DIFFERENCE:-1);
						if (interval>=0) {
							pair.node1=nodesArray[interval];
							pair.alignmentStart1=start+tmpArray[0];
							pair.alignmentEnd1=start+tmpArray[1];
							pair.node2=nodesArray[to];
							projection=edge.projectOnto_sharedSubstring(to,tmpArray,-1);
							if (projection!=-1) {
								pair.alignmentStart2=nodeTo.start+tmpArray[0];
								pair.alignmentEnd2=nodeTo.start+tmpArray[1];
								pair.orientation=edge.orientation==0;
								pair.diffs=avgDiffs;
								pair.inCanonicalForm(); 
								if ( ( Intervals.isApproximatelyContained(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end) ||  
									   Intervals.areApproximatelyIdentical(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end)
									 ) &&
								     ( Intervals.isApproximatelyContained(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end) ||  
									   Intervals.areApproximatelyIdentical(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end)
									 )
								   ) {
									// Alignments that do not satisfy the conditions
									// above can be assigned the wrong Step2 type after
									// trim().
									pair.trim();
									newEdge = new Edge(pair);
									if (edge.orientation==0) newEdge.orientation=0;
									else if (edge.orientation==1) newEdge.orientation=1;
									else newEdge.orientation=-1;
									lastNewEdgePrime++; lastNewEdge++;
									if (lastNewEdge==newEdges.length) {
										Edge[] newEdgesPrime = new Edge[newEdges.length+GROWTH_RATE];
										System.arraycopy(newEdges,0,newEdgesPrime,0,newEdges.length);
										newEdges=newEdgesPrime;
									}
									newEdges[lastNewEdge]=newEdge;
									if (nSharedSubstrings==1) edge.sharedSubstring=-1;
								}
							}
						}
					}
				}
					
				// Removing the current edge iff all edges with stored coordinates have 
				// been used (see $fixPeriodicUndersplits()$ for details). The current
				// edge is always removed if $weakOrLongPeriod=false$.
				if (!weakOrLongPeriod || lastNewEdgePrime>=0) {
					redirectedEdges++;
					if (!weakOrLongPeriod || lastNewEdgePrime+1==edge.getNAlignments() || lastNewEdgePrime+1==nStoredAlignments || edge.nDistinctTypes()==0) {
						neighbors[i][j]=null;
						for (k=0; k<nNeighbors[to]; k++) {
							if (neighbors[to][k]==edge) {
								neighbors[to][k]=null;
								break;
							}
						}
					}
					else {
						redirectedEdgesKept++;
						if (edge.orientation==0) edge.nAlignmentsForward-=lastNewEdgePrime+1;
						else if (edge.orientation==1) edge.nAlignmentsBackward-=lastNewEdgePrime+1;
						else {
							edge.nAlignmentsForward=edge.getNAlignments()-lastNewEdgePrime+1;  // Arbitrarily setting forward rather than backward alignments.
							edge.nAlignmentsBackward=0;
							edge.orientation=-1;
						}
						edge.avgDiffs=avgDiffs*edge.getNAlignments();
					}
				}
			}
			if (lastNewEdge==-1) continue;
			processedNodes++;
			
			// Adding new edges
			for (j=0; j<=lastNewEdge; j++) {
				edge=newEdges[j];
				nodeID1=edge.nodeID1;
				k=findEdge(nodeID1,edge,false);
				if (k==Math.NEGATIVE_INFINITY) {
					addEdge(nodeID1,edge);
					addEdge(edge.nodeID2,edge);
					nNewEdges++;
				}
				else {
					neighbors[nodeID1][k].addTypes(edge);
					neighbors[nodeID1][k].addDiffs(edge);
					neighbors[nodeID1][k].addOrientation(edge);
				}
			}
		}
		
		// Removing $null$ neighbors
		for (i=0; i<nNodes; i++) {
			k=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				neighbors[i][++k]=neighbors[i][j];
			}
			nNeighbors[i]=k+1;
		}
		
		System.err.println("redirectWeakSubstrings> "+processedNodes+" nodes processed, "+redirectedEdges+" edges redirected ("+redirectedEdgesKept+" kept), "+nNewEdges+" new edges created.");
	}
	
	
	/**
	 * @return a shortest interval in $containedIntervals$ that approximately contains
	 * $[start..end]$, or -1 if no such interval can be found.
	 */
	private static final int redirectWeakSubstrings_getContainedInterval(int start, int end, int[] containedIntervals, int lastContained, double minJaccard, int maxDifference) {
		int i;
		int length, minLength, minInterval, minLengthIntersection, minIntervalIntersection;
		Node node;
		
		minInterval=-1; minLength=Math.POSITIVE_INFINITY;
		minIntervalIntersection=-1; minLengthIntersection=Math.POSITIVE_INFINITY;
		for (i=0; i<=lastContained; i++) {
			node=nodesArray[containedIntervals[i]];
			if (Intervals.isApproximatelyContained(start,end,node.start,node.end)) {
				length=node.length();
				if (length<minLength) {
					minLength=length;
					minInterval=containedIntervals[i];
				}
			}
			else if (Intervals.jaccardSimilarity(start,end,node.start,node.end)>=minJaccard && Intervals.difference(start,end,node.start,node.end)<=maxDifference) {
				length=node.length();
				if (length<minLengthIntersection) {
					minLengthIntersection=length;
					minIntervalIntersection=containedIntervals[i];
				}
			}
		}
		return minInterval>=0?minInterval:minIntervalIntersection;
	}
	
	
	/**
	 * Nodes with a low coverage computed on their assigned alignments are likely to be 
	 * unreliable repeat annotations, or to be the concatenation of distinct repeats. The 
	 * procedure disconnects all such intervals.
	 */
	private static final void disconnectLowCoverageNodes() {
		final int MIN_COVERAGE = IO.minRepeatCoverage<<1;  // Arbitrary
		int i, j, k;
		int to;
		Node node;
		Edge edge;
		
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if ( ( node.type==Constants.INTERVAL_ALIGNMENT || 
				   (node.type==Constants.INTERVAL_DENSE_SUBSTRING && node.isWeak)
				 ) && 
				 node.avgCoverage<MIN_COVERAGE
			   ) {
				for (j=0; j<nNeighbors[i]; j++) {
					edge=neighbors[i][j];
					if (edge==null) continue;
					to=edge.getTo(i);
					neighbors[i][j]=null;
					for (k=0; k<nNeighbors[to]; k++) {
						if (neighbors[to][k]==edge) {
							neighbors[to][k]=null;
							break;
						}
					}
				}
			}			
		}
		
		// Removing $null$ neighbors
		for (i=0; i<nNodes; i++) {
			k=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				neighbors[i][++k]=neighbors[i][j];
			}
			nNeighbors[i]=k+1;
		}
	}

	
	/**
	 * During factorization, an alignment is assigned to a periodic substring, with a 
	 * specific readB, if it belongs to the periodic pattern of the substring, or if its 
	 * intervals in both readA and readB fall inside the intervals of the substring. The 
	 * periodic substring in readA is then merged into a periodic substring interval, 
	 * whose range on readA might become larger. An alignment that falls inside a periodic
	 * substring interval in readA might not be assigned to it, e.g. because readA and 
	 * readB might not have traced a periodic substring at all, the periodic substring 
	 * interval in readA might be contained in or straddle a dense substring, and the 
	 * alignment might have been assigned to such dense substring. This is probably the 
	 * best we can do when each read pile is factorized separately. This is a special case
	 * of the more general problem that an alignment might be assigned to a container 
	 * interval rather than to a contained interval: see $fixUndersplits()$.
	 *
	 * During interval graph construction, however, we know the factorization of both 
	 * reads of each alignment. Given an edge, let one of its nodes be "red" iff it is a
	 * nonperiodic interval, and the projection of the alignment that induces the edge 
	 * falls inside a periodic interval in the same read. The procedure redirects every
	 * edge such that both its nodes are red, or such that one node is red and the other 
	 * is periodic, to corresponding periodic intervals. In other words, if an alignment
	 * falls inside periodic intervals in both reads, it is reassigned to such intervals.
	 * If it falls in a periodic interval in just one read, its assignment is not changed.
	 *
	 * Remark: if, as a result of the reassignment, the overlap between two nodes is too
	 * short, the new edge is not created, but the old edge is removed anyway. This might 
	 * be overly aggressive.
	 *
	 * Remark: recall that periodic substrings are the only substring type that is 
	 * detected by using alignment patterns that involve readA and a specific readB, 
	 * rather than readA and all readB's. No other substring type is detected in the same 
	 * way and can thus cause the problem above systematically.
	 *
	 * Remark: an edge is redirected only if the periodic substrings that contain the
	 * alignments at its ends are both of the same type, and, if they have a long period, 
	 * only if the alignments are at least as long as the period.
	 *
	 * Remark: new edges might have unknown orientation, since the procedure unpacks all
	 * alignments that are collapsed into the same original edge. However, such new edges
	 * are periodic-periodic, so their orientation is not necessarily important. One could
	 * try and fix all unknown orientations downstream, but this is likely to be overkill 
	 * for the moment. See e.g. \cite{baaijens2017novo}, in which they start from a node 
	 * of minimal degree to assign orientations, and if there are conflicts in the 
	 * resulting orientation, they search for an orientation that leads to a minimal 
	 * amount of conflicts.
	 *
	 * Remark: if $periodicOrDense=FALSE$, a variant of this procedure is applied to dense 
	 * substrings of substring type, which are now assumed to be stored in the 
	 * $periodicSubstrings$ array of each node. This is necessary, since e.g. a weak dense 
	 * substring of substring type might coincide with a suffix of a suffix substring, so, 
	 * using just the readA pile, it is not clear whether some alignments should be 
	 * assigned to the suffix substring or to the dense substring of substring type. Note 
	 * that the factorization procedure likely assigns such alignments to the suffix 
	 * substring, since weak dense substrings are assigned just to alignments that are not
	 * already implied by non-weak substrings. This procedure reassigns such alignments to
	 * the dense substring of substring type if they can be assigned in this way in readB
	 * as well.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$. 
	 * The procedure does not create new intervals.
	 *
	 * @param tmp temporary space, of size at least two.
	 */
	private static final void fixPeriodicUndersplits(boolean periodicOrDense, int[] tmp) {
		final int MIN_NONPERIODIC_NEIGHBORS = 3;
		final int MIN_ALIGNMENT_LENGTH = Alignments.minAlignmentLength-IO.quantum;  // Arbitrary
		boolean containedInShortPeriod, containedInLongPeriod;
		int i, j, k;
		int read, to, nodeID1, nodeID2, orientation, projection, period, nPeriodicNeighbors;
		int type1, type2, type1Prime, type2Prime, nStoredAlignments;
		int nNewEdges, nNewEdgesPrime, newEdges, decomposedEdges, decomposedEdgesKept, nSharedSubstrings;
		double avgDiffs;
		Node node, nodeTo;
		Edge edge;
		AlignmentPair pair = new AlignmentPair();
		
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				// Using the $printed$ tag to mark processed and new edges.
				neighbors[i][j].printed=-1;
			}
		}
		decomposedEdges=0; decomposedEdgesKept=0; newEdges=0;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			containedInShortPeriod=(node.isContainedInterval&(1<<Constants.INTERVAL_PERIODIC))!=0;
			containedInLongPeriod=(node.isContainedInterval&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0;
			if ( ( periodicOrDense && 
				   ( ( node.type==Constants.INTERVAL_PERIODIC ||
				       (node.type==Constants.INTERVAL_ALIGNMENT && (containedInShortPeriod||containedInLongPeriod))
				     ) || 
				     node.lastPeriodicSubstring==-1
				   )
				 ) ||	 
				 ( !periodicOrDense && 
				   ( node.type==Constants.INTERVAL_DENSE_SUBSTRING ||
					 node.lastPeriodicSubstring==-1
				   )
				 )
			   ) continue;
			read=node.read;
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];	
				if (edge==null) continue;
				if (edge==null || edge.printed==1) continue;
				to=edge.getTo(i); nodeTo=nodesArray[to];
				if (nodeTo.read==read) continue;
				orientation=edge.orientation;
				if (edge.getNAlignments()>1) {
					if (orientation==2) orientation=-1;
				}
				avgDiffs=edge.avgDiffs/edge.getNAlignments();
				nSharedSubstrings=edge.sharedSubstring==-1?0:edge.getNAlignments()-(edge.containment==-1?0:1)-(edge.insertion==-1?0:1)-edge.getNOverlaps();  // Upper bound
				nStoredAlignments=edge.nStoredAlignments();
				nNewEdges=0; nNewEdgesPrime=0;
				
				// Containment sub-edge
				if (edge.containment!=-1 && edge.containment!=Constants.CONTAINMENT_IDENTICAL && edge.containment!=(i==edge.nodeID1?Constants.CONTAINMENT_ONE_IN_TWO:Constants.CONTAINMENT_TWO_IN_ONE)) {
					nodeID1=-1; nodeID2=-1;
					projection=edge.projectOnto_containment(i,tmp,-1);
					if (projection!=-1) {
						type1=node.inPeriodicSubstring(tmp[0],tmp[1]);
						if (type1==0 && periodicOrDense) type1Prime=node.intersectsPeriodicSubstring(tmp[0],tmp[1]);
						else type1Prime=type1;
						if (type1Prime>0) {
							nodeID1=getPeriodicInterval(node.start+tmp[0],node.start+tmp[1],read,i,periodicOrDense,type1Prime);
							if (nodeID1>=0) {
								period=periodicOrDense?nodesArray[nodeID1].period:-1;
								if (type1Prime==1 || (type1Prime==2 && (period<0 || tmp[1]-tmp[0]+1>=period))) {
									pair.node1=nodesArray[nodeID1];
									pair.alignmentStart1=node.start+tmp[0];
									pair.alignmentEnd1=node.start+tmp[1];
									projection=edge.projectOnto_containment(to,tmp,-1);
									if (projection!=-1) {
										containedInShortPeriod=(nodeTo.isContainedInterval&(1<<Constants.INTERVAL_PERIODIC))!=0;
										containedInLongPeriod=(nodeTo.isContainedInterval&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0;
				   						if ( ( periodicOrDense && 
											   ( nodeTo.type==Constants.INTERVAL_PERIODIC ||
												 (nodeTo.type==Constants.INTERVAL_ALIGNMENT && (containedInShortPeriod||containedInLongPeriod))
											   )
											 ) ||
											 (!periodicOrDense && nodeTo.type==Constants.INTERVAL_DENSE_SUBSTRING)
										   ) {
				   							if ( ( periodicOrDense &&
											       ( (nodeTo.type==Constants.INTERVAL_PERIODIC && nodeTo.hasLongPeriod==(type1Prime==2)) ||
												     (nodeTo.type==Constants.INTERVAL_ALIGNMENT && ((containedInShortPeriod&&(type1Prime==1))||(containedInLongPeriod&&(type1Prime==2))))
												   )
												 ) ||
												 !periodicOrDense
											   ) {
												nodeID2=to;
												pair.node2=nodesArray[nodeID2];
												pair.alignmentStart2=nodeTo.start+tmp[0];
												pair.alignmentEnd2=nodeTo.start+tmp[1];
											}
				   						}
				   						else if (periodicOrDense || nodeTo.type!=Constants.INTERVAL_PERIODIC) {
											type2=nodeTo.inPeriodicSubstring(tmp[0],tmp[1]);
											if (type2==0 && periodicOrDense) type2Prime=nodeTo.intersectsPeriodicSubstring(tmp[0],tmp[1]);
											else type2Prime=type2;
											if (type2Prime==type1Prime) {
												nodeID2=getPeriodicInterval(nodeTo.start+tmp[0],nodeTo.start+tmp[1],nodeTo.read,to,periodicOrDense,type2Prime);
												if (nodeID2>=0) {
													period=periodicOrDense?nodesArray[nodeID2].period:-1;
													if (type2Prime==2 && period>0 && tmp[1]-tmp[0]+1<period) nodeID2=-1;
													else {
														pair.node2=nodesArray[nodeID2];
														pair.alignmentStart2=nodeTo.start+tmp[0];
														pair.alignmentEnd2=nodeTo.start+tmp[1];
													}
												}
											}
										}
									}
								}
							}
						}
					}
					if ( nodeID1!=-1 && nodeID2!=-1 && 
						 ( Intervals.isApproximatelyContained(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end) ||  
					   	   Intervals.areApproximatelyIdentical(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end)
					 	 ) &&
				         ( Intervals.isApproximatelyContained(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end) ||  
					   	   Intervals.areApproximatelyIdentical(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end)
					 	 )
				       ) {
						// Alignments that do not satisfy the conditions above can be
						// assigned the wrong Step2 type after trim().
						pair.trim();
						if (pair.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && pair.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
							pair.orientation=orientation==0;
							pair.diffs=avgDiffs;
							fixPeriodicUndersplits_addEdge(pair,edge,i,to);
							nNewEdgesPrime++;
						}
						else {
							// NOP. The original edge might be disconnected as a result.
						}
						nNewEdges++;
						edge.containment=-1;
					}
				}
				
				// Insertion sub-edge
				if (edge.insertion!=-1 && edge.insertion!=Constants.INSERTION_BOTH && edge.insertion!=(i==edge.nodeID1?Constants.INSERTION_ONE_IN_TWO:Constants.INSERTION_TWO_IN_ONE)) {
					nodeID1=-1; nodeID2=-1;
					projection=edge.projectOnto_insertion(i,tmp,-1);
					if (projection!=-1) {
						type1=node.inPeriodicSubstring(tmp[0],tmp[1]);
						if (type1==0 && periodicOrDense) type1Prime=node.intersectsPeriodicSubstring(tmp[0],tmp[1]);
						else type1Prime=type1;
						if (type1Prime>0) {
							nodeID1=getPeriodicInterval(node.start+tmp[0],node.start+tmp[1],read,i,periodicOrDense,type1Prime);
							if (nodeID1>=0) {
								period=periodicOrDense?nodesArray[nodeID1].period:-1;
								if (type1Prime==1 || (type1Prime==2 && (period<0 || tmp[1]-tmp[0]+1>=period))) {
									pair.node1=nodesArray[nodeID1];
									pair.alignmentStart1=node.start+tmp[0];
									pair.alignmentEnd1=node.start+tmp[1];
									projection=edge.projectOnto_insertion(to,tmp,-1);
									if (projection!=-1) {
										containedInShortPeriod=(nodeTo.isContainedInterval&(1<<Constants.INTERVAL_PERIODIC))!=0;
										containedInLongPeriod=(nodeTo.isContainedInterval&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0;
										if ( ( periodicOrDense && 
											   ( nodeTo.type==Constants.INTERVAL_PERIODIC ||
												 (nodeTo.type==Constants.INTERVAL_ALIGNMENT && (containedInShortPeriod||containedInLongPeriod))
											   )
											 ) ||
											 (!periodicOrDense && nodeTo.type==Constants.INTERVAL_DENSE_SUBSTRING)
										   ) {
											if ( ( periodicOrDense &&
											       ( (nodeTo.type==Constants.INTERVAL_PERIODIC && nodeTo.hasLongPeriod==(type1Prime==2)) ||
												     (nodeTo.type==Constants.INTERVAL_ALIGNMENT && ((containedInShortPeriod&&(type1Prime==1))||(containedInLongPeriod&&(type1Prime==2))))
												   )
												 ) ||
												 !periodicOrDense
											   ) {
												nodeID2=to;
												pair.node2=nodesArray[nodeID2];
												pair.alignmentStart2=nodeTo.start+tmp[0];
												pair.alignmentEnd2=nodeTo.start+tmp[1];
											}
										}
				   						else if (periodicOrDense || nodeTo.type!=Constants.INTERVAL_PERIODIC) {
											type2=nodeTo.inPeriodicSubstring(tmp[0],tmp[1]);
											if (type2==0 && periodicOrDense) type2Prime=nodeTo.intersectsPeriodicSubstring(tmp[0],tmp[1]);
											else type2Prime=type2;
											if (type2Prime==type1Prime) {
												nodeID2=getPeriodicInterval(nodeTo.start+tmp[0],nodeTo.start+tmp[1],nodeTo.read,to,periodicOrDense,type2Prime);
												if (nodeID2>=0) {
													period=periodicOrDense?nodesArray[nodeID2].period:-1;
													if (type2Prime==2 && period>0 && tmp[1]-tmp[0]+1<period) nodeID2=-1;
													else {
														pair.node2=nodesArray[nodeID2];
														pair.alignmentStart2=nodeTo.start+tmp[0];
														pair.alignmentEnd2=nodeTo.start+tmp[1];
													}
												}
											}
										}
									}
								}
							}
						}
					}
					if ( nodeID1!=-1 && nodeID2!=-1 &&
						 ( Intervals.isApproximatelyContained(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end) ||  
					   	   Intervals.areApproximatelyIdentical(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end)
					 	 ) &&
				         ( Intervals.isApproximatelyContained(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end) ||  
					   	   Intervals.areApproximatelyIdentical(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end)
					 	 )
					   ) {
						pair.trim();
						if (pair.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && pair.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
							pair.orientation=orientation==0;
							pair.diffs=avgDiffs;
							fixPeriodicUndersplits_addEdge(pair,edge,i,to);
							nNewEdgesPrime++;
						}
						else {
							// NOP. The original edge might be disconnected as a result.
						}
						nNewEdges++;
						edge.insertion=-1;
					}
				}
				
				// Overlap sub-edges
				if (edge.overlap!=-1) {
					for (k=Constants.OVERLAP_PREFIX_PREFIX; k<=Constants.OVERLAP_SUFFIX_SUFFIX; k<<=1) {
						if ((edge.overlap&k)!=0) {
							nodeID1=-1; nodeID2=-1;
							projection=edge.projectOnto_overlap(i,tmp,k,-1);
							if (projection!=-1) {
								type1=node.inPeriodicSubstring(tmp[0],tmp[1]);
								if (type1==0 && periodicOrDense) type1Prime=node.intersectsPeriodicSubstring(tmp[0],tmp[1]);
								else type1Prime=type1;
								if (type1Prime>0) {
									nodeID1=getPeriodicInterval(node.start+tmp[0],node.start+tmp[1],read,i,periodicOrDense,type1Prime);
									if (nodeID1>=0) {
										period=periodicOrDense?nodesArray[nodeID1].period:-1;
										if (type1Prime==1 || (type1Prime==2 && (period<0 || tmp[1]-tmp[0]+1>=period))) {
											pair.node1=nodesArray[nodeID1];
											pair.alignmentStart1=node.start+tmp[0];
											pair.alignmentEnd1=node.start+tmp[1];
											projection=edge.projectOnto_overlap(to,tmp,k,-1);
											if (projection!=-1) {
												containedInShortPeriod=(nodeTo.isContainedInterval&(1<<Constants.INTERVAL_PERIODIC))!=0;
												containedInLongPeriod=(nodeTo.isContainedInterval&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0;
												if ( ( periodicOrDense && 
													   ( nodeTo.type==Constants.INTERVAL_PERIODIC ||
														 (nodeTo.type==Constants.INTERVAL_ALIGNMENT && (containedInShortPeriod||containedInLongPeriod))
													   )
													 ) ||
													 (!periodicOrDense && nodeTo.type==Constants.INTERVAL_DENSE_SUBSTRING)
												   ) {
													if ( ( periodicOrDense &&
													       ( (nodeTo.type==Constants.INTERVAL_PERIODIC && nodeTo.hasLongPeriod==(type1Prime==2)) ||
														     (nodeTo.type==Constants.INTERVAL_ALIGNMENT && ((containedInShortPeriod&&(type1Prime==1))||(containedInLongPeriod&&(type1Prime==2))))
														   )
														 ) ||
														 !periodicOrDense
													   ) {
														nodeID2=to;
														pair.node2=nodesArray[nodeID2];
														pair.alignmentStart2=nodeTo.start+tmp[0];
														pair.alignmentEnd2=nodeTo.start+tmp[1];
													}
												}
					   							else if (periodicOrDense || nodeTo.type!=Constants.INTERVAL_PERIODIC) {
													type2=nodeTo.inPeriodicSubstring(tmp[0],tmp[1]);
													if (type2==0 && periodicOrDense) type2Prime=nodeTo.intersectsPeriodicSubstring(tmp[0],tmp[1]);
													else type2Prime=type2;
													if (type2Prime==type1Prime) {
														nodeID2=getPeriodicInterval(nodeTo.start+tmp[0],nodeTo.start+tmp[1],nodeTo.read,to,periodicOrDense,type2Prime);
														if (nodeID2>=0) {
															period=periodicOrDense?nodesArray[nodeID2].period:-1;
															if (type2Prime==2 && period>0 && tmp[1]-tmp[0]+1<period) nodeID2=-1;
															else {
																pair.node2=nodesArray[nodeID2];
																pair.alignmentStart2=nodeTo.start+tmp[0];
																pair.alignmentEnd2=nodeTo.start+tmp[1];
															}
														}
													}
												}
											}
										}
									}
								}
							}
							if ( nodeID1!=-1 && nodeID2!=-1 &&
		   						 ( Intervals.isApproximatelyContained(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end) ||  
		   					   	   Intervals.areApproximatelyIdentical(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end)
		   					 	 ) &&
		   				         ( Intervals.isApproximatelyContained(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end) ||  
		   					   	   Intervals.areApproximatelyIdentical(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end)
		   					 	 )
							   ) {
								pair.trim();
								if (pair.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && pair.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
									pair.orientation=orientation==0;
									pair.diffs=avgDiffs;
									fixPeriodicUndersplits_addEdge(pair,edge,i,to);
									nNewEdgesPrime++;
								}
								else {
									// NOP. The original edge might be disconnected as a
									// result.
								}
								nNewEdges++;
								edge.removeOverlap(k);
							}
						}
					}
					if (edge.overlap==0) edge.overlap=-1;
				}
				
				// Shared substring sub-edges
				if (edge.sharedSubstring!=-1) {
					nodeID1=-1; nodeID2=-1;
					projection=edge.projectOnto_sharedSubstring(i,tmp,-1);  // Projecting one of the possibly many shared substring intervals (only this one has been stored).	
					if (projection!=-1) {
						type1=node.inPeriodicSubstring(tmp[0],tmp[1]);
						if (type1==0 && periodicOrDense) type1Prime=node.intersectsPeriodicSubstring(tmp[0],tmp[1]);
						else type1Prime=type1;
						if (type1Prime>0) {
							nodeID1=getPeriodicInterval(node.start+tmp[0],node.start+tmp[1],read,i,periodicOrDense,type1Prime);
							if (nodeID1>=0) {
								period=periodicOrDense?nodesArray[nodeID1].period:-1;
								if (type1Prime==1 || (type1Prime==2 && (period<0 || tmp[1]-tmp[0]+1>=period))) {
									pair.node1=nodesArray[nodeID1];
									pair.alignmentStart1=node.start+tmp[0];
									pair.alignmentEnd1=node.start+tmp[1];
									projection=edge.projectOnto_sharedSubstring(to,tmp,-1);
									if (projection!=-1) {
										containedInShortPeriod=(nodeTo.isContainedInterval&(1<<Constants.INTERVAL_PERIODIC))!=0;
										containedInLongPeriod=(nodeTo.isContainedInterval&(1<<(Constants.INTERVAL_PERIODIC+1)))!=0;
										if ( ( periodicOrDense && 
											   ( nodeTo.type==Constants.INTERVAL_PERIODIC ||
												 (nodeTo.type==Constants.INTERVAL_ALIGNMENT && (containedInShortPeriod||containedInLongPeriod))
											   )
											 ) ||
											 (!periodicOrDense && nodeTo.type==Constants.INTERVAL_DENSE_SUBSTRING)
										   ) {
											if ( ( periodicOrDense &&
											       ( nodeTo.type==Constants.INTERVAL_PERIODIC ||
												     (nodeTo.type==Constants.INTERVAL_ALIGNMENT && ((containedInShortPeriod&&(type1Prime==1))||(containedInLongPeriod&&(type1Prime==2))))
												   )
												 ) ||
												 !periodicOrDense
											   ) {
												nodeID2=to;
												pair.node2=nodesArray[nodeID2];
												pair.alignmentStart2=nodeTo.start+tmp[0];
												pair.alignmentEnd2=nodeTo.start+tmp[1];
											}
										}
			   							else if (periodicOrDense || nodeTo.type!=Constants.INTERVAL_PERIODIC) {
											type2=nodeTo.inPeriodicSubstring(tmp[0],tmp[1]);
											if (type2==0 && periodicOrDense) type2Prime=nodeTo.intersectsPeriodicSubstring(tmp[0],tmp[1]);
											else type2Prime=type2;
											if (type2Prime==type1Prime) {
												nodeID2=getPeriodicInterval(nodeTo.start+tmp[0],nodeTo.start+tmp[1],nodeTo.read,to,periodicOrDense,type2Prime);
												if (nodeID2>=0) {
													period=periodicOrDense?nodesArray[nodeID2].period:-1;
													if (type2Prime==2 && period>0 && tmp[1]-tmp[0]+1<period) nodeID2=-1;
													else {
														pair.node2=nodesArray[nodeID2];
														pair.alignmentStart2=nodeTo.start+tmp[0];
														pair.alignmentEnd2=nodeTo.start+tmp[1];
													}
												}
											}
										}
									}
								}
							}
						}
					}
					if ( nodeID1!=-1 && nodeID2!=-1 &&
						 ( Intervals.isApproximatelyContained(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end) ||  
					   	   Intervals.areApproximatelyIdentical(pair.alignmentStart1,pair.alignmentEnd1,pair.node1.start,pair.node1.end)
					 	 ) &&
				         ( Intervals.isApproximatelyContained(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end) ||  
					   	   Intervals.areApproximatelyIdentical(pair.alignmentStart2,pair.alignmentEnd2,pair.node2.start,pair.node2.end)
					 	 )
					   ) {
						pair.trim();
						if (pair.alignmentLength1()>=MIN_ALIGNMENT_LENGTH && pair.alignmentLength2()>=MIN_ALIGNMENT_LENGTH) {
							pair.orientation=orientation==0;
							pair.diffs=avgDiffs;
							fixPeriodicUndersplits_addEdge(pair,edge,i,to);
							nNewEdgesPrime++;
						}
						else {
							// NOP. The original edge might be disconnected as a result.
						}
						nNewEdges++;
						if (nSharedSubstrings==1) edge.sharedSubstring=-1;
					}
				}
				
				// Removing the current edge iff all edges with stored coordinates have 
				// been used. This is quite aggressive; one could instead remove the edge
				// iff all its alignments or step 2 types have been used (an edge might
				// have no step 2 type but some alignments, since e.g. it might have
				// multiple containment alignments, but the coordinates of just one of
				// them are stored in the edge).
				if (nNewEdges>0) {
					decomposedEdges++;
					newEdges+=nNewEdgesPrime;
					if (nNewEdges==edge.getNAlignments() || nNewEdges==nStoredAlignments || edge.nDistinctTypes()==0) {
						neighbors[i][j]=null;
						for (k=0; k<nNeighbors[to]; k++) {
							if (neighbors[to][k]==edge) {
								neighbors[to][k]=null;
								break;
							}
						}
					}
					else {
						decomposedEdgesKept++;
						if (orientation==0) edge.nAlignmentsForward-=nNewEdges;
						else if (orientation==1) edge.nAlignmentsBackward-=nNewEdges;
						else {
							edge.nAlignmentsForward=edge.getNAlignments()-nNewEdges;  // Arbitrarily setting forward rather than backward alignments.
							edge.nAlignmentsBackward=0;
							edge.orientation=-1;
						}
						edge.avgDiffs=avgDiffs*edge.getNAlignments();
					}
				}
				edge.printed=1;
			}
		}
		
		// Deleting all edges $(u,v)$ such that $v$ is the only periodic neighbor of $u$.
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].type==Constants.INTERVAL_PERIODIC || nNeighbors[i]<1+MIN_NONPERIODIC_NEIGHBORS) continue;
			nPeriodicNeighbors=0; k=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]!=null && nodesArray[neighbors[i][j].getTo(i)].type==Constants.INTERVAL_PERIODIC) {
					nPeriodicNeighbors++;
					k=j;
				}
			}
			if (nPeriodicNeighbors!=1) continue;
			edge=neighbors[i][k]; to=edge.getTo(i);
			neighbors[i][k]=null;
			for (j=0; j<nNeighbors[to]; j++) {
				if (neighbors[to][j]==edge) {
					neighbors[to][j]=null;
					break;
				}
			}
		}
		
		// Removing $null$ neighbors
		for (i=0; i<nNodes; i++) {
			k=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				neighbors[i][++k]=neighbors[i][j];
			}
			nNeighbors[i]=k+1;
		}
		
		System.err.println((periodicOrDense?"Periodic":"Substring")+" undersplit resolution: "+decomposedEdges+" edges decomposed; "+decomposedEdgesKept+" decomposed edges kept; "+newEdges+" new edges created.");
	}
	
	
	/**
	 * Removes all edges, between a short-period interval and a nonperiodic interval, that
	 * involve just a proper substring of the nonperiodic interval. Assuming that all 
	 * fixes related to periodic intervals have already been applied, such edges are more 
	 * likely to connect different clusters in the graph than to be contained in a single 
	 * cluster.
	 */
	private static final void removePeriodicNonperiodicEdges() {
		boolean found;
		int i, j, k;
		int to, projection;
		Node node, nodeTo;
		Edge edge;
		
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if (node.type!=Constants.INTERVAL_PERIODIC || node.hasLongPeriod) continue;
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (edge==null) continue;
				to=edge.getTo(i); nodeTo=nodesArray[to];
				if (nodeTo.type==Constants.INTERVAL_PERIODIC) continue;				
				projection=edge.projectOnto(to,tmpArray,true);
				found=false;
				for (k=0; k<=projection; k+=2) {
					if (Intervals.areApproximatelyIdentical(tmpArray[k],tmpArray[k+1],0,nodeTo.length()-1)) {
						found=true;
						break;
					}
				}
				if (found) continue;
				neighbors[i][j]=null;
				for (k=0; k<nNeighbors[to]; k++) {
					if (neighbors[to][k]==edge) {
						neighbors[to][k]=null;
						break;
					}
				}
			}
		}
		
		// Removing $null$ neighbors
		for (i=0; i<nNodes; i++) {
			k=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				neighbors[i][++k]=neighbors[i][j];
			}
			nNeighbors[i]=k+1;
		}
	}
	
	
	/**
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$.
	 *
	 * @param pos search in $nodesArray$ starts from $pos$ (excluded) and proceeds 
	 * backwards and forwards;
	 * @param periodicOrDense TRUE=periodic; FALSE=dense substring of substring type;
	 * @param hasLongPeriod 1=false, 2=true;
	 * @return a smallest periodic interval of $read$ that approximately contains 
	 * $[start..end]$ or is identical to it. If no such interval can be found, returns a 
	 * periodic interval with largest intersection with $[start..end]$. Returns -1 if no 
	 * periodic interval intersects $[start..end]$.
	 */
	private static final int getPeriodicInterval(int start, int end, int read, int pos, boolean periodicOrDense, int hasLongPeriod) {
		final boolean HLP = hasLongPeriod==2;
		int i, minLength, minInterval, maxInterval;
		int intersection, maxIntersection;
		Node tmpNode;
		
		minInterval=-1; minLength=Math.POSITIVE_INFINITY;
		maxInterval=-1; maxIntersection=-1;
		for (i=pos-1; i>=0; i--) {
			tmpNode=nodesArray[i];
			if (tmpNode.read!=read) break;
			if ( (periodicOrDense && (tmpNode.type!=Constants.INTERVAL_PERIODIC || tmpNode.hasLongPeriod!=HLP)) ||
				 (!periodicOrDense && tmpNode.type!=Constants.INTERVAL_DENSE_SUBSTRING)
			   ) continue;
			if ( Intervals.areApproximatelyIdentical(start,end,tmpNode.start,tmpNode.end) ||
				 Intervals.isApproximatelyContained(start,end,tmpNode.start,tmpNode.end)
			   ) {
				   if (tmpNode.length()<minLength) minLength=tmpNode.length();
				   minInterval=i;
			}
			if (periodicOrDense) {
				intersection=Intervals.intersectionLength(start,end,tmpNode.start,tmpNode.end);
				if (intersection>maxIntersection) {
					maxIntersection=intersection;
					maxInterval=i;
				}
			}
		}
		for (i=pos+1; i<nNodes; i++) {
			tmpNode=nodesArray[i];
			if (tmpNode.read!=read || tmpNode.start>=end) break;
			if ( (periodicOrDense && (tmpNode.type!=Constants.INTERVAL_PERIODIC || tmpNode.hasLongPeriod!=HLP)) ||
				 (!periodicOrDense && tmpNode.type!=Constants.INTERVAL_DENSE_SUBSTRING)
			   ) continue;
			if ( Intervals.areApproximatelyIdentical(start,end,tmpNode.start,tmpNode.end) ||
				 Intervals.isApproximatelyContained(start,end,tmpNode.start,tmpNode.end)
			   ) {
				   if (tmpNode.length()<minLength) minLength=tmpNode.length();
				   minInterval=i;
			}
			if (periodicOrDense) {
				intersection=Intervals.intersectionLength(start,end,tmpNode.start,tmpNode.end);
				if (intersection>maxIntersection) {
					maxIntersection=intersection;
					maxInterval=i;
				}
			}
		}
		if (minInterval!=-1) return minInterval;
		return maxInterval;
	}

	
	/**
	 * Sets the $periodicSubstrings$ list of every node based on its current neighbors.
	 * Any projection is allowed.
	 *
	 * @param add TRUE: adds substrings to the existing $periodicSubstrings$ list, rather 
	 * than resetting the list;
	 * @param tmp temporary space, of size at least two.
	 */
	public static final void setPeriodicSubstringsOfNodes_neighbors(boolean add, int[] tmp) {
		int i, j, k;
		int read, last;
		Node node, otherNode;
		
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if (!add) node.lastPeriodicSubstring=-1;
			if (node.type==Constants.INTERVAL_PERIODIC) continue;
			read=node.read;
			for (j=0; j<nNeighbors[i]; j++) {
				otherNode=nodesArray[neighbors[i][j].getTo(i)];
				if (otherNode.read==read || otherNode.type!=Constants.INTERVAL_PERIODIC) continue;
				last=neighbors[i][j].projectOnto(i,tmp,true);
				for (k=0; k<=last; k+=2) node.addPeriodicSubstring(tmp[k],tmp[k+1],otherNode.hasLongPeriod);
			}
		}
	}
	
	
	/**
	 * Let two periodic substrings of a non-periodic node be equivalent if they intersect 
	 * or are at most some distance apart. The procedure computes the connected components
	 * of such relation, and takes the union of all periodic substrings in the same 
	 * component.
	 *
	 * Remark: periodic substrings that are contained in one another were already merged 
	 * by $Node.addPeriodicSubstring()$.
	 *
	 * Remark: the procedure can reallocate $tmpBoolean$.
	 *
	 * @return max number of periodic substrings in a node after the procedure completes.
	 */
	public static final int mergePeriodicSubstringsOfNodes() {
		final int DISTANCE_THRESHOLD = IO.quantum;
		int i, j, k;
		int top, current, newStart, newEnd, nPeriodicSubstrings, maxPeriodicSubstrings;
		int nNodesWithPeriodicSubstrings, nNodesMerged, type;
		Node node;
		
		if (stack==null || stack.length<nNodes) stack = new int[nNodes];
		nNodesWithPeriodicSubstrings=0; nNodesMerged=0; maxPeriodicSubstrings=0; current=-1;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if (node.lastPeriodicSubstring==-1) continue;
			nNodesWithPeriodicSubstrings++;
			if (node.lastPeriodicSubstring==2) {
				if (1>maxPeriodicSubstrings) maxPeriodicSubstrings=1;				
				continue;
			}
			if (tmpBoolean.length<node.lastPeriodicSubstring+1) tmpBoolean = new boolean[node.lastPeriodicSubstring+1];
			Math.set(tmpBoolean,node.lastPeriodicSubstring,false);
			for (j=0; j<=node.lastPeriodicSubstring; j+=3) {
				if (tmpBoolean[j]) continue;
				stack[0]=j; top=0;
				newStart=node.periodicSubstrings[j];
				newEnd=node.periodicSubstrings[j+1];
				type=node.periodicSubstrings[j+2];
				while (top>=0) {
					current=stack[top--];
					for (k=0; k<=node.lastPeriodicSubstring; k+=3) {
						if (k==j || k==current || tmpBoolean[k] || node.periodicSubstrings[k+2]!=type) continue;
						if (node.periodicSubstrings[k+1]<node.periodicSubstrings[current]-DISTANCE_THRESHOLD || node.periodicSubstrings[k]>node.periodicSubstrings[current+1]+DISTANCE_THRESHOLD) continue;
						newStart=Math.min(newStart,node.periodicSubstrings[k]);
						newEnd=Math.max(newEnd,node.periodicSubstrings[k+1]);
						tmpBoolean[k]=true;
						stack[++top]=k;
					}
				}
				node.periodicSubstrings[j]=newStart;
				node.periodicSubstrings[j+1]=newEnd;
			}
			k=-1;
			for (j=0; j<=node.lastPeriodicSubstring; j+=3) {
				if (tmpBoolean[j]) continue;
				node.periodicSubstrings[++k]=node.periodicSubstrings[j];
				node.periodicSubstrings[++k]=node.periodicSubstrings[j+1];
				node.periodicSubstrings[++k]=node.periodicSubstrings[j+2];
			}
			if (k!=node.lastPeriodicSubstring) {
				node.lastPeriodicSubstring=k;
				nNodesMerged++;
			}
			nPeriodicSubstrings=(k+1)/3;
			if (nPeriodicSubstrings>maxPeriodicSubstrings) maxPeriodicSubstrings=nPeriodicSubstrings;			
		}
		System.err.println("Nodes with periodic substrings: "+nNodesWithPeriodicSubstrings+" ("+IO.getPercent(nNodesWithPeriodicSubstrings,nNodes)+"% of all nodes)");
		System.err.println("Periodic substrings merged in nodes: "+nNodesMerged+" ("+IO.getPercent(nNodesMerged,nNodesWithPeriodicSubstrings)+"%)");
		System.err.println("Max. merged periodic substrings in a node: "+maxPeriodicSubstrings);
		return maxPeriodicSubstrings;
	}	
	
	
	/**
	 * Stores in the $periodicSubstrings$ list of every nonperiodic and non-dense-
	 * substring node, the set of substrings induced by the fact that the node straddles 
	 * with or contains a dense substring of substring type in the same read.
	 *
	 * Remark: the procedure stores both weak and non-weak dense substrings. Weak dense
	 * substrings would probably be enough, since an alignment cointained in a non-weak 
	 * dense substring of substring type is not likely to be assigned to a containing 
	 * interval.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$.
	 * The procedure does not create new intervals.
	 *
	 * @param add TRUE: adds substrings to the existing $periodicSubstrings$ list, rather 
	 * than resetting the list.
	 */
	private static final void setDenseSubstringsOfNodes_sameRead(boolean add) {
		boolean added;
		int i, j;
		int read, start, end, count;

		count=0;
		for (i=0; i<nNodes; i++) {
			if (!add) nodesArray[i].lastPeriodicSubstring=-1;
			if (nodesArray[i].type==Constants.INTERVAL_DENSE_SUBSTRING || nodesArray[i].type==Constants.INTERVAL_PERIODIC) continue;
			added=false;
			read=nodesArray[i].read; start=nodesArray[i].start; end=nodesArray[i].end;
			for (j=i-1; j>=0; j--) {
				if (nodesArray[j].read!=read) break;
				if (nodesArray[j].type!=Constants.INTERVAL_DENSE_SUBSTRING) continue;
				if ( ( !Intervals.areApproximatelyIdentical(start,end,nodesArray[j].start,nodesArray[j].end) &&
					   !Intervals.isApproximatelyContained(start,end,nodesArray[j].start,nodesArray[j].end) &&
					   ( Intervals.straddles(start,end,nodesArray[j].start,nodesArray[j].end) ||
				         Intervals.isApproximatelyContained(nodesArray[j].start,nodesArray[j].end,start,end)
					   )
					 ) ||
					 Intervals.isApproximatelyContained(start,end,nodesArray[j].start,nodesArray[j].end)
				   ) {
					nodesArray[i].addPeriodicSubstring(Math.max(start,nodesArray[j].start)-start,Math.min(end,nodesArray[j].end)-start,false);
					added=true;
				}
			}
			for (j=i+1; j<nNodes; j++) {
				if (nodesArray[j].read!=read || nodesArray[j].start>=end) break;
				if (nodesArray[j].type!=Constants.INTERVAL_DENSE_SUBSTRING) continue;
				if ( ( !Intervals.areApproximatelyIdentical(start,end,nodesArray[j].start,nodesArray[j].end) &&
					   !Intervals.isApproximatelyContained(start,end,nodesArray[j].start,nodesArray[j].end) &&
					   ( Intervals.straddles(start,end,nodesArray[j].start,nodesArray[j].end) ||
				         Intervals.isApproximatelyContained(nodesArray[j].start,nodesArray[j].end,start,end)
					   )
					 ) ||
					 Intervals.isApproximatelyContained(start,end,nodesArray[j].start,nodesArray[j].end)
				   ) {
					nodesArray[i].addPeriodicSubstring(Math.max(start,nodesArray[j].start)-start,Math.min(end,nodesArray[j].end)-start,false);
					added=true;
				}
			}
			if (added) count++;
		}
		System.err.println(count+" ("+IO.getPercent(count,nNodes)+"%) non-substring-type nodes have a same-read dense substring of substring type");
	}
	
	
	/**
	 * Adds a new edge that corresponds to exactly one alignment of $source_edge$.
	 *
	 * @param pair $trim()$ is assumed to have already been called;
	 * @param source_nodeX the ID of the node of $source_edge$ that corresponds to
	 * $pair.nodeX$.
	 */
	private static final void fixPeriodicUndersplits_addEdge(AlignmentPair pair, Edge source_edge, int source_node1, int source_node2) {
		int i, j, nodeID1;
		Node target_node1, target_node2;
		Edge newEdge;
		
		target_node1=pair.node1; target_node2=pair.node2;
		pair.inCanonicalForm();
		//pair.trim();  // Assumed to have already been called
		newEdge = new Edge(pair);
		transferFlexibleOverhangs(source_edge,newEdge,source_node1,target_node1.nodeID,source_node2,target_node2.nodeID);
		newEdge.printed=1;  // To avoid processing by $fixPeriodicUndersplits()$.
		nodeID1=newEdge.nodeID1;
		j=findEdge(nodeID1,newEdge,false);
		if (j==Math.NEGATIVE_INFINITY) {
			addEdge(nodeID1,newEdge);
			addEdge(newEdge.nodeID2,newEdge);
		}
		else {
			neighbors[nodeID1][j].addTypes(newEdge);
			neighbors[nodeID1][j].addDiffs(newEdge);
			neighbors[nodeID1][j].addOrientation(newEdge);
		}
	}
	
	
	/**
	 * Translates flexibility information from $source_edge$ to $target_edge$.
	 *
	 * @param sourceNodeX,targetNodeX a map from the ID of a node in $source_edge$ to the
	 * ID of a node in $target_edge$.
	 */
	private static final void transferFlexibleOverhangs(Edge source_edge, Edge target_edge, int source_nodeA, int target_nodeA, int source_nodeB, int target_nodeB) {
		int i, j;
		int jSource, newOverhangA, newOverhangB;

		if (target_edge.overlap==-1 || source_edge.overlap==-1 || source_edge.flexibleOverlap==0 || source_edge.maxOverhangs==null) return;
		if (target_edge.maxOverhangs==null) {
			target_edge.maxOverhangs = new int[4][2];
			for (i=0; i<4; i++) {
				target_edge.maxOverhangs[i][0]=target_edge.overhangs[i][0];
				target_edge.maxOverhangs[i][1]=target_edge.overhangs[i][1];
			}
		}
		j=0;
		for (i=Constants.OVERLAP_PREFIX_PREFIX; i<=Constants.OVERLAP_SUFFIX_SUFFIX; i<<=1) {
			if ((target_edge.overlap&i)==0) {
				j++;
				continue;
			}
			if (target_edge.nodeID1==target_nodeA) {
				if (source_nodeA==source_edge.nodeID1) jSource=j;
				else jSource=Constants.oppositeDirectionOverlap(j);
			}
			else {
				if (source_nodeB==source_edge.nodeID1) jSource=j;
				else jSource=Constants.oppositeDirectionOverlap(j);
			}
			if ((source_edge.flexibleOverlap&jSource)==0) {
				j++;
				continue;
			}
			if (source_edge.nodeID1==source_nodeA) {
				newOverhangA=nodesArray[target_nodeA].length()-(nodesArray[source_nodeA].length()-source_edge.maxOverhangs[jSource][0]);
				newOverhangB=nodesArray[target_nodeB].length()-(nodesArray[source_nodeB].length()-source_edge.maxOverhangs[jSource][1]);
			}
			else {
				newOverhangA=nodesArray[target_nodeA].length()-(nodesArray[source_nodeA].length()-source_edge.maxOverhangs[jSource][1]);
				newOverhangB=nodesArray[target_nodeB].length()-(nodesArray[source_nodeB].length()-source_edge.maxOverhangs[jSource][0]);
			}
			if (target_edge.nodeID1==target_nodeA) {
				target_edge.maxOverhangs[j][0]=newOverhangA;
				target_edge.maxOverhangs[j][1]=newOverhangB;
			}
			else {
				target_edge.maxOverhangs[j][0]=newOverhangB;
				target_edge.maxOverhangs[j][1]=newOverhangA;
			}
			if (target_edge.maxOverhangs[j][0]!=target_edge.overhangs[j][0] || target_edge.maxOverhangs[j][1]!=target_edge.overhangs[j][1]) target_edge.flexibleOverlap|=j;
			j++;
		}
	}
	
	
	/**
	 * Marks as periodic all nonperiodic nodes that belong to an overlap loop with two 
	 * nodes (one of which might already be periodic), and updates all incident edges
	 * (the procedure considers overlap loops with two nodes just for simplicity).
	 *
	 * Then, the procedure marks some intervals as periodic, using:
	 * - edges that contain both overlap and shared substring alignments (see 
	 * $fixOverlapSharedSubstringEdges_*$ procedures for details);
	 * - edges with just one overlap type, but with multiple alignments, that cover most
	 * of the surface of an interval (these edges are also used in $fixMultiEdges()$).
	 *
	 * Remark: the procedure could be executed multiple times, until no new node is marked
	 * as periodic.
	 *
	 * Remark: the procedure does not assume $nodesArray$ to be sorted, and it does not 
	 * create new intervals.
	 */
	private static final void fixPeriodicLoops() {
		final int THRESHOLD = IO.quantum;
		final double JACCARD_THRESHOLD = 0.2;  // Arbitrary
		boolean hasLongPeriod;
		int i, j;
		int p, q, pPrime, qPrime, x, y, xPrime, yPrime, xDoublePrime, yDoublePrime;
		int lastLoop, previousLastLoop, from, to, period, nLoops, nTransformedEdges;
		int row, length;
		Node otherNode;
		Edge edge;
		int[] loops = stack;  // Reusing the shared temporary area
		int[] stats = new int[6];  // Indexed by $Constants.INTERVAL_*$ values.
		
		nLoops=0; lastLoop=-1;
		
		// Collecting nodes to be transformed to periodic, from overlap loops.
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].type!=Constants.INTERVAL_PERIODIC) nodesArray[i].period=-1;
		}
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=-1;  // Using this temporary field to mark processed edges.
		}
		Math.set(stats,stats.length-1,0);
		for (i=0; i<nNodes; i++) {
			if (i%100000==0) System.err.println("fixPeriodicLoops> "+i+" nodes done ("+IO.getPercent(i,nNodes)+"%)");			
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (edge.printed!=-1 || edge.overlap==-1 || (nodesArray[edge.nodeID1].type==Constants.INTERVAL_PERIODIC && nodesArray[edge.nodeID2].type==Constants.INTERVAL_PERIODIC)) continue;
				edge.printed=1;
				if ((edge.overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && (edge.overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) {
					nLoops++;
					p=edge.overhangs[Constants.overlap2id(Constants.OVERLAP_SUFFIX_PREFIX)][0];
 				   	pPrime=nodesArray[edge.nodeID1].length()-edge.overhangs[Constants.overlap2id(Constants.OVERLAP_PREFIX_SUFFIX)][0];
				   	if (nodesArray[edge.nodeID2].type==Constants.INTERVAL_PERIODIC) {
						if (p<=pPrime+THRESHOLD) {
							if (nodesArray[edge.nodeID1].period==-1) {
								nodesArray[edge.nodeID1].hasLongPeriod=nodesArray[edge.nodeID2].hasLongPeriod;
								nodesArray[edge.nodeID1].period=Math.max(0,nodesArray[edge.nodeID2].period);
							}
							else {
								if (!nodesArray[edge.nodeID2].hasLongPeriod) nodesArray[edge.nodeID1].hasLongPeriod=false;
								if (nodesArray[edge.nodeID1].period==0) nodesArray[edge.nodeID1].period=Math.max(0,nodesArray[edge.nodeID2].period);
								else if (nodesArray[edge.nodeID2].period>0) nodesArray[edge.nodeID1].period=Math.min(nodesArray[edge.nodeID1].period,nodesArray[edge.nodeID2].period);
							}
							loops[++lastLoop]=edge.nodeID1;
							stats[nodesArray[edge.nodeID1].type]++;
						}
						continue;
					}
					q=edge.overhangs[Constants.overlap2id(Constants.OVERLAP_PREFIX_SUFFIX)][1];
					qPrime=nodesArray[edge.nodeID2].length()-edge.overhangs[Constants.overlap2id(Constants.OVERLAP_SUFFIX_PREFIX)][1];
					if (nodesArray[edge.nodeID1].type==Constants.INTERVAL_PERIODIC) {
						if (q<=qPrime+THRESHOLD) {
							if (nodesArray[edge.nodeID2].period==-1) {
								nodesArray[edge.nodeID2].hasLongPeriod=nodesArray[edge.nodeID1].hasLongPeriod;
								nodesArray[edge.nodeID2].period=Math.max(0,nodesArray[edge.nodeID1].period);
							}
							else {
								if (!nodesArray[edge.nodeID1].hasLongPeriod) nodesArray[edge.nodeID2].hasLongPeriod=false;
								if (nodesArray[edge.nodeID2].period==0) nodesArray[edge.nodeID2].period=Math.max(0,nodesArray[edge.nodeID1].period);
								else if (nodesArray[edge.nodeID1].period>0) nodesArray[edge.nodeID2].period=Math.min(nodesArray[edge.nodeID2].period,nodesArray[edge.nodeID1].period);
							}
							loops[++lastLoop]=edge.nodeID2;
							stats[nodesArray[edge.nodeID2].type]++;
						}
						continue;
					}
					if (pPrime>=p && qPrime>=q && Math.abs(p,qPrime-q)<=THRESHOLD && Math.abs(pPrime,nodesArray[edge.nodeID1].length()-(qPrime-q))<=THRESHOLD && Math.abs(nodesArray[edge.nodeID2].length()-(pPrime-p),qPrime)<=THRESHOLD) {
						period=fixPeriodicLoops_getPeriod(edge,nodesArray[edge.nodeID1],nodesArray[edge.nodeID2],true);
						if (nodesArray[edge.nodeID1].period==-1) {
							nodesArray[edge.nodeID1].hasLongPeriod=true;
							nodesArray[edge.nodeID1].period=period;
						}
						else if (nodesArray[edge.nodeID1].hasLongPeriod) nodesArray[edge.nodeID1].period=Math.min(nodesArray[edge.nodeID1].period,period);
						loops[++lastLoop]=edge.nodeID1;
						stats[nodesArray[edge.nodeID1].type]++;
						if (nodesArray[edge.nodeID2].period==-1) {
							nodesArray[edge.nodeID2].hasLongPeriod=true;
							nodesArray[edge.nodeID2].period=period;
						}
						else if (nodesArray[edge.nodeID2].hasLongPeriod) nodesArray[edge.nodeID2].period=Math.min(nodesArray[edge.nodeID2].period,period);
						loops[++lastLoop]=edge.nodeID2;
						stats[nodesArray[edge.nodeID2].type]++;
					}
				}
				else if ((edge.overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0 && (edge.overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) {
					nLoops++;
					p=edge.overhangs[Constants.overlap2id(Constants.OVERLAP_SUFFIX_SUFFIX)][0];
					pPrime=nodesArray[edge.nodeID1].length()-edge.overhangs[Constants.overlap2id(Constants.OVERLAP_PREFIX_PREFIX)][0];
					if (nodesArray[edge.nodeID2].type==Constants.INTERVAL_PERIODIC) {
						if (p<=pPrime+THRESHOLD) {
							if (nodesArray[edge.nodeID1].period==-1) {
								nodesArray[edge.nodeID1].hasLongPeriod=nodesArray[edge.nodeID2].hasLongPeriod;
								nodesArray[edge.nodeID1].period=Math.max(0,nodesArray[edge.nodeID2].period);
							}
							else {
								if (!nodesArray[edge.nodeID2].hasLongPeriod) nodesArray[edge.nodeID1].hasLongPeriod=false;
								if (nodesArray[edge.nodeID1].period==0) nodesArray[edge.nodeID1].period=Math.max(0,nodesArray[edge.nodeID2].period);
								else if (nodesArray[edge.nodeID2].period>0) nodesArray[edge.nodeID1].period=Math.min(nodesArray[edge.nodeID1].period,nodesArray[edge.nodeID2].period);
							}
							loops[++lastLoop]=edge.nodeID1;
							stats[nodesArray[edge.nodeID1].type]++;
						}
						continue;
					}
					q=nodesArray[edge.nodeID2].length()-edge.overhangs[Constants.overlap2id(Constants.OVERLAP_PREFIX_PREFIX)][1];
					qPrime=edge.overhangs[Constants.overlap2id(Constants.OVERLAP_SUFFIX_SUFFIX)][1];
					if (nodesArray[edge.nodeID1].type==Constants.INTERVAL_PERIODIC) {
						if (q<=qPrime+THRESHOLD) {
							if (nodesArray[edge.nodeID2].period==-1) {
								nodesArray[edge.nodeID2].hasLongPeriod=nodesArray[edge.nodeID1].hasLongPeriod;
								nodesArray[edge.nodeID2].period=Math.max(0,nodesArray[edge.nodeID1].period);
							}
							else {
								if (!nodesArray[edge.nodeID1].hasLongPeriod) nodesArray[edge.nodeID2].hasLongPeriod=false;
								if (nodesArray[edge.nodeID2].period==0) nodesArray[edge.nodeID2].period=Math.max(0,nodesArray[edge.nodeID1].period);
								else if (nodesArray[edge.nodeID1].period>0) nodesArray[edge.nodeID2].period=Math.min(nodesArray[edge.nodeID2].period,nodesArray[edge.nodeID1].period);
							}
							loops[++lastLoop]=edge.nodeID2;
							stats[nodesArray[edge.nodeID2].type]++;
						}
						continue;
					}
					if (pPrime>=p && q>=qPrime && Math.abs(p,q-qPrime)<=THRESHOLD && Math.abs(nodesArray[edge.nodeID1].length()-(q-qPrime),pPrime)<=THRESHOLD && Math.abs(pPrime-p,qPrime)<=THRESHOLD) {
						period=fixPeriodicLoops_getPeriod(edge,nodesArray[edge.nodeID1],nodesArray[edge.nodeID2],false);
						if (nodesArray[edge.nodeID1].period==-1) {
							nodesArray[edge.nodeID1].hasLongPeriod=true;
							nodesArray[edge.nodeID1].period=period;
						}
						else if (nodesArray[edge.nodeID1].hasLongPeriod) nodesArray[edge.nodeID1].period=Math.min(nodesArray[edge.nodeID1].period,period);
						loops[++lastLoop]=edge.nodeID1;
						stats[nodesArray[edge.nodeID1].type]++;
						if (nodesArray[edge.nodeID2].period==-1) {
							nodesArray[edge.nodeID2].hasLongPeriod=true;
							nodesArray[edge.nodeID2].period=period;
						}
						else if (nodesArray[edge.nodeID2].hasLongPeriod) nodesArray[edge.nodeID2].period=Math.min(nodesArray[edge.nodeID2].period,period);
						loops[++lastLoop]=edge.nodeID2;
						stats[nodesArray[edge.nodeID2].type]++;
					}
				}
			}
		}
		System.err.println("fixPeriodicLoops> "+nLoops+" overlap loops found");
		if (nLoops>=0) {
			System.err.println("fixPeriodicLoops> "+(lastLoop+1)+" distinct non-periodic nodes transformed to periodic");
			System.err.println("fixPeriodicLoops> transformations: alignment="+stats[0]+", prefix/suffix="+(stats[1]+stats[2])+", prefixSuffix="+stats[3]+", dense="+stats[4]+", singleDeletion="+stats[5]);
		}
		
		// Collecting nodes to be transformed to periodic, from overlap+sharedSubstring
		// edges.
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=-1;
		}
		previousLastLoop=lastLoop; Math.set(stats,stats.length-1,0);
		for (i=0; i<nNodes; i++) {
			if (i%100000==0) System.err.println("fixPeriodicLoops> "+i+" nodes done ("+IO.getPercent(i,nNodes)+"%)");			
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (edge.printed!=-1 || nodesArray[edge.nodeID1].type==Constants.INTERVAL_PERIODIC || nodesArray[edge.nodeID2].type==Constants.INTERVAL_PERIODIC) continue;
				edge.printed=1;
				if (edge.overlap!=-1 && edge.sharedSubstring!=-1 && edge.sharedSubstringOverhangs!=null) {
					if ((edge.overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) lastLoop=fixOverlapSharedSubstringEdges_suffixPrefix(edge.nodeID1,nodesArray[edge.nodeID1].length(),edge.sharedSubstringOverhangs[0],edge.sharedSubstringOverhangs[1],edge.overhangs[Constants.overlap2id(Constants.OVERLAP_SUFFIX_PREFIX)][0],edge.nodeID2,nodesArray[edge.nodeID2].length(),edge.sharedSubstringOverhangs[2],edge.sharedSubstringOverhangs[3],edge.overhangs[Constants.overlap2id(Constants.OVERLAP_SUFFIX_PREFIX)][1],loops,lastLoop,THRESHOLD,stats);
					if ((edge.overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) lastLoop=fixOverlapSharedSubstringEdges_suffixPrefix(edge.nodeID2,nodesArray[edge.nodeID2].length(),edge.sharedSubstringOverhangs[2],edge.sharedSubstringOverhangs[3],edge.overhangs[Constants.overlap2id(Constants.OVERLAP_PREFIX_SUFFIX)][1],edge.nodeID1,nodesArray[edge.nodeID1].length(),edge.sharedSubstringOverhangs[0],edge.sharedSubstringOverhangs[1],edge.overhangs[Constants.overlap2id(Constants.OVERLAP_PREFIX_SUFFIX)][0],loops,lastLoop,THRESHOLD,stats);
					if ((edge.overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) lastLoop=fixOverlapSharedSubstringEdges_suffixSuffix(edge.nodeID1,nodesArray[edge.nodeID1].length(),edge.sharedSubstringOverhangs[0],edge.sharedSubstringOverhangs[1],edge.overhangs[Constants.overlap2id(Constants.OVERLAP_SUFFIX_SUFFIX)][0],edge.nodeID2,nodesArray[edge.nodeID2].length(),edge.sharedSubstringOverhangs[2],edge.sharedSubstringOverhangs[3],edge.overhangs[Constants.overlap2id(Constants.OVERLAP_SUFFIX_SUFFIX)][1],loops,lastLoop,THRESHOLD,stats);
					if ((edge.overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) lastLoop=fixOverlapSharedSubstringEdges_prefixPrefix(edge.nodeID1,nodesArray[edge.nodeID1].length(),edge.sharedSubstringOverhangs[0],edge.sharedSubstringOverhangs[1],edge.overhangs[Constants.overlap2id(Constants.OVERLAP_PREFIX_PREFIX)][0],edge.nodeID2,nodesArray[edge.nodeID2].length(),edge.sharedSubstringOverhangs[2],edge.sharedSubstringOverhangs[3],edge.overhangs[Constants.overlap2id(Constants.OVERLAP_PREFIX_PREFIX)][1],loops,lastLoop,THRESHOLD,stats);
				}
			}
		}
		System.err.println("fixPeriodicLoops> "+(lastLoop-previousLastLoop)+" non-periodic nodes transformed to periodic using overlap+sharedSubstring edges");
		if (lastLoop-previousLastLoop>0) System.err.println("fixPeriodicLoops> transformations: alignment="+stats[0]+", prefix/suffix="+(stats[1]+stats[2])+", prefixSuffix="+stats[3]+", dense="+stats[4]+", singleDeletion="+stats[5]);
		
		// Collecting nodes to be transformed to periodic, from overlap multi-edges.
		previousLastLoop=lastLoop; Math.set(stats,stats.length-1,0);
		for (i=0; i<nNodes; i++) {
			if (i%100000==0) System.err.println("fixPeriodicLoops> "+i+" nodes done ("+IO.getPercent(i,nNodes)+"%)");
			if (nodesArray[i].type==Constants.INTERVAL_PERIODIC) continue;
			length=nodesArray[i].length();
			hasLongPeriod=true; period=Math.POSITIVE_INFINITY; 
			p=-1; q=length;  // Relative positions
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (!edge.hasOneType() || !edge.isType_overlap() || edge.getNAlignments()==1 || edge.overhangs==null || !edge.overhangMatricesDiffer()) continue;
				otherNode=nodesArray[edge.getTo(i)];
				if (otherNode.type==Constants.INTERVAL_PERIODIC && !otherNode.hasLongPeriod) hasLongPeriod=false;
				if (otherNode.period>0 && otherNode.period<period) period=otherNode.period;
				row=Constants.overlap2id(edge.overlap);
				if (i==edge.nodeID1) {
					if (edge.overlap==Constants.OVERLAP_SUFFIX_PREFIX || edge.overlap==Constants.OVERLAP_SUFFIX_SUFFIX) {
						x=edge.overhangs[row][0];
						if (x<q) q=x;
					}
					else {
						x=length-edge.overhangs[row][0]-1;
						if (x>p) p=x;
					}
				}
				else {
					if (edge.overlap==Constants.OVERLAP_PREFIX_SUFFIX || edge.overlap==Constants.OVERLAP_SUFFIX_SUFFIX) {
						x=edge.overhangs[row][1];
						if (x<q) q=x;
					}
					else {
						x=length-edge.overhangs[row][1]-1;
						if (x>p) p=x;
					}
				}
			}
			if ((q-p-1)<=JACCARD_THRESHOLD*length) {
				nodesArray[i].period=period==Math.POSITIVE_INFINITY?0:period;
				nodesArray[i].hasLongPeriod=hasLongPeriod;
				loops[++lastLoop]=i;
				stats[nodesArray[i].type]++;
			}
		}
		System.err.println("fixPeriodicLoops> "+(lastLoop-previousLastLoop)+" non-periodic nodes transformed to periodic using overlap multi-edges");
		if (lastLoop-previousLastLoop>0) System.err.println("fixPeriodicLoops> transformations: alignment="+stats[0]+", prefix/suffix="+(stats[1]+stats[2])+", prefixSuffix="+stats[3]+", dense="+stats[4]+", singleDeletion="+stats[5]);
		
		// Removing duplicates
		if (lastLoop>0) {
			if (lastLoop>0) Arrays.sort(loops,0,lastLoop+1);
			j=0;
			for (i=1; i<=lastLoop; i++) {
				if (loops[i]==loops[j]) continue;
				loops[++j]=loops[i];
			}
			lastLoop=j;
		}
		System.err.println("fixPeriodicLoops> "+(lastLoop+1)+" distinct non-periodic nodes transformed to periodic in total");
		
		// Transforming nodes and incident edges
		for (i=0; i<=lastLoop; i++) nodesArray[loops[i]].type=Constants.INTERVAL_PERIODIC;
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=-1;  // Using this temporary field to mark transformed edges.
		}
		nTransformedEdges=0;
		for (i=0; i<=lastLoop; i++) {
			from=loops[i];
			for (j=0; j<nNeighbors[from]; j++) {
				if (neighbors[from][j].printed!=-1) continue;
				neighbors[from][j].printed=1;
				nTransformedEdges+=IntervalGraphStep2.removeEdgeTypes(neighbors[from][j])?1:0;
			}
		}
		System.err.println("fixPeriodicLoops> "+nTransformedEdges+" edges transformed");
	}
	
	
	/**
	 * @param forwardOrReverse the overlap loop is in the forward (TRUE) or reverse-
	 * complement (FALSE) orientation;
	 * @return an estimate of the long period induced by a two-node overlap loop.
	 */
	private static final int fixPeriodicLoops_getPeriod(Edge edge, Node node1, Node node2, boolean forwardOrReverse) {
		final int i, j;
		
		if (forwardOrReverse) {
			i=Constants.overlap2id(Constants.OVERLAP_SUFFIX_PREFIX);
			j=Constants.overlap2id(Constants.OVERLAP_PREFIX_SUFFIX);
		}
		else {
			i=Constants.overlap2id(Constants.OVERLAP_SUFFIX_SUFFIX);
			j=Constants.overlap2id(Constants.OVERLAP_PREFIX_PREFIX);
		}
		return ( node1.length()-edge.overhangs[i][0] +
				 node1.length()-edge.overhangs[j][0] +
				 node2.length()-edge.overhangs[i][1] +
				 node2.length()-edge.overhangs[j][1] )>>2;
	}
	
	
	/**
	 * Uses an edge that contains both a suffix-prefix alignment and a shared substring
	 * alignment, to infer that one of its adjacent nodes is periodic. Such node is added
	 * to $loops$.
	 *
	 * Remark: the procedure assumes that none of the nodes adjacent to the edge is 
	 * periodic.
	 *
	 * @param $node1$ is the node with a suffix overlap, $node2$ with a prefix overlap;
	 * @return the last element of $loops$ after the procedure ends.
	 */
	private static final int fixOverlapSharedSubstringEdges_suffixPrefix(int node1, int length1, int sharedSubstringOverhangLeft1, int sharedSubstringOverhangRight1, int overlapOverhang1, int node2, int length2, int sharedSubstringOverhangLeft2, int sharedSubstringOverhangRight2, int overlapOverhang2, int[] loops, int lastLoop, int threshold, int[] stats) {
		int p, q, pPrime, qPrime, x, y, xPrime, yPrime, xDoublePrime, yDoublePrime, period;
		double overlapRatio, sharedSubstringRatio;
		
		sharedSubstringRatio=((double)(length1-sharedSubstringOverhangLeft1-sharedSubstringOverhangRight1+1))/(length2-sharedSubstringOverhangLeft2-sharedSubstringOverhangRight2+1);
		p=overlapOverhang1;
	   	q=length2-overlapOverhang2;
		overlapRatio=(length1-p)/(q+1.0);
		x=sharedSubstringOverhangLeft2;
		xDoublePrime=sharedSubstringOverhangLeft1;
		yDoublePrime=length2-sharedSubstringOverhangRight2;
		if (xDoublePrime<=threshold && x<q-threshold && yDoublePrime>=q) {
			qPrime=(int)(xDoublePrime+(q-x)*sharedSubstringRatio);
			xPrime=(int)(p+x*overlapRatio);
			if (qPrime>=xPrime-threshold) {
				period=length1-Math.max(qPrime,xPrime);
				if (nodesArray[node1].period==-1) {
					nodesArray[node1].hasLongPeriod=true;
					nodesArray[node1].period=period;
				}
				else if (nodesArray[node1].hasLongPeriod) nodesArray[node1].period=Math.min(nodesArray[node1].period,period);
				loops[++lastLoop]=node1;
				stats[nodesArray[node1].type]++;
			}
		}
		y=length1-sharedSubstringOverhangRight1;
		if (yDoublePrime>=length2-threshold && y>p+threshold && xDoublePrime<=p) {
			sharedSubstringRatio=1.0/sharedSubstringRatio;
			pPrime=(int)(x+(p-xDoublePrime)*sharedSubstringRatio);
			overlapRatio=1.0/overlapRatio;
			yPrime=(int)((y-p)*overlapRatio);
			if (pPrime<=yPrime+threshold) {
				period=Math.min(pPrime,yPrime);
				if (nodesArray[node2].period==-1) {
					nodesArray[node2].hasLongPeriod=true;
					nodesArray[node2].period=period;
				}
				else if (nodesArray[node2].hasLongPeriod) nodesArray[node2].period=Math.min(nodesArray[node2].period,period);
				loops[++lastLoop]=node2;
				stats[nodesArray[node2].type]++;
			}
		}
		return lastLoop;
	}
	
	
	/**
	 * Remark: $node1$ and $node2$ have no particular meaning and could be interchanged.
	 */
	private static final int fixOverlapSharedSubstringEdges_suffixSuffix(int node1, int length1, int sharedSubstringOverhangLeft1, int sharedSubstringOverhangRight1, int overlapOverhang1, int node2, int length2, int sharedSubstringOverhangLeft2, int sharedSubstringOverhangRight2, int overlapOverhang2, int[] loops, int lastLoop, int threshold, int[] stats) {
		int p, q, pPrime, qPrime, x, y, xPrime, yPrime, xDoublePrime, yDoublePrime, period;
		double overlapRatio, sharedSubstringRatio;
		
		sharedSubstringRatio=((double)(length1-sharedSubstringOverhangLeft1-sharedSubstringOverhangRight1+1))/(length2-sharedSubstringOverhangLeft2-sharedSubstringOverhangRight2+1);
		p=overlapOverhang1;
	   	q=overlapOverhang2;
		overlapRatio=((double)(length1-p))/(length2-q);
		x=length2-sharedSubstringOverhangRight2;
		xDoublePrime=sharedSubstringOverhangLeft1;
		yDoublePrime=sharedSubstringOverhangLeft2;
		if (xDoublePrime<=threshold && x>q+threshold && yDoublePrime<=q) {
			qPrime=(int)(xDoublePrime+(x-q)*sharedSubstringRatio);
			xPrime=(int)(p+(length2-x)*overlapRatio);
			if (qPrime>=xPrime-threshold) {
				period=length1-Math.max(qPrime,xPrime);
				if (nodesArray[node1].period==-1) {
					nodesArray[node1].hasLongPeriod=true;
					nodesArray[node1].period=period;
				}
				else if (nodesArray[node1].hasLongPeriod) nodesArray[node1].period=Math.min(nodesArray[node1].period,period);
				loops[++lastLoop]=node1;
				stats[nodesArray[node1].type]++;
			}
		}
		y=length1-sharedSubstringOverhangRight1;
		if (yDoublePrime<=threshold && y>p+threshold && xDoublePrime<=p) {
			sharedSubstringRatio=1.0/sharedSubstringRatio;
			pPrime=(int)(x-(p-xDoublePrime)*sharedSubstringRatio);
			overlapRatio=1.0/overlapRatio;
			yPrime=length2-(int)((y-p)*overlapRatio);
			if (pPrime>=yPrime-threshold) {
				period=length2-Math.max(pPrime,yPrime);
				if (nodesArray[node2].period==-1) {
					nodesArray[node2].hasLongPeriod=true;
					nodesArray[node2].period=period;
				}
				else if (nodesArray[node2].hasLongPeriod) nodesArray[node2].period=Math.min(nodesArray[node2].period,period);
				loops[++lastLoop]=node2;
				stats[nodesArray[node2].type]++;
			}
		}
		return lastLoop;
	}
	
	
	/**
	 * Remark: $node1$ and $node2$ have no particular meaning and could be interchanged.
	 */
	private static final int fixOverlapSharedSubstringEdges_prefixPrefix(int node1, int length1, int sharedSubstringOverhangLeft1, int sharedSubstringOverhangRight1, int overlapOverhang1, int node2, int length2, int sharedSubstringOverhangLeft2, int sharedSubstringOverhangRight2, int overlapOverhang2, int[] loops, int lastLoop, int threshold, int[] stats) {
		int p, q, pPrime, qPrime, x, y, xPrime, yPrime, xDoublePrime, yDoublePrime, period;
		double overlapRatio, sharedSubstringRatio;
		
		sharedSubstringRatio=((double)(length1-sharedSubstringOverhangLeft1-sharedSubstringOverhangRight1+1))/(length2-sharedSubstringOverhangLeft2-sharedSubstringOverhangRight2+1);
		p=length1-overlapOverhang1;
	   	q=length2-overlapOverhang2;
		overlapRatio=((double)p)/q;
		x=sharedSubstringOverhangLeft2;
		xDoublePrime=length1-sharedSubstringOverhangRight1;
		yDoublePrime=length2-sharedSubstringOverhangRight2;
		if (xDoublePrime>=length1-threshold && x<q-threshold && yDoublePrime>=q) {
			qPrime=(int)(xDoublePrime-(q-x)*sharedSubstringRatio);
			xPrime=(int)(p-x*overlapRatio);
			if (qPrime<=xPrime+threshold) {
				period=Math.min(qPrime,xPrime);
				if (nodesArray[node1].period==-1) {
					nodesArray[node1].hasLongPeriod=true;
					nodesArray[node1].period=period;
				}
				else if (nodesArray[node1].hasLongPeriod) nodesArray[node1].period=Math.min(nodesArray[node1].period,period);
				loops[++lastLoop]=node1;
				stats[nodesArray[node1].type]++;
			}
		}
		y=sharedSubstringOverhangLeft1;
		if (yDoublePrime>=length2-threshold && y<p-threshold && xDoublePrime>=p) {
			sharedSubstringRatio=1.0/sharedSubstringRatio;
			pPrime=(int)(x+(xDoublePrime-p)*sharedSubstringRatio);
			overlapRatio=1.0/overlapRatio;
			yPrime=(int)((p-y)*overlapRatio);
			if (pPrime<=yPrime+threshold) {
				period=Math.min(pPrime,yPrime);
				if (nodesArray[node2].period==-1) {
					nodesArray[node2].hasLongPeriod=true;
					nodesArray[node2].period=period;
				}
				else if (nodesArray[node2].hasLongPeriod) nodesArray[node2].period=Math.min(nodesArray[node2].period,period);
				loops[++lastLoop]=node2;
				stats[nodesArray[node2].type]++;
			}
		}
		return lastLoop;
	}
	
	
	/**
	 * Considers every edge between a non-periodic node X and any other node, with just 
	 * one overlap type and multiple alignments, as a proof that the overlap substring of
	 * X is periodic. The procedure creates new periodic nodes from the longest such 
	 * overlaps at both the prefix and the suffix of X, and it redirects edges 
	 * accordingly.
	 *
	 * Remark: using multi-edges as a proof of periodicity might be overly aggressive.
	 * This is done also in $fixPeriodicLoops()$, which is assumed to have been 
	 * executed before this procedure. Specifically, the procedure assumes that the union
	 * of the longest multi-edge prefix and of the longest multi-edge suffix overlaps of X
	 * does not cover most of X.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted, and it might create new 
	 * intervals. At the end of the procedure, $nodesArray$ is not necessarily sorted.
	 */
	private static final void fixMultiEdges() {
		boolean hasLongPeriodPrefix, hasLongPeriodSuffix;
		int i, j, k, p, q, x;
		int row, nEdges, fromNode, length, periodPrefix, periodSuffix, lastUndersplitNode, nodeID;
		Node node, otherNode, nodePrefix, nodeSuffix;
		Edge edge;
		Edge tmpEdge = new Edge();
		AlignmentPair tmpPair = new AlignmentPair();
		int[] stats = new int[3];
		int[] undersplitNodeIsUsed = new int[] {Math.POSITIVE_INFINITY,Math.POSITIVE_INFINITY};
		Node[] undersplitNodes = new Node[2];
		
		fromNode=nNodes-1;
		for (i=0; i<nNodes; i++) {
			if (i%100000==0) System.err.println("fixMultiEdges> "+i+" nodes done ("+IO.getPercent(i,nNodes)+"%)");
			node=nodesArray[i];
			if (node.type==Constants.INTERVAL_PERIODIC) continue;
			length=node.length();
			hasLongPeriodPrefix=true; hasLongPeriodSuffix=true;
			periodPrefix=Math.POSITIVE_INFINITY; periodSuffix=Math.POSITIVE_INFINITY;
			p=-1; q=length;  // Relative positions
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (!edge.hasOneType() || !edge.isType_overlap() || edge.getNAlignments()==1 || edge.overhangs==null || !edge.overhangMatricesDiffer()) continue;
				otherNode=nodesArray[edge.getTo(i)];
				row=Constants.overlap2id(edge.overlap);
				if (i==edge.nodeID1) {
					if (edge.overlap==Constants.OVERLAP_SUFFIX_PREFIX || edge.overlap==Constants.OVERLAP_SUFFIX_SUFFIX) {
						x=edge.overhangs[row][0];
						if (x<q) q=x;
						if (otherNode.type==Constants.INTERVAL_PERIODIC && !otherNode.hasLongPeriod) hasLongPeriodSuffix=false;
						if (otherNode.period>0 && otherNode.period<periodSuffix) periodSuffix=otherNode.period;
					}
					else {
						x=length-edge.overhangs[row][0]-1;
						if (x>p) p=x;
						if (otherNode.type==Constants.INTERVAL_PERIODIC && !otherNode.hasLongPeriod) hasLongPeriodPrefix=false;
						if (otherNode.period>0 && otherNode.period<periodPrefix) periodPrefix=otherNode.period;
					}
				}
				else {
					if (edge.overlap==Constants.OVERLAP_PREFIX_SUFFIX || edge.overlap==Constants.OVERLAP_SUFFIX_SUFFIX) {
						x=edge.overhangs[row][1];
						if (x<q) q=x;
						if (otherNode.type==Constants.INTERVAL_PERIODIC && !otherNode.hasLongPeriod) hasLongPeriodSuffix=false;
						if (otherNode.period>0 && otherNode.period<periodSuffix) periodSuffix=otherNode.period;
					}
					else {
						x=length-edge.overhangs[row][1]-1;
						if (x>p) p=x;
						if (otherNode.type==Constants.INTERVAL_PERIODIC && !otherNode.hasLongPeriod) hasLongPeriodPrefix=false;
						if (otherNode.period>0 && otherNode.period<periodPrefix) periodPrefix=otherNode.period;
					}
				}
			}
			lastUndersplitNode=-1;
			if (p>=0) {
				nodeID=getPeriodicInterval(node.start,node.start+p,node.read,i,true,hasLongPeriodSuffix?2:1);
				if (nodeID>=0) nodePrefix=nodesArray[nodeID];
				else {
					nodePrefix = new Node(node.read,Constants.INTERVAL_PERIODIC,-1,node.start,node.start+p,node.isLeftMaximal,true,false,node.isContainedInterval,-1,-1,periodSuffix==Math.POSITIVE_INFINITY?0:periodSuffix,hasLongPeriodSuffix,-1,-1,0,0);
					fromNode++; expandTo(fromNode);
					nodePrefix.nodeID=fromNode; nodesArray[fromNode]=nodePrefix; nNeighbors[fromNode]=0;
				}
				undersplitNodes[++lastUndersplitNode]=nodePrefix;
			}
			if (q<length) {
				nodeID=getPeriodicInterval(node.start+q,node.end,node.read,i,true,hasLongPeriodPrefix?2:1);
				if (nodeID>=0) nodeSuffix=nodesArray[nodeID];
				else {
					nodeSuffix = new Node(node.read,Constants.INTERVAL_PERIODIC,-1,node.start+q,node.end,true,node.isRightMaximal,false,node.isContainedInterval,-1,-1,periodPrefix==Math.POSITIVE_INFINITY?0:periodPrefix,hasLongPeriodPrefix,-1,-1,0,0);
					fromNode++; expandTo(fromNode);
					nodeSuffix.nodeID=fromNode; nodesArray[fromNode]=nodeSuffix; nNeighbors[fromNode]=0;
				}
				undersplitNodes[++lastUndersplitNode]=nodeSuffix;
			}			
			for (j=0; j<nNeighbors[i]; j++) fixUndersplits_handleEdge(false,node,node.start,j,neighbors[i][j],tmpEdge,undersplitNodes,lastUndersplitNode,undersplitNodeIsUsed,0,stats,tmpPair);
			k=0;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]!=null) {
					if (k!=j) neighbors[i][k]=neighbors[i][j];
					k++;
				}
			}
			nNeighbors[i]=k;
		}
		fromNode=mergeUndersplitNodes(nNodes,fromNode);
		System.out.println("fixMultiEdges> "+stats[2]+" overlap multi-edges redirected, "+(fromNode+1-nNodes)+" periodic intervals created.");
		nNodes=fromNode+1;
		if (IO.CONSISTENCY_CHECKS) {
			nEdges=0;
			for (i=0; i<nNodes; i++) nEdges+=nNeighbors[i];
			if (nEdges%2!=0) {
				System.err.println("fixMultiEdges> ERROR: wrong number of edges at the end of the procedure.");
				System.exit(1);
			}
		}
	}

	
	
	
	
	
	
	
	// ---------------------------------- UNDERSPLITS ------------------------------------
	
	private static final int MIN_INTERVAL_LENGTH_FOR_UNDERSPLIT = Alignments.minAlignmentLength+IO.quantum;  // If too big, some undersplits are not resolved.
	
	/**
	 * A nonperiodic interval $I$ of a read is said to be an \emph{undersplit} if it 
	 * contains an instance of another repeat $R$ as a substring, and if such instance is 
	 * not assigned to an interval. Undersplits can happen because of weaknesses in the 
	 * factorization pipeline. Consider all other instances of $R$ in the read set. If 
	 * they all have been correctly factorized, $I$ will have them as neighbors of 
	 * insertion type in the interval graph, and the alignments with all such neighbors 
	 * will have similar positions in $I$. If some instances of $R$ have been undersplit, 
	 * then $I$ will have them as neighbors of overlap or shared substring type in the 
	 * interval graph, again with similar alignment positions in $I$.
	 *
	 * The procedure tries to detect such cases in the interval graph, and to fix them by 
	 * creating new intervals for the instances of $R$, and by redirecting edges from $I$ 
	 * to such new intervals. If a similar interval already exists in the graph, edges are
	 * redirected to it. The original node $I$ might become disconnected in the graph.
	 * This has similarities to some heuristic described in \cite{bejerano2004into}.
	 *
	 * By reusing existing intervals, the procedure tries also to solve the problem that 
	 * an alignment might be assigned to a container interval rather than to a contained 
	 * interval, thereby inducing wrong edges in the interval graph. This is not likely to
	 * be a problem of $Factorize.alignments2intervals()$, which assigns alignments to
	 * their shortest containing intervals, but of the previous steps, i.e. of implied 
	 * alignments.
	 *
	 * Remark: this procedure detects peaks induced by maximal intervals that align to the
	 * current interval. Such distribution might not contain peaks, even though the same
	 * distribution, but limited to a specific type of alignments (e.g. overlaps), might
	 * contain peaks.
	 *
	 * Remark: the result of this procedure does not depend on the order in which nodes 
	 * with undersplits are processed.
	 *
	 * Remark: this procedure can simplify the overlap graphs of kernels, since it can
	 * transform overlaps into containments or insertions. However, it is not clear 
	 * whether repeat clusters become better separated in the interval graph, since
	 * containment edges might also be created by this procedure.
	 *
	 * Remark: even if the overall structure of the interval graph gets significantly 
	 * simpler after this procedure, the other steps of the pipeline are still useful in
	 * theory, since this procedure acts just in the presence of clear signals.
	 *
	 * Remark: undersplits are more likely to occur in long intervals that are maximal by 
	 * containment, i.e. precisely in kernel nodes.
	 *
	 * Remark: this procedure differs from the criteria used to establish edges in e.g.
	 * $IntervalGraphStep2.inStep2()$. First, because it uses all edges to find peaks.
	 * Second, because the criteria used at construction time do not alter the topology of
	 * the graph, but they just set labels on the edges (e.g. overlap vs shared 
	 * substring).
	 *
	 * Remark: the procedure assumes that $Alignments.minAlignmentLength$ has already been
	 * set, and that the $periodicSubstrings$ list of every node has already been built.
	 * The procedure might create new intervals: in this case, at the end of the procedure
	 * $nodesArray$ is not necessarily sorted.
	 *
	 * @param isSorted TRUE iff $nodesArray$ is sorted by $read,start,nodeID$ when the
	 * procedure is called;
	 * @return TRUE iff at least one edge has been redirected.
	 */
	public static final boolean fixUndersplits(boolean isSorted) {
		final int MIN_EVENTS_PER_PEAK = 10;  // Arbitrary
		int i, j, k, p;
		int last, nUndersplitCandidates, nReusableNodes, fromNode, maxPeaks;
		int previousStats, nUndersplitNodes, nUndersplitSingletons;
		int lastEvent, nPeaks;
		Node query = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		AlignmentPair tmpPair = new AlignmentPair();
		int[] events = new int[maxAlignmentsPerRead<<1];
		int[] sortedEvents = new int[events.length];
		int[] peaks = new int[1000];  // Arbitrary upper bound on $nPeaks+2$
		int[] stats = new int[3];
		int[] undersplitNodeIsUsed, undersplitNodeIsUsedPrime;
		int[][] reusedStats;
		boolean[] peaksLeft, peaksRight, undersplitNodeIsNew;
		Node[] undersplitNodes;

		// Counting the number of undersplit candidates
		DensityEstimationTree.allocateMemory(maxAlignmentsPerRead);
		nUndersplitCandidates=0; maxPeaks=0;	
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].type==Constants.INTERVAL_PERIODIC || nodesArray[i].length()<MIN_INTERVAL_LENGTH_FOR_UNDERSPLIT) continue;
			lastEvent=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				k=neighbors[i][j].isMaximalOverlap(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
				k=neighbors[i][j].isMaximalInsertionOrContainment(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
				k=neighbors[i][j].isMaximalSharedSubstring(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
			}
			if (lastEvent==-1 || ((lastEvent+1)>>1)<MIN_EVENTS_PER_PEAK) continue;
			lastTmpPoint=-1;
			for (j=0; j<lastEvent; j+=2) {
				lastTmpPoint++;
				ensureTmpPoints(lastTmpPoint);
				tmpPoints[lastTmpPoint].position=events[j];
				tmpPoints[lastTmpPoint].mass=1;
			}
			lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);			
			nPeaks=hasPeaks(tmpPoints,lastTmpPoint,peaks,MIN_EVENTS_PER_PEAK);
			if (nPeaks>0) {
				nUndersplitCandidates++;
				if (nPeaks>maxPeaks) maxPeaks=nPeaks;
			}
		}
		System.err.println("Undersplit candidates: "+nUndersplitCandidates+" ("+IO.getPercent(nUndersplitCandidates,nNodes)+"%)");
		maxPeaks<<=1;  // Arbitrary. Reassigning edges might make the new number of peaks arbitrarily large with respect to the old number of peaks.
		expandTo(nNodes+nUndersplitCandidates-1);
		reusedStats = new int[nNodes][2];
			
		// Loading reusable nodes (using the extra space in $nodesArray$ as temporary
		// space)
		p=nNodes-1;
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].type==Constants.INTERVAL_PERIODIC || nodesArray[i].length()<MIN_INTERVAL_LENGTH_FOR_UNDERSPLIT) continue;
			lastEvent=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				k=neighbors[i][j].isMaximalOverlap(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
				k=neighbors[i][j].isMaximalInsertionOrContainment(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
				k=neighbors[i][j].isMaximalSharedSubstring(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
			}
			if (lastEvent==-1 || ((lastEvent+1)>>1)<MIN_EVENTS_PER_PEAK) continue;
			lastTmpPoint=-1;
			for (j=0; j<lastEvent; j+=2) {
				lastTmpPoint++;
				ensureTmpPoints(lastTmpPoint);
				tmpPoints[lastTmpPoint].position=events[j];
				tmpPoints[lastTmpPoint].mass=1;
			}
			lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
			nPeaks=hasPeaks(tmpPoints,lastTmpPoint,peaks,MIN_EVENTS_PER_PEAK);
			if (nPeaks>0) nodesArray[++p]=nodesArray[i];
		}
		Node.order=Node.READ_START_NODEID;
		if (!isSorted && nNodes>1) Arrays.sort(nodesArray,0,nNodes); 
		if (p+1-nNodes>1) Arrays.sort(nodesArray,nNodes,p+1);
		reusableNodes = new Node[nUndersplitCandidates];
		nReusableNodes=loadReusableNodes(nodesArray,nNodes,p,false);
		System.err.println("N. reusable nodes: "+nReusableNodes);
		if (!isSorted) {
			Node.order=Node.NODE_ID;
			if (nNodes>1) Arrays.sort(nodesArray,0,nNodes);
		}		
		
		// Fixing undersplits
		peaksLeft = new boolean[maxPeaks+2]; peaksRight = new boolean[maxPeaks+2];
		undersplitNodes = new Node[((maxPeaks+2)*(maxPeaks+1))>>1];
		undersplitNodeIsNew = new boolean[undersplitNodes.length];
		undersplitNodeIsUsed = new int[undersplitNodes.length];
		undersplitNodeIsUsedPrime = new int[undersplitNodes.length];
		stats[0]=0; stats[1]=0; stats[2]=0;
		nUndersplitNodes=0; nUndersplitSingletons=0;
		fromNode=nNodes;
		for (i=0; i<nNodes; i++) {
			if (i%10000==0) System.err.println("fixUndersplits> processing node "+i+" ("+IO.getPercent(i,nNodes)+"%)");			
			if (nodesArray[i].type==Constants.INTERVAL_PERIODIC || nodesArray[i].length()<MIN_INTERVAL_LENGTH_FOR_UNDERSPLIT) continue;
			lastEvent=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				k=neighbors[i][j].isMaximalOverlap(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
				k=neighbors[i][j].isMaximalInsertionOrContainment(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
				k=neighbors[i][j].isMaximalSharedSubstring(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
			}
			if (lastEvent==-1 || ((lastEvent+1)>>1)<MIN_EVENTS_PER_PEAK) continue;
			lastTmpPoint=-1;
			for (j=0; j<lastEvent; j+=2) {
				lastTmpPoint++;
				ensureTmpPoints(lastTmpPoint);
				tmpPoints[lastTmpPoint].position=events[j];
				tmpPoints[lastTmpPoint].mass=1;
			}
			lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
			nPeaks=hasPeaks(tmpPoints,lastTmpPoint,peaks,MIN_EVENTS_PER_PEAK);			
			if (nPeaks<1) continue;		
			peaks[0]=0; peaks[nPeaks+1]=nodesArray[i].length()-1;
			sortEvents(events,lastEvent,tmpPoints,lastTmpPoint,sortedEvents);
			previousStats=stats[2];			
			fromNode=fixUndersplits_impl(nodesArray[i],peaks,nPeaks+1,sortedEvents,(lastTmpPoint<<1)+1,reusableNodes,nReusableNodes-1,fromNode,query,stats,peaksLeft,peaksRight,undersplitNodes,undersplitNodeIsNew,undersplitNodeIsUsed,undersplitNodeIsUsedPrime,tmpPair,false)+1;
			if (stats[2]!=previousStats) {
				nUndersplitNodes++;
				if (nNeighbors[i]==0) nUndersplitSingletons++;
				reusedStats[i][1]=1;
				for (j=0; j<undersplitNodes.length; j++) {
					if (undersplitNodes[j]==null) break;
					if (!undersplitNodeIsNew[j] && undersplitNodeIsUsed[j]!=0) reusedStats[undersplitNodes[j].nodeID][0]+=undersplitNodeIsUsed[j];
				}
			}
			else {
				//
				// The following has been commented just because of too much output:
				//System.err.println("Node "+i+" has "+nPeaks+" peaks but no edge moved?! Sorted events:");
				//for (int x=0; x<=(lastTmpPoint<<1)+1; x+=2) System.err.print(sortedEvents[x]+",");
				//for (int x=1; x<=(lastTmpPoint<<1)+1; x+=2) System.err.print("("+sortedEvents[x]+"),");
				//System.err.println();
				//
			}			
		}	
		System.err.println("Undersplits fixed. Merging undersplit nodes...");
		fromNode=mergeUndersplitNodes(nNodes,fromNode-1)+1;
		System.err.println("Nodes subjected to undersplit resolution: "+nUndersplitNodes+" ("+IO.getPercent(nUndersplitNodes,nNodes)+"%).");
		System.err.println("Disconnected nodes created by undersplit resolution: "+nUndersplitSingletons+" ("+IO.getPercent(nUndersplitSingletons,nUndersplitNodes)+"%).");
		System.err.println("Nodes created by undersplit resolution: "+(fromNode-nNodes)+" ("+IO.getPercent(fromNode-nNodes,nNodes)+"%).");
		System.err.println("Edges redirected by undersplit resolution: "+stats[2]);
		//System.err.println("Old nodes actually reused as undersplit nodes:");
		j=0;
		for (i=0; i<nNodes; i++) {
			if (reusedStats[i][0]!=0) j++;
		}
		System.err.println(j+" total nodes reused");
		DensityEstimationTree.deallocateMemory();
		nNodes=fromNode;
		return nUndersplitNodes>0;
	}
	
	
	/**
	 * Checks the consistency of $neighbors[from..to]$.
	 *
	 * @param events temporary space.
	 */
	public static final void checkConsistency(int from, int to, int[] events) {
		boolean found, foundPrime;
		int i, j, k, p;
		int lastEvent, otherNode, lowerBound;
		int[] tmpArray = new int[10];

		for (i=from; i<=to; i++) {
			if (nodesArray[i]==null) {
				System.err.println("checkConsistency> ERROR: null "+i);
				System.exit(1);
			}
			if (nodesArray[i].nodeID!=i) {
				System.err.println("checkConsistency> ERROR: node "+i+" has id="+nodesArray[i].nodeID+"?!");
				System.err.println("node at position "+i+" in nodesArray: "+nodesArray[i]);
				System.err.println("node at position "+nodesArray[i].nodeID+" in nodesArray: "+nodesArray[nodesArray[i].nodeID]);
				System.exit(1);
			}
			if (nodesArray[i].length()<10) {
				System.err.println("checkConsistency> ERROR: node "+i+" has length="+nodesArray[i].length()+"?!");
				System.err.println(nodesArray[i]);
				System.exit(1);
			}
			if (nodesArray[i].type==Constants.INTERVAL_PERIODIC || nodesArray[i].length()<MIN_INTERVAL_LENGTH_FOR_UNDERSPLIT) continue;
			lastEvent=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				k=neighbors[i][j].isMaximalOverlap(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
				k=neighbors[i][j].isMaximalInsertionOrContainment(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
				k=neighbors[i][j].isMaximalSharedSubstring(i,events,lastEvent+1);
				if (k>lastEvent) lastEvent=k;
			}
			for (j=0; j<=lastEvent; j++) {
				if (events[j]<0 || events[j]>nodesArray[i].length()) {
					System.err.println("checkConsistency> ERROR: WRONG EVENT: "+events[j]+" length="+nodesArray[i].length()+" node="+i+" nodeID="+nodesArray[i].nodeID);
					System.err.println("checkConsistency> node: "+nodesArray[i]);
					System.err.print("checkConsistency> events: ");
					for (int x=0; x<=lastEvent; x++) System.err.print(events[x]+" ");
					System.err.println();
					System.err.println("checkConsistency> edges to neighbors of the node:");
					for (int x=0; x<=nNeighbors[i]; x++) System.err.println(neighbors[i][x]);
					System.exit(1);
				}
			}
			
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				if (neighbors[i][j].nDistinctTypes()==0) {
					System.err.println("checkConsistency> ERROR: the following edge has no valid step2 type?!");
					System.err.println(neighbors[i][j]);
					System.err.println("node 1: "+nodesArray[neighbors[i][j].nodeID1]);
					System.err.println("node 2: "+nodesArray[neighbors[i][j].nodeID2]);
					System.exit(1);
				}
				if (neighbors[i][j].nTypes()>neighbors[i][j].getNAlignments()) {
					System.err.println("checkConsistency> ERROR: the following edge has more alignment types than alignments?!");
					System.err.println(neighbors[i][j]);
					System.err.println("node 1: "+nodesArray[neighbors[i][j].nodeID1]);
					System.err.println("node 2: "+nodesArray[neighbors[i][j].nodeID2]);
					System.exit(1);
				}
				/*if (neighbors[i][j].avgDiffs>neighbors[i][j].getNAlignments()) {
					System.err.println("checkConsistency> ERROR: the following edge has too large avgDiffs:");
					System.err.println(neighbors[i][j]);
					System.exit(1);
				}*/
			}
			
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				if (neighbors[i][j].overhangs!=null) {
					for (int z=0; z<4; z++) {
						if ((neighbors[i][j].overhangs[z][0]!=-1 || neighbors[i][j].overhangs[z][1]!=-1) && neighbors[i][j].overlap!=-1 && (neighbors[i][j].overlap&Constants.id2overlap(z))==0) {
							System.err.println("checkConsistency> ERROR: the following edge has overhangs for an absent overlap type: "+z+" :: "+neighbors[i][j].overhangs[z][0]+","+neighbors[i][j].overhangs[z][1]);
							System.err.println(neighbors[i][j]);
							System.err.println("node 1: "+nodesArray[neighbors[i][j].nodeID1]);
							System.err.println("node 2: "+nodesArray[neighbors[i][j].nodeID2]);
							System.exit(1);
						}
						if ((neighbors[i][j].overhangs[z][0]==-1 || neighbors[i][j].overhangs[z][1]==-1) && neighbors[i][j].overlap!=-1 && (neighbors[i][j].overlap&Constants.id2overlap(z))!=0) {
							System.err.println("checkConsistency> ERROR: the following edge does not have overhangs for an overlap type: "+z+" :: "+neighbors[i][j].overhangs[z][0]+","+neighbors[i][j].overhangs[z][1]);
							System.err.println(neighbors[i][j]);
							System.err.println("node 1: "+nodesArray[neighbors[i][j].nodeID1]);
							System.err.println("node 2: "+nodesArray[neighbors[i][j].nodeID2]);
							System.exit(1);
						}
					}
				}
			}
			
			for (j=0; j<nNeighbors[i]; j++) {
				if (nodesArray[neighbors[i][j].nodeID1].read>nodesArray[neighbors[i][j].nodeID2].read) {
					System.err.println("checkConsistency> ERROR: the nodes of the following edge do not satisfy nodeID1.read<=nodeID2.read:");
					System.err.println(neighbors[i][j]);
					System.err.println("node "+neighbors[i][j].nodeID1+": "+nodesArray[neighbors[i][j].nodeID1]);
					System.err.println("node "+neighbors[i][j].nodeID2+": "+nodesArray[neighbors[i][j].nodeID2]);
					System.exit(1);
				}
				/*if ( neighbors[i][j].sharedSubstring!=-1 && neighbors[i][j].sharedSubstringOverhangs!=null && 
					 ( (neighbors[i][j].sharedSubstringOverhangs[0]==0 && neighbors[i][j].sharedSubstringOverhangs[1]==0) ||
					   (neighbors[i][j].sharedSubstringOverhangs[2]==0 && neighbors[i][j].sharedSubstringOverhangs[3]==0)
				     ) 
				   ) {
					System.err.println("checkConsistency> ERROR: node "+i+" has an edge with zero shared substring overhangs:");
					System.err.println(neighbors[i][j].sharedSubstringOverhangs[0]+","+neighbors[i][j].sharedSubstringOverhangs[1]+","+neighbors[i][j].sharedSubstringOverhangs[2]+","+neighbors[i][j].sharedSubstringOverhangs[3]);
					System.err.println("node1: "+nodesArray[neighbors[i][j].nodeID1]);
					System.err.println("node2: "+nodesArray[neighbors[i][j].nodeID2]);
					System.exit(1);
				}*/
			}
			for (j=0; j<nNeighbors[i]; j++) {
				if (nodesArray[i].type==Constants.INTERVAL_PERIODIC || nodesArray[neighbors[i][j].getTo(i)].type==Constants.INTERVAL_PERIODIC) continue;
				/*if (neighbors[i][j].containment!=-1 && neighbors[i][j].insertion!=-1) {
					System.err.println("checkConsistency> ERROR: node "+i+" has an edge with both insertion and containment:");
					System.err.println(neighbors[i][j]);
					System.err.println("node "+neighbors[i][j].nodeID1+": "+nodesArray[neighbors[i][j].nodeID1]);
					System.err.println("node "+neighbors[i][j].nodeID2+": "+nodesArray[neighbors[i][j].nodeID2]);
					System.exit(1);
				}*/
				if (neighbors[i][j].containment!=-1) {
					p=neighbors[i][j].projectOnto_containment(i,tmpArray,-1);
					if (p>=0 && (tmpArray[0]<0 || tmpArray[1]<0 || tmpArray[1]>nodesArray[i].length()+IDENTITY_THRESHOLD)) {
						System.err.println("checkConsistency> ERROR: node "+i+" has an edge with wrong containment overhangs:");
						System.err.println(neighbors[i][j].containmentOverhangLeft+", "+neighbors[i][j].containmentOverhangRight);
						System.err.println(neighbors[i][j]);
						System.exit(1);
					}
					p=neighbors[i][j].projectOnto_containment(neighbors[i][j].getTo(i),tmpArray,-1);
					if (p>=0 && (tmpArray[0]<0 || tmpArray[1]<0 || tmpArray[1]>nodesArray[neighbors[i][j].getTo(i)].length()+IDENTITY_THRESHOLD)) {
						System.err.println("checkConsistency> ERROR: node "+i+" has an edge with wrong containment overhangs:");
						System.err.println(neighbors[i][j].containmentOverhangLeft+", "+neighbors[i][j].containmentOverhangRight);
						System.err.println(neighbors[i][j]);
						System.exit(1);
					}
				}
				if (neighbors[i][j].insertion!=-1) {
					p=neighbors[i][j].projectOnto_insertion(i,tmpArray,-1);
					if (p>=0 && (tmpArray[0]<0 || tmpArray[1]<0 || tmpArray[1]>nodesArray[i].length()+IDENTITY_THRESHOLD)) {
						System.err.println("checkConsistency> ERROR: node "+i+" has an edge with wrong insertion overhangs:");
						System.err.println(neighbors[i][j].containmentOverhangLeft+", "+neighbors[i][j].containmentOverhangRight);
						System.err.println(neighbors[i][j]);
						System.exit(1);
					}
					p=neighbors[i][j].projectOnto_insertion(neighbors[i][j].getTo(i),tmpArray,-1);
					if (p>=0 && (tmpArray[0]<0 || tmpArray[1]<0 || tmpArray[1]>nodesArray[neighbors[i][j].getTo(i)].length()+IDENTITY_THRESHOLD)) {
						System.err.println("checkConsistency> ERROR: node "+i+" has an edge with wrong insertion overhangs:");
						System.err.println(neighbors[i][j].containmentOverhangLeft+", "+neighbors[i][j].containmentOverhangRight);
						System.err.println(neighbors[i][j]);
						System.exit(1);
					}
				}
				if (neighbors[i][j].overlap!=-1 && neighbors[i][j].overhangs!=null) {
					for (k=Constants.OVERLAP_PREFIX_PREFIX; k<=Constants.OVERLAP_SUFFIX_SUFFIX; k<<=1) {
						p=neighbors[i][j].projectOnto_overlap(i,tmpArray,k,-1);
						if (p>=0 && (tmpArray[0]<0 || tmpArray[1]<0 || tmpArray[1]>nodesArray[i].length()+IDENTITY_THRESHOLD)) {
							System.err.println("checkConsistency> ERROR: node "+i+" has an edge with wrong overlap overhangs: (overlap="+k+")");
							for (int x=0; x<neighbors[i][j].overhangs.length; x++) {
								for (int y=0; y<neighbors[i][j].overhangs[x].length; y++) System.err.println(neighbors[i][j].overhangs[x][y]+",");
								System.err.println();
							}
							System.err.println(neighbors[i][j]);
							System.exit(1);
						}
						p=neighbors[i][j].projectOnto_overlap(neighbors[i][j].getTo(i),tmpArray,k,-1);
						if (p>=0 && (tmpArray[0]<0 || tmpArray[1]<0 || tmpArray[1]>nodesArray[neighbors[i][j].getTo(i)].length()+IDENTITY_THRESHOLD)) {
							System.err.println("checkConsistency> ERROR: node "+i+" has an edge with wrong overlap overhangs: (overlap="+k+")");
							for (int x=0; x<neighbors[i][j].overhangs.length; x++) {
								for (int y=0; y<neighbors[i][j].overhangs[x].length; y++) System.err.println(neighbors[i][j].overhangs[x][y]+",");
								System.err.println();
							}
							System.err.println(neighbors[i][j]);
							System.exit(1);
						}
					}
				}
				if (neighbors[i][j].sharedSubstring!=-1 && neighbors[i][j].sharedSubstringOverhangs!=null) {
					p=neighbors[i][j].projectOnto_sharedSubstring(i,tmpArray,-1);
					if (p>=0 && (tmpArray[0]<0 || tmpArray[1]<0 || tmpArray[1]>nodesArray[i].length()+IDENTITY_THRESHOLD)) {
						System.err.println("checkConsistency> ERROR: node "+i+" has an edge with wrong shared substring overhangs:");
						for (int x=0; x<neighbors[i][j].sharedSubstringOverhangs.length; x++) System.err.print(neighbors[i][j].sharedSubstringOverhangs[x]+",");
						System.err.println();
						System.err.println(neighbors[i][j]);
						System.exit(1);
					}
					p=neighbors[i][j].projectOnto_sharedSubstring(neighbors[i][j].getTo(i),tmpArray,-1);
					if (p>=0 && (tmpArray[0]<0 || tmpArray[1]<0 || tmpArray[1]>nodesArray[neighbors[i][j].getTo(i)].length()+IDENTITY_THRESHOLD)) {
						System.err.println("checkConsistency> ERROR: node "+i+" has an edge with wrong shared substring overhangs:");
						for (int x=0; x<neighbors[i][j].sharedSubstringOverhangs.length; x++) System.err.print(neighbors[i][j].sharedSubstringOverhangs[x]+",");
						System.err.println();
						System.err.println(neighbors[i][j]);
						System.exit(1);
					}
				}
			}
			for (j=0; j<nNeighbors[i]; j++) {
				if ( ((nodesArray[i].type==Constants.INTERVAL_PERIODIC && !nodesArray[i].hasLongPeriod) && (nodesArray[neighbors[i][j].getTo(i)].type!=Constants.INTERVAL_PERIODIC || nodesArray[neighbors[i][j].getTo(i)].hasLongPeriod) && neighbors[i][j].sharedSubstring==-1 && neighbors[i][j].insertion==-1 ) ||
					 ((nodesArray[i].type!=Constants.INTERVAL_PERIODIC || nodesArray[i].hasLongPeriod) && (nodesArray[neighbors[i][j].getTo(i)].type==Constants.INTERVAL_PERIODIC && !nodesArray[neighbors[i][j].getTo(i)].hasLongPeriod) && neighbors[i][j].sharedSubstring==-1 && neighbors[i][j].insertion==-1 )
				   ) {
				   System.err.println("checkConsistency> ERROR: nodes "+i+" and "+neighbors[i][j].getTo(i)+" are shortPeriodPeriodic-nonPeriodicOrLongPeriod with a wrong type:");
				   System.err.println("edge: "+neighbors[i][j]);
				   System.err.println("node1: "+nodesArray[i]);
				   System.err.println("node2: "+nodesArray[neighbors[i][j].getTo(i)]);
				   System.exit(1);
				}
			}
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].flexibleOverlap!=0 && neighbors[i][j].flexibleOverlap!=-1) {
					if (neighbors[i][j].maxOverhangs==null) {
						System.err.println("checkConsistency> ERROR: the following edge has flexibleOverlap="+neighbors[i][j].flexibleOverlap+" but maxOverhangs=null:");
						System.err.println(neighbors[i][j]);
						System.exit(1);
					}
					p=-1;
					for (k=Constants.OVERLAP_PREFIX_PREFIX; k<=Constants.OVERLAP_SUFFIX_SUFFIX; k<<=1) {
						p++;
						if ((neighbors[i][j].flexibleOverlap&k)==0) continue;
						if (neighbors[i][j].maxOverhangs[p][0]<neighbors[i][j].overhangs[p][0] || neighbors[i][j].maxOverhangs[p][1]<neighbors[i][j].overhangs[p][1]) {
							System.err.println("checkConsistency> ERROR: the following edge has maxOverhangs<overhangs:");
							System.err.println(neighbors[i][j]);
							System.exit(1);
						}
					}
				}
			}
		}
		for (i=from; i<=to; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].nodeID1!=i && neighbors[i][j].nodeID2!=i) {
					System.err.println("checkConsistency> ERROR: node "+i+" has an edge not related to it: "+neighbors[i][j]);
					System.exit(1);
				}
			}
		}
		for (i=from; i<=to; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				otherNode=neighbors[i][j].getTo(i);
				for (k=j+1; k<nNeighbors[i]; k++) {
					if (neighbors[i][k].getTo(i)==otherNode && neighbors[i][k].supplement==neighbors[i][j].supplement) {
						System.err.println("checkConsistency> ERROR: node "+i+" has multiple edges to node "+otherNode+", and they are "+(neighbors[i][j]==neighbors[i][k]?"the same":"different")+" objects");
						System.err.println(neighbors[i][j]);
						System.err.println(neighbors[i][k]);
						System.err.println("node "+i+": "+nodesArray[i]);
						System.err.println("node "+otherNode+": "+nodesArray[otherNode]);
						System.exit(1);
					}
				}
				foundPrime=false; found=false;
				for (k=0; k<nNeighbors[otherNode]; k++) {
					if (neighbors[otherNode][k].getTo(otherNode)==i) {
						foundPrime=true;
						if (neighbors[otherNode][k].supplement==neighbors[i][j].supplement) {
							found=true;
							break;
						}
					}
				}
				if (!found) {
					System.err.println("checkConsistency> ERROR: node "+i+" has an edge to node "+otherNode+", but not vice versa (of the same type).  foundPrime="+foundPrime);
					System.err.println("edge: "+neighbors[i][j]);
					System.err.println("node: "+nodesArray[i]);
					System.err.println("node: "+nodesArray[otherNode]);
					System.exit(1);
				}
			}
		}
		for (i=to-1; i>=from; i--) {
			for (j=i-1; j>=0; j--) {
				if (nodesArray[i]==nodesArray[j]) {
					System.err.println("checkConsistency> ERROR: nodes "+i+" and "+j+" point the same object with nodeID="+nodesArray[i].nodeID);
					System.exit(1);
				}
			}
		}
		for (i=from; i<=to; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				lowerBound=neighbors[i][j].getNOverlaps()+(neighbors[i][j].containment==-1?0:1)+(neighbors[i][j].insertion==-1?0:1)+(neighbors[i][j].sharedSubstring==-1?0:1);
				if (neighbors[i][j].getNAlignments()<lowerBound) {
					System.err.println("checkConsistency> ERROR: the following edge has "+neighbors[i][j].getNAlignments()+" alignments but "+lowerBound+" types: "+neighbors[i][j]);
					System.err.println(nodesArray[i]);
					System.err.println(nodesArray[neighbors[i][j].getTo(i)]);
					System.exit(1);
				}
			}
		}
		for (i=from; i<=to; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				if (neighbors[i][j].containment==Constants.CONTAINMENT_IDENTICAL && neighbors[i][j].insertion!=-1) {
					System.err.println("checkConsistency> ERROR: the following edge is of type CONTAINMENT_IDENTICAL and INSERTION?! "+neighbors[i][j]);
					System.err.println(nodesArray[neighbors[i][j].nodeID1]);
					System.err.println(nodesArray[neighbors[i][j].nodeID2]);
					System.err.println("overhangs: "+neighbors[i][j].containmentOverhangLeft+","+neighbors[i][j].containmentOverhangRight);
					System.exit(1);
				}
			}
		}
	}
	
	
	/**
	 * After having transformed an edge $(U,V)$ from the current node $U$, into an edge 
	 * $(W,V)$ from the corresponding undersplit node $W$, it might happen that $(W,V)$ 
	 * occurs twice in $neighbors[V]$.
	 *
	 * Remark: the procedure assumes that at most one edge of $neighbors[V]$ is equivalent
	 * to $edge$.
	 *
	 * @param neighbor ID of node $V$;
	 * @return the position of $edge$ in $neighbors[neighbor]$ at the end of the 
	 * procedure, or -1 if $edge$ is not found.
	 */
	private static final int mergeIdenticalEdges(int neighbor, Edge edge) {	
		int i, j;
		int edgePosition, equivalentEdgePosition;	

		if (nNeighbors[neighbor]<1) return -1;
		if (nNeighbors[neighbor]==1) return neighbors[neighbor][0]==edge?0:-1;
		edgePosition=-1; equivalentEdgePosition=-1;
		for (i=0; i<nNeighbors[neighbor]; i++) {
			if (neighbors[neighbor][i]==edge) {
				edgePosition=i;
				continue;
			}
			if (neighbors[neighbor][i].equals(edge)) {				
				equivalentEdgePosition=i;
				edge.addTypes(neighbors[neighbor][i]);
				edge.addDiffs(neighbors[neighbor][i]);
				edge.addOrientation(neighbors[neighbor][i]);
			}
		}
		if (equivalentEdgePosition==-1) return edgePosition;
		j=equivalentEdgePosition;
		for (i=equivalentEdgePosition+1; i<nNeighbors[neighbor]; i++) neighbors[neighbor][j++]=neighbors[neighbor][i];
		nNeighbors[neighbor]=j;
		return edgePosition<equivalentEdgePosition?edgePosition:edgePosition-1;
	}
	
	
	/**
	 * Stores in $destination[0..(lastPoint+1)*2-1]$ the sequence of sorted and compacted 
	 * events in $points[0..lastPoint]$, in which each value is followed by a side flag:
	 * 0=left; 1=right; 2=left and right.
	 *
	 * @param points sorted and compacted.
	 */
	private static final void sortEvents(int[] events, int lastEvent, Point[] points, int lastPoint, int[] destination) {
		int i, j, k;
		
		k=-1;
		for (i=0; i<=lastPoint; i++) {
			destination[++k]=(int)points[i].position;
			destination[k+1]=-1;
			for (j=0; j<=lastEvent; j+=2) {
				if (events[j]!=destination[k]) continue;
				if (events[j+1]==0) {
					if (destination[k+1]==-1) destination[k+1]=0;
					else if (destination[k+1]==1) destination[k+1]=2;
				}
				else {
					if (destination[k+1]==-1) destination[k+1]=1;
					else if (destination[k+1]==0) destination[k+1]=2;
				}
				k++;
			}
		}
	}
	
	
	/**
	 * Fits a density estimation tree on $points[0..lastPoint]$, and finds runs of 
	 * local-maximum leaves such that each has mass at least $minEventsPerPeak$. 
	 * The procedure stores the centers of mass of such runs (if any) in $peaks[1..X]$, 
	 * sorted, where $X$ is returned in output.
	 *
	 * Remark: $minEventsPerPeak$ is typically a very weak lower bound. Considering only 
	 * peaks with very large mass would be suboptimal, since: (1) it might disregard 
	 * smaller peaks that are nonetheless clear; (2) the number of events of a peak is 
	 * related to the frequency of the corresponding repeat, and there is no clean way to 
	 * compare such frequencies.
	 *
	 * Remark: the procedure assumes that $Alignments.minAlignmentLength$ has already been
	 * set, that $lastPoint+1 \geq minEventsPerPeak$, and that $points$ is sorted and 
	 * compacted.
	 *
	 * @return the number of peaks found (zero if no peak was found).
	 */
	private static final int hasPeaks(Point[] points, int lastPoint, int[] peaks, int minEventsPerPeak) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_INTERVAL_LENGTH = IO.quantum>>1;
		final int MIN_LOCAL_MAX_DISTANCE = IO.quantum>>1;
		int i, k;
		int nLocalMaximumLeaves;
		int[] tmpPeaks;
		
		if (points[lastPoint].position-points[0].position<=IDENTITY_THRESHOLD) {
			peaks[1]=Math.round(Points.getCenterOfMass(points,0,lastPoint,false,-1,lastPoint));
			return 1;
		}
		if (Points.areUniformlyDistributed(points,0,lastPoint,true,(points[lastPoint].position-points[0].position)/Points.DEFAULT_NBINS)) return 0;
		nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(points,0,lastPoint,MIN_INTERVAL_LENGTH,MIN_LOCAL_MAX_DISTANCE,true,-1,-1,-1,true,true);
		if (nLocalMaximumLeaves==0) return 0;
		DensityEstimationTree.markRunsOfLocalMaximumLeaves(points);
		k=0;
		for (i=0; i<DensityEstimationTree.nRuns; i++) {
			if (DensityEstimationTree.getMassOfRun(i,points)>=minEventsPerPeak) k++;
		}
		if (peaks.length<k) peaks = new int[1+k];
		k=0;
		for (i=0; i<DensityEstimationTree.nRuns; i++) {
			if (DensityEstimationTree.getMassOfRun(i,points)>=minEventsPerPeak) peaks[++k]=Math.round(DensityEstimationTree.getCenterOfMassOfRun(i,points,lastPoint,IDENTITY_THRESHOLD));
		}
		return k;
	}
	
	
	/**
	 * Loads in $reusableNodes[0..X-1]$, where $X$ is returned in output, the intervals of
	 * $newNodesArray[0..first-1]$ that are approximately contained in an interval of
	 * $newNodesArray[first..last]$ (if $prefixOrSuffix=FALSE$), or that are approximately
	 * identical to a suffix of a suffix interval in $newNodesArray[first..last]$ that
	 * starts at its $firstMaximalStart$ coordinate (if $prefixOrSuffix=TRUE$). The 
	 * symmetrical holds for prefix intervals. The result is sorted by $read$, then by 
	 * $start$.
	 *
	 * Remark: reusable nodes might be periodic.
	 *
	 * Remark: the procedures in the undersplit pipeline assume that reusable nodes are
	 * not themselves undersplit.
	 *
	 * @param newNodesArray intervals $[0..first-1]$ and $[first..last]$ are assumed to be
	 * each independently sorted by $read$, then by $start$.
	 */
	private static final int loadReusableNodes(Node[] newNodesArray, int first, int last, boolean prefixOrSuffix) {
		final int RESIZE_UNIT = 100;
		int i, j, k;
		int firstIForNextJ;
		Node[] tmpArray;

		for (i=0; i<first; i++) newNodesArray[i].visited=-1;  // Using the $visited$ field to mark added nodes.
		i=0; j=first; k=-1; firstIForNextJ=-1;
		while (i<first && j<=last) {
			if (newNodesArray[i].visited!=-1 || newNodesArray[i].nodeID==newNodesArray[j].nodeID || newNodesArray[i].read<newNodesArray[j].read) {
				i++;
				continue;
			}
			if (newNodesArray[i].read>newNodesArray[j].read) {
				j++;
				if (firstIForNextJ!=-1) {
					i=firstIForNextJ;
					firstIForNextJ=-1;
				}
				continue;
			}
			if (firstIForNextJ==-1 && j<last && newNodesArray[j+1].read==newNodesArray[j].read && newNodesArray[i].end>=newNodesArray[j+1].start) firstIForNextJ=i;
			if (newNodesArray[i].start<newNodesArray[j].start-Intervals.absoluteThreshold) {
				i++;
				continue;
			}
			if (newNodesArray[i].start>=newNodesArray[j].end) {
				j++;
				if (firstIForNextJ!=-1) {
					i=firstIForNextJ;
					firstIForNextJ=-1;
				}
				continue;
			}
			if (Intervals.areApproximatelyIdentical(newNodesArray[i].start,newNodesArray[i].end,newNodesArray[j].start,newNodesArray[j].end)) {
				i++;
				continue;
			}
			if (prefixOrSuffix) {
				if (newNodesArray[j].type==Constants.INTERVAL_DENSE_SUFFIX) {
					if (Intervals.areApproximatelyIdentical(newNodesArray[i].start,newNodesArray[i].end,newNodesArray[j].firstMaximalStart,newNodesArray[j].end)) {
						k++;
						if (k==reusableNodes.length) {
							tmpArray = new Node[reusableNodes.length+RESIZE_UNIT];
							System.arraycopy(reusableNodes,0,tmpArray,0,reusableNodes.length);
							reusableNodes=tmpArray;
						}
						reusableNodes[k]=newNodesArray[i];
						newNodesArray[i].visited=0;
					}
				}
				else if (newNodesArray[j].type==Constants.INTERVAL_DENSE_PREFIX) {
					if (Intervals.areApproximatelyIdentical(newNodesArray[i].start,newNodesArray[i].end,newNodesArray[j].start,newNodesArray[j].lastMaximalEnd)) {
						k++;
						if (k==reusableNodes.length) {
							tmpArray = new Node[reusableNodes.length+RESIZE_UNIT];
							System.arraycopy(reusableNodes,0,tmpArray,0,reusableNodes.length);
							reusableNodes=tmpArray;
						}
						reusableNodes[k]=newNodesArray[i];
						newNodesArray[i].visited=0;
					}
				}
				else if (newNodesArray[j].type==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
					if ( Intervals.areApproximatelyIdentical(newNodesArray[i].start,newNodesArray[i].end,newNodesArray[j].firstMaximalStart,newNodesArray[j].end) ||
						 Intervals.areApproximatelyIdentical(newNodesArray[i].start,newNodesArray[i].end,newNodesArray[j].start,newNodesArray[j].lastMaximalEnd)
					   ) {
						k++;
						if (k==reusableNodes.length) {
							tmpArray = new Node[reusableNodes.length+RESIZE_UNIT];
							System.arraycopy(reusableNodes,0,tmpArray,0,reusableNodes.length);
							reusableNodes=tmpArray;
						}
   						reusableNodes[k]=newNodesArray[i];
   						newNodesArray[i].visited=0;	   
					}
				}
			}
			else {
				if (Intervals.isApproximatelyContained(newNodesArray[i].start,newNodesArray[i].end,newNodesArray[j].start,newNodesArray[j].end)) {
					k++;
					if (k==reusableNodes.length) {
						tmpArray = new Node[reusableNodes.length+RESIZE_UNIT];
						System.arraycopy(reusableNodes,0,tmpArray,0,reusableNodes.length);
						reusableNodes=tmpArray;
					}
					reusableNodes[k]=newNodesArray[i];
					newNodesArray[i].visited=0;
				}
			}
			i++;
		}
		return k+1;
	}
	
	
	/**
	 * @return a node in $reusableNodes[0..lastReusableNode]$ that is an approximately 
	 * identical interval to $query$; NULL if no such node can be found;
	 * @param reusableNodes sorted by $read$, then by $start$.
	 */
	private static final Node getReusableNode(Node query, Node[] reusableNodes, int lastReusableNode) {
		int i, j;
		int previousOrder;
		
		previousOrder=Node.order;
		Node.order=Node.READ_START_NODEID;
		i=Arrays.binarySearch(reusableNodes,0,lastReusableNode+1,query);
		if (i<0) {
			i=-i-1;
			i=Math.min(i,lastReusableNode);
		}
		j=i;
		while (j>=0) {
			if (reusableNodes[j].read!=query.read || reusableNodes[j].start<query.start-Intervals.absoluteThreshold) break;
			if (Intervals.areApproximatelyIdentical(reusableNodes[j].start,reusableNodes[j].end,query.start,query.end)) {
				Node.order=previousOrder;
				return reusableNodes[j];
			}
			j--;
		}
		j=i+1;
		while (j<=lastReusableNode) {
			if (reusableNodes[j].read!=query.read || reusableNodes[j].start>=query.start) break;
			if (Intervals.areApproximatelyIdentical(reusableNodes[j].start,reusableNodes[j].end,query.start,query.end)) {
				Node.order=previousOrder;
				return reusableNodes[j];
			}
			j++;
		}
		Node.order=previousOrder;
		return null;
	}
	
	
	/**
	 * Assume that the maximal events on $source$ induced by its neighbors have peaks at 
	 * relative positions $peaks[0..lastPeak]$. The procedure considers all possible 
	 * undersplit nodes that conform to the sides of $events[0..lastEvent]$, and redirects
	 * to a shortest undersplit node $X$ each edge adjacent to $source$ whose overlap with 
	 * $source$ is limited to $X$. Thus, $source$ might have no adjacent edge after the 
	 * procedure completes.
	 *
	 * Remark: if $useAllUndersplitNodes=FALSE$, only undersplit nodes to which enough 
	 * neighbors of $source$ get reassigned are actually created. The threshold is set to 
	 * a percentile of the set of numbers of assigned neighbors to each undersplit node.
	 *
	 * Remark: edges that correspond to more than one alignment are not redirected.
	 *
	 * Remark: the procedure does not necessarily build undersplit nodes, but tries first
	 * to reuse nodes in $reusableNodes$ that are sufficiently similar to the desired 
	 * intervals. For simplicity, such nodes are assumed not to be affected by undersplit.
	 *
	 * Remark: the procedure assumes that the $periodicSubstrings$ list of every node has 
	 * already been built.
	 *
	 * @param peaks sorted, relative positions, including 0 and $source.length()-1$;
	 * @param events sorted relative positions; each event is followed by a side tag: 
	 * 0=left, 1=right, 2=left and right;
	 * @param fromNode stores newly-created undersplit nodes in $nodesArray$ starting from
	 * $fromNode$ (inclusive);
	 * @param query temporary space;
	 * @param stats number of: 0=nodes created; 1=nodes reused; 2=edges redirected;
	 * @param peaksLeft,peaksRight,undersplitNodes,undersplitNodeIsNew temporary space;
	 * @param undersplitNodeIsUsed,undersplitNodeIsUsedPrime temporary space;
	 * @return the position in $nodesArray$ of the last undersplit node created by the
	 * procedure, if any.
	 */
	private static final int fixUndersplits_impl(Node source, int[] peaks, int lastPeak, int[] events, int lastEvent, Node[] reusableNodes, int lastReusableNode, int fromNode, Node query, int[] stats, boolean[] peaksLeft, boolean[] peaksRight, Node[] undersplitNodes, boolean[] undersplitNodeIsNew, int[] undersplitNodeIsUsed, int[] undersplitNodeIsUsedPrime, AlignmentPair tmpPair, boolean useAllUndersplitNodes) {
		final int GROWTH_RATE = 1000;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_NEIGHBORS_THRESHOLD = 3;  // Arbitrary
		final int NEIGHBORS_THRESHOLD_QUARTILE = 3;  // Arbitrary
		final int start = source.start;
		final int read = source.read;
		boolean isPeriodic, nodesReused;
		int i, j, k;
		int to, out, end, lastUndersplitNode, firstJForNextI, checked, nNodesWithEdges;
		int neighborsThreshold, periodicSubstringType;
		Node tmpNode;
		Edge tmpEdgePrime;
		int[] newNNeighbors;
		Node[] newNodesArray;
		Edge[][] newNeighbors;
		
		// Building all possible undersplit nodes that conform to $events$.
		Math.set(peaksLeft,lastPeak,false); Math.set(peaksRight,lastPeak,false);
		peaksRight[0]=true; peaksLeft[lastPeak]=true;
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastEvent && j<=lastPeak) {
			if (peaks[j]<events[i]-IDENTITY_THRESHOLD) {
				j++;
				continue;
			}
			if (peaks[j]>events[i]+IDENTITY_THRESHOLD) {
				i+=2;
				if (firstJForNextI!=-1) {
					j=firstJForNextI;
					firstJForNextI=-1;
				}
				continue;
			}
			if (firstJForNextI==-1 && i<lastEvent && peaks[j]>=events[i+1]-IDENTITY_THRESHOLD) firstJForNextI=j;
			if (events[i+1]==0 || events[i+1]==2) peaksLeft[j]=true;
			if (events[i+1]==1 || events[i+1]==2) peaksRight[j]=true;
			j++;
		}
		lastUndersplitNode=-1;
		query.read=source.read;
		for (i=0; i<lastPeak; i++) {
			if (!peaksRight[i]) continue;
			query.start=source.start+peaks[i];
			end=i==0?lastPeak-1:lastPeak;
			for (j=i+1; j<=end; j++) {
				if (!peaksLeft[j]) continue;
				query.end=source.start+peaks[j];
				tmpNode=getReusableNode(query,reusableNodes,lastReusableNode);
				if (tmpNode==null) {
					isPeriodic=false;
					periodicSubstringType=source.inPeriodicSubstring(peaks[i],peaks[j]);
					if (periodicSubstringType==1) {
						isPeriodic=true;
						lastUndersplitNode++;
						undersplitNodes[lastUndersplitNode] = new Node(source.read,Constants.INTERVAL_PERIODIC,-1,source.start+peaks[i],source.start+peaks[j],i==0?source.isLeftMaximal:true,j==lastPeak?source.isRightMaximal:true,false,source.isContainedInterval,-1,-1,0,false,-1,-1,0,0);
						undersplitNodeIsNew[lastUndersplitNode]=true;						
					}
					else if (periodicSubstringType==2) {
						k=getPeriodicInterval(source.start+peaks[i],source.start+peaks[j],source.read,source.nodeID,true,2);
						if (k!=-1 && nodesArray[k].period>0 && peaks[j]-peaks[i]+1>=nodesArray[k].period) {
							isPeriodic=true;
							lastUndersplitNode++;
							undersplitNodes[lastUndersplitNode] = new Node(source.read,Constants.INTERVAL_PERIODIC,-1,source.start+peaks[i],source.start+peaks[j],i==0?source.isLeftMaximal:true,j==lastPeak?source.isRightMaximal:true,false,source.isContainedInterval,-1,-1,0,true,-1,-1,0,0);
							undersplitNodeIsNew[lastUndersplitNode]=true;							
						}
					}
					if (!isPeriodic) {
						lastUndersplitNode++;
						undersplitNodes[lastUndersplitNode] = new Node(source.read,Constants.INTERVAL_ALIGNMENT,-1,source.start+peaks[i],source.start+peaks[j],i==0?source.isLeftMaximal:true,j==lastPeak?source.isRightMaximal:true,false,source.isContainedInterval,-1,-1,0,false,-1,-1,0,0);
						undersplitNodeIsNew[lastUndersplitNode]=true;
					}
					stats[0]++;
				}
				else if (tmpNode!=source) {
					lastUndersplitNode++;
					undersplitNodes[lastUndersplitNode]=tmpNode;
					undersplitNodeIsNew[lastUndersplitNode]=false;
					stats[1]++;
				}
			}
		}
		if (lastUndersplitNode==-1) return fromNode-1;
		
		// Removing duplicate potential undersplit nodes (induced by assigning the same
	 	// reusable node to multiple pairs of peaks).
		// Remark: sorting $undersplitNodes$ is not likely to be useful, since the array
		// is assumed to contain few elements.
		for (i=0; i<lastUndersplitNode; i++) {
			if (undersplitNodes[i]==null) continue;
			for (j=i+1; j<=lastUndersplitNode; j++) {
				if (undersplitNodes[j]==null) continue;
				if (undersplitNodes[j]==undersplitNodes[i]) undersplitNodes[j]=null;
			}
		}
		j=0;
		for (i=1; i<=lastUndersplitNode; i++) {
			if (undersplitNodes[i]==null) continue;
			j++;
			undersplitNodes[j]=undersplitNodes[i];
			undersplitNodeIsNew[j]=undersplitNodeIsNew[i];
		}
		lastUndersplitNode=j;
		if (lastUndersplitNode+1<undersplitNodes.length) undersplitNodes[lastUndersplitNode+1]=null;  // End marker, used by the calling procedure.
		
		// Filling $undersplitNodeIsUsed$ and estimating $neighborsThreshold$.
		Math.set(undersplitNodeIsUsed,lastUndersplitNode,0);
		for (j=0; j<nNeighbors[source.nodeID]; j++) fixUndersplits_handleEdge(true,source,start,j,neighbors[source.nodeID][j],tmpEdge,undersplitNodes,lastUndersplitNode,undersplitNodeIsUsed,-1,null,null);
		nNodesWithEdges=0; j=-1;
		for (i=0; i<=lastUndersplitNode; i++) {
			if (undersplitNodeIsUsed[i]>0) {
				nNodesWithEdges++;
				j=i;
			}
		}
		if (nNodesWithEdges==1 && undersplitNodeIsUsed[j]==nNeighbors[source.nodeID]) return fromNode-1;
		if (useAllUndersplitNodes) neighborsThreshold=1;
		else {
			if (lastUndersplitNode<3) neighborsThreshold=MIN_NEIGHBORS_THRESHOLD;
			else {
				System.arraycopy(undersplitNodeIsUsed,0,undersplitNodeIsUsedPrime,0,lastUndersplitNode+1);
				if (lastUndersplitNode>0) Arrays.sort(undersplitNodeIsUsedPrime,0,lastUndersplitNode+1);
				j=0;
				while (j<=lastUndersplitNode && undersplitNodeIsUsedPrime[j]==0) j++;
				if (j==lastUndersplitNode+1) neighborsThreshold=MIN_NEIGHBORS_THRESHOLD;
				else {
					i=(NEIGHBORS_THRESHOLD_QUARTILE*(lastUndersplitNode+1-j))>>2;
					neighborsThreshold=Math.max(MIN_NEIGHBORS_THRESHOLD,undersplitNodeIsUsedPrime[j+i]);
				}
			}
		}
		for (i=0; i<=lastUndersplitNode; i++) {
			if (undersplitNodeIsUsed[i]<neighborsThreshold) undersplitNodeIsUsed[i]=0;
		}
		
		// Adding to $nodesArray$ just undersplit nodes with enough neighbors
		out=fromNode-1; nodesReused=false;
		for (i=0; i<=lastUndersplitNode; i++) {
			if (undersplitNodeIsUsed[i]<neighborsThreshold) continue;
			if (!undersplitNodeIsNew[i]) {
				nodesReused=true;
				continue;
			}
			out++;
			expandTo(out);
			undersplitNodes[i].nodeID=out;
			nodesArray[out]=undersplitNodes[i];
			nNeighbors[out]=0;
		}
		if (out==fromNode-1 && !nodesReused) return out;		
		
		// Transforming edges from $source$
		for (j=0; j<nNeighbors[source.nodeID]; j++) fixUndersplits_handleEdge(false,source,start,j,neighbors[source.nodeID][j],tmpEdge,undersplitNodes,lastUndersplitNode,undersplitNodeIsUsed,neighborsThreshold,stats,tmpPair);
		
		// Compacting the remaining neighbors of $source$, if any.
		i=0;
		for (j=0; j<nNeighbors[source.nodeID]; j++) {
			if (neighbors[source.nodeID][j]!=null) {
				if (i!=j) neighbors[source.nodeID][i]=neighbors[source.nodeID][j];
				i++;
			}
		}
		nNeighbors[source.nodeID]=i;		
		
		return out;
	}
	
	
	/**
	 * Remark: the procedure tries to redirect edges with: (1) one alignment; (2) many
	 * alignments, but they are all overlaps of the same overlap type; (3) many 
	 * alignments, reconstructing and redirecting each one of them, one by one.
	 *
	 * Remark: the procedure does not handle the case in which alignments of different
	 * overlap types could be extracted. This case is not likely to be frequent in 
	 * practice.
	 *
	 * Remark: the procedure assumes $sourceEdge.avgDiffs$ to be equal to the sum of
	 * alignment diffs, not the average.
	 *
	 * @param mode TRUE=checks whether $edge$ can be transformed, and increments counts in
	 * $undersplitNodeIsUsed$; FALSE=applies the transformation, reading the counts in
	 * $undersplitNodeIsUsed$;
	 * @param sourceColumn used just when $mode=false$;
	 * @param neighborsThreshold used just when $mode=false$;
	 * @param stats used just when $mode=false$;
	 * @param tmpEdge temporary space;
	 * @param tmpPair temporary space, used just when $mode=false$;
	 */
	private static final void fixUndersplits_handleEdge(boolean mode, Node source, int start, int sourceColumn, Edge sourceEdge, Edge tmpEdge, Node[] undersplitNodes, int lastUndersplitNode, int[] undersplitNodeIsUsed, int neighborsThreshold, int[] stats, AlignmentPair tmpPair) {
		if (sourceEdge==null || sourceEdge.isSelfLoop() || sourceEdge.sameRead()) return;
		final int sourceID = source.nodeID;
		final int to = sourceEdge.getTo(sourceID);
		final int nAlignments = sourceEdge.getNAlignments();
		final double diffsPerAlignment = sourceEdge.avgDiffs/nAlignments;  // Arbitrary
		boolean extracted;
		int j, k, n;
		int nExtracted;
		
		// One alignment
		if (nAlignments==1) {
			k=getShortestUndersplitNode(source,sourceEdge,undersplitNodes,lastUndersplitNode);
			if (k==-1) return;
			if (mode) {
				if (sourceEdge.canBeTransformed(sourceID,undersplitNodes[k].nodeID,undersplitNodes[k].start-start,undersplitNodes[k].end-start,undersplitNodes[k].type)) undersplitNodeIsUsed[k]++;
			}
			else {
				if (undersplitNodeIsUsed[k]<neighborsThreshold) return;
				sourceEdge.clone(tmpEdge);  // Backing up the original edge
				if (sourceEdge.transform(sourceID,undersplitNodes[k].nodeID,undersplitNodes[k].start-start,undersplitNodes[k].end-start,tmpPair)) {
					// $IntervalGraphStep2.inStep2_checkEdge()$ is already called inside
					// $Edge.transform()$.
					mergeIdenticalEdges(to,sourceEdge);
					moveEdge(sourceID,sourceColumn,undersplitNodes[k].nodeID);
					stats[2]++;
				}
				else {
					// Restoring the original edge
					tmpEdge.clone(sourceEdge);
				}
			}
			return;
		}
		
		// Multiple alignments
		if ( sourceEdge.orientation==-1 || 
			 ( sourceEdge.orientation==2 && 	 
			   !(sourceEdge.getNAlignments()==2 && sourceEdge.isType_overlap() && sourceEdge.nDistinctTypes()<=2)
			 )
		   ) return;
		
		// Special case: just multiple overlaps of a single overlap type.
		if (sourceEdge.nDistinctTypes()==1 && sourceEdge.nOverlapTypes()==1) {
			tmpEdge.clear();
			tmpEdge.overlap=sourceEdge.overlap;
			if (tmpEdge.overhangs==null) tmpEdge.overhangs = new int[4][2];
			for (j=0; j<tmpEdge.overhangs.length; j++) {
				// $maxOverhangs$ stores longest overhangs = shortest overlaps.
				System.arraycopy(sourceEdge.overhangs[j],0,tmpEdge.overhangs[j],0,sourceEdge.overhangs[j].length);
			}
			fixUndersplits_handleEdge_impl_overlaps(mode,source,start,sourceColumn,sourceEdge,tmpEdge,undersplitNodes,lastUndersplitNode,undersplitNodeIsUsed,neighborsThreshold,stats,tmpPair,diffsPerAlignment);
			return;
		}
		
		// General case: trying to transform every alignment that can be extracted from
	 	// $sourceEdge$.
		nExtracted=0;
		if (sourceEdge.nAlignmentsOfType(0)==1) {  // Containment
			tmpEdge.clear();
			tmpEdge.orientation=orientationOfExtractedEdge_nonOverlap(sourceEdge);
			tmpEdge.containment=sourceEdge.containment;
			tmpEdge.containmentSubtype=sourceEdge.containmentSubtype;
			tmpEdge.containmentOverhangLeft=sourceEdge.containmentOverhangLeft;
			tmpEdge.containmentOverhangRight=sourceEdge.containmentOverhangRight;
			if (fixUndersplits_handleEdge_impl(mode,source,start,sourceColumn,sourceEdge,tmpEdge,undersplitNodes,lastUndersplitNode,undersplitNodeIsUsed,neighborsThreshold,stats,tmpPair,diffsPerAlignment)) {
				nExtracted++;
				if (!mode) {
					if (sourceEdge.orientation==2 && sourceEdge.nDistinctTypes()==2 && sourceEdge.getNAlignments()==2 && sourceEdge.isType_overlap()) sourceEdge.orientation=Constants.overlap2orientation(sourceEdge.overlap);
					if (sourceEdge.nAlignmentsForward>0) sourceEdge.nAlignmentsForward--;
					else sourceEdge.nAlignmentsBackward--;
					sourceEdge.avgDiffs-=diffsPerAlignment;
					sourceEdge.setType_noContainment();
				}
			}
		}
		if (sourceEdge.nAlignmentsOfType(1)==1) {  // Overlap
			tmpEdge.clear();
			tmpEdge.orientation=Constants.overlap2orientation(sourceEdge.overlap);
			tmpEdge.overlap=sourceEdge.overlap;
			if (tmpEdge.overhangs==null) tmpEdge.overhangs = new int[4][2];
			for (j=0; j<tmpEdge.overhangs.length; j++) {
				System.arraycopy(sourceEdge.overhangs[j],0,tmpEdge.overhangs[j],0,sourceEdge.overhangs[j].length);
			}
			if (fixUndersplits_handleEdge_impl(mode,source,start,sourceColumn,sourceEdge,tmpEdge,undersplitNodes,lastUndersplitNode,undersplitNodeIsUsed,neighborsThreshold,stats,tmpPair,diffsPerAlignment)) {
				nExtracted++;
				if (!mode) {
					if (sourceEdge.orientation==2 && sourceEdge.nDistinctTypes()==2 && sourceEdge.getNAlignments()==2) sourceEdge.orientation=orientationOfExtractedEdge_nonOverlap(sourceEdge);
					if (sourceEdge.nAlignmentsForward>0) sourceEdge.nAlignmentsForward--;
					else sourceEdge.nAlignmentsBackward--;
					sourceEdge.avgDiffs-=diffsPerAlignment;
					sourceEdge.setType_noOverlap();
				}
			}
		}
		if (sourceEdge.nAlignmentsOfType(2)==1) {  // Insertion
			tmpEdge.clear();
			tmpEdge.orientation=orientationOfExtractedEdge_nonOverlap(sourceEdge);
			tmpEdge.insertion=sourceEdge.insertion;
			tmpEdge.containmentOverhangLeft=sourceEdge.containmentOverhangLeft;
			tmpEdge.containmentOverhangRight=sourceEdge.containmentOverhangRight;			
			if (fixUndersplits_handleEdge_impl(mode,source,start,sourceColumn,sourceEdge,tmpEdge,undersplitNodes,lastUndersplitNode,undersplitNodeIsUsed,neighborsThreshold,stats,tmpPair,diffsPerAlignment)) {
				nExtracted++;
				if (!mode) {
					if (sourceEdge.orientation==2 && sourceEdge.nDistinctTypes()==2 && sourceEdge.getNAlignments()==2 && sourceEdge.isType_overlap()) sourceEdge.orientation=Constants.overlap2orientation(sourceEdge.overlap);
					if (sourceEdge.nAlignmentsForward>0) sourceEdge.nAlignmentsForward--;
					else sourceEdge.nAlignmentsBackward--;
					sourceEdge.avgDiffs-=diffsPerAlignment;
					sourceEdge.setType_noInsertion();
				}
			}
		}
		if (sourceEdge.nAlignmentsOfType(3)==1) {  // Shared substring
			tmpEdge.clear();
			tmpEdge.orientation=orientationOfExtractedEdge_nonOverlap(sourceEdge);
			tmpEdge.sharedSubstring=sourceEdge.sharedSubstring;
			if (tmpEdge.sharedSubstringOverhangs==null) tmpEdge.sharedSubstringOverhangs = new int[4];
			System.arraycopy(sourceEdge.sharedSubstringOverhangs,0,tmpEdge.sharedSubstringOverhangs,0,4);
			if (fixUndersplits_handleEdge_impl(mode,source,start,sourceColumn,sourceEdge,tmpEdge,undersplitNodes,lastUndersplitNode,undersplitNodeIsUsed,neighborsThreshold,stats,tmpPair,diffsPerAlignment)) {
				nExtracted++;
				if (!mode) {
					if (sourceEdge.orientation==2 && sourceEdge.nDistinctTypes()==2 && sourceEdge.getNAlignments()==2 && sourceEdge.isType_overlap()) sourceEdge.orientation=Constants.overlap2orientation(sourceEdge.overlap);
					if (sourceEdge.nAlignmentsForward>0) sourceEdge.nAlignmentsForward--;
					else sourceEdge.nAlignmentsBackward--;
					sourceEdge.avgDiffs-=diffsPerAlignment;
					sourceEdge.setType_noSharedSubstring();
				}
			}
		}
		if (!mode) {
			// Removing $sourceEdge$ from $neighbors$.
			if (nExtracted==nAlignments) {
				n=nNeighbors[sourceID];
				for (j=0; j<n; j++) {
					if (neighbors[sourceID][j]==sourceEdge) {
						neighbors[sourceID][j]=null;
						break;
					}
				}
				n=nNeighbors[to];
				for (j=0; j<n; j++) {
					if (neighbors[to][j]==sourceEdge) {
						for (k=j+1; k<n; k++) neighbors[to][k-1]=neighbors[to][k];
						nNeighbors[to]--;
						break;
					}
				}
			}
		}
	}
	
	
	/**
	 * @return the orientation of an edge, of type different from overlap, extracted from 
	 * $sourceEdge$: 0=same orientation; 1=opposite; -1=unknown.
	 */
	private static final int orientationOfExtractedEdge_nonOverlap(Edge sourceEdge) {
		int otherOrientation;
		
		if (sourceEdge.orientation<=1) return sourceEdge.orientation;
		if (sourceEdge.nDistinctTypes()==2 && sourceEdge.getNAlignments()==2 && sourceEdge.isType_overlap()) {
			otherOrientation=Constants.overlap2orientation(sourceEdge.overlap);
			if (otherOrientation==0) return 1;
			else if (otherOrientation==1) return 0;
		}
		return -1;
	}
	
	
	/**
	 * @return if $mode=FALSE$: TRUE iff the transformation was successful and $newEdge$ 
	 * got disconnected from $source$; if $mode=TRUE$: TRUE iff $newEdge$ can be 
	 * transformed.
	 */
	private static final boolean fixUndersplits_handleEdge_impl(boolean mode, Node source, int start, int sourceColumn, Edge sourceEdge, Edge newEdge, Node[] undersplitNodes, int lastUndersplitNode, int[] undersplitNodeIsUsed, int neighborsThreshold, int[] stats, AlignmentPair tmpPair, double diffsPerAlignment) {
		int j, k, to;
		Edge createdEdge;
		
		newEdge.nodeID1=sourceEdge.nodeID1;
		newEdge.nodeID2=sourceEdge.nodeID2;
		newEdge.on=sourceEdge.on;
		newEdge.supplement=sourceEdge.supplement;
		newEdge.isRedundant=sourceEdge.isRedundant;
		newEdge.avgDiffs=diffsPerAlignment;
		if (newEdge.orientation==0) newEdge.nAlignmentsForward=1;
		else newEdge.nAlignmentsBackward=1;
		newEdge.implied=sourceEdge.implied;
		k=getShortestUndersplitNode(source,newEdge,undersplitNodes,lastUndersplitNode);
		if (k==-1) return false;
		if (mode) {
			if (newEdge.canBeTransformed(source.nodeID,undersplitNodes[k].nodeID,undersplitNodes[k].start-start,undersplitNodes[k].end-start,undersplitNodes[k].type)) {
				undersplitNodeIsUsed[k]++;
				return true;
			}
		}
		else {
			to=newEdge.getTo(source.nodeID);
			if (undersplitNodeIsUsed[k]<neighborsThreshold) return false;
			if (newEdge.transform(source.nodeID,undersplitNodes[k].nodeID,undersplitNodes[k].start-start,undersplitNodes[k].end-start,tmpPair)) {
				// $IntervalGraphStep2.inStep2_checkEdge()$ is already called inside
				// $Edge.transform()$.
				createdEdge = new Edge();
				newEdge.clone(createdEdge);
				addEdge(to,createdEdge);
				j=mergeIdenticalEdges(to,createdEdge);
				moveEdge(to,j,undersplitNodes[k].nodeID);
				neighbors[to][j]=createdEdge;  // Undoing the effect of $moveEdge()$
				stats[2]++;				
				return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Variant of $fixUndersplits_handleEdge_impl()$ for the case in which $sourceEdge$
	 * contains just multiple overlap edges of one overlap type. If the deepest overlap in
	 * $source$ is contained in a valid undersplit node of $source$, then $sourceEdge$ is 
	 * redirected from $source$ to the undersplit node, and its overhang values are 
	 * updated.
	 *
	 * Remark: the procedure could be generalized to the case in which all overlaps in 
	 * $sourceEdge$ use the same side of $source$, but not necessarily the same side of 
	 * the other node.
	 */
	private static final boolean fixUndersplits_handleEdge_impl_overlaps(boolean mode, Node source, int start, int sourceColumn, Edge sourceEdge, Edge newEdge, Node[] undersplitNodes, int lastUndersplitNode, int[] undersplitNodeIsUsed, int neighborsThreshold, int[] stats, AlignmentPair tmpPair, double diffsPerAlignment) {
		final int oldLength = source.length();
		final int row = Constants.overlap2id(sourceEdge.overlap);
		final int to = sourceEdge.getTo(source.nodeID);
		int k, newLength;
		Edge createdEdge;
		
		newEdge.nodeID1=sourceEdge.nodeID1;
		newEdge.nodeID2=sourceEdge.nodeID2;
		newEdge.on=sourceEdge.on;
		newEdge.supplement=sourceEdge.supplement;
		newEdge.isRedundant=sourceEdge.isRedundant;
		newEdge.avgDiffs=diffsPerAlignment;
		newEdge.orientation=sourceEdge.orientation;
		if (newEdge.orientation==0) newEdge.nAlignmentsForward=1;
		else newEdge.nAlignmentsBackward=1;
		newEdge.implied=sourceEdge.implied;
		k=getShortestUndersplitNode(source,newEdge,undersplitNodes,lastUndersplitNode);
		if (k==-1) return false;
		if (mode) {
			undersplitNodeIsUsed[k]++;
			return true;
		}
		else {
			if (undersplitNodeIsUsed[k]<neighborsThreshold) return false;
			// Checking whether the redirected $newEdge$ is valid for Step 2.
			newLength=undersplitNodes[k].length();
			if (newEdge.nodeID1==source.nodeID) {
				newEdge.nodeID1=undersplitNodes[k].nodeID;
				newEdge.overhangs[row][0]=Math.max(0,newLength-(oldLength-newEdge.overhangs[row][0]));
			}
			else {
				newEdge.nodeID2=undersplitNodes[k].nodeID;
				newEdge.overhangs[row][1]=Math.max(0,newLength-(oldLength-newEdge.overhangs[row][1]));
			}
			if (IntervalGraphStep2.inStep2_checkEdge(newEdge)==-1) return false;
			// Redirecting $sourceEdge$
			if (sourceEdge.nodeID1==source.nodeID) {
				sourceEdge.nodeID1=undersplitNodes[k].nodeID;
				sourceEdge.overhangs[row][0]=Math.max(0,newLength-(oldLength-sourceEdge.overhangs[row][0]));
				if (sourceEdge.maxOverhangs!=null) sourceEdge.maxOverhangs[row][0]=Math.max(0,newLength-(oldLength-sourceEdge.maxOverhangs[row][0]));
			}
			else {
				sourceEdge.nodeID2=undersplitNodes[k].nodeID;
				sourceEdge.overhangs[row][1]=Math.max(0,newLength-(oldLength-sourceEdge.overhangs[row][1]));
				if (sourceEdge.maxOverhangs!=null) sourceEdge.maxOverhangs[row][1]=Math.max(0,newLength-(oldLength-sourceEdge.maxOverhangs[row][1]));
			}
			mergeIdenticalEdges(to,sourceEdge);
			moveEdge(source.nodeID,sourceColumn,undersplitNodes[k].nodeID);
			stats[2]++;
			return true;
		}
	}
	
	
	/**
	 * @return the position in $undersplitNodes$ of a shortest node such that edge $edge$
	 * adjacent to $oldNode$ can be transformed into an edge to the undersplit node;
	 * returns -1 if no such undersplit node can be found.
	 */
	private static final int getShortestUndersplitNode(Node oldNode, Edge edge, Node[] undersplitNodes, int lastUndersplitNode) {
		int i;
		int length, min, minLength;
		
		min=-1; minLength=Math.POSITIVE_INFINITY;
		for (i=0; i<=lastUndersplitNode; i++) {			
			if (edge.canBeTransformed(oldNode.nodeID,undersplitNodes[i].nodeID,undersplitNodes[i].start-oldNode.start,undersplitNodes[i].end-oldNode.start,undersplitNodes[i].type)) {
			    length=undersplitNodes[i].length();
			    if (length<minLength) {
				    minLength=length;
				    min=i;
			    }
		    }
		}
		return min;
	}
	
	
	/**
	 * Moves edge $neighbors[sourceNode][column]$ to row $destinationNode$ of $neighbors$.
	 *
	 * Remark: if an edge identical to $neighbors[sourceNode][sourceColumn]$ is found in
	 * $neighbors[destinationNode]$, the procedure does not call $addTypes()$ and 
	 * $addDiffs()$, since those have already been called by $mergeIdenticalEdges()$.
	 */
	private static final void moveEdge(int sourceNode, int sourceColumn, int destinationNode) {
		boolean found;
		int j;
		
		found=false;
		for (j=0; j<nNeighbors[destinationNode]; j++) {
			if (neighbors[destinationNode][j]==null) continue;
			if (neighbors[destinationNode][j].equals(neighbors[sourceNode][sourceColumn])) {
				neighbors[destinationNode][j]=neighbors[sourceNode][sourceColumn];
				found=true;
				break;
			}
		}
		if (!found) addEdge(destinationNode,neighbors[sourceNode][sourceColumn]);
		neighbors[sourceNode][sourceColumn]=null;
	}
	
	
	/**
	 * Remark: the procedure does not assume $neighbors[nodeID]$ to be sorted, and it 
	 * considers all edges.
	 *
	 * @return position in $neighbors[nodeID]$ of an edge equal to $edge$, and of the 
	 * specified type; $Math.NEGATIVE_INFINITY$ if no edge equal to $edge$ and of the
	 * specified type is found.
	 */
	private static final int findEdge(int nodeID, Edge edge, boolean supplement) {
		for (int j=0; j<nNeighbors[nodeID]; j++) {
			if (neighbors[nodeID][j]==null || neighbors[nodeID][j].supplement!=supplement) continue;
			if (neighbors[nodeID][j].equals(edge)) return j;
		}
		return Math.NEGATIVE_INFINITY;
	}
	
	
	/**
	 * @return position of the first pointer to the $edge$ object in $neighbors[nodeID]$,
	 * or -1 if no such pointer is found.
	 */
	private static final int findEdgePointer(int nodeID, Edge edge) {
		for (int j=0; j<nNeighbors[nodeID]; j++) {
			if (neighbors[nodeID][j]==edge) return j;
		}
		return -1;
	}
	
	
	/**
	 * Remark: the procedure tries to add $edge$ to the first null slot in 
	 * $neighbors[nodeID]$. If such a null slot exists, $nNeighbors[nodeID]$ is not 
	 * incremented.
	 *
	 * Remark: the procedure does not remove null neighbors of $nodeID$. Such compaction
	 * should be perfomed by the caller.
	 *
	 * @return the position in $neighbors[nodeID]$ at which $edge$ is inserted.
	 */
	private static final int addEdge(int nodeID, Edge edge) {
		final int MIN_EDGES_PER_NODE = 10;
		int i, j, last;
		Edge[] tmpEdges;
		
		// Adding $edge$ to the first null slot, if any.
		for (j=0; j<nNeighbors[nodeID]; j++) {
			if (neighbors[nodeID][j]==null) {
				neighbors[nodeID][j]=edge;
				return j;
			}
		}
		
		// Adding $edge$ to the end of the array
		nNeighbors[nodeID]++;
		if (neighbors[nodeID]==null || neighbors[nodeID].length==0) neighbors[nodeID] = new Edge[MIN_EDGES_PER_NODE];
		else if (neighbors[nodeID].length<nNeighbors[nodeID]) {
			tmpEdges = new Edge[1+(neighbors[nodeID].length<<1)];
			System.arraycopy(neighbors[nodeID],0,tmpEdges,0,neighbors[nodeID].length);
			neighbors[nodeID]=tmpEdges;
		}
		neighbors[nodeID][nNeighbors[nodeID]-1]=edge;
		
		return nNeighbors[nodeID]-1;
	}
	
	
	/**
	 * Removes duplicates among the newly-created undersplit nodes, which are assumed to 
	 * be in $nodesArray[from..to]$. Duplicates can occur because e.g. two intervals in 
	 * the same read might be straddling, and their shared part might induce two new but 
	 * identical undersplit nodes.
	 *
	 * @return the absolute position of the last node in $nodesArray$ after merging.
	 */
	private static final int mergeUndersplitNodes(int from, int to) {
		boolean isPeriodic;
		int i, j, k, h, l;
		int read, start, end, sumStart, sumEnd, nMerged, newNode, otherNode, out;
		int hole, source, tmpNumber;
		Edge oldEdge;
		Edge[] tmpNeighbors;
		int[][] mergedCoordinates;
		
		// Sorting new nodes		
		Node.order=Node.READ_START_NODEID;
		if (to+1-from>1) Arrays.sort(nodesArray,from,to+1);
		Node.order=Node.UNSORTED;
		for (i=from; i<=to; i++) {
			source=nodesArray[i].nodeID;
			if (source==i) continue;
			tmpNeighbors=neighbors[i];
			tmpNumber=nNeighbors[i];
			hole=i;
			do {
				neighbors[hole]=neighbors[source];
				nNeighbors[hole]=nNeighbors[source];
				for (j=0; j<nNeighbors[hole]; j++) {
					if (neighbors[hole][j].nodeID1==source) neighbors[hole][j].nodeID1=hole;
					else neighbors[hole][j].nodeID2=hole;
				}
				nodesArray[hole].nodeID=hole;
				hole=source;
				source=nodesArray[hole].nodeID;
			} while (source!=i);
			neighbors[hole]=tmpNeighbors;
			nNeighbors[hole]=tmpNumber;
			for (j=0; j<nNeighbors[hole]; j++) {
				if (neighbors[hole][j].nodeID1==source) neighbors[hole][j].nodeID1=hole;
				else neighbors[hole][j].nodeID2=hole;
			}
			nodesArray[hole].nodeID=hole;
		}
		
		// Merging new nodes
		mergedCoordinates = new int[to-from+1][2];
		Math.set(mergedCoordinates,-1);
		for (i=from; i<=to; i++) nodesArray[i].visited=-1;
		for (i=from; i<=to; i++) {
			if (nodesArray[i].visited!=-1) continue;
			read=nodesArray[i].read;
			start=nodesArray[i].start;
			end=nodesArray[i].end;
			sumStart=start; sumEnd=end; nMerged=1;
			for (j=i+1; j<=to; j++) {
				if (nodesArray[j].read!=read || nodesArray[j].start>=end) break;
				if (!Intervals.areApproximatelyIdentical(start,end,nodesArray[j].start,nodesArray[j].end)) continue;
				nodesArray[j].visited=i;
				nMerged++;
				sumStart+=nodesArray[j].start; sumEnd+=nodesArray[j].end;
				if (nodesArray[j].type==Constants.INTERVAL_PERIODIC) {
					if (nodesArray[i].type!=Constants.INTERVAL_PERIODIC) {
						nodesArray[i].type=Constants.INTERVAL_PERIODIC;
						nodesArray[i].period=0;
						nodesArray[i].hasLongPeriod=nodesArray[j].hasLongPeriod;
					}
					else nodesArray[i].hasLongPeriod&=nodesArray[j].hasLongPeriod;
				}
				nodesArray[i].isLeftMaximal&=nodesArray[j].isLeftMaximal;
				nodesArray[i].isRightMaximal&=nodesArray[j].isRightMaximal;
				nodesArray[i].isContainedInterval|=nodesArray[j].isContainedInterval;
			}
			mergedCoordinates[i-from][0]=sumStart/nMerged;
			mergedCoordinates[i-from][1]=sumEnd/nMerged;
		}
		
		// Updating $neighbors$
		for (i=from; i<=to; i++) {
			newNode=nodesArray[i].visited;
			if (newNode==-1) continue;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null) continue;
				otherNode=neighbors[i][j].getTo(i);
				mergeUndersplitNodes_transformEdge(neighbors[i][j],i,newNode,mergedCoordinates[newNode-from][0],mergedCoordinates[newNode-from][1]);
				k=findEdge(newNode,neighbors[i][j],false);
				if (k==Math.NEGATIVE_INFINITY) addEdge(newNode,neighbors[i][j]);
				else {
					oldEdge=neighbors[newNode][k];
					neighbors[i][j].addTypes(oldEdge);
					neighbors[i][j].addDiffs(oldEdge);
					neighbors[i][j].addOrientation(oldEdge);
					neighbors[i][j].enforceOneAlignment();
					neighbors[newNode][k]=neighbors[i][j];
					k=findEdgePointer(otherNode,oldEdge);
					if (k==-1) {
						System.err.println("mergeUndersplitNodes> ERROR: incorrect edge pointers");
						System.exit(1);
					}
					neighbors[otherNode][k]=null;
				}
			}
		}
		
		// Compacting $nodesArray$ and $neighbors$.
		j=from-1;
		for (i=from; i<=to; i++) {
			if (nodesArray[i].visited!=-1) continue;
			nodesArray[i].start=mergedCoordinates[i-from][0];
			nodesArray[i].end=mergedCoordinates[i-from][1];
			j++;
			nodesArray[j]=nodesArray[i];
			nodesArray[j].nodeID=j;
			neighbors[j]=neighbors[i];
			nNeighbors[j]=nNeighbors[i];
			for (k=0; k<nNeighbors[j]; k++) {
				if (neighbors[j][k]==null) continue;
				if (neighbors[j][k].nodeID1==i) neighbors[j][k].nodeID1=j;
				else neighbors[j][k].nodeID2=j;
			}
		}
		for (i=j+1; i<=to; i++) {
			nodesArray[i]=null; neighbors[i]=null;
			nNeighbors[i]=0;
		}
		out=j;
		mergedCoordinates=null;
		
		// Removing null neighbors and self-loops from merged undersplit nodes
		for (i=from; i<=out; i++) {
			k=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j]==null || neighbors[i][j].getTo(i)==i) continue;
				neighbors[i][++k]=neighbors[i][j];
			}
			nNeighbors[i]=k+1;
		}
		
		// Removing undersplit nodes without incident edges (first time).
		j=from-1;
		for (i=from; i<=out; i++) {
			if (nNeighbors[i]==0) continue;
			j++;
			if (j!=i) {
				nodesArray[j]=nodesArray[i];
				nodesArray[j].nodeID=j;
				nNeighbors[j]=nNeighbors[i];
				neighbors[j]=neighbors[i];
				for (k=0; k<nNeighbors[j]; k++) {
					if (neighbors[j][k].nodeID1==i) neighbors[j][k].nodeID1=j;
					else neighbors[j][k].nodeID2=j;
				}
			}
		}
		for (i=j+1; i<=out; i++) {
			nodesArray[i]=null; neighbors[i]=null;
			nNeighbors[i]=0;
		}
		out=j;
		
		// Removing null neighbors from the entire $neighbors$.
		for (i=from; i<=out; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				otherNode=neighbors[i][j].getTo(i);
				if (otherNode>=from) continue;
				nodesArray[otherNode].visited=-1;
			}
		}
		for (i=from; i<=out; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				otherNode=neighbors[i][j].getTo(i);
				if (otherNode>=from || nodesArray[otherNode].visited!=-1) continue;
				k=-1;
				for (h=0; h<nNeighbors[otherNode]; h++) {
					if (neighbors[otherNode][h]==null) continue;
					neighbors[otherNode][++k]=neighbors[otherNode][h];
				}
				nNeighbors[otherNode]=k+1;
				nodesArray[otherNode].visited=0;
			}
		}
		
		// Enforcing that edges incident to periodic nodes must belong to Step 2.
		for (i=from; i<=out; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=0;
		}
		for (i=from; i<=out; i++) {
			isPeriodic=nodesArray[i].type==Constants.INTERVAL_PERIODIC;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].printed!=0) continue;
				otherNode=neighbors[i][j].getTo(i);
				if (!isPeriodic && nodesArray[otherNode].type!=Constants.INTERVAL_PERIODIC) {
					neighbors[i][j].printed=1;
					continue;
				}
				k=IntervalGraphStep2.inStep2_checkEdge(neighbors[i][j]);
				if (k==-1) {
					neighbors[i][j].printed=-1;
					if (otherNode<from) neighbors[otherNode][findEdgePointer(otherNode,neighbors[i][j])]=null;
				}
				else neighbors[i][j].printed=1;
			}
		}
		for (i=from; i<=out; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].printed!=-1) continue;
				otherNode=neighbors[i][j].getTo(i);
				if (otherNode<from) nodesArray[otherNode].visited=-1;
			}
		}
		for (i=from; i<=out; i++) {
			k=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].printed!=-1) {
					neighbors[i][++k]=neighbors[i][j];
					continue;
				}
				otherNode=neighbors[i][j].getTo(i);
				if (otherNode>=from || nodesArray[otherNode].visited!=-1) continue;
				l=-1;
				for (h=0; h<nNeighbors[otherNode]; h++) {
					if (neighbors[otherNode][h]==null) continue;
					neighbors[otherNode][++l]=neighbors[otherNode][h];
				}
				nNeighbors[otherNode]=l+1;
				nodesArray[otherNode].visited=0;
			}
			nNeighbors[i]=k+1;
		}

		// Removing undersplit nodes without incident edges (second time).
		j=from-1;
		for (i=from; i<=out; i++) {
			if (nNeighbors[i]==0) continue;
			j++;
			if (j!=i) {
				nodesArray[j]=nodesArray[i];
				nodesArray[j].nodeID=j;
				nNeighbors[j]=nNeighbors[i];
				neighbors[j]=neighbors[i];
				for (k=0; k<nNeighbors[j]; k++) {
					if (neighbors[j][k].nodeID1==i) neighbors[j][k].nodeID1=j;
					else neighbors[j][k].nodeID2=j;
				}
			}
		}
		for (i=j+1; i<=out; i++) {
			nodesArray[i]=null; neighbors[i]=null;
			nNeighbors[i]=0;
		}
		out=j;
		
		return out;
	}
	
	
	/**
	 * Remark: the procedure assumes that $edge$ has just one type in one orientation, and
	 * that $oldNode$ and $newNode$ come from the same read.
	 */
	private static final void mergeUndersplitNodes_transformEdge(Edge edge, int oldNode, int newNode, int newNodeStart, int newNodeEnd) {
		if (!edge.hasOneType()) {
			System.err.println("mergeUndersplitNodes_transformEdge> ERROR: the edge does not have just one type in one orientation");	
			System.err.println("nAlignments="+edge.getNAlignments()+" containment="+edge.containment+" overlap="+edge.overlap+" insertion="+edge.insertion+" sharedSubstring="+edge.sharedSubstring+" orientation="+edge.orientation);
			System.err.println("edge: "+edge);
			System.err.println("node1: "+edge.nodeID1+": "+nodesArray[edge.nodeID1]);
			System.err.println("node2: "+edge.nodeID2+": "+nodesArray[edge.nodeID2]);
			System.exit(1);
		}
		
		mergeUndersplitNodes_transformEdge_impl(edge,oldNode,newNode,newNodeStart,newNodeEnd);
		if (edge.nodeID1==oldNode) edge.nodeID1=newNode;
		else edge.nodeID2=newNode;
	}
	
	
	/**
	 * After replacing $oldNode$ with $newNode$, which is either a shrinked or an expanded
	 * version of $oldNode$: (1) a containment edge can become an identity edge; (2) an 
	 * insertion edge can become a containment edge or an identity edge; (3) an overlap 
	 * edge can become a containment edge; (4) a shared substring edge can become a 
	 * containment edge or an overlap edge.
	 */
	private static final void mergeUndersplitNodes_transformEdge_impl(Edge edge, int oldNode, int newNode, int newNodeStart, int newNodeEnd) {
		int i, j;
		int otherNode, oldNodeStart, oldNodeEnd;		
		
		oldNodeStart=nodesArray[oldNode].start;
		oldNodeEnd=nodesArray[oldNode].end;
		otherNode=edge.getTo(oldNode);
		if (oldNode==edge.nodeID1) {
			if (edge.isType(Constants.CONTAINMENT_TWO_IN_ONE) || edge.isType(Constants.INSERTION_TWO_IN_ONE)) {
				if (newNodeStart>oldNodeStart) {
					if (newNodeStart<=oldNodeStart+edge.containmentOverhangLeft) edge.containmentOverhangLeft-=newNodeStart-oldNodeStart;
					else edge.containmentOverhangLeft=0;
				}
				else if (newNodeStart<oldNodeStart) edge.containmentOverhangLeft+=oldNodeStart-newNodeStart;
				if (newNodeEnd>oldNodeEnd) edge.containmentOverhangRight+=newNodeEnd-oldNodeEnd;
				else if (newNodeEnd<oldNodeEnd) {
					if (newNodeEnd>=oldNodeEnd-edge.containmentOverhangRight) edge.containmentOverhangRight-=oldNodeEnd-newNodeEnd;
					else edge.containmentOverhangRight=0;
				}
				if (edge.containmentOverhangLeft==0 && edge.containmentOverhangRight==0) {
					edge.setType_noContainment();
					edge.setType_noInsertion();
					edge.setType(Constants.CONTAINMENT_IDENTICAL,true);
				}
				else if (edge.containmentOverhangLeft==0) {
					if ( edge.isType(Constants.INSERTION_TWO_IN_ONE) &&
				         ( (edge.orientation==0 && !nodesArray[otherNode].isRightMaximal) ||
						   (edge.orientation==1 && !nodesArray[otherNode].isLeftMaximal)
						 )
				       ) {
					   edge.setType_noInsertion();
					   edge.setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
					}
					edge.containmentSubtype=Constants.CONTAINMENT_PREFIX;
				}
				else if (edge.containmentOverhangRight==0) {
					if ( edge.isType(Constants.INSERTION_TWO_IN_ONE) &&
				         ( (edge.orientation==0 && !nodesArray[otherNode].isLeftMaximal) ||
						   (edge.orientation==1 && !nodesArray[otherNode].isRightMaximal)
						 )
				       ) {
					   edge.setType_noInsertion();
					   edge.setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
					}
					edge.containmentSubtype=Constants.CONTAINMENT_SUFFIX;
				}
			}
			else if (edge.isType_overlap() && edge.overhangs!=null) {
				j=0;
				for (i=Constants.OVERLAP_PREFIX_PREFIX; i<=Constants.OVERLAP_SUFFIX_SUFFIX; i<<=1) {
					if (!edge.isType(i)) {
						j++;
						continue;
					}
					if (i==Constants.OVERLAP_PREFIX_PREFIX || i==Constants.OVERLAP_PREFIX_SUFFIX) {
						if (newNodeEnd>oldNodeEnd) edge.overhangs[j][0]+=newNodeEnd-oldNodeEnd;
						else if (newNodeEnd<oldNodeEnd) {
							if (newNodeEnd>=oldNodeEnd-edge.overhangs[j][0]) edge.overhangs[j][0]-=oldNodeEnd-newNodeEnd;
							else edge.overhangs[j][0]=0;
						}
					}
					else if (i==Constants.OVERLAP_SUFFIX_PREFIX || i==Constants.OVERLAP_SUFFIX_SUFFIX) {
						if (newNodeStart>oldNodeStart) {
							if (newNodeStart<=oldNodeStart+edge.overhangs[j][0]) edge.overhangs[j][0]-=newNodeStart-oldNodeStart;
							else edge.overhangs[j][0]=0;
						}
						else if (newNodeStart<oldNodeStart) edge.overhangs[j][0]+=oldNodeStart-newNodeStart;
					}
					if (edge.overhangs[j][0]==0) {
						edge.setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
						if (i==Constants.OVERLAP_SUFFIX_PREFIX || i==Constants.OVERLAP_PREFIX_PREFIX) {
						   	edge.containmentSubtype=Constants.CONTAINMENT_PREFIX;
						   	edge.containmentOverhangLeft=0;
						   	edge.containmentOverhangRight=edge.overhangs[j][1];
						}
						else if (i==Constants.OVERLAP_PREFIX_SUFFIX || i==Constants.OVERLAP_SUFFIX_SUFFIX) {
							edge.containmentSubtype=Constants.CONTAINMENT_SUFFIX;
							edge.containmentOverhangLeft=edge.overhangs[j][1];
							edge.containmentOverhangRight=0;
						}
						edge.setType_noOverlap();
						break;
					}
					j++;
				}
			}
			else if (edge.isType(Constants.SHARED_SUBSTRING) && edge.sharedSubstringOverhangs!=null) {
				if (newNodeStart>oldNodeStart) {
					if (newNodeStart<=oldNodeStart+edge.sharedSubstringOverhangs[0]) edge.sharedSubstringOverhangs[0]-=newNodeStart-oldNodeStart;
					else edge.sharedSubstringOverhangs[0]=0;
				}
				else if (newNodeStart<oldNodeStart) edge.sharedSubstringOverhangs[0]+=oldNodeStart-newNodeStart;
				if (newNodeEnd>oldNodeEnd) edge.sharedSubstringOverhangs[1]+=newNodeEnd-oldNodeEnd;
				else if (newNodeEnd<oldNodeEnd) {
					if (newNodeEnd>=oldNodeEnd-edge.sharedSubstringOverhangs[1]) edge.sharedSubstringOverhangs[1]-=oldNodeEnd-newNodeEnd;
					else edge.sharedSubstringOverhangs[1]=0;
				}
				if (edge.sharedSubstringOverhangs[0]==0 && edge.sharedSubstringOverhangs[1]==0) {
					edge.setType(Constants.SHARED_SUBSTRING,false);
					edge.setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
					edge.containmentOverhangLeft=edge.sharedSubstringOverhangs[2];
					edge.containmentOverhangRight=edge.sharedSubstringOverhangs[3];
					if (edge.containmentOverhangLeft<=IDENTITY_THRESHOLD) edge.containmentSubtype=Constants.CONTAINMENT_PREFIX;
					else if (edge.containmentOverhangRight<=IDENTITY_THRESHOLD) edge.containmentSubtype=Constants.CONTAINMENT_SUFFIX;
					else edge.containmentSubtype=Constants.CONTAINMENT_SUBSTRING;				
				}
				else if (edge.sharedSubstringOverhangs[0]==0) {
					if (edge.sharedSubstringOverhangs[2]<=IDENTITY_THRESHOLD) {
						edge.setType(Constants.SHARED_SUBSTRING,false);
						edge.setType(Constants.OVERLAP_PREFIX_PREFIX,true);
						if (edge.overhangs==null) edge.overhangs = new int[4][2];
						Math.set(edge.overhangs,-1);
						edge.overhangs[0][0]=edge.sharedSubstringOverhangs[1];
						edge.overhangs[0][1]=edge.sharedSubstringOverhangs[3];
					}
					else if (edge.sharedSubstringOverhangs[3]<=IDENTITY_THRESHOLD) {
						edge.setType(Constants.SHARED_SUBSTRING,false);
						edge.setType(Constants.OVERLAP_PREFIX_SUFFIX,true);
						if (edge.overhangs==null) edge.overhangs = new int[4][2];
						Math.set(edge.overhangs,-1);
						edge.overhangs[1][0]=edge.sharedSubstringOverhangs[1];
						edge.overhangs[1][1]=edge.sharedSubstringOverhangs[2];
					}
				}
				else if (edge.sharedSubstringOverhangs[1]==0) {
					if (edge.sharedSubstringOverhangs[2]<=IDENTITY_THRESHOLD) {
						edge.setType(Constants.SHARED_SUBSTRING,false);
						edge.setType(Constants.OVERLAP_SUFFIX_PREFIX,true);
						if (edge.overhangs==null) edge.overhangs = new int[4][2];
						Math.set(edge.overhangs,-1);
						edge.overhangs[2][0]=edge.sharedSubstringOverhangs[0];
						edge.overhangs[2][1]=edge.sharedSubstringOverhangs[3];
					}
					else if (edge.sharedSubstringOverhangs[3]<=IDENTITY_THRESHOLD) {
						edge.setType(Constants.SHARED_SUBSTRING,false);
						edge.setType(Constants.OVERLAP_SUFFIX_SUFFIX,true);
						if (edge.overhangs==null) edge.overhangs = new int[4][2];
						Math.set(edge.overhangs,-1);
						edge.overhangs[3][0]=edge.sharedSubstringOverhangs[0];
						edge.overhangs[3][1]=edge.sharedSubstringOverhangs[2];
					}
				}
			}
		}
		else {
			if (edge.isType(Constants.CONTAINMENT_ONE_IN_TWO) || edge.isType(Constants.INSERTION_ONE_IN_TWO)) {
				if (newNodeStart>oldNodeStart) {
					if (newNodeStart<=oldNodeStart+edge.containmentOverhangLeft) edge.containmentOverhangLeft-=newNodeStart-oldNodeStart;
					else edge.containmentOverhangLeft=0;
				}
				else if (newNodeStart<oldNodeStart) edge.containmentOverhangLeft+=oldNodeStart-newNodeStart;
				if (newNodeEnd>oldNodeEnd) edge.containmentOverhangRight+=newNodeEnd-oldNodeEnd;
				else if (newNodeEnd<oldNodeEnd) {
					if (newNodeEnd>=oldNodeEnd-edge.containmentOverhangRight) edge.containmentOverhangRight-=oldNodeEnd-newNodeEnd;
					else edge.containmentOverhangRight=0;
				}
				if (edge.containmentOverhangLeft==0 && edge.containmentOverhangRight==0) {
					edge.setType_noContainment();
					edge.setType_noInsertion();
					edge.setType(Constants.CONTAINMENT_IDENTICAL,true);
				}
				else if (edge.containmentOverhangLeft==0) {
					if ( edge.isType(Constants.INSERTION_ONE_IN_TWO) &&
				         ( (edge.orientation==0 && !nodesArray[otherNode].isRightMaximal) ||
						   (edge.orientation==1 && !nodesArray[otherNode].isLeftMaximal)
						 )
				       ) {
					   edge.setType_noInsertion();
					   edge.setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
					}
					edge.containmentSubtype=Constants.CONTAINMENT_PREFIX;
				}
				else if (edge.containmentOverhangRight==0) {
					if ( edge.isType(Constants.INSERTION_ONE_IN_TWO) &&
				         ( (edge.orientation==0 && !nodesArray[otherNode].isLeftMaximal) ||
						   (edge.orientation==1 && !nodesArray[otherNode].isRightMaximal)
						 )
				       ) {
					   edge.setType_noInsertion();
					   edge.setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
					}
					edge.containmentSubtype=Constants.CONTAINMENT_SUFFIX;
				}
			}
			else if (edge.isType_overlap() && edge.overhangs!=null) {
				j=0;
				for (i=Constants.OVERLAP_PREFIX_PREFIX; i<=Constants.OVERLAP_SUFFIX_SUFFIX; i<<=1) {
					if (!edge.isType(i)) {
						j++;
						continue;
					}
					if (i==Constants.OVERLAP_PREFIX_PREFIX || i==Constants.OVERLAP_SUFFIX_PREFIX) {
						if (newNodeEnd>oldNodeEnd) edge.overhangs[j][1]+=newNodeEnd-oldNodeEnd;
						else if (newNodeEnd<oldNodeEnd) {
							if (newNodeEnd>=oldNodeEnd-edge.overhangs[j][1]) edge.overhangs[j][1]-=oldNodeEnd-newNodeEnd;
							else edge.overhangs[j][1]=0;
						}
					}
					else if (i==Constants.OVERLAP_PREFIX_SUFFIX || i==Constants.OVERLAP_SUFFIX_SUFFIX) {
						if (newNodeStart>oldNodeStart) {
							if (newNodeStart<=oldNodeStart+edge.overhangs[j][1]) edge.overhangs[j][1]-=newNodeStart-oldNodeStart;
							else edge.overhangs[j][1]=0;
						}
						else if (newNodeStart<oldNodeStart) edge.overhangs[j][1]+=oldNodeStart-newNodeStart;
					}
					if (edge.overhangs[j][1]==0) {
						edge.setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
						if (i==Constants.OVERLAP_PREFIX_PREFIX || i==Constants.OVERLAP_PREFIX_SUFFIX) {
						   	edge.containmentSubtype=Constants.CONTAINMENT_PREFIX;
						   	edge.containmentOverhangLeft=0;
						   	edge.containmentOverhangRight=edge.overhangs[j][0];
						}
						else if (i==Constants.OVERLAP_SUFFIX_PREFIX || i==Constants.OVERLAP_SUFFIX_SUFFIX) {
							edge.containmentSubtype=Constants.CONTAINMENT_SUFFIX;
							edge.containmentOverhangLeft=edge.overhangs[j][0];
							edge.containmentOverhangRight=0;
						}
						edge.setType_noOverlap();
						break;
					}
					j++;
				}
			}
			else if (edge.isType(Constants.SHARED_SUBSTRING) && edge.sharedSubstringOverhangs!=null) {
				if (newNodeStart>oldNodeStart) {
					if (newNodeStart<=oldNodeStart+edge.sharedSubstringOverhangs[2]) edge.sharedSubstringOverhangs[2]-=newNodeStart-oldNodeStart;
					else edge.sharedSubstringOverhangs[2]=0;
				}
				else if (newNodeStart<oldNodeStart) edge.sharedSubstringOverhangs[2]+=oldNodeStart-newNodeStart;
				if (newNodeEnd>oldNodeEnd) edge.sharedSubstringOverhangs[3]+=newNodeEnd-oldNodeEnd;
				else if (newNodeEnd<oldNodeEnd) {
					if (newNodeEnd>=oldNodeEnd-edge.sharedSubstringOverhangs[3]) edge.sharedSubstringOverhangs[3]-=oldNodeEnd-newNodeEnd;
					else edge.sharedSubstringOverhangs[3]=0;
				}
				if (edge.sharedSubstringOverhangs[2]==0 && edge.sharedSubstringOverhangs[3]==0) {
					edge.setType(Constants.SHARED_SUBSTRING,false);
					edge.setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
					edge.containmentOverhangLeft=edge.sharedSubstringOverhangs[0];
					edge.containmentOverhangRight=edge.sharedSubstringOverhangs[1];
					if (edge.containmentOverhangLeft<=IDENTITY_THRESHOLD) edge.containmentSubtype=Constants.CONTAINMENT_PREFIX;
					else if (edge.containmentOverhangRight<=IDENTITY_THRESHOLD) edge.containmentSubtype=Constants.CONTAINMENT_SUFFIX;
					else edge.containmentSubtype=Constants.CONTAINMENT_SUBSTRING;
				}
				else if (edge.sharedSubstringOverhangs[2]==0) {
					if (edge.sharedSubstringOverhangs[0]<=IDENTITY_THRESHOLD) {
						edge.setType(Constants.SHARED_SUBSTRING,false);
						edge.setType(Constants.OVERLAP_PREFIX_PREFIX,true);
						if (edge.overhangs==null) edge.overhangs = new int[4][2];
						Math.set(edge.overhangs,-1);
						edge.overhangs[0][0]=edge.sharedSubstringOverhangs[1];
						edge.overhangs[0][1]=edge.sharedSubstringOverhangs[3];
					}
					else if (edge.sharedSubstringOverhangs[1]<=IDENTITY_THRESHOLD) {
						edge.setType(Constants.SHARED_SUBSTRING,false);
						edge.setType(Constants.OVERLAP_SUFFIX_PREFIX,true);
						if (edge.overhangs==null) edge.overhangs = new int[4][2];
						Math.set(edge.overhangs,-1);
						edge.overhangs[2][0]=edge.sharedSubstringOverhangs[0];
						edge.overhangs[2][1]=edge.sharedSubstringOverhangs[3];
					}
				}
				else if (edge.sharedSubstringOverhangs[3]==0) {
					if (edge.sharedSubstringOverhangs[0]<=IDENTITY_THRESHOLD) {
						edge.setType(Constants.SHARED_SUBSTRING,false);
						edge.setType(Constants.OVERLAP_PREFIX_SUFFIX,true);
						if (edge.overhangs==null) edge.overhangs = new int[4][2];
						Math.set(edge.overhangs,-1);
						edge.overhangs[1][0]=edge.sharedSubstringOverhangs[1];
						edge.overhangs[1][1]=edge.sharedSubstringOverhangs[2];
					}
					else if (edge.sharedSubstringOverhangs[1]<=IDENTITY_THRESHOLD) {
						edge.setType(Constants.SHARED_SUBSTRING,false);
						edge.setType(Constants.OVERLAP_SUFFIX_SUFFIX,true);
						if (edge.overhangs==null) edge.overhangs = new int[4][2];
						Math.set(edge.overhangs,-1);
						edge.overhangs[3][0]=edge.sharedSubstringOverhangs[0];
						edge.overhangs[3][1]=edge.sharedSubstringOverhangs[2];
					}
				}
			}
		}
	}
	
	
	/**
	 * Removes all nodes with no incident edge.
	 *
	 * Remark: the procedure prints kernel tags files as specified by $IntervalGraphStep3.
	 * printKernelTags()$ for all nodes it removes; every such node is printed as if it 
	 * had no kernel tag.
	 *
	 * Remark: the procedure does not assume $nodesArray$ to be sorted, and it does not 
	 * create new intervals.
	 *
	 * @return the number of removed nodes.
	 */
	private static final int removeDisconnectedNodes(String tagsDir, String tagsPrefix) throws IOException {
		final int ARTIFICIAL_KERNEL = Math.POSITIVE_INFINITY;
		int i, j;
		int out, nNodesPrime;
		Node node;
		
		// Assigning new IDs to connected nodes
		nNodesPrime=0;
		for (i=0; i<nNodes; i++) {
			if (nNeighbors[i]>0) nNodesPrime++;
		}
		if (nNodesPrime==nNodes) return 0;
		out=nNodes-nNodesPrime;
		j=-1;
		for (i=0; i<nNodes; i++) {
			if (nNeighbors[i]==0) continue;
			nodesArray[i].nodeID=++j;
		}
		
		// Updating edges
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=-1;
		}
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].printed==1) continue;
				neighbors[i][j].nodeID1=nodesArray[neighbors[i][j].nodeID1].nodeID;
				neighbors[i][j].nodeID2=nodesArray[neighbors[i][j].nodeID2].nodeID;
				neighbors[i][j].printed=1;
			}
		}
		
		// Printing kernel tags file
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if (nNeighbors[i]>0) {
				node.lastKernel=-1;
				continue;
			}
			node.lastKernel=0;
			if (node.kernels==null || node.kernels.length==0) node.kernels = new int[1];
			node.kernels[0]=ARTIFICIAL_KERNEL;
			if (node.kernelOrientations==null || node.kernelOrientations.length==0) node.kernelOrientations = new byte[1];
			node.kernelOrientations[0]=(byte)0;
		}
		IntervalGraphStep3.pathKernelLengths=null;
		IntervalGraphStep3.printKernelTags(tagsDir+"/"+IO.TAGS_ROOT_FILE,tagsPrefix,true,true,ARTIFICIAL_KERNEL);
		
		// Compacting $nodesArray$ and $neighbors$.
		for (i=0; i<nNodes; i++) {
			if (nNeighbors[i]==0) continue;
			j=nodesArray[i].nodeID;
			nodesArray[j]=nodesArray[i];
			nNeighbors[j]=nNeighbors[i];
			neighbors[j]=neighbors[i];
		}
		for (i=nNodesPrime; i<nNodes; i++) {
			nodesArray[i]=null; neighbors[i]=null;
			nNeighbors[i]=0;
		}
		nNodes=nNodesPrime;
		return out;
	}
	
	
	
	
	
	
	
	
	// ------------------------------ FIXING OVERLAP FORKS -------------------------------

	/**
	 * Tries to remove "overlap forks", i.e. nodes V with overlap edges, from the same end
	 * of V, to multiple neighbors that do not have overlap or containment edges with
	 * one another. The procedure tries to remove wrong forks due to undersplitting the 
	 * neighbors during factorization, i.e. by the fact that every neighbors contains a 
	 * spurious end. See procedure $isSpuriousOverlapFork$ for details. Forks that do not 
	 * conform to specific undersplit patterns are kept, and are assumed to represent 
	 * shared substrings among different repeats.
	 *
	 * Remark: we could just avoid removing wrong forks and stop a contig at a fork during
	 * kernel assembly. However, this would split the sequence of a repeat into contigs.
	 *
	 * Remark: the procedure does not necessarily completely remove a fork, but it can 
	 * reduce the number of its overlap edges.
	 *
	 * Remark: the procedure could be seen as a way to detect undersplits that could not
	 * be detected by $fixUndersplits()$, by exploiting specific signals. The problematic 
	 * neighbors might not be undersplits that could be resolved by $fixUndersplits()$.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start$, and it
	 * might create new intervals.
	 */
	private static final void fixOverlapForks() {
		final int MIN_OVERLAP_NEIGHBORS = 100;
		boolean oneFixNeighbor;
		int i, j, k;
		int to, idGenerator, previousIdGenerator, isFork, lastSuffixOverlapNeighbor, lastOverlapNeighbor;
		int nForks, nForksSuffix, nForksPrefix, nForksSuffixPrefix, nFixedNeighbors;
		int originalNNeighbors;
		Node bigIntervalTmp = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Node smallIntervalTmp = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		AlignmentPair tmpPair = new AlignmentPair();
		Edge tmpEdge = new Edge();
		int[] splittingPosition, positions;

		idGenerator=nNodes;
		splittingPosition = new int[MIN_OVERLAP_NEIGHBORS];
		positions = new int[(MIN_OVERLAP_NEIGHBORS)<<1];
		if (tmpBoolean.length<MIN_OVERLAP_NEIGHBORS+1) tmpBoolean = new boolean[MIN_OVERLAP_NEIGHBORS+1];
		nForks=0; nForksSuffix=0; nForksPrefix=0; nForksSuffixPrefix=0; nFixedNeighbors=0;
		for (i=0; i<nNodes; i++) {
			Edge.order=Edge.UNSORTED-1-i;
			if (nNeighbors[i]>1) Arrays.sort(neighbors[i],0,nNeighbors[i]);
			originalNNeighbors=nNeighbors[i];
			overlapEdgesBoundaries(i,tmpArray);
			lastSuffixOverlapNeighbor=tmpArray[0];
			lastOverlapNeighbor=tmpArray[1];
			if (splittingPosition.length<nNeighbors[i]) splittingPosition = new int[nNeighbors[i]];
			if (positions.length<nNeighbors[i]<<1) positions = new int[nNeighbors[i]<<1];
			isFork=isSpuriousOverlapFork(i,lastSuffixOverlapNeighbor,lastOverlapNeighbor,splittingPosition,positions);
			if (isFork==0) continue;
			else nForks++;
			oneFixNeighbor=false;			
			if ((isFork&1)!=0) {  // Fork on the suffix side
				nForksSuffix++;
				for (j=0; j<=lastSuffixOverlapNeighbor; j++) {
					if (splittingPosition[j]<=0) continue;
					nFixedNeighbors++;
					to=neighbors[i][j].getTo(i);
					idGenerator=fixNeighbor(to,i,neighbors[i][j].isPrefixOverlap(to),splittingPosition[j],idGenerator,bigIntervalTmp,smallIntervalTmp,tmpPair,tmpEdge,positions);
					oneFixNeighbor=true;
				}
			}
			if ((isFork&2)!=0) {  // Fork on the prefix side
				nForksPrefix++;
				for (j=lastSuffixOverlapNeighbor+1; j<=lastOverlapNeighbor; j++) {
					if (splittingPosition[j]<=0) continue;
					nFixedNeighbors++;
					to=neighbors[i][j].getTo(i);
					idGenerator=fixNeighbor(to,i,neighbors[i][j].isPrefixOverlap(to),splittingPosition[j],idGenerator,bigIntervalTmp,smallIntervalTmp,tmpPair,tmpEdge,positions);
					oneFixNeighbor=true;
				}
			}
			if ((isFork&1)!=0 && (isFork&2)!=0) nForksSuffixPrefix++;
			if (oneFixNeighbor) {
				// Removing null neighbors from the fork
				k=-1;
				for (j=0; j<nNeighbors[i]; j++) {
					if (neighbors[i][j]==null) continue;
					k++;
					if (k==j) continue;
					neighbors[i][k]=neighbors[i][j];
				}
				nNeighbors[i]=k+1;
			}
		}
		System.err.println("fixOverlapForks> Detected "+nForks+" forks ("+nForksSuffix+" suffix, "+nForksPrefix+" prefix, "+nForksSuffixPrefix+" suffix-prefix)");
		System.err.println("fixOverlapForks> Fixed neighbors: "+nFixedNeighbors+" ("+IO.getPercent(nFixedNeighbors,nNodes)+"%)");
		System.err.println("fixOverlapForks> New nodes created: "+(idGenerator-nNodes)+" ("+IO.getPercent(idGenerator-nNodes,nNodes)+"%)");
		nNodes=idGenerator;
	}


	/**
	 * Splits the neighbor $nodesArray[n]$ of fork $forkID$ into three intervals: a prefix 
	 * and a suffix, separated at position $splittingPosition$, and a new interval that
	 * coincides with the old interval. Edges adjacent to node $n$ are redirected to the
	 * new nodes their alignments fall into. The procedure tries to reuse an existing
	 * interval for the interval that does not correspond to a prefix or suffix substring,
	 * and it creates a new node only if a similar one is not found.
	 *
	 * Remark: if a new interval does not collect alignments, it is not created. In 
	 * particular, if the interval that does not correspond to a prefix or suffix 
	 * substring is too small, it is automatically not created. The information that such
	 * short segment is repetitive is preserved in the new interval that coincides with 
	 * the old interval.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read$, and the
	 * $avgDiffs$ fields of edges not to have been divided by the number of alignments per
	 * edge.
	 *
	 * @param n $nodesArray[n]$ is assumed not to be periodic;
	 * @param prefixOrSuffix $nodesArray[n]$ is considered as a prefix (TRUE) or a suffix
	 * (FALSE) interval;
	 * @param splittingPosition relative to $nodesArray[n].start$;
	 * @param bigIntervalTmp,smallIntervalTmp,tmpPair,tmpEdge temporary space;
	 * @param idGenerator new nodes get IDs starting from this value;
	 * @param tmpNullNeighbors temporary space, of size equal at least to the number of 
	 * neighbors of $n$;
	 * @return the value of $idGenerator$ after the procedure completes.
	 */
	private static final int fixNeighbor(int n, int forkID, boolean prefixOrSuffix, int splittingPosition, int idGenerator, Node bigIntervalTmp, Node smallIntervalTmp, AlignmentPair tmpPair, Edge tmpEdge, int[] tmpNullNeighbors) {
		int i, j, k, p;
		int read, smallStart, smallEnd, smallComplementStart, smallComplementEnd;
		int nTransformedAlignments, bigIntervalID, smallIntervalID, neighbor;
		int projectionStart, projectionEnd, projectionStartNeighbor, projectionEndNeighbor;
		int lastNullNeighbor, originalNNeighbors;
		Node oldInterval, bigInterval, smallInterval, newNode, tmpNode;
		Edge edge, newEdge;

		// Building big and small intervals and adding them to $nodesArray$
		oldInterval=nodesArray[n]; read=oldInterval.read;		
		if (prefixOrSuffix) {
			smallComplementStart=oldInterval.start;
			smallComplementEnd=oldInterval.start+splittingPosition-1;
			smallStart=oldInterval.start+splittingPosition;
			smallEnd=oldInterval.end;
		}
		else {
			smallStart=oldInterval.start;
			smallEnd=oldInterval.start+splittingPosition;
			smallComplementStart=oldInterval.start+splittingPosition+1;
			smallComplementEnd=oldInterval.end;
		}		
		bigIntervalTmp.clear();
		bigIntervalTmp.read=read;
		bigIntervalTmp.type=Constants.INTERVAL_ALIGNMENT;
		bigIntervalTmp.start=oldInterval.start;
		bigIntervalTmp.end=oldInterval.end;
		bigIntervalTmp.isLeftMaximal=oldInterval.isLeftMaximal;
		bigIntervalTmp.isRightMaximal=oldInterval.isRightMaximal;
		bigIntervalTmp.isContainedInterval=oldInterval.isContainedInterval;
		bigInterval=bigIntervalTmp;
		expandTo(idGenerator);
		nodesArray[idGenerator]=bigInterval; bigInterval.nodeID=idGenerator;		
		if (neighbors[idGenerator]==null || neighbors[idGenerator].length<nNeighbors[n]) neighbors[idGenerator] = new Edge[nNeighbors[n]];
		nNeighbors[idGenerator]=0;
		idGenerator++;
		smallInterval=null;
		for (i=n-1; i>=0; i--) {
			if (nodesArray[i].read!=read) break;
			if (Intervals.areApproximatelyIdentical(nodesArray[i].start,nodesArray[i].end,smallStart,smallEnd)) {
				smallInterval=nodesArray[i];
				break;
			}
		}
		if (smallInterval==null) {
			for (i=n+1; i<nNodes; i++) {
				if (nodesArray[i].read!=read) break;
				if (Intervals.areApproximatelyIdentical(nodesArray[i].start,nodesArray[i].end,smallStart,smallEnd)) {
					smallInterval=nodesArray[i];
					break;
				}
			}
		}
		if (smallInterval==null) {
			smallIntervalTmp.clear();
			smallIntervalTmp.read=read;
			smallIntervalTmp.type=Constants.INTERVAL_ALIGNMENT;
			smallIntervalTmp.start=smallStart;
			smallIntervalTmp.end=smallEnd;
			smallIntervalTmp.isLeftMaximal=prefixOrSuffix?true:oldInterval.isLeftMaximal;
			smallIntervalTmp.isRightMaximal=prefixOrSuffix?oldInterval.isRightMaximal:true;
			smallIntervalTmp.isContainedInterval=oldInterval.isContainedInterval;
			smallInterval=smallIntervalTmp;
			expandTo(idGenerator);
			nodesArray[idGenerator]=smallInterval; smallInterval.nodeID=idGenerator;
			if (neighbors[idGenerator]==null || neighbors[idGenerator].length<nNeighbors[n]) neighbors[idGenerator] = new Edge[nNeighbors[n]];
			nNeighbors[idGenerator]=0;
			idGenerator++;
		}

		// Transforming edges
		lastNullNeighbor=-1; originalNNeighbors=nNeighbors[n];
		newNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		oldInterval.clone(newNode); newNode.nodeID=oldInterval.nodeID;
		newNode.start=smallComplementStart; newNode.end=smallComplementEnd;
		if (prefixOrSuffix) newNode.isRightMaximal=true;
		else newNode.isLeftMaximal=true;
		nTransformedAlignments=0;
		for (i=0; i<nNeighbors[n]; i++) {
			edge=neighbors[n][i];
			if (edge.orientation==2 || edge.orientation==-1) {
				// Cannot decompress the alignments inside the edge
				edge.changeType_alignment(n);
				if (edge.nodeID1==n) edge.nodeID1=bigInterval.nodeID;
				else edge.nodeID2=bigInterval.nodeID;
				neighbors[n][i]=null;
				j=findEdge(bigInterval.nodeID,edge,false);
				if (j==Math.NEGATIVE_INFINITY) addEdge(bigInterval.nodeID,edge);
				else {
					neighbors[bigInterval.nodeID][j].addTypes(edge);
					neighbors[bigInterval.nodeID][j].addDiffs(edge);
					neighbors[bigInterval.nodeID][j].addOrientation(edge);
				}
				continue;
			}
			neighbor=edge.getTo(n);
			neighbors[n][i]=null;
			neighbors[neighbor][findEdge(neighbor,edge,false)]=null;
			tmpNullNeighbors[++lastNullNeighbor]=neighbor;
			if (edge.containment!=-1) {
				p=edge.projectOnto_containment(n,tmpArray,-1);
				if (p>=0) {
					projectionStart=tmpArray[0]; 
					projectionEnd=tmpArray[1];
				}
				else {
					projectionStart=0;
					projectionEnd=nodesArray[n].length()-1;
				}
				p=edge.projectOnto_containment(neighbor,tmpArray,-1);
				if (p>=0) {
					projectionStartNeighbor=tmpArray[0]; 
					projectionEndNeighbor=tmpArray[1];
				}
				else {
					projectionStartNeighbor=0;
					projectionEndNeighbor=nodesArray[neighbor].length()-1;
				}
				fixNeighbor_impl(n,nodesArray[neighbor],projectionStart,projectionEnd,projectionStartNeighbor,projectionEndNeighbor,edge.orientation,edge.avgDiffs/edge.getNAlignments(),edge.implied,smallInterval,smallComplementStart,smallComplementEnd,bigInterval,newNode,tmpPair,tmpEdge);
				nTransformedAlignments++;
			}
			if (edge.insertion!=-1) {
				p=edge.projectOnto_insertion(n,tmpArray,-1);
				if (p>=0) {
					projectionStart=tmpArray[0]; 
					projectionEnd=tmpArray[1];
				}
				else {
					projectionStart=0;
					projectionEnd=nodesArray[n].length()-1;
				}
				p=edge.projectOnto_insertion(neighbor,tmpArray,-1);
				if (p>=0) {
					projectionStartNeighbor=tmpArray[0]; 
					projectionEndNeighbor=tmpArray[1];
				}
				else {
					projectionStartNeighbor=0;
					projectionEndNeighbor=nodesArray[neighbor].length()-1;
				}
				fixNeighbor_impl(n,nodesArray[neighbor],projectionStart,projectionEnd,projectionStartNeighbor,projectionEndNeighbor,edge.orientation,edge.avgDiffs/edge.getNAlignments(),edge.implied,smallInterval,smallComplementStart,smallComplementEnd,bigInterval,newNode,tmpPair,tmpEdge);
				nTransformedAlignments++;
			}
			if (edge.overlap!=-1 && edge.overhangs!=null) {
				for (j=Constants.OVERLAP_PREFIX_PREFIX; j<=Constants.OVERLAP_SUFFIX_SUFFIX; j<<=1) {
					p=edge.projectOnto_overlap(n,tmpArray,j,-1);
					if (p<0) {
						// $edge.overhangs=null$
						edge.changeType_alignment(n);
						if (edge.nodeID1==n) edge.nodeID1=bigInterval.nodeID;
						else edge.nodeID2=bigInterval.nodeID;
						k=findEdge(bigInterval.nodeID,edge,false);
						if (k==Math.NEGATIVE_INFINITY) addEdge(bigInterval.nodeID,edge);
						else {
							neighbors[bigInterval.nodeID][k].addTypes(edge);
							neighbors[bigInterval.nodeID][k].addDiffs(edge);
							neighbors[bigInterval.nodeID][k].addOrientation(edge);
						}
						continue;
					}
					projectionStart=tmpArray[0]; projectionEnd=tmpArray[1];
					p=edge.projectOnto_overlap(neighbor,tmpArray,j,-1);
					if (p<0) {
						// $edge.overhangs=null$
						edge.changeType_alignment(n);
						if (edge.nodeID1==n) edge.nodeID1=bigInterval.nodeID;
						else edge.nodeID2=bigInterval.nodeID;
						k=findEdge(bigInterval.nodeID,edge,false);
						if (k==Math.NEGATIVE_INFINITY) addEdge(bigInterval.nodeID,edge);
						else {
							neighbors[bigInterval.nodeID][k].addTypes(edge);
							neighbors[bigInterval.nodeID][k].addDiffs(edge);
							neighbors[bigInterval.nodeID][k].addOrientation(edge);
						}
						continue;
					}
					projectionStartNeighbor=tmpArray[0]; projectionEndNeighbor=tmpArray[1];
					fixNeighbor_impl(n,nodesArray[neighbor],projectionStart,projectionEnd,projectionStartNeighbor,projectionEndNeighbor,edge.orientation,edge.avgDiffs/edge.getNAlignments(),edge.implied,smallInterval,smallComplementStart,smallComplementEnd,bigInterval,newNode,tmpPair,tmpEdge);
					nTransformedAlignments++;
				}
			}
			if (edge.sharedSubstring!=-1 && edge.sharedSubstringOverhangs!=null) {
				p=edge.projectOnto_sharedSubstring(n,tmpArray,-1);
				if (p<0) {
					// $edge.sharedSubstringOverhangs=null$
					edge.changeType_alignment(n);
					if (edge.nodeID1==n) edge.nodeID1=bigInterval.nodeID;
					else edge.nodeID2=bigInterval.nodeID;
					j=findEdge(bigInterval.nodeID,edge,false);
					if (j==Math.NEGATIVE_INFINITY) addEdge(bigInterval.nodeID,edge);
					else {
						neighbors[bigInterval.nodeID][j].addTypes(edge);
						neighbors[bigInterval.nodeID][j].addDiffs(edge);
						neighbors[bigInterval.nodeID][j].addOrientation(edge);
					}
					continue;
				}
				projectionStart=tmpArray[0]; projectionEnd=tmpArray[1];
				p=edge.projectOnto_sharedSubstring(neighbor,tmpArray,-1);
				if (p<0) {
					// $edge.sharedSubstringOverhangs=null$
					edge.changeType_alignment(n);
					if (edge.nodeID1==n) edge.nodeID1=bigInterval.nodeID;
					else edge.nodeID2=bigInterval.nodeID;
					j=findEdge(bigInterval.nodeID,edge,false);
					if (j==Math.NEGATIVE_INFINITY) addEdge(bigInterval.nodeID,edge);
					else {
						neighbors[bigInterval.nodeID][j].addTypes(edge);
						neighbors[bigInterval.nodeID][j].addDiffs(edge);
						neighbors[bigInterval.nodeID][j].addOrientation(edge);
					}
					continue;
				}
				projectionStartNeighbor=tmpArray[0]; projectionEndNeighbor=tmpArray[1];
				fixNeighbor_impl(n,nodesArray[neighbor],projectionStart,projectionEnd,projectionStartNeighbor,projectionEndNeighbor,edge.orientation,edge.avgDiffs/edge.getNAlignments(),edge.implied,smallInterval,smallComplementStart,smallComplementEnd,bigInterval,newNode,tmpPair,tmpEdge);			
				nTransformedAlignments++;
			}
			if (nTransformedAlignments<edge.getNAlignments()) {
				// Remark: the original edge might have contained just multiple overlap
				// edges, or just multiple insertion edges, only one of which can be
				// decompressed. In this case we should set the type of the edge after 
			 	// decompression to overlap or insertion, but we don't have the
				// information to do so at this point, so we set it to the generic shared
				// substring type. This problem could be solved by storing the number of
				// alignments of each type inside an edge, at the cost of increasing 
				// space.
				tmpEdge.clear();
				if (read<nodesArray[neighbor].read) {
					tmpEdge.nodeID1=bigInterval.nodeID;
					tmpEdge.nodeID2=neighbor;
				}
				else {
					tmpEdge.nodeID2=bigInterval.nodeID;
					tmpEdge.nodeID1=neighbor;
				}
				tmpEdge.sharedSubstring=Constants.SHARED_SUBSTRING;
				tmpEdge.avgDiffs=(edge.getNAlignments()-nTransformedAlignments)*edge.avgDiffs/edge.getNAlignments();
				tmpEdge.nAlignmentsForward=edge.nAlignmentsForward==0?0:edge.nAlignmentsForward-nTransformedAlignments;
				tmpEdge.nAlignmentsBackward=edge.nAlignmentsBackward==0?0:edge.nAlignmentsBackward-nTransformedAlignments;
				tmpEdge.implied=edge.implied;
				tmpEdge.orientation=edge.orientation;
				j=findEdge(bigInterval.nodeID,tmpEdge,false);
				if (j==Math.NEGATIVE_INFINITY) {
					newEdge = new Edge(); tmpEdge.clone(newEdge);
					addEdge(bigInterval.nodeID,newEdge);
					addEdge(neighbor,newEdge);
				}
				else {
					neighbors[bigInterval.nodeID][j].addTypes(tmpEdge);
					neighbors[bigInterval.nodeID][j].addDiffs(tmpEdge);
					neighbors[bigInterval.nodeID][j].addOrientation(tmpEdge);
				}
			}
		}		
		
		// Removing null edges from $neighbors$
		j=-1;
		for (i=0; i<nNeighbors[n]; i++) {
			if (neighbors[n][i]==null) continue;
			j++;
			if (j==i) continue;
			neighbors[n][j]=neighbors[n][i];
		}
		nNeighbors[n]=j+1;
		for (i=0; i<=lastNullNeighbor; i++) {
			neighbor=tmpNullNeighbors[i];
			if (neighbor==forkID) continue;			
			k=-1;
			for (j=0; j<nNeighbors[neighbor]; j++) {
				if (neighbors[neighbor][j]==null) continue;
				k++;
				if (k==j) continue;
				neighbors[neighbor][k]=neighbors[neighbor][j];
			}
			nNeighbors[neighbor]=k+1;
		}
		if (IO.CONSISTENCY_CHECKS && nNeighbors[n]>originalNNeighbors) {
			System.err.println("fixNeighbor> ERROR: nNeighbors: before="+originalNNeighbors+" after="+nNeighbors[n]);
			System.exit(1);
		}

		// Removing new nodes without neighbors
		nodesArray[n]=newNode;  // Kept even if it has no neighbor
		bigIntervalID=bigInterval.nodeID; smallIntervalID=smallInterval.nodeID;
		if (nNeighbors[bigIntervalID]==0) {
			nodesArray[bigIntervalID]=null; neighbors[bigIntervalID]=null;
			if (smallIntervalID!=bigIntervalID+1) return bigIntervalID;
			if (nNeighbors[smallIntervalID]==0) {
				nodesArray[smallIntervalID]=null; neighbors[smallIntervalID]=null;
				return bigIntervalID;
			}
			nodesArray[bigIntervalID] = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
			nodesArray[smallIntervalID].clone(nodesArray[bigIntervalID]);
			nodesArray[bigIntervalID].nodeID=bigIntervalID;
			nodesArray[smallIntervalID]=null;
			if (neighbors[bigIntervalID]==null || neighbors[bigIntervalID].length<neighbors[smallIntervalID].length) neighbors[bigIntervalID] = new Edge[neighbors[smallIntervalID].length];
			System.arraycopy(neighbors[smallIntervalID],0,neighbors[bigIntervalID],0,nNeighbors[smallIntervalID]);
			nNeighbors[bigIntervalID]=nNeighbors[smallIntervalID];
			neighbors[smallIntervalID]=null; nNeighbors[smallIntervalID]=0;
			for (i=0; i<nNeighbors[bigIntervalID]; i++) {
				if (neighbors[bigIntervalID][i].nodeID1==smallIntervalID) neighbors[bigIntervalID][i].nodeID1=bigIntervalID;
				else neighbors[bigIntervalID][i].nodeID2=bigIntervalID;
			}
			return smallIntervalID;
		}
		else {
			tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
			nodesArray[bigIntervalID].clone(tmpNode);
			tmpNode.nodeID=bigIntervalID;
			nodesArray[bigIntervalID]=tmpNode;
			if (smallIntervalID!=bigIntervalID+1) return bigIntervalID+1;
			if (nNeighbors[smallIntervalID]==0) {
				nodesArray[smallIntervalID]=null; neighbors[smallIntervalID]=null;
				return smallIntervalID;
			}
			tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
			nodesArray[smallIntervalID].clone(tmpNode);
			tmpNode.nodeID=smallIntervalID;
			nodesArray[smallIntervalID]=tmpNode;
			return smallIntervalID+1;
		}		
	}


	/**
	 * Remark: the procedure assumes that no edge between $nodeID$ and $neighbor$ exists
	 * in $neighbors$, i.e. that the original edge we want to redirect has been removed
	 * from $neighbors$.
	 *
	 * Remark: the procedure updates the global variable $Edge.idGenerator$.
	 *
	 * @param nodeID ID of the neighbor of the fork that has to be fixed;
	 * @param neighbor neighbor of $nodeID$;
	 * @param projectionStart,projectionEnd projection of the alignment on $nodeID$
	 * (relative to $nodesArray[nodeID].start$);
	 * @param projectionStartNeighbor,projectionEndNeighbor projection of the alignment on
	 * $neighbor$ (relative to $neighbor.start$);
	 * @param smallInterval node that corresponds to the small interval in the splitting
	 * of $nodeID$;
	 * @param smallStartComplement,smallEndComplement absolute;
	 * @param bigInterval node that corresponds to the big interval in the splitting
	 * of $nodeID$;
	 * @param newNode a version of $nodesArray[nodeID]$ with updated boundaries, but with
	 * the same $nodeID$ as $nodesArray[nodeID]$;
	 * @param tmpPair,tmpEdge temporary space.
	 */
	private static final void fixNeighbor_impl(int nodeID, Node neighbor, int projectionStart, int projectionEnd, int projectionStartNeighbor, int projectionEndNeighbor, int orientation, double diffs, boolean implied, Node smallInterval, int smallStartComplement, int smallEndComplement, Node bigInterval, Node newNode, AlignmentPair tmpPair, Edge tmpEdge) {
		int i;
		Edge newEdge = null;

		tmpPair.orientation=orientation==0;
		tmpPair.diffs=diffs;
		tmpPair.implied=implied;
		if (nodesArray[nodeID].read<neighbor.read) {
			tmpPair.node2=neighbor;
			tmpPair.alignmentStart2=neighbor.start+projectionStartNeighbor;
			tmpPair.alignmentEnd2=neighbor.start+projectionEndNeighbor;
		}
		else {
			tmpPair.node1=neighbor;
			tmpPair.alignmentStart1=neighbor.start+projectionStartNeighbor;
			tmpPair.alignmentEnd1=neighbor.start+projectionEndNeighbor;
		}
		if ( Intervals.isApproximatelyContained(nodesArray[nodeID].start+projectionStart,nodesArray[nodeID].start+projectionEnd,smallInterval.start,smallInterval.end) ||
			 Intervals.areApproximatelyIdentical(nodesArray[nodeID].start+projectionStart,nodesArray[nodeID].start+projectionEnd,smallInterval.start,smallInterval.end)
		   ) {
			if (nodesArray[nodeID].read<neighbor.read) {
				tmpPair.node1=smallInterval;
				tmpPair.alignmentStart1=nodesArray[nodeID].start+projectionStart>=smallInterval.start?nodesArray[nodeID].start+projectionStart:smallInterval.start;
				tmpPair.alignmentEnd1=nodesArray[nodeID].start+projectionEnd<=smallInterval.end?nodesArray[nodeID].start+projectionEnd:smallInterval.end;
			}
			else {
				tmpPair.node2=smallInterval;
				tmpPair.alignmentStart2=nodesArray[nodeID].start+projectionStart>=smallInterval.start?nodesArray[nodeID].start+projectionStart:smallInterval.start;
				tmpPair.alignmentEnd2=nodesArray[nodeID].start+projectionEnd<=smallInterval.end?nodesArray[nodeID].start+projectionEnd:smallInterval.end;
			}			
			tmpEdge.clear(); tmpEdge.set(tmpPair);
			i=findEdge(smallInterval.nodeID,tmpEdge,false);
			if (i==Math.NEGATIVE_INFINITY) {
				newEdge = new Edge(); tmpEdge.clone(newEdge);
				addEdge(smallInterval.nodeID,newEdge);
				addEdge(neighbor.nodeID,newEdge);
			}
			else {
				neighbors[smallInterval.nodeID][i].addTypes(tmpEdge);
				neighbors[smallInterval.nodeID][i].addDiffs(tmpEdge);
				neighbors[smallInterval.nodeID][i].addOrientation(tmpEdge);
			}			
		}
		else if ( Intervals.isApproximatelyContained(nodesArray[nodeID].start+projectionStart,nodesArray[nodeID].start+projectionEnd,smallStartComplement,smallEndComplement) ||
			      Intervals.areApproximatelyIdentical(nodesArray[nodeID].start+projectionStart,nodesArray[nodeID].start+projectionEnd,smallStartComplement,smallEndComplement)
		        ) {
			if (nodesArray[nodeID].read<neighbor.read) {
				tmpPair.node1=newNode;
				tmpPair.alignmentStart1=nodesArray[nodeID].start+projectionStart>=smallStartComplement?nodesArray[nodeID].start+projectionStart:smallStartComplement;
				tmpPair.alignmentEnd1=nodesArray[nodeID].start+projectionEnd<=smallEndComplement?nodesArray[nodeID].start+projectionEnd:smallEndComplement;
			}
			else {
				tmpPair.node2=newNode;
				tmpPair.alignmentStart2=nodesArray[nodeID].start+projectionStart>=smallStartComplement?nodesArray[nodeID].start+projectionStart:smallStartComplement;
				tmpPair.alignmentEnd2=nodesArray[nodeID].start+projectionEnd<=smallEndComplement?nodesArray[nodeID].start+projectionEnd:smallEndComplement;
			}			
			tmpEdge.clear(); tmpEdge.set(tmpPair);
			i=findEdge(nodeID,tmpEdge,false);
			if (i==Math.NEGATIVE_INFINITY) {
				newEdge = new Edge(); tmpEdge.clone(newEdge);
				addEdge(nodeID,newEdge);
				addEdge(neighbor.nodeID,newEdge);
			}
			else {
				neighbors[nodeID][i].addTypes(tmpEdge);
				neighbors[nodeID][i].addDiffs(tmpEdge);
				neighbors[nodeID][i].addOrientation(tmpEdge);
			}
		}
		else {
			if (nodesArray[nodeID].read<neighbor.read) {
				tmpPair.node1=bigInterval;
				tmpPair.alignmentStart1=nodesArray[nodeID].start+projectionStart;
				tmpPair.alignmentEnd1=nodesArray[nodeID].start+projectionEnd;
			}
			else {
				tmpPair.node2=bigInterval;
				tmpPair.alignmentStart2=nodesArray[nodeID].start+projectionStart;
				tmpPair.alignmentEnd2=nodesArray[nodeID].start+projectionEnd;
			}
			tmpEdge.clear(); tmpEdge.set(tmpPair);
			i=findEdge(bigInterval.nodeID,tmpEdge,false);
			if (i==Math.NEGATIVE_INFINITY) {
				newEdge = new Edge(); tmpEdge.clone(newEdge);
				addEdge(bigInterval.nodeID,newEdge);
				addEdge(neighbor.nodeID,newEdge);
			}
			else {
				neighbors[bigInterval.nodeID][i].addTypes(tmpEdge);
				neighbors[bigInterval.nodeID][i].addDiffs(tmpEdge);
				neighbors[bigInterval.nodeID][i].addOrientation(tmpEdge);
			}
		}		
	}


	/**
	 * Remark: the procedure assumes that all overlap neighbors of a node $i$ form two
	 * (possibly empty) blocks at the beginning of $neighbors[i]$: overlaps with a suffix
	 * of $i$, and overlaps with a prefix of $i$.
	 *
	 * @return out[0]=last overlap edge with a suffix of $node$; out[1]=last overlap edge.
	 */
	private static final void overlapEdgesBoundaries(int node, int[] out) {
		int j;
		int lastSuffixOverlapNeighbor=-1;
		int lastOverlapNeighbor=-1;

		for (j=0; j<nNeighbors[node]; j++) {
			if (neighbors[node][j].overlap==-1) {
				lastOverlapNeighbor=j-1;
				lastSuffixOverlapNeighbor=j-1;
				break;
			}
			if (!neighbors[node][j].isSuffixOverlap(node)) {
				lastSuffixOverlapNeighbor=j-1;
				break;
			}
		}
		for (j=j+1; j<nNeighbors[node]; j++) {
			if (neighbors[node][j].overlap==-1) {
				lastOverlapNeighbor=j-1;
				break;
			}
		}
		out[0]=lastSuffixOverlapNeighbor;
		out[1]=lastOverlapNeighbor;
	}


	/**
	 * Decides whether there are at least two neighbors of $node$ such that: (1) they both
	 * overlap with a suffix of $node$; (2) they neither overlap with one another, nor are
	 * contained in one another; (3) they are caused by the fact that they end with
	 * spurious substrings. See procedure $getSplittingPosition()$ for details.
	 * The procedure searches symmetrically for pairs of neighbors that overlap with a
	 * prefix of $node$. All neighbors of $node$ that are not compatible with at least one
	 * other neighbor of $node$ are marked in array $splittingPosition$.
	 *
	 * Remark: the procedure assumes that e.g. prefix substrings with a spurious end are 
	 * neighbors of a fork. However, it could also happen that one of them is the fork 
	 * itself, and its prefix overlaps with a number of other intervals. If the number of 
	 * such intervals is large enough, such case should have been detected and resolved by 
	 * $fixUndersplits()$. Otherwise, $fixOverlapForks_suffixPrefix()$ attempts at 
	 * resolving them.
	 * Another situation that might create a fork is when a repeat occurs in a read as 
	 * e.g. a suffix of another repeat, but it happens to have low quality, thus it is not
	 * assigned an independent interval, and all the alignments that should be assigned to
	 * it are assigned to the containing interval instead.
	 *
	 * Remark: a hub could also be the full sequence of a repeat, and e.g. many 
	 * incompatible suffix intervals might overlap with its prefix if the repeat was not 
	 * detected inside the suffix intervals. It's not clear how to solve this kind of 
	 * fork, since a clear undersplit signal might be missing from the suffix intervals.
	 *
	 * Remark: for simplicity, the procedure recognizes just a small set of
	 * incompatibility patterns, and it marks $node$ as a fork only if such patterns are
	 * found.
	 *
	 * Remark: if $node$ is of prefix type, and there are multiple incompatible overlap
	 * neighbors of prefix type whose suffixes overlap with a prefix of $node$, the fork 
	 * is correct, since it corresponds to a set of prefix substrings, each with a 
	 * distinct substring at its prefix, but sharing a central substring.
	 *
	 * Remark: the procedure considers only prefix, suffix, prefix-suffix intervals as
	 * candidate forks, but any interval that can overlap a prefix/suffix interval
	 * (i.e. substring type and alignment type) could be considered without changing the
	 * code. The code is also not dependent on the neighbors being of prefix/suffix type,
	 * however intervals of other types (e.g. substring) might be less likely to include a
	 * spurious suffix or prefix.
	 *
	 * Remark: the procedure performs a simplified version of transitive reduction, i.e.
	 * it does not check whether one node can be reached from the other by a sequence of
	 * overlap edges. This is done for simplicity. Thus, the current criterion makes it
	 * more likely to mark nodes as incompatible.
	 *
	 * Remark: the procedure assumes that all overlap neighbors of a node $i$ form two
	 * (possibly empty) blocks at the beginning of $neighbors[i]$: overlaps with a suffix
	 * of $i$, and overlaps with a prefix of $i$.
	 *
	 * @param splittingPosition for each neighbor of $node$: -3=neighbor not considered;
	 * -2=incorrect neighbor; -1=don't split; >=0: splitting position (relative to the
	 * start of the neighbor); if the procedure returns zero, the content of the array is
	 * undefined;
	 * @param positions temporary space, of size at least equal to twice the number of
	 * neighbors of each neighbor of $node$;
	 * @return whether $node$ is an overlap fork on its suffix (=1), prefix (=2), on both
	 * sides (=3), or on no side (=0).
	 */
	private static final int isSpuriousOverlapFork(int node, int lastSuffixOverlapNeighbor, int lastOverlapNeighbor, int[] splittingPosition, int[] positions) {
		boolean isPrefixOverlap1, isSuffixOverlap1, isPrefixOverlap2, isSuffixOverlap2, isNonempty2;
		boolean suffixIncompatible, prefixIncompatible, suffixCorrect, prefixCorrect;
		int j, k, h, p;
		int from, to, lastPosition, first, last;
		int neighbor1, neighbor2, start1, start2, end1, end2, orientation1, orientation2;
		int firstMaximalStart2, lastMaximalEnd2, out;
		Node n1, n2;

		n1=nodesArray[node];
		if (n1.type!=Constants.INTERVAL_DENSE_PREFIX && n1.type!=Constants.INTERVAL_DENSE_SUFFIX && n1.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX && n1.type!=Constants.INTERVAL_DENSE_SINGLEDELETION) return 0;
		suffixIncompatible=false;  // TRUE: >=2 incomp. neighbors on the suffix side.
		prefixIncompatible=false;
		suffixCorrect=true;  // TRUE: no -2 on the suffix side.
		prefixCorrect=true;
		Math.set(splittingPosition,lastOverlapNeighbor,-3);
		first=(n1.type==Constants.INTERVAL_DENSE_PREFIX||n1.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX||n1.type==Constants.INTERVAL_DENSE_SINGLEDELETION)?0:lastSuffixOverlapNeighbor+1;
		last=(n1.type==Constants.INTERVAL_DENSE_SUFFIX||n1.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX||n1.type==Constants.INTERVAL_DENSE_SINGLEDELETION)?lastOverlapNeighbor:lastSuffixOverlapNeighbor;
		for (j=first; j<=last; j++) {
			if (j>lastSuffixOverlapNeighbor) {
				from=lastSuffixOverlapNeighbor+1; 
				to=lastOverlapNeighbor;
			}
			else {
				from=0; 
				to=lastSuffixOverlapNeighbor;
			}
			if (neighbors[node][j].getNAlignments()!=1) continue;
			neighbor1=neighbors[node][j].getTo(node);
			n1=nodesArray[neighbor1];
			if (n1.type!=Constants.INTERVAL_DENSE_PREFIX && n1.type!=Constants.INTERVAL_DENSE_SUFFIX && n1.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX && n1.type!=Constants.INTERVAL_DENSE_SINGLEDELETION) continue;
			start1=n1.start; end1=n1.end;
			isPrefixOverlap1=neighbors[node][j].isPrefixOverlap(neighbor1);
			if (isPrefixOverlap1 && n1.type!=Constants.INTERVAL_DENSE_PREFIX && n1.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX && n1.type!=Constants.INTERVAL_DENSE_SINGLEDELETION) continue;
			isSuffixOverlap1=neighbors[node][j].isSuffixOverlap(neighbor1);
			if (isSuffixOverlap1 && n1.type!=Constants.INTERVAL_DENSE_SUFFIX && n1.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX && n1.type!=Constants.INTERVAL_DENSE_SINGLEDELETION) continue;
			orientation1=neighbors[node][j].orientation;
			lastPosition=-1;
			// Collecting all projections, on $neighbor1$, of other neighbors of $node$
			// that are incompatible with $neighbor1$.
			for (k=from; k<=to; k++) {
				if (k==j || neighbors[node][k].getNAlignments()!=1) continue;
				neighbor2=neighbors[node][k].getTo(node);
				n2=nodesArray[neighbor2];
				if (n2.type!=Constants.INTERVAL_DENSE_PREFIX && n2.type!=Constants.INTERVAL_DENSE_SUFFIX && n2.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX && n2.type!=Constants.INTERVAL_DENSE_SINGLEDELETION) continue;
				start2=n2.start; end2=n2.end;
				firstMaximalStart2=n2.firstMaximalStart;
				lastMaximalEnd2=n2.lastMaximalEnd;
				isPrefixOverlap2=neighbors[node][k].isPrefixOverlap(neighbor2);
				if (isPrefixOverlap2 && n2.type!=Constants.INTERVAL_DENSE_PREFIX && n2.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX && n2.type!=Constants.INTERVAL_DENSE_SINGLEDELETION) continue;
				isSuffixOverlap2=neighbors[node][k].isSuffixOverlap(neighbor2);
				if (isSuffixOverlap2 && n2.type!=Constants.INTERVAL_DENSE_SUFFIX && n2.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX && n2.type!=Constants.INTERVAL_DENSE_SINGLEDELETION) continue;
				orientation2=neighbors[node][k].orientation;
				for (h=0; h<nNeighbors[neighbor1]; h++) {
					if (neighbors[neighbor1][h].getTo(neighbor1)!=neighbor2) continue;
					if (neighborsAreCompatible(neighbors[neighbor1][h],neighbor1,neighbor2,isPrefixOverlap1,isSuffixOverlap1,isPrefixOverlap2,isSuffixOverlap2,orientation1,orientation2)) {
						// $neighbor1$ and $neighbor2$ are compatible: no need to collect
						// projections.
						break;
					}
					if (j<=lastSuffixOverlapNeighbor) suffixIncompatible=true;
					else prefixIncompatible=true;
					if (neighbors[neighbor1][h].getNAlignments()!=1) {
						// $neighbor1$ and $neighbor2$ are not compatible, but there are
						// multiple alignments between them. We don't recognize this
						// pattern.
						splittingPosition[j]=-2;
						if (j<=lastSuffixOverlapNeighbor) suffixCorrect=false;
						else prefixCorrect=false;
						break;
					}
					p=neighbors[neighbor1][h].projectOnto_sharedSubstring(neighbor2,tmpArray,-1);
					if (p<0) {
						// $neighbor1$ and $neighbor2$ are not compatible, but they do not
						// have a shared substring alignment. We don't recognize this
						// pattern.
						splittingPosition[j]=-2;
						if (j<=lastSuffixOverlapNeighbor) suffixCorrect=false;
						else prefixCorrect=false;
						break;
					}
					if (isPrefixOverlap2) isNonempty2=lastMaximalEnd2==-1||tmpArray[1]>=lastMaximalEnd2-IDENTITY_THRESHOLD;
					else isNonempty2=firstMaximalStart2==-1||tmpArray[0]<=firstMaximalStart2+IDENTITY_THRESHOLD;
					p=neighbors[neighbor1][h].projectOnto_sharedSubstring(neighbor1,tmpArray,-1);
					if (p<0) {
						// $neighbor1$ and $neighbor2$ are not compatible, but they
						// do not have a shared substring alignment. We don't
						// recognize this pattern.
						splittingPosition[j]=-2;
						if (j<=lastSuffixOverlapNeighbor) suffixCorrect=false;
						else prefixCorrect=false;
						break;
					}
					if (isPrefixOverlap1) {
						positions[++lastPosition]=tmpArray[1];
						positions[++lastPosition]=isNonempty2?0:1;
					}
					else {
						positions[++lastPosition]=tmpArray[0];
						positions[++lastPosition]=isNonempty2?0:1;
					}
					break;
				}
				if (splittingPosition[j]==-2) break;
			}
			if (splittingPosition[j]==-2) {
				if (j>lastSuffixOverlapNeighbor) return (suffixCorrect&&suffixIncompatible)?1:0;
				else {
					suffixCorrect=false;
					j=lastSuffixOverlapNeighbor+1;
					continue;
				}
			}
			if (lastPosition==-1) continue;
			splittingPosition[j]=getSplittingPosition(n1,positions,lastPosition,isPrefixOverlap1);
			if (splittingPosition[j]==-2) {
				if (j>lastSuffixOverlapNeighbor) return (suffixCorrect&&suffixIncompatible)?1:0;
				else {
					suffixCorrect=false;
					j=lastSuffixOverlapNeighbor+1;
					continue;
				}
			}
		}
		out=0;
		if (suffixCorrect&&suffixIncompatible) out=1;
		if (prefixCorrect&&prefixIncompatible) out|=2;
		return out;
	}


	/**
	 * @param neighbor1,neighbor2 IDs of neighbors of a candidate overlap fork; both
	 * are assumed to overlap with the suffix (respectively, prefix) of the candidate
	 * fork;
	 * @param edge an edge between $neighbor1$ and $neighbor2$, possibly representing
	 * more than one alignment;
	 * @param isPrefixOverlapX,isSuffixOverlapX mark whether the prefix or the suffix of
	 * neighbor X is involved in an overlap with the candidate fork;
	 * @param orientationX orientation of neighbor X with respect to the candidate fork;
	 * @return TRUE iff $edge$ represents an overlap, a containment, or an insertion.
	 */
	private static final boolean neighborsAreCompatible(Edge edge, int neighbor1, int neighbor2, boolean isPrefixOverlap1, boolean isSuffixOverlap1, boolean isPrefixOverlap2, boolean isSuffixOverlap2, int orientation1, int orientation2) {
		// Overlap 1->2
		if ( ( (isPrefixOverlap1 && edge.isSuffixOverlap(neighbor1)) ||
			   (isSuffixOverlap1 && edge.isPrefixOverlap(neighbor1))
			 ) &&
			 ( (isPrefixOverlap2 && edge.isPrefixOverlap(neighbor2)) ||
			   (isSuffixOverlap2 && edge.isSuffixOverlap(neighbor2))
			 )
		   ) return true;
		// Overlap 2->1
		if ( ( (isPrefixOverlap1 && edge.isPrefixOverlap(neighbor1)) ||
			   (isSuffixOverlap1 && edge.isSuffixOverlap(neighbor1))
		     ) &&
		     ( (isPrefixOverlap2 && edge.isSuffixOverlap(neighbor2)) ||
		       (isSuffixOverlap2 && edge.isPrefixOverlap(neighbor2))
		     )
		   ) return true;
		// Containment/insertion 1->2
		if ( ( ((edge.containment==Constants.CONTAINMENT_ONE_IN_TWO||edge.insertion==Constants.INSERTION_ONE_IN_TWO) && edge.nodeID1==neighbor1) ||
		       ((edge.containment==Constants.CONTAINMENT_TWO_IN_ONE||edge.insertion==Constants.INSERTION_TWO_IN_ONE) && edge.nodeID2==neighbor1)
		     ) &&
		     ( edge.containmentSubtype==Constants.CONTAINMENT_SUBSTRING ||
		       (edge.containmentSubtype==Constants.CONTAINMENT_PREFIX && isPrefixOverlap2) ||
		       (edge.containmentSubtype==Constants.CONTAINMENT_SUFFIX && isSuffixOverlap2)
		     ) &&
		     orientationsAreCompatible(orientation1,orientation2,edge.orientation)
		   ) return true;
		// Containment/insertion 2->1
		if ( ( ((edge.containment==Constants.CONTAINMENT_ONE_IN_TWO||edge.insertion==Constants.INSERTION_ONE_IN_TWO) && edge.nodeID1==neighbor2) ||
		       ((edge.containment==Constants.CONTAINMENT_TWO_IN_ONE||edge.insertion==Constants.INSERTION_TWO_IN_ONE) && edge.nodeID2==neighbor2)
		     ) &&
		     ( edge.containmentSubtype==Constants.CONTAINMENT_SUBSTRING ||
		       (edge.containmentSubtype==Constants.CONTAINMENT_PREFIX && isPrefixOverlap1) ||
		       (edge.containmentSubtype==Constants.CONTAINMENT_SUFFIX && isSuffixOverlap1)
		  	 ) &&
		  	 orientationsAreCompatible(orientation1,orientation2,edge.orientation)
		   ) return true;
		return false;
	}


	/**
	 * @param orientation1 orientation of an alignment between interval A and interval B;
	 * @param orientation2 orientation of an alignment between interval A and interval C;
	 * @param orientation3 orientation of an alignment between interval B and interval C;
	 * @return TRUE iff such triplet of orientations is compatible.
	 */
	private static final boolean orientationsAreCompatible(int orientation1, int orientation2, int orientation3) {
		if (orientation1==0) {
			if (orientation2==0 && (orientation3==0 || orientation3==2)) return true;
			else if (orientation2==1 && (orientation3==1 || orientation3==2)) return true;
			else if (orientation2==2 && orientation3!=-1) return true;
		}
		else if (orientation1==1) {
			if (orientation2==0 && (orientation3==1 || orientation3==2)) return true;
			else if (orientation2==1 && (orientation3==0 || orientation3==2)) return true;
			else if (orientation2==2 && orientation3!=-1) return true;
		}
		else if (orientation1==2 && orientation2!=-1 && orientation3!=-1) return true;
		return false;
	}


	/**
	 * Assume that the full copy of the repeat is of prefix type, and that we are
	 * considering neighbors of the fork that overlap with its suffix. The procedure
	 * considers the case in which the fork is a substring (possibly a suffix) of the full
	 * copy, and the neighbors of the fork are suffixes of prefixes of the full copy,
	 * possibly followed by spurious segments not in the full copy. $neighbor$ is
	 * considered correct if either: (1) it has no last maximal event far from the
	 * end (i.e. no spurious end), and no projection, from an incompatible neighbor, is
	 * followed by a maximal event in the incompatible neighbor; or (2) it has a last
	 * maximal event far from the end (i.e. a spurious end), and: (2.1) no projection,
	 * from an incompatible neighbor, that is far from the rightmost projection, is
	 * followed by a maximal event in the incompatible neighbor; (2.2) if all projections
	 * are approximately identical, they coincide with or follow the last maximal end of
	 * $neighbor$.
	 *
	 * Remark: this captures the case in which $neighbor$ is a substring of the full copy
	 * (the only incompatible neighbors are suffixes of shorter prefixes of the full copy,
	 * followed by spurious segments), and a substring followed by a spurious end (the
	 * only incompatible neighbors are substrings of the full copy, suffixes of
	 * longer prefixes followed by a spurious end, and suffixes of shorter prefixes
	 * followed by a spurious end).
	 *
	 * Remark: this does not capture the case in which there are two reference prefix 
	 * substrings with internal modules ABC and ABD, and the fork is a substring of B. 
	 * This is a real fork and is left unchanged.
	 *
	 * Remark: if an incompatible neighbor does not have any alignment with $neighbor$, it
	 * is assumed not to appear in $positions$, and it is disregarded.
	 *
	 * @param positions[0..lastPosition] all projections on $neighbor$ of its incompatible
	 * fork neighbors (relative to $neighbor.start$, not necessarily sorted); each
	 * projection is followed by a one if the corresponding position in the incompatible
	 * neighbor has no maximal event "after" it (in a suitable direction), by a zero
	 * otherwise;
	 * @param prefixOrSuffix substring type of $neighbor$ that should be considered by the
	 * procedure: TRUE=prefix, FALSE=suffix;
	 * @return -2 if $neighbor$ is incorrect; -1 if $neighbor$ is correct, but should
	 * not be split; otherwise, a non-negative splitting position relative to
	 * $neighbor.start$.
	 */
	private static final int getSplittingPosition(Node neighbor, int[] positions, int lastPosition, boolean prefixOrSuffix) {
		final int LOCAL_THRESHOLD = IDENTITY_THRESHOLD+(IDENTITY_THRESHOLD>>1);
		boolean allEmpty;
		int i;
		int min, max, average, firstMaximalStart, lastMaximalEnd, start, end;

		// Setting to -1 maximal boundaries if they are too close to string boundaries
		start=neighbor.start; end=neighbor.end;
		if (prefixOrSuffix && neighbor.lastMaximalEnd!=-1 && neighbor.lastMaximalEnd<end-LOCAL_THRESHOLD) lastMaximalEnd=neighbor.lastMaximalEnd;
		else lastMaximalEnd=-1;
		if (!prefixOrSuffix && neighbor.firstMaximalStart!=-1 && neighbor.firstMaximalStart>start+LOCAL_THRESHOLD) firstMaximalStart=neighbor.firstMaximalStart;
		else firstMaximalStart=-1;
	
		// Checking correctness
		min=Math.POSITIVE_INFINITY; max=0; average=0; allEmpty=true;
		for (i=0; i<=lastPosition; i+=2) {
			average+=positions[i];
			if (positions[i]<min) min=positions[i];
			if (positions[i]>max) max=positions[i];
			allEmpty&=positions[i+1]==1;
		}
		average/=(lastPosition+1)>>1;
		if ((prefixOrSuffix && lastMaximalEnd==-1) || (!prefixOrSuffix && firstMaximalStart==-1)) return allEmpty?-1:-2;		
		if (!allEmpty) {
			for (i=0; i<=lastPosition; i+=2) {
				if (positions[i+1]==1) continue;
				if ((prefixOrSuffix && max-positions[i]>IDENTITY_THRESHOLD) || (!prefixOrSuffix && positions[i]-min>IDENTITY_THRESHOLD)) return -2;
			}
		}
		if (max-min<=IDENTITY_THRESHOLD) {
			if (prefixOrSuffix) {
				if (start+average>lastMaximalEnd-IDENTITY_THRESHOLD) {
					if (start+average<end-IDENTITY_THRESHOLD) return average;
					else return -1;
				}
				else return -2;
			} 
			else {
				if (start+average<firstMaximalStart+IDENTITY_THRESHOLD) {
					if (average>IDENTITY_THRESHOLD) return average;
					else return -1;
				}
				else return -2;
			}
		}
		else {
			if (prefixOrSuffix) {
				if (max>lastMaximalEnd-start) {
					if (max<end-IDENTITY_THRESHOLD) return max;
					else return -1;
				}
				else return lastMaximalEnd-start;
			}
			else {
				if (min<firstMaximalStart-start) {
					if (min>IDENTITY_THRESHOLD) return min;
					else return -1;
				}
				else return firstMaximalStart-start;
			}
		}
	}

	
	/**
	 * Fixes overlap forks induced by suffix substrings that contain a spurious prefix, 
	 * i.e. by suffix substrings whose leftmost alignment was assigned by mistake. If a 
	 * suffix substring with multiple overlaps on its suffix side is such that all its 
	 * overlaps on the suffix side start close to its $firstMaximalStart$ position, and if
	 * $firstMaximalStart$ is sufficiently far from the first position, then a new 
	 * interval is created with coordinates $[firstMaximalStart..end]$, where $end$ is the
	 * last position of the suffix substring, and edges are redirected. An existing node 
	 * is reused if possible. The same holds for prefix and prefix-suffix intervals.
	 *
	 * Remark: the code is a simplified variant of $fixUndersplits()$.
	 *
	 * Remark: the procedure might create new intervals.
	 *
	 * @param isSorted TRUE iff $nodesArray$ is sorted by $read,start,nodeID$ when the
	 * procedure is called;
	 * @return TRUE iff at least one edge has been redirected.
	 */
	private static final boolean fixOverlapForks_suffixPrefix(boolean isSorted) {
		final int MIN_SPURIOUS_LENGTH = IO.quantum<<1;  // Arbitrary
		int i, j, p;
		int maxNewNodes, fromNode, nReusableNodes, nPeaks, lastEvent, previousStats;
		int nCandidatesPrefix, nCandidatesSuffix, nCandidatesPrefixSuffix;
		int nUndersplitNodes, nUndersplitSingletons;
		Node query = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Node node;
		AlignmentPair tmpPair = new AlignmentPair();
		boolean[] undersplitNodeIsNew, peaksLeft, peaksRight;
		int[] stats = new int[3];
		int[] peaks, events, undersplitNodeIsUsed, undersplitNodeIsUsedPrime;
		Node[] undersplitNodes;
		int[][] reusedStats;
		
		// Counting the number of candidate overlap forks
		nCandidatesPrefix=0; nCandidatesSuffix=0; nCandidatesPrefixSuffix=0;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if ( node.length()<MIN_INTERVAL_LENGTH_FOR_UNDERSPLIT ||
				 ( (node.type!=Constants.INTERVAL_DENSE_SUFFIX || node.firstMaximalStart==-1) && 
				   (node.type!=Constants.INTERVAL_DENSE_PREFIX || node.lastMaximalEnd==-1) &&
				   (node.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX || (node.firstMaximalStart==-1 && node.lastMaximalEnd==-1)) &&
				   (node.type!=Constants.INTERVAL_DENSE_SINGLEDELETION || (node.firstMaximalStart==-1 && node.lastMaximalEnd==-1))
				 )
			   ) continue;
			if (node.type==Constants.INTERVAL_DENSE_SUFFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_suffix(i,MIN_SPURIOUS_LENGTH)) nCandidatesSuffix++;
			}
			else if (node.type==Constants.INTERVAL_DENSE_PREFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_prefix(i,MIN_SPURIOUS_LENGTH)) nCandidatesPrefix++;
			}
			else if (node.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_suffix(i,MIN_SPURIOUS_LENGTH) || fixOverlapForks_suffixPrefix_shouldFix_prefix(i,MIN_SPURIOUS_LENGTH)) nCandidatesPrefixSuffix++;
			}
		}
		if (nCandidatesPrefix+nCandidatesSuffix+nCandidatesPrefixSuffix==0) return false;
		maxNewNodes=nCandidatesPrefix+nCandidatesSuffix+(nCandidatesPrefixSuffix<<1);
		expandTo(nNodes+maxNewNodes-1);
		System.err.println("fixUndersplits_suffixPrefix> "+(nCandidatesPrefix+nCandidatesSuffix+nCandidatesPrefixSuffix)+" candidate forks found");
		
		// Loading reusable nodes (using the extra space just allocated in $nodesArray$
		// as temporary space).
		p=nNodes-1;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if ( node.length()<MIN_INTERVAL_LENGTH_FOR_UNDERSPLIT ||
				 ( (node.type!=Constants.INTERVAL_DENSE_SUFFIX || node.firstMaximalStart==-1) && 
				   (node.type!=Constants.INTERVAL_DENSE_PREFIX || node.lastMaximalEnd==-1) &&
				   (node.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX || (node.firstMaximalStart==-1 && node.lastMaximalEnd==-1)) &&
				   (node.type!=Constants.INTERVAL_DENSE_SINGLEDELETION || (node.firstMaximalStart==-1 && node.lastMaximalEnd==-1))
				 )
			   ) continue;
			if (node.type==Constants.INTERVAL_DENSE_SUFFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_suffix(i,MIN_SPURIOUS_LENGTH)) nodesArray[++p]=nodesArray[i];
			}
			else if (node.type==Constants.INTERVAL_DENSE_PREFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_prefix(i,MIN_SPURIOUS_LENGTH)) nodesArray[++p]=nodesArray[i];
			}
			else if (node.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_suffix(i,MIN_SPURIOUS_LENGTH) || fixOverlapForks_suffixPrefix_shouldFix_prefix(i,MIN_SPURIOUS_LENGTH)) nodesArray[++p]=nodesArray[i];
			}
		}
		Node.order=Node.READ_START_NODEID;
		if (!isSorted && nNodes>1) Arrays.sort(nodesArray,0,nNodes);
		if (p+1-nNodes>1) Arrays.sort(nodesArray,nNodes,p+1);
		reusableNodes = new Node[maxNewNodes];  // Remark: reusable nodes might be more than $maxNewNodes$, although this is not likely.
		nReusableNodes=loadReusableNodes(nodesArray,nNodes,p,true);
		System.err.println("fixUndersplits_suffixPrefix> "+nReusableNodes+" reusable nodes found");
		if (!isSorted) {
			Node.order=Node.NODE_ID;
			if (nNodes>1) Arrays.sort(nodesArray,0,nNodes);
		}
		
		// Fixing undersplits
		reusedStats = new int[nNodes][2];
		peaks = new int[4]; peaksLeft = new boolean[4]; peaksRight = new boolean[4];
		events = new int[2];
		undersplitNodes = new Node[maxNewNodes];
		undersplitNodeIsNew = new boolean[undersplitNodes.length];
		undersplitNodeIsUsed = new int[undersplitNodes.length];
		undersplitNodeIsUsedPrime = new int[undersplitNodes.length];
		stats[0]=0; stats[1]=0; stats[2]=0;
		nUndersplitNodes=0; nUndersplitSingletons=0;
		fromNode=nNodes;
		for (i=0; i<nNodes; i++) {
			if (i%10000==0) System.err.println("fixUndersplits_suffixPrefix> "+i+" nodes done ("+IO.getPercent(i,nNodes)+")");
			node=nodesArray[i];
			if ( node.length()<MIN_INTERVAL_LENGTH_FOR_UNDERSPLIT ||
				 ( (node.type!=Constants.INTERVAL_DENSE_SUFFIX || node.firstMaximalStart==-1) && 
				   (node.type!=Constants.INTERVAL_DENSE_PREFIX || node.lastMaximalEnd==-1) &&
				   (node.type!=Constants.INTERVAL_DENSE_PREFIXSUFFIX || (node.firstMaximalStart==-1 && node.lastMaximalEnd==-1)) &&
				   (node.type!=Constants.INTERVAL_DENSE_SINGLEDELETION || (node.firstMaximalStart==-1 && node.lastMaximalEnd==-1))
				 )
			   ) continue;
			nPeaks=0; lastEvent=-1;
			if (node.type==Constants.INTERVAL_DENSE_SUFFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_suffix(i,MIN_SPURIOUS_LENGTH)) {
					peaks[++nPeaks]=node.firstMaximalStart;
					events[++lastEvent]=node.firstMaximalStart-node.start;
					events[++lastEvent]=1;
				}
			}
			else if (node.type==Constants.INTERVAL_DENSE_PREFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_prefix(i,MIN_SPURIOUS_LENGTH)) {
					peaks[++nPeaks]=node.lastMaximalEnd;
					events[++lastEvent]=node.lastMaximalEnd-node.start;
					events[++lastEvent]=0;
				}
			}
			else if (node.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
				if (fixOverlapForks_suffixPrefix_shouldFix_suffix(i,MIN_SPURIOUS_LENGTH)) {
					peaks[++nPeaks]=node.firstMaximalStart;
					events[++lastEvent]=node.firstMaximalStart-node.start;
					events[++lastEvent]=1;
				}
				if (fixOverlapForks_suffixPrefix_shouldFix_prefix(i,MIN_SPURIOUS_LENGTH)) {
					peaks[++nPeaks]=node.lastMaximalEnd;
					events[++lastEvent]=node.lastMaximalEnd-node.start;
					events[++lastEvent]=0;
				}
			}
			if (nPeaks==0) continue;
			peaks[0]=0; peaks[nPeaks+1]=node.length()-1;
			previousStats=stats[2];
			fromNode=fixUndersplits_impl(node,peaks,nPeaks+1,events,lastEvent,reusableNodes,nReusableNodes-1,fromNode,query,stats,peaksLeft,peaksRight,undersplitNodes,undersplitNodeIsNew,undersplitNodeIsUsed,undersplitNodeIsUsedPrime,tmpPair,true)+1;
			if (stats[2]!=previousStats) {
				nUndersplitNodes++;
				if (nNeighbors[i]==0) nUndersplitSingletons++;
				reusedStats[i][1]=1;
				for (j=0; j<undersplitNodes.length; j++) {
					if (undersplitNodes[j]==null) break;
					if (!undersplitNodeIsNew[j] && undersplitNodeIsUsed[j]!=0) reusedStats[undersplitNodes[j].nodeID][0]+=undersplitNodeIsUsed[j];
				}
			}
			else {
				//
				// The following has been commented just because of too much output:
				//System.err.println("Node "+i+" has "+nPeaks+" peaks but no edge moved?! Sorted events:");
				//for (int x=0; x<=(lastTmpPoint<<1)+1; x+=2) System.err.print(sortedEvents[x]+",");
				//for (int x=1; x<=(lastTmpPoint<<1)+1; x+=2) System.err.print("("+sortedEvents[x]+"),");
				//System.err.println();
				//
			}			
		}	
		System.err.println("fixUndersplits_suffixPrefix> Merging undersplit nodes...");
		fromNode=mergeUndersplitNodes(nNodes,fromNode-1)+1;
		System.err.println("fixUndersplits_suffixPrefix> Nodes subjected to fix: "+nUndersplitNodes+" ("+IO.getPercent(nUndersplitNodes,nNodes)+"%).");
		System.err.println("fixUndersplits_suffixPrefix> Disconnected nodes created by fix: "+nUndersplitSingletons+" ("+IO.getPercent(nUndersplitSingletons,nUndersplitNodes)+"%).");
		System.err.println("fixUndersplits_suffixPrefix> New nodes created by fix: "+(fromNode-nNodes)+" ("+IO.getPercent(fromNode-nNodes,nNodes)+"%).");
		System.err.println("fixUndersplits_suffixPrefix> Edges redirected by fix: "+stats[2]);
		//System.err.println("Old nodes actually reused as undersplit nodes:");
		j=0;
		for (i=0; i<nNodes; i++) {
			if (reusedStats[i][0]!=0) {
				j++;
				//
				// The following has been commented just because of too much output:
				//System.err.println("node "+i+": "+reusedStats[i][0]+" edges redirected to it; subjected to undersplit resolution: "+reusedStats[i][1]);
				//
			}
		}
		System.err.println("fixUndersplits_suffixPrefix> "+j+" total nodes reused");
		nNodes=fromNode;
		return nUndersplitNodes>0;
	}
	
	
	/**
	 * @param minSpuriousLength minimum length of a spurious prefix;
	 * @return TRUE iff $nodeID$ should be fixed on its suffix side, i.e. if all overlaps
	 * are only on its suffix side, if there are at least two overlaps, and if their 
	 * projections are all approximately identical to [firstMaximalStart..end].
	 */
	private static final boolean fixOverlapForks_suffixPrefix_shouldFix_suffix(int nodeID, int minSpuriousLength) {
		final int THRESHOLD = IO.quantum;
		final int firstMaximalStart = nodesArray[nodeID].firstMaximalStart;
		int j;
		int nOverhangs, overlapEnd;
		Edge edge;
		
		if (nodesArray[nodeID].firstMaximalStart==-1 || nodesArray[nodeID].firstMaximalStart<nodesArray[nodeID].start+minSpuriousLength) return false;
		nOverhangs=0;
		for (j=0; j<nNeighbors[nodeID]; j++) {
			edge=neighbors[nodeID][j];
			if (edge.overlap==-1) continue;
			overlapEnd=edge.getOverlapEnd(nodeID);
			if (overlapEnd==0 || overlapEnd==2) return false;
			nOverhangs++;
			if (edge.nodeID1==nodeID) {
				if ((edge.overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) {
					edge.projectOnto_overlap(nodeID,tmpArray,Constants.OVERLAP_SUFFIX_PREFIX,-1);
					if (Math.abs(firstMaximalStart,tmpArray[0])>THRESHOLD) return false;
				}
				if ((edge.overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) {
					edge.projectOnto_overlap(nodeID,tmpArray,Constants.OVERLAP_SUFFIX_SUFFIX,-1);
					if (Math.abs(firstMaximalStart,tmpArray[0])>THRESHOLD) return false;
				}
			}
			else {
				if ((edge.overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) {
					edge.projectOnto_overlap(nodeID,tmpArray,Constants.OVERLAP_PREFIX_SUFFIX,-1);
					if (Math.abs(firstMaximalStart,tmpArray[0])>THRESHOLD) return false;
				}
				if ((edge.overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) {
					edge.projectOnto_overlap(nodeID,tmpArray,Constants.OVERLAP_SUFFIX_SUFFIX,-1);
					if (Math.abs(firstMaximalStart,tmpArray[0])>THRESHOLD) return false;
				}
			}
		}
		return nOverhangs>1;
	}
	
	
	/**
	 * @param minSpuriousLength minimum length of a spurious prefix;
	 * @return TRUE iff $nodeID$ should be fixed on its prefix side, i.e. if all overlaps
	 * are only on its prefix side, if there are at least two overlaps, and if their 
	 * projections are all approximately identical to [start..lastMaximalEnd].
	 */
	private static final boolean fixOverlapForks_suffixPrefix_shouldFix_prefix(int nodeID, int minSpuriousLength) {
		final int THRESHOLD = IO.quantum;
		final int lastMaximalEnd = nodesArray[nodeID].lastMaximalEnd;
		int j;
		int nOverhangs, overlapEnd;
		Edge edge;
		
		if (nodesArray[nodeID].lastMaximalEnd==-1 || nodesArray[nodeID].lastMaximalEnd>nodesArray[nodeID].end-minSpuriousLength) return false;
		nOverhangs=0;
		for (j=0; j<nNeighbors[nodeID]; j++) {
			edge=neighbors[nodeID][j];
			if (edge.overlap==-1) continue;
			overlapEnd=edge.getOverlapEnd(nodeID);
			if (overlapEnd==1 || overlapEnd==2) return false;
			nOverhangs++;
			if (edge.nodeID1==nodeID) {
				if ((edge.overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) {
					edge.projectOnto_overlap(nodeID,tmpArray,Constants.OVERLAP_PREFIX_PREFIX,-1);
					if (Math.abs(lastMaximalEnd,tmpArray[1])>THRESHOLD) return false;
				}
				if ((edge.overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) {
					edge.projectOnto_overlap(nodeID,tmpArray,Constants.OVERLAP_PREFIX_SUFFIX,-1);
					if (Math.abs(lastMaximalEnd,tmpArray[1])>THRESHOLD) return false;
				}
			}
			else {
				if ((edge.overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) {
					edge.projectOnto_overlap(nodeID,tmpArray,Constants.OVERLAP_PREFIX_PREFIX,-1);
					if (Math.abs(lastMaximalEnd,tmpArray[1])>THRESHOLD) return false;
				}
				if ((edge.overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) {
					edge.projectOnto_overlap(nodeID,tmpArray,Constants.OVERLAP_PREFIX_SUFFIX,-1);
					if (Math.abs(lastMaximalEnd,tmpArray[1])>THRESHOLD) return false;
				}
			}
		}
		return nOverhangs>1;
	}
	
	
	
	
	
	
	
	
	// ---------------------------------- SERIALIZATION ----------------------------------
	
	/**
	 * Serializes the whole graph.
	 */
	public static final void serialize(String path) throws IOException {
		int i, j;
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(path),IO.BUFFER_SIZE);
		
		IO.writeInt(maxAlignmentsPerRead,out);
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=-1;
		}
		IO.writeInt(nNodes,out);
		for (i=0; i<nNodes; i++) {
			nodesArray[i].serialize(out);
			IO.writeInt(nNeighbors[i],out);
			if (nNeighbors[i]==0) continue;
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].printed==1) {
					out.write(1);
					IO.writeInt(neighbors[i][j].getTo(i),out);
					IO.writeInt(neighbors[i][j].supplement?1:0,out);
				}
				else {
					out.write(0);
					neighbors[i][j].serialize(out);
					neighbors[i][j].printed=1;
				}
			}
		}
		out.close();
	}
	
	
	/**
	 * Wraps $serializeComponents()$ to save only frequent components.
	 *
	 * Remark: the procedure prints kernel tags files as specified by $IntervalGraphStep3.
	 * printKernelTags()$ for all nodes that do not belong to any serialized component. 
	 * If a node belongs to some non-serialized components, all such components are used 
	 * as kernel tags. Otherwise, the node is printed as if it had no kernel tag.
	 *
	 * @param components if NULL, considers all components for serialization;
	 * @return the number of components that have been serialized by the procedure.
	 */
	public static final int serializeFrequentComponents(int[] components, int nComponents, boolean strict, boolean discardDisconnected, String outputDir, String suffix, String tagsDir, String tagsPrefix) throws IOException {
		final int THRESHOLD_SIZE = IO.minOccurrencesInGenome*IO.coverage;
		final int THRESHOLD_MAXIMAL = (THRESHOLD_SIZE)<<1;
		final int ARTIFICIAL_KERNEL = Math.POSITIVE_INFINITY;
		int i, j, k;
		int nFrequentComponents, nPrintedComponents;
		int nSerialized, nNotSerialized, nodesWithTags, component;
		Node node;
		boolean[] isComponentFrequent, isComponentLarge, isComponentPeriodic, printComponent;
		int[] size, nMaximal;
		
		// Collecting statistics
		size = new int[nComponents];
		nMaximal = new int[nComponents];
		getComponentStats(components,nComponents,strict,discardDisconnected,size,nMaximal);		
		isComponentFrequent = new boolean[nComponents];
		isComponentLarge = new boolean[nComponents];
		Math.set(isComponentFrequent,nComponents-1,false);
		Math.set(isComponentLarge,nComponents-1,false);
		nFrequentComponents=0;
		for (i=0; i<nComponents; i++) {
			if (nMaximal[i]>=THRESHOLD_MAXIMAL) {
				isComponentFrequent[i]=true;
				nFrequentComponents++;
			}
			if (size[i]>=THRESHOLD_SIZE) isComponentLarge[i]=true;
		}
		isComponentPeriodic = new boolean[nComponents];
		markPeriodicComponents(components,nComponents-1,isComponentPeriodic,strict,0);
		
		// Serializing components
		printComponent = new boolean[nComponents];
		Math.set(printComponent,nComponents-1,false);
		nPrintedComponents=0;
		for (i=0; i<nComponents; i++) {
			if (isComponentFrequent[i] || (isComponentPeriodic[i] && isComponentLarge[i])) {
				printComponent[i]=true;
				nPrintedComponents++;
			}
		}
		serializeComponents(components,nComponents,size,printComponent,strict,discardDisconnected,outputDir,suffix);
		
		// Printing kernel tags file
		nodesWithTags=0;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			nSerialized=0; nNotSerialized=0;
			for (j=0; j<=node.lastComponent; j++) {
				component=node.components[j];
				if (component==-1) continue;
				if (printComponent[component]) nSerialized++;
				else nNotSerialized++;
			}
			if (nSerialized!=0) {
				node.lastKernel=-1;
				continue;
			}
			nodesWithTags++;
			node.lastPathWithStart=-1; node.lastPathWithEnd=-1;
			if (nNotSerialized==0) {
				node.lastKernel=0;
				if (node.kernels==null || node.kernels.length==0) node.kernels = new int[1];
				node.kernels[0]=ARTIFICIAL_KERNEL;
				if (node.kernelOrientations==null || node.kernelOrientations.length==0) node.kernelOrientations = new byte[1];
				node.kernelOrientations[0]=(byte)0;
			}
			else {
				if (node.kernels==null || node.kernels.length<nNotSerialized) node.kernels = new int[nNotSerialized];
				k=-1;
				for (j=0; j<=node.lastComponent; j++) {
					component=node.components[j];
					if (component!=-1) node.kernels[++k]=component;
				}
				node.lastKernel=k;
				if (node.kernelOrientations==null || node.kernelOrientations.length<k+1) node.kernelOrientations = new byte[k+1];
				Math.set(node.kernelOrientations,k,(byte)0);
			}
		}
		if (nodesWithTags>0) {
			IntervalGraphStep3.pathKernelLengths=null;
			IntervalGraphStep3.printKernelTags(tagsDir+"/"+IO.TAGS_PREFIX+(tagsPrefix.length()==0?"root":tagsPrefix)+".txt",tagsPrefix,true,true,ARTIFICIAL_KERNEL);
		}
		
		return nPrintedComponents;
	}
	
	
	/**
	 * Like $serialize()$, but stores only nodes whose component is in $componentIDs$,
	 * each node in potentially multiple files, associated with all the components to 
	 * which it belongs. If $strict=false$, the procedure saves also nodes with more than 
	 * one component: the motivation for this is that a large cluster might be mostly
	 * frontier, and just few of its nodes might have a single component. Otherwise, only 
	 * nodes with a single component are stored, and each node is stored in exactly one 
	 * file.
	 *
	 * Remark: the serialized nodes in each file contain the $nodeID$ values in the 
	 * original graph, rather than in the new, smaller graph.
	 *
	 * @param componentIDs if NULL, all components are serialized;
	 * @param printComponent one flag per element in $componentIDs$; components not marked
	 * in this array are discarded; if the array is null, no component is discarded;
	 * @param componentIDs sorted set of distinct connected component IDs; only components
	 * in this set are considered;
	 * @param discardDisconnected a node that has no ON neighbor in component C is not 
	 * stored in component C.
	 */
	public static final void serializeComponents(int[] componentIDs, int nComponents, int[] componentSizes, boolean[] printComponent, boolean strict, boolean discardDisconnected, String outputDir, String suffix) throws IOException {
		final int N_OUTPUT_FILES = IO.MAX_OPEN_FILES-1;
		final int N_PRINTED_LONGS;
		boolean found;
		int i, j, k, h, n;
		int component, firstComponent, nNeighborsInComponent;
		Node neighbor;
		BufferedOutputStream[] outputFiles;
		
		// Marking nodes to be discarded in specific components
		if (discardDisconnected) {
			for (i=0; i<nNodes; i++) {
				if (nodesArray[i].lastComponent==-1) continue;
				if (strict && nodesArray[i].lastComponent>0) continue;
				if (nodesArray[i].lastComponent==0 && nodesArray[i].components[0]==-1) continue;
				for (j=0; j<=nodesArray[i].lastComponent; j++) {
					component=nodesArray[i].components[j];
					found=false;
					for (h=0; h<nNeighbors[i]; h++) {
						n=neighbors[i][h].getTo(i);
						if ( neighbors[i][h].on && 
							 (strict?nodesArray[n].lastComponent==0:true) &&
							 Math.linearSearch_unsorted(nodesArray[n].components,0,nodesArray[n].lastComponent+1,component,true)>=0
						   ) {
							found=true;
							break;
						}
					}
					if (!found) nodesArray[i].components[j]=-2-component;
				}
			}
		}
		
		// Sorting all components (some of which are negative) of all nodes
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].lastComponent<=0) continue;
			if (nodesArray[i].lastComponent>0) Arrays.sort(nodesArray[i].components,0,nodesArray[i].lastComponent+1);
		}
		
		// Allocating $printedComponents$ flags on edges
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				if (neighbors[i][j].printedComponents!=null) continue;
				k=Math.ceil(nodesArray[neighbors[i][j].nodeID1].lastComponent+1,64);
				neighbors[i][j].printedComponents = new long[k];
				Math.set(neighbors[i][j].printedComponents,k-1,0L);
			}
		}
		
		// Writing to buffers
		outputFiles = new BufferedOutputStream[nComponents];
		firstComponent=0;
		while (firstComponent<nComponents) {
			for (i=0; i<N_OUTPUT_FILES && firstComponent+i<nComponents; i++) {
				if (printComponent!=null && !printComponent[firstComponent+i]) continue;
				outputFiles[firstComponent+i] = new BufferedOutputStream(new FileOutputStream(outputDir+"/"+(componentIDs==null?firstComponent+i:componentIDs[firstComponent+i])+suffix),IO.BUFFER_SIZE);
				IO.writeInt(maxAlignmentsPerRead,outputFiles[firstComponent+i]);
				IO.writeInt(componentSizes[firstComponent+i],outputFiles[firstComponent+i]);
			}
			for (i=0; i<nNodes; i++) {
				if (nodesArray[i].lastComponent==-1) {
					if (i%10000==0) System.err.println("Processed node "+i);
					continue;
				}
				if (nodesArray[i].lastComponent==0 && nodesArray[i].components[0]==-1) {
					if (i%10000==0) System.err.println("Processed node "+i);
					continue;
				}
				if (strict && nodesArray[i].lastComponent>0) {
					if (i%10000==0) System.err.println("Processed node "+i);
					continue;
				}
				for (j=0; j<=nodesArray[i].lastComponent; j++) {
					component=nodesArray[i].components[j];
					if (component<0) continue;
					k=componentIDs==null?component:Arrays.binarySearch(componentIDs,component);
					if (k<firstComponent || k>=firstComponent+N_OUTPUT_FILES || (printComponent!=null&&!printComponent[k])) continue;
					nodesArray[i].serialize(outputFiles[k]);
					nNeighborsInComponent=0;
					for (h=0; h<nNeighbors[i]; h++) {
						neighbor=nodesArray[neighbors[i][h].getTo(i)];
						if (Arrays.binarySearch(neighbor.components,0,neighbor.lastComponent+1,component)>=0 && (strict?neighbor.lastComponent==0:true)) nNeighborsInComponent++;
					}
					IO.writeInt(nNeighborsInComponent,outputFiles[k]);					
					if (nNeighborsInComponent==0) continue;
					for (h=0; h<nNeighbors[i]; h++) {
						neighbor=nodesArray[neighbors[i][h].getTo(i)];
						if (Arrays.binarySearch(neighbor.components,0,neighbor.lastComponent+1,component)<0 || (strict?neighbor.lastComponent!=0:false)) continue;
						if (neighbors[i][h].isPrintedInComponent(component)) {
							outputFiles[k].write(1);
							IO.writeInt(neighbors[i][h].getTo(i),outputFiles[k]);
							IO.writeInt(neighbors[i][h].supplement?1:0,outputFiles[k]);
						}
						else {
							outputFiles[k].write(0);
							neighbors[i][h].serialize(outputFiles[k]);
							neighbors[i][h].setPrintedInComponent(component);
						}
					}
				}
				if (i%10000==0) System.err.println("Processed node "+i);
			}
			for (i=0; i<N_OUTPUT_FILES && firstComponent+i<nComponents; i++) {
				if (outputFiles[firstComponent+i]!=null) {
					outputFiles[firstComponent+i].close();
					outputFiles[firstComponent+i]=null;
				}
			}
			firstComponent+=N_OUTPUT_FILES;
		}
		
		// Deallocating $printedComponents$ flags on edges
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printedComponents=null;
		}
	}

	
	/**
	 * @param resetNodeIDs TRUE: the $nodeID$ fields in the serialized file are 
	 * overwritten with new values, relative to the new file;
	 * @param keepSameReadEdges FALSE: removes all edges between nodes that belong to the
	 * same read;
	 * @return if $resetNodeIDs=TRUE$, the original $nodeID$ field of the $i$-th node, for
	 * every $i$; null otherwise.
	 */
	public static final int[] deserialize(String path, boolean keepSameReadEdges, boolean resetNodeIDs) throws IOException {
		boolean supplement;
		int i, j, k;
		int deserialized, otherNode, otherNodePrime, read, nodeID;
		int[] new2old;
		Node tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		BufferedInputStream out = new BufferedInputStream(new FileInputStream(path),IO.BUFFER_SIZE);
		
		// Loading file
		maxAlignmentsPerRead=IO.readInt(out);	
		Node.order=Node.NODE_ID;
		nNodes=IO.readInt(out);
		if (nodesArray==null || nodesArray.length<nNodes) nodesArray = new Node[nNodes];
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i]==null) nodesArray[i] = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		}
		nNeighbors = new int[nNodes];
		neighbors = new Edge[nNodes][0];
		for (i=0; i<nNodes; i++) {
			nodesArray[i].deserialize(out);	
			nodeID=nodesArray[i].nodeID;
			nNeighbors[i]=IO.readInt(out);
			if (nNeighbors[i]==0) {
				System.err.println("Node "+i+" has no neighbor");
				continue;
			}
			if (nNeighbors[i]<0) {
				System.err.println("ERROR 0: negative number of neighbors?! "+nNeighbors[i]+" for node "+nodeID+", i="+i+", nNodes="+nNodes+", file="+path);
				System.exit(1);
			}
			neighbors[i] = new Edge[nNeighbors[i]];
			for (j=0; j<nNeighbors[i]; j++) {
				deserialized=out.read();
				if (deserialized==1) {
					otherNode=IO.readInt(out);
					supplement=IO.readInt(out)==1;
					tmpNode.nodeID=otherNode;
					otherNodePrime=Arrays.binarySearch(nodesArray,0,i,tmpNode);
					if (otherNodePrime<0) {
						System.err.println("deserialize> ERROR 1: node "+otherNode+" not found.");
						System.exit(1);
					}
					for (k=0; k<nNeighbors[otherNodePrime]; k++) {
						if (neighbors[otherNodePrime][k]==null) continue;
						if (neighbors[otherNodePrime][k].getTo(otherNode)==nodeID && neighbors[otherNodePrime][k].supplement==supplement) {
							neighbors[i][j]=neighbors[otherNodePrime][k];
							break;
						}
					}
					if (neighbors[i][j]==null) {
						System.err.println("deserialize> ERROR 1.5: edge not found");
						System.exit(1);
					}
				}
				else {
					neighbors[i][j] = new Edge();
					neighbors[i][j].deserialize(out);
				}
			}
		}
		out.close();
		
		// Removing same-read edges
		if (!keepSameReadEdges) {
			for (i=0; i<nNodes; i++) {
				nodeID=nodesArray[i].nodeID; read=nodesArray[i].read;
				for (j=0; j<nNeighbors[i]; j++) {
					otherNode=neighbors[i][j].getTo(nodeID);
					tmpNode.nodeID=otherNode;
					otherNodePrime=Arrays.binarySearch(nodesArray,0,nNodes,tmpNode);
					if (otherNodePrime<0) {
						System.err.println("deserialize> ERROR 2: node "+otherNode+" not found.");
						System.exit(1);
					}
					if (nodesArray[otherNodePrime].read==read) neighbors[i][j]=null;
				}
				k=-1;
				for (j=0; j<nNeighbors[i]; j++) {
					if (neighbors[i][j]!=null) neighbors[i][++k]=neighbors[i][j];
				}
				nNeighbors[i]=k+1;
			}
			// Removing nodes with no neighbor
			j=-1;
			for (i=0; i<nNodes; i++) {
				if (nNeighbors[i]==0) continue;
				j++;
				if (j!=i) {
					nodesArray[j]=nodesArray[i];
					if (neighbors[j].length<nNeighbors[i]) neighbors[j] = new Edge[nNeighbors[i]];
					System.arraycopy(neighbors[i],0,neighbors[j],0,nNeighbors[i]);
					nNeighbors[j]=nNeighbors[i];
				}
			}
			nNodes=j+1;
		}
		
		// Resetting $nodeID$ values
		return resetNodeIDs?deserialize_resetNodeIDs(tmpNode,true):null;
	}
	
	
	/**
	 * Like $deserialize()$, but does not load edges in memory.
	 *
	 * @param reuseNodes reuses the existing $Node$ objects in $nodesArray$, if any.
	 */
	public static final int[] deserialize_noEdges(String path, boolean resetNodeIDs, boolean reuseNodes) throws IOException {
		int i, j;
		int deserialized, nodeID, nn;
		Node tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		BufferedInputStream out = new BufferedInputStream(new FileInputStream(path),IO.BUFFER_SIZE);
		
		maxAlignmentsPerRead=IO.readInt(out);	
		Node.order=Node.NODE_ID;
		nNodes=IO.readInt(out);
		if (nodesArray==null) {
			nodesArray = new Node[nNodes];
			for (i=0; i<nNodes; i++) nodesArray[i] = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		}
		else if (nodesArray.length<nNodes) {
			if (reuseNodes) {
				Node[] newArray = new Node[nNodes];
				System.arraycopy(nodesArray,0,newArray,0,nodesArray.length);
				for (i=nodesArray.length; i<nNodes; i++) newArray[i] = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
				nodesArray=newArray;
			}
			else {
				nodesArray = new Node[nNodes];
				for (i=0; i<nNodes; i++) nodesArray[i] = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
			}
		}
		else if (!reuseNodes) {
			for (i=0; i<nNodes; i++) nodesArray[i] = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		}
		if (nNeighbors==null || nNodes>nNeighbors.length) nNeighbors = new int[nNodes];
		neighbors=null;
		for (i=0; i<nNodes; i++) {
			nodesArray[i].deserialize(out);	
			nodeID=nodesArray[i].nodeID;
			nNeighbors[i]=0;
			nn=IO.readInt(out);
			if (nn==0) {
				System.err.println("Node "+i+" has no neighbor");
				continue;
			}
			if (nn<0) {
				System.err.println("deserialize_noEdges> ERROR: negative number of neighbors?! "+nn+" for node "+nodeID+", i="+i+", nNodes="+nNodes+", file="+path);
				System.exit(1);
			}
			for (j=0; j<nn; j++) {
				deserialized=out.read();
				if (deserialized==1) {
					IO.readInt(out);
					IO.readInt(out);
				}
				else Edge.skip(out);
			}
		}
		out.close();
		return resetNodeIDs?deserialize_resetNodeIDs(tmpNode,false):null;
	}
	
	
	/**
	 * Overwrites the $nodeID$ fields from the serialized file with new values, relative 
	 * to the new file.
	 *
	 * @return the original $nodeID$ field of the $i$-th node, for every $i$.
	 */
	public static final int[] deserialize_resetNodeIDs(Node tmpNode, boolean edgesExist) {
		int i, j;
		int nodeID, otherNode, otherNodePrime;
		int[] new2old;
		
		new2old = new int[nNodes];
		if (edgesExist) {
			for (i=0; i<nNodes; i++) {
				for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=-1;
			}
			for (i=0; i<nNodes; i++) {
				nodeID=nodesArray[i].nodeID;
				for (j=0; j<nNeighbors[i]; j++) {
					if (neighbors[i][j].printed==1) continue;
					otherNode=neighbors[i][j].getTo(nodeID);
					tmpNode.nodeID=otherNode;
					otherNodePrime=Arrays.binarySearch(nodesArray,0,nNodes,tmpNode);
					if (otherNodePrime<0) {
						System.err.println("deserialize_resetNodeIDs> ERROR 1: node "+otherNode+" not found.");
						System.exit(1);
					}
					if (neighbors[i][j].nodeID1==otherNode) {
						neighbors[i][j].nodeID1=otherNodePrime;
						neighbors[i][j].nodeID2=i;
					}
					else {
						neighbors[i][j].nodeID2=otherNodePrime;
						neighbors[i][j].nodeID1=i;
					}
					neighbors[i][j].printed=1;
				}
			}
		}
		for (i=0; i<nNodes; i++) {
			new2old[i]=nodesArray[i].nodeID;
			nodesArray[i].nodeID=i;
		}
		return new2old;
	}
	
	
	public static final void printNew2OldArray(int[] array, String path) throws IOException {
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(path),IO.BUFFER_SIZE);
		for (int i=0; i<array.length; i++) IO.writeInt(array[i],out);
		out.close();
	}
	
	
	/**
	 * Stores in $nNeighbors[0..X]$ the list of distinct reads in the graph, where X is
	 * returned in output.
	 *
	 * @param sorted TRUE: nodes are already sorted by read.
	 */
	public static final int getDistinctReads(boolean sorted) {
		int i, j;
		int previous;
		
		for (i=0; i<nNodes; i++) nNeighbors[i]=nodesArray[i].read;
		if (!sorted) Arrays.sort(nNeighbors,0,nNodes);
		j=0; previous=nNeighbors[0];
		for (i=1; i<nNodes; i++) {
			if (nNeighbors[i]==previous) continue;
			previous=nNeighbors[i];
			nNeighbors[++j]=nNeighbors[i];
		}
		return j;
	}
	
	
	/**
	 * Uses the connected components of the graph to partition every file that describes 
	 * reads and alignments: the information about a read is written to the files of every
	 * connected component that contains at least one node in the read; the information 
	 * about an alignment is written to the files of every connected component that 
	 * contains at least one node in each read involved in the alignment.
	 * Using this procedure rather than the other $filter*()$ procedures makes the whole 
	 * pipeline faster, since this procedure requires just a few passes over the files to
	 * be filtered.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read$, and it assumes
	 * $Reads.readIDs$ and $Reads.nReads$ to be correctly initialized (the former possibly
	 * to NULL).
	 *
	 * @param components sorted; if NULL, considers all components rather than a subset;
	 * @param strict if FALSE, the procedure uses also nodes with more than one component
	 * to decide the mapping between reads and components; if TRUE, only nodes with a 
	 * single component are used.
	 */
	public static final void filterAlignmentsAndReads(int[] componentIDs, int nComponents, boolean strict, String readLengthsFilePath, String readPhredFilePath, String readShortPeriodFilePath, String alignmentFilePath, String shortPeriodFilePath, String outputDir) throws IOException {
		final int INCREMENT = 10;  // Arbitrary
		final boolean FILTER_QUALITIES = readPhredFilePath!=null && (new File(readPhredFilePath)).exists();
		int i, j, r;
		int last, read, readA, readB, component, max, nOutputFiles, firstComponent;
		final int nReads = Reads.nReads;
		Node node;
		String str1, str2;
		BufferedReader readLengthsFile, readPhredFile, readShortPeriodFile;
		BufferedReader alignmentFile, shortPeriodFile;
		BufferedWriter[] readIDsFiles, readLengthsFiles, readPhredFiles, readShortPeriodFiles;
		BufferedWriter[] alignmentFiles, shortPeriodFiles;
		int[] lastComponent2read, lastRead2component;
		int[][] component2read, read2component;
		
		// Building $read2component$.
		lastComponent2read = new int[nComponents];
		Math.set(lastComponent2read,nComponents-1,-1);
		component2read = new int[nComponents][INCREMENT];
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			last=node.lastComponent;
			if (last==-1 || (last==0 && node.components[0]==-1) || (strict && last>0)) continue;
			read=node.read;
			for (j=0; j<=last; j++) {
				component=node.components[j];
				if (componentIDs!=null) {
					component=Arrays.binarySearch(componentIDs,0,nComponents,component);
					if (component<0) continue;
				}
				if (lastComponent2read[component]!=-1 && read==component2read[component][lastComponent2read[component]]) continue;	
				lastComponent2read[component]++;
				if (lastComponent2read[component]==component2read[component].length) {
					int[] newArray = new int[component2read[component].length+INCREMENT];
					System.arraycopy(component2read[component],0,newArray,0,component2read[component].length);
					component2read[component]=newArray;
				}
				component2read[component][lastComponent2read[component]]=read;
			}
			if (i%10000==0) System.err.println("Processed node "+i);
		}
		lastRead2component = new int[nReads];
		Math.set(lastRead2component,nReads-1,-1);
		read2component = new int[nReads][INCREMENT];
		for (i=0; i<nComponents; i++) {
			component=componentIDs==null?i:componentIDs[i];
			for (j=0; j<=lastComponent2read[i]; j++) {
				read=Reads.indexOfRead(component2read[i][j]);
				lastRead2component[read]++;
				if (lastRead2component[read]==read2component[read].length) {
					int[] newArray = new int[read2component[read].length+INCREMENT];
					System.arraycopy(read2component[read],0,newArray,0,read2component[read].length);
					read2component[read]=newArray;
				}
				read2component[read][lastRead2component[read]]=component;
			}
			if (i%10==0) System.err.println("Processed component "+i);
		}
		lastComponent2read=null;
		for (i=0; i<nComponents; i++) component2read[i]=null;
		component2read=null;
		max=0;
		for (i=0; i<nReads; i++) max=Math.max(max,lastRead2component[i]);
		max++;
		if (tmpArray==null || tmpArray.length<max) tmpArray = new int[max];
		
		// Filtering read IDs, lengths and qualities.
		readIDsFiles = new BufferedWriter[nComponents];
		readLengthsFiles = new BufferedWriter[nComponents];
		if (FILTER_QUALITIES) readPhredFiles = new BufferedWriter[nComponents];
		else readPhredFiles=null;
		nOutputFiles=(IO.MAX_OPEN_FILES-2)/3; firstComponent=0;
		while (firstComponent<nComponents) {
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				readIDsFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+(firstComponent+i)+"-"+IO.READS_IDS),IO.BUFFER_SIZE);
				readLengthsFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+(firstComponent+i)+"-"+IO.READS_LENGTHS),IO.BUFFER_SIZE);
				if (FILTER_QUALITIES) readPhredFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+(firstComponent+i)+"-"+IO.READS_PHRED_PREFIX+readPhredFilePath.substring(readPhredFilePath.lastIndexOf("."))),IO.BUFFER_SIZE);
			}
			readLengthsFile = new BufferedReader(new FileReader(readLengthsFilePath),IO.BUFFER_SIZE);
			if (FILTER_QUALITIES) readPhredFile = new BufferedReader(new FileReader(readPhredFilePath),IO.BUFFER_SIZE);
			else readPhredFile=null;
			str1=readLengthsFile.readLine();
			str2=FILTER_QUALITIES?readPhredFile.readLine():null;
			for (read=0; read<nReads; read++) {
				last=lastRead2component[read];
				for (i=0; i<=last; i++) {
					component=read2component[read][i];
					if (component>=firstComponent && component<firstComponent+nOutputFiles) {
						readIDsFiles[component].write(Reads.readAtIndex(read)+"\n");
						readLengthsFiles[component].write(str1+"\n");
						if (FILTER_QUALITIES) readPhredFiles[component].write(str2+"\n");
					}
				}
				str1=readLengthsFile.readLine();
				str2=FILTER_QUALITIES?readPhredFile.readLine():null;
				if (read%10000==0) System.err.println("Processed read "+read);
			}
			readLengthsFile.close();
			if (FILTER_QUALITIES) readPhredFile.close();
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				readIDsFiles[firstComponent+i].close();
				readIDsFiles[firstComponent+i]=null;
				readLengthsFiles[firstComponent+i].close();
				readLengthsFiles[firstComponent+i]=null;
				if (FILTER_QUALITIES) {
					readPhredFiles[firstComponent+i].close();
					readPhredFiles[firstComponent+i]=null;
				}
			}
			firstComponent+=nOutputFiles;
		}
		
		// Filtering short-period tracks
		readShortPeriodFiles = new BufferedWriter[nComponents];
		nOutputFiles=IO.MAX_OPEN_FILES-1; firstComponent=0;
		while (firstComponent<nComponents) {
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) readShortPeriodFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+(firstComponent+i)+"-reads-shortPeriod.txt"),IO.BUFFER_SIZE);
			readShortPeriodFile = new BufferedReader(new FileReader(readShortPeriodFilePath),IO.BUFFER_SIZE);
			str1=readShortPeriodFile.readLine();
			if (Reads.readIDsAreCompact) {
				while (str1!=null) {
					read=Integer.parseInt(str1.substring(0,str1.indexOf(",")));
					last=lastRead2component[read-Reads.firstRead];
					for (i=0; i<=last; i++) {
						component=read2component[read-Reads.firstRead][i];
						if (component>=firstComponent && component<firstComponent+nOutputFiles) readShortPeriodFiles[component].write(str1+"\n");
					}
					str1=readShortPeriodFile.readLine();
					if (read%10000==0) System.err.println("Processed read "+read);
				}
			}
			else {
				r=0;
				while (str1!=null && r<nReads) {
					read=Integer.parseInt(str1.substring(0,str1.indexOf(",")));
					if (Reads.readIDs[r]<read) {
						r++;
						continue;
					}
					else if (Reads.readIDs[r]>read) {
						str1=readShortPeriodFile.readLine();
						continue;
					}
					last=lastRead2component[r];
					for (i=0; i<=last; i++) {
						component=read2component[r][i];
						if (component>=firstComponent && component<firstComponent+nOutputFiles) readShortPeriodFiles[component].write(str1+"\n");
					}
					str1=readShortPeriodFile.readLine();
					if (read%10000==0) System.err.println("Processed read "+read);
				}
			}
			readShortPeriodFile.close();
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				readShortPeriodFiles[firstComponent+i].close();
				readShortPeriodFiles[firstComponent+i]=null;
			}
			firstComponent+=nOutputFiles;
		}
		
		// Filtering alignments
		alignmentFiles = new BufferedWriter[nComponents];
		shortPeriodFiles = new BufferedWriter[nComponents];
		nOutputFiles=IO.MAX_OPEN_FILES-2; firstComponent=0;
		while (firstComponent<nComponents) {
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				alignmentFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+(firstComponent+i)+"-alignments.txt"),IO.BUFFER_SIZE);
				shortPeriodFiles[firstComponent+i] = new BufferedWriter(new FileWriter(outputDir+"/"+(firstComponent+i)+"-alignments-shortPeriod.txt"),IO.BUFFER_SIZE);
				alignmentFiles[firstComponent+i].write(".\n");  // Fake header 
				alignmentFiles[firstComponent+i].write(".\n");
			}
			alignmentFile = new BufferedReader(new FileReader(alignmentFilePath),IO.BUFFER_SIZE);
			shortPeriodFile = new BufferedReader(new FileReader(shortPeriodFilePath),IO.BUFFER_SIZE);
			alignmentFile.readLine();
			alignmentFile.readLine();
			str1=alignmentFile.readLine(); str2=shortPeriodFile.readLine();
			r=0;
			while (str1!=null) {
				Alignments.readAlignmentFile(str1);
				readA=Reads.indexOfRead(Alignments.readA-1);
				readB=Reads.indexOfRead(Alignments.readB-1);
				last=Math.setIntersection(read2component[readA],0,lastRead2component[readA],read2component[readB],0,lastRead2component[readB],tmpArray,0);
				for (i=0; i<=last; i++) {
					component=tmpArray[i];
					if (component>=firstComponent && component<firstComponent+nOutputFiles) {
						alignmentFiles[component].write(str1+"\n");
						shortPeriodFiles[component].write(str2+"\n");
					}
				}
				str1=alignmentFile.readLine(); str2=shortPeriodFile.readLine();
				r++;
				if (r%1000000==0) System.err.println("Processed alignment "+r);
			}
			alignmentFile.close(); shortPeriodFile.close();
			for (i=0; i<nOutputFiles && firstComponent+i<nComponents; i++) {
				alignmentFiles[firstComponent+i].close();
				alignmentFiles[firstComponent+i]=null;
				shortPeriodFiles[firstComponent+i].close();
				shortPeriodFiles[firstComponent+i]=null;
			}
			firstComponent+=nOutputFiles;
		}
		lastRead2component=null;
		for (i=0; i<nReads; i++) read2component[i]=null;
		read2component=null;
	}
	
	
	/**
	 * Stores in $outputPath,shortPeriodOutputPath$ just the lines of $inputPath,
	 * shortPeriodInputPath$ where both reads of the alignment contain at least one node 
	 * of the graph.
	 *
	 * @param lastRead the procedure assumes that $getDistinctReads()$ has already been
	 * executed, and that this is its output.
	 */
	public static final void filterLA(String inputPath, String outputPath, String shortPeriodInputPath, String shortPeriodOutputPath, int lastRead) throws IOException {
		int i, j;
		int readA, readB;
		String str1, str2;
		BufferedReader inputFile, shortPeriodInputFile;
		BufferedWriter outputFile, shortPeriodOutputFile;
		
		inputFile = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
		shortPeriodInputFile = new BufferedReader(new FileReader(shortPeriodInputPath),IO.BUFFER_SIZE);
		outputFile = new BufferedWriter(new FileWriter(outputPath),IO.BUFFER_SIZE);
		shortPeriodOutputFile = new BufferedWriter(new FileWriter(shortPeriodOutputPath),IO.BUFFER_SIZE);
		str1=inputFile.readLine(); str1=inputFile.readLine();  // Input header
		outputFile.write(".\n"); outputFile.write(".\n");  // Output header
		str1=inputFile.readLine(); str2=shortPeriodInputFile.readLine();
		i=0;  // $nNeighbors$ pointer
		j=0;  // $nodesArray$ pointer: first interval of a read.
		while (str1!=null && i<=lastRead) {
			Alignments.readAlignmentFile(str1);
			readA=Alignments.readA-1;  // Read IDs in LAshow start from one, but in all our code they consistently start from zero.
			while (i<=lastRead && nNeighbors[i]<readA) i++;
			if (i>lastRead) break;
			if (nNeighbors[i]>readA) {
				str1=inputFile.readLine(); str2=shortPeriodInputFile.readLine();
				continue;
			}
			readB=Alignments.readB-1;
			if (Arrays.binarySearch(nNeighbors,0,lastRead+1,readB)>=0) {
				outputFile.write(str1+"\n");
				shortPeriodOutputFile.write(str2+"\n");
			}
			str1=inputFile.readLine(); str2=shortPeriodInputFile.readLine();
		}
		inputFile.close(); inputFile=null; 
		shortPeriodInputFile.close(); shortPeriodInputFile=null; 
		outputFile.close(); outputFile=null;
		shortPeriodOutputFile.close(); shortPeriodOutputFile=null;
	}
	
	
	/**
	 * Stores in $outputPath$ just the lines of $inputPath$ that refer to a read that 
	 * contains at least one node of the graph.
	 *
	 * @param lastRead the procedure assumes that $getDistinctReads()$ has already been
	 * executed, and that this is its output.
	 */
	public static final void filterShortPeriodTracks(String inputPath, String outputPath, int lastRead) throws IOException {
		int i;
		int read, trackRead;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		br = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
		bw = new BufferedWriter(new FileWriter(outputPath),IO.BUFFER_SIZE);
		str=br.readLine();
		if (str==null) {
			bw.close(); br.close();
			return;
		}
		trackRead=Integer.parseInt(str.substring(0,str.indexOf(",")));
		for (i=0; i<=lastRead; i++) {
			read=IntervalGraph.nNeighbors[i];
			if (trackRead<read) {
				str=br.readLine();
				if (str!=null) trackRead=Integer.parseInt(str.substring(0,str.indexOf(",")));
			}
			else if (trackRead==read) {
				bw.write(str);
				str=br.readLine();
				if (str!=null) trackRead=Integer.parseInt(str.substring(0,str.indexOf(",")));
			}
		}
		bw.close(); br.close();
	}
	
	
	/**
	 * Stores in $outputPath$ just the lines of $inputPath$ that correspond to reads that
	 * contain at least one node of the graph.
	 *
	 * @param lastRead the procedure assumes that $getDistinctReads()$ has already been
	 * executed, and that this is its output.
	 */
	public static final void filterQualities(String inputPath, String outputPath, int lastRead) throws IOException {
		int i, j;
		String str;
		BufferedReader inputFile;
		BufferedWriter outputFile;
		
		inputFile = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
		outputFile = new BufferedWriter(new FileWriter(outputPath),IO.BUFFER_SIZE);
		i=0; str=inputFile.readLine(); j=0;
		while (str!=null && i<=lastRead) {
			if (j==nNeighbors[i]) {
				outputFile.write(str+"\n");
				i++;
			}
			str=inputFile.readLine(); j++;
		}
		inputFile.close(); inputFile=null;
		outputFile.close(); outputFile=null;
	}
	
	
	/**
	 * Stores in $outputPath$ just the lines of $inputPath$ that correspond to reads that
	 * contain at least one node of the graph.
	 *
	 * @param lastRead the procedure assumes that $getDistinctReads()$ has already been
	 * executed, and that this is its output.
	 */
	public static final void filterReadLengths(String inputPath, String outputPath, int lastRead) throws IOException {
		int i, j;
		String str;
		BufferedReader inputFile;
		BufferedWriter outputFile;
		
		inputFile = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
		outputFile = new BufferedWriter(new FileWriter(outputPath),IO.BUFFER_SIZE);
		i=0; str=inputFile.readLine(); j=0;
		while (str!=null && i<=lastRead) {
			if (j==nNeighbors[i]) {
				outputFile.write(str+"\n");
				i++;
			}
			str=inputFile.readLine(); j++;
		}
		inputFile.close(); inputFile=null;
		outputFile.close(); outputFile=null;
	}
	
	
	/**
	 * Stores in $shortPeriodPath$ a bitvector that marks with a one every alignment of 
	 * $inputPath$ that is contained in the union of a set of overlapping short-period 
	 * intervals in readA.
	 *
	 * Remark: this procedure could be merged with $buildShortPeriodTrack()$.
	 */
	public static final void buildShortPeriodBitvector(String inputPath, String shortPeriodPath) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i;
		int readA;
		String str;
		BufferedReader inputFile;
		BufferedWriter shortPeriodFile;
	
		inputFile = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
		shortPeriodFile = new BufferedWriter(new FileWriter(shortPeriodPath),IO.BUFFER_SIZE);
		inputFile.readLine(); inputFile.readLine();  // Input header
		str=inputFile.readLine(); 
		i=0;  // $nodesArray$ pointer: first interval of a read.
		while (str!=null && i<nNodes) {
			Alignments.readAlignmentFile(str);
			readA=Alignments.readA-1;  // Read IDs in LAshow start from one, but in all our code they consistently start from zero.
			while (i<nNodes && nodesArray[i].read<readA) i++;
			if (i<nNodes) shortPeriodFile.write(inShortPeriod_nodes(i,readA,Alignments.startA,Alignments.endA,IDENTITY_THRESHOLD)?"1\n":"0\n");
			str=inputFile.readLine();
		}
		inputFile.close(); inputFile=null; shortPeriodFile.close(); shortPeriodFile=null;
	}
	
	
	/**
	 * Remark: the procedure assumes that $nodesArray$ is sorted by $read,start$.
	 *
	 * @param j smallest index of an element of $nodesArray$ that occurs in a read with ID
	 * $>=read$;
	 * @param threshold to decide identity of positions;
	 * @return TRUE iff $[start..end]$ belongs to the union of a set of overlapping or 
	 * adjacent short-period interval graph nodes of $read$.
	 */
	private static final boolean inShortPeriod_nodes(int j, int read, int start, int end, int threshold) {
		boolean found;
		int k;
		int windowStart, windowEnd;
	
		found=false;
		k=j; windowStart=-1; windowEnd=-1;
		if (nodesArray[k].read==read) {
			while (k<nNodes) {
				if (nodesArray[k].read>read) break;
				if (nodesArray[k].type!=Constants.INTERVAL_PERIODIC || nodesArray[k].hasLongPeriod) {
					k++;
					continue;
				}
				if (windowStart==-1) {
					windowStart=nodesArray[k].start; windowEnd=nodesArray[k].end;
					k++;
					continue;
				}
				if (nodesArray[k].start>windowEnd+threshold) {
					if (Intervals.isApproximatelyContained(start,end,windowStart,windowEnd) || Intervals.areApproximatelyIdentical(start,end,windowStart,windowEnd)) {
						found=true;
						break;
					}
					windowStart=nodesArray[k].start; 
					windowEnd=nodesArray[k].end;
					k++;
					continue;
				}
				windowEnd=Math.max(windowEnd,nodesArray[k].end);
				k++;
			}
			if (windowStart!=-1 && (Intervals.isApproximatelyContained(start,end,windowStart,windowEnd) || Intervals.areApproximatelyIdentical(start,end,windowStart,windowEnd))) found=true;
		}
		return found;
	}

	
	/**
	 * Prints a (possibly empty) file, with every line in the following format:
	 *
	 * r,n,s_{i_1},e_{i_1},...,s_{i_n},e_{i_n}
	 *
	 * where $r$ is the ID of a read, $[s_x..e_x]$ is a maximal union of interval graph 
	 * nodes marked as short-period in $r$, $s_x$ are sorted in increasing order in each 
	 * line, and lines are sorted by $r$. Only reads with at least one short-period 
	 * interval are stored.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by read.
	 */
	public static final void buildShortPeriodTrack(String outputPath) throws IOException {
		final int CAPACITY = 100;  // Arbitrary
		int i, j, k, h;
		int currentRead, lastWindow;
		Node currentNode;
		BufferedWriter bw;
		ShortPeriodWindow tmpWindow;
		ShortPeriodWindow[] windows;
		
		windows = new ShortPeriodWindow[CAPACITY];
		for (i=0; i<windows.length; i++) windows[i] = new ShortPeriodWindow();
		bw = new BufferedWriter(new FileWriter(outputPath),IO.BUFFER_SIZE);
		currentRead=-1; lastWindow=-1;
		for (i=0; i<nNodes; i++) {
			currentNode=nodesArray[i];
			if (currentNode.read!=currentRead) {
				if (lastWindow!=-1) {
					if (lastWindow>0) {
						Arrays.sort(windows,0,lastWindow+1);
						j=0;
						for (k=1; k<=lastWindow; k++) {
							if (windows[k].start>windows[j].end+IDENTITY_THRESHOLD) {
								j++;
								if (j==windows.length) {
									ShortPeriodWindow[] newWindows = new ShortPeriodWindow[windows.length+CAPACITY];
									System.arraycopy(windows,0,newWindows,0,windows.length);
									for (h=windows.length; h<newWindows.length; h++) newWindows[h] = new ShortPeriodWindow();
									windows=newWindows;
								}
								tmpWindow=windows[j];
								windows[j]=windows[k];
								windows[k]=tmpWindow;
							}
							else windows[j].end=Math.max(windows[j].end,windows[k].end);
						}
						lastWindow=j;
					}
					bw.write(currentRead+","+(lastWindow+1)+",");
					for (j=0; j<=lastWindow; j++) bw.write(windows[j].start+","+windows[j].end+",");
					bw.write("\n");
				}
				currentRead=currentNode.read;
				if (currentNode.type==Constants.INTERVAL_PERIODIC && !currentNode.hasLongPeriod) {
					lastWindow=0;
					windows[0].start=currentNode.start;
					windows[0].end=currentNode.end;
				}
				else lastWindow=-1;
			}
			else if (currentNode.type==Constants.INTERVAL_PERIODIC && !currentNode.hasLongPeriod) {
				lastWindow++;
				if (lastWindow==windows.length) {
					ShortPeriodWindow[] newWindows = new ShortPeriodWindow[windows.length+CAPACITY];
					System.arraycopy(windows,0,newWindows,0,windows.length);
					for (h=windows.length; h<newWindows.length; h++) newWindows[h] = new ShortPeriodWindow();
					windows=newWindows;
				}
				windows[lastWindow].start=currentNode.start;
				windows[lastWindow].end=currentNode.end;
			}
		}
		// Last read
		if (lastWindow!=-1) {
			if (lastWindow>0) {
				Arrays.sort(windows,0,lastWindow+1);
				j=0;
				for (k=1; k<=lastWindow; k++) {
					if (windows[k].start>windows[j].end+IDENTITY_THRESHOLD) {
						j++;
						if (j==windows.length) {
							ShortPeriodWindow[] newWindows = new ShortPeriodWindow[windows.length+CAPACITY];
							System.arraycopy(windows,0,newWindows,0,windows.length);
							for (h=windows.length; h<newWindows.length; h++) newWindows[h] = new ShortPeriodWindow();
							windows=newWindows;
						}
						tmpWindow=windows[j];
						windows[j]=windows[k];
						windows[k]=tmpWindow;
					}
					else windows[j].end=Math.max(windows[j].end,windows[k].end);
				}
				lastWindow=j;
			}
			bw.write(currentRead+","+(lastWindow+1)+",");
			for (j=0; j<=lastWindow; j++) bw.write(windows[j].start+","+windows[j].end+",");
			bw.write("\n"); 
		}
		bw.close();
		for (i=0; i<windows.length; i++) windows[i]=null;
		windows=null;
	}
	
	
	/**
	 * @return sorted by $r,s_i$; possibly of length zero.
	 */
	public static final ShortPeriodWindow[] loadShortPeriodTrack(String inputPath) throws IOException {
		int i, j, p, q;
		int read, nWindows;
		String str;
		BufferedReader br;
		ShortPeriodWindow[] out;
		
		// Counting number of windows
		br = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
		nWindows=0;
		str=br.readLine();
		while (str!=null) {
			p=str.indexOf(","); q=str.indexOf(",",p+1);
			nWindows+=Integer.parseInt(str.substring(p+1,q));
			str=br.readLine();
		}
		br.close();
		
		// Loading windows
		out = new ShortPeriodWindow[nWindows];
		if (nWindows==0) return out;
		for (i=0; i<nWindows; i++) out[i] = new ShortPeriodWindow();
		br = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
		str=br.readLine(); j=-1;
		while (str!=null) {
			p=str.indexOf(","); q=str.indexOf(",",p+1);
			read=Integer.parseInt(str.substring(0,p));
			nWindows=Integer.parseInt(str.substring(p+1,q));
			for (i=0; i<nWindows; i++) {
				j++;
				out[j].read=read;
				p=q+1; q=str.indexOf(",",p+1);
				out[j].start=Integer.parseInt(str.substring(p,q));
				p=q+1; q=str.indexOf(",",p+1);
				out[j].end=Integer.parseInt(str.substring(p,q));
			}
			str=br.readLine();
		}
		br.close();
		return out;
	}
	
	
	/**
	 * @return for every element of $Reads.readIDs$, the position of the first element in 
	 * $shortPeriodWindows$ with the same read ID, or -1 if no such element exists.
	 */
	public static final int[] buildTrackPointers(IntervalGraph.ShortPeriodWindow[] shortPeriodWindows) {
		int i, j;
		final int nWindows = shortPeriodWindows.length;
		int[] out = new int[Reads.nReads];
	
		Math.set(out,Reads.nReads-1,-1);
		if (nWindows==0) return out;
		i=0; j=0;
		while (i<nWindows && j<Reads.nReads) {
			if (shortPeriodWindows[i].read<Reads.readIDs[j]) {
				i++;
				continue;
			}
			if (shortPeriodWindows[i].read>Reads.readIDs[j]) {
				j++;
				continue;
			}
			if (out[j]==-1) out[j]=i;
			i++;
		}
		return out;
	}
	
	
	public static class ShortPeriodWindow implements Comparable {
		public int read, start, end;
		
		public int compareTo(Object other) {
			ShortPeriodWindow otherWindow = (ShortPeriodWindow)other;
			if (read<otherWindow.read) return -1;
			else if (read>otherWindow.read) return 1;
			if (start<otherWindow.start) return -1;
			else if (start>otherWindow.start) return 1;
			return 0;
		}
	}
	
	
	/**
	 * Remark: the procedure performs a simple linear scan over $windows$, since the 
	 * number of tracks per read is likely small.
	 *
	 * @param from position in $windows$ from which the search starts;
	 * @param trackPointers output of $buildTrackPointers()$;
	 * @return TRUE iff $[start..end]$ is approximately contained in an element of 
	 * $windows$.
	 */
	public static final boolean inShortPeriodTrack(int read, int start, int end, ShortPeriodWindow[] windows, int[] trackPointers) {
		int i = trackPointers[Reads.indexOfRead(read)];
		if (i==-1) return false;
		final int nWindows = windows.length;
		
		while (i<nWindows) {
			if (windows[i].read!=read) break;
			if (Intervals.isApproximatelyContained(start,end,windows[i].start,windows[i].end) || Intervals.areApproximatelyIdentical(start,end,windows[i].start,windows[i].end)) return true;
			i++;
		}
		return false;
	}
	
	
	/**
	 * Remark: the procedure does not distinguish between left- and right-maximal, since
	 * reads are sampled from both strands and we would need to put them in a common
	 * orientation to distinguish between the two sides. Distinguishing between the two 
	 * sides is not likely to give any reliable information (not even for periodic
	 * substrings), because of long random insertions.
	 *
	 * @param components computes just the size of component IDs in $[0..nComponents-1]$; 
	 * if NULL, all components are used;
	 * @param strict counts only nodes with exactly one component tag;
	 * @param discardDisconnected a node that has no ON neighbor in component C does not 
	 * contribute to the count of component C;
	 * @param size the number of nodes in each component;
	 * @param nMaximal the number of maximal ends of intervals in each component.
	 */
	public static final void getComponentStats(int[] components, int nComponents, boolean strict, boolean discardDisconnected, int[] size, int[] nMaximal) {
		boolean found;
		int i, j, k;
		int component, neighbor;
		Node node;
		
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			if (strict && node.lastComponent!=0) continue;
			for (j=0; j<=node.lastComponent; j++) {
				component=node.components[j];
				if (component==-1) continue;
				if (discardDisconnected) {
					found=false;
					for (k=0; k<nNeighbors[i]; k++) {
						neighbor=neighbors[i][k].getTo(i);
						if ( neighbors[i][k].on && 
							 (strict?nodesArray[neighbor].lastComponent==0:true) &&
						     (components==null?true:Arrays.binarySearch(nodesArray[neighbor].components,0,nodesArray[neighbor].lastComponent+1,component)>=0)
						   ) {
							found=true;
							break;
						}
					}
					if (!found) continue;
				}
				k=components==null?component:Arrays.binarySearch(components,0,nComponents,component);
				if (k>=0) {
					size[k]++;
					if (node.isLeftMaximal) nMaximal[k]++;
					if (node.isRightMaximal) nMaximal[k]++;
				}
			}
		}
	}
	
	
	/**
	 * @return the first element of the array contains the number of distinct component
	 * tags of all nodes in the graph; the following elements contain the component tags 
	 * themselves.
	 */
	public static final int[] getComponents() {
		final int GROWTH_RATE = 100;
		int i, j, p, c;
		int lastComponent, last;
		int[] components, tmpComponents;
		
		components = new int[GROWTH_RATE];
		lastComponent=0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			last=nodesArray[i].lastComponent;
			for (j=0; j<=last; j++) {
				c=nodesArray[i].components[j];
				p=Math.linearSearch_unsorted(components,1,lastComponent+1,c,false);
				if (p<0) {
					lastComponent++;
					if (lastComponent==components.length) {
						tmpComponents = new int[components.length+GROWTH_RATE];
						System.arraycopy(components,0,tmpComponents,0,components.length);
						components=tmpComponents;
					}
					components[lastComponent]=c;
				}
			}
		}
		components[0]=lastComponent;
		return components;
	}
	
	
	/**
	 * Sets the variable $avgDiffs$ of every edge to an average, assuming it is a sum.
	 * This procedure should be called only when needed, and the status of $avgDiffs$ as
	 * a sum should be preserved throughout the code.
	 *
	 * Remark: the procedure does not assumes $nodesArray$ to be sorted, and it does not 
	 * create new intervals.
	 */
	public static final void setAvgDiffs() {
		int i, j;
		Edge edge;
		
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (edge.avgDiffs==-1) continue;
				if (edge.nAlignmentsForward>=0 && edge.nAlignmentsBackward>=0) {  // The sign is used to mark a processed edge, since each edge occurs twice in the matrix.
					edge.avgDiffs/=edge.getNAlignments();
					edge.nAlignmentsForward=-edge.nAlignmentsForward;
					edge.nAlignmentsBackward=-edge.nAlignmentsBackward;
				}
				else {
					edge.nAlignmentsForward=-edge.nAlignmentsForward;
					edge.nAlignmentsBackward=-edge.nAlignmentsBackward;
				}				
			}
		}
	}
	
	
	/**
	 * Sets $contains$, $isContained$, $hasRepeat$, $isRCPalindrome$ for every node.
	 *
	 * Remark: the procedure does not assume $nodesArray$ to be sorted, and it does not 
	 * create new intervals.
	 */
	private static final void finalizeEdgeAndNodeProperties() {
		boolean periodicI, periodicJ;
		int i, j;
		int nRepetitiveIntervals, nPalindromicIntervals;
		int nContainingIntervals, nContainedIntervals;
		int nFlexibleEdgesP_P, nFlexibleEdgesP_NP, nFlexibleEdgesNP_NP;
		IntervalGraph.Node node, nodeTo;
		IntervalGraph.Edge edge;
		
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			node.isContained=false; node.contains=false;
			node.hasRepeat=false; node.isRCPalindrome=false;
		}
		nRepetitiveIntervals=0; nPalindromicIntervals=0;
		nContainingIntervals=0; nContainedIntervals=0;
		nFlexibleEdgesP_P=0; nFlexibleEdgesP_NP=0; nFlexibleEdgesNP_NP=0;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			periodicI=node.type==Constants.INTERVAL_PERIODIC;
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				nodeTo=nodesArray[edge.getTo(i)];
				periodicJ=nodeTo.type==Constants.INTERVAL_PERIODIC;
				if (!edge.supplement && edge.flexibleOverlap!=0 && edge.flexibleOverlap!=-1) {
					if (periodicI) {
						if (periodicJ) nFlexibleEdgesP_P++;
						else nFlexibleEdgesP_NP++;
					}
					else {
						if (periodicJ) nFlexibleEdgesP_NP++;
						else nFlexibleEdgesNP_NP++;
					}
				}
				
				// Detecting whether an interval has (possibly reverse-complemented)
				// internal repeats.
				if (node.type!=Constants.INTERVAL_PERIODIC && nodeTo.type!=Constants.INTERVAL_PERIODIC) {
					// Different types
					if (edge.overlap!=-1) {
						if ( edge.containment==Constants.CONTAINMENT_ONE_IN_TWO || 
							 edge.insertion==Constants.INSERTION_ONE_IN_TWO ||
							 edge.containment==Constants.CONTAINMENT_IDENTICAL
						   ) {
						   if (!nodesArray[edge.nodeID2].hasRepeat) {
							   nodesArray[edge.nodeID2].hasRepeat=true;
							   nRepetitiveIntervals++;
						   }
						}
						if ( edge.containment==Constants.CONTAINMENT_TWO_IN_ONE || 
							 edge.insertion==Constants.INSERTION_TWO_IN_ONE ||
							 edge.containment==Constants.CONTAINMENT_IDENTICAL
						   ) {
						   if (!nodesArray[edge.nodeID1].hasRepeat) {
							   nodesArray[edge.nodeID1].hasRepeat=true;
							   nRepetitiveIntervals++;
						   }
						}
						if ( ((edge.overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 && (edge.overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) ||
							 ((edge.overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && (edge.overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0)
						   ) {
						   if (!nodesArray[edge.nodeID1].hasRepeat) {
							   nodesArray[edge.nodeID1].hasRepeat=true;
							   nRepetitiveIntervals++;
						   }
						}
						if ( ((edge.overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 && (edge.overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) ||
							 ((edge.overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && (edge.overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0)
						   ) {
						   if (!nodesArray[edge.nodeID2].hasRepeat) {
							   nodesArray[edge.nodeID2].hasRepeat=true;
							   nRepetitiveIntervals++;
						   }
						}
						// We could also handle the combination overlap+sharedSubstring,
						// but which node has an internal repeat depends on whether the
						// shared substring involves a substring inside an overlapping
						// region and a substring outside (just one node), or two
						// substrings inside the overlapping regions (both nodes).
					    // We skip this for the moment.
						//
						// We could detect borders by checking $overlap==0$. 
						// We skip this for the moment.
					}
					if (edge.insertion!=-1) {
						if (edge.insertion==Constants.INSERTION_TWO_IN_ONE && (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE || edge.sharedSubstring!=-1)) {
							if (!nodesArray[edge.nodeID1].hasRepeat) {
								nodesArray[edge.nodeID1].hasRepeat=true;
								nRepetitiveIntervals++;
							}
						}
						if (edge.insertion==Constants.INSERTION_ONE_IN_TWO && (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO || edge.sharedSubstring!=-1)) {
							if (!nodesArray[edge.nodeID2].hasRepeat) {
								nodesArray[edge.nodeID2].hasRepeat=true;
								nRepetitiveIntervals++;
							}
						}
					}
					if (edge.sharedSubstring!=-1) {
						if (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO || edge.containment==Constants.CONTAINMENT_IDENTICAL) {
							if (!nodesArray[edge.nodeID2].hasRepeat) {
								nodesArray[edge.nodeID2].hasRepeat=true;
								nRepetitiveIntervals++;
							}
						}
						if (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE || edge.containment==Constants.CONTAINMENT_IDENTICAL) {
							if (!nodesArray[edge.nodeID1].hasRepeat) {
								nodesArray[edge.nodeID1].hasRepeat=true;
								nRepetitiveIntervals++;
							}
						}
					}
					// Same type
					if (edge.containment!=-1 && edge.insertion==-1 && edge.overlap==-1 && edge.sharedSubstring==-1 && edge.getNAlignments()>1) {
						if (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE) {
							if (!nodesArray[edge.nodeID1].hasRepeat) {
								nodesArray[edge.nodeID1].hasRepeat=true;
								nRepetitiveIntervals++;
							}
						}
						else if (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO) {
							if (!nodesArray[edge.nodeID2].hasRepeat) {
								nodesArray[edge.nodeID2].hasRepeat=true;
								nRepetitiveIntervals++;
							}
						}
					}
					if (edge.containment==-1 && edge.insertion!=-1 && edge.overlap==-1 && edge.sharedSubstring==-1 && edge.getNAlignments()>1) {
						if (edge.insertion==Constants.INSERTION_TWO_IN_ONE) {
							if (!nodesArray[edge.nodeID1].hasRepeat) {
								nodesArray[edge.nodeID1].hasRepeat=true;
								nRepetitiveIntervals++;
							}
						}
						else if (edge.insertion==Constants.INSERTION_ONE_IN_TWO) {
							if (!nodesArray[edge.nodeID2].hasRepeat) {
								nodesArray[edge.nodeID2].hasRepeat=true;
								nRepetitiveIntervals++;
							}
						}
					}
					if (edge.containment==-1 && edge.insertion==-1 && edge.overlap!=-1 && Math.popcount(edge.overlap)==1 && edge.sharedSubstring==-1 && edge.getNAlignments()>1) {
						if (!nodesArray[edge.nodeID1].hasRepeat) {
							nodesArray[edge.nodeID1].hasRepeat=true;
							nRepetitiveIntervals++;
						}
						if (!nodesArray[edge.nodeID2].hasRepeat) {
							nodesArray[edge.nodeID2].hasRepeat=true;
							nRepetitiveIntervals++;
						}
					}
				}
				
				// Detecting entirely palindromic intervals
				if (edge.containment==Constants.CONTAINMENT_IDENTICAL && edge.orientation==2) {
					if (!nodesArray[edge.nodeID1].isRCPalindrome) {
						nodesArray[edge.nodeID1].isRCPalindrome=true;
						nPalindromicIntervals++;
					}
					if (!nodesArray[edge.nodeID2].isRCPalindrome) {
						nodesArray[edge.nodeID2].isRCPalindrome=true;
						nPalindromicIntervals++;
					}
				}
				
				// Marking contained/containing intervals
				if (!edge.supplement) {
					if (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO) {
						if (!nodesArray[edge.nodeID1].isContained) {
							nodesArray[edge.nodeID1].isContained=true;
							nContainedIntervals++;
						}
						if (!nodesArray[edge.nodeID2].contains) {
							nodesArray[edge.nodeID2].contains=true;
							nContainingIntervals++;
						}
					}
					else if (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE) {
						if (!nodesArray[edge.nodeID2].isContained) {
							nodesArray[edge.nodeID2].isContained=true;
							nContainedIntervals++;
						}
						if (!nodesArray[edge.nodeID1].contains) {
							nodesArray[edge.nodeID1].contains=true;
							nContainingIntervals++;
						}
					}
				}
			}
		}
		System.out.println("Intervals with internal repeats: "+nRepetitiveIntervals+" ("+IO.getPercent(nRepetitiveIntervals,nNodes)+"%)");
		System.out.println("Palindromic intervals: "+nPalindromicIntervals+" ("+IO.getPercent(nPalindromicIntervals,nNodes)+"%)");
		System.out.println("Intervals that contain other intervals: "+nContainingIntervals+" ("+IO.getPercent(nContainingIntervals,nNodes)+"%)");
		System.out.println("Intervals that are contained in other intervals: "+nContainedIntervals+" ("+IO.getPercent(nContainedIntervals,nNodes)+"%)");
		System.out.println("Flexible edges: PP="+(nFlexibleEdgesP_P/2)+" P-NP="+(nFlexibleEdgesP_NP/2)+" NP-NP="+(nFlexibleEdgesNP_NP/2));
	}
	
	
	/**
	 * Sets the $isContained$ and $contains$ flags of each node.
	 *
	 * Remark: the procedure does not use supplement edges.
	 *
	 * @param onlyOn uses only ON edges.
	 */
	public static final void markContainedContainingIntervals(boolean onlyOn) {
		int i, j;
		int nContainingIntervals, nContainedIntervals;
		IntervalGraph.Node node;
		IntervalGraph.Edge edge;
		
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			node.isContained=false;
			node.contains=false;
		}
		nContainingIntervals=0; nContainedIntervals=0;
		for (i=0; i<nNodes; i++) {
			node=nodesArray[i];
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j];
				if (edge.supplement) continue;
				if (onlyOn && !edge.on) continue;
				if (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO) {
					if (!nodesArray[edge.nodeID1].isContained) {
						nodesArray[edge.nodeID1].isContained=true;
						nContainedIntervals++;
					}
					if (!nodesArray[edge.nodeID2].contains) {
						nodesArray[edge.nodeID2].contains=true;
						nContainingIntervals++;
					}
				}
				else if (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE) {
					if (!nodesArray[edge.nodeID2].isContained) {
						nodesArray[edge.nodeID2].isContained=true;
						nContainedIntervals++;
					}
					if (!nodesArray[edge.nodeID1].contains) {
						nodesArray[edge.nodeID1].contains=true;
						nContainingIntervals++;
					}
				}
			}			
		}
		System.out.println("Intervals that contain other intervals: "+nContainingIntervals+" ("+((100.0*nContainingIntervals)/nNodes)+"%)");
		System.out.println("Intervals that are contained in other intervals: "+nContainedIntervals+" ("+((100.0*nContainedIntervals)/nNodes)+"%)");
	}
	
	
	/**
	 * Flags in $marked$ the components in $components$ that contain at least one periodic
	 * node.
	 *
	 * @param components[0..lastComponent] sorted; if NULL, all components are used;
	 * @param strict if TRUE, discards nodes with more than one component;
	 * @param periodicType 0=any period; 1=short-period only; 2=long-period only.
	 */
	public static final void markPeriodicComponents(int[] components, int lastComponent, boolean[] marked, boolean strict, int periodicType) {
		int i, j, c;
		int last;

		Math.set(marked,lastComponent,false);
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].type!=Constants.INTERVAL_PERIODIC) continue;
			if ((nodesArray[i].hasLongPeriod && periodicType==1) || (!nodesArray[i].hasLongPeriod && periodicType==2)) continue;
			last=nodesArray[i].lastComponent;
			if (last==-1 || (last==0 && nodesArray[i].components[0]==-1) || (strict && last>0)) continue;
			for (j=0; j<=last; j++) {
				if (components==null) marked[nodesArray[i].components[j]]=true;
				else {
					c=Arrays.binarySearch(components,0,lastComponent+1,nodesArray[i].components[j]);
					if (c<0) {
						System.err.println("markPeriodicComponents> ERROR: component not found among those in input: "+nodesArray[i].components[j]);
						for (int x=0; x<=lastComponent; x++) System.err.print(components[x]+", ");
						System.err.println();
						System.exit(1);
					}
					marked[c]=true;
				}
			}
		}
	}
	
	
	/**
	 * Transforms into containments all overlaps whose overhangs are too short.
	 *
	 * Remark: if an edge corresponds to multiple overlap alignments, or if an edge is 
	 * already also of overlap type, it is transformed into the containment that 
	 * corresponds to a largest overlap.
	 *
	 * Remark: the procedure does not consider edges in which one of the nodes is periodic
	 * with short period, since in this case short overhangs can be informative.
	 *
	 * Remark: the procedure uses the $printed$ flag of $Edge$ to mark edges that have
	 * already been processed (but not necessarily modified).
	 *
	 * Remark: the procedure does not assumes $nodesArray$ to be sorted, and it does not 
	 * create new intervals.
	 */
	private static final void removeShortOverlaps() {
		final int THRESHOLD = IO.quantum<<1;  // Arbitrary
		boolean edgeChanged;
		int i, j, k;
		int to, changedEdges, overlap;
		long totalEdges;
		Edge edge;
		int[][] matrix;
		
		totalEdges=0;
		for (i=0; i<nNodes; i++) {
			for (j=0; j<nNeighbors[i]; j++) neighbors[i][j].printed=0;
			totalEdges+=nNeighbors[i];
		}
		totalEdges>>=1; changedEdges=0;
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].type==Constants.INTERVAL_PERIODIC && !nodesArray[i].hasLongPeriod) continue;
			for (j=0; j<nNeighbors[i]; j++) {
				edge=neighbors[i][j]; to=edge.getTo(i);				
				if (edge.printed==1 || !edge.isType_overlap() || (nodesArray[to].type==Constants.INTERVAL_PERIODIC && !nodesArray[to].hasLongPeriod)) continue;
				edge.printed=1;
				overlap=edge.overlap;
				if (edge.maxOverhangs!=null) matrix=edge.maxOverhangs;
				else matrix=edge.overhangs;
				edgeChanged=false;
				for (k=0; k<4; k++) {
					if ((overlap&Constants.id2overlap(k))==0) continue;
					if (IO.CONSISTENCY_CHECKS) {
						if (matrix[k][0]==-1 || matrix[k][1]==-1) {
							System.err.println("removeShortOverlaps> ERROR: the following edge has no overhangs for overlap "+Constants.id2overlap(k)+": "+matrix[k][0]+","+matrix[k][1]);
							System.err.println(edge);
							System.err.println("node1: "+nodesArray[edge.nodeID1]);
							System.err.println("node2: "+nodesArray[edge.nodeID2]);
							System.exit(1);
						}
					}
					if (matrix[k][0]>THRESHOLD && matrix[k][1]>THRESHOLD) continue;
					edgeChanged=true;
					if (matrix[k][0]<=THRESHOLD) {
						if (matrix[k][1]<=THRESHOLD) {
							edge.setType_noContainment();
							edge.setType(Constants.CONTAINMENT_IDENTICAL,true);
						}
						else {
							if (!edge.isType_containment()) {	
								edge.setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
								if (k==0 || k==2) {
									edge.containmentOverhangLeft=0; 
									edge.containmentOverhangRight=edge.overhangs[k][1];
								}
								else {
									edge.containmentOverhangLeft=edge.overhangs[k][1];
									edge.containmentOverhangRight=0;
								}
							}
							else if (edge.isType(Constants.CONTAINMENT_ONE_IN_TWO)) {
								if (-edge.containmentOverhangLeft-edge.containmentOverhangRight<-edge.overhangs[k][1]) {
									if (k==0 || k==2) {
										edge.containmentOverhangLeft=0; 
										edge.containmentOverhangRight=edge.overhangs[k][1];
									}
									else {
										edge.containmentOverhangLeft=edge.overhangs[k][1];
										edge.containmentOverhangRight=0;
									}
								}
							}
							else if (edge.isType(Constants.CONTAINMENT_TWO_IN_ONE)) edge.containment=Constants.CONTAINMENT_IDENTICAL;
						}
					}
					else {
						if (!edge.isType_containment()) {	
							edge.setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
							if (k<=1) {
								edge.containmentOverhangLeft=0; 
								edge.containmentOverhangRight=edge.overhangs[k][0];
							}
							else {
								edge.containmentOverhangLeft=edge.overhangs[k][0];
								edge.containmentOverhangRight=0;
							}
						}
						else if (edge.isType(Constants.CONTAINMENT_ONE_IN_TWO)) edge.containment=Constants.CONTAINMENT_IDENTICAL;
						else if (edge.isType(Constants.CONTAINMENT_TWO_IN_ONE)) {
							if (-edge.containmentOverhangLeft-edge.containmentOverhangRight<-edge.overhangs[k][0]) {
								if (k<=1) {
									edge.containmentOverhangLeft=0; 
									edge.containmentOverhangRight=edge.overhangs[k][0];
								}
								else {
									edge.containmentOverhangLeft=edge.overhangs[k][0];
									edge.containmentOverhangRight=0;
								}
							}
						}
					}
					edge.removeOverlap(Constants.id2overlap(k));
				}
				if (edgeChanged) {
					changedEdges++;
					if (edge.overlap==0) edge.setType_noOverlap();
					// Disallowing identity and insertion (copied from $Edge.addTypes()$).
					if (edge.isType(Constants.CONTAINMENT_IDENTICAL) && edge.isType_insertion()) edge.setType_noInsertion();
				}
			}
		}
		System.err.println(changedEdges+" edges transformed from overlap to containment ("+IO.getPercent(changedEdges,totalEdges)+"%)");
	}

	
	
	
	// ------------------------ EDGES INDUCED BY ALL ALIGNMENTS  -------------------------
	
	/**
	 * Uses an alignment that is approximately identical to, or that approximately 
	 * contains, an interval in readA and an interval in readB, to create an identity edge
	 * of supplement type between such intervals. This is useful, since the only alignment
	 * that proves identity between two intervals might be assigned to a different 
	 * interval in one of the reads (for example if a substring interval is contained in a
	 * suffix interval, and the alignment falls between the two). The procedure considers 
	 * all pairs of intervals that are approximately identical to, or that are 
	 * approximately contained in, the interval of the alignment in readA and readB.
	 *
	 * Remark: for speed, the procedure considers just pairs of intervals that can belong 
	 * to kernels in Step 3.
	 *
	 * Remark: one might think of using a similar all-pairs procedure to build all the 
	 * edges of the interval graph, rather than using just the assignments of alignments 
	 * to intervals in the interval-connection file. This might make clusters harder to 
	 * detect, since e.g. it might connect, with a containment edge, a fragment of a 
	 * repeat in readA, to fragments of all repeats that contain it in readB, rather than 
	 * just to the fragment of the same repeat in readB. It would also make the graphs 
	 * significantly larger, slowing down all the procedures downstream.
	 *
	 * Remark: the procedure assumes $nodesArray$ to be sorted by $read,start,nodeID$, and
	 * it does not create new intervals. The procedure does not assume $neighbors[i]$ to 
	 * be sorted for any $i$, and it considers all edges in $neighbors[i]$.
	 *
	 * @param kernelsComputed if TRUE, assumes that kernel tags have already been assigned
	 * to nodes;
	 * @param alignmentsFilePath assumed to contain all and only the alignments between 
	 * reads such that at least one node of the graph occurs in each read; the same 
	 * alignment can occur twice (once for each read considered as readA);
	 * @return the number of edges added by the procedure.
	 */
	public static final int addIdentityEdges(boolean kernelsComputed, String alignmentsFilePath) throws IOException {
		final int THRESHOLD = IO.quantum*3;  // Arbitrary
		final int PREVIOUS_ORDER = Node.order;
		int i, j, p, q;
		int readA, readB, nAdded, nAlignments, nEdges, diffs;
		String str;
		Node tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		AlignmentPair tmpAlignmentPair = new AlignmentPair();
		BufferedReader alignmentsFile;
		
		System.err.println("Adding supplement edges of identity type...");		
		Node.order=Node.READ_START_NODEID;
		nEdges=getNEdges();
		alignmentsFile = new BufferedReader(new FileReader(alignmentsFilePath),IO.BUFFER_SIZE);
		str=alignmentsFile.readLine();
		str=alignmentsFile.readLine();  // Skipping the first two lines
		str=alignmentsFile.readLine();
		nAdded=0; nAlignments=0;
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			nAlignments++;
			readA=Alignments.readA-1;
			tmpNode.read=readA; tmpNode.start=Alignments.startA; tmpNode.id=-1;
			p=Arrays.binarySearch(nodesArray,0,nNodes,tmpNode);
			if (p<0) {
				p=-1-p;
				p=Math.min(p,nNodes-1);
			}
			for (i=p-1; i>=0; i--) {
				if (nodesArray[i].read!=readA || nodesArray[i].start<Alignments.startA-THRESHOLD) break;
				if (nodesArray[i].end<=Alignments.startA || nodesArray[i].end>Alignments.endA+THRESHOLD) continue;
				if (kernelsComputed?nodesArray[i].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[i])) continue;
				readB=Alignments.readB-1;
				tmpNode.read=readB; tmpNode.start=Alignments.startB; tmpNode.id=-1;
				q=Arrays.binarySearch(nodesArray,0,nNodes,tmpNode);
				if (q<0) {
					q=-1-q;
					q=Math.min(q,nNodes-1);
				}
				for (j=q-1; j>=0; j--) {
					if (nodesArray[j].read!=readB || nodesArray[j].start<Alignments.startB-THRESHOLD) break;
					if (j==i) continue;
					if (nodesArray[j].end<=Alignments.startB || nodesArray[j].end>Alignments.endB+THRESHOLD) continue;
					if (kernelsComputed?nodesArray[j].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[j])) continue;
					if (Intervals.areApproximatelyIdentical(nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation)) {
						diffs=Alignments.projectDiffs(Alignments.diffs,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end);
						if (addIdentityEdge(i,j,Alignments.orientation?0:1,diffs,tmpAlignmentPair)) nAdded++;
					}
				}
				for (j=q; j<nNodes; j++) {
					if (nodesArray[j].read!=readB || nodesArray[j].start>=Alignments.endB) break;
					if (j==i) continue;
					if (nodesArray[j].end>Alignments.endB+THRESHOLD) continue;
					if (kernelsComputed?nodesArray[j].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[j])) continue;
					if (Intervals.areApproximatelyIdentical(nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation)) {
						diffs=Alignments.projectDiffs(Alignments.diffs,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end);
						if (addIdentityEdge(i,j,Alignments.orientation?0:1,diffs,tmpAlignmentPair)) nAdded++;
					}
				}
			}
			for (i=p; i<nNodes; i++) {
				if (nodesArray[i].read!=readA || nodesArray[i].start>=Alignments.endA) break;
				if (nodesArray[i].end>Alignments.endA+THRESHOLD) continue;
				if (kernelsComputed?nodesArray[i].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[i])) continue;
				readB=Alignments.readB-1;
				tmpNode.read=readB; tmpNode.start=Alignments.startB; tmpNode.id=-1;
				q=Arrays.binarySearch(nodesArray,0,nNodes,tmpNode);
				if (q<0) {
					q=-1-q;
					q=Math.min(q,nNodes-1);
				}
				for (j=q-1; j>=0; j--) {
					if (nodesArray[j].read!=readB || nodesArray[j].start<Alignments.startB-THRESHOLD) break;
					if (j==i) continue;
					if (nodesArray[j].end<=Alignments.startB || nodesArray[j].end>Alignments.endB+THRESHOLD) continue;
					if (kernelsComputed?nodesArray[j].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[j])) continue;
					if (Intervals.areApproximatelyIdentical(nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation)) {
						diffs=Alignments.projectDiffs(Alignments.diffs,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end);
						if (addIdentityEdge(i,j,Alignments.orientation?0:1,diffs,tmpAlignmentPair)) nAdded++;
					}
				}
				for (j=q; j<nNodes; j++) {
					if (nodesArray[j].read!=readB || nodesArray[j].start>=Alignments.endB) break;
					if (j==i) continue;
					if (nodesArray[j].end>Alignments.endB+THRESHOLD) continue;
					if (kernelsComputed?nodesArray[j].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[j])) continue;
					if (Intervals.areApproximatelyIdentical(nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation)) {
						diffs=Alignments.projectDiffs(Alignments.diffs,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end);
						if (addIdentityEdge(i,j,Alignments.orientation?0:1,diffs,tmpAlignmentPair)) nAdded++;
					}
				}
			}
			if (nAlignments%1000000==0) System.err.println("Processed alignment "+(nAlignments/1000000)+"M");
			str=alignmentsFile.readLine();
		}
		alignmentsFile.close();
		Node.order=PREVIOUS_ORDER;
		System.err.println("Added "+nAdded+" supplement edges of identity type ("+IO.getPercent(nAdded,nEdges)+"%)");
		return nAdded;
	}
	
	
	/**
	 * Remark: the procedure does not assume $neighbors[i]$ to be sorted, and it considers
	 * all edges in $neighbors[i]$.
	 *
	 * @param orientation 0=forward; 1=RC;
	 * @param diffs estimate of the number of diffs to be assigned to the edge (derived 
	 * from the diffs of some alignment); this is a number, not a rate;
	 * @param tmpAlignmentPair temporary space;
	 * @return TRUE iff a new edge is added by the procedure.
	 */
	private static final boolean addIdentityEdge(int n1, int n2, int orientation, int diffs, AlignmentPair tmpAlignmentPair) {
		final int RESIZE_UNIT = 10;  // Arbitrary
		boolean added;
		int i, n, m;
		double normalizedDiffs;
		Edge tmpEdge;
		Edge[] tmpNeighbors;
		
		if (nodesArray[n1].type==Constants.INTERVAL_PERIODIC) {
			if (nodesArray[n2].type!=Constants.INTERVAL_PERIODIC) return false;
			if (nodesArray[n2].hasLongPeriod!=nodesArray[n1].hasLongPeriod) return false;
		} 
		else if (nodesArray[n2].type==Constants.INTERVAL_PERIODIC) return false;
		
		// Checking whether the edge already exists
		normalizedDiffs=2.0*diffs/(nodesArray[n1].end-nodesArray[n1].start+nodesArray[n2].end-nodesArray[n2].start+2);
		n=nNeighbors[n1]<nNeighbors[n2]?n1:n2;
		m=n==n1?n2:n1;		
		for (i=0; i<nNeighbors[n]; i++) {
			tmpEdge=neighbors[n][i];
			if (tmpEdge.supplement) continue;
			if (tmpEdge.getTo(n)==m) {
				if (tmpEdge.containment==Constants.CONTAINMENT_IDENTICAL) {
					added=false;
					if (tmpEdge.orientation==0 && orientation==1) {
						tmpEdge.orientation=2;
						tmpEdge.nAlignmentsBackward++;
						added=true;
					}
					else if (tmpEdge.orientation==1 && orientation==0) {
						tmpEdge.nAlignmentsForward++;
						tmpEdge.orientation=2;
						added=true;
					}
					if (added) {
						if (tmpEdge.avgDiffs!=-1) tmpEdge.avgDiffs+=normalizedDiffs;
						return true;
					}
					else return false;
				}
				else if (tmpEdge.containment!=-1 || tmpEdge.insertion!=-1) return false;
			}
		}
		for (i=0; i<nNeighbors[n]; i++) {
			tmpEdge=neighbors[n][i];
			if (!tmpEdge.supplement) continue;
			if (tmpEdge.getTo(n)==m) {
				if (tmpEdge.containment==-1 && tmpEdge.insertion==-1) {
					tmpEdge.containment=Constants.CONTAINMENT_IDENTICAL;
					if (tmpEdge.orientation==-1) {
						tmpEdge.nAlignmentsForward++;						
					}
					else if (tmpEdge.orientation==0) {
						if (orientation==0) tmpEdge.nAlignmentsForward++;
						else {
							tmpEdge.nAlignmentsBackward++;
							tmpEdge.orientation=2;
						}
					}
					else if (tmpEdge.orientation==1) {
						if (orientation==0) {
							tmpEdge.nAlignmentsForward++;
							tmpEdge.orientation=2;
						}
						else tmpEdge.nAlignmentsBackward++;
					}
					else if (tmpEdge.orientation==2) {
						if (orientation==0) tmpEdge.nAlignmentsForward++;
						else tmpEdge.nAlignmentsBackward++;
					}
					if (tmpEdge.avgDiffs!=-1) tmpEdge.avgDiffs+=normalizedDiffs;
					return true;
				}
				else return false;
			}
		}
		
		// Creating a new edge
		tmpAlignmentPair.node1=nodesArray[n1];
		tmpAlignmentPair.alignmentStart1=nodesArray[n1].start;
		tmpAlignmentPair.alignmentEnd1=nodesArray[n1].end;
		tmpAlignmentPair.node2=nodesArray[n2];
		tmpAlignmentPair.alignmentStart2=nodesArray[n2].start;
		tmpAlignmentPair.alignmentEnd2=nodesArray[n2].end;
		tmpAlignmentPair.orientation=orientation==0;
		tmpAlignmentPair.diffs=normalizedDiffs;
		tmpAlignmentPair.inCanonicalForm();
		tmpEdge = new Edge(tmpAlignmentPair);
		tmpEdge.supplement=true;
		if (neighbors[n1].length==nNeighbors[n1]) {
			tmpNeighbors = new Edge[nNeighbors[n1]+RESIZE_UNIT];
			System.arraycopy(neighbors[n1],0,tmpNeighbors,0,nNeighbors[n1]);
			neighbors[n1]=tmpNeighbors;
		}
		neighbors[n1][nNeighbors[n1]++]=tmpEdge;
		if (neighbors[n2].length==nNeighbors[n2]) {
			tmpNeighbors = new Edge[nNeighbors[n2]+RESIZE_UNIT];
			System.arraycopy(neighbors[n2],0,tmpNeighbors,0,nNeighbors[n2]);
			neighbors[n2]=tmpNeighbors;
		}
		neighbors[n2][nNeighbors[n2]++]=tmpEdge;
		return true;
	}
	
	
	/**
	 * Like $addIdentityEdges()$, but adds just proper containment edges in which the
	 * contained node is a kernel, and the container node might not be a kernel.
	 * Containment edges are added only if they satisfy $IntervalGraphStep2.inStep2()$, 
	 * i.e. if they are true containments but they do not belong to the interval graph 
	 * since the alignments that prove the containment were not assigned to the intervals.
	 * Substring intervals are not used as container nodes for non-substring intervals, 
	 * since their containment rules are not specific.
	 */
	public static final int addContainmentEdges(boolean kernelsComputed, String alignmentsFilePath) throws IOException {
		final int THRESHOLD = IO.quantum*3;  // Arbitrary
		final int PREVIOUS_ORDER = Node.order;
		int i, j, p, q;
		int readA, readB, nAdded, nAlignments, nEdges, diffs;
		String str;
		Node tmpNode = new Node(-1,-1,-1,-1,-1,false,false,false,0,-1,-1,-1,false,-1,-1,0,0);
		Edge tmpEdge = new Edge();
		AlignmentPair tmpAlignmentPair = new AlignmentPair();
		BufferedReader alignmentsFile;
		
		System.err.println("Adding supplement edges of containment type...");		
		Node.order=Node.READ_START_NODEID;
		nEdges=getNEdges();
		alignmentsFile = new BufferedReader(new FileReader(alignmentsFilePath),IO.BUFFER_SIZE);
		str=alignmentsFile.readLine();
		str=alignmentsFile.readLine();  // Skipping the first two lines
		str=alignmentsFile.readLine();
		nAdded=0; nAlignments=0;
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			nAlignments++;
			readA=Alignments.readA-1;
			tmpNode.read=readA; tmpNode.start=Alignments.startA; tmpNode.id=-1;
			p=Arrays.binarySearch(nodesArray,0,nNodes,tmpNode);
			if (p<0) {
				p=-1-p;
				p=Math.min(p,nNodes-1);
			}
			for (i=p-1; i>=0; i--) {
				if (nodesArray[i].read!=readA || nodesArray[i].start<Alignments.startA-THRESHOLD) break;
				if (nodesArray[i].end<=Alignments.startA || nodesArray[i].end>Alignments.endA+THRESHOLD) continue;
				if (kernelsComputed?nodesArray[i].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[i])) continue;
				readB=Alignments.readB-1;
				tmpNode.read=readB; tmpNode.start=Alignments.startB; tmpNode.id=-1;
				q=Arrays.binarySearch(nodesArray,0,nNodes,tmpNode);
				if (q<0) {
					q=-1-q;
					q=Math.min(q,nNodes-1);
				}
				for (j=q-1; j>=0; j--) {
					if (nodesArray[j].read!=readB) break;
					if (j==i) continue;
					if ( (nodesArray[j].type==Constants.INTERVAL_DENSE_SUBSTRING && nodesArray[i].type!=Constants.INTERVAL_DENSE_SUBSTRING) || 
						 (nodesArray[j].type==Constants.INTERVAL_PERIODIC && nodesArray[i].type!=Constants.INTERVAL_PERIODIC)
					   ) continue;
					if (nodesArray[j].end<=Alignments.startB) continue;
					if (Intervals.isApproximatelyContained(nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,tmpArray)) {
						diffs=Alignments.projectDiffs(Alignments.diffs,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,nodesArray[i].start,nodesArray[i].end,tmpArray[0],tmpArray[1]);
						if (addContainmentEdge(i,j,Alignments.orientation?0:1,diffs,tmpArray[0],tmpArray[1],tmpAlignmentPair,tmpEdge)) nAdded++;
					}
				}
				for (j=q; j<nNodes; j++) {
					if (nodesArray[j].read!=readB || nodesArray[j].start>=Alignments.endB) break;
					if (j==i) continue;
					if ( (nodesArray[j].type==Constants.INTERVAL_DENSE_SUBSTRING && nodesArray[i].type!=Constants.INTERVAL_DENSE_SUBSTRING) || 
						 (nodesArray[j].type==Constants.INTERVAL_PERIODIC && nodesArray[i].type!=Constants.INTERVAL_PERIODIC)
					   ) continue;
					if (Intervals.isApproximatelyContained(nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,tmpArray)) {
						diffs=Alignments.projectDiffs(Alignments.diffs,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,nodesArray[i].start,nodesArray[i].end,tmpArray[0],tmpArray[1]);
						if (addContainmentEdge(i,j,Alignments.orientation?0:1,diffs,tmpArray[0],tmpArray[1],tmpAlignmentPair,tmpEdge)) nAdded++;
					}
				}
			}
			for (i=p; i<nNodes; i++) {
				if (nodesArray[i].read!=readA || nodesArray[i].start>=Alignments.endA) break;
				if (nodesArray[i].end>Alignments.endA+THRESHOLD) continue;
				if (kernelsComputed?nodesArray[i].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[i])) continue;
				readB=Alignments.readB-1;
				tmpNode.read=readB; tmpNode.start=Alignments.startB; tmpNode.id=-1;
				q=Arrays.binarySearch(nodesArray,0,nNodes,tmpNode);
				if (q<0) {
					q=-1-q;
					q=Math.min(q,nNodes-1);
				}
				for (j=q-1; j>=0; j--) {
					if (nodesArray[j].read!=readB) break;
					if (j==i) continue;
					if ( (nodesArray[j].type==Constants.INTERVAL_DENSE_SUBSTRING && nodesArray[i].type!=Constants.INTERVAL_DENSE_SUBSTRING) || 
						 (nodesArray[j].type==Constants.INTERVAL_PERIODIC && nodesArray[i].type!=Constants.INTERVAL_PERIODIC)
					   ) continue;
					if (nodesArray[j].end<=Alignments.startB) continue;
					if (Intervals.isApproximatelyContained(nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,tmpArray)) {
						diffs=Alignments.projectDiffs(Alignments.diffs,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,nodesArray[i].start,nodesArray[i].end,tmpArray[0],tmpArray[1]);
						if (addContainmentEdge(i,j,Alignments.orientation?0:1,diffs,tmpArray[0],tmpArray[1],tmpAlignmentPair,tmpEdge)) nAdded++;
					}
				}
				for (j=q; j<nNodes; j++) {
					if (nodesArray[j].read!=readB || nodesArray[j].start>=Alignments.endB) break;
					if (j==i) continue;
					if ( (nodesArray[j].type==Constants.INTERVAL_DENSE_SUBSTRING && nodesArray[i].type!=Constants.INTERVAL_DENSE_SUBSTRING) || 
						 (nodesArray[j].type==Constants.INTERVAL_PERIODIC && nodesArray[i].type!=Constants.INTERVAL_PERIODIC)
					   ) continue;
					if (Intervals.isApproximatelyContained(nodesArray[i].start,nodesArray[i].end,nodesArray[j].start,nodesArray[j].end,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,tmpArray)) {
						diffs=Alignments.projectDiffs(Alignments.diffs,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,nodesArray[i].start,nodesArray[i].end,tmpArray[0],tmpArray[1]);
						if (addContainmentEdge(i,j,Alignments.orientation?0:1,diffs,tmpArray[0],tmpArray[1],tmpAlignmentPair,tmpEdge)) nAdded++;
					}
				}
			}
			if (nAlignments%1000000==0) System.err.println("Processed alignment "+(nAlignments/1000000)+"M");
			str=alignmentsFile.readLine();
		}
		alignmentsFile.close();
		Node.order=PREVIOUS_ORDER;
		System.err.println("Added "+nAdded+" supplement edges of containment type ("+IO.getPercent(nAdded,nEdges)+"%)");
		return nAdded;
	}
	
	
	/**
	 * Adds an edge in which $n1$ is contained in $n2$, iff such edge conforms to 
	 * $IntervalGraphStep2.inStep2()$.
	 *
	 * Remark: the procedure does not assume $neighbors[i]$ to be sorted, and it considers
	 * all edges in $neighbors[i]$.
	 *
	 * @param orientation 0=forward; 1=RC;
	 * @param diffs estimate of the number of diffs to be assigned to the edge (derived 
	 * from the diffs of some alignment); this is a number, not a rate;
	 * @param start2,end2 absolute positions, in the read of node $n2$, of the substring 
	 * of $n2$ to which $n1$ is projected;
	 * @param tmpAlignmentPair temporary space;
	 * @param newEdge temporary space;
	 * @return TRUE iff a new edge is added by the procedure.
	 */
	private static final boolean addContainmentEdge(int n1, int n2, int orientation, int diffs, int start2, int end2, AlignmentPair tmpAlignmentPair, Edge newEdge) {
		final int RESIZE_UNIT = 10;  // Arbitrary
		int i, n, m;
		int type;
		final double normalizedDiffs = 2.0*diffs/(nodesArray[n1].end-nodesArray[n1].start+end2-start2+2);
		Edge tmpEdge;
		Edge[] tmpNeighbors;
		
		// Creating a new edge
		tmpAlignmentPair.node1=nodesArray[n1];
		tmpAlignmentPair.alignmentStart1=nodesArray[n1].start;
		tmpAlignmentPair.alignmentEnd1=nodesArray[n1].end;
		tmpAlignmentPair.node2=nodesArray[n2];
		tmpAlignmentPair.alignmentStart2=start2;
		tmpAlignmentPair.alignmentEnd2=end2;
		tmpAlignmentPair.orientation=orientation==0;
		tmpAlignmentPair.diffs=normalizedDiffs;
		tmpAlignmentPair.inCanonicalForm();
		type=IntervalGraphStep2.inStep2(tmpAlignmentPair);
		if (!(type==Constants.CONTAINMENT_ONE_IN_TWO && tmpAlignmentPair.node1==nodesArray[n1]) && !(type==Constants.CONTAINMENT_TWO_IN_ONE && tmpAlignmentPair.node2==nodesArray[n1])) return false;
		
		// Checking whether the edge already exists
		newEdge.set(tmpAlignmentPair);
		n=nNeighbors[n1]<nNeighbors[n2]?n1:n2;
		m=n==n1?n2:n1;
		for (i=0; i<nNeighbors[n]; i++) {
			tmpEdge=neighbors[n][i];
			if (tmpEdge.supplement) continue;
			if (tmpEdge.getTo(n)==m) {
				if (tmpEdge.containment!=-1 || tmpEdge.insertion!=-1) return false;
				break;
			}
		}
		for (i=0; i<nNeighbors[n]; i++) {
			tmpEdge=neighbors[n][i];
			if (!tmpEdge.supplement) continue;
			if (tmpEdge.getTo(n)==m) {
				if (tmpEdge.containment==-1 && tmpEdge.insertion==-1) {
					tmpEdge.addTypes(newEdge);
					tmpEdge.addDiffs(newEdge);
					tmpEdge.addOrientation(newEdge);
					return true;
				}
				else return false;
			}
		}
		
		// Adding the edge
		tmpEdge = new Edge(tmpAlignmentPair);
		tmpEdge.supplement=true;
		if (neighbors[n1].length==nNeighbors[n1]) {
			tmpNeighbors = new Edge[nNeighbors[n1]+RESIZE_UNIT];
			System.arraycopy(neighbors[n1],0,tmpNeighbors,0,nNeighbors[n1]);
			neighbors[n1]=tmpNeighbors;
		}
		neighbors[n1][nNeighbors[n1]++]=tmpEdge;
		if (neighbors[n2].length==nNeighbors[n2]) {
			tmpNeighbors = new Edge[nNeighbors[n2]+RESIZE_UNIT];
			System.arraycopy(neighbors[n2],0,tmpNeighbors,0,nNeighbors[n2]);
			neighbors[n2]=tmpNeighbors;
		}
		neighbors[n2][nNeighbors[n2]++]=tmpEdge;
		return true;
	}
	
	
	/**
	 * Tries to use all alignments in $alignmentsFilePath$ to add supplement edges between 
	 * pairs of neighbors $u_1,u_2$ of the same overlap fork $v$, that overlap with $v$ on 
	 * the same side of $v$. This is useful, since an alignment that says something about 
	 * the similarity of $u_1$ and $u_2$ might not have been assigned to $u_1$ or to $u_2$
	 * during factorization, and thus such information might be missing from the interval
	 * graph.
	 *
	 * Remark: this procedure is similar to $addIdentityEdges()$. See the latter for an
	 * explanation of why it is not necessary to build the whole interval graph in the
	 * same way. In practice, this procedure makes the graphs significantly larger, and
	 * marking the edges created by this procedure as regular rather than supplement would
	 * significantly slow down the procedures downstream.
	 *
	 * Remark: for speed, the procedure considers just pairs of intervals that can belong 
	 * to kernels in Step 3.
	 *
	 * Remark: the procedure does not assume $nodesArray$ to be sorted, and it does not 
	 * create new intervals. The procedure does not assume $neighbors[i]$ to be sorted for
	 * any $i$, and it considers all edges in $neighbors[i]$.
	 *
	 * Remark: the procedure sets field $visited=1$ for all forks that it considered (but
	 * for which it did not necessarily add edges). Forks that were already marked are
	 * skipped.
	 *
	 * @param kernelsComputed if TRUE, assumes that kernel tags have already been assigned
	 * to nodes;
	 * @param alignmentsFilePath assumed to contain all and only the alignments between 
	 * reads such that at least one node of the graph occurs in each read. The same 
	 * alignment is expected to occur twice, once for each read considered as readA.
	 * @return TRUE iff the procedure created at least one new overlap edge (possibly 
	 * weak, i.e. not classified as overlap).
	 */
	public static final boolean addEdgesToOverlapNeighbors(boolean kernelsComputed, String alignmentsFilePath) throws IOException {
		final int TUPLES_UNIT = 50;  // Arbitrary
		final int MIN_ALIGNMENT_LENGTH = Alignments.minAlignmentLength>>2;  // Arbitrary
		boolean success, addEdge, out, outPrime;
		int i, j, n;
		int nPairs, firstSameReadPair;
		int alignmentRead1, alignmentRead2, alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2, alignmentReadA;
		int pairRead, pairRead1, pairRead2, blockStart, blockRead, blockRead1, blockRead2;
		int overhangStart1, overhangStart2, overhangEnd1, overhangEnd2;
		int diffs, nAdded, nEdges;
		String str;
		BufferedReader alignmentsFile;
		AlignmentPair tmpAlignmentPair = new AlignmentPair(null,null,-1,-1,-1,-1,false,-1,false);
		Edge edge, previousEdge;
		int[] tmpArray1 = new int[2];
		int[] tmpArray2 = new int[2];
		OverlapPair[] pairs;
		AlignmentTuple[] tuples, tmpTuples;
		
		// Building $pairs$.
		pairs=addEdgesToOverlapNeighbors_getPairs(kernelsComputed,tmpArray);
		if (pairs==null) return false;
		nPairs=tmpArray[0]+1;
		firstSameReadPair=nPairs;
		for (i=nPairs-1; i>=0; i--) {
			if (pairs[i].read1!=pairs[i].read2) {
				firstSameReadPair=i+1;
				break;
			}
		}
		nAdded=0; out=false; nEdges=getNEdges();
		
if (IO.SHOW_INTERACTIVE) {
System.err.println("addEdgesToOverlapNeighbors> PAIRS at the very beginning:  nPairs="+nPairs);
for (int x=0; x<nPairs; x++) System.err.println(pairs[x]);		
}
		
		// Fixing overlap pairs in different reads.
		//
		// Remark: the procedure skips just through blocks of pairs with the same
		// $read1,read2$, but it does not use the first/last positions of overhangs to 
	    // skip pairs. This is because alignments are read one by one from disk, so we
	 	// don't know the next alignment unless we load it from disk, thus we cannot
		// maintain a $firstIForNextAlignment$ variable as we scan $pairs$. Of course we
		// could look ahead in the disk buffer: we don't implement this for simplicity.
		if (firstSameReadPair>0) {
			System.err.println("Adding supplement edges between pairs of overlap neighbors (different reads)...");
			i=0; pairRead1=pairs[i].read1; pairRead2=pairs[i].read2;
			blockStart=0; blockRead1=pairRead1; blockRead2=pairRead2;
			alignmentsFile = new BufferedReader(new FileReader(alignmentsFilePath),IO.BUFFER_SIZE);
			str=alignmentsFile.readLine();
			str=alignmentsFile.readLine();  // Skipping the first two lines
			str=alignmentsFile.readLine();
			success=Alignments.readAlignmentFile(str);
			if (success) {
				if (Alignments.readA<Alignments.readB) {
					alignmentRead1=Alignments.readA-1; 
					alignmentRead2=Alignments.readB-1;
					alignmentStart1=Alignments.startA;
					alignmentEnd1=Alignments.endA;
					alignmentStart2=Alignments.startB;
					alignmentEnd2=Alignments.endB;
				}
				else {
					alignmentRead1=Alignments.readB-1; 
					alignmentRead2=Alignments.readA-1;
					alignmentStart1=Alignments.startB;
					alignmentEnd1=Alignments.endB;
					alignmentStart2=Alignments.startA;
					alignmentEnd2=Alignments.endA;
				}
			}
			else {
				alignmentRead1=-1;
				alignmentRead2=-1;
				alignmentStart1=-1;
				alignmentEnd1=-1;
				alignmentStart2=-1;
				alignmentEnd2=-1;
			}
			n=0;
			while (success) {
				if (i==firstSameReadPair || pairRead1>alignmentRead1 || (pairRead1==alignmentRead1 && pairRead2>alignmentRead2)) {
					str=alignmentsFile.readLine();
					success=Alignments.readAlignmentFile(str);
					if (!success) break;					
					if (Alignments.readA<Alignments.readB) {
						alignmentRead1=Alignments.readA-1; 
						alignmentRead2=Alignments.readB-1;
						alignmentStart1=Alignments.startA;
						alignmentEnd1=Alignments.endA;
						alignmentStart2=Alignments.startB;
						alignmentEnd2=Alignments.endB;
					}
					else {
						alignmentRead1=Alignments.readB-1; 
						alignmentRead2=Alignments.readA-1;
						alignmentStart1=Alignments.startB;
						alignmentEnd1=Alignments.endB;
						alignmentStart2=Alignments.startA;
						alignmentEnd2=Alignments.endA;
					}
					if (alignmentRead1==blockRead1 && alignmentRead2==blockRead2) {
						i=blockStart;
						pairRead1=pairs[i].read1; pairRead2=pairs[i].read2;
					}
					else {
						if (i==firstSameReadPair) break;
						blockStart=i;
						blockRead1=pairRead1; blockRead2=pairRead2;
					}
					n++;
					if (n%1000000==0) System.err.println("Processed alignment "+(n/1000000)+"M");
					continue;
				}
				if (pairRead1<alignmentRead1 || (pairRead1==alignmentRead1 && pairRead2<alignmentRead2)) {
					i++;
					if (i<firstSameReadPair) {
						pairRead1=pairs[i].read1; pairRead2=pairs[i].read2;
						blockStart=i;
						blockRead1=pairRead1; blockRead2=pairRead2;
					}
					continue;
				}
				if ( ( Intervals.areApproximatelyIdentical(alignmentStart1,alignmentEnd1,nodesArray[pairs[i].neighbor1].start,nodesArray[pairs[i].neighbor1].end) ||
					   Intervals.isApproximatelyContained(nodesArray[pairs[i].neighbor1].start,nodesArray[pairs[i].neighbor1].end,alignmentStart1,alignmentEnd1)
					 ) &&
					 ( Intervals.areApproximatelyIdentical(alignmentStart2,alignmentEnd2,nodesArray[pairs[i].neighbor2].start,nodesArray[pairs[i].neighbor2].end) ||
					   Intervals.isApproximatelyContained(nodesArray[pairs[i].neighbor2].start,nodesArray[pairs[i].neighbor2].end,alignmentStart2,alignmentEnd2)
					 )
				   ) {
	   				// Skipping pairs for which the alignment involves the entire
					// intervals of the neighbors, since such alignments should have
					// already been considered by $addIdentityEdges()$.
					i++;
					if (i<firstSameReadPair) {
						pairRead1=pairs[i].read1;
						pairRead2=pairs[i].read2;
					}
					continue;
				}
				overhangStart1=pairs[i].overhangStart(true);
				overhangEnd1=pairs[i].overhangEnd(true);
				overhangStart2=pairs[i].overhangStart(false);
				overhangEnd2=pairs[i].overhangEnd(false);
				addEdge=false;				
				if (Alignments.orientation) {
					if ( ( (pairs[i].overhang1>0 && pairs[i].overhang2>0) && 
						   ( Intervals.areApproximatelyIdentical(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,Alignments.orientation) ||
						     Intervals.isApproximatePrefixOrSuffix(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,Alignments.orientation,false) ||
						     Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,overhangStart1,overhangEnd1,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,Alignments.orientation,false)
						   ) 
					     ) ||
					     ( (pairs[i].overhang1<0 && pairs[i].overhang2<0) &&
						   ( Intervals.areApproximatelyIdentical(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,Alignments.orientation) ||
						     Intervals.isApproximatePrefixOrSuffix(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,Alignments.orientation,true) ||
						     Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,overhangStart1,overhangEnd1,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,Alignments.orientation,true)
						   )
						 ) 
					   ) addEdge=true;
				}
				else {
					if ( ( (pairs[i].overhang1>0 && pairs[i].overhang2<0) &&
						   ( Intervals.areApproximatelyIdentical(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,Alignments.orientation) ||
						     Intervals.isApproximatePrefixOrSuffix(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,Alignments.orientation,true) ||
						     Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,overhangStart1,overhangEnd1,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,Alignments.orientation,false)
						   )
					     ) ||
					     ( (pairs[i].overhang1<0 && pairs[i].overhang2>0) &&
						   ( Intervals.areApproximatelyIdentical(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,Alignments.orientation) ||
						     Intervals.isApproximatePrefixOrSuffix(overhangStart1,overhangEnd1,overhangStart2,overhangEnd2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,Alignments.orientation,false) ||
						     Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,overhangStart1,overhangEnd1,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,Alignments.orientation,true)
						   ) 
					     ) 
					   ) addEdge=true;
				}
				if (addEdge) {
					tmpAlignmentPair.node1=nodesArray[pairs[i].neighbor1];
					tmpAlignmentPair.node2=nodesArray[pairs[i].neighbor2];
					tmpAlignmentPair.alignmentStart1=Math.max(alignmentStart1,tmpAlignmentPair.node1.start);
					tmpAlignmentPair.alignmentEnd1=Math.min(alignmentEnd1,tmpAlignmentPair.node1.end);
					if (tmpAlignmentPair.alignmentEnd1-tmpAlignmentPair.alignmentStart1+1<MIN_ALIGNMENT_LENGTH) {
						i++;
						if (i<firstSameReadPair) {
							pairRead1=pairs[i].read1;
							pairRead2=pairs[i].read2;
						}
						continue;
					}
					tmpAlignmentPair.alignmentStart2=Math.max(alignmentStart2,tmpAlignmentPair.node2.start);
					tmpAlignmentPair.alignmentEnd2=Math.min(alignmentEnd2,tmpAlignmentPair.node2.end);
					if (tmpAlignmentPair.alignmentEnd2-tmpAlignmentPair.alignmentStart2+1<MIN_ALIGNMENT_LENGTH) {
						i++;
						if (i<firstSameReadPair) {
							pairRead1=pairs[i].read1;
							pairRead2=pairs[i].read2;
						}
						continue;
					}
					tmpAlignmentPair.orientation=Alignments.orientation;
					diffs=Alignments.projectDiffs(Alignments.diffs,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,tmpAlignmentPair.alignmentStart1,tmpAlignmentPair.alignmentEnd1,tmpAlignmentPair.alignmentStart2,tmpAlignmentPair.alignmentEnd2);
					tmpAlignmentPair.diffs=2.0*diffs/(tmpAlignmentPair.alignmentEnd1-tmpAlignmentPair.alignmentStart1+tmpAlignmentPair.alignmentEnd2-tmpAlignmentPair.alignmentStart2+2);					
					tmpAlignmentPair.inCanonicalForm();
					edge = new Edge(tmpAlignmentPair);
					edge.supplement=true;
					j=findEdge(tmpAlignmentPair.node1.nodeID,edge,false);
					if (j==Math.NEGATIVE_INFINITY) {
						j=findEdge(tmpAlignmentPair.node1.nodeID,edge,true);
						if (j==Math.NEGATIVE_INFINITY) {
							addEdge(tmpAlignmentPair.node1.nodeID,edge);
							addEdge(tmpAlignmentPair.node2.nodeID,edge);
							if (edge.overlap!=-1 || edge.isOverlap_weak(IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2)) out=true;
							nAdded++;
						}
						else {
							previousEdge=neighbors[tmpAlignmentPair.node1.nodeID][j];
							if (!previousEdge.contains(edge,IDENTITY_THRESHOLD)) {
								previousEdge.addTypes(edge);
								previousEdge.addDiffs(edge);
								previousEdge.addOrientation(edge);
								if (edge.overlap!=-1 || edge.isOverlap_weak(IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2)) out=true;
								nAdded++;
							}
							else {
								// The supplement edge is already contained in a
								// supplement edge.
							}
						}
					}
					else {
						previousEdge=neighbors[tmpAlignmentPair.node1.nodeID][j];
						if (!previousEdge.contains(edge,IDENTITY_THRESHOLD)) {
							j=findEdge(tmpAlignmentPair.node1.nodeID,edge,true);
							if (j==Math.NEGATIVE_INFINITY) {
								addEdge(tmpAlignmentPair.node1.nodeID,edge);
								addEdge(tmpAlignmentPair.node2.nodeID,edge);
								if (edge.overlap!=-1 || edge.isOverlap_weak(IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2)) out=true;
								nAdded++;
							}
							else {
								previousEdge=neighbors[tmpAlignmentPair.node1.nodeID][j];
								if (!previousEdge.contains(edge,IDENTITY_THRESHOLD)) {
									previousEdge.addTypes(edge);
									previousEdge.addDiffs(edge);
									previousEdge.addOrientation(edge);
									if (edge.overlap!=-1 || edge.isOverlap_weak(IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2)) out=true;
									nAdded++;
								}
								else {
									// The supplement edge is already contained in a
									// supplement edge.
								}
							}
						}
						else {
							// The supplement edge is already contained in a regular edge
						}
					}
				}
				i++;
				if (i<firstSameReadPair) {
					pairRead1=pairs[i].read1;
					pairRead2=pairs[i].read2;
				}
			}
			alignmentsFile.close(); alignmentsFile=null;
			System.err.println("Added "+nAdded+" supplement edges between pairs of overlap neighbors (different reads, "+IO.getPercent(nAdded,nEdges)+"%)");
		}
		
		// Fixing overlap pairs in the same read
		if (firstSameReadPair<nPairs) {
			System.err.println("Adding supplement edges between pairs of overlap neighbors (same read)...");
			tuples = new AlignmentTuple[TUPLES_UNIT];
			nAdded=0; j=-1;
			i=firstSameReadPair; pairRead=pairs[i].read1;
			blockStart=i; blockRead=pairRead;
			alignmentsFile = new BufferedReader(new FileReader(alignmentsFilePath),IO.BUFFER_SIZE);
			alignmentsFile.readLine();
			alignmentsFile.readLine();  // Skipping the first two lines
			str=alignmentsFile.readLine();
			success=Alignments.readAlignmentFile(str);
			while (success) {
				alignmentReadA=Alignments.readA-1;
				if (alignmentReadA<blockRead) {
					str=alignmentsFile.readLine();
					success=Alignments.readAlignmentFile(str);
					continue;
				}
				if (alignmentReadA>blockRead) {
					if (j!=-1) {
						if (j>0) Arrays.sort(tuples,0,j+1);
						outPrime=addEdgesToOverlapNeighbors_sameRead(pairs,nPairs-1,blockStart,tuples,j,tmpArray1,tmpArray2);
						if (outPrime) out=true;
						nAdded+=tmpArray[0];
						i=tmpArray[1];
					}
					else i=nextRead1(pairs,blockStart);
					if (i%1000==0) System.err.println("Processed pair "+(i%1000)+"K");
					if (i==nPairs) break;
					pairRead=pairs[i].read1;
					blockStart=i; blockRead=pairRead;
					j=-1;
					continue;
				}
				j++;
				if (j==tuples.length) {
					tmpTuples = new AlignmentTuple[tuples.length];
					System.arraycopy(tuples,0,tmpTuples,0,tuples.length);
					tuples = new AlignmentTuple[tuples.length+TUPLES_UNIT];
					System.arraycopy(tmpTuples,0,tuples,0,tmpTuples.length);
					tmpTuples=null;
				}
				if (tuples[j]==null) tuples[j] = new AlignmentTuple(Alignments.orientation,Alignments.readA,Alignments.readB,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.diffs);
				else tuples[j].set(Alignments.orientation,Alignments.readA,Alignments.readB,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.diffs);
				str=alignmentsFile.readLine();
				success=Alignments.readAlignmentFile(str);
			}
			if (j!=-1) {
				if (j>0) Arrays.sort(tuples,0,j+1);
				outPrime=addEdgesToOverlapNeighbors_sameRead(pairs,nPairs-1,blockStart,tuples,j,tmpArray1,tmpArray2);
				if (outPrime) out=true;
				nAdded+=tmpArray[0];
			}
			alignmentsFile.close(); alignmentsFile=null;
			System.err.println("Added "+nAdded+" supplement edges between pairs of overlap neighbors (same read, "+IO.getPercent(nAdded,nEdges)+"%)");
		}
		
		return out;
	}
	
	
	/**
	 * @return the first position after $blockStart$ in $pairs$ with a different value of
	 * $read1$, or $pairs.length$ if none is found.
	 */
	private static final int nextRead1(OverlapPair[] pairs, int blockStart) {
		final int blockRead = pairs[blockStart].read1;
		int i = blockStart+1;
		while (i<pairs.length && pairs[i].read1==blockRead) i++;
		return i;
	}
	
	
	/**
	 * Adds edges to pairs of neighbors in the same readA, if there are two alignments 
	 * that project their overhangs to the same substring of the same readB, or that 
	 * project one overhang to a prefix or suffix of the other in the same readB.
	 * The procedure stores in $tmpArray[0]$ the number of edges added by the procedure, 
	 * and in $tmpArray[1]$ the first position in $pairs$ after the current block.
	 *
	 * Remark: this is useful, since the pipeline assumes that same-read alignments have
	 * not been computed. However, the creation of such artificial alignments should be
	 * done sparingly (e.g. just to simplify overlap forks), since transitivity is not 
	 * true in general.
	 *
	 * Remark: the procedure does not assume $neighbors[i]$ to be sorted for any $i$, and
	 * it considers all edges in $neighbors[i]$.
	 *
	 * Remark: the procedure uses $stack$ as temporary space.
	 * 
	 * @param pairs sorted by $read1,start1$;
	 * @param blockStart all pairs with the same $read1=read2=X$ for a given $X$ are 
	 * assumed to form a compact block that starts at index $blockStart$ in $pairs$;
	 * @param tuples[0..lastTuple] all alignments with $readA=X$, sorted by 
	 * $readA,readB,startA$;
	 * @return TRUE iff the procedure created at least one new overlap edge (possibly 
	 * weak, i.e. not classified as overlap).
	 */
	private static final boolean addEdgesToOverlapNeighbors_sameRead(OverlapPair[] pairs, int lastPair, int blockStart, AlignmentTuple[] tuples, int lastTuple, int[] tmpArray1, int[] tmpArray2) {
		boolean success, addEdge, orientationI, out;
		int i, j, k, p;
		final int blockRead;
		int overhangStart1, overhangEnd1, overhangStart2, overhangEnd2;
		int projectionStartI, projectionEndI, readBI, firstIForNextP, nAdded;
		AlignmentPair tmpAlignmentPair = new AlignmentPair(null,null,-1,-1,-1,-1,false,-1,false);
		Edge edge, previousEdge;
		
		nAdded=0; out=false;
		p=blockStart; blockRead=pairs[p].read1; i=0; firstIForNextP=-1;
		while (p<=lastPair && pairs[p].read1==blockRead) {
			overhangStart1=pairs[p].overhangStart(true);
			overhangEnd1=pairs[p].overhangEnd(true);
			overhangStart2=pairs[p].overhangStart(false);
			overhangEnd2=pairs[p].overhangEnd(false);
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> CONSIDERING PAIR "+pairs[p]);	
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> overhang1: ["+overhangStart1+".."+overhangEnd1+"] overhang2: ["+overhangStart2+".."+overhangEnd2+"]");
			while (i<=lastTuple) {
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> 0  considering alignment tuple: "+tuples[i]);
				if (tuples[i].startA>=nodesArray[pairs[p].neighbor1].end) break;
				if (tuples[i].endA<nodesArray[pairs[p].neighbor1].start) {
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> FAIL 1");
					i++;
					continue;
				}
				if (firstIForNextP==-1 && p+1<=lastPair && pairs[p+1].read1==blockRead && tuples[i].endA>=nodesArray[pairs[p+1].neighbor1].start) firstIForNextP=i;
				if ( !Intervals.areApproximatelyIdentical(tuples[i].startA,tuples[i].endA,overhangStart1,overhangEnd1) &&
					 !Intervals.isApproximatelyContained(overhangStart1,overhangEnd1,tuples[i].startA,tuples[i].endA)
				   ) {
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> FAIL 2");
					i++;
					continue;
				}
				readBI=tuples[i].readB;
				orientationI=tuples[i].orientation;
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> 1");
				success=Intervals.project(overhangStart1,overhangEnd1,tuples[i].startA,tuples[i].endA,tuples[i].startB,tuples[i].endB,orientationI,stack,0);
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> 2  success="+success);
				if (!success) {
					i++;
					continue;
				}
				projectionStartI=stack[0]; projectionEndI=stack[1];
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> 3  projectionI: ["+projectionStartI+".."+projectionEndI+"]");
				for (j=i+1; j<=lastTuple; j++) {
					if (tuples[j].readB!=readBI || tuples[j].startA>=nodesArray[pairs[p].neighbor2].end) break;
					if (tuples[j].endA<nodesArray[pairs[p].neighbor2].start) continue;
					if ( !Intervals.areApproximatelyIdentical(tuples[j].startA,tuples[j].endA,overhangStart2,overhangEnd2) &&
						 !Intervals.isApproximatelyContained(overhangStart2,overhangEnd2,tuples[j].startA,tuples[j].endA)
					   ) continue;
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> 4  considering other tuple: "+tuples[j]);
					addEdge=false;
					if (tuples[j].orientation==orientationI) {
						if ( ( (pairs[p].overhang1>0 && pairs[p].overhang2>0) && 
							   ( Intervals.areApproximatelyIdentical(overhangStart2,overhangEnd2,projectionStartI,projectionEndI,tuples[j].startA,tuples[j].endA,tuples[j].startB,tuples[j].endB,tuples[j].orientation) ||
							     Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,projectionStartI,projectionEndI,tuples[j].startA,tuples[j].endA,tuples[j].startB,tuples[j].endB,tuples[j].orientation,false) ||
							     Intervals.isApproximatePrefixOrSuffix(projectionStartI,projectionEndI,overhangStart2,overhangEnd2,tuples[j].startB,tuples[j].endB,tuples[j].startA,tuples[j].endA,tuples[j].orientation,false)
							   ) 
						     ) ||
						     ( (pairs[p].overhang1<0 && pairs[p].overhang2<0) &&
							   ( Intervals.areApproximatelyIdentical(overhangStart2,overhangEnd2,projectionStartI,projectionEndI,tuples[j].startA,tuples[j].endA,tuples[j].startB,tuples[j].endB,tuples[j].orientation) ||
							     Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,projectionStartI,projectionEndI,tuples[j].startA,tuples[j].endA,tuples[j].startB,tuples[j].endB,tuples[j].orientation,true) ||
							     Intervals.isApproximatePrefixOrSuffix(projectionStartI,projectionEndI,overhangStart2,overhangEnd2,tuples[j].startB,tuples[j].endB,tuples[j].startA,tuples[j].endA,tuples[j].orientation,true)
							   )
							 ) 
						   ) addEdge=true;
					}
					else {
						if ( ( (pairs[p].overhang1>0 && pairs[p].overhang2<0) &&
							   ( Intervals.areApproximatelyIdentical(overhangStart2,overhangEnd2,projectionStartI,projectionEndI,tuples[j].startA,tuples[j].endA,tuples[j].startB,tuples[j].endB,tuples[j].orientation) ||
							     Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,projectionStartI,projectionEndI,tuples[j].startA,tuples[j].endA,tuples[j].startB,tuples[j].endB,tuples[j].orientation,false) ||
								 Intervals.isApproximatePrefixOrSuffix(projectionStartI,projectionEndI,overhangStart2,overhangEnd2,tuples[j].startB,tuples[j].endB,tuples[j].startA,tuples[j].endA,tuples[j].orientation,true)
							   )
						     ) ||
						     ( (pairs[p].overhang1<0 && pairs[p].overhang2>0) &&
							   ( Intervals.areApproximatelyIdentical(projectionStartI,projectionEndI,overhangStart2,overhangEnd2,tuples[j].startB,tuples[j].endB,tuples[j].startA,tuples[j].endA,tuples[j].orientation) ||
							     Intervals.isApproximatePrefixOrSuffix(overhangStart2,overhangEnd2,projectionStartI,projectionEndI,tuples[j].startA,tuples[j].endA,tuples[j].startB,tuples[j].endB,tuples[j].orientation,true) ||
								 Intervals.isApproximatePrefixOrSuffix(projectionStartI,projectionEndI,overhangStart2,overhangEnd2,tuples[j].startB,tuples[j].endB,tuples[j].startA,tuples[j].endA,tuples[j].orientation,false)
							   ) 
						     ) 
						   ) addEdge=true;
					}
					if (!addEdge) continue;
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> 5");
					success=Intervals.project(nodesArray[pairs[p].neighbor1].start,nodesArray[pairs[p].neighbor1].end,nodesArray[pairs[p].neighbor2].start,nodesArray[pairs[p].neighbor2].end,tuples[i].startA,tuples[i].endA,tuples[i].startB,tuples[i].endB,orientationI,tuples[j].startA,tuples[j].endA,tuples[j].startB,tuples[j].endB,tuples[j].orientation,stack);
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> 6  success="+success);
					if (!success) continue;
					tmpAlignmentPair.node1=nodesArray[pairs[p].neighbor1];
					tmpAlignmentPair.node2=nodesArray[pairs[p].neighbor2];
					tmpAlignmentPair.alignmentStart1=stack[0];
					tmpAlignmentPair.alignmentEnd1=stack[1];
					tmpAlignmentPair.alignmentStart2=stack[2];
					tmpAlignmentPair.alignmentEnd2=stack[3];
					tmpAlignmentPair.orientation=stack[4]==1;
					tmpAlignmentPair.diffs=0;  // We do not try to estimate the diffs of such artificial same-read alignment.
					tmpAlignmentPair.inCanonicalForm();
					edge = new Edge(tmpAlignmentPair);
					edge.supplement=true;
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> 7  adding edge "+edge);
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> node1: "+nodesArray[edge.nodeID1]);
if (IO.SHOW_INTERACTIVE) System.err.println("addEdgesToOverlapNeighbors_sameRead> node2: "+nodesArray[edge.nodeID2]);
					k=findEdge(tmpAlignmentPair.node1.nodeID,edge,false);
					if (k==Math.NEGATIVE_INFINITY) {
						k=findEdge(tmpAlignmentPair.node1.nodeID,edge,true);
						if (k==Math.NEGATIVE_INFINITY) {
							addEdge(tmpAlignmentPair.node1.nodeID,edge);
							addEdge(tmpAlignmentPair.node2.nodeID,edge);
							if (edge.overlap!=-1 || edge.isOverlap_weak(IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2)) out=true;
							nAdded++;
						}
						else {
							previousEdge=neighbors[tmpAlignmentPair.node1.nodeID][k];
							if (!previousEdge.contains(edge,IDENTITY_THRESHOLD)) {
								previousEdge.addTypes(edge);
								previousEdge.addDiffs(edge);
								previousEdge.addOrientation(edge);
								if (edge.overlap!=-1 || edge.isOverlap_weak(IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2)) out=true;
								nAdded++;
							}
							else {
								// The supplement edge is already contained in a
								// supplement edge.
							}
						}
					}
					else {
						previousEdge=neighbors[tmpAlignmentPair.node1.nodeID][k];
						if (!previousEdge.contains(edge,IDENTITY_THRESHOLD)) {
							k=findEdge(tmpAlignmentPair.node1.nodeID,edge,true);
							if (k==Math.NEGATIVE_INFINITY) {
								addEdge(tmpAlignmentPair.node1.nodeID,edge);
								addEdge(tmpAlignmentPair.node2.nodeID,edge);
								if (edge.overlap!=-1 || edge.isOverlap_weak(IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2)) out=true;
								nAdded++;
							}
							else {
								previousEdge=neighbors[tmpAlignmentPair.node1.nodeID][k];
								if (!previousEdge.contains(edge,IDENTITY_THRESHOLD)) {
									previousEdge.addTypes(edge);
									previousEdge.addDiffs(edge);
									previousEdge.addOrientation(edge);
									if (edge.overlap!=-1 || edge.isOverlap_weak(IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2)) out=true;
									nAdded++;
								}
								else {
									// The supplement edge is already contained in a
									// supplement edge.
								}
							}
						}
						else {
							// The supplement edge is already contained in a regular edge
						}
					}
				}
				i++;
			}
			p++;
			if (firstIForNextP!=-1) i=firstIForNextP;
			firstIForNextP=-1;
		}
		tmpArray[0]=nAdded; tmpArray[1]=p;
		return out;
	}
	
	
	/**
	 * Representation of an alignment used by $addEdgesToOverlapNeighbors*$ procedures.
	 */
	private static class AlignmentTuple implements Comparable {
		public boolean orientation;
		public int readA, readB, startA, endA, startB, endB, diffs;
		
		public AlignmentTuple(boolean o, int ra, int rb, int sa, int ea, int sb, int eb, int d) {
			set(o,ra,rb,sa,ea,sb,eb,d);
		}
		
		public void set(boolean o, int ra, int rb, int sa, int ea, int sb, int eb, int d) {
			orientation=o;
			readA=ra; 
			readB=rb;
			startA=sa;
			endA=ea;
			startB=sb;
			endB=eb;
			diffs=d;
		}
		
		/**
		 * Sorts by $readA$, $readB$, $startA$.
		 */
		public int compareTo(Object other) {
			AlignmentTuple otherTuple = (AlignmentTuple)other;
			if (readA<otherTuple.readA) return -1;
			if (readA>otherTuple.readA) return 1;
			if (readB<otherTuple.readB) return -1;
			if (readB>otherTuple.readB) return 1;
			if (startA<otherTuple.startA) return -1;
			if (startA>otherTuple.startA) return 1;
			return 0;
		}
		
		public String toString() {
			return readA+"["+startA+".."+endA+"] x "+readB+"["+startB+".."+endB+"] orientation="+orientation+" diffs="+diffs;
		}
	}
	
	
	/**
	 * Builds the sorted list of all distinct pairs of overlap neighbors that use the same
	 * side (prefix or suffix) of the same overlap fork. These are the pairs of nodes 
	 * between which additional edges should be added by $addEdgesToOverlapNeighbors()$.
	 *
	 * Remark: nodes with field $visited=1$ are not considered as overlap forks.
	 *
	 * Remark: edges incident to periodic intervals are discarded, as well as edges that
	 * correspond to more than one alignment.
	 *
	 * Remark: the output list is not sorted by the ID of the fork. See class 
	 * $OverlapPair$ for details.
	 *
	 * Remark: the procedure uses $stack$ as temporary space.
	 *
	 * Remark: the procedure does not assume $neighbors[i]$ to be sorted for any $i$, and
	 * it considers all edges in $neighbors[i]$, including supplement edges.
	 *
	 * @param kernelsComputed if TRUE, assumes that kernel tags have already been assigned
	 * to nodes;
	 * @param out output array: 0=last element in the output (-1 if the output is null);
	 * @return NULL if no node has more than one overlap neighbor on the same side.
	 */
	private static final OverlapPair[] addEdgesToOverlapNeighbors_getPairs(boolean kernelsComputed, int[] out) {
		final int GROWTH_RATE = 100;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_OVERHANG_LENGTH = IO.quantum<<1;  // Arbitrary
		int i, j, k, h, p, q;
		int nOverlapsPrefix, nOverlapsSuffix, nPairs, tmpArrayLength, overhang, to, nAlignments;
		int[] tmpArray1 = new int[GROWTH_RATE];
		int[] tmpArray2 = new int[GROWTH_RATE];
		OverlapPair tmpPair;
		OverlapPair[] pairs;
		
		// Counting the number of pairs
		tmpArrayLength=0; nPairs=0;
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].visited==1 || nodesArray[i].type==Constants.INTERVAL_PERIODIC) continue;
			if (kernelsComputed?nodesArray[i].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[i])) continue;
			nOverlapsPrefix=0; nOverlapsSuffix=0;
			for (j=0; j<nNeighbors[i]; j++) {
				to=neighbors[i][j].getTo(i);
				if (nodesArray[to].type==Constants.INTERVAL_PERIODIC) continue;
				if (kernelsComputed?nodesArray[to].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[to])) continue;
				nAlignments=(neighbors[i][j].getNAlignments())<<1;
				if (nAlignments>tmpArray1.length) {
					tmpArray1 = new int[nAlignments+GROWTH_RATE];
					tmpArray2 = new int[nAlignments+GROWTH_RATE];
				}
				overhang=neighbors[i][j].getOverhang(i,true);
				if (overhang==Math.POSITIVE_INFINITY) overhang=neighbors[i][j].getOverhang_weak(i,true,IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2);
				if (overhang!=Math.POSITIVE_INFINITY && Math.abs(overhang)>=MIN_OVERHANG_LENGTH) nOverlapsPrefix++;
				overhang=neighbors[i][j].getOverhang(i,false);
				if (overhang==Math.POSITIVE_INFINITY) overhang=neighbors[i][j].getOverhang_weak(i,false,IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2);
				if (overhang!=Math.POSITIVE_INFINITY && Math.abs(overhang)>=MIN_OVERHANG_LENGTH) nOverlapsPrefix++;
			}
			if (nOverlapsPrefix<=1 && nOverlapsSuffix<=1) continue;
			if (nOverlapsPrefix>tmpArrayLength) tmpArrayLength=nOverlapsPrefix;
			if (nOverlapsSuffix>tmpArrayLength) tmpArrayLength=nOverlapsSuffix;
			if (nOverlapsPrefix>1) nPairs+=((nOverlapsPrefix-1)*nOverlapsPrefix)>>1;
			if (nOverlapsSuffix>1) nPairs+=((nOverlapsSuffix-1)*nOverlapsSuffix)>>1;
		}
		if (nPairs==0) {
			out[0]=-1;
			return null;
		}
		if (stack.length<(tmpArrayLength<<1)) stack = new int[tmpArrayLength<<1];
		pairs = new OverlapPair[nPairs];
		
		// Loading pairs
		k=-1;
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].visited==1 || nodesArray[i].type==Constants.INTERVAL_PERIODIC) continue;
			if (kernelsComputed?nodesArray[i].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[i])) continue;
			// Prefix side
			h=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				to=neighbors[i][j].getTo(i);
				if (nodesArray[to].type==Constants.INTERVAL_PERIODIC) continue;
				if (kernelsComputed?nodesArray[to].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[to])) continue;
				overhang=neighbors[i][j].getOverhang(i,true);
				if (overhang==Math.POSITIVE_INFINITY) overhang=neighbors[i][j].getOverhang_weak(i,true,IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2);
				if (overhang==Math.POSITIVE_INFINITY || Math.abs(overhang)<MIN_OVERHANG_LENGTH) continue;
				stack[++h]=to; stack[++h]=overhang;
			}
			if (h>1) {
				for (p=0; p<h; p+=2) {
					if (stack[p]==i) continue;
					for (q=p+2; q<h; q+=2) {
						if (stack[q]==stack[p]) continue;
						k++;
						pairs[k] = new OverlapPair(true,i,stack[p],stack[p+1],stack[q],stack[q+1]);
						pairs[k].inCanonicalForm();
					}
				}
				nodesArray[i].visited=1;
			}
			// Suffix side
			h=-1;
			for (j=0; j<nNeighbors[i]; j++) {
				to=neighbors[i][j].getTo(i);
				if (nodesArray[to].type==Constants.INTERVAL_PERIODIC) continue;
				if (kernelsComputed?nodesArray[to].lastKernel==-1:!IntervalGraphStep3.canBeKernel(nodesArray[to])) continue;
				overhang=neighbors[i][j].getOverhang(i,false);
				if (overhang==Math.POSITIVE_INFINITY) overhang=neighbors[i][j].getOverhang_weak(i,false,IntervalGraphStep3.NEW_OVERLAPS_THRESHOLD,tmpArray1,tmpArray2);
				if (overhang==Math.POSITIVE_INFINITY || Math.abs(overhang)<MIN_OVERHANG_LENGTH) continue;
				stack[++h]=to; stack[++h]=overhang;
			}
			if (h>1) {
				for (p=0; p<h; p+=2) {
					if (stack[p]==i) continue;
					for (q=p+2; q<h; q+=2) {
						if (stack[q]==stack[p]) continue;
						k++;
						pairs[k] = new OverlapPair(false,i,stack[p],stack[p+1],stack[q],stack[q+1]);
						pairs[k].inCanonicalForm();
					}
				}
				nodesArray[i].visited=1;
			}
		}
		
		// Removing approximately identical pairs
		if (k>0) {
			OverlapPair.order=OverlapPair.ORDER_DISTINCT;
			Arrays.sort(pairs,0,k+1);
			j=0;
			for (i=1; i<=k; i++) {
				if ( pairs[i].neighbor1==pairs[j].neighbor2 && pairs[i].neighbor2==pairs[j].neighbor2 && 
					 ( (pairs[i].overhang1>=0 && pairs[j].overhang1>=0 && Math.abs(pairs[i].overhang1,pairs[j].overhang1)<=IDENTITY_THRESHOLD) ||
					   (pairs[i].overhang1<0 && pairs[j].overhang1<0 && Math.abs(-pairs[i].overhang1,-pairs[j].overhang1)<=IDENTITY_THRESHOLD)
					 ) &&
  					 ( (pairs[i].overhang2>=0 && pairs[j].overhang2>=0 && Math.abs(pairs[i].overhang2,pairs[j].overhang2)<=IDENTITY_THRESHOLD) ||
  					   (pairs[i].overhang2<0 && pairs[j].overhang2<0 && Math.abs(-pairs[i].overhang2,-pairs[j].overhang2)<=IDENTITY_THRESHOLD)
  					 ) 
				   ) continue;
				j++;
				tmpPair=pairs[j];
				pairs[j]=pairs[i];
				pairs[i]=tmpPair;
			}
			k=j;
			if (k>0) {
				OverlapPair.order=OverlapPair.ORDER_READS;
				Arrays.sort(pairs,0,k+1);
			}
		}
		System.err.println((k+1)+" distinct pairs of overlap neighbors loaded");
		out[0]=k;
		return pairs;
	}
	
	
	/**
	 * A pair of neighbors of an overlap fork $v$, that overlap $v$ on the same side of 
	 * $v$. The neighbors can be on the same read.
	 */
	private static class OverlapPair implements Comparable {
		public static final int ORDER_UNSORTED = 0;
		public static final int ORDER_DISTINCT = 1;
		public static final int ORDER_READS = 2;
		
		public static int order;
		
		public boolean side;  // TRUE=prefix side of the fork; FALSE=suffix side of the fork.
		public int from;  // ID of the overlap fork
		public int neighbor1, neighbor2;  // neighbor1: the one with smaller read; if the two reads are the same, the one with smaller first position of its overhang.
		public int overhang1, overhang2;  // positive (respectively, negative) if $overhangX$ is a prefix (respectively, suffix) of $neighborX$.
		public int read1, read2;
		
		
		public String toString() {
			return "neighbor1: "+nodesArray[neighbor1]+"\n neighbor2: "+nodesArray[neighbor2]+"\n overhang1="+overhang1+" overhang2="+overhang2;
		}
		
		
		public OverlapPair(boolean s, int f, int n1, int o1, int n2, int o2) {
			side=s;
			from=f;
			neighbor1=n1; neighbor2=n2;
			read1=nodesArray[neighbor1].read;
			read2=nodesArray[neighbor2].read;
			overhang1=o1; overhang2=o2;
		}
		
		
		/**
		 * @param neighborOneOrTwo TRUE=1, FALSE=2;
		 * @return the absolute first position of the overhang in the read of the selected 
		 * neighbor.
		 */
		public final int overhangStart(boolean neighborOneOrTwo) {
			if (neighborOneOrTwo) return overhang1>0?nodesArray[neighbor1].start:nodesArray[neighbor1].end+overhang1;
			else return overhang2>0?nodesArray[neighbor2].start:nodesArray[neighbor2].end+overhang2;
		}
		
		
		/**
		 * @param neighborOneOrTwo TRUE=1, FALSE=2;
		 * @return the absolute last position of the overhang in the read of the selected
		 * neighbor.
		 */
		public final int overhangEnd(boolean neighborOneOrTwo) {
			if (neighborOneOrTwo) return overhang1>0?nodesArray[neighbor1].start+overhang1-1:nodesArray[neighbor1].end;
			else return overhang2>0?nodesArray[neighbor2].start+overhang2-1:nodesArray[neighbor2].end;
		}
		
		
		/**
		 * Assigns the neighbor with smaller readID to *1 variables. If the neighbors have
		 * the same readID, *1 variables are assigned to the neighbor whose overhang has 
		 * the smaller starting position in the read.
		 */
		public final void inCanonicalForm() {			
			int tmp;
			
			if (read1<read2 || (read1==read2 && overhangStart(true)<=overhangStart(false))) return;
			tmp=read1; read1=read2; read2=tmp;
			tmp=neighbor1; neighbor1=neighbor2; neighbor2=tmp;
			tmp=overhang1; overhang1=overhang2; overhang2=tmp;
		}
		
		
		/**
		 * If $order=ORDER_READS$, puts all pairs in which $read1=read2$ in a consecutive 
		 * block at the end of the list; then, sorts each block by $read1$, $read2$, and 
		 * first position of the overhang of $neighbor1$.
		 *
		 * If $order=ORDER_DISTINCT$, sorts by $neighbor1,neighbor2,overhang1,overhang2$.
		 */
		public int compareTo(Object other) {
			OverlapPair otherPair = (OverlapPair)other;
			if (order==ORDER_READS) {
				if (read1!=read2 && otherPair.read1==otherPair.read2) return -1;
				if (read1==read2 && otherPair.read1!=otherPair.read2) return 1;
				if (read1<otherPair.read1) return -1;
				if (read1>otherPair.read1) return 1;
				if (read2<otherPair.read2) return -1;
				if (read2>otherPair.read2) return 1;
				final int start1 = overhangStart(true);
				final int otherStart1 = otherPair.overhangStart(true);
				if (start1<otherStart1) return -1;
				if (start1>otherStart1) return 1;
			}
			else if (order==ORDER_DISTINCT) {
				if (neighbor1<otherPair.neighbor1) return -1;
				if (neighbor1>otherPair.neighbor1) return 1;
				if (neighbor2<otherPair.neighbor2) return -1;
				if (neighbor2>otherPair.neighbor2) return 1;
				if (overhang1<otherPair.overhang1) return -1;
				if (overhang1>otherPair.overhang1) return 1;
				if (overhang2<otherPair.overhang2) return -1;
				if (overhang2>otherPair.overhang2) return 1;
			}
			return 0;
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


	
	// --------------------------------- DATA STRUCTURES ---------------------------------
	
	/**
	 * Represents the pairing between two intervals in distinct reads induced by an 
	 * alignment that is assigned to both.
	 */
	public static class AlignmentPair implements Comparable {
		public Node node1;  // The node with smaller value of $read$
		public Node node2;
		
		// Properties of the alignment
		public int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;  // In the forward orientation of both reads
		public boolean orientation;
		public double diffs;
		public boolean implied;  // From $Factorize.writeIntervalConnectionFile$

		
		public AlignmentPair(Node n1, Node n2, int s1, int e1, int s2, int e2, boolean o, double d, boolean m) {
			node1=n1;
			node2=n2;
			alignmentStart1=s1;
			alignmentEnd1=e1;
			alignmentStart2=s2;
			alignmentEnd2=e2;
			orientation=o;
			diffs=2.0*d/(e1-s1+e2-s2+2);
			implied=m;
		}
		
		
		public AlignmentPair() {
			node1=null;
			node2=null;
			alignmentStart1=-1;
			alignmentEnd1=-1;
			alignmentStart2=-1;
			alignmentEnd2=-1;
			orientation=false;
			diffs=-1.0;
			implied=false;
		}
		
		
		public int alignmentLength1() {
			return alignmentEnd1-alignmentStart1;
		}
		
		
		public int alignmentLength2() {
			return alignmentEnd2-alignmentStart2;
		}
		
		
		public String toString() {
			return (node1==null?"null":node1.read+",*"+node1.nodeID+"*["+node1.start+" ("+alignmentStart1+".."+alignmentEnd1+") "+node1.end+"]")+" x "+(node2==null?"null":node2.read+",*"+node2.nodeID+"*["+node2.start+" ("+alignmentStart2+".."+alignmentEnd2+") "+node2.end+"] "+(node1==null?"--":node1.type)+"::"+(node2==null?"--":node2.type))+" -- orientation="+orientation+" diffs="+diffs;
		}
		
		
		/**
		 * Creates a copy of this alignment pair, with $node1$ and $node1$ swapped.
		 * This violates the definition of $node1$ and $node2$.
		 */
		public AlignmentPair symmetrize() {
			return new AlignmentPair(node2,node1,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,orientation,diffs,implied);
		}
		
		
		/**
		 * Uses just the order between $node1$ fields.
		 */
		public int compareTo(Object other) {
			AlignmentPair otherPair = (AlignmentPair)other;
			return node1.compareTo(otherPair.node1);
		}
		
		
		/**
		 * Uses just equality between $node1$ fields.
		 */
		public boolean equals(Object other) {
			AlignmentPair otherPair = (AlignmentPair)other;
			return node1.equals(otherPair.node1);
		}
		
		
		/**
		 * Forces $alignmentStart*,alignmentEnd*$ to be inside the intervals of the
		 * corresponding nodes.
		 */
		public void trim() {
			alignmentStart1=Math.max(node1.start,alignmentStart1);
			alignmentEnd1=Math.min(node1.end,alignmentEnd1);
			alignmentStart2=Math.max(node2.start,alignmentStart2);
			alignmentEnd2=Math.min(node2.end,alignmentEnd2);
		}
		
		
		/**
		 * Enforces $node1.read<=node2.read$.
		 */
		public void inCanonicalForm() {
			if (node1.read<=node2.read) return;
			int tmpAlignmentStart1, tmpAlignmentEnd1;
			Node tmpNode;
			
			tmpNode=node1;
			tmpAlignmentStart1=alignmentStart1;
			tmpAlignmentEnd1=alignmentEnd1;
			node1=node2;
			alignmentStart1=alignmentStart2;
			alignmentEnd1=alignmentEnd2;
			node2=tmpNode; tmpNode=null;
			alignmentStart2=tmpAlignmentStart1;
			alignmentEnd2=tmpAlignmentEnd1;
		}
		
	}
	
	
	/**
	 * A node of the interval graph, i.e. a specific repeat interval of a specific read.
	 */
	public static class Node implements Comparable {
		public static final int UNSORTED = 0;
		public static final int START = 1;
		public static final int NODE_ID = 2;
		public static final int READ_START_NODEID = 3;
		public static final int KERNEL_BIDIRECTED = 4;
		public static final int BIDIRECTED = 5;
		public static int order;
		public static int hashMode = 0;  // 0=read,type,ID. 1=read,start,end.
		
		/**
		 * Properties of the interval
		 */
		public int read;  // Zero-based
		public int type, id;
		public int start, end;
		public boolean isLeftMaximal, isRightMaximal;
		public boolean isWeak;
		public int isContainedInterval;  // See $Factorize.markContainedIntervals()$ for details.
		public int period;
		public int[] periodicSubstrings;  // For each periodic substring: first and last position relative to $start$, and 1/0 flag for long period.
		public int lastPeriodicSubstring;
		public int minPrefixLength, minSuffixLength;
		public boolean hasLongPeriod;  // Just tells if the period is long or short, i.e. it is not necessarily equal to the corresponding flag produced during factorization.
		public int firstMaximalStart, lastMaximalEnd;
		public int nAssignedAlignments;
		public double avgCoverage;
		
		/**
		 * Properties of the interval graph
		 */
		public int nodeID;
		public boolean isContained, contains;
		public boolean hasRepeat;  // There is a repeat in the sequence
		public boolean isRCPalindrome;
		public int[] components;
		public int[] kernels;  // An interval can belong to multiple kernels simultaneously. All such IDs are non-negative.
		public byte[] kernelOrientations;  // 0=forward; 1=RC; 2=both.
		public int lastKernel, lastComponent;
		public int bidirectedGraphNode;  // ID of the bidirected graph node this node is projected to (-1 if it is not projected).
		public byte bidirectedGraphNodeOrientation;  // 0=forward; 1=RC; 2=both.
		public boolean isFullLength;  // TRUE iff the interval is a full-length copy of its kernel repeat
		public boolean frontier;
		public int lastPathWithStart, lastPathWithEnd;
		public int[] pathsWithStart, pathsWithEnd;  // Positive ($x$): beginning of kernel path $x$. Negative ($-1-x$): end of kernel path $x$. Both $x$ and $-1-x$ could occur in each array.
		
		/**
		 * Temporary space
		 */
		public int visited, iteration;
		public boolean isMaximal, isMinimal;
		public boolean isSingletonComponent;
		
		
		public final void set(int r, int t, int i, int s, int e, boolean lm, boolean rm, boolean w, int c, int mpl, int msl, int p, boolean hlp, int fms, int lme, int naa, double avgc) {
			// Properties of the interval
			read=r; type=t; id=i;
			start=s; end=e;
			isLeftMaximal=lm; isRightMaximal=rm;
			isWeak=w; isContainedInterval=c;
			period=p;
			lastPeriodicSubstring=-1;
			minPrefixLength=mpl;
			minSuffixLength=msl;
			if (period>0) hasLongPeriod=period>=PeriodicSubstrings.MIN_PERIOD_LONG;
			else hasLongPeriod=hlp;
			firstMaximalStart=fms;
			lastMaximalEnd=lme;
			nAssignedAlignments=naa;
			avgCoverage=avgc;
			
			// Properties of the interval graph
			nodeID=-1;
			isContained=false; contains=false;
			hasRepeat=false;
			isRCPalindrome=false;
			lastKernel=-1; lastComponent=-1;
			lastPathWithStart=-1; lastPathWithEnd=-1;
			pathsWithStart=null; pathsWithEnd=null;
			isFullLength=false; frontier=false;
			bidirectedGraphNode=-1; bidirectedGraphNodeOrientation=-1;
			kernelOrientations=null;
			
			// Temporary space
			visited=-1; isMaximal=false; isMinimal=false; isSingletonComponent=false;
		}
		
		
		public Node(int r, int t, int i, int s, int e, boolean lm, boolean rm, boolean w, int c, int mpl, int msl, int p, boolean hlp, int fms, int lme, int naa, double avgc) {
			set(r,t,i,s,e,lm,rm,w,c,mpl,msl,p,hlp,fms,lme,naa,avgc);
		}
		
		
		public Node() {
			clear();
		}
		
		
		/**
		 * Copies onto $otherNode$ all attributes of this node except for $nodeID$, which
		 * is set to -1.
		 */
		public final void clone(Node otherNode) {
			// Properties of the interval
			otherNode.read=read;
			otherNode.type=type;
			otherNode.id=id;
			otherNode.start=start;
			otherNode.end=end;
			otherNode.isLeftMaximal=isLeftMaximal;
			otherNode.isRightMaximal=isRightMaximal;
			otherNode.isWeak=isWeak;
			otherNode.isContainedInterval=isContainedInterval;
			otherNode.period=period;
			otherNode.lastPeriodicSubstring=lastPeriodicSubstring;
			if (lastPeriodicSubstring>=0) {
				otherNode.periodicSubstrings = new int[lastPeriodicSubstring+1];
				System.arraycopy(periodicSubstrings,0,otherNode.periodicSubstrings,0,lastPeriodicSubstring+1);
			}
			otherNode.minPrefixLength=minPrefixLength;
			otherNode.minSuffixLength=minSuffixLength;
			otherNode.hasLongPeriod=hasLongPeriod;
			otherNode.firstMaximalStart=firstMaximalStart;
			otherNode.lastMaximalEnd=lastMaximalEnd;
			
			// Properties of the interval graph
			otherNode.nodeID=-1;
			otherNode.isContained=isContained;
			otherNode.contains=contains;
			otherNode.hasRepeat=hasRepeat;
			otherNode.isRCPalindrome=isRCPalindrome;
			otherNode.lastComponent=lastComponent;
			if (lastComponent>=0) {
				otherNode.components = new int[lastComponent+1];
				System.arraycopy(components,0,otherNode.components,0,lastComponent+1);
			}
			otherNode.lastKernel=lastKernel;
			if (lastKernel>=0) {
				otherNode.kernels = new int[lastKernel+1];
				System.arraycopy(kernels,0,otherNode.kernels,0,lastKernel+1);
				if (kernelOrientations!=null) {
					otherNode.kernelOrientations = new byte[lastKernel+1];
					System.arraycopy(kernelOrientations,0,otherNode.kernelOrientations,0,lastKernel+1);
				}
				else otherNode.kernelOrientations=null;
			}
			else otherNode.kernelOrientations=null;
			otherNode.bidirectedGraphNode=bidirectedGraphNode; 
			otherNode.bidirectedGraphNodeOrientation=bidirectedGraphNodeOrientation;
			otherNode.lastPathWithStart=lastPathWithStart;
			if (lastPathWithStart>=0) {
				otherNode.pathsWithStart = new int[lastPathWithStart+1];
				System.arraycopy(pathsWithStart,0,otherNode.pathsWithStart,0,lastPathWithStart+1);
			}
			else otherNode.pathsWithStart=null;
			otherNode.lastPathWithEnd=lastPathWithEnd;
			if (lastPathWithEnd>=0) {
				otherNode.pathsWithEnd = new int[lastPathWithEnd+1];
				System.arraycopy(pathsWithEnd,0,otherNode.pathsWithEnd,0,lastPathWithEnd+1);
			}
			else otherNode.pathsWithEnd=null;
			otherNode.isFullLength=isFullLength;
			otherNode.frontier=frontier;
		}
		
		
		public final void clear() {
			// Properties of the interval
			read=-1;
			type=-1;
			id=-1;
			start=-1;
			end=-1;
			isLeftMaximal=false;
			isRightMaximal=false;
			isWeak=false;
			isContainedInterval=0;
			period=-1;
			lastPeriodicSubstring=-1;
			minPrefixLength=-1;
			minSuffixLength=-1;
			hasLongPeriod=false;
			firstMaximalStart=-1;
			lastMaximalEnd=-1;
			nAssignedAlignments=0;
			avgCoverage=0.0;
			
			// Properties of the interval graph
			nodeID=-1;
			isContained=false;
			contains=false;
			hasRepeat=false;
			isRCPalindrome=false;
			lastComponent=-1;
			lastKernel=-1;
			kernelOrientations=null;
			bidirectedGraphNode=-1; bidirectedGraphNodeOrientation=-1;
			lastPathWithStart=-1; lastPathWithEnd=-1;
			pathsWithStart=null; pathsWithEnd=null;
			isFullLength=false;
			frontier=false;
		}
		
		
		public final boolean isRepresentative() {
			return !isContained && isLeftMaximal && isRightMaximal;
		}
		
		
		/**
		 * Remark: making $Node$ implement $Comparable<Node>$ instead of just $Comparable$
		 * creates problems with the global $nodes$ hashmap: it returns null when queried 
		 * with an element that is already in the map, even though $equals()$ and 
		 * $hashCode()$ are implemented correctly.
		 */
		public final boolean equals(Object other) {
			Node otherNode = (Node)other;
			if (hashMode==0) return otherNode.read==read && otherNode.type==type && otherNode.id==id;
			else return otherNode.read==read && otherNode.start==start && otherNode.end==end && otherNode.type==type;
		}
		
		
		public int hashCode() {
			if (hashMode==0) return (read+"-"+type+"-"+id).hashCode();
			else return (read+"-"+start+"-"+end+"-"+type).hashCode();
		}
		
		
		public final int length() {
			return end-start+1;
		}
		
		
		public int compareTo(Object other) {
			Node otherNode = (Node)other;
			if (order==START) {
				if (start<otherNode.start) return -1;
				if (start>otherNode.start) return 1;
			}
			else if (order==NODE_ID) {
				if (nodeID<otherNode.nodeID) return -1;
				if (nodeID>otherNode.nodeID) return 1;
			}
			else if (order==READ_START_NODEID) {
				if (read<otherNode.read) return -1;
				if (read>otherNode.read) return 1;
				if (start<otherNode.start) return -1;
				if (start>otherNode.start) return 1;
				if (nodeID<otherNode.nodeID) return -1;
				if (nodeID>otherNode.nodeID) return 1;
			}
			else if (order==KERNEL_BIDIRECTED) {  // Assumes that both nodes have $lastKernel=0$.
				if (kernels[0]<otherNode.kernels[0]) return -1;
				else if (kernels[0]>otherNode.kernels[0]) return 1;
				if (bidirectedGraphNode<otherNode.bidirectedGraphNode) return -1;
				else if (bidirectedGraphNode>otherNode.bidirectedGraphNode) return 1;
			}
			else if (order==BIDIRECTED) {
				if (bidirectedGraphNode<otherNode.bidirectedGraphNode) return -1;
				else if (bidirectedGraphNode>otherNode.bidirectedGraphNode) return 1;
			}
			return 0;
		}

		
		/**
		 * @param clusterField used only if not null;
		 * @param pos position in $clusterField$; used only if the latter is not null;
		 * @param printedComponents sorted
		 */
		public final String toDot(int[] printedComponents, int[] clusterField, int pos /*int[] kernel2family*/) {
//			return nodeID+" [label=\""+length()+"-"+read+"\"];\n";


/*			if (kernels==null || kernels.length==0 || lastKernel<0) return nodeID+" [kernel=-1];\n";
			if (lastKernel>0) return nodeID+" [kernel=-2];\n";
			return nodeID+" [kernel="+kernels[0]+"];\n";
*/			

/*			
			if (components==null || components.length==0 || lastComponent<0) return nodeID+" [component=-1];\n";
			if (lastComponent>0) return nodeID+" [component=-2];\n";
			return nodeID+" [component="+components[0]+"];\n";
*/			
			
			
//			return nodeID+" [period="+period+"];\n";


//			return nodeID+(clusterField!=null?" [cluster="+clusterField[pos]+"];\n":" ;\n");
			
			
//			return nodeID+" [family="+(lastKernel==0?kernel2family[kernels[lastKernel]]:"-1")+"];\n";
			
			
			
//			return nodeID+" [bc="+betweennessCentrality+"];\n";
			
/*			int tag = 0;
			if (isMaximal && !isMinimal) tag=1;
			else if (!isMaximal && isMinimal) tag=2;
			return nodeID+" [tag="+tag+"];\n";
*/			
			
/*			
			if (components==null || components.length==0 || lastComponent<0) return nodeID+" [component=-1];\n";
			if (lastComponent>0) return nodeID+" [component=-2];\n";
			if (printedComponents!=null && printedComponents.length!=0 && Arrays.binarySearch(printedComponents,components[0])<0) return nodeID+" [component=-1];\n";
			return nodeID+" [component="+components[0]+"];\n";
*/

			
			int component;
			if (components==null || components.length==0 || lastComponent<0) component=-1;
			else if (lastComponent>0) component=-2;
			else component=components[0];
			int avgCoverageBin = -1;
			if (avgCoverage<10*IO.coverage) avgCoverageBin=1;
			else if (avgCoverage<20*IO.coverage) avgCoverageBin=2;
			else if (avgCoverage<40*IO.coverage) avgCoverageBin=3;
			else if (avgCoverage<80*IO.coverage) avgCoverageBin=4;
			else avgCoverageBin=5;
			int cluster;
			if (lastComponent<0 || (lastComponent==0 && components[0]==-1)) cluster=0;
			else if (lastComponent>0) cluster=1;
			else cluster=2+components[lastComponent];
			return nodeID+" [component="+component+",type="+type+",isContained="+(isContained?"1":"0")+",isLeftMaximal="+(isLeftMaximal?"1":"0")+",isRightMaximal="+(isRightMaximal?"1":"0")+",nKernels="+(lastKernel+1)+",read="+read+",start="+start+",end="+end+",period="+period+",weak="+isWeak+",nAssignedAlignments="+nAssignedAlignments+",avgCoverageBin="+avgCoverageBin+",avgCoverage="+IO.format(avgCoverage)+",hasLongPeriod="+hasLongPeriod+",cluster="+cluster+"];\n";
		}
		
		
		/**
		 * Remark: the procedure does not add $[first..last]$ if there is already a 
		 * substring of the same type that contains it. Otherwise, the procedure merges 
		 * $[first..last]$ with all existing substrings of the same type that are 
		 * approximately identical to it or contained in it.
		 *
		 * Remark: $periodicSubstrings$ is not assumed to be sorted.
		 *
		 * Remark: the procedure can reallocate $tmpBoolean$.
		 *
		 * @param first,last relative to $start$.
		 */
		public final void addPeriodicSubstring(int first, int last, boolean hasLongPeriod) {
			final int TYPE = hasLongPeriod?1:0;
			boolean found;
			int i, j;
			int newFirst, newLast, newLength;			
			
			if (lastPeriodicSubstring==-1) {
				if (periodicSubstrings==null) periodicSubstrings = new int[(MIN_PERIODIC_SUBSTRINGS_PER_INTERVAL)*3];
				periodicSubstrings[0]=first;
				periodicSubstrings[1]=last;
				periodicSubstrings[2]=TYPE;
				lastPeriodicSubstring=2;
				return;
			}
			
			// Adding $[newFirst..newLast]$ iff no existing substring of the same type 
			// contains it, and merging it to the first existing container of the same
			// type otherwise.
			found=false;
			for (i=0; i<=lastPeriodicSubstring; i+=3) {
				if (periodicSubstrings[i+2]==TYPE && Intervals.isApproximatelyContained(first,last,periodicSubstrings[i],periodicSubstrings[i+1])) {
					found=true;
					periodicSubstrings[i]=Math.min(periodicSubstrings[i],first);
					periodicSubstrings[i+1]=Math.max(periodicSubstrings[i+1],last);
					break;
				}
			}
			if (found) return;
			
			// Removing substrings of the same type that are identical to or contained in
			// $[first..last]$, and computing their union with $[first..last]$.
			newLength=(lastPeriodicSubstring+1)/3;
			if (newLength>tmpBoolean.length) tmpBoolean = new boolean[newLength];
			Math.set(tmpBoolean,newLength-1,false);
			found=false; newFirst=first; newLast=last;
			for (i=0; i<=lastPeriodicSubstring; i+=3) {
				if ( periodicSubstrings[i+2]==TYPE &&
					 ( Intervals.areApproximatelyIdentical(periodicSubstrings[i],periodicSubstrings[i+1],first,last) ||
					   Intervals.isApproximatelyContained(periodicSubstrings[i],periodicSubstrings[i+1],first,last)
				     )
				   ) {
					found=true; tmpBoolean[i/3]=true;
					newFirst=Math.min(periodicSubstrings[i],newFirst);
					newLast=Math.max(periodicSubstrings[i+1],newLast);
				}
			}
			if (found) {
				j=-1;
				for (i=0; i<=lastPeriodicSubstring; i+=3) {
					if (tmpBoolean[i/3]) continue;
					periodicSubstrings[++j]=periodicSubstrings[i];
					periodicSubstrings[++j]=periodicSubstrings[i+1];
					periodicSubstrings[++j]=periodicSubstrings[i+2];
				}
				lastPeriodicSubstring=j;
			}
			if (lastPeriodicSubstring+3>=periodicSubstrings.length) {
				int[] newPeriodicSubstrings = new int[periodicSubstrings.length<<1];
				System.arraycopy(periodicSubstrings,0,newPeriodicSubstrings,0,periodicSubstrings.length);
				periodicSubstrings=newPeriodicSubstrings;
			}
			periodicSubstrings[++lastPeriodicSubstring]=newFirst;
			periodicSubstrings[++lastPeriodicSubstring]=newLast;
			periodicSubstrings[++lastPeriodicSubstring]=TYPE;
		}
		
		
		/**
		 * Remark: a long-period periodic interval is returned only if there is no short-
		 * period periodic interval that contains $[start..end]$. This is because short-
		 * period periodic intervals are less likely to straddle, whereas long-period 
		 * periodic intervals are more likely to straddle. In this case, $[start..end]$ 
		 * might belong to the union of two straddling long-period intervals while not 
		 * belonging to any long-period interval.
		 *
		 * Remark: $periodicSubstrings$ is not assumed to be sorted.
		 *
		 * @param first,last relative to $start$;
		 * @param tmp temporary space, of size at least equal to the number of periodic
		 * substrings of the node;
		 * @return 0=not in a periodic substring; 1=in a periodic substring with short
		 * period; 2=in a periodic substring with long period.
		 */
		public final int inPeriodicSubstring(int first, int last) {
			int i;
			int length, lengthShort, lengthLong;
			
			lengthShort=Math.POSITIVE_INFINITY; lengthLong=Math.POSITIVE_INFINITY;
			if (lastPeriodicSubstring==-1) return 0;
			for (i=0; i<=lastPeriodicSubstring; i+=3) {
				if ( Intervals.isApproximatelyContained(first,last,periodicSubstrings[i],periodicSubstrings[i+1]) ||
					 Intervals.areApproximatelyIdentical(first,last,periodicSubstrings[i],periodicSubstrings[i+1])
				   ) {
				   length=periodicSubstrings[i+1]-periodicSubstrings[i]+1;
				   if (periodicSubstrings[i+2]==0 && length<lengthShort) lengthShort=length;
				   else if (periodicSubstrings[i+2]==1 && length<lengthLong) lengthLong=length;
				}
			}
			if (lengthShort!=Math.POSITIVE_INFINITY) return 1;
			if (lengthLong!=Math.POSITIVE_INFINITY) return 2;
			return 0;
		}
		
		
		/**
		 * Analogous to $inPeriodicSubstring()$.
		 */
		public final int intersectsPeriodicSubstring(int first, int last) {
			final double MIN_INTERSECTION_RATIO = 1.0/3;  // Arbitrary
			final int MIN_INTERSECTION_LENGTH = (int)((last-first+1)*MIN_INTERSECTION_RATIO);
			int i;
			int length, lengthShort, lengthLong;
			
			lengthShort=Math.POSITIVE_INFINITY; lengthLong=Math.POSITIVE_INFINITY;
			if (lastPeriodicSubstring==-1) return 0;
			for (i=0; i<=lastPeriodicSubstring; i+=3) {
				if (Intervals.intersectionLength(first,last,periodicSubstrings[i],periodicSubstrings[i+1])>=MIN_INTERSECTION_LENGTH) {
				   length=periodicSubstrings[i+1]-periodicSubstrings[i]+1;
				   if (periodicSubstrings[i+2]==0 && length<lengthShort) lengthShort=length;
				   else if (periodicSubstrings[i+2]==1 && length<lengthLong) lengthLong=length;
				}
			}
			if (lengthShort!=Math.POSITIVE_INFINITY) return 1;
			if (lengthLong!=Math.POSITIVE_INFINITY) return 2;
			return 0;
		}
		
		
		/**
		 * Remark: the procedure assumes that the elements of $periodicSubstrings$ do not
		 * overlap.
		 *
		 * @param tmpArray temporary space, of size at least equal to the twice the n. of
		 * periodic substrings;
		 * @return TRUE iff a substring of length $>=gapLength$ of the node does not
		 * belong to any periodic substring.
		 */
		public final boolean gapInPeriodicSubstrings(int gapLength, int[] tmpArray) {
			final int nodeLength = length();
			int i, j;
			
			if (lastPeriodicSubstring==-1) return nodeLength>=gapLength;
			j=-1;
			for (i=0; i<lastPeriodicSubstring; i+=3) {
				tmpArray[++j]=periodicSubstrings[i];
				tmpArray[++j]=periodicSubstrings[i+1];
			}
			Arrays.sort(tmpArray,0,j+1);
			if (tmpArray[0]>=gapLength) return true;
			for (i=2; i<j; i+=2) {
				if (tmpArray[i]-tmpArray[i-1]-1>=gapLength) return true;
			}
			return nodeLength-tmpArray[j]-1>=gapLength;
		}

		
		public final boolean hasPeriodicNeighbor() {
			for (int j=0; j<nNeighbors[nodeID]; j++) {
				if (nodesArray[neighbors[nodeID][j].getTo(nodeID)].type==Constants.INTERVAL_PERIODIC) return true;
			}
			return false;
		}
		
		
		/**
		 * Remark: the procedure does not assume $components$ to be sorted.
		 *
		 * @return the position of $component$ in $components$, or -1 if $component$ 
		 * cannot be found in $components$.
		 */
		public final int inComponent(int component) {
			if (components==null) return -1;
			for (int i=0; i<=lastComponent; i++) {
				if (components[i]==component) return i;
			}
			return -1;
		}
		
		
		public String toString() {
			String out;
			out="("+nodeID+") "+read+","+type+","+id+" ["+start+".."+end+"], "+isLeftMaximal+","+isRightMaximal+","+isWeak+","+isContainedInterval+",*"+hasLongPeriod+"*,"+period+","+lastPeriodicSubstring+","+minPrefixLength+","+minSuffixLength+" | "+nodeID+","+isContained+","+contains+","+hasRepeat+","+isRCPalindrome+","+isFullLength+","+frontier+" | "+lastComponent+",";
			for (int i=0; i<=lastComponent; i++) out+=components[i]+",";
			out+=" << "+lastKernel+",";
			for (int i=0; i<=lastKernel; i++) out+="("+kernels[i]+","+(kernelOrientations==null?"-":kernelOrientations[i])+"), ";
			out+=" >> ";
			out+="$$ "+lastPathWithStart+",";
			for (int i=0; i<=lastPathWithStart; i++) out+=pathsWithStart[i]+",";
			out+=" $$ ";
			out+="%% "+lastPathWithEnd+",";
			for (int i=0; i<=lastPathWithEnd; i++) out+=pathsWithEnd[i]+",";
			out+=" %%";
			out+=" | "+nAssignedAlignments+","+IO.format(avgCoverage)+" // "+bidirectedGraphNode;
			return out;
		}
		
		
		public final void serialize(BufferedOutputStream out) throws IOException {
			int i, c, mask;
			
			// Properties of the interval
			IO.writeInt(read,out);
			IO.writeInt(type,out);
			IO.writeInt(id,out);
			IO.writeInt(start,out);
			IO.writeInt(end,out);
			c=0x00000000; mask=0x00000001;
			if (isLeftMaximal) c|=mask;
			mask<<=1;
			if (isRightMaximal) c|=mask;
			mask<<=1;
			if (isWeak) c|=mask;
			mask<<=1;
			if (hasLongPeriod) c|=mask;
			IO.writeInt(c,out);
			IO.writeInt(isContainedInterval,out);
			IO.writeInt(period,out);
			IO.writeInt(lastPeriodicSubstring,out);
			if (lastPeriodicSubstring>=2) {
				for (i=0; i<=lastPeriodicSubstring; i+=3) {
					IO.writeInt(periodicSubstrings[i],out);
					IO.writeInt(periodicSubstrings[i+1],out);
					IO.writeInt(periodicSubstrings[i+2],out);
				}
			}
			IO.writeInt(minPrefixLength,out);
			IO.writeInt(minSuffixLength,out);
			IO.writeInt(nAssignedAlignments,out);
			IO.writeDouble(avgCoverage,out);
			
			// Properties of the interval graph
			IO.writeInt(nodeID,out);
			c=0x00000000; mask=0x00000001;
			if (isContained) c|=mask;
			mask<<=1;
			if (contains) c|=mask;
			mask<<=1;
			if (hasRepeat) c|=mask;
			mask<<=1;
			if (isRCPalindrome) c|=mask;
			mask<<=1;
			if (isFullLength) c|=mask;
			mask<<=1;
			if (frontier) c|=mask;
			IO.writeInt(c,out);
			IO.writeInt(lastComponent,out);
			for (i=0; i<=lastComponent; i++) IO.writeInt(components[i],out);
			IO.writeInt(lastKernel,out);
			for (i=0; i<=lastKernel; i++) IO.writeInt(kernels[i],out);
			IO.writeInt(lastPathWithStart,out);
			for (i=0; i<=lastPathWithStart; i++) IO.writeInt(pathsWithStart[i],out);
			IO.writeInt(lastPathWithEnd,out);
			for (i=0; i<=lastPathWithEnd; i++) IO.writeInt(pathsWithEnd[i],out);
		}
		
		
		public final void deserialize(BufferedInputStream in) throws IOException {
			int i, c, mask;
			
			// Properties of the interval
			read=IO.readInt(in);
			type=IO.readInt(in);
			id=IO.readInt(in);
			start=IO.readInt(in);
			end=IO.readInt(in);
			c=IO.readInt(in);
			mask=0x00000001;
			isLeftMaximal=(c&mask)!=0;
			mask<<=1;
			isRightMaximal=(c&mask)!=0;
			mask<<=1;
			isWeak=(c&mask)!=0;
			mask<<=1;
			hasLongPeriod=(c&mask)!=0;
			isContainedInterval=IO.readInt(in);
			period=IO.readInt(in);
			lastPeriodicSubstring=IO.readInt(in);
			if (lastPeriodicSubstring>=2) {
				if (periodicSubstrings==null || periodicSubstrings.length<lastPeriodicSubstring+1) periodicSubstrings = new int[lastPeriodicSubstring+1];
				for (i=0; i<=lastPeriodicSubstring; i+=3) {
					periodicSubstrings[i]=IO.readInt(in);
					periodicSubstrings[i+1]=IO.readInt(in);
					periodicSubstrings[i+2]=IO.readInt(in);
				}
			}
			minPrefixLength=IO.readInt(in);
			minSuffixLength=IO.readInt(in);
			nAssignedAlignments=IO.readInt(in);
			avgCoverage=IO.readDouble(in);
			
			// Properties of the interval graph
			nodeID=IO.readInt(in);
			c=IO.readInt(in);
			mask=0x00000001;
			isContained=(c&mask)!=0;
			mask<<=1;
			contains=(c&mask)!=0;
			mask<<=1;
			hasRepeat=(c&mask)!=0;
			mask<<=1;
			isRCPalindrome=(c&mask)!=0;
			mask<<=1;
			isFullLength=(c&mask)!=0;
			mask<<=1;
			frontier=(c&mask)!=0;
			lastComponent=IO.readInt(in);
			if (components==null || components.length<lastComponent+1) components = new int[lastComponent+1];
			for (i=0; i<=lastComponent; i++) components[i]=IO.readInt(in);
			lastKernel=IO.readInt(in);
			if (kernels==null || kernels.length<lastKernel+1) kernels = new int[lastKernel+1];
			for (i=0; i<=lastKernel; i++) kernels[i]=IO.readInt(in);
			lastPathWithStart=IO.readInt(in);
			if (lastPathWithStart>=0) {
				pathsWithStart = new int[lastPathWithStart+1];
				for (i=0; i<=lastPathWithStart; i++) pathsWithStart[i]=IO.readInt(in);
			}
			else pathsWithStart=null;
			lastPathWithEnd=IO.readInt(in);
			if (lastPathWithEnd>=0) {
				pathsWithEnd = new int[lastPathWithEnd+1];
				for (i=0; i<=lastPathWithEnd; i++) pathsWithEnd[i]=IO.readInt(in);
			}
			else pathsWithEnd=null;
		}
		
		
		/**
		 * @return TRUE iff the node has the same kernel IDs as $otherNode$, and if
		 * kernels are in compatible orientations.
		 */
		public final boolean sameKernels(IntervalGraph.Node otherNode) {
			if (lastKernel!=otherNode.lastKernel) return false;
			for (int i=0; i<=lastKernel; i++) {
				if ( kernels[i]!=otherNode.kernels[i] || 
					 ( (kernelOrientations[i]==0 && otherNode.kernelOrientations[i]==1) ||
					   (kernelOrientations[i]==1 && otherNode.kernelOrientations[i]==0)
					 )
				   ) return false;
			}
			return true;
		}
		
		
		/**
		 * @return 1=this node has some kernel IDs, and they are a proper subset of those 
		 * of $otherNode$ (taking orientations into account); 2=opposite case; 3=identical
		 * kernel IDs; 0=none of the above.
		 */
		public final int subsetOfKernels(IntervalGraph.Node otherNode) {
			final int nKernelsI = lastKernel+1;
			final int nKernelsJ = otherNode.lastKernel+1;
			if (nKernelsI==0 || nKernelsJ==0) return 0;
			boolean canBeSubsetI, canBeSubsetJ;
			int i, j;
			int intersection;
			
			canBeSubsetI=true; canBeSubsetJ=true;
			i=0; j=0; intersection=0;
			while (i<nKernelsI && j<nKernelsJ) {
				if (kernels[i]<otherNode.kernels[j]) {
					i++; 
					canBeSubsetI=false;
					if (!canBeSubsetJ) return 0;
					continue;
				}
				else if (kernels[i]>otherNode.kernels[j]) {
					j++;
					canBeSubsetJ=false;
					if (!canBeSubsetI) return 0;
					continue;
				}
				if (kernelOrientations[i]!=otherNode.kernelOrientations[j]) return 0;
				intersection++; i++; j++;
			}
			if (intersection==0) return 0;
			if (nKernelsI==nKernelsJ && intersection==nKernelsI) return 3;
			if (nKernelsI<nKernelsJ && intersection==nKernelsI) return 1;
			if (nKernelsJ<nKernelsI && intersection==nKernelsJ) return 2;
			return 0;
		}
		
		
		/**
		 * Removes from $kernels,kernelOrientations,pathsWithStart,pathsWithEnd$ every 
		 * element of $kernels$ that does not occur in $otherNode$.
		 *
		 * Remark: the procedure assumes all kernel arrays of both nodes to be sorted.
		 *
		 * @param tmp1,tmp2 of size at least twice the number of kernels in this node.
		 */
		public final void removeKernels(Node otherNode, int[] tmp1, int[] tmp2) {
			int i, j;
			int tmpValue, lastAbsent, lastAbsentPrime;
			
			if (lastKernel==-1) return;
			lastAbsent=Math.setMinus(kernels,lastKernel,otherNode.kernels,otherNode.lastKernel,tmp1);
			if (lastAbsent==-1) return;
			
			// Cleaning $kernels,kernelOrientations$.
			i=0; j=0;
			while (i<=lastKernel && j<=lastAbsent) {
				if (kernels[i]<tmp1[j]) {
					i++;
					continue;
				}
				else if (kernels[i]>tmp1[j]) {
					j++;
					continue;
				}
				else {
					kernels[i]=-1;
					i++; j++;
				}
			}
			j=-1;
			for (i=0; i<=lastKernel; i++) {
				if (kernels[i]==-1) continue;
				j++;
				kernels[j]=kernels[i];
				kernelOrientations[j]=kernelOrientations[i];
			}
			lastKernel=j;
			
			// Cleaning $pathsWithStart,pathsWithEnd$.
			if (lastPathWithStart>=0 || lastPathWithEnd>=0) {
				System.arraycopy(tmp1,0,tmp1,lastAbsent+1,lastAbsent+1);
				for (i=0; i<=lastAbsent; i++) {
					tmpValue=-tmp1[lastAbsent-i];
					tmp1[lastAbsent-i]=-tmp1[i];
					tmp1[i]=tmpValue;
				}
				lastAbsentPrime=(lastAbsent+1)<<1;
				if (lastPathWithStart>=0) {
					lastPathWithStart=Math.setMinus(pathsWithStart,lastPathWithStart,tmp1,lastAbsentPrime,tmp2);
					if (lastPathWithStart>=0) System.arraycopy(tmp2,0,pathsWithStart,0,lastPathWithStart+1);
				}
				if (lastPathWithEnd>=0) {
					lastPathWithEnd=Math.setMinus(pathsWithEnd,lastPathWithEnd,tmp1,lastAbsentPrime,tmp2);
					if (lastPathWithEnd>=0) System.arraycopy(tmp2,0,pathsWithEnd,0,lastPathWithEnd+1);
				}
			}
		}
		
	}
	
	
	/**
	 * An edge of the interval graph
	 */
	public static class Edge implements Comparable {
		public static final int UNSORTED = -1;
		public static int order = UNSORTED;
		private static int idGenerator = -1;
		
		public int nodeID1;  // The node with smaller value of $read$ (order inherited from $AlignmentPair$)
		public int length1;
		public int nodeID2;
		public int length2;
		public boolean on;  // ON/OFF switch
		public boolean supplement;  // Supplement edges are just a way to store additional information, and should not be considered as regular edges of the graph.
		
		/**
		 * Information for Step 2
		 */
		public int containment;
		public int containmentSubtype;  // Set for both strict containment and insertion.
		public int overlap;
		public int insertion;
		public int sharedSubstring;
		public boolean isRedundant;
		public int[][] overhangs;  // Row $i$ corresponds to bit $i$ in $overlap$; column $j$ corresponds to $nodeIDj$. When merging edges, this matrix stores minima. -1 values mean no overlap.
		public int[][] maxOverhangs;  // When merging edges, this matrix stores maxima. The matrix exists only if multiple alignments are collapsed into this edge, and if at least one difference in overhangs exists. The matrix is a copy of $overhangs$ everywhere there is no difference in overhangs.
		public int containmentOverhangLeft, containmentOverhangRight;  // Set for both strict containment and insertion.
		public int[] sharedSubstringOverhangs;  // $node1left,node1right,node2left,node2right$. If the edge represents multiple alignments, this contains just information about one of them. Both $nodeXleft$ and $nodeXright$ can be zero if the alignment is a containment or insertion, but the type of the intervals does not recognize it as a containment or insertion.
		public double avgDiffs;  // -1=unknown.
		public int nAlignmentsForward, nAlignmentsBackward;
		public boolean implied;  // From $Factorize.writeIntervalConnectionFile$
		public int orientation;  // 0=same orientation; 1=opposite; 2=both; -1=unknown.
		public int flexibleOverlap;  // Bit $i$ is one iff the $i$-th overlap type is valid and flexible.
		
		/**
		 * Temporary space
		 */
		public int printed;
		public long[] printedComponents;  // Refers to the $components$ field of $nodeID1$
		public boolean usedInKernelGraph;
		
		
		public Edge(AlignmentPair pair) {
			set(pair);
		}
		
		
		/**
		 * Removes all information about edge type
		 */
		public final void clearTypes() {
			setType_noContainment(); containmentSubtype=-1;
			setType_noOverlap();
			setType_noInsertion();
			setType(Constants.SHARED_SUBSTRING,false);
		}
		
		
		public final void set(AlignmentPair pair) {
			int i, t;
			
			nodeID1=pair.node1.nodeID;
			nodeID2=pair.node2.nodeID;			
			orientation=pair.orientation?0:1;
			clearTypes();
			printed=-1; nAlignmentsForward=0; nAlignmentsBackward=0; avgDiffs=0.0;
			implied=false; on=true; supplement=false;
			containmentOverhangLeft=-1; containmentOverhangRight=-1;
			
			// Assigning Step2 type
			i=-1;
			i=IntervalGraphStep2.inStep2(pair);
			if (i==-1) { 
				System.err.println("ERROR: The following edge is not in Step 2: "+pair.toString());
				System.exit(1);
			}
			if (i<=Constants.OVERLAP_SUFFIX_SUFFIX) setType(i,true);
			else if (i<=Constants.CONTAINMENT_TWO_IN_ONE) setType(i,true);
			else if (i<=Constants.INSERTION_BOTH) setType(i,true);
			else if (i==Constants.SHARED_SUBSTRING) setType(i,true);
			
			// Assigning overhangs
			t=-1;
			if (i<=Constants.OVERLAP_SUFFIX_SUFFIX) t=i;
			if (t!=-1) {
				overhangs = new int[4][2];
				setOverhangs(t,pair);
			}
			t=-1;
			if (i<=Constants.CONTAINMENT_TWO_IN_ONE) t=i;
			if (t!=-1) {
				if (t==Constants.CONTAINMENT_IDENTICAL) {
					containmentSubtype=-1;
					containmentOverhangLeft=-1;
					containmentOverhangRight=-1;
				}
				else {
					containmentSubtype=IntervalGraphStep2.inStep2_getContainmentSubtype(t,pair,tmpArray);
					containmentOverhangLeft=tmpArray[0];
					containmentOverhangRight=tmpArray[1];
				}
			}
			t=-1;
			if (i<=Constants.INSERTION_BOTH) t=i;
			if (t!=-1 && t!=Constants.INSERTION_BOTH) {
				containmentSubtype=IntervalGraphStep2.inStep2_getContainmentSubtype(t,pair,tmpArray);
				containmentOverhangLeft=tmpArray[0];
				containmentOverhangRight=tmpArray[1];
			}
			t=-1;
			if (i==Constants.SHARED_SUBSTRING) t=i;
			if (t!=-1) {
				sharedSubstringOverhangs = new int[4];
				sharedSubstringOverhangs[0]=pair.alignmentStart1>pair.node1.start?pair.alignmentStart1-pair.node1.start:0;
				sharedSubstringOverhangs[1]=pair.node1.end>pair.alignmentEnd1?pair.node1.end-pair.alignmentEnd1:0;
				sharedSubstringOverhangs[2]=pair.alignmentStart2>pair.node2.start?pair.alignmentStart2-pair.node2.start:0;
				sharedSubstringOverhangs[3]=pair.node2.end>pair.alignmentEnd2?pair.node2.end-pair.alignmentEnd2:0;
			}
			if (pair.node1.type==Constants.INTERVAL_PERIODIC && pair.node2.type!=Constants.INTERVAL_PERIODIC) pair.node2.addPeriodicSubstring(Math.max(pair.node2.start,pair.alignmentStart2)-pair.node2.start,Math.min(pair.node2.end,pair.alignmentEnd2)-pair.node2.start,pair.node1.hasLongPeriod);
			else if (pair.node1.type!=Constants.INTERVAL_PERIODIC && pair.node2.type==Constants.INTERVAL_PERIODIC) pair.node1.addPeriodicSubstring(Math.max(pair.node1.start,pair.alignmentStart1)-pair.node1.start,Math.min(pair.node1.end,pair.alignmentEnd1)-pair.node1.start,pair.node2.hasLongPeriod);
			
			// Number of alignments
			if (pair.orientation) nAlignmentsForward=1;
			else nAlignmentsBackward=1;
			avgDiffs=pair.diffs;
			implied=pair.implied;
		}
		
		
		/**
		 * Ad hoc constructor used just to add shared substring edges induced by overlaps 
		 * in the same read (see procedures $addSameReadEdges*()$).
		 *
		 * @param type type of edge determined by the shared substring procedure.
		 */
		public Edge(int n1, int n2, int type) {
			nodeID1=n1;
			nodeID2=n2;
			
			orientation=0;
			clearTypes();
			printed=-1; nAlignmentsForward=0; nAlignmentsBackward=0; avgDiffs=-1;
			on=true; supplement=false;
			
			if (type==Constants.CONTAINMENT_IDENTICAL) {
				setType(type,true);
				containmentSubtype=-1;
				containmentOverhangLeft=-1;
				containmentOverhangRight=-1;
			}
			else if (type==Constants.SHARED_SUBSTRING) {
				setType(type,true);
				sharedSubstringOverhangs=null;
			}
			else {
				System.err.println("ERROR: egde type not recognized: "+type);
				System.exit(1);
			}
			
			nAlignmentsForward=1;  // Artificial alignment, to avoid divisions by zero.
			nAlignmentsBackward=0;
		}
		
		
		/**
		 * NOP constructor
		 */
		public Edge() {
			clear();
		}
		
		
		/**
		 * Removes all properties from the edge.
		 */
		public final void clear() {
			nodeID1=-1;
			nodeID2=-1;
			on=false;
			supplement=false;
			
			clearTypes();
			isRedundant=false;
			overhangs=null; maxOverhangs=null;
			containmentOverhangLeft=-1;
			containmentOverhangRight=-1;
			sharedSubstringOverhangs=null;
			avgDiffs=-1;
			nAlignmentsForward=0;
			nAlignmentsBackward=0;
			implied=false;
			orientation=-1;
			flexibleOverlap=0;
		}
		
		
		public final int getNAlignments() {
			return nAlignmentsForward+nAlignmentsBackward;
		}
		
		
		public final int getTo(int from) {
			return from==nodeID1?nodeID2:nodeID1;
		}
		
		
		public final boolean isSelfLoop() {
			return nodeID1==nodeID2;
		}
		
		
		public final boolean sameRead() {
			return nodesArray[nodeID1].read==nodesArray[nodeID2].read;
		}
		
		
		/**
		 * @param type 0=containment; 1=overlap; 2=insertion; 3=shared substring;
		 * @return a number smaller than $Math.POSITIVE_INFINITY$ if the number of 
		 * alignments can be determined; $Math.POSITIVE_INFINITY$ otherwise.
		 */
		public final int nAlignmentsOfType(int type) {
			final int nAlignments = getNAlignments();
			int n;
			
			if (type==0) {
				if (!isType_containment()) return 0;
				if (nAlignments==1 || nAlignments==nTypes()) return 1;
				return Math.POSITIVE_INFINITY;
			}
			else if (type==1) {
				if (!isType_overlap()) return 0;
				if (nAlignments==1) return 1;
				n=nOverlapTypes();
				return nAlignments==nTypes()?n:Math.POSITIVE_INFINITY;
				// $maxOverhangs=null$ does not guarantee that there is just one overlap 
				// alignment.
			}
			else if (type==2) {
				if (!isType_insertion()) return 0;
				if (nAlignments==1) return 1;
				if (insertion==Constants.INSERTION_BOTH) return nAlignments==nTypes()?2:Math.POSITIVE_INFINITY;
				return nAlignments==nTypes()?1:Math.POSITIVE_INFINITY;
			}
			else if (type==3) {
				if (!isType_sharedSubstring()) return 0;
				if (nAlignments==1 || nAlignments==nTypes()) return 1;
				return Math.POSITIVE_INFINITY;
			}
			return 0;
		}
		
		
		/**
		 * @param side TRUE=prefix, FALSE=suffix;
		 * @return the number of distinct overlap alignment types that use side $side$ of 
		 * $node$. This is a lower bound on the number of alignments that use side $side$ 
		 * of $node$.
		 */
		public final int nOverlapAlignments(int node, boolean side) {
			if (!isType_overlap()) return 0;
			int out = 0;
			if (side) {
				if ((overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) out++;
				if ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && node==nodeID1) out++;
				if ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && node==nodeID2) out++;
			}
			else {
				if ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && node==nodeID2) out++;
				if ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && node==nodeID1) out++;
				if ((overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) out++;
			}
			return out;
		}
		
		
		private final void setOverhangs(int type, AlignmentPair pair) {
			int i, j;
			
			if (overhangs==null) overhangs = new int[4][2];
			Math.set(overhangs,-1); maxOverhangs=null;
			if (type==Constants.OVERLAP_PREFIX_PREFIX) {
				overhangs[0][0]=pair.node1.end+1>pair.alignmentEnd1?pair.node1.end+1-pair.alignmentEnd1:0;
				overhangs[0][1]=pair.node2.end+1>pair.alignmentEnd2?pair.node2.end+1-pair.alignmentEnd2:0;
			}
			else if (type==Constants.OVERLAP_PREFIX_SUFFIX) {
				overhangs[1][0]=pair.node1.end+1>pair.alignmentEnd1?pair.node1.end+1-pair.alignmentEnd1:0;
				overhangs[1][1]=pair.alignmentStart2>pair.node2.start?pair.alignmentStart2-pair.node2.start:0;
			}
			else if (type==Constants.OVERLAP_SUFFIX_PREFIX) {
				overhangs[2][0]=pair.alignmentStart1>pair.node1.start?pair.alignmentStart1-pair.node1.start:0;
				overhangs[2][1]=pair.node2.end+1>pair.alignmentEnd2?pair.node2.end+1-pair.alignmentEnd2:0;
			}
			else if (type==Constants.OVERLAP_SUFFIX_SUFFIX) {
				overhangs[3][0]=pair.alignmentStart1>pair.node1.start?pair.alignmentStart1-pair.node1.start:0;
				overhangs[3][1]=pair.alignmentStart2>pair.node2.start?pair.alignmentStart2-pair.node2.start:0;
			}
		}
		
		
		/**
		 * Stores in $overhangs$ (respectively, in $maxOverhangs$) the min (respectively, 
		 * max) between the non-null overhangs of this edge and those of $other$, and 
		 * updates $flexibleOverlap$.
		 * 
		 * Remark: $maxOverhangs$ is allocated only if strictly necessary.
		 */
		private final void copyOverhangs(Edge other) {
			if (other.overhangs==null) return;
			boolean allocateMax;
			int i, b;
			int[][] otherMax;
			
			// Allocating $overhangs$ and $maxOverhangs$, if needed.
			if (overhangs==null) {
				overhangs = new int[4][2];
				Math.set(overhangs,-1);
			}
			allocateMax=false;
			if (maxOverhangs==null) {
				if (other.maxOverhangs!=null) allocateMax=true;
				else {
					for (i=0; i<4; i++) {
						if ( (overhangs[i][0]!=-1 && other.overhangs[i][0]!=-1 && overhangs[i][0]!=other.overhangs[i][0]) ||
							 (overhangs[i][1]!=-1 && other.overhangs[i][1]!=-1 && overhangs[i][1]!=other.overhangs[i][1])
						   ) {
							allocateMax=true;
							break;
						}
					}
				}
			}
			if (allocateMax) {
				maxOverhangs = new int[4][2];
				Math.copy(overhangs,maxOverhangs);
			}
			
			// Updating $overhangs$
			for (i=0; i<4; i++) {
				if (overhangs[i][0]==-1) {
					if (other.overhangs[i][0]!=-1) overhangs[i][0]=other.overhangs[i][0];
				}
				else {
					if (other.overhangs[i][0]!=-1) overhangs[i][0]=Math.min(overhangs[i][0],other.overhangs[i][0]);
				}
				if (overhangs[i][1]==-1) {
					if (other.overhangs[i][1]!=-1) overhangs[i][1]=other.overhangs[i][1];
				}
				else {
					if (other.overhangs[i][1]!=-1) overhangs[i][1]=Math.min(overhangs[i][1],other.overhangs[i][1]);
				}
			}
			
			// Updating $maxOverhangs$
			if (maxOverhangs==null) return;
			otherMax=other.maxOverhangs!=null?other.maxOverhangs:other.overhangs;
			for (i=0; i<4; i++) {
				if (maxOverhangs[i][0]==-1) {
					if (otherMax[i][0]!=-1) maxOverhangs[i][0]=otherMax[i][0];
				}
				else {
					if (otherMax[i][0]!=-1) maxOverhangs[i][0]=Math.max(maxOverhangs[i][0],otherMax[i][0]);
				}
				if (maxOverhangs[i][1]==-1) {
					if (otherMax[i][1]!=-1) maxOverhangs[i][1]=otherMax[i][1];
				}
				else {
					if (otherMax[i][1]!=-1) maxOverhangs[i][1]=Math.max(maxOverhangs[i][1],otherMax[i][1]);
				}
			}
			
			// Updating $flexibleOverlap$
			b=Constants.OVERLAP_PREFIX_PREFIX;
			for (i=0; i<4; i++) {
				if ((overhangs[i][0]!=-1 && maxOverhangs[i][0]!=overhangs[i][0]) || (overhangs[i][1]!=-1 && maxOverhangs[i][1]!=overhangs[i][1])) flexibleOverlap|=b;
				b<<=1;
			}
		}
		
		
		/**
		 * Stores in $containmentOverhang{Left,Right}$ the overhangs that correspond to
		 * the containment or insertion, between the current edge and $other$, that 
		 * maximizes the length of the substring of the containing interval.
		 *
		 * Remark: the procedure assumes the current edge and $other$ to be of the same
		 * containment or insertion type, which is different from $CONTAINMENT_IDENTICAL$
		 * and $INSERTION_BOTH$, and to be incident to the same nodes.
		 */
		private final void copyContainmentOverhangs(Edge other) {
			if (-other.containmentOverhangLeft-other.containmentOverhangRight>-containmentOverhangLeft-containmentOverhangRight) {
				containmentOverhangLeft=other.containmentOverhangLeft;
				containmentOverhangRight=other.containmentOverhangRight;
			}
		}
		
		
		/**
		 * Stores in $sharedSubstringOverhangs$ the overhangs that correspond to the 
		 * shared substring, between the current edge and $other$, that maximizes the sum 
		 * of the lengths in both reads.
		 *
		 * Remark: the procedure assumes the current edge and $other$ to be incident to 
		 * the same nodes, and to have non-null $sharedSubstringOverhangs$.
		 */
		private final void copySharedSubstringOverhangs(Edge other) {
			int support, otherSupport;
			
			if (sharedSubstringOverhangs[0]<=0 && sharedSubstringOverhangs[1]<=0 && sharedSubstringOverhangs[2]<=0 && sharedSubstringOverhangs[3]<=0) support=0;
			else support=-sharedSubstringOverhangs[0]-sharedSubstringOverhangs[1]-sharedSubstringOverhangs[2]-sharedSubstringOverhangs[3];			
			if (other.sharedSubstringOverhangs[0]<=0 && other.sharedSubstringOverhangs[1]<=0 && other.sharedSubstringOverhangs[2]<=0 && other.sharedSubstringOverhangs[3]<=0) otherSupport=0;
			else otherSupport=-other.sharedSubstringOverhangs[0]-other.sharedSubstringOverhangs[1]-other.sharedSubstringOverhangs[2]-other.sharedSubstringOverhangs[3];
			if (otherSupport>support) {
				if (nodeID1==other.nodeID1) System.arraycopy(other.sharedSubstringOverhangs,0,sharedSubstringOverhangs,0,4);
				else {
					sharedSubstringOverhangs[0]=other.sharedSubstringOverhangs[2];
					sharedSubstringOverhangs[1]=other.sharedSubstringOverhangs[3];
					sharedSubstringOverhangs[2]=other.sharedSubstringOverhangs[0];
					sharedSubstringOverhangs[3]=other.sharedSubstringOverhangs[1];
				}
			}
		}
		
		
		/**
		 * Remark: containment (not necessarily proper) and overlap are not mutually 
		 * exclusive. If both happen, then the containing interval has an internal repeat 
		 * that occurs at one of its ends. If containment equals identity, then the 
		 * containing interval has a border. The same holds between insertion and overlap.
		 *
		 * Containment and insertion are also both possible at the same time, if 
		 * containment is not identity. In this case, the $containmentOverhang*$ variables
		 * are assumed to refer to containment.
		 *
		 * Remark: in the interval graph pipeline we don't test for strict overlap, strict
		 * insertion, etc. This can be interpreted as allowing all possible alternatives.
		 */
		public final void addTypes(Edge other) {
			int tmp;
			
			// Containment
			if (other.isType_containment() && !isType(Constants.CONTAINMENT_IDENTICAL)) {
				if (!isType_containment()) {
					copyType_containment(other);
					if (!other.isType(Constants.CONTAINMENT_IDENTICAL)) {
						containmentOverhangLeft=other.containmentOverhangLeft;
						containmentOverhangRight=other.containmentOverhangRight;
					}
				}
				else {
					if (other.isType(Constants.CONTAINMENT_IDENTICAL)) {
						setType_noContainment();
						setType(Constants.CONTAINMENT_IDENTICAL,true);
					}
					else if ( (isType(Constants.CONTAINMENT_ONE_IN_TWO) && other.isType(Constants.CONTAINMENT_TWO_IN_ONE)) ||
						      (isType(Constants.CONTAINMENT_TWO_IN_ONE) && other.isType(Constants.CONTAINMENT_ONE_IN_TWO))
					        ) {
						setType_noContainment();
						setType(Constants.CONTAINMENT_IDENTICAL,true);
					}
					else {
						// Same type
						copyContainmentOverhangs(other);
					}
				}
			}
			
			// Allowing multiple types of overlap between the same pair of nodes
			if (!isType_overlap()) {
				if (other.isType_overlap()) {
					copyType_overlap(other);
					copyOverhangs(other);
				}
			}
			else {
				if (other.isType_overlap()) {
					addType_overlap(other);
					copyOverhangs(other);
				}
			}
			
			// Insertion
			if (other.isType_insertion() && !isType(Constants.INSERTION_BOTH)) {
				if (!isType_insertion()) {
					copyType_insertion(other);
					if ((!isType_containment() || isType(Constants.CONTAINMENT_IDENTICAL)) && (!other.isType_containment() || other.isType(Constants.CONTAINMENT_IDENTICAL))) {
						containmentOverhangLeft=other.containmentOverhangLeft;
						containmentOverhangRight=other.containmentOverhangRight;
					}
				}
				else if (differentType_insertion(other)) {
					// The two intervals are approximately identical, but they have not
					// been marked as $CONTAINMENT_IDENTICAL$ by previous procedures.
					setType_noInsertion();
					setType(Constants.INSERTION_BOTH,true);
				}
				else {
					// Same type
					copyContainmentOverhangs(other);
				}
			}
			
			// Shared substring
			if (other.isType(Constants.SHARED_SUBSTRING)) {
				tmp=sharedSubstring;
				copyType_sharedSubstring(other);
				if (other.sharedSubstringOverhangs!=null) {
					if (tmp==-1 || sharedSubstringOverhangs==null) {
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						System.arraycopy(other.sharedSubstringOverhangs,0,sharedSubstringOverhangs,0,4);
					}
					else copySharedSubstringOverhangs(other);
				}
			}
			
			// Disallowing identity and insertion
			if (isType(Constants.CONTAINMENT_IDENTICAL) && isType_insertion()) setType_noInsertion();
		}

		
		/**
		 * Remark: taking the minimum rather than averaging diff values does not allow one
		 * to cluster a periodic component by just removing edges with high value.
		 *
		 * Remark: the procedure handles the case in which $avgDiffs=-1$.
		 */
		public final void addDiffs(Edge other) {
			nAlignmentsForward+=other.nAlignmentsForward;
			nAlignmentsBackward+=other.nAlignmentsBackward;
			if (avgDiffs!=-1) {
				if (other.avgDiffs!=-1) avgDiffs+=other.avgDiffs;
				else avgDiffs=-1;
			}
		}
		
		
		public final void addOrientation(Edge other) {
			if (orientation==-1 || orientation==2) return;
			if (other.orientation==-1) orientation=-1;
			else if (other.orientation!=orientation) orientation=2;
		}
		
		
		/**
		 * Assume that the edge contains exactly two alignments. The procedure arbitrarily
		 * keeps just the alignment that maximizes the common area between the intervals.
		 *
		 * Remark: the procedure assumes that $avgDiffs$ is a sum, not an avg.
		 */
		public final void enforceOneAlignment() {
			int i;
			int row, containmentSurface, insertionSurface, overlapSurface, sharedSubstringSurface;
			int overlapSurface1, overlapSurface2, overlapType1, overlapType2;
			
			avgDiffs/=2;
			flexibleOverlap=0;
			
			// The two alignments have the same type
			if (isType(Constants.INSERTION_BOTH)) {
				setType(Constants.CONTAINMENT_IDENTICAL,true);
				setType_noInsertion();
				setType_noOverlap();
				setType_noSharedSubstring();
				if (orientation==0) nAlignmentsForward=1;
				else if (orientation==1) nAlignmentsBackward=1;
				else {
					orientation=-1;
					nAlignmentsForward=1;  // Arbitrary choice between forward and backward
					nAlignmentsBackward=0;
				}
				return;
			}
			if (isType_overlap() && !isType_containment() && !isType_insertion() && !isType_sharedSubstring() && overhangs!=null) {
				overlapSurface1=-1; overlapSurface2=-1; overlapType1=-1; overlapType2=-1;
				for (i=0; i<4; i++) {
					if ((overlap&(1<<i))==0) continue;
					if (overlapSurface1==-1) {
						overlapSurface1=nodesArray[nodeID1].length()-overhangs[i][0]+nodesArray[nodeID2].length()-overhangs[i][1];
						overlapType1=1<<i;
					}
					else {
						overlapSurface2=nodesArray[nodeID1].length()-overhangs[i][0]+nodesArray[nodeID2].length()-overhangs[i][1];
						overlapType2=1<<i;
						break;
					}
				}
				if (overlapSurface1>=overlapSurface2) overlap=overlapType1;
				else overlap=overlapType2;
				orientation=Constants.overlap2orientation(overlap);
				if (orientation==0) nAlignmentsForward=1;
				else nAlignmentsBackward=1;
				return;
			}
			if (isType_sharedSubstring() && !isType_containment() && !isType_insertion() && !isType_overlap()) {
				// We do not change $sharedSubstringOverhangs$.
				if (orientation==0) nAlignmentsForward=1;
				else if (orientation==1) nAlignmentsBackward=1;
				else {
					orientation=-1;
					nAlignmentsForward=1;  // Arbitrary choice between forward and backward
					nAlignmentsBackward=0;
				}
				return;
			}
			
			// The two alignments have different types
			if (orientation==0) nAlignmentsForward=1;
			else if (orientation==1) nAlignmentsBackward=1;
			else {
				orientation=-1;
				nAlignmentsForward=1;  // Arbitrary choice between forward and backward
				nAlignmentsBackward=0;
			}
			if (isType(Constants.CONTAINMENT_IDENTICAL)) {
				setType_noInsertion();
				setType_noOverlap();
				setType_noSharedSubstring();
				return;
			}
			containmentSurface=0;
			if (isType(Constants.CONTAINMENT_ONE_IN_TWO) || isType(Constants.CONTAINMENT_TWO_IN_ONE)) containmentSurface=nodesArray[nodeID1].length()+nodesArray[nodeID2].length()-containmentOverhangLeft-containmentOverhangRight;
			insertionSurface=0;
			if (isType(Constants.INSERTION_ONE_IN_TWO) || isType(Constants.INSERTION_TWO_IN_ONE)) insertionSurface=nodesArray[nodeID1].length()+nodesArray[nodeID2].length()-containmentOverhangLeft-containmentOverhangRight;
			overlapSurface=0;
			if (isType_overlap() && overhangs!=null) {
				row=Constants.overlap2id(overlap);
				overlapSurface=nodesArray[nodeID1].length()-overhangs[row][0]+nodesArray[nodeID2].length()-overhangs[row][1];
			}
			sharedSubstringSurface=0;
			if (isType(Constants.SHARED_SUBSTRING) && sharedSubstringOverhangs!=null) sharedSubstringSurface=nodesArray[nodeID1].length()-sharedSubstringOverhangs[0]-sharedSubstringOverhangs[1]+nodesArray[nodeID2].length()-sharedSubstringOverhangs[2]-sharedSubstringOverhangs[3];
			if (containmentSurface>=insertionSurface && containmentSurface>=overlapSurface && containmentSurface>=sharedSubstringSurface) {
				setType_noInsertion();
				setType_noOverlap();
				setType_noSharedSubstring();
			}
			else if (insertionSurface>=containmentSurface && insertionSurface>=overlapSurface && insertionSurface>=sharedSubstringSurface) {
				setType_noContainment();
				setType_noOverlap();
				setType_noSharedSubstring();
			}
			else if (overlapSurface>=containmentSurface && overlapSurface>=insertionSurface && overlapSurface>=sharedSubstringSurface) {
				setType_noContainment();
				setType_noInsertion();
				setType_noSharedSubstring();
			}
			else if (sharedSubstringSurface>=containmentSurface && sharedSubstringSurface>=insertionSurface && sharedSubstringSurface>=overlapSurface) {
				setType_noContainment();
				setType_noInsertion();
				setType_noOverlap();
			}
		}
		
		
		public boolean equals(Object other) {
			Edge otherEdge = (Edge)other;
			return (otherEdge.nodeID1==nodeID1 && otherEdge.nodeID2==nodeID2) || 
				   (otherEdge.nodeID1==nodeID2 && otherEdge.nodeID2==nodeID1);
		}
		
		
		public int hashCode() {
			return (nodeID1+"-"+nodeID2).hashCode();
		}
		
		
		/**
		 * Puts all edges with $on=true$ first; if $order=UNSORTED$, returns. If 
		 * $order > UNSORTED$, further sorts by $getTo(order-(UNSORTED+1))$. If 
		 * $order < UNSORTED$, further puts all overlap edges first, then puts all edges 
		 * that overlap with a suffix of node $(UNSORTED-1)-order$ first.
		 */
		public int compareTo(Object other) {
			Edge otherEdge = (Edge)other;
			if (on) {
				if (!otherEdge.on) return -1;
			}
			else {
				if (otherEdge.on) return 1;
			}
			if (order==UNSORTED) return 0;
			if (order>UNSORTED) {
				int to = getTo(order-(UNSORTED+1));
				int otherTo = otherEdge.getTo(order-(UNSORTED+1));
				if (to<otherTo) return -1;
				if (to>otherTo) return 1;
			}
			else {
				if (overlap==-1) {
					if (otherEdge.overlap!=-1) return 1;
				}
				else {
					if (otherEdge.overlap==-1) return -1;
					int node = (UNSORTED-1)-order;
					if (isSuffixOverlap(node)) {
						if (!otherEdge.isSuffixOverlap(node)) return -1;
					}
					else {
						if (otherEdge.isSuffixOverlap(node)) return 1;
					}
				}
			}
			return 0;
		}
		
		
		public final String toDot(double threshold) {
			final double WIDTH_SCALE = 0.25;
			boolean periodicNonperiodic = false;
			if ( (nodesArray[nodeID1].type==Constants.INTERVAL_PERIODIC && nodesArray[nodeID2].type!=Constants.INTERVAL_PERIODIC) ||
				 (nodesArray[nodeID1].type!=Constants.INTERVAL_PERIODIC && nodesArray[nodeID2].type==Constants.INTERVAL_PERIODIC)
			   ) periodicNonperiodic=true;
			String color = "gray";
			if (containment!=-1) color="green";
			else if (overlap!=-1) color="blue";
			else if (insertion!=-1) color="red";
			else if (sharedSubstring!=-1) color="yellow";
			//return nodeID1+" -- "+nodeID2+" [color=\""+(periodicNonperiodic?"red":"black")+"\",penwidth=\""+(getNAlignments()*WIDTH_SCALE)+"\"];\n";
			return nodeID1+" -- "+nodeID2+" [color=\""+color+"\"];\n";
			
			
			
/*			if (nAlignmentsForward>8 || nAlignmentsBackward>8) return nodeID1+" -- "+nodeID2+" [color=red]\n;";
			else return nodeID1+" -- "+nodeID2+" [color=black];\n";
*/			
//			return isPrimary?nodeID1+" -- "+nodeID2+";\n":"";
			
//			return nodeID1+" -- "+nodeID2+" [color=\""+(isPrimary?"black":"red")+"\"];\n";
			
/*			if (avgDiffs<threshold) return nodeID1+" -- "+nodeID2+";\n";
			return "";
*/
/*
			String color = "gray";
			if (containment!=-1) color="green";
			else if (overlap!=-1) color="blue";
			else if (insertion!=-1) color="red";
			else if (sharedSubstring!=-1) color="yellow";
			return nodeID1+" -- "+nodeID2+" [color="+color+"];\n";
*/			
			
/*			
			String tagsOther = " [color=black];\n";
			String tagsContainment = " [color=gray];\n";
			String tagsInsertion = " [color=red];\n";
			String tagsHeadTail = "";
			if (overlap==Constants.OVERLAP_PREFIX_PREFIX) tagsHeadTail="arrowhead=dot,arrowtail=dot";
			else if (overlap==Constants.OVERLAP_PREFIX_SUFFIX) tagsHeadTail="arrowhead=dot,arrowtail=odot";
			else if (overlap==Constants.OVERLAP_SUFFIX_PREFIX) tagsHeadTail="arrowhead=odot,arrowtail=dot";
			else if (overlap==Constants.OVERLAP_SUFFIX_SUFFIX) tagsHeadTail="arrowhead=odot,arrowtail=odot";
			String tagsOverlap = " [color=green,penwidth=3,"+tagsHeadTail+"];\n";
			
			if (containment==CONTAINMENT_ONE_IN_TWO) return nodeID1+" -> "+nodeID2+tagsContainment;
			else if (containment==CONTAINMENT_TWO_IN_ONE) return nodeID2+" -> "+nodeID1+tagsContainment;
			else if (containment==CONTAINMENT_IDENTICAL) return nodeID1+" -- "+nodeID2+tagsContainment;
			else if (overlap!=-1) return nodeID1+" -- "+nodeID2+tagsOverlap;
			else if (insertion==INSERTION_ONE_IN_TWO) return nodeID1+" -> "+nodeID2+tagsInsertion;
			else if (insertion==INSERTION_TWO_IN_ONE) return nodeID2+" -> "+nodeID1+tagsInsertion;
			else if (insertion==INSERTION_BOTH) return nodeID1+" -- "+nodeID2+tagsInsertion;
			return nodeID1+" -- "+nodeID2+tagsOther;
*/			
			/*
			String tags = " [containment="+containment+",overlap="+overlap+",insertion="+insertion+"];\n";
			if (containment==CONTAINMENT_ONE_IN_TWO || insertion==INSERTION_ONE_IN_TWO) return nodeID1+" -> "+nodeID2+tags;
			else if (containment==CONTAINMENT_TWO_IN_ONE || insertion==INSERTION_TWO_IN_ONE) return nodeID2+" -> "+nodeID1+tags;
			return nodeID1+" -> "+nodeID2+tags+nodeID2+" -> "+nodeID1+tags;			
			*/
		}
		
		
		public String toString() {
			String out = "     nodeID1="+nodeID1+", nodeID2="+nodeID2+" orient="+orientation+" nAlignments="+getNAlignments()+" (C="+containment+",C'="+containmentSubtype+",O="+overlap+",I="+insertion+",S="+sharedSubstring+") supplement="+supplement+", on="+on+" \n";
			out+="node1: "+(nodesArray[nodeID1]==null?"NULL("+nodeID1+")\n":" length="+nodesArray[nodeID1].length()+" period="+nodesArray[nodeID1].period+" read="+nodesArray[nodeID1].read+"\n");
			out+="node2: "+(nodesArray[nodeID2]==null?"NULL("+nodeID2+")\n":" length="+nodesArray[nodeID2].length()+" period="+nodesArray[nodeID2].period+" read="+nodesArray[nodeID2].read+"\n");
			out+="overlap overhangs: \n";
			if (overhangs!=null) {
				for (int i=0; i<overhangs.length; i++) {
					for (int j=0; j<overhangs[i].length; j++) out+=overhangs[i][j]+",";
					out+="\n";
				}
			}
			out+="containmentOverhangs: "+containmentOverhangLeft+","+containmentOverhangRight+"\n";
			out+="shared substring overhangs: \n";
			if (sharedSubstringOverhangs!=null) {
				for (int i=0; i<sharedSubstringOverhangs.length; i++) out+=sharedSubstringOverhangs[i]+",";
			}
			out+="\n";
			return out;	
			
//			return nodeID1+"--"+nodeID2+" containment="+containment+", overlap="+overlap+", insertion="+insertion+", sharedSubstring="+sharedSubstring+", nAlignments="+getNAlignments()+", avgDiffs="+avgDiffs+", orientation="+orientation+", supplement="+supplement;
		}
		
		
		/**
		 * @return 0 (respectively, 1) if the edge is of overlap type and $node$ uses its 
		 * prefix (respectively, suffix). -1 if the edge is not of overlap type. 2 if 
		 * $node$ uses both ends.
		 */
		public final int getOverlapEnd(int node) {
			boolean prefixEnd = false;
			boolean suffixEnd = false;
			
			if (overlap==-1) return -1;
			if ((overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) prefixEnd=true;
			if ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) {
				if (node==nodeID1) prefixEnd=true;
				else suffixEnd=true;
			}
			if ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) {
				if (node==nodeID1) suffixEnd=true;
				else prefixEnd=true;
			}
			if ((overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) suffixEnd=true;
			if (prefixEnd) {
				if (!suffixEnd) return 0;
				else return 2;
			}
			else if (suffixEnd) return 1;
			return -1;
		}
		
		
		/**
		 * The procedure adds pairs $(p,s)$ to $out[k..X]$, where $p$ is the position of a 
		 * maximal event induced by overlap, with respect to the beginning of $node$, and 
		 * $s$ is zero (respectively, one) if the overlap lies to the left (respectively, 
		 * right) of the event in $node$. The procedure returns $X$, or $k-1$ if no 
		 * maximal event is represented by this egde.
		 *
		 * Remark: an edge can encode multiple maximal overlaps.
		 */
		public final int isMaximalOverlap(int node, int[] out, int k) {
			int last = k-1;
			
			if (overlap==-1) return last;
			if ((overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) {
				if (node==nodeID1 && nodesArray[nodeID2].isLeftMaximal) {
					out[++last]=nodesArray[nodeID1].length()-overhangs[0][0]-1;
					if (out[last]<0) out[last]=0;
					out[++last]=0;
				}
				if (node==nodeID2 && nodesArray[nodeID1].isLeftMaximal) {
					out[++last]=nodesArray[nodeID2].length()-overhangs[0][1]-1;
					if (out[last]<0) out[last]=0;
					out[++last]=0;
				}
			}
			if ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) {
				if (node==nodeID1 && nodesArray[nodeID2].isRightMaximal) {
					out[++last]=nodesArray[nodeID1].length()-overhangs[1][0]-1;
					if (out[last]<0) out[last]=0;
					out[++last]=0;
				}
				if (node==nodeID2 && nodesArray[nodeID1].isLeftMaximal) {
					out[++last]=overhangs[1][1]-1;
					if (out[last]<0) out[last]=0;
					out[++last]=1;
				}
			}
			if ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) {
				if (node==nodeID1 && nodesArray[nodeID2].isLeftMaximal) {
 					out[++last]=overhangs[2][0]-1;
					if (out[last]<0) out[last]=0;
 					out[++last]=1;
				}
				if (node==nodeID2 && nodesArray[nodeID1].isRightMaximal) {
					out[++last]=nodesArray[nodeID2].length()-overhangs[2][1]-1;
					if (out[last]<0) out[last]=0;
 					out[++last]=0;
				}
			}
			if ((overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) {
				if (node==nodeID1 && nodesArray[nodeID2].isRightMaximal) {
					out[++last]=overhangs[3][0]-1;
					if (out[last]<0) out[last]=0;
 					out[++last]=1;
				}
				if (node==nodeID2 && nodesArray[nodeID1].isRightMaximal) {
					out[++last]=overhangs[3][1]-1;
					if (out[last]<0) out[last]=0;
 					out[++last]=1;
				}
			}
			return last;
		}
		
		
		/**
		 * The procedure adds pairs $(p,s)$ to $out[k..X]$, where $p$ is the position of a 
		 * maximal event induced by insertion or containment, with respect to the 
		 * beginning of $node$, and $s$ is zero (respectively, one) if the insertion lies 
		 * to the left (respectively, right) of the event in $node$. The procedure returns
		 * $X$, or $k-1$ if no maximal event is represented by this egde.
		 */
		public final int isMaximalInsertionOrContainment(int node, int[] out, int k) {
			int last = k-1; 
			int contained = -1;
			
			if ((insertion==-1 && containment==-1) || getNAlignments()!=1) return last;
			if ((insertion==Constants.INSERTION_ONE_IN_TWO || containment==Constants.CONTAINMENT_ONE_IN_TWO) && node==nodeID2) contained=nodeID1;
			else if ((insertion==Constants.INSERTION_TWO_IN_ONE || containment==Constants.CONTAINMENT_TWO_IN_ONE) && node==nodeID1) contained=nodeID2;
			else return last;
			if (containmentSubtype==Constants.CONTAINMENT_PREFIX) {
				if ( ((orientation==0 || orientation==2) && nodesArray[contained].isRightMaximal) ||
					 ((orientation==1 || orientation==2) && nodesArray[contained].isLeftMaximal)
				   ) {
					out[++last]=nodesArray[node].length()-containmentOverhangRight-1;
					if (out[last]<0) out[last]=0;
					out[++last]=0;
				}
			}
			else if (containmentSubtype==Constants.CONTAINMENT_SUFFIX) {
				if ( ((orientation==0 || orientation==2) && nodesArray[contained].isLeftMaximal) ||
					 ((orientation==1 || orientation==2) && nodesArray[contained].isRightMaximal)
				   ) {
					out[++last]=containmentOverhangLeft;
					out[++last]=1;
				}
			}
			else {
				if (nodesArray[contained].isLeftMaximal) {
					if (orientation==0 || orientation==2) {
						out[++last]=containmentOverhangLeft;
						out[++last]=1;
					}
					if (orientation==1 || orientation==2) {
						out[++last]=nodesArray[node].length()-containmentOverhangRight-1;
						if (out[last]<0) out[last]=0;
						out[++last]=0;
					}
				}
				if (nodesArray[contained].isRightMaximal) {
					if (orientation==0 || orientation==2) {
						out[++last]=nodesArray[node].length()-containmentOverhangRight-1;
						if (out[last]<0) out[last]=0;
						out[++last]=0;
					}
					if (orientation==1 || orientation==2) {
						out[++last]=containmentOverhangLeft;
						out[++last]=1;
					}
				}
			}
			return last;
		}
		
		
		/**
		 * The procedure adds pairs $(p,s)$ to $out[k..X]$, where $p$ is the position of a 
		 * maximal event induced by a shared substring, with respect to the beginning of 
		 * $node$, and $s$ is zero (respectively, one) if the shared substring lies to the
		 * left (respectively, right) of the event in $node$. The procedure returns $X$,
		 * or $k-1$ if no maximal event is represented by this egde, or if this edge
		 * represents more than one shared substring.
		 */
		public final int isMaximalSharedSubstring(int node, int[] out, int k) {
			final int IDENTITY_THRESHOLD = IO.quantum;
			int position, last;
			
			last=k-1;
			if (sharedSubstring==-1 || sharedSubstringOverhangs==null || getNAlignments()!=1) return last;
			if (nodesArray[nodeID1].read==nodesArray[nodeID2].read) return last;
			if (node==nodeID1) {
				position=sharedSubstringOverhangs[0];
				if (position>IDENTITY_THRESHOLD && position<nodesArray[node].length()-IDENTITY_THRESHOLD) {
					out[++last]=position;
					out[++last]=1;
				}
				position=nodesArray[node].length()-sharedSubstringOverhangs[1]-1;
				if (position<0) position=0;
				if (position>IDENTITY_THRESHOLD && position<nodesArray[node].length()-IDENTITY_THRESHOLD) {
					out[++last]=position;
					out[++last]=0;
				}
			}
			else if (node==nodeID2) {
				position=sharedSubstringOverhangs[2];
				if (position>IDENTITY_THRESHOLD && position<nodesArray[node].length()-IDENTITY_THRESHOLD) {
					out[++last]=position;
					out[++last]=1;
				}
				position=nodesArray[node].length()-sharedSubstringOverhangs[3]-1;
				if (position<0) position=0;
				if (position>IDENTITY_THRESHOLD && position<nodesArray[node].length()-IDENTITY_THRESHOLD) {
					out[++last]=position;
					out[++last]=0;
				}
			}
			return last;
		}
		
		
		/**
		 * Resets to $newNode$ the one of $nodeID1,nodeID2$ that equals $oldNode \neq 
		 * newNode$, where $newNode$ is substring $[newNodeStart..newNodeEnd]$ of 
		 * $oldNode$, and updates all properties of the edge accordingly.
		 *
		 * Remark: the procedure works only if the edge is of containment, insertion, 
		 * overlap, or shared substring type (i.e. if it has a proper-substring 
		 * projection onto $oldNode$).
		 *
		 * Remark: the procedure assumes that the edge contains just one alignment.
		 *
		 * @return true iff the transformation has been performed succesfully. If the 
		 * transformation has not been performed successfully, the state of the edge is
		 * not necessarily the one it was before the procedure was called.
		 */
		public final boolean transform(int oldNode, int newNode, int newNodeStart, int newNodeEnd, AlignmentPair tmpPair) {
			boolean transformed;
			int otherNode, tmp;	

			if (oldNode==nodeID1) otherNode=nodeID2;
			else otherNode=nodeID1;
			transformed=transform_impl(oldNode,newNode,newNodeStart,newNodeEnd,tmpPair);		
			if (transformed) {
				if (oldNode==nodeID1) nodeID1=newNode;
				else nodeID2=newNode;				
				if (nodesArray[nodeID1].type!=Constants.INTERVAL_PERIODIC && nodesArray[nodeID2].type!=Constants.INTERVAL_PERIODIC) return true;
				tmp=IntervalGraphStep2.inStep2_checkEdge(this);
				if (tmp!=-1) return true;
				else return false;
			}
			else return false;
		}

		
		private final boolean transform_impl(int oldNode, int newNode, int newNodeStart, int newNodeEnd, AlignmentPair tmpPair) {
			boolean transformed;
			int otherNode, projectionStart, projectionEnd, last;
			
			if (oldNode==nodeID1) otherNode=nodeID2;
			else otherNode=nodeID1;
			last=projectOnto(oldNode,tmpArray,true);
			if (last==-1) return false;
			projectionStart=tmpArray[0]; projectionEnd=tmpArray[1];
			transformed=false;
			if (Intervals.areApproximatelyIdentical(newNodeStart,newNodeEnd,projectionStart,projectionEnd)) {
				if ( (isType(Constants.CONTAINMENT_ONE_IN_TWO) && oldNode==nodeID2) ||
					 (isType(Constants.CONTAINMENT_TWO_IN_ONE) && oldNode==nodeID1) ||
					 (isType(Constants.INSERTION_ONE_IN_TWO) && oldNode==nodeID2) ||
					 (isType(Constants.INSERTION_TWO_IN_ONE) && oldNode==nodeID1)
				   ) {
					transform_identical_containmentOrInsertion();
					transformed=true;
				}
				else if (isType_overlap()) {
					transform_identical_overlap(oldNode,newNode);
					transformed=true;
				}
				else if (isType(Constants.SHARED_SUBSTRING)) {
					transform_identical_sharedSubstring(oldNode,newNode,newNodeStart,newNodeEnd);
					transformed=true;
				}
			}
			else if (Intervals.isApproximatelyContained(projectionStart,projectionEnd,newNodeStart,newNodeEnd)) {
				if ( (isType(Constants.CONTAINMENT_ONE_IN_TWO) && oldNode==nodeID2) ||
					 (isType(Constants.CONTAINMENT_TWO_IN_ONE) && oldNode==nodeID1) ||
					 (isType(Constants.INSERTION_ONE_IN_TWO) && oldNode==nodeID2) ||
					 (isType(Constants.INSERTION_TWO_IN_ONE) && oldNode==nodeID1)
				   ) {
					transform_contained_containmentOrInsertion(oldNode,newNode,newNodeStart,newNodeEnd,projectionStart,projectionEnd);
					transformed=true;
				}
				else if (isType_overlap()) {
					transform_contained_overlap(oldNode,newNode,newNodeStart,newNodeEnd,projectionStart,projectionEnd,tmpPair);
					transformed=true;
				}
				else if (isType(Constants.SHARED_SUBSTRING)) {
					transform_contained_sharedSubstring(oldNode,newNode,newNodeStart,newNodeEnd,projectionStart,projectionEnd,tmpPair);
					transformed=true;
				}
			}
			return transformed;
		}
		
		
		/**
		 * Same as $transform()$, but does not perform the actual transformation.
		 */
		public final boolean canBeTransformed(int oldNode, int newNode, int newNodeStart, int newNodeEnd, int newNodeType) {
			boolean transformed;
			int otherNode, projectionStart, projectionEnd, last;
			
			if (oldNode==nodeID1) otherNode=nodeID2;
			else otherNode=nodeID1;
			last=projectOnto(oldNode,tmpArray,true);
			if (last==-1) return false;
			projectionStart=tmpArray[0]; projectionEnd=tmpArray[1];
			transformed=false;
			if (Intervals.areApproximatelyIdentical(newNodeStart,newNodeEnd,projectionStart,projectionEnd)) {
				if ( (containment==Constants.CONTAINMENT_ONE_IN_TWO && oldNode==nodeID2) ||
					 (containment==Constants.CONTAINMENT_TWO_IN_ONE && oldNode==nodeID1) ||
					 (insertion==Constants.INSERTION_ONE_IN_TWO && oldNode==nodeID2) ||
					 (insertion==Constants.INSERTION_TWO_IN_ONE && oldNode==nodeID1)
				   ) transformed=true;
				else if (overlap!=-1) transformed=true;
				else if (sharedSubstring!=-1) transformed=true;
			}
			else if (Intervals.isApproximatelyContained(projectionStart,projectionEnd,newNodeStart,newNodeEnd)) {
				if ( (containment==Constants.CONTAINMENT_ONE_IN_TWO && oldNode==nodeID2) ||
					 (containment==Constants.CONTAINMENT_TWO_IN_ONE && oldNode==nodeID1) ||
					 (insertion==Constants.INSERTION_ONE_IN_TWO && oldNode==nodeID2) ||
					 (insertion==Constants.INSERTION_TWO_IN_ONE && oldNode==nodeID1)
				   ) transformed=true;
				else if (overlap!=-1) transformed=true;
				else if (sharedSubstring!=-1) transformed=true;
			}
			return transformed;
		}
		
		
		/**
		 * Writes to $out$ the first and last positions of the substrings of $node$ that 
		 * are covered by the alignments that correspond to this edge (relative to the 
		 * first position of $node$).
		 *
		 * Remark: the procedure considers only containment, insertion, overlap (iff
		 * $includeOverlap=TRUE$), and shared substring types.
		 * 
		 * @return the last position of $out$ that stores a projection, or -1 if no 
		 * projection was performed.
		 */
		public final int projectOnto(int node, int[] out, boolean includeOverlap) {
			int i, last;
			
			last=-1;
			if (containment!=-1) last=projectOnto_containment(node,out,last);
			if (insertion!=-1) last=projectOnto_insertion(node,out,last);
			if (includeOverlap && overlap!=-1) {
				for (i=Constants.OVERLAP_PREFIX_PREFIX; i<=Constants.OVERLAP_SUFFIX_SUFFIX; i<<=1) {
					if ((overlap&i)!=0) last=projectOnto_overlap(node,out,i,last);
				}
			}
			if (sharedSubstring!=-1) last=projectOnto_sharedSubstring(node,out,last);
			return last;
		}
		
		
		public final int projectOnto_containment(int node, int[] out, int last) {
			int pos = last;
			int start=Math.NEGATIVE_INFINITY, end=Math.NEGATIVE_INFINITY;
			
			if (containment==Constants.CONTAINMENT_ONE_IN_TWO) {
				if (node==nodeID2) {
					start=containmentOverhangLeft;
					end=nodesArray[nodeID2].length()-containmentOverhangRight-1;
				}
				else if (node==nodeID1) {
					start=0;
					end=nodesArray[nodeID1].length()-1;
				}
			}
			else if (containment==Constants.CONTAINMENT_TWO_IN_ONE) {
				if (node==nodeID1) {
					start=containmentOverhangLeft;
					end=nodesArray[nodeID1].length()-containmentOverhangRight-1;
				}
				else if (node==nodeID2) {
					start=0; 
					end=nodesArray[nodeID2].length()-1;
				}
			}
			else if (containment==Constants.CONTAINMENT_IDENTICAL) {
				if (node==nodeID1) {
					start=0; 
					end=nodesArray[nodeID1].length()-1;
				}
				else if (node==nodeID2) {
					start=0;
					end=nodesArray[nodeID2].length()-1;
				}
			}
			if (start!=Math.NEGATIVE_INFINITY && end!=Math.NEGATIVE_INFINITY) {
				if (start>=end) {
					System.err.println("projectOnto_containment> ERROR containment="+containment+" start="+start+" end="+end+" containmentOverhangLeft="+containmentOverhangLeft+" containmentOverhangRight="+containmentOverhangRight);			
					System.err.println("edge: "+this);
					System.err.println("node1 id="+nodeID1+": "+nodesArray[nodeID1]);
					System.err.println("node1.length="+nodesArray[nodeID1].length());
					System.err.println("node2 id="+nodeID2+": "+nodesArray[nodeID2]);
					System.err.println("node2.length="+nodesArray[nodeID2].length());
					if (overhangs!=null) {
						for (int x=0; x<4; x++) System.err.println(overhangs[x][0]+","+overhangs[x][1]);
					}
					System.exit(1);
				}
				out[++pos]=Math.max(start,0); 
				out[++pos]=Math.min(end,nodesArray[node].length()-1);
			}
			return pos;
		}
		
		
		public final int projectOnto_insertion(int node, int[] out, int last) {
			int pos = last;
			int start=Math.NEGATIVE_INFINITY, end=Math.NEGATIVE_INFINITY;
			
			if (insertion==Constants.INSERTION_ONE_IN_TWO) {
				if (node==nodeID2) {
					start=containmentOverhangLeft;
					end=nodesArray[nodeID2].length()-containmentOverhangRight-1;
				}
				else if (node==nodeID1) {
					start=0; 
					end=nodesArray[nodeID1].length()-1;
				}
			}
			else if (insertion==Constants.INSERTION_TWO_IN_ONE) {
				if (node==nodeID1) {
					start=containmentOverhangLeft;
					end=nodesArray[nodeID1].length()-containmentOverhangRight-1;
				}
				else if (node==nodeID2) {
					start=0;
					end=nodesArray[nodeID2].length()-1;
				}
			}
			if (start!=Math.NEGATIVE_INFINITY && end!=Math.NEGATIVE_INFINITY) {
				if (start>=end) {
					System.err.println("projectOnto_insertion> ERROR start="+start+" end="+end+" containmentOverhangLeft="+containmentOverhangLeft+" containmentOverhangRight="+containmentOverhangRight);			
					System.err.println("node: "+node+": "+nodesArray[node]);
					System.exit(1);
				}
				out[++pos]=Math.max(start,0); 
				out[++pos]=Math.min(end,nodesArray[node].length()-1);
			}
			return pos;
		}
		
		
		public final int projectOnto_overlap(int node, int[] out, int overlapType, int last) {
			int pos = last;
			int start=Math.NEGATIVE_INFINITY, end=Math.NEGATIVE_INFINITY;
			
			if (overhangs==null) return pos;
			if (overlapType==Constants.OVERLAP_PREFIX_PREFIX) {
				if (node==nodeID1) {
					start=0;
					end=nodesArray[nodeID1].length()-overhangs[0][0]-1;
				}
				else if (node==nodeID2) {
					start=0;
					end=nodesArray[nodeID2].length()-overhangs[0][1]-1;
				}
			}
			else if (overlapType==Constants.OVERLAP_PREFIX_SUFFIX) {
				if (node==nodeID1) {
					start=0;
					end=nodesArray[nodeID1].length()-overhangs[1][0]-1;
				}
				else if (node==nodeID2) {
					start=overhangs[1][1];
					end=nodesArray[nodeID2].length()-1;
				}
			}
			else if (overlapType==Constants.OVERLAP_SUFFIX_PREFIX) {
				if (node==nodeID1) {
					start=overhangs[2][0];
					end=nodesArray[nodeID1].length()-1;
				}
				else if (node==nodeID2) {
					start=0;
					end=nodesArray[nodeID2].length()-overhangs[2][1]-1;
				}
			}
			else if (overlapType==Constants.OVERLAP_SUFFIX_SUFFIX) {
				if (node==nodeID1) {
					start=overhangs[3][0];
					end=nodesArray[nodeID1].length()-1;
				}
				else if (node==nodeID2) {
					start=overhangs[3][1];
					end=nodesArray[nodeID2].length()-1;
				}
			}
			if (start!=Math.NEGATIVE_INFINITY && end!=Math.NEGATIVE_INFINITY) {
				if (start>=end) {
					System.err.println("projectOnto_overlap> ERROR start="+start+" end="+end+" edge="+this);
					System.err.println("node1 id="+nodeID1+": "+nodesArray[nodeID1]);
					System.err.println("node1.length="+nodesArray[nodeID1].length());
					System.err.println("node2 id="+nodeID2+": "+nodesArray[nodeID2]);
					System.err.println("node2.length="+nodesArray[nodeID2].length());
					for (int x=0; x<4; x++) System.err.println(overhangs[x][0]+","+overhangs[x][1]);		
					System.exit(1);
				}
				
				out[++pos]=Math.max(start,0); 
				out[++pos]=Math.min(end,nodesArray[node].length()-1);
			}	
			return pos;
		}
		
		
		public final int projectOnto_sharedSubstring(int node, int[] out, int last) {
			int pos = last;
			int start=Math.NEGATIVE_INFINITY, end=Math.NEGATIVE_INFINITY;
			
			if (sharedSubstringOverhangs==null) return pos;
			if (sharedSubstring!=-1) {
				if (node==nodeID1) {
					start=sharedSubstringOverhangs[0];
					end=nodesArray[nodeID1].length()-sharedSubstringOverhangs[1]-1;
				}
				else if (node==nodeID2) {
					start=sharedSubstringOverhangs[2];
					end=nodesArray[nodeID2].length()-sharedSubstringOverhangs[3]-1;
				}
			}
			if (start!=Math.NEGATIVE_INFINITY && end!=Math.NEGATIVE_INFINITY) {
				if (start>=end) {
					System.err.println("projectOnto_sharedSubstring> ERROR start="+start+" end="+end+" edge="+this);	
					System.err.println(sharedSubstringOverhangs[0]+","+sharedSubstringOverhangs[1]);
					System.err.println(sharedSubstringOverhangs[2]+","+sharedSubstringOverhangs[3]);
					System.err.println("node1="+nodeID1+": "+nodesArray[nodeID1]);
					System.err.println("node2="+nodeID2+": "+nodesArray[nodeID2]);
					System.err.println("edge info: "+this);
					System.exit(1);
				}
				out[++pos]=Math.max(start,0); 
				out[++pos]=Math.min(end,nodesArray[node].length()-1);
			}
			return pos;
		}
		
		
		/**
		 * If the new interval is identical to the projection on the old interval, and if
		 * the edge was an insertion or a containment, the edge is transformed into an 
		 * identity.
		 */
		private final void transform_identical_containmentOrInsertion() {
			setType_noContainment();
			setType(Constants.CONTAINMENT_IDENTICAL,true);
			containmentSubtype=-1;
			setType_noInsertion();
			containmentOverhangLeft=-1;
			containmentOverhangRight=-1;
		}
		
		
		/**
		 * If the new interval is identical to the projection on the old interval, and if
		 * the edge was an overlap, the edge is transformed into an insertion or a
		 * containment.
		 */
		private final void transform_identical_overlap(int oldNode, int newNode) {
			int i = getOverhangID();
			
			if (oldNode==nodeID1) {
				if (isType(Constants.OVERLAP_PREFIX_PREFIX) || isType(Constants.OVERLAP_SUFFIX_PREFIX)) {
					if (IntervalGraphStep2.isContained_alignment(nodeID2,newNode,nodesArray[nodeID2].start,nodesArray[nodeID2].end-overhangs[i][1],nodesArray[newNode].start,nodesArray[newNode].end,nodesArray[nodeID2].isLeftMaximal,nodesArray[nodeID2].isRightMaximal,nodesArray[newNode].isLeftMaximal,nodesArray[newNode].isRightMaximal,orientation==0?true:false)) {
						setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
						containmentSubtype=Constants.CONTAINMENT_PREFIX;
						containmentOverhangLeft=0;
						containmentOverhangRight=overhangs[i][1];
					}
					else if ( (isType(Constants.OVERLAP_PREFIX_PREFIX) && nodesArray[nodeID1].isLeftMaximal) ||
						 	  (isType(Constants.OVERLAP_SUFFIX_PREFIX) && nodesArray[nodeID1].isRightMaximal)
					        ) {
						setType(Constants.INSERTION_ONE_IN_TWO,true);
						containmentSubtype=Constants.CONTAINMENT_PREFIX;
						containmentOverhangLeft=0;
						containmentOverhangRight=overhangs[i][1];
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=0;
						sharedSubstringOverhangs[1]=0;
						sharedSubstringOverhangs[2]=0;
						sharedSubstringOverhangs[3]=overhangs[i][1];
					}
				}
				else if (isType(Constants.OVERLAP_PREFIX_SUFFIX) || isType(Constants.OVERLAP_SUFFIX_SUFFIX)) {
					if (IntervalGraphStep2.isContained_alignment(nodeID2,newNode,nodesArray[nodeID2].start+overhangs[i][1],nodesArray[nodeID2].end,nodesArray[newNode].start,nodesArray[newNode].end,nodesArray[nodeID2].isLeftMaximal,nodesArray[nodeID2].isRightMaximal,nodesArray[newNode].isLeftMaximal,nodesArray[newNode].isRightMaximal,orientation==0?true:false)) {
						setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
						containmentSubtype=Constants.CONTAINMENT_SUFFIX;
						containmentOverhangLeft=overhangs[i][1];
						containmentOverhangRight=0;
					}
					else if ( (isType(Constants.OVERLAP_PREFIX_SUFFIX) && nodesArray[nodeID1].isLeftMaximal) ||
						 	  (isType(Constants.OVERLAP_SUFFIX_SUFFIX) && nodesArray[nodeID1].isRightMaximal)
					   		) {
						setType(Constants.INSERTION_ONE_IN_TWO,true);
						containmentSubtype=Constants.CONTAINMENT_SUFFIX;
						containmentOverhangLeft=overhangs[i][1];
						containmentOverhangRight=0;
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=0;
						sharedSubstringOverhangs[1]=0;
						sharedSubstringOverhangs[2]=overhangs[i][1];
						sharedSubstringOverhangs[3]=0;
					}
				}
			}
			else {
				if (isType(Constants.OVERLAP_PREFIX_PREFIX) || isType(Constants.OVERLAP_PREFIX_SUFFIX)) {
					if (IntervalGraphStep2.isContained_alignment(nodeID1,newNode,nodesArray[nodeID1].start,nodesArray[nodeID1].end-overhangs[i][0],nodesArray[newNode].start,nodesArray[newNode].end,nodesArray[nodeID1].isLeftMaximal,nodesArray[nodeID1].isRightMaximal,nodesArray[newNode].isLeftMaximal,nodesArray[newNode].isRightMaximal,orientation==0?true:false)) {
						setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
						containmentSubtype=Constants.CONTAINMENT_PREFIX;
						containmentOverhangLeft=0;
						containmentOverhangRight=overhangs[i][0];
					}
					else if ( (isType(Constants.OVERLAP_PREFIX_PREFIX) && nodesArray[nodeID2].isLeftMaximal) ||
						 	  (isType(Constants.OVERLAP_PREFIX_SUFFIX) && nodesArray[nodeID2].isRightMaximal)
					   		) {
						setType(Constants.INSERTION_TWO_IN_ONE,true);
						containmentSubtype=Constants.CONTAINMENT_PREFIX;
						containmentOverhangLeft=0;
						containmentOverhangRight=overhangs[i][0];
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=0;
						sharedSubstringOverhangs[1]=overhangs[i][0];
						sharedSubstringOverhangs[2]=0;
						sharedSubstringOverhangs[3]=0;
					}
				}
				else if (isType(Constants.OVERLAP_SUFFIX_PREFIX) || isType(Constants.OVERLAP_SUFFIX_SUFFIX)) {
					if (IntervalGraphStep2.isContained_alignment(nodeID1,newNode,nodesArray[nodeID1].start+overhangs[i][0],nodesArray[nodeID1].end,nodesArray[newNode].start,nodesArray[newNode].end,nodesArray[nodeID1].isLeftMaximal,nodesArray[nodeID1].isRightMaximal,nodesArray[newNode].isLeftMaximal,nodesArray[newNode].isRightMaximal,orientation==0?true:false)) {
						setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
						containmentSubtype=Constants.CONTAINMENT_SUFFIX;
						containmentOverhangLeft=overhangs[i][0];
						containmentOverhangRight=0;
					}
					else if ( (isType(Constants.OVERLAP_SUFFIX_PREFIX) && nodesArray[nodeID2].isLeftMaximal) ||
						 	  (isType(Constants.OVERLAP_SUFFIX_SUFFIX) && nodesArray[nodeID2].isRightMaximal)
					   		) {
						setType(Constants.INSERTION_TWO_IN_ONE,true);
						containmentSubtype=Constants.CONTAINMENT_SUFFIX;
						containmentOverhangLeft=overhangs[i][0];
						containmentOverhangRight=0;
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=overhangs[i][0];
						sharedSubstringOverhangs[1]=0;
						sharedSubstringOverhangs[2]=0;
						sharedSubstringOverhangs[3]=0;
					}
				}
			}
			setType_noOverlap();
		}
		
		
		/**
		 * If the new interval is identical to the projection on the old interval, and if
		 * the edge was a shared substring, the edge is transformed into an insertion or a
		 * containment.
		 */
		private final void transform_identical_sharedSubstring(int oldNode, int newNode, int newNodeStart, int newNodeEnd) {
			int otherStart, otherEnd;
			
			if (oldNode==nodeID1) {
				otherStart=sharedSubstringOverhangs[2];
				otherEnd=nodesArray[nodeID2].length()-sharedSubstringOverhangs[3]-1;
				if (IntervalGraphStep2.isContained_alignment(nodeID2,newNode,nodesArray[nodeID2].start+otherStart,nodesArray[nodeID2].start+otherEnd,nodesArray[newNode].start,nodesArray[newNode].end,nodesArray[nodeID2].isLeftMaximal,nodesArray[nodeID2].isRightMaximal,nodesArray[newNode].isLeftMaximal,nodesArray[newNode].isRightMaximal,orientation==0?true:false)) {
					setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
					if (otherStart<=IDENTITY_THRESHOLD) containmentSubtype=Constants.CONTAINMENT_PREFIX;
					else if (otherEnd>=nodesArray[nodeID2].length()-IDENTITY_THRESHOLD) containmentSubtype=Constants.CONTAINMENT_SUFFIX;
					else containmentSubtype=Constants.CONTAINMENT_SUBSTRING;
					containmentOverhangLeft=otherStart;
					containmentOverhangRight=sharedSubstringOverhangs[3];
					setType_noSharedSubstring();
				}
				else if ( (newNodeStart>IDENTITY_THRESHOLD && newNodeEnd<nodesArray[nodeID1].length()-IDENTITY_THRESHOLD) ||
						  ( newNodeStart<=IDENTITY_THRESHOLD &&
						    ( ( orientation==0 &&
						 	    ( otherEnd<nodesArray[nodeID2].length()-IDENTITY_THRESHOLD ||
							      nodesArray[nodeID1].isLeftMaximal
						        )
						      ) ||
							  ( orientation==1 &&
							    ( otherStart>IDENTITY_THRESHOLD ||
							      nodesArray[nodeID1].isLeftMaximal
						        )
						      )
						    )
						  ) ||
						  ( newNodeEnd>=nodesArray[nodeID1].length()-IDENTITY_THRESHOLD && 
						    ( ( orientation==0 &&
							    ( otherStart>IDENTITY_THRESHOLD ||
							      nodesArray[nodeID1].isRightMaximal
						        )
						      ) ||
							  ( orientation==1 &&
							    ( otherEnd<nodesArray[nodeID2].length()-IDENTITY_THRESHOLD ||
							      nodesArray[nodeID1].isRightMaximal
						        )
						      )
						    )
						  )
					    ) {
					setType(Constants.INSERTION_ONE_IN_TWO,true);
					if (otherStart<=IDENTITY_THRESHOLD) containmentSubtype=Constants.CONTAINMENT_PREFIX;
					else if (otherEnd>=nodesArray[nodeID2].length()-IDENTITY_THRESHOLD) containmentSubtype=Constants.CONTAINMENT_SUFFIX;
					else containmentSubtype=Constants.CONTAINMENT_SUBSTRING;
					containmentOverhangLeft=otherStart;
					containmentOverhangRight=sharedSubstringOverhangs[3];
					setType_noSharedSubstring();
				}
				else {
					sharedSubstringOverhangs[0]=0;
					sharedSubstringOverhangs[1]=0;
				}
			}
			else if (oldNode==nodeID2) {
				otherStart=sharedSubstringOverhangs[0];
				otherEnd=nodesArray[nodeID1].length()-sharedSubstringOverhangs[1]-1;
				if (IntervalGraphStep2.isContained_alignment(nodeID1,newNode,nodesArray[nodeID1].start+otherStart,nodesArray[nodeID1].start+otherEnd,nodesArray[newNode].start,nodesArray[newNode].end,nodesArray[nodeID1].isLeftMaximal,nodesArray[nodeID1].isRightMaximal,nodesArray[newNode].isLeftMaximal,nodesArray[newNode].isRightMaximal,orientation==0?true:false)) {
					setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
					if (otherStart<=IDENTITY_THRESHOLD) containmentSubtype=Constants.CONTAINMENT_PREFIX;
					else if (otherEnd>=nodesArray[nodeID1].length()-IDENTITY_THRESHOLD) containmentSubtype=Constants.CONTAINMENT_SUFFIX;
					else containmentSubtype=Constants.CONTAINMENT_SUBSTRING;
					containmentOverhangLeft=otherStart;
					containmentOverhangRight=sharedSubstringOverhangs[1];
					setType_noSharedSubstring();
				}
				else if ( (newNodeStart>IDENTITY_THRESHOLD && newNodeEnd<nodesArray[nodeID2].length()-IDENTITY_THRESHOLD) ||
						  ( newNodeStart<=IDENTITY_THRESHOLD && 
						    ( ( orientation==0 &&
						 	    ( otherEnd<nodesArray[nodeID1].length()-IDENTITY_THRESHOLD ||
							      nodesArray[nodeID2].isLeftMaximal
						        )
						      ) ||
							  ( orientation==1 &&
							    ( otherStart>IDENTITY_THRESHOLD ||
							      nodesArray[nodeID2].isLeftMaximal
						        )
						      )
						    )
						  ) ||
						  ( newNodeEnd>=nodesArray[nodeID2].length()-IDENTITY_THRESHOLD && 
						    ( ( orientation==0 &&
							    ( otherStart>IDENTITY_THRESHOLD ||
							      nodesArray[nodeID2].isRightMaximal
						        )
						      ) ||
							  ( orientation==1 &&
							    ( otherEnd<nodesArray[nodeID1].length()-IDENTITY_THRESHOLD ||
							      nodesArray[nodeID2].isRightMaximal
						        )
						      )
						    )
						  )
					    ) {
					setType(Constants.INSERTION_TWO_IN_ONE,true);
					if (otherStart<=IDENTITY_THRESHOLD) containmentSubtype=Constants.CONTAINMENT_PREFIX;
					else if (otherEnd>=nodesArray[nodeID1].length()-IDENTITY_THRESHOLD) containmentSubtype=Constants.CONTAINMENT_SUFFIX;
					else containmentSubtype=Constants.CONTAINMENT_SUBSTRING;
					containmentOverhangLeft=otherStart;
					containmentOverhangRight=sharedSubstringOverhangs[1];
					setType_noSharedSubstring();
				}
				else {
					sharedSubstringOverhangs[2]=0;
					sharedSubstringOverhangs[3]=0;
				}
			}
		}
		
		
		/**
		 * If the projection on the old interval is contained in the new interval, and if
		 * the edge was an insertion or a containment, the edge is transformed into an 
		 * insertion or a containment.
		 */
		private final void transform_contained_containmentOrInsertion(int oldNode, int newNode, int newNodeStart, int newNodeEnd, int projectionStart, int projectionEnd) {
			final int otherNode = oldNode==nodeID1?nodeID2:nodeID1;			
			
			if (IntervalGraphStep2.isContained_alignment(otherNode,newNode,nodesArray[otherNode].start,nodesArray[otherNode].end,nodesArray[newNode].start+(projectionStart>newNodeStart?projectionStart-newNodeStart:0),nodesArray[newNode].end-(newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0),nodesArray[otherNode].isLeftMaximal,nodesArray[otherNode].isRightMaximal,nodesArray[newNode].isLeftMaximal,nodesArray[newNode].isRightMaximal,orientation==0?true:false)) {
				setType_noInsertion();
				if (oldNode==nodeID1) setType(Constants.CONTAINMENT_TWO_IN_ONE,true);
				else setType(Constants.CONTAINMENT_ONE_IN_TWO,true);
				if (Math.abs(projectionStart,newNodeStart)<=IDENTITY_THRESHOLD) {
					containmentSubtype=Constants.CONTAINMENT_PREFIX;
					containmentOverhangLeft=0;
					containmentOverhangRight=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
				}
				else if (Math.abs(projectionEnd,newNodeEnd)<=IDENTITY_THRESHOLD) {
					containmentSubtype=Constants.CONTAINMENT_SUFFIX;
					containmentOverhangLeft=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
					containmentOverhangRight=0;
				}
				else {
					containmentSubtype=Constants.CONTAINMENT_SUBSTRING;
					containmentOverhangLeft=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
					containmentOverhangRight=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
				}
			}
			else if ( (projectionStart>newNodeStart+IDENTITY_THRESHOLD && projectionEnd<newNodeEnd-IDENTITY_THRESHOLD && (nodesArray[otherNode].isLeftMaximal||nodesArray[otherNode].isRightMaximal)) ||
					  ( orientation==0 && 
					    ( (projectionStart<=newNodeStart+IDENTITY_THRESHOLD && nodesArray[otherNode].isRightMaximal) ||
					      (projectionEnd>=newNodeEnd-IDENTITY_THRESHOLD && nodesArray[otherNode].isLeftMaximal)
					    )
					  ) ||
					  ( orientation==1 && 
	 				    ( (projectionStart<=newNodeStart+IDENTITY_THRESHOLD && nodesArray[otherNode].isLeftMaximal) ||
	 				      (projectionEnd>=newNodeEnd-IDENTITY_THRESHOLD && nodesArray[otherNode].isRightMaximal)
	 				    )
	 				  )
				    ) {
				if (oldNode==nodeID1) setType(Constants.INSERTION_TWO_IN_ONE,true);
				else setType(Constants.INSERTION_ONE_IN_TWO,true);
				setType_noContainment();
				if (Math.abs(projectionStart,newNodeStart)<=IDENTITY_THRESHOLD) {
					containmentSubtype=Constants.CONTAINMENT_PREFIX;
					containmentOverhangLeft=0;
					containmentOverhangRight=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
				}
				else if (Math.abs(projectionEnd,newNodeEnd)<=IDENTITY_THRESHOLD) {
					containmentSubtype=Constants.CONTAINMENT_SUFFIX;
					containmentOverhangLeft=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
					containmentOverhangRight=0;
				}
				else {
					containmentSubtype=Constants.CONTAINMENT_SUBSTRING;
					containmentOverhangLeft=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
					containmentOverhangRight=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
				}
			}
			else {
				setType(Constants.SHARED_SUBSTRING,true);
				setType_noContainment();
				setType_noInsertion();
				if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
				if (oldNode==nodeID1) {
					sharedSubstringOverhangs[0]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
					sharedSubstringOverhangs[1]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
					sharedSubstringOverhangs[2]=0;
					sharedSubstringOverhangs[3]=0;
				}
				else {
					sharedSubstringOverhangs[0]=0;
					sharedSubstringOverhangs[1]=0;
					sharedSubstringOverhangs[2]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
					sharedSubstringOverhangs[3]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
				}
			}			
		}
		
		
		/**
		 * If the projection on the old interval is contained in the new interval, and
		 * if the edge was an overlap, the edge is transformed into an overlap or a shared
		 * substring.
		 */
		private final void transform_contained_overlap(int oldNode, int newNode, int newNodeStart, int newNodeEnd, int projectionStart, int projectionEnd, AlignmentPair tmpPair) {
			int i, j;
			
			tmpPair.orientation=orientation==0?true:false;
			i=getOverhangID();
			if (oldNode==nodeID1) {
				tmpPair.node1=nodesArray[newNode];
				tmpPair.node2=nodesArray[nodeID2];
				if (isType(Constants.OVERLAP_PREFIX_PREFIX)) {
					tmpPair.alignmentStart1=nodesArray[newNode].start;
					tmpPair.alignmentEnd1=nodesArray[newNode].start+projectionEnd;
					tmpPair.alignmentStart2=nodesArray[nodeID2].start;
					tmpPair.alignmentEnd2=nodesArray[nodeID2].end-overhangs[i][1];
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						overhangs[i][0]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
						if ((flexibleOverlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 && maxOverhangs!=null) maxOverhangs[i][0]=nodesArray[newNode].length()-(nodesArray[oldNode].length()-maxOverhangs[i][0]);
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=0;
						sharedSubstringOverhangs[1]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
						sharedSubstringOverhangs[2]=0;
						sharedSubstringOverhangs[3]=overhangs[i][1];
						setType_noOverlap();
					}
				}
				else if (isType(Constants.OVERLAP_PREFIX_SUFFIX)) {
					tmpPair.alignmentStart1=nodesArray[newNode].start;
					tmpPair.alignmentEnd1=nodesArray[newNode].start+projectionEnd;
					tmpPair.alignmentStart2=nodesArray[nodeID2].start+overhangs[i][1];
					tmpPair.alignmentEnd2=nodesArray[nodeID2].end;
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						overhangs[i][0]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
						if ((flexibleOverlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && maxOverhangs!=null) maxOverhangs[i][0]=nodesArray[newNode].length()-(nodesArray[oldNode].length()-maxOverhangs[i][0]);
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=0;
						sharedSubstringOverhangs[1]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
						sharedSubstringOverhangs[2]=overhangs[i][1];
						sharedSubstringOverhangs[3]=0;
						setType_noOverlap();
					}
				}
				else if (isType(Constants.OVERLAP_SUFFIX_PREFIX)) {
					tmpPair.alignmentStart1=nodesArray[newNode].start+projectionStart-newNodeStart;
					tmpPair.alignmentEnd1=nodesArray[newNode].end;
					tmpPair.alignmentStart2=nodesArray[nodeID2].start;
					tmpPair.alignmentEnd2=nodesArray[nodeID2].end-overhangs[i][1];
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						overhangs[i][0]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
						if ((flexibleOverlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && maxOverhangs!=null) maxOverhangs[i][0]=nodesArray[newNode].length()-(nodesArray[oldNode].length()-maxOverhangs[i][0]);
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
						sharedSubstringOverhangs[1]=0;
						sharedSubstringOverhangs[2]=0;
						sharedSubstringOverhangs[3]=overhangs[i][1];
						setType_noOverlap();
					}
				}
				else if (isType(Constants.OVERLAP_SUFFIX_SUFFIX)) {
					tmpPair.alignmentStart1=nodesArray[newNode].start+projectionStart-newNodeStart;
					tmpPair.alignmentEnd1=nodesArray[newNode].end;
					tmpPair.alignmentStart2=nodesArray[nodeID2].start+overhangs[i][1];
					tmpPair.alignmentEnd2=nodesArray[nodeID2].end;
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						overhangs[i][0]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
						if ((flexibleOverlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0 && maxOverhangs!=null) maxOverhangs[i][0]=nodesArray[newNode].length()-(nodesArray[oldNode].length()-maxOverhangs[i][0]);
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
						sharedSubstringOverhangs[1]=0;
						sharedSubstringOverhangs[2]=overhangs[i][1];
						sharedSubstringOverhangs[3]=0;
						setType_noOverlap();
					}
				}
			}
			else if (oldNode==nodeID2) {
				tmpPair.node1=nodesArray[nodeID1];
				tmpPair.node2=nodesArray[newNode];
				if (isType(Constants.OVERLAP_PREFIX_PREFIX)) {
					tmpPair.alignmentStart1=nodesArray[nodeID1].start;
					tmpPair.alignmentEnd1=nodesArray[nodeID1].end-overhangs[i][0];
					tmpPair.alignmentStart2=nodesArray[newNode].start;
					tmpPair.alignmentEnd2=nodesArray[newNode].start+projectionEnd;
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						overhangs[i][1]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
						if ((flexibleOverlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 && maxOverhangs!=null) maxOverhangs[i][1]=nodesArray[newNode].length()-(nodesArray[oldNode].length()-maxOverhangs[i][1]);
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=0;
						sharedSubstringOverhangs[1]=overhangs[i][0];
						sharedSubstringOverhangs[2]=0;
						sharedSubstringOverhangs[3]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
						setType_noOverlap();
					}
				}
				else if (isType(Constants.OVERLAP_PREFIX_SUFFIX)) {
					tmpPair.alignmentStart1=nodesArray[nodeID1].start;
					tmpPair.alignmentEnd1=nodesArray[nodeID1].end-overhangs[i][0];
					tmpPair.alignmentStart2=nodesArray[newNode].start+projectionStart-newNodeStart;
					tmpPair.alignmentEnd2=nodesArray[newNode].end;
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						overhangs[i][1]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
						if ((flexibleOverlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && maxOverhangs!=null) maxOverhangs[i][1]=nodesArray[newNode].length()-(nodesArray[oldNode].length()-maxOverhangs[i][1]);
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=0;
						sharedSubstringOverhangs[1]=overhangs[i][0];
						sharedSubstringOverhangs[2]=projectionEnd>newNodeEnd?projectionEnd-newNodeEnd:0;
						sharedSubstringOverhangs[3]=0;
						setType_noOverlap();
					}
				}
				else if (isType(Constants.OVERLAP_SUFFIX_PREFIX)) {
					tmpPair.alignmentStart1=nodesArray[nodeID1].start+overhangs[i][0];
					tmpPair.alignmentEnd1=nodesArray[nodeID1].end;
					tmpPair.alignmentStart2=nodesArray[newNode].start;
					tmpPair.alignmentEnd2=nodesArray[newNode].start+projectionEnd;
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						overhangs[i][1]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
						if ((flexibleOverlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && maxOverhangs!=null) maxOverhangs[i][1]=nodesArray[newNode].length()-(nodesArray[oldNode].length()-maxOverhangs[i][1]);
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=overhangs[i][0];
						sharedSubstringOverhangs[1]=0;
						sharedSubstringOverhangs[2]=0;
						sharedSubstringOverhangs[3]=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
						setType_noOverlap();
					}
				}
				else if (isType(Constants.OVERLAP_SUFFIX_SUFFIX)) {
					tmpPair.alignmentStart1=nodesArray[nodeID1].start+overhangs[i][0];
					tmpPair.alignmentEnd1=nodesArray[nodeID1].end;
					tmpPair.alignmentStart2=nodesArray[newNode].start+projectionStart-newNodeStart;
					tmpPair.alignmentEnd2=nodesArray[newNode].end;
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						overhangs[i][1]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
						if ((flexibleOverlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0 && maxOverhangs!=null) maxOverhangs[i][1]=nodesArray[newNode].length()-(nodesArray[oldNode].length()-maxOverhangs[i][1]);
					}
					else {
						setType(Constants.SHARED_SUBSTRING,true);
						if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
						sharedSubstringOverhangs[0]=overhangs[i][0];
						sharedSubstringOverhangs[1]=0;
						sharedSubstringOverhangs[2]=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
						sharedSubstringOverhangs[3]=0;
						setType_noOverlap();
					}
				}
			}
		}
		
		
		/**
		 * If the projection on the old interval is contained in the new interval, and
		 * if the edge was a shared substring, the edge is transformed into a shared 
		 * substring or an overlap.
		 */
		private final void transform_contained_sharedSubstring(int oldNode, int newNode, int newNodeStart, int newNodeEnd, int projectionStart, int projectionEnd, AlignmentPair tmpPair) {
			final int newOverhangLeft, newOverhangRight;
			
			tmpPair.orientation=orientation==0?true:false;
			newOverhangLeft=projectionStart>newNodeStart?projectionStart-newNodeStart:0;
			newOverhangRight=newNodeEnd>projectionEnd?newNodeEnd-projectionEnd:0;
			if (Math.abs(projectionStart,newNodeStart)<=IDENTITY_THRESHOLD) {
				if (oldNode==nodeID1) {
					tmpPair.node1=nodesArray[newNode];
					tmpPair.node2=nodesArray[nodeID2];
					tmpPair.alignmentStart1=nodesArray[newNode].start+newOverhangLeft;
					tmpPair.alignmentEnd1=nodesArray[newNode].end-newOverhangRight;
					tmpPair.alignmentStart2=nodesArray[nodeID2].start+sharedSubstringOverhangs[2];
					tmpPair.alignmentEnd2=nodesArray[nodeID2].end-sharedSubstringOverhangs[3];
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						if (sharedSubstringOverhangs[2]<=IDENTITY_THRESHOLD) {
							setType(Constants.OVERLAP_PREFIX_PREFIX,true);
							if (overhangs==null) {
								overhangs = new int[4][2];
								Math.set(overhangs,-1);
							}
							overhangs[0][0]=newOverhangRight;
							overhangs[0][1]=sharedSubstringOverhangs[3];
							flexibleOverlap=0;
							setType_noSharedSubstring();
						}
						else if (sharedSubstringOverhangs[3]<=IDENTITY_THRESHOLD) {
							setType(Constants.OVERLAP_PREFIX_SUFFIX,true);
							if (overhangs==null) {
								overhangs = new int[4][2];
								Math.set(overhangs,-1);
							}
							overhangs[1][0]=newOverhangRight;
							overhangs[1][1]=sharedSubstringOverhangs[2];
							flexibleOverlap=0;
							setType_noSharedSubstring();
						}
					}
					else {
						sharedSubstringOverhangs[0]=newOverhangLeft;
						sharedSubstringOverhangs[1]=newOverhangRight;
					}
				}
				else if (oldNode==nodeID2) {
					tmpPair.node1=nodesArray[nodeID1];
					tmpPair.node2=nodesArray[newNode];
					tmpPair.alignmentStart1=nodesArray[nodeID1].start+sharedSubstringOverhangs[0];
					tmpPair.alignmentEnd1=nodesArray[nodeID1].end-sharedSubstringOverhangs[1];
					tmpPair.alignmentStart2=nodesArray[newNode].start+newOverhangLeft;
					tmpPair.alignmentEnd2=nodesArray[newNode].end-newOverhangRight;
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						if (sharedSubstringOverhangs[0]<=IDENTITY_THRESHOLD) {
							setType(Constants.OVERLAP_PREFIX_PREFIX,true);
							if (overhangs==null) {
								overhangs = new int[4][2];
								Math.set(overhangs,-1);
							}
							overhangs[0][0]=sharedSubstringOverhangs[1];
							overhangs[0][1]=newOverhangRight;
							flexibleOverlap=0;
							setType_noSharedSubstring();
						}
						else if (sharedSubstringOverhangs[1]<=IDENTITY_THRESHOLD) {
							setType(Constants.OVERLAP_SUFFIX_PREFIX,true);
							if (overhangs==null) {
								overhangs = new int[4][2];
								Math.set(overhangs,-1);
							}
							overhangs[2][0]=sharedSubstringOverhangs[0];
							overhangs[2][1]=newOverhangRight;
							flexibleOverlap=0;
							setType_noSharedSubstring();
						}
					}
					else {
						sharedSubstringOverhangs[2]=newOverhangLeft;
						sharedSubstringOverhangs[3]=newOverhangRight;
					}
				}
			}
			else if (Math.abs(projectionEnd,newNodeEnd)<=IDENTITY_THRESHOLD) {
				if (oldNode==nodeID1) {
					tmpPair.node1=nodesArray[newNode];
					tmpPair.node2=nodesArray[nodeID2];
					tmpPair.alignmentStart1=nodesArray[newNode].start+newOverhangLeft;
					tmpPair.alignmentEnd1=nodesArray[newNode].end-newOverhangRight;
					tmpPair.alignmentStart2=nodesArray[nodeID2].start+sharedSubstringOverhangs[2];
					tmpPair.alignmentEnd2=nodesArray[nodeID2].end-sharedSubstringOverhangs[3];
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						if (sharedSubstringOverhangs[2]<=IDENTITY_THRESHOLD) {
							setType(Constants.OVERLAP_SUFFIX_PREFIX,true);
							if (overhangs==null) {
								overhangs = new int[4][2];
								Math.set(overhangs,-1);
							}
							overhangs[2][0]=newOverhangLeft;
							overhangs[2][1]=sharedSubstringOverhangs[3];
							flexibleOverlap=0;
							setType_noSharedSubstring();
						}
						else if (sharedSubstringOverhangs[3]<=IDENTITY_THRESHOLD) {
							setType(Constants.OVERLAP_SUFFIX_SUFFIX,true);
							if (overhangs==null) {
								overhangs = new int[4][2];
								Math.set(overhangs,-1);
							}
							overhangs[3][0]=newOverhangLeft;
							overhangs[3][1]=sharedSubstringOverhangs[2];
							flexibleOverlap=0;
							setType_noSharedSubstring();
						}
					}
					else {
						sharedSubstringOverhangs[0]=newOverhangLeft;
						sharedSubstringOverhangs[1]=newOverhangRight;
					}
				}
				else if (oldNode==nodeID2) {
					tmpPair.node1=nodesArray[nodeID1];
					tmpPair.node2=nodesArray[newNode];
					tmpPair.alignmentStart1=nodesArray[nodeID1].start+sharedSubstringOverhangs[0];
					tmpPair.alignmentEnd1=nodesArray[nodeID1].end-sharedSubstringOverhangs[1];
					tmpPair.alignmentStart2=nodesArray[newNode].start+newOverhangLeft;
					tmpPair.alignmentEnd2=nodesArray[newNode].end-newOverhangRight;
					if (IntervalGraphStep2.inStep2_overlap(tmpPair)!=-1) {
						if (sharedSubstringOverhangs[0]<=IDENTITY_THRESHOLD) {
							setType(Constants.OVERLAP_PREFIX_SUFFIX,true);
							if (overhangs==null) {
								overhangs = new int[4][2];
								Math.set(overhangs,-1);
							}
							overhangs[1][0]=sharedSubstringOverhangs[1];
							overhangs[1][1]=newOverhangLeft;
							flexibleOverlap=0;
							setType_noSharedSubstring();
						}
						else if (sharedSubstringOverhangs[1]<=IDENTITY_THRESHOLD) {
							setType(Constants.OVERLAP_SUFFIX_SUFFIX,true);
							if (overhangs==null) {
								overhangs = new int[4][2];
								Math.set(overhangs,-1);
							}
							overhangs[3][0]=sharedSubstringOverhangs[0];
							overhangs[3][1]=newOverhangLeft;
							flexibleOverlap=0;
							setType_noSharedSubstring();
						}
					}
					else {
						sharedSubstringOverhangs[2]=newOverhangLeft;
						sharedSubstringOverhangs[3]=newOverhangRight;
					}
				}
			}
			else {
				if (oldNode==nodeID1) {
					sharedSubstringOverhangs[0]=newOverhangLeft;
					sharedSubstringOverhangs[1]=newOverhangRight;
				}
				else if (oldNode==nodeID2) {
					sharedSubstringOverhangs[2]=newOverhangLeft;
					sharedSubstringOverhangs[3]=newOverhangRight;
				}
			}			
		}
		
		
		public final void serialize(BufferedOutputStream out) throws IOException {
			int i, c, mask;
			
			IO.writeInt(nodeID1,out);
			IO.writeInt(nodeID2,out);
			
			// Information for Step 2
			IO.writeInt(containment,out);
			IO.writeInt(containmentSubtype,out);
			IO.writeInt(overlap,out);
			IO.writeInt(insertion,out);
			IO.writeInt(sharedSubstring,out);
			if (overhangs==null) {
				for (i=0; i<4; i++) {
					IO.writeInt(-1,out);
					IO.writeInt(-1,out);
				}
			}
			else {
				for (i=0; i<4; i++) {
					IO.writeInt(overhangs[i][0],out);
					IO.writeInt(overhangs[i][1],out);
				}
			}
			IO.writeInt(containmentOverhangLeft,out);
			IO.writeInt(containmentOverhangRight,out);
			if (sharedSubstringOverhangs==null) for (i=0; i<4; i++) IO.writeInt(-1,out);
			else {
				for (i=0; i<4; i++) IO.writeInt(sharedSubstringOverhangs[i],out);
			}
			c=0x00000000; mask=0x00000001;
			if (isRedundant) c|=mask;
			mask<<=1;
			if (on) c|=mask;
			mask<<=1;
			if (implied) c|=mask;
			mask<<=1;
			if (supplement) c|=mask;
			IO.writeInt(c,out);
			IO.writeInt(nAlignmentsForward,out);
			IO.writeInt(nAlignmentsBackward,out);
			IO.writeInt(orientation,out);
			IO.writeDouble(avgDiffs,out);
		}
		
		
		public final void deserialize(BufferedInputStream in) throws IOException {
			int i, c, mask;
			
			nodeID1=IO.readInt(in);
			nodeID2=IO.readInt(in);
			
			// Information for Step 2
			containment=IO.readInt(in);
			containmentSubtype=IO.readInt(in);
			overlap=IO.readInt(in);
			insertion=IO.readInt(in);
			sharedSubstring=IO.readInt(in);
			if (overhangs==null) overhangs = new int[4][2];
			for (i=0; i<4; i++) {
				overhangs[i][0]=IO.readInt(in);
				overhangs[i][1]=IO.readInt(in);
			}
			containmentOverhangLeft=IO.readInt(in);
			containmentOverhangRight=IO.readInt(in);
			if (sharedSubstringOverhangs==null) sharedSubstringOverhangs = new int[4];
			for (i=0; i<4; i++) sharedSubstringOverhangs[i]=IO.readInt(in);
			c=IO.readInt(in); mask=0x00000001;
			isRedundant=(c&mask)!=0;
			mask<<=1;
			on=(c&mask)!=0;
			mask<<=1;
			implied=(c&mask)!=0;
			mask<<=1;
			supplement=(c&mask)!=0;
			nAlignmentsForward=IO.readInt(in);
			nAlignmentsBackward=IO.readInt(in);
			orientation=IO.readInt(in);
			avgDiffs=IO.readDouble(in);
		}
		
		
		/**
		 * Like $deserialize()$, but just advances the buffer by one edge.
		 */
		public static final void skip(BufferedInputStream in) throws IOException {
			int i;
			
			IO.readInt(in);
			IO.readInt(in);
			
			// Information for Step 2
			IO.readInt(in);
			IO.readInt(in);
			IO.readInt(in);
			IO.readInt(in);
			IO.readInt(in);
			for (i=0; i<4; i++) { IO.readInt(in); IO.readInt(in); }
			IO.readInt(in);
			IO.readInt(in);
			for (i=0; i<4; i++) IO.readInt(in);
			IO.readInt(in);
			IO.readInt(in);
			IO.readInt(in);
			IO.readInt(in);
			IO.readDouble(in);
		}
		
		
		public final void setType(int edgeType, boolean on) {
			switch (edgeType) {
				case Constants.OVERLAP_PREFIX_PREFIX: overlap=on?edgeType:-1; break;
				case Constants.OVERLAP_PREFIX_SUFFIX: overlap=on?edgeType:-1; break;
				case Constants.OVERLAP_SUFFIX_PREFIX: overlap=on?edgeType:-1; break;
				case Constants.OVERLAP_SUFFIX_SUFFIX: overlap=on?edgeType:-1; break;
				case Constants.CONTAINMENT_IDENTICAL: containment=on?edgeType:-1; break;
				case Constants.CONTAINMENT_ONE_IN_TWO: containment=on?edgeType:-1; break;
				case Constants.CONTAINMENT_TWO_IN_ONE: containment=on?edgeType:-1; break;
				case Constants.INSERTION_ONE_IN_TWO: insertion=on?edgeType:-1; break;
				case Constants.INSERTION_TWO_IN_ONE: insertion=on?edgeType:-1; break;
				case Constants.INSERTION_BOTH: insertion=on?edgeType:-1; break;
				case Constants.SHARED_SUBSTRING: sharedSubstring=on?edgeType:-1; break;
			}
		}


		public final boolean isType(int edgeType) {
			switch (edgeType) {
				case Constants.OVERLAP_PREFIX_PREFIX: return overlap==edgeType;
				case Constants.OVERLAP_PREFIX_SUFFIX: return overlap==edgeType;
				case Constants.OVERLAP_SUFFIX_PREFIX: return overlap==edgeType;
				case Constants.OVERLAP_SUFFIX_SUFFIX: return overlap==edgeType;
				case Constants.CONTAINMENT_IDENTICAL: return containment==edgeType;
				case Constants.CONTAINMENT_ONE_IN_TWO: return containment==edgeType;
				case Constants.CONTAINMENT_TWO_IN_ONE: return containment==edgeType;
				case Constants.INSERTION_ONE_IN_TWO: return insertion==edgeType;
				case Constants.INSERTION_TWO_IN_ONE: return insertion==edgeType;
				case Constants.INSERTION_BOTH: return insertion==edgeType;
				case Constants.SHARED_SUBSTRING: return sharedSubstring==edgeType;
			}
			return false;
		}


		public final boolean isType_overlap() {
			return overlap!=-1;
		}


		public final boolean isType_containment() {
			return containment!=-1;
		}


		public final boolean isType_insertion() {
			return insertion!=-1;
		}
		
		
		public final boolean isType_sharedSubstring() {
			return sharedSubstring!=-1;
		}
		
		
		/**
		 * @return the number of distinct types of the edge.
		 */
		public int nDistinctTypes() {
			int nTypes = 0;
			if (isType_containment()) nTypes++;
			if (isType_insertion()) nTypes++;
			if (isType_overlap()) nTypes++;
			if (isType_sharedSubstring()) nTypes++;
			return nTypes;
		}
		
		
		/**
		 * @return the number of types of the edge, counting sub-types of insertion and
		 * overlap, if any.
		 */
		public int nTypes() {
			int nTypes = 0;
			if (isType_containment()) nTypes++;
			if (isType_insertion()) nTypes+=nInsertionTypes();
			if (isType_overlap()) nTypes+=nOverlapTypes();
			if (isType_sharedSubstring()) nTypes++;
			return nTypes;
		}
		
		
		/**
		 * @return the number of alignments whose coordinates are stored in the edge.
		 */
		public int nStoredAlignments() {
			int out = 0;
			if (isType_containment()) out++;
			if (isType_insertion()) out++;
			if (isType_overlap()) out+=nStoredAlignments_overlap();
			if (isType_sharedSubstring()) out++;
			return out;
		}


		public final void setType_noOverlap() {
			overlap=-1;
			flexibleOverlap=0;
			overhangs=null; maxOverhangs=null;
		}

		
		public final void setType_noSharedSubstring() {
			sharedSubstring=-1;
			sharedSubstringOverhangs=null;
		}
		

		public final void setType_noContainment() {
			containment=-1;
		}


		public final void setType_noInsertion() {
			insertion=-1;
		}


		public final void copyType_containment(Edge other) {
			containment=other.containment;
		}

		
		/**
		 * Overwrites the overlap type of this edge with the one of $other$.
		 *
		 * Remark: the procedure alters $flexibleOverlap$.
		 */
		public final void copyType_overlap(Edge other) {
			overlap=other.overlap;
		}
		
		
		/**
		 * Assumes that the edge is of overlap type.
		 */
		public final void addType_overlap(Edge other) {
			int myOverlap, otherOverlap;
			
			myOverlap=overlap;
			otherOverlap=other.overlap;
			for (int i=Constants.OVERLAP_PREFIX_PREFIX; i<=Constants.OVERLAP_SUFFIX_SUFFIX; i<<=1) {
				if ((otherOverlap&i)!=0) myOverlap|=i;
			}
			overlap=myOverlap;
		}
		
		
		/**
		 * Removes all information related to an overlap of type $overlapType$.
		 */
		public final void removeOverlap(int overlapType) {
			if (overlap==-1 || (overlap&overlapType)==0) return;
			final int row = Constants.overlap2id(overlapType);
			
			overlap&=~overlapType;
			overhangs[row][0]=-1;
			overhangs[row][1]=-1;
			if (maxOverhangs!=null) {
				maxOverhangs[row][0]=-1;
				maxOverhangs[row][1]=-1;
			}
		}
		

		public final void copyType_insertion(Edge other) {
			insertion=other.insertion;
		}


		public final void copyType_sharedSubstring(Edge other) {
			sharedSubstring=other.sharedSubstring;
		}


		public final boolean differentType_insertion(Edge other) {
			return insertion!=other.insertion;
		}
		
		
		public final int getOverhangID() {
			if (isType(Constants.OVERLAP_PREFIX_PREFIX)) return 0;
			if (isType(Constants.OVERLAP_PREFIX_SUFFIX)) return 1;
			if (isType(Constants.OVERLAP_SUFFIX_PREFIX)) return 2;
			if (isType(Constants.OVERLAP_SUFFIX_SUFFIX)) return 3;
			return -1;
		}
		
		
		/**
		 * @return TRUE iff $edge$ is an overlap that uses the prefix part of $node$.
		 */
		public final boolean isPrefixOverlap(int node) {
			if (overlap==-1 || (node!=nodeID1 && node!=nodeID2)) return false;
			return (overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 ||
				   ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && nodeID1==node) ||
				   ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && nodeID2==node);
		}
	
		
		/**
		 * @return TRUE iff $edge$ is an overlap that uses the suffix part of $node$.
		 */
		public final boolean isSuffixOverlap(int node) {
			if (overlap==-1 || (node!=nodeID1 && node!=nodeID2)) return false;
			return ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && nodeID2==node) ||
				   ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && nodeID1==node) ||
				   (overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0;
		}
		
		
		/**
		 * Remark: the procedure assumes the $components$ field of nodes to be sorted.
		 */
		public boolean isPrintedInComponent(int component) {
			int i;
			Node node;
			
			node=nodesArray[nodeID1];
			i=Arrays.binarySearch(node.components,0,node.lastComponent+1,component);
			if (i<0) {
				System.err.println("Edge.isPrintedInComponent()> ERROR: component "+component+" not found in node "+nodeID1);
				System.exit(1);
			}
			return (printedComponents[i/64]&(1L<<(i%64)))!=0;
		}
		
		
		/**
		 * Remark: the procedure assumes the $components$ field of nodes to be sorted.
		 */
		public void setPrintedInComponent(int component) {
			int i;
			Node node;
			
			node=nodesArray[nodeID1];
			i=Arrays.binarySearch(node.components,0,node.lastComponent+1,component);
			if (i<0) {
				System.err.println("Edge.setPrintedInComponent()> ERROR: component "+component+" not found in node "+nodeID1);
				System.exit(1);
			}
			printedComponents[i/64]|=1L<<(i%64);
		}
		
		
		/**
		 * @return TRUE iff the edge has just one type and one orientation, or just one
		 * type and unknown orientation (but possibly multiple alignments).
		 */
		public final boolean hasOneType() {
			int nTypes;
			
			nTypes=getNOverlaps();
			if (nTypes>1) return false;
			if ((nTypes+(containment==-1?0:1)+(insertion==-1?0:1)+(sharedSubstring==-1?0:1))>1) return false;
			if (orientation==2) return false;
			return true;
		}
		
		
		public final int getNOverlaps() {
			int i, out;
			
			if (overlap==-1) return 0;
			out=0;
			for (i=Constants.OVERLAP_PREFIX_PREFIX; i<=Constants.OVERLAP_SUFFIX_SUFFIX; i<<=1) {
				if ((overlap&i)!=0) out++;
			}
			return out;
		}
		
		
		/**
		 * Resets $otherEdge$ and copies all properties of this edge onto it (excluding 
		 * $id$, which is set to -1).
		 */
		public final void clone(Edge otherEdge) {
			int i, j;
			
			otherEdge.clear();
			otherEdge.nodeID1=nodeID1;
			otherEdge.nodeID2=nodeID2;
			otherEdge.containment=containment;
			otherEdge.containmentSubtype=containmentSubtype;
			otherEdge.overlap=overlap;
			otherEdge.insertion=insertion;
			otherEdge.sharedSubstring=sharedSubstring;
			otherEdge.isRedundant=isRedundant;
			if (overhangs==null) otherEdge.overhangs=null;
			else {
				otherEdge.overhangs = new int[4][2];
				for (i=0; i<overhangs.length; i++) System.arraycopy(overhangs[i],0,otherEdge.overhangs[i],0,overhangs[i].length);
			}
			otherEdge.containmentOverhangLeft=containmentOverhangLeft;
			otherEdge.containmentOverhangRight=containmentOverhangRight;
			if (sharedSubstringOverhangs==null) otherEdge.sharedSubstringOverhangs=null;
			else {
				otherEdge.sharedSubstringOverhangs = new int[4];
				System.arraycopy(sharedSubstringOverhangs,0,otherEdge.sharedSubstringOverhangs,0,sharedSubstringOverhangs.length);
			}
			otherEdge.on=on;
			otherEdge.avgDiffs=avgDiffs;
			otherEdge.nAlignmentsForward=nAlignmentsForward;
			otherEdge.nAlignmentsBackward=nAlignmentsBackward;
			otherEdge.implied=implied;
			otherEdge.orientation=orientation;
		}


		/**
		 * Updates the types of the edge when node $nodeID$ changes just its type to 
		 * $Constants.INTERVAL_ALIGNMENT$.
		 */
		public final void changeType_alignment(int nodeID) {
			if (nodeID==nodeID1) {
				if (containment==Constants.CONTAINMENT_TWO_IN_ONE) {
					if ( (nodesArray[nodeID2].isLeftMaximal || nodesArray[nodeID2].isRightMaximal) &&
						 (nodesArray[nodeID1].type==Constants.INTERVAL_DENSE_PREFIX || nodesArray[nodeID1].type==Constants.INTERVAL_DENSE_SUFFIX || nodesArray[nodeID1].type==Constants.INTERVAL_DENSE_PREFIXSUFFIX || nodesArray[nodeID1].type==Constants.INTERVAL_DENSE_SUBSTRING || nodesArray[nodeID1].type==Constants.INTERVAL_DENSE_SINGLEDELETION) 
					   ) {
						   containment=-1;
						   if (insertion==-1) insertion=Constants.INSERTION_TWO_IN_ONE;
						   else if (insertion==Constants.INSERTION_ONE_IN_TWO) insertion=Constants.INSERTION_BOTH;
					}
				}
				if (overlap!=-1 && (overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && nodesArray[nodeID1].isLeftMaximal) {
					overlap&=~Constants.OVERLAP_PREFIX_SUFFIX;
					flexibleOverlap&=~Constants.OVERLAP_PREFIX_SUFFIX;
				}
				if (overlap!=-1 && (overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 && nodesArray[nodeID1].isLeftMaximal) {
					overlap&=~Constants.OVERLAP_PREFIX_PREFIX;
					flexibleOverlap&=~Constants.OVERLAP_PREFIX_PREFIX;
				}
				if (overlap!=-1 && (overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && nodesArray[nodeID1].isRightMaximal) {
					overlap&=~Constants.OVERLAP_SUFFIX_PREFIX;
					flexibleOverlap&=~Constants.OVERLAP_SUFFIX_PREFIX;
				}
				if (overlap!=-1 && (overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0 && nodesArray[nodeID1].isRightMaximal) {
					overlap&=~Constants.OVERLAP_SUFFIX_SUFFIX;
					flexibleOverlap&=~Constants.OVERLAP_SUFFIX_SUFFIX;
				}
			}
			else {
				if (containment==Constants.CONTAINMENT_ONE_IN_TWO) {
					if ( (nodesArray[nodeID1].isLeftMaximal || nodesArray[nodeID1].isRightMaximal) &&
						 (nodesArray[nodeID2].type==Constants.INTERVAL_DENSE_PREFIX || nodesArray[nodeID2].type==Constants.INTERVAL_DENSE_SUFFIX || nodesArray[nodeID2].type==Constants.INTERVAL_DENSE_PREFIXSUFFIX || nodesArray[nodeID2].type==Constants.INTERVAL_DENSE_SUBSTRING || nodesArray[nodeID2].type==Constants.INTERVAL_DENSE_SINGLEDELETION) 
					   ) {
						   containment=-1;
						   if (insertion==-1) insertion=Constants.INSERTION_ONE_IN_TWO;
						   else if (insertion==Constants.INSERTION_TWO_IN_ONE) insertion=Constants.INSERTION_BOTH;
					}
				}
				if (overlap!=-1 && (overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && nodesArray[nodeID2].isRightMaximal) {
					overlap&=~Constants.OVERLAP_PREFIX_SUFFIX;
					flexibleOverlap&=~Constants.OVERLAP_PREFIX_SUFFIX;
				}
				if (overlap!=-1 && (overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 && nodesArray[nodeID2].isLeftMaximal) {
					overlap&=~Constants.OVERLAP_PREFIX_PREFIX;
					flexibleOverlap&=~Constants.OVERLAP_PREFIX_PREFIX;
				}
				if (overlap!=-1 && (overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && nodesArray[nodeID2].isLeftMaximal) {
					overlap&=~Constants.OVERLAP_SUFFIX_PREFIX;
					flexibleOverlap&=~Constants.OVERLAP_SUFFIX_PREFIX;
				}
				if (overlap!=-1 && (overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0 && nodesArray[nodeID2].isRightMaximal) {
					overlap&=~Constants.OVERLAP_SUFFIX_SUFFIX;
					flexibleOverlap&=~Constants.OVERLAP_SUFFIX_SUFFIX;
				}
			}
		}
		
		
		/**
		 * Updates the types of the edge after the type of a node has been changed to 
		 * $Constants.INTERVAL_PERIODIC$.
		 *
		 * @param tmpEdge1 temporary space;
		 * @param out the new edge after changing node type.
		 */
		public final void changeType_periodic(Edge tmpEdge1, Edge out) {
			int k;
			int row, result;
			
			out.clear();
			tmpEdge1.clear();
			tmpEdge1.nodeID1=nodeID1;
			tmpEdge1.nodeID2=nodeID2;
			out.nodeID1=nodeID1;
			out.nodeID2=nodeID2;
			
			if (containment==Constants.CONTAINMENT_IDENTICAL) {
				tmpEdge1.containment=Constants.CONTAINMENT_IDENTICAL;
				tmpEdge1.containmentOverhangLeft=0;
				tmpEdge1.containmentOverhangRight=0;
				tmpEdge1.orientation=orientation;
				result=IntervalGraphStep2.inStep2_checkEdge(tmpEdge1);
				if (result!=-1) {
					out.addTypes(tmpEdge1);
					if (orientation==0 || orientation==1) out.addOrientation(tmpEdge1);
					if (orientation==1) out.nAlignmentsBackward++;
					else out.nAlignmentsForward++;
				}
				tmpEdge1.setType_noContainment();
			}
			else if (containment==Constants.CONTAINMENT_TWO_IN_ONE || containment==Constants.CONTAINMENT_ONE_IN_TWO) {
				tmpEdge1.containment=containment;
				tmpEdge1.containmentOverhangLeft=containmentOverhangLeft;
				tmpEdge1.containmentOverhangRight=containmentOverhangRight;
				tmpEdge1.orientation=orientation;
				result=IntervalGraphStep2.inStep2_checkEdge(tmpEdge1);
				if (result!=-1) {
					out.addTypes(tmpEdge1);
					if (orientation==0 || orientation==1) out.addOrientation(tmpEdge1);
					if (orientation==1) out.nAlignmentsBackward++;
					else out.nAlignmentsForward++;
				}
				tmpEdge1.setType_noContainment();
			}
			if (overlap!=-1) {
				for (k=Constants.OVERLAP_PREFIX_PREFIX; k<=Constants.OVERLAP_SUFFIX_SUFFIX; k<<=1) {
					if ((overlap&k)!=0) {
						tmpEdge1.overlap=k;
						row=Constants.overlap2id(k);
						if (tmpEdge1.overhangs==null) {
							tmpEdge1.overhangs = new int[4][2];
							Math.set(tmpEdge1.overhangs,-1);
						}
						tmpEdge1.overhangs[row][0]=overhangs[row][0];
						tmpEdge1.overhangs[row][1]=overhangs[row][1];
						if (maxOverhangs!=null) {
							if (tmpEdge1.maxOverhangs==null) {
								tmpEdge1.maxOverhangs = new int[4][2];
								Math.set(tmpEdge1.maxOverhangs,-1);
							}
							tmpEdge1.maxOverhangs[row][0]=maxOverhangs[row][0];
							tmpEdge1.maxOverhangs[row][1]=maxOverhangs[row][1];
						}
						tmpEdge1.orientation=Constants.overlap2orientation(k);
						result=IntervalGraphStep2.inStep2_checkEdge(tmpEdge1);
						if (result!=-1) {
							out.addTypes(tmpEdge1);
							out.addOrientation(tmpEdge1);
							if (tmpEdge1.orientation==1) out.nAlignmentsBackward++;
							else out.nAlignmentsForward++;
						}
						tmpEdge1.setType_noOverlap();
					}
				}
			}
			if (insertion!=-1) {
				tmpEdge1.insertion=insertion;
				tmpEdge1.containmentOverhangLeft=containmentOverhangLeft;
				tmpEdge1.containmentOverhangRight=containmentOverhangRight;
				tmpEdge1.orientation=orientation;
				result=IntervalGraphStep2.inStep2_checkEdge(tmpEdge1);
				if (result!=-1) {
					out.addTypes(tmpEdge1);
					if (orientation==0 || orientation==1) out.addOrientation(tmpEdge1);
					if (orientation==1) out.nAlignmentsBackward++;
					else out.nAlignmentsForward++;
				}
				tmpEdge1.setType_noInsertion();
			}
			if (sharedSubstring!=-1) {
				tmpEdge1.sharedSubstring=sharedSubstring;
				if (sharedSubstringOverhangs!=null) {
					if (tmpEdge1.sharedSubstringOverhangs==null) tmpEdge1.sharedSubstringOverhangs = new int[4];
					System.arraycopy(sharedSubstringOverhangs,0,tmpEdge1.sharedSubstringOverhangs,0,4);
				}
				tmpEdge1.orientation=orientation;
				result=IntervalGraphStep2.inStep2_checkEdge(tmpEdge1);
				if (result!=-1) {
					out.addTypes(tmpEdge1);
					if (orientation==0 || orientation==1) out.addOrientation(tmpEdge1);
					if (orientation==1) out.nAlignmentsBackward++;
					else out.nAlignmentsForward++;
				}
				tmpEdge1.setType_noSharedSubstring();
			}
			
			out.avgDiffs=(avgDiffs/getNAlignments())*out.getNAlignments();
		}
		
		
		/**
		 * @return the number of distinct overlap types in this edge.
		 */
		public final int nOverlapTypes() {
			int i, nTypes;
		
			if (overlap==-1) return 0;
			nTypes=0;
			for (i=Constants.OVERLAP_PREFIX_PREFIX; i<=Constants.OVERLAP_SUFFIX_SUFFIX; i<<=1) {
				if ((overlap&i)!=0) nTypes++;
			}
			return nTypes;
		}
		
		
		/**
		 * @return the number of overlap alignments whose coordinates are stored in the 
		 * edge.
		 */
		public final int nStoredAlignments_overlap() {
			int i, j, out;
		
			if (overlap==-1) return 0;
			out=0;
			for (i=Constants.OVERLAP_PREFIX_PREFIX; i<=Constants.OVERLAP_SUFFIX_SUFFIX; i<<=1) {
				if ((overlap&i)==0) continue;
				out++;
				j=Constants.overlap2id(i);
				if ( maxOverhangs!=null && 
					 ( (overhangs[j][0]!=-1 && maxOverhangs[j][0]!=-1 && maxOverhangs[j][0]!=overhangs[j][0]) || 
					   (overhangs[j][1]!=-1 && maxOverhangs[j][1]!=-1 && maxOverhangs[j][1]!=overhangs[j][1])
					 )
				   ) out++;
			}
			return out;
		}
		
		
		public final int nInsertionTypes() {
			if (insertion==-1) return 0;
			return insertion==Constants.INSERTION_BOTH?2:1;
		}
	
	
		/**
		 * Remark: the procedure assumes the length of the output overhang to be nonzero, 
		 * if it exists.
		 *
		 * @param prefixOrSuffix TRUE=prefix, FALSE=suffix;
		 * @return the length of the overlap overhang, in $to$, between node $from$ and 
		 * $to$, where $to$ is the other node of the edge that uses side $prefixOrSuffix$ 
		 * of $from$. The returned value is positive (respectively, negative) if the 
		 * overhang substring is a prefix (respectively, a suffix) of $to$. If the edge is
		 * not an overlap, or if it does not use side $prefixOrSuffix$ of $from$, or if 
		 * $to$ overlaps side $prefixOrSuffix$ of $from$ on two sides of $to$, or if $to$ 
		 * overlaps side $prefixOrSuffix$ of $from$ multiple times on the same side of 
		 * $to$, the procedure returns $Math.POSITIVE_INFINITY$.
		 */
		public final int getOverhang(int from, boolean prefixOrSuffix) {
			int out1, out2;
			
			if (overlap==-1) return Math.POSITIVE_INFINITY;
			out1=Math.POSITIVE_INFINITY; out2=Math.POSITIVE_INFINITY;
			if (prefixOrSuffix) {
				if (from==nodeID1) {
					if ((overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 && (maxOverhangs==null || (maxOverhangs[0][0]==overhangs[0][0] && maxOverhangs[0][1]==overhangs[0][1]))) out1=-overhangs[0][1];
					if ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && (maxOverhangs==null || (maxOverhangs[1][0]==overhangs[1][0] && maxOverhangs[1][1]==overhangs[1][1]))) out2=overhangs[1][1];
				}
				else {
					if ((overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0 && (maxOverhangs==null || (maxOverhangs[0][0]==overhangs[0][0] && maxOverhangs[0][1]==overhangs[0][1]))) out1=-overhangs[0][0];
					if ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && (maxOverhangs==null || (maxOverhangs[2][0]==overhangs[2][0] && maxOverhangs[2][1]==overhangs[2][1]))) out2=overhangs[2][0];
				}
			}
			else {
				if (from==nodeID1) {
					if ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0 && (maxOverhangs==null || (maxOverhangs[2][0]==overhangs[2][0] && maxOverhangs[2][1]==overhangs[2][1]))) out1=-overhangs[2][1];
					if ((overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0 && (maxOverhangs==null || (maxOverhangs[3][0]==overhangs[3][0] && maxOverhangs[3][1]==overhangs[3][1]))) out2=overhangs[3][1];
				}
				else {
					if ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0 && (maxOverhangs==null || (maxOverhangs[1][0]==overhangs[1][0] && maxOverhangs[1][1]==overhangs[1][1]))) out1=-overhangs[1][0];
					if ((overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0 && (maxOverhangs==null || (maxOverhangs[3][0]==overhangs[3][0] && maxOverhangs[3][1]==overhangs[3][1]))) out2=overhangs[3][0];
				}
			}
			if (out1!=Math.POSITIVE_INFINITY) {
				if (out2!=Math.POSITIVE_INFINITY) return Math.POSITIVE_INFINITY;
				return out1;
			}
			else return out2;
		}
		
		
		/**
		 * Same as $getOverhang()$, but the edge is not required to be of overlap type, 
		 * i.e. the overlap we are testing for does not need to satisfy type and
		 * maximality constraints.
		 *
		 * @param tmpArray* temporary space, of size at least equal to twice the number of
		 * alignments in the edge.
		 */
		public final int getOverhang_weak(int from, boolean prefixOrSuffix, int threshold, int[] tmpArray1, int[] tmpArray2) {
			int k;
			int to, last1, last2, overlapType, length;
			int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
			Node node1, node2;

			last1=projectOnto(from,tmpArray1,false);
			if (last1==-1) return Math.POSITIVE_INFINITY;
			to=getTo(from); 
			last2=projectOnto(to,tmpArray2,false);
			if (last2==-1 || last2!=last1) return Math.POSITIVE_INFINITY;
			node1=nodesArray[from]; start1=node1.start; end1=node1.end;
			node2=nodesArray[to]; start2=node2.start; end2=node2.end;
			length=Math.POSITIVE_INFINITY;
			for (k=0; k<last1; k+=2) {
				// This works since $projectOnto()$ writes its output in the same order
				// for both nodes.
				alignmentStart1=start1+tmpArray1[k];
				alignmentEnd1=start1+tmpArray1[k+1];
				alignmentStart2=start2+tmpArray2[k];
				alignmentEnd2=start2+tmpArray2[k+1];
				overlapType=IntervalGraphStep2.isOverlap_simple_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,orientation==0,threshold);
				if (overlapType==-1) continue;
				if (overlapType==Constants.OVERLAP_PREFIX_PREFIX && prefixOrSuffix) {
					if (length!=Math.POSITIVE_INFINITY) return Math.POSITIVE_INFINITY;
					length=-(end2-alignmentEnd2);
				}
				else if (overlapType==Constants.OVERLAP_PREFIX_SUFFIX && prefixOrSuffix) {
					if (length!=Math.POSITIVE_INFINITY) return Math.POSITIVE_INFINITY;
					length=alignmentStart2-start2;
				}
				else if (overlapType==Constants.OVERLAP_SUFFIX_PREFIX && !prefixOrSuffix) {
					if (length!=Math.POSITIVE_INFINITY) return Math.POSITIVE_INFINITY;
					length=-(end2-alignmentEnd2);
				}
				else if (overlapType==Constants.OVERLAP_SUFFIX_SUFFIX && !prefixOrSuffix) {
					if (length!=Math.POSITIVE_INFINITY) return Math.POSITIVE_INFINITY;
					length=alignmentStart2-start2;
				}
			}			
			return length;
		}
		
		
		/**
		 * Like $getOverhang_weak()$, but returns true iff at least one of the non-overlap
		 * alignments encoded in the edge is an overlap (except for type and maximality
		 * constraints).
		 */
		public final boolean isOverlap_weak(int threshold, int[] tmpArray1, int[] tmpArray2) {
			int k;
			int last1, last2, overlapType;
			int start1, end1, start2, end2, alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
			Node node1, node2;

			if (orientation==2 || orientation==-1) return false;
			last1=projectOnto(nodeID1,tmpArray1,false);
			if (last1==-1) return false;
			last2=projectOnto(nodeID2,tmpArray2,false);
			if (last2==-1 || last2!=last1) return false;
			node1=nodesArray[nodeID1]; start1=node1.start; end1=node1.end;
			node2=nodesArray[nodeID2]; start2=node2.start; end2=node2.end;
			for (k=0; k<last1; k+=2) {
				// This works since $projectOnto()$ writes its output in the same order
				// for both nodes.
				alignmentStart1=start1+tmpArray1[k];
				alignmentEnd1=start1+tmpArray1[k+1];
				alignmentStart2=start2+tmpArray2[k];
				alignmentEnd2=start2+tmpArray2[k+1];
				overlapType=IntervalGraphStep2.isOverlap_simple_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,orientation==0,threshold);
				if (overlapType!=-1) return true;
			}			
			return false;
		}
		
		
		/**
		 * Checks whether the current edge contains an alignment that is identical to the
		 * one stored in $otherEdge$ (which is assumed to contain just one alignment),
		 * i.e. an alignment with the same type and with the same start and end positions 
		 * in both intervals.
		 *
		 * Remark: due to the way edges are represented, the alignment in $otherEdge$ 
		 * might occur inside the current edge, but might be invisible to this procedure.
		 *
		 * @param threshold (>=0) for deciding identity.
		 */
		public final boolean contains(Edge otherEdge, int threshold) {
			final boolean direction;
			int i, j;
			int[] tmpArray = new int[4];
			
			if (nodeID1==otherEdge.nodeID1 && nodeID2==otherEdge.nodeID2) direction=true;
			else if (nodeID1==otherEdge.nodeID2 && nodeID2==otherEdge.nodeID1) direction=false;
			else return false;
			
			if (otherEdge.containment!=-1 && containment!=-1) {
				if (otherEdge.containment==Constants.CONTAINMENT_IDENTICAL) return containment==Constants.CONTAINMENT_IDENTICAL;
				if (otherEdge.containment!=(direction?containment:Constants.oppositeDirectionContainment(containment))) return false;
				projectOnto_containment(nodeID1,tmpArray,-1);
				otherEdge.projectOnto_containment(nodeID1,tmpArray,1);
				if (Math.abs(tmpArray[0],tmpArray[2])>threshold || Math.abs(tmpArray[1],tmpArray[3])>threshold) return false;
				projectOnto_containment(nodeID2,tmpArray,-1);
				otherEdge.projectOnto_containment(nodeID2,tmpArray,1);
				return Math.abs(tmpArray[0],tmpArray[2])<=threshold && Math.abs(tmpArray[1],tmpArray[3])<=threshold;
			}
			else if (otherEdge.overlap!=-1 && overlap!=-1 && (overlap&(direction?otherEdge.overlap:Constants.oppositeDirectionOverlap(otherEdge.overlap)))!=0) {
				i=Constants.overlap2id(direction?otherEdge.overlap:Constants.oppositeDirectionOverlap(otherEdge.overlap));
				j=Constants.overlap2id(otherEdge.overlap);
				return (Math.abs(overhangs[i][0],otherEdge.overhangs[j][direction?0:1])<=threshold && Math.abs(overhangs[i][1],otherEdge.overhangs[j][direction?1:0])<=threshold) ||
					   (maxOverhangs!=null && Math.abs(maxOverhangs[i][0],otherEdge.overhangs[j][direction?0:1])<=threshold && Math.abs(maxOverhangs[i][1],otherEdge.overhangs[j][direction?1:0])<=threshold);
			}
			else if (otherEdge.insertion!=-1 && insertion!=-1 && insertion==(direction?otherEdge.insertion:Constants.oppositeDirectionInsertion(otherEdge.insertion))) {
				projectOnto_insertion(nodeID1,tmpArray,-1);
				otherEdge.projectOnto_insertion(nodeID1,tmpArray,1);
				if (Math.abs(tmpArray[0],tmpArray[2])>threshold || Math.abs(tmpArray[1],tmpArray[3])>threshold) return false;
				projectOnto_insertion(nodeID2,tmpArray,-1);
				otherEdge.projectOnto_insertion(nodeID2,tmpArray,1);
				return Math.abs(tmpArray[0],tmpArray[2])<=threshold && Math.abs(tmpArray[1],tmpArray[3])<=threshold;
			}
			else if (otherEdge.sharedSubstring!=-1 && sharedSubstring!=-1) {
				projectOnto_sharedSubstring(nodeID1,tmpArray,-1);
				otherEdge.projectOnto_sharedSubstring(nodeID1,tmpArray,1);
				if (Math.abs(tmpArray[0],tmpArray[2])>threshold || Math.abs(tmpArray[1],tmpArray[3])>threshold) return false;
				projectOnto_sharedSubstring(nodeID2,tmpArray,-1);
				otherEdge.projectOnto_sharedSubstring(nodeID2,tmpArray,1);
				return Math.abs(tmpArray[0],tmpArray[2])<=threshold && Math.abs(tmpArray[1],tmpArray[3])<=threshold;
			}
			return false;
		}
		
		
		/**
		 * @return TRUE iff $maxOverhangs$ is not null and is different from $overhangs$ 
		 * in at least one cell.
		 */
		public final boolean overhangMatricesDiffer() {
			if (overhangs==null || maxOverhangs==null) return false;
			for (int i=0; i<4; i++) {
				if (overhangs[i][0]!=maxOverhangs[i][0] || overhangs[i][1]!=maxOverhangs[i][1]) return true;
			}
			return false;
		}
		
		
		/**
		 * Draws a schematic representation of the edge.
		 *
		 * @param reference ID of the node to use as reference in the drawing;
		 * @param fromX,fromY coordinates of the pixel from which to draw $reference$; the
		 * drawing can extend to values smaller than $fromX$, but not smaller than 
		 * $fromY$;
		 * @param referencePointX,referencePointY coordinates of the top-left pixel of
		 * $reference$, drawn inside its snippet in $image$;
		 * @param otherPointX,otherPointY coordinates of the top-left pixel of the other
		 * interval, drawn inside its snippet in $image$.
		 */
		public final void draw(BufferedImage image, int reference, int fromX, int fromY, int referencePointX, int referencePointY, int otherPointX, int otherPointY, int referenceSnippetMinY, int referenceSnippetMaxY, int otherSnippetMinY, int otherSnippetMaxY) {
			final int QUANTUM = 10;
			final int INTERVAL_HEIGHT = 20;
			final int VERTICAL_SPACE = 40;
			final Color LINE_COLOR = new Color(Colors.COLOR_GUIDE);
			final int LINE_WIDTH = 4;
			final int LINE_WIDTH_SNIPPET_GUIDES = 2;
			final int TEXT_OFFSET = 4;
			final int COLOR_TEXT = Colors.COLOR_TEXT;
			final int COLOR_SNIPPET_GUIDES = Colors.COLOR_HISTOGRAM;
			final Node referenceNode = nodesArray[reference];
			final int referenceColor = Colors.type2color(referenceNode.type);
			final int referenceLength = referenceNode.end-referenceNode.start+1;
			final int otherNode = nodeID1==reference?nodeID2:nodeID1;
			final int otherColor = Colors.type2color(nodesArray[otherNode].type);
			final int otherLength = nodesArray[otherNode].end-nodesArray[otherNode].start+1;
			int x, y, p1, q1, p2, q2;
			int from, to, lastY, delta, overhangLeft, overhangRight;
			String label;
			Graphics2D graphics = image.createGraphics();
			graphics.setFont(new Font("SansSerif ",Font.PLAIN,15));
			FontMetrics fontMetrics = graphics.getFontMetrics();
			
			// Drawing reference
			for (x=0; x<=referenceLength/QUANTUM; x++) {
				for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(fromX+x,fromY+y,referenceColor);
			}
			Colors.drawTriangle(image,fromX+referenceLength/QUANTUM,fromY,INTERVAL_HEIGHT,true,referenceColor);
			label=(nAlignmentsForward>0?" nFWD="+nAlignmentsForward:"")+(nAlignmentsBackward>0?" nBKW="+nAlignmentsBackward:"")+" ORIENT="+orientation+" diffs="+IO.format(avgDiffs/getNAlignments())+(on?"":" OFF")+(supplement?" SUPP":"");
			graphics.setColor(new Color(COLOR_TEXT));
			graphics.drawString(label,fromX+TEXT_OFFSET,fromY+INTERVAL_HEIGHT-TEXT_OFFSET);
			lastY=fromY+INTERVAL_HEIGHT+VERTICAL_SPACE;
		
			// Overlap
			if (overlap!=-1) {
				if ((overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) {
					p1=fromX; 
					q1=fromX+(referenceLength-overhangs[0][reference==nodeID1?0:1])/QUANTUM;
					q2=q1;
					from=q2-otherLength/QUANTUM;
					p2=from+overhangs[0][otherNode==nodeID1?0:1]/QUANTUM;
					for (x=from; x<=q2; x++) {
						for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(x,lastY+y,otherColor);
					}
					Colors.drawTriangle(image,from,lastY,INTERVAL_HEIGHT,false,otherColor);
					label="OVERLAP PP";
					label+=" // "+overhangs[0][otherNode==nodeID1?0:1];
					label+=" // "+IO.format(((double)Math.max(q1-p1+1,q2-p2+1))/Math.min(q1-p1+1,q2-p2+1))+"x";
					graphics.setColor(new Color(COLOR_TEXT));
					graphics.drawString(label,from+TEXT_OFFSET,lastY+INTERVAL_HEIGHT-TEXT_OFFSET);
					graphics.setColor(LINE_COLOR); graphics.setStroke(new BasicStroke(LINE_WIDTH));
					graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,p2,lastY);
					graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,q2,lastY);
					lastY+=INTERVAL_HEIGHT+VERTICAL_SPACE;
					// Drawing inside the snippets
					graphics.setColor(new Color(COLOR_SNIPPET_GUIDES));
					graphics.setStroke(new BasicStroke(LINE_WIDTH_SNIPPET_GUIDES));
					graphics.drawLine(referencePointX+p1-fromX,referenceSnippetMinY,referencePointX+p1-fromX,referenceSnippetMaxY);
					graphics.drawLine(referencePointX+q1-fromX,referenceSnippetMinY,referencePointX+q1-fromX,referenceSnippetMaxY);
					graphics.drawLine(otherPointX,otherSnippetMinY,otherPointX,otherSnippetMaxY);
					graphics.drawLine(otherPointX+q2-p2,otherSnippetMinY,otherPointX+q2-p2,otherSnippetMaxY);
				}
				if ((overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) {
					if (reference==nodeID1) {
						p1=fromX; 
						q1=fromX+(referenceLength-overhangs[1][0])/QUANTUM;
						q2=q1;
						from=q2-otherLength/QUANTUM;
						p2=from+overhangs[1][1]/QUANTUM;
						to=from+otherLength/QUANTUM;
					}
					else {
						p1=fromX+overhangs[1][1]/QUANTUM;
						q1=fromX+referenceLength/QUANTUM;
						p2=p1;
						from=p2;
						q2=from+(otherLength-overhangs[1][0])/QUANTUM;
						to=from+otherLength/QUANTUM;
					}
					for (x=from; x<=to; x++) {
						for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(x,lastY+y,otherColor);
					}
					Colors.drawTriangle(image,to,lastY,INTERVAL_HEIGHT,true,otherColor);
					label="OVERLAP "+(reference==nodeID1?"PS":"SP");
					label+=" // "+overhangs[1][otherNode==nodeID1?0:1];
					label+=" // "+IO.format(((double)Math.max(q1-p1+1,q2-p2+1))/Math.min(q1-p1+1,q2-p2+1))+"x";
					graphics.setColor(new Color(COLOR_TEXT));
					graphics.drawString(label,from+TEXT_OFFSET,lastY+INTERVAL_HEIGHT-TEXT_OFFSET);
					graphics.setColor(LINE_COLOR); graphics.setStroke(new BasicStroke(LINE_WIDTH));
					graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,p2,lastY);
					graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,q2,lastY);
					lastY+=INTERVAL_HEIGHT+VERTICAL_SPACE;
					// Drawing inside the snippet
					graphics.setColor(new Color(COLOR_SNIPPET_GUIDES));
					graphics.setStroke(new BasicStroke(LINE_WIDTH_SNIPPET_GUIDES));
					graphics.drawLine(referencePointX+p1-fromX,referenceSnippetMinY,referencePointX+p1-fromX,referenceSnippetMaxY);
					graphics.drawLine(referencePointX+q1-fromX,referenceSnippetMinY,referencePointX+q1-fromX,referenceSnippetMaxY);
					graphics.drawLine(otherPointX+p2-from,otherSnippetMinY,otherPointX+p2-from,otherSnippetMaxY);
					graphics.drawLine(otherPointX+q2-from,otherSnippetMinY,otherPointX+q2-from,otherSnippetMaxY);
				}
				if ((overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) {
					if (reference==nodeID1) {
						p1=fromX+overhangs[2][0]/QUANTUM;
						q1=fromX+referenceLength/QUANTUM;
						p2=p1;
						from=p2;
						to=p2+otherLength/QUANTUM;
						q2=to-overhangs[2][1]/QUANTUM;
					}
					else {
						p1=fromX;
						q1=fromX+(referenceLength-overhangs[2][1])/QUANTUM;
						q2=q1;
						from=q2-otherLength/QUANTUM;
						p2=from+overhangs[2][0]/QUANTUM;
						to=from+otherLength/QUANTUM;
					}
					for (x=from; x<=to; x++) {
						for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(x,lastY+y,otherColor);
					}
					Colors.drawTriangle(image,to,lastY,INTERVAL_HEIGHT,true,otherColor);
					label="OVERLAP "+(reference==nodeID1?"SP":"PS");
					label+=" // "+overhangs[2][otherNode==nodeID1?0:1];
					label+=" // "+IO.format(((double)Math.max(q1-p1+1,q2-p2+1))/Math.min(q1-p1+1,q2-p2+1))+"x";
					graphics.setColor(new Color(COLOR_TEXT));
					graphics.drawString(label,from+TEXT_OFFSET,lastY+INTERVAL_HEIGHT-TEXT_OFFSET);
					graphics.setColor(LINE_COLOR); graphics.setStroke(new BasicStroke(LINE_WIDTH));
					graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,p2,lastY);
					graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,q2,lastY);
					lastY+=INTERVAL_HEIGHT+VERTICAL_SPACE;
					// Drawing inside the snippet
					graphics.setColor(new Color(COLOR_SNIPPET_GUIDES));
					graphics.setStroke(new BasicStroke(LINE_WIDTH_SNIPPET_GUIDES));
					graphics.drawLine(referencePointX+p1-fromX,referenceSnippetMinY,referencePointX+p1-fromX,referenceSnippetMaxY);
					graphics.drawLine(referencePointX+q1-fromX,referenceSnippetMinY,referencePointX+q1-fromX,referenceSnippetMaxY);
					graphics.drawLine(otherPointX+p2-from,otherSnippetMinY,otherPointX+p2-from,otherSnippetMaxY);
					graphics.drawLine(otherPointX+q2-from,otherSnippetMinY,otherPointX+q2-from,otherSnippetMaxY);
				}
				if ((overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) {
					p1=fromX+overhangs[3][reference==nodeID1?0:1]/QUANTUM;
					q1=fromX+referenceLength/QUANTUM;
					p2=p1;
					to=p2+otherLength/QUANTUM;
					q2=to-overhangs[3][otherNode==nodeID1?0:1]/QUANTUM;
					for (x=p2; x<=to; x++) {
						for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(x,lastY+y,otherColor);
					}
					Colors.drawTriangle(image,p2,lastY,INTERVAL_HEIGHT,false,otherColor);
					label="OVERLAP SS";
					label+=" // "+overhangs[3][otherNode==nodeID1?0:1];
					label+=" // "+IO.format(((double)Math.max(q1-p1+1,q2-p2+1))/Math.min(q1-p1+1,q2-p2+1))+"x";
					graphics.setColor(new Color(COLOR_TEXT));
					graphics.drawString(label,p2+TEXT_OFFSET,lastY+INTERVAL_HEIGHT-TEXT_OFFSET);
					graphics.setColor(LINE_COLOR); graphics.setStroke(new BasicStroke(LINE_WIDTH));
					graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,p2,lastY);
					graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,q2,lastY);
					lastY+=INTERVAL_HEIGHT+VERTICAL_SPACE;
					// Drawing inside the snippet
					graphics.setColor(new Color(COLOR_SNIPPET_GUIDES));
					graphics.setStroke(new BasicStroke(LINE_WIDTH_SNIPPET_GUIDES));
					graphics.drawLine(referencePointX+p1-fromX,referenceSnippetMinY,referencePointX+p1-fromX,referenceSnippetMaxY);
					graphics.drawLine(referencePointX+q1-fromX,referenceSnippetMinY,referencePointX+q1-fromX,referenceSnippetMaxY);
					graphics.drawLine(otherPointX+to-q2,otherSnippetMinY,otherPointX+to-q2,otherSnippetMaxY);
					graphics.drawLine(otherPointX+to-p2,otherSnippetMinY,otherPointX+to-p2,otherSnippetMaxY);
				}
			}
		
			// Containment/insertion
			if (containment==Constants.CONTAINMENT_IDENTICAL) {
				p1=fromX;
				q1=fromX+referenceLength/QUANTUM;
				if (referenceLength>otherLength) {
					delta=(referenceLength-otherLength)/2;
					p2=p1+delta/QUANTUM;
				}
				else {
					delta=(otherLength-referenceLength)/2;
					p2=p1-delta/QUANTUM;
				}				
				q2=p2+otherLength/QUANTUM;
				for (x=p2; x<=q2; x++) {
					for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(x,lastY+y,otherColor);
				}
				if (orientation==0) Colors.drawTriangle(image,q2,lastY,INTERVAL_HEIGHT,true,otherColor);
				else if (orientation==1) Colors.drawTriangle(image,p2,lastY,INTERVAL_HEIGHT,false,otherColor);
				label="IDENTITY";
				label+=" // "+IO.format(((double)Math.max(q1-p1+1,q2-p2+1))/Math.min(q1-p1+1,q2-p2+1))+"x";
				graphics.setColor(new Color(COLOR_TEXT));
				graphics.drawString(label,p2+TEXT_OFFSET,lastY+INTERVAL_HEIGHT-TEXT_OFFSET);
				graphics.setColor(LINE_COLOR); graphics.setStroke(new BasicStroke(LINE_WIDTH));
				graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,p2,lastY);
				graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,q2,lastY);
				lastY+=INTERVAL_HEIGHT+VERTICAL_SPACE;
				// Drawing inside the snippet
				graphics.setColor(new Color(COLOR_SNIPPET_GUIDES));
				graphics.setStroke(new BasicStroke(LINE_WIDTH_SNIPPET_GUIDES));
				graphics.drawLine(referencePointX+p1-fromX,referenceSnippetMinY,referencePointX+p1-fromX,referenceSnippetMaxY);
				graphics.drawLine(referencePointX+q1-fromX,referenceSnippetMinY,referencePointX+q1-fromX,referenceSnippetMaxY);
				graphics.drawLine(otherPointX,otherSnippetMinY,otherPointX,otherSnippetMaxY);
				graphics.drawLine(otherPointX+q2-p2,otherSnippetMinY,otherPointX+q2-p2,otherSnippetMaxY);
			}
			else if ( ((containment==Constants.CONTAINMENT_ONE_IN_TWO||insertion==Constants.INSERTION_ONE_IN_TWO) && reference==nodeID2) || 
			          ((containment==Constants.CONTAINMENT_TWO_IN_ONE||insertion==Constants.INSERTION_TWO_IN_ONE) && reference==nodeID1)
			        ) {
				p1=fromX+containmentOverhangLeft/QUANTUM;
				q1=fromX+(referenceLength-containmentOverhangRight)/QUANTUM;
				p2=p1;
				q2=p1+otherLength/QUANTUM;
				for (x=p2; x<=q2; x++) {
					for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(x,lastY+y,otherColor);
				}
				if (orientation==0) Colors.drawTriangle(image,q2,lastY,INTERVAL_HEIGHT,true,otherColor);
				else if (orientation==1) Colors.drawTriangle(image,p2,lastY,INTERVAL_HEIGHT,false,otherColor);
				if (containment==Constants.CONTAINMENT_ONE_IN_TWO || containment==Constants.CONTAINMENT_TWO_IN_ONE) label="CONTAINMENT";
				else label="INSERTION";
				label+=" // "+IO.format(((double)Math.max(q1-p1+1,q2-p2+1))/Math.min(q1-p1+1,q2-p2+1))+"x";
				graphics.setColor(new Color(COLOR_TEXT));
				graphics.drawString(label,p2+TEXT_OFFSET,lastY+INTERVAL_HEIGHT-TEXT_OFFSET);
				graphics.setColor(LINE_COLOR); graphics.setStroke(new BasicStroke(LINE_WIDTH));
				graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,p2,lastY);
				graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,q2,lastY);
				lastY+=INTERVAL_HEIGHT+VERTICAL_SPACE;
				// Drawing inside the snippet
				graphics.setColor(new Color(COLOR_SNIPPET_GUIDES));
				graphics.setStroke(new BasicStroke(LINE_WIDTH_SNIPPET_GUIDES));
				graphics.drawLine(referencePointX+p1-fromX,referenceSnippetMinY,referencePointX+p1-fromX,referenceSnippetMaxY);
				graphics.drawLine(referencePointX+q1-fromX,referenceSnippetMinY,referencePointX+q1-fromX,referenceSnippetMaxY);
				graphics.drawLine(otherPointX,otherSnippetMinY,otherPointX,otherSnippetMaxY);
				graphics.drawLine(otherPointX+q2-p2,otherSnippetMinY,otherPointX+q2-p2,otherSnippetMaxY);
			}
			else if ( ((containment==Constants.CONTAINMENT_ONE_IN_TWO||insertion==Constants.INSERTION_ONE_IN_TWO) && reference==nodeID1) || 
			          ((containment==Constants.CONTAINMENT_TWO_IN_ONE||insertion==Constants.INSERTION_TWO_IN_ONE) && reference==nodeID2)
			        ) {
				p1=fromX;
				q1=p1+referenceLength/QUANTUM;
				p2=p1;
				from=p2-containmentOverhangLeft/QUANTUM;
				q2=from+(otherLength-containmentOverhangRight)/QUANTUM;
				to=from+otherLength/QUANTUM;
				for (x=from; x<=to; x++) {
					for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(x,lastY+y,otherColor);
				}
				Colors.drawTriangle(image,to,lastY,INTERVAL_HEIGHT,true,otherColor);
				if (containment==Constants.CONTAINMENT_ONE_IN_TWO || containment==Constants.CONTAINMENT_TWO_IN_ONE) label="CONTAINMENT";
				else label="INSERTION";
				label+=" // "+containmentOverhangLeft+","+containmentOverhangRight;
				label+=" // "+IO.format(((double)Math.max(q1-p1+1,q2-p2+1))/Math.min(q1-p1+1,q2-p2+1))+"x";
				graphics.setColor(new Color(COLOR_TEXT));
				graphics.drawString(label,from+TEXT_OFFSET,lastY+INTERVAL_HEIGHT-TEXT_OFFSET);
				graphics.setColor(LINE_COLOR); graphics.setStroke(new BasicStroke(LINE_WIDTH));
				if (orientation==0) {
					graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,p2,lastY);
					graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,q2,lastY);
				}
				else {
					graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,q2,lastY);
					graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,p2,lastY);
				}
				lastY+=INTERVAL_HEIGHT+VERTICAL_SPACE;
				// Drawing inside the snippet
				graphics.setColor(new Color(COLOR_SNIPPET_GUIDES));
				graphics.setStroke(new BasicStroke(LINE_WIDTH_SNIPPET_GUIDES));
				graphics.drawLine(referencePointX+p1-fromX,referenceSnippetMinY,referencePointX+p1-fromX,referenceSnippetMaxY);
				graphics.drawLine(referencePointX+q1-fromX,referenceSnippetMinY,referencePointX+q1-fromX,referenceSnippetMaxY);
				graphics.drawLine(otherPointX+p2-from,otherSnippetMinY,otherPointX+p2-from,otherSnippetMaxY);
				graphics.drawLine(otherPointX+q2-from,otherSnippetMinY,otherPointX+q2-from,otherSnippetMaxY);
			}
		
			// Shared substring
			if (sharedSubstring!=-1) {
				overhangLeft=sharedSubstringOverhangs[reference==nodeID1?0:2];
				p1=fromX+overhangLeft/QUANTUM;
				overhangRight=sharedSubstringOverhangs[reference==nodeID1?1:3];
				q1=fromX+(referenceLength-overhangRight)/QUANTUM;
				p2=p1;
				overhangLeft=sharedSubstringOverhangs[reference==nodeID1?2:0];
				from=p2-overhangLeft/QUANTUM;
				to=from+otherLength/QUANTUM;
				overhangRight=sharedSubstringOverhangs[reference==nodeID1?3:1];
				q2=to-overhangRight/QUANTUM;
				for (x=from; x<=to; x++) {
					for (y=0; y<INTERVAL_HEIGHT; y++) image.setRGB(x,lastY+y,otherColor);
				}
				Colors.drawTriangle(image,to,lastY,INTERVAL_HEIGHT,true,otherColor);
				label="SUBSTRING // "+overhangLeft+","+overhangRight;
				label+=" // "+IO.format(((double)Math.max(q1-p1+1,q2-p2+1))/Math.min(q1-p1+1,q2-p2+1))+"x";
				graphics.setColor(new Color(COLOR_TEXT));
				graphics.drawString(label,from+TEXT_OFFSET,lastY+INTERVAL_HEIGHT-TEXT_OFFSET);
				graphics.setColor(LINE_COLOR); graphics.setStroke(new BasicStroke(LINE_WIDTH));
				if (orientation==0) {
					graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,p2,lastY);
					graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,q2,lastY);
				}
				else {
					graphics.drawLine(p1,fromY+INTERVAL_HEIGHT,q2,lastY);
					graphics.drawLine(q1,fromY+INTERVAL_HEIGHT,p2,lastY);
				}
				lastY+=INTERVAL_HEIGHT+VERTICAL_SPACE;
				// Drawing inside the snippet
				graphics.setColor(new Color(COLOR_SNIPPET_GUIDES));
				graphics.setStroke(new BasicStroke(LINE_WIDTH_SNIPPET_GUIDES));
				graphics.drawLine(referencePointX+p1-fromX,referenceSnippetMinY,referencePointX+p1-fromX,referenceSnippetMaxY);
				graphics.drawLine(referencePointX+q1-fromX,referenceSnippetMinY,referencePointX+q1-fromX,referenceSnippetMaxY);
				graphics.drawLine(otherPointX+p2-from,otherSnippetMinY,otherPointX+p2-from,otherSnippetMaxY);
				graphics.drawLine(otherPointX+q2-from,otherSnippetMinY,otherPointX+q2-from,otherSnippetMaxY);
			}
		}
 	}
	
	
	/**
	 * Computes basic statistics on the edges and on the neighbors of node $nodeID$.
	 *
	 * @param neighborType output array: 0=alignment; 1=prefix; 2=suffix; 3=prefix+suffix;
	 * 4=substring; 5=single deletion; 6=short period; 7=long period;
	 * @param edgeType output array: 0=overlap prefix; 1=overlap suffix; 2=identical; 
	 * 3=containing; 4=contained; 5=insertion containing; 6=insertion contained;
	 * 7=shared substring.
	 */
	public static final void getEdgeStats(int nodeID, int[] neighborType, int[] edgeType) {
		int j;
		Node nodeTo;
		Edge edge;
		
		Math.set(neighborType,7,0); Math.set(edgeType,7,0);
		for (j=0; j<nNeighbors[nodeID]; j++) {
			edge=neighbors[nodeID][j];
			nodeTo=nodesArray[edge.getTo(nodeID)];
			if (nodeTo.type<=Constants.INTERVAL_DENSE_SINGLEDELETION) neighborType[nodeTo.type]++;
			else if (nodeTo.hasLongPeriod) neighborType[7]++;
			else neighborType[6]++;			
			if ( edge.overlap==Constants.OVERLAP_PREFIX_PREFIX || 
				 (edge.overlap==Constants.OVERLAP_PREFIX_SUFFIX && edge.nodeID1==nodeID) ||
				 (edge.overlap==Constants.OVERLAP_SUFFIX_PREFIX && edge.nodeID2==nodeID)
			   ) edgeType[0]++;
			if ( edge.overlap==Constants.OVERLAP_SUFFIX_SUFFIX || 
				 (edge.overlap==Constants.OVERLAP_PREFIX_SUFFIX && edge.nodeID2==nodeID) ||
				 (edge.overlap==Constants.OVERLAP_SUFFIX_PREFIX && edge.nodeID1==nodeID)
			   ) edgeType[1]++;
			if (edge.containment==Constants.CONTAINMENT_IDENTICAL) edgeType[2]++;
			else if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID2==nodeID) ||
				      (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID1==nodeID)
				    ) edgeType[3]++;
			else if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==nodeID) ||
				      (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID2==nodeID)
				    ) edgeType[4]++;
			if ( (edge.insertion==Constants.INSERTION_ONE_IN_TWO && edge.nodeID2==nodeID) ||
				 (edge.insertion==Constants.INSERTION_TWO_IN_ONE && edge.nodeID1==nodeID)
			   ) edgeType[5]++;
			else if ( (edge.insertion==Constants.INSERTION_ONE_IN_TWO && edge.nodeID1==nodeID) ||
				      (edge.insertion==Constants.INSERTION_TWO_IN_ONE && edge.nodeID2==nodeID)
			        ) edgeType[6]++;
			if (edge.sharedSubstring!=-1) edgeType[7]++;
		}
	}
	
	
	/**
	 * Stores in $out[i]$ the number of times relative position $i$ of $nodeID$ is the 
	 * beginning/end of an alignment with a node, and the corresponding position in such 
	 * node coincides with its first/last position (and is maximal if $isMaximal=true$). 
	 * Positions in $out$ are binned with bin size $quantum$.
	 *
	 * @return the last nonzero cell in $out$.
	 */
	public static final int getEdgeHistogram(int nodeID, int[] out, int quantum, boolean isMaximal) {
		final int DISTANCE_THRESHOLD = IO.quantum;
		int j;
		int to, p1, q1, referenceLength, overhangLeft, overhangRight, last;
		Node node, nodeTo;
		Edge edge;
		
		node=nodesArray[nodeID];
		referenceLength=node.length();
		Math.set(out,referenceLength/quantum+1,0);
		last=-1;
		for (j=0; j<nNeighbors[nodeID]; j++) {
			edge=neighbors[nodeID][j];
			to=edge.getTo(nodeID);
			nodeTo=nodesArray[to];
			
			// Overlap
			if (edge.overlap!=-1) {
				if ((edge.overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) {
					p1=0;
					q1=(referenceLength-edge.overhangs[0][nodeID==edge.nodeID1?0:1])/quantum;
					if (!isMaximal || nodeTo.isLeftMaximal) {
						out[q1]++;
						if (q1>last) last=q1;
					}
				}
				if ((edge.overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) {
					p1=0; 
					q1=(referenceLength-edge.overhangs[1][nodeID==edge.nodeID1?0:1])/quantum;
					if (!isMaximal || nodeTo.isRightMaximal) {
						out[q1]++;
						if (q1>last) last=q1;
					}
				}
				if ((edge.overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) {
					p1=edge.overhangs[2][nodeID==edge.nodeID1?0:1]/quantum;
					q1=referenceLength/quantum;
					if (!isMaximal || nodeTo.isLeftMaximal) {
						out[p1]++;
						if (p1>last) last=p1;
					}
				}
				if ((edge.overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) {
					p1=edge.overhangs[3][nodeID==edge.nodeID1?0:1]/quantum;
					q1=referenceLength/quantum;
					if (!isMaximal || nodeTo.isRightMaximal) {
						out[p1]++;
						if (p1>last) last=p1;
					}
				}
			}
	
			// Containment/insertion
			if (edge.containment==Constants.CONTAINMENT_IDENTICAL) {
				p1=0;
				q1=referenceLength/quantum;
				if (!isMaximal) {
					out[p1]++;
					out[q1]++;
				}
			}
			else if ( ((edge.containment==Constants.CONTAINMENT_ONE_IN_TWO||edge.insertion==Constants.INSERTION_ONE_IN_TWO) && nodeID==edge.nodeID1) || 
			          ((edge.containment==Constants.CONTAINMENT_TWO_IN_ONE||edge.insertion==Constants.INSERTION_TWO_IN_ONE) && nodeID==edge.nodeID2)
			        ) {
				p1=0;
				q1=referenceLength/quantum;
				if (!isMaximal) {
					out[p1]++;
					out[q1]++;
				}
			}
			else if ( ((edge.containment==Constants.CONTAINMENT_ONE_IN_TWO||edge.insertion==Constants.INSERTION_ONE_IN_TWO) && nodeID==edge.nodeID2) || 
			          ((edge.containment==Constants.CONTAINMENT_TWO_IN_ONE||edge.insertion==Constants.INSERTION_TWO_IN_ONE) && nodeID==edge.nodeID1)
			        ) {
				p1=edge.containmentOverhangLeft/quantum;
				q1=(referenceLength-edge.containmentOverhangRight)/quantum;
				if (edge.orientation==0) {
					if (!isMaximal || nodeTo.isLeftMaximal) {
						out[p1]++;
						if (p1>last) last=p1;
					}
					if (!isMaximal || nodeTo.isRightMaximal) {
						out[q1]++;
						if (q1>last) last=q1;
					}
				}
				else if (edge.orientation==1) {
					if (!isMaximal || nodeTo.isRightMaximal) {
						out[p1]++;
						if (p1>last) last=p1;
					}
					if (!isMaximal || nodeTo.isLeftMaximal) {
						out[q1]++;
						if (q1>last) last=q1;
					}
				}
			}
	
			// Shared substring
			if (edge.sharedSubstring!=-1) {
				overhangLeft=edge.sharedSubstringOverhangs[nodeID==edge.nodeID1?0:2];
				p1=overhangLeft/quantum;
				overhangRight=edge.sharedSubstringOverhangs[nodeID==edge.nodeID1?1:3];
				q1=(referenceLength-overhangRight)/quantum;
				if (edge.orientation==0) {
					if ( !isMaximal || 
						 (overhangLeft>DISTANCE_THRESHOLD && edge.sharedSubstringOverhangs[nodeID==edge.nodeID1?2:0]<=DISTANCE_THRESHOLD && nodeTo.isLeftMaximal)
				       ) {
						out[p1]++;
						if (p1>last) last=p1;
					}
					if ( !isMaximal || 
						 (overhangRight>DISTANCE_THRESHOLD && edge.sharedSubstringOverhangs[nodeID==edge.nodeID1?3:1]<=DISTANCE_THRESHOLD && nodeTo.isRightMaximal)
				       ) {
						out[q1]++;
						if (q1>last) last=q1;
					}
				}
				else if (edge.orientation==1) {
					if ( !isMaximal || 
						 (overhangLeft>DISTANCE_THRESHOLD && edge.sharedSubstringOverhangs[nodeID==edge.nodeID1?3:1]<=DISTANCE_THRESHOLD && nodeTo.isRightMaximal)
					   ) {
						out[p1]++;
						if (p1>last) last=p1;
					}
					if ( !isMaximal || 
						 (overhangRight>DISTANCE_THRESHOLD && edge.sharedSubstringOverhangs[nodeID==edge.nodeID1?2:0]<=DISTANCE_THRESHOLD && nodeTo.isLeftMaximal)
					   ) {
						out[q1]++;
						if (q1>last) last=q1;
					}
				}
			}
		}
		
		return last;
	}
	
	
	/**
	 * Stores in $histogram$ the binned distribution of clustering coefficient values for
	 * all nodes in the graph. The number of equal-size bins equals the size of
	 * $histogram$.
	 *
	 * Remark: the procedure does not consider self-edges, and it does not consider edges
	 * with $on=false$.
	 *
	 * @param tmpArray* temporary space, of size at least equal to the max degree of a 
	 * node;
	 * @param edgesSorted TRUE iff $neighbors[i]$ is already sorted by $nodeID$ for every 
	 * node $i$.
	 */
	public static final void getClusteringCoefficientHistogram(int[] histogram, boolean edgesSorted, int[] tmpArray1, int[] tmpArray2) {
		final double HISTOGRAM_QUANTUM = 1.0/histogram.length;
		int i;
		double coefficient;

		if (!edgesSorted) {
			for (i=0; i<nNodes; i++) {
				if (nNeighbors[i]<=1) continue;
				Edge.order=i;
				Arrays.sort(neighbors,0,nNeighbors[i]+1);
			}
		}
		Math.set(histogram,histogram.length-1,0);
		for (i=0; i<nNodes; i++) {
			coefficient=getClusteringCoefficient(i,tmpArray1,tmpArray2);
			if (coefficient==-1) continue;
			histogram[coefficient==1.0?histogram.length-1:(int)(coefficient/HISTOGRAM_QUANTUM)]++;
		}
	}
	
	
	/**
	 * Remark: the procedure assumes that $neighbors[i]$ is already sorted by $nodeID$.
	 *
	 * @param tmpArray* temporary space, of size at least equal to the max degree of a 
	 * node;
	 * @return -1 if the coefficient could not be computed.
	 */
	public static final double getClusteringCoefficient(int i, int[] tmpArray1, int[] tmpArray2) {
		int j, k;
		int to, last1, last2, nEdges, other;
		double out;
		
		if (nNeighbors[i]<=1) return -1;
		last1=-1;
		for (j=0; j<nNeighbors[i]; j++) {
			if (!neighbors[i][j].on) break;
			to=neighbors[i][j].getTo(i);
			if (to==i) continue;
			tmpArray1[++last1]=to;
		}
		if (last1<=0) return -1;
		nEdges=0;		
		for (j=0; j<nNeighbors[i]; j++) {
			if (!neighbors[i][j].on) break;
			to=neighbors[i][j].getTo(i);
			if (to==i) continue;
			last2=-1;
			for (k=0; k<nNeighbors[to]; k++) {
				if (!neighbors[to][k].on) break;
				other=neighbors[to][k].getTo(to);
				if (other==i || other==j) continue;
				tmpArray2[++last2]=other;
			}			
			if (last2>=0) nEdges+=Math.setIntersectionSize(tmpArray1,0,last1,tmpArray2,0,last2,to);
		}
		out=((double)nEdges)/((last1+1)*last1);
		if (IO.CONSISTENCY_CHECKS && out>1.0) {
			System.err.println("getClusteringCoefficient> ERROR: coefficient="+out+">1.0 last1="+last1+" nEdges="+nEdges);		
			System.exit(1);
		}
		return out;
	}
	
	
	/**
	 * Stores in $histogram$ the binned distribution of degree values for all nodes in the
	 * graph. The number of equal-size bins equals the size of $histogram$. 
	 *
	 * Remark: the procedure does not consider self-loops, and it does not consider edges
	 * with $on=false$.
	 */
	public static final void getDegreeHistogram(int[] histogram, int maxDegree, boolean edgesSorted) {
		final double HISTOGRAM_QUANTUM = ((double)maxDegree)/histogram.length;
		int i;
		int nEdges;
		
		Math.set(histogram,histogram.length-1,0);
		for (i=0; i<nNodes; i++) {
			nEdges=getDegree(i,edgesSorted);
			histogram[nEdges==maxDegree?histogram.length-1:(int)(nEdges/HISTOGRAM_QUANTUM)]++;
		}
	}
	
	
	/**
	 * Remark: the procedure does not take into account OFF edges and self-loops.
	 *
	 * @param edgesSorted TRUE iff all ON edges are at the beginning of $neighbors[i]$.
	 */
	public static final int getDegree(int i, boolean edgesSorted) {
		int j;
		int nEdges;
		
		nEdges=0;
		for (j=0; j<nNeighbors[i]; j++) {
			if (!neighbors[i][j].on) {
				if (edgesSorted) break;
				else continue;
			}
			if (neighbors[i][j].getTo(i)!=i) nEdges++;
		}
		return nEdges;
	}
	
	
	public static class PriorityQueueElement implements Comparable {
		int id, size;
	
		public PriorityQueueElement(int i, int s) {
			id=i;
			size=s;
		}
	
		public boolean equals(Object other) {
			PriorityQueueElement otherElement = (PriorityQueueElement)other;
			return id==otherElement.id;
		}
	
		public int compareTo(Object other) {
			PriorityQueueElement otherElement = (PriorityQueueElement)other;
			if (size<otherElement.size) return -1;
			else if (size>otherElement.size) return 1;
			return 0;
		}
	}
	
	
	/**
	 * @param out output array, containing the number of short-period (0) and long-period
	 * (1) nodes in $nodesArray$.
	 */
	public static final void nPeriodicNodes(int[] out) {
		int i;
		int nShort, nLong;
		
		nShort=0; nLong=0;
		for (i=0; i<nNodes; i++) {
			if (nodesArray[i].type!=Constants.INTERVAL_PERIODIC) continue;
			if (nodesArray[i].hasLongPeriod) nLong++;
			else nShort++;
		}
		out[0]=nShort; out[1]=nLong;
	}
	
	
}