package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.util.PriorityQueue;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Stream;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Factorize;


/**
 * Tries to split a connected component into subcomponents, based on modularity of 
 * containment edges only. The resulting subcomponents might not be connected when
 * considering containment edges only.
 *
 * Remark: changing the clustering algorithm can change the kernels built in the following
 * steps, since some nodes that are contained and in a cluster with one setting of the 
 * parameters, might become not contained and in another cluster with another setting, and
 * thus be used as kernels.
 *
 * Remark: the graphs in input to this stage look very similar to those  in
 * \cite{estill2010taxonomy}, which are however derived from a matrix of alignment 
 * scores between all pairs of the 5' long terminal repeats of LTR transposons (i.e. not 
 * between the full copies of the entire transposons). Transposon occurrences are detected
 * by other repeat-finding programs, with the assembled genome as input. It is interesting
 * that interval graphs built from a sequenced genome and from likely true repeats look 
 * similar to the interval graphs we observe. Similar graphs appear also in  
 * \cite{novak2010graph}, but their nodes are short reads and their edges are just suffix- 
 * prefix overlaps.
 *
 * Remark: this program might print some tags files, but the intervals stored in there 
 * are not going to be merged with the intervals produced by the steps downstream (as
 * specified by $IntervalGraphStep3.printKernelTags()$). This is done to highlight
 * such intervals.
 */
public class IntervalGraphStep2 {	
	/**
	 * Parameters of the pipeline
	 */
	public static final String GRAPH_SUFFIX = ".graph";
	private static final int MIN_PEELING_SIZE = 1000;
	private static final int ALIGNMENT_OVERLAP_THRESHOLD = 5*IO.quantum;  // Used for deciding overlaps with alignment intervals. Chosen arbitrarily.
	
	/**
	 * Encoding of the recursive peeling tree
	 */
	private static PeelingNode peelingTree;
	
	/**
	 * Temporary space
	 */
	public static Point[] tmpPoints;
	public static int lastTmpPoint;
	private static int[] tmpArray = new int[4];
	// Seed set expansion variables
	private static int[] degrees, clusterField;
	private static PriorityQueue<IntervalGraph.SeedExtensionPair> queue;
	private static Stream stream;
	private static IntervalGraph.AlignmentPair tmpPair = new IntervalGraph.AlignmentPair();
	private static IntervalGraph.Edge tmpEdge = new IntervalGraph.Edge();
	
	/**
	 * Temporary space for $getComponents_clustering()$.
	 */
	private static int[] clusterSize, clusterSizePrime, clusterTmp, clusterTmpPrime;
	private static int[] lastOutNeighbor, nInNeighbors, components;
	private static double[] deltas;
	private static int[][] outNeighbors;
	
	
	public static void main(String[] args) throws IOException {
		final int MIN_SIZE_FOR_CLUSTERING = 50;  // In number of nodes. Arbitrary.
		final int MAX_SUBGRAPH_EDGES = 300000;  // Maximum size of a graph to be considered small. Arbitrary.
		final int ARTIFICIAL_KERNEL = Math.POSITIVE_INFINITY;
		
		final String GRAPH_FILE = args[0];
		final String GRAPH_FILE_PREFIX = GRAPH_FILE.substring(0,GRAPH_FILE.length()-GRAPH_SUFFIX.length());
		final String COMPONENTS_DIR = args[1];
		IO.coverage=Integer.parseInt(args[2]);
		IO.minOccurrencesInGenome=Integer.parseInt(args[3]);
		IO.minRepeatCoverage=(IO.minOccurrencesInGenome*IO.coverage)-1;
		Alignments.minAlignmentLength=Integer.parseInt(args[4]);
		final boolean PRINT_DOT = Integer.parseInt(args[5])==1;
		final String TAGS_DIR = args[6];
		final String TAGS_PREFIX = args[7];
		final int N_READS = Integer.parseInt(args[8]);
		
		final String READS_IDS = GRAPH_FILE_PREFIX+"-"+IO.READS_IDS;
		final String READS_LENGTHS = GRAPH_FILE_PREFIX+"-"+IO.READS_LENGTHS;
		final String READS_PHRED = GRAPH_FILE_PREFIX+"-"+IO.READS_PHRED_PREFIX+IO.getPhredSuffix(GRAPH_FILE_PREFIX+"-"+IO.READS_PHRED_PREFIX);
		final String READS_SHORTPERIOD = GRAPH_FILE_PREFIX+"-"+IO.READS_SHORTPERIOD;
		final String ALIGNMENTS_FILE = GRAPH_FILE_PREFIX+"-"+IO.ALIGNMENTS_FILE;
		final String ALIGNMENTS_SHORTPERIOD = GRAPH_FILE_PREFIX+"-"+IO.ALIGNMENTS_SHORTPERIOD;
		final int MIN_CLIQUE_SIZE = 5;  // Arbitrary
		
		boolean graphIsSmall;
		int i;
		int nEdges, nComponents, maxDegree;
		BufferedWriter bw;
		IntervalGraph.Node node;
		int[] components, new2old, subgraph;
		
		// Loading graph
		IO.initialize();
		System.err.println("STEP2> Loading interval graph...");
		new2old=IntervalGraph.deserialize(GRAPH_FILE,true,true);
		maxDegree=IntervalGraph.turnOffEdges(-2,true);  // All edges
		if (IO.CONSISTENCY_CHECKS) {
			System.err.println("STEP2> Consistency checks started...");
			components=IntervalGraph.getComponents();
			if (components[0]!=1) {
				System.err.println("STEP2> ERROR: the nodes in the input graph have "+components[0]+" distinct component tags.");
				System.exit(1);
			}
			IntervalGraph.checkConsistency(0,IntervalGraph.nNodes-1,new int[maxDegree*10]);
			System.err.println("STEP2> Consistency checks passed");
		}
		IntervalGraph.printNew2OldArray(new2old,GRAPH_FILE_PREFIX+".new2old");
		new2old=null;
		System.err.println("STEP2> Interval graph loaded:");
		IntervalGraph.getNodeStats();
		nEdges=IntervalGraph.getEdgeStats();
		graphIsSmall=nEdges<=MAX_SUBGRAPH_EDGES;
		System.err.println("STEP2> nNodes="+IntervalGraph.nNodes+" nEdges="+nEdges);
		IntervalGraph.printNodeStats();
		IntervalGraph.printEdgeStats(nEdges);
		System.err.println();	
		if (IntervalGraph.nNodes<MIN_SIZE_FOR_CLUSTERING) return;
		// Remark: $avgDiffs$ in the graph is currently a sum, not an avg., and nodes are
		// sorted by $read,start$.		
		
		// Clustering
		System.err.println("STEP2> Clustering...");
		subgraph = new int[IntervalGraph.nNodes];
		for (i=0; i<IntervalGraph.nNodes; i++) subgraph[i]=i;
		degrees = new int[IntervalGraph.nNodes];
		clusterField = new int[IntervalGraph.nNodes];
		stream = new Stream(256);
		nComponents=0;
		if (graphIsSmall) {
			nComponents=getComponents_clustering(subgraph,0,IntervalGraph.nNodes-1,0,degrees,clusterField,stream,tmpArray,true,true,MIN_CLIQUE_SIZE);
		}
		else {
			// ..........................
			nComponents=getComponents_clustering(subgraph,0,IntervalGraph.nNodes-1,0,degrees,clusterField,stream,tmpArray,true,false,MIN_CLIQUE_SIZE);
		}
		System.err.println("STEP2> "+nComponents+" clusters found");
		
		// Outputting
		if (PRINT_DOT) {
			System.err.println("STEP2> Printing the whole graph (with cluster tags) to <"+GRAPH_FILE_PREFIX+"-components.dot>...");
			IntervalGraph.toDot(GRAPH_FILE_PREFIX+"-components.dot",null,null,Math.POSITIVE_INFINITY);
		}
		if (nComponents>1) {
			System.err.println("STEP2> Saving all clusters...");
			IntervalGraph.serializeFrequentComponents(null,nComponents,true,true,COMPONENTS_DIR,GRAPH_SUFFIX,TAGS_DIR,TAGS_PREFIX);  // Strict, since in the following steps of the pipeline we will use just nodes with one component for building kernels.
			System.err.println("STEP2> done");
			System.err.println("STEP2> Filtering read and alignment files by component...");
			Reads.loadReadIDs(READS_IDS,N_READS);
			IntervalGraph.filterAlignmentsAndReads(null,nComponents,true,READS_LENGTHS,READS_PHRED,READS_SHORTPERIOD,ALIGNMENTS_FILE,ALIGNMENTS_SHORTPERIOD,COMPONENTS_DIR);
			System.err.println("STEP2> done");
		}		
		// Remark: $avgDiffs$ in the serialized files is a sum, not an avg., and nodes are
		// sorted by $read,start$.		
		
			
		/*
		bw = new BufferedWriter(new FileWriter(COMPONENTS_DIR+"/toSplitFurther.txt"),IO.BUFFER_SIZE);
		components=getComponents(bw);
		bw.close();
		nComponents=components.length;
		System.err.println("STEP2> "+nComponents+" clusters found ("+(peelingTree.size()-1)+" components found by peeling)");
		if (nComponents<=1) return;
		printPeelingTree(COMPONENTS_DIR+"/peelingTree.dot");
		System.err.println("STEP2> Printing the whole graph (with component tags) to <"+GRAPH_FILE_PREFIX+"-components.dot>...");
		*/
			

	
	}


	/**
	 * ----------->Tries to cluster the nodes of a subgraph that cannot be further split by peeling.
	 * Then, marks the clusters as components, propagates component tags to all nodes not 
	 * in a cluster, and disconnects nodes with more than one component tag.
	 * Breaking large connected components induced by overlap edges only, into hundreds of
	 * clusters by modularity, was already done in \cite{novak2010graph}, but from short 
	 * reads. Modularity has been used also in \cite{guo2017replong} to cluster a graph of
	 * overlap edges only, from long PacBio reads.
	 *
	 * Remark: the nodes that belong to the same cluster at the end of the procedure might
	 * not necessarily form a connected subgraph.
	 *
	 * Remark: the procedure assumes that all nodes in $subgraph[first..last]$ are 
	 * assigned one and the same component ID before the procedure is called.
	 *
	 * Remark: the procedure allocates internally some matrices of size proportional to 
	 * the number of clusters in the subgraph. This number is hard to estimate by the 
	 * caller, and it is expected to be small. The procedure also allocates matrices of
	 * size proportional to the number of nodes in the subgraph. This size is known by the
	 * caller, who could pass such matrices to the function. However, the number of
	 * matrices is large, so we make the function allocate them to reduce clutter in its 
	 * interface. 
	 *
     * Remark: if the subgraph is too large, the procedure builds maximal cliques just 
	 * from a random subset of nodes, and it does not refine clusters.
	 *
	 * Remark: when all the nodes of the subgraph are assigned to the same component,
	 * the kernels built later might assemble intervals in different clusters, merging 
	 * distinct repeats.
	 *
	 * @param firstComponent component IDs are assigned starting from $firstComponent$,
	 * included;
	 * @param degrees,clusterField temporary space, of size at least equal to the total 
	 * number of nodes in the interval graph;
	 * @param out temporary space, of size at least 2;
	 * @param stream temporary space;
	 * @param minCliqueSize min clique size for $IntervalGraph.getSeeds()$;
	 * @return 0 if $subgraph$ is not split, otherwise the number of components found.
	 */
	public static final int getComponents_clustering(int[] subgraph, int first, int last, int firstComponent, int[] degrees, int[] clusterField, Stream stream, int[] out, boolean sorted, boolean graphIsSmall, int minCliqueSize) throws IOException {
		final double SAMPLING_RATE_LARGE = 0.01;  // To limit the time of $IntervalGraph.getMaximalCliques_impl()$ in large graphs. Arbitrary.
		final double SAMPLING_RATE = graphIsSmall?-1:SAMPLING_RATE_LARGE;
		final int MAX_N_CLIQUES = 100000;  // To limit the time of $IntervalGraph.getMaximalCliques_impl()$. Arbitrary.
		final int MAX_CLIQUE_SIZE = -1;//40;  // To limit the time of $IntervalGraph.getMaximalCliques_impl()$. Arbitrary.
		final double UNDO_RATE = 0.5;  // Arbitrary
		final int GROWTH_RATE = 10;  // Arbitrary
		
		int i, j;
		int node, nNodes, nEdges, lastSeed, nClusters, nodesWithCluster;
		int[] seeds, clusters;
		IntervalGraph.Edge edge;		
		
		// Building seeds
		nNodes=last-first+1;
		if (!sorted) IntervalGraph.sortEdges(subgraph,first,last);
		nEdges=IntervalGraph.getNOnEdges(subgraph,first,last,true);		
		IntervalGraph.getSeeds(subgraph,first,last,minCliqueSize,stream,MAX_CLIQUE_SIZE,SAMPLING_RATE,MAX_N_CLIQUES,out);
		lastSeed=out[0]; nClusters=out[1];
		if (nClusters<=1) return 0;
		seeds=IntervalGraph.seeds; clusters=IntervalGraph.clusters;
		
		// Extending seeds
		if (clusterSize==null || clusterSize.length<nClusters) clusterSize = new int[nClusters];
		if (clusterSizePrime==null || clusterSizePrime.length<nClusters) clusterSizePrime = new int[nClusters];
		if (deltas==null || deltas.length<nClusters) deltas = new double[nClusters];
		System.err.println("getComponents_clustering> 2 extendSeeds started. nClusters="+nClusters);
		if (queue==null) queue = new PriorityQueue<IntervalGraph.SeedExtensionPair>(nNodes);
		IntervalGraph.extendSeeds(subgraph,first,last,lastSeed,nClusters,clusterSize,nEdges,deltas,queue,clusterField,degrees,true,-1);
		System.err.println("extendSeeds ended");

		// Refining clusters
		if (graphIsSmall) {
			if (clusterTmp==null || clusterTmp.length<nNodes) clusterTmp = new int[nNodes];
			if (clusterTmpPrime==null || clusterTmpPrime.length<nNodes) clusterTmpPrime = new int[nNodes];
			if (outNeighbors==null || outNeighbors.length<nNodes) outNeighbors = new int[nNodes][GROWTH_RATE];
			if (lastOutNeighbor==null || lastOutNeighbor.length<nNodes) lastOutNeighbor = new int[nNodes];
			if (nInNeighbors==null || nInNeighbors.length<nNodes) nInNeighbors = new int[nNodes];
			if (components==null || components.length<nNodes) components = new int[nNodes];
			nClusters=IntervalGraph.refineClusters(subgraph,first,last,nClusters,clusterSize,clusterSizePrime,nEdges,queue,deltas,clusterTmp,clusterTmpPrime,outNeighbors,lastOutNeighbor,nInNeighbors,components,clusterField,degrees,true);
		}
		else {
			nClusters=IntervalGraph.refineClusters_large(subgraph,first,last,nClusters,clusterSize,clusterSizePrime,nEdges,queue,deltas,clusterField,degrees,true);
		}
		System.err.println("refineClusters ended: "+nClusters+" final clusters assigned to the subgraph");
		if (nClusters==1) return 0;
		nodesWithCluster=0;  // Undoing clustering if the total number of nodes assigned to any cluster is smaller than the total number of nodes assigned to no cluster.
		for (i=0; i<nClusters; i++) nodesWithCluster+=clusterSize[i];
		if (nodesWithCluster<=nNodes*UNDO_RATE) return 0;
		
		// Converting $clusterField$ into the $components$ field of nodes
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			if (clusterField[node]<0) IntervalGraph.nodesArray[node].lastComponent=-1;
			else {
				IntervalGraph.nodesArray[node].lastComponent=0;
				if (IntervalGraph.nodesArray[node].components==null || IntervalGraph.nodesArray[node].components.length==0) IntervalGraph.nodesArray[node].components = new int[1];
				IntervalGraph.nodesArray[node].components[0]=firstComponent+clusterField[node];
			}
		}
		// Remark: calling $IntervalGraph.propagateComponentTags()$ is not a good idea in
		// practice, since it adds to the clusters nodes that are in the middle of
		// filaments and that it might be better to keep unassigned.
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		/*
		// Setting to OFF all containment edges adjacent to nodes with more than one
		// component.
		for (i=first; i<=last; i++) {
			node=subgraph[i];
			if (IntervalGraph.nodesArray[node].lastComponent==0) continue;
			for (j=0; j<IntervalGraph.nNeighbors[node]; j++) {
				edge=IntervalGraph.neighbors[node][j];
				if (!edge.on) continue;
				if (edge.containment==-1) continue;
				edge.on=false;
			}
		}
		
		// Writing to $toSplitFurther$
		if (!graphIsSmall && toSplitFurther!=null) {
			for (i=0; i<nClusters; i++) toSplitFurther.write((firstComponent+i)+"\n");
		}
		*/
		
		
System.err.println("getComponents_clustering> returns nClusters="+nClusters);
		return nClusters;
	}






	
	

	
	





	
	
	
	
	// --------------------------------- TYPES OF EDGE -----------------------------------
	
	/**
	 * Decides whether the edge induced by an alignment pair belongs to the Step 2 graph
	 * or not. 
	 *
	 * Remark: in general, at this step of the pipeline, alignments might contain regions 
	 * that flank intervals, e.g. if intervals have been trimmed for low quality. However,
	 * some procedures that decide edge type work correctly only if alignments are approx.
	 * contained in both intervals. For example, containment procedures do not handle the 
	 * case in which an interval is fully contained in a long alignment in readA, and the
	 * readB side of the alignment is fully contained in a readB interval. 
	 * One could rewrite all procedures to handle any configuration of alignments and 
	 * intervals: we leave this to the future due to lack of time.
	 *
	 * Remark: above: if one calls $AlignmentPair.trim()$ on an aligment that is not 
	 * approx. contained in both intervals, the procedures that decide edge type might 
	 * give wrong results.
	 *
	 * Remark: even after keeping only Step-2 edges, it is possible that distinct 
	 * repeats belong to the same connected component, since in general distinct 
	 * repeats can contain similar substrings. Distinct repeats can share a substring 
	 * because: (1) they straddle; (2) one is contained in the other; (3) they just 
	 * share a substring, and such substring has erroneously not been marked as a 
	 * distinct repeat. See procedures $inStep2_containment$, $inStep2_overlap$,
	 * $inStep2_insertion$ for details.
	 *
	 * @return -1: if the edge does not belong to the graph in Step 2. Otherwise, a
	 * suitable code denoting the type of edge.
	 */
	public static final int inStep2(IntervalGraph.AlignmentPair pair) {
		int i;
		
		i=inStep2_containment(pair);
		if (i>=0) return i;
		i=inStep2_overlap(pair);
		if (i>=0) return i;
		i=inStep2_insertion(pair);
		if (i>=0) return i;
		i=inStep2_sharedSubstring(pair);
		if (i>=0) return i;
		return -1;
	}
	
	
	/**
	 * Decides whether the edge induced by an alignment pair belongs to the Step 2 graph 
	 * or not, based on whether the two intervals can be (not necessarily proper) 
	 * substrings of one another.
	 *
	 * For each interval type, the procedure considers just: (1) all possible 
	 * substrings which can be obtained by inserting long low-quality regions into a 
	 * copy; (2) all substrings that are allowed by the replication modes of the 
	 * interval. Specifically, the procedure does not consider an interval that is 
	 * both left- and right-maximal as a candidate substring of another interval: see 
	 * procedure $isInsertion$ for details.
	 *
	 * Remark: even if we just consider these edges, a connected component could still
	 * contain distinct repeats that happen to be adjacent in the genome. In the 
	 * simple example of procedure $Factorize.alignments2intervals()$, $X$ could occur
	 * in read 2 immediately followed by a region of low quality, $Z$ could occur in 
	 * read 3 immediately preceded by a region of low quality, and $Z$ might not be 
	 * split into $X_2$ and $Y_1$ in read 3.
	 *
	 * Remark: let $W$ be a dense substring of substring type in readA, and let $\ell$
	 * be the minimum length of a substring of $W$. Copying to readB a substring $V$
	 * of $W$ can produce, in readB: (1) an alignment interval (if $|V|$ is
	 * $Alignments.minAlignmentLength$); (2) a dense substring of prefix 
	 * (respectively, suffix) type (if $V$ is a suffix (respectively, prefix) of $W$ 
	 * of length $\ell$, or if the distance between $V$ and an end of $W$ is small
	 * enough compared to $\ell$); (3) a dense substring of prefix and suffix type 
	 * (if $V$ is an infix of $W$ of length $\ell$); (4) a dense substring of 
	 * substring type (if $V$ is a long infix of $W$). This does not depend on the 
	 * left- and right-maximality of $V$ in readB.
	 *
	 * Remark: alignments that involve two entire end-to-end intervals are marked as
	 * $CONTAINMENT_IDENTICAL$ even if the types of interval differ. E.g. because of 
	 * nontransitivity, two diverged instances of a repeat might not have alignment piles 
	 * with similar shape.
	 *
	 * @return -1 if no substring relation can be detected. Otherwise, a code
	 * $CONTAINMENT_*$.
	 */
	public static final int inStep2_containment(IntervalGraph.AlignmentPair pair) {
		final int type1 = pair.node1.type;
		final int type2 = pair.node2.type;
		
		// Identity
		if ( ( (type1!=Constants.INTERVAL_PERIODIC && type2!=Constants.INTERVAL_PERIODIC) ||
			   (type1==Constants.INTERVAL_PERIODIC && type2==Constants.INTERVAL_PERIODIC && pair.node1.hasLongPeriod==pair.node2.hasLongPeriod)
			 ) &&
			 Math.abs(pair.alignmentStart1,pair.node1.start)<=IntervalGraph.IDENTITY_THRESHOLD && 
		     Math.abs(pair.alignmentEnd1-1,pair.node1.end)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(pair.alignmentStart2,pair.node2.start)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(pair.alignmentEnd2-1,pair.node2.end)<=IntervalGraph.IDENTITY_THRESHOLD
		   ) return Constants.CONTAINMENT_IDENTICAL;
		
		// Proper containment
		if (type1==Constants.INTERVAL_ALIGNMENT) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isContained_alignment_alignment(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isContained_alignment_prefixOrSuffix(pair,true);			
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isContained_alignment_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX || type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isContained_alignment_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isContained_alignment_substring(pair);
			if (type2==Constants.INTERVAL_PERIODIC) return isContained_periodic_alignment(pair);
		}
		if (type1==Constants.INTERVAL_DENSE_PREFIX) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isContained_alignment_prefixOrSuffix(pair,true);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isContained_prefix_prefixOrSuffix(pair,true);			
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isContained_prefix_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX || type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isContained_prefix_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isContained_substring_prefixOrSuffix(pair);
		}
		if (type1==Constants.INTERVAL_DENSE_SUFFIX) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isContained_alignment_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isContained_prefix_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isContained_suffix_suffix(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX || type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isContained_suffix_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isContained_substring_prefixOrSuffix(pair);
		}
		if (type1==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isContained_alignment_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isContained_prefix_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isContained_suffix_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX) return isContained_prefixSuffix_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isContained_substring_prefixOrSuffix(pair);  // Same criterion for prefixSuffix substrings and for just prefix or just suffix substrings
			if (type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isContained_prefixSuffix_singleDeletion(pair);
		}
		if (type1==Constants.INTERVAL_DENSE_SUBSTRING) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isContained_alignment_substring(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isContained_substring_prefixOrSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isContained_substring_prefixOrSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX) return isContained_substring_prefixOrSuffix(pair);  // Same criterion for prefixSuffix substrings and for just prefix or just suffix substrings
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isContained_substring_substring(pair);
			if (type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isContained_substring_singleDeletion(pair);
		}
		if (type1==Constants.INTERVAL_DENSE_SINGLEDELETION) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isContained_alignment_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isContained_prefix_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isContained_suffix_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX) return isContained_prefixSuffix_singleDeletion(pair);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isContained_substring_singleDeletion(pair);
			if (type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isContained_substring_substring(pair);  // Same criteria as for the substring-substring case
		}
		if (type1==Constants.INTERVAL_PERIODIC) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isContained_periodic_alignment(pair);
			if (type2==Constants.INTERVAL_PERIODIC) return isContained_periodic_periodic(pair);
		}
		return -1;
	}

	
	/**
	 * @param containmentType except for $CONTAINMENT_IDENTICAL$ and $INSERTION_BOTH$, can
	 * be any containment or insertion type;
	 * @param out 0: left containment/insertion overhang; 1: right containment/insertion
	 * overhang.
	 */
	public static final int inStep2_getContainmentSubtype(int containmentType, IntervalGraph.AlignmentPair pair, int[] out) {
		final int alignmentStart, alignmentEnd, start, end;
		
		if (containmentType==Constants.CONTAINMENT_ONE_IN_TWO || containmentType==Constants.INSERTION_ONE_IN_TWO) {
			alignmentStart=pair.alignmentStart2;
			alignmentEnd=pair.alignmentEnd2-1;
			start=pair.node2.start;
			end=pair.node2.end;
		}
		else if (containmentType==Constants.CONTAINMENT_TWO_IN_ONE || containmentType==Constants.INSERTION_TWO_IN_ONE) {
			alignmentStart=pair.alignmentStart1;
			alignmentEnd=pair.alignmentEnd1-1;
			start=pair.node1.start;
			end=pair.node1.end;
		}
		else return -1;
		out[0]=alignmentStart>start?alignmentStart-start:0;	
		out[1]=end>alignmentEnd?end-alignmentEnd:0;
		if (Math.abs(alignmentStart,start)<=IntervalGraph.IDENTITY_THRESHOLD) return Constants.CONTAINMENT_PREFIX;
		else if (Math.abs(alignmentEnd,end)<=IntervalGraph.IDENTITY_THRESHOLD) return Constants.CONTAINMENT_SUFFIX;
		else return Constants.CONTAINMENT_SUBSTRING;
	}

	
    /** 
	 * Decides whether the edge induced by an alignment pair belongs to the Step 2 graph 
	 * or not, based on whether the two intervals can have a suffix-prefix overlap.
	 * 
	 * @return -1 if no overlap can be detected; otherwise, a code $OVERLAP_*$.
	 */
	public static final int inStep2_overlap(IntervalGraph.AlignmentPair pair) {
		final int type1 = pair.node1.type;
		final int type2 = pair.node2.type;
		
		if (type1==Constants.INTERVAL_ALIGNMENT) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isOverlap_alignment_alignment(pair);
			if (type2==Constants.INTERVAL_PERIODIC) return isOverlap_periodic_alignment(pair);
		}
		if (type1==Constants.INTERVAL_DENSE_PREFIX) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isOverlap_prefix_alignment(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isOverlap_prefix_prefixOrSuffix(pair,true);			
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isOverlap_prefix_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX || type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isOverlap_prefixSuffix_prefixOrSuffix(pair,true);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isOverlap_substring_prefixOrSuffix(pair,true);
		}
		if (type1==Constants.INTERVAL_DENSE_SUFFIX) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isOverlap_suffix_alignment(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isOverlap_prefix_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isOverlap_suffix_suffix(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX || type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isOverlap_prefixSuffix_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isOverlap_substring_prefixOrSuffix(pair,false);
		}
		if (type1==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isOverlap_prefixSuffix_alignment(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isOverlap_prefixSuffix_prefixOrSuffix(pair,true);
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isOverlap_prefixSuffix_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX) return isOverlap_simple(pair);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isOverlap_substring_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isOverlap_simple(pair);
		}
		if (type1==Constants.INTERVAL_DENSE_SUBSTRING) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isOverlap_prefixSuffix_alignment(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isOverlap_substring_prefixOrSuffix(pair,true);
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isOverlap_substring_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX) return isOverlap_substring_prefixSuffix(pair);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return isOverlap_substring_substring(pair);
			if (type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return -1;
		}
		if (type1==Constants.INTERVAL_DENSE_SINGLEDELETION) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isOverlap_prefixSuffix_alignment(pair);
			if (type2==Constants.INTERVAL_DENSE_PREFIX) return isOverlap_prefixSuffix_prefixOrSuffix(pair,true);
			if (type2==Constants.INTERVAL_DENSE_SUFFIX) return isOverlap_prefixSuffix_prefixOrSuffix(pair,false);
			if (type2==Constants.INTERVAL_DENSE_PREFIXSUFFIX) return isOverlap_simple(pair);
			if (type2==Constants.INTERVAL_DENSE_SUBSTRING) return -1;
			if (type2==Constants.INTERVAL_DENSE_SINGLEDELETION) return isOverlap_simple(pair);
		}
		if (type1==Constants.INTERVAL_PERIODIC) {
			if (type2==Constants.INTERVAL_ALIGNMENT) return isOverlap_periodic_alignment(pair);
			if (type2==Constants.INTERVAL_PERIODIC) return isOverlap_periodic_periodic(pair);
		}
		return -1;
	}
	
	
	/**
	 * Remark: a short-period periodic interval can be a substring of a transposon, and a 
	 * diverged substring of a short-period periodic interval might become an independent 
	 * repeat that aligns just to a substring of the original periodic interval. Such
	 * relationships are not marked as insertions.
	 */
	private static final int inStep2_insertion(IntervalGraph.AlignmentPair pair) {
		final boolean hasLongPeriod1, hasLongPeriod2;
		final int orientation, type1, type2, length1, length2;
		
		orientation=pair.orientation?0:1;
		type1=pair.node1.type; type2=pair.node2.type;
		hasLongPeriod1=pair.node1.hasLongPeriod; hasLongPeriod2=pair.node2.hasLongPeriod; 
		length1=pair.node1.length(); length2=pair.node2.length();
		if (length2<length1-IntervalGraph.IDENTITY_THRESHOLD && !Constants.compareTypes(type2,type1,hasLongPeriod2,hasLongPeriod1,orientation)) return -1;
		if (length1<length2-IntervalGraph.IDENTITY_THRESHOLD && !Constants.compareTypes(type1,type2,hasLongPeriod1,hasLongPeriod2,orientation)) return -1;
		
		return isInsertion(pair);
	}
	
	
	private static final int inStep2_sharedSubstring(IntervalGraph.AlignmentPair pair) {
		return isSharedSubstring(pair);
	}
	
	
	/**
	 * Ensures that $edge$, which is assumed to correspond to just one alignment, belongs 
	 * to Step 2. If it belongs to Step 2, but its type is not correct, the procedure 
	 * fixes it. If it doesn't belong to Step 2, or if its type is correct, the procedure 
	 * returns without doing anything.
	 *
	 * @return -1: the edge does not belong to the graph in Step 2; 0: the edge is already
	 * of the correct type and is not altered; 1: the edge was of the wrong type and is 
	 * now fixed.
	 */
	public static final int inStep2_checkEdge(IntervalGraph.Edge edge) {
		int tmp, last;
		
		// Checking $edge$
		tmpPair.node1=IntervalGraph.nodesArray[edge.nodeID1];
		tmpPair.node2=IntervalGraph.nodesArray[edge.nodeID2];
		last=edge.projectOnto(edge.nodeID1,tmpArray,true);
		if (last==-1) return -1;
		tmpPair.alignmentStart1=tmpPair.node1.start+tmpArray[0];
		tmpPair.alignmentEnd1=tmpPair.node1.start+tmpArray[1];
		last=edge.projectOnto(edge.nodeID2,tmpArray,true);
		if (last==-1) return -1;
		tmpPair.alignmentStart2=tmpPair.node2.start+tmpArray[0];
		tmpPair.alignmentEnd2=tmpPair.node2.start+tmpArray[1];
		tmpPair.orientation=edge.orientation==0;
		tmpPair.diffs=edge.avgDiffs;
		tmpPair.inCanonicalForm();
		tmp=inStep2(tmpPair);
		if (tmp==-1) return -1;
		else if (tmp==edge.containment || tmp==edge.overlap || tmp==edge.insertion || tmp==edge.sharedSubstring) return 0;

		// Fixing $edge$
		tmpEdge.set(tmpPair);
		tmpEdge.clone(edge);
		return 1;
	}
	
	
	/**
	 * Checks whether each alignment type in $edge$ is correct according to Step 2. If 
	 * not, it assigns a shared substring type to the corresponding alignment.
	 *
	 * Remark: the procedure disregards $maxOverhangs$ and considers only $overhangs$,
	 * i.e. it considers just longest overlaps.
	 *
	 * @return TRUE iff $edge$ has been modified by the procedure.
	 */
	public static final boolean removeEdgeTypes(IntervalGraph.Edge edge) {
		boolean out;
		int overlap, last;
		
		out=false;
		if (edge.orientation==2 || edge.orientation==-1) {
			out=removeEdgeTypes_toSharedSubstring(edge,-3);
			return out;
		}
		tmpPair.node1=IntervalGraph.nodesArray[edge.nodeID1];
		tmpPair.node2=IntervalGraph.nodesArray[edge.nodeID2];
		tmpPair.orientation=edge.orientation==0;
		tmpPair.diffs=edge.avgDiffs;
		if (edge.containment!=-1) {
			last=edge.projectOnto_containment(edge.nodeID1,tmpArray,-1);
			if (last!=-1) {
				tmpPair.alignmentStart1=tmpArray[0];
				tmpPair.alignmentEnd1=tmpArray[1];
				last=edge.projectOnto_containment(edge.nodeID2,tmpArray,-1);
				if (last!=-1) {
					tmpPair.alignmentStart2=tmpArray[0];
					tmpPair.alignmentEnd2=tmpArray[1];
					tmpPair.inCanonicalForm();
					if (inStep2_containment(tmpPair)==-1) out|=removeEdgeTypes_toSharedSubstring(edge,-1);
				}
				else out|=removeEdgeTypes_toSharedSubstring(edge,-1);
			}
			else out|=removeEdgeTypes_toSharedSubstring(edge,-1);
		}
		if (edge.insertion!=-1) {
			if (edge.insertion==Constants.INSERTION_BOTH) out|=removeEdgeTypes_toSharedSubstring(edge,-2);
			else {
				last=edge.projectOnto_insertion(edge.nodeID1,tmpArray,-1);
				if (last!=-1) {
					tmpPair.alignmentStart1=tmpArray[0];
					tmpPair.alignmentEnd1=tmpArray[1];
					last=edge.projectOnto_insertion(edge.nodeID2,tmpArray,-1);
					if (last!=-1) {
						tmpPair.alignmentStart2=tmpArray[0];
						tmpPair.alignmentEnd2=tmpArray[1];
						tmpPair.inCanonicalForm();
						if (inStep2_insertion(tmpPair)==-1) out|=removeEdgeTypes_toSharedSubstring(edge,-2);
					}
					else out|=removeEdgeTypes_toSharedSubstring(edge,-2);
				}
				else out|=removeEdgeTypes_toSharedSubstring(edge,-2);
			}
		}
		if (edge.overlap!=-1) {
			for (overlap=Constants.OVERLAP_PREFIX_PREFIX; overlap<=Constants.OVERLAP_SUFFIX_SUFFIX; overlap<<=1) {
				if (edge.overlap==-1) break;
				if ((edge.overlap&overlap)==0) continue;
				last=edge.projectOnto_overlap(edge.nodeID1,tmpArray,overlap,-1);
				if (last!=-1) {
					tmpPair.alignmentStart1=tmpArray[0];
					tmpPair.alignmentEnd1=tmpArray[1];
					last=edge.projectOnto_overlap(edge.nodeID2,tmpArray,overlap,-1);
					if (last!=-1) {	
						tmpPair.alignmentStart2=tmpArray[0];
						tmpPair.alignmentEnd2=tmpArray[1];
						tmpPair.inCanonicalForm();
						if (inStep2_overlap(tmpPair)==-1) out|=removeEdgeTypes_toSharedSubstring(edge,overlap);
					}
					else out|=removeEdgeTypes_toSharedSubstring(edge,overlap);
				}
				else out|=removeEdgeTypes_toSharedSubstring(edge,overlap);
			}
		}
		return out;
	}
	
	
	/**
	 * Converts to a shared substring type an alignment of type $type$ inside $edge$.
	 *
	 * @param type -1=containment; -2=insertion; -3=all types; >0=overlaps;
	 * @return TRUE iff $edge$ has been modified by the procedure.
	 */
	private static final boolean removeEdgeTypes_toSharedSubstring(IntervalGraph.Edge edge, int type) {
		boolean out;
		int id;
		
		out=false;
		if (edge.sharedSubstringOverhangs==null) {
			edge.sharedSubstringOverhangs = new int[4];
			Math.set(edge.sharedSubstringOverhangs,3,0);
		}
		if (type==-1) {
			if (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO) {
				edge.sharedSubstringOverhangs[0]=0;
				edge.sharedSubstringOverhangs[1]=0;
				edge.sharedSubstringOverhangs[2]=edge.containmentOverhangLeft;
				edge.sharedSubstringOverhangs[3]=edge.containmentOverhangRight;
				edge.containment=-1;
				edge.sharedSubstring=Constants.SHARED_SUBSTRING;
				out=true;
			}
			else if (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE) {
				edge.sharedSubstringOverhangs[0]=edge.containmentOverhangLeft;
				edge.sharedSubstringOverhangs[1]=edge.containmentOverhangRight;
				edge.sharedSubstringOverhangs[2]=0;
				edge.sharedSubstringOverhangs[3]=0;
				edge.containment=-1;
				edge.sharedSubstring=Constants.SHARED_SUBSTRING;
				out=true;
			}
		}
		else if (type==-2) {
			if (edge.insertion==Constants.INSERTION_ONE_IN_TWO) {
				edge.sharedSubstringOverhangs[0]=0;
				edge.sharedSubstringOverhangs[1]=0;
				edge.sharedSubstringOverhangs[2]=edge.containmentOverhangLeft;
				edge.sharedSubstringOverhangs[3]=edge.containmentOverhangRight;
				edge.insertion=-1;
				edge.sharedSubstring=Constants.SHARED_SUBSTRING;
				out=true;
			}
			else if (edge.insertion==Constants.INSERTION_TWO_IN_ONE) {
				edge.sharedSubstringOverhangs[0]=edge.containmentOverhangLeft;
				edge.sharedSubstringOverhangs[1]=edge.containmentOverhangRight;
				edge.sharedSubstringOverhangs[2]=0;
				edge.sharedSubstringOverhangs[3]=0;
				edge.insertion=-1;
				edge.sharedSubstring=Constants.SHARED_SUBSTRING;
				out=true;
			}
		}
		else if (type==Constants.OVERLAP_PREFIX_PREFIX) {
			id=Constants.overlap2id(Constants.OVERLAP_PREFIX_PREFIX);
			edge.sharedSubstringOverhangs[0]=0;
			edge.sharedSubstringOverhangs[1]=edge.overhangs[id][0];
			edge.sharedSubstringOverhangs[2]=0;
			edge.sharedSubstringOverhangs[3]=edge.overhangs[id][1];		
			edge.overlap&=~type;
			if (edge.overlap==0) edge.overlap=-1;
			edge.sharedSubstring=Constants.SHARED_SUBSTRING;
			out=true;
		}
		else if (type==Constants.OVERLAP_PREFIX_SUFFIX) {
			id=Constants.overlap2id(Constants.OVERLAP_PREFIX_SUFFIX);
			edge.sharedSubstringOverhangs[0]=0;
			edge.sharedSubstringOverhangs[1]=edge.overhangs[id][0];
			edge.sharedSubstringOverhangs[2]=edge.overhangs[id][1];
			edge.sharedSubstringOverhangs[3]=0;
			edge.overlap&=~type;
			if (edge.overlap==0) edge.overlap=-1;
			edge.sharedSubstring=Constants.SHARED_SUBSTRING;
			out=true;
		}
		else if (type==Constants.OVERLAP_SUFFIX_PREFIX) {
			id=Constants.overlap2id(Constants.OVERLAP_SUFFIX_PREFIX);
			edge.sharedSubstringOverhangs[0]=edge.overhangs[id][0];
			edge.sharedSubstringOverhangs[1]=0;
			edge.sharedSubstringOverhangs[2]=0;
			edge.sharedSubstringOverhangs[3]=edge.overhangs[id][1];
			edge.overlap&=~type;
			if (edge.overlap==0) edge.overlap=-1;
			edge.sharedSubstring=Constants.SHARED_SUBSTRING;
			out=true;
		}
		else if (type==Constants.OVERLAP_SUFFIX_SUFFIX) {
			id=Constants.overlap2id(Constants.OVERLAP_SUFFIX_SUFFIX);
			edge.sharedSubstringOverhangs[0]=edge.overhangs[id][0];
			edge.sharedSubstringOverhangs[1]=0;
			edge.sharedSubstringOverhangs[2]=edge.overhangs[id][1];
			edge.sharedSubstringOverhangs[3]=0;
			edge.overlap&=~type;
			if (edge.overlap==0) edge.overlap=-1;
			edge.sharedSubstring=Constants.SHARED_SUBSTRING;
			out=true;
		}
		else if (type==-3) {
			if (edge.containment!=-1) {
				removeEdgeTypes_toSharedSubstring(edge,-1);
				out=true;
			}
			if (edge.insertion!=-1) {
				removeEdgeTypes_toSharedSubstring(edge,-2);
				out=true;
			}
			if (edge.overlap!=-1 && (edge.overlap&Constants.OVERLAP_PREFIX_PREFIX)!=0) {
				removeEdgeTypes_toSharedSubstring(edge,Constants.OVERLAP_PREFIX_PREFIX);
				out=true;
			}
			if (edge.overlap!=-1 && (edge.overlap&Constants.OVERLAP_PREFIX_SUFFIX)!=0) {
				removeEdgeTypes_toSharedSubstring(edge,Constants.OVERLAP_PREFIX_SUFFIX);
				out=true;
			}
			if (edge.overlap!=-1 && (edge.overlap&Constants.OVERLAP_SUFFIX_PREFIX)!=0) {
				removeEdgeTypes_toSharedSubstring(edge,Constants.OVERLAP_SUFFIX_PREFIX);
				out=true;
			}
			if (edge.overlap!=-1 && (edge.overlap&Constants.OVERLAP_SUFFIX_SUFFIX)!=0) {
				removeEdgeTypes_toSharedSubstring(edge,Constants.OVERLAP_SUFFIX_SUFFIX);
				out=true;
			}
		}
		return out;
	}
	
	
	
	
	
	
	
		
	// ----------------------------- CONTAINMENT PROCEDURES ------------------------------
	
	public static final int isContained_alignment_alignment(IntervalGraph.AlignmentPair pair) {
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final int start1, end1, start2, end2;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		start1=pair.node1.start; end1=pair.node1.end;
		start2=pair.node2.start; end2=pair.node2.end;
		alignmentStart1=pair.alignmentStart1;
		alignmentEnd1=pair.alignmentEnd1-1;
		alignmentStart2=pair.alignmentStart2;
		alignmentEnd2=pair.alignmentEnd2-1;
		isLeftMaximal1=pair.node1.isLeftMaximal;
		isRightMaximal1=pair.node1.isRightMaximal;
		isLeftMaximal2=pair.node2.isLeftMaximal;
		isRightMaximal2=pair.node2.isRightMaximal;
		
		return isContained_alignment_alignment_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,isLeftMaximal1,isRightMaximal1,isLeftMaximal2,isRightMaximal2,pair.orientation);
	}
	
	
	private static final int isContained_alignment_alignment_impl(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean isLeftMaximal1, boolean isRightMaximal1, boolean isLeftMaximal2, boolean isRightMaximal2, boolean orientation) {
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		final boolean isLessThan21 = alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD;
		
		if (isLeftMaximal1 && isRightMaximal1) {
			if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
				 isLessThan21 
			   ) {
				if ( isLeftMaximal2 && !isRightMaximal2 &&
					 ( (orientation && Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD) ||
					   (!orientation && Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD)
					 )
				   ) return Constants.CONTAINMENT_TWO_IN_ONE;
				if ( isRightMaximal2 && !isLeftMaximal2 &&
					 ( (orientation && Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD) ||
					   (!orientation && Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD) 
					 ) 
				   ) return Constants.CONTAINMENT_TWO_IN_ONE;
				if (!isLeftMaximal2 && !isRightMaximal2) return Constants.CONTAINMENT_TWO_IN_ONE;
			}
		}
		else if (isLeftMaximal1) {
			if (isLeftMaximal2 && isRightMaximal2) {
				if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
					 ( (orientation && Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD) ||
				       (!orientation && Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD)
				     ) &&
				     Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     isLessThan12
			       ) return Constants.CONTAINMENT_ONE_IN_TWO;
			}
			else if (isLeftMaximal2 && !isRightMaximal2) {
				if ( orientation && 
					 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD
				   ) {
	   				if (Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
					if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
				}
				if (!orientation) {
					if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					     Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD
				       ) {
						if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
					}
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
						 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD 
					   ) {
						if (Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
					}
				}
			}
			else if (isRightMaximal2 && !isLeftMaximal2) {
				if (orientation) {
					if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD
				       ) {
	   					if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
					}
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
						 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD 
					   ) {
						if (Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
					}
				}
				if ( !orientation && 
					 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD
				   ) {
					if (Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
   					if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
				}
			}
			else {
				if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     isLessThan21
		           ) return Constants.CONTAINMENT_TWO_IN_ONE;
			}
		}
		else if (isRightMaximal1) {
			if (isLeftMaximal2 && isRightMaximal2) {
				if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
					 ( (orientation && Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD) ||
				       (!orientation && Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD)
				     ) &&
				     Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     isLessThan12
			       ) return Constants.CONTAINMENT_ONE_IN_TWO;
			}
			else if (isLeftMaximal2 && !isRightMaximal2) {
				if (orientation) {
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					     Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD
				       ) {
						if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
					}
					if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
						 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD 
					   ) {
						if (Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
					}
				}
				if ( !orientation && 
					 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD
				   ) {
	   				if (Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
					if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
				}
			}
			else if (isRightMaximal2 && !isLeftMaximal2) {
				if ( orientation && 
					 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD
				   ) {
					if (Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
   					if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
				}
				if (!orientation) {
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD
				       ) {
	   					if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
					}
					if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
						 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD 
					   ) {
						if (Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
					}
				}
			}
			else {
				if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     isLessThan21
		           ) return Constants.CONTAINMENT_TWO_IN_ONE;
			}
		}
		else {
			if (isLeftMaximal2) {
				if (orientation) {
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					     Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD
				       ) {
	   					if (Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
						if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
					}
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			     	     Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
				 		 isLessThan12
			   		   ) return Constants.CONTAINMENT_ONE_IN_TWO;
				}
				else {
					if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					     Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD
				       ) {
	   					if (Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
						if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
					}
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			     	     Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
				 		 isLessThan12
			   		   ) return Constants.CONTAINMENT_ONE_IN_TWO;
				}
			}
			else if (isRightMaximal2) {
				if (orientation) {
					if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD
				       ) {
	   					if (Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
						if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
					}
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			     	     Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
				 		 isLessThan12
			   		   ) return Constants.CONTAINMENT_ONE_IN_TWO;
				}
				else {
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD
				       ) {
	   					if (Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
						if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
					}
					if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			     	     Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
				 		 isLessThan12
			   		   ) return Constants.CONTAINMENT_ONE_IN_TWO;
				}
			}
			else {
				if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     isLessThan12 
				   ) return Constants.CONTAINMENT_ONE_IN_TWO;
				if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
				     isLessThan21 
				   ) return Constants.CONTAINMENT_TWO_IN_ONE;
			}
		}
		return -1;
	}
	
	
	/**
	 * @param prefixOrSuffix true=prefix, false=suffix.
	 */
	private static final int isContained_alignment_prefixOrSuffix(IntervalGraph.AlignmentPair pair, boolean prefixOrSuffix) {
		final boolean isLeftMaximal1, isRightMaximal1;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id1;
		if (pair.node1.type==Constants.INTERVAL_ALIGNMENT) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id1=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id1=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		
		return isContained_alignment_prefixOrSuffix_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,isLeftMaximal1,isRightMaximal1,id1,prefixOrSuffix,pair.orientation);
	}
	
	
	private static final int isContained_alignment_prefixOrSuffix_impl(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean isLeftMaximal1, boolean isRightMaximal1, int id1, boolean prefixOrSuffix, boolean orientation) {
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 isLessThan12 
		   ) {
			if (prefixOrSuffix) {
				if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD) return id1;
				if ( alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD && 
					 ( (orientation && !isLeftMaximal1) || (!orientation && !isRightMaximal1) )
				   ) return id1;
			}
			else {
				if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD) return id1;
				if ( alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
					 ( (orientation && !isRightMaximal1) || (!orientation && !isLeftMaximal1) ) 
				   ) return id1;
			}
		}
		return -1;
	}
	
	
	private static final int isContained_alignment_prefixSuffix(IntervalGraph.AlignmentPair pair) {
		final boolean isLeftMaximal1, isRightMaximal1;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id1;
		if (pair.node1.type==Constants.INTERVAL_ALIGNMENT) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id1=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id1=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		
		return isContained_alignment_prefixSuffix_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,isLeftMaximal1,isRightMaximal1,id1);
	}
	
	
	private static final int isContained_alignment_prefixSuffix_impl(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean isLeftMaximal1, boolean isRightMaximal1, int id1) {
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 isLessThan12 ) {
	 		if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD || 
	 			 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD
	 		   ) return id1;
	 		else if (!isLeftMaximal1 || !isRightMaximal1) return id1;
		}
		return -1;
	}
	
	
	private static final int isContained_alignment_substring(IntervalGraph.AlignmentPair pair) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id1;
		if (pair.node1.type==Constants.INTERVAL_ALIGNMENT) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id1=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id1=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		
		return isContained_alignment_substring_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,id1);
	}
	
	
	private static final int isContained_alignment_substring_impl(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, int id1) {
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 isLessThan12
		   ) return id1;
		return -1;
	}
	
	
	/**
	 * @param prefixOrSuffix true=prefix, false=suffix.
	 */
	private static final int isContained_prefix_prefixOrSuffix(IntervalGraph.AlignmentPair pair, boolean prefixOrSuffix) {
		final boolean isLeftMaximal1, isLeftMaximal2;
		final boolean isRightMaximal1, isRightMaximal2;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;  // In case of prefix-prefix, $alignment{Start,End}1$ refer to the first element of $pair$.
		final int id1, id2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_PREFIX) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isLeftMaximal2=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			isRightMaximal2=pair.node2.isRightMaximal;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id1=Constants.CONTAINMENT_ONE_IN_TWO; id2=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isLeftMaximal2=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			isRightMaximal2=pair.node1.isRightMaximal;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id1=Constants.CONTAINMENT_TWO_IN_ONE; id2=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		final boolean isLessThan21 = alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD;
		
		if (prefixOrSuffix) {
			if (pair.orientation) { 
				if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD
				   ) {
					 if (Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return id1;
					 if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return id2;
				}
				if ( alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 isLessThan12 &&
					 !isLeftMaximal1
				   ) return id1;
				if ( alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 isLessThan21 &&
					 !isLeftMaximal2
				   ) return id2;
			}
		}
		else {
			if (!pair.orientation) {
				if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD 
				   ) {
					 if (Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return id1;
					 if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return id2;
				}
				if ( alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 isLessThan12 &&
					 !isLeftMaximal1
				   ) return id1;
				if ( alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
					 isLessThan21 &&
					 !isRightMaximal2
				   ) return id2;
			}
		}
		return -1;
	}
	
	
	private static final int isContained_prefix_prefixSuffix(IntervalGraph.AlignmentPair pair) {
		final boolean isLeftMaximal1;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id1;
		if (pair.node1.type==Constants.INTERVAL_DENSE_PREFIX) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id1=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id1=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 isLessThan12
		   ) {
			if (pair.orientation) {
				if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD) return id1;
				if (alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD && !isLeftMaximal1) return id1;
			}
			else {
				if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD) return id1;
				if (alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD && !isLeftMaximal1) return id1;
			}
		}
		return -1;
	}
	
	
	private static final int isContained_suffix_suffix(IntervalGraph.AlignmentPair pair) {
		final boolean isRightMaximal1, isRightMaximal2;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;  // $alignment{Start,End}1$ refers to the first element in $pair$.
		start1=pair.node1.start; start2=pair.node2.start;
		end1=pair.node1.end; end2=pair.node2.end;
		isRightMaximal1=pair.node1.isRightMaximal;
		isRightMaximal2=pair.node2.isRightMaximal;
		alignmentStart1=pair.alignmentStart1;
		alignmentEnd1=pair.alignmentEnd1-1;
		alignmentStart2=pair.alignmentStart2;
		alignmentEnd2=pair.alignmentEnd2-1;
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		final boolean isLessThan21 = alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD;
		
		if (pair.orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD 
			   ) {
		    	 if (Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan12) return Constants.CONTAINMENT_ONE_IN_TWO;
				 if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && isLessThan21) return Constants.CONTAINMENT_TWO_IN_ONE;
			}
			if ( alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 isLessThan12 &&	 
				 !isRightMaximal1
			   ) return Constants.CONTAINMENT_ONE_IN_TWO;
			if ( alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 isLessThan21 &&
				 !isRightMaximal2
			   ) return Constants.CONTAINMENT_TWO_IN_ONE;
		}
		return -1;
	}
	
	
	private static final int isContained_suffix_prefixSuffix(IntervalGraph.AlignmentPair pair) {
		final boolean isRightMaximal1;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id1;
		if (pair.node1.type==Constants.INTERVAL_DENSE_SUFFIX) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			isRightMaximal1=pair.node1.isRightMaximal;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id1=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			isRightMaximal1=pair.node2.isRightMaximal;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id1=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
			 isLessThan12
		   ) {
			if (pair.orientation) {
				if (Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD) return id1;
				if (alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD && !isRightMaximal1) return id1;
			}
			else {
				if (Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD) return id1;
				if (alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD && !isRightMaximal1) return id1;
			}
		}
		return -1;
	}
	
	
	private static final int isContained_prefixSuffix_prefixSuffix(IntervalGraph.AlignmentPair pair) {
		final boolean isLeftMaximal1, isRightMaximal1, isLeftMaximal2, isRightMaximal2;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;  // $alignment{Start,End}1$ refer to the first element in $pair$.
		start1=pair.node1.start; start2=pair.node2.start;
		end1=pair.node1.end; end2=pair.node2.end;
		isLeftMaximal1=pair.node1.isLeftMaximal;
		isRightMaximal1=pair.node1.isRightMaximal;
		isLeftMaximal2=pair.node2.isLeftMaximal;
		isRightMaximal2=pair.node2.isRightMaximal;
		alignmentStart1=pair.alignmentStart1;
		alignmentEnd1=pair.alignmentEnd1-1;
		alignmentStart2=pair.alignmentStart2;
		alignmentEnd2=pair.alignmentEnd2-1;
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		final boolean isLessThan21 = alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 isLessThan12 ) {
			if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD ||
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD 
			   ) return Constants.CONTAINMENT_ONE_IN_TWO;
			if (!isLeftMaximal1 && !isRightMaximal1) return Constants.CONTAINMENT_ONE_IN_TWO;
		}
		if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 isLessThan21 ) {
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD ||
				 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD 
			   ) return Constants.CONTAINMENT_TWO_IN_ONE;
			if (!isLeftMaximal2 && !isRightMaximal2) return Constants.CONTAINMENT_TWO_IN_ONE;
		}
		return -1;
	}
	
	
	private static final int isContained_prefixSuffix_singleDeletion(IntervalGraph.AlignmentPair pair) {
		final boolean isLeftMaximal1, isRightMaximal1;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id1;
		if (pair.node1.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id1=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id1=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 isLessThan12 ) {
			if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD ||
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD 
			   ) return id1;
			if (!isLeftMaximal1 && !isRightMaximal1) return id1;
		}
		return -1;
	}
	
	
	private static final int isContained_substring_prefixOrSuffix(IntervalGraph.AlignmentPair pair) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_SUBSTRING) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id2=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id2=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		final boolean isLessThan21 = alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 isLessThan21
		   ) return id2;
		return -1;
	}
	
	
	private static final int isContained_substring_substring(IntervalGraph.AlignmentPair pair) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		start1=pair.node1.start; start2=pair.node2.start;
		end1=pair.node1.end; end2=pair.node2.end;
		alignmentStart1=pair.alignmentStart1;
		alignmentEnd1=pair.alignmentEnd1-1;
		alignmentStart2=pair.alignmentStart2;
		alignmentEnd2=pair.alignmentEnd2-1;
		final boolean isLessThan12 = alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD;
		final boolean isLessThan21 = alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 isLessThan12
		   ) return Constants.CONTAINMENT_ONE_IN_TWO;
		if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 isLessThan21
		   ) return Constants.CONTAINMENT_TWO_IN_ONE;
		return -1;
	}
	
	
	private static final int isContained_substring_singleDeletion(IntervalGraph.AlignmentPair pair) {
		final int start1, end1, start2, end2, alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_SUBSTRING) {
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id2=Constants.CONTAINMENT_TWO_IN_ONE;
		}
		else {
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id2=Constants.CONTAINMENT_ONE_IN_TWO;
		}
		final boolean isLessThan21 = alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD;
		
		if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 isLessThan21
		   ) return id2;
		return -1;
	}
	
	
	/**
	 * @param *1 related to the non-alignment node;
	 * @param *2 related to the alignment node;
	 * @return true iff the alignment between nodes $nonAlignmentNodeID$ and 
	 * $alignmentNodeID$ can be of containment type.
	 */
	public static final boolean isContained_alignment(int nonAlignmentNodeID, int alignmentNodeID, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean isLeftMaximal1, boolean isRightMaximal1, boolean isLeftMaximal2, boolean isRightMaximal2, boolean orientation) {
		IntervalGraph.Node nonAlignmentNode = IntervalGraph.nodesArray[nonAlignmentNodeID];
		IntervalGraph.Node alignmentNode = IntervalGraph.nodesArray[alignmentNodeID];
		
		switch (IntervalGraph.nodesArray[nonAlignmentNodeID].type) {
			case Constants.INTERVAL_ALIGNMENT: return isContained_alignment_alignment_impl(alignmentNode.start,alignmentNode.end,nonAlignmentNode.start,nonAlignmentNode.end,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,isLeftMaximal2,isRightMaximal2,isLeftMaximal1,isRightMaximal1,orientation)!=-1;
			case Constants.INTERVAL_DENSE_PREFIX: return isContained_alignment_prefixOrSuffix_impl(alignmentNode.start,alignmentNode.end,nonAlignmentNode.start,nonAlignmentNode.end,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,isLeftMaximal2,isRightMaximal2,/*irrelevant*/Constants.CONTAINMENT_ONE_IN_TWO,true,orientation)!=-1;
			case Constants.INTERVAL_DENSE_SUFFIX: return isContained_alignment_prefixOrSuffix_impl(alignmentNode.start,alignmentNode.end,nonAlignmentNode.start,nonAlignmentNode.end,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,isLeftMaximal2,isRightMaximal2,/*irrelevant*/Constants.CONTAINMENT_ONE_IN_TWO,false,orientation)!=-1;
			case Constants.INTERVAL_DENSE_PREFIXSUFFIX: return isContained_alignment_prefixSuffix_impl(alignmentNode.start,alignmentNode.end,nonAlignmentNode.start,nonAlignmentNode.end,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,isLeftMaximal2,isRightMaximal2,/*irrelevant*/Constants.CONTAINMENT_ONE_IN_TWO)!=-1;
			case Constants.INTERVAL_DENSE_SUBSTRING: return isContained_alignment_substring_impl(alignmentNode.start,alignmentNode.end,nonAlignmentNode.start,nonAlignmentNode.end,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,/*irrelevant*/Constants.CONTAINMENT_ONE_IN_TWO)!=-1;
			case Constants.INTERVAL_DENSE_SINGLEDELETION: return isContained_alignment_prefixSuffix_impl(alignmentNode.start,alignmentNode.end,nonAlignmentNode.start,nonAlignmentNode.end,alignmentStart2,alignmentEnd2,alignmentStart1,alignmentEnd1,isLeftMaximal2,isRightMaximal2,/*irrelevant*/Constants.CONTAINMENT_ONE_IN_TWO)!=-1;
		}
		return false;
	}
	
	
	/**
	 * Alignments that are shorter than the period of a long-period periodic interval 
	 * might be induced by nested repeats inside the long period. Alignments that are 
	 * too long compared to the period imply that the other interval is long enough to 
	 * have been marked as periodic.
	 *
	 * Remark: only alignment intervals can be contained in a long-period interval, since
	 * other types of interval could be induced by nested repeats that are longer than a
	 * period.
	 *
	 * Remark: a short-period periodic interval can be a substring of a transposon, and a 
	 * diverged substring of a short-period periodic interval might become an independent 
	 * repeat that aligns just to a substring of the original periodic interval. Such
	 * relationships are not marked as containment.
	 */
	private static final int isContained_periodic_alignment(IntervalGraph.AlignmentPair pair) {
		final double THRESHOLD = 3.0;
		final boolean orientation, isLeftMaximal2, isRightMaximal2;
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentEnd1, alignmentStart2, alignmentEnd2;
		final int id, length2, period1;
		orientation=pair.orientation;
		if (pair.node1.type==Constants.INTERVAL_PERIODIC) {
			if (!pair.node1.hasLongPeriod) return -1;
			period1=pair.node1.period;
			if (period1<=0) return -1;
			start1=pair.node1.start; start2=pair.node2.start;
			end1=pair.node1.end; end2=pair.node2.end;
			isLeftMaximal2=pair.node2.isLeftMaximal;
			isRightMaximal2=pair.node2.isRightMaximal;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			id=Constants.CONTAINMENT_TWO_IN_ONE;
			length2=pair.node2.length();
		}
		else {
			if (!pair.node2.hasLongPeriod) return -1;
			period1=pair.node2.period;
			if (period1<=0) return -1;
			start1=pair.node2.start; start2=pair.node1.start;
			end1=pair.node2.end; end2=pair.node1.end;
			isLeftMaximal2=pair.node1.isLeftMaximal;
			isRightMaximal2=pair.node1.isRightMaximal;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			id=Constants.CONTAINMENT_ONE_IN_TWO;
			length2=pair.node1.length();
		}
		final boolean isLessThan21 = alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD;		
		
		if ( length2<=period1+IntervalGraph.IDENTITY_THRESHOLD || 
			 length2>=THRESHOLD*period1 || 
		     Math.abs(alignmentStart2,start2)>IntervalGraph.IDENTITY_THRESHOLD || 
			 Math.abs(alignmentEnd2,end2)>IntervalGraph.IDENTITY_THRESHOLD ||
			 !isLessThan21
		   ) return -1;
		if ((orientation?isLeftMaximal2:isRightMaximal2) && alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD) return -1;
		if ((orientation?isRightMaximal2:isLeftMaximal2) && alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD) return -1;
		return id;
	}
	
	
	private static final int isContained_periodic_periodic(IntervalGraph.AlignmentPair pair) {
		if (pair.node1.hasLongPeriod!=pair.node2.hasLongPeriod) return -1;
		
		return isContained_alignment_alignment(pair);
	}
	
	
	
	
	// -------------------------------- OVERLAP PROCEDURES -------------------------------
	
	/**
	 * @return the symmetrized version of $string$ if $id1=2$ and $id2=1$, $string$ itself
	 * otherwise.
	 */
	private static final int orient(int string, int id1, int id2) {
		if (id1==1) return string;
		return Constants.oppositeDirectionOverlap(string);
	}
	
	
	/**
	 * Checks just alignment start/end positions, does not use maximality.
	 */
	public static final int isOverlap_simple(IntervalGraph.AlignmentPair pair) {
		return isOverlap_simple_impl(pair.node1.start,pair.node1.end,pair.node2.start,pair.node2.end,pair.alignmentStart1,pair.alignmentEnd1-1,pair.alignmentStart2,pair.alignmentEnd2-1,pair.orientation,IntervalGraph.IDENTITY_THRESHOLD);
	}
	
	
	public static final int isOverlap_simple_impl(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean orientation, int threshold) {
		if (orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=threshold &&
				 alignmentStart1>start1+threshold &&
				 Math.abs(alignmentStart2,start2)<=threshold &&
				 alignmentEnd2<end2-threshold
			   ) return Constants.OVERLAP_SUFFIX_PREFIX;
			if ( Math.abs(alignmentStart1,start1)<=threshold &&
				 alignmentEnd1<end1-threshold &&
				 Math.abs(alignmentEnd2,end2)<=threshold &&
				 alignmentStart2>start2+threshold
			   ) return Constants.OVERLAP_PREFIX_SUFFIX;
		}
		else {
			if ( Math.abs(alignmentEnd1,end1)<=threshold &&
				 alignmentStart1>start1+threshold &&
				 Math.abs(alignmentEnd2,end2)<=threshold &&
				 alignmentStart2>start2+threshold
			   ) return Constants.OVERLAP_SUFFIX_SUFFIX;
			if ( Math.abs(alignmentStart1,start1)<=threshold &&
				 alignmentEnd1<end1-threshold &&
				 Math.abs(alignmentStart2,start2)<=threshold &&
				 alignmentEnd2<end2-threshold
			   ) return Constants.OVERLAP_PREFIX_PREFIX;
		}
		return -1;
	}
	
	
	public static final int isOverlap_alignment_alignment(IntervalGraph.AlignmentPair pair) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		start1=pair.node1.start; end1=pair.node1.end;
		start2=pair.node2.start; end2=pair.node2.end;
		alignmentStart1=pair.alignmentStart1;
		alignmentEnd1=pair.alignmentEnd1-1;
		alignmentStart2=pair.alignmentStart2;
		alignmentEnd2=pair.alignmentEnd2-1;
		isLeftMaximal1=pair.node1.isLeftMaximal;
		isRightMaximal1=pair.node1.isRightMaximal;
		isLeftMaximal2=pair.node2.isLeftMaximal;
		isRightMaximal2=pair.node2.isRightMaximal;
		
		return isOverlap_alignment_alignment_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,isLeftMaximal1,isRightMaximal1,isLeftMaximal2,isRightMaximal2,pair.orientation);
	}
	
	
	private static final int isOverlap_alignment_alignment_impl(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean isLeftMaximal1, boolean isRightMaximal1, boolean isLeftMaximal2, boolean isRightMaximal2, boolean orientation) {
		if (orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal1 && !isLeftMaximal2
			   ) return Constants.OVERLAP_SUFFIX_PREFIX;
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal1 && !isRightMaximal2
			   ) return Constants.OVERLAP_PREFIX_SUFFIX;
		}
		else {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal1 && !isRightMaximal2
			   ) return Constants.OVERLAP_SUFFIX_SUFFIX;
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal1 && !isLeftMaximal2
			   ) return Constants.OVERLAP_PREFIX_PREFIX;
		}
		return -1;
	}
	
	
	/**
	 * Remark: for a suffix-prefix overlap between prefix substring V and prefix substring
	 * W, we don't require V to be non-right-maximal, since V could be a prefix of a 
	 * master prefix substring U, and W could be a prefix of U as well but with a low
	 * quality region to its left.
	 */
	private static final int isOverlap_prefix_prefixOrSuffix(IntervalGraph.AlignmentPair pair, boolean prefixOrSuffix) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final int minPrefixLength, minSuffixLength;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		final int length1, length2, alignmentLength1, alignmentLength2;
		final int id1, id2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_PREFIX) {
			start1=pair.node1.start; end1=pair.node1.end;
			start2=pair.node2.start; end2=pair.node2.end;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			isLeftMaximal2=pair.node2.isLeftMaximal;
			isRightMaximal2=pair.node2.isRightMaximal;
			minPrefixLength=pair.node1.minPrefixLength;
			minSuffixLength=pair.node2.minSuffixLength;
			length1=pair.node1.length();
			length2=pair.node2.length();
			alignmentLength1=pair.alignmentLength1();
			alignmentLength2=pair.alignmentLength2();
			id1=1; id2=2;
		}
		else {
			start1=pair.node2.start; end1=pair.node2.end;
			start2=pair.node1.start; end2=pair.node1.end;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			isLeftMaximal2=pair.node1.isLeftMaximal;
			isRightMaximal2=pair.node1.isRightMaximal;
			minPrefixLength=pair.node2.minPrefixLength;
			minSuffixLength=pair.node1.minSuffixLength;
			length1=pair.node2.length();
			length2=pair.node1.length();
			alignmentLength1=pair.alignmentLength2();
			alignmentLength2=pair.alignmentLength1();
			id1=2; id2=1;
		}
		
		if (prefixOrSuffix && pair.orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal2
			   ) return orient(Constants.OVERLAP_SUFFIX_PREFIX,id1,id2);
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal1
			   ) return orient(Constants.OVERLAP_PREFIX_SUFFIX,id1,id2);
		}
		if (!prefixOrSuffix && !pair.orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal2
			   ) return orient(Constants.OVERLAP_SUFFIX_SUFFIX,id1,id2);
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal1
			   ) return orient(Constants.OVERLAP_PREFIX_PREFIX,id1,id2);
		}
		if (!prefixOrSuffix && pair.orientation) {
			if ( alignmentLength1<length1-IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentLength2<length2-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 (alignmentEnd1>=start1+minPrefixLength-IntervalGraph.IDENTITY_THRESHOLD || !isRightMaximal2) &&
				 (alignmentStart2<=end2-minSuffixLength+IntervalGraph.IDENTITY_THRESHOLD || !isLeftMaximal1)
			   ) return orient(Constants.OVERLAP_PREFIX_SUFFIX,id1,id2);
		}
		return -1;
	}

	
	/**
	 * Remark: a symmetrical observation to the one in $isOverlap_prefix_prefixOrSuffix$
	 * applies.
	 */
	private static final int isOverlap_suffix_suffix(IntervalGraph.AlignmentPair pair) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		start1=pair.node1.start; end1=pair.node1.end;
		start2=pair.node2.start; end2=pair.node2.end;
		alignmentStart1=pair.alignmentStart1;
		alignmentEnd1=pair.alignmentEnd1-1;
		alignmentStart2=pair.alignmentStart2;
		alignmentEnd2=pair.alignmentEnd2-1;
		isLeftMaximal1=pair.node1.isLeftMaximal;
		isRightMaximal1=pair.node1.isRightMaximal;
		isLeftMaximal2=pair.node2.isLeftMaximal;
		isRightMaximal2=pair.node2.isRightMaximal;
		
		if (pair.orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal1
			   ) return Constants.OVERLAP_SUFFIX_PREFIX;
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal2
			   ) return Constants.OVERLAP_PREFIX_SUFFIX;
		}
		return -1;
	}
	
	
	/**
	 * Remark: for a suffix-prefix overlap between V of substring type and W of prefix 
	 * type, we don't require V to be non-right-maximal and W to be non-left-maximal, 
	 * since V could be a substring of a master substring type U, and W could be a suffix 
	 * of U.
	 */
	private static final int isOverlap_substring_prefixOrSuffix(IntervalGraph.AlignmentPair pair, boolean prefixOrSuffix) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2, minLength, length2;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		final int id1, id2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_SUBSTRING) {
			start1=pair.node1.start; end1=pair.node1.end;
			start2=pair.node2.start; end2=pair.node2.end;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			isLeftMaximal2=pair.node2.isLeftMaximal;
			isRightMaximal2=pair.node2.isRightMaximal;
			id1=1; id2=2;
			minLength=pair.node1.minPrefixLength;
			length2=pair.node2.length();
		}
		else {
			start1=pair.node2.start; end1=pair.node2.end;
			start2=pair.node1.start; end2=pair.node1.end;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			isLeftMaximal2=pair.node1.isLeftMaximal;
			isRightMaximal2=pair.node1.isRightMaximal;
			id1=2; id2=1;
			minLength=pair.node2.minPrefixLength;
			length2=pair.node1.length();
		}
		
		if (length2>minLength+ALIGNMENT_OVERLAP_THRESHOLD) return -1;
		if ( prefixOrSuffix &&
		     Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
		     alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD
		   ) {
			if ( pair.orientation &&
				 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD
			   ) return orient(Constants.OVERLAP_SUFFIX_PREFIX,id1,id2);
			if ( !pair.orientation &&
				 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD
			   ) return orient(Constants.OVERLAP_PREFIX_PREFIX,id1,id2);
		}
		if ( !prefixOrSuffix &&
		     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
		     alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD
		   ) {
			if ( pair.orientation &&
				 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD
			   ) return orient(Constants.OVERLAP_PREFIX_SUFFIX,id1,id2);
			if ( !pair.orientation &&
				 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD
			   ) return orient(Constants.OVERLAP_SUFFIX_SUFFIX,id1,id2);
		}
		return -1;
	}
	
	
	/**
	 * Remark: the same observations as in $isOverlap_prefix_prefixOrSuffix$ and
	 * $isOverlap_suffix_suffix$ apply.
	 */
	private static final int isOverlap_prefixSuffix_prefixOrSuffix(IntervalGraph.AlignmentPair pair, boolean prefixOrSuffix) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		final int id1, id2;
		final int minPrefixLength, minSuffixLength;
		if (pair.node1.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX || pair.node1.type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
			start1=pair.node1.start; end1=pair.node1.end;
			start2=pair.node2.start; end2=pair.node2.end;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			isLeftMaximal2=pair.node2.isLeftMaximal;
			isRightMaximal2=pair.node2.isRightMaximal;
			id1=1; id2=2;
			minPrefixLength=pair.node1.minPrefixLength;
			minSuffixLength=pair.node1.minSuffixLength;
		}
		else {
			start1=pair.node2.start; end1=pair.node2.end;
			start2=pair.node1.start; end2=pair.node1.end;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			isLeftMaximal2=pair.node1.isLeftMaximal;
			isRightMaximal2=pair.node1.isRightMaximal;
			id1=2; id2=1;
			minPrefixLength=pair.node2.minPrefixLength;
			minSuffixLength=pair.node2.minSuffixLength;
		}
		
		if ( prefixOrSuffix &&
		     Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
		     alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD
		   ) {
			if ( pair.orientation &&
				 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 (alignmentStart1>=end1-minSuffixLength-ALIGNMENT_OVERLAP_THRESHOLD)
			   ) return orient(Constants.OVERLAP_SUFFIX_PREFIX,id1,id2);
			if ( !pair.orientation &&
				 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 (alignmentEnd1<=start1+minPrefixLength+ALIGNMENT_OVERLAP_THRESHOLD)
			   ) return orient(Constants.OVERLAP_PREFIX_PREFIX,id1,id2);
		}
		if ( !prefixOrSuffix &&
		     Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
		     alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD
		   ) {
			if ( pair.orientation &&
				 Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 (alignmentEnd1<=start1+minPrefixLength+ALIGNMENT_OVERLAP_THRESHOLD)
			   ) return orient(Constants.OVERLAP_PREFIX_SUFFIX,id1,id2);
			if ( !pair.orientation &&
				 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 (alignmentStart1>=end1-minSuffixLength-ALIGNMENT_OVERLAP_THRESHOLD)
			   ) return orient(Constants.OVERLAP_SUFFIX_SUFFIX,id1,id2);
		}
		return -1;
	}
	
	
	/**
	 * Remark: if the alignment interval is long enough, it should exhibit dense behavior,
	 * unless it corresponds to the non-prefix part of the prefix interval.
	 */
	private static final int isOverlap_prefix_alignment(IntervalGraph.AlignmentPair pair) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final int minPrefixLength;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		final int id1, id2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_PREFIX) {
			start1=pair.node1.start; end1=pair.node1.end;
			start2=pair.node2.start; end2=pair.node2.end;
			minPrefixLength=pair.node1.minPrefixLength;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			isLeftMaximal2=pair.node2.isLeftMaximal;
			isRightMaximal2=pair.node2.isRightMaximal;
			id1=1; id2=2;
		}
		else {
			start1=pair.node2.start; end1=pair.node2.end;
			start2=pair.node1.start; end2=pair.node1.end;
			minPrefixLength=pair.node2.minPrefixLength;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			isLeftMaximal2=pair.node1.isLeftMaximal;
			isRightMaximal2=pair.node1.isRightMaximal;
			id1=2; id2=1;
		}
		
		return isOverlap_prefix_alignment_impl(start1,end1,start2,end2,minPrefixLength,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,isLeftMaximal1,isRightMaximal1,isLeftMaximal2,isRightMaximal2,id1,id2,pair.orientation);
	}
	
	
	private static final int isOverlap_prefix_alignment_impl(int start1, int end1, int start2, int end2, int minPrefixLength, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean isLeftMaximal1, boolean isRightMaximal1, boolean isLeftMaximal2, boolean isRightMaximal2, int id1, int id2, boolean orientation) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		
		if (orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal2 &&
				 length2<=Alignments.minAlignmentLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_SUFFIX_PREFIX,id1,id2);
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal1 &&
				 alignmentEnd1<=start1+minPrefixLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_PREFIX_SUFFIX,id1,id2);
		}
		else {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal2 &&
				 length2<=Alignments.minAlignmentLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_SUFFIX_SUFFIX,id1,id2);
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal1 &&
				 alignmentEnd1<=start1+minPrefixLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_PREFIX_PREFIX,id1,id2);
		}
		return -1;
	}
	
	
	/**
	 * Symmetric to $isOverlap_prefix_alignment()$.
	 */
	private static final int isOverlap_suffix_alignment(IntervalGraph.AlignmentPair pair) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final int minSuffixLength;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		final int id1, id2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_SUFFIX) {
			start1=pair.node1.start; end1=pair.node1.end;
			start2=pair.node2.start; end2=pair.node2.end;
			minSuffixLength=pair.node1.minSuffixLength;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			isLeftMaximal1=pair.node1.isLeftMaximal;
			isRightMaximal1=pair.node1.isRightMaximal;
			isLeftMaximal2=pair.node2.isLeftMaximal;
			isRightMaximal2=pair.node2.isRightMaximal;
			id1=1; id2=2;
		}
		else {
			start1=pair.node2.start; end1=pair.node2.end;
			start2=pair.node1.start; end2=pair.node1.end;
			minSuffixLength=pair.node2.minSuffixLength;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			isLeftMaximal1=pair.node2.isLeftMaximal;
			isRightMaximal1=pair.node2.isRightMaximal;
			isLeftMaximal2=pair.node1.isLeftMaximal;
			isRightMaximal2=pair.node1.isRightMaximal;
			id1=2; id2=1;
		}
		
		return isOverlap_suffix_alignment_impl(start1,end1,start2,end2,minSuffixLength,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,isLeftMaximal1,isRightMaximal1,isLeftMaximal2,isRightMaximal2,id1,id2,pair.orientation);
	}
	
	
	private static final int isOverlap_suffix_alignment_impl(int start1, int end1, int start2, int end2, int minSuffixLength, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, boolean isLeftMaximal1, boolean isRightMaximal1, boolean isLeftMaximal2, boolean isRightMaximal2, int id1, int id2, boolean orientation) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		
		if (orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal1 &&
				 alignmentStart1>=end1-minSuffixLength-ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_SUFFIX_PREFIX,id1,id2);
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal2 &&
				 length2<=Alignments.minAlignmentLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_PREFIX_SUFFIX,id1,id2);
		}
		else {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 !isRightMaximal1 &&
				 alignmentStart1>=end1-minSuffixLength-ALIGNMENT_OVERLAP_THRESHOLD	 
			   ) return orient(Constants.OVERLAP_SUFFIX_SUFFIX,id1,id2);
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 !isLeftMaximal2 &&
				 length2<=Alignments.minAlignmentLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_PREFIX_PREFIX,id1,id2);
		}
		return -1;
	}
	
	
	/**
	 * Remark: if the alignment interval is long enough, it should exhibit dense behavior.
	 */
	private static final int isOverlap_prefixSuffix_alignment(IntervalGraph.AlignmentPair pair) {
		final int start1, start2, end1, end2;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final int id1, id2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
			start1=pair.node1.start; end1=pair.node1.end;
			start2=pair.node2.start; end2=pair.node2.end;
			alignmentStart1=pair.alignmentStart1;
			alignmentEnd1=pair.alignmentEnd1-1;
			alignmentStart2=pair.alignmentStart2;
			alignmentEnd2=pair.alignmentEnd2-1;
			id1=1; id2=2;
		}
		else {
			start1=pair.node2.start; end1=pair.node2.end;
			start2=pair.node1.start; end2=pair.node1.end;
			alignmentStart1=pair.alignmentStart2;
			alignmentEnd1=pair.alignmentEnd2-1;
			alignmentStart2=pair.alignmentStart1;
			alignmentEnd2=pair.alignmentEnd1-1;
			id1=2; id2=1;
		}
		
		return isOverlap_prefixSuffix_alignment_impl(start1,end1,start2,end2,alignmentStart1,alignmentEnd1,alignmentStart2,alignmentEnd2,id1,id2,pair.orientation);
	}

	
	private static final int isOverlap_prefixSuffix_alignment_impl(int start1, int end1, int start2, int end2, int alignmentStart1, int alignmentEnd1, int alignmentStart2, int alignmentEnd2, int id1, int id2, boolean orientation) {
		final int length1 = end1-start1+1;
		final int length2 = end2-start2+1;
		
		if (orientation) {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 length2<=Alignments.minAlignmentLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_SUFFIX_PREFIX,id1,id2);
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 length2<=Alignments.minAlignmentLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_PREFIX_SUFFIX,id1,id2);
		}
		else {
			if ( Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD &&
				 length2<=Alignments.minAlignmentLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_SUFFIX_SUFFIX,id1,id2);
			if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD &&
				 Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD &&
				 alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD &&
				 length2<=Alignments.minAlignmentLength+ALIGNMENT_OVERLAP_THRESHOLD
			   ) return orient(Constants.OVERLAP_PREFIX_PREFIX,id1,id2);
		}
		return -1;
	}
	
	
	private static final int isOverlap_substring_prefixSuffix(IntervalGraph.AlignmentPair pair) {
		final int minLength, length2;
		if (pair.node1.type==Constants.INTERVAL_DENSE_SUBSTRING) {
			minLength=pair.node1.minPrefixLength;
			length2=pair.node2.length();
		}
		else {
			minLength=pair.node2.minPrefixLength;
			length2=pair.node1.length();
		}
		
		if (length2>minLength+ALIGNMENT_OVERLAP_THRESHOLD) return -1;
		return isOverlap_simple(pair);
	}
	
	
	/**
	 * We use the weak property just with dense substrings of substring type, since it 
	 * does not seem necessary with the other types.
	 */
	private static final int isOverlap_substring_substring(IntervalGraph.AlignmentPair pair) {
		if (pair.node1.isWeak!=pair.node2.isWeak) return -1;
		return isOverlap_simple(pair);
	}
	
	
	/**
	 * Only long period allowed. The alignment interval must not be too long, otherwise it
	 * should have been classified as long period.
	 */
	private static final int isOverlap_periodic_alignment(IntervalGraph.AlignmentPair pair) {
		final double THRESHOLD = 3.0;
		final int length2, period1;
		if (pair.node1.type==Constants.INTERVAL_PERIODIC) {
			if (!pair.node1.hasLongPeriod) return -1;
			period1=pair.node1.period;
			if (period1<=0) return -1;
			length2=pair.node2.length();
		}
		else {
			if (!pair.node2.hasLongPeriod) return -1;
			period1=pair.node2.period;
			if (period1<=0) return -1;
			length2=pair.node1.length();
		}
		
		if (length2>=THRESHOLD*period1) return -1;
		return isOverlap_alignment_alignment(pair);
	}
	
	
	private static final int isOverlap_periodic_periodic(IntervalGraph.AlignmentPair pair) {
		if (pair.node1.hasLongPeriod!=pair.node2.hasLongPeriod) return -1;
		
		return isOverlap_alignment_alignment(pair);
	}

	
	/**
	 * Let $X$ be a repeat. It could happen that a full occurrence of $X$ is subject 
	 * to iterated insertions of other repeats, breaking $X$ into intervals that are
	 * both left- and right-maximal. Such insertions usually occur at random 
	 * positions. If many occurrences of $X$ are subject to insertions, then an
	 * occurrence of $X$ that is not subject to insertions should look like a dense 
	 * substring of substring type, and its fragments should be assigned to the 
	 * correct connected component by procedure $inStep2_containment$. This procedure 
	 * marks the cases in which the occurrence of $X$ without insertions is not 
	 * detected as a dense substring of substring type: in this case the alignments
	 * might be ignored during the factorization of the occurrence of $X$ without 
	 * insertions, and they could be eventually assigned to the full occurrence of $X$
	 * without insertions.
	 *
	 * Remark: the same happens when $X$ is subject to multiple internal deletions: in 
	 * this case the left- and right-maximal fragments are adjacent in the copy with 
	 * deletions. Having multiple internal deletions is very common for transposable 
	 * elements: see e.g. \cite{quesneville2003detection}.
	 *
	 * Remark: this is the first step of the pipeline in which we explicitly model random
	 * insertions of repeats inside others.
	 *
	 * Remark: we also support intervals that are not left-maximal or not 
	 * right-maximal.
	 *
	 * Remark: allowing in Step 2 an edge between nodes that could be related by 
	 * repeat insertion might put distinct repeats into the same connected component. 
	 * Assume e.g. that two distinct repeats $X$ and $Y$ contain a similar substring 
	 * $Z$, which also occurs by itself. Such substring should be assigned a distinct 
	 * interval in every occurrence of $X$ and $Y$. However, if this doesn't happen in
	 * just one occurrence of $X$ and in just one occurrence of $Y$, then the Step 2 
	 * graph contains an edge that connects an interval of $X$ to an interval of $Z$, 
	 * and an interval of $Z$ to an interval of $Y$.
	 *
	 * Remark: the problem described above is very likely to happen with dense 
	 * substrings, since detecting contained intervals is nontrivial.
	 *
	 * @return a code $INSERTION_*$ if the interval of one node could have been 
	 * generated by the insertion of other repeats inside the occurrence of the repeat
	 * in the other node; -1 if we cannot establish whether the alignment pair represents
	 * such an insertion or not.
	 */
	public static final int isInsertion(IntervalGraph.AlignmentPair pair) {
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final int start1, end1, start2, end2;
		final boolean isLeftMaximal1, isLeftMaximal2, isRightMaximal1, isRightMaximal2;
		final boolean oneInTwo, twoInOne;
		int out;
		start1=pair.node1.start; end1=pair.node1.end;
		start2=pair.node2.start; end2=pair.node2.end;
		alignmentStart1=pair.alignmentStart1;
		alignmentEnd1=pair.alignmentEnd1-1;
		alignmentStart2=pair.alignmentStart2;
		alignmentEnd2=pair.alignmentEnd2-1;
		isLeftMaximal1=pair.node1.isLeftMaximal;
		isRightMaximal1=pair.node1.isRightMaximal;
		isLeftMaximal2=pair.node2.isLeftMaximal;
		isRightMaximal2=pair.node2.isRightMaximal;
	
		if ( Math.abs(alignmentStart1,start1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd1,end1)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 ( ( isLeftMaximal1 && 
			     ( (pair.orientation && alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD) ||
				   (!pair.orientation && alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD && alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD)
			     )
			   ) ||
		       ( isRightMaximal1 && 
			     ( (pair.orientation && alignmentEnd2<end2-IntervalGraph.IDENTITY_THRESHOLD && alignmentStart2>=start2-IntervalGraph.IDENTITY_THRESHOLD) ||
				   (!pair.orientation && alignmentStart2>start2+IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd2<=end2+IntervalGraph.IDENTITY_THRESHOLD)
			     )	 
		       ) 
			 )
		   ) oneInTwo=true;
		else oneInTwo=false;
		if ( Math.abs(alignmentStart2,start2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 Math.abs(alignmentEnd2,end2)<=IntervalGraph.IDENTITY_THRESHOLD && 
			 ( ( isLeftMaximal2 && 
			     ( (pair.orientation && alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD) ||
				   (!pair.orientation && alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD && alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD)
			     )
			   ) ||
			   ( isRightMaximal2 && 
			     ( (pair.orientation && alignmentEnd1<end1-IntervalGraph.IDENTITY_THRESHOLD && alignmentStart1>=start1-IntervalGraph.IDENTITY_THRESHOLD) ||
				   (!pair.orientation && alignmentStart1>start1+IntervalGraph.IDENTITY_THRESHOLD && alignmentEnd1<=end1+IntervalGraph.IDENTITY_THRESHOLD)
		         )
		       )
			 ) 
		   ) twoInOne=true;
		else twoInOne=false;
		out=-1;
		if (oneInTwo && !twoInOne) out=Constants.INSERTION_ONE_IN_TWO;
		else if (!oneInTwo && twoInOne) out=Constants.INSERTION_TWO_IN_ONE;
		else if (oneInTwo && twoInOne) {
			// The two intervals are approximately identical, but they have not been
			// marked as $CONTAINMENT_IDENTICAL$ in previous steps.
		}
		return out;
	}
	
	
	/**
	 * Two intervals just share a substring.
	 */
	public static final int isSharedSubstring(IntervalGraph.AlignmentPair pair) {
		final int threshold = IntervalGraph.IDENTITY_THRESHOLD;
		final int alignmentStart1, alignmentStart2, alignmentEnd1, alignmentEnd2;
		final int start1, end1, start2, end2;
		final int type1, type2;
		int out;
		start1=pair.node1.start; end1=pair.node1.end;
		start2=pair.node2.start; end2=pair.node2.end;
		alignmentStart1=pair.alignmentStart1;
		alignmentEnd1=pair.alignmentEnd1-1;
		alignmentStart2=pair.alignmentStart2;
		alignmentEnd2=pair.alignmentEnd2-1;
		type1=pair.node1.type;
		type2=pair.node2.type;
		
		out=-1;
		if ( alignmentStart1>=start1-threshold && 
			 alignmentEnd1<=end1+threshold && 
			 alignmentStart2>=start2-threshold && 
			 alignmentEnd2<=end2+threshold ) out=Constants.SHARED_SUBSTRING;
		return out;
	}
	
	
	/**
	 * Assume that an interval A is strictly contained in another interval B in the same 
	 * read. Assume also that there is an edge (A,X) that says that A is contained in (but
	 * not identical to) interval X in another read. In Step 3, kernels will be built from 
	 * intervals that are maximal by containment, so A will not be part of a kernel. 
	 * However, if B and X are instances of distinct repeats, A is a string that occurs in
	 * two distinct repeats, so it should not be considered just as a part of X (or of B).
	 * If B and X are instances of the same repeat, then X does not contain evidence for A
	 * to be an independent repeat, but B contains such evidence (otherwise A would not
	 * have been created), so A should not be considered just as a part of X (or of B).
	 * The procedure converts all containment edges of this type into insertion edges.
	 *
	 * Remark: the observations above depend on the maximality of the start/end of A (note
	 * that A strictly contained in B does not imply that A is right- or left-maximal, due
	 * to heuristics in the aligner). The observations above hold for any type of 
	 * interval, but for simplicity the procedure considers just A of alignment interval 
	 * type, and B not short-period (it can be weak dense).
	 *
	 * Remark: a repeat contained in another might be merged with its container later on
	 * in the pipeline, e.g. in the construction of the kernel graph in Step 3. So this 
	 * heuristic does not guarantee that the two will remain distinct in the output.
	 *
	 * Remark: to decide if A is contained, it is crucial to use the $isContainedInterval$
	 * field of an interval graph node, rather than scanning the nodes of the interval 
	 * graph. This is because the procedure might be called also e.g. on a Step 2 
	 * component, which might contain only a subset of all the intervals in a read.
	 *
	 * @param supplementOnly TRUE=works only on supplement edges;
	 * @param containmentEdges temporary space, of size >=max degree of a node;
	 * @return the number of edges changed by the procedure.
	 */
	public static final int containment2insertion(boolean supplementOnly, int[] containmentEdges) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MASK = (1<<Constants.INTERVAL_ALIGNMENT)|(1<<Constants.INTERVAL_DENSE_PREFIX)|(1<<Constants.INTERVAL_DENSE_SUBSTRING)|(1<<(Constants.INTERVAL_PERIODIC+1));
		final int nNodes = IntervalGraph.nNodes;
		int i, j;
		int nNeighbors, containment, lastContainmentEdge, out;
		IntervalGraph.Node nodeI;
		IntervalGraph.Edge edge;
		
		out=0;
		for (i=0; i<nNodes; i++) {
			nodeI=IntervalGraph.nodesArray[i];
			if (nodeI.type!=Constants.INTERVAL_ALIGNMENT || (nodeI.isContainedInterval&MASK)==0) continue;
			nNeighbors=IntervalGraph.nNeighbors[i];
			lastContainmentEdge=-1;
			for (j=0; j<nNeighbors; j++) {
				if ((supplementOnly?IntervalGraph.neighbors[i][j].supplement:true) && containment2insertion_isContainmentEdge(i,j,IDENTITY_THRESHOLD)) containmentEdges[++lastContainmentEdge]=j;
			}
			if (lastContainmentEdge==-1) continue;
			for (j=0; j<=lastContainmentEdge; j++) {
				edge=IntervalGraph.neighbors[i][containmentEdges[j]];
				containment=edge.containment;
				edge.setType_noContainment();
				if (containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==i) edge.insertion=Constants.INSERTION_ONE_IN_TWO;
				else edge.insertion=Constants.INSERTION_TWO_IN_ONE;
				// Fields $containmentSubtype,containmentOverhang*$ used by insertion are
				// also used by containment, so nothing else needs to be changed.
			}
			out+=lastContainmentEdge+1;
		}
		return out;
	}
	
	
	private static final boolean containment2insertion_isContainmentEdge(int nodeID, int edgeID, int identityThreshold) {
		final IntervalGraph.Node node = IntervalGraph.nodesArray[nodeID];
		final IntervalGraph.Edge edge = IntervalGraph.neighbors[nodeID][edgeID];
		final int otherNodeType = IntervalGraph.nodesArray[edge.getTo(nodeID)].type;
		
		if ( otherNodeType>=Constants.INTERVAL_DENSE_PREFIX && otherNodeType<=Constants.INTERVAL_DENSE_SINGLEDELETION &&
			 ( (edge.nodeID1==nodeID && edge.containment==Constants.CONTAINMENT_ONE_IN_TWO) ||
			   (edge.nodeID2==nodeID && edge.containment==Constants.CONTAINMENT_TWO_IN_ONE)
			 ) &&
			 ( ( node.isRightMaximal &&
			     ( ((edge.orientation==0 || edge.orientation==2) && edge.containmentOverhangRight>identityThreshold) ||
			       ((edge.orientation==1 || edge.orientation==2) && edge.containmentOverhangLeft>identityThreshold)
			     )
			   ) ||
			   ( node.isLeftMaximal &&
 			     ( ((edge.orientation==0 || edge.orientation==2) && edge.containmentOverhangLeft>identityThreshold) ||
 			       ((edge.orientation==1 || edge.orientation==2) && edge.containmentOverhangRight>identityThreshold)
 			     )
 			   )
			 )
		   ) return true;
		else return false;
	}
	
	
	/**
	 * Like $containment2insertion()$, but for overlap edges. This is useful, since Step 3
	 * builds kernels following overlap edges. Assume e.g. that an alignment interval A is
	 * at the suffix of a suffix interval B. An occurrence X of B in another read, might 
	 * overlap with A (e.g. because X is not right-maximal), thus A and B would end up in
	 * the same kernel even though they were clearly separated in their read.
	 */
	public static final int overlap2sharedSubstring(boolean supplementOnly, int[] overlapEdges) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MASK = (1<<Constants.INTERVAL_ALIGNMENT)|(1<<Constants.INTERVAL_DENSE_PREFIX)|(1<<Constants.INTERVAL_DENSE_SUBSTRING)|(1<<(Constants.INTERVAL_PERIODIC+1));
		final int nNodes = IntervalGraph.nNodes;
		int i, j;
		int nNeighbors, overlap, lastOverlapEdge, out;
		IntervalGraph.Node nodeI;
		IntervalGraph.Edge edge;
		
		out=0;
		for (i=0; i<nNodes; i++) {
			nodeI=IntervalGraph.nodesArray[i];
			if (nodeI.type!=Constants.INTERVAL_ALIGNMENT || (nodeI.isContainedInterval&MASK)==0) continue;
			nNeighbors=IntervalGraph.nNeighbors[i];
			lastOverlapEdge=-1;
			for (j=0; j<nNeighbors; j++) {
				if ((supplementOnly?IntervalGraph.neighbors[i][j].supplement:true) && overlap2sharedSubstring_isOverlapEdge(i,j,IDENTITY_THRESHOLD)) overlapEdges[++lastOverlapEdge]=j;
			}
			if (lastOverlapEdge==-1) continue;
			for (j=0; j<=lastOverlapEdge; j++) {
				edge=IntervalGraph.neighbors[i][overlapEdges[j]];
				overlap=edge.overlap;
				if (edge.sharedSubstringOverhangs==null) edge.sharedSubstringOverhangs = new int[4];
				Math.set(edge.sharedSubstringOverhangs,3,0);
				if (overlap==Constants.OVERLAP_PREFIX_PREFIX) {
					edge.sharedSubstringOverhangs[0]=0;
					edge.sharedSubstringOverhangs[1]=edge.overhangs[0][0];
					edge.sharedSubstringOverhangs[2]=0;
					edge.sharedSubstringOverhangs[3]=edge.overhangs[0][1];
				}
				else if (overlap==Constants.OVERLAP_PREFIX_SUFFIX) {
					edge.sharedSubstringOverhangs[0]=0;
					edge.sharedSubstringOverhangs[1]=edge.overhangs[1][0];
					edge.sharedSubstringOverhangs[2]=edge.overhangs[1][1];
					edge.sharedSubstringOverhangs[3]=0;
				}
				else if (overlap==Constants.OVERLAP_SUFFIX_PREFIX) {
					edge.sharedSubstringOverhangs[0]=edge.overhangs[2][0];
					edge.sharedSubstringOverhangs[1]=0;
					edge.sharedSubstringOverhangs[2]=0;
					edge.sharedSubstringOverhangs[3]=edge.overhangs[2][1];
				}
				else if (overlap==Constants.OVERLAP_SUFFIX_SUFFIX) {
					edge.sharedSubstringOverhangs[0]=edge.overhangs[3][0];
					edge.sharedSubstringOverhangs[1]=0;
					edge.sharedSubstringOverhangs[2]=edge.overhangs[3][1];
					edge.sharedSubstringOverhangs[3]=0;
				}
				edge.setType_noOverlap();
				edge.sharedSubstring=Constants.SHARED_SUBSTRING;
			}
			out+=lastOverlapEdge+1;
		}
		return out;
	}
	
	
	private static final boolean overlap2sharedSubstring_isOverlapEdge(int nodeID, int edgeID, int identityThreshold) {
		final IntervalGraph.Node node = IntervalGraph.nodesArray[nodeID];
		final IntervalGraph.Edge edge = IntervalGraph.neighbors[nodeID][edgeID];
		final int otherNodeType = IntervalGraph.nodesArray[edge.getTo(nodeID)].type;
		
		if ( otherNodeType>=Constants.INTERVAL_DENSE_PREFIX && otherNodeType<=Constants.INTERVAL_DENSE_SINGLEDELETION &&
			 ( ( edge.nodeID1==nodeID &&
				 ( (node.isRightMaximal && (edge.overlap==Constants.OVERLAP_SUFFIX_PREFIX || edge.overlap==Constants.OVERLAP_SUFFIX_SUFFIX)) ||
				   (node.isLeftMaximal && (edge.overlap==Constants.OVERLAP_PREFIX_PREFIX || edge.overlap==Constants.OVERLAP_PREFIX_SUFFIX))
				 )
			   ) ||
			   ( edge.nodeID2==nodeID &&
   				 ( (node.isRightMaximal && (edge.overlap==Constants.OVERLAP_PREFIX_SUFFIX || edge.overlap==Constants.OVERLAP_SUFFIX_SUFFIX)) ||
   				   (node.isLeftMaximal && (edge.overlap==Constants.OVERLAP_PREFIX_PREFIX || edge.overlap==Constants.OVERLAP_SUFFIX_PREFIX))
   				 )
   			   )
			 )				   
		   ) return true;
		else return false;
	}
	
	
	
	
	
	
	
	
	// ------------------------------ COMPONENT DETECTION --------------------------------
	
	/**
	 * The procedure tries to build clusters of nodes that correspond to families of 
	 * repeats. Since such clusters are ultimately the connected components of some 
	 * relation, they are called "components" throughout.
	 *
	 * The procedure discards every edge whose type is not containment, and assumes that 
	 * the graph consists of a single connected component after such edges are removed. 
	 * In particular, the procedure discards overlap edges, since they could be induced by 
	 * repeats of type $Z = X_2 Y_1$ where $X = X_1 X_2$ and $Y = Y_1 Y_2$ are repeats: 
	 * occurrences of $Z$ without $X$ and $Y$ overlap occurrences of $X$ and of $Y$, 
	 * putting $X$ and $Y$ in the same connected component of the interval graph.
	 *
	 * Then, the procedure continues recursively, discarding the frontier of the graph
	 * (i.e. all minimal and maximal nodes by containment, as well as unary paths adjacent
	 * to such nodes), computing the connected components of the resulting graph, 
	 * propagating component tags to the frontier, and reactivating edges that connect
	 * frontier nodes with just one tag to their connected component. Nodes in the 
	 * frontier might get multiple component tags. The procedure continues recursively on 
	 * each connected component of size at least $MIN_PEELING_SIZE$ (estimated from the
	 * data).
	 *
	 * The procedure is based on peeling frontier nodes, since such nodes are likely to 
	 * connect different repeat families. Indeed, a substring that is shared by different 
	 * repeat families is likely to be rare in the read set, since we need to find 
	 * precisely the right interval in a read, and it is likely to be short, so it is not 
	 * likely to contain other intervals. Factorization might produce intervals that have 
	 * been undersplit, i.e. that are the concatenation of intervals that correspond to 
	 * repeats in different families. Such factorization errors are not likely to be 
	 * contained in other intervals, since they are likely to be rare, and they do not 
	 * have a biological meaning, thus they are not likely to have a full alignment that 
	 * proves they are entirely contained in another interval.
	 *
	 * Remark: clusters inside a connected component of the interval graph with just 
	 * containment edges could be explained as follows. The elements inside a cluster are 
	 * long fragments of a repeat, which tend to align with many other fragments of the
	 * same repeat. The elements that connect different clusters are either short 
	 * fragments of a repeat, which are shared by different repeats and align with 
	 * comparably fewer other fragments (i.e. only with the fragments that contain them), 
	 * or they are concatenations of different repeats, which align with many fragments in
	 * each cluster.
	 *
	 * Remark: maximal nodes by containment that are assigned multiple components by this
	 * procedure are mostly intervals of dense substrings of substring type, or alignment 
	 * intervals that either straddle adjacent repeats, or contain/intersect dense 
	 * substrings of substring type or other complex regions with many maximal events. 
	 * Such regions are hard to factorize and can thus contain multiple repeat units. 
	 * Minimal nodes by containment that are assigned multiple components by this
	 * procedure are obviously mostly alignment intervals.
	 *
	 * Remark: recursion continues only if a frontier is small with respect to its
	 * subgraph. The threshold is inferred from the subgraphs themselves, by first 
	 * performing a full recursion.
	 *
	 * Remark: periodic nodes might be part of the frontier, but this is not likely. See
	 * $IntervalGraph.turnOffFrontier()$ for details.
	 *
	 * Remark: a node with multiple component tags could have a component tag that never
	 * occurs as the only component tag of a node (e.g. at level $i$ of recursion some 
	 * nodes in the frontier are assigned tag $t$, but the subgraph associated with $t$
	 * gets split into connected components at recursion level $i+1$, and all nodes in
	 * $t$ are assigned the tags of such components).
	 *
	 * Remark: given a component tag $c$, the subgraph of the interval graph induced by 
	 * containment edges, and by nodes with just component tag $c$, is not necessarily
	 * connected.
	 *
	 * Remark: the procedure could be parallelized by assigning subintervals to threads.
	 *
	 * @param toSplitFurther writes to this file the IDs of components that should be 
	 * split further, if any (since they come from running a less accurate but faster 
	 * version of clustering on very large connected components);
	 * @return the sorted list of distinct tags of nodes assigned to just one component. 
	 * The tags of all nodes (including nodes assigned to more than one component) are IDs
	 * of nodes in $peelingTree$ (large components), or IDs of small components not 
	 * subject to further processing.
	 */
	private static final int[] getComponents(BufferedWriter toSplitFurther) throws IOException {
		final int MIN_CLIQUE_SIZE = 5;  // Arbitrary
		int i, j;
		int nComponents, lastLargeComponent, nNodesInLargeComponents;
		int firstComponent, nMaximal, nMinimal, maxComponents, maxDegree;
		double threshold;
		int[] subgraph, frontier, tmp, singletonComponentsStats;
		int[][] componentSize, frontierStats;
		IntervalGraph.Node node;
		PeelingNode pNode;
		PriorityQueue queue;

		subgraph = new int[IntervalGraph.nNodes];
		for (i=0; i<IntervalGraph.nNodes; i++) subgraph[i]=i;
		tmp = new int[IntervalGraph.nNodes];
		peelingTree = new PeelingNode(-1);
		
		// Removing non-containment edges
		IntervalGraph.turnOffEdges(0,false);
		for (i=0; i<IntervalGraph.nNodes; i++) IntervalGraph.nodesArray[i].frontier=false;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			node.lastComponent=0;
			if (node.components==null) node.components = new int[IntervalGraph.MIN_COMPONENTS];
			node.components[0]=0;
		}
		nComponents=1; lastLargeComponent=0;
		componentSize = new int[IntervalGraph.nNodes][4];  // Reused by all recursive calls
		componentSize[0][0]=IntervalGraph.nNodes; componentSize[0][1]=0; componentSize[0][2]=0;
		
		// Recursion 1: estimating frontier threshold.
		frontier = new int[IntervalGraph.nNodes];  // Reused by all recursive calls
		threshold=estimateFrontierThreshold(nComponents,lastLargeComponent,subgraph,frontier,tmp,componentSize,MIN_CLIQUE_SIZE);
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 6");
		IntervalGraph.turnOffEdges(0,false);  // Undoing the effect of $estimateFrontierThreshold$ on the ON status of all edges
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 7");
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			node.frontier=false;
			node.lastComponent=0;
			node.components[0]=0;
		}
		nComponents=1; lastLargeComponent=0;
		for (i=0; i<IntervalGraph.nNodes; i++) subgraph[i]=i;
		componentSize[0][0]=IntervalGraph.nNodes; componentSize[0][1]=0; componentSize[0][2]=0;
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 8");
		
		// Recursion 2: peeling
		degrees = new int[IntervalGraph.nNodes];
		clusterField = new int[IntervalGraph.nNodes];
		queue = new PriorityQueue();
		stream = new Stream(256);
		firstComponent=nComponents;
		frontierStats = new int[Constants.INTERVAL_PERIODIC+1][8];
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 11");
		nNodesInLargeComponents=0;
		for (i=0; i<=lastLargeComponent; i++) nNodesInLargeComponents+=componentSize[i][0];
		for (i=0; i<lastLargeComponent; i++) {
			pNode=peelingTree.addChild(componentSize[i][2]);
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 12 getComponents_impl("+i+") start");
			firstComponent+=getComponents_impl(subgraph,componentSize[i][1],componentSize[i+1][1]-1,firstComponent,frontier,tmp,componentSize,pNode,false,threshold,frontierStats,toSplitFurther,MIN_CLIQUE_SIZE);
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 12 getComponents_impl("+i+") end");
		}
		pNode=peelingTree.addChild(componentSize[lastLargeComponent][2]);
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 12 getComponents_impl("+lastLargeComponent+") start");
		firstComponent+=getComponents_impl(subgraph,componentSize[lastLargeComponent][1],nNodesInLargeComponents-1,firstComponent,frontier,tmp,componentSize,pNode,false,threshold,frontierStats,toSplitFurther,MIN_CLIQUE_SIZE);
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 12 getComponents_impl("+lastLargeComponent+") end");
		degrees=null; clusterField=null; queue.clear(); stream.clear(true);
if (IO.SHOW_INTERACTIVE) System.err.println("----> fuck 13");

		IntervalGraph.turnOffEdges(0,true);  // Keeps OFF containment edges to OFF
		IntervalGraph.markContainedContainingIntervals(true);
		System.err.println("STEP2> Maximal/minimal nodes with more than one component in some frontier:");
		nMaximal=0; nMinimal=0;
		for (i=0; i<frontierStats.length; i++) {
			nMaximal+=frontierStats[i][0];
			nMinimal+=frontierStats[i][4];
		}
		for (i=0; i<frontierStats.length; i++) {
			System.err.println(i+", maximal: "+IO.format((100*(double)frontierStats[i][0])/IntervalGraph.nNodes)+"% of all nodes, "+IO.format((100*(double)frontierStats[i][1])/nMaximal)+"% left- or right-max, "+IO.format((100*(double)frontierStats[i][2])/nMaximal)+"% both, "+IO.format((100*(double)frontierStats[i][3])/nMaximal)+"% neither");
			System.err.println(i+", minimal: "+IO.format((100*(double)frontierStats[i][4])/IntervalGraph.nNodes)+"% of all nodes, "+IO.format((100*(double)frontierStats[i][5])/nMinimal)+"% left- or right-max, "+IO.format((100*(double)frontierStats[i][6])/nMinimal)+"% both, "+IO.format((100*(double)frontierStats[i][7])/nMinimal)+"% neither");
		}
		System.err.println("STEP2> Found "+firstComponent+" total components. nNodes="+IntervalGraph.nNodes);
				
		// Assigning component tags to singleton components
		maxComponents=maxComponentTagsPerNode();
		maxDegree=IntervalGraph.getMaxDegree();
		singletonComponentsStats = new int[6];
		getSingletonComponents(singletonComponentsStats,new int[maxComponents*maxDegree],new int[maxComponents*maxDegree]);
		System.err.println("Singleton components: "+singletonComponentsStats[0]+"; no neighbor: "+singletonComponentsStats[1]+"; all neighbors with same components: "+singletonComponentsStats[2]+" (exactly one component: "+singletonComponentsStats[3]+"); at least one periodic neighbor: "+singletonComponentsStats[4]+"; all neighbors periodic: "+singletonComponentsStats[5]);
		cleanSingletonComponents(new int[firstComponent]);
		
		// Computing the distinct tags of nodes with just one tag
		return getComponentTags(firstComponent);
	}
	
	
	/**
	 * Extracts the distinct component tags of nodes with just one component tag. These 
	 * can include components smaller than $MIN_PEELING_SIZE$.
	 *
	 * @param nComponents total number of components;
	 * @return the sorted list of distinct tags.
	 */
	public static final int[] getComponentTags(int nComponents) {
		int i;
		int lastTag, component, nNodesMultipleComponents;
		IntervalGraph.Node node;
		int[] tags, out;
		
		lastTag=-1; nNodesMultipleComponents=0;
		tags = new int[nComponents];
		Math.set(tags,nComponents-1,-1);
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastComponent!=0) {
				if (node.lastComponent>0) nNodesMultipleComponents++;
				else {
					// A node can have no component, iff it is a singleton component with
					// no neighbor.
				}
				continue;
			}
			component=node.components[0];
			if (tags[component]==-1) tags[component]=++lastTag;
		}
		System.err.println("Nodes with multiple components: "+nNodesMultipleComponents+" ("+(((double)nNodesMultipleComponents)/IntervalGraph.nNodes)+"%)");
		out = new int[lastTag+1];
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (node.lastComponent!=0) continue;
			component=node.components[0];
			out[tags[component]]=component;
		}
		if (out.length>1) Arrays.sort(out,0,out.length);
		return out;
	}
	
	
	public static final int[] getComponentsSize(int nComponents) {
		int i, j;
		int[] out = new int[nComponents];	
		
		for (i=0; i<IntervalGraph.nNodes; i++) {
			for (j=0; j<=IntervalGraph.nodesArray[i].lastComponent; j++) out[IntervalGraph.nodesArray[i].components[j]]++;
		}
		return out;
	}
	
	
	
	
	
	/**
	 * The histogram of all ratios between the size of a subgraph built recursively by 
	 * $getComponents_impl()$, and the size of its frontier, should have a peak near one,
	 * corresponding to subgraphs that are mostly frontier and that should not be split. 
	 * A subgraph that should be split consists of a number of components, connected 
	 * by a small frontier.
	 * The procedure returns a value F that separates the peak near one from the rest of 
	 * the distribution. Only subgraphs with ratio bigger than F are assumed to contain 
	 * clusters, and are targets for peeling and splitting.
	 *
	 * Remark: using just the number of maximal/minimal nodes rather than the size of the
	 * frontier ends up removing too many ratios in practice.
	 *
	 * Remark: the procedure assumes that the initial splitting that takes place at the 
	 * root node of the recursion tree, by just removing non-containment edges, has
	 * already been performed, as well as $sortNodesByComponent()$. The procedure also 
	 * assumes $lastLargeComponent>=0$.
	 */
	private static final double estimateFrontierThreshold(int nComponents, int lastLargeComponent, int[] subgraph, int[] frontier, int[] tmp, int[][] componentSize, int minCliqueSize) throws IOException {
		final double MIN_FRONTIER_THRESHOLD = 1.0/0.8;  // Arbitrary
		int i;
		int nNodesInLargeComponents, firstComponent;
		double threshold;

		// Collecting size ratios
		if (tmpPoints==null || tmpPoints.length<subgraph.length) {
			tmpPoints = new Point[subgraph.length];
			for (i=0; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
		}
		lastTmpPoint=-1;
		firstComponent=nComponents;
		nNodesInLargeComponents=0;
		for (i=0; i<=lastLargeComponent; i++) nNodesInLargeComponents+=componentSize[i][0];
		for (i=0; i<lastLargeComponent; i++) firstComponent+=getComponents_impl(subgraph,componentSize[i][1],componentSize[i+1][1]-1,firstComponent,frontier,tmp,componentSize,null,true,-1,null,null,minCliqueSize);
		firstComponent+=getComponents_impl(subgraph,componentSize[lastLargeComponent][1],nNodesInLargeComponents-1,firstComponent,frontier,tmp,componentSize,null,true,-1,null,null,minCliqueSize);
		
		// Estimating threshold
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);
		threshold=Points.getThreshold(tmpPoints,lastTmpPoint,0.0,false,-1);
		threshold=Math.max(threshold,MIN_FRONTIER_THRESHOLD);
		System.err.println("STEP2> Frontier ratios: ");
		for (i=0; i<=lastTmpPoint; i++) System.err.println(tmpPoints[i].position+","+tmpPoints[i].mass);
		System.err.println("STEP2> Threshold="+threshold);
		return threshold;
	}

	
	/**
	 * @param frontier,tmp,componentSize temporary space reused by all recursive calls;
	 * @param root can be null;
	 * @param logFrontierRatio if true, appends to $tmpPoints$ the ratio between the size 
	 * of the subgraph and the size of its frontier, does not start seed set expansion, 
	 * and does not call $IntervalGraph.tagFrontier()$ since it might be too slow in 
	 * practice; this should not change the final choice of the threshold on frontier 
	 * ratios too much;
	 * @param threshold if $threshold>0$, stops recursion if the ratio between subgraph 
	 * size and frontier size is at most $threshold$; otherwise $threshold$ is not used;
	 * @return the total number of connected components (not necessarily large) found by 
	 * this procedure and by all its recursive calls.
	 */
	private static final int getComponents_impl(int[] subgraph, int first, int last, int firstComponent, int[] frontier, int[] tmp, int[][] componentSize, PeelingNode root, boolean logFrontierRatio, double threshold, int[][] frontierStats, BufferedWriter toSplitFurther, int minCliqueSize) throws IOException {
		int i;
		int lastFrontier, decisionNumber;
		int subgraphComponent, nComponents, totalNComponents, lastLargeComponent, nNodesInLargeComponents;
		IntervalGraph.Node tmpNode;
		PeelingNode pNode;
		
		subgraphComponent=IntervalGraph.nodesArray[subgraph[first]].components[0];
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 1");
		IntervalGraph.turnOffFrontier(subgraph,first,last,frontier,tmpArray);
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 2");
		lastFrontier=tmpArray[0];
		decisionNumber=lastFrontier-first+1;  
		if (logFrontierRatio && decisionNumber!=0) {
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=((double)(last-first+1))/decisionNumber;
			tmpPoints[lastTmpPoint].mass=1;
		}
		if ( lastFrontier==last || (threshold>0.0 && decisionNumber!=0 && ((double)(last-first+1))/decisionNumber<=threshold) ) {
			for (i=first; i<=lastFrontier; i++) IntervalGraph.nodesArray[frontier[i]].frontier=false;
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 3");
			IntervalGraph.turnOffEdges(subgraph,first,last,0);  // Resetting to ON all containment edges
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 4");
			if (logFrontierRatio) return 0;
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 5");
			return getComponents_clustering(subgraph,first,last,firstComponent,degrees,clusterField,stream,tmpArray,true,true/* SHOULD BE GRAPHISSMALL... */,minCliqueSize);
		}
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 6");
		IntervalGraph.getConnectedComponents(subgraph,first,last,firstComponent,1,false,tmpArray);
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 7");
		nComponents=tmpArray[1];
		if (nComponents==1) {
			for (i=first; i<=last; i++) {
				tmpNode=IntervalGraph.nodesArray[subgraph[i]];
				tmpNode.lastComponent=0;
				tmpNode.components[0]=subgraphComponent;
			}
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 8");
			IntervalGraph.turnOffEdges(subgraph,first,last,0);  // Resetting to ON all containment edges
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 9");
			if (logFrontierRatio) return 0;
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 10");
			return getComponents_clustering(subgraph,first,last,firstComponent,degrees,clusterField,stream,tmpArray,true,true/* SHOULD BE GRAPHISSMALL... */,minCliqueSize);
		}
		else {
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 11");
			IntervalGraph.tagFrontier(subgraph,first,last,frontier,lastFrontier,tmp,frontierStats,subgraphComponent,logFrontierRatio);
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 12");
		}
		
		// Recursion
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 13");
		lastLargeComponent=IntervalGraph.sortNodesByComponent(subgraph,first,last,firstComponent,nComponents,MIN_PEELING_SIZE,componentSize,tmp);
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 14");
		if (lastLargeComponent<=0) {
			for (i=first; i<=lastFrontier; i++) IntervalGraph.nodesArray[frontier[i]].frontier=false;
			for (i=first; i<=last; i++) {
				tmpNode=IntervalGraph.nodesArray[subgraph[i]];
				tmpNode.lastComponent=0;
				tmpNode.components[0]=subgraphComponent;
			}
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 15");
			if (last+1-first>1) Arrays.sort(subgraph,first,last+1);
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 16");
			IntervalGraph.turnOffEdges(subgraph,first,last,0);  // Resetting to ON all containment edges
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 17");
			if (logFrontierRatio) return 0;
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 18");
			return getComponents_clustering(subgraph,first,last,firstComponent,degrees,clusterField,stream,tmpArray,true,true/* SHOULD BE GRAPHISSMALL... */,minCliqueSize);
		}
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 19");
		nNodesInLargeComponents=0;
		for (i=0; i<=lastLargeComponent; i++) nNodesInLargeComponents+=componentSize[firstComponent+i][0];
		totalNComponents=nComponents;
		for (i=0; i<lastLargeComponent; i++) {
			pNode=root==null?null:root.addChild(componentSize[firstComponent+i][2]);
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 20 ----> component "+i+" start");
			totalNComponents+=getComponents_impl(subgraph,componentSize[firstComponent+i][1],componentSize[firstComponent+i+1][1]-1,firstComponent+totalNComponents,frontier,tmp,componentSize,pNode,logFrontierRatio,threshold,frontierStats,toSplitFurther,minCliqueSize);
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 20 ----> component "+i+" end");
		}
		pNode=root==null?null:root.addChild(componentSize[firstComponent+lastLargeComponent][2]);
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 20 ----> component "+lastLargeComponent+" start");
		totalNComponents+=getComponents_impl(subgraph,componentSize[firstComponent+lastLargeComponent][1],first+nNodesInLargeComponents-1,firstComponent+totalNComponents,frontier,tmp,componentSize,pNode,logFrontierRatio,threshold,frontierStats,toSplitFurther,minCliqueSize);		
if (IO.SHOW_INTERACTIVE) System.err.println("----> getComponents_impl ----> fuck 20 ----> component "+lastLargeComponent+" end");
		return totalNComponents;
	}
	
	
	
	
	
	/**
	 * Prints just the tree of recursion calls of the peeling process. This is not the 
	 * full tree of components, since components smaller than $MIN_PEELING_SIZE$ are not 
	 * part of the tree.
	 */
	private static final void printPeelingTree(String path) throws IOException {
		if (peelingTree==null) {
			System.err.println("ERROR: peelingTree=null");
			System.exit(1);
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
		bw.write("digraph G {");
		peelingTree.toDot(bw);
		bw.write("}");
		bw.close();
		System.err.println("Peeling tree printed to "+path);
	}
	
	
	public static class PeelingNode {
		private static final int MIN_CHILDREN = 10;
		public int id;
		public PeelingNode parent;
		public PeelingNode[] children;
		public int lastChild;
	
		public PeelingNode(int id) {
			this.id=id;
			parent=null;
			children = new PeelingNode[MIN_CHILDREN];
			lastChild=-1;
		}
	
		public final PeelingNode addChild(int newID) {
			PeelingNode tmpNode;
		
			lastChild++;
			if (lastChild==children.length) {
				PeelingNode[] tmp = new PeelingNode[children.length<<1];
				System.arraycopy(children,0,tmp,0,children.length);
				children=tmp;
			}
			tmpNode = new PeelingNode(newID);
			children[lastChild]=tmpNode;
			tmpNode.parent=this;
			return tmpNode;
		}
	
		public final void toDot(BufferedWriter bw) throws IOException {
			for (int i=0; i<=lastChild; i++) bw.write(id+" -> "+children[i].id+" [color=black];\n");
			for (int i=0; i<=lastChild; i++) children[i].toDot(bw);
		}
	
		public final PeelingNode find(int component) {
			if (id==component) return this;
			PeelingNode tmpNode;
			for (int i=0; i<=lastChild; i++) {
				tmpNode=children[i].find(component);
				if (tmpNode!=null) return tmpNode;
			}
			return null;
		}
	
		public final int size() {
			int out = 1;
			for (int i=0; i<=lastChild; i++) out+=children[i].size();
			return out;
		}
		
		public final int toArray(int[] array, int first) {
			array[first]=id;
			int f = first+1;
			for (int i=0; i<=lastChild; i++) f=children[i].toArray(array,f);
			return f;
		}
	}
	
	
	/**
	 * @return the maximum number of component tags assigned to a node.
	 */
	private static final int maxComponentTagsPerNode() {
		int i, last, max;
		
		max=0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			last=IntervalGraph.nodesArray[i].lastComponent;
			if (last>max) max=last;
		}
		return max+1;
	}

	
	/**
	 * Assume that $getComponents()$ has already been executed, and that some components 
	 * have just one node. The procedure propagates component tags from every 
	 * non-singleton-component node, to all singleton-component nodes that can be reached
	 * from it via singleton-component nodes only, using any type of non-supplement edge 
	 * (both ON and OFF). This propagation might disconnect the graph that is induced by 
	 * the nodes with the same component tag after the procedure ends.
	 *
	 * Remark: singleton components can occur because e.g.: (1) some alignments might not 
	 * be assigned to an interval in one or both reads; (2) some dense substrings of
	 * substring type might just share substrings with other dense substrings of substring
	 * type; (3) an interval I of a read R that corresponds to a dense substring of 
	 * substring type might have few neighbors in the interval graph, since most of the 
	 * alignments contained in I in the pile of R might form alignment intervals and get 
	 * assigned to them; (4) a node might be connected to distinct periodic components.
	 *
	 * Remark: the procedure assigns a non-periodic interval to a mostly-periodic 
	 * component if e.g. all its neighbors belong to that component. Thus, components can 
	 * contain both periodic and non-periodic intervals in general. In practice, 
	 * non-periodic intervals might not be concentrated in a single cluster of a 
	 * periodic component.
	 *
	 * Remark: singleton components with no neighbor are assigned no component tag.
	 *
	 * Remark: what was a singleton component might get transformed down the line into a 
	 * kernel of size one, if all the neighbors of the node have one and the same 
	 * component tag. Note however that kernels of size one can also be produced by
	 * iterative peeling.
	 *
	 * Remark: the procedure assumes the $isSingletonComponent$ field of every node to
	 * have already been assigned.	
	 *
	 * @param tmpComponents temporary space, of size at least equal to the number of
	 * distinct components in the interval graph.
	 */
	private static final void cleanSingletonComponents(int[] tmpComponents) {
		int i;
		
		for (i=0; i<IntervalGraph.nNodes; i++) IntervalGraph.nodesArray[i].visited=-1;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			if (IntervalGraph.nodesArray[i].isSingletonComponent) continue;
			propagateComponentTagsToSingleton(i,tmpComponents);
		}
	}
	
	
	/**
	 * Adds the component tags of $node$ to all singleton-component nodes reachable from 
	 * it via singleton-component nodes only. The procedure assumes the $visited$ field of
	 * every node to have been initialized to a negative value.
	 *
	 * Remark: the procedure does not use supplement edges.
	 *
	 * @param tmpComponents temporary space, of size at least equal to the number of
	 * distinct components in the interval graph.
	 */
	private static final void propagateComponentTagsToSingleton(int node, int[] tmpComponents) {
		int j;
		int to, top, last, next, nNeighbors;
		IntervalGraph.Node source, tmpNode;
		IntervalGraph.Edge edge;

		// Initializing the stack
		source=IntervalGraph.nodesArray[node];
		source.visited=node;
		top=-1;
		nNeighbors=IntervalGraph.nNeighbors[node];
		for (j=0; j<nNeighbors; j++) {
			edge=IntervalGraph.neighbors[node][j];
			if (edge.supplement) continue;
			to=edge.getTo(node);
			if (!IntervalGraph.nodesArray[to].isSingletonComponent) continue;
			IntervalGraph.stack[++top]=to;
			IntervalGraph.nodesArray[to].visited=node;
		}
		
		// Reachability
		while (top>=0) {
			to=IntervalGraph.stack[top--];
			tmpNode=IntervalGraph.nodesArray[to];
			last=Math.setUnion(tmpNode.components,tmpNode.lastComponent,source.components,source.lastComponent,tmpComponents);
			if (tmpNode.visited>=0 && last==tmpNode.lastComponent) continue;
			if (tmpNode.components.length<last+1) tmpNode.components = new int[last+1];
			System.arraycopy(tmpComponents,0,tmpNode.components,0,last+1);
			tmpNode.lastComponent=last;
			nNeighbors=IntervalGraph.nNeighbors[to];
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[to][j];
				if (edge.supplement) continue;
				next=edge.getTo(to);
				tmpNode=IntervalGraph.nodesArray[next];
				if (!tmpNode.isSingletonComponent || tmpNode.visited==node) continue;
				tmpNode.visited=node;
				IntervalGraph.stack[++top]=next;
			}
		}
	}
	
	
	/**
	 * Assume that $getComponents()$ has already been executed, and that some components 
	 * have just one node. The procedure sets the $isSingletonComponent$ field of every
	 * such node to TRUE, sets its $lastComponent$ field to -1, and stores the following 
	 * quantities about singleton components in the output array $stats$:
	 *
	 * 0: total number of singleton components;
	 * 1: without any incident edge;
	 * 2: with all non-singleton neighbors having the same set of component tags;
	 * 3: with all non-singleton neighbors having exactly one component tag;
	 * 4: with at least one non-singleton periodic neighbor;
	 * 5: with all non-singleton neighbors of periodic type.
	 *
	 * Remark: the procedure uses only ON edges, and does not use supplement edges.
	 *
	 * @param tmpComponentsFrom,tmpComponentsTo temporary space, of size at least equal to
	 * the maximum number of component tags per node, times the maximum number of 
	 * neighbors of a node.
	 */
	private static final void getSingletonComponents(int[] stats, int[] tmpComponentsFrom, int[] tmpComponentsTo) {
		boolean neighborsHaveSameComponents, periodicNeighbor, allPeriodicNeighbors;
		int i, j;
		int nNeighbors, nGoodNeighbors, lastComponent, last;
		int nSingleton, nSingleton_noNeighbor, nSingleton_sameComponents;
		int nSingleton_oneComponent, nSingleton_atLeastOnePeriodic, nSingleton_allPeriodic;
		IntervalGraph.Node node, neighbor;
		IntervalGraph.Edge edge;
		
		// Setting node fields
		nSingleton=0; nSingleton_noNeighbor=0;
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			j=isSingletonComponent(i);
			if (j==1) {
				node.isSingletonComponent=true;
				node.lastComponent=-1;
				nSingleton++;
			}
			else if (j==0) {
				node.isSingletonComponent=true;
				node.lastComponent=-1;
				nSingleton++;
				nSingleton_noNeighbor++;
			}
			else node.isSingletonComponent=false;
		}
		
		// Collecting statistics
		nSingleton_sameComponents=0; nSingleton_oneComponent=0;
		nSingleton_atLeastOnePeriodic=0; nSingleton_allPeriodic=0; 
		for (i=0; i<IntervalGraph.nNodes; i++) {
			node=IntervalGraph.nodesArray[i];
			if (!node.isSingletonComponent) continue;
			nNeighbors=IntervalGraph.nNeighbors[i];
			
			// Periodic statistics
			periodicNeighbor=false; allPeriodicNeighbors=true; nGoodNeighbors=0;
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[i][j];
				if (edge.supplement) continue;
				neighbor=IntervalGraph.nodesArray[edge.getTo(i)];
				if (neighbor.isSingletonComponent) continue;
				nGoodNeighbors++;
				if (neighbor.type==Constants.INTERVAL_PERIODIC) periodicNeighbor=true;
				if (neighbor.type!=Constants.INTERVAL_PERIODIC) allPeriodicNeighbors=false;
			}
			if (periodicNeighbor) nSingleton_atLeastOnePeriodic++;
			if (nGoodNeighbors>0 && allPeriodicNeighbors) nSingleton_allPeriodic++;
			
			// Component statistics
			lastComponent=-1; neighborsHaveSameComponents=true; nGoodNeighbors=0;
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[i][j];
				if (edge.supplement) continue;
				neighbor=IntervalGraph.nodesArray[edge.getTo(i)];
				if (neighbor.isSingletonComponent) continue;
				nGoodNeighbors++;
				if (lastComponent>=0 && neighbor.lastComponent!=lastComponent) neighborsHaveSameComponents=false;
				last=Math.setUnion(tmpComponentsFrom,lastComponent,neighbor.components,neighbor.lastComponent,tmpComponentsTo);
				if (lastComponent>=0 && last!=lastComponent) neighborsHaveSameComponents=false;
				lastComponent=last;
				System.arraycopy(tmpComponentsTo,0,tmpComponentsFrom,0,lastComponent+1);
			}
			if (nGoodNeighbors>0) {
				if (neighborsHaveSameComponents) nSingleton_sameComponents++;
				if (lastComponent==0) nSingleton_oneComponent++;
			}
		}
		
		// Output
		stats[0]=nSingleton;
		stats[1]=nSingleton_noNeighbor;
		stats[2]=nSingleton_sameComponents;
		stats[3]=nSingleton_oneComponent;
		stats[4]=nSingleton_atLeastOnePeriodic;
		stats[5]=nSingleton_allPeriodic;
	}
	
	
	/**
	 * @return 1 iff $nodesArray[node]$ is the only node of a singleton component by 
	 * containment, and has non-supplement, non-containment neighbors; 
	 * 0 iff $nodesArray[node]$ is the only node of a singleton component by containment, 
	 * and has no non-supplement, non-containment neighbor; 
	 * -1 iff $nodesArray[node]$ is not a singleton component by containment.
	 */
	private static final int isSingletonComponent(int node) {
		int j;
		int nNeighbors, nNeighborsPrime;
		IntervalGraph.Edge edge;
		
		nNeighbors=IntervalGraph.nNeighbors[node];
		nNeighborsPrime=0;
		for (j=0; j<nNeighbors; j++) {
			edge=IntervalGraph.neighbors[node][j];
			if (edge.supplement) continue;
			nNeighborsPrime++;
			if (edge.on && edge.containment!=-1) return -1;
		}
		return nNeighborsPrime==0?0:1;
	}
	
}