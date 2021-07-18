package de.mpi_cbg.revant.fragments;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.CharStream;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.PeriodicSubstrings;
import de.mpi_cbg.revant.intervalgraph.BidirectedGraph;
import de.mpi_cbg.revant.intervalgraph.IntervalGraph;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;


/**
 * Given a set of fragments and a reference sequence, built from the interval graph, the 
 * program tries to build one or more reference sequences that are close to the center of
 * the balls induced by the fragments in edit distance space. This is done by computing a
 * small dominator of the graph induced by fragments and overlap/containment alignments.
 *
 * Remark: the program returns with exit status 255 iff no dominator could be computed.
 * Otherwise it exits with standard exit statuses.
 *
 * Remark: fragments typically belong to highly overlapping balls, and the graph 
 * clustering procedures used elsewhere (like $IntervalGraphStep2.getComponents_
 * clustering()$) are not able to partition them in practice.
 *
 * Remark: a short-period repeat might have few fragment-fragment alignments because the
 * fragments share too many exact k-mers and the aligner caps them. The program tries to
 * check for this, using a simple threshold for "few alignments".
 *
 * Remark: the process of checking if fragments align to the reference, and of building a
 * better reference sequence, could have been done before generating the fragment files, 
 * using the same alignments used in factorization. We do it at this step for simplicity, 
 * and as an additional consistency check over the previous steps.
 */
public class FragmentsStep2 {
	/**
	 * Parameters of the pipeline
	 */
	private static final int IDENTITY_THRESHOLD = IO.quantum;
	private static final int WHITE_NODE = 0;
	private static final int GRAY_NODE = 1;
	private static final int BLACK_NODE = 2;
	private static final int GREEDY_EXTENSION_NODE = 3;
	public static final double DOMINATOR_THRESHOLD = 0.9;  // Arbitrary
	private static int MIN_N_FRAGMENTS;
	public static final int NO_DOMINATOR_FOUND = 255;
	
	/**
	 * Variables for connected dominator
	 */
	private static int[] dominator = new int[300];  // nodeID,iteration,nDominated
	private static QueueElement[] nodeQueue, edgeQueue, queueArray;
	private static QueueElement[] nodeQueueRemove, nodeQueueUpdate;
	private static QueueElement[] edgeQueueRemove, edgeQueueUpdate;
	private static int lastDominator, nodeQueueLast, edgeQueueLast;
	private static int nodeQueueRemoveLast, edgeQueueRemoveLast;
	private static int nodeQueueUpdateLast, edgeQueueUpdateLast;
	
	/**
	 * Header ID generator
	 */
	private static int header; 
	
	/**
	 * Temporary space
	 */
	private static int[] tmpArray = new int[100];
	private static int[] components;
	private static Point[] tmpPoints;
	private static boolean[] isAcyclic;
	private static boolean[][] visited;
	
	
	/**
	 * Remark: the program might build no new reference sequence, if no connected 
	 * component of the graph used to compute the dominator contains enough nodes.
	 *
	 * @param args 
	 * 0: number of fragments (including those in the old reference sequence);
	 * 1: directory containing alignments of fragments (to fragments and reference);
	 * 4: min. number of fragments that must align to a reference;
	 * 6: min length of an alignment used in the repeat inference pipeline;
	 * 7: output dir containing the new references; the new reference file might contain 
	 * more than one string;
	 * 8: output dir containing the new fragments; the new fragments file contains only 
	 * the fragments that align to one of the new references.
	 */
	public static void main(String[] args) throws IOException {
		Reads.nReads=Integer.parseInt(args[0]);
		final String ALIGNMENTS_DIR = args[1];
		final String BASIN_DESCRIPTORS_DIR = args[2];
		final String BASIN_DESCRIPTOR_ID = args[3];
		MIN_N_FRAGMENTS=Integer.parseInt(args[4]);
		final String FRAGMENTS_STRINGS_FILE = args[5];
		final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[6]);
		final String NEW_REFERENCES_DIR = args[7];
		final String NEW_FRAGMENTS_DIR = args[8];
		
		final String READ_LENGTHS_FILE = ALIGNMENTS_DIR+"/"+IO.READS_LENGTHS;
		final String FF_ALIGNMENTS_PATH = ALIGNMENTS_DIR+"/"+IO.ALIGNMENTS_FILE_STEP1;
		final String FR_ALIGNMENTS_PATH = ALIGNMENTS_DIR+"/"+IO.FRAGMENTS_REFERENCE_FILE+".txt";
		final String BASIN_DESCRIPTOR_FILE = BASIN_DESCRIPTORS_DIR+"/basin-"+BASIN_DESCRIPTOR_ID+".txt";
		final String CONNECTION_FILE = ALIGNMENTS_DIR+"/connection.txt";
		final String NEW_REFERENCES_FILE = NEW_REFERENCES_DIR+"/reference-"+BASIN_DESCRIPTOR_ID+".txt";
		final String NEW_REFERENCES_LENGTHS = NEW_REFERENCES_DIR+"/reference-"+BASIN_DESCRIPTOR_ID+"-lengths.txt";
		final String NEW_FRAGMENTS_FILE = NEW_FRAGMENTS_DIR+"/fragments-"+BASIN_DESCRIPTOR_ID+".txt";
		final String NEW_FRAGMENTS_LENGTHS = NEW_FRAGMENTS_DIR+"/fragments-"+BASIN_DESCRIPTOR_ID+"-lengths.txt";
		final String OUTPUT_GRAPH = NEW_FRAGMENTS_DIR+"/graph-"+BASIN_DESCRIPTOR_ID+".dot";  // The graph used by $connectedDominator()$.
		
		final int YIELD_THRESHOLD = 5;  // Arbitrary
		final double ERROR_RATE_QUANTILE = 0.1;  // Arbitrary
		boolean referenceIsPeriodic;
		int i, j;
		int nComponents, maxDegree, referenceLength, referencePeriod;
		double threshold;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		int[] degrees, clusterField, subgraph;
		CharStream[] nodeLabels;
		
		// Building the interval graph
		IO.initialize();
		br = new BufferedReader(new FileReader(BASIN_DESCRIPTOR_FILE));
		str=br.readLine();
		br.close(); br=null;
		IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
		referenceLength=tmpArray[1]; 
		referenceIsPeriodic=tmpArray[3]==Constants.INTERVAL_PERIODIC;
		referencePeriod=tmpArray[4];
		Alignments.minAlignmentLength=MIN_ALIGNMENT_LENGTH;
		Reads.loadReadIDs(0,Reads.nReads-1);
		Reads.loadReadLengths(READ_LENGTHS_FILE);
		IntervalGraph.maxAlignmentsPerRead=IntervalGraph.buildIntervalGraph(FF_ALIGNMENTS_PATH,CONNECTION_FILE,null,null,null,false,false,false,false,false,false,false,null,null,false);
		if (IntervalGraph.nNodes<MIN_N_FRAGMENTS) System.exit(NO_DOMINATOR_FOUND);
		if (referenceIsPeriodic && referencePeriod>0 && referencePeriod<PeriodicSubstrings.MIN_PERIOD_LONG && IntervalGraph.getNEdges()<IntervalGraph.nNodes*FragmentsStep1.N_FRAGMENTS_REFERENCE_ALIGNMENTS_THRESHOLD) {
			System.err.println("FragmentsStep2> Too few fragment-fragment alignments in short-period repeat "+BASIN_DESCRIPTOR_ID+": k-mer capping by the aligner?");
			System.exit(NO_DOMINATOR_FOUND);
		}
		IntervalGraph.sortNodesArray();
		IntervalGraph.setAvgDiffs();
		
		// Dominator
		System.err.println("FragmentsStep2> Building dominator...");
		maxDegree=connectedDominator(DOMINATOR_THRESHOLD,YIELD_THRESHOLD);
		System.err.println("FragmentsStep2> Dominator built ("+((lastDominator+1)/3)+" nodes)");
		if (lastDominator==-1) System.exit(NO_DOMINATOR_FOUND);
		threshold=getErrorRateThreshold(FF_ALIGNMENTS_PATH,ERROR_RATE_QUANTILE);
		System.err.println("FragmentsStep2> Error rate threshold: "+threshold);
		nodeLabels=loadNodeLabels(FRAGMENTS_STRINGS_FILE);
		header=Math.POSITIVE_INFINITY>>1;
		extendDominator(maxDegree,threshold,referenceLength,referenceLength==0?0:(referenceIsPeriodic?referencePeriod:-1),nodeLabels,NEW_REFERENCES_FILE,NEW_REFERENCES_LENGTHS);
		header=Math.POSITIVE_INFINITY>>1;
		printNewFragments(nodeLabels,NEW_FRAGMENTS_FILE,NEW_FRAGMENTS_LENGTHS);
		
		// Printing the graph set up by $connectedDominator()$.
		bw = new BufferedWriter(new FileWriter(OUTPUT_GRAPH));
		bw.write("digraph G {\n");
		printNodes(referenceLength,FR_ALIGNMENTS_PATH,bw);
		printEdges(false,bw);
		bw.write("}\n");
		bw.close();
		
		System.exit(0);
	}

	





	
	// ---------------------------- DOMINATOR PROCEDURES ---------------------------------
	
	/**
	 * Stores in global array $dominator[0..lastDominator]$ a connected dominating set of 
	 * vertices in the subgraph induced by containment and overlap edges only. The 
	 * procedure builds a dominator for each large connected component, and implements the 
	 * approx. O(|E|*|V|) Modified Greedy Algorithm I from \cite{guha1998approximation}.
	 *
	 * Assume that the set of fragments is generated by several variants of a repeat.
	 * Every (unknown) variant is the center of an edit distance ring with fixed radii 
	 * (10% <= error <= 15%) and arbitrary mass (proportional to the frequency of the 
	 * variant in the genome). The ring could actually be a ball if several variants are
	 * very close to one another. The true variant is the consensus of all fragments in 
	 * its ring, and it is likely different from every one of its observed fragments.
	 * Different rings are likely to share fragments.
	 *
	 * If the rings are sufficiently well separated (e.g. if the center of each ring 
	 * does not belong to any other ring), the dominator algorithm is likely to start from
	 * a fragment close to the center of a heaviest ring (since this has largest degree), 
	 * and it is likely to cross two rings when it chooses an edge instead of a node 
	 * (such edge is likely to connect to a fragment that is close to the center of the 
	 * next ring).
	 * 
	 * If two rings are highly overlapping, the dominator algorithm might start from an
	 * off-center node in their intersection (since it might have similar degree to a node
	 * close to a center), and it might then pick other off-center nodes to cover the 
	 * remaining points (since they might have similar yield as nodes close to centers).
	 * To try to avoid this, the procedure chooses a node with largest degree, among all 
	 * nodes with similar max yield. However, this is not likely to be useful when several
	 * rings overlap (e.g. in a "string" of >2 highly overlapping rings), since in this 
	 * case every off-center node might have both high yield and high degree, and the 
	 * algorithm might keep picking off-center nodes at every iteration.
	 * (One could also sort edges with similar yield by the degree of their white nodes,
	 * but this is probably not needed.)
	 *
	 * Remark: the idea of using a connected dominator comes from the fact that, in the 
	 * simple case of a single long repeat, a good reference sequence could be built by 
	 * taking few fragments connected by overlap edges (all other fragments are contained 
	 * in those of overlapping with those). One could model the problem more accurately by
	 * formulating it as a set cover, in which every fragment covers all parts of the 
	 * overall fragments surface that directly align to it; this would guarantee that 
	 * most fragments fully align to the solution, which is something a dominator does not
	 * guarantee.
	 *
	 * Remark: the procedure turns on just containment and overlap edges, and it does NOT
	 * restore the ON state of the other edges when it exits.
	 *
	 * Remark: the procedure stores in the $visited$ field of every node the color 
	 * assigned to it, and in the $iteration$ field of every node the last time its color
	 * was changed.
	 *
	 * @param totalNodesThreshold the procedure builds a dominator to which at least this 
	 * fraction of the total nodes is connected;
	 * @param yieldThreshold when the procedure chooses a node, it picks a node with 
	 * largest degree, among all nodes with yield $>=maxNodeYield-yieldThreshold$;
	 * @return the largest degree of a node, considering only edges used by the procedure.
	 */
	private static int connectedDominator(double totalNodesThreshold, int yieldThreshold) {
		final int GROWTH_RATE = 300;  // Arbitrary, multiple of 3.
		final int nNodes = IntervalGraph.nNodes;
		int i, j;
		int nNeighbors, nEdges, iteration, component;
		int grayNode, whiteNode, nBlackNodes, nWhiteNodes, maxWhiteNodes, nComponents;
		QueueElement candidate;
		IntervalGraph.Edge edge;
		int[] lastWhite = new int[nNodes];
		int[] degree = new int[nNodes];  // The ON degree of each node
		int[] componentSize, maxDegree, maxNode;
		
		// Turning on just containment and overlap edges. Removing self-loops. Sorting
		// edges by $on,getTo()$.
		nEdges=0;
		for (i=0; i<nNodes; i++) {
			nNeighbors=IntervalGraph.nNeighbors[i];
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[i][j];
				edge.on=(edge.isType_overlap()||edge.isType_containment()) && edge.getTo(i)!=i;
			}
			if (nNeighbors>1) {
				IntervalGraph.Edge.order=IntervalGraph.Edge.UNSORTED+1+i;
				Arrays.sort(IntervalGraph.neighbors[i],0,nNeighbors);
			}
			for (j=0; j<nNeighbors; j++) {
				if (!IntervalGraph.neighbors[i][j].on) break;
			}
			lastWhite[i]=j-1; degree[i]=j;
			nEdges+=j;
			IntervalGraph.nodesArray[i].visited=WHITE_NODE;
			IntervalGraph.nodesArray[i].iteration=-1;
		}
		nEdges>>=1;
		
		// Computing connected components
		nComponents=IntervalGraph.getConnectedComponents(1,true,true);
		componentSize = new int[nComponents];
		Math.set(componentSize,nComponents-1,0);
		maxDegree = new int[nComponents];
		Math.set(maxDegree,nComponents-1,0);
		maxNode = new int[nComponents];
		Math.set(maxNode,nComponents-1,-1);
		for (i=0; i<IntervalGraph.nNodes; i++) {
			if (IntervalGraph.nodesArray[i].components==null || IntervalGraph.nodesArray[i].components.length==0) continue;
			component=IntervalGraph.nodesArray[i].components[0];
			componentSize[component]++;
			j=degree[i];
			if (j>maxDegree[component]) {
				maxDegree[component]=j;
				maxNode[component]=i;
			}
		}
		System.err.println("connectedDominator> "+nComponents+" connected components");
		
		// Processing large components
		if (nEdges>0) {
			queueArray = new QueueElement[nEdges];
			nodeQueue = new QueueElement[nNodes];
			edgeQueue = new QueueElement[nEdges];
			nodeQueueRemove = new QueueElement[nNodes];
			nodeQueueUpdate = new QueueElement[nNodes];
			edgeQueueRemove = new QueueElement[nEdges];
			edgeQueueUpdate = new QueueElement[nEdges];
			lastDominator=-1;
			for (i=0; i<nComponents; i++) {
				if (componentSize[i]<MIN_N_FRAGMENTS) continue;
		
				// Adding to the dominator a node with largest ON degree
				iteration=0; nodeQueueLast=-1; edgeQueueLast=-1;
				makeBlack(maxNode[i],lastWhite,true,iteration);
				nWhiteNodes=componentSize[i]-1-maxDegree[i]; nBlackNodes=1;
				dominator[++lastDominator]=maxNode[i];
				dominator[++lastDominator]=0;
				dominator[++lastDominator]=maxDegree[i];
		
				// Iterating
				maxWhiteNodes=(int)((1.0-totalNodesThreshold)*componentSize[i]);
				while (nWhiteNodes>maxWhiteNodes) {
					if (IO.CONSISTENCY_CHECKS) {
						int sum = 0;
						for (int x=0; x<nNodes; x++) {
							if (IntervalGraph.nodesArray[x].components!=null && IntervalGraph.nodesArray[x].components[0]==i && IntervalGraph.nodesArray[x].visited==WHITE_NODE) sum++;
						}
						if (sum!=nWhiteNodes) {
							System.err.println("connectedDominator> ERROR 1: nWhiteNodes="+nWhiteNodes+" but the actual number is "+sum);
							System.exit(1);
						}
						if (nodeQueueLast>=0 && nodeQueue[0].yield>nWhiteNodes) {
							System.err.println("connectedDominator> ERROR 2: maxNodeYield="+nodeQueue[0].yield+" > nWhiteNodes="+nWhiteNodes+": "+nodeQueue[0]);
							System.exit(1);
						}
						if (edgeQueueLast>=0 && edgeQueue[0].yield>nWhiteNodes) {
							System.err.println("connectedDominator> ERROR 3: maxEdgeYield="+edgeQueue[0].yield+" > nWhiteNodes="+nWhiteNodes+": "+edgeQueue[0]);
							System.exit(1);
						}
					}
					iteration++;
					if (nodeQueueLast>=0 && (edgeQueueLast==-1 || edgeQueue[0].yield<(nodeQueue[0].yield<<1)/*Remark after Thm 2.1 in the paper*/)) {
						candidate=nodeQueue[getMaxDegreeNode(yieldThreshold,degree)];
						makeBlack(candidate.max,lastWhite,false,iteration);
						nWhiteNodes-=candidate.yield; nBlackNodes++;
						if (nBlackNodes*3>dominator.length) {
							int[] newArray = new int[dominator.length+GROWTH_RATE];
							System.arraycopy(dominator,0,newArray,0,dominator.length);
							dominator=newArray;
						}
						dominator[++lastDominator]=candidate.max;
						dominator[++lastDominator]=iteration;
						dominator[++lastDominator]=nNodes-nWhiteNodes;
						System.err.println("connectedDominator> added node "+IntervalGraph.nodesArray[candidate.max]);
					}
					else {
						candidate=edgeQueue[0];
						edgeQueueRemove[++edgeQueueRemoveLast]=candidate;
						if (IntervalGraph.nodesArray[candidate.min].visited==GRAY_NODE) {
							grayNode=candidate.min;
							whiteNode=candidate.max;
						}
						else {
							grayNode=candidate.max;
							whiteNode=candidate.min;
						}
						makeBlack(grayNode,lastWhite,false,iteration);
						makeBlack(whiteNode,lastWhite,false,iteration+1/*$iteration+1$ is used here, but not in array $dominator$*/);
						nWhiteNodes-=candidate.yield; nBlackNodes+=2;
						if (nBlackNodes*3>dominator.length) {
							int[] newArray = new int[dominator.length+GROWTH_RATE];
							System.arraycopy(dominator,0,newArray,0,dominator.length);
							dominator=newArray;
						}
						dominator[++lastDominator]=grayNode;
						dominator[++lastDominator]=iteration;
						dominator[++lastDominator]=nNodes-nWhiteNodes;
						dominator[++lastDominator]=whiteNode;
						dominator[++lastDominator]=iteration;
						dominator[++lastDominator]=nNodes-nWhiteNodes;
						iteration++;
						System.err.println("connectedDominator> added edge from "+IntervalGraph.nodesArray[grayNode]+" to "+IntervalGraph.nodesArray[whiteNode]+" yield="+candidate.yield);
					}
				}
			}
		}
		
		// Returning global max degree
		j=0;
		for (i=0; i<nComponents; i++) j=Math.max(j,maxDegree[i]);
		return j;
	}
	
	
	/**
	 * The position in $nodeQueue$ of a node with max degree, among all those with yield
	 * $>=maxYield-yieldThreshold$.
	 *
	 * @param degree array containing the degree of each node.
	 */
	private static final int getMaxDegreeNode(int yieldThreshold, int[] degree) {
		final int maxYield = nodeQueue[0].yield;
		int i;
		int deg, maxDegree, maxNode;
		
		maxNode=0; maxDegree=degree[nodeQueue[0].max];
		for (i=1; i<=nodeQueueLast; i++) {
			if (nodeQueue[i].yield<maxYield-yieldThreshold) break;
			deg=degree[nodeQueue[i].max];
			if (deg>maxDegree) { maxDegree=deg; maxNode=i; }
		}
		return maxNode;
	}
	
	
	/**
	 * Changes the color of $node$ and of its white neighbors, and updates $nodeQueue,
	 * edgeQueue$ to reflect such change.
	 *
	 * Remark: the procedure first collects all edit operations in auxiliary arrays 
	 * ${node,edge}QueueRemove,{node,edge}QueueUpdate$, and then applies all such edits in
	 * batch. In highly connected graphs, this is much faster in practice than querying
	 * priority queues.
	 *
	 * @param firstTime TRUE=this is the first time the procedure is called.
	 */
	private static final void makeBlack(int node, int[] lastWhite, boolean firstTime, int iteration) {
		int i, j, k, p;
		int neighbor, grayNode, whiteNode, lastChanged, nNeighbors;
		int grayNodePrime, nNeighborsPrime;
		IntervalGraph.Node neighborNode;
		IntervalGraph.Edge edge;

		// Changing colors
		IntervalGraph.nodesArray[node].visited=BLACK_NODE;
		IntervalGraph.nodesArray[node].iteration=iteration;
		if (tmpArray.length<lastWhite[node]+1) tmpArray = new int[lastWhite[node]+1];
		for (i=0; i<=lastWhite[node]; i++) {
			neighbor=IntervalGraph.neighbors[node][i].getTo(node);
			IntervalGraph.nodesArray[neighbor].visited=GRAY_NODE;
			IntervalGraph.nodesArray[neighbor].iteration=iteration;
			tmpArray[i]=neighbor;
			removeEdgeFromWhiteList(IntervalGraph.neighbors[node][i],neighbor,lastWhite);
		}
		lastChanged=lastWhite[node];
		
		// Updating the yield of nodes
		nodeQueueRemoveLast=-1; nodeQueueUpdateLast=-1;
		lastWhite[node]=-1;
		if (!firstTime) nodeQueueRemove[++nodeQueueRemoveLast] = new QueueElement(-1,node,-1);
		// We use the $isMaximal$ field of a node to mark if its white list has
		// already been updated.
		for (i=0; i<=lastChanged; i++) {
			// The white list of a changed node does not need to be updated, unless it is
			// adjacent to another changed node.
			IntervalGraph.nodesArray[tmpArray[i]].isMaximal=true;
		}
		for (i=0; i<=lastChanged; i++) {
			grayNode=tmpArray[i];
			nNeighbors=IntervalGraph.nNeighbors[grayNode];
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[grayNode][j];
				if (!edge.on) break;
				neighbor=edge.getTo(grayNode);
				IntervalGraph.nodesArray[neighbor].isMaximal=false;
			}
		}		
		for (i=0; i<=lastChanged; i++) {
			grayNode=tmpArray[i];
			if (!IntervalGraph.nodesArray[grayNode].isMaximal) {
				updateWhiteList(grayNode,lastWhite);
				IntervalGraph.nodesArray[grayNode].isMaximal=true;
			}
			nodeQueue[++nodeQueueLast] = new QueueElement(-1,grayNode,lastWhite[grayNode]+1);
			nNeighbors=IntervalGraph.nNeighbors[grayNode];
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[grayNode][j];
				if (!edge.on) break;
				neighbor=edge.getTo(grayNode);
				neighborNode=IntervalGraph.nodesArray[neighbor];
				if ( !neighborNode.isMaximal && 
				     ( neighborNode.visited==WHITE_NODE || (neighborNode.visited==GRAY_NODE && Arrays.binarySearch(tmpArray,0,lastChanged+1,neighbor)<0) )
				   ) {
					updateWhiteList(neighbor,lastWhite);
					neighborNode.isMaximal=true;
					if (neighborNode.visited==GRAY_NODE/*firstTime is always false here*/) nodeQueueUpdate[++nodeQueueUpdateLast] = new QueueElement(-1,neighbor,lastWhite[neighbor]+1);
				}
			}
		}
		
		// Updating the yield of edges
		edgeQueueRemoveLast=-1; edgeQueueUpdateLast=-1;
		// We use the $isMaximal$ field of a gray (respectively, white) node to mark if
		// the yield of all its adjacent edges to white (resp. gray) neighbors has
		// already been updated.
		for (i=0; i<=lastChanged; i++) {
			grayNode=tmpArray[i];
			nNeighbors=IntervalGraph.nNeighbors[grayNode];
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[grayNode][j];
				if (!edge.on) break;
				neighbor=edge.getTo(grayNode);
				neighborNode=IntervalGraph.nodesArray[neighbor];
				neighborNode.isMaximal=false;
				if (neighborNode.visited!=WHITE_NODE) continue;
				nNeighborsPrime=IntervalGraph.nNeighbors[neighbor];
				for (k=lastWhite[neighbor]+1; k<nNeighborsPrime; k++) {
					edge=IntervalGraph.neighbors[neighbor][k];
					if (!edge.on) break;
					IntervalGraph.nodesArray[edge.getTo(neighbor)].isMaximal=false;
				}
			}
		}
		// Removing black-gray edges
		if (!firstTime) {
			for (i=0; i<=lastChanged; i++) {
				grayNode=tmpArray[i];
				edgeQueueRemove[++edgeQueueRemoveLast] = new QueueElement(Math.min(node,grayNode),Math.max(node,grayNode),-1);
			}
		}
		// Updating the yield of all edges adjacent to gray neighbors of changed nodes
		if (!firstTime) {
			for (i=0; i<=lastChanged; i++) {
				grayNode=tmpArray[i];
				nNeighbors=IntervalGraph.nNeighbors[grayNode];
				for (j=0; j<nNeighbors; j++) {
					edge=IntervalGraph.neighbors[grayNode][j];
					if (!edge.on) break;
					neighbor=edge.getTo(grayNode);
					neighborNode=IntervalGraph.nodesArray[neighbor];
					if (neighborNode.visited!=GRAY_NODE) continue;
					p=Arrays.binarySearch(tmpArray,0,lastChanged+1,neighbor);
					if (p<0 || p>i) edgeQueueRemove[++edgeQueueRemoveLast] = new QueueElement(Math.min(grayNode,neighbor),Math.max(grayNode,neighbor),-1);
					if (!neighborNode.isMaximal) {
						for (k=0; k<=lastWhite[neighbor]; k++) {
							whiteNode=IntervalGraph.neighbors[neighbor][k].getTo(neighbor);	
							edgeQueueUpdate[++edgeQueueUpdateLast] = new QueueElement(Math.min(neighbor,whiteNode),Math.max(neighbor,whiteNode),getYield(neighbor,whiteNode,lastWhite));
						}
						neighborNode.isMaximal=true;
					}
				}
			}
		}
		// Updating the yield of all edges adjacent to white neighbors of changed nodes
		for (i=0; i<=lastChanged; i++) {
			grayNode=tmpArray[i];
			nNeighbors=IntervalGraph.nNeighbors[grayNode];
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[grayNode][j];
				if (!edge.on) break;
				neighbor=edge.getTo(grayNode);
				neighborNode=IntervalGraph.nodesArray[neighbor];
				if (neighborNode.visited!=WHITE_NODE) continue;
				edgeQueue[++edgeQueueLast] = new QueueElement(Math.min(grayNode,neighbor),Math.max(grayNode,neighbor),getYield(grayNode,neighbor,lastWhite));
				if (!firstTime && !neighborNode.isMaximal) {
					nNeighborsPrime=IntervalGraph.nNeighbors[neighbor];
					for (k=lastWhite[neighbor]+1; k<nNeighborsPrime; k++) {
						edge=IntervalGraph.neighbors[neighbor][k];
						if (!edge.on) break;
						grayNodePrime=edge.getTo(neighbor);
						if (grayNodePrime!=grayNode && !IntervalGraph.nodesArray[grayNodePrime].isMaximal) {
							edgeQueueUpdate[++edgeQueueUpdateLast] = new QueueElement(Math.min(neighbor,grayNodePrime),Math.max(neighbor,grayNodePrime),getYield(neighbor,grayNodePrime,lastWhite));
						}
					}
					neighborNode.isMaximal=true;
				}
			}
		}
		
		// Applying all queries in batch.
		// Remark: both node and edge queues have necessarily changed at this point.
		nodeQueueLast=sortQueue(nodeQueue,nodeQueueLast,nodeQueueRemove,nodeQueueRemoveLast,nodeQueueUpdate,nodeQueueUpdateLast);
		edgeQueueLast=sortQueue(edgeQueue,edgeQueueLast,edgeQueueRemove,edgeQueueRemoveLast,edgeQueueUpdate,edgeQueueUpdateLast);
	}
	
	
	/**
	 * Compacts all the white neighbors of $node$ to the beginning of its neighbors list.
	 */
	private static final void updateWhiteList(int node, int[] lastWhite) {
		int i, j;
		int neighbor;
		IntervalGraph.Edge tmpEdge;
		
		j=-1;
		for (i=0; i<=lastWhite[node]; i++) {
			neighbor=IntervalGraph.neighbors[node][i].getTo(node);
			if (IntervalGraph.nodesArray[neighbor].visited==WHITE_NODE) {
				j++;
				if (j!=i) {
					tmpEdge=IntervalGraph.neighbors[node][j];
					IntervalGraph.neighbors[node][j]=IntervalGraph.neighbors[node][i];
					IntervalGraph.neighbors[node][i]=tmpEdge;
				}
			}
		}
		lastWhite[node]=j;
	}
	
	
	/**
	 * Assume that the white neighbors of a node are compacted and sorted by neighbor ID. 
	 * The procedure removes from the white list the edge that is equivalent to $edge$, 
	 * and puts it just outside the white list.
	 */
	private static final void removeEdgeFromWhiteList(IntervalGraph.Edge edge, int node, int[] lastWhite) {
		int i, j;
		final int previousOrder = IntervalGraph.Edge.order;
		final int last = lastWhite[node];
		IntervalGraph.Edge toRemove;
		
		IntervalGraph.Edge.order=IntervalGraph.Edge.UNSORTED+1+node;
		i=Arrays.binarySearch(IntervalGraph.neighbors[node],0,last+1,edge);
		IntervalGraph.Edge.order=previousOrder;
		if (i<0) return;
		toRemove=IntervalGraph.neighbors[node][i];
		for (j=i; j<last; j++) IntervalGraph.neighbors[node][j]=IntervalGraph.neighbors[node][j+1];
		IntervalGraph.neighbors[node][last]=toRemove;
		lastWhite[node]--;
	}
	
	
	/**
	 * Remark: the procedure assumes $neighbors[X][0..lastWhite]$ to contain all and only
	 * the white neighbors of node $X$, sorted by $getTo(X)$.
	 *
	 * Remark: one might think of computing the yield just on containment edges; this is 
	 * not accurate, since a fragment might not be contained in any other fragment but 
	 * have e.g. a suffix-prefix overlap with two fragments in the dominator.
	 *
	 * @param node1,node2 one node is assumed to be gray, the other white;
	 * @return the yield of a gray-white edge, i.e. the size of the union of the sets of 
	 * white neighbors of $node1$ and $node2$ (this number includes the one of $node1,
	 * node2$ that is white).
	 */
	private static final int getYield(int node1, int node2, int[] lastWhite) {
		final int last1 = lastWhite[node1];
		final int last2 = lastWhite[node2];
		int i1, i2, j;
		int from1, from2;
		
		i1=0; i2=0; j=0;
		while (i1<=last1 && i2<=last2) {
			from1=IntervalGraph.neighbors[node1][i1].getTo(node1);
			from2=IntervalGraph.neighbors[node2][i2].getTo(node2);
			if (from1<from2) { j++; i1++; }
			else if (from1>from2) { j++; i2++; }
			else { j++; i1++; i2++; }
		}
		return j+(last1+1-i1)+(last2+1-i2);
	}
	
	
	/**
	 * Stores in $positive$ the list $positive \setminus negative$, sorted by yield, where
	 * yield values are taken from $update$.
	 *
	 * Remark: the procedure uses global array $queueArray$.
	 *
	 * @param positive,negative,update not assumed to be sorted by any criterion;
	 * @return the new value of $lastPositive$.
	 */
	private static final int sortQueue(QueueElement[] positive, int lastPositive, QueueElement[] negative, int lastNegative, QueueElement[] update, int lastUpdate) {
		final int newSize = lastPositive-lastNegative;
		int i, j, k;
		int length;

		QueueElement.order=QueueElement.ORDER_IDS;
		if (lastPositive>0) Arrays.sort(positive,0,lastPositive+1);

		// Applying $update$ to $positive$.
		if (lastUpdate>=0) {
			if (lastUpdate>0) Arrays.sort(update,0,lastUpdate+1); 
			i=0; j=0;
			while (i<=lastUpdate) {
				if (positive[j].compareTo(update[i])==-1) {
					j++;
					continue;
				}
				else if (positive[j].compareTo(update[i])==1) {
					if (IO.CONSISTENCY_CHECKS) {
						System.err.println("sortQueues> ERROR: element to be updated not found?!");
						System.err.println("POSITIVE: ");
						for (int x=0; x<=j; x++) System.err.print(positive[x]+", ");
						System.err.println("...");
						System.err.println("UPDATE: ");
						for (int x=0; x<=i; x++) System.err.print(update[x]+", ");
						System.err.println("...");
						System.exit(1);
					}
				}
				else {
					positive[j].yield=update[i].yield;
					j++; i++; 
				}
			}
		}
		
		// Computing $positive \setminus negative$.
		if (lastNegative>=0) {
			if (queueArray.length<newSize) queueArray = new QueueElement[newSize];
			if (lastNegative>0) Arrays.sort(negative,0,lastNegative+1);
			i=0; j=0; k=-1;
			while (i<=lastNegative) {
				if (j>lastPositive) {
					if (negative[i].min==-1) {  // Error only if it is a node
						System.err.println("sortQueues> ERROR: element to be deleted not found?!");
						System.err.println("POSITIVE: ");
						for (int x=0; x<=lastPositive; x++) System.err.print(positive[x]+", ");
						System.err.println("NEGATIVE: ");
						for (int x=0; x<=i; x++) System.err.print(negative[x]+", ");
						System.err.println("...");
						System.exit(1);
					}
					i++;
					continue;
				}
				if (positive[j].compareTo(negative[i])==-1) {
					queueArray[++k]=positive[j];
					j++;
					continue;
				}
				else if (positive[j].compareTo(negative[i])==1) {
					if (positive[j].min==-1) {  // Error only if it is a node
						System.err.println("sortQueues> ERROR: element to be deleted not found?!");
						System.err.println("POSITIVE: ");
						for (int x=0; x<=j; x++) System.err.print(positive[x]+", ");
						System.err.println("...");
						System.err.println("NEGATIVE: ");
						for (int x=0; x<=i; x++) System.err.print(negative[x]+", ");
						System.err.println("...");
						System.exit(1);
					}
					else i++;
				}
				else { j++; i++; }
			}
			length=lastPositive+1-j;
			if (length>0) {
				System.arraycopy(positive,j,queueArray,k+1,length);
				k+=length;
			}
			lastPositive=k;
			if (lastPositive>=0) System.arraycopy(queueArray,0,positive,0,lastPositive+1);
		}
		
		// Sorting $positive$ by yield.
		if (lastPositive>0) {
			QueueElement.order=QueueElement.ORDER_YIELD;
			Arrays.sort(positive,0,lastPositive+1);
		}
		return lastPositive;
	}
	
	
	private static class QueueElement implements Comparable {
		public static final int ORDER_IDS = 0;
		public static final int ORDER_YIELD = 1;
		public static int order;
		public int min, max, yield;
	
		public QueueElement(int min, int max, int y) {
			this.min=min; this.max=max; yield=y;
		}
	
		public boolean equals(Object other) {
			QueueElement otherElement = (QueueElement)other;
			return min==otherElement.min && max==otherElement.max;
		}
	
		public int compareTo(Object other) {
			QueueElement otherElement = (QueueElement)other;
			if (order==ORDER_IDS) {
				if (min<otherElement.min) return -1;
				else if (min>otherElement.min) return 1;
				if (max<otherElement.max) return -1;
				else if (max>otherElement.max) return 1;
			}
			else if (order==ORDER_YIELD) {
				if (yield<otherElement.yield) return 1;
				else if (yield>otherElement.yield) return -1;
			}
			return 0;
		}
		
		public String toString() {
			return "("+min+","+max+"|"+yield+")";
		}
	}
	
	
	
	
	
	
	
	// -------------------------- REFERENCE SEQUENCE PROCEDURES --------------------------
	
	/**
	 * Remark: due to heuristics in the aligner, the error rate between two reads might
	 * be greater than the one specified to the aligner on the command line (!). E.g. 
	 * in some fragment files, daligner called with identity 0.7 has a distribution of 
	 * error rates with most mass between 0.2 and 0.35, centered at 0.27, and stretching 
	 * from 0.16 to 0.37.
	 *
	 * @return a quantile of the distribution of error rates of fragment-fragment 
	 * alignments, which is expected to be unimodal. Computing error rates on the edges of
	 * the interval graph rather than on all input alignments should give similar results.
	 */
	private static final double getErrorRateThreshold(String ffAlignmentsFile, double quantile) throws IOException {
		int i;
		int nAlignments, lastTmpPoint;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(ffAlignmentsFile));
		br.readLine(); br.readLine();  // Skipping header
		str=br.readLine(); nAlignments=0;
		while (str!=null) {
			nAlignments++;
			str=br.readLine();
		}
		br.close();
		if (tmpPoints==null || tmpPoints.length<nAlignments) {
			tmpPoints = new Point[nAlignments];
			for (i=0; i<tmpPoints.length; i++) tmpPoints[i] = new Point();
		}
		br = new BufferedReader(new FileReader(ffAlignmentsFile));
		br.readLine(); br.readLine();  // Skipping header
		str=br.readLine(); lastTmpPoint=-1;
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			lastTmpPoint++;
			tmpPoints[lastTmpPoint].position=((double)(Alignments.diffs<<1))/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2);
			tmpPoints[lastTmpPoint].mass=1;
			str=br.readLine();
		}
		br.close();
		lastTmpPoint=Points.sortAndCompact(tmpPoints,lastTmpPoint);		
		return Points.quantile(tmpPoints,0,lastTmpPoint,true,quantile,false);
	}
	
	
	/**
	 * Grows new reference sequences from the nodes of the dominator, and prints them to
	 * $outputFile$. Specifically, the procedure builds the bidirected graph induced by
	 * dominator nodes and by all their overlap edges, and lists all paths. Every such 
	 * path is then extended greedily, using any remaining node in the interval graph, and
	 * possibly connecting a bidirected path to another bidirected path. See 
	 * $buildBidirectedGraph(),extendBidirectedPaths(),concatenateExtendedPaths()$ for 
	 * details.
	 *
	 * Remark: if the repeat is periodic, bidirected paths are not greedily extended, 
	 * since every fragment is likely to overlap with every other fragment. If the repeat 
	 * is short-period, bidirected paths are not even built, and every node in the 
	 * dominator is reported as a separate reference. If the repeat is long-period, only 
	 * dominator nodes that are shorter than the period can have incident edges in the
	 * bidirected graph.
	 *
	 * Remark: the procedure stores the ID of the ring cluster a node belongs to (if any)
	 * in its field $iteration$.
	 * 
	 * @param threshold output of $getErrorRateThreshold()$;
	 * @param referenceLength string length of the current reference; this is assumed to 
	 * be a good estimate of the string length of every non-periodic variant (it is not
	 * used if the repeat is periodic);
	 * @param period estimate of the period of the repeat: -1=non-periodic; 0=periodic, 
	 * unknown period; >0=periodic, known period.
	 */
	private static final void extendDominator(int maxDegree, double threshold, int referenceLength, int period, CharStream[] nodeLabels, String outputFile, String outputLengths) throws IOException {
		final int MAX_LENGTH = referenceLength+(IO.quantum<<2);  // Arbitrary
		final int CAPACITY = 10000;  // Arbitrary
		int i, j, k;
		int nNodes, maxNodes, nBalls, previousIteration, start;
		IntervalGraph.Node tmpNode;
		CharStream tmpStream1, tmpStream2, tmpStream3;
		boolean[] concatenated;
		CharStream[] extendedPathsLabels, bidirectedGraphStrings;
		boolean[][] pathEndExtended;
		int[][] nOverlaps, path2path;
			 
		// Partitioning the graph where the dominator alg. selected an edge over a node
		if (lastDominator==2) { nBalls=1; maxNodes=1; }
		else {
			if (tmpArray==null || tmpArray.length<lastDominator/3)  tmpArray = new int[lastDominator/3];
			nBalls=1; nNodes=1; maxNodes=1; previousIteration=dominator[1];
			for (i=3; i<lastDominator; i+=3) {
				if (dominator[i+1]==previousIteration) {
					maxNodes=Math.max(maxNodes,nNodes);
					tmpArray[nBalls-1]=previousIteration;
					nBalls++; nNodes=1;
				}
				else {
					nNodes++;
					previousIteration=dominator[i+1];
				}
			}
			maxNodes=Math.max(maxNodes,nNodes);
		}
		System.err.println("extendDominator> "+nBalls+" ring clusters found");
		nNodes=IntervalGraph.nNodes;
		if (nBalls==1) {
			for (i=0; i<nNodes; i++) {
				tmpNode=IntervalGraph.nodesArray[i];
				if (tmpNode.iteration!=-1) tmpNode.iteration=0;
			}
		}
		else {
			for (i=0; i<nNodes; i++) {
				tmpNode=IntervalGraph.nodesArray[i];
				if (tmpNode.iteration!=-1) {
					j=Arrays.binarySearch(tmpArray,0,nBalls-1,tmpNode.iteration);
					tmpNode.iteration=j>=0?j:-1-j;
				}
			}
		}
		
		// Initializing reused data structures
		extendedPathsLabels = new CharStream[maxNodes];
		for (i=0; i<maxNodes; i++) extendedPathsLabels[i] = new CharStream(CAPACITY);
		if (period==-1) {
			path2path = new int[maxNodes][4];
			pathEndExtended = new boolean[maxNodes][2];
			bidirectedGraphStrings = new CharStream[maxNodes];
			concatenated = new boolean[maxNodes];
			tmpStream1 = new CharStream(CAPACITY); 
			tmpStream2 = new CharStream(CAPACITY);
			tmpStream3 = new CharStream(CAPACITY);
		}
		else {
			path2path=null;
			pathEndExtended=null;
			bidirectedGraphStrings=null;
			concatenated=null;
			tmpStream1=null;
			tmpStream2=null;
			tmpStream3=null;
		}
		for (i=0; i<nNodes; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			tmpNode.isMaximal=false;  // End of a bidirected graph path
			tmpNode.isMinimal=false;  // End of an extended bidirected graph path
			tmpNode.isSingletonComponent=false;  // In current bidirected graph path
			tmpNode.frontier=false;  // Greedy extension of a bidirected graph path
		}
		
		// Building a distinct bidirected graph for every cluster of rings
		if (lastDominator==2) {
			if (period>=0) {  // No greedy extension
				tmpNode=IntervalGraph.nodesArray[dominator[0]];
				extendedPathsLabels[0].clear(false);
				extendedPathsLabels[0].push(nodeLabels[dominator[0]],0,tmpNode.length()-1,0);
			}
			else {
				BidirectedGraph.allocateMemory(1,1,1);
				BidirectedGraph.clear();
				BidirectedGraph.nNodes=1; BidirectedGraph.lastNeighbor[0]=-1; 
				tmpNode=IntervalGraph.nodesArray[dominator[0]];
				BidirectedGraph.nodeLength[0]=tmpNode.length();
				BidirectedGraph.intervalGraphPointers[0]=tmpNode;
				BidirectedGraph.lastAssemblyPath=0;
				BidirectedGraph.assemblyPaths = new int[1][1];
				BidirectedGraph.assemblyPaths[0][0]=0;
				BidirectedGraph.assemblyPathLengths = new int[1];
				BidirectedGraph.assemblyPathLengths[0]=1;
				BidirectedGraph.assemblyPathStringLengths = new int[1];
				BidirectedGraph.assemblyPathStringLengths[0]=tmpNode.length();
				extendBidirectedPaths(0,MAX_LENGTH,maxNodes,nodeLabels,extendedPathsLabels,path2path,pathEndExtended,bidirectedGraphStrings,tmpStream1,tmpStream2,tmpStream3);
			}
			printPaths(extendedPathsLabels,1,outputFile,outputLengths);
		}
		else {
			nOverlaps = new int[(lastDominator+1)/3][2];
			buildNOverlaps(threshold,nOverlaps);
			BidirectedGraph.allocateMemory(maxNodes,(maxNodes*maxDegree)>>1,maxDegree);
			nBalls=1; start=0; nNodes=1; previousIteration=dominator[1];
			for (i=3; i<lastDominator; i+=3) {
				if (dominator[i+1]==previousIteration) {
					if (i-1==start) {
						if (period>=0) {  // No greedy extension
							tmpNode=IntervalGraph.nodesArray[dominator[start]];
							extendedPathsLabels[0].clear(false);
							extendedPathsLabels[0].push(nodeLabels[dominator[start]],0,tmpNode.length()-1,0);
						}
						else {
							BidirectedGraph.clear();
							BidirectedGraph.nNodes=1; BidirectedGraph.lastNeighbor[0]=-1; 
							tmpNode=IntervalGraph.nodesArray[dominator[start]];
							BidirectedGraph.nodeLength[0]=tmpNode.length();
							BidirectedGraph.intervalGraphPointers[0]=tmpNode;
							BidirectedGraph.lastAssemblyPath=0;
							BidirectedGraph.assemblyPaths = new int[1][1];
							BidirectedGraph.assemblyPaths[0][0]=0;
							BidirectedGraph.assemblyPathLengths = new int[1];
							BidirectedGraph.assemblyPathLengths[0]=1;
							BidirectedGraph.assemblyPathStringLengths = new int[1];
							BidirectedGraph.assemblyPathStringLengths[0]=tmpNode.length();
							extendBidirectedPaths(0,MAX_LENGTH,maxNodes,nodeLabels,extendedPathsLabels,path2path,pathEndExtended,bidirectedGraphStrings,tmpStream1,tmpStream2,tmpStream3);
						}
						printPaths(extendedPathsLabels,1,outputFile,outputLengths);
					}
					else {
						buildBidirectedGraph(start,i-1,threshold,nOverlaps,period);
						BidirectedGraph.toDot("./bidirectedGraph-"+(nBalls-1)+".dot",null);
						if (period>=0) {  // No greedy extension nor concatenation
							for (k=0; k<BidirectedGraph.nNodes; k++) {
								extendedPathsLabels[k].clear(false);
								extendedPathsLabels[k].push(nodeLabels[BidirectedGraph.intervalGraphPointers[k].nodeID],0,BidirectedGraph.nodeLength[k]-1,0);
							}
							printPaths(extendedPathsLabels,BidirectedGraph.nNodes,outputFile,outputLengths);
						}
						else {
							if (path2path.length<BidirectedGraph.lastAssemblyPath+1) path2path = new int[BidirectedGraph.lastAssemblyPath+1][4];
							if (pathEndExtended.length<BidirectedGraph.lastAssemblyPath+1) pathEndExtended = new boolean[BidirectedGraph.lastAssemblyPath+1][2];
							if (concatenated.length<BidirectedGraph.lastAssemblyPath+1) concatenated = new boolean[BidirectedGraph.lastAssemblyPath+1];
							if (extendedPathsLabels.length<BidirectedGraph.lastAssemblyPath+1) {
								CharStream[] newStrings = new CharStream[BidirectedGraph.lastAssemblyPath+1];
								System.arraycopy(extendedPathsLabels,0,newStrings,0,extendedPathsLabels.length);
								for (k=extendedPathsLabels.length; k<newStrings.length; k++) newStrings[k] = new CharStream(CAPACITY);
								extendedPathsLabels=newStrings; 
							}
							if (bidirectedGraphStrings.length<BidirectedGraph.lastAssemblyPath+1) {
								CharStream[] newStrings = new CharStream[BidirectedGraph.lastAssemblyPath+1];
								System.arraycopy(bidirectedGraphStrings,0,newStrings,0,bidirectedGraphStrings.length);
								for (k=bidirectedGraphStrings.length; k<newStrings.length; k++) newStrings[k] = new CharStream(CAPACITY);
								bidirectedGraphStrings=newStrings; 
							}
							extendBidirectedPaths(nBalls-1,MAX_LENGTH,maxNodes,nodeLabels,extendedPathsLabels,path2path,pathEndExtended,bidirectedGraphStrings,tmpStream1,tmpStream2,tmpStream3);
							if (BidirectedGraph.lastAssemblyPath==0) printPaths(extendedPathsLabels,1,outputFile,outputLengths);
							else concatenateExtendedPaths(extendedPathsLabels,path2path,outputFile,outputLengths,concatenated,tmpStream1);
						}
					}
					nBalls++; start=i; nNodes=1;
				}
				else previousIteration=dominator[i+1];
			}
			if (lastDominator==start) {
				if (period>=0) {  // No greedy extension
					tmpNode=IntervalGraph.nodesArray[dominator[start]];
					extendedPathsLabels[0].clear(false);
					extendedPathsLabels[0].push(nodeLabels[dominator[start]],0,tmpNode.length()-1,0);
				}
				else {
					BidirectedGraph.allocateMemory(1,1,1);
					BidirectedGraph.clear();
					BidirectedGraph.nNodes=1; BidirectedGraph.lastNeighbor[0]=-1; 
					tmpNode=IntervalGraph.nodesArray[dominator[start]];
					BidirectedGraph.nodeLength[0]=tmpNode.length();
					BidirectedGraph.intervalGraphPointers[0]=tmpNode;
					BidirectedGraph.lastAssemblyPath=0;
					BidirectedGraph.assemblyPaths = new int[1][1];
					BidirectedGraph.assemblyPaths[0][0]=0;
					BidirectedGraph.assemblyPathLengths = new int[1];
					BidirectedGraph.assemblyPathLengths[0]=1;
					BidirectedGraph.assemblyPathStringLengths = new int[1];
					BidirectedGraph.assemblyPathStringLengths[0]=tmpNode.length();
					extendBidirectedPaths(0,MAX_LENGTH,maxNodes,nodeLabels,extendedPathsLabels,path2path,pathEndExtended,bidirectedGraphStrings,tmpStream1,tmpStream2,tmpStream3);
				}
				printPaths(extendedPathsLabels,1,outputFile,outputLengths);
			}
			else {
				buildBidirectedGraph(start,lastDominator,threshold,nOverlaps,period);
				BidirectedGraph.toDot("./bidirectedGraph-"+(nBalls-1)+".dot",null);
				if (period>=0) {  // No greedy extension nor concatenation
					for (k=0; k<BidirectedGraph.nNodes; k++) {
						extendedPathsLabels[k].clear(false);
						extendedPathsLabels[k].push(nodeLabels[BidirectedGraph.intervalGraphPointers[k].nodeID],0,BidirectedGraph.nodeLength[k]-1,0);
					}
					printPaths(extendedPathsLabels,BidirectedGraph.nNodes,outputFile,outputLengths);
				}
				else {
					if (path2path.length<BidirectedGraph.lastAssemblyPath+1) path2path = new int[BidirectedGraph.lastAssemblyPath+1][4];
					if (pathEndExtended.length<BidirectedGraph.lastAssemblyPath+1) pathEndExtended = new boolean[BidirectedGraph.lastAssemblyPath+1][2];
					if (concatenated.length<BidirectedGraph.lastAssemblyPath+1) concatenated = new boolean[BidirectedGraph.lastAssemblyPath+1];
					if (extendedPathsLabels.length<BidirectedGraph.lastAssemblyPath+1) {
						CharStream[] newStrings = new CharStream[BidirectedGraph.lastAssemblyPath+1];
						System.arraycopy(extendedPathsLabels,0,newStrings,0,extendedPathsLabels.length);
						for (k=extendedPathsLabels.length; k<newStrings.length; k++) newStrings[k] = new CharStream(CAPACITY);
						extendedPathsLabels=newStrings; 
					}
					if (bidirectedGraphStrings.length<BidirectedGraph.lastAssemblyPath+1) {
						CharStream[] newStrings = new CharStream[BidirectedGraph.lastAssemblyPath+1];
						System.arraycopy(bidirectedGraphStrings,0,newStrings,0,bidirectedGraphStrings.length);
						for (k=bidirectedGraphStrings.length; k<newStrings.length; k++) newStrings[k] = new CharStream(CAPACITY);
						bidirectedGraphStrings=newStrings; 
					}
					extendBidirectedPaths(nBalls-1,MAX_LENGTH,maxNodes,nodeLabels,extendedPathsLabels,path2path,pathEndExtended,bidirectedGraphStrings,tmpStream1,tmpStream2,tmpStream3);
					if (BidirectedGraph.lastAssemblyPath==0) printPaths(extendedPathsLabels,1,outputFile,outputLengths);
					else concatenateExtendedPaths(extendedPathsLabels,path2path,outputFile,outputLengths,concatenated,tmpStream1);
				}
			}
		}
		
		// Deallocating reused data structures
		for (i=0; i<maxNodes; i++) {
			extendedPathsLabels[i].deallocate(); extendedPathsLabels[i]=null;
		}
		extendedPathsLabels=null;
		if (period==-1) {
			for (i=0; i<maxNodes; i++) bidirectedGraphStrings[i]=null;
			bidirectedGraphStrings=null;
			tmpStream1.deallocate(); tmpStream1=null;
			tmpStream2.deallocate(); tmpStream2=null;
			tmpStream3.deallocate(); tmpStream3=null;
		}
	}
	
	
	/**
	 * Builds the bidirected graph induced by all nodes in $dominator[start..end]$ and by
	 * the overlap edges between them, and prints all paths. An overlap edge between 
	 * dominator nodes is not used if it has error rate $>threshold$ and if there are 
	 * other overlap edges, on the same end of a node, with error rate $<=threshold$. This
	 * is done to avoid assembling the centers of different rings that are very close to 
	 * one another.
	 *
	 * Remark: the procedure does not merge a dominator node that is contained in another
	 * dominator node, since both their sequences are necessary to cover the specificied 
	 * fraction of all nodes of the graph. This is why containment edges are not used by
	 * the procedure.
	 *
	 * Remark: extending every node of the dominator greedily, in isolation, without first
	 * trying to assemble dominator nodes, likely avoids merging the dominator nodes at 
	 * the center of very close rings, but it is probably too fragile: by picking the 
	 * absolute best at every step we might not sucessfully merge dominator nodes that
	 * should be merged.
	 *
	 * Remark: if the repeat is short-period no edge is built, since every fragment is
	 * likely to overlap with every other; if the repeat is long-period, only fragments
	 * shorter than the period can be incident to overlap edges.
	 *
	 * @param nOverlaps output of $buildNOverlaps()$;
	 * @param period -1=non-periodic; 0=periodic, unknown period; >0=periodic, known 
	 * period.
	 */
	private static final void buildBidirectedGraph(int start, int end, double threshold, int[][] nOverlaps, int period) {
		final int PERIOD_THRESHOLD = period-IO.quantum;  // Arbitrary
		boolean overlapAlignments;
		int i, j, k;
		int from, to, nodeID, nNeighbors, node1, node2, node1Prime, node2Prime, component, nComponents;
		IntervalGraph.Node tmpNode;
		IntervalGraph.Edge intervalEdge;
		BidirectedGraph.BidirectedEdge bidirectedEdge;
	
		BidirectedGraph.clear();
		BidirectedGraph.nNodes=(end-start+1)/3;
		for (i=start; i<end; i+=3) {
			j=(i-start)/3;
			tmpNode=IntervalGraph.nodesArray[dominator[i]];
			BidirectedGraph.nodeLength[j]=tmpNode.length();
			BidirectedGraph.intervalGraphPointers[j]=tmpNode;
		}
		if (period>=0 && period<PeriodicSubstrings.MIN_PERIOD_LONG) {
			for (i=0; i<BidirectedGraph.nNodes; i++) BidirectedGraph.lastNeighbor[i]=-1;
			BidirectedGraph.printAllPaths_trivial(0);
			return;
		}
		for (i=start; i<end; i+=3) {
			nodeID=dominator[i];
			nNeighbors=IntervalGraph.nNeighbors[nodeID];
			for (j=0; j<nNeighbors; j++) {
				intervalEdge=IntervalGraph.neighbors[nodeID][j];
				if (!intervalEdge.on) break;
				intervalEdge.printed=-1;  // Using the $printed$ field of each interval graph edge to mark whether an edge has already been used by this procedure.
			}
		}
		for (i=start; i<end; i+=3) {
			from=dominator[i];
			if (period>=0 && IntervalGraph.nodesArray[from].length()>=PERIOD_THRESHOLD) continue;
			nNeighbors=IntervalGraph.nNeighbors[from];
			for (j=0; j<nNeighbors; j++) {
				intervalEdge=IntervalGraph.neighbors[from][j];
				if (!intervalEdge.on) break;
				if (intervalEdge.printed==1) continue;
				intervalEdge.printed=1;
				if (!intervalEdge.isType_overlap()) continue;
				to=intervalEdge.getTo(from);
				if (period>=0 && IntervalGraph.nodesArray[to].length()>=PERIOD_THRESHOLD) continue;
				if (from==intervalEdge.nodeID1) {
					node1=(i-start)/3;
					node1Prime=i/3;
					node2=buildBidirectedGraph_findNode(to,start,end);
					if (node2!=-1) {
						node2=(node2-start)/3;
						node2Prime=node2/3;
					}
					else node2Prime=-1;
				}
				else {
					node1=buildBidirectedGraph_findNode(to,start,end);
					if (node1!=-1) {
						node1=(node1-start)/3;
						node1Prime=node1/3;
					}
					else node1Prime=-1;
					node2=(i-start)/3;
					node2Prime=i/3;
				}
				if (node1==-1 || node2==-1) continue;
				overlapAlignments=intervalEdge.nAlignmentsOfType(1)!=Math.POSITIVE_INFINITY;
				for (k=Constants.OVERLAP_PREFIX_PREFIX; k<=Constants.OVERLAP_SUFFIX_SUFFIX; k<<=1) {
					if ((intervalEdge.overlap&k)==0) continue;
					if ( intervalEdge.avgDiffs>threshold &&
						 ( (k==Constants.OVERLAP_PREFIX_PREFIX && (nOverlaps[node1Prime][0]>1 || nOverlaps[node2Prime][0]>1)) ||
						   (k==Constants.OVERLAP_PREFIX_SUFFIX && (nOverlaps[node1Prime][0]>1 || nOverlaps[node2Prime][1]>1)) ||
						   (k==Constants.OVERLAP_SUFFIX_PREFIX && (nOverlaps[node1Prime][1]>1 || nOverlaps[node2Prime][0]>1)) ||
						   (k==Constants.OVERLAP_SUFFIX_SUFFIX && (nOverlaps[node1Prime][1]>1 || nOverlaps[node2Prime][1]>1))
						 )
					   ) continue;
					bidirectedEdge=BidirectedGraph.getEdge();
					bidirectedEdge.node1=node1; 
					bidirectedEdge.node2=node2; 
					bidirectedEdge.type=k;
					bidirectedEdge.setOverhangs(k,intervalEdge.overhangs);
					bidirectedEdge.nAlignments=overlapAlignments?1:2;  // Assigning 2 alignments is arbitrary.
					BidirectedGraph.addEdge(bidirectedEdge);
				}
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			System.err.print("buildBidirectedGraph> Consistency checks started... ");
			BidirectedGraph.checkConsistency();
			System.err.println("OK");
		}
		
		// Printing all paths of all acyclic connected components
		if (components==null || BidirectedGraph.nNodes>components.length) components = new int[BidirectedGraph.nNodes];
		BidirectedGraph.printAllPaths_init();
		if (BidirectedGraph.nOneEnd==0) {
			for (i=0; i<BidirectedGraph.nNodes; i++) components[i]=0;
		}
		else {
			nComponents=BidirectedGraph.getConnectedComponents(components);
			if (isAcyclic==null || nComponents>isAcyclic.length) isAcyclic = new boolean[nComponents];
			Math.set(isAcyclic,nComponents-1,true);
			if (visited==null || BidirectedGraph.nNodes>visited.length) visited = new boolean[BidirectedGraph.nNodes][2];
			Math.set(visited,BidirectedGraph.nNodes-1,false);
			for (i=0; i<BidirectedGraph.nOneEnd; i++) {
				from=BidirectedGraph.oneEnd[i]>=0?BidirectedGraph.oneEnd[i]:-1-BidirectedGraph.oneEnd[i];
				if (!isAcyclic[components[from]]) continue;
				component=components[from];
				if (BidirectedGraph.cycleFrom(from,BidirectedGraph.oneEnd[i]>=0?0:1,visited)) isAcyclic[component]=false;
			}
			for (i=0; i<BidirectedGraph.nNodes; i++) components[i]=isAcyclic[components[i]]?1:0;
		}
		BidirectedGraph.printAllPaths(0,true,components);
	}
	
	
	/**
	 * @return the position of $node$ in $dominator[start..end]$, or -1 if not found. 
	 * This is just a linear scan, since $end-start+1$ is assumed to be short.
	 */
	private static final int buildBidirectedGraph_findNode(int node, int start, int end) {
		for (int i=start; i<end; i+=3) {
			if (dominator[i]==node) return i;
		}
		return -1;
	}
	
	
	/**
	 * Stores in $out[i]$ the number of overlap edges with error rate $<=threshold$ 
	 * that use the prefix (out[i][0]) or the suffix (out[i][1]) of the $i$-th node in
	 * $dominator$. $out$ is assumed to be already big enough.
	 */
	private static final void buildNOverlaps(double threshold, int[][] out) {
		int i, j;
		int node, nNeighbors;
		IntervalGraph.Edge edge;
		
		Math.set(out,0);
		for (i=0; i<lastDominator; i+=3) {
			node=dominator[i];
			nNeighbors=IntervalGraph.nNeighbors[node];
			for (j=0; j<nNeighbors; j++) {
				edge=IntervalGraph.neighbors[node][j];
				if (!edge.on) break;
				if (!edge.isType_overlap() || edge.avgDiffs>threshold) continue;
				if (edge.isPrefixOverlap(node)) out[i/3][0]++;
				if (edge.isSuffixOverlap(node)) out[i/3][1]++;
			}
		}
	}
	
	
	/**
	 * Tries to extend every path of the bidirected graph of ring cluster $ballID$, using
	 * best-scoring edges in the interval graph (see $extendBidirectedPaths_*()$ for 
	 * details). The strings that label each extended path are stored in 
	 * $extendedPathsLabels$. For every extended path, a row of $path2path$ tells whether
	 * the extended path should be concatenated with another extended path, in the 
	 * following format:
	 *
	 * pathLeft,overlapLeft,pathRight,overlapRight
	 *
	 * where $path*$ is the ID of a path ($x$ if path $x$ should be concatenated from its
	 * prefix, $-1-x$ if $x$ should be concatenated from its suffix), and $overlap*$ is 
	 * the number of characters of $path*$ that already belong to the current path. 
	 * $path*$ is $Math.POSITIVE_INFINITY$ if no concatenation should be done.
	 *
	 * @param maxLength max allowed string length of an extended path;
	 * @param maxNodes upper bound on the number of nodes in the ring cluster;
	 * @param nodeLabels the label of every interval graph node;
	 * @param pathEndExtended,bidirectedGraphStrings temporary space, of size at least
	 * equal to the number of paths in the bidirected graph;
	 * @param tmpStream* temporary space.
	 */
	private static final void extendBidirectedPaths(int ballID, int maxLength, int maxNodes, CharStream[] nodeLabels, CharStream[] extendedPathsLabels, int[][] path2path, boolean[][] pathEndExtended, CharStream[] bidirectedGraphStrings, CharStream tmpStream1, CharStream tmpStream2, CharStream tmpStream3) throws IOException {
		boolean nodeOrientation;
		int i, j;
		int node, nodeEnd, nodesInPath, last, length;
		final int lastPath = BidirectedGraph.lastAssemblyPath;
		IntervalGraph.Node tmpNode;


StringBuilder sb = new StringBuilder(10000);

		
		// Initializing data structures
		if (tmpArray.length<maxNodes*3) tmpArray = new int[maxNodes*3];
		Math.set(pathEndExtended,lastPath,false);
		Math.set(path2path,lastPath,Math.POSITIVE_INFINITY);
		for (i=0; i<=lastPath; i++) {
			node=BidirectedGraph.assemblyPaths[i][0];
			BidirectedGraph.intervalGraphPointers[node>=0?node:-1-node].isMaximal=true;
			node=BidirectedGraph.assemblyPaths[i][BidirectedGraph.assemblyPathLengths[i]-1];
			BidirectedGraph.intervalGraphPointers[node>=0?node:-1-node].isMaximal=true;
		}
		for (i=0; i<BidirectedGraph.nNodes; i++) bidirectedGraphStrings[i] = nodeLabels[BidirectedGraph.intervalGraphPointers[i].nodeID];
		for (i=0; i<=lastPath; i++) extendedPathsLabels[i].clear(false);
		
		// Extending
		for (i=0; i<=lastPath; i++) {
			nodesInPath=BidirectedGraph.assemblyPathLengths[i];
			for (j=0; j<nodesInPath; j++) {
				node=BidirectedGraph.assemblyPaths[i][j];
				if (node<0) node=-1-node;
				BidirectedGraph.intervalGraphPointers[node].isSingletonComponent=true;
			}
			if (nodesInPath==1) {
				tmpStream2.clear(false);
				tmpNode=BidirectedGraph.intervalGraphPointers[BidirectedGraph.assemblyPaths[i][0]];
				node=tmpNode.nodeID; nodeOrientation=true;
				length=maxLength-tmpNode.length();
				if (!pathEndExtended[i][0] && !pathEndExtended[i][1]) {
					last=extendBidirectedPaths_undirected(node,ballID,length,tmpArray);
					System.err.print("extendBidirectedPaths_undirected(): ");
					for (int x=0; x<=last; x++) System.err.print(tmpArray[x]+",");
					System.err.println();
					if (last!=-1) {
						getSequence_undirected(node,tmpArray,last,nodeLabels,tmpStream2);
						node=tmpArray[last-2];
						if (node<0) {
							node=-1-node;
							nodeOrientation=false;
						}
						else nodeOrientation=true;
						for (j=0; j<last; j+=3) length-=tmpArray[j+1]+tmpArray[j+2];
					}
					else tmpStream2.push(nodeLabels[node],0,tmpNode.length()-1,0);
				}
				tmpStream1.clear(false);
				if (!pathEndExtended[i][0]) {
					last=extendBidirectedPaths_directed(node,nodeOrientation?0:1,ballID,length,tmpArray,pathEndExtended,path2path[i],true);					
					System.err.print("extendBidirectedPaths_directed(), left: ");
					for (int x=0; x<=last; x++) System.err.print(tmpArray[x]+",");
					System.err.println();
					if (last!=-1) {
						for (j=0; j<last; j+=2) length-=tmpArray[j+1];
						IntervalGraph.nodesArray[tmpArray[last-1]>=0?tmpArray[last-1]:-1-tmpArray[last-1]].isMinimal=true;
						getSequence_directed(node,nodeOrientation,nodeOrientation?0:1,false,tmpArray,last,nodeLabels,tmpStream1);
					}
					else IntervalGraph.nodesArray[node].isMinimal=true;
					pathEndExtended[i][0]=true;
				}
				tmpStream3.clear(false);
				if (!pathEndExtended[i][1]) {
					last=extendBidirectedPaths_directed(node,nodeOrientation?1:0,ballID,length,tmpArray,pathEndExtended,path2path[i],false);
					System.err.print("extendBidirectedPaths_directed(), right: ");
					for (int x=0; x<=last; x++) System.err.print(tmpArray[x]+",");
					System.err.println();
					if (last!=-1) {
						for (j=0; j<last; j+=2) length-=tmpArray[j+1];
						IntervalGraph.nodesArray[tmpArray[last-1]>=0?tmpArray[last-1]:-1-tmpArray[last-1]].isMinimal=true;
						getSequence_directed(node,nodeOrientation,nodeOrientation?1:0,false,tmpArray,last,nodeLabels,tmpStream3);
					}
					else IntervalGraph.nodesArray[node].isMinimal=true;
					pathEndExtended[i][1]=true;
				}
				extendedPathsLabels[i].clear(false);
				if (tmpStream1.nCharacters()>0) extendedPathsLabels[i].push(tmpStream1,0,tmpStream1.nCharacters()-1,0);
				if (tmpStream2.nCharacters()>0) extendedPathsLabels[i].push(tmpStream2,0,tmpStream2.nCharacters()-1,0);
				if (tmpStream3.nCharacters()>0) extendedPathsLabels[i].push(tmpStream3,0,tmpStream3.nCharacters()-1,0);
				
				
System.out.println("Ball "+ballID+" :: Bidirected path "+i+" (singleton) :: "+BidirectedGraph.intervalGraphPointers[BidirectedGraph.assemblyPaths[i][0]].length()+" -> "+(tmpStream1.nCharacters()+tmpStream2.nCharacters()+tmpStream3.nCharacters())+" :: pathEndExtended="+pathEndExtended[i][0]+","+pathEndExtended[i][1]);
sb.delete(0,sb.length());
if (tmpStream1.nCharacters()>0) tmpStream1.toString(sb);
if (tmpStream2.nCharacters()>0) tmpStream2.toString(sb);
if (tmpStream3.nCharacters()>0) tmpStream3.toString(sb);
System.out.println(sb);
					
					
					
			}
			else {
				length=maxLength-BidirectedGraph.assemblyPathStringLengths[i];
				tmpStream1.clear(false);
				if (!pathEndExtended[i][0]) {
					j=BidirectedGraph.assemblyPaths[i][0];
					if (j>=0) {
						nodeEnd=0;
						nodeOrientation=true;
					}
					else {
						j=-1-j;
						nodeEnd=1;
						nodeOrientation=false;
					}
					node=BidirectedGraph.intervalGraphPointers[j].nodeID;
					last=extendBidirectedPaths_directed(node,nodeEnd,ballID,length,tmpArray,pathEndExtended,path2path[i],true);					
					System.err.print("extendBidirectedPaths_directed(), left: ");
					for (int x=0; x<=last; x++) System.err.print(tmpArray[x]+",");
					System.err.println();
					if (last!=-1) {
						for (j=0; j<last; j+=2) length-=tmpArray[j+1];
						IntervalGraph.nodesArray[tmpArray[last-1]>=0?tmpArray[last-1]:-1-tmpArray[last-1]].isMinimal=true;
						getSequence_directed(node,nodeOrientation,nodeEnd,false,tmpArray,last,nodeLabels,tmpStream1);
					}
					else IntervalGraph.nodesArray[node].isMinimal=true;
					pathEndExtended[i][0]=true;
				}
				tmpStream3.clear(false);
				if (!pathEndExtended[i][1]) {
					j=BidirectedGraph.assemblyPaths[i][nodesInPath-1];
					if (j>=0) {
						nodeEnd=1;
						nodeOrientation=true;
					}
					else {
						j=-1-j;
						nodeEnd=0;
						nodeOrientation=false;
					}
					node=BidirectedGraph.intervalGraphPointers[j].nodeID;
					last=extendBidirectedPaths_directed(node,nodeEnd,ballID,length,tmpArray,pathEndExtended,path2path[i],false);
					System.err.print("extendBidirectedPaths_directed(), right: ");
					for (int x=0; x<=last; x++) System.err.print(tmpArray[x]+",");
					System.err.println();
					if (last!=-1) {
						for (j=0; j<last; j+=2) length-=tmpArray[j+1];
						IntervalGraph.nodesArray[tmpArray[last-1]>=0?tmpArray[last-1]:-1-tmpArray[last-1]].isMinimal=true;
						getSequence_directed(node,nodeOrientation,nodeEnd,false,tmpArray,last,nodeLabels,tmpStream3);
					}
					else IntervalGraph.nodesArray[node].isMinimal=true;
					pathEndExtended[i][1]=true;
				}
				BidirectedGraph.printPath(i,bidirectedGraphStrings,null,null,tmpStream2,0,-1);
				extendedPathsLabels[i].clear(false);
				if (tmpStream1.nCharacters()>0) extendedPathsLabels[i].push(tmpStream1,0,tmpStream1.nCharacters()-1,0);
				if (tmpStream2.nCharacters()>0) extendedPathsLabels[i].push(tmpStream2,0,tmpStream2.nCharacters()-1,0);
				if (tmpStream3.nCharacters()>0) extendedPathsLabels[i].push(tmpStream3,0,tmpStream3.nCharacters()-1,0);
			
			


System.out.println("Ball "+ballID+" :: Bidirected path "+i+" :: "+BidirectedGraph.assemblyPathStringLengths[i]+" -> "+(tmpStream1.nCharacters()+tmpStream2.nCharacters()+tmpStream3.nCharacters())+" :: pathEndExtended="+pathEndExtended[i][0]+","+pathEndExtended[i][1]);
sb.delete(0,sb.length());
if (tmpStream1.nCharacters()>0) tmpStream1.toString(sb);
if (tmpStream2.nCharacters()>0) tmpStream2.toString(sb);
if (tmpStream3.nCharacters()>0) tmpStream3.toString(sb);
System.out.println(sb);




				
			}
			for (j=0; j<nodesInPath; j++) {
				node=BidirectedGraph.assemblyPaths[i][j];
				if (node<0) node=-1-node;
				BidirectedGraph.intervalGraphPointers[node].isSingletonComponent=false;
			}
		}
	}
	
	
	/**
	 * Stores in $out$ the sequence of nodes obtained by extending $node$ greedily by 
	 * error rate, from its end $nodeEnd$, using just overlap edges or prefix/suffix
	 * containment edges. The procedure stops if the node of a bidirected path, or a node 
	 * at the end of an extended bidirected path, is met.
	 *
	 * Remark: greedy is motivated by the fact that the sequence of the initial node is 
	 * assumed to be close to the center of a ring. Greedy extension is not likely to 
	 * merge paths that belong to the centers of distinct (highly overlapping) rings.
	 *
	 * @param node assumed to be the first/last node of a bidirected path that has yet to
	 * be extended;
	 * @param nodeEnd 0=prefix, 1=suffix;
	 * @param maxLength max number of characters that can be added to $node$;
	 * @param out output array: list of pairs $(nodeID,overhang)$, where $overhang$ is the
	 * number of extra characters added by $nodeID$, and $nodeID$ is negative if it should
	 * be taken in reverse orientation, positive otherwise; $out$ is assumed to be big
	 * enough to contain the output;
   	 * @param pathEndExtended flags indicating that the first/last node of a bidirected
   	 * path has already been extended; the procedure might set one entry to TRUE;
	 * @param path2path row corresponding to the bidirected path to which $node$ belongs;
	 * @param path2pathSide tells whether $node$ is at the beginning (true) or at the end 
	 * (false) of its bidirected path;
	 * @return the last element of $out$.
	 */
	private static final int extendBidirectedPaths_directed(int node, int nodeEnd, int ballID, int maxLength, int[] out, boolean[][] pathEndExtended, int[] path2path, boolean path2pathSide) {
		boolean minType;  // TRUE=overlap; FALSE=contaiment.
		int i;
		int to, last, nNeighbors, minEdge, overhang, length, toEnd, pointer;
		double minDiffs;
		IntervalGraph.Node toNode;
		IntervalGraph.Edge edge;
		
		last=-1; length=0;
		while (true) {
			nNeighbors=IntervalGraph.nNeighbors[node];
			minDiffs=Math.POSITIVE_INFINITY; minEdge=-1; minType=false;
			for (i=0; i<nNeighbors; i++) {
				edge=IntervalGraph.neighbors[node][i];
				if (!edge.on) break;
				if (edge.orientation==-1 || edge.orientation==2 || edge.nTypes()>1 || edge.avgDiffs>=minDiffs) continue;
				to=edge.getTo(node);
				if (IntervalGraph.nodesArray[to].iteration!=ballID) continue;
				if (edge.getOverlapEnd(node)==nodeEnd) {
					minDiffs=edge.avgDiffs;
					minEdge=i;
					minType=true;
				}
				else if ( ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==node) ||
					  	    (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID2==node)
					      ) && 
				          ( ( nodeEnd==0 &&
						      ( (edge.orientation==0 && edge.containmentSubtype==Constants.CONTAINMENT_SUFFIX) ||
						  	    (edge.orientation==1 && edge.containmentSubtype==Constants.CONTAINMENT_PREFIX)
						      )
						    ) ||
						    ( nodeEnd==1 &&
						      ( (edge.orientation==0 && edge.containmentSubtype==Constants.CONTAINMENT_PREFIX) ||
							    (edge.orientation==1 && edge.containmentSubtype==Constants.CONTAINMENT_SUFFIX)
						      )
						    )
					      )
					    ) {
					minDiffs=edge.avgDiffs;
					minEdge=i;
					minType=false;
				}
			}
			if (minEdge==-1) break;
			edge=IntervalGraph.neighbors[node][minEdge];
			if (minType) {
				overhang=Math.abs(edge.getOverhang(node,nodeEnd==0));
if (overhang<0) {
	System.err.println("extendBidirectedPaths_directed> 1  negative overhang: "+overhang+" :: edge.getOverhang()="+edge.getOverhang(node,nodeEnd==0));
	System.err.println("extendBidirectedPaths_directed> edge: "+edge);
	System.exit(1);
}
			}
			else {
				overhang=edge.containmentSubtype==Constants.CONTAINMENT_PREFIX?edge.containmentOverhangRight:edge.containmentOverhangLeft;
if (overhang<0) {
	System.err.println("extendBidirectedPaths_directed> 2  negative overhang: "+overhang);
	System.err.println("extendBidirectedPaths_directed> edge: "+edge);
	System.exit(1);
}
			}
			to=edge.getTo(node);
			toNode=IntervalGraph.nodesArray[to];
			toEnd=edge.getOverlapEnd(to);
			if (toNode.isSingletonComponent/*=in the current bidirected path*/|| toNode.isMinimal/*=end of some extended path*/|| length+overhang>maxLength) break;
			if (toNode.isMaximal/*=end of a bidirected path*/) {
				pointer=markPathEnd(toNode,toEnd,pathEndExtended);
				if (pointer!=Math.POSITIVE_INFINITY) {
					if (edge.orientation==0) out[++last]=to;
					else {
						out[++last]=-1-to;
						nodeEnd=1-nodeEnd;
					}
					out[++last]=overhang;
					length+=overhang;
					if (path2pathSide) {
						path2path[0]=pointer;
						path2path[1]=toNode.length()-overhang;
					}
					else {
						path2path[2]=pointer;
						path2path[3]=toNode.length()-overhang;
					}
					toNode.frontier=true;
				}
				break;
			}
			if (edge.orientation==0) out[++last]=to;
			else {
				out[++last]=-1-to;
				nodeEnd=1-nodeEnd;
			}
			out[++last]=overhang;
			length+=overhang;
			toNode.frontier=true;
			node=to;
		}
		return last;
	}
	
	
	/**
	 * Variant that checks, at each step, if the greedy choice is a strict containment
	 * (extending in any direction, and possibly in both directions), and that stops as 
	 * soon as the greedy choice is not a containment.
	 * 
	 * @param maxLength max number of characters that can be added to $node$;
	 * @param out output array: list of triplets $(nodeID,overhangLeft,overhangRight)$, 
	 * where $overhang*$ is the number of extra characters added by $nodeID$, and $nodeID$
	 * is negative if it should be taken in reverse orientation, positive otherwise;
	 * @return the last element of $out$.
	 */
	private static final int extendBidirectedPaths_undirected(int node, int ballID, int maxLength, int[] out) {
		boolean minType;  // TRUE=overlap; FALSE=contaiment.
		int i;
		int to, last, nNeighbors, minEdge, length;
		double minDiffs;
		IntervalGraph.Node toNode;
		IntervalGraph.Edge edge;
		
		last=-1; length=0;
		while (true) {
			nNeighbors=IntervalGraph.nNeighbors[node];
			minDiffs=Math.POSITIVE_INFINITY; minEdge=-1; minType=false;
			for (i=0; i<nNeighbors; i++) {
				edge=IntervalGraph.neighbors[node][i];
				if (!edge.on) break;
				if (edge.orientation==-1 || edge.orientation==2 || edge.nTypes()>1 || edge.avgDiffs>=minDiffs) continue;
				to=edge.getTo(node);
				if (IntervalGraph.nodesArray[to].iteration!=ballID) continue;
				if (edge.isType_overlap()) {
					minDiffs=edge.avgDiffs;
					minEdge=i;
					minType=true;
				}
				else if ( (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO && edge.nodeID1==node) ||
					  	  (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE && edge.nodeID2==node)
					    ) {
					minDiffs=edge.avgDiffs;
					minEdge=i;
					minType=false;
				}
			}
			if (minEdge==-1 || minType) break;
			edge=IntervalGraph.neighbors[node][minEdge];
			to=edge.getTo(node);
			if (to==node) break;
			toNode=IntervalGraph.nodesArray[to];
			if (toNode.isMaximal/*=end of a bidirected path*/|| toNode.isMinimal/*=end of some extended path*/| length+edge.containmentOverhangLeft+edge.containmentOverhangRight>maxLength) break;
			out[++last]=edge.orientation==0?to:-1-to;
			out[++last]=edge.containmentOverhangLeft;
			out[++last]=edge.containmentOverhangRight;
			toNode.frontier=true;
			node=to;
		}
		return last;
	}
	
	
	/**
	 * Sets to TRUE the entry of $pathEndExtended$ that corresponds to $node$ taken from
	 * the end that is opposite to $nodeEnd$, if any.
	 * 
	 * @param node an end of a bidirected path;
	 * @param nodeEnd end of $node$ reached by the caller (0=prefix, 1=suffix);
	 * @return $Math.POSITIVE_INFINITY$ if no such entry of $pathEndExtended$ exists, or 
	 * if it exists but it is already set to TRUE; otherwise, if the entry belongs to path
	 * $i$: value $i$ if it is the beginning of the path, $-1-i$ if it is the end.
	 */
	private static final int markPathEnd(IntervalGraph.Node node, int nodeEnd, boolean[][] pathEndExtended) {
		boolean bidirectedOrientation;
		int i;
		int bidirectedNode;
		final int lastPath = BidirectedGraph.lastAssemblyPath;
		
		for (i=0; i<=lastPath; i++) {
			// First node
			bidirectedNode=BidirectedGraph.assemblyPaths[i][0];
			if (bidirectedNode>=0) bidirectedOrientation=true;
			else {
				bidirectedNode=-1-bidirectedNode;
				bidirectedOrientation=false;
			}
			if (BidirectedGraph.intervalGraphPointers[bidirectedNode]==node) {
				if ((bidirectedOrientation && nodeEnd==0) || (!bidirectedOrientation && nodeEnd==1)) {
					if (pathEndExtended[i][0]) return Math.POSITIVE_INFINITY;
					else {
						pathEndExtended[i][0]=true;
						return i;
					}
				}
			}
			// Last node
			bidirectedNode=BidirectedGraph.assemblyPaths[i][BidirectedGraph.assemblyPathLengths[i]-1];
			if (bidirectedNode>=0) bidirectedOrientation=true;
			else {
				bidirectedNode=-1-bidirectedNode;
				bidirectedOrientation=false;
			}
			if (BidirectedGraph.intervalGraphPointers[bidirectedNode]==node) {
				if ((bidirectedOrientation && nodeEnd==1) || (!bidirectedOrientation && nodeEnd==0)) {
					if (pathEndExtended[i][1]) return Math.POSITIVE_INFINITY;
					else {
						pathEndExtended[i][1]=true;
						return -1-i;
					}
				}
			}
		}
		return Math.POSITIVE_INFINITY;
	}
	
	
	/**
	 * @return the string label of every node of the interval graph.
	 */
	private static final CharStream[] loadNodeLabels(String path) throws IOException {
		final int CAPACITY = 20000;  // Arbitrary
		final int nNodes = IntervalGraph.nNodes;
		final int nReads = Reads.nReads;
		int i;
		int currentFragment, currentLength;
		String str;
		IntervalGraph.Node currentNode;
		CharStream tmpStream;
		BufferedReader br;
		CharStream[] nodeLabels;
		
		nodeLabels = new CharStream[nNodes];
		for (i=0; i<nNodes; i++) nodeLabels[i] = new CharStream(Reads.getReadLength(IntervalGraph.nodesArray[i].read));
		tmpStream = new CharStream(CAPACITY);
		i=0;
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();  // Header
		currentFragment=0; currentLength=Reads.getReadLength(0);
		while (str!=null && i<nNodes) {
			tmpStream.ensureSize(currentLength);
			tmpStream.load(currentLength,br);
			currentNode=IntervalGraph.nodesArray[i];
			if (currentNode.read==currentFragment) nodeLabels[i++].push(tmpStream,currentNode.start,currentNode.end,0);
			str=br.readLine();  // Next header
			currentFragment++;
			if (currentFragment<nReads) currentLength=Reads.getReadLength(currentFragment);
		}
		br.close(); br=null;
		tmpStream.deallocate(); tmpStream=null;
		return nodeLabels;
	}
	
	
	/**
	 * Appends to $stream$ the characters from the sequence of nodes generated by 
	 * $extendBidirectedPaths_directed()$.
	 *
	 * @param extendedNode initial node from which extension started;
	 * @param extendedNodeOrientation orientation of $extendedNode$ WRT a path that is 
	 * implicit from the context;
	 * @param extendedNodeEnd 0=prefix, 1=suffix;
	 * @param useExtendedNode if TRUE, the procedure adds the sequence of $extendedNode$
	 * to $stream$;
	 * @param stream the initial characters in the stream are deleted.
	 */
	private static final void getSequence_directed(int extendedNode, boolean extendedNodeOrientation, int extendedNodeEnd, boolean useExtendedNode, int[] array, int last, CharStream[] nodeLabels, CharStream stream) {
		int i;
		int node, nodeLength, length;
		
		// Allocating
		length=useExtendedNode?IntervalGraph.nodesArray[extendedNode].length():0;
		for (i=0; i<last; i+=2) length+=array[i+1];
		stream.ensureSize(length);
		stream.clear(false);
		
		if (extendedNodeOrientation) {
			if (extendedNodeEnd==0) {
				// Left overhangs
				for (i=last-1; i>=0; i-=2) {
					node=array[i];
					if (node>=0) stream.push(nodeLabels[node],0,array[i+1]-1,0);
					else {
						node=-1-node;
						nodeLength=IntervalGraph.nodesArray[node].length();
						stream.push(nodeLabels[node],nodeLength-array[i+1],nodeLength-1,1);
					}
				}
				if (useExtendedNode) {
					nodeLength=IntervalGraph.nodesArray[extendedNode].length();
					stream.push(nodeLabels[extendedNode],0,nodeLength-1,0);
				}
			}
			else {
				// Right overhangs
				if (useExtendedNode) {
					nodeLength=IntervalGraph.nodesArray[extendedNode].length();
					stream.push(nodeLabels[extendedNode],0,nodeLength-1,0);
				}
				for (i=0; i<last; i+=2) {
					node=array[i];
					if (node>=0) {
						nodeLength=IntervalGraph.nodesArray[node].length();
						stream.push(nodeLabels[node],nodeLength-array[i+1],nodeLength-1,0);
					}
					else {
						node=-1-node;
						stream.push(nodeLabels[node],0,array[i+1]-1,1);
					}
				}
			}
		}
		else {
			if (extendedNodeEnd==1) {
				// Left overhangs
				for (i=last-1; i>=0; i-=2) {
					node=array[i];
					if (node>=0) {
						nodeLength=IntervalGraph.nodesArray[node].length();
						stream.push(nodeLabels[node],nodeLength-array[i+1],nodeLength-1,1);
					}
					else {
						node=-1-node;
						stream.push(nodeLabels[node],0,array[i+1]-1,0);
					}
				}
				if (useExtendedNode) {
					nodeLength=IntervalGraph.nodesArray[extendedNode].length();
					stream.push(nodeLabels[extendedNode],0,nodeLength-1,1);
				}
			}
			else {
				// Right overhangs
				if (useExtendedNode) {
					nodeLength=IntervalGraph.nodesArray[extendedNode].length();
					stream.push(nodeLabels[extendedNode],0,nodeLength-1,1);
				}
				for (i=0; i<last; i+=2) {
					node=array[i];
					if (node>=0) stream.push(nodeLabels[node],0,array[i+1]-1,1);
					else {
						node=-1-node;
						nodeLength=IntervalGraph.nodesArray[node].length();
						stream.push(nodeLabels[node],nodeLength-array[i+1],nodeLength-1,0);
					}
				}
			}
		}
	}
	
	
	/**
	 * Stores in $stream$ the characters from the sequence of nodes generated by 
	 * $extendBidirectedPaths_undirected()$.
	 *
	 * @param extendedNode initial node from which extension started, assumed to be in the 
	 * forward orientation; its sequence is included in the output;
	 * @param stream the initial characters in the stream are deleted.
	 */
	private static final void getSequence_undirected(int extendedNode, int[] array, int last, CharStream[] nodeLabels, CharStream stream) {
		int i;
		int node, nodeLength, length;
		
		// Allocating
		length=IntervalGraph.nodesArray[extendedNode].length();
		for (i=0; i<last; i+=3) length+=array[i+1]+array[i+2];
		stream.ensureSize(length);
		stream.clear(false);
		
		// Left overhangs
		for (i=last-2; i>=0; i-=3) {
			node=array[i];
			if (node>=0) {
				if (array[i+1]>0) stream.push(nodeLabels[node],0,array[i+1]-1,0);
			}
			else {
				node=-1-node;
				if (array[i+2]>0) {
					nodeLength=IntervalGraph.nodesArray[node].length();
					stream.push(nodeLabels[node],nodeLength-array[i+2],nodeLength-1,1);
				}
			}
		}
		// Initial node
		nodeLength=IntervalGraph.nodesArray[extendedNode].length();
		stream.push(nodeLabels[extendedNode],0,nodeLength-1,0);
		// Right overhangs
		for (i=0; i<last; i+=3) {
			node=array[i];
			if (node>=0) {
				if (array[i+2]>0) {
					nodeLength=IntervalGraph.nodesArray[node].length();
					stream.push(nodeLabels[node],nodeLength-array[i+2],nodeLength-1,0);
				}
			}
			else {
				node=-1-node;
				if (array[i+1]>0) stream.push(nodeLabels[node],0,array[i+1]-1,1);
			}
		}
	}
	
	
	/**
	 * Concatenates the extended paths built by $extendBidirectedPaths()$ as specified by
	 * $path2path$, and writes the resulting strings to $outputFile$. The procedure 
	 * assumes that every extended path belongs to at most one concatenation chain.
	 *
	 * Remark: one could try to concatenate bidirected paths more aggressively, e.g. one
	 * could build a bidirected graph of bidirected paths, using all edges in the interval
	 * graph, since a bidirected path might be just a small alternative of a full repeat.
	 * We don't do this, since we consider more prudent to extend each bidirected path
	 * greedily.
	 *
	 * Remark: some of the resulting extended paths might be very similar to one another;
	 * however, there is no good reason for discarding any of them at this stage, since 
	 * they are all needed to cover the specified fraction of fragments (e.g. the 
	 * fragments in the symmetric difference of their rings).
	 *
	 * @param concatenated temporary space, of size at least equal to the number of 
	 * extended paths;
	 * @param tmpStream temporary space, resized if necessary.
	 */
	private static final void concatenateExtendedPaths(CharStream[] extendedPathsLabels, int[][] path2path, String outputFile, String outputLengths, boolean[] concatenated, CharStream tmpStream) throws IOException {
		int i, j;
		int jSide, top, overlap, pathID, pathEnd, length;
		final int lastPath = BidirectedGraph.lastAssemblyPath;
		BufferedWriter bw1, bw2;
		
		if (tmpArray==null || tmpArray.length<(lastPath+1)<<1) tmpArray = new int[(lastPath+1)<<1];
		Math.set(concatenated,lastPath,false);
		bw1 = new BufferedWriter(new FileWriter(outputFile,true/*Appending*/));
		bw2 = new BufferedWriter(new FileWriter(outputLengths,true/*Appending*/));
		for (i=0; i<=lastPath; i++) {
			if (concatenated[i]) continue;
			tmpStream.clear(false);
			// Concatenating from the left of the current path
			top=-1; j=i; jSide=0;
			while (path2path[j][jSide]!=Math.POSITIVE_INFINITY) {
				tmpArray[++top]=path2path[j][jSide];
				tmpArray[++top]=path2path[j][jSide+1];
				System.err.println("concatenateExtendedPaths> 1  concatenating path "+path2path[j][jSide]+" to path "+i);
				j=path2path[j][jSide];
				if (j<0) {
					j=-1-j;
					jSide=0;
				}
				else jSide=2;
				concatenated[j]=true;
			}
			for (j=top-1; j>=0; j-=2) {
				if (tmpArray[j]>=0) {
					pathEnd=0;
					pathID=tmpArray[j];
				}
				else {
					pathEnd=1;
					pathID=-1-tmpArray[j];
				}
				overlap=tmpArray[j+1];
				if (pathEnd==1) tmpStream.push(extendedPathsLabels[pathID],0,extendedPathsLabels[pathID].nCharacters()-1-overlap,0);
				else tmpStream.push(extendedPathsLabels[pathID],overlap,extendedPathsLabels[pathID].nCharacters()-1,1);
			}
			// Concatenating the current path
			tmpStream.push(extendedPathsLabels[i],0,extendedPathsLabels[i].nCharacters()-1,0);
			concatenated[i]=true;
			// Concatenating from the right of the current path
			top=-1; j=i; jSide=2;
			while (path2path[j][jSide]!=Math.POSITIVE_INFINITY) {
				tmpArray[++top]=path2path[j][jSide];
				tmpArray[++top]=path2path[j][jSide+1];
				System.err.println("concatenateExtendedPaths> 2  concatenating path "+path2path[j][jSide]+" to path "+i);
				j=path2path[j][jSide];
				if (j<0) {
					j=-1-j;
					jSide=0;
				}
				else jSide=2;
				concatenated[j]=true;
			}
			for (j=0; j<top; j+=2) {
				if (tmpArray[j]>=0) {
					pathEnd=0;
					pathID=tmpArray[j];
				}
				else {
					pathEnd=1;
					pathID=-1-tmpArray[j];
				}
				overlap=tmpArray[j+1];
				if (pathEnd==0) tmpStream.push(extendedPathsLabels[pathID],overlap,extendedPathsLabels[pathID].nCharacters()-1,0);
				else tmpStream.push(extendedPathsLabels[pathID],0,extendedPathsLabels[pathID].nCharacters()-1-overlap,1);
			}
			// Writing to file
			length=tmpStream.nCharacters();
			IO.writeFakeHeader(header++/*Unlikely to collide with any meaningful header*/,length,null,bw1);			
			tmpStream.print(bw1,true,0);
			bw1.write('\n');
			bw2.write(length+"\n");
		}
		bw1.close(); bw2.close();
	}
	
	
	private static final void printPaths(CharStream[] pathLabels, int nPaths, String outputFile, String outputLengths) throws IOException {
		int i;
		int length;
		BufferedWriter bw1, bw2;
		
		bw1 = new BufferedWriter(new FileWriter(outputFile,true/*Appending*/));
		bw2 = new BufferedWriter(new FileWriter(outputLengths,true/*Appending*/));
		for (i=0; i<nPaths; i++) {
			length=pathLabels[i].nCharacters();
			IO.writeFakeHeader(header++/*Unlikely to collide with any meaningful header*/,length,null,bw1);
			pathLabels[i].print(bw1,true,0);
			bw1.write('\n');
			bw2.write(length+"\n");
		}
		bw1.close(); bw2.close();
	}
	
	
	
	
	
	
	
	
	// ----------------------------- PRINTING PROCEDURES ---------------------------------
	
	/**
	 * For every node of the interval graph that is connected to some dominator node and
	 * that is not used to build a new reference, the procedure writes its string label to 
	 * $newFragmentsFile$.
	 */
	private static final void printNewFragments(CharStream[] nodeLabels, String newFragmentsFile, String newLengthsFile) throws IOException {
		final int nNodes = IntervalGraph.nNodes;
		int i;
		IntervalGraph.Node tmpNode;
		BufferedWriter bw1, bw2;
		
		bw1 = new BufferedWriter(new FileWriter(newFragmentsFile));
		bw2 = new BufferedWriter(new FileWriter(newLengthsFile));
		for (i=0; i<nNodes; i++) {
			tmpNode=IntervalGraph.nodesArray[i];
			if (tmpNode.visited==WHITE_NODE/*Never touched by $connectedDominator()$*/ || tmpNode.visited==BLACK_NODE/*Dominator*/ || tmpNode.frontier/*Extension of a dominator*/) continue;
			IO.writeFakeHeader(header++/*Unlikely to collide with any meaningful header*/,tmpNode.length(),null,bw1);
			nodeLabels[i].print(bw1,true,0);
			bw1.write('\n');
			bw2.write(tmpNode.length()+"\n");
		}
		bw1.close(); bw2.close();
	}
	
	
	/**
	 * Node fields in the output:
	 *
	 * referenceBlock: ID that encodes all blocks of the reference to which the fragment 
	 *                 aligns (all alignments are used, not just those that cover the 
	 *                 whole fragment);
	 * component: connected component;
	 * isReference: 1 if the node belongs to the basin's reference;
	 * dominatorTag: type of node assigned by the dominator procedure.
	 *
	 * Remark: in the output, node IDs start from one.
	 */
	private static final void printNodes(int referenceLength, String fragmentsReferenceFile, BufferedWriter bw) throws IOException {
		final int N_REFERENCE_BLOCKS = 4;  // Arbitrary
		final int MIN_BLOCK_INTERSECTION = IO.quantum<<1;  // Arbitrary
		final int BLOCK_LENGTH = referenceLength/N_REFERENCE_BLOCKS;
		boolean currentReference;
		int i;
		int currentRead;
		String str;
		BufferedReader br;
		IntervalGraph.Node node, tmpNode;
		boolean[] printed = new boolean[Reads.nReads];
		boolean[] currentBlocks = new boolean[N_REFERENCE_BLOCKS];
		
		Math.set(printed,IntervalGraph.nNodes-1,false);
		
		// Printing fragments that appear in $fragmentsReferenceFile$.
		tmpNode = new IntervalGraph.Node();
		tmpNode.start=-1; tmpNode.end=-1;
		br = new BufferedReader(new FileReader(fragmentsReferenceFile),IO.BUFFER_SIZE);
		currentRead=-1; currentReference=false;
		Math.set(currentBlocks,N_REFERENCE_BLOCKS-1,false);
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			if (Alignments.readA-1!=currentRead) {
				if (currentRead!=-1) {
					tmpNode.read=currentRead;
					i=-1-Arrays.binarySearch(IntervalGraph.nodesArray,0,IntervalGraph.nNodes,tmpNode);
					if (i<IntervalGraph.nNodes && IntervalGraph.nodesArray[i].read==currentRead) {
						node=IntervalGraph.nodesArray[i];
						bw.write((i+1)+" [read=\""+node.read+"\",start=\""+node.start+"\",end=\""+node.end+"\",referenceBlock=\""+blocks2id(currentBlocks,N_REFERENCE_BLOCKS)+"\",component=\""+(node.lastComponent==0?node.components[0]:-1)+"\",dominatorTag=\""+node.visited+"\",isReference=\""+(currentReference?"1":"0")+"\",extension=\""+(node.frontier?GREEDY_EXTENSION_NODE:node.visited)+"\"];\n");
						printed[i]=true;
					}
				}
				currentRead=Alignments.readA-1;
				currentReference=Alignments.diffs<=FragmentsStep1.TOO_FEW_DIFFS;
				Math.set(currentBlocks,N_REFERENCE_BLOCKS-1,false);
				printNodes_markBlocks(Alignments.startB,Alignments.endB,BLOCK_LENGTH,MIN_BLOCK_INTERSECTION,currentBlocks);
				str=br.readLine();
				continue;
			}
			currentReference|=Alignments.diffs<=FragmentsStep1.TOO_FEW_DIFFS;
			printNodes_markBlocks(Alignments.startB,Alignments.endB,BLOCK_LENGTH,MIN_BLOCK_INTERSECTION,currentBlocks);
			str=br.readLine();
		}
		br.close();
		if (currentRead!=-1) {
			tmpNode.read=currentRead;
			i=-1-Arrays.binarySearch(IntervalGraph.nodesArray,0,IntervalGraph.nNodes,tmpNode);
			if (i<IntervalGraph.nNodes && IntervalGraph.nodesArray[i].read==currentRead) {
				node=IntervalGraph.nodesArray[i];
				bw.write((i+1)+" [read=\""+node.read+"\",start=\""+node.start+"\",end=\""+node.end+"\",referenceBlock=\""+blocks2id(currentBlocks,N_REFERENCE_BLOCKS)+"\",component=\""+(node.lastComponent==0?node.components[0]:-1)+"\",dominatorTag=\""+node.visited+"\",isReference=\""+(currentReference?"1":"0")+"\",extension=\""+(node.frontier?GREEDY_EXTENSION_NODE:node.visited)+"\"];\n");
				printed[i]=true;
			}
		}
		
		// Printing remaining fragments
		for (i=0; i<IntervalGraph.nNodes; i++) {
			if (!printed[i]) {
				node=IntervalGraph.nodesArray[i];
				bw.write((i+1)+" [read=\""+node.read+"\",start=\""+node.start+"\",end=\""+node.end+"\",referenceBlock=\"0\",component=\""+(node.lastComponent==0?node.components[0]:-1)+"\",dominatorTag=\""+node.visited+"\",isReference=\"0\",extension=\""+(node.frontier?GREEDY_EXTENSION_NODE:node.visited)+"\"];\n");
			}
		}
	}
	
	
	private static final void printNodes_markBlocks(int start, int end, int blockLength, int minIntersection, boolean[] blocks) {
		if (blockLength==0) return;
		int firstBlock = start/blockLength;
		if (firstBlock>=blocks.length) return;
		if ((firstBlock+1)*blockLength<start+minIntersection) firstBlock++;
		int lastBlock = end/blockLength;
		if (lastBlock>=blocks.length) lastBlock=blocks.length-1;
		if (end<lastBlock*blockLength+minIntersection) lastBlock--;
		for (int i=firstBlock; i<=lastBlock; i++) blocks[i]=true;
	}
	
	
	/**
	 * Maps the configuration in $blocks$ to an integer identifier.
	 */
	private static final int blocks2id(boolean[] blocks, int nBlocks) {
		int mask = 0x000000001;
		int out = 0;
		for (int i=0; i<nBlocks; i++) {
			if (blocks[i]) out|=mask;
			mask<<=1;
		}
		return out;
	}
	
	
	/**
	 * Remark: the procedure prints ON edges only.
	 * Remark: node IDs start from one.
	 *
	 * @param dominatorOnly TRUE=prints only edges between nodes with $visited=2$.
	 */
	private static final void printEdges(boolean dominatorOnly, BufferedWriter bw) throws IOException {
		final double WIDTH_SCALE = 10;
		final String[] TYPE2COLOR = new String[] { 
			// Same as in $IntervalGraph$:
			// 0=containment, 1=overlap, 2=insertion, 3=sharedSubstring. (4=unknown)
			"yellow","green","red","blue","black" 
		};
		int i, j;
		int type;
		double identity;
		IntervalGraph.Edge edge;
		
		for (i=0; i<IntervalGraph.nNodes; i++) {
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
				if (!IntervalGraph.neighbors[i][j].on) break;
				IntervalGraph.neighbors[i][j].printed=-1;
			}
		}
		for (i=0; i<IntervalGraph.nNodes; i++) {
			if (dominatorOnly && IntervalGraph.nodesArray[i].visited!=BLACK_NODE) continue;
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
				edge=IntervalGraph.neighbors[i][j];
				if (!edge.on) break;
				if (edge.printed==1) continue;
				if (dominatorOnly && IntervalGraph.nodesArray[edge.getTo(i)].visited!=BLACK_NODE) continue;
				edge.printed=1;
				if (edge.nDistinctTypes()!=1) type=4;
				else {
					if (edge.isType_containment()) type=0;
					else if (edge.isType_overlap()) type=1;
					else if (edge.isType_insertion()) type=2;
					else if (edge.isType_sharedSubstring()) type=3;
					else type=4;
				}
				identity=1.0-edge.avgDiffs;
				if (type==0) {
					if (edge.containment==Constants.CONTAINMENT_ONE_IN_TWO) bw.write((edge.nodeID1+1)+" -> "+(edge.nodeID2+1));
					else if (edge.containment==Constants.CONTAINMENT_TWO_IN_ONE) bw.write((edge.nodeID2+1)+" -> "+(edge.nodeID1+1));
					else if (edge.containment==Constants.CONTAINMENT_IDENTICAL) bw.write((edge.nodeID1+1)+" -- "+(edge.nodeID2+1));
				}
				else bw.write((edge.nodeID1+1)+" -- "+(edge.nodeID2+1));
				bw.write(" [type=\""+type+"\",color=\""+TYPE2COLOR[type]+"\",weight=\""+(identity*WIDTH_SCALE)+"\"];\n");
			}
		}
		for (i=0; i<IntervalGraph.nNodes; i++) {
			for (j=0; j<IntervalGraph.nNeighbors[i]; j++) IntervalGraph.neighbors[i][j].printed=-1;
		}
	}
	
	
	/**
	 * Remark: node IDs start from one.
	 */
	private static final void printDominatorSubgraph(boolean onlyOverlap, BufferedWriter bw) throws IOException {
		final double NODE_SCALING = 0.01;
		final int nNodes = IntervalGraph.nNodes;
		int i, j;
		IntervalGraph.Node node;
		
		bw.write("digraph G {\n");
		for (i=0; i<lastDominator; i+=3) {
			node=IntervalGraph.nodesArray[dominator[i]];
			bw.write((dominator[i]+1)+" [iteration=\""+dominator[i+1]+"\",nDominated=\""+(100*((double)dominator[i+2])/nNodes)+"\",size=\""+(node.length()*NODE_SCALING)+"\"];\n");
		}
		if (onlyOverlap) {
			for (i=0; i<nNodes; i++) {
				for (j=0; j<IntervalGraph.nNeighbors[i]; j++) {
					if (!IntervalGraph.neighbors[i][j].isType_overlap()) IntervalGraph.neighbors[i][j].on=false;
				}
			}
		}
		printEdges(true,bw);
		bw.write("}\n");
	}

	
	/**
	 * (This procedure is unused)
	 *
	 * Builds $distanceMatrix$ (cell $(i,j)$ is the min error rate observed in any
	 * alignment between fragments $i,j$, of any alignment type) and $identityHistogram$ 
	 * (histogram of all observed identity values over all alignments).
	 */
	private static final void buildDistanceMatrix(String fragmentFragmentAlignmentsPath) throws IOException {
		final int N_HISTOGRAM_BINS = 100;
		final double MIN_IDENTITY = 0.6;
		final double IDENTITY_QUANTUM = (1.0-MIN_IDENTITY)/N_HISTOGRAM_BINS;
		int i, j;
		int lengthA, lengthB;
		double errorRate, identity;
		String str;
		BufferedReader alignmentsFile;
		int[] identityHistogram = new int[N_HISTOGRAM_BINS];
		double[][] distanceMatrix;
		
		distanceMatrix = new double[Reads.nReads][Reads.nReads];
		for (i=0; i<distanceMatrix.length; i++) {
			for (j=0; j<distanceMatrix[i].length; j++) distanceMatrix[i][j]=1.0;
			distanceMatrix[i][i]=0;
		}
		alignmentsFile = new BufferedReader(new FileReader(fragmentFragmentAlignmentsPath),IO.BUFFER_SIZE);
		alignmentsFile.readLine(); alignmentsFile.readLine();  // Skipping header
		str=alignmentsFile.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			lengthA=Reads.getReadLength(Alignments.readA-1);
			lengthB=Reads.getReadLength(Alignments.readB-1);
			errorRate=2*((double)Alignments.diffs)/(lengthA+lengthB);
			identity=1.0-errorRate;
			if (identity<MIN_IDENTITY) {
				System.err.println("buildMatrix> ERROR: too low identity "+identity+": "+str);
				System.exit(1);
			}
			if (errorRate<distanceMatrix[Alignments.readA-1][Alignments.readB-1]) {
				distanceMatrix[Alignments.readA-1][Alignments.readB-1]=errorRate;
				distanceMatrix[Alignments.readB-1][Alignments.readA-1]=errorRate;
			}
			identityHistogram[(int)Math.floor((identity-MIN_IDENTITY)/IDENTITY_QUANTUM)]++;			
			str=alignmentsFile.readLine();
		}
		alignmentsFile.close();
	}


}