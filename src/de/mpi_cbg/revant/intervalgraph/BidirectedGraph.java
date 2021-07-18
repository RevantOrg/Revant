package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.util.PriorityQueue;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.CharStream;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.DAG;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Colors;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Factorize;


public class BidirectedGraph {
	/**
	 * The bidirected graph
	 */
	private static final int NEIGHBORS_UNIT = 10;  // Arbitrarily chosen
	public static BidirectedEdge[][] neighbors;
	public static int[] lastNeighbor;
	public static int[] nodeLength;
	public static IntervalGraph.Node[] intervalGraphPointers;
	public static int nNodes;
	public static boolean[] removed;  // TRUE: the node has been removed.
	public static int nRemoved;
	public static int maxNodes, maxEdges, maxDegree;
	
	/**
	 * All assembly paths in the bidirected graph
	 */
	public static int lastAssemblyPath;
	public static int[][] assemblyPaths;
	public static int[] assemblyPathLengths;
	public static int[] assemblyPathStringLengths;
	public static int[][] node2paths;
	public static int[] node2paths_last;
	
	/**
	 * Temporary space for bidirected graph procedures
	 */
	private static boolean[][] isActive;
	public static int[] oneEnd;
	public static int nOneEnd;
	private static int[] reachedFrom;
	private static int[] stack;
	private static BidirectedEdge[] edgePool;
	private static PriorityEnd[] priorityEndPool;
	private static int lastUsedEdge;
	private static int lastPriorityEnd;
	public static int[] tmp;
	public static int[][] shortestDistance;
	public static PriorityQueue<PriorityEnd> queue;
	private static BidirectedEdge[][] newNeighbors;
	private static int[] newLastNeighbor;
	public static boolean[] oneEndFlags;
	private static CharStream charStack;
	private static boolean[] allPathsPrintedFrom;
	private static CharStream[] strings;
	
	/**
	 * Temporary space for DAG procedures
	 */
	private static int[][] original2dag;
	private static int[] dag2original;
	private static int[][] inNeighbors, outNeighbors;
	private static int[] lastInNeighbor, lastInNeighborPrime, lastOutNeighbor;
	private static int[] components;
	private static int[] minimalVertices;
	private static int[] nVerticesPerComponent;
	private static int[] sorted2original, original2sorted;
	private static int[] nPaths;
	private static int[][] inOverhangs;
	private static int[][] pathLengths;
	
	
	// ------------------------- BASIC PROCEDURES ON THE GRAPH ---------------------------
	
	/**
	 * @param maxN max number of nodes in a bidirected graph;
	 * @param maxE max number of edges in a bidirected graph (and in a DAG);
	 * @param maxD max degree in a bidirected graph (and in a DAG).
	 */
	public static final void allocateMemory(int maxN, int maxE, int maxD) {
		int i, maxDagNodes;
		
		// The bidirected graph
		maxNodes=maxN; maxEdges=maxE; maxDegree=maxD;
		neighbors = new BidirectedEdge[maxNodes][NEIGHBORS_UNIT];
		lastNeighbor = new int[maxNodes];
		Math.set(lastNeighbor,lastNeighbor.length-1,-1);
		newNeighbors = new BidirectedEdge[maxNodes][NEIGHBORS_UNIT];
		newLastNeighbor = new int[maxNodes];
		nNodes=0;
		removed = new boolean[maxNodes];
		nRemoved=0;
		nodeLength = new int[maxNodes];
		intervalGraphPointers = new IntervalGraph.Node[maxNodes];
		
		// Temporary space for bidirected graph procedures
		isActive = new boolean[maxNodes][2];
		oneEnd = new int[maxNodes];
		reachedFrom = new int[maxNodes];
		original2dag = new int[maxNodes][2];
		dag2original = new int[maxNodes<<1];
		stack = new int[maxNodes*3];  // *3 needed by $printAllPathsFrom()$.
		edgePool = new BidirectedEdge[maxEdges];  // Automatically enlarged if needed
		for (i=0; i<maxEdges; i++) edgePool[i] = new BidirectedEdge();
		lastUsedEdge=-1;
		tmp = new int[10];
		shortestDistance = new int[maxNodes][2];
		queue = new PriorityQueue<PriorityEnd>();
		priorityEndPool = new PriorityEnd[maxNodes];  // We don't set it to $maxNodes*2$ or $maxNodes*2+k$ for some constant $k$, because it is automatically enlarged.
		for (i=0; i<maxNodes; i++) priorityEndPool[i] = new PriorityEnd();
		lastPriorityEnd=-1;
		oneEndFlags = new boolean[maxNodes];
		
		// Temporary space for DAG procedures
		maxDagNodes=maxNodes<<1;
		inNeighbors = new int[maxDagNodes][NEIGHBORS_UNIT];
		lastInNeighbor = new int[maxDagNodes];
		lastInNeighborPrime = new int[maxDagNodes];
		outNeighbors = new int[maxDagNodes][NEIGHBORS_UNIT];
		lastOutNeighbor = new int[maxDagNodes];
		components = new int[maxDagNodes];
		minimalVertices = new int[maxDagNodes];
		nVerticesPerComponent = new int[1];  // Grows dynamically
		sorted2original = new int[maxDagNodes];
		original2sorted = new int[maxDagNodes];
		nPaths = new int[maxDagNodes];
		inOverhangs = new int[maxDagNodes][NEIGHBORS_UNIT];
		pathLengths = new int[maxDagNodes][2];
	}
	
	
	public static final void deallocateMemory() {
		int i;
		
		// The bidirected graph
		neighbors=null;
		newNeighbors=null;
		lastNeighbor=null;
		newLastNeighbor=null;
		nNodes=0;
		removed=null;
		nRemoved=0;
		nodeLength=null;
		for (i=0; i<intervalGraphPointers.length; i++) intervalGraphPointers[i]=null;
		intervalGraphPointers=null;
		
		// Temporary space for bidirected graph procedures
		isActive=null;
		oneEnd=null;
		reachedFrom=null;
		original2dag=null;
		stack=null;
		for (i=0; i<edgePool.length; i++) edgePool[i]=null;
		edgePool=null;
		for (i=0; i<priorityEndPool.length; i++) priorityEndPool[i]=null;
		priorityEndPool=null;
		
		// Temporary space for DAG procedures
		inNeighbors=null;
		lastInNeighbor=null;
		outNeighbors=null;
		lastOutNeighbor=null;
		components=null;
		minimalVertices=null;
		nVerticesPerComponent=null;
		sorted2original=null;
		original2sorted=null;
		nPaths=null;
	}
	
	
	public static final void clear() {
		int i;
		
		Math.set(lastNeighbor,lastNeighbor.length-1,-1);
		Math.set(nodeLength,nodeLength.length-1,0);
		Math.set(removed,removed.length-1,false);
		nRemoved=0;
		for (i=0; i<intervalGraphPointers.length; i++) intervalGraphPointers[i]=null;
		nNodes=0;
		nOneEnd=0;
		lastUsedEdge=-1;
		lastAssemblyPath=-1;
	}
	
	
	public static final int getMaxDegree() {
		int i, j, max;
		
		max=0;
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			j=lastNeighbor[i];
			if (j>max) max=j;
		}
		return max+1;
	}
	
	
	/**
	 * Sorts each array $neighbors[i]$ according to the one of $node1$ and $node2$ that is 
	 * different from $i$.
	 */
	public static final void sortNeighbors() {
		for (int i=0; i<nNodes; i++) {
			if (removed[i] || lastNeighbor[i]<=0) continue;
			BidirectedEdge.order=i;
			Arrays.sort(neighbors[i],0,lastNeighbor[i]+1);
		}
	}
	
	
	/**
	 * Gets an edge from the pool
	 */
	public static final BidirectedEdge getEdge() {
		lastUsedEdge++;
		if (lastUsedEdge==edgePool.length) {
			BidirectedEdge[] edgePoolPrime = new BidirectedEdge[edgePool.length<<1];
			System.arraycopy(edgePool,0,edgePoolPrime,0,lastUsedEdge);
			for (int i=lastUsedEdge; i<edgePoolPrime.length; i++) edgePoolPrime[i] = new BidirectedEdge();
			edgePool=edgePoolPrime;
		}
		return edgePool[lastUsedEdge];
	}
	
	
	public static final void resetEdgePool() {
		lastUsedEdge=-1;
	}
	
	
	/**
	 * Remark: the procedure assumes each edge to be stored among the neighbors of both of
	 * its incident nodes.
	 *
	 * @return -1 iff there is no edge, of any type, between nodes $x$ and $y$; the 
	 * position of the edge in $neighbors[x]$ otherwise.
	 */
	public static final int findEdge(int x, int y) {
		for (int i=0; i<=lastNeighbor[x]; i++) {
			if ( (neighbors[x][i].node1==x && neighbors[x][i].node2==y) ||
				 (neighbors[x][i].node2==x && neighbors[x][i].node1==y)
			   ) return i;
		}
		return -1;
	}
	
	
	/**
	 * @param xSide TRUE=the overlap involves a prefix of $x$; FALSE=the overlap involves 
	 * a suffix of $x$.
	 * @return TRUE iff there is an edge between $node$ and $neighbor$ that uses the 
	 * specified sides in the overlap.
	 */
	private static final boolean findEdge(int node, boolean nodeSide, int neighbor, boolean neighborSide) {
		final int last = lastNeighbor[node];
		int i;
		BidirectedEdge edge;
		
		for (i=0; i<=last; i++) {
			edge=neighbors[node][i];
			if (edge.node1==node && edge.node2==neighbor) {
				if (edge.type==Constants.sidesToOverlap(nodeSide,neighborSide)) return true;
				break;
			}
			if (edge.node1==neighbor && edge.node2==node) {
				if (edge.type==Constants.sidesToOverlap(neighborSide,nodeSide)) return true;
				break;
			}
		}
		return false;
	}
	
	
	public static final void addEdge(BidirectedEdge edge) {
		addEdge(edge,neighbors,lastNeighbor,nodeLength);
	}
	
	
	/**
	 * Adds to the input matrices a pointer to $edge$, without creating duplicates.
	 */
	private static final void addEdge(BidirectedEdge edge, BidirectedEdge[][] neighbors, int[] lastNeighbor, int[] nodeLength) {
		int foundX, foundY;
		int i;
		final int x = edge.node1;
		final int y = edge.node2;
		final int type = edge.type;
		BidirectedEdge[] tmpArray;
		
		// Searching the neighbors of $node1$ for $edge$.
		foundX=-1;
		for (i=0; i<=lastNeighbor[x]; i++) {
			if ( (neighbors[x][i].node1==x && neighbors[x][i].node2==y && neighbors[x][i].type==type) ||
				 (neighbors[x][i].node2==x && neighbors[x][i].node1==y && neighbors[x][i].type==Constants.oppositeDirectionOverlap(type))
			   ) {
			   foundX=i;
			   break;
			}
		}
		// Searching the neighbors of $node2$ for $edge$.
		foundY=-1;
		for (i=0; i<=lastNeighbor[y]; i++) {
			if ( (neighbors[y][i].node1==x && neighbors[y][i].node2==y && neighbors[y][i].type==type) ||
				 (neighbors[y][i].node2==x && neighbors[y][i].node1==y && neighbors[y][i].type==Constants.oppositeDirectionOverlap(type))
			   ) {
			   foundY=i;
			   break;
			}
		}
		// Adding $edge$
		if (foundX==-1) {
			if (foundY==-1) {
				lastNeighbor[x]++;
				if (lastNeighbor[x]==neighbors[x].length) {
					tmpArray = new BidirectedEdge[neighbors[x].length+NEIGHBORS_UNIT];
					System.arraycopy(neighbors[x],0,tmpArray,0,neighbors[x].length);
					neighbors[x]=tmpArray;
				}
				neighbors[x][lastNeighbor[x]]=edge;
				lastNeighbor[y]++;
				if (lastNeighbor[y]==neighbors[y].length) {
					tmpArray = new BidirectedEdge[neighbors[y].length+NEIGHBORS_UNIT];
					System.arraycopy(neighbors[y],0,tmpArray,0,neighbors[y].length);
					neighbors[y]=tmpArray;
				}
				neighbors[y][lastNeighbor[y]]=edge;
			}
			else {
				neighbors[y][foundY].copyOverhangs(edge,nodeLength);
				neighbors[y][foundY].nAlignments+=edge.nAlignments;
				lastNeighbor[x]++;
				if (lastNeighbor[x]==neighbors[x].length) {
					tmpArray = new BidirectedEdge[neighbors[x].length+NEIGHBORS_UNIT];
					System.arraycopy(neighbors[x],0,tmpArray,0,neighbors[x].length);
					neighbors[x]=tmpArray;
				}
				neighbors[x][lastNeighbor[x]]=neighbors[y][foundY];
			}
		}
		else {
			neighbors[x][foundX].copyOverhangs(edge,nodeLength);
			neighbors[x][foundX].nAlignments+=edge.nAlignments;
			if (foundY==-1) {
				lastNeighbor[y]++;
				if (lastNeighbor[y]==neighbors[y].length) {
					tmpArray = new BidirectedEdge[neighbors[y].length+NEIGHBORS_UNIT];
					System.arraycopy(neighbors[y],0,tmpArray,0,neighbors[y].length);
					neighbors[y]=tmpArray;
				}
				neighbors[y][lastNeighbor[y]]=neighbors[x][foundX];
			}
			else {
				if (neighbors[x][foundX]!=neighbors[y][foundY]) {
					System.err.println("ERROR: edge duplication in BidirectedGraph.addEdge()");
					System.exit(1);
				}
			}
		}
	}
	
	
	/**
	 * Sets to FALSE the $marked$ flag of all edges adjacent to a non-removed node (if 
	 * $skipRemoved=TRUE$), or of every edge (if $skipRemoved=FALSE$).
	 */
	public static final void cleanMarkedFlags(boolean skipRemoved) {
		int i, j;
		
		if (skipRemoved) {
			for (i=0; i<nNodes; i++) {
				if (removed[i]) continue;
				for (j=0; j<=lastNeighbor[i]; j++) neighbors[i][j].marked=false;
			}
		}
		else {
			for (i=0; i<nNodes; i++) {
				for (j=0; j<=lastNeighbor[i]; j++) neighbors[i][j].marked=false;
			}
		}
	}
	
	
	/**
	 * Deletes all marked edges incident to $node$, and removes the neighbors of $node$ 
	 * that get disconnected after such deletion.
	 *
	 * @param tmpArray temporary space, of size at least equal to the number of edges
	 * incident to $node$;
	 * @return the number of distinct nodes removed by the procedure.
	 */
	public static final int removeMarkedEdges(int node, int[] tmpArray) {
		int i, j, k;
		int neighbor, lastMarkedNeighbor, previous, out;
		BidirectedEdge tmpEdge;
		
		// Removing marked edges from $node$.
		j=-1; lastMarkedNeighbor=-1;
		for (i=0; i<=lastNeighbor[node]; i++) {
			if (!neighbors[node][i].marked) {
				j++;
				tmpEdge=neighbors[node][j];
				neighbors[node][j]=neighbors[node][i];
				neighbors[node][i]=tmpEdge;
			}
			else tmpArray[++lastMarkedNeighbor]=neighbors[node][i].getTo(node);
		}
		lastNeighbor[node]=j;
		if (lastMarkedNeighbor>0) {
			Arrays.sort(tmpArray,0,lastMarkedNeighbor+1);
			j=0; previous=tmpArray[0];
			for (k=1; k<=lastMarkedNeighbor; k++) {
				if (tmpArray[k]!=previous) {
					tmpArray[++j]=tmpArray[k];
					previous=tmpArray[k];
				}
			}
			lastMarkedNeighbor=j;
		}
		
		// Removing marked edges from the neighbors of $node$.
		for (k=0; k<=lastMarkedNeighbor; k++) {
			neighbor=tmpArray[k];
			j=-1;
			for (i=0; i<=lastNeighbor[neighbor]; i++) {
				if (!neighbors[neighbor][i].marked) {
					j++;
					tmpEdge=neighbors[neighbor][j];
					neighbors[neighbor][j]=neighbors[neighbor][i];
					neighbors[neighbor][i]=tmpEdge;
				}
			}
			lastNeighbor[neighbor]=j;
		}
		
		// Removing neighbors of $node$ that get disconnected.
		out=0;
		for (k=0; k<=lastMarkedNeighbor; k++) {
			neighbor=tmpArray[k];
			if (lastNeighbor[neighbor]==-1) {
				removeNode(neighbor);
				out++;
			}
		}
		return out;
	}
	
	
	/**
	 * Remark: the procedure preserves node IDs. It just marks removed nodes in $removed$
	 * and deletes every connection to them from $neighbors$.
	 */
	public static final void removeNode(int node) {
		removed[node]=true;
		disconnectNode(node);
		nRemoved++;
	}
	
	
	public static final void disconnectNode(int node) {
		int i, j;
		int to, neighbor, last;
		
		for (i=0; i<=lastNeighbor[node]; i++) {
			neighbor=neighbors[node][i].getTo(node);
			last=-1;
			for (j=0; j<=lastNeighbor[neighbor]; j++) {
				to=neighbors[neighbor][j].getTo(neighbor);
				if (to!=node) neighbors[neighbor][++last]=neighbors[neighbor][j];
			}
			lastNeighbor[neighbor]=last;
		}
		lastNeighbor[node]=-1;
	}
	
	
	/**
	 * @return TRUE iff the graph contains no edge.
	 */
	public static final boolean onlyDisconnectedNodes() {
		int i;
		
		if (nNodes==1) return true;
		for (i=0; i<nNodes; i++) {
			if (lastNeighbor[i]>=0) return false;
		}
		return true;
	}
	
	
	/**
	 * Marks as removed all nodes with no neighbor, unless the graph has just one node.
	 * 
	 * Remark: this procedure is currently never used.
	 */
	public static final void removeDisconnectedNodes() {
		int i;
		
		if (nNodes==1) return;
		for (i=0; i<nNodes; i++) {
			if (lastNeighbor[i]==-1) {
				removed[i]=true;
				nRemoved++;
			}
		}
	}

	
	/**
	 * Remark: the procedure assumes $neighbors$ to be sorted by destination node.
	 * Different edges might point to the same node.
	 *
	 * @return TRUE iff $node1$ and $node2$ are connected to the same nodes with the same
	 * (if $orientation=0,2$) or complementary (if $orientation=1$) egde types.
	 */
	public static final boolean sameNeighbors(int node1, int node2, int orientation, BidirectedEdge tmpEdge) {
		boolean found;
		int i, j, jPrime;
		int to, previousOrder;
		
		if (lastNeighbor[node1]!=lastNeighbor[node2]) return false;
		for (i=0; i<=lastNeighbor[node1]; i++) {
			tmpEdge.clone(neighbors[node1][i]);
			if (tmpEdge.node1==node1) tmpEdge.node1=node2;
			else tmpEdge.node2=node2;
			previousOrder=BidirectedEdge.order;
			BidirectedEdge.order=node2;
			j=Arrays.binarySearch(neighbors[node2],0,lastNeighbor[node2]+1,tmpEdge);
			BidirectedEdge.order=previousOrder;
			if (j<0) return false;
			if (sameType(node1,neighbors[node1][i],node2,neighbors[node2][j],orientation)) continue;
			to=neighbors[node1][i].getTo(node1);
			found=false;
			jPrime=j-1;
			while (jPrime>=0) {
				if (neighbors[node2][jPrime].getTo(node2)!=to) break;
				if (sameType(node1,neighbors[node1][i],node2,neighbors[node2][jPrime],orientation)) {
					found=true;
					break;
				}
				jPrime--;
			}
			if (found) continue;
			jPrime=j+1;
			while (jPrime<=lastNeighbor[node2]) {
				if (neighbors[node2][jPrime].getTo(node2)!=to) break;
				if (sameType(node1,neighbors[node1][i],node2,neighbors[node2][jPrime],orientation)) {
					found=true;
					break;
				}
				jPrime++;
			}
			if (!found) return false;
		}
		return true;
	}
	
	
	/**
	 * @return true iff $edge1$ and $edge2$ have the same type when seen from $from1$ and 
	 * $from2$, respectively (if $orientation=0,2$), or if the type of $edge2$ equals the 
	 * type of $edge1$ after swapping prefix with suffix in $from1$ (if $orientation=1$).
	 */
	private static final boolean sameType(int from1, BidirectedEdge edge1, int from2, BidirectedEdge edge2, int orientation) {
		final int type1 = edge1.node1==from1?edge1.type:Constants.oppositeDirectionOverlap(edge1.type);
		final int type2 = edge2.node1==from2?edge2.type:Constants.oppositeDirectionOverlap(edge2.type);
		if (orientation==0 || orientation==2) return type1==type2;
		else if (orientation==1) return type2==Constants.reverseComplementOverlap(type1,0);
		return false;
	}
	
	
	/**
	 * Remark: the procedure assumes that an edge incident to $node*$ is stored in both
	 * lists $neighbors[node*]$.
	 *
	 * @param overlapType one of the $OVERLAP_*$ codes.
	 */
	public static final boolean areOverlapping(int node1, int node2, int overlapType) {
		final int oppositeType = Constants.oppositeDirectionOverlap(overlapType);
		final int last = lastNeighbor[node1];
		int i;
		BidirectedEdge edge;
		
		for (i=0; i<=last; i++) {
			edge=neighbors[node1][i];
			if ( (edge.node1==node1 && edge.node2==node2 && edge.type==overlapType) ||
				 (edge.node1==node2 && edge.node2==node1 && edge.type==oppositeType)
			   ) return true;
		}
		return false;
	}
	
	
	/**
	 * @param neighbor* IDs in the bidirected graph of two neighbors of a fork;
	 * @param overhangSide* TRUE=the overhang with the fork is a prefix of the interval; 
	 * FALSE=the overhang with the fork is a suffix;
	 * @return TRUE iff the bidirected graph already contains an overlap, between the two
	 * neighbors of the fork, that involves the side of the overhang with the fork of one
	 * neighbor, and the side opposite to the overhang with the fork of the other 
	 * neighbor.
	 */
	public static final boolean neighborsOfForkOverlap(int neighbor1, boolean overhangSide1, int neighbor2, boolean overhangSide2) {
		if (overhangSide1) {
			if (overhangSide2) {
				if (areOverlapping(neighbor1,neighbor2,Constants.OVERLAP_PREFIX_SUFFIX) || areOverlapping(neighbor2,neighbor1,Constants.OVERLAP_PREFIX_SUFFIX)) return true;
			}
			else {
				if (areOverlapping(neighbor1,neighbor2,Constants.OVERLAP_PREFIX_PREFIX) || areOverlapping(neighbor2,neighbor1,Constants.OVERLAP_SUFFIX_SUFFIX)) return true;
			}
		}
		else {
			if (overhangSide2) {
				if (areOverlapping(neighbor1,neighbor2,Constants.OVERLAP_SUFFIX_SUFFIX) || areOverlapping(neighbor2,neighbor1,Constants.OVERLAP_PREFIX_PREFIX)) return true;
			}
			else {
				if (areOverlapping(neighbor1,neighbor2,Constants.OVERLAP_SUFFIX_PREFIX) || areOverlapping(neighbor2,neighbor1,Constants.OVERLAP_SUFFIX_PREFIX)) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Checks whether every neighbor $v$ on $side1$ of $node1$ (excluding possibly 
	 * $node1$ and $node2$) is also a neighbor of $node2$ on its $side2$, and if the side 
	 * of $v$ is the same in both cases.
	 *
	 * Remark: the procedure could list the neighbors that do not satisfy this property,
	 * so that the caller might check e.g. containment in $node2$ using the interval 
	 * graph. We omit this feature since it is not needed at the moment.
	 * 
	 * @param side* TRUE=prefix; FALSE=suffix.
	 */
	public static final boolean neighborsAreSubset(int node1, boolean side1, int node2, boolean side2) {
		int i, j;
		final int last = lastNeighbor[node1];
		BidirectedEdge edge;
		
		for (i=0; i<=last; i++) {
			edge=neighbors[node1][i];
			if (side1) {
				if (edge.node1==node1) {
					if (edge.type==Constants.OVERLAP_PREFIX_PREFIX && edge.node2!=node2 && edge.node2!=node1 && !findEdge(node2,side2,edge.node2,true)) return false;
					else if (edge.type==Constants.OVERLAP_PREFIX_SUFFIX && edge.node2!=node2 && edge.node2!=node1 && !findEdge(node2,side2,edge.node2,false)) return false;
				}
				else {
					if (edge.type==Constants.OVERLAP_PREFIX_PREFIX && edge.node1!=node2 && edge.node1!=node1 && !findEdge(node2,side2,edge.node1,true)) return false;
					else if (edge.type==Constants.OVERLAP_SUFFIX_PREFIX && edge.node1!=node2 && edge.node1!=node1 && !findEdge(node2,side2,edge.node1,false)) return false;
				}
			}
			else {
				if (edge.node1==node1) {
					if (edge.type==Constants.OVERLAP_SUFFIX_PREFIX && edge.node2!=node2 && edge.node2!=node1 && !findEdge(node2,side2,edge.node2,true)) return false;
					else if (edge.type==Constants.OVERLAP_SUFFIX_SUFFIX && edge.node2!=node2 && edge.node2!=node1 && !findEdge(node2,side2,edge.node2,false)) return false;
				}
				else {
					if (edge.type==Constants.OVERLAP_PREFIX_SUFFIX && edge.node1!=node2 && edge.node1!=node1 && !findEdge(node2,side2,edge.node1,true)) return false;
					else if (edge.type==Constants.OVERLAP_SUFFIX_SUFFIX && edge.node1!=node2 && edge.node1!=node1 && !findEdge(node2,side2,edge.node1,false)) return false;
				}
			}
		}
		return true;
	}
	
	
	/**
	 * @param outNodes output array: stores the (not necessarily distinct) interval graph
	 * nodes that correspond to the neighbors of $node$, in no specific order; the array 
	 * is assumed to be of sufficient size;
	 * @return the number of edges incident to $node$ that overlap with its prefix (if
	 * $side=TRUE$) or with its suffix (if $side=FALSE$).
	 */
	public static final int edgesOnSide(int node, boolean side, IntervalGraph.Node[] outNodes) {
		final int last = lastNeighbor[node];
		int i, out;
		BidirectedEdge edge;

		out=0;
		for (i=0; i<=last; i++) {
			edge=neighbors[node][i];
			if ( edge.node1==node && 
				 ( side && (edge.type==Constants.OVERLAP_PREFIX_PREFIX || edge.type==Constants.OVERLAP_PREFIX_SUFFIX) ) ||
				 ( !side && (edge.type==Constants.OVERLAP_SUFFIX_PREFIX || edge.type==Constants.OVERLAP_SUFFIX_SUFFIX) )
			   ) {
			    outNodes[out]=intervalGraphPointers[edge.node2];
				out++;
			}
			else if ( edge.node2==node && 
				 	  ( side && (edge.type==Constants.OVERLAP_PREFIX_PREFIX || edge.type==Constants.OVERLAP_SUFFIX_PREFIX) ) ||
				      ( !side && (edge.type==Constants.OVERLAP_PREFIX_SUFFIX || edge.type==Constants.OVERLAP_SUFFIX_SUFFIX) )
			   		) {
			    outNodes[out]=intervalGraphPointers[edge.node1];
				out++;
			}
		}
		return out;
	}
	
	
	/**
	 * Breaks (at a shortest overlap) a cycle in the bidirected graph, such that no node 
	 * in the cycle has neighbors outside the cycle.
	 *
	 * Remark: such cycles can exist in the graph (and they can even have just two nodes),
	 * even after the cycle-removal procedures of $IntervalGraph.buildIntervalGraph()$.
	 *
	 * Remark: the procedure uses global array $components[]$ as temporary space.
	 */
	public static final void breakSimpleCycles() {
		boolean found;
		int i, j, k;
		int neighbor, node, nextNode, minNode, minEdge, overlap, minOverlap;
		BidirectedEdge edge;
		
		Math.set(components,nNodes-1,0);
		for (i=0; i<nNodes; i++) {
			if (components[i]!=0 || !breakSimpleCycles_useNode(i)) {
				components[i]=1;
				continue;
			}
			
			// Finding a simple cycle from $i$.
			components[i]=2;
			j=neighbors[i][0].getTo(i);
			if (breakSimpleCycles_useNode(j)) {
				node=j;
				overlap=neighbors[i][0].getOverlapLength();
				minNode=i; minEdge=0; minOverlap=overlap;
			}
			else {
				node=-1;
				minNode=-1; minEdge=-1; minOverlap=Math.POSITIVE_INFINITY;
			}
			while (node>=0) {
				components[node]=2;
				nextNode=-1;
				for (j=0; j<=lastNeighbor[node]; j++) {
					neighbor=neighbors[node][j].getTo(node);
					if (neighbor==i) {
						overlap=neighbors[node][j].getOverlapLength();
						if (overlap<minOverlap) {
							minNode=node; minEdge=j; minOverlap=overlap;
						}
					}
					else if (components[neighbor]==0 && breakSimpleCycles_useNode(neighbor)) {
						nextNode=neighbor;
						overlap=neighbors[node][j].getOverlapLength();
						if (overlap<minOverlap) {
							minNode=node; minEdge=j; minOverlap=overlap;
						}
						break;
					}
				}
				node=nextNode;
			}
			found=components[neighbors[i][1].getTo(i)]==2;
			
			// Cleaning marks in $components$.
			node=i;
			while (node>=0) {
				components[node]=1;
				nextNode=-1;
				for (j=0; j<=lastNeighbor[node]; j++) {
					neighbor=neighbors[node][j].getTo(node);
					if (components[neighbor]==2) {
						components[neighbor]=1;
						nextNode=neighbor;
						break;
					}
				}
				node=nextNode;
			}
			
			// Disconnecting the cycle
			if (found) {
				edge=neighbors[minNode][minEdge];
				for (j=minEdge; j<lastNeighbor[minNode]; j++) neighbors[minNode][j]=neighbors[minNode][j+1];
				neighbors[minNode][lastNeighbor[minNode]]=edge;
				lastNeighbor[minNode]--;
				neighbor=edge.getTo(minNode);
				for (j=0; j<=lastNeighbor[neighbor]; j++) {
					if (neighbors[neighbor][j]==edge) {
						for (k=j; k<lastNeighbor[neighbor]; k++) neighbors[neighbor][k]=neighbors[neighbor][k+1];
						neighbors[neighbor][lastNeighbor[neighbor]]=edge;
						lastNeighbor[neighbor]--;
						break;
					}
				}
			}
		}
	}
	
	
	/**
	 * @return TRUE iff node $i$ has exactly one edge that uses its prefix, and exactly 
	 * one edge that uses its suffix.
	 */
	private static final boolean breakSimpleCycles_useNode(int i) {
		boolean foundPrefix, foundSuffix;
		int j;
		
		if (lastNeighbor[i]!=1) return false;
		foundPrefix=false; foundSuffix=false;
		for (j=0; j<=lastNeighbor[i]; j++) {
			if ( (i==neighbors[i][j].node1 && (neighbors[i][j].type==Constants.OVERLAP_PREFIX_PREFIX || neighbors[i][j].type==Constants.OVERLAP_PREFIX_SUFFIX)) ||
				 (i==neighbors[i][j].node2 && (neighbors[i][j].type==Constants.OVERLAP_PREFIX_PREFIX || neighbors[i][j].type==Constants.OVERLAP_SUFFIX_PREFIX))
			   ) {
			    if (foundPrefix) return false;
				foundPrefix=true;
			}
			if ( (i==neighbors[i][j].node1 && (neighbors[i][j].type==Constants.OVERLAP_SUFFIX_PREFIX || neighbors[i][j].type==Constants.OVERLAP_SUFFIX_SUFFIX)) ||
				 (i==neighbors[i][j].node2 && (neighbors[i][j].type==Constants.OVERLAP_PREFIX_SUFFIX || neighbors[i][j].type==Constants.OVERLAP_SUFFIX_SUFFIX))
			   ) {
				if (foundSuffix) return false;
				foundSuffix=true;
			}
		}
		return foundPrefix && foundSuffix;
	}
	
	
	/**
	 * @return TRUE iff the graph contains at least one fork, i.e. a node with two edges
	 * on the same side.
	 */
	public static final boolean containsFork() {
		if (nNodes<=2) return false;
		for (int i=0; i<nNodes; i++) {
			if (isFork(i)!=0) return true;
		}
		return false;
	}
	
	
	/**
	 * @return 0=node $i$ is not a fork on any side; 1=fork only on prefix side; 
	 * 2=fork only on suffix side; 3=fork on both sides.
	 */
	private static final int isFork(int i) {
		int j;
		int nPrefix, nSuffix;
		final int last = lastNeighbor[i];
		BidirectedEdge edge;
		
		if (last<=0) return 0;
		nPrefix=0; nSuffix=0;
		for (j=0; j<=last; j++) {
			edge=neighbors[i][j];
			if (i==edge.node1) {
				if (edge.type==Constants.OVERLAP_PREFIX_PREFIX || edge.type==Constants.OVERLAP_PREFIX_SUFFIX) nPrefix++;
				else nSuffix++;
			}
			else {
				if (edge.type==Constants.OVERLAP_PREFIX_PREFIX || edge.type==Constants.OVERLAP_SUFFIX_PREFIX) nPrefix++;
				else nSuffix++;
			}
		}
		if (nPrefix>1) return nSuffix>1?3:1;
		else if (nSuffix>1) return 2;
		else return 0;
	}
	
	
	/**
	 * Sets $components[i]$ to the ID of the connected component of the graph to which
	 * node $i$ belongs (-1 if node $i$ was removed).
	 *
	 * @return number of connected components.
	 */
	public static final int getConnectedComponents(int[] components) {
		int i, j;
		int from, to, top, lastComponent;
		
		Math.set(components,nNodes-1,-1);
		lastComponent=-1;
		for (i=0; i<nNodes; i++) {
			if (removed[i] || components[i]!=-1) continue;
			lastComponent++; 
			components[i]=lastComponent; top=0; stack[0]=i;
			while (top>=0) {
				from=stack[top--];
				for (j=0; j<=lastNeighbor[from]; j++) {
					to=neighbors[from][j].getTo(from);
					if (components[to]!=-1) continue;
					components[to]=lastComponent;
					stack[++top]=to;
				}
			}
		}
		return lastComponent+1;
	}
	
	
	/**
	 * Basic consistency checks
	 */
	public static final void checkConsistency() {
		int i, j;
		int last;
		
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			last=lastNeighbor[i];
			for (j=0; j<=last; j++) {
				if (neighbors[i][j].node1!=i && neighbors[i][j].node2!=i) {
					System.err.println("BidirectedGraph.checkConsistency> ERROR: the following edge does not contain adjacent node "+i+": "+neighbors[i][j]);
					System.exit(1);
				}
				if (neighbors[i][j].node1<0 || neighbors[i][j].node2<0) {
					System.err.println("BidirectedGraph.checkConsistency> ERROR: the following edge contains a negative adjacent node: "+neighbors[i][j]);
					System.exit(1);
				}
			}
		}
	}




	
	
	
	
	// ---------------------------- CONTRACTION PROCEDURES -------------------------------
	
	/**
	 * Creates the bidirected graph induced by the current graph, and by a surjective 
	 * mapping of the non-removed existing nodes, to a smaller set of new nodes. The edges 
	 * of each node in the new graph are sorted by the ID of the destination node.
	 * Removed nodes are discarded.
	 *
	 * @param newNodes read, start, end;
	 * @param substringOrIdentical TRUE: every node in the current graph is a substring of
	 * a node in the new graph; FALSE: every node in the current graph is approximately
	 * identical to a node in the new graph;
	 * @param node2newNode 0: ID of the new node, or -1 if the node does not get 
	 * projected; 1 (2): if $substringOrIdentical=true$, relative start (end) inside the 
	 * new node; if $substringOrIdentical=false$, column one is the orientation of the 
	 * node with respect to its new node.
	 */
	public static final void contract(int[][] newNodes, int lastNewNode, IntervalGraph.Node[] newPointers, int[][] node2newNode, boolean substringOrIdentical) {
		int i, j;
		BidirectedEdge newEdge;
		int[] tmpArray;
		BidirectedEdge[][] tmpMatrix;
		
		cleanMarkedFlags(true);
		Math.set(newLastNeighbor,lastNewNode,-1);
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			for (j=0; j<=lastNeighbor[i]; j++) {
				if (neighbors[i][j].marked) continue;
				newEdge=substringOrIdentical?getContractedEdge_substring(neighbors[i][j],newNodes,node2newNode):getContractedEdge_identical(neighbors[i][j],newNodes,node2newNode);
				if (newEdge.node1!=newEdge.node2) {  // No self-loops
					addEdge(newEdge,newNeighbors,newLastNeighbor,nodeLength);
				}
				neighbors[i][j].marked=true;
			}
		}
		Math.set(removed,lastNewNode,false); nRemoved=0;
		tmpArray=lastNeighbor; lastNeighbor=newLastNeighbor; newLastNeighbor=tmpArray;
		tmpMatrix=neighbors; neighbors=newNeighbors; newNeighbors=tmpMatrix;
		nNodes=lastNewNode+1;
		sortNeighbors();
		for (i=0; i<=lastNewNode; i++) nodeLength[i]=newNodes[i][2]-newNodes[i][1]+1;
		for (i=0; i<=lastNewNode; i++) intervalGraphPointers[i]=newPointers[i];
	}
	
	
	/**
	 * Returns the edge in the contracted bidirected graph induced by edge $oldEdge$ in 
	 * the original graph, and by a surjective mapping of existing nodes to a smaller set 
	 * of new nodes, such that every existing node is a substring of a new node.
	 */
	private static final BidirectedEdge getContractedEdge_substring(BidirectedEdge oldEdge, int[][] newNodes, int[][] node2newNode) {
		final int node1 = node2newNode[oldEdge.node1][0];
		final int node1Length = newNodes[node1][2]-newNodes[node1][1]+1;
		final int node1Start = node2newNode[oldEdge.node1][1];
		final int node1End = node2newNode[oldEdge.node1][2];
		final int node2 = node2newNode[oldEdge.node2][0];
		final int node2Length = newNodes[node2][2]-newNodes[node2][1]+1;
		final int node2Start = node2newNode[oldEdge.node2][1];
		final int node2End = node2newNode[oldEdge.node2][2];
		int overhang1=-1, overhang2=-1;
		BidirectedEdge newEdge = getEdge();
		
		if (oldEdge.type==Constants.OVERLAP_PREFIX_PREFIX) {
			overhang1=oldEdge.overhang1+node1Length-node1End;
			overhang2=oldEdge.overhang2+node2Length-node2End;
		}
		else if (oldEdge.type==Constants.OVERLAP_PREFIX_SUFFIX) {
			overhang1=oldEdge.overhang1+node1Length-node1End;
			overhang2=oldEdge.overhang2+node2Start;
		}
		else if (oldEdge.type==Constants.OVERLAP_SUFFIX_PREFIX) {
			overhang1=oldEdge.overhang1+node1Start;
			overhang2=oldEdge.overhang2+node2Length-node2End;
		}
		else if (oldEdge.type==Constants.OVERLAP_SUFFIX_SUFFIX) {
			overhang1=oldEdge.overhang1+node1Start;
			overhang2=oldEdge.overhang2+node2Start;
		}
		newEdge.node1=node1; newEdge.node2=node2; newEdge.type=oldEdge.type;
		newEdge.overhang1=overhang1; newEdge.overhang2=overhang2;
		newEdge.nAlignments=oldEdge.nAlignments;
		return newEdge;
	}
	
	
	/**
	 * Returns the edge in the contracted bidirected graph induced by edge $oldEdge$ in 
	 * the original graph, and by a surjective mapping of existing nodes to a smaller set 
	 * of new nodes, such that every existing node is approximately identical to a new 
	 * node.
	 *
	 * Remark: the procedure uses only columns 0 and 1 of $node2newNode$, where column 1
	 * is the orientation with respect to the representative in column 0.
	 */
	private static final BidirectedEdge getContractedEdge_identical(BidirectedEdge oldEdge, int[][] newNodes, int[][] node2newNode) {
		final int node1 = node2newNode[oldEdge.node1][0];
		final int node2 = node2newNode[oldEdge.node2][0];
		final int oldLength1 = nodeLength[oldEdge.node1];
		final int oldLength2 = nodeLength[oldEdge.node2];
		final int newLength1 = newNodes[node1][2]-newNodes[node1][1]+1;
		final int newLength2 = newNodes[node2][2]-newNodes[node2][1]+1;
		int rc = -1;
		if (node2newNode[oldEdge.node1][1]==1 && (node2newNode[oldEdge.node2][1]==0 || node2newNode[oldEdge.node2][1]==2)) rc=0;
		else if ((node2newNode[oldEdge.node1][1]==0 || node2newNode[oldEdge.node1][1]==2) && node2newNode[oldEdge.node2][1]==1) rc=1;
		
		BidirectedEdge newEdge = getEdge();
		newEdge.node1=node1; newEdge.node2=node2;
		newEdge.type=rc==-1?oldEdge.type:Constants.reverseComplementOverlap(oldEdge.type,rc);
		newEdge.overhang1=Math.round(oldEdge.overhang1*(((double)newLength1)/oldLength1));
		newEdge.overhang2=Math.round(oldEdge.overhang2*(((double)newLength2)/oldLength2));
		newEdge.nAlignments=oldEdge.nAlignments;
		return newEdge;
	}
	
	
	/**
	 * Compacts the graph, deleting all nodes that were marked as removed.
	 */
	public static final void deleteRemoved() {
		int i, j, k;
		BidirectedEdge tmpEdge;
		BidirectedEdge[] tmpNeighbors;
		if (nRemoved==0) return;
		
		j=-1;
		for (i=0; i<nNodes; i++) {
			if (!removed[i]) stack[i]=++j;  // Using $stack$ as temporary space.
		}
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			k=-1;
			for (j=0; j<=lastNeighbor[i]; j++) {
				if (removed[neighbors[i][j].getTo(i)]) continue;
				neighbors[i][j].node1=stack[neighbors[i][j].node1];
				neighbors[i][j].node2=stack[neighbors[i][j].node2];
				k++;
				if (k!=j) {
					tmpEdge=neighbors[i][k];
					neighbors[i][k]=neighbors[i][j];
					neighbors[i][j]=tmpEdge;
				}
			}
			lastNeighbor[i]=k;
		}
		j=-1;
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			j++;
			lastNeighbor[j]=lastNeighbor[i];
			tmpNeighbors=neighbors[j];
			neighbors[j]=neighbors[i];
			neighbors[i]=tmpNeighbors;
			nodeLength[j]=nodeLength[i];
			intervalGraphPointers[j]=intervalGraphPointers[i];
		}
		nNodes-=nRemoved; nRemoved=0;
		Math.set(removed,nNodes-1,false);
	}
	
	
	
	
	
	
	
	
	// ------------------------------ NUMBER OF ASSEMBLIES -------------------------------
	
	/**
	 * Computes the number of valid walks, in the transitively-reduced graph, from a node 
	 * whose incident edges use all just one end, to a different node whose incident edges
	 * use all just one end (or to the other end of the same node). This is an estimate of
	 * the number of distinct assemblies. See $getReachableSubgraph()$ for a definition of
	 * valid walks. The procedure also estimates the min/max string length of a valid 
	 * walk.
	 *
	 * Remark: the procedure does not work on bidirected graphs that contain a directed 
	 * cycle. This might be too conservative, since cycles might not affect the paths 
	 * between two one-end nodes (or between any pair of one-end nodes).
	 *
	 * Remark: the procedure handles the case in which the bidirected graph has multiple
	 * connected components, possibly of size one.
	 *
	 * Remark: this procedure takes $O(T+mE+m^2)$ time, where $T$ is the time for 
	 * transitive reduction, $E$ is the number of edges in the transitively-reduced graph,
	 * and $m$ is the number of nodes with incident edges on just one end. 
	 * If the graph is a kernel this might be too much, since $m$ could be proportional to
	 * the number of fragments of a single repeat created by insertion, and to the number 
	 * of distinct repeats that share substrings.
	 *
	 * Remark: the main loop of this procedure could be parallelized.
	 *
	 * @param nPaths temporary space, with at least $nNodes$ cells;
	 * @param pathLengths temporary space, with at least $nNodes$ rows and 4 columns;
	 * @param nComp number of connected components of the bidirected graph;
	 * @param comp for every node $i$ of the bidirected graph, the ID of the connected 
	 * component of the bidirected graph it belongs to;
	 * @param compPaths output array: for every connected component, the sequence of nodes
	 * in a longest assembly path;
	 * @param compLengths output array: number of nodes in each path of $compPaths$;
	 * @param compStringLengths output array: string length of each path in $compPaths$;
	 * @param out output array; $out[0]=-1$ if the graph contains a directed cycle, the 
	 * number of assemblies otherwise; $out[1]$ and $out[2]$: min/max string length of an 
	 * assembly; $out[3]$: number of transitively-reduced edges; $out[4]$: number of 
	 * one-end nodes; out[5]: number of nodes without neighbors.
	 */
	public static final void estimateAssemblies(int[] nPaths, int[][] pathLengths, int nComp, int[] comp, int[][] compPaths, int[] compLengths, int[] compStringLengths, int[] out) {
		boolean oneEndReached;
		int i, j, w;
		int source, sourceEnd, sourceComp, nVertices, isDAG, sourceLength, length;
		int minLength, maxLength, maxLengthPrime, maxEnd,vertex, vertexOrientation, totalNPaths, nComponents;
		final int nOneEnd, nSingletons;
		
		if (nNodes==1) {
			out[0]=1;
			out[1]=nodeLength[0]; out[2]=nodeLength[0];
			out[3]=0;
			out[4]=0;
			out[5]=1;
			return;
		}
		getNodesWithOneEnd(oneEnd,oneEndFlags,tmp);
		nOneEnd=tmp[0]; nSingletons=tmp[1]; minLength=tmp[2]; maxLength=tmp[3];
		if (nOneEnd==0) {
			if (nSingletons==nNodes) {
				out[0]=nSingletons;
				out[1]=minLength; out[2]=maxLength;
				out[3]=0;
				out[4]=0;
				out[5]=nSingletons;
			}
			else {
				out[0]=-1;
				out[1]=-1; out[2]=-1;
				out[3]=-1;
				out[4]=0;
				out[5]=nSingletons;
			}
			return;
		}
		out[3]=transitiveReduction(shortestDistance,isActive,queue);
		out[4]=nOneEnd;
		out[5]=nSingletons;
		
		// Initializing reused data structures
		Math.set(reachedFrom,nNodes-1,-1);
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			for (j=0; j<=lastNeighbor[i]; j++) {
				neighbors[i][j].traversed12=false;
				neighbors[i][j].traversed21=false;
			}
		}
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			original2dag[i][0]=-1;
			original2dag[i][1]=-1;
		}
		Math.set(lastInNeighbor,(nNodes<<1)-1,-1);
		Math.set(lastOutNeighbor,(nNodes<<1)-1,-1);
		Math.set(nPaths,nNodes-1,0);
		for (i=0; i<nNodes; i++) { 
			if (removed[i]) continue;
			pathLengths[i][0]=0; pathLengths[i][1]=0; pathLengths[i][2]=-1; pathLengths[i][3]=-1;
		}
		Math.set(compLengths,nComp-1,0); Math.set(compStringLengths,nComp-1,0);
		
		// Counting paths
		totalNPaths=0;		
		for (i=0; i<nOneEnd; i++) {
			source=oneEnd[i]>=0?oneEnd[i]:-1-oneEnd[i];
			sourceEnd=oneEnd[i]>=0?1:0;
			oneEndReached=getReachableSubgraph(source,sourceEnd,true/*Sorting required by $getReachableDAG()$*/,reachedFrom,lastInNeighborPrime/*Used as temporary space*/,oneEndFlags);
			if (!oneEndReached) {
				getReachableSubgraph_clean(source,true,reachedFrom);
				out[0]=-1;
				out[1]=-1; out[2]=-1;
				return;
			}
			nVertices=getReachableDAG(source,sourceEnd,original2dag,dag2original,inOverhangs);
			for (j=0; j<nVertices; j++) lastInNeighbor[j]++;
			for (j=0; j<nVertices; j++) lastOutNeighbor[j]++;
			nComponents=DAG.getConnectedComponents(nVertices,inNeighbors,lastInNeighbor,outNeighbors,lastOutNeighbor,components,stack);
			if (nComponents<=0 || nComponents>nVertices) {
				IO.printCriticalErr("Error in estimateAssemblies(): wrong number of connected components.");
				System.exit(1);
			}		
			System.arraycopy(lastInNeighbor,0,lastInNeighborPrime,0,nVertices);  // Since the following procedure modifies $lastInNeighbor$.			
			if (nComponents>nVerticesPerComponent.length) nVerticesPerComponent = new int[nComponents<<1];
			isDAG=DAG.topologicalSort(nVertices,inNeighbors,lastInNeighborPrime,outNeighbors,lastOutNeighbor,components,1,minimalVertices,-1,nVerticesPerComponent,sorted2original,original2sorted);
			if (isDAG!=0) {
				getReachableDAG_clean(source,sourceEnd,original2dag,dag2original,nVertices);
				getReachableSubgraph_clean(source,true,reachedFrom);
				out[0]=-1; 
				out[1]=-1; out[2]=-1;
				return;
			}
			DAG.nPaths(original2dag[source][sourceEnd],inNeighbors,lastInNeighbor,sorted2original,original2sorted,nVertices,nPaths);
			for (j=0; j<nOneEnd; j++) {  // The loop is correct also for $j=i$.
				vertex=oneEnd[j]>=0?oneEnd[j]:-1-oneEnd[j];
				vertexOrientation=oneEnd[j]>=0?0:1;
				if (original2dag[vertex][vertexOrientation]>=0) totalNPaths+=nPaths[original2dag[vertex][vertexOrientation]];
			}
			DAG.nPaths_clean(original2dag[source][sourceEnd],sorted2original,original2sorted,nVertices,nPaths);
			DAG.pathWeights(original2dag[source][sourceEnd],inNeighbors,lastInNeighbor,outNeighbors,lastOutNeighbor,sorted2original,original2sorted,nVertices,inOverhangs,pathLengths);
			sourceLength=intervalGraphPointers[source].length();
			maxLengthPrime=0; maxEnd=-1;
			for (j=0; j<nOneEnd; j++) {  // The loop is correct also for $j=i$.
				vertex=oneEnd[j]>=0?oneEnd[j]:-1-oneEnd[j];
				vertexOrientation=oneEnd[j]>=0?0:1;
				if (original2dag[vertex][vertexOrientation]==-1) continue;
				w=pathLengths[original2dag[vertex][vertexOrientation]][0];
				if (IO.CONSISTENCY_CHECKS && w==0) {
					System.err.println("ERROR: DAG vertex "+original2dag[vertex][vertexOrientation]+" is reachable from DAG vertex "+original2dag[source][sourceEnd]+" but has min path length zero?!");
					System.exit(1);
				}
				length=sourceLength+w;
				if (length<minLength) minLength=length;
				length=sourceLength+pathLengths[original2dag[vertex][vertexOrientation]][1];
				if (length>maxLengthPrime) {
					maxLengthPrime=length;
					maxEnd=j;
				}
			}
			if (maxLengthPrime>maxLength) maxLength=maxLengthPrime;
			sourceComp=comp[source];
			if (maxLengthPrime>compStringLengths[sourceComp]) {
				compStringLengths[sourceComp]=maxLengthPrime;
				vertex=oneEnd[maxEnd]>=0?oneEnd[maxEnd]:-1-oneEnd[maxEnd];
				vertexOrientation=oneEnd[maxEnd]>=0?0:1;
				DAG.getLongestPath(original2dag[vertex][vertexOrientation],pathLengths,compPaths,compLengths,sourceComp,true);
				for (j=0; j<compLengths[sourceComp]; j++) compPaths[sourceComp][j]=dag2original[compPaths[sourceComp][j]];
			}
			DAG.pathWeights_clean(original2dag[source][sourceEnd],sorted2original,original2sorted,nVertices,pathLengths);
			getReachableDAG_clean(source,sourceEnd,original2dag,dag2original,nVertices);
			getReachableSubgraph_clean(source,true,reachedFrom);
		}
		if (nSingletons>0) {
			for (i=0; i<nNodes; i++) {
				if (removed[i] || lastNeighbor[i]>=0) continue;
				sourceComp=comp[i];
				compLengths[sourceComp]=1;
				compStringLengths[sourceComp]=nodeLength[i];
				if (compPaths[sourceComp].length<1) compPaths[sourceComp] = new int[1];
				compPaths[sourceComp][0]=i;
				totalNPaths++;
			}
		}
		if (totalNPaths==0) {
			out[0]=-1;
			out[1]=-1; out[2]=-1;
			return;
		}
		if (IO.CONSISTENCY_CHECKS) {
			if (minLength==Math.POSITIVE_INFINITY || (totalNPaths-nSingletons)%2!=0) {
				System.err.println("estimateAssemblies> ERROR with the following bidirected graph:");
				for (i=0; i<nNodes; i++) {
					if (removed[i]) continue;
					System.err.print(i+": ");
					for (j=0; j<=lastNeighbor[i]; j++) System.err.print(neighbors[i][j]+" ");
					System.err.println();
				}
				System.err.println("totalNPaths="+totalNPaths);
				System.err.println("oneEnd:");
				for (i=0; i<nOneEnd; i++) System.err.print(oneEnd[i]+",");
				System.err.println();
				System.exit(1);
			}
			checkLongestPaths(nComp,compPaths,compLengths);
		}
		out[0]=totalNPaths>>1; out[1]=minLength; out[2]=maxLength;
	}
	
	
	/**
	 * Checks whether longest paths of different connected components share nodes.
	 */
	private static final void checkLongestPaths(int nComp, int[][] compPaths, int[] compLengths) {
		int i, j;
		int max, lastI, lastJ;
		
		max=0;
		for (i=0; i<nComp; i++) {
			if (compLengths[i]>max) max=compLengths[i];
		}
		max*=3;
		if (max>stack.length) stack = new int[max];
		for (i=0; i<nComp; i++) {
			System.arraycopy(compPaths[i],0,stack,0,compLengths[i]);
			lastI=Math.makePositive(stack,0,compLengths[i]-1);
			for (j=i+1; j<nComp; j++) {
				System.arraycopy(compPaths[j],0,stack,max,compLengths[j]);
				lastJ=Math.makePositive(stack,max,max+compLengths[j]-1);
				if (Math.nonemptyIntersection(stack,0,lastI,stack,max,lastJ)) {
					System.err.println("estimateAssemblies> ERROR: the longest paths of connected components "+i+" and "+j+" share nodes?!");
					System.exit(1);
				}
			}
		}
	}
	
	
	/**
	 * Stores in $oneEnd$ the (not necessarily sorted) list of all nodes such that, for 
	 * every node, all the edges it is incident to involve just one end of the node. 
	 * $oneEnd[i]$ equals the node ID if its only end is a prefix; it equals minus one 
	 * minus the node ID if its only end is a suffix.
	 *
	 * Remark: $transitiveReduction()$ does not change if a node is one-end or not.
	 *
	 * @param out output array, of size at least equal to the number of nodes in the 
	 * bidirected graph; the procedure sets to TRUE all and only the nodes with one end;
	 * @param stats output array: 0=number of elements in $oneEnd$; 1=number of nodes with
	 * no neighbor; 2(3)=min(max) string length of a node with no neighbor.
	 */
	public static final void getNodesWithOneEnd(int[] oneEnd, boolean[] out, int[] stats) {
		boolean prefixEnd, suffixEnd;
		int i, j;
		int lastOneEnd, type, nSingletons, length, minLength, maxLength;
		
		Math.set(out,nNodes-1,false);
		lastOneEnd=-1; nSingletons=0; minLength=Math.POSITIVE_INFINITY; maxLength=0;
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			if (lastNeighbor[i]==-1) {
				nSingletons++;
				length=nodeLength[i];
				if (length<minLength) minLength=length;
				if (length>maxLength) maxLength=length;
				continue;
			}
			prefixEnd=false; suffixEnd=false;
			for (j=0; j<=lastNeighbor[i]; j++) {
				type=neighbors[i][j].type;
				if (neighbors[i][j].node1==i) {
					if (type==Constants.OVERLAP_PREFIX_PREFIX || type==Constants.OVERLAP_PREFIX_SUFFIX) prefixEnd=true;
					else suffixEnd=true;
				}
				else {
					if (type==Constants.OVERLAP_SUFFIX_PREFIX || type==Constants.OVERLAP_PREFIX_PREFIX) prefixEnd=true;
					else suffixEnd=true;
				}
			}
			if (prefixEnd) {
				if (suffixEnd) continue;
				oneEnd[++lastOneEnd]=i;
				out[i]=true;
			}
			else {
				if (!suffixEnd) continue;
				oneEnd[++lastOneEnd]=-1-i;
				out[i]=true;
			}
		}
		stats[0]=lastOneEnd+1; stats[1]=nSingletons; stats[2]=minLength; stats[3]=maxLength;
	}
	
	
	/**
	 * @return -2: no end; -1: two ends; 0: one end, prefix used in alignments; 
	 * 1: one end, suffix used in alignments.
	 */
	public static final int isOneEnd(int i) {
		boolean prefixEnd, suffixEnd;
		int j, type;
		
		if (lastNeighbor[i]==-1) return -1;
		prefixEnd=false; suffixEnd=false;
		for (j=0; j<=lastNeighbor[i]; j++) {
			type=neighbors[i][j].type;
			if (neighbors[i][j].node1==i) {
				if (type==Constants.OVERLAP_PREFIX_PREFIX || type==Constants.OVERLAP_PREFIX_SUFFIX) prefixEnd=true;
				else suffixEnd=true;
			}
			else {
				if (type==Constants.OVERLAP_SUFFIX_PREFIX || type==Constants.OVERLAP_PREFIX_PREFIX) prefixEnd=true;
				else suffixEnd=true;
			}
		}
		if (prefixEnd) {
			if (suffixEnd) return -1;
			else return 0;
		}
		else if (suffixEnd) return 1;
		else return -2;
	}
	
	
	/**
	 * A walk in a bidirected graph is valid if it exits every node from the opposite 
	 * side from which it enters it, and if it uses every edge at most twice in two
	 * opposite orientations. In other words, we allow, in the reconstructed string, both 
	 * $W$ and the reverse of $W$ to occur, where $W$ is the string that results from 
	 * overlapping node $i$ with node $j$ for an edge $(i,j)$. The procedure marks all 
	 * nodes and edges in the bidirected graph that are reachable from $source$ with a 
	 * valid walk.
	 *
	 * Remark: the procedure avoids transitively-reduced edges.
	 *
	 * Remark: the reachable subgraph can contain cycles, and this could happen even if 
	 * every edge could be traversed just once.
	 *
	 * Remark: since an edge or a node can be traversed in both directions from the same 
	 * source, via different walks or even in the same walk, the algorithm can push
	 * the same node twice on the DFS stack, if it is traversed in a new direction.
	 *
	 * Remark: this procedure takes $O(E)$ time, where $E$ is the number of edges in the
	 * graph.
	 *
	 * Remark: the procedure assumes the $inSubgraph*$ marks of every edge to be 
	 * initialized to false. 
	 *
	 * @param reachedFrom the procedure assumes $reachedFrom$ to be initialized to all -1;
	 * for each reachable node, the procedure stores in $reachedFrom$ the end from which 
	 * it can be reached from $source$: -1: unreachable, 0=prefix, 1=suffix, 2=both;
	 * @param sourceEnd the side from which $source$ was reached before the procedure 
	 * starts: 0=prefix, 1=suffix;
	 * @param sort if TRUE, sorts $neighbors[i]$ by $-1-i$ for every reachable $i$;
	 * @param tmpArray (used only if $sort=TRUE$) temporary space, of size at least equal 
	 * to the total number of nodes in the bidirected graph;
	 * @param oneEndFlags input array: marks all and only the nodes that use just one end 
	 * in the bidirected graph;
	 * @return TRUE iff at least one node with one end, and different from $source$, is 
	 * reached during the traversal, or if $source$ has one end but can be reached from 
	 * both ends.
	 */	
	private static final boolean getReachableSubgraph(int source, int sourceEnd, boolean sort, int[] reachedFrom, int[] tmpArray, boolean[] oneEndFlags) {
		boolean pushed, out;
		int i, j;
		int top, from, to, type, last, previous;
		
		// Reachability
		out=false;
		reachedFrom[source]=sourceEnd;
		stack[0]=source; top=0; last=-1;
		while (top>=0) {
			from=stack[top--];
			if (from!=source && oneEndFlags[from]) out=true;
			if (sort) tmpArray[++last]=from;
			for (j=0; j<=lastNeighbor[from]; j++) {
				if (neighbors[from][j].isRedundant || neighbors[from][j].inSubgraph(from)) continue;
				to=neighbors[from][j].getTo(from);
				type=neighbors[from][j].getType(from);
				pushed=false;
				if (reachedFrom[from]==0 || reachedFrom[from]==2) {
					if (type==Constants.OVERLAP_SUFFIX_PREFIX) {
						neighbors[from][j].setSubgraph(from,true);
						if (reachedFrom[to]==-1) {
							reachedFrom[to]=0;
							pushed=true;
						}
						else if (reachedFrom[to]==1) {
							reachedFrom[to]=2;
							pushed=true;
						}
					}
					else if (type==Constants.OVERLAP_SUFFIX_SUFFIX) {
						neighbors[from][j].setSubgraph(from,true);
						if (reachedFrom[to]==-1) {
							reachedFrom[to]=1;
							pushed=true;
						}
						else if (reachedFrom[to]==0) {
							reachedFrom[to]=2;
							pushed=true;
						}
					}
				}
				if (reachedFrom[from]==1 || reachedFrom[from]==2) {
					if (type==Constants.OVERLAP_PREFIX_PREFIX) {
						neighbors[from][j].setSubgraph(from,true);
						if (reachedFrom[to]==-1) {
							reachedFrom[to]=0;
							pushed=true;
						}
						else if (reachedFrom[to]==1) {
							reachedFrom[to]=2;
							pushed=true;
						}
					}
					else if (type==Constants.OVERLAP_PREFIX_SUFFIX) {
						neighbors[from][j].setSubgraph(from,true);
						if (reachedFrom[to]==-1) {
							reachedFrom[to]=1;
							pushed=true;
						}
						else if (reachedFrom[to]==0) {
							reachedFrom[to]=2;
							pushed=true;
						}
					}
				}
				if (pushed) stack[++top]=to;
			}
		}
		if (!out && oneEndFlags[source] && reachedFrom[source]==2) out=true;
		
		// Sorting
		if (sort) {
			if (last>0) Arrays.sort(tmpArray,0,last+1);
			previous=tmpArray[0];
			if (lastNeighbor[previous]>0) {
				BidirectedEdge.order=-1-previous;
				Arrays.sort(neighbors[previous],0,lastNeighbor[previous]+1);
			}
			for (i=1; i<=last; i++) {
				if (tmpArray[i]==previous) continue;
				previous=tmpArray[i];
				if (lastNeighbor[previous]>0) {
					BidirectedEdge.order=-1-previous;
					Arrays.sort(neighbors[previous],0,lastNeighbor[previous]+1);
				}
			}
		}
		
		return out;
	}
	
	
	/**
	 * Resets to all -1 the data structures used by $getReachableSubgraph()$, in time 
	 * proportional to the size of the reachable subgraph.
	 *
	 * Remark: the procedure implicitly avoids transitively-reduced edges, since 
	 * $getReachableSubgraph()$ avoided them.
	 *
	 * @param sorted TRUE iff $neighbors[i]$ is sorted by $-1-i$.
	 */
	private static final void getReachableSubgraph_clean(int source, boolean sorted, int[] reachedFrom) {
		int i, top, from, to;
		
		reachedFrom[source]=-1;
		top=0;
		stack[top]=source;
		while (top>=0) {
			from=stack[top--];
			for (i=0; i<=lastNeighbor[from]; i++) {
				if (!neighbors[from][i].inSubgraph(from)) {
					if (sorted) break;
					else continue;
				}
				neighbors[from][i].setSubgraph(from,false);
				to=neighbors[from][i].getTo(from);
				stack[++top]=to;
				reachedFrom[to]=-1;
			}
		}
	}
	
	
	/**
	 * Stores in $inNeighbors$, $lastInNeighbor$, $outNeighbors$, $lastOutNeighbor$ the
	 * DAG that corresponds to the subgraph of the bidirected graph that is reachable from 
	 * $source$ via valid walks. The DAG is the union of all such valid (directed) walks.
	 * Every node in the bidirected graph is mapped to at most two vertices of the DAG, 
	 * that correspond to the two directions in which the node can be traversed by a walk.
	 *
	 * Remark: the DAG can contain directed cycles, e.g. when the same node of the 
	 * bidirected graph is used twice inside each of two valid walks, with different 
	 * orientations, but such orientations appear in different order in the two walks.
	 *
	 * Remark: the procedure takes linear time in the size of the reachable subgraph of 
	 * the bidirected graph from $source$.
	 *
	 * Remark: the procedure assumes that $getReachableSubgraph()$ has already been
	 * executed, with its sorting option activated (this implies that the procedure
	 * automatically discards transitively-reduced edges); that $original2dag$, 
	 * $lastInNeighbor$ and $lastOutNeighbor$ are initialized to all -1; and that the 
	 * $traversed*$ marks of every edge are initialized to false.
	 *
	 * @param sourceEnd the side from which $source$ was reached before the procedure 
	 * starts: 0=prefix, 1=suffix;
	 * @param original2dag output matrix: for every node $i$ in the bidirected graph, the 
	 * ID of the DAG vertex that corresponds to $i$ in the forward (column 0) or reverse 
	 * (column 1) orientation, if any (IDs in the DAG start from zero);
	 * @param dag2original output matrix: for every vertex in the DAG, the ID of the 
	 * corresponding node $i$ of the bidirected graph, in the forward ($i$) or reverse
	 * ($-1-i$) orientation;
	 * @param inOverhangs output matrix: for every incoming arc to a node, the overhang in
	 * the destination node;
	 * @return the number of vertices in the reachable DAG.
	 */
	private static final int getReachableDAG(int source, int sourceEnd, int[][] original2dag, int[] dag2original, int[][] inOverhangs) {
		int i;
		int top, from, to, dagFrom, dagTo, type, lastDagNode;
		int[] tmpArray;
			
		lastDagNode=0;  // $source$ in its only direction.
		original2dag[source][sourceEnd]=0;
		dag2original[0]=sourceEnd==0?source:-1-source;
		top=0;
		stack[top]=source;
		while (top>=0) {
			from=stack[top--];
			for (i=0; i<=lastNeighbor[from]; i++) {
				if (!neighbors[from][i].inSubgraph(from)) break;
				if (neighbors[from][i].isTraversed(from)) continue;
				neighbors[from][i].setTraversed(from,true);
				to=neighbors[from][i].getTo(from);
				type=neighbors[from][i].getType(from);
				if (type==Constants.OVERLAP_SUFFIX_PREFIX) {
					if (original2dag[from][0]==-1) {
						original2dag[from][0]=++lastDagNode;
						dag2original[lastDagNode]=from;
					}
					if (original2dag[to][0]==-1) {
						original2dag[to][0]=++lastDagNode;
						dag2original[lastDagNode]=to;
					}
					dagFrom=original2dag[from][0];
					dagTo=original2dag[to][0];
				}
				else if (type==Constants.OVERLAP_SUFFIX_SUFFIX) {
					if (original2dag[from][0]==-1) {
						original2dag[from][0]=++lastDagNode;
						dag2original[lastDagNode]=from;
					}
					if (original2dag[to][1]==-1) {
						original2dag[to][1]=++lastDagNode;
						dag2original[lastDagNode]=-1-to;
					}
					dagFrom=original2dag[from][0];
					dagTo=original2dag[to][1];
				}
				else if (type==Constants.OVERLAP_PREFIX_SUFFIX) {
					if (original2dag[from][1]==-1) {
						original2dag[from][1]=++lastDagNode;
						dag2original[lastDagNode]=-1-from;
					}
					if (original2dag[to][1]==-1) {
						original2dag[to][1]=++lastDagNode;
						dag2original[lastDagNode]=-1-to;
					}
					dagFrom=original2dag[from][1];
					dagTo=original2dag[to][1];
				}
				else {
					if (original2dag[from][1]==-1) {
						original2dag[from][1]=++lastDagNode;
						dag2original[lastDagNode]=-1-from;
					}
					if (original2dag[to][0]==-1) {
						original2dag[to][0]=++lastDagNode;
						dag2original[lastDagNode]=to;
					}
					dagFrom=original2dag[from][1];
					dagTo=original2dag[to][0];
				}
				lastOutNeighbor[dagFrom]++;
				if (lastOutNeighbor[dagFrom]==outNeighbors[dagFrom].length) {
					tmpArray = new int[outNeighbors[dagFrom].length+NEIGHBORS_UNIT];
					System.arraycopy(outNeighbors[dagFrom],0,tmpArray,0,outNeighbors[dagFrom].length);
					outNeighbors[dagFrom]=tmpArray;
				}
				outNeighbors[dagFrom][lastOutNeighbor[dagFrom]]=dagTo;
				lastInNeighbor[dagTo]++;
				if (lastInNeighbor[dagTo]==inNeighbors[dagTo].length) {
					tmpArray = new int[inNeighbors[dagTo].length+NEIGHBORS_UNIT];
					System.arraycopy(inNeighbors[dagTo],0,tmpArray,0,inNeighbors[dagTo].length);
					inNeighbors[dagTo]=tmpArray;
				}
				inNeighbors[dagTo][lastInNeighbor[dagTo]]=dagFrom;
				if (lastInNeighbor[dagTo]==inOverhangs[dagTo].length) {
					tmpArray = new int[inOverhangs[dagTo].length+NEIGHBORS_UNIT];
					System.arraycopy(inOverhangs[dagTo],0,tmpArray,0,inOverhangs[dagTo].length);
					inOverhangs[dagTo]=tmpArray;
				}
				inOverhangs[dagTo][lastInNeighbor[dagTo]]=from==neighbors[from][i].node1?neighbors[from][i].overhang2:neighbors[from][i].overhang1;
				stack[++top]=to;
			}
		}
		return lastDagNode+1;
	}
	
	
	/**
	 * Resets to all -1 the data structures used by $getReachableDAG()$, in time 
	 * proportional to the size of the reachable subgraph.
	 *
	 * Remark: the procedure assumes $neighbors[i]$ to have all the edges that belong to
	 * the reachable bidirected graph clustered at the beginning of the array.
	 *
	 * @param nVertices number of vertices in the DAG reachable from $source$.
	 */
	private static final void getReachableDAG_clean(int source, int sourceEnd, int[][] original2dag, int[] dag2original, int nVertices) {
		int i;
		int top, from, to, vertex;
		
		Math.set(lastInNeighbor,nVertices-1,-1);
		Math.set(lastOutNeighbor,nVertices-1,-1);
		vertex=original2dag[source][0];
		if (vertex>=0) dag2original[vertex]=-1;
		original2dag[source][0]=-1;
		vertex=original2dag[source][1];
		if (vertex>=0) dag2original[vertex]=-1;
		original2dag[source][1]=-1;
		top=0;
		stack[top]=source;
		while (top>=0) {
			from=stack[top--];
			for (i=0; i<=lastNeighbor[from]; i++) {
				if (!neighbors[from][i].inSubgraph(from)) break;
				if (!neighbors[from][i].isTraversed(from)) continue;
				neighbors[from][i].setTraversed(from,false);
				to=neighbors[from][i].getTo(from);
				stack[++top]=to;
				vertex=original2dag[to][0];
				if (vertex>=0) dag2original[vertex]=-1;
				original2dag[to][0]=-1;
				vertex=original2dag[to][1];
				if (vertex>=0) dag2original[vertex]=-1;
				original2dag[to][1]=-1;
			}
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// ------------------------------ TRANSITIVE REDUCTION -------------------------------
	
	/**
	 * Uses a simple $O(E^2)$ algorithm to mark as redundant edges implied by other edges, 
	 * where $E$ is the number of edges. An edge $(v,w)$ with a given overlap type is 
	 * redundant iff there is a valid walk, of length greater than one, that connects $v$ 
	 * and $w$ with the same ends as $(v,w)$ (see $reduce_impl()$ for details). Such walk 
	 * is detected by running reachability from every neighbor of $v$ different from $w$.
	 *
	 * Remark: usually an edge $u->w$ is called transitive if there is a node $v$ and 
	 * edges $u->v$ and $v->w$. Hopefully our approach marks more edges as redundant.
	 *
	 * Remark: the procedure does not require the neighbors of a node to be sorted.
	 *
	 * Remark: the procedure might be very slow for repeats, since $E$ might be 
	 * proportional to the number of occurrences of a repeat.
	 *
	 * Remark: this procedure could be parallelized.
	 *
	 * Remark: in \cite{baaijens2017novo} they define "double transitivity", in which 
	 * $u->w$ is double-transitive if $u->v$ and $v->w$ are transitive. Our procedure
	 * marks as redundant double-transitive edges as well.
	 *
	 * @param shortestDistance temporary space, with at least $nNodes$ rows and two 
	 * columns;
	 * @param isActive temporary space, with at least $nNodes$ rows and two columns;
	 * @param queue temporary space;
	 * @return the number of edges marked as redundant.
	 */
	public static final int transitiveReduction(int[][] shortestDistance, boolean[][] isActive, PriorityQueue<PriorityEnd> queue) {
		int i, j, v, w;
		int to, type, forkType, out;
		
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			for (j=0; j<=lastNeighbor[i]; j++) neighbors[i][j].isRedundant=false;
		}
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			isActive[i][0]=false;
			isActive[i][1]=false; 
		}
		Math.set(shortestDistance,nNodes-1,Math.POSITIVE_INFINITY);
		queue.clear();
		out=0;
		for (v=0; v<nNodes; v++) {
			if (removed[v] || lastNeighbor[v]<=0) continue;
			forkType=isFork(v);
			// First direction (starting from the prefix of $v$).
			if (forkType==2 || forkType==3) {
				reduce(v,0,shortestDistance,isActive,queue);
				for (i=0; i<=lastNeighbor[v]; i++) {
					to=neighbors[v][i].getTo(v);
					type=neighbors[v][i].getType(v);
					if ((type==Constants.OVERLAP_SUFFIX_PREFIX && !isActive[to][0]) || (type==Constants.OVERLAP_SUFFIX_SUFFIX && !isActive[to][1])) {
						neighbors[v][i].isRedundant=true;
						out++;
					}
				}
				for (i=0; i<=lastNeighbor[v]; i++) {
					to=neighbors[v][i].getTo(v);
					isActive[to][0]=false; isActive[to][1]=false;
				}
			}
			// Second direction (starting from the suffix of $v$).
			if (forkType==1 || forkType==3) {
				reduce(v,1,shortestDistance,isActive,queue);
				for (i=0; i<=lastNeighbor[v]; i++) {
					to=neighbors[v][i].getTo(v);
					type=neighbors[v][i].getType(v);
					if ((type==Constants.OVERLAP_PREFIX_PREFIX && !isActive[to][0]) || (type==Constants.OVERLAP_PREFIX_SUFFIX && !isActive[to][1])) {
						neighbors[v][i].isRedundant=true;
						out++;
					}
				}
				for (i=0; i<=lastNeighbor[v]; i++) {
					to=neighbors[v][i].getTo(v);
					isActive[to][0]=false; isActive[to][1]=false;
				}
			}
		}
		return out;
	}
	
	
	/**
	 * Tries to reduce all edges that are incident to $source$ on the end opposite to 
	 * $sourceEnd$ (0=prefix, 1=suffix).
	 *
	 * Remark: the procedure can reduce also self-loops.
	 *
	 * @param shortestDistance assumed to be initialized to all $Math.POSITIVE_INFINITY$;
	 * the procedure leaves it in this state when it ends;
	 * @param isActive[i][j] temporary space, with at least $nNodes$ rows and two columns;
	 * represents an edge to $i$ from $source$, connected to the prefix ($j=0$) or the 
	 * suffix ($j=1$) of $i$; the procedure sets some TRUE cells to FALSE;
	 * @param queue temporary space, assumed to be empty.	
	 */
	private static final void reduce(int source, int sourceEnd, int[][] shortestDistance, boolean[][] isActive, PriorityQueue<PriorityEnd> queue) {
		int i;
		int w, type, nActive, maxActive, delta, largestOffset, offset;
		
		// Marking active edges (including self-loops).
		nActive=0; largestOffset=0;
		for (i=0; i<=lastNeighbor[source]; i++) {
			w=neighbors[source][i].getTo(source);
			type=neighbors[source][i].getType(source);
			if (sourceEnd==0 && (type==Constants.OVERLAP_SUFFIX_PREFIX||type==Constants.OVERLAP_SUFFIX_SUFFIX)) {
				isActive[w][type==Constants.OVERLAP_SUFFIX_PREFIX?0:1]=true;
				nActive++;
				offset=neighbors[source][i].getOverhang(source);
				if (offset>largestOffset) largestOffset=offset;
			}
			else if (sourceEnd==1 && (type==Constants.OVERLAP_PREFIX_PREFIX||type==Constants.OVERLAP_PREFIX_SUFFIX)) {
				isActive[w][type==Constants.OVERLAP_PREFIX_PREFIX?0:1]=true;
				nActive++;
				offset=neighbors[source][i].getOverhang(source);
				if (offset>largestOffset) largestOffset=offset;
			}
		}
		if (nActive==0) return;
		
		// Reducing active edges
		for (i=0; i<=lastNeighbor[source]; i++) {
			if (neighbors[source][i].isRedundant) continue;
			w=neighbors[source][i].getTo(source);
			if (w==source) continue;  // Self-loops should not be used to start reachability
			type=neighbors[source][i].getType(source);
			delta=(isActive[w][0]?1:0)+(isActive[w][1]?1:0);
			if (delta==0) continue;  // It is useless to start reachability from a vertex that has already been reached from all possible ends.
			maxActive=nActive-delta;
			if (maxActive==0) continue;  // Nothing to reduce, except for edges to $w$.
			if (sourceEnd==0) {
				if (type==Constants.OVERLAP_SUFFIX_PREFIX && isActive[w][0]) nActive-=reduce_impl(source,w,0,neighbors[source][i].getOverhang(source),largestOffset,shortestDistance,isActive,queue);
				else if (type==Constants.OVERLAP_SUFFIX_SUFFIX && isActive[w][1]) nActive-=reduce_impl(source,w,1,neighbors[source][i].getOverhang(source),largestOffset,shortestDistance,isActive,queue);
			}
			else {
				if (type==Constants.OVERLAP_PREFIX_PREFIX && isActive[w][0]) nActive-=reduce_impl(source,w,0,neighbors[source][i].getOverhang(source),largestOffset,shortestDistance,isActive,queue);
				else if (type==Constants.OVERLAP_PREFIX_SUFFIX && isActive[w][1]) nActive-=reduce_impl(source,w,1,neighbors[source][i].getOverhang(source),largestOffset,shortestDistance,isActive,queue);
			}
			nActive+=delta;
		}
	}
	
	
	/**
	 * Marks as reduced an edge $e=(u,v)$ of the bidirected graph, iff there is a path 
	 * from a node $w$ such that: (1) the first edge $e'=(u,w)$ of the path uses the same 
	 * end of $u$ as $e$; (2) the last edge of the path uses the same end of $v$ as $e$;
	 * (3) the path traverses each node from opposite ends; (4) the offset between the 
	 * starting position of $u$ and the starting position of every interval in the path is
	 * at most the offset between $u$ and $v$.
	 *
	 * Remark: $u$ itself is not considered reachable, since this could mark as reduced 
	 * some edges from $u$.
	 *
	 * @param transitiveReductionSource ID of node $u$;
	 * @param source ID of node $w$;
	 * @param sourceEnd the side from which $w$ was reached before the procedure starts: 
	 * 0=prefix, 1=suffix;
	 * @param offset the offset between $transitiveReductionSource$ and $source$;
	 * @param largestOffset the largest offset between $transitiveReductionSource$ and any
	 * of its neighbors that use the same end of $transitiveReductionSource$ as $e$; 
	 * used for stopping reachability;
	 * @param shortestDistance assumed to be initialized to all $Math.POSITIVE_INFINITY$;
	 * the procedure leaves it in this state when it ends;
	 * @param isActive[i][j] represents an edge to $i$ from $u$, connected to the 
	 * prefix ($j=0$) or the suffix ($j=1$) of $i$;
	 * @param queue temporary space, assumed to be empty;
	 * @return the number of cells $(i,j)$ such that $isActive[i][j]$ is set to FALSE by 
	 * the procedure.
	 */
	private static final int reduce_impl(int transitiveReductionSource, int source, int sourceEnd, int offset, int largestOffset, int[][] shortestDistance, boolean[][] isActive, PriorityQueue<PriorityEnd> queue) {
		final int THRESHOLD = IO.quantum;
		int i;
		int from, to, currentStart, currentEnd, newStart, type, nDeactivated, lastReached, overhang;
		BidirectedEdge edge;
		PriorityEnd currentPriorityEnd, newPriorityEnd;
		PriorityEnd tmpPriorityEnd = getPriorityEnd();
		
		// Computing the earliest start of every node end reachable from $source$.
		shortestDistance[source][sourceEnd]=offset;
		newPriorityEnd=getPriorityEnd(); newPriorityEnd.set(source,sourceEnd,offset);
		queue.add(newPriorityEnd);
		stack[0]=source; stack[1]=sourceEnd; lastReached=1;
		while (!queue.isEmpty()) {
			currentPriorityEnd=queue.poll();
			currentStart=currentPriorityEnd.earliestStart;
			if (currentStart>largestOffset+THRESHOLD) break;
			from=currentPriorityEnd.node; currentEnd=currentPriorityEnd.end;
			for (i=0; i<=lastNeighbor[from]; i++) {
				edge=neighbors[from][i]; to=edge.getTo(from);
				if (to==source) continue;  // $to$ can be $transitiveReductionSource$ because we handle self-loops in $transitiveReductionSource$.
				type=edge.getType(from);
				if ((currentEnd==0 && type==Constants.OVERLAP_SUFFIX_PREFIX) || (currentEnd==1 && type==Constants.OVERLAP_PREFIX_PREFIX)) {
					newStart=currentStart+edge.getOverhang(from);
					if (newStart<shortestDistance[to][0]) {
						shortestDistance[to][0]=newStart;
						tmpPriorityEnd.node=to; tmpPriorityEnd.end=0;
						if (queue.contains(tmpPriorityEnd)) queue.remove(tmpPriorityEnd);
						else {
							stack[++lastReached]=to;
							stack[++lastReached]=0;
						}
						if (to!=transitiveReductionSource) {
							newPriorityEnd=getPriorityEnd();
							newPriorityEnd.set(to,0,newStart);
							queue.add(newPriorityEnd);
						}
					}
				}
				else if ((currentEnd==0 && type==Constants.OVERLAP_SUFFIX_SUFFIX) || (currentEnd==1 && type==Constants.OVERLAP_PREFIX_SUFFIX)) {
					newStart=currentStart+edge.getOverhang(from);
					if (newStart<shortestDistance[to][1]) {
						shortestDistance[to][1]=newStart;
						tmpPriorityEnd.node=to; tmpPriorityEnd.end=1;
						if (queue.contains(tmpPriorityEnd)) queue.remove(tmpPriorityEnd);
						else {
							stack[++lastReached]=to;
							stack[++lastReached]=1;
						}
						if (to!=transitiveReductionSource) {
							newPriorityEnd=getPriorityEnd();
							newPriorityEnd.set(to,1,newStart);
							queue.add(newPriorityEnd);
						}
					}
				}
			}
		}
		
		// Deactivating neighbors of $transitiveReductionSource$.
		nDeactivated=0;
		for (i=0; i<=lastNeighbor[transitiveReductionSource]; i++) {
			edge=neighbors[transitiveReductionSource][i];
			to=edge.getTo(transitiveReductionSource);
			if (to==source) continue;
			overhang=edge.getOverhang(transitiveReductionSource);
			if (isActive[to][0] && shortestDistance[to][0]<=overhang+THRESHOLD) {
				isActive[to][0]=false;
				nDeactivated++;
			}
			if (isActive[to][1] && shortestDistance[to][1]<=overhang+THRESHOLD) {
				isActive[to][1]=false;
				nDeactivated++;
			}
		}
		
		// Cleaning up for next call
		queue.clear();
		resetPriorityEndPool();  // The next call can reuse all elements in the pool.
		for (i=0; i<=lastReached; i+=2) shortestDistance[stack[i]][stack[i+1]]=Math.POSITIVE_INFINITY;
		
		return nDeactivated;
	}
	
	
	private static final PriorityEnd getPriorityEnd() {
		lastPriorityEnd++;
		if (lastPriorityEnd==priorityEndPool.length) {
			PriorityEnd[] priorityEndPoolPrime = new PriorityEnd[priorityEndPool.length<<1];
			System.arraycopy(priorityEndPool,0,priorityEndPoolPrime,0,lastPriorityEnd);
			for (int i=lastPriorityEnd; i<priorityEndPoolPrime.length; i++) priorityEndPoolPrime[i] = new PriorityEnd();
			priorityEndPool=priorityEndPoolPrime;
		}
		return priorityEndPool[lastPriorityEnd];
	}
	
	
	private static final void resetPriorityEndPool() {
		lastPriorityEnd=-1;
	}
	
	
	private static class PriorityEnd implements Comparable {
		public int node;
		public int end;
		public int earliestStart;
		
		public void set(int n, int e, int s) {
			node=n;
			end=e;
			earliestStart=s;
		}
		
		public boolean equals(Object other) {
			PriorityEnd otherEnd = (PriorityEnd)other;
			return node==otherEnd.node && end==otherEnd.end;
		} 
		
		public int compareTo(Object other) {
			PriorityEnd otherEnd = (PriorityEnd)other;
			if (earliestStart<otherEnd.earliestStart) return -1;
			else if (earliestStart>otherEnd.earliestStart) return 1;
			return 0;
		}
	}

		
	
	
	
	
	

	// --------------------- PRINTING AND SERIALIZATION PROCEDURES -----------------------
	
	/**
	 * Makes sure that all data structures required by $printAllPaths()$ are initialized.
	 */
	public static final void printAllPaths_init() {
		if (oneEnd==null || oneEnd.length<nNodes) oneEnd = new int[nNodes];
		if (oneEndFlags==null || oneEndFlags.length<nNodes) oneEndFlags = new boolean[nNodes];
		if (tmp==null) tmp = new int[4];
		getNodesWithOneEnd(oneEnd,oneEndFlags,tmp);
		nOneEnd=tmp[0];
		if (IO.SHOW_INTERACTIVE) System.err.println("nOneEnd="+nOneEnd+" nSingletons="+tmp[1]);
		if (stack==null || stack.length<nNodes*3) stack = new int[nNodes*3]; 
	}
	
	
	/**
	 * Prints the string label of every path from a one-end node to another one-end node, 
	 * and of every node with no neighbor, assuming that the bidirected graph has no 
	 * directed cycle.
	 *
	 * Remark: this procedure assembles all repeats that have been collapsed to the same
	 * kernel. This is important for reconstructing very long or very fragmented TEs, or 
	 * duplicated regions near the centromeres/telomeres (hundreds of kb), copied from 
	 * other chromosomes or from the same chromosome \cite{guy2000genomic}, which might be
	 * visible in the factorization. Such regions could contain other repeats (e.g. TEs), 
	 * but the latter should be assigned to different intervals and hopefully belong to 
	 * different kernels.
	 *
	 * Remark: the procedure also assumes that no node is marked as removed.
	 *
	 * @param inputPath file containing the strings that correspond to every node in the 
	 * bidirected graph, one string per line. Recall that nodes in the bidirected graph 
	 * are not necessarily identical to the nodes in the interval graph stored in 
	 * $intervalGraphPointers$: they could be e.g. slightly different substrings of the 
	 * same reads.
	 */
	public static final void printAllPaths(String inputPath, String outputPath) throws IOException {
		final int CAPACITY = 4096;  // In characters. Arbitrary.
		int i, c, pathID;
		BufferedReader br;
		BufferedWriter bw;
		CharStream[] strings;
		
		printAllPaths_init();
		strings=loadNodeLabels(inputPath);
		if (allPathsPrintedFrom==null || nNodes>allPathsPrintedFrom.length) allPathsPrintedFrom = new boolean[nNodes];
		Math.set(allPathsPrintedFrom,nNodes-1,false);
		if (charStack==null) charStack = new CharStream((CAPACITY)<<1);
		bw = new BufferedWriter(new FileWriter(outputPath),IO.BUFFER_SIZE);
		pathID=0;
		for (i=0; i<nOneEnd-1; i++) {  // No need to use the last one-end node as source
			pathID=printAllPathsFrom(oneEnd[i]>=0?oneEnd[i]:-1-oneEnd[i],oneEnd[i]>=0?1:0,oneEndFlags,allPathsPrintedFrom,strings,charStack,bw,pathID);
		}
		for (i=0; i<nNodes; i++) {
			if (lastNeighbor[i]==-1) {
				IO.writeFakeHeader(pathID++,strings[i].nCharacters(),null,bw);
				strings[i].print(bw,true,0);
				bw.write('\n');
			}
		}
		strings=null;
		bw.close();
	}
	
	
	/**
	 * Variant that stores in $assemblyPaths$ the distinct sequences of node IDs of every 
	 * path, and in $node2paths[v]$ the sorted list of all distinct path IDs that use node
	 * $v$. Every row of $assemblyPaths$ stores $x>=0$ (respectively, $-1-x < 0$) iff the 
	 * prefix (respectively, suffix) of node $x$ is connected to the (possibly absent) 
	 * previous node of the sequence. Rows of $assemblyPaths$ are not sorted in any
	 * specific way.
	 *
	 * @param firstTag path IDs in $node2paths$ start from this value;
	 * @param alreadyInitialized TRUE iff $printAllPaths_init()$ has already been called;
	 * @param flags array with one cell per node; 1=print all paths from the node; 
	 * 0=print a path that consists just of the node. The array is ignored if null.
	 */
	public static final void printAllPaths(int firstTag, boolean alreadyInitialized, int[] flags) {
		boolean equal;
		int i, j, k;
		int pathID, tmp, lengthI, lengthJ, nSingletons;
		int[] tmpArray;

		if (!alreadyInitialized) printAllPaths_init();
		
		// Adding paths
		if (allPathsPrintedFrom==null || nNodes>allPathsPrintedFrom.length) allPathsPrintedFrom = new boolean[nNodes];
		Math.set(allPathsPrintedFrom,nNodes-1,false);
		pathID=0;
		for (i=0; i<nOneEnd-1; i++) {  // No need to use the last one-end node as source
			if (flags!=null && flags[i]==0) continue;
			try { pathID=printAllPathsFrom(oneEnd[i]>=0?oneEnd[i]:-1-oneEnd[i],oneEnd[i]>=0?1:0,oneEndFlags,allPathsPrintedFrom,null,null,null,pathID); }
			catch (IOException e) { /* Never happens */ }
		}
		nSingletons=0;
		for (i=0; i<nNodes; i++) {
			if (lastNeighbor[i]==-1 || (flags!=null && flags[i]==0)) nSingletons++;
		}
		ensureAssemblyPaths(pathID+nSingletons);
		for (i=0; i<nNodes; i++) {
			if (lastNeighbor[i]==-1 || (flags!=null && flags[i]==0)) {
				if (assemblyPaths[pathID]==null || assemblyPaths[pathID].length==0) assemblyPaths[pathID] = new int[1];
				assemblyPaths[pathID][0]=i;
				assemblyPathLengths[pathID]=1;
				assemblyPathStringLengths[pathID]=nodeLength[i];
				pathID++;
			}
		}
		lastAssemblyPath=pathID-1;
		
		// Removing duplicated paths
		for (i=0; i<lastAssemblyPath; i++) {
			lengthI=assemblyPathLengths[i];
			if (lengthI<=1) continue;
			for (j=i+1; j<=lastAssemblyPath; j++) {
				lengthJ=assemblyPathLengths[j];
				if (lengthJ<=1 || lengthJ!=lengthI) continue;
				equal=true;
				for (k=0; k<lengthI; k++) {
					if (assemblyPaths[i][k]!=assemblyPaths[j][k]) {
						equal=false;
						break;
					}
				}
				if (equal) {
					assemblyPathLengths[j]=0;
					continue;
				}
				equal=true;
				for (k=0; k<lengthI; k++) {
					if (assemblyPaths[i][k]!=-1-assemblyPaths[j][lengthJ-1-k]) {
						equal=false;
						break;
					}
				}
				if (equal) assemblyPathLengths[j]=0;
			}
		}
		j=-1;
		for (i=0; i<=lastAssemblyPath; i++) {
			if (assemblyPathLengths[i]==0) continue;
			j++;
			if (j!=i) {
				tmpArray=assemblyPaths[j];
				assemblyPaths[j]=assemblyPaths[i];
				assemblyPaths[i]=tmpArray;
				tmp=assemblyPathLengths[j];
				assemblyPathLengths[j]=assemblyPathLengths[i];
				assemblyPathLengths[i]=tmp;
				tmp=assemblyPathStringLengths[j];
				assemblyPathStringLengths[j]=assemblyPathStringLengths[i];
				assemblyPathStringLengths[i]=tmp;
			}
		}
		lastAssemblyPath=j;
		
		// Assigning path tags to bidirected graph nodes
		printAllPaths_initNode2paths(firstTag);
	}
	
	
	/**
	 * @param firstTag path IDs in $node2paths$ start from this value.
	 */
	private static final void printAllPaths_initNode2paths(int firstTag) {
		final int GROWTH_RATE = 4;  // Arbitrary
		int i, j, k;
		int tmp;
		
		if (node2paths==null || node2paths.length<nNodes) {
			node2paths = new int[nNodes][GROWTH_RATE];
			node2paths_last = new int[nNodes];
		}
		Math.set(node2paths_last,nNodes-1,-1);
		for (i=0; i<=lastAssemblyPath; i++) {
			for (j=0; j<assemblyPathLengths[i]; j++) {
				k=assemblyPaths[i][j];
				if (k<0) k=-1-k;
				node2paths_last[k]++;
				if (node2paths[k]==null || node2paths_last[k]==node2paths[k].length) {
					int[] newArray = new int[node2paths[k].length+GROWTH_RATE];
					System.arraycopy(node2paths[k],0,newArray,0,node2paths[k].length);
					node2paths[k]=newArray;
				}
				node2paths[k][node2paths_last[k]]=firstTag+i;
			}
		}
		for (i=0; i<nNodes; i++) {
			if (node2paths_last[i]<=0) continue; 
			Arrays.sort(node2paths[i],0,node2paths_last[i]+1);
			k=0; tmp=node2paths[i][0];
			for (j=1; j<=node2paths_last[i]; j++) {
				if (node2paths[i][j]==tmp) continue;
				node2paths[i][++k]=node2paths[i][j];
				tmp=node2paths[i][j];
			}
			node2paths_last[i]=k;
		}
	}
	
	
	/**
	 * Variant that ignores edges and creates one distinct path per node.
	 *
	 * @param firstTag path IDs in $node2paths$ start from this value.
	 */
	public static final void printAllPaths_trivial(int firstTag) {
		int i;
		
		// Adding paths
		ensureAssemblyPaths(nNodes);
		for (i=0; i<nNodes; i++) {
			if (assemblyPaths[i]==null || assemblyPaths[i].length==0) assemblyPaths[i] = new int[1];
			assemblyPaths[i][0]=i;
			assemblyPathLengths[i]=1;
			assemblyPathStringLengths[i]=nodeLength[i];
		}
		lastAssemblyPath=nNodes-1;
		
		// Assigning path tags to bidirected graph nodes
		printAllPaths_initNode2paths(firstTag);
	}
	
	
	/**
	 * Variant of $printAllPaths()$ that stores in $assemblyPaths$ just one longest path
	 * per connected component of the bidirected graph.
	 *
	 * @param firstTag path IDs in $node2paths$ start from this value;
	 * @param componentPaths a longest path per connected component, as produced by
	 * $estimateAssemblies()$; all such paths are distinct;
	 * @param node2component component ID of every bidirected graph node, as produced by
	 * $getConnectedComponents()$.
	 */
	public static final void printLongestPaths(int firstTag, int nComponents, int[][] componentPaths, int[] componentLengths, int[] componentStringLengths, int[] node2component) {
		int i, length;
		
		printAllPaths_init();
		
		// Adding paths
		ensureAssemblyPaths(nComponents);
		System.arraycopy(componentStringLengths,0,assemblyPathStringLengths,0,nComponents);
		for (i=0; i<nComponents; i++) {
			length=componentLengths[i];
			if (assemblyPaths[i]==null || assemblyPaths[i].length<length) assemblyPaths[i] = new int[length];
			System.arraycopy(componentPaths[i],0,assemblyPaths[i],0,length);
			assemblyPathLengths[i]=length;
		}
		lastAssemblyPath=nComponents-1;
		
		// Assigning path tags to bidirected graph nodes
		if (node2paths==null || node2paths.length<nNodes) {
			node2paths = new int[nNodes][1];
			node2paths_last = new int[nNodes];
		}
		Math.set(node2paths_last,nNodes-1,0);
		for (i=0; i<nNodes; i++) {
			if (node2paths[i]==null || node2paths[i].length==0) node2paths[i] = new int[1];
			node2paths[i][0]=firstTag+node2component[i];
		}
	}
	
	
	/**
	 * Stores in $orientation$ the orientation of every bidirected graph node with respect
	 * to the longest paths computed by $printLongestPaths()$. 
	 *
	 * Remark: the procedure works also when the bidirected graph contains a cycle.
	 *
	 * @param tmpBoolean temporary space of size at least $nNodes$.
	 */
	public static final void propagateOrientationFromLongestPath(byte[] orientation, boolean[] tmpBoolean) {
		final int GROWTH_RATE = 100;  // Arbitrary
		byte oldOrientation, newOrientation;
		int i, j, k;
		int top, length, vertex, neighbor, type;
		BidirectedEdge edge;
		
		// Marking all nodes in a longest path
		Math.set(tmpBoolean,nNodes-1,false); Math.set(orientation,nNodes-1,(byte)(-1));
		for (i=0; i<=lastAssemblyPath; i++) {
			length=assemblyPathLengths[i];
			for (j=0; j<length; j++) {
				vertex=assemblyPaths[i][j];
				if (vertex<0) {
					vertex=-1-vertex;
					orientation[vertex]=1;
				}
				else orientation[vertex]=0;
				tmpBoolean[vertex]=true;
			}
		}
		
		// Propagating
		if (stack==null || stack.length<nNodes) stack = new int[nNodes];
		for (i=0; i<=lastAssemblyPath; i++) {
			length=assemblyPathLengths[i];
			for (j=0; j<length; j++) {
				stack[0]=assemblyPaths[i][j]; top=0;
				while (top>=0) {
					vertex=stack[top--];
					if (vertex<0) vertex=-1-vertex;
					for (k=0; k<=lastNeighbor[vertex]; k++) {
						edge=neighbors[vertex][k];
						neighbor=edge.getTo(vertex);
						if (tmpBoolean[neighbor]) continue;
						type=edge.getType(vertex);
						newOrientation=-1;
						if (orientation[vertex]==0) {
							switch (type) {
								case Constants.OVERLAP_PREFIX_PREFIX: newOrientation=1; break;
								case Constants.OVERLAP_PREFIX_SUFFIX: newOrientation=0; break;
								case Constants.OVERLAP_SUFFIX_PREFIX: newOrientation=0; break;
								case Constants.OVERLAP_SUFFIX_SUFFIX: newOrientation=1; break;
							}
						}
						else {
							switch (type) {
								case Constants.OVERLAP_PREFIX_PREFIX: newOrientation=0; break;
								case Constants.OVERLAP_PREFIX_SUFFIX: newOrientation=1; break;
								case Constants.OVERLAP_SUFFIX_PREFIX: newOrientation=1; break;
								case Constants.OVERLAP_SUFFIX_SUFFIX: newOrientation=0; break;
							}
						}
						oldOrientation=orientation[neighbor];
						if (oldOrientation==-1) orientation[neighbor]=newOrientation;
						else orientation[neighbor]=IntervalGraphStep3.addOrientation[oldOrientation][newOrientation];
						if (orientation[neighbor]!=oldOrientation) {
							top++;
							if (top>=stack.length) {
								int[] newStack = new int[stack.length+GROWTH_RATE];
								System.arraycopy(stack,0,newStack,0,stack.length);
								stack=newStack;
							}
							stack[top]=neighbor;
						}
					}
				}
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<nNodes; i++) {
				if (orientation[i]==-1) {
					System.err.println("propagateLongestPathOrientations> ERROR: the "+i+"-th node of the bidirected graph does not have an orientation?!");
					System.err.println("lastAssemblyPath="+lastAssemblyPath+" node lastNeighbor="+lastNeighbor[i]);
					System.err.println("assembly paths:");
					for (int x=0; x<=lastAssemblyPath; x++) {
						for (int y=0; y<assemblyPathLengths[x]; y++) System.err.print(assemblyPaths[x][y]+",");
						System.err.println();
					}
					System.exit(1);
				}
			}
		}
	}
	
	
	/**
	 * Initializes assembly path data structures when the bidirected graph contains just 
	 * one node.
	 */
	public static final void printAllPaths_singleton(int firstTag) {
		lastAssemblyPath=0;
		if (assemblyPaths==null || assemblyPaths.length==0) assemblyPaths = new int[1][1];
		if (assemblyPaths[0]==null || assemblyPaths[0].length==0) assemblyPaths[0] = new int[1];
		assemblyPaths[0][0]=0;
		if (assemblyPathLengths==null || assemblyPathLengths.length==0) assemblyPathLengths = new int[1];
		assemblyPathLengths[0]=1;
		if (assemblyPathStringLengths==null || assemblyPathStringLengths.length==0) assemblyPathStringLengths = new int[1];
		assemblyPathStringLengths[0]=nodeLength[0];
		if (node2paths==null || node2paths.length==0) node2paths = new int[1][1];
		if (node2paths[0]==null || node2paths[0].length==0) node2paths[0] = new int[1];
		node2paths[0][0]=firstTag;
		if (node2paths_last==null || node2paths_last.length==0) node2paths_last = new int[1];
		node2paths_last[0]=0;
	}
	
	
	/**
	 * Computes the number of assembly paths that have the first position (if 
	 * $startOrEnd=TRUE$) or the last position (if $startOrEnd=FALSE$) of $node$ at one of
	 * their ends. An assembly path is counted twice if the requested end of $node$ occurs
	 * both at the beginning and at the end of a path.
	 *
	 * @param firstTag path IDs in $node2paths$ start from this value.
	 */
	public static final int nAssemblyPaths(int node, int firstTag, boolean startOrEnd) {
		int i;
		int path, last, out;
		
		out=0;
		last=node2paths_last[node];
		if (startOrEnd) {
			for (i=0; i<=last; i++) {
				path=node2paths[node][i]-firstTag;
				if (assemblyPaths[path][0]==node) out++;
				if (assemblyPaths[path][assemblyPathLengths[path]-1]==-1-node) out++;
			}
		}
		else {
			for (i=0; i<=last; i++) {
				path=node2paths[node][i]-firstTag;
				if (assemblyPaths[path][0]==-1-node) out++;
				if (assemblyPaths[path][assemblyPathLengths[path]-1]==node) out++;
			}
		}
		return out;
	}
	
	
	/**
	 * Like $nAssemblyPaths()$, but stores in $out$ the sorted IDs of all paths.
	 * A path ID $x>=0$ means that the requested end of $node$ appears at the beginning of
	 * path $x$; a path ID $x=-1-y < 0$ means that the requested end of $node$ occurs at 
	 * the end of path $y$. Beginning and end of a path are based on the sequences of node
	 * IDs in $assemblyPaths$, and are arbitrary with respect to the string a path 
	 * represents.
	 *
	 * @param firstTag path IDs in $node2paths$ start from this value;
	 * @param out assumed to be large enough to contain all path IDs.
	 */
	public static final void getAssemblyPaths(int node, int firstTag, boolean startOrEnd, int[] out) {
		int i, j;
		int path, last;
		
		last=node2paths_last[node];	
		j=-1;
		if (startOrEnd) {
			for (i=0; i<=last; i++) {
				path=node2paths[node][i]-firstTag;
				if (assemblyPaths[path][0]==node) out[++j]=firstTag+path;
				if (assemblyPaths[path][assemblyPathLengths[path]-1]==-1-node) out[++j]=-1-(firstTag+path);
			}
		}
		else {
			for (i=0; i<=last; i++) {
				path=node2paths[node][i]-firstTag;
				if (assemblyPaths[path][0]==-1-node) out[++j]=firstTag+path;
				if (assemblyPaths[path][assemblyPathLengths[path]-1]==node) out[++j]=-1-(firstTag+path);
			}
		}
		if (j>0) Arrays.sort(out,0,j+1);
	}
	
	
	/**
	 * Stores in $out[i]$ the orientation in which $node$ is used in path 
	 * $node2paths[node][i]$: 0=forward, 1=RC, 2=both.
	 *
	 * @param firstTag path IDs in $node2paths$ start from this value.
	 */
	public static final void getAssemblyPathOrientations(int node, int firstTag, byte[] out) {
		byte orientation;
		int i, j;
		int path, last, length;
		
		last=node2paths_last[node];
		for (i=0; i<=last; i++) {
			path=node2paths[node][i]-firstTag;
			orientation=-1;
			length=assemblyPathLengths[path];
			for (j=0; j<length; j++) {
				if (assemblyPaths[path][j]==node) {
					if (orientation==-1) orientation=0;
					else if (orientation==0) { /* NOP */ }
					else if (orientation==1) orientation=2;
				}
				else if (assemblyPaths[path][j]==-1-node) {
					if (orientation==-1) orientation=1;
					else if (orientation==0) orientation=2;
					else if (orientation==1) { /* NOP */ }
				}
			}
			out[i]=orientation;
		}
	}
	
	
	/**
	 * Prints to stream $bw$ the string label of all paths from $source$, which is assumed
	 * to be a one-end node, to every other one-end node. If $bw=null$, the procedure
	 * instead lists the sequence of nodes in every path, in the global data structure
	 * $assemblyPaths$ (such sequences might not be all distinct). Every such sequence
	 * stores $x>=0$ (respectively, $-1-x < 0$) if the prefix (respectively, suffix) of 
	 * node $x$ is connected to the (possibly absent) previous node of the sequence.
	 *
	 * Remark: the procedure prints twice every path $source-...-RC(source)$ that uses the 
	 * same edge from $source$ (once in the forward and once in the RC orientation). Every
	 * other path $source-...-RC(source)$ is printed exactly once.
	 *
	 * Remark: the procedure uses global variable $stack$.
	 *
	 * @param sourceEnd the side of $source$ not involved in overlaps;
	 * @param oneEndFlags for each node, tells whether it is used from just one side;
	 * @param allPathsPrintedFrom for each node, tells whether the procedure has been 
	 * already executed from that node (the procedure marks cell $source$ at the end);
	 * @param strings (used only if $bw!=null$) strings that correspond to every node in 
	 * the bidirected graph;
	 * @param charStack (used only if $bw!=null$) temporary space;
	 * @param pathID first ID to assign to a path;
	 * @return first available path ID after the procedure completes.
	 */
	private static final int printAllPathsFrom(int source, int sourceEnd, boolean[] oneEndFlags, boolean[] allPathsPrintedFrom, CharStream[] strings, CharStream charStack, BufferedWriter bw, int pathID) throws IOException {
		final boolean printStrings = bw!=null;
		boolean pushed, doContinue;
		int i, j, k;
		int top, type, first, last, from, fromEnd, to, toEnd, toOverhang, toLength;
		final int sourceLength;
		BidirectedEdge edge;
		
		sourceLength=nodeLength[source];
		if (!printStrings) ensureAssemblyPaths(pathID+1);
		stack[0]=sourceLength;  // String delta
		stack[1]=-1;  // Last visited neighbor
		stack[2]=sourceEnd==0?source:-1-source;  // Positive=prefix, negative=suffix.
		if (printStrings) charStack.push(strings[source],0,sourceLength-1,sourceEnd);
		top=2;
		while (top>=0) {
			fromEnd=stack[top]>=0?0:1;
			from=stack[top]>=0?stack[top]:-1-stack[top];
			if ((from!=source || fromEnd!=sourceEnd) && oneEndFlags[from]) {
				if (!allPathsPrintedFrom[from]) {
					if (printStrings) {
						IO.writeFakeHeader(pathID,charStack.nCharacters(),null,bw);
						charStack.print(bw,true,0);
						bw.write('\n');
					}
					else addAssemblyPath(pathID,top);
					pathID++;
				}
				if (printStrings) charStack.pop(stack[top-2],false);
				top-=3;
				continue;
			}
			pushed=false;
			for (j=stack[top-1]+1; j<=lastNeighbor[from]; j++) {
				edge=neighbors[from][j];
				if (edge.isRedundant) continue;
				to=edge.getTo(from);
				if (to==source) {
					doContinue=false;
					for (k=0; k<=lastNeighbor[source]; k++) {
						if (neighbors[source][k]==edge) {
							if (k<stack[1]) doContinue=true;
							break;
						}
					}
					if (doContinue) continue;
				}
				toEnd=-1;
				type=edge.getType(from);
				if (fromEnd==0) {
					if (type==Constants.OVERLAP_SUFFIX_PREFIX) toEnd=0;
					else if (type==Constants.OVERLAP_SUFFIX_SUFFIX) toEnd=1;
				}
				else {
					if (type==Constants.OVERLAP_PREFIX_PREFIX) toEnd=0;
					else if (type==Constants.OVERLAP_PREFIX_SUFFIX) toEnd=1;
				}
				if (toEnd==-1) continue;
				stack[top-1]=j;
				toOverhang=edge.getOverhang(to);
				if (top+3>=stack.length) {
					int[] newStack = new int[stack.length<<1];
					System.arraycopy(stack,0,newStack,0,top+1);
					stack=newStack;
				}
				stack[++top]=toOverhang;
				stack[++top]=-1;
				stack[++top]=toEnd==0?to:-1-to;
				if (printStrings) {
					toLength=strings[to].nCharacters();
					if (toEnd==0) {
						first=toLength-toOverhang;
						last=toLength-1;
					}
					else {
						first=0;
						last=toOverhang-1;
					}
					charStack.push(strings[to],first,last,toEnd);
				}
				pushed=true;
				break;
			}
			if (!pushed) {
				if (printStrings) charStack.pop(stack[top-2],false);
				top-=3;
			}
		}
		allPathsPrintedFrom[source]=true;
		return pathID;
	}
	
	
	/**
	 * @param visited temporary matrix with at least $nNodes$ rows and 2 columns; the
	 * procedure assumes it to be all zeros at the beginning, and resets it at the end;
	 * @return TRUE iff there is a path from end $sourceEnd$ of $source$, that uses the
	 * same end of the same node more than once.
	 */
	public static final boolean cycleFrom(int source, int sourceEnd, boolean[][] visited) {
		boolean pushed;
		int i, j, k;
		int top, from, fromEnd, to, toEnd, toPrime, type;
		BidirectedEdge edge;
		
		stack[0]=-1;  // Last visited neighbor
		stack[1]=sourceEnd==0?-1-source:source;  // Positive=prefix, negative=suffix.
		visited[source][sourceEnd]=true;
		top=1;
		while (top>0) {
			fromEnd=stack[top]>=0?0:1;
			from=stack[top]>=0?stack[top]:-1-stack[top];
			pushed=false;
			for (j=stack[top-1]+1; j<=lastNeighbor[from]; j++) {
				edge=neighbors[from][j];
				if (edge.isRedundant) continue;				
				to=edge.getTo(from);
				toEnd=-1;
				type=edge.getType(from);
				if (fromEnd==0) {
					if (type==Constants.OVERLAP_SUFFIX_PREFIX) toEnd=0;
					else if (type==Constants.OVERLAP_SUFFIX_SUFFIX) toEnd=1;
				}
				else {
					if (type==Constants.OVERLAP_PREFIX_PREFIX) toEnd=0;
					else if (type==Constants.OVERLAP_PREFIX_SUFFIX) toEnd=1;
				}
				if (toEnd==-1) continue;
				if (visited[to][toEnd]) return true;
				toPrime=toEnd==0?to:-1-to;
				stack[top-1]=j;
				if (top+2>=stack.length) {
					int[] newStack = new int[stack.length<<1];
					System.arraycopy(stack,0,newStack,0,top+1);
					stack=newStack;
				}
				stack[++top]=-1;
				stack[++top]=toPrime;
				visited[to][toEnd]=true;
				pushed=true;
				break;
			}
			if (!pushed) {
				top-=2;
				visited[from][fromEnd]=false;
			}
		}
		return false;
	}
	
	
	private static final void ensureAssemblyPaths(int nPaths) {
		final int GROWTH_RATE = 16;  // Arbitrary
		
		if (assemblyPaths!=null && nPaths<=assemblyPaths.length) return;
		int[][] newAssemblyPaths = new int[nPaths+GROWTH_RATE][0];
		if (assemblyPaths!=null) System.arraycopy(assemblyPaths,0,newAssemblyPaths,0,assemblyPaths.length);
		assemblyPaths=newAssemblyPaths;
		int[] newAssemblyPathLengths = new int[nPaths+GROWTH_RATE];
		if (assemblyPathLengths!=null) System.arraycopy(assemblyPathLengths,0,newAssemblyPathLengths,0,assemblyPathLengths.length);
		assemblyPathLengths=newAssemblyPathLengths;
		int[] newAssemblyPathStringLengths = new int[nPaths+GROWTH_RATE];
		if (assemblyPathStringLengths!=null) System.arraycopy(assemblyPathStringLengths,0,newAssemblyPathStringLengths,0,assemblyPathStringLengths.length);
		assemblyPathStringLengths=newAssemblyPathStringLengths;
	}
	
	
	/**
	 * Remark: the procedure assumes that global variable $stack$ is the DFS stack used by
	 * $printAllPathsFrom()$.
	 *
	 * @param top last used element in $stack$.
	 */
	private static final void addAssemblyPath(int pathID, int top) {
		final int length = (top+1)/3;
		int i;
		
		ensureAssemblyPaths(pathID+1);
		if (assemblyPaths[pathID]==null || assemblyPaths[pathID].length<length) assemblyPaths[pathID] = new int[length];
		for (i=0; i<length; i++) assemblyPaths[pathID][i]=stack[2+3*i];
		assemblyPathLengths[pathID]=length;
		assemblyPathStringLengths[pathID]=0;
		for (i=0; i<length; i++) assemblyPathStringLengths[pathID]+=stack[3*i];
	}
	
	
	/**
	 * Loads the string labels of all nodes.
	 */
	private static final CharStream[] loadNodeLabels(String inputPath) throws IOException {
		final int CAPACITY = 4096;  // In characters. Arbitrary.
		int i, c;
		BufferedReader br;
		
		if (strings==null) {
			strings = new CharStream[nNodes];
			for (i=0; i<nNodes; i++) strings[i] = new CharStream(CAPACITY);
		}
		else if (nNodes>strings.length) {
			CharStream[] newStrings = new CharStream[nNodes];
			System.arraycopy(strings,0,newStrings,0,strings.length);
			for (i=strings.length; i<nNodes; i++) newStrings[i] = new CharStream(CAPACITY);
			strings=newStrings;
		}
		br = new BufferedReader(new FileReader(inputPath));
		for (i=0; i<nNodes; i++) {
			strings[i].load(nodeLength[i],br);
			c=br.read();
			if (IO.CONSISTENCY_CHECKS && c!=-1 && c!='\n') {
				System.err.println("loadNodeLabels> ERROR: wrong input file format.");
				System.exit(1);
			}
		}
		br.close();
		return strings;
	}
	
	
	/**
	 * Prints the string label of the $pathID$-th path in $assemblyPaths$.
	 *
	 * @param strings the strings that correspond to every node in the bidirected graph;
	 * if NULL, they are loaded from $stringsFile$;
	 * @param stringsFile one string per line;
	 * @param outputString path of the output file; if NULL, the procedure does not create
	 * an output file;
	 * @param str temporary buffer for concatenating strings; at the end of the procedure
	 * it contains the full string of the path;
	 * @param prefixLength if >0, prints just a prefix of this length; 
	 * @return length of the printed string.
	 */
	public static final int printPath(int pathID, CharStream[] strings, String stringsFile, String outputString, CharStream str, int prefixLength, int header) throws IOException {
		final int CAPACITY = 4096;  // In characters. Arbitrary.
		int i, j;
		int end, from, fromEnd, to, toEnd, toOverhang, toLength, first, last, nCharacters;
		BidirectedEdge edge;
		BufferedWriter bw;
		
		if (strings==null) strings=loadNodeLabels(stringsFile);
		str.ensureSize(CAPACITY);
		str.clear(false);
		fromEnd=assemblyPaths[pathID][0];
		if (fromEnd>=0) { from=fromEnd; end=0; }
		else { from=-1-fromEnd; end=1; }
		str.push(strings[from],0,strings[from].nCharacters()-1,end);
		for (i=1; i<assemblyPathLengths[pathID]; i++) {
			if (prefixLength>0 && str.nCharacters()>=prefixLength) break;
			toEnd=assemblyPaths[pathID][i];
			if (toEnd>=0) { to=toEnd; end=0; }
			else { to=-1-toEnd; end=1; }
			for (j=0; j<=lastNeighbor[from]; j++) {
				edge=neighbors[from][j];
				if (edge.getTo(from)!=to) continue;
				toOverhang=edge.getOverhang(to);
				toLength=strings[to].nCharacters();
				if (end==0) {
					first=toLength-toOverhang;
					last=toLength-1;
				}
				else {
					first=0;
					last=toOverhang-1;
				}
				str.push(strings[to],first,last,end);
				break;
			}
			from=to;
		}
		if (outputString!=null) {
			bw = new BufferedWriter(new FileWriter(outputString),IO.BUFFER_SIZE);
			nCharacters=str.nCharacters();
			IO.writeFakeHeader(header,((prefixLength>0)&&(prefixLength<nCharacters))?prefixLength:nCharacters,null,bw);
			str.print(bw,true,prefixLength);
			bw.write('\n');
			bw.close();
		}
		return prefixLength>0?prefixLength:str.nCharacters();
	}
	
	
	/**
	 * Does not store $intervalGraphPointers$ and $removed$.
	 */
	public static final void serialize(String path) throws IOException {
		int i, j, k;
		int otherNode;
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(path),IO.BUFFER_SIZE);
		
		// Writing nodes and arcs
		IO.writeInt(nNodes,out);
		for (i=0; i<nNodes; i++) IO.writeInt(nodeLength[i],out);
		cleanMarkedFlags(false);
		for (i=0; i<nNodes; i++) {
			IO.writeInt(lastNeighbor[i]+1,out);
			for (j=0; j<=lastNeighbor[i]; j++) {
				if (neighbors[i][j].marked) {
					IO.writeInt(1,out);
					otherNode=neighbors[i][j].getTo(i);
					IO.writeInt(otherNode,out);
					for (k=0; k<=lastNeighbor[otherNode]; k++) {
						if (neighbors[otherNode][k]==neighbors[i][j]) {
							IO.writeInt(k,out);
							break;
						}
					}
				}
				else {
					IO.writeInt(0,out);
					neighbors[i][j].serialize(out);
					neighbors[i][j].marked=true;
				}
			}
		}
		
		// Writing assembly paths
		IO.writeInt(lastAssemblyPath+1,out);
		for (i=0; i<=lastAssemblyPath; i++) {
			IO.writeInt(assemblyPathLengths[i],out);
			IO.writeInt(assemblyPathStringLengths[i],out);
			for (j=0; j<assemblyPathLengths[i]; j++) IO.writeInt(assemblyPaths[i][j],out);
		}
		
		out.close();
	}
	
	
	/**
	 * Does not initialize $intervalGraphPointers$, and sets $removed$ to all zeros.
	 */
	public static final void deserialize(String path) throws IOException {
		int i, j, k;
		int otherNode, stored;
		BufferedInputStream in = new BufferedInputStream(new FileInputStream(path),IO.BUFFER_SIZE);
		
		// Reading nodes and arcs
		nNodes=IO.readInt(in);
		nodeLength = new int[nNodes];
		for (i=0; i<nNodes; i++) nodeLength[i]=IO.readInt(in);
		neighbors = new BidirectedEdge[nNodes][0];
		lastNeighbor = new int[nNodes];
		for (i=0; i<nNodes; i++) {
			lastNeighbor[i]=IO.readInt(in)-1;
			neighbors[i] = new BidirectedEdge[lastNeighbor[i]+1];
			for (j=0; j<=lastNeighbor[i]; j++) {
				stored=IO.readInt(in);
				if (stored==1) {
					otherNode=IO.readInt(in);
					k=IO.readInt(in);
					neighbors[i][j]=neighbors[otherNode][k];
				}
				else {
					neighbors[i][j] = new BidirectedEdge();
					neighbors[i][j].deserialize(in);
				}
			}
		}
		removed = new boolean[nNodes];
		Math.set(removed,nNodes-1,false);
		
		// Reading assembly paths
		lastAssemblyPath=IO.readInt(in)-1;
		assemblyPathLengths = new int[lastAssemblyPath+1];
		assemblyPathStringLengths = new int[lastAssemblyPath+1];
		assemblyPaths = new int[lastAssemblyPath+1][0];
		for (i=0; i<=lastAssemblyPath; i++) {
			assemblyPathLengths[i]=IO.readInt(in);
			assemblyPathStringLengths[i]=IO.readInt(in);
			assemblyPaths[i] = new int[assemblyPathLengths[i]];
			for (j=0; j<assemblyPathLengths[i]; j++) assemblyPaths[i][j]=IO.readInt(in);
		}
		
		in.close();
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static final void print() {
		int i, j;
		
		System.err.println("Bidirected graph (nNodes="+nNodes+", nRemoved="+nRemoved+"):");
		for (i=0; i<nNodes; i++) {
			System.err.print(i+" ("+(removed[i]?"1":"0")+","+nodeLength[i]+"): ");
			for (j=0; j<=lastNeighbor[i]; j++) System.err.print((neighbors[i][j]==null?"null":neighbors[i][j].getTo(i))+",");
			System.err.println();
		}
	}
	
	
	/**
	 * Prints the sorted and compacted list of all overhangs of node $node$.
	 *
	 * @param out temporary space, of size at least equal to the number of neighbors of
	 * $node$.
	 */
	public static final void printOverhangs(int node, Point[] out) {
		int i, last;
		
		last=lastNeighbor[node];
		for (i=0; i<=last; i++) {
			out[i].position=neighbors[node][i].getOverhang(node);
			out[i].mass=1;
		}
		last=Points.sortAndCompact(out,last);
		for (i=0; i<=last; i++) System.out.println(out[i]);
	}
	
	
	/**
	 * Remark: the procedure does not print nodes that are marked as removed.
	 *
	 * @param matrix with at least 4 columns: 0=readID, 1=start, 2=end, 3=type; 
	 * if NULL, $intervalGraphPointer$ is used instead.
	 */
	public static final void toDot(String path, int[][] matrix) throws IOException {
		int i, j, t;
		String type, fillColor;
		BufferedWriter bw = new BufferedWriter(new FileWriter(path),IO.BUFFER_SIZE);
		
		cleanMarkedFlags(true);
		bw.write("digraph G {\n");
		bw.write("bgcolor=\""+Colors.COLOR_DOT_BACKGROUND+"\";\n");
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			t=matrix!=null?matrix[i][3]:intervalGraphPointers[i].type;
			if (t==Constants.INTERVAL_DENSE_PREFIX) type="P // ";
			else if (t==Constants.INTERVAL_DENSE_SUFFIX) type="S // ";
			else if (t==Constants.INTERVAL_DENSE_PREFIXSUFFIX) type="PS // ";
			else if (t==Constants.INTERVAL_PERIODIC && intervalGraphPointers[i].hasLongPeriod) type="L // ";
			else type="";
			fillColor=Colors.type2color_dot(t);
			bw.write(i+" [shape=\"plaintext\",label=\""+type+(matrix!=null?matrix[i][0]:intervalGraphPointers[i].read)+" ["+(matrix!=null?matrix[i][1]:intervalGraphPointers[i].start)+".."+(matrix!=null?matrix[i][2]:intervalGraphPointers[i].end)+"] // "+intervalGraphPointers[i].nodeID+(intervalGraphPointers[i].isWeak?",W":"")+"\",style=\"filled\",color=\""+fillColor+"\",fontcolor=\""+Colors.COLOR_DOT_TEXT+"\",fontname=\"Helvetica\",fontsize=\"15\"];\n");
		}
		for (i=0; i<nNodes; i++) {
			if (removed[i]) continue;
			for (j=0; j<=lastNeighbor[i]; j++) {
				if (!neighbors[i][j].marked) {
					bw.write(neighbors[i][j].toDot());
					neighbors[i][j].marked=true;
				}
			}
		}
		bw.write("}\n");
		bw.close();
		cleanMarkedFlags(true);
	}
	

	
	
	
	
	
		
	// -------------------------------- DATA STRUCTURES ----------------------------------
	
	/**
	 * An edge of the bidirected graph
	 */
	public static class BidirectedEdge implements Comparable {
		/**
		 * Sorting order:
		 *
		 * $node1$: edges are sorted by $node2$;
		 * $node2$: edges are sorted by $node1$;
		 * $-1-Z$: edges with $inSubgraphXY=true$, where $nodeX=Z$, precede the others.
 		 */
		public static int order;
		
		/**
		 * Properties of an edge
		 */
		public int node1, node2;
		public int type;  // Exactly one of the $OVERLAP_*$ constants
		public int overhang1, overhang2;  // Minima
		public int maxOverhang1, maxOverhang2;  // Maxima
		public boolean isRedundant;  // For transitive reduction
		public int nAlignments;
		
		/**
		 * Temporary space for subgraph construction/traversal
		 */
		public boolean inSubgraph12, inSubgraph21;
		public boolean traversed12, traversed21;
		public boolean marked;
		
		
		/**
		 * @return average of the lengths in the two nodes.
		 */
		public int getOverlapLength() {
			return (nodeLength[node1]-overhang1+nodeLength[node2]-overhang2+2)>>1;
		}
		
		
		public void clone(BidirectedEdge from) {
			node1=from.node1;
			node2=from.node2;
			type=from.type;
			overhang1=from.overhang1;
			overhang2=from.overhang2;
			maxOverhang1=from.maxOverhang1;
			maxOverhang2=from.maxOverhang2;
			isRedundant=from.isRedundant;
			nAlignments=from.nAlignments;
		}
		
		
		public int compareTo(Object other) {
			BidirectedEdge otherEdge = (BidirectedEdge)other;
			if (order>=0) {
				int n = order==node1?node2:node1;
				int otherN = order==otherEdge.node1?otherEdge.node2:otherEdge.node1;
				if (n<otherN) return -1;
				else if (n>otherN) return 1;
				return 0;
			}
			else {
				boolean s = order==-1-node1?inSubgraph12:inSubgraph21;
				boolean otherS = order==-1-otherEdge.node1?otherEdge.inSubgraph12:otherEdge.inSubgraph21;
				if (s && !otherS) return -1;
				else if (!s && otherS) return 1;
				return 0;
			}
		}
		
		
		public final int getTo(int from) {
			return from==node1?node2:node1;
		}
		
		
		/**
		 * @return the overlap type as seen when going from node $from$ to the other node
		 * (recall that $type$ is the overlap seen when going from $node1$ to $node2$).
		 */
		public final int getType(int from) {
			if (from==node1) return type;
			else if (from==node2) return Constants.oppositeDirectionOverlap(type);
			return -1;
		}

		
		public final int getOverhang(int node) {
			return node==node1?overhang1:overhang2;
		}
		
		
		/**
		 * @param overhangs rows: overhang types; column 0: node1; column 1: node2.
		 */
		public final void setOverhangs(int type, int[][] overhangs) {
			if (type==Constants.OVERLAP_PREFIX_PREFIX) {
				overhang1=overhangs[0][0];
				overhang2=overhangs[0][1];
			}
			else if (type==Constants.OVERLAP_PREFIX_SUFFIX) {
				overhang1=overhangs[1][0];
				overhang2=overhangs[1][1];
			}
			else if (type==Constants.OVERLAP_SUFFIX_PREFIX) {
				overhang1=overhangs[2][0];
				overhang2=overhangs[2][1];
			}
			else if (type==Constants.OVERLAP_SUFFIX_SUFFIX) {
				overhang1=overhangs[3][0];
				overhang2=overhangs[3][1];
			}
		}
		
		
		/**
		 * Copies the overhangs of $other$ iff the length of the string that results from
		 * the overlap is smaller than the length of the string that results from the 
		 * current overlap.
		 *
		 * Remark: the procedure assumes that the type of $other$ and of this edge are
		 * identical up to swapping $node1$ and $node2$.
		 */
		private final void copyOverhangs(BidirectedEdge other, int[] nodeLength) {
			final double length = overhang1+((nodeLength[node1]-overhang1+nodeLength[node2]-overhang2)/2.0)+overhang2;
			final double otherLength = other.overhang1+((nodeLength[other.node1]-other.overhang1+nodeLength[other.node2]-other.overhang2)/2.0)+other.overhang2;
			if (otherLength<length) {
				overhang1=other.node1==node1?other.overhang1:other.overhang2;
				overhang2=other.node2==node2?other.overhang2:other.overhang1;
			}
		}
		
		
		/**
		 * Stores in $interval$ the relative interval of node $nodeID$ that is not covered
		 * by the overlap, in the forward direction of node $nodeID$.
		 */
		public final void getNonoverlappingInterval(int nodeID, int[] interval) {
			final int length = nodeLength[nodeID];
			
			if (nodeID==node1) {
				if (type==Constants.OVERLAP_PREFIX_PREFIX || type==Constants.OVERLAP_PREFIX_SUFFIX) {
					interval[0]=length-overhang1;
					interval[1]=length-1;
				}
				else if (type==Constants.OVERLAP_SUFFIX_PREFIX || type==Constants.OVERLAP_SUFFIX_SUFFIX) {
					interval[0]=0;
					interval[1]=overhang1-1;
				}
			}
			else if (nodeID==node2) {
				if (type==Constants.OVERLAP_PREFIX_PREFIX || type==Constants.OVERLAP_SUFFIX_PREFIX) {
					interval[0]=length-overhang2;
					interval[1]=length-1;
				}
				else if (type==Constants.OVERLAP_PREFIX_SUFFIX || type==Constants.OVERLAP_SUFFIX_SUFFIX) {
					interval[0]=0;
					interval[1]=overhang2-1;
				}
			}
		}
		
		
		// ------------------------ TEMPORARY SPACE PROCEDURES ---------------------------
		
		/**
		 * Tells whether the edge has been traversed in the direction $(from,w)$, where 
		 * $w$ is the node different from $from$, by the subgraph construction procedure.
		 */
		public final boolean inSubgraph(int from) {
			return from==node1?inSubgraph12:inSubgraph21;
		}
		
		
		public final void setSubgraph(int from, boolean b) {
			if (from==node1) inSubgraph12=b;
			else if (from==node2) inSubgraph21=b;
		}
		
		
		/**
		 * Tells whether the edge has been traversed in the direction that starts from
		 * $from$ and ends in the node different from $from$.
		 */
		public final boolean isTraversed(int from) {
			return from==node1?traversed12:traversed21;
		}
		
		
		public final void setTraversed(int from, boolean b) {
			if (from==node1) traversed12=b;
			else if (from==node2) traversed21=b;
		}
		
		
		// -------------------- PRINTING AND SERIALIZATION PROCEDURES --------------------
		
		private final void serialize(BufferedOutputStream out) throws IOException {
			IO.writeInt(node1,out);
			IO.writeInt(node2,out);
			IO.writeInt(type,out);
			IO.writeInt(overhang1,out);
			IO.writeInt(overhang2,out);
			IO.writeInt(maxOverhang1,out);
			IO.writeInt(maxOverhang2,out);
			IO.writeInt(isRedundant?1:0,out);
			IO.writeInt(nAlignments,out);
		}
		
		
		private final void deserialize(BufferedInputStream in) throws IOException {
			node1=IO.readInt(in);
			node2=IO.readInt(in);
			type=IO.readInt(in);
			overhang1=IO.readInt(in);
			overhang2=IO.readInt(in);
			maxOverhang1=IO.readInt(in);
			maxOverhang2=IO.readInt(in);
			isRedundant=IO.readInt(in)==1;
			nAlignments=IO.readInt(in);
		}
		
		
		public String toString() {
			return "("+node1+","+overhang1+"), ("+node2+","+overhang2+"), "+type+", "+isRedundant+" | nAlignments="+nAlignments+", marked="+marked+" | inSubgraph12="+inSubgraph12+" inSubgraph21="+inSubgraph21;
		}
		
		
		public String toDot() {
			final double WIDTH_SCALE = 4.0;
			final int PEN_WIDTH = Math.min((int)(nAlignments*WIDTH_SCALE),10);  // Arbitrary. Avoids printing huge edges.
			return node1+" -> "+node2+" [dir=\"both\",arrowtail=\""+getArrow(node1)+"\",arrowhead=\""+getArrow(node2)+"\",style=\""+(isRedundant?"dotted":"bold")+"\",label=\""+intervalGraphPointers[node1].nodeID+"="+overhang1+","+intervalGraphPointers[node2].nodeID+"="+overhang2+"\",penwidth=\""+PEN_WIDTH+"\",color=\""+Colors.COLOR_DOT_EDGE+"\",fontcolor=\""+Colors.COLOR_DOT_TEXT+"\",fontname=\"Helvetica\",fontsize=\"8\"];\n";
		}
		
		
		/**
		 * @param node false=node1, true=node2.
		 */
		private final String getArrow(int node) {
			if (node==node1) {
				switch (type) {
					case Constants.OVERLAP_PREFIX_PREFIX: return "tee";
					case Constants.OVERLAP_PREFIX_SUFFIX: return "tee";
					case Constants.OVERLAP_SUFFIX_PREFIX: return "inv";
					case Constants.OVERLAP_SUFFIX_SUFFIX: return "inv";
				}
			}
			else {
				switch (type) {
					case Constants.OVERLAP_PREFIX_PREFIX: return "tee";
					case Constants.OVERLAP_PREFIX_SUFFIX: return "inv";
					case Constants.OVERLAP_SUFFIX_PREFIX: return "tee";
					case Constants.OVERLAP_SUFFIX_SUFFIX: return "inv";
				}
			}
			return "";
		}
	}
	
	
}