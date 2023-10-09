package de.mpi_cbg.revant.apps;

import java.util.Arrays;
import java.util.PriorityQueue;
import java.io.*;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;

/**
 * Prints a simplified assembly graph from the alignments before and after filtering.
 * This should help checking whether the filtering induces a visible decrease in
 * complexity while at the same time preserving linearity.
 */
public class BuildAssemblyGraph {
	/**
	 * Simplistic encoding of an assembly graph. Every read is a node, every alignment is
	 * an edge. See $addEdge()$ for details.
	 */
	private static int[][] neighbors;
	private static int[] lastNeighbor;
	
	/**
	 * Remark: the program prints only large components.
	 *
	 * @param args
	 * 2: 0=loose, 1=tight;
	 * 3: 2=use only suffix-prefix alignments; 1=use every alignment, except containment/
	 * identity; 0=use every alignment.
	 */
	public static void main(String[] args) throws IOException {
		final String INPUT_DIR = args[0];
		final int N_READS = Integer.parseInt(args[1]);
		final String FILTERING_MODE = args[2];
		final int ALIGNMENT_TYPE = Integer.parseInt(args[3]);
		final double MAX_ERROR_RATE = Double.parseDouble(args[4]);
		final String OUTPUT_DIR = args[5];
		
		final String ALIGNMENTS_FILE = INPUT_DIR+"/LAshow-reads-reads.txt";
		final String BITVECTOR_FILE = INPUT_DIR+"/LAshow-reads-reads.txt.mode"+FILTERING_MODE+".bitvector";
		final String TANDEM_BITVECTOR_FILE = INPUT_DIR+"/LAshow-reads-reads.txt.tandem.bitvector";
		final String FULLY_UNIQUE_FILE = INPUT_DIR+"/reads-fullyUnique-new.txt";
		final String READS_TRANSLATED_FILE = INPUT_DIR+"/reads-translated-disambiguated.txt";
		final String BOUNDARIES_FILE = INPUT_DIR+"/reads-translated-boundaries-new.txt";
		final String ALPHABET_FILE = INPUT_DIR+"/alphabet-cleaned.txt";
		final String READ_IDS_FILE = INPUT_DIR+"/reads-ids.txt";
		final String READ_LENGTHS_FILE = INPUT_DIR+"/reads-lengths.txt";
			
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_COMPONENT_SIZE = 10;  // Arbitrary
        final int MAX_DISTANCE = 10;  // Arbitrary, in edges.
		
		boolean keep, hasPrefix, hasSuffix, hasPrefixPrime, hasSuffixPrime, hasLoop;
		int i, j, k;
		int type, nAlignments, idGenerator, top, node, neighbor, nNeighbors, nComponents, size, nextFullyUnique;
        int prefixOrSuffixA, prefixOrSuffixB, maxNeighbors, nSingletons, nTips, nLoops;
		double errorRate;
		String str1, str2, str3;
		RepeatAlphabet.Character tmpChar;
        Edge tmpEdge;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		boolean[] containsUnique;
		int[] component, componentSize, stack;
		BufferedWriter[] bws;
        Edge[] edges;
		
		// Building the graph
		System.err.println("Building the graph...");
		tmpChar = new RepeatAlphabet.Character();
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		Reads.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		neighbors = new int[N_READS][0];
		lastNeighbor = new int[N_READS];
		Math.set(lastNeighbor,N_READS-1,-1);
		br1 = new BufferedReader(new FileReader(ALIGNMENTS_FILE));
		br2 = new BufferedReader(new FileReader(BITVECTOR_FILE));
		br3 = new BufferedReader(new FileReader(TANDEM_BITVECTOR_FILE));
		str1=br1.readLine(); str1=br1.readLine();  // Skipping header
		str1=br1.readLine(); 
		str2=br2.readLine(); str3=br3.readLine(); nAlignments=0;
		while (str1!=null) {
			nAlignments++;
			Alignments.readAlignmentFile(str1);
			errorRate=((double)(Alignments.diffs<<1))/(Alignments.endA-Alignments.startA+1+Alignments.endB-Alignments.startB+1);
			if (errorRate>MAX_ERROR_RATE) {
				str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
				continue;
			}
			type=Alignments.readAlignmentFile_getType(IDENTITY_THRESHOLD);
			if ((ALIGNMENT_TYPE==2 && type>3) || (ALIGNMENT_TYPE==1 && type>4)) {
				str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
				continue;
			}
			keep=str2.equalsIgnoreCase("1")&&str3.equalsIgnoreCase("1");
			addEdge(Alignments.readA-1,Alignments.readB-1,keep,type);
			addEdge(Alignments.readB-1,Alignments.readA-1,keep,Alignments.transposeType[type]);
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
			if (nAlignments%10000==0) System.err.println("Processed "+nAlignments+" alignments");
		}
		br1.close(); br2.close(); br3.close();
        
        // Removing duplicated edges
		System.err.println("Removing duplicated edges...");
        maxNeighbors=0;
        for (i=0; i<N_READS; i++) {  
            if (lastNeighbor[i]>maxNeighbors) maxNeighbors=lastNeighbor[i];
        }
        maxNeighbors=(maxNeighbors+1)>>1;
        edges = new Edge[maxNeighbors];
        for (i=0; i<maxNeighbors; i++) edges[i] = new Edge();
		for (i=0; i<N_READS; i++) {
			if (lastNeighbor[i]<=1) continue;
            nNeighbors=(lastNeighbor[i]+1)>>1;
            for (j=0; j<=lastNeighbor[i]; j+=2) {
                k=j>>1;
                edges[k].neighbor=neighbors[i][j];
                edges[k].type=neighbors[i][j+1];
            }
			Arrays.parallelSort(edges,0,nNeighbors);
			k=0;
			for (j=1; j<nNeighbors; j++) {
				if (edges[j].neighbor!=edges[k].neighbor || edges[j].type!=edges[k].type) {
                    k++; tmpEdge=edges[k]; edges[k]=edges[j]; edges[j]=tmpEdge;
                }
			}
            lastNeighbor[i]=-1;
            for (j=0; j<=k; j++) {
                neighbors[i][++lastNeighbor[i]]=edges[j].neighbor;
                neighbors[i][++lastNeighbor[i]]=edges[j].type;
            }
		}
        
        // Removing spurious edges
        getCutVertices(N_READS);
        simplify(MAX_DISTANCE,maxNeighbors,N_READS);
        
        // Computing number of nodes with only-prefix or only-suffix alignments
        nSingletons=0; nTips=0; nLoops=0;
        for (i=0; i<N_READS; i++) {
            if (lastNeighbor[i]==-1) { nSingletons++; continue; }
            hasPrefix=false; hasSuffix=false; hasLoop=false;
            hasPrefixPrime=false; hasSuffixPrime=false;
            for (j=0; j<=lastNeighbor[i]; j+=2) {
                if (j>0 && neighbors[i][j]!=neighbors[i][j-2]) {
                    if (hasPrefixPrime && hasSuffixPrime) hasLoop=true;
                    hasPrefixPrime=false; hasSuffixPrime=false;
                }
                if (neighbors[i][j+1]<=1) { hasPrefix=true; hasPrefixPrime=true; }
                else if (neighbors[i][j+1]>=2 && neighbors[i][j+1]<=3) { hasSuffix=true; hasSuffixPrime=true; }
            }
            if (hasPrefixPrime && hasSuffixPrime) hasLoop=true;
            if (hasPrefix!=hasSuffix) nTips++;
            if (hasLoop) nLoops++;
        }
        System.err.println("Reads with no alignment: "+nSingletons+" ("+((100.0*nSingletons)/N_READS)+"%)");
        System.err.println("Reads with only prefix or only suffix alignments: "+nTips+" ("+((100.0*nTips)/N_READS)+"%)");
        System.err.println("Reads with suffix and prefix alignments to the same other read: "+nLoops+" ("+((100.0*nLoops)/N_READS)+"%)");
		
		// Building $containsUnique$.		
		RepeatAlphabet.deserializeAlphabet(ALPHABET_FILE,2);
		containsUnique = new boolean[N_READS];
		Arrays.fill(containsUnique,false);
		br1 = new BufferedReader(new FileReader(READS_TRANSLATED_FILE));
		br2 = new BufferedReader(new FileReader(BOUNDARIES_FILE));
		str1=br1.readLine(); str2=br2.readLine(); j=0;
		while (str1!=null) {
			if (str1.length()!=0) {
				RepeatAlphabet.loadBoundaries(str2);
				containsUnique[j]=RepeatAlphabet.containsUnique(str1,RepeatAlphabet.boundaries,Reads.getReadLength(j),tmpChar);
			}
			j++; str1=br1.readLine(); str2=br2.readLine();
		}
		br1.close(); br2.close();
		
		System.err.println("Computing weakly-connected components of the input...");
		component = new int[N_READS];
		Math.set(component,N_READS-1,-1);
		stack = new int[N_READS]; idGenerator=-1;
		for (i=0; i<N_READS; i++) {
			if (component[i]!=-1) continue;
			idGenerator++;
			component[i]=idGenerator;
			top=0; stack[0]=i; size=1;
			while (top>=0) {
				node=stack[top--];
				for (j=0; j<=lastNeighbor[node]; j+=2) {
					neighbor=neighbors[node][j];
					if (neighbor<0) neighbor=-1-neighbor;
					if (component[neighbor]!=-1) continue;
					component[neighbor]=idGenerator;
					stack[++top]=neighbor;
					size++;
				}
			}
		}
		nComponents=idGenerator+1;
		System.err.println("DONE. "+nComponents+" components.");
		componentSize = new int[nComponents];
		Math.set(componentSize,nComponents-1,0);
		for (i=0; i<N_READS; i++) componentSize[component[i]]++;
        printHistogram(componentSize,nComponents,N_READS);
		j=-1;
		for (i=0; i<nComponents; i++) {
			if (componentSize[i]>=MIN_COMPONENT_SIZE) componentSize[++j]=i;
		}
		nComponents=j+1;
		System.err.println("Printing "+nComponents+" weakly-connected components of size at least "+MIN_COMPONENT_SIZE);
		bws = new BufferedWriter[nComponents];
		for (i=0; i<nComponents; i++) {
			bws[i] = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/all-component-"+componentSize[i]+".dot"));
			bws[i].write("graph G {\n");
		}
		br1 = new BufferedReader(new FileReader(FULLY_UNIQUE_FILE));
		str1=br1.readLine();
		nextFullyUnique=str1!=null?Integer.parseInt(str1):Math.POSITIVE_INFINITY;
		for (i=0; i<N_READS; i++) {
			k=Arrays.binarySearch(componentSize,0,nComponents,component[i]);
			if (k>=0) {
				while (nextFullyUnique<i) {
					str1=br1.readLine();
					nextFullyUnique=str1!=null?Integer.parseInt(str1):Math.POSITIVE_INFINITY;
				}
				if (i==nextFullyUnique) {
					bws[k].write(i+" [uniqueStatus=\"2\"];\n");
					str1=br1.readLine();
					nextFullyUnique=str1!=null?Integer.parseInt(str1):Math.POSITIVE_INFINITY;
				}
				else bws[k].write(i+" [uniqueStatus=\""+(containsUnique[i]?1:0)+"\"];\n");
				for (j=0; j<=lastNeighbor[i]; j+=2) {
					neighbor=neighbors[i][j]>=0?neighbors[i][j]:-1-neighbors[i][j];
					if (i<neighbor) bws[k].write(i+" -- "+neighbor+";\n");
				}
			}
		}
		for (i=0; i<nComponents; i++) { bws[i].write("}"); bws[i].close(); }
		br1.close();
		
		System.err.println("Computing weakly-connected components after filters...");
		Math.set(component,N_READS-1,-1);
		idGenerator=-1;
		for (i=0; i<N_READS; i++) {
			if (component[i]!=-1) continue;
			idGenerator++;
			component[i]=idGenerator;
			top=0; stack[0]=i; size=1;
			while (top>=0) {
				node=stack[top--];
				for (j=0; j<=lastNeighbor[node]; j+=2) {
					neighbor=neighbors[node][j];
					if (neighbor<0 || component[neighbor]!=-1) continue;
					component[neighbor]=idGenerator;
					stack[++top]=neighbor;
					size++;
				}
			}
		}
		nComponents=idGenerator+1;
		System.err.println("DONE. "+nComponents+" components.");
		componentSize = new int[nComponents];
		Math.set(componentSize,nComponents-1,0);
		for (i=0; i<N_READS; i++) componentSize[component[i]]++;
        printHistogram(componentSize,nComponents,N_READS);
		j=-1;
		for (i=0; i<nComponents; i++) {
			if (componentSize[i]>=MIN_COMPONENT_SIZE) componentSize[++j]=i;
		}
		nComponents=j+1;
		System.err.println("Printing "+nComponents+" components of size at least "+MIN_COMPONENT_SIZE);
		if (nComponents>bws.length) bws = new BufferedWriter[nComponents];
		for (i=0; i<nComponents; i++) {
			bws[i] = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/filtered-component-"+componentSize[i]+".dot"));
			bws[i].write("graph G {\n");
		}
		br1 = new BufferedReader(new FileReader(FULLY_UNIQUE_FILE));
		str1=br1.readLine();
		nextFullyUnique=str1!=null?Integer.parseInt(str1):Math.POSITIVE_INFINITY;
		for (i=0; i<N_READS; i++) {
			k=Arrays.binarySearch(componentSize,0,nComponents,component[i]);
			if (k>=0) {
				while (nextFullyUnique<i) {
					str1=br1.readLine();
					nextFullyUnique=str1!=null?Integer.parseInt(str1):Math.POSITIVE_INFINITY;
				}
				if (i==nextFullyUnique) {
					bws[k].write(i+" [uniqueStatus=\"2\"];\n");
					str1=br1.readLine();
					nextFullyUnique=str1!=null?Integer.parseInt(str1):Math.POSITIVE_INFINITY;
				}
				else bws[k].write(i+" [uniqueStatus=\""+(containsUnique[i]?1:0)+"\"];\n");
				for (j=0; j<=lastNeighbor[i]; j+=2) {
					if (neighbors[i][j]>=0 && i<neighbors[i][j]) bws[k].write(i+" -- "+neighbors[i][j]+";\n");
				}
			}
		}
		for (i=0; i<nComponents; i++) { bws[i].write("}"); bws[i].close(); }
		br1.close();
	}
    
    
    private static final void printHistogram(int[] componentSize, int nComponents, int nNodes) {
        int i;
        double count;
        int[] tmpArray;
        
        System.err.println("Number of nodes in components of size >=:");
        tmpArray = new int[nComponents];
        System.arraycopy(componentSize,0,tmpArray,0,nComponents);
        Arrays.parallelSort(tmpArray);
        count=tmpArray[nComponents-1];
        for (i=nComponents-2; i>=0; i--) {
            if (tmpArray[i]!=tmpArray[i+1]) System.err.println(tmpArray[i+1]+","+(count/nNodes));
            count+=tmpArray[i];
        }
        System.err.println(tmpArray[0]+","+(count/nNodes));
    }
    
	
	/**
	 * Adds $to$ to $from$.
     * 
     * @param alignmentType as in $Alignments.readAlignmentFile_getType()$, expressed as 
     * if $from$ were readA and $to$ were readB.
	 */
	private static final void addEdge(int from, int to, boolean keep, int alignmentType) {
		final int CAPACITY = 4;  // Arbitrary
		
		if (neighbors[from]==null || neighbors[from].length==0) neighbors[from] = new int[CAPACITY];
		else if (lastNeighbor[from]+2>=neighbors[from].length) {
			int[] newArray = new int[neighbors[from].length<<1];
			System.arraycopy(neighbors[from],0,newArray,0,lastNeighbor[from]+1);
			neighbors[from]=newArray;
		}
		neighbors[from][++lastNeighbor[from]]=keep?to:-1-to;
        neighbors[from][++lastNeighbor[from]]=alignmentType;
	}
	
    
    private static class Edge implements Comparable {
        public int neighbor, type;
        
        public Edge() { }
        
        public void set(int n, int t) {
            this.neighbor=n; this.type=t;
        }
        
        public int compareTo(Object other) {
            Edge otherEdge = (Edge)other;
            if (neighbor<otherEdge.neighbor) return -1;
            else if (neighbor>otherEdge.neighbor) return 1;
            if (type<otherEdge.type) return -1;
            else if (type>otherEdge.type) return 1;
            return 0;
        }
        
        public boolean equals(Object other) {
            Edge otherEdge = (Edge)other;
            return neighbor==otherEdge.neighbor && type==otherEdge.type;
        }
    }
    
    
    
    
    
    
    
    
    // --------------------------- SIMPLIFICATION PROCEDURES -----------------------------
    
    /**
     * Sorted list of all the cutvertices of the undirected graph induced by the edges
     * that are kept after filtering.
     */
    private static int[] cutVertices;
    private static int lastCutVertex;
    
    /**
     * Reused temporary space, of size equal to the total number of nodes in the graph, to
     * store the shortest path length.
     */
    private static int[] distance;
    
    /**
     * Reused temporary space, at least as big as the max number of neighbors of a node.
     */
    private static int[] tmpArray, component, componentSize, stack, lastNeighborPrime;
    private static int[][] neighborsPrime;
    private static PriorityQueue<Node> tmpQueue;
    
    /**
     * Reused temporary space
     */
    private static Node[] nodePool;
    private static int lastNodeInPool;
    private static int[] edgesToRemove;
    private static int edgesToRemove_last;
    
    
    /**
     * Stores in global variable $cutVertices$ the sorted list of nodes of the graph whose
     * removal increases the number of connected components. Uses the Hopcroft-Tarjan 
     * algorithm, which is sequential and considers the graph as undirected.
     */
    private static final void getCutVertices(int nReads) {
        int i, p;
        int top, nodeID, lastChild, reachableTime, parentPointer, neighbor, pushedChildren, isCutvertex;
        int[] stack, times;
        
        times = new int[nReads];
        Math.set(times,nReads-1,0);
        cutVertices = new int[100];  // Arbitrary
        lastCutVertex=-1;
        stack = new int[500];  // Arbitrary, multiple of 5.
        for (i=0; i<nReads; i++) {
            if (times[i]!=0 || lastNeighbor[i]<=1) continue;
            times[i]=1; pushedChildren=0;
            top=4; stack[0]=i; stack[1]=-2; stack[2]=times[i]; stack[3]=0; stack[4]=-1;
            while (top>=0) {
                p=top;
                nodeID=stack[top-4];
                lastChild=stack[top-3];
                reachableTime=stack[top-2];
                isCutvertex=stack[top-1];
                parentPointer=stack[top];
                lastChild+=2;
                if (lastChild<=lastNeighbor[nodeID]) {
                    neighbor=neighbors[nodeID][lastChild];
                    if (neighbor<0 || (parentPointer!=-1 && neighbor==stack[parentPointer-4])) {
                        stack[p-3]=lastChild;
                        continue;
                    }
                    if (times[neighbor]!=0) reachableTime=Math.min(reachableTime,times[neighbor]);
                    else {
                        times[neighbor]=times[nodeID]+1;
                        if (nodeID==i) pushedChildren++;
                        if (top+5>=stack.length) {
                            int[] newArray = new int[stack.length<<1];
                            System.arraycopy(stack,0,newArray,0,stack.length);
                            stack=newArray;
                        }
                        stack[++top]=neighbor;
                        stack[++top]=-2;
                        stack[++top]=times[neighbor];
                        stack[++top]=0;
                        stack[++top]=p;
                    }
                    stack[p-3]=lastChild; stack[p-2]=reachableTime;
                }
                else {
                    if (isCutvertex==1) {
                        lastCutVertex++;
                        if (lastCutVertex==cutVertices.length) {
                            int[] newArray = new int[cutVertices.length<<1];
                            System.arraycopy(cutVertices,0,newArray,0,cutVertices.length);
                            cutVertices=newArray;
                        }
                        cutVertices[lastCutVertex]=nodeID;
                    }
                    if (parentPointer!=-1 && parentPointer!=4) {
                        if (reachableTime>=times[stack[parentPointer-4]]) stack[parentPointer-1]=1;
                        stack[parentPointer-2]=Math.min(stack[parentPointer-2],reachableTime);
                    }
                    top-=5;
                }
            }
            if (pushedChildren>1) {
                lastCutVertex++;
                if (lastCutVertex==cutVertices.length) {
                    int[] newArray = new int[cutVertices.length<<1];
                    System.arraycopy(cutVertices,0,newArray,0,cutVertices.length);
                    cutVertices=newArray;
                }
                cutVertices[lastCutVertex]=i;
            }
        }
        if (lastCutVertex>0) Arrays.parallelSort(cutVertices,0,lastCutVertex+1);
        times=null; stack=null;
        System.err.print("Found "+(lastCutVertex+1)+" cut vertices: ");
        for (i=0; i<=lastCutVertex; i++) System.err.print(cutVertices[i]+" ");
        System.err.println();
    }
    
    
    /**
     * Removes spurious edges that connect distinct regions of the genome (see 
     * $isNodeAnomalous()$ for details).
     *
     * Remark: the procedure assumes $neighbors$ to be sorted, and $getCutVertices()$ to
     * have already been executed.
     *
     * @param maxDistance see $isNodeAnomalous()$;
     * @param maxNeighbors upper bound on the number of neighbors of a node.
     */
    private static final void simplify(int maxDistance, int maxNeighbors, int nReads) {
        int i, j;
        int nAnomalousNodes, nPairs, first;
        int[] anomalousNodes;
        Pair[] pairs;
        
        // Allocating reused space
        distance = new int[nReads];
        Math.set(distance,nReads-1,Math.POSITIVE_INFINITY);
        nodePool = new Node[maxNeighbors];  // Arbitrary
        for (i=0; i<nodePool.length; i++) nodePool[i] = new Node();
        lastNodeInPool=-1;
        edgesToRemove = new int[20];  // Arbitrary
        tmpQueue = new PriorityQueue<Node>();
        component = new int[maxNeighbors];
        componentSize = new int[maxNeighbors];
        stack = new int[maxNeighbors];
        tmpArray = new int[maxNeighbors];
        neighborsPrime = new int[maxNeighbors][2];  // Arbitrary
        lastNeighborPrime = new int[maxNeighbors];
        anomalousNodes = new int[100];  // Arbitrary
        
        // Collecting edges to remove
        nAnomalousNodes=0; edgesToRemove_last=-1;
        j=0;
        for (i=0; i<nReads; i++) {
            if (j>lastCutVertex || i<cutVertices[j]) {
System.err.println("isNodeAnomalous? "+i);
                if (isNodeAnomalous(i,maxDistance)) {
                    nAnomalousNodes++;
                    if (nAnomalousNodes>anomalousNodes.length) {
                        int[] newArray = new int[nAnomalousNodes<<1];
                        System.arraycopy(anomalousNodes,0,newArray,0,anomalousNodes.length);
                        anomalousNodes=newArray;
                    }
                    anomalousNodes[nAnomalousNodes-1]=i;
                }
            }
            else if (i==cutVertices[j]) j++;
        }
        System.err.print("Found "+nAnomalousNodes+" nodes with spurious edges: ");
        for (i=0; i<nAnomalousNodes; i++) System.err.print(anomalousNodes[i]+" ");
        System.err.println();
        
        // Deallocating reused space
        distance=null;
        for (i=0; i<nodePool.length; i++) nodePool[i]=null;
        nodePool=null;
        tmpQueue.clear(); tmpQueue=null;
        component=null; componentSize=null; stack=null; tmpArray=null;
        for (i=0; i<neighborsPrime.length; i++) neighborsPrime[i]=null;
        neighborsPrime=null;
        
        // Removing edges
        if (edgesToRemove_last==-1) return;
        pairs = new Pair[edgesToRemove_last];
        for (i=0; i<edgesToRemove_last; i+=2) pairs[i>>1] = new Pair(edgesToRemove[i],edgesToRemove[i+1]);
        edgesToRemove=null; nPairs=(edgesToRemove_last+1)>>1;
        if (nPairs>0) {
            Arrays.parallelSort(pairs,0,nPairs);
            j=0;
            for (i=1; i<nPairs; i++) {
                if (!pairs[i].equals(pairs[j])) pairs[++j]=pairs[i];
            }
            nPairs=j+1;
        }
        first=0;
        for (i=1; i<nPairs; i++) {
            if (pairs[i].from!=pairs[first].from) {
                simplify_impl(pairs,first,i-1);
                first=i;
            }
        }
        simplify_impl(pairs,first,nPairs-1);
        System.err.println("Removed "+(nPairs>>1)+" spurious edges");
    }
    
    
    /**
     * Removes from $neighbors$ all the edges in $pairs[start..end]$. 
     * Remark: $neighbors$ is assumed to be sorted.
     */
    private static final void simplify_impl(Pair[] pairs, int start, int end) {
        int i, p, q;
        final int row = pairs[start].from;
        final int last = lastNeighbor[row];
        
        p=start; q=0;
        while (p<=end && q<last) {
            if (neighbors[row][q]>pairs[p].to) { p++; continue; }
            else if (neighbors[row][q]<pairs[p].to) { q+=2; continue; }
            neighbors[row][q]=Math.POSITIVE_INFINITY;
            neighbors[row][q+1]=Math.POSITIVE_INFINITY;
            p++; q+=2;
        }
        i=-1;
        for (q=0; q<=last; q++) {
            if (neighbors[row][q]!=Math.POSITIVE_INFINITY) neighbors[row][++i]=neighbors[row][q];
        }
        lastNeighbor[row]=i;
    }
    
    
    /**
     * In a perfect assembly graph, removing a node would keep its neighbors connected.
     * The procedure returns true iff, after removing $nodeID$ from the graph, its 
     * neighbors are not all mutually reachable via undirected paths through edges that 
     * are kept after filtering. In this case the procedure adds to $edgesToRemove$ every 
     * edge (of any overlap type) from $nodeID$ to a neighbor that does not belong to the
     * largest cluster of reachable neighbors.
     *
     * @param maxDistance two neighbors are considered unreachable if the shortest
     * undirected path between them is longer than this.
     */
    private static final boolean isNodeAnomalous(int nodeID, int maxDistance) {
        int i, j;
        int last, top, lastComponent, currentNode, currentNeighbor, maxSize, maxComponent;
        
        // Collecting the distinct neighbors of $nodeID$ that are different from $nodeID$.
        j=-1;
        last=lastNeighbor[nodeID];
        for (i=0; i<=last; i+=2) {
            if (neighbors[nodeID][i]>=0) tmpArray[++j]=neighbors[nodeID][i];
        }
        if (j>0) Arrays.parallelSort(tmpArray,0,j+1);
        last=j;        
        j=0;
        while (j<=last && tmpArray[j]==nodeID) j++;
        if (j>last) return false;
        else if (j!=0) {
            for (i=j; i<=last; i++) tmpArray[i-j]=tmpArray[i];
            last-=j;
        }
        j=0;
        for (i=1; i<=last; i++) {
            if (tmpArray[i]!=nodeID && tmpArray[i]!=tmpArray[j]) tmpArray[++j]=tmpArray[i];
        }
        last=j;
        
boolean fabio = nodeID==171 || nodeID==142 || nodeID==621;
if (fabio) {
    System.err.print("isNodeAnomalous> 0  neighbors of node "+nodeID+": ");
    for (int x=0; x<=last; x++) System.err.print(tmpArray[x]+" ");
    System.err.println();
}        
        
        // Removing $nodeID$ from the graph and clustering distinct neighbors by mutual
        // reachability via bidirected walks.
        Math.set(lastNeighborPrime,last,-1);
        for (i=0; i<last; i++) {
            for (j=i+1; j<=last; j++) {
if (fabio) System.err.println("isNodeAnomalous> 1  finding shortest path between "+tmpArray[i]+" and "+tmpArray[j]);                
                if (shortestPath(tmpArray[i],tmpArray[j],nodeID,maxDistance)) {
if (fabio) System.err.println("isNodeAnomalous> 1  shortest path found between "+tmpArray[i]+" and "+tmpArray[j]);
                    lastNeighborPrime[i]++;
                    if (lastNeighborPrime[i]==neighborsPrime[i].length) {
                        int[] newArray = new int[neighborsPrime[i].length<<1];
                        System.arraycopy(neighborsPrime[i],0,newArray,0,neighborsPrime[i].length);
                        neighborsPrime[i]=newArray;
                    }
                    neighborsPrime[i][lastNeighborPrime[i]]=j;
                    lastNeighborPrime[j]++;
                    if (lastNeighborPrime[j]==neighborsPrime[j].length) {
                        int[] newArray = new int[neighborsPrime[j].length<<1];
                        System.arraycopy(neighborsPrime[j],0,newArray,0,neighborsPrime[j].length);
                        neighborsPrime[j]=newArray;
                    }
                    neighborsPrime[j][lastNeighborPrime[j]]=i;
                }
            }
        }
        Math.set(component,last,-1); Math.set(componentSize,last,0);
        j=-1; lastComponent=-1;
        for (i=0; i<=last; i++) {
            if (component[i]!=-1) continue;
            component[i]=++lastComponent;
            stack[0]=i; top=0;
            while (top>=0) {
                currentNode=stack[top--];
                componentSize[lastComponent]++;
                for (j=0; j<=lastNeighborPrime[currentNode]; j++) {
                    currentNeighbor=neighborsPrime[currentNode][j];
                    if (component[currentNeighbor]==-1) {
                        component[currentNeighbor]=lastComponent;
                        stack[++top]=currentNeighbor;
                    }
                }
            }
        }
        if (lastComponent==0) return false;
        
        // Keeping only edges that connect $nodeID$ to the largest component.
        maxSize=0; maxComponent=-1;
        for (i=0; i<=lastComponent; i++) {
            if (componentSize[i]>maxSize) { maxSize=componentSize[i]; maxComponent=i; }
        }
        j=-1;
        for (i=0; i<=last; i++) {
            if (component[i]!=maxComponent) {
                addEdgeToRemove(nodeID,tmpArray[i]);
                addEdgeToRemove(tmpArray[i],nodeID);
            }
        }
        return true;
    }
    
    
    /**
     * Remark: the procedure assumes that global variable $distance$ is initialized to 
     * infinity, and it resets it to this state at the end. The procedure assumes that 
     * global variable $tmpQueue$ is empty, and it resets it to this state at the end.
     *
     * Remark: considering the graph as undirected rather than as bidirected makes the 
     * computation of pairwise reachability much faster in practice, when $maxDistance$ is
     * a number of edges rather than bps.
     *
     * @param excludeNode paths that use this node are not considered;
     * @param maxDistance stop the search when only distances greater than this are left;
     * @return TRUE iff there is an undirected path $from--to$ with length $<=maxDistance$
     * and that uses only edges that are kept after filtering.
     */
    private static final boolean shortestPath(int from, int to, int excludeNode, int maxDistance) {
        boolean out, found;
        int i;
        int nodeID, nodeDistance, last, neighbor, edgeType;
        Node currentNode, neighborNode;
        
boolean fabio = from==621 && to==1119;        
        
        distance[from]=0;
        currentNode=getNodeFromPool();
        currentNode.set(from); tmpQueue.add(currentNode);
        while (!tmpQueue.isEmpty()) {
            currentNode=tmpQueue.poll();
            nodeID=currentNode.id; 
            nodeDistance=distance[nodeID];
            if (nodeDistance>maxDistance) break;
            last=lastNeighbor[nodeID];
            
            
if (fabio) System.err.println("Considering node "+nodeID+" with distance "+nodeDistance+" and "+(last+1)+" neighbors");
            
            found=false;
            for (i=0; i<=last; i+=2) {
                neighbor=neighbors[nodeID][i];
                if (neighbor<0 || neighbor==excludeNode) continue;
                if (nodeDistance+1<distance[neighbor]) {
                    neighborNode=getNodeFromPool();
                    neighborNode.set(neighbor);
                    distance[neighbor]=nodeDistance+1;
                    if (neighbor==to && distance[neighbor]<=maxDistance) { 
if (fabio) System.err.println(to+" is a neighbor!");                        
                        found=true; break; 
                    }
                    tmpQueue.remove(neighborNode); tmpQueue.add(neighborNode);
if (fabio) System.err.println("Added node "+neighbor+" with distance "+distance[neighbor]);
                }
            }
            if (found) break;
        }
if (from==621 && to==1119) System.err.println("shortestPath> distance between "+from+" and "+to+" is "+distance[to]);        
        out=distance[to]<=maxDistance;
        for (i=0; i<=lastNodeInPool; i++) distance[nodePool[i].id]=Math.POSITIVE_INFINITY;
        tmpQueue.clear(); lastNodeInPool=-1;
        return out;
    }

    
    private static final Node getNodeFromPool() {
        int i;
        
        lastNodeInPool++;
        if (lastNodeInPool==nodePool.length) {
            Node[] newArray = new Node[nodePool.length<<1];
            System.arraycopy(nodePool,0,newArray,0,nodePool.length);
            for (i=nodePool.length; i<newArray.length; i++) newArray[i] = new Node();
            nodePool=newArray;
        }
        return nodePool[lastNodeInPool];
    }
    
    
    private static final void addEdgeToRemove(int from, int to) {
        if (edgesToRemove_last+2>=edgesToRemove.length) {
            int[] newArray = new int[edgesToRemove.length<<1];
            System.arraycopy(edgesToRemove,0,newArray,0,edgesToRemove.length);
            edgesToRemove=newArray;
        }
        edgesToRemove[++edgesToRemove_last]=from;
        edgesToRemove[++edgesToRemove_last]=to;
    }
    
    
    private static class Node implements Comparable {
        public int id;
        
        public Node() { this.id=Math.POSITIVE_INFINITY; }
        
        public void set(int i) { this.id=i; }
        
        public int compareTo(Object other) {
            Node otherNode = (Node)other;
            final int d1 = distance[id];
            final int d2 = distance[otherNode.id];
            if (d1<d2) return -1;
            else if (d1>d2) return 1;
            return 0;
        }
        
        public boolean equals(Object other) {
            Node otherNode = (Node)other;
            return id==otherNode.id;
        }
    }
    
    
    private static class Pair implements Comparable {
        public int from, to;
        
        public Pair(int f, int t) {
            this.from=f; this.to=t;
        }
        
        public int compareTo(Object other) {
            Pair otherPair = (Pair)other;
            if (from<otherPair.from) return -1;
            else if (from>otherPair.from) return 1;
            if (to<otherPair.to) return -1;
            else if (to>otherPair.to) return 1;
            return 0;
        }
        
        public boolean equals(Object other) {
            Pair otherPair = (Pair)other;
            return from==otherPair.from && to==otherPair.to;
        }
    }
    
    
}