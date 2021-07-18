package de.mpi_cbg.revant.util;

import java.io.IOException;
import java.io.BufferedWriter;
import de.mpi_cbg.revant.util.Math;


public class DAG {

	/**
	 * Stores in $components[i]$ the ID of the connected component of vertex $i$ in the
	 * DAG represented by matrices $inNeighbors$ and $outNeighbors$. Connected component
	 * IDs are between zero and the number of connected components minus one.
	 *
	 * @param stack temporary space of size at least $nVertices$;
	 * @return the number of connected components.
	 */
	public static final int getConnectedComponents(int nVertices, int[][] inNeighbors, int[] nInNeighbors, int[][] outNeighbors, int[] nOutNeighbors, int[] components, int[] stack) {
		int i, j;
		int lastComponent, lastVertex, from, to;
		
		Math.set(components,nVertices-1,-1);
		lastComponent=-1;
		for (i=0; i<nVertices; i++) {
			if (components[i]!=-1) continue;
			lastComponent++;
			components[i]=lastComponent;
			lastVertex=0;
			stack[lastVertex]=i;
			while (lastVertex>=0) {
				from=stack[lastVertex--];
				for (j=0; j<nInNeighbors[from]; j++) {
					to=inNeighbors[from][j];
					if (components[to]!=-1) {
						if (components[to]!=lastComponent) {
							IO.printCriticalErr("getConnectedComponents> ERROR 1");
							System.exit(1);
						}
						else continue;
					}
					components[to]=lastComponent;
					stack[++lastVertex]=to;
				}
				for (j=0; j<nOutNeighbors[from]; j++) {
					to=outNeighbors[from][j];
					if (components[to]!=-1) {
						if (components[to]!=lastComponent) {
							IO.printCriticalErr("getConnectedComponents> ERROR 2");
							System.exit(1);
						}
						else continue;
					}
					components[to]=lastComponent;
					stack[++lastVertex]=to;
				}
			}
		}
		return lastComponent+1;
	}


	/**
	 * Stores in $componentSize$ the number of vertices in each connected component.
	 *
	 * @param components the connected component of each vertex;
	 * @return the largest size of a connected component.
	 */
	public static final int largestComponent(int[] components, int nVertices, int[] componentSize, int nComponents) {
		int i, maxSize;

		Math.set(componentSize,nComponents-1,0);
		for (i=0; i<nVertices; i++) {
			if (components[i]==-1) {
				IO.printCriticalErr("Error: a component is -1");
				System.exit(1);
			}
			componentSize[components[i]]++;
		}
		maxSize=0;
		for (i=0; i<nComponents; i++) {
			if (componentSize[i]>maxSize) maxSize=componentSize[i];
		}
		return maxSize;
	}


	/**
	 * Stores in $longestPaths[i]$ the length of a longest path in connected component $i$
	 * for all connected components.
	 *
	 * @param sorted2original output of procedure $topologicalSort$, in which every
	 * connected component is a consecutive block;
	 * @param components the connected component of each vertex;
	 * @param pathLengths temporary space with at least $nVertices$ elements, used to store
	 * the length of a longest path from a minimal vertex to each vertex in a connected
	 * component;
	 * @return the length of a longest path in any connected component.
	 */
	public static final int longestPathInComponent(int nVertices, int[] sorted2original, int[] nInNeighbors, int[][] inNeighbors, int[] components, Point[] longestPaths, int[] pathLengths) {
		int i, j;
		int vertex, currentComponent, currentMax, max, length, neighborLength, nIn;
		
		Math.set(pathLengths,nVertices-1,0);
		currentComponent=-1; currentMax=0; max=0;
		for (i=0; i<nVertices; i++) {
			vertex=sorted2original[i];
			if (components[vertex]!=currentComponent && currentComponent!=-1) {
				longestPaths[currentComponent].position=currentMax;
				longestPaths[currentComponent].mass=1;
				if (currentMax>max) max=currentMax;
				currentMax=0;
			}
			currentComponent=components[vertex];
			nIn=nInNeighbors[vertex];
			length=0;
			for (j=0; j<nIn; j++) {
				neighborLength=pathLengths[inNeighbors[vertex][j]];
				if (neighborLength>length) length=neighborLength;
			}
			length++;
			pathLengths[vertex]=length;
			if (length>currentMax) currentMax=length;
		}
		if (currentComponent!=-1) {
			longestPaths[currentComponent].position=currentMax;
			longestPaths[currentComponent].mass=1;
			if (currentMax>max) max=currentMax;
		}
		return max;
	}

	
	/**
	 * Stores in $minimalVertices$ the minimal vertices in the partial order of the DAG, 
	 * sorted by connected component.
	 *
	 * @param nVerticesPerComponent temporary space, of size at least $nComponents$, used 
	 * to store the number of minimal vertices in each connected component;
	 * @return the total number of minimal vertices.
	 */
	public static final int getMinimalVertices(int nVertices, int[][] inNeighbors, int[] nInNeighbors, int[] components, int nComponents, int[] minimalVertices, int[] nVerticesPerComponent) {
		int i, lastMinimal;
		
		Math.set(nVerticesPerComponent,nComponents-1,0);
		for (i=0; i<nVertices; i++) {
			if (nInNeighbors[i]==0) nVerticesPerComponent[components[i]]++;
		}
		for (i=1; i<nComponents; i++) nVerticesPerComponent[i]+=nVerticesPerComponent[i-1];
		lastMinimal=nVerticesPerComponent[nComponents-1]-1;
		for (i=nComponents-1; i>=1; i--) nVerticesPerComponent[i]=nVerticesPerComponent[i-1];
		nVerticesPerComponent[0]=0;
		for (i=0; i<nVertices; i++) {
			if (nInNeighbors[i]==0) minimalVertices[nVerticesPerComponent[components[i]]++]=i;
		}
		return lastMinimal+1;
	}


	/**
	 * Stores in $sorted2original$ the topologically sorted list of vertex IDs of the DAG
	 * of $nVertices$ vertices described by the matrices $inNeighbors$ and $outNeighbors$.
	 * $sorted2original$ consists of blocks such that all vertices inside a block belong
	 * to the same connected component. However, the minimal vertices of a connected
	 * component are not necessarily at the beginning of its block. The procedure stores
	 * in $original2sorted[i]$ the position in $sorted2original$ of the vertex with ID 
	 * $i$.
	 *
	 * Remark: $nInNeighbors$ and $minimalVertices$ are modified by the procedure.
	 * In particular, at the end of the procedure, $minimalVertices$ DOES NOT necessarily 
	 * contain the set of all minimal vertices.
	 *
	 * @param components connected component ID of each vertex;
	 * @param minimalVertices temporary space of size at least $nVertices$;
	 * @param lastMinimal if >=0, the procedure assumes that minimal vertices have already
	 * been computed by the caller, and stored in $minimalVertices[0..lastMinimal]$;
	 * @param nVerticesPerComponent temporary space of size at least $nComponents$, used 
	 * to store the number of minimal vertices in each connected component;
	 * @param sorted2original, original2sorted temporary space of size at least 
	 * $nVertices$;
	 * @return 0 if the graph is a DAG, $i+1$ if node $i$ is involved in a directed cycle.
	 */
	public static final int topologicalSort(int nVertices, int[][] inNeighbors, int[] nInNeighbors, int[][] outNeighbors, int[] nOutNeighbors, int[] components, int nComponents, int[] minimalVertices, int lastMinimal, int[] nVerticesPerComponent, int[] sorted2original, int[] original2sorted) {
		int i, j, k;
		int lastOut, component;

		// Finding minimal vertices in the partial order of the DAG, and sorting them by
		// connected component.
		if (lastMinimal<0) {
			lastMinimal=getMinimalVertices(nVertices,inNeighbors,nInNeighbors,components,nComponents,minimalVertices,nVerticesPerComponent)-1;
			if (lastMinimal==-1) return -1;
		}

		// Sorting topologically
		lastOut=-1;
		while (lastMinimal>=0) {
			i=minimalVertices[lastMinimal--];
			sorted2original[++lastOut]=i;
			original2sorted[i]=lastOut;
			for (j=0; j<nOutNeighbors[i]; j++) {
				k=outNeighbors[i][j];
				nInNeighbors[k]--;
				if (nInNeighbors[k]==0) minimalVertices[++lastMinimal]=k;
			}
		}

		// Checking whether the graph has a directed cycle
		for (i=0; i<nVertices; i++) {
			if (nInNeighbors[i]>0) return i+1;
		}
		return 0;
	}
	
	
	/**
	 * Stores in $nPaths[i]$ the number of directed paths from $source$ to $i$, for every 
	 * $i$.
	 *
	 * @param nPaths assumed initialized to all zeros.
	 */
	public static final void nPaths(int source, int[][] inNeighbors, int[] nInNeighbors, int[] sorted2original, int[] original2sorted, int nVertices, int[] nPaths) {
		int i, j, s, v;
		
		nPaths[source]=1;
		s=original2sorted[source];
		for (i=s+1; i<nVertices; i++) {
			v=sorted2original[i];
			for (j=0; j<nInNeighbors[v]; j++) nPaths[v]+=nPaths[inNeighbors[v][j]];
		}
	}
	
	
	public static final void nPaths_clean(int source, int[] sorted2original, int[] original2sorted, int nVertices, int[] nPaths) {
		for (int i=original2sorted[source]; i<nVertices; i++) nPaths[sorted2original[i]]=0;
	}
	
	
	/**
	 * Given a topologically-sorted DAG, the procedure stores in $weights[i][0]$ 
	 * (respectively, $weights[i][1]$) the minimum (respectively, maximum) weight of a 
	 * path from $source$ to vertex $i$, for every $i$. Column 2 (respectively, 3) stores 
	 * the ID of the in-neighbor in a min-weight (respectively, max-weight) path.
	 *
	 * @param inWeights positive weight of every incoming arc to every vertex;
	 * @param weights assumed to be initialized to all zeros.
	 */
	public static final void pathWeights(int source, int[][] inNeighbors, int[] nInNeighbors, int[][] outNeighbors, int[] nOutNeighbors, int[] sorted2original, int[] original2sorted, int nVertices, int[][] inWeights, int[][] weights) {
		int i, j;
		int v, neighbor, wMin, wMax;
		
		
		for (i=original2sorted[source]+1; i<nVertices; i++) {
			v=sorted2original[i];
			for (j=0; j<nInNeighbors[v]; j++) {
				neighbor=inNeighbors[v][j];
				if (neighbor!=source && weights[neighbor][1]==0) continue;  // Unreachable from $source$.
				wMin=weights[neighbor][0]+inWeights[v][j];
				wMax=weights[neighbor][1]+inWeights[v][j];
				if (weights[v][0]==0 || wMin<weights[v][0]) {
					weights[v][0]=wMin;
					weights[v][2]=neighbor;
				}
				if (weights[v][1]==0 || wMax>weights[v][1]) {
					weights[v][1]=wMax;
					weights[v][3]=neighbor;
				}
			}
		}
	}
	
	
	/**
	 * Stores in $out[outRow]$ the nodes of the longest path to $destination$ that is 
	 * encoded in $weights$.
	 *
	 * @param longestOrShortest if FALSE, stores a shortest path;
	 * @param weights as defined in $pathWeights()$.
	 */
	public static final void getLongestPath(int destination, int[][] weights, int[][] out, int[] outLength, int outRow, boolean longestOrShortest) {
		final int COLUMN = longestOrShortest?3:2;
		int i, j;
		int length;
		
		i=destination; length=0;
		while (i!=-1) {
			length++;
			i=weights[i][COLUMN];
		}
		if (length>out[outRow].length) out[outRow] = new int[length];
		i=destination; j=length-1;
		while (i!=-1) {
			out[outRow][j--]=i;
			i=weights[i][COLUMN];
		}
		outLength[outRow]=length;
	}
	
	
	/**
	 * Resets to zero all cells of $weights$ altered by $pathWeights()$.
	 */
	public static final void pathWeights_clean(int source, int[] sorted2original, int[] original2sorted, int nVertices, int[][] weights) {
		int i, j;
		
		for (i=original2sorted[source]+1; i<nVertices; i++) {
			j=sorted2original[i];
			weights[j][0]=0; weights[j][1]=0; weights[j][2]=-1; weights[j][3]=-1;
		}
	}
	
	
	/**
	 * Positive arcs are dotted, negative arcs are continuous lines.
	 *
	 * @param componentID if nonnegative, prints just the component with the specified ID.
	 */
	public static final void printDOT(int nVertices, int[][] outNeighbors, int[] nOutNeighbors, int[] components, int componentID, BufferedWriter bw) throws IOException {
		int i, j, vertex;
		String style;
		
		bw.write("digraph G {\n");
		bw.write("node [shape=point];\n");
		for (i=0; i<nVertices; i++) {
			if (componentID!=-1 && components[i]!=componentID) continue;
			for (j=0; j<nOutNeighbors[i]; j++) {
				vertex=outNeighbors[i][j];
				if (vertex>=0) style=" [style=dotted];\n";
				else {
					vertex=-1-vertex;
					style=";\n";
				}
				bw.write(i+" -> "+vertex+style);
			}
		}
		bw.write("}\n");
	}
	
	
	/**
	 * Sets to $-1-i$ the elements $i$ of $inNeighbors$ and $outNeighbors$ that correspond 
	 * to non-redundant arcs by transitivity. This very simple algorithm computes the 
	 * longest path from each vertex $v$ to every one of its out-neighbors $w$, and marks 
	 * $(v,w)$ as non-redundant iff the length of a longest path from $v$ to $w$ is one.
	 *
	 * Remark: the procedure assumes $topologicalSort()$ to have already been executed.
	 *
	 * @param components the connected component of each vertex;
	 * @param isActive temporary space with at least $nVertices$ elements, initialized to 
	 * false;
	 * @param longestPath temporary space with at least $nVertices$ elements, initialized 
	 * to -1.
	 */
	public static final void transitiveReduction(int nVertices, int[][] inNeighbors, int[] nInNeighbors, int[][] outNeighbors, int[] nOutNeighbors, int[] sorted2original, int[] original2sorted, int[] components, boolean[] isActive, int[] longestPath) {
		int i, j;
		int source, s, component, activeMet, vertex, inNeighbor, length, maxLength;
		
		for (source=0; source<nVertices; source++) {
			// Computing longest paths
			component=components[source];
			for (i=0; i<nOutNeighbors[source]; i++) {
				isActive[outNeighbors[source][i]]=true;
			}
			activeMet=0; longestPath[source]=0;
			s=original2sorted[source];
			for (i=s+1; i<nVertices; i++) {
				vertex=sorted2original[i];
				if (components[vertex]!=component) break;
				maxLength=-1;
				for (j=0; j<nInNeighbors[vertex]; j++) {
					inNeighbor=inNeighbors[vertex][j];
					if (inNeighbor<0) inNeighbor=-1-inNeighbor;
					length=longestPath[inNeighbor];
					if (length>maxLength) maxLength=length;					
				}
				if (maxLength==-1) continue;  // Unreachable from $source$
				longestPath[vertex]=maxLength+1;
				if (isActive[vertex]) activeMet++;
				if (activeMet==nOutNeighbors[source]) break;
			}
			
			// Marking non-redundant arcs
			for (i=0; i<nOutNeighbors[source]; i++) {
				vertex=outNeighbors[source][i];
				if (longestPath[vertex]>1) continue;
				if (IO.CONSISTENCY_CHECKS && longestPath[vertex]<1) {
					System.err.println("DAG.transitiveReduction> ERROR: wrong longest path.");
					System.exit(1);
				}
				outNeighbors[source][i]=-1-outNeighbors[source][i];
				for (j=0; j<nInNeighbors[vertex]; j++) {
					if (inNeighbors[vertex][j]==source) {
						inNeighbors[vertex][j]=-1-inNeighbors[vertex][j];
						break;
					}
				}
			}
			
			// Cleaning up for next iteration
			activeMet=0;
			longestPath[source]=-1;
			for (i=s+1; i<nVertices; i++) {
				vertex=sorted2original[i];
				if (components[vertex]!=component) break;
				longestPath[vertex]=-1;
				if (isActive[vertex]) activeMet++;
				if (activeMet==nOutNeighbors[source]) break;
			}
			for (i=0; i<nOutNeighbors[source]; i++) {
				vertex=outNeighbors[source][i];
				if (vertex<0) vertex=-1-vertex;
				isActive[vertex]=false;
			}
		}
	}
	
	
	public static final int nArcs(int nVertices, int[][] outNeighbors, int[] nOutNeighbors) {
		int out = 0;
		for (int i=0; i<nVertices; i++) out+=nOutNeighbors[i];
		return out;
	}
	
	
	/**
	 * Disregards positive entries in $outNeighbors$.
	 */
	public static final int nNegativeArcs(int nVertices, int[][] outNeighbors, int[] nOutNeighbors) {
		int i, j, out;
		
		out=0;
		for (i=0; i<nVertices; i++) {
			for (j=0; j<nOutNeighbors[i]; j++) {
				if (outNeighbors[i][j]<0) out++;
			}
		}
		return out;
	}
	

}