package de.mpi_cbg.revant.apps;

import java.util.Arrays;
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
	 * an edge.
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
		
		boolean keep, hasPrefix, hasSuffix, hasPrefixPrime, hasSuffixPrime, hasLoop;
		int i, j, k;
		int type, nAlignments, idGenerator, top, node, neighbor, nNeighbors, nComponents, size, nextFullyUnique;
        int prefixOrSuffixA, prefixOrSuffixB, max, nSingletons, nTips, nLoops;
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
			if ((ALIGNMENT_TYPE==2 && type!=0) || (ALIGNMENT_TYPE==1 && type==2)) {
				str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
				continue;
			}
			keep=str2.equalsIgnoreCase("1")&&str3.equalsIgnoreCase("1");
            if (type==0) {
                prefixOrSuffixA=Alignments.startA<=IDENTITY_THRESHOLD?0:1;
                prefixOrSuffixB=Alignments.startB<=IDENTITY_THRESHOLD?0:1;
            }
            else { prefixOrSuffixA=2; prefixOrSuffixB=2; }
            
            
            
            
if (Alignments.readA==1110) {
    String type2 = "";
    if (Alignments.startA<=100) type2="PREFIX";
    else if (Alignments.endA>=Reads.getReadLength(Alignments.readA-1)-100) type2="SUFFIX";
    String keep1 = str2.equalsIgnoreCase("1")?"KEEP-N":"-";
    String keep2 = str3.equalsIgnoreCase("1")?"KEEP-T":"-";
    System.err.println((keep?"KEEP ":"DELETE ")+" "+keep1+" "+keep2+" "+type2+" "+str1);
}
            
            
            
			addEdge(Alignments.readA-1,Alignments.readB-1,keep,prefixOrSuffixA);
			addEdge(Alignments.readB-1,Alignments.readA-1,keep,prefixOrSuffixB);
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
			if (nAlignments%10000==0) System.err.println("Processed "+nAlignments+" alignments");
		}
		br1.close(); br2.close(); br3.close();
        
        // Removing duplicated edges
		System.err.println("Removing duplicated edges...");
        max=0;
        for (i=0; i<N_READS; i++) {  
            if (lastNeighbor[i]>max) max=lastNeighbor[i];
        }
        max=(max+1)>>1;
        edges = new Edge[max];
        for (i=0; i<max; i++) edges[i] = new Edge();
		for (i=0; i<N_READS; i++) {
			if (lastNeighbor[i]<=1) continue;
            nNeighbors=(lastNeighbor[i]+1)>>1;
            for (j=0; j<=lastNeighbor[i]; j+=2) {
                k=j>>1;
                edges[k].neighbor=neighbors[i][j];
                edges[k].prefixOrSuffix=neighbors[i][j+1];
            }
			Arrays.sort(edges,0,nNeighbors);
			k=0;
			for (j=1; j<nNeighbors; j++) {
				if (edges[j].neighbor!=edges[k].neighbor || edges[j].prefixOrSuffix!=edges[k].prefixOrSuffix) {
                    k++; tmpEdge=edges[k]; edges[k]=edges[j]; edges[j]=tmpEdge;
                }
			}
            lastNeighbor[i]=-1;
            for (j=0; j<=k; j++) {
                neighbors[i][++lastNeighbor[i]]=edges[j].neighbor;
                neighbors[i][++lastNeighbor[i]]=edges[j].prefixOrSuffix;
            }
		}
        
        // Computing number of nodes with only-prefix or only-suffix alignments
        nSingletons=0; nTips=0; nLoops=0;
        for (i=0; i<N_READS; i++) {
            if (lastNeighbor[i]==-1) { nSingletons++; continue; }
            hasPrefix=false; hasSuffix=false; hasLoop=false;
            hasPrefixPrime=false; hasSuffixPrime=false;
            for (j=0; j<=lastNeighbor[i]; j+=2) {
                if (neighbors[i][j+1]==0) hasPrefix=true;
                else if (neighbors[i][j+1]==1) hasSuffix=true;
                if (j>0 && neighbors[i][j]!=neighbors[i][j-2]) {
                    if (hasPrefixPrime && hasSuffixPrime) hasLoop=true;
                    hasPrefixPrime=false; hasSuffixPrime=false;
                }
                if (neighbors[i][j+1]==0) hasPrefixPrime=true;
                else if (neighbors[i][j+1]==1) hasSuffixPrime=true;
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
		
		// Printing input connected components
		System.err.println("Computing input connected components...");
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
		System.err.println("Printing "+nComponents+" components of size at least "+MIN_COMPONENT_SIZE);
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
		
		// Printing filtered connected components
		System.err.println("Computing filtered connected components...");
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
        Arrays.sort(tmpArray);
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
     * @param prefixOrSuffixFrom the alignment uses a prefix (0), suffix (1), or neither 
     * (2), of $from$.
	 */
	private static final void addEdge(int from, int to, boolean keep, int prefixOrSuffixFrom) {
		final int CAPACITY = 4;  // Arbitrary
		
		if (neighbors[from]==null || neighbors[from].length==0) neighbors[from] = new int[CAPACITY];
		else if (lastNeighbor[from]+2>=neighbors[from].length) {
			int[] newArray = new int[neighbors[from].length<<1];
			System.arraycopy(neighbors[from],0,newArray,0,lastNeighbor[from]+1);
			neighbors[from]=newArray;
		}
		neighbors[from][++lastNeighbor[from]]=keep?to:-1-to;
        neighbors[from][++lastNeighbor[from]]=prefixOrSuffixFrom;
	}
	
    
    private static class Edge implements Comparable {
        public int neighbor, prefixOrSuffix;
        
        public Edge() { }
        
        public void set(int n, int p) {
            this.neighbor=n; this.prefixOrSuffix=p;
        }
        
        public int compareTo(Object other) {
            Edge otherEdge = (Edge)other;
            if (neighbor<otherEdge.neighbor) return -1;
            else if (neighbor>otherEdge.neighbor) return 1;
            if (prefixOrSuffix<otherEdge.prefixOrSuffix) return -1;
            else if (prefixOrSuffix>otherEdge.prefixOrSuffix) return 1;
            return 0;
        }
        
        public boolean equals(Object other) {
            Edge otherEdge = (Edge)other;
            return neighbor==otherEdge.neighbor && prefixOrSuffix==otherEdge.prefixOrSuffix;
        }
    }
    
}