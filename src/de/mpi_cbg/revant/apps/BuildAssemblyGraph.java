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
		
		boolean keep;
		int i, j, k;
		int type, nAlignments, idGenerator, top, node, neighbor, nComponents, size, nextFullyUnique;
		double errorRate;
		String str1, str2, str3;
		RepeatAlphabet.Character tmpChar;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		boolean[] containsUnique;
		int[] component, componentSize, stack;
		BufferedWriter[] bws;
		
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
			addEdge(Alignments.readA-1,Alignments.readB-1,keep);
			addEdge(Alignments.readB-1,Alignments.readA-1,keep);
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
			if (nAlignments%10000==0) System.err.println("Processed "+nAlignments+" alignments");
		}
		br1.close(); br2.close(); br3.close();
		System.err.println("Removing duplicated edges...");
		for (i=0; i<N_READS; i++) {
			if (lastNeighbor[i]<=0) continue;
			Arrays.sort(neighbors[i],0,lastNeighbor[i]+1);
			k=0;
			for (j=1; j<=lastNeighbor[i]; j++) {
				if (neighbors[i][j]!=neighbors[i][k]) neighbors[i][++k]=neighbors[i][j];
			}
			lastNeighbor[i]=k;
		}
		
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
				for (j=0; j<=lastNeighbor[node]; j++) {
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
				for (j=0; j<=lastNeighbor[i]; j++) {
					neighbor=neighbors[i][j]>=0?neighbors[i][j]:-1-neighbors[i][j];
					if (i<neighbor) bws[k].write(i+" -- "+neighbor+";\n");
				}
			}
		}
		for (i=0; i<nComponents; i++) { bws[i].write("}"); bws[i].close(); }
		br1.close();
		
		// Outputting filtered connected components
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
				for (j=0; j<=lastNeighbor[node]; j++) {
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
				for (j=0; j<=lastNeighbor[i]; j++) {
					if (neighbors[i][j]>=0 && i<neighbors[i][j]) bws[k].write(i+" -- "+neighbors[i][j]+";\n");
				}
			}
		}
		for (i=0; i<nComponents; i++) { bws[i].write("}"); bws[i].close(); }
		br1.close();
	}
	
	
	/**
	 * Adds $to$ to $from$.
	 */
	private static final void addEdge(int from, int to, boolean keep) {
		final int CAPACITY = 2;  // Arbitrary
		
		lastNeighbor[from]++;
		if (neighbors[from]==null || neighbors[from].length==0) neighbors[from] = new int[CAPACITY];
		else if (lastNeighbor[from]==neighbors[from].length) {
			int[] newArray = new int[neighbors[from].length<<1];
			System.arraycopy(neighbors[from],0,newArray,0,neighbors[from].length);
			neighbors[from]=newArray;
		}
		neighbors[from][lastNeighbor[from]]=keep?to:-1-to;
	}
	
}