package de.mpi_cbg.revant.biology;

import java.io.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignments;


public class ClusterRepbaseGraph {
	
	private static int[] componentsInfo;


	/**
	
	   java ClusterRepbaseGraph ../repbase/repbaseTree.txt ../repbase/repbaseConcatenation.ref 1 1 ../lashow -lashow.txt ../coordinates -hub.txt ../componentsInfo.txt 1> hubGraph.dot
	
	 * java ClusterRepbaseGraph ../repbase/repbaseTree.txt ../repbase/repbaseConcatenation.ref 0 1287 ../lashow -lashow.txt ../coordinates -family.txt ../componentsInfo.txt 1> clusterRepbaseGraph.dot 2> ./edgeHistograms.txt
	 */
	public static void main(String[] args) throws IOException {
		final String TREE_FILE = args[0];
		final String FASTA_FILE = args[1];
		final int FIRST_FILE = Integer.parseInt(args[2]);
		final int LAST_FILE = Integer.parseInt(args[3]);
		final String ALIGNMENT_FILES_PATH = args[4];
		final String ALIGNMENT_SUFFIX = args[5];
		final String FRAGMENT_FILES_PATH = args[6];
		final String FRAGMENT_SUFFIX = args[7];
		final String COMPONENTS_INFO_FILE = args[8];
		int i, maxDepth;
		StringBuilder[] levels;
		
		IO.initialize();
		maxDepth=RepBase.loadRepbase(TREE_FILE,FASTA_FILE,false);
		System.out.println("digraph G {");
		System.out.println("rankdir=LR;");
		
		// Printing RepBase nodes on distinct levels for each depth
		levels = new StringBuilder[maxDepth+1];
		for (i=0; i<=maxDepth; i++) levels[i] = new StringBuilder();
		RepBase.root.toDot(levels);
		for (i=0; i<=maxDepth; i++) {
			System.out.println("{ rank=same;");
			System.out.println(levels[i]);
			System.out.println("}");
		}
		RepBase.root.toDot();
		
		// Printing clusters->RepBase arcs
		loadComponentsInfo(COMPONENTS_INFO_FILE);
		printClusterRepbaseGraph(ALIGNMENT_FILES_PATH,FIRST_FILE,LAST_FILE,ALIGNMENT_SUFFIX);
		printRepbaseHeatmaps(FIRST_FILE,LAST_FILE,FRAGMENT_FILES_PATH,FRAGMENT_SUFFIX,ALIGNMENT_FILES_PATH,ALIGNMENT_SUFFIX);
		
		System.out.println("}");
	}
	
	
	private static final void loadComponentsInfo(String path) throws IOException {
		int i, nComponents;
		String str;
		BufferedReader br;
		
		// Counting components
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		nComponents=0;
		str=br.readLine();
		while (str!=null) {
			nComponents++;
			str=br.readLine();
		}
		br.close();
		
		// Loading file
		componentsInfo = new int[nComponents];
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		str=br.readLine();
		i=-1;
		while (str!=null) {
			componentsInfo[++i]=Integer.parseInt(str.trim());
			str=br.readLine();
		}
		br.close();
	}
	

	private static final String type2color(int type) {
		switch(type) {
			case 1: return "orange";
			case 2: return "red";
			case 3: return "yellow";
			case 4: return "green";
			default: return "white";
		}
	}
	
	
	/**
	 * Prints to the standard output a DOT file that connects cluster IDs $[firstFile..
	 * lastFile]$ to RepBase leaves, using the alignment files stored in 
	 * $aligmentFilesPath$ (one file per cluster). Prints to the standard error the
	 * coverage histogram of each edge.
	 *
	 * Remark: the procedure assumes that $RepBase.loadRepbase()$ has already been 
	 * executed.
	 */
	private static final void printClusterRepbaseGraph(String aligmentFilesPath, int firstFile, int lastFile, String suffix) throws IOException {
		final double MAX_PEN_WIDTH = 10.0;
		final int MAX_NEIGHBORS = 10;
		final int MAX_LENGTH = 50000;
		int i, j, k;
		int idGenerator, found, lastNeighbor, descendant;
		double totalCount;
		String str, path;
		BufferedReader alignmentsFile;
		File file;
		RepBase.RepbaseTreeNode neighbor, tmpNode;
		boolean[] isActive;
		RepBase.RepbaseTreeNode[] neighbors;
		int[][] histograms = new int[MAX_NEIGHBORS][MAX_LENGTH];
		
		// Putting all cluster nodes in the same DOT subgraph
		idGenerator=RepBase.nNodes-1;
		System.out.println("{ rank=sink; \n");
		for (i=firstFile; i<=lastFile; i++) {
			path=aligmentFilesPath+"/"+i+suffix;
			file = new File(path);
			if (!file.exists()) continue;
			idGenerator++;
			System.out.println(idGenerator+" [color=\""+type2color(componentsInfo[i])+"\",style=\"filled\"];\n");
		}
		System.out.println("}");
		
		// Connecting clusters to RepBase leaves
		neighbors = new RepBase.RepbaseTreeNode[RepBase.string2node.length];
		isActive = new boolean[RepBase.string2node.length];
		idGenerator=RepBase.nNodes-1;
		for (i=firstFile; i<=lastFile; i++) {
			Math.set(histograms,0);
			path=aligmentFilesPath+"/"+i+suffix;
			file = new File(path);
			if (!file.exists()) continue;
			if (file.length()==0L) {
				idGenerator++;
				continue;
			}
			alignmentsFile = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
			str=alignmentsFile.readLine();
			str=alignmentsFile.readLine();  // Skipping the first two lines
			str=alignmentsFile.readLine();
			lastNeighbor=-1;
			while (str!=null) {
				Alignments.readAlignmentFile(str);
				neighbor=RepBase.string2node[Alignments.readB-1];
				found=-1;
				for (j=0; j<=lastNeighbor; j++) {
					if (neighbors[j].id==neighbor.id) {
						found=j;
						break;
					}
				}
				if (found!=-1) {
					neighbors[found].count++;
					for (j=Alignments.startB; j<=Alignments.endB; j++) histograms[found][j]++;
				}
				else {
					neighbors[++lastNeighbor]=neighbor;
					neighbor.count=1;
					for (j=Alignments.startB; j<=Alignments.endB; j++) histograms[lastNeighbor][j]++;
				}
				str=alignmentsFile.readLine();
			}
			alignmentsFile.close();
			idGenerator++;
			if (lastNeighbor==-1) continue;
			
			// Printing histograms
			for (j=0; j<=lastNeighbor; j++) {
				System.err.println(idGenerator+" -> "+neighbors[j].label);
				for (k=0; k<neighbors[j].lengthMax; k++) System.err.print(histograms[j][k]+",");
				System.err.println();
			}
			
			// Removing arcs to nodes that have just one descendant in $neighbors$
			for (j=0; j<=lastNeighbor; j++) isActive[j]=true;
			for (j=0; j<=lastNeighbor; j++) {
				descendant=-1;
				for (k=0; k<=lastNeighbor; k++) {
					if (k==j) continue;
					tmpNode=neighbors[k].parent;
					while (tmpNode!=null) {
						if (tmpNode.id==neighbors[j].id) {
							if (descendant==-1) descendant=k;
							else if (descendant>=0 && descendant!=k) descendant=-2; 
							break;
						}
						tmpNode=tmpNode.parent;
					}
				}
				if (descendant>=0) {
					isActive[j]=false;
					neighbors[descendant].count+=neighbors[j].count;
				}
			}
			
			// Assigning weights to arcs
			totalCount=0.0;
			for (j=0; j<=lastNeighbor; j++) {
				if (!isActive[j]) continue;
				totalCount+=neighbors[j].count;
			}
			for (j=0; j<=lastNeighbor; j++) {
				if (!isActive[j]) continue;
				System.out.println(idGenerator+" -> "+neighbors[j].id+" [color=\"red\",minlen=\"10\",penwidth=\""+IO.format((MAX_PEN_WIDTH*neighbors[j].count)/totalCount)+"\"];\n");
			}
		}
	}


	/**
	 * Prints, in the current directory, a heatmap on the RepBase topology induced by each
	 * interval type. Each interval of a given type contributes 1/N to each of the N
	 * distinct RepBase nodes it aligns to, independent of the number of alignments to 
	 * each.
	 */
	private static final void printRepbaseHeatmaps(int firstFile, int lastFile, String fragmentFilesPath, String fragmentSuffix, String aligmentFilesPath, String alignmentSuffix) throws IOException {
		int i, j, k, h;
		int readA, previousReadA, readB, previousReadB, lastNode;
		String str, path;
		BufferedReader fragmentsFile, alignmentsFile;
		BufferedWriter bw;
		File file;
		double[] tmp;
		int[][] types;  // column 0: type of each interval; 1: n. of components to which the interval belongs.
		String[] tokens;
		RepBase.RepbaseTreeNode[] nodes;
		
		// Collecting type counts
		RepBase.root.initTypeCounts();
		nodes = new RepBase.RepbaseTreeNode[RepBase.nNodes];
		for (i=firstFile; i<=lastFile; i++) {			
			// Loading interval types
			path=fragmentFilesPath+"/"+i+fragmentSuffix;
			file = new File(path);
			if (!file.exists()) continue;
			fragmentsFile = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
			str=fragmentsFile.readLine();
			j=0;
			while (str!=null) {
				j++;
				str=fragmentsFile.readLine();
			}
			fragmentsFile.close();
			types = new int[j][2];
			fragmentsFile = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
			str=fragmentsFile.readLine();
			j=-1;
			while (str!=null) {
				tokens=str.split(",");
				j++;
				types[j][0]=Integer.parseInt(tokens[3]);
				h=0;
				for (k=4; k<tokens.length; k++) {
					if (tokens[k].length()>0) h++;
				}
				types[j][1]=h==0?1:h;
				str=fragmentsFile.readLine();
			}
			fragmentsFile.close();
			
			// Reading alignments
			path=aligmentFilesPath+"/"+i+alignmentSuffix;
			file = new File(path);
			if (!file.exists() || file.length()==0L) continue;
			alignmentsFile = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
			str=alignmentsFile.readLine();
			str=alignmentsFile.readLine();  // Skipping the first two lines
			str=alignmentsFile.readLine();
			lastNode=-1; previousReadA=-1; previousReadB=-1;
			while (str!=null) {
				Alignments.readAlignmentFile(str);
				readA=Alignments.readA;
				readB=Alignments.readB;
				if (readA!=previousReadA) {
					for (j=0; j<=lastNode; j++) nodes[j].typeCounts[types[previousReadA-1][0]]+=1.0/((lastNode+1)*types[previousReadA-1][1]);
					lastNode=0;
					nodes[lastNode]=RepBase.string2node[readB-1];
					previousReadA=readA;
				}
				else if (readB!=previousReadB) {
					nodes[++lastNode]=RepBase.string2node[readB-1];
					previousReadB=readB;
				}
				str=alignmentsFile.readLine();
			}
			for (j=0; j<=lastNode; j++) nodes[j].typeCounts[types[previousReadA-1][0]]+=1.0/((lastNode+1)*types[previousReadA-1][1]);
			alignmentsFile.close();
		}
		
		// Normalizing counts
		tmp = new double[Constants.INTERVAL_PERIODIC+1];
		RepBase.root.getMaxCounts(tmp);
		RepBase.root.divideCounts(tmp);
		
		// Printing
		for (i=0; i<=Constants.INTERVAL_PERIODIC; i++) {
			if (i==Constants.INTERVAL_DENSE_SUFFIX) continue;
			bw = new BufferedWriter(new FileWriter("./repbaseHeatmap-type"+i+".dot"),IO.BUFFER_SIZE);
			bw.write("digraph G {\n");
			RepBase.root.mark(i);
			RepBase.root.toDot(i,bw);
			bw.write("}\n");
			bw.close();
		}
	}


}