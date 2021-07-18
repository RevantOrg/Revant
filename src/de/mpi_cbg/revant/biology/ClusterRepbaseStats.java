package de.mpi_cbg.revant.biology;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.Alignments;


/**
 * Computes statistics on the purity of clusters with respect to the RepBase tree.
 */
public class ClusterRepbaseStats {
	
	private static File[] files;
	private static int lastFile;
	private static int[] countHistogram;
	private static int lastCount;
	private static double[] entropies;
	private static int lastEntropy;
	
	
	/**
	 * java ClusterRepbaseStats ./lashow ./repbase/repbaseTree.txt ./repbase/repbaseConcatenation.ref ./stats/countHistogram.txt ./stats/entropies.txt 0.2 0.8
	 */
	public static void main(String[] args) throws IOException {
		final String ALIGNMENTS_DIR = args[0];
		final String REPBASE_TREE_FILE = args[1];
		final String REPBASE_FASTA_FILE = args[2];
		final String COUNT_FILE = args[3];
		final String ENTROPY_FILE = args[4];
		final double ENTROPY_MIN = Double.parseDouble(args[5]);
		final double ENTROPY_MAX = Double.parseDouble(args[6]);
		int i, maxDepth;
		
		IO.initialize();
		maxDepth=RepBase.loadRepbase(REPBASE_TREE_FILE,REPBASE_FASTA_FILE,false);
		for (i=0; i<=maxDepth; i++) {
			getStats(ALIGNMENTS_DIR,i,ENTROPY_MIN,ENTROPY_MAX);
			printStats(COUNT_FILE+"-"+i,ENTROPY_FILE+"-"+i);
		}
	}
	
	
	/**
	 * Prints to stdout clusters with entropy in $[entropyMin..entropyMax)$.
	 */
	private static final void getStats(String aligmentsDir, int treeDepth, double entropyMin, double entropyMax) throws IOException {
		final int MAX_NODES_PER_CLUSTER = 100;
		int i, j;
		int found, lastNeighbor;
		double entropy;
		String str;
		BufferedReader alignmentsFile;
		RepBase.RepbaseTreeNode neighbor;
		RepBase.RepbaseTreeNode[] neighbors;
		
		// Counting
		countHistogram = new int[MAX_NODES_PER_CLUSTER]; lastCount=-1;
		neighbors = new RepBase.RepbaseTreeNode[RepBase.string2node.length];
		files = (new File(aligmentsDir)).listFiles(); lastFile=-1;
		entropies = new double[files.length]; lastEntropy=-1;
		for (i=0; i<files.length; i++) {
			if (files[i].isDirectory() || files[i].length()==0L) continue;			
			lastFile++;
			alignmentsFile = new BufferedReader(new FileReader(files[i].getPath()));
			str=alignmentsFile.readLine();
			str=alignmentsFile.readLine();  // Skipping the first two lines
			str=alignmentsFile.readLine();
			lastNeighbor=-1;
			while (str!=null) {
				Alignments.readAlignmentFile(str);
				neighbor=RepBase.string2node[Alignments.readB-1];
				while (neighbor.depth>treeDepth) neighbor=neighbor.parent;
				found=-1;
				for (j=0; j<=lastNeighbor; j++) {
					if (neighbors[j].id==neighbor.id) {
						found=j;
						break;
					}
				}
				if (found!=-1) neighbors[found].count++;
				else {
					neighbors[++lastNeighbor]=neighbor;
					neighbor.count=1;
				}
				str=alignmentsFile.readLine();
			}
			alignmentsFile.close();
			countHistogram[lastNeighbor+1]++;
			if (lastNeighbor+1>lastCount) lastCount=lastNeighbor+1;
			if (lastNeighbor>=0) {
				entropy=entropy(neighbors,lastNeighbor);
				entropies[++lastEntropy]=entropy;
				if (entropy>=entropyMin && entropy<entropyMax) System.out.println("treeDepth "+treeDepth+"> entropy of file "+files[i].getName()+" is "+entropy+" in ["+entropyMin+".."+entropyMax+")");
			}
		}
		
		// Removing empty files from $files$
		i=-1;
		for (j=0; j<files.length; j++) {
			if (files[j].isDirectory() || files[j].length()==0L) continue;
			if (j!=i) files[++i]=files[j];
		}
	}
	
	
	private static final double entropy(RepBase.RepbaseTreeNode[] neighbors, int lastNeighbor) {
		if (lastNeighbor==0) return 0;
		int i;
		double n, p, entropy;		
		
		n=0.0;
		for (i=0; i<=lastNeighbor; i++) n+=neighbors[i].count;
		entropy=0.0;
		for (i=0; i<=lastNeighbor; i++) {
			p=neighbors[i].count/n;
			entropy+=-p*Math.log(p);
		}
		return entropy;
	}
	
	
	private static final void printStats(String countFile, String entropyFile) throws IOException {
		int i;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(countFile),IO.BUFFER_SIZE);
		for (i=0; i<=lastCount; i++) bw.write(i+","+countHistogram[i]+"\n");
		bw.close();
		
		bw = new BufferedWriter(new FileWriter(entropyFile),IO.BUFFER_SIZE);
		Arrays.sort(entropies);
		for (i=0; i<entropies.length; i++) bw.write(entropies[i]+"\n");
		bw.close();
	}

}