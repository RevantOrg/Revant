package de.mpi_cbg.revant.intervalgraph;

import java.io.*;
import java.util.Arrays;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.Leaf;
import de.mpi_cbg.revant.util.Leaves;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Factorize;
import de.mpi_cbg.revant.factorize.PeriodicSubstrings;


/**
 * The program returns with exit status 255 iff all path kernels in input should be kept
 * in output (but it does not perform any copying itself), otherwise it exits with
 * standard exit statuses.
 *
 * Remark: fields $K,P$ of the basin descriptor headers printed by this program carry 
 * no information.
 * 
 * Remark: one might discard arcs of the kernel graph that connect path kernels in the 
 * same cluster, since redundancy between such path kernels should have already been 
 * detected by $IntervalGraphStep3$. However, in practice there is still redundancy, and
 * the bigger basins produced by $IntervalGraphStep3$ allow to detect it, so the program
 * takes the aggressive approach of not removing same-cluster arcs. 
 *
 * Remark: this program does not consider cases in which a single kernel should have been
 * assembled across distinct clusters. This is by design, and it might be fixed in the 
 * following steps of the pipeline.
 *
 * Remark: every $tags-X-step4.txt$ file built by this program contains a single number 
 * $X$ in the $prefix$ field of its rows, and its kernel tag IDs have no relation to the 
 * kernel tags IDs of previous steps. Thus, the file should be interpreted using the 
 * $X-newKernel2oldKernel.txt$ file printed by this program.
 */
public class IntervalGraphStep4 {
	
	public static final int KEEP_ALL_PATH_KERNELS = 255;
	
	/**
	 * Assigns to each kernel ID, the ID of the basin it corresponds to.
	 */
	public static int[][] kernel2descriptor;
	
	/**
	 * Data structures used by $newTag2oldTags(),loadKernel2Descriptor()$.
	 */
	private static int nComponents, nRows;
	private static int[] components;  // Distinct component IDs
	private static int[] first;  // First position of $components[i]$ in the sorted list.
	private static int[] newTags;
	private static int[][] oldTags;
	
	/**
	 * Temporary space
	 */
	private static int[] tmpIO = new int[10];
	
	
	
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		final int COMPONENT_ID = Integer.parseInt(args[1]);
		final int N_ALIGNMENTS = Integer.parseInt(args[2]);
		final int N_READS = Integer.parseInt(args[3]);
		Alignments.minAlignmentLength=Integer.parseInt(args[4]);
		final String QUALITY_THRESHOLDS_FILE = args[5];
		IO.coverage=Integer.parseInt(args[6]);
		IO.minOccurrencesInGenome=Integer.parseInt(args[7]);
		final double MIN_PERIODIC_FRACTION = 0.9;  // Arbitrary
		
		final int MIN_NODES_PER_KERNEL = IO.minRepeatCoverage+1;
		final String DESCRIPTORS_DIR = STEP1_DIR+"/"+IO.TAGS_DIR;
		final String OUTPUT_DIR = DESCRIPTORS_DIR+"/"+IO.STEP4_DIR;
		final String ALIGNMENTS_FILE = STEP1_DIR+"/"+COMPONENT_ID+"-"+IO.ALIGNMENTS_FILE;
		final String SHORTPERIOD_BITVECTOR_FILE = STEP1_DIR+"/"+COMPONENT_ID+"-"+IO.ALIGNMENTS_SHORTPERIOD;
		final String READ_IDS_FILE = STEP1_DIR+"/"+COMPONENT_ID+"-"+IO.READS_IDS;
		final String READ_LENGTHS_FILE = STEP1_DIR+"/"+COMPONENT_ID+"-"+IO.READS_LENGTHS;
		final String QUALITIES_FILE = STEP1_DIR+"/"+COMPONENT_ID+"-"+IO.READS_PHRED_PREFIX+IO.getPhredSuffix(STEP1_DIR+"/"+COMPONENT_ID+"-"+IO.READS_PHRED_PREFIX);
		final String SHORTPERIOD_TRACKS_FILE = STEP1_DIR+"/"+COMPONENT_ID+"-"+IO.READS_SHORTPERIOD;
		final String GRAPH_FILE = STEP1_DIR+"/"+COMPONENT_ID+".graph";
		
		int i;
		int out, minFrequency, nRepeats;
		byte[] kernelLabels;
		int[] lastNode2kernel;
		int[][] nodes2kernel, newKernelFrequency;
		IntervalGraph.Node[] nodes;
		
		// Building the kernel graph
		IO.initialize();
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		IntervalGraphStep3.maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadQualities(QUALITY_THRESHOLDS_FILE,(new File(QUALITIES_FILE)).exists()?QUALITIES_FILE:null);
		out=loadNodes(DESCRIPTORS_DIR,IO.BASIN_DESCRIPTOR_PREFIX+"-"+COMPONENT_ID+"-",GRAPH_FILE,true,true,false);
		if (out==-1 || out==-2) System.exit(KEEP_ALL_PATH_KERNELS);
		IntervalGraph.nPeriodicNodes(tmpIO);
		if (tmpIO[0]+tmpIO[1]>=IntervalGraph.nNodes*MIN_PERIODIC_FRACTION) {
			// Conservatively keeping all periodic kernels; this is because we want to
			// keep periodic kernels from different Step2 clusters distinct (otherwise the
			// Step2 clustering itself might have been useless).
			System.exit(KEEP_ALL_PATH_KERNELS);
		}
		System.err.println("buildKernelGraph> "+out+" initial path kernels marked as periodic");
		IntervalGraphStep3.shortPeriodAlignments=IntervalGraphStep3.loadShortPeriodBitvector(N_ALIGNMENTS,ALIGNMENTS_FILE,SHORTPERIOD_BITVECTOR_FILE);
		IntervalGraphStep3.shortPeriodWindows=IntervalGraph.loadShortPeriodTrack(SHORTPERIOD_TRACKS_FILE);
		IntervalGraphStep3.trackPointers=IntervalGraph.buildTrackPointers(IntervalGraphStep3.shortPeriodWindows);
		IntervalGraphStep3.buildPathKernelLongestNode();
		nodes2kernel = new int[IntervalGraphStep3.nNodesInKernel][0];
		lastNode2kernel = new int[IntervalGraphStep3.nNodesInKernel];
		nodes = new IntervalGraph.Node[IntervalGraphStep3.nNodesInKernel];
		for (i=0; i<IntervalGraphStep3.nNodesInKernel; i++) nodes[i]=IntervalGraph.nodesArray[IntervalGraphStep3.nodesInKernel[i]];
		IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
		if (IntervalGraphStep3.nNodesInKernel>1) Arrays.sort(nodes,0,IntervalGraphStep3.nNodesInKernel);
		IntervalGraphStep3.buildNodes2Kernel(nodes,IntervalGraphStep3.nNodesInKernel,ALIGNMENTS_FILE,nodes2kernel,lastNode2kernel,N_ALIGNMENTS);		
		IntervalGraphStep3.nPathKernelPeriodic=IntervalGraphStep3.buildKernelGraph_main(N_ALIGNMENTS,ALIGNMENTS_FILE,nodes,IntervalGraphStep3.nNodesInKernel,nodes2kernel,lastNode2kernel,out);
		IntervalGraphStep3.shortPeriodAlignments=null; 
		IntervalGraphStep3.shortPeriodWindows=null;
		IntervalGraphStep3.trackPointers=null;
		IntervalGraphStep3.basinStats(IntervalGraphStep3.tmpIO,true);
		IntervalGraphStep3.updateNKernelNodes();
if (IO.SHOW_INTERACTIVE) IntervalGraphStep3.printKernelGraph("./kernelGraph-before.dot",null,IntervalGraphStep3.kernelFrequency,kernel2descriptor);
		
		// Finding redundant path kernels using the kernel graph
		minFrequency=IntervalGraphStep3.getKernelFrequencyThreshold();
		System.err.println("STEP4> minFrequency="+minFrequency+"  nodes in the interval graph: "+IntervalGraph.nNodes);
		IntervalGraphStep3.buildFrequencyGraph_prefixSuffix(); 
		IntervalGraphStep3.buildFrequencyGraph_full();
		newKernelFrequency = new int[IntervalGraphStep3.nPathKernels][4];
		kernelLabels = new byte[IntervalGraphStep3.nPathKernels];
		nRepeats=IntervalGraphStep3.markRepeats(minFrequency,kernelLabels,newKernelFrequency); // It might happen that no path kernel gets a positive number in $kernelLabels$, and that some interval graph nodes get no kernel tag. 
		IntervalGraphStep3.deallocateFrequencyGraphs();
		System.err.println("STEP4> found "+nRepeats+" independent repeats (out of "+IntervalGraphStep3.nPathKernels+" path kernels)");
		if (IntervalGraphStep3.nPathKernelPeriodic>0 || IntervalGraphStep3.nPathKernelCyclic>0) {
			// Conservatively keeping all periodic kernels; this is because we want to
			// keep periodic kernels from different Step2 clusters distinct (otherwise the
			// Step2 clustering itself might have been useless).
			nRepeats+=IntervalGraphStep3.markRepeats_pathKernelPeriodic(kernelLabels,true);
		}
		if (nRepeats==0) System.exit(KEEP_ALL_PATH_KERNELS);  // Conservatively keeping all existing kernels
		IntervalGraphStep3.updateBidirectedGraphNodePointer();
		IntervalGraphStep3.basinStats(IntervalGraphStep3.tmpIO,false);
		nRepeats=IntervalGraphStep3.discardRareKernels(MIN_NODES_PER_KERNEL);
		if (nRepeats==0) System.exit(KEEP_ALL_PATH_KERNELS);  // Conservatively keeping all existing kernels
		IntervalGraphStep3.basinStats(IntervalGraphStep3.tmpIO,false);
		IntervalGraphStep3.getFrequencyPeriodic();
if (IO.SHOW_INTERACTIVE) IntervalGraphStep3.printKernelGraph("./kernelGraph-after.dot",kernelLabels,newKernelFrequency!=null?newKernelFrequency:IntervalGraphStep3.kernelFrequency,kernel2descriptor);
		IntervalGraphStep3.printKernelTags(OUTPUT_DIR+"/"+IO.TAGS_PREFIX+COMPONENT_ID+"-step4.txt",COMPONENT_ID+"",true,false,-1);		
		IntervalGraphStep3.printBasinDescriptors(kernelLabels,OUTPUT_DIR,IO.BASIN_DESCRIPTOR_PREFIX,kernel2descriptor,false);
		Points.deallocateMemory();
		printKernel2descriptor(kernelLabels,COMPONENT_ID,OUTPUT_DIR);
		System.exit(0);
	}

	
	/**
	 * Initializes the data structures of $IntervalGraph$ and $IntervalGraphStep3$ that 
	 * are necessary for building the kernel graph.
	 *
	 * Remark: it is important to load field $bidirectedGraphNode$ of every node as well,
	 * since it is used for building the Step4 kernel graph (see e.g.
	 * $buildBidirectedGraph_impl_allFlags(), markRepeats_markRedundantKernels()$ in 
	 * $IntervalGraphStep3$).
	 *
	 * Remark: the procedure loads all nodes in the graph of a connected component, 
	 * rather than just the nodes in the basin descriptor files. This is done in order to
	 * print tags correctly in $IntervalGraphStep3.printKernelTags()$, i.e. to print also 
	 * nodes that do not belong to any Step 3 basin, and to possibly merge them with other 
	 * intervals. Nodes that were written by Step 2 to the tags file of a connected 
	 * component, and never merged with Step 3 intervals, are now printed in output and 
	 * they might be merged with Step 3 intervals.
	 *
	 * Remark: some nodes in a basin descriptor might not be present in the graph of the
	 * connected component.
	 *
	 * Remark: kernels with unknown period are classified as short-period. 
	 *
	 * @param graphFile nodes in the file are assumed to be sorted by $read,start,nodeID$;
	 * @param abortIfOneKernel the procedure stops immediately if there is just one path 
	 * kernel;
	 * @param loadNodesInKernel if TRUE, initializes also $nodesInKernel,blockKernels$ in
	 * $IntervalGraphStep3$;
	 * @param multipleCalls does not reuse $Node$ objects in $IntervalGraph.nodesArray$,
	 * if any, so that the nodes produced by multiple calls are distinct objects;
	 * @return -1 if there is just one path kernel; -2 if there is no path kernel;
	 * otherwise, the number of path kernels marked in $pathKernelPeriodic$.
	 */
	public static final int loadNodes(String descriptorsDir, String descriptorPrefix, String graphFile, boolean abortIfOneKernel, boolean loadNodesInKernel, boolean multipleCalls) throws IOException {
		final int PREFIX_LENGTH = descriptorPrefix.length();
		int i, j, k, p, q, r;
		int last, nNodesInKernel, nPathKernels, lastNewNode;
		int currentKernel, lastKernel, out;
		IntervalGraph.Node currentNode;
		String str;
		BufferedReader br;
		File directory, file;
		int[] tmpArray = new int[20];
		String[] files;
		
		// Computing the number of path kernels
		directory = new File(descriptorsDir);
		files=directory.list(); nPathKernels=0;
		for (i=0; i<files.length; i++) {
			if (files[i].length()<=PREFIX_LENGTH || !files[i].substring(0,PREFIX_LENGTH).equalsIgnoreCase(descriptorPrefix)) continue;
			p=files[i].indexOf("-"); if (p<0) continue;
			q=files[i].indexOf("-",p+1); if (q<0) continue;
			r=files[i].indexOf("-",q+1); if (r<0) continue;
			nPathKernels++;
		}
		if (nPathKernels==0) return -2;
		if (nPathKernels==1 && abortIfOneKernel) return -1;
		IntervalGraphStep3.nPathKernels=nPathKernels;
		if (IntervalGraphStep3.pathKernelLengths==null || nPathKernels>IntervalGraphStep3.pathKernelLengths.length) IntervalGraphStep3.pathKernelLengths = new int[nPathKernels];
		if (IntervalGraphStep3.pathKernelPeriodic==null || nPathKernels>IntervalGraphStep3.pathKernelPeriodic.length) IntervalGraphStep3.pathKernelPeriodic = new byte[nPathKernels];
		if (IntervalGraphStep3.kernelSize==null || nPathKernels>IntervalGraphStep3.kernelSize.length) IntervalGraphStep3.kernelSize = new int[nPathKernels];
		if (IntervalGraphStep3.kernelFrequency==null || nPathKernels>IntervalGraphStep3.kernelFrequency.length) IntervalGraphStep3.kernelFrequency = new int[nPathKernels][4];
		if (IntervalGraphStep3.nKernelNodes==null || nPathKernels>IntervalGraphStep3.nKernelNodes.length) IntervalGraphStep3.nKernelNodes = new int[nPathKernels];
		if (IntervalGraphStep3.pathKernel2Kernel==null || nPathKernels>IntervalGraphStep3.pathKernel2Kernel.length) IntervalGraphStep3.pathKernel2Kernel = new int[nPathKernels];
		if (IntervalGraph.stack==null || nPathKernels>IntervalGraph.stack.length) IntervalGraph.stack = new int[nPathKernels];
		if (IntervalGraphStep3.tmpPoints==null || nPathKernels>IntervalGraphStep3.tmpPoints.length) {
			IntervalGraphStep3.tmpPoints = new Point[nPathKernels];
			for (i=0; i<IntervalGraphStep3.tmpPoints.length; i++) IntervalGraphStep3.tmpPoints[i] = new Point();
		}
		if (IntervalGraphStep3.tmpLeaves==null || nPathKernels>IntervalGraphStep3.tmpLeaves.length) {
			IntervalGraphStep3.tmpLeaves = new Leaf[nPathKernels];
			for (i=0; i<IntervalGraphStep3.tmpLeaves.length; i++) IntervalGraphStep3.tmpLeaves[i] = new Leaf();
		}
		for (i=0; i<nPathKernels; i++) IntervalGraphStep3.pathKernel2Kernel[i]=i;
		Points.allocateMemory(nPathKernels);
		Leaves.allocateMemory(nPathKernels);
		IntervalGraphStep3.frequencyPairs = new IntervalGraphStep3.FrequencyPair[100];  // Arbitrary, needed by $getFrequency$.
		for (i=0; i<IntervalGraphStep3.frequencyPairs.length; i++) IntervalGraphStep3.frequencyPairs[i] = new IntervalGraphStep3.FrequencyPair();
		//System.err.println("Step4.loadNodes> nPathKernels="+nPathKernels);		
		
		// Building $kernel2descriptor$.
		if (kernel2descriptor==null || nPathKernels>kernel2descriptor.length) kernel2descriptor = new int[nPathKernels][3];
		lastKernel=-1;
		for (i=0; i<files.length; i++) {
			if (files[i].length()<=PREFIX_LENGTH || !files[i].substring(0,PREFIX_LENGTH).equalsIgnoreCase(descriptorPrefix)) continue;
			p=files[i].indexOf("-"); if (p<0) continue;
			q=files[i].indexOf("-",p+1); if (q<0) continue;
			r=files[i].indexOf("-",q+1); if (r<0) continue;
			lastKernel++;
			kernel2descriptor[lastKernel][0]=Integer.parseInt(files[i].substring(p+1,q));
			kernel2descriptor[lastKernel][1]=Integer.parseInt(files[i].substring(q+1,r));
			kernel2descriptor[lastKernel][2]=Integer.parseInt(files[i].substring(r+1,files[i].indexOf(".",r+1)));
		}
		
		// Building $IntervalGraph.nodesArray$.
		IntervalGraph.deserialize_noEdges(graphFile,true,!multipleCalls);
		lastNewNode=IntervalGraph.nNodes-1;
		for (i=0; i<=lastNewNode; i++) {
			currentNode=IntervalGraph.nodesArray[i];
			currentNode.lastKernel=-1; currentNode.lastPathWithStart=-1; currentNode.lastPathWithEnd=-1;
			currentNode.bidirectedGraphNode=-1;
		}
		currentKernel=-1; nNodesInKernel=0; out=0;
		for (i=0; i<files.length; i++) {
			if (files[i].length()<=PREFIX_LENGTH || !files[i].substring(0,PREFIX_LENGTH).equalsIgnoreCase(descriptorPrefix)) continue;
			p=files[i].indexOf("-"); if (p<0) continue;
			q=files[i].indexOf("-",p+1); if (q<0) continue;
			r=files[i].indexOf("-",q+1); if (r<0) continue;
			currentKernel++;
			br = new BufferedReader(new FileReader(descriptorsDir+"/"+files[i]));
			str=br.readLine();
			IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
			IntervalGraphStep3.pathKernelLengths[currentKernel]=tmpArray[1];
			if (tmpArray[3]==Constants.INTERVAL_PERIODIC) {
				out++;
				if (tmpArray[4]>=PeriodicSubstrings.MIN_PERIOD_LONG) IntervalGraphStep3.pathKernelPeriodic[currentKernel]=(byte)2;
				else {
					// Unknown periods are classified as short
					IntervalGraphStep3.pathKernelPeriodic[currentKernel]=(byte)1;
				}
			}
			else IntervalGraphStep3.pathKernelPeriodic[currentKernel]=(byte)0;
			j=0; str=br.readLine();
			if (str!=null) {
				IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
				if (tmpArray[6]>=0) nNodesInKernel++;
			}
			while (str!=null) {
				if (j>=IntervalGraph.nNodes) {
					// Node not in the graph
					lastNewNode=loadNodes_nodeNotFound(lastNewNode,currentKernel,tmpArray);
					str=br.readLine(); 
					if (str!=null) {
						IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
						if (tmpArray[6]>=0) nNodesInKernel++;
					}
					continue;
				}
				currentNode=IntervalGraph.nodesArray[j];
				if (currentNode.read<tmpArray[0] || (currentNode.read==tmpArray[0] && currentNode.start<tmpArray[1])) {
					j++;
					continue;
				}
				if (currentNode.read>tmpArray[0] || (currentNode.read==tmpArray[0] && currentNode.start>tmpArray[1])) {
					// Node not in the graph
					lastNewNode=loadNodes_nodeNotFound(lastNewNode,currentKernel,tmpArray);
					str=br.readLine(); 
					if (str!=null) {
						IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
						if (tmpArray[6]>=0) nNodesInKernel++;
					}
					continue;
				}
				k=j;
				while (currentNode.read==tmpArray[0] && currentNode.start==tmpArray[1] && (currentNode.end!=tmpArray[2] || currentNode.type!=tmpArray[7])) {
					k++;
					if (k==IntervalGraph.nNodes) break;
					currentNode=IntervalGraph.nodesArray[k];
				}
				if (currentNode.read!=tmpArray[0] || currentNode.start!=tmpArray[1] || currentNode.end!=tmpArray[2] || currentNode.type!=tmpArray[7]) {
					// Node not in the graph
					lastNewNode=loadNodes_nodeNotFound(lastNewNode,currentKernel,tmpArray);
					str=br.readLine();
					if (str!=null) {
						IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
						if (tmpArray[6]>=0) nNodesInKernel++;
					}
					continue;
				}
				loadNodes_editNode(currentNode,currentKernel,tmpArray);
				str=br.readLine();
				if (str!=null) {
					IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
					if (tmpArray[6]>=0) nNodesInKernel++;
				}
			}
			br.close();
		}
		if (lastNewNode>=IntervalGraph.nNodes) {
			IntervalGraph.Node.order=IntervalGraph.Node.READ_START_NODEID;
			Arrays.sort(IntervalGraph.nodesArray,0,lastNewNode+1);
			for (i=0; i<=lastNewNode; i++) IntervalGraph.nodesArray[i].nodeID=i;
			IntervalGraph.nNodes=lastNewNode+1;
		}
		for (i=0; i<IntervalGraph.nNodes; i++) {
			currentNode=IntervalGraph.nodesArray[i];
			if (currentNode.lastPathWithStart>0) Arrays.sort(currentNode.pathsWithStart,0,currentNode.lastPathWithStart+1);
			if (currentNode.lastPathWithEnd>0) Arrays.sort(currentNode.pathsWithEnd,0,currentNode.lastPathWithEnd+1);
			// $currentNode.kernels$ is already sorted by construction.
		}
		System.err.println("loadNodes> nNodes="+IntervalGraph.nNodes+" nNodesInKernel="+nNodesInKernel);
		
		// Building $IntervalGraphStep3.nodesInKernel$.
		if (loadNodesInKernel) {
			IntervalGraphStep3.nNodesInKernel=nNodesInKernel;
			if (IntervalGraphStep3.nodesInKernel==null || nNodesInKernel>IntervalGraphStep3.nodesInKernel.length) IntervalGraphStep3.nodesInKernel = new int[nNodesInKernel];
			if (IntervalGraphStep3.blockKernels==null || nNodesInKernel>IntervalGraphStep3.blockKernels.length) IntervalGraphStep3.blockKernels = new int[nNodesInKernel];
			last=-1; currentKernel=-1;
			for (i=0; i<files.length; i++) {
				if (files[i].length()<=PREFIX_LENGTH || !files[i].substring(0,PREFIX_LENGTH).equalsIgnoreCase(descriptorPrefix)) continue;
				p=files[i].indexOf("-"); if (p<0) continue;
				q=files[i].indexOf("-",p+1); if (q<0) continue;
				r=files[i].indexOf("-",q+1); if (r<0) continue;
				currentKernel++;
				br = new BufferedReader(new FileReader(descriptorsDir+"/"+files[i]));
				str=br.readLine();  // Skipping header
				str=br.readLine(); j=0;
				if (str!=null) IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
				while (str!=null) {
					if (tmpArray[6]<0) {
						str=br.readLine();
						if (str!=null) IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
						continue;
					}
					currentNode=IntervalGraph.nodesArray[j];
					if (currentNode.read<tmpArray[0] || currentNode.start<tmpArray[1]) {
						j++;
						continue;
					}
					k=j;
					while (currentNode.end!=tmpArray[2] || currentNode.type!=tmpArray[7]) {
						k++;
						currentNode=IntervalGraph.nodesArray[k];
					}
					last++;
					IntervalGraphStep3.nodesInKernel[last]=k;
					IntervalGraphStep3.blockKernels[last]=currentKernel;
					str=br.readLine();
					if (str!=null) IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
				}
				br.close();
			}
		}
		
		// Setting $IntervalGraphStep3.nPathKernelCyclic$.
		j=0;
		for (i=0; i<IntervalGraphStep3.nPathKernels; i++) {
			if (IntervalGraphStep3.pathKernelLengths[i]==0) j++;
		}
		IntervalGraphStep3.nPathKernelCyclic=j;
		
		return out;
	}
	
	
	/**
	 * Adds kernel $currentKernel$ to $currentNode$ using the information in $tmpArray$,
	 * which is assumed to be the output of $IntervalGraphStep3.
	 * readBasinDescriptorRecord()$.
	 *
	 * Remark: the procedure assumes that Step2 clusters do not share nodes. Even in this
	 * case, however, the same interval of the same type might occur in basin descriptors
	 * from different Step2 clusters, because of changes in the topology applied in Step3.
	 */
	private static final void loadNodes_editNode(IntervalGraph.Node currentNode, int currentKernel, int[] tmpArray) {
		final int GROWTH_RATE = 5;  // Arbitrary
		
		if (currentNode.lastKernel==-1) {
			if (currentNode.kernels==null || currentNode.kernels.length==0) currentNode.kernels = new int[GROWTH_RATE];
			if (currentNode.kernelOrientations==null || currentNode.kernelOrientations.length==0) currentNode.kernelOrientations = new byte[GROWTH_RATE];
			currentNode.kernels[0]=currentKernel;
			currentNode.lastKernel=0;
			currentNode.kernelOrientations[0]=(byte)tmpArray[3];
			if (tmpArray[4]!=0) {
				if (currentNode.pathsWithStart==null || currentNode.pathsWithStart.length==0) currentNode.pathsWithStart = new int[GROWTH_RATE];
				currentNode.pathsWithStart[0]=tmpArray[4]==1?currentKernel:-1-currentKernel;
				currentNode.lastPathWithStart=0;
			}
			else currentNode.lastPathWithStart=-1;
			if (tmpArray[5]!=0) {
				if (currentNode.pathsWithEnd==null || currentNode.pathsWithEnd.length==0) currentNode.pathsWithEnd = new int[GROWTH_RATE];
				currentNode.pathsWithEnd[0]=tmpArray[5]==1?currentKernel:-1-currentKernel;
				currentNode.lastPathWithEnd=0;
			}
			else currentNode.lastPathWithEnd=-1;
			currentNode.bidirectedGraphNode=tmpArray[6];
		}
		else {
			if (IO.CONSISTENCY_CHECKS) {
				final int componentID = kernel2descriptor[currentKernel][0];
				final int clusterID = kernel2descriptor[currentKernel][1];
				int i, kernel;
				for (i=0; i<=currentNode.lastKernel; i++) {
					kernel=currentNode.kernels[i];
					if (kernel2descriptor[kernel][0]!=componentID || kernel2descriptor[kernel][1]!=clusterID) {
/*						System.err.println("loadNodes_editNode> WARNING: the following node belongs both to cluster "+componentID+"-"+clusterID+" and to cluster "+kernel2descriptor[kernel][0]+"-"+kernel2descriptor[kernel][1]);
						System.err.println(currentNode);
*/					}					
				}
			}
			currentNode.lastKernel++;
			if (currentNode.lastKernel==currentNode.kernels.length) {
				int[] newKernels = new int[currentNode.kernels.length+GROWTH_RATE];
				System.arraycopy(currentNode.kernels,0,newKernels,0,currentNode.kernels.length);
				currentNode.kernels=newKernels;
				byte[] newOrientations = new byte[currentNode.kernelOrientations.length+GROWTH_RATE];
				System.arraycopy(currentNode.kernelOrientations,0,newOrientations,0,currentNode.kernelOrientations.length);
				currentNode.kernelOrientations=newOrientations;
			}
			currentNode.kernels[currentNode.lastKernel]=currentKernel;
			currentNode.kernelOrientations[currentNode.lastKernel]=(byte)tmpArray[3];
			if (tmpArray[4]!=0) {
				currentNode.lastPathWithStart++;
				if (currentNode.pathsWithStart==null) currentNode.pathsWithStart = new int[GROWTH_RATE];
				else if (currentNode.lastPathWithStart==currentNode.pathsWithStart.length) {
					int[] newArray = new int[currentNode.pathsWithStart.length+GROWTH_RATE];
					System.arraycopy(currentNode.pathsWithStart,0,newArray,0,currentNode.pathsWithStart.length);
					currentNode.pathsWithStart=newArray;
				}
				currentNode.pathsWithStart[currentNode.lastPathWithStart]=tmpArray[4]==1?currentKernel:-1-currentKernel;
			}
			if (tmpArray[5]!=0) {
				currentNode.lastPathWithEnd++;
				if (currentNode.pathsWithEnd==null) currentNode.pathsWithEnd = new int[GROWTH_RATE];
				else if (currentNode.lastPathWithEnd==currentNode.pathsWithEnd.length) {
					int[] newArray = new int[currentNode.pathsWithEnd.length+GROWTH_RATE];
					System.arraycopy(currentNode.pathsWithEnd,0,newArray,0,currentNode.pathsWithEnd.length);
					currentNode.pathsWithEnd=newArray;
				}
				currentNode.pathsWithEnd[currentNode.lastPathWithEnd]=tmpArray[5]==1?currentKernel:-1-currentKernel;
			}
			if (tmpArray[6]>=0) {
				if (IO.CONSISTENCY_CHECKS) {
					if (currentNode.bidirectedGraphNode==-1) {
						final int componentID = kernel2descriptor[currentKernel][0];
						final int clusterID = kernel2descriptor[currentKernel][1];
						int i, kernel;
						for (i=0; i<=currentNode.lastKernel; i++) {
							kernel=currentNode.kernels[i];
							if (kernel!=currentKernel && kernel2descriptor[kernel][0]==componentID && kernel2descriptor[kernel][1]==clusterID) {
								System.err.println("loadNodes_editNode> ERROR: the following node belongs to the reference of kernel "+componentID+"-"+clusterID+"-"+kernel2descriptor[currentKernel][2]+", but it does not belong to the reference of kernel "+componentID+"-"+clusterID+"-"+kernel2descriptor[kernel][2]+":");
								System.err.println(currentNode);
								System.exit(1);
							}
						}
					}
					else if (tmpArray[6]!=currentNode.bidirectedGraphNode) {
						System.err.println("loadNodes_editNode> ERROR: the following node is mapped to different bidirected graph nodes by different kernels:  ");
						System.err.println("by kernel "+currentKernel+" -> "+tmpArray[6]);
						System.err.println("by some existing kernel -> "+currentNode.bidirectedGraphNode);
						System.err.println(currentNode);
						System.exit(1);
					}
				}
				currentNode.bidirectedGraphNode=tmpArray[6];
			}
		}
	}
	
	
	/**
	 * Adds a new node to the end of $IntervalGraph.nodesArray$, using the information in
	 * $tmpArray$, which is assumed to be the output of $IntervalGraphStep3.
	 * readBasinDescriptorRecord()$.
	 *
	 * @return the new value of $lastNewNode$ after the addition.
	 */
	private static final int loadNodes_nodeNotFound(int lastNewNode, int currentKernel, int[] tmpArray) {
		final int GROWTH_RATE = 5;  // Arbitrary
		final int GROWTH_RATE_NODESARRAY = 100;  // Arbitrary
		int i;
		IntervalGraph.Node currentNode;
		
		// Editing an existing node, if any.
		for (i=IntervalGraph.nNodes; i<=lastNewNode; i++) {
			currentNode=IntervalGraph.nodesArray[i];
			if (currentNode.read==tmpArray[0] && currentNode.start==tmpArray[1] && currentNode.end==tmpArray[2] && currentNode.type==tmpArray[7]) {
				loadNodes_editNode(currentNode,currentKernel,tmpArray);
				return lastNewNode;
			}
		}
		
		// Adding a new node
		lastNewNode++;
		if (lastNewNode==IntervalGraph.nodesArray.length) {
			IntervalGraph.Node[] newArray = new IntervalGraph.Node[IntervalGraph.nodesArray.length+GROWTH_RATE_NODESARRAY];
			System.arraycopy(IntervalGraph.nodesArray,0,newArray,0,lastNewNode);
			IntervalGraph.nodesArray=newArray;
		}
		currentNode = new IntervalGraph.Node(tmpArray[0],tmpArray[7],-1,tmpArray[1],tmpArray[2],tmpArray[9]==1,tmpArray[10]==1,false,0,0,0,tmpArray[8],false,0,0,0,0.0);
		currentNode.lastKernel=0;
		currentNode.kernels = new int[GROWTH_RATE];
		currentNode.kernels[0]=currentKernel;
		currentNode.kernelOrientations = new byte[GROWTH_RATE];
		currentNode.kernelOrientations[0]=(byte)tmpArray[3];
		if (tmpArray[4]!=0) {
			currentNode.pathsWithStart = new int[GROWTH_RATE];
			currentNode.pathsWithStart[0]=tmpArray[4]==1?currentKernel:-1-currentKernel;
			currentNode.lastPathWithStart=0;
		}
		else currentNode.lastPathWithStart=-1;
		if (tmpArray[5]!=0) {
			currentNode.pathsWithEnd = new int[GROWTH_RATE];
			currentNode.pathsWithEnd[0]=tmpArray[5]==1?currentKernel:-1-currentKernel;
			currentNode.lastPathWithEnd=0;
		}
		else currentNode.lastPathWithEnd=-1;
		currentNode.bidirectedGraphNode=tmpArray[6];
		IntervalGraph.nodesArray[lastNewNode]=currentNode;
		return lastNewNode;
	}
	
	
	/**
	 * Writes the map $i->kernel2descriptor[i]$ for all $i$ selected in $kernelLabels$.
	 */
	private static final void printKernel2descriptor(byte[] kernelLabels, int componentID, String outputDir) throws IOException {
		int i;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(outputDir+"/"+componentID+"-"+IO.KERNEL2DESCRIPTOR_FILE),IO.BUFFER_SIZE);
		for (i=0; i<IntervalGraphStep3.nPathKernels; i++) {
			if ((kernelLabels[i]>0 || kernelLabels[i]==-3) && IntervalGraphStep3.kernelSize[i]>0) bw.write(i+","+IO.printRow(kernel2descriptor[i])+"\n");
		}
		bw.close();
	}
	
	
	/**
	 * Loads in $components,first,newTags,oldTags$ all kernel2descriptor files in 
	 * $filesDir$, sorted by component and new tag.
	 *
	 * @param filesDir if NULL, the procedure just initializes all data structures to NULL
	 * and exits.
	 */
	public static final void loadKernel2Descriptor(String filesDir) throws IOException {
		final int SUFFIX_LENGTH = IO.KERNEL2DESCRIPTOR_FILE.length();
		int i, j, p, q;
		int length, newTag, component, old1, old2;
		String str;
		BufferedReader br;
		File directory;
		String[] files;
		Kernel2DescriptorRow[] rows;
		
		if (filesDir==null) {
			nRows=0; nComponents=0; components=null; newTags=null; oldTags=null;
			return;
		}
		
		// Computing the number of entries
		directory = new File(filesDir);
		files=directory.list(); nRows=0;
		for (i=0; i<files.length; i++) {
			length=files[i].length();
			if (length<=SUFFIX_LENGTH || !files[i].substring(length-SUFFIX_LENGTH).equalsIgnoreCase(IO.KERNEL2DESCRIPTOR_FILE)) continue;
			br = new BufferedReader(new FileReader(filesDir+"/"+files[i]));
			str=br.readLine();
			while (str!=null) {
				nRows++;
				str=br.readLine();
			}
			br.close();
		}
		if (nRows==0) {
			nComponents=0; components=null; newTags=null; oldTags=null;
			return;
		}
		
		// Loading entries
		rows = new Kernel2DescriptorRow[nRows];
		j=-1;
		for (i=0; i<files.length; i++) {
			length=files[i].length();
			if (length<=SUFFIX_LENGTH || !files[i].substring(length-SUFFIX_LENGTH).equalsIgnoreCase(IO.KERNEL2DESCRIPTOR_FILE)) continue;
			br = new BufferedReader(new FileReader(filesDir+"/"+files[i]));
			str=br.readLine();
			while (str!=null) {
				if (str.length()==0) {
					str=br.readLine();
					continue;
				}
				j++;
				p=str.indexOf(",");
				newTag=Integer.parseInt(str.substring(0,p));
				q=str.indexOf("-",p+1);
				component=Integer.parseInt(str.substring(p+1,q));
				p=q; q=str.indexOf("-",p+1);
				old1=Integer.parseInt(str.substring(p+1,q));
				p=q;
				old2=Integer.parseInt(str.substring(p+1));
				rows[j] = new Kernel2DescriptorRow(newTag,component,old1,old2);
				str=br.readLine();
			}
			br.close();
		}
		
		// Building the output
		if (nRows>1) Arrays.sort(rows);
		newTags = new int[nRows]; oldTags = new int[nRows][2];
		for (i=0; i<nRows; i++) {
			newTags[i]=rows[i].newTag;
			oldTags[i][0]=rows[i].old1; oldTags[i][1]=rows[i].old2;
		}
		nComponents=1;
		for (i=1; i<nRows; i++) {
			if (rows[i].component!=rows[i-1].component) nComponents++;
		}
		components = new int[nComponents]; first = new int[nComponents];
		j=0; components[0]=rows[0].component; first[0]=0;
		for (i=1; i<nRows; i++) {
			component=rows[i].component;
			if (component!=rows[i-1].component) {
				j++; components[j]=component; first[j]=i;
			}
		}
	}
	
	
	private static class Kernel2DescriptorRow implements Comparable {
		public int newTag, component, old1, old2;
		
		public Kernel2DescriptorRow(int n, int c, int o1, int o2) {
			newTag=n; component=c; old1=o1; old2=o2;
		}
		
		/**
		 * Sorts by $component,newTag$.
		 */
		public int compareTo(Object other) {
			Kernel2DescriptorRow otherRow = (Kernel2DescriptorRow)other;
			if (component<otherRow.component) return -1;
			else if (component>otherRow.component) return 1;
			if (newTag<otherRow.newTag) return -1;
			else if (newTag>otherRow.newTag) return 1;
			return 0;
		}
	}
	
	
	/**
	 * Remark: the procedure assumes that $loadKernel2Descriptor()$ has already been 
	 * executed.
	 *
	 * @return an array containing the two old tags associated with $component,newTag$, or
	 * NULL if $component$ is not found.
	 */
	public static final int[] newTag2oldTags(int component, int newTag) {
		if (nRows==0) return null;
		int i;
		
		i=Arrays.binarySearch(components,0,nComponents,component);
		if (i<0) return null;
		i=Arrays.binarySearch(newTags,first[i],i==nComponents-1?nRows:first[i+1],newTag);
		if (i<0) {
			System.err.println("newTag2oldTags> ERROR: tag "+newTag+" not found in component "+component);
			System.exit(1);
		}
		return oldTags[i];
	}
	
	
	/**
	 * @return the value of $oldTags[][0]$ shared by all the tags of $component$ in 
	 * $tags[firstTag..lastTag]$, or -1 if such tags do not share the same value.
	 */
	public static final int newTag2oldTag(int component, int[] tags, int firstTag, int lastTag) {
		if (nRows==0) return -1;
		int i, j;
		int c, from, to, out;
		
		c=Arrays.binarySearch(components,0,nComponents,component);
		if (c<0) return -1;
		from=first[c];
		to=c==nComponents-1?nRows:first[c+1];
		out=-1;
		for (i=firstTag; i<=lastTag; i++) {
			j=Arrays.binarySearch(newTags,from,to,tags[i]);
			if (j<0) {
				System.err.println("newTag2oldTags> ERROR: tag "+tags[i]+" not found in component "+component);
				System.exit(1);
			}
			if (out==-1) out=oldTags[j][0];
			else if (out!=oldTags[j][0]) return -1;
		}
		return out;
	}
	
	
	
	
	
	
	
	
	
	
	// --------------------------------- UNUSED CODE -------------------------------------
	
	/**
	 * (This procedure is not currently used)
	 * Removes all arcs $(u,v)$ of the kernel graph such that $u$ and $v$ have the same
	 * value of $kernel2descriptor[0]$. We assume that there is no redundancy between 
	 * path kernels in the same cluster (this should have been detected by 
	 * $IntervalGraphStep3$).
	 *
	 * @return the number of arcs removed by the procedure.
	 */
	private static final int buildKernelGraph_removeSameClusterArcs() {
		byte tmpB;
		int i, j, k, h, p, q;
		int neighbor, last, out, tmpI, cluster;
		final int nPathKernels = IntervalGraphStep3.nPathKernels;
		final int[] lastInNeighbor = IntervalGraphStep3.lastInNeighbor;
		final int[] lastOutNeighbor = IntervalGraphStep3.lastOutNeighbor;
		final int[][] inNeighbors = IntervalGraphStep3.inNeighbors;
		final int[][] outNeighbors = IntervalGraphStep3.outNeighbors;
		final byte[][] inNeighborFlags = IntervalGraphStep3.inNeighborFlags;
		
		out=0;
		for (i=0; i<nPathKernels; i++) {
			cluster=kernel2descriptor[i][0];
			last=lastInNeighbor[i];
			for (j=0; j<=last; j++) {
				neighbor=inNeighbors[i][j];
				if (kernel2descriptor[neighbor][0]==cluster) {
					out++;
					inNeighbors[i][j]=-1-neighbor;
					for (k=0; k<=lastOutNeighbor[neighbor]; k++) {
						if (outNeighbors[neighbor][k]==i) {
							outNeighbors[neighbor][k]=-1-outNeighbors[neighbor][k];
							break;
						}
					}
				}
			}
		}
		for (i=0; i<nPathKernels; i++) {
			k=-1;
			for (j=0; j<=lastInNeighbor[i]; j++) {
				if (inNeighbors[i][j]<0) continue;
				k++;
				if (k!=j) {
					tmpI=inNeighbors[i][k];
					inNeighbors[i][k]=inNeighbors[i][j];
					inNeighbors[i][j]=tmpI;
					p=k<<2; q=j<<2;
					for (h=0; h<=3; h++) {
						tmpB=inNeighborFlags[i][p+h];
						inNeighborFlags[i][p+h]=inNeighborFlags[i][q+h];
						inNeighborFlags[i][q+h]=tmpB;
					}
				}
			}
			lastInNeighbor[i]=k;
		}
		for (i=0; i<nPathKernels; i++) {
			k=-1;
			for (j=0; j<=lastOutNeighbor[i]; j++) {
				if (outNeighbors[i][j]<0) continue;
				k++;
				if (k!=j) {
					tmpI=outNeighbors[i][k];
					outNeighbors[i][k]=outNeighbors[i][j];
					outNeighbors[i][j]=tmpI;
				}
			}
			lastOutNeighbor[i]=k;
		}
		return out;
	}
	
}