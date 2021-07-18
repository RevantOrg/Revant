package de.mpi_cbg.revant.intervalgraph;

import java.io.IOException;
import java.io.FileWriter;
import java.io.BufferedWriter;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.factorize.Alignments;


/**
 * Discards every edge whose type is not containment, and computes the connected 
 * components of the resulting graph. In particular, the procedure discards 
 * overlap edges, since they could be induced by repeats of type $Z = X_2 Y_1$ where 
 * $X = X_1 X_2$ and $Y = Y_1 Y_2$ are repeats: occurrences of $Z$ without $X$ and $Y$
 * overlap occurrences of $X$ and of $Y$, putting $X$ and $Y$ in the same connected 
 * component of the interval graph.
 *
 * Remark: each output file contains all edges that connect nodes in that component.
 * Edges that connect nodes in different components are discarded, and will thus be 
 * invisible to the following steps of the pipeline.
 */
public class IntervalGraphStep2_0 {
	/**
	 * Parameters of the pipeline
	 */
	private static int MIN_PEELING_SIZE = 1000;  // Minimum size of a component for it to be peeled. Will be estimated from the data.
	
	
	/**
	 * java -Xms2G -Xmx2G IntervalGraphStep2_0 ~/Desktop/components/0.graph ~/Desktop/components/0-step2_0 1 3 1000
	 */
	public static void main(String[] args) throws IOException {
		final String GRAPH_SUFFIX = ".graph";
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
		
		int i, lastPoint, nEdges, nComponents, nPrintedComponents;
		BufferedWriter bw;
		int[] size, nMaximal, components, new2old;
		IO.initialize();

		System.err.println("STEP2_0> Loading interval graph...");
		new2old=IntervalGraph.deserialize(GRAPH_FILE,true,true);
		if (IO.CONSISTENCY_CHECKS) {
			System.err.println("STEP2_0> Consistency checks started...");
			components=IntervalGraph.getComponents();
			if (components[0]!=1) {
				System.err.println("STEP2_0> ERROR: the nodes in the input graph have "+components[0]+" component tags.");
				System.exit(1);
			}
			IntervalGraph.checkConsistency(0,IntervalGraph.nNodes-1,new int[IntervalGraph.getMaxDegree()*10]);
			System.err.println("STEP2_0> Consistency checks passed");
		}
		IntervalGraph.printNew2OldArray(new2old,GRAPH_FILE_PREFIX+".new2old");
		new2old=null;
		System.err.println("STEP2_0> Interval graph loaded:");
		IntervalGraph.getNodeStats();
		nEdges=IntervalGraph.getEdgeStats();
		System.err.println("STEP2_0> nNodes="+IntervalGraph.nNodes+" nEdges="+nEdges);
		IntervalGraph.printNodeStats();
		IntervalGraph.printEdgeStats(nEdges);
		System.err.println();
		
		System.err.println("STEP2_0> Computing connected components of ON edges...");
		IntervalGraph.turnOffEdges(0,false);
		nComponents=IntervalGraph.getConnectedComponents(1,true,true);
		System.err.println("STEP2_0> "+nComponents+" components created by removing just non-containment edges");
		if (nComponents>1) MIN_PEELING_SIZE=estimateComponentSizeThreshold(nComponents,MIN_PEELING_SIZE);
		bw = new BufferedWriter(new FileWriter(GRAPH_FILE_PREFIX+".minPeelingSize"));
		bw.write(MIN_PEELING_SIZE+"\n");
		bw.close();
		System.err.println("Histogram of component sizes:");
		size = new int[nComponents];
		nMaximal = new int[nComponents];
		IntervalGraph.getComponentStats(null,nComponents,false,false,size,nMaximal);		
		IntervalGraph.ensureTmpPoints(nComponents);
		for (i=0; i<nComponents; i++) {
			IntervalGraph.tmpPoints[i].position=size[i];
			IntervalGraph.tmpPoints[i].mass=1;
		}
		lastPoint=Points.sortAndCompact(IntervalGraph.tmpPoints,nComponents-1);
		for (i=0; i<=lastPoint; i++) System.err.println(IntervalGraph.tmpPoints[i]);
		
		System.err.println("STEP2_0> Saving all frequent components...");
		nPrintedComponents=IntervalGraph.serializeFrequentComponents(components,nComponents,false,false,COMPONENTS_DIR,".graph",TAGS_DIR,TAGS_PREFIX);
		System.err.println("STEP2_0> "+nPrintedComponents+" components saved (out of total "+nComponents+")");
		if (PRINT_DOT) {
			System.err.println("STEP2_0> Printing the whole graph (with component tags) to <"+GRAPH_FILE_PREFIX+"-components.dot>...");
			IntervalGraph.toDot(GRAPH_FILE_PREFIX+"-components.dot",null,null,Math.POSITIVE_INFINITY);
		}
	}
	
	
	/**
	 * Assume that all edges of type different from containment have been removed from the
	 * interval graph, and that all nodes in the interval graph have been assigned exactly
	 * one connected component. The histogram of sizes of all connected components should 
	 * have a peak near one, corresponding to components that are too small to be split
	 * further. The procedure returns a value S that separates the peak near one from the 
	 * rest of the distribution. Only components of size bigger than S will be targets for
	 * peeling and splitting.
	 *
	 * @param maxSize only the region [0..maxSize] of the histogram of sizes is used.
	 * This is necessary, since there can be a large difference between the largest and
	 * the smallest size, and in this case the estimation of the threshold can be
	 * inaccurate. It is easy for the caller to come up with a loose estimate of 
	 * $maxSize$.
	 */
	private static final int estimateComponentSizeThreshold(int nComponents, int maxSize) {
		int i, threshold;
		int[] sizes;
		
		// Collecting component sizes
		sizes=IntervalGraphStep2.getComponentsSize(nComponents);
		if (IntervalGraphStep2.tmpPoints==null || IntervalGraphStep2.tmpPoints.length<nComponents) {
			IntervalGraphStep2.tmpPoints = new Point[nComponents];
			for (i=0; i<IntervalGraphStep2.tmpPoints.length; i++) IntervalGraphStep2.tmpPoints[i] = new Point();
		}
		IntervalGraphStep2.lastTmpPoint=-1;
		for (i=0; i<nComponents; i++) {
			if (sizes[i]>maxSize) continue;
			IntervalGraphStep2.lastTmpPoint++;
			IntervalGraphStep2.tmpPoints[IntervalGraphStep2.lastTmpPoint].position=sizes[i];
			IntervalGraphStep2.tmpPoints[IntervalGraphStep2.lastTmpPoint].mass=1;
		}
		sizes=null;
		
		// Estimating threshold
		IntervalGraphStep2.lastTmpPoint=Points.sortAndCompact(IntervalGraphStep2.tmpPoints,IntervalGraphStep2.lastTmpPoint);	
		threshold=(int)Points.getThreshold(IntervalGraphStep2.tmpPoints,IntervalGraphStep2.lastTmpPoint,0.0,true,-1);
		System.err.println("STEP2_0> Components sizes: ");
		for (i=0; i<=IntervalGraphStep2.lastTmpPoint; i++) System.err.println(IntervalGraphStep2.tmpPoints[i].position+","+IntervalGraphStep2.tmpPoints[i].mass);
		System.err.println("STEP2_0> Threshold="+threshold);
		return threshold;
	}
	
}