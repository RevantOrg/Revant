package de.mpi_cbg.revant.intervalgraph;

import java.io.IOException;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Factorize;


/**
 * Uses all edges induced by all alignments to take connected components of the interval  
 * graph. All edges are stored in the output files.  
 * 
 * Remark: clustering by connected components is also done in \cite{pertea2003tigr}, as  
 * well as processing the resulting components in parallel downstream.
 *
 * Remark: in every serialized file in output, nodes are sorted by $read,start$, and
 * the $avgDiffs$ field of edges is a sum, not an avg.
 *
 * Remark: this program might print some tags files, but the intervals stored in there 
 * are not going to be merged with the intervals produced by the steps downstream (as
 * specified by $IntervalGraphStep3.printKernelTags()$). This is done to highlight
 * such intervals.
 */
public class IntervalGraphStep1 {	
	
	public static void main(String[] args) throws IOException {
		final String INPUT_DIR = args[0];
		final int N_READS = Integer.parseInt(args[1]);
		IO.coverage=Integer.parseInt(args[2]);
		IO.minOccurrencesInGenome=Integer.parseInt(args[3]);
		Alignments.minAlignmentLength=Integer.parseInt(args[4]);
		final boolean DISCARD_CONTAINED_IN_SHORT_PERIOD = Integer.parseInt(args[5])==1;
		final boolean DISCARD_CONTAINED_IN_LONG_PERIOD = Integer.parseInt(args[6])==1;
		final boolean DISCARD_PERIOD_IN_PERIOD = Integer.parseInt(args[7])==1;
		final boolean DISCARD_DENSE_IN_DENSE = Integer.parseInt(args[8])==1;
		final boolean FIX_DANGLING = Integer.parseInt(args[9])==1;
		final boolean FIX_UNASSIGNED = Integer.parseInt(args[10])==1;
		final boolean ADD_SAME_READ_EDGES = Integer.parseInt(args[11])==1;
		
		IO.minRepeatCoverage=(IO.minOccurrencesInGenome*IO.coverage)-1;
		final String ALIGNMENTS_FILE = INPUT_DIR+"/"+IO.ALIGNMENTS_FILE_STEP1;
		final String CONNECTION_FILE = INPUT_DIR+"/"+IO.CONNECTION_FILE;
		final String COMPONENTS_DIR = INPUT_DIR+"/"+IO.STEP1_DIR;
		final String PERIODIC_INTERVALS = INPUT_DIR+"/"+IO.PERIODIC_INTERVALS;
		final String DENSE_INTERVALS = INPUT_DIR+"/"+IO.DENSE_INTERVALS;
		final String ALIGNMENT_INTERVALS = INPUT_DIR+"/"+IO.ALIGNMENT_INTERVALS;
		final String TAGS_DIR = COMPONENTS_DIR+"/"+IO.TAGS_DIR;
		final String TAGS_PREFIX = "";
		final String ALIGNMENTS_SHORTPERIOD = INPUT_DIR+"/"+IO.ALIGNMENTS_SHORTPERIOD;
		final String READS_SHORTPERIOD = INPUT_DIR+"/"+IO.READS_SHORTPERIOD;
		final String READS_LENGTHS = INPUT_DIR+"/"+IO.READS_LENGTHS;
		final String READS_IDS = INPUT_DIR+"/"+IO.READS_IDS;
		final String READS_PHRED = INPUT_DIR+"/"+IO.READS_PHRED_PREFIX+IO.getPhredSuffix(INPUT_DIR+"/"+IO.READS_PHRED_PREFIX);
		
		int i;
		int lastPoint, nComponents, nPrintedComponents;
		int[] size, nMaximal;
		
		IO.initialize();
		Reads.nReads=N_READS;
		Reads.loadReadIDs(READS_IDS,N_READS);
		System.err.println("STEP1> Building interval graph...");
		IntervalGraph.maxAlignmentsPerRead=IntervalGraph.buildIntervalGraph(ALIGNMENTS_FILE,CONNECTION_FILE,PERIODIC_INTERVALS,DENSE_INTERVALS,ALIGNMENT_INTERVALS,DISCARD_CONTAINED_IN_SHORT_PERIOD,DISCARD_CONTAINED_IN_LONG_PERIOD,DISCARD_PERIOD_IN_PERIOD,DISCARD_DENSE_IN_DENSE,FIX_DANGLING,FIX_UNASSIGNED,ADD_SAME_READ_EDGES,TAGS_DIR,TAGS_PREFIX,true);
		System.err.println("STEP1> Interval graph built. maxAlignmentsPerRead="+IntervalGraph.maxAlignmentsPerRead);
		IntervalGraph.getNodeStats();
		IntervalGraph.printNodeStats();
		System.err.println("STEP1> Computing connected components...");
		IntervalGraph.turnOffEdges(-2,false);
		nComponents=IntervalGraph.getConnectedComponents(1,true,true);
		System.err.println("STEP1> "+nComponents+" connected components found");
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
		System.err.println("STEP1> Saving all frequent components...");
		nPrintedComponents=IntervalGraph.serializeFrequentComponents(null,nComponents,false,false,COMPONENTS_DIR,".graph",TAGS_DIR,TAGS_PREFIX);
		System.err.println("STEP1> "+nPrintedComponents+" components saved");
		
		System.err.println("STEP1> Building short-period bitvector and tracks...");
		IntervalGraph.buildShortPeriodBitvector(ALIGNMENTS_FILE,ALIGNMENTS_SHORTPERIOD);
		IntervalGraph.buildShortPeriodTrack(READS_SHORTPERIOD);
		System.err.println("STEP1> done");
		
		System.err.println("STEP1> Filtering read and alignment files by component...");
		IntervalGraph.filterAlignmentsAndReads(null,nComponents,false,READS_LENGTHS,READS_PHRED,READS_SHORTPERIOD,ALIGNMENTS_FILE,ALIGNMENTS_SHORTPERIOD,COMPONENTS_DIR);
		System.err.println("STEP1> done");
	}
	
}