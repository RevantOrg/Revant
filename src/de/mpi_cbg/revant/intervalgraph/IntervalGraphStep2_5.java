package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.util.PriorityQueue;
import java.io.IOException;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Stream;
import de.mpi_cbg.revant.factorize.Alignments;


/**
 * Applies again the clustering procedure of $IntervalGraphStep2$, to a Step 2 component 
 * that results from the simplified splitting of a very large connected component.
 * The component in input might not be connected when considering containment edges only.
 */
public class IntervalGraphStep2_5 {
	
	/**
	 * java -Xms2G -Xmx4G IntervalGraphStep2_5 /Users/ramseysnow/Desktop/SchmidteaPolychroa_650Mb_PacBioSequel/step1/0-step2_0/92-step2/4.graph /Users/ramseysnow/Desktop/SchmidteaPolychroa_650Mb_PacBioSequel/step1/0-step2_0/92-step2/4-step2_5 1 3 1000
	 */
	public static void main(String[] args) throws IOException {
		final String GRAPH_SUFFIX = ".graph";
		final String GRAPH_FILE=args[0];
		final String GRAPH_FILE_PREFIX = GRAPH_FILE.substring(0,GRAPH_FILE.length()-GRAPH_SUFFIX.length());
		final String COMPONENTS_DIR = args[1];
		IO.coverage=Integer.parseInt(args[2]);
		IO.minOccurrencesInGenome=Integer.parseInt(args[3]);
		IO.minRepeatCoverage=(IO.minOccurrencesInGenome*IO.coverage)-1;
		Alignments.minAlignmentLength=Integer.parseInt(args[4]);
		final String TAGS_DIR = args[5];
		final String TAGS_PREFIX = args[6];
		
		final int MIN_COMPONENT_SIZE = 3;  // Arbitrary lower bound for clustering a component
		final int MIN_CLIQUE_SIZE = 5;  // Arbitrary
		int i, nEdges, nComponents, nComponentsPrime, firstComponent;
		PriorityQueue queue;
		Stream stream;
		int[] tmp, components, new2old, subgraph, degrees, clusterField, out;
		int[][] componentSize;
		IO.initialize();

		System.err.println("STEP2.5> Loading interval graph...");
		new2old=IntervalGraph.deserialize(GRAPH_FILE,true,true);
		if (IO.CONSISTENCY_CHECKS) {
			System.err.println("STEP2.5> Consistency checks started...");
			IntervalGraph.checkConsistency(0,IntervalGraph.nNodes-1,new int[IntervalGraph.getMaxDegree()*10]);
			System.err.println("STEP2.5> Consistency checks passed");
		}
		IntervalGraph.printNew2OldArray(new2old,GRAPH_FILE_PREFIX+".new2old");
		new2old=null;
		IntervalGraph.turnOffEdges(0,false);  // Keeping only containment edges
		System.err.println("STEP2.5> Interval graph loaded:");
		IntervalGraph.getNodeStats();
		nEdges=IntervalGraph.getEdgeStats();
		System.err.println("STEP2.5> nNodes="+IntervalGraph.nNodes+" nEdges="+nEdges);
		System.err.println();
		IntervalGraph.printNodeStats();
		IntervalGraph.printEdgeStats(nEdges);
		System.err.println();
		if (IntervalGraph.nNodes<=MIN_COMPONENT_SIZE) return;
		
		System.err.println("STEP2.5> Computing connected components...");
		nComponents=IntervalGraph.getConnectedComponents(1,true,true);
		System.err.println("STEP2.5> "+nComponents+" connected components found");
		subgraph = new int[IntervalGraph.nNodes];
		for (i=0; i<IntervalGraph.nNodes; i++) subgraph[i]=i;
		degrees = new int[IntervalGraph.nNodes];
		clusterField = new int[IntervalGraph.nNodes];
		queue = new PriorityQueue();
		stream = new Stream(256);
		out = new int[2];
		firstComponent=nComponents;
		if (nComponents==1) {
			System.err.println("STEP2.5> Clustering the only connected component...");
			nComponentsPrime=IntervalGraphStep2.getComponents_clustering(subgraph,0,IntervalGraph.nNodes-1,firstComponent,degrees,clusterField,stream,out,true/* EDGES SORTED????? */,true/* GRAPHISSMALL????? */,MIN_CLIQUE_SIZE);
			if (nComponentsPrime>1) {
				System.err.println("STEP2.5> "+nComponentsPrime+" clusters found inside the only connected component");
				System.err.println("STEP2.5> Saving all frequent clusters and components...");
				components=IntervalGraphStep2.getComponentTags(nComponents+nComponentsPrime);
				IntervalGraph.serializeFrequentComponents(components,components.length,true,true,COMPONENTS_DIR,".graph",TAGS_DIR,TAGS_PREFIX);  // Strict, since we will use just nodes with one component for building kernels downstream.				
			}
			else System.err.println("STEP2.5> Just one cluster found inside the only connected component");
		}
		else {
			componentSize = new int[nComponents][4];
			componentSize[0][0]=IntervalGraph.nNodes; componentSize[0][1]=0; componentSize[0][2]=0;
			tmp = new int[IntervalGraph.nNodes];
			IntervalGraph.sortNodesByComponent(subgraph,0,IntervalGraph.nNodes-1,0,nComponents,1,componentSize,tmp);
			tmp=null;
			for (i=0; i<nComponents-1; i++) {  // Keeping also components of size one
				if (componentSize[i][0]<=MIN_COMPONENT_SIZE) continue;
				System.err.println("STEP2.5> Clustering component "+i+"...");
				nComponentsPrime=IntervalGraphStep2.getComponents_clustering(subgraph,componentSize[i][1],componentSize[i+1][1]-1,firstComponent,degrees,clusterField,stream,out,true/* EDGES SORTED?????? */,true/* GRAPHISSMALL????? */,MIN_CLIQUE_SIZE);
				if (nComponentsPrime>1) {
					System.err.println("STEP2.5> "+nComponentsPrime+" clusters found inside component "+i);
					firstComponent+=nComponentsPrime;
				}
				else System.err.println("STEP2.5> Just one cluster found inside component "+i);
				queue.clear(); stream.clear(false);
			}
			if (componentSize[nComponents-1][0]>MIN_COMPONENT_SIZE) {
				System.err.println("STEP2.5> Clustering component "+(nComponents-1)+"...");
				nComponentsPrime=IntervalGraphStep2.getComponents_clustering(subgraph,componentSize[nComponents-1][1],IntervalGraph.nNodes-1,firstComponent,degrees,clusterField,stream,out,true/* EDGES SORTED?????? */,true/* GRAPHISSMALL????? */,MIN_CLIQUE_SIZE);
				if (nComponentsPrime>1) {
					System.err.println("STEP2.5> "+nComponentsPrime+" clusters found inside component "+(nComponents-1));
					firstComponent+=nComponentsPrime;
				}
				else System.err.println("STEP2.5> Just one cluster found inside component "+(nComponents-1));
			}
			System.err.println("STEP2.5> Saving all frequent clusters and components...");
			components=IntervalGraphStep2.getComponentTags(firstComponent);
			IntervalGraph.serializeFrequentComponents(components,components.length,true,true,COMPONENTS_DIR,".graph",TAGS_DIR,TAGS_PREFIX);  // Strict, since we will use just nodes with one component for building kernels downstream.
		}
		queue=null; stream.clear(true); stream=null;		
	}

}