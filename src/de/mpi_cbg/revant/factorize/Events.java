package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import java.util.Random;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Clusters;
import de.mpi_cbg.revant.util.Leaf;

/**
 * Machinery for clustering arbitrary $Event$ objects.
 * In practice only maximal events are used.
 */
public class Events {
	/**
	 * Parameters of the pipeline
	 */
	public static int nKmeansIterations = 0;	
	public static double minMassRatio = 0.1;  // A run of local-maximum leaves whose mass is smaller than $minMassRatio$ times the maximum mass of a run is not assigned a cluster center
	public static double minKeptMassRatio = 0.5;  // A warning is raised if less than this fraction of events does not belong to output clusters
	public static int minIntervalLength = IO.quantum>>1;  // Intervals of a density estimation tree are not split if they are this length or smaller
	public static int minLocalMaxDistance = minIntervalLength;  // Minimum distance between two local maxima of a density estimation tree for them to be considered distinct peaks
	public static int smallLeafDistance = minIntervalLength>>2;  // Arbitrary

	/**
	 * Input data structures
	 */
	public static Event[] events;
	public static int lastEvent;
	public static int nDeletedEvents;  // Incremented by procedure $delete$

	/**
	 * Output data structures
	 */
	public static int[] nClusters;
	public static int[][] clusterSizes;
	public static double[][] clusterMeans;
	public static boolean[] keep;
	public static int bestIteration;

	/**
	 * Temporary space
	 */
	private static double[] averageDistances;  // The average distance of a point from each cluster
	private static double[] clusterVariances;
	private static double[][] clusterMassInStandardDeviation;  // The number of parentheses at distance at most one (0) and two (1) standard deviations from the mean.
	private static int[] nOpen, nClosed;
	private static Event query;


	public static final void allocateMemory(int maxEvents, int maxClusters) {
		// Input data structures
		events = new Event[maxEvents];
		for (int i=0; i<maxEvents; i++) events[i] = new Event();

		// Output data structures
		int iter = nKmeansIterations==0?1:nKmeansIterations;
		nClusters = new int[maxClusters];
		clusterSizes = new int[iter][maxClusters];
		clusterMeans = new double[iter][maxClusters];
		keep = new boolean[maxClusters];

		// Temporary space
		averageDistances = new double[maxClusters];
		clusterVariances = new double[maxClusters];
		clusterMassInStandardDeviation = new double[maxClusters][2];
		nOpen = new int[maxClusters];
		nClosed = new int[maxClusters];
		query = new Event();
	}


	/**
	 * Tries to cluster events, and to keep only clusters that resemble Gaussians.
	 * Specifically, assume that the set of events consists of short substrings with a
	 * high density of events, separated by longer substrings with a low density of
	 * events. The procedure fits a density estimation tree on the set of all events and,
	 * if $nKmeansIterations>0$, it uses it to randomly initialize $nKmeansIterations$ 
	 * iterations of k-means and to pick the best. If $nKmeansIterations=0$, the density 
	 * estimation tree is used directly to output a set of clusters, without running 
	 * k-means.
	 *
	 * Remark: the procedure does not assume $events$ to be sorted. It has the side effect
	 * of sorting and compacting $events$ (see procedure $Points.sortAndCompact$).
	 *
	 * Remark: this procedure does not use any knowledge of the domain, since clusters of
	 * events can have arbitrary properties. For example, a cluster could contain just
	 * open or just closed parentheses, but also a mixture of the two, or neither of the
	 * two since its events could come from a different source. Two consecutive clusters
	 * could both contain just open or just closed parentheses, i.e. they do not
	 * necessarily represent an open and a closed parenthesis.
	 *
	 * @return the k-means iteration that gives the best clustering (zero if
	 * $nIterations=0$); -1 if no clustering can be performed.
	 */
	public static final int cluster() {
		final int MAX_LARGEST_MASS = 20;  // Arbitrary
		final int MIN_HOLE_LENGTH = IO.quantum;  // Arbitrary
		int i;
		int nLocalMaximumLeaves, nKeptClusters, maxMass, largestMass;
		double length, massInKeptClusters, totalMass;
		double silhouette, bestSilhouette;

		totalMass=lastEvent+1;
		lastEvent=Points.sortAndCompact(events,lastEvent);
		for (i=0; i<=lastEvent; i++) events[i].isUsedForClustering=true;


if (IO.SHOW_STD_ERR_PRIME) IO.printErr("EVENTS:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=lastEvent; x++) IO.printErr(events[x]); }


		if (lastEvent==0) {
			clusterSizes[0][0]=events[lastEvent].getMass();
			clusterMeans[0][0]=events[lastEvent].position;
			nClusters[0]=1;
			bestIteration=0;
			events[0].cluster=0;
		}
		else {
			// Fitting density estimation tree
			largestMass=Points.largestMass(events,lastEvent,-1);
			if (largestMass>0) largestMass=Math.min(largestMass,MAX_LARGEST_MASS);
if (IO.SHOW_STD_ERR_PRIME) IO.printErr("Events.cluster> largestMass="+largestMass);
			nLocalMaximumLeaves=DensityEstimationTree.buildDensityEstimationTree(events,0,lastEvent,minIntervalLength,minLocalMaxDistance,true,largestMass,MIN_HOLE_LENGTH,smallLeafDistance,false,false);
			if (DensityEstimationTree.lastLeaf==0) {
				length=events[DensityEstimationTree.leaves[0].lastPoint].position-events[DensityEstimationTree.leaves[0].firstPoint].position+1;
				if (length<Intervals.jaccardThreshold*ReadA.readLength) {
					DensityEstimationTree.leaves[0].isLocalMaximum=true;
					nLocalMaximumLeaves=1;
				}
			}
			if (nLocalMaximumLeaves==0) {
if (IO.SHOW_STD_ERR_PRIME) IO.printErr(":'-(");
				return -1;
			}

if (IO.SHOW_STD_ERR_PRIME) IO.printErr("LEAVES OF THE D.E.T. ON EVENTS:");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr("firstPoint="+DensityEstimationTree.leaves[x].firstPoint+" lastPoint="+DensityEstimationTree.leaves[x].lastPoint+" isLocalMaximum="+DensityEstimationTree.leaves[x].isLocalMaximum+" marked="+DensityEstimationTree.leaves[x].marked+" VALUE="+DensityEstimationTree.leaves[x].value); }


			DensityEstimationTree.markRunsOfLocalMaximumLeaves(events);


if (IO.SHOW_STD_ERR_PRIME) IO.printErr("LEAVES OF THE D.E.T. ON EVENTS (after markRuns...):");
if (IO.SHOW_STD_ERR_PRIME) { for (int x=0; x<=DensityEstimationTree.lastLeaf; x++) IO.printErr("firstPoint="+DensityEstimationTree.leaves[x].firstPoint+" lastPoint="+DensityEstimationTree.leaves[x].lastPoint+" isLocalMaximum="+DensityEstimationTree.leaves[x].isLocalMaximum+" marked="+DensityEstimationTree.leaves[x].marked+" VALUE="+DensityEstimationTree.leaves[x].value); }


			// Clustering
			if (nKmeansIterations==0) {
				getRoughClustering();
				bestIteration=0;
			}
			else {
				DensityEstimationTree.discardLocalMinima(events,0,lastEvent);
				maxMass=DensityEstimationTree.getMassOfHeaviestRun(events);
				bestIteration=-1; bestSilhouette=-2;
				for (i=0; i<nKmeansIterations; i++) {
					initializeKmeans(i,maxMass);
					Clusters.kmeans(nClusters[i],events,lastEvent,clusterMeans[i],clusterSizes[i]);
					silhouette=Clusters.getSilhouette(events,lastEvent,nClusters[i],clusterSizes[i],averageDistances);
					if (silhouette>bestSilhouette) {
						bestSilhouette=silhouette;
						bestIteration=i;
					}
				}

				// Filtering clusters
				nKeptClusters=Clusters.filterClusters(events,lastEvent,nClusters[bestIteration],clusterSizes[bestIteration],clusterMeans[bestIteration],clusterVariances,clusterMassInStandardDeviation,keep);
				if (nKeptClusters==0) return -1;
				massInKeptClusters=0;
				for (i=0; i<nClusters[bestIteration]; i++) {
					if (keep[i]) massInKeptClusters+=clusterSizes[bestIteration][i];
				}
				if (massInKeptClusters<minKeptMassRatio*totalMass) { if (IO.SHOW_STD_ERR) IO.printErr("Warning: just the "+(massInKeptClusters/totalMass)+" fraction of all events to be clustered belong to the final clusters we keep."); }
			}
		}

		// Adding cluster centers to $Factors.splits$.
		Math.set(nOpen,nClusters[bestIteration]-1,0);
		Math.set(nClosed,nClusters[bestIteration]-1,0);
		for (i=0; i<=lastEvent; i++) {
			if (!events[i].isUsedForClustering) continue;
			nOpen[events[i].cluster]+=events[i].nOpen;
			nClosed[events[i].cluster]+=events[i].nClosed;
		}
		for (i=0; i<nClusters[bestIteration]; i++) {
			if (nKmeansIterations>0 && !keep[i]) continue;
			Factors.lastSplit++;
			Factors.splits[Factors.lastSplit].clear();
			Factors.splits[Factors.lastSplit].position=(int)clusterMeans[bestIteration][i];
			Factors.splits[Factors.lastSplit].nOpen=nOpen[i];
			Factors.splits[Factors.lastSplit].nClosed=nClosed[i];
			// $nUnknown$ is left to zero, since no event has unknown direction.
if (IO.SHOW_STD_ERR) IO.printErr("Events.cluster> adds split "+Factors.splits[Factors.lastSplit].position);

		}
		return bestIteration;
	}


	/**
	 * Let $\ell_i$ and $\ell_j$ be the first local-maximum leaf and the last local-
	 * maximum leaf of a run of local-maximum leaves detected by procedure
	 * $DensityEstimationTree.buildDensityEstimationTree$. This procedure assigns a
	 * cluster center to the center of mass of all events between the first event of
	 * $\ell_i$ and the last event of $\ell_j$, inclusive.
	 *
	 * If the run contains exactly two local maxima, the cluster center is either assigned
	 * to the midpoint of the leftmost local-minimum leaf between such local maxima, or to
	 * the midpoint between the last position of the first local maximum and the first
	 * position of the second local maximum, if no local-minimum leaf is found between the
	 * local maxima.
	 *
	 * The special treatment for runs with exactly two local maxima is motivated by the
	 * fact that two adjacent repeats could produce two, very close density peaks at their
	 * interface, separated by a very steep local minimum (e.g. if the two repeats never
	 * occur together elsewhere). Assigning the cluster mean based on the center of mass
	 * would move the mean closer to the repeat with higher frequency, which is not
	 * necessarily correct.
	 *
	 * Remark: the procedure assumes that procedure
	 * $DensityEstimationTree.markRunsOfLocalMaximumLeaves$ has already been executed.
	 *
	 * Remark: the output is stored in iteration zero of the clustering data structures.
	 */
	private static final void getRoughClustering() {
		boolean found;
		int i, j;
		int lastLeaf, first, last, nPoints;
		int firstLocalMaximum, lastLocalMaximum, lastLocalMinimum;
		double sum;
		Leaf[] leaves;

		leaves=DensityEstimationTree.leaves;
		lastLeaf=DensityEstimationTree.lastLeaf;
		firstLocalMaximum=0;
		lastLocalMaximum=DensityEstimationTree.lastLocalMaximum;
		lastLocalMinimum=DensityEstimationTree.lastLocalMinimum;
		nClusters[0]=0;
		for (i=0; i<=lastLocalMaximum; i++) {
			if (!leaves[i].marked) continue;

			// Computing mass and center of mass of the run
			first=leaves[firstLocalMaximum].firstPoint;
			last=leaves[i].lastPoint;
if (IO.SHOW_STD_ERR) IO.printErr("%%%%% COMPUTING THE CENTER OF MASS OF EVENTS ["+first+".."+last+"]");
			sum=0.0; nPoints=0;
			for (j=first; j<=last; j++) {
				events[j].cluster=nClusters[0];				
				sum+=events[j].position*events[j].getMass();
				nPoints+=events[j].getMass();
			}
			clusterSizes[0][nClusters[0]]=nPoints;

			// Computing cluster mean
			if (i==firstLocalMaximum+1) {
				found=false;
				for (j=lastLocalMaximum+1; j<=lastLocalMinimum; j++) {
					if ( leaves[j].firstPoint>=leaves[firstLocalMaximum].lastPoint &&
						 leaves[j].lastPoint<=leaves[i].firstPoint 
					   ) {
						clusterMeans[0][nClusters[0]]=(events[leaves[j].firstPoint].position+events[leaves[j].lastPoint].position)/2;
						found=true;
						break;
					}
				}
				if (!found) clusterMeans[0][nClusters[0]]=(events[leaves[firstLocalMaximum].lastPoint].position+events[leaves[i].firstPoint].position)/2;
			}
			else clusterMeans[0][nClusters[0]]=sum/nPoints;

if (IO.SHOW_STD_ERR) IO.printErr("center of mass="+clusterMeans[0][nClusters[0]]);

			nClusters[0]++;
			firstLocalMaximum=i+1;
		}
		
		// Marking all events outside the runs as not being used for clustering
		for (i=0; i<=lastEvent; i++) {
			if (events[i].cluster==-1) events[i].isUsedForClustering=false;
		}
	}


	/**
	 * Creates one cluster center per run of local-maximum leaves. The center of the
	 * cluster is set to a random event inside the run. A cluster center is not created if
	 * the mass of events in the run is less than $minMassRatio$ times the maximum mass of
	 * events in any run.
	 *
	 * @param iteration row of $clusterMeans$ and $nClusters$ to be used;
	 * @param maxMass maximum mass of events in a run of local-maximum leaves.
	 */
	private static final void initializeKmeans(int iteration, int maxMass) {
		int i, k;
		int lastLeaf, firstLocalMaximum, first, last;
		double mass;
		Leaf[] leaves;

		lastLeaf=DensityEstimationTree.lastLeaf;
		leaves=DensityEstimationTree.leaves;
		nClusters[iteration]=0;
		firstLocalMaximum=0;
		k=0;  // Pointer to events
		for (i=0; i<=lastLeaf; i++) {
			if (!leaves[i].isLocalMaximum) break;
			if (!leaves[i].marked) continue;
			first=leaves[firstLocalMaximum].firstPoint;
			last=leaves[i].lastPoint;
			mass=0.0;
			while (k<first) k++;
			while (k<=last) mass+=events[k++].getMass();
			if (mass<maxMass*minMassRatio) continue;
			clusterMeans[iteration][nClusters[iteration]++]=events[first+Math.random.nextInt(last-first+1)].position;
			firstLocalMaximum=i+1;
		}
	}


	/**
	 * Marks all events in $events[first..lastEvent]$ as not to be deleted.
	 */
	public static final void undelete(int first) {
		for (int i=first; i<=lastEvent; i++) events[i].deleted=false;
	}


	/**
	 * Sorts $events[first..lastEvent]$ by $order$.
	 */
	public static final void sort(int first, int order) {
		if (Event.order!=order) {
			Event.order=order;
			if (lastEvent>0) Arrays.sort(events,first,lastEvent+1);
			if (first!=0) Event.order=Event.UNSORTED;
		}
	}


	/**
	 * Sets $deleted=true$ for all events in $events[first..lastEvent]$ that belong to
	 * interval $readA[from..to]$, and increments $nDeletedEvents$ by the number of
	 * deleted events.
	 *
	 * Remark: the procedure assumes $Events.events[first..lastEvent]$ to be sorted by
	 * position.
	 */
	public static final void delete(int first, int from, int to) {
		int start, end;

		query.position=from;
		start=Arrays.binarySearch(events,first,lastEvent+1,query);
		if (start<0) start=-start-1;
		query.position=to;
		end=Arrays.binarySearch(events,first,lastEvent+1,query);
		if (end<0) {
			end=-end-1;
			end=Math.min(end,lastEvent);
		}
		for (int i=start; i<=end; i++) events[i].deleted=true;
		nDeletedEvents+=end-start+1;
	}


	/**
	 * Pushes out of visibility all events in $events[first..lastEvent]$ that are marked 
	 * as deleted.
	 */
	public static final void removeDeletedEvents(int first) {
		int i, j;
		Event tmp;

		j=first-1;
		for (i=first; i<=lastEvent; i++) {
			if (!events[i].deleted) {
				j++;
				if (j!=i) {
					tmp=events[j];
					events[j]=events[i];
					events[i]=tmp;
				}
			}
		}
		lastEvent=j;
	}
	
	
	public static final void ensureSpace_events(int nEvents) {
		final double GROWTH_RATE = 1.5;  // Arbitrary
		if (events.length>=nEvents) return;
		
		Event[] newEvents = new Event[(int)(nEvents*GROWTH_RATE)];
		System.arraycopy(events,0,newEvents,0,events.length);
		for (int i=events.length; i<newEvents.length; i++) newEvents[i] = new Event();
		events=newEvents;
	}

}