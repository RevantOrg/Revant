package de.mpi_cbg.revant.util;

import de.mpi_cbg.revant.util.Math;


public class Clusters {

	/**
	 * Uses k-means to group into $k$ clusters all points with $isUsedForClustering=true$.
	 * Assumes that $clusterMeans$ contains the initial position of the clusters.
	 */
	public static final void kmeans(int k, Point[] points, int lastPoint, double[] clusterMeans, int[] clusterSizes) {
		boolean assignmentsChanged;
		int c, p, minCluster, mass;
		double distance, minDistance;

		for (p=0; p<=lastPoint; p++) points[p].cluster=-1;
		while (true) {
			// Assigning points to clusters
			assignmentsChanged=false;
			for (p=0; p<=lastPoint; p++) {
				if (!points[p].isUsedForClustering) continue;
				minDistance=Math.POSITIVE_INFINITY; minCluster=-1;
				for (c=0; c<k; c++) {
					if (clusterMeans[c]==-1) continue;  // Empty clusters cannot be refilled
					distance=points[p].position-clusterMeans[c];
					if (distance<0) distance=-distance;
					if (distance<minDistance) {
						minDistance=distance;
						minCluster=c;
					}
				}
				if (points[p].cluster!=minCluster) {
					assignmentsChanged=true;
					points[p].cluster=minCluster;
				}
			}
			if (!assignmentsChanged) break;  // Convergence
			// Computing new cluster values
			Math.set(clusterMeans,k-1,0);
			Math.set(clusterSizes,k-1,0);
			for (p=0; p<=lastPoint; p++) {
				if (!points[p].isUsedForClustering) continue;
				mass=points[p].getMass();
				clusterMeans[points[p].cluster]+=points[p].position*mass;
				clusterSizes[points[p].cluster]+=mass;
			}
			for (c=0; c<k; c++) {
				if (clusterSizes[c]!=0) clusterMeans[c]/=clusterSizes[c];
				else clusterMeans[c]=-1;
			}
		}
	}


	/**
	 * @return the average silhouette over all marked points (the higher the better).
	 * Remark: this is not necessarily always positive, even though the clustering is
	 * performed by k-means. If we divided the sum of all intra-cluster distances by the
	 * size of the cluster, we would obtain an always positive measure under k-means.
	 * @param averageDistances temporary space with at least $lastPoint+1$ elements.
	 */
	public static final double getSilhouette(Point[] points, int lastPoint, int nClusters, int[] clusterSizes, double[] averageDistances) {
		int c, p, pPrime, totalMass;
		double minDistance, quality;

		quality=0; totalMass=0;
		for (p=0; p<=lastPoint; p++) {
			if (!points[p].isUsedForClustering) continue;
			Math.set(averageDistances,nClusters-1,0);
			for (pPrime=0; pPrime<=lastPoint; pPrime++) {
				if (!points[pPrime].isUsedForClustering) continue;
				if (pPrime==p) continue;
				averageDistances[points[pPrime].cluster]+=Math.abs(points[p].position-points[pPrime].position)*points[pPrime].getMass();
			}
			minDistance=Math.POSITIVE_INFINITY;
			for (c=0; c<nClusters; c++) {
				if (clusterSizes[c]==0) continue;
				if (c==points[p].cluster) {
					if (clusterSizes[c]>1) averageDistances[c]/=clusterSizes[c]-1;
					continue;
				}
				averageDistances[c]/=clusterSizes[c];
				if (averageDistances[c]<minDistance) minDistance=averageDistances[c];
			}
			if (clusterSizes[points[p].cluster]>1) {  // If there is no other point in the cluster of a point $p$, the silhouette of $p$ is not defined.
				quality+=points[p].getMass()*(minDistance-averageDistances[points[p].cluster])/Math.max(minDistance,averageDistances[points[p].cluster]);
				totalMass+=points[p].getMass();
			}
		}
		return quality/totalMass;
	}


	/**
	 * Keeps only clusters that are sufficiently similar to a Gaussian, i.e. that
	 * have small enough variance and large enough mass at distance at most one standard
	 * deviation from the mean. Clusters that do not conform to this model are very likely
	 * the union of at least two peaks, and cutting the read at their center of mass is
	 * not necessarily correct. Such clusters are discarded (although it would be better
	 * to try and recluster their points). All other clusters are flagged in bitvector
	 * $keep$.
	 *
	 * @param clusterVariances,clusterMassInStandardDeviation Temporary space;
	 * @return the number of clusters that are kept.
	 */
	public static final int filterClusters(Point[] points, int lastPoint, int nClusters, int[] clusterSizes, double[] clusterMeans, double[] clusterVariances, double[][] clusterMassInStandardDeviation, boolean[] keep) {
		int c, p, nKept, minCluster;
		double distance, minDistance;
		double threshold1, threshold2;

		// Assigning again every point to its closest cluster center. This is necessary,
		// since points are overwritten by every k-means iteration.
		for (p=0; p<=lastPoint; p++) {
			if (!points[p].isUsedForClustering) {
				points[p].cluster=-1;
				continue;
			}
			minCluster=-1; minDistance=Math.POSITIVE_INFINITY;
			for (c=0; c<nClusters; c++) {
				distance=Math.abs(points[p].position-clusterMeans[c]);
				if (distance<minDistance) {
					minDistance=distance;
					minCluster=c;
				}
			}
			points[p].cluster=minCluster;
		}

		// Computing, for each cluster, the variance and the number of points at distance
		// from the mean at most equal to one and two standard deviations.
		Math.set(clusterVariances,nClusters-1,0.0);
		for (c=0; c<nClusters; c++) {
			clusterMassInStandardDeviation[c][0]=0.0;
			clusterMassInStandardDeviation[c][1]=0.0;
		}
		for (p=0; p<=lastPoint; p++) {
			if (!points[p].isUsedForClustering) continue;
			clusterVariances[points[p].cluster]+=(points[p].position-clusterMeans[points[p].cluster])*(points[p].position-clusterMeans[points[p].cluster])*points[p].getMass();
		}
		for (c=0; c<nClusters; c++) {
			if (clusterSizes[c]!=0) clusterVariances[c]=Math.sqrt(clusterVariances[c]/clusterSizes[c]);
		}
		for (p=0; p<=lastPoint; p++) {
			if (!points[p].isUsedForClustering) continue;
			c=points[p].cluster;
			if (Math.abs(points[p].position-clusterMeans[c])<=clusterVariances[c]) clusterMassInStandardDeviation[c][0]+=points[p].getMass();
			if (Math.abs(points[p].position-clusterMeans[c])<=2*clusterVariances[c]) clusterMassInStandardDeviation[c][1]+=points[p].getMass();
		}

		// Filtering
		threshold1=Math.floor(Math.FRACTION_IN_STANDARD_DEVIATIONS[0]*10)/10;  // We try to be a little more permissive than Gaussian
		threshold2=Math.floor(Math.FRACTION_IN_STANDARD_DEVIATIONS[1]*10)/10;  // We try to be a little more permissive than Gaussian
		for (c=0; c<nClusters; c++) keep[c]=false;
		nKept=0;
		for (c=0; c<nClusters; c++) {
			if (clusterSizes[c]==0) continue;
			if (Math.round(clusterMassInStandardDeviation[c][1]/clusterSizes[c],1)<threshold2) {
				if (IO.SHOW_STD_ERR) IO.printErr("Warning: fraction of mass in two standard deviations of cluster "+clusterMeans[c]+" is too small ("+(clusterMassInStandardDeviation[c][1]/clusterSizes[c])+" rather than "+threshold2+"). Cluster discarded.");
				continue;
			}
			if (Math.round(clusterMassInStandardDeviation[c][0]/clusterSizes[c],1)<threshold1) {
				if (IO.SHOW_STD_ERR) IO.printErr("Warning: fraction of mass in one standard deviation of cluster "+clusterMeans[c]+" is too small ("+(clusterMassInStandardDeviation[c][0]/clusterSizes[c])+" rather than "+threshold1+"). Cluster discarded.");
				continue;
			}
			keep[c]=true;
			nKept++;
		}
		return nKept;
	}

}