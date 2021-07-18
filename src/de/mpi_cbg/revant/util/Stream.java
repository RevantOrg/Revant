package de.mpi_cbg.revant.util;

import de.mpi_cbg.revant.util.Math;


/**
 * Expanding and contracting stack of integers
 */
public class Stream {

	private final int INTS_PER_REGION;
	private final int LOG2_INTS_PER_REGION;
	private final int MASK;

	private int[][] regions;
	private int topRegion, topPointer;  // Last used element in the stack
	private int nElements;  // Total number of elements in the stream


	/**
	 * @param intsPerRegion must be a power of two.
	 */
	public Stream(int intsPerRegion) {
		INTS_PER_REGION=intsPerRegion;
		LOG2_INTS_PER_REGION=Math.log2(intsPerRegion);
		final int THIRTYTWO_MINUS_LOG2_INTS_PER_REGION=32-LOG2_INTS_PER_REGION;
		MASK=0xFFFFFFFF>>>THIRTYTWO_MINUS_LOG2_INTS_PER_REGION;
		regions = new int[1][intsPerRegion];
		topRegion=0; topPointer=-1;
		nElements=0;
	}


	public void clear(boolean deallocate) {
		if (deallocate) regions = new int[1][INTS_PER_REGION];
		topRegion=0; topPointer=-1;
		nElements=0;
	}


	public void deallocate() {
		int nRegions = regions.length;
		for (int i=0; i<nRegions; i++) regions[i]=null;
		regions=null;
	}


	public final int nElements() {
		return nElements;
	}


	/**
	 * Appends $value$ to the stack, possibly expanding the stack. Expansion implies:
	 * (1) doubling the size of $regions$; (2) copying $regions.length$ pointers (using
	 * the internalized procedure $System.arraycopy$); (3) increasing the number of
	 * allocated bits by one region.
	 */
	public final void push(int value) {
		if (topPointer==INTS_PER_REGION-1) {
			int nRegions = regions.length;
			if (topRegion==nRegions-1) {
				int[][] newRegions = new int[nRegions<<1][0];
				System.arraycopy(regions,0,newRegions,0,nRegions);
				regions=newRegions;
			}
			topRegion++;
			regions[topRegion] = new int[INTS_PER_REGION];
			topPointer=0;
		}
		else topPointer++;
		regions[topRegion][topPointer]=value;
		nElements++;
	}


	/**
	 * Removes the last element from the top of the stack, possibly contracting the stack.
	 * The last unused region (if any) is immediately set to NULL, but not necessarily
	 * immediately deallocated by the garbage collector.
	 */
	public final void pop() {
		nElements--;
		if (topPointer==0) {
			regions[topRegion]=null;
			topRegion--;
			topPointer=INTS_PER_REGION-1;
		}
		else topPointer--;
	}


	/**
	 * Reads the $i$th element in the stack. The procedure assumes that $i$ is a valid
	 * index: no explicit check is performed.
	 */
	public final int getElementAt(int i) {
		int pointer = i&MASK;
		i>>>=LOG2_INTS_PER_REGION;
		return regions[i][pointer];
	}

}