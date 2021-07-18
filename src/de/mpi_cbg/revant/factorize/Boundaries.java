package de.mpi_cbg.revant.factorize;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Histograms;


public class Boundaries {
	/**
	 * Upper bound on the part of a factor to be used in boundary refinement
	 */
	public static int boundaryRefinementLength;
	
	
	/**
	 * Moves the boundary of a factor to the leftmost (if $direction=-1$) or to the
	 * rightmost (if $direction=1$) position of maximum absolute second derivative in
	 * $histogram$ inside a neighborhood of the boundary (see procedures
	 * $getFirstPositionForBoundaryRefinement$, $getLastPositionForBoundaryRefinement$,
	 * $Histograms.getMaxSecondDerivative$).
	 * The procedure does not alter the boundary if this makes a factor shorter than
	 * $Factors.minFactorLength$.
	 *
	 * If $compareToPrevious=TRUE$, the boundary is moved only if the new value is smaller
	 * (if $direction=-1$) or larger (if $direction=1$) than the existing value. If
	 * $compareToPrevious=TRUE$ and $direction=0$, the previous value of the boundary is
	 * disregarded. If $direction=0$, the boundary is set to the center of mass of all
	 * positions with maximum absolute second derivative.
	 *
	 * @param start first position of the first factor in the sequence of factors;
	 * @param end last position of the last factor in the sequence of factors.
	 */
	public static final void refineBoundary(Factor factor, Factor previousFactor, Factor nextFactor, int start, int end, boolean isLeftBoundary, double[] histogram, int direction, boolean compareToPrevious, boolean isQuality) {
		final int STEP = IO.quantum;
		int firstPosition, newFirstPosition, lastPosition, newLastPosition, newBoundary;
		
		if (isLeftBoundary) {
			firstPosition=getFirstPositionForBoundaryRefinement(factor,previousFactor,true,histogram,isQuality);
			if (firstPosition==end+1) return;
			lastPosition=getLastPositionForBoundaryRefinement(factor,nextFactor,true,histogram,isQuality);
			if (lastPosition==start-1) return;
			if (lastPosition-firstPosition+1<IO.quantum) return;
			newBoundary=getMaxSecondDerivative(firstPosition,lastPosition,histogram,STEP,direction,isQuality);
			if ( compareToPrevious &&
			     ( (direction==-1 && newBoundary>=factor.firstPosition) ||
			       (direction==1 && newBoundary<=factor.firstPosition) )
			   ) return;
			// Refusing to shorten a factor below $Factors.minFactorLength$
			newFirstPosition=newBoundary+1;
			if (factor.lastPosition-newFirstPosition+1<Factors.minFactorLength) return;
			if (previousFactor!=null) {
				newLastPosition=newBoundary;
				if (newLastPosition-previousFactor.firstPosition+1<Factors.minFactorLength) return;
			}
			// Updating the boundary
			factor.firstPosition=newBoundary+1;
			factor.length=factor.lastPosition-factor.firstPosition+1;
			if (previousFactor!=null) {
				previousFactor.lastPosition=newBoundary;
				previousFactor.length=previousFactor.lastPosition-previousFactor.firstPosition+1;
			}
		}
		else {
			firstPosition=getFirstPositionForBoundaryRefinement(factor,previousFactor,false,histogram,isQuality);
			if (firstPosition==end+1) return;
			lastPosition=getLastPositionForBoundaryRefinement(factor,nextFactor,false,histogram,isQuality);
			if (lastPosition==start-1) return;
			if (lastPosition-firstPosition+1<IO.quantum) return;
			newBoundary=getMaxSecondDerivative(firstPosition,lastPosition,histogram,STEP,direction,isQuality);
			if ( compareToPrevious &&
			     ( (direction==-1 && newBoundary>=factor.lastPosition) ||
			       (direction==1 && newBoundary<=factor.lastPosition) )
			   ) return;
			// Refusing to shorten a factor below $Factors.minFactorLength$
			newLastPosition=newBoundary;
			if (newLastPosition-factor.firstPosition+1<Factors.minFactorLength) return;
			if (nextFactor!=null) {
				newFirstPosition=newBoundary+1;
				if (nextFactor.lastPosition-newFirstPosition+1<Factors.minFactorLength) return;
			}
			// Updating the boundary
			factor.lastPosition=newBoundary;
			factor.length=factor.lastPosition-factor.firstPosition+1;
			if (nextFactor!=null) {
				nextFactor.firstPosition=newBoundary+1;
				nextFactor.length=nextFactor.lastPosition-nextFactor.firstPosition+1;
			}
		}
	}
	
	
	/**
	 * Wraps $Histograms.getMaxSecondDerivative$ to work on $Reads.qualities$.
	 */
	private static final int getMaxSecondDerivative(int firstPosition, int lastPosition, double[] f, int step, int direction, boolean isQuality) {
		if (!isQuality) return Histograms.getMaxSecondDerivative(firstPosition,lastPosition,f,ReadA.readLength,step,direction);
		int i = Histograms.getMaxSecondDerivative(firstPosition/Reads.QUALITY_SPACING,lastPosition/Reads.QUALITY_SPACING,f,Reads.getQualityArrayLength(ReadA.id),Math.ceil(step,Reads.QUALITY_SPACING),direction);
		return Reads.QUALITY_SPACING*i+(Reads.QUALITY_SPACING>>1);  // A primitive way of smoothing
	}
	
	
	/**
	 * The refinement neighborhood of the left (respectively, right) boundary of a factor
	 * starts $boundaryRefinementLength$ positions before the first (respectively,
	 * last) position of the factor. If $histogram$ is approximately constant inside the
	 * interval between such starting position and the boundary, the refinement
	 * neighborhood starts at the boundary (otherwise, procedure $refineBoundary$ could
	 * move the boundary to the starting position just because it is the leftmost position
	 * in the neighborhood).
	 *
	 * @return the first position of $readA$ to be considered by $refineBoundary$ on
	 * $histogram$. Such position belongs to interval $[factors[0].firstPosition ..
	 * factors[lastFactor].lastPosition+1]$, so it's not necessarily a valid index in
	 * $ReadA$.
	 */
	private static final int getFirstPositionForBoundaryRefinement(Factor factor, Factor previousFactor, boolean isLeftBoundary, double[] histogram, boolean isQuality) {
		final int BIN_LENGTH = 1;
		int firstPosition;

		if (isLeftBoundary) {
			if (previousFactor==null) firstPosition=factor.firstPosition;
			else {
				firstPosition=previousFactor.lastPosition-(boundaryRefinementLength<previousFactor.length?boundaryRefinementLength:previousFactor.length>>1);
				if (Histograms.isConstant(histogram,isQuality?firstPosition/Reads.QUALITY_SPACING:firstPosition,isQuality?previousFactor.lastPosition/Reads.QUALITY_SPACING:previousFactor.lastPosition,BIN_LENGTH)) firstPosition=factor.firstPosition;
			}
		}
		else {
			firstPosition=factor.lastPosition-(boundaryRefinementLength<factor.length?boundaryRefinementLength:factor.length>>1);
			if (Histograms.isConstant(histogram,isQuality?firstPosition/Reads.QUALITY_SPACING:firstPosition,isQuality?factor.lastPosition/Reads.QUALITY_SPACING:factor.lastPosition,BIN_LENGTH)) firstPosition=factor.lastPosition+1;
		}
		return firstPosition;
	}
	
	
	/**
	 * The refinement neighborhood of the left (respectively, right) boundary of a factor
	 * ends $boundaryRefinementLength$ positions after the first (respectively, last)
	 * position of the factor. If $histogram$ is approximately constant inside the
	 * interval between the boundary and such ending position, the refinement neighborhood
	 * ends at the boundary (otherwise, procedure $refineBoundary$ could move the boundary
	 * to the ending position just because it is the rightmost position in the
	 * neighborhood).
	 *
	 * @return the last position of $readA$ to be considered by $refineBoundary$ on
	 * $histogram$. Such position belongs to interval $[factors[0].firstPosition-1 ..
	 * factors[lastFactor].lastPosition]$, so it's not necessarily a valid index in
	 * $ReadA$.
	 */
	private static final int getLastPositionForBoundaryRefinement(Factor factor, Factor nextFactor, boolean isLeftBoundary, double[] histogram, boolean isQuality) {
		final int BIN_LENGTH = 1;
		int lastPosition;

		if (isLeftBoundary) {
			lastPosition=factor.firstPosition+(boundaryRefinementLength<factor.length?boundaryRefinementLength:factor.length>>1);
			if (Histograms.isConstant(histogram,isQuality?factor.firstPosition/Reads.QUALITY_SPACING:factor.firstPosition,isQuality?lastPosition/Reads.QUALITY_SPACING:lastPosition,BIN_LENGTH)) lastPosition=factor.firstPosition-1;
		}
		else {
			if (nextFactor==null) lastPosition=factor.lastPosition;
			else {
				lastPosition=nextFactor.firstPosition+(boundaryRefinementLength<nextFactor.length?boundaryRefinementLength:nextFactor.length>>1);
				if (Histograms.isConstant(histogram,isQuality?nextFactor.firstPosition/Reads.QUALITY_SPACING:nextFactor.firstPosition,isQuality?lastPosition/Reads.QUALITY_SPACING:lastPosition,BIN_LENGTH)) lastPosition=factor.lastPosition;
			}
		}
		return lastPosition;
	}
	

}