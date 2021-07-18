package de.mpi_cbg.revant.util;

import java.io.*;
import de.mpi_cbg.revant.util.Math;


/**
 * Expanding stack of two-bit DNA characters.
 */
public class CharStream {

	private static final char[] ALPHABET = new char[] { 'a','c','g','t','A','C','G','T' };
	private static final int ALPHABET_LENGTH = 8;
	private static final int HALF_ALPHABET_LENGTH = 4;
	private static final int CHARS_PER_INT = 32;
	private static final int INTS_PER_REGION = 128;  // Must be a power of two
	private static final int LOG2_INTS_PER_REGION = Math.log2(INTS_PER_REGION);
	private static final int CHARS_PER_REGION = INTS_PER_REGION*CHARS_PER_INT;
	private static final int MASK = (~0)>>>(32-LOG2_INTS_PER_REGION);
	private static final int GROWTH_RATE = 1000;  // Arbitrary

	private long[][] regions;
	private int topRegion;  // Region of the last used integer in the stack
	private int topPointer;  // Inside $topRegion$.
	private int topChar;  // Inside $regions[topRegion][topPointer]$.
	private int nCharacters;  // Total number of characters in the stream


	/**
	 * @param initialCapacity in characters.
	 */
	public CharStream(int initialCapacity) {
		regions = new long[Math.max(1,Math.ceil(initialCapacity,CHARS_PER_REGION))][INTS_PER_REGION];
		topRegion=0; topPointer=0; topChar=-1;
		nCharacters=0;
	}


	public void clear(boolean deallocate) {
		if (deallocate) regions = new long[1][INTS_PER_REGION];
		topRegion=0; topPointer=0; topChar=-1;
		nCharacters=0;
	}


	public void deallocate() {
		int nRegions = regions.length;
		for (int i=0; i<nRegions; i++) regions[i]=null;
		regions=null;
	}


	public final int nCharacters() {
		return nCharacters;
	}
	
	
	/**
	 * Reads the integer encoding of the $j$-th character in the stack. This is an index 
	 * in $ALPHABET$, not a character. The procedure assumes that $j$ is a valid position:
	 * no explicit check is performed.
	 */
	public final int getCharacterAt(int j) {
		int i = j/CHARS_PER_INT;
		final long value = regions[i>>>LOG2_INTS_PER_REGION][i&MASK];
		final int shift = (j%CHARS_PER_INT)<<1;
		return (int)((value&((0x3L)<<shift))>>>shift);
	}

	
	/**
	 * Appends $source[first..last]$, in the forward ($direction=0$) or RC ($direction=1$)
	 * direction, possibly expanding the stack.
	 */
	public final void push(CharStream source, int first, int last, int direction) {
		ensureSize(nCharacters+(last-first+1));
		if (direction==0) {
			for (int i=first; i<=last; i++) push(source.getCharacterAt(i),true);
		}
		else {
			for (int i=last; i>=first; i--) push(3-source.getCharacterAt(i),true);
		}
	}
	

	/**
	 * Appends $character$ to the stack, possibly expanding the stack.
	 *
	 * @param bigEnough TRUE if there is enough space for one more character.
	 */
	public final void push(int character, boolean bigEnough) {
		if (topChar==CHARS_PER_INT-1) {
			if (topPointer==INTS_PER_REGION-1) {
				if (!bigEnough) ensureSize(nCharacters+GROWTH_RATE);
				topRegion++;
				topPointer=0;
			}
			else topPointer++;
			topChar=0;
		}
		else topChar++;
		final int shift = topChar<<1;
		regions[topRegion][topPointer]&=~((0x3L)<<shift);
		regions[topRegion][topPointer]|=((long)character)<<shift;
		nCharacters++;
	}

	
	/**
	 * @param newSize in characters.
	 */
	public final void ensureSize(int newSize) {
		final int nRegions = regions.length;
		if (newSize<=nRegions*CHARS_PER_REGION) return;
		int nNewRegions = Math.ceil(newSize,CHARS_PER_REGION);
		long[][] newRegions = new long[nNewRegions][0];
		System.arraycopy(regions,0,newRegions,0,nRegions);
		for (int i=nRegions; i<nNewRegions; i++) newRegions[i] = new long[INTS_PER_REGION];
		regions=newRegions;
	}

	
	/**
	 * Removes the last $nChars$ characters from the stack (no check is performed).
	 *
	 * @param contract TRUE: deallocates unused space.
	 */
	public final void pop(int nChars, boolean contract) {
		final int previousTopRegion = topRegion;
		nCharacters-=nChars;
		int i = (nCharacters-1)/CHARS_PER_INT;
		topPointer=i&MASK;
		topRegion=i>>>LOG2_INTS_PER_REGION;
		topChar=(nCharacters-1)%CHARS_PER_INT;
		if (contract && topRegion<previousTopRegion) {
			long[][] newRegions = new long[topRegion+1][0];
			System.arraycopy(regions,0,newRegions,0,topRegion+1);
			regions=newRegions;
		}
	}


	/**
	 * Removes the last character from the stack.
	 *
	 * @param contract TRUE: deallocates unused space.
	 */
	public final void pop(boolean contract) {
		final int previousTopRegion = topRegion;
		nCharacters--;
		if (topChar==0) {
			if (topPointer==0) {
				topRegion--;
				topPointer=INTS_PER_REGION-1;
			}
			else topPointer--;
			topChar=CHARS_PER_INT-1;
		}
		else topChar--;
		if (contract && topRegion<previousTopRegion) {
			long[][] newRegions = new long[topRegion+1][0];
			System.arraycopy(regions,0,newRegions,0,topRegion+1);
			regions=newRegions;
		}
	}
	
	
	/**
	 * @param direction TRUE=forward; FALSE=RC;
	 * @param prefixLength if >0, prints just a prefix of this length.
	 */
	public final void print(BufferedWriter bw, boolean direction, int prefixLength) throws IOException {
		if (direction) {
			final int LAST = ((prefixLength>0)&&(prefixLength<nCharacters))?prefixLength-1:nCharacters-1;
			for (int i=0; i<=LAST; i++) bw.write(ALPHABET[getCharacterAt(i)]);
		}
		else {
			final int LAST = ((prefixLength>0)&&(prefixLength<nCharacters))?nCharacters-prefixLength:0;
			for (int i=nCharacters-1; i>=LAST; i--) bw.write(ALPHABET[3-getCharacterAt(i)]);
		}
	}
	
	
	/**
	 * Loads $nCharacters$ from the current position of $br$. Characters not in $ALPHABET$ 
	 * are not loaded. Capitalization in the input is ignored.
	 */
	public final void load(int nCharacters, BufferedReader br) throws IOException {
		int i, j, value;
		
		ensureSize(nCharacters);
		clear(false);
		for (i=0; i<nCharacters; i++) {
			value=br.read();
			for (j=0; j<ALPHABET_LENGTH; j++) {
				if (ALPHABET[j]==value) {
					push(j%HALF_ALPHABET_LENGTH,false);
					break;
				}
			}
		}		
	}
	
	
	public void toString(StringBuilder sb) {
		final int LAST = nCharacters-1;
		for (int i=0; i<=LAST; i++) sb.append(ALPHABET[getCharacterAt(i)]);
	}

}