package de.mpi_cbg.revant.apps;

import java.io.*;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Intervals;


/**
 * Basic tools for handling an alphabet of repeat units.
 */
public class RepeatAlphabet {
	public static final int UNIQUE = -1;
	public static final char SEPARATOR_MINOR = ',';
	public static final char SEPARATOR_MAJOR = ':';
	
	/**
	 * Length of every repeat
	 */
	public static int[] repeatLengths;
	public static boolean[] isPeriodic;
	
	/**
	 * The recoded alphabet
	 */
	public static Character[] alphabet;
	public static int lastUnique, lastPeriodic, lastAlphabet;
	public static int maxOpenLength_unique;
	public static long[] alphabetCount;
	
	/**
	 * All alignments of a given readA
	 */
	private static AlignmentRow[] alignments;
	private static int lastAlignment;
	
	/**
	 * The recoded sequence of a given readA
	 */
	public static Block[] sequence;
	public static int lastInSequence;
	
	/**
	 * Length of the recoded sequence of every read
	 */
	public static int[] sequenceLengths;
	
	/**
	 * Sorted lists of read IDs: with a blue interval; fully non-repetitive; fully 
	 * contained in a single repeat block; translated into more than one block.
	 */
	private static int[] blueIntervals_reads;
	private static int[] fullyUnique, fullyContained, translated;
	
	/**
	 * For every read in $blueIntervals_reads$: tuples (position,length,nHaplotypes).
	 *
	 * Remark: could be an array of bytes, since there are fewer than 128 blocks per read 
	 * in practice: we keep it as an array of 32-bit integers to be fully generic.
	 */
	private static int[][] blueIntervals;
	private static int blueIntervals_last;
	
	/**
	 * For every read in $translated$: block boundaries; flags marking non-repetitive
	 * blocks; full translation; tandem intervals.
	 */ 
	private static int[][] boundaries_all;
	private static byte[][] isBlockUnique_all;
	private static int[][][] translation_all;
	private static int[][] tandems;
	private static int[] lastTandem;
	
	/**
	 * Data structures for fixing inaccurate periodic endpoints.
	 */
	private static Spacer[] spacers;
	public static int lastSpacer, nRigidSpacers, nBridgingSpacers, nInactiveSpacers;
	private static double[][] spacerNeighbors;
	private static int[] lastSpacerNeighbor;
	
	/**
	 * Temporary space used by procedure $recodeRead()$.
	 */
	private static int[] periodicIntervals, points;
	private static int lastPoint;
	private static Block newBlock;
	private static Character tmpCharacter;
	
	/**
	 * Temporary space encoding an undirected graph and its connected components (see
	 * procedure $compactAlphabet()$).
	 */
	private static int[][] neighbors;
	private static int[] lastNeighbor;
	private static int[] stack, stack2, connectedComponent, connectedComponentSize;
	private static int nComponents;
	private static Character[] mergedCharacters;
	
	/**
	 * Global variables used by $fixPeriodicEndpoints_updateTranslation*()$ and 
	 * $fixPeriodicEndpoints_collectCharacterInstances()$.
	 */
	private static int lastLeft, blockCursor, currentBoundary, nBoundariesWritten;
	
	/**
	 * Temporary space
	 */
	private static String[][] blocks;
	private static int[] lastInBlock, lastInBlock_int;
	public static int[] boundaries;
	private static int[][] intBlocks;
	private static boolean[] isBlockUnique;
	private static Pair[] tandemIntervals;
	private static Kmer[] kmerPool;
	private static int lastKmerPool;
	private static boolean[] tmpBoolean;
	private static Character[] leftCharacters, rightCharacters;
	
	
	/**
	 * Like $Reads.loadReadLengths()$, but for repeats.
	 */
	public static final int loadRepeatLengths(String repeatLengthsFile, int nRepeats) throws IOException {
		int i;
		int length, maxLength;
		BufferedReader br;

		repeatLengths = new int[nRepeats];
		br = new BufferedReader(new FileReader(repeatLengthsFile));
		maxLength=0;
		for (i=0; i<nRepeats; i++) {
			length=Integer.parseInt(br.readLine());
			repeatLengths[i]=length;
			maxLength=Math.max(maxLength,length);
		}
		br.close();
		return maxLength;
	}
	
	
	public static final void loadIsPeriodic(String isPeriodicFile, int nRepeats) throws IOException {
		isPeriodic = new boolean[nRepeats];
		BufferedReader br = new BufferedReader(new FileReader(isPeriodicFile));
		for (int i=0; i<nRepeats; i++) isPeriodic[i]=Integer.parseInt(br.readLine())==1;
		br.close();
	}
	
	
	
	
	// ------------ ALPHABET CONSTRUCTION AND READ TRANSLATION PROCEDURES ----------------
	
	/**
	 * Recodes every read, and collects in $outputFile$ every character (not necessarily 
	 * distinct) with discriminative power.
	 *
	 * Remark: the aligner might not find every occurrence of every repeat in every read.
	 * Such missing occurrences will become non-repetitive sequence, which will induce
	 * several spurious alignments downstream.
	 *
	 * Remark: the error rate in read-repeat alignments should not be too high compared to
	 * the one used in read-read alignments, otherwise the following could happen. Assume
	 * that there are many variants of the same repeat, such that each variant occurs in
	 * a single locus of the genome, and the variants can be discriminated by read-read 
	 * alignments. Every such variant would get tagged with the same repeat ID, and any 
	 * alignment strictly contained inside such variants would be filtered out downstream 
	 * because assumed to be repeat-induced -- even though such alignments come from the 
	 * same locus of the genome. In practice this might not be a problem, because of the 
	 * context around the variants in the genome.
	 *
	 * @param alignmentsFile alignments between reads (readA) and repeats (readB);
	 * @param maxError alignments with error rate greater than this are discarded;
	 * @param minAlignmentLength alignments shorter than this on either side are 
	 * discarded;
	 * @param uniqueFile the procedure writes just $maxOpenLength_unique$ there.
	 */
	public static final void collectCharacterInstances(String alignmentsFile, double maxError, int minAlignmentLength, int distanceThreshold, int lengthThreshold, String outputFile, String uniqueFile) throws IOException {
		final int ALPHABET_CAPACITY = 100000;  // Arbitrary
		final int ALIGNMENTS_CAPACITY = 100000;  // Arbitrary
		final int SEQUENCE_CAPACITY = 1000000;  // Arbitrary
		final int MAX_SEQUENCE_LENGTH = 1000;  // Arbitrary
		int i;
		int row, readA, previousReadA;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		// Allocating memory
		if (alignments==null || alignments.length<ALIGNMENTS_CAPACITY) alignments = new AlignmentRow[ALIGNMENTS_CAPACITY];
		for (i=0; i<alignments.length; i++) {	
			if (alignments[i]==null) alignments[i] = new AlignmentRow();
		}
		if (sequence==null || sequence.length<SEQUENCE_CAPACITY) sequence = new Block[SEQUENCE_CAPACITY];
		for (i=0; i<sequence.length; i++) {
			if (sequence[i]==null) sequence[i] = new Block();
		}
		if (points==null || points.length<(ALIGNMENTS_CAPACITY)<<1) points = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (periodicIntervals==null || periodicIntervals.length<(ALIGNMENTS_CAPACITY)<<1) periodicIntervals = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (sequenceLengths==null || sequenceLengths.length<MAX_SEQUENCE_LENGTH) sequenceLengths = new int[MAX_SEQUENCE_LENGTH];
		Math.set(sequenceLengths,MAX_SEQUENCE_LENGTH-1,0);
		if (tmpCharacter==null) tmpCharacter = new Character();
		if (stack==null || stack.length<ALIGNMENTS_CAPACITY) stack = new int[ALIGNMENTS_CAPACITY];
		
		// Collecting characters from all recoded reads
		maxOpenLength_unique=0;
		bw = new BufferedWriter(new FileWriter(outputFile));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); 
		previousReadA=-1; lastAlignment=-1; row=0;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			if ( (2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2)>maxError ||
				 Math.max(Alignments.endA-Alignments.startA,Alignments.endB-Alignments.startB)+1<minAlignmentLength
			   ) {
				str=br.readLine(); row++;
				continue;
			}
			readA=Alignments.readA-1;
			if (previousReadA==-1 || readA!=previousReadA) {
				if (previousReadA!=-1) {
					cleanAlignments(distanceThreshold);
					if (lastAlignment!=-1) {
						recodeRead(distanceThreshold);
						if (lastInSequence>=0) addCharacterInstances(bw);
						sequenceLengths[lastInSequence+1<MAX_SEQUENCE_LENGTH?lastInSequence+1:MAX_SEQUENCE_LENGTH-1]++;
					}
					else sequenceLengths[0]++;
				}
				previousReadA=readA; lastAlignment=0;
				alignments[0].set(readA,Math.max(Alignments.startA,0),Math.min(Alignments.endA,Reads.getReadLength(readA)-1),Alignments.readB-1,Math.max(Alignments.startB,0),Math.min(Alignments.endB,repeatLengths[Alignments.readB-1]-1),Alignments.orientation,Alignments.diffs);
			}
			else {
				lastAlignment++;
				if (lastAlignment==alignments.length) {
					AlignmentRow[] newAlignments = new AlignmentRow[alignments.length<<1];
					System.arraycopy(alignments,0,newAlignments,0,alignments.length);
					for (i=alignments.length; i<newAlignments.length; i++) newAlignments[i] = new AlignmentRow();
					alignments=newAlignments;
				}
				alignments[lastAlignment].set(readA,Math.max(Alignments.startA,0),Math.min(Alignments.endA,Reads.getReadLength(readA)-1),Alignments.readB-1,Math.max(Alignments.startB,0),Math.min(Alignments.endB,repeatLengths[Alignments.readB-1]-1),Alignments.orientation,Alignments.diffs);
			}
			str=br.readLine(); row++;
		}
		br.close();
		if (previousReadA!=-1) {
			cleanAlignments(distanceThreshold);
			if (lastAlignment!=-1) {
				recodeRead(distanceThreshold);
				if (lastInSequence>=0) addCharacterInstances(bw);
				sequenceLengths[lastInSequence+1<MAX_SEQUENCE_LENGTH?lastInSequence+1:MAX_SEQUENCE_LENGTH-1]++;
			}
			else sequenceLengths[0]++;
		}
		bw.close();
		bw = new BufferedWriter(new FileWriter(uniqueFile));
		bw.write(maxOpenLength_unique+"\n");
		bw.close();
	}
	
	
	/**
	 * Stores in global variable $sequence$ a recoding of the read based on repeat 
	 * alignments. The procedure collects the endpoints of every maximal range of 
	 * overlapping periodic alignments, as well as the endpoints of every non-periodic 
	 * alignment that is not contained in a periodic range (non-periodic alignments might
	 * straddle one another and periodic ranges). Points are clustered, and the resulting
	 * blocks of the read become elements of $sequence$.
	 *
	 * Remark: the procedure assumes that reads do not contain long low-quality regions.
	 *
	 * Remark: endpoints are clustered by connected components based on distance 
	 * $distanceThreshold$. If two repeats overlap by at most this much, their overlapping 
	 * endpoints are merged into one, so the procedure creates blocks of length greater
	 * than $distanceThreshold$.
	 *
	 * Remark: assume that two consecutive real boundaries between (possibly overlapping)
	 * repeats are such that the set of all endpoints between them forms a dense region. 
	 * The procedure might discard all the points in such a dense region, and collapse the 
	 * two repeats into a single block. The correct way of solving this problem would be
	 * to perform density estimation on the distribution of endpoints and to select peaks.
	 * We don't do this for simplicity.
	 *
	 * Remark: it might happen that some alignments are not used for building blocks, and
	 * that the region of the read they cover might become part of a non-repetitive block.
	 * One might think of recording how such a non-repetitive block is covered by
	 * alignments (e.g. fully covered, partially covered, prefix, suffix, substring) since
	 * it might be a useful feature for discrimination. We don't pursue this for 
	 * simplicity.
	 *
	 * Remark: the characters in blocks of $sequence$ are not objects in $alphabet$, and 
	 * the procedure does not look them up in $alphabet$: this is left to the caller.
	 *
	 * Remark: the procedure uses global arrays $alignments,periodicIntervals,points,
	 * sequence,stack$. $stack$ is assumed to be non-null.
	 */
	public static final void recodeRead(int distanceThreshold) {
		final int CLUSTERING_DISTANCE = distanceThreshold;
		final int MAX_DENSE_LENGTH = (CLUSTERING_DISTANCE)<<1;  // Arbitrary
		int i, j, k;
		int lastPeriodicInterval, currentStart, currentEnd, first, firstZero, lastStack;
		int firstJForNextI, inPeriodic;
		int startA, endA, newLength, component, nComponents;
		final int lengthA = Reads.getReadLength(alignments[0].readA);
		
		if (periodicIntervals.length<(lastAlignment+1)<<1) periodicIntervals = new int[(lastAlignment+1)<<1];
		if (points.length<(lastAlignment+1)<<1) points = new int[(lastAlignment+1)<<1];
		AlignmentRow.order=AlignmentRow.ORDER_STARTA;
		if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
		if (stack.length<lastAlignment+1) stack = new int[lastAlignment+1];
		
		// Building maximal periodic intervals
		lastPeriodicInterval=-1; currentStart=-1; currentEnd=-1; lastStack=-1;
		for (i=0; i<=lastAlignment; i++) {
			if (!isPeriodic[alignments[i].readB]) continue;
			if (currentStart==-1) {
				currentStart=alignments[i].startA; currentEnd=alignments[i].endA;
				lastStack=0; stack[0]=alignments[i].readB;
				continue;
			}
			if (alignments[i].startA<currentEnd-distanceThreshold) {
				currentEnd=Math.max(currentEnd,alignments[i].endA);
				stack[++lastStack]=alignments[i].readB;
				continue;
			}
			if (lastStack>0) {
				Arrays.sort(stack,0,lastStack+1);
				j=0;
				for (k=1; k<=lastStack; k++) {
					if (stack[k]!=stack[j]) stack[++j]=stack[k];
				}
				lastStack=j;
			}
			if (alignments[i].startA>currentEnd+distanceThreshold || Arrays.binarySearch(stack,0,lastStack+1,alignments[i].readB)<0) {
				periodicIntervals[++lastPeriodicInterval]=currentStart;
				periodicIntervals[++lastPeriodicInterval]=currentEnd;
				currentStart=alignments[i].startA; currentEnd=alignments[i].endA; 
				lastStack=0; stack[0]=alignments[i].readB;
			}
			else currentEnd=Math.max(currentEnd,alignments[i].endA);
		}
		if (currentStart!=-1) {
			periodicIntervals[++lastPeriodicInterval]=currentStart;
			periodicIntervals[++lastPeriodicInterval]=currentEnd;
		}
		
		// Collecting readA endpoints of all alignments
		lastPoint=-1;
		for (i=0; i<=lastAlignment; i++) {
			points[++lastPoint]=alignments[i].startA;
			points[++lastPoint]=alignments[i].endA;
		}
		if (lastPoint>0) Arrays.sort(points,0,lastPoint+1);
		initializeGraph(lastPoint+1);		
		
		// Marking points: (0) outside periodic intervals; (1) inside a periodic interval
		// and close to its endpoints; (-1) inside a periodic interval and far from its
		// endpoints. $connectedComponent$ is used as temporary space.
		Math.set(connectedComponent,lastPoint,0);
		i=0; j=0; firstJForNextI=-1;
		while (i<=lastPoint) {
			if (j>lastPeriodicInterval || periodicIntervals[j]-distanceThreshold>points[i]) {
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (periodicIntervals[j+1]+distanceThreshold<points[i]) {
				j+=2;
				continue;
			}
			if (firstJForNextI==-1 && i<lastPoint && periodicIntervals[j+1]+distanceThreshold>=points[i+1]) firstJForNextI=j;
			if (points[i]>periodicIntervals[j]+distanceThreshold && points[i]<periodicIntervals[j+1]-distanceThreshold) {
				if (connectedComponent[i]!=1) connectedComponent[i]=-1;
			}
			else connectedComponent[i]=1;
			j+=2;
		}
		j=-1;
		for (i=0; i<=lastPoint; i++) {
			if (connectedComponent[i]==-1) continue;
			j++;
			points[j]=points[i]; connectedComponent[j]=connectedComponent[i];
		}
		lastPoint=j;
		
		// Removing points outside periodic intervals if they are in a dense region
		firstZero=-1;
		for (i=0; i<=lastPoint; i++) {
			if (connectedComponent[i]==1) {
				if (firstZero!=-1) {
					first=firstZero;
					for (j=firstZero+1; j<i; j++) {
						if (points[j]-points[j-1]>CLUSTERING_DISTANCE) {
							if (points[j-1]-points[first]+1>MAX_DENSE_LENGTH) {
								for (k=first; k<=j-1; k++) connectedComponent[k]=-1;
							}
							first=j;
						}
					}
					if (points[i-1]-points[first]+1>MAX_DENSE_LENGTH) {
						for (k=first; k<=i-1; k++) connectedComponent[k]=-1;
					}
					firstZero=-1;
				}
			}
			else if (firstZero==-1) firstZero=i;
		}
		if (firstZero!=-1) {
			first=firstZero;
			for (j=firstZero+1; j<=lastPoint; j++) {
				if (points[j]-points[j-1]>CLUSTERING_DISTANCE) {
					if (points[j-1]-points[first]+1>MAX_DENSE_LENGTH) {
						for (k=first; k<=j-1; k++) connectedComponent[k]=-1;
					}
					first=j;
				}
			}
			if (points[lastPoint]-points[first]+1>MAX_DENSE_LENGTH) {
				for (k=first; k<=lastPoint; k++) connectedComponent[k]=-1;
			}
		}
		j=-1;
		for (i=0; i<=lastPoint; i++) {
			if (connectedComponent[i]==-1) continue;
			points[++j]=points[i];
		}
		lastPoint=j;
		
		// Clustering all surviving points
		for (i=0; i<lastPoint; i++) {
			for (j=i+1; j<=lastPoint; j++) {
				if (points[j]>points[i]+CLUSTERING_DISTANCE) break;
				addEdge(i,j); addEdge(j,i);
			}
		}
		nComponents=getConnectedComponent(lastPoint+1);
		if (stack.length<(nComponents<<1)) stack = new int[nComponents<<1];
		for (i=0; i<(nComponents<<1); i++) stack[i]=0;
		for (i=0; i<=lastPoint; i++) {
			component=connectedComponent[i]<<1;
			stack[component]+=points[i];
			stack[component+1]++;
		}
		for (i=0; i<nComponents; i++) points[i]=stack[i<<1]/stack[(i<<1)+1];
		lastPoint=nComponents-1;
		
		// Creating the sequence
		if (newBlock==null) newBlock = new Block();
		lastInSequence=-1;
		// First non-repetitive block (if any).
		i=points[0];
		if (i>distanceThreshold) {
			lastInSequence=0;
			sequence[0].setUnique(i,true,i>=lengthA-distanceThreshold);
		}
		// Middle blocks
		if (stack==null || stack.length<lastAlignment+1) stack = new int[lastAlignment+1];
		i=1; j=0; firstJForNextI=-1;
		while (i<=lastPoint) {
			startA=points[i-1]; endA=points[i];
			if (j>lastAlignment || alignments[j].startA>=endA) {
				lastInSequence++;
				if (newBlock.lastCharacter==-1) newBlock.setUnique(endA-startA+1,startA<=distanceThreshold,endA>=lengthA-distanceThreshold);
				else {
					if (stack.length<newBlock.lastCharacter+1) stack = new int[newBlock.lastCharacter+1];
					newBlock.cleanCharacters(endA-startA+1,lastInSequence==0||(i==lastPoint&&lengthA-points[lastPoint]<=distanceThreshold),stack);
				}
				if (lastInSequence==sequence.length) {
					Block[] newSequence = new Block[sequence.length<<1];
					System.arraycopy(sequence,0,newSequence,0,sequence.length);
					for (k=sequence.length; k<newSequence.length; k++) newSequence[k] = new Block();
					sequence=newSequence;
				}
				sequence[lastInSequence].copyFrom(newBlock);
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				newBlock.reset();
				continue;
			}
			else if (alignments[j].endA<=startA) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastPoint && alignments[j].endA>endA) firstJForNextI=j;
			if (Intervals.areApproximatelyIdentical(alignments[j].startA,alignments[j].endA,startA,endA)) newBlock.addCharacter(alignments[j],distanceThreshold,alignments[j].startA,alignments[j].endA,tmpCharacter);
			else if ( (isPeriodic[alignments[j].readB] && Intervals.isApproximatelyContained(alignments[j].startA,alignments[j].endA,startA,endA)) ||
				      (!isPeriodic[alignments[j].readB] && Intervals.isApproximatelyContained(startA,endA,alignments[j].startA,alignments[j].endA))
				    ) newBlock.addCharacter(alignments[j],distanceThreshold,startA,endA,tmpCharacter);
			j++;
		}
		// Last non-repetitive block (if any).
		i=points[lastPoint];
		if (lengthA-i>distanceThreshold) {
			lastInSequence++;
			if (lastInSequence==sequence.length) {
				Block[] newSequence = new Block[sequence.length<<1];
				System.arraycopy(sequence,0,newSequence,0,sequence.length);
				for (k=sequence.length; k<newSequence.length; k++) newSequence[k] = new Block();
				sequence=newSequence;
			}
			sequence[lastInSequence].setUnique(lengthA-i,i<=distanceThreshold,true);
		}
	}
	
	
	/**
	 * Appends to $bw$ the characters of $sequence$.
	 */
	private static final void addCharacterInstances(BufferedWriter bw) throws IOException {	
		int i, j;
		
		for (i=0; i<=lastInSequence; i++) {
			if (sequence[i].isUnique()) {
				if (!sequence[i].characters[0].openStart && !sequence[i].characters[0].openEnd) {
					bw.write(sequence[i].characters[0].toString());
					bw.newLine();
				}
				else maxOpenLength_unique=Math.max(maxOpenLength_unique,sequence[i].getLength());
			}
			else {
				for (j=0; j<=sequence[i].lastCharacter; j++) {
					bw.write(sequence[i].characters[j].toString());
					bw.newLine();
				}
			}
		}
	}
	
	
	/**
	 * Quantizes the endpoints of every repeat, and removes character instances that are
	 * less specific than others (this is done just to reduce the size of the alphabet, to
	 * speed up the procedures downstream).
	 *
	 * Remark: the procedure compacts only characters in the same orientation.
	 *
	 * Remark: the procedure does not explicitly remove characters with no discriminative
	 * power, since it assumes they have not been added to $alphabet$ in ther first place.
	 * 
	 * Remark: assume that the repeat database contains a repeat that occurs in full just 
	 * once in the genome, and that is longer than every read. At the end of the procedure 
	 * this corresponds to an open prefix character and an open suffix character. A 
	 * similar observation holds also for a repeat with a long substring that occurs just 
	 * once in the genome, and for a long short-period repeat that occurs just once in the
	 * genome.
	 */
	public static final void compactInstances() {
		final int QUANTUM = IO.quantum;
		int i, j, k;
		int first;
		Character tmpChar;
		
		System.err.println("Quantizing characters... ");
		for (i=0; i<=lastAlphabet; i++) alphabet[i].quantize(QUANTUM);
		if (lastUnique>0) Arrays.sort(alphabet,0,lastUnique+1);
		if (lastPeriodic>lastUnique+1) Arrays.sort(alphabet,lastUnique+1,lastPeriodic+1);
		if (lastAlphabet>lastPeriodic+1) Arrays.sort(alphabet,lastPeriodic+1,lastAlphabet+1);
		j=0;
		for (i=1; i<=lastAlphabet; i++) {
			if ( alphabet[i].repeat!=alphabet[j].repeat || alphabet[i].orientation!=alphabet[j].orientation || 
				 alphabet[i].start!=alphabet[j].start || alphabet[i].end!=alphabet[j].end || alphabet[i].length!=alphabet[j].length ||
			     alphabet[i].openStart!=alphabet[j].openStart || alphabet[i].openEnd!=alphabet[j].openEnd
			   ) {
				j++;
				tmpChar=alphabet[j];
				alphabet[j]=alphabet[i];
				alphabet[i]=tmpChar;
			}
		}
		lastAlphabet=j;
		lastUnique=-1; i=0;
		while (i<=lastAlphabet && alphabet[i].repeat==UNIQUE) lastUnique=i++;
		lastPeriodic=lastUnique;
		while (i<=lastAlphabet && alphabet[i].start==-1) lastPeriodic=i++;
		System.err.println("DONE");
		
		System.err.println("Discarding repetitive characters that are implied by other repetitive characters... ("+(lastUnique+1)+" unique characters)");
		for (i=lastUnique+1; i<=lastAlphabet; i++) alphabet[i].flag=1;
		for (i=lastUnique+1; i<=lastAlphabet; i++) {
			if (i%1000==0) System.err.println("Processed "+i+" characters");
			for (j=i+1; j<=lastAlphabet; j++) {
				if ( alphabet[j].repeat!=alphabet[i].repeat || alphabet[j].orientation!=alphabet[i].orientation || 
					 alphabet[j].implies_tooFarAfter(alphabet[i])
				   ) break;
				if (alphabet[j].flag==0) continue;
				if (alphabet[i].implies(alphabet[j])) alphabet[j].flag=0;
			}
			for (j=i-1; j>=0; j--) {
				if (alphabet[j].repeat!=alphabet[i].repeat || alphabet[j].orientation!=alphabet[i].orientation) break;
				if (alphabet[j].flag==0) continue;
				if (alphabet[i].implies(alphabet[j])) alphabet[j].flag=0;
			}
		}
		j=lastUnique;
		for (i=lastUnique+1; i<=lastAlphabet; i++) {
			if (alphabet[i].flag==0) continue;
			j++;
			tmpChar=alphabet[j];
			alphabet[j]=alphabet[i];
			alphabet[i]=tmpChar;
		}
		lastAlphabet=j;
		lastPeriodic=lastUnique; i=lastUnique+1;
		while (i<=lastAlphabet && alphabet[i].start==-1) i++;
		lastPeriodic=i-1;
		System.err.println("DONE  "+(lastAlphabet+1)+" characters after filtering.");
	}
	
	
	/**
	 * Ensures that the reverse-complement of every character in the alphabet, is also in
	 * the alphabet.
	 *
	 * Remark: the reverse-complement of a character might be implied by a character that
	 * is already in the alphabet: in this case the new character is not added. Or, it
	 * might imply characters that are already in the alphabet: this case is not handled
	 * and it is left to the caller to run $compactInstances()$ again.
	 */
	public static final void closeAlphabetByRC() {
		final int QUANTUM = 1000;  // Arbitrary
		boolean openStart, openEnd, foundImplied;
		int i, j, k;
		int from, newFrom, repeat, start, end, length, found, nAdded, lastPeriodicNew;
		Character tmpCharacter;
		Character[] newAlphabet;
		
		// Unique
		newAlphabet = new Character[lastAlphabet+1];
		System.arraycopy(alphabet,0,newAlphabet,0,lastUnique+1);
		for (i=lastUnique+1; i<=lastAlphabet; i++) alphabet[i].flag=0;
		System.err.println("Adding reverse-complement characters... ");
		nAdded=0;
		
		// Periodic
		tmpCharacter = new Character();
		k=lastUnique; from=lastUnique+1; newFrom=from;
		for (i=lastUnique+1; i<=lastPeriodic; i++) {
			if (i%1000==0) System.err.println("Processed "+i+" characters");
			repeat=alphabet[i].repeat;
			openStart=alphabet[i].openStart; openEnd=alphabet[i].openEnd;		
			if (repeat!=alphabet[from].repeat) {
				if (k>newFrom) Arrays.sort(newAlphabet,newFrom,k+1);
				from=i; newFrom=k+1;
			}
			k++;
			if (k==newAlphabet.length) {
				Character[] newArray = new Character[newAlphabet.length+QUANTUM];
				System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
				newAlphabet=newArray;
			}
			newAlphabet[k]=alphabet[i];
			if (alphabet[i].flag==1) continue;
			tmpCharacter.copyFrom(alphabet[i]); tmpCharacter.reverseComplement();
			length=alphabet[i].length;
			if (alphabet[i].orientation) {
				found=-1; foundImplied=false;
				for (j=i+1; j<=lastPeriodic; j++) {
					if (alphabet[j].repeat!=repeat) break;
					if (!alphabet[j].orientation) {
						if (alphabet[j].length==length && alphabet[j].openStart==openStart && alphabet[j].openEnd==openEnd) {
							found=j;
							break;
						}
						if (alphabet[j].implies(tmpCharacter)) foundImplied=true;
					}
				}
				if (found==-1) {
					if (!foundImplied) {
						k++;
						if (k==newAlphabet.length) {
							Character[] newArray = new Character[newAlphabet.length+QUANTUM];
							System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
							newAlphabet=newArray;
						}
						newAlphabet[k] = new Character(repeat,false,-1,-1,length,openStart,openEnd);
						nAdded++;
					}
				}
				else alphabet[found].flag=1;
			}
			else {
				found=-1; foundImplied=false;
				for (j=i-1; j>lastUnique; j--) {
					if (alphabet[j].repeat!=repeat) break;
					if (alphabet[j].orientation) {
						if (alphabet[j].length==length && alphabet[j].openStart==openStart && alphabet[j].openEnd==openEnd) {
							found=j;
							break;
						}
						if (alphabet[j].implies(tmpCharacter)) foundImplied=true;
					}
				}
				if (found==-1 && !foundImplied) {
					k++;
					if (k==newAlphabet.length) {
						Character[] newArray = new Character[newAlphabet.length+QUANTUM];
						System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
						newAlphabet=newArray;
					}
					newAlphabet[k] = new Character(repeat,true,-1,-1,length,openStart,openEnd);
					nAdded++;
				}
			}
		}
		if (k>newFrom) Arrays.sort(newAlphabet,newFrom,k+1);
		lastPeriodicNew=k;
		
		// Nonperiodic
		from=lastPeriodic+1; newFrom=k+1;
		for (i=lastPeriodic+1; i<=lastAlphabet; i++) {
			if (i%1000==0) System.err.println("Processed "+i+" characters");
			repeat=alphabet[i].repeat;
			openStart=alphabet[i].openStart; openEnd=alphabet[i].openEnd;			
			if (repeat!=alphabet[from].repeat) {
				if (k>newFrom) Arrays.sort(newAlphabet,newFrom,k+1);
				from=i; newFrom=k+1;
			}			
			k++;
			if (k==newAlphabet.length) {
				Character[] newArray = new Character[newAlphabet.length+QUANTUM];
				System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
				newAlphabet=newArray;
			}
			newAlphabet[k]=alphabet[i];
			if (alphabet[i].flag==1) continue;
			tmpCharacter.copyFrom(alphabet[i]); tmpCharacter.reverseComplement();
			start=alphabet[i].start; end=alphabet[i].end;		
			if (alphabet[i].orientation) {
				found=-1; foundImplied=false;
				for (j=i+1; j<=lastAlphabet; j++) {
					if (alphabet[j].repeat!=repeat) break;
					if (!alphabet[j].orientation) {
						if (alphabet[j].start==start && alphabet[j].end==end && alphabet[j].openStart==openStart && alphabet[j].openEnd==openEnd) {
							found=j;
							break;
						}
						if (alphabet[j].implies(tmpCharacter)) foundImplied=true;
					}
				}				
				if (found==-1) {	
					if (!foundImplied) {
						k++;
						if (k==newAlphabet.length) {
							Character[] newArray = new Character[newAlphabet.length+QUANTUM];
							System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
							newAlphabet=newArray;
						}
						newAlphabet[k] = new Character(repeat,false,start,end,0,openStart,openEnd);
						nAdded++;
					}
				}
				else alphabet[found].flag=1;
			}
			else {
				found=-1; foundImplied=false;
				for (j=i-1; j>lastUnique; j--) {
					if (alphabet[j].repeat!=repeat) break;
					if (alphabet[j].orientation) {
						if (alphabet[j].start==start && alphabet[j].end==end && alphabet[j].openStart==openStart && alphabet[j].openEnd==openEnd) {
							found=j;
							break;
						}
						if (alphabet[j].implies(tmpCharacter)) foundImplied=true;
					}
				}
				if (found==-1 && !foundImplied) {
					k++;
					if (k==newAlphabet.length) {
						Character[] newArray = new Character[newAlphabet.length+QUANTUM];
						System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
						newAlphabet=newArray;
					}
					newAlphabet[k] = new Character(repeat,true,start,end,0,openStart,openEnd);
					nAdded++;
				}
			}
		}
		if (k>newFrom) Arrays.sort(newAlphabet,newFrom,k+1);
		lastAlphabet=k; lastPeriodic=lastPeriodicNew;
		alphabet=newAlphabet;
		System.err.println("DONE "+nAdded+" reverse-complement characters added ("+(lastAlphabet+1)+" total characters)");
	}
	
	
	public static final void serializeAlphabet(String path) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
		bw.write(lastUnique+(SEPARATOR_MINOR+"")+lastPeriodic+(SEPARATOR_MINOR+"")+lastAlphabet+(SEPARATOR_MINOR+"")+maxOpenLength_unique+"\n");
		for (int i=0; i<=lastAlphabet; i++) bw.write(alphabet[i].toString()+"\n");
		bw.close();
	}
	
	
	/**
	 * @param lastElement loads just the elements up to this (included): 0=lastUnique; 
	 * 1=lastPeriodic; 2=lastAlphabet. The $alphabet$ array is set to be big enough to 
	 * contain all characters, but the entries after $lastElement$ are undefined.
	 */
	public static final void deserializeAlphabet(String path, int lastElement) throws IOException {
		int i, p, q;
		int last;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();
		if (str==null) {
			lastUnique=-1; lastPeriodic=-1; lastAlphabet=-1;
			br.close();
			return;
		}
		p=str.indexOf(SEPARATOR_MINOR+"");
		lastUnique=Integer.parseInt(str.substring(0,p));
		p++; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
		lastPeriodic=Integer.parseInt(str.substring(p,q));
		p=q+1; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
		lastAlphabet=Integer.parseInt(str.substring(p,q));
		maxOpenLength_unique=Integer.parseInt(str.substring(q+1));
		alphabet = new Character[lastAlphabet+1];
		switch (lastElement) {
			case 0: last=lastUnique; break;
			case 1: last=lastPeriodic; break;
			default: last=lastAlphabet;
		}
		for (i=0; i<=last; i++) {
			alphabet[i] = new Character();
			alphabet[i].deserialize(br.readLine());
		}
		br.close();
	}
	
	
	/**
	 * Initializes the global data structures that describe a temporary graph.
	 */
	private static final void initializeGraph(int nNodes) {
		final int NEIGHBORS_CAPACITY = 10;  // Arbitrary
		int i;
		
		if (neighbors==null) neighbors = new int[nNodes][NEIGHBORS_CAPACITY];
		else if (neighbors.length<nNodes) {
			int[][] newNeighbors = new int[nNodes][0];
			for (i=0; i<neighbors.length; i++) newNeighbors[i]=neighbors[i];
			for (i=neighbors.length; i<newNeighbors.length; i++) newNeighbors[i] = new int[NEIGHBORS_CAPACITY];
			neighbors=newNeighbors;
		}
		if (lastNeighbor==null || lastNeighbor.length<nNodes) lastNeighbor = new int[nNodes];
		Math.set(lastNeighbor,nNodes-1,-1);
		if (connectedComponent==null || connectedComponent.length<nNodes) connectedComponent = new int[nNodes];
		if (stack==null || stack.length<nNodes) stack = new int[nNodes];
	}
	
	
	/**
	 * To the global data structures that describe a temporary graph.
	 */
	private static final void addEdge(int i, int j) {
		lastNeighbor[i]++;
		if (lastNeighbor[i]==neighbors[i].length) {
			int[] newArray = new int[neighbors[i].length<<1];
			System.arraycopy(neighbors[i],0,newArray,0,neighbors[i].length);
			neighbors[i]=newArray;
		}
		neighbors[i][lastNeighbor[i]]=j;
	}
	
	
	/**
	 * @return the number of connected components of the temporary graph. 
	 * Component IDs are marked in global array $connectedComponents$.
	 */
	private static final int getConnectedComponent(int nNodes) {
		int i, j;
		int last, idGenerator, top, currentNode, currentComponent, neighbor;
		int nComponents;
	
		// Computing components
		idGenerator=-1;
		Math.set(connectedComponent,0,nNodes-1,-1);
		for (i=0; i<nNodes; i++) {
			if (connectedComponent[i]!=-1) continue;
			currentComponent=++idGenerator;
			connectedComponent[i]=currentComponent;
			stack[0]=i; top=0;
			while (top>=0) {
				currentNode=stack[top--];
				last=lastNeighbor[currentNode];
				for (j=0; j<=last; j++) {
					neighbor=neighbors[currentNode][j];
					if (connectedComponent[neighbor]==-1) {
						connectedComponent[neighbor]=currentComponent;
						stack[++top]=neighbor;
					}
					else if (connectedComponent[neighbor]!=currentComponent) {
						System.err.println("getConnectedComponent> ERROR: different connected components "+currentComponent+" :: "+connectedComponent[neighbor]);
						System.exit(1);
					}
				}
			}
		}
		nComponents=idGenerator+1;
	
		// Computing component size
		connectedComponentSize = new int[nComponents];
		Math.set(connectedComponentSize,0,nComponents-1,0);
		for (i=0; i<nNodes; i++) connectedComponentSize[connectedComponent[i]]++;
	
		return nComponents;
	}
	
	
	/**
	 * To STDERR
	 */
	public static final void printSequenceLengths() {
		int i, j;
		
		for (j=sequenceLengths.length-1; j>=0; j--) {
			if (sequenceLengths[j]!=0) break;
		}
		for (i=0; i<=j; i++) System.err.println(i+","+RepeatAlphabet.sequenceLengths[i]);
	}
	
	
	/**
	 * Like $collectCharacterInstances()$, but looks up in the alphabet every character of
	 * every block of a recoded read.
	 *
	 * @param lastTranslatedRead the index in $Reads.readIDs$ of the last read that has 
	 * already been translated (-1 if no read has been translated yet);
	 * @param isLastChunk TRUE iff this is the last chunk of alignments to be processed;
	 * @param outputFile_characters translation of every read (one per line); a read has
	 * an empty line if it contains no repeat;
	 * @param outputFile_boundaries character boundaries (one read per line);
	 * @param outputFile_fullyUniqueReads list of IDs of reads that contain no repeat
	 * (zero-based);
	 * @param outputFile_oneRepeatReads list of IDs of reads that are fully contained in a
	 * single repeat (zero-based).
	 */
	public static final void translateReads(String alignmentsFile, int lastTranslatedRead, double maxError, int minAlignmentLength, int quantum, boolean isLastChunk, String outputFile_characters, String outputFile_boundaries, String outputFile_histogram, String outputFile_fullyUniqueReads, String outputFile_oneRepeatReads) throws IOException {
		final int ALIGNMENTS_CAPACITY = 100000;  // Arbitrary
		final int SEQUENCE_CAPACITY = 1000000;  // Arbitrary
		final int MAX_HISTOGRAM_LENGTH = 1000;  // Arbitrary
		int i, j;
		int row, readA, previousReadA;
		String str;
		BufferedReader br;
		BufferedWriter bw1, bw2, bw3, bw4;
		long[] histogram;
		
		// Allocating memory
		if (alignments==null || alignments.length<ALIGNMENTS_CAPACITY) alignments = new AlignmentRow[ALIGNMENTS_CAPACITY];
		for (i=0; i<alignments.length; i++) {	
			if (alignments[i]==null) alignments[i] = new AlignmentRow();
		}
		if (sequence==null || sequence.length<SEQUENCE_CAPACITY) sequence = new Block[SEQUENCE_CAPACITY];
		for (i=0; i<sequence.length; i++) {
			if (sequence[i]==null) sequence[i] = new Block();
		}
		if (points==null || points.length<(ALIGNMENTS_CAPACITY)<<1) points = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (periodicIntervals==null || periodicIntervals.length<(ALIGNMENTS_CAPACITY)<<1) periodicIntervals = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (newBlock==null) newBlock = new Block();
		if (tmpCharacter==null) tmpCharacter = new Character();
		histogram = new long[MAX_HISTOGRAM_LENGTH+1];
		Math.set(histogram,MAX_HISTOGRAM_LENGTH,0L);
		if (stack==null || stack.length<ALIGNMENTS_CAPACITY) stack = new int[ALIGNMENTS_CAPACITY];
		
		// Translating every read using the alphabet
		bw1 = new BufferedWriter(new FileWriter(outputFile_characters));
		bw2 = new BufferedWriter(new FileWriter(outputFile_boundaries));
		bw3 = new BufferedWriter(new FileWriter(outputFile_fullyUniqueReads));
		bw4 = new BufferedWriter(new FileWriter(outputFile_oneRepeatReads));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); 
		previousReadA=-1; lastAlignment=-1; row=0; j=lastTranslatedRead+1;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			if ( (2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2)>maxError ||
				 Math.max(Alignments.endA-Alignments.startA,Alignments.endB-Alignments.startB)+1<minAlignmentLength
			   ) {
				str=br.readLine(); row++;
				continue;
			}
			readA=Alignments.readA-1;
			if (previousReadA==-1 || readA!=previousReadA) {
				if (previousReadA!=-1) {
					while (Reads.readIDs[j]<previousReadA) {
						bw1.newLine(); bw2.newLine(); histogram[0]++;
						bw3.write(Reads.readIDs[j]+"\n");
						j++;
					}
					cleanAlignments(quantum);
					if (lastAlignment!=-1) {
						recodeRead(quantum);
						if (lastInSequence==-1) {
							bw1.newLine(); bw2.newLine(); histogram[0]++;
							bw3.write(previousReadA+"\n");
						}
						else if (lastInSequence==0) {
							translateRead(bw1,bw2,quantum);
							histogram[1]++;
							if (sequence[0].isUnique()) bw3.write(previousReadA+"\n");
							else bw4.write(previousReadA+"\n");
						}
						else {
							translateRead(bw1,bw2,quantum);
							histogram[lastInSequence+1]++;							
						}						
					}
					else { 
						bw1.newLine(); bw2.newLine(); histogram[0]++; 
						bw3.write(previousReadA+"\n");
					}
					j++;
				}
				previousReadA=readA; lastAlignment=0;
				alignments[0].set(readA,Math.max(Alignments.startA,0),Math.min(Alignments.endA,Reads.getReadLength(readA)-1),Alignments.readB-1,Math.max(Alignments.startB,0),Math.min(Alignments.endB,repeatLengths[Alignments.readB-1]-1),Alignments.orientation,Alignments.diffs);
			}
			else {
				lastAlignment++;
				if (lastAlignment==alignments.length) {
					AlignmentRow[] newAlignments = new AlignmentRow[alignments.length<<1];
					System.arraycopy(alignments,0,newAlignments,0,alignments.length);
					for (i=alignments.length; i<newAlignments.length; i++) newAlignments[i] = new AlignmentRow();
					alignments=newAlignments;
				}
				alignments[lastAlignment].set(readA,Math.max(Alignments.startA,0),Math.min(Alignments.endA,Reads.getReadLength(readA)-1),Alignments.readB-1,Math.max(Alignments.startB,0),Math.min(Alignments.endB,repeatLengths[Alignments.readB-1]-1),Alignments.orientation,Alignments.diffs);
			}
			str=br.readLine(); row++;
		}
		br.close();
		if (previousReadA!=-1) {
			while (Reads.readIDs[j]<previousReadA) {
				bw1.newLine(); bw2.newLine(); histogram[0]++;
				bw3.write(Reads.readIDs[j]+"\n");
				j++;
			}
			cleanAlignments(quantum);
			if (lastAlignment!=-1) {
				recodeRead(quantum);
				if (lastInSequence==-1) { 
					bw1.newLine(); bw2.newLine(); histogram[0]++; 
					bw3.write(previousReadA+"\n");
				}
				else if (lastInSequence==0) {
					translateRead(bw1,bw2,quantum);
					histogram[1]++; 
					if (sequence[0].isUnique()) bw3.write(previousReadA+"\n");
					else bw4.write(previousReadA+"\n");
				}
				else {
					translateRead(bw1,bw2,quantum);
					histogram[lastInSequence+1]++;
				}
			}
			else { 
				bw1.newLine(); bw2.newLine(); histogram[0]++;
				bw3.write(previousReadA+"\n");
			}
			j++;
		}
		if (isLastChunk) {
			previousReadA++;
			while (previousReadA<=Reads.lastRead) {
				bw1.newLine(); bw2.newLine(); histogram[0]++;
				bw3.write(previousReadA+"\n");
				previousReadA++;
			}
		}
		bw1.close(); bw2.close(); bw3.close(); bw4.close();
		bw1 = new BufferedWriter(new FileWriter(outputFile_histogram));
		for (i=0; i<=MAX_HISTOGRAM_LENGTH; i++) bw1.write(histogram[i]+"\n");
		bw1.close();
	}
	
	
	private static final void translateRead(BufferedWriter bw1, BufferedWriter bw2, int quantum) throws IOException {
		int i, j, k;
		
		for (i=0; i<=sequence[0].lastCharacter; i++) sequence[0].characters[i].quantize(quantum);
		if (sequence[0].isUnique()) translate_unique(sequence[0].characters[0],bw1);
		else if (sequence[0].characters[0].start==-1) translate_periodic(sequence[0],bw1);
		else translate(sequence[0],quantum,bw1);
		k=0;
		while (k<=lastPoint && points[k]<=quantum) k++;
		for (i=1; i<=lastInSequence; i++) {
			bw1.write(SEPARATOR_MAJOR+"");
			bw2.write((i>1?(SEPARATOR_MINOR+""):"")+points[k++]);
			for (j=0; j<=sequence[i].lastCharacter; j++) sequence[i].characters[j].quantize(quantum);
			if (sequence[i].isUnique()) translate_unique(sequence[i].characters[0],bw1);
			else if (sequence[i].characters[0].start==-1) translate_periodic(sequence[i],bw1);
			else translate(sequence[i],quantum,bw1);
		}
		bw1.newLine(); bw2.newLine();
	}

	
	/**
	 * Writes: $X$ if $character$ is closed and is identical to character $X$; $-1-X$ if
	 * $character$ is half-open and could be an instance of any closed block $>=X$;
	 * $lastAlphabet+1$ if $character$ is half-open but longer than any closed block.
	 *
	 * @param character assumed to be quantized.
	 */
	private static final void translate_unique(Character character, BufferedWriter bw) throws IOException {
		int i;
		final int length = character.length;
		
		i=Arrays.binarySearch(alphabet,0,lastUnique+1,character);
		if (character.isOpen()) {
			if (i<0) i=-1-i;
			if (i==lastUnique+1) bw.write((lastAlphabet+1)+"");
			else {
				while (i>=0 && alphabet[i].getLength()>=length) i--;
				i++;
				bw.write((-1-i)+"");
			}
		}
		else {
			if (i<0) {
				System.err.println("translate_unique> ERROR: closed unique character not found in the alphabet");
				System.err.println("query: "+character);
				System.err.println("candidate in alphabet: "+alphabet[-1-i]);
				System.exit(1);
			}
			bw.write(i+"");
		}
	}
	
	
	/**
	 * Same logic as $translateRead_unique()$, but taking orientation and open status of
	 * endpoints into account.
	 */
	private static final void translate_periodic(Block block, BufferedWriter bw) throws IOException {
		final int CAPACITY = 100;  // Arbitrary
		boolean orientation, openStart, openEnd;
		int i, j, k;
		int repeat, length, last;
		final int lastCharacter = block.lastCharacter;
		Character character;		
		
		// Collecting characters
		if (stack==null) stack = new int[CAPACITY];
		last=-1;
		for (j=0; j<=lastCharacter; j++) {
			character=block.characters[j];
			repeat=character.repeat; orientation=character.orientation; length=character.length;
			openStart=character.openStart; openEnd=character.openEnd;
			i=Arrays.binarySearch(alphabet,lastUnique+1,lastPeriodic+1,character);
			if (openStart || openEnd) {
				if (i<0) i=-1-i;
				if (alphabet[i].repeat!=repeat || alphabet[i].orientation!=orientation) i--;
				k=i;
				while (k>lastUnique && alphabet[k].repeat==repeat && alphabet[k].orientation==orientation && alphabet[k].length>=length) k--;
				k++;
				while (k<=lastPeriodic && alphabet[k].repeat==repeat && alphabet[k].orientation==orientation && !alphabet[k].implies(character)) k++;
				i=k;
				if (alphabet[i].repeat!=repeat || alphabet[i].orientation!=orientation || !alphabet[i].implies(character)) {
					System.err.println("translate_periodic> ERROR: open periodic repeat not found in the alphabet");
					System.err.println("query: "+character);
					System.err.println("first candidate in alphabet: "+alphabet[i]);
					System.exit(1);
				}
				last=appendToStack(-1-i,last);
			}
			else {
				if (i<0) {
					System.err.println("translate_periodic> ERROR: closed periodic character not found in the alphabet:");
					System.err.println("query: "+character);
					System.err.println("candidate in alphabet: "+alphabet[-1-i]);
					System.err.println("characters in the block:");
					for (int x=0; x<=lastCharacter; x++) System.err.println(block.characters[x]+"  isPeriodic="+isPeriodic[block.characters[x].repeat]);
					System.exit(1);
				}
				last=appendToStack(i,last);
			}
		}
		
		// Removing duplicates
		if (last>0) Arrays.sort(stack,0,last+1);
		j=stack[0]; bw.write(stack[0]+"");
		for (i=1; i<=last; i++) {
			if (stack[i]==j) continue;
			bw.write((SEPARATOR_MINOR+"")+stack[i]);
			j=stack[i];
		}
	}
	
	
	/**
	 * @param last last element currently in the stack;
	 * @return the new value of $last$.
	 */
	private static final int appendToStack(int value, int last) {
		last++;
		if (last==stack.length) {
			int[] newArray = new int[stack.length<<1];
			System.arraycopy(stack,0,newArray,0,stack.length);
			stack=newArray;
		}
		stack[last]=value;
		return last;
	}
	
	
	/**
	 * For every nonperiodic character $X$ in the block: if $X$ is half-open, prints the 
	 * sorted list of all characters in $alphabet$ that are similar to or imply $X$. If 
	 * $X$ is closed, prints the closed character in $alphabet$ that is identical to $X$.
	 *
	 * @param block its characters are assumed to be quantized.
	 */
	private static final void translate(Block block, int quantum, BufferedWriter bw) throws IOException {
		final int CAPACITY = 100;  // Arbitrary
		boolean orientation, found;
		int i, j, k;
		int last, repeat, length;
		final int lastCharacter = block.lastCharacter;
		Character character;
		
		// Collecting characters
		if (stack==null) stack = new int[CAPACITY];
		last=-1;
		for (k=0; k<=lastCharacter; k++) {
			character=block.characters[k];
			repeat=character.repeat; orientation=character.orientation; length=character.length;
			i=Arrays.binarySearch(alphabet,lastPeriodic+1,lastAlphabet+1,character);
			if (character.isOpen()) {
				if (i<0) i=-1-i;
				found=false;
				for (j=i; j<=lastAlphabet; j++) {
					if (alphabet[j].repeat!=repeat || alphabet[j].orientation!=orientation) break;
					if (alphabet[j].implies(character)) {
					   found=true;
					   last=appendToStack(j,last);
					}
				}
				for (j=i-1; j>lastPeriodic; j--) {
					if (alphabet[j].repeat!=repeat || alphabet[j].orientation!=orientation) break;
					if (alphabet[j].implies(character)) {
						found=true;
						last=appendToStack(j,last);
					}
				}
				if (!found) {
					System.err.println("translate> ERROR: open nonperiodic character not found in the alphabet:");
					System.err.println("query: "+character);
					System.err.println("candidate in alphabet: "+alphabet[Math.min(i,lastAlphabet)]);
					System.exit(1);
				}
			}
			else {
				if (i<0) {
					System.err.println("translate> ERROR: closed nonperiodic character not found in the alphabet:");
					System.err.println("query: "+character);
					System.err.println("candidate in alphabet: "+alphabet[-1-i]);
					System.exit(1);
				}
				last=appendToStack(i,last);
			}
		}
		
		// Removing duplicates
		if (last>0) Arrays.sort(stack,0,last+1);
		j=stack[0]; bw.write(stack[0]+"");
		for (i=1; i<=last; i++) {
			if (stack[i]==j) continue;
			bw.write((SEPARATOR_MINOR+"")+stack[i]);
			j=stack[i];
		}
	}
	
	
	/**
	 * Adds to $characterCount$ every character instance in a row $str$ of the translated 
	 * reads file.
	 * 
	 * @param characterCount one cell per character of the alphabet, plus one.
	 */
	public static final void incrementCharacterCounts(String str, long[] characterCount) {
		boolean orientation;
		int i, j, k;
		int c, to, nBlocks, last, repeat;
		
		if (str.length()==0) return;
		nBlocks=loadBlocks(str);
		for (i=0; i<nBlocks; i++) {
			last=lastInBlock[i];
			for (j=0; j<=last; j++) {
				c=Integer.parseInt(blocks[i][j]);
				if (c<0) {
					c=-1-c;
					if (c<=lastUnique) {
						for (k=c; k<=lastUnique; k++) characterCount[k]++;
					}
					else {
						repeat=alphabet[c].repeat; orientation=alphabet[c].orientation;
						for (k=c; k<=lastPeriodic; k++) {
							if (alphabet[k].repeat!=repeat || alphabet[k].orientation!=orientation) break;
							if (alphabet[k].implies(alphabet[c])) characterCount[k]++;
						}
					}
				}
				else characterCount[c]++;
			}
		}
	}
	
	
	/**
	 * Let $alphabet[x_0],...,alphabet[x_n]$ be all and only the characters that 
	 * correspond to the same substring of the same repeat, regardless of orientation and 
	 * open status. The procedure resets $characterCount[x_i]$ to $\sum_{j=0}^{n} 
	 * characterCount[x_j]$ for every $i$.
	 *
	 * Remark: discarding open status is necessary, because $Character.implies()$ works 
	 * only in one orientation, so the same substring of the same repeat might have 
	 * different open status in different orientations.
	 *
	 * @param marked temporary space, of size at least $lastAlphabet+1$.
	 */
	public static final void symmetrizeCharacterCounts(long[] characterCount, boolean[] marked) {
		int i, j;
		int repeat, start, end, length;
		long sum;
		
		Math.set(marked,lastAlphabet,false);
		
		// Periodic
		for (i=lastUnique+1; i<lastPeriodic; i++) {
			if (marked[i]) continue;
			repeat=alphabet[i].repeat; length=alphabet[i].length; sum=characterCount[i];
			for (j=i+1; j<=lastPeriodic; j++) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].length==length) sum+=characterCount[j];
			}
			for (j=i-1; j>lastUnique; j--) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].length==length) sum+=characterCount[j];
			}
			characterCount[i]=sum;
			for (j=i+1; j<=lastPeriodic; j++) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].length==length) {
					characterCount[j]=sum;
					marked[j]=true;
				}
			}
			for (j=i-1; j>lastUnique; j--) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].length==length) {
					characterCount[j]=sum;
					marked[j]=true;
				}
			}
		}
		
		// Nonperiodic
		for (i=lastPeriodic+1; i<lastAlphabet; i++) {
			if (marked[i]) continue;
			repeat=alphabet[i].repeat; start=alphabet[i].start; end=alphabet[i].end;
			sum=characterCount[i];
			for (j=i+1; j<=lastAlphabet; j++) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].start==start && alphabet[j].end==end) sum+=characterCount[j];
			}
			for (j=i-1; j>lastPeriodic; j--) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].start==start && alphabet[j].end==end) sum+=characterCount[j];
			}
			characterCount[i]=sum;
			for (j=i+1; j<=lastAlphabet; j++) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].start==start && alphabet[j].end==end) {
					characterCount[j]=sum;
					marked[j]=true;
				}
			}
			for (j=i-1; j>lastPeriodic; j--) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].start==start && alphabet[j].end==end) {
					characterCount[j]=sum;
					marked[j]=true;
				}
			}
		}
	}
	
	
	/**
	 * Builds a histogram of symmetrized $characterCount$ values. Rows: counts. Columns:
	 * 0: unique, open; 
	 * 1: unique, closed;
	 * 2: periodic;
	 * 3: nonperiodic.
	 *
	 * @param marked the same array used by $symmetrizeCharacterCounts()$, after that
	 * procedure completes.
	 */
	public static final long[][] getCharacterHistogram(long[] characterCount, boolean[] marked, int maxFrequency) throws IOException {
		int i;
		long count;
		long[][] characterHistogram;

		// Unique
		characterHistogram = new long[maxFrequency+1][4];
		Math.set(characterHistogram,0);
		for (i=0; i<=lastUnique; i++) {
			count=characterCount[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][1/*All closed*/]++;
		}
		count=characterCount[lastAlphabet+1];
		if (count>maxFrequency) count=maxFrequency;
		characterHistogram[(int)count][0/*Open*/]++;
		// Periodic
		for (i=lastUnique+1; i<=lastPeriodic; i++) {
			if (marked[i]) continue;
			count=characterCount[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][2]++;
		}		
		// Nonperiodic
		for (i=lastPeriodic+1; i<=lastAlphabet; i++) {
			if (marked[i]) continue;
			count=characterCount[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][3]++;
		}
		return characterHistogram;
	}

	
	/**
	 * Loads the recoded read $str$ in global variables $blocks,lastInBlock$.
	 *
	 * @return the number of blocks in $blocks$.
	 */
	private static final int loadBlocks(String str) {
		int i;
		int nBlocks;
		String[] tokens;
		
		tokens=str.split(SEPARATOR_MAJOR+"");
		nBlocks=tokens.length;
		if (blocks==null) blocks = new String[nBlocks][0];
		else if (blocks.length<nBlocks) {
			String[][] newArray = new String[nBlocks][0];
			System.arraycopy(blocks,0,newArray,0,blocks.length);
			blocks=newArray;
		}
		if (lastInBlock==null || lastInBlock.length<nBlocks) lastInBlock = new int[nBlocks];
		Math.set(lastInBlock,nBlocks-1,-1);
		for (i=0; i<nBlocks; i++) {
			if (tokens[i].indexOf(SEPARATOR_MINOR+"")>=0) {
				blocks[i]=tokens[i].split(SEPARATOR_MINOR+"");
				lastInBlock[i]=blocks[i].length-1;
			}
			else {
				if (blocks[i].length==0) blocks[i] = new String[1];
				blocks[i][0]=tokens[i];
				lastInBlock[i]=0;
			}
		}
		return nBlocks;
	}
	
	
	/**
	 * Transforms string matrix $blocks$ into integer matrix $intBlocks$, expanding every 
	 * negative ID into a list of positive IDs, and builds array $isBlockUnique$.
	 *
	 * Remark: the procedure follows the logic of $incrementCharacterCounts()$.
	 *
	 * @param boundaries of each block;
	 * @param tmpChar temporary space.
	 */
	private static final void loadIntBlocks(int nBlocks, int[] boundaries, int readLength, Character tmpChar) {
		final int QUANTUM = IO.quantum;
		boolean orientation, unique, found;
		int i, j, k, c;
		int to, last, value, nElements, repeat;
		
		// Allocating memory
		if (intBlocks==null) intBlocks = new int[nBlocks][0];
		else if (intBlocks.length<nBlocks) {
			int[][] newArray = new int[nBlocks][0];
			System.arraycopy(intBlocks,0,newArray,0,intBlocks.length);
			intBlocks=newArray;
		}
		if (lastInBlock_int==null || lastInBlock_int.length<nBlocks) lastInBlock_int = new int[nBlocks];
		if (isBlockUnique==null || isBlockUnique.length<nBlocks) isBlockUnique = new boolean[nBlocks];
		Math.set(isBlockUnique,nBlocks-1,false);
		
		// Building arrays
		for (i=0; i<nBlocks; i++) {
			last=lastInBlock[i];
			nElements=0;
			for (j=0; j<=last; j++) {
				value=Integer.parseInt(blocks[i][j]);
				found=false;
				if (value>=0) { found=true; nElements++; }
				else {
					value=-1-value;
					if (value<=lastUnique) { found=true; nElements+=lastUnique+1-value; }
					else {
						tmpChar.copyFrom(alphabet[value]);
						tmpChar.length=(i==nBlocks-1?readLength:boundaries[i])-(i==0?0:boundaries[i-1]+1);
						if (i==0) tmpChar.openStart=true;
						if (i==nBlocks-1) tmpChar.openEnd=true;
						repeat=alphabet[value].repeat; orientation=alphabet[value].orientation;
						tmpChar.quantize(QUANTUM);
						for (k=value; k<=lastPeriodic; k++) {
							if (alphabet[k].repeat!=repeat || alphabet[k].orientation!=orientation) break;
							if (alphabet[k].implies(tmpChar)) { found=true; nElements++; }
						}
						if (!found) {
							// Looking for a shorter version, since the length of 
							// character $value$ was not computed in the same way as we do
							// here.
							tmpChar.length-=QUANTUM;
							for (k=value; k<=lastPeriodic; k++) {
								if (alphabet[k].repeat!=repeat || alphabet[k].orientation!=orientation) break;
								if (alphabet[k].implies(tmpChar)) { found=true; nElements++; }
							}
						}
					}
				}
				if (!found) {
					System.err.println("loadIntBlocks> ERROR: character ID "+blocks[i][j]+" in a translated read is not implied by any character in the alphabet.");
					System.err.println("tmpChar: "+tmpChar);
					System.err.println("alphabet["+value+"]: "+alphabet[value]);
					System.err.println("readLength="+readLength);
					System.err.println("blocks: ");
					for (int x=0; x<nBlocks; x++) {
						System.err.print(x+": ");
						for (int y=0; y<=lastInBlock[x]; y++) System.err.print(blocks[x][y]+",");
						System.err.println();
					}
					System.err.print("boundaries: ");
					for (int x=0; x<nBlocks-1; x++) System.err.print(boundaries[x]+",");
					System.err.println();
					System.exit(1);
				}
			}
			if (intBlocks[i].length<nElements) intBlocks[i] = new int[nElements];
			k=-1; unique=false;
			for (j=0; j<=last; j++) {
				value=Integer.parseInt(blocks[i][j]);
				if (value>=0) {
					intBlocks[i][++k]=value;
					if (value==lastAlphabet+1 || value<=lastUnique) unique=true;
				}
				else {
					value=-1-value;
					if (value<=lastUnique) {
						for (c=value; c<=lastUnique; c++) intBlocks[i][++k]=c;
						unique=true;
					}
					else {
						tmpChar.copyFrom(alphabet[value]);
						tmpChar.length=(i==nBlocks-1?readLength:boundaries[i])-(i==0?0:boundaries[i-1]+1);
						if (i==0) tmpChar.openStart=true;
						if (i==nBlocks-1) tmpChar.openEnd=true;
						repeat=alphabet[value].repeat; orientation=alphabet[value].orientation;
						tmpChar.quantize(QUANTUM);
						found=false;
						for (c=value; c<=lastPeriodic; c++) {
							if (alphabet[c].repeat!=repeat || alphabet[c].orientation!=orientation) break;
							if (alphabet[c].implies(tmpChar)) { found=true; intBlocks[i][++k]=c; }
						}
						if (!found) {
							// Looking for a shorter version, since the length of 
							// character $value$ was not computed in the same way as we do
							// here.
							tmpChar.length-=QUANTUM;
							for (c=value; c<=lastPeriodic; c++) {
								if (alphabet[c].repeat!=repeat || alphabet[c].orientation!=orientation) break;
								if (alphabet[c].implies(tmpChar)) { found=true; intBlocks[i][++k]=c; }
							}
						}
					}
				}
			}
			lastInBlock_int[i]=k; isBlockUnique[i]=unique;
		}
	}
	
	
	/**
	 * Stores in global variable $boundaries[0..X]$ the contents of $str$, where $X$ is
	 * returned in output.
	 */
	public static final int loadBoundaries(String str) {
		int i;
		int nBoundaries;
		String[] tokens;
		
		if (str.length()==0) nBoundaries=0;
		else if (str.indexOf(SEPARATOR_MINOR+"")>=0) {
			tokens=str.split(SEPARATOR_MINOR+"");
			nBoundaries=tokens.length;
			if (boundaries==null || boundaries.length<nBoundaries) boundaries = new int[nBoundaries];
			for (i=0; i<nBoundaries; i++) boundaries[i]=Integer.parseInt(tokens[i]);
		}
		else {
			nBoundaries=1;
			if (boundaries==null || boundaries.length<nBoundaries) boundaries = new int[nBoundaries];
			boundaries[0]=Integer.parseInt(str);
		}
		return nBoundaries-1;
	}



	
	// ----------------------- ALPHABET CLEANING PROCEDURES ------------------------------

	public static final void loadAlphabetCount(String file, int alphabetSize) throws IOException {
		int i;
		BufferedReader br;

		alphabetCount = new long[alphabetSize+1];
		br = new BufferedReader(new FileReader(file));
		for (i=0; i<alphabetSize+1; i++) alphabetCount[i]=Long.parseLong(br.readLine());
		br.close();
	}
	
	
	/**
	 * Appends to $bw$ the new unique closed characters that result from removing every 
	 * character with count smaller than $minCount$ from $read2characters$, and updates 
	 * $maxOpenLength_unique$ as well.
	 *
	 * Remark: the procedure assumes that global variables $alphabet,alphabetCount$ have
	 * already been initialized, and it uses global temporary variable $boundaries$.
	 *
	 * @param keepPeriodic TRUE=does not filter out rare periodic characters; periodic
	 * characters are often supported by a strong signal (of several overlapping read-
	 * repeat alignments), and removing them and replacing them with non-repetitive
	 * characters does not allow to remove several spurious alignments downstream;
	 * @param read2boundaries boundaries of the characters in $read2characters$;
	 * @param tmpChar temporary space.
	 */
	public static final void cleanTranslatedRead_collectCharacterInstances(String read2characters, String read2boundaries, int readLength, int minCount, boolean keepPeriodic, int quantum, BufferedWriter bw, Character tmpChar) throws IOException {
		boolean isUnique;
		int i, j, k;
		int c, length, first, nBlocks, nBoundaries;
		String[] tokens;
		
		if (read2characters.length()==0) return;
		nBoundaries=loadBoundaries(read2boundaries)+1;
		nBlocks=loadBlocks(read2characters);
		removeRareCharacters(nBlocks,minCount,lastUnique,lastPeriodic,lastAlphabet,keepPeriodic);
		first=-1;
		for (i=0; i<nBlocks; i++) {
			if (lastInBlock[i]==-1) isUnique=true;
			else {
				j=Integer.parseInt(blocks[i][0]);
				if (j<0) j=-1-j;
				isUnique=j<=lastUnique||j==lastAlphabet+1;
			}
			if (!isUnique) {
				if (first!=-1) {
					tmpChar.repeat=UNIQUE;
					tmpChar.orientation=true;
					tmpChar.start=-1; tmpChar.end=-1;
					length=boundaries[i-1]-(first==0?0:boundaries[first-1]);
					tmpChar.length=length;
					tmpChar.openStart=first==0; tmpChar.openEnd=false;
					tmpChar.quantize(quantum);
					j=Arrays.binarySearch(alphabet,0,lastUnique+1,tmpChar);
					if (tmpChar.isOpen()) {
						if (j<0) j=-1-j;
						if (j==lastUnique+1) maxOpenLength_unique=Math.max(maxOpenLength_unique,length);
					}
					else if (j<0) bw.write(tmpChar.toString()+"\n");
				}
				first=-1;
			}
			else if (first==-1) first=i;
		}
		if (first!=-1) {
			tmpChar.repeat=UNIQUE;
			tmpChar.orientation=true;
			tmpChar.start=-1; tmpChar.end=-1;
			length=readLength-(first==0?0:boundaries[first-1]);
			tmpChar.length=length;
			tmpChar.openStart=first==0; tmpChar.openEnd=true;
			tmpChar.quantize(quantum);
			j=Arrays.binarySearch(alphabet,0,lastUnique+1,tmpChar);
			if (tmpChar.isOpen()) {
				if (j<0) j=-1-j;
				if (j==lastUnique+1) maxOpenLength_unique=Math.max(maxOpenLength_unique,length);
			}
			else if (j<0) bw.write(tmpChar.toString()+"\n");
		}
	}

	
	/**
	 * Removes from $blocks$ all non-unique characters with count smaller than $minCount$.
	 *
	 * Remark: the procedure assumes that global variable $alphabetCount$ has already been
	 * initialized.
	 *
	 * @param keepPeriodic TRUE=does not remove rare periodic characters.
	 */
	private static final void removeRareCharacters(int nBlocks, int minCount, int lastUnique, int lastPeriodic, int lastAlphabet, boolean keepPeriodic) {
		boolean found, deleted, orientation;
		int i, j, k, h;
		int c, last, repeat;
		
		for (i=0; i<nBlocks; i++) {
			last=lastInBlock[i];
			k=-1;
			for (j=0; j<=last; j++) {
				c=Integer.parseInt(blocks[i][j]);
				deleted=false;
				if (c<0) {
					c=-1-c;
					if (c>lastUnique) {
						if ((keepPeriodic?c>lastPeriodic:true) && alphabetCount[c]<minCount) {
							repeat=alphabet[c].repeat; orientation=alphabet[c].orientation;
							h=c+1; found=false;
							while (h<=lastPeriodic) {
								if (alphabet[h].repeat!=repeat || alphabet[h].orientation!=orientation) break;
								if (alphabetCount[h]>=minCount && alphabet[h].implies(alphabet[c])) {
									found=true;
									break;
								}
								h++;
							}
							if (!found) deleted=true;
							else blocks[i][j]=h+"";
						}
					}
					else {
						// NOP: we don't filter out unique characters, and non-periodic
						// repeats cannot be negative.
					}
				}
				else if (c>lastUnique && (keepPeriodic?c>lastPeriodic:true) && alphabetCount[c]<minCount) deleted=true;
				if (!deleted) blocks[i][++k]=blocks[i][j];
			}
			lastInBlock[i]=k;
		}
	}
	
	
	/**
	 * Removes from $alphabet$ all characters with $alphabetCount < minCount$, and adds to
	 * $alphabet$ all the new unique characters in $newCharactersFile$.
	 *
	 * Remark: $alphabetCount$ is not valid after the procedure completes.
	 *
	 * @param keepPeriodic TRUE=does not remove rare periodic characters;
	 * @return an array that maps every non-unique character in the original $alphabet$ to
	 * its position in the new $alphabet$ (or to -1 if the character does not appear in 
	 * the new alphabet). Positions are relative to the corresponding values of 
	 * $lastUnique+1$.
	 */
	public static final int[] cleanTranslatedRead_updateAlphabet(int nNewCharacters, String newCharactersFile, int minCount, boolean keepPeriodic) throws IOException {
		int i, j;
		int lastPeriodicPrime;
		String str;
		Character tmpChar;
		BufferedReader br;
		int[] old2new;
		Character[] newAlphabet;
		
		// Removing rare characters
		old2new = new int[lastAlphabet-lastUnique];
		Math.set(old2new,old2new.length-1,-1);
		if (keepPeriodic) {
			for (i=lastUnique+1; i<=lastPeriodic; i++) old2new[i-lastUnique-1]=i-lastUnique-1;
			j=lastPeriodic; lastPeriodicPrime=lastPeriodic; 
		}
		else { j=lastUnique; lastPeriodicPrime=-1; }
		for (i=keepPeriodic?lastPeriodic+1:lastUnique+1; i<=lastAlphabet; i++) {
			if (alphabetCount[i]<minCount) continue;
			j++;
			tmpChar=alphabet[j];
			alphabet[j]=alphabet[i];
			alphabet[i]=tmpChar;
			if (alphabet[j].start==-1) lastPeriodicPrime=j;
			old2new[i-lastUnique-1]=j-lastUnique-1;
		}
		lastAlphabet=j;
		lastPeriodic=lastPeriodicPrime!=-1?lastPeriodicPrime:lastUnique;
		
		// Adding new unique characters
		if (nNewCharacters>0) {
			if (alphabet.length>=lastAlphabet+1+nNewCharacters) {
				for (i=lastAlphabet; i>lastUnique; i--) {
					alphabet[i+nNewCharacters]=alphabet[i];
					alphabet[i]=null;
				}
				br = new BufferedReader(new FileReader(newCharactersFile));
				for (i=0; i<nNewCharacters; i++) {
					alphabet[lastUnique+1+i] = new Character();
					alphabet[lastUnique+1+i].deserialize(br.readLine());
				}
				br.close();
				lastUnique+=nNewCharacters; lastPeriodic+=nNewCharacters; lastAlphabet+=nNewCharacters;
				Arrays.sort(alphabet,0,lastUnique+1);
			}
			else {
				newAlphabet = new Character[lastAlphabet+1+nNewCharacters];
				System.arraycopy(alphabet,0,newAlphabet,0,lastUnique+1);
				br = new BufferedReader(new FileReader(newCharactersFile));
				for (i=0; i<nNewCharacters; i++) {
					newAlphabet[lastUnique+1+i] = new Character();
					newAlphabet[lastUnique+1+i].deserialize(br.readLine());
				}
				br.close();
				i=lastUnique+1+nNewCharacters;
				Arrays.sort(newAlphabet,0,i);
				System.arraycopy(alphabet,lastUnique+1,newAlphabet,i,lastPeriodic-lastUnique);
				i=lastPeriodic+nNewCharacters;
				System.arraycopy(alphabet,lastPeriodic+1,newAlphabet,i,lastAlphabet-lastPeriodic);
				alphabet=newAlphabet;
				lastUnique+=nNewCharacters; lastPeriodic+=nNewCharacters; lastAlphabet+=nNewCharacters;
			}
		}
		
		return old2new;
	}
	
	
	/**
	 * Updates the translation of a read into characters, assuming that $alphabetCount$ 
	 * refers to the OLD alphabet before $cleanTranslatedRead_updateAlphabet()$.
	 * 
	 * Remark: real repeats might be inadvertedly discarded by this process, since the
	 * frequency threshold is not perfect. Discarded repeat characters are transformed 
	 * into non-repetitive: this completely discards the information that the substring 
	 * of the read is similar to a repeat (this information might have been particularly 
	 * useful if the discarded repeat was periodic), and it might affect alignment 
	 * filtering downstream (any alignment that contains this new nonrepetitive character
	 * might be kept, precisely because of this character). We keep things as they are for
	 * simplicity.
	 *	
	 * @param read2characters_old old translation;
	 * @param read2boundaries_old old block boundaries of the translation;
	 * @param newAlphabet obtained from the old alphabet by running $cleanTranslatedRead_
	 * updateAlphabet()$;
	 * @param old2new the output of $cleanTranslatedRead_updateAlphabet()$;
	 * @param read2characters_new,read2boundaries_new output files;
	 * @return the number of blocks in the new translation.
	 */
	public static final int cleanTranslatedRead_updateTranslation(String read2characters_old, String read2boundaries_old, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, int[] old2new, int readLength, int minCount, boolean keepPeriodic, int quantum, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, Character tmpChar) throws IOException {
		boolean isUnique;
		int i, j;
		int c, length, first, last, nBlocks, nBoundaries, nAppendedBlocks;
		String[] tokens;
		
		if (read2characters_old.length()==0) return 0;
		nBoundaries=loadBoundaries(read2boundaries_old)+1;
		nBlocks=loadBlocks(read2characters_old);
		alphabet=oldAlphabet; lastUnique=lastUnique_old; lastPeriodic=lastPeriodic_old; lastAlphabet=lastAlphabet_old;
		removeRareCharacters(nBlocks,minCount,lastUnique_old,lastPeriodic_old,lastAlphabet_old,keepPeriodic);
		alphabet=newAlphabet; lastUnique=lastUnique_new; lastPeriodic=lastPeriodic_new; lastAlphabet=lastAlphabet_new;
		first=-1; nAppendedBlocks=0;
		for (i=0; i<nBlocks; i++) {
			if (lastInBlock[i]==-1) isUnique=true;
			else {
				j=Integer.parseInt(blocks[i][0]);
				if (j<0) j=-1-j;
				isUnique=j<=lastUnique_old||j==lastAlphabet_old+1;
			}
			if (!isUnique) {
				if (first!=-1) {
					tmpChar.repeat=UNIQUE;
					tmpChar.orientation=true;
					tmpChar.start=-1; tmpChar.end=-1;
					length=boundaries[i-1]-(first==0?0:boundaries[first-1]);
					tmpChar.length=length;
					tmpChar.openStart=first==0; tmpChar.openEnd=false;
					tmpChar.quantize(quantum);
					if (nAppendedBlocks>0) read2characters_new.write(SEPARATOR_MAJOR+"");				
					translate_unique(tmpChar,read2characters_new);
					if (first>0) {
						if (nAppendedBlocks>1) read2boundaries_new.write(SEPARATOR_MINOR+"");
						read2boundaries_new.write(boundaries[first-1]+"");
					}
					nAppendedBlocks++;
				}
				first=-1;
				if (nAppendedBlocks>0) read2characters_new.write(SEPARATOR_MAJOR+"");
				last=lastInBlock[i];
				for (j=0; j<=last; j++) {
					c=Integer.parseInt(blocks[i][j]);
					if (c<0) {
						c=-1-c;
						if (c<=lastUnique_old) {
							if (j>0) read2characters_new.write(SEPARATOR_MINOR+"");
							translate_unique(oldAlphabet[c],read2characters_new);
						}
						else {
							c=-1-(lastUnique+1+old2new[c-lastUnique_old-1]);
							read2characters_new.write((j>0?SEPARATOR_MINOR+"":"")+c);
						}
					}
					else if (c==lastAlphabet_old+1) {
						tmpChar.repeat=UNIQUE;
						tmpChar.orientation=true;
						tmpChar.start=-1; tmpChar.end=-1;
						length=(i==nBlocks-1?readLength:boundaries[i])-(i==0?0:boundaries[i-1]);
						tmpChar.length=length;
						tmpChar.openStart=i==0;
						tmpChar.openEnd=i==nBlocks-1;
						tmpChar.quantize(quantum);
						if (j>0) read2characters_new.write(SEPARATOR_MINOR+"");
						translate_unique(tmpChar,read2characters_new);
					}
					else if (c<=lastUnique_old) {
						if (j>0) read2characters_new.write(SEPARATOR_MINOR+"");
						translate_unique(oldAlphabet[c],read2characters_new);
					}
					else {
						c=lastUnique+1+old2new[c-lastUnique_old-1];
						read2characters_new.write((j>0?SEPARATOR_MINOR+"":"")+c);
					}
				}
				if (i>0) {
					if (nAppendedBlocks>1) read2boundaries_new.write(SEPARATOR_MINOR+"");
					read2boundaries_new.write(boundaries[i-1]+"");
				}
				nAppendedBlocks++;
			}
			else if (first==-1) first=i;
		}
		if (first!=-1) {
			if (nAppendedBlocks>0) {
				tmpChar.repeat=UNIQUE;
				tmpChar.orientation=true;
				tmpChar.start=-1; tmpChar.end=-1;
				length=readLength-(first==0?0:boundaries[first-1]);
				tmpChar.length=length;
				tmpChar.openStart=first==0; tmpChar.openEnd=true;
				tmpChar.quantize(quantum);
				if (nAppendedBlocks>0) read2characters_new.write(SEPARATOR_MAJOR+"");
				translate_unique(tmpChar,read2characters_new);
				if (first>0) {
					if (nAppendedBlocks>1) read2boundaries_new.write(SEPARATOR_MINOR+"");
					read2boundaries_new.write(boundaries[first-1]+"");
				}
				nAppendedBlocks++;
			}
			else {
				// Reads with a single unique block are recoded as having no block
				return 0;
			}
		}
		return nAppendedBlocks;
	}

	

	
	// ------------------------------ K-MER PROCEDURES -----------------------------------
	
	/**
	 * Tells whether a one-mer that is marked as unique based just on frequency, should be 
	 * used to mark a read block as a unique address on the genome. M=multiple characters 
	 * occur in the read block; E=the block is the first/last one of the read; 
	 * oneMerFilter=selection criterion.
     *
	 * M  E   oneMerFilter
	 *        0  1  2  3
	 * 0  0	  1  1  1  1
     * 0  1   1  1  1  0
	 * 1  0   1  0  1  1
	 * 1  1   1  0  0  0
	 */
	private static final boolean[][][] ONEMER_TABLE = new boolean[][][] {
		{ {true,true,true,true}, {true,true,true,false} },
		{ {true,false,true,true}, {true,false,false,false} }
	};
	
	
	/**
	 * If $newKmers$ is not null, the procedure uses every length-k window that satisfies 
	 * the conditions in $uniqueMode,multiMode$ (see procedure $isValidWindow()$ for 
	 * details) to add a k-mer to $newKmers$. If a block contains multiple characters, 
	 * every character can be used to build a k-mer. K-mers are canonized before being 
	 * added to $newKmers$ (see $Kmer.canonize()$). Array $avoidedIntervals$ contains 
	 * tuples (position,length,nHaplotypes) sorted by position: if a window contains one 
	 * such interval, it is not used to build k-mers.
	 *
	 * If $newKmers$ is null, the procedure checks instead if a window contains a k-mer in 
	 * $oldKmers$, and if so it appends a (position,length,nHaplotypes) tuple to 
	 * $avoidedIntervals$ ($nHaplotypes$ is decided according to $haplotypeCoverage$).
	 * The new value of $lastAvoidedInterval$ is returned in output.
     *
	 * Remark: one-mers collected by this procedure might have a different (and even 
	 * smaller) count than the one produced by the $getCharacterHistogram()$ pipeline, 
	 * because of the constraints and of canonization. E.g. a one-mer might be rare in 
	 * this procedure but frequent in $getCharacterHistogram()$, if the latter counted all
	 * the occurrences of the character but the former counts just closed occurrences.
	 *
	 * Remark: often several characters map to the same block in practice, and most of 
	 * the k-mers (k>=2) that result from this are rare.
	 *
	 * Remark: the procedure enumerates k-mers everywhere, including strictly inside 
	 * tandems. Substrings that would generate too many k-mers are skipped to avoid large
	 * slowdowns.
	 *
	 * Remark: the procedure uses global variables $blocks,intBlocks,stack$.
	 * 
	 * @param oneMerFilter if $newKmers$ is null, use the criterion specified in 
	 * $ONEMER_TABLE[][][oneMerFilter]$ to decide whether a one-mer should be considered a 
	 * unique address on the genome (for $k>1$ we run no filter, since we assume that 
	 * contexts longer than one are enough to disambiguate; alternatively, one could try 
	 * to avoid unique k-mer occurrences that cover the first/last block of a read, as 
	 * well, but this is likely overengineering);
	 * @param haplotypeCoverage coverage of one haplotype;
	 * @param tmpKmer temporary space;
	 * @param tmpArray2 temporary space, of size at least k;
	 * @param tmpArray3 temporary space, of size at least 2k;
	 * @param tmpMap temporary hashmap, used only if $newKmers$ is not null.
	 */
	public static final int getKmers(String str, int k, int uniqueMode, int multiMode, int oneMerFilter, HashMap<Kmer,Kmer> newKmers, HashMap<Kmer,Kmer> oldKmers, int[] avoidedIntervals, int lastAvoidedInterval, int haplotypeCoverage, int readLength, int[] boundaries, Kmer tmpKmer, int[] tmpArray2, int[] tmpArray3, HashMap<Kmer,Kmer> tmpMap, Character tmpChar) {
		final int MAX_KMERS_TO_ENUMERATE = 500000;  // Arbitrary, just for speedup.
		int i, j, p;
		int nBlocks, sum, start, end, nHaplotypes, nKmers, out;
		Kmer key, value;
		Iterator<Kmer> iterator;
		
		// Loading blocks
		out=lastAvoidedInterval;
		if (str.length()<(k<<1)-1) return out;
		i=-1; nBlocks=0;
		while (true) {
			i=str.indexOf(SEPARATOR_MAJOR+"",i+1);
			if (i<0) break;
			else nBlocks++;
		}
		nBlocks++;
		if (nBlocks<k) return out;
		loadBlocks(str); loadIntBlocks(nBlocks,boundaries,readLength,tmpChar);
		sum=0;
		for (i=0; i<nBlocks; i++) sum+=lastInBlock_int[i]+1;
		if (stack==null || stack.length<sum*3) stack = new int[sum*3];
		
		// Processing every k-mer in the read
		if (newKmers==null) {
			j=0;
			for (i=0; i<=nBlocks-k; i++) {
				while (j<lastAvoidedInterval && avoidedIntervals[j]<i) j+=3;
				if (j<lastAvoidedInterval && avoidedIntervals[j]+avoidedIntervals[j+1]-1<=i+k-1) continue;
				if (!isValidWindow(i,k,nBlocks,uniqueMode,multiMode,readLength)) continue;
				nKmers=lastInBlock_int[i]+1;
				for (p=i+1; p<=i+k-1; p++) nKmers*=lastInBlock_int[p]+1;
				if (nKmers<0 || nKmers>MAX_KMERS_TO_ENUMERATE) continue;
				nHaplotypes=getKmers_impl(i,k,null,oldKmers,haplotypeCoverage,tmpKmer,stack,tmpArray2,tmpArray3);
				if (nHaplotypes!=-1 && (k>1 || ONEMER_TABLE[lastInBlock_int[i]>0?1:0][i==0||i==nBlocks-1?1:0][oneMerFilter])) { 
					avoidedIntervals[++out]=i; avoidedIntervals[++out]=k; avoidedIntervals[++out]=nHaplotypes; 
				}
			}
		}
		else {
			tmpMap.clear(); lastKmerPool=-1;
			j=0;
			for (i=0; i<=nBlocks-k; i++) {
				while (j<lastAvoidedInterval && avoidedIntervals[j]<i) j+=3;
				if (j<lastAvoidedInterval && avoidedIntervals[j]+avoidedIntervals[j+1]-1<=i+k-1) continue;			
				if (!isValidWindow(i,k,nBlocks,uniqueMode,multiMode,readLength)) continue;			
				nKmers=lastInBlock_int[i]+1;
				for (p=i+1; p<=i+k-1; p++) nKmers*=lastInBlock_int[p]+1;
				if (nKmers<0 || nKmers>MAX_KMERS_TO_ENUMERATE) continue;				
				getKmers_impl(i,k,tmpMap,oldKmers,haplotypeCoverage,tmpKmer,stack,tmpArray2,tmpArray3);
			}
			iterator=tmpMap.keySet().iterator();
			while (iterator.hasNext()) {
				key=iterator.next();
				value=newKmers.get(key);
				if (value==null) {
					value = new Kmer(key,k);
					value.count=key.count; value.sameReadCount=(int)key.count;
					newKmers.put(value,value);
				}
				else {
					value.sameReadCount=Math.max(value.sameReadCount,(int)key.count);
					value.count+=key.count;
				}
			}
		}
		
		return out;
	}

	
	/**
	 * Tells whether window $blocks[first..first+k-1]$ satisfies the following conditions.
	 *
	 * Remark: the procedure assumes that global variable $boundaries$ has been loaded 
	 * with the boundaries of the current blocks.
	 *
	 * @param uniqueMode blocks with unique characters:
	 * 0: are allowed; 
	 * 1: are allowed everywhere, except in the first and last block of the window;
	 * 2: are not allowed;
	 * @param multiMode blocks with multiple characters:
	 * 0: are allowed;
	 * 1: are not allowed iff they are the first/last one of the read; endblocks are more
	 *    likely to match several characters because they contain just a fraction of a 
	 *    character; allowing an endblock if it matches just one character is useful, 
	 *    since it might be the only way to detect e.g. a transposon that is longer than 
	 *    every read and that occurs just once in the genome, or an extremely long 
	 *    satellite that occurs just once in the genome.
	 * 2: are not allowed.
	 */
	private static final boolean isValidWindow(int first, int k, int nBlocks, int uniqueMode, int multiMode, int readLength) {
		final int MIN_BLOCK_LENGTH = IO.quantum<<1;  // Arbitrary
		int i;
		int start, end;
		
		if (uniqueMode==1) {
			if (isBlockUnique[first] || isBlockUnique[first+k-1]) return false;
		}
		else if (uniqueMode==2) {
			for (i=0; i<=k-1; i++) {
				if (isBlockUnique[first+i]) return false;
			}
		}
		if (multiMode==1) {
			if ((first==0 && lastInBlock_int[first]>0) || (first+k-1==nBlocks-1 && lastInBlock_int[first+k-1]>0)) return false;
		}
		else if (multiMode==2) {
			for (i=0; i<=k-1; i++) {
				if (lastInBlock_int[first+i]>0) return false;
			}
		}
		
		// Avoiding the whole window if it contains a short block
		if (nBlocks>1) {
			if (first==0) { start=0; end=boundaries[0]; }
			else if (first==nBlocks-1) { start=boundaries[nBlocks-2]; end=readLength-1; }
			else { start=boundaries[first-1]; end=boundaries[first]; }
			if (end-start+1<MIN_BLOCK_LENGTH) return false;
			if (first+k-1==0) { start=0; end=boundaries[0]; }
			else if (first+k-1==nBlocks-1) { start=boundaries[nBlocks-2]; end=readLength-1; }
			else { start=boundaries[first+k-2]; end=boundaries[first+k-1]; }
			if (end-start+1<MIN_BLOCK_LENGTH) return false;
			for (i=1; i<=k-2; i++) {
				if (boundaries[i]-boundaries[i-1]<MIN_BLOCK_LENGTH) return false;
			}
		}
		else { /* NOP: the only block of the read is assumed to be long. */ }
		
		return true;
	}
	
	
	public static final void kmerPool_init(int k) {
		kmerPool = new Kmer[100];  // Arbitrary
		for (int i=0; i<kmerPool.length; i++) {
			kmerPool[i] = new Kmer();
			kmerPool[i].sequence = new int[k];
		}
	}
	
	
	private static final Kmer kmerPool_allocate(int k) {
		final int GROWTH_RATE = 100;  // Arbitrary
		final int length = kmerPool.length;
		lastKmerPool++;
		if (lastKmerPool==length) {
			Kmer[] newArray = new Kmer[length+GROWTH_RATE];
			System.arraycopy(kmerPool,0,newArray,0,length);
			for (int i=length; i<newArray.length; i++) {
				newArray[i] = new Kmer();
				newArray[i].sequence = new int[k];
			}
			kmerPool=newArray;
		}
		return kmerPool[lastKmerPool];
	}
	
	
	/**
	 * If $newKmers$ is not null, adds every possible canonized instance of 
	 * $intBlocks[first..first+k-1]$ to $newKmers$. Otherwise, checks whether one of the 
	 * possible canonized instances of $intBlocks[first..first+k-1]$ belongs to
	 * $oldKmers$, and if so returns an estimate of the number of haplotypes in which it
	 * occurs (returns -1 otherwise).
	 *
	 * Remark: the objects added to $newKmers$ come from $kmerPool$, which is assumed to
	 * be already initialized.
	 *
	 * @param key temporary space;
	 * @param tmpArray1 temporary space, of size at least equal to 3 times the number of 
	 * elements in $intBlocks$;
	 * @param tmpArray2 temporary space, of size at least k;
	 * @param tmpArray3 temporary space, of size at least 2k.
	 */
	private static final int getKmers_impl(int first, int k, HashMap<Kmer,Kmer> newKmers, HashMap<Kmer,Kmer> oldKmers, int haplotypeCoverage, Kmer key, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) {
		int i;
		int top1, top2, row, column, lastChild, length;
		Kmer value, newKey;
		
		top1=-1; top2=-1;
		tmpArray1[++top1]=-1; tmpArray1[++top1]=-1; tmpArray1[++top1]=first-1;
		while (top1>=0) {
			row=tmpArray1[top1]; column=tmpArray1[top1-1]; lastChild=tmpArray1[top1-2];
			if (row==first+k-1) {
				key.set(tmpArray2,0,k,true); 
				key.canonize(k,tmpArray3);				
				if (newKmers!=null) {
					value=newKmers.get(key);
					if (value==null) {
						newKey=kmerPool_allocate(k);
						newKey.set(key.sequence,0,k,true); 
						newKey.count=1;
						newKmers.put(newKey,newKey);
					}
					else value.count++;
				}
				else {
					value=oldKmers.get(key);
					if (value!=null) return (int)(value.count/haplotypeCoverage);
				}
			}
			if (row==first+k-1 || lastChild==lastInBlock_int[row+1]) { top1-=3; top2--; }
			else {
				lastChild++;
				tmpArray1[top1-2]=lastChild;
				tmpArray1[++top1]=-1; tmpArray1[++top1]=lastChild; tmpArray1[++top1]=row+1;
				tmpArray2[++top2]=intBlocks[row+1][lastChild];
			}
		}
		return -1;
	}
	
	
	/**
	 * For every k-mer in $kmers$, the procedure adds both its canonized prefix and suffix
	 * to $kMinusOneMers$, setting their ${previous,next}Character$ fields.
	 *
	 * @param k length of a k-mer (not of a (k-1)-mer);
	 * @param kMinusOneMer temporary space;
	 * @param tmpArray temporary space of size at least 2(k-1).
	 */
	public static final void getKMinusOneMers(HashMap<Kmer,Kmer> kmers, HashMap<Kmer,Kmer> kMinusOneMers, int k, Kmer kMinusOneMer, int[] tmpArray) {
		boolean sameOrientation;
		int c;
		Kmer kmer, oldKey, newKey;
		Iterator<Kmer> iterator;
		
		kMinusOneMers.clear();
		iterator=kmers.keySet().iterator();
		while (iterator.hasNext()) {
			kmer=iterator.next();
			// Prefix
			System.arraycopy(kmer.sequence,0,kMinusOneMer.sequence,0,kmer.sequence.length-1);
			sameOrientation=kMinusOneMer.canonize(k-1,tmpArray);
			oldKey=kMinusOneMers.get(kMinusOneMer);
			c=canonizeCharacter(kmer.sequence[k-1],sameOrientation);
			if (sameOrientation) {
				if (oldKey==null) {
					newKey = new Kmer(kMinusOneMer,k-1);
					newKey.previousCharacter=-1; newKey.nextCharacter=c;
					kMinusOneMers.put(newKey,newKey);
				}
				else {
					if (oldKey.nextCharacter==-1) oldKey.nextCharacter=c;
					else if (oldKey.nextCharacter>=0 && oldKey.nextCharacter!=c) oldKey.nextCharacter=-2;			
				}
			}
			else {
				if (oldKey==null) {
					newKey = new Kmer(kMinusOneMer,k-1);
					newKey.previousCharacter=c; newKey.nextCharacter=-1;
					kMinusOneMers.put(newKey,newKey);
				}
				else {
					if (oldKey.previousCharacter==-1) oldKey.previousCharacter=c;
					else if (oldKey.previousCharacter>=0 && oldKey.previousCharacter!=c) oldKey.previousCharacter=-2;			
				}
			}
			// Suffix
			System.arraycopy(kmer.sequence,1,kMinusOneMer.sequence,0,kmer.sequence.length-1);
			sameOrientation=kMinusOneMer.canonize(k-1,tmpArray);
			oldKey=kMinusOneMers.get(kMinusOneMer);
			c=canonizeCharacter(kmer.sequence[0],sameOrientation);
			if (sameOrientation) {
				if (oldKey==null) {
					newKey = new Kmer(kMinusOneMer,k-1);
					newKey.previousCharacter=c; newKey.nextCharacter=-1;
					kMinusOneMers.put(newKey,newKey);
				}
				else {
					if (oldKey.previousCharacter==-1) oldKey.previousCharacter=c;
					else if (oldKey.previousCharacter>=0 && oldKey.previousCharacter!=c) oldKey.previousCharacter=-2;			
				}
			}
			else {
				if (oldKey==null) {
					newKey = new Kmer(kMinusOneMer,k-1);
					newKey.previousCharacter=-1; newKey.nextCharacter=c;
					kMinusOneMers.put(newKey,newKey);
				}
				else {
					if (oldKey.nextCharacter==-1) oldKey.nextCharacter=c;
					else if (oldKey.nextCharacter>=0 && oldKey.nextCharacter!=c) oldKey.nextCharacter=-2;			
				}
			}
		}
	}
	
	
	/**
	 * Tries to disambiguate the first and the last block of $str$ using the surrounding
	 * context of length $k$, and appends the updated $str$ to $bw$. This is useful, since
	 * ambiguous endblocks are not included in the unique intervals that are built 
	 * downstream, and thus a suffix-prefix alignment is more likely to be filtered out by
	 * $filterAlignments_*()$ if its endblock is ambiguous.
	 *
	 * Remark: this procedure might disambiguate just few endblocks in practice. This is
	 * expected, since we need a strong signal from a length-k context to proceed. Even
	 * few disambiguated endblocks might have a large effect on the number of filtered
	 * alignments, in theory.
	 *
	 * Remark: we are aware that we are treating blocks that contain several characters
	 * differently, depending on whether they occur in the middle or at the end of a read.
	 * One might want to try and disambiguate blocks in the middle of a read, as well. We
	 * don't do it for simplicity.
	 *
	 * Remark: multiple characters in the same first/last block might canonize to the same
	 * character, since we wrote all the characters that imply the block, and the alphabet
	 * might contain several open/closed versions of the same substring of the same repeat
	 * in the same orientation. If all the characters in a block canonize to the same
	 * character, the procedure selects just one of them (which one is irrelevant for 
	 * computing unique substrings and for filtering alignments downstream).
	 *
	 * Remark: for simplicity the procedure loads all the blocks of $str$, even though 
	 * just the first and the last $k+1$ blocks would suffice.
	 *	
	 * @param kmers a set of k-mers, with their $*Character$ fields correctly set; this
 	 * might be just a subset of all k-mers (e.g. only those with large frequency);
 	 * @param tightMode TRUE=accept a prediction only if every k-mer of $intBlocks[first..
 	 * first+k-1]$ is either adjacent to no character or to a single character, and if
 	 * such a character is the same across all k-mers; FALSE=accept a prediction if all
 	 * k-mers that predict a single character (possibly just one k-mer) agree, even though
 	 * other k-mers might be adjacent to multiple characters;
 	 * @param context temporary space;
 	 * @param tmpArray1 temporary space, of size at least equal to 3 times the number of 
 	 * elements in $intBlocks$;
 	 * @param tmpArray2 temporary space, of size at least k;
 	 * @param tmpArray3 temporary space, of size at least 2k;
	 * @param out adds to cell 0 the number of blocks disambiguated by the procedure 
	 * (0,1,2), and to cell 1 the total number of ambiguous blocks (0,1,2);
	 * @param ambiguityHistogram adds to cell $i$ of row 0 the number of endblocks that 
	 * contain $i>=1$ characters, for every read with at least two blocks; adds to cell 
	 * $i$ of row 1 the number of internal blocks that contain $i>=2$ characters, for 
	 * every read with at least 3 blocks.
	 */
	public static final void fixEndBlocks(String str, int k, HashMap<Kmer,Kmer> kmers, boolean tightMode, Kmer context, int[] boundaries, int readLength, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3, Character tmpChar, BufferedWriter bw, int[] out, int[][] ambiguityHistogram) throws IOException {
		int i, j, c, d, p, q;
		int nBlocks, sum, fixedFirst, fixedLast;
		
		if (str.length()==0) { bw.newLine(); return; }
		
		// Loading blocks
		nBlocks=loadBlocks(str); loadIntBlocks(nBlocks,boundaries,readLength,tmpChar);
		if (IO.CONSISTENCY_CHECKS) {
			for (i=1; i<nBlocks-1; i++) {
				for (j=0; j<=lastInBlock_int[i]; j++) {
					if (alphabet[intBlocks[i][j]].isOpen()) {
						System.err.println("fixEndBlocks> ERROR: interior block "+i+" contains open characters?!");
						for (j=0; j<=lastInBlock_int[i]; j++) System.err.println(alphabet[intBlocks[i][j]]);
						System.err.println("Translated read: "+str);
						System.err.println();
						System.exit(1);
					}
				}
			}
		}
		if (nBlocks<2) {
			bw.write(str); bw.newLine();
			return;
		}
		ambiguityHistogram[0][lastInBlock_int[0]+1>=ambiguityHistogram[0].length?ambiguityHistogram[0].length-1:lastInBlock_int[0]+1]++;
		ambiguityHistogram[0][lastInBlock_int[nBlocks-1]+1>=ambiguityHistogram[0].length?ambiguityHistogram[0].length-1:lastInBlock_int[nBlocks-1]+1]++;
		for (i=1; i<nBlocks-1; i++) {
			if (lastInBlock_int[i]>0) ambiguityHistogram[1][lastInBlock_int[i]+1>=ambiguityHistogram[1].length?ambiguityHistogram[1].length-1:lastInBlock_int[i]+1]++;
		}
		if (lastInBlock_int[0]==0 && lastInBlock_int[nBlocks-1]==0) {
			bw.write(str); bw.newLine();
			return;
		}
		if (lastInBlock_int[0]>0) out[1]++;
		if (lastInBlock_int[nBlocks-1]>0) out[1]++;
		
		// Non-kmer-based fixes
		p=str.indexOf(SEPARATOR_MAJOR+""); q=str.lastIndexOf(SEPARATOR_MAJOR+"");
		fixedFirst=-1; fixedLast=-1;
		if (lastInBlock_int[0]>0) fixedFirst=fixEndBlocks_impl_basic(0);
		if (lastInBlock_int[nBlocks-1]>0) fixedLast=fixEndBlocks_impl_basic(nBlocks-1);
		if (fixedFirst>=0 && fixedLast>=0) {
			out[0]+=2;
			bw.write(fixedFirst+str.substring(p,q+1)+fixedLast); bw.newLine();
			return;
		}
		
		// Kmer-based fixes
		if (nBlocks<k+1) {
			if (fixedFirst>=0) {
				out[0]++;
				bw.write(fixedFirst+str.substring(p)); 
			}
			else if (fixedLast>=0) {
				out[0]++;
				bw.write(str.substring(0,q+1)+fixedLast);
			}
			else bw.write(str);
			bw.newLine();
			return;
		}
		sum=0;
		for (i=0; i<nBlocks; i++) sum+=lastInBlock_int[i]+1;
		if (stack==null || stack.length<sum*3) stack = new int[sum*3];
		if (lastInBlock_int[0]>0) {
			if (fixedFirst>=0) {
				out[0]++;
				bw.write(fixedFirst+str.substring(p,q+1));
			}
			else {
				c=fixEndBlocks_impl_kmer(1,nBlocks,k,kmers,tightMode,context,tmpArray1,tmpArray2,tmpArray3,str);
				if (c>=0) { 
					out[0]++;
					bw.write(c+str.substring(p,q+1)); 
				}
				else bw.write(str.substring(0,q+1));
			}
		}
		else bw.write(str.substring(0,q+1));
		if (lastInBlock_int[nBlocks-1]>0) {
			if (fixedLast>=0) {
				out[0]++;
				bw.write(fixedLast+"");
			}
			else {
				c=fixEndBlocks_impl_kmer(nBlocks-k-1,nBlocks,k,kmers,tightMode,context,tmpArray1,tmpArray2,tmpArray3,str);
				if (c>=0) { 
					out[0]++;
					bw.write(c+"");
				}
				else bw.write(str.substring(q+1));
			}
		}
		else bw.write(str.substring(q+1));
		bw.newLine();
	}
	
	
	/**
	 * If block $blockID$ contains multiple characters, but they are all canonized in the
	 * same way, the procedure selects just one of the most open out of them.
	 *
	 * @return -1 if the block contains distinct characters after canonization; >=0 the 
	 * selected character if the block contains a single character after canonization.
	 */
	private static final int fixEndBlocks_impl_basic(int blockID) {
		int i, c, cPrime;
		int count, canonized, canonizedCount, nonCanonized;
		
		canonized=-1; canonizedCount=-1; nonCanonized=-1;
		for (i=0; i<=lastInBlock_int[blockID]; i++) {
			c=intBlocks[blockID][i]; cPrime=canonizeCharacter(c,true);
			count=(alphabet[c].openStart?1:0)+(alphabet[c].openEnd?1:0);
			if (canonized==-1) { canonized=cPrime; canonizedCount=count; nonCanonized=c; }
			else if (cPrime==canonized) {
				if (count>canonizedCount) { canonizedCount=count; nonCanonized=c; }
			}
			else return -1;
		}
		return nonCanonized;
	}
	
	
	/**
	 * Uses every k-mer that can be built from $intBlocks[first..first+k-1]$ and that 
	 * occurs in $kmers$, as a context for disambiguating the first (if $first=1$) or the
	 * last (if $first=nBlocks-k$) block of $intBlocks$.
	 *
	 * @return -1: impossible to perform a prediction; -2: a prediction is possible, but
	 * the only predicted character cannot be found in the block to disambiguate; >=0: the
	 * character that results from disambiguating the block.
	 */
	private static final int fixEndBlocks_impl_kmer(int first, int nBlocks, int k, HashMap<Kmer,Kmer> kmers, boolean tightMode, Kmer context, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3, String str) {
		boolean sameOrientation;
		int i, c;
		int top1, top2, count, row, column, lastChild, predictedCharacter, selectedCharacter, selectedCount;
		Kmer key;
		
		// Collecting the single character predicted by every possible context
		predictedCharacter=-1;  // In the forward orientation
		top1=-1; top2=-1;
		tmpArray1[++top1]=-1; tmpArray1[++top1]=-1; tmpArray1[++top1]=first-1;
		while (top1>=0) {
			row=tmpArray1[top1]; column=tmpArray1[top1-1]; lastChild=tmpArray1[top1-2];
			if (row==first+k-1) {
				context.set(tmpArray2,0,k,true); 
				sameOrientation=context.canonize(k,tmpArray3);
				key=kmers.get(context);
				if (key!=null) {
					if (first==1) {
						if (sameOrientation) {
							if (tightMode && key.previousCharacter==-2) return -1;
							if (key.previousCharacter>=0) {
								if (predictedCharacter==-1) predictedCharacter=key.previousCharacter;
								else if (predictedCharacter!=key.previousCharacter) return -1;
							}
						}
						else {
							if (tightMode && key.nextCharacter==-2) return -1;
							if (key.nextCharacter>=0) {
								c=canonizeCharacter(key.nextCharacter,false);
								if (predictedCharacter==-1) predictedCharacter=c;
								else if (predictedCharacter!=c) return -1;
							}
						}
					}
					else {
						if (sameOrientation) {
							if (tightMode && key.nextCharacter==-2) return -1;
							if (key.nextCharacter>=0) {
								if (predictedCharacter==-1) predictedCharacter=key.nextCharacter;
								else if (predictedCharacter!=key.nextCharacter) return -1;
							}
						}
						else {
							if (tightMode && key.previousCharacter==-2) return -1;
							if (key.previousCharacter>=0) {
								c=canonizeCharacter(key.previousCharacter,false);
								if (predictedCharacter==-1) predictedCharacter=c;
								else if (predictedCharacter!=c) return -1;
							}
						}
					}
				}
			}
			if (row==first+k-1 || lastChild==lastInBlock_int[row+1]) { top1-=3; top2--; }
			else {
				lastChild++;
				tmpArray1[top1-2]=lastChild;
				tmpArray1[++top1]=-1; tmpArray1[++top1]=lastChild; tmpArray1[++top1]=row+1;
				tmpArray2[++top2]=intBlocks[row+1][lastChild];
			}
		}
		if (predictedCharacter==-1) return -1;
		
		// Choosing the most open among all the characters in the block that canonize to
		// the character predicted by the context.
		selectedCharacter=-1; selectedCount=-1;
		if (first==1) {
			for (i=0; i<=lastInBlock_int[0]; i++) {
				if (canonizeCharacter(intBlocks[0][i],true)==predictedCharacter) {
					count=(alphabet[intBlocks[0][i]].openStart?1:0)+(alphabet[intBlocks[0][i]].openEnd?1:0);
					if (selectedCharacter==-1 || (intBlocks[0][i]!=selectedCharacter && count>selectedCount)) {
						selectedCharacter=intBlocks[0][i];
						selectedCount=count;
					}
				}
			}
		}
		else {
			for (i=0; i<=lastInBlock_int[nBlocks-1]; i++) {
				if (canonizeCharacter(intBlocks[nBlocks-1][i],true)==predictedCharacter) {
					count=(alphabet[intBlocks[nBlocks-1][i]].openStart?1:0)+(alphabet[intBlocks[nBlocks-1][i]].openEnd?1:0);
					if (selectedCharacter==-1 || (intBlocks[nBlocks-1][i]!=selectedCharacter && count>selectedCount)) {
						selectedCharacter=intBlocks[nBlocks-1][i];
						selectedCount=count;
					}
				}
			}
		}
		return selectedCharacter==-1?-2:selectedCharacter;
	}
	

	/**
	 * Writes to $outputFile$ the tandem intervals of every translated read in
	 * $translatedFile$, as a sequence of pairs $firstBlock,lastBlock$.
	 *
	 * Remark: the procedure assumes that array $repeatLength$ has already been loaded.
	 *
	 * @return the total number of tandem intervals found.
	 */
	public static final long getTandemIntervals(String translatedFile, String boundariesFile, String readLengthsFile, String outputFile) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i;
		int nBlocks, lastTandem, readLength;
		long out;
		String str1, str2, str3;
		Character tmpChar;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		
		tmpChar = new Character();
		br1 = new BufferedReader(new FileReader(translatedFile));
		br2 = new BufferedReader(new FileReader(boundariesFile));
		br3 = new BufferedReader(new FileReader(readLengthsFile));
		bw = new BufferedWriter(new FileWriter(outputFile));
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); out=0;
		while (str1!=null) {
			if (str1.length()==0 || str1.indexOf(SEPARATOR_MAJOR+"")<0) {
				bw.newLine();
				str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
				continue;
			}
			nBlocks=loadBlocks(str1); 
			loadBoundaries(str2);
			readLength=Integer.parseInt(str3);
			loadIntBlocks(nBlocks,boundaries,readLength,tmpChar);
			lastTandem=getTandemIntervals_impl(nBlocks,boundaries,readLength,IDENTITY_THRESHOLD,tmpChar);
			if (lastTandem==-1) {
				bw.newLine();
				str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
				continue;
			}
			out+=lastTandem+1;
			bw.write(tandemIntervals[0].from+""+SEPARATOR_MINOR+""+tandemIntervals[0].to);
			for (i=1; i<=lastTandem; i++) bw.write(SEPARATOR_MINOR+""+tandemIntervals[i].from+""+SEPARATOR_MINOR+""+tandemIntervals[i].to);
			bw.newLine();
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
		}
		br1.close(); br2.close(); br3.close(); bw.close();
		return out;
	}
	
	
	/**
	 * Stores in global variable $tandemIntervals$ the set of all maximal tandem intervals
	 * of $intBlocks$, where a tandem interval is a maximal substring such that all its
	 * blocks have a character in common. The procedure searches also for maximal 
	 * intervals that all share an inexact character, with some limited tolerance. Tandem 
	 * intervals that intersect are merged (adjacent tandem intervals are not merged).
	 *
	 * Remark: this definition of tandem is likely too simple for real data, since the
	 * aligner might fail to align some repeat to some tandem unit, or the alignment might
	 * be a bit off and give rise to a different character in our alphabet.
	 *
	 * Remark: the procedure assumes that array $repeatLength$ has already been loaded, 
	 * and it uses global arrays $stack,stack2$.
	 *
	 * @param nBlocks assumed >1;
	 * @param tmpChar temporary space;
	 * @return the last element in $tandemIntervals$.
	 */
	private static final int getTandemIntervals_impl(int nBlocks, int[] boundaries, int readLength, int distanceThreshold, Character tmpChar) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		boolean found, processFirstBlock, processLastBlock;
		int i, j, k, h, c, d;
		int max, last1, last2, lastInterval, tandemLength, from, to, toPrime, length, repeatLength;
		Pair tmpInterval;
		int[] tmpArray;
		
		// Allocating memory
		if (tandemIntervals==null) tandemIntervals = new Pair[100];  // Arbitrary
		max=0;
		for (i=0; i<nBlocks; i++) max=Math.max(max,lastInBlock_int[i]+1);
		if (stack==null || stack.length<max) stack = new int[max];
		if (stack2==null || stack2.length<max) stack2 = new int[max];
		
		// Collecting intervals that share a character
		processFirstBlock=true; processLastBlock=true; lastInterval=-1;
		for (i=0; i<nBlocks-1; i++) {
			System.arraycopy(intBlocks[i],0,stack,0,lastInBlock_int[i]+1);
			last1=lastInBlock_int[i]; tandemLength=1;
			for (j=i+1; j<nBlocks; j++) {
				last2=Math.setIntersection(intBlocks[j],0,lastInBlock_int[j],stack,0,last1,stack2,0);
				if (last2==-1) break;
				tandemLength++;
				tmpArray=stack; stack=stack2; stack2=tmpArray;
				last1=last2;
			}
			if (tandemLength>1) {
				found=false;
				if (i>0) {
					for (j=0; j<=lastInBlock_int[i-1]; j++) {
						if (intBlocks[i-1][j]>lastAlphabet || alphabet[intBlocks[i-1][j]].repeat==-1) continue;
						tmpChar.copyFrom(alphabet[intBlocks[i-1][j]]);
						repeatLength=repeatLengths[tmpChar.repeat];
						length=boundaries[i-1]-(i==1?0:boundaries[i-2]);
						if (intBlocks[i-1][j]<=lastPeriodic) tmpChar.length=length;
						else {
							if (tmpChar.orientation) tmpChar.start=tmpChar.end-length+1>0?tmpChar.end-length+1:0;
							else tmpChar.end=tmpChar.start+length-1<repeatLength-1?tmpChar.start+length-1:repeatLength-1;
						}
						for (k=0; k<=last1; k++) {
							if (tmpChar.isSuffixOf(alphabet[stack[k]],distanceThreshold)) {
								found=true;
								break;
							}
						}
						if (found) break;
					}
				}
				from=found?i-1:i;
				found=false;
				if (i+tandemLength<=nBlocks-1) {
					for (j=0; j<=lastInBlock_int[i+tandemLength]; j++) {
						if (intBlocks[i+tandemLength][j]>lastAlphabet || alphabet[intBlocks[i+tandemLength][j]].repeat==-1) continue;
						tmpChar.copyFrom(alphabet[intBlocks[i+tandemLength][j]]);
						tmpChar.openEnd=true;
						repeatLength=repeatLengths[tmpChar.repeat];
						length=(i+tandemLength==nBlocks-1?readLength:boundaries[i+tandemLength])-boundaries[i+tandemLength-1];
						if (intBlocks[i+tandemLength][j]<=lastPeriodic) tmpChar.length=length;
						else {
							if (tmpChar.orientation) tmpChar.end=tmpChar.start+length-1<repeatLength-1?tmpChar.start+length-1:repeatLength-1;
							else tmpChar.start=tmpChar.end-length+1>0?tmpChar.end-length+1:0;
						}
						for (k=0; k<=last1; k++) {
							if (tmpChar.isPrefixOf(alphabet[stack[k]],distanceThreshold)) {
								found=true;
								break;
							}
						}
						if (found) break;
					}
				}
				to=found?i+tandemLength:i+tandemLength-1;
				lastInterval++;
				if (lastInterval>=tandemIntervals.length) {
					Pair[] newArray = new Pair[tandemIntervals.length<<1];
					System.arraycopy(tandemIntervals,0,newArray,0,tandemIntervals.length);
					tandemIntervals=newArray;
				}
				if (tandemIntervals[lastInterval]==null) tandemIntervals[lastInterval] = new Pair(from,to);
				else tandemIntervals[lastInterval].set(from,to);
				if (from<=1) processFirstBlock=false;
				if (to>=nBlocks-2) processLastBlock=false;
			}
		}
		
		// Adding interval $[0..1]$ if it is a tandem.
		if (processFirstBlock) {
			found=false;
			for (i=0; i<=lastInBlock_int[0]; i++) {
				if (intBlocks[0][i]>lastAlphabet || alphabet[intBlocks[0][i]].repeat==-1) continue;
				tmpChar.copyFrom(alphabet[intBlocks[0][i]]);
				repeatLength=repeatLengths[tmpChar.repeat];
				length=boundaries[0];
				if (intBlocks[0][i]<=lastPeriodic) tmpChar.length=length;
				else {
					if (tmpChar.orientation) tmpChar.start=tmpChar.end-length+1>0?tmpChar.end-length+1:0;
					else tmpChar.end=tmpChar.start+length-1<repeatLength-1?tmpChar.start+length-1:repeatLength-1;
				}
				for (j=0; j<=lastInBlock_int[1]; j++) {
					if (intBlocks[1][j]>lastAlphabet || alphabet[intBlocks[1][j]].repeat==-1) continue;
					if (tmpChar.isSuffixOf(alphabet[intBlocks[1][j]],distanceThreshold)) {
						found=true;
						break;
					}
				}
				if (found) break;
			}
			if (found) {
				lastInterval++;
				if (lastInterval>=tandemIntervals.length) {
					Pair[] newArray = new Pair[tandemIntervals.length<<1];
					System.arraycopy(tandemIntervals,0,newArray,0,tandemIntervals.length);
					tandemIntervals=newArray;
				}
				if (tandemIntervals[lastInterval]==null) tandemIntervals[lastInterval] = new Pair(0,1);
				else tandemIntervals[lastInterval].set(0,1);
			}
		}
	
		// Adding interval $[nBlocks-2..nBlocks-1]$ if it is a tandem.
		if (processLastBlock) {
			found=false;
			for (i=0; i<=lastInBlock_int[nBlocks-1]; i++) {
				if (intBlocks[nBlocks-1][i]>lastAlphabet || alphabet[intBlocks[nBlocks-1][i]].repeat==-1) continue;
				tmpChar.copyFrom(alphabet[intBlocks[nBlocks-1][i]]);
				repeatLength=repeatLengths[tmpChar.repeat];
				length=readLength-boundaries[nBlocks-2];
				if (intBlocks[nBlocks-1][i]<=lastPeriodic) tmpChar.length=length;
				else {
					if (tmpChar.orientation) tmpChar.end=tmpChar.start+length-1<repeatLength-1?tmpChar.start+length-1:repeatLength-1;
					else tmpChar.start=tmpChar.end-length+1>0?tmpChar.end-length+1:0;
				}
				for (j=0; j<=lastInBlock_int[nBlocks-2]; j++) {
					if (intBlocks[nBlocks-2][j]>lastAlphabet || alphabet[intBlocks[nBlocks-2][j]].repeat==-1) continue;
					if (tmpChar.isPrefixOf(alphabet[intBlocks[nBlocks-2][j]],distanceThreshold)) {
						found=true;
						break;
					}
				}
				if (found) break;
			}
			if (found) {
				lastInterval++;
				if (lastInterval>=tandemIntervals.length) {
					Pair[] newArray = new Pair[tandemIntervals.length<<1];
					System.arraycopy(tandemIntervals,0,newArray,0,tandemIntervals.length);
					tandemIntervals=newArray;
				}
				if (tandemIntervals[lastInterval]==null) tandemIntervals[lastInterval] = new Pair(nBlocks-2,nBlocks-1);
				else tandemIntervals[lastInterval].set(nBlocks-2,nBlocks-1);
			}
		}
		
		// Collecting intervals that share similar characters
		for (i=0; i<nBlocks-1; i++) {
			to=0;
			for (j=0; j<=lastInBlock_int[i]; j++) {
				c=intBlocks[i][j];
				if (c<=lastUnique || c>lastAlphabet) continue;
				toPrime=0;
				for (k=i+1; k<=nBlocks-1; k++) {
					found=false;
					for (h=0; h<=lastInBlock_int[k]; h++) {
						d=intBlocks[k][h];
						if (d<=lastUnique || d>lastAlphabet) continue;
						if ( alphabet[d].repeat==alphabet[c].repeat && alphabet[d].orientation==alphabet[c].orientation &&
							 ( (alphabet[c].repeat<=lastPeriodic && Math.abs(alphabet[c].length,alphabet[d].length)<=IDENTITY_THRESHOLD) ||
							   ( alphabet[c].repeat>lastPeriodic && 
								 ( (alphabet[c].start==alphabet[d].start && Math.abs(alphabet[c].end,alphabet[d].end)<=IDENTITY_THRESHOLD) ||
								   (alphabet[c].end==alphabet[d].end && Math.abs(alphabet[c].start,alphabet[d].start)<=IDENTITY_THRESHOLD)
								 )
							   )
							 )
						   ) {
							found=true;
							break;
						}
					}
					if (!found) break;
					else toPrime=k;
				}
				to=Math.max(to,toPrime);
			}
			if (to!=0) {
				lastInterval++;
				if (lastInterval>=tandemIntervals.length) {
					Pair[] newArray = new Pair[tandemIntervals.length<<1];
					System.arraycopy(tandemIntervals,0,newArray,0,tandemIntervals.length);
					tandemIntervals=newArray;
				}
				if (tandemIntervals[lastInterval]==null) tandemIntervals[lastInterval] = new Pair(i,to);
				else tandemIntervals[lastInterval].set(i,to);
			}
		}
		
		// Merging overlapping intervals
		if (lastInterval<=0) return lastInterval;
		Arrays.sort(tandemIntervals,0,lastInterval+1);
		j=0;
		for (i=0; i<=lastInterval; i++) {
			if (tandemIntervals[i].from<=tandemIntervals[j].to) tandemIntervals[j].to=tandemIntervals[i].to;
			else {
				j++;
				tmpInterval=tandemIntervals[j]; tandemIntervals[j]=tandemIntervals[i]; tandemIntervals[i]=tmpInterval;
			}
		}
		return j;
	}
	
	
	private static class Pair implements Comparable {
		public int from, to;
		
		public Pair(int f, int t) {
			set(f,t);
		}
		
		public void set(int f, int t) {
			this.from=f; this.to=t;
		}
		
		public int compareTo(Object other) {
			Pair otherPair = (Pair)other;
			if (from<otherPair.from) return -1;
			else if (from>otherPair.from) return 1;
			return 0;
		}
	}
	
	
	/**
	 * Initializes global data structures $tandems,lastTandem$ from the content of 
	 * $inputFile$.
	 */
	public static final void loadTandemIntervals(String inputFile, int nReads) throws IOException {
		int i, j;
		int length;
		String str;
		BufferedReader br;
		String[] tokens;
		
		tandems = new int[nReads][0];
		lastTandem = new int[nReads];
		Math.set(lastTandem,nReads-1,-1);
		br = new BufferedReader(new FileReader(inputFile));
		str=br.readLine(); i=-1;
		while (str!=null) {
			i++;
			if (str.length()==0) {
				str=br.readLine();
				continue;
			}
			tokens=str.split(SEPARATOR_MINOR+"");
			length=tokens.length;
			tandems[i] = new int[length];
			for (j=0; j<length; j++) tandems[i][j]=Integer.parseInt(tokens[j]);
			lastTandem[i]=length-1;
			str=br.readLine();
		}
		br.close();
	}
	
	
	/**
	 * Returns the smallest position in $alphabet$ of a character with the same
	 * repeat, start/end (or length), and orientation, as $alphabet[characterID]$. This is
	 * done to merge all instances of the same substring of the same repeat, regardless of
	 * the open status of their endpoints.
	 *
	 * @param sameOrientation TRUE=keeps the same orientation; FALSE=reverse-
	 * complements the orientation of the character.
	 */
	private static final int canonizeCharacter(int characterID, boolean sameOrientation) {
		boolean orientation;
		int i;
		int last, repeat, length, start, end;
		
		if (characterID==lastAlphabet+1 || characterID<=lastUnique) {
			// Unique characters have no orientation and are all closed
			return characterID;
		}
		else if (characterID<=lastPeriodic) {
			repeat=alphabet[characterID].repeat;
			orientation=alphabet[characterID].orientation;
			length=alphabet[characterID].length;
			last=-1;
			if (orientation) {
				if (sameOrientation) {
					for (i=characterID-1; i>lastUnique; i--) {
						if (alphabet[i].repeat!=repeat) break;
						if (alphabet[i].length==length) last=i;
					}
					return last==-1?characterID:last;
				}
				else {
					for (i=characterID+1; i<=lastPeriodic; i++) {
						if (alphabet[i].repeat!=repeat) break;
						if (!alphabet[i].orientation && alphabet[i].length==length) {
							last=i;
							break;
						}
					}
					if (last==-1) {
						System.err.println("canonizeCharacter> ERROR: reverse-complement of the following character not found in the alphabet: "+alphabet[characterID]);
						System.exit(1);
					}
					return last;
				}
			}
			else {
				if (sameOrientation) {
					for (i=characterID-1; i>lastUnique; i--) {
						if (alphabet[i].repeat!=repeat || alphabet[i].orientation) break;
						if (alphabet[i].length==length) last=i;
					}
					return last==-1?characterID:last;
				}
				else {
					for (i=characterID-1; i>lastUnique; i--) {
						if (alphabet[i].repeat!=repeat) break;
						if (alphabet[i].orientation && alphabet[i].length==length) last=i;
					}
					if (last==-1) {
						System.err.println("canonizeCharacter> ERROR: reverse-complement of the following character not found in the alphabet: "+alphabet[characterID]);
						System.exit(1);
					}
					return last;
				}
			}
		}
		else {
			repeat=alphabet[characterID].repeat;
			orientation=alphabet[characterID].orientation;
			start=alphabet[characterID].start; end=alphabet[characterID].end;
			last=-1;
			if (orientation) {
				if (sameOrientation) {
					for (i=characterID-1; i>lastPeriodic; i--) {
						if (alphabet[i].repeat!=repeat) break;
						if (alphabet[i].start==start && alphabet[i].end==end) last=i;
					}
					return last==-1?characterID:last;
				}
				else {
					for (i=characterID+1; i<=lastAlphabet; i++) {
						if (alphabet[i].repeat!=repeat) break;
						if (!alphabet[i].orientation && alphabet[i].start==start && alphabet[i].end==end) {
							last=i;
							break;
						}
					}
					if (last==-1) {
						System.err.println("canonizeCharacter> ERROR: reverse-complement of the following character not found in the alphabet: "+alphabet[characterID]);
						System.exit(1);
					}
					return last;
				}
			}
			else {
				if (sameOrientation) {
					for (i=characterID-1; i>lastPeriodic; i--) {
						if (alphabet[i].repeat!=repeat || alphabet[i].orientation) break;
						if (alphabet[i].start==start && alphabet[i].end==end) last=i;
					}
					return last==-1?characterID:last;
				}
				else {
					for (i=characterID-1; i>lastPeriodic; i--) {
						if (alphabet[i].repeat!=repeat) break;
						if (alphabet[i].orientation && alphabet[i].start==start && alphabet[i].end==end) last=i;
					}
					if (last==-1) {
						System.err.println("canonizeCharacter> ERROR: reverse-complement of the following character not found in the alphabet: "+alphabet[characterID]);
						System.exit(1);
					}
					return last;
				}
			}
		}
	}
	
	
	public static class Kmer implements Comparable {
		public int[] sequence;  // Positions in $alphabet$.
		public long count;
		public int sameReadCount;  // Max frequency inside a read
		
		/**
		 * -1: unknown; -2: more than one character. >=0: canonized index in $alphabet$
		 * of the only character that precedes/follows in the recoded reads.
		 */
		public int previousCharacter, nextCharacter;
		
		public Kmer() { }
		
		public Kmer(Kmer otherKmer, int k) {
			set(otherKmer.sequence,0,k,false);
		}
		
		public Kmer(String str, int k) {
			sequence = new int[k];
			String[] tokens = str.split(",");
			for (int i=0; i<k; i++) sequence[i]=Integer.parseInt(tokens[i]);
			count=Long.parseLong(tokens[k]);
			if (tokens.length>k+1) {
				previousCharacter=Integer.parseInt(tokens[k+1]);
				nextCharacter=Integer.parseInt(tokens[k+2]);
			}
			else { previousCharacter=-1; nextCharacter=-1; }
		}
		
		/**
		 * @param reuseSequence TRUE=reuses the existing array if long enough.
		 */
		public void set(int[] fromArray, int first, int k, boolean reuseSequence) {
			if (sequence==null || sequence.length<k || !reuseSequence) sequence = new int[k];
			System.arraycopy(fromArray,first,sequence,0,k);
			count=0; previousCharacter=-1; nextCharacter=-1;
		}
		
		/**
		 * Resets $sequence$ to the lexicographically smaller between $sequence$ with
		 * every character canonized, and the reverse of $sequence$ with the complement of 
		 * every character canonized.
		 * 
		 * @tmpArray temporary space of size at least 2k;
		 * @return TRUE iff canonization did not change the orientation of $sequence$.
		 */
		public boolean canonize(int k, int[] tmpArray) {
			boolean smaller;
			int i;
			
			for (i=0; i<k; i++) tmpArray[i]=canonizeCharacter(sequence[i],true);
			for (i=0; i<k; i++) tmpArray[k+i]=canonizeCharacter(sequence[k-1-i],false);
			smaller=true;
			for (i=0; i<k; i++) {
				if (tmpArray[i]<tmpArray[k+i]) break;
				else if (tmpArray[i]>tmpArray[k+i]) {
					smaller=false; break;
				}
			}
			System.arraycopy(tmpArray,smaller?0:k,sequence,0,k);
			return smaller;
		}
		
		public boolean equals(Object other) {
			Kmer otherKmer = (Kmer)other;
			final int k = sequence.length;
			for (int i=0; i<k; i++) {
				if (sequence[i]!=otherKmer.sequence[i]) return false;
			}
			return true;
		}
		
		public int compareTo(Object other) {
			Kmer otherKmer = (Kmer)other;
			final int k = sequence.length;
			for (int i=0; i<k; i++) {
				if (sequence[i]<otherKmer.sequence[i]) return -1;
				else if (sequence[i]>otherKmer.sequence[i]) return 1;
			}
			return 0;
		}
		
		public int hashCode() {
			return Arrays.hashCode(sequence);
		}
		
		public String toString() {
			final int k = sequence.length;
			String out = sequence[0]+"";
			for (int i=1; i<k; i++) out+=","+sequence[i];
			return out;
		}
	}
	
	
	
	
	// ------------------------- UNIQUE INTERVALS PROCEDURES -----------------------------
	
	/**
	 * Initializes global variables $fullyUnique,fullyContained$ to the sorted list of IDs
	 * of reads that are fully nonrepetitive and fully contained in a single repeat block,
	 * respectively.
	 */
	public static final void loadReadsFully(String fullyUniqueFile, int nFullyUnique, String fullyContainedFile, int nFullyContained) throws IOException {
		int i;
		BufferedReader br;
		
		fullyUnique = new int[nFullyUnique];
		br = new BufferedReader(new FileReader(fullyUniqueFile));
		for (i=0; i<nFullyUnique; i++) fullyUnique[i]=Integer.parseInt(br.readLine());
		br.close();
		fullyContained = new int[nFullyContained];
		br = new BufferedReader(new FileReader(fullyContainedFile));
		for (i=0; i<nFullyContained; i++) fullyContained[i]=Integer.parseInt(br.readLine());
		br.close();
	}
	
	
	/**
	 * Loads $translated$ and, for every read in it, the row of the boundaries file in
	 * $boundaries_all$.
	 *
	 * Remark: the procedure assumes that the read length of every read has already been 
	 * loaded in $Reads.readLengths$.
	 *
	 * @param loadIsBlockUnique loads the bitvector $isBlockUnique_all$ for every read in 
	 * $translated$;
	 * @param loadTranslation loads in $translation_all$ the translation of every read 
	 * in $translated$;
	 * @return the max length of a block in $translation_all$, if this data structure is 
	 * loaded; not defined otherwise.
	 */
	public static final int loadAllBoundaries(String translatedFile, boolean loadIsBlockUnique, boolean loadTranslation, String boundariesFile) throws IOException {
		int i, j, r;
		int read, cell, mask, block, nBlocks, nBytes, nReads, nBoundaries, readLength, max;
		String str;
		Character tmpChar;
		BufferedReader br;
		String[] tokens;
		
		// Computing the set of translated reads that are not fully contained in a repeat
		br = new BufferedReader(new FileReader(translatedFile));
		str=br.readLine();
		nReads=0;
		while (str!=null) {
			if (str.length()!=0 && str.indexOf(SEPARATOR_MAJOR+"")>=0) nReads++;
			str=br.readLine();
		}
		br.close();
		translated = new int[nReads];
		br = new BufferedReader(new FileReader(translatedFile));
		str=br.readLine(); r=0; i=-1;
		while (str!=null) {
			if (str.length()!=0 && str.indexOf(SEPARATOR_MAJOR+"")>=0) translated[++i]=r;
			str=br.readLine(); r++;
		}
		br.close();
		
		// Loading $boundaries_all$.
		boundaries_all = new int[nReads][0];
		br = new BufferedReader(new FileReader(boundariesFile));
		str=br.readLine(); r=-1;
		while (str!=null) {
			if (str.length()==0) {
				str=br.readLine();
				continue;
			}
			r++;
			if (str.indexOf(SEPARATOR_MINOR+"")>=0) {
				tokens=str.split(SEPARATOR_MINOR+"");
				nBoundaries=tokens.length;
				boundaries_all[r] = new int[nBoundaries];
				for (i=0; i<nBoundaries; i++) boundaries_all[r][i]=Integer.parseInt(tokens[i]);
			}
			else {
				boundaries_all[r] = new int[1];
				boundaries_all[r][0]=Integer.parseInt(str);
			}
			str=br.readLine();
		}
		br.close();
		
		// Loading $isBlockUnique_all$.
		tmpChar=null;
		if (loadIsBlockUnique) {
			tmpChar = new Character();
			isBlockUnique_all = new byte[nReads][0];
			br = new BufferedReader(new FileReader(translatedFile));
			str=br.readLine(); r=-1;
			while (str!=null) {
				if (str.length()==0 || str.indexOf(SEPARATOR_MAJOR+"")<0) {
					str=br.readLine();
					continue;
				}
				r++;
				nBlocks=loadBlocks(str);
				loadIntBlocks(nBlocks,boundaries_all[r],Reads.getReadLength(translated[r]),tmpChar);
				nBytes=Math.ceil(nBlocks,8);
				isBlockUnique_all[r] = new byte[nBytes];
				for (i=0; i<nBytes; i++) {
					cell=0; mask=1;
					for (j=0; j<8; j++) {
						block=i*8+j;
						if (block==nBlocks) break;
						if (isBlockUnique[block]) cell|=mask;
						mask<<=1;
					}
					isBlockUnique_all[r][i]=(byte)cell;
				}
				str=br.readLine();
			}
			br.close();
		}
		
		// Loading $translation_all$.
		max=0;
		if (loadTranslation) {
			if (tmpChar==null) tmpChar = new Character();
			translation_all = new int[nReads][0][0];
			br = new BufferedReader(new FileReader(translatedFile));
			str=br.readLine(); r=-1;
			while (str!=null) {
				if (str.length()==0 || str.indexOf(SEPARATOR_MAJOR+"")<0) {
					str=br.readLine();
					continue;
				}
				r++;
				nBlocks=loadBlocks(str);
				loadIntBlocks(nBlocks,boundaries_all[r],Reads.getReadLength(translated[r]),tmpChar);
				translation_all[r] = new int[nBlocks][0];
				for (i=0; i<nBlocks; i++) {
					translation_all[r][i] = new int[lastInBlock_int[i]+1];
					System.arraycopy(intBlocks[i],0,translation_all[r][i],0,lastInBlock_int[i]+1);
					max=Math.max(max,lastInBlock_int[i]+1);
				}
				str=br.readLine();
			}
			br.close();
		}
		return max;
	}
	
	
	/**
	 * Loads all the \emph{blue intervals}, i.e. sequences of repeat and unique characters
	 * that are likely to occur just once in the genome, from all reads that contain one.
	 */
	public static final void loadBlueIntervals(String intervalsFile) throws IOException {
		final int GROWTH_RATE = 1000;  // Arbitrary
		int i;
		int read;
		String str;
		BufferedReader br;
		String[] tokens;
		
		blueIntervals = new int[GROWTH_RATE][0];
		blueIntervals_reads = new int[GROWTH_RATE];
		br = new BufferedReader(new FileReader(intervalsFile));
		str=br.readLine(); read=-1; blueIntervals_last=-1;
		while (str!=null) {
			read++;
			if (str.length()==0) {
				str=br.readLine();
				continue;
			}
			blueIntervals_last++;
			if (blueIntervals_last==blueIntervals.length) {
				int[][] newArray = new int[blueIntervals.length+GROWTH_RATE][0];
				System.arraycopy(blueIntervals,0,newArray,0,blueIntervals.length);
				blueIntervals=newArray;
				int[] newArray2 = new int[blueIntervals.length+GROWTH_RATE];
				System.arraycopy(blueIntervals_reads,0,newArray2,0,blueIntervals_reads.length);
				blueIntervals_reads=newArray2;
			}
			blueIntervals_reads[blueIntervals_last]=read;
			tokens=str.split(SEPARATOR_MINOR+"");
			blueIntervals[blueIntervals_last] = new int[tokens.length];
			for (i=0; i<tokens.length; i++) blueIntervals[blueIntervals_last][i]=Integer.parseInt(tokens[i]);
			str=br.readLine();
		}
		br.close();
	}
	
	
	/**
	 * Writes to $outputFile$ a zero for every alignment of $alignmentsFile$ that belongs 
	 * to a \emph{red region} on both readA and readB, i.e. to a region that fully belongs
	 * to a repeat character, or to a sequence of repeat characters, that is likely to 
	 * have more than H occurrences in the genome, where H is the number of haplotypes.
	 * The intervals of the alignment in the two reads might cover mismatching sequences 
	 * of boundaries and different characters, but we can safely discard the alignment 
	 * anyway, since it just encodes a similarity between substrings of (possibly 
	 * different) repeats.
	 *
	 * Remark: the ad hoc handling of suffix-prefix overlaps in $filterAlignments_tight()$
	 * is not needed here, since we keep an alignment just if it contains a unique 
	 * interval in either read.
	 *
	 * Remark: every alignment between intervals at frequency <=H should be kept in the
	 * assembly graph. In heterozygous regions such alignments might create multiple
	 * alternative paths, and these should be resolved downstream. One might label every
	 * alignment with the estimated number of haplotypes it involves, as a clue for the
	 * assembler, but this might require extracting many more shortest unique intervals
	 * (since an interval that occurs once in X haplotypes might be contained in an 
	 * interval that occurs once in < X haplotypes). We don't do this for simplicity.
	 *
	 * Remark: the procedure filters out some alignments even when the less stringent 
	 * settings possible are used. Specifically, it always removes alignments that are 
	 * fully inside a repeat in both reads (or in just one read for the $_tight$ version),
	 * and alignments such that all their contained k-mers (k>1) start or end with a 
	 * unique block that is too short to be considered a unique signature.
	 *
	 * Remark: the procedure does not need $alphabet$, but it needs the following arrays: 
	 * isBlockUnique_all | fullyUnique, fullyContained, boundaries_all, blueIntervals, 
	 * blueIntervals_reads. 
	 *
	 * @param alignmentsFile output of LAshow, assumed to be sorted by readA;
	 * @param minIntersection_nonrepetitive min. length of a non-repetitive substring of 
	 * the alignment, for the alignment not to be considered red; this should not be too 
	 * small, since short non-repetitive regions might not address a unique locus of the 
	 * genome, and they might even be occurrences of short repeats that were not aligned 
	 * to the repeat database because of heuristics of the aligner;
	 * @param out output array containing the number of alignments for each type (columns)
	 * specified in $Alignments.readAlignmentFile_getType()$; row 0: all alignments in 
	 * input; row 1: all alignments kept in output; row 2: only input alignments that
	 * intersect a non-repetitive region (these are kept in output).
	 */
	public static final void filterAlignments_loose(String alignmentsFile, String outputFile, int minIntersection_nonrepetitive, int minBlueIntervalLength, long[][] out) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int p;
		int row, readA, readB, type, value, isRepetitive;
		int lastFullyContained, lastFullyUnique, lastTranslated, lastBlueInterval;
		final int nFullyContained = fullyContained.length;
		final int nFullyUnique = fullyUnique.length;
		final int nTranslated = translated.length;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(outputFile));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); row=0; 
		lastFullyContained=0; lastFullyUnique=0; lastTranslated=0; lastBlueInterval=0;
		while (str!=null)  {
			if (row%1000000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			type=Alignments.readAlignmentFile_getType(IDENTITY_THRESHOLD);
			out[0][type]++;	
			// Processing readA
			readA=Alignments.readA-1;
			while (lastFullyUnique<nFullyUnique && fullyUnique[lastFullyUnique]<readA) lastFullyUnique++;
			if (lastFullyUnique<nFullyUnique && fullyUnique[lastFullyUnique]==readA) {
				bw.write("1\n"); str=br.readLine(); row++;
				out[1][type]++; out[2][type]++;
				continue;
			}
			isRepetitive=-3;
			while (lastFullyContained<nFullyContained && fullyContained[lastFullyContained]<readA) lastFullyContained++;
			if (lastFullyContained<nFullyContained && fullyContained[lastFullyContained]==readA) isRepetitive=0;
			else {
				while (lastTranslated<nTranslated && translated[lastTranslated]<readA) lastTranslated++;
				while (lastBlueInterval<=blueIntervals_last && blueIntervals_reads[lastBlueInterval]<readA) lastBlueInterval++;
				isRepetitive=inRedRegion(readA,Alignments.startA,Alignments.endA,lastTranslated,lastBlueInterval<=blueIntervals_last&&blueIntervals_reads[lastBlueInterval]==readA?lastBlueInterval:-1,-1,Reads.getReadLength(readA),minIntersection_nonrepetitive,minBlueIntervalLength);
			}
			if (isRepetitive!=0) {
				bw.write("1\n"); str=br.readLine(); row++;
				out[1][type]++; 
				if (isRepetitive==-1) out[2][type]++;
				continue;
			}
			// Processing readB
			readB=Alignments.readB-1;
			if (readInArray(readB,fullyUnique,nFullyUnique-1,lastFullyUnique)>=0) {
				bw.write("1\n"); str=br.readLine(); row++;
				out[1][type]++; out[2][type]++;
				continue;
			}
			else if (readInArray(readB,fullyContained,nFullyContained-1,lastFullyContained)>=0) {
				bw.write("0\n"); str=br.readLine(); row++;
				continue;
			}
			p=readInArray(readB,translated,nTranslated-1,lastTranslated);
			if (p<0) {
				System.err.println("filterAlignments_loose> ERROR: read "+readB+" is neither translated, nor fully unique, nor fully contained in a repeat.");
				System.exit(1);
			}
			isRepetitive=inRedRegion(readB,Alignments.startB,Alignments.endB,p,-2,lastBlueInterval,Reads.getReadLength(readB),minIntersection_nonrepetitive,minBlueIntervalLength);
			if (isRepetitive==0) bw.write("0\n"); 
			else {
				bw.write("1\n");
				out[1][type]++;
				if (isRepetitive==-1) out[2][type]++;
			}
			str=br.readLine(); row++;
		}
		br.close(); bw.close();
	}
	
	
	/**
	 * Tells whether interval $readID[intervalStart..intervalEnd]$ fully belongs to a 
	 * single repeat character, or to a sequence of repeat characters, that is likely to
	 * have more than H occurrences in the genome, where H is the number of haplotypes.
	 *
	 * Remark: the procedure needs the following arrays: 
	 * isBlockUnique_all, boundaries_all, blueIntervals, blueIntervals_reads.
	 * 
	 * @param readID assumed to contain more than one block;
	 * @param boundariesAllID position of the read in $boundaries_all,isBlockUnique_all,
	 * translation_all$;
	 * @param blueIntervalsID position of the read in $blueIntervals$; or -1 if the 
	 * read does not occur in $blueIntervals$; or -2 if the ID of the read in 
	 * $blueIntervals$ is unknown;
	 * @param blueIntervalsStart a position in $blueIntervals$ from which to start
	 * the search when $blueIntervalsID=-2$;
	 * @param minIntersection_nonrepetitive min. length of a non-repetitive substring of 
	 * the alignment, for the alignment to be considered non-repetitive;
	 * @return interval $readID[intervalStart..intervalEnd]$:
	 * -1: belongs to or straddles a non-repetitive region;
	 * -2: contains a sequence of $>=minBlueIntervalLength$ repeat characters that likely
	 *     occurs <=H times in the genome;
	 *  0: none of the above is true.
	 */
	private static final int inRedRegion(int readID, int intervalStart, int intervalEnd, int boundariesAllID, int blueIntervalsID, int blueIntervalsStart, int readLength, int minIntersection_nonrepetitive, int minBlueIntervalLength) {
		int i, j;
		int mask, cell, start, end, blockStart, blockEnd, firstBlock, lastBlock;
		final int nBlocks = boundaries_all[boundariesAllID].length+1;
		final int nBytes = Math.ceil(nBlocks,8);
		
		// Checking the nonrepetitive blocks of the read, if any.
		cell=isBlockUnique_all[boundariesAllID][0]; 
		if (cell!=0) {
			mask=1; blockStart=0; blockEnd=boundaries_all[boundariesAllID][0];
			if ( (cell&mask)!=0 && 
				 Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive &&
				 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd)
				 )
			   ) return -1;
			mask<<=1;
			for (j=1; j<8; j++) {
				if (j==nBlocks-1) break;
				blockStart=boundaries_all[boundariesAllID][j-1];
				blockEnd=boundaries_all[boundariesAllID][j];
				if ( (cell&mask)!=0 && 
					 Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive &&
					 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) || 
					   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
					   !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd)
					 )
				   ) return -1;
				mask<<=1;
			}
		}
		for (i=1; i<nBytes; i++) {
			cell=isBlockUnique_all[boundariesAllID][i]; 
			if (cell==0) continue;
			mask=1;
			for (j=0; j<8; j++) {
				if ((i<<3)+j==nBlocks-1) break;
				blockStart=boundaries_all[boundariesAllID][(i<<3)+j-1];
				blockEnd=boundaries_all[boundariesAllID][(i<<3)+j];
				if ( (cell&mask)!=0 && 
					 Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive &&
					 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
					   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
					   !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd)
					 )
				   ) return -1;
				mask<<=1;
			}
		}
		cell=isBlockUnique_all[boundariesAllID][nBytes-1]; 
		if (cell!=0) {
			mask=1<<((nBlocks%8)-1);
			blockStart=boundaries_all[boundariesAllID][nBlocks-2];
			blockEnd=readLength-1;
			if ( (cell&mask)!=0 &&
				 Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive &&
				 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd)
				 )
			   ) return -1;
		}
		
		// Checking the repetitive blocks of the read, if any.
		if (blueIntervalsID==-2) blueIntervalsID=readInArray(readID,blueIntervals_reads,blueIntervals_reads.length-1,blueIntervalsStart);
		if (blueIntervalsID!=-1) {
			for (i=0; i<blueIntervals[blueIntervalsID].length; i+=3) {
				if (blueIntervals[blueIntervalsID][i+1]<minBlueIntervalLength) continue;		
				firstBlock=blueIntervals[blueIntervalsID][i];
				start=firstBlock==0?0:boundaries_all[boundariesAllID][firstBlock-1];
				lastBlock=firstBlock+blueIntervals[blueIntervalsID][i+1]-1;
				end=lastBlock==nBlocks-1?readLength-1:boundaries_all[boundariesAllID][lastBlock];
				if ( Intervals.areApproximatelyIdentical(start,end,intervalStart,intervalEnd) ||
					 Intervals.isApproximatelyContained(start,end,intervalStart,intervalEnd)
				   ) return -2;
			}
		}
		
		return 0;
	}
	
	
	/**
	 * @param array sorted;
	 * @param position an arbitrary position in $array$ from which to start the search
	 * (might be greater than $last$);
	 * @return the position of $read$ in $array[0..last]$, or -1 if it does not occur.
	 */
	private static final int readInArray(int read, int[] array, int last, int position) {
		final int MAX_DIFF = 100;  // Arbitrary
		int i;
		int value;
		
		if (position<=last) {
			value=array[position];
			if (value==read) return position;
			else if (read<value) {
				if (read>=value-MAX_DIFF) {
					for (i=position-1; i>=0; i--) {
						if (array[i]<read) break;
						if (array[i]==read) return i;
					}
				}
				else {
					i=Arrays.binarySearch(array,0,position,read);
					if (i>=0) return i;
				}
			}
			else {
				if (read<=value+MAX_DIFF) {
					for (i=position+1; i<=last; i++) {
						if (array[i]>read) break;
						if (array[i]==read) return i;
					}
				}
				else {
					i=Arrays.binarySearch(array,position+1,last+1,read);
					if (i>=0) return i;
				}
			}
		}
		else {
			i=Arrays.binarySearch(array,0,last+1,read);
			if (i>=0) return i;
		}
		return -1;
	}
	
	
	/**
	 * Writes to $outputFile$ a one for every alignment of $alignmentsFile$ that contains,
	 * in both readA and readB, a region that is likely to have at most H occurrences in
	 * the genome, where H is the number of haplotypes. If $mode=TRUE$, the procedure 
	 * additionally requires that the intervals of the alignment in the two reads cover
	 * matching sequences of boundaries and matching characters.
	 *
	 * Remark: the procedure needs the following arrays: 
	 * translation_all, alphabet | fullyUnique, fullyContained, boundaries_all, 
	 * blueIntervals, blueIntervals_reads.
	 *
	 * @param suffixPrefixMode TRUE = a suffix-prefix alignment is kept even if it 
	 * contains a unique interval in just one read; this is allowed when the alignment
	 * straddles a unique interval on the other read, on the side that is opposite to the
	 * end of the read. Alignments kept by this criterion are mostly noisy in practice, so
	 * it seems safe to disable it.
	 * @param out row 2 stores the number of input alignments that overlap a non-
	 * repetitive region on both reads (these are kept in output).
	 */
	public static final void filterAlignments_tight(String alignmentsFile, String outputFile, boolean mode, boolean suffixPrefixMode, int minIntersection_nonrepetitive, int minIntersection_repetitive, int minBlueIntervalLength, long[][] out) throws IOException {
		final int DISTANCE_THRESHOLD = IO.quantum;
		boolean orientation, overlapsUniqueA, straddlesLeftA, straddlesRightA;
		int p, q;
		int row, readA, readB, startA, endA, startB, endB, lengthA, lengthB, type;
		int lastFullyContained, lastFullyUnique, lastTranslated, lastBlueInterval, readAInTranslated;
		final int nFullyContained = fullyContained.length;
		final int nFullyUnique = fullyUnique.length;
		final int nTranslated = translated.length;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(outputFile));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); row=0; 
		lastFullyContained=0; lastFullyUnique=0; lastTranslated=0; lastBlueInterval=0;
		while (str!=null)  {
			if (row%1000000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			type=Alignments.readAlignmentFile_getType(DISTANCE_THRESHOLD);
			out[0][type]++;
			// Processing readA
			readA=Alignments.readA-1; readB=Alignments.readB-1; orientation=Alignments.orientation;			
			startA=Alignments.startA; endA=Alignments.endA;
			startB=Alignments.startB; endB=Alignments.endB;
			lengthA=Reads.getReadLength(readA); lengthB=Reads.getReadLength(readB);
			readAInTranslated=-1;
			while (lastFullyUnique<nFullyUnique && fullyUnique[lastFullyUnique]<readA) lastFullyUnique++;
			if (lastFullyUnique<nFullyUnique && fullyUnique[lastFullyUnique]==readA) {
				overlapsUniqueA=true; straddlesLeftA=false; straddlesRightA=false;
			}
			else {
				while (lastFullyContained<nFullyContained && fullyContained[lastFullyContained]<readA) lastFullyContained++;
				if (lastFullyContained<nFullyContained && fullyContained[lastFullyContained]==readA) {
					bw.write("0\n"); str=br.readLine(); row++;					
					continue;
				}
				while (lastTranslated<nTranslated && translated[lastTranslated]<readA) lastTranslated++;
				if (translated[lastTranslated]!=readA) {
					System.err.println("filterAlignments_tight> ERROR: readA not found in translated: "+readA+" :: "+translated[lastTranslated]);
					System.exit(1);
				}
				readAInTranslated=lastTranslated;
				while (lastBlueInterval<=blueIntervals_last && blueIntervals_reads[lastBlueInterval]<readA) lastBlueInterval++;
				q=inBlueRegion(readA,startA,endA,lastTranslated,(lastBlueInterval<=blueIntervals_last&&blueIntervals_reads[lastBlueInterval]==readA)?lastBlueInterval:-1,-1,lengthA,minIntersection_nonrepetitive,minIntersection_repetitive,minBlueIntervalLength);
				if (q==-1) {
					bw.write("0\n"); str=br.readLine(); row++;
					continue;
				}
				overlapsUniqueA=q==0||q==1; straddlesLeftA=q==3||q==5; straddlesRightA=q==4||q==5;
				if ( !overlapsUniqueA && q!=2 && 
					 ( !suffixPrefixMode || 
					   ( !(startA<=DISTANCE_THRESHOLD && straddlesRightA) && 
						 !(endA>=lengthA-DISTANCE_THRESHOLD && straddlesLeftA)
					   )
					 )
				   ) {
   					bw.write("0\n"); str=br.readLine(); row++;
   					continue;
				}
			}
			// Processing readB
			if (readInArray(readB,fullyUnique,nFullyUnique-1,lastFullyUnique)>=0) {
				if (mode) {
					if (overlapsUniqueA) {
						bw.write("1\n");
						out[1][type]++; out[2][type]++;
					}
					else bw.write("0\n");
				}
				else {
					bw.write("1\n");
					out[1][type]++;
					if (overlapsUniqueA) out[2][type]++;
				}
				str=br.readLine(); row++;
				continue;
			}
			else if (readInArray(readB,fullyContained,nFullyContained-1,lastFullyContained)>=0) {
				bw.write("0\n"); str=br.readLine(); row++;
				continue;
			}
			p=readInArray(readB,translated,nTranslated-1,lastTranslated);
			q=inBlueRegion(readB,startB,endB,p,-2,lastBlueInterval,Reads.getReadLength(readB),minIntersection_nonrepetitive,minIntersection_repetitive,minBlueIntervalLength);
			if (q==-1) bw.write("0\n");
			else if (q==0 || q==1) {
				if (mode) {
					if (overlapsUniqueA) {
						bw.write("1\n");
						out[1][type]++; out[2][type]++;
					}
					else bw.write("0\n");
				}
				else {
					bw.write("1\n");
					out[1][type]++;
					if (overlapsUniqueA) out[2][type]++;
				}
			}
			else if (q==2) {
				if (overlapsUniqueA) {
					if (mode) bw.write("0\n");
					else {
						bw.write("1\n");
						out[1][type]++;
					}
				}
				else if (straddlesLeftA) {
					if ( suffixPrefixMode && 
						 ( (orientation && startB<=DISTANCE_THRESHOLD) || 
						   (!orientation && endB>=lengthB-DISTANCE_THRESHOLD)
						 )
					   ) {
						if (mode) {
							if (sameFactorization(readAInTranslated,startA,endA,p,startB,endB,orientation)) {
								bw.write("1\n");
								out[1][type]++;
							}
							else bw.write("0\n");
						}
						else {
							bw.write("1\n");
							out[1][type]++;
						}
					}
					else bw.write("0\n");
				}
				else if (straddlesRightA) {
					if ( suffixPrefixMode && 
					     ( (orientation && endB>=lengthB-DISTANCE_THRESHOLD) || 
					       (!orientation && startB<=DISTANCE_THRESHOLD)
						 )
					   ) {
						if (mode) {
							if (sameFactorization(readAInTranslated,startA,endA,p,startB,endB,orientation)) {
								bw.write("1\n");
								out[1][type]++;
							}
							else bw.write("0\n");
						}
						else {
							bw.write("1\n");
							out[1][type]++;
						}
					}
					else bw.write("0\n");
				}
				else {
					if (mode) {
						if (sameFactorization(readAInTranslated,startA,endA,p,startB,endB,orientation)) {
							bw.write("1\n");
							out[1][type]++;
						}
						else bw.write("0\n");
					}
					else {
						bw.write("1\n");
						out[1][type]++;
					}
				}
			}
			else if (q==3 || q==5) {
				if (overlapsUniqueA || straddlesLeftA || straddlesRightA) bw.write("0\n");
				else {
					if ( suffixPrefixMode && 
					     ( (orientation && startA<=DISTANCE_THRESHOLD) || 
					       (!orientation && endA>=lengthA-DISTANCE_THRESHOLD)
						 )
					   ) {
						if (mode) {
							if (sameFactorization(readAInTranslated,startA,endA,p,startB,endB,orientation)) {
								bw.write("1\n");
								out[1][type]++;
							}
							else bw.write("0\n");
						}
						else {
							bw.write("1\n");
							out[1][type]++;
						}
					}
					else bw.write("0\n");
				}
			}
			else if (q==4 || q==5) {
				if (overlapsUniqueA || straddlesLeftA || straddlesRightA) bw.write("0\n");
				else {
					if ( suffixPrefixMode && 
					     ( (orientation && endA>=lengthA-DISTANCE_THRESHOLD) || 
					       (!orientation && startA<=DISTANCE_THRESHOLD)
					     )
					   ) {
						if (mode) {
							if (sameFactorization(readAInTranslated,startA,endA,p,startB,endB,orientation)) {
								bw.write("1\n");
								out[1][type]++;
							}
							else bw.write("0\n");
						}
						else {
							bw.write("1\n");
							out[1][type]++;
						}
					}
					else bw.write("0\n");
				}
			}
			str=br.readLine(); row++;
		}
		br.close(); bw.close();
	}
	
	
	/**
	 * The dual of $inRedRegion()$. The procedure needs the following arrays: 
	 * boundaries_all, translation_all, blueIntervals, blueIntervals_reads.
	 *
	 * @param minBlueIntervalLength in blocks, not characters;
	 * @return interval $readID[intervalStart..intervalEnd]$:
	 * 0: is fully contained in a non-repetitive region, or fully contains a non-
	 *    repetitive region, and the overlap with the non-repetitive region has length 
	 *    $>=minIntersection_nonrepetitive$;
	 * 1: straddles a non-repetitive region by $>=minIntersection_nonrepetitive$ bps;
	 * 2: contains a sequence of $>=minBlueIntervalLength$ repeat characters that is 
	 *    likely to occur <=H times in the genome, where H is the number of haplotypes;
	 * 3: straddles, but does not fully contain, a sequence in point (2), on the left side
	 *    of the interval;
	 * 4: straddles, but does not fully contain, a sequence in point (2), on the right
	 *    side of the interval;
	 * 5: both (3) and (4) are true;
	 * -1: none of the above is true.
	 */
	private static final int inBlueRegion(int readID, int intervalStart, int intervalEnd, int boundariesAllID, int blueIntervalsID, int blueIntervalsStart, int readLength, int minIntersection_nonrepetitive, int minIntersection_repetitive, int minBlueIntervalLength) {
		boolean straddlesLeft, straddlesRight;
		int i, j;
		int start, end, blockStart, blockEnd, firstBlock, lastBlock;
		final int nBlocks = boundaries_all[boundariesAllID].length+1;
		final int nBytes = Math.ceil(nBlocks,8);
		
		// Checking the nonrepetitive blocks of the read, if any.
		blockStart=0; blockEnd=boundaries_all[boundariesAllID][0];
		if (translation_all[boundariesAllID][0][0]<=lastUnique || translation_all[boundariesAllID][0][0]==lastAlphabet+1) {
			if (Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive) {
				if ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
				     Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) || 
					 Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd)
				   ) return 0;
				else return 1;
			}
		}
		for (i=1; i<nBlocks-1; i++) {
			blockStart=boundaries_all[boundariesAllID][i-1];
			blockEnd=boundaries_all[boundariesAllID][i];
			if (translation_all[boundariesAllID][i][0]<=lastUnique || translation_all[boundariesAllID][i][0]==lastAlphabet+1) {
				if (Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive) {
					if ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
					     Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
						 Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd)
					   ) return 0;
					else return 1;
				}
			}
		}
		blockStart=boundaries_all[boundariesAllID][nBlocks-2];
		blockEnd=readLength-1;		
		if (translation_all[boundariesAllID][nBlocks-1][0]<=lastUnique || translation_all[boundariesAllID][nBlocks-1][0]==lastAlphabet+1) {
			if (Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive) {
				if ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
				     Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
					 Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd)
				   ) return 0;
				else return 1;
			}
		}
		
		// Checking the repetitive blocks of the read, if any.
		if (blueIntervalsID==-2) blueIntervalsID=readInArray(readID,blueIntervals_reads,blueIntervals_reads.length-1,blueIntervalsStart);
		if (blueIntervalsID!=-1) {
			straddlesLeft=false; straddlesRight=false;
			for (i=0; i<blueIntervals[blueIntervalsID].length; i+=3) {
				if (blueIntervals[blueIntervalsID][i+1]<minBlueIntervalLength) continue;
				firstBlock=blueIntervals[blueIntervalsID][i];
				start=firstBlock==0?0:boundaries_all[boundariesAllID][firstBlock-1];
				lastBlock=firstBlock+blueIntervals[blueIntervalsID][i+1]-1;
				end=lastBlock==nBlocks-1?readLength-1:boundaries_all[boundariesAllID][lastBlock];
				if ( Intervals.areApproximatelyIdentical(start,end,intervalStart,intervalEnd) ||
					 Intervals.isApproximatelyContained(start,end,intervalStart,intervalEnd)
				   ) return 2;
				else if (intervalEnd>=start+minIntersection_repetitive && intervalEnd<end && intervalStart<start) straddlesRight=true;
				else if (end>=intervalStart+minIntersection_repetitive && end<intervalEnd && start<intervalStart) straddlesLeft=true;
			}
			if (straddlesLeft) {
				if (straddlesRight) return 5;
				else return 3;
			}
			else if (straddlesRight) return 4;
		}
		
		return -1;
	}
	
	
	/**
	 * @return TRUE iff the specified interval contains at least 2 non-unique blocks.
	 */
	private static final boolean containsTwoBlocks(int boundariesAllID, int intervalStart, int intervalEnd, int readLength) {
		boolean isUnique;
		int i, c;
		int nContained, blockStart, blockEnd;
		final int nBlocks = boundaries_all[boundariesAllID].length+1;
		
		c=translation_all[boundariesAllID][0][0];
		isUnique=c<=lastUnique||c==lastAlphabet+1;
		blockStart=0;
		blockEnd=boundaries_all[boundariesAllID][0];
		if ( !isUnique &&
			 Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) && 
			 !Intervals.areApproximatelyIdentical(blockStart,blockEnd,intervalStart,intervalEnd) 
		   ) nContained=1;
		else nContained=0;
		for (i=1; i<nBlocks-1; i++) {
			c=translation_all[boundariesAllID][i][0];
			isUnique=c<=lastUnique||c==lastAlphabet+1;
			blockStart=boundaries_all[boundariesAllID][i-1];
			blockEnd=boundaries_all[boundariesAllID][i];
			if ( !isUnique &&
				 Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) && 
				 !Intervals.areApproximatelyIdentical(blockStart,blockEnd,intervalStart,intervalEnd) 
			   ) nContained++;
		}
		c=translation_all[boundariesAllID][nBlocks-1][0];
		isUnique=c<=lastUnique||c==lastAlphabet+1;
		blockStart=boundaries_all[boundariesAllID][nBlocks-2];
		blockEnd=readLength-1;
		if ( !isUnique &&
			 Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) && 
			 !Intervals.areApproximatelyIdentical(blockStart,blockEnd,intervalStart,intervalEnd) 
		   ) nContained++;
		return nContained>=2;
	}
	
	
	/**
	 * @param read* index in $translated_all$;
	 * @return TRUE iff $[startA..endA]$ intersects the same number of blocks as
	 * $[startB..endB]$, with boundaries at similar positions, and with at least one 
	 * matching character per block.
	 */
	private static final boolean sameFactorization(int readA, int startA, int endA, int readB, int startB, int endB, boolean orientation) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		int i;
		final int nBlocksA = translation_all[readA].length;
		final int nBlocksB = translation_all[readB].length;
		int nBlocks, blockEnd, firstBlockA, firstBlockB, lastBlockA, lastBlockB;
		double ratio;
		
		// Computing first and last block of readA and readB
		firstBlockA=-1;
		for (i=0; i<=nBlocksA-2; i++) {
			blockEnd=boundaries_all[readA][i];
			if (startA<blockEnd-IDENTITY_THRESHOLD) {
				firstBlockA=i;
				break;
			}
		}
		if (firstBlockA==-1) firstBlockA=nBlocksA-1;		
		lastBlockA=-1;
		for (i=firstBlockA; i<=nBlocksA-2; i++) {
			blockEnd=boundaries_all[readA][i];
			if (endA<=blockEnd+IDENTITY_THRESHOLD) {
				lastBlockA=i;
				break;
			}
		}
		if (lastBlockA==-1) lastBlockA=nBlocksA-1;
		firstBlockB=-1;
		for (i=0; i<=nBlocksB-2; i++) {
			blockEnd=boundaries_all[readB][i];
			if (startA<blockEnd-IDENTITY_THRESHOLD) {
				firstBlockB=i;
				break;
			}
		}
		if (firstBlockB==-1) firstBlockB=nBlocksB-1;		
		lastBlockB=-1;	
		for (i=firstBlockB; i<=nBlocksB-2; i++) {
			blockEnd=boundaries_all[readB][i];
			if (endA<=blockEnd+IDENTITY_THRESHOLD) {
				lastBlockB=i;
				break;
			}
		}
		if (lastBlockB==-1) lastBlockB=nBlocksB-1;
		
		// Finding matching blocks
		nBlocks=lastBlockA-firstBlockA+1;
		if (nBlocks!=lastBlockB-firstBlockB+1) return false;
		ratio=((double)(endB-startB+1))/(endA-startA+1);
		if (orientation) {
			for (i=0; i<nBlocks-1; i++) {
				if ( Math.abs((boundaries_all[readA][firstBlockA+i]-startA)*ratio-(boundaries_all[readB][firstBlockB+i]-startB))>IDENTITY_THRESHOLD ||
					 !nonemptyIntersection(readA,firstBlockA+i,readB,firstBlockB+i,true)
				   ) return false;
			}
			if (!nonemptyIntersection(readA,lastBlockA,readB,lastBlockB,true)) return false;
		}
		else {
			for (i=0; i<nBlocks-1; i++) {
				if ( Math.abs((boundaries_all[readA][firstBlockA+i]-startA)*ratio-(endB-boundaries_all[readB][lastBlockB-i-1]))>IDENTITY_THRESHOLD ||
					 !nonemptyIntersection(readA,firstBlockA+i,readB,lastBlockB-i,false)
				   ) return false;
			}
			if (!nonemptyIntersection(readA,lastBlockA,readB,firstBlockB,false)) return false;
		}
		return true;
	}
	
	
	/**
	 * Remark: the procedure uses global array $stack$ as temporary space.
	 *
	 * @param read* index in $translated_all$;
	 * @return TRUE iff array $translation_all[readA][blockA]$ has at least one character 
	 * in common with the characters (if $mode=TRUE$) or the reverse-complemented 
	 * characters (if $mode=FALSE$) of array $translation_all[readB][blockB]$; characters
	 * are canonized with $canonizeCharacter()$ before being compared.
	 */
	private static final boolean nonemptyIntersection(int readA, int blockA, int readB, int blockB, boolean mode) {
		int i;
		final int lengthA = translation_all[readA][blockA].length;
		final int lengthB = translation_all[readB][blockB].length;
		
		if (stack==null || stack.length<lengthA+lengthB) stack = new int[lengthA+lengthB];
		for (i=0; i<lengthA; i++) stack[i]=canonizeCharacter(translation_all[readA][blockA][i],true);
		for (i=0; i<lengthB; i++) stack[lengthA+i]=canonizeCharacter(translation_all[readB][blockB][i],mode);
		if (lengthA>1) Arrays.sort(stack,0,lengthA);
		if (lengthB>1) Arrays.sort(stack,lengthA,lengthA+lengthB);
		return Math.nonemptyIntersection(stack,0,lengthA-1,stack,lengthA,lengthA+lengthB-1);
	}
	
	
	/**
	 * @return TRUE iff the translated read in $str$ contains a non-repetitive character.
	 */
	public static final boolean containsUnique(String str, int[] boundaries, int readLength, Character tmpChar) {
		final int nBlocks = loadBlocks(str); 
		loadIntBlocks(nBlocks,boundaries,readLength,tmpChar);
		for (int i=0; i<nBlocks; i++) {
			if (isBlockUnique[i]) return true;
		}
		return false;
	}
	
	
	/**
	 * Let a tandem be a maximal sequence of adjacent blocks with characters in common.
	 * The procedure writes a zero for every alignment that should not be trusted because
	 * it either (1) is strictly contained in a tandem, or (2) is identical to a tandem, 
	 * but the tandem is not a blue interval, or (3) straddles or contains a tandem, does
	 * not contain a nonrepetitive character, and all blue intervals in the alignment fall
	 * inside a tandem.
	 *
	 * Remark: every block inside a tandem might be tagged with several characters, and it
	 * might happen that a sequence of such characters occurs just once inside the tandem,
	 * and that it is considered blue based on its total frequency in the read set. In 
	 * practice such a sequence is likely noise.
	 *
	 * @param bothReads discards an alignment if the conditions above hold on both reads
	 * (TRUE) or on just one read (FALSE);
	 * @param minIntersection_nonrepetitive min. length of a non-repetitive substring of 
	 * the alignment, for the alignment not to be considered red; this should not be too 
	 * small, since short non-repetitive regions might not address a unique locus of the 
	 * genome, and they might even be occurrences of short repeats that were not aligned 
	 * to the repeat database because of heuristics of the aligner;
	 * @param out output array containing the number of alignments for each type (columns)
	 * specified in $Alignments.readAlignmentFile_getType()$; row 0: all alignments in 
	 * input; row 1: all alignments kept in output.
	 */
	public static final void filterAlignments_tandem(String alignmentsFile, boolean bothReads, int minIntersection_nonrepetitive, int minAlignmentLength_readRepeat, String outputFile, long[][] out) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int TANDEM_DISTANCE_THRESHOLD = minAlignmentLength_readRepeat;  // Arbitrary
		boolean found, containedA, identicalA, containedB, identicalB;
		int i, p;
		int length, lengthA, lengthB, row, type, nBlocks, readA, readB, startA, endA, startB, endB;
		int lastTranslated, lastBlueInterval, blueIntervalA, blueIntervalB;
		int firstTandemPosition, lastTandemPosition, firstTandemBlockA, lastTandemBlockA, firstTandemBlockB, lastTandemBlockB;
		final int nTranslated = translated.length;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(outputFile));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); row=0; 
		lastTranslated=0; lastBlueInterval=0;
		while (str!=null)  {
			if (row%1000000==0) System.err.println("Processed "+row+" alignments (tandem)");
			Alignments.readAlignmentFile(str);
			type=Alignments.readAlignmentFile_getType(IDENTITY_THRESHOLD);
			out[0][type]++;
			readA=Alignments.readA-1; readB=Alignments.readB-1;
			if ( (bothReads && (lastTandem[readA]==-1 || lastTandem[readB]==-1)) ||
				 (!bothReads && (lastTandem[readA]==-1 && lastTandem[readB]==-1))
			   ) {
				out[1][type]++;
				bw.write("1\n"); str=br.readLine(); row++;
				continue;
			}
			// Processing readA
			lengthA=Reads.getReadLength(readA);
			while (lastTranslated<nTranslated && translated[lastTranslated]<readA) lastTranslated++;
			if (lastTranslated<nTranslated && translated[lastTranslated]==readA) {
				startA=Alignments.startA; endA=Alignments.endA;
				while (lastBlueInterval<=blueIntervals_last && blueIntervals_reads[lastBlueInterval]<readA) lastBlueInterval++;
				if (lastBlueInterval>blueIntervals_last || blueIntervals_reads[lastBlueInterval]!=readA) blueIntervalA=-1;
				else blueIntervalA=lastBlueInterval;
				nBlocks=boundaries_all[lastTranslated].length+1;
				containedA=false; identicalA=false; firstTandemBlockA=-1; lastTandemBlockA=-1;
				for (i=0; i<lastTandem[readA]; i+=2) {
					if (tandems[readA][i]==0) firstTandemPosition=0;
					else firstTandemPosition=boundaries_all[lastTranslated][tandems[readA][i]-1];
					if (tandems[readA][i+1]==nBlocks-1) lastTandemPosition=lengthA-1;
					else lastTandemPosition=boundaries_all[lastTranslated][tandems[readA][i+1]];
					if (Intervals.areApproximatelyIdentical(startA,endA,firstTandemPosition,lastTandemPosition)) identicalA=true;
					else if (Intervals.isApproximatelyContained(startA,endA,firstTandemPosition,lastTandemPosition)) containedA=true;
					if (identicalA || containedA) {
						firstTandemBlockA=tandems[readA][i]; 
						lastTandemBlockA=tandems[readA][i+1];
						break;
					}
				}
				if (!bothReads) {
					if ( containedA || 
					     (identicalA && !filterAlignments_tandem_isBlue(blueIntervalA,firstTandemBlockA,lastTandemBlockA)) || 
					     (!identicalA && filterAlignments_tandem_allContainedBlueInTandem(readA,lastTranslated,blueIntervalA,startA,endA,nBlocks,lengthA,TANDEM_DISTANCE_THRESHOLD) && !filterAlignments_tandem_containsUnique(lastTranslated,startA,endA,nBlocks,lengthA,minIntersection_nonrepetitive))
					   ) {
						bw.write("0\n"); str=br.readLine(); row++;
						continue;
					}
				}
				else {
					if ( (identicalA && filterAlignments_tandem_isBlue(blueIntervalA,firstTandemBlockA,lastTandemBlockA)) || 
					     (!identicalA && !containedA && !filterAlignments_tandem_allContainedBlueInTandem(readA,lastTranslated,blueIntervalA,startA,endA,nBlocks,lengthA,TANDEM_DISTANCE_THRESHOLD))
					   ) {
						out[1][type]++;
						bw.write("1\n"); str=br.readLine(); row++;
						continue;
					}
				}
			}
			else if (bothReads) {
				out[1][type]++;
				bw.write("1\n"); str=br.readLine(); row++;
				continue;
			}
			// Processing readB
			lengthB=Reads.getReadLength(readB);
			p=readInArray(readB,translated,nTranslated-1,lastTranslated);
			if (p>=0) {
				blueIntervalB=readInArray(readB,blueIntervals_reads,blueIntervals_reads.length-1,lastBlueInterval);
				startB=Alignments.startB; endB=Alignments.endB;
				nBlocks=boundaries_all[p].length+1;
				containedB=false; identicalB=false; firstTandemBlockB=-1; lastTandemBlockB=-1;
				for (i=0; i<lastTandem[readB]; i+=2) {
					if (tandems[readB][i]==0) firstTandemPosition=0;
					else firstTandemPosition=boundaries_all[p][tandems[readB][i]-1];
					if (tandems[readB][i+1]==nBlocks-1) lastTandemPosition=lengthB-1;
					else lastTandemPosition=boundaries_all[p][tandems[readB][i+1]];
					if (Intervals.areApproximatelyIdentical(startB,endB,firstTandemPosition,lastTandemPosition)) identicalB=true;
					else if (Intervals.isApproximatelyContained(startB,endB,firstTandemPosition,lastTandemPosition)) containedB=true;
					if (identicalB || containedB) {
						firstTandemBlockB=tandems[readB][i];
						lastTandemBlockB=tandems[readB][i+1];
						break;
					}
				}
				if (!bothReads) {
					if ( containedB || 
					     (identicalB && !filterAlignments_tandem_isBlue(blueIntervalB,firstTandemBlockB,lastTandemBlockB)) || 
					     (!identicalB && filterAlignments_tandem_allContainedBlueInTandem(readB,p,blueIntervalB,startB,endB,nBlocks,lengthB,TANDEM_DISTANCE_THRESHOLD) && !filterAlignments_tandem_containsUnique(p,startB,endB,nBlocks,lengthB,minIntersection_nonrepetitive))
					   ) {
						bw.write("0\n"); str=br.readLine(); row++;
						continue;
					}
				}
				else {
					if ( (identicalB && filterAlignments_tandem_isBlue(blueIntervalB,firstTandemBlockB,lastTandemBlockB)) || 
					     (!identicalB && !containedB && !filterAlignments_tandem_allContainedBlueInTandem(readB,p,blueIntervalB,startB,endB,nBlocks,lengthB,TANDEM_DISTANCE_THRESHOLD))
					   ) {
						out[1][type]++;
						bw.write("1\n"); str=br.readLine(); row++;
						continue;
					}
				}
			}
			else if (bothReads) {
				out[1][type]++;
				bw.write("1\n"); str=br.readLine(); row++;
				continue;
			}
			if (bothReads) bw.write("0\n");
			else bw.write("1\n");
			str=br.readLine(); row++;
		}
		br.close(); bw.close();
	}
	
	
	/**
	 * @return TRUE iff $[firstTandemBlock..lastTandemBlock]$ is in 
	 * $blueIntervals[blueIntervalsID]$.
	 */
	private static final boolean filterAlignments_tandem_isBlue(int blueIntervalsID, int firstTandemBlock, int lastTandemBlock) {
		int i;
		int firstBlock, lastBlock, lastBlueInterval;
		
		if (blueIntervalsID==-1) return false;
		lastBlueInterval=blueIntervals[blueIntervalsID].length-1;
		for (i=0; i<lastBlueInterval; i+=3) {
			firstBlock=blueIntervals[blueIntervalsID][i];
			if (firstBlock>firstTandemBlock) break;
			else if (firstBlock==firstTandemBlock) {
				lastBlock=firstBlock+blueIntervals[blueIntervalsID][i+1]-1;
				if (lastBlock==lastTandemBlock) return true;
			}
		}
		return false;
	}

	
	/**
	 * @return TRUE iff $[intervalStart..intervalEnd]$ strictly contains blue intervals,
	 * and all blue intervals strictly contained in the interval, are also strictly 
	 * contained in a tandem. So, if $[intervalStart..intervalEnd]$ contains a tandem, and
	 * the tandem coincides with a blue interval, the procedure returns FALSE. Blue 
	 * intervals that coincide with $[intervalStart..intervalEnd]$ are not considered.
	 * @param distanceThreshold a contained blue interval that coincides with a tandem 
	 * interval, but that is at most this far from the read start/end, is considered
	 * strictly contained in the tandem interval.
	 */
	private static final boolean filterAlignments_tandem_allContainedBlueInTandem(int readID, int boundariesAllID, int blueIntervalsID, int intervalStart, int intervalEnd, int nBlocks, int readLength, int distanceThreshold) {
		boolean containsBlue, found;
		int i, j;
		int firstBlock, lastBlock, firstPosition, lastPosition, lastBlueInterval, lastTandemInterval;
		
		if (blueIntervalsID==-1) return false;
		lastBlueInterval=blueIntervals[blueIntervalsID].length-1;
		lastTandemInterval=lastTandem[readID];
		if (lastTandemInterval==-1) return false;
		j=0; containsBlue=false;
		for (i=0; i<lastBlueInterval; i+=3) {
			firstBlock=blueIntervals[blueIntervalsID][i];
			if (firstBlock==0) firstPosition=0;
			else firstPosition=boundaries_all[boundariesAllID][firstBlock-1];
			if (firstPosition>=intervalEnd) break;
			lastBlock=firstBlock+blueIntervals[blueIntervalsID][i+1]-1;
			if (lastBlock==nBlocks-1) lastPosition=Reads.getReadLength(readID)-1;
			else lastPosition=boundaries_all[boundariesAllID][lastBlock];
			if (lastPosition<=intervalStart) continue;
			if ( Intervals.isApproximatelyContained(firstPosition,lastPosition,intervalStart,intervalEnd) &&
				 !Intervals.areApproximatelyIdentical(firstPosition,lastPosition,intervalStart,intervalEnd)
			   ) {
				containsBlue=true; found=false;
				while (j<lastTandemInterval) {
					if (tandems[readID][j+1]<lastBlock) {
						j+=2;
						continue;
					}
					if (tandems[readID][j]>firstBlock) break;
					if ( tandems[readID][j]<firstBlock || tandems[readID][j+1]>lastBlock ||
						 (tandems[readID][j]==firstBlock && tandems[readID][j+1]==lastBlock && (firstPosition<distanceThreshold || lastPosition>readLength-distanceThreshold))
					   ) {
						found=true;
						break;
					}
					j+=2;
				}
				if (!found) return false;
			}
		}
		return containsBlue;
	}
	
	
	/**
	 * Remark: this procedure assumes that $translation_all$ has already been loaded.
	 *
	 * @return TRUE iff $[intervalStart..intervalEnd]$ contains a non-repetitive block.
	 */
	private static final boolean filterAlignments_tandem_containsUnique(int boundariesAllID, int intervalStart, int intervalEnd, int nBlocks, int readLength, int minIntersection_nonrepetitive) {
		int i, j, c;
		int start, end, last;
		
		for (i=0; i<=nBlocks-1; i++) {
			last=translation_all[boundariesAllID][i].length-1;
			for (j=0; j<=last; j++) {
				c=translation_all[boundariesAllID][i][j];
				if (c>lastUnique && c<=lastAlphabet) break;
				start=i==0?0:boundaries_all[boundariesAllID][i-1]+1;
				end=i==nBlocks-1?readLength-1:boundaries_all[boundariesAllID][i];
				if (end-start+1>=minIntersection_nonrepetitive && Intervals.isApproximatelyContained(start,end,intervalStart,intervalEnd)) return true;
			}
		}
		return false;
	}
	
	
	
	
	
	
	
	
	// ------------------------- READ BREAKING AT LOW QUALITY ----------------------------
	
	/**
	 * Assume that the input reads are CLR and that they contain long low-quality regions.
	 * The procedure translates read-read and read-repeat alignments based on the output 
	 * of $Reads.breakReads()$. This transformation allows handling low-quality regions 
	 * correctly while translating the reads in the alphabet of repeats: specifically, the
	 * whole translation pipeline can be run as a black box on the alignments of broken
	 * reads.
	 *
	 * Remark: assume that we just run the repeat translation pipeline as-is on reads with 
	 * long low-quality regions. If a low-quality region $L$ falls inside a repetitive 
	 * region $RLS$, the repetitive substrings $R$ and $S$ would be modeled as closed 
	 * characters and, since low-quality regions are randomly distributed, such characters
	 * would likely have low frequency in the dataset and they would likely be filtered 
	 * out of the alphabet. $R$ and $S$ would thus be converted into non-repetitive 
	 * substrings, and some repeat-induced alignments that contain them might be preserved
	 * because of this.
	 * 
	 * Remark: one might think of modeling $R,S$ above as open characters. Open characters 
	 * might be compressed away during alphabet construction, so $R,S$ might get replaced 
	 * with closed characters of much bigger length. The translated read might then become
	 * $UVW$, where $U,V$ are closed characters (possibly the same closed character) with 
	 * $|U|>|R|$, $|W|>|S|$, and where $V$ is non-repetitive. Any k-mer built on such a 
	 * translation would not reflect the real structure of the read. One might assume that
	 * such a k-mer would likely be globally infrequent and would thus be removed, but it 
	 * might be kept if it happens to be frequent.
	 *
	 * Remark: breaking the read at $L$ above would still allow to model $R,S$ as half-
	 * open, while at the same time taking advantage of the existing determinization 
	 * procedures if $R,S$ happen to match several repeats, and forbidding any k-mer to 
	 * contain $L$. 
	 *
	 * Remark: breaking reads is acceptable in this context, since we don't need their
	 * long-range information (as we do in assembly): we just need to translate reads in 
	 * the alphabet of repeats and to mark unique substrings in the recoded alphabet, and
	 * such unique substrings cannot contain a low-quality region anyway.
	 *
	 * Remark: "short" low-quality intervals pose the same problems as "long" low-quality 
	 * intervals, since the quality track comes from observing low coverage, and possibly
	 * alignments breaks, in both cases. It makes sense to break reads even at short low-
	 * quality intervals.
	 *
	 * Remark: the procedure discards any alignment that contains a low-quality region in
	 * readA or readB, and it trims the readA and the readB side of an alignment whose 
	 * suffix/prefix overlaps a low-quality region (the readA side and the readB side are
	 * trimmed independently).
	 *
	 * Remark: the translated alignments can have any length. It is left to procedures
	 * downstream to keep only long-enough alignments.
	 *
	 * Remark: the procedure assumes that $Reads.breakReads_old2new$ has already been 
	 * loaded, and that $Reads.readLengths$ has been already initialized with the lengths
	 * of the old reads.
	 *
	 * @param inputFile assumed to be sorted by $readA,readB,orientation,startA,startB$;
	 * @param translateB if FALSE, the readB side of every alignment is kept intact (this
	 * is useful if $inputFile$ contains read-repeat alignments);
	 * @param addHeader adds the two header lines to the output;
	 * @param outputFile sorted in the same order as $inputFile$.
	 */
	public static final void breakReads_translateAlignments(String inputFile, boolean translateB, boolean addHeader, String outputFile) throws IOException {
		final int ALIGNMENTS_CAPACITY = 100;  // Arbitrary
		boolean found;
		int i;
		int last, end, intersection, maxIntersection, readA, readB, currentReadA;
		int readA_new, startA_new, endA_new, readB_new, startB_new, endB_new, diffs_new;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		AlignmentRow.order=AlignmentRow.ORDER_READA_READB_ORIENTATION_STARTA_STARTB_ENDA_ENDB;
		if (alignments==null || alignments.length<ALIGNMENTS_CAPACITY) alignments = new AlignmentRow[ALIGNMENTS_CAPACITY];
		for (i=0; i<alignments.length; i++) {	
			if (alignments[i]==null) alignments[i] = new AlignmentRow();
		}
		bw = new BufferedWriter(new FileWriter(outputFile));
		if (addHeader) { bw.newLine(); bw.newLine(); }
		br = new BufferedReader(new FileReader(inputFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); currentReadA=-1;
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			readA=Alignments.readA-1;
			if (readA!=currentReadA) {
				if (currentReadA!=-1) {
					if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
					for (i=0; i<=lastAlignment; i++) bw.write((alignments[i].readA+1)+"  "+(alignments[i].readB+1)+"  "+(alignments[i].orientation?'n':'c')+"  ["+alignments[i].startA+".. "+alignments[i].endA+"] x ["+(alignments[i].orientation?alignments[i].startB+".. "+alignments[i].endB:alignments[i].endB+".. "+alignments[i].startB)+"] ( "+alignments[i].diffs+" diffs)\n");
				}
				currentReadA=readA; lastAlignment=-1;
			}
			readB=Alignments.readB-1;
			if (Reads.breakReads_containsLowQuality(readA,Alignments.startA,Alignments.endA) || (translateB && Reads.breakReads_containsLowQuality(readB,Alignments.startB,Alignments.endB))) {
				str=br.readLine();
				continue;
			}
			maxIntersection=0; startA_new=-1; endA_new=-1; readA_new=-1;
			for (i=0; i<Reads.last_old2new[readA]; i+=3) {
				if (Reads.breakReads_old2new[readA][i]>Alignments.endA) break;
				if (Reads.breakReads_old2new[readA][i+1]<Alignments.startA) continue;
				intersection=Intervals.intersectionLength(Alignments.startA,Alignments.endA,Reads.breakReads_old2new[readA][i],Reads.breakReads_old2new[readA][i+1]);
				if (intersection>maxIntersection) {
					maxIntersection=intersection;
					startA_new=Reads.breakReads_old2new[readA][i];
					endA_new=Reads.breakReads_old2new[readA][i+1];
					readA_new=Reads.breakReads_old2new[readA][i+2];
				}
			}
			endA_new=Math.min(endA_new,Alignments.endA)-startA_new;
			startA_new=Math.max(startA_new,Alignments.startA)-startA_new;
			if (translateB) {
				maxIntersection=0; startB_new=-1; endB_new=-1; readB_new=-1;
				for (i=0; i<Reads.last_old2new[readB]; i+=3) {
					if (Reads.breakReads_old2new[readB][i]>Alignments.endB) break;
					if (Reads.breakReads_old2new[readB][i+1]<Alignments.startB) continue;
					intersection=Intervals.intersectionLength(Alignments.startB,Alignments.endB,Reads.breakReads_old2new[readB][i],Reads.breakReads_old2new[readB][i+1]);
					if (intersection>maxIntersection) {
						maxIntersection=intersection;
						startB_new=Reads.breakReads_old2new[readB][i];
						endB_new=Reads.breakReads_old2new[readB][i+1];
						readB_new=Reads.breakReads_old2new[readB][i+2];
					}
				}
				endB_new=Math.min(endB_new,Alignments.endB)-startB_new;
				startB_new=Math.max(startB_new,Alignments.startB)-startB_new;
			}
			else { startB_new=Alignments.startB; endB_new=Alignments.endB; readB_new=readB; }
			diffs_new=(int)((endA_new-startA_new+1+endB_new-startB_new+1)*((double)Alignments.diffs)/(Alignments.endA-Alignments.startA+1+Alignments.endB-Alignments.startB+1));
			lastAlignment++;
			if (lastAlignment==alignments.length) {
				AlignmentRow[] newAlignments = new AlignmentRow[alignments.length<<1];
				System.arraycopy(alignments,0,newAlignments,0,alignments.length);
				for (i=alignments.length; i<newAlignments.length; i++) newAlignments[i] = new AlignmentRow();
				alignments=newAlignments;
			}
			alignments[lastAlignment].set(readA_new,startA_new,endA_new,readB_new,startB_new,endB_new,Alignments.orientation,diffs_new);
			str=br.readLine();
		}
		br.close();
		if (currentReadA!=-1) {
			if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
			for (i=0; i<=lastAlignment; i++) bw.write((alignments[i].readA+1)+"  "+(alignments[i].readB+1)+"  "+(alignments[i].orientation?'n':'c')+"  ["+alignments[i].startA+".. "+alignments[i].endA+"] x ["+(alignments[i].orientation?alignments[i].startB+".. "+alignments[i].endB:alignments[i].endB+".. "+alignments[i].startB)+"] ( "+alignments[i].diffs+" diffs)\n");
		}
		bw.close();
	}
	
	
	/**
	 * Splits $alignmentsFile_new,bitvectorFile_new$, which are the full alignments file
	 * and the full bitvector of the broken reads, into chunks that correspond to those of
	 * the alignments file of the unbroken reads.
	 *
	 * Remark: the procedure assumes that $Reads.breakReads_old2new$ is already loaded.
	 *
	 * @param lastReadA_old one-based, like in LAshow;
	 * @param bitvectorFile_new ignored if NULL.
	 */
	public static final void breakReads_splitAlignments(int[] lastReadA_old, int nChunks_old, String alignmentsFile_new, String bitvectorFile_new, String alignmentsPrefix_new, String bitvectorPrefix_new) throws IOException {
		int i, j;
		int readA_new, lastReadA_new;
		String str1, str2;
		BufferedReader br1, br2;
		BufferedWriter bw1, bw2;
		
		br1 = new BufferedReader(new FileReader(alignmentsFile_new));
		if (bitvectorFile_new!=null) br2 = new BufferedReader(new FileReader(bitvectorFile_new));
		else br2=null;
		str1=br1.readLine(); str1=br1.readLine();  // Skipping header
		str1=br1.readLine(); j=0;
		if (bitvectorFile_new!=null) str2=br2.readLine();
		else str2=null;
		for (i=0; i<nChunks_old; i++) {
			lastReadA_new=Reads.breakReads_old2new[lastReadA_old[i]-1][Reads.last_old2new[lastReadA_old[i]-1]];
			bw1 = new BufferedWriter(new FileWriter(alignmentsPrefix_new+i+".txt"));
			bw1.newLine(); bw1.newLine();
			if (bitvectorFile_new!=null) bw2 = new BufferedWriter(new FileWriter(bitvectorPrefix_new+i+".txt"));
			else bw2=null;
			readA_new=Alignments.readAlignmentFile_readA(str1)-1;
			while (readA_new<=lastReadA_new) {
				bw1.write(str1); bw1.newLine(); 
				j++;
				if (j%10000000==0) System.err.println("Processed "+j+" alignments");
				if (bitvectorFile_new!=null) { bw2.write(str2); bw2.newLine(); str2=br2.readLine(); }
				else str2=null;
				str1=br1.readLine();
				if (str1==null) break;
				Alignments.readAlignmentFile_readA(str1);
				readA_new=Alignments.readA-1;
			}
			bw1.close();
			if (bitvectorFile_new!=null) bw2.close();
		}
		br1.close(); 
		if (bitvectorFile_new!=null) br2.close();
	}
	
	
	/**
	 * Given bitvectors, that filter the new alignments created by $breakIntervals_
	 * translateAlignments()$, the procedure builds corresponding bitvectors that filter
	 * the old alignments.
	 * 
	 * Remark: the procedure assumes that $Reads.breakReads_new2old$ has already been
	 * loaded, and it uses global array $alignments$.
	 *
	 * @param *bitvectorFile_old output files.
	 */
	public static final void breakReads_translateBitvector(String bitvectorFile_new, String tandemBitvectorFile_new, String alignmentsFile_new, String bitvectorFile_old, String tandemBitvectorFile_old, String alignmentsFile_old) throws IOException {
		final int ALIGNMENTS_CAPACITY = 100;  // Arbitrary
		int i;
		int currentReadA, oldReadA, oldFirstA, oldReadB, oldFirstB;
		String str1, str2, str3, str4;
		BufferedReader br1, br2, br3, br4;
		BufferedWriter bw1, bw2;
		AlignmentRow tmpAlignment = new AlignmentRow();
		
		if (alignments==null || alignments.length<ALIGNMENTS_CAPACITY) alignments = new AlignmentRow[ALIGNMENTS_CAPACITY];
		for (i=0; i<alignments.length; i++) alignments[i] = new AlignmentRow();
		AlignmentRow.order=AlignmentRow.ORDER_READA_READB_ORIENTATION_STARTA_STARTB_ENDA_ENDB;
		bw1 = new BufferedWriter(new FileWriter(bitvectorFile_old));
		bw2 = new BufferedWriter(new FileWriter(tandemBitvectorFile_old));
		br1 = new BufferedReader(new FileReader(alignmentsFile_new));
		str1=br1.readLine(); str1=br1.readLine();  // Skipping header
		str1=br1.readLine();
		br2 = new BufferedReader(new FileReader(bitvectorFile_new));
		str2=br2.readLine();
		br3 = new BufferedReader(new FileReader(alignmentsFile_old));
		str3=br3.readLine(); str3=br3.readLine();  // Skipping header
		str3=br3.readLine();
		br4 = new BufferedReader(new FileReader(tandemBitvectorFile_new));
		str4=br4.readLine();
		currentReadA=-1; lastAlignment=-1;
		while (str1!=null) {
			Alignments.readAlignmentFile_readA(str1);
			oldReadA=Reads.breakReads_new2old[Alignments.readA-1][0]; 
			if (oldReadA!=currentReadA) {
				if (currentReadA!=-1) {
					if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
					// This procedure calls $Alignments.readAlignmentFile()$, so we should
					// call it again afterwards.
					str3=breakReads_translateBitvector_impl(currentReadA,br3,str3,bw1,bw2,tmpAlignment);
				}
				Alignments.readAlignmentFile(str1);
				oldFirstA=Reads.breakReads_new2old[Alignments.readA-1][1];
				oldReadB=Reads.breakReads_new2old[Alignments.readB-1][0]; 
				oldFirstB=Reads.breakReads_new2old[Alignments.readB-1][1];
				currentReadA=oldReadA; lastAlignment=0; 
				alignments[lastAlignment].set(oldReadA,oldFirstA+Alignments.startA,oldFirstA+Alignments.endA,oldReadB,oldFirstB+Alignments.startB,oldFirstB+Alignments.endB,Alignments.orientation,Alignments.diffs);
				alignments[lastAlignment].flag=Integer.parseInt(str2.trim())!=0;
				alignments[lastAlignment].flag2=Integer.parseInt(str4.trim())!=0;
			}
			else {
				lastAlignment++;
				if (lastAlignment==alignments.length) {
					AlignmentRow[] newAlignments = new AlignmentRow[alignments.length<<1];
					System.arraycopy(alignments,0,newAlignments,0,alignments.length);
					for (i=alignments.length; i<newAlignments.length; i++) newAlignments[i] = new AlignmentRow();
					alignments=newAlignments;
				}
				Alignments.readAlignmentFile(str1);
				oldFirstA=Reads.breakReads_new2old[Alignments.readA-1][1];
				oldReadB=Reads.breakReads_new2old[Alignments.readB-1][0]; 
				oldFirstB=Reads.breakReads_new2old[Alignments.readB-1][1];
				alignments[lastAlignment].set(oldReadA,oldFirstA+Alignments.startA,oldFirstA+Alignments.endA,oldReadB,oldFirstB+Alignments.startB,oldFirstB+Alignments.endB,Alignments.orientation,Alignments.diffs);
				alignments[lastAlignment].flag=Integer.parseInt(str2.trim())!=0;
				alignments[lastAlignment].flag2=Integer.parseInt(str4.trim())!=0;
			}
			str1=br1.readLine(); str2=br2.readLine(); str4=br4.readLine();
		}
		if (currentReadA!=-1) {
			if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
			breakReads_translateBitvector_impl(currentReadA,br3,str3,bw1,bw2,tmpAlignment);
		}
		br1.close(); br2.close(); br3.close(); br4.close(); bw1.close(); bw2.close();
	}
	
	
	/**
	 * Remark: an old alignment might not be identical to the translation of its 
	 * corresponding new alignment, because of trimming. Thus, the procedure has to do a
	 * binary search inside $alignments$ and to return an element that maximizes the 
	 * intersection with the old alignment.
	 *
	 * @param br old alignments file;
	 * @param firstString the first string from $br$ to be processed;
	 * @param bw* output files (old bitvectors): 1=unique substrings; 2=tandem.
	 * @param tmpAlignment temporary space;
	 * @return the new value of $firstString$ after the procedure completes.
	 */
	private static final String breakReads_translateBitvector_impl(int currentReadA, BufferedReader br, String firstString, BufferedWriter bw1, BufferedWriter bw2, AlignmentRow tmpAlignment) throws IOException {
		int i, p;
		int readA, readB, intersection, maxIntersection, maxAlignment;
		String str;
		
		str=firstString;
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			readA=Alignments.readA-1;
			if (readA>currentReadA) return str;
			if (readA<currentReadA) {
				bw1.write("0\n"); bw2.write("1\n");
				str=br.readLine();
				continue;
			}
			readB=Alignments.readB-1;
			if (Reads.breakReads_containsLowQuality(readA,Alignments.startA,Alignments.endA) || Reads.breakReads_containsLowQuality(readB,Alignments.startB,Alignments.endB)) {
				bw1.write("0\n"); bw2.write("1\n");
				str=br.readLine();
				continue;
			}
			tmpAlignment.set(readA,Alignments.startA,Alignments.endA,readB,Alignments.startB,Alignments.endB,Alignments.orientation,Alignments.diffs);
			p=Arrays.binarySearch(alignments,0,lastAlignment+1,tmpAlignment);
			if (p>=0) {
				bw1.write(alignments[p].flag?"1\n":"0\n");
				bw2.write(alignments[p].flag2?"1\n":"0\n");
				str=br.readLine();
				continue;
			}
			p=-1-p; maxIntersection=0; maxAlignment=-1;
			for (i=p-1; i>=0; i--) {
				if (alignments[i].readB!=tmpAlignment.readB || alignments[i].orientation!=tmpAlignment.orientation) break;
				intersection=Intervals.intersectionLength(alignments[i].startA,alignments[i].endA,Alignments.startA,Alignments.endA)+Intervals.intersectionLength(alignments[i].startB,alignments[i].endB,Alignments.startB,Alignments.endB);
				if (intersection>maxIntersection) {
					maxIntersection=intersection;
					maxAlignment=i;
				}
			}
			for (i=p; i<=lastAlignment; i++) {
				if (alignments[i].readB!=tmpAlignment.readB || alignments[i].orientation!=tmpAlignment.orientation || alignments[i].startA>=tmpAlignment.endA) break;
				intersection=Intervals.intersectionLength(alignments[i].startA,alignments[i].endA,Alignments.startA,Alignments.endA)+Intervals.intersectionLength(alignments[i].startB,alignments[i].endB,Alignments.startB,Alignments.endB);
				if (intersection>maxIntersection) {
					maxIntersection=intersection;
					maxAlignment=i;
				}
			}
			if (maxAlignment==-1) {
				System.err.println("breakIntervals_translateFilter_impl> ERROR: no new alignment intersects the old alignment "+str);
				System.exit(1);
			}
			bw1.write(alignments[maxAlignment].flag?"1\n":"0\n");
			bw2.write(alignments[maxAlignment].flag2?"1\n":"0\n");
			str=br.readLine();
		}
		return null;
	}
	
	
	/**
	 * Enforces that the result of $fixEndBlocks()$, executed on the broken reads, is
	 * consistent on the two sides of every long low-quality region.
	 * 
	 * Remark: if a long low-quality region replaces a substring without preserving the
	 * substring's length, no constraint can be enforced.
	 *
	 * Remark: the procedure assumes that $Reads.readLengths$ contains the length of every
	 * new read, and that $Reads.breakReads_new2old$ has already been loaded.
	 *
	 * @param mode every long low-quality region: TRUE=replaces a substring with a random 
	 * substring of approx. the same length; FALSE=is an insertion;
	 * @param translatedFile the translated reads before disambiguation;
	 * @param translatedFile_disambiguated_new (output file) updated version of 
	 * $translatedFile_disambiguated$.
	 */
	public static final void breakReads_checkDisambiguation(boolean mode, int lengthThreshold, String translatedFile, String translatedFile_disambiguated, String boundariesFile, String translatedFile_disambiguated_new) throws IOException {
		final int CAPACITY = 10;  // Arbitrary
		boolean success, leftUnique, rightUnique;
		int r;
		int oldRead, currentOldRead, nBlocks, leftLength, nAttempted, nRolledBack;
		int last1, last2, leftEnd1_last, leftEnd2_last, rightEnd1_last, rightEnd2_last;
		String str1, str2, str3;
		StringBuilder leftEnd1_str, leftEnd2_str;
		Character tmpChar;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		int[] leftEnd1, leftEnd2, rightEnd1, rightEnd2;
		int[] tmpArray, tmpArray1, tmpArray2;
		String[] tokens;
		
		tmpChar = new Character();
		leftEnd1 = new int[CAPACITY]; leftEnd2 = new int[CAPACITY];
		rightEnd1 = new int[CAPACITY]; rightEnd2 = new int[CAPACITY];
		tmpArray1 = new int[CAPACITY]; tmpArray2 = new int[CAPACITY];
		leftEnd1_str = new StringBuilder(); leftEnd2_str = new StringBuilder();
		bw = new BufferedWriter(new FileWriter(translatedFile_disambiguated_new));
		br1 = new BufferedReader(new FileReader(translatedFile));
		br2 = new BufferedReader(new FileReader(translatedFile_disambiguated));
		br3 = new BufferedReader(new FileReader(boundariesFile));
		nAttempted=0; nRolledBack=0;
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); r=-1;
		currentOldRead=-1; leftEnd1_last=-1; leftEnd2_last=-1; rightEnd1_last=-1; rightEnd2_last=-1; leftLength=0;
		leftUnique=false; rightUnique=false;
		while (str1!=null) {
			r++;
			oldRead=Reads.breakReads_new2old[r][0]; 
			if (oldRead!=currentOldRead) {
				if (leftEnd2_str.length()!=0) bw.write(SEPARATOR_MAJOR+leftEnd2_str.toString()+"\n");
				currentOldRead=oldRead;
				nBlocks=loadBlocks(str1);
				if (nBlocks<2) {
					bw.write(str2); bw.newLine();
					leftEnd1_str.delete(0,leftEnd1_str.length());
					leftEnd2_str.delete(0,leftEnd2_str.length());
					leftEnd1_last=-1; leftEnd2_last=-1; leftLength=0; leftUnique=false/*Irrelevant*/;
				}
				else {
					bw.write(str2.substring(0,str2.lastIndexOf(SEPARATOR_MAJOR+"")));
					loadBoundaries(str3);
					loadIntBlocks(nBlocks,boundaries,Reads.readLengths[r],tmpChar);
					if (lastInBlock_int[nBlocks-1]>=leftEnd1.length) leftEnd1 = new int[lastInBlock_int[nBlocks-1]+1];
					System.arraycopy(intBlocks[nBlocks-1],0,leftEnd1,0,lastInBlock_int[nBlocks-1]+1);
					leftEnd1_last=lastInBlock_int[nBlocks-1];
					leftEnd1_str.delete(0,leftEnd1_str.length());
					leftEnd1_str.append(str1.lastIndexOf(SEPARATOR_MAJOR+"")+1);
					nBlocks=loadBlocks(str2); loadIntBlocks(nBlocks,boundaries,Reads.readLengths[r],tmpChar);
					if (lastInBlock_int[nBlocks-1]>=leftEnd2.length) leftEnd2 = new int[lastInBlock_int[nBlocks-1]+1];
					System.arraycopy(intBlocks[nBlocks-1],0,leftEnd2,0,lastInBlock_int[nBlocks-1]+1);
					leftEnd2_last=lastInBlock_int[nBlocks-1];
					leftEnd2_str.delete(0,leftEnd2_str.length());
					leftEnd2_str.append(str2.lastIndexOf(SEPARATOR_MAJOR+"")+1);
					leftLength=Reads.readLengths[r]-boundaries[nBlocks-2]-1;
					leftUnique=isBlockUnique[nBlocks-1];
				}
				str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
				continue;
			}
			nBlocks=loadBlocks(str1); 
			if (nBlocks<2) {
				if (leftEnd2_str.length()!=0) bw.write(SEPARATOR_MAJOR+leftEnd2_str.toString()+"\n");
				bw.write(str2); bw.newLine();
				leftEnd1_str.delete(0,leftEnd1_str.length());
				leftEnd2_str.delete(0,leftEnd2_str.length());
				leftEnd1_last=-1; leftEnd2_last=-1; leftLength=0; leftUnique=false/*Irrelevant*/;
				str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
				continue;
			}
			loadBoundaries(str3);
			loadIntBlocks(nBlocks,boundaries,Reads.readLengths[r],tmpChar);
			if (lastInBlock_int[0]>=rightEnd1.length) rightEnd1 = new int[lastInBlock_int[0]+1];
			System.arraycopy(intBlocks[0],0,rightEnd1,0,lastInBlock_int[0]+1);
			rightEnd1_last=lastInBlock_int[0];
			if (lastInBlock_int[nBlocks-1]>=tmpArray1.length) tmpArray1 = new int[lastInBlock_int[nBlocks-1]+1];
			System.arraycopy(intBlocks[nBlocks-1],0,tmpArray1,0,lastInBlock_int[nBlocks-1]+1);
			last1=lastInBlock_int[nBlocks-1];
			nBlocks=loadBlocks(str2); loadIntBlocks(nBlocks,boundaries,Reads.readLengths[r],tmpChar);
			if (lastInBlock_int[0]>=rightEnd2.length) rightEnd2 = new int[lastInBlock_int[0]+1];
			System.arraycopy(intBlocks[0],0,rightEnd2,0,lastInBlock_int[0]+1);
			rightEnd2_last=lastInBlock_int[0];
			if (lastInBlock_int[nBlocks-1]>=tmpArray2.length) tmpArray2 = new int[lastInBlock_int[nBlocks-1]+1];
			System.arraycopy(intBlocks[nBlocks-1],0,tmpArray2,0,lastInBlock_int[nBlocks-1]+1);
			rightUnique=isBlockUnique[0];
			last2=lastInBlock_int[nBlocks-1];
			if ((leftEnd1_last>0 && leftEnd2_last==0) || (rightEnd1_last>0 && rightEnd2_last==0 && leftEnd2_last>=0)) {
				nAttempted++;
				success=breakReads_checkDisambiguation_impl(mode,mode?Reads.breakReads_new2old[r][1]-Reads.breakReads_new2old[r-1][2]-1:0,leftUnique,leftLength,rightUnique,boundaries[0],lengthThreshold,leftEnd1,leftEnd1_last,leftEnd1_str,leftEnd2_str,rightEnd1,rightEnd1_last,str1,str2,leftEnd2,leftEnd2_last,rightEnd2,rightEnd2_last,bw);
				if (!success) {
					nRolledBack++;
					System.err.println("breakReads_checkDisambiguation> rolled back: "+str1+" <- "+str2);
				}
			}
			else {
				if (leftEnd2_str.length()!=0) bw.write(SEPARATOR_MAJOR+leftEnd2_str.toString()+"\n");
				bw.write(str2.substring(0,str2.lastIndexOf(SEPARATOR_MAJOR+"")));
			}
			tmpArray=leftEnd1; leftEnd1=tmpArray1; tmpArray1=tmpArray;
			leftEnd1_last=last1;
			leftEnd1_str.delete(0,leftEnd1_str.length());
			leftEnd1_str.append(str1.lastIndexOf(SEPARATOR_MAJOR+"")+1);
			tmpArray=leftEnd2; leftEnd2=tmpArray2; tmpArray2=tmpArray;
			leftEnd2_last=last2;
			leftEnd2_str.delete(0,leftEnd2_str.length());
			leftEnd2_str.append(str2.lastIndexOf(SEPARATOR_MAJOR+"")+1);
			leftLength=Reads.readLengths[r]-boundaries[nBlocks-2]-1;
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
		}
		if (leftEnd2_str.length()!=0) bw.write(SEPARATOR_MAJOR+leftEnd2_str.toString()+"\n");
		br1.close(); br2.close(); br3.close(); bw.close();
		System.err.println("Rolled back "+nRolledBack+" pairs out of "+nAttempted+" attempted ("+(((double)(nRolledBack*100))/nAttempted)+"%)");
	}
	
	
	/**
	 * The disambiguation of the left and right side of a low-quality region is rolled 
	 * back if it uses open characters, or if it uses characters whose length is 
	 * incompatible with the old read. If just one side was disambiguated, the procedure
	 * tries to set the other side to the same character.
	 *
	 * Remark: the procedure assumes that either the left or the right side were
	 * disambiguated, and that the side that was not disambiguated contains at least one
	 * character.
	 *
	 * Remark: if a low-quality regions is a replacement, the procedure assumes that only
	 * the character on the left/right side of it can occur inside the replaced interval.
	 * This is done for simplicity.
	 *
	 * Remark: unique and periodic characters are treated the same way as nonperiodic.
	 *
	 * @param lowQualityLength length of the low-quality interval between the new reads
	 * (discarded if low-quality intervals are assumed to be insertions);
	 * @param lengthThreshold tolerance in length comparisons;
	 * @param leftLength,rightEnd length of the left and right blocks of the new reads;
	 * @param leftUnique,rightUnique tell whether the left or right block is non-
	 * repetitive;
	 * @return FALSE iff the disambiguation is rejected.
	 */
	private static final boolean breakReads_checkDisambiguation_impl(boolean mode, int lowQualityLength, boolean leftUnique, int leftLength, boolean rightUnique, int rightLength, int lengthThreshold, int[] leftEnd1, int leftEnd1_last, StringBuilder leftEnd1_str, StringBuilder leftEnd2_str, int[] rightEnd1, int rightEnd1_last, String str1, String str2, int[] leftEnd2, int leftEnd2_last, int[] rightEnd2, int rightEnd2_last, BufferedWriter bw) throws IOException {
		int i;
		int leftCharacter, leftCharacterLength, rightCharacter, rightCharacterLength;
		final int availableLength = leftLength+(mode?lowQualityLength:0)+rightLength;

		if (leftUnique || rightUnique) {
			// Disambiguation accepted
			bw.write(SEPARATOR_MAJOR+""+leftEnd2_str.toString()); bw.newLine();
			bw.write(str2.substring(0,str2.lastIndexOf(SEPARATOR_MAJOR+"")));
			return true;
		}
		if (leftEnd2_last==0) {
			leftCharacter=leftEnd2[0];
			if (alphabet[leftCharacter].openEnd) {
				// Disambiguation rejected
				bw.write(SEPARATOR_MAJOR+leftEnd1_str.toString()); bw.newLine();
				bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
				return false;
			}
			leftCharacterLength=alphabet[leftCharacter].getLength();
			if (rightEnd2_last==0) rightCharacter=rightEnd2[0];
			else {
				rightCharacter=-1;
				for (i=0; i<=rightEnd1_last; i++) {
					if (rightEnd1[i]==leftCharacter) {
						rightCharacter=leftCharacter;
						break;
					}
				}
			}
			if (rightCharacter==leftCharacter) {
				if ( (mode && (Math.abs(leftCharacterLength,availableLength)<=lengthThreshold || (leftCharacterLength<<1)<=availableLength+lengthThreshold)) ||
				     (!mode && Math.abs(leftCharacterLength,availableLength)<=lengthThreshold)
				   ) {
					// Disambiguation accepted
					bw.write(SEPARATOR_MAJOR+""+rightCharacter); bw.newLine();
					bw.write(rightCharacter+str1.substring(str1.indexOf(SEPARATOR_MAJOR+""),str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return true;
				}
				else {
					// Disambiguation rejected
					bw.write(SEPARATOR_MAJOR+leftEnd1_str.toString()); bw.newLine();
					bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return false;
				}
			}
			else if (rightCharacter!=-1) {
				if (alphabet[rightCharacter].openStart) {
					// Disambiguation rejected
					bw.write(SEPARATOR_MAJOR+""+leftCharacter); bw.newLine();
					bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return false;
				}
				rightCharacterLength=alphabet[rightCharacter].getLength();
				if ( (mode && leftCharacterLength+rightCharacterLength<=availableLength+lengthThreshold) ||
				     (!mode && Math.abs(leftCharacterLength+rightCharacterLength,availableLength)<=lengthThreshold)
				   ) {
					// Disambiguation accepted
					bw.write(SEPARATOR_MAJOR+""+leftCharacter); bw.newLine();
					bw.write(rightCharacter+str1.substring(str1.indexOf(SEPARATOR_MAJOR+""),str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return true;
				}
				else {
					// Disambiguation rejected
					bw.write(SEPARATOR_MAJOR+leftEnd1_str.toString()); bw.newLine();
					bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return false;
				}
			}
			else {
				// Cannot disambiguate the right end
				bw.write(SEPARATOR_MAJOR+""+leftEnd1_str.toString()); bw.newLine();
				bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
				return false;
			}
		}
		else if (rightEnd2_last==0) {
			rightCharacter=rightEnd2[0];
			if (alphabet[rightCharacter].openStart) {
				// Disambiguation rejected
				bw.write(SEPARATOR_MAJOR+""+leftEnd1_str.toString()); bw.newLine();
				bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
				return false;
			}
			rightCharacterLength=alphabet[rightCharacter].getLength();
			if (leftEnd2_last==0) leftCharacter=leftEnd2[0];
			else {
				leftCharacter=-1;
				for (i=0; i<=leftEnd1_last; i++) {
					if (leftEnd1[i]==rightCharacter) {
						leftCharacter=rightCharacter;
						break;
					}
				}
			}
			if (leftCharacter==rightCharacter) {
				if ( (mode && (Math.abs(rightCharacterLength,availableLength)<=lengthThreshold || (rightCharacterLength<<1)<=availableLength+lengthThreshold)) ||
				     (!mode && Math.abs(rightCharacterLength,availableLength)<=lengthThreshold)
				   ) {
					// Disambiguation accepted
					bw.write(SEPARATOR_MAJOR+""+rightCharacter); bw.newLine();
					bw.write(rightCharacter+str1.substring(str1.indexOf(SEPARATOR_MAJOR+""),str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return true;
				}
				else {
					// Disambiguation rejected
					bw.write(SEPARATOR_MAJOR+""+leftEnd1_str.toString()); bw.newLine();
					bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return false;
				}
			}
			else if (leftCharacter!=-1) {
				if (alphabet[leftCharacter].openEnd) {
					// Disambiguation rejected
					bw.write(SEPARATOR_MAJOR+""+leftEnd1_str.toString()); bw.newLine();
					bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return false;
				}
				leftCharacterLength=alphabet[leftCharacter].getLength();
				if ( (mode && leftCharacterLength+rightCharacterLength<=availableLength+lengthThreshold) ||
				     (!mode && Math.abs(leftCharacterLength+rightCharacterLength,availableLength)<=lengthThreshold)
				   ) {
					// Disambiguation accepted
					bw.write(SEPARATOR_MAJOR+""+leftCharacter); bw.newLine();
					bw.write(rightCharacter+str1.substring(str1.indexOf(SEPARATOR_MAJOR+""),str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return true;
				}
				else {
					// Disambiguation rejected
					bw.write(SEPARATOR_MAJOR+leftEnd1_str.toString()); bw.newLine();
					bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
					return false;
				}
			}
			else {
				// Cannot disambiguate the left end
				bw.write(SEPARATOR_MAJOR+""+leftEnd1_str.toString()); bw.newLine();
				bw.write(str1.substring(0,str1.lastIndexOf(SEPARATOR_MAJOR+"")));
				return false;
			}
		}
		else return false;  // Never happens
	}
	
	
	
	
	
	
	
	
	// ------------------------ INACCURATE PERIODIC ENDPOINTS ----------------------------
	
	/**
	 * The aligner used to map the repeat library to the reads might terminate alignments
	 * early for periodic repeats, leaving gaps (regions not covered by any alignment) 
	 * between pairs of occurrences of periodic repeats, and between pairs of occurrences
	 * of a periodic and a non-periodic repeat, which are in reality consecutive in the 
	 * genome. We observed gaps of length <=500 bps (and even reaching some kbps) with 
	 * daligner on pbsim-simulated CLR reads from a simulated, fully-repetitive genome.
	 * The position and length of such gaps is not concordant across reads, and this
	 * breaks down our entire approach based on unique k-mers on the repeat alphabet.
	 * The approach breaks down in different ways: since several such gaps are considered
	 * unique signatures, the overlap graph becomes too complex. If we ask not to use any 
	 * non-repetitive block, then very few k-mers are considered a signature, since they 
	 * are too rare, and the overlap graph becomes too fragmented. We try to fix this by 
	 * closing the gaps in a concordant way across reads, propagating arbitrary decisions
	 * greedily using the read-read alignments that we already have.
	 *
	 * Remark: since this follows every alignment, it makes 1-mers and 2-mers (on the 
	 * repeat alphabet), that come from different regions of the genome, more similar, and 
	 * thus it reduces the number of unique 1-mers and 2-mers downstream, shifting the
	 * burden of unicity to longer k-mers.
	 *
	 * This procedure loads in $spacers$:
	 *
	 * 1. every pair of adjacent periodic blocks (both of them might contain non-periodic 
	 * repeats as well);
	 * 2.1 every short non-repetitive block, between two periodic blocks, such that one of 
	 * them contains non-periodic repeats;
	 * 2.2 every short non-repetitive block, between two periodic blocks, such that none
	 * of them contains non-periodic repeats;
	 * 3. every short non-repetitive block, between a periodic and a non-periodic block.
	 *
	 * An element of type 1, 2.1, 3 is called \emph{rigid}, since there is just one way of 
	 * assigning a breakpoint to it. Short non-repetitive blocks at the beginning/end of a
	 * read are loaded like 2.2.
	 *
	 * Remark: we need to load a spacer of type 3 in read X, because another read Y might 
	 * cover, with its beginning or end, some part of the spacer in read X, and some part 
	 * of the periodic block next to it in read X. The spacer in read Y might align only
	 * to spacers of type 3, in which case it would become disconnected in the spacers 
	 * graph, and the breakpoint decision taken in the connected component of the spacer 
	 * in read X would not propagate to the spacer in read Y.
	 * Another reason for loading spacers of type 3 is that the the non-repetitive part of
	 * the spacer might come from a periodic repeat that is completely different from the
	 * adjacent one, or from a non-periodic repeat, and using the spacers graph we could
	 * take a consistent decision on how to handle its breakpoints.
	 *
	 * Remark: an obvious limitation of this approach is the dependence on periodic 
	 * repeats. I.e. if read X contains a spacer of type 3, and read Y contains just its
	 * nonperiodic-nonrepetitive side because read Y ends, the instance on read Y does
	 * not become a spacer and, if we set the spacer to periodic in read X, we will never
	 * be able to set it to periodic in read Y. This would require propagating repeat 
	 * tags, and it's probably too laborious to implement.
	 *
	 * Remark: the procedure is sequential just for simplicity.
	 *
	 * Remark: the procedure assumes that $Reads.readLengths$ has already been loaded.
	 *
	 * Remark: the procedure needs the following arrays:
	 * translation_all, isBlockUnique_all, boundaries_all.
	 *
	 * @param maxBlockLength max length of a block in $translation_all$.
	 */
	public static final void loadSpacers(int maxSpacerLength, int maxBlockLength) {
		final int LAST_THREE_BITS = 0x00000007;
		boolean foundPeriodic, foundNonperiodic;
		int i, j, k, c;
		int cell, mask, length, nBlocks, nSpacers, readLength;
		final int nTranslatedReads = translated.length;
		boolean[] isBlockPeriodic, isBlockNonperiodic;
		int[] lengthHistogram;
		
		spacers = new Spacer[100];  // Arbitrary
		isBlockPeriodic = new boolean[100];  // Arbitrary
		isBlockNonperiodic = new boolean[100];  // Arbitrary
		nSpacers=0; nRigidSpacers=0;
		lengthHistogram = new int[11];  // Arbitrary, multiples of $IO.quantum$.
		Math.set(lengthHistogram,lengthHistogram.length-1,0);
		lastSpacer=-1;
		for (i=0; i<nTranslatedReads; i++) {
			nBlocks=translation_all[i].length;
			if (nBlocks<2) continue;
			if (isBlockPeriodic.length<nBlocks) isBlockPeriodic = new boolean[nBlocks];
			if (isBlockNonperiodic.length<nBlocks) isBlockNonperiodic = new boolean[nBlocks];
			for (j=0; j<nBlocks; j++) {
				foundPeriodic=false; foundNonperiodic=false;
				length=translation_all[i][j].length;
				for (k=0; k<length; k++) {
					c=translation_all[i][j][k];
					if (c>lastUnique && c<=lastPeriodic) {
						foundPeriodic=true;
						if (foundNonperiodic) break;
					}
					else if (c>lastPeriodic) {
						foundNonperiodic=true;
						if (foundPeriodic) break;
					}
				}
				isBlockPeriodic[j]=foundPeriodic; isBlockNonperiodic[j]=foundNonperiodic;
			}
			// First block
			cell=isBlockUnique_all[i][0]; mask=1;
			if ((cell&mask)!=0 && isBlockPeriodic[1] && !isBlockNonperiodic[1]) {
				length=boundaries_all[i][0];
				lengthHistogram[length/IO.quantum<lengthHistogram.length?length/IO.quantum:lengthHistogram.length-1]++;
				if (length<=maxSpacerLength) {
					lastSpacer++; nSpacers++;
					if (lastSpacer==spacers.length) {
						Spacer[] newArray = new Spacer[spacers.length<<1];
						System.arraycopy(spacers,0,newArray,0,spacers.length);
						spacers=newArray;
					}
					spacers[lastSpacer] = new Spacer(translated[i],0,boundaries_all[i][0],false,false,0);
				}
			}
			// Intermediate blocks		
			for (j=1; j<nBlocks-1; j++) {
				cell=isBlockUnique_all[i][j>>3]; mask=1<<(j&LAST_THREE_BITS);
				if ((cell&mask)!=0) {
					if (isBlockPeriodic[j-1] && isBlockPeriodic[j+1] && !(isBlockNonperiodic[j-1] && isBlockNonperiodic[j+1])) {
						length=boundaries_all[i][j]-boundaries_all[i][j-1];
						lengthHistogram[length/IO.quantum<lengthHistogram.length?length/IO.quantum:lengthHistogram.length-1]++;
						if (length<=maxSpacerLength) {
							lastSpacer++; nSpacers++;
							if (lastSpacer==spacers.length) {
								Spacer[] newArray = new Spacer[spacers.length<<1];
								System.arraycopy(spacers,0,newArray,0,spacers.length);
								spacers=newArray;
							}
							spacers[lastSpacer] = new Spacer(translated[i],boundaries_all[i][j-1],boundaries_all[i][j],isBlockNonperiodic[j-1],isBlockNonperiodic[j+1],j);
						}
					}
					else if ((isBlockNonperiodic[j-1] && (isBlockPeriodic[j+1]&&!isBlockNonperiodic[j+1])) || ((isBlockPeriodic[j-1]&&!isBlockNonperiodic[j-1]) && isBlockNonperiodic[j+1])) {
						length=boundaries_all[i][j]-boundaries_all[i][j-1];
						lengthHistogram[length/IO.quantum<lengthHistogram.length?length/IO.quantum:lengthHistogram.length-1]++;
						if (length<=maxSpacerLength) {
							lastSpacer++; nSpacers++;
							if (lastSpacer==spacers.length) {
								Spacer[] newArray = new Spacer[spacers.length<<1];
								System.arraycopy(spacers,0,newArray,0,spacers.length);
								spacers=newArray;
							}
							spacers[lastSpacer] = new Spacer(translated[i],boundaries_all[i][j-1],boundaries_all[i][j],isBlockNonperiodic[j-1],isBlockNonperiodic[j+1],j);
						}
					}
				}
				else if (isBlockPeriodic[j] && isBlockPeriodic[j-1]) {
					lastSpacer++; nRigidSpacers++;
					if (lastSpacer==spacers.length) {
						Spacer[] newArray = new Spacer[spacers.length<<1];
						System.arraycopy(spacers,0,newArray,0,spacers.length);
						spacers=newArray;
					}
					spacers[lastSpacer] = new Spacer(translated[i],boundaries_all[i][j-1],boundaries_all[i][j-1],true,true,-1);
				}
			}
			// Last block
			cell=isBlockUnique_all[i][(nBlocks-1)>>3]; mask=1<<((nBlocks-1)&LAST_THREE_BITS);
			if ((cell&mask)!=0 && isBlockPeriodic[nBlocks-2] && !isBlockNonperiodic[nBlocks-2]) {
				readLength=Reads.getReadLength(translated[i]);
				length=readLength-boundaries_all[i][nBlocks-2];
				lengthHistogram[length/IO.quantum<lengthHistogram.length?length/IO.quantum:lengthHistogram.length-1]++;
				if (length<=maxSpacerLength) {
					lastSpacer++; nSpacers++;
					if (lastSpacer==spacers.length) {
						Spacer[] newArray = new Spacer[spacers.length<<1];
						System.arraycopy(spacers,0,newArray,0,spacers.length);
						spacers=newArray;
					}
					spacers[lastSpacer] = new Spacer(translated[i],boundaries_all[i][nBlocks-2],readLength-1,false,false,nBlocks-1);
				}
			}
			else if (isBlockPeriodic[nBlocks-1] && isBlockPeriodic[nBlocks-2]) {
				lastSpacer++; nRigidSpacers++;
				if (lastSpacer==spacers.length) {
					Spacer[] newArray = new Spacer[spacers.length<<1];
					System.arraycopy(spacers,0,newArray,0,spacers.length);
					spacers=newArray;
				}
				spacers[lastSpacer] = new Spacer(translated[i],boundaries_all[i][nBlocks-2],boundaries_all[i][nBlocks-2],true,true,-1);
			}
		}
		System.err.println("Loaded "+nSpacers+" spacers ("+(((double)nSpacers)/nTranslatedReads)+" per translated read) and "+nRigidSpacers+" rigid spacers ("+(((double)nRigidSpacers)/nTranslatedReads)+" per translated read).");
		System.err.println("Histogram of all observed spacer lengths:");
		for (i=0; i<lengthHistogram.length; i++) System.err.println((i*IO.quantum)+": "+lengthHistogram[i]);		
	}
	
	
	/**
	 * @param array* assumed to contain valid positions in $alphabet$;
	 * @param tmpArray* temporary space, assumed to be large enough;
	 * @return if $arrays1,array2$ both contain some periodic character with the same 
	 * periodic repeat in the same orientation, all such distinct $(repeat,orientation)$
	 * pairs are stored in $tmpArray3[0..X]$ in sorted order, and $X$ is returned in
	 * output (a repeat $X$ in RC orientation is endoded as $-1-X$).
	 */
	private static final int samePeriod(int[] array1, int length1, int[] array2, int length2, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) {
		int i, j, c;
		int last1, last2;
		Character character;
		
		last1=-1;
		for (i=0; i<length1; i++) {
			c=array1[i];
			if (c>lastUnique && c<=lastPeriodic) {
				character=alphabet[c];
				tmpArray1[++last1]=character.orientation?character.repeat:-1-character.repeat;
			}
		}
		if (last1==-1) return -1;
		if (last1>0) {
			Arrays.sort(tmpArray1,0,last1+1);
			j=0;
			for (i=1; i<=last1; i++) {
				if (tmpArray1[i]==tmpArray1[j]) continue;
				j++;
				tmpArray1[j]=tmpArray1[i];
			}
			last1=j;
		}
		last2=-1;
		for (i=0; i<length2; i++) {
			c=array2[i];
			if (c>lastUnique && c<=lastPeriodic) {
				character=alphabet[c];
				tmpArray2[++last2]=character.orientation?character.repeat:-1-character.repeat;
			}
		}
		if (last2==-1) return -1;
		if (last2>0) {
			Arrays.sort(tmpArray2,0,last2+1);
			j=0;
			for (i=1; i<=last2; i++) {
				if (tmpArray2[i]==tmpArray2[j]) continue;
				j++;
				tmpArray2[j]=tmpArray2[i];
			}
			last2=j;
		}
		return Math.setIntersection(tmpArray1,0,last1,tmpArray2,0,last2,tmpArray3,0);
	}
	
	
	/**
	 * Variant of the procedure above.
	 */
	private static final boolean samePeriod(String[] array1, int length1, String[] array2, int length2, Character[] alphabet, int lastUnique, int lastPeriodic, int[] tmpArray1, int[] tmpArray2) {
		int i, c;
		int last1, last2;
		Character character;
		
		last1=-1;
		for (i=0; i<length1; i++) {
			c=Integer.parseInt(array1[i]);
			if (c<0) c=-1-c;			
			if (c>lastUnique && c<=lastPeriodic) {
				character=alphabet[c];
				tmpArray1[++last1]=character.orientation?character.repeat:-1-character.repeat;
			}
		}
		if (last1==-1) return false;
		if (last1>0) Arrays.sort(tmpArray1,0,last1+1);
		last2=-1;
		for (i=0; i<length2; i++) {
			c=Integer.parseInt(array2[i]);
			if (c<0) c=-1-c;
			if (c>lastUnique && c<=lastPeriodic) {
				character=alphabet[c];
				tmpArray2[++last2]=character.orientation?character.repeat:-1-character.repeat;
			}
		}
		if (last2==-1) return false;
		if (last2>0) Arrays.sort(tmpArray2,0,last2+1);
		return Math.nonemptyIntersection(tmpArray1,0,last1,tmpArray2,0,last2);
	}
	
	
	/**
	 * Creates an edge between two spacers iff at least one of them is not rigid, and 
	 * there is a read-read alignment in $alignmentsFile$ that makes the two spacers
	 * intersect when projected to the same read. An edge $(x,y)$ is represented with two
	 * integers in $spacerNeighbors[x]$: the first is either $y$ or $-1-y$, depending on
	 * the orientation of the alignment; the second is the quantity to add to $y.first$ to
	 * get the projection of $x.first$ on $y$'s read (might be negative). The entry in 
	 * $spacerNeighbors[y]$ for the same edge might have a different absolute value.
	 *
	 * Remark: the procedure initializes the field $breakpoint$ of every spacer, and marks
	 * bridging and inactive spacers. A spacer is \emph{bridging} if its two breakpoints 
	 * should be deleted, i.e. iff: (1) it has a periodic block to its left and to its 
	 * right; (2) both periodic blocks contain the same repeat R in the same orientation; 
	 * (3) the spacer (and some context to its left and right) aligns to a string that is
	 * fully contained inside a periodic block of another read; (4) such block also 
	 * contains R in the same orientation. A spacer is \emph{inactive} if its two 
	 * breakpoints are likely correct and should not be edited, i.e. iff it aligns to a 
	 * string that is fully contained inside a non-periodic repeat in another read. Note 
	 * that, given the two breakpoints of a spacer, there are just 3 possible operations: 
	 * delete both (bridging), keep both (inactive), and merge them into one (active 
	 * spacer). We don't consider a fourth operation of keeping both breakpoints but 
	 * changing their positions.
	 *
	 * Remark: some inactive spacers could be detected by just looking at the original
	 * version of the read translation file (before rare characters got removed). However,
	 * looking at all read-read alignments is needed anyway for spacers shorter than 
	 * $minAlignmentLength$, and for spacers where a non-periodic alignment never occured.
	 *
	 * Remark: the procedure is sequential just for simplicity.
	 *
	 * Remark: $spacerNeighbors[i]$ is sorted by decreasing alignment similarity for every
	 * $i$. The second integer is infinity when projecting to a rigid spacer.
	 *
	 * @return the number of spacers whose $breakpoint$ field has already been assigned.
	 */
	public static final int loadSpacerNeighbors(String alignmentsFile, int minAlignmentLength_readRepeat, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MIN_INTERSECTION = (3*IDENTITY_THRESHOLD)>>1;  // Arbitrary
		final int PERIODIC_CONTEXT = minAlignmentLength_readRepeat>>1;  // Arbitrary
		final int CAPACITY = 6;  // Arbitrary, multiple of 3.
		boolean readB_fullyUnique, readB_fullyContained;
		int i, j, k, p;
		int firstSpacer, readA, readB, nEdges, row, max, last, fromA, toA, fromB, toB, readA_length, readA_translatedIndex, readB_translatedIndex;
		int nSingletonSpacers_rigid, nSingletonSpacers_nonRigid_all, nSingletonSpacers_nonRigid_bridging;
		String str;
		BufferedReader br;
		Spacer tmpSpacer = new Spacer();
		Edge tmpEdge;
		Edge[] edges;
		
		// Loading edges
		for (i=0; i<=lastSpacer; i++) spacers[i].breakpoint=-1;
		spacerNeighbors = new double[lastSpacer+1][CAPACITY];
		lastSpacerNeighbor = new int[lastSpacer+1];
		Math.set(lastSpacerNeighbor,lastSpacer,-1);
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); firstSpacer=0; nEdges=0; row=0;
		while (str!=null) {
			row++;
			if (row%1000000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			readA=Alignments.readA-1;
			if (readA<spacers[firstSpacer].read) {
				str=br.readLine();
				continue;
			}
			while (firstSpacer<=lastSpacer && spacers[firstSpacer].read<readA) firstSpacer++;
			if (firstSpacer>lastSpacer) break;
			if (readA<spacers[firstSpacer].read) {
				str=br.readLine();
				continue;
			}
			readA_length=Reads.getReadLength(readA);
			readA_translatedIndex=Arrays.binarySearch(translated,readA);
			readB=Alignments.readB-1;
			readB_fullyUnique=Arrays.binarySearch(fullyUnique,readB)>=0;
			readB_fullyContained=Arrays.binarySearch(fullyContained,readB)>=0;
			if (!readB_fullyUnique && !readB_fullyContained) readB_translatedIndex=Arrays.binarySearch(translated,readB);
			else readB_translatedIndex=-1;
			for (i=firstSpacer; i<=lastSpacer; i++) {
				if (spacers[i].read!=readA) break;
				if ( (spacers[i].first==spacers[i].last && spacers[i].first>Alignments.startA+IDENTITY_THRESHOLD && spacers[i].first<Alignments.endA-IDENTITY_THRESHOLD) ||
					 ( spacers[i].first!=spacers[i].last &&
					   Intervals.isApproximatelyContained(spacers[i].first,spacers[i].last,Alignments.startA,Alignments.endA) &&
					   !Intervals.areApproximatelyIdentical(spacers[i].first,spacers[i].last,Alignments.startA,Alignments.endA)
					 )
				   ) {
					p=readHasSpacer(readB,firstSpacer,tmpSpacer);
					if (p>=0) {
						for (j=p; j<=lastSpacer; j++) {
							if (spacers[j].read!=readB) break;
							if ( (spacers[j].first==spacers[j].last && spacers[j].first>Alignments.startB+IDENTITY_THRESHOLD && spacers[j].first<Alignments.endB-IDENTITY_THRESHOLD) ||
								 ( spacers[j].first!=spacers[j].last &&
								   Intervals.isApproximatelyContained(spacers[j].first,spacers[j].last,Alignments.startB,Alignments.endB) &&
								   !Intervals.areApproximatelyIdentical(spacers[j].first,spacers[j].last,Alignments.startB,Alignments.endB)
								 )
							   ) nEdges+=loadSpacerNeighbors_impl(i,j,MIN_INTERSECTION,IDENTITY_THRESHOLD)?1:0;
						}
					}
				}
				
				// Initializing the $breakpoint$ field.				
				if (readB_fullyUnique) continue;			
				if ( spacers[i].first!=spacers[i].last && !spacers[i].isRigid() &&
					 Intervals.isApproximatelyContained(spacers[i].first,spacers[i].last,Alignments.startA,Alignments.endA) &&
					 !Intervals.areApproximatelyIdentical(spacers[i].first,spacers[i].last,Alignments.startA,Alignments.endA)
				   ) {
					last=translation_all[readA_translatedIndex].length-1;
					if (spacers[i].blockID>0 && spacers[i].blockID<last) p=samePeriod(translation_all[readA_translatedIndex][spacers[i].blockID-1],translation_all[readA_translatedIndex][spacers[i].blockID-1].length,translation_all[readA_translatedIndex][spacers[i].blockID+1],translation_all[readA_translatedIndex][spacers[i].blockID+1].length,tmpArray1,tmpArray2,tmpArray3);
					else p=-1;
				   	if (p>=0) {
						if (readB_fullyContained) {
						   	// Assuming that readB is fully periodic, for simplicity.
							spacers[i].breakpoint=Math.POSITIVE_INFINITY;
						}
						else {
							fromA=boundaries_all[readA_translatedIndex][spacers[i].blockID-1]-PERIODIC_CONTEXT;
							if (spacers[i].blockID>=2 && fromA<boundaries_all[readA_translatedIndex][spacers[i].blockID-2]) fromA=boundaries_all[readA_translatedIndex][spacers[i].blockID-2];
							else if (spacers[i].blockID==1 && fromA<0) fromA=0;
							toA=boundaries_all[readA_translatedIndex][spacers[i].blockID]+PERIODIC_CONTEXT;
							if (spacers[i].blockID<last-1 && toA>boundaries_all[readA_translatedIndex][spacers[i].blockID+1]) toA=boundaries_all[readA_translatedIndex][spacers[i].blockID+1];
							else if (spacers[i].blockID==last-1 && toA>readA_length-1) toA=readA_length-1;
							Alignments.projectIntersection(fromA,toA,tmpArray1);
							fromB=tmpArray1[0]; toB=tmpArray1[1];
							if (toB<=boundaries_all[readB_translatedIndex][0]+IDENTITY_THRESHOLD && blockContainsRepeat(readB_translatedIndex,0,tmpArray3,p,Alignments.orientation)) {
								spacers[i].breakpoint=Math.POSITIVE_INFINITY;
							}
							else {
								last=translation_all[readB_translatedIndex].length-2;
								for (j=1; j<=last; j++) {
									if ( ( Intervals.isApproximatelyContained(fromB,toB,boundaries_all[readB_translatedIndex][j-1],boundaries_all[readB_translatedIndex][j]) ||
										   Intervals.areApproximatelyIdentical(fromB,toB,boundaries_all[readB_translatedIndex][j-1],boundaries_all[readB_translatedIndex][j])
										 ) && blockContainsRepeat(readB_translatedIndex,j,tmpArray3,p,Alignments.orientation)
									   ) {
									   	spacers[i].breakpoint=Math.POSITIVE_INFINITY;
										break;
									}
								}
								if (fromB>=boundaries_all[readB_translatedIndex][last]-IDENTITY_THRESHOLD && blockContainsRepeat(readB_translatedIndex,last+1,tmpArray3,p,Alignments.orientation)) {
									spacers[i].breakpoint=Math.POSITIVE_INFINITY;
								}
							}
						}
					}
					if (spacers[i].breakpoint==-1 && !readB_fullyContained) { 
						// Marking as nonperiodic only if not already marked as periodic.
						// Remark: we assume that $readB_fullyContained$ means that readB
						// is fully periodic, for simplicity.
						Alignments.projectIntersection(spacers[i].first,spacers[i].last,tmpArray1);
						fromB=tmpArray1[0]; toB=tmpArray1[1];
						if (toB<=boundaries_all[readB_translatedIndex][0]+IDENTITY_THRESHOLD && !isBlockPeriodic(readB_translatedIndex,0) && isBlockNonperiodic(readB_translatedIndex,0)) spacers[i].breakpoint=Math.POSITIVE_INFINITY-1;
						else {
							last=translation_all[readB_translatedIndex].length-2;
							for (j=1; j<=last; j++) {
								if ( ( Intervals.isApproximatelyContained(fromB,toB,boundaries_all[readB_translatedIndex][j-1],boundaries_all[readB_translatedIndex][j]) ||
									   Intervals.areApproximatelyIdentical(fromB,toB,boundaries_all[readB_translatedIndex][j-1],boundaries_all[readB_translatedIndex][j])
									 ) && !isBlockPeriodic(readB_translatedIndex,j) && isBlockNonperiodic(readB_translatedIndex,j)
								   ) {
								   	spacers[i].breakpoint=Math.POSITIVE_INFINITY-1;
									break;
								}
							}
							if (fromB>=boundaries_all[readB_translatedIndex][last]-IDENTITY_THRESHOLD && !isBlockPeriodic(readB_translatedIndex,last+1) && isBlockNonperiodic(readB_translatedIndex,last+1)) spacers[i].breakpoint=Math.POSITIVE_INFINITY-1;
						}
					}
				}
			}
			str=br.readLine();
		}
		br.close();
		System.err.println("Loaded "+nEdges+" spacer edges");
		
		// Removing duplicates and sorting edges
		max=0;
		for (i=0; i<=lastSpacer; i++) max=Math.max(max,lastSpacerNeighbor[i]);
		max+=1;
		edges = new Edge[max];
		for (i=0; i<max; i++) edges[i] = new Edge();
		for (i=0; i<=lastSpacer; i++) {
			last=lastSpacerNeighbor[i];
			if (last==0) continue;
			nEdges=(last+1)/3;
			for (j=0; j<last; j+=3) edges[j/3].set(spacerNeighbors[i][j],spacerNeighbors[i][j+1],spacerNeighbors[i][j+2]);
			if (nEdges>1) {
				Edge.order=Edge.ORDER_NEIGHBOR;
				Arrays.sort(edges,0,nEdges);
				k=0;
				for (j=1; j<nEdges; j++) {
					if (edges[j].neighbor!=edges[k].neighbor) {
						k++;
						tmpEdge=edges[k]; edges[k]=edges[j]; edges[j]=tmpEdge;
					}
					else edges[k].diffs=Math.min(edges[k].diffs,edges[j].diffs);
				}
				nEdges=k+1;
			}
			if (nEdges>1) {
				Edge.order=Edge.ORDER_DIFFS;
				Arrays.sort(edges,0,nEdges);
			}
			k=-1;
			for (j=0; j<nEdges; j++) {
				spacerNeighbors[i][++k]=edges[j].neighbor;
				spacerNeighbors[i][++k]=edges[j].offset;
			}
			lastSpacerNeighbor[i]=k;
		}
		
		// Computing statistics
		nBridgingSpacers=0; nInactiveSpacers=0;
		nSingletonSpacers_rigid=0; nSingletonSpacers_nonRigid_all=0; nSingletonSpacers_nonRigid_bridging=0;
		for (i=0; i<=lastSpacer; i++) {
			if (spacers[i].breakpoint==Math.POSITIVE_INFINITY) nBridgingSpacers++;
			if (spacers[i].breakpoint==Math.POSITIVE_INFINITY-1) nInactiveSpacers++;
			if (lastSpacerNeighbor[i]>=0) continue;
			if (spacers[i].isRigid()) nSingletonSpacers_rigid++;
			else {
				nSingletonSpacers_nonRigid_all++;
				if (spacers[i].breakpoint==Math.POSITIVE_INFINITY) nSingletonSpacers_nonRigid_bridging++;
			}
		}
		System.err.println("Total singleton spacers: "+(nSingletonSpacers_rigid+nSingletonSpacers_nonRigid_all)+" ("+((100.0*(nSingletonSpacers_rigid+nSingletonSpacers_nonRigid_all))/(lastSpacer+1))+"%)");
		System.err.println("Rigid singleton spacers: "+nSingletonSpacers_rigid+" ("+((100.0*nSingletonSpacers_rigid)/(nSingletonSpacers_rigid+nSingletonSpacers_nonRigid_all))+"%)");
		System.err.println("Non-rigid singleton spacers: "+nSingletonSpacers_nonRigid_all+" ("+((100.0*nSingletonSpacers_nonRigid_all)/(nSingletonSpacers_rigid+nSingletonSpacers_nonRigid_all))+"%)");
		System.err.println("Non-rigid singleton spacers that are bridging: "+nSingletonSpacers_nonRigid_bridging+" ("+((100.0*nSingletonSpacers_nonRigid_bridging)/nSingletonSpacers_nonRigid_all)+"%)");
		return nBridgingSpacers;
	}
	
	
	/**
	 * @param readID row of $translation_all$.
	 */
	private static final boolean isBlockPeriodic(int readID, int blockID) {
		final int last = translation_all[readID][blockID].length-1;
		int i, c;
		
		for (i=0; i<=last; i++) {
			c=translation_all[readID][blockID][i];
			if (c>lastUnique && c<=lastPeriodic) return true;
		}
		return false;
	}
	
	
	/**
	 * @param readID row of $translation_all$;
	 * @param repeatIDs assumed to be sorted; a repeat $X$ in RC orientation is encoded as
	 * $-1-X$;
	 * @return TRUE iff the block contains a character with the same repeat ID and the 
	 * same orientation (if $sameOrientation=TRUE$) or with opposite orientation (if 
	 * $sameOrientation=FALSE$) in $repeatIDs[0..last]$.
	 */
	private static final boolean blockContainsRepeat(int readID, int blockID, int[] repeatIDs, int last, boolean sameOrientation) {
		final int lastPrime = translation_all[readID][blockID].length-1;
		int i, c;
		
		for (i=0; i<=lastPrime; i++) {
			c=translation_all[readID][blockID][i];
			if (c==lastAlphabet+1) continue;
			if (sameOrientation) c=alphabet[c].orientation?alphabet[c].repeat:-1-alphabet[c].repeat;
			else c=alphabet[c].orientation?-1-alphabet[c].repeat:alphabet[c].repeat;
			if (Arrays.binarySearch(repeatIDs,0,last+1,c)>=0) return true;
		}
		return false;
	}
	
	
	/**
	 * @param readID row of $translation_all$.
	 */
	private static final boolean isBlockNonperiodic(int readID, int blockID) {
		final int last = translation_all[readID][blockID].length-1;
		int i, c;
		
		for (i=0; i<=last; i++) {
			c=translation_all[readID][blockID][i];
			if (c>lastPeriodic) return true;
		}
		return false;
	}
	
	
	private static final boolean loadSpacerNeighbors_impl(int spacerAID, int spacerBID, int minIntersection, int identityThreshold) {
		boolean rigidA, rigidB, addEdge;
		int startA, endA, startB, endB, offsetAB, offsetBA;
		double ratio;
		final Spacer spacerA = spacers[spacerAID];
		final Spacer spacerB = spacers[spacerBID];
		
		// Evaluating conditions
		rigidA=spacerA.isRigid(); rigidB=spacerB.isRigid();
		if (rigidA && rigidB) return false;
		if (rigidA) {
			ratio=((double)(Alignments.endB-Alignments.startB+1))/(Alignments.endA-Alignments.startA+1);
			if (spacerA.first==spacerA.last || spacerA.rigidLeft) {
				if (Alignments.orientation) startB=Alignments.startB+(int)((spacerA.first-Alignments.startA)*ratio);
				else startB=Alignments.endB-(int)((spacerA.first-Alignments.startA)*ratio);
				addEdge=startB>=spacerB.first-identityThreshold && startB<=spacerB.last+identityThreshold;
				offsetAB=startB-spacerB.first; offsetBA=Math.POSITIVE_INFINITY;
			}
			else {
				if (Alignments.orientation) endB=Alignments.startB+(int)((spacerA.last-Alignments.startA)*ratio);
				else endB=Alignments.endB-(int)((spacerA.last-Alignments.startA)*ratio);
				addEdge=endB>=spacerB.first-identityThreshold && endB<=spacerB.last+identityThreshold;
				offsetAB=endB-spacerB.first; offsetBA=Math.POSITIVE_INFINITY;
			}
		}
		else if (rigidB) {
			ratio=((double)(Alignments.endA-Alignments.startA+1))/(Alignments.endB-Alignments.startB+1);
			if (spacerB.first==spacerB.last || spacerB.rigidLeft) {
				if (Alignments.orientation) startA=Alignments.startA+(int)((spacerB.first-Alignments.startB)*ratio);
				else startA=Alignments.endA-(int)((spacerB.first-Alignments.startB)*ratio);
				addEdge=startA>=spacerA.first-identityThreshold && startA<=spacerA.last+identityThreshold;
				offsetBA=startA-spacerA.first; offsetAB=Math.POSITIVE_INFINITY;
			}
			else {
				if (Alignments.orientation) endA=Alignments.startA+(int)((spacerB.last-Alignments.startB)*ratio);
				else endA=Alignments.endA-(int)((spacerB.last-Alignments.startB)*ratio);
				addEdge=endA>=spacerA.first-identityThreshold && endA<=spacerA.last+identityThreshold;
				offsetBA=endA-spacerA.first; offsetAB=Math.POSITIVE_INFINITY;
			}
		}
		else {
			addEdge=false;
			ratio=((double)(Alignments.endB-Alignments.startB+1))/(Alignments.endA-Alignments.startA+1);
			if (Alignments.orientation) {
				startB=Alignments.startB+(int)((spacerA.first-Alignments.startA)*ratio);
				endB=Alignments.startB+(int)((spacerA.last-Alignments.startA)*ratio);
			}
			else {
				startB=Alignments.endB-(int)((spacerA.last-Alignments.startA)*ratio);
				endB=Alignments.endB-(int)((spacerA.first-Alignments.startA)*ratio);
			}
			if (Intervals.intersectionLength(startB,endB,spacerB.first,spacerB.last)>=minIntersection) addEdge=true;
			offsetAB=startB-spacerB.first;
			ratio=((double)(Alignments.endA-Alignments.startA+1))/(Alignments.endB-Alignments.startB+1);
			if (Alignments.orientation) {
				startA=Alignments.startA+(int)((spacerB.first-Alignments.startB)*ratio);
				endA=Alignments.startA+(int)((spacerB.last-Alignments.startB)*ratio);
			}
			else {
				startA=Alignments.endA-(int)((spacerB.last-Alignments.startB)*ratio);
				endA=Alignments.endA-(int)((spacerB.first-Alignments.startB)*ratio);
			}
			if (Intervals.intersectionLength(startA,endA,spacerA.first,spacerA.last)>=minIntersection) addEdge=true;
			offsetBA=startA-spacerA.first;
		}
		if (!addEdge) return false;
		
		// Adding edges
		if (lastSpacerNeighbor[spacerAID]+2>=spacerNeighbors[spacerAID].length) {
			double[] newArray = new double[spacerNeighbors[spacerAID].length<<1];
			System.arraycopy(spacerNeighbors[spacerAID],0,newArray,0,spacerNeighbors[spacerAID].length);
			spacerNeighbors[spacerAID]=newArray;
		}
		ratio=((double)(Alignments.diffs<<1))/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2);
		spacerNeighbors[spacerAID][++lastSpacerNeighbor[spacerAID]]=Alignments.orientation?spacerBID:-1-spacerBID;
		spacerNeighbors[spacerAID][++lastSpacerNeighbor[spacerAID]]=offsetAB;
		spacerNeighbors[spacerAID][++lastSpacerNeighbor[spacerAID]]=ratio;
		if (lastSpacerNeighbor[spacerBID]+2>=spacerNeighbors[spacerBID].length) {
			double[] newArray = new double[spacerNeighbors[spacerBID].length<<1];
			System.arraycopy(spacerNeighbors[spacerBID],0,newArray,0,spacerNeighbors[spacerBID].length);
			spacerNeighbors[spacerBID]=newArray;
		}
		spacerNeighbors[spacerBID][++lastSpacerNeighbor[spacerBID]]=Alignments.orientation?spacerAID:-1-spacerAID;
		spacerNeighbors[spacerBID][++lastSpacerNeighbor[spacerBID]]=offsetBA;
		spacerNeighbors[spacerBID][++lastSpacerNeighbor[spacerBID]]=ratio;
		return true;
	}
	
	
	/**
	 * In DOT format
	 */
	public static final void printSpacerNeighbors(String path) throws IOException {
		int i, j;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(path));
		bw.write("digraph G {"); bw.newLine();
		for (i=0; i<=lastSpacer; i++) {
			bw.write(i+" [isRigid=\""+spacers[i].isRigid()+"\", rigidLeft=\""+spacers[i].rigidLeft+"\", rigidRight=\""+spacers[i].rigidRight+"\", read=\""+spacers[i].read+"\", first=\""+spacers[i].first+"\", last=\""+spacers[i].last+"\", breakpoint=\""+spacers[i].breakpoint+"\"];");
			bw.newLine();
			for (j=0; j<=lastSpacerNeighbor[i]; j+=2) {
				if (spacerNeighbors[i][j]>=0) bw.write(i+" -> "+(int)spacerNeighbors[i][j]+" [orientation=\"true\", shift=\""+(int)spacerNeighbors[i][j+1]+"\"];");
				else bw.write(i+" -> "+(int)(-1-spacerNeighbors[i][j])+" [orientation=\"false\", shift=\""+(int)spacerNeighbors[i][j+1]+"\"];");
				bw.newLine();
			}
		}
		bw.write("}"); bw.newLine();
		bw.close();
	}
	
	
	private static class Edge implements Comparable {
		public static final int ORDER_NEIGHBOR = 0;
		public static final int ORDER_DIFFS = 1;
		public static int order;
		
		public double neighbor, offset, diffs;
		
		public Edge() { }
		
		public Edge(double n, double o, double d) {
			set(n,o,d);
		}
		
		public void set(double n, double o, double d) {
			this.neighbor=n; this.offset=o; this.diffs=d;
		}
		
		public int compareTo(Object other) {
			Edge otherEdge = (Edge)other;
			if (order==ORDER_NEIGHBOR) {
				if (neighbor<otherEdge.neighbor) return -1;
				else if (neighbor>otherEdge.neighbor) return 1;
			}
			else if (order==ORDER_DIFFS) {
				if (diffs<otherEdge.diffs) return -1;
				else if (diffs>otherEdge.diffs) return 1;
			}
			return 0;
		}
	}
	
	
	/**
	 * Identical to $readInArray()$. 
	 *
	 * @param position a position in $spacers$ from which to start the search; must be the
	 * first element of a read block; might be greater than $lastSpacer$;
	 * @return the first position in $spacers$ that belongs to $read$, or -1 if $read$ has
	 * no spacer.
	 */
	private static final int readHasSpacer(int read, int position, Spacer tmpSpacer) {
		final int MAX_DIFF = 100;  // Arbitrary
		int i, j;
		int value, out;
		
		if (position<=lastSpacer) {
			value=spacers[position].read;
			if (value==read) return position;
			else if (read<value) {
				if (read>=value-MAX_DIFF) {
					for (i=position-1; i>=0; i--) {
						if (spacers[i].read<read) break;
						if (spacers[i].read==read) {
							out=i;
							for (j=i-1; j>=0; j--) {
								if (spacers[j].read!=read) break;
								out=j;
							}
							return out;
						}
					}
				}
				else {
					tmpSpacer.read=read;
					i=Arrays.binarySearch(spacers,0,position,tmpSpacer);
					if (i>=0) {
						out=i;
						for (j=i-1; j>=0; j--) {
							if (spacers[j].read!=read) break;
							out=j;
						}
						return out;
					}
				}
			}
			else {
				if (read<=value+MAX_DIFF) {
					for (i=position+1; i<=lastSpacer; i++) {
						if (spacers[i].read>read) break;
						else if (spacers[i].read==read) return i;
					}
				}
				else {
					tmpSpacer.read=read;
					i=Arrays.binarySearch(spacers,position+1,lastSpacer+1,tmpSpacer);
					if (i>=0) {
						out=i;
						for (j=i-1; j>=0; j--) {
							if (spacers[j].read!=read) break;
							out=j;
						}
						return out;
					}
				}
			}
		}
		else {
			tmpSpacer.read=read;
			i=Arrays.binarySearch(spacers,0,lastSpacer+1,tmpSpacer);
			if (i>=0) {
				out=i;
				for (j=i-1; j>=0; j--) {
					if (spacers[j].read!=read) break;
					out=j;
				}
				return out;
			}
		}
		return -1;
	}
	
	
	/**
	 * Assigns a breakpoint to every spacer, by first propagating the breakpoints of rigid
	 * spacers depth-first, and then by propagating arbitrary decisions from every spacer
	 * that was not reached from a rigid spacer.
	 *
	 * Remark: one could alternatively propagate breadth-first. The traversal order is not
	 * likely to matter in practice. Depth-first makes sense since edges are sorted by
	 * decreasing similarity (but one could put all edges in a priority queue instead). 
	 * Doing a global alignment of all the reads that have overlapping spacers is
	 * impractical.
	 *
	 * Remark: before propagating normal spacers, the procedure propagates bridging
	 * spacers and inactive spacers. This is because, if a spacer is bridging (i.e. fully
	 * belonging to a periodic repeat), it is likely that every other spacer that aligns 
	 * to it is bridging. And if a spacer is inactive, it is likely that every spacer that
	 * aligns to it has both breakpoints at approx. the same positions and is inactive as
	 * well.
	 *
	 * Remark: the procedure assumes that $Reads.readLengths$ has already been loaded.
	 *
	 * @param nAlreadyAssigned number of spacers whose $breakpoint$ field has already been 
	 * decided by $loadSpacerNeighbors()$.
	 */
	public static final void assignBreakpoints() throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		boolean orientation, propagate;
		int i, j;
		int top, currentSpacer, neighbor, offset, breakpoint, readLength, nAssigned, nAlreadyAssigned;
		
		// Propagating from bridging spacers
		if (stack==null || stack.length<lastSpacer+1) stack = new int[lastSpacer+1];
		nAlreadyAssigned=0;
		if (nBridgingSpacers!=0) {
			for (i=0; i<=lastSpacer; i++) {
				if (spacers[i].breakpoint!=Math.POSITIVE_INFINITY) continue;
				top=0; stack[0]=i;
				while (top>=0) {
					currentSpacer=stack[top--];
					for (j=0; j<=lastSpacerNeighbor[currentSpacer]; j+=2) {
						neighbor=(int)spacerNeighbors[currentSpacer][j];
						if (neighbor<0) neighbor=-1-neighbor;
						if (spacers[neighbor].breakpoint==Math.POSITIVE_INFINITY-2) continue;
						spacers[neighbor].breakpoint=Math.POSITIVE_INFINITY-2;
						stack[++top]=neighbor;
					}
				}
			}
			nBridgingSpacers=0;
			for (i=0; i<=lastSpacer; i++) {
				if (spacers[i].breakpoint==Math.POSITIVE_INFINITY-2) {
					spacers[i].breakpoint=Math.POSITIVE_INFINITY;
					nBridgingSpacers++;
				}
				else if (spacers[i].breakpoint==Math.POSITIVE_INFINITY) nBridgingSpacers++;
			}
			nAlreadyAssigned=nBridgingSpacers;
		}
		
		// Propagating from inactive spacers
		if (nInactiveSpacers!=0) {
			for (i=0; i<=lastSpacer; i++) {
				if (spacers[i].breakpoint!=Math.POSITIVE_INFINITY-1) continue;
				top=0; stack[0]=i;
				while (top>=0) {
					currentSpacer=stack[top--];
					for (j=0; j<=lastSpacerNeighbor[currentSpacer]; j+=2) {
						neighbor=(int)spacerNeighbors[currentSpacer][j];
						if (neighbor<0) neighbor=-1-neighbor;
						if (spacers[neighbor].breakpoint==Math.POSITIVE_INFINITY-2 || spacers[neighbor].breakpoint==Math.POSITIVE_INFINITY) continue;
						spacers[neighbor].breakpoint=Math.POSITIVE_INFINITY-2;
						stack[++top]=neighbor;
					}
				}
			}
			nInactiveSpacers=0;
			for (i=0; i<=lastSpacer; i++) {
				if (spacers[i].breakpoint==Math.POSITIVE_INFINITY-2) {
					spacers[i].breakpoint=Math.POSITIVE_INFINITY-1;
					nInactiveSpacers++;
				}
				else if (spacers[i].breakpoint==Math.POSITIVE_INFINITY-1) nInactiveSpacers++;
			}
			nAlreadyAssigned+=nInactiveSpacers;
		}
		
		// Deactivating every spacer with no neighbor
		nAssigned=0;
		for (i=0; i<=lastSpacer; i++) {
			if (spacers[i].breakpoint!=-1) continue;
			if (lastSpacerNeighbor[i]==-1) {
				spacers[i].breakpoint=Math.POSITIVE_INFINITY-1;
				nAssigned++;
			}
		}
		
		// Propagating from rigid spacers
		if (nRigidSpacers!=0) {
			for (i=0; i<=lastSpacer; i++) {
				if (!spacers[i].isRigid() || spacers[i].breakpoint!=-1) continue;
				spacers[i].setBreakpoint(); nAssigned++;
				top=0; stack[0]=i;
				while (top>=0) {
					currentSpacer=stack[top--];
					for (j=0; j<=lastSpacerNeighbor[currentSpacer]; j+=2) {
						neighbor=(int)spacerNeighbors[currentSpacer][j];
						if (neighbor<0) {
							neighbor=-1-neighbor;
							orientation=false;
						}
						else orientation=true;
						if (spacers[neighbor].isRigid() || spacers[neighbor].breakpoint!=-1) continue;						
						offset=(int)spacerNeighbors[currentSpacer][j+1];
						breakpoint=spacers[neighbor].first+offset+(orientation?spacers[currentSpacer].breakpoint-spacers[currentSpacer].first:spacers[currentSpacer].last-spacers[currentSpacer].breakpoint);
						if (breakpoint>=spacers[neighbor].first-IDENTITY_THRESHOLD && breakpoint<=spacers[neighbor].last+IDENTITY_THRESHOLD) {
							breakpoint=Math.max(breakpoint,spacers[neighbor].first);
							breakpoint=Math.min(breakpoint,spacers[neighbor].last);
							spacers[neighbor].breakpoint=breakpoint;
							stack[++top]=neighbor;
							nAssigned++;
						}
					}
				}
			}
		}
		
		// Propagating from non-rigid spacers that are strictly inside a read
		if (nAssigned+nAlreadyAssigned<lastSpacer+1) {
			for (i=0; i<=lastSpacer; i++) {
				if (spacers[i].breakpoint!=-1) continue;
				readLength=Reads.getReadLength(spacers[i].read);
				if (spacers[i].first==0 || spacers[i].last==readLength-1) continue;
				spacers[i].setBreakpoint(); nAssigned++;
				top=0; stack[0]=i;
				while (top>=0) {
					currentSpacer=stack[top--];
					for (j=0; j<=lastSpacerNeighbor[currentSpacer]; j+=2) {
						neighbor=(int)spacerNeighbors[currentSpacer][j];
						if (neighbor<0) {
							neighbor=-1-neighbor;
							orientation=false;
						}
						else orientation=true;
						if (spacers[neighbor].breakpoint!=-1) continue;
						offset=(int)spacerNeighbors[currentSpacer][j+1];
						breakpoint=spacers[neighbor].first+offset+(orientation?spacers[currentSpacer].breakpoint-spacers[currentSpacer].first:spacers[currentSpacer].last-spacers[currentSpacer].breakpoint);
						if (breakpoint>=spacers[neighbor].first-IDENTITY_THRESHOLD && breakpoint<=spacers[neighbor].last+IDENTITY_THRESHOLD) {
							breakpoint=Math.max(breakpoint,spacers[neighbor].first);
							breakpoint=Math.min(breakpoint,spacers[neighbor].last);
							spacers[neighbor].breakpoint=breakpoint;
							stack[++top]=neighbor;
							nAssigned++;
						}
					}
				}
			}
		}
		
		// Propagating from the remaining non-rigid spacers
		if (nAssigned+nAlreadyAssigned<lastSpacer+1) {
			for (i=0; i<=lastSpacer; i++) {
				if (spacers[i].breakpoint!=-1) continue;
				spacers[i].setBreakpoint(); nAssigned++;
				top=0; stack[0]=i;
				while (top>=0) {
					currentSpacer=stack[top--];
					for (j=0; j<=lastSpacerNeighbor[currentSpacer]; j+=2) {
						neighbor=(int)spacerNeighbors[currentSpacer][j];
						if (neighbor<0) {
							neighbor=-1-neighbor;
							orientation=false;
						}
						else orientation=true;
						if (spacers[neighbor].breakpoint!=-1) continue;
						offset=(int)spacerNeighbors[currentSpacer][j+1];
						breakpoint=spacers[neighbor].first+offset+(orientation?spacers[currentSpacer].breakpoint-spacers[currentSpacer].first:spacers[currentSpacer].last-spacers[currentSpacer].breakpoint);
						if (breakpoint>=spacers[neighbor].first-IDENTITY_THRESHOLD && breakpoint<=spacers[neighbor].last+IDENTITY_THRESHOLD) {
							breakpoint=Math.max(breakpoint,spacers[neighbor].first);
							breakpoint=Math.min(breakpoint,spacers[neighbor].last);
							spacers[neighbor].breakpoint=breakpoint;
							stack[++top]=neighbor;
							nAssigned++;
						}
					}
				}
			}
		}
		
		if (nAssigned+nAlreadyAssigned!=lastSpacer+1) {
			System.err.println("fixSpacers> ERROR: assigned "+(nAssigned+nAlreadyAssigned)+" spacers out of "+(lastSpacer+1));
			System.err.println("Spacers that were not assigned:");
			for (int x=0; x<=lastSpacer; x++) {
				if (spacers[x].breakpoint==-1) System.err.println(spacers[x]);
			}
			System.exit(1);
		}
	}
	
	
	/**
	 * Stores $spacers$ to multiple files, each containing reads up to $lastRead[i]$
	 * (included, zero-based).
	 */
	public static final void serializeSpacers(String prefix, int[] lastRead) throws IOException {
		int i, j;
		BufferedWriter bw;
		
		j=0;
		bw = new BufferedWriter(new FileWriter(prefix+j+".txt"));
		for (i=0; i<=lastSpacer; i++) {
			while (spacers[i].read>lastRead[j]) {
				bw.close();
				j++;
				bw = new BufferedWriter(new FileWriter(prefix+j+".txt"));
			}
			spacers[i].serialize(bw);
		}
		bw.close();
	}
	
	
	public static final void deserializeSpacers(String spacersFile, int nSpacers) throws IOException {
		int i;
		String str;
		BufferedReader br;
		
		spacers = new Spacer[nSpacers];
		br = new BufferedReader(new FileReader(spacersFile));
		str=br.readLine(); i=0;
		while (str!=null) {
			spacers[i] = new Spacer();
			spacers[i].deserialize(str);
			str=br.readLine(); i++;
		}
		br.close();
		lastSpacer=nSpacers-1;
	}
	
	
	/**
	 * Alters an existing read translation using spacer breakpoints. Every new character 
	 * instance that is not already in $alphabet$ is written to $bw$. Every character in
	 * $alphabet$ that is used in the new translation is marked in $used$.
	 *
	 * Remark: short non-periodic blocks between a periodic and a non-periodic block are 
	 * merged with their adjacent periodic block.
	 *
	 * Remark: non-periodic repeats do not change after this procedure completes.
	 *
	 * @param used same size as $alphabet$;
	 * @param spacersCursor current position in $spacers$;
	 * @param isBlock* temporary space, with a number of cells at least equal to the
	 * number of blocks in the translation;
	 * @param tmpArray* temporary space, with a number of cells at least equal to the max
	 * number of elements in a block of the translation;
	 * @return the new value of $spacersCursor$.
	 */
	public static final int fixPeriodicEndpoints_collectCharacterInstances(int readID, int spacersCursor, String read2characters, String read2boundaries, int readLength, int maxSpacerLength, BufferedWriter bw, boolean[] used, boolean[] isBlockPeriodic, boolean[] isBlockNonperiodic, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) throws IOException {
		final int CAPACITY = 10;  // Arbitrary
		final int QUANTUM = IO.quantum;
		boolean found, foundPeriodic, foundNonperiodic, isUnique;
		int i, j, c;
		int length, nBlocks, last;
		
		if (read2characters.length()==0) return spacersCursor;
		if (tmpCharacter==null) tmpCharacter = new Character();
		if (leftCharacters==null) {
			leftCharacters = new Character[CAPACITY];
			for (i=0; i<leftCharacters.length; i++) leftCharacters[i] = new Character();
		}
		loadBoundaries(read2boundaries);
		nBlocks=loadBlocks(read2characters);
		if (tmpBoolean==null || tmpBoolean.length<nBlocks) tmpBoolean = new boolean[nBlocks];
		Math.set(tmpBoolean,nBlocks-1,false);
		loadIntBlocks(nBlocks,boundaries,readLength,tmpCharacter);
		if (nBlocks==1) {
			for (j=0; j<=lastInBlock_int[0]; j++) used[intBlocks[0][j]]=true;
			return spacersCursor;
		}
		for (i=0; i<nBlocks; i++) {
			foundPeriodic=false; foundNonperiodic=false;
			last=lastInBlock_int[i];
			for (j=0; j<=last; j++) {
				c=intBlocks[i][j];
				if (c>lastUnique && c<=lastPeriodic) {
					foundPeriodic=true;
					if (foundNonperiodic) break;
				}
				else if (c>lastPeriodic) {
					foundNonperiodic=true;
					if (foundPeriodic) break;
				}
			}
			isBlockPeriodic[i]=foundPeriodic; isBlockNonperiodic[i]=foundNonperiodic;
		}
		
		// First block. Remark: if the block is unique and shrinks, it does not create a
		// new unique character in the alphabet.
		lastLeft=-1;
		if (isBlockUnique[0] && isBlockPeriodic[1] && !isBlockNonperiodic[1] && boundaries[0]<=maxSpacerLength) {
			while (spacersCursor<=lastSpacer && spacers[spacersCursor].read<readID) spacersCursor++;
			if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>0) {
				System.err.println("fixPeriodicEndpoints_collectCharacterInstances> ERROR (1): spacer not found: "+readID+"[0.."+boundaries[0]+"]");
				System.exit(1);
			}
			if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY) {
				length=spacers[spacersCursor].last-spacers[spacersCursor].first+1;
				leftCharacters_load(false,1,length,true,1==nBlocks-1);
				tmpBoolean[0]=true;
			}
			else if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) { /* NOP */ }
			else {
				length=spacers[spacersCursor].breakpoint>QUANTUM?spacers[spacersCursor].last-spacers[spacersCursor].breakpoint:boundaries[0];
				leftCharacters_load(false,1,length,spacers[spacersCursor].breakpoint<=QUANTUM,1==nBlocks-1);
				if (spacers[spacersCursor].breakpoint<=QUANTUM) tmpBoolean[0]=true;
			}
			i=2;
		}
		else i=1;
		
		// Intermediary blocks
		while (i<=nBlocks-2) {
			if (!isBlockUnique[i] || boundaries[i]-boundaries[i-1]>maxSpacerLength) {
				if (lastLeft!=-1) {
					leftCharacters_setOpen(false,false);
					leftCharacters_clear(used,bw,QUANTUM);
					tmpBoolean[i-1]=true;
				}
				i+=(isBlockUnique[i]||isBlockPeriodic[i+1])?2:1;
				continue;
			}
			if (isBlockPeriodic[i-1] && isBlockPeriodic[i+1]) {
				if (isBlockNonperiodic[i-1] && isBlockNonperiodic[i+1]) {
					if (lastLeft!=-1) {
						leftCharacters_setOpen(false,false);
						leftCharacters_clear(used,bw,QUANTUM);
						tmpBoolean[i-1]=true;
					}
					i+=2;
				}
				else if (isBlockNonperiodic[i-1] || isBlockNonperiodic[i+1]) {
					while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[i-1])) spacersCursor++;
					if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>boundaries[i-1]) {
						System.err.println("fixPeriodicEndpoints_collectCharacterInstances> ERROR (2): spacer not found: "+readID+"["+boundaries[i-1]+".."+boundaries[i]+"]");
						System.exit(1);
					}
					if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
						if (lastLeft!=-1) {
							leftCharacters_setOpen(false,false);
							leftCharacters_clear(used,bw,QUANTUM);
							tmpBoolean[i-1]=true;
						}
					}
					else {
						length=spacers[spacersCursor].breakpoint-spacers[spacersCursor].first;
						if (lastLeft!=-1) {
							leftCharacters_addLength(length);
							leftCharacters_setOpen(false,false);
						}
						else leftCharacters_load(false,i-1,length,i-1==0,false);
						leftCharacters_clear(used,bw,QUANTUM);
						length=spacers[spacersCursor].last-spacers[spacersCursor].breakpoint;
						leftCharacters_load(false,i+1,length,false,i+1==nBlocks-1);
						tmpBoolean[i-1]=true; tmpBoolean[i]=true;
					}
					i+=2;
				}
				else {
					while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[i-1])) spacersCursor++;
					if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>boundaries[i-1]) {
						System.err.println("fixPeriodicEndpoints_collectCharacterInstances> ERROR (2): spacer not found: "+readID+"["+boundaries[i-1]+".."+boundaries[i]+"]");
						System.exit(1);
					}
					if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY) {
						leftCharacters_bridge(i,nBlocks,readLength,tmpArray1,tmpArray2,tmpArray3);
						tmpBoolean[i-1]=true; tmpBoolean[i]=true;
					}
					else if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
						if (lastLeft!=-1) {
							leftCharacters_setOpen(false,false);
							leftCharacters_clear(used,bw,QUANTUM);
							tmpBoolean[i-1]=true;
						}
					}
					else {
						length=spacers[spacersCursor].breakpoint-spacers[spacersCursor].first;
						if (lastLeft!=-1) leftCharacters_addLength(length);
						else leftCharacters_load(false,i-1,length,i-1==0,false);
						leftCharacters_clear(used,bw,QUANTUM);
						length=spacers[spacersCursor].last-spacers[spacersCursor].breakpoint;
						leftCharacters_load(false,i+1,length,false,i+1==nBlocks-1);
						tmpBoolean[i-1]=true; tmpBoolean[i]=true;
					}
					i+=2;
				}
			}
			else if (isBlockPeriodic[i-1]) {
				if (isBlockNonperiodic[i-1]) {
					if (lastLeft!=-1) {
						leftCharacters_setOpen(false,false);
						leftCharacters_clear(used,bw,QUANTUM);
						tmpBoolean[i-1]=true;
					}
					i+=2;
				}
				else {
					while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[i-1])) spacersCursor++;
					if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>boundaries[i-1]) {
						System.err.println("fixPeriodicEndpoints_collectCharacterInstances> ERROR (3): spacer not found: "+readID+"["+boundaries[i-1]+".."+boundaries[i]+"]");
						System.exit(1);
					}
					if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
						if (lastLeft!=-1) {
							leftCharacters_setOpen(false,false);
							leftCharacters_clear(used,bw,QUANTUM);
							tmpBoolean[i-1]=true;
						}
					}
					else {
						length=boundaries[i]-boundaries[i-1];
						if (lastLeft==-1) leftCharacters_load(false,i-1,length,i-1==0,false);
						else {
							leftCharacters_addLength(length);
							leftCharacters_setOpen(false,false);
						}
						leftCharacters_clear(used,bw,QUANTUM);
						tmpBoolean[i-1]=true; tmpBoolean[i]=true;
					}
					i+=2;
				}
			}
			else if (isBlockPeriodic[i+1]) {
				if (isBlockNonperiodic[i+1]) i+=2;
				else {
					while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[i-1])) spacersCursor++;
					if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>boundaries[i-1]) {
						System.err.println("fixPeriodicEndpoints_collectCharacterInstances> ERROR (4): spacer not found: "+readID+"["+boundaries[i-1]+".."+boundaries[i]+"]");
						System.exit(1);
					}
					if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) { /* NOP */ }
					else {
						length=boundaries[i]-boundaries[i-1];
						leftCharacters_load(false,i+1,length,false,i+1==nBlocks-1);
						tmpBoolean[i]=true; tmpBoolean[i+1]=true;
					}
					i+=2;
				}
			}
			else i+=2;
		}
		
		// Last block. Remark: if the block is unique and shrinks, it does not create a
		// new unique character in the alphabet.
		if (i==nBlocks-1) {
			if (isBlockUnique[nBlocks-1] && isBlockPeriodic[nBlocks-2] && !isBlockNonperiodic[nBlocks-2] && readLength-boundaries[nBlocks-2]<=maxSpacerLength) {
				while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[nBlocks-2])) spacersCursor++;
				if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>boundaries[nBlocks-2]) {
					System.err.println("fixPeriodicEndpoints_collectCharacterInstances> ERROR (3): spacer not found: "+readID+"["+boundaries[nBlocks-2]+".."+(readLength-1)+"]");
					System.exit(1);
				}
				if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY) {
					length=readLength-boundaries[nBlocks-2];
					if (lastLeft==-1) leftCharacters_load(false,nBlocks-2,length,false,true);
					else leftCharacters_addLength(length);
					leftCharacters_clear(used,bw,QUANTUM);
					tmpBoolean[nBlocks-2]=true; tmpBoolean[nBlocks-1]=true;
				}
				else if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
					if (lastLeft!=-1) {
						leftCharacters_setOpen(false,false);
						leftCharacters_clear(used,bw,QUANTUM);
						tmpBoolean[nBlocks-2]=true;
					}
				}
				else {
					length=spacers[spacersCursor].breakpoint<readLength-QUANTUM?spacers[spacersCursor].breakpoint-spacers[spacersCursor].first:readLength-boundaries[nBlocks-2];
					if (lastLeft==-1) leftCharacters_load(false,nBlocks-2,length,false,spacers[spacersCursor].breakpoint>=readLength-QUANTUM);
					else leftCharacters_addLength(length);
					leftCharacters_clear(used,bw,QUANTUM);
					tmpBoolean[nBlocks-2]=true; 
					if (spacers[spacersCursor].breakpoint>=readLength-QUANTUM) tmpBoolean[nBlocks-1]=true;
				}
			}
			else {
				if (lastLeft!=-1) {
					leftCharacters_setOpen(false,false);
					leftCharacters_clear(used,bw,QUANTUM);
					tmpBoolean[nBlocks-2]=true;
				}
			}
		}
		else {
			if (lastLeft!=-1) {
				leftCharacters_clear(used,bw,QUANTUM);
				tmpBoolean[nBlocks-1]=true;
			}
		}
		
		// Marking as used every character in every other block
		for (i=0; i<nBlocks; i++) {
			if (tmpBoolean[i]) continue;
			for (j=0; j<=lastInBlock_int[i]; j++) used[intBlocks[i][j]]=true;
		}
		
		return spacersCursor;
	}
	
	
	/**
	 * Adds to $leftCharacters$ every character in block $intBlockID$. If the character is
	 * periodic, resets its open status and its length.
	 * Remark: lengths are not quantized.
	 *
	 * @param lengthMode TRUE=resets all lengths to $length$; FALSE: adds $length$ to 
	 * every existing length.
	 */
	private static final void leftCharacters_load(boolean lengthMode, int intBlockID, int length, boolean openStart, boolean openEnd) {
		boolean isPeriodic;
		int i, j, c;
		final int last = lastInBlock_int[intBlockID];
		
		for (i=0; i<=last; i++) {
			c=intBlocks[intBlockID][i];
			isPeriodic=c>lastUnique&&c<=lastPeriodic;
			lastLeft++;
			if (lastLeft==leftCharacters.length) {
				Character[] newArray = new Character[leftCharacters.length<<1];
				System.arraycopy(leftCharacters,0,newArray,0,leftCharacters.length);
				for (j=leftCharacters.length; j<newArray.length; j++) newArray[j] = new Character();
				leftCharacters=newArray;
			}
			leftCharacters[lastLeft].copyFrom(alphabet[c]);
			if (isPeriodic) {
				if (lengthMode) leftCharacters[lastLeft].length=length;
				else leftCharacters[lastLeft].length+=length;
				if (leftCharacters[lastLeft].orientation) {
					leftCharacters[lastLeft].openStart=openStart;
					leftCharacters[lastLeft].openEnd=openEnd;
				}
				else {
					leftCharacters[lastLeft].openStart=openEnd;
					leftCharacters[lastLeft].openEnd=openStart;
				}
			}
		}
	}
	
	
	/**
	 * Adds $length$ to every periodic character.
	 */
	private static final void leftCharacters_addLength(int length) {
		for (int i=0; i<=lastLeft; i++) {
			if (leftCharacters[i].repeat!=UNIQUE && leftCharacters[i].start==-1) leftCharacters[i].length+=length;
		}
	}
	
	
	/**
	 * @param openStatus TRUE=open, FALSE=closed;
	 * @param startOrEnd sets to $openStatus$ this end of every periodic character
	 * (TRUE=start, FALSE=end).
	 */
	private static final void leftCharacters_setOpen(boolean startOrEnd, boolean openStatus) {
		for (int i=0; i<=lastLeft; i++) {
			if (leftCharacters[i].repeat!=UNIQUE && leftCharacters[i].start==-1) {
				if (leftCharacters[i].orientation) {
					if (startOrEnd) leftCharacters[i].openStart=openStatus;
					else leftCharacters[i].openEnd=openStatus;
				}
				else {
					if (startOrEnd) leftCharacters[i].openEnd=openStatus;
					else leftCharacters[i].openStart=openStatus;
				}
			}
		}
	}
	
	
	/**
	 * Empties $leftCharacters$, marking in $used$ every element already in $alphabet$, 
	 * and writing to $bw$ every element not in $alphabet$.
	 */
	private static final void leftCharacters_clear(boolean[] used, BufferedWriter bw, int quantum) throws IOException {
		int i, p;
		
		for (i=0; i<=lastLeft; i++) {
			if (leftCharacters[i].repeat==UNIQUE) {
				p=Arrays.binarySearch(alphabet,0,lastUnique+1,leftCharacters[i]);
				if (p>=0) used[p]=true;
				else {
					System.err.println("leftCharacters_clear> ERROR: new unique character not in the alphabet: "+leftCharacters[i]);
					System.exit(1);
				}
			}
			else if (leftCharacters[i].start==-1) {
				leftCharacters[i].quantize(quantum);
				p=Arrays.binarySearch(alphabet,lastUnique+1,lastPeriodic+1,leftCharacters[i]);
				if (p>=0) used[p]=true;
				else bw.write(leftCharacters[i].toString()+"\n");
			}
			else {
				p=Arrays.binarySearch(alphabet,lastPeriodic+1,lastAlphabet+1,leftCharacters[i]);
				if (p>=0) used[p]=true;
				else {
					System.err.println("leftCharacters_clear> ERROR: new non-periodic character not in the alphabet: "+leftCharacters[i]);
					System.err.println("lastLeft="+lastLeft+", blockCursor="+blockCursor+", currentBoundary="+currentBoundary+", nBoundariesWritten="+nBoundariesWritten);
					System.err.println("leftCharacters:");
					for (int x=0; x<=lastLeft; x++) System.err.println(leftCharacters[x]);
					System.exit(1);
				}
			}
		}
		lastLeft=-1;
	}
	
	
	/**
	 * Assume that $intBlockID$ is a unique block, and that both $intBlockID-1$ and 
	 * $intBlockID+1$ are periodic. If $leftCharacters$ is not empty, the procedure 
	 * extends every periodic character in it up to the end of $intBlockID+1$, and adds
	 * an instance of every periodic repeat that occurs in $intBlockID+1$ but does not
	 * occur in $intBlockID-1$. If $leftCharacters$ is empty, it is initialized in the
	 * same way (only using periodic repeats).
	 *
	 * @param tmpArray* temporary space, of size at least equal to the max number of 
	 * characters in a block.
	 */
	private static final void leftCharacters_bridge(int intBlockID, int nBlocks, int readLength, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) {
		final boolean openStart = intBlockID-1==0;
		final boolean openEnd = intBlockID+1==nBlocks-1;
		final int length = (intBlockID+1==nBlocks-1?readLength:boundaries[intBlockID+1])-boundaries[intBlockID-1];
		boolean found, foundOpen, orientation;
		int i, j, c, p;
		int last, last1, last2, last3, foundLength, repeat;
		Character character;
		
		j=-1; last=lastInBlock_int[intBlockID-1];
		for (i=0; i<=last; i++) {
			c=intBlocks[intBlockID-1][i];
			if (c<=lastUnique || c>lastPeriodic) continue;
			character=alphabet[c];
			tmpArray1[++j]=character.orientation?character.repeat:-1-character.repeat;
		}
		last1=j;
		if (last1>0) Arrays.sort(tmpArray1,0,last1+1);
		j=0;
		for (i=1; i<=last1; i++) {
			if (tmpArray1[i]==tmpArray1[j]) continue;
			tmpArray1[++j]=tmpArray1[i];
		}
		last1=j;
		j=-1; last=lastInBlock_int[intBlockID+1];
		for (i=0; i<=last; i++) {
			c=intBlocks[intBlockID+1][i];
			if (c<=lastUnique || c>lastPeriodic) continue;
			character=alphabet[c];
			tmpArray2[++j]=character.orientation?character.repeat:-1-character.repeat;
		}
		last2=j;
		if (last2>0) Arrays.sort(tmpArray2,0,last2+1);
		j=0;
		for (i=1; i<=last2; i++) {
			if (tmpArray2[i]==tmpArray2[j]) continue;
			tmpArray2[++j]=tmpArray2[i];
		}
		last2=j;
		last3=Math.setUnion(tmpArray1,last1,tmpArray2,last2,tmpArray3);
		if (lastLeft==-1) {
			for (i=0; i<=last3; i++) {
				lastLeft++;
				if (lastLeft==leftCharacters.length) {
					Character[] newArray = new Character[leftCharacters.length<<1];
					System.arraycopy(leftCharacters,0,newArray,0,leftCharacters.length);
					for (j=leftCharacters.length; j<newArray.length; j++) newArray[j] = new Character();
					leftCharacters=newArray;
				}
				leftCharacters[lastLeft].start=-1; leftCharacters[lastLeft].end=-1;
				if (tmpArray3[i]<0) {
					leftCharacters[lastLeft].orientation=false;
					leftCharacters[lastLeft].repeat=-1-tmpArray3[i];
				}
				else {
					leftCharacters[lastLeft].orientation=true;
					leftCharacters[lastLeft].repeat=tmpArray3[i];
				}
				leftCharacters[lastLeft].length=(intBlockID-1==0?boundaries[intBlockID-1]:boundaries[intBlockID-1]-boundaries[intBlockID-2])+length;
				if (leftCharacters[lastLeft].orientation) {
					leftCharacters[lastLeft].openStart=openStart;
					leftCharacters[lastLeft].openEnd=openEnd;
				}
				else {
					leftCharacters[lastLeft].openStart=openEnd;
					leftCharacters[lastLeft].openEnd=openStart;
				}
			}
		}
		else {
			leftCharacters_addLength(length);
			foundOpen=false; foundLength=0;
			for (i=0; i<=lastLeft; i++) {
				if (leftCharacters[i].repeat!=UNIQUE && leftCharacters[i].start==-1) {
					if (leftCharacters[i].openStart || leftCharacters[i].openEnd) foundOpen=true;
					foundLength=Math.max(foundLength,leftCharacters[i].length);
				}
			}
			leftCharacters_setOpen(false,openEnd);
			p=lastLeft;
			for (i=0; i<=last3; i++) {
				if (tmpArray3[i]>=0) { repeat=tmpArray3[i]; orientation=true; }
				else { repeat=-1-tmpArray3[i]; orientation=false; }
				found=false;
				for (j=0; j<=p; j++) {
					if (leftCharacters[j].repeat==repeat && leftCharacters[j].orientation==orientation) {
						found=true;
						break;
					}
				}
				if (!found) {
					lastLeft++;
					if (lastLeft==leftCharacters.length) {
						Character[] newArray = new Character[leftCharacters.length<<1];
						System.arraycopy(leftCharacters,0,newArray,0,leftCharacters.length);
						for (j=leftCharacters.length; j<newArray.length; j++) newArray[j] = new Character();
						leftCharacters=newArray;
					}
					leftCharacters[lastLeft].repeat=repeat;
					leftCharacters[lastLeft].start=-1; leftCharacters[lastLeft].end=-1;
					leftCharacters[lastLeft].length=foundLength;
					if (orientation) {
						leftCharacters[lastLeft].orientation=true;
						leftCharacters[lastLeft].openStart=foundOpen;
						leftCharacters[lastLeft].openEnd=openEnd;
					}
					else {
						leftCharacters[lastLeft].orientation=false;
						leftCharacters[lastLeft].openStart=openEnd;
						leftCharacters[lastLeft].openEnd=foundOpen;
					}
				}
			}
		}
	}
	
	
	/**
	 * Like $fixPeriodicEndpoints_collectCharacterInstances()$, but looks up in the new
	 * alphabet every existing and new character induced by fixing spacers.
	 *
	 * @param out the procedure cumulates the number of spacers fixed at every length 
	 * (multiple of $IO.quantum$);
	 * @param tmpArray* temporary space, with a number of cells at least equal to the max
	 * number of elements in a block of the translation.
	 */
	public static final int fixPeriodicEndpoints_updateTranslation(int readID, int readLength, int spacersCursor, int maxSpacerLength, String read2characters_old, String read2boundaries_old, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, BufferedWriter fullyContained_new, int[] out, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) throws IOException {
		final int CAPACITY = 10;  // Arbitrary
		final int QUANTUM = IO.quantum;
		boolean isUnique, isBlockPeriodicLeft, isBlockPeriodicRight, isBlockNonperiodicLeft, isBlockNonperiodicRight;
		int i, c, p;
		int nBlocks, length;
	
		if (read2characters_old.length()==0) {
			read2characters_new.newLine(); read2boundaries_new.newLine();
			return spacersCursor;
		}
		if (tmpCharacter==null) tmpCharacter = new Character();
		if (leftCharacters==null) {
			leftCharacters = new Character[CAPACITY];
			for (i=0; i<leftCharacters.length; i++) leftCharacters[i] = new Character();
		}
		loadBoundaries(read2boundaries_old);
		nBlocks=loadBlocks(read2characters_old);
		
		// First block
		spacersCursor=fixPeriodicEndpoints_updateTranslation_firstBlock(nBlocks,readID,readLength,maxSpacerLength,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,spacersCursor,QUANTUM,read2characters_new,read2boundaries_new,out,tmpArray1);
		if (nBlocks==1) return spacersCursor;
		else if (nBlocks==2) return fixPeriodicEndpoints_updateTranslation_secondBlock(readID,readLength,spacersCursor,maxSpacerLength,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,QUANTUM,read2characters_new,read2boundaries_new,fullyContained_new,out,tmpArray1);
		
		// Intermediate blocks
		while (blockCursor<=nBlocks-2) {
			p=Integer.parseInt(blocks[blockCursor][0]);
			if (p<0) p=-1-p;
			isUnique=p<=lastUnique_old||p==lastAlphabet_old+1;
			isBlockPeriodicLeft=false; isBlockNonperiodicLeft=false;
			for (i=0; i<=lastInBlock[blockCursor-1]; i++) {
				c=Integer.parseInt(blocks[blockCursor-1][i]);
				if (c<0) c=-1-c;
				if (c>lastUnique_old && c<=lastPeriodic_old) isBlockPeriodicLeft=true;
				else isBlockNonperiodicLeft=true;
			}
			isBlockPeriodicRight=false; isBlockNonperiodicRight=false;
			for (i=0; i<=lastInBlock[blockCursor+1]; i++) {
				c=Integer.parseInt(blocks[blockCursor+1][i]);
				if (c<0) c=-1-c;
				if (c>lastUnique_old && c<=lastPeriodic_old) isBlockPeriodicRight=true;
				else isBlockNonperiodicRight=true;
			}
			if (!isUnique || boundaries[blockCursor]-boundaries[blockCursor-1]>maxSpacerLength) fixPeriodicEndpoints_updateTranslation_noSpacer((isUnique||isBlockPeriodicRight)?2:1,readID,readLength,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray1);
			else if (isBlockPeriodicLeft && isBlockPeriodicRight) {
				if (isBlockNonperiodicLeft && isBlockNonperiodicRight) fixPeriodicEndpoints_updateTranslation_noSpacer(2,readID,readLength,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray1);
				else spacersCursor=fixPeriodicEndpoints_updateTranslation_spacer(readID,readLength,spacersCursor,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,tmpArray1,tmpArray2,tmpArray3,out);
			}
			else if (isBlockPeriodicLeft && isBlockNonperiodicRight) {
				if (isBlockNonperiodicLeft) fixPeriodicEndpoints_updateTranslation_noSpacer(1,readID,readLength,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray1);
				else {
					while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[blockCursor-1])) spacersCursor++;
					if (spacersCursor>lastSpacer || spacers[spacersCursor].first>boundaries[blockCursor-1]) {
						System.err.println("fixPeriodicEndpoints_updateTranslation> ERROR: spacer not found: "+readID+"["+boundaries[blockCursor-1]+".."+boundaries[blockCursor]+"]");
						System.exit(1);
					}
					if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
						fixPeriodicEndpoints_updateTranslation_noSpacer(2,readID,readLength,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray1);
					}
					else {
						if (currentBoundary!=-1) {
							read2characters_new.write(SEPARATOR_MAJOR+"");
							read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+currentBoundary);
							nBoundariesWritten++;
						}
						length=blockCursor-1==0?boundaries[blockCursor-1]:boundaries[blockCursor-1]-boundaries[blockCursor-2];
						if (lastLeft==-1) leftCharacters_load_prime(true,blockCursor-1,length,blockCursor-1==0,false,oldAlphabet,lastUnique_old,lastPeriodic_old);
						leftCharacters_addLength(boundaries[blockCursor]-boundaries[blockCursor-1]);
						leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,QUANTUM,tmpArray1);
						currentBoundary=boundaries[blockCursor];
						blockCursor+=2;
					}
				}
			}
			else if (isBlockPeriodicRight && isBlockNonperiodicLeft) {
				if (isBlockNonperiodicRight) fixPeriodicEndpoints_updateTranslation_noSpacer(2,readID,readLength,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray1);
				else {
					while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[blockCursor-1])) spacersCursor++;
					if (spacersCursor>lastSpacer || spacers[spacersCursor].first>boundaries[blockCursor-1]) {
						System.err.println("fixPeriodicEndpoints_updateTranslation> ERROR: spacer not found: "+readID+"["+boundaries[blockCursor-1]+".."+boundaries[blockCursor]+"]");
						System.exit(1);
					}
					if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) fixPeriodicEndpoints_updateTranslation_noSpacer(2,readID,readLength,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray1);
					else {
						if (currentBoundary!=-1) {
							read2characters_new.write(SEPARATOR_MAJOR+"");
							read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+currentBoundary);
							nBoundariesWritten++;
						}
						length=blockCursor-1==0?boundaries[blockCursor-1]:boundaries[blockCursor-1]-boundaries[blockCursor-2];
						writeBlock(blockCursor-1,length,blockCursor-1==0,false,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray1);
						currentBoundary=boundaries[blockCursor-1];
						length=(blockCursor+1==nBlocks-1?readLength:boundaries[blockCursor+1])-boundaries[blockCursor];
						leftCharacters_load_prime(true,blockCursor+1,length,false,blockCursor+1==nBlocks-1,oldAlphabet,lastUnique_old,lastPeriodic_old);
						leftCharacters_addLength(boundaries[blockCursor]-boundaries[blockCursor-1]);
						blockCursor+=2;
					}
				}
			}
			else fixPeriodicEndpoints_updateTranslation_noSpacer(2,readID,readLength,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray1);
		}
		
		// Last block
		if (blockCursor==nBlocks-1) spacersCursor=fixPeriodicEndpoints_updateTranslation_lastBlock(readID,readLength,nBlocks,spacersCursor,maxSpacerLength,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray1);
		else if (blockCursor==nBlocks) fixPeriodicEndpoints_updateTranslation_afterLastBlock(readLength,nBlocks,QUANTUM,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,tmpArray1);
		else { /* NOP: this case is already fully handled before. */ }
		read2characters_new.newLine(); read2boundaries_new.newLine();
		if (nBoundariesWritten==0) fullyContained_new.write(readID+"\n");
		
		return spacersCursor;
	}
	
	
	/**
	 * @return the new value of $spacersCursor$ after the procedure completes.
	 */
	private static final int fixPeriodicEndpoints_updateTranslation_spacer(int readID, int readLength, int spacersCursor, int nBlocks, int quantum, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3, int[] out) throws IOException {
		int length;
		
		while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[blockCursor-1])) spacersCursor++;
		if (spacersCursor>lastSpacer || spacers[spacersCursor].first>boundaries[blockCursor-1]) {
			System.err.println("fixPeriodicEndpoints_updateTranslation_spacer> ERROR: spacer not found: "+readID+"["+boundaries[blockCursor-1]+".."+boundaries[blockCursor]+"]");
			System.exit(1);
		}
		if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY) {
			fixPeriodicEndpoints_updateTranslation_bridge(readLength,nBlocks,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,tmpArray1,tmpArray2,tmpArray3);
			return spacersCursor;
		}
		else if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
			fixPeriodicEndpoints_updateTranslation_noSpacer(2,readID,readLength,nBlocks,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,out,tmpArray3);
			return spacersCursor;
		}
		if (currentBoundary!=-1) {
			read2characters_new.write(SEPARATOR_MAJOR+"");
			read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+currentBoundary);
			nBoundariesWritten++;
		}
		out[Math.min((boundaries[blockCursor]-boundaries[blockCursor-1])/IO.quantum,out.length-1)]++;
		length=blockCursor-1==0?boundaries[blockCursor-1]:boundaries[blockCursor-1]-boundaries[blockCursor-2];
		if (lastLeft==-1) leftCharacters_load_prime(true,blockCursor-1,length,blockCursor-1==0,false,oldAlphabet,lastUnique_old,lastPeriodic_old);
		leftCharacters_addLength(spacers[spacersCursor].breakpoint-spacers[spacersCursor].first);
		leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray1);
		length=boundaries[blockCursor+1]-boundaries[blockCursor];
		leftCharacters_load_prime(true,blockCursor+1,length,false,blockCursor+1==nBlocks-1,oldAlphabet,lastUnique_old,lastPeriodic_old);
		leftCharacters_addLength(spacers[spacersCursor].last-spacers[spacersCursor].breakpoint);
		leftCharacters_setOpen(true,false);
		currentBoundary=spacers[spacersCursor].breakpoint;
		blockCursor+=2;
		return spacersCursor;
	}
	
	
	/**
	 * @param shift advance $blockCursor$ by {1,2,3}.
	 */
	private static final void fixPeriodicEndpoints_updateTranslation_noSpacer(int shift, int readID, int readLength, int nBlocks, int quantum, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, int[] out, int[] tmpArray) throws IOException {
		int length;

		if (currentBoundary!=-1) {
			read2characters_new.write(SEPARATOR_MAJOR+"");
			read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+currentBoundary);
			nBoundariesWritten++;
		}
		if (lastLeft!=-1) leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
		else writeBlock(blockCursor-1,blockCursor-1==0?boundaries[blockCursor-1]:boundaries[blockCursor-1]-boundaries[blockCursor-2],blockCursor-1==0,false,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
		if (shift==1) currentBoundary=boundaries[blockCursor-1];
		else {
			read2characters_new.write(SEPARATOR_MAJOR+"");
			read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+boundaries[blockCursor-1]);
			nBoundariesWritten++;
			writeBlock(blockCursor,boundaries[blockCursor]-boundaries[blockCursor-1],false,false,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
			if (shift==2) currentBoundary=boundaries[blockCursor];
			else {
				read2characters_new.write(SEPARATOR_MAJOR+"");
				read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+boundaries[blockCursor]);
				nBoundariesWritten++;
				length=(blockCursor+1==nBlocks-1?readLength:boundaries[blockCursor+1])-boundaries[blockCursor];
				writeBlock(blockCursor+1,length,false,blockCursor+1==nBlocks-1,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
				if (blockCursor+1<nBlocks-1) currentBoundary=boundaries[blockCursor+1];
			}
		}
		blockCursor+=shift;
	}
	
	
	/**
	 * Identical to $leftCharacters_bridge()$.
	 */
	private static final void fixPeriodicEndpoints_updateTranslation_bridge(int readLength, int nBlocks, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3) {
		final boolean openStart = blockCursor-1==0;
		final boolean openEnd = blockCursor+1==nBlocks-1;
		final int length = (blockCursor+1==nBlocks-1?readLength:boundaries[blockCursor+1])-boundaries[blockCursor-1];
		boolean found, foundOpen, orientation;
		int i, j, c, p;
		int last, last1, last2, last3, foundLength, repeat;
		Character character;
		
		j=-1; last=lastInBlock[blockCursor-1];
		for (i=0; i<=last; i++) {
			c=Integer.parseInt(blocks[blockCursor-1][i]);
			if (c<0) c=-1-c;
			if (c<=lastUnique_old || c>lastPeriodic_old) continue;
			character=oldAlphabet[c];
			tmpArray1[++j]=character.orientation?character.repeat:-1-character.repeat;
		}
		last1=j;
		if (last1>0) Arrays.sort(tmpArray1,0,last1+1);
		j=0;
		for (i=1; i<=last1; i++) {
			if (tmpArray1[i]==tmpArray1[j]) continue;
			tmpArray1[++j]=tmpArray1[i];
		}
		last1=j;
		j=-1; last=lastInBlock[blockCursor+1];
		for (i=0; i<=last; i++) {
			c=Integer.parseInt(blocks[blockCursor+1][i]);
			if (c<0) c=-1-c;
			if (c<=lastUnique_old || c>lastPeriodic_old) continue;
			character=oldAlphabet[c];
			tmpArray2[++j]=character.orientation?character.repeat:-1-character.repeat;
		}
		last2=j;
		if (last2>0) Arrays.sort(tmpArray2,0,last2+1);
		j=0;
		for (i=1; i<=last2; i++) {
			if (tmpArray2[i]==tmpArray2[j]) continue;
			tmpArray2[++j]=tmpArray2[i];
		}
		last2=j;
		last3=Math.setUnion(tmpArray1,last1,tmpArray2,last2,tmpArray3);
		if (lastLeft==-1) {
			for (i=0; i<=last3; i++) {
				lastLeft++;
				if (lastLeft==leftCharacters.length) {
					Character[] newArray = new Character[leftCharacters.length<<1];
					System.arraycopy(leftCharacters,0,newArray,0,leftCharacters.length);
					for (j=leftCharacters.length; j<newArray.length; j++) newArray[j] = new Character();
					leftCharacters=newArray;
				}
				leftCharacters[lastLeft].start=-1; leftCharacters[lastLeft].end=-1;
				if (tmpArray3[i]<0) {
					leftCharacters[lastLeft].orientation=false;
					leftCharacters[lastLeft].repeat=-1-tmpArray3[i];
				}
				else {
					leftCharacters[lastLeft].orientation=true;
					leftCharacters[lastLeft].repeat=tmpArray3[i];
				}
				leftCharacters[lastLeft].length=(blockCursor-1==0?boundaries[blockCursor-1]:boundaries[blockCursor-1]-boundaries[blockCursor-2])+length;
				if (leftCharacters[lastLeft].orientation) {
					leftCharacters[lastLeft].openStart=openStart;
					leftCharacters[lastLeft].openEnd=openEnd;
				}
				else {
					leftCharacters[lastLeft].openStart=openEnd;
					leftCharacters[lastLeft].openEnd=openStart;
				}
			}
		}
		else {
			leftCharacters_addLength(length);
			foundOpen=false; foundLength=0;
			for (i=0; i<=lastLeft; i++) {
				if (leftCharacters[i].repeat!=UNIQUE && leftCharacters[i].start==-1) {
					if (leftCharacters[i].openStart || leftCharacters[i].openEnd) foundOpen=true;
					foundLength=Math.max(foundLength,leftCharacters[i].length);
				}
			}
			leftCharacters_setOpen(false,openEnd);
			p=lastLeft;
			for (i=0; i<=last3; i++) {
				if (tmpArray3[i]>=0) { repeat=tmpArray3[i]; orientation=true; }
				else { repeat=-1-tmpArray3[i]; orientation=false; }
				found=false;
				for (j=0; j<=p; j++) {
					if (leftCharacters[j].repeat==repeat && leftCharacters[j].orientation==orientation) {
						found=true;
						break;
					}
				}
				if (!found) {
					lastLeft++;
					if (lastLeft==leftCharacters.length) {
						Character[] newArray = new Character[leftCharacters.length<<1];
						System.arraycopy(leftCharacters,0,newArray,0,leftCharacters.length);
						for (j=leftCharacters.length; j<newArray.length; j++) newArray[j] = new Character();
						leftCharacters=newArray;
					}
					leftCharacters[lastLeft].repeat=repeat;
					leftCharacters[lastLeft].start=-1; leftCharacters[lastLeft].end=-1;
					leftCharacters[lastLeft].length=foundLength;
					if (orientation) {
						leftCharacters[lastLeft].orientation=true;
						leftCharacters[lastLeft].openStart=foundOpen;
						leftCharacters[lastLeft].openEnd=openEnd;
					}
					else {
						leftCharacters[lastLeft].orientation=false;
						leftCharacters[lastLeft].openStart=openEnd;
						leftCharacters[lastLeft].openEnd=foundOpen;
					}
				}
			}			
		}
		blockCursor+=2;
	}
	
	
	/**
	 * @return the new value of $spacersCursor$.
	 */
	private static final int fixPeriodicEndpoints_updateTranslation_firstBlock(int nBlocks, int readID, int readLength, int maxSpacerLength, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, int spacersCursor, int quantum, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, int[] out, int[] tmpArray) throws IOException {
		boolean isUnique, foundPeriodic, foundNonperiodic;
		int i, c, p;
		int length, last;
		
		nBoundariesWritten=0;
		
		// Single block
		if (nBlocks==1) {
			writeBlock(0,readLength,true,true,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
			read2characters_new.newLine(); read2boundaries_new.newLine();
			return spacersCursor;
		}
		
		// Multiple blocks
		p=Integer.parseInt(blocks[0][0]);
		if (p<0) p=-1-p;
		isUnique=p<=lastUnique_old||p==lastAlphabet_old+1;
		foundPeriodic=false; foundNonperiodic=false;
		for (i=0; i<=lastInBlock[1]; i++) {
			c=Integer.parseInt(blocks[1][i]);
			if (c<0) c=-1-c;
			if (c>lastUnique_old && c<=lastPeriodic_old) foundPeriodic=true;
			else foundNonperiodic=true;
		}
		if (isUnique && boundaries[0]<=maxSpacerLength && foundPeriodic && !foundNonperiodic) {
			// The first block is a spacer
			while (spacersCursor<=lastSpacer && spacers[spacersCursor].read<readID) spacersCursor++;
			if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>0) {
				System.err.println("fixPeriodicEndpoints_updateTranslation> ERROR (1): spacer not found: "+readID+"[0.."+boundaries[0]+"]");
				System.exit(1);
			}
			if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY) {
				currentBoundary=-1; length=boundaries[0];
				leftCharacters_load_prime(true,1,(nBlocks==2?readLength:boundaries[1])-boundaries[0],true,nBlocks==2,oldAlphabet,lastUnique_old,lastPeriodic_old);
				leftCharacters_addLength(length);
			}
			else if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
				// The first block is an inactive spacer
				writeBlock(0,boundaries[0],true,false,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
				currentBoundary=boundaries[0];
				lastLeft=-1;
			}
			else {
				if (spacers[spacersCursor].breakpoint>quantum) {
					tmpCharacter.repeat=UNIQUE; tmpCharacter.orientation=true;
					tmpCharacter.start=-1; tmpCharacter.end=-1;
					tmpCharacter.openStart=true; tmpCharacter.openEnd=false;
					tmpCharacter.length=spacers[spacersCursor].breakpoint;
					tmpCharacter.quantize(quantum);
					p=fixPeriodicEndpoints_lookupUnique(tmpCharacter,newAlphabet,lastUnique_new,lastAlphabet_new);
					read2characters_new.write(p+"");
					currentBoundary=spacers[spacersCursor].breakpoint;
					length=spacers[spacersCursor].last-spacers[spacersCursor].breakpoint;
				}
				else { currentBoundary=-1; length=boundaries[0]; }
				leftCharacters_load_prime(true,1,(nBlocks==2?readLength:boundaries[1])-boundaries[0],spacers[spacersCursor].breakpoint<=quantum,nBlocks==2,oldAlphabet,lastUnique_old,lastPeriodic_old);
				leftCharacters_addLength(length);
			}
			out[Math.min(boundaries[0]/IO.quantum,out.length-1)]++;
			blockCursor=2;
		}
		else { 
			// The first block is not a spacer
			blockCursor=1; lastLeft=-1; currentBoundary=-1;
		}
		return spacersCursor;
	}
	
	
	/**
	 * Handles the second block when the read contains exactly two blocks.
	 *
 	 * @return the new value of $spacersCursor$.
	 */
	private static final int fixPeriodicEndpoints_updateTranslation_secondBlock(int readID, int readLength, int spacersCursor, int maxSpacerLength, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, int quantum, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, BufferedWriter fullyContained_new, int[] out, int[] tmpArray) throws IOException {
		boolean isUnique, foundPeriodic, foundNonperiodic;
		int i, c, p;

		p=Integer.parseInt(blocks[1][0]);
		if (p<0) p=-1-p;
		isUnique=p<=lastUnique_old||p==lastAlphabet_old+1;
		foundPeriodic=false; foundNonperiodic=false;
		for (i=0; i<=lastInBlock[0]; i++) {
			c=Integer.parseInt(blocks[0][i]);
			if (c<0) c=-1-c;
			if (c>lastUnique_old && c<=lastPeriodic_old) foundPeriodic=true;
			else foundNonperiodic=true;
		}
		if (isUnique && readLength-boundaries[0]<=maxSpacerLength && foundPeriodic && !foundNonperiodic) {
			// The second block is a spacer
			while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[0])) spacersCursor++;
			if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>boundaries[0]) {
				System.err.println("fixPeriodicEndpoints_updateTranslation> ERROR (2): spacer not found: "+readID+"["+boundaries[0]+".."+(readLength-1)+"]");
				System.exit(1);
			}
			if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY) {
				writeBlock(0,readLength,true,true,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
				fullyContained_new.write(readID+"\n");
			}
			else if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
				// The second block is an inactive spacer
				writeBlock(0,boundaries[0],true,false,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
				read2characters_new.write(SEPARATOR_MAJOR+"");
				read2boundaries_new.write(boundaries[0]+"");
				nBoundariesWritten++;
				writeBlock(1,readLength-boundaries[0],false,true,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
			}
			else {
				if (readLength-spacers[spacersCursor].breakpoint<quantum) {
					writeBlock(0,readLength,true,true,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
					fullyContained_new.write(readID+"\n");
				}
				else {
					leftCharacters_load_prime(false,0,spacers[spacersCursor].breakpoint-spacers[spacersCursor].first,true,false,oldAlphabet,lastUnique_old,lastPeriodic_old);
					leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
					read2characters_new.write(SEPARATOR_MAJOR+"");
					read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+spacers[spacersCursor].breakpoint);
					nBoundariesWritten++;
					tmpCharacter.repeat=UNIQUE; tmpCharacter.orientation=true;
					tmpCharacter.start=-1; tmpCharacter.end=-1;
					tmpCharacter.openStart=false; tmpCharacter.openEnd=true;
					tmpCharacter.length=readLength-spacers[spacersCursor].breakpoint;
					tmpCharacter.quantize(quantum);
					p=fixPeriodicEndpoints_lookupUnique(tmpCharacter,newAlphabet,lastUnique_new,lastAlphabet_new);
					read2characters_new.write(p+"");
				}
			}
			out[Math.min((readLength-boundaries[0])/IO.quantum,out.length-1)]++;
		}
		else {
			// The second block is not a spacer
			if (blockCursor==1) {
				writeBlock(0,boundaries[0],true,false,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
				read2characters_new.write(SEPARATOR_MAJOR+"");
				read2boundaries_new.write(boundaries[0]+"");
				nBoundariesWritten++;
				writeBlock(1,readLength-boundaries[0],false,true,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
			}
			else {
				fixPeriodicEndpoints_updateTranslation_afterLastBlock(readLength,2,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,read2boundaries_new,tmpArray);
				fullyContained_new.write(readID+"\n");
			}
		}
		read2characters_new.newLine(); read2boundaries_new.newLine();
		return spacersCursor;
	}

	
	/**
	 * Handles the case where the block cursor is on the last block.
	 *
	 * @return the new value of $spacersCursor$.
	 */
	private static final int fixPeriodicEndpoints_updateTranslation_lastBlock(int readID, int readLength, int nBlocks, int spacersCursor, int maxSpacerLength, int quantum, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, int[] out, int[] tmpArray) throws IOException {
		boolean isUnique, foundPeriodic, foundNonperiodic;
		int i, c, p;
		int length;
		
		if (currentBoundary!=-1) {
			read2characters_new.write(SEPARATOR_MAJOR+"");
			read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+currentBoundary);
			nBoundariesWritten++;
		}
		p=Integer.parseInt(blocks[nBlocks-1][0]);
		if (p<0) p=-1-p;
		isUnique=p<=lastUnique_old||p==lastAlphabet_old+1;
		foundPeriodic=false; foundNonperiodic=false;
		for (i=0; i<=lastInBlock[nBlocks-2]; i++) {
			c=Integer.parseInt(blocks[nBlocks-2][i]);
			if (c<0) c=-1-c;
			if (c>lastUnique_old && c<=lastPeriodic_old) foundPeriodic=true;
			else foundNonperiodic=true;
		}
		if (isUnique && readLength-boundaries[nBlocks-2]<=maxSpacerLength && foundPeriodic && !foundNonperiodic) {
			// The last block is a spacer
			while (spacersCursor<=lastSpacer && (spacers[spacersCursor].read<readID || spacers[spacersCursor].first<boundaries[nBlocks-2])) spacersCursor++;
			if (spacersCursor>lastSpacer || spacers[spacersCursor].read>readID || spacers[spacersCursor].first>boundaries[nBlocks-2]) {
				System.err.println("fixPeriodicEndpoints_updateTranslation> ERROR (6): spacer not found: "+readID+"["+boundaries[nBlocks-2]+".."+(readLength-1)+"]");
				System.err.println("fixPeriodicEndpoints_updateTranslation> spacers[spacersCursor]="+spacers[spacersCursor]);
				System.exit(1);
			}
			if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY) {
				length=readLength-boundaries[nBlocks-2];
				if (lastLeft!=-1) {
					leftCharacters_addLength(length);
					leftCharacters_setOpen(false,true);
				}
				else leftCharacters_load_prime(false,nBlocks-2,length,false,true,oldAlphabet,lastUnique_old,lastPeriodic_old);
				leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
			}
			else if (spacers[spacersCursor].breakpoint==Math.POSITIVE_INFINITY-1) {
				// The last block is an inactive spacer
				if (lastLeft!=-1) leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
				else writeBlock(nBlocks-2,boundaries[nBlocks-2]-boundaries[nBlocks-3],false,false,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
				read2characters_new.write(SEPARATOR_MAJOR+"");
				read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+boundaries[nBlocks-2]);
				nBoundariesWritten++;
				writeBlock(nBlocks-1,readLength-boundaries[nBlocks-2],false,true,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
			}
			else {
				if (readLength-spacers[spacersCursor].breakpoint>quantum) {
					length=spacers[spacersCursor].breakpoint-boundaries[nBlocks-2];
					if (lastLeft!=-1) {
						leftCharacters_addLength(length);
						leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
					}
					else {
						leftCharacters_load_prime(false,nBlocks-2,length,false,false,oldAlphabet,lastUnique_old,lastPeriodic_old);
						leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
					}
					read2characters_new.write(SEPARATOR_MAJOR+"");
					read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+spacers[spacersCursor].breakpoint);
					nBoundariesWritten++;
					tmpCharacter.repeat=UNIQUE; tmpCharacter.orientation=true;
					tmpCharacter.start=-1; tmpCharacter.end=-1;
					tmpCharacter.openStart=false; tmpCharacter.openEnd=true;
					tmpCharacter.length=readLength-spacers[spacersCursor].breakpoint;
					tmpCharacter.quantize(quantum);
					p=fixPeriodicEndpoints_lookupUnique(tmpCharacter,newAlphabet,lastUnique_new,lastAlphabet_new);
					read2characters_new.write(p+"");
				}
				else {
					length=readLength-boundaries[nBlocks-2];
					if (lastLeft!=-1) {
						leftCharacters_addLength(length);
						leftCharacters_setOpen(false,true);
					}
					else leftCharacters_load_prime(false,nBlocks-2,length,false,true,oldAlphabet,lastUnique_old,lastPeriodic_old);
					leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
				}
			}
			out[Math.min((readLength-boundaries[nBlocks-2])/IO.quantum,out.length-1)]++;
		}
		else {
			// The last block is not a spacer
			if (lastLeft!=-1) leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
			else writeBlock(nBlocks-2,boundaries[nBlocks-2]-boundaries[nBlocks-3],false,false,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
			read2characters_new.write(SEPARATOR_MAJOR+"");
			read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+boundaries[nBlocks-2]);
			nBoundariesWritten++;
			writeBlock(nBlocks-1,readLength-boundaries[nBlocks-2],false,true,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
		}
		return spacersCursor;
	}
	
	
	/**
	 * Handles the case where the block cursor is immediately after the last block.
	 */
	private static final void fixPeriodicEndpoints_updateTranslation_afterLastBlock(int readLength, int nBlocks, int quantum, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, int[] tmpArray) throws IOException {
		if (currentBoundary!=-1) {
			read2characters_new.write(SEPARATOR_MAJOR+"");
			read2boundaries_new.write((nBoundariesWritten>0?SEPARATOR_MINOR+"":"")+currentBoundary);
			nBoundariesWritten++;
		}
		if (lastLeft!=-1) leftCharacters_clear_prime(newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,quantum,tmpArray);
		else writeBlock(nBlocks-1,readLength-boundaries[nBlocks-2],false,true,quantum,oldAlphabet,lastUnique_old,lastPeriodic_old,lastAlphabet_old,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,read2characters_new,tmpArray);
	}
	
	
	/**
	 * Variant of $leftCharacters_load()$.
	 */
	private static final void leftCharacters_load_prime(boolean lengthMode, int blockID, int length, boolean openStart, boolean openEnd, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old) {
		boolean isPeriodic;
		int i, j, c;
		final int last = lastInBlock[blockID];
		
		for (i=0; i<=last; i++) {
			c=Integer.parseInt(blocks[blockID][i]);
			if (c<0) c=-1-c;
			isPeriodic=c>lastUnique_old&&c<=lastPeriodic_old;
			lastLeft++;
			if (lastLeft==leftCharacters.length) {
				Character[] newArray = new Character[leftCharacters.length<<1];
				System.arraycopy(leftCharacters,0,newArray,0,leftCharacters.length);
				for (j=leftCharacters.length; j<newArray.length; j++) newArray[j] = new Character();
				leftCharacters=newArray;
			}
			leftCharacters[lastLeft].copyFrom(oldAlphabet[c]);
			if (isPeriodic) {
				if (lengthMode) leftCharacters[lastLeft].length=length;
				else leftCharacters[lastLeft].length+=length;
				if (leftCharacters[lastLeft].orientation) {
					leftCharacters[lastLeft].openStart=openStart;
					leftCharacters[lastLeft].openEnd=openEnd;
				}
				else {
					leftCharacters[lastLeft].openStart=openEnd;
					leftCharacters[lastLeft].openEnd=openStart;
				}
			}
		}
	}
	
	
	/**
	 * Variant of $leftCharacters_clear()$.
	 */
	private static final void leftCharacters_clear_prime(Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, BufferedWriter read2characters_new, int quantum, int[] tmpArray) throws IOException {
		int i, p;
		
		p=-1;
		for (i=0; i<=lastLeft; i++) {
			if (leftCharacters[i].repeat!=UNIQUE && leftCharacters[i].start==-1) leftCharacters[i].quantize(quantum);
			p=fixPeriodicEndpoints_updateTranslation_impl(leftCharacters[i],newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,quantum,tmpArray,p);
		}
		if (p>0) Arrays.sort(tmpArray,0,p+1);
		read2characters_new.write(tmpArray[0]+"");
		for (i=1; i<=p; i++) {
			if (tmpArray[i]!=tmpArray[i-1]) read2characters_new.write((SEPARATOR_MINOR+"")+tmpArray[i]);
		}
		lastLeft=-1;
	}
	
	
	private static final void writeBlock(int blockID, int length, boolean openStart, boolean openEnd, int quantum, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, BufferedWriter read2characters_new, int[] tmpArray) throws IOException {
		int i, c, p;
		final int last = lastInBlock[blockID];
		
		p=-1;
		for (i=0; i<=last; i++) {			
			c=Integer.parseInt(blocks[blockID][i]);
			if (c<0) c=-1-c;
			if (c==lastAlphabet_old+1) tmpArray[++p]=lastAlphabet_new+1;
			else {
				tmpCharacter.copyFrom(oldAlphabet[c]);		
				if (c>lastUnique_old && c<=lastPeriodic_old) {
					tmpCharacter.length=length;
					if (tmpCharacter.orientation) {
						tmpCharacter.openStart=openStart;
						tmpCharacter.openEnd=openEnd;
					}
					else {
						tmpCharacter.openStart=openEnd;
						tmpCharacter.openEnd=openStart;
					}
					tmpCharacter.quantize(quantum);
				}
				p=fixPeriodicEndpoints_updateTranslation_impl(tmpCharacter,newAlphabet,lastUnique_new,lastPeriodic_new,lastAlphabet_new,quantum,tmpArray,p);
			}
		}
		if (p>0) Arrays.sort(tmpArray,0,p+1);
		read2characters_new.write(tmpArray[0]+"");
		for (i=1; i<=p; i++) {
			if (tmpArray[i]!=tmpArray[i-1]) read2characters_new.write((SEPARATOR_MINOR+"")+tmpArray[i]);
		}
	}
	
	
	/**
	 * @param quantum if $character$ is not found in $newAlphabet$, the procedure tries
	 * again by changing the length of $character$ by $quantum$. This is useful since, due
	 * to quantization, even the length of a character that was not altered by 
	 * $fixPeriodicEndpoints_updateTranslation()$, where the length is derived from the 
	 * boundaries file, might be slightly different from the length of the character
	 * pointed to by the block in the translation file.
	 * @param out the procedure appends characters to $out[last+1..]$;
	 * @return the new value of $last$ after appending.
	 */
	private static final int fixPeriodicEndpoints_updateTranslation_impl(Character character, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, int quantum, int[] out, int last) {
		int p;
		
		if (character.repeat==UNIQUE) {
			p=Arrays.binarySearch(newAlphabet,0,lastUnique_new+1,character);
			if (p<0) {
				character.length+=quantum;
				p=Arrays.binarySearch(newAlphabet,0,lastUnique_new+1,character);
				character.length-=quantum;
			}
			if (p<0) {
				character.length-=quantum;
				p=Arrays.binarySearch(newAlphabet,0,lastUnique_new+1,character);
				character.length+=quantum;
			}
			if (p<0) {
				System.err.println("fixPeriodicEndpoints_updateTranslation_impl> ERROR: unique character not found in the new alphabet\n query: "+character+"\n first candidate in the new alphabet: "+newAlphabet[p]);
				System.exit(1);
			}
			out[++last]=p;
		}
		else if (character.start==-1) {
			p=fixPeriodicEndpoints_lookupPeriodic(character,newAlphabet,lastUnique_new,lastPeriodic_new);
			if (p<0) {
				character.length+=quantum;
				p=fixPeriodicEndpoints_lookupPeriodic(character,newAlphabet,lastUnique_new,lastPeriodic_new);
				character.length-=quantum;
			}
			if (p<0) {
				character.length-=quantum;
				p=fixPeriodicEndpoints_lookupPeriodic(character,newAlphabet,lastUnique_new,lastPeriodic_new);
				character.length+=quantum;
			}
			if (p<0) {
				System.err.println("fixPeriodicEndpoints_updateTranslation_impl> ERROR: periodic character not found in the new alphabet\n query: "+character+"\n first candidate in the new alphabet: "+newAlphabet[p]);
				System.exit(1);
			}
			out[++last]=p;
		}
		else {
			p=Arrays.binarySearch(newAlphabet,lastPeriodic_new+1,lastAlphabet_new+1,character);
			if (p<0) {
				System.err.println("fixPeriodicEndpoints_updateTranslation_impl> ERROR: nonperiodic character not found in the new alphabet\n query: "+character+"\n first candidate in the new alphabet: "+newAlphabet[-1-p]);
				System.exit(1);
			}
			out[++last]=p;
		}
		return last;
	}
	
	
	/**
	 * @return if $character$ occurs in $newAlphabet$, the position of it; otherwise, the 
	 * position of the smallest element of $newAlphabet$ that implies $character$; -1
	 * if no element of $newAlphabet$ implies $character$.
	 */
	private static final int fixPeriodicEndpoints_lookupPeriodic(Character character, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new) {
		boolean orientation;
		int k, p;
		int repeat, length;
		
		p=Arrays.binarySearch(newAlphabet,lastUnique_new+1,lastPeriodic_new+1,character);
		if (p>=0) return p;
		repeat=character.repeat; orientation=character.orientation; length=character.length;
		p=-1-p;
		if (newAlphabet[p].repeat!=repeat || newAlphabet[p].orientation!=orientation) p--;
		k=p;
		while (k>lastUnique_new && newAlphabet[k].repeat==repeat && newAlphabet[k].orientation==orientation && newAlphabet[k].length>=length) k--;
		k++;
		while (k<=lastPeriodic_new && newAlphabet[k].repeat==repeat && newAlphabet[k].orientation==orientation && !newAlphabet[k].implies(character)) k++;
		p=k;
		return p;
	}
	
	
	/**
	 * @return if $character$ occurs in $newAlphabet$, the position of it; otherwise, the 
	 * position of the smallest element of $newAlphabet$ that implies $character$; -1
	 * if no element of $newAlphabet$ implies $character$.
	 */
	private static final int fixPeriodicEndpoints_lookupUnique(Character character, Character[] newAlphabet, int lastUnique_new, int lastAlphabet_new) {
		int p;
		
		p=Arrays.binarySearch(newAlphabet,0,lastUnique_new+1,character);
		if (p>=0) return p;
		p=-1-p;
		if (newAlphabet[p].repeat!=UNIQUE) return lastAlphabet_new+1;
		else return p;
	}
	
	
	
	
	// ---------------------------------- WOBBLING ---------------------------------------
	
	/**
	 * Sets $out[i]=true$ iff $alphabet[i]$ is a non-periodic character that has to be
	 * wobbled (every periodic character must be wobbled, marking those is not needed).
	 *
	 * @param read2characters a row of the translated reads file;
	 * @param tmpArray temporary space, of size at least equal to the number of blocks.
	 */
	public static final void wobble_markAlphabet(String read2characters, boolean[] out, int[] tmpArray) throws IOException {
		int i, j, k, c;
		int last, nBlocks, nPeriodicBlocks;
		
		// Marking blocks to be wobbled
		if (read2characters.length()==0) return;
		nBlocks=loadBlocks(read2characters);
		if (nBlocks==1) return;
		Math.set(tmpArray,nBlocks-1,0);
		nPeriodicBlocks=0;
		for (i=0; i<nBlocks; i++) {
			last=lastInBlock[i];
			for (j=0; j<=last; j++) {
				c=Integer.parseInt(blocks[i][j]);
				if (c<0) c=-1-c;
				if (c>lastUnique && c<=lastPeriodic) {
					tmpArray[i]=1;
					nPeriodicBlocks++;
					break;
				}
			}
		}
		if (nPeriodicBlocks==0) return;
		if (tmpArray[0]==1 && tmpArray[1]==0) tmpArray[1]=2;
		for (i=1; i<nBlocks-1; i++) {
			if (tmpArray[i]!=1) continue;
			if (tmpArray[i-1]==0) tmpArray[i-1]=2;
			if (tmpArray[i+1]==0) tmpArray[i+1]=2;
		}
		if (tmpArray[nBlocks-1]==1 && tmpArray[nBlocks-2]==0) tmpArray[nBlocks-2]=2;
		
		// Marking every non-periodic character as to be wobbled.
		for (i=0; i<nBlocks; i++) {
			if (tmpArray[i]==1) {
				for (j=0; j<=lastInBlock[i]; j++) {
					c=Integer.parseInt(blocks[i][j]);
					if (c>=0 && c>lastPeriodic) out[c]=true;
				}
			}
			else if (tmpArray[i]==2) {
				last=-1;
				for (j=0; j<=lastInBlock[i]; j++) {
					c=Integer.parseInt(blocks[i][j]);
					if (c<0) {  // Must be a unique character
						c=-1-c;
						if (IO.CONSISTENCY_CHECKS && c>lastUnique) {
							System.err.println("wobble_markAlphabet> ERROR: negative character is not unique: "+(-1-c)+" lastUnique="+lastUnique);
							System.exit(1);
						}
						for (k=c; k<=lastUnique; k++) out[k]=true;
					}
					else out[c]=true;
				}
			}
		}
	}
	
	
	/**
	 * Builds the new alphabet that results from wobbling the characters of the current 
	 * alphabet that are marked in $flags$, as well as every periodic character. Each 
	 * section of the new alphabet (unique, periodic, non-periodic) might not be sorted 
	 * and might contain duplicates.
	 *
	 * @param nFlags number of TRUE elements in $flags$, excluding periodic characters;
	 * @param out output array: contains the values of $lastUnique,lastPeriodic,
	 * lastAlphabet$ for the new alphabet.
	 */
	public static final Character[] wobble_extendAlphabet(boolean[] flags, int nFlags, int quantum_wobble, int quantum_alphabet, int[] out) throws IOException {
		final int MAX_NEWCHARS_PER_CHAR = ((quantum_wobble<<1)/quantum_alphabet)*((quantum_wobble<<1)/quantum_alphabet);
		int i;
		int lastUnique_new, lastPeriodic_new, lastAlphabet_new;
		Character tmpCharacter;
		Character[] alphabet_new;
		
		tmpCharacter = new Character();
		alphabet_new = new Character[lastAlphabet+(nFlags+lastPeriodic-lastUnique)*MAX_NEWCHARS_PER_CHAR];
		lastUnique_new=-1;
		for (i=0; i<=lastUnique; i++) {
			alphabet_new[++lastUnique_new]=alphabet[i];
			if (flags[i]) lastUnique_new=wobble_extendAlphabet_impl(i,quantum_wobble,quantum_alphabet,alphabet_new,lastUnique_new,tmpCharacter);
		}
		lastPeriodic_new=lastUnique_new;
		for (i=lastUnique+1; i<=lastPeriodic; i++) {
			alphabet_new[++lastPeriodic_new]=alphabet[i];
			lastPeriodic_new=wobble_extendAlphabet_impl(i,quantum_wobble,quantum_alphabet,alphabet_new,lastPeriodic_new,tmpCharacter);
		}
		lastAlphabet_new=lastPeriodic_new;
		for (i=lastPeriodic+1; i<=lastAlphabet; i++) {
			alphabet_new[++lastAlphabet_new]=alphabet[i];
			if (flags[i]) lastAlphabet_new=wobble_extendAlphabet_impl(i,quantum_wobble,quantum_alphabet,alphabet_new,lastAlphabet_new,tmpCharacter);
		}
		out[0]=lastUnique_new; out[1]=lastPeriodic_new; out[2]=lastAlphabet_new;
		return alphabet_new;
	}
	
	
	/**
	 * Adds to $alphabet_new[lastCharacter_new+1..]$ every character that is produced by
	 * wobbling $alphabet[c]$ and that is not already in $alphabet$.
	 *
	 * Remark: running this procedure with different values of $c$ might introduce 
	 * duplicates in $alphabet_new$.
	 *
	 * @param alphabet_new assumed to be large enough;
	 * @return the new value of $lastCharacter_new$ after the procedure completes.
	 */
	private static final int wobble_extendAlphabet_impl(int c, int quantum_wobble, int quantum_alphabet, Character[] alphabet_new, int lastCharacter_new, Character tmpCharacter) {
		int i, j;
		int last, start;
		final int startPrime = alphabet[c].start;
		final int length = alphabet[c].getLength();
		Character newCharacter;
		
		if (c<=lastPeriodic) {
			for (i=quantum_alphabet; i<=quantum_wobble; i+=quantum_alphabet) {
				tmpCharacter.copyFrom(alphabet[c]);
				tmpCharacter.length=length-i;
				if (!wobble_find(tmpCharacter,c,false)) {
					newCharacter = new Character();
					newCharacter.copyFrom(tmpCharacter);
					alphabet_new[++lastCharacter_new]=newCharacter;
				}
				tmpCharacter.length=length+i;
				if (!wobble_find(tmpCharacter,c,true)) {
					newCharacter = new Character();
					newCharacter.copyFrom(tmpCharacter);
					alphabet_new[++lastCharacter_new]=newCharacter;
				}
			}
		}
		else {
			for (i=quantum_alphabet; i<=quantum_wobble; i+=quantum_alphabet) {
				start=startPrime-i;
				for (j=quantum_alphabet; j<=quantum_wobble; j+=quantum_alphabet) {
					tmpCharacter.copyFrom(alphabet[c]);
					tmpCharacter.start=start;
					tmpCharacter.end=tmpCharacter.start+(length-j)-1;
					if (!wobble_find(tmpCharacter,c,false)) {
						newCharacter = new Character();
						newCharacter.copyFrom(tmpCharacter);
						alphabet_new[++lastCharacter_new]=newCharacter;
					}
					tmpCharacter.end=tmpCharacter.start+(length+j)-1;
					if (!wobble_find(tmpCharacter,c,false)) {
						newCharacter = new Character();
						newCharacter.copyFrom(tmpCharacter);
						alphabet_new[++lastCharacter_new]=newCharacter;
					}
				}
				start=startPrime+i;
				for (j=quantum_alphabet; j<=quantum_wobble; j+=quantum_alphabet) {
					tmpCharacter.copyFrom(alphabet[c]);
					tmpCharacter.start=start;
					tmpCharacter.end=tmpCharacter.start+(length-j)-1;
					if (!wobble_find(tmpCharacter,c,true)) {
						newCharacter = new Character();
						newCharacter.copyFrom(tmpCharacter);
						alphabet_new[++lastCharacter_new]=newCharacter;
					}
					tmpCharacter.end=tmpCharacter.start+(length+j)-1;
					if (!wobble_find(tmpCharacter,c,true)) {
						newCharacter = new Character();
						newCharacter.copyFrom(tmpCharacter);
						alphabet_new[++lastCharacter_new]=newCharacter;
					}
				}
			}
		}
		return lastCharacter_new;
	}
	
	
	/**
	 * @param c searches for an exact occurrence of $character$ starting from position $c$
	 * in $alphabet$ (excluded);
	 * @param direction TRUE=forward, FALSE=backward.
	 */
	private static final boolean wobble_find(Character character, int c, boolean direction) {
		int i;
		int from, to;
		final boolean orientation = character.orientation;
		final int repeat = character.repeat;
		final int start = character.start;
		
		if (direction) {
			from=c+1;
			if (c<=lastUnique) to=lastUnique;
			else if (c<=lastPeriodic) to=lastPeriodic;
			else to=lastAlphabet;
			for (i=from; i<=to; i++) {
				if (alphabet[i].repeat!=repeat || alphabet[i].orientation!=orientation || alphabet[i].start>start) break;
				if (alphabet[i].equals(character)) return true;
			}
		}
		else {
			to=c-1;
			if (c<=lastUnique) from=0;
			else if (c<=lastPeriodic) from=lastUnique+1;
			else from=lastPeriodic+1;
			for (i=to; i>=from; i--) {
				if (alphabet[i].repeat!=repeat || alphabet[i].orientation!=orientation || alphabet[i].start<start) break;
				if (alphabet[i].equals(character)) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Stores in $outputFile$ the map oldCharacter -> newCharacter.
	 */
	public static final void wobble_buildOld2New(Character[] alphabet_old, int lastAlphabet_old, Character[] alphabet_new, int lastAlphabet_new, String outputFile) throws IOException {
		int i, j, k;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(outputFile));
		i=0; j=0;
		while (i<=lastAlphabet_old) {
			k=alphabet_new[j].compareTo(alphabet_old[i]);
			if (k<0) { j++; continue; }
			else if (k==0) { bw.write(j+""); i++; j++; }
			else {  
				System.err.println("wobble_buildOld2New> ERROR: old character not found in the new alphabet: "+alphabet_old[i]);
				System.exit(1);
			}
		}
		bw.close();
	}
	
	
	/**
	 * Rewrites the translation of a read so that periodic blocks, and blocks adjacent to
	 * periodic blocks, \emph{wobble}, i.e. they contain, in addition to their original
	 * characters in $alphabet$, other similar characters whose length is $<=wobbleLength$
	 * bps away from the original. Reads with just one block are not altered.
	 *
	 * Remark: wobbling is designed to increase the number of edges in a highly
	 * disconnected overlap graph where the endpoints of repeat occurrences are uncertain.
	 * However, such increase in frequency may make some k-mers be classified as repeats 
	 * rather than as unique addresses on the genome, and this might \emph{remove} some
	 * edges from the overlap graph.
	 *
	 * Remark: this procedure might put multiple unique characters inside a block that
	 * contains a single unique character. Thus, throughout the code, no test for the
	 * non-repetitiveness of a block must require a block to contain a single character.
	 *
	 * @param read2characters row of the translated reads file;
	 * @param alphabet_new the new alphabet that contains the results of wobbling; the 
	 * procedure translates $read2characters$ to this alphabet;
	 * @param old2new map $alphabet -> alphabet_new$;
	 * @param output array; the procedure cumulates to $out[0]$ the number of blocks to 
	 * which wobbling was applied, and to $out[1]$ the total number of blocks in the read.
	 */
	public static final void wobble(String read2characters, int wobbleLength, int[] old2new, Character[] alphabet_new, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, BufferedWriter bw, int[] tmpArray1, int[] tmpArray2, int[] out) throws IOException {
		int i, j, c;
		int last, nBlocks, nPeriodicBlocks;
		
		// Marking blocks to be wobbled
		if (read2characters.length()==0) {
			bw.newLine();
			return;
		}
		nBlocks=loadBlocks(read2characters);
		out[1]+=nBlocks;
		Math.set(tmpArray1,nBlocks-1,0);
		if (nBlocks>1) {
			nPeriodicBlocks=0;
			for (i=0; i<nBlocks; i++) {
				last=lastInBlock[i];
				for (j=0; j<=last; j++) {
					c=Integer.parseInt(blocks[i][j]);
					if (c<0) c=-1-c;
					if (c>lastUnique && c<=lastPeriodic) {
						tmpArray1[i]=1;
						nPeriodicBlocks++;
						break;
					}
				}
			}
			if (nPeriodicBlocks>0) {
				if (tmpArray1[0]==1 && tmpArray1[1]==0) tmpArray1[1]=2;
				for (i=1; i<nBlocks-1; i++) {
					if (tmpArray1[i]!=1) continue;
					if (tmpArray1[i-1]==0) tmpArray1[i-1]=2;
					if (tmpArray1[i+1]==0) tmpArray1[i+1]=2;
				}
				if (tmpArray1[nBlocks-1]==1 && tmpArray1[nBlocks-2]==0) tmpArray1[nBlocks-2]=2;
			}
		}
		
		// Wobbling blocks
		for (i=0; i<nBlocks; i++) {
			if (i>0) bw.write(SEPARATOR_MAJOR+"");
			if (tmpArray1[i]==0) {
				c=Integer.parseInt(blocks[i][0]);
				if (c>=0) bw.write(old2new[c]+"");
				else bw.write((-1-old2new[-1-c])+"");
				for (j=1; j<=lastInBlock[i]; j++) {
					c=Integer.parseInt(blocks[i][j]);
					if (c>=0) bw.write(SEPARATOR_MINOR+""+old2new[c]);
					else bw.write(SEPARATOR_MINOR+""+(-1-old2new[-1-c]));
				}
			}
			else {
				last=-1;
				for (j=0; j<=lastInBlock[i]; j++) last=wobble_impl(Integer.parseInt(blocks[i][j]),wobbleLength,old2new,alphabet_new,lastUnique_new,lastPeriodic_new,lastAlphabet_new,tmpArray2,last);
				if (last>0) Arrays.sort(tmpArray2,0,last+1);
				j=0; bw.write(tmpArray2[0]+"");
				for (j=1; j<=last; j++) {
					if (tmpArray2[j]!=tmpArray2[j-1]) bw.write(SEPARATOR_MINOR+""+tmpArray2[j]);
				}
				out[0]++;
			}
		}
		bw.newLine();
	}
	
	
	/**
	 * Let $c$ be the ID of a character in $alphabet$. The procedure translates $c$ to 
	 * $alphabet_new$, appends it to $out[outLast+1..]$, and then appends all the other 
	 * characters of $alphabet_new$ that are the wobble of $alphabet[c]$. Characters are 
	 * appended in no particular order.
	 * 
	 * @param c if negative, the characters added to $outLast$ are negative;
	 * @return the new value of $outLast$ after the procedure completes.
	 */
	private static final int wobble_impl(int c, int quantum, int[] old2new, Character[] alphabet_new, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, int[] out, int outLast) {
		int i, d, w;
		int cPrime, from, to;
		
		cPrime=old2new[c<0?-1-c:c];
		if (cPrime<=lastUnique_new) { from=0; to=lastUnique_new; }
		else if (cPrime>lastUnique_new && cPrime<=lastPeriodic_new) { from=lastUnique_new+1; to=lastPeriodic_new; }
		else if (cPrime<=lastAlphabet_new) { from=lastPeriodic_new+1; to=lastAlphabet_new; }
		else {
			out[++outLast]=c<0?-1-(lastAlphabet_new+1):lastAlphabet_new+1;
			return outLast;
		}
		d=cPrime-1;
		while (d>=from) {
			w=isWobbleOf(d,cPrime,quantum,alphabet_new);
			if (w==-1) break;
			else if (w==1) out[++outLast]=c<0?-1-d:d;
			d--;
		}
		out[++outLast]=c<0?-1-cPrime:cPrime;
		d=cPrime+1;
		while (d<=to) {
			w=isWobbleOf(d,cPrime,quantum,alphabet_new);
			if (w==-1) break;
			else if (w==1) out[++outLast]=c<0?-1-d:d;
			d++;
		}
		return outLast;
	}
	
	
	/**
	 * @return 1=$alphabet[x]$ is identical to $alphabet[y]$ in everything except the 
	 * length, which differs by at most $quantum$, and the starting position, which 
	 * differs by at most $quantum$; -1=$alphabet[x]$ is out of the range of candidates 
	 * that could satisfy (1); 0=none of the above.
	 */
	private static final int isWobbleOf(int x, int y, int quantum, Character[] alphabet) {
		if ( alphabet[x].repeat!=alphabet[y].repeat || alphabet[x].orientation!=alphabet[y].orientation || 
			 alphabet[x].start<alphabet[y].start-quantum || alphabet[x].start>alphabet[y].start+quantum ) return -1;
		if (alphabet[x].openStart!=alphabet[y].openStart || alphabet[x].openEnd!=alphabet[y].openEnd) return 0;
		int lengthX, lengthY;
		if (alphabet[x].start==-1) {
			lengthX=alphabet[x].length;
			lengthY=alphabet[y].length;
		}
		else {
			lengthX=alphabet[x].end-alphabet[x].start+1;
			lengthY=alphabet[y].end-alphabet[y].start+1;
		}
		return Math.abs(lengthX,lengthY)<=quantum?1:0;
	}
	
	
	public static class Spacer implements Comparable {
		boolean rigidLeft, rigidRight;
		public int read, first, last, breakpoint;
		
		/**
		 * Index of the non-repetitive block (if any) in the translation of $read$.
		 */
		public int blockID;
		
		public Spacer() { }
		
		public Spacer(int r, int f, int l, boolean fl, boolean fr, int b) {
			this.read=r; this.first=f; this.last=l;
			this.rigidLeft=fl; this.rigidRight=fr;
			this.blockID=b;
			breakpoint=-1;
		}
		
		public boolean isRigid() { return rigidLeft||rigidRight; }
		
		public void setBreakpoint() {
			if (first==last) breakpoint=first;
			else if (rigidLeft) breakpoint=first;
			else if (rigidRight) breakpoint=last;
			else breakpoint=(first+last)>>1;  // Arbitrary
		}
		
		/**
		 * By $read$ only.
		 */
		public int compareTo(Object other) {
			Spacer otherSpacer = (Spacer)other;
			if (read<otherSpacer.read) return -1;
			else if (read>otherSpacer.read) return 1;
			return 0;
		}
		
		public String toString() { return read+"["+first+".."+breakpoint+".."+last+"] ("+blockID+")"; }
		
		public void serialize(BufferedWriter bw) throws IOException {
			bw.write(read+""+SEPARATOR_MINOR);
			bw.write(first+""+SEPARATOR_MINOR);
			bw.write(last+""+SEPARATOR_MINOR);
			bw.write(breakpoint+""+SEPARATOR_MINOR);
			bw.write(blockID+""+SEPARATOR_MINOR);
			bw.write((rigidLeft?"1":"0")+SEPARATOR_MINOR);
			bw.write((rigidRight?"1":"0")+(SEPARATOR_MINOR+"\n"));
		}
		
		public void deserialize(String str) throws IOException {
			String[] tokens = str.split(SEPARATOR_MINOR+"");
			read=Integer.parseInt(tokens[0]);
			first=Integer.parseInt(tokens[1]);
			last=Integer.parseInt(tokens[2]);
			breakpoint=Integer.parseInt(tokens[3]);
			blockID=Integer.parseInt(tokens[4]);
			rigidLeft=Integer.parseInt(tokens[5])==1;
			rigidRight=Integer.parseInt(tokens[6])==1;
		}
	}
	






	
	// ------------------------------ DATA STRUCTURES ------------------------------------
	
	/**
	 * A block of a recoded read, which is a set of characters.
	 */
	public static class Block {
		private static final int CAPACITY = 10;  // Arbitrary
				
		public Character[] characters;
		public int lastCharacter;
		
		/**
		 * An empty block
		 */
		public Block() {
			characters = new Character[CAPACITY];
			for (int i=0; i<characters.length; i++) characters[i] = new Character();
			reset();
		}
		
		public void reset() { lastCharacter=-1; }
		
		/**
		 * A non-repetitive block of given length.
		 */
		public Block(int length, boolean openStart, boolean openEnd) {
			characters = new Character[CAPACITY];
			setUnique(length,openStart,openEnd);
		}
		
		/**
		 * Sets the current block to a non-repetitive of given length.
		 */
		public final void setUnique(int length, boolean openStart, boolean openEnd) {
			characters[0] = new Character(UNIQUE,true,-1,-1,length,openStart,openEnd);
			lastCharacter=0;
		}
		
		public final boolean isUnique() {
			return lastCharacter==0 && characters[0].repeat==UNIQUE;
		}
		
		private final void enlarge() {
			Character[] newCharacters = new Character[characters.length<<1];
			System.arraycopy(characters,0,newCharacters,0,characters.length);
			for (int i=characters.length; i<newCharacters.length; i++) newCharacters[i] = new Character();
			characters=newCharacters;
		}
		
		/**
		 * Makes this block contain the same data as $otherBlock$ (but all 
		 * pointers stay distinct). Field $id$ is not changed.
		 */
		public final void copyFrom(Block otherBlock) {
			int i;
			
			if (characters==null) {
				characters = new Character[otherBlock.lastCharacter+1];
				for (i=0; i<characters.length; i++) characters[i] = new Character();
			}
			else if (characters.length<otherBlock.lastCharacter+1) {
				Character[] newCharacters = new Character[otherBlock.lastCharacter+1];
				System.arraycopy(characters,0,newCharacters,0,characters.length);
				for (i=characters.length; i<newCharacters.length; i++) newCharacters[i] = new Character();
				characters=newCharacters;
			}
			lastCharacter=otherBlock.lastCharacter;
			for (i=0; i<=lastCharacter; i++) characters[i].copyFrom(otherBlock.characters[i]);
		}
		
		
		public String toString() {
			String out = (isUnique()?"0":"1")+(SEPARATOR_MINOR+"")+(characters[0].start==-1?"0":"1")+(SEPARATOR_MINOR+"")+lastCharacter+(SEPARATOR_MINOR+"");
			out+=characters[0].toString();
			for (int i=1; i<=lastCharacter; i++) out+=SEPARATOR_MINOR+""+characters[i].toString();
			return out;
		}
		
		
		/**
		 * Symmetrical to $toString()$.
		 */
		public void deserialize(String str) {
			int i, p, q;
			
			p=str.indexOf(SEPARATOR_MINOR+"");  // Skipping
			p=str.indexOf(SEPARATOR_MINOR+"",p+1);  // Skipping
			q=str.indexOf(SEPARATOR_MINOR+"",p+1);
			lastCharacter=Integer.parseInt(str.substring(p+1,q));
			if (characters==null) {
				characters = new Character[lastCharacter+1];
				for (i=0; i<=lastCharacter; i++) characters[i] = new Character();
			}
			else if (characters.length<lastCharacter+1) {
				Character[] newArray = new Character[lastCharacter+1];
				System.arraycopy(characters,0,newArray,0,characters.length);
				for (i=characters.length; i<=lastCharacter; i++) newArray[i] = new Character();
				characters=newArray;
			}
			p=q+1;
			for (i=0; i<=lastCharacter; i++) p=characters[i].deserialize(str);
		}
		
		
		/**
		 * @return the average length of all characters.
		 */
		public int getLength() {
			int i, out;
			
			out=0;
			for (i=0; i<=lastCharacter; i++) out+=characters[i].getLength();
			return out/(lastCharacter+1);
		}
		
		
		/**
		 * Adds an alignment to the current block. If the block already contains an 
		 * alignment to the same repeat, the new alignment replaces the old one based
		 * on a canonical (and arbitrary) order. The alignment might contain the block and
		 * be much larger than it.
		 *
		 * Remark: characters might not be sorted after the procedure completes.
		 *
		 * @param alignmentStartA,alignmentEndA substring of $[alignment.startA..
		 * alignment.endA]$ that corresponds to this block of readA (the alignment might
		 * span multiple blocks);
		 * @param tmpCharacter temporary space.
		 */
		public final void addCharacter(AlignmentRow alignment, int distanceThreshold, int alignmentStartA, int alignmentEndA, Character tmpCharacter) {
			int i;
			int found, alignmentStartB, alignmentEndB;
			final double ratio = ((double)(alignment.endB-alignment.startB+1))/(alignment.endA-alignment.startA+1);
			
			alignmentStartB=-1; alignmentEndB=-1;
			tmpCharacter.repeat=alignment.readB; tmpCharacter.orientation=alignment.orientation;
			if (isPeriodic[alignment.readB]) {
				if (alignment.orientation) {
					tmpCharacter.openStart=alignmentStartA<=distanceThreshold;
					tmpCharacter.openEnd=alignmentEndA>=Reads.getReadLength(alignment.readA)-distanceThreshold;
				}
				else {
					tmpCharacter.openEnd=alignmentStartA<=distanceThreshold;
					tmpCharacter.openStart=alignmentEndA>=Reads.getReadLength(alignment.readA)-distanceThreshold;
				}
				tmpCharacter.start=-1; tmpCharacter.end=-1;
				tmpCharacter.length=alignmentEndA-alignmentStartA+1; // A is correct here
			}
			else {
				if (alignment.orientation) {
					alignmentStartB=alignment.startB+(alignmentStartA>alignment.startA?(int)((alignmentStartA-alignment.startA)*ratio):0);
					alignmentEndB=alignment.endB-(alignment.endA>alignmentEndA?(int)((alignment.endA-alignmentEndA)*ratio):0);
					tmpCharacter.openStart=alignmentStartA<=distanceThreshold && alignmentStartB>distanceThreshold;
					tmpCharacter.openEnd=alignmentEndA>=Reads.getReadLength(alignment.readA)-distanceThreshold && alignmentEndB<repeatLengths[alignment.readB]-distanceThreshold;
				}
				else {
					alignmentStartB=alignment.startB+(alignment.endA>alignmentEndA?(int)((alignment.endA-alignmentEndA)*ratio):0);
					alignmentEndB=alignment.endB-(alignmentStartA>alignment.startA?(int)((alignmentStartA-alignment.startA)*ratio):0);
					tmpCharacter.openEnd=alignmentStartA<=distanceThreshold && alignmentEndB<repeatLengths[alignment.readB]-distanceThreshold;
					tmpCharacter.openStart=alignmentEndA>=Reads.getReadLength(alignment.readA)-distanceThreshold && alignmentStartB>distanceThreshold;
				}
				tmpCharacter.start=alignmentStartB; tmpCharacter.end=alignmentEndB;
				tmpCharacter.length=0;
			}
			
			found=-1;
			for (i=0; i<=lastCharacter; i++) {
				if (characters[i].repeat==tmpCharacter.repeat && characters[i].orientation==tmpCharacter.orientation) {
					found=i;
					break;
				}
			}
			if (found==-1) {
				lastCharacter++;
				if (lastCharacter==characters.length) enlarge();	
				characters[lastCharacter].copyFrom(tmpCharacter);
			}
			else {
				if (isPeriodic[alignment.readB]) characters[found].length=Math.max(characters[found].length,tmpCharacter.length);
				else if ( alignmentStartB<characters[found].start || 
						  (alignmentStartB==characters[found].start && alignmentEndB<characters[found].end)
			        	) {
					characters[found].start=alignmentStartB;
					characters[found].end=alignmentEndB;
				}
				if (tmpCharacter.openStart) characters[found].openStart=true;
				if (tmpCharacter.openEnd) characters[found].openEnd=true;
			}
		}
		
		
		/**
		 * (1) Sorts $characters$. (2) If the block contains a periodic repeat, removes 
		 * all non-periodic repeats from the block and sets the length of every periodic 
		 * repeat to $periodicLength$. (3) If the block is not the first or last of its 
		 * read, all open marks are removed. Then, if the block contains several 
		 * characters from the same periodic repeat in the same orientation, only a single
		 * one (closed if possible) is kept. If the block contains characters from the 
		 * same nonperiodic repeat, in the same orientation, and with similar substrings, 
		 * only a single one (closed if possible) is kept.
		 *
		 * @param stack temporary space, of size at least $lastCharacter+1$.
		 */
		public final void cleanCharacters(int periodicLength, boolean firstOrLast, int[] stack) {
			final int DISTANCE_THRESHOLD = 100;  // Arbitrary
			boolean foundPeriodic, foundNonperiodic;
			int i, j;
			int top, selected, count, countSelected, currentCharacter;
			Character tmpCharacter;
			
			// Constraint 2
			foundPeriodic=false; foundNonperiodic=false;
			for (i=0; i<=lastCharacter; i++) {
				if (isPeriodic[characters[i].repeat]) {
					foundPeriodic=true;
					characters[i].length=periodicLength;
				}
				else foundNonperiodic=true;
			}
			if (foundPeriodic && foundNonperiodic) {
				j=-1;
				for (i=0; i<=lastCharacter; i++) {
					if (!isPeriodic[characters[i].repeat]) continue;
					j++;
					tmpCharacter=characters[j];
					characters[j]=characters[i];
					characters[i]=tmpCharacter;
				}
				lastCharacter=j;
			}
			if (lastCharacter>0) Arrays.sort(characters,0,lastCharacter+1);
			
			// Constraint 3 (using temporary field $Character.flag$).
			if (!firstOrLast) {
				for (i=0; i<=lastCharacter; i++) { characters[i].openStart=false; characters[i].openEnd=false; }
			}
			if (lastCharacter>0) {
				for (i=0; i<=lastCharacter; i++) characters[i].flag=-1;
				for (i=0; i<=lastCharacter; i++) {
					if (characters[i].flag!=-1) continue;
					characters[i].flag=i; top=0; stack[0]=i;
					if (characters[i].isOpen()) { selected=-1; countSelected=0; }
					else { 
						selected=i;
						countSelected=(characters[i].openStart?0:1)+(characters[i].openEnd?0:1);
					}
					while (top>=0) {
						currentCharacter=stack[top--];
						for (j=currentCharacter+1; j<=lastCharacter; j++) {
							if (characters[j].repeat!=characters[currentCharacter].repeat || characters[j].orientation!=characters[currentCharacter].orientation) break;
							if (characters[j].flag!=-1) continue;
							if ( characters[j].start==-1 ||
						   	 	 ( characters[j].start<=characters[currentCharacter].start+DISTANCE_THRESHOLD &&
						           Math.abs(characters[j].end,characters[currentCharacter].end)<=DISTANCE_THRESHOLD
						         )
							   ) {
								characters[j].flag=i;
								count=(characters[j].openStart?0:1)+(characters[j].openEnd?0:1);
								if ( selected==-1 || count>countSelected || 
									 ( count==selected &&
									   ( (characters[j].start==-1 && characters[j].length<characters[selected].length) ||
									     (characters[j].start!=-1 && characters[j].end-characters[j].start<characters[selected].end-characters[selected].start)
									   )
									 )
								   ) {
									selected=j; countSelected=count;
									continue;
								}
								stack[++top]=j;
							}
						}
					}
					if (selected==-1 || selected==i) {
						// If all characters are open, we leave them flagged by $i$.
					}
					else {
						for (j=i; j<=lastCharacter; j++) {
							if (characters[j].repeat!=characters[i].repeat || characters[j].orientation!=characters[i].orientation) break;
							if (characters[j].flag==i) characters[j].flag=selected;
						}
					}
				}
				j=-1;
				for (i=0; i<=lastCharacter; i++) {
					if (characters[i].flag!=i) continue;
					j++;
					tmpCharacter=characters[j];
					characters[j]=characters[i];
					characters[i]=tmpCharacter;
				}
				lastCharacter=j;
			}
		}		
		
	}
	
	
	/**
	 * A substring of a repeat in a specific orientation and with open/closed endpoints.
	 */
	public static class Character implements Comparable {
		public int repeat;
		public boolean orientation;
		public int start, end;
		
		/**
		 * Tell whether the repeat might continue after its start/end.
		 *
		 * Remark: these flags refer to the positions written in $start,end$, so they do
		 * not need to be reversed when reverse-complementing the character.
		 */
		public boolean openStart, openEnd;
		
		public int length;  // >0: unique or periodic; 0: repetitive non-periodic.
		
		public int flag;  // Temporary space
		
		
		public Character() { reset(); }
		
		
		public void reset() { 
			repeat=-1; orientation=false; start=-1; end=-1; length=-1;
			openStart=true; openEnd=true; flag=-1;
		}
		
		
		public Character(int r, boolean o, int s, int e, int l, boolean os, boolean oe) {
			repeat=r; orientation=o; start=s; end=e; length=l;
			openStart=os; openEnd=oe;
		}
		
		
		public int compareTo(Object other) {
			Character otherCharacter = (Character)other;
			
			if (repeat==UNIQUE && otherCharacter.repeat!=UNIQUE) return -1;
			else if (repeat!=UNIQUE && otherCharacter.repeat==UNIQUE) return 1;
			if (start==-1 && otherCharacter.start!=-1) return -1;
			else if (start!=-1 && otherCharacter.start==-1) return 1;
			
			if (repeat<otherCharacter.repeat) return -1;
			else if (repeat>otherCharacter.repeat) return 1;
			if (orientation && !otherCharacter.orientation) return -1;
			else if (!orientation && otherCharacter.orientation) return 1;
			if (start<otherCharacter.start) return -1;
			else if (start>otherCharacter.start) return 1;
			if (end<otherCharacter.end) return -1;
			else if (end>otherCharacter.end) return 1;
			if (length<otherCharacter.length) return -1;
			else if (length>otherCharacter.length) return 1;
			if (openStart && !otherCharacter.openStart) return -1;
			else if (!openStart && otherCharacter.openStart) return 1;
			if (openEnd && !otherCharacter.openEnd) return -1;
			else if (!openEnd && otherCharacter.openEnd) return 1;
			return 0;
		}
		
		
		public final void copyFrom(Character otherCharacter) {
			repeat=otherCharacter.repeat;
			orientation=otherCharacter.orientation;
			start=otherCharacter.start; end=otherCharacter.end;
			length=otherCharacter.length;
			openStart=otherCharacter.openStart; openEnd=otherCharacter.openEnd;
		}
		
		
		public boolean equals(Object other) {
			Character otherCharacter = (Character)other;
			return repeat==otherCharacter.repeat && orientation==otherCharacter.orientation && start==otherCharacter.start && end==otherCharacter.end && length==otherCharacter.length && openStart==otherCharacter.openStart && openEnd==otherCharacter.openEnd;
		}
		
		
		public String toString() {
			return (repeat==UNIQUE?"0":"1")+SEPARATOR_MINOR+(start==-1?"0":"1")+(SEPARATOR_MINOR+"")+repeat+(SEPARATOR_MINOR+"")+(orientation?"1":"0")+(SEPARATOR_MINOR+"")+start+(SEPARATOR_MINOR+"")+end+(SEPARATOR_MINOR+"")+length+(SEPARATOR_MINOR+"")+(openStart?"1":"0")+(SEPARATOR_MINOR+"")+(openEnd?"1":"0");
		}
		
		
		/**
		 * Symmetrical to $toString()$.
		 *
		 * @param p the first position of $str$ to be read;
		 * @return the value of $p$ when the procedure completes.
		 */
		public int deserialize(String str) {
			int p, q;
			
			p=str.indexOf(SEPARATOR_MINOR+"");  // Skipping
			p=str.indexOf(SEPARATOR_MINOR+"",p+1)+1;  // Skipping
			q=str.indexOf(SEPARATOR_MINOR+"",p+1);
			repeat=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
			orientation=Integer.parseInt(str.substring(p,q))==1;
			p=q+1; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
			start=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
			end=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
			length=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
			openStart=Integer.parseInt(str.substring(p,q))==1;
			p=q+1; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
			if (q>=0) {
				openEnd=Integer.parseInt(str.substring(p,q))==1;
				return q+1;
			}
			else {
				openEnd=Integer.parseInt(str.substring(p))==1;
				return str.length();
			}
		}
		
		
		public final int getLength() {
			return (start==-1||end==-1)?length:end-start+1;
		}
		
		
		/**
		 * Tells whether this character is more specific than, or as specific as, 
		 * $otherCharacter$ (a character implies itself). 
		 *
		 * Remark: an open character is just a character that might continue on its open
		 * sides. After $compactInstances()$, only open characters that cannot be implied
		 * by any other character in the alphabet survive (e.g. long satellites). The open
		 * status of a character should not affect anything downstream: open characters
		 * are just characters like any other.
		 *
		 * Remark: both characters are assumed to be quantized, and repetitive with the 
		 * same repeat and in the same orientation. Unique characters are not assumed to
		 * imply one another, since fully-open and half-open unique characters are 
		 * assumed to have been discarded.
		 */
		public final boolean implies(Character otherCharacter) {
			if (start==-1) {  // Periodic
				if (length<otherCharacter.length) return false;
				else if (length==otherCharacter.length) {
					if ((openStart && !otherCharacter.openStart) || (openEnd && !otherCharacter.openEnd)) return false;
					else return true;
				}
				else {
					if (!openStart && !openEnd && !otherCharacter.openStart && !otherCharacter.openEnd) return false;
					else if ((openStart && !otherCharacter.openStart) || (openEnd && !otherCharacter.openEnd)) return false;
					else return true;
				}
			}
			else {  // Nonperiodic
				if (otherCharacter.start==start && otherCharacter.end==end) {
					if ((openStart && !otherCharacter.openStart) || (openEnd && !otherCharacter.openEnd)) return false;
					else return true;
				}
				else if (Intervals.isContained(otherCharacter.start,otherCharacter.end,start,end)) {
					if ( (start==otherCharacter.start && openStart && !otherCharacter.openStart) ||
						 (otherCharacter.start>start && !otherCharacter.openStart) ||
						 (end==otherCharacter.end && openEnd && !otherCharacter.openEnd) ||
						 (otherCharacter.end<end && !otherCharacter.openEnd)
					   ) return false;
					else return true;
				}
				else return false;
			}
		}
		
		
		/**
		 * Tells whether neither this character nor any character that follows in the 
		 * alphabet can be implied by $otherCharacter$.
		 * 
		 * Remark: the procedure assumes that the characters have the same repeats.
		 */
		public final boolean implies_tooFarAfter(Character otherCharacter) {
			if (start==-1) {  // Periodic
				return length>otherCharacter.length;
			}
			else {  // Nonperiodic
				return start>=otherCharacter.end;
			}
		}
		
		
		public final boolean sameOpen(Character otherCharacter) {
			return openStart==otherCharacter.openStart && openEnd==otherCharacter.openEnd;
		}
		
		
		public final boolean isOpen() { return openStart || openEnd; }
		
		
		public final void quantize(int quantum) {
			if (start!=-1) start=Math.round(start,quantum)*quantum;
			if (end!=-1) end=Math.round(end,quantum)*quantum; 
			length=Math.round(length,quantum)*quantum;
		}
		
		
		/**
		 * Not necessarily proper
		 */
		public final boolean isSuffixOf(Character otherCharacter, int threshold) {
			if (repeat!=otherCharacter.repeat || repeat==UNIQUE || orientation!=otherCharacter.orientation) return false;
			if (start==-1) return length<=otherCharacter.length+threshold;
			else return Math.abs(end,otherCharacter.end)<=threshold && start>=otherCharacter.start-threshold;
		}
		
		
		/**
		 * Not necessarily proper
		 */
		public final boolean isPrefixOf(Character otherCharacter, int threshold) {
			if (repeat!=otherCharacter.repeat || repeat==UNIQUE || orientation!=otherCharacter.orientation) return false;
			if (start==-1) return length<=otherCharacter.length+threshold;
			else return Math.abs(start,otherCharacter.start)<=threshold && end<=otherCharacter.end+threshold;
		}

		
		public final void reverseComplement() {
			if (repeat==UNIQUE) return;
			orientation=!orientation;
		}
	}
	
	
	
	
	
	
	
	
	// ----------------------------- ALIGNMENT PROCEDURES --------------------------------
	
	/**
	 * Basic filters on repeat alignments.
	 *
	 * Remark: non-periodic alignments that straddle other (periodic or non-periodic)
	 * alignments (possibly with the same readB) are kept, so their endpoints might create
	 * blocks when the read is recoded downstream. This is allowed, because some repeats 
	 * are seen to straddle in practice. Keeping all such alignments might create regions
	 * with many events and unclear signal: these are filtered out by $recodeRead()$.
	 *
	 * @param alignments all alignments of a given readA.
	 */
	private static final void cleanAlignments(int distanceThreshold) {
		int i, j;
		int startA, endA;
		AlignmentRow tmpAlignment;
		
		// Removing exact duplicates
		AlignmentRow.order=AlignmentRow.ORDER_READB_ORIENTATION_STARTA_ENDA_STARTB_ENDB;
		Arrays.sort(alignments,0,lastAlignment+1);
		j=0;
		for (i=1; i<=lastAlignment; i++) {
			if (!alignments[i].equals(alignments[j])) {
				j++;
				if (j!=i) {
					tmpAlignment=alignments[j];
					alignments[j]=alignments[i];
					alignments[i]=tmpAlignment;
				}
			}
		}
		lastAlignment=j;
		
		// If an alignment is strictly contained in another one in readA, we discard it
		// because it is less specific (regardless of what happens in readB).
		AlignmentRow.order=AlignmentRow.ORDER_STARTA;
		Arrays.sort(alignments,0,lastAlignment+1);
		for (i=0; i<=lastAlignment; i++) alignments[i].flag=true;
		for (i=0; i<=lastAlignment; i++) {
			if (!alignments[i].flag) continue;
			startA=alignments[i].startA; endA=alignments[i].endA;
			for (j=i+1; j<=lastAlignment; j++) {
				if (alignments[j].startA>=endA) break;
				if ( !Intervals.areApproximatelyIdentical(alignments[j].startA,alignments[j].endA,startA,endA) &&
					 Intervals.isApproximatelyContained(alignments[j].startA,alignments[j].endA,startA,endA)
				   ) alignments[j].flag=false;
			}
			for (j=i-1; j>=0; j--) {
				if (alignments[j].startA<startA-distanceThreshold) break;
				if ( !Intervals.areApproximatelyIdentical(alignments[j].startA,alignments[j].endA,startA,endA) &&
					 Intervals.isApproximatelyContained(alignments[j].startA,alignments[j].endA,startA,endA)
				   ) alignments[j].flag=false;
			}
		}
		j=-1;
		for (i=0; i<=lastAlignment; i++) {
			if (!alignments[i].flag) continue;
			j++;
			tmpAlignment=alignments[j];
			alignments[j]=alignments[i];
			alignments[i]=tmpAlignment;
		}
		lastAlignment=j;
	}
	
	
	/**
	 * A row of the alignments file
	 */
	private static class AlignmentRow implements Comparable {
		public static final int ORDER_READB_ORIENTATION_STARTA_ENDA_STARTB_ENDB = 0;
		public static final int ORDER_STARTA = 1;
		public static final int ORDER_READA_READB_ORIENTATION_STARTA_STARTB_ENDA_ENDB = 2;
		public static int order;
		
		public int readA, readB, startA, endA, startB, endB, diffs;
		public boolean orientation;
		
		/**
		 * Temporary space
		 */
		public boolean flag, flag2;
		
		public AlignmentRow() {
			this.readA=-1; this.startA=-1; this.endA=-1;
			this.readB=-1; this.startB=-1; this.endB=-1;
			this.orientation=false; this.diffs=-1;
		}
		
		public void set(int readA, int startA, int endA, int readB, int startB, int endB, boolean orientation, int diffs) {
			this.readA=readA; this.startA=startA; this.endA=endA;
			this.readB=readB; this.startB=startB; this.endB=endB;
			this.orientation=orientation; this.diffs=diffs;
		}
		
		public void copyFrom(AlignmentRow otherAlignment) {
			readA=otherAlignment.readA; startA=otherAlignment.startA; endA=otherAlignment.endA;
			readB=otherAlignment.readB; startB=otherAlignment.startB; endB=otherAlignment.endB;
			orientation=otherAlignment.orientation;
			diffs=Math.min(diffs,otherAlignment.diffs);
		}
		
		public double getRatio() {
			return ((double)(endA-startA))/(endB-startB);
		}
		
		public boolean equals(Object other) {
			AlignmentRow otherRow = (AlignmentRow)other;
			return otherRow.readA==readA && otherRow.readB==readB && 
				   otherRow.startA==startA && otherRow.endA==endA && 
				   otherRow.startB==startB && otherRow.endB==endB && 
				   otherRow.orientation==orientation;
		}
		
		public int compareTo(Object other) {
			AlignmentRow otherAlignment = (AlignmentRow)other;
			
			if (order==ORDER_READB_ORIENTATION_STARTA_ENDA_STARTB_ENDB) {
				if (readB<otherAlignment.readB) return -1;
				else if (readB>otherAlignment.readB) return 1;
				if (orientation && !otherAlignment.orientation) return -1;
				else if (!orientation && otherAlignment.orientation) return 1;
				if (startA<otherAlignment.startA) return -1;
				else if (startA>otherAlignment.startA) return 1;
				if (endA<otherAlignment.endA) return -1;
				else if (endA>otherAlignment.endA) return 1;
				if (startB<otherAlignment.startB) return -1;
				else if (startB>otherAlignment.startB) return 1;
				if (endB<otherAlignment.endB) return -1;
				else if (endB>otherAlignment.endB) return 1;
			}
			else if (order==ORDER_STARTA) {
				if (startA<otherAlignment.startA) return -1;
				else if (startA>otherAlignment.startA) return 1;
			}
			else if (order==ORDER_READA_READB_ORIENTATION_STARTA_STARTB_ENDA_ENDB) {
				if (readA<otherAlignment.readA) return -1;
				else if (readA>otherAlignment.readA) return 1;
				if (readB<otherAlignment.readB) return -1;
				else if (readB>otherAlignment.readB) return 1;
				if (orientation && !otherAlignment.orientation) return -1;
				else if (!orientation && otherAlignment.orientation) return 1;
				if (startA<otherAlignment.startA) return -1;
				else if (startA>otherAlignment.startA) return 1;
				if (startB<otherAlignment.startB) return -1;
				else if (startB>otherAlignment.startB) return 1;
				if (endA<otherAlignment.endA) return -1;
				else if (endA>otherAlignment.endA) return 1;
				if (endB<otherAlignment.endB) return -1;
				else if (endB>otherAlignment.endB) return 1;
			}
			return 0;
		}
		
		public String toString() {
			return readA+"["+startA+".."+endA+"] x "+readB+"["+startB+".."+endB+"] :: "+orientation;
		}	
	}

}