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
	 * Temporary space
	 */
	private static String[][] blocks;
	private static int[] lastInBlock, lastInBlock_int;
	private static int[] boundaries;
	private static int[][] intBlocks;
	private static boolean[] isBlockUnique;
	private static Pair[] tandemIntervals;
	private static Kmer[] kmerPool;
	private static int lastKmerPool;
	
	
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
	 * sequence$.
	 */
	public static final void recodeRead(int distanceThreshold) {
		final int CLUSTERING_DISTANCE = distanceThreshold;
		final int MAX_DENSE_LENGTH = (CLUSTERING_DISTANCE)<<1;  // Arbitrary
		int i, j, k;
		int lastPeriodicInterval, currentStart, currentEnd, first, firstZero;
		int firstJForNextI, inPeriodic;
		int startA, endA, newLength, component, nComponents;
		final int lengthA = Reads.getReadLength(alignments[0].readA);
		
		if (periodicIntervals.length<(lastAlignment+1)<<1) periodicIntervals = new int[(lastAlignment+1)<<1];
		if (points.length<(lastAlignment+1)<<1) points = new int[(lastAlignment+1)<<1];
		AlignmentRow.order=AlignmentRow.ORDER_STARTA;
		if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
		

boolean fabio = false;
for (i=0; i<=lastAlignment; i++) {
	if (alignments[i].readA==779 && alignments[i].readB==1016) {
		fabio=true;
		break;
	}
}
if (fabio) {
	System.err.println("VITTU> 1  alignments:");
	for (i=0; i<=lastAlignment; i++) System.err.println(alignments[i]);
}
		
		
		
		// Building maximal periodic intervals
		lastPeriodicInterval=-1; currentStart=-1; currentEnd=-1;
		for (i=0; i<=lastAlignment; i++) {
			if (!isPeriodic[alignments[i].readB]) continue;
			if (currentStart==-1) {
				currentStart=alignments[i].startA; currentEnd=alignments[i].endA;
				continue;
			}
			if (alignments[i].startA>=currentEnd-distanceThreshold) {
				periodicIntervals[++lastPeriodicInterval]=currentStart;
				periodicIntervals[++lastPeriodicInterval]=currentEnd;
				currentStart=alignments[i].startA; currentEnd=alignments[i].endA;
				continue;
			}
			currentEnd=Math.max(currentEnd,alignments[i].endA);
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

		
if (fabio) {
	System.err.println("VITTU> 2  points:");
	for (i=0; i<=lastPoint; i++) System.err.println(points[i]);
}				
		
		
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
			if ( Intervals.areApproximatelyIdentical(alignments[j].startA,alignments[j].endA,startA,endA) || 
				 (isPeriodic[alignments[j].readB] && Intervals.isApproximatelyContained(alignments[j].startA,alignments[j].endA,startA,endA))
			   ) newBlock.addCharacter(alignments[j],distanceThreshold,alignments[j].startA,alignments[j].endA,tmpCharacter);
			else if (!isPeriodic[alignments[j].readB] && Intervals.isApproximatelyContained(startA,endA,alignments[j].startA,alignments[j].endA)) newBlock.addCharacter(alignments[j],distanceThreshold,startA,endA,tmpCharacter);
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
	public static final void compactInstances(int distanceThreshold, int lengthThreshold) {
		final int QUANTUM = IO.quantum;
		int i, j, k;
		int first;
		Character tmpChar;
		
		System.err.println("Quantizing characters... ");
		for (i=0; i<=lastAlphabet; i++) alphabet[i].quantize(QUANTUM);
		Arrays.sort(alphabet,0,lastAlphabet+1);
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
				if (alphabet[i].implies(alphabet[j],distanceThreshold)) alphabet[j].flag=0;
			}
			for (j=i-1; j>=0; j--) {
				if (alphabet[j].repeat!=alphabet[i].repeat || alphabet[j].orientation!=alphabet[i].orientation) break;
				if (alphabet[j].flag==0) continue;
				if (alphabet[i].implies(alphabet[j],distanceThreshold)) alphabet[j].flag=0;
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
	 * the alphabet. Reverse-complement pairs have symmetrical open status.
	 */
	public static final void closeAlphabetByRC() {
		final int QUANTUM = 1000;  // Arbitrary
		boolean openStart, openEnd;
		int i, j, k;
		int from, newFrom, repeat, start, end, length, found, nAdded, lastPeriodicNew;
		Character[] newAlphabet;
		
		// Unique
		newAlphabet = new Character[lastAlphabet+1];
		System.arraycopy(alphabet,0,newAlphabet,0,lastUnique+1);
		for (i=lastUnique+1; i<=lastAlphabet; i++) alphabet[i].flag=0;
		System.err.println("Adding reverse-complement characters... ");
		nAdded=0;
		
		// Periodic
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
			length=alphabet[i].length;
			if (alphabet[i].orientation) {
				found=-1;
				for (j=i+1; j<=lastPeriodic; j++) {
					if (alphabet[j].repeat!=repeat) break;
					if (!alphabet[j].orientation && alphabet[j].length==length && alphabet[j].openStart==openEnd && alphabet[j].openEnd==openStart) {
						found=j;
						break;
					}
				}
				if (found==-1) {
					k++;
					if (k==newAlphabet.length) {
						Character[] newArray = new Character[newAlphabet.length+QUANTUM];
						System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
						newAlphabet=newArray;
					}
					newAlphabet[k] = new Character(repeat,false,-1,-1,length,openEnd,openStart);
					nAdded++;
				}
				else alphabet[found].flag=1;
			}
			else {
				found=-1;
				for (j=i-1; j>lastUnique; j--) {
					if (alphabet[j].repeat!=repeat) break;
					if (alphabet[j].orientation && alphabet[j].length==length && alphabet[j].openStart==openEnd && alphabet[j].openEnd==openStart) {
						found=j;
						break;
					}
				}
				if (found==-1) {
					k++;
					if (k==newAlphabet.length) {
						Character[] newArray = new Character[newAlphabet.length+QUANTUM];
						System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
						newAlphabet=newArray;
					}
					newAlphabet[k] = new Character(repeat,true,-1,-1,length,openEnd,openStart);
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
			start=alphabet[i].start; end=alphabet[i].end;
			if (alphabet[i].orientation) {
				found=-1;
				for (j=i+1; j<=lastAlphabet; j++) {
					if (alphabet[j].repeat!=repeat) break;
					if (!alphabet[j].orientation && alphabet[j].start==start && alphabet[j].end==end && alphabet[j].openStart==openEnd && alphabet[j].openEnd==openStart) {
						found=j;
						break;
					}
				}
				if (found==-1) {
					k++;
					if (k==newAlphabet.length) {
						Character[] newArray = new Character[newAlphabet.length+QUANTUM];
						System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
						newAlphabet=newArray;
					}
					newAlphabet[k] = new Character(repeat,false,start,end,0,openEnd,openStart);
					nAdded++;
				}
				else alphabet[found].flag=1;
			}
			else {
				found=-1;
				for (j=i-1; j>lastUnique; j--) {
					if (alphabet[j].repeat!=repeat) break;
					if (alphabet[j].orientation && alphabet[j].start==start && alphabet[j].end==end && alphabet[j].openStart==openEnd && alphabet[j].openEnd==openStart) {
						found=j;
						break;
					}
				}
				if (found==-1) {
					k++;
					if (k==newAlphabet.length) {
						Character[] newArray = new Character[newAlphabet.length+QUANTUM];
						System.arraycopy(newAlphabet,0,newArray,0,newAlphabet.length);
						newAlphabet=newArray;
					}
					newAlphabet[k] = new Character(repeat,true,start,end,0,openEnd,openStart);
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
				while (k<=lastPeriodic && alphabet[k].repeat==repeat && alphabet[k].orientation==orientation && !alphabet[k].implies(character,-1)) k++;
				i=k;
				if (alphabet[i].repeat!=repeat || alphabet[i].orientation!=orientation || !alphabet[i].implies(character,-1)) {
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
					if (alphabet[j].implies(character,quantum)) {
					   found=true;
					   last=appendToStack(j,last);
					}
				}
				for (j=i-1; j>lastPeriodic; j--) {
					if (alphabet[j].repeat!=repeat || alphabet[j].orientation!=orientation) break;
					if (alphabet[j].implies(character,quantum)) {
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
							if (alphabet[k].implies(alphabet[c],-1)) characterCount[k]++;
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
	 */
	private static final void loadIntBlocks(int nBlocks) {
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
						repeat=alphabet[value].repeat; orientation=alphabet[value].orientation;
						for (k=value; k<=lastPeriodic; k++) {
							if (alphabet[k].repeat!=repeat || alphabet[k].orientation!=orientation) break;
							if (alphabet[k].implies(alphabet[value],-1)) { found=true; nElements++; }
						}
					}
				}
				if (!found) {
					System.err.println("loadIntBlocks> ERROR: character ID "+blocks[i][j]+" in a translated read is not implied by any character in the alphabet.");
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
						repeat=alphabet[value].repeat; orientation=alphabet[value].orientation;
						for (c=value; c<=lastPeriodic; c++) {
							if (alphabet[c].repeat!=repeat || alphabet[c].orientation!=orientation) break;
							if (alphabet[c].implies(alphabet[value],-1)) intBlocks[i][++k]=c;
						}
					}
				}
			}
			lastInBlock_int[i]=k; isBlockUnique[i]=unique;
		}
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
		if (read2boundaries.length()==0) nBoundaries=0;
		else if (read2boundaries.indexOf(SEPARATOR_MINOR+"")>=0) {
			tokens=read2boundaries.split(SEPARATOR_MINOR+"");
			nBoundaries=tokens.length;
			if (boundaries==null || boundaries.length<nBoundaries) boundaries = new int[nBoundaries];
			for (i=0; i<nBoundaries; i++) boundaries[i]=Integer.parseInt(tokens[i]);
		}
		else {
			nBoundaries=1;
			if (boundaries==null || boundaries.length<nBoundaries) boundaries = new int[nBoundaries];
			boundaries[0]=Integer.parseInt(read2boundaries);
		}
		nBlocks=loadBlocks(read2characters);
		removeRareCharacters(nBlocks,minCount,lastUnique,lastPeriodic,lastAlphabet,keepPeriodic);
		first=-1;
		for (i=0; i<nBlocks; i++) {
			if (lastInBlock[i]==-1) isUnique=true;
			else if (lastInBlock[i]>0) isUnique=false;
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
								if (alphabetCount[h]>=minCount && alphabet[h].implies(alphabet[c],-1)) {
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
			for (i=lastUnique+1; i<=lastPeriodic; i++) old2new[i-lastUnique+1]=i-lastUnique+1;
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
		if (read2boundaries_old.length()==0) nBoundaries=0;
		else if (read2boundaries_old.indexOf(SEPARATOR_MINOR+"")>=0) {
			tokens=read2boundaries_old.split(SEPARATOR_MINOR+"");
			nBoundaries=tokens.length;
			if (boundaries==null || boundaries.length<nBoundaries) boundaries = new int[nBoundaries];
			for (i=0; i<nBoundaries; i++) boundaries[i]=Integer.parseInt(tokens[i]);
		}
		else {
			nBoundaries=1;
			if (boundaries==null || boundaries.length<nBoundaries) boundaries = new int[nBoundaries];
			boundaries[0]=Integer.parseInt(read2boundaries_old);
		}
		nBlocks=loadBlocks(read2characters_old);
		alphabet=oldAlphabet; lastUnique=lastUnique_old; lastPeriodic=lastPeriodic_old; lastAlphabet=lastAlphabet_old;
		removeRareCharacters(nBlocks,minCount,lastUnique_old,lastPeriodic_old,lastAlphabet_old,keepPeriodic);
		alphabet=newAlphabet; lastUnique=lastUnique_new; lastPeriodic=lastPeriodic_new; lastAlphabet=lastAlphabet_new;
		first=-1; nAppendedBlocks=0;
		for (i=0; i<nBlocks; i++) {
			if (lastInBlock[i]==-1) isUnique=true;
			else if (lastInBlock[i]>0) isUnique=false;
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
	 * Remark: in the latter case above, 1-mers in $oldKmers$ are not used if they occur
	 * in an endblock of the read, or in a block with multiple characters (for $k>1$ we do
	 * not do this, since we assume that longer contexts are enough to disambiguate).
     *
	 * Remark: 1-mers collected by this procedure might have a different (and even 
	 * smaller) count than the one produced by the $getCharacterHistogram()$ pipeline, 
	 * because of the constraints and of canonization. E.g. a 1-mer might be rare in this 
	 * procedure but frequent in $getCharacterHistogram()$, if the latter counted all the
	 * occurrences of the character but the former counts just closed occurrences.
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
	 * @param haplotypeCoverage coverage of one haplotype;
	 * @param tmpKmer temporary space;
	 * @param tmpArray2 temporary space, of size at least k;
	 * @param tmpArray3 temporary space, of size at least 2k;
	 * @param tmpMap temporary hashmap, used only if $newKmers$ is not null.
	 */
	public static final int getKmers(String str, int k, int uniqueMode, int multiMode, HashMap<Kmer,Kmer> newKmers, HashMap<Kmer,Kmer> oldKmers, int[] avoidedIntervals, int lastAvoidedInterval, int haplotypeCoverage, Kmer tmpKmer, int[] tmpArray2, int[] tmpArray3, HashMap<Kmer,Kmer> tmpMap) {
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
		loadBlocks(str); loadIntBlocks(nBlocks);
		sum=0;
		for (i=0; i<nBlocks; i++) sum+=lastInBlock_int[i]+1;
		if (stack==null || stack.length<sum*3) stack = new int[sum*3];
		
		// Processing every k-mer in the read
		if (newKmers==null) {
			j=0;
			for (i=0; i<=nBlocks-k; i++) {
				while (j<lastAvoidedInterval && avoidedIntervals[j]<i) j+=3;
				if (j<lastAvoidedInterval && avoidedIntervals[j]+avoidedIntervals[j+1]-1<=i+k-1) continue;
				if (!isValidWindow(i,k,nBlocks,uniqueMode,multiMode)) continue;
				nKmers=lastInBlock_int[i]+1;
				for (p=i+1; p<=i+k-1; p++) nKmers*=lastInBlock_int[p]+1;
				if (nKmers<0 || nKmers>MAX_KMERS_TO_ENUMERATE) continue;
				nHaplotypes=getKmers_impl(i,k,null,oldKmers,haplotypeCoverage,tmpKmer,stack,tmpArray2,tmpArray3);
				if (nHaplotypes!=-1 && (k>1?true:i>0&&i<nBlocks-1&&lastInBlock_int[i]==0)) { avoidedIntervals[++out]=i; avoidedIntervals[++out]=k; avoidedIntervals[++out]=nHaplotypes; }
			}
		}
		else {
			tmpMap.clear(); lastKmerPool=-1;
			j=0;
			for (i=0; i<=nBlocks-k; i++) {
				while (j<lastAvoidedInterval && avoidedIntervals[j]<i) j+=3;
				if (j<lastAvoidedInterval && avoidedIntervals[j]+avoidedIntervals[j+1]-1<=i+k-1) continue;
				if (!isValidWindow(i,k,nBlocks,uniqueMode,multiMode)) continue;
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
	 * Tells whether window $blocks[first..first+k-1]$ satisfies the following:
	 *
	 * @param uniqueMode blocks with unique characters:
	 * 0: are allowed; 
	 * 1: are allowed everywhere, except in the first and last block of the window;
	 * 2: are not allowed;
	 * @param multiMode blocks with multiple characters:
	 * 0: are allowed; 
	 * 1: are not allowed if they are the first/last one (endblocks are more likely to
	 *    match several characters because they contain just a fraction of a character);
	 * 2: are not allowed.
	 */
	private static final boolean isValidWindow(int first, int k, int nBlocks, int uniqueMode, int multiMode) {
		int i;
		
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
				key.set(tmpArray2,0,k,true); key.canonize(k,tmpArray3);
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
	public static final void fixEndBlocks(String str, int k, HashMap<Kmer,Kmer> kmers, boolean tightMode, Kmer context, int[] tmpArray1, int[] tmpArray2, int[] tmpArray3, BufferedWriter bw, int[] out, int[][] ambiguityHistogram) throws IOException {
		int i, j, c, d, p, q;
		int nBlocks, sum, fixedFirst, fixedLast;
		
		if (str.length()==0) { bw.newLine(); return; }
		
		// Loading blocks
		nBlocks=loadBlocks(str); loadIntBlocks(nBlocks);
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
	 * $inputFile$, as a sequence of pairs $firstBlock,lastBlock$.
	 *
	 * @return the total number of tandem intervals found.
	 */
	public static final long getTandemIntervals(String inputFile, String outputFile) throws IOException {
		int i;
		int nBlocks, lastTandem;
		long out;
		String str;
		BufferedReader br;
		BufferedWriter bw;
		
		br = new BufferedReader(new FileReader(inputFile));
		bw = new BufferedWriter(new FileWriter(outputFile));
		str=br.readLine(); out=0;
		while (str!=null) {
			if (str.length()==0 || str.indexOf(SEPARATOR_MAJOR+"")<0) {
				bw.newLine();
				str=br.readLine();
				continue;
			}
			nBlocks=loadBlocks(str); loadIntBlocks(nBlocks);
			lastTandem=getTandemIntervals_impl(nBlocks);
			if (lastTandem==-1) {
				bw.newLine();
				str=br.readLine();
				continue;
			}
			out+=lastTandem+1;
			bw.write(tandemIntervals[0].from+""+SEPARATOR_MINOR+""+tandemIntervals[0].to);
			for (i=1; i<=lastTandem; i++) bw.write(SEPARATOR_MINOR+""+tandemIntervals[i].from+""+SEPARATOR_MINOR+""+tandemIntervals[i].to);
			bw.newLine();
			str=br.readLine();
		}
		br.close(); bw.close();
		return out;
	}
	
	
	/**
	 * Stores in global variable $tandemIntervals$ the set of all maximal tandem intervals
	 * of $intBlocks$, where a tandem interval is a substring such that all its blocks
	 * have a character in common. Tandem intervals that intersect are merged.
	 *
	 * Remark: this definition of tandem is likely too simple for real data, since the
	 * aligner might fail to align some repeat to some tandem unit, or the alignment might
	 * be a bit off and give rise to a different character in our alphabet. 
	 *
	 * Remark: the procedure uses global arrays $stack,stack2$.
	 *
	 * @return the last element in $tandemIntervals$.
	 */
	private static final int getTandemIntervals_impl(int nBlocks) {
		int i, j;
		int max, last1, last2, lastInterval, tandemLength;
		Pair tmpInterval;
		int[] tmpArray;
		
		// Allocating memory
		if (tandemIntervals==null) tandemIntervals = new Pair[100];  // Arbitrary
		max=0;
		for (i=0; i<nBlocks; i++) max=Math.max(max,lastInBlock_int[i]+1);
		if (stack==null || stack.length<max) stack = new int[max];
		if (stack2==null || stack2.length<max) stack2 = new int[max];
		
		// Collecting intervals
		lastInterval=-1;
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
				lastInterval++;
				if (lastInterval>=tandemIntervals.length) {
					Pair[] newArray = new Pair[tandemIntervals.length<<1];
					System.arraycopy(tandemIntervals,0,newArray,0,tandemIntervals.length);
					tandemIntervals=newArray;
				}
				if (tandemIntervals[lastInterval]==null) tandemIntervals[lastInterval] = new Pair(i,i+tandemLength-1);
				else tandemIntervals[lastInterval].set(i,i+tandemLength-1);
			}
		}
		if (lastInterval<=0) return lastInterval;
		
		// Merging overlapping intervals
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
	 * @param loadIsBlockUnique loads the bitvector $isBlockUnique_all$ for every read in 
	 * $translated$;
	 * @param loadTranslation loads in $translation_all$ the translation of every read 
	 * in $translated$.
	 */
	public static final void loadAllBoundaries(String translatedFile, boolean loadIsBlockUnique, boolean loadTranslation, String boundariesFile) throws IOException {
		int i, j, r;
		int read, cell, mask, block, nBlocks, nBytes, nReads, nBoundaries;
		String str;
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
		if (loadIsBlockUnique) {
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
				loadIntBlocks(nBlocks);
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
		if (loadTranslation) {
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
				loadIntBlocks(nBlocks);
				translation_all[r] = new int[nBlocks][0];
				for (i=0; i<nBlocks; i++) {
					translation_all[r][i] = new int[lastInBlock_int[i]+1];
					System.arraycopy(intBlocks[i],0,translation_all[r][i],0,lastInBlock_int[i]+1);
				}
				str=br.readLine();
			}
			br.close();
		}
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
	public static final void filterAlignments_loose(String alignmentsFile, String outputFile, int minIntersection_nonrepetitive, long[][] out) throws IOException {
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
				isRepetitive=inRedRegion(readA,Alignments.startA,Alignments.endA,lastTranslated,lastBlueInterval<=blueIntervals_last&&blueIntervals_reads[lastBlueInterval]==readA?lastBlueInterval:-1,-1,Reads.getReadLength(readA),minIntersection_nonrepetitive);
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
			isRepetitive=inRedRegion(readB,Alignments.startB,Alignments.endB,p,-2,lastBlueInterval,Reads.getReadLength(readB),minIntersection_nonrepetitive);
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
	 * -2: contains a sequence of repeat characters that likely occurs <=H times in the 
	 *     genome;
	 *  0: none of the above is true.
	 */
	private static final int inRedRegion(int readID, int intervalStart, int intervalEnd, int boundariesAllID, int blueIntervalsID, int blueIntervalsStart, int readLength, int minIntersection_nonrepetitive) {
		int i, j;
		int mask, cell, start, end, blockStart, blockEnd, firstBlock, lastBlock;
		final int nBlocks = boundaries_all[boundariesAllID].length+1;
		final int nBytes = Math.ceil(nBlocks,8);
		
		// Checking the nonrepetitive blocks of the read, if any.
		cell=isBlockUnique_all[boundariesAllID][0]; 
		if (cell!=0) {
			mask=1; blockStart=0; blockEnd=boundaries_all[boundariesAllID][0];
			if ( (cell&mask)!=0 && 
				 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
				      Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive
				   )
				 )
			   ) return -1;
			mask<<=1;
			for (j=1; j<8; j++) {
				if (j==nBlocks-1) break;
				blockStart=boundaries_all[boundariesAllID][j-1];
				blockEnd=boundaries_all[boundariesAllID][j];
				if ( (cell&mask)!=0 && 
					 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) || 
					   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
					   ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
						  Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive
					   )
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
					 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
					   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
					   ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
						  Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive
					   )
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
				 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
					  Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive
				   )
				 )
			   ) return -1;
		}
		
		// Checking the repetitive blocks of the read, if any.
		if (blueIntervalsID==-2) blueIntervalsID=readInArray(readID,blueIntervals_reads,blueIntervals_reads.length-1,blueIntervalsStart);
		if (blueIntervalsID!=-1) {
			for (i=0; i<blueIntervals[blueIntervalsID].length; i+=3) {
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
	public static final void filterAlignments_tight(String alignmentsFile, String outputFile, boolean mode, boolean suffixPrefixMode, int minIntersection_nonrepetitive, int minIntersection_repetitive, long[][] out) throws IOException {
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
				q=inBlueRegion(readA,startA,endA,lastTranslated,(lastBlueInterval<=blueIntervals_last&&blueIntervals_reads[lastBlueInterval]==readA)?lastBlueInterval:-1,-1,lengthA,minIntersection_nonrepetitive,minIntersection_repetitive);
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
			q=inBlueRegion(readB,startB,endB,p,-2,lastBlueInterval,Reads.getReadLength(readB),minIntersection_nonrepetitive,minIntersection_repetitive);
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
	 * @return interval $readID[intervalStart..intervalEnd]$:
	 * 0: fully belongs to a non-repetitive region;
	 * 1: straddles a non-repetitive region;
	 * 2: contains a sequence of repeat characters that is likely to occur <=H times in 
	 *    the genome, where H is the number of haplotypes;
	 * 3: straddles, but does not fully contain, a sequence in point (2), on the left side
	 *    of the interval;
	 * 4: straddles, but does not fully contain, a sequence in point (2), on the right
	 *    side of the interval;
	 * 5: both (3) and (4) are true;
	 * -1: none of the above is true.
	 */
	private static final int inBlueRegion(int readID, int intervalStart, int intervalEnd, int boundariesAllID, int blueIntervalsID, int blueIntervalsStart, int readLength, int minIntersection_nonrepetitive, int minIntersection_repetitive) {
		boolean straddlesLeft, straddlesRight;
		int i, j;
		int start, end, blockStart, blockEnd, firstBlock, lastBlock;
		final int nBlocks = boundaries_all[boundariesAllID].length+1;
		final int nBytes = Math.ceil(nBlocks,8);
		
		// Checking the nonrepetitive blocks of the read, if any.
		blockStart=0; blockEnd=boundaries_all[boundariesAllID][0];
		if (translation_all[boundariesAllID][0].length==1 && (translation_all[boundariesAllID][0][0]<=lastUnique || translation_all[boundariesAllID][0][0]==lastAlphabet+1)) {
			if ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
			     Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd)
			   ) return 0;
			else if ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
			           Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive
			   		) return 1;
		}
		for (i=1; i<nBlocks-1; i++) {
			blockStart=boundaries_all[boundariesAllID][i-1];
			blockEnd=boundaries_all[boundariesAllID][i];
			if (translation_all[boundariesAllID][i].length==1 && (translation_all[boundariesAllID][i][0]<=lastUnique || translation_all[boundariesAllID][i][0]==lastAlphabet+1)) {
				if ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
				     Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd)
				   ) return 0;
				else if ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
				      	   Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive
				   		) return 1;
			}
		}
		blockStart=boundaries_all[boundariesAllID][nBlocks-2];
		blockEnd=readLength-1;
		if (translation_all[boundariesAllID][nBlocks-1].length==1 && (translation_all[boundariesAllID][nBlocks-1][0]<=lastUnique || translation_all[boundariesAllID][nBlocks-1][0]==lastAlphabet+1)) {
			if ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
			     Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd)
			   ) return 0;
			else if ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
			           Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection_nonrepetitive
			        ) return 1;
		}
		
		// Checking the repetitive blocks of the read, if any.
		if (blueIntervalsID==-2) blueIntervalsID=readInArray(readID,blueIntervals_reads,blueIntervals_reads.length-1,blueIntervalsStart);
		if (blueIntervalsID!=-1) {
			straddlesLeft=false; straddlesRight=false;
			for (i=0; i<blueIntervals[blueIntervalsID].length; i+=3) {
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
	public static final boolean containsUnique(String str) {
		final int nBlocks = loadBlocks(str); 
		loadIntBlocks(nBlocks);
		for (int i=0; i<nBlocks; i++) {
			if (isBlockUnique[i]) return true;
		}
		return false;
	}
	
	
	/**
	 * Let a tandem be a maximal sequence of adjacent blocks with characters in common.
	 * The procedure writes a zero for every alignment that should not be trusted because
	 * it either (1) is strictly contained in a tandem, or (2) is identical to a tandem, 
	 * but the tandem is not a blue interval, or (3) straddles or contains a tandem, and 
	 * all blue intervals in the alignment fall inside a tandem.
	 *
	 * Remark: every block inside a tandem might be tagged with several characters, and it
	 * might happen that a sequence of such characters occurs just once inside the tandem,
	 * and that it is considered blue based on its total frequency in the read set. In 
	 * practice such a sequence is likely noise.
	 *
	 * @param bothReads discards an alignment if the conditions above hold on both reads
	 * (TRUE) or on just one read (FALSE);
	 * @param out output array containing the number of alignments for each type (columns)
	 * specified in $Alignments.readAlignmentFile_getType()$; row 0: all alignments in 
	 * input; row 1: all alignments kept in output.
	 */
	public static final void filterAlignments_tandem(String alignmentsFile, boolean bothReads, String outputFile, long[][] out) throws IOException {
		final int IDENTITY_THRESHOLD = IO.quantum;
		boolean found, containedA, identicalA, containedB, identicalB;
		int i, p;
		int length, row, type, nBlocks, readA, readB, startA, endA, startB, endB;
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
					if (tandems[readA][i+1]==nBlocks-1) lastTandemPosition=Reads.getReadLength(readA)-1;
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
					if (containedA || (identicalA && !filterAlignments_tandem_isBlue(blueIntervalA,firstTandemBlockA,lastTandemBlockA)) || (!identicalA && filterAlignments_tandem_allContainedBlueInTandem(readA,lastTranslated,blueIntervalA,startA,endA,nBlocks))) {
						bw.write("0\n"); str=br.readLine(); row++;
						continue;
					}
				}
				else {
					if ((identicalA && filterAlignments_tandem_isBlue(blueIntervalA,firstTandemBlockA,lastTandemBlockA)) || (!identicalA && !containedA && !filterAlignments_tandem_allContainedBlueInTandem(readA,lastTranslated,blueIntervalA,startA,endA,nBlocks))) {
if (readA==940 && readB==1108) System.err.println("filterAlignments_tandem> 2  WHAT??? WE KEPT ALIGNMENT "+str);
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
			p=readInArray(readB,translated,nTranslated-1,lastTranslated);
			if (p>=0) {
				blueIntervalB=readInArray(readB,blueIntervals_reads,blueIntervals_reads.length-1,lastBlueInterval);
				startB=Alignments.startB; endB=Alignments.endB;
				nBlocks=boundaries_all[p].length+1;
				containedB=false; identicalB=false; firstTandemBlockB=-1; lastTandemBlockB=-1;
				for (i=0; i<lastTandem[readB]; i+=2) {
					if (tandems[readB][i]==0) firstTandemPosition=0;
					else firstTandemPosition=boundaries_all[p][tandems[readB][i]-1];
					if (tandems[readB][i+1]==nBlocks-1) lastTandemPosition=Reads.getReadLength(readB)-1;
					else lastTandemPosition=boundaries_all[p][tandems[readB][i+1]];
if (readA==940 && readB==1108) System.err.println("filterAlignments_tandem> 3.5  ["+startB+".."+endB+"] :: ["+firstTandemPosition+".."+lastTandemPosition+"]");				
					if (Intervals.areApproximatelyIdentical(startB,endB,firstTandemPosition,lastTandemPosition)) identicalB=true;
					else if (Intervals.isApproximatelyContained(startB,endB,firstTandemPosition,lastTandemPosition)) containedB=true;
					if (identicalB || containedB) {
						firstTandemBlockB=tandems[readB][i];
						lastTandemBlockB=tandems[readB][i+1];
						break;
					}
				}
				if (!bothReads) {
					if (containedB || (identicalB && !filterAlignments_tandem_isBlue(blueIntervalB,firstTandemBlockB,lastTandemBlockB)) || (!identicalB && filterAlignments_tandem_allContainedBlueInTandem(readB,p,blueIntervalB,startB,endB,nBlocks))) {
						bw.write("0\n"); str=br.readLine(); row++;
						continue;
					}
				}
				else {
					if ((identicalB && filterAlignments_tandem_isBlue(blueIntervalB,firstTandemBlockB,lastTandemBlockB)) || (!identicalB && !containedB && !filterAlignments_tandem_allContainedBlueInTandem(readB,p,blueIntervalB,startB,endB,nBlocks))) {
if (readA==940 && readB==1108) {
	System.err.println("filterAlignments_tandem> 4.1  WHAT??? WE KEPT ALIGNMENT "+str);
	System.err.println("tandems of readA:");
	for (int x=0; x<=lastTandem[readA]; x++) System.err.print(tandems[readA][x]+",");
	System.err.println();
	System.err.println("tandems of readB:");
	for (int x=0; x<=lastTandem[readB]; x++) System.err.print(tandems[readB][x]+",");
	System.err.println();
//							System.err.println("containedA="+containedA+" identicalA="+identicalA+" firstTandemBlockB="+firstTandemBlockB+" lastTandemBlockB="+lastTandemBlockB);
}
						out[1][type]++;
						bw.write("1\n"); str=br.readLine(); row++;
						continue;
					}
				}
			}
			else if (bothReads) {
if (readA==940 && readB==1108) System.err.println("filterAlignments_tandem> 3  WHAT??? WE KEPT ALIGNMENT "+str);
				out[1][type]++;
				bw.write("1\n"); str=br.readLine(); row++;
				continue;
			}
			if (bothReads) bw.write("0\n");
			else {
				bw.write("1\n");
if (readA==940 && readB==1108) System.err.println("filterAlignments_tandem> 5  WHAT??? WE KEPT ALIGNMENT "+str);		
			}
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
	 */
	private static final boolean filterAlignments_tandem_allContainedBlueInTandem(int readID, int boundariesAllID, int blueIntervalsID, int intervalStart, int intervalEnd, int nBlocks) {
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
					if (tandems[readID][j]<firstBlock || tandems[readID][j+1]>lastBlock) {
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
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		int[] leftEnd1, leftEnd2, rightEnd1, rightEnd2;
		int[] tmpArray, tmpArray1, tmpArray2;
		String[] tokens;
		
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
					loadIntBlocks(nBlocks);
					if (lastInBlock_int[nBlocks-1]>=leftEnd1.length) leftEnd1 = new int[lastInBlock_int[nBlocks-1]+1];
					System.arraycopy(intBlocks[nBlocks-1],0,leftEnd1,0,lastInBlock_int[nBlocks-1]+1);
					leftEnd1_last=lastInBlock_int[nBlocks-1];
					leftEnd1_str.delete(0,leftEnd1_str.length());
					leftEnd1_str.append(str1.lastIndexOf(SEPARATOR_MAJOR+"")+1);
					nBlocks=loadBlocks(str2); loadIntBlocks(nBlocks);
					if (lastInBlock_int[nBlocks-1]>=leftEnd2.length) leftEnd2 = new int[lastInBlock_int[nBlocks-1]+1];
					System.arraycopy(intBlocks[nBlocks-1],0,leftEnd2,0,lastInBlock_int[nBlocks-1]+1);
					leftEnd2_last=lastInBlock_int[nBlocks-1];
					leftEnd2_str.delete(0,leftEnd2_str.length());
					leftEnd2_str.append(str2.lastIndexOf(SEPARATOR_MAJOR+"")+1);
					leftLength=Reads.readLengths[r]-Integer.parseInt(str3.substring(str3.lastIndexOf(",")+1).trim())-1;
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
			loadIntBlocks(nBlocks);
			if (lastInBlock_int[0]>=rightEnd1.length) rightEnd1 = new int[lastInBlock_int[0]+1];
			System.arraycopy(intBlocks[0],0,rightEnd1,0,lastInBlock_int[0]+1);
			rightEnd1_last=lastInBlock_int[0];
			if (lastInBlock_int[nBlocks-1]>=tmpArray1.length) tmpArray1 = new int[lastInBlock_int[nBlocks-1]+1];
			System.arraycopy(intBlocks[nBlocks-1],0,tmpArray1,0,lastInBlock_int[nBlocks-1]+1);
			last1=lastInBlock_int[nBlocks-1];
			nBlocks=loadBlocks(str2); loadIntBlocks(nBlocks);
			if (lastInBlock_int[0]>=rightEnd2.length) rightEnd2 = new int[lastInBlock_int[0]+1];
			System.arraycopy(intBlocks[0],0,rightEnd2,0,lastInBlock_int[0]+1);
			rightEnd2_last=lastInBlock_int[0];
			if (lastInBlock_int[nBlocks-1]>=tmpArray2.length) tmpArray2 = new int[lastInBlock_int[nBlocks-1]+1];
			System.arraycopy(intBlocks[nBlocks-1],0,tmpArray2,0,lastInBlock_int[nBlocks-1]+1);
			rightUnique=isBlockUnique[0];
			last2=lastInBlock_int[nBlocks-1];
			if ((leftEnd1_last>0 && leftEnd2_last==0) || (rightEnd1_last>0 && rightEnd2_last==0 && leftEnd2_last>=0)) {
				nAttempted++;
				success=breakReads_checkDisambiguation_impl(mode,mode?Reads.breakReads_new2old[r][1]-Reads.breakReads_new2old[r-1][2]-1:0,leftUnique,leftLength,rightUnique,Integer.parseInt(str3.substring(0,str3.indexOf(","))),lengthThreshold,leftEnd1,leftEnd1_last,leftEnd1_str,leftEnd2_str,rightEnd1,rightEnd1_last,str1,str2,leftEnd2,leftEnd2_last,rightEnd2,rightEnd2_last,bw);
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
			leftLength=Reads.readLengths[r]-Integer.parseInt(str3.substring(str3.lastIndexOf(",")+1).trim());
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
			for (int i=1; i<=lastCharacter; i++) out+=SEPARATOR_MINOR+characters[i].toString();
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
		public boolean openStart, openEnd;  // The repeat might continue after start/end
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
		public final boolean implies(Character otherCharacter, int quantum) {
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