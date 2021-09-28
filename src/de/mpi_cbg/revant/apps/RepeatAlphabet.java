package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.*;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Intervals;


/**
 * Basic tools for handling an alphabet of repeat units.
 */
public class RepeatAlphabet {
	/**
	 * Marks non-repetitive blocks in a read
	 */
	public static final int NON_REPETITIVE = -1;
	
	/**
	 * Length of every repeat
	 */
	public static int[] repeatLengths;
	public static boolean[] isPeriodic;
	
	/**
	 * The recoded alphabet
	 */
	public static Character[] alphabet;
	public static int lastAlphabet;
	
	/**
	 * All alignments of a given readA
	 */
	private static AlignmentRow[] alignments;
	private static int lastAlignment;
	
	/**
	 * The recoded sequence of a given readA
	 */
	public static Character[] sequence;
	public static int lastInSequence;
	
	/**
	 * Length of the recoded sequence of every read
	 */
	private static int[] sequenceLengths;
	private static int lastSequenceLength;	
	
	/**
	 * Temporary space used by procedure $recode()$.
	 */
	private static int[] periodicIntervals, points;
	private static Character newCharacter;
	private static Tuple tmpTuple;
	
	/**
	 * Temporary space encoding an undirected graph and its connected components (see
	 * procedure $compactAlphabet()$).
	 */
	private static int[][] neighbors;
	private static int[] lastNeighbor;
	private static int[] stack, connectedComponent, connectedComponentSize;
	private static int nComponents;
	private static Character[] mergedChars;
	
	
	/**
	 * Loads read lengths, and repeat lengths and periodicity.
	 */
	public static final void initialize(int nReads, String readIDsFile, String readLengthsFile, int nRepeats, String repeatLengthsFile, String isPeriodicFile) throws IOException {
		Reads.nReads=nReads;
		Reads.maxReadLength=Reads.loadReadLengths(readLengthsFile);
		Reads.loadReadIDs(readIDsFile,Reads.nReads);
		loadRepeatLengths(repeatLengthsFile,nRepeats);
		loadIsPeriodic(isPeriodicFile,nRepeats);
	}
	
	
	/**
	 * Like $Reads.loadReadLengths()$, but for repeats.
	 */
	private static final int loadRepeatLengths(String repeatLengthsFile, int nRepeats) throws IOException {
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
	
	
	private static final void loadIsPeriodic(String isPeriodicFile, int nRepeats) throws IOException {
		isPeriodic = new boolean[nRepeats];
		BufferedReader br = new BufferedReader(new FileReader(isPeriodicFile));
		for (int i=0; i<nRepeats; i++) isPeriodic[i]=Integer.parseInt(br.readLine())==1;
		br.close();
	}
	
	
	/**
	 * Recodes every read, and collects in $alphabet[0..lastAlphabet]$ every distinct 
	 * character with discriminative power.
	 *
	 * Remark: the procedure assumes that $initialize()$ has already been called.
	 *
	 * @param alignmentsFile alignments between reads (readA) and repeats (readB);
	 * @param maxError alignments with error rate greater than this are discarded;
	 * @param read*File of reads, not of repeats.
	 */
	public static final void buildAlphabet(String alignmentsFile, double maxError, int distanceThreshold, int lengthThreshold) throws IOException {
		final int ALPHABET_CAPACITY = 100000;  // Arbitrary
		final int ALIGNMENTS_CAPACITY = 100000;  // Arbitrary
		final int SEQUENCE_CAPACITY = 1000000;  // Arbitrary
		boolean found;
		int i, j, k;
		int row, readA, previousReadA;
		String str;
		BufferedReader br;
		Character tmpCharacter;
		
		// Allocating memory
		if (alphabet==null || alphabet.length<ALPHABET_CAPACITY) alphabet = new Character[ALPHABET_CAPACITY];
		lastAlphabet=-1;
		if (alignments==null || alignments.length<ALIGNMENTS_CAPACITY) alignments = new AlignmentRow[ALIGNMENTS_CAPACITY];
		for (i=0; i<alignments.length; i++) {	
			if (alignments[i]==null) alignments[i] = new AlignmentRow();
		}
		if (sequence==null || sequence.length<SEQUENCE_CAPACITY) sequence = new Character[SEQUENCE_CAPACITY];
		for (i=0; i<sequence.length; i++) {
			if (sequence[i]==null) sequence[i] = new Character();
		}
		if (points==null || points.length<(ALIGNMENTS_CAPACITY)<<1) points = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (periodicIntervals==null || periodicIntervals.length<(ALIGNMENTS_CAPACITY)<<1) periodicIntervals = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (sequenceLengths==null || sequenceLengths.length<Reads.nReads) sequenceLengths = new int[Reads.nReads];
		if (newCharacter==null) newCharacter = new Character();
		if (tmpTuple==null) tmpTuple = new Tuple();
		
		// Collecting all characters (not necessarily distinct) from all recoded reads.
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); 
		previousReadA=-1; lastAlignment=-1; row=0;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments, "+(lastAlphabet+1)+" characters collected.");
			Alignments.readAlignmentFile(str);
			if ((2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.endB+2)>maxError) {
				str=br.readLine();
				continue;
			}
			readA=Alignments.readA-1;
			if (previousReadA==-1 || readA!=previousReadA) {
				if (previousReadA!=-1) {
					cleanAlignments(distanceThreshold);
					recodeRead(distanceThreshold);
					if (lastInSequence>0) addCharacters2Alphabet();
					sequenceLengths[previousReadA]=lastInSequence+1;
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
			recodeRead(distanceThreshold);
			if (lastInSequence>0) addCharacters2Alphabet();
			sequenceLengths[previousReadA]=lastInSequence+1;
		}
		
		// Sorting and compacting the alphabet
		compactAlphabet(distanceThreshold,lengthThreshold);
	}
	
	
	/**
	 * Stores in global variable $sequence$ a recoding of the read based on repeat 
	 * alignments. The procedure collects the endpoints of every maximal range of 
	 * overlapping periodic alignments, and of every non-periodic alignment that is not 
	 * contained in a periodic range and that does not straddle another non-periodic 
	 * alignment. Every other alignment is discarded. Points are clustered, and the 
	 * resulting chunks of the read become characters of $sequence$.
	 *
	 * Remark: the procedure assumes that reads do not contain long low-quality regions.
	 *
	 * Remark: endpoint clustering is very simplistic for now, and should be improved.
	 *
	 * Remark: the characters of $sequence$ are not objects in $alphabet$, and the 
	 * procedure does not look them up in $alphabet$: this is left to the caller.
	 *
	 * Remark: the procedure uses global arrays $alignments,periodicIntervals,points,
	 * sequence$.
	 */
	public static final void recodeRead(int distanceThreshold) {
		int i, j, k;
		int lastPeriodicInterval, currentStart, currentEnd;
		int firstJForNextI, inPeriodic, lastPoint, previousPoint, numerator, denominator;
		int startA, endA, newLength;
		final int lengthA = Reads.getReadLength(alignments[0].readA);
		
		if (periodicIntervals.length<(lastAlignment+1)<<1) periodicIntervals = new int[(lastAlignment+1)<<1];
		if (points.length<(lastAlignment+1)<<1) points = new int[(lastAlignment+1)<<1];
		AlignmentRow.order=AlignmentRow.ORDER_STARTA;
		Arrays.sort(alignments,0,lastAlignment+1);
		
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
		
		// Clustering readA endpoints of alignments
		i=0; j=0; firstJForNextI=-1; inPeriodic=-1; lastPoint=-1;
		while (i<=lastAlignment) {
			if (j>lastPeriodicInterval || periodicIntervals[j]>=alignments[i].endA-distanceThreshold) {
				if (inPeriodic==-1) {
					points[++lastPoint]=alignments[i].startA;
					points[++lastPoint]=alignments[i].endA;
				}
				else if (isPeriodic[alignments[i].readB]) {
					if (Math.abs(alignments[i].startA,periodicIntervals[inPeriodic])<=distanceThreshold) points[++lastPoint]=alignments[i].startA;
					if (Math.abs(alignments[i].endA,periodicIntervals[inPeriodic+1])<=distanceThreshold) points[++lastPoint]=alignments[i].endA;
				}
				i++; inPeriodic=-1;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (periodicIntervals[j+1]<alignments[i].startA-distanceThreshold) {
				j+=2;
				continue;
			}
			if (firstJForNextI==-1 && i<lastAlignment && periodicIntervals[j+1]>=alignments[i+1].startA) firstJForNextI=j;
			if ( Intervals.isApproximatelyContained(alignments[i].startA,alignments[i].endA,periodicIntervals[j],periodicIntervals[j+1]) ||
				 Intervals.areApproximatelyIdentical(alignments[i].startA,alignments[i].endA,periodicIntervals[j],periodicIntervals[j+1])
			   ) inPeriodic=j;
			j+=2;
		}
		if (inPeriodic==-1) {
			points[++lastPoint]=alignments[i].startA;
			points[++lastPoint]=alignments[i].endA;
		}
		else if (isPeriodic[alignments[i].readB]) {
			if (Math.abs(alignments[i].startA,periodicIntervals[inPeriodic])<=distanceThreshold) points[++lastPoint]=alignments[i].startA;
			if (Math.abs(alignments[i].endA,periodicIntervals[inPeriodic+1])<=distanceThreshold) points[++lastPoint]=alignments[i].endA;
		}
		Arrays.sort(points,0,lastPoint+1);
		j=-1; previousPoint=points[0]; numerator=points[0]; denominator=1;
		for (i=1; i<=lastPoint; i++) {
			if (points[i]<=previousPoint+distanceThreshold) {
				numerator+=points[i];
				denominator++;
				continue;
			}
			points[++j]=numerator/denominator;
			previousPoint=points[i]; numerator=points[i]; denominator=1;
		}
		points[++j]=numerator/denominator;
		lastPoint=j;
		
		// Creating the sequence
		lastInSequence=-1;
		// First non-repetitive character (if any).
		i=points[0];
		if (i>distanceThreshold) {
			lastInSequence=0;
			sequence[0].setNonrepetitive(i,true,i>=lengthA-distanceThreshold);
		}
		// Middle characters
		i=1; j=0; firstJForNextI=-1;
		while (i<=lastPoint) {
			startA=points[i-1]; endA=points[i];
			if (j>lastAlignment || alignments[j].startA>=endA) {
				if (newCharacter.lastTuple==-1) newCharacter.setNonrepetitive(endA-startA+1,startA<=distanceThreshold,endA>=lengthA-distanceThreshold);
				else newCharacter.cleanTuples(endA-startA+1);
				lastInSequence++;
				if (lastInSequence==sequence.length) {
					Character[] newSequence = new Character[sequence.length<<1];
					System.arraycopy(sequence,0,newSequence,0,sequence.length);
					for (k=sequence.length; k<newSequence.length; k++) newSequence[k] = new Character();
					sequence=newSequence;
				}
				sequence[lastInSequence].copyFrom(newCharacter);
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				newCharacter.reset();
				continue;
			}
			else if (alignments[j].endA<=startA) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastPoint && alignments[j].endA>endA) firstJForNextI=j;
			if ( Intervals.areApproximatelyIdentical(alignments[j].startA,alignments[j].endA,startA,endA) || 
				 Intervals.isApproximatelyContained(alignments[j].startA,alignments[j].endA,startA,endA)
			   ) newCharacter.addTuple(alignments[j],distanceThreshold,tmpTuple);
			j++;
		}
		// Last non-repetitive character (if any).
		i=points[lastPoint];
		if (lengthA-i>distanceThreshold) {
			lastInSequence++;
			if (lastInSequence==sequence.length) {
				Character[] newSequence = new Character[sequence.length<<1];
				System.arraycopy(sequence,0,newSequence,0,sequence.length);
				for (k=sequence.length; k<newSequence.length; k++) newSequence[k] = new Character();
				sequence=newSequence;
			}
			sequence[lastInSequence].setNonrepetitive(lengthA-i,i<=distanceThreshold,true);
		}
	}
	
	
	/**
	 * Appends all the characters of $sequence$ to the end of $alphabet$.
	 */
	private static final void addCharacters2Alphabet() {
		final int GROWTH_RATE = 100000;  // Arbitrary
		int i, j;
		final int newLength = lastAlphabet+lastInSequence+2;
		
		if (newLength>alphabet.length) {
			Character[] newAlphabet = new Character[newLength+GROWTH_RATE];
			System.arraycopy(alphabet,0,newAlphabet,0,alphabet.length);
			for (i=alphabet.length; i<newAlphabet.length; i++) newAlphabet[i] = new Character();
			alphabet=newAlphabet;
		}
		j=lastAlphabet;
		for (i=0; i<=lastInSequence; i++) alphabet[++j].copyFrom(sequence[i]);
		lastAlphabet=j;
	}
	
	
	/**
	 * Removes characters that can be considered substrings of others, and merges all 
	 * surviving characters by connected component.
	 *
	 * Remark: the procedure does not explicitly remove characters with no discriminative
	 * power, since it assumes they have not been added to $alphabet$.
	 */
	private static final void compactAlphabet(int distanceThreshold, int lengthThreshold) {
		int i, j;
		Character tmpChar;
		
		// Discarding characters that are implied by other characters
		Arrays.sort(alphabet,0,lastAlphabet+1);
		for (i=0; i<=lastAlphabet; i++) alphabet[i].id=1;
		for (i=0; i<=lastAlphabet; i++) {
			if (alphabet[i].isNonrepetitive()) continue;
			for (j=i+1; j<=lastAlphabet; j++) {
				if (!alphabet[j].sameRepeat(alphabet[i])) break;
				if (alphabet[j].isNonrepetitive() || alphabet[j].id==0) continue;
				if (alphabet[i].implies(alphabet[j],distanceThreshold,lengthThreshold)) alphabet[j].id=0;
			}
			for (j=i-1; j>=0; j--) {
				if (!alphabet[j].sameRepeat(alphabet[i])) break;
				if (alphabet[j].isNonrepetitive() || alphabet[j].id==0) continue;
				if (alphabet[i].implies(alphabet[j],distanceThreshold,lengthThreshold)) alphabet[j].id=0;
			}
		}
		j=-1;
		for (i=0; i<=lastAlphabet; i++) {
			if (alphabet[i].id==0) continue;
			j++;
			tmpChar=alphabet[j];
			alphabet[j]=alphabet[i];
			alphabet[i]=tmpChar;
		}
		lastAlphabet=j;
		
		// Merging surviving characters
		initializeGraph(lastAlphabet+1);
		for (i=1; i<lastAlphabet; i++) {
			for (j=i+1; j<=lastAlphabet; j++) {
				if (!alphabet[j].sameRepeat(alphabet[i])) break;
				if (alphabet[j].isSimilar(alphabet[i],distanceThreshold,lengthThreshold) && alphabet[j].sameOpen(alphabet[i])) {
					addEdge(i,j); addEdge(j,i);
				}
			}
		}
		nComponents=getConnectedComponent(lastAlphabet+1);
		if (mergedChars==null) {
			mergedChars = new Character[nComponents];
			for (i=0; i<mergedChars.length; i++) mergedChars[i] = new Character();
		}
		else if (mergedChars.length<nComponents) {
			Character[] newArray = new Character[nComponents];
			System.arraycopy(mergedChars,0,newArray,0,mergedChars.length);
			for (i=0; i<mergedChars.length; i++) newArray[i].reset();
			for (i=mergedChars.length; i<newArray.length; i++) newArray[i] = new Character();
			mergedChars=newArray;
		}
		for (i=0; i<=lastAlphabet; i++) {
			tmpChar=mergedChars[connectedComponent[i]];
			if (tmpChar.lastTuple==-1) tmpChar.copyFrom(alphabet[i]);
			else tmpChar.addTupleEnds(alphabet[i]);
		}
		for (i=0; i<nComponents; i++) mergedChars[i].normalizeTupleEnds(connectedComponentSize[i]);
		for (i=0; i<nComponents; i++) alphabet[i].copyFrom(mergedChars[i]);
		lastAlphabet=nComponents-1;
		
		// Sorting the results of the merge
		Arrays.sort(alphabet,0,lastAlphabet+1);
		for (i=0; i<=lastAlphabet; i++) alphabet[i].id=i;
	}
	
	
	public static final void serializeAlphabet(String path) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
		for (int i=0; i<=lastAlphabet; i++) bw.write(alphabet[i].toString()+"\n");
		bw.close();
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
		for (i=0; i<nNodes; i++) connectedComponent[i]=-1;
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
		for (i=0; i<nComponents; i++) connectedComponentSize[i]=0;
		for (i=0; i<nNodes; i++) connectedComponentSize[connectedComponent[i]]++;
	
		return nComponents;
	}
	
	
	/**
	 * A character of the recoded alphabet, which is a set of substrings of repeats in 
	 * specific orientations and with open/closed endpoints.
	 */
	private static class Character implements Comparable {
		private static final int CAPACITY = 10;  // Arbitrary
		
		public int id;
		public Tuple[] tuples;
		public int lastTuple;
		
		
		/**
		 * An empty character
		 */
		public Character() {
			tuples = new Tuple[CAPACITY];
			reset();
		}
		
		
		public void reset() { id=-1; lastTuple=-1; }
		
		
		/**
		 * A non-repetitive character of given length.
		 */
		public Character(int length, boolean openStart, boolean openEnd) {
			id=-1;
			tuples = new Tuple[CAPACITY];
			setNonrepetitive(length,openStart,openEnd);
		}
		
		
		/**
		 * Sets the current character to a non-repetitive of given length.
		 */
		public final void setNonrepetitive(int length, boolean openStart, boolean openEnd) {
			tuples[0] = new Tuple(NON_REPETITIVE,true,-1,-1,length,openStart,openEnd);
			lastTuple=0;
		}
		
		
		public final boolean isNonrepetitive() {
			return lastTuple==0 && tuples[0].repeat==NON_REPETITIVE;
		}
		
		
		private final void enlarge() {
			Tuple[] newTuples = new Tuple[tuples.length<<1];
			System.arraycopy(tuples,0,newTuples,0,tuples.length);
			for (int i=tuples.length; i<newTuples.length; i++) newTuples[i] = new Tuple();
			tuples=newTuples;
		}
		
		
		/**
		 * Makes this character contain the same data as $otherCharacter$ (but all 
		 * pointers stay distinct). Field $id$ is not changed.
		 */
		public final void copyFrom(Character otherCharacter) {
			int i;
			
			if (tuples==null) tuples = new Tuple[otherCharacter.lastTuple+1];			
			else if (tuples.length<otherCharacter.lastTuple+1) {
				Tuple[] newTuples = new Tuple[otherCharacter.lastTuple+1];
				System.arraycopy(tuples,0,newTuples,0,tuples.length);
				for (i=tuples.length; i<newTuples.length; i++) tuples[i] = new Tuple();
			}
			lastTuple=otherCharacter.lastTuple;
			for (i=0; i<=lastTuple; i++) tuples[i].copyFrom(otherCharacter.tuples[i]);
		}
		
		
		/**
		 * Sorting just by $repeat,orientation$ of the tuples.
		 */
		public int compareTo(Object other) {
			int i;
			int last;
			
			Character otherCharacter = (Character)other;
			if (lastTuple<otherCharacter.lastTuple) return -1;
			else if (lastTuple>otherCharacter.lastTuple) return 1;
			for (i=0; i<=lastTuple; i++) {
				if (tuples[i].repeat<otherCharacter.tuples[i].repeat) return -1;
				else if (tuples[i].repeat>otherCharacter.tuples[i].repeat) return 1;
				if (tuples[i].orientation && !otherCharacter.tuples[i].orientation) return -1;
				else if (!tuples[i].orientation && otherCharacter.tuples[i].orientation) return 1;
			}
			return 0;
		}
		
		
		public boolean equals(Object other) {
			Character otherCharacter = (Character)other;
			if (lastTuple!=otherCharacter.lastTuple) return false;
			for (int i=0; i<=lastTuple; i++) {
				if (!tuples[i].equals(otherCharacter.tuples[i])) return false;
			}
			return true;
		}
		
		
		/**
		 * @return TRUE iff $otherCharacter$ is identical to this character,
		 * when considering just the $repeat$ and $orientation$ fields of all tuples.
		 */
		public final boolean sameRepeat(Character otherCharacter) {
			if (lastTuple!=otherCharacter.lastTuple) return false;
			for (int i=0; i<=lastTuple; i++) {
				if (tuples[i].repeat!=otherCharacter.tuples[i].repeat || tuples[i].orientation!=otherCharacter.tuples[i].orientation) return false;
			}
			return true;
		}
		
		
		/**
		 * Remark: the procedure assumes that the two characters are similar 
		 * according to $sameRepeat()$.
		 */
		public final boolean isSimilar(Character otherCharacter, int distanceThreshold, int lengthThreshold) {
			for (int i=0; i<=lastTuple; i++) {
				if ( Math.abs(tuples[i].start-otherCharacter.tuples[i].start)>distanceThreshold || 
					 Math.abs(tuples[i].end-otherCharacter.tuples[i].end)>distanceThreshold || 
					 Math.abs(tuples[i].length-otherCharacter.tuples[i].length)>lengthThreshold
				   ) return false;
			}
			return true;
		}
		
		
		public final boolean sameOpen(Character otherCharacter) {
			if (lastTuple!=otherCharacter.lastTuple) return false;
			for (int i=0; i<=lastTuple; i++) {
				if (tuples[i].openStart!=otherCharacter.tuples[i].openStart || tuples[i].openEnd!=otherCharacter.tuples[i].openEnd) return false;
			}
			return true;			
		}
			
		
		
		/**
		 * Remark: the procedure assumes that the two characters are similar 
		 * acccording to $sameReadB()$.
		 *
		 * @return the Euclidean distance between the feature vectors of the two
		 * characters.
		 */
		public final double getDistance(Character otherCharacter) {
			int i;
			double out;
			
			out=0;
			for (i=0; i<=lastTuple; i++) {
				out+=Math.pow(tuples[i].start-otherCharacter.tuples[i].start,2);
				out+=Math.pow(tuples[i].end-otherCharacter.tuples[i].end,2);
				out+=Math.pow(tuples[i].length-otherCharacter.tuples[i].length,2);
			}
			return Math.sqrt(out);
		}
		
		
		public int hashCode() {
			String out = id+",";
			for (int i=0; i<=lastTuple; i++) out+=tuples[i].toString()+",";
			return out.hashCode();
		}
		
		
		public String toString() {
			String out = id+","+(lastTuple+1)+",";
			for (int i=0; i<=lastTuple; i++) out+=tuples[i].toString()+",";
			return out;
		}
		
		
		/**
		 * @return the average length of all tuples.
		 */
		public int getLength() {
			int i, out;
			
			out=0;
			for (i=0; i<=lastTuple; i++) out+=tuples[i].getLength();
			return out/(lastTuple+1);
		}
		
		
		/**
		 * Adds an alignment to the current character. If the character already contains 
		 * an alignment to the same repeat, the new alignment replaces the old one based
		 * on a canonical (and arbitrary) order.
		 *
		 * Remark: tuples might not be sorted after the procedure completes.
		 *
		 * @param tmpTuple temporary space.
		 */
		public final void addTuple(AlignmentRow alignment, int distanceThreshold, Tuple tmpTuple) {
			int i;
			int found;
			
			if (alignment.orientation) {
				tmpTuple.openStart=alignment.startA<=distanceThreshold;
				tmpTuple.openEnd=alignment.endA>=repeatLengths[alignment.readB]-distanceThreshold;
			}
			else {
				tmpTuple.openEnd=alignment.startA<=distanceThreshold;
				tmpTuple.openStart=alignment.endA>=repeatLengths[alignment.readB]-distanceThreshold;
			}
			tmpTuple.repeat=alignment.readB; tmpTuple.orientation=alignment.orientation;
			if (!isPeriodic[alignment.readB]) {
				tmpTuple.start=alignment.startB;
				tmpTuple.end=alignment.endB;
				tmpTuple.length=0;
			}
			else {
				tmpTuple.start=-1; tmpTuple.end=-1;
				tmpTuple.length=alignment.endA-alignment.startA+1; // A is correct here
			}
			found=-1;
			for (i=0; i<=lastTuple; i++) {
				if (tuples[i].repeat==tmpTuple.repeat && tuples[i].orientation==tmpTuple.orientation) {
					found=i;
					break;
				}
			}
			if (found==-1) {
				lastTuple++;
				if (lastTuple==tuples.length) enlarge();	
				tuples[lastTuple].copyFrom(tmpTuple);
			}
			else {
				if (isPeriodic[alignment.readB]) tuples[found].length=Math.max(tuples[found].length,tmpTuple.length);
				else if ( alignment.startB<tuples[found].start || 
						  (alignment.startB==tuples[found].start && alignment.endB<tuples[found].end)
			        	) {
					tuples[found].start=alignment.startB;
					tuples[found].end=alignment.endB;
				}
				if (tmpTuple.openStart) tuples[found].openStart=true;
				if (tmpTuple.openEnd) tuples[found].openEnd=true;
			}
		}
		
		
		/**
		 * (1) Sorts $tuples$. (2) Makes open marks homogeneous across all tuples. (3) If 
		 * the character contains a periodic repeat, removes all non-periodic repeats from
		 * the character and sets the length of every periodic repeat to $periodicLength$.
		 */
		public final void cleanTuples(int periodicLength) {
			boolean foundPeriodic, foundNonperiodic, openLeft, openRight;
			int i, j;
			Tuple tmpTuple;
			
			// Enforcing periodic/nonperiodic criteria
			foundPeriodic=false; foundNonperiodic=false;
			for (i=0; i<=lastTuple; i++) {
				if (isPeriodic[tuples[i].repeat]) {
					foundPeriodic=true;
					tuples[i].length=periodicLength;
				}
				else foundNonperiodic=true;
			}
			if (!foundPeriodic) {
				Arrays.sort(tuples,0,lastTuple+1);
				return;
			}
			if (foundNonperiodic) {
				j=-1;
				for (i=0; i<=lastTuple; i++) {
					if (!isPeriodic[tuples[i].repeat]) continue;
					j++;
					tmpTuple=tuples[j];
					tuples[j]=tuples[i];
					tuples[i]=tmpTuple;
				}
				lastTuple=j;
			}
			
			// Enforcing open criteria
			openLeft=false; openRight=false;
			for (i=0; i<=lastTuple; i++) {
				if ((tuples[i].orientation && tuples[i].openStart) ||  (!tuples[i].orientation && tuples[i].openEnd)) openLeft=true;
				if ((tuples[i].orientation && tuples[i].openEnd) ||  (!tuples[i].orientation && tuples[i].openStart)) openRight=true;				
			}
			for (i=0; i<=lastTuple; i++) {
				if (tuples[i].orientation) { tuples[i].openStart=openLeft; tuples[i].openEnd=openRight; }
				else { tuples[i].openStart=openRight; tuples[i].openEnd=openLeft; }
			}			
			
			Arrays.sort(tuples,0,lastTuple+1);
		}
		
		
		public final boolean openStart() { return tuples[0].openStart; }
		
		
		public final boolean openEnd() { return tuples[0].openEnd; }
		
		
		/**
		 * Adds the tuple ends of $otherCharacter$ to those of the current character,
		 * assuming that $otherCharacter$ satisfies $sameRepeat()$ and $isSimilar()$.
		 */
		public final void addTupleEnds(Character otherCharacter) {
			for (int i=0; i<=lastTuple; i++) {
				tuples[i].start+=otherCharacter.tuples[i].start;
				tuples[i].end+=otherCharacter.tuples[i].end;
				tuples[i].length+=otherCharacter.tuples[i].length;
			}
		}
		
		
		/**
		 * Divides every tuple end by $denominator$.
		 */
		public final void normalizeTupleEnds(int denominator) {
			for (int i=0; i<=lastTuple; i++) {
				tuples[i].start/=denominator; tuples[i].end/=denominator;
				tuples[i].length/=denominator;
			}
		}
		
		
		/**
		 * Remark: the two characters are assumed both to be repetitive and to satisfy
		 * $sameRepeat()$.
		 *
		 * @return TRUE iff no tuple of $otherCharacter$ can add information WRT the
		 * corresponding tuple of this character.
		 */
		public final boolean implies(Character otherCharacter, int distanceThreshold, int lengthThreshold) {
			for (int i=0; i<=lastTuple; i++) {
				if (tuples[i].start==-1) {  // Periodic
					if (Math.abs(tuples[i].length,otherCharacter.tuples[i].length)<=lengthThreshold) return false;
					else if (tuples[i].length<otherCharacter.tuples[i].length) return false;
					else if ( (!otherCharacter.tuples[i].openStart && !otherCharacter.tuples[i].openEnd) ||
						      (!otherCharacter.tuples[i].openStart && tuples[i].openStart) ||
					          (!otherCharacter.tuples[i].openEnd && tuples[i].openEnd)
							) return false;
				}
				else {  // Nonperiodic
					if (Intervals.areApproximatelyIdentical(otherCharacter.tuples[i].start,otherCharacter.tuples[i].end,tuples[i].start,tuples[i].end)) return false;
					else if (Intervals.isApproximatelyContained(otherCharacter.tuples[i].start,otherCharacter.tuples[i].end,tuples[i].start,tuples[i].end)) {
						if ( (Math.abs(tuples[i].start,otherCharacter.tuples[i].start)<=distanceThreshold && !tuples[i].openStart && otherCharacter.tuples[i].openStart) ||
							 (otherCharacter.tuples[i].start>tuples[i].start+distanceThreshold && !otherCharacter.tuples[i].openStart) ||
							 (Math.abs(tuples[i].end,otherCharacter.tuples[i].end)<=distanceThreshold && !tuples[i].openEnd && otherCharacter.tuples[i].openEnd) ||
							 (otherCharacter.tuples[i].end<tuples[i].end-distanceThreshold && !otherCharacter.tuples[i].openEnd)
						   ) return false;
					}
					else return false;
				}
			}
			return true;
		}
		
	}
	
	
	/**
	 * A substring of a repeat in a specific orientation
	 */
	private static class Tuple implements Comparable {
		public int repeat;
		public boolean orientation;
		public int start, end;
		public boolean openStart, openEnd;  // The repeat might continue beyond start/end
		public int length;  // >0: unique or periodic; 0: repetitive non-periodic.
		
		public Tuple() {
			repeat=-1; orientation=false; start=-1; end=-1; length=-1;
			openStart=true; openEnd=true;
		}
		
		public Tuple(int r, boolean o, int s, int e, int l, boolean os, boolean oe) {
			repeat=r; orientation=o; start=s; end=e; length=l;
			openStart=os; openEnd=oe;
		}
		
		public int compareTo(Object other) {
			Tuple otherTuple = (Tuple)other;
			if (repeat<otherTuple.repeat) return -1;
			else if (repeat>otherTuple.repeat) return 1;
			if (orientation && !otherTuple.orientation) return -1;
			else if (!orientation && otherTuple.orientation) return 1;
			if (start<otherTuple.start) return -1;
			else if (start>otherTuple.start) return 1;
			if (end<otherTuple.end) return -1;
			else if (end>otherTuple.end) return 1;
			if (length<otherTuple.length) return -1;
			else if (length>otherTuple.length) return 1;
			return 0;
		}
		
		public void copyFrom(Tuple otherTuple) {
			repeat=otherTuple.repeat;
			orientation=otherTuple.orientation;
			start=otherTuple.start; end=otherTuple.end;
			length=otherTuple.length;
			openStart=otherTuple.openStart; openEnd=otherTuple.openEnd;
		}
		
		public boolean equals(Object other) {
			Tuple otherTuple = (Tuple)other;
			return repeat==otherTuple.repeat && orientation==otherTuple.orientation && start==otherTuple.start && end==otherTuple.end && length==otherTuple.length && openStart==otherTuple.openStart && openEnd==otherTuple.openEnd;
		}
		
		public String toString() {
			return repeat+","+(orientation?"1":"0")+","+start+","+end+","+length+","+(openStart?"1":"0")+","+(openEnd?"1":"0");
		}
		
		public int getLength() {
			return (start==-1||end==-1)?length:end-start+1;
		}
	}
	
	
	
	
	
	
	
	
	// ----------------------------- ALIGNMENT PROCEDURES --------------------------------
	
	/**
	 * Basic filters on repeat alignments.
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
		
		// Removing non-periodic alignments that straddle other non-periodic alignments
		// (possibly with the same readB): these do not provide a clear signal.
		for (i=0; i<=lastAlignment; i++) alignments[i].flag=true;
		for (i=0; i<=lastAlignment; i++) {
			if (isPeriodic[alignments[i].readB]) continue;
			startA=alignments[i].startA; endA=alignments[i].endA;
			for (j=i+1; j<=lastAlignment; j++) {
				if (alignments[j].startA>=endA-distanceThreshold) break;
				if (isPeriodic[alignments[j].readB]) continue;
				if (Intervals.areApproximatelyIdentical(alignments[j].startA,alignments[j].endA,startA,endA)) continue;
				alignments[i].flag=false; alignments[j].flag=false;
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
		public static int order;
		
		public int readA, readB, startA, endA, startB, endB, diffs;
		public boolean orientation;
		
		/**
		 * Temporary space
		 */
		public boolean flag;
		
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
			return 0;
		}
		
		public String toString() {
			return readA+"["+startA+".."+endA+"] x "+readB+"["+startB+".."+endB+"] :: "+orientation;
		}	
	}

}