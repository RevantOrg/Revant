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
	public static final int NON_REPETITIVE = -1;
	public static final String SEPARATOR_MINOR = ",";
	public static final String SEPARATOR_MAJOR = "|";
	
	/**
	 * Length of every repeat
	 */
	public static int[] repeatLengths;
	public static boolean[] isPeriodic;
	
	/**
	 * The recoded alphabet
	 */
	public static Character[] alphabet;
	public static int lastNonrepetitive, lastPeriodic, lastAlphabet;
	public static int maxOpenLength_nonperiodic;
	
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
	public static int[] sequenceLengths;
	
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
	
	
	/**
	 * Recodes every read, and collects in $outputFile$ every character (not necessarily 
	 * distinct) with discriminative power.
	 *
	 * @param alignmentsFile alignments between reads (readA) and repeats (readB);
	 * @param maxError alignments with error rate greater than this are discarded.
	 */
	public static final void collectCharacterInstances(String alignmentsFile, double maxError, int distanceThreshold, int lengthThreshold, String outputFile) throws IOException {
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
		if (sequence==null || sequence.length<SEQUENCE_CAPACITY) sequence = new Character[SEQUENCE_CAPACITY];
		for (i=0; i<sequence.length; i++) {
			if (sequence[i]==null) sequence[i] = new Character();
		}
		if (points==null || points.length<(ALIGNMENTS_CAPACITY)<<1) points = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (periodicIntervals==null || periodicIntervals.length<(ALIGNMENTS_CAPACITY)<<1) periodicIntervals = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (sequenceLengths==null || sequenceLengths.length<MAX_SEQUENCE_LENGTH) sequenceLengths = new int[MAX_SEQUENCE_LENGTH];
		Math.set(sequenceLengths,MAX_SEQUENCE_LENGTH-1,0);
		if (newCharacter==null) newCharacter = new Character();
		if (tmpTuple==null) tmpTuple = new Tuple();
		
		// Collecting characters from all recoded reads
		bw = new BufferedWriter(new FileWriter(outputFile));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); 
		previousReadA=-1; lastAlignment=-1; row=0;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments, "+(lastAlphabet+1)+" characters collected.");
			Alignments.readAlignmentFile(str);
			if ((2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.endB+2)>maxError) {
				str=br.readLine(); row++;
				continue;
			}
			readA=Alignments.readA-1;
			if (previousReadA==-1 || readA!=previousReadA) {
				if (previousReadA!=-1) {
					cleanAlignments(distanceThreshold);
					if (lastAlignment!=-1) {
						recodeRead(distanceThreshold);
						if (lastInSequence>0) addCharacterInstances(bw);
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
				if (lastInSequence>0) addCharacterInstances(bw);
				sequenceLengths[lastInSequence+1<MAX_SEQUENCE_LENGTH?lastInSequence+1:MAX_SEQUENCE_LENGTH-1]++;
			}
			else sequenceLengths[0]++;
		}
		bw.close();
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
	 * Remark: the characters of $sequence$ are not objects in $alphabet$, and the 
	 * procedure does not look them up in $alphabet$: this is left to the caller.
	 *
	 * Remark: the procedure uses global arrays $alignments,periodicIntervals,points,
	 * sequence$.
	 */
	public static final void recodeRead(int distanceThreshold) {
		int i, j, k;
		int lastPeriodicInterval, currentStart, currentEnd;
		int firstJForNextI, inPeriodic, lastPoint;
		int startA, endA, newLength, component, nComponents;
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
		Arrays.sort(points,0,lastPoint+1);
		initializeGraph(lastPoint+1);
		for (i=0; i<lastPoint; i++) {
			for (j=i+1; j<=lastPoint; j++) {
				if (points[j]>points[i]+distanceThreshold) break;
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
	 * Appends the characters of $sequence$ to the end of $alphabet$.
	 */
	private static final void addCharacterInstances(BufferedWriter bw) throws IOException {	
		for (int i=0; i<=lastInSequence; i++) {
			if (sequence[i].isNonrepetitive()) {
				if (!sequence[i].openStart() && !sequence[i].openEnd()) {
					bw.write(sequence[i].toString());
					bw.newLine();
				}
				else maxOpenLength_nonperiodic=Math.max(maxOpenLength_nonperiodic,sequence[i].getLength());
			}
			else {
				bw.write(sequence[i].toString());
				bw.newLine();
			}
		}
	}
	
	
	/**
	 * Removes character instances that can be considered substrings of others, and merges
	 * all surviving character instances by connected component.
	 *
	 * Remark: the procedure does not explicitly remove characters with no discriminative
	 * power, since it assumes they have not been added to $alphabet$.
	 */
	public static final void compactInstances(int distanceThreshold, int lengthThreshold) {
		int i, j, k;
		Character tmpChar;
		Tuple tupleI;
		
		for (k=0; k<=lastAlphabet; k++) {
			if (!alphabet[k].isNonrepetitive()) break;
		}
		System.err.println("Discarding repetitive characters that are implied by other repetitive characters... ("+k+" non-repetitive characters)");
		for (i=k; i<=lastAlphabet; i++) alphabet[i].id=1;
		for (i=k; i<=lastAlphabet; i++) {
			if (i%1000==0) System.err.println("Processed "+i+" characters");
			tupleI=alphabet[i].tuples[0];
			for (j=i+1; j<=lastAlphabet; j++) {
				if (!alphabet[j].sameRepeat(alphabet[i]) || alphabet[j].tuples[0].implies_tooFarAfter(tupleI,distanceThreshold,lengthThreshold)) break;
				if (alphabet[j].id==0) continue;
				if (alphabet[i].implies(alphabet[j],distanceThreshold,lengthThreshold)) alphabet[j].id=0;
			}
			for (j=i-1; j>=0; j--) {
				if (!alphabet[j].sameRepeat(alphabet[i]) || alphabet[j].tuples[0].implies_tooFarBefore(tupleI,distanceThreshold,lengthThreshold)) break;
				if (alphabet[j].id==0) continue;
				if (alphabet[i].implies(alphabet[j],distanceThreshold,lengthThreshold)) alphabet[j].id=0;
			}
		}
		j=k-1;
		for (i=k; i<=lastAlphabet; i++) {
			if (alphabet[i].id==0) continue;
			j++;
			tmpChar=alphabet[j];
			alphabet[j]=alphabet[i];
			alphabet[i]=tmpChar;
		}
		lastAlphabet=j;
		System.err.println("DONE  "+(lastAlphabet+1)+" characters after filtering.");	
		
		System.err.println("Merging surviving characters... ");
		initializeGraph(lastAlphabet+1);
		for (i=0; i<lastAlphabet; i++) {
			if (i%1000==0) System.err.println("Processed "+i+" characters");
			tupleI=alphabet[i].tuples[0];
			for (j=i+1; j<=lastAlphabet; j++) {
				if (!alphabet[j].sameRepeat(alphabet[i]) || alphabet[j].tuples[0].tooFarAfter(tupleI,distanceThreshold,lengthThreshold)) break;
				if (alphabet[j].isSimilar(alphabet[i],distanceThreshold,lengthThreshold) && alphabet[j].sameOpen(alphabet[i])) {
					addEdge(i,j); addEdge(j,i);
				}
			}
		}
		System.err.print("Graph built. Computing connected components... ");
		nComponents=getConnectedComponent(lastAlphabet+1);
		System.err.println("DONE");
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
		System.err.print("DONE  "+(lastAlphabet+1)+" characters after the merge.");
		
		System.err.print("Sorting the results of the merge... ");
		Arrays.sort(alphabet,0,lastAlphabet+1);
		lastNonrepetitive=-1; lastPeriodic=-1;
		i=0;
		while (i<=lastAlphabet && alphabet[i].isNonrepetitive()) {
			alphabet[i].id=i;
			lastNonrepetitive=i;
			i++;
		}
		while (i<=lastAlphabet && alphabet[i].tuples[0].start==-1) {
			alphabet[i].id=i;
			lastPeriodic=i;
			i++;
		}
		while (i<=lastAlphabet) {
			alphabet[i].id=i;
			i++;
		}
		System.err.println("DONE");
	}
	
	
	public static final void serializeAlphabet(String path) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
		bw.write(lastNonrepetitive+SEPARATOR_MINOR+lastPeriodic+SEPARATOR_MINOR+lastAlphabet+SEPARATOR_MINOR+maxOpenLength_nonperiodic+"\n");
		for (int i=0; i<=lastAlphabet; i++) bw.write(alphabet[i].toString()+"\n");
		bw.close();
	}
	
	
	public static final void deserializeAlphabet(String path) throws IOException {
		int i, p, q;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();
		p=str.indexOf(SEPARATOR_MINOR);
		lastNonrepetitive=Integer.parseInt(str.substring(0,p));
		p++; q=str.indexOf(SEPARATOR_MINOR,p+1);
		lastPeriodic=Integer.parseInt(str.substring(p,q));
		p=q+1; q=str.indexOf(SEPARATOR_MINOR,p+1);
		lastAlphabet=Integer.parseInt(str.substring(p,q));
		maxOpenLength_nonperiodic=Integer.parseInt(str.substring(q+1));
		alphabet = new Character[lastAlphabet+1];
		for (i=0; i<=lastAlphabet; i++) {
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
	 * Like $buildAlphabet()$, but looks up in the alphabet every character of a recoded 
	 * read.
	 *
	 * Remark: this procedure should be parallelized.
	 *
	 * @param outputFile contains the translation of one read per line.
	 */
	public static final void translateReads(String alignmentsFile, double maxError, int distanceThreshold, int lengthThreshold, String outputFile) throws IOException {
		final int ALIGNMENTS_CAPACITY = 100000;  // Arbitrary
		final int SEQUENCE_CAPACITY = 1000000;  // Arbitrary
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
		if (sequence==null || sequence.length<SEQUENCE_CAPACITY) sequence = new Character[SEQUENCE_CAPACITY];
		for (i=0; i<sequence.length; i++) {
			if (sequence[i]==null) sequence[i] = new Character();
		}
		if (points==null || points.length<(ALIGNMENTS_CAPACITY)<<1) points = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (periodicIntervals==null || periodicIntervals.length<(ALIGNMENTS_CAPACITY)<<1) periodicIntervals = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (newCharacter==null) newCharacter = new Character();
		if (tmpTuple==null) tmpTuple = new Tuple();
		
		// Translating every read using the alphabet
		bw = new BufferedWriter(new FileWriter(outputFile));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); 
		previousReadA=-1; lastAlignment=-1; row=0;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			if ((2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.endB+2)>maxError) {
				str=br.readLine(); row++;
				continue;
			}
			readA=Alignments.readA-1;
			if (previousReadA==-1 || readA!=previousReadA) {
				if (previousReadA!=-1) {
					cleanAlignments(distanceThreshold);
					if (lastAlignment!=-1) {
						recodeRead(distanceThreshold);
						if (lastInSequence<=0) bw.newLine();
						else translateRead(bw,distanceThreshold,lengthThreshold);
					}
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
				if (lastInSequence<=0) bw.newLine();
				else translateRead(bw,distanceThreshold,lengthThreshold);
			}
		}
		bw.close();
	}
	
	
	private static final void translateRead(BufferedWriter bw, int distanceThreshold, int lengthThreshold) throws IOException {
		int i;
		
		if (sequence[0].isNonrepetitive()) translate_nonrepetitive(sequence[0],lengthThreshold,bw);
		else if (sequence[0].tuples[0].start==-1) translate_periodic(sequence[0],lengthThreshold,bw);
		else translate(sequence[0],distanceThreshold,bw);
		for (i=1; i<=lastInSequence; i++) {
			bw.write(SEPARATOR_MAJOR);
			if (sequence[i].isNonrepetitive()) translate_nonrepetitive(sequence[i],lengthThreshold,bw);
			else if (sequence[i].tuples[0].start==-1) translate_periodic(sequence[i],lengthThreshold,bw);
			else translate(sequence[i],distanceThreshold,bw);
		}
		bw.newLine();
	}

	
	/**
	 * Writes: $X$ if $character$ is closed and has length most similar to $X$; $-1-X$ if
	 * $character$ is half-open and could be an instance of any closed character $>=X$;
	 * $lastAlphabet+1$ if $character$ is half-open but longer than any closed character.
	 */
	private static final void translate_nonrepetitive(Character character, int lengthThreshold, BufferedWriter bw) throws IOException {
		int i;
		final int length = character.tuples[0].length;
		
		i=Arrays.binarySearch(alphabet,0,lastNonrepetitive+1,character);
		if (i<0) i=-1-i;
		if (character.isOpen()) {
			if (i==lastNonrepetitive+1) bw.write((lastAlphabet+1)+"");
			else {
				while (i>=0 && alphabet[i].tuples[0].getLength()>=length-lengthThreshold) i--;
				i++;
				bw.write((-1-i)+"");
			}
		}
		else {
			if (i==0) bw.write("0");
			else if (i<=lastNonrepetitive) {
				if (alphabet[i].tuples[0].length-length<length-alphabet[i-1].tuples[0].length) bw.write(i+"");
				else bw.write((i-1)+"");
			}
			else bw.write((i-1)+"");
		}
	}
	
	
	/**
	 * Same logic as $translateRead_nonrepetitive()$.
	 */
	private static final void translate_periodic(Character character, int lengthThreshold, BufferedWriter bw) throws IOException {
		int i;
		final int length = character.getLength();
		
		i=Arrays.binarySearch(alphabet,lastNonrepetitive+1,lastPeriodic+1,character);
		if (i<0) i=-1-i;
		if (character.isOpen()) {
			if (!alphabet[i].sameRepeat(character)) i--;
			while (i>lastNonrepetitive && alphabet[i].sameRepeat(character) && alphabet[i].getLength()>=length-lengthThreshold) i--;
			i++;
			if (!alphabet[i].sameRepeat(character)) {
				System.err.println("translateRead_periodic> ERROR: open periodic character not found in the alphabet");
				System.err.println("query: "+character);
				System.err.println("first candidate in alphabet: "+alphabet[i]);
				System.exit(1);
			}
			bw.write((-1-i)+"");
		}
		else {
			if (i<=lastPeriodic && alphabet[i].sameRepeat(character)) {
				if (i-1>lastNonrepetitive && alphabet[i-1].sameRepeat(character)) {
					if (alphabet[i].getLength()-length<length-alphabet[i-1].getLength()) bw.write(i+"");
					else bw.write((i-1)+"");
				}
				else bw.write(i+"");
			}
			else {
				if (i-1>lastNonrepetitive && alphabet[i-1].sameRepeat(character)) bw.write((i-1)+"");
				else {
					System.err.println("translateRead_periodic> ERROR: closed periodic character not found in the alphabet: "+character);
					System.exit(1);
				}
			}
		}
	}
	
	
	/**
	 * For nonperiodic repeats. If $character$ is half-open, prints the sorted list of all
	 * characters in $alphabet$ that are similar to or imply $character$. If $character$ 
	 * is closed, prints the sorted list of all closed characters in $alphabet$ that are 
	 * similar to $character$.
	 */
	private static final void translate(Character character, int distanceThreshold, BufferedWriter bw) throws IOException {
		final int CAPACITY = 100;  // Arbitrary
		int i, j;
		int lastCharacter;
		
		if (stack==null) stack = new int[CAPACITY];
		lastCharacter=-1;
		i=Arrays.binarySearch(alphabet,lastPeriodic+1,lastAlphabet+1,character);
		if (i<0) i=-1-i;
		if (character.isOpen()) {
			for (j=i+1; j<=lastAlphabet; j++) {
				if (!alphabet[j].sameRepeat(character) || alphabet[j].tuples[0].implies_tooFarAfter(character.tuples[0],distanceThreshold,Math.POSITIVE_INFINITY)) break;
				if (alphabet[j].implies(character,distanceThreshold,Math.POSITIVE_INFINITY) || alphabet[j].isSimilar(character,distanceThreshold,Math.POSITIVE_INFINITY)) {
					lastCharacter++;
					if (lastCharacter==stack.length) {
						int[] newArray = new int[stack.length<<1];
						System.arraycopy(stack,0,newArray,0,stack.length);
						stack=newArray;
					}
					stack[lastCharacter]=j;
				}
			}
			for (j=i-1; j>lastPeriodic; j--) {
				if (!alphabet[j].sameRepeat(character) || alphabet[j].tuples[0].implies_tooFarBefore(character.tuples[0],distanceThreshold,Math.POSITIVE_INFINITY)) break;
				if (alphabet[j].implies(character,distanceThreshold,Math.POSITIVE_INFINITY) || alphabet[j].isSimilar(character,distanceThreshold,Math.POSITIVE_INFINITY)) {
					lastCharacter++;
					if (lastCharacter==stack.length) {
						int[] newArray = new int[stack.length<<1];
						System.arraycopy(stack,0,newArray,0,stack.length);
						stack=newArray;
					}
					stack[lastCharacter]=j;
				}
			}
			if (lastCharacter>0) Arrays.sort(stack,0,lastCharacter+1);
			bw.write(stack[0]+"");
			for (j=1; j<=lastCharacter; j++) bw.write(SEPARATOR_MINOR+stack[j]);
		}
		else {
			lastCharacter=-1;
			for (j=i+1; j<=lastAlphabet; j++) {
				if (!alphabet[j].sameRepeat(character) || alphabet[j].tuples[0].implies_tooFarAfter(character.tuples[0],distanceThreshold,Math.POSITIVE_INFINITY)) break;
				if (alphabet[j].isSimilar(character,distanceThreshold,Math.POSITIVE_INFINITY) && alphabet[j].sameOpen(character)) {
					lastCharacter++;
					if (lastCharacter==stack.length) {
						int[] newArray = new int[stack.length<<1];
						System.arraycopy(stack,0,newArray,0,stack.length);
						stack=newArray;
					}
					stack[lastCharacter]=j;
				}
			}
			for (j=i-1; j>lastPeriodic; j--) {
				if (!alphabet[j].sameRepeat(character) || alphabet[j].tuples[0].implies_tooFarBefore(character.tuples[0],distanceThreshold,Math.POSITIVE_INFINITY)) break;
				if (alphabet[j].isSimilar(character,distanceThreshold,Math.POSITIVE_INFINITY) && alphabet[j].sameOpen(character)) {
					lastCharacter++;
					if (lastCharacter==stack.length) {
						int[] newArray = new int[stack.length<<1];
						System.arraycopy(stack,0,newArray,0,stack.length);
						stack=newArray;
					}
					stack[lastCharacter]=j;
				}
			}
			if (lastCharacter>0) Arrays.sort(stack,0,lastCharacter+1);
			bw.write(stack[0]+"");
			for (j=1; j<=lastCharacter; j++) bw.write(SEPARATOR_MINOR+stack[j]);
		}
	}
	
	
	
	
	
	
	
	
	
	/**
	 * A character of the recoded alphabet, which is a set of substrings of repeats in 
	 * specific orientations and with open/closed endpoints.
	 */
	public static class Character implements Comparable {
		private static final int CAPACITY = 10;  // Arbitrary
		public static int order;
		
		public int id;
		public Tuple[] tuples;
		public int lastTuple;
		
		
		/**
		 * An empty character
		 */
		public Character() {
			tuples = new Tuple[CAPACITY];
			for (int i=0; i<tuples.length; i++) tuples[i] = new Tuple();
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
			
			if (tuples==null) {
				tuples = new Tuple[otherCharacter.lastTuple+1];
				for (i=0; i<tuples.length; i++) tuples[i] = new Tuple();
			}
			else if (tuples.length<otherCharacter.lastTuple+1) {
				Tuple[] newTuples = new Tuple[otherCharacter.lastTuple+1];
				System.arraycopy(tuples,0,newTuples,0,tuples.length);
				for (i=tuples.length; i<newTuples.length; i++) newTuples[i] = new Tuple();
				tuples=newTuples;
			}
			lastTuple=otherCharacter.lastTuple;
			for (i=0; i<=lastTuple; i++) tuples[i].copyFrom(otherCharacter.tuples[i]);
		}
		
		
		/**
		 * First sorts by nonrepetitive, periodic, nonperiodic; then by number of tuples;
		 * then lexicographically on the tuples.
		 */
		public int compareTo(Object other) {
			final boolean isNonrepetitive = isNonrepetitive();
			final boolean isPeriodic = tuples[0].start==-1;
			final Character otherCharacter = (Character)other;
			final boolean otherCharacterIsNonrepetitive = otherCharacter.isNonrepetitive();
			final boolean otherCharacterIsPeriodic = otherCharacter.tuples[0].start==-1;
			int i, j;
			
			if (isNonrepetitive && !otherCharacterIsNonrepetitive) return -1;
			else if (!isNonrepetitive && otherCharacterIsNonrepetitive) return 1;
			if (isPeriodic && !otherCharacterIsPeriodic) return -1;
			else if (!isPeriodic && otherCharacterIsPeriodic) return 1;
			if (lastTuple<otherCharacter.lastTuple) return -1;
			else if (lastTuple>otherCharacter.lastTuple) return 1;
			for (i=0; i<=lastTuple; i++) {
				j=tuples[i].compareTo(otherCharacter.tuples[i]);
				if (j!=0) return j;
			}
			return 0;
		}
		
		
		/**
		 * Sorting strings lexicographically induces the same order as in $compareTo()$.
		 */
		public String toString() {
			String out = (isNonrepetitive()?"0":"1")+SEPARATOR_MINOR+(tuples[0].start==-1?"0":"1")+SEPARATOR_MINOR+lastTuple+SEPARATOR_MINOR;
			out+=tuples[0].toString();
			for (int i=1; i<=lastTuple; i++) out+=SEPARATOR_MINOR+tuples[i].toString();
			return out;
		}
		
		
		/**
		 * Symmetrical to $toString()$.
		 */
		public void deserialize(String str) {
			int i, p, q;
			
			id=-1;
			p=str.indexOf(SEPARATOR_MINOR);  // Skipping
			p=str.indexOf(SEPARATOR_MINOR,p+1);  // Skipping
			q=str.indexOf(SEPARATOR_MINOR,p+1);
			lastTuple=Integer.parseInt(str.substring(p+1,q));
			if (tuples==null) {
				tuples = new Tuple[lastTuple+1];
				for (i=0; i<=lastTuple; i++) tuples[i] = new Tuple();
			}
			else if (tuples.length<lastTuple+1) {
				Tuple[] newArray = new Tuple[lastTuple+1];
				System.arraycopy(tuples,0,newArray,0,tuples.length);
				for (i=tuples.length; i<=lastTuple; i++) newArray[i] = new Tuple();
				tuples=newArray;
			}
			p=q+1;
			for (i=0; i<=lastTuple; i++) p=tuples[i].deserialize(str,p);
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
				tmpTuple.openEnd=alignment.endA>=Reads.getReadLength(alignment.readA)-distanceThreshold;
			}
			else {
				tmpTuple.openEnd=alignment.startA<=distanceThreshold;
				tmpTuple.openStart=alignment.endA>=Reads.getReadLength(alignment.readA)-distanceThreshold;
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
		
		
		public final boolean isOpen() { return tuples[0].openStart || tuples[0].openEnd; }

		
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
		 * Remark: when changing this procedure, update $implies_tooFar*()$ as well.
		 *
		 * Remark: the two characters are assumed both to be repetitive and to satisfy
		 * $sameRepeat()$.
		 *
		 * Remark: the procedure assumes that both characters are nonrepetitive. 
		 * Nonrepetitive characters are not assumed to imply one another, since fully-open
		 * and half-open nonrepetitive characters are assumed to have been discarded.
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
	public static class Tuple implements Comparable {
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
		
		/**
		 * Assumes the same order as in $compareTo()$.
		 */
		public final boolean tooFarAfter(Tuple otherTuple, int distanceThreshold, int lengthThreshold) {
			return start>otherTuple.start+distanceThreshold || end>otherTuple.end+distanceThreshold || length>otherTuple.length+lengthThreshold;
		}
		
		public final void copyFrom(Tuple otherTuple) {
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
		
		/**
		 * Symmetrical to $toString()$.
		 *
		 * @param p the first position of $str$ to be read;
		 * @return the value of $p$ when the procedure completes.
		 */
		public int deserialize(String str, int p) {
			int q;
			
			q=str.indexOf(SEPARATOR_MINOR,p+1);
			repeat=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SEPARATOR_MINOR,p+1);
			orientation=Integer.parseInt(str.substring(p,q))==1;
			p=q+1; q=str.indexOf(SEPARATOR_MINOR,p+1);
			start=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SEPARATOR_MINOR,p+1);
			end=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SEPARATOR_MINOR,p+1);
			length=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(SEPARATOR_MINOR,p+1);
			openStart=Integer.parseInt(str.substring(p,q))==1;
			p=q+1; q=str.indexOf(SEPARATOR_MINOR,p+1);
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
		 * Assumes that the container characters have the same repeats, and uses the same
		 * criteria as $Character.implies()$.
		 */
		public final boolean implies_tooFarAfter(Tuple otherTuple, int distanceThreshold, int lengthThreshold) {
			if (start==-1) {  // Periodic
				return length>otherTuple.length+lengthThreshold;
			}
			else {  // Nonperiodic
				return start>=otherTuple.end-distanceThreshold;
			}
		}
		
		public final boolean implies_tooFarBefore(Tuple otherTuple, int distanceThreshold, int lengthThreshold) {
			if (start==-1) {  // Periodic
				return false;
			}
			else {  // Nonperiodic
				return start<otherTuple.start-distanceThreshold;
			}
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
// ----> This is not satisfactory, since repeats can straddle one another...
// Either do like Arne, explicitly modeling this as smaller characters.
// Or find maximal ranges of straddling alignments like periodic, and assign them a specific type.
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