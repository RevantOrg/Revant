package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Arrays;
import java.util.Hashtable;

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
	public static int maxOpenLength_nonperiodic;
	
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
	 * Temporary space used by procedure $recode()$.
	 */
	private static int[] periodicIntervals, points;
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
		maxOpenLength_nonperiodic=0;
		bw = new BufferedWriter(new FileWriter(outputFile));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); 
		previousReadA=-1; lastAlignment=-1; row=0;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			if ((2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2)>maxError) {
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
	 * resulting blocks of the read become elements of $sequence$.
	 *
	 * Remark: the procedure assumes that reads do not contain long low-quality regions.
	 *
	 * Remark: the characters in blocks of $sequence$ are not objects in $alphabet$, and 
	 * the procedure does not look them up in $alphabet$: this is left to the caller.
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
		if (newBlock==null) newBlock = new Block();
		lastInSequence=-1;
		// First non-repetitive block (if any).
		i=points[0];
		if (i>distanceThreshold) {
			lastInSequence=0;
			sequence[0].setUnique(i,true,i>=lengthA-distanceThreshold);
		}
		// Middle blocks
		i=1; j=0; firstJForNextI=-1;
		while (i<=lastPoint) {
			startA=points[i-1]; endA=points[i];
			if (j>lastAlignment || alignments[j].startA>=endA) {
				if (newBlock.lastCharacter==-1) newBlock.setUnique(endA-startA+1,startA<=distanceThreshold,endA>=lengthA-distanceThreshold);
				else newBlock.cleanCharacters(endA-startA+1);
				lastInSequence++;
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
				 Intervals.isApproximatelyContained(alignments[j].startA,alignments[j].endA,startA,endA)
			   ) newBlock.addCharacter(alignments[j],distanceThreshold,tmpCharacter);
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
				if (!sequence[i].openStart() && !sequence[i].openEnd()) {
					bw.write(sequence[i].characters[0].toString());
					bw.newLine();
				}
				else maxOpenLength_nonperiodic=Math.max(maxOpenLength_nonperiodic,sequence[i].getLength());
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
	 * Quantizes the endpoints of every repeat and removes character instances that can be
	 * considered substrings of others.
	 *
	 * Remark: the procedure does not explicitly remove characters with no discriminative
	 * power, since it assumes they have not been added to $alphabet$ in ther first place.
	 */
	public static final void compactInstances(int distanceThreshold, int lengthThreshold) {
		final int QUANTUM = IO.quantum;
		int i, j, k;
		int first;
		Character tmpChar;
		
		System.err.println("Quantizing characters... ");
		for (i=0; i<=lastAlphabet; i++) {
boolean fabio = alphabet[i].repeat==1114 && alphabet[i].start==0 && alphabet[i].end==835;
if (fabio) System.err.println("compactInstances> 1  original character:  "+alphabet[i]);
			alphabet[i].quantize(QUANTUM);
if (fabio) System.err.println("compactInstances> 1  quantized character: "+alphabet[i]);
		}
		Arrays.sort(alphabet,0,lastAlphabet+1);
		j=0;
		for (i=1; i<=lastAlphabet; i++) {
boolean fabio = (alphabet[i].repeat==1114 && alphabet[i].start==0 && alphabet[i].end==800) ||
	            (alphabet[i].repeat==1114 && alphabet[i].start==0 && alphabet[i].end==1400);
			if ( alphabet[i].repeat!=alphabet[j].repeat || alphabet[i].orientation!=alphabet[j].orientation || 
				 alphabet[i].start!=alphabet[j].start || alphabet[i].end!=alphabet[j].end || alphabet[i].length!=alphabet[j].length ||
			     alphabet[i].openStart!=alphabet[j].openStart || alphabet[i].openEnd!=alphabet[j].openEnd
			   ) {
if (fabio) System.err.println("compactInstances> 2  character added to quantized alphabet: "+alphabet[i]);				   
				j++;
				tmpChar=alphabet[j];
				alphabet[j]=alphabet[i];
				alphabet[i]=tmpChar;
			}
else {
	if (fabio) System.err.println("compactInstances> 3  character NOT ADDED to quantized alphabet: "+alphabet[i]);
}			
		}
		lastAlphabet=j;
		lastUnique=-1; i=0;
		while (i<=lastAlphabet && alphabet[i].repeat==UNIQUE) {
			alphabet[i].id=i;
			lastUnique=i;
			i++;
		}
		lastPeriodic=lastUnique;
		while (i<=lastAlphabet && alphabet[i].start==-1) {
			alphabet[i].id=i;
			lastPeriodic=i;
			i++;
		}
		while (i<=lastAlphabet) {
			alphabet[i].id=i;
			i++;
		}
		System.err.println("DONE");
		
		System.err.println("Discarding repetitive characters that are implied by other repetitive characters... ("+(lastUnique+1)+" non-repetitive characters)");
		for (i=lastUnique+1; i<=lastAlphabet; i++) alphabet[i].id=1;
		for (i=lastUnique+1; i<=lastAlphabet; i++) {
			if (i%1000==0) System.err.println("Processed "+i+" characters");
			for (j=i+1; j<=lastAlphabet; j++) {
				if ( alphabet[j].repeat!=alphabet[i].repeat || alphabet[j].orientation!=alphabet[i].orientation || 
					 alphabet[j].implies_tooFarAfter(alphabet[i])
				   ) break;
				if (alphabet[j].id==0) continue;
				if (alphabet[i].implies(alphabet[j],distanceThreshold)) {
					alphabet[j].id=0;
if ( (alphabet[j].repeat==1114 && alphabet[j].start==0 && alphabet[j].end==800) ||
	 (alphabet[i].repeat==1114 && alphabet[i].start==0 && alphabet[i].end==1400)
   ) System.err.println("compactInstances> 4  character "+alphabet[j]+" is implied by character "+alphabet[i]);
				}
			}
			for (j=i-1; j>=0; j--) {
				if (alphabet[j].repeat!=alphabet[i].repeat || alphabet[j].orientation!=alphabet[i].orientation) break;
				if (alphabet[j].id==0) continue;
				if (alphabet[i].implies(alphabet[j],distanceThreshold)) {
					alphabet[j].id=0;
if ( (alphabet[j].repeat==1114 && alphabet[j].start==0 && alphabet[j].end==800) ||
	 (alphabet[i].repeat==1114 && alphabet[i].start==0 && alphabet[i].end==1400)
   ) System.err.println("compactInstances> 5  character "+alphabet[j]+" is implied by character "+alphabet[i]);
				}
			}
		}
		j=lastUnique;
		for (i=lastUnique+1; i<=lastAlphabet; i++) {
			if (alphabet[i].id==0) continue;
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
	
	
	public static final void serializeAlphabet(String path) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
		bw.write(lastUnique+(SEPARATOR_MINOR+"")+lastPeriodic+(SEPARATOR_MINOR+"")+lastAlphabet+(SEPARATOR_MINOR+"")+maxOpenLength_nonperiodic+"\n");
		for (int i=0; i<=lastAlphabet; i++) bw.write(alphabet[i].toString()+"\n");
		bw.close();
	}
	
	
	public static final void deserializeAlphabet(String path) throws IOException {
		int i, p, q;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();
		p=str.indexOf(SEPARATOR_MINOR+"");
		lastUnique=Integer.parseInt(str.substring(0,p));
		p++; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
		lastPeriodic=Integer.parseInt(str.substring(p,q));
		p=q+1; q=str.indexOf(SEPARATOR_MINOR+"",p+1);
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
	 * Like $collectCharacterInstances()$, but looks up in the alphabet every character of
	 * every block of a recoded read.
	 *
	 * @param lastTranslatedRead the index in $Reads.readIDs$ of the last read that has 
	 * already been translated (-1 if no read has been translated yet);
	 * @param outputFile contains the translation of one read per line.
	 */
	public static final void translateReads(String alignmentsFile, int lastTranslatedRead, double maxError, int quantum, String outputFile) throws IOException {
		final int ALIGNMENTS_CAPACITY = 100000;  // Arbitrary
		final int SEQUENCE_CAPACITY = 1000000;  // Arbitrary
		int i, j;
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
		if (newBlock==null) newBlock = new Block();
		if (tmpCharacter==null) tmpCharacter = new Character();
		
		// Translating every read using the alphabet
		bw = new BufferedWriter(new FileWriter(outputFile));
		br = new BufferedReader(new FileReader(alignmentsFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); 
		previousReadA=-1; lastAlignment=-1; row=0; j=lastTranslatedRead+1;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			if ((2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2)>maxError) {
				str=br.readLine(); row++;
				continue;
			}
			readA=Alignments.readA-1;
			if (previousReadA==-1 || readA!=previousReadA) {
				if (previousReadA!=-1) {
					while (Reads.readIDs[j]<previousReadA) {
						bw.newLine();
						j++;
					}
					cleanAlignments(quantum);
					if (lastAlignment!=-1) {
						recodeRead(quantum);
						if (lastInSequence<=0) bw.newLine();
						else {
System.err.println("Translating read "+previousReadA);							
							translateRead(bw,quantum);
						}
					}
					else bw.newLine();
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
				bw.newLine();
				j++;
			}
			cleanAlignments(quantum);
			if (lastAlignment!=-1) {
				recodeRead(quantum);
				if (lastInSequence<=0) bw.newLine();
				else {
System.err.println("Translating read "+previousReadA);
					translateRead(bw,quantum);
				}
			}
			else bw.newLine();
			j++;
		}
		bw.close();
	}
	
	
	private static final void translateRead(BufferedWriter bw, int quantum) throws IOException {
		int i, j;
				
		for (i=0; i<=sequence[0].lastCharacter; i++) sequence[0].characters[i].quantize(quantum);
		if (sequence[0].isUnique()) translate_unique(sequence[0],bw);
		else if (sequence[0].characters[0].start==-1) translate_periodic(sequence[0],bw);
		else translate(sequence[0],quantum,bw);
		for (i=1; i<=lastInSequence; i++) {
			bw.write(SEPARATOR_MAJOR+"");
			for (j=0; j<=sequence[i].lastCharacter; j++) sequence[i].characters[j].quantize(quantum);
			if (sequence[i].isUnique()) translate_unique(sequence[i],bw);
			else if (sequence[i].characters[0].start==-1) translate_periodic(sequence[i],bw);
			else translate(sequence[i],quantum,bw);
		}
		bw.newLine();
	}

	
	/**
	 * Writes: $X$ if $block$ is closed and is identical to character $X$; $-1-X$ if
	 * $block$ is half-open and could be an instance of any closed block $>=X$;
	 * $lastAlphabet+1$ if $block$ is half-open but longer than any closed block.
	 *
	 * @param block its character is assumed to be quantized.
	 */
	private static final void translate_unique(Block block, BufferedWriter bw) throws IOException {
		int i;
		final Character character = block.characters[0];
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
				System.err.println("translateRead_unique> ERROR: closed unique character not found in the alphabet");
				System.err.println("query: "+character);
				System.err.println("candidate in alphabet: "+alphabet[-1-i]);
				System.exit(1);
			}
			bw.write(i+"");
		}
	}
	
	
	/**
	 * Same logic as $translateRead_unique()$, but taking orientation into account.
	 */
	private static final void translate_periodic(Block block, BufferedWriter bw) throws IOException {
		final int CAPACITY = 100;  // Arbitrary
		boolean orientation;
		int i, j;
		int repeat, length, last;
		final int lastCharacter = block.lastCharacter;
		Character character;
		
		// Collecting characters
		if (stack==null) stack = new int[CAPACITY];
		last=-1;
		for (j=0; j<=lastCharacter; j++) {
			character=block.characters[j];
			repeat=character.repeat; orientation=character.orientation; length=character.length;
			i=Arrays.binarySearch(alphabet,lastUnique+1,lastPeriodic+1,character);
			if (character.isOpen()) {
				if (i<0) i=-1-i;
				if (alphabet[i].repeat!=repeat || alphabet[i].orientation!=orientation) i--;
				while (i>lastUnique && alphabet[i].repeat==repeat && alphabet[i].orientation==orientation && alphabet[i].length>=length) i--;
				i++;
				if (alphabet[i].repeat!=repeat || alphabet[i].orientation!=orientation) {
					System.err.println("translateRead_periodic> ERROR: open periodic repeat not found in the alphabet");
					System.err.println("query: "+character);
					System.err.println("first candidate in alphabet: "+alphabet[i]);
					System.exit(1);
				}
				last=appendToStack(-1-i,last);
			}
			else {
				if (i<0) {
					System.err.println("translateRead_periodic> ERROR: closed periodic character not found in the alphabet:");
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
					if (alphabet[j].repeat!=repeat || alphabet[j].orientation!=orientation || alphabet[j].implies_tooFarAfter(character)) break;
					if (alphabet[j].implies(character,quantum) || alphabet[j].isSimilar(character)) {
					   found=true;
					   last=appendToStack(j,last);
					}
				}
				for (j=i-1; j>lastPeriodic; j--) {
					if (alphabet[j].repeat!=repeat || alphabet[j].orientation!=orientation) break;
					if (alphabet[j].implies(character,quantum) || alphabet[j].isSimilar(character)) {
						found=true;
						last=appendToStack(j,last);
					}
				}
				if (!found) {
					System.err.println("translateRead> ERROR: open nonperiodic character not found in the alphabet:");
					System.err.println("query: "+character);
					System.err.println("candidate in alphabet: "+alphabet[Math.min(i,lastAlphabet)]);
					System.exit(1);
				}
			}
			else {
				if (i<0) {
					System.err.println("translateRead> ERROR: closed nonperiodic character not found in the alphabet:");
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
	 *
	 * @param characterCount one cell per character of the alphabet, plus one; assumed to 
	 * be all zeros.
	 */
	public static final void incrementCharacterCounts(String str, long[] characterCount) {
		int i, j, k;
		int to, c;
		String[] tokens;
		String[][] blocks;
		
		// Building $characterCount$.
		if (str.length()==0) return;
		tokens=str.split(SEPARATOR_MAJOR+"");
		blocks = new String[tokens.length][0];
		for (i=0; i<=tokens.length-1; i++) {
			if (tokens[i].indexOf(SEPARATOR_MINOR+"")>=0) blocks[i]=tokens[i].split(SEPARATOR_MINOR+"");
			else {
				blocks[i] = new String[1];
				blocks[i][0]=tokens[i];
			}
		}
		for (i=0; i<blocks.length; i++) {
			for (j=0; j<blocks[i].length; j++) {
				c=Integer.parseInt(blocks[i][j]);
				if (c<0) {
					c=-1-c;
					to=c<=lastUnique?lastUnique:lastPeriodic;
					for (k=c; k<=to; k++) characterCount[k]++;
				}
				else characterCount[c]++;
			}
		}
	}
	
	
	/**
	 * Adds to $kmers$ every k-mer of the translated read $str$. If a block contains 
	 * multiple characters, every character is used to build a distinct kmer ("or").
     *
	 * @param sb temporary space, assumed to be empty.
	 */
	public static final void getKmers(String str, int k, Hashtable<String,Long> kmers, StringBuilder sb) {
		int i;
		int nBlocks;
		String[] tokens;
		String[][] blocks;
		
		if (str.length()<(k<<1)-1) return;
		i=-1; nBlocks=0;
		while (true) {
			i=str.indexOf(SEPARATOR_MAJOR+"",i+1);
			if (i<0) break;
			else nBlocks++;
		}
		nBlocks++;
		if (nBlocks<k) return;
		tokens=str.split(SEPARATOR_MAJOR+"");
		blocks = new String[nBlocks][0];
		for (i=0; i<=tokens.length-1; i++) {
			if (tokens[i].indexOf(SEPARATOR_MINOR+"")>=0) blocks[i]=tokens[i].split(SEPARATOR_MINOR+"");
			else {
				blocks[i] = new String[1];
				blocks[i][0]=tokens[i];
			}
		}
		if (stack==null || stack.length<(k+1)*3) stack = new int[(k+1)*3];
		for (i=0; i<=nBlocks-k; i++) loadTranslatedRead_impl(blocks,i,k,kmers,sb);
	}
	
	
	/**
	 * Adds $blocks[first..first+k-1]$ to $kmers$.
	 *
	 * @param sb temporary space.
	 */
	private static final void loadTranslatedRead_impl(String[][] blocks, int first, int k, Hashtable<String,Long> kmers, StringBuilder sb) {
		int i;
		int top, row, column, lastChild, length;
		Long value;
		String key;
		
		top=-1;
		stack[++top]=-1; stack[++top]=-1; stack[++top]=first-1;
		while (top>=0) {
			row=stack[top]; column=stack[top-1]; lastChild=stack[top-2];
			if (row==first+k-1) {
				key=sb.toString(); value=kmers.get(key);
				if (value==null) kmers.put(key,Long.valueOf(1));
				else kmers.put(key,Long.valueOf(value.longValue()+1));
			}
			if (row==first+k-1 || lastChild==blocks[row+1].length-1) {
				top-=3;
				i=sb.length()-1;
				while (i>=0 && sb.charAt(i)!=SEPARATOR_MAJOR) i--;
				i++;
				sb.delete(i,sb.length());
			}
			else {
				lastChild++;
				stack[top-2]=lastChild;
				stack[++top]=-1; stack[++top]=lastChild; stack[++top]=row+1;
				sb.append((row<first?"":(SEPARATOR_MAJOR+""))+blocks[row+1][lastChild]);
			}
		}
	}
	
	
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
		 * on a canonical (and arbitrary) order.
		 *
		 * Remark: characters might not be sorted after the procedure completes.
		 *
		 * @param tmpCharacter temporary space.
		 */
		public final void addCharacter(AlignmentRow alignment, int distanceThreshold, Character tmpCharacter) {
			int i;
			int found;
			
			if (alignment.orientation) {
				tmpCharacter.openStart=alignment.startA<=distanceThreshold;
				tmpCharacter.openEnd=alignment.endA>=Reads.getReadLength(alignment.readA)-distanceThreshold;
			}
			else {
				tmpCharacter.openEnd=alignment.startA<=distanceThreshold;
				tmpCharacter.openStart=alignment.endA>=Reads.getReadLength(alignment.readA)-distanceThreshold;
			}
			tmpCharacter.repeat=alignment.readB; tmpCharacter.orientation=alignment.orientation;
			if (!isPeriodic[alignment.readB]) {
				tmpCharacter.start=alignment.startB;
				tmpCharacter.end=alignment.endB;
				tmpCharacter.length=0;
			}
			else {
				tmpCharacter.start=-1; tmpCharacter.end=-1;
				tmpCharacter.length=alignment.endA-alignment.startA+1; // A is correct here
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
				else if ( alignment.startB<characters[found].start || 
						  (alignment.startB==characters[found].start && alignment.endB<characters[found].end)
			        	) {
					characters[found].start=alignment.startB;
					characters[found].end=alignment.endB;
				}
				if (tmpCharacter.openStart) characters[found].openStart=true;
				if (tmpCharacter.openEnd) characters[found].openEnd=true;
			}
		}
		
		
		/**
		 * (1) Sorts $characters$. (2) Makes open marks homogeneous across all characters.
		 * (3) If the block contains a periodic repeat, removes all non-periodic repeats 
		 * from the block and sets the length of every periodic repeat to
		 * $periodicLength$.
		 */
		public final void cleanCharacters(int periodicLength) {
			boolean foundPeriodic, foundNonperiodic, openLeft, openRight;
			int i, j;
			Character tmpCharacter;
			
			// Enforcing periodic/nonperiodic criteria
			foundPeriodic=false; foundNonperiodic=false;
			for (i=0; i<=lastCharacter; i++) {
				if (isPeriodic[characters[i].repeat]) {
					foundPeriodic=true;
					characters[i].length=periodicLength;
				}
				else foundNonperiodic=true;
			}
			if (!foundPeriodic) {
				Arrays.sort(characters,0,lastCharacter+1);
				return;
			}
			if (foundNonperiodic) {
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
			
			// Enforcing open criteria
			openLeft=false; openRight=false;
			for (i=0; i<=lastCharacter; i++) {
				if ((characters[i].orientation && characters[i].openStart) ||  (!characters[i].orientation && characters[i].openEnd)) openLeft=true;
				if ((characters[i].orientation && characters[i].openEnd) ||  (!characters[i].orientation && characters[i].openStart)) openRight=true;				
			}
			for (i=0; i<=lastCharacter; i++) {
				if (characters[i].orientation) { characters[i].openStart=openLeft; characters[i].openEnd=openRight; }
				else { characters[i].openStart=openRight; characters[i].openEnd=openLeft; }
			}			
			
			Arrays.sort(characters,0,lastCharacter+1);
		}
		
		
		public final boolean openStart() { return characters[0].openStart; }
		
		
		public final boolean openEnd() { return characters[0].openEnd; }
	}
	
	
	/**
	 * A substring of a repeat in a specific orientation and with open/closed endpoints.
	 */
	public static class Character implements Comparable {
		public int id;
		public int repeat;
		public boolean orientation;
		public int start, end;
		public boolean openStart, openEnd;  // The repeat might continue beyond start/end
		public int length;  // >0: unique or periodic; 0: repetitive non-periodic.
		
		public Character() { reset(); }
		
		public void reset() { 
			id=-1;
			repeat=-1; orientation=false; start=-1; end=-1; length=-1;
			openStart=true; openEnd=true;
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
		 * Both characters are assumed to be quantized, and repetitive with the same 
		 * repeat. Unique characters are not assumed to imply one another, since fully-
		 * open and half-open unique characters are assumed to have been discarded.
		 */
		public final boolean implies(Character otherCharacter, int quantum) {
			if (start==-1) {  // Periodic
				if (length!=otherCharacter.length) return false;
				else if (length<otherCharacter.length) return false;
				else if ( (!otherCharacter.openStart && !otherCharacter.openEnd) ||
					      (!otherCharacter.openStart && openStart) ||
				          (!otherCharacter.openEnd && openEnd)
						) return false;
			}
			else {  // Nonperiodic
				if (otherCharacter.start==start && otherCharacter.end==end) {
					if (Math.abs(getLength(),repeatLengths[repeat])<=quantum) {
						if ( (!openStart && !openEnd && (otherCharacter.openStart || otherCharacter.openEnd)) ||
							 ((!openStart || !openEnd) && otherCharacter.openStart && otherCharacter.openEnd)
						   ) return true;
						else return false;
					}
					else return false;
				}
				else if (Intervals.isContained(otherCharacter.start,otherCharacter.end,start,end)) {
					if ( (start==otherCharacter.start && openStart && !otherCharacter.openStart) ||
						 (otherCharacter.start>start && !otherCharacter.openStart) ||
						 (end==otherCharacter.end && openEnd && !otherCharacter.openEnd) ||
						 (otherCharacter.end<end && !otherCharacter.openEnd)
					   ) return false;
				}
				else return false;
			}
			return true;
		}
		
		
		/**
		 * Assumes that the characters have the same repeats, and uses the same criteria 
		 * as $implies()$.
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
		
		
		/**
		 * Assumes that the characters have the same repeat and are quantized.
		 */
		public final boolean isSimilar(Character otherCharacter) {
			return start==otherCharacter.start && end==otherCharacter.end && length==otherCharacter.length && 
				   openStart==otherCharacter.openStart && openEnd==otherCharacter.openEnd;
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
// Either do like Arne, explicitly modeling this as smaller blocks.
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