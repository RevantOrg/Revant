package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.HashMap;

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
	private static boolean[] isBlockUnique, isBlockOpen;
	
	
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
	 * @param alignmentsFile alignments between reads (readA) and repeats (readB);
	 * @param maxError alignments with error rate greater than this are discarded;
	 * @param uniqueFile the procedure writes just $maxOpenLength_unique$ there.
	 */
	public static final void collectCharacterInstances(String alignmentsFile, double maxError, int distanceThreshold, int lengthThreshold, String outputFile, String uniqueFile) throws IOException {
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
		bw = new BufferedWriter(new FileWriter(uniqueFile));
		bw.write(maxOpenLength_unique+"\n");
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
		int firstJForNextI, inPeriodic;
		int startA, endA, newLength, component, nComponents;
		final int lengthA = Reads.getReadLength(alignments[0].readA);
		
		if (periodicIntervals.length<(lastAlignment+1)<<1) periodicIntervals = new int[(lastAlignment+1)<<1];
		if (points.length<(lastAlignment+1)<<1) points = new int[(lastAlignment+1)<<1];
		AlignmentRow.order=AlignmentRow.ORDER_STARTA;
		if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
		
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
		if (lastPoint>0) Arrays.sort(points,0,lastPoint+1);
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
	 * Quantizes the endpoints of every repeat and removes character instances that can be
	 * considered substrings of others.
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
	 * the alphabet. Reverse-complement pairs do not need to have the same open
	 * information.
	 */
	public static final void closeAlphabetByRC() {
		final int QUANTUM = 1000;  // Arbitrary
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
					if (!alphabet[j].orientation && alphabet[j].length==length) {
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
					newAlphabet[k] = new Character(repeat,false,-1,-1,length,alphabet[i].openEnd,alphabet[i].openStart);
					nAdded++;
				}
				else alphabet[found].flag=1;
			}
			else {
				found=-1;
				for (j=i-1; j>lastUnique; j--) {
					if (alphabet[j].repeat!=repeat) break;
					if (alphabet[j].orientation && alphabet[j].length==length) {
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
					newAlphabet[k] = new Character(repeat,true,-1,-1,length,alphabet[i].openEnd,alphabet[i].openStart);
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
					if (!alphabet[j].orientation && alphabet[j].start==start && alphabet[j].end==end) {
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
					newAlphabet[k] = new Character(repeat,false,start,end,0,alphabet[i].openEnd,alphabet[i].openStart);
					nAdded++;
				}
				else alphabet[found].flag=1;
			}
			else {
				found=-1;
				for (j=i-1; j>lastUnique; j--) {
					if (alphabet[j].repeat!=repeat) break;
					if (alphabet[j].orientation && alphabet[j].start==start && alphabet[j].end==end) {
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
					newAlphabet[k] = new Character(repeat,true,start,end,0,alphabet[i].openEnd,alphabet[i].openStart);
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
	 * @param outputFile_characters translation of every read (one per line); a read has
	 * an empty line if it contains no repeat, or if it is fully contained in a single
	 * repeat;
	 * @param outputFile_boundaries character boundaries (one read per line);
	 * @param outputFile_fullyUniqueReads list of IDs of reads that contain no repeat;
	 * @param outputFile_oneRepeatReads list of IDs of reads that are fully contained in a
	 * single repeat.
	 */
	public static final void translateReads(String alignmentsFile, int lastTranslatedRead, double maxError, int quantum, String outputFile_characters, String outputFile_boundaries, String outputFile_histogram, String outputFile_fullyUniqueReads, String outputFile_oneRepeatReads) throws IOException {
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
			if ((2.0*Alignments.diffs)/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2)>maxError) {
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
							bw1.newLine(); bw2.newLine(); histogram[0]++;
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
					bw1.newLine(); bw2.newLine(); histogram[0]++; 
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
				System.err.println("translateRead_unique> ERROR: closed unique character not found in the alphabet");
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
				while (k<i && !alphabet[k].implies(character,-1)) k++;
				i=k;
				if (alphabet[i].repeat!=repeat || alphabet[i].orientation!=orientation || !alphabet[i].implies(character,-1)) {
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
	 * Resets $characterCount[x]$ and $characterCount[y]$ to the sum of their values for 
	 * every $x,y$ that are the reverse-complement of each other.
	 *
	 * @param marked temporary space, of size at least $lastAlphabet+1$.
	 */
	public static final void symmetrizeCharacterCounts(long[] characterCount, boolean[] marked) {
		boolean orientation, openStart, openEnd, isHalfOpen;
		int i, j;
		int repeat, start, end, length;
		
		Math.set(marked,lastAlphabet,false);
		
		// Periodic
		for (i=lastUnique+1; i<lastPeriodic; i++) {
			if (marked[i]) continue;
			repeat=alphabet[i].repeat; orientation=alphabet[i].orientation;
			length=alphabet[i].length; 
			openStart=alphabet[i].openStart; openEnd=alphabet[i].openEnd;
			isHalfOpen=openStart!=openEnd;
			for (j=i+1; j<=lastPeriodic; j++) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].length==length && alphabet[j].orientation!=orientation) {
					if ( (isHalfOpen && alphabet[j].openStart!=openStart && alphabet[j].openEnd!=openEnd) ||
						 (!isHalfOpen && alphabet[j].openStart==openStart && alphabet[j].openEnd==openEnd)
					   ) {
						characterCount[i]+=characterCount[j];
						characterCount[j]=characterCount[i];
						marked[j]=true;
						break;
					}
				}
			}
		}
		
		// Nonperiodic
		for (i=lastPeriodic+1; i<lastAlphabet; i++) {
			if (marked[i]) continue;
			repeat=alphabet[i].repeat; orientation=alphabet[i].orientation;
			start=alphabet[i].start; end=alphabet[i].end;
			openStart=alphabet[i].openStart; openEnd=alphabet[i].openEnd;
			isHalfOpen=openStart!=openEnd;
			for (j=i+1; j<=lastAlphabet; j++) {
				if (alphabet[j].repeat!=repeat) break;
				if (marked[j]) continue;
				if (alphabet[j].start==start && alphabet[j].end==end && alphabet[j].orientation!=orientation) {
					if ( (isHalfOpen && alphabet[j].openStart!=openStart && alphabet[j].openEnd!=openEnd) ||
						 (!isHalfOpen && alphabet[j].openStart==openStart && alphabet[j].openEnd==openEnd)
					   ) {
						characterCount[i]+=characterCount[j];
						characterCount[j]=characterCount[i];
						marked[j]=true;
						break;
					}
				}
			}
		}
	}
	
	
	/**
	 * Builds a histogram of symmetrized $characterCount$ values. Rows: counts. Columns:
	 * 0: unique, open; 
	 * 1: unique, closed;
	 * 2: periodic, open;
	 * 3: periodic, closed; 
	 * 4: nonperiodic, open;
	 * 5: nonperiodic, closed.
	 *
	 * @param marked the same array used by $symmetrizeCharacterCounts()$, after that
	 * procedure completes.
	 */
	public static final long[][] getCharacterHistogram(long[] characterCount, boolean[] marked, int maxFrequency) throws IOException {
		int i;
		long count;
		long[][] characterHistogram;

		// Unique
		characterHistogram = new long[maxFrequency+1][3<<1];
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
			characterHistogram[(int)count][alphabet[i].isOpen()?2:3]++;
		}		
		// Nonperiodic
		for (i=lastPeriodic+1; i<=lastAlphabet; i++) {
			if (marked[i]) continue;
			count=characterCount[i];
			if (count>maxFrequency) count=maxFrequency;
			characterHistogram[(int)count][alphabet[i].isOpen()?4:5]++;
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
	 * negative ID into a list of positive IDs, and builds arrays $isBlock{Unique,Open}$.
	 *
	 * Remark: the procedure follows the logic of $incrementCharacterCounts()$.
	 */
	private static final void loadIntBlocks(int nBlocks) {
		boolean orientation, unique;
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
		if (isBlockOpen==null || isBlockOpen.length<nBlocks) isBlockOpen = new boolean[nBlocks];
		Math.set(isBlockOpen,nBlocks-1,false);
		
		// Building arrays
		isBlockOpen[0]=true; isBlockOpen[nBlocks-1]=true;
		for (i=0; i<nBlocks; i++) {
			last=lastInBlock[i];
			nElements=0;
			for (j=0; j<=last; j++) {
				value=Integer.parseInt(blocks[i][j]);
				if (value>=0) nElements++;
				else {
					value=-1-value;
					if (value<=lastUnique) nElements+=lastUnique+1-value;
					else {
						repeat=alphabet[value].repeat; orientation=alphabet[value].orientation;
						for (k=value; k<=lastPeriodic; k++) {
							if (alphabet[k].repeat!=repeat || alphabet[k].orientation!=orientation) break;
							if (alphabet[k].implies(alphabet[value],-1)) nElements++;
						}
					}
				}
			}
			if (intBlocks[i].length<nElements) intBlocks[i] = new int[nElements];
			k=-1; unique=false;
			for (j=0; j<=last; j++) {
				value=Integer.parseInt(blocks[i][j]);
				if (value>=0) {
					intBlocks[i][++k]=value;
					if (value==lastAlphabet+1 || alphabet[value].isOpen()) isBlockOpen[i]=true;
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
					isBlockOpen[i]=true;
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
	 * @param read2boundaries boundaries of the characters in $read2characters$;
	 * @param tmpChar temporary space.
	 */
	public static final void cleanTranslatedRead_collectCharacterInstances(String read2characters, String read2boundaries, int readLength, int minCount, int quantum, BufferedWriter bw, Character tmpChar) throws IOException {
		int i, j, k;
		int c, length, first, nBlocks, nBoundaries;
		String[] tokens;
		
		if (read2characters.length()==0) return;
		if (read2boundaries.indexOf(SEPARATOR_MINOR+"")>=0) {
			tokens=read2boundaries.split(SEPARATOR_MINOR+"");
			nBoundaries=tokens.length;
			if (boundaries==null || boundaries.length<nBoundaries) boundaries = new int[nBoundaries];
			for (i=0; i<nBoundaries; i++) boundaries[i]=Integer.parseInt(tokens[i]);
		}
		else {
			nBoundaries=1;
			boundaries[0]=Integer.parseInt(read2boundaries);
		}
		nBlocks=loadBlocks(read2characters);
		removeRareCharacters(nBlocks,minCount,lastUnique,lastPeriodic,lastAlphabet);
		first=-1;
		for (i=0; i<nBlocks; i++) {
			if (lastInBlock[i]!=-1) {
				if (first!=-1) {
					tmpChar.repeat=UNIQUE;
					tmpChar.orientation=false;
					tmpChar.start=-1; tmpChar.end=-1;
					length=(i==nBlocks-1?readLength:boundaries[i])-(first==0?0:boundaries[first-1]);
					tmpChar.length=length;
					tmpChar.openStart=first==0;
					tmpChar.openEnd=i==nBlocks-1;
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
			tmpChar.orientation=false;
			tmpChar.start=-1; tmpChar.end=-1;
			length=readLength-(first==0?0:boundaries[first-1]);
			tmpChar.length=length;
			tmpChar.openStart=first==0; 
			tmpChar.openEnd=true;
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
	 * Removes from $blocks$ all characters with count smaller than $minCount$.
	 *
	 * Remark: the procedure assumes that global variable $alphabetCount$ has already been
	 * initialized.
	 */
	private static final void removeRareCharacters(int nBlocks, int minCount, int lastUnique, int lastPeriodic, int lastAlphabet) {
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
						if (alphabetCount[c]<minCount) {
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
				else if (alphabetCount[c]<minCount) deleted=true;
				if (!deleted) blocks[i][++k]=blocks[i][j];
			}
			lastInBlock[i]=k;
		}
	}
	
	
	/**
	 * Removes from $alphabet$ all characters with $alphabetCount < minCount$, and adds to
	 * $alphabet$ all new unique characters in $newCharactersFile$.
	 *
	 * Remark: $alphabetCount$ is not valid after the procedure completes.
	 *
	 * @return an array that maps every non-unique character in the original $alphabet$ to
	 * its position in the new $alphabet$. Positions are relative to the corresponding
	 * values of $lastUnique+1$.
	 */
	public static final int[] cleanTranslatedRead_updateAlphabet(int nNewCharacters, String newCharactersFile, int minCount) throws IOException {
		int i, j;
		String str;
		Character tmpChar;
		BufferedReader br;
		int[] old2new;
		Character[] newAlphabet;
		
		// Removing rare characters
		old2new = new int[lastAlphabet-lastUnique];
		Math.set(old2new,old2new.length-1,-1);
		j=lastUnique;
		for (i=lastUnique+1; i<=lastAlphabet; i++) {
			if (alphabetCount[i]<minCount) continue;
			j++;
			tmpChar=alphabet[j];
			alphabet[j]=alphabet[i];
			alphabet[i]=tmpChar;
			if (alphabet[j].start==-1) lastPeriodic=j;
			old2new[i-lastUnique-1]=j-lastUnique-1;
		}
		lastAlphabet=j;
		
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
	 * @param read2characters_old old translation;
	 * @param read2boundaries_old old block boundaries of the translation;
	 * @param newAlphabet obtained from the old alphabet by running $cleanTranslatedRead_
	 * updateAlphabet()$;
	 * @param old2new the output of $cleanTranslatedRead_updateAlphabet()$;
	 * @param read2characters_new,read2boundaries_new output files;
	 * @return the number of blocks in the new translation.
	 */
	public static final int cleanTranslatedRead_updateTranslation(String read2characters_old, String read2boundaries_old, Character[] oldAlphabet, int lastUnique_old, int lastPeriodic_old, int lastAlphabet_old, Character[] newAlphabet, int lastUnique_new, int lastPeriodic_new, int lastAlphabet_new, int[] old2new, int readLength, int minCount, int quantum, BufferedWriter read2characters_new, BufferedWriter read2boundaries_new, Character tmpChar) throws IOException {
		int i, j;
		int c, length, first, last, nBlocks, nBoundaries, nAppendedBlocks;
		String[] tokens;
		
		if (read2characters_old.length()==0) return 0;
		if (read2boundaries_old.indexOf(SEPARATOR_MINOR+"")>=0) {
			tokens=read2boundaries_old.split(SEPARATOR_MINOR+"");
			nBoundaries=tokens.length;
			if (boundaries==null || boundaries.length<nBoundaries) boundaries = new int[nBoundaries];
			for (i=0; i<nBoundaries; i++) boundaries[i]=Integer.parseInt(tokens[i]);
		}
		else {
			nBoundaries=1;
			boundaries[0]=Integer.parseInt(read2boundaries_old);
		}
		nBlocks=loadBlocks(read2characters_old);
		alphabet=oldAlphabet; lastUnique=lastUnique_old; lastPeriodic=lastPeriodic_old; lastAlphabet=lastAlphabet_old;
		removeRareCharacters(nBlocks,minCount,lastUnique_old,lastPeriodic_old,lastAlphabet_old);
		alphabet=newAlphabet; lastUnique=lastUnique_new; lastPeriodic=lastPeriodic_new; lastAlphabet=lastAlphabet_new;
		first=-1; nAppendedBlocks=0;
		for (i=0; i<nBlocks; i++) {
			if (lastInBlock[i]!=-1) {
				if (first!=-1) {
					tmpChar.repeat=UNIQUE;
					tmpChar.orientation=false;
					tmpChar.start=-1; tmpChar.end=-1;
					length=(i==nBlocks-1?readLength:boundaries[i])-(first==0?0:boundaries[first-1]);
					tmpChar.length=length;
					tmpChar.openStart=first==0;
					tmpChar.openEnd=i==nBlocks-1;
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
						tmpChar.orientation=false;
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
				tmpChar.orientation=false;
				tmpChar.start=-1; tmpChar.end=-1;
				length=readLength-(first==0?0:boundaries[first-1]);
				tmpChar.length=length;
				tmpChar.openStart=first==0;
				tmpChar.openEnd=true;
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
				// Reads with a single unique block are recoded to have no block
				return 0;
			}
		}
		return nAppendedBlocks;
	}

	

	
	// ------------------------------ K-MER PROCEDURES -----------------------------------
	
	/**
	 * If $newKmers$ is not null, the procedure uses every length-k window that satisfies 
	 * the conditions in $uniqueMode,openMode,multiMode$ (see procedure $isValidKmer()$ 
	 * for details) to add a k-mer to $newKmers$. If a block contains multiple characters, 
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
	 * Remark: 1-mers collected by this procedure might have a different (and even 
	 * smaller) count than the one produced by the $getCharacterHistogram()$ pipeline, 
	 * because of the constraints and of canonization. E.g. a 1-mer might be rare in this 
	 * procedure but frequent in $getCharacterHistogram()$, if the latter counted all the
	 * occurrences of the character but the former counts just closed occurrences.
	 *
	 * Remark: often several characters map to the same block in practice, and most of 
	 * the k-mers (k>=2) that result from this are rare.
	 *
	 * Remark: the procedure uses global variables $blocks,intBlocks,stack$.
	 * 
	 * @param haplotypeCoverage coverage of one haplotype;
	 * @param tmpKmer temporary space;
	 * @param tmpArray2 temporary space, of size at least k;
	 * @param tmpArray3 temporary space, of size at least 2k.
	 */
	public static final int getKmers(String str, int k, int uniqueMode, int openMode, int multiMode, HashMap<Kmer,Kmer> newKmers, HashMap<Kmer,Kmer> oldKmers, int[] avoidedIntervals, int lastAvoidedInterval, int haplotypeCoverage, Kmer tmpKmer, int[] tmpArray2, int[] tmpArray3) {
		int i, j;
		int nBlocks, sum, start, end, nHaplotypes, out;
		
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
		for (i=0; i<nBlocks; i++) sum+=lastInBlock[i]+1;
		if (stack==null || stack.length<sum*3) stack = new int[sum*3];
		
		// Loading k-mers
		j=0;
		for (i=0; i<=nBlocks-k; i++) {
			while (j<lastAvoidedInterval && avoidedIntervals[j]<i) j+=3;
			if (j<lastAvoidedInterval && avoidedIntervals[j]+avoidedIntervals[j+1]-1<=i+k-1) continue;
			if (!isValidWindow(i,k,uniqueMode,openMode,multiMode)) continue;
			nHaplotypes=getKmers_impl(i,k,newKmers,oldKmers,haplotypeCoverage,tmpKmer,stack,tmpArray2,tmpArray3);
			if (newKmers==null && nHaplotypes!=-1) { avoidedIntervals[++out]=i; avoidedIntervals[++out]=k; avoidedIntervals[++out]=nHaplotypes; }
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
	 * @param openMode open blocks:
	 * 0: are allowed; 
	 * 1: are not allowed;
	 * @param multiMode blocks with multiple characters:
	 * 0: are allowed; 
	 * 1: are not allowed if they are open (open blocks are more likely to match several 
	 *    characters because they contain just a fraction of a character);
	 * 2: are not allowed.
	 */
	private static final boolean isValidWindow(int first, int k, int uniqueMode, int openMode, int multiMode) {
		int i;
		
		if (uniqueMode==1) {
			if (isBlockUnique[first] || isBlockUnique[first+k-1]) return false;
		}
		else if (uniqueMode==2) {
			for (i=0; i<=k-1; i++) {
				if (isBlockUnique[first+i]) return false;
			}
		}
		if (openMode==1) {
			for (i=0; i<=k-1; i++) {
				if (isBlockOpen[first+i]) return false;
			}
		}
		if (multiMode==1) {
			for (i=0; i<=k-1; i++) {
				if (lastInBlock_int[first+i]>0 && isBlockOpen[first+i]) return false;
			}
		}
		else if (multiMode==2) {
			for (i=0; i<=k-1; i++) {
				if (lastInBlock_int[first+i]>0) return false;
			}
		}
		return true;
	}
	
	
	/**
	 * If $newKmers$ is not null, adds every possible canonized instance of 
	 * $intBlocks[first..first+k-1]$ to $newKmers$. Otherwise, checks whether one of the 
	 * possible canonized instances of $intBlocks[first..first+k-1]$ belongs to
	 * $oldKmers$, and if so returns an estimate of the number of haplotypes in which it
	 * occurs (returns -1 otherwise).
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
				key.set(tmpArray2,0,k); key.canonize(k,tmpArray3);
				if (newKmers!=null) {
					value=newKmers.get(key);
					if (value==null) {
						newKey = new Kmer(key,k); newKey.count=1;
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
		
		public Kmer() { }
		
		public Kmer(Kmer otherKmer, int k) {
			sequence = new int[k];
			System.arraycopy(otherKmer.sequence,0,sequence,0,k);
			count=0;
		}
		
		public Kmer(String str, int k) {
			sequence = new int[k];
			String[] tokens = str.split(",");
			for (int i=0; i<k; i++) sequence[i]=Integer.parseInt(tokens[i]);
			count=Long.parseLong(tokens[k]);
		}
		
		public void set(int[] fromArray, int first, int k) {
			sequence = new int[k];
			System.arraycopy(fromArray,first,sequence,0,k);
			count=0;
		}
		
		/**
		 * Resets $sequence$ to the lexicographically smaller between $sequence$ with
		 * every character canonized, and the reverse of $sequence$ with the complement of 
		 * every character canonized.
		 * 
		 * @tmpArray temporary space of size at least 2k.
		 */
		public void canonize(int k, int[] tmpArray) {
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
	 * Tuples (position,length,nHaplotypes). Could be an array of bytes, since there are
	 * fewer than 128 blocks per read in practice: we keep it as an array of 32-bit 
	 * integers to be fully generic.
	 */
	private static int[][] blueIntervals;
	private static int blueIntervals_last;
	
	/**
	 * Sorted lists of read IDs: with a blue interval; fully non-repetitive; fully 
	 * contained in a single repeat block; translated into several blocks.
	 */
	private static int[] blueIntervals_reads;
	private static int[] fullyUnique, fullyContained, translated;
	
	/**
	 * For every read in $translated$: block boundaries; flags marking non-repetitive
	 * blocks; full translation.
	 */ 
	private static int[][] boundaries_all;
	private static byte[][] isBlockUnique_all;
	private static int[][][] translation_all;
	
	
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
	 * Loads the entire boundaries file in $boundaries_all$, and stores in 
	 * $isBlockUnique_all$ a bitmap for every block of every read.
	 *
	 * @param loadIsBlockUnique loads the bitvector $isBlockUnique_all$;
	 * @param loadTranslation loads in $translation_all$ the translation of every read 
	 * that has one.
	 */
	public static final void loadAllBoundaries(String translatedFile, boolean loadIsBlockUnique, boolean loadTranslation, String boundariesFile) throws IOException {
		int i, j, r;
		int read, cell, mask, nBlocks, nBytes, nReads, nBoundaries;
		String str;
		BufferedReader br;
		String[] tokens;
		
		// Computing the set of translated reads
		br = new BufferedReader(new FileReader(translatedFile));
		str=br.readLine();
		nReads=0;
		while (str!=null) {
			if (str.length()!=0) nReads++;
			str=br.readLine();
		}
		br.close();
		translated = new int[nReads];
		br = new BufferedReader(new FileReader(translatedFile));
		str=br.readLine(); r=0; i=-1;
		while (str!=null) {
			if (str.length()!=0) translated[++i]=r;
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
				if (str.length()==0) {
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
						if (isBlockUnique[i*8+j]) cell|=mask;
						mask<<=1;
					}
					isBlockUnique_all[r][i]=(byte)mask;
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
				if (str.length()==0) {
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
		
		br = new BufferedReader(new FileReader(intervalsFile));
		str=br.readLine(); read=-1; blueIntervals_last=-1;
		while (str!=null) {
			read++;
			if (str.length()==0) {
				str=br.readLine();
				continue;
			}
			blueIntervals_last++;
			if (blueIntervals_last>blueIntervals.length) {
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
	 * to a repeat, or to a sequence of repeats, that is likely to occur multiple times in 
	 * the genome. The intervals of the alignment in the two reads might cover mismatching 
	 * sequences of boundaries and different characters, but we can safely discard the 
	 * alignment anyway, since it just encodes a similarity between substrings of 
	 * (possibly different) repeats.
	 *
	 * Remark: the procedure does not need $alphabet$, but it needs the following arrays: 
	 * isBlockUnique_all | fullyUnique, fullyContained, boundaries_all, blueIntervals, blueIntervals_reads. 
	 *
	 * @param alignmentsFile output of LAshow, assumed to be sorted by readA;
	 * @param minIntersection min. length of a non-repetitive substring of the alignment,
	 * for the alignment not to be considered red.
	 */
	public static final void filterAlignments_loose(String alignmentsFile, String outputFile, int minIntersection) throws IOException {
		boolean isRepetitive;
		int p;
		int row, readA, readB;
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
			if (row%100000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			// Processing readA
			readA=Alignments.readA-1;
			while (lastFullyUnique<nFullyUnique && fullyUnique[lastFullyUnique]<readA) lastFullyUnique++;
			if (lastFullyUnique<nFullyUnique && fullyUnique[lastFullyUnique]==readA) {
				bw.write("1\n"); str=br.readLine(); row++;
				continue;
			}
			isRepetitive=false;
			while (lastFullyContained<nFullyContained && fullyContained[lastFullyContained]<readA) lastFullyContained++;
			if (lastFullyContained<nFullyContained && fullyContained[lastFullyContained]==readA) isRepetitive=true;
			else {
				while (lastTranslated<nTranslated && translated[lastTranslated]<readA) lastTranslated++;
				while (lastBlueInterval<=blueIntervals_last && blueIntervals_reads[lastBlueInterval]<readA) lastBlueInterval++;
				if (inRedRegion(readA,Alignments.startA,Alignments.endA,lastTranslated,lastBlueInterval<=blueIntervals_last&&blueIntervals_reads[lastBlueInterval]==readA?lastBlueInterval:-1,-1,Reads.getReadLength(readA),minIntersection)) isRepetitive=true;
			}
			if (!isRepetitive) {
				bw.write("1\n"); str=br.readLine(); row++;
				continue;
			}
			// Processing readB
			readB=Alignments.readB-1;
			if (readInArray(readB,fullyUnique,nFullyUnique-1,lastFullyUnique)>=0) {
				bw.write("1\n"); str=br.readLine(); row++;
				continue;
			}
			else if (readInArray(readB,fullyContained,nFullyContained-1,lastFullyContained)>=0) {
				bw.write("0\n"); str=br.readLine(); row++;
				continue;
			}
			p=readInArray(readB,translated,nTranslated-1,lastTranslated);
			isRepetitive=inRedRegion(readB,Alignments.startB,Alignments.endB,p,-2,lastBlueInterval,Reads.getReadLength(readB),minIntersection);			
			bw.write(isRepetitive?"0\n":"1\n"); 
			str=br.readLine(); row++;
		}
		br.close(); bw.close();
	}
	
	
	/**
	 * Tells whether interval $readID[intervalStart..intervalEnd]$ fully belongs to a 
	 * single repeat character, or to a sequence of repeat characters, that is likely to 
	 * occur multiple times in the genome.
	 *
	 * Remark: the procedure needs the following arrays: 
	 * isBlockUnique_all, boundaries_all, blueIntervals, blueIntervals_reads.
	 * 
	 * @param readID assumed to contain more than one block;
	 * @param boundariesAllID position of the read in $boundaries_all,isBlockUnique_all,
	 * tranlsation_all$;
	 * @param blueIntervalsID position of the read in $blueIntervals$; or -1 if the 
	 * read does not occur in $blueIntervals$; or -2 if the ID of the read in 
	 * $blueIntervals$ is unknown;
	 * @param blueIntervalsStart a position in $blueIntervals$ from which to start
	 * the search when $blueIntervalsID=-2$;
	 * @param minIntersection min. length of a non-repetitive substring of the alignment,
	 * for the alignment to be considered non-repetitive.
	 */
	private static final boolean inRedRegion(int readID, int intervalStart, int intervalEnd, int boundariesAllID, int blueIntervalsID, int blueIntervalsStart, int readLength, int minIntersection) {
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
				      Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection
				   )
				 )
			   ) return false;
			mask<<=1;
			for (j=1; j<8; j++) {
				if (j==nBlocks-1) break;
				blockStart=boundaries_all[boundariesAllID][j-1];
				blockEnd=boundaries_all[boundariesAllID][j];
				if ( (cell&mask)!=0 && 
					 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) || 
					   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
					   ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
						  Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection
					   )
					 ) 
				   ) return false;
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
						  Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection
					   )
					 )
				   ) return false;
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
					  Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection
				   )
				 )
			   ) return false;
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
				   ) return false;
			}
		}
		
		return true;
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
	 * Writes to $outputFile$ a one for every alignment of $alignmentsFile$ that contains
	 * a unique region in both readA and readB. If $mode=TRUE$, the procedure additionally
	 * requires that the intervals of the alignment in the two reads cover matching 
	 * sequences of boundaries and matching characters.
	 *
	 * Remark: the procedure needs the following arrays: 
	 * translation_all, alphabet | fullyUnique, fullyContained, boundaries_all, blueIntervals, blueIntervals_reads.
	 */
	public static final void filterAlignments_tight(String alignmentsFile, String outputFile, boolean mode, int minIntersection) throws IOException {
		boolean isFullyUniqueA;
		int p;
		int row, readA, readB, startA, endA, startB, endB;
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
			if (row%100000==0) System.err.println("Processed "+row+" alignments");
			Alignments.readAlignmentFile(str);
			// Processing readA
			readA=Alignments.readA-1; readB=Alignments.readB-1;
			startA=Alignments.startA; endA=Alignments.endA;
			startB=Alignments.startB; endB=Alignments.endB;
			while (lastFullyUnique<nFullyUnique && fullyUnique[lastFullyUnique]<readA) lastFullyUnique++;
			if (lastFullyUnique<nFullyUnique && fullyUnique[lastFullyUnique]==readA) {
				isFullyUniqueA=true;
				if (mode) {
					if (readInArray(readB,fullyUnique,nFullyUnique-1,lastFullyUnique)>=0) bw.write("1\n");
					else bw.write("0\n");
					str=br.readLine(); row++;
					continue;
				}
			}
			else {
				isFullyUniqueA=false;
				while (lastFullyContained<nFullyContained && fullyContained[lastFullyContained]<readA) lastFullyContained++;
				if (lastFullyContained<nFullyContained && fullyContained[lastFullyContained]==readA) {
					bw.write("0\n"); str=br.readLine(); row++;
					continue;
				}
				while (lastTranslated<nTranslated && translated[lastTranslated]<readA) lastTranslated++;
				while (lastBlueInterval<=blueIntervals_last && blueIntervals_reads[lastBlueInterval]<readA) lastBlueInterval++;
				if (!inBlueRegion(readA,startA,endA,lastTranslated,lastBlueInterval<=blueIntervals_last&&blueIntervals_reads[lastBlueInterval]==readA?lastBlueInterval:-1,-1,Reads.getReadLength(readA),minIntersection)) {
					bw.write("0\n"); str=br.readLine(); row++;
					continue;
				}
			}
			// Processing readB
			if (readInArray(readB,fullyUnique,nFullyUnique-1,lastFullyUnique)>=0) {
				if (mode) bw.write(isFullyUniqueA?"1\n":"0\n");
				else bw.write("1\n");
				str=br.readLine(); row++;
				continue;
			}
			else if (readInArray(readB,fullyContained,nFullyContained-1,lastFullyContained)>=0) {
				bw.write("0\n"); str=br.readLine(); row++;
				continue;
			}
			p=readInArray(readB,translated,nTranslated-1,lastTranslated);
			if (inBlueRegion(readB,startB,endB,p,-2,lastBlueInterval,Reads.getReadLength(readB),minIntersection)) {
				if (mode) {
					if (sameFactorization(readA,startA,endA,readB,startB,endB,Alignments.orientation)) bw.write("1\n");
					else bw.write("0\n");
				}
				else bw.write("1\n");
			}
			else bw.write("0\n"); 
			str=br.readLine(); row++;
		}
		br.close(); bw.close();
	}
	
	
	/**
	 * The dual of $inRedRegion()$: tells whether interval $readID[intervalStart..
	 * intervalEnd]$ fully belongs to a non-repetitive region, or straddles a non-
	 * repetitive region, or contains a sequence of repeat characters that is likely to 
	 * occur just once in the genome.
	 *
	 * Remark: the procedure needs the following arrays: 
	 * boundaries_all, translation_all, blueIntervals, blueIntervals_reads.
	 */
	private static final boolean inBlueRegion(int readID, int intervalStart, int intervalEnd, int boundariesAllID, int blueIntervalsID, int blueIntervalsStart, int readLength, int minIntersection) {
		int i, j;
		int start, end, blockStart, blockEnd, firstBlock, lastBlock;
		final int nBlocks = boundaries_all[boundariesAllID].length+1;
		final int nBytes = Math.ceil(nBlocks,8);
		
		// Checking the nonrepetitive blocks of the read, if any.
		blockStart=0; blockEnd=boundaries_all[boundariesAllID][0];
		if ( translation_all[boundariesAllID][0].length==1 && (translation_all[boundariesAllID][0][0]<=lastUnique || translation_all[boundariesAllID][0][0]==lastAlphabet+1) &&
			 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
			   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
			   ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
			      Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection
			   )
			 )
		   ) return true;
		for (i=1; i<nBlocks-1; i++) {
			blockStart=boundaries_all[boundariesAllID][i-1];
			blockEnd=boundaries_all[boundariesAllID][i];
			if ( translation_all[boundariesAllID][i].length==1 && (translation_all[boundariesAllID][i][0]<=lastUnique || translation_all[boundariesAllID][i][0]==lastAlphabet+1) &&
				 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
				   ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
				      Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection
				   )
				 )
			   ) return true;
		}
		blockStart=boundaries_all[boundariesAllID][nBlocks-2];
		blockEnd=readLength-1;
		if ( translation_all[boundariesAllID][nBlocks-1].length==1 && (translation_all[boundariesAllID][nBlocks-1][0]<=lastUnique || translation_all[boundariesAllID][nBlocks-1][0]==lastAlphabet+1) &&
			 ( Intervals.areApproximatelyIdentical(intervalStart,intervalEnd,blockStart,blockEnd) ||
			   Intervals.isApproximatelyContained(intervalStart,intervalEnd,blockStart,blockEnd) ||
			   ( !Intervals.isApproximatelyContained(blockStart,blockEnd,intervalStart,intervalEnd) &&
			      Intervals.intersectionLength(intervalStart,intervalEnd,blockStart,blockEnd)>=minIntersection
			   )
			 )
		   ) return true;
		
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
				   ) return true;
			}
		}
		
		return false;
	}
	
	
	/**
	 * @param read* index in $translated_all$;
	 * @return TRUE iff $readA[startA..endA]$ intersects the same number of blocks as
	 * $read[startB..endB]$, with boundaries at similar positions, and with at least one 
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
					 !Math.nonemptyIntersection(translation_all[readA][firstBlockA+i],0,translation_all[readA][firstBlockA+i].length-1,translation_all[readB][firstBlockB+i],0,translation_all[readB][firstBlockB+i].length-1)
				   ) return false;
			}
			if (!Math.nonemptyIntersection(translation_all[readA][lastBlockA],0,translation_all[readA][lastBlockA].length-1,translation_all[readB][lastBlockB],0,translation_all[readB][lastBlockB].length-1)) return false;
		}
		else {
			for (i=0; i<nBlocks-1; i++) {
				if ( Math.abs((boundaries_all[readA][firstBlockA+i]-startA)*ratio-(endB-boundaries_all[readB][lastBlockB-i-1]))>IDENTITY_THRESHOLD ||
					 !nonemptyIntersectionRC(readA,firstBlockA+i,readB,lastBlockB-i)
				   ) return false;
			}
			if (!nonemptyIntersectionRC(readA,lastBlockA,readB,firstBlockB)) return false;
		}
		return true;
	}
	
	
	/**
	 * Remark: the procedure uses global array $stack$ as temporary space.
	 *
	 * @return TRUE iff array $translation_all[readA][blockA]$ has at least one character 
	 * in common with the reverse-complemented characters of array 
	 * $translation_all[readB][blockB]$.
	 */
	private static final boolean nonemptyIntersectionRC(int readA, int blockA, int readB, int blockB) {
		int i;
		final int lengthB = translation_all[readB][blockB].length;
		
		if (stack==null || stack.length<lengthB) stack = new int[lengthB];
		for (i=0; i<lengthB; i++) stack[i]=canonizeCharacter(translation_all[readB][blockB][i],false);
		if (lengthB>1) Arrays.sort(stack,0,lengthB);
		return Math.nonemptyIntersection(translation_all[readA][blockA],0,translation_all[readA][blockA].length-1,stack,0,lengthB-1);
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
		public int repeat;
		public boolean orientation;
		public int start, end;
		public boolean openStart, openEnd;  // The repeat might continue beyond start/end
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