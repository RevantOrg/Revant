package de.mpi_cbg.revant.apps;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.imageio.*;
import java.text.*;

import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.Intervals;

/**
 * 
 *
 * @param args 
 */
public class FilterAlignments {
	
	/**
	 * Marks periodic repeats
	 */
	private static final boolean[] isPeriodic;
	
	
	
	private static final int DISTANCE_THRESHOLD = 20;  // Arbitrary
	private static final int LENGTH_THRESHOLD = 100;  // Arbitrary
	private static final String NON_REPETITIVE = "nonrepetitive";
	private static final String PAF_SEPARATOR = "\t";
	private static final int HEADER_ROWS = 100;  // in pixels
	private static int rFrom, gFrom, bFrom, rTo, gTo, bTo;
	private static Character[] alphabet;
	private static int lastAlphabet;
	private static int[] tmpArray;
	private static Character[] sequence;
	private static AlignmentRow[] alignments;
	
	
	private static int[] sequenceLengths;
	private static int lastSequenceLength;
	
	private static int K_VALUE;
	private static long idGenerator;
	
	/**
	 * Every distinct k-mer in a recoded read
	 */
	private static KmerTree kmerSet;
	
	private static int[] kmerCountHistogram, edgeCountHistogram, degreeHistogram;
	private static Kmer[] stack;
	
	
	private static final int ALPHABET_CAPACITY = 100000;  // Arbitrary
	private static final int SEQUENCE_CAPACITY = 1000000;  // Arbitrary
	private static final int ALIGNMENTS_CAPACITY = 100000;  // Arbitrary
	private static final int STACK_CAPACITY = 10000;
	
	
	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		K_VALUE=Integer.parseInt(args[1]);
		
		int i;
		int previous;
		long count;
		
		Reads.readLengths = new int[N_READS];
		Math.set(Reads.readLengths,N_READS-1,READ_LENGTH);
		
		
		System.err.println("Building the alphabet...");
		buildAlphabet(INPUT_FILE);
		System.err.println("Alphabet built. Size: "+(lastAlphabet+1));
		System.err.println("Histogram of lengths of recoded sequences:");
		Arrays.sort(sequenceLengths,0,lastSequenceLength+1);
		previous=sequenceLengths[0]; count=1;
		for (i=1; i<=lastSequenceLength; i++) {
			if (sequenceLengths[i]==previous) count++;
			else {
				System.err.println(previous+","+count);
				previous=sequenceLengths[i]; count=1;
			}
		}
		System.err.println(previous+","+count);
		
		System.err.println("Building k-mers...");
		buildKmers(INPUT_FILE);
		System.err.println(idGenerator+" k-mers built");
		kmerCountHistogram = new int[1001];
		for (i=0; i<kmerCountHistogram.length; i++) kmerCountHistogram[i]=0;
		edgeCountHistogram = new int[1001];
		for (i=0; i<edgeCountHistogram.length; i++) edgeCountHistogram[i]=0;
		degreeHistogram = new int[1001];
		for (i=0; i<degreeHistogram.length; i++) degreeHistogram[i]=0;
		kmerSet.kmerFrequencyStats();
		System.err.println("Histogram of node counts:");
		for (i=0; i<kmerCountHistogram.length; i++) System.err.println(i+","+kmerCountHistogram[i]);
		System.err.println("Histogram of edge counts:");
		for (i=0; i<edgeCountHistogram.length; i++) System.err.println(i+","+edgeCountHistogram[i]);
		System.err.println("Histogram of node degrees:");
		for (i=0; i<degreeHistogram.length; i++) System.err.println(i+","+degreeHistogram[i]);
		
		System.err.println("Printing connected components of the DBG...");
		stack = new Kmer[STACK_CAPACITY];
		idGenerator=0;
		kmerSet.kmerCleanConnectedComponents();
		kmerSet.kmerConnectedComponents();
		System.err.println((idGenerator+1)+" connected components built");
		

System.err.println("digraph G {");
kmerSet.kmerCompressUnaryPaths();
kmerSet.kmerToDot();
System.err.println("}");

	}
	
	
	
	
	
	
	
	
	
	
	
	// ----------------------------- ALIGNMENT PROCEDURES --------------------------------
	
	/**
	 * [UP TO DATE] Simple filters on repeat mappings.
	 */
	private static final int clean(AlignmentRow[] alignments, int lastAlignment) {
		int i, j;
		int length, maxLength, maxAlignment, intervalStart, intervalEnd;
		int startA, endA;
		String intervalReadB;
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
				if (alignments[j].startA<startA-DISTANCE_THRESHOLD) break;
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
		
		// If multiple substrings of the same repeat overlap in the read, and if the
		// repeat is not periodic, we keep only a longest alignment. We keep all
		// alignments if the repeat is periodic, so that they can be merged downstream.
		AlignmentRow.order=AlignmentRow.ORDER_READB_ORIENTATION_STARTA_ENDA_STARTB_ENDB;
		Arrays.sort(alignments,0,lastAlignment+1);
		for (i=0; i<=lastAlignment; i++) alignments[i].flag=false;
		intervalReadB=alignments[0].readB;
		intervalStart=alignments[0].startA;
		intervalEnd=alignments[0].endA;
		maxAlignment=0; maxLength=intervalEnd-intervalStart+1;
		for (i=1; i<=lastAlignment; i++) {
			if (alignments[i].readB!=intervalReadB || alignments[i].startA>intervalEnd-DISTANCE_THRESHOLD) {
				alignments[maxAlignment].flag=true;
				intervalReadB=alignments[i].readB;
				intervalStart=alignments[i].startA;
				intervalEnd=alignments[i].endA;
				maxAlignment=i; maxLength=intervalEnd-intervalStart+1;
				continue;
			}
			if (alignments[i].endA>intervalEnd) intervalEnd=alignments[i].endA;
			length=alignments[i].endA-alignments[i].startA+1;
			if (length>maxLength) {
				maxLength=length;
				maxAlignment=i;
			}
		}
		alignments[maxAlignment].flag=true;
		j=-1;
		for (i=0; i<=lastAlignment; i++) {
			if (!isPeriodic[alignments[i].readB] && !alignments[i].flag) continue;
			j++;
			tmpAlignment=alignments[j];
			alignments[j]=alignments[i];
			alignments[i]=tmpAlignment;
		}
		lastAlignment=j;
		
		return lastAlignment;
	}
	
	
	/**
	 * [UP TO DATE] A row of the alignments file.
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
			this.lengthA=-1; this.lengthB=-1;
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

	
	
	
	// ----------------------------- ALPHABET PROCEDURES ---------------------------------
	
	/**
	 * [UP TO DATE] Recodes every read and collects in $alphabet$ all distinct repeat 
	 * intervals.
	 *
	 * @param inputFile alignments between reads (readA) and repeats (readB);
	 * @param maxError alignments with error rate greater than this are discarded.
	 */
	private static final void buildAlphabet(String inputFile, double maxError) throws IOException {
		final int SEQUENCE_LENGTH_CAPACITY = 1000000;  // Arbitrary
		boolean found;
		int i, j, k;
		int lastAlignment, lastInSequence, row, readA, previousReadA;
		String str;
		BufferedReader br;
		Character tmpCharacter;
		
		// Allocating memory
		if (alphabet==null || alphabet.length<ALPHABET_CAPACITY) alphabet = new Character[ALPHABET_CAPACITY];
		lastAlphabet=-1;
		if (tmpArray==null || tmpArray.length<(ALIGNMENTS_CAPACITY)<<1) tmpArray = new int[(ALIGNMENTS_CAPACITY)<<1];
		if (sequence==null || sequence.length<SEQUENCE_CAPACITY) sequence = new Character[SEQUENCE_CAPACITY];
		if (alignments==null || alignments.length<ALIGNMENTS_CAPACITY) alignments = new AlignmentRow[ALIGNMENTS_CAPACITY];
		for (i=0; i<alignments.length; i++) {	
			if (alignments[i]==null) alignments[i] = new AlignmentRow();
		}
		
		// Collecting all characters (not necessarily distinct), from all recoded
		// sequences.
		sequenceLengths = new int[SEQUENCE_LENGTH_CAPACITY];
		lastSequenceLength=-1;
		br = new BufferedReader(new FileReader(inputFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); previousReadA=null; lastAlignment=-1; row=0;
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
					lastAlignment=clean(alignments,lastAlignment);
					lastInSequence=recode(alignments,lastAlignment,sequence,tmpArray);
					addCharacters2Alphabet(sequence,lastInSequence);
					addSequenceLength(lastInSequence+1);
				}
				previousReadA=readA; lastAlignment=0;
				alignments[0].set(readA,Alignments.startA,Alignments.endA,Alignments.readB-1,Alignments.startB,Alignments.endB,Alignments.orientation,Alignments.diffs);
			}
			else {
				lastAlignment++;
				if (lastAlignment==alignments.length) {
					AlignmentRow[] newAlignments = new AlignmentRow[alignments.length<<1];
					System.arraycopy(alignments,0,newAlignments,0,alignments.length);
					for (i=alignments.length; i<newAlignments.length; i++) newAlignments[i] = new AlignmentRow();
					alignments=newAlignments;
				}
				alignments[lastAlignment].set(readA,Alignments.startA,Alignments.endA,Alignments.readB-1,Alignments.startB,Alignments.endB,Alignments.orientation,Alignments.diffs);
			}
			str=br.readLine(); row++;
		}
		br.close();
		if (previousReadA!=null) {
			lastAlignment=clean(alignments,lastAlignment);
			lastInSequence=recode(alignments,lastAlignment,sequence,tmpArray);
			addCharacters2Alphabet(sequence,lastInSequence);
			addSequenceLength(lastInSequence+1);
		}
		
		// Sorting and compacting the alphabet
		Arrays.sort(alphabet,0,lastAlphabet+1);
		k=0;
		for (i=1; i<=lastAlphabet; i++) {
			found=false;
			for (j=k; j>=0; j--) {
				if (!alphabet[i].sameReadB(alphabet[j])) break;
				if (alphabet[i].isSimilar(alphabet[j])) {
					found=true;
					break;
				}
			}
			if (!found) {
				k++;
				tmpCharacter=alphabet[k];
				alphabet[k]=alphabet[i];
				alphabet[i]=tmpCharacter;
			}
		}
		lastAlphabet=k;
		
		// Assigning IDs
		for (i=0; i<=lastAlphabet; i++) alphabet[i].id=i;
	}
	
	
	private static final void addSequenceLength(int length) {
		lastSequenceLength++;
		if (lastSequenceLength==sequenceLengths.length) {
			int[] newSequenceLengths = new int[sequenceLengths.length<<1];
			System.arraycopy(sequenceLengths,0,newSequenceLengths,0,sequenceLengths.length);
			sequenceLengths=newSequenceLengths;
		}
		sequenceLengths[lastSequenceLength]=length;
	}
	
	
	/**
	 * Stores in $sequence$ a recoding of the read based on repeat mappings. The endpoints
	 * of all alignments with a non-periodic repeat are collected, as well as the
	 * endpoints of maximal ranges of overlapping alignments with periodic repeats. Points
	 * are clustered, and the resulting chunks of the read become characters of 
	 * $sequence$.
	 *
	 * Remark: the characters in $sequence$ are not elements of $alphabet$, and the 
	 * procedure does not look them up in $alphabet$: this is left to the caller.
	 *
	 * @param alignments all the alignments of a given readA;
	 * @param sequence assumed to be long enough;
	 * @param periodicIntervals temporary space, with at least $2(lastAlignment+1)$ cells;
	 * @param points temporary space, with at least $2(lastAlignment+1)$ cells;
	 * @return the last element in $sequence$.
	 */
	private static final int recode(AlignmentRow[] alignments, int lastAlignment, Character[] sequence, int[] periodicIntervals, int[] points) {
		int i, j;
		int lastPeriodicInterval, currentStart, currentEnd;
		int firstJForNextI, inPeriodic, lastPoint, previousPoint, numerator, denominator;
		
		int startA, endA, newLength;
		
		final int lengthA = Reads.getReadLength(alignments[0].readA);
		Character character;
		Character newCharacter = new Character();
		
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
			if (alignments[i].startA>=currentEnd-DISTANCE_THRESHOLD) {
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
		
		// Collecting readA endpoints of alignments
		i=0; j=0; firstJForNextI=-1; inPeriodic=-1; lastPoint=-1;
		while (i<=lastAlignment) {
			if (j>lastPeriodicInterval || periodicIntervals[j]>=alignments[i].endA-DISTANCE_THRESHOLD) {
				if (inPeriodic!=-1) {
					points[++lastPoint]=alignments[i].startA;
					points[++lastPoint]=alignments[i].endA;
				}
				else {
					if (Math.abs(alignments[i].startA,periodicIntervals[inPeriodic])<=DISTANCE_THRESHOLD) points[++lastPoint]=alignments[i].startA;
					if (Math.abs(alignments[i].endA,periodicIntervals[inPeriodic+1])<=DISTANCE_THRESHOLD) points[++lastPoint]=alignments[i].endA;
				}
				i++; inPeriodic=-1;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				continue;
			}
			if (periodicIntervals[j+1]<alignments[i].startA-DISTANCE_THRESHOLD) {
				j+=2;
				continue;
			}
			if (firstJForNextI==-1 && i<lastAlignment && periodicIntervals[j+1]>=alignments[i+1].startA) firstJForNextI=j;
			if ( Intervals.isApproximatelyContained(alignments[i].startA,alignments[i].endA,periodicIntervals[j],periodicIntervals[j+1]) ||
				 Intervals.areApproximatelyIdentical(alignments[i].startA,alignments[i].endA,periodicIntervals[j],periodicIntervals[j+1])
			   ) inPeriodic=j;
		}
		if (inPeriodic!=-1) {
			points[++lastPoint]=alignments[i].startA;
			points[++lastPoint]=alignments[i].endA;
		}
		else {
			if (Math.abs(alignments[i].startA,periodicIntervals[inPeriodic])<=DISTANCE_THRESHOLD) points[++lastPoint]=alignments[i].startA;
			if (Math.abs(alignments[i].endA,periodicIntervals[inPeriodic+1])<=DISTANCE_THRESHOLD) points[++lastPoint]=alignments[i].endA;
		}
		Arrays.sort(points,0,lastPoint+1);
		j=-1; previousPoint=points[0]; numerator=points[0]; denominator=1;
		for (i=1; i<=lastPoint; i++) {
			if (points[i]<=previousPoint+DISTANCE_THRESHOLD) {
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
		out=-1;
		// First non-repetitive character
		i=points[0];
		if (i>DISTANCE_THRESHOLD) {
			sequence[0].setNonrepetitive(i);
			sequence[0].openStart=true;
			sequence[0].openEnd=i>=lengthA-DISTANCE_THRESHOLD;
		}
		// Middle characters
		i=1; j=0; firstJForNextI=-1; out=-1;
		while (i<=lastPoint) {
			startA=points[i-1]; endA=points[i];
			if (j>lastAlignment || alignments[j].startA>=endA) {
				if (newCharacter.lastTuple==-1) newCharacter.setNonrepetitive(endA-startA+1);
				else Arrays.sort(newCharacter.tuples,0,newCharacter.lastTuple+1);
				out++;
				sequence[out].copyFrom(newCharacter);
				sequence[out].openStart=startA<=DISTANCE_THRESHOLD;
				sequence[out].openEnd=endA>=lengthA-DISTANCE_THRESHOLD;
				i++;
				if (firstJForNextI!=-1) j=firstJForNextI;
				firstJForNextI=-1;
				newCharacter.clear();
				continue;
			}
			else if (alignments[j].endA<=startA) {
				j++;
				continue;
			}
			if (firstJForNextI==-1 && i<lastPoint && alignments[j].endA>=endA) firstJForNextI=j;
			if ( Intervals.areApproximatelyIdentical(alignments[j].startA,alignments[j].endA,startA,endA) || 
				 Intervals.isApproximatelyContained(alignments[j].startA,alignments[j].endA,startA,endA)
			   ) newCharacter.add(alignments[j],startA,endA);
			j++;
		}
		// Last non-repetitive character
		i=points[lastPoint];
		if (lengthA-i>DISTANCE_THRESHOLD) {
			out++;
			sequence[out].setNonrepetitive(lengthA-i);
			sequence[out].openStart=i<=DISTANCE_THRESHOLD;
			sequence[out].openEnd=true;
		}
------>
		return out;
	}
	
	
	/**
	 * Appends all the characters in $sequence$ to the end of $alphabet$.
	 */
	private static final void addCharacters2Alphabet(Character[] sequence, int lastCharacter) {
		if (lastCharacter<2) return;
		
		final int newLength = lastAlphabet+1+lastCharacter+1-2;
		if (newLength>alphabet.length) {
			Character[] newAlphabet = new Character[(alphabet.length<<1)>newLength?(alphabet.length<<1):newLength<<1];
			System.arraycopy(alphabet,0,newAlphabet,0,alphabet.length);
			alphabet=newAlphabet;
		}
		System.arraycopy(sequence,1,alphabet,lastAlphabet+1,lastCharacter+1-2);
		lastAlphabet+=lastCharacter+1-2;
	}
	
	
	/**
	 * [ALMOST UP TO DATE] A character of the recoded alphabet, which is a set of substrings of 
	 * repeats in specific orientations, and with open/closed endpoints in readA.
	 */
	private static class Character implements Comparable {
		private static final int CAPACITY = 10;  // Arbitrary
		
		public int id;
		public Tuple[] tuples;
		public int lastTuple;
		public boolean openStart;  // TRUE: readA begins close to the start of the char
		public boolean openEnd;  // TRUE: readA ends close to the end of the char
		
		/**
		 * An empty character
		 */
		public Character() {
			tuples = new Tuple[CAPACITY];
			clear();
		}
		
		public void clear() {
			id=-1; lastTuple=-1; openStart=false; openEnd=false;
		}
		
		/**
		 * A non-repetitive character of given length.
		 */
		public Character(int length, boolean openStart, boolean openEnd) {
			id=-1;
			tuples = new Tuple[CAPACITY];
			setNonrepetitive(length);
			this.openStart=openStart; this.openEnd=openEnd;
		}
		
		/**
		 * Sets the current character to a non-repetitive of given length.
		 */
		public final void setNonrepetitive(int length) {
			tuples[0] = new Tuple(NON_REPETITIVE,true,-1,-1,length);
			lastTuple=0;
		}
		
		public final boolean isNonrepetitive() {
			return lastTuple==0 && tuples[0].readB.equals(NON_REPETITIVE);
		}
		
		private final void enlarge() {
			Tuple[] newTuples = new Tuple[tuples.length<<1];
			System.arraycopy(tuples,0,newTuples,0,tuples.length);
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
			openStart=otherCharacter.openStart;
			openEnd=otherCharacter.openEnd;
		}
		
		/**
		 * Sorting just by $readB,orientation$ of the tuples.
		 */
		public int compareTo(Object other) {
			int i, j;
			int last;
			
			Character otherCharacter = (Character)other;
			if (lastTuple<otherCharacter.lastTuple) return -1;
			else if (lastTuple>otherCharacter.lastTuple) return 1;
			for (i=0; i<=lastTuple; i++) {
				j=tuples[i].readB.compareTo(otherCharacter.tuples[i].readB);
				if (j!=0) return j;
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
		 * when considering only the $readB$ and $orientation$ fields of all tuples.
		 */
		public final boolean sameReadB(Character otherCharacter) {
			if (lastTuple!=otherCharacter.lastTuple) return false;
			for (int i=0; i<=lastTuple; i++) {
				if (!tuples[i].readB.equalsIgnoreCase(otherCharacter.tuples[i].readB) || tuples[i].orientation!=otherCharacter.tuples[i].orientation) return false;
			}
			return true;
		}
		
		/**
		 * Remark: the procedure assumes that the two characters are similar 
		 * according to $sameReadB()$.
		 */
		public final boolean isSimilar(Character otherCharacter) {
			for (int i=0; i<=lastTuple; i++) {
				if ( Math.abs(tuples[i].startB-otherCharacter.tuples[i].startB)>DISTANCE_THRESHOLD || 
					 Math.abs(tuples[i].endB-otherCharacter.tuples[i].endB)>DISTANCE_THRESHOLD || 
					 Math.abs(tuples[i].lengthB-otherCharacter.tuples[i].lengthB)>LENGTH_THRESHOLD
				   ) return false;
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
				out+=Math.pow(tuples[i].startB-otherCharacter.tuples[i].startB,2);
				out+=Math.pow(tuples[i].endB-otherCharacter.tuples[i].endB,2);
				out+=Math.pow(tuples[i].lengthB-otherCharacter.tuples[i].lengthB,2);
			}
			return Math.sqrt(out);
		}
		
		public int hashCode() {
			String out = id+",";
			for (int i=0; i<=lastTuple; i++) out+=tuples[i].toString()+",";
			return out.hashCode();
		}
		
		public String toString() {
			if (lastTuple==0 && tuples[0].readB==NON_REPETITIVE) return id+": "+NON_REPETITIVE+","+tuples[0].lengthB;
			String out = id+": ";
			for (int i=0; i<=lastTuple; i++) out+=tuples[i].readB+(tuples[i].orientation?"+":"-")+"["+tuples[i].startB+".."+tuples[i].endB+"], ";
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
		
		
-------> up to date until here		
		/**
		 * Appends the alignment of a non-short-period readB to the current 
		 * character, which is assumed to correspond to interval $[startA..
		 * endA]$ in readA.
		 *
		 * Remark: tuples might not be sorted after the procedure completes.
		 */
		public final void add(AlignmentRow alignment, int startA, int endA) {
			int startB, endB;
			double ratio;
			
			if (isIdentical(startA,endA,alignment.startA,alignment.endA)) {
				lastTuple++;
				if (lastTuple==tuples.length) enlarge();	
				tuples[lastTuple] = new Tuple(alignment.readB,alignment.orientation,alignment.startB,alignment.endB,0);
			}
			else if (isContained(startA,endA,alignment.startA,alignment.endA)) {
				lastTuple++;
				if (lastTuple==tuples.length) enlarge();
				ratio=alignment.getRatio();
				if (alignment.orientation) {
					startB=(int)(alignment.startB+(startA-alignment.startA)/ratio);
					endB=(int)(alignment.endB-(alignment.endA-endA)/ratio);
				}
				else {
					startB=(int)(alignment.startB+(alignment.endA-endA)/ratio);
					endB=(int)(alignment.endB-(startA-alignment.startA)/ratio);
				}
				tuples[lastTuple] = new Tuple(alignment.readB,alignment.orientation,startB,endB,0);
			}
			else {
				// NOP
			}
		}
	}
	
	
	/**
	 * [UP TO DATE] A substring of a readB in a specific orientation.
	 */
	private static class Tuple implements Comparable {
		public int readB;
		public boolean orientation;
		public int startB;
		public int endB;
		public int lengthB;  // >0: non-periodic or short-period; 0: others.
		
		public Tuple() {
			readB=-1; orientation=false; startB=-1; endB=-1; lengthB=-1;
		}
		
		public Tuple(String rb, boolean o, int sb, int eb, int lb) {
			this.readB=rb;
			this.orientation=o;
			this.startB=sb;
			this.endB=eb;
			this.lengthB=lb;
		}
		
		public int compareTo(Object other) {
			Tuple otherTuple = (Tuple)other;
			int i = readB.compareTo(otherTuple.readB);
			if (i!=0) return i;
			if (orientation && !otherTuple.orientation) return -1;
			else if (!orientation && otherTuple.orientation) return 1;
			if (startB<otherTuple.startB) return -1;
			else if (startB>otherTuple.startB) return 1;
			if (endB<otherTuple.endB) return -1;
			else if (endB>otherTuple.endB) return 1;
			if (lengthB<otherTuple.lengthB) return -1;
			else if (lengthB>otherTuple.lengthB) return 1;
			return 0;
		}
		
		public void copyFrom(Tuple otherTuple) {
			readB=otherTuple.readB;
			orientation=otherTuple.orientation;
			startB=otherTuple.startB;
			endB=otherTuple.endB;
			lengthB=otherTuple.lengthB;
		}
		
		public boolean equals(Object other) {
			Tuple otherTuple = (Tuple)other;
			return readB.equalsIgnoreCase(otherTuple.readB) && orientation==otherTuple.orientation && startB==otherTuple.startB && endB==otherTuple.endB && lengthB==otherTuple.lengthB;
		}
		
		public String toString() {
			return readB+","+orientation+","+startB+","+endB+","+lengthB;
		}
		
		public int getLength() {
			return (startB==-1||endB==-1)?lengthB:endB-startB+1;
		}
	}

	
	
		
	// ------------------------- K-MER PROCEDURES ------------------------------
	
	/**
	 * A simple trie that collects all k-mers.
	 */
	private static class KmerTree {
		public static final int CAPACITY = 2;  // Arbitrary
		
		public int character;  // Position in $alphabet$.
		public KmerTree parent;
		public KmerTree[] children;
		public int lastChild;
		public Kmer kmer;
		
		
		public KmerTree(int c) {
			this.character=c;
			parent=null;
			children = new KmerTree[CAPACITY];
			lastChild=-1;
			kmer=null;
		}
		
		
		/**
		 * @return the new descendant created, or the existing descendant of the
		 * node, with label equal to string $sequence[from..to]$.
		 */
		public KmerTree add(Character[] sequence, int from, int to) {
			int i;
			KmerTree currentNode;
			
			currentNode=this;
			for (i=from; i<=to; i++) currentNode=currentNode.addChild(sequence[i].id);
			return currentNode;
		}
		
		
		/**
		 * @return the new child created, or the existing child of the node with
		 * label $c$.
		 */
		public KmerTree addChild(int c) {
			for (int i=0; i<=lastChild; i++) {
				if (children[i].character==c) return children[i];
			}
			lastChild++;
			if (lastChild==children.length) {
				KmerTree[] newChildren = new KmerTree[children.length<<1];
				System.arraycopy(children,0,newChildren,0,children.length);
				children=newChildren;
			}
			children[lastChild] = new KmerTree(c);
			children[lastChild].parent=this;
			return children[lastChild];
		}
		
		
		public String getLabel(boolean simple) {
			KmerTree currentNode;
			String out;
			
			out=""; currentNode=this;
			while (currentNode!=kmerSet) {
				if (simple) out=alphabet[currentNode.character].id+"_"+out;
				else out=alphabet[currentNode.character].toString()+"\n"+out;
				currentNode=currentNode.parent;
			}
			return out;
		}
		
		
		public long getLengthBps() {
			KmerTree currentNode;
			long out;
			
			out=0; currentNode=this;
			while (currentNode!=kmerSet) {
				out+=alphabet[currentNode.character].getLength();
				currentNode=currentNode.parent;
			}
			return out;
		}
		
		
		public long getLengthBps_first() {
			KmerTree currentNode;
			
			currentNode=this;
			while (currentNode.parent!=kmerSet) currentNode=currentNode.parent;
			return alphabet[currentNode.character].getLength();
		}
		
		
		public long getLengthBps_last() {
			return alphabet[this.character].getLength();
		}
		
		
		public final void kmerFrequencyStats() {
			if (lastChild==-1) {
				kmer.frequencyStats();
				return;
			}
			for (int i=0; i<=lastChild; i++) children[i].kmerFrequencyStats();
		}
		
		
		public final void kmerCleanConnectedComponents() {
			if (lastChild==-1) {
				kmer.cleanConnectedComponent();
				return;
			}
			for (int i=0; i<=lastChild; i++) children[i].kmerCleanConnectedComponents();
		}
		
		
		public final void kmerConnectedComponents() {
			if (lastChild==-1) {
				kmer.getConnectedComponent();
				return;
			}
			for (int i=0; i<=lastChild; i++) children[i].kmerConnectedComponents();
		}
		
		
		public final void kmerToDot() {
			if (lastChild==-1) {
				kmer.toDot();
				return;
			}
			for (int i=0; i<=lastChild; i++) children[i].kmerToDot();
		}
		
		
		public final void kmerCompressUnaryPaths() {
			if (lastChild==-1) {
				kmer.compressUnaryPaths();
				return;
			}
			for (int i=0; i<=lastChild; i++) children[i].kmerCompressUnaryPaths();
		}
		
	}
	
	
	/**
	 * A node of the de Bruijn graph.
	 */
	private static class Kmer {
		public static final int CAPACITY = 2;  // Arbitrary
		private static final int NODE_SIZE = 10;
		private static final double WEIGHT_RATE = 0.01;
		
		public long id;
		public KmerTree node;
		public long count;
		public Kmer[] leftExtensions, rightExtensions;
		public long[] leftExtensionCounts, rightExtensionCounts;
		public int lastLeftExtension, lastRightExtension;
		public int connectedComponent;
		
		
		public Kmer() {
			id=-1; count=0; lastLeftExtension=-1; lastRightExtension=-1;
			leftExtensions = new Kmer[CAPACITY];
			rightExtensions = new Kmer[CAPACITY];
			rightExtensionCounts = new long[CAPACITY];
			leftExtensionCounts = new long[CAPACITY];
			connectedComponent=-1;
			node=null;
		}
		
		
		public void addRightExtension(Kmer neighbor) {
			for (int i=0; i<=lastRightExtension; i++) {
				if (rightExtensions[i].equals(neighbor)) {
					rightExtensionCounts[i]++;
					return;
				}
			}
			lastRightExtension++;
			if (lastRightExtension==rightExtensions.length) {
				int newLength = rightExtensions.length<<1;
				Kmer[] newArray = new Kmer[newLength];
				System.arraycopy(rightExtensions,0,newArray,0,rightExtensions.length);
				rightExtensions=newArray;
				long[] newArray2 = new long[newLength];
				System.arraycopy(rightExtensionCounts,0,newArray2,0,rightExtensionCounts.length);
				rightExtensionCounts=newArray2;
			}
			rightExtensions[lastRightExtension]=neighbor;
			rightExtensionCounts[lastRightExtension]=1;
		}
		
		
		public void addLeftExtension(Kmer neighbor) {
			for (int i=0; i<=lastLeftExtension; i++) {
				if (leftExtensions[i].equals(neighbor)) {
					leftExtensionCounts[i]++;
					return;
				}
			}
			lastLeftExtension++;
			if (lastLeftExtension==leftExtensions.length) {
				int newLength = leftExtensions.length<<1;
				Kmer[] newArray = new Kmer[newLength];
				System.arraycopy(leftExtensions,0,newArray,0,leftExtensions.length);
				leftExtensions=newArray;
				long[] newArray2 = new long[newLength];
				System.arraycopy(leftExtensionCounts,0,newArray2,0,leftExtensionCounts.length);
				leftExtensionCounts=newArray2;
			}
			leftExtensions[lastLeftExtension]=neighbor;
			leftExtensionCounts[lastLeftExtension]=1;
		}
		
		
		public void frequencyStats() {
			final int degree = lastLeftExtension+lastRightExtension+2;
			int i;
			
			kmerCountHistogram[(int)(count>=kmerCountHistogram.length?kmerCountHistogram.length-1:count)]++;
			for (i=0; i<=lastLeftExtension; i++) edgeCountHistogram[(int)(leftExtensionCounts[i]>=edgeCountHistogram.length?kmerCountHistogram.length-1:leftExtensionCounts[i])]++;
			for (i=0; i<=lastRightExtension; i++) edgeCountHistogram[(int)(rightExtensionCounts[i]>=edgeCountHistogram.length?edgeCountHistogram.length-1:rightExtensionCounts[i])]++;
			degreeHistogram[degree>=degreeHistogram.length?degreeHistogram.length-1:degree]++;
			
			
if (degree==29) {
	System.err.println("k-mer with degree "+degree+": ");
	System.err.println(node.getLabel(false));
	System.err.println("right neighbors:  "+(lastRightExtension+1));
	for (int x=0; x<=lastRightExtension; x++) {
		System.err.println(rightExtensions[x].node.getLabel(false));
		System.err.println();
	}
	System.err.println("left neighbors:  "+(lastLeftExtension+1));
	for (int x=0; x<=lastLeftExtension; x++) {
		System.err.println(leftExtensions[x].node.getLabel(false));
		System.err.println();
	}
}
			
		}
		
		
		public void cleanConnectedComponent() {
			connectedComponent=-1;
		}
		
		
		/**
		 * Sets the $connectedComponent$ field of every k-mer that is reachable 
		 * from this k-mer in the de Bruijn graph, and that does not already 
		 * have a connected component.
		 */
		public void getConnectedComponent() {
			if (connectedComponent>=0) return;
			final int newComponent = (int)(++idGenerator);
			int i;
			int top;
			Kmer currentKmer, neighbor;			
			
			this.connectedComponent=newComponent;
			stack[0]=this; top=0;
			while (top>=0) {
				currentKmer=stack[top--];
				ensureStack(top+1+(currentKmer.lastLeftExtension+1)+(currentKmer.lastRightExtension+1));
				for (i=0; i<=currentKmer.lastLeftExtension; i++) {
					neighbor=currentKmer.leftExtensions[i];
					if (neighbor.connectedComponent==-1) {
						neighbor.connectedComponent=newComponent;
						stack[++top]=neighbor;
					}
					else if (neighbor.connectedComponent!=newComponent) {
						System.err.println("getConnectedComponent> ERROR: different connected components "+newComponent+" :: "+neighbor.connectedComponent);
						System.exit(1);
					}
				}
				for (i=0; i<=currentKmer.lastRightExtension; i++) {
					neighbor=currentKmer.rightExtensions[i];
					if (neighbor.connectedComponent==-1) {
						neighbor.connectedComponent=newComponent;
						stack[++top]=neighbor;
					}
					else if (neighbor.connectedComponent!=newComponent) {
						System.err.println("getConnectedComponent> ERROR: different connected components "+newComponent+" :: "+neighbor.connectedComponent);
						System.exit(1);
					}
				}
			}
		}
		
		/**
		 *
		 */
		public void toDot() {
			final int degree = lastLeftExtension+lastRightExtension+2;
			if (degree==0) return;
			
			System.out.println(node.getLabel(true)+" [size=\""+((int)Math.log(NODE_SIZE*count))+"\"];");
			for (int i=0; i<=lastRightExtension; i++) System.out.println(node.getLabel(true)+" -> "+rightExtensions[i].node.getLabel(true)+"[weight=\""+Math.max(rightExtensionCounts[i]*WEIGHT_RATE,1.0)+"\"];");
			for (int i=0; i<=lastLeftExtension; i++) System.out.println(leftExtensions[i].node.getLabel(true)+" -> "+node.getLabel(true)+"[weight=\""+Math.max(leftExtensionCounts[i]*WEIGHT_RATE,1.0)+"\"];");
		}
		
		
		/**
		 * For simplicity it does not compress simple cycles.
		 */
		public void compressUnaryPaths() {
			boolean inSimpleCycle;
			int i, p, q;
			long pathLength;
			Kmer currentKmer, previousKmer, leftKmer, rightKmer;
			
			if (lastLeftExtension!=0 || lastRightExtension!=0 || leftExtensions[0]==rightExtensions[0]) return;
			rightKmer=rightExtensions[0];
			pathLength=node.getLengthBps();
			
//System.err.println("compressUnaryPaths> CALLED FROM "+this+" ----------------");			
			
			
			// Check for presence of a simple cycle
			currentKmer=this; inSimpleCycle=false;
			while (currentKmer.lastLeftExtension==0 && currentKmer.lastRightExtension==0) {
				currentKmer=currentKmer.leftExtensions[0];
				if (currentKmer==this) {
					inSimpleCycle=true;
					break;
				}
			}
			if (inSimpleCycle) return;
			
			// Left traversal
			currentKmer=this; previousKmer=null; leftKmer=leftExtensions[0];
			while (currentKmer.lastLeftExtension==0 && currentKmer.lastRightExtension==0) {
				
/*				
System.err.println("compressUnaryPaths> currentKmer: "+currentKmer);
System.err.println("compressUnaryPaths> previousKmer: "+previousKmer);
System.err.println("compressUnaryPaths> right-extensions: ");
for (int x=0; x<=currentKmer.lastRightExtension; x++) System.err.println(x+": "+currentKmer.rightExtensions[x]);
System.err.println("compressUnaryPaths> left-extensions: ");
for (int x=0; x<=currentKmer.lastLeftExtension; x++) System.err.println(x+": "+currentKmer.leftExtensions[x]);
*/
				
				currentKmer.lastLeftExtension=-1;
				currentKmer.lastRightExtension=-1;				
				previousKmer=currentKmer;
				currentKmer=leftKmer;
				pathLength+=currentKmer.node.getLengthBps_first();
				leftKmer=currentKmer.leftExtensions[0];
			}
			p=-1;
			for (i=0; i<=currentKmer.lastRightExtension; i++) {
				if (currentKmer.rightExtensions[i]==previousKmer) {
					p=i;
					break;
				}
			}
			leftKmer=currentKmer;
			pathLength-=leftKmer.node.getLengthBps_first();
			
			// Right traversal
			currentKmer=rightKmer; previousKmer=this; rightKmer=currentKmer.rightExtensions[0];
			pathLength+=currentKmer.node.getLengthBps_last();
			while (currentKmer.lastLeftExtension==0 && currentKmer.lastRightExtension==0) {
				currentKmer.lastLeftExtension=-1;
				currentKmer.lastRightExtension=-1;				
				previousKmer=currentKmer;
				currentKmer=rightKmer;
				pathLength+=currentKmer.node.getLengthBps_last();
				rightKmer=currentKmer.rightExtensions[0];
			}
			q=-1;
			for (i=0; i<=currentKmer.lastLeftExtension; i++) {
				if (currentKmer.leftExtensions[i]==previousKmer) {
					q=i;
					break;
				}
			}
			rightKmer=currentKmer;
			pathLength-=rightKmer.node.getLengthBps_last();
			
			// Connecting the two ends of the unary path
			leftKmer.rightExtensions[p]=this;
			this.lastLeftExtension=0;
			this.leftExtensions[0]=leftKmer;
			this.lastRightExtension=0;
			this.rightExtensions[0]=rightKmer;
			rightKmer.leftExtensions[q]=this;
			
			System.err.println("unary path length bps: "+pathLength);
		}
	}
	
	
	/**
	 * Recodes every read and builds the set of distinct k-mers.
	 */
	private static final void buildKmers(String inputFile) throws IOException {
		final int SEQUENCE_LENGTH_CAPACITY = 1000000;  // Arbitrary
		boolean found, orientation;
		int i, j;
		int lastAlignment, startA, endA, startB, endB, diffs, lengthA, lengthB;
		int mask, lastInSequence, row, previous, count;
		String str, readA, previousReadA, readB;
		BufferedReader br;
		Character tmpCharacter;
		String[] tokens;
				
		// Collecting all characters (not necessarily distinct), from all
		// recoded sequences.
		kmerSet = new KmerTree(0); idGenerator=0;
		br = new BufferedReader(new FileReader(inputFile));
		str=br.readLine(); previousReadA=null; lastAlignment=-1; row=0;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments, "+idGenerator+" k-mers collected.");
			tokens=str.split(PAF_SEPARATOR);
			readA=tokens[0]; 
			readB=tokens[5];
			orientation=tokens[4].charAt(0)=='+';
			startA=Integer.parseInt(tokens[2]);
			endA=Integer.parseInt(tokens[3]);
			startB=Integer.parseInt(tokens[7]);
			endB=Integer.parseInt(tokens[8]);
			diffs=Integer.parseInt(tokens[10])-Integer.parseInt(tokens[9]);
			lengthA=Integer.parseInt(tokens[1]);
			lengthB=Integer.parseInt(tokens[6]);
			if (previousReadA==null || !readA.equalsIgnoreCase(previousReadA)) {
				if (previousReadA!=null) {
					lastAlignment=clean(alignments,lastAlignment);
					lastInSequence=recode(alignments,lastAlignment,sequence,tmpArray);
					
					
Character[] backup = new Character[lastInSequence+1];
System.arraycopy(sequence,0,backup,0,lastInSequence+1);
boolean fabio = toAlphabet(sequence,lastInSequence);
if (fabio) {
	System.err.println("FUCK! FOUND THE TRIMER IN SEQUENCE "+previousReadA);
	for (int x=0; x<=lastInSequence; x++) System.err.println("---> "+x+": "+backup[x]+"  ---> "+sequence[x]);
	System.err.println();
}



					extractKmers(sequence,lastInSequence);
				}
				previousReadA=readA; lastAlignment=0;
				alignments[0].set(readA,startA,endA,readB,startB,endB,orientation,diffs,lengthA,lengthB);
			}
			else {
				lastAlignment++;
				if (lastAlignment==alignments.length) {
					AlignmentRow[] newAlignments = new AlignmentRow[alignments.length<<1];
					System.arraycopy(alignments,0,newAlignments,0,alignments.length);
					for (i=alignments.length; i<newAlignments.length; i++) newAlignments[i] = new AlignmentRow();
					alignments=newAlignments;
				}
				alignments[lastAlignment].set(readA,startA,endA,readB,startB,endB,orientation,diffs,lengthA,lengthB);
			}
			str=br.readLine(); row++;
		}
		br.close();
		if (previousReadA!=null) {
			lastAlignment=clean(alignments,lastAlignment);
			lastInSequence=recode(alignments,lastAlignment,sequence,tmpArray);
boolean fabio = toAlphabet(sequence,lastInSequence);
if (fabio) System.err.println("FUCK! FOUND THE TRIMER IN SEQUENCE "+previousReadA);
			extractKmers(sequence,lastInSequence);
		}
	}
	
	
	/**
	 * Assumes that $alphabet$ has already been built, and replaces every 
	 * character in $sequence$ with the result of its lookup in $alphabet$.
	 */
	private static final boolean toAlphabet(Character[] sequence, int lastCharacter) {
		int i;
		Character newCharacter;

		if (lastCharacter<2) return false;
		for (i=1; i<lastCharacter; i++) {
			newCharacter=inAlphabet(sequence[i]);
			if (newCharacter==null) {
				System.err.println("toAlphabet> ERROR: character not found in the alphabet: "+sequence[i]);
				System.exit(1);
			}
			sequence[i]=newCharacter;
		}
		
		
		
		
		
		
		
		
boolean found=false;		
for (int x=0; x<lastCharacter-4; x++) {
/*	if ( sequence[x].lastTuple==0 && sequence[x].tuples[0].readB.equalsIgnoreCase("HERVH") && sequence[x].tuples[0].orientation && sequence[x].tuples[0].startB==25 && sequence[x].tuples[0].endB==3965 &&
		 sequence[x+1].isNonrepetitive() && sequence[x+1].tuples[0].lengthB==1278 &&
		 sequence[x+2].lastTuple==0 && sequence[x+2].tuples[0].readB.equalsIgnoreCase("LTR7A") && sequence[x+2].tuples[0].orientation && sequence[x+2].tuples[0].startB==3 && sequence[x+2].tuples[0].endB==430
	   ) {
*/
	if (sequence[x].id==44 && sequence[x+1].id==48793 && sequence[x+2].id==105 && sequence[x+3].id==48593 && sequence[x+4].id==39) {	
		   System.err.println("FUCK! I FOUND THE 5-mer AT POSITION "+x+" IN A SEQUENCE");
		   found=true;
	   }
}		
return found;
		
		
	}
	
	
	/**
	 * Remark: the procedure assumes $alphabet$ to be sorted.
	 *
	 * @return an element of $alphabet$ that is most similar to $query$, or NULL
	 * if no such element exists.
	 */
	private static final Character inAlphabet(Character query) {
		int i, p;
		double distance, minDistance;
		Character character, minCharacter;
		
		p=Arrays.binarySearch(alphabet,0,lastAlphabet+1,query);
		if (p<0) p=-1-p;
		minCharacter=null; minDistance=Integer.MAX_VALUE;
		for (i=p-1; i>=0; i--) {
			character=alphabet[i];
			if (!character.sameReadB(query)) break;
			if (character.isSimilar(query)) {
				distance=character.getDistance(query);
				if (distance<minDistance) {
					minDistance=distance;
					minCharacter=character;
				}	
			}
		}
		for (i=p; i<=lastAlphabet; i++) {
			character=alphabet[i];
			if (!character.sameReadB(query)) break;
			if (character.isSimilar(query)) {
				distance=character.getDistance(query);
				if (distance<minDistance) {
					minDistance=distance;
					minCharacter=character;
				}	
			}
		}
		return minCharacter;
	}
	
	
	/**
	 * Adds to $kmerSet$ all substrings of length $K_VALUE$ of $sequence$.
	 */
	private static final void extractKmers(Character[] sequence, int lastCharacter) {
		int i;
		Character newCharacter;
		KmerTree node;
		Kmer currentKmer, previousKmer;

		if (lastCharacter+1-2<K_VALUE) return;
		previousKmer=null;
		for (i=1; i<=lastCharacter-K_VALUE; i++) {
			node=kmerSet.add(sequence,i,i+K_VALUE-1);
			if (node.kmer==null) {
				node.kmer = new Kmer();
				node.kmer.id=++idGenerator;
				node.kmer.node=node;
			}
			currentKmer=node.kmer;
			currentKmer.count++;
			if (previousKmer!=null) {
				previousKmer.addRightExtension(currentKmer);
				currentKmer.addLeftExtension(previousKmer);
			}
			previousKmer=currentKmer;
		}
	}
	
	
	private static final void ensureStack(int newSize) {
		if (stack.length>=newSize) return;
		Kmer[] newStack = new Kmer[(stack.length<<1)<newSize?newSize:stack.length<<1];
		System.arraycopy(stack,0,newStack,0,stack.length);
		stack=newStack;
	}
	
	
	
	
	// -------------------------- DRAWING PROCEDURES ---------------------------
	
	/**
	 * Draws a cleaned alignment of every read.
	 */
	private static final void draw(String inputFile, boolean transparent) throws IOException {
		boolean orientation;
		int i;
		int lastAlignment, startA, endA, startB, endB, diffs, lengthA, lengthB;
		int row;
		String str, readA, previousReadA, readB;
		BufferedReader br;
		String[] tokens;
		
		// Allocating memory
		if (alignments==null || alignments.length<ALIGNMENTS_CAPACITY) alignments = new AlignmentRow[ALIGNMENTS_CAPACITY];
		for (i=0; i<alignments.length; i++) alignments[i] = new AlignmentRow();
		
		// Drawing
		br = new BufferedReader(new FileReader(inputFile));
		str=br.readLine(); previousReadA=null; lastAlignment=-1; row=0;
		while (str!=null)  {
			if (row%100000==0) System.err.println("Processed "+row+" alignments, alphabet="+(lastAlphabet+1));
			tokens=str.split(PAF_SEPARATOR);
			readA=tokens[0]; 
			readB=tokens[5];
			orientation=tokens[4].charAt(0)=='+';
			startA=Integer.parseInt(tokens[2]);
			endA=Integer.parseInt(tokens[3]);
			startB=Integer.parseInt(tokens[7]);
			endB=Integer.parseInt(tokens[8]);
			diffs=Integer.parseInt(tokens[10])-Integer.parseInt(tokens[9]);
			lengthA=Integer.parseInt(tokens[1]);
			lengthB=Integer.parseInt(tokens[6]);
			if (previousReadA==null || !readA.equalsIgnoreCase(previousReadA)) {
				if (previousReadA!=null) {
					lastAlignment=clean(alignments,lastAlignment);
					draw_impl(alignments,lastAlignment,transparent);
				}
				previousReadA=readA; lastAlignment=0;
				alignments[0].set(readA,startA,endA,readB,startB,endB,orientation,diffs,lengthA,lengthB);
			}
			else {
				lastAlignment++;
				if (lastAlignment==alignments.length) {
					AlignmentRow[] newAlignments = new AlignmentRow[alignments.length<<1];
					System.arraycopy(alignments,0,newAlignments,0,alignments.length);
					for (i=alignments.length; i<newAlignments.length; i++) newAlignments[i] = new AlignmentRow();
					alignments=newAlignments;
				}
				alignments[lastAlignment].set(readA,startA,endA,readB,startB,endB,orientation,diffs,lengthA,lengthB);
			}
			str=br.readLine(); row++;
		}
		br.close();
		if (previousReadA!=null) {
			lastAlignment=clean(alignments,lastAlignment);
			draw_impl(alignments,lastAlignment,transparent);
		}
	}
	
	
	private static final void draw_impl(AlignmentRow[] alignments, int lastAlignment, boolean transparent) throws IOException {
		final int QUANTUM = 1;
		final int COLOR_BACKGROUND = 0x00050A0D;
		final int N_COLUMNS = alignments[0].lengthA/QUANTUM;
		final int N_ROWS = HEADER_ROWS+lastAlignment+1;
		final int HALF_MAX_HEIGHT = N_ROWS>>1;
		final int TRUNCATION_TAG_WIDTH_PIXELS = 3;
		final int TRUNCATION_MULTIPLE= 1000;
		final int TRANSPARENT_BORDER_PIXELS = 2;
		final int COLOR_REFERENCE = 0x00244254;
		final int COLOR_REFERENCE_HIGHLIGHT = 0x00509FB5;
		final String DEFAULT_FONT = "Barlow";
		int i, j, x, y;
		int firstX, lastX, leftTruncation, rightTruncation, aLength, height, fromY;
		String label;
		BufferedImage image;
		final NumberFormat formatter = NumberFormat.getInstance();
		formatter.setMaximumFractionDigits(2);
		
		image = new BufferedImage(N_COLUMNS,N_ROWS,BufferedImage.TYPE_INT_RGB);
		for (x=0; x<N_COLUMNS; x++) {
			for (y=0; y<N_ROWS; y++) setRGB(image,x,y,COLOR_BACKGROUND);
		}
		fromY=0;
		for (j=0; j<=lastAlignment; j++) {
			// Drawing truncated triangles
			leftTruncation=alignments[j].orientation?alignments[j].startB:alignments[j].lengthB-alignments[j].endB;
			rightTruncation=alignments[j].orientation?alignments[j].lengthB-alignments[j].endB:alignments[j].startB;
			firstX=alignments[j].startA-leftTruncation;
			lastX=alignments[j].endA+rightTruncation;
			aLength=lastX-firstX+1;
			for (x=alignments[j].startA; x<=alignments[j].endA; x++) {
				height=Math.round(HALF_MAX_HEIGHT*(alignments[j].orientation?lastX-x:x-firstX)/aLength);
				if ( (leftTruncation>DISTANCE_THRESHOLD && ((x-alignments[j].startA)/QUANTUM)<=TRUNCATION_TAG_WIDTH_PIXELS) || 
					 (rightTruncation>DISTANCE_THRESHOLD && ((alignments[j].endA-x)/QUANTUM)<=TRUNCATION_TAG_WIDTH_PIXELS)
				   ) {
					for (y=fromY+HALF_MAX_HEIGHT-height; y<=fromY+HALF_MAX_HEIGHT+height; y++) setRGB(image,x/QUANTUM,y,COLOR_REFERENCE_HIGHLIGHT);
				}
				else {
					if (transparent) {
						if (x<=alignments[j].startA+TRANSPARENT_BORDER_PIXELS || x>=alignments[j].endA-TRANSPARENT_BORDER_PIXELS) {
							for (y=fromY+HALF_MAX_HEIGHT-height; y<=fromY+HALF_MAX_HEIGHT+height; y++) setRGB(image,x/QUANTUM,y,COLOR_REFERENCE);
						} 
						else {
							for (y=fromY+HALF_MAX_HEIGHT-height; y<=fromY+HALF_MAX_HEIGHT-height+TRANSPARENT_BORDER_PIXELS; y++) setRGB(image,x/QUANTUM,y,COLOR_REFERENCE);
							for (y=fromY+HALF_MAX_HEIGHT+height-TRANSPARENT_BORDER_PIXELS; y<=fromY+HALF_MAX_HEIGHT+height; y++) setRGB(image,x/QUANTUM,y,COLOR_REFERENCE);
						}
					}
					else {
						for (y=fromY+HALF_MAX_HEIGHT-height; y<=fromY+HALF_MAX_HEIGHT+height; y++) setRGB(image,x/QUANTUM,y,COLOR_REFERENCE);
					}
				}
			}

			// Drawing text
			final Graphics2D g2d = image.createGraphics();
			g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
			g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
			g2d.setColor(new Color(COLOR_REFERENCE_HIGHLIGHT));
			final FontMetrics fontMetrics = g2d.getFontMetrics();
			label=""+alignments[j].readB.toUpperCase();
			x=((alignments[j].startA+alignments[j].endA)/2)/QUANTUM-(int)Math.round(((double)fontMetrics.stringWidth(label))/2);
			y=fromY+HALF_MAX_HEIGHT+fontMetrics.getAscent()/2;
			g2d.drawString(label,x,y);
			if (leftTruncation>DISTANCE_THRESHOLD) {
				label=""+formatter.format(((double)leftTruncation)/TRUNCATION_MULTIPLE);
				x=alignments[j].startA/QUANTUM+TRUNCATION_TAG_WIDTH_PIXELS+TRUNCATION_TAG_WIDTH_PIXELS;
				y=fromY+HALF_MAX_HEIGHT+fontMetrics.getAscent()/2;
				g2d.drawString(label,x,y);
			}
			if (rightTruncation>DISTANCE_THRESHOLD) {
				label=""+formatter.format(((double)rightTruncation)/TRUNCATION_MULTIPLE);
				x=alignments[j].endA/QUANTUM-TRUNCATION_TAG_WIDTH_PIXELS-TRUNCATION_TAG_WIDTH_PIXELS-fontMetrics.stringWidth(label);
				y=fromY+HALF_MAX_HEIGHT+fontMetrics.getAscent()/2;
				g2d.drawString(label,x,y);
			}
		}
		
		// Writing to disk
		ImageIO.write(image,"png",new File(alignments[0].readA.replace('/','-')+".png"));
	}
	
	
	private static final void setRGB(BufferedImage image, int x, int y, int color) {
		if (x<0 || y<0 || x>=image.getWidth() || y>=image.getHeight()) return;
		image.setRGB(x,y,color);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 * Adds $query to $alphabet$ and sorts $alphabet$.
	 */
/*	private static final void addToAlphabet(Character query) {
		lastAlphabet++;
		if (lastAlphabet==alphabet.length) {
			Character[] newAlphabet = new Character[alphabet.length<<1];
			System.arraycopy(alphabet,0,newAlphabet,0,alphabet.length);
			alphabet=newAlphabet;
		}
		alphabet[lastAlphabet]=query;
		Arrays.sort(alphabet,0,lastAlphabet+1);
	}
*/	
	
	
	

}