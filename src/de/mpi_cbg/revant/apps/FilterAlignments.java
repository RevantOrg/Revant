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
	 * Marks periodic repeats: 0=noperiodic; 1=short-period; 2=long-period.
	 */
	private static final byte[] isPeriodic;
	
	
	
	private static final int DISTANCE_THRESHOLD = 20;  // Arbitrary
	private static final int LENGTH_THRESHOLD = 100;  // Arbitrary
	private static final String NON_REPETITIVE = "nonrepetitive";
	private static final String PAF_SEPARATOR = "\t";
	private static final int HEADER_ROWS = 100;  // in pixels
	private static int rFrom, gFrom, bFrom, rTo, gTo, bTo;
	
	
	
	
	
	private static int K_VALUE;
	private static long idGenerator;
	
	/**
	 * Every distinct k-mer in a recoded read
	 */
	private static KmerTree kmerSet;
	
	private static int[] kmerCountHistogram, edgeCountHistogram, degreeHistogram;
	private static Kmer[] stack;
	
	
	
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
	 * [UP TO DATE] Basic filters on repeat alignments.
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
			if (isPeriodic[alignments[i].readB]!=0) continue;
			startA=alignments[i].startA; endA=alignments[i].endA;
			for (j=i+1; j<=lastAlignment; j++) {
				if (alignments[j].startA>=endA-distanceThreshold) break;
				if (isPeriodic[alignments[j].readB]!=0) continue;
				if (Intervals.areApproximatelyIdentical(alignments[j].startA,alignments[j].endA,startA,endA) continue;
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
	 * The recoded alphabet
	 */
	private static Character[] alphabet;
	private static int lastAlphabet;
	
	/**
	 * All alignments of a given readA
	 */
	private static AlignmentRow[] alignments;
	private static int lastAlignment;
	
	/**
	 * The recoded sequence of a given readA
	 */
	private static Character[] sequence;
	private static int lastInSequence;
	
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
	 * Recodes every read and collects in $alphabet$ all distinct characters.
	 *
	 * @param inputFile alignments between reads (readA) and repeats (readB);
	 * @param maxError alignments with error rate greater than this are discarded.
	 */
	private static final void buildAlphabet(String inputFile, double maxError, int distanceThreshold) throws IOException {
		final int ALPHABET_CAPACITY = 100000;  // Arbitrary
		final int ALIGNMENTS_CAPACITY = 100000;  // Arbitrary
		final int SEQUENCE_CAPACITY = 1000000;  // Arbitrary
		final int SEQUENCE_LENGTH_CAPACITY = 1000000;  // Arbitrary
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
		if (sequenceLengths==null || sequenceLengths.length<SEQUENCE_LENGTH_CAPACITY) sequenceLengths = new int[SEQUENCE_LENGTH_CAPACITY];
		if (newCharacter==null) newCharacter = new Character();
		if (tmpTuple==null) tmpTuple = new Tuple();
		
		// Collecting all characters (not necessarily distinct) from all recoded reads.
		br = new BufferedReader(new FileReader(inputFile));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine(); 
		previousReadA=null; lastAlignment=-1; row=0; lastSequenceLength=-1;
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
					cleanAlignments();
					recodeRead(distanceThreshold);
					if (lastInSequence>0 || (lastInSequence==0 && !sequence[0].isNonrepetitive())) addCharacters2Alphabet();
					addSequenceLength(lastInSequence+1);
				}
				previousReadA=readA; lastAlignment=0;
				alignments[0].set(readA,Math.max(Alignments.startA,0),Math.min(Alignments.endA,Reads.getReadLength(readA)-1),Alignments.readB-1,Math.max(Alignments.startB,0),Math.min(Alignments.endB,Reads.getReadLength(readB)-1),Alignments.orientation,Alignments.diffs);
			}
			else {
				lastAlignment++;
				if (lastAlignment==alignments.length) {
					AlignmentRow[] newAlignments = new AlignmentRow[alignments.length<<1];
					System.arraycopy(alignments,0,newAlignments,0,alignments.length);
					for (i=alignments.length; i<newAlignments.length; i++) newAlignments[i] = new AlignmentRow();
					alignments=newAlignments;
				}
				alignments[lastAlignment].set(readA,Math.max(Alignments.startA,0),Math.min(Alignments.endA,Reads.getReadLength(readA)-1),Alignments.readB-1,Math.max(Alignments.startB,0),Math.min(Alignments.endB,Reads.getReadLength(readB)-1),Alignments.orientation,Alignments.diffs);
			}
			str=br.readLine(); row++;
		}
		br.close();
		if (previousReadA!=null) {
			cleanAlignments();
			recodeRead(distanceThrehsold);
			if (lastInSequence>0 || (lastInSequence==0 && !sequence[0].isNonrepetitive())) addCharacters2Alphabet();
			addSequenceLength(lastInSequence+1);
		}
		
		// Sorting and compacting the alphabet
		--------------->
		
		// Assigning IDs
		for (i=0; i<=lastAlphabet; i++) alphabet[i].id=i;
	}
	
	
	private static final void addSequenceLength(int length) {
		final int GROWTH_RATE = 100000;  // Arbitrary
		
		lastSequenceLength++;
		if (lastSequenceLength==sequenceLengths.length) {
			int[] newSequenceLengths = new int[sequenceLengths.length+GROWTH_RATE];
			System.arraycopy(sequenceLengths,0,newSequenceLengths,0,sequenceLengths.length);
			sequenceLengths=newSequenceLengths;
		}
		sequenceLengths[lastSequenceLength]=length;
	}
	
	
	/**
	 * [UP TO DATE] Stores in global variable $sequence$ a recoding of the read based on repeat 
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
	 * Remark: the characters of $sequence$ are not elements of $alphabet$, and the 
	 * procedure does not look them up in $alphabet$: this is left to the caller.
	 *
	 * Remark: the procedure uses global arrays $alignments,periodicIntervals,points,
	 * sequence$.
	 */
	private static final void recodeRead(int distanceThreshold) {
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
			if (isPeriodic[alignments[i].readB]==0) continue;
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
				else if (isPeriodic[alignments[i].readB]!=0) {
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
		else if (isPeriodic[alignments[i].readB]!=0) {
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
		i=1; j=0; firstJForNextI=-1; out=-1;
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
	 * [UP TO DATE] Appends all the characters of $sequence$ to the end of $alphabet$.
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
	 * Temporary space encoding an undirected graph and its connected 
	 * components.
	 */
	private static int[][] neighbors;
	private static int[] lastNeighbor;
	private static int[] stack, connectedComponent, connectedComponentSize;
	private static int nComponents;
	private static Char[] mergedChars;
	
	
	
	/**
	 *
	 *
	 */
	private static final void compactAlphabet(int distanceThreshold, int lengthThreshold) {
		int i, j, k;
		int isClosedI, isClosedJ;
		Character tmpChar;
		
		// Discarding characters implied by other characters
		Arrays.sort(alphabet,0,lastAlphabet+1);
		for (i=0; i<=lastAlphabet; i++) alphabet[i].id=1;
		for (i=0; i<=lastAlphabet; i++) {
			if (alphabet[i].isNonrepetitive()) continue;
			isClosedI=alphabet[i].isClosed();b
			for (j=i+1; j<=lastAlphabet; j++) {
				if (!alphabet[j].sameRepeat(alphabet[i])) break;
				if (alphabet[j].isNonrepetitive()) continue;
				isClosedJ=alphabet[j].isClosed();
				if (isClosedJ==2 || (isClosedJ==1 && isClosedI==0)) continue;
				if (alphabet[i].implies(alphabet[j])) alphabet[j].id=0;
			}
			for (j=i-1; j>=0; j--) {
				if (!alphabet[j].sameRepeat(alphabet[i])) break;
				if (alphabet[j].isNonrepetitive()) continue;
				isClosedJ=alphabet[j].isClosed();
				if (isClosedJ==2 || (isClosedJ==1 && isClosedI==0)) continue;
				if (alphabet[i].implies(alphabet[j])) alphabet[j].id=0;
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
		
		
		
		
		
		// Paritioning closed and open characters
		j=-1;
		for (i=0; i<=lastAlphabet; i++) {
			if (alphabet[i].isClosed()==2) {
				j++;
				tmpChar=alphabet[j];
				alphabet[j]=alphabet[i];
				alphabet[i]=tmpChar;
			}
		}
		lastClosed=j;
		
		// Merging characters with closed ends
		ensureGraph(lastClosed+1);
		Arrays.sort(alphabet,0,lastClosed+1);
		k=0;
		for (i=1; i<lastClosed; i++) {
			for (j=i+1; j<=lastClosed; j++) {
				if (!alphabet[j].sameRepeat(alphabet[i])) break;
				if (alphabet[j].isSimilar(alphabet[i],distanceThreshold,lengthThreshold)) {
					addEdge(i,j); addEdge(j,i);
				}
			}
		}
		nComponents=getConnectedComponent(lastClosed+1);
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
		for (i=0; i<=lastClosed; i++) {
			tmpChar=mergedChars[connectedComponent[i]];
			if (tmpChar.lastTuple==-1) tmpChar.copyFrom(alphabet[i]);
			else tmpChar.addTupleEnds(alphabet[i]);
		}
		for (i=0; i<nComponents; i++) mergedChars[i].normalizeTupleEnds(connectedComponentSize[i]);
		for (i=0; i<nComponents; i++) alphabet[i].copyFrom(mergedChars[i]);
		lastClosedPrime=nComponents-1;
		
		// 
		Arrays.sort(alphabet,lastClosed+1,lastAlphabet+1);
		
		
		
		
		// Merging remaining open characters
		
		
		
		
	}
	
	
	
	private static final void ensureGraph(int nNodes) {
		final int NEIGHBORS_CAPACITY = 10;  // Arbitrary
		
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
	 * [UP TO DATE] A character of the recoded alphabet, which is a set of substrings of 
	 * repeats in specific orientations and with open/closed endpoints.
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
		
		public void reset() {
			id=-1; lastTuple=-1;
		}
		
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
		public final boolean isSimilar(Character otherCharacter, int distanceThrehsold, int lengthThreshold) {
			for (int i=0; i<=lastTuple; i++) {
				if ( Math.abs(tuples[i].start-otherCharacter.tuples[i].start)>distanceThreshold || 
					 Math.abs(tuples[i].end-otherCharacter.tuples[i].end)>distanceThreshold || 
					 Math.abs(tuples[i].length-otherCharacter.tuples[i].length)>lengthThreshold
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
			if (lastTuple==0 && tuples[0].read==NON_REPETITIVE) return id+": "+NON_REPETITIVE+","+tuples[0].length;
			String out = id+": ";
			for (int i=0; i<=lastTuple; i++) out+=tuples[i].read+(tuples[i].orientation?"+":"-")+"["+tuples[i].start+".."+tuples[i].end+"], ";
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
				tmpTuple.openEnd=alignment.endA>=Reads.getReadLength(alignment.readB)-distanceThreshold;
			}
			else {
				tmpTuple.openEnd=alignment.startA<=distanceThreshold;
				tmpTuple.openStart=alignment.endA>=Reads.getReadLength(alignment.readB)-distanceThreshold;
			}
			tmpTuple.repeat=alignment.readB; tmpTuple.orientation=alignment.orientation;
			if (isPeriodic[alignment.readB]==0) {
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
				if (tuples[i].repeat==tmpTuple.read && tuples[i].orientation==tmpTuple.orientation) {
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
				if (isPeriodic[alignment.readB]!=0) tuples[found].length=Math.max(tuples[found].length,tmpTuple.length);
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
		public static final void cleanTuples(int periodicLength) {
			boolean foundPeriodic, foundNonperiodic, openLeft, openRight;
			int i, j;
			Tuple tmpTuple;
			
			// Enforcing periodic/nonperiodic criteria
			foundPeriodic=false; foundNonperiodic=false;
			for (i=0; i<=lastTuple; i++) {
				if (isPeriodic[tuples[i].repeat]!=0) {
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
					if (isPeriodic[tuples[i].repeat]==0) continue;
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
				if (tuples[i].orientation) { openStart=openLeft; openEnd=openRight; }
				else { openStart=openRight; openEnd=openLeft; }
			}			
			
			Arrays.sort(tuples,0,lastTuple+1);
		}
		
		
		/**
		 * @return 0=fully open; 1=half-open; 2=fully closed.
		 */
		public final int isClosed() {
			if (tuples[0].openStart && tuples[0].openEnd) return 0;
			else if (tuples[0].openStart || tuples[0].openEnd) return 1;
			else return 2;
		}
		
		
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
		 * Remark: the characters are assumed both to be repetitive and to satisfy
		 * $sameRepeat()$.
		 *
		 * @return TRUE iff no tuple of $otherCharacter$ adds information WRT the 
		 * corresponding tuple of this character.
		 */
		public static final boolean implies(Character otherCharacter, int distanceThreshold) {
			int i;

			for (i=0; i<=lastTuple; i++) {
				if (tuples[i].start==-1) {  // Periodic
					if (Math.abs(tuples[i].length,otherCharacter.tuples[i].length)<distanceThreshold) {
						if ( (tuples[i].openStart && !otherCharacter.tuples[i].openStart) ||
						     (tuples[i].openEnd && !otherCharacter.tuples[i].openEnd)
						   ) return false;
					}
					else if (tuples[i].length<otherCharacter.tuples[i].length) return false;
					else if (!otherCharacter.tuples[i].openStart && !otherCharacter.tuples[i].openEnd) return false;
				}
				else {  // Nonperiodic
					if (Intervals.areApproximatelyIdentical(otherCharacter.tuples[i].start,otherCharacter.tuples[i].end,tuples[i].start,tuples[i].end)) {
						if ( (tuples[i].openStart && !otherCharacter.tuples[i].openStart) ||
						     (tuples[i].openEnd && !otherCharacter.tuples[i].openEnd) ||
							 (tuples[i].openStart && otherCharacter.tuples[i].openStart && tuples[i].openEnd && otherCharacter.tuples[i].openEnd)
						   ) return false;
					}
					else if (Intervals.isApproximatelyContained(otherCharacter.tuples[i].start,otherCharacter.tuples[i].end,tuples[i].start,tuples[i].end)) {
						if ( (Math.abs(tuples[i].start,otherCharacter.tuples[i].start)<=distanceThreshold && tuples[i].openStart && !otherCharacter.tuples[i].openStart) ||
							 (otherCharacter.tuples[i].start>tuples[i].start+distanceThreshold && !otherCharacter.tuples[i].openStart) ||
							 (Math.abs(tuples[i].end,otherCharacter.tuples[i].end)<=distanceThreshold && tuples[i].openEnd && !otherCharacter.tuples[i].openEnd) ||
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
	 * [UP TO DATE] A substring of a repeat in a specific orientation.
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
			return repeat+","+(orientation?"FWD":"REV")+"["+start+".."+end+"] (L="+length+")";
		}
		
		public int getLength() {
			return (start==-1||end==-1)?length:end-start+1;
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