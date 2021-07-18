package de.mpi_cbg.revant.factorize;

import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Histograms;


public class Alignments {
	
	public static final String BPS_LABEL = "bps,";
	public static final int BPS_LABEL_LENGTH = BPS_LABEL.length();
	
	/**
	 * Parameters of the pipeline
	 */
	public static final double MAX_ALIGNMENT_ERROR = 0.3;
	public static final double MIN_ALIGNMENT_ERROR = 0.2;  // Alignments of identical substrings have error rate between $MIN_ALIGNMENT_ERROR$ and $MAX_ALIGNMENT_ERROR$ in a typical PacBio run.
	public static int maxAlignmentsPerRead;  // Upper bound on the number of alignments per readA
	public static int maxAlignmentsPerPair;  // Upper bound on the number of alignments per (readA,readB) pair.
	public static int maxReadsPerRead;  // Upper bound on the number of distinct readB per readA
	public static int maxOccurrences;  // Upper bound on the number of occurrences of a dense or periodic substring in the pile of a readA
	public static final boolean AGGRESSIVE_MAXIMALITY = false;  // Affects $isLeftMaximal$ and $isRightMaximal$
	public static int minAlignmentLength;  // The aligner returns alignments of at least this length
	public static int maximalityTolerance = IO.quantum<<1;

	/**
	 * Data structures
	 */
	public static int nAlignments;
	public static int[][] alignments;
	public static boolean[] marked;
	public static int[] tmp;  // $tmp[0]$: maximal; $tmp[1]$: B-maximal.
	
	
	/**
	 * Temporary space loaded by procedure $readAlignmentFile$.
	 */
	public static boolean orientation;
	public static int readA, readB, startA, endA, startB, endB, diffs;	
	private static StringBuilder sb = new StringBuilder();


	public static final void loadAlignments(String path) throws IOException {
		int i, p, q;
		int comma, previousReadA, previousReadB, readLength;
		int first, firstPrime, nAlignmentsPerRead, nAlignmentsPerPair, nReadB, nOccurrences;
		int mask;
		String str;
		BufferedReader br;

		tmp = new int[2];
		alignments = new int[nAlignments][8];
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		str=br.readLine(); str=br.readLine();  // Skipping the first two lines
		str=br.readLine();
		i=0; 
		previousReadA=-1; previousReadB=-1; 
		maxAlignmentsPerRead=0; first=0; 
		nReadB=0; maxReadsPerRead=0;
		maxOccurrences=0; nOccurrences=0;
		maxAlignmentsPerRead=0; nAlignmentsPerPair=0; firstPrime=0;
		while (str!=null) {
			readAlignmentFile(str);
			alignments[i][0]=readA;
			alignments[i][1]=readB;
			alignments[i][2]=orientation?1:0;
			alignments[i][3]=startA<0?0:startA;
			readLength=Reads.getReadLength(readA-1);
			alignments[i][4]=endA-1>=readLength?readLength-1:endA-1;  // The interval excludes position $endA$.
			alignments[i][5]=startB;  // $startB$ and $endB$ are already reversed by $readAlignmentFile$, if necessary.
			if (readB-1<Reads.readLengths.length) {
				readLength=Reads.getReadLength(readB-1);
				alignments[i][6]=endB-1>=readLength?readLength-1:endB-1;  // The interval excludes position $endB$.
			}
			else alignments[i][6]=endB-1;
			alignments[i][7]=diffs;

			// Updating stats
			if (readA!=previousReadA) {
				nAlignmentsPerRead=i-first;
				if (nAlignmentsPerRead>maxAlignmentsPerRead) maxAlignmentsPerRead=nAlignmentsPerRead;
				first=i; 
				if (nReadB>maxReadsPerRead) maxReadsPerRead=nReadB;
				nReadB=1;
				nOccurrences/=minAlignmentLength;
				if (nOccurrences>maxOccurrences) maxOccurrences=nOccurrences;
				nOccurrences=readB-1<Reads.readLengths.length?Reads.getReadLength(readB-1):Reads.maxReadLength;
				nAlignmentsPerPair=i-firstPrime;
				if (nAlignmentsPerPair>maxAlignmentsPerPair) maxAlignmentsPerPair=nAlignmentsPerPair;
				firstPrime=i;
			}
			else if (readB!=previousReadB) {
				nAlignmentsPerPair=i-firstPrime;
				if (nAlignmentsPerPair>maxAlignmentsPerPair) maxAlignmentsPerPair=nAlignmentsPerPair;
				firstPrime=i;
				nReadB++;
				nOccurrences+=readB-1<Reads.readLengths.length?Reads.getReadLength(readB-1):Reads.maxReadLength;
			}
			previousReadA=readA; previousReadB=readB;
			
			i++;
			str=br.readLine();
		}
		// Updating stats at the end of the file
		nAlignmentsPerRead=i-first;
		if (nAlignmentsPerRead>maxAlignmentsPerRead) maxAlignmentsPerRead=nAlignmentsPerRead;
		if (nReadB>maxReadsPerRead) maxReadsPerRead=nReadB;
		nOccurrences/=minAlignmentLength;
		if (nOccurrences>maxOccurrences) maxOccurrences=nOccurrences;
		nAlignmentsPerPair=i-firstPrime;
		if (nAlignmentsPerPair>maxAlignmentsPerPair) maxAlignmentsPerPair=nAlignmentsPerPair;
		maxOccurrences>>=1;
		
		br.close();
		marked = new boolean[nAlignments];
		//IO.printErr("maxAlignmentsPerRead="+maxAlignmentsPerRead+" maxAlignmentsPerPair="+maxAlignmentsPerPair);
	}
	
	
	/**
	 * Loads row $str$ of the alignment file into corresponding global variables
	 */
	public static final boolean readAlignmentFile(String str1) {
		if (str1==null || str1.length()==0) return false;
		int p, q;
		int tmpStartB, tmpEndB;
			
		p=0;
		while (p<str1.length() && str1.charAt(p)==' ') p++;
		q=p+1;
		while (q<str1.length() && str1.charAt(q)!=' ') q++;
		readA=intSubstring(str1,p,q);
		p=q+1;
		while (p<str1.length() && str1.charAt(p)==' ') p++;
		q=p+1;
		while (q<str1.length() && str1.charAt(q)!=' ') q++;
		readB=intSubstring(str1,p,q);
		p=q+1;
		while (p<str1.length() && str1.charAt(p)==' ') p++;
		orientation=str1.charAt(p)=='n';
		p=nextParenthesis(true,p,str1);
		q=str1.indexOf(".",p+1);
		startA=intSubstring(str1,p+1,q);
		if (startA<0) startA=0;
		p=q+1;
		while (p<str1.length() && str1.charAt(p)=='.') p++;
		q=nextParenthesis(false,p+1,str1);	
		endA=intSubstring(str1,p,q);
		p=nextParenthesis(true,q+1,str1);
		q=str1.indexOf(".",p+1);
		tmpStartB=intSubstring(str1,p+1,q);
		p=q+1;
		while (p<str1.length() && str1.charAt(p)=='.') p++;
		q=nextParenthesis(false,p+1,str1);
		tmpEndB=intSubstring(str1,p,q);
		startB=orientation?tmpStartB:tmpEndB;
		endB=orientation?tmpEndB:tmpStartB;
		p=str1.indexOf(BPS_LABEL,q+1);  // "New" format.
		if (p>=0) p+=BPS_LABEL_LENGTH-1;
		else p=str1.indexOf("(",q+1);  // "Old" format.
		q=str1.indexOf("d",p+1);
		diffs=intSubstring(str1,p+1,q);
		
		return true;
	}
	
	
	/**
	 * Variant of $readAlignmentFile()$ that just loads $readA$.
	 *
	 * @return first position in $str1$ after $readA$ (or -1 if $str1$ is emtpy).
	 */
	public static final int readAlignmentFile_readA(String str1) {
		if (str1==null || str1.length()==0) return -1;
		int p, q;
			
		p=0;
		while (p<str1.length() && str1.charAt(p)==' ') p++;
		q=p+1;
		while (q<str1.length() && str1.charAt(q)!=' ') q++;
		readA=intSubstring(str1,p,q);
		return q;
	}
	
	
	/**
	 * Loads the parts of $str1$ that were not loaded by $readAlignmentFile_readA()$.
	 *
	 * @param from output of $readAlignmentFile_readA()$.
	 */
	public static final boolean readAlignmentFile_rest(String str1, int from) {
		int p, q;
		int tmpStartB, tmpEndB;
		
		p=from+1;
		while (p<str1.length() && str1.charAt(p)==' ') p++;
		q=p+1;
		while (q<str1.length() && str1.charAt(q)!=' ') q++;
		readB=intSubstring(str1,p,q);
		p=q+1;
		while (p<str1.length() && str1.charAt(p)==' ') p++;
		orientation=str1.charAt(p)=='n';
		p=nextParenthesis(true,p,str1);
		q=str1.indexOf(".",p+1);
		startA=intSubstring(str1,p+1,q);
		if (startA<0) startA=0;
		p=q+1;
		while (p<str1.length() && str1.charAt(p)=='.') p++;
		q=nextParenthesis(false,p+1,str1);	
		endA=intSubstring(str1,p,q);
		p=nextParenthesis(true,q+1,str1);
		q=str1.indexOf(".",p+1);
		tmpStartB=intSubstring(str1,p+1,q);
		p=q+1;
		while (p<str1.length() && str1.charAt(p)=='.') p++;
		q=nextParenthesis(false,p+1,str1);
		tmpEndB=intSubstring(str1,p,q);
		startB=orientation?tmpStartB:tmpEndB;
		endB=orientation?tmpEndB:tmpStartB;
		p=str1.indexOf(BPS_LABEL,q+1);  // "New" format.
		if (p>=0) p+=BPS_LABEL_LENGTH-1;
		else p=str1.indexOf("(",q+1);  // "Old" format.
		q=str1.indexOf("d",p+1);
		diffs=intSubstring(str1,p+1,q);
		
		return true;
	}
	
	
	/**
	 * Stores in $array$ the global variables written by $readAlignmentsFile()$.
	 */
	public static final void readAlignmentFile_store(int[] array) {
		array[0]=orientation?0:1;
		array[1]=readA;
		array[2]=readB;
		array[3]=startA;
		array[4]=endA;
		array[5]=startB;
		array[6]=endB;
		array[7]=diffs;
	}
	
	
	/**
	 * Loads from $array$ the global variables written by $readAlignmentsFile()$.
	 */
	public static final void readAlignmentFile_load(int[] array) {
		orientation=array[0]==0;
		readA=array[1];
		readB=array[2];
		startA=array[3];
		endA=array[4];
		startB=array[5];
		endB=array[6];
		diffs=array[7];
	}
	
	
	/**
	 * @return TRUE iff $array$ encodes the same information as the global variables 
	 * written by $readAlignmentsFile()$.
	 */
	public static final boolean readAlignmentFile_equals(int[] array) {
		return orientation==(array[0]==0) && readA==array[1] && readB==array[2] && startA==array[3] && endA==array[4] && startB==array[5] && endB==array[6] && diffs==array[7];
	}
	
	
	/**
	 * Flips the roles of readA and readB in the global variables written by
	 * $readAlignmentsFile()$.
	 */
	public static final void readAlignmentFile_flip() {
		int tmp;
		
		tmp=readA; readA=readB; readB=tmp;
		tmp=startA; startA=startB; startB=tmp;
		tmp=endA; endA=endB; endB=tmp;
	}
	
	
	public static final int nextParenthesis(boolean open, int p, String str) {
		final int k = str.indexOf(open?"[":"]",p);
		final int kPrime = str.indexOf(open?"<":">",p);
		if (k>=0 && kPrime>=0) return k<kPrime?k:kPrime;
		if (k>=0) return k;
		if (kPrime>=0) return kPrime;
		return -1;
	}


	/**
	 * Converts a substring of the concatenation file into an integer.
	 * Indices work as in $String.substring$.
	 */
	public static final int intSubstring(String str, int start, int end) {
		char c;
		int i;
		
		sb.delete(0,sb.length());
		for (i=start; i<end; i++) {
			c=str.charAt(i);
			if (c!=',' && c!=' ') sb.append(c);
		}
		return Integer.parseInt(sb.toString());
	}


	/**
	 * Tests whether row $alignmentID$ of matrix $alignments$ is a left-maximal alignment,
	 * i.e. whether we can prove, using just $alignments$, that the corresponding
	 * substring of $readA$ occurs in $readB$ with a different left-context of high
	 * quality. The procedure assumes that the alignment algorithm returns maximal
	 * alignments.
	 *
	 * If $AGGRESSIVE_MAXIMALITY=true$, the procedure skips substrings of low quality
	 * immediately to the left of the alignment in $readA$ and $readB$, and it searches
	 * for an alignment that is close to the left boundary of the low quality substring to
	 * prove that the left context of the alignment is similar in the two reads.
	 * In practice returning 1 if we don't find an alignment with such features is not
	 * a reliable proof of maximality (such alignment might not have been returned by the
	 * aligner, or it might start too far from the boundary), and it is better to set
	 * $AGGRESSIVE_MAXIMALITY=false$.
	 *
	 * Remark: the procedure stores its output in array $tmp$: $tmp[0]=1$: the alignment
	 * is left-maximal; 0: the alignment is not left-maximal; -1: cannot be decided.
	 * $tmp[1]$ contains the output for B-maximality.
	 */
	public static final void isLeftMaximal(int alignmentID) {
		final int randomInsertionWindow = 1;
		boolean readAHasInsertion, readBHasInsertion;
		int i, a, b, aPrime, bPrime, distance, distancePrime;
		int readA = alignments[alignmentID][0]-1;
		int readB = alignments[alignmentID][1]-1;
		boolean orientation = alignments[alignmentID][2]==1;
		int startA = alignments[alignmentID][3];
		int startB = alignments[alignmentID][5];
		int endB = alignments[alignmentID][6];

		isLeftMaximal_impl(readA,startA,readB,startB);
		if (!AGGRESSIVE_MAXIMALITY) return;

		readAHasInsertion=Reads.hasLowQuality(readA,startA-randomInsertionWindow*Reads.QUALITY_SPACING,startA-1,true);
		readBHasInsertion=Reads.hasLowQuality(readB,startB-randomInsertionWindow*Reads.QUALITY_SPACING,startB-1,true);
		if ( (readAHasInsertion && (startA-randomInsertionWindow*Reads.QUALITY_SPACING)<minAlignmentLength) ||
			 (readBHasInsertion && (startB-randomInsertionWindow*Reads.QUALITY_SPACING)<minAlignmentLength) ) {
			tmp[0]=-1;
			return;
		}

		// Finding the rightmost $b$ such that
		// $readA[b-randomInsertionWindow*Reads.QUALITY_SPACING+1..b]$ has high quality,
		// and the leftmost $a$ such that $readA[a+1..b]$ has high quality, and such that
		// $b-a \geq minAlignmentLength$. Doing the same for $readB$.
		if (readAHasInsertion) {
			b=skipLeft(readA,true,startA-(randomInsertionWindow<<1)*Reads.QUALITY_SPACING,false);
			if (b==-1) {
				tmp[0]=-1;
				return;
			}
			b+=randomInsertionWindow*Reads.QUALITY_SPACING-1;
		}
		else b=startA-1;
		if (b<minAlignmentLength) {
			tmp[0]=-1;
			return;
		}
		a=skipLeft(readA,true,b-(randomInsertionWindow<<1)*Reads.QUALITY_SPACING,true);
		if (a==-1) a=0;
		else a+=randomInsertionWindow*Reads.QUALITY_SPACING;
		if (b-a<minAlignmentLength) {
			tmp[0]=-1;
			return;
		}
		if (readBHasInsertion) {
			bPrime=skipLeft(readB,true,startB-(randomInsertionWindow<<1)*Reads.QUALITY_SPACING,false);
			if (bPrime==-1) {
				tmp[0]=-1;
				return;
			}
			bPrime+=randomInsertionWindow*Reads.QUALITY_SPACING-1;
		}
		else bPrime=startB-1;
		if (bPrime<minAlignmentLength) {
			tmp[0]=-1;
			return;
		}
		aPrime=skipLeft(readB,true,bPrime-(randomInsertionWindow<<1)*Reads.QUALITY_SPACING,true);
		if (aPrime==-1) aPrime=0;
		else aPrime+=randomInsertionWindow*Reads.QUALITY_SPACING;
		if (bPrime-aPrime<minAlignmentLength) {
			tmp[0]=-1;
			return;
		}

		// Finding an alignment between $readA[a..b]$ and $readB[aPrime..bPrime]$, in the
		// same orientation.
		i=alignmentID-1;
		while ( i>=0 &&
		        alignments[i][0]==alignments[alignmentID][0] &&
		        alignments[i][1]==alignments[alignmentID][1] &&
		        alignments[i][2]==alignments[alignmentID][2] ) {
		    distance=b-alignments[i][4];
		    if (distance<0) distance=-distance;
		    distancePrime=bPrime-alignments[i][6];
		    if (distancePrime<0) distancePrime=-distancePrime;
		    if ( ( distance<=maximalityTolerance &&
		           Intervals.isApproximatelyContained(alignments[i][3],alignments[i][4],a,b)
		     	 ) &&
		     	 ( distancePrime<=maximalityTolerance &&
		     	   Intervals.isApproximatelyContained(alignments[i][5],alignments[i][6],aPrime,bPrime)
		    	 )
		       ) {
		    	tmp[0]=0;
		    	return;  // The left-context in $readA$ is similar to the left-context in $readB$
		    }
		    i--;
		}
		tmp[0]=1;
	}
	
	
	public static final void isLeftMaximal_impl(int readA, int startA, int readB, int startB) {
		if (Reads.isLeftMaximal(startB,readB,true)) tmp[1]=1;
		else tmp[1]=0;
		if (tmp[1]!=1 || !Reads.isLeftMaximal(startA,readA,true)) tmp[0]=0;
		else tmp[0]=1;
	}	


	/**
	 * Tests whether row $alignmentID$ of matrix $alignments$ is a right-maximal alignment
	 * i.e. whether we can prove, using just $alignments$, that the corresponding
	 * substring of $readA$ occurs in $readB$ with a different right-context of high
	 * quality. The procedure assumes that the alignment algorithm returns maximal
	 * alignments.
	 *
	 * If $AGGRESSIVE_MAXIMALITY=true$, the procedure skips substrings of low quality
	 * immediately to the right of the alignment in $readA$ and $readB$, and it searches
	 * for an alignment that is close to the right boundary of the low quality substring
	 * to prove that the right context of the alignment is similar in the two reads.
	 * In practice returning 1 if we don't find an alignment with such features is not
	 * a reliable proof of maximality (such alignment might not have been returned by the
	 * aligner, or it might start too far from the boundary), and it is better to set
	 * $AGGRESSIVE_MAXIMALITY=false$.
	 *
	 * Remark: the procedure stores its output in array $tmp$: $tmp[0]=1$: the alignment
	 * is right-maximal; 0: the alignment is not right-maximal; -1: cannot be decided.
	 * $tmp[1]$ contains the output for B-maximality.
	 */
	public static final void isRightMaximal(int alignmentID) {
		final int randomInsertionWindow = 1;
		boolean readAHasInsertion, readBHasInsertion;
		int i, a, b, aPrime, bPrime, distance, distancePrime;
		int readA = alignments[alignmentID][0]-1;
		int readB = alignments[alignmentID][1]-1;
		boolean orientation = alignments[alignmentID][2]==1;
		int endA = alignments[alignmentID][4];
		int startB = alignments[alignmentID][5];
		int endB = alignments[alignmentID][6];
		
		isRightMaximal_impl(readA,endA,readB,endB);
		if (!AGGRESSIVE_MAXIMALITY) return;

		readAHasInsertion=Reads.hasLowQuality(readA,endA+1,endA+randomInsertionWindow*Reads.QUALITY_SPACING,true);
		readBHasInsertion=Reads.hasLowQuality(readB,endB+1,endB+randomInsertionWindow*Reads.QUALITY_SPACING,true);		
		if ( (readAHasInsertion && Reads.getReadLength(readA)-(endA+randomInsertionWindow*Reads.QUALITY_SPACING)-1<minAlignmentLength) ||
			 (readBHasInsertion && Reads.getReadLength(readB)-(endB+randomInsertionWindow*Reads.QUALITY_SPACING)-1<minAlignmentLength) ) {
			tmp[0]=-1;
			return;
		}

		// Finding the leftmost $a$ such that $readA[endA+1..a-1]$ has low quality,
		// and the rightmost $b$ such that $readA[a+1..b]$ has high quality, and such that
		// $b-a \geq minAlignmentLength$. Doing the same for $readB$.
		if (readAHasInsertion) {
			a=skipRight(readA,true,endA+randomInsertionWindow*Reads.QUALITY_SPACING,false);
			if (a==-1) {
				tmp[0]=-1;
				return;
			}
		}
		else a=endA+1;
		if (Reads.getReadLength(readA)-a<minAlignmentLength) {
			tmp[0]=-1;
			return;
		}
		b=skipRight(readA,true,a+randomInsertionWindow*Reads.QUALITY_SPACING,true);
		if (b==-1) b=Reads.getReadLength(readA)-1;
		else b--;
		if (b-a+1<minAlignmentLength) {
			tmp[0]=-1;
			return;
		}
		if (readBHasInsertion) {
			aPrime=skipRight(readB,true,endB+randomInsertionWindow*Reads.QUALITY_SPACING,false);
			if (aPrime==-1) {
				tmp[0]=-1;
				return;
			}
		}
		else aPrime=endB+1;
		if (Reads.getReadLength(readB)-aPrime<minAlignmentLength) {
			tmp[0]=-1;
			return;
		}
		bPrime=skipRight(readB,true,aPrime+randomInsertionWindow*Reads.QUALITY_SPACING,true);
		if (bPrime==-1) bPrime=Reads.getReadLength(readB)-1;
		else bPrime--;
		if (bPrime-aPrime+1<minAlignmentLength) {
			tmp[0]=-1;
			return;
		}

		// Finding an alignment between $readA[a..b]$ and $readB[aPrime..bPrime]$, in the
		// same orientation.
		i=alignmentID+1;
		while ( i<nAlignments &&
		        alignments[i][0]==alignments[alignmentID][0] &&
		        alignments[i][1]==alignments[alignmentID][1] &&
		        alignments[i][2]==alignments[alignmentID][2] ) {
		    distance=a-alignments[i][3];
		    if (distance<0) distance=-distance;
		    distancePrime=aPrime-alignments[i][5];
		    if (distancePrime<0) distancePrime=-distancePrime;
		    if ( ( distance<=maximalityTolerance &&
		           Intervals.isApproximatelyContained(alignments[i][3],alignments[i][4],a,b)
		     	 ) &&
		     	 ( distancePrime<=maximalityTolerance &&
		     	   Intervals.isApproximatelyContained(alignments[i][5],alignments[i][6],aPrime,bPrime)
		    	 )
		       ) {
		    	tmp[0]=0;
		    	return;  // The left-context in $readA$ is similar to the left-context in $readB$
		    }
		    i++;
		}
		tmp[0]=1;
	}
	
	
	public static final void isRightMaximal_impl(int readA, int endA, int readB, int endB) {
		if (Reads.isRightMaximal(endB,readB,true)) tmp[1]=1;
		else tmp[1]=0;
		if (tmp[1]!=1 || !Reads.isRightMaximal(endA,readA,true)) tmp[0]=0;
		else tmp[0]=1;
	}
	
	
	/**
	 * Wraps $Histograms.skipLeft$ to work on $Reads.qualities$.
	 */
	private static final int skipLeft(int readID, boolean orientation, int start, boolean sign) {
		final int WINDOW = 2;
		
		int i = Histograms.skipLeft(Reads.getQualityArray(readID),Reads.getQualityArrayLength(readID),orientation,start/Reads.QUALITY_SPACING,WINDOW,Reads.MIN_RANDOM_QUALITY_SCORE,sign);
		if (i==-1) return -1;
		return i*Reads.QUALITY_SPACING;
	}
	
	
	/**
	 * Wraps $Histograms.skipRight$ to work on $Reads.qualities$.
	 */
	private static final int skipRight(int readID, boolean orientation, int start, boolean sign) {
		final int WINDOW = 2;
		
		int i = Histograms.skipRight(Reads.getQualityArray(readID),Reads.getQualityArrayLength(readID),orientation,start/Reads.QUALITY_SPACING+1,WINDOW,Reads.MIN_RANDOM_QUALITY_SCORE,sign);
		if (i==-1) return -1;
		return i*Reads.QUALITY_SPACING;
	}


	/**
	 * Let $[a..b]$, $[x1..x2]$ and $[y1..y2]$ be intervals in $readA$. The procedure
	 * returns the distance between the first position of $[y1..y2] \intersect [a..b]$ and
	 * the first position of $[x1..x2] \intersect [a..b]$, where intersections are assumed
	 * to be nonempty.
	 *
	 * @param alignment1 id of the alignment in $alignments$ with interval in $readA$
	 * equal to $[x1..x2]$;
	 * @param alignment2 id of the alignment in $alignments$ with interval in $readA$
	 * equal to $[y1..y2]$.
	 */
	public static final int getAShift(int alignment1, int alignment2, int a) {
		int f1 = alignments[alignment1][3]>a?alignments[alignment1][3]:a;
		int f2 = alignments[alignment2][3]>a?alignments[alignment2][3]:a;
		return f2-f1;
	}


	/**
	 * Let $[a..b]$ be an interval in $readA$, and let $([x1..x2],[x1',x2'])$ and
	 * $([y1..y2],[y1'..y2'])$ be alignments between $readA$ and $readB$. The procedure
	 * estimates the distance between the first position of $[y1'..y2']$ that aligns to
	 * $[a..b]$ and the first position of $[x1'..x2']$ that aligns to $[a..b]$.
	 * The procedure assumes that $[x1..x2] \intersect [a..b] \neq \emptyset$ and
	 * $[y1..y2] \intersect [a..b] \neq \emptyset$.
	 *
	 * @param alignment1 id of the alignment in $alignments$ with interval in $readA$
	 * equal to $[x1..x2]$;
	 * @param alignment2 id of the alignment in $alignments$ with interval in $readA$
	 * equal to $[y1..y2]$.
	 */
	public static final int getBShift(int alignment1, int alignment2, int a) {
		int f1, f2;
		double dilation;

		dilation=((double)(alignments[alignment1][6]-alignments[alignment1][5]))/(alignments[alignment1][4]-alignments[alignment1][3]);
		f1=alignments[alignment1][3]>a?alignments[alignment1][3]:a;
		f1=alignments[alignment1][5]+(int)((f1-alignments[alignment1][3])*dilation);
		dilation=((double)(alignments[alignment2][6]-alignments[alignment2][5]))/(alignments[alignment2][4]-alignments[alignment2][3]);
		f2=alignments[alignment2][3]>a?alignments[alignment2][3]:a;
		f2=alignments[alignment2][5]+(int)((f2-alignments[alignment2][3])*dilation);
		return f2-f1;
	}


	/**
	 * Let $[a..b]$ be an interval in $readA$, and let $([x1..x2],[x1',x2'])$ be an
	 * alignment between $readA$ and $readB$. The procedure estimates the first position
	 * $f$ in $[x1'..x2']$ that aligns to $[a..b]$. The procedure assumes that $[x1..x2]
	 * \intersect [a..b] \neq \emptyset$.
	 *
	 * @param alignment id of the alignment in $alignments$ with interval in $readA$
	 * equal to $[x1..x2]$.
	 */
	public static final int getF(int alignment, int a) {
		int f, readBLength;
		double dilation;

		dilation=((double)(alignments[alignment][6]-alignments[alignment][5]))/(alignments[alignment][4]-alignments[alignment][3]);
		f=alignments[alignment][3]>a?alignments[alignment][3]:a;
		f=alignments[alignment][5]+(int)((f-alignments[alignment][3])*dilation);
		readBLength=Reads.getReadLength(alignments[alignment][1]-1);
		if (f>=readBLength) f=readBLength-1;
		return f;
	}


	/**
	 * Let $[a..b]$ be an interval in $readA$, and let $([x1..x2],[x1',x2'])$ be an
	 * alignment between $readA$ and $readB$. The procedure estimates the last position
	 * $g$ in $[x1'..x2']$ that aligns to $[a..b]$. The procedure assumes that $[x1..x2]
	 * \intersect [a..b] \neq \emptyset$.
	 *
	 * @param alignment id of the alignment in $alignments$ with interval in $readA$
	 * equal to $[x1..x2]$.
	 */
	public static final int getG(int alignment, int b) {
		int g, readBLength;
		double dilation;

		dilation=((double)(alignments[alignment][6]-alignments[alignment][5]))/(alignments[alignment][4]-alignments[alignment][3]);
		g=alignments[alignment][4]<b?alignments[alignment][4]:b;
		g=alignments[alignment][6]-(int)((alignments[alignment][4]-g)*dilation);
		if (g<0) g=0;
		return g;
	}


	public static final int getALength(int id) {
		return alignments[id][4]-alignments[id][3]+1;
	}


	public static final int getBLength(int id) {
		return alignments[id][6]-alignments[id][5]+1;
	}


	public static final double getAvgDiffs(int id) {
		return ((double)(alignments[id][7]<<1))/(getALength(id)+getBLength(id));
	}
	
	
	/**
	 * Projects the number of differences in $[oldStartA..oldEndA] x [oldStartB..oldEndB]$
	 * onto $[newStartA..newEndA] x [newStartB..newEndB]$.
	 */
	public static final int projectDiffs(int diffs, int oldStartA, int oldEndA, int oldStartB, int oldEndB, int newStartA, int newEndA, int newStartB, int newEndB) {
		return (int)((((double)(diffs))/(oldEndA-oldStartA+oldEndB-oldStartB+2))*(newEndA-newStartA+newEndB-newStartB+2));
	}
	
	
	/**
	 * Builds the coverage histogram of each readA from an alignments file.
	 * The first number of each histogram is the read ID (starting from zero).
	 *
	 * @param inputPath output of LAshow;
	 * @param readLengths length of each readA in the alignments file;
	 * @param tmpArray temporary space, of size at least equal to the length of the 
	 * longest read in the file.
	 */
	public static final void getCoverageHistograms(String inputPath, String outputPath, int[] readLengths, int[] tmpArray) throws IOException {
		int i;
		int previousReadA;
		String str;
		BufferedReader br;
		BufferedWriter bw;

		br = new BufferedReader(new FileReader(inputPath),IO.BUFFER_SIZE);
		bw = new BufferedWriter(new FileWriter(outputPath),IO.BUFFER_SIZE);
		str=br.readLine(); str=br.readLine();  // Skipping the first two lines
		str=br.readLine();
		previousReadA=-1;
		while (str!=null) {
			readAlignmentFile(str);
			readA--;
			if (readA!=previousReadA) {
				if (previousReadA!=-1) {
					bw.write(previousReadA+",");
					for (i=0; i<readLengths[previousReadA]; i++) bw.write(tmpArray[i]+",");
					bw.write("\n");
				}
				Math.set(tmpArray,readLengths[readA]-1,0);
				previousReadA=readA;
			}
			if (startA<0) startA=0;
			endA--;  // The interval excludes position $endA$.
			if (endA>=readLengths[readA]) endA=readLengths[readA]-1;
			for (i=startA; i<=endA; i++) tmpArray[i]++;
			str=br.readLine();
		}
		bw.write(readA+",");
		for (i=0; i<readLengths[readA]; i++) bw.write(tmpArray[i]+",");
		bw.write("\n");
		br.close(); bw.close();
	}

}