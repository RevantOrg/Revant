package de.mpi_cbg.revant.factorize;

import java.util.Arrays;
import java.io.*;
import java.nio.charset.Charset;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Histograms;


public class Reads {
	/**
	 * Quality score constants, according to DBdump and DAmar. 
	 *
	 * Remark: quality scores are the opposite of quality. High score means high average
	 * dissimilarity between a B read and the A read, among a specific set of alignments, 
	 * i.e. low quality of the corresponding interval in the A read.
	 */
	public static final byte MIN_QUALITY_SCORE = 1;
	public static final byte MAX_QUALITY_SCORE = 51;
	public static final byte N_QUALITY_SCORES = MAX_QUALITY_SCORE-MIN_QUALITY_SCORE+1;
	public static final byte ASCII_FROM1 = 97;
	public static final byte ASCII_TO1 = 122;
	public static final byte ASCII_FROM2 = 65;
	public static final byte ASCII_TO2 = 89;
	public static final byte DAZZLER_THRESHOLD = 26;
	public static final byte DAMAR_OFFSET = 33;
	
	/**
	 * Parameters of the pipeline
	 */
	public static int maxReadLength;
	public static int minLongInsertionLength;
	public static final double MAX_QUALITY_STD = ((double)(MAX_QUALITY_SCORE-MIN_QUALITY_SCORE))/10;  // Regression tree intervals on a quality histogram are not split if their standard deviation is at most this much.
	public static byte MIN_RANDOM_QUALITY_SCORE;  // Minimum quality score for a substring to be considered a random insertion
	public static byte MAX_HIGH_QUALITY_SCORE;  // Maximum quality score for a substring to be considered of high quality
	public static int QUALITY_SPACING;  // DBdump returns one quality score for $QUALITY_SPACING$ consecutive positions
	public static double QUALITY_FRACTION = 0.5;  // Arbitrary

	/**
	 * Data structures
	 */
	public static int nReads;
	public static int[] readLengths;
	public static double[][] qualities;
	public static boolean virtualQualities;
	public static int[] readIDs = null;  // If not null: (1) contains the sorted set of IDs of reads in a subset of all reads; (2) $readLengths$ and $qualities$ refer just to the reads in $readIDs$, in order; (3) global variable $nReads$ equals the length of $readIDs$. 
	public static boolean readIDsAreCompact;  // TRUE=$readIDs$ contains a compact range.
	public static int firstRead, lastRead;
	
	/**
	 * Quality gradient variables
	 */
	public static final int GRADIENT_FROM = 0x0000FF00;
	public static final int GRADIENT_TO = 0x00FF0000;
	private static int rFrom, gFrom, bFrom, rTo, gTo, bTo;
	private static int rQuantum, gQuantum, bQuantum;

	/**
	 * Temporary space
	 */
	private static int[] tmp = new int[2];

	
	/**
	 * Loads a sorted set of $nr$ read IDs, which identify the reads in a subset of 
	 * the entire read set.
	 */
	public static final void loadReadIDs(String path, int nr) throws IOException {
		int i;
		BufferedReader br;
		
		nReads=nr;
		readIDs = new int[nReads];
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		for (i=0; i<nReads; i++) readIDs[i]=Integer.parseInt(br.readLine());
		br.close();
		firstRead=readIDs[0]; lastRead=readIDs[nReads-1];
		readIDsAreCompact=true;
		for (i=1; i<nReads; i++) {
			if (readIDs[i]!=readIDs[i-1]+1) {
				readIDsAreCompact=false;
				break;
			}
		}
	}
	
	
	/**
	 * Sets $readIDs$ to the compact range $[from..to]$ (zero-based).
	 */	
	public static final void loadReadIDs(int from, int to) {
		readIDs=null;
		nReads=to-from+1;
		firstRead=from; lastRead=to;
		readIDsAreCompact=true;
	}

	
	/**
	 * @return the position of $read$ in $readIDs$.
	 */
	public static final int indexOfRead(int read) {
		if (readIDsAreCompact) return read-firstRead;
		else return Arrays.binarySearch(readIDs,0,nReads,read);
	}

	
	/**
	 * @return $readIDs[i]$.
	 */
	public static final int readAtIndex(int i) {
		if (readIDsAreCompact) return firstRead+i;
		else return readIDs[i];
	}

	
	/**
	 * Remark: $path$ might contain the lengths of just a subset of size $nReads$ of all 
	 * reads.
	 *
	 * @return max length of a read.
	 */
	public static final int loadReadLengths(String path) throws IOException {
		int i;
		int length, maxLength;
		BufferedReader br;

		readLengths = new int[nReads];
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		maxLength=0;
		for (i=0; i<nReads; i++) {
			length=Integer.parseInt(br.readLine());
			readLengths[i]=length;
			maxLength=Math.max(maxLength,length);
		}
		br.close();
		return maxLength;
	}
	
	
	public static final int getReadLength(int readID) {
		if (readIDsAreCompact) return readLengths[readID-firstRead];
		int i = Arrays.binarySearch(readIDs,0,readIDs.length,readID);
		if (i<0) {
			System.err.println("Reads.getReadLength> ERROR: read ID not found in readIDs?! "+readID);
			System.exit(1);
		}
		return readLengths[i];
	}
	
	
	/**
	 * @return the sum of the length of every read.
	 */
	public static final long getTotalReadLength() {
		long out = 0;
		for (int i=0; i<nReads; i++) out+=readLengths[i];
		return out;
	}
	
	
	/**
	 * @return the total number of basepairs with quality score $<=maxScore$.
	 */
	public static final long getSurfaceAtQuality(byte maxScore) {
		int i, j;
		int length;
		long out;
		
		if (virtualQualities) return maxScore>=MAX_HIGH_QUALITY_SCORE?getTotalReadLength():0;
		out=0;
		for (i=0; i<nReads; i++) {
			length=qualities[i].length;
			for (j=0; j<length; j++) {
				if (qualities[i][j]<=maxScore) out+=QUALITY_SPACING;
			}
		}
		return out;
	}
	
	
	/**
	 * Stores in $out[0..X]$ the sorted sequence of maximal intervals of read $readID$ 
	 * with quality score $<=maxScore$, where $X$ is returned in output.
	 * 
	 * @param out assumed to be large enough.
	 */
	public static final int getQualityAsTrack(int readID, byte maxScore, int[] out) {
		int i, j;
		int currentStart, currentEnd, length;
		double[] qualityArray;
		
		qualityArray=getQualityArray(readID);
		length=qualityArray.length;
		currentStart=-1; currentEnd=-1;
		j=-1;
		for (i=0; i<length; i++) {
			if (qualityArray[i]<=maxScore) {
				if (currentStart==-1) currentStart=i;
				currentEnd=i;
			}
			else {
				if (currentStart!=-1) {
					out[++j]=currentStart*QUALITY_SPACING; 
					out[++j]=(currentEnd+1)*QUALITY_SPACING-1;
				}
				currentStart=-1; currentEnd=-1;
			}
		}
		if (currentStart!=-1) {
			out[++j]=currentStart*QUALITY_SPACING; 
			out[++j]=(currentEnd+1)*QUALITY_SPACING-1;
		}
		return j;
	}


	/**
	 * Loads the intrinsic quality values produced by DBdump.
	 *
	 * Remark: computing quality values for repeats is nontrivial in practice, since low 
	 * intrinsic quality just means few alignments covering a region: this might happen 
	 * either because the region is noise, or because the region is real but it was masked
	 * with a repeat track or a dust track. Using such masks might be unavoidable to get 
	 * the alignments in reasonable time and using a reasonable amount of disk space (for 
	 * LAS files), even when computing just 1x-versus-10x alignments. But this might mask
	 * entire repeats, and not only hide them from inference, but make the inference
	 * process assume that those repeats are instead regions of low quality. Another issue
	 * might be k-mer masking performed internally by the aligner, due to too many 
	 * occurrences of the same repeat: this might also remove several alignments and 
	 * affect intrinsic quality. The min length of an alignment should also be small 
	 * enough to allow seeding from short non-masked regions or from boundaries between
	 * adjacent masked regions.
	 *
	 * Remark: DASqv and DBdump produce one quality score for $QUALITY_SPACING$ 
	 * consecutive characters in a read. This score is actually a measure of error rather
	 * than of quality, i.e. high score means low quality.
	 *
	 * Remark: making a mistake by assigning low quality to a region with high quality
	 * makes the ends of repeat intervals non-maximal, and thus increases the size of
	 * clusters in the interval graph and the complexity of kernel assembly graphs (since
	 * it increases the probability that an edge is marked as containment or overlap).
	 * Assigning high quality to a region with low quality should ultimately yield smaller
	 * clusters and simpler (but more numerous) kernel graphs, oversplitting repeats.
	 *
	 * @param qualitiesFile if NULL, the procedure sets every cell of every read to the
	 * $MAX_HIGH_QUALITY_SCORE$ value stored in $qualityThresholdsFile$. $qualitiesFile$ 
	 * might contain the qualities of just a subset of size $nReads$ of all reads.
	 */
	public static final void loadQualities(String qualityThresholdsFile, String qualitiesFile) throws IOException {
		byte mode;
		int i, j;
		int readID;
		String str;
		BufferedReader br;
		byte[] tmpBytes;

		// Loading quality thresholds
		loadQualities_thresholds(qualityThresholdsFile);

		// Loading read qualities
		if (qualitiesFile!=null) {
			virtualQualities=false;
			qualities = new double[nReads][1];
			mode=filename2qualityMode(qualitiesFile);
			br = new BufferedReader(new FileReader(qualitiesFile),IO.BUFFER_SIZE);
			for (i=0; i<nReads; i++) {
				str=br.readLine();
				readID=readIDsAreCompact?firstRead+i:readIDs[i];
				if (str==null || str.length()<Math.ceil(getReadLength(readID),QUALITY_SPACING)) {
					IO.printCriticalErr("Malformed quality file: read "+readID+", length="+getReadLength(readID)+" qualityLength="+(str==null?"null":str.length()));
					System.exit(1);
				}
				tmpBytes=str.getBytes(Charset.forName("US-ASCII"));
				qualities[i] = new double[tmpBytes.length];
				for (j=0; j<tmpBytes.length; j++) qualities[i][j]=ascii2quality(tmpBytes[j],mode);
			}
			br.close();
		}
		else {
			virtualQualities=true;
			qualities = new double[1][Math.ceil(maxReadLength,QUALITY_SPACING)];
			for (j=0; j<qualities[0].length; j++) qualities[0][j]=MAX_HIGH_QUALITY_SCORE;
		}
		
		// Loading colors
		loadQualities_colors();
	}
	
	
	public static final void loadQualities_thresholds(String qualityThresholdsFile) throws IOException {
		int i;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(qualityThresholdsFile),IO.BUFFER_SIZE);
		str=br.readLine();
		i=str.indexOf(" ");
		MIN_RANDOM_QUALITY_SCORE=Byte.parseByte(str.substring(0,i));
		str=br.readLine();
		i=str.indexOf(" ");
		MAX_HIGH_QUALITY_SCORE=Byte.parseByte(str.substring(0,i));
		str=br.readLine();
		i=str.indexOf(" ");
		QUALITY_SPACING=Integer.parseInt(str.substring(0,i));
		br.close();
	}
	
	
	private static final void loadQualities_colors() {
		int mask = 0x000000FF;
		
		bFrom=GRADIENT_FROM&mask;
		mask<<=8;
		gFrom=(GRADIENT_FROM&mask)>>8;
		mask<<=8;
		rFrom=(GRADIENT_FROM&mask)>>16;
		mask=0x000000FF;
		bTo=GRADIENT_TO&mask;
		mask<<=8;
		gTo=(GRADIENT_TO&mask)>>8;
		mask<<=8;
		rTo=(GRADIENT_TO&mask)>>16;
		rQuantum=(rTo-rFrom+1)/N_QUALITY_SCORES;
		gQuantum=(gTo-gFrom+1)/N_QUALITY_SCORES;
		bQuantum=(bTo-bFrom+1)/N_QUALITY_SCORES;
	}

	
	/**
	 * Transforms a character into a quality value in [1..51].
	 *
	 * @param mode 0=uses the convention of DASqv and DBdump; 1=uses the convention of
	 * DAmar;
	 * @return -1 if the character is not the valid encoding of a quality value.
	 */																	 
	private static final byte ascii2quality(byte c, byte mode) {
		if (mode==0) {
			if (c>=ASCII_FROM1 && c<=ASCII_TO1) return (byte)(c-ASCII_FROM1+1);
			else if (c>=ASCII_FROM2 && c<=ASCII_TO2) return (byte)((ASCII_TO1-ASCII_FROM1+1)+(c-ASCII_FROM2)+1);
			else return -1;
		}
		else if (mode==1) {
			final byte value = (byte)(c-DAMAR_OFFSET);
			return value>MAX_QUALITY_SCORE?MAX_QUALITY_SCORE:value;  // Zero (=no covering alignment) is mapped to 56.
		}
		else return -1;
	}
	
	
	/**
	 * Inverse of $ascii2quality()$.
	 */
	private static final byte quality2ascii(byte value, byte mode) {
		if (mode==0) {
			if (value<=DAZZLER_THRESHOLD) return (byte)(ASCII_FROM1-1+value);
			else if (value<=MAX_QUALITY_SCORE) return (byte)(ASCII_FROM2+value-DAZZLER_THRESHOLD-1);
			else return ASCII_TO2;
		}
		else if (mode==1) return (byte)(value+DAMAR_OFFSET);
		else return -1;
	}
	
	
	/**
	 * Infers whether a quality file is in the DBDUMP or DAMAR format, from its 
	 * extension.
	 */
	private static final byte filename2qualityMode(String filename) {
		final int p = filename.lastIndexOf('.');
		if (p<0) return 0;
		if (filename.substring(p).equalsIgnoreCase(IO.DBDUMP_QUALITY_EXTENSION)) return 0;
		else if (filename.substring(p).equalsIgnoreCase(IO.DAMAR_QUALITY_EXTENSION)) return 1;
		return 0;
	}
	
	
	/**
	 * Remark: the procedure uses gradient variables initialized in $loadQualities()$.
	 */
	public static final int quality2color(int quality) {
		final int OFFSET = quality-MIN_QUALITY_SCORE+1;
		final int MASK = 0x000000FF;
		int out;
		
		out=0x00000000;
		out=(bFrom+OFFSET*bQuantum)&MASK;
		out|=((gFrom+OFFSET*gQuantum)&MASK)<<8;
		out|=((rFrom+OFFSET*rQuantum)&MASK)<<16;
		return out;
	}

	
	/**
	 * Converts an interval $[start..end]$ in a read, into an interval in the quality 
	 * track of the read.
	 */
	public static final void readEnds2qualityEnds(int readID, int start, int end, boolean orientation, int[] tmpArray) {
		final int length = getReadLength(readID);
		int startPrime, endPrime;
		
		if (start<0) start=0;
		else if (start>=length) start=length-1;
		if (end<start) end=start;
		else if (end>=length) end=length-1;
		if (!orientation) {
			startPrime=length-end-1;
			endPrime=length-start-1;
			start=startPrime; end=endPrime;
		}
		start=Math.round(start,QUALITY_SPACING);
		end=Math.round(end,QUALITY_SPACING)-1;
		if (end<start) end=start;
		tmpArray[0]=start; tmpArray[1]=end;
	}


	/**
	 * @return true if at least $QUALITY_FRACTION$ fraction of the substring has quality
	 * score at least $MIN_RANDOM_QUALITY_SCORE$, or if the average quality score in the 
	 * substring is at least $MIN_RANDOM_QUALITY_SCORE$.
	 */
	public static final boolean hasLowQuality(int readID, int start, int end, boolean orientation) {
		readEnds2qualityEnds(readID,start,end,orientation,tmp);
		start=tmp[0]; end=tmp[1];
		return Histograms.getFraction(getQualityArray(readID),start,end,MIN_RANDOM_QUALITY_SCORE,true)>=QUALITY_FRACTION ||
			   Histograms.getAverage(getQualityArray(readID),start,end)>=MIN_RANDOM_QUALITY_SCORE;
	}


	/**
	 * @return true iff the specified substring of the read has low quality according to 
	 * procedure $hasLowQuality()$, and if it does not contain a long substring with 
	 * average quality score at most $MAX_HIGH_QUALITY_SCORE$.
	 */
	public static final boolean isRandomInsertion(int readID, int start, int end, boolean orientation) {
		final int WINDOW = 4;  // Multiple of $Reads.QUALITY_SPACING$. Arbitrary.
		
        if (!hasLowQuality(readID,start,end,true)) return false;
		readEnds2qualityEnds(readID,start,end,orientation,tmp);
		start=tmp[0]; end=tmp[1];
	    return end-start+1<=WINDOW*Reads.QUALITY_SPACING?true:!Histograms.hasSubstringWithAverage(getQualityArray(readID),start,end,MAX_HIGH_QUALITY_SCORE+1,false,WINDOW,true,-1);
	}


	/**
	 * @param orientation TRUE=forward, FALSE=reverse.
	 */
	public static final double getAverageQuality(int readID, int start, int end, boolean orientation) {
		readEnds2qualityEnds(readID,start,end,orientation,tmp);
		start=tmp[0]; end=tmp[1];
		return Histograms.getAverage(getQualityArray(readID),start,end);
	}

	
	public static final boolean isLeftMaximal(int start, int readID, boolean orientation) {
		final int WINDOW = 1;  // Multiple of $Reads.QUALITY_SPACING$. Arbitrary.
		return start>IO.quantum && !hasLowQuality(readID,start-WINDOW*Reads.QUALITY_SPACING,start-1,orientation);
	}
	

	public static final boolean isRightMaximal(int end, int readID, boolean orientation) {
		final int WINDOW = 1;  // Multiple of $Reads.QUALITY_SPACING$. Arbitrary.
		return getReadLength(readID)-end-1>IO.quantum && !hasLowQuality(readID,end+1,end+WINDOW*Reads.QUALITY_SPACING,orientation);
	}
	
	
	public static final boolean isLeftRightMaximal(int pos, int readID, boolean orientation) {
		final int WINDOW = 1;  // Multiple of $Reads.QUALITY_SPACING$. Arbitrary.
		return pos>IO.quantum && getReadLength(readID)-pos-1>IO.quantum && !hasLowQuality(readID,pos-WINDOW*Reads.QUALITY_SPACING,pos+WINDOW*Reads.QUALITY_SPACING,orientation);
	}
	
	
	/**
	 * Removes low-quality regions from the start/end of an interval of readA.
	 *
	 * @param range input interval;
	 * @param window a multiple of $QUALITY_SPACING$ over which average quality is 
	 * computed at the beginning/end of the interval;
	 * @return the new interval, stored again in $range$, or pair (-1,-1) if no high-
	 * quality window can be found. The procedure does not check whether the output is a 
	 * valid interval, i.e. the output might have start > end.
	 */
	public static final void trim(int[] range, int window, boolean trimLeft, boolean trimRight) {
		final int readID = ReadA.id;
		readEnds2qualityEnds(readID,range[0],range[1],true,tmp);
		final int start = tmp[0]; 
		final int end = tmp[1];
		int p;
		double quality;
		
		if (trimLeft) {
			quality=Histograms.getAverage(getQualityArray(readID),start,start+window-1);
			if (quality>=MIN_RANDOM_QUALITY_SCORE) {
				p=Histograms.skipRight(getQualityArray(readID),getQualityArrayLength(readID),true,start,window,MIN_RANDOM_QUALITY_SCORE-1,false);
				if (p==-1) {
					range[0]=-1; range[1]=-1;
					return;
				}
				range[0]=p*QUALITY_SPACING;
			}
		}
		if (trimRight) {
			quality=Histograms.getAverage(getQualityArray(readID),end-window+1,end);
			if (quality>=MIN_RANDOM_QUALITY_SCORE) {
				p=Histograms.skipLeft(getQualityArray(readID),getQualityArrayLength(readID),true,end-window+1,window,MIN_RANDOM_QUALITY_SCORE-1,false);
				if (p==-1) {
					range[0]=-1; range[1]=-1;
					return;
				}
				range[1]=(p+window)*QUALITY_SPACING;
			}
		}
	}
	
	
	public static final double[] getQualityArray(int readID) {
		if (virtualQualities) return qualities[0];
		if (!readIDsAreCompact) {
			int i = Arrays.binarySearch(readIDs,0,readIDs.length,readID);
			if (i<0) {
				System.err.println("getQualityArray> ERROR: readID="+readID+" not found?!");
				System.exit(1);
			}
			return qualities[i];
		}
		else return qualities[readID-firstRead];
	}
	
	
	public static final int getQualityArrayLength(int readID) {
		return Math.ceil(getReadLength(readID),QUALITY_SPACING);
	}
	
	
	/**
	 * Stores in $out[0..X]$ the sorted sequence of maximal disjoint intervals with 
	 * quality>=minLowQuality in substring $[from..to]$ or $readID$; every interval is 
	 * represented as a pair of positions in the read (i.e. not in the quality array). 
	 * X is returned in output.
	 *
	 * Remark: the procedure assumes tha QUALITY_SPACING has already been initialized.
	 *
	 * @param out assumed to be large enough.
	 */
	public static final int getLowQualityIntervals(int readID, byte minLowQuality, int from, int to, int[] out) {
		int i;
		int start, end, last;
		final double[] qualityArray = getQualityArray(readID);
		
		start=-1; end=-1; last=-1;
		for (i=0; i<qualityArray.length; i++) {
			if (qualityArray[i]<minLowQuality) {
				if (start!=-1 && Intervals.intersect(start,end,from,to)) { 
					out[++last]=Math.max(start,from);
					out[++last]=Math.min(end,to); 
				}
				start=-1; end=-1;
			}
			else {
				if (start==-1) start=i*QUALITY_SPACING;
				end=(i+1)*QUALITY_SPACING-1;
			}
		}
		if (start!=-1 && Intervals.intersect(start,end,from,to)) { 
			out[++last]=Math.max(start,from);
			out[++last]=Math.min(end,to);
		}
		return last;
	}

	
	/**
	 * Prints the intersection and symmetric difference of all values $>=threshold$ in two
	 * PHRED files (which are assumed to contain one sequence of PHRED characters per
	 * line).
	 *
	 * @param threshold in [1..51];
	 * @param verbose TRUE=outputs to STDERR all pairs of reads with intersection<=50%.
	 */
	public static final void comparePhreds(String phredFile1, String phredFile2, int threshold, boolean verbose) throws IOException {
		final double VERBOSE_INTERSECTION = 0.5;  // Arbitrary
		byte mode1, mode2;
		int i;
		int value1, value2, min;
		long oneMinusTwo, oneIntersectTwo, twoMinusOne, union, totalLength;
		long localOneMinusTwo, localOneIntersectTwo, localTwoMinusOne, localUnion, minLocalIntersection;
		String str1, str2, minStr1, minStr2;
		BufferedReader br1, br2;
		byte[] tmpBytes1, tmpBytes2;
		
		mode1=filename2qualityMode(phredFile1);
		mode2=filename2qualityMode(phredFile2);
		br1 = new BufferedReader(new FileReader(phredFile1));
		br2 = new BufferedReader(new FileReader(phredFile2));
		str1=br1.readLine(); str2=br2.readLine();
		oneMinusTwo=0L; oneIntersectTwo=0L; twoMinusOne=0L; totalLength=0L;
		minLocalIntersection=Math.POSITIVE_INFINITY; minStr1=null; minStr2=null;
		while (str1!=null && str2!=null) {
			if (str1.length()!=str2.length())  {
				System.err.println("WARNING> strings of different lengths:  "+str1.length()+"!="+str2.length());
				System.err.println(str1);
				System.err.println(str2);
				min=Math.min(str1.length(),str2.length());
				if (str1.length()!=min) str1=str1.substring(0,min);
				if (str2.length()!=min) str2=str2.substring(0,min);
			}
			tmpBytes1=str1.getBytes(Charset.forName("US-ASCII"));
			tmpBytes2=str2.getBytes(Charset.forName("US-ASCII"));
			totalLength+=tmpBytes1.length;
			localOneMinusTwo=0L; localOneIntersectTwo=0L; localTwoMinusOne=0L;
			for (i=0; i<tmpBytes1.length; i++) {
				value1=ascii2quality(tmpBytes1[i],mode1);
				value2=ascii2quality(tmpBytes2[i],mode2);
				if (value1>=threshold) {
					if (value2>=threshold) localOneIntersectTwo++;
					else localOneMinusTwo++;
				}
				else if (value2>=threshold) localTwoMinusOne++;
			}
			localUnion=localOneMinusTwo+localOneIntersectTwo+localTwoMinusOne;
			if (verbose && localUnion>0 && localOneIntersectTwo<=localUnion*VERBOSE_INTERSECTION) {
				System.err.println("INTERSECTION<="+VERBOSE_INTERSECTION);
				System.err.println("str1: "+str1);
				System.err.println("str2: "+str2);
				System.err.println();
			}
			if (localUnion>0 && localOneIntersectTwo<minLocalIntersection) {
				minLocalIntersection=localOneIntersectTwo;
				minStr1=str1; minStr2=str2;
			}
			oneMinusTwo+=localOneMinusTwo; oneIntersectTwo+=localOneIntersectTwo; twoMinusOne+=localTwoMinusOne;			
			str1=br1.readLine(); str2=br2.readLine();
		}
		br1.close(); br2.close();
		union=oneMinusTwo+oneIntersectTwo+twoMinusOne;
		System.out.println("PHRED scores >="+threshold+": union="+IO.getPercent(union,totalLength)+"% of all values");
		System.out.println("Intersection="+IO.getPercent(oneIntersectTwo,union)+"% of union");
		System.out.println("Quality1 \\setminus Quality2="+IO.getPercent(oneMinusTwo,union)+"% of union");
		System.out.println("Quality2 \\setminus Quality1="+IO.getPercent(twoMinusOne,union)+"% of union");
		System.out.println();
		System.out.println("Read with minimum intersection: ");
		System.err.println(minStr1);
		System.err.println(minStr2);
	}
	
	
	/**
	 * Stores in $histogram$ the number of times each quality value occurs in 
	 * $qualityFile$.
	 */
	public static final void qualityHistogram(String qualityFile, long[] histogram) throws IOException {
		final byte mode = filename2qualityMode(qualityFile);
		final int histogramLength = histogram.length;
		int i;
		int length, value;
		String str;
		BufferedReader br;
		byte[] tmpBytes;
		
		Math.set(histogram,histogram.length-1,0);
		br = new BufferedReader(new FileReader(qualityFile));
		str=br.readLine();
		while (str!=null) {
			tmpBytes=str.getBytes(Charset.forName("US-ASCII"));
			length=tmpBytes.length;
			for (i=0; i<length; i++) {
				value=ascii2quality(tmpBytes[i],mode);
				histogram[value>=histogramLength?histogramLength-1:value]++;
			}
			str=br.readLine();
		}
		br.close();
	}
	
	
	/**
	 * Prints the intersection and symmetric difference of sets Q and T, where Q is the
	 * set of all maximal intervals of consecutive values $>=qualityThreshold$ in 
	 * $phredFile$, and T is the set built from $trackFile$ as follows: every maximal 
	 * sequence of intervals at distance $<=distanceThreshold$ from one another is merged
	 * to a single interval.
	 *
	 * Remark: the procedure assumes that QUALITY_SPACING has already been initialized.
	 *
	 * @param phredFile assumed to contain one sequence of PHRED characters per line;
	 * @param trackFile assumed to contain a sorted sequence of intervals per line (or an
	 * empty line if no track is present in a read);
	 * @param qualityThreshold in [1..51].
	 */
	public static final void comparePhred2track(String phredFile, String trackFile, int qualityThreshold, int distanceThreshold) throws IOException {
		final int MAX_READ_LENGTH = 1000000;
		byte mode;
		int i;
		int sum, value, start, end, currentStart, currentEnd, lastQualityInterval, lastTrackInterval, lastDifferenceInterval;
		long qualityMinusTrack, trackMinusQuality, intersection, totalLength, qualityUnion;
		String str1, str2;
		BufferedReader br1, br2;
		byte[] tmpBytes;
		int[] qualityIntervals = new int[Math.ceil(MAX_READ_LENGTH,QUALITY_SPACING)];
		int[] trackIntervals = new int[MAX_READ_LENGTH];
		int[] differenceIntervals = new int[MAX_READ_LENGTH*3];
		String[] tokens;
		
		mode=filename2qualityMode(phredFile);
		br1 = new BufferedReader(new FileReader(phredFile));
		br2 = new BufferedReader(new FileReader(trackFile));
		qualityMinusTrack=0L; trackMinusQuality=0L; intersection=0L; totalLength=0L;
		str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null) {
			// Loading low-quality intervals
			tmpBytes=str1.getBytes(Charset.forName("US-ASCII"));
			totalLength+=tmpBytes.length*QUALITY_SPACING;
			lastQualityInterval=-1; currentStart=-1; currentEnd=-1;
			for (i=0; i<tmpBytes.length; i++) {
				value=ascii2quality(tmpBytes[i],mode);
				if (value>=qualityThreshold) {
					if (currentStart==-1) {
						currentStart=i*QUALITY_SPACING;
						currentEnd=(i+1)*QUALITY_SPACING-1;
					}
					else currentEnd=(i+1)*QUALITY_SPACING-1;
				}
				else if (currentStart!=-1) {
					qualityIntervals[++lastQualityInterval]=currentStart;
					qualityIntervals[++lastQualityInterval]=currentEnd;
					currentStart=-1; currentEnd=-1;
				}
			}
			if (currentStart!=-1) {
				qualityIntervals[++lastQualityInterval]=currentStart;
				qualityIntervals[++lastQualityInterval]=currentEnd;
			}			
			// Loading track intervals
			if (str2.length()==0) lastTrackInterval=-1;
			else {
				lastTrackInterval=-1; currentStart=-1; currentEnd=-1;
				tokens=str2.split(" ");
				for (i=0; i<tokens.length; i+=2) {
					start=Integer.parseInt(tokens[i]);
					end=Integer.parseInt(tokens[i+1]);
					if (currentStart==-1) { currentStart=start; currentEnd=end; }
					else {
						if (start<=currentEnd+distanceThreshold) currentEnd=end;
						else {
							trackIntervals[++lastTrackInterval]=currentStart;
							trackIntervals[++lastTrackInterval]=currentEnd;
							currentStart=start; currentEnd=end;
						}
					}
				}
				if (currentStart!=-1) {
					trackIntervals[++lastTrackInterval]=currentStart;
					trackIntervals[++lastTrackInterval]=currentEnd;
				}
			}			
			// Comparing the two sets of intervals
			intersection+=Intervals.intersectionLength(qualityIntervals,lastQualityInterval,trackIntervals,lastTrackInterval,1,true);
			lastDifferenceInterval=Intervals.difference(qualityIntervals,lastQualityInterval,trackIntervals,lastTrackInterval,differenceIntervals,0,-1,1);
			sum=0;
			for (i=0; i<lastDifferenceInterval; i+=3) sum+=differenceIntervals[i+2]-differenceIntervals[i+1]+1;
			qualityMinusTrack+=sum;
			lastDifferenceInterval=Intervals.difference(trackIntervals,lastTrackInterval,qualityIntervals,lastQualityInterval,differenceIntervals,0,-1,1);
			sum=0;
			for (i=0; i<lastDifferenceInterval; i+=3) sum+=differenceIntervals[i+2]-differenceIntervals[i+1]+1;
			trackMinusQuality+=sum;
			str1=br1.readLine(); str2=br2.readLine(); 
		}
		qualityUnion=qualityMinusTrack+intersection;
		System.out.println("PHRED scores >="+qualityThreshold+": union="+IO.getPercent(qualityUnion,totalLength)+"% of all basepairs");
		System.out.println("LowQuality \\intersect Track = "+IO.getPercent(intersection,qualityUnion)+"% of all low-quality");
		System.out.println("Track \\setminus LowQuality = "+IO.getPercent(trackMinusQuality,intersection+trackMinusQuality)+"% of all track");
	}
	
	
	/**
	 * Prints the intersection of sets T1 and T2, where TX is the set built from 
	 * $trackFileX$ as follows: every maximal sequence of intervals at distance $<=
	 * distanceThresholdX$ from one another is merged into a single interval.
	 *
	 * @param trackFileX assumed to contain a sorted sequence of intervals per line (or an
	 * empty line if no track is present in a read).
	 */
	public static final void compareTrack2track(String trackFile1, String trackFile2, int distanceThreshold1, int distanceThreshold2) throws IOException {
		final int MAX_READ_LENGTH = 1000000;
		int i;
		int sum, start, end, currentStart, currentEnd, lastInterval1, lastInterval2, lastDifferenceInterval, readID;
		long oneMinusTwo, twoMinusOne, intersection, union1;
		String str1, str2;
		BufferedReader br1, br2;
		int[] intervals1 = new int[MAX_READ_LENGTH];
		int[] intervals2 = new int[MAX_READ_LENGTH];
		int[] differenceIntervals = new int[MAX_READ_LENGTH*3];
		String[] tokens;
		
		br1 = new BufferedReader(new FileReader(trackFile1));
		br2 = new BufferedReader(new FileReader(trackFile2));
		oneMinusTwo=0L; twoMinusOne=0L; intersection=0L; readID=0;
		str1=br1.readLine(); str2=br2.readLine();
		while (str1!=null) {
			readID++;  // 1-based
			// Loading track 1 intervals
			if (str1.length()==0) lastInterval1=-1;
			else {
				lastInterval1=-1; currentStart=-1; currentEnd=-1;
				tokens=str1.split(" ");
				for (i=0; i<tokens.length; i+=2) {
					start=Integer.parseInt(tokens[i]);
					end=Integer.parseInt(tokens[i+1]);
					if (currentStart==-1) { currentStart=start; currentEnd=end; }
					else {
						if (start<=currentEnd+distanceThreshold1) currentEnd=end;
						else {
							intervals1[++lastInterval1]=currentStart;
							intervals1[++lastInterval1]=currentEnd;
							currentStart=start; currentEnd=end;
						}
					}
				}
				if (currentStart!=-1) {
					intervals1[++lastInterval1]=currentStart;
					intervals2[++lastInterval1]=currentEnd;
				}
			}
			// Loading track 2 intervals
			if (str2.length()==0) lastInterval2=-1;
			else {
				lastInterval2=-1; currentStart=-1; currentEnd=-1;
				tokens=str2.split(" ");
				for (i=0; i<tokens.length; i+=2) {
					start=Integer.parseInt(tokens[i]);
					end=Integer.parseInt(tokens[i+1]);
					if (currentStart==-1) { currentStart=start; currentEnd=end; }
					else {
						if (start<=currentEnd+distanceThreshold2) currentEnd=end;
						else {
							intervals2[++lastInterval2]=currentStart;
							intervals2[++lastInterval2]=currentEnd;
							currentStart=start; currentEnd=end;
						}
					}
				}
				if (currentStart!=-1) {
					intervals2[++lastInterval2]=currentStart;
					intervals2[++lastInterval2]=currentEnd;
				}
			}			
			// Comparing the two sets of intervals
			intersection+=Intervals.intersectionLength(intervals1,lastInterval1,intervals2,lastInterval2,1,true);
			lastDifferenceInterval=Intervals.difference(intervals1,lastInterval1,intervals2,lastInterval2,differenceIntervals,0,-1,1);
			sum=0;
			for (i=0; i<lastDifferenceInterval; i+=3) sum+=differenceIntervals[i+2]-differenceIntervals[i+1]+1;
			oneMinusTwo+=sum;
			lastDifferenceInterval=Intervals.difference(intervals2,lastInterval2,intervals1,lastInterval1,differenceIntervals,0,-1,1);
			sum=0;
			for (i=0; i<lastDifferenceInterval; i+=3) sum+=differenceIntervals[i+2]-differenceIntervals[i+1]+1;
			twoMinusOne+=sum;
			str1=br1.readLine(); str2=br2.readLine(); 
		}
		System.out.println("Track1 \\intersect Track2 = "+IO.getPercent(intersection,oneMinusTwo+intersection)+"% of all Track1 and "+IO.getPercent(intersection,twoMinusOne+intersection)+"% of all Track2");
	}
	
	
	/**
	 * Given the qualities in $phredFile$, the procedure writes to $outputFile$ a copy
     * of them such that every value $>=minLowQuality$ is replaced with 
	 * $MIN_RANDOM_QUALITY_SCORE-1$ iff: (1) it belongs to the tandem track; (2) it does
	 * not intersect any interval in the (non-generalized, raw) dust track.
	 *
	 * Remark: the proceure assumes $QUALITY_SPACING$ to have already been initialized.
	 *
	 * @param minBadQuality a value, not a character;
	 * @param dustFile,tandemFile assumed to contain a sorted sequence of intervals per 
	 * line (or an empty line if no track is present in a read).
	 */
	public static final void addTracksToPhred(String phredFile, byte minLowQuality, String qualityThresholdsFile, String dustFile, String tandemFile, String outputFile) throws IOException {
		final int TANDEM_THRESHOLD = IO.quantum;  // Arbitrary
		final int MIN_TANDEM_INTERSECTION = (int)(QUALITY_SPACING*0.6);  // Arbitrary
		final int MAX_READ_LENGTH = 1000000;  // Arbitrary
		boolean intersectsDust, inTandem;
		byte mode, value, newGoodQuality;
		int i, j, j2, j3;
		int readID, start, end, length, lastDustInterval, lastTandemInterval, firstJ3ForNextI;
		long nLowQuality, nLowQualityInDust, nLowQualityInTandem, nConverted;
		String str1, str2, str3;
		BufferedReader br1, br2, br3;
		BufferedWriter bw;
		byte[] qualities;
		int[] dustIntervals = new int[MAX_READ_LENGTH];
		int[] tandemIntervals = new int[MAX_READ_LENGTH];
		String[] tokens;
		
		mode=filename2qualityMode(phredFile);
		loadQualities_thresholds(qualityThresholdsFile);
		newGoodQuality=quality2ascii((byte)(MIN_RANDOM_QUALITY_SCORE-1),mode);
		br1 = new BufferedReader(new FileReader(phredFile));
		br2 = new BufferedReader(new FileReader(dustFile));
		br3 = new BufferedReader(new FileReader(tandemFile));
		bw = new BufferedWriter(new FileWriter(outputFile));
		str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine(); 
		nLowQuality=0L; nLowQualityInDust=0L; nLowQualityInTandem=0L; nConverted=0L; 
		readID=0;
		while (str1!=null) {
			readID++;  // 1-based
			qualities=str1.getBytes(Charset.forName("US-ASCII"));
			lastDustInterval=-1;
			if (str2.length()>0) {
				tokens=str2.split(" ");
				for (i=0; i<tokens.length; i++) dustIntervals[++lastDustInterval]=Integer.parseInt(tokens[i]);
			}
			lastTandemInterval=-1;
			if (str3.length()>0) {
				tokens=str3.split(" ");
				for (i=0; i<tokens.length; i++) tandemIntervals[++lastTandemInterval]=Integer.parseInt(tokens[i]);
			}
			if (lastTandemInterval!=-1) {
				start=tandemIntervals[0]; end=tandemIntervals[1];
				j=-1;
				for (i=2; i<lastTandemInterval; i+=2) {
					if (tandemIntervals[i]<=end+TANDEM_THRESHOLD) end=tandemIntervals[i+1];
					else { 
						tandemIntervals[++j]=start; tandemIntervals[++j]=end; 
						start=tandemIntervals[i]; end=tandemIntervals[i+1];
					}
				}
				tandemIntervals[++j]=start; tandemIntervals[++j]=end;
				lastTandemInterval=j;
			}
			length=qualities.length;
			j2=0; j3=0; firstJ3ForNextI=-1;
			for (i=0; i<length; i++) {
				value=ascii2quality(qualities[i],mode);
				if (value<minLowQuality) {
					bw.write(str1.charAt(i));
					continue;
				}
				nLowQuality++;
				start=i*QUALITY_SPACING; end=(i+1)*QUALITY_SPACING-1;
				intersectsDust=false;
				while (j2<lastDustInterval && dustIntervals[j2+1]<start) j2+=2;
				if (j2<lastDustInterval && dustIntervals[j2]<end) intersectsDust=true;
				if (intersectsDust) nLowQualityInDust++;
				inTandem=false;
				while (j3<lastTandemInterval && tandemIntervals[j3+1]<start) j3+=2;
				firstJ3ForNextI=-1;
				while (j3<lastTandemInterval) {
					if (tandemIntervals[j3]>=end) break;
					if (firstJ3ForNextI==-1 && i<length-1 && tandemIntervals[j3+1]>=(i+1)*QUALITY_SPACING) firstJ3ForNextI=j3;
					if ( Intervals.isContained(start,end,tandemIntervals[j3],tandemIntervals[j3+1]) || 
						 Intervals.intersectionLength(start,end,tandemIntervals[j3],tandemIntervals[j3+1])>=MIN_TANDEM_INTERSECTION
					   ) {
						inTandem=true;
						break;
					}
					j3+=2;
				}
				if (inTandem) {
					nLowQualityInTandem++;
					if (intersectsDust) bw.write(str1.charAt(i));
					else {
						bw.write(newGoodQuality);
						nConverted++;
					}
				}
				else bw.write(str1.charAt(i));
				if (firstJ3ForNextI!=-1) j3=firstJ3ForNextI;
			}
			bw.write("\n");
			str1=br1.readLine(); str2=br2.readLine(); str3=br3.readLine();
		}
		br1.close(); br2.close(); br3.close(); bw.close();
		System.err.println("addTracksToPhred> Converted "+nConverted+" low-quality values out of "+nLowQuality+" ("+IO.getPercent(nConverted,nLowQuality)+"%)");
		System.err.println("addTracksToPhred> Low-quality values that are also dust: "+nLowQualityInDust+" ("+IO.getPercent(nLowQualityInDust,nLowQuality)+"%)");
		System.err.println("addTracksToPhred> Low-quality values that are also tandem: "+nLowQualityInTandem+" ("+IO.getPercent(nLowQualityInTandem,nLowQuality)+"%)");
	}

}