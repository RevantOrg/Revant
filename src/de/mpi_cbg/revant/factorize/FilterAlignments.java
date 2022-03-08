package de.mpi_cbg.revant.factorize;

import java.io.*;
import java.util.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;

/**
 * An aligner typically returns alignments with error rate inside a wide range. Even when  
 * explicitly asked to just report alignments with error <=E, an aligner might return 
 * alignments with error way above E: for example, we observed that the alignments 
 * returned by daligner when called with max error rate 0.25, have avg. error rate in a 
 * Gaussian that peaks at 0.25, but that has approx. half of its mass in [0.25..0.35].
 * The same applies to min alignment length: in particular, the length of an alignment 
 * might be shorter than what was asked to the aligner, in just one read.
 * This program discards every alignment with avg error > E or with min length over the 
 * two reads < L.
 *
 * Remark: it is not yet clear what is the optimal setting of E for repeat inference. 
 * Setting E to the avg. error rate can discard half of all alignments, so some repeats 
 * might not get detected; on the other hand, the repeats we detect would likely contain 
 * only fragments that come from the same repeat variant, and their consensus sequences 
 * could be very accurate; the graphs used in repeat inference should also be simpler. 
 * Setting E too large would create more complex graphs that connect several variants of 
 * the same repeat, the fragments assigned to a repeat might come from different variants,
 * and the consensus might be less accurate; on the other hand, most of the repeats should
 * be detected.
 */
public class FilterAlignments  {
	
	/**
	 * Repeat track, if any.
	 */
	private static int[] lastRepeatTrack;
	private static int[][] repeatTrack;

	/**
	 * Temporary space
	 */
	private static int[] tmpArray1, tmpArray2, tmpArray3, tmpArray4, tmpArray5;
	private static long tmpSurface;
	private static String tmpString;
	private static boolean VERBOSE;

	
	/**
	 * Remark: the program can work on the alignments of a single block.
	 *
	 * @param args 
	 *  4: number of reads in the read lengths file;
	 *  5: max length of a read in the read lengths file;
	 *  8: directory that contains all repeat track files (use "null" to disregard);
	 *  9: common prefix of all repeat track files (use "null" to disregard); 
	 *     every file is assumed to contain a distinct repeat track; only alignments that
	 *     map a large-enough contiguous interval outside any repeat track, to a large-
	 *     enough contiguous interval outside any repeat track, are kept (details in 
	 *     $hasNonrepeatSurface()$);
	 * 11: first readID in the block (zero-based);
	 * 12: last readID in the block (zero-based);
	 * 13: number of rows in the alignments file (not the n. of distinct alignments);
	 * 15: file containing the quality of each read (use "null" to disregard).
	 */
	public static void main(String[] args) throws IOException {
		final double MAX_ERROR = Double.parseDouble(args[0]);
		final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[1]);
		final String INPUT_FILE = args[2];
		final String OUTPUT_FILE = args[3];
		final int N_READS = Integer.parseInt(args[4]);
		final int MAX_READ_LENGTH = Integer.parseInt(args[5]);
		final String READ_LENGTHS_FILE = args[6];
		final String READ_IDS_FILE = args[7];
		final String REPEAT_TRACK_DIR = args[8];
		final String REPEAT_TRACK_PREFIX = args[9];
		final byte REPEAT_TRACK_FORMAT = Byte.parseByte(args[10]);
		final int FIRST_READ = Integer.parseInt(args[11]);
		final int LAST_READ = Integer.parseInt(args[12]);
		final int N_ALIGNMENTS = Integer.parseInt(args[13]);
		final String QUALITY_THRESHOLDS_FILE = args[14];
		final String QUALITIES_FILE = args[15].equalsIgnoreCase("null")?null:args[15];
		VERBOSE=Integer.parseInt(args[16])==1;
		
		final boolean REPEAT_TRACK_EXISTS = !REPEAT_TRACK_DIR.equalsIgnoreCase("null") && !REPEAT_TRACK_PREFIX.equalsIgnoreCase("null");
		final int BUFFER_SIZE = 1000000;
		final int HISTOGRAM_BINS = 100;
		final int LENGTH_HISTOGRAM_MIN = 500;
		final int LENGTH_HISTOGRAM_MAX = (MIN_ALIGNMENT_LENGTH)<<1;
		final int LENGTH_HISTOGRAM_QUANTUM = (LENGTH_HISTOGRAM_MAX-LENGTH_HISTOGRAM_MIN+1)/HISTOGRAM_BINS;
		final double ERROR_HISTOGRAM_MAX = 0.4;
		final double ERROR_HISTOGRAM_QUANTUM = ERROR_HISTOGRAM_MAX/HISTOGRAM_BINS;
		final int MIN_INTERSECTION = IO.quantum;
		final double INTERSECTION_HISTOGRAM_QUANTUM = 1.0/HISTOGRAM_BINS;
		final int MIN_REPEAT_TRACK_INTERSECTION = IO.quantum<<1;
		final boolean CHECK_SYMMETRY_AGGRESSIVE = false;
		int i, p;
		int length, maxLength, readA, readB, previousReadA, previousReadB, last1, last2, diffs;
		int nAlignments, nAlignmentsKept1, nAlignmentsKept2, nAlignmentsWithIntersection, nAsymmetricDiffs;
		long surface, alignmentsSurface, intersectionA, intersectionB, alignmentsIntersection;
		double error, ratio;
		String str, key;
		BufferedReader br;
		BufferedWriter bw;
		ValuePair value;
		HashMap<String,ValuePair> map;
		final int[] errorHistogram = new int[HISTOGRAM_BINS];
		final int[] lengthHistogram = new int[HISTOGRAM_BINS];
		final int[] intersectionHistogram = new int[HISTOGRAM_BINS];
		
		IO.initialize();
		Reads.nReads=N_READS;
		Reads.maxReadLength=MAX_READ_LENGTH;
		Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.loadReadIDs(READ_IDS_FILE,N_READS);
		Reads.loadQualities(QUALITY_THRESHOLDS_FILE,QUALITIES_FILE);
		Math.set(errorHistogram,HISTOGRAM_BINS-1,0);
		Math.set(lengthHistogram,HISTOGRAM_BINS-1,0);
		Math.set(intersectionHistogram,HISTOGRAM_BINS-1,0);
		tmpArray5 = new int[Reads.maxReadLength];
		if (REPEAT_TRACK_EXISTS) {
			maxLength=loadRepeatTracks(REPEAT_TRACK_DIR,REPEAT_TRACK_PREFIX,REPEAT_TRACK_FORMAT,FIRST_READ,LAST_READ,MIN_REPEAT_TRACK_INTERSECTION);
			maxLength*=3;
			tmpArray1 = new int[maxLength]; 
			tmpArray2 = new int[maxLength];
			tmpArray3 = new int[maxLength];
			tmpArray4 = new int[maxLength];
		}
		else tmpArray3 = new int[100];  // Arbitrary
		
		System.err.println("FilterAlignments> Making sure that diff values are symmetric...");
		map = new HashMap<String,ValuePair>(N_ALIGNMENTS>>1);
		br = new BufferedReader(new FileReader(INPUT_FILE),BUFFER_SIZE);
		nAlignments=0; nAsymmetricDiffs=0;
		str=br.readLine(); str=br.readLine();  // Header
		str=br.readLine();
		while (str!=null) {
			nAlignments++;
			Alignments.readAlignmentFile(str);
			key=IO.getCanonicalForm(Alignments.readA,Alignments.readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
			value=map.get(key);
			if (value!=null) {
				value.count++;
				if (Alignments.diffs!=value.diffs) nAsymmetricDiffs++;
				if (Alignments.diffs<value.diffs) value.diffs=Alignments.diffs;
			}
			else map.put(key,new ValuePair(1,Alignments.diffs));
			str=br.readLine();
		}
		br.close();
		if (IO.CONSISTENCY_CHECKS) {
			if (nAlignments%2!=0) {
				System.err.println("FilterAlignments> ERROR: the input file contains an odd number of alignments: "+nAlignments);
				System.exit(1);
			}
			Set<Map.Entry<String,ValuePair>> set = map.entrySet();
			Iterator<Map.Entry<String,ValuePair>> iterator = set.iterator();
			Map.Entry<String,ValuePair> next;
			while (iterator.hasNext()) {
				next=iterator.next();
				if (next.getValue().count!=2) {
					System.err.println("FilterAlignments> ERROR: this alignment is present "+next.getValue().count+" times: "+next.getKey());
					System.exit(1);
				}
			}
		}
		nAlignments>>=1;
		System.err.println("FilterAlignments> The input file contains "+nAsymmetricDiffs+" distinct alignments with asymmetric diffs ("+IO.getPercent(nAsymmetricDiffs,nAlignments)+"% of all distinct alignments)");
		
		System.err.println("FilterAlignments> Filtering alignments...");
		br = new BufferedReader(new FileReader(INPUT_FILE),BUFFER_SIZE);
		bw = new BufferedWriter(new FileWriter(OUTPUT_FILE),BUFFER_SIZE);
		last1=-1; last2=-1; previousReadA=-1; previousReadB=-1;
		nAlignmentsKept1=0; nAlignmentsKept2=0; nAlignmentsWithIntersection=0;
		alignmentsSurface=0L; alignmentsIntersection=0L;
		str=br.readLine(); bw.write(str); bw.newLine();  // Header
		str=br.readLine(); bw.write(str); bw.newLine();
		str=br.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			// Filtering by alignment error rate and length
			surface=Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2;
			key=IO.getCanonicalForm(Alignments.readA,Alignments.readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
			diffs=map.get(key).diffs;
			error=((double)(diffs<<1))/surface;
			if (error>=ERROR_HISTOGRAM_MAX) errorHistogram[HISTOGRAM_BINS-1]++;
			else errorHistogram[(int)(error/ERROR_HISTOGRAM_QUANTUM)]++;
			if (error>MAX_ERROR) {
				str=br.readLine();
				continue;
			}
			length=Math.min(Alignments.endA-Alignments.startA,Alignments.endB-Alignments.startB)+1;
			if (length>=LENGTH_HISTOGRAM_MAX) lengthHistogram[HISTOGRAM_BINS-1]++;
			else if (length<=LENGTH_HISTOGRAM_MIN) lengthHistogram[0]++;
			else lengthHistogram[(length-LENGTH_HISTOGRAM_MIN)/LENGTH_HISTOGRAM_QUANTUM]++;
			if (length<MIN_ALIGNMENT_LENGTH) {
				str=br.readLine();
				continue;
			}
			nAlignmentsKept1++; alignmentsSurface+=surface;
			// Filtering by repeat track, if any.
			if (REPEAT_TRACK_EXISTS) {
				readA=Alignments.readA-1;
				if (readA!=previousReadA) {
					last1=getNonrepeatSurfaceA(readA,FIRST_READ,MIN_ALIGNMENT_LENGTH);
					previousReadA=readA;
				}
				readB=Alignments.readB-1;
				if (readB!=previousReadB) {
					last2=getNonrepeatSurfaceB(readB,FIRST_READ,MIN_ALIGNMENT_LENGTH);
					previousReadB=readB;
				}
				tmpArray3[0]=Alignments.startA; tmpArray3[1]=Alignments.endA;
				intersectionA=Intervals.intersectionLength(tmpArray3,1,repeatTrack[readA-FIRST_READ],lastRepeatTrack[readA-FIRST_READ],1,true);
				tmpArray3[0]=Alignments.startB; tmpArray3[1]=Alignments.endB;
				intersectionB=Intervals.intersectionLength(tmpArray3,1,repeatTrack[readB-FIRST_READ],lastRepeatTrack[readB-FIRST_READ],1,true);
				if (intersectionA>=MIN_INTERSECTION || intersectionB>=MIN_INTERSECTION) nAlignmentsWithIntersection++;
				alignmentsIntersection+=intersectionA+intersectionB;
				ratio=((double)(intersectionA+intersectionB))/surface;
				intersectionHistogram[Math.min((int)(ratio/INTERSECTION_HISTOGRAM_QUANTUM),HISTOGRAM_BINS-1)]++;
				if (CHECK_SYMMETRY_AGGRESSIVE && last1!=-1 && last2!=-1) {
					boolean result1 = hasNonrepeatSurface(tmpArray1,last1,tmpArray2,last2,FIRST_READ,LAST_READ,MIN_ALIGNMENT_LENGTH);
					Alignments.readAlignmentFile_flip();
					boolean result2 = hasNonrepeatSurface(tmpArray2,last2,tmpArray1,last1,FIRST_READ,LAST_READ,MIN_ALIGNMENT_LENGTH);
					if (result1!=result2) {
						System.err.println("FilterAlignments> ERROR: hasNonrepeatSurface() is asymmetric for "+str);
						System.exit(1);
					}
					Alignments.readAlignmentFile_flip();
				}
				if (last1==-1 || last2==-1 || !hasNonrepeatSurface(tmpArray1,last1,tmpArray2,last2,FIRST_READ,LAST_READ,MIN_ALIGNMENT_LENGTH)) {
					str=br.readLine();
					continue;
				}
			}			
			nAlignmentsKept2++;
			if (Alignments.diffs==diffs) bw.write(str); 
			else {
				diffsNumberInString(str,tmpArray3);
				for (i=0; i<tmpArray3[0]; i++) bw.write(str.charAt(i));
				bw.write(diffs+"");
				length=str.length();
				for (i=tmpArray3[1]+1; i<length; i++) bw.write(str.charAt(i));
			}
			bw.newLine();			
			str=br.readLine();
		}
		br.close(); bw.close();
		map.clear(); map=null;
		if (IO.CONSISTENCY_CHECKS) {
			if (nAlignmentsKept1%2!=0) {
				System.err.println("FilterAlignments> ERROR: nAlignmentsKept1="+nAlignmentsKept1+" is odd.");
				System.exit(1);
			}
			if (REPEAT_TRACK_EXISTS && nAlignmentsKept2%2!=0) {
				System.err.println("FilterAlignments> ERROR: nAlignmentsKept2="+nAlignmentsKept2+" is odd.");
				System.exit(1);
			}
		}
		nAlignmentsKept1>>=1; nAlignmentsKept2>>=1; nAlignmentsWithIntersection>>=1;
		System.out.println("FilterAlignments> Kept "+nAlignmentsKept2+" out of "+nAlignmentsKept1+" passing basic filters ("+IO.getPercent(nAlignmentsKept2,nAlignmentsKept1)+"%) and out of "+nAlignments+" total ("+IO.getPercent(nAlignmentsKept2,nAlignments)+"%)");
		if (REPEAT_TRACK_EXISTS) {
			System.out.println("FilterAlignments> Alignments that intersect the repeat track: "+nAlignmentsWithIntersection+" out of "+nAlignmentsKept1+" passing basic filters ("+IO.getPercent(nAlignmentsWithIntersection,nAlignmentsKept1)+"%)");
			System.out.println("FilterAlignments> Total intersection surface: "+alignmentsIntersection+" out of "+alignmentsSurface+" total alignment surface ("+IO.getPercent(alignmentsIntersection,alignmentsSurface)+"%)");
			System.out.println("FilterAlignments> Number of alignments with intersection:");
			for (i=0; i<HISTOGRAM_BINS; i++) System.out.println(IO.format(i*INTERSECTION_HISTOGRAM_QUANTUM)+","+intersectionHistogram[i]);
		}
		System.out.println("Histogram of avg. alignment error rates:");
		for (i=0; i<HISTOGRAM_BINS; i++) System.out.println(IO.format(i*ERROR_HISTOGRAM_QUANTUM)+","+(errorHistogram[i]>>1));
		System.out.println("Histogram of min alignment lengths:");
		for (i=0; i<HISTOGRAM_BINS; i++) System.out.println((LENGTH_HISTOGRAM_MIN+i*LENGTH_HISTOGRAM_QUANTUM)+","+(lengthHistogram[i]>>1));
		
		// Printing union of all repeat tracks
		if (REPEAT_TRACK_EXISTS) printRepeatTrack(REPEAT_TRACK_DIR+"/"+REPEAT_TRACK_PREFIX+"-union.txt",N_READS,REPEAT_TRACK_FORMAT);
	}	
	
	
	private static class ValuePair {
		public int count, diffs;
		
		public ValuePair(int c, int d) {
			this.count=c; this.diffs=d;
		}
	}
	
	
	/**
	 * @return interval [X..Y] such that str[X..Y] is the maximal string of numbers before 
	 * "diffs".
	 */
	private static final void diffsNumberInString(String str, int[] out) {
		final int ASCII_CHARACTER_ZERO = 48;
		final int ASCII_CHARACTER_NINE = 57;
		int c, i, p;
		
		p=str.indexOf(" diffs");
		for (i=p-1; i>=0; i--) {
			c=str.charAt(i);
			if (c<ASCII_CHARACTER_ZERO || c>ASCII_CHARACTER_NINE) break;
		}
		out[0]=i+1; out[1]=p-1;
	}

	
	




	
	// -------------------------------- REPEAT TRACKS ------------------------------------
	
	/**
	 * Interprets every file in $directory/prefix*$ as a distinct repeat track, and stores
	 * their union in global variables $repeatTrack,lastRepeatTrack$.
	 *
	 * Remark: the procedure computes statistics on every track compared to the union of 
	 * all previous tracks in lexicographic order.
	 *
	 * Remark: the procedure assumes that quality tracks have already been loaded.
	 *
	 * @param format format of the repeat track file: see $loadRepeatTrack()$ for details;
	 * @param minIntersection if two tracks intersect by at least this quantity in a read,
	 * an error occurs;
	 * @return the number of integers in a longest repeat track after the merge.
	 */
	private static final int loadRepeatTracks(String directory, String prefix, byte format, int firstRead, int lastRead, int minIntersection) throws IOException {
		final int INCREMENT = 10;  // Arbitrary, even.
		final int N_READS = lastRead-firstRead+1;
		final byte QUALITY_THRESHOLD = Reads.MAX_HIGH_QUALITY_SCORE;
		int i, j;
		int nTracks, track1, track2, length, maxLength, intersectionLength, last;
		long previousSurface, newSurface, totalHighQualitySurface, totalReadLength;
		File file;
		long[] trackSurface, highQualityIntersection;
		long[][] totalIntersection;
		int[][] lastRepeatTracks;
		int[][][] repeatTracks;
		String[] list, repeatTrackNames;
		
		// Loading all repeat tracks
		file = new File(directory);
		list=file.list();
		nTracks=0;
		for (i=0; i<list.length; i++) {
			if (list[i].startsWith(prefix)) nTracks++;
		}
		if (nTracks==0) {
			System.err.println("loadRepeatTracks> ERROR: no repeat track file found in <"+directory+">.");
			System.exit(1);
		}
		Arrays.sort(list);
		repeatTracks = new int[nTracks][N_READS][INCREMENT];
		lastRepeatTracks = new int[nTracks][N_READS];
		repeatTrackNames = new String[nTracks];
		trackSurface = new long[nTracks];
		Math.set(trackSurface,nTracks-1,0);
		j=-1; maxLength=0;
		for (i=0; i<list.length; i++) {
			if (!list[i].startsWith(prefix)) continue;
			j++;
			length=loadRepeatTrack(directory+"/"+list[i],format,firstRead,lastRead,repeatTracks[j],lastRepeatTracks[j]);
			maxLength=Math.max(maxLength,length);
			repeatTrackNames[j]=tmpString;
			trackSurface[j]=tmpSurface;
		}
		totalHighQualitySurface=Reads.getSurfaceAtQuality(QUALITY_THRESHOLD);
		
		// Checking if two tracks intersect
		if (IO.CONSISTENCY_CHECKS) {
			highQualityIntersection = new long[nTracks];
			Math.set(highQualityIntersection,nTracks-1,0);
			for (i=0; i<N_READS; i++) {
				j=Reads.getQualityAsTrack(firstRead+i,Reads.MAX_HIGH_QUALITY_SCORE,tmpArray5);
				for (track1=0; track1<nTracks; track1++) {
					intersectionLength=Intervals.intersectionLength(repeatTracks[track1][i],lastRepeatTracks[track1][i],tmpArray5,j,1,true);
					highQualityIntersection[track1]+=intersectionLength;
				}
			}
			totalIntersection = new long[nTracks][0];
			for (track1=0; track1<nTracks; track1++) {
				totalIntersection[track1] = new long[nTracks-track1-1];
				for (track2=track1+1; track2<nTracks; track2++) {
					totalIntersection[track1][track2-track1-1]=0;
					for (i=0; i<N_READS; i++) {
						intersectionLength=Intervals.intersectionLength(repeatTracks[track1][i],lastRepeatTracks[track1][i],repeatTracks[track2][i],lastRepeatTracks[track2][i],minIntersection,true);
						if (intersectionLength>0 && VERBOSE) {
							System.err.println("WARNING> Repeat tracks "+repeatTrackNames[track1]+" and "+repeatTrackNames[track2]+" intersect in read "+(firstRead+i)+" (zero-based) by >="+minIntersection+" bps:");
							System.err.print(repeatTrackNames[track1]+": ");
							for (j=0; j<lastRepeatTracks[track1][i]; j+=2) System.err.print("["+repeatTracks[track1][i][j]+".."+repeatTracks[track1][i][j+1]+"], ");
							System.err.println();
							System.err.print(repeatTrackNames[track2]+": ");
							for (j=0; j<lastRepeatTracks[track2][i]; j+=2) System.err.print("["+repeatTracks[track2][i][j]+".."+repeatTracks[track2][i][j+1]+"], ");
							System.err.println();
						}
						totalIntersection[track1][track2-track1-1]+=intersectionLength;
					}
				}
			}
			if (nTracks>1) {
				System.err.println("loadRepeatTracks> Track intersection statistics:");
				System.err.println("T1, T2, |T1 \\setminus T2|/|T1|, |T1 \\intersect T2|/|T1 \\union T2|, |T2 \\setminus T1|/|T2|");
				for (i=0; i<nTracks; i++) {
					for (j=i+1; j<nTracks; j++) System.err.println(repeatTrackNames[i]+", "+repeatTrackNames[j]+", "+IO.getPercent(trackSurface[i]-totalIntersection[i][j-i-1],trackSurface[i])+"%, "+IO.getPercent(totalIntersection[i][j-i-1],trackSurface[i]+trackSurface[j]-totalIntersection[i][j-i-1])+"%, "+IO.getPercent(trackSurface[j]-totalIntersection[i][j-i-1],trackSurface[j])+"%");
				}
				System.err.println();
				System.err.println("loadRepeatTracks> How each track intersects high-quality regions (Q):");
				System.err.println("Q, T, |Q \\setminus T|/|Q|, |Q \\intersect T|/|Q \\union T|, |T \\setminus Q|/|T|");
				for (i=0; i<nTracks; i++) System.err.println("Q, "+repeatTrackNames[i]+", "+IO.getPercent(totalHighQualitySurface-highQualityIntersection[i],totalHighQualitySurface)+"%, "+IO.getPercent(highQualityIntersection[i],totalHighQualitySurface+trackSurface[i]-highQualityIntersection[i])+"%, "+IO.getPercent(trackSurface[i]-highQualityIntersection[i],trackSurface[i])+"%");
				System.err.println();
			}
		}

		// Merging all tracks
		System.err.println("loadRepeatTracks> Intersection of a track (T) with all previous tracks (P):");
		System.err.println("P, T, |P \\setminus T|/|P|, |P \\intersect T|/|P \\union T|, |T \\setminus P|/|T|");
		if (tmpArray1==null || tmpArray1.length<maxLength<<1) tmpArray1 = new int[maxLength<<1];
		repeatTrack = new int[N_READS][INCREMENT];
		lastRepeatTrack = new int[N_READS];
		Math.set(lastRepeatTrack,N_READS-1,-1);
		totalIntersection = new long[nTracks][1];
		Math.set(totalIntersection,0);
		previousSurface=0;
		for (track1=0; track1<nTracks; track1++) {
			newSurface=0;
			for (i=0; i<N_READS; i++) {
				if (track1>0) {
					intersectionLength=Intervals.intersectionLength(repeatTrack[i],lastRepeatTrack[i],repeatTracks[track1][i],lastRepeatTracks[track1][i],minIntersection,true);
					totalIntersection[track1][0]+=intersectionLength;
				}
				last=Intervals.union(repeatTrack[i],lastRepeatTrack[i],repeatTracks[track1][i],lastRepeatTracks[track1][i],tmpArray1);
				lastRepeatTrack[i]=last;
				if (repeatTrack[i].length<last+1) repeatTrack[i] = new int[last+1];
				System.arraycopy(tmpArray1,0,repeatTrack[i],0,last+1);
				for (j=0; j<last; j+=2) newSurface+=repeatTrack[i][j+1]-repeatTrack[i][j]+1;
			}
			if (track1>0) {
				System.err.println("..., "+repeatTrackNames[track1]+", "+IO.getPercent(previousSurface-totalIntersection[track1][0],previousSurface)+"%, "+IO.getPercent(totalIntersection[track1][0],previousSurface+trackSurface[track1]-totalIntersection[track1][0])+"%, "+IO.getPercent(trackSurface[track1]-totalIntersection[track1][0],trackSurface[track1])+"%");
				System.err.println("..., "+repeatTrackNames[track1]+", "+(previousSurface-totalIntersection[track1][0])+", "+totalIntersection[track1][0]+", "+(trackSurface[track1]-totalIntersection[track1][0]));
			}
			else System.err.println("..., "+repeatTrackNames[track1]+", 0, 0, "+trackSurface[track1]);
			previousSurface=newSurface;
		}
		totalReadLength=Reads.getTotalReadLength();
		System.err.println("loadRepeatTracks> Total repeat surface: "+previousSurface+" ("+IO.getPercent(previousSurface,totalReadLength)+"% of total read surface)");
		if (IO.CONSISTENCY_CHECKS) {
			highQualityIntersection[0]=0;
			for (i=0; i<N_READS; i++) {
				j=Reads.getQualityAsTrack(firstRead+i,Reads.MAX_HIGH_QUALITY_SCORE,tmpArray5);
				highQualityIntersection[0]+=Intervals.intersectionLength(repeatTrack[i],lastRepeatTrack[i],tmpArray5,j,1,true);
			}
			System.err.println("loadRepeatTracks> How the union of all tracks intersects high-quality regions (Q):");
			System.err.println("Q, U, |Q \\setminus U|/|Q|, |Q \\intersect U|/|Q \\union U|, |U \\setminus Q|/|U|");
			System.err.println("Q, U, "+IO.getPercent(totalHighQualitySurface-highQualityIntersection[0],totalHighQualitySurface)+"%, "+IO.getPercent(highQualityIntersection[0],totalHighQualitySurface+previousSurface-highQualityIntersection[0])+"%, "+IO.getPercent(previousSurface-highQualityIntersection[0],previousSurface)+"%");
		}		
		System.err.println("loadRepeatTracks> Total high-quality surface: "+IO.getPercent(totalHighQualitySurface,totalReadLength)+"% of total read surface");
		System.err.println();
		
		// Deallocating
		for (i=0; i<nTracks; i++) {
			for (j=0; j<N_READS; j++) repeatTracks[i][j]=null;
			lastRepeatTracks[i]=null;
		}
		
		return maxLength;
	}
	
	
	/**
	 * Loads in matrices $repeatTrack,lastRepeatTrack$ the set of repeat tracks of all 
	 * reads with ID in $[firstRead..lastRead]$. Loads in global variable $tmpString$ the 
	 * string ID of the track, and in $tmpSurface$ the total surface of the track.
	 *
	 * @param format 0=$inputFile$ is in DBdump format; 1=every row encodes a list of 
	 * intervals, with format $read start_1 end_1 ... start_n end_n$ (all numbers are 
	 * zero-based); 2=every row encodes a single interval, with format $read start end$
	 * (all numbers are zero-based);
	 * @param firstRead,lastRead zero-based;
	 * @param repeatTrack,lastRepeatTrack assumed to already contain a large enough number
	 * of rows;
	 * @return the number of integers in a longest repeat track.
	 */
	public static final int loadRepeatTrack(String inputFile, byte format, int firstRead, int lastRead, int[][] repeatTrack, int[] lastRepeatTrack) throws IOException {
		final int INCREMENT = 10;  // Arbitrary, even.
		final int N_READS = lastRead-firstRead+1;
		char c;
		int i, j, p;
		int read, nElements, last, out;
		String str;
		BufferedReader br;
		String[] tokens;
		
		// Loading the file
		tmpString=null;
		Math.set(lastRepeatTrack,N_READS-1,-1);
		br = new BufferedReader(new FileReader(inputFile));
		str=br.readLine(); out=0;
		if (format==0) {
			while (str!=null) {
				c=str.charAt(0);
				if (c=='@') {
					tmpString=str.substring(str.lastIndexOf(" ")+1).trim();
					str=br.readLine();
					continue;
				}
				if (c!='R') {
					str=br.readLine();
					continue;
				}
				p=str.indexOf(" ");			
				read=Integer.parseInt(str.substring(p+1))-1-firstRead;
				do {
					str=br.readLine();
					c=str.charAt(0);
				} while (c!='T');
				tokens=str.split(" ");
				nElements=tokens.length-2;
				if (repeatTrack[read].length<nElements) repeatTrack[read] = new int[nElements];
				for (i=0; i<nElements; i++) repeatTrack[read][i]=Integer.parseInt(tokens[2+i]);
				lastRepeatTrack[read]=nElements-1;
				out=Math.max(out,nElements);
				str=br.readLine();
			}
		}
		else if (format==1) {
			while (str!=null) {
				tokens=str.split(" ");
				read=Integer.parseInt(tokens[0])-firstRead;
				nElements=tokens.length-1;
				if (repeatTrack[read].length<nElements) repeatTrack[read] = new int[nElements];
				for (i=0; i<nElements; i++) repeatTrack[read][i]=Integer.parseInt(tokens[1+i]);
				lastRepeatTrack[read]=nElements-1;
				str=br.readLine();
			}
		}
		else if (format==2) {
			while (str!=null) {
				tokens=str.split(" ");
				read=Integer.parseInt(tokens[0])-firstRead;
				if (lastRepeatTrack[read]+2>=repeatTrack[read].length) {
					int[] newTrack = new int[Math.max(repeatTrack[read].length<<1,lastRepeatTrack[read]+3)];
					System.arraycopy(repeatTrack[read],0,newTrack,0,repeatTrack[read].length);
					repeatTrack[read]=newTrack;
				}
				repeatTrack[read][++lastRepeatTrack[read]]=Integer.parseInt(tokens[1]);
				repeatTrack[read][++lastRepeatTrack[read]]=Integer.parseInt(tokens[2]);
				str=br.readLine();
			}
		}
		br.close();
		
		// Making sure intervals do not overlap, and computing the total surface.
		tmpSurface=0;
		for (i=firstRead; i<=lastRead; i++) {
			last=lastRepeatTrack[i-firstRead];
			for (j=2; j<last; j+=2) repeatTrack[i-firstRead][j]=Math.max(repeatTrack[i-firstRead][j],repeatTrack[i-firstRead][j-1]+1);
			for (j=0; j<last; j+=2) tmpSurface+=repeatTrack[i-firstRead][j+1]-repeatTrack[i-firstRead][j]+1;
		}
		
		return out;
	}	
	
	
	/**
	 * Loads in $tmpArray1[0..X]$ the sorted set of intervals, each of length 
	 * $>=minIntervalLength$, that do not belong to the repeat track of readA, where X is 
	 * returned in output. Every interval is a triplet $(-1,start,end)$.
	 *
	 * Remark: the procedure assumes that read lengths have already been loaded.
	 *
	 * Remark: the procedure uses global arrays $tmpArray{1,3}$, which are assumed to be 
	 * large enough.
	 *
	 * @param readA zero-based.
	 */
	private static final int getNonrepeatSurfaceA(int readA, int firstRead, int minIntervalLength) {
		final int readAPrime = readA-firstRead;
		final int lastTrackA = lastRepeatTrack[readAPrime];
		int i;
		
		if (lastTrackA==-1) {
			tmpArray1[1]=0; tmpArray1[2]=Reads.getReadLength(readA)-1; 
			return 2;
		}
		else {
			tmpArray3[0]=0; tmpArray3[1]=Reads.getReadLength(readA)-1;
			return Intervals.difference(tmpArray3,1,repeatTrack[readAPrime],lastTrackA,tmpArray1,0,-1,minIntervalLength);
		}
	}
	
	
	/**
	 * Like $getNonrepeatSurfaceA()$, but stores the result in $tmpArray2$.
	 */
	private static final int getNonrepeatSurfaceB(int readB, int firstRead, int minIntervalLength) {
		final int readBPrime = readB-firstRead;
		final int lastTrackB = lastRepeatTrack[readBPrime];
		int i;
		
		if (lastTrackB==-1) {
			tmpArray2[1]=0; tmpArray2[2]=Reads.getReadLength(readB)-1; 
			return 2;
		}
		else {
			tmpArray3[0]=0; tmpArray3[1]=Reads.getReadLength(readB)-1;
			return Intervals.difference(tmpArray3,1,repeatTrack[readBPrime],lastTrackB,tmpArray2,0,-1,minIntervalLength);
		}
	}
	
	
	/**
	 * Remark: the procedure uses global arrays $tmpArray{3,4}$, which are assumed to be 
	 * large enough.
	 *
	 * @param intervalsX,lastIntervalX output of $getNonrepeatSurfaceX()$;
	 * @param firstRead,lastRead zero-based;
	 * @return TRUE iff the alignment just loaded by $Alignments.readAlignmentFile()$ maps
	 * a non-repeat substring of readA of length $>=minIntervalLength$, to a non-repeat
	 * substring of readB of length $>=minIntervalLength$ (or vice versa).
	 */
	private static final boolean hasNonrepeatSurface(int[] intervalsA, int lastIntervalA, int[] intervalsB, int lastIntervalB, int firstRead, int lastRead, int minIntervalLength) {
		final int readA = Alignments.readA-1;
		final int readAPrime = readA-firstRead;
		final int lastTrackA = lastRepeatTrack[readAPrime];
		final int readB = Alignments.readB-1;
		final int readBPrime = readB-firstRead;
		final int lastTrackB = lastRepeatTrack[readBPrime];
		int i;
		int from, to, last4;
		
		for (i=0; i<lastIntervalA; i+=3) {
			from=Math.max(intervalsA[i+1],Alignments.startA);
			to=Math.min(intervalsA[i+2],Alignments.endA);
			if (to-from+1<minIntervalLength) continue;
			Intervals.project(intervalsA[i+1],intervalsA[i+2],Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.orientation,tmpArray3,0);
			last4=Intervals.difference(tmpArray3,1,repeatTrack[readBPrime],lastTrackB,tmpArray4,0,-1,minIntervalLength);
			if (last4!=-1) return true;
		}
		for (i=0; i<lastIntervalB; i+=3) {
			from=Math.max(intervalsB[i+1],Alignments.startB);
			to=Math.min(intervalsB[i+2],Alignments.endB);
			if (to-from+1<minIntervalLength) continue;
			Intervals.project(intervalsB[i+1],intervalsB[i+2],Alignments.startB,Alignments.endB,Alignments.startA,Alignments.endA,Alignments.orientation,tmpArray3,0);
			last4=Intervals.difference(tmpArray3,1,repeatTrack[readAPrime],lastTrackA,tmpArray4,0,-1,minIntervalLength);
			if (last4!=-1) return true;
		}
		return false;
	}
	
	
	/**
	 * @param format see $loadRepeatTrack()$.
	 */
	private static final void printRepeatTrack(String file, int nReads, byte format) throws IOException {
		int i, j;
		int last;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(file));
		if (format==0) {
			for (i=0; i<nReads; i++) {
				bw.write("R "+(Reads.firstRead+i+1)+"\n");  // Read IDs are 1-based
				last=lastRepeatTrack[i];
				bw.write("T0 "+((last+1)>>1));
				for (j=0; j<=last; j++) {
					bw.write(" "+repeatTrack[i][j]);  // Coordinates are 0-based.
				}
				bw.write("\n");
			}
		}
		else if (format==1) {
			for (i=0; i<nReads; i++) {
				bw.write((Reads.firstRead+i)+" ");
				last=lastRepeatTrack[i];
				for (j=0; j<=last; j++) bw.write(repeatTrack[i][j]+" ");
				bw.write("\n");
			}
		}
		else if (format==2) {
			for (i=0; i<nReads; i++) {
				last=lastRepeatTrack[i];
				for (j=0; j<last; j+=2) bw.write((Reads.firstRead+i)+" "+repeatTrack[i][j]+" "+repeatTrack[i][j+1]+"\n");
			}
		}
		bw.close();
	}
	
}


