package de.mpi_cbg.revant.fragments;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.PeriodicSubstrings;
import de.mpi_cbg.revant.factorize.Intervals;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;


/** 
 * Prints a report with the following lines: 
 * 
 * id,ar,nr,nf,ni,p,d,l
 *
 * id: ID of a basin;
 * ar: 1 iff enough fragments fully align to the reference (so the reference can be used
 *     as a guide sequence to compute the consensus);
 * nr: fraction of fragments that are not covered with alignments to the reference;
 * nf: fraction of NR that are not covered with alignments to other fragments;
 * ni: fraction of NF that are not covered with alignments to other intervals of the basin
 *     (these might be longer than fragments), considering the input alignments used in 
 *     the whole repeat analysis pipeline;
 * p: fraction of fragments that are likely permutations of the reference;
 * d: fraction of fragments that are likely deletions of the reference;
 * l: fraction of fragments that are too long WRT the reference.
 *  
 * The procedure writes also the artificial connection file that will be used by 
 * $FragmentsStep2$, discarding and trimming some fragments: see 
 * $writeArtificialConnectionFile()$ for details.
 *
 * Remark: NI is useful, since fragments are arbitrary substrings of the original  
 * intervals produced by repeat analysis (e.g. they might be just very-high-quality
 * substrings of them), and this could create issues with the min alignment length used
 * for fragment-fragment alignments. Some intervals might not even produce any fragment
 * at all. Moreover, the exact settings of the aligner used to build the alignments used 
 * for repeat analysis might be unknown. And even if they are known, the heuristics of the
 * aligner might produce different results when invoked on fragments rather than on full 
 * reads, since fragments can be significantly shorter than reads.
 */
public class FragmentsStep1 {
	/**
	 * Global constants
	 */
	private static final int IDENTITY_THRESHOLD = IO.quantum;
	private static final int SURFACE_THRESHOLD = IO.quantum<<1;  // Arbitrary	
	public static final int TOO_FEW_DIFFS = 1;  // Arbitrary
	public static final int N_FRAGMENTS_REFERENCE_ALIGNMENTS_THRESHOLD = 2;  // Arbitrary
	
	/**
	 * All read-read alignments, sorted by $readA,startA$.
	 */
	private static FFAlignment[] reads2reads;
	private static int N_RR_ALIGNMENTS;
	
	/**
	 * Number of fragments in the current module
	 */
	private static int N_FRAGMENTS;
	
	/**
	 * $read,start,end$ coordinates of each fragment in the current module 
	 * (rows=coordinate types, columns=fragments).
	 */
	private static int[][] fragmentCoordinates;
	
	/**
	 * Number of fragment-fragment alignments in the current module
	 */
	private static int N_FF_ALIGNMENTS;
	
	/**
	 * All fragment-fragment alignments in the current module, sorted by $readA,startA$.
	 */
	private static FFAlignment[] fragments2fragments;
	
	/**
	 * First element of $fragments2fragments$ whose readA is equal to a given fragment.
	 */
	private static int[] fragmentFirst;
	
	/**
	 * $read,start,end$ coordinates of each interval in the basin of the current module 
	 * (rows=coordinate types, columns=intervals).
	 */
	private static int[][] intervalCoordinates;
	
	/**
	 * Number of intervals in the basin of the current module
	 */
	private static int N_INTERVALS;
	
	/**
	 * Temporary space
	 */
	private static FFAlignment tmpAlignment = new FFAlignment();
	private static int[] tmpArray = new int[100];  // Arbitrary (even length).
	private static int[] tmpArray1 = new int[100];  // Arbitrary (even length).
	private static int[] tmpArray2 = new int[100];  // Arbitrary (even length).
	private static int[] tmpArray3 = new int[100];  // Arbitrary (even length).
	private static Pair[] tmpPairs = new Pair[100];  // Arbitrary
	
	
	/**
	 * @param args 
	 * 0: Step1 directory;
	 * 1: min length of an alignment between fragments and reference or other fragments;
	 * 2: (0/1) consider alignments of fragments to reference with very small diffs as 
	 * errors.
	 */
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		Alignments.minAlignmentLength=Integer.parseInt(args[1]);
		final boolean SMALL_DIFFS_ARE_ERRORS = Integer.parseInt(args[2])==1;
		
		final String RR_ALIGNMENTS_FILE = STEP1_DIR+"/../"+IO.ALIGNMENTS_FILE_STEP1;
		final String STEP4_DIR = STEP1_DIR+"/"+IO.TAGS_DIR+"/"+IO.STEP4_DIR;
		final String STEP5_DIR = STEP4_DIR+"/"+IO.STEP5_DIR;
		final String FRAGMENTS_COORDINATES_DIR = STEP5_DIR+"/"+IO.FRAGMENTS_LABEL;
		final String FRAGMENTS_ALIGNMENTS_DIR = STEP5_DIR+"/fragments-strings-alignments";
		final String BASINS_DIR = STEP4_DIR;
		
		final int MIN_INTERSECTION_LENGTH = Alignments.minAlignmentLength>>1;  // Arbitrary
		final int LENGTH_THRESHOLD = IO.quantum<<2;  // The same as in $IntervalGraphStep3.getBasins_ensureKernelLengths()$.
		final int GROWTH_RATE = 100;  // Arbitrary
		boolean currentReference, isCyclicAndNonperiodic;
		int i, j;
		int currentRead, currentSurface, length, last;
		int ffUncoveredSurface, rrUncoveredSurface, lastTmpWindow, nAlignments;
		int REFERENCE_LENGTH, BASIN_TYPE, PATH_ID, PERIOD, PERIOD_PRIME;
		String str, id, FRAGMENTS_LENGTHS_FILE, FF_ALIGNMENTS_FILE, FR_ALIGNMENTS_FILE, FRAGMENTS_COORDINATES_FILE, BASINS_FILE, CONNECTION_FILE;
		BufferedReader br;
		BufferedWriter statsFile;
		File file;
		boolean[] fragmentFound = new boolean[100];  // Arbitrary, will be resized.
		int[] stats = new int[6];  
		// 0=not fully aligning to reference; 
		// 1=not fully aligning to reference, and not fully aligning to fragments;
		// 2=like [1], and moreover not fully aligning to fragments using read-read
		// alignments.
		// 3=likely permutations of reference;
		// 4=likely deletions of reference;
		// 5=too long WRT reference.
		int[] tmpIO = new int[2];
		int[] frUncoveredHistogram = new int[100];  // Bins are multiples of $IO.quantum$
		int[] ffUncoveredHistogram = new int[100];  // Bins are multiples of $IO.quantum$
		int[] rrUncoveredHistogram = new int[100];  // Bins are multiples of $IO.quantum$
		String[] list;
		PermutationWindow[] tmpWindows = new PermutationWindow[100];  // Arbitrary
		int[][] newFragmentCoordinates = new int[100][2];  // Arbitrary, will be resized.
		
		System.err.print("Loading read-read alignments... ");
		br = new BufferedReader(new FileReader(RR_ALIGNMENTS_FILE));
		i=0; str=br.readLine();
		while (str!=null) {
			i++;
			str=br.readLine();
		}
		br.close();
		N_RR_ALIGNMENTS=i;
		IO.initialize();
		for (i=0; i<tmpPairs.length; i++) tmpPairs[i] = new Pair();
		loadReads2reads(RR_ALIGNMENTS_FILE);
		System.err.println("done");
		Math.set(frUncoveredHistogram,frUncoveredHistogram.length-1,0);
		Math.set(ffUncoveredHistogram,ffUncoveredHistogram.length-1,0);
		Math.set(rrUncoveredHistogram,rrUncoveredHistogram.length-1,0);
		for (i=0; i<tmpWindows.length; i++) tmpWindows[i] = new PermutationWindow();
		statsFile = new BufferedWriter(new FileWriter(FRAGMENTS_ALIGNMENTS_DIR+"/step1-stats.txt"));
		file = new File(FRAGMENTS_ALIGNMENTS_DIR);
		list=file.list();
		for (i=0; i<list.length; i++) {
			if (list[i].length()<IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH || !list[i].substring(0,IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH).equalsIgnoreCase(IO.MODULE_ALIGNMENTS_DIR_PREFIX)) continue;
			
			// Initializing data structures for the current repeat module
			id=list[i].substring(IO.MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH);
			System.err.print("Loading data structures for repeat "+id+"... ");
			BASINS_FILE=BASINS_DIR+"/"+IO.BASIN_PREFIX+id+".txt";
			br = new BufferedReader(new FileReader(BASINS_FILE));
			str=br.readLine(); br.close();
			IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
			REFERENCE_LENGTH=tmpArray[1]; PATH_ID=tmpArray[2]; BASIN_TYPE=tmpArray[3]; N_INTERVALS=tmpArray[5];
			PERIOD=tmpArray[4]>=PeriodicSubstrings.MIN_PERIOD_LONG?tmpArray[4]:-1;
			PERIOD_PRIME=tmpArray[4];
			isCyclicAndNonperiodic=(PATH_ID==-1)&&(BASIN_TYPE!=Constants.INTERVAL_PERIODIC);
			loadIntervalCoordinates(BASINS_FILE);
			FRAGMENTS_LENGTHS_FILE=FRAGMENTS_ALIGNMENTS_DIR+"/"+list[i]+"/reads-lengths.txt";
			br = new BufferedReader(new FileReader(FRAGMENTS_LENGTHS_FILE));
			N_FRAGMENTS=0; str=br.readLine();
			while (str!=null) {
				N_FRAGMENTS++;
				str=br.readLine();
			}
			br.close();
			Reads.nReads=N_FRAGMENTS;
			Reads.loadReadLengths(FRAGMENTS_LENGTHS_FILE);
			Reads.loadReadIDs(0,N_FRAGMENTS-1);
			if (N_FRAGMENTS>newFragmentCoordinates.length) newFragmentCoordinates = new int[N_FRAGMENTS][2];
			Math.set(stats,stats.length-1,0);
			if (BASIN_TYPE!=Constants.INTERVAL_PERIODIC && PATH_ID!=-1) {
				for (j=0; j<N_FRAGMENTS; j++) {
					length=Reads.getReadLength(j);
					if (length>REFERENCE_LENGTH+LENGTH_THRESHOLD) {
						System.err.println("WARNING: Fragment "+j+" has length "+length+" > referenceLength="+REFERENCE_LENGTH+" (this is not necessarily an error)");
						newFragmentCoordinates[j][0]=-1; newFragmentCoordinates[j][1]=-1;
						stats[5]++;
					}
					else { newFragmentCoordinates[j][0]=0; newFragmentCoordinates[j][1]=0; }
				}
			}
			FF_ALIGNMENTS_FILE=FRAGMENTS_ALIGNMENTS_DIR+"/"+list[i]+"/LAshow.txt";
			br = new BufferedReader(new FileReader(FF_ALIGNMENTS_FILE));
			str=br.readLine(); str=br.readLine();  // Skipping header
			N_FF_ALIGNMENTS=0; str=br.readLine();
			while (str!=null) {
				N_FF_ALIGNMENTS++;
				str=br.readLine();
			}
			br.close();
			loadFragments2Fragments(FF_ALIGNMENTS_FILE);
			FR_ALIGNMENTS_FILE=FRAGMENTS_ALIGNMENTS_DIR+"/"+list[i]+"/"+IO.FRAGMENTS_REFERENCE_FILE+".txt";
			FRAGMENTS_COORDINATES_FILE=FRAGMENTS_COORDINATES_DIR+"/fragments-"+id+".txt";
			loadFragmentCoordinates(FRAGMENTS_COORDINATES_FILE);
			if (N_FRAGMENTS>fragmentFound.length) fragmentFound = new boolean[N_FRAGMENTS];
			Math.set(fragmentFound,N_FRAGMENTS-1,false);
			CONNECTION_FILE=FRAGMENTS_ALIGNMENTS_DIR+"/"+list[i]+"/connection.txt";
			System.err.println("done");
			
			// Checking fragments that have some alignment to the reference
			br = new BufferedReader(new FileReader(FR_ALIGNMENTS_FILE),IO.BUFFER_SIZE);
			currentRead=-1; currentSurface=-1; currentReference=false; lastTmpWindow=-1;
			str=br.readLine(); str=br.readLine();  // Skipping header
			str=br.readLine(); nAlignments=0;
			while (str!=null) {
				nAlignments++;
				Alignments.readAlignmentFile(str);
				if (Alignments.readA-1!=currentRead) {
					if (currentRead!=-1) {
						fragmentFound[currentRead]=true;
						if (newFragmentCoordinates[currentRead][0]!=-1) {
							checkSurfaces(currentRead,currentReference,isCyclicAndNonperiodic,BASIN_TYPE,PERIOD,REFERENCE_LENGTH,MIN_INTERSECTION_LENGTH,tmpWindows,lastTmpWindow,newFragmentCoordinates,stats,frUncoveredHistogram,ffUncoveredHistogram,rrUncoveredHistogram,list[i]);
							if ((BASIN_TYPE!=Constants.INTERVAL_PERIODIC || PERIOD>=PeriodicSubstrings.MIN_PERIOD_LONG) && newFragmentCoordinates[currentRead][0]!=-1 && lastTmpWindow>0) {
								if (isPermutation(tmpWindows,lastTmpWindow,PERIOD,tmpIO)) {
									newFragmentCoordinates[currentRead][0]=-1;
									newFragmentCoordinates[currentRead][1]=-1;
									stats[3]++;
								}
								else {
									if (isDeletion(tmpWindows,lastTmpWindow,PERIOD)) {
										newFragmentCoordinates[currentRead][0]=-1;
										newFragmentCoordinates[currentRead][1]=-1;
										stats[4]++;
									}
								}
							}
						}
					}
					currentRead=Alignments.readA-1; 
					currentReference=Alignments.diffs<=TOO_FEW_DIFFS;
					if (currentReference) {
						if (SMALL_DIFFS_ARE_ERRORS) {
							System.err.println("ERROR: Fragment "+(Alignments.readA-1)+" is too similar to the reference ("+Alignments.diffs+" diffs):");
							System.err.println("   "+str);
							System.exit(1);
						}
						else {
							System.err.println("WARNING: Fragment "+(Alignments.readA-1)+" is too similar to the reference ("+Alignments.diffs+" diffs):");
							System.err.println("   "+str);
						}
					}
					lastTmpWindow=0;
					tmpWindows[0].set(Alignments.readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.diffs);
					str=br.readLine();
					continue;
				}
				if (Alignments.diffs<=TOO_FEW_DIFFS) {
					if (SMALL_DIFFS_ARE_ERRORS) {
						System.err.println("ERROR: Fragment "+(Alignments.readA-1)+" is too similar to the reference ("+Alignments.diffs+" diffs):");
						System.err.println("   "+str);
						System.exit(1);
					}
					else {
						System.err.println("WARNING: Fragment "+(Alignments.readA-1)+" is too similar to the reference ("+Alignments.diffs+" diffs):");
						System.err.println("   "+str);
					}
					currentReference=true;
				}
				lastTmpWindow++;
				if (lastTmpWindow==tmpWindows.length) {
					PermutationWindow[] newArray = new PermutationWindow[tmpWindows.length+GROWTH_RATE];
					System.arraycopy(tmpWindows,0,newArray,0,tmpWindows.length);
					for (j=tmpWindows.length; j<newArray.length; j++) newArray[j] = new PermutationWindow();
					tmpWindows=newArray;
				}
				tmpWindows[lastTmpWindow].set(Alignments.readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB,Alignments.diffs);
				str=br.readLine();
			}
			br.close();
			if (currentRead!=-1) {
				fragmentFound[currentRead]=true;
				if (newFragmentCoordinates[currentRead][0]!=-1) {
					checkSurfaces(currentRead,currentReference,isCyclicAndNonperiodic,BASIN_TYPE,PERIOD,REFERENCE_LENGTH,MIN_INTERSECTION_LENGTH,tmpWindows,lastTmpWindow,newFragmentCoordinates,stats,frUncoveredHistogram,ffUncoveredHistogram,rrUncoveredHistogram,list[i]);
					if ((BASIN_TYPE!=Constants.INTERVAL_PERIODIC || PERIOD>=PeriodicSubstrings.MIN_PERIOD_LONG) && newFragmentCoordinates[currentRead][0]!=-1 && lastTmpWindow>0) {
						if (isPermutation(tmpWindows,lastTmpWindow,PERIOD,tmpIO)) {
							newFragmentCoordinates[currentRead][0]=-1;
							newFragmentCoordinates[currentRead][1]=-1;
							stats[3]++;
						}
						else {
							if (isDeletion(tmpWindows,lastTmpWindow,PERIOD)) {
								newFragmentCoordinates[currentRead][0]=-1;
								newFragmentCoordinates[currentRead][1]=-1;
								stats[4]++;
							}
						}
					}
				}
			}
			if ( BASIN_TYPE==Constants.INTERVAL_PERIODIC && PERIOD_PRIME>0 && (PERIOD_PRIME)<PeriodicSubstrings.MIN_PERIOD_LONG && 
			     nAlignments<N_FRAGMENTS*N_FRAGMENTS_REFERENCE_ALIGNMENTS_THRESHOLD
			   )  System.err.println("WARNING: Too few fragment-reference alignments in short-period repeat "+id+": k-mer capping by the aligner?");
			
			// Checking fragments with no alignment to the reference
			for (currentRead=0; currentRead<N_FRAGMENTS; currentRead++) {
				if (fragmentFound[currentRead]) continue;
				stats[0]++;
				currentSurface=Reads.getReadLength(currentRead);
				System.err.println("Fragment "+currentRead+" does not fully align to the reference for "+currentSurface+" bp");
				currentSurface/=IO.quantum;
				currentSurface=Math.min(currentSurface,frUncoveredHistogram.length-1);
				frUncoveredHistogram[currentSurface]++;
				ffUncoveredSurface=uncoveredSurface_fragment2fragment(currentRead,null,-1,newFragmentCoordinates[currentRead]);
				if (ffUncoveredSurface>=SURFACE_THRESHOLD) {
					stats[1]++;
					System.err.println("   Fragment "+currentRead+" does not fully align to the reference nor to other fragments for "+ffUncoveredSurface+" bp");
					ffUncoveredSurface/=IO.quantum;
					ffUncoveredSurface=Math.min(ffUncoveredSurface,ffUncoveredHistogram.length-1);
					ffUncoveredHistogram[ffUncoveredSurface]++;
					rrUncoveredSurface=uncoveredSurface_read2read(currentRead,MIN_INTERSECTION_LENGTH);
					if (rrUncoveredSurface>=SURFACE_THRESHOLD && !(BASIN_TYPE==Constants.INTERVAL_PERIODIC && PERIOD==-1 && Reads.getReadLength(currentRead)>REFERENCE_LENGTH)) {
						stats[2]++;
						System.err.println("      Fragment "+currentRead+" does not fully align to the reference, nor to other fragments, nor to basin intervals using read-read alignments, for "+rrUncoveredSurface+" bp ("+list[i]+").");
						rrUncoveredSurface/=IO.quantum;
						rrUncoveredSurface=Math.min(rrUncoveredSurface,rrUncoveredHistogram.length-1);
						rrUncoveredHistogram[rrUncoveredSurface]++;
					}
				}
			}
			statsFile.write(id+","+((stats[0]+stats[3]+stats[4]+stats[5]>(1.0-FragmentsStep2.DOMINATOR_THRESHOLD)*N_FRAGMENTS)?"0":"1")+","+IO.format(100*((double)stats[0])/N_FRAGMENTS)+","+(stats[0]!=0?IO.format(100*((double)stats[1])/stats[0]):"0")+","+(stats[1]!=0?IO.format(100*((double)stats[2])/stats[1]):"0")+","+(stats[3]!=0?IO.format(100*((double)stats[3])/N_FRAGMENTS):"0")+","+(stats[4]!=0?IO.format(100*((double)stats[4])/N_FRAGMENTS):"0")+","+(stats[5]!=0?IO.format(100*((double)stats[5])/N_FRAGMENTS):"0")+"\n");
			writeArtificialConnectionFile(FF_ALIGNMENTS_FILE,CONNECTION_FILE,newFragmentCoordinates);			
		}
		statsFile.close();
		
		// Printing global stats
		System.out.println("Histogram of unaligned FR surface over all modules:");
		last=frUncoveredHistogram.length;
		for (i=frUncoveredHistogram.length-1; i>=0; i--) {
			if (frUncoveredHistogram[i]!=0) {
				last=i;
				break;
			}
		}
		if (last<frUncoveredHistogram.length) {
			for (i=0; i<=last; i++) System.out.println(i*IO.quantum+","+frUncoveredHistogram[i]);
		}
		else System.out.println("<empty>");
		System.out.println("Histogram of unaligned FF surface over all modules:");
		last=ffUncoveredHistogram.length;
		for (i=ffUncoveredHistogram.length-1; i>=0; i--) {
			if (ffUncoveredHistogram[i]!=0) {
				last=i;
				break;
			}
		}
		if (last<ffUncoveredHistogram.length) {
			for (i=0; i<=last; i++) System.out.println(i*IO.quantum+","+ffUncoveredHistogram[i]);
		}
		else System.out.println("<empty>");
		System.out.println("Histogram of unaligned RR surface over all modules:");
		last=rrUncoveredHistogram.length;
		for (i=rrUncoveredHistogram.length-1; i>=0; i--) {
			if (rrUncoveredHistogram[i]!=0) {
				last=i;
				break;
			}
		}
		if (last<rrUncoveredHistogram.length) {
			for (i=0; i<=last; i++) System.out.println(i*IO.quantum+","+rrUncoveredHistogram[i]);
		}
		else System.out.println("<empty>");
	}
	
	
	/**
	 * Remark: the procedure uses global variables $tmpArray{1,2,3}$ as temporary space.
	 *
	 * @param isCyclicAndNonperiodic if the kernel is cyclic and nonperiodic, the 
	 * reference chosen by $getRepresentativeInterval_nonperiodic()$ is just a longest 
	 * one, and this might not be a node of the kernel graph; the latter nodes, being 
	 * maximal by containment, might contain arbitrarily long spurious substrings, which 
	 * might not align to the newly chosen reference;
	 * @param basinType,period,referenceLength it is possible that a short-period kernel 
	 * contains fragments from long-period or non-periodic intervals, that such fragments 
	 * are longer than the reference, and that they have regions that do not align to the 
	 * reference nor to any other interval in the kernel: see e.g. $IntervalGraphStep3.
	 * handleShortPeriodGraph()$.
	 * @param windows not assumed to be sorted in any way; the procedure might change the
	 * order.
	 */
	private static final void checkSurfaces(int currentRead, boolean currentReference, boolean isCyclicAndNonperiodic, int basinType, int period, int referenceLength, int minIntersectionLength, PermutationWindow[] windows, int lastWindow, int[][] newFragmentCoordinates, int[] stats, int[] frUncoveredHistogram, int[] ffUncoveredHistogram, int[] rrUncoveredHistogram, String fileName) {
		int i;
		int last, currentSurface, currentStart, currentEnd;
		int ffUncoveredSurface, rrUncoveredSurface;
		
		// Computing the total surface covered by alignments to the reference, in all
		// orientations.
		PermutationWindow.order=PermutationWindow.ORDER_STARTA;
		if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
		currentStart=-1; currentEnd=-1; last=-1;
		for (i=0; i<=lastWindow; i++) {
			if (currentStart==-1) { 
				currentStart=windows[i].startA; currentEnd=windows[i].endA;
				continue;
			}
			if (windows[i].startA>currentEnd+IDENTITY_THRESHOLD) {
				if (last+2>=tmpArray1.length) {
					int[] newArray = new int[tmpArray1.length<<1];
					System.arraycopy(tmpArray1,0,newArray,0,last+1);
					tmpArray1=newArray;
				}
				tmpArray1[++last]=currentStart; tmpArray1[++last]=currentEnd;
				currentStart=windows[i].startA; currentEnd=windows[i].endA;
			}
			else currentEnd=Math.max(currentEnd,windows[i].endA);
		}
		if (currentStart!=-1) {
			if (last+2>=tmpArray1.length) {
				int[] newArray = new int[tmpArray1.length<<1];
				System.arraycopy(tmpArray1,0,newArray,0,last+1);
				tmpArray1=newArray;
			}
			tmpArray1[++last]=currentStart; tmpArray1[++last]=currentEnd;
		}
		currentSurface=Reads.getReadLength(currentRead);
		for (i=0; i<last; i+=2) currentSurface-=tmpArray1[i+1]-tmpArray1[i]+1;
		
		// Checking other covered surfaces
		if (currentSurface>=SURFACE_THRESHOLD && !currentReference && !isCyclicAndNonperiodic) {
			stats[0]++;
			System.err.println("Fragment "+currentRead+" does not fully align to the reference for "+currentSurface+" bp");
			currentSurface/=IO.quantum;
			currentSurface=Math.min(currentSurface,frUncoveredHistogram.length-1);
			frUncoveredHistogram[currentSurface]++;
			ffUncoveredSurface=uncoveredSurface_fragment2fragment(currentRead,tmpArray1,last,newFragmentCoordinates[currentRead]);
			// Remark: it can happen that $ffUncoveredSurface$ is small, but the output of
			// $uncoveredSurface_read2read()$ is large, if the fragment is a substring
			// around the reference (this might happen even if $currentReference=false$).
			if (ffUncoveredSurface>=SURFACE_THRESHOLD) {
				stats[1]++;
				System.err.println("   Fragment "+currentRead+" does not fully align to the reference nor to other fragments for "+ffUncoveredSurface+" bp");
				ffUncoveredSurface/=IO.quantum;
				ffUncoveredSurface=Math.min(ffUncoveredSurface,ffUncoveredHistogram.length-1);
				ffUncoveredHistogram[ffUncoveredSurface]++;
				rrUncoveredSurface=uncoveredSurface_read2read(currentRead,minIntersectionLength);
				if (rrUncoveredSurface>=SURFACE_THRESHOLD && !(basinType==Constants.INTERVAL_PERIODIC && period<PeriodicSubstrings.MIN_PERIOD_LONG && Reads.getReadLength(currentRead)>referenceLength)) {
					stats[2]++;
					System.err.println("      Fragment "+currentRead+" does not fully align to the reference, nor to other fragments, nor to basin intervals using read-read alignments, for "+rrUncoveredSurface+" bp ("+fileName+").");
					rrUncoveredSurface/=IO.quantum;
					rrUncoveredSurface=Math.min(rrUncoveredSurface,rrUncoveredHistogram.length-1);
					rrUncoveredHistogram[rrUncoveredSurface]++;
				}
			}
		}
		else {
			// Done just to update $newFragmentCoordinates$.
			uncoveredSurface_fragment2fragment(currentRead,tmpArray1,last,newFragmentCoordinates[currentRead]);
		}
	}
	
	
	/**
	 * Assume that $windows$ contains all and only the aligments of a given fragment $A$ 
	 * to the reference $B$. The procedure tries to detect if $A$ is a permutation of $B$,
	 * by finding a pair of not necessarily disjoint substrings $A[i..j],A[i'..j']$ of 
	 * $A$, with $i'>i$, that are aligned to substrings $B[p..q],B[p'..q']$ of $B$, with 
	 * $p' <= p$, and such that $A[i..j],A[i'..j']$ are never aligned to any pair of 
	 * substrings $B[x..y],B[x'..y']$ with $x'>x$.
	 * This is essentially identical to $IntervalGraphStep3.buildKernelGraph_
	 * removePermutations()$.
	 *
	 * Remark: for simplicity, the procedure does not check for permutations with RC, i.e. 
 	 * it considers alignments in each orientation separately.
	 *
	 * Remark: the procedure uses alignments of fragments to the reference, so it does not
	 * work for any fragment that does not align to the reference. One could think of 
	 * looking for peaks in the alignments of a fragment to every other fragment instead.
	 * However, since a kernel might contain other repeats, such peaks are not guaranted 
	 * to appear only for permutations/deletions.
	 *
	 * Remark: the procedure handles the case in which $B[p..q]$ and $B[p'..q']$ coincide.
	 * Remark: the procedure should not be called for short-period kernels.
	 *
	 * @param windows not assumed to be sorted in any way; the procedure might change the
	 * order;
	 * @param period the period of the long-period kernel, or -1 if it does not have a 
	 * long period;
	 * @param tmp temporary space with at least two cells.
	 */
	private static final boolean isPermutation(PermutationWindow[] windows, int lastWindow, int period, int[] tmp) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int DISTANCE_THRESHOLD = IO.quantum;  // Arbitrary
		boolean orientationI, found;
		int i, j, iPrime, jPrime;
		int startAI, endAI, startBI, endBI, startAJ, endAJ, startBJ, endBJ;
		int startIPrime, endIPrime, startJPrime, endJPrime;
		
		PermutationWindow.order=PermutationWindow.ORDER_READB_ORIENTATION_STARTA_STARTB;
		if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
		for (i=0; i<lastWindow; i++) {
			orientationI=windows[i].orientation;
			startAI=windows[i].startA; endAI=windows[i].endA;
			if (period>0 && endAI-startAI+1>period+IDENTITY_THRESHOLD) continue;
			startBI=windows[i].startB; endBI=windows[i].endB;
			for (j=i+1; j<=lastWindow; j++) {
				if (windows[j].orientation!=orientationI) break;
				startAJ=windows[j].startA;
				if (period>0 && startAJ>=startAI+period-IDENTITY_THRESHOLD) break;
				if (startAJ<=startAI+DISTANCE_THRESHOLD) continue;
				endAJ=windows[j].endA;
				if (period>0 && endAJ-startAJ+1>period+IDENTITY_THRESHOLD) continue;
				if (period>0 && endAJ-startAI+1>period+IDENTITY_THRESHOLD) continue;
				if (Intervals.isApproximatelyContained(startAJ,endAJ,startAI,endAI) || Intervals.isApproximatelyContained(startAI,endAI,startAJ,endAJ)) continue;
				startBJ=windows[j].startB; endBJ=windows[j].endB;
				if (Intervals.isApproximatelyContained(startBJ,endBJ,startBI,endBI) || Intervals.isApproximatelyContained(startBI,endBI,startBJ,endBJ)) continue;
				if (period>0 && Intervals.unionLength(startBI,endBI,startBJ,endBJ)>period+IDENTITY_THRESHOLD) continue;
				if ((orientationI && startBJ>startBI+DISTANCE_THRESHOLD) || (!orientationI && endBJ<endBI-DISTANCE_THRESHOLD)) continue;
				found=false;
				iPrime=i+1;
				while (iPrime<=lastWindow && windows[iPrime].orientation==orientationI && windows[iPrime].startA==startAI) iPrime++;
				iPrime--;
				while (iPrime>=0) {
					if (iPrime==i)  {
						iPrime--;
						continue;
					}
					if (windows[iPrime].orientation!=orientationI) break;
					if (!Intervals.isApproximatelyContained(startAI,endAI,windows[iPrime].startA,windows[iPrime].endA)) {
						iPrime--;
						continue;
					}
					Intervals.project(startAI,endAI,windows[iPrime].startA,windows[iPrime].endA,windows[iPrime].startB,windows[iPrime].endB,windows[iPrime].orientation,tmp,0);
					startIPrime=tmp[0]; endIPrime=tmp[1];
					jPrime=j+1;
					while (jPrime<=lastWindow && windows[jPrime].orientation==orientationI && windows[jPrime].startA==startAJ) jPrime++;
					jPrime--;
					while (jPrime>=0) {
						if (jPrime==j)  {
							jPrime--;
							continue;
						}
						if (windows[jPrime].orientation!=orientationI) break;
						if (!Intervals.isApproximatelyContained(startAJ,endAJ,windows[jPrime].startA,windows[jPrime].endA)) {
							jPrime--;
							continue;
						}
						Intervals.project(startAJ,endAJ,windows[jPrime].startA,windows[jPrime].endA,windows[jPrime].startB,windows[jPrime].endB,windows[jPrime].orientation,tmp,0);
						startJPrime=tmp[0]; endJPrime=tmp[1];
						if ( (orientationI && startJPrime>startIPrime+DISTANCE_THRESHOLD && (period>0?Intervals.unionLength(startIPrime,endIPrime,startJPrime,endJPrime)<=period+IDENTITY_THRESHOLD:true)) || 
							 (!orientationI && endJPrime<endIPrime-DISTANCE_THRESHOLD && (period>0?Intervals.unionLength(startIPrime,endIPrime,startJPrime,endJPrime)<=period+IDENTITY_THRESHOLD:true))
						   ) {
							found=true;
							break;
						}
						jPrime--;
					}
					if (found) break;
					iPrime--;
				}
				if (!found) return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Consider all alignments in a given orientation. The procedure returns TRUE if there
	 * is a large gap in the readB intervals of two alignments, such that: (1) the gap 
	 * maps to approx. a single position P in readA; (2) there is no alignment, in the 
	 * same orientation, that contains P in readA.
	 * This is similar to $IntervalGraphStep3.buildKernelGraph_removeDeletions()$.
	 *
	 * Remark: the procedure uses $tmpArray$ as temporary space.
	 * 
	 * @param windows not assumed to be sorted in any way; the procedure might change the
	 * order;
	 * @param period <=0: the kernel is non-periodic; >0: the kernel has this long period.
	 */
	private static final boolean isDeletion(PermutationWindow[] windows, int lastWindow, int period) {
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int MAX_GAP = IO.quantum<<1;  // Arbitrary
		final int VALIDATION_DISTANCE = Alignments.minAlignmentLength>>1;  // Arbitrary
		boolean found, intervalOrientation;
		int i, j;
		int intervalStartB, intervalEndB, intervalEndA, lastDeletion, firstJForNextI;
		
		// Finding candidates
		PermutationWindow.order=PermutationWindow.ORDER_READB_ORIENTATION_STARTB;
		if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
		intervalOrientation=windows[0].orientation;
		intervalStartB=windows[0].startB; intervalEndB=windows[0].endB;
		intervalEndA=intervalOrientation?windows[0].endA:windows[0].startA;
		lastDeletion=-1;
		for (i=1; i<=lastWindow; i++) {
			if (windows[i].orientation!=intervalOrientation) {
				intervalOrientation=windows[i].orientation;
				intervalStartB=windows[i].startB; intervalEndB=windows[i].endB; 
				intervalEndA=intervalOrientation?windows[i].endA:windows[i].startA;
				continue;
			}
			if (windows[i].startB<=intervalEndB+MAX_GAP) {
				intervalEndB=Math.max(intervalEndB,windows[i].endB);
				intervalEndA=intervalOrientation?Math.max(intervalEndA,windows[i].endA):Math.min(intervalEndA,windows[i].startA);
			}
			else {
				if (Math.abs(intervalEndA,intervalOrientation?windows[i].startA:windows[i].endA)<=IDENTITY_THRESHOLD) {
					if (period>0) {
						if (windows[i].endB-intervalStartB+1<=period+IDENTITY_THRESHOLD) {
							if (lastDeletion+2>=tmpArray.length) {
								int[] newArray = new int[tmpArray.length<<1];
								System.arraycopy(tmpArray,0,newArray,0,tmpArray.length);
								tmpArray=newArray;
							}
							tmpArray[++lastDeletion]=intervalOrientation?0:1;
							tmpArray[++lastDeletion]=intervalEndA;
						}
					}
					else {
						if (lastDeletion+2>=tmpArray.length) {
							int[] newArray = new int[tmpArray.length<<1];
							System.arraycopy(tmpArray,0,newArray,0,tmpArray.length);
							tmpArray=newArray;
						}
						tmpArray[++lastDeletion]=intervalOrientation?0:1;
						tmpArray[++lastDeletion]=intervalEndA;
					}
				}
				intervalStartB=windows[i].startB; intervalEndB=windows[i].endB;
				intervalEndA=intervalOrientation?windows[i].endA:windows[i].startA;
			}
		}
		if (lastDeletion==-1) return false;
		
		// Validating candidates
		PermutationWindow.order=PermutationWindow.ORDER_READB_ORIENTATION_STARTA_STARTB;
		if (lastWindow>0) Arrays.sort(windows,0,lastWindow+1);
		j=0; firstJForNextI=-1;
		for (i=0; i<lastDeletion; i+=2) {
			found=false;
			while (j<=lastWindow) {				
				if (firstJForNextI==-1 && i<lastDeletion-1 && windows[j].orientation==(tmpArray[i+2]==0) && tmpArray[i+3]>=windows[j].startA && tmpArray[i+3]<=windows[j].endA) firstJForNextI=j;
				if ((!windows[j].orientation && tmpArray[i]==0) || windows[j].startA>tmpArray[i+1]) break;
				if (tmpArray[i+1]>=windows[j].startA+VALIDATION_DISTANCE && tmpArray[i+1]<=windows[j].endA-VALIDATION_DISTANCE) {
					found=true;
					break;
				}
				j++;
			}
			if (!found) return true;
			if (firstJForNextI!=-1) {
				j=firstJForNextI;
				firstJForNextI=-1;
			}
		}
		return false;
	}
	
	
	public static class PermutationWindow implements Comparable {
		public static final int ORDER_READB_ORIENTATION_STARTA_STARTB = 0;
		public static final int ORDER_READB_ORIENTATION_STARTB = 1;
		public static final int ORDER_STARTA = 2;
		public static final int ORDER_READB_STARTB = 3;
		public static int order;
		public boolean orientation;
		public int readB, startA, endA, startB, endB, diffs;
		
		public PermutationWindow() {
			this.readB=-1;
			this.orientation=false;
			this.startA=-1;
			this.endA=-1;
			this.startB=-1;
			this.endB=-1;
			this.diffs=0;
		}
		
		public PermutationWindow(int rb, boolean or, int s, int e, int sb, int eb, int d) {
			set(rb,or,s,e,sb,eb,d);
		}
		
		public void set(int rb, boolean or, int s, int e, int sb, int eb, int d) {
			readB=rb;
			orientation=or;
			startA=s;
			endA=e;
			startB=sb;
			endB=eb;
			diffs=d;
		}
		
		public int compareTo(Object other) {
			PermutationWindow otherWindow = (PermutationWindow)other;
			if (order==ORDER_READB_ORIENTATION_STARTA_STARTB) {
				if (readB<otherWindow.readB) return -1;
				else if (readB>otherWindow.readB) return 1;
				if (orientation && !otherWindow.orientation) return -1;
				else if (!orientation && otherWindow.orientation) return 1;
				if (startA<otherWindow.startA) return -1;
				else if (startA>otherWindow.startA) return 1;
				if (startB<otherWindow.startB) return -1;
				else if (startB>otherWindow.startB) return 1;
			}
			else if (order==ORDER_READB_ORIENTATION_STARTB) {
				if (readB<otherWindow.readB) return -1;
				else if (readB>otherWindow.readB) return 1;
				if (orientation && !otherWindow.orientation) return -1;
				else if (!orientation && otherWindow.orientation) return 1;
				if (startB<otherWindow.startB) return -1;
				else if (startB>otherWindow.startB) return 1;
			}
			else if (order==ORDER_STARTA) {
				if (startA<otherWindow.startA) return -1;
				else if (startA>otherWindow.startA) return 1;
			}
			else if (order==ORDER_READB_STARTB) {
				if (readB<otherWindow.readB) return -1;
				else if (readB>otherWindow.readB) return 1;
				if (startB<otherWindow.startB) return -1;
				else if (startB>otherWindow.startB) return 1;
			}
			return 0;
		}
		
		public boolean equals(Object other) {
			PermutationWindow otherWindow = (PermutationWindow)other;
			return orientation==otherWindow.orientation && startA==otherWindow.startA && endA==otherWindow.endA && startB==otherWindow.startB && endB==otherWindow.endB;
		}
		
		/**
		 * @param from first position of $array$.
		 */
		public void writeTo(int[] array, int from) {
			array[from+0]=orientation?0:1;
			array[from+1]=startA;
			array[from+2]=endA;
			array[from+3]=startB;
			array[from+4]=endB;
		}
		
		public String toString() {
			return (orientation?0:1)+"["+startA+".."+endA+"] x ["+startB+".."+endB+"]";
		}
	}
	
	
	/**
	 * Writes the artificial connection file that will be used by $FragmentsPrintGraph$ to
	 * build a fragment-fragment alignment graph. Fragment intervals with just a prefix or 
	 * suffix that does not align to any other fragment are trimmed. Fragments with 
	 * internal regions that do not align to any other fragment are discarded.
	 */
	private static final void writeArtificialConnectionFile(String alignmentsFile, String connectionFile, int[][] newFragmentCoordinates) throws IOException {
		int readA;
		String str;
		BufferedWriter bw;
		BufferedReader br;
		
		bw = new BufferedWriter(new FileWriter(connectionFile));
		br = new BufferedReader(new FileReader(alignmentsFile),IO.BUFFER_SIZE);
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			readA=Alignments.readA-1;
			if (newFragmentCoordinates[readA][0]==-1) bw.write("\n");
			else bw.write(Constants.INTERVAL_ALIGNMENT+","+readA+","+newFragmentCoordinates[readA][0]+","+newFragmentCoordinates[readA][1]+",0,0,0,0,0,1,1\n");
			str=br.readLine();
		}
		bw.close();
	}
	
	
	
	
	
	
	
	
	// ----------------------------- FRAGMENTS PROCEDURES --------------------------------
	
	private static final void loadFragmentCoordinates(String path) throws IOException {
		int i, p, q;
		String str;
		BufferedReader br;

		fragmentCoordinates = new int[3][N_FRAGMENTS];
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();
		for (i=0; i<N_FRAGMENTS; i++) {
			p=str.indexOf(",");
			fragmentCoordinates[0][i]=Integer.parseInt(str.substring(0,p));
			p++; q=str.indexOf(",",p+1);
			fragmentCoordinates[1][i]=Integer.parseInt(str.substring(p,q));
			p=q+1;
			fragmentCoordinates[2][i]=Integer.parseInt(str.substring(p));
			str=br.readLine();
		}
		br.close();
	}
	
	
	/**
	 * Loads all alignments in $fragments2fragments,fragmentFirst$.
	 */
	private static final void loadFragments2Fragments(String path) throws IOException {
		int i;
		int current, previous, previousFirst;
		String str;
		BufferedReader br;
		
		// Loading $fragments2fragments$.
		fragments2fragments = new FFAlignment[N_FF_ALIGNMENTS];
		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		for (i=0; i<N_FF_ALIGNMENTS; i++) {
			Alignments.readAlignmentFile(str);
			fragments2fragments[i] = new FFAlignment(Alignments.readA-1,Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
			str=br.readLine();
		}
		br.close();
		if (N_FF_ALIGNMENTS>1) {
			FFAlignment.order=FFAlignment.ORDER_READA_STARTA;
			Arrays.sort(fragments2fragments);
		}
		
		// Loading $fragmentFirst$.
		fragmentFirst = new int[N_FRAGMENTS];
		Math.set(fragmentFirst,N_FRAGMENTS-1,-1);
		previous=-1; previousFirst=-1;
		for (i=0; i<N_FF_ALIGNMENTS; i++) {
			current=fragments2fragments[i].readA;
			if (current!=previous) {
				if (previous!=-1) fragmentFirst[previous]=previousFirst;
				previous=current; previousFirst=i;
			}
		}
		if (previous!=-1) fragmentFirst[previous]=previousFirst;
	}
	
	
	/**
	 * Remark: the procedure uses global variables $tmpArray{2,3}$ as temporary space.
	 *
	 * @param fragment2reference list of fragment-reference alignments (ignored if null);
	 * @param out output array: if $fragment$ contains exactly one contiguous region 
	 * covered by fragment-{fragment,reference} alignments (in any orientation), the array
	 * stores the first and last such covered positions; otherwise, it stores (-1,-1);
	 * @return the surface of $fragment$ that is not covered by any fragment-{fragment,
	 * reference} alignment.
	 */
	private static final int uncoveredSurface_fragment2fragment(int fragment, int[] fragment2reference, int lastFragment2reference, int[] out) {
		int i;
		int startA, endA, currentStart, currentEnd, currentSurface, last, length;
		final int fragmentLength = Reads.getReadLength(fragment);
		
		if (fragmentFirst[fragment]==-1) {
			out[0]=-1; out[1]=-1;
			return fragmentLength;
		}
		
		// Building fragment-fragment track
		currentStart=-1; currentEnd=-1; last=-1;
		for (i=fragmentFirst[fragment]; i<N_FF_ALIGNMENTS; i++) {
			if (fragments2fragments[i].readA!=fragment) break;
			startA=fragments2fragments[i].startA; endA=fragments2fragments[i].endA;
			if (currentStart==-1) {
				currentStart=startA; currentEnd=endA;
				continue;
			}
			if (startA>currentEnd+IDENTITY_THRESHOLD) {
				if (last+2>=tmpArray3.length) {
					int[] newArray = new int[tmpArray3.length<<1];
					System.arraycopy(tmpArray3,0,newArray,0,last+1);
					tmpArray3=newArray;
				}
				tmpArray3[++last]=currentStart; tmpArray3[++last]=currentEnd;
				currentStart=startA; currentEnd=endA;
			}
			else currentEnd=Math.max(currentEnd,endA);
		}
		if (currentStart!=-1) {
			if (last+2>=tmpArray3.length) {
				int[] newArray = new int[tmpArray3.length<<1];
				System.arraycopy(tmpArray3,0,newArray,0,last+1);
				tmpArray3=newArray;
			}
			tmpArray3[++last]=currentStart; tmpArray3[++last]=currentEnd;
		}
		
		// Merging with fragment-reference track, if any.
		if (fragment2reference!=null && lastFragment2reference>=0) {
			length=last+lastFragment2reference+2;
			if (tmpArray2.length<length) tmpArray2 = new int[length];
			last=Intervals.union(tmpArray3,last,fragment2reference,lastFragment2reference,tmpArray2);
		}
		else {
			int[] array = tmpArray2;
			tmpArray2=tmpArray3;
			tmpArray3=array;
		}
		
		// Computing statistics
		currentSurface=fragmentLength;
		for (i=0; i<last; i+=2) currentSurface-=tmpArray2[i+1]-tmpArray2[i]+1;
		if (last!=1) { out[0]=-1; out[1]=-1; }
		else { out[0]=tmpArray2[0]; out[1]=tmpArray2[1]; }
		return currentSurface;
	}

	
	
	
	


	
	// --------------------------- ALL ALIGNMENTS PROCEDURES -----------------------------
	
	private static final void loadIntervalCoordinates(String path) throws IOException {
		int i, p, q;
		String str;
		BufferedReader br;

		intervalCoordinates = new int[3][N_INTERVALS];
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();  // Skipping header
		str=br.readLine();
		for (i=0; i<N_INTERVALS; i++) {
			p=str.indexOf(",");
			intervalCoordinates[0][i]=Integer.parseInt(str.substring(0,p));
			p++; q=str.indexOf(",",p+1);
			intervalCoordinates[1][i]=Integer.parseInt(str.substring(p,q));
			p=q+1; q=str.indexOf(",",p+1);
			intervalCoordinates[2][i]=Integer.parseInt(str.substring(p,q));
			str=br.readLine();
		}
		br.close();
	}
	
	
	private static final void loadReads2reads(String path) throws IOException {
		int i;
		String str;
		BufferedReader br;

		reads2reads = new FFAlignment[N_RR_ALIGNMENTS];
		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		for (i=0; i<N_RR_ALIGNMENTS; i++) {
			Alignments.readAlignmentFile(str);
			reads2reads[i] = new FFAlignment(Alignments.readA-1,Alignments.readB-1,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
			str=br.readLine();
		}
		br.close();
		if (N_RR_ALIGNMENTS>1) {
			FFAlignment.order=FFAlignment.ORDER_READA_STARTA;
			Arrays.sort(reads2reads);  // Sorting by $readA,startA$.
		}
	}
	
	
	/**
	 * Remark: it can happen that a fragment has a long substring that is not covered by
	 * alignments between intervals. This might be caused by e.g. $IntervalGraphStep3.
	 * checkKernelsWithReads()$ and $IntervalGraphStep3.getBasins_ensureKernelLengths()$.
	 *
	 * Remark: the procedure uses global variables $tmpAlignment,tmpArray,tmpPairs$, which
	 * are assumed to be already initialized and large enough.
	 * 
	 * @param minIntersectionLength intersections shorter than this, between the readB 
	 * side of an alignment and another fragment, are not used;
	 * @return the surface of $fragment$ that is not covered by any fragment-interval 
	 * alignment induced by a read-read alignment.
	 */
	private static final int uncoveredSurface_read2read(int fragment, int minIntersectionLength) {
		int i, j;
		int last, lastPair, startA, endA, currentStart, currentEnd, currentSurface;
		final int fragmentRead  = fragmentCoordinates[0][fragment];
		final int fragmentStart  = fragmentCoordinates[1][fragment];
		final int fragmentEnd  = fragmentCoordinates[2][fragment];
		final int fragmentLength = fragmentEnd-fragmentStart+1;
		
		// Collecting all induced projections of an interval onto $fragment$.
		tmpAlignment.readA=fragmentRead; tmpAlignment.startA=fragmentStart;
		i=Arrays.binarySearch(reads2reads,0,N_RR_ALIGNMENTS,tmpAlignment);
		if (i<0) i=-i-1;
		if (i==N_RR_ALIGNMENTS || reads2reads[i].readA>fragmentRead) i--;
		while (i>=0 && reads2reads[i].readA==fragmentRead) i--;
		i++; lastPair=-1;	
		while (i<N_RR_ALIGNMENTS) {
			if (reads2reads[i].readA!=fragmentRead || reads2reads[i].startA>=fragmentEnd) break;
			if (reads2reads[i].endA<=fragmentStart) {
				i++;
				continue;
			}
			last=intersectsIntervals(reads2reads[i].readB,reads2reads[i].startB,reads2reads[i].endB,minIntersectionLength);
			for (j=0; j<last; j+=2) lastPair=projectOntoReadA(i,tmpArray[j],tmpArray[j+1],lastPair);
			i++;
		}
		if (lastPair>0) Arrays.sort(tmpPairs,0,lastPair+1);
		
		// Computing uncovered surface
		currentStart=-1; currentEnd=-1; currentSurface=fragmentLength;
		for (i=0; i<=lastPair; i++) {
			startA=tmpPairs[i].start; endA=tmpPairs[i].end;
			if (currentStart==-1) {
				currentStart=startA; currentEnd=endA;
				continue;
			}
			if (startA>currentEnd+IDENTITY_THRESHOLD) {
				currentSurface-=currentEnd-currentStart+1;
				currentStart=startA; currentEnd=endA;
			}
			else currentEnd=Math.max(currentEnd,endA);
		}
		if (currentStart!=-1) currentSurface-=currentEnd-currentStart+1;
		return currentSurface;
	}
	
	
	/**
	 * @param out output array, storing all intervals that result from the intersection of
	 * $[start..end]$ with an element of $fragmentCoordinates$;
	 * @param minIntersectionLength only intersections of at least this length are 
	 * reported;
	 * @return the last element in $out$, or -1 if $[start..end]$ does not intersect with
	 * any interval.
	 */
	private static final int intersectsFragments(int read, int start, int end, int minIntersectionLength, int[] out) {
		int i;
		int last, from, to;
	
		last=-1;
		i=Arrays.binarySearch(fragmentCoordinates[0],read);
		if (i<0) return last;
		while (i>=0 && fragmentCoordinates[0][i]==read) i--;
		i++;
		while (i<N_FRAGMENTS) {
			if (fragmentCoordinates[0][i]!=read) break;
			from=Math.max(start,fragmentCoordinates[1][i]);
			to=Math.min(end,fragmentCoordinates[2][i]);
			if (to-from+1>=minIntersectionLength) { out[++last]=from; out[++last]=to; }
			i++;
		}
		return last;
	}
	
	
	/**
	 * Stores in global variable $tmpArray$ all intervals that result from the 
	 * intersection of $[start..end]$ with an element of $intervalCoordinates$.
	 *
	 * @param minIntersectionLength only intersections of at least this length are 
	 * reported;
	 * @return the last element in $out$, or -1 if $[start..end]$ does not intersect with
	 * any interval.
	 */
	private static final int intersectsIntervals(int read, int start, int end, int minIntersectionLength) {
		final int GROWTH_RATE = 100;  // Arbitrary
		int i;
		int last, from, to;
	
		last=-1;
		i=Arrays.binarySearch(intervalCoordinates[0],read);
		if (i<0) return last;
		while (i>=0 && intervalCoordinates[0][i]==read) i--;
		i++;
		while (i<N_INTERVALS) {
			if (intervalCoordinates[0][i]!=read) break;
			from=Math.max(start,intervalCoordinates[1][i]);
			to=Math.min(end,intervalCoordinates[2][i]);
			if (to-from+1>=minIntersectionLength) { 
				if (last+2>=tmpArray.length) {
					int[] newArray = new int[tmpArray.length+GROWTH_RATE];
					System.arraycopy(tmpArray,0,newArray,0,last+1);
					tmpArray=newArray;
				}
				tmpArray[++last]=from; tmpArray[++last]=to;
			}
			i++;
		}
		return last;
	}
	
	
	/**
	 * @param startB,endB assumed to be included in the readB side of the alignment;
	 * @param out appends to the end of the global variable $tmpPairs$ a new pair that 
	 * represents the projection;
	 * @return the new value of $last$.
	 */
	private static final int projectOntoReadA(int alignmentID, int startB, int endB, int last) {
		final int GROWTH_RATE = 1000;  // Arbitrary
		final int alignmentStartA = reads2reads[alignmentID].startA;
		final int alignmentEndA = reads2reads[alignmentID].endA;
		final int alignmentStartB = reads2reads[alignmentID].startB;
		final int alignmentEndB = reads2reads[alignmentID].endB;
		final double ratio = ((double)(alignmentEndA-alignmentStartA+1))/(alignmentEndB-alignmentStartB+1);
		int startA, endA;
		
		if (reads2reads[alignmentID].orientation) {
			startA=alignmentStartA+(int)((startB-alignmentStartB)*ratio);
			endA=alignmentStartA+(int)((endB-alignmentStartB)*ratio);
		}
		else {
			startA=alignmentStartA+(int)((alignmentEndB-endB)*ratio);
			endA=alignmentStartA+(int)((alignmentEndB-startB)*ratio);
		}
		last++;
		if (last==tmpPairs.length) {
			Pair[] newOut = new Pair[tmpPairs.length+GROWTH_RATE];
			System.arraycopy(tmpPairs,0,newOut,0,tmpPairs.length);
			for (int i=tmpPairs.length; i<newOut.length; i++) newOut[i] = new Pair();
			tmpPairs=newOut;
		}
		tmpPairs[last].start=startA; tmpPairs[last].end=endA;
		return last;
	}

	
	
	
	


	
	// ------------------------------- DATA STRUCTURES -----------------------------------
	
	public static class FFAlignment implements Comparable {
		public static final byte ORDER_READA_STARTA = 0;
		public static final byte ORDER_READB_STARTB = 1;
		public static final byte ORDER_ID = 2;
		public static final byte ORDER_NNEIGHBORS = 3;
		
		public static byte order;
		public boolean orientation;
		public int readA, readB, startA, endA, startB, endB;
		public int id, score, nNeighbors;  // Used by $FragmentsStep3$.
		
		public FFAlignment() {
			orientation=false;
			readA=-1; readB=-1; startA=-1; endA=-1; startB=-1; endB=-1;
			id=-1; score=-1;
		}	
		
		public FFAlignment(int readA, int readB, boolean orientation, int startA, int endA, int startB, int endB) {
			set(readA,readB,orientation,startA,endA,startB,endB);
		}
		
		public void set(int readA, int readB, boolean orientation, int startA, int endA, int startB, int endB) {
			this.readA=readA; this.readB=readB; this.orientation=orientation;
			this.startA=startA; this.endA=endA; this.startB=startB; this.endB=endB;
		}
		
		public int compareTo(Object other) {
			FFAlignment otherAlignment = (FFAlignment)other;
			if (order==ORDER_READA_STARTA) {
				if (readA<otherAlignment.readA) return -1;
				else if (readA>otherAlignment.readA) return 1;
				else if (startA<otherAlignment.startA) return -1;
				else if (startA>otherAlignment.startA) return 1;
			}
			else if (order==ORDER_READB_STARTB) {
				if (readB<otherAlignment.readB) return -1;
				else if (readB>otherAlignment.readB) return 1;
				else if (startB<otherAlignment.startB) return -1;
				else if (startB>otherAlignment.startB) return 1;
			}
			else if (order==ORDER_ID) {
				if (id<otherAlignment.id) return -1;
				else if (id>otherAlignment.id) return 1;
			}
			else if (order==ORDER_NNEIGHBORS) {
				if (nNeighbors<otherAlignment.nNeighbors) return -1;
				else if (nNeighbors>otherAlignment.nNeighbors) return 1;
			}
			return 0;
		}
		
		public boolean equals(Object other) {
			FFAlignment otherAlignment = (FFAlignment)other;
			return readA==otherAlignment.readA && startA==otherAlignment.startA && endA==otherAlignment.endA && readB==otherAlignment.readB && startB==otherAlignment.startB && endB==otherAlignment.endB && orientation==otherAlignment.orientation;
		}
		
		public String toString() {
			return id+": "+readA+"["+startA+".."+endA+"] x "+readB+"["+startB+".."+endB+"] "+orientation+" // "+score;
		}
		
		public void inCanonicalForm() {
			if (readA<readB) return;
			int tmp;
			tmp=readA; readA=readB; readB=tmp;
			tmp=startA; startA=startB; startB=tmp;
			tmp=endA; endA=endB; endB=tmp;
		}
		
		public int hashCode() {
			return (readA+"-"+startA+"-"+endA+"-"+readB+"-"+startB+"-"+endB+"-"+(orientation?"1":"0")).hashCode();
		}
	}
	
	
	public static class Pair implements Comparable {
		public int start, end;
		
		public Pair() {
			start=-1; end=-1;
		}
		
		public Pair(int s, int e) {
			this.start=s; this.end=e;
		}
		
		public int compareTo(Object other) {
			Pair otherPair = (Pair)other;
			if (start<otherPair.start) return -1;
			else if (start>otherPair.start) return 1;
			return 0;
		}
		
		public String toString() {
			return "("+start+","+end+")";
		}
	}

}