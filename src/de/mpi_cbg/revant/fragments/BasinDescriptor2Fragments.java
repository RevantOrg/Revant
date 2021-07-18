package de.mpi_cbg.revant.fragments;

import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.PeriodicSubstrings;
import de.mpi_cbg.revant.intervalgraph.IntervalGraph;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;


/**
 * Prints a description of the set of fragments to be used for computing the consensus of
 * a path kernel specified in a descriptor file built by
 * $IntervalGraphStep3.printBasinDescriptors()$.
 * 
 * - Descriptors containing too few intervals are not printed. 
 * - Intervals used for building the reference sequence of the descriptor (and all 
 * intervals contained in them) can be omitted. Note that this can generate fragments that
 * \emph{contain} the reference as a substring.
 * - Intervals in the same read are merged according to the criteria in 
 * $IntervalGraphStep3.printKernelTags_shouldBeMerged()$. Those criteria might create 
 * merges that are longer than the reference sequence of the descriptor (e.g. if several 
 * intervals have a suffix-prefix overlap).
 * - Only the maximal disjoint high-quality subintervals of the merges are printed.
 */
public class BasinDescriptor2Fragments {
	
	private static int[][] forbidden = new int[20][3];  // $read,start,end$
	private static int lastForbidden;
	private static int[] tmpArray = new int[200];
	
	
	/**
	 * @param args 
	 * 0: Step1 directory;
	 * 1: min. number of fragments for a descriptor file to be output at all; 
	 * 2: min alignment length used during repeat inference;
	 * 3: (0/1) outputs only fragments that belong to exactly one kernel;
	 * 4: (0/1) prints also the fragments used for building the reference.
	 */
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		final int MIN_N_FRAGMENTS = Integer.parseInt(args[1]);
		final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[2]);
		final boolean ONE_KERNEL_ONLY = Integer.parseInt(args[3])==1;
		final boolean PRINT_REFERENCE_FRAGMENTS = Integer.parseInt(args[4])==1;
		
		final String STEP4_DIR = STEP1_DIR+"/finalOutput/step4";
		final String STEP5_DIR = STEP4_DIR+"/step5";
		final String BASIN_DESCRIPTOR_FILES = STEP5_DIR+"/descriptors.txt";
		final String INPUT_DIR = STEP4_DIR;
		final String OUTPUT_DIR = STEP5_DIR+"/fragments";
		final String READ_LENGTHS_FILE = STEP1_DIR+"/../reads-lengths.txt";
		final String READ_IDS_FILE = STEP1_DIR+"/../reads-ids.txt";
		final String QUALITY_THRESHOLDS_FILE = STEP1_DIR+"/../qualityThresholds.txt";
		String QUALITIES_FILE = STEP1_DIR+"/../"+IO.READS_PHRED_PREFIX+IO.getPhredSuffix(STEP1_DIR+"/../"+IO.READS_PHRED_PREFIX);
		File tmpFile = new File(QUALITIES_FILE);
		if (!tmpFile.exists()) QUALITIES_FILE=null;
		tmpFile=null;
		
		final int QUALITY_WINDOW = MIN_ALIGNMENT_LENGTH>>1;  // Arbitrary
		final int IDENTITY_THRESHOLD = IO.quantum;
		final int BASIN_LABEL_LENGTH = IO.BASIN_DESCRIPTOR_PREFIX.length()+1;
		boolean merged;
		int i;
		int read, start, end, currentRead, maxReadLength;
		int bidirectedGraphNode, type, period, nFragments, otherKernels, lastWindow;
		String str, file, fileID;
		BufferedReader br, sequenceInfoFiles;
		BufferedWriter bw;
		IntervalGraph.Node tmpNode = new IntervalGraph.Node();
		int[] pathKernelLengths = new int[1];  // Artificial array of size one
		byte[] pathKernelPeriodic = new byte[1];  // Artificial array of size one
		IntervalGraph.Node[] windows = null;
		
		IO.initialize();
		Alignments.minAlignmentLength=MIN_ALIGNMENT_LENGTH;
		br = new BufferedReader(new FileReader(READ_LENGTHS_FILE));
		i=0; str=br.readLine();
		while (str!=null) {
			i++;
			str=br.readLine();
		}
		br.close();
		Reads.nReads=i;
		Reads.loadReadIDs(READ_IDS_FILE,i);
		maxReadLength=Reads.loadReadLengths(READ_LENGTHS_FILE);
		Reads.maxReadLength=maxReadLength;
		Reads.loadQualities(QUALITY_THRESHOLDS_FILE,QUALITIES_FILE);
		sequenceInfoFiles = new BufferedReader(new FileReader(BASIN_DESCRIPTOR_FILES));
		file=sequenceInfoFiles.readLine();
		while (file!=null) {
			System.out.println("Printing fragments of file "+file);
			
			// Checking the number of fragments
			br = new BufferedReader(new FileReader(INPUT_DIR+"/"+file.trim()));
			str=br.readLine();
			IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
			nFragments=tmpArray[5];
			if (nFragments<MIN_N_FRAGMENTS) {
				br.close();
				file=sequenceInfoFiles.readLine();
				continue;
			}
			type=tmpArray[3]; period=tmpArray[4];
			pathKernelLengths[0]=tmpArray[1];
			pathKernelPeriodic[0]=tmpArray[3]==Constants.INTERVAL_PERIODIC?(period<PeriodicSubstrings.MIN_PERIOD_LONG?(byte)1:(byte)2):(byte)0;
			
			// Building maximal intervals of overlapping/adjacent forbidden fragments
			lastForbidden=-1;
			if (!PRINT_REFERENCE_FRAGMENTS) {
				str=br.readLine();
				while (str!=null) {
					IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
					start=tmpArray[1]; end=tmpArray[2]; bidirectedGraphNode=tmpArray[6]; type=tmpArray[7]; period=tmpArray[8];
					if (!isForbidden(bidirectedGraphNode,start,end,type,period)) {
						str=br.readLine();
						continue;
					}
					read=tmpArray[0]; 
					if (lastForbidden==-1 || read!=forbidden[lastForbidden][0] || start>forbidden[lastForbidden][2]+IDENTITY_THRESHOLD) {
						lastForbidden++;
						if (lastForbidden>=forbidden.length) {
							int[][] newMatrix = new int[forbidden.length<<1][3];
							Math.copy(forbidden,newMatrix);
							forbidden=newMatrix;
						}
						forbidden[lastForbidden][0]=read;
						forbidden[lastForbidden][1]=start;
						forbidden[lastForbidden][2]=end;
					}
					else forbidden[lastForbidden][2]=Math.max(forbidden[lastForbidden][2],end);
					str=br.readLine();
				}
			}
			br.close();
			
			// Writing strings
			if (windows==null) {
				windows = new IntervalGraph.Node[nFragments];
				for (i=0; i<nFragments; i++) windows[i] = new IntervalGraph.Node();
			}
			else if (windows.length<nFragments) {
				IntervalGraph.Node[] newWindows = new IntervalGraph.Node[nFragments];
				System.arraycopy(windows,0,newWindows,0,windows.length);
				for (i=windows.length; i<nFragments; i++) newWindows[i] = new IntervalGraph.Node();
				windows=newWindows;
			}
			lastWindow=-1;
			fileID=file.substring(BASIN_LABEL_LENGTH,file.indexOf("."));
			bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+IO.FRAGMENTS_LABEL+"-"+fileID+".txt"));
			br = new BufferedReader(new FileReader(INPUT_DIR+"/"+file.trim()));
			currentRead=-1; nFragments=0;
			str=br.readLine(); str=br.readLine();  // Skipping header
			while (str!=null) {
				IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
				read=tmpArray[0]; start=tmpArray[1]; end=tmpArray[2];
				bidirectedGraphNode=tmpArray[6]; type=tmpArray[7]; period=tmpArray[8];
				otherKernels=tmpArray[11];
				if (ONE_KERNEL_ONLY && otherKernels==1) {
					str=br.readLine();
					continue;
				}
				if (!PRINT_REFERENCE_FRAGMENTS && isForbidden(bidirectedGraphNode,start,end,type,period)) {
					str=br.readLine();
					continue;
				}
				basinDescriptorRecord2Node(tmpNode,tmpArray,0);
				if (currentRead==-1) {
					currentRead=read; lastWindow=0;
					if (windows[0]==null) windows[0] = new IntervalGraph.Node();
					IntervalGraphStep3.printKernelTags_clone(tmpNode,windows[0]);
				}
				else if (read!=currentRead) {
					for (i=0; i<=lastWindow; i++) nFragments+=printQualityIntervals(windows[i].read,windows[i].start,windows[i].end,QUALITY_WINDOW,Reads.MAX_HIGH_QUALITY_SCORE,bw);
					currentRead=read; lastWindow=0;
					if (windows[0]==null) windows[0] = new IntervalGraph.Node();
					IntervalGraphStep3.printKernelTags_clone(tmpNode,windows[0]);
				}
				else {
					merged=false;
					for (i=lastWindow; i>=0; i--) {
						if (IntervalGraphStep3.printKernelTags_shouldBeMerged(tmpNode,windows[i],IDENTITY_THRESHOLD,pathKernelLengths,pathKernelPeriodic,false)) {
							IntervalGraphStep3.printKernelTags_merge(tmpNode,windows[i],false,IDENTITY_THRESHOLD);
							merged=true;
							// Merging with multiple windows is allowed
						}
					}
					if (!merged) IntervalGraphStep3.printKernelTags_clone(tmpNode,windows[++lastWindow]);
				}
				str=br.readLine();
			}
			for (i=0; i<=lastWindow; i++) nFragments+=printQualityIntervals(windows[i].read,windows[i].start,windows[i].end,QUALITY_WINDOW,Reads.MAX_HIGH_QUALITY_SCORE,bw);
			bw.close();
			if (nFragments<MIN_N_FRAGMENTS) {
				(new File(OUTPUT_DIR+"/"+IO.FRAGMENTS_LABEL+"-"+fileID+".txt")).delete();
				(new File(OUTPUT_DIR+"/"+IO.REFERENCE_LABEL+"-"+fileID+".txt")).delete();
			}
			file=sequenceInfoFiles.readLine();
		}
		sequenceInfoFiles.close();
	}
	
	
	/**
	 * @param node not NULL;
	 * @param record output of $IntervalGraphStep3.readBasinDescriptorRecord()$;
	 * @param kernel artificial kernel ID assigned to the node.
	 */
	private static final void basinDescriptorRecord2Node(IntervalGraph.Node node, int[] record, int kernel) {
		node.read=record[0]; node.start=record[1]; node.end=record[2];
		node.lastKernel=0;
		if (node.kernels==null) node.kernels = new int[1];
		node.kernels[0]=kernel;
		if (node.kernelOrientations==null) node.kernelOrientations = new byte[1];
		node.kernelOrientations[0]=(byte)record[3];
		final int atStart = record[4];
		final int atEnd = record[5];
		if (atStart==1) {
			node.lastPathWithStart=0;
			if (node.pathsWithStart==null) node.pathsWithStart = new int[1];
			node.pathsWithStart[0]=kernel;
		}
		else if (atStart==2) {
			node.lastPathWithStart=0;
			if (node.pathsWithStart==null) node.pathsWithStart = new int[1];
			node.pathsWithStart[0]=-1-kernel;
		}
		else node.lastPathWithStart=-1;
		if (atEnd==1) {
			node.lastPathWithEnd=0;
			if (node.pathsWithEnd==null) node.pathsWithEnd = new int[1];
			node.pathsWithEnd[0]=kernel;
		}
		else if (atEnd==2) {
			node.lastPathWithEnd=0;
			if (node.pathsWithEnd==null) node.pathsWithEnd = new int[1];
			node.pathsWithEnd[0]=-1-kernel;
		}
		else node.lastPathWithEnd=-1;
	}
	
	
	/**
	 * Prints to $bw$ all the maximal disjoint intervals of $read[start..end]$ such that 
	 * all the windows of size $window$ of an interval have average quality at most 
	 * $maxHighQualityScore$. Parts in $forbidden$ are not printed.
	 *
	 * @return number of intervals printed.
	 */
	private static final int printQualityIntervals(int read, int start, int end, int window, int maxHighQualityScore, BufferedWriter bw) throws IOException {
		final int UPPER_BOUND = (maxHighQualityScore*window)/Reads.QUALITY_SPACING;
		int i, j, w, s, e;
		int last, intervalStart, intervalEnd, out;
		double qualitySum;
		final double[] histogram = Reads.getQualityArray(read);
		
		Reads.readEnds2qualityEnds(read,start,end,true,tmpArray);
		s=tmpArray[0]; e=tmpArray[1]; w=window/Reads.QUALITY_SPACING;
		intervalStart=-1; intervalEnd=-1; i=s; out=0;
		qualitySum=0.0;
		for (j=0; j<w; j++) qualitySum+=histogram[i+j];
		while (i+w-1<=e) {
			if (qualitySum>UPPER_BOUND) {
				if (intervalStart!=-1) {
					last=subtractForbidden(read,Math.max(start,intervalStart*Reads.QUALITY_SPACING),Math.min(end,(intervalEnd+1)*Reads.QUALITY_SPACING-1),tmpArray);
					for (j=0; j<last; j+=2) {
						if (tmpArray[j+1]-tmpArray[j]+1>=window) {
							bw.write(read+","); bw.write(tmpArray[j]+","); bw.write(tmpArray[j+1]+"\n");
							out++;
						}
					}
					intervalStart=-1;
					i=intervalEnd+1;
					if (i+w-1<=e) {
						qualitySum=0.0;
						for (j=0; j<w; j++) qualitySum+=histogram[i+j];
					}
				}
				else {
					i++; j=i+w-1; 
					if (j<=e) qualitySum+=histogram[j]-histogram[i-1];
				}
				continue;
			}
			if (intervalStart==-1) intervalStart=i;
			intervalEnd=i+w-1;
			i++; j=i+w-1; 
			if (j<=e) qualitySum+=histogram[j]-histogram[i-1];
		}
		if (intervalStart!=-1) {
			last=subtractForbidden(read,Math.max(start,intervalStart*Reads.QUALITY_SPACING),Math.min(end,(intervalEnd+1)*Reads.QUALITY_SPACING-1),tmpArray);
			for (j=0; j<last; j+=2) {
				if (tmpArray[j+1]-tmpArray[j]+1>=window) {
					bw.write(read+","); bw.write(tmpArray[j]+","); bw.write(tmpArray[j+1]+"\n");
					out++;
				}
			}
		}
		return out;
	}
	
	
	private static final boolean isForbidden(int bidirectedGraphNode, int start, int end, int type, int period) {
		return bidirectedGraphNode>=0 && !(type==Constants.INTERVAL_PERIODIC && period>=PeriodicSubstrings.MIN_PERIOD_LONG && end-start+1>=period<<1);
	}
	
	
	/**
	 * Remark: the procedure assumes the intervals in $forbidden$ to be disjoint, and 
	 * $tmpArray$ to be large enough.
	 *
	 * @param tmpArray output array, containing the sorted list of intervals that results 
	 * from $read[start..end] \setminus forbidden$;
	 * @return the last element of $tmpArray$; -1 if the subtraction is empty.
	 */
	private static final int subtractForbidden(int read, int start, int end, int[] tmpArray) {
		int i;
		int last, intersectionStart, intersectionEnd;
		
		last=-1;
		for (i=0; i<=lastForbidden; i++) {
			if (forbidden[i][0]>read) break;
			if (forbidden[i][0]<read) continue;
			if (forbidden[i][1]>end) break;
			if (forbidden[i][2]<start) continue;
			intersectionStart=Math.max(forbidden[i][1],start);
			intersectionEnd=Math.min(forbidden[i][2],end);
			if (intersectionStart==start && intersectionEnd==end) return -1;
			if (intersectionStart>start) {
				if (last==-1) tmpArray[++last]=start;
				tmpArray[++last]=intersectionStart-1;
			}
			if (intersectionEnd<end) tmpArray[++last]=forbidden[i][2]+1;
		}
		if (last!=-1) {
			if ((last+1)%2!=0) tmpArray[++last]=end;
		} 
		else { tmpArray[0]=start; tmpArray[1]=end; last=1; }
		return last;
	}
	
	
}