package de.mpi_cbg.revant.fragments;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.CharStream;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.PeriodicSubstrings;
import de.mpi_cbg.revant.intervalgraph.BidirectedGraph;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;


/**
 * Prints the reference sequence of a path kernel described in the first line of a 
 * descriptor file built by $IntervalGraphStep3.printBasinDescriptors()$.
 * Descriptors containing too few fragments are not printed.
 * 
 * Remark: the procedure assumes that the $*-kernel*-labels-*.txt$ file, that contains
 * the labels of all nodes of the specified bidirected graph, has already been built.
 */
public class BasinDescriptor2Reference {
	
	/**
	 * @param args 
	 * 0: Step1 directory; 
	 * 1: min number of fragments for a descriptor file to be output at all;
	 * 2: min alignment length used during repeat inference.
	 */
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		final int MIN_N_FRAGMENTS = Integer.parseInt(args[1]);
		final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[2]);
		
		final String STEP3_DIR = STEP1_DIR+"/finalOutput";
		final String STEP4_DIR = STEP3_DIR+"/step4";
		final String STEP5_DIR = STEP4_DIR+"/step5";
		final String BASIN_DESCRIPTOR_FILES = STEP5_DIR+"/descriptors.txt";
		final String INPUT_DIR = STEP4_DIR;
		final String INPUT_DIR_BIDIRECTED = STEP3_DIR;  // We assume that this directory containing alternative descriptor files whose $K,P$ fields of the header should be used (discarded if null).
		final String OUTPUT_DIR = STEP5_DIR+"/references-strings";
		final String READS_DIR = STEP5_DIR+"/reads";  // We assume that this directory contains the read strings files of every read in a basin descriptor.
		
		final int BASIN_LABEL_LENGTH = IO.BASIN_DESCRIPTOR_PREFIX.length()+1;
		final boolean SAME_INPUT_DIR = INPUT_DIR.equalsIgnoreCase(INPUT_DIR_BIDIRECTED);
		boolean found, hasLongPeriod;
		int i, p, q;
		int length, pathLength, pathID, type, period, maxRead, maxStart, maxEnd, nFragments;
		int kernelID, header, lastTuple;
		String file, fileID, str, componentID, clusterID, pathKernelID, kernelStrings, readStringFile;
		FileReader fr;
		BufferedReader br, brPrime, sequenceInfoFiles;
		CharStream stream;
		BufferedWriter bw;
		int[] out = new int[5];
		int[] tmpArray = new int[20];
		StatsTuple[] tuples = new StatsTuple[100];  // Arbitrary capacity
		
		Alignments.minAlignmentLength=MIN_ALIGNMENT_LENGTH;
		header=Math.POSITIVE_INFINITY>>1;
		sequenceInfoFiles = new BufferedReader(new FileReader(BASIN_DESCRIPTOR_FILES));
		file=sequenceInfoFiles.readLine(); lastTuple=-1;
		while (file!=null) {		
			p=BASIN_LABEL_LENGTH;
			q=file.indexOf("-",p);
			componentID=file.substring(p,q);
			p=q; q=file.indexOf("-",p+1);
			clusterID=file.substring(p+1,q);
			p=q; q=file.indexOf(".",p+1);
			pathKernelID=file.substring(p+1,q);
			fileID=file.substring(BASIN_LABEL_LENGTH,q);
			br = new BufferedReader(new FileReader(INPUT_DIR+"/"+file.trim()));
			str=br.readLine();
			IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
			kernelID=tmpArray[0]; pathLength=tmpArray[1]; pathID=tmpArray[2];
			type=tmpArray[3]; period=tmpArray[4]; nFragments=tmpArray[5];
			if (nFragments<MIN_N_FRAGMENTS) {
				br.close();
				file=sequenceInfoFiles.readLine();
				continue;
			}
			if (!SAME_INPUT_DIR) {
				brPrime = new BufferedReader(new FileReader(INPUT_DIR_BIDIRECTED+"/"+file.trim()));
				str=brPrime.readLine();
				brPrime.close();
				IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
				kernelID=tmpArray[0]; pathID=tmpArray[2];
			}
			System.out.println("Printing reference of file "+file);
			lastTuple++;
			if (lastTuple==tuples.length) {
				StatsTuple[] newTuples = new StatsTuple[tuples.length<<1];
				System.arraycopy(tuples,0,newTuples,0,tuples.length);
				tuples=newTuples;
			}
			tuples[lastTuple] = new StatsTuple(type,pathLength,period,fileID);
			hasLongPeriod=type==Constants.INTERVAL_PERIODIC && period>=PeriodicSubstrings.MIN_PERIOD_LONG;
			System.out.print("   Loading bidirected graph... ");
			BidirectedGraph.deserialize(STEP1_DIR+"/"+componentID+"-clusters/"+clusterID+"-kernel"+kernelID+(pathLength==0?"-"+IO.CYCLIC_LABEL:"")+".bdgraph");
			System.out.println("done. nNodes="+BidirectedGraph.nNodes+" lastAssemblyPath="+BidirectedGraph.lastAssemblyPath+" assemblyPaths.length="+BidirectedGraph.assemblyPaths.length);
			BidirectedGraph.print();
			if (pathLength==0) {
				System.out.print("   Printing a representative of the cyclic graph... ");
				if (type==Constants.INTERVAL_PERIODIC) getRepresentativeInterval_periodic(period,br,out,tmpArray);
				else getRepresentativeInterval_nonperiodic(br,out);
				br.close();
				maxRead=out[0]; maxStart=out[1]; maxEnd=out[2];
				readStringFile=READS_DIR+"/read"+(maxRead+1)+".txt";
				br = new BufferedReader(new FileReader(readStringFile));
				str=br.readLine();
				br.close();
				bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+IO.REFERENCE_LABEL+"-"+fileID+".txt"),IO.BUFFER_SIZE);
				IO.writeFakeHeader(header++/*Unlikely to collide with any meaningful header*/,maxEnd-maxStart+1,null,bw);
				bw.write(str.substring(maxStart,maxEnd+1));
				bw.write('\n');
				bw.close();
				str=null;
			}
			else {
				br.close();
				System.out.println("   Printing the "+pathID+"-th path... period="+period+" hasLongPeriod="+hasLongPeriod);
				found=false; kernelStrings=null;
				try { 
					kernelStrings=STEP1_DIR+"/"+componentID+"-clusters/"+clusterID+"-kernel"+kernelID+"-"+IO.ONE_NODE_ONLY_LABEL+"-"+IO.LABELS_LABEL+".txt";
					fr = new FileReader(kernelStrings); 
					found=true;
				}
				catch (FileNotFoundException e) { }
				if (!found) {
					try {
						kernelStrings=STEP1_DIR+"/"+componentID+"-clusters/"+clusterID+"-kernel"+kernelID+"-"+IO.LABELS_LABEL+".txt";
						fr = new FileReader(kernelStrings);
						found=true;
					}
					catch (FileNotFoundException e) { }
				}
				if (!found) {
					try {
						kernelStrings=STEP1_DIR+"/"+componentID+"-clusters/"+clusterID+"-kernel"+kernelID+"-"+IO.TOO_MANY_PATHS_LABEL+"-"+IO.LABELS_LABEL+".txt";
						fr = new FileReader(kernelStrings);
						found=true;
					}
					catch (FileNotFoundException e) { }
				}
				if (!found) {
					System.err.println("ERROR> labels not found for descriptor "+file);
					System.exit(1);
				}
				stream = new CharStream(4096);  // Arbitrary capacity
				length=BidirectedGraph.printPath(pathID,null,kernelStrings,OUTPUT_DIR+"/"+IO.REFERENCE_LABEL+"-"+fileID+".txt",stream,(hasLongPeriod&&(period<pathLength))?period:0,header++/*Unlikely to collide with any meaningful header*/);
				stream.deallocate(); stream=null;
				if (length!=pathLength && !hasLongPeriod) {
					System.err.println("ERROR: expectedLength="+pathLength+" != actualLength="+length);
					System.exit(1);
				}
			}
			file=sequenceInfoFiles.readLine();
		}
		
		// Printing stats
		if (lastTuple>0) Arrays.sort(tuples,0,lastTuple+1);
		System.err.println("Statistics on basin descriptor files:  [length,period,type,id]");
		for (i=0; i<=lastTuple; i++) System.err.println(tuples[i].length+","+tuples[i].period+","+tuples[i].type+","+tuples[i].id);
	}
	
	
	/**
	 * @param br assumed to be positioned on the first line of the basin descriptor file;
	 * @param out output array: 0=read ID, 1=absolute start, 2=absolute end of a longest 
	 * interval.
	 */
	private static final void getRepresentativeInterval_nonperiodic(BufferedReader br, int[] out) throws IOException {
		int i, p, q;
		int read, start, end, length, maxLength, maxRead, maxStart, maxEnd;
		String str;
		
		maxRead=-1; maxStart=-1; maxEnd=-1; maxLength=0;
		str=br.readLine();
		while (str!=null) {
			p=str.indexOf(","); 
			read=Integer.parseInt(str.substring(0,p));
			q=str.indexOf(",",p+1);
			start=Integer.parseInt(str.substring(p+1,q));
			p=q; q=str.indexOf(",",p+1);
			end=Integer.parseInt(str.substring(p+1,q));
			length=end-start+1;
			if (length>maxLength) {
				maxLength=length; maxRead=read; maxStart=start; maxEnd=end;
			}
			str=br.readLine();
		}
		out[0]=maxRead; out[1]=maxStart; out[2]=maxEnd;
	}
	
	
	/**
	 * @param br assumed to be positioned on the first line of the basin descriptor file;
	 * @param out output array: 0=read ID; 1=absolute start; 2=absolute end; if 
	 * $period$ is short, the returned substring is a longest short-period interval in the
	 * basin; otherwise, it is a unit of a long-period interval in the basin. The 
	 * procedure tries to take such a unit from an interval whose start or end is adjacent 
	 * to a high-quality region of a read, since it is more likely that the resulting 
	 * substring is a full copy of the period. If this is not possible, the procedure 
	 * returns a random substring of length equal to the period, and this might be a 
	 * rotation of the actual period.
	 *
	 * Remark: alternatively, we could return a longest merge of overlapping intervals in 
	 * the same read, as done in $BasinDescriptor2Fragments$, or one of the intervals 
	 * produced by $BasinDescriptor2Fragments$. We don't do that to keep the code simple.
	 */
	private static final void getRepresentativeInterval_periodic(int period, BufferedReader br, int[] out, int[] tmpArray) throws IOException {
		boolean ms, me;
		int i, p, q;
		int prd, type, read, start, end, length, maxLength, maxRead, maxStart, maxEnd;
		String str;
		
		maxRead=-1; maxStart=-1; maxEnd=-1; maxLength=0;
		str=br.readLine();
		while (str!=null) {
			IntervalGraphStep3.readBasinDescriptorRecord(str,tmpArray);
			read=tmpArray[0]; start=tmpArray[1]; end=tmpArray[2];
			type=tmpArray[7]; prd=tmpArray[8];
			if (type==-1) {
				System.err.println("ERROR: line with no type?! "+str);
				System.exit(1);
			}
			if (type!=Constants.INTERVAL_PERIODIC || (prd<PeriodicSubstrings.MIN_PERIOD_LONG && period>=PeriodicSubstrings.MIN_PERIOD_LONG) || (prd>=PeriodicSubstrings.MIN_PERIOD_LONG && period<PeriodicSubstrings.MIN_PERIOD_LONG)) {
				str=br.readLine();
				continue;
			}
			length=end-start+1; ms=tmpArray[9]==1; me=tmpArray[10]==1;
			if (period>=PeriodicSubstrings.MIN_PERIOD_LONG) {
				out[0]=read;
				if (ms || !me) {
					out[1]=start;
					out[2]=start+period-1;
					return;
				}
				else if (me) {
					out[1]=end-period+1;
					out[2]=end;
					return;
				}
				else if (length>maxLength) {
					maxLength=length; maxRead=read; maxStart=start; maxEnd=end;
				}
			}
			else {
				if (length>maxLength) {
					maxLength=length; maxRead=read; maxStart=start; maxEnd=end;
				}
			}
			str=br.readLine();
		}
		out[0]=maxRead; out[1]=maxStart;
		out[2]=period>=PeriodicSubstrings.MIN_PERIOD_LONG?maxStart+period-1:maxEnd;
	}
	
	
	private static class StatsTuple implements Comparable {
		public int type, length, period;
		public String id;
		
		public StatsTuple(int t, int l, int p, String d) {
			type=t; length=l; period=p; id=d;
		}
		
		public int compareTo(Object other) {
			StatsTuple otherTuple = (StatsTuple)other;
			if (length<otherTuple.length) return 1;
			else if (length>otherTuple.length) return -1;
			if (period<otherTuple.period) return 1;
			else if (period>otherTuple.period) return -1;
			return 0;
		}
	}
	
}