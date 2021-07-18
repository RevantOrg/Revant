package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class IntervalGraphStep2_getParallelWork {

	
	/**
	 * args[0]: Step 1 dir.
	 * args[1]: nThreads.
	 */
	public static void main(String[] args) throws IOException {
		String FOLDER = args[0];
		int nThreads = Integer.parseInt(args[1]);
		
		final String SUFFIX = ".graph";
		final String SUFFIX_MPS = ".minPeelingSize";
		final String STEP_2_0_SUFFIX = "-step2_0";
		int i, j, nFiles, lastFile, currentThread, mps, filesPerThread;
		long workPerThread, currentWork;
		File folder, subfolder, file;
		BufferedReader br;
		BufferedWriter bw;
		int[] minPeelingSize;
		long[] lengths;
		String[] pathsPrime, paths, subpaths;
		
		// Counting the number of files
		nFiles=0;
		folder = new File(FOLDER);
		pathsPrime=folder.list();
		if (pathsPrime.length>1) Arrays.sort(pathsPrime);
		for (i=0; i<pathsPrime.length; i++) {
			subfolder = new File(FOLDER+"/"+pathsPrime[i]);
			if (!subfolder.isDirectory()) continue;			
			subpaths=subfolder.list();
			for (j=0; j<subpaths.length; j++) {
				file = new File(FOLDER+"/"+pathsPrime[i]+"/"+subpaths[j]);
				if (!file.isFile() || subpaths[j].length()<SUFFIX.length() || !subpaths[j].substring(subpaths[j].length()-SUFFIX.length(),subpaths[j].length()).equals(SUFFIX)) continue;
				nFiles++;
			}
		}
		paths = new String[nFiles]; lengths = new long[nFiles]; minPeelingSize = new int[nFiles];
		
		// Loading file lengths
		lastFile=-1; workPerThread=0;
		for (i=0; i<pathsPrime.length; i++) {
			subfolder = new File(FOLDER+"/"+pathsPrime[i]);
			if (!subfolder.isDirectory()) continue;
			file = new File(FOLDER+"/"+pathsPrime[i].substring(0,pathsPrime[i].indexOf(STEP_2_0_SUFFIX))+SUFFIX_MPS);
			br = new BufferedReader(new FileReader(file.getAbsolutePath()));
			mps=Integer.parseInt(br.readLine());
			br.close();
			subpaths=subfolder.list();
			if (subpaths.length>1) Arrays.sort(subpaths);
			for (j=0; j<subpaths.length; j++) {
				file = new File(FOLDER+"/"+pathsPrime[i]+"/"+subpaths[j]);
				if (!file.isFile() || subpaths[j].length()<SUFFIX.length() || !subpaths[j].substring(subpaths[j].length()-SUFFIX.length(),subpaths[j].length()).equals(SUFFIX)) continue;
				lastFile++;
				paths[lastFile]=FOLDER+"/"+pathsPrime[i]+"/"+subpaths[j];
				lengths[lastFile]=file.length();
				minPeelingSize[lastFile]=mps;
				workPerThread+=lengths[lastFile];
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<=lastFile; i++) {
				file = new File(paths[i]);
				if (!file.exists()) {
					System.err.println("IntervalGraphStep2_getParallelWork> ERROR: the following file does not exist: "+paths[i]);
					System.exit(1);
				}
			}
		}
		workPerThread/=nThreads;
		System.out.println("workPerThread="+workPerThread);
		
		// Splitting work into approx. equal parts (using file size as a proxy for
		// work).
/*		currentThread=0; currentWork=0L;
		bw = new BufferedWriter(new FileWriter(FOLDER+"/step2-thread"+currentThread+".txt"));
		for (i=0; i<nFiles; i++) {
			bw.write(paths[i]+"!"+minPeelingSize[i]+"\n");
			currentWork+=lengths[i];
			if (currentWork>=workPerThread) {
				bw.close();
				currentThread++; currentWork=0L;
				bw = new BufferedWriter(new FileWriter(FOLDER+"/step2-thread"+currentThread+".txt"));
			}
		}
		bw.close();
*/		
		
		// Splitting work into approx. equal parts (using file number as a proxy for
		// work).
		filesPerThread=(lastFile+1)/nThreads;
		for (i=0; i<nThreads; i++) {
			bw = new BufferedWriter(new FileWriter(FOLDER+"/step2-thread"+i+".txt"));
			for (j=i*filesPerThread; j<Math.min((i+1)*filesPerThread,lastFile+1); j++) bw.write(paths[j]+"!"+minPeelingSize[j]+"\n");
			bw.close();
		}
		
	}

}