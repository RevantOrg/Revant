package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class IntervalGraphStep2_0_getParallelWork {

	
	/**
	 * args[0]: Step 1 dir.
	 * args[1]: nThreads.
	 */
	public static void main(String[] args) throws IOException {
		String FOLDER = args[0];
		int nThreads = Integer.parseInt(args[1]);
		
		final String SUFFIX = ".graph";
		int i, j, lastFile, currentThread, filesPerThread;
		long workPerThread, currentWork;
		File folder, file;
		BufferedWriter bw;
		long[] lengths;
		String[] paths;
		
		// Loading graph files and lengths
		folder = new File(FOLDER);
		paths=folder.list();
		if (paths.length>1) Arrays.sort(paths);
		lengths = new long[paths.length];
		for (i=0; i<paths.length; i++) {
			file = new File(FOLDER+"/"+paths[i]);
			lengths[i]=file.length();
		}
		lastFile=-1; workPerThread=0;
		for (i=0; i<paths.length; i++) {
			file = new File(FOLDER+"/"+paths[i]);
			if (!file.isFile() || !paths[i].substring(paths[i].length()-SUFFIX.length(),paths[i].length()).equals(SUFFIX)) continue;
			lastFile++;
			paths[lastFile]=FOLDER+"/"+paths[i];
			lengths[lastFile]=lengths[i];
			workPerThread+=lengths[i];
		}
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<=lastFile; i++) {
				file = new File(paths[i]);
				if (!file.exists()) {
					System.err.println("IntervalGraphStep2_0_getParallelWork> ERROR: the following file does not exist: "+paths[i]);
					System.exit(1);
				}
			}
		}
		workPerThread/=nThreads;
		System.out.println("workPerThread="+workPerThread);
		
		// Splitting work into approx. equal parts (using file size as a proxy for
		// work).
/*		currentThread=0; currentWork=0L;
		bw = new BufferedWriter(new FileWriter(FOLDER+"/step2_0-thread"+currentThread+".txt"));
		for (i=0; i<=lastFile; i++) {
			bw.write(paths[i]+"\n");
			currentWork+=lengths[i];
			if (currentWork>=workPerThread) {
				bw.close();
				currentThread++; currentWork=0L;
				bw = new BufferedWriter(new FileWriter(FOLDER+"/step2_0-thread"+currentThread+".txt"));
			}
		}
		bw.close();
*/		
		
		// Splitting work into approx. equal parts (using file number as a proxy for
		// work).
		filesPerThread=(lastFile+1)/nThreads;
		for (i=0; i<nThreads; i++) {
			bw = new BufferedWriter(new FileWriter(FOLDER+"/step2_0-thread"+i+".txt"));
			for (j=i*filesPerThread; j<Math.min((i+1)*filesPerThread,lastFile+1); j++) bw.write(paths[j]+"\n");
			bw.close();
		}
		
	}

}