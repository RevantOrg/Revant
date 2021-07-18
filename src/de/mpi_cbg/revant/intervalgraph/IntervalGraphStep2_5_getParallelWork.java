package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class IntervalGraphStep2_5_getParallelWork {


	public static void main(String[] args) throws IOException {
		String ROOT_FOLDER = args[0];
		final int N_THREADS = Integer.parseInt(args[1]);
		
		final String SUFFIX_GRAPH = ".graph";
		final String SUFFIX_STEP2_0 = "-step2_0";
		final String SUFFIX_STEP2 = "-step2";
		int i, j, k;
		int lastFile, currentThread, filesPerThread;
		long workPerThread, currentWork;
		String step1path, step20path, step2path;
		File file, step1dir, step1file, step20file, step20dir, step2file, toSplitFurtherFile;
		BufferedWriter bw;
		long[] lengths;
		String[] step1paths, step20paths, step2paths, toSplitFurtherIDs, files, filesPrime;
		
		// Loading graph files
		step1dir = new File(ROOT_FOLDER);
		step1paths=step1dir.list();
		files = new String[step1paths.length];
		lastFile=-1;
		for (i=0; i<step1paths.length; i++) {
			step1path=ROOT_FOLDER+"/"+step1paths[i];
			step1dir = new File(step1path);
			if (!step1dir.isDirectory() || step1paths[i].length()<=SUFFIX_STEP2_0.length() || !step1paths[i].substring(step1paths[i].length()-SUFFIX_STEP2_0.length(),step1paths[i].length()).equals(SUFFIX_STEP2_0)) continue;
			step20paths=step1dir.list();
			for (j=0; j<step20paths.length; j++) {
				step20path=step1path+"/"+step20paths[j];
				step20dir = new File(step20path);
				if (!step20dir.isDirectory() || step20paths[j].length()<=SUFFIX_STEP2.length() || !step20paths[j].substring(step20paths[j].length()-SUFFIX_STEP2.length(),step20paths[j].length()).equals(SUFFIX_STEP2)) continue;			
				// Clusters from large components produced by Step 2
				toSplitFurtherFile = new File(step20path+"/toSplitFurther.txt");
				if (toSplitFurtherFile.length()==0L) continue;
				toSplitFurtherIDs=readToSplitFurtherFile(step20path+"/toSplitFurther.txt");
				for (k=0; k<toSplitFurtherIDs.length; k++) {
					lastFile++;
					if (lastFile>=files.length) {
						filesPrime = new String[files.length<<1];
						System.arraycopy(files,0,filesPrime,0,files.length);
						files=filesPrime;
					}
					files[lastFile]=step20path+"/"+toSplitFurtherIDs[k]+SUFFIX_GRAPH;
				}
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<=lastFile; i++) {
				file = new File(files[i]);
				if (!file.exists()) {
					System.err.println("IntervalGraphStep2_5_getParallelWork> ERROR: the following file does not exist: "+files[i]);
					System.exit(1);
				}
			}
		}
		if (lastFile>0) Arrays.sort(files,0,lastFile+1);
		
		// Adding file sizes
		lengths = new long[lastFile+1];
		workPerThread=0L;
		for (i=0; i<=lastFile; i++) {
			file = new File(files[i]);
			lengths[i]=file.length();
			workPerThread+=lengths[i];
		}
		workPerThread/=N_THREADS;
		
		// Splitting work into almost equal parts
/*		currentThread=0; currentWork=0L;
		bw = new BufferedWriter(new FileWriter(ROOT_FOLDER+"/step2_5-thread"+currentThread+".txt"));
		for (i=0; i<=lastFile; i++) {
			bw.write(files[i]+"\n");
			currentWork+=lengths[i];
			if (currentWork>=workPerThread) {
				bw.close();
				currentThread++; currentWork=0L;
				bw = new BufferedWriter(new FileWriter(ROOT_FOLDER+"/step2_5-thread"+currentThread+".txt"));
			}
		}
		bw.close();
*/		
		
		// Splitting work into approx. equal parts (using file number as a proxy for
		// work).
		filesPerThread=(lastFile+1)/N_THREADS;
		for (i=0; i<N_THREADS; i++) {
			bw = new BufferedWriter(new FileWriter(ROOT_FOLDER+"/step2_5-thread"+i+".txt"));
			for (j=i*filesPerThread; j<Math.min((i+1)*filesPerThread,lastFile+1); j++) bw.write(files[j]+"\n");
			bw.close();
		}
		
		
	}
	
	
	public static final String[] readToSplitFurtherFile(String path) throws IOException {
		int i, count;
		String str;
		BufferedReader br;
		String[] out;
		
		br = new BufferedReader(new FileReader(path));
		count=0; str=br.readLine();
		while (str!=null) {
			count++;
			str=br.readLine();
		}
		br.close();
		out = new String[count];
		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); i=-1;
		while (str!=null) {
			out[++i]=str;
			str=br.readLine();
		}
		br.close();
		return out;
	}
	
}