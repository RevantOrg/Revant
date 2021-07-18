package de.mpi_cbg.revant.intervalgraph;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;


public class IntervalGraphStep3_getParallelWork {
	
	public final static String SUFFIX_STEP20 = "-step2_0";
	public final static String SUFFIX_STEP2 = "-step2";
	public final static String SUFFIX_STEP25 = "-step2_5";
	public final static String SUFFIX_GRAPH = ".graph";
	

	public static void main(String[] args) throws IOException {
		String ROOT_FOLDER = args[0];
		final int N_THREADS = Integer.parseInt(args[1]);

		boolean found20, found2, found25;
		int i, j, k, h;
		int lastFile, currentThread, filesPerThread;
		long workPerThread, currentWork;
		String folderPath, step2FolderPath, step25FolderPath, step25FolderPath_nested;
		String step1FilePath, step20FilePath, step2FilePath, step25FilePath, step25FilePath_nested;
		File folder, step1Folder, step2Folder, step25Folder;
		File step1File, step20File, step2File, step25File;
		File file;
		BufferedWriter bw;
		long[] lengths;
		String[] step1Paths, step20Paths, step2Paths, step25Paths;
		String[] paths, pathsPrime;
		
		step1Folder = new File(ROOT_FOLDER);
		step1Paths=step1Folder.list();
		paths = new String[step1Paths.length];
		lastFile=-1;
		for (i=0; i<step1Paths.length; i++) {
			folderPath=ROOT_FOLDER+"/"+step1Paths[i];
			folder = new File(folderPath);
			if (!folder.isDirectory()) continue;
			if (step1Paths[i].length()>=SUFFIX_STEP20.length() && step1Paths[i].substring(step1Paths[i].length()-SUFFIX_STEP20.length(),step1Paths[i].length()).equals(SUFFIX_STEP20)) {
				// Entering the Step-2.0 folder
				step20Paths=folder.list();
				found20=false;
				for (j=0; j<step20Paths.length; j++) {
					step20FilePath=folderPath+"/"+step20Paths[j];
					step20File = new File(step20FilePath);
					if (!step20File.isFile() || step20FilePath.length()<=SUFFIX_GRAPH.length() || !step20FilePath.substring(step20FilePath.length()-SUFFIX_GRAPH.length(),step20FilePath.length()).equals(SUFFIX_GRAPH)) continue;
					step2FolderPath=folderPath+"/"+step20Paths[j].substring(0,step20Paths[j].length()-SUFFIX_GRAPH.length())+SUFFIX_STEP2;
					step2Folder = new File(step2FolderPath);
				
					// This Step-2.0 graph doesn't have a Step-2 folder: using the
					// Step-2.0 graph.
					if (!step2Folder.exists() || !step2Folder.isDirectory()) {
						lastFile++;
						if (lastFile>=paths.length) {
							pathsPrime = new String[paths.length<<1];
							System.arraycopy(paths,0,pathsPrime,0,paths.length);
							paths=pathsPrime;
						}
						paths[lastFile]=step20FilePath;
						found20=true;
						continue;
					}
					else {
						// Entering the Step-2 folder
						found2=false;
						step2Paths=step2Folder.list();
						for (k=0; k<step2Paths.length; k++) {
							step2FilePath=step2FolderPath+"/"+step2Paths[k];
							step2File = new File(step2FilePath);
							if (!step2File.isFile() || step2FilePath.length()<=SUFFIX_GRAPH.length() || !step2FilePath.substring(step2FilePath.length()-SUFFIX_GRAPH.length(),step2FilePath.length()).equals(SUFFIX_GRAPH)) continue;
							step25FolderPath=step2FolderPath+"/"+step2Paths[k].substring(0,step2Paths[k].length()-SUFFIX_GRAPH.length())+SUFFIX_STEP25;
							step25Folder = new File(step25FolderPath);
			
							// This Step-2 graph has no Step-2.5 folder: using the
							// Step-2 graph.
							if (!step25Folder.exists() || !step25Folder.isDirectory()) {
								lastFile++;
								if (lastFile>=paths.length) {
									pathsPrime = new String[paths.length<<1];
									System.arraycopy(paths,0,pathsPrime,0,paths.length);
									paths=pathsPrime;
								}
								paths[lastFile]=step2FilePath;
								found2=true; found20=true;
								continue;
							}
			
							// Entering the Step-2.5 folder
							found25=false;
							step25Paths=step25Folder.list();
							for (h=0; h<step25Paths.length; h++) {
								step25FilePath=step25FolderPath+"/"+step25Paths[h];
								step25File = new File(step25FilePath);
								if (!step25File.isFile() || step25FilePath.length()<=SUFFIX_GRAPH.length() || !step25FilePath.substring(step25FilePath.length()-SUFFIX_GRAPH.length(),step25FilePath.length()).equals(SUFFIX_GRAPH)) continue;
								lastFile++;
								if (lastFile>=paths.length) {
									pathsPrime = new String[paths.length<<1];
									System.arraycopy(paths,0,pathsPrime,0,paths.length);
									paths=pathsPrime;
								}
								paths[lastFile]=step25FilePath;
								found25=true; found2=true; found20=true;
							}
							if (!found25) {
								// The Step-2.5 folder is empty: using the Step-2
								// graph.
								lastFile++;
								if (lastFile>=paths.length) {
									pathsPrime = new String[paths.length<<1];
									System.arraycopy(paths,0,pathsPrime,0,paths.length);
									paths=pathsPrime;
								}
								paths[lastFile]=step2FilePath;
								found2=true; found20=true;
							}
						}
						if (!found2) {
							// The Step-2 folder is empty: using the Step-2.0 graph.
							lastFile++;
							if (lastFile>=paths.length) {
								pathsPrime = new String[paths.length<<1];
								System.arraycopy(paths,0,pathsPrime,0,paths.length);
								paths=pathsPrime;
							}
							paths[lastFile]=step20FilePath;
							found20=true;
						}
					}
				}
				if (!found20) {
					// The Step-2.0 folder is empty: using the Step-1 graph.
					lastFile++;
					if (lastFile>=paths.length) {
						pathsPrime = new String[paths.length<<1];
						System.arraycopy(paths,0,pathsPrime,0,paths.length);
						paths=pathsPrime;
					}
					paths[lastFile]=ROOT_FOLDER+"/"+step1Paths[i].substring(0,step1Paths[i].indexOf(SUFFIX_STEP20))+SUFFIX_GRAPH;
				}
			}
		}
		if (IO.CONSISTENCY_CHECKS) {
			for (i=0; i<=lastFile; i++) {
				file = new File(paths[i]);
				if (!file.exists()) {
					System.err.println("IntervalGraphStep3_getParallelWork> ERROR: the following file does not exist: "+paths[i]);
					System.exit(1);
				}
			}
		}
		if (lastFile>0) Arrays.sort(paths,0,lastFile+1);
		
		// Adding file sizes
		lengths = new long[lastFile+1];
		workPerThread=0L;
		for (i=0; i<=lastFile; i++) {
			file = new File(paths[i]);
			lengths[i]=file.length();
			workPerThread+=lengths[i];
		}
		workPerThread/=N_THREADS;
		
		// Splitting work into almost equal parts
/*		currentThread=0; currentWork=0L;
		bw = new BufferedWriter(new FileWriter(ROOT_FOLDER+"/step3-thread"+currentThread+".txt"));
		for (i=0; i<=lastFile; i++) {
			bw.write(paths[i]+"\n");
			currentWork+=lengths[i];
			if (currentWork>=workPerThread) {
				bw.close();
				currentThread++; currentWork=0L;
				bw = new BufferedWriter(new FileWriter(ROOT_FOLDER+"/step3-thread"+currentThread+".txt"));
			}
		}
		bw.close();
*/
		
		// Splitting work into approx. equal parts (using file number as a proxy for
		// work).
		filesPerThread=(lastFile+1)/N_THREADS;
		for (i=0; i<N_THREADS; i++) {
			bw = new BufferedWriter(new FileWriter(ROOT_FOLDER+"/step3-thread"+i+".txt"));
			for (j=i*filesPerThread; j<Math.min((i+1)*filesPerThread,lastFile+1); j++) bw.write(paths[j]+"\n");
			bw.close();
		}		
		
		
	}
	
	
	/**
	 * Converts a path of a $.graph$ file in the file system hierarchy, into a unique ID
	 * that keeps only the numbers of each step.
	 *
	 * @param lastCharacterInPrefix last position of a prefix of $path$ that should not
	 * be used by the procedure.
	 */
	public static final String getFileID(String path, int lastCharacterInPrefix) {
		final char SEPARATOR = '/';
		
		int i, p, q;
		int last;
		String suffix, out;
		
		out=""; suffix=null;
		last=path.length()-1;
		while (last>lastCharacterInPrefix) {
			if (path.substring(last-SUFFIX_GRAPH.length()+1,last+1).equalsIgnoreCase(SUFFIX_GRAPH)) suffix=SUFFIX_GRAPH;
			else if (path.substring(last-SUFFIX_STEP25.length()+1,last+1).equalsIgnoreCase(SUFFIX_STEP25)) suffix=SUFFIX_STEP25;
			else if (path.substring(last-SUFFIX_STEP2.length()+1,last+1).equalsIgnoreCase(SUFFIX_STEP2)) suffix=SUFFIX_STEP2;
			else if (path.substring(last-SUFFIX_STEP20.length()+1,last+1).equalsIgnoreCase(SUFFIX_STEP20)) suffix=SUFFIX_STEP20;
			else {
				System.err.println("getFileID> ERROR");
				System.exit(1);
			}
			q=last-suffix.length();
			p=-1;
			for (i=q; i>=0; i--) {
				if (path.charAt(i)==SEPARATOR) {
					p=i+1;
					break;
				}
			}
			if (out.length()==0) out=path.substring(p,q+1);
			else out=path.substring(p,q+1)+"-"+out;
			last=p-2;
		}
		return out;
	}
	
	
}