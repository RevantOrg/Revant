package de.mpi_cbg.revant.consensus;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Alignments;


/**
 * Simulates how Pippel is building the repeat track, i.e. it just creates a track by 
 * merging all the readA intervals of all block-consensi alignment files.
 */
public class PippelBuildRepeatTrackOfDatabase {

	/**
	 * @param args 0=directory containing one subdirectory per module; every subdirectory
	 * contains a block-consensus alignments file. 
	 */
	public static void main(String[] args) throws IOException {
		final String ROOT_DIR = args[0];
		final int FIRST_READ = Integer.parseInt(args[1]);
		final int LAST_READ = Integer.parseInt(args[2]);
		final double MAX_ERROR = Double.parseDouble(args[3]);
		
		final int N_READS = LAST_READ-FIRST_READ+1;
		final String PREFIX = IO.MODULE_ALIGNMENTS_DIR_PREFIX;
		
		int i, j;
		int last;
		File file;
		int[] lastRepeatTrack;
		String[] list;
		int[][] repeatTrack;
		
		// Loading and merging tracks
		repeatTrack = new int[N_READS][0];
		lastRepeatTrack = new int[N_READS];
		Math.set(lastRepeatTrack,N_READS-1,-1);
		file = new File(ROOT_DIR);
		list=file.list();
		for (i=0; i<list.length; i++) {
			if (!list[i].startsWith(PREFIX)) continue;
			System.err.println("Loading "+list[i]+"...");
			alignments2repeatTrack(ROOT_DIR+"/"+list[i]+"/LAshow-block-consensus.txt",MAX_ERROR,FIRST_READ,repeatTrack,lastRepeatTrack);
		}
		ConsensusStep1.compactTracks(N_READS,FIRST_READ,repeatTrack,lastRepeatTrack,false);
		
		// Outputting
		for (i=0; i<N_READS; i++) {
			last=lastRepeatTrack[i];
			if (last==-1) continue;
			System.out.print((FIRST_READ+i)+" ");
			for (j=0; j<last; j+=2) System.out.print(repeatTrack[i][j]+" "+repeatTrack[i][j+1]+" ");
			System.out.println();
		}
	}
	
	
	/**
	 * Appends the contents of $path$ to $repeatTrack,lastRepeatTrack$.
	 *
	 * @param path assumed to contain block-consensus alignments, i.e. readA is a read ID
	 * and readB is a consensus ID;
	 * @param maxError only alignments with error rate at most this are used.
	 */
	private static final void alignments2repeatTrack(String path, double maxError, int firstRead, int[][] repeatTrack, int[] lastRepeatTrack) throws IOException {
		int last, length, read, start, end, currentRead, currentStart, currentEnd;
		double error;
		String str;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); str=br.readLine();  // Skipping header
		str=br.readLine();
		currentRead=-1; currentStart=-1; currentEnd=-1;
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			error=((double)(Alignments.diffs<<1))/(Alignments.endA-Alignments.startA+Alignments.endB-Alignments.startB+2);
			if (error>maxError) {
				str=br.readLine();
				continue;
			}
			read=Alignments.readA-1-firstRead; start=Alignments.startA; end=Alignments.endA;
			if (read!=currentRead) {
				if (currentRead!=-1) {
					last=lastRepeatTrack[currentRead];
					length=repeatTrack[currentRead].length;
					if (last+2>=length) {
						int[] newTrack = new int[Math.max(length<<1,last+3)];
						System.arraycopy(repeatTrack[currentRead],0,newTrack,0,length);
						repeatTrack[currentRead]=newTrack;
					}
					repeatTrack[currentRead][last+1]=currentStart;
					repeatTrack[currentRead][last+2]=currentEnd;
					lastRepeatTrack[currentRead]+=2;
				}
				currentRead=read; currentStart=start; currentEnd=end;
			}
			else if (start>currentEnd) {
				last=lastRepeatTrack[currentRead];
				length=repeatTrack[currentRead].length;
				if (last+2>=length) {
					int[] newTrack = new int[Math.max(length<<1,last+3)];
					System.arraycopy(repeatTrack[currentRead],0,newTrack,0,length);
					repeatTrack[currentRead]=newTrack;
				}
				repeatTrack[currentRead][last+1]=currentStart;
				repeatTrack[currentRead][last+2]=currentEnd;
				lastRepeatTrack[currentRead]+=2;
				currentStart=start; currentEnd=end;
			}
			else currentEnd=Math.max(currentEnd,end);
			str=br.readLine();
		}
		br.close();
		if (currentRead!=-1) {
			last=lastRepeatTrack[currentRead];
			length=repeatTrack[currentRead].length;
			if (last+2>=length) {
				int[] newTrack = new int[Math.max(length<<1,last+3)];
				System.arraycopy(repeatTrack[currentRead],0,newTrack,0,length);
				repeatTrack[currentRead]=newTrack;
			}
			repeatTrack[currentRead][last+1]=currentStart;
			repeatTrack[currentRead][last+2]=currentEnd;
			lastRepeatTrack[currentRead]+=2;
		}
	}
	
}