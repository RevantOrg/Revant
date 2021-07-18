package de.mpi_cbg.revant.fragments;

import java.util.Arrays;
import java.util.HashSet;
import java.io.*;
import java.awt.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.intervalgraph.IntervalGraph;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;


public class PippelFilterAlignments {

	/**
	 * @param args 
	 * 0: root dir of Pippel's fragment-reference alignments to be filtered;
	 * 1: Fabio's directory from which fragments and references were copied; it must 
	 *    contain a $basins$ subfolder with basin descriptors.
	 */
	public static void main(String[] args) throws IOException {
		final String FRAGMENTS_ALIGNMENTS_DIR = args[0];
		final String FABIO_ROOT_DIR = args[1];
		
		final String REFERENCE_LENGTHS_DIR = FABIO_ROOT_DIR+"/references-strings-new/lengths";
		final String BASINS_DIR = FABIO_ROOT_DIR+"/basins";
		final int CAPACITY = 100;  // Arbitrary
		final int POOL_CAPACITY = 1000;  // Arbitrary
		final int MIN_CONFLICT_INTERSECTION = IO.quantum<<1;  // Arbitrary
		
		int i, p;
		int length, nReferences, basinType;
		String str, id;
		String BASINS_FILE, REFERENCE_LENGTHS_FILE, FR_ALIGNMENTS_FILE;
		File file;
		BufferedReader br;
		BufferedWriter bw;
		int[] referenceLengths = new int[CAPACITY];
		FragmentsStep3.stack = new int[POOL_CAPACITY];
		int[] tmpArray = new int[10];
		String[] list;
		FragmentsStep1.PermutationWindow[] tmpWindows = new FragmentsStep1.PermutationWindow[CAPACITY];
		FragmentsStep3.edgePool = new IntervalGraph.Edge[POOL_CAPACITY];
		FragmentsStep3.alignments = new FragmentsStep1.FFAlignment[CAPACITY];
		
		IO.initialize();
		for (i=0; i<FragmentsStep3.alignments.length; i++) FragmentsStep3.alignments[i] = new FragmentsStep1.FFAlignment();
		for (i=0; i<tmpWindows.length; i++) tmpWindows[i] = new FragmentsStep1.PermutationWindow();
		for (i=0; i<FragmentsStep3.edgePool.length; i++) FragmentsStep3.edgePool[i] = new IntervalGraph.Edge();
		file = new File(FRAGMENTS_ALIGNMENTS_DIR);
		list=file.list();
		for (i=0; i<list.length; i++) {
			p=list[i].indexOf("_");
			if (p<0) continue;
			p=list[i].indexOf("_",p+1);
			if (p<0) continue;
			if (!((new File(FRAGMENTS_ALIGNMENTS_DIR+"/"+list[i])).isDirectory())) continue;
			id=list[i].replace("_","-");
			System.err.print("Loading data structures for repeat "+id+"... ");
			
			// Basin descriptor
			BASINS_FILE=BASINS_DIR+"/"+IO.BASIN_PREFIX+id+".txt";
			br = new BufferedReader(new FileReader(BASINS_FILE));
			str=br.readLine(); br.close();
			IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
			basinType=tmpArray[3];
			
			// Loading reference lengths
			REFERENCE_LENGTHS_FILE=REFERENCE_LENGTHS_DIR+"/reference-"+id+"-lengths.txt";
			br = new BufferedReader(new FileReader(REFERENCE_LENGTHS_FILE));
			nReferences=0; str=br.readLine();
			while (str!=null) {
				length=Integer.parseInt(str);
				nReferences++;
				if (nReferences==referenceLengths.length) {
					int[] newLengths = new int[referenceLengths.length<<1];
					System.arraycopy(referenceLengths,0,newLengths,0,referenceLengths.length);
					referenceLengths=newLengths;
				}
				referenceLengths[nReferences]=length;
				str=br.readLine();
			}
			br.close();
			
			// Building conflict graphs
			FR_ALIGNMENTS_FILE=FRAGMENTS_ALIGNMENTS_DIR+"/"+list[i]+"/"+IO.FRAGMENTS_REFERENCE_FILE+".txt";
			FragmentsStep3.filterAlignments(basinType,FR_ALIGNMENTS_FILE,MIN_CONFLICT_INTERSECTION,referenceLengths,true,null,false);
		}
		
		// Symmetrizing bitvectors
		buildSymmetricBitvectors(FRAGMENTS_ALIGNMENTS_DIR);
	}
	
	
	/**
	 * Uses the bitvectors built for the fragment-reference alignments file, to build the 
	 * bitvectors of every corresponding reference-fragment alignment file.
	 *
	 * @param path root directory of fragments-reference alignments.
	 */
	private static final void buildSymmetricBitvectors(String path) throws IOException {
		int i, p;
		int nOnes;
		String str1, str2;
		BufferedReader br1, br2;
		BufferedWriter bw;
		File file;
		FragmentsStep1.FFAlignment tmpAlignment, newAlignment;
		HashSet<FragmentsStep1.FFAlignment> set;
		String[] list;
		
		tmpAlignment = new FragmentsStep1.FFAlignment();
		set = new HashSet<FragmentsStep1.FFAlignment>(10000);
		file = new File(path);
		list=file.list();
		for (i=0; i<list.length; i++) {
			p=list[i].indexOf("_");
			if (p<0) continue;
			p=list[i].indexOf("_",p+1);
			if (p<0) continue;
			if (!((new File(path+"/"+list[i])).isDirectory())) continue;
			if (!file.exists()) {
				System.err.println("WARNING: file not found: "+path+"/"+list[i]+"/"+IO.REFERENCE_FRAGMENTS_FILE+".txt");
				continue;
			}
			// Loading fragment-reference alignments
			set.clear();
			br1 = new BufferedReader(new FileReader(path+"/"+list[i]+"/"+IO.FRAGMENTS_REFERENCE_FILE+".txt"));
			br2 = new BufferedReader(new FileReader(path+"/"+list[i]+"/"+IO.FRAGMENTS_REFERENCE_FILE+"-keep.txt"));
			str1=br1.readLine(); str1=br1.readLine();  // Skipping header
			str1=br1.readLine();
			str2=br2.readLine();
			while (str2!=null) {
				if (str2.charAt(0)=='1') {
					Alignments.readAlignmentFile(str1);
					newAlignment = new FragmentsStep1.FFAlignment(Alignments.readA,Alignments.readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
					newAlignment.inCanonicalForm();
					set.add(newAlignment);
				}
				str1=br1.readLine(); str2=br2.readLine();
			}
			br1.close(); br2.close();
			
			// Filtering reference-fragment alignments
			br1 = new BufferedReader(new FileReader(path+"/"+list[i]+"/"+IO.REFERENCE_FRAGMENTS_FILE+".txt"));
			bw = new BufferedWriter(new FileWriter(path+"/"+list[i]+"/"+IO.REFERENCE_FRAGMENTS_FILE+"-keep.txt"));
			str1=br1.readLine(); str1=br1.readLine();  // Skipping header
			str1=br1.readLine(); nOnes=0;
			while (str1!=null) {
				Alignments.readAlignmentFile(str1);
				tmpAlignment.set(Alignments.readA,Alignments.readB,Alignments.orientation,Alignments.startA,Alignments.endA,Alignments.startB,Alignments.endB);
				tmpAlignment.inCanonicalForm();
				if (set.contains(tmpAlignment)) {
					nOnes++;
					bw.write("1\n");
				}
				else bw.write("0\n");
				str1=br1.readLine();
			}
			br1.close(); bw.close();
			if (nOnes!=set.size()) {
				System.err.println("ERROR: written "+nOnes+" ones instead of "+set.size()+", module "+list[i]);
				System.exit(1);
			}
		}
		set.clear();
	}

}