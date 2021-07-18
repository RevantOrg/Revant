package de.mpi_cbg.revant.factorize;

import java.util.*;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;

/**
 * 
 */
public class PAF2Repeats {

	
	public static void main(String[] args) throws IOException {
		final String REFERENCE_REPEATS_FILE = args[0];
		final String ALIGNMENTS_FILE = args[1];
		
		final double SURFACE_THRESHOLD = 0.9;  // Arbitrary
		
		int p;
		int start, previousStart, previousEnd, previousLength, previousSurface, withAlignment, withLowSurface;
		String str, key, reference, previousReference;
		HashMap<String,Integer> sequenceNames;
		BufferedReader br;
		BufferedWriter bw;
		String[] tokens;
		
		IO.initialize();
		System.out.print("Loading names of reference repeats... ");
		sequenceNames = new HashMap<String,Integer>(1000);
		br = new BufferedReader(new FileReader(REFERENCE_REPEATS_FILE));
		str=br.readLine();
		while (str!=null) {
			if (str.length()>0 && str.charAt(0)=='>') sequenceNames.put(str.substring(1),Integer.valueOf(0));
			str=br.readLine();
		}
		br.close();
		System.out.println(" DONE");
		
		System.out.println("Checking reference->mine alignments... ");
		br = new BufferedReader(new FileReader(ALIGNMENTS_FILE));
		previousReference=null; previousStart=-1; previousEnd=-1; 
		previousLength=0; previousSurface=0; withLowSurface=0;
		str=br.readLine();
		while (str!=null) {
			tokens=str.split("\t");
			reference=tokens[0];
			if (sequenceNames.containsKey(reference)) sequenceNames.put(reference,Integer.valueOf(1));
			if (previousReference==null) {
				previousReference=reference;
				previousLength=Integer.parseInt(tokens[1]);
				previousStart=Integer.parseInt(tokens[2]);
				previousEnd=Integer.parseInt(tokens[3]);
				previousSurface=0;
			}
			else if (!reference.equals(previousReference)) {
				previousSurface+=previousEnd-previousStart+1;
				if (previousSurface<previousLength*SURFACE_THRESHOLD) {
					withLowSurface++;
					System.out.println("WARNING: reference <"+previousReference+"> is covered just by "+IO.getPercent(previousSurface,previousLength)+"%");				
				}
				previousReference=reference;
				previousLength=Integer.parseInt(tokens[1]);
				previousStart=Integer.parseInt(tokens[2]);
				previousEnd=Integer.parseInt(tokens[3]);
				previousSurface=0;
			}
			else {
				start=Integer.parseInt(tokens[2]);
				if (start<=previousEnd) previousEnd=Math.max(previousEnd,Integer.parseInt(tokens[3]));
				else {
					previousSurface+=previousEnd-previousStart+1;
					previousStart=start;
					previousEnd=Integer.parseInt(tokens[3]);
				}
			}
			str=br.readLine();
		}
		previousSurface+=previousEnd-previousStart+1;
		if (previousSurface<previousLength*SURFACE_THRESHOLD) {
			withLowSurface++;
			System.out.println("WARNING: reference <"+previousReference+"> is covered just by "+IO.getPercent(previousSurface,previousLength)+"%");
		}
		br.close();
		
		// Counting the number of references that have an alignment
		Set<Map.Entry<String,Integer>> set = sequenceNames.entrySet();
		Iterator<Map.Entry<String,Integer>> iterator = set.iterator();
		withAlignment=0;
		while (iterator.hasNext()) {
			Map.Entry<String,Integer> entry = iterator.next();
			if (entry.getValue().intValue()==1) withAlignment++;
			else System.out.println(entry.getKey()+" has no alignment to our repeats");
		}
		System.out.println(withAlignment+" out of "+sequenceNames.size()+" reference repeats have at least one alignment with my repeats ("+IO.getPercent(withAlignment,sequenceNames.size())+"%)");
		System.out.println(withLowSurface+" out of "+sequenceNames.size()+" reference repeats have low covered surface ("+IO.getPercent(withLowSurface,sequenceNames.size())+"%)");
	}
	

}