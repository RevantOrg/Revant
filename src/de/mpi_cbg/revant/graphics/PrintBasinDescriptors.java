package de.mpi_cbg.revant.graphics;

import java.awt.image.*;
import java.io.*;
import javax.imageio.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;


/**
 * For every basin descriptor in a connected component of the global interval graph, the  
 * procedure draws a graphical representation of every node in the kernel of which the 
 * descriptor is a path. This is useful for debugging $IntervalGraphStep4$.
 */
public class PrintBasinDescriptors {
	
	private static final String[] SUFFIXES = new String[] {"-"+IO.LABEL_COORDINATES_LABEL+".txt","-"+IO.LABEL_COORDINATES_LABEL+"-"+IO.CYCLIC_LABEL+".txt","-"+IO.ONE_NODE_ONLY_LABEL+"-"+IO.LABEL_COORDINATES_LABEL+".txt","-"+IO.TOO_MANY_PATHS_LABEL+"-"+IO.LABEL_COORDINATES_LABEL+".txt"};
	
	/**
	 * 
	 */
	public static void main(String[] args) throws Exception {
		final String COMPONENT_ID = args[0];
		final String STEP1_DIR = args[1];
		final String READ_LENGTHS_FILE = args[2];  // All reads in the dataset
		final int N_READS = Integer.parseInt(args[3]);  // All reads in the dataset
		final int MAX_READ_LENGTH = Integer.parseInt(args[4]);
		final String INPUT_PICTURES_DIR = args[5];
		final String OUTPUT_PICTURES_DIR = args[6];
		
		final String DESCRIPTORS_DIR = STEP1_DIR+"/finalOutput";
		final String DESCRIPTOR_PREFIX = "basin-"+COMPONENT_ID+"-";
		final int X_CONTEXT = 10;
		final int Y_CONTEXT = -1;
		int i, j, p, q;
		int kernelID, nodeID, read, start, end;
		String str, fileName, clusterID;
		FileReader fr;
		BufferedReader br;
		BufferedImage snapshot;
		File directory, file;
		int[] tmpArray = new int[20];
		String[] files;
			
		Reads.nReads=N_READS;
		Reads.maxReadLength=MAX_READ_LENGTH;
		Reads.loadReadLengths(READ_LENGTHS_FILE);
		directory = new File(DESCRIPTORS_DIR);
		files=directory.list();
		for (i=0; i<files.length; i++) {
			if (files[i].indexOf(DESCRIPTOR_PREFIX)==-1) continue;
			p=DESCRIPTOR_PREFIX.length(); q=files[i].indexOf("-",p+1);
			clusterID=files[i].substring(p,q);
			br = new BufferedReader(new FileReader(DESCRIPTORS_DIR+"/"+files[i]));
			str=br.readLine();
			br.close();
			IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
			kernelID=tmpArray[0];
			fileName=STEP1_DIR+"/"+COMPONENT_ID+"-clusters/"+clusterID+"-kernel"+kernelID;
			fr=null;
			for (j=0; j<SUFFIXES.length; j++) {
				try { 
					fr = new FileReader(fileName+SUFFIXES[j]);
					break;
				}
				catch (FileNotFoundException e) { }
			}
			if (fr==null) {
				System.err.println("WARNING: kernel <"+fileName+"> not found.");
				continue;
			}
			br = new BufferedReader(fr);
			str=br.readLine(); nodeID=0;
			while (str!=null) {
				p=str.indexOf(",");
				read=Integer.parseInt(str.substring(0,p));
				p++; q=str.indexOf(",",p);
				start=Integer.parseInt(str.substring(p,q));
				p=q+1;
				end=Integer.parseInt(str.substring(p));
				snapshot=PrintReadAlignments.getIntervalContext(INPUT_PICTURES_DIR+"/read"+(read+1)+".png",read,start,end,-1,X_CONTEXT,Y_CONTEXT,tmpArray);
				if (snapshot!=null) {
					fileName=OUTPUT_PICTURES_DIR+"/"+files[i].substring(0,files[i].indexOf("."))+"-node"+nodeID+".png";
					ImageIO.write(snapshot,"png",new File(fileName));
					System.out.println("Printed file "+fileName);
				}
				str=br.readLine(); nodeID++;
			}
		}
	}

}