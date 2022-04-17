package de.mpi_cbg.revant.biology;

import java.util.Arrays;
import java.io.*;
import de.mpi_cbg.revant.util.IO;


/**
 * Removes from a Repbase file all characters that are not in the DNA alphabet. Pads with
 * multiple copies every satellite that is shorter than the min length of an alignment, 
 * as well as every non-satellite whose padding would make an alignment contain at least a
 * given number of copies of the repeat (such repeats become satellites). Non-satellites
 * for which such padding is not possible are discarded. Makes all headers compatible with
 * PacBio conventions. Outputs an array of repeat lengths, and a bitvector that marks
 * repeats that are periodic after the filtering.
 */
public class FilterRepbase {
	
	private static final String SATELLITE_LABEL = "sat";
	
	/**
	 * @param args
	 *  
	 * 2: if a non-satellite repeat is shorter than minAlignmentLength, and if at least 
	 * this number of copies >1 fit into an alignment, the repeat is padded with 
	 * copies of itself and considered periodic, otherwise it is discarded; use any number
	 * <=1 to discard every non-satellite repeat that is shorter than minAlignmentLength;
	 *
	 * 3: discard every repeat that does not contain this keyword (use NULL to disable the
	 * filter).
	 */
 	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
		final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[1]);
		final double MIN_COPIES_IN_ALIGNMENT = Double.parseDouble(args[2]);
		final String KEYWORD = args[3].equalsIgnoreCase("null")?null:args[3].toLowerCase();
		final String OUTPUT_FILE = args[4];
		final String OUTPUT_FILE_ISPERIODIC = args[5];
		final String OUTPUT_FILE_LENGTHS = args[6];
		
		boolean isPeriodic;
		int idGenerator, length;
		String str, header;
		StringBuilder buffer1, buffer2;
		BufferedReader br;
		BufferedWriter bw1, bw2, bw3;
		
		buffer1 = new StringBuilder();
		buffer2 = new StringBuilder();
		br = new BufferedReader(new FileReader(INPUT_FILE));
		bw1 = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		bw2 = new BufferedWriter(new FileWriter(OUTPUT_FILE_ISPERIODIC));
		bw3 = new BufferedWriter(new FileWriter(OUTPUT_FILE_LENGTHS));
		idGenerator=-1; header="";
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)=='>') {
				if (KEYWORD==null || header.indexOf(KEYWORD)>=0) {
					IO.removeNonDNACharacters(buffer1);
					length=buffer1.length();
					if (length!=0) {
						isPeriodic=header.indexOf(SATELLITE_LABEL)>=0;
						if (length<MIN_ALIGNMENT_LENGTH) {
							if (isPeriodic || (MIN_COPIES_IN_ALIGNMENT>1.0 && MIN_ALIGNMENT_LENGTH>=length*MIN_COPIES_IN_ALIGNMENT)) {
								isPeriodic=true;
								buffer2.delete(0,buffer2.length());
								buffer2.append(buffer1);
								while (buffer2.length()<MIN_ALIGNMENT_LENGTH) buffer2.append(buffer1);
								IO.writeFakeHeader(++idGenerator,buffer2.length(),header,bw1);
								bw1.write(buffer2+"\n");
								bw3.write(buffer2.length()+"\n");
							}
							else {
								buffer1.delete(0,length);
								header=str.substring(1).trim().toLowerCase();
								str=br.readLine();
								continue;
							}
						}
						else {
							IO.writeFakeHeader(++idGenerator,length,header,bw1);
							bw1.write(buffer1+"\n");
							bw3.write(length+"\n");
						}
						bw2.write((isPeriodic?"1":"0")+"\n");
						buffer1.delete(0,length);
					}
				}
				header=str.substring(1).trim().toLowerCase();
			}
			else buffer1.append(str.toLowerCase());
			str=br.readLine();
		}
		br.close();
		if (KEYWORD==null || header.indexOf(KEYWORD)>=0) {
			IO.removeNonDNACharacters(buffer1);
			length=buffer1.length();
			isPeriodic=header.indexOf(SATELLITE_LABEL)>=0;
			if (length<MIN_ALIGNMENT_LENGTH) {
				if (isPeriodic || (MIN_COPIES_IN_ALIGNMENT>1.0 && MIN_ALIGNMENT_LENGTH>=length*MIN_COPIES_IN_ALIGNMENT)) {
					isPeriodic=true;
					buffer2.delete(0,buffer2.length());
					buffer2.append(buffer1);
					while (buffer2.length()<MIN_ALIGNMENT_LENGTH) buffer2.append(buffer1);
					IO.writeFakeHeader(++idGenerator,buffer2.length(),header,bw1);
					bw1.write(buffer2+"\n");
					bw3.write(buffer2.length()+"\n");
				}
				else { /*NOP*/ }
			}
			else {
				IO.writeFakeHeader(++idGenerator,length,header,bw1);
				bw1.write(buffer1+"\n");
				bw3.write(length+"\n");
			}
			bw2.write((isPeriodic?"1":"0")+"\n");
		}
		bw1.close(); bw2.close(); bw3.close();
 	}
	
}