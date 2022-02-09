package de.mpi_cbg.revant.apps;

import java.io.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;


/**
 * Splits the set of all character instances into approx. equal-sized parts, aligned with
 * repeat boundaries.
 */
public class SplitCharacterInstances {

	public static void main(String[] args) throws IOException {
		final int N_INSTANCES = Integer.parseInt(args[0]);
		final int N_PIECES = Integer.parseInt(args[1]);
		final String INPUT_FILE = args[2];
		final String OUTPUT_PREFIX = args[3];
		
		int fileID, nInstancesInFile, lastUnique, lastPeriodic;
		final int quantum = N_INSTANCES/N_PIECES;
		String str;
		BufferedReader br;
		BufferedWriter outputFile, headerFile;
		RepeatAlphabet.Character previousCharacter, character;
		
		fileID=0; nInstancesInFile=0;
		outputFile = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+fileID+".txt"));
		headerFile = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+fileID+"-header.txt"));
		lastUnique=-1; lastPeriodic=-1;
		previousCharacter = new RepeatAlphabet.Character();
		character = new RepeatAlphabet.Character();
		br = new BufferedReader(new FileReader(INPUT_FILE),IO.BUFFER_SIZE);
		str=br.readLine();
		while (str!=null) {
			character.deserialize(str);
			if (character.repeat!=previousCharacter.repeat && nInstancesInFile>=quantum) {
				outputFile.close();
				headerFile.write(lastUnique+(RepeatAlphabet.SEPARATOR_MINOR+"")+lastPeriodic+(RepeatAlphabet.SEPARATOR_MINOR+"")+(nInstancesInFile-1)+(RepeatAlphabet.SEPARATOR_MINOR+"")+"-1");
				headerFile.newLine(); headerFile.close();
				System.out.println(nInstancesInFile+" instances");
				fileID++;
				outputFile = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+fileID+".txt"));
				headerFile = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+fileID+"-header.txt"));
				nInstancesInFile=0; lastUnique=-1; lastPeriodic=-1;
			}
			outputFile.write(str); outputFile.newLine(); nInstancesInFile++;
			if (character.repeat==RepeatAlphabet.UNIQUE) lastUnique=nInstancesInFile-1;
			else if (character.start==-1) lastPeriodic=nInstancesInFile-1;
			previousCharacter.copyFrom(character);
			str=br.readLine();
		}
		br.close(); outputFile.close();
		headerFile.write(lastUnique+(RepeatAlphabet.SEPARATOR_MINOR+"")+lastPeriodic+(RepeatAlphabet.SEPARATOR_MINOR+"")+(nInstancesInFile-1)+(RepeatAlphabet.SEPARATOR_MINOR+"")+"-1");
		headerFile.newLine(); headerFile.close();
		fileID++;
		while (fileID<N_PIECES) {
			outputFile = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+fileID+".txt"));
			outputFile.close();
			headerFile = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+fileID+"-header.txt"));
			headerFile.write("-1"+RepeatAlphabet.SEPARATOR_MINOR+"-1"+RepeatAlphabet.SEPARATOR_MINOR+"-1"+RepeatAlphabet.SEPARATOR_MINOR+"-1");
			headerFile.close();
			fileID++;
		}
		System.out.println(nInstancesInFile+" instances");
	}

}