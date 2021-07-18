package de.mpi_cbg.revant.factorize;

import java.io.*;
import java.util.*;
import java.awt.image.*;
import javax.imageio.*;

/**
 *
 */
public class StudyPeriodic {
	
	/**
     * Plotting parameters
	 */
	private static final int QUANTUM_ROWS = 100;
	private static final int QUANTUM_COLUMNS = 1;
	private static final int COLOR_A = 0x00000000; //0x00FF0000;
	private static final int COLOR_T = 0x00000000; //0x00FF8000;
	private static final int COLOR_C = 0x00000000; //0x00004080;
	private static final int COLOR_G = 0xFFFFFFFF; //0x000080FF;
	
	/**
	 * IO parameters
	 */
	private static String tableFile, stringsFile;
	private static int[][] table;
	private static int nRows;
	private static final int nColumns = 6;


	public static void main(String[] args) throws Exception {
		int i, j;
		int minA, maxA, lengthB, maxLengthB;
		int color;
		String str;
		String[] tokens;
		BufferedReader br;
		BufferedImage image;
		tableFile=args[0];
		stringsFile=args[1];
		
		// Loading table
		nRows=0;
		br = new BufferedReader(new FileReader(tableFile));
		str=br.readLine();
		while (str!=null) {
			nRows++;
			str=br.readLine();
		}
		br.close();
		table = new int[nRows][nColumns];
		br = new BufferedReader(new FileReader(tableFile));
		str=br.readLine();
		i=0;
		while (str!=null) {
			tokens=str.split(",");
			for (j=0; j<nColumns; j++) table[i][j]=Integer.parseInt(tokens[j]);
			i++;
			str=br.readLine();
		}
		br.close();
		
		// Plotting the union of all readA intervals
		minA=Integer.MAX_VALUE; maxA=-1; maxLengthB=0;
		for (i=0; i<nRows; i++) {
			if (table[i][1]<minA) minA=table[i][1];
			if (table[i][2]>maxA) maxA=table[i][2];
			lengthB=table[i][5]-table[i][4]+1;
			if (lengthB>maxLengthB) maxLengthB=lengthB;
		}
		System.err.println("ReadA interval of length "+(maxA-minA+1));
		image = new BufferedImage(7*QUANTUM_ROWS,QUANTUM_COLUMNS*(maxA-minA+1),BufferedImage.TYPE_INT_RGB);
		str=getString(table[0][0],minA,maxA);
		plotString(str,image,QUANTUM_ROWS,QUANTUM_COLUMNS);
		ImageIO.write(image,"png",new File("studyPeriodic-readA.png"));
		
		// Plotting another row
		i=nRows-1;
		image = new BufferedImage(7*QUANTUM_ROWS,QUANTUM_COLUMNS*(table[i][5]-table[i][4]+1),BufferedImage.TYPE_INT_RGB);
		str=getString(table[i][3],table[i][4],table[i][5]);
		plotString(str,image,QUANTUM_ROWS,QUANTUM_COLUMNS);
		ImageIO.write(image,"png",new File("studyPeriodic-readB.png"));
				
				
		/*		
		// Plotting each readB interval
		image = new BufferedImage(QUANTUM_ROWS*nRows,QUANTUM_COLUMNS*maxLengthB,BufferedImage.TYPE_INT_RGB);
		for (i=0; i<nRows; i++) {
			str=getString(table[i][3],table[i][4],table[i][5]);
			plotString(str,image,QUANTUM_ROWS*i,QUANTUM_ROWS,QUANTUM_COLUMNS);
			System.err.println("Plotted row "+i+" out of "+nRows);
		}
		ImageIO.write(image,"png",new File("studyPeriodic-readBs.png"));
		*/
	}



	private static final String getString(int readID, int start, int end) throws IOException {
		int i, j, offset;
		String str, out;
		BufferedReader br;
		
		br = new BufferedReader(new FileReader(stringsFile));
		str=br.readLine();
		i=-1; out=null;
		while (str!=null) {
			if (str.charAt(0)!='S') {
				str=br.readLine();
				continue;
			}
			i++;
			if (i!=readID) {
				str=br.readLine();
				continue;
			}
			offset=str.indexOf(" ",2)+1;
			out=str.substring(offset+start,offset+end+1);
			break;
		}
		br.close();
		return out;
	}
	
	
	private static final void plotString(String str, BufferedImage image, int quantumRows, int quantumColumns) {
		int i, j, k;
		int length, row;
		final int color = 0xFFFFFFFF;
		
		length=str.length();
		for (i=0; i<length; i++) {
			row=-1;
			switch (str.charAt(i)) {
				case 'a': row=0; break;
				case 'c': row=2; break;
				case 'g': row=4; break;
				case 't': row=6; break;
			}
			row*=quantumRows;
			for (j=row; j<row+quantumRows; j++) {
				for (k=0; k<quantumColumns; k++) image.setRGB(j,quantumColumns*i+k,color);
			}
		}
	}
	
	/*
	private static final void plotString(String str, BufferedImage image, int row, int quantumRows, int quantumColumns) {
		int i, j, k;
		int length, color;
		
		length=str.length();
		for (i=0; i<length; i++) {
			color=-1;
			switch (str.charAt(i)) {
				case 'a': color=COLOR_A; break;
				case 'c': color=COLOR_C; break;
				case 'g': color=COLOR_G; break;
				case 't': color=COLOR_T; break;
			}
			for (j=row; j<row+quantumRows; j++) {
				for (k=0; k<quantumColumns; k++) image.setRGB(j,quantumColumns*i+k,color);
			}
		}
	}
	*/

}