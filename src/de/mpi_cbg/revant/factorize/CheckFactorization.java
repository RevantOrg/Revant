package de.mpi_cbg.revant.factorize;

import java.io.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Constants;


/**
 * Checks whether the intervals assigned to alignments in the interval-connection file
 * also occur in the interval-* files.
 */
public class CheckFactorization {


	private static int[][] alignments, dense, periodic;
	

	/**
	 * java CheckFactorization ~/Desktop/planarian/LAshow-3.txt ~/Desktop/planarian/connection-3.txt ~/Desktop/planarian/intervals-alignments-3.txt ~/Desktop/planarian/intervals-dense-3.txt ~/Desktop/planarian/intervals-periodic-3.txt
	 */
	public static void main(String[] args) throws IOException {
		final String LASHOW_FILE_PATH = args[0];
		final String CONNECTION_FILE_PATH = args[1];
		final String ALIGNMENTS_FILE_PATH = args[2];
		final String DENSE_FILE_PATH = args[3];
		final String PERIODIC_FILE_PATH = args[4];
		
		int i, p, inconsistent;
		int nAlignments, nDense, nPeriodic;
		String str, strPrime;
		BufferedReader alignmentsFile, denseFile, periodicFile, lashowFile, connectionFile;
		int[] tmp = new int[100];
		boolean[] alignmentsUsed, denseUsed, periodicUsed;
		
		// Loading alignments
		alignmentsFile = new BufferedReader(new FileReader(ALIGNMENTS_FILE_PATH),IO.BUFFER_SIZE);
		nAlignments=0;
		str=alignmentsFile.readLine();
		while (str!=null) {
			nAlignments++;
			str=alignmentsFile.readLine();
		}
		alignmentsFile.close();
		alignments = new int[nAlignments][4];  // read,id,start,end
		alignmentsUsed = new boolean[nAlignments];
		for (i=0; i<nAlignments; i++) alignmentsUsed[i]=false;
		alignmentsFile = new BufferedReader(new FileReader(ALIGNMENTS_FILE_PATH),IO.BUFFER_SIZE);
		i=0;
		str=alignmentsFile.readLine();
		while (str!=null) {
			p=str.indexOf(",");
			alignments[i][0]=Integer.parseInt(str.substring(0,p));
			AlignmentInterval.readAlignmentsFile(str,p+1,tmp,0);
			alignments[i][1]=tmp[0];
			alignments[i][2]=tmp[1];
			alignments[i][3]=tmp[2];
			i++;
			str=alignmentsFile.readLine();
		}
		alignmentsFile.close();
		
		// Loading dense
		denseFile = new BufferedReader(new FileReader(DENSE_FILE_PATH),IO.BUFFER_SIZE);
		nDense=0;
		str=denseFile.readLine();
		while (str!=null) {
			nDense++;
			str=denseFile.readLine();
		}
		denseFile.close();
		dense = new int[nDense][5];  // read,id,start,end,type
		denseUsed = new boolean[nDense];
		for (i=0; i<nDense; i++) denseUsed[i]=false;
		denseFile = new BufferedReader(new FileReader(DENSE_FILE_PATH),IO.BUFFER_SIZE);
		i=0;
		str=denseFile.readLine();
		while (str!=null) {
			p=str.indexOf(",");
			dense[i][0]=Integer.parseInt(str.substring(0,p));
			DenseSubstring.readDenseSubstringsFile(str,p+1,tmp,0);
			dense[i][1]=tmp[0];
			dense[i][2]=tmp[2];
			dense[i][3]=tmp[3];
			dense[i][4]=DenseSubstring.getType(tmp);
			i++;
			str=denseFile.readLine();
		}
		denseFile.close();
		
		// Loading periodic
		periodicFile = new BufferedReader(new FileReader(PERIODIC_FILE_PATH),IO.BUFFER_SIZE);
		nPeriodic=0;
		str=periodicFile.readLine();
		while (str!=null) {
			nPeriodic++;
			str=periodicFile.readLine();
		}
		periodicFile.close();
		periodic = new int[nPeriodic][4];  // read,id,start,end
		periodicUsed = new boolean[nPeriodic];
		for (i=0; i<nPeriodic; i++) periodicUsed[i]=false;
		periodicFile = new BufferedReader(new FileReader(PERIODIC_FILE_PATH),IO.BUFFER_SIZE);
		i=0;
		str=periodicFile.readLine();
		while (str!=null) {
			p=str.indexOf(",");
			periodic[i][0]=Integer.parseInt(str.substring(0,p));
			PeriodicSubstringInterval.readPeriodicSubstringsFile(str,p+1,tmp,0);
			periodic[i][1]=tmp[0];
			periodic[i][2]=tmp[1];
			periodic[i][3]=tmp[2];
			i++;
			str=periodicFile.readLine();
		}
		periodicFile.close();
		
		// Checking connection file
		lashowFile = new BufferedReader(new FileReader(LASHOW_FILE_PATH),IO.BUFFER_SIZE);
		connectionFile = new BufferedReader(new FileReader(CONNECTION_FILE_PATH),IO.BUFFER_SIZE);
		str=lashowFile.readLine(); str=lashowFile.readLine();  // Skipping header
		str=lashowFile.readLine();
		strPrime=connectionFile.readLine();
		i=0;
		while (str!=null && strPrime!=null) {
			if (strPrime.length()!=0) {
				Alignments.readAlignmentFile(str);
				Factorize.readIntervalConnectionFile(strPrime);
				if ( (Factorize.type==Constants.INTERVAL_ALIGNMENT && !find(alignments,Alignments.readA-1,Factorize.id,Factorize.start,Factorize.end,-1,alignmentsUsed)) ||
					 (Factorize.type>=Constants.INTERVAL_DENSE_PREFIX && Factorize.type<=Constants.INTERVAL_DENSE_SINGLEDELETION && !find(dense,Alignments.readA-1,Factorize.id,Factorize.start,Factorize.end,Factorize.type,denseUsed)) ||
					 (Factorize.type==Constants.INTERVAL_PERIODIC && !find(periodic,Alignments.readA-1,Factorize.id,Factorize.start,Factorize.end,-1,periodicUsed))
				   ) {
					System.err.println("Line "+(i+1)+" of the connection file not found (starting from one): read="+(Alignments.readA-1)+" type="+Factorize.type+" id="+Factorize.id+" start="+Factorize.start+" end="+Factorize.end);
					System.exit(1);
				}
			}
			str=lashowFile.readLine();
			strPrime=connectionFile.readLine();
			i++;
		}
		System.out.println("Connection file is consistent with interval files");
		System.out.println();
		
		// Checking $*used$ bitvectors
		inconsistent=0;
		for (i=0; i<nAlignments; i++) {
			if (!alignmentsUsed[i]) {
				inconsistent++;
				System.err.println("Alignment interval "+i+" never used in the connection file: read="+alignments[i][0]+" id="+alignments[i][1]+" start="+alignments[i][2]+" end="+alignments[i][3]);
			}
		}
		for (i=0; i<nDense; i++) {
			if (!denseUsed[i]) {
				inconsistent++;
				System.err.println("Dense interval "+i+" never used in the connection file: read="+dense[i][0]+" id="+dense[i][1]+" start="+dense[i][2]+" end="+dense[i][3]);
			}
		}
		for (i=0; i<nPeriodic; i++) {
			if (!periodicUsed[i]) {
				inconsistent++;
				System.err.println("Periodic interval "+i+" never used in the connection file: read="+periodic[i][0]+" id="+periodic[i][1]+" start="+periodic[i][2]+" end="+periodic[i][3]);
			}
		}
		System.out.println(inconsistent+" inconsistencies found between interval files and connection file");
	}
	
	
	private static final boolean find(int[][] matrix, int read, int id, int start, int end, int type, boolean[] usedMarks) {
		for (int i=0; i<matrix.length; i++) {
			if (matrix[i][0]<read) continue;
			if (matrix[i][0]>read) break;
			if (matrix[i][1]==id && matrix[i][2]==start && matrix[i][3]==end && (matrix==dense?matrix[i][4]==type:true)) {
				usedMarks[i]=true;
				return true;
			}
		}
		return false;
	}

}