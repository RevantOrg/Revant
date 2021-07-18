package de.mpi_cbg.revant.util;

import java.util.Arrays;
import java.util.HashMap;
import java.io.*;
import java.util.zip.Deflater;
import de.mpi_cbg.revant.util.Point;
import de.mpi_cbg.revant.util.Points;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;


/**
 * Constants about known repeats
 */
public class Biology {
	/**
	 * Nucleosome
	 */
	public static final int NUCLEOSOME_LENGTH = 146;
	public static final int LINKER_MIN_LENGTH = 10;
	public static final int LINKER_MAX_LENGTH = 80;
	
	/**
	 * Centromere
	 */
	public static final int CENPB_BOX_LENGTH = 17;
	
	/**
	 * Telomere
	 */
	public static final int TELOMERE_PERIOD = 6;
	public static final int TELOMERE_MAX_UNITS = 4000;
	public static final int SUBTELOMERIC_REGION_MAX_LENGTH = 300000;
	
	/**
	 * Microsatellites
	 */
	public static final int MICROSATELLITE_MIN_PERIOD = 2;
	public static final int MICROSATELLITE_MAX_PERIOD = 6;
	public static final int MICROSATELLITE_MAX_LENGTH = 1000;
	
	/**
	 * Minisatellites
	 */
	public static final int MINISATELLITE_MIN_PERIOD = 10;
	public static final int MINISATELLITE_MAX_PERIOD = 60;
	public static final int MINISATELLITE_MIN_LENGTH = 1000;
	public static final int MINISATELLITE_MAX_LENGTH = 100000;
	
	/**
	 * Satellites
	 */
	public static final int SATELLITE_MIN_PERIOD = 150;
	public static final int SATELLITE_MAX_PERIOD = 180;
	public static int[] OTHER_SATELLITE_PERIODS;  // From the literature
	public static final int SATELLITE_MIN_LENGTH = 2000;
	public static final int SATELLITE_MAX_LENGTH = 5000000;  // "Several Mbp"...
	public static final int SATELLITE_HIGHERORDERBLOCK_LENGTH = 2000;
	
	/**
	 * Transposable elements
	 */
	public static final int MIR_LENGTH = 262;
	public static final int ALU_LENGTH = 300;
	public static final int MITE_MAX_LENGTH = 500;
	public static final int HAT_LENGTH = 3000;
	public static final int LINE2_LENGTH = 3300;
	public static final int LINE1_LENGTH = 6000;
	public static final int LINE1_MAX_ADJACENT_LENGTH = 1000;  // "Several kb"...
	public static final int LONG_TERMINAL_REPEATS_MIN_LENGTH = 200;
	public static final int LONG_TERMINAL_REPEATS_MAX_LENGTH = 3000;
	public static final int TERMINAL_INVERTED_REPEAT_MIN_LENGTH = 70;
	public static final int TERMINAL_INVERTED_REPEAT_MAX_LENGTH = 400;
	
	public static int maxPeriod;  // Largest period observed in the literature
	
	
	public static final void initialize(String directory) throws IOException {
		loadPeriods(directory);
		//loadRepbase(directory);
		
		maxPeriod=Math.max(TELOMERE_PERIOD,Math.max(MICROSATELLITE_MAX_PERIOD,Math.max(MINISATELLITE_MAX_PERIOD,Math.max(SATELLITE_MAX_PERIOD,OTHER_SATELLITE_PERIODS[OTHER_SATELLITE_PERIODS.length-1]))));
		System.err.println("maxPeriod="+maxPeriod);
	}
	
	
	/**
	 * @return true iff $p$ can be considered a period of a known repeat.
	 */
	public static final boolean isKnownPeriod(int p) {
		int i;
		
		if ( isPeriod(p,MICROSATELLITE_MIN_PERIOD,MICROSATELLITE_MAX_PERIOD) ||
			 isPeriod(p,MINISATELLITE_MIN_PERIOD,MINISATELLITE_MAX_PERIOD) ||
			 isPeriod(p,SATELLITE_MIN_PERIOD,SATELLITE_MAX_PERIOD) ||
			 isPeriod(p,NUCLEOSOME_LENGTH,NUCLEOSOME_LENGTH) ||
			 isPeriod(p,NUCLEOSOME_LENGTH+LINKER_MIN_LENGTH,NUCLEOSOME_LENGTH+LINKER_MAX_LENGTH) ||
		     isPeriod(p,TELOMERE_PERIOD,TELOMERE_PERIOD) ||
			 isPeriod(p,CENPB_BOX_LENGTH,CENPB_BOX_LENGTH)
		   ) return true;
		for (i=0; i<OTHER_SATELLITE_PERIODS.length; i++) {
			if (isPeriod(p,OTHER_SATELLITE_PERIODS[i],OTHER_SATELLITE_PERIODS[i])) return true;
		}
		return false;
	}
	
	
	/**
	 * @return true iff $low*i<=p<=high*i$ for some integer $i>=1$.
	 */
	private static final boolean isPeriod(int p, int low, int high) {
		int i, min, max;
		
		i=1;
		do {
			min=low*i; max=high*i;
			if (p>=min && p<=max) return true;
			i++;
		} while (min<=p);
		return false;
	}

	
	
	
	
	
	
	
	
	// --------------------------------- KNOWN PERIODS -----------------------------------
	private static int tmpPeriods;
	
	
	/**
	 * Initializes $OTHER_SATELLITE_PERIODS$ to the sorted list of all distinct periods 
	 * observed in the literature.
	 *
	 * Remark: ideally we should load just the periods that have been observed in close
	 * relatives of the genome in input. We could achieve this by adding a filter to
	 * $loadPeriods_meltersEtAl$.
	 */
	private static final void loadPeriods(String directory) throws IOException {
		int i;
		int lastPeriod1, lastPeriod2, lastPeriodMerge;
		Point[] periods1, periods2, merge;
		
		periods1=loadPeriods_cunial(directory); lastPeriod1=tmpPeriods;
		periods2=loadPeriods_meltersEtAl(directory); lastPeriod2=tmpPeriods;
		merge = new Point[lastPeriod1+1+lastPeriod2+1];
		for (i=0; i<merge.length; i++) merge[i] = new Point();
		lastPeriodMerge=Points.merge(periods1,lastPeriod1,periods2,lastPeriod2,merge);
		OTHER_SATELLITE_PERIODS = new int[lastPeriodMerge+1];
		for (i=0; i<=lastPeriodMerge; i++) OTHER_SATELLITE_PERIODS[i]=(int)merge[i].position;
		
if (IO.SHOW_STD_ERR) {
	IO.printErr("Known periods loaded: ");
	for (i=0; i<OTHER_SATELLITE_PERIODS.length; i++) IO.printErr(OTHER_SATELLITE_PERIODS[i]+", ");
	System.err.println();
}
	}
	
	
	/**
	 * Loads the file of periods collected manually by Fabio Cunial from the literature.
	 *
	 * Remark: the last period in the output is stored in $tmpPeriods$.
	 */
	private static final Point[] loadPeriods_cunial(String directory) throws IOException {
		final String PATH = directory+"/cunial.txt";
		int i, nPeriods, lastPeriod;
		String str;
		Point[] periodPoints;
		BufferedReader br;
		
		// Counting the number of periods
		br = new BufferedReader(new FileReader(PATH),IO.BUFFER_SIZE);
		nPeriods=0;
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)!='#') nPeriods++;
			str=br.readLine();
		}
		br.close();
		periodPoints = new Point[nPeriods];
		for (i=0; i<nPeriods; i++) periodPoints[i] = new Point();
		
		// Loading periods
		br = new BufferedReader(new FileReader(PATH),IO.BUFFER_SIZE);
		lastPeriod=-1;
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)=='#') {
				str=br.readLine();
				continue;
			}
			lastPeriod++;
			periodPoints[lastPeriod].position=Integer.parseInt(str.trim());
			periodPoints[lastPeriod].mass=1;
			str=br.readLine();
		}
		br.close();
		
		// Returning sorted and compacted periods
		lastPeriod=Points.sortAndCompact(periodPoints,lastPeriod);
		tmpPeriods=lastPeriod;
		return periodPoints;
	}
	
	
	/**
	 * Loads the supplementary material of:
	 *
	 * Melters, Daniel P., et al. "Comparative analysis of tandem repeats from 
	 * hundreds of species reveals unique insights into centromere evolution." Genome 
	 * biology 14.1 (2013): R10.
	 *
	 * Remark: the last period in the output is stored in $tmpPeriods$.
	 */ 
	private static final Point[] loadPeriods_meltersEtAl(String directory) throws IOException {
		final String PATH1 = directory+"/melters1.txt";
		final String PATH2 = directory+"/melters2.txt";
		int i, nPeriods, lastPeriod;
		String str;
		String[] tokens1, tokens2;
		Point[] periodPoints;
		BufferedReader br;
		
		// Computing an upper bound on the number of distinct periods in the two files
		br = new BufferedReader(new FileReader(PATH1),IO.BUFFER_SIZE);
		nPeriods=0;
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)!='#') nPeriods++;
			str=br.readLine();
		}
		br.close();
		br = new BufferedReader(new FileReader(PATH2),IO.BUFFER_SIZE);
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)=='#') {
				str=br.readLine();
				continue;
			}
			tokens1=str.split("\t");		
			if (tokens1[4].indexOf(",")<0) nPeriods++;
			else {
				tokens2=tokens1[4].split(",");
				nPeriods+=tokens2.length;
			}
			if (tokens1[5].indexOf(",")<0) nPeriods++;
			else {
				tokens2=tokens1[5].split(",");
				nPeriods+=tokens2.length;
			}
			str=br.readLine();
		}
		br.close();
		periodPoints = new Point[nPeriods];
		for (i=0; i<nPeriods; i++) periodPoints[i] = new Point();
		
		// Loading all periods
		lastPeriod=-1;
		br = new BufferedReader(new FileReader(PATH1),IO.BUFFER_SIZE);
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)=='#') {
				str=br.readLine();
				continue;
			}
			tokens1=str.split("\t");
			lastPeriod++;
			periodPoints[lastPeriod].position=Integer.parseInt(tokens1[6].trim());
			periodPoints[lastPeriod].mass=1;
			str=br.readLine();
		}
		br.close();
		br = new BufferedReader(new FileReader(PATH2),IO.BUFFER_SIZE);
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)=='#') {
				str=br.readLine();
				continue;
			}
			tokens1=str.split("\t");			
			if (tokens1[4].indexOf(",")<0) {
				lastPeriod++;
				periodPoints[lastPeriod].position=Integer.parseInt(tokens1[4].trim());
				periodPoints[lastPeriod].mass=1;
			}
			else {
				tokens2=tokens1[4].split(",");
				for (i=0; i<tokens2.length; i++) {
					lastPeriod++;
					periodPoints[lastPeriod].position=Integer.parseInt(tokens2[i].trim());
					periodPoints[lastPeriod].mass=1;
				}
			}
			if (tokens1[5].indexOf(",")<0) {
				lastPeriod++;
				periodPoints[lastPeriod].position=Integer.parseInt(tokens1[5].trim());
				periodPoints[lastPeriod].mass=1;
			}
			else {
				tokens2=tokens1[5].split(",");
				for (i=0; i<tokens2.length; i++) {
					lastPeriod++;
					periodPoints[lastPeriod].position=Integer.parseInt(tokens2[i].trim());
					periodPoints[lastPeriod].mass=1;
				}
			}
			str=br.readLine();
		}
		br.close();
		
		// Returning the set of distinct periods
		lastPeriod=Points.sortAndCompact(periodPoints,lastPeriod);
		tmpPeriods=lastPeriod;
		return periodPoints;
	}
	
	
	
	
	

}