package de.mpi_cbg.revant.util;

import java.util.Arrays;
import java.util.Random;
import java.text.NumberFormat;
import java.util.Locale;
import java.io.BufferedWriter;
import java.io.BufferedOutputStream;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;


public class IO {
	
	public static final boolean SHOW_STD_OUT = false;
	public static final boolean SHOW_STD_ERR = false;
	public static final boolean SHOW_STD_ERR_PRIME = false;
	public static final boolean SHOW_INTERACTIVE = false;
	public static final boolean SHOW_INTERACTIVE_2 = false;
	public static final boolean CONSISTENCY_CHECKS = true;
	public static final int MAX_OPEN_FILES = 1024;
	public static final int BUFFER_SIZE = 1024;
	
	public static final String STEP1_DIR = "step1";
	public static final String TAGS_DIR = "finalOutput";
	public static final String STEP4_DIR = "step4";
	public static final String STEP5_DIR = "step5";
	
	public static final String ALIGNMENTS_FILE_STEP1 = "LAshow.txt";	
	public static final String ALIGNMENTS_FILE = "alignments.txt";
	public static final String ALIGNMENTS_SHORTPERIOD = "alignments-shortPeriod.txt";
	
	public static final String CONNECTION_FILE = "connection-all.txt";
	public static final String PERIODIC_INTERVALS = "intervals-periodic-all.txt";
	public static final String DENSE_INTERVALS = "intervals-dense-all.txt";
	public static final String ALIGNMENT_INTERVALS = "intervals-alignments-all.txt";
	
	public static final String READS_IDS = "reads-ids.txt";
	public static final String READS_LENGTHS = "reads-lengths.txt";
	public static final String READS_PHRED_PREFIX = "reads-phred";
	public static final String READS_PHRED_SUFFIX_DAMAR = ".damar";
	public static final String READS_PHRED_SUFFIX_DBDUMP = ".dbdump";
	public static final String READS_SHORTPERIOD = "reads-shortPeriod.txt";
	
	public static final String KERNEL2DESCRIPTOR_FILE = "newKernel2oldKernel.txt";
	public static final String BASIN_DESCRIPTOR_PREFIX = "basin";
	public static final String TAGS_PREFIX = "tags-";
	public static final String TAGS_ROOT_FILE = TAGS_PREFIX+"root.txt";
	
	public static final String LABELS_LABEL = "labels";
	public static final String LABEL_COORDINATES_LABEL = "labelCoordinates";
	public static final String ONE_NODE_ONLY_LABEL = "oneNodeOnly";
	public static final String CYCLIC_LABEL = "cyclic";
	public static final String TOO_MANY_PATHS_LABEL = "tooManyPaths";
	
	public static final String REFERENCE_LABEL = "reference";
	public static final String FRAGMENTS_LABEL = "fragments";
	public static final String MODULE_ALIGNMENTS_DIR_PREFIX = "test-basin-fragments-";
	public static final int MODULE_ALIGNMENTS_DIR_PREFIX_LENGTH = MODULE_ALIGNMENTS_DIR_PREFIX.length();
	public static final String BASIN_PREFIX = "basin-";
	public static final int BASIN_PREFIX_LENGTH = BASIN_PREFIX.length();
	
	public static final String STEP5_DESCRIPTORS_FILE = "descriptors.txt";
	public static final String STEP5_ADJACENCY_FILE = "adjacencies.txt";
	
	public static final String DBDUMP_QUALITY_EXTENSION = ".dbdump";
	public static final String DAMAR_QUALITY_EXTENSION = ".damar";
	
	public static final String FRAGMENTS_REFERENCE_FILE = "LAshow-fragments-reference";
	public static final String REFERENCE_FRAGMENTS_FILE = "LAshow-reference-fragments";
	
	public static final char[] DNA_ALPHABET_LOWERCASE = new char[] {'a','c','g','t'};


	/**
	 * Global parameters of the pipeline
	 */
	public static int minOccurrencesInGenome;  // Minimum number of occurrences in the genome for a substring to be considered a repeat.
	public static int coverage;  // Coverage of the genome
	public static int minCoverage;  // Minimum coverage for a substring not to be considered a random insertion
	public static int minRepeatCoverage;  // Minimum coverage for a substring to be considered a repeat
	public static int quantum = 100;
	public static double insertionRate = 0.1;
	public static double deletionRate = 0.05;
	public static int maxDeletionLength = 15;  // Maximum length of a deletion in a read
	public static int maxInsertionLength = 20;  // Maximum length of a short insertion in a read
	public static double maxDeltaStd = 0.02;  // Regression tree intervals on function $\delta$ are not split if their standard deviation is at most this much
	public static boolean includeSelfAlignments;


	private static NumberFormat formatter;
	private static StringBuilder sb;
	
	
	
	public static final void initialize() {
		formatter = NumberFormat.getNumberInstance(Locale.ENGLISH);
		//formatter.setMinimumFractionDigits(4);
		formatter.setMaximumFractionDigits(4);
		formatter.setGroupingUsed(false);
		sb = new StringBuilder();
	}

	
	public static final String getPhredSuffix(String prefix) throws IOException {
		File file = new File(prefix+READS_PHRED_SUFFIX_DAMAR);
		if (file.exists()) return READS_PHRED_SUFFIX_DAMAR;
		else return READS_PHRED_SUFFIX_DBDUMP;
	}

	
	public static final String format(double number) {
		return formatter.format(number);
	}
	
	
	public static final void printOut(Object obj) {
		System.out.println(obj.toString());
	}
	
	
	public static final void printErr(Object obj) {
		System.err.println(obj.toString());
	}
	
	
	public static final void printCriticalErr(Object obj) {
		System.err.println(obj.toString());
	}
	
	
	public static final String printRow(int[] ids) {
		final int length = ids.length;
		if (length==0) return "";
		sb.delete(0,sb.length());
		for (int i=0; i<ids.length-1; i++) {
			sb.append(ids[i]);
			sb.append('-');
		}
		sb.append(ids[length-1]);
		return sb.toString();
	}
	
	
	public static final void writeInt(int n, BufferedOutputStream out) throws IOException {
		int mask = 0x000000FF;
		int shift = 0;
		for (int i=0; i<4; i++) {
			out.write((n&mask)>>>shift);
			mask<<=8;
			shift+=8;
		}
	}
	
	
	public static final int readInt(BufferedInputStream in) throws IOException {
		final int mask = 0x000000FF;
		int out = 0x00000000;
		int shift = 0;
		for (int i=0; i<4; i++) {
			out|=(in.read()&mask)<<shift;
			shift+=8;
		}
		return out;
	}
	
	
	public static final void writeLong(long n, BufferedOutputStream out) throws IOException {
		long mask = 0x00000000000000FF;
		int shift = 0;
		for (int i=0; i<8; i++) {
			out.write((int)((n&mask)>>>shift));
			mask<<=8;
			shift+=8;
		}
	}
	
	
	public static final long readLong(BufferedInputStream in) throws IOException {
		final long mask = 0x00000000000000FF;
		long out = 0x0000000000000000;
		int shift = 0;
		long c;
		for (int i=0; i<8; i++) {
			c=in.read();
			out|=(c&mask)<<shift;
			shift+=8;
		}
		return out;
	}
	
	
	public static final void writeDouble(double n, BufferedOutputStream out) throws IOException {
		writeLong(Double.doubleToLongBits(n),out);
	}
	
	
	public static final double readDouble(BufferedInputStream in) throws IOException {
		return Double.longBitsToDouble(readLong(in));
	}


	public static final String getPercent(long numerator, long denominator) {
		return denominator==0?"0":format(100*(((double)numerator)/denominator));
	}
	
	
	/**
	 * Remark: alignment intervals are in the forward orientation of both reads.
	 *
	 * @return the canonical concatenation of the coordinates that define an alignment.
	 */
	public static final String getCanonicalForm(int readA, int readB, boolean orientation, int startA, int endA, int startB, int endB) {
		if (readA<readB) return readA+","+startA+","+endA+","+readB+","+startB+","+endB+","+(orientation?"1":"0");
		else return readB+","+startB+","+endB+","+readA+","+startA+","+endA+","+(orientation?"1":"0");
	}
	
	
	/**
	 * Prints a header in PacBio format that works with DALIGNER programs.
	 *
	 * @param suffix arbitrary suffix, can be null.
	 */
	public static final void writeFakeHeader(int id, int stringLength, String suffix, BufferedWriter bw) throws IOException {
		bw.write(">U0/"+id+"/0_"+stringLength+(suffix!=null?" "+suffix:"")+"\n");
	}
	
	
	/**
	 * @param sb assumed to contain only lowercase characters.
	 */
	public static final void removeNonDNACharacters(StringBuilder sb) {
		char c;
		int i;
		final int alphabetLength = DNA_ALPHABET_LOWERCASE.length;
		final int length = sb.length();
		
		sb.ensureCapacity(length<<1);
		for (i=0; i<length; i++) {
			c=sb.charAt(i);
			if (Arrays.binarySearch(DNA_ALPHABET_LOWERCASE,0,alphabetLength,c)>=0) sb.append(c);
		}
		sb.delete(0,length);
	}

}