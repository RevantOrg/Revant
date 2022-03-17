package de.mpi_cbg.revant.biology;

import java.util.Arrays;
import java.util.Vector;
import java.util.Random;
import java.io.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.imageio.*;
import java.text.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Colors;


/**
 * Simulates a repetitive genome through several epochs. In every epoch a new repeat is 
 * generated, either by creating a random repeat from scratch, or by perturbing an 
 * instance of an existing repeat. In the latter case, the two repeats are connected by an
 * edge in an evolutionary tree. Every repeat is born with a type (identical to the types 
 * used in factorization, e.g. prefix, substring, short-period, long-period), that decides
 * how to create its instances when it replicates; with an insertion preference, that 
 * tells how likely the repeat is to insert into every one of the existing evolutionary 
 * trees; and with a desired frequency in the genome. The new repeat is replicated at 
 * random in the current genome, until its desired frequency or the total repeat bps of 
 * the genome are reached. The length of repeats and instances are decided according to 
 * the min. alignment length specified by the user, so that repeat instances can be 
 * detected using this setting of the aligner. The simulator can use a Repbase file given
 * in input, and it classifies Repbase repeats into satellites and non-satellites.
 */
public class GenomeSimulator {
	
	public static final char[] DNA_ALPHABET = new char[] {'a','c','g','t'};
	
	/**
	 * Repeats. All unique substrings of the genome are assumed to belong to the same 
	 * artificial repeat, which forms its own evolutionary tree with just one node.
	 */
	public static final int N_REPEAT_TYPES = Constants.INTERVAL_UNKNOWN+1;
	private static int DELTA;  // For single-deletion repeat type
	private static Repeat[] repeats;
	private static int lastRepeat, lastTree;
	private static Repeat unique;
	
	/**
	 * Root of the linked list that represents the genome. $root.next$ points to the
	 * leftmost real instance.
	 */
	private static RepeatInstance root;
	
	/**
	 * Row $i$ stores all the repeat instances that belong to the $i$-th tree and that 
	 * have length $>=2*minAlignmentLength$, in no particular order.
	 */
	private static RepeatInstance[][] byTree;
	private static int[] lastByTree;
	
	/**
	 * Temporary array: cell $i$ contains the number of columns in all rows of $byTree$ 
	 * that are $<=i$ (including row zero).
	 */
	private static int[] cumulativeByTree;
	
	/**
	 * Temporary space used by $getMultiInstanceString()$.
	 */
	private static StringBuilder sbLeft, sbRight, forReverse;
	
	/**
	 * Drawing constants
	 */
	private static final int ROWS_PER_LINE = 100;
	private static final int ROWS_PER_LINE_OVER_TWO = ROWS_PER_LINE>>1;
	private static final int TRUNCATION_THRESHOLD = 100;
	private static final int TRUNCATION_TAG_WIDTH_PIXELS = 3;
	private static final int TRUNCATION_MULTIPLE = 1000;
	public static final String DEFAULT_FONT = "Barlow";	
	
	
	/**
	 *  0: random seed, for reproducing a specific simulation (-1 to discard it);
	 *  1: number of bps in non-repetitive sequence (an int, so at most 2 billion);
	 *  2: ratio between total repetitive bps and total unique bps in the genome;  
	 *  3: directory containing the files that describe the repeat model;
   	 *  4: (output file) the simulated genome in FASTA format ("null" to discard it);
   	 *  5: (output file) random reads from the genome, sampled uniformly and with no error
   	 *     ("null" to discard it);
   	 *  6: desired coverage of the simulated genome, for producing the reads above.
	 *
	 *  7: (output file) the set of repeats created by the simulation, in FASTA format 
	 *     ("null" to discard it);
	 *  8: (output file) the evolutionary tree of simulated repeats, in DOT format ("null" 
	 *     to discard it);
	 *
	 *  9: (output file) an image of the genome (use "null" to discard it);
	 * 10: pixel width of the output image;
	 * 
	 * 11: (output file) the de Bruijn graph of the repeat annotation, in DOT format 
	 *     ("null" to discard it);
	 * 12: order k of the de Bruijn graph;
	 * 13: read length to be used by the de Bruijn graph; we assume to observe every 
	 *     substring of k and k+1 consecutive repeat annotations of total length <= this;
	 * 14: how to handle unique sequence in the alphabet of the DBG: 1=a unique sequence 
	 *     does not match any other unique sequence, not even one of the same length; 
	 *     0=unique sequences of the same length are assumed to match;
	 * 15: distance threshold for building the alphabet of the DBG.
	 */
	public static void main(String[] args) throws IOException {
		final long RANDOM_SEED = Long.parseLong(args[0]);
		final int N_UNIQUE_BPS = Integer.parseInt(args[1]);
		final double REPEAT_BPS_OVER_UNIQUE_BPS = Double.parseDouble(args[2]);
		final String REPEAT_MODEL_DIR = args[3];
		final String GENOME_SEQUENCE_FILE = args[4];
		final String READS_FILE = args[5];
		final int READS_COVERAGE = Integer.parseInt(args[6]);
		final String OUTPUT_FASTA = args[7];
		final String OUTPUT_DOT = args[8];
		final String OUTPUT_IMAGE = args[9];
		final int N_OUTPUT_COLUMNS = Integer.parseInt(args[10]);
		final String OUTPUT_DBG = args[11];
		final int DBG_K = Integer.parseInt(args[12]);
		final int DBG_READ_LENGTH = Integer.parseInt(args[13]);
		final boolean DBG_UNIQUE_MODE = Integer.parseInt(args[14])==1;
		final int DBG_DISTANCE_THRESHOLD = Integer.parseInt(args[15]);
		
		final int STRING_CAPACITY = 1000000;  // Arbitrary
		final long N_REPEAT_BPS = (long)(REPEAT_BPS_OVER_UNIQUE_BPS*N_UNIQUE_BPS);
		RepeatModel model = new RepeatModel(REPEAT_MODEL_DIR);
		StringBuilder sb = new StringBuilder(STRING_CAPACITY);
		Random random = RANDOM_SEED==-1?new Random():new Random(RANDOM_SEED);
		BufferedWriter bw;
		
		int i;
		long nRepeatBps;
		String parentString, repeatString;
		
		// Building the genome
		initializeGenome(N_UNIQUE_BPS,model,sb,random);
		nRepeatBps=0;
		while (nRepeatBps<N_REPEAT_BPS) {
			nRepeatBps+=epoch(model,sb,random);
			System.err.println("Epoch completed, "+nRepeatBps+" repeat bps, "+N_UNIQUE_BPS+" unique bps, "+lastRepeat+" distinct repeats.");
		}
		
		// Printing output
		if (!OUTPUT_FASTA.equalsIgnoreCase("null")) {
			bw = new BufferedWriter(new FileWriter(OUTPUT_FASTA));
			for (i=0; i<=lastRepeat; i++) bw.write(repeats[i].toString(i,model.minAlignmentLength));
			bw.close();
		}
		if (!OUTPUT_DOT.equalsIgnoreCase("null")) {
			bw = new BufferedWriter(new FileWriter(OUTPUT_DOT));
			bw.write("digraph G { \n");
			for (i=0; i<=lastRepeat; i++) {
				if (repeats[i].parent!=null) {
					parentString=repeats[i].parent.id;
					parentString=parentString.replace(' ','_');
					parentString=parentString.replace('/','_');
					parentString=parentString.replace(',','_');
					parentString=parentString.replace('-','_');
					repeatString=repeats[i].id;
					repeatString=repeatString.replace(' ','_');
					repeatString=repeatString.replace('/','_');
					repeatString=repeatString.replace(',','_');
					repeatString=repeatString.replace('-','_');
					bw.write(parentString+" -> "+repeatString+";  \n");
				}
				else {
					repeatString=repeats[i].id;
					repeatString=repeatString.replace(' ','_');
					repeatString=repeatString.replace('/','_');
					repeatString=repeatString.replace(',','_');
					repeatString=repeatString.replace('-','_');
					bw.write(repeatString+";  \n");
				}
			}
			bw.write("}");
			bw.close();
		}
		if (!OUTPUT_IMAGE.equalsIgnoreCase("null")) drawGenome(N_OUTPUT_COLUMNS,OUTPUT_IMAGE);
		if (!OUTPUT_DBG.equalsIgnoreCase("null")) buildDBG(DBG_K,DBG_READ_LENGTH,model.minAlignmentLength,DBG_UNIQUE_MODE,DBG_DISTANCE_THRESHOLD,OUTPUT_DBG,true);
		if (!GENOME_SEQUENCE_FILE.equalsIgnoreCase("null")) printGenome(GENOME_SEQUENCE_FILE);
		if (!READS_FILE.equalsIgnoreCase("null")) printReads(DBG_READ_LENGTH,READS_COVERAGE,READS_FILE,random);
	}
	
	
	
	
	
	
	
	
	// ------------------------------------ GENOME ---------------------------------------
	
	/**
	 * Sets the genome to a non-repetitive random sequence of length $nUniqueBps >=
	 * 2*minAlignmentLength$.
	 */
	private static final void initializeGenome(int nUniqueBps, RepeatModel model, StringBuilder sb, Random random) {
		final int CAPACITY_ROWS = 100;  // Arbitrary
		final int CAPACITY_COLUMNS = 1000;  // Arbitrary
		RepeatInstance instance;
		
		// Repeats
		unique = new Repeat();
		lastRepeat=-1;
		unique.initialize(nUniqueBps,model,sb,random);
		repeats = new Repeat[CAPACITY_ROWS];
		repeats[0]=unique;
		lastTree=0; unique.tree=0;
		
		// Instances
		root = new RepeatInstance();
		instance = new RepeatInstance(unique.sequence,unique,true,0,unique.sequenceLength-1);
		root.previous=null; root.next=instance; 
		instance.previous=root; instance.next=null;
		byTree = new RepeatInstance[CAPACITY_ROWS][0];
		byTree[0] = new RepeatInstance[CAPACITY_COLUMNS];
		byTree[0][0]=instance;
		lastByTree = new int[CAPACITY_ROWS];
		lastByTree[0]=0;
		cumulativeByTree = new int[CAPACITY_ROWS];
	}
	
	
	/**
	 * If Repbase strings were provided but they have all been used up by the simulation,
	 * the program reverts to the non-Repbase case.
	 *
	 * @return the number of new repetitive basepairs created by the procedure.
	 */
	private static final long epoch(RepeatModel model, StringBuilder sb, Random random) {
		final int MAX_NO_CONTRIBUTION_ITERATIONS = 10;  // Arbitrary
		final int CAPACITY = 1000;  // Arbitrary
		boolean isSatellite, isUsed;
		int i, j;
		int repeatTree, insertionTree, noContribution, id, length;
		long newBps, totalNewBps;
		Repeat repeat;
		RepeatInstance instance;
		
		// Creating a new repeat
		repeat = new Repeat();
		if (lastRepeat==0 || random.nextDouble()>model.fromInstanceProb) {
			// New repeat from scratch
			if (model.hasUnusedRepbaseStrings()) {
				do {
					isSatellite=model.getType(random)>=Constants.INTERVAL_PERIODIC;
					if (isSatellite) {
						id=random.nextInt(model.lastRepbaseSat+1);
						isUsed=model.repbaseSatIsUsed[id];
						model.repbaseSatIsUsed[id]=true;
					}
					else {
						id=random.nextInt(model.lastRepbase+1);
						isUsed=model.repbaseIsUsed[id];
						model.repbaseIsUsed[id]=true;
					}
				} while (isUsed);
				repeat.initialize(model,isSatellite,id,random);
			}
			else repeat.initialize(model,sb,random);
			repeatTree=-1;
		}
		else {
			// New repeat from one or more instances
			if (random.nextDouble()<=model.fromInstanceProb_multiInstance) {
				instance=getRandomInstance(model.minAlignmentLength<<1/*Arbitrary*/,random);
				repeat.initialize(instance,model,random);
				repeatTree=-1;
			}
			else {
				instance=getRandomInstance(Math.POSITIVE_INFINITY,random);
				repeat.initialize(instance,model,sb,random);
				repeatTree=instance.repeat.tree;
			}
		}
		// $lastRepeat$ is incremented by $Repeat.initialize$.
		if (lastRepeat==repeats.length) {
			Repeat[] newRepeats = new Repeat[repeats.length<<1];
			System.arraycopy(repeats,0,newRepeats,0,repeats.length);
			repeats=newRepeats;
		}
		repeats[lastRepeat]=repeat;
		if (repeatTree==-1) {
			lastTree++;
			if (lastTree==byTree.length) {
				RepeatInstance[][] newByTree = new RepeatInstance[byTree.length<<1][0];
				for (i=0; i<byTree.length; i++) newByTree[i]=byTree[i];
				byTree=newByTree;
				int[] newLastByTree = new int[byTree.length];
				System.arraycopy(lastByTree,0,newLastByTree,0,lastByTree.length);
				lastByTree=newLastByTree;
				int[] newCumulativeByTree = new int[byTree.length];
				System.arraycopy(cumulativeByTree,0,newCumulativeByTree,0,cumulativeByTree.length);
				cumulativeByTree=newCumulativeByTree;
			}	
			byTree[lastTree] = new RepeatInstance[CAPACITY];
			lastByTree[lastTree]=-1;
			repeatTree=lastTree;
		}
		repeat.tree=repeatTree;
		
		// Replication wave	
		totalNewBps=0; noContribution=0;
		for (i=0; i<repeat.frequency; i++) {
			newBps=epoch_impl(repeat,model,sb,random);
			if (newBps==0) {
				noContribution++;
				if (noContribution==MAX_NO_CONTRIBUTION_ITERATIONS || treesWithLongInstances()==0) break;
			}
			else {
				noContribution=0;
				totalNewBps+=newBps;
			}
			System.err.println("Replication wave "+i+" out of "+repeat.frequency);
		}
		return totalNewBps;
	}
	
	
	/**
	 * Creates a new instance of $repeat$ and inserts it into an existing instance of a
	 * repeat (possibly of the same repeat, or of the artificial repeat that represents 
	 * unique sequence).
	 *
	 * @return the number of new basepairs created, or zero if insertion did not succeed:
	 * this happens iff no instance of the destination repeat is long enough. In this case
	 * the procedure fails instead of trying with another destination repeat.
	 */
	private static final int epoch_impl(Repeat repeat, RepeatModel model, StringBuilder sb, Random random) {
		int i, j, p;
		int toTree;
		final int fromTree = repeat.tree;
		RepeatInstance fromInstance, toInstance, prefix, suffix;
		
		do { fromInstance=repeat.getInstance(model,sb,random); }
		while (fromInstance.sequenceLength<model.minAlignmentLength);
		if (lastTree<=1) {
			toTree=0;
			if (lastByTree[0]==-1) return 0;
			i=lastByTree[0]==0?0:random.nextInt(lastByTree[0]+1);
			lastByTree[fromTree]++;
			if (lastByTree[fromTree]==byTree[fromTree].length) {
				RepeatInstance[] newByTree = new RepeatInstance[byTree[fromTree].length<<1];
				System.arraycopy(byTree[fromTree],0,newByTree,0,byTree[fromTree].length);
				byTree[fromTree]=newByTree;
			}
			byTree[fromTree][lastByTree[fromTree]]=fromInstance;
		}
		else {
			toTree=Arrays.binarySearch(repeat.insertionProbCumulative,random.nextDouble());
			if (toTree<0) toTree=-1-toTree;
			if (lastByTree[toTree]==-1) return 0;
			i=lastByTree[toTree]==0?0:random.nextInt(lastByTree[toTree]+1);
			lastByTree[fromTree]++;
			if (lastByTree[fromTree]==byTree[fromTree].length) {
				RepeatInstance[] newByTree = new RepeatInstance[byTree[fromTree].length<<1];
				System.arraycopy(byTree[fromTree],0,newByTree,0,byTree[fromTree].length);
				byTree[fromTree]=newByTree;
			}
			byTree[fromTree][lastByTree[fromTree]]=fromInstance;
		}
		toInstance=byTree[toTree][i];
		p=model.minAlignmentLength+random.nextInt(Math.max(toInstance.sequenceLength-(model.minAlignmentLength<<1),1));
		prefix=toInstance.getPrefix(p);
		suffix=toInstance.getSuffix(toInstance.sequenceLength-p);
		prefix.previous=toInstance.previous;
		toInstance.previous.next=prefix;
		prefix.next=fromInstance;
		fromInstance.previous=prefix;
		fromInstance.next=suffix;
		suffix.previous=fromInstance;
		suffix.next=toInstance.next;
		if (toInstance.next!=null) toInstance.next.previous=suffix;
		if (prefix.sequenceLength>=model.minAlignmentLength<<1) {
			byTree[toTree][i]=prefix;
			if (suffix.sequenceLength>=model.minAlignmentLength<<1) {
				lastByTree[toTree]++;
				if (lastByTree[toTree]==byTree[toTree].length) {
					RepeatInstance[] newByTree = new RepeatInstance[byTree[toTree].length<<1];
					System.arraycopy(byTree[toTree],0,newByTree,0,byTree[toTree].length);
					byTree[toTree]=newByTree;
				}
				byTree[toTree][lastByTree[toTree]]=suffix;
			}
		}
		else if (suffix.sequenceLength>=model.minAlignmentLength<<1) byTree[toTree][i]=suffix;
		else {
			for (j=i; j<lastByTree[toTree]; j++) byTree[toTree][j]=byTree[toTree][j+1];
			lastByTree[toTree]--;
		}
		return fromInstance.sequenceLength;
	}
	
	
	/**
	 * @return an instance of length at most $maxLength$, chosen at random among those in 
	 * $byTree$.
	 */
	private static final RepeatInstance getRandomInstance(int maxLength, Random random) {
		int i;
		int repeatTree;
		RepeatInstance instance;
		
		do {
			cumulativeByTree[0]=0;
			cumulativeByTree[1]=lastByTree[1]+1;
			for (i=2; i<=lastTree; i++) cumulativeByTree[i]=cumulativeByTree[i-1]+lastByTree[i]+1;
			i=1+random.nextInt(cumulativeByTree[lastTree]);
			repeatTree=Arrays.binarySearch(cumulativeByTree,1,lastTree+1,i);
			if (repeatTree<0) repeatTree=-1-repeatTree;
			while (repeatTree>0 && lastByTree[repeatTree]==-1) repeatTree--;
			instance=byTree[repeatTree][i-cumulativeByTree[repeatTree-1]-1];
		}
		while (instance.sequenceLength>maxLength);
		return instance;
	}
	
	
	/**
	 * Returns a substring of length $length$ of the genome, grown alternatively from the 
	 * left and from the right side of $seed$. Any instance can be incorporated in the 
	 * process, including short and periodic.
	 *
	 * Remark: the procedure uses global arrays $sbLeft,sbRight,forReverse$.
	 *
	 * @param length assumed to be bigger than $seed$.
	 */
	private static final String getMultiInstanceString(RepeatInstance seed, int length) {
		boolean direction;
		int currentLength, delta;
		RepeatInstance instanceLeft, instanceRight;
		
		if (sbLeft==null) sbLeft = new StringBuilder();
		if (sbRight==null) sbRight = new StringBuilder();
		if (forReverse==null) forReverse = new StringBuilder();
		sbLeft.delete(0,sbLeft.length());
		sbRight.delete(0,sbRight.length());
		instanceLeft=seed; instanceRight=seed;
		currentLength=seed.sequenceLength;
		direction=false;
		while (currentLength<length && (instanceLeft!=null || instanceRight!=null)) {
			if (direction) {
				instanceRight=instanceRight.next;
				if (instanceRight!=null) {
					if (currentLength+instanceRight.sequenceLength<=length) {
						sbRight.append(instanceRight.sequence);
						currentLength+=instanceRight.sequenceLength;
					}
					else {
						sbRight.append(instanceRight.sequence.substring(0,length-currentLength));
						currentLength=length;
					}
				}
			}
			else {
				instanceLeft=instanceLeft.previous;
				if (instanceLeft!=null) {
					forReverse.delete(0,forReverse.length());
					if (currentLength+instanceLeft.sequenceLength<=length) {
						forReverse.append(instanceLeft.sequence);
						currentLength+=instanceLeft.sequenceLength;
					}
					else {
						forReverse.append(instanceLeft.sequence.substring(instanceLeft.sequenceLength-(length-currentLength)));
						currentLength=length;
					}
					forReverse.reverse();
					sbLeft.append(forReverse);
				}
			}
			direction=!direction;
		}
		return (sbLeft.length()>0?sbLeft.reverse().toString():"")+seed.sequence+(sbRight.length()>0?sbRight.toString():"");		
	}
	
	
	/**
	 * @return the number of nonempty rows in $lastByTree$.
	 */
	private static final int treesWithLongInstances() {
		int i, out;
		
		out=0;
		for (i=0; i<=lastTree; i++) {
			if (lastByTree[i]>=0) out++;
		}
		return out;
	}
	
	
	/**
	 * @param mode TRUE: length in basepairs; FALSE: length in number of intances.
	 */
	private static final long getGenomeLength(boolean mode) {
		long out;
		RepeatInstance currentInstance;
		
		out=0; currentInstance=root.next;
		if (mode) {
			do {
				out+=currentInstance.sequenceLength;
				currentInstance=currentInstance.next;
			} while (currentInstance!=null);
		}
		else {
			do {
				out++;
				currentInstance=currentInstance.next;
			} while (currentInstance!=null);
		}
		return out;
	}
	
	
	/**
	 * @param nColumns pixel width of the image (one pixel corresponds to one basepair).
	 */
	private static final void drawGenome(int nColumns, String outputFile) throws IOException {
		int x, y;
		int fromX, fromY;
		final long genomeLength = getGenomeLength(true);
		final long nInstances = getGenomeLength(false);
		final int nLines = 1+(int)(genomeLength/nColumns);
		final int nRows = nLines*ROWS_PER_LINE;
		RepeatInstance currentInstance;
		BufferedImage image;
		Graphics2D g2d;
		int[] tmpArray = new int[2];
		int[] pixels;
		
		image = new BufferedImage(nColumns,nRows,BufferedImage.TYPE_INT_RGB);
		pixels=((DataBufferInt)image.getRaster().getDataBuffer()).getData();
		Arrays.fill(pixels,Colors.COLOR_BACKGROUND);
		System.err.println("Drawing genome...  ("+nInstances+" blocks, "+genomeLength+" bps)");
		currentInstance=root;
		fromX=0; fromY=0;
		while (currentInstance.next!=null) {
			currentInstance=currentInstance.next;
			currentInstance.draw(fromX,fromY,image,nColumns,tmpArray);
			fromX=tmpArray[0]; fromY=tmpArray[1];
		}
		ImageIO.write(image,"png",new File(outputFile));
	}
	
	
	/**
	 * As a single string.
	 */
	private static final void printGenome(String outputFile) throws IOException {
		final long genomeLength = getGenomeLength(true);
		RepeatInstance currentInstance;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(outputFile));
		bw.write(">U0/1/0_"+(genomeLength-1)+"/ \n");
		currentInstance=root;
		while (currentInstance.next!=null) {
			currentInstance=currentInstance.next;
			bw.write(currentInstance.sequence);
		}
		bw.write("\n"); bw.close();
	}
	
	
	/**
	 * Randomly sampled, no error.
	 */
	private static final void printReads(int readLength, int coverage, String outputFile, Random random) throws IOException {
		int i, j;
		final int genomeLength = (int)getGenomeLength(true);
		final int nReads = (int)((coverage*genomeLength)/readLength);
		RepeatInstance currentInstance;
		StringBuilder sb;
		BufferedWriter bw;
		
		// Building the genome sequence
		sb = new StringBuilder(genomeLength);
		currentInstance=root;
		while (currentInstance.next!=null) {
			currentInstance=currentInstance.next;
			sb.append(currentInstance.sequence);
		}
		
		// Sampling reads
		bw = new BufferedWriter(new FileWriter(outputFile));
		for (i=0; i<nReads; i++) {
			j=random.nextInt(genomeLength-readLength+1);
			bw.write(">U0/"+i+"/0_"+readLength+"/ \n");
			bw.write(sb.substring(j,j+readLength));
			bw.write("\n");
		}
		bw.close();
	}


	
	
	
	

		
	// --------------------------- REPEATS AND REPEAT INSTANCES --------------------------
	
	private static final void reverseComplement(StringBuilder sb) {
		char c, d;
		int i;
		final int length = sb.length();
		final int halfLength = length/2;

		for (i=0; i<halfLength; i++) {
			c=sb.charAt(i); d=sb.charAt(length-1-i);
			switch (d) {
				case 'a': sb.setCharAt(i,'t'); break;
				case 'c': sb.setCharAt(i,'g'); break;
				case 'g': sb.setCharAt(i,'c'); break;
				case 't': sb.setCharAt(i,'a'); break;
			}
			switch (c) {
				case 'a': sb.setCharAt(length-1-i,'t'); break;
				case 'c': sb.setCharAt(length-1-i,'g'); break;
				case 'g': sb.setCharAt(length-1-i,'c'); break;
				case 't': sb.setCharAt(length-1-i,'a'); break;
			}
		}
		if (length%2!=0) {
			i=halfLength;
			switch (sb.charAt(i)) {
				case 'a': sb.setCharAt(i,'t'); break;
				case 'c': sb.setCharAt(i,'g'); break;
				case 'g': sb.setCharAt(i,'c'); break;
				case 't': sb.setCharAt(i,'a'); break;
			}
		}
	}
	
	
	public static class RepeatInstance implements Comparable {
		public static final byte ORDER_ID = 0;
		public static final byte ORDER_FEATURES = 1;
		public static byte order;
		
		public int id, characterID;
		public String sequence;
		public int sequenceLength;
		public Repeat repeat;
		public boolean orientation;
		public RepeatInstance previous, next;  // In genome order
		
		/**
		 * The substring of the repeat of which this instance is a copy. 
		 * $repeatStart <= repeatEnd$ regardless of $orientation$.
		 * Not used if the repeat type is single deletion or periodic.
		 */
		public int repeatStart, repeatEnd;
		
		
		public RepeatInstance() {
			sequence=null;
			sequenceLength=0;
			repeat=null;
			orientation=false;
			repeatStart=-1; repeatEnd=-1;
		}
		
		
		public RepeatInstance(String s, Repeat r, boolean o, int rs, int re) {
			sequence=s;
			sequenceLength=sequence.length();
			repeat=r;
			orientation=o;
			repeatStart=rs;
			repeatEnd=re;
		}
		
		
		/**
		 * Makes this instance a clone of $otherInstance$, excluding fields $previous,
		 * next$.
		 *
		 * @param mode TRUE: simple clone; FALSE: this instance becomes an instance of
		 * the same length of the artificial unique-sequence repeat.
		 */
		public RepeatInstance(RepeatInstance otherInstance, boolean mode) {
			if (mode) {
				this.sequence=otherInstance.sequence;
				this.sequenceLength=otherInstance.sequenceLength;
				this.repeat=otherInstance.repeat;
				this.orientation=otherInstance.orientation;
				this.previous=null; this.next=null;
				this.repeatStart=otherInstance.repeatStart;
				this.repeatEnd=otherInstance.repeatEnd;
			}
			else {
				this.sequence=otherInstance.sequence;
				this.sequenceLength=otherInstance.sequenceLength;
				this.repeat=unique;
				this.orientation=true;
				this.previous=null; this.next=null;
				this.repeatStart=-1;
				this.repeatEnd=-1;
			}
		}
		
		
		public RepeatInstance getPrefix(int prefixLength) {
			RepeatInstance out = new RepeatInstance();
			out.sequence=sequence.substring(0,prefixLength);
			out.sequenceLength=prefixLength;
			out.repeat=repeat;
			out.orientation=orientation;
			if (orientation) {
				out.repeatStart=repeatStart;
				out.repeatEnd=repeatStart+prefixLength-1;
			}
			else {
				out.repeatEnd=repeatEnd;
				out.repeatStart=repeatStart+sequenceLength-prefixLength-1;
			}
			out.previous=previous;
			out.next=null;
			return out;
		}
		
		
		public RepeatInstance getSuffix(int suffixLength) {
			RepeatInstance out = new RepeatInstance();
			out.sequence=sequence.substring(sequenceLength-suffixLength);
			out.sequenceLength=suffixLength;
			out.repeat=repeat;
			out.orientation=orientation;
			if (orientation) {
				out.repeatStart=repeatEnd-suffixLength+1;
				out.repeatEnd=repeatEnd;
			}
			else {
				out.repeatEnd=repeatEnd-(sequenceLength-suffixLength);
				out.repeatStart=repeatStart;
			}
			out.previous=null;
			out.next=next;
			return out;
		}
		
		
		/**
		 * Draws the repeat instance starting from point $(fromX,toX)$ of $image$, and 
		 * handling overflows by continuing to draw on the following rows.
		 *
		 * @param nColumns total number of columns in $image$;
		 * @return out output array: 0=the new value of fromX; 1=the new value of fromY
		 * after the procedure completes.
		 */
		public final void draw(int fromX, int fromY, BufferedImage image, int nColumns, int[] out) throws IOException {
			boolean isOverflow;
			int overflow, missingLength, windowStart, windowEnd;
			int repeatStartPrime, repeatEndPrime;
			
			if (repeat.type==-1) {
				missingLength=sequenceLength;
				while (true) {
					overflow=Math.max(fromX+missingLength-nColumns,0);
					if (overflow>0) {
						fromX=0; fromY+=ROWS_PER_LINE;
						missingLength=overflow;
					}
					else {
						fromX+=missingLength;
						break;
					}
				}
				out[0]=fromX; out[1]=fromY;
			}
			else if (repeat.type<Constants.INTERVAL_PERIODIC) {
				if (repeat.type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
					// Arbitrary choice, just for display purposes.
					repeatStartPrime=(repeat.sequenceLength-sequenceLength)>>1; 
					repeatEndPrime=repeat.sequenceLength-repeatStartPrime;
				}
				else { repeatStartPrime=repeatStart; repeatEndPrime=repeatEnd; }
				missingLength=sequenceLength; windowStart=repeatStartPrime; isOverflow=false;
				while (true) {
					overflow=Math.max(fromX+missingLength-nColumns,0);
					windowEnd=repeatEndPrime-overflow;
					drawSubstring(repeat,windowStart,windowEnd,orientation,image,fromX,fromY,isOverflow);
					if (overflow>0) {
						fromX=0; fromY+=ROWS_PER_LINE; isOverflow=true;
						missingLength-=windowEnd-windowStart+1;
						windowStart=windowEnd+1;
					}
					else {
						fromX+=windowEnd-windowStart+1; isOverflow=false;
						missingLength=0;
						break;
					}
				}
				out[0]=fromX; out[1]=fromY;
			}
			else {
				missingLength=sequenceLength; isOverflow=false;
				while (true) {
					overflow=Math.max(fromX+missingLength-nColumns,0);
					if (overflow>0) {
						drawPeriodicSubstring(repeat,sequenceLength,orientation,image,fromX,nColumns-1,fromY,isOverflow,true);
						fromX=0; fromY+=ROWS_PER_LINE; isOverflow=true;
						missingLength=overflow;
					}
					else {
						drawPeriodicSubstring(repeat,sequenceLength,orientation,image,fromX,fromX+missingLength-1,fromY,isOverflow,false);
						fromX+=missingLength;
						break;
					}
				}
				out[0]=fromX; out[1]=fromY;
			}
		}
		
		
		public String toString() {
			return repeat.type+">"+repeat.id+"["+repeatStart+".."+repeatEnd+"](length="+sequenceLength+"), orient="+orientation;
		}
		
		
		/**
		 * $ORDER_FEATURES$ has the following meaning:
		 *
		 * For non-periodic and non-unique: 
		 * repeat.id, repeatStart, repeatEnd.
		 *
		 * For all other types: 
		 * repeat.id, sequenceLength.
		 *
		 * Remark: orientation is not used on purpose.
		 */
		public int compareTo(Object other) {
			RepeatInstance otherInstance = (RepeatInstance)other;
			
			if (order==ORDER_ID) {
				if (id<otherInstance.id) return -1;
				else if (id>otherInstance.id) return 1;
			}
			else if (order==ORDER_FEATURES) {
				int i = repeat.id.compareTo(otherInstance.repeat.id);
				if (i<0) return -1;
				else if (i>0) return 1;
				if (repeat.type>=0 && repeat.type<=Constants.INTERVAL_DENSE_SINGLEDELETION) {
					if (repeatStart<otherInstance.repeatStart) return -1;
					else if (repeatStart>otherInstance.repeatStart) return 1;
					if (repeatEnd<otherInstance.repeatEnd) return -1;
					else if (repeatEnd>otherInstance.repeatEnd) return 1;
				}
				else {
					if (sequenceLength<otherInstance.sequenceLength) return -1;
					else if (sequenceLength>otherInstance.sequenceLength) return 1;
				}
			}
			return 0;
		}
		
		
		/**
		 * Similar to $compareTo()$.
		 * 
		 * @param uniqueMode TRUE: a unique sequence does not match any other unique
		 * sequence, not even one of the same length; FALSE: unique sequences of the same
		 * length are assumed to match.
		 */
		public boolean isApproximatelyIdentical(RepeatInstance otherInstance, int distanceThreshold, boolean uniqueMode) {
			if (!repeat.id.equalsIgnoreCase(otherInstance.repeat.id)) return false;
			if (repeat.type==-1) {
				if (uniqueMode) return false;
				else return Math.abs(sequenceLength,otherInstance.sequenceLength)<=distanceThreshold<<1;
			}
			else if (repeat.type<=Constants.INTERVAL_DENSE_SINGLEDELETION) {
				return Math.abs(repeatStart,otherInstance.repeatStart)<=distanceThreshold && 
					   Math.abs(repeatEnd,otherInstance.repeatEnd)<=distanceThreshold;
			}
			else return Math.abs(sequenceLength,otherInstance.sequenceLength)<=distanceThreshold<<1;
		}
	}
	
	
	/**
	 * Draws the substring $[repeatStart..repeatEnd]$ of a non-periodic repeat, 
	 * starting from point $(fromX,toX)$ of $image$, assuming that the drawing does 
	 * not overflow the image.
	 *
	 * @param repeatLength length of the full repeat;
	 * @param isOverflow TRUE iff the substring is the result of splitting a repeat
	 * instance because it overflowed; in this case, left-truncation information is 
	 * not printed. 
	 */
	private static final void drawSubstring(Repeat repeat, int repeatStart, int repeatEnd, boolean orientation, BufferedImage image, int fromX, int fromY, boolean isOverflow) throws IOException {
		int x, y;
		int height, firstX, lastX, color;
		final int sequenceLength = repeatEnd-repeatStart+1;
		final int toX = fromX+sequenceLength-1;
		final int leftTruncation = orientation?repeatStart:repeat.sequenceLength-repeatEnd;
		final int rightTruncation = orientation?repeat.sequenceLength-repeatEnd:repeatStart;
		final int colorPlain = Colors.type2color(repeat.type);
		final int colorHighlight = Colors.COLOR_TEXT;
		String label;
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
		g2d.setColor(new Color(colorHighlight));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
		final NumberFormat formatter = NumberFormat.getInstance();
		formatter.setMaximumFractionDigits(2);
		
		// Drawing truncated triangle
		firstX=fromX-(orientation?repeatStart:repeat.sequenceLength-repeatEnd);
		lastX=fromX+sequenceLength+(orientation?repeat.sequenceLength-repeatEnd:repeatStart);
		for (x=fromX; x<=toX; x++) {
			height=Math.round(ROWS_PER_LINE_OVER_TWO*(orientation?lastX-x:x-firstX)/repeat.sequenceLength);
			if ( (!isOverflow && leftTruncation>TRUNCATION_THRESHOLD && (x-fromX)<=TRUNCATION_TAG_WIDTH_PIXELS) || 
				 (rightTruncation>TRUNCATION_THRESHOLD && (fromX+sequenceLength-x)<=TRUNCATION_TAG_WIDTH_PIXELS)
			   ) color=colorHighlight;
			else color=colorPlain;
			for (y=fromY+ROWS_PER_LINE_OVER_TWO-height; y<=fromY+ROWS_PER_LINE_OVER_TWO+height; y++) Colors.setRGB(image,x,y,color);
		}
		
		// Drawing text
		label="ID="+repeat.id+", T="+Colors.type2string(repeat.type)+", L="+sequenceLength;
		x=(fromX+toX)/2-Math.round(fontMetrics.stringWidth(label),2);
		y=fromY+ROWS_PER_LINE_OVER_TWO+fontMetrics.getAscent()/2;
		g2d.drawString(label,x,y);
		if (!isOverflow && leftTruncation>TRUNCATION_THRESHOLD) {
			label=""+formatter.format(((double)leftTruncation)/TRUNCATION_MULTIPLE);
			x=fromX+TRUNCATION_TAG_WIDTH_PIXELS+TRUNCATION_TAG_WIDTH_PIXELS;
			y=fromY+ROWS_PER_LINE_OVER_TWO+fontMetrics.getAscent()/2;
			g2d.drawString(label,x,y);
		}
		if (rightTruncation>TRUNCATION_THRESHOLD) {
			label=""+formatter.format(((double)rightTruncation)/TRUNCATION_MULTIPLE);
			x=toX-TRUNCATION_TAG_WIDTH_PIXELS-TRUNCATION_TAG_WIDTH_PIXELS-fontMetrics.stringWidth(label);
			y=fromY+ROWS_PER_LINE_OVER_TWO+fontMetrics.getAscent()/2;
			g2d.drawString(label,x,y);
		}
	}
	
	
	/**
	 * Draws a periodic repeat over points $(fromX..toX,fromY..)$ of $image$, assuming 
	 * that the drawing does not overflow the image.
	 */
	private static final void drawPeriodicSubstring(Repeat repeat, int sequenceLength, boolean orientation, BufferedImage image, int fromX, int toX, int fromY, boolean isOverflow, boolean willOverflow) {
		int x, y;
		final int colorPlain = Colors.type2color(repeat.type);
		final int colorHighlight = Colors.COLOR_TEXT;
		int color;
		String label;
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
		g2d.setColor(new Color(colorHighlight));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
	
		// Drawing rectangle
		for (x=fromX; x<=toX; x++) {
			if ( (!isOverflow && (x-fromX)<=TRUNCATION_TAG_WIDTH_PIXELS) || 
				 (!willOverflow && toX-x<=TRUNCATION_TAG_WIDTH_PIXELS)
			   ) color=colorHighlight;
			else color=colorPlain;
			for (y=fromY; y<=fromY+ROWS_PER_LINE; y++) Colors.setRGB(image,x,y,color);
		}
	
		// Drawing text
		label="ID="+repeat.id+", O="+(orientation?"FWD":"RC")+", P="+repeat.sequenceLength+", L="+sequenceLength;
		x=(fromX+toX)/2-Math.round(fontMetrics.stringWidth(label),2);
		y=fromY+ROWS_PER_LINE_OVER_TWO+fontMetrics.getAscent()/2;
		g2d.drawString(label,x,y);
	}
	
	
	private static class Repeat {
		public int tree;
		public String id;
		public Repeat parent;
		public int type;  // -1 for unique sequence
		public String sequence;  // Of the period if periodic
		public int sequenceLength;
		public double[] insertionProb;  // In every tree of past repeats
		public double[] insertionProbCumulative;
		public int frequency;
		
		
		/**
		 * Creates a new unique region of given $length$ according to $model$.
		 *
		 * Remark: the procedure does not set field $tree$.
		 */
		public final void initialize(int length, RepeatModel model, StringBuilder sb, Random random) {
			id=(++lastRepeat)+"";
			parent=null;
			type=-1;
			sequenceLength=length;
			model.getUniqueString(sequenceLength,sb,random);
			sequence=sb.toString();
			insertionProb=null;
			insertionProbCumulative=null;
			frequency=1;  // Never repeats
		}
		
		
		/**
		 * Creates a new random repeat from scratch, according to $model$.
		 *
		 * Remark: the procedure does not set field $tree$.
		 */											  
		public final void initialize(RepeatModel model, StringBuilder sb, Random random) {
			int i;
			
			id=(++lastRepeat)+"";
			parent=null;
			type=model.getType(random);
			if (type>=Constants.INTERVAL_DENSE_PREFIX && type<=Constants.INTERVAL_DENSE_SUBSTRING) {
				do { sequenceLength=model.getLength(random); }
				while (sequenceLength<model.minAlignmentLength<<1);
				model.getString(sequenceLength,sb,random);
				sequence=sb.toString();
			}
			else if (type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
				do { sequenceLength=model.getLength(random); }
				while (sequenceLength<DELTA+(model.minAlignmentLength<<1));
				model.getString(sequenceLength,sb,random);
				sequence=sb.toString();
			}
			else if (type==Constants.INTERVAL_PERIODIC) {
				sequenceLength=model.getShortPeriod(random);
				model.getShortPeriodString(sequenceLength,sb,random);
				sequence=sb.toString();
			}
			else {
				sequenceLength=model.getLength(random);
				model.getString(sequenceLength,sb,random);
				sequence=sb.toString();
			}
			insertionProb=Math.sampleFromSimplex(lastTree+1,random);
			insertionProbCumulative = new double[lastTree+1];
			insertionProbCumulative[0]=insertionProb[0];
			for (i=1; i<=lastTree; i++) insertionProbCumulative[i]=insertionProb[i]+insertionProbCumulative[i-1];
			frequency=type==Constants.INTERVAL_PERIODIC?model.getFrequencyPeriodic(random):model.getFrequency(random);
		}
		
		
		/**
		 * Creates a new repeat from an instance of a past one. The new repeat inherits 
		 * the type of the past repeat (if length allows), a perturbed version of its 
		 * insertion preferences, and the exact sequence of the instance, but it gets a 
		 * new frequency.
		 *
		 * Remark: the procedure does not set field $tree$.
		 *
		 * @param instance assumed to be of length at least $minAlignmentLength$.
		 */
		public final void initialize(RepeatInstance instance, RepeatModel model, StringBuilder sb, Random random) {
			int i, period;
			
			id=(++lastRepeat)+"";
			parent=instance.repeat;
			if (parent.type<Constants.INTERVAL_PERIODIC) {
				sequence=""+instance.sequence;
				sequenceLength=instance.sequenceLength;
				if (parent.type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
					if (sequenceLength>=DELTA+(model.minAlignmentLength<<1)) type=parent.type;
					else if (sequenceLength>=model.minAlignmentLength<<1) type=Constants.INTERVAL_DENSE_SUBSTRING;
					else type=Constants.INTERVAL_ALIGNMENT;
				}
				else if (parent.type>=Constants.INTERVAL_DENSE_PREFIX && parent.type<=Constants.INTERVAL_DENSE_SUBSTRING) {
					if (sequenceLength>=model.minAlignmentLength<<1) type=parent.type;
					else type=Constants.INTERVAL_ALIGNMENT;
				}
				else type=parent.type;
			}
			else {
				type=parent.type;
				if (instance.sequenceLength<=parent.sequenceLength) {
					period=instance.sequenceLength;
					sequence=instance.sequence;
					sequenceLength=period;
				}
				else {
					do { period=model.perturbPeriod(parent.sequenceLength,random); }
					while (period>instance.sequenceLength);
					i=random.nextInt(instance.sequenceLength-period);
					sequence=instance.sequence.substring(i,i+period);
					sequenceLength=period;
				}
			}
			insertionProb=model.perturbInsertionProb(parent.insertionProb,random);
			insertionProbCumulative = new double[insertionProb.length];
			insertionProbCumulative[0]=insertionProb[0];
			for (i=1; i<insertionProb.length; i++) insertionProbCumulative[i]=insertionProb[i]+insertionProbCumulative[i-1];
			frequency=model.getFrequency(random);
			System.err.println("Created a new repeat from an instance of element "+parent.id);
		}
		
		
		/**
		 * Creates a new non-periodic repeat, by taking the concatenation of several 
		 * existing instances, starting from $instance$ and growing it. The new repeat has
		 * a random non-periodic type and new insertion preferences and frequency.
		 *
		 * Remark: the procedure does not set field $tree$.
		 *
		 * @param instance assumed to be of length $>=minAlignmentLength$.
		 */
		public final void initialize(RepeatInstance instance, RepeatModel model, Random random) {
			final int minLength = instance.sequenceLength+model.minAlignmentLength;
			int i;
			
			id=(++lastRepeat)+"";
			parent=null;
			while (true) {
				type=model.getType(random);
				if (type==Constants.INTERVAL_ALIGNMENT) {
					do { sequenceLength=model.getLength(random); }
					while (sequenceLength<minLength);
					sequence=getMultiInstanceString(instance,sequenceLength);
					break;
				}
				else if (type>=Constants.INTERVAL_DENSE_PREFIX && type<=Constants.INTERVAL_DENSE_SUBSTRING) {
					do { sequenceLength=model.getLength(random); }
					while (sequenceLength<Math.max(model.minAlignmentLength<<1,minLength));
					sequence=getMultiInstanceString(instance,sequenceLength);
					break;
				}
				else if (type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
					do { sequenceLength=model.getLength(random); }
					while (sequenceLength<Math.max(DELTA+(model.minAlignmentLength<<1),minLength));
					sequence=getMultiInstanceString(instance,sequenceLength);
					break;
				}
			}
			insertionProb=Math.sampleFromSimplex(lastTree+1,random);
			insertionProbCumulative = new double[lastTree+1];
			insertionProbCumulative[0]=insertionProb[0];
			for (i=1; i<=lastTree; i++) insertionProbCumulative[i]=insertionProb[i]+insertionProbCumulative[i-1];
			frequency=model.getFrequency(random);
			System.err.println("Created a new multi-instance repeat of length "+sequenceLength+" from seed instance "+instance);
		}
		
		
		/**
		 * Creates a new random repeat from a Repbase entry.
		 *
		 * Remark: the procedure does not set field $tree$.
		 */											  
		public final void initialize(RepeatModel model, boolean isSatellite, int repbaseID, Random random) {
			int i;
			
			lastRepeat++;
			id=isSatellite?model.repbaseSatIDs[repbaseID]:model.repbaseIDs[repbaseID];
			parent=null;
			if (isSatellite) {
				sequence=model.repbaseSat[repbaseID];
				sequenceLength=sequence.length();
				type=Constants.INTERVAL_PERIODIC;
			}
			else {
				sequence=model.repbase[repbaseID];
				sequenceLength=sequence.length();
				do { type=model.getType(random); }
				while ( type==Constants.INTERVAL_PERIODIC || 
					    (type==Constants.INTERVAL_DENSE_SINGLEDELETION && sequenceLength<DELTA+(model.minAlignmentLength<<1)) ||
						(type>=Constants.INTERVAL_DENSE_PREFIX && type<=Constants.INTERVAL_DENSE_SUBSTRING && sequenceLength<model.minAlignmentLength<<1)
					  );
				
			}
			insertionProb=Math.sampleFromSimplex(lastTree+1,random);
			insertionProbCumulative = new double[lastTree+1];
			insertionProbCumulative[0]=insertionProb[0];
			for (i=1; i<=lastTree; i++) insertionProbCumulative[i]=insertionProb[i]+insertionProbCumulative[i-1];
			frequency=type==Constants.INTERVAL_PERIODIC?model.getFrequencyPeriodic(random):model.getFrequency(random);
			System.err.println("Created Repbase repeat "+id);
		}
		
		
		/**
		 * Remark: (1) the prefix-suffix type has equal prob. of generating a prefix or a 
		 * suffix, and uniform prob. of generating every length >=minAlignmentLength;
		 * (2) substring type has uniform prob. of generating every substring length >=
		 * minAlignmentLength, and uniform prob. of first position; (3) single-deletion
		 * ensures that a prefix and a suffix of length >=minAlignmentLength are always 
		 * present; (4) short-period can have a perturbed phase; (5) long-period can have 
		 * an integer multiple of the period, always in phase.
		 */
		public final RepeatInstance getInstance(RepeatModel model, StringBuilder sb, Random random) {
			boolean orientation;
			int i, p;
			int length, phase, repeatStart, repeatEnd;
			String prefix, suffix;
			
			repeatStart=-1; repeatEnd=-1;
			sb.delete(0,sb.length());
			if (type==Constants.INTERVAL_ALIGNMENT) {
				repeatStart=0; repeatEnd=sequenceLength-1;
				sb.append(sequence);
			}
			else if (type==Constants.INTERVAL_DENSE_PREFIX) {
				length=model.minAlignmentLength+random.nextInt(sequence.length()-model.minAlignmentLength);
				repeatStart=0; repeatEnd=length-1;
				sb.append(sequence.substring(repeatStart,repeatEnd+1));
			}
			else if (type==Constants.INTERVAL_DENSE_SUFFIX) {
				length=model.minAlignmentLength+random.nextInt(sequence.length()-model.minAlignmentLength);
				repeatStart=sequence.length()-length;
				repeatEnd=sequence.length()-1;
				sb.append(sequence.substring(repeatStart,repeatEnd+1));
			}
			else if (type==Constants.INTERVAL_DENSE_PREFIXSUFFIX) {
				length=model.minAlignmentLength+random.nextInt(sequence.length()-model.minAlignmentLength);
				if (random.nextBoolean()) {
					repeatStart=0; repeatEnd=length-1;
				}
				else {
					repeatStart=sequence.length()-length;
					repeatEnd=sequence.length()-1;
				}
				sb.append(sequence.substring(repeatStart,repeatEnd+1));
			}
			else if (type==Constants.INTERVAL_DENSE_SUBSTRING) {
				length=model.minAlignmentLength+random.nextInt(sequence.length()-model.minAlignmentLength);
				repeatStart=random.nextInt(sequence.length()-length);
				repeatEnd=repeatStart+length-1;
				sb.append(sequence.substring(repeatStart,repeatEnd+1));
			}
			else if (type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
				length=random.nextInt(sequence.length()-(model.minAlignmentLength<<1));							
				p=model.minAlignmentLength+random.nextInt((int)Math.max(sequence.length()-(model.minAlignmentLength<<1)-length,1));
				sb.append(sequence.substring(0,p));
				sb.append(sequence.substring(p+length));
				repeatStart=p; repeatEnd=p+length;
			}
			else if (type==Constants.INTERVAL_PERIODIC) {
				length=model.getLengthShortPeriod(random);
				phase=random.nextInt((int)(model.maxPhaseDifference*sequence.length()));
				suffix=sequence.substring(phase);
				prefix=sequence.substring(0,sequence.length()-phase);
				do { sb.append(suffix); sb.append(prefix); }
				while (sb.length()<length);
				if (sb.length()>length) sb.delete(length,sb.length());
				repeatStart=0; repeatEnd=sequenceLength-1;
			}
			else if (type==Constants.INTERVAL_PERIODIC+1) {
				length=model.getLengthLongPeriod(random);
				for (i=0; i<length; i++) sb.append(sequence);
				repeatStart=0; repeatEnd=sequenceLength-1;
			}
			if (random.nextDouble()<=model.rcProbability) {
				orientation=false;
				reverseComplement(sb);
			}
			else orientation=true;
			model.addReplicationNoise(sb,random);
			return new RepeatInstance(sb.toString(),this,orientation,repeatStart,repeatEnd);
		}
		
		
		/**
		 * FASTA format. Replicates the sequence if it is shorter than 
		 * $minAlignmentLength$.
		 */
		public String toString(int newID, int minAlignmentLength) {
			final int length = insertionProb==null?0:insertionProb.length;
			int i;
			
			String out = ">U0/"+newID+"/0_"+(sequenceLength-1)+"/"+id+"/"+type+"/"+frequency+"/";
			if (length>0) out+="P(0)="+insertionProb[0];
			for (i=1; i<length; i++) out+=",P("+i+")="+insertionProb[i];
			out+="\n";
			i=0;
			do { out+=sequence; i+=sequenceLength; } while (i<minAlignmentLength);
			out+="\n";
			return out;
		}
	}
	
	
	/**
	 * Interface to all the properties of repeats and of repeat instances that the user 
	 * specifies in several config files.
	 */
	private static class RepeatModel {
		public int minAlignmentLength;
		private int minShortPeriod, maxShortPeriod, shortPeriodRange;
		
		/**
		 * Probability of generating a type. Types are as in $Constants.java$, but 7 
		 * means long-period instead of unknown.
		 */
		private double[] typeProb, typeProbCumulative;
		
		/**
		 * Probability of creating a character
		 */
		private double[] charProb, charProbCumulative;
		private double[] charProbShortPeriod, charProbShortPeriodCumulative;
		private double[] charProbUnique, charProbUniqueCumulative;
		
		/**
		 * Probability of creating a repeat length
		 */
		private double[] lengthProb;  // Prob. of block i*mAL+1..(i+1)*mAL
		private double[] lengthProbCumulative;
		private double[] lengthProbShortPeriod;  // Prob. of block i*mAL+1..(i+1)*mAL
		private double[] lengthProbShortPeriodCumulative;
		private double[] lengthProbLongPeriod;  // Prob. of a tandem of i instances
		private double[] lengthProbLongPeriodCumulative;
		
		/**
		 * Probability of creating a new repeat from an existing instance
		 */
		private double fromInstanceProb, fromInstanceProb_multiInstance;
		private double maxPeriodDifference, maxPhaseDifference;
		private double insertionProbPerturbation;
		
		/**
		 * Probability of editing a sequence
		 */
		private double rcProbability;
		private double mismatchProb, insertionProb, deletionProb;
		private double mPlusI, mPlusIPlusD;
		private int maxIndelLength;
		
		/**
		 * Probability of a repeat having a given frequency in the genome.
		 */
		private double[] frequencyProb;  // Prob. of block 10^i+1..10^(i+1)
		private double[] frequencyProbCumulative;
		private double[] frequencyProbPeriodic;  // Prob. of frequency $i$.
		private double[] frequencyProbPeriodicCumulative;
		
		/**
		 * Repbase file
	     */
		private String[] repbase, repbaseSat;
		private int lastRepbase, lastRepbaseSat;
		private String[] repbaseIDs, repbaseSatIDs;  // IDs in the file
		
		/**
		 * TRUE iff a string has been already used to create a repeat.
		 */
		private boolean[] repbaseIsUsed, repbaseSatIsUsed;
		
		
		public RepeatModel(String configDir) throws IOException {
			int i, p;
			int length;
			String str, repbaseFile;
			BufferedReader br;
			String[] tokens;
			
			br = new BufferedReader(new FileReader(configDir+"/repeatModel.txt"));
			str=br.readLine();
			p=str.indexOf("/");
			minAlignmentLength=Integer.parseInt(p>=0?str.substring(0,p).trim():str.trim());
			DELTA=minAlignmentLength>>2;  // Arbitrary
			str=br.readLine();
			p=str.indexOf("/");
			minShortPeriod=Integer.parseInt(p>=0?str.substring(0,p).trim():str.trim());
			str=br.readLine();
			p=str.indexOf("/");
			maxShortPeriod=Integer.parseInt(p>=0?str.substring(0,p).trim():str.trim());
			shortPeriodRange=maxShortPeriod-minShortPeriod+1;
			
			// Loading $typeProb*$.
			typeProb = new double[N_REPEAT_TYPES];
			for (i=0; i<N_REPEAT_TYPES; i++) {
				str=br.readLine();
				p=str.indexOf("/");
				typeProb[i]=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			}
			typeProbCumulative = new double[N_REPEAT_TYPES];
			typeProbCumulative[0]=typeProb[0];
			for (i=1; i<N_REPEAT_TYPES; i++) typeProbCumulative[i]=typeProb[i]+typeProbCumulative[i-1];
			
			// Loading $charProb*$.
			charProb = new double[DNA_ALPHABET.length];
			for (i=0; i<DNA_ALPHABET.length; i++) {
				str=br.readLine();
				p=str.indexOf("/");
				charProb[i]=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			}
			charProbCumulative = new double[DNA_ALPHABET.length];
			charProbCumulative[0]=charProb[0];
			for (i=1; i<DNA_ALPHABET.length; i++) charProbCumulative[i]=charProb[i]+charProbCumulative[i-1];
			
			// Loading $charProbShortPeriod*$.
			charProbShortPeriod = new double[DNA_ALPHABET.length];
			for (i=0; i<DNA_ALPHABET.length; i++) {
				str=br.readLine();
				p=str.indexOf("/");
				charProbShortPeriod[i]=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			}
			charProbShortPeriodCumulative = new double[DNA_ALPHABET.length];
			charProbShortPeriodCumulative[0]=charProbShortPeriod[0];
			for (i=1; i<DNA_ALPHABET.length; i++) charProbShortPeriodCumulative[i]=charProbShortPeriod[i]+charProbShortPeriodCumulative[i-1];
			
			// Loading $charProbUnique*$.
			charProbUnique = new double[DNA_ALPHABET.length];
			for (i=0; i<DNA_ALPHABET.length; i++) {
				str=br.readLine();
				p=str.indexOf("/");
				charProbUnique[i]=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			}
			charProbUniqueCumulative = new double[DNA_ALPHABET.length];
			charProbUniqueCumulative[0]=charProbUnique[0];
			for (i=1; i<DNA_ALPHABET.length; i++) charProbUniqueCumulative[i]=charProbUnique[i]+charProbUniqueCumulative[i-1];
			
			// Loading from-instance probabilities
			str=br.readLine();
			p=str.indexOf("/");
			fromInstanceProb=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			str=br.readLine();
			p=str.indexOf("/");
			fromInstanceProb_multiInstance=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			str=br.readLine();
			p=str.indexOf("/");
			maxPeriodDifference=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			
			// Loading $maxPhaseDifference$.
			str=br.readLine();
			p=str.indexOf("/");
			maxPhaseDifference=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			
			// Loading $insertionProbPerturbation$.
			str=br.readLine();
			p=str.indexOf("/");
			insertionProbPerturbation=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			
			// Loading $rcProbability$.
			str=br.readLine();
			p=str.indexOf("/");
			rcProbability=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			
			// Loading noise probabilities.
			str=br.readLine();
			p=str.indexOf("/");
			mismatchProb=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			str=br.readLine();
			p=str.indexOf("/");
			insertionProb=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			mPlusI=mismatchProb+insertionProb;
			str=br.readLine();
			p=str.indexOf("/");
			deletionProb=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			mPlusIPlusD=mPlusI+deletionProb;
			str=br.readLine();
			p=str.indexOf("/");
			maxIndelLength=Integer.parseInt(p>=0?str.substring(0,p).trim():str.trim());
			
			br.close();
			
			// Loading $lengthProb*$.
			br = new BufferedReader(new FileReader(configDir+"/lengthProb.txt"));
			str=br.readLine();  // Skipping header
			str=br.readLine();
			br.close();
			tokens=str.split(" ");
			length=tokens.length;
			lengthProb = new double[length];
			for (i=0; i<length; i++) lengthProb[i]=Double.parseDouble(tokens[i]);
			tokens=null;
			lengthProbCumulative = new double[length];
			lengthProbCumulative[0]=lengthProb[0];
			for (i=1; i<length; i++) lengthProbCumulative[i]=lengthProb[i]+lengthProbCumulative[i-1];
			
			// Loading $lengthProbShortPeriod*$.
			br = new BufferedReader(new FileReader(configDir+"/lengthProbShortPeriod.txt"));
			str=br.readLine();  // Skipping header
			str=br.readLine();
			br.close();
			tokens=str.split(" ");
			length=tokens.length;
			lengthProbShortPeriod = new double[length];
			for (i=0; i<length; i++) lengthProbShortPeriod[i]=Double.parseDouble(tokens[i]);
			tokens=null;
			lengthProbShortPeriodCumulative = new double[length];
			lengthProbShortPeriodCumulative[0]=lengthProbShortPeriod[0];
			for (i=1; i<length; i++) lengthProbShortPeriodCumulative[i]=lengthProbShortPeriod[i]+lengthProbShortPeriodCumulative[i-1];
			
			// Loading $lengthProbLongPeriod*$.
			br = new BufferedReader(new FileReader(configDir+"/lengthProbLongPeriod.txt"));
			str=br.readLine();  // Skipping header
			str=br.readLine();
			br.close();
			tokens=str.split(" ");
			length=tokens.length;
			lengthProbLongPeriod = new double[length];
			for (i=0; i<length; i++) lengthProbLongPeriod[i]=Double.parseDouble(tokens[i]);
			tokens=null;
			lengthProbLongPeriodCumulative = new double[length];
			lengthProbLongPeriodCumulative[0]=lengthProbLongPeriod[0];
			for (i=1; i<length; i++) lengthProbLongPeriodCumulative[i]=lengthProbLongPeriod[i]+lengthProbLongPeriodCumulative[i-1];
			
			// Loading $frequencyProb*$.
			br = new BufferedReader(new FileReader(configDir+"/frequencyProb.txt"));
			str=br.readLine();  // Skipping header
			str=br.readLine();
			br.close();
			tokens=str.split(" ");
			length=tokens.length;
			frequencyProb = new double[length];
			for (i=0; i<length; i++) frequencyProb[i]=Double.parseDouble(tokens[i]);
			tokens=null;
			frequencyProbCumulative = new double[length];
			frequencyProbCumulative[0]=frequencyProb[0];
			for (i=1; i<length; i++) frequencyProbCumulative[i]=frequencyProb[i]+frequencyProbCumulative[i-1];
			
			// Loading $frequencyProbPeriodic*$.
			br = new BufferedReader(new FileReader(configDir+"/frequencyProbPeriodic.txt"));
			str=br.readLine();  // Skipping header
			str=br.readLine();
			br.close();
			tokens=str.split(" ");
			length=tokens.length;
			frequencyProbPeriodic = new double[length];
			for (i=0; i<length; i++) frequencyProbPeriodic[i]=Double.parseDouble(tokens[i]);
			tokens=null;
			frequencyProbPeriodicCumulative = new double[length];
			frequencyProbPeriodicCumulative[0]=frequencyProbPeriodic[0];
			for (i=1; i<length; i++) frequencyProbPeriodicCumulative[i]=frequencyProbPeriodic[i]+frequencyProbPeriodicCumulative[i-1];
			
			// Loading Repbase strings, if any.
			repbaseFile=configDir+"/repbase.fa";
			if (new File(repbaseFile).exists()) loadRepbase(repbaseFile,minAlignmentLength);
			else { lastRepbase=-1; lastRepbaseSat=-1; }
		}
		
		
		public final int getType(Random random) {
			double p = random.nextDouble();
			int out = Arrays.binarySearch(typeProbCumulative,p);
			
System.err.println("getType> choosing type "+out+" since p="+p);			
			
			return out<0?-out-1:out;
		}
		
		
		/**
		 * @return the length of a non-periodic repeat (>=minAlignmentLength).
		 */
		public final int getLength(Random random) {
			int i = Arrays.binarySearch(lengthProbCumulative,random.nextDouble());
			if (i<0) i=-i-1;
			return i==0?minAlignmentLength:i*minAlignmentLength+1+random.nextInt(minAlignmentLength);
		}
		
		
		/**
		 * @return the length of a short-period instance (>=minAlignmentLength).
		 */
		public final int getLengthShortPeriod(Random random) {
			int i = Arrays.binarySearch(lengthProbShortPeriodCumulative,random.nextDouble());
			if (i<0) i=-i-1;
			return i==0?minAlignmentLength:i*minAlignmentLength+1+random.nextInt(minAlignmentLength);
		}
		
		
		/**
		 * @return the number of period units in a long-period instance (>=1).
		 */
		public final int getLengthLongPeriod(Random random) {
			int i = Arrays.binarySearch(lengthProbLongPeriodCumulative,random.nextDouble());
			return 1+(i<0?-i-1:i);
		}
		
		
		public final int getShortPeriod(Random random) {
			return minShortPeriod+random.nextInt(shortPeriodRange);
		}
		
		
		/**
		 * Loads in $out$ a random string of $length$ basepairs.
		 */
		public final void getString(int length, StringBuilder out, Random random) {
			int i, j;
			
			out.delete(0,out.length());
			out.ensureCapacity(length);
			for (i=0; i<length; i++) {
				j=Arrays.binarySearch(charProbCumulative,random.nextDouble());
				out.append(DNA_ALPHABET[j<0?-j-1:j]);
			}
		}
		
		
		public final void getShortPeriodString(int length, StringBuilder out, Random random) {
			int i, j;
			
			out.delete(0,out.length());
			out.ensureCapacity(length);
			for (i=0; i<length; i++) {
				j=Arrays.binarySearch(charProbShortPeriodCumulative,random.nextDouble());
				out.append(DNA_ALPHABET[j<0?-j-1:j]);
			}
		}
		
		
		public final void getUniqueString(int length, StringBuilder out, Random random) {
			int i, j;
			
			out.delete(0,out.length());
			out.ensureCapacity(length);
			for (i=0; i<length; i++) {
				j=Arrays.binarySearch(charProbUniqueCumulative,random.nextDouble());
				out.append(DNA_ALPHABET[j<0?-j-1:j]);
			}
		}
		
		
		/**
		 * @return a value which is at most $maxPeriodDifference*oldPeriod$ from 
		 * $oldPeriod$.
		 */
		public final int perturbPeriod(int oldPeriod, Random random) {
			final int delta = (int)Math.ceil(oldPeriod*maxPeriodDifference);
			return delta==0?oldPeriod:oldPeriod-delta+random.nextInt(delta<<1);
		}
		
		
		/**
		 * @return a randomly perturbed version of $insertionProb$.
		 */
		public final double[] perturbInsertionProb(double[] insertionProb, Random random) {
			return Math.perturbInSimplex(insertionProb,insertionProbPerturbation,random);
		}
		
		
		/**
		 * @return the frequency of a repeat (>=2); every frequency inside the same bin of 
		 * $frequencyProb$ is considered equally likely.
		 */
		public final int getFrequency(Random random) {
			int i = Arrays.binarySearch(frequencyProbCumulative,random.nextDouble());
			if (i<0) i=-i-1;
			final int from = 10^(i)+1;
			final int to = 10^(i+1);
			return from+random.nextInt(to-from+1);
		}
		
		
		/**
		 * @return the frequency of a periodic repeat (>=1).
		 */
		public final int getFrequencyPeriodic(Random random) {
			int i = Arrays.binarySearch(frequencyProbPeriodicCumulative,random.nextDouble());
			if (i<0) i=-i-1;
			return 1+i;
		}
		
		
		public final void addReplicationNoise(StringBuilder sb, Random random) {
			int i, p;
			int c, cPrime, nBases;
			final int length = sb.length();
			double prob;
		
			p=0;
			while (p<length) {
				c=Arrays.binarySearch(DNA_ALPHABET,0,DNA_ALPHABET.length,sb.charAt(p));
				prob=random.nextDouble();
				if (prob<=mismatchProb) {
					do { cPrime=random.nextInt(4); }
					while (cPrime==c);
					sb.append(DNA_ALPHABET[cPrime]);
					p++;
				}
				else if (prob<=mPlusI) {
					nBases=1+random.nextInt(maxIndelLength);
					for (i=0; i<nBases; i++) {
						cPrime=random.nextInt(4);
						sb.append(DNA_ALPHABET[cPrime]);
					}
				}
				else if (prob<=mPlusIPlusD) {
					nBases=1+random.nextInt(maxIndelLength);
					p+=nBases;
				}
				else {
					sb.append(DNA_ALPHABET[c]);
					p++;
				}
			}
			// Insertion after the last character
			prob=random.nextDouble();
			if (prob>mismatchProb && prob<=mPlusI) {
				nBases=1+random.nextInt(maxIndelLength);
				for (i=0; i<nBases; i++) {
					cPrime=random.nextInt(4);
					sb.append(DNA_ALPHABET[cPrime]);
				}
			}
			sb.delete(0,length);
		}
		
		
		/**
		 * Initializes all global Repbase variables using $file$. A string is considered
		 * satellite or non-satellite based on its Repbase label. Non-satellite strings 
		 * shorter than $minAlignmentLength$ are not loaded. Satellite strings of any 
		 * length are assumed to be periods (possibly long).
		 *
		 * @param file characters not in {A,C,G,T} are discarded.
		 */
		private final void loadRepbase(String file, int minAlignmentLength) throws IOException {
			final int CAPACITY = 100;  // Arbitrary
			final String SAT_LABEL_1 = "satellite";
			final String SAT_LABEL_2 = "sat";
			boolean currentIsSat;
			String str;
			StringBuilder buffer;
			BufferedReader br;
			
			buffer = new StringBuilder();
			repbase = new String[CAPACITY]; 
			repbaseSat = new String[CAPACITY];
			repbaseIDs = new String[CAPACITY];
			repbaseSatIDs = new String[CAPACITY];
			lastRepbase=-1; lastRepbaseSat=-1;
			br = new BufferedReader(new FileReader(file));
			str=br.readLine(); currentIsSat=false;
			while (str!=null) {
				if (str.length()!=0 && str.charAt(0)=='>') {
					cleanBuffer(buffer);
					if (buffer.length()!=0) {
						if (currentIsSat) repbaseSat[lastRepbaseSat]=buffer.toString();
						else {
							if (buffer.length()>=minAlignmentLength) repbase[lastRepbase]=buffer.toString();
							else lastRepbase--;
						}
						buffer.delete(0,buffer.length());
					}
					str=str.toLowerCase();
					if (str.indexOf(SAT_LABEL_1)>=0 || str.indexOf(SAT_LABEL_2)>=0) {
						lastRepbaseSat++;
						if (lastRepbaseSat==repbaseSat.length) {
							String[] newArray = new String[repbaseSat.length<<1];
							System.arraycopy(repbaseSat,0,newArray,0,repbaseSat.length);
							repbaseSat=newArray;
						}
						if (lastRepbaseSat==repbaseSatIDs.length) {
							String[] newArray = new String[repbaseSatIDs.length<<1];
							System.arraycopy(repbaseSatIDs,0,newArray,0,repbaseSatIDs.length);
							repbaseSatIDs=newArray;
						}
						repbaseSatIDs[lastRepbaseSat]=str.substring(1).replaceAll("\t",",");
						currentIsSat=true;
					}
					else {
						lastRepbase++;
						if (lastRepbase==repbase.length) {
							String[] newArray = new String[repbase.length<<1];
							System.arraycopy(repbase,0,newArray,0,repbase.length);
							repbase=newArray;
						}
						if (lastRepbase==repbaseIDs.length) {
							String[] newArray = new String[repbaseIDs.length<<1];
							System.arraycopy(repbaseIDs,0,newArray,0,repbaseIDs.length);
							repbaseIDs=newArray;
						}
						repbaseIDs[lastRepbase]=str.substring(1).replaceAll("\t",",");
						currentIsSat=false;
					}
				}	
				else buffer.append(str.trim().toLowerCase());
				str=br.readLine();
			}
			br.close();
			cleanBuffer(buffer);
			if (buffer.length()!=0) {
				if (currentIsSat) repbaseSat[lastRepbaseSat]=buffer.toString();
				else {
					if (buffer.length()>=minAlignmentLength) repbase[lastRepbase]=buffer.toString();
					else lastRepbase--;
				}
			}
			repbaseIsUsed = new boolean[lastRepbase+1];
			Math.set(repbaseIsUsed,lastRepbase,false);
			repbaseSatIsUsed = new boolean[lastRepbaseSat+1];
			Math.set(repbaseSatIsUsed,lastRepbaseSat,false);
			System.err.println("Repbase file loaded: "+(lastRepbase+1)+" non-satellites, "+(lastRepbaseSat+1)+" satellites.");
		}
		
		
		/**
		 * Removes non-DNA characters from $sb$.
		 */
		private static final void cleanBuffer(StringBuilder sb) {
			int i = 0;
			while (i<sb.length()) {
				if (Arrays.binarySearch(DNA_ALPHABET,0,DNA_ALPHABET.length,sb.charAt(i))<0) sb.deleteCharAt(i);
				else i++;  
			}
		}
		
		
		/**
		 * @return TRUE iff the model contains some Repbase strings.
		 */
		public final boolean hasRepbaseStrings() { return lastRepbase+lastRepbaseSat>0; }		
		
		
		/**
		 * @return TRUE iff some Repbase strings have not been used in the simulation yet.
		 */
		public final boolean hasUnusedRepbaseStrings() { 
			int i, count;
			
			count=0;
			for (i=0; i<=lastRepbase; i++) {
				if (!repbaseIsUsed[i]) count++;
			}
			for (i=0; i<=lastRepbaseSat; i++) {
				if (!repbaseSatIsUsed[i]) count++;
			}
			return count>0;
		}
		
	}
	
	
	
	
	
	
	
	
	// ------------------------------- DBG PROCEDURES ------------------------------------
	
	/**
	 * Builds a de Bruijn graph of repeat instances, i.e. a bidirected graph whose nodes
	 * are k-mer endpoints. A character is a distinct substring of a distinct repeat in a 
	 * distinct orientation. A k-mer is a sequence of $order$ adjacent characters that 
	 * fully occur inside a read of length $readLength$. A (k+1)-mer is defined in a 
	 * similar way, i.e. it must fully occur inside a read.
	 *	
	 * Remark: one could also include the last and the first cropped annotation in a read,
	 * if we observe at least $minAlignmentLength$ bps of them: we don't do this in order 
	 * to be more conservative in the result.
	 *
	 * Remark: frequencies of k-mers and (k+1)-mers are not recorded.
	 *
	 * @param minAlignmentLength repeat instances shorter than this are transformed into 
	 * unique sequence;
	 * @param uniqueMode TRUE: a unique sequence does not match any other unique sequence, 
	 * not even one of the same length; FALSE: unique sequences of the same length are 
	 * assumed to match;
	 * @param outputFile prints to this file a DOT representation of the graph (discarded
	 * if NULL).
	 */
	private static final void buildDBG(int order, int readLength, int minAlignmentLength, boolean uniqueMode, int distanceThreshold, String outputFile, boolean verbose) throws IOException {
		final int CAPACITY = 1000;  // Arbitrary
		final int EDGE_CAPACITY = 2;
		boolean orientation;
		int i, j, k, c;
		int last, length, from, to, idGenerator, nCharacters, nKmers, nEnds, nEdges;
		int currentKmerEnd, previousKmerEnd;
		RepeatInstance currentInstance;
		KmerNode kmerSet, currentNode;
		BufferedWriter bw;
		int[] stack, lastNeighbor;
		RepeatInstance[] instances;
		int[][] neighbors;
		
		// Building alphabet and recoded sequence
		System.err.print("Building alphabet... ");
		instances = new RepeatInstance[CAPACITY];
		i=-1; currentInstance=root.next;
		while (currentInstance!=null) {
			i++;
			if (i==instances.length) {
				RepeatInstance[] newArray = new RepeatInstance[instances.length<<1];
				System.arraycopy(instances,0,newArray,0,instances.length);
				instances=newArray;
			}
			instances[i] = new RepeatInstance(currentInstance,currentInstance.sequenceLength>=minAlignmentLength);
			instances[i].id=i;
			currentInstance=currentInstance.next;
		}
		last=i;
		RepeatInstance.order=RepeatInstance.ORDER_FEATURES;
		Arrays.sort(instances,0,last+1);
		neighbors = new int[last+1][EDGE_CAPACITY];
		lastNeighbor = new int[last+1];
		Math.set(lastNeighbor,last,-1);
		for (i=0; i<last; i++) {
			for (j=i+1; j<=last; j++) {
				if (!instances[j].repeat.id.equalsIgnoreCase(instances[i].repeat.id)) break;
				if (!instances[j].isApproximatelyIdentical(instances[i],distanceThreshold,uniqueMode)) continue;
				lastNeighbor[i]++;
				if (lastNeighbor[i]==neighbors[i].length) {
					int[] newArray = new int[neighbors[i].length<<1];
					System.arraycopy(neighbors[i],0,newArray,0,neighbors[i].length);
					neighbors[i]=newArray;
				}
				neighbors[i][lastNeighbor[i]]=j;
				lastNeighbor[j]++;
				if (lastNeighbor[j]==neighbors[j].length) {
					int[] newArray = new int[neighbors[j].length<<1];
					System.arraycopy(neighbors[j],0,newArray,0,neighbors[j].length);
					neighbors[j]=newArray;
				}
				neighbors[j][lastNeighbor[j]]=i;
			}
		}
		stack = new int[last+1];
		nCharacters=getCharacters(instances,last,neighbors,lastNeighbor,stack);
		System.err.println("DONE, "+(nCharacters<<1)+" characters.");
		if (verbose) {
			for (int x=0; x<=last; x++) System.err.println("Character "+instances[x].characterID+" => repeat instance "+instances[x]);
		}
		
		// Building all k-mers
		System.err.print("Building "+order+"-mers... ");
		RepeatInstance.order=RepeatInstance.ORDER_ID;
		Arrays.sort(instances,0,last+1);
		idGenerator=-1;
		kmerSet = new KmerNode(-1);
		for (i=0; i<=last-order+1; i++) {
			to=i-1; length=0;
			do { length+=instances[++to].sequenceLength; } 
			while (to<last && length<readLength);
			if (length>readLength) to--;
			for (j=i; j<=to-order+1; j++) {
				// Forward orientation
				currentNode=kmerSet;
				for (k=0; k<order; k++) {
					currentNode=currentNode.getChild(instances[j+k].characterID);
					if (currentNode==null) break;
				}
				if (currentNode==null) {
					// RC orientation
					currentNode=kmerSet;
					for (k=0; k<order; k++) {
						c=instances[j+order-1-k].characterID;
						if (c%2==0) c++;
						else c--;
						currentNode=currentNode.getChild(c);
						if (currentNode==null) break;
					}
				}
				if (currentNode==null) {
					currentNode=kmerSet;
					for (k=0; k<order; k++) currentNode=currentNode.addChild(instances[j+k].characterID);
					currentNode.kmer=++idGenerator;
					if (verbose) {
						System.err.print("KMER "+currentNode.kmer+" = ");
						System.err.print(instances[j].characterID+"");
						for (k=1; k<order; k++) System.err.print("-"+instances[j+k].characterID);
						System.err.println();
					}
				}
			}
		}
		nKmers=idGenerator+1;
		System.err.println("DONE: "+nKmers+" distinct "+order+"-mers (merged with their RC).");
		
		// Building the DBG
		System.err.print("Building DBG... ");
		nEnds=nKmers<<1;
		if (neighbors.length<nEnds) {
			int[][] newArray = new int[nEnds][0];
			for (i=0; i<neighbors.length; i++) newArray[i]=neighbors[i];
			for (i=neighbors.length; i<nEnds; i++) newArray[i] = new int[EDGE_CAPACITY];
			neighbors=newArray;
		}
		if (lastNeighbor.length<nEnds) lastNeighbor = new int[nEnds];
		for (i=0; i<nEnds; i++) lastNeighbor[i]=-1;
		// Edges that correspond to k-mers
		for (i=0; i<nKmers; i++) {
			from=i<<1; to=(i<<1)+1;
			lastNeighbor[from]++;
			if (lastNeighbor[from]==neighbors[from].length) {
				int[] newArray = new int[neighbors[from].length<<1];
				System.arraycopy(neighbors[from],0,newArray,0,neighbors[from].length);
				neighbors[from]=newArray;
			}
			neighbors[from][lastNeighbor[from]]=to;
			lastNeighbor[to]++;
			if (lastNeighbor[to]==neighbors[to].length) {
				int[] newArray = new int[neighbors[to].length<<1];
				System.arraycopy(neighbors[to],0,newArray,0,neighbors[to].length);
				neighbors[to]=newArray;
			}
			neighbors[to][lastNeighbor[to]]=from;
		}
		// Edges that correspond to (k+1)-mers
		for (i=0; i<=last-order; i++) {
			to=i-1; length=0;
			do { length+=instances[++to].sequenceLength; } 
			while (to<last && length<readLength);
			if (length>readLength) to--;
			previousKmerEnd=Math.POSITIVE_INFINITY;
			for (j=i; j<=to-order+1; j++) {
				// Forward orientation
				orientation=true;
				currentNode=kmerSet;
				for (k=0; k<order; k++) {
					currentNode=currentNode.getChild(instances[j+k].characterID);
					if (currentNode==null) break;
				}
				if (currentNode==null) {
					// RC orientation
					orientation=false;
					currentNode=kmerSet;
					for (k=0; k<order; k++) {
						c=instances[j+order-1-k].characterID;
						if (c%2==0) c++;
						else c--;
						currentNode=currentNode.getChild(c);
					}
				}
				currentKmerEnd=orientation?(currentNode.kmer<<1):(currentNode.kmer<<1)+1;
				if (previousKmerEnd!=Math.POSITIVE_INFINITY) {
					lastNeighbor[previousKmerEnd]++;
					if (lastNeighbor[previousKmerEnd]==neighbors[previousKmerEnd].length) {
						int[] newArray = new int[neighbors[previousKmerEnd].length<<1];
						System.arraycopy(neighbors[previousKmerEnd],0,newArray,0,neighbors[previousKmerEnd].length);
						neighbors[previousKmerEnd]=newArray;
					}
					neighbors[previousKmerEnd][lastNeighbor[previousKmerEnd]]=currentKmerEnd;
					lastNeighbor[currentKmerEnd]++;
					if (lastNeighbor[currentKmerEnd]==neighbors[currentKmerEnd].length) {
						int[] newArray = new int[neighbors[currentKmerEnd].length<<1];
						System.arraycopy(neighbors[currentKmerEnd],0,newArray,0,neighbors[currentKmerEnd].length);
						neighbors[currentKmerEnd]=newArray;
					}
					neighbors[currentKmerEnd][lastNeighbor[currentKmerEnd]]=previousKmerEnd;
				}
				previousKmerEnd=orientation?(currentNode.kmer<<1)+1:(currentNode.kmer<<1);
			}
		}
		// Removing duplicated edges
		nEdges=0;
		for (i=0; i<nEnds; i++) {
			if (lastNeighbor[i]<=0) {
				nEdges+=lastNeighbor[i]+1;
				continue;
			}
			Arrays.sort(neighbors[i],0,lastNeighbor[i]+1);
			k=0;
			for (j=1; j<=lastNeighbor[i]; j++) {
				if (neighbors[i][j]!=neighbors[i][k]) neighbors[i][++k]=neighbors[i][j];
			}
			lastNeighbor[i]=k;
			nEdges+=k+1;
		}
		nEdges>>=1;
		System.err.println("DONE: "+nEnds+" "+order+"-mer endpoints, "+nEdges+" edges.");
		
		// Printing the DBG
		if (outputFile!=null) {
			bw = new BufferedWriter(new FileWriter(outputFile));
			bw.write("graph G {\n");
			for (i=0; i<nEnds; i++) {
				for (j=0; j<=lastNeighbor[i]; j++) {
					if (neighbors[i][j]>i) bw.write(i+" -- "+neighbors[i][j]+";\n");
				}
			}
			bw.write("}\n");
			bw.close();
		}
	}
	
	
	/**
	 * Sets to $2x$ (respectively, $2x+1$) field $characterID$ of every instance in 
	 * $instances$, if its node in the pairwise similarity graph belongs to connected 
	 * component $x$ and is in forward (respectively, RC) orientation.
	 *
	 * @param stack temporary space, of size at least equal to $nNodes$;
	 * @return number of connected components, i.e. half the number of characters in the
	 * alphabet.
	 */
	private static final int getCharacters(RepeatInstance[] instances, int lastInstance, int[][] neighbors, int[] lastNeighbor, int[] stack) {
		int i, j;
		int from, to, top, lastComponent;
		
		// Computing connected components
		for (i=0; i<=lastInstance; i++) instances[i].characterID=-1;
		lastComponent=-1;
		for (i=0; i<=lastInstance; i++) {
			if (instances[i].characterID!=-1) continue;
			lastComponent++; 
			instances[i].characterID=lastComponent; top=0; stack[0]=i;
			while (top>=0) {
				from=stack[top--];
				for (j=0; j<=lastNeighbor[from]; j++) {
					to=neighbors[from][j];
					if (instances[to].characterID!=-1) continue;
					instances[to].characterID=lastComponent;
					stack[++top]=to;
				}
			}
		}
		
		// Mapping different orientations to different characters
		for (i=0; i<=lastInstance; i++) {
			if (instances[i].orientation) instances[i].characterID<<=1;
			else instances[i].characterID=(instances[i].characterID<<1)+1;
		}
		
		return lastComponent+1;
	}
	
	
	/**
	 * A simple trie for collecting k-mers.
	 */
	private static class KmerNode {
		public static final int CAPACITY = 2;  // Arbitrary
		public int character;  // Label of the edge from the parent
		public KmerNode parent;
		public KmerNode[] children;
		public int lastChild;
		public int kmer;  // ID of the k-mer
	
		public KmerNode(int c) {
			this.character=c;
			parent=null;
			children = new KmerNode[CAPACITY];
			lastChild=-1;
			kmer=-1;
		}
	
		/**
		 * @return the new child created, or the existing child of the node with
		 * label $c$.
		 */
		public KmerNode addChild(int c) {
			for (int i=0; i<=lastChild; i++) {
				if (children[i].character==c) return children[i];
			}
			lastChild++;
			if (lastChild==children.length) {
				KmerNode[] newChildren = new KmerNode[children.length<<1];
				System.arraycopy(children,0,newChildren,0,children.length);
				children=newChildren;
			}
			children[lastChild] = new KmerNode(c);
			children[lastChild].parent=this;
			return children[lastChild];
		}
	
		/**
		 * @return the child of the node with label $c$, or NULL if it does not exist.
		 */
		public KmerNode getChild(int c) {
			for (int i=0; i<=lastChild; i++) {
				if (children[i].character==c) return children[i];
			}
			return null;
		}
	}
	
}