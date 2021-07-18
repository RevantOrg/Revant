package de.mpi_cbg.revant.biology;

import java.util.Arrays;
import java.util.Vector;
import java.util.Random;
import java.io.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.Constants;


public class GenomeSimulator {
	public static final char[] DNA_ALPHABET = new char[] {'a','c','g','t'};
	
	/**
	 * Repeats
	 */
	public static final int N_REPEAT_TYPES = Constants.INTERVAL_UNKNOWN+1;
	private static Repeat[] repeats;
	private static int lastRepeat, lastComponent;
	private static Repeat unique;
	
	/**
	 * Root of the linked list that represents the genome. $root.next$ points to the
	 * leftmost real instance.
	 */
	private static RepeatInstance root;
	
	/**
	 * Row $i$ stores all the repeat instances that belong to the $i$-th component, in no
	 * particular order.
	 */
	private static RepeatInstance[][] byComponent;
	private static int[] lastByComponent;
	
	/**
	 * Temporary array: cell $i$ contains the numer of instances in all components $<=i$
	 * (excluding component zero).
	 */
	private static int[] cumulativeByComponent;
	
	
	
	
	
	/**
	 *
	 * @param args
	 * 0=number of unique basepairs (an integer, so at most 2 billion);
	 */
	public static void main(String[] args) throws IOException {
		final int N_UNIQUE_BPS = Integer.parseInt(args[0]);
		final double REPEAT_BPS_OVER_UNIQUE_BPS = Double.parseDouble(args[1]);
		final long RANDOM_SEED = Long.parseLong(args[2]);  // -1 to discard it
		final String REPEAT_MODEL_DIR = args[3];
		
		final int STRING_CAPACITY = 100000;  // Arbitrary
		final long N_REPEAT_BPS = (long)(REPEAT_BPS_OVER_UNIQUE_BPS*N_UNIQUE_BPS);
		RepeatModel model = new RepeatModel(REPEAT_MODEL_DIR);
		StringBuilder sb = new StringBuilder(STRING_CAPACITY);
		Random random = RANDOM_SEED==-1?new Random():new Random(RANDOM_SEED);
		
		long nRepeatBps;
		
		initializeGenome(N_UNIQUE_BPS,model,sb,random);
		nRepeatBps=0;
		while (nRepeatBps<N_REPEAT_BPS) {
			
			
			
		}
	
	
	}
	
	
	
	
	
	// ------------------------------------ GENOME ---------------------------------------
	
	/**
	 * Sets the genome to a non-repetitive random sequence of length $nUniqueBps$.
	 */
	private static final void initializeGenome(int nUniqueBps, RepeatModel model, StringBuilder sb, Random random) {
		final int CAPACITY_ROWS = 100;  // Arbitrary
		final int CAPACITY_COLUMNS = 1000;  // Arbitrary
		RepeatInstance instance;
		
		// Repeats
		lastRepeat=-1; lastComponent=-1;
		unique = new Repeat();
		unique.initialize(nUniqueBps,model,sb,random);
		repeats = new Repeat[CAPACITY_ROWS];
		repeats[0]=unique;
		
		// Instances
		root = new RepeatInstance();
		instance = new RepeatInstance(unique.sequence,unique,true,0,unique.sequenceLength-1);
		root.previous=null; root.next=instance; 
		instance.previous=root; instance.next=null;
		byComponent = new RepeatInstance[CAPACITY_ROWS][0];
		byComponent[0] = new RepeatInstance[CAPACITY_COLUMNS];
		byComponent[0][1]=instance;
		lastByComponent = new int[CAPACITY_ROWS];
		lastByComponent[0]=0;
		cumulativeByComponent = new int[CAPACITY_ROWS];
	}
	
	
	/**
	 *
	 * @return the number of new repetitive basepairs created.
	 */
	private static final long epoch(RepeatModel model, StringBuilder sb, Random random) {
		final int CAPACITY = 1000;  // Arbitrary
		boolean newComponent;
		int i, j;
		int insertionComponent;
		long newBps;
		Repeat repeat;
		
		repeat = new Repeat();
		if (lastRepeat==0 || random.nextDouble()>model.fromInstanceProb) {
			// New repeat from scratch
			repeat.initialize(model,sb,random);
			newComponent=true;
		}
		else {
			// New repeat from instance chosen uniformly at random
			cumulativeByComponent[0]=0;
			cumulativeByComponent[1]=lastByComponent[1]+1;
			for (i=2; i<=lastComponent; i++) cumulativeByComponent[i]=cumulativeByComponent[i-1]+lastByComponent[i]+1;
			i=1+random.nextInt(cumulativeByComponent[lastComponent]);
			j=Arrays.binarySearch(cumulativeByComponent,1,lastComponent+1,i);
			if (j<0) j=-j-1;
			repeat.initialize(byComponent[j][i-cumulativeByComponent[j-1]-1],model,sb,random);
			newComponent=false;
		}
		lastRepeat++;
		if (lastRepeat==repeats.length) {
			Repeat[] newRepeats = new Repeat[repeats.length<<1];
			System.arraycopy(repeats,0,newRepeats,0,repeats.length);
			repeats=newRepeats;
		}
		repeats[lastRepeat]=repeat;
		
		// Replication wave
		if (newComponent && lastComponent==byComponent.length) {
			RepeatInstance[][] newByComponent = new RepeatInstance[byComponent.length<<1][0];
			for (i=0; i<byComponent.length; i++) newByComponent[i] = byComponent[i];
			byComponent=newByComponent;
			byComponent[lastComponent] = new RepeatInstance[CAPACITY];
			int[] newLastByComponent = new int[byComponent.length];
			System.arraycopy(lastByComponent,0,newLastByComponent,0,byComponent.length);
			lastByComponent=newLastByComponent;
			lastByComponent[lastComponent]=-1;
		}
		newBps=0;
		for (i=0; i<repeat.frequency; i++) {
//			----------------->
			
			
			
			
		}
		
		
		return 0;
	}
	
	
	
	
	
	
	
	
	
	
	// --------------------------- REPEATS AND REPEAT INSTANCES --------------------------
	
	public static class RepeatInstance {
		public String sequence;
		public int sequenceLength;
		public Repeat repeat;
		public boolean orientation;
		public int repeatStart, repeatEnd;  // Not useful if the repeat type is single deletion or periodic.
		public RepeatInstance previous, next;  // In genome order
		
		
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
	}
	
	
	private static class Repeat {
		public int component;  
		public int id;
		public Repeat parent;
		public int type;  // -1 for unique sequence
		public String sequence;  // Of the period if periodic
		public int sequenceLength;
		public double[] insertionProb;  // In every component of past repeats
		public double[] insertionProbCumulative;
		public int frequency;
		
		
		/**
		 * Creates a new unique region of given $length$ according to $model$.
		 */
		public final void initialize(int length, RepeatModel model, StringBuilder sb, Random random) {
			id=++lastRepeat;
			parent=null;
			component=++lastComponent;
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
		 */											  
		public final void initialize(RepeatModel model, StringBuilder sb, Random random) {
			int i;
			
			id=++lastRepeat;
			parent=null;
			component=++lastComponent;
			type=model.getType(random);
			if (type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
				do { sequenceLength=model.getLength(random); }
				while (sequenceLength<(model.minAlignmentLength)<<1);
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
			insertionProb=Math.sampleFromSimplex(lastComponent,random);
			insertionProbCumulative[0]=insertionProb[0];
			for (i=1; i<=lastComponent; i++) insertionProbCumulative[i]=insertionProb[i]+insertionProbCumulative[i-1];
			frequency=model.getFrequency(random);
		}
		
		
		/**
		 * Creates a new repeat from an instance of a past one. The new repeat inherits 
		 * the type of the past repeat, a perturbed version of its insertion preferences,
		 * and the exact sequence of the instance, but it gets a new frequency.
		 *
		 * @param instance assumed to be of length at least $minAlignmentLength$.
		 */
		public final void initialize(RepeatInstance instance, RepeatModel model, StringBuilder sb, Random random) {
			int i, period;
			
			id=++lastRepeat;
			parent=instance.repeat;
			component=parent.component;
			type=parent.type;
			if (type<6) {
				sequence=""+instance.sequence;
				sequenceLength=instance.sequenceLength;
			}
			else {
				do { period=model.perturbPeriod(parent.sequenceLength,random); }
				while (period>instance.sequenceLength);
				i=random.nextInt(instance.sequenceLength-period);
				sequence=instance.sequence.substring(i,i+period);
				sequenceLength=period;
			}
			insertionProb=model.perturbInsertionProb(parent.insertionProb,random);
			insertionProbCumulative[0]=insertionProb[0];
			for (i=1; i<=lastComponent; i++) insertionProbCumulative[i]=insertionProb[i]+insertionProbCumulative[i-1];
			frequency=model.getFrequency(random);
		}
		
		
		/**
		 * Remark: (1) prefix-suffix type has equal prob. of generating a prefix or a 
		 * suffix, and uniform prob. of generating every length >=minAlignmentLength;
		 * (2) substring type has uniform prob. of generating every substring length >=
		 * minAlignmentLength, and uniform prob. of first position; (3) single-deletion
		 * ensures that a prefix and a suffix of >=minAlignmentLength are always present;
		 * (4) short-period can have a perturbed phase; (5) long-period can have an 
		 * integer multiple of the period, always in phase.
		 */
		public final RepeatInstance getInstance(RepeatModel model, StringBuilder sb, Random random) {
			boolean orientation;
			int i, p;
			int length, phase, repeatStart, repeatEnd;
			String prefix, suffix;
			
			repeatStart=-1; repeatEnd=-1;
			sb.delete(0,sb.length());
			if (type==Constants.INTERVAL_ALIGNMENT) {
				repeatStart=0; repeatEnd=sequenceLength;
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
				length=random.nextInt(sequence.length()-(model.minAlignmentLength)<<1);
				p=model.minAlignmentLength+random.nextInt(sequence.length()-((model.minAlignmentLength)<<1)-length);
				sb.append(sequence.substring(0,p));
				sb.append(sequence.substring(p+length));
				repeatStart=0; repeatEnd=sequenceLength-1;
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
			addReplicationNoise(sb,model,random);
			return new RepeatInstance(sb.toString(),this,orientation,repeatStart,repeatEnd);
		}
		
		
		public final void insertInComponent(RepeatModel model, StringBuilder sb, Random random) {
			int i;
			int insertionComponent;
			RepeatInstance instance;
			
			instance=getInstance(model,sb,random);
			insertionComponent=Arrays.binarySearch(insertionProbCumulative,random.nextDouble());
			if (insertionComponent<0) insertionComponent=-insertionComponent-1;
			i=random.nextInt(lastByComponent[insertionComponent]+1);
			
			
			
			
			
			
		}
		
		
		public String toString() {
			final int length = insertionProb.length;
			int i;
			
			String out = ">"+id+"/"+type+"/";
			if (length>0) out+="P(0)="+insertionProb[0];
			for (i=1; i<length; i++) out+=",P("+i+")="+insertionProb[i];
			out+="/0_"+(sequenceLength-1)+"\n";
			out+=sequence+"\n";
			return out;
		}
	}
	
	
	private static final void reverseComplement(StringBuilder sb) {
		char c, d;
		int i;
		final int length = sb.length();
		
		for (i=0; i<length/2; i++) {
			c=sb.charAt(i); d=sb.charAt(length-i);
			switch (d) {
				case 'a': sb.setCharAt(i,'t'); break;
				case 'c': sb.setCharAt(i,'g'); break;
				case 'g': sb.setCharAt(i,'c'); break;
				case 't': sb.setCharAt(i,'a'); break;
			}
			switch (c) {
				case 'a': sb.setCharAt(length-i,'t'); break;
				case 'c': sb.setCharAt(length-i,'g'); break;
				case 'g': sb.setCharAt(length-i,'c'); break;
				case 't': sb.setCharAt(length-i,'a'); break;
			}
		}
		if (length%2!=0) {
			i=length/2;
			switch (sb.charAt(i)) {
				case 'a': sb.setCharAt(i,'t'); break;
				case 'c': sb.setCharAt(i,'g'); break;
				case 'g': sb.setCharAt(i,'c'); break;
				case 't': sb.setCharAt(i,'a'); break;
			}
		}
	}
	
	
	private static final void addReplicationNoise(StringBuilder sb, RepeatModel model, Random random) {
		int i, p;
		int c, cPrime, nBases;
		final int length = sb.length();
		double prob;
		
		p=-1;
		while (p<length) {
			prob=random.nextDouble();
			if (p==-1) {
				if (prob>model.mismatchProb && prob<=model.mPlusI) {
					nBases=1+random.nextInt(model.maxIndelLength);
					for (i=0; i<nBases; i++) {
						cPrime=random.nextInt(4);
						sb.append(DNA_ALPHABET[cPrime]);
					}
				}
				else p++;
			}
			else {
				c=sb.charAt(p);
				if (prob<=model.mismatchProb) {
					do { cPrime=random.nextInt(4); }
					while (cPrime==c);
					sb.append(DNA_ALPHABET[cPrime]);
					p++;
				}
				else if (prob<=model.mPlusI) {
					nBases=1+random.nextInt(model.maxIndelLength);
					for (i=0; i<nBases; i++) {
						cPrime=random.nextInt(4);
						sb.append(DNA_ALPHABET[cPrime]);
					}
				}
				else if (prob<=model.mPlusIPlusD) {
					nBases=1+random.nextInt(model.maxIndelLength);
					p+=nBases;
				}
				else {
					sb.append(c);
					p++;
				}
			}
		}
		// Insertion after the last character
		prob=random.nextDouble();
		if (prob>model.mismatchProb && prob<=model.mPlusI) {
			nBases=1+random.nextInt(model.maxIndelLength);
			for (i=0; i<nBases; i++) {
				cPrime=random.nextInt(4);
				sb.append(DNA_ALPHABET[cPrime]);
			}
		}
		sb.delete(0,length);
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
		 * Character generation probabilities
		 */
		private double[] charProb, charProbCumulative;
		private double[] charProbShortPeriod, charProbShortPeriodCumulative;
		private double[] charProbUnique, charProbUniqueCumulative;
		
		/**
		 * Length generation probabilities
		 */
		private double[] lengthProb;  // Prob. of block i*mAL+1..(i+1)*mAL
		private double[] lengthProbCumulative;
		private double[] lengthProbShortPeriod;  // Prob. of block i*mAL+1..(i+1)*mAL
		private double[] lengthProbShortPeriodCumulative;
		private double[] lengthProbLongPeriod;  // Prob. of a tandem of 1..X instances
		private double[] lengthProbLongPeriodCumulative;
		
		/**
		 * From-instance probabilities
		 */
		private double fromInstanceProb;
		private double maxPeriodDifference, maxPhaseDifference;
		private double insertionProbPerturbation;
		
		/**
		 * Sequence probabilities
		 */
		private double rcProbability;
		private double mismatchProb, insertionProb, deletionProb;
		private double mPlusI, mPlusIPlusD;
		private int maxIndelLength;
		
		/**
		 * Frequency probabilities
		 */
		private double[] frequencyProb;  // Prob. of block 10^i+1..10^(i+1)
		private double[] frequencyProbCumulative;

		
		public RepeatModel(String configDir) throws IOException {
			int i, p;
			int length;
			String str;
			BufferedReader br;
			String[] tokens;
			
			br = new BufferedReader(new FileReader(configDir+"/repeatModel.txt"));
			str=br.readLine();
			p=str.indexOf("/");
			minAlignmentLength=Integer.parseInt(p>=0?str.substring(0,p).trim():str.trim());
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
			for (i=0; i<DNA_ALPHABET.length; i++) {
				str=br.readLine();
				p=str.indexOf("/");
				charProb[i]=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			}
			charProbCumulative = new double[DNA_ALPHABET.length];
			charProbCumulative[0]=charProb[0];
			for (i=1; i<DNA_ALPHABET.length; i++) charProbCumulative[i]=charProb[i]+charProbCumulative[i-1];
			
			// Loading $charProbShortPeriod*$.
			for (i=0; i<DNA_ALPHABET.length; i++) {
				str=br.readLine();
				p=str.indexOf("/");
				charProbShortPeriod[i]=Double.parseDouble(p>=0?str.substring(0,p).trim():str.trim());
			}
			charProbShortPeriodCumulative = new double[DNA_ALPHABET.length];
			charProbShortPeriodCumulative[0]=charProbShortPeriod[0];
			for (i=1; i<DNA_ALPHABET.length; i++) charProbShortPeriodCumulative[i]=charProbShortPeriod[i]+charProbShortPeriodCumulative[i-1];
			
			// Loading $charProbUnique*$.
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
			tokens[i]=null;
			lengthProbCumulative = new double[length];
			lengthProbCumulative[0]=lengthProb[0];
			for (i=0; i<length; i++) lengthProbCumulative[i]=lengthProb[i]+lengthProbCumulative[i-1];
			
			// Loading $lengthProbShortPeriod*$.
			br = new BufferedReader(new FileReader(configDir+"/lengthProbShortPeriod.txt"));
			str=br.readLine();  // Skipping header
			str=br.readLine();
			br.close();
			tokens=str.split(" ");
			length=tokens.length;
			lengthProbShortPeriod = new double[length];
			for (i=0; i<length; i++) lengthProbShortPeriod[i]=Double.parseDouble(tokens[i]);
			tokens[i]=null;
			lengthProbShortPeriodCumulative = new double[length];
			lengthProbShortPeriodCumulative[0]=lengthProbShortPeriod[0];
			for (i=0; i<length; i++) lengthProbShortPeriodCumulative[i]=lengthProbShortPeriod[i]+lengthProbShortPeriodCumulative[i-1];
			
			// Loading $lengthProbLongPeriod*$.
			br = new BufferedReader(new FileReader(configDir+"/lengthProbLongPeriod.txt"));
			str=br.readLine();  // Skipping header
			str=br.readLine();
			br.close();
			tokens=str.split(" ");
			length=tokens.length;
			lengthProbLongPeriod = new double[length];
			for (i=0; i<length; i++) lengthProbLongPeriod[i]=Double.parseDouble(tokens[i]);
			tokens[i]=null;
			lengthProbLongPeriodCumulative = new double[length];
			lengthProbLongPeriodCumulative[0]=lengthProbLongPeriod[0];
			for (i=0; i<length; i++) lengthProbLongPeriodCumulative[i]=lengthProbLongPeriod[i]+lengthProbLongPeriodCumulative[i-1];
			
			// Loading $frequencyProb*$.
			br = new BufferedReader(new FileReader(configDir+"/frequencyProb.txt"));
			str=br.readLine();  // Skipping header
			str=br.readLine();
			br.close();
			tokens=str.split(" ");
			length=tokens.length;
			frequencyProb = new double[length];
			for (i=0; i<length; i++) frequencyProb[i]=Double.parseDouble(tokens[i]);
			tokens[i]=null;
			frequencyProbCumulative = new double[length];
			frequencyProbCumulative[0]=frequencyProb[0];
			for (i=0; i<length; i++) frequencyProbCumulative[i]=frequencyProb[i]+frequencyProbCumulative[i-1];
		}
		
		
		public final int getType(Random random) {
			int out = Arrays.binarySearch(typeProbCumulative,random.nextDouble());
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
		 * Loads int $out$ a random string of $length$ basepairs.
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
		 * @return a value which is at most $maxPeriodDifference*period$ from $period$.
		 */
		public final int perturbPeriod(int oldPeriod, Random random) {
			final int delta = (int)Math.ceil(oldPeriod*maxPeriodDifference);
			return oldPeriod-delta+random.nextInt(delta<<1);
		}
		
		
		/**
		 * @return a randomly perturbed version of $insertionProb$.
		 */
		public final double[] perturbInsertionProb(double[] insertionProb, Random random) {
			return Math.perturbInSimplex(insertionProb,insertionProbPerturbation,random);
		}
		
		
		/**
		 * @return the frequency of a repeat (>=2); every element inside the same bin of 
		 * $frequencyProb$ is considered equally likely.
		 */
		public final int getFrequency(Random random) {
			int i = Arrays.binarySearch(frequencyProbCumulative,random.nextDouble());
			if (i<0) i=-i-1;
			final int from = 10^(i)+1;
			final int to = 10^(i+1);
			return from+random.nextInt(to-from+1);
		}
	}

}