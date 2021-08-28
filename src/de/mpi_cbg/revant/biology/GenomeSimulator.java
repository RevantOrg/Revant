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


public class GenomeSimulator {
	
	public static final char[] DNA_ALPHABET = new char[] {'a','c','g','t'};
	
	/**
	 * Repeats. All unique substrings of the genome are assumed to belong to the same 
	 * artificial repeat, which forms its own evolutionary tree with just one node.
	 */
	public static final int N_REPEAT_TYPES = Constants.INTERVAL_UNKNOWN+1;
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
	 * Drawing constants
	 */
	private static final int ROWS_PER_LINE = 100;
	private static final int ROWS_PER_LINE_OVER_TWO = ROWS_PER_LINE>>1;
	private static final int TRUNCATION_THRESHOLD = 100;
	private static final int TRUNCATION_TAG_WIDTH_PIXELS = 3;
	private static final int TRUNCATION_MULTIPLE = 1000;
	public static final String DEFAULT_FONT = "Barlow";
	
	
	
	/**
	 *
	 * @param args
	 * 0=number of unique basepairs (an int, so at most 2 billion);
	 */
	public static void main(String[] args) throws IOException {
		final int N_UNIQUE_BPS = Integer.parseInt(args[0]);
		final double REPEAT_BPS_OVER_UNIQUE_BPS = Double.parseDouble(args[1]);
		final long RANDOM_SEED = Long.parseLong(args[2]);  // -1 to discard it
		final String REPEAT_MODEL_DIR = args[3];
		final int N_OUTPUT_COLUMNS = Integer.parseInt(args[4]);
		final String OUTPUT_IMAGE = args[5];
		
		final int STRING_CAPACITY = 100000;  // Arbitrary
		final long N_REPEAT_BPS = (long)(REPEAT_BPS_OVER_UNIQUE_BPS*N_UNIQUE_BPS);
		RepeatModel model = new RepeatModel(REPEAT_MODEL_DIR);
		StringBuilder sb = new StringBuilder(STRING_CAPACITY);
		Random random = RANDOM_SEED==-1?new Random():new Random(RANDOM_SEED);
		
		long nRepeatBps;
		
		initializeGenome(N_UNIQUE_BPS,model,sb,random);
		nRepeatBps=0;
		while (nRepeatBps<N_REPEAT_BPS) {
			nRepeatBps+=epoch(model,sb,random);
			System.err.println("Epoch completed, "+nRepeatBps+" repeat bps, "+N_UNIQUE_BPS+" unique bps.");
		}
		drawGenome(N_OUTPUT_COLUMNS,OUTPUT_IMAGE);
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
		unique.initialize(nUniqueBps,model,sb,random);
		repeats = new Repeat[CAPACITY_ROWS];
		repeats[0]=unique;
		lastRepeat=0; lastTree=0; unique.tree=0;
		
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
	 *
	 * @return the number of new repetitive basepairs created by the procedure.
	 */
	private static final long epoch(RepeatModel model, StringBuilder sb, Random random) {
		final int CAPACITY = 1000;  // Arbitrary
		int i, j;
		int repeatTree, insertionTree;
		long newBps, totalNewBps;
		Repeat repeat;		
		
		// Creating a new repeat
		repeat = new Repeat();
		if (lastRepeat==0 || random.nextDouble()>model.fromInstanceProb) {
			// New repeat from scratch
			repeat.initialize(model,sb,random);
			repeatTree=-1;
		}
		else {
			// New repeat from an instance chosen uniformly at random
			cumulativeByTree[0]=0;
			cumulativeByTree[1]=lastByTree[1]+1;
			for (i=2; i<=lastTree; i++) cumulativeByTree[i]=cumulativeByTree[i-1]+lastByTree[i]+1;
			i=1+random.nextInt(cumulativeByTree[lastTree]);
			repeatTree=Arrays.binarySearch(cumulativeByTree,1,lastTree+1,i);
			if (repeatTree<0) repeatTree=-1-repeatTree;
			repeat.initialize(byTree[repeatTree][i-cumulativeByTree[repeatTree-1]-1],model,sb,random);
		}
		lastRepeat++;
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
				System.arraycopy(lastByTree,0,newLastByTree,0,byTree.length);
				lastByTree=newLastByTree;
			}	
			byTree[lastTree] = new RepeatInstance[CAPACITY];
			lastByTree[lastTree]=-1;
			repeatTree=lastTree;
		}
		repeat.tree=repeatTree;
		
		// Replication wave	
		totalNewBps=0;
		for (i=0; i<repeat.frequency; i++) {
			newBps=epoch_impl(repeat,model,sb,random);
			if (newBps==0 && treesWithLongInstances()==0) break;
			totalNewBps+=newBps;
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
	 * In basepairs
	 */
	private static final long getGenomeLength() {
		long out;
		RepeatInstance currentInstance;
		
		out=0;
		currentInstance=root.next;
		do {
			out+=currentInstance.sequenceLength;
			currentInstance=currentInstance.next;
		} while (currentInstance!=null);
		return out;
	}
	
	
	/**
	 * @param nColumns pixel width of the image (one pixel corresponds to one basepair).
	 */
	private static final void drawGenome(int nColumns, String outputFile) throws IOException {
		int x, y;
		int fromX, fromY;
		final long genomeLength = getGenomeLength();
		final int nLines = 1+(int)(genomeLength/nColumns);
		final int nRows = nLines*ROWS_PER_LINE;
		RepeatInstance currentInstance;
		BufferedImage image;
		int[] tmpArray = new int[2];

		image = new BufferedImage(nColumns,nRows,BufferedImage.TYPE_INT_RGB);
		for (x=0; x<nColumns; x++) {
			for (y=0; y<nRows; y++) Colors.setRGB(image,x,y,Colors.COLOR_BACKGROUND);
		}
		currentInstance=root;
		fromX=0; fromY=0;
		while (currentInstance.next!=null) {
			currentInstance=currentInstance.next;
			currentInstance.draw(fromX,fromY,image,nColumns,tmpArray);
			fromX=tmpArray[0]; fromY=tmpArray[1];
		}
		ImageIO.write(image,"png",new File(outputFile));
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
	
	
	public static class RepeatInstance {
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
				missingLength=sequenceLength; windowStart=repeatStart; isOverflow=false;
				while (true) {
					overflow=Math.max(fromX+missingLength-nColumns,0);
					windowEnd=repeatEnd-overflow;
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
		public int id;
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
			id=++lastRepeat;
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
			final int DELTA = model.minAlignmentLength>>2;
			int i;
			
			id=++lastRepeat;
			parent=null;
			type=model.getType(random);
			if (type>=Constants.INTERVAL_DENSE_PREFIX && type<=Constants.INTERVAL_DENSE_SUBSTRING) {
				do { sequenceLength=model.getLength(random); }
				while (sequenceLength<model.minAlignmentLength<<1);
				model.getString(sequenceLength,sb,random);
				sequence=sb.toString();
			}
			else if (type<=Constants.INTERVAL_DENSE_SINGLEDELETION) {
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
			frequency=model.getFrequency(random);
		}
		
		
		/**
		 * Creates a new repeat from an instance of a past one. The new repeat inherits 
		 * the type of the past repeat, a perturbed version of its insertion preferences,
		 * and the exact sequence of the instance, but it gets a new frequency.
		 *
		 * Remark: the procedure does not set field $tree$.
		 *
		 * @param instance assumed to be of length at least $minAlignmentLength$.
		 */
		public final void initialize(RepeatInstance instance, RepeatModel model, StringBuilder sb, Random random) {
			int i, period;
			
			id=++lastRepeat;
			parent=instance.repeat;
			type=parent.type;
			if (type<Constants.INTERVAL_PERIODIC) {
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
			insertionProbCumulative = new double[insertionProb.length];
			insertionProbCumulative[0]=insertionProb[0];
			for (i=1; i<insertionProb.length; i++) insertionProbCumulative[i]=insertionProb[i]+insertionProbCumulative[i-1];
			frequency=model.getFrequency(random);
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
System.err.println("PERKELE> 1  length="+length+"  sequence.length()="+sequence.length());
				p=model.minAlignmentLength+random.nextInt((int)Math.max(sequence.length()-(model.minAlignmentLength<<1)-length,1));
System.err.println("PERKELE> 2  p="+p);
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
			model.addReplicationNoise(sb,random);
			return new RepeatInstance(sb.toString(),this,orientation,repeatStart,repeatEnd);
		}
		
		
		/**
		 * FASTA format
		 */
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
		private double fromInstanceProb;
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
	}

}