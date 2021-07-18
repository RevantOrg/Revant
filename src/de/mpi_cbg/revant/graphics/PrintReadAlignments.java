package de.mpi_cbg.revant.graphics;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.imageio.*;
import java.text.*;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.DensityEstimationTree;
import de.mpi_cbg.revant.util.Histograms;
import de.mpi_cbg.revant.util.Colors;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignment;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Factorize;
import de.mpi_cbg.revant.factorize.Reads;
import de.mpi_cbg.revant.factorize.ReadA;
import de.mpi_cbg.revant.factorize.PeriodicSubstringInterval;
import de.mpi_cbg.revant.factorize.DenseSubstring;
import de.mpi_cbg.revant.factorize.AlignmentInterval;
import de.mpi_cbg.revant.intervalgraph.IntervalGraph;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep4;


/**
 *
 */
public class PrintReadAlignments {
	
	private static final int QUANTUM = 10;
	
	public static final int INTERVAL_HEIGHT = 20;
	
	private static final int MAXIMALITY_RADIUS = 2;  // in pixels
	private static final int B_MAXIMALITY_RADIUS = 0;  // in pixels
	private static final int SPACE_BEFORE_INTERVAL = 40;
	public static final int SMALL_SPACE = 5;
	public static final int HISTOGRAM_ROWS = 100;  // in pixels
	public static final int HISTOGRAM_THICKNESS = 4;  // in pixels
	private static final int MAX_INTERVALS_PER_READ = 30000;
	private static final int READ_AHEAD_LIMIT = MAX_INTERVALS_PER_READ*100;  // Characters
	private static final int TRACK_HEIGHT = QUANTUM;  // In pixels
	private static final int TRACK_SPACE = 2;  // In pixels
	public static final String DEFAULT_FONT = "Barlow";
	
	// Subcoverage plots
	private static final int SUBCOVERAGE_HEIGHT = 15;
	private static final int SUBCOVERAGE_SAMPLING = 8;
	private static final int SUBCOVERAGE_THICKNESS = 1;
	private static final int SUBCOVERAGE_COLOR = 0x002C2C2C;

	// Data structures for plotting intervals
	private static int[] tmpArray;
	private static int[][] periodicIntervals, denseIntervals, alignmentIntervals;
	private static int lastIntervalPeriodic, lastIntervalDense, lastIntervalAlignment, lastKernelTag;
	private static DisplayInterval[] intervals;
	private static DisplayAlignment[] alignments;
	private static int[] subcoverage;
	private static int[][] kernelTags;
	private static String[] kernelIDs;
	private static String[] referenceNames, alienReferenceNames;
	private static int[] referenceLengths, alienReferenceLengths;


	public static void main(String[] args) throws Exception {
		int i, j, k, x, y;
		int lastDisplayAlignment, lastInterval, previousReadA;
		int nRows, nColumns, fromY, guideLimit;
		double maxErrorRate;
		String str;
		String referenceAlignmentsFile, referenceNamesFile, referenceLengthsFile;
		String alienReferenceAlignmentsFile, alienReferenceNamesFile, alienReferenceLengthsFile;
		String auxiliaryCoverageFile, kernel2descriptorDir;
		DisplayInterval tmpInterval = new DisplayInterval();
		DisplayInterval previousInterval, currentInterval;
		BufferedImage image;
		BufferedReader denseSubstringIntervals, periodicIntervals, alignmentIntervals, intervalConnectionFile;
		BufferedReader lowComplexityTrackFile, tandemTrackFile, kernelTags, referenceAlignments, alienReferenceAlignments;
		
		// Allocating
		denseSubstringIntervals=args[18].equalsIgnoreCase("null")?null:new BufferedReader(new FileReader(args[18]));
		periodicIntervals=args[19].equalsIgnoreCase("null")?null:new BufferedReader(new FileReader(args[19]));
		alignmentIntervals=args[20].equalsIgnoreCase("null")?null:new BufferedReader(new FileReader(args[20]));
		intervalConnectionFile=args[22].equalsIgnoreCase("null")?null:new BufferedReader(new FileReader(args[22]));
		referenceAlignmentsFile=args[23].equalsIgnoreCase("null")?null:args[23];
		referenceLengthsFile=args[24].equalsIgnoreCase("null")?null:args[24];
		referenceNamesFile=args[25].equalsIgnoreCase("null")?null:args[25];
		alienReferenceAlignmentsFile=args[26].equalsIgnoreCase("null")?null:args[26];
		alienReferenceLengthsFile=args[27].equalsIgnoreCase("null")?null:args[27];
		alienReferenceNamesFile=args[28].equalsIgnoreCase("null")?null:args[28];
		auxiliaryCoverageFile=args[29].equalsIgnoreCase("null")?null:args[29];
		lowComplexityTrackFile=args[30].equalsIgnoreCase("null")?null:new BufferedReader(new FileReader(args[30]));
		tandemTrackFile=args[31].equalsIgnoreCase("null")?null:new BufferedReader(new FileReader(args[31]));
		kernelTags=args[32].equalsIgnoreCase("null")?null:new BufferedReader(new FileReader(args[32]));
		maxErrorRate=Double.parseDouble(args[33]);
		final boolean TRANSPARENT = Integer.parseInt(args[34])==1;
		kernel2descriptorDir=args[35].equalsIgnoreCase("null")?null:args[35];
		final boolean PRINT_MAXIMAL_ENDS = true;
		initIntervals();
		args[18]="null"; args[19]="null"; args[20]="null"; args[22]="null";  // To avoid overwriting interval files from $Factorize.loadInput()$.
		Factorize.loadInput(args);
		ReadA.allocateMemory(Reads.maxReadLength,Alignments.maxAlignmentsPerRead);
		initAlignments();
		DensityEstimationTree.allocateMemory(Alignments.maxAlignmentsPerRead);
		initColors();
		referenceNames=referenceNamesFile==null?null:loadReferenceNames(referenceNamesFile);
		referenceLengths=referenceLengthsFile==null?null:loadReferenceLengths(referenceLengthsFile);
		alienReferenceNames=alienReferenceNamesFile==null?null:loadReferenceNames(alienReferenceNamesFile);
		alienReferenceLengths=alienReferenceLengthsFile==null?null:loadReferenceLengths(alienReferenceLengthsFile);
		if (referenceAlignmentsFile!=null) {
			referenceAlignments = new BufferedReader(new FileReader(referenceAlignmentsFile),IO.BUFFER_SIZE);
			str=referenceAlignments.readLine();  // Skipping the header 
			str=referenceAlignments.readLine();
			referenceAlignments.mark(READ_AHEAD_LIMIT);
		}
		else referenceAlignments=null;
		if (alienReferenceAlignmentsFile!=null) {
			alienReferenceAlignments = new BufferedReader(new FileReader(alienReferenceAlignmentsFile),IO.BUFFER_SIZE);
			str=alienReferenceAlignments.readLine();  // Skipping the header 
			str=alienReferenceAlignments.readLine();
			alienReferenceAlignments.mark(READ_AHEAD_LIMIT);
		}
		else alienReferenceAlignments=null;
		IntervalGraphStep4.loadKernel2Descriptor(kernel2descriptorDir);
		
		// Scanning alignments
		i=0; previousReadA=-1;
		while (i<Alignments.nAlignments) {
			// Loading intervals and alignments
			ReadA.initialize(i,false,false);
			if (lowComplexityTrackFile!=null || tandemTrackFile!=null) {
				while (previousReadA<ReadA.id-1) {
					if (lowComplexityTrackFile!=null) lowComplexityTrackFile.readLine();
					if (tandemTrackFile!=null) tandemTrackFile.readLine();
					previousReadA++;
				}
			}
			previousReadA=ReadA.id;
			subcoverage = new int[Reads.getReadLength(ReadA.id)/QUANTUM+1];
			IO.printOut("Printing read "+ReadA.id+" L="+Reads.getReadLength(ReadA.id)+", A="+(ReadA.lastSortedAlignment+1));
			

			/*
int SELECTED_READ = 1465282;  // 0-based
if (ReadA.id>SELECTED_READ) break;
if (ReadA.id<SELECTED_READ) {
	i+=ReadA.lastSortedAlignment+1;
	for (j=0; j<=ReadA.lastSortedAlignment; j++) intervalConnectionFile.readLine();	
	continue;
}
*/

			ReadA.getCoverageHistogram(false);			
			lastInterval=loadIntervals(ReadA.id,denseSubstringIntervals,periodicIntervals,alignmentIntervals,kernelTags);
			if (lastInterval==-1) {
				i+=ReadA.lastSortedAlignment+1;
				for (j=0; j<=ReadA.lastSortedAlignment; j++) intervalConnectionFile.readLine();	
				continue;
			}
			DisplayInterval.order=DisplayInterval.ORDER_TYPE_START;
			if (lastInterval>0) Arrays.sort(intervals,0,lastInterval+1);
			lastDisplayAlignment=-1;
			for (j=0; j<=ReadA.lastSortedAlignment; j++) {
				lastDisplayAlignment++;
				alignments[lastDisplayAlignment].alignment=ReadA.sortedAlignments[lastDisplayAlignment];
				str=intervalConnectionFile.readLine();
				if (str==null || str.length()==0) {
					alignments[lastDisplayAlignment].interval=null;
					continue;
				}
				Factorize.readIntervalConnectionFile(str);
				tmpInterval.type=Factorize.type;
				tmpInterval.start=Factorize.start;
				tmpInterval.end=Factorize.end;
				k=Arrays.binarySearch(intervals,0,lastInterval+1,tmpInterval);
				if (k<0) {
					System.err.println("ERROR: the following alignment is assigned to an interval that is not present in the interval files?!");
					System.err.println("Alignment: "+ReadA.sortedAlignments[lastDisplayAlignment].toStringBoundaries());
					System.err.println("Interval: "+tmpInterval.type+" ["+tmpInterval.start+".."+tmpInterval.end+"]");
					System.err.println("connection file string: "+str);
					System.exit(1);
				}
				else alignments[lastDisplayAlignment].interval=intervals[k];
			}
			
			// Allocating the image
			image = new BufferedImage(1,1,BufferedImage.TYPE_INT_RGB);
			nRows = HISTOGRAM_ROWS+
				    2*SMALL_SPACE+getNRows_qualities()+
				    2*SMALL_SPACE+getNRows_kernelTags()+
				    getNRows_classicalView(lastDisplayAlignment,lastInterval,image)+
					getNRows_explodedView(lastDisplayAlignment,lastInterval,image)+
					SPACE_BEFORE_INTERVAL;
			nColumns=Reads.getReadLength(ReadA.id)/QUANTUM+1;
			nRows=Math.min(nRows,Integer.MAX_VALUE/nColumns);  // Limitation of $BufferedImage$.
			image = new BufferedImage(nColumns,nRows,BufferedImage.TYPE_INT_RGB);
			for (x=0; x<nColumns; x++) {
				for (y=0; y<nRows; y++) Colors.setRGB(image,x,y,Colors.COLOR_BACKGROUND);
			}
			
			// Drawing
			if (alienReferenceAlignmentsFile!=null && alienReferenceLengthsFile!=null) drawReference(image,0,HISTOGRAM_ROWS,alienReferenceAlignments,ReadA.id,alienReferenceNames,alienReferenceLengths,Colors.COLOR_REFERENCE_ALIEN,Colors.COLOR_REFERENCE_ALIEN_HIGHLIGHT,Alignments.minAlignmentLength,maxErrorRate,TRANSPARENT);
			if (referenceAlignmentsFile!=null && referenceLengthsFile!=null) drawReference(image,0,HISTOGRAM_ROWS,referenceAlignments,ReadA.id,referenceNames,referenceLengths,Colors.COLOR_REFERENCE,Colors.COLOR_REFERENCE_HIGHLIGHT,Alignments.minAlignmentLength,maxErrorRate,TRANSPARENT);
			if (lowComplexityTrackFile!=null) drawTrack(0,lowComplexityTrackFile,image,HISTOGRAM_ROWS+2*SMALL_SPACE-TRACK_HEIGHT-TRACK_SPACE,nRows,nColumns);
			if (tandemTrackFile!=null) drawTrack(1,tandemTrackFile,image,HISTOGRAM_ROWS+2*SMALL_SPACE-TRACK_HEIGHT-TRACK_SPACE,nRows,nColumns);
			if (auxiliaryCoverageFile!=null) fromY=drawAuxiliaryCoverageHistogram(auxiliaryCoverageFile,image,0,HISTOGRAM_ROWS-1);
			fromY=drawCoverageHistogram(image,0,HISTOGRAM_ROWS-1);
			fromY=drawQualities(image,HISTOGRAM_ROWS+2*SMALL_SPACE,nRows,nColumns);
			guideLimit=fromY;
			fromY=drawKernelTags(lastKernelTag,image,nRows,nColumns,fromY+2*SMALL_SPACE,guideLimit);
			fromY=drawClassicalView(lastDisplayAlignment,lastInterval,image,nRows,nColumns,PRINT_MAXIMAL_ENDS,fromY+2*SMALL_SPACE,guideLimit);
			drawExplodedView(lastDisplayAlignment,image,nRows,nColumns,PRINT_MAXIMAL_ENDS,fromY,guideLimit);
			
			// Writing to disk
			ImageIO.write(image,"png",new File((lastInterval>=0?"":"notFactorized-")+"read"+(ReadA.id+1)+".png"));
			i+=ReadA.lastSortedAlignment+1;
		}
		if (denseSubstringIntervals!=null) denseSubstringIntervals.close();
		if (periodicIntervals!=null) periodicIntervals.close();
		if (alignmentIntervals!=null) alignmentIntervals.close();
		if (intervalConnectionFile!=null) intervalConnectionFile.close();
		if (kernelTags!=null) kernelTags.close();
		if (referenceAlignments!=null) referenceAlignments.close();
		if (alienReferenceAlignments!=null) alienReferenceAlignments.close();
	}


	private static final void initIntervals() {
		periodicIntervals = new int[MAX_INTERVALS_PER_READ][11];
		denseIntervals = new int[MAX_INTERVALS_PER_READ][14];
		alignmentIntervals = new int[MAX_INTERVALS_PER_READ][8];
		lastIntervalPeriodic=-1; lastIntervalDense=-1; lastIntervalAlignment=-1;
		tmpArray = new int[30];
		intervals = new DisplayInterval[MAX_INTERVALS_PER_READ];
		for (int i=0; i<MAX_INTERVALS_PER_READ; i++) intervals[i] = new DisplayInterval();
		kernelTags = new int[MAX_INTERVALS_PER_READ][0];
		kernelIDs = new String[MAX_INTERVALS_PER_READ];
	}
	
	
	private static final void initAlignments() {
		alignments = new DisplayAlignment[ReadA.sortedAlignments.length];
		for (int i=0; i<alignments.length; i++) alignments[i] = new DisplayAlignment();
	}


	/**
	 * @return the last loaded interval in $intervals$.
	 */
	private static final int loadIntervals(int read, BufferedReader denseBuffer, BufferedReader periodicBuffer, BufferedReader alignmentsBuffer, BufferedReader kernelTagsBuffer) throws IOException {
		int i, j;
		
		// Loading $*Intervals$ arrays.
		if (periodicBuffer!=null) periodicBuffer.mark(READ_AHEAD_LIMIT);
		if (denseBuffer!=null) denseBuffer.mark(READ_AHEAD_LIMIT);
		if (alignmentsBuffer!=null) alignmentsBuffer.mark(READ_AHEAD_LIMIT);
		if (periodicBuffer!=null) lastIntervalPeriodic=PeriodicSubstringInterval.loadPeriodicIntervals(read,periodicBuffer,READ_AHEAD_LIMIT,periodicIntervals,tmpArray);
		else lastIntervalPeriodic=-1;
		if (denseBuffer!=null) lastIntervalDense=DenseSubstring.loadDenseIntervals(read,denseBuffer,READ_AHEAD_LIMIT,denseIntervals,tmpArray);
		else lastIntervalDense=-1;
		if (alignmentsBuffer!=null) lastIntervalAlignment=AlignmentInterval.loadAlignmentIntervals(read,alignmentsBuffer,READ_AHEAD_LIMIT,alignmentIntervals,tmpArray);
		else lastIntervalAlignment=-1;
		
		// Loading kernel tags
		if (kernelTagsBuffer!=null) {
			kernelTagsBuffer.mark(READ_AHEAD_LIMIT);
			lastKernelTag=IntervalGraphStep3.loadKernelTags(read,kernelTagsBuffer,READ_AHEAD_LIMIT,kernelTags,kernelIDs);
		}
		else lastKernelTag=-1;
	
		// Loading $intervals$.
		j=-1;
		for (i=0; i<=lastIntervalPeriodic; i++) intervals[++j].initFromPeriodic(periodicIntervals[i]);
		for (i=0; i<=lastIntervalDense; i++) intervals[++j].initFromDense(denseIntervals[i]);
		for (i=0; i<=lastIntervalAlignment; i++) intervals[++j].initFromAlignment(alignmentIntervals[i]);
		return j;
	}
	
	
	private static final int drawKernelTags(int lastKernelTag, BufferedImage image, int nRows, int nColumns, int from, int guideLimit) {
		final int MASK = 0x00FFFFFF;
		boolean found;
		int i, j, k;
		int fromY, nextFromY, fromX, toX, color;
		
		fromY=from;
		for (i=0; i<=lastKernelTag; i++) {
			j=from; fromX=kernelTags[i][0]/QUANTUM; toX=kernelTags[i][1]/QUANTUM;
			while (true) {
				found=false;
				for (k=fromX; k<=toX; k++) {
					color=Colors.getRGB(image,k,j);
					if ( (color&MASK)==(Colors.COLOR_KERNELTAG&MASK) || (color&MASK)==(Colors.COLOR_KERNELTAG_HIGHLIGHT&MASK) || 
						 (color&MASK)==(Colors.COLOR_KERNELTAG_FULL&MASK) || (color&MASK)==(Colors.COLOR_KERNELTAG_FULL_HIGHLIGHT&MASK) || 
						 (color&MASK)==(Colors.COLOR_KERNELTAG_NOKERNEL_HIGHLIGHT&MASK) ||
					     Colors.isTypeColor(color)
					   ) {
						found=true;
						break;
					}
				}
				if (found) j+=INTERVAL_HEIGHT+QUANTUM/2;
				else break;
			}
			nextFromY=drawKernelTag(i,image,j,true,guideLimit)+QUANTUM/2;
			fromY=Math.max(fromY,nextFromY);
		}
		return fromY;
	}
	
	
	/**
	 * Remark: if an interval has multiple kernel tags, the procedure prints the smallest 
	 * difference between the length of the interval, and the length of a kernel that does
	 * not coincide with the whole interval.
	 *
	 * @param fromY the first Y coordinate available for drawing;
	 * @param guideLimit smallest Y coordinate up to which guides should be printed;
	 * @return the first Y coordinate available for drawing after the procedure 
	 * completes.
	 */
	private static final int drawKernelTag(int tag, BufferedImage image, int fromY, boolean printGuides, int guideLimit) {
		final int MASK = 0x00FFFFFF;
		final int OFFSET = 5;
		final int TRUNCATION_TAG_WIDTH_PIXELS = 3;
		final int UNASSIGNED_TAG_WIDTH_PIXELS = 2;
		final int LENGTH_THRESHOLD = IO.quantum<<1;  // Arbitrary
		boolean assignedToKernel;
		int i, x, y;
		int fromS, toS, fromE, toE, fromL, length, minLength, kernel, kernelLength, color;
		final int start = kernelTags[tag][0]; 
		final int end = kernelTags[tag][1];
		final boolean inReference = kernelTags[tag][2]==1;
		final int lastKernel = kernelTags[tag][3];
		final int lastPathWithStart = kernelTags[tag][4];
		final int lastPathWithEnd = kernelTags[tag][5];
		final int firstKernel = lastKernel==-1?-1:kernelTags[tag][6];
		final int firstKernelOrientation = lastKernel==-1?-1:kernelTags[tag][6+(lastKernel+1)];
		final int firstKernelLength = lastKernel==-1?-1:kernelTags[tag][6+((lastKernel+1)<<1)];
		final int colorPlain, colorHighlight;
		String label, lengthLabel;
		
		// Deciding text
		final NumberFormat formatter = NumberFormat.getInstance();
		formatter.setMaximumFractionDigits(2);
		if (kernelIDs[tag].length()==0) {
			// Tag with no component, printed by Step1. The first and only kernel, if any,
			// is a component ID.
			label=firstKernel==-1?"":(firstKernel+"");
			assignedToKernel=false;
			colorPlain=Colors.COLOR_KERNELTAG_NOKERNEL;
			colorHighlight=Colors.COLOR_KERNELTAG_NOKERNEL_HIGHLIGHT;
		}
		else if (lastKernel==-1) {
			label=kernelIDs[tag]+"";
			assignedToKernel=false;
			colorPlain=Colors.COLOR_KERNELTAG_NOKERNEL;
			colorHighlight=Colors.COLOR_KERNELTAG_NOKERNEL_HIGHLIGHT;
		}
		else {
			length=end-start+1;
			if (lastKernel==0) {
				if (lastPathWithStart!=-1 && lastPathWithEnd!=-1 && Math.abs(firstKernelLength,length)<=LENGTH_THRESHOLD) lengthLabel="(WHOLE)";
				else if (firstKernelLength>0 && firstKernelLength>=length) lengthLabel="("+formatter.format(((double)(firstKernelLength-length))/1000)+")";
				else lengthLabel="";
				label=getKernelString(kernelIDs[tag],firstKernel)+"   "+lengthLabel;
			}
			else {
				minLength=Math.POSITIVE_INFINITY;
				fromS=6+((lastKernel+1)*3); toS=fromS+lastPathWithStart;
				fromE=toS+1; toE=fromE+lastPathWithEnd;
				fromL=6+((lastKernel+1)<<1);
				for (i=0; i<=lastKernel; i++) {
					kernel=kernelTags[tag][6+i];
					if ( (lastPathWithStart==-1 || (Arrays.binarySearch(kernelTags[tag],fromS,toS+1,kernel)<0 && Arrays.binarySearch(kernelTags[tag],fromS,toS+1,-1-kernel)<0)) ||
					     (lastPathWithEnd==-1 || (Arrays.binarySearch(kernelTags[tag],fromE,toE+1,kernel)<0 && Arrays.binarySearch(kernelTags[tag],fromE,toE+1,-1-kernel)<0))
					   ) {
						kernelLength=kernelTags[tag][fromL+i];
						if (kernelLength>=length && kernelLength<minLength) minLength=kernelLength;
					}
				}
				if (minLength!=Math.POSITIVE_INFINITY) lengthLabel="("+formatter.format(((double)(minLength-length))/1000)+")";
				else lengthLabel="";
				label=getKernelString_multipleKernels(kernelIDs[tag],kernelTags[tag],6,6+lastKernel)+"+   "+lengthLabel;
			}
			assignedToKernel=true;
			colorPlain=inReference?Colors.COLOR_KERNELTAG_FULL:Colors.COLOR_KERNELTAG;
			colorHighlight=inReference?Colors.COLOR_KERNELTAG_FULL_HIGHLIGHT:Colors.COLOR_KERNELTAG_HIGHLIGHT;
		}
		
		// Guides
		if (printGuides) {
			x=start;
			for (y=fromY-1; y>=guideLimit; y--) {
				if ((Colors.getRGB(image,x/QUANTUM,y)&MASK)==(Colors.COLOR_BACKGROUND&MASK) || (Colors.getRGB(image,x/QUANTUM,y)&MASK)==(Colors.COLOR_BACKGROUND_LIGHT&MASK)) Colors.setRGB(image,x/QUANTUM,y,Colors.COLOR_GUIDE);
			}
			x=end;
			for (y=fromY-1; y>=guideLimit; y--) {
				if ((Colors.getRGB(image,x/QUANTUM,y)&MASK)==(Colors.COLOR_BACKGROUND&MASK) || (Colors.getRGB(image,x/QUANTUM,y)&MASK)==(Colors.COLOR_BACKGROUND_LIGHT&MASK)) Colors.setRGB(image,x/QUANTUM,y,Colors.COLOR_GUIDE);
			}
		}
		
		// Interval
		if (assignedToKernel) {
			for (x=start; x<=end; x++) {
				if ( (lastPathWithStart==-1 && ((x-start)/QUANTUM)<=TRUNCATION_TAG_WIDTH_PIXELS) || 
					 (lastPathWithEnd==-1 && ((end-x)/QUANTUM)<=TRUNCATION_TAG_WIDTH_PIXELS)
				   ) color=colorHighlight;
				else color=colorPlain;
				for (y=fromY; y<fromY+INTERVAL_HEIGHT; y++) Colors.setRGB(image,x/QUANTUM,y,color);
			}
		}
		else {
			color=colorPlain;
			for (x=start; x<=end; x++) {
				for (y=fromY; y<fromY+INTERVAL_HEIGHT; y++) Colors.setRGB(image,x/QUANTUM,y,color);
			}
			color=colorHighlight;
			for (y=fromY; y<fromY+INTERVAL_HEIGHT; y++) {
				for (x=0; x<UNASSIGNED_TAG_WIDTH_PIXELS; x++) Colors.setRGB(image,(start/QUANTUM)+x,y,color);
			}
			for (y=fromY; y<fromY+INTERVAL_HEIGHT; y++) {
				for (x=0; x<UNASSIGNED_TAG_WIDTH_PIXELS; x++) Colors.setRGB(image,(end/QUANTUM)-x,y,color);
			}
			for (x=start; x<=end; x++) {
				for (y=0; y<UNASSIGNED_TAG_WIDTH_PIXELS; y++) Colors.setRGB(image,x/QUANTUM,fromY+y,color);
			}
			for (x=start; x<=end; x++) {
				for (y=0; y<UNASSIGNED_TAG_WIDTH_PIXELS; y++) Colors.setRGB(image,x/QUANTUM,fromY+INTERVAL_HEIGHT-1-y,color);
			}
		}
		
		// Drawing text
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
		g2d.setColor(new Color(colorHighlight));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
		x=((start+end)/2)/QUANTUM-Math.round(fontMetrics.stringWidth(label),2);
		y=fromY+INTERVAL_HEIGHT-(INTERVAL_HEIGHT-fontMetrics.getAscent());
		g2d.drawString(label,x,y);
		
		// Orientation
		if (assignedToKernel && lastKernel==0) {
			if (firstKernelOrientation==0 || firstKernelOrientation==2) drawDirectionalTriangle(true,end/QUANTUM-TRUNCATION_TAG_WIDTH_PIXELS-1,fromY+(INTERVAL_HEIGHT>>1),colorHighlight,image);
			if (firstKernelOrientation==1 || firstKernelOrientation==2) drawDirectionalTriangle(false,start/QUANTUM+TRUNCATION_TAG_WIDTH_PIXELS+1,fromY+(INTERVAL_HEIGHT>>1),colorHighlight,image);
		}
	
		return fromY+INTERVAL_HEIGHT;
	}
	
	
	/**
	 * @return the string that represents a kernel tag, taking the kernel2descriptor files
	 * built in Step 4 into account, if any.
	 */
	private static final String getKernelString(String prefix, int kernelTag) {
		final int prefixLength = prefix.length();
		final String prefixPrime = prefixLength==0?"":prefix+".";
		int[] oldTags;
		
		if (prefixLength==0 || prefix.indexOf(".")>=0) return prefixPrime+kernelTag;
		oldTags=IntervalGraphStep4.newTag2oldTags(Integer.parseInt(prefix),kernelTag);
		if (oldTags==null) return prefixPrime+kernelTag;
		else return prefixPrime+oldTags[0]+"."+oldTags[1];
	}
	
	
	private static final String getKernelString_multipleKernels(String prefix, int[] tags, int firstTag, int lastTag) {
		final int prefixLength = prefix.length();
		final String prefixPrime = prefixLength==0?"":prefix+".";
		int oldTag;
		
		if (prefixLength==0 || prefix.indexOf(".")>=0) return prefixPrime;
		oldTag=IntervalGraphStep4.newTag2oldTag(Integer.parseInt(prefix),tags,firstTag,lastTag);
		if (oldTag==-1) return prefixPrime;
		else return prefixPrime+oldTag+".";
	}
	
	
	private static final int getNRows_kernelTags() {
		return (lastKernelTag+1)*(INTERVAL_HEIGHT*2);
	}
	
	
	private static final int drawClassicalView(int lastDisplayAlignment, int lastDisplayInterval, BufferedImage image, int nRows, int nColumns, boolean printMaximalEnds, int from, int guideLimit) {
		final int MASK = 0x00FFFFFF;
		boolean found;
		int i, j, k, fromY, previousFromY, nextFromY, fromX, toX;
		
		fromY=from;
		
		// Intervals
		DisplayInterval.order=DisplayInterval.ORDER_START;
		if (lastDisplayInterval>0) Arrays.sort(intervals,0,lastDisplayInterval+1);
		for (i=0; i<=lastDisplayInterval; i++) {
			j=from; fromX=intervals[i].start/QUANTUM; toX=intervals[i].end/QUANTUM;
			while (true) {
				found=false;
				for (k=fromX; k<=toX; k++) {
					if (Colors.isTypeColor(Colors.getRGB(image,k,j))) {
						found=true;
						break;
					}
				}
				if (found) j+=intervals[i].getDrawingHeight()+QUANTUM/2;
				else break;
			}
			nextFromY=intervals[i].draw(image,j,false,guideLimit)+QUANTUM/2;
			fromY=Math.max(fromY,nextFromY);
		}
		
		// Title
		fromY=drawHeader("INPUT",image,nRows,nColumns,fromY);
		
		// Alignments
		drawDotMatrix(image,nRows,nColumns,from,fromY-1,DOT_GRID_SPACING,true);
		DisplayAlignment.order=DisplayAlignment.ORDER_STARTA_ENDA;
		if (lastDisplayAlignment>0) Arrays.sort(alignments,0,lastDisplayAlignment+1);
		previousFromY=fromY;
		for (i=0; i<=lastDisplayAlignment; i++) fromY=alignments[i].draw(image,fromY,printMaximalEnds,nRows,nColumns,true);
		for (i=0; i<nColumns; i++) {
			for (j=previousFromY; j<fromY; j++) {
				if ((Colors.getRGB(image,i,j)&MASK)==(Colors.COLOR_BACKGROUND&MASK)) Colors.setRGB(image,i,j,Colors.COLOR_BACKGROUND_LIGHT);
			}
		}
		drawDotMatrix(image,nRows,nColumns,previousFromY,fromY-1,DOT_GRID_SPACING,false);
		
		return fromY;
	}
	
	
	private static final int getNRows_classicalView(int lastDisplayAlignment, int lastDisplayInterval, BufferedImage image) {
		return getNRows_drawHeader(image)+(lastDisplayAlignment+1)+(lastDisplayInterval+1)*(INTERVAL_HEIGHT*2);
	}
	
	
	private static final int drawHeader(String label, BufferedImage image, int nRows, int nColumns, int from) {
		final int X_OFFSET = DOT_GRID_SPACING;  // In pixels
		int fromY;
		
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.BOLD,32));
		g2d.setColor(new Color(Colors.COLOR_TEXT_DARK));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
		fromY=from+(SPACE_BEFORE_INTERVAL/2)+fontMetrics.getAscent();
		g2d.drawString(label,X_OFFSET,fromY);
		return fromY+(SPACE_BEFORE_INTERVAL/2);
	}
	
	
	private static final int getNRows_drawHeader(BufferedImage image) {
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,26));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
		return (SPACE_BEFORE_INTERVAL/2)+fontMetrics.getAscent()+(SPACE_BEFORE_INTERVAL/2);
	}
	
	
	/**
	 * This visualization was partially suggested by Gene Myers, but the exploded view 
	 * <a href="https://en.wikipedia.org/wiki/Exploded-view_drawing">dates back to the
	 * Renaissance</a>.
	 */
	private static final int drawExplodedView(int lastDisplayAlignment, BufferedImage image, int nRows, int nColumns, boolean printMaximalEnds, int from, int guideLimit) {
		int i, fromY;
		DisplayInterval currentInterval, previousInterval;
		
		fromY=drawHeader("OUTPUT",image,nRows,nColumns,from);
		DisplayAlignment.order=DisplayAlignment.ORDER_EXPLODED;
		if (lastDisplayAlignment>0) Arrays.sort(alignments,0,lastDisplayAlignment+1);
		previousInterval=null;
		for (i=0; i<=lastDisplayAlignment; i++) {
			currentInterval=alignments[i].interval;
			if (currentInterval==null) {
				if (previousInterval!=null) {
					drawBoundingBox(i-1,fromY-1,image,nRows,nColumns);
					fromY+=(SPACE_BEFORE_INTERVAL*3)>>1;
				}
				previousInterval=null;
				fromY=alignments[i].draw(image,fromY,printMaximalEnds,nRows,nColumns,true);
				continue;
			}
			if (currentInterval!=previousInterval) {
				if (previousInterval!=null) {		
					drawBoundingBox(i-1,fromY-1,image,nRows,nColumns);
					fromY+=(SPACE_BEFORE_INTERVAL*3)>>1;
				}
				fromY=currentInterval.draw(image,fromY,true,guideLimit);
				previousInterval=currentInterval;
				fromY+=SMALL_SPACE;
			}
			fromY=alignments[i].draw(image,fromY,printMaximalEnds,nRows,nColumns,false);
		}
		if (alignments[i-1].interval!=null) drawBoundingBox(i-1,fromY-1,image,nRows,nColumns);
		drawDotMatrix(image,nRows,nColumns,from,fromY-1,DOT_GRID_SPACING,true);
		
		return fromY;
	}
	
	
	/**
	 * @param lastAlignment the last element in $alignments$ that is assigned to a given
	 * interval;
	 * @param lastY the Y coordinate of the last element,
	 */
	private static final void drawBoundingBox(int lastAlignment, int lastY, BufferedImage image, int nRows, int nColumns) {
		final int OFFSET_PIXELS = 10;
		int i, j, k, x, y, p;
		int minX, maxX, minY, color, colorDarker, from, to, first, last, firstAlignment, maxSubcoverage;
		int nLeftMaximal, nRightMaximal, nLeftMaximalB, nRightMaximalB, nLeftRightMaximal, maxWidth;
		double numerator;
		String label;
		Integer readB;
		HashSet<Integer> distinctReadB = new HashSet<Integer>();
		HashSet<Integer> readBWithMaximal = new HashSet<Integer>();
		
		minX=alignments[lastAlignment].alignment.startA(); 
		maxX=alignments[lastAlignment].alignment.endA();
		minY=lastY;
		nLeftMaximal=0; nRightMaximal=0; nLeftRightMaximal=0; nLeftMaximalB=0; nRightMaximalB=0;
		numerator=0.0;
		i=lastAlignment;
		while (i>=0 && alignments[i].interval==alignments[lastAlignment].interval) {
			minX=Math.min(minX,alignments[i].alignment.startA());
			maxX=Math.max(maxX,alignments[i].alignment.endA());
			numerator+=alignments[i].alignment.errorRate();
			minY--;
			if (alignments[i].alignment.isLeftMaximal==1) nLeftMaximal++;
			if (alignments[i].alignment.isRightMaximal==1) nRightMaximal++;
			if (alignments[i].alignment.isLeftMaximal==1 && alignments[i].alignment.isRightMaximal==1) nLeftRightMaximal++;
			if (alignments[i].alignment.isLeftMaximalB==1) nLeftMaximalB++;
			if (alignments[i].alignment.isRightMaximalB==1) nRightMaximalB++;
			readB=Integer.valueOf(Alignments.alignments[alignments[i].alignment.id][1]);
			distinctReadB.add(readB);
			if (alignments[i].alignment.isLeftMaximalB==1 || alignments[i].alignment.isRightMaximalB==1) readBWithMaximal.add(readB);
			i--;
		}
		firstAlignment=i+1;
		
		// Sub-coverage histogram
		for (i=0; i<=(maxX-minX)/QUANTUM; i++) subcoverage[i]=0;
		i=lastAlignment;
		while (i>=0 && alignments[i].interval==alignments[lastAlignment].interval) {
			for (j=(alignments[i].alignment.startA()-minX)/QUANTUM; j<=(alignments[i].alignment.endA()-minX)/QUANTUM; j++) subcoverage[j]++;
			i--;
		}
		maxSubcoverage=0;
		for (i=0; i<=(maxX-minX)/QUANTUM; i++) maxSubcoverage=Math.max(maxSubcoverage,subcoverage[i]);
		
		// Drawing bounding box
		minY-=SMALL_SPACE;
		color=Colors.type2color(alignments[lastAlignment].interval.type);
		first=Math.max(minY,0); 
		last=Math.min(lastY+OFFSET_PIXELS+SUBCOVERAGE_HEIGHT,nRows);
		from=Math.max(minX/QUANTUM-OFFSET_PIXELS,0); 
		to=Math.min(maxX/QUANTUM+OFFSET_PIXELS,nColumns-1);
		for (x=from; x<=to; x++) Colors.setRGB(image,x,first,color);
		for (x=from; x<=to; x++) Colors.setRGB(image,x,last,color);
		first=Math.max(minX/QUANTUM-OFFSET_PIXELS,0);
		last=Math.min(maxX/QUANTUM+OFFSET_PIXELS,nColumns-1);
		from=Math.max(minY,0); 
		to=Math.min(lastY+OFFSET_PIXELS+SUBCOVERAGE_HEIGHT,nRows-1);
		for (y=from; y<=to; y++) Colors.setRGB(image,first,y,color);
		for (y=from; y<=to; y++) Colors.setRGB(image,last,y,color);
		
		// Drawing sub-coverage histogram
		for (i=0; i<=(maxX-minX)/QUANTUM; i+=SUBCOVERAGE_SAMPLING) {
			for (k=i; k<=i+SUBCOVERAGE_THICKNESS; k++) {
				for (j=0; j<=(subcoverage[i]*SUBCOVERAGE_HEIGHT)/maxSubcoverage; j++) Colors.setRGB(image,Math.min(minX/QUANTUM+k,nColumns-1),Math.min(lastY+OFFSET_PIXELS+SUBCOVERAGE_HEIGHT-1-j,nRows-1),color);
				Colors.setRGB(image,Math.min(minX/QUANTUM+k,nColumns-1),Math.min(lastY+OFFSET_PIXELS+SUBCOVERAGE_HEIGHT-1-(subcoverage[i]*SUBCOVERAGE_HEIGHT)/maxSubcoverage,nRows-1),Colors.COLOR_MAXIMAL_DESATURATED);
			}
		}
		
		// Drawing text
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
		g2d.setColor(new Color(Colors.COLOR_TEXT_BOUNDINGBOX));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
		final NumberFormat formatter = NumberFormat.getInstance();
		formatter.setMaximumFractionDigits(2);
		
		// Deciding where to print the alignment stats
		label="A="+(lastAlignment-firstAlignment+1);
		maxWidth=fontMetrics.stringWidth(label);
		maxWidth=Math.max(maxWidth,getNColumns_drawSquares(nLeftMaximal,nLeftMaximalB,image));
		maxWidth=Math.max(maxWidth,getNColumns_drawSquares(nRightMaximal,nRightMaximalB,image));
		maxWidth=Math.max(maxWidth,getNColumns_drawSquares(nLeftRightMaximal,0,image));
		label="B="+distinctReadB.size()+","+readBWithMaximal.size();
		maxWidth=Math.max(maxWidth,fontMetrics.stringWidth(label));
		from=maxX/QUANTUM+OFFSET_PIXELS+OFFSET_PIXELS;
		if (from+maxWidth<nColumns) x=from;
		else {
			from=minX/QUANTUM-OFFSET_PIXELS-OFFSET_PIXELS-maxWidth;
			if (from>=0) x=from;
			else {
				from=minX/QUANTUM;
				if (from+maxWidth<nColumns) x=from;
				else return;  // Stats not printable anywhere
			}
		}
		
		// N. alignments
		label="A="+(lastAlignment-firstAlignment+1);
		from=Math.max(minY+fontMetrics.getAscent(),0); 
		y=from;
		g2d.drawString(label,x,y);
		// Maximality
		if (nLeftMaximal>0 || nLeftMaximalB>0) y=drawSquares(true,false,nLeftMaximal,nLeftMaximalB,image,x,y,nRows,nColumns);
		if (nRightMaximal>0 || nRightMaximalB>0) y=drawSquares(false,true,nRightMaximal,nRightMaximalB,image,x,y,nRows,nColumns);
		if (nLeftRightMaximal>0) y=drawSquares(true,true,nLeftRightMaximal,-1,image,x,y,nRows,nColumns);
		// Distinct B-reads
		y+=fontMetrics.getAscent();
		g2d.drawString("B="+distinctReadB.size()+","+readBWithMaximal.size(),x,y);
		y+=fontMetrics.getAscent();
		g2d.drawString("E="+formatter.format(numerator/(lastAlignment-firstAlignment+1)),x,y);
	}
	
	
	private static final int getNRows_explodedView(int lastDisplayAlignment, int lastDisplayInterval, BufferedImage image) {
		return getNRows_drawHeader(image)+(lastDisplayAlignment+1)+(lastDisplayInterval+1)*(INTERVAL_HEIGHT+SPACE_BEFORE_INTERVAL+SMALL_SPACE+SUBCOVERAGE_HEIGHT+(SPACE_BEFORE_INTERVAL*3)>>1);
	}
	
	
	/**
	 * @param to the last row (included) to be used to draw the coverage histogram.
	 */
	private static final int drawCoverageHistogram(BufferedImage image, int from, int to) {
		final int readLength = Reads.getReadLength(ReadA.id);
		final int max = (int)Histograms.max(ReadA.coverage,0,readLength-1);
		final double resolution = ((double)(to-from))/max;
		int i, j, p;
		int last;
		double average;

		for (i=0; i<readLength; i+=QUANTUM) {
			average=0.0;
			last=Math.min(i+QUANTUM-1,readLength-1);
			for (j=i; j<=last; j++) average+=ReadA.coverage[j];
			average=resolution*average/(last-i+1);
			p=to-(int)average;
			for (j=p; j<=p+(HISTOGRAM_THICKNESS>>1) && j<=to; j++) Colors.setRGB(image,i/QUANTUM,j,Colors.COLOR_HISTOGRAM);
			for (j=p-1; j>=p-(HISTOGRAM_THICKNESS>>1) && j>=0; j--) Colors.setRGB(image,i/QUANTUM,j,Colors.COLOR_HISTOGRAM);
		}
		return to+1;
	}
	
	
	/**
	 * A version of $drawCoverageHistogram()$ that reads the coverage from a file.
	 *
	 * @param path text file with a row of comma-separated coverages for every read.
	 */
	private static final int drawAuxiliaryCoverageHistogram(String path, BufferedImage image, int from, int to) throws IOException {
		final int readLength = Reads.getReadLength(ReadA.id);
		int i, j, p;
		int max, last;
		double average, resolution;
		String str;
		BufferedReader br;
		double[] coverage;
		String[] tokens;
		
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		str=br.readLine();
		while (str!=null) {
			p=str.indexOf(",");
			if (Integer.parseInt(str.substring(0,p))==ReadA.id) {
				tokens=str.split(",");
				coverage = new double[tokens.length-1];
				for (i=1; i<tokens.length; i++) coverage[i-1]=Integer.parseInt(tokens[i]);
				System.err.println("tokens length="+(tokens.length-1)+" readLength="+readLength);		
				max=(int)Histograms.max(coverage,0,readLength-1);
				resolution=((double)(to-from))/max;
				for (i=0; i<readLength; i+=QUANTUM) {
					average=0.0;
					last=Math.min(i+QUANTUM-1,readLength-1);
					for (j=i; j<=last; j++) average+=coverage[j];
					average=resolution*average/(last-i+1);
					p=to-(int)average;
					for (j=p; j<=p+(HISTOGRAM_THICKNESS>>1) && j<=to; j++) Colors.setRGB(image,i/QUANTUM,j,Colors.COLOR_HISTOGRAM_AUXILIARY);
					for (j=p-1; j>=p-(HISTOGRAM_THICKNESS>>1) && j>=0; j--) Colors.setRGB(image,i/QUANTUM,j,Colors.COLOR_HISTOGRAM_AUXILIARY);
				}
				br.close();
				return to+1;
			}
			str=br.readLine();
		}
		br.close();
		return to+1;
	}	

	
	private static class DisplayAlignment implements Comparable {
		public static final int ORDER_STARTA_ENDA = 1;
		public static final int ORDER_EXPLODED = 2;
		public static int order;
		
		public Alignment alignment;
		public DisplayInterval interval;
		
		
		public int compareTo(Object other) {
			int p, length, otherLength, hashCode, otherHashCode;
			DisplayAlignment otherAlignment = (DisplayAlignment)other;
			
			if (order==ORDER_STARTA_ENDA) {
				if (alignment.startA()<otherAlignment.alignment.startA()) return -1;
				else if (alignment.startA()>otherAlignment.alignment.startA()) return 1;
				if (alignment.endA()<otherAlignment.alignment.endA()) return -1;
				else if (alignment.endA()>otherAlignment.alignment.endA()) return 1;
				else return 0;
			}
			else if (order==ORDER_EXPLODED) {
				// Putting all unassigned alignments at the end, sorted by first position.
				if (interval==null) {
					if (otherAlignment.interval==null) return compareStart(otherAlignment.alignment);
					else return 1;
				}
				if (otherAlignment.interval==null) return -1;
			
				// Sorting implied alignments by the first position of their implying 
				// interval.
				if (interval.start<otherAlignment.interval.start) return -1;
				else if (interval.start>otherAlignment.interval.start) return 1;
				if (interval.end<otherAlignment.interval.end) return -1;
				else if (interval.end>otherAlignment.interval.end) return 1;
				
				// Sorting implied alignments by the hashcode of their implying interval
				hashCode=interval.hashCode();
				otherHashCode=otherAlignment.interval.hashCode();
				if (hashCode<otherHashCode) return -1;
				else if (hashCode>otherHashCode) return 1;
			
				// Sorting alignments implied by the same interval
				if ( interval.type==Constants.INTERVAL_ALIGNMENT ||
					 interval.type==Constants.INTERVAL_DENSE_SUBSTRING
				   ) {
					p=compareMaximality(otherAlignment.alignment);
					if (p!=0) return p;
					return compareLengths(otherAlignment.alignment);
				}
				else if ( interval.type==Constants.INTERVAL_DENSE_PREFIX ||
					      interval.type==Constants.INTERVAL_DENSE_SUFFIX ||
						  interval.type==Constants.INTERVAL_DENSE_PREFIXSUFFIX || 
						  interval.type==Constants.INTERVAL_DENSE_SINGLEDELETION ||
						  interval.type==Constants.INTERVAL_PERIODIC  
				        ) {
					p=comparePrefixSuffix(otherAlignment);
					if (p!=0) return p;
					return compareLengths(otherAlignment.alignment);
				}
				else return compareStart(otherAlignment.alignment);
			}
			else return 0;
		}
		
		
		/**
		 * Puts maximal alignments at the beginning.
		 */
		private final int compareMaximality(Alignment otherAlignment) {
			boolean isMaximal = alignment.isLeftMaximal==1 || alignment.isRightMaximal==1;
			boolean isMaximalBoth = alignment.isLeftMaximal==1 && alignment.isRightMaximal==1;
			boolean otherIsMaximal = otherAlignment.isLeftMaximal==1 || otherAlignment.isRightMaximal==1;
			boolean otherIsMaximalBoth = otherAlignment.isLeftMaximal==1 && otherAlignment.isRightMaximal==1;
			
			if (isMaximalBoth) {
				if (!otherIsMaximalBoth) return -1;
				return 0;
			}
			else if (isMaximal) {
				if (otherIsMaximalBoth) return 1;
				else if (otherIsMaximal) return 0;
				else return -1;
			}
			else {
				if (otherIsMaximalBoth || otherIsMaximal) return 1;
				else return 0;
			}
		}
		
		
		/**
		 * Puts longer alignments at the beginning.
		 */
		private final int compareLengths(Alignment otherAlignment) {
			final int length = Alignments.alignments[alignment.id][4]-Alignments.alignments[alignment.id][3]+1;
			final int otherLength = Alignments.alignments[otherAlignment.id][4]-Alignments.alignments[otherAlignment.id][3]+1;
			if (length>otherLength) return -1;
			else if (length<otherLength) return 1;
			return 0;
		}
		
		
		/**
		 * Sorts alignments by first position.
		 */
		private final int compareStart(Alignment otherAlignment) {
			final int start = Alignments.alignments[alignment.id][3];
			final int otherStart = Alignments.alignments[otherAlignment.id][3];
			if (start<otherStart) return -1;
			else if (start>otherStart) return 1;
			return 0;
		}
		
		
		/**
		 * Puts prefix alignments at the beginning, followed by suffix alignments, 
		 * followed by substring alignments.
		 */
		private final int comparePrefixSuffix(DisplayAlignment otherAlignment) {
			final int DISTANCE_THRESHOLD = IO.quantum;
			final int startDistance = Math.abs(Alignments.alignments[alignment.id][3],interval.start);
			final int endDistance = Math.abs(Alignments.alignments[alignment.id][4],interval.end);
			final int otherStartDistance = Math.abs(Alignments.alignments[otherAlignment.alignment.id][3],otherAlignment.interval.start);
			final int otherEndDistance = Math.abs(Alignments.alignments[otherAlignment.alignment.id][4],otherAlignment.interval.end);
			
			if (startDistance<=DISTANCE_THRESHOLD && endDistance<=DISTANCE_THRESHOLD) {
				if (otherStartDistance<=DISTANCE_THRESHOLD && otherEndDistance<=DISTANCE_THRESHOLD) return 0;
				else return -1;
			}
			else if (startDistance<=DISTANCE_THRESHOLD && endDistance>DISTANCE_THRESHOLD) {
				if (otherStartDistance<=DISTANCE_THRESHOLD && otherEndDistance<=DISTANCE_THRESHOLD) return 1;
				if (otherStartDistance<=DISTANCE_THRESHOLD && otherEndDistance>DISTANCE_THRESHOLD) return 0;
				else return -1;
			}
			else if (startDistance>DISTANCE_THRESHOLD && endDistance<=DISTANCE_THRESHOLD) {
				if (otherStartDistance<=DISTANCE_THRESHOLD && otherEndDistance<=DISTANCE_THRESHOLD) return 1;
				else if (otherStartDistance<=DISTANCE_THRESHOLD && otherEndDistance>DISTANCE_THRESHOLD) return 1;
				else if (otherStartDistance>DISTANCE_THRESHOLD && otherEndDistance<=DISTANCE_THRESHOLD) return 0;
				else return -1;
			}
			else {
				if (otherStartDistance<=DISTANCE_THRESHOLD && otherEndDistance<=DISTANCE_THRESHOLD) return 1;
				else if (otherStartDistance<=DISTANCE_THRESHOLD && otherEndDistance>DISTANCE_THRESHOLD) return 1;
				else if (otherStartDistance>DISTANCE_THRESHOLD && otherEndDistance<=DISTANCE_THRESHOLD) return 1;
				else return 0;
			}
		}
		
		
		/**
		 * @param fromY the first Y coordinate available fro drawing;
		 * @param printMaximalEnds TRUE: highlights the maximal start/end of the alignment;
		 * @return the first Y coordinate available for drawing after the procedure completes.
		 */
		public final int draw(BufferedImage image, int fromY, boolean printMaximalEnds, int nRows, int nColumns, boolean saturated) {
			int i, x, y;
			int start, end, xPrime, yPrime, color, radius;

			// Printing alignment
			start=Alignments.alignments[alignment.id][3];
			start-=1;
			if (start<0) start=0;
			end=Alignments.alignments[alignment.id][4];
			end-=1;
			color=diffs2color(Alignments.getAvgDiffs(alignment.id),saturated);
			for (i=start; i<=end; i++) Colors.setRGB(image,i/QUANTUM,fromY,color);
		
			// Printing maximal ends
			color=saturated?Colors.COLOR_MAXIMAL:Colors.COLOR_MAXIMAL_DESATURATED;
			if (printMaximalEnds) {
				if (alignment.isLeftMaximal==1 || alignment.isLeftMaximalB==1) {
					radius=alignment.isLeftMaximal==1?MAXIMALITY_RADIUS:B_MAXIMALITY_RADIUS;
					x=start/QUANTUM;
					y=fromY;
					for (xPrime=x-radius; xPrime<=x+radius; xPrime++) {
						if (xPrime<0 || xPrime>=nColumns) continue;
						for (yPrime=y-radius; yPrime<=y+radius; yPrime++) {
							if (yPrime<0 || yPrime>=nRows) continue;
							Colors.setRGB(image,xPrime,yPrime,color);
						}
					}
				}
				if (alignment.isRightMaximal==1 || alignment.isRightMaximalB==1) {
					radius=alignment.isRightMaximal==1?MAXIMALITY_RADIUS:B_MAXIMALITY_RADIUS;
					x=end/QUANTUM;
					y=fromY;
					for (xPrime=x-radius; xPrime<=x+radius; xPrime++) {
						if (xPrime<0 || xPrime>=nColumns) continue;
						for (yPrime=y-radius; yPrime<=y+radius; yPrime++) {
							if (yPrime<0 || yPrime>=nRows) continue;
							Colors.setRGB(image,xPrime,yPrime,color);
						}
					}
				}
			}
		
			return fromY+1;		
		}
		
	}
	
	
	/**
	 * @param fromY the first Y coordinate available for drawing;
	 * @return the first Y coordinate available for drawing after the procedure completes.
	 */
	private static final int drawQualities(BufferedImage image, int fromY, int nRows, int nColumns) {
		final int READ_ID = ReadA.id;
		final int N_BINS = Reads.getQualityArrayLength(READ_ID);
		final int MAX_Y = Math.min(fromY+QUANTUM,nRows);
		int i, x, y;
		int maxX, color;
		final double[] array = Reads.getQualityArray(READ_ID);
		
		for (i=0; i<N_BINS; i++) {
			color=Reads.quality2color((int)array[i]);
			maxX=Math.min((i+1)*Reads.QUALITY_SPACING/QUANTUM,nColumns);
			for (x=i*Reads.QUALITY_SPACING/QUANTUM; x<maxX; x++) {
				for (y=fromY; y<MAX_Y; y++) Colors.setRGB(image,x,y,color);
			}
		}
		
		return MAX_Y;
	}
	
	
	private static final int getNRows_qualities() {
		return QUANTUM;
	}
	
	
	/**
	 * @param trackID 0=low-complexity; 1=tandem;
	 * @param fromY the first Y coordinate available for drawing;
	 * @return the first Y coordinate available for drawing after the procedure completes.
	 */
	private static final int drawTrack(int trackID, BufferedReader trackFile, BufferedImage image, int fromY, int nRows, int nColumns) throws IOException {
		final int TRACK_OFFSET = Reads.QUALITY_SPACING;  // In string positions
		final int TRACK_RADIUS_LARGE = 4;  // In pixels
		final int TRACK_RADIUS_SMALL = 2;  // In pixels
		int i, x, y;
		int xPrime, start, end, color;
		String str;
		String[] tokens;
		
	    str=trackFile.readLine();
	   	if (str==null || str.length()==0) return fromY;
		tokens=str.split(" ");
		switch (trackID) {
			case 0: color=Colors.COLOR_TRACK_LOWCOMPLEXITY; break;
			case 1: color=Colors.COLOR_TRACK_TANDEM; break;
			default: color=Colors.COLOR_GUIDE;
		}
		for (i=0; i<tokens.length; i+=2) {
			start=Integer.parseInt(tokens[i]);
			end=Integer.parseInt(tokens[i+1]);
			// Start
			x=start/QUANTUM;
			for (xPrime=x-TRACK_RADIUS_LARGE; xPrime<=x+TRACK_RADIUS_LARGE; xPrime++) {
				if (xPrime<0 || xPrime>=nColumns) continue;
				for (y=fromY+(TRACK_HEIGHT>>1)-TRACK_RADIUS_LARGE; y<=fromY+(TRACK_HEIGHT>>1)+TRACK_RADIUS_LARGE; y++) {
					if (y<0 || y>=nRows) continue;
					Colors.setRGB(image,xPrime,y,color);
				}
			}
			// Middle
			for (x=start+TRACK_OFFSET; x<=end; x+=TRACK_OFFSET) {
				for (xPrime=x/QUANTUM-TRACK_RADIUS_SMALL; xPrime<=x/QUANTUM+TRACK_RADIUS_SMALL; xPrime++) {
					if (xPrime<0 || xPrime>=nColumns) continue;
					for (y=fromY+(TRACK_HEIGHT>>1)-TRACK_RADIUS_SMALL; y<=fromY+(TRACK_HEIGHT>>1)+TRACK_RADIUS_SMALL; y++) {
						if (y<0 || y>=nRows) continue;
						Colors.setRGB(image,xPrime,y,color);
					}
				}
			}
			// End
			x=end/QUANTUM;
			for (xPrime=x-TRACK_RADIUS_LARGE; xPrime<=x+TRACK_RADIUS_LARGE; xPrime++) {
				if (xPrime<0 || xPrime>=nColumns) continue;
				for (y=fromY+(TRACK_HEIGHT>>1)-TRACK_RADIUS_LARGE; y<=fromY+(TRACK_HEIGHT>>1)+TRACK_RADIUS_LARGE; y++) {
					if (y<0 || y>=nRows) continue;
					Colors.setRGB(image,xPrime,y,color);
				}
			}
		}
		return fromY+TRACK_HEIGHT;
	}
	
	
	/**
	 * Draws all alignments between reference repeats of length $referenceLengths[]$ and 
	 * read $readID$, as pieces of a horizontal directed triangle that represents the
	 * reference.
	 *
	 * @param br alignments file; the procedure starts reading from the current position
	 * of $br$;
	 * @param minAlignmentLength alignments shorter than this are not drawn;
	 * @param maxErrorRate alignments with error rate bigger than this are not drawn.
	 */
	private static final void drawReference(BufferedImage image, int fromY, int nRows, BufferedReader br, int readID, String[] referenceNames, int[] referenceLengths, int colorPlain, int colorHighlight, int minAlignmentLength, double maxErrorRate, boolean transparent) throws IOException {
		final int HALF_MAX_HEIGHT = nRows>>1;
		final int TRUNCATION_THRESHOLD = 100;
		final int TRUNCATION_TAG_WIDTH_PIXELS = 3;
		final int TRUNCATION_MULTIPLE= 1000;
		final int TRANSPARENT_BORDER_PIXELS = 2;
		int x, y, firstX, lastX, height, color;
		int leftTruncation, rightTruncation, referenceLength;
		double aLength, errorRate;
		String str, label;
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
		g2d.setColor(new Color(colorHighlight));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
		final NumberFormat formatter = NumberFormat.getInstance();
		formatter.setMaximumFractionDigits(2);
		
		str=br.readLine();
		while (str!=null) {
			Alignments.readAlignmentFile(str);
			if (Alignments.readA-1<readID) {
				str=br.readLine();
				continue;
			}
			else if (Alignments.readA-1>readID) {
				br.reset();
				return;
			}
			br.mark(READ_AHEAD_LIMIT);
			
			// Drawing truncated triangle
			errorRate=((double)(Alignments.diffs<<1))/(Alignments.endB-Alignments.startB+Alignments.endA-Alignments.startA+2);
			if (errorRate>maxErrorRate) {
				str=br.readLine();
				continue;
			}
			referenceLength=referenceLengths[Alignments.readB-1];
			if (referenceLength<minAlignmentLength) {
				str=br.readLine();
				continue;
			}
			// Correcting for wrong coordinates by daligner
			if (Alignments.startB<0) Alignments.startB=0;
			if (Alignments.endB>=referenceLength) Alignments.endB=referenceLength;
			firstX=Alignments.startA-(Alignments.orientation?Alignments.startB:referenceLength-Alignments.endB);
			lastX=Alignments.endA+(Alignments.orientation?referenceLength-Alignments.endB:Alignments.startB);
			leftTruncation=Alignments.orientation?Alignments.startB:referenceLength-Alignments.endB;
			rightTruncation=Alignments.orientation?referenceLength-Alignments.endB:Alignments.startB;
			aLength=lastX-firstX+1;
			if (aLength<minAlignmentLength) {
				str=br.readLine();
				continue;
			}
			for (x=Alignments.startA; x<=Alignments.endA; x++) {
				height=Math.round(HALF_MAX_HEIGHT*(Alignments.orientation?lastX-x:x-firstX)/aLength);
				if ( (leftTruncation>TRUNCATION_THRESHOLD && ((x-Alignments.startA)/QUANTUM)<=TRUNCATION_TAG_WIDTH_PIXELS) || 
					 (rightTruncation>TRUNCATION_THRESHOLD && ((Alignments.endA-x)/QUANTUM)<=TRUNCATION_TAG_WIDTH_PIXELS)
				   ) {
					color=colorHighlight;
					for (y=fromY+HALF_MAX_HEIGHT-height; y<=fromY+HALF_MAX_HEIGHT+height; y++) Colors.setRGB(image,x/QUANTUM,y,color);
				}
				else {
					color=colorPlain;
					if (transparent) {
						if (x<=Alignments.startA+TRANSPARENT_BORDER_PIXELS || x>=Alignments.endA-TRANSPARENT_BORDER_PIXELS) {
							for (y=fromY+HALF_MAX_HEIGHT-height; y<=fromY+HALF_MAX_HEIGHT+height; y++) Colors.setRGB(image,x/QUANTUM,y,color);
						} 
						else {
							for (y=fromY+HALF_MAX_HEIGHT-height; y<=fromY+HALF_MAX_HEIGHT-height+TRANSPARENT_BORDER_PIXELS; y++) Colors.setRGB(image,x/QUANTUM,y,color);
							for (y=fromY+HALF_MAX_HEIGHT+height-TRANSPARENT_BORDER_PIXELS; y<=fromY+HALF_MAX_HEIGHT+height; y++) Colors.setRGB(image,x/QUANTUM,y,color);
						}
					}
					else {
						for (y=fromY+HALF_MAX_HEIGHT-height; y<=fromY+HALF_MAX_HEIGHT+height; y++) Colors.setRGB(image,x/QUANTUM,y,color);
					}
				}
			}
			
			// Drawing text
			label=""+referenceNames[Alignments.readB-1].toUpperCase();
			x=((Alignments.startA+Alignments.endA)/2)/QUANTUM-Math.round(fontMetrics.stringWidth(label),2);
			y=fromY+HALF_MAX_HEIGHT+fontMetrics.getAscent()/2;
			g2d.drawString(label,x,y);
			if (leftTruncation>TRUNCATION_THRESHOLD) {
				label=""+formatter.format(((double)leftTruncation)/TRUNCATION_MULTIPLE);
				x=Alignments.startA/QUANTUM+TRUNCATION_TAG_WIDTH_PIXELS+TRUNCATION_TAG_WIDTH_PIXELS;
				y=fromY+HALF_MAX_HEIGHT+fontMetrics.getAscent()/2;
				g2d.drawString(label,x,y);
			}
			if (rightTruncation>TRUNCATION_THRESHOLD) {
				label=""+formatter.format(((double)rightTruncation)/TRUNCATION_MULTIPLE);
				x=Alignments.endA/QUANTUM-TRUNCATION_TAG_WIDTH_PIXELS-TRUNCATION_TAG_WIDTH_PIXELS-fontMetrics.stringWidth(label);
				y=fromY+HALF_MAX_HEIGHT+fontMetrics.getAscent()/2;
				g2d.drawString(label,x,y);
			}
			
			// Next iteration
			str=br.readLine();
		}
		br.close();
	}
	
	
	private static final String[] loadReferenceNames(String file) throws Exception {
		int i;
		String str;
		BufferedReader br;
		String[] names;
		
		// Computing the number of names
		br = new BufferedReader(new FileReader(file));
		i=0; str=br.readLine();
		while (str!=null) {
			i++;
			str=br.readLine();
		}
		br.close();
		
		// Loading names
		names = new String[i];
		br = new BufferedReader(new FileReader(file));
		i=-1; str=br.readLine();
		while (str!=null) {
			names[++i]=str.trim();
			str=br.readLine();
		}
		br.close();
		
		return names;
	}
	
	
	private static final int[] loadReferenceLengths(String file) throws Exception {
		int i;
		String str;
		BufferedReader br;
		int[] lengths;
		
		// Computing the number of lengths
		br = new BufferedReader(new FileReader(file));
		i=0; str=br.readLine();
		while (str!=null) {
			i++;
			str=br.readLine();
		}
		br.close();
		
		// Loading lengths
		lengths = new int[i];
		br = new BufferedReader(new FileReader(file));
		i=-1; str=br.readLine();
		while (str!=null) {
			lengths[++i]=Integer.parseInt(str.trim());
			str=br.readLine();
		}
		br.close();
		
		return lengths;
	}
	
	
	private static class DisplayInterval implements Comparable {
		public static final int ORDER_START = 0;
		public static final int ORDER_TYPE_START = 1;
		public static int order;
		
		public int type;
		public int start, end;
		public int period;
		public boolean hasLongPeriod, equalsOnePeriod;
		public boolean isWeak;
		public int minPrefixLength, minSuffixLength;
		
		
		public final void initFromPeriodic(int[] buffer) {
			type=Constants.INTERVAL_PERIODIC;
			start=buffer[1];
			end=buffer[2];
			period=buffer[6];
			hasLongPeriod=buffer[7]==1;
			equalsOnePeriod=buffer[8]==1;
			isWeak=false;
		}
		
			 
		public final void initFromDense(int[] buffer) {
			type=buffer[1];
			start=buffer[2];
			end=buffer[3];
			isWeak=buffer[6]==1;
			hasLongPeriod=false;
			minPrefixLength=buffer[8];
			minSuffixLength=buffer[9];
		}
		
		
		public final void initFromAlignment(int[] buffer) {
			type=Constants.INTERVAL_ALIGNMENT;
			start=buffer[1];
			end=buffer[2];
			hasLongPeriod=false;
			isWeak=false;
		}
		
		
		public int compareTo(Object other) {
			DisplayInterval otherInterval = (DisplayInterval)other;
			if (order==ORDER_START) {
				if (start<otherInterval.start) return -1;
				else if (start>otherInterval.start) return 1;
				if (end<otherInterval.end) return -1;
				else if (end>otherInterval.end) return 1;
			}
			else if (order==ORDER_TYPE_START) {
				if (type<otherInterval.type) return -1;
				else if (type>otherInterval.type) return 1;
				if (start<otherInterval.start) return -1;
				else if (start>otherInterval.start) return 1;
				if (end<otherInterval.end) return -1;
				else if (end>otherInterval.end) return 1;
			}
			return 0;
		}
		
		
		/**
		 * @param fromY the first Y coordinate available for drawing;
		 * @param guideLimit smallest Y coordinate up to which guides should be printed;
		 * @return the first Y coordinate available for drawing after the procedure 
		 * completes.
		 */
		public final int draw(BufferedImage image, int fromY, boolean printGuides, int guideLimit) {
			final int MASK = 0x00FFFFFF;
			final int OFFSET = 5;
			int x, y;
			int color;
			String label;
			
			// Guides
			if (printGuides) {
				x=start;
				for (y=fromY-1; y>=guideLimit; y--) {
					if ((Colors.getRGB(image,x/QUANTUM,y)&MASK)==(Colors.COLOR_BACKGROUND&MASK) || (Colors.getRGB(image,x/QUANTUM,y)&MASK)==(Colors.COLOR_BACKGROUND_LIGHT&MASK)) Colors.setRGB(image,x/QUANTUM,y,Colors.COLOR_GUIDE);
				}
				x=end;
				for (y=fromY-1; y>=guideLimit; y--) {
					if ((Colors.getRGB(image,x/QUANTUM,y)&MASK)==(Colors.COLOR_BACKGROUND&MASK) || (Colors.getRGB(image,x/QUANTUM,y)&MASK)==(Colors.COLOR_BACKGROUND_LIGHT&MASK)) Colors.setRGB(image,x/QUANTUM,y,Colors.COLOR_GUIDE);
				}
			}
			
			// Interval
			color=Colors.type2color(type);
			for (x=start; x<=end; x++) {
				for (y=fromY; y<fromY+INTERVAL_HEIGHT; y++) Colors.setRGB(image,x/QUANTUM,y,color);
			}
			
			// Text
			final Graphics2D g2d = image.createGraphics();
			g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
			g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
			g2d.setColor(new Color(Colors.COLOR_TEXT));
			final FontMetrics fontMetrics = g2d.getFontMetrics();
			final NumberFormat formatter = NumberFormat.getInstance();
			formatter.setMaximumFractionDigits(2);
			label=""+type2text(type)+(type==Constants.INTERVAL_PERIODIC&&period>0?" "+period:"")+((type==Constants.INTERVAL_DENSE_PREFIX||type==Constants.INTERVAL_DENSE_SUFFIX||type==Constants.INTERVAL_DENSE_SUBSTRING)&&isWeak?" // W":"")+(type==Constants.INTERVAL_PERIODIC&&hasLongPeriod?" // L"+(equalsOnePeriod?"P":""):"");
			x=start/QUANTUM+OFFSET;
			y=fromY+INTERVAL_HEIGHT-(INTERVAL_HEIGHT-fontMetrics.getAscent());
			g2d.drawString(label,x,y);
			label=""+formatter.format(((double)(end-start+1))/1000);
			x=end/QUANTUM-OFFSET-fontMetrics.stringWidth(label);
			y=fromY+INTERVAL_HEIGHT-(INTERVAL_HEIGHT-fontMetrics.getAscent());
			g2d.drawString(label,x,y);
			
			// Directional square
			if (minPrefixLength>0 && type==Constants.INTERVAL_DENSE_PREFIX || type==Constants.INTERVAL_DENSE_PREFIXSUFFIX || type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
				drawDirectionalSquare((start+minPrefixLength)/QUANTUM,fromY+(INTERVAL_HEIGHT>>1),true,type,image,fromY+INTERVAL_HEIGHT,Math.min(image.getWidth(),(end+1)/QUANTUM));
			}
			if (minSuffixLength>0 && type==Constants.INTERVAL_DENSE_SUFFIX || type==Constants.INTERVAL_DENSE_PREFIXSUFFIX || type==Constants.INTERVAL_DENSE_SINGLEDELETION) {
				drawDirectionalSquare((end-minSuffixLength)/QUANTUM,fromY+(INTERVAL_HEIGHT>>1),false,type,image,fromY+INTERVAL_HEIGHT,Math.min(image.getWidth(),(end+1)/QUANTUM));
			}
			if (minPrefixLength>0 && type==Constants.INTERVAL_DENSE_SUBSTRING) {
				x=repositionInterval(start/QUANTUM,(end+1)/QUANTUM,fromY+(INTERVAL_HEIGHT>>1),minPrefixLength/QUANTUM,type,image);
				if (x!=-1) {
					drawDirectionalSquare(x,fromY+(INTERVAL_HEIGHT>>1),false,type,image,fromY+INTERVAL_HEIGHT,Math.min(image.getWidth(),(end+1)/QUANTUM));
					drawDirectionalSquare(x+minPrefixLength/QUANTUM,fromY+(INTERVAL_HEIGHT>>1),true,type,image,fromY+INTERVAL_HEIGHT,Math.min(image.getWidth(),(end+1)/QUANTUM));
				}
			}
			if (period>0 && type==Constants.INTERVAL_PERIODIC) {
				x=repositionInterval(start/QUANTUM,(end+1)/QUANTUM,fromY+(INTERVAL_HEIGHT>>1),period/QUANTUM,type,image);
				if (x!=-1) {
					drawDirectionalSquare(x,fromY+(INTERVAL_HEIGHT>>1),false,type,image,fromY+INTERVAL_HEIGHT,Math.min(image.getWidth(),(end+1)/QUANTUM));
					drawDirectionalSquare(x+period/QUANTUM,fromY+(INTERVAL_HEIGHT>>1),true,type,image,fromY+INTERVAL_HEIGHT,Math.min(image.getWidth(),(end+1)/QUANTUM));
				}
			}
				
			return fromY+INTERVAL_HEIGHT;
		}
		
		
		public final int getDrawingHeight() {
			return INTERVAL_HEIGHT;
		}
		
	}




	// --------------------------------- COLORS ------------------------------------------
	
	private static final int DOT_GRID_SPACING = 20;
	
	
	private static int rFrom, rTo, bFrom, bTo, gFrom, gTo;
	private static int rFromDesaturated, rToDesaturated, bFromDesaturated, bToDesaturated, gFromDesaturated, gToDesaturated;
	
	
	public static final void initColors() {
		int mask;
		
		// Saturated
		mask=0x000000FF;
		bFrom=Colors.GRADIENT_FROM&mask;
		mask<<=8;
		gFrom=(Colors.GRADIENT_FROM&mask)>>8;
		mask<<=8;
		rFrom=(Colors.GRADIENT_FROM&mask)>>16;
		mask=0x000000FF;
		bTo=Colors.GRADIENT_TO&mask;
		mask<<=8;
		gTo=(Colors.GRADIENT_TO&mask)>>8;
		mask<<=8;
		rTo=(Colors.GRADIENT_TO&mask)>>16;
		
		// Desaturated
		mask=0x000000FF;
		bFromDesaturated=Colors.GRADIENT_FROM_DESATURATED&mask;
		mask<<=8;
		gFromDesaturated=(Colors.GRADIENT_FROM_DESATURATED&mask)>>8;
		mask<<=8;
		rFromDesaturated=(Colors.GRADIENT_FROM_DESATURATED&mask)>>16;
		mask=0x000000FF;
		bToDesaturated=Colors.GRADIENT_TO_DESATURATED&mask;
		mask<<=8;
		gToDesaturated=(Colors.GRADIENT_TO_DESATURATED&mask)>>8;
		mask<<=8;
		rToDesaturated=(Colors.GRADIENT_TO_DESATURATED&mask)>>16;
	}
	
	
	public static final int diffs2color(double avgDiffs, boolean saturated) {
		final int MASK = 0x000000FF;
		int rValue, gValue, bValue, out;
		double tmp;
		
		tmp=avgDiffs/Alignments.MAX_ALIGNMENT_ERROR;
		if (tmp>1.0) tmp=1.0;
		if (saturated) {
			rValue=rFrom+Math.round(tmp*(rTo-rFrom+1));
			gValue=gFrom+Math.round(tmp*(gTo-gFrom+1));
			bValue=bFrom+Math.round(tmp*(bTo-bFrom+1));
		}
		else {
			rValue=rFromDesaturated+Math.round(tmp*(rToDesaturated-rFromDesaturated+1));
			gValue=gFromDesaturated+Math.round(tmp*(gToDesaturated-gFromDesaturated+1));
			bValue=bFromDesaturated+Math.round(tmp*(bToDesaturated-bFromDesaturated+1));
		}
		out=0x00000000;
		out|=bValue&MASK;
		out|=(gValue&MASK)<<8;
		out|=(rValue&MASK)<<16;
		if (out==0) out=Colors.GRADIENT_FROM;
		return out;
	}
	
	
	public static final String type2text(int type) {
		switch (type) {
			case Constants.INTERVAL_ALIGNMENT: return "WHOLE";
			case Constants.INTERVAL_DENSE_PREFIX: return "PREFIX";
			case Constants.INTERVAL_DENSE_SUFFIX: return "SUFFIX";
			case Constants.INTERVAL_DENSE_PREFIXSUFFIX: return "PREFIX+SUFFIX";
			case Constants.INTERVAL_DENSE_SUBSTRING: return "SUBSTRING";
			case Constants.INTERVAL_DENSE_SINGLEDELETION: return "DELETION";
			case Constants.INTERVAL_PERIODIC: return "PERIODIC";
			default: return "?";
		}
	}
	
	
	private static final int DOT_SIZE_PIXELS = 2;
	
	
	/**
	 *
	 * @param darkOrLight TRUE=dark, FALSE=light.
	 */
	private static final void drawDotMatrix(BufferedImage image, int nRows, int nColumns, int fromRow, int toRow, int period, boolean darkOrLight) {
		for (int i=fromRow; i<=toRow; i+=period) {
			for (int j=0; j<=nColumns; j+=period) drawDot(j,i,image,nRows,nColumns,darkOrLight);
		}
	}
	
	
	/**
	 *
	 * @param x,y center of the dot;
	 * @param darkOrLight TRUE=dark, FALSE=light.
	 */
	private static final void drawDot(int x, int y, BufferedImage image, int nRows, int nColumns, boolean darkOrLight) {		
		// First column
		Colors.colorIfBackground(x-1,y-1,darkOrLight?0x001A1A1A:0x00282D2F,image,nRows,nColumns);
		Colors.colorIfBackground(x-1,y,darkOrLight?0x00222222:0x002D3235,image,nRows,nColumns);
		Colors.colorIfBackground(x-1,y+1,darkOrLight?0x00181818:0x002F3437,image,nRows,nColumns);
		
		// Second column
		Colors.colorIfBackground(x,y-1,darkOrLight?0x001C1C1C:0x002B3032,image,nRows,nColumns);
		Colors.colorIfBackground(x,y,darkOrLight?0x00262526:0x0032373A,image,nRows,nColumns);
		Colors.colorIfBackground(x,y+1,darkOrLight?0x001C1C1C:0x0035393D,image,nRows,nColumns);
		
		// Third column
		Colors.colorIfBackground(x+1,y-1,darkOrLight?0x00101010:0x002B3032,image,nRows,nColumns);
		Colors.colorIfBackground(x+1,y,darkOrLight?0x00151515:0x0032373A,image,nRows,nColumns);
		Colors.colorIfBackground(x+1,y+1,darkOrLight?0x000F0F0F:0x0035393D,image,nRows,nColumns);
	}


	/**
	 * Draws the full/empty squares and the maximality statistics of a set of alignments.
	 *
	 * @param fullX TRUE: square X is full; FALSE: empty;
	 * @param valueX the number assigned to square X.
	 */
	private static final int drawSquares(boolean fullLeft, boolean fullRight, int valueLeft, int valueRight, BufferedImage image, int fromX, int fromY, int nRows, int nColumns) {
		final int COLOR = Colors.COLOR_TEXT_BOUNDINGBOX;
		final int SIDE = (MAXIMALITY_RADIUS)<<2;
		final int HSPACE = SIDE>>1;
		int x, y, yPrime;
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
		
		if (fromX+SIDE+HSPACE+SIDE+HSPACE>=nColumns) return fromY;
		
		// Text
		if (valueLeft>0 || valueRight>0) {
			String label;
			if (valueLeft!=-1 && valueRight==-1) label=valueLeft+"";
			else if (valueLeft==-1 && valueRight!=-1) label=valueRight+"";
			else if (valueLeft!=-1 && valueRight!=-1) label=valueLeft+","+valueRight;
			else label="";
			g2d.setColor(new Color(Colors.COLOR_TEXT_BOUNDINGBOX));
			y=fromY+fontMetrics.getAscent();
			g2d.drawString(label,fromX+SIDE+HSPACE+SIDE+HSPACE,y);
		}
		
		// Squares
		yPrime=fromY+fontMetrics.getAscent()-SIDE-SIDE/3;
		if (fullLeft) {
			for (x=fromX; x<fromX+SIDE; x++) {
				for (y=yPrime; y<yPrime+SIDE; y++) Colors.setRGB(image,x,y,COLOR);
			}
		}
		else {
			for (x=fromX; x<fromX+SIDE; x++) Colors.setRGB(image,x,yPrime,COLOR);
			for (x=fromX; x<fromX+SIDE; x++) Colors.setRGB(image,x,yPrime+SIDE-1,COLOR);
			for (y=yPrime; y<yPrime+SIDE; y++) Colors.setRGB(image,fromX,y,COLOR);
			for (y=yPrime; y<yPrime+SIDE; y++) Colors.setRGB(image,fromX+SIDE-1,y,COLOR);
		}
		if (fullRight) {
			for (x=fromX+SIDE+HSPACE; x<fromX+SIDE+HSPACE+SIDE; x++) {
				for (y=yPrime; y<yPrime+SIDE; y++) Colors.setRGB(image,x,y,COLOR);
			}
		}
		else {
			for (x=fromX+SIDE+HSPACE; x<fromX+SIDE+HSPACE+SIDE; x++) Colors.setRGB(image,x,yPrime,COLOR);
			for (x=fromX+SIDE+HSPACE; x<fromX+SIDE+HSPACE+SIDE; x++) Colors.setRGB(image,x,yPrime+SIDE-1,COLOR);
			for (y=yPrime; y<yPrime+SIDE; y++) Colors.setRGB(image,fromX+SIDE+HSPACE,y,COLOR);
			for (y=yPrime; y<yPrime+SIDE; y++) Colors.setRGB(image,fromX+SIDE+HSPACE+SIDE-1,y,COLOR);
		}
		
		return fromY+fontMetrics.getAscent();
	}
	
	
	private static final int getNColumns_drawSquares(int valueLeft, int valueRight, BufferedImage image) {
		final int SIDE = (MAXIMALITY_RADIUS)<<2;
		final int HSPACE = SIDE>>1;
		final Graphics2D g2d = image.createGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
		final FontMetrics fontMetrics = g2d.getFontMetrics();
		
		String label;
		if (valueLeft!=-1 && valueRight==-1) label=valueLeft+"";
		else if (valueLeft==-1 && valueRight!=-1) label=valueRight+"";
		else if (valueLeft!=-1 && valueRight!=-1) label=valueLeft+","+valueRight;
		else label="";
		return SIDE+HSPACE+SIDE+HSPACE+fontMetrics.stringWidth(label);
	}

	
	/**
	 * Side of the square drawn by $drawDirectionalSquare()$, in pixels.
	 */
	private static final int SIDE_PIXELS = 3;


	/**
	 * Tries to position a range of length $lengthPixels$ inside the colored header 
	 * $[startPixels..endPixels]$ of an interval of type $type$, so that it does not 
	 * intersect with existing annotation on the interval.
	 *
	 * @param yPixels Y coordinate of the center of the interval, in pixels;
	 * @return the first X coordinate, in pixels, assigned to the beginning of the 
	 * interval, or -1 if the range cannot be positioned.
	 */
	private static final int repositionInterval(int startPixels, int endPixels, int yPixels, int lengthPixels, int type, BufferedImage image) {
		final int MASK = 0x00FFFFFF;
		final int PROBE_SIDE = 10;
		final int expectedColor = Colors.type2color(type);
		boolean found;
		int i, j, x;
		int leftmost, rightmost;
		
		// Finding the leftmost left coordinate
		leftmost=-1;
		for (x=startPixels+2*SIDE_PIXELS; x<=endPixels-2*SIDE_PIXELS; x++) {
			found=false;
			for (i=Math.max(0,x-PROBE_SIDE); i<=Math.min(endPixels,x+PROBE_SIDE); i++) {
				for (j=yPixels-SIDE_PIXELS; j<=yPixels+SIDE_PIXELS; j++) {
					if ((Colors.getRGB(image,i,j)&MASK)!=expectedColor) {
						found=true;
						break;
					}
				}
				if (found) break;
			}
			if (found) continue;
			leftmost=x;
			break;
		}
		
		// Finding the rightmost right coordinate
		rightmost=-1;
		for (x=endPixels-2*SIDE_PIXELS; x>=startPixels+2*SIDE_PIXELS; x--) {
			found=false;
			for (i=Math.max(0,x-PROBE_SIDE); i<=Math.min(endPixels,x+PROBE_SIDE); i++) {
				for (j=yPixels-SIDE_PIXELS; j<=yPixels+SIDE_PIXELS; j++) {
					if ((Colors.getRGB(image,i,j)&MASK)!=expectedColor) {
						found=true;
						break;
					}
				}
				if (found) break;
			}
			if (found) continue;
			rightmost=x;
			break;
		}

		if (leftmost!=-1 && rightmost!=-1 && leftmost+lengthPixels-1<=rightmost) return leftmost;
		else return -1;
	}


	/**
	 * Draws a square with a little line to the left/right to indicate direction.
	 *
	 * @param x,y center of the square;
	 * @param direction TRUE=left, FALSE=right.
	 */
	private static final void drawDirectionalSquare(int x, int y, boolean direction, int intervalType, BufferedImage image, int nRows, int nColumns) {
		final int EXPECTED_COLOR = Colors.type2color(intervalType);
		final int COLOR = Colors.makeLighter(EXPECTED_COLOR);
		final int LINE_MULTIPLE = 4;
		final int MASK = 0x00FFFFFF;
		int i, j;
		int fromX, fromY, toX, toY;
		
		fromX=Math.max(0,x-SIDE_PIXELS); toX=Math.min(x+SIDE_PIXELS,nColumns-1);
		fromY=Math.max(0,y-SIDE_PIXELS); toY=Math.min(y+SIDE_PIXELS,nRows-1);
		for (i=fromX; i<=toX; i++) {
			for (j=fromY; j<=toY; j++) {
				if ((Colors.getRGB(image,i,j)&MASK)==EXPECTED_COLOR) Colors.setRGB(image,i,j,COLOR);
			}
		}
		if (direction) {
			fromX=Math.max(0,x-LINE_MULTIPLE*SIDE_PIXELS); 
			toX=Math.max(0,x-SIDE_PIXELS);
		}
		else {
			fromX=Math.min(x+SIDE_PIXELS,nColumns-1);
			toX=Math.min(x+LINE_MULTIPLE*SIDE_PIXELS,nColumns-1); 
		}
		for (i=fromX; i<=toX; i++) {
			if ((Colors.getRGB(image,i,y)&MASK)==EXPECTED_COLOR) Colors.setRGB(image,i,y,COLOR);
		}
	}
	
	
	/**
	 * Side of the triangle drawn by $drawDirectionalTriangle()$, in pixels.
	 */
	private static final int TRIANGLE_SIDE_PIXELS = (INTERVAL_HEIGHT/3)<<1;
	
	
	/**
	 * @param direction TRUE=left-to-right, FALSE=right-to-left;
	 * @param xCoordinate X position of the head of the horizontal triangle (in pixels);
	 * @param fromY coordinate of the center of the triangle.
	 */
	private static final void drawDirectionalTriangle(boolean direction, int xCoordinate, int fromY, int color, BufferedImage image) {
		int x, y, h;
		final int fromX, toX, length;
		final int HALF_HEIGHT = TRIANGLE_SIDE_PIXELS>>1;
		
		if (direction) {
			toX=xCoordinate;
			fromX=toX-TRIANGLE_SIDE_PIXELS;
		}
		else {
			fromX=xCoordinate;
			toX=fromX+TRIANGLE_SIDE_PIXELS;
		}
		length=toX-fromX+1;
		for (x=fromX; x<=toX; x++) {
			h=Math.round((HALF_HEIGHT*(direction?toX-x:x-fromX))/length);
			for (y=fromY-h; y<=fromY+h; y++) Colors.setRGB(image,x,y,color);
		}
	}
	
	




	
	// ----------------------------- EDGE DRAWING PROCEDURES -----------------------------
	
	/**
	 * Returns a snapshot of interval $[intervalFirst..intervalLast]$ of type 
	 * $intervalType$ taken from the figure in $inputFile$.
	 *
	 * Remark: the procedure assumes $Reads.readLengths$ to have already been loaded.
	 *
	 * @param intervalType if -1, the snapshot contains the whole column $[intervalFirst..
	 * intervalLast]$ of the figure, with no additional formatting;
	 * @param intervalFirst,intervalLast absolute coordinates in the read;
	 * @param xContext number of positions in the read to include in the horizontal 
	 * dimension as a context;
	 * @param yContext number of alignments to include in the vertical dimension as a 
	 * context;
	 * @param out output array; if not null, it contains: 0=the horizontal and 1=the 
	 * vertical coordinate, in the output image, of the top-left point of the interval;
	 * 2=the vertical coordinate, in the output image, of the of the bottom-most point of
	 * the interval box;
	 * @return NULL if the interval could not be found.
	 */
	public static final BufferedImage getIntervalContext(String inputFile, int read, int intervalFirst, int intervalLast, int intervalType, int xContext, int yContext, int[] out) throws IOException {
		final int MASK = 0x00FFFFFF;
		final double FRACTION = 0.99;  // Arbitrary
		final int VERTICAL_SPACE = 10;  // Arbitrary. In pixels.
		final int DARKER_ITERATIONS = 3;
		final int INTERVAL_COLOR = Colors.type2color(intervalType);
		boolean found;
		int i, x, y;
		int inputImageWidth, inputImageHeight, firstY, lastY, sum;
		int fromX, toX, toY1, toY1Prime, fromY2, toY2;
		BufferedImage inputImage, outputImage;
		Color color;
		
		inputImage=ImageIO.read(new File(inputFile));
		inputImageWidth=inputImage.getWidth();
		inputImageHeight=inputImage.getHeight();
		fromX=Math.max((intervalFirst-xContext)/QUANTUM,0);	
		toX=Math.min((intervalLast+xContext)/QUANTUM,Reads.getReadLength(read)/QUANTUM-1);
		if (out!=null) out[0]=intervalFirst/QUANTUM-fromX;
		
		// Finding header
		toY1=-1;
		for (y=0; y<inputImageHeight; y++) {
			found=false;
			for (x=0; x<inputImageWidth; x++) {			
				if ((inputImage.getRGB(x,y)&MASK)==(Colors.COLOR_BACKGROUND_LIGHT&MASK)) {
					found=true;
					break;
				}
			}
			if (found) {
				toY1=y-1;
				break;
			}
		}
		
		// Finding interval
		if (intervalType==-1) {
			fromY2=toY1+1;
			firstY=fromY2+1;
			lastY=inputImageHeight-2;
			toY2=inputImageHeight-1;
		}
		else {
			firstY=Math.POSITIVE_INFINITY; lastY=-1;
			for (y=inputImageHeight-1; y>toY1; y--) {
				sum=0;
				for (x=intervalFirst/QUANTUM; x<=intervalLast/QUANTUM; x++) {
					if ((inputImage.getRGB(x,y)&MASK)==(INTERVAL_COLOR&MASK)) sum++;
				}
				if (sum<FRACTION*(intervalLast-intervalFirst)/QUANTUM) continue;
				if (lastY==-1) lastY=y;
				if (y<firstY) firstY=y;
			}
			if (lastY==-1 || firstY==Math.POSITIVE_INFINITY) {
				System.err.println("WARNING: cannot find interval in image?! read="+read+" ["+intervalFirst+".."+intervalLast+"] type="+intervalType);
				return null;
			}
			toY2=Math.min(lastY+yContext,inputImageHeight-1);
			toY1Prime=-1;
			for (y=firstY-1; y>toY1; y--) {
				found=false;
				for (x=0; x<inputImageWidth; x++) {			
					if ((inputImage.getRGB(x,y)&MASK)==(Colors.COLOR_BACKGROUND_LIGHT&MASK)) {
						found=true;
						break;
					}
				}
				if (found) {
					toY1Prime=y;
					break;
				}
			}
			fromY2=Math.max(firstY-yContext,toY1Prime+1);
		}
		
		// Darkening output rectangles
		if (intervalType!=-1) {
			for (y=0; y<=toY1; y++) {
				for (x=fromX; x<=toX; x++) {
					color = new Color(inputImage.getRGB(x,y));
					for (i=0; i<DARKER_ITERATIONS; i++) color=color.darker();
					inputImage.setRGB(x,y,color.getRGB());
				}
			}
			for (y=fromY2; y<firstY; y++) {
				for (x=fromX; x<=toX; x++) {
					color = new Color(inputImage.getRGB(x,y));
					for (i=0; i<DARKER_ITERATIONS; i++) color=color.darker();
					inputImage.setRGB(x,y,color.getRGB());
				}
			}
			for (y=lastY+1; y<=toY2; y++) {
				for (x=fromX; x<=toX; x++) {
					color = new Color(inputImage.getRGB(x,y));
					for (i=0; i<DARKER_ITERATIONS; i++) color=color.darker();
					inputImage.setRGB(x,y,color.getRGB());
				}
			}
		}
		
		// Building output image
		outputImage = new BufferedImage(toX-fromX+1,toY1+1+VERTICAL_SPACE+(toY2-fromY2+1),BufferedImage.TYPE_INT_RGB);
		for (y=0; y<=toY1; y++) {
			for (x=fromX; x<=toX; x++) outputImage.setRGB(x-fromX,y,inputImage.getRGB(x,y));
		}
		for (y=toY1+1; y<=toY1+VERTICAL_SPACE; y++) {
			for (x=fromX; x<=toX; x++) outputImage.setRGB(x-fromX,y,Colors.COLOR_BACKGROUND);
		}
		for (y=fromY2; y<=toY2; y++) {
			for (x=fromX; x<=toX; x++) outputImage.setRGB(x-fromX,(toY1+VERTICAL_SPACE+1)+(y-fromY2),inputImage.getRGB(x,y));
		}
		if (out!=null) {
			out[1]=toY1+VERTICAL_SPACE+firstY-fromY2;
			out[2]=toY1+VERTICAL_SPACE+lastY-fromY2;
		}
		return outputImage;
	}
	
	
	/**
	 * Draws details on every edge incident to the interval identified by 
	 * $read,type,start,end$, in one output file per edge. The procedure assumes that 
	 * pictures of the factorization of all related reads are already located in directory
	 * $snippetsPath$.
	 *
	 * @param quantum number of bases per pixel in the image;
	 * @param histogram temporary space, of size at least equal to the length of a longest
	 * read;
	 * @param degreeHistogram,clusteringCoefficientHistogram histograms computed on all 
	 * nodes in the graph;
	 * @param tmpArray* temporary space, of size at least equal to the max degree of a 
	 * node in the graph;
	 * @param sameKernel draws only edges that are of overlap type and are incident to 
	 * nodes in the same kernel.
	 */
	public static final void drawEdges(int read, int type, int start, int end, String snippetsPath, String outputPrefix, String outputExtension, int quantum, int[] histogram, int[] degreeHistogram, int[] clusteringCoefficientHistogram, int[] tmpArray1, int[] tmpArray2, boolean sameKernel) throws IOException {
		final int MIN_ALIGNMENT_LENGTH = 1000;
		final int X_CONTEXT = 2000;
		final int Y_CONTEXT = 400;
		final int SNIPPETS_Y_BORDER = 200;
		final int SNIPPETS_X_BORDER = 50;
		final int HEADER_SIZE = (INTERVAL_HEIGHT)<<1;
		final int SEPARATOR_SIZE = 10;
		final int HISTOGRAM_HEIGHT = (HISTOGRAM_ROWS)<<1;
		final int HISTOGRAM_SPACE = SMALL_SPACE;
		final int SNIPPET_SQUARE_SIZE = 40;
		final int HISTOGRAMS_SQUARE_SIZE = SNIPPET_SQUARE_SIZE/2;
		final int HISTOGRAMS_PANEL_WIDTH = 386;  // In pixels
		final int HISTOGRAMS_PANEL_HEIGHT = 100;  // In pixels
		final int HISTOGRAMS_PANEL_VERTICAL_OFFSET = 80;  // In pixels
		final int HISTOGRAMS_PANEL_HORIZONTAL_OFFSET = 220;  // In pixels
		final int HISTOGRAMS_PANEL_HORIZONTAL_BORDER = 10;  // In pixels
		final int HISTOGRAMS_PANEL_VERTICAL_BORDER = 10;  // In pixels
		int i, j, x, y;
		int nNeighbors, nRows, nColumns, middleColumns, lastY, fromX, fromY;
		int referenceLength, nRowsReference, nColumnsReference;
		int referencePointX, referencePointY, otherPointX, otherPointY;
		int lastCell, degree, otherDegree, offset, sideFromY, sideToY;
		int kernel;
		double clusteringCoefficient, otherClusteringCoefficient;
		String label;
		IntervalGraph.Node neighbor;
		IntervalGraph.Edge edge;
		Graphics2D graphics;
		FontMetrics fontMetrics;
		BufferedImage snippetLeft, snippetRight, out;
		int[] coordinates = new int[3];
		int[] neighborType1 = new int[8];
		int[] neighborType2 = new int[8];
		int[] edgeType1 = new int[8];
		int[] edgeType2 = new int[8];
		
		for (i=0; i<IntervalGraph.nNodes; i++) {
			if (IntervalGraph.nodesArray[i].read!=read || IntervalGraph.nodesArray[i].type!=type || IntervalGraph.nodesArray[i].start!=start || IntervalGraph.nodesArray[i].end!=end) continue;
			referenceLength=IntervalGraph.nodesArray[i].length();
			nNeighbors=IntervalGraph.nNeighbors[i];
			snippetLeft=getIntervalContext(snippetsPath+"/read"+(read+1)+".png",read,start,end,type,X_CONTEXT,Y_CONTEXT,coordinates);
			if (snippetLeft==null) {
				if (type==Constants.INTERVAL_PERIODIC) snippetLeft=getIntervalContext(snippetsPath+"/read"+(read+1)+".png",read,start,end,Constants.INTERVAL_ALIGNMENT,X_CONTEXT,Y_CONTEXT,coordinates);
				if (snippetLeft==null) return;
			}
			referencePointX=SNIPPETS_X_BORDER+coordinates[0]; 
			referencePointY=SNIPPETS_Y_BORDER+coordinates[1];
			sideFromY=coordinates[1];
			sideToY=coordinates[2];
			nRowsReference=snippetLeft.getHeight();
			nColumnsReference=snippetLeft.getWidth();
			lastCell=IntervalGraph.getEdgeHistogram(i,histogram,quantum,false);
			IntervalGraph.getEdgeStats(i,neighborType1,edgeType1);
			degree=IntervalGraph.getDegree(i,true);
			clusteringCoefficient=IntervalGraph.getClusteringCoefficient(i,tmpArray1,tmpArray2);
			kernel=IntervalGraph.nodesArray[i].lastKernel==0?IntervalGraph.nodesArray[i].kernels[0]:-1;
			
			// Building images
			for (j=0; j<nNeighbors; j++) {
				neighbor=IntervalGraph.nodesArray[IntervalGraph.neighbors[i][j].getTo(IntervalGraph.nodesArray[i].nodeID)];
				if (sameKernel && (IntervalGraph.neighbors[i][j].overlap==-1 || neighbor.lastKernel!=0 || neighbor.kernels[0]!=kernel)) continue;
				
				// Allocating the image
				snippetRight=getIntervalContext(snippetsPath+"/read"+(neighbor.read+1)+".png",neighbor.read,neighbor.start,neighbor.end,neighbor.type,X_CONTEXT,Y_CONTEXT,coordinates);
				if (snippetRight==null) {
					if (neighbor.type==Constants.INTERVAL_PERIODIC) snippetRight=getIntervalContext(snippetsPath+"/read"+(neighbor.read+1)+".png",neighbor.read,neighbor.start,neighbor.end,Constants.INTERVAL_ALIGNMENT,X_CONTEXT,Y_CONTEXT,coordinates);
					if (snippetRight==null) continue;
				}
				nRows=SNIPPETS_Y_BORDER+Math.max(nRowsReference,snippetRight.getHeight())+SNIPPETS_X_BORDER;
				middleColumns=(referenceLength+2*(neighbor.length()-MIN_ALIGNMENT_LENGTH))/quantum;
				otherPointX=SNIPPETS_X_BORDER+nColumnsReference+middleColumns+coordinates[0];
				otherPointY=SNIPPETS_Y_BORDER+coordinates[1];
				nColumns=SNIPPETS_X_BORDER+nColumnsReference+middleColumns+snippetRight.getWidth()+SNIPPETS_X_BORDER;
				out = new BufferedImage(nColumns,nRows,BufferedImage.TYPE_INT_RGB);
				for (y=0; y<nRows; y++) {
					for (x=0; x<nColumns; x++) out.setRGB(x,y,Colors.COLOR_BACKGROUND);
				}
				graphics=out.createGraphics();
				graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
				graphics.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
				fontMetrics=graphics.getFontMetrics();
				
				// Filling the image: left snippet.
				drawFrame(graphics,SNIPPETS_X_BORDER,SNIPPETS_Y_BORDER,nColumnsReference,nRowsReference,SNIPPET_SQUARE_SIZE,true,true,false,-1,sideFromY,sideToY,false);
				label="READ "+read+" ["+start+".."+end+"] // NODE_ID "+IntervalGraph.nodesArray[i].nodeID;
				graphics.setColor(new Color(Colors.COLOR_TEXT));
				graphics.drawString(label,SNIPPETS_X_BORDER+nColumnsReference-fontMetrics.stringWidth(label),SNIPPETS_Y_BORDER-fontMetrics.getAscent());
				for (y=0; y<nRowsReference; y++) {
					for (x=0; x<nColumnsReference; x++) out.setRGB(SNIPPETS_X_BORDER+x,SNIPPETS_Y_BORDER+y,snippetLeft.getRGB(x,y));
				}
				// Right snippet
				drawFrame(graphics,SNIPPETS_X_BORDER+nColumnsReference+middleColumns,SNIPPETS_Y_BORDER,snippetRight.getWidth(),snippetRight.getHeight(),SNIPPET_SQUARE_SIZE,true,true,false,-1,coordinates[1],coordinates[2],true);
				label="NODE_ID "+neighbor.nodeID+" // READ "+neighbor.read+" ["+neighbor.start+".."+neighbor.end+"]";
				graphics.setColor(new Color(Colors.COLOR_TEXT));
				graphics.drawString(label,SNIPPETS_X_BORDER+nColumnsReference+middleColumns,SNIPPETS_Y_BORDER-fontMetrics.getAscent());
				for (y=0; y<snippetRight.getHeight(); y++) {
					for (x=0; x<snippetRight.getWidth(); x++) out.setRGB(SNIPPETS_X_BORDER+nColumnsReference+middleColumns+x,SNIPPETS_Y_BORDER+y,snippetRight.getRGB(x,y));
				}
				// Edge details in the middle
				fromX=(middleColumns-referenceLength/quantum)>>1;
				fromY=Math.min(nColumnsReference,snippetRight.getHeight())/3;
				IntervalGraph.neighbors[i][j].draw(out,IntervalGraph.nodesArray[i].nodeID,SNIPPETS_X_BORDER+nColumnsReference+fromX,SNIPPETS_Y_BORDER+fromY,referencePointX,referencePointY,otherPointX,otherPointY,SNIPPETS_Y_BORDER,SNIPPETS_Y_BORDER+nRowsReference,SNIPPETS_Y_BORDER,SNIPPETS_Y_BORDER+snippetRight.getHeight());
				drawEdgeHistogram(histogram,lastCell,SNIPPETS_X_BORDER+nColumnsReference+fromX,SNIPPETS_Y_BORDER+fromY-HISTOGRAM_SPACE,HISTOGRAM_HEIGHT,IntervalGraph.nNeighbors[i],quantum,out);
				
				// Histograms of left snippet
				drawFrame(graphics,SNIPPETS_X_BORDER+nColumnsReference-HISTOGRAMS_PANEL_WIDTH,SNIPPETS_Y_BORDER-HISTOGRAMS_PANEL_VERTICAL_OFFSET-HISTOGRAMS_PANEL_HEIGHT,HISTOGRAMS_PANEL_WIDTH,HISTOGRAMS_PANEL_HEIGHT,HISTOGRAMS_SQUARE_SIZE,true,false,true,Colors.COLOR_BACKGROUND,-1,-1,false);
				drawEdgeStats(neighborType1,edgeType1,out,SNIPPETS_X_BORDER+nColumnsReference-HISTOGRAMS_PANEL_WIDTH+HISTOGRAMS_PANEL_HORIZONTAL_BORDER,SNIPPETS_Y_BORDER-HISTOGRAMS_PANEL_VERTICAL_OFFSET-HISTOGRAMS_PANEL_VERTICAL_BORDER,HISTOGRAMS_PANEL_HEIGHT-HISTOGRAMS_PANEL_VERTICAL_BORDER*2);
				drawDegreeHistograms(out,graphics,fontMetrics,SNIPPETS_X_BORDER+nColumnsReference-HISTOGRAMS_PANEL_WIDTH+HISTOGRAMS_PANEL_HORIZONTAL_OFFSET,SNIPPETS_Y_BORDER-HISTOGRAMS_PANEL_VERTICAL_OFFSET-HISTOGRAMS_PANEL_VERTICAL_BORDER,degreeHistogram,clusteringCoefficientHistogram,HISTOGRAMS_PANEL_HEIGHT-HISTOGRAMS_PANEL_VERTICAL_BORDER*2,degree,tmpArray1.length,clusteringCoefficient);
				
				// Histograms of right snippet
				IntervalGraph.getEdgeStats(neighbor.nodeID,neighborType2,edgeType2);
				otherDegree=IntervalGraph.getDegree(neighbor.nodeID,true);
				otherClusteringCoefficient=IntervalGraph.getClusteringCoefficient(neighbor.nodeID,tmpArray1,tmpArray2);
				offset=SNIPPETS_X_BORDER+nColumnsReference+middleColumns;
				drawFrame(graphics,offset,SNIPPETS_Y_BORDER-HISTOGRAMS_PANEL_VERTICAL_OFFSET-HISTOGRAMS_PANEL_HEIGHT,HISTOGRAMS_PANEL_WIDTH,HISTOGRAMS_PANEL_HEIGHT,HISTOGRAMS_SQUARE_SIZE,true,false,true,Colors.COLOR_BACKGROUND,-1,-1,false);
				drawEdgeStats(neighborType2,edgeType2,out,offset+HISTOGRAMS_PANEL_HORIZONTAL_BORDER,SNIPPETS_Y_BORDER-HISTOGRAMS_PANEL_VERTICAL_OFFSET-HISTOGRAMS_PANEL_VERTICAL_BORDER,HISTOGRAMS_PANEL_HEIGHT-HISTOGRAMS_PANEL_VERTICAL_BORDER*2);
				drawDegreeHistograms(out,graphics,fontMetrics,offset+HISTOGRAMS_PANEL_HORIZONTAL_OFFSET,SNIPPETS_Y_BORDER-HISTOGRAMS_PANEL_VERTICAL_OFFSET-HISTOGRAMS_PANEL_VERTICAL_BORDER,degreeHistogram,clusteringCoefficientHistogram,HISTOGRAMS_PANEL_HEIGHT-HISTOGRAMS_PANEL_VERTICAL_BORDER*2,otherDegree,tmpArray1.length,otherClusteringCoefficient);
				
				ImageIO.write(out,outputExtension,new File(outputPrefix+(IntervalGraph.neighbors[i][j].supplement?"suppl-":"")+read+"-"+neighbor.read+"-"+i+"-"+j+"."+outputExtension));
			}
			break;
		}
	}
	
	
	/**
	 * Draws a rectangular frame, to the right and down WRT point $(x,y)$.
	 * 
	 * @param squareSize size of the squares at the corners of the frame;
	 * @param horizontalLines TRUE=draw the horizontal lines of the frame;
	 * @param sideFromY,sideToY (relative to $y$) if non-negative, the vertical interval 
	 * $[fromY..toY]$ is highlighted, with a left or right border as decided by 
	 * $leftOrRight$.
	 */
	private static final void drawFrame(Graphics2D graphics, int x, int y, int width, int height, int squareSize, boolean horizontalLines, boolean verticalLines, boolean fill, int fillColor, int sideFromY, int sideToY, boolean leftOrRight) {
		final int STROKE = 4;  // In pixels
		final double STROKE_RATE = 1.2;
		final int STROKE_OFFSET = (int)(STROKE*STROKE_RATE);
		final int SQUARE_SIZE = Math.min(squareSize,width/6);  // In pixels
		final int COLOR_HIGHLIGHT = Colors.COLOR_MAXIMAL_DESATURATED;
		final int x1 = x;
		final int y1 = y;
		final int x2 = x1+width;
		final int y2 = y1;
		final int x3 = x2;
		final int y3 = y+height;
		final int x4 = x1;
		final int y4 = y3;
		
		graphics.setColor(new Color(Colors.FRAME_COLOR)); 
		graphics.setStroke(new BasicStroke(STROKE));
		if (horizontalLines) {
			graphics.drawLine(x1,y1,x2,y2);
			graphics.drawLine(x3,y3,x4,y4);
		}
		if (verticalLines) {
			graphics.drawLine(x1,y1,x4,y4);
			graphics.drawLine(x2,y2,x3,y3);
		}
		graphics.fillRect(x1-STROKE_OFFSET,y1-STROKE_OFFSET,SQUARE_SIZE,SQUARE_SIZE);
		graphics.fillRect(x2+STROKE_OFFSET-SQUARE_SIZE,y2-STROKE_OFFSET,SQUARE_SIZE,SQUARE_SIZE);
		graphics.fillRect(x3+STROKE_OFFSET-SQUARE_SIZE,y3+STROKE_OFFSET-SQUARE_SIZE,SQUARE_SIZE,SQUARE_SIZE);
		graphics.fillRect(x4-STROKE_OFFSET,y4+STROKE_OFFSET-SQUARE_SIZE,SQUARE_SIZE,SQUARE_SIZE);
		if (sideFromY>=0) {
			graphics.setColor(new Color(COLOR_HIGHLIGHT)); 
			graphics.fillRect(leftOrRight?x1-STROKE_OFFSET:x2,y+sideFromY,STROKE_OFFSET,sideToY-sideFromY);
		}
		if (fill) {
			graphics.setColor(new Color(fillColor));
			graphics.fillRect(x,y,width,height);
		}
	}

	
	/**
	 * Prints $out[0..last]$ starting at pixel $(fromX,fromY)$ of $image$ and going right 
	 * and up. Every position of $out$ is drawn with a single horizontal pixel, but with
	 * possibly multiple vertical pixels.
	 *
	 * @param height maximum height of the histogram (in pixels).
	 */
	private static final void drawEdgeHistogram(int[] out, int last, int fromX, int fromY, int height, int max, int quantum, BufferedImage image) {
		int i, j;
		int p;
		
		for (i=0; i<=last; i++) {
			p=(out[i]*height)/max;
			for (j=0; j<p; j++) Colors.setRGB(image,fromX+i,fromY-j,Colors.COLOR_HISTOGRAM);
			if (out[i]>0) Colors.setRGB(image,fromX+i,fromY-p,Colors.COLOR_MAXIMAL_DESATURATED);
		}
	}
	
	
	/**
	 * Draws the statistics on the edges incident to a given node, produced by 
	 * $IntervalGraph.getEdgeStats()$.
	 *
	 * @param x,y histograms are drawn to the right and up WRT this point;
	 * @param maxHeight maximum height of the histograms (in pixels).
	 */
	private static final void drawEdgeStats(int[] neighborType, int[] edgeType, BufferedImage image, int x, int y, int maxHeight) {
		final int COLOR_HIGHLIGHT = Colors.COLOR_MAXIMAL_DESATURATED;
		final int COLOR_EDGETYPE = Colors.FRAME_COLOR;
		final int HISTOGRAM_THICKNESS = 5;  // In pixels
		final int TIP_THICKNESS = 2;  // In pixels
		final int DELTA_X = HISTOGRAM_THICKNESS*2;  // In pixels
		int p, q;
		int max, fromX, fromY, toX, toY, color;
		
		// Drawing $neighborType$
		max=Math.max(neighborType,7);
		fromY=y; toY=y-maxHeight*neighborType[0]/max;
		fromX=x; toX=x+HISTOGRAM_THICKNESS;
		color=Colors.type2color(Constants.INTERVAL_ALIGNMENT);
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,color);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		toY=y-maxHeight*neighborType[1]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		color=Colors.type2color(Constants.INTERVAL_DENSE_PREFIX);
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,color);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		toY=y-maxHeight*neighborType[2]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		color=Colors.type2color(Constants.INTERVAL_DENSE_SUFFIX);
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,color);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		toY=y-maxHeight*neighborType[3]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		color=Colors.type2color(Constants.INTERVAL_DENSE_PREFIXSUFFIX);
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,color);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		toY=y-maxHeight*neighborType[4]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		color=Colors.type2color(Constants.INTERVAL_DENSE_SUBSTRING);
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,color);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		toY=y-maxHeight*neighborType[5]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		color=Colors.type2color(Constants.INTERVAL_DENSE_SINGLEDELETION);
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,color);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		toY=y-maxHeight*neighborType[6]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		color=Colors.type2color(Constants.INTERVAL_PERIODIC);
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,color);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		toY=y-maxHeight*neighborType[7]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		color=Colors.type2color(Constants.INTERVAL_PERIODIC);
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,color);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		
		// Drawing $edgeType$
		max=Math.max(edgeType,7);
		fromY=y; toY=y-maxHeight*edgeType[0]/max;
		fromX+=DELTA_X*2; toX=fromX+HISTOGRAM_THICKNESS;
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_EDGETYPE);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		fromY=y; toY=y-maxHeight*edgeType[1]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_EDGETYPE);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		
		fromY=y; toY=y-maxHeight*edgeType[2]/max;
		fromX+=DELTA_X*2; toX=fromX+HISTOGRAM_THICKNESS;
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_EDGETYPE);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		fromY=y; toY=y-maxHeight*edgeType[3]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_EDGETYPE);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		fromY=y; toY=y-maxHeight*edgeType[4]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_EDGETYPE);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		
		fromY=y; toY=y-maxHeight*edgeType[5]/max;
		fromX+=DELTA_X*2; toX=fromX+HISTOGRAM_THICKNESS;
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_EDGETYPE);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		fromY=y; toY=y-maxHeight*edgeType[6]/max;
		fromX+=DELTA_X; toX=fromX+HISTOGRAM_THICKNESS;
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_EDGETYPE);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
		
		fromY=y; toY=y-maxHeight*edgeType[7]/max;
		fromX+=DELTA_X*2; toX=fromX+HISTOGRAM_THICKNESS;
		for (p=fromX; p<=toX; p++) {
			for (q=fromY; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_EDGETYPE);
		}
		for (p=fromX; p<=toX; p++) {
			for (q=0; q<TIP_THICKNESS; q++) Colors.setRGB(image,p,toY-q,COLOR_HIGHLIGHT);
		}
	}

	
	/**
	 * Draws the histograms produced by $IntervalGraph.getDegreeHistogram()$ and
	 * $IntervalGraph.getClusteringCoefficientHistogram()$.
	 *
	 * @param x,y histograms are drawn to the right and up WRT this point;
	 * @param maxHeight maximum height of the histograms (in pixels).
	 */
	private static final void drawDegreeHistograms(BufferedImage image, Graphics graphics, FontMetrics fontMetrics, int x, int y, int[] degreeHistogram, int[] clusteringCoefficientHistogram, int maxHeight, int selectedDegree, int maxDegree, double selectedClusteringCoefficient) {
		final int COLOR_HISTOGRAM_LOCAL = Colors.COLOR_REFERENCE;
		final int COLOR_HISTOGRAM_HIGHLIGHT = Colors.COLOR_HISTOGRAM;
		final int HISTOGRAM_THICKNESS = 1;  // In pixels
		final int DELTA_X = HISTOGRAM_THICKNESS;  // In pixels
		final int HORIZONTAL_SPACE = 15;  // In pixels
		final double DEGREE_QUANTUM = ((double)maxDegree)/degreeHistogram.length;
		final double CLUSTERING_QUANTUM = 1.0/clusteringCoefficientHistogram.length;
		int i, p, q;
		int last, max, fromX, fromXPrime, toX, toY, color;
		String label;
		
		// Degree
		last=degreeHistogram.length-1;
		max=Math.max(degreeHistogram,last);
		fromX=x;
		for (i=0; i<=last; i++) {
			toY=y-maxHeight*degreeHistogram[i]/max;
			toX=fromX+HISTOGRAM_THICKNESS;
			for (p=fromX; p<=toX; p++) {
				for (q=y; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_HISTOGRAM_LOCAL);
			}
			fromX+=DELTA_X;
		}
		i=(int)(selectedDegree/DEGREE_QUANTUM);
		if (i>=degreeHistogram.length) i=degreeHistogram.length-1;
		toY=y-maxHeight*degreeHistogram[i]/max;
		fromXPrime=x+i*DELTA_X;
		toX=fromXPrime+HISTOGRAM_THICKNESS;
		for (p=fromXPrime; p<=toX; p++) {
			for (q=y; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_HISTOGRAM_HIGHLIGHT);
		}
		graphics.setFont(new Font(DEFAULT_FONT,Font.PLAIN,15));
		graphics.setColor(new Color(COLOR_HISTOGRAM_HIGHLIGHT));
		label=""+selectedDegree;
		graphics.drawString(label,x+degreeHistogram.length*DELTA_X/2-fontMetrics.stringWidth(label)/2,y-maxHeight+fontMetrics.getAscent()/2);
		
		// Clustering coefficient
		last=clusteringCoefficientHistogram.length-1;
		max=Math.max(clusteringCoefficientHistogram,last);
		fromX+=HORIZONTAL_SPACE;
		fromXPrime=fromX;
		for (i=0; i<=last; i++) {
			toY=y-maxHeight*clusteringCoefficientHistogram[i]/max;
			toX=fromX+HISTOGRAM_THICKNESS;
			for (p=fromX; p<=toX; p++) {
				for (q=y; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_HISTOGRAM_LOCAL);
			}
			fromX+=DELTA_X;
		}
		if (selectedClusteringCoefficient>=0) {
			i=(int)(selectedClusteringCoefficient/CLUSTERING_QUANTUM);
			if (i>=clusteringCoefficientHistogram.length) i=clusteringCoefficientHistogram.length-1;
			toY=y-maxHeight*clusteringCoefficientHistogram[i]/max;
			fromXPrime+=i*DELTA_X;
			toX=fromXPrime+HISTOGRAM_THICKNESS;
			for (p=fromXPrime; p<=toX; p++) {
				for (q=y; q>=toY; q--) Colors.setRGB(image,p,q,COLOR_HISTOGRAM_HIGHLIGHT);
			}
		}
	}
	

}