package de.mpi_cbg.revant.consensus;

import java.io.*;
import java.util.*;
import java.text.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Colors;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Factorize;
import de.mpi_cbg.revant.factorize.ReadA;
import de.mpi_cbg.revant.factorize.Alignments;


/**
 * Used in the experiment in which we mapped my lungfish modules to 1x of the reads.
 * Builds the suffix trie of adjacent occurrences of repeats (no difference between a full
 * copy and a substring of a repeat).
 */
public class MappingExcerciseUnique {
	
	static NumberFormat formatter = NumberFormat.getInstance();
	
	/**
	 * java MappingExcerciseUnique nAlignments rootDir 
	 */
	public static void main(String[] args) throws Exception {
		int N_ALIGNMENTS = Integer.parseInt(args[0]);
		String ROOT_DIR = args[1];
		int N_BLOCKS = 120;
		int N_READS = 1868764;
		int MAX_READ_LENGTH = 500000;  // Arbitrary
		String LENGTHS_FILE = ROOT_DIR+"/readLengths.txt";
		int MAX_ALIGNMENTS_PER_READ = 10000;  // Arbitrary
		String ALIGNMENTS_FILE = ROOT_DIR+"/alignments-all-clean.txt";
		String QUALITY_THRESHOLDS_FILE = ROOT_DIR+"/qualityThresholds.txt";
		final int IDENTITY_THRESHOLD = IO.quantum;
		
		int i, j, k;
		int top, currentAlignment, currentTag, currentEnd, start, tag, readID;
		Node root, currentNode;
		int[] alignmentStack = new int[MAX_ALIGNMENTS_PER_READ];
		Node[] nodeStack = new Node[MAX_ALIGNMENTS_PER_READ];
		formatter.setMaximumFractionDigits(2);
		
		Factorize.loadInput(new String[] {N_READS+"",MAX_READ_LENGTH+"",LENGTHS_FILE,N_ALIGNMENTS+"",ALIGNMENTS_FILE,QUALITY_THRESHOLDS_FILE,"null","1","3","0","1000","1000","0","null","null","null","null","null","null","null","0","null"});
		ReadA.allocateMemory(MAX_READ_LENGTH,MAX_ALIGNMENTS_PER_READ);
		root = new Node(-1); i=0;
		while (i<Alignments.nAlignments) {
			ReadA.initialize(i,false,false);
			readID=ReadA.id;
			for (j=0; j<=ReadA.lastSortedAlignment; j++) {
				currentTag=Alignments.alignments[j][1]-1;
				currentNode=root.addChild(currentTag,readID);
				top=0; alignmentStack[0]=j; nodeStack[0]=currentNode;
				while (top>=0) {
					currentAlignment=alignmentStack[top];
					currentNode=nodeStack[top];
					top--;
					currentEnd=Alignments.alignments[currentAlignment][4];
					for (k=currentAlignment+1; k<=ReadA.lastSortedAlignment; k++) {
						start=Alignments.alignments[k][3];
						if (start>currentEnd+IDENTITY_THRESHOLD) break;
						if (start<currentEnd-IDENTITY_THRESHOLD) continue;
						tag=Alignments.alignments[k][1]-1;
						top++;
						alignmentStack[top]=k;
						nodeStack[top]=currentNode.addChild(tag,readID);
					}
				}
			}
			i=ReadA.lastAlignment+1;
		}
		
		// Printing
		System.out.println("digraph G {\n");
		System.out.println("bgcolor=\""+Colors.COLOR_DOT_BACKGROUND+"\";\n");
		root.toDot("");
		System.out.println("}");
	}

	
	private static class Node implements Comparable {
		public static final int ORDER_ID = 0;
		public static final int ORDER_COUNT = 1;
		public static int order;
		public int id;
		public Node[] children;
		public int lastChild;
		public int count;
		public String dotPrefix;
		public int previousRead;
		
		// Temporary space
		private static Node tmpNode = new Node(-1);
		private static int dotID;
		
		public Node(int id) {
			this.id=id;
			count=0;
			children=null;
			lastChild=-1;
			previousRead=-1;
		}
		
		
		public int compareTo(Object other) {
			Node otherNode = (Node)other;
			if (order==ORDER_ID) {
				if (id<otherNode.id) return -1;
				else if (id>otherNode.id) return 1;
			}
			else if (order==ORDER_COUNT) {  // Non-increasing
				if (count>otherNode.count) return -1;
				else if (count<otherNode.count) return 1;
			}
			return 0;
		}
		
		
		public Node addChild(int id, int read) {
			final int CAPACITY = 10;
			int i;
			final int previousOrder = Node.order;
			Node newChild;
			
			// Checking if the child exists
			tmpNode.id=id;
			if (lastChild==0) {
				if (children[0].id==id) {
					if (children[0].previousRead!=read) {
						children[0].previousRead=read;
						children[0].count++;
					}
					return children[0];
				}
			}
			else if (lastChild>0) {
				Node.order=Node.ORDER_ID;
				i=Arrays.binarySearch(children,0,lastChild+1,tmpNode);
				Node.order=previousOrder;
				if (i>=0) {
					if (children[i].previousRead!=read) {
						children[i].previousRead=read;
						children[i].count++;
					}
					return children[i];
				}
			}
			
			// Adding a new child
			lastChild++;
			if (children==null || children.length==lastChild) {
				Node[] newChildren = new Node[(children==null||children.length==0)?CAPACITY:children.length<<1];
				if (children!=null && children.length>0) System.arraycopy(children,0,newChildren,0,children.length);
				children=newChildren;
			}
			newChild = new Node(id);
			newChild.previousRead=read;
			newChild.count=1;
			children[lastChild]=newChild;
			if (lastChild>0) {
				Node.order=ORDER_ID;
				Arrays.sort(children,0,lastChild+1);
				Node.order=previousOrder;
			}
			return newChild;
		}
		
		
		/**
		 * @return the number of leaves in the subtree.
		 */
		public int nLeaves() {
			int i, out;
			
			if (lastChild==-1) return 1;
			out=0;
			for (i=0; i<=lastChild; i++) out+=children[i].nLeaves();
			return out;
		}
		
		
		/**
		 * @return the sum of the count values of all leaves in the subtree. 
		 * The procedure sets to this value the count field of every internal node.
		 */
		public int getFrequency() {
			int i;
			
			if (lastChild==-1) return count;
			for (i=0; i<=lastChild; i++) count+=children[i].getFrequency();
			return count;
		}

		
		public void toDot(String dotPrefix) {
			final double PEN_SCALE = 10.0;
			final double NODE_SCALE = 10.0;
			final double FONT_SCALE = 100.0;
			final int oldID = dotID;
			final String fillColor = Colors.type2color_dot(Constants.INTERVAL_ALIGNMENT);
			final int nodeSize = Math.max(1,(int)((1.0/Math.sqrt(count))*NODE_SCALE));
			final int fontSize = Math.max(1,(int)((1.0/Math.sqrt(count))*FONT_SCALE));
			int edgeSize;
			
			System.out.println(dotID+" [shape=\"circle\",style=\"filled\",fixedsize=\"true\",width="+nodeSize+",fontsize=\""+fontSize+"\",label=\""+dotPrefix+" \n "+count+"\",color=\""+fillColor+"\",fontcolor=\""+Colors.COLOR_DOT_TEXT+"\",fontname=\"Helvetica\"];\n");
			this.dotPrefix=dotPrefix;
			for (int i=0; i<=lastChild; i++) {
				dotID++;
				edgeSize=Math.max(1,(int)((1.0/Math.sqrt(children[i].count))*PEN_SCALE));
				System.out.println(oldID+" -> "+dotID+" [dir=\"both\",arrowhead=\"none\",arrowtail=\"none\",penwidth=\""+edgeSize+"\",color=\""+Colors.COLOR_DOT_EDGE+"\"];\n");
				children[i].toDot((dotPrefix.length()==0?"":dotPrefix+"--")+children[i].id);
			}
		}
	}

}