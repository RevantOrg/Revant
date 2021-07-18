package de.mpi_cbg.revant.graphics;

import java.util.Arrays;
import java.io.*;
import java.text.*;
import de.mpi_cbg.revant.util.Math;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Colors;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.factorize.Alignments;
import de.mpi_cbg.revant.factorize.Factorize;
import de.mpi_cbg.revant.factorize.PeriodicSubstrings;
import de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3;


/**
 * Represents a set of basin descriptors as a trie with frequencies.
 * The program prints to STDOUT in DOT format. A good layout is produced by: 
 * circo -O -Teps output.dot
 */
public class PrintBasinDescriptorTree {
	
	static NumberFormat formatter = NumberFormat.getInstance();
	
	/**
	 * @param args
	 * 2: 0=all basin descriptors in a directory; frequency based on basin descriptor 
	 * intervals. The frequency of an internal node may not be correct, since it is just 
	 * the sum of the frequencies of all its leaves (plus the number of intervals that 
	 * could not be assigned to any node below it), but the same interval might have been 
	 * assigned to multiple leaves.
	 * 1=only basin descriptors selected by $IntervalGraphStepStep5$; frequency based on 
	 * number of fragments.
	 */
	public static void main(String[] args) throws IOException {
		final String STEP1_DIR = args[0];
		final int MIN_ALIGNMENT_LENGTH = Integer.parseInt(args[1]);
		final boolean MODE = Integer.parseInt(args[2])==1;
		
		final String INPUT_DIR = STEP1_DIR+"/"+IO.TAGS_DIR+"/"+IO.STEP4_DIR;
		final String DESCRIPTORS_FILE = INPUT_DIR+"/"+IO.STEP5_DIR+"/"+IO.STEP5_DESCRIPTORS_FILE;
		final String FRAGMENTS_DIR = INPUT_DIR+"/"+IO.STEP5_DIR+"/"+IO.FRAGMENTS_LABEL; 
		final int BASIN_LABEL_LENGTH = IO.BASIN_PREFIX.length();
		final String SUFFIX = ".txt";
		final int SUFFIX_LENGTH = SUFFIX.length();
		final int TAGS_PREFIX_LENGTH = IO.TAGS_PREFIX.length();
		
		boolean isCyclic;
		int i, j, p, q, r, s;
		int nLeaves, lastLeaf, count, type, frequency, period;
		String str, str2;
		Node currentNode;
		Node root = new Node(0);
		BufferedReader br, br2;
		File directory;
		int[] tmpArray = new int[20];
		String[] files, tokens;
		Node[] leaves;
		formatter.setMaximumFractionDigits(2);
		Alignments.minAlignmentLength=MIN_ALIGNMENT_LENGTH;
		
		if (!MODE) {
			// Adding tags from $basin-*-*-*.txt$ files.
			directory = new File(INPUT_DIR);
			files=directory.list();
			directory=null;
			for (i=0; i<files.length; i++) {
				if (files[i]==null || files[i].length()<BASIN_LABEL_LENGTH || !files[i].substring(0,BASIN_LABEL_LENGTH).equals(IO.BASIN_PREFIX) || !files[i].substring(files[i].length()-SUFFIX_LENGTH).equals(SUFFIX)) continue;
				br = new BufferedReader(new FileReader(INPUT_DIR+"/"+files[i]));
				str=br.readLine();
				if (str==null) {
					System.err.println("ERROR: the following file is empty: "+files[i]);
					continue;
				}
				IntervalGraphStep3.readBasinDescriptorHeader(str,tmpArray);
				isCyclic=tmpArray[1]==0; type=tmpArray[3]; frequency=tmpArray[5]; period=tmpArray[4];
				count=0; str=br.readLine();
				while (str!=null && str.length()>0) {
					count++;
					str=br.readLine();
				}
				br.close();
				if (IO.CONSISTENCY_CHECKS && count!=frequency) {
					System.err.println("ERROR: inconsistent frequency in the following file ("+frequency+"!="+count+"): "+files[i]);
					continue;
				}
				// Adding to the tree
				str=files[i].substring(0,files[i].indexOf("."));
				tokens=str.split("-");
				currentNode=root;
				for (j=1; j<tokens.length; j++) currentNode=currentNode.addChild(Integer.parseInt(tokens[j]));
				currentNode.type=isCyclic?-1:type;
				currentNode.count=count;
				currentNode.period=period;
			}
			// Adding tags from $tags-root.txt$, $tags-*.txt$ and $tags-*-*.txt$ files.
			for (i=0; i<files.length; i++) {
				if (files[i]==null || files[i].length()==0 || !files[i].substring(0,TAGS_PREFIX_LENGTH).equals(IO.TAGS_PREFIX) || !files[i].substring(files[i].length()-SUFFIX_LENGTH).equals(SUFFIX)) continue;
				p=files[i].indexOf("-"); q=files[i].indexOf("-",p+1);
				if (q>=0 && files[i].indexOf("-",q+1)>=0) continue;
				br = new BufferedReader(new FileReader(INPUT_DIR+"/"+files[i]));
				str=br.readLine();
				if (str==null) {
					System.err.println("ERROR: the following file is empty: "+files[i]);
					continue;
				}
				// Computing frequency
				count=0;
				while (str!=null) {
					count++;
					str=br.readLine();
				}
				br.close();
				// Adding frequency to the tree
				str=files[i].substring(0,files[i].indexOf("."));
				tokens=str.split("-");
				currentNode=root;
				if (!tokens[1].equalsIgnoreCase("root") && (tokens.length<3 || !tokens[2].equalsIgnoreCase("step4"))) {
					for (j=1; j<tokens.length; j++) currentNode=currentNode.addChild(Integer.parseInt(tokens[j]));
				}
				currentNode.noKernelCount+=count;
			}
			files=null;
			// Highlighting basins selected in Step 5
			br = new BufferedReader(new FileReader(DESCRIPTORS_FILE));
			str=br.readLine();
			while (str!=null) {
				p=str.indexOf("-"); q=str.indexOf("-",p+1); r=str.indexOf("-",q+1); s=str.indexOf(".",r+1);
				currentNode=root;
				currentNode=currentNode.addChild(Integer.parseInt(str.substring(p+1,q)));
				currentNode=currentNode.addChild(Integer.parseInt(str.substring(q+1,r)));
				currentNode=currentNode.addChild(Integer.parseInt(str.substring(r+1,s)));
				currentNode.selected=true;
				str=br.readLine();
			}
			br.close();
		}
		else {
			br = new BufferedReader(new FileReader(DESCRIPTORS_FILE));
			str=br.readLine();
			while (str!=null) {
				br2 = new BufferedReader(new FileReader(INPUT_DIR+"/"+str.trim()));
				str2=br2.readLine();
				br2.close();
				if (str2==null) {
					System.err.println("ERROR: the following file is empty: "+str);
					continue;
				}
				IntervalGraphStep3.readBasinDescriptorHeader(str2,tmpArray);
				isCyclic=tmpArray[1]==0; type=tmpArray[3]; frequency=tmpArray[5]; period=tmpArray[4];
				str2=str.substring(IO.BASIN_PREFIX_LENGTH,str.indexOf("."));
				tokens=str2.split("-");
				// Computing frequency
				try { br2 = new BufferedReader(new FileReader(FRAGMENTS_DIR+"/fragments-"+str2.trim()+".txt")); }
				catch (IOException e) {
					// The basin might not have any fragment
					str=br.readLine();
					continue;
				}
				count=0; str2=br2.readLine();
				while (str2!=null && str2.length()>0) {
					count++;
					str2=br2.readLine();
				}
				br2.close();
				// Adding to the tree
				currentNode=root;
				for (j=1; j<tokens.length; j++) currentNode=currentNode.addChild(Integer.parseInt(tokens[j]));
				currentNode.type=isCyclic?-1:type;
				currentNode.count=count;
				currentNode.period=period;
				str=br.readLine();
			}
			br.close();
		}
		
		// Printing the tree
		nLeaves=root.nLeaves();
		System.err.println("There are "+nLeaves+" distinct kernel tags");
		root.getFrequency(true); root.getFrequency(false); 
		root.removeUnaryPaths();
		System.out.println("digraph G {\n");
		System.out.println("bgcolor=\""+Colors.COLOR_DOT_BACKGROUND+"\";\n");
		Node.dotID=0;
		root.toDot("");
		System.out.println("}");
		leaves = new Node[nLeaves];
		lastLeaf=root.collectLeaves(leaves,0);
		if (lastLeaf!=nLeaves-1) {
			System.err.println("ERROR in collectLeaves()");
			System.exit(1);
		}
		Node.order=Node.ORDER_COUNT;
		Arrays.sort(leaves,0,nLeaves);
		System.err.println("Leaves in order of frequency:");
		for (i=0; i<nLeaves; i++) System.err.println(leaves[i].count+","+leaves[i].dotPrefix);
	}
	
	
	private static class Node implements Comparable {
		public static final int ORDER_ID = 0;
		public static final int ORDER_COUNT = 1;
		public static int order;
		public boolean selected;
		public int type;  // -1=cyclic; -2=unknown.
		public int id;
		public int count, noKernelCount, period;
		public Node[] children;
		public int lastChild;
		public String dotPrefix;
		
		// Temporary space
		private static Node tmpNode = new Node(-1);
		private static int dotID;
		
		
		public Node(int id) {
			this.id=id;
			count=0; noKernelCount=0;
			children=null;
			lastChild=-1;
			type=-2;
			period=0;
			selected=false;
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
		
		
		public Node addChild(int id) {
			final int CAPACITY = 10;
			int i;
			final int previousOrder = Node.order;
			Node newChild;
			
			// Checking if the child exists
			tmpNode.id=id;
			if (lastChild==0) {
				if (children[0].id==id) return children[0];
			}
			else if (lastChild>0) {
				Node.order=Node.ORDER_ID;
				i=Arrays.binarySearch(children,0,lastChild+1,tmpNode);
				Node.order=previousOrder;
				if (i>=0) return children[i];
			}
			
			// Adding a new child
			lastChild++;
			if (children==null || children.length==lastChild) {
				Node[] newChildren = new Node[(children==null||children.length==0)?CAPACITY:children.length<<1];
				if (children!=null && children.length>0) System.arraycopy(children,0,newChildren,0,children.length);
				children=newChildren;
			}
			newChild = new Node(id);
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
		 * Stores in $out[from..X]$ the leaves in the subtree of this node, where $X$ is
		 * returned in output.
		 */
		public int collectLeaves(Node[] out, int from) {
			int i, last;
			
			if (lastChild==-1) {
				out[from]=this;
				return from;
			}
			last=from-1;
			for (i=0; i<=lastChild; i++) {
				last=children[i].collectLeaves(out,from);
				from=last+1;
			}
			return last;
		}
		
		
		/**
		 * @param which TRUE=count, FALSE=noKernelCount;
		 * @return the sum of the count values of all leaves in the subtree. 
		 * The procedure sets to this value the count field of every internal node.
		 */
		public int getFrequency(boolean which) {
			int i;
			
			if (which) {
				if (lastChild==-1) return count;
				for (i=0; i<=lastChild; i++) count+=children[i].getFrequency(true);
				return count;
			}
			else {
				if (lastChild==-1) return noKernelCount;
				for (i=0; i<=lastChild; i++) noKernelCount+=children[i].getFrequency(false);
				return noKernelCount;
			}
		}
		
		
		public void removeUnaryPaths() {
			Node currentNode = this;
			while (currentNode.lastChild==0) currentNode=currentNode.children[0];
			type=currentNode.type;
			count=currentNode.count;
			noKernelCount=currentNode.noKernelCount;
			period=currentNode.period;
			children=currentNode.children;
			lastChild=currentNode.lastChild;
			for (int i=0; i<=lastChild; i++) children[i].removeUnaryPaths();
		}

		
		/**
		 * White = cyclic.
		 * Green border = selected by Step 5.
		 */
		public void toDot(String dotPrefix) {
			final double PEN_SCALE = 5.0;
			final double NODE_SCALE = 1.0;
			final double FONT_SCALE = 10.0;
			final int oldID = dotID;
			final int nodeSize = Math.max(1,(int)(Math.sqrt(count)*NODE_SCALE));
			final int fontSize = Math.max(1,(int)(Math.sqrt(count)*FONT_SCALE));
			int edgeSize, penwidth;
			String fillColor, borderColor;
			
			switch (type) {
				case -1: fillColor="#FFFFFF"; break;  // Cyclic
				case -2: fillColor=Colors.COLOR_DOT_EDGE; break;  // Unknown
				case Constants.INTERVAL_PERIODIC: fillColor=period>=PeriodicSubstrings.MIN_PERIOD_LONG?Colors.COLOR_DOT_LONGPERIOD:Colors.type2color_dot(Constants.INTERVAL_PERIODIC); break;
				default: fillColor=Colors.type2color_dot(type); break;
			}
			// Remark: the period is not printed if it is zero, since it means that the
			// previous steps of the pipeline could not estimate it.
			penwidth=Math.max(1,(int)(Math.sqrt(count)*PEN_SCALE));
			System.out.println(dotID+" [shape=\"circle\",style=\"filled\",fixedsize=\"true\",width="+nodeSize+",fontsize=\""+fontSize+"\",label=\""+dotPrefix+" \n f="+formatter.format(((double)count)/1000)+"k"+(noKernelCount==0?"":"\n("+formatter.format(((double)noKernelCount)/1000)+"k)")+(type==Constants.INTERVAL_PERIODIC&&period>0?"\n p="+period:"")+"\",fillcolor=\""+fillColor+"\""+(selected?",penwidth=\""+penwidth+"\",color=\""+Colors.COLOR_DOT_REFERENCE_ALIEN_HIGHLIGHT+"\"":"")+",fontcolor=\""+Colors.COLOR_DOT_TEXT+"\",fontname=\"Arial\"];\n");
			this.dotPrefix=dotPrefix;
			for (int i=0; i<=lastChild; i++) {
				dotID++;
				edgeSize=Math.max(1,(int)(Math.sqrt(children[i].count)*PEN_SCALE));
				System.out.println(oldID+" -> "+dotID+" [dir=\"both\",arrowhead=\"none\",arrowtail=\"none\",penwidth=\""+edgeSize+"\",color=\""+Colors.COLOR_DOT_EDGE+"\"];\n");
				children[i].toDot((dotPrefix.length()==0?"":dotPrefix+".")+children[i].id);
			}
		}
	}
	
}