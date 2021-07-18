package de.mpi_cbg.revant.biology;

import java.util.Arrays;
import java.util.HashMap;
import java.io.*;
import java.util.zip.Deflater;
import de.mpi_cbg.revant.util.IO;
import de.mpi_cbg.revant.util.Constants;
import de.mpi_cbg.revant.util.Math;


/**
 * Basic functions to deal with a RepBase dataset
 */
public class RepBase {

	private static final String TRANSPOSABLE_ELEMENT_LABEL = "transposable-element";
	private static final String OTHER_LABEL = "other";
	private static final String REST_LABEL = "rest";
	private static final int MAX_STRINGS = 100000;
	
	/**
	 * Main data structures
	 */
	public static RepbaseTreeNode root;
	public static int nNodes;  // Number of nodes under $root$, included.
	public static RepbaseTreeNode[] string2node;
	
	/**
	 * Temporary space
	 */
	private static double[] tmp;
	
	
	/**
	 * Remark: elements in $FASTA_FILE$ that are not present in $treeFile$ are added as
	 * new children of the root of $TREE_FILE$.
	 *
	 * @return the maximum depth of a node.
	 */
	public static final int loadRepbase(String treeFile, String fastaFile, boolean compress) throws IOException {
		int nStrings, maxDepth;
		HashMap<String,RepbaseTreeNode> map;
		
		// Building
		map = new HashMap<String,RepbaseTreeNode>();
		root=loadRepbaseTree(treeFile,map);
		System.err.println("RepBase> Tree loaded");
		nStrings=loadRepbaseFasta(fastaFile,map,compress);
		string2node=loadString2node(nStrings,fastaFile,map);
		System.err.println("RepBase> Strings loaded");
		map.clear(); map=null;
		root.depth=0; maxDepth=root.setDepth();
		nNodes=root.setID(0)+1;
		
		// Computing statistics
		tmp = new double[MAX_STRINGS];
		computeStats(root,true);
		computeStats(root,false);
		System.err.println("RepBase> Statistics computed");
		return maxDepth;
	}
	
	
	/**
	 * Builds a representation of the RepBase tree using the HTML of the multiple choice
	 * box at <http://www.girinst.org/repbase/update/browse.php>. I couldn't find an
	 * explicit representation of the hierarchy in the RepBase download files.
	 *
	 * @return the root of the tree.
	 */
	public static final RepbaseTreeNode loadRepbaseTree(String path, HashMap<String,RepbaseTreeNode> map) throws IOException {
		int i, j, k, h;
		int depth;
		String str, label;
		BufferedReader br;
		RepbaseTreeNode root, currentNode, tmpNode;
		
		// Loading the tree from $path$
		root = new RepbaseTreeNode("#",-1);
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		str=br.readLine();
		currentNode=root;
		while (str!=null) {
			i=str.indexOf(">")+1; j=str.indexOf("<",i+1)-1;
			depth=0; h=i-1;
			for (k=i; k<=j; k++) {
				if (str.charAt(k)==';') {
					depth++;
					h=k;
				}
			}
			label=str.substring(h+1,j+1).trim().toLowerCase().replaceAll(" ","-");
			tmpNode = new RepbaseTreeNode(label,depth);
			map.put(label,tmpNode);
			while (currentNode.depth>=depth) currentNode=currentNode.parent;
			currentNode.addChild(tmpNode);
			currentNode=tmpNode;
			str=br.readLine();
		}
		br.close();
		
		// Adding artificial node $OTHER_LABEL$
		tmpNode = new RepbaseTreeNode(OTHER_LABEL,1);
		root.addChild(tmpNode);
		map.put(OTHER_LABEL,tmpNode);
		return root;
	}
	
	
	/**
	 * Loads in the tree a set of properties computed from the RepBase FASTA file $path$.
	 *
	 * @return the number of strings in $path$.
	 */
	public static final int loadRepbaseFasta(String path, HashMap<String,RepbaseTreeNode> map, boolean compress) throws IOException {
		int i, j;
		int denominator, nBytes, nStrings;
		double numerator;
		byte[] uncompressed, compressed;
		String str;
		StringBuilder buffer;
		BufferedReader br;
		RepbaseTreeNode node, tmpNode;
		Deflater compressor;
		
		buffer = new StringBuilder();
		node=null;
		nStrings=0;
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		str=br.readLine();
		while (str!=null) {
			if (str.length()>0 && str.charAt(0)=='>') {
				nStrings++;
				if (buffer.length()>0) {
					node.addLength(buffer.length());
					numerator=buffer.length();
					if (compress) {
						uncompressed=buffer.toString().getBytes();
						compressed = new byte[uncompressed.length];
						compressor = new Deflater(Deflater.BEST_COMPRESSION);
						compressor.setInput(uncompressed);
						compressor.finish();
						denominator=compressor.deflate(compressed);
						compressor.end();
					}
					else denominator=1;
					node.addCompressibility(numerator/denominator);
					buffer.delete(0,buffer.length());
				}
				i=str.indexOf("\t");
				j=str.lastIndexOf("\t");
				if (i==-1||j==-1) {
					// Anomalous header: trying to match a prefix to the tree.
					tmpNode=root.search(str.trim().toLowerCase().replaceAll(" ","-"));
					if (tmpNode!=null) node=tmpNode;
					else {
						node = new RepbaseTreeNode(str,1);
						root.addChild(node);
						map.put(str,node);						
						System.err.println("UNKNOWN TYPE ADDED TO THE TREE: "+str);
						//node=map.get(OTHER_LABEL);
					}
				}
				else {
					str=str.substring(i+1,j).trim().toLowerCase().replaceAll(" ","-");
					node=map.get(str);
					if (node==null) {  // Correct header but element not in the tree
						node = new RepbaseTreeNode(str,1);
						root.addChild(node);
						map.put(str,node);
					}
				}
			}
			else buffer.append(str.trim());
			str=br.readLine();
		}
		br.close();
		node.addLength(buffer.length());
		numerator=buffer.length();
		if (compress) {
			uncompressed=buffer.toString().getBytes();
			compressed = new byte[uncompressed.length];
			compressor = new Deflater(Deflater.BEST_COMPRESSION);
			compressor.setInput(uncompressed);
			compressor.finish();
			denominator=compressor.deflate(compressed);
			compressor.end();
		}
		else denominator=1;
		node.addCompressibility(numerator/denominator);
		return nStrings;
	}
	
	
	/**
	 * Returns a map from each string in file $path$ to the corresponding leaf in the 
	 * RepBase tree.
	 */
	private static final RepbaseTreeNode[] loadString2node(int nStrings, String path, HashMap<String,RepbaseTreeNode> map) throws IOException {
		int i, j, r;
		String str;
		StringBuilder buffer;
		BufferedReader br;
		RepbaseTreeNode node, tmpNode;
		RepbaseTreeNode[] out;
		
		out = new RepbaseTreeNode[nStrings];		
		buffer = new StringBuilder();
		node=null;
		br = new BufferedReader(new FileReader(path),IO.BUFFER_SIZE);
		str=br.readLine(); r=-1;
		while (str!=null) {
			if (str.length()>0 && str.charAt(0)=='>') {
				i=str.indexOf("\t");
				j=str.lastIndexOf("\t");
				if (i==-1||j==-1) {
					// Anomalous header: trying to match a prefix to the tree.
					tmpNode=root.search(str.trim().toLowerCase().replaceAll(" ","-"));
					if (tmpNode!=null) node=tmpNode;
					else node=map.get(str);
				}
				else {
					str=str.substring(i+1,j).trim().toLowerCase().replaceAll(" ","-");
					node=map.get(str);
				}
				out[++r]=node;
			}
			str=br.readLine();
		}
		br.close();
		return out;
	}
	
	
	/**
	 * Computes statistics for all nodes in the tree.
	 *
	 * @param type TRUE=lengths, FALSE=compressibilities.
	 */
	private static final double[] computeStats(RepbaseTreeNode node, boolean type) {
		int i, quarter, last;
		double[] from1, from2;
		
		// Base case
		if (node.lastChild==-1) return node.computeStats(type);
		
		// Recursion
		last=-1;
		from1=computeStats(node.children[0],type);
		for (i=1; i<=node.lastChild; i++) {
			from2=computeStats(node.children[i],type);
			last=Math.setUnion(from1,from1.length-1,from2,from2.length-1,tmp);
			from1=Arrays.copyOf(tmp,last+1);
		}
		if (type) node.lastLength=last==-1?0:last+1;
		else node.lastCompressibility=last==-1?0:last+1;
		if (last==-1) {
			if (type) {
				node.lengthMin=-1;
				node.lengthMax=-1;
				node.length1q=-1;
				node.length2q=-1;
				node.length3q=-1;
			}
			else {
				node.compMin=-1;
				node.compMax=-1;
				node.comp1q=-1;
				node.comp2q=-1;
				node.comp3q=-1;
			}
			return new double[0];
		}
		else if (last==0) {
			if (type) {
				node.lengthMin=tmp[0];
				node.lengthMax=tmp[0];
				node.length1q=tmp[0];
				node.length2q=tmp[0];
				node.length3q=tmp[0];
			}
			else {
				node.compMin=tmp[0];
				node.compMax=tmp[0];
				node.comp1q=tmp[0];
				node.comp2q=tmp[0];
				node.comp3q=tmp[0];
			}
			return Arrays.copyOf(tmp,last+1);
		}
		else if (last==1) {
			if (type) {
				node.lengthMin=tmp[0];
				node.lengthMax=tmp[1];
				node.length1q=tmp[0];
				node.length2q=tmp[0];
				node.length3q=tmp[1];
			}
			else {
				node.compMin=tmp[0];
				node.compMax=tmp[1];
				node.comp1q=tmp[0];
				node.comp2q=tmp[0];
				node.comp3q=tmp[1];
			}
			return Arrays.copyOf(tmp,last+1);
		}
		else if (last==2) {
			if (type) {
				node.lengthMin=tmp[0];
				node.lengthMax=tmp[2];
				node.length1q=tmp[0];
				node.length2q=tmp[1];
				node.length3q=tmp[2];
			}
			else {
				node.compMin=tmp[0];
				node.compMax=tmp[2];
				node.comp1q=tmp[0];
				node.comp2q=tmp[1];
				node.comp3q=tmp[2];
			}
			return Arrays.copyOf(tmp,last+1);
		}
		quarter=(last+1)>>2;
		if (type) {
			node.lengthMin=tmp[0];
			node.lengthMax=tmp[last];
			node.length1q=tmp[quarter-1];
			node.length2q=tmp[(quarter<<1)-1];
			node.length3q=tmp[3*quarter-1];
		}
		else {
			node.compMin=tmp[0];
			node.compMax=tmp[last];
			node.comp1q=tmp[quarter-1];
			node.comp2q=tmp[(quarter<<1)-1];
			node.comp3q=tmp[3*quarter-1];
		}
		return Arrays.copyOf(tmp,last+1);
	}
	
	
	public static class RepbaseTreeNode {
		final static int MAX_CHILDREN = 1000;
		final static int MAX_COLOR = 255;
		public RepbaseTreeNode parent;
		public RepbaseTreeNode[] children;
		public int lastChild;
		public int depth;
		public int count;
		public int id;
		
		/**
		 * Biological properties
		 */	
		public String label;
		public double[] lengths, compressibilities;
		public int lastLength, lastCompressibility;
		
		/**
		 * Aggregate statistics
		 */
		public double lengthMin, length1q, length2q, length3q, lengthMax;
		public double compMin, comp1q, comp2q, comp3q, compMax;
		
		/**
		 * Temporary space
		 */
		public double[] typeCounts;
		private boolean marked;
		
		
		public RepbaseTreeNode(String label, int depth) {
			children = new RepbaseTreeNode[MAX_CHILDREN];
			lastChild=-1;
			this.label=label;
			this.depth=depth;
			
			// Biological properties
			lengths = new double[MAX_STRINGS];
			lastLength=-1;
			compressibilities = new double[MAX_STRINGS];
			lastCompressibility=-1;
			
			// Temporary space
			typeCounts = new double[Constants.INTERVAL_PERIODIC+1];
		}
		
		
		public final void addChild(RepbaseTreeNode node) {
			children[++lastChild]=node;
			node.parent=this;
		}
		
		
		public final void addLength(int length) {
			lengths[++lastLength]=length;
		}
		
		
		public final void addCompressibility(double compressibility) {
			compressibilities[++lastCompressibility]=compressibility;
		}
		
		
		/**
		 * @param type TRUE=lengths, FALSE=compressibilities.
		 */
		public final double[] computeStats(boolean type) {
			int quarter;
			
			if (type) {
				if (lastLength==-1) {
					lengthMin=-1;
					lengthMax=-1;
					length1q=-1;
					length2q=-1;
					length3q=-1;
					return new double[0];
				}
				else if (lastLength==0) {
					lengthMin=lengths[0];
					lengthMax=lengths[0];
					length1q=lengths[0];
					length2q=lengths[0];
					length3q=lengths[0];
					return Arrays.copyOf(lengths,lastLength+1);
				}
				else if (lastLength==1) {
					Arrays.sort(lengths,0,lastLength+1);
					lengthMin=lengths[0];
					lengthMax=lengths[1];
					length1q=lengths[0];
					length2q=lengths[0];
					length3q=lengths[1];
					return Arrays.copyOf(lengths,lastLength+1);
				}
				else if (lastLength==2) {
					Arrays.sort(lengths,0,lastLength+1);
					lengthMin=lengths[0];
					lengthMax=lengths[2];
					length1q=lengths[0];
					length2q=lengths[1];
					length3q=lengths[2];
					return Arrays.copyOf(lengths,lastLength+1);
				}
				Arrays.sort(lengths,0,lastLength+1);
				lengthMin=lengths[0];
				lengthMax=lengths[lastLength];
				quarter=(lastLength+1)>>2;
				length1q=lengths[quarter-1];
				length2q=lengths[(quarter<<1)-1];
				length3q=lengths[3*quarter-1];
				if (lastLength==lengths.length-1) return lengths;
				return Arrays.copyOf(lengths,lastLength+1);
			}
			else {
				if (lastCompressibility==-1) {
					compMin=-1;
					compMax=-1;
					comp1q=-1;
					comp2q=-1;
					comp3q=-1;
					return new double[0];
				}
				else if (lastCompressibility==0) {
					compMin=compressibilities[0];
					compMax=compressibilities[0];
					comp1q=compressibilities[0];
					comp2q=compressibilities[0];
					comp3q=compressibilities[0];
					return Arrays.copyOf(compressibilities,lastCompressibility+1);
				}
				else if (lastCompressibility==1) {
					Arrays.sort(lengths,0,lastLength+1);
					compMin=compressibilities[0];
					compMax=compressibilities[1];
					comp1q=compressibilities[0];
					comp2q=compressibilities[0];
					comp3q=compressibilities[1];
					return Arrays.copyOf(compressibilities,lastCompressibility+1);
				}
				else if (lastCompressibility==2) {
					Arrays.sort(lengths,0,lastLength+1);
					compMin=compressibilities[0];
					compMax=compressibilities[2];
					comp1q=compressibilities[0];
					comp2q=compressibilities[1];
					comp3q=compressibilities[2];
					return Arrays.copyOf(compressibilities,lastCompressibility+1);
				}
				Arrays.sort(compressibilities,0,lastCompressibility+1);
				compMin=compressibilities[0];
				compMax=compressibilities[lastCompressibility];
				quarter=(lastCompressibility+1)>>2;
				comp1q=compressibilities[quarter-1];
				comp2q=compressibilities[(quarter<<1)-1];
				comp3q=compressibilities[3*quarter-1];
				if (lastCompressibility==compressibilities.length-1) return compressibilities;
				return Arrays.copyOf(compressibilities,lastCompressibility+1);
			}
		}
		
		
		/**
		 * @return the maximum value of $depth$ of a node in the subtree rooted at this 
		 * node.
		 */
		private final int setDepth() {
			int i, d, max;
			
			max=depth;
			for (i=0; i<=lastChild; i++) {
				children[i].depth=depth+1;
				d=children[i].setDepth();
				if (d>max) max=d;
			}
			return max;
		}
		
		
		/**
		 * Assigns IDs to all nodes in the subtree rooted at this node (included), 
		 * starting from $id$, which is assigned to this node.
		 *
		 * @return the last ID assigned.
		 */
		public final int setID(int id) {
			int i, d;
			
			this.id=id;
			d=id;
			for (i=0; i<=lastChild; i++) d=children[i].setID(d+1);
			return d;
		}
		
		
		/**
		 * @return the label of a node that matches a prefix of $str$.
		 */
		public final RepbaseTreeNode search(String str) {
			int i, length;
			RepbaseTreeNode tmpNode;
			
			for (i=0; i<=lastChild; i++) {
				tmpNode=children[i].search(str);
				if (tmpNode!=null) return tmpNode;
			}
			length=label.length();
			if (length<str.length() && str.indexOf(label)>=0) return this;
			return null;
		}
		
		
		public final void initTypeCounts() {
			int i;
			
			Math.set(typeCounts,typeCounts.length-1,0.0);
			for (i=0; i<=lastChild; i++) children[i].initTypeCounts();
		}
		
		
		/**
		 * Remark: the procedure also sums prefix and suffix types for each node, storing 
		 * the sum at the prefix index of $typeCounts$.
		 */
		public final void getMaxCounts(double[] out) {
			int i;
			
			typeCounts[Constants.INTERVAL_DENSE_PREFIX]+=typeCounts[Constants.INTERVAL_DENSE_SUFFIX];
			typeCounts[Constants.INTERVAL_DENSE_SUFFIX]=0.0;
			for (i=0; i<typeCounts.length; i++) {
				if (i==Constants.INTERVAL_DENSE_SUFFIX) continue;
				if (typeCounts[i]>out[i]) out[i]=typeCounts[i];
			}
			for (i=0; i<=lastChild; i++) children[i].getMaxCounts(out);
		}
		
		
		/**
		 * Remark: does not consider the position of $typeCounts$ that corresponds to
		 * suffix substrings.
		 */
		public final void divideCounts(double[] denominator) {
			int i;
			
			for (i=0; i<typeCounts.length; i++) {
				if (i==Constants.INTERVAL_DENSE_SUFFIX) continue;
				typeCounts[i]/=denominator[i];
			}
			for (i=0; i<=lastChild; i++) children[i].divideCounts(denominator);
		}
		
		
		/**
		 * Appends the descriptor of each node to the buffer that corresponds to its 
		 * $depth$ value.
		 */
		public final void toDot(StringBuilder[] levels) {
			levels[depth].append(id+" [label=\""+label+"\"];\n");
			for (int i=0; i<=lastChild; i++) children[i].toDot(levels);
		}
		
		
		/**
		 * Prints the topology of the tree to the standard output
		 */
		public final void toDot() {
			for (int i=0; i<=lastChild; i++) {
				System.out.println(id+" -> "+children[i].id+";");
				children[i].toDot();
			}
		}
		
		
		/**
		 * Prints a heatmap on the topology of the tree, induced by normalized counts of
		 * type $type$. Only the subgraph of the tree induced by marked nodes is printed.
		 */
		public final void toDot(int type, BufferedWriter bw) throws IOException {
			if (!marked) return;
			bw.write(id+" [label=\""+label+"\",color=\"#"+Integer.toHexString((int)(MAX_COLOR*typeCounts[type]))+"0000\",style=\"filled\",fontcolor=\"white\"];\n");
			for (int i=0; i<=lastChild; i++) {
				if (!children[i].marked) continue;
				bw.write(id+" -> "+children[i].id+";\n");
				children[i].toDot(type,bw);
			}
		}
		
		
		public final void mark(int type) {
			marked=false;
			for (int i=0; i<=lastChild; i++) {
				children[i].mark(type);
				marked|=children[i].marked;
			}
			if (typeCounts[type]>0) marked=true;
		}
		
	}
	
	
	/**
	 * Expands just the subtree of transposable elements
	 */
	public static final void printTree(RepbaseTreeNode root, boolean printLabel) throws IOException {
		int i, j;
		String label;
		RepbaseTreeNode node;
		BufferedWriter bw, restBw;
		
		restBw = new BufferedWriter(new FileWriter(REST_LABEL+".txt"),IO.BUFFER_SIZE);
		for (i=0; i<=root.lastChild; i++) {
			node=root.children[i];
			label=node.label;
			if (label.equals(TRANSPOSABLE_ELEMENT_LABEL)) {
				for (j=0; j<=node.lastChild; j++) {
					bw = new BufferedWriter(new FileWriter(node.children[j].label+".txt"),IO.BUFFER_SIZE);
					printLeaves(node.children[j],printLabel,bw);
					bw.close();
				}
			}
			else if (node.lastChild>=0 || node.label==OTHER_LABEL) {
				bw = new BufferedWriter(new FileWriter(label+".txt"),IO.BUFFER_SIZE);
				printLeaves(node,printLabel,bw);
				bw.close();
			}
			else printLeaves(node,printLabel,restBw);
		}
		restBw.close();
	}
	
	
	/**
	 * Prints all leaves in the subtree rooted at $node$
	 */
	public static final void printLeaves(RepbaseTreeNode node, boolean printLabel, BufferedWriter bw) throws IOException {
		if (node.lastChild==-1 && node.length2q!=-1 && node.comp2q!=-1) {
			if (node.lastLength!=node.lastCompressibility) {
				System.err.println("Error");
				System.exit(1);
			}
			bw.write((printLabel?node.label+",":"")+(node.lastLength+1)+","+node.length2q+","+node.comp2q+"\n");
		}
		for (int i=0; i<=node.lastChild; i++) printLeaves(node.children[i],printLabel,bw);
	}
	
}