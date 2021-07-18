package de.mpi_cbg.revant.graphics;

import java.io.*;
import de.mpi_cbg.revant.util.IO;


public class PrintComponentTree {
	
	
	/**
	 * java PrintComponentTree ~/Desktop/components 1 ~/Desktop/components/componentTree.dot 2
	 */
	public static void main(String[] args) throws IOException {
		final String ROOT_FOLDER = args[0];
		final int N_STEP1_COMPONENTS = Integer.parseInt(args[1]);
		final String TREE_FILE = args[2];
		final int MIN_SIZE = Integer.parseInt(args[3]);
		int i, j, p, q;
		String str, fileName;
		ComponentTreeNode root, step1node, step2node, kernelNode;
		BufferedReader br;
		BufferedWriter bw;
		File[] files;
		
		// Loading component tree
		root = new ComponentTreeNode(0);
		for (i=0; i<N_STEP1_COMPONENTS; i++) {
			step1node=root.addChild(i);
			files = new File(ROOT_FOLDER+"/"+i+"-step2").listFiles();
			for (j=0; j<files.length; j++) {
				fileName=files[j].getName();
				q=fileName.indexOf(".repeat");
				if (q<0) continue;
				p=fileName.indexOf("-");
				step2node=step1node.addChild(Integer.parseInt(fileName.substring(0,p)));
				kernelNode=step2node.addChild(Integer.parseInt(fileName.substring(p+1,q)));
				br = new BufferedReader(new FileReader(files[j].getAbsolutePath()),IO.BUFFER_SIZE);
				str=br.readLine();
				kernelNode.size=0;
				while (str!=null) {
					kernelNode.size++;
					str=br.readLine();
				}				
				br.close();
				//System.err.println("loaded file <"+files[j].getAbsolutePath()+">, size="+kernelNode.size);
			}
		}
		
		// Printing to DOT
		root.getSize(); root.getDepth(0);
		System.err.println("maxSize="+root.size);
		bw = new BufferedWriter(new FileWriter(TREE_FILE),IO.BUFFER_SIZE);
		bw.write("digraph G {\n"); 
		root.print(bw,root.size,MIN_SIZE); 
		bw.write("}\n");
		bw.close();
	}


	private static class ComponentTreeNode {
		private static final int MIN_CHILDREN = 10;
		
		public int id, lastChild;
		public int size, depth;
		public ComponentTreeNode parent;
		public ComponentTreeNode[] children;
		
		
		public ComponentTreeNode(int id) {
			this.id=id;
			children = new ComponentTreeNode[MIN_CHILDREN];
			lastChild=-1;
			size=0;
			parent=null;
		}
		
		
		public final ComponentTreeNode addChild(int childID) {
			int i;
			ComponentTreeNode tmpNode;
			ComponentTreeNode[] tmp;
			
			for (i=0; i<=lastChild; i++) {
				if (children[i].id==childID) return children[i];
			}
			tmpNode = new ComponentTreeNode(childID);
			tmpNode.parent=this;
			if (lastChild+1>children.length-1) {
				tmp = new ComponentTreeNode[children.length<<1];
				System.arraycopy(children,0,tmp,0,lastChild+1);
				children=tmp;
			}
			children[++lastChild]=tmpNode;
			return tmpNode;
		}
		
		
		/**
		 * Sets $size$ to the sum of $size$ values at the leaves.
		 */
		public final void getSize() {
			int i;
			
			if (lastChild==-1) return;
			size=0;
			for (i=0; i<=lastChild; i++) {
				children[i].getSize();
				size+=children[i].size;
			}
		}
		
		
		/**
		 * Assigns a tree depth to each node
		 */
		public final void getDepth(int d) {
			depth=d;
			for (int i=0; i<=lastChild; i++) children[i].getDepth(d+1);
		}
		
		
		public String toString() {
			String out;
			ComponentTreeNode node;
			
			out="n"; node=this;
			while (node!=null) {
				out=node.id+"_"+out;
				node=node.parent;
			}
			out="n"+out;
			return out;
		}
		
		
		private boolean shouldBePrinted(int minSize) {
			if (size<minSize) return false;
			if ((depth==1 || depth==2) && lastChild==0) return false;
			return true;
		}

		
		public void print(BufferedWriter bw, int maxSize, int minSize) throws IOException {
			final double MAX_WIDTH = 1000;
			int i;
			double nodeSize;
			String stringID, sizeStr;
			String[] tokens;
			
			// Printing node
			if (!shouldBePrinted(minSize)) return;
			stringID=toString();
			nodeSize=(MAX_WIDTH*size)/maxSize;
			sizeStr=IO.format(nodeSize);
			if (sizeStr.indexOf(",")>=0) {
				tokens=sizeStr.split(",");
				sizeStr="";
				for (i=0; i<tokens.length; i++) sizeStr+=tokens[i];
			}
			bw.write(stringID+" [shape=ellipse,width="+sizeStr+",height="+sizeStr+"];\n");
			
			// Printing edges
			for (i=0; i<=lastChild; i++) {
				if (children[i].shouldBePrinted(minSize)) bw.write(stringID+" -> "+children[i].toString()+";\n");
			}
			for (i=0; i<=lastChild; i++) {
				if (children[i].shouldBePrinted(minSize)) children[i].print(bw,maxSize,minSize);
			}
		}
	}
	

}