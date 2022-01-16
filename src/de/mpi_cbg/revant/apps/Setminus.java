package de.mpi_cbg.revant.apps;

import java.io.*;

/**
 * Given two files containing non-negative sorted integers, one per line, the program 
 * prints to STDOUT the difference between the first file and the second file.
 *
 * Remark: the Unix $comm$ program expects the files to be sorted lexicographically, not
 * numerically. $grep$ using one file as a pattern set is too slow in practice.
 */
public class Setminus {
	
	public static void main(String[] args) throws IOException {
		final String FILE_1 = args[0];
		final String FILE_2 = args[1];
		int number1, number2;
		String str1, str2;
		BufferedReader br1, br2;
		
		br1 = new BufferedReader(new FileReader(FILE_1));
		br2 = new BufferedReader(new FileReader(FILE_2));
		str1=br1.readLine(); 
		if (str1==null) return;
		number1=Integer.parseInt(str1);
		str2=br2.readLine();
		number2=-1;
		if (str2!=null) number2=Integer.parseInt(str2);
		while (str1!=null && str2!=null) {
			if (number2<number1) {
				str2=br2.readLine();
				if (str2!=null) number2=Integer.parseInt(str2);
			}
			else if (number2>number1) {
				System.out.println(str1);
				str1=br1.readLine();
				if (str1!=null) number1=Integer.parseInt(str1);
			}
			else {
				str1=br1.readLine();
				if (str1!=null) number1=Integer.parseInt(str1);
				str2=br2.readLine();
				if (str2!=null) number2=Integer.parseInt(str2);
			}
		}
		while (str1!=null) {
			System.out.println(str1);
			str1=br1.readLine();
		}
		br1.close(); br2.close();
	}

}