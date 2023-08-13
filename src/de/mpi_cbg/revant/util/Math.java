package de.mpi_cbg.revant.util;

import java.util.Random;
import java.util.Arrays;


/**
 * Wraps $java.lang.Math$, adding a few custom constants and implementations.
 */
public class Math {
	/**
	 * Global symbols for infinity
	 */
	public static final int POSITIVE_INFINITY = Integer.MAX_VALUE;
	public static final int NEGATIVE_INFINITY = Integer.MIN_VALUE;

	/**
	 * Row i \in [0..99]: i+1 degrees of freedom. Columns: chi square with cumulative
	 * probability at most X. 0: X=0.9. 1: X=0.95. 2: X=0.975. 3: X=0.99. 4: X=0.999.
	 *
	 * Source: <https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test>.
	 */
	public static final double[][] CHI_SQUARE = new double[][] {
		{2.706,3.841,5.024,6.635,10.828},
		{4.605,5.991,7.378,9.21,13.816},
		{6.251,7.815,9.348,11.345,16.266},
		{7.779,9.488,11.143,13.277,18.467},
		{9.236,11.07,12.833,15.086,20.515},
		{10.645,12.592,14.449,16.812,22.458},
		{12.017,14.067,16.013,18.475,24.322},
		{13.362,15.507,17.535,20.09,26.125},
		{14.684,16.919,19.023,21.666,27.877},
		{15.987,18.307,20.483,23.209,29.588},
		{17.275,19.675,21.92,24.725,31.264},
		{18.549,21.026,23.337,26.217,32.91},
		{19.812,22.362,24.736,27.688,34.528},
		{21.064,23.685,26.119,29.141,36.123},
		{22.307,24.996,27.488,30.578,37.697},
		{23.542,26.296,28.845,32,39.252},
		{24.769,27.587,30.191,33.409,40.79},
		{25.989,28.869,31.526,34.805,42.312},
		{27.204,30.144,32.852,36.191,43.82},
		{28.412,31.41,34.17,37.566,45.315},
		{29.615,32.671,35.479,38.932,46.797},
		{30.813,33.924,36.781,40.289,48.268},
		{32.007,35.172,38.076,41.638,49.728},
		{33.196,36.415,39.364,42.98,51.179},
		{34.382,37.652,40.646,44.314,52.62},
		{35.563,38.885,41.923,45.642,54.052},
		{36.741,40.113,43.195,46.963,55.476},
		{37.916,41.337,44.461,48.278,56.892},
		{39.087,42.557,45.722,49.588,58.301},
		{40.256,43.773,46.979,50.892,59.703},
		{41.422,44.985,48.232,52.191,61.098},
		{42.585,46.194,49.48,53.486,62.487},
		{43.745,47.4,50.725,54.776,63.87},
		{44.903,48.602,51.966,56.061,65.247},
		{46.059,49.802,53.203,57.342,66.619},
		{47.212,50.998,54.437,58.619,67.985},
		{48.363,52.192,55.668,59.893,69.347},
		{49.513,53.384,56.896,61.162,70.703},
		{50.66,54.572,58.12,62.428,72.055},
		{51.805,55.758,59.342,63.691,73.402},
		{52.949,56.942,60.561,64.95,74.745},
		{54.09,58.124,61.777,66.206,76.084},
		{55.23,59.304,62.99,67.459,77.419},
		{56.369,60.481,64.201,68.71,78.75},
		{57.505,61.656,65.41,69.957,80.077},
		{58.641,62.83,66.617,71.201,81.4},
		{59.774,64.001,67.821,72.443,82.72},
		{60.907,65.171,69.023,73.683,84.037},
		{62.038,66.339,70.222,74.919,85.351},
		{63.167,67.505,71.42,76.154,86.661},
		{64.295,68.669,72.616,77.386,87.968},
		{65.422,69.832,73.81,78.616,89.272},
		{66.548,70.993,75.002,79.843,90.573},
		{67.673,72.153,76.192,81.069,91.872},
		{68.796,73.311,77.38,82.292,93.168},
		{69.919,74.468,78.567,83.513,94.461},
		{71.04,75.624,79.752,84.733,95.751},
		{72.16,76.778,80.936,85.95,97.039},
		{73.279,77.931,82.117,87.166,98.324},
		{74.397,79.082,83.298,88.379,99.607},
		{75.514,80.232,84.476,89.591,100.888},
		{76.63,81.381,85.654,90.802,102.166},
		{77.745,82.529,86.83,92.01,103.442},
		{78.86,83.675,88.004,93.217,104.716},
		{79.973,84.821,89.177,94.422,105.988},
		{81.085,85.965,90.349,95.626,107.258},
		{82.197,87.108,91.519,96.828,108.526},
		{83.308,88.25,92.689,98.028,109.791},
		{84.418,89.391,93.856,99.228,111.055},
		{85.527,90.531,95.023,100.425,112.317},
		{86.635,91.67,96.189,101.621,113.577},
		{87.743,92.808,97.353,102.816,114.835},
		{88.85,93.945,98.516,104.01,116.092},
		{89.956,95.081,99.678,105.202,117.346},
		{91.061,96.217,100.839,106.393,118.599},
		{92.166,97.351,101.999,107.583,119.85},
		{93.27,98.484,103.158,108.771,121.1},
		{94.374,99.617,104.316,109.958,122.348},
		{95.476,100.749,105.473,111.144,123.594},
		{96.578,101.879,106.629,112.329,124.839},
		{97.68,103.01,107.783,113.512,126.083},
		{98.78,104.139,108.937,114.695,127.324},
		{99.88,105.267,110.09,115.876,128.565},
		{100.98,106.395,111.242,117.057,129.804},
		{102.079,107.522,112.393,118.236,131.041},
		{103.177,108.648,113.544,119.414,132.277},
		{104.275,109.773,114.693,120.591,133.512},
		{105.372,110.898,115.841,121.767,134.746},
		{106.469,112.022,116.989,122.942,135.978},
		{107.565,113.145,118.136,124.116,137.208},
		{108.661,114.268,119.282,125.289,138.438},
		{109.756,115.39,120.427,126.462,139.666},
		{110.85,116.511,121.571,127.633,140.893},
		{111.944,117.632,122.715,128.803,142.119},
		{113.038,118.752,123.858,129.973,143.344},
		{114.131,119.871,125,131.141,144.567},
		{115.223,120.99,126.141,132.309,145.789},
		{116.315,122.108,127.282,133.476,147.01},
		{117.407,123.225,128.422,134.642,148.23},
		{118.498,124.342,129.561,135.807,149.449}
	};
	
	
	/**
	 * Rows: number of trials. Columns: significance level 0.10, 0.05, 0.02, 0.01. 
	 * Cell $(i,j)$: maximum absolute difference with $i+1$ trials and significance $j$,
	 * for a continuous distribution.
	 *
	 * Source: Patrick D. T. O'Connor and Andre Kleyner, "Practical Reliability 
	 * Engineering", Fifth Edition. Wiley 2012.
	 *<http://onlinelibrary.wiley.com/store/10.1002/9781119961260.app3/asset/app3.pdf?v=1&t=j3tyc55h&s=3defbd3058d1520b202a2d56d7d71ba0b25e611f>.
	 */
	public static final double[][] KOLMOGOROV_SMIRNOV_ABSOLUTE = new double[][] {
		{0.95000,0.97500,0.99000,0.99500},
		{0.77639,0.84189,0.90000,0.92929},
		{0.63604,0.70760,0.78456,0.82900},
		{0.56522,0.62394,0.68887,0.73424},
		{0.50945,0.56328,0.62718,0.66853},
		{0.46799,0.51926,0.57741,0.61661},
		{0.43607,0.48342,0.53844,0.57581},
		{0.40962,0.45427,0.50654,0.54179},
		{0.38746,0.43001,0.47960,0.51332},
		{0.36866,0.40925,0.45662,0.48893},
		{0.35242,0.39122,0.43670,0.46770},
		{0.33815,0.37543,0.41918,0.44905},
		{0.32549,0.36143,0.40362,0.43247},
		{0.31417,0.34890,0.38970,0.41762},
		{0.30397,0.33760,0.37713,0.40420},
		{0.29472,0.32733,0.36571,0.39201},
		{0.28627,0.31796,0.35528,0.38086},
		{0.27851,0.30936,0.34569,0.37062},
		{0.27136,0.30143,0.33685,0.36117},
		{0.26473,0.29408,0.32866,0.35241},
		{0.25858,0.28724,0.32104,0.34427},
		{0.25283,0.28087,0.31394,0.33666},
		{0.24746,0.27490,0.30728,0.32954},
		{0.24242,0.26931,0.30104,0.32286},
		{0.23768,0.26404,0.29516,0.31657},
		{0.23320,0.25907,0.28962,0.31064},
		{0.22898,0.25438,0.28438,0.30502},
		{0.22497,0.24993,0.27942,0.29971},
		{0.22117,0.24571,0.27471,0.29466},
		{0.21756,0.24170,0.27023,0.28987},
		{0.21412,0.23788,0.26596,0.28530},
		{0.21085,0.23424,0.26189,0.28094},
		{0.20771,0.23076,0.25801,0.27677},
		{0.20472,0.22743,0.25429,0.27279},
		{0.20185,0.22425,0.26073,0.26897},
		{0.19910,0.22119,0.24732,0.26532},
		{0.19646,0.21826,0.24404,0.26180},
		{0.19392,0.21544,0.24089,0.25843},
		{0.19148,0.21273,0.23786,0.25518},
		{0.18913,0.21012,0.23494,0.25205}
	};
	
	
	/**
	 * Constants used to approximate table $KOLMOGOROV_SMIRNOV_ABSOLUTE$ for a number of 
	 * samples greater than 40, for each significance level.
	 *
	 * Source: Patrick D. T. O'Connor and Andre Kleyner, "Practical Reliability 
	 * Engineering", Fifth Edition. Wiley 2012.
	 *<http://onlinelibrary.wiley.com/store/10.1002/9781119961260.app3/asset/app3.pdf?v=1&t=j3tyc55h&s=3defbd3058d1520b202a2d56d7d71ba0b25e611f>.
	 */
	public static final double[] KOLMOGOROV_SMIRNOV_ABSOLUTE_LARGE = new double[] {1.22,1.36,1.51,1.63};
	
	
	/**
	 * Shared source of randomness
	 */
	public static final Random random = new Random(0);  // The seed is fixed, for reproducibility.
	
	
	/**
	 * Mass of the normal distribution inside positive integer multiples of the standard 
	 * deviation.
	 */
	public static final double[] FRACTION_IN_STANDARD_DEVIATIONS = new double[] {
		0.682689492137,
		0.954499736104,
		0.997300203937,
		0.999999426697
	};
	
	
	/**
	 * @param n positive;
	 * @param d positive;
	 * @return ceil(n/d).
	 */
	public static final int ceil(int n, int d) {
		return 1+(n-1)/d;
	}
	
	
	public static final long ceilLong(long n, long d) {
		return 1+(n-1)/d;
	}
	
	
	/**
	 * @return $x-y$ if $x \geq y$, 0 otherwise.
	 */
	private static final int doz(int x, int y) {
		final int d = x-y;
		return d & (bitwiseEquivalence(d,(x^y)&(d^x))>>31);
	}
	
	
	private static final int bitwiseEquivalence(int x, int y) {
		return ~(x^y);
	}
	
	
	public static final int min(int x, int y) {
		return x-doz(x,y);
	}
	
	
	public static final int max(int x, int y) {
		return y+doz(x,y);
	}

	
	private static final long dozLong(long x, long y) {
		final long d = x-y;
		return d & (bitwiseEquivalenceLong(d,(x^y)&(d^x))>>63);
	}
	
	
	private static final long bitwiseEquivalenceLong(long x, long y) {
		return ~(x^y);
	}
	
	
	public static final long minLong(long x, long y) {
		return x-dozLong(x,y);
	}
	
	
	public static final long maxLong(long x, long y) {
		return y+dozLong(x,y);
	}

	
	/**
	 * @param x nonnegative;
	 * @param y nonnegative;
	 * @return $|x-y|$.
	 */
	public static final int abs(int x, int y) {
		return max(x,y)-min(x,y);
	}
	
	
	public static final int round(int n, int d) {
		return (n+(d>>1))/d;
	}
	
	
	/**
	 * @return $value$, rounded to the $nDigits$-th decimal digit.
	 */
	public static final double round(double value, int nDigits) {
		double scale = 1.0;
		for (int i=0; i<nDigits; i++) scale*=10.0;
		return Math.round(value*scale)/scale;
	}
	
	
	/**
	 * Stores in a prefix of $to$ the union of the sorted sets $from1[0..last1]$ and 
	 * $from2[0..last2]$.
	 *
	 * @return the last element of $to$.
	 */
	public static final int setUnion(int[] from1, int last1, int[] from2, int last2, int[] to) {
		int i1, i2, j;
		
		i1=0; i2=0; j=0;
		while (i1<=last1 && i2<=last2) {
			if (from1[i1]<from2[i2]) to[j++]=from1[i1++];
			else if (from1[i1]>from2[i2]) to[j++]=from2[i2++];
			else { 
				to[j++]=from1[i1];
				i1++; i2++;
			}
		}
		while (i1<=last1) to[j++]=from1[i1++];
		while (i2<=last2) to[j++]=from2[i2++];
		return j-1;
	}
	
	
	/**
	 * Identical to the above, but for doubles.
	 */
	public static final int setUnion(double[] from1, int last1, double[] from2, int last2, double[] to) {
		int i1, i2, j;
		
		i1=0; i2=0; j=0;
		while (i1<=last1 && i2<=last2) {
			if (from1[i1]<from2[i2]) to[j++]=from1[i1++];
			else if (from1[i1]>from2[i2]) to[j++]=from2[i2++];
			else { 
				to[j++]=from1[i1];
				i1++; i2++;
			}
		}
		while (i1<=last1) to[j++]=from1[i1++];
		while (i2<=last2) to[j++]=from2[i2++];
		return j-1;
	}
	
	
	/**
	 * Stores in $y[fromY..]$ the intersection of the sorted sets $x1[from1..last1]$ and 
	 * $x2[from2..last2]$.
	 *
	 * @return the last element of $y$.
	 */
	public static final int setIntersection(int[] x1, int from1, int last1, int[] x2, int from2, int last2, int[] y, int fromY) {
		int i1, i2, j;
		
		i1=from1; i2=from2; j=fromY;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]<x2[i2]) i1++;
			else if (x1[i1]>x2[i2]) i2++;
			else {
				y[j++]=x1[i1];
				i1++; i2++;
			}
		}
		return j-1;
	}
	
	
	/**
	 * A variant of $setIntersection()$ that just measures intersection size.
	 *
	 * @param exclude integer $exclude$ does not contribute to the count.
	 */
	public static final int setIntersectionSize(int[] x1, int from1, int last1, int[] x2, int from2, int last2, int exclude) {
		int i1, i2, size;
		
		i1=from1; i2=from2; size=0;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]==exclude) {
				i1++;
				continue;
			}
			if (x2[i2]==exclude) {
				i2++;
				continue;
			}
			if (x1[i1]<x2[i2]) i1++;
			else if (x1[i1]>x2[i2]) i2++;
			else {
				size++;
				i1++; i2++;
			}
		}
		return size;
	}
	
	
	public static final boolean setIdentity(int[] x1, int from1, int last1, int[] x2, int from2, int last2) {
		final int length1 = last1-from1+1;
		if (length1!=last2-from2+1) return false;
		for (int i=0; i<length1; i++) {
			if (x1[from1+i]!=x2[from2+i]) return false;
		}
		return true;
	}
	
	
	/**
	 * A variant of $setIntersectionSize()$ that just tests for nonempty intersection.
	 */
	public static final boolean nonemptyIntersection(int[] x1, int from1, int last1, int[] x2, int from2, int last2) {
		int i1, i2;
		
		i1=from1; i2=from2;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]<x2[i2]) i1++;
			else if (x1[i1]>x2[i2]) i2++;
			else return true;
		}
		return false;
	}
	
	
	/**
	 * A variant of $nonemptyIntersection()$ in which element $exclude$ does not 
	 * contribute to the intersection.
	 */
	public static final boolean nonemptyIntersection(int[] x1, int from1, int last1, int[] x2, int from2, int last2, int exclude) {
		int i1, i2;
		
		i1=from1; i2=from2;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]==exclude) {
				i1++;
				continue;
			}
			if (x2[i2]==exclude) {
				i2++;
				continue;
			}
			if (x1[i1]<x2[i2]) i1++;
			else if (x1[i1]>x2[i2]) i2++;
			else return true;
		}
		return false;
	}
	
	
	/**
	 * Like $setIntersectionSize()$, but considers only the elements of $x2$ that are 
	 * flagged in $include2$.
	 */
	public static final int setIntersectionSize(int[] x1, int from1, int last1, int[] x2, int from2, int last2, int exclude, boolean[] include2) {
		int i1, i2, size;
		
		i1=from1; i2=from2; size=0;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]==exclude) {
				i1++;
				continue;
			}
			if (x2[i2]==exclude || !include2[i2]) {
				i2++;
				continue;
			}
			if (x1[i1]<x2[i2]) i1++;
			else if (x1[i1]>x2[i2]) i2++;
			else {
				size++;
				i1++; i2++;
			}
		}
		return size;
	}
	
	
	/**
	 * @return the position of the last element of the sorted set $x2[from2..last2]$ that 
	 * occurs also in the sorted set $x1[from1..last1]$; $from2-1$ if no such element 
	 * exists.
	 */
	public static final int lastInIntersection(int[] x1, int from1, int last1, int[] x2, int from2, int last2) {
		int i1, i2, last;
		
		i1=from1; i2=from2; last=from2-1;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]<x2[i2]) i1++;
			else if (x1[i1]>x2[i2]) i2++;
			else {
				last=i2;
				i1++; i2++;
			}
		}
		return last;
	}
	
	
	/**
	 * Stores in a prefix of $to$ the subtraction of the sorted set $from2[0..last2]$ from 
	 * the sorted set $from1[0..last1]$.
	 *
	 * @return the last element of $to$.
	 */
	public static final int setMinus(int[] from1, int last1, int[] from2, int last2, int[] to) {
		int i1, i2, j;
		
		i1=0; i2=0; j=0;
		while (i1<=last1 && i2<=last2) {
			if (from1[i1]<from2[i2]) to[j++]=from1[i1++];
			else if (from1[i1]>from2[i2]) i2++;
			else { i1++; i2++; }
		}
		while (i1<=last1) to[j++]=from1[i1++];
		return j-1;
	}
	
	
	/**
	 * A variant of $setMinus()$ that just measures the size of the difference.
	 */
	public static final int setMinusSize(int[] from1, int last1, int[] from2, int last2) {
		int i1, i2, out;
		
		i1=0; i2=0; out=0;
		while (i1<=last1 && i2<=last2) {
			if (from1[i1]<from2[i2]) out++;
			else if (from1[i1]>from2[i2]) i2++;
			else { i1++; i2++; }
		}
		while (i1<=last1) out++;
		return out;
	}
	
	
	/**
	 * A variant of $setMinus()$ that just tests for nonempty difference.
	 */
	public static final boolean nonemptySetMinus(int[] from1, int last1, int[] from2, int last2) {
		int i1, i2;
		
		i1=0; i2=0;
		while (i1<=last1 && i2<=last2) {
			if (from1[i1]<from2[i2]) return true;
			else if (from1[i1]>from2[i2]) i2++;
			else { i1++; i2++; }
		}
		return i1<=last1;
	}

	
	/**
	 * Takes columns $c1$ and $c2$ from $from$ and copies them, swapped, into $to$.
	 */
	public static final void swapColumns(int[][] from, int c1, int c2, int[][] to) {
		for (int i=0; i<from.length; i++) {
			to[i][c1]=from[i][c2]; 
			to[i][c2]=from[i][c1];
		}
	}
	
	
	/**
	 * Sets all elements of $matrix$ to the same value $c$.
	 */
	public static final void set(int[][] matrix, int c) {
		for (int i=0; i<matrix.length; i++) Arrays.fill(matrix[i],c);
	}
	
	public static final void set(int[][] matrix, int lastRow, int c) {
		for (int i=0; i<=lastRow; i++) Arrays.fill(matrix[i],c);
	}
	
	
	/**
	 * Sets all elements of $matrix$ to the same value $c$.
	 */
	public static final void set(boolean[][] matrix, boolean c) {
		for (int i=0; i<matrix.length; i++) Arrays.fill(matrix[i],c);
	}
	
	
	public static final void set(boolean[][] matrix, int lastRow, boolean c) {
		for (int i=0; i<=lastRow; i++) Arrays.fill(matrix[i],c);
	}
	
	
	
	public static final void set(double[][] matrix, double c) {
		for (int i=0; i<matrix.length; i++) Arrays.fill(matrix[i],c);
	}
	
	
	/**
	 * Sets all elements of $matrix[0..last]$ to the same value $c$.
	 */
	public static final void set(boolean[] matrix, int last, boolean c) {
		Arrays.fill(matrix,0,last+1,c);
	}
	
	
	/**
	 * Sets all elements of $matrix[0..last]$ to the same value $c$.
	 */
	public static final void set(int[] matrix, int last, int c) {
		Arrays.fill(matrix,0,last+1,c);
	}
	
	
	/**
	 * Sets all elements of $matrix[0..last]$ to the same value $c$.
	 */
	public static final void set(byte[] matrix, int last, byte c) {
		Arrays.fill(matrix,0,last+1,c);
	}
	
	
	/**
	 * Sets all elements of $matrix[first..last]$ to the same value $c$.
	 */
	public static final void set(boolean[] matrix, int first, int last, boolean c) {
		Arrays.fill(matrix,first,last+1,c);
	}
	
	
	/**
	 * Sets all elements of $matrix[first..last]$ to the same value $c$.
	 */
	public static final void set(int[] matrix, int first, int last, int c) {
		Arrays.fill(matrix,first,last+1,c);
	}
	
	/**
	 * Sets all elements of $matrix[0..last]$ to the same value $c$.
	 */
	public static final void set(long[] matrix, int last, long c) {
		Arrays.fill(matrix,0,last+1,c);
	}
	
	
	/**
	 * Sets all elements of $matrix[0..last]$ to the same value $c$.
	 */
	public static final void set(double[] matrix, int last, double c) {
		Arrays.fill(matrix,0,last+1,c);
	}
	
	
	public static final void set(long[][] matrix, long c) {
		for (int i=0; i<matrix.length; i++) Arrays.fill(matrix[i],c);
	}
	
	
	/**
	 * Let $source[0..sourceLast]$ and $subset[0..subsetLast]$ be two sorted lists of
	 * distinct numbers, such that every element in $subset[0..subsetLast]$ is contained 
	 * in $source[0..sourceLast]$. The procedure sets $flags[i]$ to TRUE iff $source[i]$
	 * is present in $subset[0..subsetLast]$.
	 *
	 * @return the total number of flags set in $flags$ after the procedure completes.
	 */
	public static final int set(int[] source, int sourceLast, int[] subset, int subsetLast, boolean[] flags) {
		int i1, i2, count;
		
		i1=0; i2=0; count=0;
		while (i1<=sourceLast && i2<=subsetLast) {
			if (source[i1]<subset[i2]) {
				if (flags[i1]) count++;
				i1++;
			}
			else if (source[i1]>subset[i2]) i2++;
			else { 
				flags[i1]=true;
				count++;
				i1++; i2++;
			}
		}
		return count;
	}
	
	
	/**
	 * Checks whether all elements of $matrix$ are equal to the same value $c$.
	 */
	public static final boolean equals(int[][] matrix, int c) {
		int i, j;
		
		for (i=0; i<matrix.length; i++) {
			for (j=0; j<matrix[i].length; j++) {
				if (matrix[i][j]!=c) return false;
			}
		}
		return true;
	}
	
	
	/**
	 * @param array assumed to be sorted;
	 * @param fromIndex,toIndex,output same conventions as $Arrays.binarySearch$;
	 * @param negative considers both $key$ and $-1-key$ as matches.
	 */
	public static final int linearSearch_sorted(int[] array, int fromIndex, int toIndex, int key, boolean negative) {
		for (int i=fromIndex; i<toIndex; i++) {
			if (array[i]==key || (negative&&array[i]==-1-key)) return i;
			if (array[i]>key) return -1-i;
		}
		return -1-toIndex;
	}
	
	
	/**
	 * @param array NOT assumed to be sorted;
	 * @param fromIndex,toIndex same conventions as $Arrays.binarySearch$;
	 * @param negative considers both $key$ and $-1-key$ as matches;
	 * @return $i>=0$ if found, -1 if not found.
	 */
	public static final int linearSearch_unsorted(int[] array, int fromIndex, int toIndex, int key, boolean negative) {
		for (int i=fromIndex; i<toIndex; i++) {
			if (array[i]==key || (negative&&array[i]==-1-key)) return i;
		}
		return -1;
	}
	
	
	/**
	 * @return $\sum_{i=1}^n i$
	 */
	public static final int summation(int n) {
		return (n*(n+1))>>1;
	}
	
	
	/**
	 * @return if $X$ is the number represented by the string $c.W * 10^{-i}$, where $c$ 
     * is a character in $[1..9]$, $W$ is a sequence of characters, and $i>0$, the 
	 * procedure returns the number represented by the string $ceil(Y) * 10^{-i}$, where 
	 * $Y=c.W$.
	 */
	public static final double ceilPrime(double x) {
		int i, j;
		
		i=0;
		while (x<1) {
			x*=10;
			i++;
		}
		x=ceil(x);
		for (j=0; j<i; j++) x/=10;
		return x;
	}
	
	
	/**
     * @return $\ceil{\log_{2}(x)}$
     */
    public static final int log2(int x) {
    	return 32-Integer.numberOfLeadingZeros(x-1);
    }
	
	
	public static final int max(int[] array, int last) {
		int out = NEGATIVE_INFINITY;
		for (int i=0; i<=last; i++) {
			if (array[i]>out) out=array[i];
		}
		return out;
	}
	
	
	/**
	 * Uses two's complement. Taken from:
	 * http://graphics.stanford.edu/~seander/bithacks.html
	 */
	public static final int abs(int x) {
		final int MASK = x>>31;
		return (x+MASK)^MASK;
	}
	
	
	/**
	 * Copies all cells of $from$ into $to$, which is assumed to be large enough but
	 * possibly bigger.
	 */
	public static final void copy(int[][] from, int[][] to) {
		int i, j;
		int nColumns;
		final int nRows = from.length;
		
		for (i=0; i<nRows; i++) {
			nColumns=from[i].length;
			for (j=0; j<nColumns; j++) to[i][j]=from[i][j];
		}
	}
	
	
	/**
	 * @return the smallest, strictly-positive value in $array[0..last]$; -1 if no such
	 * value can be found.
	 */
	public static final int minPositive(int[] array, int last) {
		int i, out;
		
		out=POSITIVE_INFINITY;
		for (i=0; i<=last; i++) {
			if (array[i]<=0) continue;
			out=min(out,array[i]);
		}
		return out==POSITIVE_INFINITY?-1:out;
	}
	
	
	/**
	 * The following procedures just wrap $java.lang.Math$.
	 */
	public static final double max(double x, double y) {
		return java.lang.Math.max(x,y);
	}
	
	public static final double min(double x, double y) {
		return java.lang.Math.min(x,y);
	}
	
	public static final double abs(double x) {
		return java.lang.Math.abs(x);
	}
	
	public static final double sqrt(double x) {
		return java.lang.Math.sqrt(x);
	}

	public static final double floor(double x) {
		return java.lang.Math.floor(x);
	}
	
	public static final double ceil(double x) {
		return java.lang.Math.ceil(x);
	}
	
	public static final int round(double x) {
		return (int)java.lang.Math.round(x);
	}
	
	public static final double pow(double x, double y) {
		return java.lang.Math.pow(x,y);
	}
	
	public static final double exp(double x) {
		return java.lang.Math.exp(x);
	}
	
	public static final double log(double x) {
		return java.lang.Math.log(x);
	}
	
	public static final int popcount(int x) {
		return java.lang.Integer.bitCount(x);
	}
	
	public static final int numberOfTrailingZeros(int x) {
		return java.lang.Integer.numberOfTrailingZeros(x);
	}
	
	
	public static final void ensureSpace(int[][] matrix, int row, int nColumns) {
		final double GROWTH_RATE = 1.5;  // Arbitrary
		if (matrix[row].length>=nColumns) return;
		
		int[] tmpRow = new int[(int)(nColumns*GROWTH_RATE)];
		set(tmpRow,tmpRow.length-1,0);
		System.arraycopy(matrix[row],0,tmpRow,0,matrix[row].length);
		matrix[row]=tmpRow;
	}
	
	
	/**
	 * Transforms every element $x < 0$ of $array[first..last]$ into $-1-x$, sorts, and 
	 * removes duplicates.
	 *
	 * @return the last element of $array$ after the procedure completes.
	 */
	public static final int makePositive(int[] array, int first, int last) {
		int i, j;
		int previous;
		
		if (last<first) return last;
		for (i=first; i<=last; i++) {
			if (array[i]<0) array[i]=-1-array[i];
		}
		if (last>first) Arrays.sort(array,first,last+1);
		j=first; previous=array[first];
		for (i=first+1; i<=last; i++) {
			if (array[i]==previous) continue;
			j++;
			if (j!=i) array[j]=array[i];
			previous=array[i];
		}
		return j;
	}
	
	
	/**
	 * @return a uniformly distributed vector from the simplex of given $dimension$.
	 */
	public static final double[] sampleFromSimplex(int dimension, Random random) {
		int i, j;
		double[] out;
		
		out = new double[dimension];
		if (dimension==1) out[0]=1.0;
		else if (dimension==2) {
			out[0]=random.nextDouble();
			out[1]=1.0-out[0];
		}
		else {
			out[0]=0;
			for (i=1; i<dimension; i++) out[i]=random.nextDouble();
			if (dimension>2) Arrays.sort(out,1,dimension);
			for (i=1; i<dimension; i++) out[i-1]=out[i]-out[i-1];
			out[dimension-1]=1.0-out[dimension-2];
		}
		return out;
	}
	
	
	/**
	 * Perturbs every element of a simplex vector $vector$ by at most $2*absPerturbation$,
	 * while still keeping it in the simplex.
	 */
	public static final double[] perturbInSimplex(double[] vector, double absPerturbation, Random random) {
		final int length = vector.length;
		int i;
		double[] out = new double[length];
		
		out[0]=vector[0];
		for (i=1; i<length; i++) out[i]=out[i-1]+vector[i];
		out[0]+=-absPerturbation+random.nextDouble()*(2*absPerturbation);
		for (i=1; i<length-1; i++) {
			do { out[i]+=-absPerturbation+random.nextDouble()*(2*absPerturbation); }
			while (out[i]<out[i-1]);
		}
		out[length-1]=1.0;
		for (i=length-1; i>0; i--) out[i]-=out[i-1];
		return out;
	}
	
	
	/**
	 * Reverses $array[first..last]$ (endpoints included).
	 */
	public static final void reverse(int[] array, int first, int last) {
		int a = first;
		int b = last;
		final int mid = (last-first+1)>>1;
		for (int i=0; i<mid; i++) {
			array[a]^=array[b]; array[b]^=array[a]; array[a]^=array[b];
			a++; b--;
		}
	}
    
}