package de.mpi_cbg.revant.util;

import java.awt.image.*;


public class Colors {
	
	public static final int COLOR_PLAIN = 0x00999999;
	
	public static final int COLOR_BACKGROUND = 0x00050A0D;
	public static final int COLOR_BACKGROUND_LIGHT = 0x001B2125; //121619;
	
	public static final int GRADIENT_FROM = 0x00E98532; // Light=low error
	public static final int GRADIENT_TO = 0x00982425;  // Dark=high error
	public static final int GRADIENT_FROM_DESATURATED = 0x00C7C7C7;  // Light=low error
	public static final int GRADIENT_TO_DESATURATED = 0x00636363; //2B2B2B;  // Dark=high error
	
	public static final int COLOR_MAXIMAL = 0x00FEED50; //00FF00;
	public static final int COLOR_MAXIMAL_DESATURATED = 0x00EFEFEF;
	
	public static final int COLOR_GUIDE = 0x00363535; //233E4F;
	
	public static final int COLOR_TEXT = 0x00D9D7D2;
	public static final int COLOR_TEXT_DARK = COLOR_BACKGROUND_LIGHT;
	public static final int COLOR_TEXT_BOUNDINGBOX = 0x00737373;
	
	public static final int COLOR_HISTOGRAM = 0x00509FB5; //3E8AA4; //262626;
	public static final int COLOR_HISTOGRAM_AUXILIARY = 0x2C5B68;
	public static final int COLOR_REFERENCE = 0x00244254;
	public static final int COLOR_REFERENCE_HIGHLIGHT = COLOR_HISTOGRAM; //447B9D;
	public static final int COLOR_REFERENCE_ALIEN = 0x002F8C6D;
	public static final int COLOR_REFERENCE_ALIEN_HIGHLIGHT = 0x0000E89B;
	public static final String COLOR_DOT_REFERENCE_ALIEN_HIGHLIGHT = "#00E89B";
	
	public static final int COLOR_TRACK_LOWCOMPLEXITY = 0x00878787;  //Reads.GRADIENT_TO;  // Same as low-quality
	
	public static final String COLOR_DOT_BACKGROUND = "#050A0D";
	public static final String COLOR_DOT_EDGE = "#636363";
	public static final String COLOR_DOT_TEXT = "#D9D7D2";
	public static final String COLOR_DOT_LONGPERIOD = "#244254";
	
	public static final int COLOR_KERNELTAG = COLOR_REFERENCE;
	public static final int COLOR_KERNELTAG_HIGHLIGHT = COLOR_REFERENCE_HIGHLIGHT;
	public static final int COLOR_KERNELTAG_FULL = COLOR_REFERENCE_ALIEN;
	public static final int COLOR_KERNELTAG_FULL_HIGHLIGHT = COLOR_REFERENCE_ALIEN_HIGHLIGHT;
	public static final int COLOR_KERNELTAG_NOKERNEL = COLOR_BACKGROUND;
	public static final int COLOR_KERNELTAG_NOKERNEL_HIGHLIGHT = COLOR_TEXT_BOUNDINGBOX;
	
	public static final int COLOR_TRACK_TANDEM = type2color(Constants.INTERVAL_PERIODIC);
	
	public static final int FRAME_COLOR = 0x00666666;
	
	
	public static final int type2color(int type) {
		switch (type) {
			case Constants.INTERVAL_ALIGNMENT: return 0x005B201D; //BC6D13;
			case Constants.INTERVAL_DENSE_PREFIX: return 0x00A6231A;
			case Constants.INTERVAL_DENSE_SUFFIX: return 0x00A6231A;
			case Constants.INTERVAL_DENSE_PREFIXSUFFIX: return 0x00A6231A;
			case Constants.INTERVAL_DENSE_SUBSTRING: return 0x00B66B15;
			case Constants.INTERVAL_DENSE_SINGLEDELETION: return 0x00A6231A;
			case Constants.INTERVAL_PERIODIC: return 0x00319FC7;
			default: return 0x00FFFFFF;
		}
	}
	
	
	
	public static final String type2color_string(int type) {
		switch (type) {
			case Constants.INTERVAL_ALIGNMENT: return "yellow";
			case Constants.INTERVAL_DENSE_PREFIX: return "red";
			case Constants.INTERVAL_DENSE_SUFFIX: return "orange";
			case Constants.INTERVAL_DENSE_PREFIXSUFFIX: return "purple";
			case Constants.INTERVAL_DENSE_SUBSTRING: return "blue";
			case Constants.INTERVAL_DENSE_SINGLEDELETION: return "green";
			case Constants.INTERVAL_PERIODIC: return "pink";
			default: return "white";
		}
	}
	
	
	
	/*
	public static final int type2color(int type) {
		switch (type) {
			case Constants.INTERVAL_ALIGNMENT: return 0x00FFE606; //FFFB00;
			case INTERVAL_DENSE_PREFIX: return 0x00BC0D01; //FF2600;
			case INTERVAL_DENSE_SUFFIX: return 0x00FE9A20; //FF9300;
			case INTERVAL_DENSE_PREFIXSUFFIX: return 0x00BC0186; //942192;
			case INTERVAL_DENSE_SUBSTRING: return 0x003E8AA4; //0x0029A3BC; //0433FF;
			case INTERVAL_DENSE_SINGLEDELETION: return 0x009501BC; //00F900;
			case Constants.INTERVAL_PERIODIC: return 0x0001BC23; //FF40FF;
			default: return 0x00FFFFFF;
		}
	}
	*/
	
	public static final int type2color_lighter(int type) {
		switch (type) {
			case Constants.INTERVAL_DENSE_PREFIX: return 0x00FF1A0E;
			case Constants.INTERVAL_DENSE_SUFFIX: return 0x00FFB560;
			case Constants.INTERVAL_DENSE_PREFIXSUFFIX: return 0x00FF4ACC;
			case Constants.INTERVAL_DENSE_SUBSTRING: return 0x0047D1FF;
			case Constants.INTERVAL_DENSE_SINGLEDELETION: return 0x00D24BFF;
			case Constants.INTERVAL_PERIODIC: return 0x001FFF3D;
			default: return 0x00FFFFFF;
		}
	}
	
	
	/*
	 * Not maintained, probably obsolete.
	 */
	public static final int type2color_darker(int type) {
		switch (type) {
			case Constants.INTERVAL_ALIGNMENT: return 0x00454300;
			case Constants.INTERVAL_DENSE_PREFIX: return 0x00420900;
			case Constants.INTERVAL_DENSE_SUFFIX: return 0x00492900;
			case Constants.INTERVAL_DENSE_PREFIXSUFFIX: return 0x00380C37;
			case Constants.INTERVAL_DENSE_SUBSTRING: return 0x0000104A;
			case Constants.INTERVAL_DENSE_SINGLEDELETION: return 0x00003700;
			case Constants.INTERVAL_PERIODIC: return 0x003E003F;
			default: return 0x00000000;
		}
	}
	
	
	public static final String type2color_dot(int type) {
		switch (type) {
			case Constants.INTERVAL_ALIGNMENT: return "#5B201D"; //BC6D13;
			case Constants.INTERVAL_DENSE_PREFIX: return "#A6231A";
			case Constants.INTERVAL_DENSE_SUFFIX: return "#A6231A";
			case Constants.INTERVAL_DENSE_PREFIXSUFFIX: return "#A6231A";
			case Constants.INTERVAL_DENSE_SUBSTRING: return "#B66B15";
			case Constants.INTERVAL_DENSE_SINGLEDELETION: return "#A6231A";
			case Constants.INTERVAL_PERIODIC: return "#319FC7";
			default: return "#FFFFFF";
		}
	}
	
	
	public static final boolean isTypeColor(int color) {
		final int MASK = 0x00FFFFFF;
		
		if ((color&MASK)==(type2color(Constants.INTERVAL_ALIGNMENT)&MASK)) return true;
		if ((color&MASK)==(type2color(Constants.INTERVAL_DENSE_PREFIX)&MASK)) return true;
		if ((color&MASK)==(type2color(Constants.INTERVAL_DENSE_SUFFIX)&MASK)) return true;
		if ((color&MASK)==(type2color(Constants.INTERVAL_DENSE_PREFIXSUFFIX)&MASK)) return true;
		if ((color&MASK)==(type2color(Constants.INTERVAL_DENSE_SUBSTRING)&MASK)) return true;
		if ((color&MASK)==(type2color(Constants.INTERVAL_DENSE_SINGLEDELETION)&MASK)) return true;
		if ((color&MASK)==(type2color(Constants.INTERVAL_PERIODIC)&MASK)) return true;
		return false;
	}
	
	
	public static final void colorIfBackground(int x, int y, int color, BufferedImage image, int nRows, int nColumns) {
		final int MASK = 0x00FFFFFF;
		
		if (x<0 || x>=nColumns || y<0 || y>=nRows) return;
		int previousColor = getRGB(image,x,y);
		if ((previousColor&MASK)==(COLOR_BACKGROUND&MASK) || (previousColor&MASK)==(COLOR_BACKGROUND_LIGHT&MASK)) setRGB(image,x,y,color);
	}
	
	
	public static final int makeLighter(int color) {
		final int DELTA = 50;
		final int MASK = 0x0000FF;
		
		int red = (color&0x00FF0000)>>16;
		red=(red+DELTA)&MASK;
		int green = (color&0x0000FF00)>>8;
		green=(green+DELTA)&MASK;
		int blue = color&0x000000FF;
		blue=(blue+DELTA)&MASK;
		return 0x00000000|(red<<16)|(green<<8)|blue;
	}
	
	
	public static final int makeDarker(int color) {
		final int DELTA = 50;
		final int MASK = 0x0000FF;
		
		int red = (color&0x00FF0000)>>16;
		red=red>DELTA?(red-DELTA)&MASK:0;
		int green = (color&0x0000FF00)>>8;
		green=green>DELTA?(green-DELTA)&MASK:0;
		int blue = color&0x000000FF;
		blue=blue>DELTA?(blue-DELTA)&MASK:0;
		return 0x00000000|(red<<16)|(green<<8)|blue;
	}
	
	
	public static final void setRGB(BufferedImage image, int x, int y, int color) {
		if (x<0 || y<0 || x>=image.getWidth() || y>=image.getHeight()) return;
		image.setRGB(x,y,color);
	}
	
	
	public static final int getRGB(BufferedImage image, int x, int y) {
		final int COLOR_BLACK = 0x00000000;
		
		if (x<0 || y<0 || x>=image.getWidth() || y>=image.getHeight()) return COLOR_BLACK;
		return image.getRGB(x,y);
	}
	
	
	/**
	 * Draws a triangle that points to the right or to the left.
	 *
	 * @param fromX,fromY top corner of the triangle;
	 * @param height height of the triangle;
	 * @param direction TRUE=pointing right; FALSE=pointing left.
	 */
	public static final void drawTriangle(BufferedImage image, int fromX, int fromY, int height, boolean direction, int color) {
		int x, y;
		
		if (direction) {
			for (x=fromX; x<fromX+height/2; x++) {
				for (y=fromY+(x-fromX); y<fromY+height-(x-fromX); y++) image.setRGB(x,y,color);
			}
		}
		else {
			for (x=fromX; x>fromX-height/2; x--) {
				for (y=fromY+(fromX-x); y<fromY+height-(fromX-x); y++) image.setRGB(x,y,color);
			}
		}
	}
	
	
}