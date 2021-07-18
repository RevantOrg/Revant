package de.mpi_cbg.revant.factorize;

import de.mpi_cbg.revant.util.Point;


public class Event extends Point {

	public int nOpen, nClosed;
	public boolean deleted;


	public void clear() {
		super.clear();
		nOpen=0;
		nClosed=0;
		deleted=false;
	}


	public int getMass() {
		return nOpen+nClosed;
	}


	public void sum(Point otherPoint) {
		super.sum(otherPoint);
		Event otherEvent = (Event)otherPoint;
		nOpen+=otherEvent.nOpen;
		nClosed+=otherEvent.nClosed;
	}


	public void copyFrom(Point otherPoint) {
		super.copyFrom(otherPoint);
		Event otherEvent = (Event)otherPoint;
		nOpen=otherEvent.nOpen;
		nClosed=otherEvent.nClosed;
	}


	public String toString() {
		return super.toString()+" nOpen="+nOpen+" nClosed="+nClosed;
	}

}