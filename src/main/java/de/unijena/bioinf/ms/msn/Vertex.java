/*
 *  This file is part of the ftree-msn library for analyzing and MSn data
 *
 *  Copyright (C) 2021 Sebastian Böcker and Friedrich-Schiller University, Jena.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with the ftree-msn library. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.ms.msn;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.ArrayList;
import java.util.LinkedList;

/**
 * @author Kerstin Scheubert
 */
public class Vertex implements Comparable<Vertex> {
	private List<Edge> in;
	private List<Edge> out;
	private Molecule label;
	private Peak peak;
	private int colour;
	private int number;
	private boolean PMD;
	//private double minGain;
	//private double maxGain;
	private ColourSetsMap sets;
	private List<Edge> optTree;


	public Vertex(Molecule label, Peak peak, int colour, int number, boolean isPMD) {
		this.label = label;
		this.peak = peak;
		this.colour = colour;
		this.number = number;
		this.PMD = isPMD;
		in = new ArrayList<Edge>(20);
		out = new ArrayList<Edge>(20);
		this.sets = null;
		this.optTree = null;
	}

	public Vertex(Molecule label, boolean isPMD) {
		this(label, null, 0, 0, isPMD);
	}

	public Vertex(Molecule label){
		this(label, false);
	}

	public Molecule getLabel(){
		return label;
	}

	public Peak getPeak(){
		return peak;
	}

	public void setPeak(Peak peak){
		this.peak = peak;
	}

	public int getColour(){
		return colour;
	}

	public void setColour(int colour) {
		this.colour = colour;
	}

	public int getNumber(){
		return number;
	}

	public void setNumber(int number){
		this.number = number;
	}

	public double getScore() {
		return label.getScore();
	}

	public void setScore(double score) {
		label.setScore(score);
	}

	public void addScore(double score) {
		label.addScore(score);
	}

	public String getCmlId(){
		return "ion_" + label;
	}

	public boolean isPMD(){
		return PMD;
	}

	public void initSets(int numberOfColours){
		sets = new ColourSetsMap(numberOfColours, colour);
	}


	/*  public double getMinGain(){
    return minGain;
  }

  public double getMaxGain(){
    return maxGain;
  }

  public boolean getMark(){
    return mark;
  }

  public void setMark(boolean mark){
    this.mark = mark;
  }

  public int getLevel(){
    return level;
  }

  public void setLevel(int level){
    this.level = level;
  } */

	public ColourSetsMap getSets() {
		return sets;
	}

	public void setSets(ColourSetsMap sets) {
		this.sets = sets;
	}

	public void addOutEdge(Edge e){
		out.add(e);
	}

	public void addInEdge(Edge e){
		in.add(e);
	}

	public List<Edge> outEdges(){
		return out;
	}

	public List<Edge> inEdges(){
		return in;
	}

	public int compareTo(Vertex another){
		return getColour() - another.getColour();		
	}

	/*public Object clone(){
		Vertex result = new Vertex(label);
  	result.in = (List<Edge>) in.clone();
  	result.out = (List<Edge>) out.clone();
		result.colour = colour;
		result.number = number;
		result.score = score;
	  result.PMD = PMD;
	  result.isBT = isBT;
	  result.minGain = minGain;
	  result.maxGain = maxGain;
	  result.level = level;
	  result.sets = null;
	  result.setsBT = null;
		return (Object) result;
	}*/

	public void delete(){
		List<Edge> neighEdges;
		boolean done;

		for (int i = 0; i < in.size(); ++i){
			in.get(i).setDeleted(true);
			neighEdges = in.get(i).from().out;
			done = false;
			for (int j = 0; j < neighEdges.size() && !done; ++j){
				if (neighEdges.get(j).to() == this){
					neighEdges.remove(j);
					done = true;
				}
			}
		}
		for (int i = 0; i < out.size(); ++i){
			out.get(i).setDeleted(true);
			neighEdges = out.get(i).to().in;
			done = false;
			for (int j = 0; j < neighEdges.size() && !done; ++j){
				if (neighEdges.get(j).from() == this){
					neighEdges.remove(j);
					done = true;
				}
			}
		}

	}

	/*  public void rankOutEdges(){
  	Collections.sort(out);
  	int i = 0, n = out.size();
  	if (n == 0) return;
  	for (Edge e : out){
  		e.score *= 1-(i/n);
  		++i;
  	}
  } */

	/*  public void calcGains(int[] amountOfColour){
    minGain = 0;
    maxGain = 0;
    int numberOfColours = amountOfColour.length;
    int[] amountOfColourLocal = new int[numberOfColours];
    double maxIn = 0, minIn = Double.POSITIVE_INFINITY, score;
    Edge edge;

    for (Iterator<Edge> eIter = in.iterator(); eIter.hasNext();){
      score = eIter.next().score;
      maxIn = Math.max(score, maxIn);
      minIn = Math.min(score,minIn);
    }
    minGain += minIn;
    maxGain += maxIn;

    double[] mins = new double[numberOfColours], maxs = new double[numberOfColours];
    for (int i=0; i < mins.length; ++i){
      mins[i] = Double.POSITIVE_INFINITY;
    }

    for (Iterator<Edge> eIter = out.iterator(); eIter.hasNext();){
      ++amountOfColourLocal[eIter.next().toColour];
    }

    for (Iterator<Edge> eIter = out.iterator(); eIter.hasNext();){
      edge = eIter.next();
      if (amountOfColourLocal[edge.toColour] == amountOfColour[edge.toColour]
       && mins[edge.toColour] > edge.score) {
        mins[edge.toColour] = edge.score;
      }
      if (maxs[edge.toColour] < edge.score) {
        maxs[edge.toColour] = edge.score;
      }
    }

    for (int i=0; i < mins.length; ++i){
      if (mins[i] != Double.POSITIVE_INFINITY){
        minGain += mins[i];
      }
    }
    for (int i=0; i < maxs.length; ++i){
      maxGain += maxs[i];
    }

  } */

	public void backtrace(boolean noTransitiveClosureScore){

		long bestSet = sets.setWithMaxScore();
		if (bestSet == 0) { //no good set could be found

			optTree = new LinkedList<Edge>();
			return; // use empty tree
		}

		optTree = backtrace(bestSet, noTransitiveClosureScore);	  
	}

	private List<Edge> backtrace(long s, boolean noTransitiveClosureScore){

		int sizeS = BitSet.size(s);
		if (sizeS == 1) {
			if (s != BitSet.singleton(this.colour)){
				System.out.println("Assertion failed: sizeS=1, S != c(v)");
			}
			return new LinkedList<Edge>();
		}

		double actScore = sets.get(s, sizeS);
		if (Double.isNaN(actScore)){
			System.out.println("Error actScore = NaN");
		}

		// first part of recurrence: out edges
		for(Edge vu : out){
			Vertex u = vu.to();
			if (BitSet.contains(u.colour, s)){
				long sMinusV = BitSet.remove(this.colour, s);
				Double newScore = u.sets.get(sMinusV, sizeS-1);				
				double colourscore=(noTransitiveClosureScore)?0:peak.scoreColours(sMinusV);	      
				if (newScore != null && Math.abs(newScore.doubleValue() + vu.getScore() + colourscore - actScore) <= 1e-10){
					List<Edge> tree = u.backtrace(sMinusV,noTransitiveClosureScore);
					tree.add(vu);
					return tree;
				}
			}
		}

		//second part of recurrence: merging
		for (int sizeOfS1 = (int)Math.ceil((sizeS+1)/2.0); sizeOfS1 <= sizeS-1; ++sizeOfS1){
			DoubleIterator vIter = sets.values(sizeOfS1);
			for (LongIterator kIter = sets.keys(sizeOfS1); kIter.hasNext();){ // here no -1 cause its the method sets handles it.
				long S1 = kIter.next(), S2 = BitSet.add(this.colour, BitSet.setMinus(s, S1));
				double scoreS1 = vIter.next(), scoreS2 = sets.get(S2, sizeS-sizeOfS1+1);
				if (!Double.isNaN(scoreS2)){ // S2 might not be found
					double diff = Math.abs(scoreS2 + scoreS1 - actScore);
					if (diff <= 1e-10){ // take care of rounding errors
						List<Edge> tree = backtrace(S1,noTransitiveClosureScore);
						tree.addAll(backtrace(S2,noTransitiveClosureScore));
						return tree;
					}
				}
			}
		}
		System.out.println("Problem: reached end of bt function!");
		System.out.println(BitSet.toString(s)+" "+toString()+" "+actScore);
		System.out.println();
		return new LinkedList<Edge>();
	}

	public List<Edge> optimalTree(boolean noTransitiveClosureScore){

		if (optTree == null) backtrace(noTransitiveClosureScore);
		return optTree;
	}

	public void setOptimalTree(List<Edge> tree){
		optTree = tree;
	}

	public String toString(){
		return label.toString();
	}

	public String dotString() throws ElementNotFoundException{
		DecimalFormat df = new DecimalFormat("0.000");
		DecimalFormat one = new DecimalFormat("0.0");
		return label.toString()+"[label=\""+label.toString()+"\\nMass: "+df.format(label.getMass())+"\\nInt: "+df.format(peak.getIntensity())+"\\nppm: "+df.format((Math.abs(peak.getMass()-label.getMass())-Element.PROTON_MASS)/peak.getMass()*1e6)+"\\nCE:"+peak.getLowestEnergy()+" "+peak.getHighestEnergy()+"\"]";
		//return label.toString()+"[label=\""+label.toString()+"\\nMass: "+df.format(label.getMass())+"\\nDBE: "+one.format(label.getDBE())+"\\nDev: "+df.format(Math.abs(peak.getMass()-label.getMass())-Element.PROTON_MASS)+"\\nCE:"+peak.getLowestEnergy()+" "+peak.getHighestEnergy()+"\"]";
	}

	public boolean equals(Object o) {
		if (o == null || !(o instanceof Vertex)) return false;
		return label.equals(((Vertex) o).label);
	}

	public String dotString(int[] ces, int charge) {
		DecimalFormat df = new DecimalFormat("0.000");		
		StringBuffer result = new StringBuffer(label.toString());
		label.addElement("H", charge);
		double ionmass = label.getMass()-charge*(Element.HYDROGEN_MASS-Element.PROTON_MASS);		
		result.append("[label=\""+label.toString()+"\\nMass: "+df.format(ionmass)+"\\nInt: "+df.format(peak.getIntensity())+"\\nppm: "+df.format(Math.abs(peak.getMass()-ionmass)/peak.getMass()*1e6)+"\\nCE:"+ces[peak.getLowestEnergy()]+"-"+ces[peak.getHighestEnergy()]);		
		result.append(spectraString(df,charge));
		result.append("\"]");
		//result.append("[label=\""+label.toString()+"\"]");
		label.addElement("H", -1*charge);
		return result.toString();
		//return label.toString()+"[label=\""+label.toString()+"\\nMass: "+df.format(label.getMass())+"\\nDBE: "+one.format(label.getDBE())+"\\nDev: "+df.format(Math.abs(peak.getMass()-label.getMass())-Element.PROTON_MASS)+"\\nCE:"+ces[peak.getLowestEnergy()]+" "+ces[peak.getHighestEnergy()]+"\"]";
		
//		DecimalFormat df = new DecimalFormat("0");		
//		StringBuffer result = new StringBuffer(label.toString());
//		label.addElement("H", charge);
//		double ionmass = label.getMass()-charge*(Element.HYDROGEN_MASS-Element.PROTON_MASS);		
//		result.append("[label=\""+df.format(ionmass));		
//		result.append("\"]");
//		label.addElement("H", -1*charge);
//		return result.toString();
	}
	
	private String spectraString(DecimalFormat df, int charge){
		if(peak.getSpectra().isEmpty())return "";
		StringBuffer result=new StringBuffer("\\nSpectra: ");		
		for(Spectrum s:peak.getSpectra()){						
			result.append("("+s.getMSLevel()+"/"+df.format(s.getParentPeak().getMass())+")");
			if(s!=peak.getSpectra().get(peak.getSpectra().size()-1))result.append("\\n");
		}		
		return result.toString();
	}

	
	public static Comparator<Vertex> massComparator(){
		return new Comparator<Vertex>(){
			public int compare(Vertex v1, Vertex v2){
				if (v1.getPeak().getMass() < v2.getPeak().getMass()) return 1;
				else if (v1.getPeak().getMass() == v2.getPeak().getMass()) return 0;
				return -1;
			}      
		};
	}
}
