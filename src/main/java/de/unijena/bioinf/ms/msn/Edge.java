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
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * @author Kerstin Scheubert
 */
public class Edge implements Comparable<Edge> {

	  private static List<Set<Molecule>> sureLosses;
	  private static List<Molecule> radicalLosses;
	  private static List<Molecule> forbiddenLosses;
	private Vertex from;	
	private Vertex to;
	private double score;
	private boolean deleted;

	public static void init(Parameters param, int multiplier, ElementTable et) throws ElementNotFoundException{
	    String losses = param.neutralLossList;
		sureLosses = new ArrayList< Set<Molecule> >();
		Set<Molecule> firstSet = new TreeSet<Molecule>();
		String[] array = losses.split(" ");
		for (String loss : array){
			firstSet.add(new Molecule(loss, et));			
		}
		sureLosses.add(firstSet);
		for (int mult = 1; mult < multiplier; ++mult){
			Set<Molecule> newSet = new TreeSet<Molecule> ();
			for (Molecule augend : sureLosses.get(sureLosses.size()-1)){
				for (Molecule addend : firstSet){
					newSet.add(augend.add(addend));
				}
			}
			sureLosses.add(newSet);
		}
		   
	    String radicalString = param.radicalLossList;
	    String[] radicalArray = radicalString.split(" ");
	    radicalLosses = new ArrayList<Molecule>(radicalArray.length);
	    for (String loss : radicalArray){
	      radicalLosses.add(new Molecule(loss, et));
	    }

	    String forbiddenString = param.forbiddenLossList;
	    String[] forbiddenArray = forbiddenString.split(" ");
	    forbiddenLosses = new ArrayList<Molecule>(forbiddenArray.length);
	    for (String loss : forbiddenArray){
	      forbiddenLosses.add(new Molecule(loss, et));
	    }
	}

	public Edge(Vertex from, Vertex to, double parentmass, Parameters param) throws ElementNotFoundException {
		this.from = from;
		this.to =to;
		score = to.getScore();
		Molecule toL = to.getLabel();
		Molecule fromL = from.getLabel();
		Molecule loss = fromL.subtract(toL);

		// Score according to H/C-Ratio
		if (param.hcRatioScoring) {
			double tohcRatio = ((double) toL.getElementAbundance("H"))/toL.getElementAbundance("C");    
			double fromhcRatio = ((double) fromL.getElementAbundance("H"))/fromL.getElementAbundance("C");    
			double quot = MathUtils.pdf(tohcRatio, param.hcAverage, param.hcDev)/MathUtils.pdf(fromhcRatio, param.hcAverage, param.hcDev);
			if (quot<1) logAndAddScore(quot);      
		}

		// score according to hetro atom C ratio
		if (param.heteroScoring){
			int fromHetero = fromL.getElementAbundance("N")+fromL.getElementAbundance("O")+fromL.getElementAbundance("P")+fromL.getElementAbundance("S");
			double fromHeteroRatio = ((double) fromHetero/fromL.getElementAbundance("C"));
			int toHetero = toL.getElementAbundance("N")+toL.getElementAbundance("O")+toL.getElementAbundance("P")+toL.getElementAbundance("S");
			double toHeteroRatio = ((double) toHetero/toL.getElementAbundance("C"));
			double quot = MathUtils.pdf(toHeteroRatio, param.heteroAverage, param.heteroDev)/MathUtils.pdf(fromHeteroRatio, param.heteroAverage, param.heteroDev);
			if (quot<1) logAndAddScore(quot);      
		}

		double toDBE = toL.getDBE();
		double fromDBE = fromL.getDBE();
		// Score according to DBE distribution
		if (param.DBEDependentScoring){
			double quot = MathUtils.pdf(toDBE, param.DBEAverage, param.DBEStdDev)/MathUtils.pdf(fromDBE, param.DBEAverage, param.DBEStdDev);
			if (quot<1) logAndAddScore(quot);      
		}

		// DBE should not be <0
		/*if (toDBE < 0) {
    	if (fromDBE > 0) {
    		score *= (0.45-0.1*toDBE);
    	} else if (toDBE < fromDBE) {
    		score *= (0.45-0.1*toDBE)/(0.45-0.1*fromDBE);    		
    	}
    }*/    	

		// Score with fraction of the parentmass.		
		logAndAddScore(1-(loss.getMass()/parentmass));

		// Collision Energy Scoring
		int toLE = to.getPeak().getLowestEnergy();
		int fromLE = from.getPeak().getLowestEnergy(), fromHE = from.getPeak().getHighestEnergy();
		if (toLE < fromLE){
			// Child appeared before parent. Unlikely!
			logAndAddScore(0.1);
		} else if (fromHE+1 < toLE){ // Neither peak can be found in a spectrum between them, thus, cannot be direct descendance
			logAndAddScore(0.1);
		} else if (fromHE < toLE){ // There is no spectrum with both peaks. The transition is too seamless.
			logAndAddScore(0.8);
		}

		boolean lossIsKnown = false;
		int counter = 0;
		for (Set<Molecule> set : sureLosses){
			++counter;
			if (set.contains(loss)) {
				logAndAddScore(10/counter);
				lossIsKnown = true;
				break;
			}
		}

		/*if (loss.size()==1 && loss.iterator().next().equals("H")){ // Bestrafung von Hs
      score = 0;
    }*/

		String mayNotBeSingle = "CN";
		if (loss.size()==1 && mayNotBeSingle.contains(loss.iterator().next())){
			logAndAddScore(0.0001);
		}
		
	    if (forbiddenLosses.contains(loss)){
	        logAndAddScore(0.001);
	      }
	      
	      if (!loss.isNeutral()){
	        if (radicalLosses.contains(loss)){
	          // TODO: Radical Scorings here
	          lossIsKnown = true;
	          logAndAddScore(0.9);
	        } else {
	          logAndAddScore(0.001);
	        }
	      }
	      if (lossIsKnown){
	  			
	  		} else if (!loss.hasPositiveDBE()) { //The loss should obey the seniorRule: DBE >= 0
	  			logAndAddScore(0.25);
	  		} /*else if (fragments != null && fragments.containsKey(noH)) {
	  			//System.out.println("Fragment "+noH);
	  			score *= (1+1000*fragments.get(noH).getScore());
	  		}*/

	  		// Loss has to be neutral, therefore has integer DBE or is one of the rare radicals!		
	  		/*if (!lossIsKnown && (loss.getDBE()*2)%2 != 0){
	  			score *= 0.25;
	  		}*/

	}

	public Edge(Vertex from, Vertex to, double score) {
		this.from = from;
		this.to = to;
		this.score = score;
	}

	public void setFrom(Vertex from) {
		this.from = from;
	}

	public Vertex from() {
		return from;
	}

	public void setTo(Vertex to) {
		this.to = to;
	}

	public Vertex to() {
		return to;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public void logAndAddScore(double score) {
		this.score += Math.log(score);
		//this.score *= score;
	}

	public void setDeleted(boolean deleted) {
		this.deleted = deleted;
	}

	public boolean isDeleted() {
		return deleted;
	}

	public int fromColour(){
		return from().getColour();
	}

	public int fromNumber(){
		return from().getNumber();
	}

	public int toColour(){
		return to().getColour();
	}

	public int toNumber(){
		return to().getNumber();
	}

	public String toString(){
		DecimalFormat df = new DecimalFormat("0.000");
		return from().getLabel().toString()+" -> "+to().getLabel().toString()+" "+df.format(score);
	}

	public String dotString() throws ElementNotFoundException{
		DecimalFormat df = new DecimalFormat("0.000");
		DecimalFormat one = new DecimalFormat("0.0");		
		Molecule loss = from().getLabel().subtract(to().getLabel());
		String result = "";
		result += from().getLabel().toString()+" -> "+to().getLabel().toString();
		result += " [label=\""+loss.toString()+"\\nMass: "+df.format(loss.getMass())+"\\nDBE: "+one.format(loss.getDBE())+"\\nScore: "+df.format(score)+"\"]";
		return result;
	}

	public int compareTo (Edge another){
		if (score > another.score) return -1;
		if (score == another.score) return 0;
		return 1;
	}

	public String getCmlId(){
		String id = from().getLabel().subtract(to().getLabel()).toString();
		id = "ion_"+id;
		return id;
	}

	public static double getMaxEdgeScore(Parameters param){

		// find max Int in all spectra
		// transform to relInt should've happened already
		double maxInt = 0.0;
		for (Spectrum s : param.input.getSpectra()){
			for (Peak p : s.getPeaks()){
				if (p.getRelIntensity() > maxInt){
					maxInt = p.getRelIntensity();
				}
			}
		}

		// Crazy hacking
		// Generate a peak with maximum intensity, which is absolutely accurate
		Peak p = new Peak(param.et.getAll()[0].getMass(), maxInt, 0);
		p.setRelIntensity(maxInt);
		Molecule m = new Molecule(param.et.getAll()[0].getName(), param.et);
		List<Molecule> d = new ArrayList<Molecule>(1);
		d.add(m);
		p.setDecompositions(d);
		List<Peak> pl = new ArrayList<Peak>(0);
		Spectrum s = new Spectrum(0, pl);
		// Some more work here when isotopes available
		p.scoreDecompositions(param, s, 0);

		// Generate an edge which has a sureLoss, and to and from are equal (=> no H/C Scoring etc.)
		List< Set<Molecule> > realSureLosses = sureLosses;
		sureLosses = new ArrayList< Set<Molecule> >(1);
		Set<Molecule> set = new TreeSet<Molecule>();
		set.add(new Molecule(param.et));
		sureLosses.add(set);
		Vertex u = new Vertex(m, p, 0, 0, false), v = new Vertex(m, p, 0, 0, false);
		Edge e = new Edge(u, v, param.input.getMass(), param);
		sureLosses = realSureLosses;
		return e.getScore();
	}

	public String dotString(int[] ces, int charge) {
		DecimalFormat df = new DecimalFormat("0.000");
		DecimalFormat one = new DecimalFormat("0.0");   
		Molecule loss = from().getLabel().subtract(to().getLabel());
		String result = "";
		result += from().getLabel().toString()+" -> "+to().getLabel().toString();
		result += " [label=\""+loss.toString()+(loss.isNeutral()?"":" radical ")+"\\nMass: "+df.format(loss.getMass())+"\\nDBE: "+one.format(loss.getDBE())+"\\nScore: "+df.format(score)+"\"]";		
		return result;
//		String result = "";
//		result += from().getLabel().toString()+" -> "+to().getLabel().toString();		
//		return result;
		
	}

}
