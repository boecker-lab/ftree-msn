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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Collections;

/**
 * @author Kerstin Scheubert
 */
public class Peak implements Comparable<Peak>, java.io.Serializable {

	private static final long serialVersionUID = 5798716547077002160L;
	private double mass;
	private double intensity;
	private double relIntensityDB;
	private double relIntensity;
	private int lowestEnergy;
	private int highestEnergy;
	private boolean overlaps;

	private List<Peak> isotopePeaks;
	private List<Molecule> decompositions;
	
	private List<Spectrum> spectra;
	/*private double[][] ms3score;
	private HashMap<Peak,HashMap<Peak,Double>> ms3scoreremaining;*/
	public HashMap<Integer,HashMap<Integer,Double>> ms3score;
	private int colour=-1;


	public Peak(double mass, double intensity, double relIntensityDB, int energy){
		this.mass = mass;
		this.intensity = intensity;
		this.relIntensityDB = relIntensityDB;
		relIntensity = 0.0;
		lowestEnergy = energy;
		highestEnergy = energy;
		decompositions = null;
		spectra=new ArrayList<Spectrum>();
	}

	public Peak(double mass, double intensity, int energy){
		this(mass, intensity, 0.0, energy);		
	}

	public double getMass(){
		return mass;
	}

	public double getIntensity(){
		return intensity;
	}
	public double getRelIntensityDB(){
		return relIntensityDB;
	}

	public double getRelIntensity(){
		if (relIntensity == 0.0) return relIntensityDB;
		return relIntensity;
	}

	public void setMass(double mass){
		this.mass = mass;	
	}

	public void setIntensity(double intensity){
		this.intensity = intensity;
	}

	public void setRelIntensityDB(double relIntensityDB){
		this.relIntensityDB = relIntensityDB;
	}

	public void setRelIntensity(double relIntensity){
		this.relIntensity = relIntensity;
	}

	public int getHighestEnergy() {
		return highestEnergy;
	}

	public void setHighestEnergy(int highestEnergy) {
		this.highestEnergy = highestEnergy;
	}

	public int getLowestEnergy() {
		return lowestEnergy;
	}

	public void setLowestEnergy(int lowestEnergy) {
		this.lowestEnergy = lowestEnergy;
	}

	public void setEnergy(int energy){
		this.lowestEnergy = energy;
		this.highestEnergy = energy;		
	}
	
	public void addSpectrum(Spectrum s){
		if(!spectra.contains(s))spectra.add(s);
	}
	
	public void addSpectra(List<Spectrum> ls){
		for(Spectrum s:ls){
			if(!spectra.contains(s))spectra.add(s);
		}
	}
	
	public void setSpectra(List<Spectrum> spectra){
		this.spectra=spectra;
	}
	
	public List<Spectrum> getSpectra(){
		return spectra;
	}
	
	/*public void setMS3Score(double[][] ms3score){
		this.ms3score=ms3score;
	}
	
	public double[][] getMS3Score(){
		return ms3score;
	}
	
	public void setMS3ScoreRemaining(HashMap<Peak,HashMap<Peak,Double>> ms3scoreremaining){
		this.ms3scoreremaining=ms3scoreremaining;
	}
	
	public HashMap<Peak,HashMap<Peak,Double>> getMS3ScoreRemaining(){
		return ms3scoreremaining;
	}*/
	
	public void setMS3Score(HashMap<Integer,HashMap<Integer,Double>> ms3score){
		this.ms3score=ms3score;
	}
	
	public HashMap<Integer,HashMap<Integer,Double>> getMS3ScoreRemaining(){
		return ms3score;
	}
	
	public void setColour(int colour){
		this.colour=colour;
	}
	
	public int getColour(){
		return colour;
	}

	public boolean hasIsotopes(){
		return !(isotopePeaks == null || isotopePeaks.isEmpty());
	}

	public List<Peak> getIsotopes(){
		return isotopePeaks;
	}

	void setIsotopes(List<Peak> p){
		isotopePeaks = p;
	}

	public void addIsotope(Peak p){
		if (isotopePeaks == null){
			isotopePeaks = new ArrayList<Peak>(4);
			isotopePeaks.add(this);
		}
		isotopePeaks.add(p);
		Collections.sort(isotopePeaks);
	}

	public boolean isEnergyAdjacent(Peak another){
		return (highestEnergy+1 >= another.lowestEnergy && lowestEnergy-1 <= another.highestEnergy);
	}

	public boolean overlaps() {
		return overlaps;
	}

	public void setOverlaps(boolean overlaps) {
		this.overlaps = overlaps;
	}

	public int compareTo(Peak p){
		if (this.mass < p.mass) return -1;
		else if (this.mass == p.mass) return 0;
		return 1;
	}

	public void decompose(Decomposer decomp, Parameters param, int mode) {
		if (decompositions == null) {
			decompositions = decomp.decompose(mass, mode);
			Collections.sort(decompositions);
		}

	}

	public List<Molecule> getDecompositions() {
		return decompositions;
	}

	public void setDecompositions(List<Molecule> decompositions) {
		this.decompositions = decompositions;
	}

	public void scoreDecompositions(Parameters param, Spectrum spectrum, int charge) {
		for (Molecule m : decompositions){
			// vertex/molecule scoring here!			
			m.addScore(getRelIntensity()*0.001); // 0-Model: Exponential-Distribution of noise peaks
			//m.logAndAddScore(getRelIntensity()); // 0-Model: pareto distribution of noise peaks
			//m.logAndaddScore(MathUtils.pdf(m.getMass(),this.getMass(),param.getErrorForMass(getMass())/param.massDeviationPenalty)); // Calculates the abs. error from the relative
			m.logAndAddScore(MathUtils.erfc(Math.abs(this.getMass()-m.getMass()-Element.PROTON_MASS*charge)*param.massDeviationPenalty/param.getErrorForMass(getMass())/Math.sqrt(2)));
		}
		if (param.isotopes && this.hasIsotopes()){
			scoreIsotopeDist(param, spectrum, charge);
		}
	}

	public void scoreIsotopeDist(Parameters param, Spectrum spectrum, int charge){
		if (!this.hasIsotopes()) return;

		System.out.println("Scoring isotopes of peak: "+this.getMass());
		List<Peak> iso =  this.getIsotopes();

		IsotopeScorer scorer;
		scorer = new IsotopeScorer(iso, param, spectrum);
		/*if (iso.size() <= 2) {
      scorer = new IsotopeScorer(iso, param, spectrum);
    } else {
      scorer = new HAdductIsotopeScorer(iso, param, spectrum);
      double intSum = 0.0;
      for (Peak p : iso){
      	intSum += p.getRelIntensity();
      }
      for (Peak p : iso){
      	p.setRelIntensityDB(p.getRelIntensity()/intSum);
      }
      System.out.println("Using HadductScorer");
    }*/

		for (Molecule m : decompositions){    
			// Calculate hypothetical Iso pattern
			Molecule ion = new Molecule(m);
			try {
				ion.addElement("H", charge);
			} catch (ElementNotFoundException e) {
				System.out.println("Hydrogen not in element table, aborting isotope scoring!");
				return;
			}
			// Decide whether to score missing last peak or not
			IsotopeDistribution candidate = 
				(scorer instanceof HAdductIsotopeScorer)?
						ion.calculateIsotopeDistribution(iso.size(),0):
							ion.calculateIsotopeDistribution(iso.size(),0);
						// we want to score no additional peak
						double isoScore = scorer.getScore(candidate.getPeaks());
						m.setIsotopeScore(isoScore);	
						// Deferred to MS2Analyzer to check plausibility of isotopes
						//m.addScore(20*isoScore+param.isobonus);
		}
	}

	public void filterDecompositions() {
		Molecule decomposition = null;
		for (Iterator<Molecule> decompIter = decompositions.iterator(); decompIter.hasNext();){
			decomposition = decompIter.next();
			if (!decomposition.hasPositiveDBE()){
				decompIter.remove();
			}
		}    
	}

	public void mergeDecompositions(Peak mergePeak, Parameters param) {
		if (decompositions == null || decompositions.isEmpty()){
			decompositions = mergePeak.decompositions;
			return;
		}

		decompositions.addAll(mergePeak.decompositions);
		Collections.sort(decompositions);
		Iterator<Molecule> iter = decompositions.iterator();
		Molecule prev = null, curr = null; //previous and current Molecules
		while(iter.hasNext()){
			curr = iter.next();
			if(curr.equals(prev)){
				prev.setScore(Math.max(prev.getScore(), curr.getScore()));
				iter.remove();
			} else {
				prev = curr; //change previous only if current is not removed
			}
		}
	}

	public void checkDecompositions(Parameters param, int charge){
		if (decompositions == null){
			return;
		}

		double error = param.getErrorForMass(getMass());
		if (this.overlaps){
			error += 10*1e-6*this.getMass();
		}
		for (Iterator<Molecule> iter = decompositions.iterator(); iter.hasNext();){
			if (Math.abs(iter.next().getMass()-this.getMass()+charge*Element.PROTON_MASS) > error){
				iter.remove();
			}
		}    
	}
	public String toString(){
		return Double.toString(mass)+" "+Double.toString(relIntensity);		
	}

	public static Comparator<Peak> intensityComparator(){
		return new Comparator<Peak>(){
			public int compare(Peak p1, Peak p2){
				if (p1.getIntensity() < p2.getIntensity()) return 1;
				else if (p1.getIntensity() == p2.getIntensity()) return 0;
				return -1;
			}      
		};
	}

	public static Comparator<Peak> relIntensityComparator(){
		return new Comparator<Peak>(){
			public int compare(Peak p1, Peak p2){
				if (p1.getRelIntensity() < p2.getRelIntensity()) return 1;
				else if (p1.getRelIntensity() == p2.getRelIntensity()) return 0;
				return -1;
			}      
		};
	}
	
	public static Comparator<Peak> massComparator(){
		return new Comparator<Peak>(){
			public int compare(Peak p1, Peak p2){
				if (p1.getMass() < p2.getMass()) return 1;
				else if (p1.getMass() == p2.getMass()) return 0;
				return -1;
			}      
		};
	}
	
	public double scoreColours(long colorset){
		/*if(ms3score==null||currcol>ms3score.length)return 0;//TODO: throw error
		double sum=0;
		String s=Long.toBinaryString(colorset);
		int l=s.length();	
		for(int j=l-1;j>=0;j--){
			if(s.charAt(j)=='1'){
				if(ms3score[currcol].length<l-1-j){sum+=0;}//TODO: throw error
				else{sum+=ms3score[currcol][l-1-j];}
			}
		}*/
		
		if(ms3score==null)return 0;
		double sum=0;
		String s=Long.toBinaryString(colorset);
		int l=s.length();	
		for(int j=l-1;j>=0;j--){
			if(s.charAt(j)=='1'){
				if(ms3score.get(colour)!=null&&ms3score.get(colour).get(l-1-j)!=null)sum+=ms3score.get(colour).get(l-1-j);		
			}
		}	
		return sum;		
	}
	
	public double scoreColours(int colour){
		if(ms3score==null)return 0;
		double sum=0;
		if(ms3score.get(this.colour)!=null&&ms3score.get(this.colour).get(colour)!=null)sum=ms3score.get(this.colour).get(colour);		
		//else if(ms3score.get(colour)!=null&&ms3score.get(colour).get(this.colour)!=null)sum=ms3score.get(colour).get(this.colour);
		return sum;		
	}
	
	public double scoreColours(Peak p2){		
		return scoreColours(p2.getColour());
	}

}
