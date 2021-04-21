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

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

/**
 * @author Kerstin Scheubert
 */
public class Compound implements java.io.Serializable {
  
  private static final long serialVersionUID = 2765192746354199635L;
  private Molecule formula;
  private String name;
  private int charge;
  private double mass;
  private double focusedMass;
  private String instrument;
  private Molecule formula2;
  
  private List<Spectrum> spectra;
  private List<Peak> ms1peaks;
  
  public Compound(Molecule formula, String name, int charge, double mass,
                  double focusedMass, String instrument, List<Spectrum> spectra){
  this.formula = formula;
  this.name = name;
  this.charge = charge;
  this.mass = mass;
  this.focusedMass = focusedMass;
  this.instrument = instrument;
  this.spectra = spectra; 
  }
  
  public Compound(String name, int charge, double focusedMass, String instrument,
                  List<Spectrum> spectra){
    this(null, name, charge, focusedMass, focusedMass, instrument, spectra);
  }
  
  public Compound(Molecule formula, String name, int charge, double mass,
      						double focusedMass, String instrument, List<Spectrum> spectra, List<Peak> ms1peaks) {
  	this(formula, name, charge, mass, focusedMass, instrument, spectra);
  	this.ms1peaks = ms1peaks;
	}

  public Molecule getFormula(){
    return formula;
  }
  
  public String getName(){
    return name;
  }
  
  public int getCharge(){
    return charge;
  }
  
  public double getMass(){
    return mass;
  }
  
  public double getFocusedMass(){
    return focusedMass;
  }
  
  public String getInstrument(){
    return instrument;
  }
  
  public List<Spectrum> getSpectra(){
    return spectra;
  }
  
  public void setFormula(Molecule formula){
  this.formula = formula;
  }
  public void setName(String name){
  this.name = name;
  }
  public void setCharge(int charge){
  this.charge = charge;
  }
  public void setMass(double mass){
  this.mass = mass;
  }
  public void setFocusedMass(double focusedMass){
  this.focusedMass = focusedMass;
  }
  public void setSpectra(List<Spectrum> spectra){
  this.spectra = spectra;
  }
  
  public List<Peak> getMS1Peaks() {
		return ms1peaks;
	}

	public void setMS1Peaks(List<Peak> ms1peaks) {
		this.ms1peaks = ms1peaks;
	}

	public void calcRelIntensities(){
    for (Spectrum spectrum : spectra){
      spectrum.calcRelIntensities();
    }
  }
  
  public void calculateSNR(){
    for (Spectrum spectrum : spectra){
      spectrum.calculateSNR();
    }
  }
  
  public void calculateRaw(){
    for (Spectrum spectrum : spectra){
      spectrum.calculateRaw();
    }
  } 

  public void copyIntensities(){
    for (Spectrum spectrum : spectra){
      spectrum.copyIntensities();
    }
  }

  public void normIntensities(int max){
    for (Spectrum spectrum : spectra){    	
      spectrum.normIntensities(max);
    }
    Spectrum ms1scan = new Spectrum(0, ms1peaks);
    ms1scan.copyIntensities();
    ms1scan.normIntensities(max);

  }
  public void normIntensities2(int desiredMax){
	  double max = 0.0;
	  for(Spectrum s:spectra){
		  for (Peak p : s.getPeaks()){
			  if (p.getIntensity() > max){
				  max = p.getRelIntensity();
			  }
		  }
	  }
	  double scaling = desiredMax/max;
	  for(Spectrum s:spectra){
		  for (Peak p : s.getPeaks()){
			  p.setRelIntensity(p.getRelIntensity()*scaling);
		  }
		  s.findMinMaxInt();
	  }
	  
}

	public void findIsotopePeaks(Parameters param){
    for (Spectrum spectrum : spectra){
      spectrum.findIsotopePeaks(param);
    }
	}
	
  public void strictFilter(double threshold){
    for (Spectrum spectrum : spectra){
      spectrum.strictFilter(threshold);
    }
  } 
  
  public List<Peak> mergePeaks(Parameters param){
    /*if (spectra.size() == 1){
      Collections.sort(spectra.get(0).getPeaks());
      return spectra.get(0).getPeaks();
    }*/
    List<Peak> result = new LinkedList<Peak>();
    for (Spectrum spectrum : getSpectra()){
      result.addAll(spectrum.getPeaks());      
    }
    Collections.sort(result);
    ListIterator<Peak> left = result.listIterator(), right = null;
    while (left.hasNext()){
      Peak leftPeak = left.next(), mergePeak = leftPeak; //
      right = result.listIterator(left.nextIndex()); //right: one beyond the peaks to merge.
      while (right.hasNext() && Math.abs(mergePeak.getMass()-(mergePeak=right.next()).getMass()) < param.mergeThreshold){
      }
      int numberMerges = right.previousIndex()-left.previousIndex()-1; // one merge less than there are Peaks to be merged
      // if right cannot be moved any more, move right behind the end of the list, thus inc number merges
      // TODO: Solves the Phlorizin bug
      if (!right.hasNext() && result.size() > 2 && Math.abs(result.get(result.size()-1).getMass()-result.get(result.size()-2).getMass())< param.mergeThreshold){
        ++numberMerges;
      }
      for(int i=0; i<numberMerges && left.hasNext(); ++i){
        mergePeak = left.next();
        //System.out.println("Merging: "+leftPeak.getMass()+" "+mergePeak.getMass()+" "+Math.abs(leftPeak.getMass()-mergePeak.getMass())+" "+param.getErrorForMass(leftPeak.getMass()));
        if (param.intenseMerge) {
          leftPeak.setMass(leftPeak.getIntensity() > mergePeak.getIntensity()?leftPeak.getMass():mergePeak.getMass());
        } else {
          leftPeak.setMass((leftPeak.getMass()+mergePeak.getMass())/2);
        }
        leftPeak.setHighestEnergy(Math.max(mergePeak.getHighestEnergy(), leftPeak.getHighestEnergy()));
        leftPeak.setLowestEnergy(Math.min(mergePeak.getLowestEnergy(), leftPeak.getLowestEnergy()));
        leftPeak.setIntensity(Math.max(mergePeak.getIntensity(), leftPeak.getIntensity()));
        leftPeak.setRelIntensity(Math.max(mergePeak.getRelIntensity(), leftPeak.getRelIntensity()));
        leftPeak.setRelIntensityDB(Math.max(mergePeak.getRelIntensityDB(), leftPeak.getRelIntensityDB()));
        leftPeak.addSpectra(mergePeak.getSpectra());
        leftPeak.mergeDecompositions(mergePeak, param);        
        left.remove();
      }
    }
    
    for (Peak p:result){
      p.checkDecompositions(param, charge);
    }
    return result;
  }
/*  public List<Peak> mergePeaks(double threshold, boolean intenseMerge){
    int numberOfPeaks = 0, index = 0, indexSmallestPeak = 0; // index of the smallest peak belonging in the merge group
    boolean mergedsth = true;
    int refPeakIndex= 0, mergeRange = 0;
    Peak refPeak = null, mergePeak = null;
    for (Iterator<Spectrum> iter = spectra.iterator(); iter.hasNext();){
      numberOfPeaks += iter.next().getPeaks().size();
    }
    Peak[] peaks = new Peak[numberOfPeaks];
    for (Iterator<Spectrum> iter = spectra.iterator(); iter.hasNext();){
      for (Iterator<Peak> peakIter = iter.next().getPeaks().iterator(); peakIter.hasNext();){
        peaks[index] = peakIter.next();
        ++index;
      }
    }
    Arrays.sort(peaks);

    double massUntilMerge = peaks[0].getMass(); //the highest mass we currently merge to.
    for (index = 0; index <= peaks.length; ++index){ // for the last merge, equality is necessary
      if (index == peaks.length || peaks[index].getMass() > massUntilMerge + threshold){
        // next mass is not in merge range, merge the smaller ones
        mergedsth = true;
        refPeakIndex = indexSmallestPeak;
        while (mergedsth){
          mergedsth = false; // not yet anything merged in this loop
          refPeak = null;
          while (refPeak == null){
            refPeak = peaks[refPeakIndex];
            ++refPeakIndex;
            // Restart at the beginning of the merge range, til no more peaks are merged. 
            if (refPeakIndex >= index) refPeakIndex = indexSmallestPeak;
          }
          for (mergeRange = indexSmallestPeak; mergeRange < index; ++mergeRange){
            mergePeak = peaks[mergeRange];
            if(mergePeak != null && mergePeak != refPeak && refPeak.isEnergyAdjacent(mergePeak)){
              
              //refPeak.setMass(mergePeak.getIntensity()*mergePeak.getMass() + refPeak.getIntensity()*refPeak.getMass()/mergePeak.getIntensity()+refPeak.getIntensity());
              if (intenseMerge) {
              	refPeak.setMass(refPeak.getIntensity() > mergePeak.getIntensity()?refPeak.getMass():mergePeak.getMass());
              } else {
              	refPeak.setMass((refPeak.getMass()+mergePeak.getMass())/2);
              }
              refPeak.setHighestEnergy(Math.max(mergePeak.getHighestEnergy(), refPeak.getHighestEnergy()));
              refPeak.setLowestEnergy(Math.min(mergePeak.getLowestEnergy(), refPeak.getLowestEnergy()));
              refPeak.setIntensity(Math.max(mergePeak.getIntensity(), refPeak.getIntensity()));
              refPeak.setRelIntensity(Math.max(mergePeak.getRelIntensity(), refPeak.getRelIntensity()));
              refPeak.setRelIntensityDB(Math.max(mergePeak.getRelIntensityDB(), refPeak.getRelIntensityDB()));
              peaks[mergeRange] = null;
              
              mergedsth = true;
            }
          }
          // Restart at the beginning of the merge range, til no more peaks are merged. 
          if (refPeakIndex >= index) refPeakIndex = indexSmallestPeak;
        }
        // set the smallest peak in the new merge range to be this peak
        indexSmallestPeak = index;
      }
      if (index != peaks.length) massUntilMerge = peaks[index].getMass();
    }
    
    List<Peak> mergedPeaks = new ArrayList<Peak>();
    for (index = 0; index < peaks.length; ++index){
      if (peaks[index] != null) mergedPeaks.add(peaks[index]);
    }
    
    return mergedPeaks;
  } */

  public void decomposePeaks(Decomposer decomp, Parameters param) {
    for (Spectrum s : spectra){
      s.decomposePeaks(decomp, param, charge);
    }
    if (ms1peaks != null && !ms1peaks.isEmpty()){
    	ms1peaks.get(0).decompose(decomp, param, charge);
    }
  }

  public void scoreDecompositions(Parameters param) {
    for (Spectrum s : spectra){
      s.scoreDecompositions(param, charge);
    }
    if (ms1peaks != null && !ms1peaks.isEmpty()){
    	//ms1peaks.get(0).setRelIntensity(ms1peaks.get(0).getIntensity());
    	ms1peaks.get(0).scoreDecompositions(param, new Spectrum(0, ms1peaks), charge);
    }
  }

  public void filterDecompositions() {
    for (Spectrum s : spectra){
      s.filterDecompositions();
    }
    if (ms1peaks != null && !ms1peaks.isEmpty()){
    	ms1peaks.get(0).filterDecompositions();
    }
  }

	public Molecule getFormula2() {
		return formula2;
	}

	public void setFormula2(Molecule formula2) {
		this.formula2 = formula2;
	}
}
